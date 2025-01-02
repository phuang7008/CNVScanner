"Manages integration with AWS S3, such as for file transfer and bucket/key operations"
import os
import subprocess
import sys
from pathlib import Path

import boto3
from boto3.s3.transfer import S3Transfer
from boto3.s3.transfer import TransferConfig
from botocore.config import Config
from botocore.exceptions import ClientError


def get_bucket_keyname(s3uri):
    """ Given an s3uri in the form 's3://bucket/example/object_key.ext',
    or s3://bucket/example/ (similar to directory),
    RETURNS: bucketname, keyname, and filename OR bucketname, keyname
    """
    if s3uri.endswith("/"):
        s3_bucket = s3uri.strip("s3://").split("/")[0]
        key = s3uri.replace("s3://" + s3_bucket + "/", "")
        return s3_bucket, key
    else:
        s3_bucket = s3uri.strip("s3://").split("/")[0]
        key = "/".join(s3uri.strip("s3://").split("/")[1:])
        filename = key.split("/")[-1]
        return s3_bucket, key, filename


def get_bucket_location(s3_bucket):
    """ Given s3 bucketname,
    RETURNS: s3 bucket location 
    """
    s3_client = boto3.client(service_name="s3")
    response = s3_client.get_bucket_location(Bucket=s3_bucket)
    location = response["LocationConstraint"]

    # for buckets created in the US Standard region, us-east-1
    # the value of LocationConstraint will be null
    default_location = "us-east-1"

    if location:
        return location
    else:
        return default_location


def s3_upload_dir_recursively(s3_bucket, s3_key, local_dir_path):
    """ Uploads directory from local absolute path to s3 object key.
    Upload is recursive (s3_obj/file.ext).
    RETURNS: None
    """
    # ensure there is a trailing '/' in s3_object so upload
    # groups object_keys into one object (similar to directory)
    if not s3_key.endswith("/"):
        s3_key += "/"
    local_dir_path = Path(local_dir_path)
    if local_dir_path.is_dir():
        for p in local_dir_path.rglob("*"):
            if Path(p).is_file():
                filepath = str(Path(*p.parts[1:]))
                mode = "upload"
                final_s3_key = s3_key + filepath
                print("Sending %s to transfer manager" % (filepath), flush=True)
                s3_upload_download_file(mode, s3_bucket, final_s3_key, str(p))
        return
    else:
        raise ValueError("{} must be a valid local directory".format(local_dir_path))


def s3_upload_download_file(mode, s3_bucket, key, local_absolutepath):
    """ Uploads and downloads objects from S3. Key must be file; 
    directory uploads/downloads not handled.
    RETURNS: A transfer manager based on parameters provided
    """
    s3_client = boto3.client(
        service_name="s3", region_name=get_bucket_location(s3_bucket)
    )
    config = TransferConfig(
        multipart_threshold=256 * 1024 * 1024,
        max_concurrency=10,
        max_io_queue=1000,
        io_chunksize=2 * 1024 * 1024,
    )
    transfer = S3Transfer(s3_client, config)

    if mode == "download":
        print(
            "Downloading %s from %s to %s... " % (key, s3_bucket, local_absolutepath),
            flush=True,
        )
        transfer.download_file(s3_bucket, key, local_absolutepath)
        print("Download of %s complete." % (key), flush=True)
        return transfer

    elif mode == "upload":
        # upload with sse
        print(
            "Uploading %s to %s in %s bucket (using sse AES256)... "
            % (local_absolutepath, key, s3_bucket),
            flush=True,
        )
        transfer.upload_file(
            local_absolutepath,
            s3_bucket,
            key,
            extra_args={"ServerSideEncryption": "AES256"},
        )
        print("Upload of %s complete." % (local_absolutepath), flush=True)
        return transfer


def check_object_validity(s3_bucket, key):
    """ heads metadata for s3 object without returning
    object itself -- useful to check if object exists
    in s3
    RETURNS: bool is file does/does not exist
    """
    s3_client = boto3.client(service_name="s3")

    try:
        s3_client.head_object(Bucket=s3_bucket, Key=key)
        print("Object '%s' found." % (key), flush=True)
        return True
    except ClientError as e:
        error_code = e.response["Error"]["Code"]
        if error_code == "404":
            print("Object '%s' does not exist in bucket." % (key), flush=True)
        return False


def run_cmd(cmd):
    """ Run executable from python script. 
    """
    if type(cmd) is list:
        shell = False
        executable = None
        print(subprocess.list2cmdline(cmd), flush=True)
    else:
        shell = True
        executable = "/bin/bash"
        print(cmd, flush=True)
    p = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=shell,
        executable=executable,
    )
    out, err = p.communicate()
    print(out, err, flush=True)
    if p.returncode != 0:
        sys.exit(p.returncode)


def get_index_file(file_for_index, file_for_index_local_copy, func1=get_bucket_keyname, func2=s3_upload_download_file, func3=run_cmd):
    """ Downloads an index file (.bai/.crai/.fai)
    or runs samtools to generate index.
    """
    object_prefix = os.path.splitext(file_for_index)[0]
    s3_bucket, object_key, basename = get_bucket_keyname(object_prefix)

    if file_for_index.endswith(".bam"):
        index_file = object_key + ".bam.bai"
        local_index_file = basename + ".bam.bai"
    elif file_for_index.endswith(".cram"):
        index_file = object_key + ".cram.crai"
        local_index_file = basename + ".cram.crai"
    elif file_for_index.endswith(".fa"):
        index_file = object_key + ".fa.fai"
        local_index_file = basename + ".fa.fai"
    elif file_for_index.endswith(".vcf"):
        index_file = object_key + ".vcf.idx"
        local_index_file = basename + ".vcf.idx"

    if check_object_validity(s3_bucket, index_file):
        mode = "download"
        currentDir = os.getcwd()
        local_path = "/".join([currentDir, local_index_file])
        s3_upload_download_file(mode, s3_bucket, index_file, local_path)
    else:
        if file_for_index_local_copy.endswith(".bam") or file_for_index_local_copy.endswith(".cram"):
            samtools_arg = ["samtools", "index", file_for_index_local_copy]
            print("\nIndexing %s..." %(file_for_index_local_copy), flush=True)
            run_cmd(samtools_arg)
        elif file_for_index_local_copy.endswith(".fa"):
            samtools_arg = ["samtools", "faidx", file_for_index_local_copy]
            print("\nIndexing %s..." %(file_for_index_local_copy), flush=True)
            run_cmd(samtools_arg)


def conditionally_download_alignment(input_file_path, func=get_index_file):
    """ Downloads alignment file (or fasta) if valid s3uri, or returns un-altered
    local path. Contains additional logic to download index file to 
    instance.
    RETURNS: input path
    """
    if input_file_path.startswith("s3://"):
        # break into bucketname and keyname
        s3_bucket, key, filename = get_bucket_keyname(input_file_path)
        # Download file
        # S3Transfer internally verifies object exists in s3 as FILE
        mode = "download"
        currentDir = os.getcwd()
        local_path = "/".join([currentDir, filename])
        s3_upload_download_file(mode, s3_bucket, key, local_path)
        
        get_index_file(input_file_path, local_path)
        return local_path

    return input_file_path


def conditionally_download_s3uri(input_file_path):
    """ Downloads input file if valid s3uri, or returns
    un-altered local path.
    RETURNS: input path
    """
    if input_file_path.startswith("s3://"):
        input_uri = input_file_path
        # break into bucketname and keyname
        s3_bucket, key, filename = get_bucket_keyname(input_uri)
        # Download file
        # S3Transfer internally verifies object exists in s3 as FILE
        mode = "download"
        currentDir = os.getcwd()
        local_abspath = "/".join([currentDir, filename])
        s3_upload_download_file(mode, s3_bucket, key, local_abspath)
        
        return local_abspath

    return input_file_path
