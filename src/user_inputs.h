/*
 * =====================================================================================
 *
 *      Filename:       user_inputs.h
 *
 *      Description:    For the functionalities of user inputs
 *
 *      Version:        1.0
 *      Created:        03/16/2020 04:45:04 PM
 *      Revision:       none
 *      Compiler:       gcc
 *
 *      Author:         Peiming (Peter) Huang (phuang@bcm.edu)
 *      Company:        Baylor College of Medicine
 *
 * =====================================================================================
 */
#ifndef USER_INPUTS_H
#define USER_INPUTS_H

#include <ctype.h>		// for isdigit()
#include <errno.h>
#include <getopt.h>
#include <stdbool.h>	// for bool return type
#include "terms.h"

/**
 * produces the usage information for end users
 */
void usage();

/**
 * process user input options and ensure they are correct!
 * @param argv[]: an array of strings that are used to decide user input options
 */
void processUserOptions(User_Input *user_inputs, int argc, char *argv[]);

/**
 * Output the User Input Options to the end user so he/she can double-check if all the options are correct!
 * @param user_inputs: a variable contains all the user input options
 */
void outputUserInputOptions(User_Input *user_inputs);

/**
 * this method is used to initialize the user input structure
 * @return returns an instance of User_Input upon the memory allocation success. On failure, the entire program will exit!
 */
User_Input * userInputInit();

/**
 * This method is used to clean up the memory after the User_Input variable is done!
 * @param User_Input: a user defined struct
 */
void userInputDestroy(User_Input *user_inputs);

void setupOutputReportFiles(User_Input *user_inputs);

void generateFileName(char *output_dir, char *base_name, char **file_in, char *string_to_append);

#endif //USER_INPUTS_H
