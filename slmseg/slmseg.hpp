/******************************************************************************
    SLMSuite - C++ wrapper
    Copyright (C) 2016  Valerio Orlandini, Alberto Magi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/


#ifndef SLMSEG_HPP
#define SLMSEG_HPP


#include <iostream>
#include <fstream>
#include <vector>

#include <RInside.h>

#include "cctable.hpp"


#define TO_FVEC Rcpp::as< std::vector<double> >
#define TO_FLOAT Rcpp::as< double >


typedef Rcpp::NumericVector Vector;
typedef Rcpp::NumericMatrix Matrix;


// double is a floating point value - float or double
class SLMSeg
{
public:
    SLMSeg(double omega, double eta, double stepeta, double fw);

    //bool set_variables(double omega = 0.3, double eta = 0.0001, unsigned int stepeta = 1000000, unsigned int fw = 0)
                       //unsigned int stepeta = 1000000, unsigned int fw = 1)
    bool set_variables()
    {
        R_.assign(omega_, "omega");
        R_.parseEvalQ("omega<-as.double(omega)");
        R_.assign(eta_, "eta");
        R_.parseEvalQ("eta<-as.double(eta)");
        R_.assign(stepeta_, "stepeta");
        R_.parseEvalQ("stepeta<-as.integer(stepeta)");
        R_.assign(fw_, "fw");
        R_.parseEvalQ("fw<-as.integer(fw)");
        R_.assign(mw_, "mw");
        R_.parseEvalQ("mw<-as.integer(mw)");

        return true;
    }

    bool load_signal_file(const char *file_name, char* chromosome)
    {
        CcTable<unsigned int> *signal_file = new CcTable<unsigned int>();

        if(!signal_file->load_table(file_name))
        {
            delete signal_file;

            return false;
        }

        for (unsigned int r = 0; r < signal_file->nrow(); r++)
        {
            bool dummy;
            char* chr = (char*) calloc(100, sizeof(char));

            try
            {
                strcpy(chr, signal_file->get_cell(r, 0, dummy).c_str());
            }
            catch (...)
            {
                strcpy(chr, "0");
            }

            if (strcmp(chr, chromosome) == 0)
            {
                pos_data_.push_back(std::stoi(signal_file->get_cell(r, 1, dummy)));
                signal_data_.push_back(std::stof(signal_file->get_cell(r, 2, dummy)));
            }
            free(chr);
        }

        delete signal_file;

        return true;
    }

    bool load_data(const std::vector<double> &signal, const std::vector<unsigned int> &pos)
    {
        signal_data_ = signal;


        if (pos.size() == signal_data_.size())
        {
            pos_data_ = pos;
        }
        else
        {
            std::cerr << "Positions vector must have the same size of the Signals one\n";

            return false;
        }


        return true;
    }

    bool param_est_seq()
    {
        Vector m_signal;

        for (unsigned int i = 0; i < signal_data_.size(); i++)
        {
            m_signal.push_back(signal_data_.at(i));
        }

        R_["data_matrix"] = m_signal;
        R_.parseEvalQ("data_matrix<-rbind(data_matrix)");

        R_.parseEval("param_list<-ParamEstSeq(data_matrix, as.double(omega))");
        Vector param_vect = R_.parseEval("as.vector(unlist(param_list))");
        mi_ = param_vect.at(0);
        smu_ = param_vect.at(1);
        sepsilon_ = param_vect.at(2);

        R_.assign(mi_, "mi");
        R_.parseEvalQ("mi<-as.double(mi)");
        R_.assign(smu_, "smu");
        R_.parseEvalQ("smu<-as.double(smu)");
        R_.assign(sepsilon_, "sepsilon");
        R_.parseEvalQ("sepsilon<-as.double(sepsilon)");

        return true;
    }

    void muk_est()
    {
        muk_ = TO_FVEC(R_.parseEval("MukEst(data_matrix, mw)"));

        R_.assign(muk_, "muk");
        R_.parseEvalQ("muk<-rbind(muk)");
    }

    void joint_seg()
    {
        total_pred_break_ = TO_FVEC(R_.parseEval("JointSeg(data_matrix, eta, omega, muk, mi, smu, sepsilon)"));

        R_.assign(total_pred_break_, "total_pred_break");
    }

    void joint_seg_in()
    {
        Vector m_pos;

        for (unsigned int i = 0; i < pos_data_.size(); i++)
        {
            m_pos.push_back(pos_data_.at(i));
        }

        R_["pos"] = m_pos;

        total_pred_break_ = TO_FVEC(R_.parseEval("JointSegIn(data_matrix, eta, omega, muk, mi, smu, sepsilon, pos, stepeta)"));

        R_.assign(total_pred_break_, "total_pred_break");
    }

    void filter_seg()
    {
        total_pred_break_filtered_ = TO_FVEC(R_.parseEval("FilterSeg(total_pred_break, fw)"));

        R_.assign(total_pred_break_filtered_, "total_pred_break_filtered");
    }

    void seg_results()
    {
        data_seg_ = TO_FVEC(R_.parseEval("SegResults(data_matrix, total_pred_break_filtered)"));
    }

    void SLM();

    void HSLM();

    double mw()
    {
        return mw_;
    }

    double mi()
    {
        return mi_;
    }

    double smu()
    {
        return smu_;
    }

    double sepsilon()
    {
        return sepsilon_;
    }

    double omega()
    {
        return omega_;
    }

    double eta()
    {
        return eta_;
    }

    double stepeta()
    {
        return stepeta_;
    }

    double fw()
    {
        return fw_;
    }

    std::vector<double> data_seg()
    {
        return data_seg_;
    }

    std::vector<unsigned int> pos_data()
    {
	return pos_data_;
    }

    std::vector<double> muk()
    {
        return muk_;
    }

    std::vector<double> total_pred_break()
    {
        return total_pred_break_;
    }

    std::vector<double> total_pred_break_filtered()
    {
        return total_pred_break_filtered_;
    }

    void erase_all_vector_data();
private:
    RInside R_;

    std::vector<double> signal_data_;
    std::vector<unsigned int> pos_data_;

    double mw_;
    double mi_;
    double smu_;
    double sepsilon_;
    double omega_;
    double fw_;
    double stepeta_;
    double eta_;

    std::vector<double> total_pred_break_;
    std::vector<double> total_pred_break_filtered_;
    std::vector<double> muk_;
    std::vector<double> data_seg_;
};


#endif // SLMSEG_HPP
