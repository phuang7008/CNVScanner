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

#include "slmseg.hpp"

SLMSeg::SLMSeg(double omega, double eta, double stepeta, double fw) {
    R_.parseEvalQ("library(SLMSeg)");

    omega_ = omega;
    eta_ = eta;
    stepeta_ = stepeta;
    fw_ = fw;

    mi_ = 0.0;
    smu_ = 0.0;
    sepsilon_ = 0.0;
    mw_ = 1;

    set_variables();
}

void SLMSeg::SLM() {
    std::cout << "Parameter estimation...\n";
    param_est_seq();

    std::cout << "Muk estimation...\n";
    muk_est();

    std::cout << "Joint segmentation...\n";
    joint_seg();

    std::cout << "Filtering predictions...\n";
    filter_seg();

    std::cout << "Producing results...\n";
    seg_results();
}

void SLMSeg::HSLM() {
    std::cout << "Parameter estimation...\n";
    param_est_seq();

    std::cout << "Muk estimation...\n";
    muk_est();

    std::cout << "Joint segmentation...\n";
    joint_seg_in();

    std::cout << "Filtering predictions...\n";
    filter_seg();

    std::cout << "Producing results...\n";
    seg_results();
}

void SLMSeg::erase_all_vector_data() {
    signal_data_.clear();
    pos_data_.clear();
}
