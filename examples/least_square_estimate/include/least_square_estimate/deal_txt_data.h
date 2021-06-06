#ifndef DEAL_TXT_DATA_H_
#define DEAL_TXT_DATA_H_

#include <fstream>  
#include <string>  
#include <iostream>  
#include <vector>
#include <string>
#include <cstdlib>
#include <Eigen/Dense>

using namespace Eigen;

namespace deal_txt
{
int get_matrix_size(std::string file_path, unsigned int &row, unsigned int &col);

int get_data_matrix(std::string file_path, unsigned int col, MatrixXd &DATA);

}

#endif