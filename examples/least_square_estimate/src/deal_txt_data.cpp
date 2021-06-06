#include <least_square_estimate/deal_txt_data.h>

namespace deal_txt
{
int get_matrix_size(std::string file_path, unsigned int &row, unsigned int &col)
{
    std::ifstream file_hd;
    file_hd.open(file_path);
    std::string line;  
    
    if (!file_hd)
    {  
        std::cout << "No such file" << std::endl;
        return 0;
    }  
    else
    {  
        unsigned int j = 0;
        unsigned int row_ = 0;
        unsigned int col_ = 0;
        getline(file_hd, line);
        row_++;
        while (line[j]!='\0')
        {
            if (line[j]==' ')
            {
                col_++;
            }
            j++;
        }
        
        while (getline(file_hd, line))
        {
            row_++;
        }
        
        row = row_;
        col = col_;

        file_hd.close();
    }

    return 0;
}

int get_data_matrix(std::string file_path, unsigned int col, MatrixXd &DATA)
{
    std::ifstream file_hd;
    file_hd.open(file_path);
    std::string line;  
    
    if (!file_hd)
    {  
        std::cout << "No such file" << std::endl;
        return 0;
    }  
    else
    {  
        unsigned int i = 0;
        while (getline(file_hd, line))
        {
            std::vector<std::string> data;
            data.resize(col);
            unsigned int j = 0;
            unsigned int k = 0;
            while (line[j]!='\0')
            {
                if (line[j]!=' ')
                {
                    data[k].push_back(line[j]);
                }
                else
                {
                    k++;
                }
                
                j++;
            }

            for (unsigned int j=0; j<col; j++)
            {
                DATA(i,j) = atof(data[j].c_str());
            }

            i++;
        }

        file_hd.close();
    }
}

}