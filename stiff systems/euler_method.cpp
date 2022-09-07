#include <iostream>
#include <vector>
#include <fstream>

/** More on that on https://dspace.spbu.ru/bitstream/11701/15202/1/method_part_1.pdf, part 9
 */

std::vector<double> split_string_into_array(std::string str) {
    std::vector<double> arr;
    while (str.length() != 0) {
        int ind_of_space = str.find(' ');
        arr.push_back(stod(str.substr(0, ind_of_space)));
        str = str.substr(ind_of_space + 1, str.length() - ind_of_space);
        if (str.find(' ') == -1) {
            if (str[0] != ' ')
                arr.push_back(stod(str));
            break;
        } else if (str.length() == 2 && str[0] == '-') {
            arr.push_back(stod(str));
            break;
        }
    }
    return arr;
}

std::vector<double> euler_method(std::vector<std::vector<double>> a_matrix, double h, std::vector <double> temp_vector){
    for(int i = 0; i < 2; ++i){
        temp_vector[i] *= (1 + (a_matrix[0][i] + a_matrix[1][i])*h);
    }
    return  temp_vector;
}


int main() {
    std::cout << "Computational workshop\n Kizeev Danil\n";
    std::ifstream file("D:\\C++ projects\\comp methods\\stiff systems\\data.txt");
    std::cout << "Opening file with start data: \n";
    if (!file) {
        std::cout << "File not found";
        exit(0);
    } else std::cout << "File has opened successfully\n";
    int n = 2;
    std::vector<std::vector<double>> a_matrix(n, std::vector<double>(n));
    std::string string_curr;
    for(int i = 0; i < n; ++i){
        getline(file, string_curr);
        std::vector<double> vect_temp = split_string_into_array(string_curr);
        for(int j = 0; j <  n; ++j){
            a_matrix[i][j] = vect_temp[j];
        }
    }
    std::cout << "\nMain matrix: " << std::endl;
    for(int i = 0; i < n; ++i){         //main Matrix
        for(int j = 0; j < n; ++j){
            std::cout << a_matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::vector<double> init_conditions(2);
    double h;
    std::cout << "Input h, pls: ";std::cin >> h;
    std::cout << "Input initial conditions: "; std::cin >> init_conditions[0] >> init_conditions[1];
    std::vector <double> result_matrix(2);
    result_matrix = euler_method(a_matrix, h, init_conditions);
    return 0;
}