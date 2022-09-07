#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "gauss.cpp"


void transforming(std::vector<std::vector<double>>&, std::vector<double>&);
double infinity_norm_of_vector(std::vector<double>);
double infinity_norm_of_matrix(std::vector<std::vector<double>>);
std::vector<double> simple_iteration_method(const std::vector<std::vector<double>>& , std::vector<double> , int);
std::vector<double> operator-(std::vector<double>, std::vector<double>);
std::vector<double> seidel_method(std::vector<std::vector<double>>, std::vector<double> , int);
std::vector<double> top_relaxation_method(std::vector<std::vector<double>>,  std::vector<double>, int, double);



void transforming(std::vector<std::vector<double>> &A_matrix, std::vector<double> &b_vector){
    int n = b_vector.size();
    for(int i = 0; i < n; ++i){
        double element_on_the_main_diag = A_matrix[i][i];
        for(int j = 0; j < n; ++j){
            if(i == j)
                A_matrix[i][j] = 0;
            else if (fabs(element_on_the_main_diag ) >= 0.000001)
                A_matrix[i][j] /= -element_on_the_main_diag;
        }
        if (fabs(element_on_the_main_diag ) >= 0.000001)
            b_vector[i] /= element_on_the_main_diag;
    }
}

double infinity_norm_of_vector(std::vector<double> vector) {
    double max_el = fabs(vector[0]);
    for(int i = 1; i < vector.size(); ++i){
       if (fabs(vector[i]) > max_el){
           max_el = fabs(vector[i]);
       }
    }
    return max_el;
}

double infinity_norm_of_matrix(std::vector<std::vector<double>> matrix){
    double max_sum = 0;
    for(int i = 0; i < matrix[0].size(); ++i){
        double temp_sum = 0;
        for(int j = 0; j < matrix[0].size(); ++j){
            temp_sum += fabs(matrix[i][j]);
        }
        if (temp_sum > max_sum)
            max_sum = temp_sum;
    }
    return max_sum;
}

std::vector<double> simple_iteration_method(const std::vector<std::vector<double>>& H_matrix, std::vector<double> g_vector, int iteration){
    int n = g_vector.size();
    std::vector<double> x_vector_current(n, 0);
    for(int k = 0; k < iteration; ++k){
        std::vector<double> x_vector_temp = x_vector_current;
        for(int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < n; ++j) {
                sum += H_matrix[i][j] * x_vector_current[j];
            }
            x_vector_temp[i] = sum + g_vector[i];
        }
        x_vector_current = x_vector_temp;
    }
    return x_vector_current;
}

std::vector<double> operator-(std::vector<double> v1, std::vector<double> v2){
    int size1 = v1.size(), size2 = v2.size();
    for(int i = size1; i < size2; ++i){
        v1.push_back(0);
    }
    for(int i = size2; i < size1; ++i){
        v2.push_back(0);
    }
    for(int i = 0; i < v1.size(); ++i){
        v1[i] -= v2[i];
    }
    return v1;
}

std::vector<double> seidel_method(std::vector<std::vector<double>> H_matrix, std::vector<double> g_vector, int iteration){
    int n = g_vector.size();
    std::vector<double> x_vector_current(n, 0);     //current vector is for the next iteration(k+1)
    for(int k = 0; k < iteration; ++k){
        std::vector<double> x_vector_temp = x_vector_current;       //temp is for current iteration(k)
        for(int i = 0; i < n; ++i) {
            double sum_first = 0, sum_second = 0;   //sum_first up to i - 1 el, meanwhile sum_second from i to n element
            for(int j = 0; j < i - 1; ++j){
                sum_first += H_matrix[i][j]*x_vector_temp[j];
            }
            for (int j = i; j < n; ++j) {
                sum_second += H_matrix[i][j] * x_vector_current[j];
            }
            x_vector_temp[i] = sum_first + sum_second + g_vector[i];
        }
        x_vector_current = x_vector_temp;
    }
    return x_vector_current;
}


std::vector<double> top_relaxation_method(std::vector<std::vector<double>> H_matrix, std::vector<double> g_vector, int iteration, double q){
    int n = g_vector.size();
    std::vector<double> x_vector_current(n, 0);     //current vector is for the previous iteration(k-1)
    for(int k = 0; k < iteration; ++k){
        std::vector<double> x_vector_temp = x_vector_current;       //temp is for current iteration(k)
        std::vector<double> x_vector_temp_new(n);
        for(int i = 0; i < n; ++i) {
            double sum_first = 0, sum_second = 0;   //sum_first up to i - 1 el, meanwhile sum_second from i to n element
            for(int j = 0; j < i - 1; ++j){
                sum_first += H_matrix[i][j]*x_vector_temp[j];
            }
            for (int j = i + 1; j < n; ++j) {
                sum_second += H_matrix[i][j] * x_vector_current[j];
            }
            x_vector_temp_new[i] = x_vector_current[i] + q*(sum_first + sum_second - x_vector_current[i] + g_vector[i]);
        }
        x_vector_current = x_vector_temp_new;
    }
    return x_vector_current;
}

int main() {

    std::ifstream file("C:\\Games\\iteration method\\data.txt");
    if(!file){
        std::cout << "File not found";
    }
    else
        std::cout << "File has opened successfully\n";

    int n;
    std::cout << "Enter number of rows and cols: " << std::endl;
    std::cin >> n;
    std::vector<std::vector<double>> arr(n, std::vector<double> (n));
    std::cout << "Enter A matrix values: " << std::endl;
    std::string string;
    //std::getline(file, string);     //first getline for reading /n symbol, which has been left by std::cout
    std::vector<double> values_temp;
    for(int i = 0; i < n; ++i){
        std::getline(file, string);
        values_temp = split_string_into_array(string);
        for(int j = 0; j < n; ++j){
            arr[i][j] = values_temp[j];
        }
    }
    std::vector<double> b_vector;
    std::cout << "Enter b-vector:" << std::endl;
    getline(file, string);
    b_vector = split_string_into_array(string);
    std::vector <double> accurate_answer = Gauss(b_vector, arr, n, n);
    transforming(arr, b_vector);            //"A" matrix is "H" now, meanwhile "b_vector" is "g" vector now
    std::cout << "Infinity norm for H matrix is: " << infinity_norm_of_matrix(arr) << std::endl;
    std::cout << "Accurate answer is: ";        //Gauss answer
    for(auto el : accurate_answer)
        std::cout << el << " ";
    std::cout << std::endl;
    std::vector<double> approximate_answer = simple_iteration_method(arr, b_vector, 7);

    std::cout << "A priori error estimate of ||x^(7) - x*||_{infinity} "
                 "= ||H||^k||x^(0)|| + ||H||^k/(1 - ||H||)*||g|| = " << pow(infinity_norm_of_matrix(arr), 7)/
                 (1 - infinity_norm_of_matrix(arr)) * infinity_norm_of_vector(b_vector) << std::endl;

    std::cout << "Approximate value of x^(7) is: ";
    for(auto el : approximate_answer)
        std::cout << el << " ";
    std::cout << std::endl;


    std::cout << "Actual error is: " << infinity_norm_of_vector(approximate_answer - accurate_answer) << std::endl;
    std::cout << "A posterior error estimate of ||x^(7) - x*||_{infinity}: "
                << infinity_norm_of_matrix(arr)/(1 - infinity_norm_of_matrix(arr))
                *infinity_norm_of_vector(approximate_answer - simple_iteration_method(arr, b_vector, 6));


    std::vector<double> approximate_answer_of_two_hundred = simple_iteration_method(arr, b_vector, 200);
    std::cout << "\nApproximate value of x^(200) is ";
    for(auto el : approximate_answer_of_two_hundred)
        std::cout << el << " ";
    std::cout << std::endl;


    std::vector<double> seidel_method_answer = seidel_method(arr, b_vector, 7);
    std::cout << "Seidel's method result for x^(7) is: ";
    for(auto el: seidel_method_answer)
        std::cout << el << " ";
    std::cout << std::endl;


    //top relaxation method
    double spectral_radius = infinity_norm_of_vector(simple_iteration_method(arr, b_vector, 10)-
            simple_iteration_method(arr, b_vector, 9))/
                    infinity_norm_of_vector(simple_iteration_method(arr, b_vector, 8) -
                    simple_iteration_method(arr, b_vector, 7));
    std::cout << "Spectral radius is: " << spectral_radius << std::endl;
    double q_opt = 2/(1 + sqrt(1 - pow(spectral_radius, 2))), q_1 = q_opt - 0.1, q_2 = q_opt + 0.1;
    std::vector<double> answer_top_relax_first = top_relaxation_method(arr, b_vector, 7, q_opt),
        answer_top_relax_second = top_relaxation_method(arr, b_vector, 7, q_1),
        answer_top_relax_third = top_relaxation_method(arr, b_vector, 7, q_2);
    for(int i = 0; i < n; ++i){
        std::cout << "Top relaxation method, result with q_opt: ";
        for(int j = 0; j < answer_top_relax_first.size(); ++j){
            std::cout << answer_top_relax_first[j] << " ";
        }
        std::cout << std::endl << "Top relaxation method, result with q_1: ";
        for(int j = 0; j < answer_top_relax_second.size(); ++j){
            std::cout << answer_top_relax_second[j] << " ";
        }
        std::cout << std::endl << "Top relaxation method, result with q_2: ";
        for(int j = 0; j < answer_top_relax_third.size(); ++j){
            std::cout << answer_top_relax_third[j] << " ";
        }
        std::cout << std::endl;

    }
    file.close();
}