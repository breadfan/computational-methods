#include <iostream>
#include <vector>
#include <cmath>

/* Program implements simple Gauss method for solving system of algebraic equations
 */

std::vector<double> Gauss(std::vector<double>, std::vector<std::vector<double>>, int, int);
std::vector<double> split_string_into_array(std::string);
void swap_vector_rows(std::vector< std::vector <double> >, int, int);
std::vector<std::vector<double> > reverse_matrix(std::vector< std::vector<double>> , int);
void residual_counting(std::vector< std::vector<double>> &, std::vector<double>&, std::vector<double>&);


std::vector< std::vector <double> > reverse_matrix(std::vector< std::vector < double > > arr, int size_vector){
    std::vector< std::vector<double> > result;
    std::vector<double> unity_vector;
    for(int i = 0; i < size_vector; ++i){           //creating unity b-vectors
        if(!i){
            for (int j = 0; j < size_vector; ++j){
                if (i == j){
                    unity_vector.push_back(1);
                }
                else{
                    unity_vector.push_back(0);
                }
            }
        }
        else{
            unity_vector[i] = 1;
            unity_vector[i - 1] = 0;
        }
        result.push_back(Gauss(unity_vector, arr, size_vector, size_vector));
    }
    return result;
}

std::vector<double> split_string_into_array(std::string str){
    std::vector<double> arr;
    while(str.length() != 0){
        int ind_of_space = str.find(' ');
        arr.push_back(stod(str.substr(0, ind_of_space)));
        str = str.substr(ind_of_space + 1, str.length() - ind_of_space);
        if(str.length() == 1){
            if(str[0] != ' ')
                arr.push_back(stod(str));
            break;
        }
        else if(str.length() == 2 && str[0] == '-'){
            arr.push_back(stod(str));
            break;
        }
    }
    return arr;
}

void swap_vector_rows(std::vector<std::vector <double> > arr, int first_ind, int second_ind){
    for(std::size_t i = 0; i < arr[first_ind].size(); ++i)
    {
        int temp = arr[first_ind][i];
        arr[first_ind][i] = arr[second_ind][i];
        arr[second_ind][i] = temp;
    }
}

std::vector<double> Gauss(std::vector<double> answer, std::vector<std::vector<double> > arr, int n, int m){
        //forward step//
    for(int j = 0; j < m; ++j){     //begin on columns
        double max_in_column = fabs(arr[0][j]);
        int ind_for_max_in_col = j;
        for(int row = 0; row < n; ++row){
            if(fabs(arr[row][j]) > max_in_column){
                max_in_column = arr[row][j];
                ind_for_max_in_col = row;
            }
        }
        if(ind_for_max_in_col != j) swap_vector_rows(arr, j, ind_for_max_in_col);
        for (int i = n - 1; i > j; --i){        //starting from bottom board
            double coefficient_after_division = arr[i][j]/arr[j][j];
            for(int k = m - 1; k > j; --k){     //starting from right border
                arr[i][k] -= arr[j][k]*coefficient_after_division;          //that's absolutely right
            }
            answer[i] -= answer[j]*coefficient_after_division;
        }
    }
        //reverse step//
    std::vector<double> result_vector;
    for(int i = n - 1; i > -1; --i){
        int id = m - 1;
        while(id > i)
            answer[i] -= arr[i][id--];
        double cur_el = answer[i]/arr[i][i];
        for(int row = n - 2; row > -1; --row){
            arr[row][i] *= cur_el;
        }
        result_vector.insert(result_vector.begin(), cur_el);
        //result_vector.push_back(cur_el);
    }
    return result_vector;
}

void residual_counting(std::vector<std::vector<double>> &arr, std::vector<double> &result, std::vector<double> &b_vector){
    std::vector<double> residue_vector;
    int n = b_vector.size();
    for(int i = 0; i < n; ++i) {
        double temp_sum = 0;
        for (int j = 0; j < n; ++j) {
            temp_sum += result[j]*arr[i][j];
        }
        residue_vector.push_back(b_vector[i] - temp_sum);
    }
    for(double el: residue_vector)
        std::cout << el << " ";
}


int main()
{
    int n, m;
    std::cout << "Enter number of rows and cols: " << std::endl;
    std::cin >> n >> m;
    std::vector< std::vector<double> > arr(n, std::vector<double> (m));
    std::cout << "Enter A matrix values: " << std::endl;
    std::string string;
    std::getline(std::cin, string);     //first getline for reading /n symbol, which has been left by std::cout
    std::vector<double> values_temp;
    for(int i = 0; i < n; ++i){
        std::getline(std::cin, string);
        values_temp = split_string_into_array(string);
        for(int j = 0; j < m; ++j){
            arr[i][j] = values_temp[j];
        }
    }
    std::vector<double> b_vector;
    std::cout << "Enter b-vector:" << std::endl;
    getline(std::cin, string);
    b_vector = split_string_into_array(string);
    std::vector <double> result = Gauss(b_vector, arr, n, m);
    std::cout << "Your answer is: ";
    for(auto el: result) {
        std::cout << el << " ";
    }
    std::cout << std::endl;
    residual_counting(arr, result, b_vector);
    std::cout << "\nNow we change our matrix: a_{11} *= 10^(-8), so the new residue vector is: ";
    arr[0][0] *= pow(10, -8);
    residual_counting(arr, result, b_vector);
    std::cout << "\nReverse matrix for A-matrix: " << std::endl;
    std::vector<std::vector<double>> _reverse_matrix = reverse_matrix(arr, n);
    for(int j = n - 1; j > -1; --j){
        for(int i = 0; i < n; ++i){
            std::cout << _reverse_matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}