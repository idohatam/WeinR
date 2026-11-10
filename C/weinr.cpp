#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List sanger_to_phred(CharacterVector r_character_vector) {
    std::vector<std::vector<int>> result;
    std::vector<std::string> c_string_vector;
    c_string_vector = Rcpp::as<std::vector<std::string>>(r_character_vector);
    std::vector<int> c_integer_vector;
    for(int i = 0; i < int(c_string_vector.size()); i++) {
        for(int j = 0; j < int(c_string_vector[i].size()); j++) {
            c_integer_vector.push_back(int(c_string_vector[i][j] - 33));
        }
        result.push_back(c_integer_vector);
        c_integer_vector.clear();
    }
    return Rcpp::wrap(result);
}

// [[Rcpp::export]]
List sanger_to_probability(CharacterVector r_character_vector) {
    int phred;
    double p;
    double lookup_table[94];
    for(int i = 0; i < 94; i ++) {
        lookup_table[i] = pow(10, -i / 10.0);
    }
    std::vector<std::vector<double>> result;
    std::vector<std::string> c_string_vector;
    c_string_vector = Rcpp::as<std::vector<std::string>>(r_character_vector);
    std::vector<double> c_double_vector;
    for(int i = 0; i < int(c_string_vector.size()); i++) {
        for(int j = 0; j < int(c_string_vector[i].size()); j++) {
        phred = int(c_string_vector[i][j] - 33);
        //p = pow(10, -phred / 10.0);
        p = lookup_table[phred];
        c_double_vector.push_back(p);
        }
        result.push_back(c_double_vector);
        c_double_vector.clear();
    }
    return Rcpp::wrap(result);
}

// [[Rcpp::export]]
NumericVector calculate_per_read_avg_q_score(List r_p_vectors) {
    std::vector<double> result;
    std::vector<std::vector<double>> c_p_vectors;
    c_p_vectors = Rcpp::as<std::vector<std::vector<double>>>(r_p_vectors);
    
    for(int i = 0; i < int(c_p_vectors.size()); i++) {
        int count = 0;
        double sum = 0;
        for(int j = 0; j < int(c_p_vectors[i].size()); j++) {
            count = count + 1;
            sum = sum + c_p_vectors[i][j];
        }
        result.push_back(-10 * log10(sum / count));
    }
    
    return Rcpp::wrap(result);
    
}

// [[Rcpp::export]]
List tabulate_qpp_mean(List r_p_vectors, const int MAX_WIDTH) {
    std::vector<std::vector<double>> result;
    std::vector<double> qpp_count(MAX_WIDTH);
    std::vector<double> qpp_sum(MAX_WIDTH);
    std::vector<std::vector<double>> c_p_vectors;
    c_p_vectors = Rcpp::as<std::vector<std::vector<double>>>(r_p_vectors);
    for(int i = 0; i < int(c_p_vectors.size()); i++) {
        for(int j = 0; j < int(c_p_vectors[i].size()); j++) {
        qpp_count[j] += 1;
        qpp_sum[j] = qpp_sum[j] + c_p_vectors[i][j];
        }
    }
    result.push_back(qpp_count);
    result.push_back(qpp_sum);
    return Rcpp::wrap(result);
}

// [[Rcpp::export]]
IntegerMatrix tabulate_qpp_median(IntegerMatrix frequency_table, CharacterVector r_character_vector, const int BIN_WIDTH) {
    std::vector<std::string> c_string_vector;
    c_string_vector = Rcpp::as<std::vector<std::string>>(r_character_vector);
    for(int i = 0; i < int(c_string_vector.size()); i++) {
        int bin;
        for(int j = 0; j < int(c_string_vector[i].size()); j++) {
            bin = int(int(c_string_vector[i][j] - 33) / BIN_WIDTH);
            frequency_table(bin, j) = frequency_table(bin, j) + 1;
        }
    }
    return frequency_table;
}

// [[Rcpp::export]]
NumericVector calculate_qpp_mean(IntegerVector r_qpp_count, NumericVector r_qpp_sum, const int MAX_WIDTH) {
    std::vector<int> c_qpp_count(MAX_WIDTH);
    std::vector<double> c_qpp_sum(MAX_WIDTH);
    std::vector<double> result(MAX_WIDTH);
    c_qpp_count = Rcpp::as<std::vector<int>>(r_qpp_count);
    c_qpp_sum = Rcpp::as<std::vector<double>>(r_qpp_sum);
    for(int i = 0; i < MAX_WIDTH; i++) {
        if(c_qpp_count[i] == 0) {
            continue;
        }
        result[i] =  -10 * log10(c_qpp_sum[i] / c_qpp_count[i]);
    }
    return Rcpp::wrap(result);
}

// [[Rcpp::export]]
NumericVector calculate_qpp_median(IntegerMatrix frequency_table, const int MAX_WIDTH, const int BINS, const int BIN_WIDTH) {
    std::vector<double> result(MAX_WIDTH);
    for(int j = 0; j < MAX_WIDTH; j++) {
        int total_frequency = 0;
        int target_frequency = 0;
        int cumulative_frequency = 0;
        for(int i = 0; i < BINS; i++) {
            total_frequency = total_frequency + frequency_table(i, j);
        }
        if(total_frequency == 0) {
            continue;
        }
        target_frequency = total_frequency / 2;
        for(int i = 0; i < BINS; i++) {
            cumulative_frequency = cumulative_frequency + frequency_table(i, j);
            if(cumulative_frequency > target_frequency) {
                int lower_bound = i * BIN_WIDTH;
                int upper_bound = i * BIN_WIDTH + BIN_WIDTH;
                double median = (lower_bound + upper_bound) / 2.0; // mathematically incorrect but good enough
                result[j] = median;
                break;
            }
        }
    }
    return Rcpp::wrap(result);
}