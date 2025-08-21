#pragma once
#include "moses.h"
#include "firm.h"
#include "market.h"

namespace moses {

// Sum across a vector
inline double Moses::sum(std::vector<double>& vec) const {
	return std::accumulate(vec.begin(), vec.end(), 0.0);
}

// Sum across an array
template<std::size_t N>
inline double Moses::sum(const std::array<double, N>& ar) const {
	    return std::accumulate(ar.begin(), ar.end(), 0.0);
}

// Sum across a vector per market
inline double Moses::sum(const std::vector<double>& vec, const int market_id) const {
	if (vec.size() != firms.size()) {
		throw std::runtime_error("Vector doesn't match in sum()");
	}
	
	double sum = 0.0;
	for (size_t i = 0; i < vec.size(); i++) {
		if (firms[i]->market_id == market_id) {
			sum += vec[i];
		}
	}
	return sum;
}

// Sum one variable across firms
template<typename T>
inline double Moses::sum(T Firm::*member) const {
    double sum = 0.0;
    for (const auto& firm : firms) {
        sum += firm.get()->*member;
    }
    return sum;
}	

// Sum one variable across markets
template<typename T>
inline double Moses::sum(T Market::*member) const {
    double sum = 0.0;
    for (const auto& mkt : markets) {
        sum += mkt.get()->*member;
    }
    return sum;
}

// Sum one variable across firms in a sector
template<typename T>
inline double Moses::sum(T Firm::*member, const int market_id) const {
    double sum = 0.0;
    for (const auto& firm : firms) {
		if (market_id == firm->market_id) {
			sum += firm.get()->*member;
		}
    }
    return sum;
}

// Weighted average sum across firms
template<typename WeightT, typename ValueT>
inline double Moses::sum_avg(WeightT Firm::*weight, ValueT Firm::*value) const {
	double sum_products = 0.0;
	double sum_weights = 0.0;

	for (const auto& firm : firms) {
		sum_products += (firm.get()->*weight) * (firm.get()->*value);
		sum_weights += firm.get()->*weight;
	}

	if (sum_weights == 0.0) {
		throw std::runtime_error("Division by 0 in sum_avg()");
	}

	return sum_products / sum_weights;
}

// Weighted average sum across firms per market
template<typename WeightT, typename ValueT>
inline double Moses::sum_avg(WeightT Firm::*weight, ValueT Firm::*value, const int market_id) const {
	double sum_products = 0.0;
	double sum_weights = 0.0;

	for (const auto& firm : firms) {
		if (firm->market_id == market_id) {
			sum_products += (firm.get()->*weight) * (firm.get()->*value);
			sum_weights += firm.get()->*weight;
		}
	}

	if (sum_weights == 0.0) {
		throw std::runtime_error("Division by 0 in sum_avg()");
	}

	return sum_products / sum_weights;
}
// Sum the multiplication across two vectors
inline double Moses::sum_mult(const std::vector<double>& v1, const std::vector<double>& v2) const {
	if (v1.size() != v2.size()) {
		throw std::runtime_error("Vectors do not match in sum_mult()");
	}			

	double sum = 0.0;
	for (size_t i = 0; i < v1.size(); i++) {
		sum += v1[i] * v2[i];
	}
	return sum;
}

// Sum the multiplication between two variables across firms
template<typename T1, typename T2>
inline double Moses::sum_mult(T1 Firm::*member1, T2 Firm::*member2, const int market_id) const {
    double sum = 0.0;
    for (const auto& firm : firms) {
		if (market_id == firm->market_id) {
			sum += firm.get()->*member1 * firm.get()->*member2;
		}
    }
    return sum;
}

// Sum the multiplication across an array and a vector
template<std::size_t N>
inline double Moses::sum_mult(const std::array<double, N>& ar, const std::vector<double> vec) const {
	if (ar.size() != vec.size()) {
		throw std::runtime_error("Vectors do not match in sum_mult()");
	}			

	double sum = 0.0;
	for (size_t i = 0; i < vec.size(); i++) {
		sum += ar[i] * vec[i];
	}
	return sum;
}

// Sum the multiplication between a variable across firms per market, and a specific variable
template<typename T1>
inline double Moses::sum_mult(T1 Firm::*member, const double var, const int market_id) const {
	double sum = 0.0;
	for (const auto& firm : firms) {
		if (market_id == firm->market_id) {
			sum += firm.get()->*member * var;
		}
	}
	return sum;
}

// Sum the ratio between two variables across an array and a vector
template<std::size_t N>
inline double Moses::sum_div(const std::array<double, N>& ar, const std::vector<double>& vec) const {
	if (ar.size() != vec.size()) {
		throw std::runtime_error("Vectors do not match in sum_mult()");
	}	

	double sum = 0.0;
	for (size_t i = 0; i < vec.size(); i++) {
		if (vec[i] == 0.0) {
			throw std::runtime_error("Division by 0 in sum_div()");
		}
		sum += ar[i] / vec[i];
	}
	return sum;
}

// Sum the ratio between two variables across two arrays
template<std::size_t N1, std::size_t N2>
inline double Moses::sum_div(const std::array<double, N1>& ar1, const std::array<double, N2>& ar2) const {
	if (ar1.size() != ar2.size()) {
		throw std::runtime_error("Vectors do not match in sum_mult()");
	}	

	double sum = 0.0;
	for (size_t i = 0; i < ar1.size(); i++) {
		if (ar2[i] == 0.0) {
			throw std::runtime_error("Division by 0 in sum_div()");
		}
		sum += ar1[i] / ar2[i];
	}
	return sum;
}

} // namespace moses
