#ifndef HVDC_COMMON_VECTOR_HPP
#define HVDC_COMMON_VECTOR_HPP

#include "common/Types.hpp"
#include <vector>
#include <numeric>
#include <cmath>
#include <stdexcept>

namespace hvdc {

class Vector {
public:
    Vector() = default;
    
    explicit Vector(Index size, Real value = 0.0)
        : data_(size, value), size_(size) {}
    
    explicit Vector(const Vec& data)
        : data_(data), size_(static_cast<Index>(data.size())) {}
    
    Index size() const { return size_; }
    
    Real& operator[](Index i) { return data_[i]; }
    const Real& operator[](Index i) const { return data_[i]; }
    
    Real& at(Index i) {
        if (i < 0 || i >= size_) throw std::out_of_range("Vector index out of range");
        return data_[i];
    }
    const Real& at(Index i) const {
        if (i < 0 || i >= size_) throw std::out_of_range("Vector index out of range");
        return data_[i];
    }
    
    Vec& raw() { return data_; }
    const Vec& raw() const { return data_; }
    
    void resize(Index size, Real value = 0.0) {
        data_.resize(size, value);
        size_ = size;
    }
    
    void assign(Index size, Real value = 0.0) {
        data_.assign(size, value);
        size_ = size;
    }
    
    void fill(Real value) {
        std::fill(data_.begin(), data_.end(), value);
    }
    
    void zero() { fill(0.0); }
    
    Real norm() const {
        Real sum = 0.0;
        for (Index i = 0; i < size_; ++i) sum += data_[i] * data_[i];
        return std::sqrt(sum);
    }
    
    Real norm_inf() const {
        Real max_val = 0.0;
        for (Index i = 0; i < size_; ++i) {
            Real abs_v = std::fabs(data_[i]);
            if (abs_v > max_val) max_val = abs_v;
        }
        return max_val;
    }
    
    Real sum() const {
        Real s = 0.0;
        for (Index i = 0; i < size_; ++i) s += data_[i];
        return s;
    }
    
    Real dot(const Vector& other) const {
        if (size_ != other.size_) throw std::runtime_error("Vector size mismatch in dot");
        Real sum = 0.0;
        for (Index i = 0; i < size_; ++i) sum += data_[i] * other.data_[i];
        return sum;
    }
    
    void axpy(Real alpha, const Vector& x) {
        if (size_ != x.size_) throw std::runtime_error("Vector size mismatch in axpy");
        for (Index i = 0; i < size_; ++i) data_[i] += alpha * x.data_[i];
    }
    
    void scale(Real alpha) {
        for (Index i = 0; i < size_; ++i) data_[i] *= alpha;
    }
    
    Vector& operator+=(const Vector& rhs) {
        if (size_ != rhs.size_) throw std::runtime_error("Vector size mismatch in +=");
        for (Index i = 0; i < size_; ++i) data_[i] += rhs.data_[i];
        return *this;
    }
    
    Vector& operator-=(const Vector& rhs) {
        if (size_ != rhs.size_) throw std::runtime_error("Vector size mismatch in -=");
        for (Index i = 0; i < size_; ++i) data_[i] -= rhs.data_[i];
        return *this;
    }
    
    Vector& operator*=(Real alpha) {
        scale(alpha);
        return *this;
    }
    
    Vector operator+(const Vector& rhs) const {
        Vector result(*this);
        result += rhs;
        return result;
    }
    
    Vector operator-(const Vector& rhs) const {
        Vector result(*this);
        result -= rhs;
        return result;
    }
    
    Vector operator*(Real alpha) const {
        Vector result(*this);
        result *= alpha;
        return result;
    }
    
    void set(Index idx, Real val) { data_[idx] = val; }
    void add(Index idx, Real val) { data_[idx] += val; }
    Real get(Index idx) const { return data_[idx]; }
    
    void copy_from(const Vector& src) {
        if (size_ != src.size_) resize(src.size_);
        data_ = src.data_;
    }
    
    void swap(Vector& other) noexcept {
        std::swap(data_, other.data_);
        std::swap(size_, other.size_);
    }

private:
    Vec data_;
    Index size_ = 0;
};

inline Vector operator*(Real alpha, const Vector& v) {
    return v * alpha;
}

} // namespace hvdc

#endif // HVDC_COMMON_VECTOR_HPP
