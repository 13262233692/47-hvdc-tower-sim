#ifndef HVDC_COMMON_SPARSE_MATRIX_HPP
#define HVDC_COMMON_SPARSE_MATRIX_HPP

#include "common/Types.hpp"
#include "common/Vector.hpp"
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <stdexcept>

namespace hvdc {

class SparseMatrix {
public:
    using Triplet = std::tuple<Index, Index, Real>;
    
    SparseMatrix() = default;
    
    SparseMatrix(Index nrows, Index ncols)
        : nrows_(nrows), ncols_(ncols) {
        row_ptr_.resize(nrows + 1, 0);
    }
    
    SparseMatrix(Index nrows, Index ncols, Index nnz_estimate)
        : nrows_(nrows), ncols_(ncols) {
        row_ptr_.resize(nrows + 1, 0);
        col_idx_.reserve(nnz_estimate);
        values_.reserve(nnz_estimate);
    }
    
    static SparseMatrix from_triplets(Index nrows, Index ncols, const std::vector<Triplet>& triplets);
    
    Index rows() const { return nrows_; }
    Index cols() const { return ncols_; }
    Index nnz() const { return static_cast<Index>(values_.size()); }
    Index num_nonzeros() const { return nnz(); }
    
    void reserve(Index nnz) {
        col_idx_.reserve(nnz);
        values_.reserve(nnz);
    }
    
    void insert(Index row, Index col, Real value) {
        if (row < 0 || row >= nrows_ || col < 0 || col >= ncols_) {
            throw std::out_of_range("SparseMatrix insert index out of range");
        }
        triplets_backup_.emplace_back(row, col, value);
        dirty_ = true;
    }
    
    void add(Index row, Index col, Real value) {
        triplets_backup_.emplace_back(row, col, value);
        dirty_ = true;
    }
    
    void set(Index row, Index col, Real value) {
        UInt64 key = (static_cast<UInt64>(static_cast<UInt32>(row)) << 32) | static_cast<UInt64>(static_cast<UInt32>(col));
        triplets_set_[key] = value;
        dirty_ = true;
    }
    
    void finalize();
    
    Real operator()(Index row, Index col) const {
        if (dirty_) const_cast<SparseMatrix*>(this)->finalize();
        Index start = row_ptr_[row];
        Index end = row_ptr_[row + 1];
        for (Index i = start; i < end; ++i) {
            if (col_idx_[i] == col) return values_[i];
        }
        return 0.0;
    }
    
    void multiply(const Vector& x, Vector& y) const;
    
    Vector multiply(const Vector& x) const {
        Vector y(nrows_);
        multiply(x, y);
        return y;
    }
    
    void diagonal(Vector& diag) const {
        diag.resize(nrows_);
        diag.zero();
        for (Index i = 0; i < nrows_; ++i) {
            diag[i] = (*this)(i, i);
        }
    }
    
    Vector diagonal() const {
        Vector d(nrows_);
        diagonal(d);
        return d;
    }
    
    void scale_rows(const Vector& scale_factors) {
        if (dirty_) finalize();
        for (Index r = 0; r < nrows_; ++r) {
            Real s = scale_factors[r];
            Index start = row_ptr_[r];
            Index end = row_ptr_[r + 1];
            for (Index i = start; i < end; ++i) values_[i] *= s;
        }
    }
    
    void clear() {
        nrows_ = ncols_ = 0;
        row_ptr_.clear();
        col_idx_.clear();
        values_.clear();
        triplets_backup_.clear();
        triplets_set_.clear();
        dirty_ = false;
    }
    
    void reset(Index nrows, Index ncols) {
        nrows_ = nrows;
        ncols_ = ncols;
        row_ptr_.assign(nrows + 1, 0);
        col_idx_.clear();
        values_.clear();
        triplets_backup_.clear();
        triplets_set_.clear();
        dirty_ = false;
    }
    
    const IndexVec& row_ptr() const {
        if (dirty_) const_cast<SparseMatrix*>(this)->finalize();
        return row_ptr_;
    }
    
    const IndexVec& col_idx() const {
        if (dirty_) const_cast<SparseMatrix*>(this)->finalize();
        return col_idx_;
    }
    
    const Vec& values() const {
        if (dirty_) const_cast<SparseMatrix*>(this)->finalize();
        return values_;
    }
    
    IndexVec& row_ptr() {
        if (dirty_) finalize();
        return row_ptr_;
    }
    
    IndexVec& col_idx() {
        if (dirty_) finalize();
        return col_idx_;
    }
    
    Vec& values() {
        if (dirty_) finalize();
        return values_;
    }

private:
    Index nrows_ = 0;
    Index ncols_ = 0;
    IndexVec row_ptr_;
    IndexVec col_idx_;
    Vec values_;
    
    std::vector<Triplet> triplets_backup_;
    std::unordered_map<UInt64, Real> triplets_set_;
    bool dirty_ = false;
};

} // namespace hvdc

#endif // HVDC_COMMON_SPARSE_MATRIX_HPP
