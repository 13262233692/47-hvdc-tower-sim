#include "common/SparseMatrix.hpp"
#include <algorithm>
#include <map>

namespace hvdc {

SparseMatrix SparseMatrix::from_triplets(Index nrows, Index ncols, const std::vector<Triplet>& triplets) {
    SparseMatrix mat(nrows, ncols);
    
    struct RowCol {
        Index row, col;
        bool operator<(const RowCol& other) const {
            if (row != other.row) return row < other.row;
            return col < other.col;
        }
    };
    
    std::map<RowCol, Real> aggregated;
    for (const auto& t : triplets) {
        aggregated[{std::get<0>(t), std::get<1>(t)}] += std::get<2>(t);
    }
    
    mat.row_ptr_.resize(nrows + 1, 0);
    mat.col_idx_.reserve(aggregated.size());
    mat.values_.reserve(aggregated.size());
    
    for (const auto& [key, val] : aggregated) {
        mat.col_idx_.push_back(key.col);
        mat.values_.push_back(val);
        mat.row_ptr_[key.row + 1] = static_cast<Index>(mat.values_.size());
    }
    
    for (Index i = 1; i <= nrows; ++i) {
        if (mat.row_ptr_[i] == 0) {
            mat.row_ptr_[i] = mat.row_ptr_[i - 1];
        }
    }
    
    mat.dirty_ = false;
    return mat;
}

void SparseMatrix::finalize() {
    if (!dirty_) return;
    
    struct RC {
        Index row, col;
        Real val;
        bool operator<(const RC& other) const {
            if (row != other.row) return row < other.row;
            return col < other.col;
        }
    };
    
    std::vector<RC> entries;
    
    for (const auto& t : triplets_backup_) {
        entries.push_back({std::get<0>(t), std::get<1>(t), std::get<2>(t)});
    }
    
    for (const auto& [key, val] : triplets_set_) {
        Index r = static_cast<Index>(key >> 32);
        Index c = static_cast<Index>(key & 0xFFFFFFFFULL);
        entries.push_back({r, c, val});
    }
    
    std::sort(entries.begin(), entries.end());
    
    col_idx_.clear();
    values_.clear();
    row_ptr_.assign(nrows_ + 1, 0);
    
    Index cur_row = -1;
    Index cur_col = -1;
    
    for (const auto& e : entries) {
        while (cur_row < e.row) {
            cur_row++;
            row_ptr_[cur_row + 1] = static_cast<Index>(col_idx_.size());
            cur_col = -1;
        }
        
        if (e.col == cur_col && !col_idx_.empty()) {
            values_.back() += e.val;
        } else {
            col_idx_.push_back(e.col);
            values_.push_back(e.val);
            cur_col = e.col;
        }
    }
    
    while (cur_row < nrows_) {
        cur_row++;
        row_ptr_[cur_row] = static_cast<Index>(col_idx_.size());
    }
    row_ptr_[nrows_] = static_cast<Index>(col_idx_.size());
    
    triplets_backup_.clear();
    triplets_set_.clear();
    dirty_ = false;
}

void SparseMatrix::multiply(const Vector& x, Vector& y) const {
    if (ncols_ != x.size()) {
        throw std::runtime_error("SparseMatrix multiply: x size mismatch");
    }
    if (y.size() != nrows_) {
        y.resize(nrows_);
    }
    y.zero();
    
    const Index* rp = row_ptr_.data();
    const Index* ci = col_idx_.data();
    const Real* v = values_.data();
    const Real* xp = x.raw().data();
    Real* yp = y.raw().data();
    
    for (Index r = 0; r < nrows_; ++r) {
        Real sum = 0.0;
        for (Index i = rp[r]; i < rp[r + 1]; ++i) {
            sum += v[i] * xp[ci[i]];
        }
        yp[r] = sum;
    }
}

} // namespace hvdc
