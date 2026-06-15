#ifndef HVDC_COMMON_PETSC_MANAGER_HPP
#define HVDC_COMMON_PETSC_MANAGER_HPP

#include "common/Types.hpp"
#include "common/SparseMatrix.hpp"
#include "common/Vector.hpp"
#include <string>
#include <memory>
#include <stdexcept>

#ifdef HVDC_USE_PETSC
#include "petsc.h"
#include "petscksp.h"
#endif

namespace hvdc {

class PETScMatrix;
class PETScVector;
class PETScSolver;

class PETScManager {
public:
    static PETScManager& instance();
    
    void initialize(int* argc, char*** argv);
    void finalize();
    
    bool is_initialized() const { return initialized_; }
    bool has_petsc() const {
#ifdef HVDC_USE_PETSC
        return true;
#else
        return false;
#endif
    }
    
    void set_options(const std::string& options_file);
    void set_option(const std::string& name, const std::string& value);
    
    std::string get_version() const;

private:
    PETScManager() = default;
    ~PETScManager();
    
    PETScManager(const PETScManager&) = delete;
    PETScManager& operator=(const PETScManager&) = delete;
    
    bool initialized_ = false;
};

class PETScVector {
public:
    PETScVector() = default;
    explicit PETScVector(Index size);
    PETScVector(Index local_size, Index global_size);
    ~PETScVector();
    
    PETScVector(const PETScVector&) = delete;
    PETScVector& operator=(const PETScVector&) = delete;
    PETScVector(PETScVector&& other) noexcept;
    PETScVector& operator=(PETScVector&& other) noexcept;
    
    void create(Index size);
    void create(Index local_size, Index global_size);
    void destroy();
    
    void from_vector(const Vector& vec);
    void to_vector(Vector& vec) const;
    
    void set(Index idx, Real val);
    void add(Index idx, Real val);
    Real get(Index idx) const;
    
    void set_all(Real val);
    void zero();
    
    void assembly_begin();
    void assembly_end();
    
    Index size() const { return size_; }
    Index local_size() const { return local_size_; }
    
    Real norm() const;
    Real dot(const PETScVector& other) const;
    void axpy(Real alpha, const PETScVector& x);
    void scale(Real alpha);
    void copy_from(const PETScVector& src);
    
    void swap(PETScVector& other) noexcept;
    
#ifdef HVDC_USE_PETSC
    Vec& petsc_vec() { return vec_; }
    const Vec& petsc_vec() const { return vec_; }
#endif

private:
    Index size_ = 0;
    Index local_size_ = 0;
#ifdef HVDC_USE_PETSC
    Vec vec_ = nullptr;
#endif
};

class PETScMatrix {
public:
    PETScMatrix() = default;
    PETScMatrix(Index nrows, Index ncols, Index nnz_per_row = 10);
    PETScMatrix(Index local_rows, Index global_rows, Index global_cols, Index nnz_per_row = 10);
    ~PETScMatrix();
    
    PETScMatrix(const PETScMatrix&) = delete;
    PETScMatrix& operator=(const PETScMatrix&) = delete;
    PETScMatrix(PETScMatrix&& other) noexcept;
    PETScMatrix& operator=(PETScMatrix&& other) noexcept;
    
    void create(Index nrows, Index ncols, Index nnz_per_row = 10);
    void create(Index local_rows, Index global_rows, Index global_cols, Index nnz_per_row = 10);
    void destroy();
    
    void from_sparse_matrix(const SparseMatrix& mat);
    void to_sparse_matrix(SparseMatrix& mat) const;
    
    void set(Index row, Index col, Real val);
    void add(Index row, Index col, Real val);
    
    void set_row(Index row, const IndexVec& cols, const Vec& vals);
    void add_row(Index row, const IndexVec& cols, const Vec& vals);
    
    void zero();
    void zero_rows(const IndexVec& rows, Real diag_value = 0.0);
    void shift_diagonal(Real alpha);
    
    void assembly_begin();
    void assembly_end();
    
    void multiply(const PETScVector& x, PETScVector& y) const;
    
    Index rows() const { return nrows_; }
    Index cols() const { return ncols_; }
    Index local_rows() const { return local_rows_; }
    Index nnz() const;
    
    void set_preallocation_nnz(const IndexVec& d_nnz, const IndexVec& o_nnz);
    
#ifdef HVDC_USE_PETSC
    Mat& petsc_mat() { return mat_; }
    const Mat& petsc_mat() const { return mat_; }
#endif

private:
    Index nrows_ = 0;
    Index ncols_ = 0;
    Index local_rows_ = 0;
#ifdef HVDC_USE_PETSC
    Mat mat_ = nullptr;
#endif
};

class PETScSolver {
public:
    PETScSolver() = default;
    ~PETScSolver();
    
    PETScSolver(const PETScSolver&) = delete;
    PETScSolver& operator=(const PETScSolver&) = delete;
    
    void create();
    void destroy();
    
    void set_operator(const PETScMatrix& A);
    void set_preconditioner(const std::string& type = "bjacobi");
    void set_ksp_type(const std::string& type = "gmres");
    
    void set_tolerances(Real rtol = 1e-5, Real atol = 1e-50, 
                        Real dtol = 1e5, Index max_iters = 10000);
    
    bool solve(const PETScVector& b, PETScVector& x);
    
    Index iterations() const { return iterations_; }
    Real residual_norm() const { return residual_; }
    bool converged() const { return converged_; }
    std::string converged_reason() const { return reason_; }
    
    void set_from_options();
    
#ifdef HVDC_USE_PETSC
    KSP& petsc_ksp() { return ksp_; }
    PC& petsc_pc() { return pc_; }
#endif

private:
    bool created_ = false;
    Index iterations_ = 0;
    Real residual_ = 0.0;
    bool converged_ = false;
    std::string reason_;
#ifdef HVDC_USE_PETSC
    KSP ksp_ = nullptr;
    PC pc_ = nullptr;
#endif
};

} // namespace hvdc

#endif // HVDC_COMMON_PETSC_MANAGER_HPP
