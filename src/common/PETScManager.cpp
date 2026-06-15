#include "common/PETScManager.hpp"
#include <utility>
#include <algorithm>

namespace hvdc {

PETScManager& PETScManager::instance() {
    static PETScManager inst;
    return inst;
}

PETScManager::~PETScManager() {
    if (initialized_) finalize();
}

void PETScManager::initialize(int* argc, char*** argv) {
    if (initialized_) return;
#ifdef HVDC_USE_PETSC
    PetscInitialize(argc, argv, nullptr, nullptr);
#endif
    initialized_ = true;
}

void PETScManager::finalize() {
    if (!initialized_) return;
#ifdef HVDC_USE_PETSC
    PetscFinalize();
#endif
    initialized_ = false;
}

void PETScManager::set_options(const std::string& options_file) {
#ifdef HVDC_USE_PETSC
    if (!options_file.empty()) {
        PetscOptionsInsertFile(PETSC_COMM_WORLD, nullptr, options_file.c_str(), PETSC_FALSE);
    }
#endif
}

void PETScManager::set_option(const std::string& name, const std::string& value) {
#ifdef HVDC_USE_PETSC
    PetscOptionsSetValue(nullptr, name.c_str(), value.c_str());
#endif
}

std::string PETScManager::get_version() const {
#ifdef HVDC_USE_PETSC
    char version[256];
    PetscGetVersion(version, 256);
    return std::string(version);
#else
    return "PETSc not available";
#endif
}

// ===== PETScVector =====

PETScVector::PETScVector(Index size) {
    create(size);
}

PETScVector::PETScVector(Index local_size, Index global_size) {
    create(local_size, global_size);
}

PETScVector::~PETScVector() { destroy(); }

PETScVector::PETScVector(PETScVector&& other) noexcept
    : size_(other.size_), local_size_(other.local_size_)
#ifdef HVDC_USE_PETSC
    , vec_(other.vec_)
#endif
{
    other.size_ = 0;
    other.local_size_ = 0;
#ifdef HVDC_USE_PETSC
    other.vec_ = nullptr;
#endif
}

PETScVector& PETScVector::operator=(PETScVector&& other) noexcept {
    if (this != &other) {
        destroy();
        size_ = other.size_;
        local_size_ = other.local_size_;
#ifdef HVDC_USE_PETSC
        vec_ = other.vec_;
        other.vec_ = nullptr;
#endif
        other.size_ = 0;
        other.local_size_ = 0;
    }
    return *this;
}

void PETScVector::create(Index size) {
    destroy();
    size_ = size;
    local_size_ = PETSC_DECIDE;
#ifdef HVDC_USE_PETSC
    VecCreate(PETSC_COMM_WORLD, &vec_);
    VecSetSizes(vec_, PETSC_DECIDE, size);
    VecSetFromOptions(vec_);
#endif
}

void PETScVector::create(Index local_size, Index global_size) {
    destroy();
    size_ = global_size;
    local_size_ = local_size;
#ifdef HVDC_USE_PETSC
    VecCreate(PETSC_COMM_WORLD, &vec_);
    VecSetSizes(vec_, local_size, global_size);
    VecSetFromOptions(vec_);
#endif
}

void PETScVector::destroy() {
#ifdef HVDC_USE_PETSC
    if (vec_) {
        VecDestroy(&vec_);
        vec_ = nullptr;
    }
#endif
    size_ = 0;
    local_size_ = 0;
}

void PETScVector::from_vector(const Vector& vec) {
    if (size_ != vec.size()) create(vec.size());
#ifdef HVDC_USE_PETSC
    VecSet(vec_, 0.0);
    for (Index i = 0; i < vec.size(); ++i) {
        VecSetValue(vec_, i, vec[i], INSERT_VALUES);
    }
    VecAssemblyBegin(vec_);
    VecAssemblyEnd(vec_);
#endif
}

void PETScVector::to_vector(Vector& vec) const {
    vec.resize(size_);
#ifdef HVDC_USE_PETSC
    for (Index i = 0; i < size_; ++i) {
        PetscScalar val;
        VecGetValues(vec_, 1, &i, &val);
        vec[i] = PetscRealPart(val);
    }
#endif
}

void PETScVector::set(Index idx, Real val) {
#ifdef HVDC_USE_PETSC
    VecSetValue(vec_, idx, val, INSERT_VALUES);
#endif
}

void PETScVector::add(Index idx, Real val) {
#ifdef HVDC_USE_PETSC
    VecSetValue(vec_, idx, val, ADD_VALUES);
#endif
}

Real PETScVector::get(Index idx) const {
#ifdef HVDC_USE_PETSC
    PetscScalar val;
    VecGetValues(vec_, 1, &idx, &val);
    return PetscRealPart(val);
#else
    return 0.0;
#endif
}

void PETScVector::set_all(Real val) {
#ifdef HVDC_USE_PETSC
    VecSet(vec_, val);
#endif
}

void PETScVector::zero() { set_all(0.0); }

void PETScVector::assembly_begin() {
#ifdef HVDC_USE_PETSC
    VecAssemblyBegin(vec_);
#endif
}

void PETScVector::assembly_end() {
#ifdef HVDC_USE_PETSC
    VecAssemblyEnd(vec_);
#endif
}

Real PETScVector::norm() const {
#ifdef HVDC_USE_PETSC
    PetscReal n;
    VecNorm(vec_, NORM_2, &n);
    return n;
#else
    return 0.0;
#endif
}

Real PETScVector::dot(const PETScVector& other) const {
#ifdef HVDC_USE_PETSC
    PetscScalar d;
    VecDot(vec_, other.vec_, &d);
    return PetscRealPart(d);
#else
    return 0.0;
#endif
}

void PETScVector::axpy(Real alpha, const PETScVector& x) {
#ifdef HVDC_USE_PETSC
    VecAXPY(vec_, alpha, x.vec_);
#endif
}

void PETScVector::scale(Real alpha) {
#ifdef HVDC_USE_PETSC
    VecScale(vec_, alpha);
#endif
}

void PETScVector::copy_from(const PETScVector& src) {
#ifdef HVDC_USE_PETSC
    VecCopy(src.vec_, vec_);
#endif
}

void PETScVector::swap(PETScVector& other) noexcept {
    std::swap(size_, other.size_);
    std::swap(local_size_, other.local_size_);
#ifdef HVDC_USE_PETSC
    std::swap(vec_, other.vec_);
#endif
}

// ===== PETScMatrix =====

PETScMatrix::PETScMatrix(Index nrows, Index ncols, Index nnz_per_row) {
    create(nrows, ncols, nnz_per_row);
}

PETScMatrix::PETScMatrix(Index local_rows, Index global_rows, Index global_cols, Index nnz_per_row) {
    create(local_rows, global_rows, global_cols, nnz_per_row);
}

PETScMatrix::~PETScMatrix() { destroy(); }

PETScMatrix::PETScMatrix(PETScMatrix&& other) noexcept
    : nrows_(other.nrows_), ncols_(other.ncols_), local_rows_(other.local_rows_)
#ifdef HVDC_USE_PETSC
    , mat_(other.mat_)
#endif
{
    other.nrows_ = other.ncols_ = other.local_rows_ = 0;
#ifdef HVDC_USE_PETSC
    other.mat_ = nullptr;
#endif
}

PETScMatrix& PETScMatrix::operator=(PETScMatrix&& other) noexcept {
    if (this != &other) {
        destroy();
        nrows_ = other.nrows_;
        ncols_ = other.ncols_;
        local_rows_ = other.local_rows_;
#ifdef HVDC_USE_PETSC
        mat_ = other.mat_;
        other.mat_ = nullptr;
#endif
        other.nrows_ = other.ncols_ = other.local_rows_ = 0;
    }
    return *this;
}

void PETScMatrix::create(Index nrows, Index ncols, Index nnz_per_row) {
    destroy();
    nrows_ = nrows;
    ncols_ = ncols;
    local_rows_ = PETSC_DECIDE;
#ifdef HVDC_USE_PETSC
    MatCreate(PETSC_COMM_WORLD, &mat_);
    MatSetSizes(mat_, PETSC_DECIDE, PETSC_DECIDE, nrows, ncols);
    MatSetFromOptions(mat_);
    MatSeqAIJSetPreallocation(mat_, nnz_per_row, nullptr);
    MatMPIAIJSetPreallocation(mat_, nnz_per_row, nullptr, nnz_per_row, nullptr);
#endif
}

void PETScMatrix::create(Index local_rows, Index global_rows, Index global_cols, Index nnz_per_row) {
    destroy();
    nrows_ = global_rows;
    ncols_ = global_cols;
    local_rows_ = local_rows;
#ifdef HVDC_USE_PETSC
    MatCreate(PETSC_COMM_WORLD, &mat_);
    MatSetSizes(mat_, local_rows, PETSC_DECIDE, global_rows, global_cols);
    MatSetFromOptions(mat_);
#endif
}

void PETScMatrix::destroy() {
#ifdef HVDC_USE_PETSC
    if (mat_) {
        MatDestroy(&mat_);
        mat_ = nullptr;
    }
#endif
    nrows_ = ncols_ = local_rows_ = 0;
}

void PETScMatrix::from_sparse_matrix(const SparseMatrix& m) {
    if (nrows_ != m.rows() || ncols_ != m.cols()) {
        create(m.rows(), m.cols(), 20);
    }
    zero();
    
    const auto& rp = m.row_ptr();
    const auto& ci = m.col_idx();
    const auto& v = m.values();
    
    for (Index r = 0; r < m.rows(); ++r) {
        for (Index i = rp[r]; i < rp[r + 1]; ++i) {
            add(r, ci[i], v[i]);
        }
    }
    assembly_begin();
    assembly_end();
}

void PETScMatrix::to_sparse_matrix(SparseMatrix& mat) const {
    mat.reset(nrows_, ncols_);
#ifdef HVDC_USE_PETSC
    const PetscInt* cols;
    const PetscScalar* vals;
    PetscInt ncols_r;
    
    for (PetscInt r = 0; r < nrows_; ++r) {
        MatGetRow(mat_, r, &ncols_r, &cols, &vals);
        for (PetscInt j = 0; j < ncols_r; ++j) {
            mat.add(r, cols[j], PetscRealPart(vals[j]));
        }
        MatRestoreRow(mat_, r, &ncols_r, &cols, &vals);
    }
#endif
}

void PETScMatrix::set(Index row, Index col, Real val) {
#ifdef HVDC_USE_PETSC
    MatSetValue(mat_, row, col, val, INSERT_VALUES);
#endif
}

void PETScMatrix::add(Index row, Index col, Real val) {
#ifdef HVDC_USE_PETSC
    MatSetValue(mat_, row, col, val, ADD_VALUES);
#endif
}

void PETScMatrix::set_row(Index row, const IndexVec& cols, const Vec& vals) {
#ifdef HVDC_USE_PETSC
    MatSetValues(mat_, 1, &row, static_cast<PetscInt>(cols.size()), 
                 cols.data(), vals.data(), INSERT_VALUES);
#endif
}

void PETScMatrix::add_row(Index row, const IndexVec& cols, const Vec& vals) {
#ifdef HVDC_USE_PETSC
    MatSetValues(mat_, 1, &row, static_cast<PetscInt>(cols.size()), 
                 cols.data(), vals.data(), ADD_VALUES);
#endif
}

void PETScMatrix::zero() {
#ifdef HVDC_USE_PETSC
    MatZeroEntries(mat_);
#endif
}

void PETScMatrix::zero_rows(const IndexVec& rows, Real diag_value) {
#ifdef HVDC_USE_PETSC
    IS is_rows;
    ISCreateGeneral(PETSC_COMM_WORLD, static_cast<PetscInt>(rows.size()), 
                     rows.data(), PETSC_COPY_VALUES, &is_rows);
    MatZeroRowsIS(mat_, is_rows, diag_value, nullptr, nullptr);
    ISDestroy(&is_rows);
#endif
}

void PETScMatrix::shift_diagonal(Real alpha) {
#ifdef HVDC_USE_PETSC
    MatShift(mat_, alpha);
#endif
}

void PETScMatrix::assembly_begin() {
#ifdef HVDC_USE_PETSC
    MatAssemblyBegin(mat_, MAT_FINAL_ASSEMBLY);
#endif
}

void PETScMatrix::assembly_end() {
#ifdef HVDC_USE_PETSC
    MatAssemblyEnd(mat_, MAT_FINAL_ASSEMBLY);
#endif
}

void PETScMatrix::multiply(const PETScVector& x, PETScVector& y) const {
#ifdef HVDC_USE_PETSC
    MatMult(mat_, x.petsc_vec(), y.petsc_vec());
#endif
}

Index PETScMatrix::nnz() const {
#ifdef HVDC_USE_PETSC
    PetscInt n;
    MatGetInfo(mat_, MAT_GLOBAL_SUM_NONZEROS, &n);
    return n;
#else
    return 0;
#endif
}

void PETScMatrix::set_preallocation_nnz(const IndexVec& d_nnz, const IndexVec& o_nnz) {
#ifdef HVDC_USE_PETSC
    MatMPIAIJSetPreallocation(mat_, 0, d_nnz.empty() ? nullptr : d_nnz.data(),
                              0, o_nnz.empty() ? nullptr : o_nnz.data());
#endif
}

// ===== PETScSolver =====

PETScSolver::~PETScSolver() { destroy(); }

void PETScSolver::create() {
    if (created_) return;
#ifdef HVDC_USE_PETSC
    KSPCreate(PETSC_COMM_WORLD, &ksp_);
    KSPGetPC(ksp_, &pc_);
#endif
    created_ = true;
}

void PETScSolver::destroy() {
#ifdef HVDC_USE_PETSC
    if (ksp_) {
        KSPDestroy(&ksp_);
        ksp_ = nullptr;
    }
#endif
    created_ = false;
}

void PETScSolver::set_operator(const PETScMatrix& A) {
    if (!created_) create();
#ifdef HVDC_USE_PETSC
    KSPSetOperators(ksp_, A.petsc_mat(), A.petsc_mat());
#endif
}

void PETScSolver::set_preconditioner(const std::string& type) {
    if (!created_) create();
#ifdef HVDC_USE_PETSC
    if (type == "jacobi")     PCSetType(pc_, PCJACOBI);
    else if (type == "bjacobi") PCSetType(pc_, PCBJACOBI);
    else if (type == "ilu")     PCSetType(pc_, PCILU);
    else if (type == "icc")     PCSetType(pc_, PCICC);
    else if (type == "gamg")    PCSetType(pc_, PCGAMG);
    else if (type == "hypre")   PCSetType(pc_, PCHYPRE);
    else if (type == "ml")      PCSetType(pc_, PCML);
    else                        PCSetType(pc_, type.c_str());
#endif
}

void PETScSolver::set_ksp_type(const std::string& type) {
    if (!created_) create();
#ifdef HVDC_USE_PETSC
    if (type == "cg")       KSPSetType(ksp_, KSPCG);
    else if (type == "gmres")  KSPSetType(ksp_, KSPGMRES);
    else if (type == "bicg")   KSPSetType(ksp_, KSPBICG);
    else if (type == "cgs")    KSPSetType(ksp_, KSPCGS);
    else if (type == "tfqmr")  KSPSetType(ksp_, KSPTFQMR);
    else if (type == "richardson") KSPSetType(ksp_, KSPRICHARDSON);
    else if (type == "preonly")    KSPSetType(ksp_, KSPPREONLY);
    else                           KSPSetType(ksp_, type.c_str());
#endif
}

void PETScSolver::set_tolerances(Real rtol, Real atol, Real dtol, Index max_iters) {
    if (!created_) create();
#ifdef HVDC_USE_PETSC
    KSPSetTolerances(ksp_, rtol, atol, dtol, max_iters);
#endif
}

bool PETScSolver::solve(const PETScVector& b, PETScVector& x) {
    if (!created_) create();
    converged_ = false;
    iterations_ = 0;
    residual_ = 0.0;
    reason_ = "not_run";
#ifdef HVDC_USE_PETSC
    KSPSolve(ksp_, b.petsc_vec(), x.petsc_vec());
    
    KSPConvergedReason r;
    KSPGetConvergedReason(ksp_, &r);
    KSPGetIterationNumber(ksp_, &iterations_);
    KSPGetResidualNorm(ksp_, &residual_);
    
    if (r > 0) {
        converged_ = true;
        switch (r) {
            case KSP_CONVERGED_RTOL: reason_ = "rtol"; break;
            case KSP_CONVERGED_ATOL: reason_ = "atol"; break;
            case KSP_CONVERGED_ITS:  reason_ = "its"; break;
            default: reason_ = "positive";
        }
    } else {
        switch (r) {
            case KSP_DIVERGED_NULL:       reason_ = "null"; break;
            case KSP_DIVERGED_ITS:        reason_ = "max_iters"; break;
            case KSP_DIVERGED_DTOL:       reason_ = "dtol"; break;
            case KSP_DIVERGED_BREAKDOWN:  reason_ = "breakdown"; break;
            case KSP_DIVERGED_BREAKDOWN_BICG: reason_ = "bicg_breakdown"; break;
            default: reason_ = "diverged";
        }
    }
#endif
    return converged_;
}

void PETScSolver::set_from_options() {
    if (!created_) create();
#ifdef HVDC_USE_PETSC
    KSPSetFromOptions(ksp_);
#endif
}

} // namespace hvdc
