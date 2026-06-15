#include "fem/ParallelAssembler.hpp"
#include "common/MPIManager.hpp"
#include <cassert>

namespace hvdc {
namespace fem {

ParallelAssembler::ParallelAssembler() = default;
ParallelAssembler::~ParallelAssembler() { clear(); }

void ParallelAssembler::clear() {
    dof_exchange_.reset();
    node_exchange_.reset();
    is_parallel_ = false;
    model_ = nullptr;
    assembler_ = nullptr;
}

void ParallelAssembler::build_dof_neighbors(
    std::vector<NeighborCommInfo>& neighbors) const
{
    neighbors.clear();
    neighbors.reserve(partition_.neighbor_ranks.size());
    
    for (std::size_t i = 0; i < partition_.neighbor_ranks.size(); ++i) {
        NeighborCommInfo info;
        info.rank = partition_.neighbor_ranks[i];
        info.send_indices = partition_.send_dof_indices[i];
        info.recv_indices = partition_.recv_dof_indices[i];
        info.send_count = static_cast<Index>(info.send_indices.size());
        info.recv_count = static_cast<Index>(info.recv_indices.size());
        neighbors.push_back(info);
    }
}

void ParallelAssembler::build_node_neighbors(
    std::vector<NeighborCommInfo>& neighbors) const
{
    neighbors.clear();
    neighbors.reserve(partition_.neighbor_ranks.size());
    
    for (std::size_t i = 0; i < partition_.neighbor_ranks.size(); ++i) {
        NeighborCommInfo info;
        info.rank = partition_.neighbor_ranks[i];
        info.send_indices = partition_.send_node_indices[i];
        info.recv_indices = partition_.recv_node_indices[i];
        info.send_count = static_cast<Index>(info.send_indices.size());
        info.recv_count = static_cast<Index>(info.recv_indices.size());
        neighbors.push_back(info);
    }
}

void ParallelAssembler::initialize(
    FEModel* model,
    Assembler* assembler,
    const FEMPartitionInfo& partition)
{
    clear();
    model_ = model;
    assembler_ = assembler;
    partition_ = partition;
    
    MPIManager& mpi = MPIManager::instance();
    is_parallel_ = mpi.has_mpi() && mpi.size() > 1 && !partition.neighbor_ranks.empty();
    
    if (!is_parallel_) return;
    
    std::vector<NeighborCommInfo> dof_neighbors;
    build_dof_neighbors(dof_neighbors);
    
    std::vector<NeighborCommInfo> node_neighbors;
    build_node_neighbors(node_neighbors);
    
    dof_exchange_ = std::make_unique<GhostExchange<Real>>();
    dof_exchange_->initialize(dof_neighbors);
    
    node_exchange_ = std::make_unique<GhostExchange<std::array<Real, 3>>>();
    node_exchange_->initialize(node_neighbors);
}

void ParallelAssembler::exchange_displacement(Vector& displacement) {
    if (!is_parallel_) return;
    std::vector<Real> data(displacement.size());
    for (Index i = 0; i < displacement.size(); ++i) {
        data[i] = displacement[i];
    }
    dof_exchange_->exchange(data);
    for (Index i = 0; i < displacement.size(); ++i) {
        displacement[i] = data[i];
    }
}

void ParallelAssembler::exchange_force(Vector& force) {
    exchange_displacement(force);
}

void ParallelAssembler::exchange_solution(Vector& solution) {
    exchange_displacement(solution);
}

void ParallelAssembler::start_displacement_exchange(Vector& displacement) {
    if (!is_parallel_) return;
    
    dof_scratch_.resize(displacement.size());
    for (Index i = 0; i < displacement.size(); ++i) {
        dof_scratch_[i] = displacement[i];
    }
    dof_exchange_->start_exchange(dof_scratch_);
}

void ParallelAssembler::finish_displacement_exchange(Vector& displacement) {
    if (!is_parallel_) return;
    
    dof_exchange_->finish_exchange(dof_scratch_);
    for (Index i = 0; i < displacement.size(); ++i) {
        displacement[i] = dof_scratch_[i];
    }
}

void ParallelAssembler::sum_assembled_matrix(SparseMatrix& K) {
    if (!is_parallel_) return;
    
    MPIManager& mpi = MPIManager::instance();
#ifdef HVDC_USE_MPI
    IndexVec local_rows, local_cols;
    Vec local_vals;
    
    K.to_triplets(local_rows, local_cols, local_vals);
    
    int local_nnz = static_cast<int>(local_vals.size());
    std::vector<int> all_nnz(mpi.size());
    MPI_Allgather(&local_nnz, 1, MPI_INT, all_nnz.data(), 1, MPI_INT, mpi.comm());
    
    std::vector<int> displs(mpi.size(), 0);
    int total_nnz = 0;
    for (int r = 0; r < mpi.size(); ++r) {
        displs[r] = total_nnz;
        total_nnz += all_nnz[r];
    }
    
    std::vector<Real> all_vals(total_nnz);
    std::vector<int> all_rows(total_nnz);
    std::vector<int> all_cols(total_nnz);
    
    std::vector<int> local_rows_int(local_rows.begin(), local_rows.end());
    std::vector<int> local_cols_int(local_cols.begin(), local_cols.end());
    
    MPI_Allgatherv(local_vals.data(), local_nnz, MPI_DOUBLE,
                   all_vals.data(), all_nnz.data(), displs.data(), MPI_DOUBLE,
                   mpi.comm());
    MPI_Allgatherv(local_rows_int.data(), local_nnz, MPI_INT,
                   all_rows.data(), all_nnz.data(), displs.data(), MPI_INT,
                   mpi.comm());
    MPI_Allgatherv(local_cols_int.data(), local_nnz, MPI_INT,
                   all_cols.data(), all_nnz.data(), displs.data(), MPI_INT,
                   mpi.comm());
    
    K.reset(K.rows(), K.cols());
    for (int i = 0; i < total_nnz; ++i) {
        K.add(static_cast<Index>(all_rows[i]), 
              static_cast<Index>(all_cols[i]), 
              all_vals[i]);
    }
    K.finalize();
#endif
}

void ParallelAssembler::sum_assembled_vector(Vector& F) {
    if (!is_parallel_) return;
    
    MPIManager& mpi = MPIManager::instance();
    std::vector<Real> local(F.size());
    for (Index i = 0; i < F.size(); ++i) {
        local[i] = F[i];
    }
    
    std::vector<Real> global(F.size(), 0.0);
    mpi.allreduce_sum(local, global);
    
    for (Index i = 0; i < F.size(); ++i) {
        F[i] = global[i];
    }
}

void ParallelAssembler::assemble_tangent_stiffness_parallel(
    const Vector& displacement,
    AssemblyResult& result,
    bool include_geometric,
    bool include_mass)
{
    if (!is_parallel_ || !assembler_) {
        if (assembler_) {
            assembler_->assemble_tangent_stiffness(
                displacement, result, include_geometric, include_mass, false);
        }
        return;
    }
    
    assembler_->assemble_tangent_stiffness(
        displacement, result, include_geometric, include_mass, false);
    
    sum_assembled_vector(result.F_int);
    sum_assembled_vector(result.F_ext);
    sum_assembled_vector(result.Residual);
}

void ParallelAssembler::assemble_internal_forces_parallel(
    const Vector& displacement,
    Vector& F_int)
{
    if (!is_parallel_ || !assembler_) {
        if (assembler_) {
            assembler_->assemble_internal_forces(displacement, F_int, false);
        }
        return;
    }
    
    assembler_->assemble_internal_forces(displacement, F_int, false);
    sum_assembled_vector(F_int);
}

} // namespace fem
} // namespace hvdc
