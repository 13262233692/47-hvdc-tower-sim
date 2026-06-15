#include "fvm/ParallelExchange.hpp"
#include "fvm/Grid.hpp"
#include <cassert>

namespace hvdc {
namespace fvm {

FVMParallelExchange::FVMParallelExchange() = default;
FVMParallelExchange::~FVMParallelExchange() { clear(); }

void FVMParallelExchange::clear() {
    scalar_exchange_.reset();
    vector_exchange_.reset();
    face_scalar_exchange_.reset();
    face_vector_exchange_.reset();
    is_parallel_ = false;
    grid_ = nullptr;
}

void FVMParallelExchange::build_neighbor_info(
    const std::vector<int>& ranks,
    const std::vector<IndexVec>& send_idx,
    const std::vector<IndexVec>& recv_idx,
    std::vector<NeighborCommInfo>& neighbors) const
{
    neighbors.clear();
    neighbors.reserve(ranks.size());
    
    for (std::size_t i = 0; i < ranks.size(); ++i) {
        NeighborCommInfo info;
        info.rank = ranks[i];
        info.send_indices = send_idx[i];
        info.recv_indices = recv_idx[i];
        info.send_count = static_cast<Index>(send_idx[i].size());
        info.recv_count = static_cast<Index>(recv_idx[i].size());
        neighbors.push_back(info);
    }
}

void FVMParallelExchange::initialize(
    const Grid* grid,
    const ParallelPartitionInfo& partition)
{
    clear();
    grid_ = grid;
    partition_ = partition;
    
    MPIManager& mpi = MPIManager::instance();
    is_parallel_ = mpi.has_mpi() && mpi.size() > 1 && !partition.neighbor_ranks.empty();
    
    if (!is_parallel_) return;
    
    std::vector<NeighborCommInfo> cell_neighbors;
    build_neighbor_info(
        partition.neighbor_ranks,
        partition.send_cell_indices,
        partition.recv_cell_indices,
        cell_neighbors
    );
    
    std::vector<NeighborCommInfo> face_neighbors;
    build_neighbor_info(
        partition.neighbor_ranks,
        partition.send_face_indices,
        partition.recv_face_indices,
        face_neighbors
    );
    
    scalar_exchange_ = std::make_unique<GhostExchange<Real>>();
    scalar_exchange_->initialize(cell_neighbors);
    
    vector_exchange_ = std::make_unique<GhostExchange<std::array<Real, 3>>>();
    vector_exchange_->initialize(cell_neighbors);
    
    face_scalar_exchange_ = std::make_unique<GhostExchange<Real>>();
    face_scalar_exchange_->initialize(face_neighbors);
    
    face_vector_exchange_ = std::make_unique<GhostExchange<std::array<Real, 3>>>();
    face_vector_exchange_->initialize(face_neighbors);
}

void FVMParallelExchange::exchange_scalar_field(ScalarField& field) {
    if (!is_parallel_) return;
    scalar_exchange_->exchange(field.raw());
}

void FVMParallelExchange::exchange_vector_field(VectorField& field) {
    if (!is_parallel_) return;
    
    std::vector<std::array<Real, 3>> vec_data(field.size());
    for (Index i = 0; i < field.size(); ++i) {
        const Vec3& v = field[i];
        vec_data[i][0] = v[0];
        vec_data[i][1] = v[1];
        vec_data[i][2] = v[2];
    }
    
    vector_exchange_->exchange(vec_data);
    
    for (Index i = 0; i < field.size(); ++i) {
        Vec3& v = field[i];
        v[0] = vec_data[i][0];
        v[1] = vec_data[i][1];
        v[2] = vec_data[i][2];
    }
}

void FVMParallelExchange::exchange_velocity_field(VelocityField& field) {
    exchange_vector_field(static_cast<VectorField&>(field));
}

void FVMParallelExchange::exchange_pressure_field(PressureField& field) {
    exchange_scalar_field(static_cast<ScalarField&>(field));
}

void FVMParallelExchange::start_scalar_exchange(ScalarField& field) {
    if (!is_parallel_) return;
    scalar_exchange_->start_exchange(field.raw());
}

void FVMParallelExchange::finish_scalar_exchange(ScalarField& field) {
    if (!is_parallel_) return;
    scalar_exchange_->finish_exchange(field.raw());
}

void FVMParallelExchange::start_vector_exchange(VectorField& field) {
    if (!is_parallel_) return;
    
    vec3_scratch_.resize(field.size());
    for (Index i = 0; i < field.size(); ++i) {
        const Vec3& v = field[i];
        vec3_scratch_[i][0] = v[0];
        vec3_scratch_[i][1] = v[1];
        vec3_scratch_[i][2] = v[2];
    }
    
    vector_exchange_->start_exchange(vec3_scratch_);
}

void FVMParallelExchange::finish_vector_exchange(VectorField& field) {
    if (!is_parallel_) return;
    
    vector_exchange_->finish_exchange(vec3_scratch_);
    
    for (Index i = 0; i < field.size(); ++i) {
        Vec3& v = field[i];
        v[0] = vec3_scratch_[i][0];
        v[1] = vec3_scratch_[i][1];
        v[2] = vec3_scratch_[i][2];
    }
}

void FVMParallelExchange::exchange_face_data(std::vector<Real>& face_data) {
    if (!is_parallel_) return;
    face_scalar_exchange_->exchange(face_data);
}

void FVMParallelExchange::exchange_face_velocity(std::vector<Vec3>& face_vel) {
    if (!is_parallel_) return;
    
    std::vector<std::array<Real, 3>> vec_data(face_vel.size());
    for (std::size_t i = 0; i < face_vel.size(); ++i) {
        vec_data[i][0] = face_vel[i][0];
        vec_data[i][1] = face_vel[i][1];
        vec_data[i][2] = face_vel[i][2];
    }
    
    face_vector_exchange_->exchange(vec_data);
    
    for (std::size_t i = 0; i < face_vel.size(); ++i) {
        face_vel[i][0] = vec_data[i][0];
        face_vel[i][1] = vec_data[i][1];
        face_vel[i][2] = vec_data[i][2];
    }
}

} // namespace fvm
} // namespace hvdc
