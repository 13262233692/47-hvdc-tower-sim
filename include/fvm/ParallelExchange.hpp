#ifndef HVDC_FVM_PARALLEL_EXCHANGE_HPP
#define HVDC_FVM_PARALLEL_EXCHANGE_HPP

#include "common/Types.hpp"
#include "common/GhostExchange.hpp"
#include "fvm/Field.hpp"
#include <vector>
#include <memory>

namespace hvdc {
namespace fvm {

class Grid;
class BoundaryConditionManager;

struct ParallelPartitionInfo {
    Index local_num_cells = 0;
    Index local_num_faces = 0;
    Index global_num_cells = 0;
    Index global_num_faces = 0;
    
    Index local_cell_start = 0;
    Index local_face_start = 0;
    
    Index num_ghost_cells = 0;
    Index num_ghost_faces = 0;
    
    std::vector<int> neighbor_ranks;
    std::vector<IndexVec> send_cell_indices;
    std::vector<IndexVec> recv_cell_indices;
    std::vector<IndexVec> send_face_indices;
    std::vector<IndexVec> recv_face_indices;
};

class FVMParallelExchange {
public:
    FVMParallelExchange();
    ~FVMParallelExchange();
    
    void initialize(const Grid* grid, 
                    const ParallelPartitionInfo& partition);
    
    void exchange_scalar_field(ScalarField& field);
    void exchange_vector_field(VectorField& field);
    void exchange_velocity_field(VelocityField& field);
    void exchange_pressure_field(PressureField& field);
    
    void start_scalar_exchange(ScalarField& field);
    void finish_scalar_exchange(ScalarField& field);
    
    void start_vector_exchange(VectorField& field);
    void finish_vector_exchange(VectorField& field);
    
    void exchange_face_data(std::vector<Real>& face_data);
    void exchange_face_velocity(std::vector<Vec3>& face_vel);
    
    bool is_parallel() const { return is_parallel_; }
    Index num_ghost_cells() const { return partition_.num_ghost_cells; }
    Index num_ghost_faces() const { return partition_.num_ghost_faces; }
    
    const ParallelPartitionInfo& partition_info() const { return partition_; }
    
    void clear();

private:
    bool is_parallel_ = false;
    const Grid* grid_ = nullptr;
    ParallelPartitionInfo partition_;
    
    std::unique_ptr<GhostExchange<Real>> scalar_exchange_;
    std::unique_ptr<GhostExchange<std::array<Real, 3>>> vector_exchange_;
    std::unique_ptr<GhostExchange<Real>> face_scalar_exchange_;
    std::unique_ptr<GhostExchange<std::array<Real, 3>>> face_vector_exchange_;
    
    mutable std::vector<std::array<Real, 3>> vec3_scratch_;
    mutable std::vector<std::array<Real, 3>> face_vec3_scratch_;
    
    void build_neighbor_info(
        const std::vector<int>& ranks,
        const std::vector<IndexVec>& send_idx,
        const std::vector<IndexVec>& recv_idx,
        std::vector<NeighborCommInfo>& neighbors) const;
};

} // namespace fvm
} // namespace hvdc

#endif // HVDC_FVM_PARALLEL_EXCHANGE_HPP
