#ifndef HVDC_FEM_PARALLEL_ASSEMBLER_HPP
#define HVDC_FEM_PARALLEL_ASSEMBLER_HPP

#include "common/Types.hpp"
#include "common/GhostExchange.hpp"
#include "fem/FEModel.hpp"
#include "fem/Assembler.hpp"
#include <vector>
#include <memory>

namespace hvdc {
namespace fem {

struct FEMPartitionInfo {
    Index local_num_nodes = 0;
    Index local_num_elements = 0;
    Index global_num_nodes = 0;
    Index global_num_elements = 0;
    
    Index local_dof_start = 0;
    Index local_num_dofs = 0;
    Index global_num_dofs = 0;
    
    Index num_ghost_nodes = 0;
    Index num_ghost_dofs = 0;
    
    std::vector<int> neighbor_ranks;
    std::vector<IndexVec> send_node_indices;
    std::vector<IndexVec> recv_node_indices;
    std::vector<IndexVec> send_dof_indices;
    std::vector<IndexVec> recv_dof_indices;
    
    IndexVec local_element_ids;
    IndexVec local_node_ids;
};

class ParallelAssembler {
public:
    ParallelAssembler();
    ~ParallelAssembler();
    
    void initialize(FEModel* model, 
                    Assembler* assembler,
                    const FEMPartitionInfo& partition);
    
    void assemble_tangent_stiffness_parallel(
        const Vector& displacement,
        AssemblyResult& result,
        bool include_geometric = true,
        bool include_mass = false);
    
    void assemble_internal_forces_parallel(
        const Vector& displacement,
        Vector& F_int);
    
    void exchange_displacement(Vector& displacement);
    void exchange_force(Vector& force);
    void exchange_solution(Vector& solution);
    
    void start_displacement_exchange(Vector& displacement);
    void finish_displacement_exchange(Vector& displacement);
    
    void sum_assembled_matrix(SparseMatrix& K);
    void sum_assembled_vector(Vector& F);
    
    bool is_parallel() const { return is_parallel_; }
    const FEMPartitionInfo& partition_info() const { return partition_; }
    
    Index local_dof_start() const { return partition_.local_dof_start; }
    Index local_num_dofs() const { return partition_.local_num_dofs; }
    Index global_num_dofs() const { return partition_.global_num_dofs; }
    
    void clear();

private:
    bool is_parallel_ = false;
    FEModel* model_ = nullptr;
    Assembler* assembler_ = nullptr;
    FEMPartitionInfo partition_;
    
    std::unique_ptr<GhostExchange<Real>> dof_exchange_;
    std::unique_ptr<GhostExchange<std::array<Real, 3>>> node_exchange_;
    
    mutable std::vector<Real> dof_scratch_;
    mutable std::vector<std::array<Real, 3>> node_scratch_;
    
    void build_dof_neighbors(std::vector<NeighborCommInfo>& neighbors) const;
    void build_node_neighbors(std::vector<NeighborCommInfo>& neighbors) const;
};

} // namespace fem
} // namespace hvdc

#endif // HVDC_FEM_PARALLEL_ASSEMBLER_HPP
