#ifndef HVDC_COMMON_GHOST_EXCHANGE_HPP
#define HVDC_COMMON_GHOST_EXCHANGE_HPP

#include "common/Types.hpp"
#include "common/MPIManager.hpp"
#include <vector>
#include <array>
#include <cstring>
#include <stdexcept>

#ifdef HVDC_USE_MPI
#include <mpi.h>
#endif

namespace hvdc {

enum class ExchangePattern : UInt8 {
    NearestNeighbor = 0,
    AllToAll = 1,
    Ring = 2
};

struct NeighborCommInfo {
    int rank = -1;
    Index send_count = 0;
    Index recv_count = 0;
    IndexVec send_indices;
    IndexVec recv_indices;
};

template <typename T>
class GhostExchange {
public:
    GhostExchange();
    ~GhostExchange();
    
    void initialize(const std::vector<NeighborCommInfo>& neighbors);
    
    void start_exchange(std::vector<T>& data);
    void finish_exchange(std::vector<T>& data);
    
    void exchange(std::vector<T>& data);
    
    Index num_neighbors() const { return neighbors_.size(); }
    const std::vector<NeighborCommInfo>& neighbors() const { return neighbors_; }
    
    bool is_active() const { return active_; }
    
    void clear();

private:
    bool active_ = false;
    std::vector<NeighborCommInfo> neighbors_;
    
    std::vector<std::vector<T>> send_buffers_;
    std::vector<std::vector<T>> recv_buffers_;
    
#ifdef HVDC_USE_MPI
    std::vector<MPI_Request> send_requests_;
    std::vector<MPI_Request> recv_requests_;
#endif
    
    bool exchanging_ = false;
    
    void pack_send_data(const std::vector<T>& data);
    void unpack_recv_data(std::vector<T>& data);
    
#ifdef HVDC_USE_MPI
    static MPI_Datatype mpi_type();
#endif
};

} // namespace hvdc

#endif // HVDC_COMMON_GHOST_EXCHANGE_HPP
