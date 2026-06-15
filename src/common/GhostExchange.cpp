#include "common/GhostExchange.hpp"
#include <cassert>

namespace hvdc {

namespace detail {

template <typename T>
struct MPITypeMap {
    static constexpr bool is_specialized = false;
};

template <> struct MPITypeMap<Real> {
    static constexpr bool is_specialized = true;
#ifdef HVDC_USE_MPI
    static MPI_Datatype type() { return MPI_DOUBLE; }
#endif
};

template <> struct MPITypeMap<int> {
    static constexpr bool is_specialized = true;
#ifdef HVDC_USE_MPI
    static MPI_Datatype type() { return MPI_INT; }
#endif
};

template <> struct MPITypeMap<Index> {
    static constexpr bool is_specialized = true;
#ifdef HVDC_USE_MPI
    static MPI_Datatype type() { 
        return sizeof(Index) == sizeof(int) ? MPI_INT : MPI_LONG_LONG; 
    }
#endif
};

template <typename T, std::size_t N>
struct MPITypeMap<std::array<T, N>> {
    static constexpr bool is_specialized = true;
#ifdef HVDC_USE_MPI
    static MPI_Datatype type() {
        static MPI_Datatype arr_type = MPI_DATATYPE_NULL;
        static bool initialized = false;
        if (!initialized) {
            MPI_Type_contiguous(static_cast<int>(N), MPITypeMap<T>::type(), &arr_type);
            MPI_Type_commit(&arr_type);
            initialized = true;
        }
        return arr_type;
    }
#endif
};

} // namespace detail

template <typename T>
GhostExchange<T>::GhostExchange() = default;

template <typename T>
GhostExchange<T>::~GhostExchange() {
    clear();
}

template <typename T>
void GhostExchange<T>::clear() {
#ifdef HVDC_USE_MPI
    if (exchanging_) {
        for (auto& req : send_requests_) {
            if (req != MPI_REQUEST_NULL) {
                MPI_Status status;
                MPI_Wait(&req, &status);
            }
        }
        for (auto& req : recv_requests_) {
            if (req != MPI_REQUEST_NULL) {
                MPI_Status status;
                MPI_Wait(&req, &status);
            }
        }
    }
#endif
    neighbors_.clear();
    send_buffers_.clear();
    recv_buffers_.clear();
#ifdef HVDC_USE_MPI
    send_requests_.clear();
    recv_requests_.clear();
#endif
    active_ = false;
    exchanging_ = false;
}

template <typename T>
void GhostExchange<T>::initialize(const std::vector<NeighborCommInfo>& neighbors) {
    clear();
    neighbors_ = neighbors;
    
    std::size_t n = neighbors_.size();
    send_buffers_.resize(n);
    recv_buffers_.resize(n);
#ifdef HVDC_USE_MPI
    send_requests_.resize(n, MPI_REQUEST_NULL);
    recv_requests_.resize(n, MPI_REQUEST_NULL);
#endif
    
    for (std::size_t i = 0; i < n; ++i) {
        send_buffers_[i].resize(neighbors_[i].send_count);
        recv_buffers_[i].resize(neighbors_[i].recv_count);
    }
    
    active_ = !neighbors_.empty();
}

template <typename T>
void GhostExchange<T>::pack_send_data(const std::vector<T>& data) {
    for (std::size_t i = 0; i < neighbors_.size(); ++i) {
        const auto& info = neighbors_[i];
        auto& buf = send_buffers_[i];
        assert(info.send_count == static_cast<Index>(info.send_indices.size()));
        assert(buf.size() == static_cast<std::size_t>(info.send_count));
        
        for (Index j = 0; j < info.send_count; ++j) {
            Index idx = info.send_indices[j];
            assert(idx < static_cast<Index>(data.size()));
            buf[j] = data[idx];
        }
    }
}

template <typename T>
void GhostExchange<T>::unpack_recv_data(std::vector<T>& data) {
    for (std::size_t i = 0; i < neighbors_.size(); ++i) {
        const auto& info = neighbors_[i];
        const auto& buf = recv_buffers_[i];
        assert(info.recv_count == static_cast<Index>(info.recv_indices.size()));
        assert(buf.size() == static_cast<std::size_t>(info.recv_count));
        
        for (Index j = 0; j < info.recv_count; ++j) {
            Index idx = info.recv_indices[j];
            assert(idx < static_cast<Index>(data.size()));
            data[idx] = buf[j];
        }
    }
}

template <typename T>
void GhostExchange<T>::start_exchange(std::vector<T>& data) {
    if (!active_) return;
    
    MPIManager& mpi = MPIManager::instance();
    if (mpi.size() <= 1) return;
    
    pack_send_data(data);
    
#ifdef HVDC_USE_MPI
    MPI_Datatype dtype = detail::MPITypeMap<T>::type();
    
    for (std::size_t i = 0; i < neighbors_.size(); ++i) {
        const auto& info = neighbors_[i];
        
        if (info.recv_count > 0) {
            int count = static_cast<int>(recv_buffers_[i].size());
            MPI_Irecv(
                recv_buffers_[i].data(),
                count,
                dtype,
                info.rank,
                0,
                mpi.comm(),
                &recv_requests_[i]
            );
        } else {
            recv_requests_[i] = MPI_REQUEST_NULL;
        }
        
        if (info.send_count > 0) {
            int count = static_cast<int>(send_buffers_[i].size());
            MPI_Isend(
                send_buffers_[i].data(),
                count,
                dtype,
                info.rank,
                0,
                mpi.comm(),
                &send_requests_[i]
            );
        } else {
            send_requests_[i] = MPI_REQUEST_NULL;
        }
    }
#endif
    
    exchanging_ = true;
}

template <typename T>
void GhostExchange<T>::finish_exchange(std::vector<T>& data) {
    if (!active_ || !exchanging_) return;
    
    MPIManager& mpi = MPIManager::instance();
    if (mpi.size() <= 1) {
        exchanging_ = false;
        return;
    }
    
#ifdef HVDC_USE_MPI
    for (std::size_t i = 0; i < neighbors_.size(); ++i) {
        if (recv_requests_[i] != MPI_REQUEST_NULL) {
            MPI_Status status;
            MPI_Wait(&recv_requests_[i], &status);
            recv_requests_[i] = MPI_REQUEST_NULL;
        }
    }
    
    for (std::size_t i = 0; i < neighbors_.size(); ++i) {
        if (send_requests_[i] != MPI_REQUEST_NULL) {
            MPI_Status status;
            MPI_Wait(&send_requests_[i], &status);
            send_requests_[i] = MPI_REQUEST_NULL;
        }
    }
#endif
    
    unpack_recv_data(data);
    exchanging_ = false;
}

template <typename T>
void GhostExchange<T>::exchange(std::vector<T>& data) {
    start_exchange(data);
    finish_exchange(data);
}

template class GhostExchange<Real>;
template class GhostExchange<Index>;
template class GhostExchange<int>;
template class GhostExchange<std::array<Real, 3>>;
template class GhostExchange<std::array<Real, 6>>;

} // namespace hvdc
