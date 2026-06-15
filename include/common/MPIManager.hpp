#ifndef HVDC_COMMON_MPI_MANAGER_HPP
#define HVDC_COMMON_MPI_MANAGER_HPP

#include "common/Types.hpp"
#include <vector>
#include <string>
#include <stdexcept>

#ifdef HVDC_USE_MPI
#include <mpi.h>
#endif

namespace hvdc {

class MPIManager {
public:
    static MPIManager& instance();
    
    void initialize(int* argc, char*** argv);
    void finalize();
    
    bool is_initialized() const { return initialized_; }
    bool has_mpi() const {
#ifdef HVDC_USE_MPI
        return true;
#else
        return false;
#endif
    }
    
    int rank() const { return rank_; }
    int size() const { return size_; }
    bool is_root() const { return rank_ == 0; }
    
#ifdef HVDC_USE_MPI
    MPI_Comm comm() const { return comm_; }
#endif
    
    void barrier() const;
    
    template <typename T>
    T broadcast(T value, int root = 0) const;
    
    template <typename T>
    void broadcast(std::vector<T>& data, int root = 0) const;
    
    template <typename T>
    void broadcast(T* data, std::size_t count, int root = 0) const;
    
    template <typename T>
    T reduce_sum(T value, int root = 0) const;
    
    template <typename T>
    void reduce_sum(const std::vector<T>& send, std::vector<T>& recv, int root = 0) const;
    
    template <typename T>
    T allreduce_sum(T value) const;
    
    template <typename T>
    void allreduce_sum(const std::vector<T>& send, std::vector<T>& recv) const;
    
    template <typename T>
    T allreduce_max(T value) const;
    
    template <typename T>
    T allreduce_min(T value) const;
    
    void partition_range(Index total, Index& local_start, Index& local_end) const;
    
    void scatter_indices(const IndexVec& global_indices, IndexVec& local_indices,
                         std::vector<IndexVec>& send_counts) const;
    
    std::string processor_name() const { return proc_name_; }
    
    double wtime() const;

private:
    MPIManager() = default;
    ~MPIManager();
    
    MPIManager(const MPIManager&) = delete;
    MPIManager& operator=(const MPIManager&) = delete;
    
    bool initialized_ = false;
    bool finalize_called_ = false;
    int rank_ = 0;
    int size_ = 1;
    std::string proc_name_ = "main";
    
#ifdef HVDC_USE_MPI
    MPI_Comm comm_ = MPI_COMM_WORLD;
#endif
};

} // namespace hvdc

#endif // HVDC_COMMON_MPI_MANAGER_HPP
