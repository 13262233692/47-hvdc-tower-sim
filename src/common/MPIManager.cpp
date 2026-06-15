#include "common/MPIManager.hpp"
#include <iostream>
#include <cstring>

namespace hvdc {

MPIManager& MPIManager::instance() {
    static MPIManager inst;
    return inst;
}

MPIManager::~MPIManager() {
    if (initialized_ && !finalize_called_) {
        finalize();
    }
}

void MPIManager::initialize(int* argc, char*** argv) {
    if (initialized_) return;
    
#ifdef HVDC_USE_MPI
    int provided;
    MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &provided);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &size_);
    
    char name[MPI_MAX_PROCESSOR_NAME];
    int len;
    MPI_Get_processor_name(name, &len);
    proc_name_ = std::string(name, len);
    
    initialized_ = true;
#else
    rank_ = 0;
    size_ = 1;
    proc_name_ = "serial";
    initialized_ = true;
#endif
}

void MPIManager::finalize() {
    if (!initialized_ || finalize_called_) return;
    
#ifdef HVDC_USE_MPI
    MPI_Finalize();
#endif
    finalize_called_ = true;
}

void MPIManager::barrier() const {
#ifdef HVDC_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

double MPIManager::wtime() const {
#ifdef HVDC_USE_MPI
    return MPI_Wtime();
#else
    return static_cast<double>(clock()) / CLOCKS_PER_SEC;
#endif
}

void MPIManager::partition_range(Index total, Index& local_start, Index& local_end) const {
    Index per_proc = total / size_;
    Index remainder = total % size_;
    
    if (rank_ < remainder) {
        local_start = rank_ * (per_proc + 1);
        local_end = local_start + per_proc + 1;
    } else {
        local_start = remainder * (per_proc + 1) + (rank_ - remainder) * per_proc;
        local_end = local_start + per_proc;
    }
}

void MPIManager::scatter_indices(const IndexVec& global_indices, IndexVec& local_indices,
                                  std::vector<IndexVec>& send_counts) const {
    Index total = static_cast<Index>(global_indices.size());
    send_counts.resize(size_);
    
    for (int r = 0; r < size_; ++r) {
        Index s, e;
        Index rt = r;
        Index per = total / size_;
        Index rem = total % size_;
        if (rt < rem) {
            s = rt * (per + 1);
            e = s + per + 1;
        } else {
            s = rem * (per + 1) + (rt - rem) * per;
            e = s + per;
        }
        send_counts[r].assign(global_indices.begin() + s, global_indices.begin() + e);
    }
    
    local_indices = send_counts[rank_];
}

} // namespace hvdc
