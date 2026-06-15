#include "common/Timer.hpp"
#include "common/MPIManager.hpp"
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace hvdc {

Timer::Timer() {
    tic_time_ = Clock::now();
}

void Timer::start(const std::string& name) {
    auto& e = entries_[name];
    e.start = Clock::now();
    e.running = true;
}

void Timer::stop(const std::string& name) {
    auto it = entries_.find(name);
    if (it == entries_.end()) return;
    auto& e = it->second;
    if (!e.running) return;
    
    auto now = Clock::now();
    std::chrono::duration<double> diff = now - e.start;
    e.total += diff.count();
    e.calls++;
    e.running = false;
}

void Timer::reset(const std::string& name) {
    auto it = entries_.find(name);
    if (it != entries_.end()) {
        it->second.total = 0.0;
        it->second.calls = 0;
        it->second.running = false;
    }
}

void Timer::reset_all() {
    entries_.clear();
}

Real Timer::elapsed(const std::string& name) const {
    auto it = entries_.find(name);
    if (it == entries_.end()) return 0.0;
    return it->second.total;
}

Index Timer::count(const std::string& name) const {
    auto it = entries_.find(name);
    if (it == entries_.end()) return 0;
    return it->second.calls;
}

Real Timer::average(const std::string& name) const {
    auto it = entries_.find(name);
    if (it == entries_.end() || it->second.calls == 0) return 0.0;
    return it->second.total / static_cast<Real>(it->second.calls);
}

std::vector<std::pair<std::string, Real>> Timer::sorted_times() const {
    std::vector<std::pair<std::string, Real>> result;
    result.reserve(entries_.size());
    for (const auto& [name, e] : entries_) {
        result.emplace_back(name, e.total);
    }
    std::sort(result.begin(), result.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    return result;
}

void Timer::print_summary() const {
    auto& mpi = MPIManager::instance();
    if (!mpi.is_root()) return;
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Timing Summary" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << std::left << std::setw(30) << "Region" 
              << std::right << std::setw(15) << "Total (s)"
              << std::setw(12) << "Calls"
              << std::setw(15) << "Avg (ms)" << std::endl;
    std::cout << std::string(72, '-') << std::endl;
    
    Real total_time = 0.0;
    for (const auto& [name, e] : entries_) total_time += e.total;
    
    auto sorted = sorted_times();
    for (const auto& [name, time] : sorted) {
        Index calls = count(name);
        Real avg = (calls > 0) ? (time / calls) * 1000.0 : 0.0;
        std::cout << std::left << std::setw(30) << name.substr(0, 28)
                  << std::right << std::setw(15) << std::fixed << std::setprecision(4) << time
                  << std::setw(12) << calls
                  << std::setw(15) << std::fixed << std::setprecision(3) << avg << std::endl;
    }
    
    std::cout << std::string(72, '-') << std::endl;
    std::cout << std::left << std::setw(30) << "Total Accounted"
              << std::right << std::setw(15) << std::fixed << std::setprecision(4) << total_time << std::endl;
    std::cout << "========================================\n" << std::endl;
}

void Timer::tic() {
    tic_time_ = Clock::now();
}

Real Timer::toc() const {
    auto now = Clock::now();
    std::chrono::duration<double> diff = now - tic_time_;
    return diff.count();
}

} // namespace hvdc
