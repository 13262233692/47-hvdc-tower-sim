#ifndef HVDC_COMMON_TIMER_HPP
#define HVDC_COMMON_TIMER_HPP

#include "common/Types.hpp"
#include <string>
#include <chrono>
#include <unordered_map>
#include <vector>

namespace hvdc {

class Timer {
public:
    Timer();
    
    void register_section(const std::string& name) { (void)name; }
    void start(const std::string& name);
    void stop(const std::string& name);
    void reset(const std::string& name);
    void reset_all();
    
    Real elapsed(const std::string& name) const;
    Index count(const std::string& name) const;
    Real average(const std::string& name) const;
    
    std::vector<std::pair<std::string, Real>> sorted_times() const;
    
    void print_summary() const;
    
    void tic();
    Real toc() const;

private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = Clock::time_point;
    
    struct TimeEntry {
        TimePoint start;
        Real total = 0.0;
        Index calls = 0;
        bool running = false;
    };
    
    mutable std::unordered_map<std::string, TimeEntry> entries_;
    TimePoint tic_time_;
};

} // namespace hvdc

#endif // HVDC_COMMON_TIMER_HPP
