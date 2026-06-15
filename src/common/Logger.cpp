#include "common/Logger.hpp"
#include "common/MPIManager.hpp"
#include <chrono>
#include <iomanip>
#include <ctime>

namespace hvdc {

Logger& Logger::instance() {
    static Logger inst;
    return inst;
}

void Logger::set_log_file(const std::string& filename) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (log_file_.is_open()) log_file_.close();
    log_file_.open(filename, std::ios::out | std::ios::app);
}

void Logger::set_rank_filter(int rank) {
    rank_filter_ = rank;
}

bool Logger::should_log(LogLevel level) const {
    if (level < min_level_) return false;
    if (rank_filter_ >= 0) {
        auto& mpi = MPIManager::instance();
        if (mpi.rank() != rank_filter_) return false;
    }
    return true;
}

void Logger::log(LogLevel level, const std::string& message,
                 const char* file, int line, const char* func) {
    if (!should_log(level)) return;
    
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto now = std::chrono::system_clock::now();
    std::time_t now_t = std::chrono::system_clock::to_time_t(now);
    std::tm now_tm = *std::localtime(&now_t);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;
    
    auto& mpi = MPIManager::instance();
    
    std::ostringstream oss;
    oss << "[" << std::put_time(&now_tm, "%Y-%m-%d %H:%M:%S")
        << "." << std::setw(3) << std::setfill('0') << ms.count()
        << "]"
        << "[" << std::setw(5) << level_to_string(level) << "]"
        << "[R" << mpi.rank() << "]";
    
    if (file) {
        std::string f = file;
        size_t pos = f.find_last_of("/\\");
        if (pos != std::string::npos) f = f.substr(pos + 1);
        oss << "[" << f << ":" << line;
        if (func) oss << " " << func;
        oss << "]";
    }
    
    oss << " " << message << std::endl;
    
    if (console_output_) {
        if (level >= LogLevel::Error) {
            std::cerr << oss.str();
        } else {
            std::cout << oss.str();
        }
    }
    
    if (log_file_.is_open()) {
        log_file_ << oss.str();
        log_file_.flush();
    }
}

std::string level_to_string(LogLevel level) {
    switch (level) {
        case LogLevel::Trace:   return "TRACE";
        case LogLevel::Debug:   return "DEBUG";
        case LogLevel::Info:    return "INFO";
        case LogLevel::Warning: return "WARN";
        case LogLevel::Error:   return "ERROR";
        case LogLevel::Fatal:   return "FATAL";
        default:                return "UNKNOWN";
    }
}

LogLevel string_to_level(const std::string& str) {
    std::string s = str;
    for (auto& c : s) c = std::toupper(c);
    if (s == "TRACE")   return LogLevel::Trace;
    if (s == "DEBUG")   return LogLevel::Debug;
    if (s == "INFO")    return LogLevel::Info;
    if (s == "WARN" || s == "WARNING") return LogLevel::Warning;
    if (s == "ERROR")   return LogLevel::Error;
    if (s == "FATAL")   return LogLevel::Fatal;
    if (s == "SILENT" || s == "OFF") return LogLevel::Silent;
    return LogLevel::Info;
}

} // namespace hvdc
