#ifndef HVDC_COMMON_LOGGER_HPP
#define HVDC_COMMON_LOGGER_HPP

#include "common/Types.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <mutex>
#include <memory>

namespace hvdc {

enum class LogLevel : UInt8 {
    Trace = 0,
    Debug = 1,
    Info = 2,
    Warning = 3,
    Error = 4,
    Fatal = 5,
    Silent = 6
};

class Logger {
public:
    static Logger& instance();
    
    void set_level(LogLevel level) { min_level_ = level; }
    LogLevel level() const { return min_level_; }
    
    void set_log_file(const std::string& filename);
    void set_console_output(bool enable) { console_output_ = enable; }
    void set_rank_filter(int rank);
    
    void log(LogLevel level, const std::string& message,
             const char* file = nullptr, int line = 0, const char* func = nullptr);
    
    bool should_log(LogLevel level) const;
    
    LogLevel min_level_ = LogLevel::Info;
    bool console_output_ = true;
    std::ofstream log_file_;
    int rank_filter_ = -1;
    mutable std::mutex mutex_;
};

std::string level_to_string(LogLevel level);
LogLevel string_to_level(const std::string& str);

#define HVDC_LOG_TRACE(msg) \
    do { \
        std::ostringstream _oss_; _oss_ << msg; \
        hvdc::Logger::instance().log(hvdc::LogLevel::Trace, _oss_.str(), __FILE__, __LINE__, __func__); \
    } while(0)

#define HVDC_LOG_DEBUG(msg) \
    do { \
        std::ostringstream _oss_; _oss_ << msg; \
        hvdc::Logger::instance().log(hvdc::LogLevel::Debug, _oss_.str(), __FILE__, __LINE__, __func__); \
    } while(0)

#define HVDC_LOG_INFO(msg) \
    do { \
        std::ostringstream _oss_; _oss_ << msg; \
        hvdc::Logger::instance().log(hvdc::LogLevel::Info, _oss_.str(), __FILE__, __LINE__, __func__); \
    } while(0)

#define HVDC_LOG_WARNING(msg) \
    do { \
        std::ostringstream _oss_; _oss_ << msg; \
        hvdc::Logger::instance().log(hvdc::LogLevel::Warning, _oss_.str(), __FILE__, __LINE__, __func__); \
    } while(0)

#define HVDC_LOG_ERROR(msg) \
    do { \
        std::ostringstream _oss_; _oss_ << msg; \
        hvdc::Logger::instance().log(hvdc::LogLevel::Error, _oss_.str(), __FILE__, __LINE__, __func__); \
    } while(0)

#define HVDC_LOG_FATAL(msg) \
    do { \
        std::ostringstream _oss_; _oss_ << msg; \
        hvdc::Logger::instance().log(hvdc::LogLevel::Fatal, _oss_.str(), __FILE__, __LINE__, __func__); \
        std::abort(); \
    } while(0)

#define HVDC_CHECK(cond, msg) \
    do { \
        if (!(cond)) { \
            HVDC_LOG_FATAL("Check failed: " #cond " - " << msg); \
        } \
    } while(0)

} // namespace hvdc

#endif // HVDC_COMMON_LOGGER_HPP
