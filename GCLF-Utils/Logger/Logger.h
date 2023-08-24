#pragma once
#include <cstdio>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
// spdlog
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/stopwatch.h"

namespace GCLF
{
namespace Utils 
{

/****************************/
/******** Logger ************/
/****************************/

class Logger
{
public:
  static std::shared_ptr<spdlog::logger> dev_logger;
  static std::shared_ptr<spdlog::logger> user_logger;
  static std::shared_ptr<spdlog::stopwatch> sw;

  static void InitLogger(
    spdlog::level::level_enum console_level = spdlog::level::trace,
    bool file_log = false,
    spdlog::level::level_enum file_level = spdlog::level::trace,
    const std::string& file_path = "logs/log.txt"
  );

  static void updateFileLog(bool file_log, spdlog::level::level_enum file_level, const std::string& file_path);

  static void StopProgram();
};

/****************************/
/******** Trace *************/
/****************************/

#if defined(TRACE_FUNC)
#define REP_FUNC {Logger::dev_logger->trace(__FUNCTION__);}
#else
#define REP_FUNC
#endif

/*************************/
/******** Time ***********/
/*************************/

#define SW_RESET Logger::sw->reset()
#define SW_COUNT Logger::sw->elapsed().count()

/****************************/
/******** Assert ************/
/****************************/

class AssertFailExcept
{
public:
};

#define GCLF_DEBUG
#if defined(GCLF_DEBUG)
#define ASSERT(condition, msg, ...) \
if(!(condition)) {\
  Logger::dev_logger->critical("in file {} line {}: " msg , __FILE__, __LINE__, ##__VA_ARGS__);\
  throw AssertFailExcept();}
#else
#define ASSERT(condition, msg, ...)
#endif

/*****************************/
/******** Breakout ***********/
/*****************************/

class BreakoutExcept
{
public:
};

#define BREAKOUT Logger::dev_logger->debug("break out."); throw BreakoutExcept()
#define UNUSED(var) (void)var


/*****************************/
/********  Memory  ***********/
/*****************************/

/** Copy from VolumeMesher::makePolyhedralMesh.cpp **/
#ifdef _MSC_VER
#include <windows.h>
#include <psapi.h>

// To ensure correct resolution of symbols, add Psapi.lib to TARGETLIBS
// and compile with -DPSAPI_VERSION=1

inline double getPeakMegabytesUsed()
{
  HANDLE hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, GetCurrentProcessId());
  if (NULL == hProcess) return 0;

  PROCESS_MEMORY_COUNTERS pmc;
  double mem = 0;
  if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)))
  {
    mem = pmc.PeakWorkingSetSize / 1048576.0;
  }

  CloseHandle(hProcess);
  return mem;
}

inline double getMegabytesUsed()
{
  HANDLE hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ, FALSE, GetCurrentProcessId());
  if (NULL == hProcess) return 0;

  PROCESS_MEMORY_COUNTERS pmc;
  double mem = 0;
  if (GetProcessMemoryInfo(hProcess, &pmc, sizeof(pmc)))
  {
    mem = pmc.WorkingSetSize / 1048576.0;
  }

  CloseHandle(hProcess);
  return mem;
}
#else
inline double getPeakMegabytesUsed() { return 0.0; }
inline double getMegabytesUsed() { return 0.0; }
#endif
}// namespace Utils
}// namespace GCLF