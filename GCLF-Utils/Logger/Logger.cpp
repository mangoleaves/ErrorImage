#include "Logger.h"

namespace GCLF
{
namespace Utils
{
bool logger_initialized = false;
std::shared_ptr<spdlog::logger> Logger::dev_logger = nullptr;
std::shared_ptr<spdlog::logger> Logger::user_logger = nullptr;
std::shared_ptr<spdlog::stopwatch> Logger::sw = nullptr;

void Logger::InitLogger(
  spdlog::level::level_enum console_level,
  bool file_log,
  spdlog::level::level_enum file_level,
  const std::string& file_path)
{
  if (logger_initialized)
  {
    dev_logger->warn("Repeat initializing logger.");
    return;
  }
  
  auto dev_console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_st>();
  dev_console_sink->set_level(console_level);
  auto user_console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_st>();
  user_console_sink->set_level(console_level);

  if (file_log)
  {
    auto dev_file_sink = std::make_shared<spdlog::sinks::basic_file_sink_st>(file_path, true);
    dev_file_sink->set_level(file_level);
    auto user_file_sink = std::make_shared<spdlog::sinks::basic_file_sink_st>(file_path, true);
    user_file_sink->set_level(file_level);
    dev_logger = std::make_shared<spdlog::logger>(
      spdlog::logger(
        "console and file logger",
        { dev_console_sink, dev_file_sink }
    ));
    user_logger = std::make_shared<spdlog::logger>(
      spdlog::logger(
        "console and file logger",
        { user_console_sink, user_file_sink }
    ));
  }
  else
  {
    dev_logger = std::make_shared<spdlog::logger>(
      "console logger", dev_console_sink);
    user_logger = std::make_shared<spdlog::logger>(
      "console logger", user_console_sink);
  }
  dev_logger->set_pattern("%^[%H:%M:%S][dev][%l]%v%$");
  user_logger->set_pattern("%^[%H:%M:%S][user][%l]%v%$");
  dev_logger->set_level(spdlog::level::trace);
  user_logger->set_level(spdlog::level::trace);

  sw = std::make_shared<spdlog::stopwatch>();
  dev_logger->info("create logger.");
  user_logger->info("create logger.");

  logger_initialized = true;
}

void Logger::updateFileLog(bool file_log, spdlog::level::level_enum file_level, const std::string& file_path)
{
  if (dev_logger->sinks().size() == 2)
    dev_logger->sinks().pop_back();
  if (user_logger->sinks().size() == 2)
    user_logger->sinks().pop_back();

  if (file_log)
  {
    auto dev_file_sink = std::make_shared<spdlog::sinks::basic_file_sink_st>(file_path, true);
    dev_file_sink->set_level(file_level);
    dev_logger->sinks().push_back(dev_file_sink);

    auto user_file_sink = std::make_shared<spdlog::sinks::basic_file_sink_st>(file_path, true);
    user_file_sink->set_level(file_level);
    user_logger->sinks().push_back(user_file_sink);
  }

  dev_logger->set_pattern("%^[%H:%M:%S][dev][%l]%v%$");
  user_logger->set_pattern("%^[%H:%M:%S][user][%l]%v%$");
}

void Logger::StopProgram()
{
  dev_logger->critical("trigger stopping program.");
  exit(1);
}
}// namespace Utils
}// namespace GCLF