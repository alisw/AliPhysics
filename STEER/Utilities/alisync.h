#ifndef STEER_Utilities_alisync_H
#define STEER_Utilities_alisync_H
#include <string>
#include <sys/types.h>
#include <time.h>
#include <algorithm>
#include <functional>


// Helpers to trim whitespace from strings. trim from start.
static inline std::string &ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
  return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
  return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
  return ltrim(rtrim(s));
}

struct MatchPathSeparator
{
  bool operator()( char ch ) const
  {
    return ch == '/';
  }
};

struct NotMatchPathSeparator
{
  bool operator()( char ch ) const
  {
    return ch != '/';
  }
};

static std::string trimTrailingSlashes(const std::string &s) {
  if (s == "/")
    return s;

  std::string o = s;
  // Find the last slash on the right boundary of the string
  // - if found trim all the trailing slashes
  // - if not found contribute to memory churn.
  std::string::reverse_iterator lsi = std::find_if(o.rbegin(),
                                                   o.rend(),
                                                   NotMatchPathSeparator());
  if (lsi != o.rbegin())
    o.erase(lsi.base());
  return o;
}

std::string normalizePath(const std::string &s) {
  std::string str = s;
  size_t index = 0;
  while (true) {
    index = str.find("/./", index);
    if (index == std::string::npos) break;

    /* Make the replacement. */
    str.erase(index, 2);
  }
  return str;
}


struct JobInfo {
  std::string filename;
  time_t startTime;
  pid_t pid;
  int exitCode;
  int retries;
  bool running;
};

struct SamePidAs {
  SamePidAs(pid_t pid) : pid(pid) {}
  pid_t pid;
  bool operator()(const JobInfo &info) { return info.pid == this->pid;}
};

#endif // STEER_Utilities_alisync_H
