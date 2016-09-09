#ifndef STEER_UTILITIES_HELPERS_H
#define STEER_UTILITIES_HELPERS_H

#include <string>
#include <vector>
#include <cstdio>

struct MergeInput {
  MergeInput(const std::string &f, int order)
  : filename(f),
    order(order)
  {}
  std::string filename;
  size_t cost;
  int order;
};

struct MergeInputComparator {
  bool operator()(const MergeInput &a, const MergeInput &b) {
    return a.order < b.order;
  }
};

inline bool ends_with(std::string const & value, std::string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

struct MergeJob {
  MergeJob (const std::string &namePrefix, int maxJobs)
  :output(namePrefix, 0, namePrefix.size() - (ends_with(namePrefix, ".root") ? 5 : 0))
  {
    // In case we only have one file, simply add back the .root extension.
    // In case we have more than one file, add an incremental index to the filename.
    if (maxJobs == 1)
    {
      output += ".root";
      return;
    }
    static int count = 1;
    char buf[16];
    snprintf(buf, 10, "%i.root", count);
    // Add numeric suffix when splitting.
    output += buf;
    count++;
  }
  std::string output;
  std::vector<size_t> inputs;
};

#endif
