// Needed to have working asserts also when building in release mode
#ifdef NDEBUG
#undef NDEBUG
#endif

#include "helpers.h"
#include <cassert>
#include <iostream>

int main() {
  assert(ends_with("foobar", "bar"));
  assert(!ends_with("foobar", "foo"));

  MergeInput a("file0.root", 0);
  MergeInput b("file1.root", 1);
  assert(a.filename == "file0.root");
  assert(a.order == 0);

  MergeInputComparator comparator;
  assert(comparator(a, b) == true);

  MergeJob job1("foo.root", 1);
  assert(job1.output == "foo.root");
  MergeJob job2("bar", 1);
  assert(job2.output == "bar.root");
  MergeJob job3("bar", 2);
  assert(job3.output == "bar1.root");
  MergeJob job4("bar.root", 2);
  assert(job4.output == "bar2.root");
}
