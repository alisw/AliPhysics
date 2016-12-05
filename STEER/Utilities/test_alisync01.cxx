#include "alisync.h"
#include <algorithm>
#include <string>
#include <cassert>
#include <iostream>
#include <cstddef>

int
main(int argc, char **argv) {
  assert(trimTrailingSlashes(std::string(".")) == ".");
  assert(trimTrailingSlashes(std::string("/")) == "/");
  assert(trimTrailingSlashes(std::string("/a")) == "/a");
  assert(trimTrailingSlashes(std::string("foobar")) == "foobar");
  assert(trimTrailingSlashes(std::string("b/")) == "b");
  assert(trimTrailingSlashes(std::string("b/a")) == "b/a");
  assert(trimTrailingSlashes(std::string("/b/a")) == "/b/a");
  assert(trimTrailingSlashes(std::string("/b/a/")) == "/b/a");
  assert(trimTrailingSlashes(std::string("/b/a/")) == "/b/a");
  assert(normalizePath("/a/./b") == "/a/b");
  assert(normalizePath("/a/b") == "/a/b");
  assert(normalizePath("/a//b") == "/a//b");
  assert(normalizePath("/././a///b.s") == "/a///b.s");
}
