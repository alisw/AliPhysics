#ifndef ALIFILEUTILITIES_H
#define ALIFILEUTILITIES_H

#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>

// put a standard header
// History: 25.4.2016 -- initial version (sandro.wenzel@cern.ch)

// A namespace encapsulating some standard static "filesystem" like functions
// The purpose of this class is to prevent (very expensive) calls
// to the OS via system or gSystem->Exec("..")
namespace AliFileUtilities {

// removes a specifc file specified by filename
// returns true of successfull (false if failure)
inline bool RemoveLocalFile(const char *filename) {
  return (std::remove(filename) == 0);
}

// brief: removes (multiple) files matching a pattern from the filesystem
// returns the number of successfull remove operations
// example RemoveLocalFile("/tmp/all*.root")
// note that a list of filenames will not work
size_t RemoveLocalFiles(const char * /*pattern*/);

// counts the local files matching a pattern
size_t CountLocalFiles(const char * /*pattern*/);

// retrieve the local files matching a pattern in a vevtor
void GetLocalFiles(const char *, std::vector<std::string> &);

// remove a directory given by paramer name and all its including things
// mimics behaviour of std::remove_all available in C++17
// returns the number of objects deleted
size_t Remove_All(const char * /*name*/);

}

#endif
