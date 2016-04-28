#include "AliFileUtilities.h"

#include "TFileCollection.h"
#include "TFileInfo.h"
#include "THashList.h"
#include <cstdio>
#include <vector>
#include <string>
#include <cstring>
#include <glob.h>
#include <cstdlib>
#include <iostream>
#include <dirent.h>
#include <sys/stat.h>
#include <cassert>

namespace AliFileUtilities {

// a little helper function to expand filepatterns with wildcards
std::vector<std::string> glob(const std::string &pat) {
  using namespace std;
  glob_t glob_result;
  glob(pat.c_str(), GLOB_TILDE, NULL, &glob_result);
  vector<string> ret;
  for (unsigned int i = 0; i < glob_result.gl_pathc; ++i) {
    ret.push_back(string(glob_result.gl_pathv[i]));
  }
  globfree(&glob_result);
  return ret;
}

size_t RemoveLocalFiles(const char *pattern) {
  std::vector<std::string> matchingfiles(glob(pattern));
  if (matchingfiles.size() == 0) {
    // no matching file found
    return 0;
  }

  size_t counter = 0;
  for (auto v : matchingfiles) {
    if(RemoveLocalFile(v.c_str()))
      counter++;
  }
  return counter;
}

size_t CountLocalFiles(const char *pattern) {
  std::vector<std::string> matchingfiles(glob(pattern));
  return matchingfiles.size();
}

bool isDotOrDotDot(const char *name) {
  if (strcmp(".", name) == 0 || strcmp("..", name) == 0)
    return true;
  return false;
}

bool IsDirectory(const char *path) {
  struct stat st_buf;
  assert(!isDotOrDotDot(path));
  stat(path, &st_buf);
  return S_ISDIR(st_buf.st_mode) && !S_ISLNK(st_buf.st_mode);
}

// internal implementation function
size_t remove_all_aux(const char *name) {
  size_t counter = 0;
  DIR *dir = opendir(name);
  if (dir!=nullptr){
    struct dirent *ent;
    ent = readdir(dir);

    while (ent != nullptr) {
      char *entry_name = ent->d_name;
      // we do not list nor follow "." or ".."
      if (isDotOrDotDot(entry_name)) {
        ent = readdir(dir);
        continue;
      }
      const size_t S = strlen(name) + strlen(entry_name) + 2;
      char next[S];
      sprintf(next, "%s/%s", name, entry_name);
      // if this is a directory and not a symbolic link; we recurse
      if (IsDirectory(next)) {
        counter += remove_all_aux(next);
      }
      // 0 is success
      if (std::remove(next) == 0){
        counter++;
      }
      ent = readdir(dir);
    }
    closedir(dir);
  }  
  return counter;
}

size_t Remove_All(const char *name) {
  if (!IsDirectory(name))
    return 0;
  size_t c = remove_all_aux(name);
  // 0 is success
  if (std::remove(name) == 0) {
    c++;
  }
  return c;
}
}
