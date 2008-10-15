#include "LHAPDF/FortranWrappers.h"
#include <string>
#include <cstring>
#include <cassert>

using namespace std;

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "5.x.x"
#endif


extern "C" {

  #define fgetlhapdfversion FC_FUNC(getlhapdfversion, GETLHAPDFVERSION)
  void fgetlhapdfversion(char* fversion, int length) {
    string version = PACKAGE_VERSION;
    strncpy(fversion, version.c_str(), length);
    for (size_t i = strlen(fversion); i < (unsigned) length; ++i) {
      fversion[i] = ' ';
    }
  }

}
