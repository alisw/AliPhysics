#include "LHAPDF/FortranWrappers.h"
#include "LHAPDF/LHAPDFConfig.h"
#include <string>
#include <cstring>
#include <cassert>

using namespace std;

#ifndef LHAPDF_VERSION
#define LHAPDF_VERSION "5.x.x"
#endif


extern "C" {

  #define fgetlhapdfversion FC_FUNC(getlhapdfversion, GETLHAPDFVERSION)
  void fgetlhapdfversion(char* fversion, int length) {
    string version = LHAPDF_VERSION;
    strncpy(fversion, version.c_str(), length);
    for (size_t i = strlen(fversion); i < (unsigned) length; ++i) {
      fversion[i] = ' ';
    }
  }

}
