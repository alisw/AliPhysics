#include "LHAPDF/FortranWrappers.h"
#include <iostream>
#include <string>
#include <cstring>
#include <cassert>
#include <unistd.h>
#include "binreloc.h"

using namespace std;


#define SIZE 499


extern "C" {


  #define fgetdirpath FC_FUNC(getdirpath, GETDIRPATH)
  void fgetdirpath(char*, int);


  #define fgetprefixpath FC_FUNC(getprefixpath, GETPREFIXPATH)
  void fgetprefixpath(char* prefixpath, int length) {
    BrInitError error;
    br_init_lib(&error);
    string prefixdir = br_find_prefix(DEFAULTPREFIXPATH);
    string test1 = prefixdir + "/share/lhapdf";
    if (access(test1.c_str(), R_OK) != 0) {
      prefixdir = DEFAULTPREFIXPATH;
    }
    assert(prefixdir.length() <= (unsigned) length);
    strncpy(prefixpath, prefixdir.c_str(), length);
    // Replace null-terminated string convention with Fortran "trailing spaces" convention:
    for (size_t i = strlen(prefixpath); i < (unsigned) length; ++i) {
      prefixpath[i] = ' ';
    }
  }


  #define fgetindexpath FC_FUNC(getindexpath, GETINDEXPATH)
  void fgetindexpath(char* indexpath, int length) {
    char tmp[SIZE+1];
    tmp[SIZE] = '\0';
    fgetdirpath(tmp, SIZE);
    //for (size_t i = 0; i < SIZE; ++i) {
    //  tmp[i] = ' ';
    //}
    for (int i = SIZE-1; i >= 0; --i) {
      if (tmp[i] != ' ') break;
      tmp[i] = '\0';
    }
    string try1(tmp), try2(tmp);
    try1 += "/PDFsets.index";
    try2 += ".index";
    if (access(try1.c_str(), R_OK) == 0) {
      assert(try1.length() <= (unsigned) length);
      strncpy(indexpath, try1.c_str(), length);
    } else {
      assert(try2.length() <= (unsigned) length);
      strncpy(indexpath, try2.c_str(), length);
    }
    // Replace null-terminated string convention with Fortran "trailing spaces" convention:
    for (size_t i = strlen(indexpath); i < (unsigned) length; ++i) {
      indexpath[i] = ' ';
    }
  }


  #define fgetdatapath FC_FUNC(getdatapath, GETDATAPATH)
  void fgetdatapath(char* datapath, int length) {
    BrInitError error;
    br_init_lib(&error);
    string sharedir = br_find_data_dir(DEFAULTLHAPATH);
    string tmp = sharedir + "/lhapdf/PDFsets";
    string test1 = tmp + "/cteq6.LHpdf";
    if (access(test1.c_str(), R_OK) != 0) {
      tmp = string(DEFAULTLHAPATH) + "/lhapdf/PDFsets";
    }
    assert(tmp.length() <= (unsigned) length);
    strncpy(datapath, tmp.c_str(), length);
    // Replace null-terminated string convention with Fortran "trailing spaces" convention:
    for (size_t i = strlen(datapath); i < (unsigned) length; ++i) {
      datapath[i] = ' ';
    }
  }


}
