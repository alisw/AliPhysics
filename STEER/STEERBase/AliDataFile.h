#ifndef ALIDATAFILE_H
#define ALIDATAFILE_H

#include "TFile.h"

class AliDataFile {
  public:
    virtual ~AliDataFile() {};
    static std::string GetFileName(const std::string &url);
    static std::string GetFileNameOADB(const std::string &url);
    static TFile *Open(const std::string &url, Option_t *opts="") {
      return TFile::Open(GetFileName(url).c_str(), opts);
    }
    static TFile *OpenOADB(const std::string &url, Option_t *opts="") {
      return TFile::Open(GetFileNameOADB(url).c_str(), opts);
    }

  /// \cond CLASSIMP
  ClassDef(AliDataFile, 0);
  /// \endcond
};

#endif
