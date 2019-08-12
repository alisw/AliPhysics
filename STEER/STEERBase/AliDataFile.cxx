#include "AliDataFile.h"
#include "AliLog.h"
#include "TSystem.h"
#include <iostream>

std::string AliDataFile::GetFileName(const std::string &url) {

  // Possible sources have the following priority:
  //
  // $ALICE_DATA/<url>
  // $ALICE_PHYSICS/<url>
  // $ALICE_ROOT/<url>
  // $ALICE_ROOT/../../../../data/analysis/YYYY/VVVV/<url>
  // /cvmfs/alice.cern.ch/data/prod/v5-XX-YY-01/<url>
  // /cvmfs/alice.cern.ch/data/analysis/YYYY/VVVV/<url>

  std::string buf, ver, year;
  const char *env;

  if ((env = gSystem->Getenv("ALIPHYSICS_VERSION"))) {
    buf = env;
    size_t cut = buf.rfind('-'); // remove build number
    size_t us = buf.find('_');
    if (cut != std::string::npos) {
      if (us != std::string::npos && us < cut) cut = us;
      ver = buf.substr(0, cut);
      if (buf.length() >= 12 && buf.substr(0, 4) == "vAN-") year = buf.substr(4, 4);
    }
  }

  if ((env = gSystem->Getenv("ALICE_DATA"))) {
    buf = std::string(env) + "/" + url;
    if (!gSystem->AccessPathName(buf.c_str())) {
      AliDebugClassF(2, "Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str());
      return buf;
    }
  }

  if ((env = gSystem->Getenv("ALICE_PHYSICS"))) {
    buf = std::string(env) + "/" + url;
    if (!gSystem->AccessPathName(buf.c_str())) {
      AliDebugClassF(2, "Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str());
      return buf;
    }
  }

  if ((env = gSystem->Getenv("ALICE_ROOT"))) {
    buf = std::string(env) + "/" + url;
    if (!gSystem->AccessPathName(buf.c_str())) {
      AliDebugClassF(2, "Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str());
      return buf;
    }
  }

  if (!ver.empty() && year.empty()) {
    // Production tag
    if ((env = gSystem->Getenv("ALICE_ROOT"))) {
      buf = std::string(env) + "/../../../../data/prod/" + "/" + ver + "/" + url;
      if (!gSystem->AccessPathName(buf.c_str())) {
        AliDebugClassF(2, "Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str());
        return buf;
      }
    }
    buf = "/cvmfs/alice.cern.ch/data/prod/" + ver + "/" + url;
    if (!gSystem->AccessPathName(buf.c_str())) {
      AliDebugClassF(2, "Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str());
      return buf;
    }
  }

  if (!ver.empty() && !year.empty()) {
    // Analysis tag
    if ((env = gSystem->Getenv("ALICE_ROOT"))) {
      buf = std::string(env) + "/../../../../data/analysis/" + year + "/" + ver + "/" + url;
      if (!gSystem->AccessPathName(buf.c_str())) {
        AliDebugClassF(2, "Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str());
        return buf;
      }
    }
    buf = "/cvmfs/alice.cern.ch/data/analysis/" + year + "/" + ver + "/" + url;
    if (!gSystem->AccessPathName(buf.c_str())) {
      AliDebugClassF(2, "Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str());
      return buf;
    }
  }

  AliErrorClassF("Data file not found on any location (use environment "
                 "variable $ALICE_DATA to manually specify one): \"%s\"", url.c_str());
  return "";
}

std::string AliDataFile::GetFileNameOADB(const std::string &url) {
  std::string buf;
  const char *env;

  if ((env = gSystem->Getenv("OADB_PATH"))) {
    buf = std::string(env) + "/" + url;
    if (!gSystem->AccessPathName(buf.c_str())) {
      AliDebugClassF(2, "Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str());
      return buf;
    }
  }

  return GetFileName("OADB/" + url);
}
