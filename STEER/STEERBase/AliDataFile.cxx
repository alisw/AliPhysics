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
  // /cvmfs/alice.cern.ch/data/analysis/YYYY/VVVV/<url>
  // root://eospublic.cern.ch//eos/experiment/alice/analysis-data/<url>

  std::string buf, ver, year;
  const char *env;

  if (env = gSystem->Getenv("ALIPHYSICS_VERSION")) {
    buf = env;
    if (buf.length() >= 12) {
      ver = buf.substr(0, 12);
      year = buf.substr(4, 4);
    }
  }

  if (env = gSystem->Getenv("ALICE_DATA")) {
    buf = std::string(env) + "/" + url;
    if (!gSystem->AccessPathName(buf.c_str())) {
      AliDebugClass(2, TString::Format("Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str()));
      return buf;
    }
  }

  if (env = gSystem->Getenv("ALICE_PHYSICS")) {
    buf = std::string(env) + "/" + url;
    if (!gSystem->AccessPathName(buf.c_str())) {
      AliDebugClass(2, TString::Format("Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str()));
      return buf;
    }
  }

  if (env = gSystem->Getenv("ALICE_ROOT")) {
    buf = std::string(env) + "/" + url;
    if (!gSystem->AccessPathName(buf.c_str())) {
      AliDebugClass(2, TString::Format("Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str()));
      return buf;
    }
  }

  if (!ver.empty() && !year.empty()) {
    if (env = gSystem->Getenv("ALICE_ROOT")) {
      buf = std::string(env) + "/../../../../data/analysis/" + year + "/" + ver + "/" + url;
      if (!gSystem->AccessPathName(buf.c_str())) {
        AliDebugClass(2, TString::Format("Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str()));
        return buf;
      }
    }
    buf = "/cvmfs/alice.cern.ch/data/analysis/" + year + "/" + ver + "/" + url;
    if (!gSystem->AccessPathName(buf.c_str())) {
      AliDebugClass(2, TString::Format("Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str()));
      return buf;
    }
  }

  buf = "root://eospublic.cern.ch//eos/experiment/alice/analysis-data/";
  buf += url;
  return buf;
}

std::string AliDataFile::GetFileNameOADB(const std::string &url) {
  std::string buf;
  const char *env;

  if (env = gSystem->Getenv("OADB_PATH")) {
    buf = std::string(env) + "/" + url;
    if (!gSystem->AccessPathName(buf.c_str())) {
      AliDebugClass(2, TString::Format("Using data file \"%s\" from \"%s\"", url.c_str(), buf.c_str()));
      return buf;
    }
  }

  return GetFileName("OADB/" + url);
}
