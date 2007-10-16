// $Id$
// Original: AliHLTLog.cxx,v 1.1 2004/05/14 09:37:22 loizides 

// lagacy logging methods for HLT TPC code
// or                                                                        //
// refer to README to build package                                          //
// or                                                                        //
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt                          //

#include "AliHLTTPCLog.h"

AliHLTTPCLog::TLogLevel AliHLTTPCLog::fgLevel=AliHLTTPCLog::kNone;

const char* AliHLTTPCLog::kEnd = "";
const char* AliHLTTPCLog::kPrec = "";
const char* AliHLTTPCLog::kHex = "";
const char* AliHLTTPCLog::kDec = "";

const char* AliHLTTPCLog::fgKeyOrigin ="__origin";
const char* AliHLTTPCLog::fgKeyKeyword="__key";
const char* AliHLTTPCLog::fgKeyMessage="__message";

stringstream AliHLTTPCLog::fgStream;

AliHLTLogging AliHLTTPCLog::fgHLTLogging;

const char* AliHLTTPCLog::Flush()
{
  // see header file for class documentation
  int severity=0;
  string origin("");
  string keyword("");
  string iter;
  string message("");
  int scanStatus=0;
  fgStream >> severity;
  while (!fgStream.eof()) {
    fgStream >> iter;
    if (scanStatus==0 && iter.compare(fgKeyOrigin)==0) {
      // idicate scan of origin message
      scanStatus=1;
      continue;
    } else if (scanStatus==1 && iter.compare(fgKeyKeyword)==0) {
      // idicate scan of keyword message
      scanStatus=2;
      continue;
    } else if (scanStatus==2 && iter.compare(fgKeyMessage)==0) {
      scanStatus=3;
      continue;
    }
    switch (scanStatus) {
    case 1:
      if (!origin.empty()) origin+=" ";
      origin+=iter;
      break;
    case 2:
      if (!keyword.empty()) keyword+=" ";
      keyword+=iter;
      break;
    default:
      // if we have come here already for the first word, we don't
      // expect origin and keyword any more
      scanStatus=3;
      if (!message.empty()) message+=" ";
      message+=iter;
    }
  }
  
  // flush the string stream and send out through the logging system
  switch (severity) {
  case kDebug:
    fgHLTLogging.Logging(kHLTLogDebug, origin.c_str(), keyword.c_str(), message.c_str());
    break;
  case kInformational:
    fgHLTLogging.Logging(kHLTLogInfo, origin.c_str(), keyword.c_str(), message.c_str());
    break;
  case kWarning:
    fgHLTLogging.Logging(kHLTLogWarning, origin.c_str(), keyword.c_str(), message.c_str());
    break;
  case kError:
    fgHLTLogging.Logging(kHLTLogError, origin.c_str(), keyword.c_str(), message.c_str());
    break;
  case kFatal:
  case kPrimary:
    fgHLTLogging.Logging(kHLTLogFatal, origin.c_str(), keyword.c_str(), message.c_str());
    break;
  }
  fgStream.clear();
  string empty("");
  fgStream.str(empty);
  return "";
}

