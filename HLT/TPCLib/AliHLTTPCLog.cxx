// $Id$
// Original: AliHLTLog.cxx,v 1.1 2004/05/14 09:37:22 loizides 

#ifndef use_logging

#include "AliHLTTPCLogging.h"

AliHLTTPCLog::TLogLevel AliHLTTPCLog::fgLevel=AliHLTTPCLog::kNone;

const char* AliHLTTPCLog::kEnd = "";
const char* AliHLTTPCLog::kPrec = "";
const char* AliHLTTPCLog::kHex = "";
const char* AliHLTTPCLog::kDec = "";
// const std::ios_base::fmtflags AliHLTTPCLog::kHex = std::ios_base::hex;
// const std::ios_base::fmtflags AliHLTTPCLog::kDec = std::ios_base::dec;


#endif
