// @(#) $Id$

#ifndef ALIHLTTPCLOGGING_H
#define ALIHLTTPCLOGGING_H

#if 1

#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCStandardIncludes.h"

#ifdef use_logging
#include <MLUCLog.hpp>
#include <MLUCLogServer.hpp>

typedef MLUCLog AliHLTTPCLog;
typedef MLUCLogServer AliHLTTPCLogServer;
typedef MLUCDevNullLogServer AliHLTTPCDevNullLogServer;
typedef MLUCStdoutLogServer AliHLTTPCStdoutLogServer;
typedef MLUCStderrLogServer AliHLTTPCStderrLogServer;
typedef MLUCStreamLogServer AliHLTTPCStreamLogServer;

#else
#include "AliHLTTPCLog.h"
#endif

#if 1
#else

class AliHLTTPCLog{
  public:
  enum TLogLevel { kNone = 0, kDebug= 0x01, kInformational = 0x02, kWarning = 0x04, kError = 0x08 , kFatal = 0x10, kPrimary = 0x80, kAll = 0x9F };
  enum TLogCmd { kEnd, kPrec, kHex, kDec };
};

#if __GNUC__>=3
#define LOG( lvl, origin, keyword ) std::cerr<<"["<<origin<<": "<<keyword<<"] "
#define ENDLOG std::endl
#else
#define LOG( lvl, origin, keyword ) cerr<<"["<<origin<<": "<<keyword<<"] "
#define ENDLOG endl
#endif

#endif


#else  // 0
#include "AliHLTLog.hpp"

typedef AliHLTLog AliHLTTPCLog;

#endif // 0


#endif // ALIHLTTPCLOGGING_H


