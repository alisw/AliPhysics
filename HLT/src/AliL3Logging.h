// @(#) $Id$

#ifndef ALIL3LOGGING_H
#define ALIL3LOGGING_H

#include "AliL3RootTypes.h"
#include "AliL3StandardIncludes.h"

#ifdef use_logging
#include <MLUCLog.hpp>
#include <MLUCLogServer.hpp>

typedef MLUCLog AliL3Log;
typedef MLUCLogServer AliL3LogServer;
typedef MLUCDevNullLogServer AliL3DevNullLogServer;
typedef MLUCStdoutLogServer AliL3StdoutLogServer;
typedef MLUCStderrLogServer AliL3StderrLogServer;
typedef MLUCStreamLogServer AliL3StreamLogServer;

#else

class AliL3Log{
  public:
  enum TLogLevel { kNone = 0, kDebug= 0x01, kInformational = 0x02, kWarning = 0x04, kError = 0x08 , kFatal = 0x10, kPrimary = 0x80, kAll = 0x9F };
  enum TLogCmd { kEnd, kPrec, kHex, kDec };
};

#if __GNUC__ == 3
#define LOG( lvl, origin, keyword ) std::cerr<<"["<<origin<<": "<<keyword<<"] "
#define ENDLOG std::endl
#else
#define LOG( lvl, origin, keyword ) cerr<<"["<<origin<<": "<<keyword<<"] "
#define ENDLOG endl
#endif

#endif
#endif // ALIL3LOGGING_H


