#ifndef ALIL3LOGGING_H
#define ALIL3LOGGING_H

#define use_logging

#ifdef use_logging

#include "/heim/franken/lib/MLUC/include/MLUCLog.hpp"
#include "/heim/franken/lib/MLUC/include/MLUCLogServer.hpp"

typedef MLUCLog AliL3Log;
typedef MLUCLogServer AliL3LogServer;
typedef MLUCDevNullLogServer AliL3DevNullLogServer;
typedef MLUCStdoutLogServer AliL3StdoutLogServer;
typedef MLUCStderrLogServer AliL3StderrLogServer;
typedef MLUCStreamLogServer AliL3StreamLogServer;

#else
#include <iostream.h>
class AliL3Log{
  public:
  enum TLogLevel { kNone = 0, kDebug= 0x01, kInformational = 0x02, kWarning = 0x04, kError = 0x08 , kFatal = 0x10, kPrimary = 0x80, kAll = 0x9F };
  enum TLogCmd { kEnd, kPrec, kHex, kDec };
};

#define LOG( lvl, origin, keyword ) cerr

#define ENDLOG endl

#endif
#endif // ALIL3LOGGING_H


