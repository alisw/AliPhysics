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

#include "AliL3Log.h"
#endif /* use_logging */ 
#endif /* ALIL3LOGGING_H */


