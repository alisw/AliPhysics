// @(#) $Id$

#ifndef ALIL3LOGGING_H
#define ALIL3LOGGING_H

#include "AliHLTRootTypes.h"
#include "AliHLTStandardIncludes.h"

#ifdef use_logging
#include <MLUCLog.hpp>
#include <MLUCLogServer.hpp>

typedef MLUCLog AliHLTLog;
typedef MLUCLogServer AliHLTLogServer;
typedef MLUCDevNullLogServer AliHLTDevNullLogServer;
typedef MLUCStdoutLogServer AliHLTStdoutLogServer;
typedef MLUCStderrLogServer AliHLTStderrLogServer;
typedef MLUCStreamLogServer AliHLTStreamLogServer;

#else

#include "AliHLTLog.h"
#endif /* use_logging */ 
#endif /* ALIL3LOGGING_H */
