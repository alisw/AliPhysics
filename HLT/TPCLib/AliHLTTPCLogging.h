// @(#) $Id$
// Original: AliHLTLogging.h,v 1.11 2004/10/09 11:28:31 loizides 

#ifndef ALIHLTTPCLOGGING_H
#define ALIHLTTPCLOGGING_H

#include "AliHLTTPCRootTypes.h"
#include "AliHLTStdIncludes.h"

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
#endif /* use_logging */ 
#endif /* ALIHLTTPCLOGGING_H */
