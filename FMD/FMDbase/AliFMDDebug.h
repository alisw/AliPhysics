#ifndef ALIFMDDebug_H
#define ALIFMDDebug_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMD.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 17:59:37 2006
    @brief   Declaration of AliFMD detector driver 
*/
//____________________________________________________________________
//
// Some more clever declarations of Debug macros 
//
#include <AliLog.h>                        // ALILOG_H
#ifdef LOG_NO_DEBUG
#define AliFMDDebug(N, A) 
#else 
/** @defn AliFMDDebug 
    @param N Debug level - always evaluated 
    @param A Argument (including paranthesis) to Form - the message to
    print.  Note, that @a A should contain balanced paranthesis, like 
    @verbatim 
      AliFMDDebug(1, ("Failed to decode line %d of %s", line, filename));
    @endverbatim 
    The point is, if the current log level isn't high enough, as
    returned by the AliLog object, then we do not want to evalute the
    call to Form, since that is an expensive call.  We should always
    put macros like this into a @c do ... @c while loop, since that
    makes sure that evaluations are local, and that we can safely put
    a @c ; after the macro call.  Note, that @c do ... @c while loop
    and the call with extra paranthis, are an old tricks used by many
    C coders (see for example Bison, the Linux kernel, and the like). 
*/
#define AliFMDDebug(N, A) \
  do { \
    if (!AliLog::IsDebugEnabled() || \
      AliLog::GetDebugLevel(MODULENAME(), ClassName()) < N)  break; \
    AliLog::Debug(N, Form A, MODULENAME(), ClassName(), FUNCTIONNAME(), \
	  	  __FILE__, __LINE__); } while (false)
#endif

#endif
//
// EOF
//
