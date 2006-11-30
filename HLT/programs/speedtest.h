// $Id$

#include <iostream.h>
#define ALIL3LOGGING_H
class AliHLTLog{
  public:
  enum TLogLevel { kNone = 0, kDebug= 0x01, kInformational = 0x02, kWarning = 0x04, kError = 0x08 , kFatal = 0x10, kPrimary = 0x80, kAll = 0x9F };
  enum TLogCmd { kEnd, kPrec, kHex, kDec };
};
#define LOG( lvl, origin, keyword ) cerr
#define ENDLOG endl
#define no_root

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "AliHLTConfMapPoint.cxx"

double CpuTime()
{
  //Return the Cputime in seconds.

  return (double)(clock()) / CLOCKS_PER_SEC;
}


