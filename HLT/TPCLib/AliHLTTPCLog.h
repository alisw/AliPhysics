// @(#) $Id$
// Original: AliHLTLog.h,v 1.2 2004/06/11 16:06:33 loizides Exp $

#ifndef ALIHLTTPCLOG_H
#define ALIHLTTPCLOG_H

#ifndef __CINT__
#include <bits/ios_base.h>
#endif

class AliHLTTPCLog {
  public:
  enum TLogLevel { kNone = 0, kDebug= 0x01, kInformational = 0x02, kWarning = 0x04, kError = 0x08 , kFatal = 0x10, kPrimary = 0x80, kAll = 0x9F };
  //enum TLogCmd { kEnd, kPrec, kHex=std::ios_base::hex, kDec=std::ios_base::dec };
  static const char* kEnd;
  static const char* kPrec;
  static const char* kHex;
  static const char* kDec;
/*   static const std::ios_base::fmtflags kHex; */
/*   static const std::ios_base::fmtflags kDec; */
  static TLogLevel fgLevel;
};

#if __GNUC__ >= 3
#define LOG( lvl, origin, keyword ) \
 if (lvl>=AliHLTTPCLog::fgLevel) std::cerr<<"["<<origin<<": "<<keyword<<"] "
#define ENDLOG std::endl
#else
#error old gcc!
#define LOG( lvl, origin, keyword ) \
 if (lvl>=AliHLTTPCLog::fgLevel) cerr<<"["<<origin<<": "<<keyword<<"] "
#define ENDLOG endl
#endif /* __GNUC__ */
#endif /* ALIHLTTPCLOG_H */
