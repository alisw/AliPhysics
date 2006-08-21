// @(#) $Id$
// Original: AliL3Log.h,v 1.2 2004/06/11 16:06:33 loizides Exp $

#ifndef ALIHLTTPCLOG_H
#define ALIHLTTPCLOG_H

class AliHLTTPCLog {
  public:
  enum TLogLevel { kNone = 0, kDebug= 0x01, kInformational = 0x02, kWarning = 0x04, kError = 0x08 , kFatal = 0x10, kPrimary = 0x80, kAll = 0x9F };
  enum TLogCmd { kEnd, kPrec, kHex, kDec };
  static TLogLevel fgLevel;
};

#if __GNUC__ == 3
#define LOG( lvl, origin, keyword ) \
 if (lvl>=AliHLTTPCLog::fgLevel) std::cerr<<"["<<origin<<": "<<keyword<<"] "
#define ENDLOG std::endl
#else
#define LOG( lvl, origin, keyword ) \
 if (lvl>=AliHLTTPCLog::fgLevel) cerr<<"["<<origin<<": "<<keyword<<"] "
#define ENDLOG endl
#endif /* __GNUC__ */
#endif /* ALIHLTTPCLOG_H */
