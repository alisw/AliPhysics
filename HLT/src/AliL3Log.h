// @(#) $Id$

#ifndef ALIL3LOG_H
#define ALIL3LOG_H

class AliL3Log{
  public:
  enum TLogLevel { kNone = 0, kDebug= 0x01, kInformational = 0x02, kWarning = 0x04, kError = 0x08 , kFatal = 0x10, kPrimary = 0x80, kAll = 0x9F };
  enum TLogCmd { kEnd, kPrec, kHex, kDec };
  static TLogLevel fgLevel;
};

#if __GNUC__ == 3
#define LOG( lvl, origin, keyword ) \
 if (lvl>=AliL3Log::fgLevel) std::cerr<<"["<<origin<<": "<<keyword<<"] "
//#define LOG( lvl, origin, keyword ) std::cerr<<"["<<origin<<": "<<keyword<<"] "
#define ENDLOG std::endl
#else
//#define LOG( lvl, origin, keyword ) cerr<<"["<<origin<<": "<<keyword<<"] "
//#define ENDLOG endl
#define LOG( lvl, origin, keyword ) \
 if (lvl>=AliL3Log::fgLevel) cerr<<"["<<origin<<": "<<keyword<<"] "
#define ENDLOG endl
#endif /* __GNUC__ */
#endif /* ALIL3LOG_H */
