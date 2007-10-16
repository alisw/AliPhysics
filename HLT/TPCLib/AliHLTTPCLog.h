// @(#) $Id$
// Original: AliHLTLog.h,v 1.2 2004/06/11 16:06:33 loizides Exp $

#ifndef ALIHLTTPCLOG_H
#define ALIHLTTPCLOG_H

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <sstream>
#include <iostream>
#include "AliHLTLogging.h"

using namespace std;

/**
 * @class AliHLTTPCLog
 * This class implements the old HLT TPC logging mechanism.
 * Logging messages are now forwarded to the HLT logging system
 * \em Note: the old LOG and ENDLOG macros should be used any longer,
 * use the HLT logging macros or AliRoot logging macros instead. 
 * @see AliHLTLogging
 */
class AliHLTTPCLog  {
  public:
  enum TLogLevel { kNone = 0, kDebug= 0x01, kInformational = 0x02, kWarning = 0x04, kError = 0x08 , kFatal = 0x10, kPrimary = 0x80, kAll = 0x9F };

  /** not used */
  static const char* kEnd;                                         //! transient
  /** not used */
  static const char* kPrec;                                        //! transient
  /** stream manipulator for hex output, but empty in the implementation */
  static const char* kHex;                                         //! transient
  /** stream manipulator for decimal output, but empty in the implementation */
  static const char* kDec;                                         //! transient

  /** the logging filter */
  static TLogLevel fgLevel;                                        // see above

  /** key to indicate the origin part */
  static const char* fgKeyOrigin;                                  //! transient
  /** key to indicate the keyword part */
  static const char* fgKeyKeyword;                                 //! transient
  /** key to indicate the message part */
  static const char* fgKeyMessage;                                 //! transient

  /** a stringstream to receive the output */
  static stringstream fgStream;                                    // see above
  /** HLT logging instance */
  static AliHLTLogging fgHLTLogging;                               // see above

  /**
   * Flush the stringstream and print output to the HLT logging system.
   * The attributes are set before the message is streamed into the
   * stringstream.<br>
   * The LOG macro sets the attributes from the macro arguments and provides
   * the stringstream.<br>
   * The ENDLOG macro calls the Flush method after the message was streamed
   * into the stringstream.
   */
  static const char* Flush();

 private:
  /** copy constructor prohibited */
  AliHLTTPCLog(const AliHLTTPCLog&);
  /** assignment operator prohibited */
  AliHLTTPCLog& operator=(const AliHLTTPCLog&);

};

/** LOG macro to be used by the TPC code 
 * \em Note: this macro should be used any longer 
 */
#define LOG( lvl, origin, keyword ) \
 if (lvl>=AliHLTTPCLog::fgLevel) AliHLTTPCLog::fgStream << lvl \
                           << " " << AliHLTTPCLog::fgKeyOrigin  << " " << origin \
                           << " " << AliHLTTPCLog::fgKeyKeyword << " " << keyword \
			   << " " << AliHLTTPCLog::fgKeyMessage << " "

/** ENDLOG macro calls the Flush method 
 * \em Note: this macro should be used any longer 
 */
#define ENDLOG AliHLTTPCLog::Flush()

#endif /* ALIHLTTPCLOG_H */
