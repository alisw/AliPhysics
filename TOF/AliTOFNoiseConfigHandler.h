#ifndef ALITOFNOISECONFIGHANDLER_H
#define ALITOFNOISECONFIGHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * * See cxx source for full Copyright notice */
/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The SAX XML file handler used by the TOF noiseDA                      //
//  to get the necessary flags to run (e.g. debug flag)                   //
//                                                                        //
//  Chiara.Zampolli (Chiara.Zampolli@cern.ch)                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
class TString;

class AliTOFNoiseConfigHandler : public TObject {

public:

  AliTOFNoiseConfigHandler();
  AliTOFNoiseConfigHandler(const AliTOFNoiseConfigHandler &sh);
  virtual ~AliTOFNoiseConfigHandler();
  AliTOFNoiseConfigHandler &operator=(const AliTOFNoiseConfigHandler &sh);

  // functions to interface to TSAXHandler
  void          OnStartDocument();
  void          OnEndDocument();
  void          OnStartElement(const char *name, const TList *attributes);
  void          OnEndElement(const char *name);
  void          OnCharacters(const char *name);
  void          OnComment(const char *name);
  void          OnWarning(const char *name);
  void          OnError(const char *name);
  void          OnFatalError(const char *name);
  void          OnCdataBlock(const char *name, Int_t len);

  Int_t GetDebugFlag() const {return fDebugFlag;}

 private:
  Int_t fDebugFlag;          // debug flag: 0-->off, 1-->on

  ClassDef(AliTOFNoiseConfigHandler,0);    // The XML file handler for the OCDB
};
#endif

