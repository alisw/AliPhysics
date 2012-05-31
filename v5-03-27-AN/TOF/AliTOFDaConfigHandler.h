#ifndef ALITOFDACONFIGHANDLER_H
#define ALITOFDACONFIGHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * * See cxx source for full Copyright notice */
/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The SAX XML file handler used by the TOF DA for PHYSICS runs          //
//  to get the necessary flags to run (e.g. debug flag)                   //
//                                                                        //
//  Chiara.Zampolli (Chiara.Zampolli@cern.ch)                             //
//  Roberto Preghenella (R+) (preghenella@bo.infn.it)                     //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
class TString;

class AliTOFDaConfigHandler : public TObject {

public:

  AliTOFDaConfigHandler();
  AliTOFDaConfigHandler(const AliTOFDaConfigHandler &sh);
  virtual ~AliTOFDaConfigHandler();
  AliTOFDaConfigHandler &operator=(const AliTOFDaConfigHandler &sh);

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

  Int_t GetMeanMultiplicity() const {return fMeanMultiplicity;};
  Int_t GetMaxHits() const {return fMaxHits;};

 private:

  Int_t fMeanMultiplicity; // mean multiplicity
  Int_t fMaxHits; // max hits 

  ClassDef(AliTOFDaConfigHandler,1);  
};
#endif

