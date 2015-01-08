#ifndef ALIRAWEVENTTAG_H
#define ALIRAWEVENTTAG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// The AliRawEventTag class handles the raw-data event-oriented tag         //
// information. One object for each raw-data event is stored in a ROOT      //
// tree inside the file controled by AliTagDB class.                        //
// For the moment the tag information includes the raw-data event header +  //
// the raw-data file GUID and the event index.                              //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TString.h>

class AliRawEventHeaderBase;

class AliRawEventTag: public TObject {
 public:

  AliRawEventTag();
  virtual ~AliRawEventTag() {};

  void             SetHeader(AliRawEventHeaderBase *header) { fHeader = header; }
  void             SetGUID(const char *guid) { fGUID = guid; }
  void             SetEventNumber(Int_t event) { fEvent = event; }

  
  AliRawEventHeaderBase *GetHeader() const { return fHeader; }
  const char *     GetGUID() const { return fGUID.Data(); }
  Int_t            GetEventNumber() const { return fEvent; }

 private:

   AliRawEventTag(const AliRawEventTag& tag);
   AliRawEventTag& operator = (const AliRawEventTag& tag);

   AliRawEventHeaderBase *fHeader;    // raw data event header
   TString                fGUID;      // GUID of the raw data file
   Int_t                  fEvent;     // raw data event number inside the file

   ClassDef(AliRawEventTag,1)  // Raw data event tag
};

#endif
