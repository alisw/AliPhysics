#ifndef ALISTREAM_H
#define ALISTREAM_H
/* Copyright(c) 1998-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////////
//
//  Class to handle files on IO
//                  
//  Author: Jiri Chudoba (CERN), 2001
//
////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TNamed.h"

// --- AliRoot header files ---
//#include <TString.h>

#include "TObjArray.h"
class TFile;

class TString;

class AliStream: public TNamed {

public:
  AliStream();
  AliStream(const char* foldername, Option_t *optioneventfoldername);
  virtual ~AliStream();

  void       AddFile(const char *fileName);
  Bool_t     NextEventInStream();
  Bool_t     OpenNextFile();//returns kFALSE in case of failure
  Bool_t     ImportgAlice();
  void       ChangeMode(Option_t* option);     // reset READ or UPDATE mode
 
  const TString& GetFolderName() const{return fEventFolderName;}
  Int_t GetNInputFiles() const {return fFileNames->GetLast()+1;}
  TString GetFileName(Int_t order) const;
  void SetFolderName(const TString name) { fEventFolderName = name ; }    
private:  
  Int_t      fLastEventSerialNr;
  Int_t      fLastEventNr;
  Int_t      fCurrentFileIndex;
  Int_t      fEvents;                //! nr. of events in the current file
  TString    fMode;                  // = 0 for READONLY, = 1 for READWRITE
  TObjArray* fFileNames;
  
  TString fEventFolderName; //Name of the folder where data for this stram will be mounted
  
  ClassDef(AliStream,1)
};

#endif // ALISTREAM_H
