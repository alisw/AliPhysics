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
#include "TObjArray.h"
#include "TString.h"

class TFile;

// --- AliRoot header files ---

class AliStream: public TNamed {

public:
  AliStream();
  AliStream(Option_t *option);
  AliStream(const AliStream& str);
  AliStream & operator=(const AliStream& str) 
    {str.Copy(*this); return (*this);}
  virtual ~AliStream();
  void AddFile(const char *fileName);
  Bool_t NextEventInStream(Int_t &eventNr);
  Bool_t OpenNextFile();
  Bool_t ImportgAlice();
  TFile* CurrentFile() { return fCurrentFile;}
  void ChangeMode(Option_t* option);     // reset READ or UPDATE mode
  Int_t GetNInputFiles() const {return fFileNames->GetLast()+1;}          
  TString GetFileName(const Int_t order) const; 
  
private:  
  void Copy(AliStream& str) const;

  Int_t fLastEventSerialNr;     // Serial number of the last event
  Int_t fLastEventNr;           // Number of the last event
  Int_t fCurrentFileIndex;      // Index of the current file
  Int_t fEvents;                //! nr. of events in the current file
  TString fMode;                // = 0 for READONLY, = 1 for READWRITE
  TFile *fCurrentFile;          //! pointer to current open file
  TObjArray * fFileNames;       // storage for TStrings with file names
  
  ClassDef(AliStream,1)
};

#endif // ALISTREAM_H
