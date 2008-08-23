#ifndef ALIITSCOMPRESSRAWDATASDD_H
#define ALIITSCOMPRESSRAWDATASDD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include<TObject.h>
#include<TString.h>
///////////////////////////////////////////////////////////////////
//                                                               //
// Class to decode the SDD Raw Data from the CarlosRX format to  //
// a compressed format consisting in a word of 32 bit per cell   //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

class AliITSCompressRawDataSDD : public TObject {

 public:
  AliITSCompressRawDataSDD(TString filename);
  void SetEventRange(Int_t first, Int_t last){
    fEventRange=kTRUE;
    fFirstEvent=first;
    fLastEvent=last;
  }
  void Compress();

 protected:
  TString fNameFile;    // name of the raw data file
  Bool_t  fEventRange;  // flag to select a range of events
  Int_t   fFirstEvent;  // first event (used only if fEventRange==kTRUE)
  Int_t   fLastEvent;  // first event (used only if fEventRange==kTRUE)

  ClassDef(AliITSCompressRawDataSDD, 0)
};

#endif
