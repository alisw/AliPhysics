#ifndef ALITRDQARECPOINTS_H
#define ALITRDQARECPOINTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDqaRecPoints.h 23387 2008-01-17 17:25:16Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  QA of black events                                                    //
//                                                                        //
//  Author:                                                               //
//    Sylwester Radomski (radomski@physi.uni-heidelberg.de)               //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TH1;
class TH1D;
class TH3D;
class TObjArray;

class AliTRDqaRecPoints : public TObject {

 public:
  
  AliTRDqaRecPoints();
  AliTRDqaRecPoints(const AliTRDqaRecPoints &qa);
  ~AliTRDqaRecPoints() {}
  AliTRDqaRecPoints& operator = (const AliTRDqaRecPoints& /*qa*/) { return *this; };
 
  void Init();
  void Reset() {}
  void AddEvent(TTree *tree);
  void Process(const char* filename);
    
  void SetNPad(Int_t nPad) {fnPad = nPad;}
  void CreateRef(Int_t ref) {fRef = ref;}

 private:
  
  Int_t fnEvents;         // number of events processed  
  TObjArray *fHist;       // histograms

  TH2D *fRefHist[540];    // reference histograms
  
  Int_t fnPad;            // something
  Int_t fRef;             // something else

  ClassDef(AliTRDqaRecPoints,0) // QA for black events  

};
#endif
