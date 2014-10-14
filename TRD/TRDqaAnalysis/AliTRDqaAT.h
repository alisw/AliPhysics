#ifndef ALITRDQAAT_H
#define ALITRDQAAT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id: AliTRDqaAT.h  $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  TRD analysis tools                                                    //
//                                                                        //
//  Author:                                                               //
//    Sylwester Radomski (radomski@physi.uni-heidelberg.de)               //              
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"  

class TH1D; 
class AliESDtrack;
class AliExternalTrackParam;

class AliTRDqaAT : public TObject {

 public:

  AliTRDqaAT();
  virtual ~AliTRDqaAT() {}
   
  static Int_t GetSector(const Double_t alpha);
  static Int_t GetStack(const AliExternalTrackParam *paramOut);
  static void  BuildRatio(TH1D *ratio, TH1D *histN, TH1D *histD);
  static void  FillStatus(TH1D *fStatusHist, UInt_t status);

  static void  PrintPID(const AliESDtrack *track);

  ClassDef(AliTRDqaAT, 0); // TRD analysis tools
};
#endif // ALITRDQAAT_H
