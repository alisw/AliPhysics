#ifndef ALIITSMAPSDD_H
#define ALIITSMAPSDD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class for SDD maps used to correct for                        //
// voltage divider shape and doping fluctuations                 //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSsegmentationSDD.h"
#include<TNamed.h>
#include "AliLog.h"
class TH1F;
class TH2F;

class AliITSMapSDD : public TNamed {

 public:
  AliITSMapSDD();
  AliITSMapSDD(Char_t *mapname);
  virtual ~AliITSMapSDD(){};

  void SetMap(TH2F* hmap);
  Bool_t CheckBounds(Int_t iAn, Int_t iTb) const {
    if(iAn<0 || iAn>=fgkNAnodPts || iTb<0 || iTb >= fgkNDrifPts){ 
      AliWarning(Form("Cell out of bounds, anode=%d time-bin=%d",iAn,iTb));
      return kFALSE;
    }
    return kTRUE;
  }
  void SetCellContent(Int_t iAn, Int_t iTb, Float_t devMicron){
    if(CheckBounds(iAn,iTb)) fMap[iAn][iTb]=devMicron;
  }

  Float_t GetCellContent(Int_t iAn, Int_t iTb) const {
    if(CheckBounds(iAn,iTb)) return fMap[iAn][iTb];
    else return 0.;
  }
  Float_t GetCorrection(Float_t z, Float_t x, AliITSsegmentationSDD *seg);
  static Int_t GetNBinsAnode() {return fgkNAnodPts;}
  static Int_t GetNBinsDrift() {return fgkNDrifPts;}

  TH2F* GetMapHisto() const;
  TH1F* GetResidualDistr(Float_t dmin=-300., Float_t dmax=300.) const;

 protected:
  static const Int_t fgkNAnodPts = 256; // number of map points along anodes
  static const Int_t fgkNDrifPts = 72; // number of map points along anodes
  Float_t fMap[fgkNAnodPts][fgkNDrifPts];   // map of deviations

  ClassDef(AliITSMapSDD,1);
};
#endif
