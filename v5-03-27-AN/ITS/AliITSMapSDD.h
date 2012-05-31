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
#include "AliITSCorrMapSDD.h"
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
    if(CheckBounds(iAn,iTb)) fMap[iAn][iTb]=(Short_t)(devMicron*10.+0.5);
  }

  Float_t GetCellContent(Int_t iAn, Int_t iTb) const {
    if(CheckBounds(iAn,iTb)) return (Float_t)fMap[iAn][iTb]/10.;
    else return 0.;
  }
  Float_t GetCorrection(Float_t z, Float_t x, AliITSsegmentationSDD *seg);
  static Int_t GetNBinsAnode() {return fgkNAnodPts;}
  static Int_t GetNBinsDrift() {return fgkNDrifPts;}
  AliITSCorrMapSDD* ConvertToNewFormat() const;

  TH2F* GetMapHisto() const;
  TH1F* GetResidualDistr(Float_t dmin=-300., Float_t dmax=300.) const;

 protected:
  static const Int_t fgkNAnodPts = 256; // number of map points along anodes
  static const Int_t fgkNDrifPts = 72; // number of map points along anodes
  Short_t fMap[fgkNAnodPts][fgkNDrifPts];   // map of deviations
                                            // stored as Short_t: integer 
                                            // values from -32000 to 32000
                                            // in the range -3.2 - 3.2 mm

  ClassDef(AliITSMapSDD,2);
};
#endif
