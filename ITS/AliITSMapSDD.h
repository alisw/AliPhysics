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

#include<TNamed.h>

class TH1F;
class TH2F;

class AliITSMapSDD : public TNamed {

 public:
  AliITSMapSDD();
  AliITSMapSDD(Char_t *mapname);
  virtual ~AliITSMapSDD(){};

  void SetMap(TH2F* hmap);
  void SetCellContent(Int_t iAn, Int_t iTb, Float_t devMicron){
    fMap[iAn][iTb]=devMicron;
  }
  Float_t GetCellContent(Int_t iAn, Int_t iTb) const {return fMap[iAn][iTb];}
  TH2F* GetMapHisto() const;
  TH1F* GetResidualDistr(Float_t dmin=-300., Float_t dmax=300.) const;

 protected:
  static const Int_t fgkNAnodPts = 256; // number of map points along anodes
  static const Int_t fgkNDrifPts = 72; // number of map points along anodes
  Float_t fMap[fgkNAnodPts][fgkNDrifPts];   // map of deviations

  ClassDef(AliITSMapSDD,1);
};
#endif
