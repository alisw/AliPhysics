#ifndef ALIITSCORRMAPSDD_H
#define ALIITSCORRMAPSDD_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Mother class for SDD maps used to correct for                 //
// voltage divider shape and doping fluctuations                 //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include "AliITSsegmentationSDD.h"
#include<TNamed.h>
#include "AliLog.h"
class TH1F;
class TH2F;

class AliITSCorrMapSDD : public TNamed {

 public:
  // maps produced from residuals stores Xtrue-Xmeasured in drift lenght.
  // Since the map computes xtrue-xmeas in drift coordinate
  // and the difference is added to measured local coordinate, we have
  // Left:  Xmeas_loc_corr = Xmeas_loc + (Xtrue-Xmeas)_loc = Xmeas_loc - (Xtrue-Xmeas)_drift
  // Right: Xmeas_loc_corr = Xmeas_loc + (Xtrue-Xmeas)_loc = Xmeas_loc + (Xtrue-Xmeas)_drift
  // hence, for the left side the sign of the correction need to inverted
  enum {kLeftInvBit = BIT(14)};   
  AliITSCorrMapSDD();
  AliITSCorrMapSDD(Char_t *mapname);
  virtual ~AliITSCorrMapSDD(){};
  //
  void   SetInversionBit(Bool_t v=kTRUE)           {SetBit(kLeftInvBit,v);}
  Bool_t GetInversionBit()               const     {return TestBit(kLeftInvBit);}
  //
  Int_t GetNBinsAnode() const {return fNAnodePts;}
  Int_t GetNBinsDrift() const {return fNDriftPts;}
  void SetNBinsAnode(Int_t nbins) {
    if(nbins<=kMaxNAnodePts) fNAnodePts=nbins;
    else AliError(Form("Max. number of anode bins = %d",kMaxNAnodePts));
  }
  void SetNBinsDrift(Int_t nbins) {
    if(nbins<=kMaxNDriftPts) fNDriftPts=nbins;
    else AliError(Form("Max. number of drift bins = %d",kMaxNDriftPts));
  }

  Bool_t CheckAnodeBounds(Int_t iAn) const {
    if(iAn<0 || iAn>=fNAnodePts)return kFALSE;
    else return kTRUE;
  }
  Bool_t CheckDriftBounds(Int_t iTb) const {
    if(iTb<0 || iTb >= fNDriftPts)return kFALSE;
    else return kTRUE;
  }

  virtual void Set1DMap(TH1F* /*hmap*/){
    AliError("Not implemented");
  }
  virtual void Set2DMap(TH2F* /*hmap*/){
    AliError("Not implemented");
  }

  virtual void ResetMap(){
    AliError("Not implemented");
  }
  virtual void SetCellContent(Int_t /*iAn*/, Int_t /*iTb*/, Float_t /*devMicron*/){
    AliError("Not implemented");
  }
  virtual Float_t GetCellContent(Int_t /*iAn*/, Int_t /*iTb*/) const {
    AliError("Not implemented");
    return -99999.;
  }

  void    ComputeGridPoints(Float_t z, Float_t x, AliITSsegmentationSDD *seg, Bool_t isReco=kTRUE);
  Float_t GetCorrection(Float_t z, Float_t x, AliITSsegmentationSDD *seg);
  Float_t GetShiftForSimulation(Float_t z, Float_t x, AliITSsegmentationSDD *seg);
  TH2F* GetMapHisto() const;
  TH1F* GetMapProfile() const;
  TH1F* GetResidualDistr(Float_t dmin=-300., Float_t dmax=300.) const;


 protected:

  enum {kMaxNAnodePts=256};// max number of map points along anodes
  enum {kMaxNDriftPts=291};// max number of map points along drift

  static const Int_t fgkNAnodePtsDefault; // default value for fNAnodePts
  static const Int_t fgkNDriftPtsDefault; // default value for fNDriftPts
  Int_t fNAnodePts; // number of map points along anodes
  Int_t fNDriftPts; // number of map points along anodes

  Float_t fXt1;   // true coordinate in lower grid point
  Float_t fXt2;   // true coordinate in upper grid point
  Float_t fXm1;   // measured coordinate in lower grid point
  Float_t fXm2;   // measured coordinate in upper grid point
  Float_t fDrLen; // drift length

  ClassDef(AliITSCorrMapSDD,2);
};
#endif

