#ifndef ALITOFRESPONSEPARAMS_H
#define ALITOFRESPONSEPARAMS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/* $Id$ */

// *
// *
// *
// * this class defines the TOF object to be stored
// * in OCDB in order to have TOF response correction
// * and actual resolution
// * 
// *
// *
// *

class TGraph;

#include "TObject.h"
#include "AliPID.h"

class AliTOFResponseParams :
public TObject
{

 public:

  AliTOFResponseParams(); // default constructor
  AliTOFResponseParams(Int_t *nPoints); // standard constructor
  virtual ~AliTOFResponseParams(); // default destructor
  AliTOFResponseParams(const AliTOFResponseParams &source); // copy constructor
  AliTOFResponseParams &operator=(const AliTOFResponseParams &source); // operator=

  Int_t GetNPoints(Int_t ipart) const {return ipart < AliPID::kSPECIES ? fNPoints[ipart] : 0;}; // getter
  Double_t GetP(Int_t ipart, Int_t ipoint) const {return ipart < AliPID::kSPECIES && ipoint < fNPoints[ipart] ? fP[ipart][ipoint] : 0.;}; // getter
  Double_t GetTExpCorr(Int_t ipart, Int_t ipoint) const {return ipart < AliPID::kSPECIES && ipoint < fNPoints[ipart] ? fTExpCorr[ipart][ipoint] : 0.;}; // getter

  void SetP(Int_t ipart, Double_t *value) {if (ipart < AliPID::kSPECIES) for (Int_t ipoint = 0; ipoint < fNPoints[ipart]; ipoint++) fP[ipart][ipoint] = value[ipoint];}; // setter
  void SetTExpCorr(Int_t ipart, Double_t *value) {if (ipart < AliPID::kSPECIES) for (Int_t ipoint = 0; ipoint < fNPoints[ipart]; ipoint++) fTExpCorr[ipart][ipoint] = value[ipoint];}; // setter

  Double_t EvalTExpCorr(Int_t ipart, Double_t p); // eval corr
  TGraph *DrawGraph(Int_t ipart, Option_t* option = ""); // draw

 private:

  static const Int_t fgkMaxPoints = 20; // max number of points
  Int_t fNPoints[AliPID::kSPECIES]; // number of points
  Double_t fP[AliPID::kSPECIES][fgkMaxPoints]; // average momentum (GeV/c)
  Double_t fTExpCorr[AliPID::kSPECIES][fgkMaxPoints]; // expected time correction (ps)

  ClassDef(AliTOFResponseParams, 1);
};

#endif /* ALITOFRESPONSEPARAMS_H */
