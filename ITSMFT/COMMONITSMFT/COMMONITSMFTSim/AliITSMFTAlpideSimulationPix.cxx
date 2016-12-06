
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TF1.h>
#include <TF2.h>
#include <TRandom.h>
#include <TLorentzVector.h>

#include "AliITSMFTAlpideSimulationPix.h"
#include "AliITSMFTSimuClusterShaper.h"
#include "AliITSMFTSDigit.h"
#include "AliITSMFTChip.h"
#include "AliITSMFTGeomTGeo.h"
#include "AliITSMFTSimuParam.h"
#include "AliITSMFTDigitPix.h"
#include "AliLog.h"

ClassImp(AliITSMFTAlpideSimulationPix)

//______________________________________________________________________
AliITSMFTAlpideSimulationPix::AliITSMFTAlpideSimulationPix() {}


//______________________________________________________________________
AliITSMFTAlpideSimulationPix::AliITSMFTAlpideSimulationPix(AliITSMFTSimuParam* sim, AliITSMFTSensMap* map)
:AliITSMFTSimulation(sim, map) {
  Init();
}


//______________________________________________________________________
AliITSMFTAlpideSimulationPix::~AliITSMFTAlpideSimulationPix() {}


//______________________________________________________________________
void AliITSMFTAlpideSimulationPix::Init() {}


//______________________________________________________________________
void AliITSMFTAlpideSimulationPix::SDigitiseChip(TClonesArray *sdarray) {
  if (fChip->GetNHits()) GenerateCluster();
  if (!fSensMap->GetEntries()) return;
  WriteSDigits(sdarray);
  ClearMap();
}


//______________________________________________________________________
void AliITSMFTAlpideSimulationPix::FrompListToDigits(TObjArray *detDigits) {
  int nsd = fSensMap->GetEntries();
  if (!nsd) return; // nothing to digitize

  UInt_t row,col;
  Int_t iCycle, modId = fChip->GetIndex();
  static AliITSMFTDigitPix dig;

  for (int i = 0; i < nsd; ++i) {
      AliITSMFTSDigit* sd = (AliITSMFTSDigit*) fSensMap->At(i); // ordered in index
      if (fSensMap->IsDisabled(sd)) continue;

      fSensMap->GetMapIndex(sd->GetUniqueID(),col,row,iCycle);
      dig.SetCoord1(col);
      dig.SetCoord2(row);

      //I.B.
      Double_t sig = sd->GetSumSignal();
      dig.SetROCycle(iCycle);
      dig.SetSignal((Int_t)sig);
      dig.SetSignalPix((Int_t)sig);
      int ntr = sd->GetNTracks();
      for (int j=0;j<ntr;j++) {
          dig.SetTrack(j,sd->GetTrack(j));
          dig.SetHit(j,sd->GetHit(j));
      }
      for (int j=ntr;j<AliITSMFTDigitPix::GetNTracks();j++) {
          dig.SetTrack(j,-3);
          dig.SetHit(j,-1);
      }
      //I.B.
      
      Int_t branch = AliITSMFTAux::kChipTypePix;
      TClonesArray &ldigits = *((TClonesArray*) detDigits->At(branch));
      int nd = ldigits.GetEntriesFast();
      new (ldigits[nd]) AliITSMFTDigitPix(dig);
  }
}


//______________________________________________________________________
void AliITSMFTAlpideSimulationPix::WriteSDigits(TClonesArray *sdarray) {
  //  This function adds each S-Digit to pList
  int nsd = fSensMap->GetEntries();

  for (int i = 0; i < nsd; ++i) {
    AliITSMFTSDigit* sd = (AliITSMFTSDigit*)fSensMap->At(i); // ordered in index
    if (fSensMap->IsDisabled(sd)) continue;
    new ((*sdarray)[sdarray->GetEntriesFast()]) AliITSMFTSDigit(*sd);
  }
  return;
}


//______________________________________________________________________
Bool_t AliITSMFTAlpideSimulationPix::AddSDigitsToChip(TSeqCollection *pItemArr, Int_t mask) {
  //    pItemArr  Array of AliITSpListItems (SDigits).
  //    mask    Track number off set value
  Int_t nItems = pItemArr->GetEntries();

  for( Int_t i=0; i<nItems; i++ ) {
    AliITSMFTSDigit * pItem = (AliITSMFTSDigit *)(pItemArr->At( i ));
    if(pItem->GetChip() != int(fChip->GetIndex()) ) AliFatal(Form("SDigits chip %d != current chip %d: exit", pItem->GetChip(),fChip->GetIndex()));

    AliITSMFTSDigit* oldItem = (AliITSMFTSDigit*)fSensMap->GetItem(pItem);
    if (!oldItem) oldItem = (AliITSMFTSDigit*)fSensMap->RegisterItem( new(fSensMap->GetFree()) AliITSMFTSDigit(*pItem) );
  }
  return true;
}


//______________________________________________________________________
void AliITSMFTAlpideSimulationPix::FinishSDigitiseChip(TObjArray *detDigits) {
  //  This function calls SDigitsToDigits which creates Digits from SDigits
  FrompListToDigits(detDigits);
  ClearMap();
}


//______________________________________________________________________
void AliITSMFTAlpideSimulationPix::DigitiseChip(TObjArray *detDigits) {
  GenerateCluster();
  FinishSDigitiseChip(detDigits);
}


//______________________________________________________________________
void AliITSMFTAlpideSimulationPix::SetResponseParam(AliITSMFTParamList* resp) {
  fResponseParam = resp;
}


//______________________________________________________________________
Double_t AliITSMFTAlpideSimulationPix::ACSFromBetaGamma(Double_t x, Double_t theta) const {
  TF1 *acs = new TF1("acs", "[0]*((1+TMath::Power(x, 2))/TMath::Power(x, 2))*(0.5*TMath::Log([1]*TMath::Power(x, 2)) - (TMath::Power(x, 2)/(1+TMath::Power(x, 2))) - [2]*TMath::Log(x))", 0, 10000);
  Double_t par1 = fSimuParam->GetACSFromBGPar1();
  Double_t par2 = fSimuParam->GetACSFromBGPar2();
  Double_t par3 = fSimuParam->GetACSFromBGPar3();
  acs->SetParameter(0, par1);
  acs->SetParameter(1, par2);
  acs->SetParameter(2, par3);
  Double_t val = acs->Eval(x)/fabs(cos(theta));
  delete acs;
  return val;
}


//______________________________________________________________________
Int_t AliITSMFTAlpideSimulationPix::GetPixelPositionResponse(Int_t idPadX, Int_t idPadZ, Float_t locx, Float_t locz, Double_t acs) const {
  Float_t centerX, centerZ;
  fSeg->DetToLocal(idPadX, idPadZ, centerX, centerZ);

  Double_t Dx = locx-centerX;
  Double_t Dy = locz-centerZ;
  Double_t sigma = 0.001; // = 10 um
  Double_t offc  = acs; // WARNING: this is just temporary! (a function for this is ready but need further testing)

  TF2 *respf = new TF2("respf", "([1]-1)*(1-TMath::Gaus(x,0,[0])*TMath::Gaus(y,0,[0]))+1",
                       -fSeg->GetPitchX()/2, fSeg->GetPitchX()/2, -fSeg->GetPitchZ()/2, fSeg->GetPitchZ()/2);
  respf->SetParameter(0, sigma);
  respf->SetParameter(1, offc);
  Int_t cs = (Int_t) round(respf->Eval(Dx, Dy));
  delete respf;
  return cs;
}


//______________________________________________________________________
Int_t AliITSMFTAlpideSimulationPix::CSSampleFromLandau(Double_t mpv, Double_t w) const {
  TF1 *landauDistr = new TF1("landauDistr","TMath::Landau(x,[0],[1])", 0, 20);
  landauDistr->SetParameter(0, mpv);
  landauDistr->SetParameter(1, w);

  // Generation according to the Landau distribution defined above
  Double_t fmax = landauDistr->GetMaximum();
  Double_t x1 = gRandom->Uniform(0, 20);
  Double_t y1 = gRandom->Uniform(0, fmax);
  while (y1 > landauDistr->Eval(x1)) {
    x1 = gRandom->Uniform(0, 20);
    y1 = gRandom->Uniform(0, fmax);
  }
  Int_t cs = (Int_t) round(x1);
  delete landauDistr;
  return cs;
}


//______________________________________________________________________
Double_t AliITSMFTAlpideSimulationPix::ComputeIncidenceAngle(TLorentzVector dir) const {
  Double_t glob[3], loc[3];
  glob[0] = dir.Px()/dir.P();
  glob[1] = dir.Py()/dir.P();
  glob[2] = dir.Pz()/dir.P();

  AliITSMFTGeomTGeo* geomTG = fChip->GetITSMFTGeomTGeo();
  geomTG->GlobalToLocalVect(fChip->GetIndex(), glob, loc);

  TVector3 pdirection(loc[0], loc[1], loc[2]);
  TVector3 normal(0., -1., 0.);

  return pdirection.Angle(normal);
}


//______________________________________________________________________
void AliITSMFTAlpideSimulationPix::GenerateCluster() {
  TObjArray *hits = fChip->GetHits();
  Int_t nhits = hits->GetEntriesFast();
  if (nhits <= 0) return;

  Float_t px, py, pz, x,  y,  z;
  Double_t etot, beta, gamma, acs, theta;
  Double_t tof, x0, x1, y0, y1, z0, z1, el, de;
  UInt_t cs, nrows, ncols;
  Int_t pid, ix, iz, nx, nz, cx, cz, idtrack;

  for (Int_t h=0; h < nhits; ++h) {
    if (!fChip->LineSegmentL(h, x0, x1, y0, y1, z0, z1, de, tof, idtrack)) continue;

    AliITSMFTHit *hit = fChip->GetHit(h);
    hit->GetMomentum(px, py, pz);

    // local coordinates
    /* I.B.
    x = (x0+x1)/2.0;
    y = (y0+y1)/2.0;
    z = (z0+z1)/2.0;
    */
    x = x0 + 0.5*x1;
    y = y0 + 0.5*y1;
    z = z0 + 0.5*z1;

    etot = hit->GetTotalEnergy();
    pid = hit->GetPID();

    TLorentzVector momen;
    momen.SetPxPyPzE(px, py, pz, etot);
    beta = momen.Beta();
    if (beta > 0.99999) return; //I.B. FIXME
    gamma = momen.Gamma();
    if (beta*gamma < 0.1) return; //I.B. FIXME
    theta = ComputeIncidenceAngle(momen);

    // Get the pixel ID
    if(!fSeg->LocalToDet(x,z,ix,iz)) return;

    acs = ACSFromBetaGamma(beta*gamma, theta);
    cs = GetPixelPositionResponse(ix, iz, x, z, acs);
    UInt_t *cshape = new UInt_t[cs];
    AliITSMFTSimuClusterShaper *csManager = new AliITSMFTSimuClusterShaper(cs);
    csManager->FillClusterRandomly();
    cshape = csManager->GetShape();
    nrows = csManager->GetNRows();
    ncols = csManager->GetNCols();
    cx = gRandom->Integer(ncols);
    cz = gRandom->Integer(nrows);

    AliDebug(10,Form("_/_/_/_/_/_/_/_/_/_/_/_/_/_/"));
    AliDebug(10,Form("_/_/_/ pALPIDE debug  _/_/_/"));
    AliDebug(10,Form("_/_/_/_/_/_/_/_/_/_/_/_/_/_/"));
    AliDebug(10,Form(" Beta*Gamma: %f", beta*gamma));
    AliDebug(10,Form("        ACS: %f", acs));
    AliDebug(10,Form("         CS: %d", cs));
    AliDebug(10,Form("      Shape: %s", csManager->ShapeSting(cs, cshape).c_str()));
    AliDebug(10,Form("     Center: %d, %d", cx, cz));
    AliDebug(10,Form("_/_/_/_/_/_/_/_/_/_/_/_/_/_/"));

    for (Int_t ipix = 0; ipix < cs; ++ipix) {
      UInt_t r = (Int_t) cshape[ipix] / nrows;
      UInt_t c = cshape[ipix] % nrows;
      nx = ix - cx + c;
      nz = iz - cz + r;
      CreateDigi(nz, nx, idtrack, h);
    }
    
    delete[] cshape;
    delete csManager;
  }
}


//______________________________________________________________________
void AliITSMFTAlpideSimulationPix::CreateDigi(UInt_t col, UInt_t row, Int_t track, Int_t hit) {
  Int_t index = fSensMap->GetIndex(col, row, 0);
  Int_t chip  = fChip->GetIndex();

  fSensMap->RegisterItem(new (fSensMap->GetFree()) AliITSMFTSDigit(track, hit, chip, index, 0.1, 0));
}
