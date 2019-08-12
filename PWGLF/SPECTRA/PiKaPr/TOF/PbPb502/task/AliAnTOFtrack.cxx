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

////////////////////////////////////////////////////////////////////////////
///                                                                       //
///                                                                       //
/// Class container for identified charged hadron spectra: TOF            //
///                                                                       //
///                                                                       //
/// Authors:                                                              //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                        //
////////////////////////////////////////////////////////////////////////////

#define LOG_NO_INFO
#define LOG_NO_DEBUG
#include "AliAnTOFtrack.h"
#include "TError.h"
#include <iostream>

//________________________________________________________________________
AliAnTOFtrack::AliAnTOFtrack()
    : fTrkMask(0)
    , fTPCPIDMask(0)
    , fTrkCutMask(0)
    , fDCAXY(-999)
    , fDCAZ(-999)
    , fLength(-999)
    , fLengthRatio(-999)
    , fTOFTime(-999)
    , fTOFMismatchTime(-999)
    , fT0TrkTime(-999)
    , fTOFchan(-999)
    , fEta(-999)
    , fPhi(-999)
    , fPt(-999)
    , fPTPC(-999)
    , fNTOFClusters(-999)
    , fTPCSignal(-999)
// fT0TrkSigma(-999)
{ // standard constructur which should be used

#ifndef LOG_NO_INFO
  ::Info("AliAnTOFtrack::AliAnTOFtrack", "**** CONSTRUCTOR CALLED ****");
#endif

  for (Int_t i = 0; i < kExpSpecies; i++) {
    fTOFExpTime[i] = -999;
    fTOFExpSigma[i] = -999;
  }

#ifndef LOG_NO_INFO
  ::Info("AliAnTOFtrack::AliAnTOFtrack", "**** END OF CONSTRUCTOR ****");
#endif
}

//________________________________________________________________________
AliAnTOFtrack::~AliAnTOFtrack()
{ //Destructor
#ifndef LOG_NO_INFO
  ::Info("AliAnTOFtrack::~AliAnTOFtrack", "**** DESTRUCTOR CALLED ****");
#endif

#ifndef LOG_NO_INFO
  ::Info("AliAnTOFtrack::~AliAnTOFtrack", "**** END OF DESTRUCTOR ****");
#endif
}

///////////////////
///Cut functions///
///////////////////

//________________________________________________________________________
Bool_t AliAnTOFtrack::PassStdCut() const
{
  for (Int_t i = 0; i < nCuts; i++)
    if (!GetMaskBit(fTrkCutMask, CutStdIndexInMask[i]))
      return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::PassCut(const Int_t cut) const
{
  if (cut < -1 || cut >= nCutVars)
    ::Fatal("AliAnTOFtrack::PassCut", "requested cut is out of bound");
  //
  //Always apply standard cuts except if requiring one different cut
  if (cut == -1)
    return PassStdCut();
  ;

  Int_t cutindex = -1;
  Int_t sum = 0;
  for (Int_t i = 0; i < nCuts; i++) { //Identify if the cut Type is of the Loose kind
    if (cut == sum) {
      cutindex = i;
      break;
    }
    sum += CutIndex[i] - 1; //One cut is not necessary
  }

  //Finally check the cuts
  for (Int_t i = 0; i < nCuts; i++) {
    if (i == cutindex)
      continue;
    else if (!GetMaskBit(fTrkCutMask, CutStdIndexInMask[i]))
      return kFALSE;
  }
  if (cutindex == -1)
    if (!GetMaskBit(fTrkCutMask, cut))
      return kFALSE; //Cut type must be of the Tight kind

  return kTRUE;
}

/////////////////////////
//////PID functions//////
/////////////////////////

//________________________________________________________________________
Float_t AliAnTOFtrack::GetDeltaT(const UInt_t id) const
{
  if (id > kExpSpecies)
    ::Fatal(" AliAnTOFtrack::GetDeltaT", "Index required is out of bound");
  return fTOFTime - fTOFExpTime[id] - fT0TrkTime;
}

//________________________________________________________________________
Float_t AliAnTOFtrack::GetDeltaSigma(const UInt_t id, const UInt_t hypo) const
{
  if (id > kExpSpecies || hypo > kExpSpecies)
    ::Fatal(" AliAnTOFtrack::GetDeltaSigma", "Index required is out of bound");
  return GetDeltaT(hypo) / fTOFExpSigma[id];
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsTPCElectron() const
{
  if (GetMaskBit(fTPCPIDMask, kIsTPCElectron))
    return kTRUE; //1.5 sigma cut for Electrons in TPC
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsTPCPiKP(const UInt_t i) const
{
  if (i >= 3)
    ::Fatal("AliAnTOFtrack::IsTPCPiKP", "Wrong index required");
  if (GetMaskBit(fTPCPIDMask, kIsTPCPion + i))
    return kTRUE; //5 sigma cut for Pi/K/P in TPC
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsTPCPiKP() const
{
  for (Int_t i = 0; i < 3; i++)
    if (IsTPCPiKP(i))
      return kTRUE; //5 sigma cut for Pi/K/P in TPC
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::ConsistentTPCTOF() const
{
  if (!IsTPCPiKP())
    return kFALSE;
  for (Int_t i = 0; i < 3; i++) {
    if (GetMaskBit(fTPCPIDMask, kIsTPCPion + i) && TMath::Abs(GetDeltaSigma(kpi + i, kpi + i)) < 5)
      return kTRUE;
  }
  return kFALSE;
}

///////////////////////
////Track functions////
///////////////////////

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsNegative() const
{
  if (GetMaskBit(fTrkMask, kNegTrk))
    return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Double_t AliAnTOFtrack::GetGammaBeta() const
{
  Double_t beta = GetBeta();
  if (beta == 1.)
    return -999.;
  return beta / (TMath::Sqrt(1. - beta * beta));
}

//________________________________________________________________________
Double_t AliAnTOFtrack::GetTOFMass() const
{
  Double_t gb = GetGammaBeta();
  return gb <= 0 ? gb : GetMomentum() / gb;
}

///////////////////////////
///DCA Utility functions///
///////////////////////////

//________________________________________________________________________
void AliAnTOFtrack::ComputeDCABin(const Double_t dca, const Bool_t xy)
{
  if (xy)
    fDCAXY = dca;
  else
    fDCAZ = dca;
};

//________________________________________________________________________
void AliAnTOFtrack::ComputeDCABin(const Double_t dcaxy, const Double_t dcaz)
{
  ComputeDCABin(dcaxy, (Bool_t)kTRUE);
  ComputeDCABin(dcaz, (Bool_t)kFALSE);
};

//////////////////////////
///T0 Utility functions///
//////////////////////////

//________________________________________________________________________
Double_t AliAnTOFtrack::GetT0Resolution(const Double_t TOFsigma) const
{
  return TMath::Sqrt(TMath::Power(fTOFExpSigma[0], 2) - TMath::Power(TOFsigma, 2) - TMath::Power(15. / GetMomentum(), 2));
}

////////////////////////
///T0 Usage functions///
////////////////////////

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0TOF(const Bool_t exclusive) const
{
  if (GetMaskBit(fTrkMask, kT0_0)) {
    if (exclusive && !GetMaskBit(fTrkMask, kT0_1) && !GetMaskBit(fTrkMask, kT0_2))
      return kTRUE;
    else if (exclusive)
      return kFALSE;
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0A(const Bool_t exclusive) const
{
  if (GetMaskBit(fTrkMask, kT0_1)) {
    if (exclusive && !GetMaskBit(fTrkMask, kT0_0) && !GetMaskBit(fTrkMask, kT0_2))
      return kTRUE;
    else if (exclusive)
      return kFALSE;
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0C(const Bool_t exclusive) const
{
  if (GetMaskBit(fTrkMask, kT0_2)) {
    if (exclusive && !GetMaskBit(fTrkMask, kT0_0) && !GetMaskBit(fTrkMask, kT0_1))
      return kTRUE;
    else if (exclusive)
      return kFALSE;
    return kTRUE;
  }
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0Fill() const
{
  if (!GetMaskBit(fTrkMask, kT0_0) && !GetMaskBit(fTrkMask, kT0_1) && !GetMaskBit(fTrkMask, kT0_2))
    return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0TOF_T0A() const
{
  if (IsT0TOF() && IsT0A())
    return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsOrT0TOF_T0A() const
{
  if (IsT0TOF() || IsT0A())
    return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0TOF_T0C() const
{
  if (IsT0TOF() && IsT0C())
    return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsOrT0TOF_T0C() const
{
  if (IsT0TOF() || IsT0C())
    return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0A_T0C() const
{
  if (IsT0A() && IsT0C())
    return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsOrT0A_T0C() const
{
  if (IsT0A() || IsT0C())
    return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsT0TOF_T0A_T0C() const
{
  if (IsT0TOF() && IsT0A() && IsT0C())
    return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnTOFtrack::IsOrT0TOF_T0A_T0C() const
{
  if (IsT0TOF() || IsT0A() || IsT0C())
    return kTRUE;
  return kFALSE;
}

//________________________________________________________________________
void AliAnTOFtrack::Print() const
{
#define PrintVar(Var) \
  cout << #Var << " " << Var << endl;
  PrintVar(fTrkMask);
  PrintVar(fTPCPIDMask);
  PrintVar(fTrkCutMask);
  PrintVar(fDCAXY);
  PrintVar(fDCAZ);
  PrintVar(fLength);
  PrintVar(fLengthRatio);
  PrintVar(fTOFTime);
  PrintVar(fTOFMismatchTime);
  for (Int_t i = 0; i < kExpSpecies; i++) {
    PrintVar(fTOFExpTime[i]);
    PrintVar(fTOFExpSigma[i]);
  }
  PrintVar(fT0TrkTime);
  PrintVar(fTOFchan);
  PrintVar(fEta);
  PrintVar(fPhi);
  PrintVar(fPt);
  PrintVar(fPTPC);
  PrintVar(fNTOFClusters);
  PrintVar(fTPCSignal);
#undef PrintVar
};
