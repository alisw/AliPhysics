/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
//                Dielectron PID                                  //
//                                                                       //
//                                                                       //
/*
Detailed description


*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TF1.h>
#include <TGraph.h>

#include <AliVTrack.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliESDtrack.h>
#include <AliESDpid.h>
#include <AliAODpidUtil.h>
#include <AliAODPid.h>

#include "AliDielectronVarManager.h"

#include "AliDielectronPID.h"

ClassImp(AliDielectronPID)

TGraph *AliDielectronPID::fgFitCorr=0x0;
Double_t AliDielectronPID::fgCorr=0.0;

AliDielectronPID::AliDielectronPID() :
  AliAnalysisCuts(),
  fNcuts(0),
  fESDpid(0x0),
  fAODpidUtil(0x0)
{
  //
  // Default Constructor
  //
  for (Int_t icut=0; icut<kNmaxPID; ++icut){
    fDetType[icut]=kTPC;
    fPartType[icut]=AliPID::kPion;
    fNsigmaLow[icut]=0;
    fNsigmaUp[icut]=0;
    fPmin[icut]=0;
    fPmax[icut]=0;
    fExclude[icut]=kFALSE;
    fFunUpperCut[icut]=0x0;
    fFunLowerCut[icut]=0x0;
    fRequirePIDbit[icut]=0;
  }
}

//______________________________________________
AliDielectronPID::AliDielectronPID(const char* name, const char* title) :
  AliAnalysisCuts(name, title),
  fNcuts(0),
  fESDpid(0x0),
  fAODpidUtil(0x0)
{
  //
  // Named Constructor
  //
  for (Int_t icut=0; icut<kNmaxPID; ++icut){
    fDetType[icut]=kTPC;
    fPartType[icut]=AliPID::kPion;
    fNsigmaLow[icut]=0;
    fNsigmaUp[icut]=0;
    fPmin[icut]=0;
    fPmax[icut]=0;
    fExclude[icut]=kFALSE;
    fFunUpperCut[icut]=0x0;
    fFunLowerCut[icut]=0x0;
    fRequirePIDbit[icut]=0;
  }
}

//______________________________________________
AliDielectronPID::~AliDielectronPID()
{
  //
  // Default Destructor
  //
}

//______________________________________________
void AliDielectronPID::AddCut(DetType det, AliPID::EParticleType type, Double_t nSigmaLow, Double_t nSigmaUp/*=-99999.*/,
                              Double_t pMin/*=0*/, Double_t pMax/*=0*/, Bool_t exclude/*=kFALSE*/,
                              UInt_t pidBitType/*=AliDielectronPID::kRequire*/)
{
  //
  // Add a pid nsigma cut
  // use response of detector 'det' in the momentum range ('pMin') to ['pMax']
  // use a sigma band betwee 'nSigmaLow' and 'nSigmaUp'
  // if nSigmaUp==-99999. then nSigmaLow will be uesd as a symmetric band +- nSigmaLow
  // specify whether to 'exclude' the given band
  //

  if (fNcuts>=kNmaxPID){
    AliError(Form("only %d pid cut ranges allowed",kNmaxPID));
    return;
  }
  if (TMath::Abs(nSigmaUp+99999.)<1e-20){
    nSigmaUp=TMath::Abs(nSigmaLow);
    nSigmaLow=-1.*nSigmaUp;
  }
  fDetType[fNcuts]=det;
  fPartType[fNcuts]=type;
  fNsigmaLow[fNcuts]=nSigmaLow;
  fNsigmaUp[fNcuts]=nSigmaUp;
  fPmin[fNcuts]=pMin;
  fPmax[fNcuts]=pMax;
  fExclude[fNcuts]=exclude;
  fRequirePIDbit[fNcuts]=pidBitType;
  ++fNcuts;
  
}

//______________________________________________
void AliDielectronPID::AddCut(DetType det, AliPID::EParticleType type, Double_t nSigmaLow, TF1 * const funUp,
                              Double_t pMin/*=0*/, Double_t pMax/*=0*/, Bool_t exclude/*=kFALSE*/,
                              UInt_t pidBitType/*=AliDielectronPID::kRequire*/)
{
  //
  // cut using a TF1 as upper cut
  //
  if (funUp==0x0){
    AliError("A valid function is required for the upper cut. Not adding the cut!");
    return;
  }
  fFunUpperCut[fNcuts]=funUp;
  AddCut(det,type,nSigmaLow,0.,pMin,pMax,exclude,pidBitType);
}

//______________________________________________
void AliDielectronPID::AddCut(DetType det, AliPID::EParticleType type, TF1 * const funLow, Double_t nSigmaUp,
                              Double_t pMin/*=0*/, Double_t pMax/*=0*/, Bool_t exclude/*=kFALSE*/,
                              UInt_t pidBitType/*=AliDielectronPID::kRequire*/)
{
  //
  // cut using a TF1 as lower cut
  //
  if (funLow==0x0){
    AliError("A valid function is required for the lower cut. Not adding the cut!");
    return;
  }
  fFunLowerCut[fNcuts]=funLow;
  AddCut(det,type,0.,nSigmaUp,pMin,pMax,exclude,pidBitType);
}

//______________________________________________
void AliDielectronPID::AddCut(DetType det, AliPID::EParticleType type, TF1 * const funLow, TF1 * const funUp,
                              Double_t pMin/*=0*/, Double_t pMax/*=0*/, Bool_t exclude/*=kFALSE*/,
                              UInt_t pidBitType/*=AliDielectronPID::kRequire*/)
{
  //
  // cut using a TF1 as lower and upper cut
  //
  if ( (funUp==0x0) || (funLow==0x0) ){
    AliError("A valid function is required for upper and lower cut. Not adding the cut!");
    return;
  }
  fFunUpperCut[fNcuts]=funUp;
  fFunLowerCut[fNcuts]=funLow;
  AddCut(det,type,0.,0.,pMin,pMax,exclude,pidBitType);
}

//______________________________________________
Bool_t AliDielectronPID::IsSelected(TObject* track)
{
  //
  // perform PID cuts
  //

  //loop over all cuts
  AliVTrack *part=static_cast<AliVTrack*>(track);
  //TODO: Which momentum to use?
  //      Different momenta for different detectors?
  Double_t mom=part->P();
  
  Bool_t selected=kFALSE;
  fESDpid=AliDielectronVarManager::GetESDpid();
  fAODpidUtil=AliDielectronVarManager::GetAODpidUtil();
  
  for (UChar_t icut=0; icut<fNcuts; ++icut){
    Double_t pMin=fPmin[icut];
    Double_t pMax=fPmax[icut];

    // test momentum range. In case pMin==pMax use all momenta
    if ( (TMath::Abs(pMin-pMax)>1e-20) && (mom<=pMin || mom>pMax) ) continue;


    switch (fDetType[icut]){
    case kITS:
      selected = IsSelectedITS(part,icut);
      break;
    case kTPC:
      selected = IsSelectedTPC(part,icut);
      break;
    case kTRD:
      selected = IsSelectedTRD(part,icut);
      break;
    case kTOF:
      selected = IsSelectedTOF(part,icut);
      break;
    }
    if (!selected) return kFALSE;
  }

  return selected;
}

//______________________________________________
Bool_t AliDielectronPID::IsSelectedITS(AliVTrack * const part, Int_t icut)
{
  //
  // ITS part of the PID check
  // Don't accept the track if there was no pid bit set
  //
  Float_t numberOfSigmas=-1000.;
  
  if (fRequirePIDbit[icut]==AliDielectronPID::kRequire&&!(part->GetStatus()&AliESDtrack::kITSpid)) return kFALSE;
  if (fRequirePIDbit[icut]==AliDielectronPID::kIfAvailable&&!(part->GetStatus()&AliESDtrack::kITSpid)) return kTRUE;

  Double_t mom=part->P();
  
  if (part->IsA()==AliESDtrack::Class()){
    // ESD case in case the PID bit is not set, don't use this track!
    AliESDtrack *track=static_cast<AliESDtrack*>(part);
    numberOfSigmas=fESDpid->NumberOfSigmasITS(track, fPartType[icut]);
  }else if(part->IsA()==AliAODTrack::Class()){
    // AOD case
    // FIXME: Is there a place to check whether the PID is was set in ESD???
    AliAODTrack *track=static_cast<AliAODTrack*>(part);
    numberOfSigmas=fAODpidUtil->NumberOfSigmasITS(track, fPartType[icut]);
  }
  
  // test if we are supposed to use a function for the cut
  if (fFunUpperCut[icut]) fNsigmaUp[icut] =fFunUpperCut[icut]->Eval(mom);
  if (fFunLowerCut[icut]) fNsigmaLow[icut]=fFunLowerCut[icut]->Eval(mom);
  
  Bool_t selected=((numberOfSigmas>=fNsigmaLow[icut])&&(numberOfSigmas<=fNsigmaUp[icut]))^fExclude[icut];
  return selected;
}

//______________________________________________
Bool_t AliDielectronPID::IsSelectedTPC(AliVTrack * const part, Int_t icut)
{
  //
  // TPC part of the PID check
  // Don't accept the track if there was no pid bit set
  //
  Float_t numberOfSigmas=-1000.;
  
  if (fRequirePIDbit[icut]==AliDielectronPID::kRequire&&!(part->GetStatus()&AliESDtrack::kTPCpid)) return kFALSE;
  if (fRequirePIDbit[icut]==AliDielectronPID::kIfAvailable&&!(part->GetStatus()&AliESDtrack::kTPCpid)) return kTRUE;

  Double_t mom=part->P();
  
  if (part->IsA()==AliESDtrack::Class()){
    // ESD case in case the PID bit is not set, don't use this track!
    AliESDtrack *track=static_cast<AliESDtrack*>(part);
    const AliExternalTrackParam *in = track->GetInnerParam();
    if (in) mom = in->GetP();
    numberOfSigmas=fESDpid->NumberOfSigmasTPC(track, fPartType[icut]);
  }else if(part->IsA()==AliAODTrack::Class()){
    // AOD case
    // FIXME: Is there a place to check whether the PID is was set in ESD???
    AliAODTrack *track=static_cast<AliAODTrack*>(part);
    const AliAODPid *pidObj = track->GetDetPid();
    if (pidObj) mom = pidObj->GetTPCmomentum();
    numberOfSigmas=fAODpidUtil->NumberOfSigmasTPC(track, fPartType[icut]);
  }
  if (fPartType[icut]==AliPID::kElectron){
    numberOfSigmas-=fgCorr;
  }
  
  // test if we are supposed to use a function for the cut
  if (fFunUpperCut[icut]) fNsigmaUp[icut] =fFunUpperCut[icut]->Eval(mom);
  if (fFunLowerCut[icut]) fNsigmaLow[icut]=fFunLowerCut[icut]->Eval(mom);
  
  Bool_t selected=((numberOfSigmas>=fNsigmaLow[icut])&&(numberOfSigmas<=fNsigmaUp[icut]))^fExclude[icut];
  return selected;
}

//______________________________________________
Bool_t AliDielectronPID::IsSelectedTRD(AliVTrack * const /*part*/, Int_t /*icut*/)
{
  //   
  // TRD part of the pid check
  //
  return kFALSE;
}

//______________________________________________
Bool_t AliDielectronPID::IsSelectedTOF(AliVTrack * const part, Int_t icut)
{
  //
  // TOF part of the PID check
  // Don't accept the track if there was no pid bit set
  //
  Float_t numberOfSigmas=-1000.;
  
  if (fRequirePIDbit[icut]==AliDielectronPID::kRequire&&!(part->GetStatus()&AliESDtrack::kTOFpid)) return kFALSE;
  if (fRequirePIDbit[icut]==AliDielectronPID::kIfAvailable&&!(part->GetStatus()&AliESDtrack::kTOFpid)) return kTRUE;

  if (part->IsA()==AliESDtrack::Class()){
    // ESD case in case the PID bit is not set, don't use this track!
    AliESDtrack *track=static_cast<AliESDtrack*>(part);
    numberOfSigmas=fESDpid->NumberOfSigmasTOF(track, fPartType[icut], fESDpid->GetTOFResponse().GetTimeZero());
  }else if(part->IsA()==AliAODTrack::Class()){
    // AOD case
    // FIXME: Is there a place to check whether the PID is was set in ESD???
    AliAODTrack *track=static_cast<AliAODTrack*>(part);
    numberOfSigmas=fAODpidUtil->NumberOfSigmasTOF(track, fPartType[icut]);
  }
  
  Bool_t selected=((numberOfSigmas>=fNsigmaLow[icut])&&(numberOfSigmas<=fNsigmaUp[icut]))^fExclude[icut];
  return selected;
}

//______________________________________________
void AliDielectronPID::SetDefaults(Int_t def){
  //
  // initialise default pid strategies
  //

  if (def==0){
    // 2sigma bands TPC:
    // - include e
    // - exclude mu,K,pi,p
    // -complete p range
    AddCut(kTPC,AliPID::kElectron,2);
    AddCut(kTPC,AliPID::kMuon,-2.,2.,0.,0.,kTRUE);
    AddCut(kTPC,AliPID::kPion,-2.,2.,0.,0.,kTRUE);
    AddCut(kTPC,AliPID::kKaon,-2.,2.,0.,0.,kTRUE);
    AddCut(kTPC,AliPID::kProton,-2.,2.,0.,0.,kTRUE);
  } else if (def==1) {
    // 2sigma bands TPC:
    // - include e 0<p<inf
    // - exclude mu,K,pi,p 0<p<2

    AddCut(kTPC,AliPID::kElectron,2.);
    AddCut(kTPC,AliPID::kMuon,-2.,2.,0.,2.,kTRUE);
    AddCut(kTPC,AliPID::kPion,-2.,2.,0.,2.,kTRUE);
    AddCut(kTPC,AliPID::kKaon,-2.,2.,0.,2.,kTRUE);
    AddCut(kTPC,AliPID::kProton,-2.,2.,0.,2.,kTRUE);
    
  } else if (def==2) {
    // include 2sigma e TPC
    // 3sigma bands TOF
    // - exclude K,P
    AddCut(kTPC,AliPID::kElectron,-3.,3.);
    AddCut(kTPC,AliPID::kPion,-3.,3.,0.,0.,kTRUE);
    AddCut(kTPC,AliPID::kProton,-3.,3.,0.,0.,kTRUE);

  } else if (def==3) {
    // include 2sigma e TPC
    // 3sigma bands TOF
    // - exclude K 0<p<1
    // - exclude P 0<p<2
    AddCut(kTPC,AliPID::kElectron,2);
    AddCut(kTOF,AliPID::kKaon,-3.,3.,0.,1.,kTRUE);
    AddCut(kTOF,AliPID::kProton,-6.,6.,0.,1.,kTRUE);
    AddCut(kTOF,AliPID::kProton,-3.,3.,1.,2.,kTRUE);
    
  } else if (def==4) {
    // include 2sigma e TPC
    // 3sigma bands TOF
    // - exclude K 0<p<1
    // - exclude P 0<p<2
    AddCut(kTPC,AliPID::kElectron,2.);
    AddCut(kTOF,AliPID::kKaon,-3.,3.,0.,1.,kTRUE);
    AddCut(kTOF,AliPID::kProton,-6.,6.,0.,1.,kTRUE);
    AddCut(kTOF,AliPID::kProton,-3.,3.,1.,2.,kTRUE);
  } else if (def==5) {
    AddCut(kTPC,AliPID::kElectron,-0.5,3);
    AddCut(kTOF,AliPID::kElectron,-3,3,0,1.5);
  } else if (def==6) {
    // lower cut TPC: parametrisation by HFE
    // upper cut TPC: 3 sigma
    // TOF ele band 3sigma 0<p<1.5GeV
    TF1 *lowerCut=new TF1("lowerCut", "[0] * TMath::Exp([1]*x)", 0, 100);
    lowerCut->SetParameters(-2.7,-0.4357);
    AddCut(kTPC,AliPID::kElectron,lowerCut,3.);
    AddCut(kTOF,AliPID::kElectron,-3,3,0,1.5);
  } else if (def==7) {
    // wide TPC cut
    // TOF ele band 3sigma 0<p<1.5GeV
    AddCut(kTPC,AliPID::kElectron,10.);
    AddCut(kTOF,AliPID::kElectron,-3,3,0,1.5);
  } else if (def==8) {
    // TOF 5 sigma inclusion if TOFpid available
    // this should reduce K,p,Pi to a large extent
    AddCut(kTOF,AliPID::kElectron,-5,5,0,200,kFALSE,AliDielectronPID::kIfAvailable);
  } else if (def==9) {
    // lower cut TPC: parametrisation by HFE
    // upper cut TPC: 3 sigma
    // TOF 5 sigma inclusion if TOFpid available
    // this should reduce K,p,Pi to a large extent
    TF1 *lowerCut=new TF1("lowerCut", "[0] * TMath::Exp([1]*x)", 0, 100);
    lowerCut->SetParameters(-2.65,-0.6757);
    AddCut(kTPC,AliPID::kElectron,lowerCut,4.);
    AddCut(kTOF,AliPID::kElectron,-5,5,0,200,kFALSE,AliDielectronPID::kIfAvailable);
  } else if (def==10) {
    AddCut(kTOF,AliPID::kElectron,-5,5,0,200,kFALSE,AliDielectronPID::kIfAvailable);
    AddCut(kTPC,AliPID::kElectron,3.);
    AddCut(kTPC,AliPID::kPion,-3.,3.,0.,0.,kTRUE);
    AddCut(kTPC,AliPID::kProton,-3.,3.,0.,0.,kTRUE);
    
  }
}

//______________________________________________
void AliDielectronPID::SetCorrVal(Double_t run)
{
  //
  // set correction value for run
  //
  fgCorr=0;
  if (!fgFitCorr) return;
  fgCorr=fgFitCorr->Eval(run);
  if (run<fgFitCorr->GetX()[0]) fgCorr=fgFitCorr->GetY()[0];
  if (run>fgFitCorr->GetX()[fgFitCorr->GetN()-1]) fgCorr=fgFitCorr->GetY()[fgFitCorr->GetN()-1];
//   printf("Corr: %f\n",fgCorr);
}

