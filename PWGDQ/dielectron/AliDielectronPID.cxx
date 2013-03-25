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
#include <AliVCluster.h>
#include <AliLog.h>
#include <AliExternalTrackParam.h>
#include <AliPIDResponse.h>
#include <AliTRDPIDResponse.h>
#include <AliESDtrack.h> //!!!!! Remove once Eta correction is treated in the tender
#include <AliAODTrack.h>
#include <AliAODPid.h>

#include "AliDielectronVarManager.h"
#include "AliDielectronVarCuts.h"

#include "AliDielectronPID.h"

ClassImp(AliDielectronPID)

TGraph *AliDielectronPID::fgFitCorr=0x0;
Double_t AliDielectronPID::fgCorr=0.0;
Double_t AliDielectronPID::fgCorrdEdx=1.0;
TF1 *AliDielectronPID::fgFunEtaCorr=0x0;
TGraph *AliDielectronPID::fgdEdxRunCorr=0x0;

AliDielectronPID::AliDielectronPID() :
  AliAnalysisCuts(),
  fNcuts(0),
  fPIDResponse(0x0)
{
  //
  // Default Constructor
  //
  for (Int_t icut=0; icut<kNmaxPID; ++icut){
    fDetType[icut]=kTPC;
    fPartType[icut]=AliPID::kPion;
    fNsigmaLow[icut]=0;
    fNsigmaUp[icut]=0;
    fmin[icut]=0;
    fmax[icut]=0;
    fExclude[icut]=kFALSE;
    fFunUpperCut[icut]=0x0;
    fFunLowerCut[icut]=0x0;
    fRequirePIDbit[icut]=0;
    fActiveCuts[icut]=-1;
    fSigmaFunLow[icut]=0;
    fSigmaFunUp[icut]=0;
    fFunSigma[icut]=0x0;
    fVarCuts[icut]=0x0;
  }
}

//______________________________________________
AliDielectronPID::AliDielectronPID(const char* name, const char* title) :
  AliAnalysisCuts(name, title),
  fNcuts(0),
  fPIDResponse(0x0)
{
  //
  // Named Constructor
  //
  for (Int_t icut=0; icut<kNmaxPID; ++icut){
    fDetType[icut]=kTPC;
    fPartType[icut]=AliPID::kPion;
    fNsigmaLow[icut]=0;
    fNsigmaUp[icut]=0;
    fmin[icut]=0;
    fmax[icut]=0;
    fExclude[icut]=kFALSE;
    fFunUpperCut[icut]=0x0;
    fFunLowerCut[icut]=0x0;
    fRequirePIDbit[icut]=0;
    fActiveCuts[icut]=0;
    fSigmaFunLow[icut]=0;
    fSigmaFunUp[icut]=0;
    fFunSigma[icut]=0x0;
    fVarCuts[icut]=0x0;
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
                              Double_t min/*=0*/, Double_t max/*=0*/, Bool_t exclude/*=kFALSE*/,
                              UInt_t pidBitType/*=AliDielectronPID::kRequire*/, Int_t var/*=-1*/)
{
  //
  // Add a pid nsigma cut
  // use response of detector 'det' in the range ['min'] to ['max'] for var
  // use a sigma band between 'nSigmaLow' and 'nSigmaUp'
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
  fmin[fNcuts]=min;
  fmax[fNcuts]=max;
  fExclude[fNcuts]=exclude;
  fRequirePIDbit[fNcuts]=pidBitType;
  fActiveCuts[fNcuts]=(var==-1 ? AliDielectronVarManager::kP : var);

  AliDebug(1,Form("Add PID cut %d: sigma [% .1f,% .1f] \t cut [% .1f,% .f] \t var %d->%s \n",
		  fNcuts,nSigmaLow,nSigmaUp,min,max,fActiveCuts[fNcuts],AliDielectronVarManager::GetValueName(fActiveCuts[fNcuts])));
  
  ++fNcuts;

}

//______________________________________________
void AliDielectronPID::AddCut(DetType det, AliPID::EParticleType type, Double_t nSigmaLow, TF1 * const funUp,
                              Double_t min/*=0*/, Double_t max/*=0*/, Bool_t exclude/*=kFALSE*/,
                              UInt_t pidBitType/*=AliDielectronPID::kRequire*/, Int_t var/*=-1*/)
{
  //
  // cut using a TF1 as upper cut
  //
  if (funUp==0x0){
    AliError("A valid function is required for the upper cut. Not adding the cut!");
    return;
  }
  fFunUpperCut[fNcuts]=funUp;
  AddCut(det,type,nSigmaLow,0.,min,max,exclude,pidBitType,var);
}

//______________________________________________
void AliDielectronPID::AddCut(DetType det, AliPID::EParticleType type, TF1 * const funLow, Double_t nSigmaUp,
                              Double_t min/*=0*/, Double_t max/*=0*/, Bool_t exclude/*=kFALSE*/,
                              UInt_t pidBitType/*=AliDielectronPID::kRequire*/, Int_t var/*=-1*/)
{
  //
  // cut using a TF1 as lower cut
  //
  if (funLow==0x0){
    AliError("A valid function is required for the lower cut. Not adding the cut!");
    return;
  }
  fFunLowerCut[fNcuts]=funLow;
  AddCut(det,type,0.,nSigmaUp,min,max,exclude,pidBitType,var);
}

//______________________________________________
void AliDielectronPID::AddCut(DetType det, AliPID::EParticleType type, TF1 * const funLow, TF1 * const funUp,
                              Double_t min/*=0*/, Double_t max/*=0*/, Bool_t exclude/*=kFALSE*/,
                              UInt_t pidBitType/*=AliDielectronPID::kRequire*/, Int_t var/*=-1*/)
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
  AddCut(det,type,0.,0.,min,max,exclude,pidBitType,var);
}

//______________________________________________
void AliDielectronPID::AddCut(DetType det, AliPID::EParticleType type, Double_t nSigmaLow, Double_t nSigmaUp,
                              Double_t min, Double_t max, Bool_t exclude,
                              UInt_t pidBitType, TF1 * const funSigma)
{
  //
  // cut using a TF1 as lower cut
  //
  if (funSigma==0x0){
    AliError("A valid function is required for the sigma cut. Not adding the cut!");
    return;
  }
  fFunSigma[fNcuts]=funSigma;
  fSigmaFunLow[fNcuts]=min;
  fSigmaFunUp[fNcuts]=max;

  AddCut(det,type,nSigmaLow,nSigmaUp,0,0,exclude,pidBitType,-1);
}

//______________________________________________
void AliDielectronPID::AddCut(DetType det, AliPID::EParticleType type, Double_t nSigmaLow, Double_t nSigmaUp,
                              AliDielectronVarCuts *var, Bool_t exclude/*=kFALSE*/,
                              UInt_t pidBitType/*=AliDielectronPID::kRequire*/)
{
  //
  // Add a pid nsigma cut
  // use response of detector 'det' in the ranges for variables defined in var
  // use a sigma band between 'nSigmaLow' and 'nSigmaUp'
  // if nSigmaUp==-99999. then nSigmaLow will be uesd as a symmetric band +- nSigmaLow
  // specify whether to 'exclude' the given band
  //
  if(!var) return;
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
  fExclude[fNcuts]=exclude;
  fRequirePIDbit[fNcuts]=pidBitType;
  fVarCuts[fNcuts]=var;

  AliDebug(1,Form("Add PID cut %d: sigma [% .1f,% .1f] \n",
		  fNcuts,nSigmaLow,nSigmaUp));
  
  ++fNcuts;

}

//______________________________________________
Bool_t AliDielectronPID::IsSelected(TObject* track)
{
  //
  // perform PID cuts
  //

  //loop over all cuts
  AliVTrack *part=static_cast<AliVTrack*>(track);
  AliESDtrack *esdTrack=0x0;
  AliAODTrack *aodTrack=0x0;
  Double_t origdEdx=-1;
  
  // apply ETa correction, remove once this is in the tender
  if( (part->IsA() == AliESDtrack::Class()) ){
    esdTrack=static_cast<AliESDtrack*>(part);
    origdEdx=esdTrack->GetTPCsignal();
    esdTrack->SetTPCsignal(origdEdx/GetEtaCorr(esdTrack)/fgCorrdEdx,esdTrack->GetTPCsignalSigma(),esdTrack->GetTPCsignalN());
  } else if ( (part->IsA() == AliAODTrack::Class()) ){
    aodTrack=static_cast<AliAODTrack*>(track);
    AliAODPid *pid=const_cast<AliAODPid*>(aodTrack->GetDetPid());
    if (pid){
      origdEdx=pid->GetTPCsignal();
      pid->SetTPCsignal(origdEdx/GetEtaCorr(aodTrack)/fgCorrdEdx);
    }
  }
  
  //Fill values
  Double_t values[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(track,values);

  Bool_t selected=kFALSE;
  fPIDResponse=AliDielectronVarManager::GetPIDResponse();
  for (UChar_t icut=0; icut<fNcuts; ++icut){
    Double_t min=fmin[icut];
    Double_t max=fmax[icut];
    Double_t val=values[fActiveCuts[icut]];

    // test var range. In case min==max do not cut
    if ( fVarCuts[icut] ) {
      if ( !fVarCuts[icut]->IsSelected(part) ) continue;
    }
    else if ( ( (TMath::Abs(min-max)>1e-20) && (val<min || val>=max) ) ) {
	continue;
    }

    // check if fFunSigma is set, then check if 'part' is in sigma range of the function
    if(fFunSigma[icut]){
        val= fPIDResponse->NumberOfSigmasTPC(part, fPartType[icut]);
        if (fPartType[icut]==AliPID::kElectron){
            val-=fgCorr;
        }
        min= fFunSigma[icut]->Eval(part->GetTPCmomentum())+fSigmaFunLow[icut];
        max= fFunSigma[icut]->Eval(part->GetTPCmomentum())+fSigmaFunUp[icut];
        if(val<min || val>=max) continue;
    }

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
    case kTRDeleEff:
      selected = IsSelectedTRDeleEff(part,icut);
      break;
    case kTRDeleEff2D:
      selected = IsSelectedTRDeleEff(part,icut,AliTRDPIDResponse::kLQ2D);
      break;
    case kTOF:
      selected = IsSelectedTOF(part,icut);
      break;
    case kEMCAL:
      selected = IsSelectedEMCAL(part,icut);
      break;
    }
    if (!selected) {
      if (esdTrack) esdTrack->SetTPCsignal(origdEdx,esdTrack->GetTPCsignalSigma(),esdTrack->GetTPCsignalN());
      else if (aodTrack){
        AliAODPid *pid=const_cast<AliAODPid*>(aodTrack->GetDetPid());
        if (pid) pid->SetTPCsignal(origdEdx);
      }

      return kFALSE;
    }
  }

  if (esdTrack) esdTrack->SetTPCsignal(origdEdx,esdTrack->GetTPCsignalSigma(),esdTrack->GetTPCsignalN());
  else if (aodTrack){
    AliAODPid *pid=const_cast<AliAODPid*>(aodTrack->GetDetPid());
    if (pid) pid->SetTPCsignal(origdEdx);
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
  AliPIDResponse::EDetPidStatus pidStatus = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,part);
  if (fRequirePIDbit[icut]==AliDielectronPID::kRequire&&(pidStatus!=AliPIDResponse::kDetPidOk)) return kFALSE;
  if (fRequirePIDbit[icut]==AliDielectronPID::kIfAvailable&&(pidStatus!=AliPIDResponse::kDetPidOk)) return kTRUE;

  Double_t mom=part->P();
  
  Float_t numberOfSigmas=fPIDResponse->NumberOfSigmasITS(part, fPartType[icut]);
  
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
  AliPIDResponse::EDetPidStatus pidStatus = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,part);
  if (fRequirePIDbit[icut]==AliDielectronPID::kRequire&&(pidStatus!=AliPIDResponse::kDetPidOk)) return kFALSE;
  if (fRequirePIDbit[icut]==AliDielectronPID::kIfAvailable&&(pidStatus!=AliPIDResponse::kDetPidOk)) return kTRUE;

  Double_t mom=part->GetTPCmomentum();
  
  Float_t numberOfSigmas=fPIDResponse->NumberOfSigmasTPC(part, fPartType[icut]);

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
Bool_t AliDielectronPID::IsSelectedTRD(AliVTrack * const part, Int_t icut)
{
  //   
  // TRD part of the pid check
  // the TRD checks on the probabilities.
  //

  AliPIDResponse::EDetPidStatus pidStatus = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTRD,part);
  if (fRequirePIDbit[icut]==AliDielectronPID::kRequire&&(pidStatus!=AliPIDResponse::kDetPidOk)) return kFALSE;
  if (fRequirePIDbit[icut]==AliDielectronPID::kIfAvailable&&(pidStatus!=AliPIDResponse::kDetPidOk)) return kTRUE;

  if (fRequirePIDbit[icut]==AliDielectronPID::kIfAvailable && (part->GetTRDntrackletsPID()<4)) return kTRUE;

  Double_t p[AliPID::kSPECIES]={0.};
  fPIDResponse->ComputeTRDProbability(part,AliPID::kSPECIES,p);
  Float_t particleProb=p[fPartType[icut]];

  Bool_t selected=((particleProb>=fNsigmaLow[icut])&&(particleProb<=fNsigmaUp[icut]))^fExclude[icut];
  return selected;
}

//______________________________________________
Bool_t AliDielectronPID::IsSelectedTRDeleEff(AliVTrack * const part, Int_t icut, AliTRDPIDResponse::ETRDPIDMethod PIDmethod)
{
  //
  // TRD part of the pid check using electron efficiency requirement
  // in this case the upper limit as well as the particle specie is ignored 
  //   and the lower limit regarded as the requested electron efficiency
  //

  AliPIDResponse::EDetPidStatus pidStatus = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTRD,part);
  if (fRequirePIDbit[icut]==AliDielectronPID::kRequire&&(pidStatus!=AliPIDResponse::kDetPidOk)) return kFALSE;
  if (fRequirePIDbit[icut]==AliDielectronPID::kIfAvailable&&(pidStatus!=AliPIDResponse::kDetPidOk)) return kTRUE;

  Double_t centrality = -1.;
  if(part->IsA() == AliESDtrack::Class())
    centrality=(const_cast<AliESDEvent*>( (static_cast<const AliESDtrack*>(part))->GetESDEvent()) )->GetCentrality()->GetCentralityPercentile("V0M");
  if(part->IsA() == AliAODTrack::Class())
    centrality=(const_cast<AliAODEvent*>( (static_cast<const AliAODTrack*>(part))->GetAODEvent()) )->GetCentrality()->GetCentralityPercentile("V0M");

  Bool_t selected=fPIDResponse->IdentifiedAsElectronTRD(part,fNsigmaLow[icut], centrality, PIDmethod)^fExclude[icut];
  return selected;
}

//______________________________________________
Bool_t AliDielectronPID::IsSelectedTOF(AliVTrack * const part, Int_t icut)
{
  //
  // TOF part of the PID check
  // Don't accept the track if there was no pid bit set
  //
  AliPIDResponse::EDetPidStatus pidStatus = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,part);
  if (fRequirePIDbit[icut]==AliDielectronPID::kRequire&&(pidStatus!=AliPIDResponse::kDetPidOk)) return kFALSE;
  if (fRequirePIDbit[icut]==AliDielectronPID::kIfAvailable&&(pidStatus!=AliPIDResponse::kDetPidOk)) return kTRUE;

  Float_t numberOfSigmas=fPIDResponse->NumberOfSigmasTOF(part, fPartType[icut]);
  
  Bool_t selected=((numberOfSigmas>=fNsigmaLow[icut])&&(numberOfSigmas<=fNsigmaUp[icut]))^fExclude[icut];
  return selected;
}

//______________________________________________
Bool_t AliDielectronPID::IsSelectedEMCAL(AliVTrack * const part, Int_t icut)
{
  //
  // emcal pid selecttion
  //

  //TODO: correct way to check for emcal pid?
  Float_t numberOfSigmas=fPIDResponse->NumberOfSigmasEMCAL(part, fPartType[icut]);

  Bool_t hasPID=numberOfSigmas>-998.;

  if (fRequirePIDbit[icut]==AliDielectronPID::kRequire&&!hasPID) return kFALSE;
  if (fRequirePIDbit[icut]==AliDielectronPID::kIfAvailable&&!hasPID) return kTRUE;


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

  } else if (def==3 || def==4) { // def3 and def4 are the same??
    // include 2sigma e TPC
    // 3sigma bands TOF
    // - exclude K 0<p<1
    // - exclude P 0<p<2
    AddCut(kTPC,AliPID::kElectron,2);
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
    
  } else if (def==11) {
    // lower cut TPC: parametrisation by HFE
    // only use from period d on !!
    // upper cut TPC: 3 sigma
    // TOF ele band 3sigma 0<p<2.0GeV
    TF1 *lowerCut=new TF1("lowerCut", "[0] * TMath::Exp([1]*x)+[2]", 0, 100);
    lowerCut->SetParameters(-3.718,-0.4,0.27);
    AddCut(kTPC,AliPID::kElectron,lowerCut,3.);
    AddCut(kTOF,AliPID::kElectron,-3,3,0,2.);
  } else if (def==12) {
    // lower cut TPC: parametrisation by HFE
    // only use from period d on !!
    // upper cut TPC: 3 sigma
    // TOF 5 sigma inclusion if TOFpid available
    // this should reduce K,p,Pi to a large extent
    TF1 *lowerCut=new TF1("lowerCut", "[0] * TMath::Exp([1]*x)+[2]", 0, 100);
    lowerCut->SetParameters(-3.718,-0.4,0.27);
    AddCut(kTPC,AliPID::kElectron,lowerCut,4.);
    AddCut(kTOF,AliPID::kElectron,-5,5,0,200,kFALSE,AliDielectronPID::kIfAvailable);
  }
  else if (def==13) {
    // TPC electron inclusion
    // TOF electron inclusion if available
    AddCut(kTOF,AliPID::kElectron,-4.,4.,0,200,kFALSE,AliDielectronPID::kIfAvailable);
    AddCut(kTPC,AliPID::kElectron,-3.5,3.5);
  }
  else if (def==14) {
    // TRD 1D 90% elec eff, 4-6 tracklets
    // TPC electron inclusion
    // TOF electron inclusion if available
    AddCut(AliDielectronPID::kTRDeleEff,AliPID::kElectron,.9,1.,3.5,6.,kFALSE,
	   AliDielectronPID::kIfAvailable,AliDielectronVarManager::kTRDpidQuality);
    AddCut(kTOF,AliPID::kElectron,-4.,4.,0,200,kFALSE,AliDielectronPID::kIfAvailable);
    AddCut(kTPC,AliPID::kElectron,-3.5,3.5);
  }
  else if (def==15) {
    // TRD 1D 90% elec eff, 4-6 tracklets, chi2 < 2
    // TPC electron inclusion
    // TOF electron inclusion if available
    AddCut(AliDielectronPID::kTRDeleEff,AliPID::kElectron,.9,1.,3.5,6.,kFALSE,
	   AliDielectronPID::kIfAvailable,AliDielectronVarManager::kTRDpidQuality);
    AddCut(AliDielectronPID::kTRDeleEff,AliPID::kElectron,.9,1.,0.,2.,kFALSE,
	   AliDielectronPID::kIfAvailable,AliDielectronVarManager::kTRDchi2);
    AddCut(kTOF,AliPID::kElectron,-4.,4.,0,200,kFALSE,AliDielectronPID::kIfAvailable);
    AddCut(kTPC,AliPID::kElectron,-3.5,3.5);
  }

}

//______________________________________________
void AliDielectronPID::SetCorrVal(Double_t run)
{
  //
  // set correction value for run
  //
  fgCorr=0.;
  fgCorrdEdx=1.;
  
  if (fgFitCorr){
    fgCorr=fgFitCorr->Eval(run);
    if (run<fgFitCorr->GetX()[0]) fgCorr=fgFitCorr->GetY()[0];
    if (run>fgFitCorr->GetX()[fgFitCorr->GetN()-1]) fgCorr=fgFitCorr->GetY()[fgFitCorr->GetN()-1];
  }

  if (fgdEdxRunCorr){
    fgCorrdEdx=fgdEdxRunCorr->Eval(run);
    if (run<fgdEdxRunCorr->GetX()[0]) fgCorrdEdx=fgdEdxRunCorr->GetY()[0];
    if (run>fgdEdxRunCorr->GetX()[fgFitCorr->GetN()-1]) fgCorrdEdx=fgdEdxRunCorr->GetY()[fgdEdxRunCorr->GetN()-1];
  }
}

//______________________________________________
Double_t AliDielectronPID::GetEtaCorr(const AliVTrack *track)
{
  //
  // return eta correction
  //
  if (!fgFunEtaCorr) return 1;
  return fgFunEtaCorr->Eval(track->Eta());
}
