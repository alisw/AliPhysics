/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
// \class for cuts on reconstructed c-deuteron -> dKpi
// \First implementation by copying the AliRDHFCutsXicTopKpi class
// \inheriting from this class
// \author Author: J. Norman (jaime.norman@cern.ch)
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliPIDResponse.h>

#include "AliRDHFCutsCdeuterontodKpi.h"
#include "AliRDHFCutsXictopKpi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliRDHFCuts.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliKFParticle.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliRDHFCutsCdeuterontodKpi);
/// \endcond

//--------------------------------------------------------------------------
AliRDHFCutsCdeuterontodKpi::AliRDHFCutsCdeuterontodKpi(const char* name) : 
AliRDHFCutsXictopKpi(name)
{
  //
  // Default Constructor
  //
//  Int_t nvars=13;
//  SetNVars(nvars);
//  TString varNames[13]={"inv. mass [GeV]",
//			"pTK [GeV/c]",
//			"pTP [GeV/c]",
//			"d0K [cm]   lower limit!",
//			"d0Pi [cm]  lower limit!",
//			"dist12 (cm)",
//			"sigmavert (cm)",
//			"dist prim-sec (cm)",
//			"pM=Max{pT1,pT2,pT3} (GeV/c)",
//			"cosThetaPoint",
//			"Sum d0^2 (cm^2)",
//			"dca cut (cm)",
//			"cut on pTpion [GeV/c]"};
//  Bool_t isUpperCut[13]={kTRUE,
//			 kFALSE,
//			 kFALSE,
//			 kFALSE,
//			 kFALSE,
//			 kFALSE,
//			 kTRUE,
//			 kFALSE,
//			 kFALSE,
//			 kFALSE,
//			 kFALSE,
//			 kTRUE,
//			 kFALSE
//			 };
//  SetVarNames(nvars,varNames,isUpperCut);
//  Bool_t forOpt[13]={kFALSE,
//		     kTRUE,
//		     kTRUE,
//		     kFALSE,
//		     kFALSE,
//		     kFALSE,
//		     kFALSE,
//		     kTRUE,
//		     kFALSE,
//		     kFALSE,
//		     kFALSE,
//		     kFALSE,
//		     kTRUE};
//  SetVarsForOpt(4,forOpt);
//  Float_t limits[2]={0,999999999.};
//  SetPtBins(2,limits);
//  for (Int_t ispecies=0;ispecies<AliPID::kSPECIES;++ispecies)
//      fPIDThreshold[ispecies]=0.;
}
//--------------------------------------------------------------------------
AliRDHFCutsCdeuterontodKpi::AliRDHFCutsCdeuterontodKpi(const AliRDHFCutsCdeuterontodKpi &source) :
  AliRDHFCutsXictopKpi(source)
{
  //
  // Copy constructor
  //
//  if (source.fPidObjprot) fPidObjprot = new AliAODPidHF(*(source.fPidObjprot));
//  else fPidObjprot = new AliAODPidHF();
//  if (source.fPidObjpion) fPidObjpion = new AliAODPidHF(*(source.fPidObjpion));
//  else fPidObjpion = new AliAODPidHF();
//  memcpy(fPIDThreshold,source.fPIDThreshold,AliPID::kSPECIES*sizeof(Double_t));
}
//--------------------------------------------------------------------------
AliRDHFCutsCdeuterontodKpi &AliRDHFCutsCdeuterontodKpi::operator=(const AliRDHFCutsCdeuterontodKpi &source)
{
  //
  // assignment operator
  //
  if(this != &source) {
    AliRDHFCuts::operator=(source);
  }
    
  return *this;
}
//---------------------------------------------------------------------------
AliRDHFCutsCdeuterontodKpi::~AliRDHFCutsCdeuterontodKpi() {
 //
 //  // Default Destructor
 //   

}

//---------------------------------------------------------------------------
void AliRDHFCutsCdeuterontodKpi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters, AliAODEvent *aod) {
  // 
  // Fills in vars the values of the variables 
  //

  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsCdeuterontodKpi::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoDecayHF3Prong *dd = (AliAODRecoDecayHF3Prong*)d;

  Int_t iter=-1;
  if(fVarsForOpt[0]){
    iter++;
    vars[iter]=dd->InvMassCdeuterondKpi();
  }
  if(fVarsForOpt[1]){
    iter++;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==321) {
	vars[iter]=dd->PtProng(iprong);
      }
    }
  }
  if(fVarsForOpt[2]){
    iter++;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==1000010020) {
	vars[iter]=dd->PtProng(iprong);
      }
    }
  }
  if(fVarsForOpt[3]){
    iter++;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==1000010020) {
	vars[iter]=dd->Getd0Prong(iprong);
      }
    }
  }
  if(fVarsForOpt[4]){
    iter++;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==211) {
	vars[iter]=dd->Getd0Prong(iprong);
      }
    }
  }
  if(fVarsForOpt[5]){
    iter++;
    vars[iter]=dd->GetDist12toPrim();
  }
  if(fVarsForOpt[6]){
    iter++;
    vars[iter]=dd->GetSigmaVert(aod);
  }
  if(fVarsForOpt[7]){
    iter++;
    vars[iter] = dd->DecayLength();
  }
  if(fVarsForOpt[8]){
    iter++;
    Float_t ptmax=0;
    for(Int_t i=0;i<3;i++){
      if(dd->PtProng(i)>ptmax)ptmax=dd->PtProng(i);
    }
    vars[iter]=ptmax;
  }
  if(fVarsForOpt[9]){
    iter++;
    vars[iter]=dd->CosPointingAngle();
  }
  if(fVarsForOpt[10]){
    iter++;
    vars[iter]=dd->Getd0Prong(0)*dd->Getd0Prong(0)+dd->Getd0Prong(1)*dd->Getd0Prong(1)+dd->Getd0Prong(2)*dd->Getd0Prong(2);
  }
  if(fVarsForOpt[11]){
    iter++;
    vars[iter]=dd->GetDCA();
  }
  if(fVarsForOpt[12]){
    iter++;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==211) {
	vars[iter]=dd->PtProng(iprong);
      }
    }
  }

  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsCdeuterontodKpi::IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent *aod) {
  //
  // Apply selection
  //

  //Printf("---- WE ARE IN THE C-DEUTERON CLASS - WELL DONE ---");
  if(!fCutsRD){
    AliError("Cut matrice not inizialized. Exit...\n");
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF3Prong* d=(AliAODRecoDecayHF3Prong*)obj;

  if(!d){
    AliError("AliAODRecoDecayHF3Prong null \n");
    return 0;
  }

  //
  // NB: pdg code of Lc... shall we change it ???
  //
  if(fKeepSignalMC) if(IsSignalMC(d,aod,4122)) return 3;

  Int_t returnvalue=3;
  Int_t returnvaluePID=3;

  if(d->Pt()<fMinPtCand) return 0;
  if(d->Pt()>fMaxPtCand) return 0;
  //Printf("1: Pt cut passed");

  if(fUseTrackSelectionWithFilterBits && d->HasBadDaughters()) return 0;
  //Printf("2: track selection filter bit and bad daughter passed");


  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d,aod)) return 0;
  }
  //Printf("3: selection on daughter tracks passed");


  // PID selection
  if(selectionLevel==AliRDHFCuts::kAll ||
     selectionLevel==AliRDHFCuts::kCandidate|| 
     selectionLevel==AliRDHFCuts::kPID) {
    switch (fPIDStrategy) {
      Printf("check PID strategy = %i",fPIDStrategy);
    case kNSigma:
      returnvaluePID = IsSelectedPID(d);
      break;
    case kNSigmaMin:
      returnvaluePID = IsSelectedPID(d);
      break;
    case kNSigmaPbPb:
      returnvaluePID = IsSelectedNSigmaPbPb(d);
      break;
    case kCombined:
      returnvaluePID = IsSelectedCombinedPID(d);
      break;
    case kCombinedSoft:
      returnvaluePID = IsSelectedCombinedPIDSoft(d);
      break;
    case kNSigmaStrong:
      returnvaluePID = IsSelectedPIDStrong(d);
      break;
    case kCombinedpPb:
      returnvaluePID = IsSelectedCombinedPIDpPb(d);
      break;
    case kCombinedpPb2:
      returnvaluePID = IsSelectedCombinedPIDpPb2(d);
      break;
    case kCombinedProb:
      returnvaluePID = IsSelectedCombinedPIDProb(d);
      break;
    }
    fIsSelectedPID=returnvaluePID;
    //Printf("is selected PID = %i",fIsSelectedPID);
  }
  //  if(fUsePID || selectionLevel==AliRDHFCuts::kPID) returnvaluePID = IsSelectedCombinedPID(d);   // to test!!
  if(returnvaluePID==0) return 0;
  //Printf("4: passed PID selection");




  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {

    Double_t pt=d->Pt();
    //Printf("\tcandidate selection");
    //Printf("\tpt = %f",pt);
    
    Int_t ptbin=PtBin(pt);
    
    Double_t mCdeuterondKpi=0.,mCdeuteronpiKd=0.;
    Int_t okCdeuterondKpi=1,okCdeuteronpiKd=1;

    Double_t mLcPDG = 3.226;

    mCdeuterondKpi=d->InvMassCdeuterondKpi();
    mCdeuteronpiKd=d->InvMassCdeuteronpiKd();
    //Printf("\tmass c deuteron dKpi = %f",mCdeuterondKpi);
    //Printf("\tmass c deuteron piKd = %f",mCdeuteronpiKd);

    // Comparison with Lc mass from PDG.
    // From the cutobject, only the upper limit is set.
    // The lower limit is hardocoded, in order to discard candidates with 
    // mass lower than the Lc.
    if(TMath::Abs(mCdeuterondKpi-mLcPDG) > fCutsRD[GetGlobalIndex(0,ptbin)]) okCdeuterondKpi = 0;
    if(TMath::Abs(mCdeuteronpiKd-mLcPDG) > fCutsRD[GetGlobalIndex(0,ptbin)]) okCdeuteronpiKd = 0;
    if(!okCdeuterondKpi && !okCdeuteronpiKd) return 0;
    //Printf("5: passed invariant mass selection");

  switch (fCutsStrategy) {

    case kStandard:
    if(TMath::Abs(d->PtProng(1)) < fCutsRD[GetGlobalIndex(1,ptbin)] || TMath::Abs(d->Getd0Prong(1))<fCutsRD[GetGlobalIndex(3,ptbin)]) return 0;//Kaon
    if(d->Pt()>=3. && d->PProng(1)<0.55) return 0;
    if(fUseSpecialCut) {
      if(TMath::Abs(d->PtProng(0)) < TMath::Abs(d->PtProng(2)) )okCdeuterondKpi=0;
      if(TMath::Abs(d->PtProng(2)) < TMath::Abs(d->PtProng(0)) )okCdeuteronpiKd=0;
    }
    if((TMath::Abs(d->PtProng(0)) < fCutsRD[GetGlobalIndex(2,ptbin)]) || (TMath::Abs(d->PtProng(2)) < fCutsRD[GetGlobalIndex(12,ptbin)])) okCdeuterondKpi=0;
    if((TMath::Abs(d->PtProng(2)) < fCutsRD[GetGlobalIndex(2,ptbin)]) || (TMath::Abs(d->PtProng(0)) < fCutsRD[GetGlobalIndex(12,ptbin)]))okCdeuteronpiKd=0;
    if(!okCdeuterondKpi && !okCdeuteronpiKd) return 0;
    //2track cuts
      // postponed
    //if(d->GetDist12toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]|| d->GetDist23toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]) return 0;
    //if(d->GetDist12toPrim()>0.5) return 0;
    //if(d->GetDist23toPrim()>0.5) return 0;
    if(fUseImpParProdCorrCut){
      if(d->Getd0Prong(0)*d->Getd0Prong(1)<0. && d->Getd0Prong(2)*d->Getd0Prong(1)<0.) return 0;
    }
    //sec vert
    if(d->DecayLength()<fCutsRD[GetGlobalIndex(7,ptbin)]) return 0;
    if(d->DecayLength()>0.5) return 0;

  //  Double_t sumd0s=d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
  //  if(sumd0s<fCutsRD[GetGlobalIndex(10,ptbin)]) return 0;
    if((d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2))<fCutsRD[GetGlobalIndex(10,ptbin)]) return 0;
    
    if(TMath::Abs(d->PtProng(0))<fCutsRD[GetGlobalIndex(8,ptbin)] && TMath::Abs(d->PtProng(1))<fCutsRD[GetGlobalIndex(8,ptbin)] && TMath::Abs(d->PtProng(2))<fCutsRD[GetGlobalIndex(8,ptbin)]) return 0;
    if(d->CosPointingAngle()< fCutsRD[GetGlobalIndex(9,ptbin)]) return 0;
    if(d->GetSigmaVert(aod)>fCutsRD[GetGlobalIndex(6,ptbin)]) return 0;
    
    //DCA
    for(Int_t i=0;i<3;i++) if(d->GetDCA(i)>fCutsRD[GetGlobalIndex(11,ptbin)]) return 0;

    // dist12 and dist23
    if(d->GetDist12toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]|| d->GetDist23toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]) return 0;
    if(d->GetDist12toPrim()>0.5) return 0;
    if(d->GetDist23toPrim()>0.5) return 0;

    break;


   case kKF:
    Int_t pdgs[3]={0,321,0};
    Bool_t constraint=kFALSE;
    if(fCutsRD[GetGlobalIndex(1,ptbin)]>0.) constraint=kTRUE;
    Double_t field=aod->GetMagneticField();
    if (returnvaluePID==1 || returnvaluePID==3){

      pdgs[0]=2122;pdgs[2]=211;
      AliKFParticle *lc1=ReconstructKF(d,pdgs,field,constraint);
      if(!lc1){
	okCdeuterondKpi=0;
      }else{
	if(lc1->GetChi2()/lc1->GetNDF()>fCutsRD[GetGlobalIndex(2,ptbin)]) okCdeuterondKpi=0;
      }
    } else if(returnvaluePID>=2){

      pdgs[0]=211;pdgs[2]=1000010020;
      AliKFParticle *lc2=ReconstructKF(d,pdgs,field,constraint);
      if(!lc2){ 
	okCdeuteronpiKd=0;
      }else{
	if(lc2->GetChi2()/lc2->GetNDF()>fCutsRD[GetGlobalIndex(2,ptbin)])okCdeuteronpiKd=0; 
      }
    }
    break;

   }

    if(okCdeuterondKpi) returnvalue=1; //cuts passed as Xic->pKpi
    if(okCdeuteronpiKd) returnvalue=2; //cuts passed as Xic->piKp
    if(okCdeuterondKpi && okCdeuteronpiKd) returnvalue=3; //cuts passed as both pKpi and piKp
   
  }


  Int_t returnvalueTot=CombinePIDCuts(returnvalue,returnvaluePID);
  //Printf("\tcuts return value = %i",returnvalue);
  //Printf("\tPID  return value = %i",returnvaluePID);
  //Printf("\tfinal cuts return value = %i",returnvalueTot);
  return returnvalueTot;
}



Int_t AliRDHFCutsCdeuterontodKpi::IsSelectedPID(AliAODRecoDecayHF* obj) {

  //Printf(" -- IN THE PID SELECTION FUNCTION - SICK");

    if(!fUsePID || !obj) return 3;
    Int_t okLcpKpi=0,okLcpiKp=0;
    Int_t returnvalue=0;
    Bool_t isMC=fPidHF->GetMC();
    Bool_t ispion0=kTRUE,ispion2=kTRUE;
    Bool_t isdeuteron0=kFALSE,isdeuteron2=kFALSE;
    Bool_t iskaon1=kFALSE;
    if(isMC) {
     fPidObjprot->SetMC(kTRUE);
     fPidObjpion->SetMC(kTRUE);
    }

   if(fPidObjprot->GetPidResponse()==0x0){
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
      AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
      fPidObjprot->SetPidResponse(pidResp);
    }
    if(fPidObjpion->GetPidResponse()==0x0){
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
      AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
      fPidObjpion->SetPidResponse(pidResp);
    }
    if(fPidHF->GetPidResponse()==0x0){
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
      AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
      fPidHF->SetPidResponse(pidResp);
    }

    for(Int_t i=0;i<3;i++){
     AliAODTrack *track=(AliAODTrack*)obj->GetDaughter(i);
     if(!track) return 0;
     if(i==1) {
      // identify kaon
      //if(track->P()<0.55){
      // fPidHF->SetTOF(kFALSE);
      // fPidHF->SetTOFdecide(kFALSE);
      //}
      Int_t isKaon=0;
      isKaon=fPIDStrategy==kNSigmaMin?fPidHF->MatchTPCTOFMin(track,3):fPidHF->MakeRawPid(track,3);
      if(isKaon>=1) iskaon1=kTRUE;
    //  if(track->P()<0.55){
    //   fPidHF->SetTOF(kTRUE);
    //   fPidHF->SetTOFdecide(kTRUE);
    //  }


      if(isKaon>=1) iskaon1=kTRUE;

      if(!iskaon1) return 0;

     }else{
     //pion or deuteron
    //if(track->P()<1.){
    //  fPidObjprot->SetTOF(kFALSE);
    //  fPidObjprot->SetTOFdecide(kFALSE);
    // }

     Int_t isDeuteron=0;
     isDeuteron=fPIDStrategy==kNSigmaMin?fPidObjprot->MatchTPCTOFMin(track,5):fPidObjprot->MakeRawPid(track,5);
     Int_t isPion=0;
     isPion=fPIDStrategy==kNSigmaMin?fPidObjpion->MatchTPCTOFMin(track,2):fPidObjpion->MakeRawPid(track,2);

    // if(track->P()<1.){
    //  fPidObjprot->SetTOF(kTRUE);
    //  fPidObjprot->SetTOFdecide(kTRUE);
    // }


     if(i==0) {
      if(isPion<0) ispion0=kFALSE;
      if(isDeuteron>=1) isdeuteron0=kTRUE;

     }
      if(!ispion0 && !isdeuteron0) return 0;
     if(i==2) {
      if(isPion<0) ispion2=kFALSE;
      if(isDeuteron>=1) isdeuteron2=kTRUE;
     }

    }
   }

    //Printf("isdeuteron0 = %i, isdeuteron2 = %i, ispion0 = %i, iskaon2 = %i", isdeuteron0, isdeuteron2, ispion0, ispion2);
    if(ispion2 && isdeuteron0 && iskaon1) okLcpKpi=1;
    if(ispion0 && isdeuteron2 && iskaon1) okLcpiKp=1;
    if(okLcpKpi) returnvalue=1; //cuts passed as Lc->pKpi
    if(okLcpiKp) returnvalue=2; //cuts passed as Lc->piKp
    if(okLcpKpi && okLcpiKp) returnvalue=3; //cuts passed as both pKpi and piKp

 return returnvalue;
}
