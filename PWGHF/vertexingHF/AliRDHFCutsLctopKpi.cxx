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
//
// Class for cuts on AOD reconstructed Lc->pKpi
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <Riostream.h>
#include <AliAnalysisManager.h>
#include <AliInputEventHandler.h>
#include <AliPIDResponse.h>

#include "AliRDHFCutsLctopKpi.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliRDHFCuts.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliKFParticle.h"
#include "AliESDVertex.h"

using std::cout;
using std::endl;

ClassImp(AliRDHFCutsLctopKpi)

//--------------------------------------------------------------------------
AliRDHFCutsLctopKpi::AliRDHFCutsLctopKpi(const char* name) : 
AliRDHFCuts(name),
fPidObjprot(0),
fPidObjpion(0),
fUseImpParProdCorrCut(kFALSE),
fPIDStrategy(kNSigma),
fCutsStrategy(kStandard),
fUseSpecialCut(kFALSE)
{
  //
  // Default Constructor
  //
  Int_t nvars=13;
  SetNVars(nvars);
  TString varNames[13]={"inv. mass [GeV]",
			"pTK [GeV/c]",
			"pTP [GeV/c]",
			"d0K [cm]   lower limit!",
			"d0Pi [cm]  lower limit!",
			"dist12 (cm)",
			"sigmavert (cm)",
			"dist prim-sec (cm)",
			"pM=Max{pT1,pT2,pT3} (GeV/c)",
			"cosThetaPoint",
			"Sum d0^2 (cm^2)",
			"dca cut (cm)",
			"cut on pTpion [GeV/c]"};
  Bool_t isUpperCut[13]={kTRUE,
			 kFALSE,
			 kFALSE,
			 kFALSE,
			 kFALSE,
			 kFALSE,
			 kTRUE,
			 kFALSE,
			 kFALSE,
			 kFALSE,
			 kFALSE,
			 kTRUE,
			 kFALSE
			 };
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[13]={kFALSE,
		     kTRUE,
		     kTRUE,
		     kFALSE,
		     kFALSE,
		     kFALSE,
		     kFALSE,
		     kTRUE,
		     kFALSE,
		     kFALSE,
		     kFALSE,
		     kFALSE,
		     kTRUE};
  SetVarsForOpt(4,forOpt);
  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);
  for (Int_t ispecies=0;ispecies<AliPID::kSPECIES;++ispecies)
      fPIDThreshold[ispecies]=0.;
}
//--------------------------------------------------------------------------
AliRDHFCutsLctopKpi::AliRDHFCutsLctopKpi(const AliRDHFCutsLctopKpi &source) :
  AliRDHFCuts(source),
  fPidObjprot(0x0),
  fPidObjpion(0x0),
  fUseImpParProdCorrCut(source.fUseImpParProdCorrCut),
  fPIDStrategy(source.fPIDStrategy),
  fCutsStrategy(source.fCutsStrategy),
  fUseSpecialCut(source.fUseSpecialCut)
{
  //
  // Copy constructor
  //
  if (source.fPidObjprot) fPidObjprot = new AliAODPidHF(*(source.fPidObjprot));
  else fPidObjprot = new AliAODPidHF();
  if (source.fPidObjpion) fPidObjpion = new AliAODPidHF(*(source.fPidObjpion));
  else fPidObjpion = new AliAODPidHF();
  memcpy(fPIDThreshold,source.fPIDThreshold,AliPID::kSPECIES*sizeof(Double_t));
}
//--------------------------------------------------------------------------
AliRDHFCutsLctopKpi &AliRDHFCutsLctopKpi::operator=(const AliRDHFCutsLctopKpi &source)
{
  //
  // assignment operator
  //
  if(this != &source) {
    
    AliRDHFCuts::operator=(source);
    delete fPidObjprot;
    fPidObjprot = new AliAODPidHF(*(source.fPidObjprot));
    delete fPidObjpion;
    fPidObjpion = new AliAODPidHF(*(source.fPidObjpion));
    fPIDStrategy=source.fPIDStrategy;
    fCutsStrategy=source.fCutsStrategy;
    memcpy(fPIDThreshold,source.fPIDThreshold,AliPID::kSPECIES*sizeof(Double_t));
  }
    
  return *this;
}
//---------------------------------------------------------------------------
AliRDHFCutsLctopKpi::~AliRDHFCutsLctopKpi() {
 //
 //  // Default Destructor
 //   
 if(fPidObjpion){
  delete fPidObjpion;
  fPidObjpion=0;
 }
 if(fPidObjprot){
  delete fPidObjprot;
  fPidObjprot=0;
 }

}

//---------------------------------------------------------------------------
void AliRDHFCutsLctopKpi::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters, AliAODEvent *aod) {
  // 
  // Fills in vars the values of the variables 
  //

  if(nvars!=fnVarsForOpt) {
    printf("AliRDHFCutsLctopKpi::GetCutsVarsForOpt: wrong number of variables\n");
    return;
  }

  AliAODRecoDecayHF3Prong *dd = (AliAODRecoDecayHF3Prong*)d;

    Int_t iter=-1;
  if(fVarsForOpt[0]){
    iter++;
    vars[iter]=dd->InvMassLcpKpi();
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
      if(TMath::Abs(pdgdaughters[iprong])==2212) {
	vars[iter]=dd->PtProng(iprong);
      }
    }
  }
  if(fVarsForOpt[3]){
    iter++;
    for(Int_t iprong=0;iprong<3;iprong++){
      if(TMath::Abs(pdgdaughters[iprong])==2212) {
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
Int_t AliRDHFCutsLctopKpi::IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent *aod) {
  //
  // Apply selection
  //

  if(!fCutsRD){
    cout<<"Cut matrice not inizialized. Exit..."<<endl;
    return 0;
  }
  //PrintAll();
  AliAODRecoDecayHF3Prong* d=(AliAODRecoDecayHF3Prong*)obj;

  if(!d){
    cout<<"AliAODRecoDecayHF3Prong null"<<endl;
    return 0;
  }


  if(fKeepSignalMC) if(IsSignalMC(d,aod,4122)) return 3;

  Int_t returnvalue=3;
  Int_t returnvaluePID=3;

  if(d->Pt()<fMinPtCand) return 0;
  if(d->Pt()>fMaxPtCand) return 0;

  if(fUseTrackSelectionWithFilterBits && d->HasBadDaughters()) return 0;

  if(selectionLevel==AliRDHFCuts::kAll ||
     selectionLevel==AliRDHFCuts::kCandidate|| 
     selectionLevel==AliRDHFCuts::kPID) {
     switch (fPIDStrategy) {
      case kNSigma:
       returnvaluePID = IsSelectedPID(d);
    
      break;
      case kCombined:
       returnvaluePID = IsSelectedCombinedPID(d);
      break;
      case kCombinedSoft:
       returnvaluePID = IsSelectedCombinedPIDSoft(d);
      break;
      case kNSigmaStrong:
       returnvaluePID = IsSelectedPIDStrong(d);
     }
     fIsSelectedPID=returnvaluePID;
  }
  //  if(fUsePID || selectionLevel==AliRDHFCuts::kPID) returnvaluePID = IsSelectedCombinedPID(d);   // to test!!
  if(returnvaluePID==0) return 0;

  // selection on candidate
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kCandidate) {

    Double_t pt=d->Pt();
    
    Int_t ptbin=PtBin(pt);
    
    Double_t mLcpKpi,mLcpiKp;
    Int_t okLcpKpi=1,okLcpiKp=1;

    Double_t mLcPDG = TDatabasePDG::Instance()->GetParticle(4122)->Mass();

    mLcpKpi=d->InvMassLcpKpi();
    mLcpiKp=d->InvMassLcpiKp();

    if(TMath::Abs(mLcpKpi-mLcPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) okLcpKpi = 0;
    if(TMath::Abs(mLcpiKp-mLcPDG)>fCutsRD[GetGlobalIndex(0,ptbin)]) okLcpiKp = 0;
    if(!okLcpKpi && !okLcpiKp) return 0;

  switch (fCutsStrategy) {

    case kStandard:
    if(TMath::Abs(d->PtProng(1)) < fCutsRD[GetGlobalIndex(1,ptbin)] || TMath::Abs(d->Getd0Prong(1))<fCutsRD[GetGlobalIndex(3,ptbin)]) return 0;//Kaon
    if(d->Pt()>=3. && d->PProng(1)<0.55) return 0;
    if(fUseSpecialCut) {
      if(TMath::Abs(d->PtProng(0)) < TMath::Abs(d->PtProng(2)) )okLcpKpi=0;
      if(TMath::Abs(d->PtProng(2)) < TMath::Abs(d->PtProng(0)) )okLcpiKp=0;
    }
    if((TMath::Abs(d->PtProng(0)) < fCutsRD[GetGlobalIndex(2,ptbin)]) || (TMath::Abs(d->PtProng(2)) < fCutsRD[GetGlobalIndex(12,ptbin)])) okLcpKpi=0;
    if((TMath::Abs(d->PtProng(2)) < fCutsRD[GetGlobalIndex(2,ptbin)]) || (TMath::Abs(d->PtProng(0)) < fCutsRD[GetGlobalIndex(12,ptbin)]))okLcpiKp=0;
    if(!okLcpKpi && !okLcpiKp) return 0;
    //2track cuts
    if(d->GetDist12toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]|| d->GetDist23toPrim()<fCutsRD[GetGlobalIndex(5,ptbin)]) return 0;
    if(d->GetDist12toPrim()>0.5) return 0;
    if(d->GetDist23toPrim()>0.5) return 0;
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
	okLcpKpi=0;
      }else{
	if(lc1->GetChi2()/lc1->GetNDF()>fCutsRD[GetGlobalIndex(2,ptbin)]) okLcpKpi=0;
      }
    } else if(returnvaluePID>=2){

      pdgs[0]=211;pdgs[2]=2212;
      AliKFParticle *lc2=ReconstructKF(d,pdgs,field,constraint);
      if(!lc2){ 
	okLcpiKp=0;
      }else{
	if(lc2->GetChi2()/lc2->GetNDF()>fCutsRD[GetGlobalIndex(2,ptbin)])okLcpiKp=0; 
      }
    }
    break;

   }

    if(okLcpKpi) returnvalue=1; //cuts passed as Lc->pKpi
    if(okLcpiKp) returnvalue=2; //cuts passed as Lc->piKp
    if(okLcpKpi && okLcpiKp) returnvalue=3; //cuts passed as both pKpi and piKp
   
  }

  // selection on daughter tracks 
  if(selectionLevel==AliRDHFCuts::kAll || 
     selectionLevel==AliRDHFCuts::kTracks) {
    if(!AreDaughtersSelected(d)) return 0;
  }


  Int_t returnvalueTot=CombinePIDCuts(returnvalue,returnvaluePID);
  return returnvalueTot;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctopKpi::IsSelectedPID(AliAODRecoDecayHF* obj) {


    if(!fUsePID || !obj) return 3;
    Int_t okLcpKpi=0,okLcpiKp=0;
    Int_t returnvalue=0;
    Bool_t isPeriodd=fPidHF->GetOnePad();
    Bool_t isMC=fPidHF->GetMC();
    Bool_t ispion0=kTRUE,ispion2=kTRUE;
    Bool_t isproton0=kFALSE,isproton2=kFALSE;
    Bool_t iskaon1=kFALSE;
    if(isPeriodd) {
     fPidObjprot->SetOnePad(kTRUE);
     fPidObjpion->SetOnePad(kTRUE);
    }
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
     // identify kaon
     if(track->P()<0.55){
      fPidHF->SetTOF(kFALSE);
      fPidHF->SetTOFdecide(kFALSE);
     }
     if(i==1) {
      Int_t isKaon=fPidHF->MakeRawPid(track,3); 
      if(isKaon>=1) iskaon1=kTRUE;
      if(track->P()<0.55){
       fPidHF->SetTOF(kTRUE);
       fPidHF->SetTOFdecide(kTRUE);
      }
      
      if(!iskaon1) return 0;
     
     }else{
     //pion or proton
    if(track->P()<1.){
      fPidObjprot->SetTOF(kFALSE);
      fPidObjprot->SetTOFdecide(kFALSE);
     }
      
     Int_t isProton=fPidObjprot->MakeRawPid(track,4);
     

     Int_t isPion=fPidObjpion->MakeRawPid(track,2);
     
     if(track->P()<1.){
      fPidObjprot->SetTOF(kTRUE);
      fPidObjprot->SetTOFdecide(kTRUE);
     }
     

     if(i==0) {
      if(isPion<0) ispion0=kFALSE;
      if(isProton>=1) isproton0=kTRUE;

     }
      if(!ispion0 && !isproton0) return 0;
     if(i==2) {
      if(isPion<0) ispion2=kFALSE;
      if(isProton>=1) isproton2=kTRUE;
     }

    }
   }

    if(ispion2 && isproton0 && iskaon1) okLcpKpi=1;
    if(ispion0 && isproton2 && iskaon1) okLcpiKp=1;
    if(okLcpKpi) returnvalue=1; //cuts passed as Lc->pKpi
    if(okLcpiKp) returnvalue=2; //cuts passed as Lc->piKp
    if(okLcpKpi && okLcpiKp) returnvalue=3; //cuts passed as both pKpi and piKp

 return returnvalue;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctopKpi::IsSelectedCombinedPID(AliAODRecoDecayHF* obj) {

  //  Printf(" -------- IsSelectedCombinedPID --------------");

    
    if(!fUsePID || !obj) {return 3;}
    Int_t okLcpKpi=0,okLcpiKp=0;
    Int_t returnvalue=0;
    Bool_t isPeriodd=fPidHF->GetOnePad();
    Bool_t isMC=fPidHF->GetMC();

    if(isPeriodd) {
	    fPidObjprot->SetOnePad(kTRUE);
	    fPidObjpion->SetOnePad(kTRUE);
    }
    if(isMC) {
	    fPidObjprot->SetMC(kTRUE);
	    fPidObjpion->SetMC(kTRUE);
    }

    AliVTrack *track0=dynamic_cast<AliVTrack*>(obj->GetDaughter(0));
    AliVTrack *track1=dynamic_cast<AliVTrack*>(obj->GetDaughter(1));
    AliVTrack *track2=dynamic_cast<AliVTrack*>(obj->GetDaughter(2));
    if (!track0 || !track1 || !track2) return 0;
    Double_t prob0[AliPID::kSPECIES];
    Double_t prob1[AliPID::kSPECIES];
    Double_t prob2[AliPID::kSPECIES];
    if(obj->Pt()<3. && track0->P()<1.) fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC);
    fPidHF->GetPidCombined()->ComputeProbabilities(track0,fPidHF->GetPidResponse(),prob0);
    if(obj->Pt()<3. && track0->P()<1.) fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);

   if(obj->Pt()<3. && track1->P()<0.55) fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC);
    fPidHF->GetPidCombined()->ComputeProbabilities(track1,fPidHF->GetPidResponse(),prob1);
   if(obj->Pt()<3. && track1->P()<0.55) fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);

    if(obj->Pt()<3. && track2->P()<1.) fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC);
    fPidHF->GetPidCombined()->ComputeProbabilities(track2,fPidHF->GetPidResponse(),prob2);
   if(obj->Pt()<3. && track2->P()<1.) fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);

    if(fPIDThreshold[AliPID::kPion]>0. && fPIDThreshold[AliPID::kKaon]>0. && fPIDThreshold[AliPID::kProton]>0.){
    okLcpiKp=  (prob0[AliPID::kPion  ]>fPIDThreshold[AliPID::kPion  ])
             &&(prob1[AliPID::kKaon  ]>fPIDThreshold[AliPID::kKaon  ])
             &&(prob2[AliPID::kProton]>fPIDThreshold[AliPID::kProton]);
    okLcpKpi=  (prob0[AliPID::kProton]>fPIDThreshold[AliPID::kProton])
             &&(prob1[AliPID::kKaon  ]>fPIDThreshold[AliPID::kKaon  ])
             &&(prob2[AliPID::kPion  ]>fPIDThreshold[AliPID::kPion  ]);
   }else{ 
		    //pion or proton
		    
		    
    if(TMath::MaxElement(AliPID::kSPECIES,prob1) == prob1[AliPID::kKaon]){
    if(TMath::MaxElement(AliPID::kSPECIES,prob0) == prob0[AliPID::kProton] && TMath::MaxElement(AliPID::kSPECIES,prob2) == prob2[AliPID::kPion]) okLcpKpi = 1;  
    if(TMath::MaxElement(AliPID::kSPECIES,prob2) == prob2[AliPID::kProton] && TMath::MaxElement(AliPID::kSPECIES,prob0) == prob0[AliPID::kPion]) okLcpiKp = 1; 
	    }
	   }
    
    if(okLcpKpi) returnvalue=1; //cuts passed as Lc->pKpi
    if(okLcpiKp) returnvalue=2; //cuts passed as Lc->piKp
    if(okLcpKpi && okLcpiKp) returnvalue=3; //cuts passed as both pKpi and piKp
    
    return returnvalue;
}
//-----------------------
Int_t AliRDHFCutsLctopKpi::CombinePIDCuts(Int_t returnvalue, Int_t returnvaluePID) const {

 Int_t returnvalueTot=0;
 Int_t okLcpKpi=0,okLcpiKp=0;
 if(returnvaluePID==1){
   if(returnvalue==1 || returnvalue==3) okLcpKpi=1;
 }
 if(returnvaluePID==2){
   if(returnvalue>=2) okLcpiKp=1;
 }
 if(returnvaluePID==3 && returnvalue>0){
  if(returnvalue==1 || returnvalue==3) okLcpKpi=1;
  if(returnvalue>=2) okLcpiKp=1;
 } 

 if(okLcpKpi) returnvalueTot=1; //cuts passed as Lc->pKpi
 if(okLcpiKp) returnvalueTot=2; //cuts passed as Lc->piKp
 if(okLcpKpi && okLcpiKp) returnvalueTot=3; //cuts passed as both pKpi and piKp
 return returnvalueTot;
}
//----------------------------------
void AliRDHFCutsLctopKpi::SetStandardCutsPP2010() {

 SetName("LctopKpiProdCuts");
 SetTitle("Production cuts for Lc analysis");

 AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");
 esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
 esdTrackCuts->SetRequireTPCRefit(kTRUE);
 esdTrackCuts->SetMinNClustersTPC(70);
 esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                          AliESDtrackCuts::kAny);
 esdTrackCuts->SetRequireITSRefit(kTRUE);
 esdTrackCuts->SetMinNClustersITS(4);
 esdTrackCuts->SetMinDCAToVertexXY(0.);
 esdTrackCuts->SetEtaRange(-0.8,0.8);
 esdTrackCuts->SetPtRange(0.3,1.e10);
 AddTrackCuts(esdTrackCuts);

 const Int_t nptbins=4;
 const Int_t nvars=13;
 Float_t* ptbins;
 ptbins=new Float_t[nptbins+1];
 
 ptbins[0]=0.;
 ptbins[1]=2.;
 ptbins[2]=3.;
 ptbins[3]=4.;
 ptbins[4]=9999.;

 SetGlobalIndex(nvars,nptbins);
 SetPtBins(nptbins+1,ptbins);

 Float_t** prodcutsval;
 prodcutsval=new Float_t*[nvars];
 for(Int_t iv=0;iv<nvars;iv++){
  prodcutsval[iv]=new Float_t[nptbins];
 }

 for(Int_t ipt=0;ipt<nptbins;ipt++){
  prodcutsval[0][ipt]=0.18;
  prodcutsval[1][ipt]=0.4;
  prodcutsval[2][ipt]=0.5;
  prodcutsval[3][ipt]=0.;
  prodcutsval[4][ipt]=0.;
  prodcutsval[5][ipt]=0.01;
  prodcutsval[6][ipt]=0.06;
  prodcutsval[7][ipt]=0.005;
  prodcutsval[8][ipt]=0.;
  prodcutsval[9][ipt]=0.;
  prodcutsval[10][ipt]=0.;
  prodcutsval[11][ipt]=0.05;
  prodcutsval[12][ipt]=0.4;
 }
 SetCuts(nvars,nptbins,prodcutsval);

 AliAODPidHF* pidObjK=new AliAODPidHF();
 Double_t sigmasK[5]={3.,1.,1.,3.,2.};
 pidObjK->SetSigma(sigmasK);
 pidObjK->SetAsym(kTRUE);
 pidObjK->SetMatch(1);
 pidObjK->SetTPC(kTRUE);
 pidObjK->SetTOF(kTRUE);
 pidObjK->SetITS(kTRUE);
 Double_t plimK[2]={0.5,0.8};
 pidObjK->SetPLimit(plimK,2);
 pidObjK->SetTOFdecide(kTRUE);

 SetPidHF(pidObjK);

 AliAODPidHF* pidObjpi=new AliAODPidHF();
 pidObjpi->SetTPC(kTRUE);
 Double_t sigmaspi[5]={3.,0.,0.,0.,0.};
 pidObjpi->SetSigma(sigmaspi);
 pidObjpi->SetTOFdecide(kTRUE);
 SetPidpion(pidObjpi);

 AliAODPidHF* pidObjp=new AliAODPidHF();
 Double_t sigmasp[5]={3.,1.,1.,3.,2.};
 pidObjp->SetSigma(sigmasp);
 pidObjp->SetAsym(kTRUE);
 pidObjp->SetMatch(1);
 pidObjp->SetTPC(kTRUE);
 pidObjp->SetTOF(kTRUE);
 pidObjp->SetITS(kTRUE);
 Double_t plimp[2]={1.,2.};
 pidObjp->SetPLimit(plimp,2);
 pidObjp->SetTOFdecide(kTRUE);

 SetPidprot(pidObjp);

 SetUsePID(kTRUE);

 PrintAll();

 for(Int_t iiv=0;iiv<nvars;iiv++){
  delete [] prodcutsval[iiv];
 }
 delete [] prodcutsval;
 prodcutsval=NULL;
 delete [] ptbins;
 ptbins=NULL;

 delete pidObjp;
 pidObjp=NULL;

 return;
}
//------------------
void AliRDHFCutsLctopKpi::SetStandardCutsPbPb2010() {

 SetName("LctopKpiProdCuts");
 SetTitle("Production cuts for Lc analysis");

 AliESDtrackCuts *esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts","default");

 esdTrackCuts->SetRequireTPCRefit(kTRUE);
 esdTrackCuts->SetMinNClustersTPC(70);
 esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                          AliESDtrackCuts::kAny);
 esdTrackCuts->SetRequireITSRefit(kTRUE);
 esdTrackCuts->SetMinNClustersITS(4);
 esdTrackCuts->SetMinDCAToVertexXYPtDep("0.0100*TMath::Max(0.,(1-TMath::Floor(TMath::Abs(pt)/2.)))");
 esdTrackCuts->SetEtaRange(-0.8,0.8);
 esdTrackCuts->SetMaxDCAToVertexXY(1.);
 esdTrackCuts->SetMaxDCAToVertexZ(1.);
 esdTrackCuts->SetPtRange(0.8,1.e10);
 AddTrackCuts(esdTrackCuts);

 const Int_t nptbins=4;
 const Int_t nvars=13;
 Float_t* ptbins;
 ptbins=new Float_t[nptbins+1];
 
 ptbins[0]=0.;
 ptbins[1]=2.;
 ptbins[2]=3.;
 ptbins[3]=4.;
 ptbins[4]=9999.;

 SetGlobalIndex(nvars,nptbins);
 SetPtBins(nptbins+1,ptbins);

 Float_t** prodcutsval;
 prodcutsval=new Float_t*[nvars];
 for(Int_t iv=0;iv<nvars;iv++){
  prodcutsval[iv]=new Float_t[nptbins];
 }

 for(Int_t ipt=0;ipt<nptbins;ipt++){
  prodcutsval[0][ipt]=0.13;
  prodcutsval[1][ipt]=0.5;
  prodcutsval[2][ipt]=0.6;
  prodcutsval[3][ipt]=0.;
  prodcutsval[4][ipt]=0.;
  prodcutsval[5][ipt]=0.01;
  prodcutsval[6][ipt]=0.04;
  prodcutsval[7][ipt]=0.006;
  prodcutsval[8][ipt]=0.8;
  prodcutsval[9][ipt]=0.3;
  prodcutsval[10][ipt]=0.;
  prodcutsval[11][ipt]=0.05;
  prodcutsval[12][ipt]=0.4;
 }
 SetCuts(nvars,nptbins,prodcutsval);

 AliAODPidHF* pidObjK=new AliAODPidHF();
 Double_t sigmasK[5]={3.,1.,1.,3.,2.};
 pidObjK->SetSigma(sigmasK);
 pidObjK->SetAsym(kTRUE);
 pidObjK->SetMatch(1);
 pidObjK->SetTPC(kTRUE);
 pidObjK->SetTOF(kTRUE);
 pidObjK->SetITS(kTRUE);
 Double_t plimK[2]={0.5,0.8};
 pidObjK->SetPLimit(plimK,2);

 SetPidHF(pidObjK);

 AliAODPidHF* pidObjpi=new AliAODPidHF();
 pidObjpi->SetTPC(kTRUE);
 Double_t sigmaspi[5]={3.,0.,0.,0.,0.};
 pidObjpi->SetSigma(sigmaspi);
 SetPidpion(pidObjpi);

 AliAODPidHF* pidObjp=new AliAODPidHF();
 Double_t sigmasp[5]={3.,1.,1.,3.,2.};
 pidObjp->SetSigma(sigmasp);
 pidObjp->SetAsym(kTRUE);
 pidObjp->SetMatch(1);
 pidObjp->SetTPC(kTRUE);
 pidObjp->SetTOF(kTRUE);
 pidObjp->SetITS(kTRUE);
 Double_t plimp[2]={1.,2.};
 pidObjp->SetPLimit(plimp,2);

 SetPidprot(pidObjp);

 SetUsePID(kTRUE);

 PrintAll();

 for(Int_t iiv=0;iiv<nvars;iiv++){
  delete [] prodcutsval[iiv];
 }
 delete [] prodcutsval;
 prodcutsval=NULL;
 delete [] ptbins;
 ptbins=NULL;

 delete pidObjp;
 pidObjp=NULL;

 return;
}
//------------------
AliKFParticle* AliRDHFCutsLctopKpi::ReconstructKF(AliAODRecoDecayHF3Prong *d,Int_t *pdgs,Double_t field,Bool_t constraint) const{
  // Method to construct the KF particle from the candidate

 const Int_t nprongs=d->GetNProngs();
 if(nprongs<=0) return 0x0;

 Int_t iprongs[nprongs];
 for(Int_t i=0;i<nprongs;i++) iprongs[i]=i;

 Double_t mass[2]={0.,0.};

 AliKFParticle *decay=d->ApplyVertexingKF(iprongs,nprongs,pdgs,constraint,field,mass);
 if(!decay) return 0x0;
 AliESDVertex *vertexESD = new AliESDVertex(decay->Parameters(),
                                            decay->CovarianceMatrix(),
					    decay->GetChi2(),
					    nprongs);
 Double_t pos[3],cov[6],chi2perNDF;
 vertexESD->GetXYZ(pos);
 vertexESD->GetCovMatrix(cov);
 chi2perNDF = vertexESD->GetChi2toNDF();
 delete vertexESD; vertexESD=NULL;
 AliAODVertex *vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);
 d->SetSecondaryVtx(vertexAOD);
 return decay;
}

//------------------
void AliRDHFCutsLctopKpi::SetStandardCutsPbPb2011() {

  // Default 2010 PbPb cut object
  SetStandardCutsPbPb2010();

  //
  // Enable all 2011 PbPb run triggers
  //  
  SetTriggerClass("");
  ResetMaskAndEnableMBTrigger();
  EnableCentralTrigger();
  EnableSemiCentralTrigger();
}
//-----------------

Bool_t AliRDHFCutsLctopKpi::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
  //
  //  // Checking if Lc is in fiducial acceptance region 
  //    //
  //
  if(fMaxRapidityCand>-998.){
    if(TMath::Abs(y) > fMaxRapidityCand) return kFALSE;
    else return kTRUE;
  }

  if(pt > 5.) {
    // applying cut for pt > 5 GeV
   AliDebug(2,Form("pt of Lc = %f (> 5), cutting at |y| < 0.8",pt));
   if (TMath::Abs(y) > 0.8) return kFALSE;
  
  } else {
   // appliying smooth cut for pt < 5 GeV
   Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5;
   Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;
  AliDebug(2,Form("pt of Lc = %f (< 5), cutting  according to the fiducial zone [%f, %f]\n",pt,minFiducialY,maxFiducialY));
  if (y < minFiducialY || y > maxFiducialY) return kFALSE;
 }
  //
  return kTRUE;
}
//--------------------------------------------------------
Int_t AliRDHFCutsLctopKpi::IsSelectedCombinedPIDSoft(AliAODRecoDecayHF* obj) {
 if(!fUsePID || !obj) {return 3;}
 Int_t okLcpKpi=0,okLcpiKp=0;
 Int_t returnvalue=0;
 
 AliVTrack *track0=dynamic_cast<AliVTrack*>(obj->GetDaughter(0));
 AliVTrack *track1=dynamic_cast<AliVTrack*>(obj->GetDaughter(1));
 AliVTrack *track2=dynamic_cast<AliVTrack*>(obj->GetDaughter(2));
 if (!track0 || !track1 || !track2) return 0;
 Double_t prob0[AliPID::kSPECIES];
 Double_t prob1[AliPID::kSPECIES];
 Double_t prob2[AliPID::kSPECIES];

 Bool_t isTOF0=fPidHF->CheckTOFPIDStatus((AliAODTrack*)obj->GetDaughter(0));
 Bool_t isTOF1=fPidHF->CheckTOFPIDStatus((AliAODTrack*)obj->GetDaughter(1));
 Bool_t isTOF2=fPidHF->CheckTOFPIDStatus((AliAODTrack*)obj->GetDaughter(2));

Bool_t isK1=kFALSE;
 if(isTOF1){ //kaon
  if(track1->P()<1.8) {
    fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);
    if(obj->Pt()<3. && track1->P()<0.55) fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC);
    fPidHF->GetPidCombined()->ComputeProbabilities(track1,fPidHF->GetPidResponse(),prob1);
   
  }else{
    AliAODTrack *trackaod1=(AliAODTrack*)(obj->GetDaughter(1));
    if(trackaod1->P()<0.55){
      fPidHF->SetTOF(kFALSE);
      fPidHF->SetTOFdecide(kFALSE);
     }
    Int_t isKaon=fPidHF->MakeRawPid(trackaod1,3);
    if(isKaon>=1) isK1=kTRUE;
    if(trackaod1->P()<0.55){
      fPidHF->SetTOF(kTRUE);
      fPidHF->SetTOFdecide(kTRUE);
     }
  }
 }else{
  if(track1->P()<0.8){
    fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC);
    fPidHF->GetPidCombined()->ComputeProbabilities(track1,fPidHF->GetPidResponse(),prob0);
  }else{
    AliAODTrack *trackaod1=(AliAODTrack*)(obj->GetDaughter(1));
     if(trackaod1->P()<0.55){
      fPidHF->SetTOF(kFALSE);
      fPidHF->SetTOFdecide(kFALSE);
     }
    Int_t isKaon=fPidHF->MakeRawPid(trackaod1,3);
    if(isKaon>=1) isK1=kTRUE;
     if(trackaod1->P()<0.55){
      fPidHF->SetTOF(kTRUE);
      fPidHF->SetTOFdecide(kTRUE);
     }
  }
 }

 Bool_t ispi0=kFALSE;
 Bool_t isp0=kFALSE;
 Bool_t ispi2=kFALSE;
 Bool_t isp2=kFALSE;

 if(isTOF0){ //proton
  if(track0->P()<2.2) {
    fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);
    if(obj->Pt()<3. && track0->P()<1.) fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC);
    fPidHF->GetPidCombined()->ComputeProbabilities(track0,fPidHF->GetPidResponse(),prob0);
  }else{
   AliAODTrack *trackaod0=(AliAODTrack*)(obj->GetDaughter(0));
   if(trackaod0->P()<1.){
      fPidObjprot->SetTOF(kFALSE);
      fPidObjprot->SetTOFdecide(kFALSE);
   }
   Int_t isProton=fPidObjprot->MakeRawPid(trackaod0,4);
   if(isProton>=1) isp0=kTRUE;
   if(trackaod0->P()<1.){
      fPidObjprot->SetTOF(kTRUE);
      fPidObjprot->SetTOFdecide(kTRUE);
   }
  }
 }else{
   if(track0->P()<1.2){
    fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC);
    fPidHF->GetPidCombined()->ComputeProbabilities(track0,fPidHF->GetPidResponse(),prob0);
   }else{
    AliAODTrack *trackaod0=(AliAODTrack*)(obj->GetDaughter(0));
    if(trackaod0->P()<1.){
      fPidObjprot->SetTOF(kFALSE);
      fPidObjprot->SetTOFdecide(kFALSE);
    }
    Int_t isProton=fPidObjprot->MakeRawPid(trackaod0,4);
    if(isProton>=1) isp0=kTRUE;
    if(trackaod0->P()<1.){
      fPidObjprot->SetTOF(kTRUE);
      fPidObjprot->SetTOFdecide(kTRUE);
    }
   }
 }

 if(isTOF2){ //proton
  if(track2->P()<2.2) {
    fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC|AliPIDResponse::kDetTOF);
    if(obj->Pt()<3. && track2->P()<1.) fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC);
    fPidHF->GetPidCombined()->ComputeProbabilities(track2,fPidHF->GetPidResponse(),prob2);
  }else{
    AliAODTrack *trackaod2=(AliAODTrack*)(obj->GetDaughter(2));
    if(trackaod2->P()<1.){
      fPidObjprot->SetTOF(kFALSE);
      fPidObjprot->SetTOFdecide(kFALSE);
    }
   Int_t isProton=fPidObjprot->MakeRawPid(trackaod2,4);
   if(isProton>=1) isp2=kTRUE;
   if(trackaod2->P()<1.){
      fPidObjprot->SetTOF(kTRUE);
      fPidObjprot->SetTOFdecide(kTRUE);
    }
  }
 }else{
   if(track2->P()<1.2){
    fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC);
    fPidHF->GetPidCombined()->ComputeProbabilities(track2,fPidHF->GetPidResponse(),prob2);
   }else{
    AliAODTrack *trackaod2=(AliAODTrack*)(obj->GetDaughter(2));
     if(trackaod2->P()<1.){
      fPidObjprot->SetTOF(kFALSE);
      fPidObjprot->SetTOFdecide(kFALSE);
    }
    Int_t isProton=fPidObjprot->MakeRawPid(trackaod2,4);
    if(isProton>=1) isp2=kTRUE;
    if(trackaod2->P()<1.){
      fPidObjprot->SetTOF(kTRUE);
      fPidObjprot->SetTOFdecide(kTRUE);
    }
   }
 }
  AliAODTrack *trackaod2=(AliAODTrack*)(obj->GetDaughter(2));
  if(fPidObjpion->MakeRawPid(trackaod2,2)>=1)ispi2=kTRUE;
  AliAODTrack *trackaod0=(AliAODTrack*)(obj->GetDaughter(2));
  if(fPidObjpion->MakeRawPid(trackaod0,2)>=1)ispi0=kTRUE;

  if(TMath::MaxElement(AliPID::kSPECIES,prob1) == prob1[AliPID::kKaon]){

    if(TMath::MaxElement(AliPID::kSPECIES,prob0) == prob0[AliPID::kProton] && TMath::MaxElement(AliPID::kSPECIES,prob2) == prob2[AliPID::kPion]) okLcpKpi = 1;
    if(TMath::MaxElement(AliPID::kSPECIES,prob2) == prob2[AliPID::kProton] && TMath::MaxElement(AliPID::kSPECIES,prob0) == prob0[AliPID::kPion]) okLcpiKp = 1;
    }

   if(!isK1 && TMath::MaxElement(AliPID::kSPECIES,prob1) == prob1[AliPID::kKaon]) isK1=kTRUE;
    if(!ispi0 && TMath::MaxElement(AliPID::kSPECIES,prob0) == prob0[AliPID::kPion]) ispi0=kTRUE;
    if(!ispi2 && TMath::MaxElement(AliPID::kSPECIES,prob2) == prob2[AliPID::kPion]) ispi2=kTRUE;
    if(!isp0 && TMath::MaxElement(AliPID::kSPECIES,prob0) == prob0[AliPID::kProton]) isp0=kTRUE;
    if(!isp2 && TMath::MaxElement(AliPID::kSPECIES,prob2) == prob2[AliPID::kProton]) isp2=kTRUE;
    if(isK1 && ispi0 && isp2) okLcpiKp = 1;
    if(isK1 && isp0 && ispi2) okLcpKpi = 1;

   if(okLcpKpi) returnvalue=1; //cuts passed as Lc->pKpi
    if(okLcpiKp) returnvalue=2; //cuts passed as Lc->piKp
    if(okLcpKpi && okLcpiKp) returnvalue=3; //cuts passed as both pKpi and piKp

    return returnvalue;
 

}
//----------------------------------------------------------
Int_t AliRDHFCutsLctopKpi::IsSelectedPIDStrong(AliAODRecoDecayHF* obj) {


    if(!fUsePID || !obj) return 3;
    Int_t okLcpKpi=0,okLcpiKp=0;
    Int_t returnvalue=0;
    Bool_t isPeriodd=fPidHF->GetOnePad();
    Bool_t isMC=fPidHF->GetMC();
    Bool_t ispion0=kTRUE,ispion2=kTRUE;
    Bool_t isproton0=kFALSE,isproton2=kFALSE;
    Bool_t iskaon1=kFALSE;
    if(isPeriodd) {
     fPidObjprot->SetOnePad(kTRUE);
     fPidObjpion->SetOnePad(kTRUE);
    }
    if(isMC) {
     fPidObjprot->SetMC(kTRUE);
     fPidObjpion->SetMC(kTRUE);
    }

    for(Int_t i=0;i<3;i++){
     AliAODTrack *track=(AliAODTrack*)obj->GetDaughter(i);
     if(!track) return 0;
     // identify kaon
     if(i==1) {
      Int_t isKaon=fPidHF->MakeRawPid(track,3);
      if(isKaon>=1) {
       iskaon1=kTRUE;
      if(fPidHF->MakeRawPid(track,2)>=1) iskaon1=kFALSE;
      }
      if(!iskaon1) return 0;
     
     }else{
     //pion or proton
     
     Int_t isProton=fPidObjprot->MakeRawPid(track,4);
     if(isProton>=1){
      if(fPidHF->MakeRawPid(track,2)>=1) isProton=-1;
      if(fPidHF->MakeRawPid(track,3)>=1) isProton=-1;
     }

     Int_t isPion=fPidObjpion->MakeRawPid(track,2);
     if(fPidHF->MakeRawPid(track,3)>=1) isPion=-1;
     if(fPidObjprot->MakeRawPid(track,4)>=1) isPion=-1;


     if(i==0) {
      if(isPion<0) ispion0=kFALSE;
      if(isProton>=1) isproton0=kTRUE;

     }
      if(!ispion0 && !isproton0) return 0;
     if(i==2) {
      if(isPion<0) ispion2=kFALSE;
      if(isProton>=1) isproton2=kTRUE;
     }

    }
   }

    if(ispion2 && isproton0 && iskaon1) okLcpKpi=1;
    if(ispion0 && isproton2 && iskaon1) okLcpiKp=1;
    if(okLcpKpi) returnvalue=1; //cuts passed as Lc->pKpi
    if(okLcpiKp) returnvalue=2; //cuts passed as Lc->piKp
    if(okLcpKpi && okLcpiKp) returnvalue=3; //cuts passed as both pKpi and piKp

 return returnvalue;
}
