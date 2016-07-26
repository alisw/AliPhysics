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
// Class for cuts on AOD reconstructed LC->e+Lambda
//
// Modified by Y.S Watanabe - wyosuke@cns.s.u-tokyo.ac.jp
//
/////////////////////////////////////////////////////////////

#include <Riostream.h>

#include <TDatabasePDG.h>
#include <TMath.h>
#include <TVector3.h>

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliRDHFCutsLctoeleLambdafromAODtracks.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliESDv0.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliRDHFCutsLctoeleLambdafromAODtracks);
/// \endcond

//--------------------------------------------------------------------------
AliRDHFCutsLctoeleLambdafromAODtracks::AliRDHFCutsLctoeleLambdafromAODtracks(const char* name) :
AliRDHFCuts(name),
  fPIDStrategy(kNSigmaCuts),
  fCombinedPIDThreshold(0.),
  fUseLambdaPID(kFALSE),
  fPidObjProton(0),
  fPidObjPion(0),
  fUseOnTheFlyV0(kFALSE),
  fUseV0Topology(0),
  fBzkG(0),
  fProdTrackTPCNclsPIDMin(0),
  fProdTrackTPCNclsRatioMin(0.0),
  fProdUseAODFilterBit(kTRUE),
  fProdAODFilterBit(4),
  fProdRejectTrackWithShared(kFALSE),
  fProdV0KinkRejection(kTRUE),
  fProdV0MassTolLambda(0.01),
  fProdV0MassTolLambdaRough(0.01),
  fProdV0PtMin(0.5),
  fProdV0PtMax(9999.),
  fProdV0CosPointingAngleToPrimVtxMin(0.99),
  fProdV0DcaDaughtersMax(1.5),
  fProdV0DaughterEtaRange(0.8),
  fProdV0DaughterPtMin(0.0),
  fProdV0DaughterTPCClusterMin(70),
  fProdV0DaughterTPCCrossRatioMin(0.0),
  fProdRfidMinV0(0.6),
  fProdRfidMaxV0(100.0),
  fProdDcaV0ToPrimVertexMin(0.),
  fProdDcaV0PrToPrimVertexMin(0.),
  fProdDcaV0PiToPrimVertexMin(0.),
  fProdV0ProperDecayLengthMax(99999.),
  fProdMassRejK0s(0.),
  fProdV0EtaMin(-9999.),
  fProdV0EtaMax(9999.),
  fProdV0RapMin(-9999.),
  fProdV0RapMax(9999.),
	fProdRoughMassTol(0.25),
	fProdRoughPtMin(0.0),
	fExcludePionTPC(kFALSE),
	fExcludeProtonTPC(kFALSE),
	fExcludeKaonTPC(kFALSE),
	fExcludenSigmaPionTPC(3.),
	fExcludenSigmaProtonTPC(3.),
	fExcludenSigmaKaonTPC(3.),
	fSigmaElectronTPCMin(-9999.),
	fSigmaElectronTPCPtDepPar0(-9999.),
	fSigmaElectronTPCPtDepPar1(-9999.),
	fSigmaElectronTPCPtDepPar2(0.),
	fSigmaElectronTPCMax(9999.),
	fSigmaElectronTOFMin(-9999.),
	fSigmaElectronTOFMax(9999.),
	fConversionMassMax(-1.)
{
  //
  // Default Constructor
  //
  for(Int_t i=0;i<3;i++){
    fPrimVert[i] = 0.;
  }

  const Int_t nvars=4;
  SetNVars(nvars);
  TString varNames[nvars]={
			   "InvMass [GeV/c2]", //0
			   "cos(Opening angle) [cos(rad)]", //1
			   "dPhiS", //2
			   "dEtaS" //3
  };

  Bool_t isUpperCut[nvars]={
			    kTRUE, //0
			    kFALSE, //1
			    kFALSE, //2
			    kFALSE //3
  };
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[nvars]={
			kTRUE, //0
			kTRUE, //1
			kTRUE, //2
			kTRUE, //3
  };
  SetVarsForOpt(nvars,forOpt);

  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);
}
//--------------------------------------------------------------------------
AliRDHFCutsLctoeleLambdafromAODtracks::AliRDHFCutsLctoeleLambdafromAODtracks(const AliRDHFCutsLctoeleLambdafromAODtracks &source) :
  AliRDHFCuts(source),
  fPIDStrategy(source.fPIDStrategy),
  fCombinedPIDThreshold(source.fCombinedPIDThreshold),
  fUseLambdaPID(source.fUseLambdaPID),
  fPidObjProton(source.fPidObjProton),
  fPidObjPion(source.fPidObjPion),
  fUseOnTheFlyV0(source.fUseOnTheFlyV0),
  fUseV0Topology(source.fUseV0Topology),
  fBzkG(source.fBzkG),
  fProdTrackTPCNclsPIDMin(source.fProdTrackTPCNclsPIDMin),
  fProdTrackTPCNclsRatioMin(source.fProdTrackTPCNclsRatioMin),
  fProdUseAODFilterBit(source.fProdUseAODFilterBit),
  fProdAODFilterBit(source.fProdAODFilterBit),
  fProdRejectTrackWithShared(source.fProdRejectTrackWithShared),
  fProdV0KinkRejection(source.fProdV0KinkRejection),
  fProdV0MassTolLambda(source.fProdV0MassTolLambda),
  fProdV0MassTolLambdaRough(source.fProdV0MassTolLambdaRough),
  fProdV0PtMin(source.fProdV0PtMin),
  fProdV0PtMax(source.fProdV0PtMax),
  fProdV0CosPointingAngleToPrimVtxMin(source.fProdV0CosPointingAngleToPrimVtxMin),
  fProdV0DcaDaughtersMax(source.fProdV0DcaDaughtersMax),
  fProdV0DaughterEtaRange(source.fProdV0DaughterEtaRange),
  fProdV0DaughterPtMin(source.fProdV0DaughterPtMin),
  fProdV0DaughterTPCClusterMin(source.fProdV0DaughterTPCClusterMin),
  fProdV0DaughterTPCCrossRatioMin(source.fProdV0DaughterTPCCrossRatioMin),
  fProdRfidMinV0(source.fProdRfidMinV0),
  fProdRfidMaxV0(source.fProdRfidMaxV0),
  fProdDcaV0ToPrimVertexMin(source.fProdDcaV0ToPrimVertexMin),
  fProdDcaV0PrToPrimVertexMin(source.fProdDcaV0PrToPrimVertexMin),
  fProdDcaV0PiToPrimVertexMin(source.fProdDcaV0PiToPrimVertexMin),
  fProdV0ProperDecayLengthMax(source.fProdV0ProperDecayLengthMax),
  fProdMassRejK0s(source.fProdMassRejK0s),
  fProdV0EtaMin(source.fProdV0EtaMin),
  fProdV0EtaMax(source.fProdV0EtaMax),
  fProdV0RapMin(source.fProdV0RapMin),
  fProdV0RapMax(source.fProdV0RapMax),
  fProdRoughMassTol(source.fProdRoughMassTol),
  fProdRoughPtMin(source.fProdRoughPtMin),
	fExcludePionTPC(source.fExcludePionTPC),
	fExcludeProtonTPC(source.fExcludeProtonTPC),
	fExcludeKaonTPC(source.fExcludeKaonTPC),
	fExcludenSigmaPionTPC(source.fExcludenSigmaPionTPC),
	fExcludenSigmaProtonTPC(source.fExcludenSigmaProtonTPC),
	fExcludenSigmaKaonTPC(source.fExcludenSigmaKaonTPC),
	fSigmaElectronTPCMin(source.fSigmaElectronTPCMin),
	fSigmaElectronTPCPtDepPar0(source.fSigmaElectronTPCPtDepPar0),
	fSigmaElectronTPCPtDepPar1(source.fSigmaElectronTPCPtDepPar1),
	fSigmaElectronTPCPtDepPar2(source.fSigmaElectronTPCPtDepPar2),
	fSigmaElectronTPCMax(source.fSigmaElectronTPCMax),
	fSigmaElectronTOFMin(source.fSigmaElectronTOFMin),
	fSigmaElectronTOFMax(source.fSigmaElectronTOFMax),
	fConversionMassMax(source.fConversionMassMax)
{
  //
  // Copy constructor
  //
  for(Int_t i=0;i<3;i++){
    fPrimVert[i] = source.fPrimVert[i];
  }
}
//--------------------------------------------------------------------------
AliRDHFCutsLctoeleLambdafromAODtracks &AliRDHFCutsLctoeleLambdafromAODtracks::operator=(const AliRDHFCutsLctoeleLambdafromAODtracks &source)
{
  //
  // assignment operator
  //

  if (this != &source) {
    AliRDHFCuts::operator=(source);
  }

  fPIDStrategy = source.fPIDStrategy;
  fCombinedPIDThreshold = source.fCombinedPIDThreshold;
  fUseLambdaPID = source.fUseLambdaPID;
  fPidObjProton = source.fPidObjProton;
  fPidObjPion = source.fPidObjPion;
  fUseOnTheFlyV0 = source.fUseOnTheFlyV0;
  fUseV0Topology = source.fUseV0Topology;
  fBzkG = source.fBzkG;
  fProdTrackTPCNclsPIDMin = source.fProdTrackTPCNclsPIDMin;
  fProdTrackTPCNclsRatioMin = source.fProdTrackTPCNclsRatioMin;
  fProdUseAODFilterBit = source.fProdUseAODFilterBit;
  fProdAODFilterBit = source.fProdAODFilterBit;
  fProdRejectTrackWithShared = source.fProdRejectTrackWithShared;
  fProdV0KinkRejection = source.fProdV0KinkRejection;
  fProdV0MassTolLambda = source.fProdV0MassTolLambda;
  fProdV0MassTolLambdaRough = source.fProdV0MassTolLambdaRough;
  fProdV0PtMin = source.fProdV0PtMin;
  fProdV0PtMax = source.fProdV0PtMax;
  fProdV0CosPointingAngleToPrimVtxMin = source.fProdV0CosPointingAngleToPrimVtxMin;
  fProdV0DcaDaughtersMax=source.fProdV0DcaDaughtersMax;
  fProdV0DaughterEtaRange=source.fProdV0DaughterEtaRange;
  fProdV0DaughterPtMin=source.fProdV0DaughterPtMin;
  fProdV0DaughterTPCClusterMin=source.fProdV0DaughterTPCClusterMin;
  fProdV0DaughterTPCCrossRatioMin=source.fProdV0DaughterTPCCrossRatioMin;
  fProdRfidMinV0=source.fProdRfidMinV0;
  fProdRfidMaxV0=source.fProdRfidMaxV0;
  fProdDcaV0ToPrimVertexMin=source.fProdDcaV0ToPrimVertexMin;
  fProdDcaV0PrToPrimVertexMin=source.fProdDcaV0PrToPrimVertexMin;
  fProdDcaV0PiToPrimVertexMin=source.fProdDcaV0PiToPrimVertexMin;
  fProdV0ProperDecayLengthMax=source.fProdV0ProperDecayLengthMax;
  fProdMassRejK0s=source.fProdMassRejK0s;
  fProdV0EtaMin = source.fProdV0EtaMin;
  fProdV0EtaMax = source.fProdV0EtaMax;
  fProdV0RapMin = source.fProdV0RapMin;
  fProdV0RapMax = source.fProdV0RapMax;
  fProdRoughMassTol = source.fProdRoughMassTol;
  fProdRoughPtMin = source.fProdRoughPtMin;
	fExcludePionTPC = source.fExcludePionTPC;
	fExcludeProtonTPC = source.fExcludeProtonTPC;
	fExcludeKaonTPC = source.fExcludeKaonTPC;
	fExcludenSigmaPionTPC = source.fExcludenSigmaPionTPC;
	fExcludenSigmaProtonTPC = source.fExcludenSigmaProtonTPC;
	fExcludenSigmaKaonTPC = source.fExcludenSigmaKaonTPC;
	fSigmaElectronTPCMin = source.fSigmaElectronTPCMin;
	fSigmaElectronTPCPtDepPar0 = source.fSigmaElectronTPCPtDepPar0;
	fSigmaElectronTPCPtDepPar1 = source.fSigmaElectronTPCPtDepPar1;
	fSigmaElectronTPCPtDepPar2 = source.fSigmaElectronTPCPtDepPar2;
	fSigmaElectronTPCMax = source.fSigmaElectronTPCMax;
	fSigmaElectronTOFMin = source.fSigmaElectronTOFMin;
	fSigmaElectronTOFMax = source.fSigmaElectronTOFMax;
	fConversionMassMax = source.fConversionMassMax;

  for(Int_t i=0;i<3;i++){
    fPrimVert[i] = source.fPrimVert[i];
  }
  
  return *this;
}

//---------------------------------------------------------------------------
AliRDHFCutsLctoeleLambdafromAODtracks::~AliRDHFCutsLctoeleLambdafromAODtracks() {
  //
  //  Default Destructor
  //
  
}

//---------------------------------------------------------------------------
void AliRDHFCutsLctoeleLambdafromAODtracks::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
  //
  // Fills in vars the values of the variables
  //
  
  if (pdgdaughters[0]==-9999) return; // dummy

  AliAODRecoCascadeHF* dd=(AliAODRecoCascadeHF*)d;
  if(!dd){
    AliError("No AliAODRecoCascadeHF object found\n");
    return;
  }
  
  if (nvars!=fnVarsForOpt) {
    AliError("Cut object seems to have the wrong number of variables\n");
    return;
  }

  //Double_t ptD=d->Pt();
  //Int_t ptbin=PtBin(ptD);
  Int_t iter=-1;


  if(fVarsForOpt[0]){
    iter++;
		UInt_t pdgdg[2]={11,3122};
    vars[iter]= dd->InvMass(2,pdgdg);
  }
  if(fVarsForOpt[1]){
    iter++;
		Double_t v0px = dd->PxProng(1);
		Double_t v0py = dd->PyProng(1);
		Double_t v0pz = dd->PzProng(1);
		Double_t epx = dd->PxProng(0);
		Double_t epy = dd->PyProng(0);
		Double_t epz = dd->PzProng(0);
    vars[iter]=  (v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz);
  }
  
  return;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctoeleLambdafromAODtracks::IsSelected(TObject* obj, Int_t selectionLevel) {
  //
  // Apply selection
  //

  if (!fCutsRD) {
    AliFatal("Cut matrix not inizialized. Exit...");
    return 0;
  }

  AliAODRecoCascadeHF* d=(AliAODRecoCascadeHF*)obj;
  if(!d){
    AliDebug(2,"AliAODRecoCascadeHF null");
    return 0;
  }

  Double_t ptD=d->Pt();
  if(ptD<fMinPtCand) return 0;
  if(ptD>fMaxPtCand) return 0;

  if (selectionLevel==AliRDHFCuts::kAll ||
      selectionLevel==AliRDHFCuts::kTracks) {
    //Performed in production stage
  }

  Int_t returnvalueCuts=1;
  // selection on candidate
  if (selectionLevel==AliRDHFCuts::kAll ||
      selectionLevel==AliRDHFCuts::kCandidate) {
    
    Double_t pt=d->Pt();
    Int_t ptbin=PtBin(pt);
    if (ptbin==-1) {
      return 0;
    }
    Bool_t okcand=kTRUE;
    
    Double_t mlamPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
		UInt_t pdgdg[2]={11,3122};
		Double_t InvMassEleLambda = d->InvMass(2,pdgdg);

		Double_t v0px = d->PxProng(1);
		Double_t v0py = d->PyProng(1);
		Double_t v0pz = d->PzProng(1);
		Double_t epx = d->PxProng(0);
		Double_t epy = d->PyProng(0);
		Double_t epz = d->PzProng(0);
		Double_t cosoa = (v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz);

    Double_t dphis_e_pr, detas_e_pr, dphis_e_pi, detas_e_pi;
    dphis_e_pr = 9999.;
    detas_e_pr = 9999.;
    dphis_e_pi = 9999.;
    detas_e_pi = 9999.;
    if(fCutsRD[GetGlobalIndex(2,ptbin)]>0 || fCutsRD[GetGlobalIndex(3,ptbin)]>0){
      AliAODTrack *trke = (AliAODTrack*)d->GetDaughter(0);
      AliAODv0 *v0 = (AliAODv0*)d->GetDaughter(1);
      if(trke && v0){
        Bool_t isparticle = kTRUE;
        if(TMath::Abs(v0->MassAntiLambda()-mlamPDG)<fProdV0MassTolLambdaRough) isparticle = kFALSE;
        AliAODTrack *cprtrk = 0;
        AliAODTrack *cpitrk = 0;
        if(isparticle){
          cprtrk = (AliAODTrack*)v0->GetDaughter(0);
          cpitrk = (AliAODTrack*)v0->GetDaughter(1);
        }else{
          cprtrk = (AliAODTrack*)v0->GetDaughter(1);
          cpitrk = (AliAODTrack*)v0->GetDaughter(0);
        }
        if(cprtrk && cpitrk)
          GetdPhiSdEtaSR125(trke,cprtrk,cpitrk,fBzkG,fPrimVert, dphis_e_pr,detas_e_pr,dphis_e_pi,detas_e_pi);
      }
    }
    
    if(InvMassEleLambda > fCutsRD[GetGlobalIndex(0,ptbin)])
      {
	okcand = kFALSE;
      }
    if(cosoa < fCutsRD[GetGlobalIndex(1,ptbin)])
      {
	okcand = kFALSE;
      }
    if(fabs(dphis_e_pr) < fCutsRD[GetGlobalIndex(2,ptbin)] && fabs(detas_e_pr) < fCutsRD[GetGlobalIndex(3,ptbin)])
      {
	okcand = kFALSE;
      }
    if(fabs(dphis_e_pi) < fCutsRD[GetGlobalIndex(2,ptbin)] && fabs(detas_e_pi) < fCutsRD[GetGlobalIndex(3,ptbin)])
      {
	okcand = kFALSE;
      }
    
    if(!okcand)  return 0;
    returnvalueCuts = 1;
  }
  
  Int_t returnvaluePID=1;
  if(selectionLevel==AliRDHFCuts::kAll ||
     selectionLevel==AliRDHFCuts::kCandidate|| 
     selectionLevel==AliRDHFCuts::kPID) {
		//Not used
  }
  
  Int_t returnvalue = 0;
  if(returnvalueCuts==1 && returnvaluePID==1) returnvalue=1;
  
  return returnvalue;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctoeleLambdafromAODtracks::IsSelectedPID(AliAODRecoDecayHF* obj) 
{
  //
  // IsSelectedPID (not used)
  //

  if(!fUsePID || !obj) return 1;

  AliAODRecoCascadeHF* dd=(AliAODRecoCascadeHF*)obj;
  AliAODTrack *part = dd->GetBachelor();
  
  Int_t returnvalue=1;

  if(fPidHF->GetPidResponse()==0x0){
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
    AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
    fPidHF->SetPidResponse(pidResp);
  }

  Int_t isElectron=fPidHF->MakeRawPid(part,0); 
  if(isElectron<1) returnvalue = 0;
  
  return returnvalue;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctoeleLambdafromAODtracks::IsSelectedCombinedPID(AliAODRecoDecayHF* obj) {
  //
  // IsSelectedCombinedPID (not used)
  //
    
  if(!fUsePID || !obj) {return 1;}

  AliAODRecoCascadeHF* dd=(AliAODRecoCascadeHF*)obj;
  AliAODTrack *part = dd->GetBachelor();
  if(!part) return 0;

  Int_t returnvalue=1;

  return returnvalue;
}

//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::SingleTrkCuts(AliAODTrack *trk, AliAODTrack *trkpid, AliAODVertex *primVert)
{
  //
  // Single Track Cut to be applied before object creation
  //

  //if(trk->GetStatus()&AliESDtrack::kITSpureSA) return kFALSE;
  //if(!(trk->GetStatus()&AliESDtrack::kITSin)) return kFALSE;
  //if(fProdUseAODFilterBit && !trk->TestFilterMask(BIT(fProdAODFilterBit))) return kFALSE;
  if(fTrackCuts){
    if(fProdAODFilterBit==7){
      Float_t ptmin, ptmax, etamin, etamax;
      fTrackCuts->GetPtRange(ptmin,ptmax);
      fTrackCuts->GetEtaRange(etamin,etamax);
      if(trk->Pt()<ptmin || trk->Pt()>ptmax) return kFALSE;
      if(trk->Eta()<etamin || trk->Eta()>etamax) return kFALSE;
    }else{
      Double_t pos[3]; primVert->GetXYZ(pos);
      Double_t cov[6]; primVert->GetCovarianceMatrix(cov);
      const AliESDVertex vESD(pos,cov,100.,100);
      if(!IsDaughterSelected(trk,&vESD,fTrackCuts)) return kFALSE;
    }
  }

	if(trkpid->GetTPCsignalN()<fProdTrackTPCNclsPIDMin) return kFALSE;
	if(trk->GetTPCNclsF()>0){
		Float_t tpcratio = (Float_t)trk->GetTPCncls()/(Float_t)trk->GetTPCNclsF();
		if(tpcratio<fProdTrackTPCNclsRatioMin) return kFALSE;
	}

  if(fProdRejectTrackWithShared){
    const TBits sharedMap = trk->GetTPCSharedMap();
    if((sharedMap.CountBits()) >= 1){
      return kFALSE;
    }
  }

  if(fUsePID)
  {
    if(fPidHF->GetPidResponse()==0x0){
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
      AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
      fPidHF->SetPidResponse(pidResp);
    }

    switch(fPIDStrategy){
      case kNSigmaCuts:
        return IsSelectedeID(trkpid);
        break;
      case kNSigmaCustomizedCuts:
        return IsSelectedCustomizedeID(trkpid);
        break;
      case kNSigmaCustomizedPtDepCuts:
        return IsSelectedCustomizedPtDepeID(trk,trkpid);
        break;
      case kCombinedCuts:
        return IsSelectedCombinedeID(trkpid);
        break;
    }

  }

  return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::SingleTrkCutsNoPID(AliAODTrack *trk, AliAODTrack *trkpid, AliAODVertex *primVert)
{
  //
  // Single Track Cut to be applied before object creation
  //

  //if(trk->GetStatus()&AliESDtrack::kITSpureSA) return kFALSE;
  //if(!(trk->GetStatus()&AliESDtrack::kITSin)) return kFALSE;
  //if(fProdUseAODFilterBit && !trk->TestFilterMask(BIT(4))) return kFALSE;
  if(fTrackCuts){
    if(fProdAODFilterBit==7){
      Float_t ptmin, ptmax, etamin, etamax;
      fTrackCuts->GetPtRange(ptmin,ptmax);
      fTrackCuts->GetEtaRange(etamin,etamax);
      if(trk->Pt()<ptmin || trk->Pt()>ptmax) return kFALSE;
      if(trk->Eta()<etamin || trk->Eta()>etamax) return kFALSE;
    }else{
      Double_t pos[3]; primVert->GetXYZ(pos);
      Double_t cov[6]; primVert->GetCovarianceMatrix(cov);
      const AliESDVertex vESD(pos,cov,100.,100);
      if(!IsDaughterSelected(trk,&vESD,fTrackCuts)) return kFALSE;
    }
  }

	if(trkpid->GetTPCsignalN()<fProdTrackTPCNclsPIDMin) return kFALSE;
	if(trk->GetTPCNclsF()>0){
		Float_t tpcratio = (Float_t)trk->GetTPCncls()/(Float_t)trk->GetTPCNclsF();
		if(tpcratio<fProdTrackTPCNclsRatioMin) return kFALSE;
	}

  if(fProdRejectTrackWithShared){
    const TBits sharedMap = trk->GetTPCSharedMap();
    if((sharedMap.CountBits()) >= 1){
      return kFALSE;
    }
  }

  return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::IsSelectedeID(AliAODTrack *trk)
{
  //
  // electron ID first shot
  //

	Double_t nSigmaPion = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kPion);
	if(fExcludePionTPC){
		if(TMath::Abs(nSigmaPion)<fExcludenSigmaPionTPC){
			return kFALSE;
		}
	}

	Double_t nSigmaProton = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kProton);
	if(fExcludeProtonTPC){
		if(TMath::Abs(nSigmaProton)<fExcludenSigmaProtonTPC){
			return kFALSE;
		}
	}

	Double_t nSigmaKaon = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kKaon);
	if(fExcludeKaonTPC){
		if(TMath::Abs(nSigmaKaon)<fExcludenSigmaKaonTPC){
			return kFALSE;
		}
	}

  Int_t isElectron=fPidHF->MakeRawPid(trk,0); 
  if(isElectron<1) return kFALSE;

	return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::IsSelectedCustomizedeID(AliAODTrack *trk)
{
  //
  // electron ID first shot
  //

	Double_t nSigmaPion = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kPion);
	if(fExcludePionTPC){
		if(TMath::Abs(nSigmaPion)<fExcludenSigmaPionTPC){
			return kFALSE;
		}
	}

	Double_t nSigmaProton = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kProton);
	if(fExcludeProtonTPC){
		if(TMath::Abs(nSigmaProton)<fExcludenSigmaProtonTPC){
			return kFALSE;
		}
	}

	Double_t nSigmaKaon = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kKaon);
	if(fExcludeKaonTPC){
		if(TMath::Abs(nSigmaKaon)<fExcludenSigmaKaonTPC){
			return kFALSE;
		}
	}

	Double_t nSigmaTPCele = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kElectron);
	Double_t nSigmaTOFele = fPidHF->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kElectron);

	if(nSigmaTPCele<fSigmaElectronTPCMin) return kFALSE;
	if(nSigmaTPCele>fSigmaElectronTPCMax) return kFALSE;
	if(nSigmaTOFele<fSigmaElectronTOFMin) return kFALSE;
	if(nSigmaTOFele>fSigmaElectronTOFMax) return kFALSE;

	return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::IsSelectedCustomizedPtDepeID(AliAODTrack *trk, AliAODTrack *trkpid)
{
  //
  // electron ID first shot
  //

	Double_t nSigmaPion = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kPion);
	if(fExcludePionTPC){
		if(TMath::Abs(nSigmaPion)<fExcludenSigmaPionTPC){
			return kFALSE;
		}
	}

	Double_t nSigmaProton = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kProton);
	if(fExcludeProtonTPC){
		if(TMath::Abs(nSigmaProton)<fExcludenSigmaProtonTPC){
			return kFALSE;
		}
	}

	Double_t nSigmaKaon = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kKaon);
	if(fExcludeKaonTPC){
		if(TMath::Abs(nSigmaKaon)<fExcludenSigmaKaonTPC){
			return kFALSE;
		}
	}

	Double_t nSigmaTPCele = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kElectron);
	Double_t nSigmaTOFele = fPidHF->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kElectron);

  if(fabs(fSigmaElectronTOFMin)<999.|| fabs(fSigmaElectronTOFMax)<999.)
  {
    if(nSigmaTOFele<fSigmaElectronTOFMin) return kFALSE;
    if(nSigmaTOFele>fSigmaElectronTOFMax) return kFALSE;
  }

	Double_t pte = trk->Pt();
	Double_t nsigmamin = fSigmaElectronTPCPtDepPar0+fSigmaElectronTPCPtDepPar1*pte+fSigmaElectronTPCPtDepPar2*pte*pte;
	if(pte>5.) nsigmamin = fSigmaElectronTPCPtDepPar0+fSigmaElectronTPCPtDepPar1*5.+fSigmaElectronTPCPtDepPar2*25.;
	if(nSigmaTPCele<nsigmamin) return kFALSE;
	if(nSigmaTPCele>fSigmaElectronTPCMax) return kFALSE;

	return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::IsSelectedCombinedeID(AliAODTrack *trk)
{
  //
  // electron ID Basyian not implemented
  //

	return kTRUE;
}

//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::SingleV0Cuts(AliAODv0 *v0, AliAODVertex *primVert)
{
  //
  // Single V0 Cut to be applied before object creation
  //

  Bool_t onFlyV0 = v0->GetOnFlyStatus(); // on-the-flight V0s
  if ( onFlyV0 && !fUseOnTheFlyV0 ) return kFALSE;
  if ( !onFlyV0 && fUseOnTheFlyV0 ) return kFALSE;
	if(!v0) return kFALSE;
	if(!(v0->GetSecondaryVtx())) return kFALSE;

  AliAODTrack *cptrack =  (AliAODTrack*)(v0->GetDaughter(0));
  AliAODTrack *cntrack =  (AliAODTrack*)(v0->GetDaughter(1));
  if(!cptrack || !cntrack) return kFALSE;
//	if(cptrack->Charge()<0 && cntrack->Charge()>0){
//		//In case sign is wrong
//		cptrack =  (AliAODTrack*)(v0->GetDaughter(1));
//		cntrack =  (AliAODTrack*)(v0->GetDaughter(0));
//	}

  if ( cptrack->Charge() == cntrack->Charge() ) return kFALSE;
  if(!(cptrack->GetStatus() & AliESDtrack::kTPCrefit) ||
     !(cntrack->GetStatus() & AliESDtrack::kTPCrefit)) return kFALSE;
  AliAODVertex *maybeKinkPos = (AliAODVertex*)cptrack->GetProdVertex();
  AliAODVertex *maybeKinkNeg = (AliAODVertex*)cntrack->GetProdVertex();
  if (fProdV0KinkRejection && (maybeKinkPos->GetType()==AliAODVertex::kKink || maybeKinkNeg->GetType()==AliAODVertex::kKink)) 
    return kFALSE;

  if ( ( ( cptrack->GetTPCClusterInfo(2,1) ) < (Float_t)fProdV0DaughterTPCClusterMin ) || 
       ( ( cntrack->GetTPCClusterInfo(2,1) ) < (Float_t)fProdV0DaughterTPCClusterMin) ) return kFALSE;


	if((cptrack->GetTPCNclsF()>0)&&(cptrack->GetTPCClusterInfo(2,1)/cptrack->GetTPCNclsF()<fProdV0DaughterTPCCrossRatioMin)) return kFALSE;
	if((cntrack->GetTPCNclsF()>0)&&(cntrack->GetTPCClusterInfo(2,1)/cntrack->GetTPCNclsF()<fProdV0DaughterTPCCrossRatioMin)) return kFALSE;

  Double_t massLambda = v0->MassLambda();
  Double_t massAntiLambda = v0->MassAntiLambda();
  Double_t mlamPDG   = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  Double_t mk0sPDG   = TDatabasePDG::Instance()->GetParticle(310)->Mass();
  if(!(fabs(massAntiLambda-mlamPDG)<fProdV0MassTolLambdaRough) && !(fabs(massLambda-mlamPDG)<fProdV0MassTolLambdaRough)) return kFALSE;

  Double_t massK0s = v0->MassK0Short();
  if(TMath::Abs(massK0s-mk0sPDG)<fProdMassRejK0s)
    return kFALSE;

  Bool_t isparticle = kTRUE;
  if(TMath::Abs(massAntiLambda-mlamPDG)<fProdV0MassTolLambdaRough) isparticle = kFALSE;


  if(TMath::Abs(v0->DcaV0Daughters())>fProdV0DcaDaughtersMax) return kFALSE;
  Double_t posVtx[3] = {0.,0.,0.};
  primVert->GetXYZ(posVtx);
  Double_t cospav0 = v0->CosPointingAngle(posVtx); 
  if(cospav0<fProdV0CosPointingAngleToPrimVtxMin) return kFALSE;
  if(v0->Pt()<fProdV0PtMin) return kFALSE;
  if(v0->Pt()>fProdV0PtMax) return kFALSE;
  if(fabs(cptrack->Eta())>fProdV0DaughterEtaRange) return kFALSE;
  if(fabs(cntrack->Eta())>fProdV0DaughterEtaRange) return kFALSE;
  if(cptrack->Pt()<fProdV0DaughterPtMin) return kFALSE;
  if(cntrack->Pt()<fProdV0DaughterPtMin) return kFALSE;

	Double_t lDcaV0ToPrimVertex = v0->DcaV0ToPrimVertex();
	Double_t lDcaPosToPrimVertex = v0->DcaPosToPrimVertex();
	Double_t lDcaNegToPrimVertex = v0->DcaNegToPrimVertex();
	if(lDcaV0ToPrimVertex < fProdDcaV0ToPrimVertexMin) return kFALSE;
  if(isparticle){
    if(lDcaPosToPrimVertex < fProdDcaV0PrToPrimVertexMin) return kFALSE;
    if(lDcaNegToPrimVertex < fProdDcaV0PiToPrimVertexMin) return kFALSE;
  }else{
    if(lDcaPosToPrimVertex < fProdDcaV0PiToPrimVertexMin) return kFALSE;
    if(lDcaNegToPrimVertex < fProdDcaV0PrToPrimVertexMin) return kFALSE;
  }
  Double_t lPosV0[3];
  lPosV0[0] = v0->DecayVertexV0X();
  lPosV0[1] = v0->DecayVertexV0Y();
  lPosV0[2] = v0->DecayVertexV0Z();
  Double_t decayvertV0 = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
  
  if(decayvertV0<fProdRfidMinV0 || decayvertV0>fProdRfidMaxV0) return kFALSE;

  Double_t ptotlam = TMath::Sqrt(pow(v0->Px(),2)+pow(v0->Py(),2)+pow(v0->Pz(),2));
  Double_t properdl = v0->DecayLengthV0(posVtx)*mlamPDG/ptotlam;
  if(properdl>fProdV0ProperDecayLengthMax) return kFALSE;


  if(fUseLambdaPID)
    {
			if(fPidObjProton->GetPidResponse()==0x0){
				AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
				AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
				AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
				fPidObjProton->SetPidResponse(pidResp);
			}
			if(fPidObjPion->GetPidResponse()==0x0){
				AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
				AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
				AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
				fPidObjPion->SetPidResponse(pidResp);
			}

      Int_t isProton= 1;
      Int_t isPion= 1;
      Double_t nsigmatpc_proton = fPidObjProton->GetSigma(0);
      Double_t nsigmatpc_pion = fPidObjPion->GetSigma(0);
      Double_t nsigmatof_proton = fPidObjProton->GetSigma(3);
      Double_t nsigmatof_pion = fPidObjPion->GetSigma(3);

      if(isparticle){
        Double_t nSigmaTPCpr = fPidObjProton->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kProton);
        Double_t nSigmaTPCpi = fPidObjPion->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kPion);
        Double_t nSigmaTOFpr = fPidObjProton->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kProton);
        Double_t nSigmaTOFpi = fPidObjPion->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kPion);
        if(fabs(nSigmaTPCpr)>nsigmatpc_proton) isProton = 0;
        if(fabs(nSigmaTPCpi)>nsigmatpc_pion) isPion = 0;
        if(nsigmatof_proton>0.01 && nsigmatof_proton<999.){
          if(fabs(nSigmaTOFpr)>nsigmatof_proton) isProton = 0;
        }
        if(nsigmatof_pion>0.01 && nsigmatof_pion<999.){
          if(fabs(nSigmaTOFpi)>nsigmatof_pion) isPion = 0;
        }
      }else{
        Double_t nSigmaTPCpr = fPidObjProton->GetPidResponse()->NumberOfSigmasTPC(cntrack,AliPID::kProton);
        Double_t nSigmaTPCpi = fPidObjPion->GetPidResponse()->NumberOfSigmasTPC(cptrack,AliPID::kPion);
        Double_t nSigmaTOFpr = fPidObjProton->GetPidResponse()->NumberOfSigmasTOF(cntrack,AliPID::kProton);
        Double_t nSigmaTOFpi = fPidObjPion->GetPidResponse()->NumberOfSigmasTOF(cptrack,AliPID::kPion);
        if(fabs(nSigmaTPCpr)>nsigmatpc_proton) isProton = 0;
        if(fabs(nSigmaTPCpi)>nsigmatpc_pion) isPion = 0;
        if(nsigmatof_proton>0.01 && nsigmatof_proton<999.){
          if(fabs(nSigmaTOFpr)>nsigmatof_proton) isProton = 0;
        }
        if(nsigmatof_pion>0.01 && nsigmatof_pion<999.){
          if(fabs(nSigmaTOFpi)>nsigmatof_pion) isPion = 0;
        }
      }
      if(isProton<1) return kFALSE;
      if(isPion<1) return kFALSE;
    }

	Double_t RapLambda = v0->RapLambda();
	if(RapLambda<fProdV0RapMin || RapLambda>fProdV0RapMax) return kFALSE;

	Double_t EtaLambda = v0->PseudoRapV0();
	if(EtaLambda<fProdV0EtaMin || EtaLambda>fProdV0EtaMax) return kFALSE;

  if(fProdRejectTrackWithShared){
    const TBits sharedMap1 = cptrack->GetTPCSharedMap();
    const TBits sharedMap2 = cntrack->GetTPCSharedMap();
    if((sharedMap1.CountBits() >= 1) ||
        (sharedMap2.CountBits() >= 1))
    {
      return kFALSE;
    }
  }

  if(fUseV0Topology>0){
    Double_t dphiprlam = 0.;
    if(isparticle){
      TVector3 v3v0pr(v0->MomPosX(),v0->MomPosY(),v0->MomPosZ());
      TVector3 v3lam(v0->Px(),v0->Py(),v0->Pz());
      dphiprlam = v3v0pr.DeltaPhi(v3lam);
    }else{
      TVector3 v3v0pr(v0->MomNegX(),v0->MomNegY(),v0->MomNegZ());
      TVector3 v3lam(v0->Px(),v0->Py(),v0->Pz());
      dphiprlam = -1.*v3v0pr.DeltaPhi(v3lam);
    }
    if(fUseV0Topology==1){
      if(dphiprlam>0) return kFALSE;
    }else if(fUseV0Topology==2){
      if(dphiprlam<0) return kFALSE;
    }
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::SelectWithRoughCuts(AliAODv0 *v0, AliAODTrack *part)
{
  //
  // Mass and pT Cut to be applied before object creation
	// Not used now
  //
	if(!v0) return kFALSE;
	if(!part) return kFALSE;

	return kTRUE;

}

//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::IsPeakRegion(AliAODv0 *v0)
{
  Double_t mLPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  Double_t massL = v0->MassLambda();
  if(TMath::Abs(massL-mLPDG)>fProdV0MassTolLambdaRough)
		massL = v0->MassAntiLambda();

  if(TMath::Abs(massL-mLPDG)<fProdV0MassTolLambda)
		return kTRUE;
	return kFALSE;
}

//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::IsPeakRegion(TLorentzVector *v0)
{
  Double_t mLPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  Double_t massL = v0->M();
  if(TMath::Abs(massL-mLPDG)<fProdV0MassTolLambda)
		return kTRUE;
	return kFALSE;
}

//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::IsSideBand(AliAODv0 *v0)
{
  Double_t mLPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  Double_t massL = v0->MassLambda();
  if(TMath::Abs(massL-mLPDG)>fProdV0MassTolLambdaRough)
		massL = v0->MassAntiLambda();

	Bool_t issideband = kFALSE;
  if((massL-mLPDG)>fProdV0MassTolLambdaRough-fProdV0MassTolLambda) issideband = kTRUE;
  if((massL-mLPDG)<-fProdV0MassTolLambdaRough+fProdV0MassTolLambda) issideband = kTRUE;
	return issideband;
}

//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::IsSideBand(TLorentzVector *v0)
{
  Double_t mLPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  Double_t massL = v0->M();
	Bool_t issideband = kFALSE;
  if((massL-mLPDG)>fProdV0MassTolLambdaRough-fProdV0MassTolLambda) issideband = kTRUE;
  if((massL-mLPDG)<-fProdV0MassTolLambdaRough+fProdV0MassTolLambda) issideband = kTRUE;
	return issideband;
}

//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::TagConversions(AliAODTrack *etrk, Int_t *id2index, AliAODEvent *evt, Int_t ntrk, Double_t &minmass)
{
  //
  // Tag conversion electron tracks
  //
	if(fConversionMassMax<0.) return kFALSE;

	Bool_t isconv = kFALSE;
	minmass = 9999.;
	Int_t trkid = etrk->GetID();

	Double_t px1 = etrk->Px();
	Double_t py1 = etrk->Py();
	Double_t pz1 = etrk->Pz();
	Double_t E1 = sqrt(px1*px1+py1*py1+pz1*pz1+0.000511*0.000511);

	for(Int_t it=0;it<ntrk;it++){
		AliAODTrack *trk2 = (AliAODTrack*) evt->GetTrack(it);
		if(!trk2) continue;
		Int_t trkid2 = trk2->GetID();
		if(abs(trkid)==abs(trkid2)) continue;
		if(etrk->Charge()*trk2->Charge()>0) continue;
		if(fProdAODFilterBit==7){
      if(!trk2->TestFilterBit(BIT(fProdAODFilterBit))) continue;
    }else{
      if(!trk2->TestFilterMask(BIT(fProdAODFilterBit))) continue;
    }

    Double_t nSigmaTPCele = 9999.;
    if(fProdAODFilterBit==7){
      if(-trkid2-1>=19000) continue;
      if(-trkid2-1<0) continue;
      Int_t index = id2index[-trkid2-1];
      AliAODTrack *partpid = (AliAODTrack*)evt->GetTrack(index);
      nSigmaTPCele = fPidHF->GetPidResponse()->NumberOfSigmasTPC(partpid,AliPID::kElectron);
    }else{
      nSigmaTPCele = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trk2,AliPID::kElectron);
    }
    if(fabs(nSigmaTPCele)>5.) continue;

		Double_t px2 = trk2->Px();
		Double_t py2 = trk2->Py();
		Double_t pz2 = trk2->Pz();
		Double_t E2 = sqrt(px2*px2+py2*py2+pz2*pz2+0.000511*0.000511);

		Double_t mass = sqrt(pow(E1+E2,2)-pow(px1+px2,2)-pow(py1+py2,2)-pow(pz1+pz2,2));
    //Double_t mass =  CalculatePhotonMass(etrk,trk2);
		if(mass<minmass) minmass = mass;
	}

	if(minmass<fConversionMassMax) isconv = kTRUE;

  return isconv;
}

//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::TagConversionsSameSign(AliAODTrack *etrk, Int_t *id2index, AliAODEvent *evt, Int_t ntrk, Double_t &minmass)
{
  //
  // Tag conversion electron tracks
  //
	if(fConversionMassMax<0.) return kFALSE;

	Bool_t isconv = kFALSE;
	minmass = 9999.;

	Int_t trkid = etrk->GetID();
	Double_t px1 = etrk->Px();
	Double_t py1 = etrk->Py();
	Double_t pz1 = etrk->Pz();
	Double_t E1 = sqrt(px1*px1+py1*py1+pz1*pz1+0.000511*0.000511);

	for(Int_t it=0;it<ntrk;it++){
		AliAODTrack *trk2 = (AliAODTrack*) evt->GetTrack(it);
		if(!trk2) continue;
		Int_t trkid2 = trk2->GetID();
		if(abs(trkid)==abs(trkid2)) continue;
		if(etrk->Charge()*trk2->Charge()<0) continue;
		if(fProdAODFilterBit==7){
      if(!trk2->TestFilterBit(BIT(fProdAODFilterBit))) continue;
    }else{
      if(!trk2->TestFilterMask(BIT(fProdAODFilterBit))) continue;
    }

    Double_t nSigmaTPCele = 9999.;
    if(fProdAODFilterBit==7){
      if(-trkid2-1>=19000) continue;
      if(-trkid2-1<0) continue;
      Int_t index = id2index[-trkid2-1];
      AliAODTrack *partpid = (AliAODTrack*)evt->GetTrack(index);
      nSigmaTPCele = fPidHF->GetPidResponse()->NumberOfSigmasTPC(partpid,AliPID::kElectron);
    }else{
      nSigmaTPCele = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trk2,AliPID::kElectron);
    }
    if(fabs(nSigmaTPCele)>5.) continue;

		Double_t px2 = trk2->Px();
		Double_t py2 = trk2->Py();
		Double_t pz2 = trk2->Pz();
		Double_t E2 = sqrt(px2*px2+py2*py2+pz2*pz2+0.000511*0.000511);

		Double_t mass = sqrt(pow(E1+E2,2)-pow(px1+px2,2)-pow(py1+py2,2)-pow(pz1+pz2,2));
		if(mass<minmass) minmass = mass;
	}

	if(minmass<fConversionMassMax) isconv = kTRUE;

  return isconv;
}

//________________________________________________________________________
void AliRDHFCutsLctoeleLambdafromAODtracks::SetSftPosR125(AliAODTrack *track,Double_t bfield,Double_t priVtx[3],Double_t *XSftR125)
{
  //
  // Sets the spatial position of the track at the radius R=1.25m in the shifted coordinate system
  //
  
  // Initialize the array to something indicating there was no propagation
  XSftR125[0]=-9999.;
  XSftR125[1]=-9999.;
  XSftR125[2]=-9999.;
   // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);
  
  // The global position of the the track
  Double_t xyz[3]={-9999.,-9999.,-9999.};  

  // The radius we want to propagate to, squared
  const Float_t RSquaredWanted(125.*125.);


  // Propagation is done in local x of the track
  for (Float_t x = 58.;x<247.;x+=1.){
    // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
    // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
    // the track is straight, i.e. has inifinite pt and doesn't get bent. 
    // If the track's momentum is smaller than infinite, it will develop a y-component,
    // which adds to the global radius
    // We don't change the propagation steps to not mess up things!

    // Stop if the propagation was not succesful. This can happen for low pt tracks
    // that don't reach outer radii
    if(!etp.PropagateTo(x,(Float_t)bfield))break;
    etp.GetXYZ(xyz); // GetXYZ returns global coordinates

    // Calculate the shifted radius we are at, squared. 
    // Compare squared radii for faster code
    Float_t shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
                                 + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);

    // Roughly reached the radius we want
    if(shiftedRadiusSquared > RSquaredWanted){
      
      // Bigger loop has bad precision, we're nearly one centimeter too far, 
      // go back in small steps.
      while (shiftedRadiusSquared>RSquaredWanted){
	// Propagate a mm inwards
	x-=.1;
	if(!etp.PropagateTo(x,bfield)){
	  // Propagation failed but we're already with a
	  // cm precision at R=1.25m so we only break the 
	  // inner loop
	  break;
	}
	// Get the global position
	etp.GetXYZ(xyz);
	// Calculate shifted radius, squared
	shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
	                     + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);
      }
      // We reached R=1.25m with a precission of a cm to a mm,
      // set the spatial position
      XSftR125[0]=xyz[0]-priVtx[0];
      XSftR125[1]=xyz[1]-priVtx[1];
      XSftR125[2]=xyz[2]-priVtx[2];
      // Done
      return;
    } // End of if roughly reached radius
  } // End of coarse propagation loop
}

//________________________________________________________________________
void AliRDHFCutsLctoeleLambdafromAODtracks::SetSftPosR(AliAODTrack *track,Double_t bfield,Double_t R, Double_t priVtx[3],Double_t *XSftR)
{
  //
  // Sets the spatial position of the track at the radius R in the shifted coordinate system
  //
  
  // Initialize the array to something indicating there was no propagation
  XSftR[0]=-9999.;
  XSftR[1]=-9999.;
  XSftR[2]=-9999.;
  if(R<3.5) return;

   // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);
  
  // The global position of the the track
  Double_t xyz[3]={-9999.,-9999.,-9999.};  

  // The radius we want to propagate to, squared
  const Float_t RSquaredWanted(R*R);


  // Propagation is done in local x of the track
  for (Float_t x = 2.;x<247.;x+=1.){
    // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
    // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
    // the track is straight, i.e. has inifinite pt and doesn't get bent. 
    // If the track's momentum is smaller than infinite, it will develop a y-component,
    // which adds to the global radius
    // We don't change the propagation steps to not mess up things!

    // Stop if the propagation was not succesful. This can happen for low pt tracks
    // that don't reach outer radii
    if(!etp.PropagateTo(x,(Float_t)bfield))break;
    etp.GetXYZ(xyz); // GetXYZ returns global coordinates

    // Calculate the shifted radius we are at, squared. 
    // Compare squared radii for faster code
    Float_t shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
                                 + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);

    // Roughly reached the radius we want
    if(shiftedRadiusSquared > RSquaredWanted){
      
      // Bigger loop has bad precision, we're nearly one centimeter too far, 
      // go back in small steps.
      while (shiftedRadiusSquared>RSquaredWanted){
	// Propagate a mm inwards
	x-=.1;
	if(!etp.PropagateTo(x,bfield)){
	  // Propagation failed but we're already with a
	  // cm precision at R=1.25m so we only break the 
	  // inner loop
	  break;
	}
	// Get the global position
	etp.GetXYZ(xyz);
	// Calculate shifted radius, squared
	shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
	                     + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);
      }
      // We reached R=1.25m with a precission of a cm to a mm,
      // set the spatial position
      XSftR[0]=xyz[0]-priVtx[0];
      XSftR[1]=xyz[1]-priVtx[1];
      XSftR[2]=xyz[2]-priVtx[2];
      // Done
      return;
    } // End of if roughly reached radius
  } // End of coarse propagation loop
}

//________________________________________________________________________
Double_t AliRDHFCutsLctoeleLambdafromAODtracks::dEtaSR125(Double_t *postrack1,Double_t *postrack2)
{
  //
  // Returns the pseudorapidity star difference
  //

  // It is important to keep the calculations easy and separated.
  // The calculation of EtaS is straight forward, one just has to
  // do it step by step to not get confused.
  Double_t ThetaS1 = TMath::Pi()/2. - TMath::ATan(postrack1[2]/125.);
  Double_t ThetaS2 = TMath::Pi()/2. - TMath::ATan(postrack2[2]/125.);
  Double_t EtaS1 =  -TMath::Log( TMath::Tan(ThetaS1/2.) );
  Double_t EtaS2 =  -TMath::Log( TMath::Tan(ThetaS2/2.) );

  return EtaS1-EtaS2;
}

//________________________________________________________________________
Double_t AliRDHFCutsLctoeleLambdafromAODtracks::dPhiSR125(Double_t *postrack1,Double_t *postrack2)
{
  //
  // returns delta phi star at R=1.25m
  // position5 at R=1.2m is stored as second radius
  //
  //
  Double_t distSft= TMath::Sqrt(TMath::Power(postrack1[0] - postrack2[0],2)
				     +TMath::Power(postrack1[1] - postrack2[1],2));
  return 2.0 * TMath::ATan(distSft/2./(125.));
}

//________________________________________________________________________
Double_t AliRDHFCutsLctoeleLambdafromAODtracks::GetdPhiSdEtaSR125(AliAODTrack *tracke, AliAODTrack *trackp,
    AliAODTrack *trackn, Double_t bfield, Double_t priVtx[3], Double_t &dPhiS_ep, Double_t &dEtaS_ep, 
    Double_t &dPhiS_en, Double_t &dEtaS_en)
{
  //
  // Returns dPhi and dEta at R 125
  //
  Double_t XSftR125_e[3];
  SetSftPosR125(tracke,bfield,priVtx, XSftR125_e);
  Double_t XSftR125_p[3];
  SetSftPosR125(trackp,bfield,priVtx, XSftR125_p);
  Double_t XSftR125_n[3];
  SetSftPosR125(trackn,bfield,priVtx, XSftR125_n);
  dPhiS_ep = dPhiSR125(XSftR125_e,XSftR125_p);
  dEtaS_ep = dEtaSR125(XSftR125_e,XSftR125_p);
  dPhiS_en = dPhiSR125(XSftR125_e,XSftR125_n);
  dEtaS_en = dEtaSR125(XSftR125_e,XSftR125_n);
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctoeleLambdafromAODtracks::IsSelected(TLorentzVector* vtrk, TLorentzVector *vv0, Double_t *cutvars, Int_t selectionLevel) 
{
  //
  // Apply selection on mixed event tracks
  //

  if (!fCutsRD) {
    AliFatal("Cut matrix not inizialized. Exit...");
    return 0;
  }

  Double_t ptD=sqrt(pow(vtrk->Px()+vv0->Px(),2)+pow(vtrk->Px()+vv0->Py(),2));
  if(ptD<fMinPtCand) return 0;
  if(ptD>fMaxPtCand) return 0;

  Double_t pt=ptD;
  Int_t ptbin=PtBin(pt);
  if (ptbin==-1) {
    return 0;
  }

  if (selectionLevel==AliRDHFCuts::kAll ||
      selectionLevel==AliRDHFCuts::kTracks) {
    //Performed in production stage
  }

  Int_t returnvalueCuts=1;
  // selection on candidate
  if (selectionLevel==AliRDHFCuts::kAll ||
      selectionLevel==AliRDHFCuts::kCandidate) {
    
    Bool_t okcand=kTRUE;
    
    Double_t mlcPDG =  TDatabasePDG::Instance()->GetParticle(4122)->Mass();
    Double_t melePDG =  TDatabasePDG::Instance()->GetParticle(11)->Mass();
    Double_t mlamPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
		Double_t v0px = vv0->Px();
		Double_t v0py = vv0->Py();
		Double_t v0pz = vv0->Pz();
		Double_t epx = vtrk->Px();
		Double_t epy = vtrk->Py();
		Double_t epz = vtrk->Pz();
		Double_t cosoa = (v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz);

    TLorentzVector vele, vlam,vlc;
    vele.SetXYZM(epx,epy,epz,melePDG);
    vlam.SetXYZM(v0px,v0py,v0pz,mlamPDG);
    vlc = vele + vlam;

    Double_t dphis_e_pr = 9999.;
    Double_t detas_e_pr = 9999.;
    Double_t dphis_e_pi = 9999.;
    Double_t detas_e_pi = 9999.;
    if(fCutsRD[GetGlobalIndex(2,ptbin)]>0 || fCutsRD[GetGlobalIndex(3,ptbin)]>0){
      Double_t xyzR125_e[3], xyzR125_pr[3], xyzR125_pi[3];
      xyzR125_e[0] = cutvars[0];
      xyzR125_e[1] = cutvars[1];
      xyzR125_e[2] = cutvars[2];
      xyzR125_pr[0] = cutvars[3];
      xyzR125_pr[1] = cutvars[4];
      xyzR125_pr[2] = cutvars[5];
      xyzR125_pi[0] = cutvars[6];
      xyzR125_pi[1] = cutvars[7];
      xyzR125_pi[2] = cutvars[8];
      dphis_e_pr = dPhiSR125(xyzR125_e,xyzR125_pr);
      detas_e_pr = dPhiSR125(xyzR125_e,xyzR125_pr);
      dphis_e_pi = dPhiSR125(xyzR125_e,xyzR125_pi);
      detas_e_pi = dPhiSR125(xyzR125_e,xyzR125_pi);
    }
    
    if(vlc.M() > fCutsRD[GetGlobalIndex(0,ptbin)])
    {
      okcand = kFALSE;
    }
    if(cosoa < fCutsRD[GetGlobalIndex(1,ptbin)])
    {
      okcand = kFALSE;
    }
    if(fabs(dphis_e_pr) < fCutsRD[GetGlobalIndex(2,ptbin)] && fabs(detas_e_pr) < fCutsRD[GetGlobalIndex(3,ptbin)])
    {
      okcand = kFALSE;
    }
    if(fabs(dphis_e_pi) < fCutsRD[GetGlobalIndex(2,ptbin)] && fabs(detas_e_pi) < fCutsRD[GetGlobalIndex(3,ptbin)])
    {
      okcand = kFALSE;
    }

    if(!okcand)  return 0;
    returnvalueCuts = 1;
  }
  
  Int_t returnvaluePID=1;
  if(selectionLevel==AliRDHFCuts::kAll ||
     selectionLevel==AliRDHFCuts::kCandidate|| 
     selectionLevel==AliRDHFCuts::kPID) {
  }
  
  Int_t returnvalue = 0;
  if(returnvalueCuts==1 && returnvaluePID==1) returnvalue=1;
  
  return returnvalue;
}
//---------------------------------------------------------------------------
Double_t AliRDHFCutsLctoeleLambdafromAODtracks::CalculatePhotonMass(AliAODTrack* trk1, AliAODTrack *trk2) 
{
  //
  // Calculate invariant mass of two tracks at the points where they are parallel
  //
  Double_t pol = 1.;
  if(fBzkG<0.) pol = -1.;

  Double_t xyz[3];

  Double_t mass = 9999.;

  Double_t Rcurv1 = trk1->Pt()/0.3/(0.1*fabs(fBzkG))*100.;
  Double_t charge1 = (Double_t) trk1->Charge();
  Double_t ux1 = trk1->Px()/trk1->Pt();
  Double_t uy1 = trk1->Py()/trk1->Pt();
  trk1->GetXYZ(xyz);
  Double_t x1 = xyz[0];
  Double_t y1 = xyz[1];
  Double_t xc1 = x1 + charge1*Rcurv1*(uy1)*pol;
  Double_t yc1 = y1 + charge1*Rcurv1*(-1.*ux1)*pol;

  Double_t Rcurv2 = trk2->Pt()/0.3/(0.1*fabs(fBzkG))*100.;
  Double_t charge2 = (Double_t) trk2->Charge();
  Double_t ux2 = trk2->Px()/trk2->Pt();
  Double_t uy2 = trk2->Py()/trk2->Pt();
  trk2->GetXYZ(xyz);
  Double_t x2 = xyz[0];
  Double_t y2 = xyz[1];
  Double_t xc2 = x2 + charge2*Rcurv2*(uy2)*pol;
  Double_t yc2 = y2 + charge2*Rcurv2*(-1.*ux2)*pol;

  if(fabs(sqrt((xc2-xc1)*(xc2-xc1)+(yc2-yc1)*(yc2-yc1))-(Rcurv1+Rcurv2))>5.) return 9999.;//not photon candidate


  Double_t dxc = xc2-xc1;
  Double_t dyc = yc2-yc1;
  Double_t nux = dyc/sqrt(dxc*dxc+dyc*dyc);
  Double_t nuy = -1.*dxc/sqrt(dxc*dxc+dyc*dyc);
  Double_t vertx1 = xc1 + (xc2-xc1)*Rcurv1/sqrt(dxc*dxc+dyc*dyc);
  Double_t verty1 = yc1 + (yc2-yc1)*Rcurv1/sqrt(dxc*dxc+dyc*dyc);
  Double_t vertx2 = xc2 + (xc1-xc2)*Rcurv2/sqrt(dxc*dxc+dyc*dyc);
  Double_t verty2 = yc2 + (yc1-yc2)*Rcurv2/sqrt(dxc*dxc+dyc*dyc);
  Double_t vertxmid = (vertx1+vertx2)/2.;
  Double_t vertymid = (verty1+verty2)/2.;

  if(vertxmid*nux+vertymid*nuy<0.) return 9999.;
  
  Double_t px1_new = trk1->Pt()*nux;
  Double_t py1_new = trk1->Pt()*nuy;
  Double_t pz1_new = trk1->Pz();
  Double_t E1_new = sqrt(px1_new*px1_new+py1_new*py1_new+pz1_new*pz1_new+0.000511);
  Double_t px2_new = trk2->Pt()*nux;
  Double_t py2_new = trk2->Pt()*nuy;
  Double_t pz2_new = trk2->Pz();
  Double_t E2_new = sqrt(px2_new*px2_new+py2_new*py2_new+pz2_new*pz2_new+0.000511);


  mass = sqrt(pow(E1_new+E2_new,2)-pow(px1_new+px2_new,2)-pow(py1_new+py2_new,2)-pow(pz1_new+pz2_new,2));

  //check
//  cout<<mass<<endl;
//  cout<<"Track 1: "<<charge1<<" "<<trk1->Px()<<" "<<trk1->Py()<<" "<<trk1->Pz()<<endl;
//  cout<<"Track 2: "<<charge2<<" "<<trk2->Px()<<" "<<trk2->Py()<<" "<<trk2->Pz()<<endl;
//  cout<<"Circle info 1: "<<Rcurv1<<" "<<xc1<<" "<<yc1<<endl;
//  cout<<"Circle info 2: "<<Rcurv2<<" "<<xc2<<" "<<yc2<<endl;
//  cout<<"New momenta 1: "<<px1_new<<" "<<py1_new<<" "<<pz1_new<<endl;
//  cout<<"New momenta 2: "<<px2_new<<" "<<py2_new<<" "<<pz2_new<<endl;
//  AliExternalTrackParam etp1;
//  etp1.CopyFromVTrack(trk1);
//  if(!etp1.PropagateTo(80.,(Float_t)fBzkG)) return 9999.;
//  etp1.GetXYZ(xyz); 
//  cout<<"At 80  cm 1: "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" -> "<<pow(xyz[0]-xc1,2)+pow(xyz[1]-yc1,2)-Rcurv1*Rcurv1<<endl;
//  if(!etp1.PropagateTo(120.,(Float_t)fBzkG)) return 9999.;
//  etp1.GetXYZ(xyz); 
//  cout<<"At 120 cm 1: "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" -> "<<pow(xyz[0]-xc1,2)+pow(xyz[1]-yc1,2)-Rcurv1*Rcurv1<<endl;
//  AliExternalTrackParam etp2;
//  etp2.CopyFromVTrack(trk2);
//  if(!etp2.PropagateTo(80.,(Float_t)fBzkG)) return 9999.;
//  etp2.GetXYZ(xyz); 
//  cout<<"At 80  cm 2: "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" -> "<<pow(xyz[0]-xc2,2)+pow(xyz[1]-yc2,2)-Rcurv2*Rcurv2<<endl;
//  if(!etp2.PropagateTo(120.,(Float_t)fBzkG)) return 9999.;
//  etp2.GetXYZ(xyz); 
//  cout<<"At 120 cm 2: "<<xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<" -> "<<pow(xyz[0]-xc2,2)+pow(xyz[1]-yc2,2)-Rcurv2*Rcurv2<<endl;
//  cout<<endl;

  return mass;
}

//________________________________________________________________________
Double_t AliRDHFCutsLctoeleLambdafromAODtracks::DeltaPhi(AliAODv0 *v0, AliAODTrack *trk)
{
  //
  // Calculate Delta phi
  //
  Double_t phiv = v0->Phi();
  Double_t phie = trk->Phi();
  Double_t dphi = phiv - phie;
  if(dphi<-M_PI) dphi += 2 *M_PI;
  if(dphi> M_PI) dphi -= 2 *M_PI;
  return dphi;
}
//________________________________________________________________________
Double_t AliRDHFCutsLctoeleLambdafromAODtracks::DeltaEta(AliAODv0 *v0, AliAODTrack *trk)
{
  //
  // Calculate Delta Eta
  //
  Double_t etav = v0->Eta();
  Double_t etae = trk->Eta();
  Double_t deta = etav - etae;
  return deta;
}
