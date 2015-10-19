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
  fProdTrackTPCNclsPIDMin(0),
  fProdTrackTPCNclsRatioMin(0.0),
  fProdUseAODFilterBit(kTRUE),
  fProdV0MassTolLambda(0.01),
  fProdV0MassTolLambdaRough(0.01),
  fProdV0PtMin(0.5),
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
	fSigmaElectronTOFMax(9999.)
{
  //
  // Default Constructor
  //

  const Int_t nvars=2;
  SetNVars(nvars);
  TString varNames[nvars]={
			   "InvMass [GeV/c2]", //0
			   "cos(Opening angle) [cos(rad)]" //1
  };

  Bool_t isUpperCut[nvars]={
			    kTRUE, //0
			    kFALSE //1
  };
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[nvars]={
			kTRUE, //0
			kTRUE //1
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
  fProdTrackTPCNclsPIDMin(source.fProdTrackTPCNclsPIDMin),
  fProdTrackTPCNclsRatioMin(source.fProdTrackTPCNclsRatioMin),
  fProdUseAODFilterBit(source.fProdUseAODFilterBit),
  fProdV0MassTolLambda(source.fProdV0MassTolLambda),
  fProdV0MassTolLambdaRough(source.fProdV0MassTolLambdaRough),
  fProdV0PtMin(source.fProdV0PtMin),
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
	fSigmaElectronTOFMax(source.fSigmaElectronTOFMax)
{
  //
  // Copy constructor
  //
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
  fProdTrackTPCNclsPIDMin = source.fProdTrackTPCNclsPIDMin;
  fProdTrackTPCNclsRatioMin = source.fProdTrackTPCNclsRatioMin;
  fProdUseAODFilterBit = source.fProdUseAODFilterBit;
  fProdV0MassTolLambda = source.fProdV0MassTolLambda;
  fProdV0MassTolLambdaRough = source.fProdV0MassTolLambdaRough;
  fProdV0PtMin = source.fProdV0PtMin;
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
    
    if(InvMassEleLambda > fCutsRD[GetGlobalIndex(0,ptbin)])
      {
	okcand = kFALSE;
      }
    if(cosoa < fCutsRD[GetGlobalIndex(1,ptbin)])
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
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::SingleTrkCuts(AliAODTrack *trk, AliAODVertex *primVert)
{
  //
  // Single Track Cut to be applied before object creation
  //

  if(trk->GetStatus()&AliESDtrack::kITSpureSA) return kFALSE;
  if(!(trk->GetStatus()&AliESDtrack::kITSin)) return kFALSE;
  if(fProdUseAODFilterBit && !trk->TestFilterMask(BIT(4))) return kFALSE;
	Double_t pos[3]; primVert->GetXYZ(pos);
	Double_t cov[6]; primVert->GetCovarianceMatrix(cov);
	const AliESDVertex vESD(pos,cov,100.,100);
	if(fTrackCuts&&!IsDaughterSelected(trk,&vESD,fTrackCuts)) return kFALSE;

	if(trk->GetTPCsignalN()<fProdTrackTPCNclsPIDMin) return kFALSE;
	if(trk->GetTPCNclsF()>0){
		Float_t tpcratio = (Float_t)trk->GetTPCncls()/(Float_t)trk->GetTPCNclsF();
		if(tpcratio<fProdTrackTPCNclsRatioMin) return kFALSE;
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
						return IsSelectedeID(trk);
						break;
					case kNSigmaCustomizedCuts:
						return IsSelectedCustomizedeID(trk);
						break;
					case kNSigmaCustomizedPtDepCuts:
						return IsSelectedCustomizedPtDepeID(trk);
						break;
					case kCombinedCuts:
						return IsSelectedCombinedeID(trk);
						break;
			}

    }

  return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::SingleTrkCutsNoPID(AliAODTrack *trk, AliAODVertex *primVert)
{
  //
  // Single Track Cut to be applied before object creation
  //

  if(trk->GetStatus()&AliESDtrack::kITSpureSA) return kFALSE;
  if(!(trk->GetStatus()&AliESDtrack::kITSin)) return kFALSE;
  if(fProdUseAODFilterBit && !trk->TestFilterMask(BIT(4))) return kFALSE;
	Double_t pos[3]; primVert->GetXYZ(pos);
	Double_t cov[6]; primVert->GetCovarianceMatrix(cov);
	const AliESDVertex vESD(pos,cov,100.,100);
	if(fTrackCuts&&!IsDaughterSelected(trk,&vESD,fTrackCuts)) return kFALSE;

	if(trk->GetTPCsignalN()<fProdTrackTPCNclsPIDMin) return kFALSE;
	if(trk->GetTPCNclsF()>0){
		Float_t tpcratio = (Float_t)trk->GetTPCncls()/(Float_t)trk->GetTPCNclsF();
		if(tpcratio<fProdTrackTPCNclsRatioMin) return kFALSE;
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
Bool_t AliRDHFCutsLctoeleLambdafromAODtracks::IsSelectedCustomizedPtDepeID(AliAODTrack *trk)
{
  //
  // electron ID first shot
  //

	Double_t nSigmaTPCele = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trk,AliPID::kElectron);
	Double_t nSigmaTOFele = fPidHF->GetPidResponse()->NumberOfSigmasTOF(trk,AliPID::kElectron);

	if(nSigmaTOFele<fSigmaElectronTOFMin) return kFALSE;
	if(nSigmaTOFele>fSigmaElectronTOFMax) return kFALSE;

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
	if(!v0) return kFALSE;
	if(!(v0->GetSecondaryVtx())) return kFALSE;

  AliAODTrack *cptrack =  (AliAODTrack*)(v0->GetDaughter(0));
  AliAODTrack *cntrack =  (AliAODTrack*)(v0->GetDaughter(1));
  if(!cptrack || !cntrack) return kFALSE;
	if(cptrack->Charge()<0 && cntrack->Charge()>0){
		//In case sign is wrong
		cptrack =  (AliAODTrack*)(v0->GetDaughter(1));
		cntrack =  (AliAODTrack*)(v0->GetDaughter(0));
	}

  if ( cptrack->Charge() == cntrack->Charge() ) return kFALSE;
  if(!(cptrack->GetStatus() & AliESDtrack::kTPCrefit) ||
     !(cntrack->GetStatus() & AliESDtrack::kTPCrefit)) return kFALSE;
  AliAODVertex *maybeKinkPos = (AliAODVertex*)cptrack->GetProdVertex();
  AliAODVertex *maybeKinkNeg = (AliAODVertex*)cntrack->GetProdVertex();
  if (maybeKinkPos->GetType()==AliAODVertex::kKink || maybeKinkNeg->GetType()==AliAODVertex::kKink) 
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

      Int_t isProton= -9999;
      Int_t isPion= -9999;
			if(isparticle){
				isProton=fPidObjProton->MakeRawPid(cptrack,4); 
				isPion=fPidObjPion->MakeRawPid(cntrack,2); 
			}else{
				isProton=fPidObjProton->MakeRawPid(cntrack,4); 
				isPion=fPidObjPion->MakeRawPid(cptrack,2); 
			}
      if(isProton<1) return kFALSE;
      if(isPion<1) return kFALSE;
    }

	Double_t RapLambda = v0->RapLambda();
	if(RapLambda<fProdV0RapMin || RapLambda>fProdV0RapMax) return kFALSE;

	Double_t EtaLambda = v0->PseudoRapV0();
	if(EtaLambda<fProdV0EtaMin || EtaLambda>fProdV0EtaMax) return kFALSE;

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

