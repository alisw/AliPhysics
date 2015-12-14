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
// Class for cuts on AOD reconstructed Xic->e+Xi
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
#include "AliRDHFCutsXictoeleXifromAODtracks.h"
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
ClassImp(AliRDHFCutsXictoeleXifromAODtracks);
/// \endcond


//--------------------------------------------------------------------------
AliRDHFCutsXictoeleXifromAODtracks::AliRDHFCutsXictoeleXifromAODtracks(const char* name) :
AliRDHFCuts(name),
  fPIDStrategy(kNSigmaCuts),
  fCombinedPIDThreshold(0.),
  fUseCascadePID(kFALSE),
  fPidObjCascPi(0),
  fPidObjCascPr(0),
  fProdTrackTPCNclsPIDMin(0),
  fProdTrackTPCNclsRatioMin(0.0),
  fProdUseAODFilterBit(kTRUE),
  fProdMassTolLambda(0.010),
  fProdMassTolXiRough(0.020),
  fProdMassTolXi(0.008),
  fProdMassRejOmega(0.008),
  fProdRfidMinV0(0.6),
  fProdRfidMaxV0(100.0),
  fProdRfidMinXi(0.6),
  fProdRfidMaxXi(100.0),
  fProdCascProperDecayLengthMax(100.0),
  fProdDcaXiDaughtersMax(100.),
  fProdDcaV0DaughtersMax(100.),
  fProdDcaBachToPrimVertexMin(0.),
  fProdDcaV0ToPrimVertexMin(0.),
  fProdDcaV0PrToPrimVertexMin(0.),
  fProdDcaV0PiToPrimVertexMin(0.),
  fProdXiCosineOfPoiningAngleMin(-1.),
  fProdV0CosineOfPoiningAngleXiMin(-1.),
  fProdCascNTPCClustersMin(0.0),
  fProdCascEtaMin(-9999.),
  fProdCascEtaMax(9999.),
  fProdCascRapMin(-9999.),
  fProdCascRapMax(9999.),
  fProdCascPtMin(0.),
  fProdCascPtMax(9999.),
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
			kTRUE, //1
  };
  SetVarsForOpt(nvars,forOpt);

  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);
}
//--------------------------------------------------------------------------
AliRDHFCutsXictoeleXifromAODtracks::AliRDHFCutsXictoeleXifromAODtracks(const AliRDHFCutsXictoeleXifromAODtracks &source) :
  AliRDHFCuts(source),
  fPIDStrategy(source.fPIDStrategy),
  fCombinedPIDThreshold(source.fCombinedPIDThreshold),
  fUseCascadePID(source.fUseCascadePID),
  fPidObjCascPi(source.fPidObjCascPi),
  fPidObjCascPr(source.fPidObjCascPr),
  fProdTrackTPCNclsPIDMin(source.fProdTrackTPCNclsPIDMin),
  fProdTrackTPCNclsRatioMin(source.fProdTrackTPCNclsRatioMin),
  fProdUseAODFilterBit(source.fProdUseAODFilterBit),
  fProdMassTolLambda(source.fProdMassTolLambda),
  fProdMassTolXiRough(source.fProdMassTolXiRough),
  fProdMassTolXi(source.fProdMassTolXi),
  fProdMassRejOmega(source.fProdMassRejOmega),
  fProdRfidMinV0(source.fProdRfidMinV0),
  fProdRfidMaxV0(source.fProdRfidMaxV0),
  fProdRfidMinXi(source.fProdRfidMinXi),
  fProdRfidMaxXi(source.fProdRfidMaxXi),
  fProdCascProperDecayLengthMax(source.fProdCascProperDecayLengthMax),
  fProdDcaXiDaughtersMax(source.fProdDcaXiDaughtersMax),
  fProdDcaV0DaughtersMax(source.fProdDcaV0DaughtersMax),
  fProdDcaBachToPrimVertexMin(source.fProdDcaBachToPrimVertexMin),
  fProdDcaV0ToPrimVertexMin(source.fProdDcaV0ToPrimVertexMin),
  fProdDcaV0PrToPrimVertexMin(source.fProdDcaV0PrToPrimVertexMin),
  fProdDcaV0PiToPrimVertexMin(source.fProdDcaV0PiToPrimVertexMin),
  fProdXiCosineOfPoiningAngleMin(source.fProdXiCosineOfPoiningAngleMin),
  fProdV0CosineOfPoiningAngleXiMin(source.fProdV0CosineOfPoiningAngleXiMin),
  fProdCascNTPCClustersMin(source.fProdCascNTPCClustersMin),
  fProdCascEtaMin(source.fProdCascEtaMin),
  fProdCascEtaMax(source.fProdCascEtaMax),
  fProdCascRapMin(source.fProdCascRapMin),
  fProdCascRapMax(source.fProdCascRapMax),
  fProdCascPtMin(source.fProdCascPtMin),
  fProdCascPtMax(source.fProdCascPtMax),
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
AliRDHFCutsXictoeleXifromAODtracks &AliRDHFCutsXictoeleXifromAODtracks::operator=(const AliRDHFCutsXictoeleXifromAODtracks &source)
{
  //
  // assignment operator
  //

  if (this != &source) {
    AliRDHFCuts::operator=(source);
  }

  fPIDStrategy = source.fPIDStrategy;
  fCombinedPIDThreshold = source.fCombinedPIDThreshold;
  fUseCascadePID = source.fUseCascadePID;
  fPidObjCascPi = source.fPidObjCascPi;
  fPidObjCascPr = source.fPidObjCascPr;
  fProdTrackTPCNclsPIDMin = source.fProdTrackTPCNclsPIDMin;
  fProdTrackTPCNclsRatioMin = source.fProdTrackTPCNclsRatioMin;
  fProdUseAODFilterBit = source.fProdUseAODFilterBit;
  fProdMassTolLambda = source.fProdMassTolLambda;
  fProdMassTolXiRough = source.fProdMassTolXiRough;
  fProdMassTolXi = source.fProdMassTolXi;
  fProdMassRejOmega = source.fProdMassRejOmega;
  fProdRfidMinV0 = source.fProdRfidMinV0;
  fProdRfidMaxV0 = source.fProdRfidMaxV0;
  fProdRfidMinXi = source.fProdRfidMinXi;
  fProdRfidMaxXi = source.fProdRfidMaxXi;
  fProdCascProperDecayLengthMax = source.fProdCascProperDecayLengthMax;
  fProdDcaXiDaughtersMax = source.fProdDcaXiDaughtersMax;
  fProdDcaV0DaughtersMax = source.fProdDcaV0DaughtersMax;
  fProdDcaBachToPrimVertexMin = source.fProdDcaBachToPrimVertexMin;
  fProdDcaV0ToPrimVertexMin = source.fProdDcaV0ToPrimVertexMin;
  fProdDcaV0PrToPrimVertexMin = source.fProdDcaV0PrToPrimVertexMin;
  fProdDcaV0PiToPrimVertexMin = source.fProdDcaV0PiToPrimVertexMin;
  fProdXiCosineOfPoiningAngleMin = source.fProdXiCosineOfPoiningAngleMin;
  fProdV0CosineOfPoiningAngleXiMin = source.fProdV0CosineOfPoiningAngleXiMin;
  fProdCascNTPCClustersMin = source.fProdCascNTPCClustersMin;
  fProdCascEtaMin = source.fProdCascEtaMin;
  fProdCascEtaMax = source.fProdCascEtaMax;
  fProdCascRapMin = source.fProdCascRapMin;
  fProdCascRapMax = source.fProdCascRapMax;
  fProdCascPtMin = source.fProdCascPtMin;
  fProdCascPtMax = source.fProdCascPtMax;
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
AliRDHFCutsXictoeleXifromAODtracks::~AliRDHFCutsXictoeleXifromAODtracks() {
  //
  //  Default Destructor
  //
  
}

//---------------------------------------------------------------------------
void AliRDHFCutsXictoeleXifromAODtracks::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
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
		UInt_t pdgdg[2]={11,3312};
    vars[iter]= dd->InvMass(2,pdgdg);
  }
  if(fVarsForOpt[1]){
    iter++;
		Double_t xipx = dd->PxProng(1);
		Double_t xipy = dd->PyProng(1);
		Double_t xipz = dd->PzProng(1);
		Double_t epx = dd->PxProng(0);
		Double_t epy = dd->PyProng(0);
		Double_t epz = dd->PzProng(0);
    vars[iter]=  (xipx*epx+xipy*epy+xipz*epz)/sqrt(xipx*xipx+xipy*xipy+xipz*xipz)/sqrt(epx*epx+epy*epy+epz*epz);
  }
  
  return;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsXictoeleXifromAODtracks::IsSelected(TObject* obj, Int_t selectionLevel) {
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
    
    Double_t mxiPDG =  TDatabasePDG::Instance()->GetParticle(3312)->Mass();
		UInt_t pdgdg[2]={11,3312};
		Double_t InvMassEleXi = d->InvMass(2,pdgdg);
		Double_t xipx = d->PxProng(1);
		Double_t xipy = d->PyProng(1);
		Double_t xipz = d->PzProng(1);
		Double_t epx = d->PxProng(0);
		Double_t epy = d->PyProng(0);
		Double_t epz = d->PzProng(0);
		Double_t cosoa = (xipx*epx+xipy*epy+xipz*epz)/sqrt(xipx*xipx+xipy*xipy+xipz*xipz)/sqrt(epx*epx+epy*epy+epz*epz);
    
    if(InvMassEleXi > fCutsRD[GetGlobalIndex(0,ptbin)])
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
		//not used. applied at single level
  }
  
  Int_t returnvalue = 0;
  if(returnvalueCuts==1 && returnvaluePID==1) returnvalue=1;
  
  return returnvalue;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsXictoeleXifromAODtracks::IsSelectedPID(AliAODRecoDecayHF* obj) 
{
  //
  // IsSelectedPID
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
Int_t AliRDHFCutsXictoeleXifromAODtracks::IsSelectedCombinedPID(AliAODRecoDecayHF* obj) {
  //
  // IsSelectedCombinedPID
  //
    
  if(!fUsePID || !obj) {return 1;}

  AliAODRecoCascadeHF* dd=(AliAODRecoCascadeHF*)obj;
  AliAODTrack *part = dd->GetBachelor();
  if(!part) return 0;

  Int_t returnvalue=1;

  return returnvalue;
}

//________________________________________________________________________
Bool_t AliRDHFCutsXictoeleXifromAODtracks::SingleTrkCuts(AliAODTrack *trk, AliAODVertex *primVert)
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
Bool_t AliRDHFCutsXictoeleXifromAODtracks::SingleTrkCutsNoPID(AliAODTrack *trk, AliAODVertex *primVert)
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
Bool_t AliRDHFCutsXictoeleXifromAODtracks::IsSelectedeID(AliAODTrack *trk)
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
Bool_t AliRDHFCutsXictoeleXifromAODtracks::IsSelectedCustomizedeID(AliAODTrack *trk)
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
	if(fabs(fSigmaElectronTOFMin)<10.&&fabs(fSigmaElectronTOFMax)<10.){
		if(nSigmaTOFele<fSigmaElectronTOFMin) return kFALSE;
		if(nSigmaTOFele>fSigmaElectronTOFMax) return kFALSE;
	}

	return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsXictoeleXifromAODtracks::IsSelectedCustomizedPtDepeID(AliAODTrack *trk)
{
  //
  // electron ID pt dependent
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
Bool_t AliRDHFCutsXictoeleXifromAODtracks::IsSelectedCombinedeID(AliAODTrack *trk)
{
  //
  // electron ID Basyian not implemented
  //

	return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsXictoeleXifromAODtracks::SingleCascadeCuts(AliAODcascade *casc,Double_t *primvert)
{
  //
  //  Single Cascade Cut
  //
	
  if(!casc) return kFALSE;
	if(!(casc->GetSecondaryVtx())) return kFALSE;
	if(!(casc->GetDecayVertexXi())) return kFALSE;

  AliAODTrack *ptrack = (AliAODTrack*) (casc->GetDaughter(0));
  AliAODTrack *ntrack = (AliAODTrack*) (casc->GetDaughter(1));
  AliAODTrack *btrack = (AliAODTrack*) (casc->GetDecayVertexXi()->GetDaughter(0));
  
  if(!ptrack||!ntrack||!btrack) return kFALSE;

	if(ptrack->Charge()<0 && ntrack->Charge()>0){
		ptrack = (AliAODTrack*) (casc->GetDaughter(1));
		ntrack = (AliAODTrack*) (casc->GetDaughter(0));
	}

  if(ptrack->GetTPCClusterInfo(2,1)<fProdCascNTPCClustersMin) return kFALSE;
  if(ntrack->GetTPCClusterInfo(2,1)<fProdCascNTPCClustersMin) return kFALSE;
  if(btrack->GetTPCClusterInfo(2,1)<fProdCascNTPCClustersMin) return kFALSE;

  Double_t pxxi = casc->MomXiX();
  Double_t pyxi = casc->MomXiY();
  Double_t ptxi = sqrt(pxxi*pxxi+pyxi*pyxi);
  if(ptxi<fProdCascPtMin || ptxi>fProdCascPtMax) return kFALSE;

  Double_t mLPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  Double_t mxiPDG =  TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  Double_t momegaPDG =  TDatabasePDG::Instance()->GetParticle(3334)->Mass();
  
  Double_t massLambda = casc->MassLambda();
  Double_t massAntiLambda = casc->MassAntiLambda();
  if(TMath::Abs(massLambda-mLPDG)>fProdMassTolLambda && TMath::Abs(massAntiLambda-mLPDG)>fProdMassTolLambda) 
    return kFALSE;

  Bool_t isparticle = kTRUE;
  //if(TMath::Abs(massAntiLambda-mLPDG)<fProdMassTolLambda) isparticle = kFALSE;
	if(casc->ChargeXi()>0) isparticle = kFALSE;
  
  Double_t massXi = casc->MassXi();
  if(TMath::Abs(massXi-mxiPDG)>fProdMassTolXiRough)
    return kFALSE;

  Double_t massOmega = casc->MassOmega();
  if(TMath::Abs(massOmega-momegaPDG)<fProdMassRejOmega)
    return kFALSE;

  
  Double_t lPosXi[3];
  lPosXi[0] = casc->DecayVertexXiX();
  lPosXi[1] = casc->DecayVertexXiY();
  lPosXi[2] = casc->DecayVertexXiZ();
  Double_t decayvertXi = TMath::Sqrt(lPosXi[0]*lPosXi[0]+lPosXi[1]*lPosXi[1]);
  Double_t lPosV0[3];
  lPosV0[0] = casc->DecayVertexV0X();
  lPosV0[1] = casc->DecayVertexV0Y();
  lPosV0[2] = casc->DecayVertexV0Z();
  Double_t decayvertV0 = TMath::Sqrt(lPosV0[0]*lPosV0[0]+lPosV0[1]*lPosV0[1]);
  
  if(decayvertV0<fProdRfidMinV0 || decayvertV0>fProdRfidMaxV0) return kFALSE;
  if(decayvertXi<fProdRfidMinXi || decayvertXi>fProdRfidMaxXi) return kFALSE;

  Double_t ptotxi = TMath::Sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)+pow(casc->MomXiZ(),2));
  Double_t properdl = casc->DecayLengthXi(primvert[0],primvert[1],primvert[2])*mxiPDG/ptotxi;
  if(properdl>fProdCascProperDecayLengthMax) return kFALSE;

	Double_t lDcaXiDaughters = casc->DcaXiDaughters();
	Double_t lDcaV0Daughters = casc->DcaV0Daughters();
	if(lDcaXiDaughters > fProdDcaXiDaughtersMax) return kFALSE;
	if(lDcaV0Daughters > fProdDcaV0DaughtersMax) return kFALSE;

	Double_t lDcaBachToPrimVertex = casc->DcaBachToPrimVertex();
	Double_t lDcaV0ToPrimVertex = casc->DcaV0ToPrimVertex();
	Double_t lDcaPosToPrimVertex = casc->DcaPosToPrimVertex();
	Double_t lDcaNegToPrimVertex = casc->DcaNegToPrimVertex();
	if(lDcaBachToPrimVertex < fProdDcaBachToPrimVertexMin) return kFALSE;
	if(lDcaV0ToPrimVertex < fProdDcaV0ToPrimVertexMin) return kFALSE;
  if(isparticle){
    if(lDcaPosToPrimVertex < fProdDcaV0PrToPrimVertexMin) return kFALSE;
    if(lDcaNegToPrimVertex < fProdDcaV0PiToPrimVertexMin) return kFALSE;
  }else{
    if(lDcaPosToPrimVertex < fProdDcaV0PiToPrimVertexMin) return kFALSE;
    if(lDcaNegToPrimVertex < fProdDcaV0PrToPrimVertexMin) return kFALSE;
  }

	Double_t lXiCosineOfPointingAngle = casc->CosPointingAngleXi(primvert[0],primvert[1],primvert[2]);
	Double_t lV0CosineOfPointingAngleXi = casc->CosPointingAngle(lPosXi);

	if(lXiCosineOfPointingAngle < fProdXiCosineOfPoiningAngleMin) return kFALSE;
	if(lV0CosineOfPointingAngleXi < fProdV0CosineOfPoiningAngleXiMin) return kFALSE;

  if(fUseCascadePID)
  {
    if(fPidObjCascPi->GetPidResponse()==0x0){
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
      AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
      fPidObjCascPi->SetPidResponse(pidResp);
    }
    if(fPidObjCascPr->GetPidResponse()==0x0){
      AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
      AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
      fPidObjCascPr->SetPidResponse(pidResp);
    }
    if(isparticle){
      Int_t isProton = -9999;
      Int_t isPion1 = -9999;
      Int_t isPion2 = -9999;
      Double_t nsigmatpc_proton = fPidObjCascPr->GetSigma(0);
      Double_t nsigmatpc_pion = fPidObjCascPi->GetSigma(0);
      Double_t nSigmaTPCpr = fPidObjCascPr->GetPidResponse()->NumberOfSigmasTPC(ptrack,AliPID::kProton);
      Double_t nSigmaTPCpi1 = fPidObjCascPi->GetPidResponse()->NumberOfSigmasTPC(ntrack,AliPID::kPion);
      Double_t nSigmaTPCpi2 = fPidObjCascPi->GetPidResponse()->NumberOfSigmasTPC(btrack,AliPID::kPion);
      if(fabs(nSigmaTPCpr)<nsigmatpc_proton) isProton = 1;
      if(fabs(nSigmaTPCpi1)<nsigmatpc_pion) isPion1 = 1;
      if(fabs(nSigmaTPCpi2)<nsigmatpc_pion) isPion2 = 1;
      if(isProton<1) return kFALSE;
      if(isPion1<1) return kFALSE;
      if(isPion2<1) return kFALSE;
    }else{
      Int_t isProton = -9999;
      Int_t isPion1 = -9999;
      Int_t isPion2 = -9999;
      Double_t nsigmatpc_proton = fPidObjCascPr->GetSigma(0);
      Double_t nsigmatpc_pion = fPidObjCascPi->GetSigma(0);
      Double_t nSigmaTPCpr = fPidObjCascPr->GetPidResponse()->NumberOfSigmasTPC(ntrack,AliPID::kProton);
      Double_t nSigmaTPCpi1 = fPidObjCascPi->GetPidResponse()->NumberOfSigmasTPC(ptrack,AliPID::kPion);
      Double_t nSigmaTPCpi2 = fPidObjCascPi->GetPidResponse()->NumberOfSigmasTPC(btrack,AliPID::kPion);
      if(fabs(nSigmaTPCpr)<nsigmatpc_proton) isProton = 1;
      if(fabs(nSigmaTPCpi1)<nsigmatpc_pion) isPion1 = 1;
      if(fabs(nSigmaTPCpi2)<nsigmatpc_pion) isPion2 = 1;
      if(isProton<1) return kFALSE;
      if(isPion1<1) return kFALSE;
      if(isPion2<1) return kFALSE;
    }
  }

	Double_t RapXi = casc->RapXi();
	if(RapXi<fProdCascRapMin || RapXi>fProdCascRapMax) return kFALSE;

	Double_t EtaXi = 0.5*TMath::Log((ptotxi+casc->MomXiZ())/(ptotxi-casc->MomXiZ()));
	if(EtaXi<fProdCascEtaMin || EtaXi>fProdCascEtaMax) return kFALSE;
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRDHFCutsXictoeleXifromAODtracks::SelectWithRoughCuts(AliAODcascade *casc, AliAODTrack *part)
{
  //
  // Mass and pT Cut to be applied before object creation
	// Not used now
  //
	if(!casc) return kFALSE;
	if(!part) return kFALSE;

	return kTRUE;
}

//________________________________________________________________________
Bool_t AliRDHFCutsXictoeleXifromAODtracks::IsPeakRegion(AliAODcascade *casc)
{
  Double_t mxiPDG =  TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  Double_t massXi = casc->MassXi();
  if(TMath::Abs(massXi-mxiPDG)<fProdMassTolXi)
		return kTRUE;
	return kFALSE;
}

//________________________________________________________________________
Bool_t AliRDHFCutsXictoeleXifromAODtracks::IsPeakRegion(TLorentzVector *casc)
{
  Double_t mxiPDG =  TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  Double_t massXi = casc->M();
  if(TMath::Abs(massXi-mxiPDG)<fProdMassTolXi)
		return kTRUE;
	return kFALSE;
}

//________________________________________________________________________
Bool_t AliRDHFCutsXictoeleXifromAODtracks::IsSideBand(AliAODcascade *casc)
{
  Double_t mxiPDG =  TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  Double_t massXi = casc->MassXi();
	Bool_t issideband = kFALSE;
  if((massXi-mxiPDG)>fProdMassTolXiRough-fProdMassTolXi) issideband = kTRUE;
  if((massXi-mxiPDG)<-fProdMassTolXiRough+fProdMassTolXi) issideband = kTRUE;
	return issideband;
}

//________________________________________________________________________
Bool_t AliRDHFCutsXictoeleXifromAODtracks::IsSideBand(TLorentzVector *casc)
{
  Double_t mxiPDG =  TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  Double_t massXi = casc->M();
	Bool_t issideband = kFALSE;
  if((massXi-mxiPDG)>fProdMassTolXiRough-fProdMassTolXi) issideband = kTRUE;
  if((massXi-mxiPDG)<-fProdMassTolXiRough+fProdMassTolXi) issideband = kTRUE;
	return issideband;
}

