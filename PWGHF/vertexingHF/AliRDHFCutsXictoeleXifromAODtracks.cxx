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
#include <TVector3.h>

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
  fUseV0Topology(0),
  fBzkG(0),
  fProdTrackTPCNclsPIDMin(0),
  fProdTrackTPCNclsRatioMin(0.0),
  fProdUseAODFilterBit(kTRUE),
  fProdAODFilterBit(4),
  fProdRejectTrackWithShared(kFALSE),
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
  fProdCascCutMinNCrossedRowsTPC(0.0),
  fProdCascratioCrossedRowsOverFindableClusterTPC(0.0),
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
	fSigmaElectronTOFMax(9999.),
	fConversionMassMax(-1.),
	fEleXiMassMax(2.5)
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
			kTRUE //3
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
  fUseV0Topology(source.fUseV0Topology),
  fBzkG(source.fBzkG),
  fProdTrackTPCNclsPIDMin(source.fProdTrackTPCNclsPIDMin),
  fProdTrackTPCNclsRatioMin(source.fProdTrackTPCNclsRatioMin),
  fProdUseAODFilterBit(source.fProdUseAODFilterBit),
  fProdAODFilterBit(source.fProdAODFilterBit),
  fProdRejectTrackWithShared(source.fProdRejectTrackWithShared),
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
 fProdCascCutMinNCrossedRowsTPC(source.fProdCascCutMinNCrossedRowsTPC),
  fProdCascratioCrossedRowsOverFindableClusterTPC(source.fProdCascratioCrossedRowsOverFindableClusterTPC),     	
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
	fSigmaElectronTOFMax(source.fSigmaElectronTOFMax),
	fConversionMassMax(source.fConversionMassMax),
	fEleXiMassMax(source.fEleXiMassMax)
{
  //
  // Copy constructor
  //
  for(Int_t i=0;i<3;i++){
    fPrimVert[i] = source.fPrimVert[i];
  }
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
  fUseV0Topology = source.fUseV0Topology;
  fBzkG = source.fBzkG;
  fProdTrackTPCNclsPIDMin = source.fProdTrackTPCNclsPIDMin;
  fProdTrackTPCNclsRatioMin = source.fProdTrackTPCNclsRatioMin;
  fProdUseAODFilterBit = source.fProdUseAODFilterBit;
  fProdAODFilterBit = source.fProdAODFilterBit;
  fProdRejectTrackWithShared = source.fProdRejectTrackWithShared;
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
  fProdCascCutMinNCrossedRowsTPC = source.fProdCascCutMinNCrossedRowsTPC;
  fProdCascratioCrossedRowsOverFindableClusterTPC = source.fProdCascratioCrossedRowsOverFindableClusterTPC; 
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
	fConversionMassMax = source.fConversionMassMax;
	fEleXiMassMax = source.fEleXiMassMax;
  
  for(Int_t i=0;i<3;i++){
    fPrimVert[i] = source.fPrimVert[i];
  }
  
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

    Double_t dphis_e_pr, detas_e_pr, dphis_e_pi, detas_e_pi, dphis_e_bach, detas_e_bach;
    dphis_e_pr = 9999.;
    detas_e_pr = 9999.;
    dphis_e_pi = 9999.;
    detas_e_pi = 9999.;
    dphis_e_bach = 9999.;
    detas_e_bach = 9999.;
    if(fCutsRD[GetGlobalIndex(2,ptbin)]>0 || fCutsRD[GetGlobalIndex(3,ptbin)]>0){
      AliAODTrack *trke = (AliAODTrack*)d->GetDaughter(0);
      AliAODcascade *casc = (AliAODcascade*)d->GetDaughter(1);
      if(trke && casc){
        AliAODTrack *cprtrk = 0;
        AliAODTrack *cpitrk = 0;
        AliAODTrack *cbtrk = 0;
        if(casc->ChargeXi()<0){
          cprtrk = (AliAODTrack*)casc->GetDaughter(0);
          cpitrk = (AliAODTrack*)casc->GetDaughter(1);
          cbtrk = (AliAODTrack*)casc->GetDecayVertexXi()->GetDaughter(0);
        }else{
          cprtrk = (AliAODTrack*)casc->GetDaughter(1);
          cpitrk = (AliAODTrack*)casc->GetDaughter(0);
          cbtrk = (AliAODTrack*)casc->GetDecayVertexXi()->GetDaughter(0);
        }
        if(cprtrk && cpitrk && cbtrk)
          GetdPhiSdEtaSR125(trke,cprtrk,cpitrk,cbtrk,fBzkG,fPrimVert, dphis_e_pr,detas_e_pr,dphis_e_pi,detas_e_pi,dphis_e_bach,detas_e_bach);
      }
    }
    
    if(InvMassEleXi > fCutsRD[GetGlobalIndex(0,ptbin)])
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
    if(fabs(dphis_e_bach) < fCutsRD[GetGlobalIndex(2,ptbin)] && fabs(detas_e_bach) < fCutsRD[GetGlobalIndex(3,ptbin)])
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
Bool_t AliRDHFCutsXictoeleXifromAODtracks::SingleTrkCuts(AliAODTrack *trk, AliAODTrack *trkpid, AliAODVertex *primVert)
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
Bool_t AliRDHFCutsXictoeleXifromAODtracks::SingleTrkCutsNoPID(AliAODTrack *trk, AliAODTrack *trkpid, AliAODVertex *primVert)
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
Bool_t AliRDHFCutsXictoeleXifromAODtracks::IsSelectedCustomizedPtDepeID(AliAODTrack *trk, AliAODTrack *trkpid)
{
  //
  // electron ID pt dependent
  //

	Double_t nSigmaTPCele = fPidHF->GetPidResponse()->NumberOfSigmasTPC(trkpid,AliPID::kElectron);
	Double_t nSigmaTOFele = fPidHF->GetPidResponse()->NumberOfSigmasTOF(trkpid,AliPID::kElectron);

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

	//=============== add the CrossedRows and Ratio cut on the daughters =============

	float nCrossedRowsTPCp = ptrack -> GetTPCCrossedRows();  
    float nCrossedRowsTPCn = ntrack -> GetTPCCrossedRows();  
    float nCrossedRowsTPCb = btrack -> GetTPCCrossedRows();
    float ratioCrossedRowsOverFindableClusterTPCp = 1.0;
    float ratioCrossedRowsOverFindableClusterTPCn = 1.0;
    float ratioCrossedRowsOverFindableClusterTPCb = 1.0;
    if (ptrack->GetTPCNclsF()>0){
        ratioCrossedRowsOverFindableClusterTPCp = nCrossedRowsTPCp/ptrack->GetTPCNclsF();
    }
 
    if (ntrack->GetTPCNclsF()>0){
        ratioCrossedRowsOverFindableClusterTPCn = nCrossedRowsTPCn/ntrack->GetTPCNclsF();
    }
 
    if (btrack->GetTPCNclsF()>0){
        ratioCrossedRowsOverFindableClusterTPCb = nCrossedRowsTPCb/btrack->GetTPCNclsF();
    }
     
    if(nCrossedRowsTPCp < fProdCascCutMinNCrossedRowsTPC) return kFALSE;
    if(nCrossedRowsTPCn < fProdCascCutMinNCrossedRowsTPC) return kFALSE;
    if(nCrossedRowsTPCb < fProdCascCutMinNCrossedRowsTPC) return kFALSE;
   
    if(ratioCrossedRowsOverFindableClusterTPCp < fProdCascratioCrossedRowsOverFindableClusterTPC) return kFALSE;
    if(ratioCrossedRowsOverFindableClusterTPCn < fProdCascratioCrossedRowsOverFindableClusterTPC) return kFALSE;
    if(ratioCrossedRowsOverFindableClusterTPCb < fProdCascratioCrossedRowsOverFindableClusterTPC) return kFALSE;

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

  if(fProdRejectTrackWithShared){
    const TBits sharedMap1 = ptrack->GetTPCSharedMap();
    const TBits sharedMap2 = ntrack->GetTPCSharedMap();
    const TBits sharedMap3 = btrack->GetTPCSharedMap();
    if((sharedMap1.CountBits() >= 1) || (sharedMap2.CountBits() >= 1) ||
        (sharedMap3.CountBits() >= 1))
    {
      return kFALSE;
    }
  }

  if(fUseV0Topology>0){
    Double_t dphiprlam = 0.;
    if(isparticle){
      TVector3 v3v0pr(casc->MomPosX(),casc->MomPosY(),casc->MomPosZ());
      TVector3 v3lam(casc->Px(),casc->Py(),casc->Pz());
      dphiprlam = v3v0pr.DeltaPhi(v3lam);
    }else{
      TVector3 v3v0pr(casc->MomNegX(),casc->MomNegY(),casc->MomNegZ());
      TVector3 v3lam(casc->Px(),casc->Py(),casc->Pz());
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

//________________________________________________________________________
Bool_t AliRDHFCutsXictoeleXifromAODtracks::TagConversions(AliAODTrack *etrk, Int_t *id2index, AliAODEvent *evt, Int_t ntrk, Double_t &minmass)
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
		if(mass<minmass) minmass = mass;
	}

	if(minmass<fConversionMassMax) isconv = kTRUE;

  return isconv;
}
//________________________________________________________________________
Bool_t AliRDHFCutsXictoeleXifromAODtracks::TagConversionsSameSign(AliAODTrack *etrk, Int_t *id2index, AliAODEvent *evt, Int_t ntrk, Double_t &minmass)
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
void AliRDHFCutsXictoeleXifromAODtracks::SetSftPosR125(AliAODTrack *track,Double_t bfield,Double_t priVtx[3],Double_t *XSftR125)
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
Double_t AliRDHFCutsXictoeleXifromAODtracks::dEtaSR125(Double_t *postrack1,Double_t *postrack2)
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
Double_t AliRDHFCutsXictoeleXifromAODtracks::dPhiSR125(Double_t *postrack1,Double_t *postrack2)
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
Double_t AliRDHFCutsXictoeleXifromAODtracks::GetdPhiSdEtaSR125(AliAODTrack *tracke, AliAODTrack *trackp,AliAODTrack *trackb,
    AliAODTrack *trackn, Double_t bfield, Double_t priVtx[3], Double_t &dPhiS_ep, Double_t &dEtaS_ep, 
    Double_t &dPhiS_en, Double_t &dEtaS_en, Double_t &dPhiS_eb, Double_t &dEtaS_eb)
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
  Double_t XSftR125_b[3];
  SetSftPosR125(trackb,bfield,priVtx, XSftR125_b);
  dPhiS_ep = dPhiSR125(XSftR125_e,XSftR125_p);
  dEtaS_ep = dEtaSR125(XSftR125_e,XSftR125_p);
  dPhiS_en = dPhiSR125(XSftR125_e,XSftR125_n);
  dEtaS_en = dEtaSR125(XSftR125_e,XSftR125_n);
  dPhiS_eb = dPhiSR125(XSftR125_e,XSftR125_b);
  dEtaS_eb = dEtaSR125(XSftR125_e,XSftR125_b);
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsXictoeleXifromAODtracks::IsSelected(TLorentzVector* vtrk, TLorentzVector *vcasc, Double_t *cutvars, Int_t selectionLevel) 
{
  //
  // Apply selection on mixed event tracks
  //

  if (!fCutsRD) {
    AliFatal("Cut matrix not inizialized. Exit...");
    return 0;
  }

  Double_t ptD=sqrt(pow(vtrk->Px()+vcasc->Px(),2)+pow(vtrk->Px()+vcasc->Py(),2));
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
    
    Double_t mlcPDG =  TDatabasePDG::Instance()->GetParticle(4132)->Mass();
    Double_t melePDG =  TDatabasePDG::Instance()->GetParticle(11)->Mass();
    Double_t mxiPDG =  TDatabasePDG::Instance()->GetParticle(3312)->Mass();
		Double_t v0px = vcasc->Px();
		Double_t v0py = vcasc->Py();
		Double_t v0pz = vcasc->Pz();
		Double_t epx = vtrk->Px();
		Double_t epy = vtrk->Py();
		Double_t epz = vtrk->Pz();
		Double_t cosoa = (v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz);

    TLorentzVector vele, vxi,vlc;
    vele.SetXYZM(epx,epy,epz,melePDG);
    vxi.SetXYZM(v0px,v0py,v0pz,mxiPDG);
    vlc = vele + vxi;

    Double_t dphis_e_pr = 9999.;
    Double_t detas_e_pr = 9999.;
    Double_t dphis_e_pi = 9999.;
    Double_t detas_e_pi = 9999.;
    Double_t dphis_e_bach = 9999.;
    Double_t detas_e_bach = 9999.;
    if(fCutsRD[GetGlobalIndex(2,ptbin)]>0 || fCutsRD[GetGlobalIndex(3,ptbin)]>0){
      Double_t xyzR125_e[3], xyzR125_pr[3], xyzR125_pi[3], xyzR125_bach[3];
      xyzR125_e[0] = cutvars[0];
      xyzR125_e[1] = cutvars[1];
      xyzR125_e[2] = cutvars[2];
      xyzR125_pr[0] = cutvars[3];
      xyzR125_pr[1] = cutvars[4];
      xyzR125_pr[2] = cutvars[5];
      xyzR125_pi[0] = cutvars[6];
      xyzR125_pi[1] = cutvars[7];
      xyzR125_pi[2] = cutvars[8];
      xyzR125_bach[0] = cutvars[9];
      xyzR125_bach[1] = cutvars[10];
      xyzR125_bach[2] = cutvars[11];
      dphis_e_pr = dPhiSR125(xyzR125_e,xyzR125_pr);
      detas_e_pr = dPhiSR125(xyzR125_e,xyzR125_pr);
      dphis_e_pi = dPhiSR125(xyzR125_e,xyzR125_pi);
      detas_e_pi = dPhiSR125(xyzR125_e,xyzR125_pi);
      dphis_e_bach = dPhiSR125(xyzR125_e,xyzR125_bach);
      detas_e_bach = dPhiSR125(xyzR125_e,xyzR125_bach);
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
    if(fabs(dphis_e_bach) < fCutsRD[GetGlobalIndex(2,ptbin)] && fabs(detas_e_bach) < fCutsRD[GetGlobalIndex(3,ptbin)])
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
//________________________________________________________________________
Double_t AliRDHFCutsXictoeleXifromAODtracks::DeltaPhi(AliAODcascade *casc, AliAODTrack *trk)
{
  //
  // Calculate Delta phi
  //
  Double_t phiv = atan2(casc->MomXiY(),casc->MomXiX());
  if(phiv<-M_PI) phiv += 2 *M_PI;
  if(phiv> M_PI) phiv -= 2 *M_PI;
  Double_t phie = trk->Phi();
  if(phie<-M_PI) phie += 2 *M_PI;
  if(phie> M_PI) phie -= 2 *M_PI;
  Double_t dphi = phiv - phie;
  if(dphi<-M_PI) dphi += 2 *M_PI;
  if(dphi> M_PI) dphi -= 2 *M_PI;
  return dphi;
}
//________________________________________________________________________
Double_t AliRDHFCutsXictoeleXifromAODtracks::DeltaEta(AliAODcascade *casc, AliAODTrack *trk)
{
  //
  // Calculate Delta Eta
  //
  Double_t ptotxi = TMath::Sqrt(pow(casc->MomXiX(),2)+pow(casc->MomXiY(),2)+pow(casc->MomXiZ(),2));
	Double_t etav = 0.5*TMath::Log((ptotxi+casc->MomXiZ())/(ptotxi-casc->MomXiZ()));
  Double_t etae = trk->Eta();
  Double_t deta = etav - etae;
  return deta;
}

//________________________________________________________________________
Double_t AliRDHFCutsXictoeleXifromAODtracks::CosOpeningAngle(AliAODcascade *casc, AliAODTrack *trk)
{
  //
  // Calculate Cosine of opening angle
  //

  Double_t xipx = casc->MomXiX();
  Double_t xipy = casc->MomXiY();
  Double_t xipz = casc->MomXiZ();
  Double_t epx = trk->Px();
  Double_t epy = trk->Py();
  Double_t epz = trk->Pz();
  Double_t cosoa = (xipx*epx+xipy*epy+xipz*epz)/sqrt(xipx*xipx+xipy*xipy+xipz*xipz)/sqrt(epx*epx+epy*epy+epz*epz);
  return cosoa;
}

