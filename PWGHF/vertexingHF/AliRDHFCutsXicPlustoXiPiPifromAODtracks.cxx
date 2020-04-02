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
// Class for cuts on AOD reconstructed Lc->p+K0s
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
#include "AliRDHFCutsXicPlustoXiPiPifromAODtracks.h"
#include "AliAODRecoCascadeHF3Prong.h"
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
ClassImp(AliRDHFCutsXicPlustoXiPiPifromAODtracks);
/// \endcond


//--------------------------------------------------------------------------
AliRDHFCutsXicPlustoXiPiPifromAODtracks::AliRDHFCutsXicPlustoXiPiPifromAODtracks(const char* name) :
  AliRDHFCuts(name),
  fPIDStrategy(kNSigmaCuts),
  fCombinedPIDThreshold(0.),
  fUseCascadePID(kFALSE),
  fPidObjCascPi(0),
  fPidObjCascPr(0),
  fProdTrackPtMin(0.),
  fProdTrackEtaRange(9999.),
  fProdUseAODFilterBit(kTRUE),
  fProdMassTolLambda(0.010),
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
  fProdCascNTPCCrossedRowsMin(0.0),
  fProdCascNTPCCrossedOverFindableRatioMin(0.0),
  fProdLikeSignDcaMax(2.0),
  fProdRoughMassTol(0.25),
  fProdRoughPtMin(0.0)
{
  //
  // Default Constructor
  //

  const Int_t nvars=14;
  SetNVars(nvars);
  TString varNames[nvars]={"Xic inv. mass [GeV/c2]",                   //  0
			   "Xic Pt [GeV/c]", //1
			   "Xi mass Tolrelance [GeV/c2]", //2
			   "Lambda mass Tolrelance [GeV/c2]", //3
			   "Max DCA pi-pi [cm]", //4
			   "Max DCA pi-casc [cm]",//5
			   "Max d0 pi [cm]",//6
			   "Max d0 Xi [cm]",//7
			   "Min d0 Xi-Bach [cm]",//8
			   "Min d0 Xi-V0 [cm]",//9
			   "Min Xic cosPA ",//10
			   "Min DecayLengthXY ",//11
			   "Min Bachelor pT",//12
			   "Max DecayLengthXY "//13
  };

  Bool_t isUpperCut[nvars]={kTRUE,  //  0
			    kFALSE, //1
			    kTRUE, //2
			    kTRUE, //3
			    kTRUE, //4
			    kTRUE, //5
			    kTRUE, //6
			    kTRUE, //7
			    kFALSE, //8
			    kFALSE, //9
			    kFALSE, //10
			    kFALSE,//11
			    kFALSE, //12
			    kTRUE //13
  };
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[nvars]={kFALSE, //  0
			kFALSE, //1
			kTRUE, //2
			kTRUE, //3
			kTRUE, //4
			kTRUE, //5
			kTRUE, //6
			kTRUE, //7
			kTRUE, //8
			kTRUE, //9
			kTRUE, //10
			kTRUE, //11
			kTRUE, //12
			kTRUE  //13
  };
  SetVarsForOpt(nvars,forOpt);

  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);
}

//--------------------------------------------------------------------------
AliRDHFCutsXicPlustoXiPiPifromAODtracks::AliRDHFCutsXicPlustoXiPiPifromAODtracks(const AliRDHFCutsXicPlustoXiPiPifromAODtracks &source) :
  AliRDHFCuts(source),
  fPIDStrategy(source.fPIDStrategy),
  fCombinedPIDThreshold(source.fCombinedPIDThreshold),
  fUseCascadePID(source.fUseCascadePID),
  fPidObjCascPi(source.fPidObjCascPi),
  fPidObjCascPr(source.fPidObjCascPr),
  fProdTrackPtMin(source.fProdTrackPtMin),
  fProdTrackEtaRange(source.fProdTrackEtaRange),
  fProdUseAODFilterBit(source.fProdUseAODFilterBit),
  fProdMassTolLambda(source.fProdMassTolLambda),
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
  fProdCascNTPCCrossedRowsMin(source.fProdCascNTPCCrossedRowsMin),
  fProdCascNTPCCrossedOverFindableRatioMin(source.fProdCascNTPCCrossedRowsMin),
  fProdLikeSignDcaMax(source.fProdLikeSignDcaMax),
  fProdRoughMassTol(source.fProdRoughMassTol),
  fProdRoughPtMin(source.fProdRoughPtMin)
{
  //
  // Copy constructor
  //
}

//--------------------------------------------------------------------------
AliRDHFCutsXicPlustoXiPiPifromAODtracks &AliRDHFCutsXicPlustoXiPiPifromAODtracks::operator=(const AliRDHFCutsXicPlustoXiPiPifromAODtracks &source)
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
  fProdUseAODFilterBit = source.fProdUseAODFilterBit;
  fProdTrackPtMin = source.fProdTrackPtMin;
  fProdTrackEtaRange = source.fProdTrackEtaRange;
  fProdMassTolLambda = source.fProdMassTolLambda;
  fProdMassTolXi = source.fProdMassTolXi;
  fProdMassRejOmega = source.fProdMassRejOmega;
  fProdRfidMinV0 = source.fProdRfidMinV0;
  fProdRfidMaxV0 = source.fProdRfidMaxV0;
  fProdRfidMinXi = source.fProdRfidMinXi;
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
  fProdCascNTPCCrossedRowsMin = source.fProdCascNTPCCrossedRowsMin;
  fProdCascNTPCCrossedOverFindableRatioMin = source.fProdCascNTPCCrossedRowsMin;
  fProdLikeSignDcaMax = source.fProdLikeSignDcaMax;
  fProdRoughMassTol = source.fProdRoughMassTol;
  fProdRoughPtMin = source.fProdRoughPtMin;

  return *this;
}

//---------------------------------------------------------------------------
AliRDHFCutsXicPlustoXiPiPifromAODtracks::~AliRDHFCutsXicPlustoXiPiPifromAODtracks() {
  //
  //  Default Destructor
  //
}

//---------------------------------------------------------------------------
void AliRDHFCutsXicPlustoXiPiPifromAODtracks::GetCutVarsForOpt(AliAODRecoDecayHF *d, Float_t *vars, Int_t nvars, Int_t *pdgdaughters) {
  //
  // Fills in vars the values of the variables
  //

  if (pdgdaughters[0]==-9999) return; // dummy

  AliAODRecoCascadeHF3Prong* dd=(AliAODRecoCascadeHF3Prong*)d;
  if(!dd){
    AliDebug(2," No AliAODRecoCascadeHF3Prong object found\n");
    return;
  }

  if (nvars!=fnVarsForOpt) {
    AliError("AliRDHFCutsXicPlustoXiPiPifromAODtracks wrong number of variables\n");
    return;
  }

  //Double_t ptD=d->Pt();
  //Int_t ptbin=PtBin(ptD);
  Int_t iter=-1;

  if(fVarsForOpt[0]){
    iter++;
    Double_t mxicPDG =  TDatabasePDG::Instance()->GetParticle(4232)->Mass();
    vars[iter]= TMath::Abs(dd->InvMassPiXiPi()-mxicPDG) ;
  }
  if(fVarsForOpt[1]){
    iter++;
    vars[iter]= dd->Pt();
  }
  if(fVarsForOpt[2]){
    iter++;
    Double_t mxiPDG =  TDatabasePDG::Instance()->GetParticle(3312)->Mass();
    vars[iter]= TMath::Abs(dd->CascMassXi()-mxiPDG);
  }
  if(fVarsForOpt[3]){
    iter++;
    Double_t mLPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
    vars[iter]= TMath::Abs(dd->CascMassLambda()-mLPDG);
  }
  Double_t dca[3];
  dd->GetDCAs(dca);
  if(fVarsForOpt[4]){
    iter++;
    vars[iter]= dca[2];
  }
  if(fVarsForOpt[5]){
    iter++;
    vars[iter]= dca[0];
  }
  if(fVarsForOpt[6]){
    iter++;
    vars[iter]= dd->Getd0Prong(0);
  }
  if(fVarsForOpt[7]){
    iter++;
    vars[iter]= dd->Getd0Prong(1);
  }
  if(fVarsForOpt[8]){
    iter++;
    vars[iter]= dd->CascDcaBachToPrimVertex();
  }
  if(fVarsForOpt[9]){
    iter++;
    vars[iter]= dd->CascDcaV0ToPrimVertex();
  }
  if(fVarsForOpt[10]){
    iter++;
    vars[iter]= dd->XicCosPointingAngle();
  }
  if(fVarsForOpt[11]){
    iter++;
    vars[iter]= dd->DecayLengthXY();
  }
  if(fVarsForOpt[12]){
    iter++;
    vars[iter]= dd->PtProng(0);
  }
  if(fVarsForOpt[13]){
    iter++;
    vars[iter]= dd->DecayLengthXY();
  }
  return;
}
//---------------------------------------------------------------------------
Int_t AliRDHFCutsXicPlustoXiPiPifromAODtracks::IsSelected(TObject* obj,Int_t selectionLevel) 
{
  //
  // Apply selection
  //

  if (!fCutsRD) {
    AliFatal("Cut matrice not inizialized. Exit...");
    return 0;
  }

  AliAODRecoCascadeHF3Prong* d=(AliAODRecoCascadeHF3Prong*)obj;
  if(!d){
    AliDebug(2," No AliAODRecoCascadeHF3Prong object found\n");
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

    Double_t mLPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
    Double_t mxiPDG =  TDatabasePDG::Instance()->GetParticle(3312)->Mass();
    Double_t mxicPDG =  TDatabasePDG::Instance()->GetParticle(4232)->Mass();
    if(TMath::Abs(d->InvMassPiXiPi()-mxicPDG) > fCutsRD[GetGlobalIndex(0,ptbin)])
      {
	okcand = kFALSE;
      }
    if(d->Pt()< fCutsRD[GetGlobalIndex(1,ptbin)])
      {
	okcand = kFALSE;
      }
    if(TMath::Abs(d->CascMassXi()-mxiPDG) > fCutsRD[GetGlobalIndex(2,ptbin)])
      {
	okcand = kFALSE;
      }
    if((TMath::Abs(d->CascMassLambda()-mLPDG) > fCutsRD[GetGlobalIndex(3,ptbin)]) &&(TMath::Abs(d->CascMassAntiLambda()-mLPDG) > fCutsRD[GetGlobalIndex(3,ptbin)]) )
      {
	okcand = kFALSE;
      }
    Double_t dca[3];
    d->GetDCAs(dca);
    if(TMath::Abs(dca[2]) > fCutsRD[GetGlobalIndex(4,ptbin)])
      {
	okcand = kFALSE;
      }
    if((TMath::Abs(dca[0]) > fCutsRD[GetGlobalIndex(5,ptbin)]) && (TMath::Abs(dca[1]) > fCutsRD[GetGlobalIndex(5,ptbin)]) ) 
      {
	okcand = kFALSE;
      }
    if((TMath::Abs(d->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(6,ptbin)]) && (TMath::Abs(d->Getd0Prong(2)) > fCutsRD[GetGlobalIndex(6,ptbin)]) ) 
      {
	okcand = kFALSE;
      }
    if((TMath::Abs(d->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(7,ptbin)])) 
      {
	okcand = kFALSE;
      }
    if((TMath::Abs(d->CascDcaBachToPrimVertex()) < fCutsRD[GetGlobalIndex(8,ptbin)])) 
      {
	okcand = kFALSE;
      }
    if((TMath::Abs(d->CascDcaV0ToPrimVertex()) < fCutsRD[GetGlobalIndex(9,ptbin)])) 
      {
	okcand = kFALSE;
      }
    if( d->XicCosPointingAngle() < fCutsRD[GetGlobalIndex(10,ptbin)]) 
      {
	okcand = kFALSE;
      }
    if( d->DecayLengthXY() < fCutsRD[GetGlobalIndex(11,ptbin)]) 
      {
	okcand = kFALSE;
      }
    if( d->DecayLengthXY() > fCutsRD[GetGlobalIndex(13,ptbin)]) 
      {
	okcand = kFALSE;
      }
    if( d->PtProng(0) < fCutsRD[GetGlobalIndex(12,ptbin)] || d->PtProng(2) < fCutsRD[GetGlobalIndex(12,ptbin)]  ) 
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

    switch(fPIDStrategy){
    case kNSigmaCuts:
      returnvaluePID = IsSelectedPID(d);
      break;
    case kCombinedCuts:
      returnvaluePID = IsSelectedCombinedPID(d);
      break;
    }
  }

  Int_t returnvalue = 0;
  if(returnvalueCuts==1 && returnvaluePID==1) returnvalue=1;

  return returnvalue;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsXicPlustoXiPiPifromAODtracks::IsSelectedPID(AliAODRecoDecayHF* obj) 
{
  //
  //  PID selection
  //

  if(!fUsePID || !obj) return 1;

  AliAODRecoCascadeHF3Prong* dd=(AliAODRecoCascadeHF3Prong*)obj;
  AliAODTrack *part1 = dd->GetBachelor1();
  AliAODTrack *part2 = dd->GetBachelor2();

  Int_t returnvalue=1;

  if(fPidHF->GetPidResponse()==0x0){
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
    AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
    fPidHF->SetPidResponse(pidResp);
  }

  Int_t isPion1=fPidHF->MakeRawPid(part1,2); 
  Int_t isPion2=fPidHF->MakeRawPid(part2,2); 

  if(isPion1<1) returnvalue = 0;
  if(isPion2<1) returnvalue = 0;

  return returnvalue;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsXicPlustoXiPiPifromAODtracks::IsSelectedCombinedPID(AliAODRecoDecayHF* obj) {
  //
  //  Combined PID selection
  //
    
  if(!fUsePID || !obj) {return 1;}

  AliAODRecoCascadeHF3Prong* dd=(AliAODRecoCascadeHF3Prong*)obj;
  AliAODTrack *part1 = dd->GetBachelor1();
  AliAODTrack *part2 = dd->GetBachelor2();
  if(!part1||!part2) return 0;

  Int_t returnvalue=1;
  Double_t probPion1 = GetPionProbabilityTPCTOF(part1);
  Double_t probPion2 = GetPionProbabilityTPCTOF(part2);
  if(probPion1<fCombinedPIDThreshold) returnvalue = 0;
  if(probPion2<fCombinedPIDThreshold) returnvalue = 0;
  return returnvalue;
}

//________________________________________________________________________
Double_t AliRDHFCutsXicPlustoXiPiPifromAODtracks::GetPionProbabilityTPCTOF(AliAODTrack* trk) 
{
  //
  //  Get Pion Probablility
  //
  //fPidHF->GetPidCombined()->SetDefaultTPCPriors();
	if(!fPidHF->GetUseCombined()) return -9999.;
  fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
  Double_t prob1[AliPID::kSPECIES];
  UInt_t detUsed1 = fPidHF->GetPidCombined()->ComputeProbabilities(trk, fPidHF->GetPidResponse(), prob1);
  if (detUsed1 != (UInt_t)fPidHF->GetPidCombined()->GetDetectorMask() ) 
    { 
      fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC);
      detUsed1 = fPidHF->GetPidCombined()->ComputeProbabilities(trk, fPidHF->GetPidResponse(), prob1);
    }
  return prob1[AliPID::kPion];
}

//________________________________________________________________________
Bool_t AliRDHFCutsXicPlustoXiPiPifromAODtracks::SingleTrkCuts(AliAODTrack *trk)
{
  //
  //  Single Track Cut
  //
  
  if(trk->GetStatus()&AliESDtrack::kITSpureSA) return kFALSE;
  if(!(trk->GetStatus()&AliESDtrack::kITSin)) return kFALSE;

  if(fProdUseAODFilterBit && !trk->TestFilterMask(BIT(4))) return kFALSE;
  //	if(!fAnalCuts->IsDaughterSelected(trk,fV1,esdTrackCuts)) return kFALSE;
  if(fabs(trk->Eta())>fProdTrackEtaRange) return kFALSE;
  if(trk->Pt()<fProdTrackPtMin) return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRDHFCutsXicPlustoXiPiPifromAODtracks::SingleCascadeCuts(AliAODcascade *casc,Double_t *primvert)
{
  //
  //  Single Cascade Cut
  //
	
  if(!casc) return kFALSE;

  AliAODTrack *ptrack = (AliAODTrack*) (casc->GetDaughter(0));
  AliAODTrack *ntrack = (AliAODTrack*) (casc->GetDaughter(1));
  AliAODTrack *btrack = (AliAODTrack*) (casc->GetDecayVertexXi()->GetDaughter(0));
  
  if(!ptrack||!ntrack||!btrack) return kFALSE;

  //obsolete
  //if(ptrack->GetTPCClusterInfo(2,1)<fProdCascNTPCClustersMin) return kFALSE;
  //if(ntrack->GetTPCClusterInfo(2,1)<fProdCascNTPCClustersMin) return kFALSE;
  //if(btrack->GetTPCClusterInfo(2,1)<fProdCascNTPCClustersMin) return kFALSE;

   if ( ptrack->GetTPCNCrossedRows()< fProdCascNTPCCrossedRowsMin ) return kFALSE;
  if ( ntrack->GetTPCNCrossedRows()< fProdCascNTPCCrossedRowsMin ) return kFALSE;
  if ( btrack->GetTPCNCrossedRows()< fProdCascNTPCCrossedRowsMin ) return kFALSE;
  if ( ptrack->GetTPCNclsF()==0 || ntrack->GetTPCNclsF()==0 || btrack->GetTPCNclsF()==0) return kFALSE;
  if ( static_cast<Double_t>(ptrack->GetTPCNCrossedRows())/static_cast<Double_t>(ptrack->GetTPCNclsF()) <= fProdCascNTPCCrossedOverFindableRatioMin ) return kFALSE;
  if ( static_cast<Double_t>(ntrack->GetTPCNCrossedRows())/static_cast<Double_t>(ntrack->GetTPCNclsF()) <= fProdCascNTPCCrossedOverFindableRatioMin ) return kFALSE;
  if ( static_cast<Double_t>(btrack->GetTPCNCrossedRows())/static_cast<Double_t>(btrack->GetTPCNclsF()) <= fProdCascNTPCCrossedOverFindableRatioMin ) return kFALSE;
  

  Double_t mLPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  Double_t mxiPDG =  TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  Double_t momegaPDG =  TDatabasePDG::Instance()->GetParticle(3334)->Mass();
  
  Double_t massLambda = casc->MassLambda();
  Double_t massAntiLambda = casc->MassAntiLambda();
  if(TMath::Abs(massLambda-mLPDG)>fProdMassTolLambda && TMath::Abs(massAntiLambda-mLPDG)>fProdMassTolLambda) 
    return kFALSE;

  Bool_t isparticle = kTRUE;
  if(TMath::Abs(massAntiLambda-mLPDG)<fProdMassTolLambda) isparticle = kFALSE;
  
  Double_t massXi = casc->MassXi();
  if(TMath::Abs(massXi-mxiPDG)>fProdMassTolXi)
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
      Int_t isProton=fPidObjCascPr->MakeRawPid(ptrack,4); 
      Int_t isPion1 =fPidObjCascPi->MakeRawPid(btrack,2); 
      Int_t isPion2 =fPidObjCascPi->MakeRawPid(ntrack,2); 
      if(isProton<1) return kFALSE;
      if(isPion1<1) return kFALSE;
      if(isPion2<1) return kFALSE;
    }else{
      Int_t isProton=fPidObjCascPr->MakeRawPid(ntrack,4); 
      Int_t isPion1 =fPidObjCascPi->MakeRawPid(btrack,2); 
      Int_t isPion2 =fPidObjCascPi->MakeRawPid(ptrack,2); 
      if(isProton<1) return kFALSE;
      if(isPion1<1) return kFALSE;
      if(isPion2<1) return kFALSE;
    }
  }

  
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsXicPlustoXiPiPifromAODtracks::SingleCascadeCutsRef(AliAODcascade *casc,Double_t *primvert)
{
  //
  //  Single Cascade Cut (without Xi mass selection)
  //  Kinematical cut is applied to compare with cascade analysis note
  //
	
  if(!casc) return kFALSE;

  AliAODTrack *ptrack = (AliAODTrack*) (casc->GetDaughter(0));
  AliAODTrack *ntrack = (AliAODTrack*) (casc->GetDaughter(1));
  AliAODTrack *btrack = (AliAODTrack*) (casc->GetDecayVertexXi()->GetDaughter(0));
  
  if(!ptrack||!ntrack||!btrack) return kFALSE;

  //Kinematical cut
  //if(TMath::Abs(ptrack->Eta())>0.8) return kFALSE;
  //if(TMath::Abs(ntrack->Eta())>0.8) return kFALSE;
  //if(TMath::Abs(btrack->Eta())>0.8) return kFALSE;
  //if(casc->RapXi()<-0.5) return kFALSE;
  //if(casc->RapXi()>0.0) return kFALSE;

  if(ptrack->GetTPCClusterInfo(2,1)<fProdCascNTPCClustersMin) return kFALSE;
  if(ntrack->GetTPCClusterInfo(2,1)<fProdCascNTPCClustersMin) return kFALSE;
  if(btrack->GetTPCClusterInfo(2,1)<fProdCascNTPCClustersMin) return kFALSE;

  Double_t mLPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  Double_t momegaPDG =  TDatabasePDG::Instance()->GetParticle(3334)->Mass();
  Double_t mxiPDG =  TDatabasePDG::Instance()->GetParticle(3312)->Mass();

  Double_t massXi = casc->MassXi();
  Double_t massOmega = casc->MassOmega();
  if(TMath::Abs(massOmega-momegaPDG)<fProdMassRejOmega)
    return kFALSE;

  
  Double_t massLambda = casc->MassLambda();
  Double_t massAntiLambda = casc->MassAntiLambda();
  if(TMath::Abs(massLambda-mLPDG)>fProdMassTolLambda && TMath::Abs(massAntiLambda-mLPDG)>fProdMassTolLambda) 
    return kFALSE;

  Bool_t isparticle = kTRUE;
  if(TMath::Abs(massAntiLambda-mLPDG)<fProdMassTolLambda) isparticle = kFALSE;
  
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
      Int_t isProton=fPidObjCascPr->MakeRawPid(ptrack,4); 
      Int_t isPion1 =fPidObjCascPi->MakeRawPid(btrack,2); 
      Int_t isPion2 =fPidObjCascPi->MakeRawPid(ntrack,2); 
      if(isProton<1) return kFALSE;
      if(isPion1<1) return kFALSE;
      if(isPion2<1) return kFALSE;
    }else{
      Int_t isProton=fPidObjCascPr->MakeRawPid(ntrack,4); 
      Int_t isPion1 =fPidObjCascPi->MakeRawPid(btrack,2); 
      Int_t isPion2 =fPidObjCascPi->MakeRawPid(ptrack,2); 
      if(isProton<1) return kFALSE;
      if(isPion1<1) return kFALSE;
      if(isPion2<1) return kFALSE;
    }
  }
  
  return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsXicPlustoXiPiPifromAODtracks::SelectWithRoughCuts(AliAODcascade *casc, AliAODTrack *part1, AliAODTrack *part2)
{
  //
  //  Select With Rough mass and pT cut before object creation
  //

  //Constants
  Double_t mpiPDG =  TDatabasePDG::Instance()->GetParticle(211)->Mass();
  Double_t mxiPDG =  TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  Double_t mxicPDG =  TDatabasePDG::Instance()->GetParticle(4232)->Mass();

  Double_t pxpi1_init = part1->Px();
  Double_t pypi1_init = part1->Py();
  Double_t pzpi1_init = part1->Pz();
  Double_t Epi1_init = sqrt(pxpi1_init*pxpi1_init+pypi1_init*pypi1_init+pzpi1_init*pzpi1_init+mpiPDG*mpiPDG);
  Double_t pxpi2_init = part2->Px();
  Double_t pypi2_init = part2->Py();
  Double_t pzpi2_init = part2->Pz();
  Double_t Epi2_init = sqrt(pxpi2_init*pxpi2_init+pypi2_init*pypi2_init+pzpi2_init*pzpi2_init+mpiPDG*mpiPDG);
  Double_t pxcasc_init = casc->MomXiX();
  Double_t pycasc_init = casc->MomXiY();
  Double_t pzcasc_init = casc->MomXiZ();
  Double_t Ecasc_init = sqrt(pxcasc_init*pxcasc_init+pycasc_init*pycasc_init+pzcasc_init*pzcasc_init+mxiPDG*mxiPDG);
  Double_t pxxic_init = pxpi1_init+pxpi2_init+pxcasc_init;
  Double_t pyxic_init = pypi1_init+pypi2_init+pycasc_init;
  Double_t pzxic_init = pzpi1_init+pzpi2_init+pzcasc_init;
  Double_t Exic_init = Epi1_init+Epi2_init+Ecasc_init;
  Double_t xicmass_init = sqrt(Exic_init*Exic_init-pxxic_init*pxxic_init-pyxic_init*pyxic_init-pzxic_init*pzxic_init);
  
  if(xicmass_init<mxicPDG-fProdRoughMassTol || xicmass_init>mxicPDG+fProdRoughMassTol) return kFALSE;
  if(sqrt(pxxic_init*pxxic_init+pyxic_init*pyxic_init)<fProdRoughPtMin) return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRDHFCutsXicPlustoXiPiPifromAODtracks::IsInFiducialAcceptance(Double_t pt, Double_t y) const
{
  if(pt > 5.) {
    // applying cut for pt > 5 GeV
    AliDebug(2,Form("pt of Xic = %f (> 5), cutting at |y| < 0.8",pt));
    if (TMath::Abs(y) > 0.8) return kFALSE;
  } else {
    // appliying smooth cut for pt < 5 GeV
    Double_t maxFiducialY = -0.2/15*pt*pt+1.9/15*pt+0.5;
    Double_t minFiducialY = 0.2/15*pt*pt-1.9/15*pt-0.5;
    AliDebug(2,Form("pt of Xic = %f (< 5), cutting  according to the fiducial zone [%f, %f]\n",pt,minFiducialY,maxFiducialY));
    if (y < minFiducialY || y > maxFiducialY) return kFALSE;
  }
  //
  return kTRUE;
}
