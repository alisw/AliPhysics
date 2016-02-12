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
#include <TLorentzVector.h>
#include <TVector3.h>

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliRDHFCutsLctopK0sfromAODtracks.h"
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
ClassImp(AliRDHFCutsLctopK0sfromAODtracks);
/// \endcond


//--------------------------------------------------------------------------
AliRDHFCutsLctopK0sfromAODtracks::AliRDHFCutsLctopK0sfromAODtracks(const char* name) :
AliRDHFCuts(name),
  fPIDStrategy(kNSigmaCuts),
  fCombinedPIDThreshold(0.),
  fUseOnTheFlyV0(kFALSE),
  fBzkG(0),
  fProdTrackTPCNclsPIDMin(0),
  fProdTrackTPCNclsRatioMin(0.0),
  fProdUseAODFilterBit(kTRUE),
  fProdAODFilterBit(4),
  fProdRejectTrackWithShared(kFALSE),
  fProdV0MassTolK0s(0.01),
  fProdV0MassRejLambda(0.00),
  fProdV0MassRejPhoton(0.00),
  fProdV0PtMin(0.5),
  fProdV0CosPointingAngleToPrimVtxMin(0.99),
  fProdV0DcaDaughtersMax(1.5),
  fProdV0DaughterEtaRange(0.8),
  fProdV0DaughterPtMin(0.0),
  fProdV0DaughterTPCClusterMin(70),
  fProdV0EtaMin(-9999.),
  fProdV0EtaMax(9999.),
  fProdV0RapMin(-9999.),
  fProdV0RapMax(9999.),
	fProdRoughMassTol(0.25),
	fProdRoughPtMin(0.0),
	fNWeightingProtonBinLimits(0),
	fWeightingProtonBins(0),
	fNWeightingK0sBinLimits(0),
	fWeightingK0sBins(0),
	fNWeightingBins(0),
	fWeight_p0(0),
	fWeight_p1(0),
	fWeight_p2(0),
	fWeight_p3(0),
	fTagV0MassTol(0)
{
  //
  // Default Constructor
  //
  for(Int_t i=0;i<3;i++){
    fPrimVert[i] = 0.;
  }

  const Int_t nvars=9;
  SetNVars(nvars);
  TString varNames[nvars]={"Lc inv. mass [GeV/c2]",                   //  0
			   "Opening angle [rad]",//1
			   "n#sigma_{TPC} max",//2
			   "n#sigma_{TOF} min",//3
			   "Decay Length min [cm]",//4
			   "cos #theta min",//5
			   "cos #theta max",//6
			   "Proton d0 max",//7
			   "V0 d0 max"//8
  };

  Bool_t isUpperCut[nvars]={kTRUE,  //  0
			    kFALSE, //1: opening angle
			    kTRUE, //2: nsigma_tpc max
			    kFALSE, //3: nsigma tof min
			    kFALSE, //4: decay length min
			    kFALSE, //5: cos the min
			    kTRUE, //6: cos the max
			    kTRUE, //7: 
			    kTRUE //8:
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
			kTRUE //8
  };
  SetVarsForOpt(nvars,forOpt);

  Float_t limits[2]={0,999999999.};
  SetPtBins(2,limits);
}
//--------------------------------------------------------------------------
AliRDHFCutsLctopK0sfromAODtracks::AliRDHFCutsLctopK0sfromAODtracks(const AliRDHFCutsLctopK0sfromAODtracks &source) :
  AliRDHFCuts(source),
  fPIDStrategy(source.fPIDStrategy),
  fCombinedPIDThreshold(source.fCombinedPIDThreshold),
  fUseOnTheFlyV0(source.fUseOnTheFlyV0),
  fBzkG(source.fBzkG),
  fProdTrackTPCNclsPIDMin(source.fProdTrackTPCNclsPIDMin),
  fProdTrackTPCNclsRatioMin(source.fProdTrackTPCNclsRatioMin),
  fProdUseAODFilterBit(source.fProdUseAODFilterBit),
  fProdAODFilterBit(source.fProdAODFilterBit),
  fProdRejectTrackWithShared(source.fProdRejectTrackWithShared),
  fProdV0MassTolK0s(source.fProdV0MassTolK0s),
  fProdV0MassRejLambda(source.fProdV0MassRejLambda),
  fProdV0MassRejPhoton(source.fProdV0MassRejPhoton),
  fProdV0PtMin(source.fProdV0PtMin),
  fProdV0CosPointingAngleToPrimVtxMin(source.fProdV0CosPointingAngleToPrimVtxMin),
  fProdV0DcaDaughtersMax(source.fProdV0DcaDaughtersMax),
  fProdV0DaughterEtaRange(source.fProdV0DaughterEtaRange),
  fProdV0DaughterPtMin(source.fProdV0DaughterPtMin),
  fProdV0DaughterTPCClusterMin(source.fProdV0DaughterTPCClusterMin),
  fProdV0EtaMin(source.fProdV0EtaMin),
  fProdV0EtaMax(source.fProdV0EtaMax),
  fProdV0RapMin(source.fProdV0RapMin),
  fProdV0RapMax(source.fProdV0RapMax),
  fProdRoughMassTol(source.fProdRoughMassTol),
  fProdRoughPtMin(source.fProdRoughPtMin),
  fNWeightingProtonBinLimits(source.fNWeightingProtonBinLimits),
  fWeightingProtonBins(NULL),
  fNWeightingK0sBinLimits(source.fNWeightingK0sBinLimits),
  fWeightingK0sBins(NULL),
  fNWeightingBins(source.fNWeightingBins),
  fWeight_p0(NULL),
  fWeight_p1(NULL),
  fWeight_p2(NULL),
  fWeight_p3(NULL),
	fTagV0MassTol(source.fTagV0MassTol)
{
  //
  // Copy constructor
  //
  for(Int_t i=0;i<3;i++){
    fPrimVert[i] = source.fPrimVert[i];
  }
	fWeightingProtonBins = new Double_t[fNWeightingProtonBinLimits];
	fWeightingK0sBins = new Double_t[fNWeightingK0sBinLimits];
	fWeight_p0 = new Double_t[fNWeightingBins];
	fWeight_p1 = new Double_t[fNWeightingBins];
	fWeight_p2 = new Double_t[fNWeightingBins];
	fWeight_p3 = new Double_t[fNWeightingBins];
	for(Int_t i=0;i<fNWeightingProtonBinLimits;i++){
		fWeightingProtonBins[i] = source.fWeightingProtonBins[i];
	}
	for(Int_t i=0;i<fNWeightingK0sBinLimits;i++){
		fWeightingK0sBins[i] = source.fWeightingK0sBins[i];
	}
	for(Int_t i=0;i<fNWeightingBins;i++){
		fWeight_p0[i] = source.fWeight_p0[i];
		fWeight_p1[i] = source.fWeight_p1[i];
		fWeight_p2[i] = source.fWeight_p2[i];
		fWeight_p3[i] = source.fWeight_p3[i];
	}
}
//--------------------------------------------------------------------------
AliRDHFCutsLctopK0sfromAODtracks &AliRDHFCutsLctopK0sfromAODtracks::operator=(const AliRDHFCutsLctopK0sfromAODtracks &source)
{
  //
  // assignment operator
  //

  if (this != &source) {
    AliRDHFCuts::operator=(source);
  }

  fPIDStrategy = source.fPIDStrategy;
  fCombinedPIDThreshold = source.fCombinedPIDThreshold;
  fUseOnTheFlyV0 = source.fUseOnTheFlyV0;
  fBzkG = source.fBzkG;
  fProdTrackTPCNclsPIDMin = source.fProdTrackTPCNclsPIDMin;
  fProdTrackTPCNclsRatioMin = source.fProdTrackTPCNclsRatioMin;
  fProdUseAODFilterBit = source.fProdUseAODFilterBit;
  fProdAODFilterBit = source.fProdAODFilterBit;
  fProdRejectTrackWithShared = source.fProdRejectTrackWithShared;
  fProdV0MassTolK0s = source.fProdV0MassTolK0s;
  fProdV0MassRejLambda = source.fProdV0MassRejLambda;
  fProdV0MassRejPhoton = source.fProdV0MassRejPhoton;
  fProdV0PtMin = source.fProdV0PtMin;
  fProdV0CosPointingAngleToPrimVtxMin = source.fProdV0CosPointingAngleToPrimVtxMin;
  fProdV0DcaDaughtersMax=source.fProdV0DcaDaughtersMax;
  fProdV0DaughterEtaRange=source.fProdV0DaughterEtaRange;
  fProdV0DaughterPtMin=source.fProdV0DaughterPtMin;
  fProdV0DaughterTPCClusterMin=source.fProdV0DaughterTPCClusterMin;
  fProdV0EtaMin = source.fProdV0EtaMin;
  fProdV0EtaMax = source.fProdV0EtaMax;
  fProdV0RapMin = source.fProdV0RapMin;
  fProdV0RapMax = source.fProdV0RapMax;
  fProdRoughMassTol = source.fProdRoughMassTol;
  fProdRoughPtMin = source.fProdRoughPtMin;
  fNWeightingBins = source.fNWeightingBins;
  fNWeightingProtonBinLimits = source.fNWeightingProtonBinLimits;
  fNWeightingK0sBinLimits = source.fNWeightingK0sBinLimits;
  fTagV0MassTol = source.fTagV0MassTol;

  for(Int_t i=0;i<3;i++){
    fPrimVert[i] = source.fPrimVert[i];
  }

	fWeightingProtonBins = new Double_t[fNWeightingProtonBinLimits];
	fWeightingK0sBins = new Double_t[fNWeightingK0sBinLimits];
	fWeight_p0 = new Double_t[fNWeightingBins];
	fWeight_p1 = new Double_t[fNWeightingBins];
	fWeight_p2 = new Double_t[fNWeightingBins];
	fWeight_p3 = new Double_t[fNWeightingBins];
	for(Int_t i=0;i<fNWeightingProtonBinLimits;i++){
		fWeightingProtonBins[i] = source.fWeightingProtonBins[i];
	}
	for(Int_t i=0;i<fNWeightingK0sBinLimits;i++){
		fWeightingK0sBins[i] = source.fWeightingK0sBins[i];
	}
	for(Int_t i=0;i<fNWeightingBins;i++){
		fWeight_p0[i] = source.fWeight_p0[i];
		fWeight_p1[i] = source.fWeight_p1[i];
		fWeight_p2[i] = source.fWeight_p2[i];
		fWeight_p3[i] = source.fWeight_p3[i];
	}
  
  
  return *this;
}

//---------------------------------------------------------------------------
AliRDHFCutsLctopK0sfromAODtracks::~AliRDHFCutsLctopK0sfromAODtracks() {
  //
  //  Default Destructor
  //

	delete [] fWeightingProtonBins;
	delete [] fWeightingK0sBins;
	delete [] fWeight_p0;
	delete [] fWeight_p1;
	delete [] fWeight_p2;
	delete [] fWeight_p3;
  
}

//---------------------------------------------------------------------------
void AliRDHFCutsLctopK0sfromAODtracks::GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
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
    Double_t mlcPDG =  TDatabasePDG::Instance()->GetParticle(4122)->Mass();
    vars[iter]= TMath::Abs(dd->InvMassLctoK0sP()-mlcPDG) ;
  }
  if(fVarsForOpt[1]){
    iter++;
    vars[iter]= dd->Pt();
  }
  
  return;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctopK0sfromAODtracks::IsSelected(TObject* obj, Int_t selectionLevel) {
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

  Double_t pt=d->Pt();
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
    Double_t mprPDG =  TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    Double_t mk0PDG =  TDatabasePDG::Instance()->GetParticle(310)->Mass();
		Double_t v0px = d->PxProng(1);
		Double_t v0py = d->PyProng(1);
		Double_t v0pz = d->PzProng(1);
		Double_t epx = d->PxProng(0);
		Double_t epy = d->PyProng(0);
		Double_t epz = d->PzProng(0);
		Double_t cosoa = (v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz);

    TLorentzVector vpr, vk0s,vlc;
    vpr.SetXYZM(epx,epy,epz,mprPDG);
    vk0s.SetXYZM(v0px,v0py,v0pz,mk0PDG);
    vlc = vpr + vk0s;
    TVector3 vboost = vlc.BoostVector();
    vpr.Boost(-vboost);
    Double_t bachcosthe = cos(vpr.Angle(vlc.Vect()));
    
    if(TMath::Abs(d->InvMassLctoK0sP()-mlcPDG) > fCutsRD[GetGlobalIndex(0,ptbin)])
      {
	okcand = kFALSE;
      }
    if(cosoa < fCutsRD[GetGlobalIndex(1,ptbin)])
      {
	okcand = kFALSE;
      }
    if(CalculateLcCosPAXY(d)*d->DecayLengthXY() < fCutsRD[GetGlobalIndex(4,ptbin)])
      {
	okcand = kFALSE;
      }
    if(bachcosthe < fCutsRD[GetGlobalIndex(5,ptbin)])
      {
	okcand = kFALSE;
      }
    if(bachcosthe > fCutsRD[GetGlobalIndex(6,ptbin)])
      {
	okcand = kFALSE;
      }
    if(fabs(d->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(7,ptbin)])
      {
	okcand = kFALSE;
      }
    if(fabs(d->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(8,ptbin)])
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

    AliAODTrack *part = (AliAODTrack*)d->GetSecondaryVtx()->GetDaughter(0);

    Double_t nSigmaTPCpr = fPidHF->GetPidResponse()->NumberOfSigmasTPC(part,AliPID::kProton);
    if(nSigmaTPCpr>fCutsRD[GetGlobalIndex(2,ptbin)]){
      returnvaluePID = -1;
    }

    Double_t nSigmaTOFpr = fPidHF->GetPidResponse()->NumberOfSigmasTOF(part,AliPID::kProton);
    if(nSigmaTOFpr<fCutsRD[GetGlobalIndex(3,ptbin)]){
      returnvaluePID = -1;
    }

  }
  
  Int_t returnvalue = 0;
  if(returnvalueCuts==1 && returnvaluePID==1) returnvalue=1;
  
  return returnvalue;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctopK0sfromAODtracks::IsSelected(TLorentzVector* vtrk, TLorentzVector *vv0, Double_t *cutvars, Int_t selectionLevel) {
  //
  // Apply selection on mixed event tracks
  //

  if (!fCutsRD) {
    AliFatal("Cut matrix not inizialized. Exit...");
    return 0;
  }

  Double_t ptD=cutvars[1];
  if(ptD<fMinPtCand) return 0;
  if(ptD>fMaxPtCand) return 0;

  Double_t pt=cutvars[1];
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
    Double_t mprPDG =  TDatabasePDG::Instance()->GetParticle(2212)->Mass();
    Double_t mk0PDG =  TDatabasePDG::Instance()->GetParticle(310)->Mass();
		Double_t v0px = vv0->Px();
		Double_t v0py = vv0->Py();
		Double_t v0pz = vv0->Pz();
		Double_t epx = vtrk->Px();
		Double_t epy = vtrk->Py();
		Double_t epz = vtrk->Pz();
		Double_t cosoa = (v0px*epx+v0py*epy+v0pz*epz)/sqrt(v0px*v0px+v0py*v0py+v0pz*v0pz)/sqrt(epx*epx+epy*epy+epz*epz);

    TLorentzVector vpr, vk0s,vlc;
    vpr.SetXYZM(epx,epy,epz,mprPDG);
    vk0s.SetXYZM(v0px,v0py,v0pz,mk0PDG);
    vlc = vpr + vk0s;
    TVector3 vboost = vlc.BoostVector();
    vpr.Boost(-vboost);
    Double_t bachcosthe = cos(vpr.Angle(vlc.Vect()));
    
    if(TMath::Abs(cutvars[0]-mlcPDG) > fCutsRD[GetGlobalIndex(0,ptbin)])
      {
	okcand = kFALSE;
      }
    if(cosoa < fCutsRD[GetGlobalIndex(1,ptbin)])
      {
	okcand = kFALSE;
      }
    if(cutvars[4] < fCutsRD[GetGlobalIndex(4,ptbin)])
      {
	okcand = kFALSE;
      }
    if(bachcosthe < fCutsRD[GetGlobalIndex(5,ptbin)])
      {
	okcand = kFALSE;
      }
    if(bachcosthe > fCutsRD[GetGlobalIndex(6,ptbin)])
      {
	okcand = kFALSE;
      }
    if(cutvars[5] > fCutsRD[GetGlobalIndex(7,ptbin)])
      {
	okcand = kFALSE;
      }
    if(cutvars[6] > fCutsRD[GetGlobalIndex(8,ptbin)])
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

    Double_t nSigmaTPCpr = cutvars[2];
    if(fabs(nSigmaTPCpr)>fCutsRD[GetGlobalIndex(2,ptbin)]){
      returnvaluePID = -1;
    }

    Double_t nSigmaTOFpr = cutvars[3];
    if(fabs(nSigmaTOFpr)<fCutsRD[GetGlobalIndex(3,ptbin)]){
      returnvaluePID = -1;
    }

  }
  
  Int_t returnvalue = 0;
  if(returnvalueCuts==1 && returnvaluePID==1) returnvalue=1;
  
  return returnvalue;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctopK0sfromAODtracks::IsSelectedPID(AliAODRecoDecayHF* obj) 
{
  //
  // IsSelectedPID ( not used)
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

  Int_t isProton=fPidHF->MakeRawPid(part,4); 
  if(isProton<1) returnvalue = 0;
  
  return returnvalue;
}

//---------------------------------------------------------------------------
Bool_t AliRDHFCutsLctopK0sfromAODtracks::IsSelectedProtonID(AliAODTrack* part) 
{
  //
  // IsSelectedProtonID
  //

  if(fPidHF->GetPidResponse()==0x0){
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
    AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
    fPidHF->SetPidResponse(pidResp);
  }

  //Int_t isProton=fPidHF->MakeRawPid(part,4); 
  //if(isProton<1) return kFALSE;
  Double_t nsigmatpc_proton = fPidHF->GetSigma(0);
  Double_t nsigmatof_proton = fPidHF->GetSigma(3);
  Double_t nSigmaTPCpr = fPidHF->GetPidResponse()->NumberOfSigmasTPC(part,AliPID::kProton);
  Double_t nSigmaTOFpr = fPidHF->GetPidResponse()->NumberOfSigmasTOF(part,AliPID::kProton);

	if(fabs(nSigmaTPCpr)<nsigmatpc_proton && fabs(nSigmaTOFpr)<nsigmatof_proton){
		return kTRUE;
	}
  
  return kFALSE;
}
//---------------------------------------------------------------------------
Bool_t AliRDHFCutsLctopK0sfromAODtracks::IsSelectedKaonID(AliAODTrack* part) 
{
  //
  // IsSelectedKaonID
  //

  if(fPidHF->GetPidResponse()==0x0){
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
    AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
    fPidHF->SetPidResponse(pidResp);
  }

  //Int_t isKaon=fPidHF->MakeRawPid(part,3); 
  //if(isKaon<1) return kFALSE;
  Double_t nsigmatpc_kaon = fPidHF->GetSigma(0);
  Double_t nsigmatof_kaon = fPidHF->GetSigma(3);
  Double_t nSigmaTPCka = fPidHF->GetPidResponse()->NumberOfSigmasTPC(part,AliPID::kKaon);
  Double_t nSigmaTOFka = fPidHF->GetPidResponse()->NumberOfSigmasTOF(part,AliPID::kKaon);

	if(fabs(nSigmaTPCka)<nsigmatpc_kaon && fabs(nSigmaTOFka)<nsigmatof_kaon){
		return kTRUE;
	}
  
  return kFALSE;
}

//---------------------------------------------------------------------------
Int_t AliRDHFCutsLctopK0sfromAODtracks::IsSelectedCombinedPID(AliAODRecoDecayHF* obj) {
  //
  // IsSelectedCombinedPID
  //
    
  if(!fUsePID || !obj) {return 1;}

  AliAODRecoCascadeHF* dd=(AliAODRecoCascadeHF*)obj;
  AliAODTrack *part = dd->GetBachelor();
  if(!part) return 0;

  Int_t returnvalue=1;
  Double_t probProton = GetProtonProbabilityTPCTOF(part);
  if(probProton<fCombinedPIDThreshold) returnvalue = 0;

  return returnvalue;
}

//________________________________________________________________________
Double_t AliRDHFCutsLctopK0sfromAODtracks::GetProtonProbabilityTPCTOF(AliAODTrack* trk) 
{
  //
  // Get Proton probability
  //
	if(!fPidHF->GetUseCombined()) return -9999.;
  fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);
  Double_t prob1[AliPID::kSPECIES];
  UInt_t detUsed1 = fPidHF->GetPidCombined()->ComputeProbabilities(trk, fPidHF->GetPidResponse(), prob1);
  if (detUsed1 != (UInt_t)fPidHF->GetPidCombined()->GetDetectorMask() ) 
    { 
      fPidHF->GetPidCombined()->SetDetectorMask(AliPIDResponse::kDetTPC);
      detUsed1 = fPidHF->GetPidCombined()->ComputeProbabilities(trk, fPidHF->GetPidResponse(), prob1);
    }
  return prob1[AliPID::kProton];
}

//________________________________________________________________________
Bool_t AliRDHFCutsLctopK0sfromAODtracks::SingleTrkCuts(AliAODTrack *trk, AliAODTrack *trkpid, AliAODVertex *primVert)
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
					return IsSelectedProtonID(trkpid);
					break;
			}
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
Bool_t AliRDHFCutsLctopK0sfromAODtracks::SingleKaonCuts(AliAODTrack *trk, AliAODVertex *primVert)
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
					return IsSelectedKaonID(trk);
					break;
			}
    }

  return kTRUE;
}
//________________________________________________________________________
Double_t AliRDHFCutsLctopK0sfromAODtracks::CalculateLcCosPAXY(AliAODRecoDecayHF *d) 
{
  AliAODRecoCascadeHF* lcobj=(AliAODRecoCascadeHF*)d;
  if(!lcobj){
    AliError("No AliAODRecoCascadeHF object found\n");
    return -9999.;
  }
  Double_t dvertx = lcobj->GetSecondaryVtx()->GetX()-lcobj->GetOwnPrimaryVtx()->GetX();
  Double_t dverty = lcobj->GetSecondaryVtx()->GetY()-lcobj->GetOwnPrimaryVtx()->GetY();
  Double_t px = lcobj->Px();
  Double_t py = lcobj->Py();
  Double_t inner = px*dvertx + py*dverty;
  if(inner<0.) return  -1.;
  return 1.;
}
//________________________________________________________________________
Bool_t AliRDHFCutsLctopK0sfromAODtracks::SingleV0Cuts(AliAODv0 *v0, AliAODVertex *primVert)
{
  //
  // Single V0 Cut to be applied before object creation
  //

  Bool_t onFlyV0 = v0->GetOnFlyStatus(); // on-the-flight V0s
  if ( onFlyV0 && !fUseOnTheFlyV0 ) return kFALSE;

  AliAODTrack *cptrack =  (AliAODTrack*)(v0->GetDaughter(0));
  AliAODTrack *cntrack =  (AliAODTrack*)(v0->GetDaughter(1));
  if(!cptrack || !cntrack) return kFALSE;
  if ( cptrack->Charge() == cntrack->Charge() ) return kFALSE;
  if(!(cptrack->GetStatus() & AliESDtrack::kTPCrefit) ||
     !(cntrack->GetStatus() & AliESDtrack::kTPCrefit)) return kFALSE;
  AliAODVertex *maybeKinkPos = (AliAODVertex*)cptrack->GetProdVertex();
  AliAODVertex *maybeKinkNeg = (AliAODVertex*)cntrack->GetProdVertex();
  if (maybeKinkPos->GetType()==AliAODVertex::kKink || maybeKinkNeg->GetType()==AliAODVertex::kKink) 
    return kFALSE;

  if ( ( ( cptrack->GetTPCClusterInfo(2,1) ) < (Float_t)fProdV0DaughterTPCClusterMin ) || 
       ( ( cntrack->GetTPCClusterInfo(2,1) ) < (Float_t)fProdV0DaughterTPCClusterMin) ) return kFALSE;

  Double_t massK0S = v0->MassK0Short();
  Double_t mk0sPDG   = TDatabasePDG::Instance()->GetParticle(310)->Mass();
  if(fabs(massK0S-mk0sPDG)>fProdV0MassTolK0s) return kFALSE;

  Double_t massLambda = v0->MassLambda();
  Double_t massAntiLambda = v0->MassAntiLambda();
  Double_t mlamPDG   = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  if((fabs(massAntiLambda-mlamPDG)<fProdV0MassRejLambda) || (fabs(massLambda-mlamPDG)<fProdV0MassRejLambda)) return kFALSE;

  Double_t pxe1 = v0->MomPosX();
  Double_t pye1 = v0->MomPosY();
  Double_t pze1 = v0->MomPosZ();
  Double_t Ee1 = sqrt(pxe1*pxe1+pye1*pye1+pze1*pze1+0.000511*0.000511);
  Double_t pxe2 = v0->MomNegX();
  Double_t pye2 = v0->MomNegY();
  Double_t pze2 = v0->MomNegZ();
  Double_t Ee2 = sqrt(pxe2*pxe2+pye2*pye2+pze2*pze2+0.000511*0.000511);
	Double_t mphoton = sqrt(pow(Ee1+Ee2,2)-pow(pxe1+pxe2,2)-pow(pye1+pye2,2)-pow(pze1+pze2,2));
	if(mphoton<fProdV0MassRejPhoton) return kFALSE;


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

	Double_t RapK0s = v0->RapK0Short();
	if(RapK0s<fProdV0RapMin || RapK0s>fProdV0RapMax) return kFALSE;

	Double_t EtaK0s = v0->PseudoRapV0();
	if(EtaK0s<fProdV0EtaMin || EtaK0s>fProdV0EtaMax) return kFALSE;


  if(fProdRejectTrackWithShared){
    const TBits sharedMap1 = cptrack->GetTPCSharedMap();
    const TBits sharedMap2 = cntrack->GetTPCSharedMap();
    if((sharedMap1.CountBits() >= 1) || (sharedMap2.CountBits() >= 1))
    {
      return kFALSE;
    }
  }


  return kTRUE;
}

//________________________________________________________________________
Bool_t AliRDHFCutsLctopK0sfromAODtracks::SelectWithRoughCuts(AliAODv0 *v0, AliAODTrack *part)
{
  //
  // Mass and pT Cut to be applied before object creation
  //

  Double_t mprPDG =  TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  Double_t mLcPDG =  TDatabasePDG::Instance()->GetParticle(4122)->Mass();

  Double_t pxpr_init = part->Px();
  Double_t pypr_init = part->Py();
  Double_t pzpr_init = part->Pz();
  Double_t Epr_init = sqrt(pxpr_init*pxpr_init+pypr_init*pypr_init+pzpr_init*pzpr_init+mprPDG*mprPDG);

  Double_t pxv0_init = v0->Px();
  Double_t pyv0_init = v0->Py();
  Double_t pzv0_init = v0->Pz();
  Double_t massv0_init = v0->MassK0Short();
  Double_t Ev0_init = sqrt(pxv0_init*pxv0_init+pyv0_init*pyv0_init+pzv0_init*pzv0_init+massv0_init*massv0_init);

  Double_t pxlc_init = pxpr_init+pxv0_init;
  Double_t pylc_init = pypr_init+pyv0_init;
  Double_t pzlc_init = pzpr_init+pzv0_init;
  Double_t Elc_init = Epr_init+Ev0_init;
  Double_t lcmass_init = sqrt(Elc_init*Elc_init-pxlc_init*pxlc_init-pylc_init*pylc_init-pzlc_init*pzlc_init);

  if(lcmass_init<mLcPDG-fProdRoughMassTol || lcmass_init>mLcPDG+fProdRoughMassTol) return kFALSE;
  if(sqrt(pxlc_init*pxlc_init+pylc_init*pylc_init)<fProdRoughPtMin) return kFALSE;

  return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsLctopK0sfromAODtracks::SelectWithRoughCutsWS(AliAODTrack *vka, AliAODTrack *part)
{
  //
  // Mass and pT Cut to be applied before object creation
  //

  Double_t mprPDG =  TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  Double_t mLcPDG =  TDatabasePDG::Instance()->GetParticle(4122)->Mass();

  Double_t pxpr_init = part->Px();
  Double_t pypr_init = part->Py();
  Double_t pzpr_init = part->Pz();
  Double_t Epr_init = sqrt(pxpr_init*pxpr_init+pypr_init*pypr_init+pzpr_init*pzpr_init+mprPDG*mprPDG);

  Double_t pxv0_init = vka->Px();
  Double_t pyv0_init = vka->Py();
  Double_t pzv0_init = vka->Pz();
  Double_t massv0_init = 0.493677;
  Double_t Ev0_init = sqrt(pxv0_init*pxv0_init+pyv0_init*pyv0_init+pzv0_init*pzv0_init+massv0_init*massv0_init);

  Double_t pxlc_init = pxpr_init+pxv0_init;
  Double_t pylc_init = pypr_init+pyv0_init;
  Double_t pzlc_init = pzpr_init+pzv0_init;
  Double_t Elc_init = Epr_init+Ev0_init;
  Double_t lcmass_init = sqrt(Elc_init*Elc_init-pxlc_init*pxlc_init-pylc_init*pylc_init-pzlc_init*pzlc_init);

  if(lcmass_init<mLcPDG-fProdRoughMassTol || lcmass_init>mLcPDG+fProdRoughMassTol) return kFALSE;
  if(sqrt(pxlc_init*pxlc_init+pylc_init*pylc_init)<fProdRoughPtMin) return kFALSE;

  return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsLctopK0sfromAODtracks::SelectWithRoughCuts(TLorentzVector *v0, TLorentzVector *part)
{
  //
  // Mass and pT Cut to be applied before object creation
  //

  Double_t mprPDG =  TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  Double_t mLcPDG =  TDatabasePDG::Instance()->GetParticle(4122)->Mass();

  Double_t pxpr_init = part->Px();
  Double_t pypr_init = part->Py();
  Double_t pzpr_init = part->Pz();
  Double_t Epr_init = sqrt(pxpr_init*pxpr_init+pypr_init*pypr_init+pzpr_init*pzpr_init+mprPDG*mprPDG);

  Double_t pxv0_init = v0->Px();
  Double_t pyv0_init = v0->Py();
  Double_t pzv0_init = v0->Pz();
  Double_t massv0_init = v0->M();
  Double_t Ev0_init = sqrt(pxv0_init*pxv0_init+pyv0_init*pyv0_init+pzv0_init*pzv0_init+massv0_init*massv0_init);

  Double_t pxlc_init = pxpr_init+pxv0_init;
  Double_t pylc_init = pypr_init+pyv0_init;
  Double_t pzlc_init = pzpr_init+pzv0_init;
  Double_t Elc_init = Epr_init+Ev0_init;
  Double_t lcmass_init = sqrt(Elc_init*Elc_init-pxlc_init*pxlc_init-pylc_init*pylc_init-pzlc_init*pzlc_init);

  if(lcmass_init<mLcPDG-fProdRoughMassTol || lcmass_init>mLcPDG+fProdRoughMassTol) return kFALSE;
  if(sqrt(pxlc_init*pxlc_init+pylc_init*pylc_init)<fProdRoughPtMin) return kFALSE;

  return kTRUE;
}
//________________________________________________________________________
Bool_t AliRDHFCutsLctopK0sfromAODtracks::SelectWithRoughCutsWS(TLorentzVector *vka, TLorentzVector *part)
{
  //
  // Mass and pT Cut to be applied before object creation
  //

  Double_t mprPDG =  TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  Double_t mLcPDG =  TDatabasePDG::Instance()->GetParticle(4122)->Mass();

  Double_t pxpr_init = part->Px();
  Double_t pypr_init = part->Py();
  Double_t pzpr_init = part->Pz();
  Double_t Epr_init = sqrt(pxpr_init*pxpr_init+pypr_init*pypr_init+pzpr_init*pzpr_init+mprPDG*mprPDG);

  Double_t pxv0_init = vka->Px();
  Double_t pyv0_init = vka->Py();
  Double_t pzv0_init = vka->Pz();
  Double_t massv0_init = 0.493677;
  Double_t Ev0_init = sqrt(pxv0_init*pxv0_init+pyv0_init*pyv0_init+pzv0_init*pzv0_init+massv0_init*massv0_init);

  Double_t pxlc_init = pxpr_init+pxv0_init;
  Double_t pylc_init = pypr_init+pyv0_init;
  Double_t pzlc_init = pzpr_init+pzv0_init;
  Double_t Elc_init = Epr_init+Ev0_init;
  Double_t lcmass_init = sqrt(Elc_init*Elc_init-pxlc_init*pxlc_init-pylc_init*pylc_init-pzlc_init*pzlc_init);

  if(lcmass_init<mLcPDG-fProdRoughMassTol || lcmass_init>mLcPDG+fProdRoughMassTol) return kFALSE;
  if(sqrt(pxlc_init*pxlc_init+pylc_init*pylc_init)<fProdRoughPtMin) return kFALSE;

  return kTRUE;
}
//________________________________________________________________________
void AliRDHFCutsLctopK0sfromAODtracks::SetMixingWeights(Int_t nbinpr, Double_t *binspr, Int_t nbink0s, Double_t *binsk0s, Double_t *p0val, Double_t *p1val, Double_t *p2val, Double_t *p3val)
{
  //
  // Set weighting factor for mixing
  //
	fNWeightingProtonBinLimits = nbinpr+1;
	fNWeightingK0sBinLimits = nbink0s+1;
	fNWeightingBins = nbinpr*nbink0s;

	fWeightingProtonBins = new Double_t[fNWeightingProtonBinLimits];
	fWeightingK0sBins = new Double_t[fNWeightingK0sBinLimits];
	fWeight_p0 = new Double_t[fNWeightingBins];
	fWeight_p1 = new Double_t[fNWeightingBins];
	fWeight_p2 = new Double_t[fNWeightingBins];
	fWeight_p3 = new Double_t[fNWeightingBins];

	for(Int_t i=0;i<fNWeightingProtonBinLimits;i++){
		fWeightingProtonBins[i] = binspr[i];
	}
	for(Int_t i=0;i<fNWeightingK0sBinLimits;i++){
		fWeightingK0sBins[i] = binsk0s[i];
	}
	for(Int_t i=0;i<fNWeightingBins;i++){
		fWeight_p0[i] = p0val[i];
		fWeight_p1[i] = p1val[i];
		fWeight_p2[i] = p2val[i];
		fWeight_p3[i] = p3val[i];
	}
	return;
}
//________________________________________________________________________
Double_t AliRDHFCutsLctopK0sfromAODtracks::GetMixingWeight(Double_t dphi, Double_t deta, Double_t pt_pr, Double_t pt_k0s)
{
  //
  // Get weighting factor for mixing
  //
	if(fNWeightingBins==0) return 1.;
	if(dphi>M_PI/2.) return 1.;//does not support away side

	Int_t ibin_pr = -9999;
	for(Int_t i=0;i<fNWeightingProtonBinLimits-1;i++){
		if(fWeightingProtonBins[i]<pt_pr && fWeightingProtonBins[i+1]>pt_pr){
			ibin_pr = i;
			break;
		}
	}
	Int_t ibin_k0s = -9999;
	for(Int_t i=0;i<fNWeightingK0sBinLimits-1;i++){
		if(fWeightingK0sBins[i]<pt_k0s && fWeightingK0sBins[i+1]>pt_k0s){
			ibin_k0s = i;
			break;
		}
	}
	if(ibin_pr<0 || ibin_pr > fNWeightingProtonBinLimits-1) return 1.;
	if(ibin_k0s<0 || ibin_k0s > fNWeightingK0sBinLimits-1) return 1.;

	Int_t ibin_comb = ibin_pr*(fNWeightingK0sBinLimits-1)+ibin_k0s;
	
	Double_t p0 = fWeight_p0[ibin_comb];
	Double_t p1 = fWeight_p1[ibin_comb];
	Double_t p2 = fWeight_p2[ibin_comb];
	Double_t p3 = fWeight_p3[ibin_comb];
	Double_t weight = p0 + p1 *TMath::Gaus(dphi,0.,p2)*TMath::Gaus(deta,0.,p3);

	return weight;
}
//________________________________________________________________________
Bool_t AliRDHFCutsLctopK0sfromAODtracks::TagV0(AliAODTrack *ptrk, AliAODEvent *evt, Int_t ntrk, Double_t &minmass)
{
  //
  // Tag conversion electron tracks
  //
	if(fTagV0MassTol<0.) return kFALSE;

  Double_t mprPDG =  TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  Double_t mpiPDG =  TDatabasePDG::Instance()->GetParticle(211)->Mass();
	Double_t mlamPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();

	Bool_t isv0 = kFALSE;
	Bool_t islam = kFALSE;
	minmass = 9999.;

	Int_t trkid = ptrk->GetID();
	Double_t px1 = ptrk->Px();
	Double_t py1 = ptrk->Py();
	Double_t pz1 = ptrk->Pz();
	Double_t Epr1 = sqrt(px1*px1+py1*py1+pz1*pz1+mprPDG*mprPDG);

	for(Int_t it=0;it<ntrk;it++){
		AliAODTrack *trk2 = (AliAODTrack*) evt->GetTrack(it);
		if(!trk2) continue;
		Int_t trkid2 = trk2->GetID();
		if(trkid==trkid2) continue;
		if(ptrk->Charge()*trk2->Charge()>0) continue;
		if(!ptrk->TestFilterMask(BIT(4))) continue;

		Double_t px2 = trk2->Px();
		Double_t py2 = trk2->Py();
		Double_t pz2 = trk2->Pz();
		Double_t E2 = sqrt(px2*px2+py2*py2+pz2*pz2+mpiPDG*mpiPDG);

		Double_t mass_lam = sqrt(pow(Epr1+E2,2)-pow(px1+px2,2)-pow(py1+py2,2)-pow(pz1+pz2,2));
		Double_t dlam = mass_lam-mlamPDG;
		if(fabs(dlam)<fabs(minmass)){
			minmass = dlam;
			islam = kTRUE;
		}
	}

	if(fabs(minmass)<fTagV0MassTol) isv0 = kTRUE;

	if(islam) minmass += mlamPDG;

  return isv0;
}

//________________________________________________________________________
Bool_t AliRDHFCutsLctopK0sfromAODtracks::TagV0SameSign(AliAODTrack *ptrk, AliAODEvent *evt, Int_t ntrk, Double_t &minmass)
{
  //
  // Tag conversion electron tracks
  //
	if(fTagV0MassTol<0.) return kFALSE;

  Double_t mprPDG =  TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  Double_t mpiPDG =  TDatabasePDG::Instance()->GetParticle(211)->Mass();
	Double_t mlamPDG =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();

	Bool_t isv0 = kFALSE;
	Bool_t islam = kFALSE;
	minmass = 9999.;

	Int_t trkid = ptrk->GetID();
	Double_t px1 = ptrk->Px();
	Double_t py1 = ptrk->Py();
	Double_t pz1 = ptrk->Pz();
	Double_t Epr1 = sqrt(px1*px1+py1*py1+pz1*pz1+mprPDG*mprPDG);

	for(Int_t it=0;it<ntrk;it++){
		AliAODTrack *trk2 = (AliAODTrack*) evt->GetTrack(it);
		if(!trk2) continue;
		Int_t trkid2 = trk2->GetID();
		if(trkid==trkid2) continue;
		if(ptrk->Charge()*trk2->Charge()<0) continue;
		if(!ptrk->TestFilterMask(BIT(4))) continue;

		Double_t px2 = trk2->Px();
		Double_t py2 = trk2->Py();
		Double_t pz2 = trk2->Pz();
		Double_t E2 = sqrt(px2*px2+py2*py2+pz2*pz2+mpiPDG*mpiPDG);

		Double_t mass_lam = sqrt(pow(Epr1+E2,2)-pow(px1+px2,2)-pow(py1+py2,2)-pow(pz1+pz2,2));
		Double_t dlam = mass_lam-mlamPDG;
		if(fabs(dlam)<fabs(minmass)){
			minmass = dlam;
			islam = kTRUE;
		}
	}

	if(fabs(minmass)<fTagV0MassTol) isv0 = kTRUE;

	if(islam) minmass += mlamPDG;

  return isv0;
}
