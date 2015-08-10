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
  fProdTrackPtMin(0.3),
  fProdTrackEtaRange(0.8),
  fProdUseAODFilterBit(kTRUE),
  fProdV0MassTolK0s(0.01),
  fProdV0PtMin(0.5),
  fProdV0CosPointingAngleToPrimVtxMin(0.99),
  fProdV0DcaDaughtersMax(1.5),
  fProdV0DaughterEtaRange(0.8),
  fProdV0DaughterPtMin(0.0),
  fProdV0DaughterTPCClusterMin(70),
fProdRoughMassTol(0.25),
fProdRoughPtMin(0.0)
{
  //
  // Default Constructor
  //

  const Int_t nvars=7;
  SetNVars(nvars);
  TString varNames[nvars]={"Lc inv. mass [GeV/c2]",                   //  0
			   "Lc pT [GeV/c]", //1
			   "Bachelor pT [GeV/c]", //2
			   "Bachelor d0 [cm]", //3
			   "V0 d0 [cm]", //4
			   "K0s mass [GeV/c2]", //5
			   "Decay Length XY [cm]" //6
  };

  Bool_t isUpperCut[nvars]={kTRUE,  //  0
			    kFALSE, //1
			    kFALSE, //2
			    kTRUE, //3
			    kTRUE, //4
			    kTRUE, //5
			    kFALSE //6
  };
  SetVarNames(nvars,varNames,isUpperCut);
  Bool_t forOpt[nvars]={kFALSE, //  0
			kFALSE, //1
			kTRUE, //2
			kTRUE, //3
			kTRUE, //4
			kTRUE, //5
			kTRUE //6
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
  fProdTrackPtMin(source.fProdTrackPtMin),
  fProdTrackEtaRange(source.fProdTrackEtaRange),
  fProdUseAODFilterBit(source.fProdUseAODFilterBit),
  fProdV0MassTolK0s(source.fProdV0MassTolK0s),
  fProdV0PtMin(source.fProdV0PtMin),
  fProdV0CosPointingAngleToPrimVtxMin(source.fProdV0CosPointingAngleToPrimVtxMin),
  fProdV0DcaDaughtersMax(source.fProdV0DcaDaughtersMax),
  fProdV0DaughterEtaRange(source.fProdV0DaughterEtaRange),
  fProdV0DaughterPtMin(source.fProdV0DaughterPtMin),
  fProdV0DaughterTPCClusterMin(source.fProdV0DaughterTPCClusterMin),
  fProdRoughMassTol(source.fProdRoughMassTol),
  fProdRoughPtMin(source.fProdRoughPtMin)
{
  //
  // Copy constructor
  //
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
  fProdTrackPtMin = source.fProdTrackPtMin;
  fProdTrackEtaRange = source.fProdTrackEtaRange;
  fProdUseAODFilterBit = source.fProdUseAODFilterBit;
  fProdV0MassTolK0s = source.fProdV0MassTolK0s;
  fProdV0PtMin = source.fProdV0PtMin;
  fProdV0CosPointingAngleToPrimVtxMin = source.fProdV0CosPointingAngleToPrimVtxMin;
  fProdV0DcaDaughtersMax=source.fProdV0DcaDaughtersMax;
  fProdV0DaughterEtaRange=source.fProdV0DaughterEtaRange;
  fProdV0DaughterPtMin=source.fProdV0DaughterPtMin;
  fProdV0DaughterTPCClusterMin=source.fProdV0DaughterTPCClusterMin;
  fProdRoughMassTol = source.fProdRoughMassTol;
  fProdRoughPtMin = source.fProdRoughPtMin;
  
  
  return *this;
}

//---------------------------------------------------------------------------
AliRDHFCutsLctopK0sfromAODtracks::~AliRDHFCutsLctopK0sfromAODtracks() {
  //
  //  Default Destructor
  //
  
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
  if(fVarsForOpt[2]){
    iter++;
    AliAODTrack *part = dd->GetBachelor();
    vars[iter]= part->Pt();
  }
  if(fVarsForOpt[3]){
    iter++;
    vars[iter]= dd->Getd0Prong(0);
  }
  if(fVarsForOpt[4]){
    iter++;
    vars[iter]= dd->Getd0Prong(1);
  }
  if(fVarsForOpt[5]){
    iter++;
    AliAODv0 *v0 = dd->Getv0();
    vars[iter]= v0->MassK0Short();
  }
  if(fVarsForOpt[6]){
    iter++;
    vars[iter]= dd->DecayLengthXY();
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
    
    Double_t mlcPDG =  TDatabasePDG::Instance()->GetParticle(4122)->Mass();
    Double_t mk0sPDG =  TDatabasePDG::Instance()->GetParticle(310)->Mass();
    AliAODTrack *part = d->GetBachelor();
    AliAODv0 *v0 = d->Getv0();
    
    if(TMath::Abs(d->InvMassLctoK0sP()-mlcPDG) > fCutsRD[GetGlobalIndex(0,ptbin)])
      {
	okcand = kFALSE;
      }
    if(d->Pt() < fCutsRD[GetGlobalIndex(1,ptbin)])
      {
	okcand = kFALSE;
      }
    if(part->Pt() < fCutsRD[GetGlobalIndex(2,ptbin)])
      {
	okcand = kFALSE;
      }
    if(TMath::Abs(d->Getd0Prong(0)) > fCutsRD[GetGlobalIndex(3,ptbin)])
      {
	okcand = kFALSE;
      }
    if(TMath::Abs(d->Getd0Prong(1)) > fCutsRD[GetGlobalIndex(4,ptbin)])
      {
	okcand = kFALSE;
      }
    if(TMath::Abs(v0->MassK0Short()-mk0sPDG) > fCutsRD[GetGlobalIndex(5,ptbin)])
      {
	okcand = kFALSE;
      }
    Double_t lccospaxy = CalculateLcCosPAXY(d);
    if(d->DecayLengthXY() * lccospaxy < fCutsRD[GetGlobalIndex(6,ptbin)])
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
Int_t AliRDHFCutsLctopK0sfromAODtracks::IsSelectedPID(AliAODRecoDecayHF* obj) 
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

  Int_t isProton=fPidHF->MakeRawPid(part,4); 
  if(isProton<1) returnvalue = 0;
  
  return returnvalue;
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
Bool_t AliRDHFCutsLctopK0sfromAODtracks::SingleTrkCuts(AliAODTrack *trk)
{
  //
  // Single Track Cut to be applied before object creation
  //

  if(trk->GetStatus()&AliESDtrack::kITSpureSA) return kFALSE;
  if(!(trk->GetStatus()&AliESDtrack::kITSin)) return kFALSE;
  if(fProdUseAODFilterBit && !trk->TestFilterMask(BIT(4))) return kFALSE;
  if(fabs(trk->Eta())>fProdTrackEtaRange) return kFALSE;
  if(trk->Pt()<fProdTrackPtMin) return kFALSE;

  if(fUsePID)
    {
      if(fPidHF->GetPidResponse()==0x0){
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
	AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
	fPidHF->SetPidResponse(pidResp);
      }

      Int_t isProton=fPidHF->MakeRawPid(trk,4); 
      if(isProton<1) return kFALSE;
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
  if(TMath::Sqrt(dvertx*dvertx+dverty*dverty)==0) return -9999.;
  if(TMath::Sqrt(px*px+py*py)==0) return -9999.;
  Double_t cospaxy = (dvertx*px+dverty*py)/TMath::Sqrt(dvertx*dvertx+dverty*dverty)/TMath::Sqrt(px*px+py*py);
  return (Float_t)cospaxy;
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
