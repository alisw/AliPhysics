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
//                                                                       //
//     Analysis task for computing single particle efficiencies          //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TH1D.h>
#include <TList.h>
#include <THashList.h>
#include <TString.h>

#include <AliPID.h>
#include <AliCFContainer.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliESDv0.h>
#include <AliESDv0Cuts.h>
#include <AliESDv0KineCuts.h>
#include <AliKFVertex.h>
#include <AliESDfriendTrack.h>
#include <AliESDfriend.h>
#include <AliTRDseedV1.h>
#include <AliTRDcluster.h>
#include <AliTRDtrackV1.h>

#include "AliDielectron.h"
#include "AliDielectronHistos.h"
#include "AliDielectronPair.h"
#include "AliDielectronCF.h"
#include "AliDielectronMC.h"
#include "AliDielectronHistos.h"
#include "AliAnalysisTaskSingleParticle.h"

#include <iostream>
using namespace std;

ClassImp(AliAnalysisTaskSingleParticle)

//_________________________________________________________________________________
AliAnalysisTaskSingleParticle::AliAnalysisTaskSingleParticle() :
  AliAnalysisTaskSE(),
  fCfContainer(0x0),
  fHistos(0x0),
  fHistogramList(),
  fSelectPhysics(kTRUE),
  fTriggerMask(AliVEvent::kMB),
  fRejectPileup(kFALSE),
  fFillTRDfriendPH(kFALSE),
  fEventFilter(0x0),
  fTrackFilter(0x0),
  fPairFilter(0x0),
  fV0Cuts(0x0),
  fLambdaFilter(0x0),
  fK0sFilter(0x0),
  fV0KineCuts(0x0),
  fCFNVars(0),
  fEventStat(0x0)
{
  //
  // Constructor
  //
  for(Int_t i=0; i<kNMaxDimensions; ++i) {
    fCFVarsEnabled[i] = -1; fCFVarsNbins[i] = 0;
    fCFVarRanges[i][0] = 0.0; fCFVarRanges[i][1] = -1.0;
    fCFVarBins[i] = "";
  }
}

//_________________________________________________________________________________
AliAnalysisTaskSingleParticle::AliAnalysisTaskSingleParticle(const char *name) :
  AliAnalysisTaskSE(name),
  fCfContainer(0x0),
  fHistos(0x0),
  fHistogramList(),
  fSelectPhysics(kTRUE),
  fTriggerMask(AliVEvent::kMB),
  fRejectPileup(kFALSE),
  fFillTRDfriendPH(kFALSE),
  fEventFilter(0x0),
  fTrackFilter(0x0),
  fPairFilter(0x0),
  fV0Cuts(0x0),
  fLambdaFilter(0x0),
  fK0sFilter(0x0),
  fV0KineCuts(0x0),
  fCFNVars(0),
  fEventStat(0x0)
{
  //
  // Constructor
  //
  for(Int_t i=0; i<kNMaxDimensions; ++i) {
    fCFVarsEnabled[i] = -1; fCFVarsNbins[i] = 0;
    fCFVarRanges[i][0] = 0.0; fCFVarRanges[i][1] = -1.0;
    fCFVarBins[i] = "";
  }
  fHistogramList.SetName("QAhistos");
  fHistogramList.SetOwner();
  DefineInput(0, TChain::Class());
  DefineOutput(1, AliCFContainer::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TH1D::Class());
}


//____________________________________________________________________
AliAnalysisTaskSingleParticle::~AliAnalysisTaskSingleParticle()
{
  //
  // Destructor
  // 
  fHistogramList.SetOwner(kFALSE);
  if(fHistos) delete fHistos;
  delete fCfContainer;
}


//____________________________________________________________________
void AliAnalysisTaskSingleParticle::AddCFVar(Int_t var, Int_t nbins, Double_t lowLim, Double_t highLim) {
  //
  // Configure variables for the CF container
  //
  if(fCFNVars>=kNMaxDimensions) return;
  fCFVarsEnabled[fCFNVars] = var;
  fCFVarsNbins[fCFNVars] = nbins;
  fCFVarRanges[fCFNVars][0] = lowLim;
  fCFVarRanges[fCFNVars][1] = highLim;
  ++fCFNVars;
}


//____________________________________________________________________
void AliAnalysisTaskSingleParticle::AddCFVar(Int_t var, const Char_t* bins) {
  //
  // Configure variables for the CF container
  //
  if(fCFNVars>=kNMaxDimensions) return;
  fCFVarsEnabled[fCFNVars] = var;
  fCFVarBins[fCFNVars] = bins;
  ++fCFNVars;
}


//_________________________________________________________________________________
void AliAnalysisTaskSingleParticle::UserCreateOutputObjects()
{
  //
  // Initialize CF container
  //
  if(fCfContainer) return;
  if(!fHistogramList.IsEmpty()) return;
  
  if (!fEventStat){
    fEventStat=new TH1D("hEventStat","Event statistics",4,-0.5,3.5);
    fEventStat->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(3,"After event filter");
    fEventStat->GetXaxis()->SetBinLabel(4,"After pileup rejection");
  }
  
  //Bool_t hasMC = AliDielectronMC::Instance()->HasMC();
  
  if(!fCfContainer) {
    for(Int_t ivar=0; ivar<fCFNVars; ++ivar) {
      if(fCFVarBins[ivar][0]!='\0') {
        TObjArray* arr = fCFVarBins[ivar].Tokenize(",");
        fCFVarsNbins[ivar] = arr->GetEntries()-1;
      }
    }
    fCFVarsNbins[fCFNVars] = AliPID::kSPECIES;
    //
    // If fFillTRDfriendPH is toggled then add some predefined variables to the CF
    // Hardcoded variables and binning !! TODO: Implement them in the VarManager ??
    // This container is filled for every TRD time bin, and every tracklet from a track
    Double_t phLimits[2] = {0.0, 2000.}; Int_t phBins = 200;
    Double_t detLimits[2] = {-0.5, 539.5}; Int_t detBins = 540;
    const Int_t kNTrdMomBins = 17;
    Double_t trdMomLimits[kNTrdMomBins+1] = {0.0, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
    Int_t trdNtb = 30; Double_t trdTbLims[2] = {0.5, 30.5};
    Int_t trdQtotNbins = 200; Double_t trdQtotLims[2] = {0.0, 10000.};
    fCFVarsNbins[fCFNVars+1] = kNTrdMomBins;
    fCFVarsNbins[fCFNVars+2] = detBins;
    fCFVarsNbins[fCFNVars+3] = trdNtb;
    fCFVarsNbins[fCFNVars+4] = phBins;
    fCFVarsNbins[fCFNVars+5] = trdQtotNbins;
    
    if(fCFNVars>0 || fFillTRDfriendPH) {
      fCfContainer = new AliCFContainer("SingleParticle", "Single Particle CF", 1, (fFillTRDfriendPH ? fCFNVars+6 : fCFNVars+1), fCFVarsNbins);
      for(Int_t ivar=0; ivar<fCFNVars; ++ivar) {
        if(fCFVarBins[ivar][0]=='\0')
	  fCfContainer->SetBinLimits(ivar, fCFVarRanges[ivar][0], fCFVarRanges[ivar][1]);
        else {
	  TObjArray* arr = fCFVarBins[ivar].Tokenize(",");
	  if(arr->GetEntries()-1>0) {
	    Double_t* binLims = new Double_t[arr->GetEntries()];
	    for(Int_t ib=0;ib<arr->GetEntries();++ib) {
	      TString binStr = arr->At(ib)->GetName();
	      binLims[ib] = binStr.Atof();
	    }
	    fCfContainer->SetBinLimits(ivar, binLims);
	  }
        }
        fCfContainer->SetVarTitle(ivar, AliDielectronVarManager::GetValueName(fCFVarsEnabled[ivar]));
      }
      fCfContainer->SetBinLimits(fCFNVars, -0.5+AliPID::kElectron, 0.5+AliPID::kProton);
      fCfContainer->SetVarTitle(fCFNVars, "specie");
      if(fFillTRDfriendPH) {
	fCfContainer->SetBinLimits(fCFNVars+1, trdMomLimits);
	fCfContainer->SetVarTitle(fCFNVars+1, "TRD_trackletMom");
	fCfContainer->SetBinLimits(fCFNVars+2, detLimits[0], detLimits[1]);
	fCfContainer->SetVarTitle(fCFNVars+2, "TRD_detector");
	fCfContainer->SetBinLimits(fCFNVars+3, trdTbLims[0], trdTbLims[1]);
	fCfContainer->SetVarTitle(fCFNVars+3, "TRD_tb");
	fCfContainer->SetBinLimits(fCFNVars+4, phLimits[0], phLimits[1]);
	fCfContainer->SetVarTitle(fCFNVars+4, "TRD_PH");
	fCfContainer->SetBinLimits(fCFNVars+5, trdQtotLims[0], trdQtotLims[1]);
	fCfContainer->SetVarTitle(fCFNVars+5, "TRD_Qtot");
      }
      fCfContainer->SetStepTitle(0, "Tracking");
    }
  }
  
  if(fHistos) {
    fHistogramList.Add(const_cast<THashList*>(fHistos->GetHistogramList()));
  }
  
  if(fCfContainer) PostData(1, fCfContainer);
  if(fHistos) PostData(2, &fHistogramList);
  PostData(3, fEventStat);
}

//_________________________________________________________________________________
void AliAnalysisTaskSingleParticle::UserExec(Option_t *)
{
  //
  // Main loop. Called for every event
  //

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  //Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  //Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if (!inputHandler) return;
  
  if ( inputHandler->GetPIDResponse() ){
    AliDielectronVarManager::SetPIDResponse( inputHandler->GetPIDResponse() );
  } else {
    AliFatal("This task needs the PID response attached to the input event handler!");   
  } 
  // Was event selected ?
  UInt_t isSelected = AliVEvent::kAny;
  if( fSelectPhysics && inputHandler && inputHandler->GetEventSelection() ) {
    isSelected = inputHandler->IsEventSelected();
    isSelected&=fTriggerMask;
  }
  
  //Before physics selection
  fEventStat->Fill(0);
  if (isSelected==0) {
    PostData(3,fEventStat);
    return;
  }
  //after physics selection
  fEventStat->Fill(1);

  //event filter
  if (fEventFilter) {
    if (!fEventFilter->IsSelected(InputEvent())) {PostData(3, fEventStat); return;}
  }
  fEventStat->Fill(2);
  
  //pileup
  if (fRejectPileup){
    if (InputEvent()->IsPileupFromSPD(3,0.8,3.,2.,5.)) {PostData(3, fEventStat); return;}
  }
  fEventStat->Fill(3);
  
  if(!fCfContainer) {PostData(3, fEventStat); return;}
  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );

  AliDielectronPID::SetCorrVal((Double_t)InputEvent()->GetRunNumber());
  
  AliESDEvent* esd = (AliESDEvent*)InputEvent();
  AliDielectronVarManager::SetEvent(esd);
  
  Double_t valuesPos[AliDielectronVarManager::kNMaxValues];
  Double_t valuesNeg[AliDielectronVarManager::kNMaxValues];
  Double_t valuesPair[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::Fill(esd, valuesPos);
  AliDielectronVarManager::Fill(esd, valuesNeg);
  AliDielectronVarManager::Fill(esd, valuesPair);
  
  if(fHistos && fHistos->GetHistogramList()->FindObject("Event"))
    fHistos->FillClass("Event", AliDielectronVarManager::kNMaxValues, valuesPos);
  
  const AliESDVertex *primaryVertex = esd->GetPrimaryVertex();
  AliKFVertex primaryVertexKF(*primaryVertex);
    
  fV0KineCuts->SetEvent(esd);
  fV0KineCuts->SetPrimaryVertex(&primaryVertexKF);
  
  Int_t pdgV0=0; Int_t pdgP=0; Int_t pdgN=0;
  for(Int_t iV0=0; iV0<InputEvent()->GetNumberOfV0s(); ++iV0) {   // loop over V0s
    //cout << "iV0 = " << iV0 << endl;
    AliESDv0 *v0 = esd->GetV0(iV0);
       
    AliESDtrack* legPos = esd->GetTrack(v0->GetPindex());
    AliESDtrack* legNeg = esd->GetTrack(v0->GetNindex());
 
    if(legPos->GetSign() == legNeg->GetSign()) {
      //cout << "V0 rejected: same sign legs" << endl;
      continue;
    }

    //To avoid ghosts test of trunk for the refit
    if(!(legPos->GetStatus() & AliESDtrack::kTPCrefit)) {
      //cout << "pos leg rejected: no TPCrefit" << endl;
      continue;
    }
    if(!(legNeg->GetStatus() & AliESDtrack::kTPCrefit)) {
      //cout << "neg leg rejected: no TPCrefit" << endl;
      continue;
    }
    
    Bool_t v0ChargesAreCorrect = (legPos->GetSign()==+1 ? kTRUE : kFALSE);
    legPos = (!v0ChargesAreCorrect ? esd->GetTrack(v0->GetNindex()) : legPos);
    legNeg = (!v0ChargesAreCorrect ? esd->GetTrack(v0->GetPindex()) : legNeg);
    
    Bool_t good = fV0KineCuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
    if(!good) {
      //cout << "V0KineCuts rejected V0" << endl;
      continue;
    }
    
    if(!fFillTRDfriendPH) {
      //cout << "No TRD friend filling" << endl;
      if(pdgV0!=22) continue;
      if(TMath::Abs(pdgN)!=11) continue;
      if(TMath::Abs(pdgP)!=11) continue;
    }
    /*
    // V0 standard cuts
    if(fV0Cuts) {
      TList v0List;
      v0List.Add(v0);
      v0List.Add(legPos); v0List.Add(legNeg);
      v0List.Add(const_cast<AliESDVertex*>(primaryVertex));
      if(!fV0Cuts->IsSelected(&v0List)) continue;
    }*/
    
    // additional track filtering
    if(fTrackFilter && !fTrackFilter->IsSelected(legPos)) {
      //cout << "pos leg rejected: track filter" << endl;
      continue;
    }
    if(fTrackFilter && !fTrackFilter->IsSelected(legNeg)) {
      //cout << "neg leg rejected: track filter" << endl;
      continue;
    }
    
    AliKFParticle* posKF = (v0ChargesAreCorrect ? new AliKFParticle(*(v0->GetParamP()),pdgP) : new AliKFParticle(*(v0->GetParamN()),pdgN));
    AliKFParticle* negKF = (v0ChargesAreCorrect ? new AliKFParticle(*(v0->GetParamN()),pdgN) : new AliKFParticle(*(v0->GetParamP()),pdgP));
    AliDielectronPair pair;
    pair.SetTracks(posKF, negKF, legPos, legNeg);
    pair.SetType(1);
    if(fPairFilter && !fPairFilter->IsSelected(&pair)) {
      //cout << "diele pair rejected: pair filter" << endl;
      continue;
    }
    
    // additional filtering on the KF pair (mass, chi2, etc.)
    // Gamma conversion inclusion cuts
    /*
    AliKFParticle* posKFGammaEle = (v0ChargesAreCorrect ? new AliKFParticle(*(v0->GetParamP()),-11) : new AliKFParticle(*(v0->GetParamN()),-11));
    AliKFParticle* negKFGammaEle = (v0ChargesAreCorrect ? new AliKFParticle(*(v0->GetParamN()),+11) : new AliKFParticle(*(v0->GetParamP()),+11));
    AliDielectronPair gammaPair;
    gammaPair.SetTracks(posKFGammaEle, negKFGammaEle, legPos, legNeg);
    gammaPair.SetType(1);
    if(fPairFilter && !fPairFilter->IsSelected(&gammaPair)) continue;
    */
    /*
    // Lambda exclusion cuts
    AliKFParticle* posKFLambdaPro = (v0ChargesAreCorrect ? new AliKFParticle(*(v0->GetParamP()),2212) : new AliKFParticle(*(v0->GetParamN()),2212));
    AliKFParticle* negKFLambdaPio = (v0ChargesAreCorrect ? new AliKFParticle(*(v0->GetParamN()),-211) : new AliKFParticle(*(v0->GetParamP()),-211));
    AliDielectronPair lambdaPair;
    lambdaPair.SetTracks(posKFLambdaPro, negKFLambdaPio, legPos, legNeg);
    lambdaPair.SetType(1);
    if(fLambdaFilter && fLambdaFilter->IsSelected(&lambdaPair)) continue;
    // Anti-Lambda exclusion cuts
    AliKFParticle* posKFALambdaPio = (v0ChargesAreCorrect ? new AliKFParticle(*(v0->GetParamP()),211) : new AliKFParticle(*(v0->GetParamN()),211));
    AliKFParticle* negKFALambdaPro = (v0ChargesAreCorrect ? new AliKFParticle(*(v0->GetParamN()),-2212) : new AliKFParticle(*(v0->GetParamP()),-2212));
    AliDielectronPair alambdaPair;
    alambdaPair.SetTracks(posKFALambdaPio, negKFALambdaPro, legPos, legNeg);
    alambdaPair.SetType(1);
    if(fLambdaFilter && fLambdaFilter->IsSelected(&alambdaPair)) continue;
    // K0s exclusion cuts
    AliKFParticle* posKFK0sPio = (v0ChargesAreCorrect ? new AliKFParticle(*(v0->GetParamP()),211) : new AliKFParticle(*(v0->GetParamN()),211));
    AliKFParticle* negKFK0sPio = (v0ChargesAreCorrect ? new AliKFParticle(*(v0->GetParamN()),-211) : new AliKFParticle(*(v0->GetParamP()),-211));
    AliDielectronPair k0sPair;
    k0sPair.SetTracks(posKFK0sPio, negKFK0sPio, legPos, legNeg);
    k0sPair.SetType(1);
    if(fK0sFilter && fK0sFilter->IsSelected(&k0sPair)) continue;
    */
    
    //  Fill histograms and the CF container
    AliDielectronVarManager::Fill(legPos, valuesPos);
    AliDielectronVarManager::Fill(legNeg, valuesNeg);
    
    if(fHistos && fHistos->GetHistogramList()->FindObject("Track_Pos"))
      fHistos->FillClass("Track_Pos", AliDielectronVarManager::kNMaxValues, valuesPos);
    if(fHistos && fHistos->GetHistogramList()->FindObject("Track_Neg"))
      fHistos->FillClass("Track_Neg", AliDielectronVarManager::kNMaxValues, valuesNeg);
    
    AliDielectronVarManager::Fill(&pair, valuesPair);
    if(fHistos && fHistos->GetHistogramList()->FindObject(Form("Pair_%s",AliDielectron::PairClassName(1))))
      fHistos->FillClass(Form("Pair_%s",AliDielectron::PairClassName(1)), AliDielectronVarManager::kNMaxValues, valuesPair);
    
    if(valuesPos[AliDielectronVarManager::kPOut]>=0.5)
      FillContainer(0, valuesPos, (v0ChargesAreCorrect ? v0->GetPindex() : v0->GetNindex()), (v0ChargesAreCorrect ? pdgP : pdgN));
    if(valuesNeg[AliDielectronVarManager::kPOut]>=0.5)
      FillContainer(0, valuesNeg, (v0ChargesAreCorrect ? v0->GetNindex() : v0->GetPindex()), (v0ChargesAreCorrect ? pdgN : pdgP));
  }   // end loop over online V0s

  //delete gammaPair;

  PostData(1, fCfContainer);
  PostData(2, &fHistogramList);
  PostData(3, fEventStat);
}

//_________________________________________________________________________________
void AliAnalysisTaskSingleParticle::FinishTaskOutput()
{
  //
  //
  //
}


//_________________________________________________________________________________
void AliAnalysisTaskSingleParticle::FillContainer(Int_t step, const Double_t* values, Int_t trackId, Int_t pdg) {
  //
  // Fill the CF container
  //
  Double_t valuesCF[kNMaxDimensions];
  for(Int_t i=0; i<fCFNVars; ++i) valuesCF[i] = values[fCFVarsEnabled[i]];
  
  Double_t pid = -1.0;
  if(TMath::Abs(pdg)==11)   pid = AliPID::kElectron;
  if(TMath::Abs(pdg)==13)   pid = AliPID::kMuon;
  if(TMath::Abs(pdg)==211)  pid = AliPID::kPion;
  if(TMath::Abs(pdg)==321)  pid = AliPID::kKaon;
  if(TMath::Abs(pdg)==2212) pid = AliPID::kProton;
  valuesCF[fCFNVars] = pid;
  
  if(fFillTRDfriendPH) {
    AliESDfriendTrack* friendTrk = fESDfriend->GetTrack(trackId);
    if(!friendTrk) return;
    if(values[AliDielectronVarManager::kTRDntracklets]>0.001) {
      TObject* o=0x0; Int_t ical = 0; 
      AliTRDtrackV1* trdTrack=0x0;
      while(1) {
	o = friendTrk->GetCalibObject(ical++);
	if(!o) break;
	TString objName = o->IsA()->GetName();
	if(!objName.Contains("AliTRDtrackV1")) continue;
	trdTrack = dynamic_cast<AliTRDtrackV1*>(o);
	break;
      }
      if(!trdTrack) return;
      
      Double_t charge = 0.0;
      for(Int_t ipl = 0; ipl < 6; ipl++) {  // loop over TRD layers
	AliTRDseedV1 *tracklet = trdTrack->GetTracklet(ipl);
	if(!tracklet) continue;
	if(!tracklet->IsOK()) continue;
	valuesCF[fCFNVars+1] = tracklet->GetMomentum();
	valuesCF[fCFNVars+2] = tracklet->GetDetector();
	charge = 0.0;
	for(Int_t itb = 0; itb < AliTRDseedV1::kNtb; itb++){
	  AliTRDcluster *c = tracklet->GetClusters(itb);
	  AliTRDcluster *c1 = tracklet->GetClusters(itb + AliTRDseedV1::kNtb);  // second pad row
	  if(c) charge += TMath::Abs(c->GetQ()); //
	  if(c1) charge += TMath::Abs(c1->GetQ());
	}
	valuesCF[fCFNVars+5] = charge;
	
	for(Int_t itb = 0; itb < AliTRDseedV1::kNtb; itb++){
	  AliTRDcluster *c = tracklet->GetClusters(itb);
	  AliTRDcluster *c1 = tracklet->GetClusters(itb + AliTRDseedV1::kNtb);  // second pad row
	  if(!(c || c1)) continue;
	  AliTRDcluster *cptr = (c ? c : c1);
	  Int_t tcal = cptr->GetLocalTimeBin();
	  Float_t sig = 0;
	  if(c) sig += TMath::Abs(c->GetQ()); //
	  if(c1) sig += TMath::Abs(c1->GetQ());
	  valuesCF[fCFNVars+3] = tcal;
	  valuesCF[fCFNVars+4] = sig;
	  
	  fCfContainer->Fill(valuesCF, step);
	}        // end loop over time bins
      }       // end loop over tracklets
    }      // end if track has TRD tracklets
  }     // end if fill TRD friend PH
  if(!fFillTRDfriendPH)
    fCfContainer->Fill(valuesCF, step);
}

