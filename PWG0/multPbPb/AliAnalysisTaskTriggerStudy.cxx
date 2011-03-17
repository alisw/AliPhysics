// AliAnalysisTaskTriggerStudy

// Author: Michele Floris, CERN
// TODO:
// - Add chi2/cluster plot for primary, secondaries and fakes


#include "AliAnalysisTaskTriggerStudy.h"
#include "AliESDInputHandler.h"
#include "AliHistoListWrapper.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TH1I.h"
#include "TH3D.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliESDCentrality.h"

#include <iostream>
#include "AliTriggerAnalysis.h"
#include "AliMultiplicity.h"
#include "TFile.h"
#include "AliLog.h"
#include "AliESDtrackCuts.h"
#include "AliESDVZERO.h"
#include "TH2F.h"
#include "AliESDUtils.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"

using namespace std;

ClassImp(AliAnalysisTaskTriggerStudy)

//const char * AliAnalysisTaskTriggerStudy::kVDNames[] = {"C0SM1","C0SM2","C0VBA","C0VBC","C0OM2"};       
const char * AliAnalysisTaskTriggerStudy::kVDNames[] = {"V0AND","V0OR","NTRACKS", "NTRACKS ESD"};//,"C0OM2"};       

AliAnalysisTaskTriggerStudy::AliAnalysisTaskTriggerStudy()
: AliAnalysisTaskSE("TaskTriggerStudy"),
  fESD(0),fHistoList(0),fIsMC(0),fTriggerAnalysis(0),fHistoSuffix(""),fNTrackletsCut(1000000),fNTrackletsCutKine(100),fRejectBGWithV0(0)
{
  // constructor

  DefineOutput(1, AliHistoListWrapper::Class());

}
AliAnalysisTaskTriggerStudy::AliAnalysisTaskTriggerStudy(const char * name)
  : AliAnalysisTaskSE(name),
    fESD(0),fHistoList(0),fIsMC(0),fTriggerAnalysis(0),fHistoSuffix(""),fNTrackletsCut(1000000),fNTrackletsCutKine(100),fRejectBGWithV0(0)
{
  //
  // Standard constructur which should be used
  //

  DefineOutput(1, AliHistoListWrapper::Class());

}

AliAnalysisTaskTriggerStudy::AliAnalysisTaskTriggerStudy(const AliAnalysisTaskTriggerStudy& obj) : 
  AliAnalysisTaskSE(obj) ,fESD (0), fIsMC(0), fTriggerAnalysis(0),fHistoSuffix(""),fNTrackletsCut(1000000),fNTrackletsCutKine(100),fRejectBGWithV0(0)
{
  //copy ctor
  fESD = obj.fESD ;
  fHistoList = obj.fHistoList;
  fTriggerAnalysis = obj.fTriggerAnalysis;
  fHistoSuffix = obj.fHistoSuffix;
  fNTrackletsCut = obj.fNTrackletsCut;
  fNTrackletsCutKine = obj.fNTrackletsCutKine;
  fRejectBGWithV0 = obj.fRejectBGWithV0;
}

AliAnalysisTaskTriggerStudy::~AliAnalysisTaskTriggerStudy(){
  // destructor

  if(!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    if(fHistoList) {
      delete fHistoList;
      fHistoList = 0;
    }
    if(fTriggerAnalysis) {
      delete fTriggerAnalysis;
      fHistoList = 0;
    }
  }
  // Histo list should not be destroyed: fListWrapper is owner!

}
void AliAnalysisTaskTriggerStudy::UserCreateOutputObjects()
{
  // Called once
  fHistoList = new AliHistoListWrapper("histoList","histogram list for trigger studies");
  fTriggerAnalysis = new AliTriggerAnalysis();
  if (fIsMC) fTriggerAnalysis->SetAnalyzeMC(1);
}


void AliAnalysisTaskTriggerStudy::UserExec(Option_t *)
{
  // User code

  // FIXME: make sure you have the right cuts here

  /* PostData(0) is taken care of by AliAnalysisTaskSE */
  PostData(1,fHistoList);

  fESD = dynamic_cast<AliESDEvent*>(fInputEvent);
  if (strcmp(fESD->ClassName(),"AliESDEvent")) {
    AliFatal("Not processing ESDs");
  }

  // get the multiplicity object
  const AliMultiplicity* mult = fESD->GetMultiplicity();
  Int_t ntracklets = mult->GetNumberOfTracklets();
  // Get Number of tracks
  Int_t ntracks    = AliESDtrackCuts::GetReferenceMultiplicity(fESD,kTRUE); // tpc only
  
  // Get V0 Multiplicity
  AliESDVZERO* esdV0 = fESD->GetVZEROData();
  Float_t multV0A=esdV0->GetMTotV0A();
  Float_t multV0C=esdV0->GetMTotV0C();
  Float_t dummy = 0;  
  Float_t multV0 =  AliESDUtils::GetCorrV0(fESD, dummy);

  // Get number of clusters in layer 1
  Float_t outerLayerSPD = mult->GetNumberOfITSClusters(1);  
  Float_t innerLayerSPD = mult->GetNumberOfITSClusters(0);  
  Float_t totalClusSPD = outerLayerSPD+innerLayerSPD;

  if ( !fIsMC &&(!(fESD->IsTriggerClassFired("CMBS2A-B-NOPF-ALL")|| fESD->IsTriggerClassFired("CMBS2C-B-NOPF-ALL") || fESD->IsTriggerClassFired("CMBAC-B-NOPF-ALL")) )) return;
  GetHistoSPD1  ("All", "All events before any selection")->Fill(outerLayerSPD);
  GetHistoV0M   ("All", "All events before any selection")->Fill(multV0);

  // V0 triggers
  Bool_t c0v0A       = fTriggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kV0A);
  Bool_t c0v0C       = fTriggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kV0C);
  Bool_t v0AHW     = (fTriggerAnalysis->V0Trigger(fESD, AliTriggerAnalysis::kASide, kTRUE) == AliTriggerAnalysis::kV0BB);// should replay hw trigger
  Bool_t v0CHW     = (fTriggerAnalysis->V0Trigger(fESD, AliTriggerAnalysis::kCSide, kTRUE) == AliTriggerAnalysis::kV0BB);// should replay hw trigger

  Bool_t v0ABG = fTriggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kV0ABG);
  Bool_t v0CBG = fTriggerAnalysis->IsOfflineTriggerFired(fESD, AliTriggerAnalysis::kV0CBG);
  Bool_t v0BG  = v0ABG || v0CBG;

  // At least one track in eta < 0.8 with pt > 0.5
  // MC Checks
  Bool_t atLeast1Track = kFALSE;
  Bool_t isSD = kFALSE;
  Bool_t isND = kFALSE;

  if(fIsMC) {
    if (!fMCEvent) {
      AliError("No MC info found");
    } else {
      Int_t nMCTracks = fMCEvent->GetNumberOfTracks();
      Int_t nPhysicalPrimaries = 0;
      for (Int_t ipart=0; ipart<nMCTracks; ipart++) { 	
	AliMCParticle *mcPart  = (AliMCParticle*)fMCEvent->GetTrack(ipart);
	// We don't care about neutrals and non-physical primaries
	if(mcPart->Charge() == 0) continue;

	//check if current particle is a physical primary
	if(!fMCEvent->IsPhysicalPrimary(ipart)) continue;
	if(mcPart->Pt()<0.5) continue;
	if(TMath::Abs(mcPart->Eta())>0.8) continue;
	atLeast1Track = kTRUE;
	break;
      }
      
      AliGenPythiaEventHeader * headPy  = 0;
      AliGenDPMjetEventHeader * headPho = 0;
      AliGenEventHeader * htmp = fMCEvent->GenEventHeader();
      if(!htmp) {
	AliFatal("Cannot Get MC Header!!");
	return;
      }
      if( TString(htmp->IsA()->GetName()) == "AliGenPythiaEventHeader") {
	headPy =  (AliGenPythiaEventHeader*) htmp;
      } else if (TString(htmp->IsA()->GetName()) == "AliGenDPMjetEventHeader") {
	headPho = (AliGenDPMjetEventHeader*) htmp;
      } else {
	AliFatal("Unknown header");	
      }
      if(headPy)   {
	//    cout << "Process: " << headPy->ProcessType() << endl;
	if(headPy->ProcessType() == 92 || headPy->ProcessType() == 93) {
	  isSD = kTRUE; // is single difractive
	}
	if(headPy->ProcessType() != 92 && headPy->ProcessType() != 93 && headPy->ProcessType() != 94) {     
	  isND = kTRUE; // is non-diffractive
	}
	
      } else if (headPho) {
	if(headPho->ProcessType() == 5 || headPho->ProcessType() == 6 ) {
	  isSD = kTRUE;
	}       
	if(headPho->ProcessType() != 5 && headPho->ProcessType() != 6  && headPho->ProcessType() != 7 ) {
	  isND = kTRUE;
	}       
      }
      
    }
  }
  
  static AliESDtrackCuts * cuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
  // cuts->SetPtRange(0.5,10000);
  // cuts->SetEtaRange(-0.8, 0.8);
  Int_t ntracksLoop = fESD->GetNumberOfTracks();

  Bool_t atLeast1TrackESD = kFALSE;
  for (Int_t iTrack = 0; iTrack<ntracksLoop; iTrack++) {    
    AliESDtrack * track = dynamic_cast<AliESDtrack*>(fESD->GetTrack(iTrack));
    // for each track
    
    // track quality cuts
    if(!cuts->AcceptTrack(track)) continue;
    if(track->Pt()<0.5) continue;
    if(TMath::Abs(track->Eta())>0.8) continue;

    atLeast1TrackESD = kTRUE;
    break;
  }



  Bool_t vdArray[kNVDEntries];
  vdArray[kVDV0AND]    = c0v0A && c0v0C;
  vdArray[kVDV0OR]     = c0v0A || c0v0C;
  vdArray[kVDNTRACKS]  = atLeast1Track;
  //  vdArray[kVDNTRACKSESD]  = atLeast1TrackESD;

  FillTriggerOverlaps("All", "All Events",vdArray);
  if(!isSD)   FillTriggerOverlaps("NSD", "NSD Events",vdArray);


  // Physics selection
  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  if(!isSelected) return;


  GetHistoSPD1  ("AllPSNoPrino", "All events after physsel (no ZDC time)")->Fill(outerLayerSPD);
  GetHistoV0M   ("AllPSNoPrino", "All events after physsel (no ZDC time)")->Fill(multV0);


  // Francesco's cuts
  // const AliESDVertex * vtxESDTPC= fESD->GetPrimaryVertexTPC(); 
  // if(vtxESDTPC->GetNContributors()<1) return;
  // if (vtxESDTPC->GetNContributors()<(-10.+0.25*fESD->GetMultiplicity()->GetNumberOfITSClusters(0)))     return;
  // const AliESDVertex * vtxESDSPD= fESD->GetPrimaryVertexSPD(); 
  // Float_t tpcContr=vtxESDTPC->GetNContributors();

  

  // GetT0 Stuff
  const Double32_t *meanT0   = fESD->GetT0TOF();
  const Double32_t  meanT0A  = 0.001* meanT0[1];
  const Double32_t  meanT0C  = 0.001* meanT0[2];
  const Double32_t  meanT0AC = 0.001* meanT0[0];
  Double32_t  T0Vertex = fESD->GetT0zVertex();
  //  Double32_t  *ampT0 =Esdevent ->GetT0amplitude();

//   cut1yesF = ( (meanC < 95. && meanA < 95.) && (meanC < -2.) ) && francescoscut
// cut1notF = ( (meanC < 95. && meanA < 95.) && (meanC < -2.) ) && ! francescoscut
// cut2 = ( (meanC < 95. && meanA < 95.) && ( (meanC-meanA) <=-0.7) && meanC > -2) )
// cut3 = ( (meanC < 95. && meanA < 95.) && ( (meanC-meanA) < 0.7 && (meanC-meanA) > -0.7 ) )
//   cut4 = ( (meanC < 95. && meanA < 95.) && (meanA < -2.)

  Bool_t cut1T0 =  ( (meanT0C < 95. && meanT0A < 95.) && (meanT0C < -2.) );
  Bool_t cut2T0 = ( (meanT0C < 95. && meanT0A < 95.) && ( (meanT0C-meanT0A) <=-0.7) && meanT0C > -2) ;

  if(ntracklets > fNTrackletsCut) return;

  // Reset histo suffix and fill reference histograms without any suffix
  fHistoSuffix = "";

  // Fast or in the outer layer  
  Int_t nFastOrOnline  = fTriggerAnalysis->SPDFiredChips(fESD, 1, 0, 2); // offline
  Int_t nFastOrOffline = fTriggerAnalysis->SPDFiredChips(fESD, 0, 0, 2); // online
  
  Bool_t c0sm1 = nFastOrOffline >= 1;
  Bool_t c0sm2 = nFastOrOffline >= 2;
  Bool_t c0sm3 = nFastOrOffline >= 3;
  Bool_t c0sm4 = nFastOrOffline >= 4;
  Bool_t c0sm5 = nFastOrOffline >= 5;
  

  // TOF triggers 
  // FIXME: move to triggeranalysis?
  AliESDHeader*h = fESD->GetHeader(); // taken the header from AliESDEvent 
  Bool_t c0OM2 = h->IsTriggerInputFired("0OM2"); // thr >= 2 (input 19)
  Bool_t c0OM3 = h->IsTriggerInputFired("0OM3"); // thr >= 3 (input 20)

  // ZDC triggers
  Bool_t zdcA    = kFALSE;
  Bool_t zdcC    = kFALSE;
  Bool_t zdcBar  = kFALSE;
  Bool_t zdcTime = fTriggerAnalysis->ZDCTimeTrigger(fESD);

  if (!fIsMC) {
    // If it's data, we use the TDCs
    zdcA   = fTriggerAnalysis->ZDCTDCTrigger(fESD, AliTriggerAnalysis::kASide, kTRUE, kFALSE) ; 
    zdcC   = fTriggerAnalysis->ZDCTDCTrigger(fESD, AliTriggerAnalysis::kCSide, kTRUE, kFALSE) ;			
    zdcBar = fTriggerAnalysis->ZDCTDCTrigger(fESD, AliTriggerAnalysis::kCentralBarrel) ;				
  } else {
    // If it's MC, we use the energy
    Double_t minEnergy = 0;
    AliESDZDC *esdZDC = fESD->GetESDZDC();    
    Double_t  zNCEnergy = esdZDC->GetZDCN1Energy();
    Double_t  zPCEnergy = esdZDC->GetZDCP1Energy();
    Double_t  zNAEnergy = esdZDC->GetZDCN2Energy();
    Double_t  zPAEnergy = esdZDC->GetZDCP2Energy();
    // zdcA = (zNAEnergy>minEnergy || zPAEnergy>minEnergy);
    // zdcC = (zNCEnergy>minEnergy || zPCEnergy>minEnergy);
    zdcA = (zNAEnergy>minEnergy);
    zdcC = (zNCEnergy>minEnergy);
  }


  // Some macros for the online triggers
  Bool_t cMBS2A = c0sm2 && c0v0A;
  Bool_t cMBS2C = c0sm2 && c0v0C;
  Bool_t cMBAC  = c0v0A && c0v0C;
  




  // Plots requested by Jurgen on 18/11/2010 + later updates (including plots for the note)
  // FIXME: will skip everything else


  if(zdcTime) {
    GetHistoSPD1  ("ZDCTIME", "ZDC Time Cut")->Fill(outerLayerSPD);
    GetHistoTracks("ZDCTIME", "ZDC Time Cut")->Fill(ntracks);
    GetHistoV0M   ("ZDCTIME", "ZDC Time Cut")->Fill(multV0);
  }

  if(zdcTime && ntracks > 1) {
    GetHistoSPD1  ("ZDCTIME1TRACK", "ZDC Time Cut & 1 Track")->Fill(outerLayerSPD);
    GetHistoTracks("ZDCTIME1TRACK", "ZDC Time Cut & 1 Track")->Fill(ntracks);
    GetHistoV0M   ("ZDCTIME1TRACK", "ZDC Time Cut & 1 Track")->Fill(multV0);
  }

  // GetHistoSPD1  ("PhysSel", "All events after physics selection and Francesco's cut")->Fill(outerLayerSPD);
  // GetHistoTracks("PhysSel", "All events after physics selection and Francesco's cut")->Fill(ntracks);
  // GetHistoV0M   ("PhysSel", "All events after physics selection and Francesco's cut")->Fill(multV0);
  if(c0v0A && c0v0C) {
    GetHistoSPD1  ("V0AND", "V0A & V0C")->Fill(outerLayerSPD);
    GetHistoTracks("V0AND", "V0A & V0C")->Fill(ntracks);
    GetHistoV0M   ("V0AND", "V0A & V0C")->Fill(multV0);
  }
  if((c0v0A && !c0v0C) || (!c0v0A && c0v0C)) {
    GetHistoSPD1  ("V0ONLYONE", "(V0A & !V0C) || (!V0A & V0C)")->Fill(outerLayerSPD);
    GetHistoTracks("V0ONLYONE", "(V0A & !V0C) || (!V0A & V0C)")->Fill(ntracks);
    GetHistoV0M   ("V0ONLYONE", "(V0A & !V0C) || (!V0A & V0C)")->Fill(multV0);
  }
  if(zdcA && zdcC) {
    GetHistoSPD1  ("ZDCAND", "ZDCA & ZDCC")->Fill(outerLayerSPD);
    GetHistoTracks("ZDCAND", "ZDCA & ZDCC")->Fill(ntracks);
    GetHistoV0M   ("ZDCAND", "ZDCA & ZDCC")->Fill(multV0);
  }
  if((c0v0A && c0v0C) && !(zdcA && zdcC)) {
    GetHistoSPD1  ("V0ANDNOTZDCAND", "(V0A & V0C) & !(ZDCA & ZDCC)")->Fill(outerLayerSPD);
    GetHistoTracks("V0ANDNOTZDCAND", "(V0A & V0C) & !(ZDCA & ZDCC)")->Fill(ntracks);
    GetHistoV0M   ("V0ANDNOTZDCAND", "(V0A & V0C) & !(ZDCA & ZDCC)")->Fill(multV0);
  }
  if((c0v0A && c0v0C) && !zdcA && !zdcC) {
    GetHistoSPD1  ("V0ANDNOTZDC", "(V0A & V0C) & !ZDCA & !ZDCC)")->Fill(outerLayerSPD);
    GetHistoTracks("V0ANDNOTZDC", "(V0A & V0C) & !ZDCA & !ZDCC)")->Fill(ntracks);
    GetHistoV0M   ("V0ANDNOTZDC", "(V0A & V0C) & !ZDCA & !ZDCC)")->Fill(multV0);
}
  if((c0v0A && c0v0C) && ((!zdcA && zdcC) || (zdcA && !zdcC))) {
    GetHistoSPD1  ("V0ANDONEZDC", "(V0A & V0C) &  ((!ZDCA && ZDCC) || (ZDCA && !ZDCC)))")->Fill(outerLayerSPD);
    GetHistoTracks("V0ANDONEZDC", "(V0A & V0C) &  ((!ZDCA && ZDCC) || (ZDCA && !ZDCC)))")->Fill(ntracks);
    GetHistoV0M   ("V0ANDONEZDC", "(V0A & V0C) &  ((!ZDCA && ZDCC) || (ZDCA && !ZDCC)))")->Fill(multV0);
  }

  if(((c0v0A && !c0v0C) || (!c0v0A && c0v0C)) && (zdcA && zdcC)) {
    GetHistoSPD1  ("V0ONEZDCBOTH", "((V0A && !V0C) || (!V0A && V0C)) && (ZDCA && ZDCC)")->Fill(outerLayerSPD);
    GetHistoTracks("V0ONEZDCBOTH", "((V0A && !V0C) || (!V0A && V0C)) && (ZDCA && ZDCC)")->Fill(ntracks);
    GetHistoV0M   ("V0ONEZDCBOTH", "((V0A && !V0C) || (!V0A && V0C)) && (ZDCA && ZDCC)")->Fill(multV0);
  }

  if((c0v0A && c0v0C) && (zdcA && zdcC)) {
    GetHistoSPD1  ("V0ANDZDCAND", "(V0A & V0C) & (ZDCA & ZDCC)")->Fill(outerLayerSPD);
    GetHistoTracks("V0ANDZDCAND", "(V0A & V0C) & (ZDCA & ZDCC)")->Fill(ntracks);
    GetHistoV0M   ("V0ANDZDCAND", "(V0A & V0C) & (ZDCA & ZDCC)")->Fill(multV0);
  }

  if((c0v0A && zdcA && !zdcC && !c0v0C) || (c0v0C && zdcC && !zdcA && !c0v0A)) {
    GetHistoSPD1  ("OneSided", "(V0A & ZDCA & !ZDCC & !V0C) || (V0C & ZDCC & !ZDCA & !V0A)")->Fill(outerLayerSPD);
    GetHistoTracks("OneSided", "(V0A & ZDCA & !ZDCC & !V0C) || (V0C & ZDCC & !ZDCA & !V0A)")->Fill(ntracks);
    GetHistoV0M   ("OneSided", "(V0A & ZDCA & !ZDCC & !V0C) || (V0C & ZDCC & !ZDCA & !V0A)")->Fill(multV0);
  }

  // GetHistoCorrelationSPDTPCVz("All", "After physics selection and Francesco's cut")->Fill(vtxESDSPD->GetZ(),vtxESDTPC->GetZ());
  // if(cut1T0)   GetHistoCorrelationSPDTPCVz("Cut1T0", "T0 Cut 1")->Fill(vtxESDSPD->GetZ(),vtxESDTPC->GetZ());
  // if(cut2T0)   GetHistoCorrelationSPDTPCVz("Cut2T0", "T0 Cut 2")->Fill(vtxESDSPD->GetZ(),vtxESDTPC->GetZ());

  // GetHistoCorrelationContrTPCSPDCls("All", "After physics selection and Francesco's cut")->Fill(totalClusSPD,tpcContr);
  // if(cut1T0)   GetHistoCorrelationContrTPCSPDCls("Cut1T0", "T0 Cut 1")->Fill(totalClusSPD,tpcContr);
  // if(cut2T0)   GetHistoCorrelationContrTPCSPDCls("Cut2T0", "T0 Cut 2")->Fill(totalClusSPD,tpcContr);

  // GetHistoCorrelationTrackletsSPDCls("All", "After physics selection and Francesco's cut")->Fill(totalClusSPD,ntracklets);
  // if(cut1T0)   GetHistoCorrelationTrackletsSPDCls("Cut1T0", "T0 Cut 1")->Fill(totalClusSPD,ntracklets);
  // if(cut2T0)   GetHistoCorrelationTrackletsSPDCls("Cut2T0", "T0 Cut 2")->Fill(totalClusSPD,ntracklets);


  return; // FIXME


  // Reject background
  if (v0BG && fRejectBGWithV0) {
    cout << "Rejection BG" << endl;
    
    return;
  }
  // Fill global histos
  GetHistoTracklets("all","All events")->Fill(ntracklets);
  FillTriggerOverlaps("All", "All Events",vdArray);
  
  // Fill some combination of trigger classes
  Bool_t cmbs1aOnline = fESD->IsTriggerClassFired("CMBS1A-B-NOPF-ALL");
  Bool_t cmbs1cOnline = fESD->IsTriggerClassFired("CMBS1C-B-NOPF-ALL");
  Bool_t cmbacOnline  = fESD->IsTriggerClassFired("CMBAC-B-NOPF-ALL");
  
  Bool_t twoOutOfThree = kFALSE;
  if (cmbs1aOnline || cmbs1cOnline ||cmbacOnline) twoOutOfThree = kTRUE;
  if (twoOutOfThree) GetHistoTracklets("All_TwoOutOfThree" ,"Events 2-out-of-3 online" )->Fill(ntracklets);


  // loop over trigger classes in the event
  TObjArray * tokens = 0;
  if(fIsMC) {
    // in case of montecarlo I override the trigger class
    tokens = new TObjArray;
    tokens->SetOwner();
    //    tokens->Add(new TObjString("CINT1B-ABCE-NOPF-ALL")); 
    tokens->Add(new TObjString("MC")); 
  }
  else {  
    TString trgClasses = fESD->GetFiredTriggerClasses();
    tokens = trgClasses.Tokenize("  ");
  }
  TIterator  * iter = (TIterator*) tokens->MakeIterator();
  
  TString classes = fESD->GetFiredTriggerClasses();

  // if (classes.Contains("SMH")) {
  //    tokens->Print();
  //    cout << classes.Data() << endl;
  // }
  //  iter->Reset();
  //Int_t itoken = 0;
  TObjString * tok=0;
  while((tok = (TObjString*) iter->Next())){
    // clean up trigger name
    TString trg = tok->GetString();
    trg.Strip(TString::kTrailing, ' ');
    trg.Strip(TString::kLeading, ' ');

    fHistoSuffix = "_";
    fHistoSuffix += trg;

    //    cout << itoken++ << " " << trg.Data() << endl;
    // continue;
    //    if (trg.Contains("SMH")) cout << itoken++ << " " << trg.Data() << endl;
    
    // Fill tracklets 
    GetHistoTracklets("All" ,"Events no offline trigger" )->Fill(ntracklets);    


    // Fill histograms mismatchs
    // TODO: check mismatch trigger class 
    if(nFastOrOffline != nFastOrOnline) {
      GetHistoTracklets("mismatchingFastOr", "Events where fast or offline differs from fast-or online")->Fill(ntracklets);
    }
    
    if (c0v0A != v0AHW){
      GetHistoTracklets("mismatchingV0A", "Events where V0A offline differs from V0A online")->Fill(ntracklets);
    }
    
    if (c0v0C != v0CHW){
      GetHistoTracklets("mismatchingV0C", "Events where V0C offline differs from V0C online")->Fill(ntracklets);
    }    
    
    // Fill a tracklet histo for each trigger type
    if(c0sm1)  GetHistoTracklets("c0sm1" ,"Events were trigger c0sm1 fired" )->Fill(ntracklets);
    if(c0sm2)  GetHistoTracklets("c0sm2" ,"Events were trigger c0sm2 fired" )->Fill(ntracklets);
    if(c0sm3)  GetHistoTracklets("c0sm3" ,"Events were trigger c0sm3 fired" )->Fill(ntracklets);
    if(c0sm4)  GetHistoTracklets("c0sm4" ,"Events were trigger c0sm4 fired" )->Fill(ntracklets);
    if(c0sm5)  GetHistoTracklets("c0sm5" ,"Events were trigger c0sm5 fired" )->Fill(ntracklets);
    if(c0OM2)  GetHistoTracklets("c0OM2" ,"Events were trigger c0OM2 fired" )->Fill(ntracklets);
    if(c0OM3)  GetHistoTracklets("c0OM3" ,"Events were trigger c0OM3 fired" )->Fill(ntracklets);
    if(c0v0A)  GetHistoTracklets("c0v0A" ,"Events were trigger c0v0A fired" )->Fill(ntracklets);
    if(c0v0C)  GetHistoTracklets("c0v0C" ,"Events were trigger c0v0C fired" )->Fill(ntracklets);
    if(cMBS2A) GetHistoTracklets("cMBS2A","Events were trigger cMBS2A fired")->Fill(ntracklets);
    if(cMBS2C) GetHistoTracklets("cMBS2C","Events were trigger cMBS2C fired")->Fill(ntracklets);
    if(cMBAC ) GetHistoTracklets("cMBAC", "Events were trigger cMBAC  fired")->Fill(ntracklets);
    if(zdcA )  GetHistoTracklets("cZDCA", "Events were trigger cZDCA  fired")->Fill(ntracklets);
    if(zdcC )  GetHistoTracklets("cZDCC", "Events were trigger cZDCC  fired")->Fill(ntracklets);    
    if(!zdcA && !zdcC )  GetHistoTracklets("NoZDC", "Events were zdc trigger don't  fired")->Fill(ntracklets);
    if( (zdcA && !zdcC) || (!zdcA && zdcC) )  GetHistoTracklets("OneSideZDC", "Events were only 1 ZDC  fired")->Fill(ntracklets);
    if( (zdcA && zdcC) )  GetHistoTracklets("TwoSideZDC", "Events were both  ZDC  fired")->Fill(ntracklets);

    if(twoOutOfThree) {
      if(!zdcA && !zdcC )  GetHistoTracklets("NoZDC_TwoOutOfThree", "Events were zdc trigger don't  fired")->Fill(ntracklets);
      if( (zdcA && !zdcC) || (!zdcA && zdcC) )  GetHistoTracklets("OneSideZDC_TwoOutOfThree", "Events were only 1 ZDC  fired")->Fill(ntracklets);
      if( (zdcA && zdcC) )  GetHistoTracklets("TwoSideZDC_TwoOutOfThree", "Events were both  ZDC  fired")->Fill(ntracklets);
    }
    if(zdcBar) GetHistoTracklets("cZDCBar","Events were trigger cZDCB  fired")->Fill(ntracklets);
    //  if() GetHistoTracklets("","Events were trigger  fired");
    
    // Fill trigger overlaps
    FillTriggerOverlaps("All", "All Events in trigger class",vdArray);


    // Fill kinematic variables for peripheral events
    static AliESDtrackCuts * cuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(); // FIXME: change this ?
    if (ntracklets < fNTrackletsCutKine) {
      // 1. Loop on tracklets
      for(Int_t itracklet = 0; itracklet < ntracklets; itracklet++){      
	GetHistoEta("All", Form("Tracklet #eta distribution, ntracklets < %d",fNTrackletsCutKine))->Fill(mult->GetEta(itracklet));
      }
      // 2. loop on tracks
      TObjArray * goodTracks = cuts->GetAcceptedTracks(fESD);
      TIterator * iterTracks = goodTracks->MakeIterator();
      AliESDtrack * esdTrack = 0;
      while ((esdTrack = (AliESDtrack*) iterTracks->Next())) {
	GetHistoPt("All", Form("Tracklet p_{T} distribution, ntracklets < %d",fNTrackletsCutKine))->Fill(esdTrack->Pt());	
      }
      delete goodTracks;
    }
    

  }
  delete tokens;
    
    // if (fIsMC) {
    

  //   if (!fMCEvent) {
  //     AliError("No MC info found");
  //   } else {
      
  //     //loop on the MC event
  //     //      Int_t nMCTracks = fMCEvent->GetNumberOfTracks();
  //     Int_t offset    = fMCEvent->GetPrimaryOffset();
  //     Int_t nMCTracks = fMCEvent->GetNumberOfPrimaries()+offset;
  //     for (Int_t ipart=offset; ipart<nMCTracks; ipart++) { 
	
  // 	AliMCParticle *mcPart  = (AliMCParticle*)fMCEvent->GetTrack(ipart);
	
  // 	// We don't care about neutrals and non-physical primaries
  // 	if(mcPart->Charge() == 0) continue;

  // PHYSICAL PRIMARY
  // 	// Get MC vertex
    // 	TArrayF   vertex;
  // 	fMCEvent->GenEventHeader()->PrimaryVertex(vertex);
  // 	Float_t zv = vertex[2];
  // 	//	Float_t zv = vtxESD->GetZ();
  // 	// Fill generated histo
  // 	hTracks[AliAnalysisMultPbTrackHistoManager::kHistoGen]->Fill(mcPart->Pt(),mcPart->Eta(),zv);
	
  //     }
  //   }
  // }
  



}

void   AliAnalysisTaskTriggerStudy::Terminate(Option_t *){
  // terminate
  // Save output in a more friendly format
  fHistoList = dynamic_cast<AliHistoListWrapper*> (GetOutputData(1));
  if (!fHistoList){
    Printf("ERROR: fHistoList not available");
    return;
  }
  TFile * f = new TFile("trigger_study.root", "recreate");
  fHistoList->GetList()->Write();
  f->Close();

}

TH1 *   AliAnalysisTaskTriggerStudy::GetHistoTracklets(const char * name, const char * title){
  // Book histo of events vs ntracklets, if needed

  TString hname = "hTracklets_";
  hname+=name;  
  hname+=fHistoSuffix;
  TH1 * h = (TH1*) fHistoList->GetList()->FindObject(hname.Data());
  
  if(!h) {
    AliInfo(Form("Booking histo %s",hname.Data()));
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    // h = new TH1F (hname.Data(), title, 1000, -0.5, 10000-0.5);
    h = new TH1F (hname.Data(), title, 4000, -0.5, 4000-0.5);
    h->Sumw2();
    h->SetXTitle("ntracklets");
    fHistoList->GetList()->Add(h);
    TH1::AddDirectory(oldStatus);
  }
  return h;
}


TH1 *   AliAnalysisTaskTriggerStudy::GetHistoSPD1(const char * name, const char * title){
  // Book histo of events vs ntracklets, if needed

  TString hname = "hSPD1_";
  hname+=name;  
  hname+=fHistoSuffix;
  TH1 * h = (TH1*) fHistoList->GetList()->FindObject(hname.Data());
  
  // const Int_t nbin = 118;
  // Float_t bins[nbin] = {0};
  // bins[0] = 0;
  // for(Int_t ibin = 1; ibin <= nbin; ibin++){
  //   if (ibin < 100) 
  //     bins[ibin] = bins [ibin-1]+1; 
  //   else if (ibin < 1000) 
  //     bins[ibin] = bins [ibin-1]+100; 
  //   else if (ibin < 10000) 
  //     bins[ibin] = bins [ibin-1]+1000; 
  // }
  

  if(!h) {
    AliInfo(Form("Booking histo %s",hname.Data()));
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    //    h = new TH1F (hname.Data(), title, nbin, bins);
    //h = new TH1F (hname.Data(), title, 100, -0.5, 100-0.5);
    h = new TH1F (hname.Data(), title, 10000, -0.5, 10000-0.5);
    h->Sumw2();
    h->SetXTitle("SPD1");
    fHistoList->GetList()->Add(h);
    TH1::AddDirectory(oldStatus);
  }
  return h;
}

TH1 *   AliAnalysisTaskTriggerStudy::GetHistoV0M(const char * name, const char * title){
  // Book histo of events vs ntracklets, if needed

  TString hname = "hV0M_";
  hname+=name;  
  hname+=fHistoSuffix;
  TH1 * h = (TH1*) fHistoList->GetList()->FindObject(hname.Data());
  
  if(!h) {
    AliInfo(Form("Booking histo %s",hname.Data()));
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    // h = new TH1F (hname.Data(), title, 1000, -0.5, 10000-0.5);
    h = new TH1F (hname.Data(), title, 500, -0.5, 500-0.5);
    h->Sumw2();
    h->SetXTitle("V0M");
    fHistoList->GetList()->Add(h);
    TH1::AddDirectory(oldStatus);
  }
  return h;
}


TH1 *   AliAnalysisTaskTriggerStudy::GetHistoTracks(const char * name, const char * title){
  // Book histo of events vs ntracklets, if needed

  TString hname = "hTPCTracks_";
  hname+=name;  
  hname+=fHistoSuffix;
  TH1 * h = (TH1*) fHistoList->GetList()->FindObject(hname.Data());
  
  if(!h) {
    AliInfo(Form("Booking histo %s",hname.Data()));
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    // h = new TH1F (hname.Data(), title, 1000, -0.5, 10000-0.5);
    h = new TH1F (hname.Data(), title, 50, -0.5, 50-0.5);
    h->Sumw2();
    h->SetXTitle("ntracks");
    fHistoList->GetList()->Add(h);
    TH1::AddDirectory(oldStatus);
  }
  return h;
}


TH1 *   AliAnalysisTaskTriggerStudy::GetHistoEta(const char * name, const char * title){
  // Book histo of events vs ntracklets, if needed

  TString hname = "hEta_";
  hname+=name;  
  hname+=fHistoSuffix;
  TH1 * h = (TH1*) fHistoList->GetList()->FindObject(hname.Data());
  
  if(!h) {
    AliInfo(Form("Booking histo %s",hname.Data()));
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    // h = new TH1F (hname.Data(), title, 1000, -0.5, 10000-0.5);
    h = new TH1F (hname.Data(), title, 20, -2, 2);
    h->Sumw2();
    h->SetXTitle("#eta");
    fHistoList->GetList()->Add(h);
    TH1::AddDirectory(oldStatus);
  }
  return h;
}

TH1 *   AliAnalysisTaskTriggerStudy::GetHistoPt(const char * name, const char * title){
  // Book histo of pt distribution, if needed

  TString hname = "hPt_";
  hname+=name;  
  hname+=fHistoSuffix;
  TH1 * h = (TH1*) fHistoList->GetList()->FindObject(hname.Data());
  
  if(!h) {
    AliInfo(Form("Booking histo %s",hname.Data()));
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    // h = new TH1F (hname.Data(), title, 1000, -0.5, 10000-0.5);
    h = new TH1F (hname.Data(), title, 100, -0.5, 1-0.5);
    h->Sumw2();
    h->SetXTitle("p_{T}");
    fHistoList->GetList()->Add(h);
    TH1::AddDirectory(oldStatus);
  }
  return h;
}

TH1 *   AliAnalysisTaskTriggerStudy::GetHistoCorrelationSPDTPCVz(const char * name, const char * title){
  // Book histo of events vs ntracklets, if needed

  TString hname = "hTPCvsSPD_";
  hname+=name;  
  hname+=fHistoSuffix;
  TH1 * h = (TH1*) fHistoList->GetList()->FindObject(hname.Data());
  
  if(!h) {
    AliInfo(Form("Booking histo %s",hname.Data()));
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    // h = new TH1F (hname.Data(), title, 1000, -0.5, 10000-0.5);
    h = new TH2F (hname.Data(), title, 80, -20, 20,  80, -20, 20);
    h->Sumw2();
    h->SetXTitle("SPD Vz");
    h->SetYTitle("TPC Vz");
    fHistoList->GetList()->Add(h);
    TH1::AddDirectory(oldStatus);
  }
  return h;
}

TH1 *   AliAnalysisTaskTriggerStudy::GetHistoCorrelationContrTPCSPDCls(const char * name, const char * title){
  // Book histo of events vs ntracklets, if needed

  TString hname = "hContrTPCvsSPDCls_";
  hname+=name;  
  hname+=fHistoSuffix;
  TH1 * h = (TH1*) fHistoList->GetList()->FindObject(hname.Data());
  
  if(!h) {
    AliInfo(Form("Booking histo %s",hname.Data()));
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    // h = new TH1F (hname.Data(), title, 1000, -0.5, 10000-0.5);
    h = new TH2F (hname.Data(), title, 1000, 0, 9000,  1000, 0, 5000);
    h->Sumw2();
    h->SetXTitle("SPD Clusters");
    h->SetYTitle("TPC V Contributors");
    fHistoList->GetList()->Add(h);
    TH1::AddDirectory(oldStatus);
  }
  return h;
}
TH1 *   AliAnalysisTaskTriggerStudy::GetHistoCorrelationTrackletsSPDCls(const char * name, const char * title){
  // Book histo of events vs ntracklets, if needed

  TString hname = "hTrackletsvsSPDCls_";
  hname+=name;  
  hname+=fHistoSuffix;
  TH1 * h = (TH1*) fHistoList->GetList()->FindObject(hname.Data());
  
  if(!h) {
    AliInfo(Form("Booking histo %s",hname.Data()));
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    // h = new TH1F (hname.Data(), title, 1000, -0.5, 10000-0.5);
    h = new TH2F (hname.Data(), title, 1000, 0, 9000,  1000, 0, 5000);
    h->Sumw2();
    h->SetXTitle("SPD Clusters");
    h->SetYTitle("N tracklets");
    fHistoList->GetList()->Add(h);
    TH1::AddDirectory(oldStatus);
  }
  return h;
}


void AliAnalysisTaskTriggerStudy::FillTriggerOverlaps (const char * name, const char * title, Bool_t * vdArray){
  //Fills a histo with the different trigger statistics in a venn like diagramm. Books it if needed.

  // Get or book histo
  TString hname = "hTrigStat_";
  hname+=name;  
  hname+=fHistoSuffix;
  TH1 * h = (TH1*) fHistoList->GetList()->FindObject(hname.Data());
  
  if(!h) {
    AliInfo(Form("Booking histo %s",hname.Data()));
    Bool_t oldStatus = TH1::AddDirectoryStatus();
    TH1::AddDirectory(kFALSE);
    Int_t nbins = 0;
    for(Int_t ientry = 0; ientry < kNVDEntries; ientry++){
      nbins = nbins | (1<<ientry);
    }
    //    cout << "NBINS " << nbins << endl;
    nbins = nbins+1;
    h = new TH1I (hname, title, nbins, -0.5, nbins-0.5);
    fHistoList->GetList()->Add(h);
    TH1::AddDirectory(oldStatus);
  
    // we look at the combinations of n triggers
    // We set a bit for each trigger to fill the diagram
    // This is much simpler and faster than any recursive function
    h->GetXaxis()->SetBinLabel(1,"NONE"); 
    for(Int_t ibin = 1; ibin < nbins; ibin++){
      TString binname = "";
      Bool_t first = kTRUE;
      for(Int_t ivdentry = 0; ivdentry < kNVDEntries; ivdentry++){
	if (ibin & (1<<ivdentry)) {
	  if(!first) binname += " & ";
	  binname += kVDNames[ivdentry];
	  first=kFALSE;
	}
      }
      h->GetXaxis()->SetBinLabel(ibin+1,binname.Data());
    }
    
  }

  UInt_t mask = 0;
  for(Int_t ivdentry = 0; ivdentry < kNVDEntries; ivdentry++){
    if(vdArray[ivdentry]) {
      mask  = mask | (1<<ivdentry);
      //      cout << " 1 "   ;
    } //else cout << " 0 ";
  }
  //  cout << hex << " = " << mask << endl;
  
  h->Fill(mask);

}
