
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>

#include "AliExternalTrackParam.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "TLorentzVector.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDv0.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskCorrPbPbMC.h"

//For MC event
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"

#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1I.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TMath.h"
#include "TList.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDatabasePDG.h"




using namespace std;
using std::cout;
using std::endl;

class AliAnalysisTaskCorrPbPbMC;
ClassImp(AliAnalysisTaskCorrPbPbMC)

//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskCorrPbPbMC::AliAnalysisTaskCorrPbPbMC():
  AliAnalysisTaskSE(),
  fAODeventCuts(),
  fESDevent(0),
  fAODevent(0),
  fMCevent(0),
  fMCstack(0),
  fInputEvent(0),
  fPIDResponse(0),
  fUtils(0),
  fOutputList(0),
  fQAList(0),
  fTreeEvent(0),
  fESDtrackCuts(0),
  fESDtrackCuts_primary(0),
  //fTrigger(AliVEvent::kINT7),
  fMultLow(0),
  fMultHigh(100),
  hNumberOfEvents(0),
  hNumberOfKaonPlus(0),
  hNumberOfKaonMinus(0),
  hNumberOfPionPlus(0),
  hNumberOfPionMinus(0),
  hNumberOfProtonPlus(0),
  hNumberOfProtonMinus(0),
  fTreeVariableCentrality(0),
  fNoGenKaonPlus_ptmax2(0),
  fNoGenKaonMinus_ptmax2(0),
  fNoGenKaonPlus_ptmax3(0),
  fNoGenKaonMinus_ptmax3(0),
  fNoGenPionPlus_ptmax2(0),
  fNoGenPionMinus_ptmax2(0),
  fNoGenPionPlus_ptmax3(0),
  fNoGenPionMinus_ptmax3(0),
  fNoGenProtonPlus_ptmax2(0),
  fNoGenProtonMinus_ptmax2(0),
  fNoGenProtonPlus_ptmax3(0),
  fNoGenProtonMinus_ptmax3(0),
  fNoKaonPlus_ptmax2(0),
  fNoKaonMinus_ptmax2(0),
  fNoKaonPlus_ptmax3(0),
  fNoKaonMinus_ptmax3(0),
  fNoPionPlus_ptmax2(0),
  fNoPionMinus_ptmax2(0),
  fNoPionPlus_ptmax3(0),
  fNoPionMinus_ptmax3(0),
  fNoProtonPlus_ptmax2(0),
  fNoProtonMinus_ptmax2(0),
  fNoProtonPlus_ptmax3(0),
  fNoProtonMinus_ptmax3(0),
  hist_KaonPlusWithoutPdg(0),
  hist_KaonPlusWithPdg(0),
  hist_KaonMinusWithoutPdg(0),
  hist_KaonMinusWithPdg(0),
  hist_PionPlusWithoutPdg(0),
  hist_PionPlusWithPdg(0),
  hist_PionMinusWithoutPdg(0),
  hist_PionMinusWithPdg(0),
  hist_ProtonPlusWithoutPdg(0),
  hist_ProtonPlusWithPdg(0),
  hist_ProtonMinusWithoutPdg(0),
  hist_ProtonMinusWithPdg(0),
  hist_GenKaonPlus(0),
  hist_GenKaonMinus(0),
  hist_GenProtonPlus(0),
  hist_GenProtonMinus(0),
  hist_GenPionPlus(0),
  hist_GenPionMinus(0)
{}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskCorrPbPbMC::AliAnalysisTaskCorrPbPbMC(const char *name):
  AliAnalysisTaskSE(name),
  fAODeventCuts(),
  fESDevent(0),
  fAODevent(0),
  fMCevent(0),
  fMCstack(0),
  fInputEvent(0),
  fPIDResponse(0),
  fUtils(0),
  fOutputList(0),
  fQAList(0),
  fTreeEvent(0),
  fESDtrackCuts(0),
  fESDtrackCuts_primary(0),
  //fTrigger(AliVEvent::kINT7),
  fMultLow(0),
  fMultHigh(100),
  hNumberOfEvents(0),
  hNumberOfKaonPlus(0),
  hNumberOfKaonMinus(0),
  hNumberOfPionPlus(0),
  hNumberOfPionMinus(0),
  hNumberOfProtonPlus(0),
  hNumberOfProtonMinus(0),
  fTreeVariableCentrality(0),
  fNoGenKaonPlus_ptmax2(0),
  fNoGenKaonMinus_ptmax2(0),
  fNoGenKaonPlus_ptmax3(0),
  fNoGenKaonMinus_ptmax3(0),
  fNoGenPionPlus_ptmax2(0),
  fNoGenPionMinus_ptmax2(0),
  fNoGenPionPlus_ptmax3(0),
  fNoGenPionMinus_ptmax3(0),
  fNoGenProtonPlus_ptmax2(0),
  fNoGenProtonMinus_ptmax2(0),
  fNoGenProtonPlus_ptmax3(0),
  fNoGenProtonMinus_ptmax3(0),
  fNoKaonPlus_ptmax2(0),
  fNoKaonMinus_ptmax2(0),
  fNoKaonPlus_ptmax3(0),
  fNoKaonMinus_ptmax3(0),
  fNoPionPlus_ptmax2(0),
  fNoPionMinus_ptmax2(0),
  fNoPionPlus_ptmax3(0),
  fNoPionMinus_ptmax3(0),
  fNoProtonPlus_ptmax2(0),
  fNoProtonMinus_ptmax2(0),
  fNoProtonPlus_ptmax3(0),
  fNoProtonMinus_ptmax3(0),
  hist_KaonPlusWithoutPdg(0),
  hist_KaonPlusWithPdg(0),
  hist_KaonMinusWithoutPdg(0),
  hist_KaonMinusWithPdg(0),
  hist_PionPlusWithoutPdg(0),
  hist_PionPlusWithPdg(0),
  hist_PionMinusWithoutPdg(0),
  hist_PionMinusWithPdg(0),
  hist_ProtonPlusWithoutPdg(0),
  hist_ProtonPlusWithPdg(0),
  hist_ProtonMinusWithoutPdg(0),
  hist_ProtonMinusWithPdg(0),
  hist_GenKaonPlus(0),
  hist_GenKaonMinus(0),
  hist_GenProtonPlus(0),
  hist_GenProtonMinus(0),
  hist_GenPionPlus(0),
  hist_GenPionMinus(0)
{
  fUtils = new AliAnalysisUtils();
  DefineInput (0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TTree::Class());
}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskCorrPbPbMC::~AliAnalysisTaskCorrPbPbMC()  {

  if (fOutputList){
    delete fOutputList;
    fOutputList = 0x0;
  }
  if (fQAList){
    delete fQAList;
    fQAList = 0x0;
  }
  if (fTreeEvent){
    delete fTreeEvent;
    fTreeEvent = 0x0;
  }
  if(fESDtrackCuts) {
    delete fESDtrackCuts;
    fESDtrackCuts=0x0;
  }


}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskCorrPbPbMC::UserCreateOutputObjects()  {
    
    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner(kTRUE);
    fQAList     -> SetOwner(kTRUE);

    OpenFile(1);
    OpenFile(2);
    OpenFile(3);

    /*
    //QA Plots of Event Selection
    fESDeventCuts.AddQAplotsToList(fQAList,kTRUE);
    fESDeventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
    */
    
    //Event Counter
    hNumberOfEvents = new TH1D ("hNumberOfEvents","",20,0,20);
    // hNumberOfEvents -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    

    //Number of Kaon finally getting selected with specified cuts event wise
    hNumberOfKaonPlus     = new TH1D ("hNumberOfKaonPlus","",3000,0,3000);
    hNumberOfKaonMinus = new TH1D ("hNumberOfKaonMinus","",3000,0,3000);
    // hNumberOfKaonPlus     -> Sumw2();
    // hNumberOfKaonMinus -> Sumw2();
    fOutputList -> Add(hNumberOfKaonPlus);
    fOutputList -> Add(hNumberOfKaonMinus);

    //Number of Pion finally getting selected with specified cuts event wise
    hNumberOfPionPlus     = new TH1D ("hNumberOfPionPlus","",3000,0,3000);
    hNumberOfPionMinus = new TH1D ("hNumberOfPionMinus","",3000,0,3000);
    // hNumberOfPionPlus     -> Sumw2();
    // hNumberOfPionMinus -> Sumw2();
    fOutputList -> Add(hNumberOfPionPlus);
    fOutputList -> Add(hNumberOfPionMinus);

    //Number of Proton finally getting selected with specified cuts event wise
    hNumberOfProtonPlus     = new TH1D ("hNumberOfProtonPlus","",3000,0,3000);
    hNumberOfProtonMinus = new TH1D ("hNumberOfProtonMinus","",3000,0,3000);
    // hNumberOfProtonPlus     -> Sumw2();
    // hNumberOfProtonMinus -> Sumw2();
    fOutputList -> Add(hNumberOfProtonPlus);
    fOutputList -> Add(hNumberOfProtonMinus);

    
    //distributions for calculating Kaon purity
    hist_KaonPlusWithoutPdg = new TH1F("hist_KaonPlusWithoutPdg"," hist_KaonPlusWithoutPdg",300,0,3.0);
    hist_KaonPlusWithPdg = new TH1F("hist_KaonPlusWithPdg"," hist_KaonPlusWithPdg",300,0,3.0);
    hist_KaonMinusWithoutPdg = new TH1F("hist_KaonMinusWithoutPdg"," hist_KaonMinusWithoutPdg",300,0,3.0);
    hist_KaonMinusWithPdg = new TH1F("hist_KaonMinusWithPdg"," hist_KaonMinusWithPdg",300,0,3.0);
    fOutputList->Add(hist_KaonPlusWithoutPdg);
    fOutputList->Add(hist_KaonPlusWithPdg);
    fOutputList->Add(hist_KaonMinusWithoutPdg);
    fOutputList->Add(hist_KaonMinusWithPdg);

    //distributions for calculating Pion purity
    hist_PionPlusWithoutPdg = new TH1F("hist_PionPlusWithoutPdg"," hist_PionPlusWithoutPdg",300,0,3.0);
    hist_PionPlusWithPdg = new TH1F("hist_PionPlusWithPdg"," hist_PionPlusWithPdg",300,0,3.0);
    hist_PionMinusWithoutPdg = new TH1F("hist_PionMinusWithoutPdg"," hist_PionMinusWithoutPdg",300,0,3.0);
    hist_PionMinusWithPdg = new TH1F("hist_PionMinusWithPdg"," hist_PionMinusWithPdg",300,0,3.0);
    fOutputList->Add(hist_PionPlusWithoutPdg);
    fOutputList->Add(hist_PionPlusWithPdg);
    fOutputList->Add(hist_PionMinusWithoutPdg);
    fOutputList->Add(hist_PionMinusWithPdg);

    //distributions for calculating Proton purity
    hist_ProtonPlusWithoutPdg = new TH1F("hist_ProtonPlusWithoutPdg"," hist_ProtonPlusWithoutPdg",300,0,3.0);
    hist_ProtonPlusWithPdg = new TH1F("hist_ProtonPlusWithPdg"," hist_ProtonPlusWithPdg",300,0,3.0);
    hist_ProtonMinusWithoutPdg = new TH1F("hist_ProtonMinusWithoutPdg"," hist_ProtonMinusWithoutPdg",300,0,3.0);
    hist_ProtonMinusWithPdg = new TH1F("hist_ProtonMinusWithPdg"," hist_ProtonMinusWithPdg",300,0,3.0);
    fOutputList->Add(hist_ProtonPlusWithoutPdg);
    fOutputList->Add(hist_ProtonPlusWithPdg);
    fOutputList->Add(hist_ProtonMinusWithoutPdg);
    fOutputList->Add(hist_ProtonMinusWithPdg);

    //Generated distributions: for efficiency
    hist_GenKaonPlus = new TH1F("hist_GenKaonPlus"," hist_GenKaonPlus",300,0,3.0);
    hist_GenKaonMinus = new TH1F("hist_GenKaonMinus"," hist_GenKaonMinus",300,0,3.0);
    hist_GenPionPlus = new TH1F("hist_GenPionPlus"," hist_GenPionPlus",300,0,3.0);
    hist_GenPionMinus = new TH1F("hist_GenPionMinus"," hist_GenPionMinus",300,0,3.0);
    hist_GenProtonPlus = new TH1F("hist_GenProtonPlus"," hist_GenProtonPlus",300,0,3.0);
    hist_GenProtonMinus = new TH1F("hist_GenProtonMinus"," hist_GenProtonMinus",300,0,3.0);
    fOutputList->Add(hist_GenKaonPlus);
    fOutputList->Add(hist_GenKaonMinus);
    fOutputList->Add(hist_GenPionPlus);
    fOutputList->Add(hist_GenPionMinus);
    fOutputList->Add(hist_GenProtonPlus);
    fOutputList->Add(hist_GenProtonMinus);
  
    //TTree object to store variables
    fTreeEvent = new TTree("fTreeEvent","Event Tree");
    fTreeEvent->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
    //Generated
    fTreeEvent->Branch("fNoGenKaonPlus_ptmax2", &fNoGenKaonPlus_ptmax2, "fNoGenKaonPlus_ptmax2/F");
    fTreeEvent->Branch("fNoGenKaonMinus_ptmax2", &fNoGenKaonMinus_ptmax2, "fNoGenKaonMinus_ptmax2/F");
    fTreeEvent->Branch("fNoGenPionPlus_ptmax2", &fNoGenPionPlus_ptmax2, "fNoGenPionPlus_ptmax2/F");
    fTreeEvent->Branch("fNoGenPionMinus_ptmax2", &fNoGenPionMinus_ptmax2, "fNoGenPionMinus_ptmax2/F");
    fTreeEvent->Branch("fNoGenProtonPlus_ptmax2", &fNoGenProtonPlus_ptmax2, "fNoGenProtonPlus_ptmax2/F");
    fTreeEvent->Branch("fNoGenProtonMinus_ptmax2", &fNoGenProtonMinus_ptmax2, "fNoGenProtonMinus_ptmax2/F");

    fTreeEvent->Branch("fNoGenKaonPlus_ptmax3", &fNoGenKaonPlus_ptmax3, "fNoGenKaonPlus_ptmax3/F");
    fTreeEvent->Branch("fNoGenKaonMinus_ptmax3", &fNoGenKaonMinus_ptmax3, "fNoGenKaonMinus_ptmax3/F");
    fTreeEvent->Branch("fNoGenPionPlus_ptmax3", &fNoGenPionPlus_ptmax3, "fNoGenPionPlus_ptmax3/F");
    fTreeEvent->Branch("fNoGenPionMinus_ptmax3", &fNoGenPionMinus_ptmax3, "fNoGenPionMinus_ptmax3/F");
    fTreeEvent->Branch("fNoGenProtonPlus_ptmax3", &fNoGenProtonPlus_ptmax3, "fNoGenProtonPlus_ptmax3/F");
    fTreeEvent->Branch("fNoGenProtonMinus_ptmax3", &fNoGenProtonMinus_ptmax3, "fNoGenProtonMinus_ptmax3/F");
   
    //reconstructed
    fTreeEvent->Branch("fNoKaonPlus_ptmax2", &fNoKaonPlus_ptmax2, "fNoKaonPlus_ptmax2/F");
    fTreeEvent->Branch("fNoKaonMinus_ptmax2", &fNoKaonMinus_ptmax2, "fNoKaonMinus_ptmax2/F");
    fTreeEvent->Branch("fNoPionPlus_ptmax2", &fNoPionPlus_ptmax2, "fNoPionPlus_ptmax2/F");
    fTreeEvent->Branch("fNoPionMinus_ptmax2", &fNoPionMinus_ptmax2, "fNoPionMinus_ptmax2/F");
    fTreeEvent->Branch("fNoProtonPlus_ptmax2", &fNoProtonPlus_ptmax2, "fNoProtonPlus_ptmax2/F");
    fTreeEvent->Branch("fNoProtonMinus_ptmax2", &fNoProtonMinus_ptmax2, "fNoProtonMinus_ptmax2/F");

    fTreeEvent->Branch("fNoKaonPlus_ptmax3", &fNoKaonPlus_ptmax3, "fNoKaonPlus_ptmax3/F");
    fTreeEvent->Branch("fNoKaonMinus_ptmax3", &fNoKaonMinus_ptmax3, "fNoKaonMinus_ptmax3/F");
    fTreeEvent->Branch("fNoPionPlus_ptmax3", &fNoPionPlus_ptmax3, "fNoPionPlus_ptmax3/F");
    fTreeEvent->Branch("fNoPionMinus_ptmax3", &fNoPionMinus_ptmax3, "fNoPionMinus_ptmax3/F");
    fTreeEvent->Branch("fNoProtonPlus_ptmax3", &fNoProtonPlus_ptmax3, "fNoProtonPlus_ptmax3/F");
    fTreeEvent->Branch("fNoProtonMinus_ptmax3", &fNoProtonMinus_ptmax3, "fNoProtonMinus_ptmax3/F");
   
    /*
    
    //Track Cuts Objects
    if(!fESDtrackCuts )
      {
	fESDtrackCuts = new AliESDtrackCuts();
	fESDtrackCuts->SetMinNClustersTPC(70);
	fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
	fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
	fESDtrackCuts->SetRequireTPCRefit(kTRUE);
	fESDtrackCuts->SetPtRange(0.15,1e10);
	fESDtrackCuts->SetEtaRange(-0.8,0.8);
      }
    if(!fESDtrackCuts_primary)
      {
	fESDtrackCuts_primary = new AliESDtrackCuts();
	fESDtrackCuts_primary = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE, 1);
	fESDtrackCuts_primary->SetPtRange(0.3,5.0);
	fESDtrackCuts_primary->SetEtaRange(-0.8,+0.8);
      }
    */
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeEvent);
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskCorrPbPbMC::UserExec(Option_t *)  {
  
    //Get Input Event
    if ( !GetEvent ()) return;
    cout<<"*********************** Found AOD event !!! ******************************"<<endl;

    
    //Get multiplicity percentile
    Float_t lV0M;
    AliMultSelection *MultSelection = (AliMultSelection*) fAODevent -> FindListObject("MultSelection");
    if( !MultSelection)
      {
	return;
      }
    else
      {
	lV0M = MultSelection->GetMultiplicityPercentile("V0M");
	cout<<"V0M: "<<lV0M<<endl;
      }

    fMCevent = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMCevent) {
      Printf("ERROR: Could not retrieve MC event \n");
      cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeEvent);
      return;
    }
  
    // fMCstack = fMCevent->Stack();
    // if (!fMCstack) {
    //   Printf("ERROR: Could not retrieve MC stack \n");
    //   cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    //   PostData(1, fOutputList);
    //   PostData(2, fQAList);
    //   PostData(3, fTreeEvent);
    //   return;
    // }

    cout<<"*********************** Found MC event !!! ******************************"<<endl;


    //"fTreeEvent" Tree Variable
    fTreeVariableCentrality = lV0M;

    
    //Initialize generated number 0f K+, K-, pi+, pi-, p and p-bar per event
    Int_t no_GenKaonPlus_perevent = 0;
    Int_t no_GenKaonMinus_perevent = 0;
    Int_t no_GenProtonPlus_perevent = 0;
    Int_t no_GenProtonMinus_perevent = 0;
    Int_t no_GenPionPlus_perevent = 0;
    Int_t no_GenPionMinus_perevent = 0;

    Int_t no_GenKaonPlus_perevent_ptmax2 = 0;
    Int_t no_GenKaonMinus_perevent_ptmax2 = 0;
    Int_t no_GenProtonPlus_perevent_ptmax2 = 0;
    Int_t no_GenProtonMinus_perevent_ptmax2 = 0;
    Int_t no_GenPionPlus_perevent_ptmax2 = 0;
    Int_t no_GenPionMinus_perevent_ptmax2 = 0;

    
    //Loop on generated MC tracks
    Int_t noGenMCtracks = fMCevent->GetNumberOfTracks();
    cout<<"No of generated MC tracks: "<<noGenMCtracks<<endl;
    for(Int_t itr_mcgen=0; itr_mcgen < noGenMCtracks; itr_mcgen++)
      {
	AliAODMCParticle *mcGenTrack = (AliAODMCParticle*) fMCevent->GetTrack(itr_mcgen);
	if (!mcGenTrack)
	  {
	    cout<<"Could not find track in MC generated loop !!!"<<endl;
	    continue;
	  }
	//cout<<"Found generated MC track !!"<<endl;

	//Out of bunch pileup event removal
	if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(itr_mcgen,fMCevent))
	  {
	    cout<<"Track belongs to out of bunch pileup events. Removed !!!!"<<endl;
	    continue;
	  }
      

	if(!mcGenTrack->IsPhysicalPrimary())
	  continue;

	Double_t TrackGen_eta = mcGenTrack->Eta();
	Double_t TrackGen_pt = mcGenTrack->Pt();
	Double_t TrackGen_charge = mcGenTrack->Charge();
	Double_t TrackGen_PID = mcGenTrack->GetPdgCode();

	if (!TrackGen_PID) continue;
	if (TrackGen_pt < 0.2) continue;
	if (TrackGen_pt > 3.0) continue;
	if (TMath::Abs(TrackGen_eta) > 0.8) continue;

	//Kaon
	if (TMath::Abs(TrackGen_PID) == 321)
	  {
	    if(TrackGen_charge > 0)  //K+
	      {
		no_GenKaonPlus_perevent += 1;
		if(TrackGen_pt < 2.0)
		  no_GenKaonPlus_perevent_ptmax2 += 1;
		hist_GenKaonPlus->Fill(TrackGen_pt);
	      }
	    if(TrackGen_charge < 0)  //K-
	      {
		no_GenKaonMinus_perevent += 1;
		if(TrackGen_pt < 2.0)
		  no_GenKaonMinus_perevent_ptmax2 += 1;
		hist_GenKaonMinus->Fill(TrackGen_pt);
	      }
	  }

	//Proton
	if (TMath::Abs(TrackGen_PID) == 2212)
	  {
	    if(TrackGen_charge > 0)  //p
	      {
		no_GenProtonPlus_perevent += 1;
		if(TrackGen_pt < 2.0)
		  no_GenProtonPlus_perevent_ptmax2 += 1;
		hist_GenProtonPlus->Fill(TrackGen_pt);
	      }
	    if(TrackGen_charge < 0)  //pbar
	      {
		no_GenProtonMinus_perevent += 1;
		if(TrackGen_pt < 2.0)
		  no_GenProtonMinus_perevent_ptmax2 += 1;
		hist_GenProtonMinus->Fill(TrackGen_pt);
	      }
	  }

	//Pion
	if (TMath::Abs(TrackGen_PID) == 211)
	  {
	    if(TrackGen_charge > 0)  //pi+
	      {
		no_GenPionPlus_perevent += 1;
		if(TrackGen_pt < 2.0)
		  no_GenPionPlus_perevent_ptmax2 += 1;
		hist_GenPionPlus->Fill(TrackGen_pt);
	      }
	    if(TrackGen_charge < 0)  //pi-
	      {
		no_GenPionMinus_perevent += 1;
		if(TrackGen_pt < 2.0)
		  no_GenPionMinus_perevent_ptmax2 += 1;
		hist_GenPionMinus->Fill(TrackGen_pt);
	      }
	  }
      }
    //end generated track loop

    //Initialize number 0f K+, K-, pi+, pi-, p and p-bar per event
    Int_t no_KaonPlus_perevent = 0;
    Int_t no_KaonMinus_perevent = 0;
    Int_t no_ProtonPlus_perevent = 0;
    Int_t no_ProtonMinus_perevent = 0;
    Int_t no_PionPlus_perevent = 0;
    Int_t no_PionMinus_perevent = 0;

    Int_t no_KaonPlus_perevent_ptmax2 = 0;
    Int_t no_KaonMinus_perevent_ptmax2 = 0;
    Int_t no_ProtonPlus_perevent_ptmax2 = 0;
    Int_t no_ProtonMinus_perevent_ptmax2 = 0;
    Int_t no_PionPlus_perevent_ptmax2 = 0;
    Int_t no_PionMinus_perevent_ptmax2 = 0;


    //Loop on reconstructed MC tracks
    Int_t noRecMCtracks = fAODevent->GetNumberOfTracks();
    cout<<"No of reconstructed MC tracks: "<<noRecMCtracks<<endl;
    for(Int_t itr_mcrec=0; itr_mcrec < noRecMCtracks; itr_mcrec++)
      {
	
	AliVTrack   *track = (AliVTrack*)fAODevent->GetTrack(itr_mcrec);
	if(!track)      continue;
	AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(track);
	if(!aodtrack)      continue;

	Int_t aodtrk_label = TMath::Abs(aodtrack->GetLabel());
	AliAODMCParticle *mcRecTrack = (AliAODMCParticle*) fMCevent->GetTrack(aodtrk_label);
	if (!mcRecTrack)
	  {
	    cout<<"Could not find MC generted track for the AODtrack Label !!!"<<endl;
	    continue;
	  }
	cout<<"Found reconstructed MC track (matched with generated) !!"<<endl;

	if(!mcRecTrack->IsPhysicalPrimary())
	  continue;

	if(!aodtrack->TestFilterBit(96))  continue;

	Double_t Track_pt = mcRecTrack->Pt();
	Double_t Track_charge = mcRecTrack->Charge();
	Double_t Track_PID = mcRecTrack->GetPdgCode();
	Double_t Track_eta = mcRecTrack->Eta();

	if (Track_pt < 0.2) continue;
	if (Track_pt > 3.0) continue;
	if (TMath::Abs(Track_eta) > 0.8) continue;

	Bool_t IsKaon = KaonSelector(track);
	Bool_t IsPion = PionSelector(track);
	Bool_t IsProton = ProtonSelector(track);
	if (!IsKaon && !IsPion && !IsProton) continue;

	
	Int_t flag = 0;
	if(IsKaon) flag+=1;
	if(IsPion) flag+=1;
	if(IsProton) flag+=1;
	cout<<"Particle identified as more than on PID: flag= "<<flag<<endl;
	if(flag>1) continue;
	
	

	if (Track_charge > 0 && IsKaon)   //K+
	  {
	    no_KaonPlus_perevent += 1;
	    if(Track_pt < 2.0)
	      no_KaonPlus_perevent_ptmax2 += 1;

	    hist_KaonPlusWithoutPdg->Fill(Track_pt);
	    if (TMath::Abs(Track_PID) == 321)
	      hist_KaonPlusWithPdg->Fill(Track_pt);
	  }
	  

	if (Track_charge < 0 && IsKaon)   //K-
	  {
	    no_KaonMinus_perevent += 1;
	    if(Track_pt < 2.0)
	      no_KaonMinus_perevent_ptmax2 += 1;

	    hist_KaonMinusWithoutPdg->Fill(Track_pt);
	    if (TMath::Abs(Track_PID) == 321)
	      hist_KaonMinusWithPdg->Fill(Track_pt);
	  }

	if (Track_charge > 0 && IsProton)   //proton
	  {
	    no_ProtonPlus_perevent += 1;
	    if(Track_pt < 2.0)
	      no_ProtonPlus_perevent_ptmax2 += 1;

	    hist_ProtonPlusWithoutPdg->Fill(Track_pt);
	    if (TMath::Abs(Track_PID) == 2212)
	      hist_ProtonPlusWithPdg->Fill(Track_pt);
	  }
	  

	if (Track_charge < 0 && IsProton)   //anti-proton
	  {
	    no_ProtonMinus_perevent += 1;
	    if(Track_pt < 2.0)
	      no_ProtonMinus_perevent_ptmax2 += 1;

	    hist_ProtonMinusWithoutPdg->Fill(Track_pt);
	    if (TMath::Abs(Track_PID) == 2212)
	      hist_ProtonMinusWithPdg->Fill(Track_pt);
	  }

	if (Track_charge > 0 && IsPion)   //pi+
	  {
	    no_PionPlus_perevent += 1;
	    if(Track_pt < 2.0)
	      no_PionPlus_perevent_ptmax2 += 1;

	    hist_PionPlusWithoutPdg->Fill(Track_pt);
	    if (TMath::Abs(Track_PID) == 211)
	      hist_PionPlusWithPdg->Fill(Track_pt);
	  }
	  

	if (Track_charge < 0 && IsPion)   //pi-
	  {
	    no_PionMinus_perevent += 1;
	    if(Track_pt < 2.0)
	      no_PionMinus_perevent_ptmax2 += 1;

	    hist_PionMinusWithoutPdg->Fill(Track_pt);
	    if (TMath::Abs(Track_PID) == 211)
	      hist_PionMinusWithPdg->Fill(Track_pt);
	  }
     }
    //end reconstructed track loop
    

    
    //Kaon
    hNumberOfKaonPlus->Fill(no_KaonPlus_perevent);
    hNumberOfKaonMinus->Fill(no_KaonMinus_perevent);
    //Pion
    hNumberOfPionPlus->Fill(no_PionPlus_perevent);
    hNumberOfPionMinus->Fill(no_PionMinus_perevent);
    //Proton
    hNumberOfProtonPlus->Fill(no_ProtonPlus_perevent);
    hNumberOfProtonMinus->Fill(no_ProtonMinus_perevent);

    //++++++++++++++++++++++++++++++++++++++++++++
    //Generated
    //Kaon
    fNoGenKaonPlus_ptmax2 = no_GenKaonPlus_perevent_ptmax2;
    fNoGenKaonMinus_ptmax2 = no_GenKaonMinus_perevent_ptmax2;
    fNoGenKaonPlus_ptmax3 = no_GenKaonPlus_perevent;
    fNoGenKaonMinus_ptmax3 = no_GenKaonMinus_perevent;
    //Pion
    fNoGenPionPlus_ptmax2 = no_GenPionPlus_perevent_ptmax2;
    fNoGenPionMinus_ptmax2 = no_GenPionMinus_perevent_ptmax2;
    fNoGenPionPlus_ptmax3 = no_GenPionPlus_perevent;
    fNoGenPionMinus_ptmax3 = no_GenPionMinus_perevent;
    //Proton
    fNoGenProtonPlus_ptmax2 = no_GenProtonPlus_perevent_ptmax2;
    fNoGenProtonMinus_ptmax2 = no_GenProtonMinus_perevent_ptmax2;
    fNoGenProtonPlus_ptmax3 = no_GenProtonPlus_perevent;
    fNoGenProtonMinus_ptmax3 = no_GenProtonMinus_perevent;
    //++++++++++++++++++++++++++++++++++++++++++++
    //Reconstructed
    //Kaon
    fNoKaonPlus_ptmax2 = no_KaonPlus_perevent_ptmax2;
    fNoKaonMinus_ptmax2 = no_KaonMinus_perevent_ptmax2;
    fNoKaonPlus_ptmax3 = no_KaonPlus_perevent;
    fNoKaonMinus_ptmax3 = no_KaonMinus_perevent;
    //Pion
    fNoPionPlus_ptmax2 = no_PionPlus_perevent_ptmax2;
    fNoPionMinus_ptmax2 = no_PionMinus_perevent_ptmax2;
    fNoPionPlus_ptmax3 = no_PionPlus_perevent;
    fNoPionMinus_ptmax3 = no_PionMinus_perevent;
    //Proton
    fNoProtonPlus_ptmax2 = no_ProtonPlus_perevent_ptmax2;
    fNoProtonMinus_ptmax2 = no_ProtonMinus_perevent_ptmax2;
    fNoProtonPlus_ptmax3 = no_ProtonPlus_perevent;
    fNoProtonMinus_ptmax3 = no_ProtonMinus_perevent;
    fTreeEvent->Fill();
    
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeEvent);
}    
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPbPbMC::GetEvent ()  //event cuts copied from my code written earlier 

{
 
  //Get Input Event
  fAODevent = dynamic_cast <AliAODEvent*>(InputEvent());
  if (!fAODevent)
    {
      AliWarning("ERROR: lAODevent not available from InputEvent(), trying with AODEvent() \n");

      fAODevent  = AODEvent();
      if(!fAODevent)
	{
	  AliWarning("ERROR: lAODevent not available from AODEvent() Aborting event!");
	  PostData(1, fOutputList);
	  PostData(2, fQAList);
	  PostData(3, fTreeEvent);
	  return kFALSE;
	}
    }

  fInputEvent = dynamic_cast <AliVEvent*>(InputEvent());
  if(!fInputEvent)
    {
      AliWarning("ERROR: fInputEvent (AliVEvent) not available \n");
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeEvent);
      return kFALSE;
    }
  
  hNumberOfEvents -> Fill(0.5);

    
  //Standard Event Cuts
  if (!fAODeventCuts.AcceptEvent(fInputEvent)) {
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeEvent);
    return kFALSE;
  }
  hNumberOfEvents -> Fill(1.5);
  
  
  //Reject Events with Incomplete DAQ
  if (fAODevent->IsIncompleteDAQ())
    {
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeEvent);
      return kFALSE;
    }
  hNumberOfEvents -> Fill(2.5);
        
  //V0 Timing Decision
  AliVVZERO *vzeroData = fInputEvent->GetVZEROData();
  if (!(vzeroData->GetV0ADecision()) || !(vzeroData->GetV0CDecision()))
    {
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeEvent);
      return kFALSE;
    }
  hNumberOfEvents -> Fill(3.5);
  
        
  //Pileup Rejection
  Int_t nClustersLayer0 = fInputEvent->GetNumberOfITSClusters(0);
  Int_t nClustersLayer1 = fInputEvent->GetNumberOfITSClusters(1);
  Int_t nTracklets      = fInputEvent->GetMultiplicity()->GetNumberOfTracklets();
  if ((nClustersLayer0 + nClustersLayer1) > 65.0 + (Double_t)nTracklets*4.0)
    {
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeEvent);
      return kFALSE;
    }
  hNumberOfEvents -> Fill(4.5);

    
  //Primary Vertex Tracks
  AliAODVertex *vertex_tracks = (AliAODVertex*) fAODevent->GetPrimaryVertexTracks();
  if (!vertex_tracks)
    {
      PostData(1, fOutputList);
      PostData(2, fQAList);
      PostData(3, fTreeEvent);
      return kFALSE;
    }
  hNumberOfEvents -> Fill(5.5);
  
        
    //Vertex Contributors Tracks
    if ( vertex_tracks->GetNContributors() < 1 )
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(6.5);

    
        
    //Primary Vertex SPD
    AliAODVertex *vertex_SPD = (AliAODVertex*) fAODevent->GetPrimaryVertexSPD();
    if (!vertex_SPD)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(7.5);
        
    //Vertex Contributors SPD
    if ( vertex_SPD->GetNContributors() < 1 )
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(8.5);

    
    //SPD Pile-up in Mult Bins
    if (fAODevent->IsPileupFromSPDInMultBins())
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(9.5);

    
    
    ////tigger/////////////                                                                                                                     
    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
    if ( !isSelected)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(10.5);
    
        
    //Cut on Z-Vertex Resolution
    if (TMath::Abs(vertex_SPD->GetZ() - vertex_tracks->GetZ()) > 0.3)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(11.5);

    //Primary Vertex Selection
    if ( vertex_tracks->GetZ() < -10.0 || vertex_tracks->GetZ() > +10.0)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(12.5);
               
    //Multiplicity
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fAODevent->FindListObject("MultSelection");
    if( !multiplicitySelection)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      } 
    hNumberOfEvents -> Fill(13.5);
                
    //Selection of Multiplicity Range
    Double_t mult_percentile = multiplicitySelection->GetMultiplicityPercentile("V0M");
    if (mult_percentile < 0.0 || mult_percentile > 100.0)
      {
	PostData(1, fOutputList);
	PostData(2, fQAList);
	PostData(3, fTreeEvent);
	return kFALSE;
      }
    hNumberOfEvents -> Fill(14.5);
    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    return kTRUE;
}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPbPbMC::PassedTrackQualityCuts (AliAODTrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);

    //select global tracks: ITS+TPC loose DCA
    // if(!track->TestFilterMask(BIT(4)))
    //   return passedTrkSelection;
    
    //Track Selection Cuts
    Int_t nTPCcluster = track->GetTPCNcls();  //TPC cluster cut
    if (nTPCcluster < 70)
      return passedTrkSelection;

    Double_t chi2TPCperClstr = track->GetTPCchi2perCluster();
    if(chi2TPCperClstr > 4)
      return passedTrkSelection;

    
    if(!(AliAODTrack::kTPCrefit)) //TPC refit
      return passedTrkSelection;

    if(track->GetKinkIndex(0) > 0) //No kink daughters
      return passedTrkSelection;

    if(!(AliAODTrack::kITSrefit)) //ITS refit
      return passedTrkSelection;
    
    
    //Kinematic cuts
    Double_t eta = track->Eta(); //Eta cut
    if(TMath::Abs(eta) > 0.8)
      return passedTrkSelection;

    Double_t pt = track->Pt(); //Pt cut
    if(pt < 0.15)
      return passedTrkSelection;

    
    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
 //_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPbPbMC::KaonSelector(AliVTrack *track)  {
 
  Double_t p[3];
  track->PxPyPz(p);

  Double_t Track_pt;
  Track_pt=TMath::Sqrt((p[0]*p[0])+(p[1]*p[1]));
  
  Double_t fTPCnSigmaPion = 0.0;
  Double_t fTPCnSigmaProton = 0.0;
  Double_t fTPCnSigmaKaon = 0.0;
  Double_t fTOFnSigmaPion = 0.0;
  Double_t fTOFnSigmaProton = 0.0;
  Double_t fTOFnSigmaKaon = 0.0;

  //TPC nsigma
  fTPCnSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  fTPCnSigmaProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  fTPCnSigmaKaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  //TOF nsigma
  fTOFnSigmaPion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
  fTOFnSigmaProton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
  fTOFnSigmaKaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);

  Double_t fTPCplusTOFnSigmaKaon = sqrt(TMath::Power(fTPCnSigmaKaon, 2.0) + TMath::Power(fTOFnSigmaKaon, 2.0));
  Double_t fTPCplusTOFnSigmaPion = sqrt(TMath::Power(fTPCnSigmaPion, 2.0) + TMath::Power(fTOFnSigmaPion, 2.0));
  Double_t fTPCplusTOFnSigmaProton = sqrt(TMath::Power(fTPCnSigmaProton, 2.0) + TMath::Power(fTOFnSigmaProton, 2.0));

 
  //Selection for pT < 0.5 : TPC only
  if( Track_pt < 0.5 )
    {
      /*
      //rejection
      Int_t flag = 0;
      if (TMath::Abs(fTPCnSigmaPion) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCnSigmaProton) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCnSigmaKaon) < 3.0)
	flag += 1;

      if (flag > 1) return kFALSE;

      if (fTPCnSigmaKaon > fTPCnSigmaProton) return kFALSE;
      if (fTPCnSigmaKaon > fTPCnSigmaPion) return kFALSE;
      */
      
      //acception
      if(TMath::Abs(fTPCnSigmaKaon) < 2.0)
	return kTRUE;
      else
	return kFALSE;
    }

  //Selection for pT > 0.5 : TOF + TPC
  if( Track_pt >= 0.5 )
    {
      //rejection
      Int_t flag = 0;
      if (TMath::Abs(fTPCplusTOFnSigmaPion) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCplusTOFnSigmaProton) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCplusTOFnSigmaKaon) < 3.0)
	flag += 1;

      if (flag > 1) return kFALSE;

      if (fTPCplusTOFnSigmaKaon > fTPCplusTOFnSigmaProton) return kFALSE;
      if (fTPCplusTOFnSigmaKaon > fTPCplusTOFnSigmaPion) return kFALSE;

      //acception
      if (fTPCplusTOFnSigmaKaon < 2.0)
	return kTRUE;
      else
	return kFALSE;
    }
  
  
  
}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPbPbMC::PionSelector(AliVTrack *track)  {
  
  Double_t p[3];
  track->PxPyPz(p);

  Double_t Track_pt;
  Track_pt=TMath::Sqrt((p[0]*p[0])+(p[1]*p[1]));
  
  Double_t fTPCnSigmaPion = 0.0;
  Double_t fTPCnSigmaProton = 0.0;
  Double_t fTPCnSigmaKaon = 0.0;
  Double_t fTOFnSigmaPion = 0.0;
  Double_t fTOFnSigmaProton = 0.0;
  Double_t fTOFnSigmaKaon = 0.0;

  //TPC nsigma
  fTPCnSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  fTPCnSigmaProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  fTPCnSigmaKaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  //TOF nsigma
  fTOFnSigmaPion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
  fTOFnSigmaProton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
  fTOFnSigmaKaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);

  Double_t fTPCplusTOFnSigmaPion = sqrt(TMath::Power(fTPCnSigmaPion, 2.0) + TMath::Power(fTOFnSigmaPion, 2.0));
  Double_t fTPCplusTOFnSigmaKaon = sqrt(TMath::Power(fTPCnSigmaKaon, 2.0) + TMath::Power(fTOFnSigmaKaon, 2.0));
  Double_t fTPCplusTOFnSigmaProton = sqrt(TMath::Power(fTPCnSigmaProton, 2.0) + TMath::Power(fTOFnSigmaProton, 2.0));

  


  //Selection for pT < 0.5 : TPC only
  if( Track_pt < 0.5 )
    {
      /*
      //rejection
      
      Int_t flag = 0;
      if (TMath::Abs(fTPCnSigmaPion) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCnSigmaProton) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCnSigmaKaon) < 3.0)
	flag += 1;

      if (flag > 1) return kFALSE;

      if (fTPCnSigmaPion > fTPCnSigmaProton) return kFALSE;
      if (fTPCnSigmaPion > fTPCnSigmaKaon) return kFALSE;
      */
      
      //acception
      
      if(TMath::Abs(fTPCnSigmaPion) < 2.0)
	return kTRUE;
      else
	return kFALSE;
    }

  //Selection for pT > 0.5 : TOF + TPC
  if( Track_pt >= 0.5 )
    {
      //rejection
      
      Int_t flag = 0;
      if (TMath::Abs(fTPCplusTOFnSigmaPion) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCplusTOFnSigmaProton) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCplusTOFnSigmaKaon) < 3.0)
	flag += 1;

      if (flag > 1) return kFALSE;

      if (fTPCplusTOFnSigmaPion > fTPCplusTOFnSigmaProton) return kFALSE;
      if (fTPCplusTOFnSigmaPion > fTPCplusTOFnSigmaKaon) return kFALSE;

      //acception
      
      if (fTPCplusTOFnSigmaPion < 2.0)
	return kTRUE;
      else
	return kFALSE;
    }
  
}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPbPbMC::ProtonSelector(AliVTrack *track)  {
  
  Double_t p[3];
  track->PxPyPz(p);

  Double_t Track_pt;
  Track_pt=TMath::Sqrt((p[0]*p[0])+(p[1]*p[1]));
  
  Double_t fTPCnSigmaPion = 0.0;
  Double_t fTPCnSigmaProton = 0.0;
  Double_t fTPCnSigmaKaon = 0.0;
  Double_t fTOFnSigmaPion = 0.0;
  Double_t fTOFnSigmaProton = 0.0;
  Double_t fTOFnSigmaKaon = 0.0;

  //TPC nsigma
  fTPCnSigmaPion = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
  fTPCnSigmaProton = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  fTPCnSigmaKaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon);
  //TOF nsigma
  fTOFnSigmaPion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
  fTOFnSigmaProton = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
  fTOFnSigmaKaon = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);

  Double_t fTPCplusTOFnSigmaPion = sqrt(TMath::Power(fTPCnSigmaPion, 2.0) + TMath::Power(fTOFnSigmaPion, 2.0));
  Double_t fTPCplusTOFnSigmaKaon = sqrt(TMath::Power(fTPCnSigmaKaon, 2.0) + TMath::Power(fTOFnSigmaKaon, 2.0));
  Double_t fTPCplusTOFnSigmaProton = sqrt(TMath::Power(fTPCnSigmaProton, 2.0) + TMath::Power(fTOFnSigmaProton, 2.0));

  
  //Selection for pT < 0.5 : TPC only
  if( Track_pt < 0.6 )
    {
      //rejection

      /*
      Int_t flag = 0;
      if (TMath::Abs(fTPCnSigmaPion) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCnSigmaProton) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCnSigmaKaon) < 3.0)
	flag += 1;

      if (flag > 1) return kFALSE;

      if (fTPCnSigmaProton > fTPCnSigmaPion) return kFALSE;
      if (fTPCnSigmaProton > fTPCnSigmaKaon) return kFALSE;
      */
      
      //acception

      if(TMath::Abs(fTPCnSigmaProton) < 2.0)
	return kTRUE;
      else
	return kFALSE;
    }

  //Selection for pT > 0.5 : TOF + TPC
  if( Track_pt >= 0.6 )
    {
      //rejection

      
      Int_t flag = 0;
      if (TMath::Abs(fTPCplusTOFnSigmaPion) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCplusTOFnSigmaProton) < 3.0)
	flag += 1;
      if (TMath::Abs(fTPCplusTOFnSigmaKaon) < 3.0)
	flag += 1;

      if (flag > 1) return kFALSE;

      if (fTPCplusTOFnSigmaProton > fTPCplusTOFnSigmaPion) return kFALSE;
      if (fTPCplusTOFnSigmaProton > fTPCplusTOFnSigmaKaon) return kFALSE;
      

      //acception
      
      if (fTPCplusTOFnSigmaProton < 2.0)
	return kTRUE;
      else
	return kFALSE;
    }
  
}


//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPbPbMC::PassedPIDSelection (AliAODTrack *track, AliPID::EParticleType type)  {
    
    //Initialization
    Bool_t passedPIDSelection=(kFALSE);
    
    //TPC Particle Identification
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,type);
    if (nsigmaTPC < -4.0) return passedPIDSelection;
    if (nsigmaTPC > +4.0) return passedPIDSelection;

    passedPIDSelection = kTRUE;
    return passedPIDSelection;
}
 //_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPbPbMC::PassedSingleParticlePileUpCuts(AliAODTrack *track)
{
  Bool_t passedTrackPileupCut = (kTRUE);
  if (!(track->HasPointOnITSLayer(1)) && !(track->HasPointOnITSLayer(4)) && !(track->HasPointOnITSLayer(5)) && !(track->GetTOFBunchCrossing() == 0))
    passedTrackPileupCut = (kFALSE);
  return passedTrackPileupCut;
}
 //_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskCorrPbPbMC::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
