#include <Riostream.h>
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenHepMCEventHeader.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"

#include "AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"



using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen)

AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen::AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen(): AliAnalysisTaskSE(), fTreeEvent(0), fListHist(0), fPIDResponse(0), fESDtrackCuts(0), fEventCuts(0), fTreeVariableCentrality(0), fNoGenKaonPlus_ptmax2(0), fNoGenKaonMinus_ptmax2(0), fNoGenKaonPlus_ptmax3(0), fNoGenKaonMinus_ptmax3(0), fNoGenPionPlus_ptmax2(0), fNoGenPionMinus_ptmax2(0), fNoGenPionPlus_ptmax3(0), fNoGenPionMinus_ptmax3(0), fNoGenProtonPlus_ptmax2(0), fNoGenProtonMinus_ptmax2(0), fNoGenProtonPlus_ptmax3(0), fNoGenProtonMinus_ptmax3(0), fvertex(0), hist2D_pt_gen_centrality(0), hist2D_pt_rec_centrality(0), hist_centrality_beforecut(0), fNch_eta0pt5(0),fMCchoice(0), fNoResoChoice(0), fNpart_1(0),fNpart_2(0)
{
  
}

AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen::AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen(const char *name): AliAnalysisTaskSE(name), fTreeEvent(0), fListHist(0), fPIDResponse(0), fESDtrackCuts(0), fEventCuts(0), fTreeVariableCentrality(0), fNoGenKaonPlus_ptmax2(0), fNoGenKaonMinus_ptmax2(0), fNoGenKaonPlus_ptmax3(0), fNoGenKaonMinus_ptmax3(0), fNoGenPionPlus_ptmax2(0), fNoGenPionMinus_ptmax2(0), fNoGenPionPlus_ptmax3(0), fNoGenPionMinus_ptmax3(0), fNoGenProtonPlus_ptmax2(0), fNoGenProtonMinus_ptmax2(0), fNoGenProtonPlus_ptmax3(0), fNoGenProtonMinus_ptmax3(0), fvertex(0), hist2D_pt_gen_centrality(0), hist2D_pt_rec_centrality(0), hist_centrality_beforecut(0), fNch_eta0pt5(0), fNoResoChoice(0), fMCchoice(0), fNpart_1(0), fNpart_2(0)
{

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
    
}

AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen::~AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen()
{
  //------------------------------------------------
  // DESTRUCTOR
  //------------------------------------------------
    
  if (fTreeEvent){
    delete fTreeEvent;
    fTreeEvent = 0x0;
  }
  if (fListHist){
    delete fListHist;
    fListHist = 0x0;
  }
  if(fESDtrackCuts) {
    delete fESDtrackCuts;
    fESDtrackCuts=0x0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen::UserCreateOutputObjects()
{
  
  OpenFile(1);
  fTreeEvent = new TTree("fTreeEvent","Event");
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

  PostData(1,fTreeEvent);
   
  
  OpenFile(2);
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);

  hist2D_pt_gen_centrality=new TH2F("hist2D_pt_gen_centrality","hist2D_pt_gen_centrality",800,0,8,100,0,100);
  hist2D_pt_rec_centrality=new TH2F("hist2D_pt_rec_centrality","hist2D_pt_rec_centrality",800,0,8,100,0,100);
  hist_centrality_beforecut=new TH1D("hist_centrality_beforecut","hist_centrality_beforecut",100,0,100);

  fListHist->Add(hist2D_pt_gen_centrality);
  fListHist->Add(hist2D_pt_rec_centrality);
  fListHist->Add(hist_centrality_beforecut);
  
  PostData(2, fListHist); 
  
}


//________________________________________________________________________
void AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen::UserExec(Option_t *)
{



  //cout<<"************************************* Event *************************************************"<<endl;
  // Main loop                                                                                                                              
  // Called for each event                                                                                                                  
  AliESDEvent *lESDevent = 0x0;
  AliMCEvent  *lMCevent  = 0x0;
  AliStack    *lMCstack  = 0x0;

  
  lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!lESDevent) {
    AliWarning("ERROR: lESDevent not available \n");
    return;
  }
  
  lMCevent = dynamic_cast<AliMCEvent *>(MCEvent());
  if (!lMCevent) {
    Printf("ERROR: Could not retrieve MC event \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }

  lMCstack = lMCevent->Stack();
  if (!lMCstack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }


  //////centrality selection/////////                                                                                                         

   
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

  if (eventHandler) lMCevent = eventHandler->MCEvent();
  AliGenEventHeader* genHeader = 0x0;
  TString genheaderName;
  if (lMCevent){
    genHeader = lMCevent->GenEventHeader();
    genheaderName = genHeader->GetName();
    if(!genHeader){ Printf(" Error::marsland: Event generator header not available!!!\n"); return; }
  }


  Double_t impParArr[10] = {0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5};
  AliGenHijingEventHeader *lHIJINGHeader = 0x0;  // event header for HIJING
  AliGenHepMCEventHeader *lHepMCHeader = 0x0;    // event header for EPOS
   
  Float_t fMCImpactParameter=0.0,fCentImpBin=-10.0;
  Float_t fNHardScatters=0.0;
  Float_t fNProjectileParticipants=0.0;
  Float_t fNTargetParticipants=0.0;
  Float_t fNNColl=0.0;
  Float_t fNNwColl=0.0;
  Float_t fNwNColl=0.0;
  Float_t fNwNwColl=0.0;
  Float_t fNpart=0.0;

    
  //For EPOS
    
  if(fMCchoice==1)
    {
      lHepMCHeader = (AliGenHepMCEventHeader*)genHeader;
      fNHardScatters = lHepMCHeader->Ncoll_hard(); // Number of hard scatterings
      fNProjectileParticipants = lHepMCHeader->Npart_proj(); // Number of projectile participants
      fNTargetParticipants     = lHepMCHeader->Npart_targ(); // Number of target participants
      fNNColl   = lHepMCHeader->Ncoll(); // Number of NN (nucleon-nucleon) collisions
      fNNwColl  = lHepMCHeader->N_Nwounded_collisions(); // Number of N-Nwounded collisions
      fNwNColl  = lHepMCHeader->Nwounded_N_collisions(); // Number of Nwounded-N collisons
      fNwNwColl = lHepMCHeader->Nwounded_Nwounded_collisions();// Number of Nwounded-Nwounded collisions
        
      fMCImpactParameter = lHepMCHeader->impact_parameter(); //impact parameter
      if (fMCImpactParameter>=impParArr[0] && fMCImpactParameter<impParArr[1]) fCentImpBin=2.5;
      if (fMCImpactParameter>=impParArr[1] && fMCImpactParameter<impParArr[2]) fCentImpBin=7.5;
      if (fMCImpactParameter>=impParArr[2] && fMCImpactParameter<impParArr[3]) fCentImpBin=15.;
      if (fMCImpactParameter>=impParArr[3] && fMCImpactParameter<impParArr[4]) fCentImpBin=25.;
      if (fMCImpactParameter>=impParArr[4] && fMCImpactParameter<impParArr[5]) fCentImpBin=35.;
      if (fMCImpactParameter>=impParArr[5] && fMCImpactParameter<impParArr[6]) fCentImpBin=45.;
      if (fMCImpactParameter>=impParArr[6] && fMCImpactParameter<impParArr[7]) fCentImpBin=55.;
      if (fMCImpactParameter>=impParArr[7] && fMCImpactParameter<impParArr[8]) fCentImpBin=65.;
      if (fMCImpactParameter>=impParArr[8] && fMCImpactParameter<impParArr[9]) fCentImpBin=75.;
      if (fMCImpactParameter<impParArr[0]  || fMCImpactParameter>impParArr[9]) fCentImpBin=-10.;
    }
    
    
  //For HIJING

  if(fMCchoice==2)
    {
      lHIJINGHeader = (AliGenHijingEventHeader*) genHeader;
      fNHardScatters = lHIJINGHeader->HardScatters();
      fNProjectileParticipants = lHIJINGHeader->ProjectileParticipants();
      fNTargetParticipants     = lHIJINGHeader->TargetParticipants();
      fNpart = lHIJINGHeader->GetTrueNPart();
      fNNColl   = lHIJINGHeader->NN();
      fNNwColl  = lHIJINGHeader->NNw();
      fNwNColl  = lHIJINGHeader->NwN();
      fNwNwColl = lHIJINGHeader->NwNw();
      fMCImpactParameter = lHIJINGHeader->ImpactParameter();
      if (fMCImpactParameter>=impParArr[0] && fMCImpactParameter<impParArr[1]) fCentImpBin=2.5;
      if (fMCImpactParameter>=impParArr[1] && fMCImpactParameter<impParArr[2]) fCentImpBin=7.5;
      if (fMCImpactParameter>=impParArr[2] && fMCImpactParameter<impParArr[3]) fCentImpBin=15.;
      if (fMCImpactParameter>=impParArr[3] && fMCImpactParameter<impParArr[4]) fCentImpBin=25.;
      if (fMCImpactParameter>=impParArr[4] && fMCImpactParameter<impParArr[5]) fCentImpBin=35.;
      if (fMCImpactParameter>=impParArr[5] && fMCImpactParameter<impParArr[6]) fCentImpBin=45.;
      if (fMCImpactParameter>=impParArr[6] && fMCImpactParameter<impParArr[7]) fCentImpBin=55.;
      if (fMCImpactParameter>=impParArr[7] && fMCImpactParameter<impParArr[8]) fCentImpBin=65.;
      if (fMCImpactParameter>=impParArr[8] && fMCImpactParameter<impParArr[9]) fCentImpBin=75.;
      if (fMCImpactParameter<impParArr[0]  || fMCImpactParameter>impParArr[9]) fCentImpBin=-10.;
    }

  if(fCentImpBin < 0.0)
    {
      PostData(1,fTreeEvent);
      return;
    }
     
    
  AliMCParticle *trackMCgen;

  Int_t trkCharge=0;
  float trkPt=0.0;
  float_t trkEta=0.0;
  Int_t trkPdgcode=0;
  float Q1[2]={0.0,0.0};
  float Q2[2]={0.0,0.0};
  float Q3[2]={0.0,0.0};
  float Q4[2]={0.0,0.0};
  float Nch[2]={0.0,0.0};
  float Nch_eta0pt5=0.0;
      
 

  Float_t no_GenKaonPlus_perevent = 0;
  Float_t no_GenKaonMinus_perevent = 0;
  Float_t no_GenProtonPlus_perevent = 0;
  Float_t no_GenProtonMinus_perevent = 0;
  Float_t no_GenPionPlus_perevent = 0;
  Float_t no_GenPionMinus_perevent = 0;

  Float_t no_GenKaonPlus_perevent_ptmax2 = 0;
  Float_t no_GenKaonMinus_perevent_ptmax2 = 0;
  Float_t no_GenProtonPlus_perevent_ptmax2 = 0;
  Float_t no_GenProtonMinus_perevent_ptmax2 = 0;
  Float_t no_GenPionPlus_perevent_ptmax2 = 0;
  Float_t no_GenPionMinus_perevent_ptmax2 = 0;

  for (Int_t iTrack = 0; iTrack < lMCevent->GetNumberOfTracks(); iTrack++)
    {
      trackMCgen = (AliMCParticle *)lMCevent->GetTrack(iTrack);
      if (!trackMCgen) continue;
      // 
      // Accepty only primary particles
      if(fNoResoChoice==1)
	{
	  if (!lMCstack->IsPhysicalPrimary(iTrack)) continue;
	}
      //
      //cout<<iTrack<<"\t"<<"Pt: "<<trackMCgen->Pt()<<"\t"<<"Eta: "<<trackMCgen->Eta()<<"\t"<<"Charge: "<<trackMCgen->Charge()<<endl;

      trkCharge=1.0/3.0*trackMCgen->Charge();
      trkPt=trackMCgen->Pt();
      trkEta=trackMCgen->Eta();
      trkPdgcode=trackMCgen->PdgCode();

      if(TMath::Abs(trkEta)>0.8)continue;
      if(TMath::Abs(trkCharge)!=1)continue;


      //Kaon
      if (TMath::Abs(trkPdgcode) == 321)
	{
	  if(trkCharge > 0)  //K+
	    {
	      if(trkPt > 0.4 && trkPt < 1.6)
		no_GenKaonPlus_perevent += 1;
	      if(trkPt < 2.0)
		no_GenKaonPlus_perevent_ptmax2 += 1;
	    }
	  if(trkCharge < 0)  //K-
	    {
	      if(trkPt > 0.4 && trkPt < 1.6)
		no_GenKaonMinus_perevent += 1;
	      if(trkPt < 2.0)
		no_GenKaonMinus_perevent_ptmax2 += 1;
	    }
	}

	//Proton
	if (TMath::Abs(trkPdgcode) == 2212)
	  {
	    if(trkCharge > 0)  //p
	      {
		if(trkPt > 0.4)
		  {
		    if(trkPt > 0.4 && trkPt < 1.6)
		      no_GenProtonPlus_perevent += 1;
		    if(trkPt < 2.0)
		      no_GenProtonPlus_perevent_ptmax2 += 1;
		  }
	      }
	    if(trkCharge < 0)  //pbar
	      {
		if(trkPt > 0.4)
		  {
		    if(trkPt > 0.4 && trkPt < 1.6)
		      no_GenProtonMinus_perevent += 1;
		    if(trkPt < 2.0)
		      no_GenProtonMinus_perevent_ptmax2 += 1;
		  }
	      }
	  }

	//Pion
	if (TMath::Abs(trkPdgcode) == 211)
	  {
	    if(trkCharge > 0)  //pi+
	      {
		if(trkPt > 0.4 && trkPt < 1.6)
		  no_GenPionPlus_perevent += 1;
		if(trkPt < 2.0)
		  no_GenPionPlus_perevent_ptmax2 += 1;
	      }
	    if(trkCharge < 0)  //pi-
	      {
		if(trkPt > 0.4 && trkPt < 1.6)
		  no_GenPionMinus_perevent += 1;
		if(trkPt < 2.0)
		  no_GenPionMinus_perevent_ptmax2 += 1;
	      }
	  }
	  
	  
	  
    }//end track loop

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  // cout<<"Value of centrality:"<<fCentImpBin<<endl;
  // cout<<"No of particaipating nucleons: "<<fNpart<<endl;
  // cout<<"No of projectile + target nucleons: "<<fNProjectileParticipants<<" + "<<fNTargetParticipants<<endl;
  // cout<<"No of collisions: "<<fNNColl<<endl;

  if(fMCchoice==2)
    {
      fNpart_1=fNpart;
    }
  else
    fNpart_1=0;
  //fNpart_1=fNpart;
  fNpart_2=fNProjectileParticipants+fNTargetParticipants;

  //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  //Tree variables
      
  fTreeVariableCentrality=fCentImpBin;

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
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

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  
  fTreeEvent->Fill();
  PostData(1,fTreeEvent);
      

}
//________________________________________________________________________
void AliAnalysisTaskNetPiKProtFluc_MCnoESD_gen::Terminate(Option_t *)
{
  
}
//----------------------------------------------------------------------------
