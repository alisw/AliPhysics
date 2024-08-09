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

#include "AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"



using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2)

AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2::AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2(): AliAnalysisTaskSE(), fTreeEvent(0), fListHist(0), fPIDResponse(0), fESDtrackCuts(0), fEventCuts(0), fTreeVariableCentrality(0), fPtsum_hadrons_less0(0), fPtsum_hadrons_greaterEtaMin(0), fNsum_hadrons_less0(0), fNsum_hadrons_greaterEtaMin(0), fNsum_pions_less0(0), fNsum_kaons_less0(0), fNsum_protons_less0(0), hist2D_pt_gen_centrality(0), hist2D_pt_rec_centrality(0), hist_centrality_beforecut(0), fNch_eta0pt5(0), fMCchoice(0), fNpart_1(0), fNpart_2(0), fEtaLeftCut(0), fEtaRightCut(0)
{
  for(int i=0; i<20; i++)
    {
      fPt_no_hadron[i] = 0;
      fPt_no_pion[i] = 0;
      fPt_no_kaon[i] = 0;
      fPt_no_proton[i] = 0;
    }
}

AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2::AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2(const char *name): AliAnalysisTaskSE(name), fTreeEvent(0), fListHist(0), fPIDResponse(0), fESDtrackCuts(0), fEventCuts(0), fTreeVariableCentrality(0), fPtsum_hadrons_less0(0), fPtsum_hadrons_greaterEtaMin(0), fNsum_hadrons_less0(0), fNsum_hadrons_greaterEtaMin(0), fNsum_pions_less0(0), fNsum_kaons_less0(0), fNsum_protons_less0(0), hist2D_pt_gen_centrality(0), hist2D_pt_rec_centrality(0), hist_centrality_beforecut(0), fNch_eta0pt5(0), fMCchoice(0), fNpart_1(0), fNpart_2(0), fEtaLeftCut(0), fEtaRightCut(0)
{
  for(int i=0; i<20; i++)
    {
      fPt_no_hadron[i] = 0;
      fPt_no_pion[i] = 0;
      fPt_no_kaon[i] = 0;
      fPt_no_proton[i] = 0;
    }
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
    
}

AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2::~AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2()
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
void AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2::UserCreateOutputObjects()
{
  
  OpenFile(1);
  fTreeEvent = new TTree("fTreeEvent","Event");
  fTreeEvent->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
  fTreeEvent->Branch("fPtsum_hadrons_less0",&fPtsum_hadrons_less0,"fPtsum_hadrons_less0/F");
  fTreeEvent->Branch("fPtsum_hadrons_greaterEtaMin",&fPtsum_hadrons_greaterEtaMin,"fPtsum_hadrons_greaterEtaMin/F");
  fTreeEvent->Branch("fNsum_hadrons_less0",&fNsum_hadrons_less0,"fNsum_hadrons_less0/F");
  fTreeEvent->Branch("fNsum_hadrons_greaterEtaMin",&fNsum_hadrons_greaterEtaMin,"fNsum_hadrons_greaterEtaMin/F");
  fTreeEvent->Branch("fNsum_pions_less0",&fNsum_pions_less0,"fNsum_pions_less0/F");
  fTreeEvent->Branch("fNsum_kaons_less0",&fNsum_kaons_less0,"fNsum_kaons_less0/F");
  fTreeEvent->Branch("fNsum_protons_less0",&fNsum_protons_less0,"fNsum_protons_less0/F");
  fTreeEvent->Branch("fPt_no_hadron",&fPt_no_hadron,"fPt_no_hadron[20]/F");
  fTreeEvent->Branch("fPt_no_pion",&fPt_no_pion,"fPt_no_pion[20]/F");
  fTreeEvent->Branch("fPt_no_kaon",&fPt_no_kaon,"fPt_no_kaon[20]/F");
  fTreeEvent->Branch("fPt_no_proton",&fPt_no_proton,"fPt_no_proton[20]/F");


  PostData(1, fTreeEvent);
  
  
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
void AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2::UserExec(Option_t *)
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
      
  double binsarray[21]={0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.5,4.0,5.0,6.0,8.0,10.0};
  
  TH1D *fPt_profile = new TH1D("fPt_profile","fPt_profile", 20, binsarray);
  Double_t pT_sum_etaLess0 = 0.0;
  Double_t N_sum_etaLess0 = 0.0;
  Double_t pT_sum_etaGreaterEtamin = 0.0;
  Double_t N_sum_etaGreaterEtamin = 0.0;


  TH1D *fPt_profile_pion = new TH1D("fPt_profile_pion","fPt_profile_pion", 20, binsarray);
  TH1D *fPt_profile_kaon = new TH1D("fPt_profile_kaon","fPt_profile_kaon", 20, binsarray);
  TH1D *fPt_profile_proton = new TH1D("fPt_profile_proton","fPt_profile_proton", 20, binsarray);
  Double_t N_sumPion_etaLess0 = 0.0;
  Double_t N_sumKaon_etaLess0 = 0.0;
  Double_t N_sumProton_etaLess0 = 0.0; 

  for (Int_t iTrack = 0; iTrack < lMCevent->GetNumberOfTracks(); iTrack++)
    {
      trackMCgen = (AliMCParticle *)lMCevent->GetTrack(iTrack);
      if (!trackMCgen) continue;
      // 
      // Accepty only primary particles
      if (!lMCstack->IsPhysicalPrimary(iTrack)) continue;
      //
      //cout<<iTrack<<"\t"<<"Pt: "<<trackMCgen->Pt()<<"\t"<<"Eta: "<<trackMCgen->Eta()<<"\t"<<"Charge: "<<trackMCgen->Charge()<<endl;

      trkCharge=1.0/3.0*trackMCgen->Charge();
      trkPt=trackMCgen->Pt();
      trkEta=trackMCgen->Eta();
      trkPdgcode=trackMCgen->PdgCode();

      if(TMath::Abs(trkEta)>0.8)continue;
      if(TMath::Abs(trkCharge)!=1)continue;
      if(trkPt<0.2)continue;
      if(trkPt>10.0)continue;

      hist2D_pt_gen_centrality->Fill(trkPt,fCentImpBin);

      if(TMath::Abs(trkCharge) > 0)
	{
	  if(trkEta < fEtaLeftCut)
	    {
	      fPt_profile->Fill(trkPt);
	      pT_sum_etaLess0 += trkPt;
	      N_sum_etaLess0 += 1.0;

	      if(TMath::Abs(trkPdgcode)==211)
		{
		  fPt_profile_pion->Fill(trkPt);
		  N_sumPion_etaLess0 += 1.0;
		}
	      if(TMath::Abs(trkPdgcode)==321)
		{
		  fPt_profile_kaon->Fill(trkPt);
		  N_sumKaon_etaLess0 += 1.0;
		}
	      if(TMath::Abs(trkPdgcode)==2212)
		{
		  fPt_profile_proton->Fill(trkPt);
		  N_sumProton_etaLess0 += 1.0;
		}

	    }
	    
	  if(trkEta > fEtaRightCut)
	    {
	      pT_sum_etaGreaterEtamin += trkPt;
	      N_sum_etaGreaterEtamin += 1.0;
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
  hist_centrality_beforecut->Fill(fCentImpBin);   
 
  fTreeVariableCentrality=fCentImpBin;
  fPtsum_hadrons_less0=pT_sum_etaLess0;
  fPtsum_hadrons_greaterEtaMin=pT_sum_etaGreaterEtamin;
  fNsum_hadrons_less0=N_sum_etaLess0;
  fNsum_hadrons_greaterEtaMin=N_sum_etaGreaterEtamin;
  fNsum_pions_less0=N_sumPion_etaLess0;
  fNsum_kaons_less0=N_sumKaon_etaLess0;
  fNsum_protons_less0=N_sumProton_etaLess0;
  
  for(int i=0; i<20; i++)
    {
      fPt_no_hadron[i]=fPt_profile->GetBinContent(i+1);
      fPt_no_pion[i]=fPt_profile_pion->GetBinContent(i+1);
      fPt_no_kaon[i]=fPt_profile_kaon->GetBinContent(i+1);
      fPt_no_proton[i]=fPt_profile_proton->GetBinContent(i+1);
    }
  fTreeEvent->Fill();
    

  fPt_profile->Delete();
  fPt_profile_pion->Delete();
  fPt_profile_kaon->Delete();
  fPt_profile_proton->Delete();
  
  PostData(1,fTreeEvent);

}
//________________________________________________________________________
void AliAnalysisTaskDiffPtFluc_MCnoESD_gen_v2::Terminate(Option_t *)
{
  
}
//----------------------------------------------------------------------------
