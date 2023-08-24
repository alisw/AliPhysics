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

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"

#include "AliAnalysisTaskResonanceVsMultiplicityMC.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"



using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskResonanceVsMultiplicityMC)

AliAnalysisTaskResonanceVsMultiplicityMC::AliAnalysisTaskResonanceVsMultiplicityMC(): AliAnalysisTaskSE(), fTreeEvent(0), fListHist(0), fPIDResponse(0), fESDtrackCuts(0), fEventCuts(0), fTreeVariableCentrality(0), fvertex(0), Profile_mean_term1(0),Profile_var_term1(0),Profile_var_term2(0),Profile_skewness_term1(0),Profile_skewness_term2(0),Profile_skewness_term3(0),Profile_kurtosis_term1(0),Profile_kurtosis_term2(0),Profile_kurtosis_term3(0),Profile_kurtosis_term4(0), hist2D_pt_gen_centrality(0), hist2D_pt_rec_centrality(0), hist_centrality_beforecut(0), fEta_max(0), fDCAxy_max(0), fDCAz_max(0), fchi2TPC_max(0), fchi2ITS_max(0), fNCR_min(0), fVertexZ_max(0), fTreeName(0)
{
  for(int i=0;i<4;i++)
    {
      fQ1_gen[i]=-999;
      fQ2_gen[i]=-999;
      fQ3_gen[i]=-999;
      fQ4_gen[i]=-999;
      fNch_gen[i]=-999;
      fQ1_rec[i]=-999;
      fQ2_rec[i]=-999;
      fQ3_rec[i]=-999;
      fQ4_rec[i]=-999;
      fNch_rec[i]=-999;
    }
}

AliAnalysisTaskResonanceVsMultiplicityMC::AliAnalysisTaskResonanceVsMultiplicityMC(const char *name): AliAnalysisTaskSE(name), fTreeEvent(0), fListHist(0), fPIDResponse(0), fESDtrackCuts(0), fEventCuts(0), fTreeVariableCentrality(0), fvertex(0),Profile_mean_term1(0),Profile_var_term1(0),Profile_var_term2(0),Profile_skewness_term1(0),Profile_skewness_term2(0),Profile_skewness_term3(0),Profile_kurtosis_term1(0),Profile_kurtosis_term2(0),Profile_kurtosis_term3(0),Profile_kurtosis_term4(0), hist2D_pt_gen_centrality(0), hist2D_pt_rec_centrality(0), hist_centrality_beforecut(0), fEta_max(0), fDCAxy_max(0), fDCAz_max(0), fchi2TPC_max(0), fchi2ITS_max(0), fNCR_min(0), fVertexZ_max(0), fTreeName(0)
{
  for(int i=0;i<4;i++)
    {
      fQ1_gen[i]=-999;
      fQ2_gen[i]=-999;
      fQ3_gen[i]=-999;
      fQ4_gen[i]=-999;
      fNch_gen[i]=-999;
      fQ1_rec[i]=-999;
      fQ2_rec[i]=-999;
      fQ3_rec[i]=-999;
      fQ4_rec[i]=-999;
      fNch_rec[i]=-999;
    }

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
    
}

AliAnalysisTaskResonanceVsMultiplicityMC::~AliAnalysisTaskResonanceVsMultiplicityMC()
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
void AliAnalysisTaskResonanceVsMultiplicityMC::UserCreateOutputObjects()
{

  //------------------------------------------------
  // Particle Identification Setup
  //------------------------------------------------
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  
   //------------------------------------------------
  // track cut
  //------------------------------------------------
  if(!fESDtrackCuts )
    {
      fESDtrackCuts = new AliESDtrackCuts();
      fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
      fESDtrackCuts->SetPtRange(0.15,8.0);
      fESDtrackCuts->SetEtaRange(-1.0*fEta_max,fEta_max);

      //++++++++++++++Track parameters to be varied for systematics++++++++++++
      fESDtrackCuts->SetMaxDCAToVertexXY(fDCAxy_max);
      fESDtrackCuts->SetMaxDCAToVertexZ(fDCAz_max);
      fESDtrackCuts->SetMaxChi2PerClusterTPC(fchi2TPC_max);
      fESDtrackCuts->SetMaxChi2PerClusterITS(fchi2ITS_max);
      fESDtrackCuts->SetMinNCrossedRowsTPC(fNCR_min);
      //fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    }

  
  OpenFile(1);
  fTreeEvent = new TTree(fTreeName,"Event");
  fTreeEvent->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
  fTreeEvent->Branch("fvertex",&fvertex,"fvertex/F");
  fTreeEvent->Branch("fNch_gen", &fNch_gen, "fNch_gen[4]/F");
  fTreeEvent->Branch("fQ1_gen", &fQ1_gen, "fQ1_gen[4]/F");
  fTreeEvent->Branch("fQ2_gen", &fQ2_gen, "fQ2_gen[4]/F");
  fTreeEvent->Branch("fQ3_gen", &fQ3_gen, "fQ3_gen[4]/F");
  fTreeEvent->Branch("fQ4_gen", &fQ4_gen, "fQ4_gen[4]/F");
  fTreeEvent->Branch("fNch_rec", &fNch_rec, "fNch_rec[4]/F");
  fTreeEvent->Branch("fQ1_rec", &fQ1_rec, "fQ1_rec[4]/F");
  fTreeEvent->Branch("fQ2_rec", &fQ2_rec, "fQ2_rec[4]/F");
  fTreeEvent->Branch("fQ3_rec", &fQ3_rec, "fQ3_rec[4]/F");
  fTreeEvent->Branch("fQ4_rec", &fQ4_rec, "fQ4_rec[4]/F");
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
void AliAnalysisTaskResonanceVsMultiplicityMC::UserExec(Option_t *)
{

  cout<<"************************************* Event *************************************************"<<endl;
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
  

    ////tigger/////////////
  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isSelected = 0;
  isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
  
  if ( ! isSelected)
    {
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      return;
    }
  
  // primary vertex
  //
  const AliESDVertex *vertex = lESDevent->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1)
    {
      // SPD vertex
      vertex = lESDevent->GetPrimaryVertexSPD();
      if (vertex->GetNContributors() < 1)
	vertex = 0x0;
    }
  if (!vertex)
    {
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      return;
    }  
  if (TMath::Abs(vertex->GetZ()) > fVertexZ_max)
    {
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      return;
    }


  //////centrality selection/////////
  Float_t lV0M;
  Int_t lEvSelCode = 100;
  AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
  if( !MultSelection)
    {
      AliWarning("AliMultSelection object not found!");
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      return;
    }
  else
    {
      lV0M = MultSelection->GetMultiplicityPercentile("V0M"); 
    }

  hist_centrality_beforecut->Fill(lV0M);
  
  Bool_t EventAccepted;
  EventAccepted = fEventCuts.AcceptEvent(lESDevent);
  if (!EventAccepted)
    {
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      return;
    }
  
  //Filling tree variables
  fTreeVariableCentrality=lV0M;
  fvertex=vertex->GetZ();

  

  //////////loop for generated////////////////
   cout<<"*********************************** Generated ********************************************"<<endl;

   //Same bunch pileup event in MC 
   //if(AliAnalysisUtils::IsSameBunchPileupInGeneratedEvent(lMCevent))
   //  continue;


   float Q1_gen[4]={0.0,0.0,0.0,0.0};
   float Q2_gen[4]={0.0,0.0,0.0,0.0};
   float Q3_gen[4]={0.0,0.0,0.0,0.0};
   float Q4_gen[4]={0.0,0.0,0.0,0.0};
   float Nch_gen[4]={0.0,0.0,0.0,0.0};

   cout<<"Generated No of tracks:*****************"<<lMCstack->GetNtrack()<<endl;
   for (Int_t ipart = 0; ipart < lMCstack->GetNtrack(); ipart++)
     {
       // This is the begining of the loop on tracks
       TParticle* particle = lMCstack->Particle(ipart);
       if(!particle) continue;
       if(!particle->GetPDG()) continue;
       Double_t lThisCharge = particle->GetPDG()->Charge()/3.0;
       // if(TMath::Abs(lThisCharge)<0.001) continue;
       if(! (lMCstack->IsPhysicalPrimary(ipart)) ) continue;

       //Out of bunch pileup event removal
       if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(ipart,lMCevent))
	 continue;
      

       Float_t trk_charge_gen = (particle->GetPDG()->Charge())/3.0;
       Float_t trk_PID_gen = particle->GetPdgCode();
       Float_t trk_eta_gen = particle->Eta();
       Float_t trk_pt_gen = particle->Pt();
       Float_t trk_VertexZ_gen = particle->Vz();


       //ACCEPTANCE CUTS : particles are rejected
       if(TMath::Abs(trk_VertexZ_gen)>fVertexZ_max)continue;
       if(TMath::Abs(trk_charge_gen)!=1)continue;
       if(TMath::Abs(trk_eta_gen)>fEta_max)continue;

       hist2D_pt_gen_centrality->Fill(trk_pt_gen,lV0M);

       //set1: 0.2<pT<2.0
       if(trk_pt_gen > 0.2 && trk_pt_gen <= 2.0)
	 {
	   //cout<<trk_charge_gen<<"\t"<<trk_VertexZ_gen<<"\t"<<trk_eta_gen<<"\t"<<trk_pt_gen<<endl;
	   Q1_gen[0]=Q1_gen[0]+TMath::Power(trk_pt_gen, 1.0);
	   Q2_gen[0]=Q2_gen[0]+TMath::Power(trk_pt_gen, 2.0);
	   Q3_gen[0]=Q3_gen[0]+TMath::Power(trk_pt_gen, 3.0);
	   Q4_gen[0]=Q4_gen[0]+TMath::Power(trk_pt_gen, 4.0);
	   Nch_gen[0]+=1.0;
	 }

       //set2: 0.2<pT<3.0
       if(trk_pt_gen > 0.2 && trk_pt_gen <= 3.0)
	 {
	   //cout<<trk_charge_gen<<"\t"<<trk_VertexZ_gen<<"\t"<<trk_eta_gen<<"\t"<<trk_pt_gen<<endl;
	   Q1_gen[1]=Q1_gen[1]+TMath::Power(trk_pt_gen, 1.0);
	   Q2_gen[1]=Q2_gen[1]+TMath::Power(trk_pt_gen, 2.0);
	   Q3_gen[1]=Q3_gen[1]+TMath::Power(trk_pt_gen, 3.0);
	   Q4_gen[1]=Q4_gen[1]+TMath::Power(trk_pt_gen, 4.0);
	   Nch_gen[1]+=1.0;

	 }

       //set3: 0.5<pT<2.0
       if(trk_pt_gen > 0.5 && trk_pt_gen <= 2.0)
	 {
	   //cout<<trk_charge_gen<<"\t"<<trk_VertexZ_gen<<"\t"<<trk_eta_gen<<"\t"<<trk_pt_gen<<endl;
	   Q1_gen[2]=Q1_gen[2]+TMath::Power(trk_pt_gen, 1.0);
	   Q2_gen[2]=Q2_gen[2]+TMath::Power(trk_pt_gen, 2.0);
	   Q3_gen[2]=Q3_gen[2]+TMath::Power(trk_pt_gen, 3.0);
	   Q4_gen[2]=Q4_gen[2]+TMath::Power(trk_pt_gen, 4.0);
	   Nch_gen[2]+=1.0;
	 }

       //set4: 0.5<pT<3.0
       if(trk_pt_gen > 0.5 && trk_pt_gen <= 3.0)
	 {
	   //cout<<trk_charge_gen<<"\t"<<trk_VertexZ_gen<<"\t"<<trk_eta_gen<<"\t"<<trk_pt_gen<<endl;
	   Q1_gen[3]=Q1_gen[3]+TMath::Power(trk_pt_gen, 1.0);
	   Q2_gen[3]=Q2_gen[3]+TMath::Power(trk_pt_gen, 2.0);
	   Q3_gen[3]=Q3_gen[3]+TMath::Power(trk_pt_gen, 3.0);
	   Q4_gen[3]=Q4_gen[3]+TMath::Power(trk_pt_gen, 4.0);
	   Nch_gen[3]+=1.0;
	 }
     }
   //end generated track loop

   for(int i1=0;i1<4;i1++)
     {
       fQ1_gen[i1]=Q1_gen[i1];
       fQ2_gen[i1]=Q2_gen[i1];
       fQ3_gen[i1]=Q3_gen[i1];
       fQ4_gen[i1]=Q4_gen[i1];
       fNch_gen[i1]=Nch_gen[i1];
     }
  
  


   ///////loop for reconstructed////////////////////////
   cout<<"*********************************** Reconstructed ********************************************"<<endl; 
   Int_t ntracks = lESDevent->GetNumberOfTracks();
    
   float Q1_rec[4]={0.0,0.0,0.0,0.0};
   float Q2_rec[4]={0.0,0.0,0.0,0.0};
   float Q3_rec[4]={0.0,0.0,0.0,0.0};
   float Q4_rec[4]={0.0,0.0,0.0,0.0};
   float Nch_rec[4]={0.0,0.0,0.0,0.0};


   for(Int_t itr = 0; itr < ntracks; itr++)
     {
       AliESDtrack *track = lESDevent->GetTrack(itr);
       if(!track)
	 {
	   AliWarning("ERROR: Could not retrieve one of the ESD tracks ...");
	   continue;
	 }
       if(fESDtrackCuts->AcceptTrack(track))
	 {

	   TParticle *particle_rec = lMCstack->Particle(TMath::Abs(track->GetLabel()));
      
	   Float_t trk_charge_rec = (particle_rec->GetPDG()->Charge())/3.0;
	   Float_t trk_PID_rec = particle_rec->GetPdgCode();
	   Float_t trk_eta_rec = particle_rec->Eta();
	   Float_t trk_pt_rec = particle_rec->Pt();
	   Float_t trk_VertexZ_rec = particle_rec->Vz();

	   //ACCEPTANCE CUTS : particles are rejected
	   //if(TMath::Abs(trk_VertexZ_rec)>10.0)continue;
	   if(TMath::Abs(trk_charge_rec)!=1)continue;
	   if(TMath::Abs(trk_eta_rec)>fEta_max)continue;
	  

	   hist2D_pt_rec_centrality->Fill(trk_pt_rec,lV0M);

	   //set1: 0.2<pT<2.0
	   if(trk_pt_rec > 0.2 && trk_pt_rec <= 2.0)
	     {
	       //cout<<trk_charge_rec<<"\t"<<trk_VertexZ_rec<<"\t"<<trk_eta_rec<<"\t"<<trk_pt_rec<<endl;
	       Q1_rec[0]=Q1_rec[0]+TMath::Power(trk_pt_rec, 1.0);
	       Q2_rec[0]=Q2_rec[0]+TMath::Power(trk_pt_rec, 2.0);
	       Q3_rec[0]=Q3_rec[0]+TMath::Power(trk_pt_rec, 3.0);
	       Q4_rec[0]=Q4_rec[0]+TMath::Power(trk_pt_rec, 4.0);
	       Nch_rec[0]+=1.0;
	     }

	   //set2: 0.2<pT<3.0
	   if(trk_pt_rec > 0.2 && trk_pt_rec <= 3.0)
	     {
	       //cout<<trk_charge_rec<<"\t"<<trk_VertexZ_rec<<"\t"<<trk_eta_rec<<"\t"<<trk_pt_rec<<endl;
	       Q1_rec[1]=Q1_rec[1]+TMath::Power(trk_pt_rec, 1.0);
	       Q2_rec[1]=Q2_rec[1]+TMath::Power(trk_pt_rec, 2.0);
	       Q3_rec[1]=Q3_rec[1]+TMath::Power(trk_pt_rec, 3.0);
	       Q4_rec[1]=Q4_rec[1]+TMath::Power(trk_pt_rec, 4.0);
	       Nch_rec[1]+=1.0;
	     }

	   //set3: 0.5<pT<2.0
	   if(trk_pt_rec > 0.5 && trk_pt_rec <= 2.0)
	     {
	       //cout<<trk_charge_rec<<"\t"<<trk_VertexZ_rec<<"\t"<<trk_eta_rec<<"\t"<<trk_pt_rec<<endl;
	       Q1_rec[2]=Q1_rec[2]+TMath::Power(trk_pt_rec, 1.0);
	       Q2_rec[2]=Q2_rec[2]+TMath::Power(trk_pt_rec, 2.0);
	       Q3_rec[2]=Q3_rec[2]+TMath::Power(trk_pt_rec, 3.0);
	       Q4_rec[2]=Q4_rec[2]+TMath::Power(trk_pt_rec, 4.0);
	       Nch_rec[2]+=1.0;
	     }

	   //set4: 0.5<pT<3.0
	   if(trk_pt_rec > 0.5 && trk_pt_rec <= 3.0)
	     {
	       //cout<<trk_charge_rec<<"\t"<<trk_VertexZ_rec<<"\t"<<trk_eta_rec<<"\t"<<trk_pt_rec<<endl;
	       Q1_rec[3]=Q1_rec[3]+TMath::Power(trk_pt_rec, 1.0);
	       Q2_rec[3]=Q2_rec[3]+TMath::Power(trk_pt_rec, 2.0);
	       Q3_rec[3]=Q3_rec[3]+TMath::Power(trk_pt_rec, 3.0);
	       Q4_rec[3]=Q4_rec[3]+TMath::Power(trk_pt_rec, 4.0);
	       Nch_rec[3]+=1.0;
	     }
      
	 }
     }//end reconstructed track loop

   for(int j1=0;j1<4;j1++)
     {
       fQ1_rec[j1]=Q1_rec[j1];
       fQ2_rec[j1]=Q2_rec[j1];
       fQ3_rec[j1]=Q3_rec[j1];
       fQ4_rec[j1]=Q4_rec[j1];
       fNch_rec[j1]=Nch_rec[j1];
     }


   fTreeEvent->Fill();
   // Post output data.
   PostData(1, fTreeEvent); 
   PostData(2, fListHist);
}
//________________________________________________________________________
void AliAnalysisTaskResonanceVsMultiplicityMC::Terminate(Option_t *)
{
  
}
//----------------------------------------------------------------------------
Double_t AliAnalysisTaskResonanceVsMultiplicityMC::MyRapidity(Double_t rE, Double_t rPz) const
{
  // Local calculation for rapidity
  Double_t ReturnValue = -100;
  if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ){ 
    ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
  }
  return ReturnValue;
} 
//----------------------------------------------------------------------------
void AliAnalysisTaskResonanceVsMultiplicityMC::SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts){
  if(fESDtrackCuts){ 
    delete fESDtrackCuts;
    fESDtrackCuts = 0x0;
  }
  fESDtrackCuts = trackcuts;
}
//----------------------------------------------------------------------------
Double_t AliAnalysisTaskResonanceVsMultiplicityMC::GetTOFBeta(AliVTrack *vtrack)
{
  AliESDtrack *esdtrack  = dynamic_cast<AliESDtrack*>(vtrack);
  if(!esdtrack) return -1;
  const Double_t c = 2.99792457999999984e-02; 
  Double_t p = esdtrack->GetTPCmomentum();
  Double_t l = esdtrack->GetIntegratedLength();
  Double_t trackT0 = fPIDResponse->GetTOFResponse().GetStartTime(p);
  Double_t timeTOF = esdtrack->GetTOFsignal()- trackT0;
  Double_t mass_square=(p*p)*(TMath::Power(c*timeTOF/l,2.0)-1);
  return mass_square;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskResonanceVsMultiplicityMC::MatchTOF(AliVTrack *vtrack)
{
  if (!vtrack) {
    AliWarning("NULL argument: impossible to check status");
    return kFALSE;
  }
  if (!(vtrack->GetStatus() & AliESDtrack::kTOFout)) return kFALSE;
  if (!(vtrack->GetStatus() & AliESDtrack::kTIME  )) return kFALSE;

  // if (!(vtrack->GetStatus() & AliESDtrack::kTOFpid)) return kFALSE;
  float probMis = fPIDResponse->GetTOFMismatchProbability(vtrack);
  if(probMis>0.01) return kFALSE;
  
  AliESDtrack *esdtrack  = dynamic_cast<AliESDtrack*>(vtrack);
  Double_t l = esdtrack->GetIntegratedLength();
  if(l<350) return kFALSE;
  
  return kTRUE;
}
