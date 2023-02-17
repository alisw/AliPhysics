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

#include "AliAnalysisTaskResonanceVsMultiplicity_MCnoESD.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"



using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskResonanceVsMultiplicity_MCnoESD)

AliAnalysisTaskResonanceVsMultiplicity_MCnoESD::AliAnalysisTaskResonanceVsMultiplicity_MCnoESD(): AliAnalysisTaskSE(), fTreeEvent(0), fListHist(0), fPIDResponse(0), fESDtrackCuts(0), fEventCuts(0), fTreeVariableCentrality(0), fvertex(0), Profile_mean_term1(0),Profile_var_term1(0),Profile_var_term2(0),Profile_skewness_term1(0),Profile_skewness_term2(0),Profile_skewness_term3(0),Profile_kurtosis_term1(0),Profile_kurtosis_term2(0),Profile_kurtosis_term3(0),Profile_kurtosis_term4(0), hist2D_pt_gen_centrality(0), hist2D_pt_rec_centrality(0), hist_centrality_beforecut(0), fNch_eta0pt5(0)
{
  for(int i=0;i<2;i++)
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

AliAnalysisTaskResonanceVsMultiplicity_MCnoESD::AliAnalysisTaskResonanceVsMultiplicity_MCnoESD(const char *name): AliAnalysisTaskSE(name), fTreeEvent(0), fListHist(0), fPIDResponse(0), fESDtrackCuts(0), fEventCuts(0), fTreeVariableCentrality(0), fvertex(0),Profile_mean_term1(0),Profile_var_term1(0),Profile_var_term2(0),Profile_skewness_term1(0),Profile_skewness_term2(0),Profile_skewness_term3(0),Profile_kurtosis_term1(0),Profile_kurtosis_term2(0),Profile_kurtosis_term3(0),Profile_kurtosis_term4(0), hist2D_pt_gen_centrality(0), hist2D_pt_rec_centrality(0), hist_centrality_beforecut(0), fNch_eta0pt5(0)
{
  for(int i=0;i<2;i++)
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

AliAnalysisTaskResonanceVsMultiplicity_MCnoESD::~AliAnalysisTaskResonanceVsMultiplicity_MCnoESD()
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
void AliAnalysisTaskResonanceVsMultiplicity_MCnoESD::UserCreateOutputObjects()
{

  //------------------------------------------------
  // Particle Identification Setup
  //------------------------------------------------
  //AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  //AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  //fPIDResponse = inputHandler->GetPIDResponse();
  
   //------------------------------------------------
  // track cut
  //------------------------------------------------

  
  OpenFile(1);
  fTreeEvent = new TTree("fTreeEvent","Event");
  fTreeEvent->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
  //fTreeEvent->Branch("fvertex",&fvertex,"fvertex/F");
  fTreeEvent->Branch("fNch_eta0pt5",&fNch_eta0pt5,"fNch_eta0pt5/F");
  fTreeEvent->Branch("fNch_gen", &fNch_gen, "fNch_gen[2]/F");
  fTreeEvent->Branch("fQ1_gen", &fQ1_gen, "fQ1_gen[2]/F");
  fTreeEvent->Branch("fQ2_gen", &fQ2_gen, "fQ2_gen[2]/F");
  fTreeEvent->Branch("fQ3_gen", &fQ3_gen, "fQ3_gen[2]/F");
  fTreeEvent->Branch("fQ4_gen", &fQ4_gen, "fQ4_gen[2]/F");
  // fTreeEvent->Branch("fNch_rec", &fNch_rec, "fNch_rec[2]/F");
  // fTreeEvent->Branch("fQ1_rec", &fQ1_rec, "fQ1_rec[2]/F");
  // fTreeEvent->Branch("fQ2_rec", &fQ2_rec, "fQ2_rec[2]/F");
  // fTreeEvent->Branch("fQ3_rec", &fQ3_rec, "fQ3_rec[2]/F");
  // fTreeEvent->Branch("fQ4_rec", &fQ4_rec, "fQ4_rec[2]/F");
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
void AliAnalysisTaskResonanceVsMultiplicity_MCnoESD::UserExec(Option_t *)
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
    Float_t fMCImpactParameter=0.0,fCentImpBin=-10.0;
    

  
    lHIJINGHeader = (AliGenHijingEventHeader*) genHeader;
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

      AliMCParticle *trackMCgen;

      Int_t trackcharge=0;
      float Track_pt=0.0;
      float_t Track_eta=0.0;
      float Q1[2]={0.0,0.0};
      float Q2[2]={0.0,0.0};
      float Q3[2]={0.0,0.0};
      float Q4[2]={0.0,0.0};
      float Nch[2]={0.0,0.0};
      float Nch_eta0pt5=0.0;

      if(fCentImpBin < 0.0)
	{
	  PostData(1,fTreeEvent);
	  return;
	}
      

      for (Int_t iTrack = 0; iTrack < lMCevent->GetNumberOfTracks(); iTrack++)
	{
	  trackMCgen = (AliMCParticle *)lMCevent->GetTrack(iTrack);
	  if (!trackMCgen) continue;
	  // 
	  // Accepty only primary particles
	  if (!lMCstack->IsPhysicalPrimary(iTrack)) continue;
	  //
	  //cout<<iTrack<<"\t"<<"Pt: "<<trackMCgen->Pt()<<"\t"<<"Eta: "<<trackMCgen->Eta()<<"\t"<<"Charge: "<<trackMCgen->Charge()<<endl;

	  trackcharge=1.0/3.0*trackMCgen->Charge();
	  Track_pt=trackMCgen->Pt();
	  Track_eta=trackMCgen->Eta();

	  if(TMath::Abs(Track_eta)>0.8)continue;
	  if(TMath::Abs(trackcharge)!=1)continue;

	  if(TMath::Abs(Track_eta)<0.5)
	    {
	      Nch_eta0pt5+=1.0;
	    }
      
	  if(Track_pt>0.2 && Track_pt<=2.0)
	    {
	      Q1[0]=Q1[0]+TMath::Power(Track_pt, 1.0);
	      Q2[0]=Q2[0]+TMath::Power(Track_pt, 2.0);
	      Q3[0]=Q3[0]+TMath::Power(Track_pt, 3.0);
	      Q4[0]=Q4[0]+TMath::Power(Track_pt, 4.0);
	      Nch[0]=Nch[0]+1;
	    }
	  if(Track_pt>0.2 && Track_pt<=3.0)
	    {
	      Q1[1]=Q1[1]+TMath::Power(Track_pt, 1.0);
	      Q2[1]=Q2[1]+TMath::Power(Track_pt, 2.0);
	      Q3[1]=Q3[1]+TMath::Power(Track_pt, 3.0);
	      Q4[1]=Q4[1]+TMath::Power(Track_pt, 4.0);
	      Nch[1]=Nch[1]+1;
	    }
	}

  

      cout<<"Value of centrality:"<<fCentImpBin<<endl;
      fTreeVariableCentrality=fCentImpBin;
      fNch_eta0pt5=Nch_eta0pt5;
      fQ1_gen[0]=Q1[0];
      fQ2_gen[0]=Q2[0];
      fQ3_gen[0]=Q3[0];
      fQ4_gen[0]=Q4[0];
      fNch_gen[0]=Nch[0];
      fQ1_gen[1]=Q1[1];
      fQ2_gen[1]=Q2[1];
      fQ3_gen[1]=Q3[1];
      fQ4_gen[1]=Q4[1];
      fNch_gen[1]=Nch[1];
  
      fTreeEvent->Fill();
      PostData(1,fTreeEvent);
      

}
//________________________________________________________________________
void AliAnalysisTaskResonanceVsMultiplicity_MCnoESD::Terminate(Option_t *)
{
  
}
//----------------------------------------------------------------------------
