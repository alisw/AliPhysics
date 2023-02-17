
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
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskCorrPP.h"

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

class AliAnalysisTaskCorrPP;
ClassImp(AliAnalysisTaskCorrPP)

//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskCorrPP::AliAnalysisTaskCorrPP():
  AliAnalysisTaskSE(),
  fAODeventCuts(),
  fESDevent(0),
  fAODevent(0),
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
  hNumberOfCascades(0),
  hNumberOfXi(0),
  hNumberOfAntiXi(0),
  hNumberOfKaonPlus(0),
  hNumberOfKaonMinus(0),
  histMassXi_vs_Pt_beforeMasscut(0),
  histMassAntiXi_vs_Pt_beforeMasscut(0),
  histMassXi_vs_Pt(0),
  histMassAntiXi_vs_Pt(0),
  fTreeVariableCentrality(0),
  fNoXi_ptmax2(0),
  fNoAntiXi_ptmax2(0),
  fNoXi_ptmax3(0),
  fNoAntiXi_ptmax3(0),
  fNoKaonPlus_ptmax2(0),
  fNoKaonMinus_ptmax2(0),
  fNoKaonPlus_ptmax3(0),
  fNoKaonMinus_ptmax3(0)
{}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskCorrPP::AliAnalysisTaskCorrPP(const char *name):
  AliAnalysisTaskSE(name),
  fAODeventCuts(),
  fESDevent(0),
  fAODevent(0),
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
  hNumberOfCascades(0),
  hNumberOfXi(0),
  hNumberOfAntiXi(0),
  hNumberOfKaonPlus(0),
  hNumberOfKaonMinus(0),
  histMassXi_vs_Pt_beforeMasscut(0),
  histMassAntiXi_vs_Pt_beforeMasscut(0),
  histMassXi_vs_Pt(0),
  histMassAntiXi_vs_Pt(0),
  fTreeVariableCentrality(0),
  fNoXi_ptmax2(0),
  fNoAntiXi_ptmax2(0),
  fNoXi_ptmax3(0),
  fNoAntiXi_ptmax3(0),
  fNoKaonPlus_ptmax2(0),
  fNoKaonMinus_ptmax2(0),
  fNoKaonPlus_ptmax3(0),
  fNoKaonMinus_ptmax3(0)
{
  fUtils = new AliAnalysisUtils();
  DefineInput (0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TTree::Class());
}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskCorrPP::~AliAnalysisTaskCorrPP()  {

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
void AliAnalysisTaskCorrPP::UserCreateOutputObjects()  {
    
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

    //Cascade counter
    hNumberOfCascades = new TH1D("hNumberOfCascades","hNumberOfCascades",20,0,20);
    fOutputList -> Add(hNumberOfCascades);
   
    
    //Number of xi finally getting selected with specified cuts
    hNumberOfXi     = new TH1D ("hNumberOfXi","",20,0,20);
    hNumberOfAntiXi = new TH1D ("hNumberOfAntiXi","",20,0,20);
    // hNumberOfXi     -> Sumw2();
    // hNumberOfAntiXi -> Sumw2();
    fOutputList -> Add(hNumberOfXi);
    fOutputList -> Add(hNumberOfAntiXi);

    //Number of Kaon finally getting selected with specified cuts
    hNumberOfKaonPlus     = new TH1D ("hNumberOfKaonPlus","",100,0,100);
    hNumberOfKaonMinus = new TH1D ("hNumberOfKaonMinus","",100,0,100);
    // hNumberOfKaonPlus     -> Sumw2();
    // hNumberOfKaonMinus -> Sumw2();
    fOutputList -> Add(hNumberOfKaonPlus);
    fOutputList -> Add(hNumberOfKaonMinus);

    //Invariant-Mass Plots if needed
    histMassXi_vs_Pt_beforeMasscut     = new TH2D ("histMassXi_vs_Pt_beforeMasscut","histMassXi_vs_Pt_beforeMasscut",100,0,10,500,1.30,1.35);
    histMassAntiXi_vs_Pt_beforeMasscut = new TH2D ("histMassAntiXi_vs_Pt_beforeMasscut","histMassAntiXi_vs_Pt_beforeMasscut",100,0,10,500,1.30,1.35);
    histMassXi_vs_Pt     = new TH2D ("histMassXi_vs_Pt","histMassXi_vs_Pt",100,0,10,500,1.30,1.35);
    histMassAntiXi_vs_Pt = new TH2D ("histMassAntiXi_vs_Pt","histMassAntiXi_vs_Pt",100,0,10,500,1.30,1.35);
    // hMassXi_vs_P     -> Sumw2();
    // hMassAntiXi_vs_P -> Sumw2();
    fOutputList -> Add(histMassXi_vs_Pt_beforeMasscut);
    fOutputList -> Add(histMassAntiXi_vs_Pt_beforeMasscut); 
    fOutputList -> Add(histMassXi_vs_Pt);
    fOutputList -> Add(histMassAntiXi_vs_Pt);

    //TTree object to store variables
    fTreeEvent = new TTree("fTreeEvent","Event Tree");
    fTreeEvent->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
    fTreeEvent->Branch("fNoXi_ptmax2", &fNoXi_ptmax2, "fNoXi_ptmax2/F");
    fTreeEvent->Branch("fNoAntiXi_ptmax2", &fNoAntiXi_ptmax2, "fNoAntiXi_ptmax2/F");
    fTreeEvent->Branch("fNoXi_ptmax3", &fNoXi_ptmax3, "fNoXi_ptmax3/F");
    fTreeEvent->Branch("fNoAntiXi_ptmax3", &fNoAntiXi_ptmax3, "fNoAntiXi_ptmax3/F");
    fTreeEvent->Branch("fNoKaonPlus_ptmax2", &fNoKaonPlus_ptmax2, "fNoKaonPlus_ptmax2/F");
    fTreeEvent->Branch("fNoKaonMinus_ptmax2", &fNoKaonMinus_ptmax2, "fNoKaonMinus_ptmax2/F");
    fTreeEvent->Branch("fNoKaonPlus_ptmax3", &fNoKaonPlus_ptmax3, "fNoKaonPlus_ptmax3/F");
    fTreeEvent->Branch("fNoKaonMinus_ptmax3", &fNoKaonMinus_ptmax3, "fNoKaonMinus_ptmax3/F");
   
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
void AliAnalysisTaskCorrPP::UserExec(Option_t *)  {
  
    //Get Input Event
    if ( !GetEvent ()) return;
    cout<<"*********************** Found event !!! ******************************"<<endl;

    
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

   
    //Initialize number 0f K+, K- per event
    Int_t no_KaonPlus_perevent = 0;
    Int_t no_KaonMinus_perevent = 0;
    Int_t no_KaonPlus_perevent_ptmax2 = 0;
    Int_t no_KaonMinus_perevent_ptmax2 = 0;

    //Loop over primary tracks to select kaon
    for( Int_t itr = 0; itr < fAODevent->GetNumberOfTracks(); itr++)
      {

	AliVTrack   *track = (AliVTrack*)fAODevent->GetTrack(itr);
	if(!track)      continue;
	AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(track);
	if(!aodtrack)      continue;

	
	if(!aodtrack->TestFilterBit(96))  continue;

	
	Double_t trkPt = aodtrack->Pt();
	Double_t trkPhi = aodtrack->Phi();
	Double_t trkEta = aodtrack->Eta();
	Double_t trkCharge = aodtrack->Charge();
	Double_t trkDCAxy = aodtrack->DCA();
	Double_t trkDCAz = aodtrack->ZAtDCA();
	Double_t trkTPCNCls = aodtrack->GetTPCNcls();
	Double_t trkChi2PerNDF = aodtrack->Chi2perNDF();
	

	Double_t p[3];
	track->PxPyPz(p);
	Double_t Track_pt;
	Track_pt=TMath::Sqrt((p[0]*p[0])+(p[1]*p[1]));

	if (Track_pt < 0.3) continue;
	if (Track_pt > 3.0) continue;
	if (TMath::Abs(trkEta) > 0.8) continue;
	
	Bool_t IsKaon = KaonSelector(track);
	if (!IsKaon) continue;

	cout<<"Found Kaon particle !!!"<<endl;

	if (track->Charge() > 0)   //K+
	  {
	    no_KaonPlus_perevent += 1;
	    if (Track_pt < 2.0)
	      no_KaonPlus_perevent_ptmax2 += 1;
	  }

	if (track->Charge() < 0)   //K-
	  {
	    no_KaonMinus_perevent += 1;
	    if (Track_pt < 2.0)
	      no_KaonMinus_perevent_ptmax2 += 1;
	  }
      }
      //end primary track loop
    

    //Primary Vertex
    AliAODVertex *primaryVertex = (AliAODVertex*) fAODevent->GetPrimaryVertex();
    Double_t vx = primaryVertex->GetX();
    Double_t vy = primaryVertex->GetY();
    Double_t vz = primaryVertex->GetZ();
    Double_t prim_vertex[3] = {-100, -100, -100};
    primaryVertex->GetXYZ(prim_vertex);
   
    //Initialize number 0f Xi, anti-Xi per event
    Int_t no_Xi_perevent = 0;
    Int_t no_antiXi_perevent = 0;
    Int_t no_Xi_perevent_ptmax2 = 0;
    Int_t no_antiXi_perevent_ptmax2 = 0;
    Int_t no_Xi_perevent_ptmax3 = 0;
    Int_t no_antiXi_perevent_ptmax3 = 0;

    
    //Loop over Reconstructed Cascades
    for (Int_t icasc=0 ; icasc<fAODevent->GetNumberOfCascades(); icasc++)  {
        
      //Get Reconstructed Cascade
      
      AliAODcascade *cascade = fAODevent->GetCascade(icasc);
      if (!cascade)
	{
	  continue;
	}
      //cout<<"*****************Cascade particle found !!!**********************************"<<endl;

       hNumberOfCascades->Fill(0.5);
	    
        
      //Get Decay Daughters
      AliAODTrack *positive_track = dynamic_cast<AliAODTrack*>( cascade->GetDaughter(0) );
      AliAODTrack *negative_track = dynamic_cast<AliAODTrack*>( cascade->GetDaughter(1) );
      AliAODTrack *bachelor_track = dynamic_cast<AliAODTrack*>( cascade->GetDecayVertexXi()->GetDaughter(0) );
      if (!positive_track) continue;
      if (!negative_track) continue;
      if (!bachelor_track) continue;

      cout<<"**************Daughter tracks of Cascade particle found !!!**************************"<<endl;

      hNumberOfCascades->Fill(1.5);
       
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //Track properties
      Int_t nTPCcluster_posTrk = positive_track->GetTPCNcls();
      Int_t nTPCcluster_negTrk = negative_track->GetTPCNcls();
      Int_t nTPCcluster_bacTrk = bachelor_track->GetTPCNcls();
      cout<<"TPC cluster: "<<nTPCcluster_posTrk<<"\t"<<nTPCcluster_negTrk<<"\t"<<nTPCcluster_bacTrk<<endl;

      Double_t chi2perTPCcluster_posTrk = positive_track->GetTPCchi2perCluster();
      Double_t chi2perTPCcluster_negTrk = negative_track->GetTPCchi2perCluster();
      Double_t chi2perTPCcluster_bacTrk = bachelor_track->GetTPCchi2perCluster();
      cout<<"Chi2 per TPC cluster: "<<chi2perTPCcluster_posTrk<<"\t"<<chi2perTPCcluster_negTrk<<"\t"<<chi2perTPCcluster_bacTrk<<endl;

      Double_t Eta_posTrk = positive_track->Eta();
      Double_t Eta_negTrk = negative_track->Eta();
      Double_t Eta_bacTrk = bachelor_track->Eta();
      cout<<"Eta: "<<Eta_posTrk<<"\t"<<Eta_negTrk<<"\t"<<Eta_bacTrk<<endl;

      Double_t Pt_posTrk = positive_track->Pt();
      Double_t Pt_negTrk = negative_track->Pt();
      Double_t Pt_bacTrk = bachelor_track->Pt();
      cout<<"Pt: "<<Pt_posTrk<<"\t"<<Pt_negTrk<<"\t"<<Pt_bacTrk<<endl;
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
      
      //Track Quality Cuts
      Bool_t Cut1 = kFALSE;
      if (!PassedTrackQualityCuts(positive_track)) continue;
      if (!PassedTrackQualityCuts(negative_track)) continue;
      if (!PassedTrackQualityCuts(bachelor_track)) continue;
      cout<<"Daughter tracks of cascade passed track quality cuts !!!!!"<<endl;

       hNumberOfCascades->Fill(2.5);
      
      //Track pile up cuts
      if (!PassedSingleParticlePileUpCuts(positive_track)) continue;
      if (!PassedSingleParticlePileUpCuts(negative_track)) continue;
      if (!PassedSingleParticlePileUpCuts(bachelor_track)) continue;
      cout<<"Passed pileup cuts for track !!!"<<endl;

      hNumberOfCascades->Fill(3.5);

      //++++++++++++++++++Topological variables+++++++++++++++++++++++++++++++++++
      //DCA of daughter tracks
      
      Double_t lDcaBachTrackToPrimVertexXi = cascade->DcaBachToPrimVertex();
      Double_t lDcaPosTrackToPrimVertexXi = cascade->DcaPosToPrimVertex();
      Double_t lDcaNegTrackToPrimVertexXi = cascade->DcaNegToPrimVertex();
      cout<<"DCAtoPV_bachelor: "<<lDcaBachTrackToPrimVertexXi<<endl;
      cout<<"DCAtoPV_positive: "<<lDcaPosTrackToPrimVertexXi<<endl;
      cout<<"DCAtoPV_negative: "<<lDcaNegTrackToPrimVertexXi<<endl;


      //Variables of V0
      
      Double_t lV0CosPointingAngleXi = cascade->CosPointingAngle(prim_vertex);
      Double_t lDcaV0DaughtersXi = cascade->DcaV0Daughters();
      Double_t lDcaV0ParticleToPrimVertexXi = cascade->DcaV0ToPrimVertex();
      //V0 Radius
      Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
      cascade->GetXYZ( lPosV0Xi );
      Double_t lV0TransRadiusXi = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
      cout<<"cosPointingAngleV0: "<<lV0CosPointingAngleXi<<endl;
      cout<<"DCA_V0_daughter: "<<lDcaV0DaughtersXi<<endl;
      cout<<"DCA_V0_toPV: "<<lDcaV0ParticleToPrimVertexXi<<endl;
      cout<<"V0 trans radius: "<<lV0TransRadiusXi<<endl;

      //Pt of V0
      Double_t lBMom[3], lNMom[3], lPMom[3];
      bachelor_track->GetPxPyPz( lBMom );
      positive_track->GetPxPyPz( lPMom );
      negative_track->GetPxPyPz( lNMom );
      Float_t lV0Pt = TMath::Sqrt(  TMath::Power( lNMom[0]+lPMom[0] , 2)
                                  + TMath::Power( lNMom[1]+lPMom[1] , 2) );
       

      // Variables of Xi
	
      Double_t lCosPointingAngleXi = cascade->CosPointingAngleXi(vx,vy,vz);
      Double_t lDcaXiDaughters = cascade->DcaXiDaughters();
      //Cascade Radius
      Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
      lPosXi[0] = cascade->DecayVertexXiX();
      lPosXi[1] = cascade->DecayVertexXiY();
      lPosXi[2] = cascade->DecayVertexXiZ();
      Double_t lTransRaiusXi = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
      cout<<"cosPointingAngleXi: "<<lCosPointingAngleXi<<endl;
      cout<<"DCA_Xidaughters: "<<lDcaXiDaughters<<endl;
      cout<<"Xi Trans Radius: "<< lTransRaiusXi <<endl;
      
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      //DCA cut of daughters to PV
      if (!PassedDaughterTrackDCAtoVertexSelectionCuts(cascade)) continue;
      cout<<"**************Daughter tracks passed DCA to PV cut !!!###############################"<<endl;

       hNumberOfCascades->Fill(4.5);

      //Extra: ptcut on Lambda
      if (lV0Pt < 0.3) continue;

      hNumberOfCascades->Fill(5.5);
     

      //Topological cuts to select cascade particle
      if (!PassedCascadeSelectionCuts(cascade)) continue;
      cout<<"Passed CascadeSelections cut !!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
     
      //Mass and Charge
      cout<<"Mass of cascade: "<<cascade->MassXi()<<"\t"<<"Charge: "<<cascade->ChargeXi()<<endl;
      
      	
      Bool_t isXi=(kFALSE),isAntiXi=(kFALSE);
      if (IsXiCandidate(cascade,positive_track,negative_track,bachelor_track))
	isXi     = kTRUE;
      if (IsAntiXiCandidate(cascade,positive_track,negative_track,bachelor_track))
	isAntiXi = kTRUE;
      if ((!isXi)&&(!isAntiXi)) continue;

      Double_t massXi = cascade->MassXi();
      Double_t ptXi = cascade->Pt();

      //Xi Selection.....filling invariant mass as function of momentum for xi
      if (isXi)
	{
	  histMassXi_vs_Pt->Fill (ptXi,massXi);
	  no_Xi_perevent+=1;

	  if(ptXi < 2.0)
	    no_Xi_perevent_ptmax2 +=1;
	  if(ptXi < 3.0)
	    no_Xi_perevent_ptmax3 +=1;

	  cout<<"************** Found Xi particle !!! *******:"<<no_Xi_perevent<<endl;
	}
      //AntiXi Selection
      if (isAntiXi)
	{
	  histMassAntiXi_vs_Pt->Fill (ptXi,massXi);
	  no_antiXi_perevent+=1;

	  if(ptXi < 2.0)
	    no_antiXi_perevent_ptmax2 +=1;
	  if(ptXi < 3.0)
	    no_antiXi_perevent_ptmax3 +=1;

	  cout<<"************** Found AntiXi particle !!! *******:"<<no_antiXi_perevent<<endl;
	}
      if((isXi) || (isAntiXi))
      cout<<"**** Found Xi particle !!! *******"<<endl;
      
    }
    //end cascade loop

    
    // cout<<"************** Xi ****************: "<<no_Xi_perevent<<endl;
    // cout<<"************** AntiXi ****************: "<<no_antiXi_perevent<<endl;
     

    hNumberOfXi->Fill(no_Xi_perevent);
    hNumberOfAntiXi->Fill(no_antiXi_perevent);
    hNumberOfKaonPlus->Fill(no_KaonPlus_perevent);
    hNumberOfKaonMinus->Fill(no_KaonMinus_perevent);

     
    //"fTreeEvent" Tree Variable
    fTreeVariableCentrality = lV0M;
    fNoXi_ptmax2 = no_Xi_perevent_ptmax2 ;
    fNoAntiXi_ptmax2 = no_antiXi_perevent_ptmax2 ;
    fNoXi_ptmax3 = no_Xi_perevent_ptmax3 ;
    fNoAntiXi_ptmax3 = no_antiXi_perevent_ptmax3 ;
    fNoKaonPlus_ptmax2 = no_KaonPlus_perevent_ptmax2 ;
    fNoKaonMinus_ptmax2 = no_KaonMinus_perevent_ptmax2 ;
    fNoKaonPlus_ptmax3 = no_KaonPlus_perevent ;
    fNoKaonMinus_ptmax3 = no_KaonMinus_perevent ;
 
    fTreeEvent->Fill();
  
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeEvent);
}    
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPP::GetEvent ()  //event cuts copied from my code written earlier 

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
Bool_t AliAnalysisTaskCorrPP::PassedTrackQualityCuts (AliAODTrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);

    //select global tracks: ITS+TPC loose DCA
    // if(!track->TestFilterMask(BIT(4)))
    //   return passedTrkSelection;
    
    //Track Selection Cuts
    Int_t nTPCcluster = track->GetTPCNcls();  //TPC cluster cut
    if (nTPCcluster <= 70)
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
    if(pt < 0.2)
      return passedTrkSelection;

    
    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
 //______________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPP::PassedDaughterTrackDCAtoVertexSelectionCuts(AliAODcascade *cascade)
{
  //Initialization
  Bool_t passedDauTrackDCAtoVertexSelection=(kFALSE);

  Double_t lDcaBachTrackToPrimVertexXi = -1.;
  Double_t lDcaPosTrackToPrimVertexXi  = -1.;
  Double_t lDcaNegTrackToPrimVertexXi  = -1.;

  //DCAs for cuts
  lDcaBachTrackToPrimVertexXi = cascade->DcaBachToPrimVertex();
  lDcaPosTrackToPrimVertexXi = cascade->DcaPosToPrimVertex();
  lDcaNegTrackToPrimVertexXi = cascade->DcaNegToPrimVertex();
  
  if (lDcaBachTrackToPrimVertexXi <= 0.06) return passedDauTrackDCAtoVertexSelection;
  if (lDcaPosTrackToPrimVertexXi <= 0.06) return passedDauTrackDCAtoVertexSelection;
  if (lDcaNegTrackToPrimVertexXi <= 0.06) return passedDauTrackDCAtoVertexSelection;
  
  passedDauTrackDCAtoVertexSelection = kTRUE;
  return passedDauTrackDCAtoVertexSelection;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPP::PassedCascadeSelectionCuts (AliAODcascade *cascade)  {

    //Initialization
    Bool_t passedCascadeSelection=(kFALSE);

    //Primary Vertex
    AliAODVertex *primaryVertex = (AliAODVertex*) fAODevent->GetPrimaryVertex();
    Double_t vx = primaryVertex->GetX();
    Double_t vy = primaryVertex->GetY();
    Double_t vz = primaryVertex->GetZ();

    Double_t primVtx[3] = {-100.0, -100.0, -100.0};
    primaryVertex->GetXYZ(primVtx);

    //Topological Selections for Lambda   //taken from a paper searched in google....please check the cut values

    if (cascade->CosPointingAngle(primVtx) <= 0.979) return passedCascadeSelection;

    hNumberOfCascades->Fill(6.5);
      
    if (cascade->DcaV0Daughters() >= 1.0) return passedCascadeSelection;

    hNumberOfCascades->Fill(7.5);
      

    Double_t lDcaV0ParticleToPrimVertexXi = -1.;
    lDcaV0ParticleToPrimVertexXi = cascade->DcaV0ToPrimVertex();
    if (lDcaV0ParticleToPrimVertexXi <= 0.08) return passedCascadeSelection;

    hNumberOfCascades->Fill(8.5);
      

    //V0 Radius
    Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
    cascade->GetXYZ( lPosV0Xi );
    Double_t radiusV0 = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
    Double_t rV01=1.7, rV02=200;
    if (radiusV0 <= rV01) return passedCascadeSelection;
    if (radiusV0 >= rV02) return passedCascadeSelection;

    hNumberOfCascades->Fill(9.5);

    //Topological Selections for Xi 
    if (cascade->CosPointingAngleXi(vx,vy,vz) <= 0.98) return passedCascadeSelection;

    hNumberOfCascades->Fill(10.5);
      
    if (cascade->DcaXiDaughters() >= 1.5) return passedCascadeSelection;

    hNumberOfCascades->Fill(11.5);
      
  //Cascade Radius
    Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
    lPosXi[0] = cascade->DecayVertexXiX();
    lPosXi[1] = cascade->DecayVertexXiY();
    lPosXi[2] = cascade->DecayVertexXiZ();
    Double_t rC = TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
    Double_t r01=0.8, r02=200;  //cascade radius to be seen from note
    if (rC <= r01) return passedCascadeSelection;
    if (rC >= r02) return passedCascadeSelection;

    hNumberOfCascades->Fill(12.5);
      
    passedCascadeSelection = kTRUE;
    return passedCascadeSelection;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPP::IsXiCandidate (AliAODcascade *casc, AliAODTrack *pos, AliAODTrack *neg, AliAODTrack *bac)  {
    
    //PID Daughters
    Bool_t passedPID=(kFALSE);
    if (PassedPIDSelection(pos,AliPID::kProton) && PassedPIDSelection(neg,AliPID::kPion) && PassedPIDSelection(bac,AliPID::kPion)) passedPID=kTRUE;
    if (!passedPID) return kFALSE;

    hNumberOfCascades->Fill(13.5);
   
    //Mass of V0(daughter of Xi)
    Double_t massV0 = casc->MassLambda();   

    Double_t massLambda_PDG=TDatabasePDG::Instance()->GetParticle(3122)->Mass();
    if (massV0 <= massLambda_PDG-0.006) return kFALSE;
    if (massV0 >= massLambda_PDG+0.006) return kFALSE;

    hNumberOfCascades->Fill(14.5);

    //Mass of Xi if bachelor is misidentified as Kaon, Omega rejection
    Double_t massOmega = casc->MassOmega();
    if (massOmega > 1.667 && massOmega < 1.677)
      return kFALSE;

    hNumberOfCascades->Fill(15.5);
 

    //Fill histogram before masscut
    histMassXi_vs_Pt_beforeMasscut->Fill(casc->Pt(),casc->MassXi());
    
    //Mass of Xi
    Double_t mass = casc->MassXi();   
 
    //Mass Selection
    Double_t massXi_PDG=TDatabasePDG::Instance()->GetParticle(3312)->Mass();
    if (mass <= massXi_PDG-0.005) return kFALSE;
    if (mass >= massXi_PDG+0.005) return kFALSE;

    hNumberOfCascades->Fill(16.5);
    
    //Assignments
    // m = mass;
    // momentum = Ptot.Vect();

    return kTRUE;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPP::IsAntiXiCandidate (AliAODcascade *casc, AliAODTrack *pos, AliAODTrack *neg, AliAODTrack *bac)  {
    
    //PID Daughters
    Bool_t passedPID=(kFALSE);
    if (PassedPIDSelection(pos,AliPID::kPion) && PassedPIDSelection(neg,AliPID::kProton) && PassedPIDSelection(bac,AliPID::kPion)) passedPID=kTRUE;
    if (!passedPID) return kFALSE;

    hNumberOfCascades->Fill(13.5);
    
    //Mass of V0(daughter of Xi)
    Double_t massV0 = casc->MassAntiLambda(); 
    
    Double_t massLambda_PDG=TDatabasePDG::Instance()->GetParticle(-3122)->Mass();
    if (massV0 <= massLambda_PDG-0.006) return kFALSE;
    if (massV0 >= massLambda_PDG+0.006) return kFALSE;

    hNumberOfCascades->Fill(14.5);

    //Mass of Xi if bachelor is misidentified as Kaon, Omega rejection
    Double_t massOmega = casc->MassOmega();
    if (massOmega > 1.667 && massOmega < 1.677)
      return kFALSE;

    hNumberOfCascades->Fill(15.5);
 

    //Fill histogram before masscut
    histMassAntiXi_vs_Pt_beforeMasscut->Fill(casc->Pt(),casc->MassXi());

    
    //Mass of Xi
    Double_t mass = casc->MassXi();

    //Mass Selection
    Double_t massXi_PDG=TDatabasePDG::Instance()->GetParticle(-3312)->Mass();
    if (mass <= massXi_PDG-0.005) return kFALSE;
    if (mass >= massXi_PDG+0.005) return kFALSE;

    hNumberOfCascades->Fill(16.5);
    
    //Assignments
    // m = mass;
    // momentum = Ptot.Vect();

    return kTRUE;
}
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPP::KaonSelector(AliVTrack *track)  {
  /*
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

  Int_t flag = 0;
  if (TMath::Abs(fTPCnSigmaPion) < 3.0)
    flag += 1;
  if (TMath::Abs(fTPCnSigmaProton) < 3.0)
    flag += 1;
  if (TMath::Abs(fTPCnSigmaKaon) < 3.0)
    flag += 1;

  if (flag > 1) return kFALSE;


  //Selection for pT < 0.5 : TPC only
  if( Track_pt < 0.5 )
    {
      if(TMath::Abs(fTPCnSigmaKaon) < 2.0)
	return kTRUE;
      else
	return kFALSE;
    }

  //Selection for pT > 0.5 : TOF + TPC
  if( Track_pt >= 0.5 )
    {
      if (fTPCplusTOFnSigmaKaon < 2.0)
	return kTRUE;
      else
	return kFALSE;
    }
  */

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
Bool_t AliAnalysisTaskCorrPP::PassedPIDSelection (AliAODTrack *track, AliPID::EParticleType type)  {
    
    //Initialization
    Bool_t passedPIDSelection=(kFALSE);
    
    //TPC Particle Identification
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,type);
    if (nsigmaTPC <= -4.0) return passedPIDSelection;
    if (nsigmaTPC >= +4.0) return passedPIDSelection;

    passedPIDSelection = kTRUE;
    return passedPIDSelection;
}
 //_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskCorrPP::PassedSingleParticlePileUpCuts(AliAODTrack *track)
{
  Bool_t passedTrackPileupCut = (kFALSE);

  /*
  if (!(track->HasPointOnITSLayer(1)) && !(track->HasPointOnITSLayer(4)) && !(track->HasPointOnITSLayer(5)) && !(track->GetTOFBunchCrossing() == 0))
    passedTrackPileupCut = (kFALSE);
  */

  Int_t flag = 0;
  if(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))
    flag+=1;
  if(track->HasPointOnITSLayer(4) || track->HasPointOnITSLayer(5))
    flag+=1;
  if(track->GetTOFBunchCrossing() == 0)
    flag+=1;

  if(flag != 0)
    passedTrackPileupCut = (kTRUE);
  
  return passedTrackPileupCut;
}
 //_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskCorrPP::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
