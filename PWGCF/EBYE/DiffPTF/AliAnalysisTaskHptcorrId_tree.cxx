
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
#include "AliAnalysisTaskHptcorrId_tree.h"

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

class AliAnalysisTaskHptcorrId_tree;
ClassImp(AliAnalysisTaskHptcorrId_tree)

//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskHptcorrId_tree::AliAnalysisTaskHptcorrId_tree():
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
  fTreeVariableCentrality(0),
  fQ1_prot(0),
  fQ2_prot(0),
  fQ3_prot(0),
  fQ4_prot(0),
  fNch_prot(0),
  fQ1_pi(0),
  fQ2_pi(0),
  fQ3_pi(0),
  fQ4_pi(0),
  fNch_pi(0),
  fQ1_K(0),
  fQ2_K(0),
  fQ3_K(0),
  fQ4_K(0),
  fNch_K(0),
  fP1_prot(0),
  fP2_prot(0),
  fP3_prot(0),
  fP4_prot(0),
  fP1_pi(0),
  fP2_pi(0),
  fP3_pi(0),
  fP4_pi(0),
  fP1_K(0),
  fP2_K(0),
  fP3_K(0),
  fP4_K(0),
  fW1_prot(0),
  fW2_prot(0),
  fW3_prot(0),
  fW4_prot(0),
  fW1_pi(0),
  fW2_pi(0),
  fW3_pi(0),
  fW4_pi(0),
  fW1_K(0),
  fW2_K(0),
  fW3_K(0),
  fW4_K(0),
  fPtMax(0),
  fPtMin(0),
  fListTRKCorr(0), 
  fHistMCEffKaonPlus(0),
  fHistMCEffKaonMinus(0),
  fHistMCEffPionPlus(0),
  fHistMCEffPionMinus(0),
  fHistMCEffProtonPlus(0),
  fHistMCEffProtonMinus(0),
  fVertexZMax(0),
  fFBNo(0),
  fChi2TPC(0),
  fChi2ITS(0),
  fPIDnSigmaPionCut(0),
  fPIDnSigmaKaonCut(0),
  fPIDnSigmaProtonCut(0),
  fTPCcrossedrows(0),
  fEtaMax(0),
  fPileupCutVal(0),
  fCentralityEstimator_flag(0),
  fTreeName(0)
{
  for(int i=0; i<9; i++)
    {
      fEffProtonPlus[i] = NULL;
      fEffProtonMinus[i] = NULL;
      fEffPionPlus[i] = NULL;
      fEffPionMinus[i] = NULL;
      fEffKaonPlus[i] = NULL;
      fEffKaonMinus[i] = NULL;
    }
}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskHptcorrId_tree::AliAnalysisTaskHptcorrId_tree(const char *name):
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
  fTreeVariableCentrality(0),
  fQ1_prot(0),
  fQ2_prot(0),
  fQ3_prot(0),
  fQ4_prot(0),
  fNch_prot(0),
  fQ1_pi(0),
  fQ2_pi(0),
  fQ3_pi(0),
  fQ4_pi(0),
  fNch_pi(0),
  fQ1_K(0),
  fQ2_K(0),
  fQ3_K(0),
  fQ4_K(0),
  fNch_K(0),
  fP1_prot(0),
  fP2_prot(0),
  fP3_prot(0),
  fP4_prot(0),
  fP1_pi(0),
  fP2_pi(0),
  fP3_pi(0),
  fP4_pi(0),
  fP1_K(0),
  fP2_K(0),
  fP3_K(0),
  fP4_K(0),
  fW1_prot(0),
  fW2_prot(0),
  fW3_prot(0),
  fW4_prot(0),
  fW1_pi(0),
  fW2_pi(0),
  fW3_pi(0),
  fW4_pi(0),
  fW1_K(0),
  fW2_K(0),
  fW3_K(0),
  fW4_K(0),
  fPtMax(0),
  fPtMin(0),
  fListTRKCorr(0), 
  fHistMCEffKaonPlus(0),
  fHistMCEffKaonMinus(0),
  fHistMCEffPionPlus(0),
  fHistMCEffPionMinus(0),
  fHistMCEffProtonPlus(0),
  fHistMCEffProtonMinus(0),
  fVertexZMax(0),
  fFBNo(0),
  fChi2TPC(0),
  fChi2ITS(0),
  fPIDnSigmaPionCut(0),
  fPIDnSigmaKaonCut(0),
  fPIDnSigmaProtonCut(0),
  fTPCcrossedrows(0),
  fEtaMax(0),
  fPileupCutVal(0),
  fCentralityEstimator_flag(0),
  fTreeName(0)
{
  for(int i=0; i<9; i++)
    {
      fEffProtonPlus[i] = NULL;
      fEffProtonMinus[i] = NULL;
      fEffPionPlus[i] = NULL;
      fEffPionMinus[i] = NULL;
      fEffKaonPlus[i] = NULL;
      fEffKaonMinus[i] = NULL;
    }
  
  fUtils = new AliAnalysisUtils();
  DefineInput (0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TTree::Class());
}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskHptcorrId_tree::~AliAnalysisTaskHptcorrId_tree()  {

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
void AliAnalysisTaskHptcorrId_tree::UserCreateOutputObjects()  {
    
    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner(kTRUE);
    fQAList     -> SetOwner(kTRUE);

    OpenFile(1);
    OpenFile(2);
    OpenFile(3);

    
    //QA Plots of Event Selection
    //fAODeventCuts.AddQAplotsToList(fQAList,kTRUE);
    fAODeventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE, fPileupCutVal);
    
    
    //Event Counter
    hNumberOfEvents = new TH1D ("hNumberOfEvents","",20,0,20);
    // hNumberOfEvents -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    
    
    //TTree object to store variables
    fTreeEvent = new TTree(fTreeName,"Event Tree");
    fTreeEvent->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
    //reconstructed
    fTreeEvent->Branch("fNch_prot", &fNch_prot, "fNch_prot/F");
    fTreeEvent->Branch("fQ1_prot", &fQ1_prot, "fQ1_prot/F");
    fTreeEvent->Branch("fQ2_prot", &fQ2_prot, "fQ2_prot/F");
    fTreeEvent->Branch("fQ3_prot", &fQ3_prot, "fQ3_prot/F");
    fTreeEvent->Branch("fQ4_prot", &fQ4_prot, "fQ4_prot/F");
    fTreeEvent->Branch("fNch_pi", &fNch_pi, "fNch_pi/F");
    fTreeEvent->Branch("fQ1_pi", &fQ1_pi, "fQ1_pi/F");
    fTreeEvent->Branch("fQ2_pi", &fQ2_pi, "fQ2_pi/F");
    fTreeEvent->Branch("fQ3_pi", &fQ3_pi, "fQ3_pi/F");
    fTreeEvent->Branch("fQ4_pi", &fQ4_pi, "fQ4_pi/F");
    fTreeEvent->Branch("fNch_K", &fNch_K, "fNch_K/F");
    fTreeEvent->Branch("fQ1_K", &fQ1_K, "fQ1_K/F");
    fTreeEvent->Branch("fQ2_K", &fQ2_K, "fQ2_K/F");
    fTreeEvent->Branch("fQ3_K", &fQ3_K, "fQ3_K/F");
    fTreeEvent->Branch("fQ4_K", &fQ4_K, "fQ4_K/F");

    //corrected
    fTreeEvent->Branch("fP1_prot", &fP1_prot, "fP1_prot/F");
    fTreeEvent->Branch("fP2_prot", &fP2_prot, "fP2_prot/F");
    fTreeEvent->Branch("fP3_prot", &fP3_prot, "fP3_prot/F");
    fTreeEvent->Branch("fP4_prot", &fP4_prot, "fP4_prot/F");
    fTreeEvent->Branch("fP1_pi", &fP1_pi, "fP1_pi/F");
    fTreeEvent->Branch("fP2_pi", &fP2_pi, "fP2_pi/F");
    fTreeEvent->Branch("fP3_pi", &fP3_pi, "fP3_pi/F");
    fTreeEvent->Branch("fP4_pi", &fP4_pi, "fP4_pi/F");
    fTreeEvent->Branch("fP1_K", &fP1_K, "fP1_K/F");
    fTreeEvent->Branch("fP2_K", &fP2_K, "fP2_K/F");
    fTreeEvent->Branch("fP3_K", &fP3_K, "fP3_K/F");
    fTreeEvent->Branch("fP4_K", &fP4_K, "fP4_K/F");
    fTreeEvent->Branch("fW1_prot", &fW1_prot, "fW1_prot/F");
    fTreeEvent->Branch("fW2_prot", &fW2_prot, "fW2_prot/F");
    fTreeEvent->Branch("fW3_prot", &fW3_prot, "fW3_prot/F");
    fTreeEvent->Branch("fW4_prot", &fW4_prot, "fW4_prot/F");
    fTreeEvent->Branch("fW1_pi", &fW1_pi, "fW1_pi/F");
    fTreeEvent->Branch("fW2_pi", &fW2_pi, "fW2_pi/F");
    fTreeEvent->Branch("fW3_pi", &fW3_pi, "fW3_pi/F");
    fTreeEvent->Branch("fW4_pi", &fW4_pi, "fW4_pi/F");
    fTreeEvent->Branch("fW1_K", &fW1_K, "fW1_K/F");
    fTreeEvent->Branch("fW2_K", &fW2_K, "fW2_K/F");
    fTreeEvent->Branch("fW3_K", &fW3_K, "fW3_K/F");
    fTreeEvent->Branch("fW4_K", &fW4_K, "fW4_K/F");
    
    
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeEvent);
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskHptcorrId_tree::UserExec(Option_t *)  {
  
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
	if (fCentralityEstimator_flag == 0)
	  lV0M = MultSelection->GetMultiplicityPercentile("V0M");
	else if (fCentralityEstimator_flag == 1)
	  lV0M = MultSelection->GetMultiplicityPercentile("CL0");
	else if (fCentralityEstimator_flag == 2)
	  lV0M = MultSelection->GetMultiplicityPercentile("CL1");
	else if (fCentralityEstimator_flag == 3)
	  lV0M = MultSelection->GetMultiplicityPercentile("CL2");

	//cout<<"V0M: "<<lV0M<<endl;
      }

    
    //"fTreeEvent" Tree Variable
    fTreeVariableCentrality = lV0M;

    //Reconstructed
    Float_t Q1[3]={0.0,0.0,0.0};
    Float_t Q2[3]={0.0,0.0,0.0};
    Float_t Q3[3]={0.0,0.0,0.0};
    Float_t Q4[3]={0.0,0.0,0.0};
    Float_t Nch[3]={0.0,0.0,0.0};

    //Corrected
    //+++++++++ numerator +++++++++++++++++
    Float_t P1[3]={0.0,0.0,0.0};
    Float_t P2[3]={0.0,0.0,0.0};
    Float_t P3[3]={0.0,0.0,0.0};
    Float_t P4[3]={0.0,0.0,0.0};
    //+++++++++ denominator +++++++++++++++++
    Float_t W1[3]={0.0,0.0,0.0};
    Float_t W2[3]={0.0,0.0,0.0};
    Float_t W3[3]={0.0,0.0,0.0};
    Float_t W4[3]={0.0,0.0,0.0};
    
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    Int_t ptBinNo = 0;
    Double_t BinCont = 0;
    Double_t EffWgt = 0;

    if(fListTRKCorr) GetMCEffCorrectionHist();

    Int_t centrality_bin = 0;
    if(lV0M > 0.0 && lV0M < 5.0)
      centrality_bin = 0;
    else if (lV0M > 5.0 && lV0M < 10.0)
      centrality_bin = 1;
    else if (lV0M > 10.0 && lV0M < 20.0)
      centrality_bin = 2;
    else if (lV0M > 20.0 && lV0M < 30.0)
      centrality_bin = 3;
    else if (lV0M > 30.0 && lV0M < 40.0)
      centrality_bin = 4;
    else if (lV0M > 40.0 && lV0M < 50.0)
      centrality_bin = 5;
    else if (lV0M > 50.0 && lV0M < 60.0)
      centrality_bin = 6;
    else if (lV0M > 60.0 && lV0M < 70.0)
      centrality_bin = 7;
    else if (lV0M > 70.0 && lV0M < 90.0)
      centrality_bin = 8;




    //Loop on reconstructed tracks
    
    for(Int_t itr=0; itr < fAODevent->GetNumberOfTracks(); itr++)
      {
	
	AliVTrack   *track = (AliVTrack*)fAODevent->GetTrack(itr);
	if(!track)      continue;
	AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(track);
	if(!aodtrack)      continue;

	Double_t trkPt = aodtrack->Pt();
	Double_t trkPhi = aodtrack->Phi();
	Double_t trkEta = aodtrack->Eta();
	Double_t trkCharge = aodtrack->Charge();
	Double_t trkDCAxy = aodtrack->DCA();
	Double_t trkDCAz = aodtrack->ZAtDCA();
	Double_t trkTPCNCls = aodtrack->GetTPCNcls();
	Double_t trkChi2PerNDF = aodtrack->Chi2perNDF();
	Double_t trkITSchi2 = aodtrack->GetITSchi2();
	Int_t trkITSNcls = aodtrack->GetITSNcls();
	Double_t trkITSchi2perNcls = trkITSchi2/trkITSNcls;
	Double_t trkTPCchi2perNcls = aodtrack->GetTPCchi2perCluster();
	Double_t trkTPCcrossedrows = aodtrack->GetTPCCrossedRows();


	//Track selectionL FilterBit 96
	//if(!aodtrack->TestFilterBit(96))  continue;
	if(!aodtrack->TestFilterBit(fFBNo))  continue;

	//cuts on TPCchi2perClstr and ITSchi2perClstr
	if (trkTPCcrossedrows < fTPCcrossedrows) continue;
	if (trkTPCchi2perNcls > fChi2TPC) continue;
	if (trkITSchi2perNcls > fChi2ITS) continue;


	//Kinematic cuts on pT and Eta
	if (TMath::Abs(trkEta) > fEtaMax) continue;
	if (trkPt < 0.2) continue;
	if (trkPt > 3.0) continue;
	if(TMath::Abs(trkCharge)!=1)continue;

	//PID selection
	Bool_t IsKaon = KaonSelector(track, fPIDnSigmaKaonCut);
	Bool_t IsPion = PionSelector(track, fPIDnSigmaPionCut);
	Bool_t IsProton = ProtonSelector(track, fPIDnSigmaProtonCut);
	if (!IsKaon && !IsPion && !IsProton) continue;

	//Check if a particle is selected as more than one type
	Int_t flag = 0;
	if(IsKaon) flag+=1;
	if(IsPion) flag+=1;
	if(IsProton) flag+=1;
	if(flag>1)
	  {
	    cout<<"Particle identified as more than one PID: flag= "<<flag<<endl;
	    continue;
	  }

	//cout<<fPtMin<<"\t++++++++++++++++++++++\t"<<fPtMax<<endl;

	if (IsProton && trkPt > fPtMin && trkPt < fPtMax)
	  {
	    if (fEffProtonPlus[centrality_bin])
	      {
		ptBinNo = fEffProtonPlus[centrality_bin]->FindBin(trkPt);
		BinCont = fEffProtonPlus[centrality_bin]->GetBinContent(ptBinNo);
		if(BinCont!=0) EffWgt = 1.0/BinCont;
		
		Q1[0]=Q1[0]+TMath::Power(trkPt, 1.0);
		Q2[0]=Q2[0]+TMath::Power(trkPt, 2.0);
		Q3[0]=Q3[0]+TMath::Power(trkPt, 3.0);
		Q4[0]=Q4[0]+TMath::Power(trkPt, 4.0);
		Nch[0]=Nch[0]+1;

		P1[0]=P1[0]+TMath::Power(trkPt*EffWgt, 1.0);
		P2[0]=P2[0]+TMath::Power(trkPt*EffWgt, 2.0);
		P3[0]=P3[0]+TMath::Power(trkPt*EffWgt, 3.0);
		P4[0]=P4[0]+TMath::Power(trkPt*EffWgt, 4.0);

		W1[0]=W1[0]+TMath::Power(EffWgt, 1.0);
		W2[0]=W2[0]+TMath::Power(EffWgt, 2.0);
		W3[0]=W3[0]+TMath::Power(EffWgt, 3.0);
		W4[0]=W4[0]+TMath::Power(EffWgt, 4.0);
	      }
	  }
	if (IsPion && trkPt > 0.2 && trkPt < 3.0)
	  {
	    if (fEffPionPlus[centrality_bin])
	      {
		ptBinNo = fEffPionPlus[centrality_bin]->FindBin(trkPt);
		BinCont = fEffPionPlus[centrality_bin]->GetBinContent(ptBinNo);
		if(BinCont!=0) EffWgt = 1.0/BinCont;
		
		Q1[1]=Q1[1]+TMath::Power(trkPt, 1.0);
		Q2[1]=Q2[1]+TMath::Power(trkPt, 2.0);
		Q3[1]=Q3[1]+TMath::Power(trkPt, 3.0);
		Q4[1]=Q4[1]+TMath::Power(trkPt, 4.0);
		Nch[1]=Nch[1]+1;

		P1[1]=P1[1]+TMath::Power(trkPt*EffWgt, 1.0);
		P2[1]=P2[1]+TMath::Power(trkPt*EffWgt, 2.0);
		P3[1]=P3[1]+TMath::Power(trkPt*EffWgt, 3.0);
		P4[1]=P4[1]+TMath::Power(trkPt*EffWgt, 4.0);

		W1[1]=W1[1]+TMath::Power(EffWgt, 1.0);
		W2[1]=W2[1]+TMath::Power(EffWgt, 2.0);
		W3[1]=W3[1]+TMath::Power(EffWgt, 3.0);
		W4[1]=W4[1]+TMath::Power(EffWgt, 4.0);
	      }
	  }
	if (IsKaon && trkPt > 0.2 && trkPt < 3.0)
	  {
	    if (fEffKaonPlus[centrality_bin])
	      {
		ptBinNo = fEffKaonPlus[centrality_bin]->FindBin(trkPt);
		BinCont = fEffKaonPlus[centrality_bin]->GetBinContent(ptBinNo);
		if(BinCont!=0) EffWgt = 1.0/BinCont; 
		    
		Q1[2]=Q1[2]+TMath::Power(trkPt, 1.0);
		Q2[2]=Q2[2]+TMath::Power(trkPt, 2.0);
		Q3[2]=Q3[2]+TMath::Power(trkPt, 3.0);
		Q4[2]=Q4[2]+TMath::Power(trkPt, 4.0);
		Nch[2]=Nch[2]+1;

		P1[2]=P1[2]+TMath::Power(trkPt*EffWgt, 1.0);
		P2[2]=P2[2]+TMath::Power(trkPt*EffWgt, 2.0);
		P3[2]=P3[2]+TMath::Power(trkPt*EffWgt, 3.0);
		P4[2]=P4[2]+TMath::Power(trkPt*EffWgt, 4.0);

		W1[2]=W1[2]+TMath::Power(EffWgt, 1.0);
		W2[2]=W2[2]+TMath::Power(EffWgt, 2.0);
		W3[2]=W3[2]+TMath::Power(EffWgt, 3.0);
		W4[2]=W4[2]+TMath::Power(EffWgt, 4.0);
	      }
	  }
      }      
    //end reconstructed track loop
    
    fNch_prot=Nch[0];
    fQ1_prot=Q1[0];
    fQ2_prot=Q2[0];
    fQ3_prot=Q3[0];
    fQ4_prot=Q4[0];
    fNch_pi=Nch[1];
    fQ1_pi=Q1[1];
    fQ2_pi=Q2[1];
    fQ3_pi=Q3[1];
    fQ4_pi=Q4[1];
    fNch_K=Nch[2];
    fQ1_K=Q1[2];
    fQ2_K=Q2[2];
    fQ3_K=Q3[2];
    fQ4_K=Q4[2];

    fP1_prot=P1[0];
    fP2_prot=P2[0];
    fP3_prot=P3[0];
    fP4_prot=P4[0];
    fP1_pi=P1[1];
    fP2_pi=P2[1];
    fP3_pi=P3[1];
    fP4_pi=P4[1];
    fP1_K=P1[2];
    fP2_K=P2[2];
    fP3_K=P3[2];
    fP4_K=P4[2];

    fW1_prot=W1[0];
    fW2_prot=W2[0];
    fW3_prot=W3[0];
    fW4_prot=W4[0];
    fW1_pi=W1[1];
    fW2_pi=W2[1];
    fW3_pi=W3[1];
    fW4_pi=W4[1];
    fW1_K=W1[2];
    fW2_K=W2[2];
    fW3_K=W3[2];
    fW4_K=W4[2];
    
    
    fTreeEvent->Fill();

    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeEvent);
}    
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHptcorrId_tree::GetEvent ()  //event cuts copied from my code written earlier 

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
    if ( vertex_tracks->GetZ() < -fVertexZMax || vertex_tracks->GetZ() > +fVertexZMax)
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
    if (mult_percentile < 0.0 || mult_percentile > 90.0)
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
Bool_t AliAnalysisTaskHptcorrId_tree::PassedTrackQualityCuts (AliAODTrack *track)  {
    
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
Bool_t AliAnalysisTaskHptcorrId_tree::KaonSelector(AliVTrack *track, Double_t nSigmaCut)  {
 
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
      if(TMath::Abs(fTPCnSigmaKaon) < nSigmaCut)
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
      if (fTPCplusTOFnSigmaKaon < nSigmaCut)
	return kTRUE;
      else
	return kFALSE;
    }
  
  
  
}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHptcorrId_tree::PionSelector(AliVTrack *track, Double_t nSigmaCut)  {
  
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
      
      if(TMath::Abs(fTPCnSigmaPion) < nSigmaCut)
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
      
      if (fTPCplusTOFnSigmaPion < nSigmaCut)
	return kTRUE;
      else
	return kFALSE;
    }
  
}

//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHptcorrId_tree::ProtonSelector(AliVTrack *track, Double_t nSigmaCut)  {
  
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

      if(TMath::Abs(fTPCnSigmaProton) < nSigmaCut)
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
      
      if (fTPCplusTOFnSigmaProton < nSigmaCut)
	return kTRUE;
      else
	return kFALSE;
    }
  
}


//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHptcorrId_tree::PassedPIDSelection (AliAODTrack *track, AliPID::EParticleType type)  {
    
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
Bool_t AliAnalysisTaskHptcorrId_tree::PassedSingleParticlePileUpCuts(AliAODTrack *track)
{
  Bool_t passedTrackPileupCut = (kTRUE);
  if (!(track->HasPointOnITSLayer(1)) && !(track->HasPointOnITSLayer(4)) && !(track->HasPointOnITSLayer(5)) && !(track->GetTOFBunchCrossing() == 0))
    passedTrackPileupCut = (kFALSE);
  return passedTrackPileupCut;
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskHptcorrId_tree::GetMCEffCorrectionHist()
{
  if(fListTRKCorr)
    {
      fHistMCEffKaonPlus = (TH1D*) fListTRKCorr->FindObject("histKaonPlusEff");
      fHistMCEffKaonMinus = (TH1D*) fListTRKCorr->FindObject("histKaonMinusEff");

      fHistMCEffPionPlus = (TH1D*) fListTRKCorr->FindObject("histPionPlusEff");
      fHistMCEffPionMinus = (TH1D*) fListTRKCorr->FindObject("histPionMinusEff");

      fHistMCEffProtonPlus = (TH1D*) fListTRKCorr->FindObject("histProtonPlusEff");
      fHistMCEffProtonMinus = (TH1D*) fListTRKCorr->FindObject("histProtonMinusEff");

      for(int i=0; i<9; i++)
	{
	  fEffProtonPlus[i] = (TH1D*)fListTRKCorr->FindObject(Form("EffProtonPlus%d",i));
	  fEffProtonMinus[i] = (TH1D*)fListTRKCorr->FindObject(Form("EffProtonMinus%d",i));
	  fEffPionPlus[i] = (TH1D*)fListTRKCorr->FindObject(Form("EffPionPlus%d",i));
	  fEffPionMinus[i] = (TH1D*)fListTRKCorr->FindObject(Form("EffPionMinus%d",i));
	  fEffKaonPlus[i] = (TH1D*)fListTRKCorr->FindObject(Form("EffKaonPlus%d",i));
	  fEffKaonMinus[i] = (TH1D*)fListTRKCorr->FindObject(Form("EffKaonMinus%d",i));

	}
    }

}
 //_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskHptcorrId_tree::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
