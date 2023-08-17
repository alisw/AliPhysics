
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
#include "AliAnalysisTaskDiffPtFluc_PIDhadrons.h"

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

class AliAnalysisTaskDiffPtFluc_PIDhadrons;
ClassImp(AliAnalysisTaskDiffPtFluc_PIDhadrons)

//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskDiffPtFluc_PIDhadrons::AliAnalysisTaskDiffPtFluc_PIDhadrons():
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
  hNumberOfKaonEtaLess0(0),
  hNumberOfPionEtaLess0(0),
  hNumberOfProtonEtaLess0(0),
  fTreeVariableCentrality(0),
  fMeanPt_less0(0),
  fMeanPt_greaterEtaMin(0),
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
  fPIDnSigmaCut(0),
  fTPCcrossedrows(0),
  fEtaMin(0),
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
   for(int i=0; i<14; i++)
    {
      fPt_factor[i] = 0;
      fPt_factor_pion[i] = 0;
      fPt_factor_kaon[i] = 0;
      fPt_factor_proton[i] = 0;
    }
}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskDiffPtFluc_PIDhadrons::AliAnalysisTaskDiffPtFluc_PIDhadrons(const char *name):
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
  hNumberOfKaonEtaLess0(0),
  hNumberOfPionEtaLess0(0),
  hNumberOfProtonEtaLess0(0),
  fTreeVariableCentrality(0),
  fMeanPt_less0(0),
  fMeanPt_greaterEtaMin(0),
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
  fPIDnSigmaCut(0),
  fTPCcrossedrows(0),
  fEtaMin(0),
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
   for(int i=0; i<14; i++)
    {
      fPt_factor[i] = 0;
      fPt_factor_pion[i] = 0;
      fPt_factor_kaon[i] = 0;
      fPt_factor_proton[i] = 0;
    }
  
  fUtils = new AliAnalysisUtils();
  DefineInput (0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TTree::Class());
}
//_____________________________________________________________________________________________________________________________________
AliAnalysisTaskDiffPtFluc_PIDhadrons::~AliAnalysisTaskDiffPtFluc_PIDhadrons()  {

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
void AliAnalysisTaskDiffPtFluc_PIDhadrons::UserCreateOutputObjects()  {
    
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
    fAODeventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);
    
    
    //Event Counter
    hNumberOfEvents = new TH1D ("hNumberOfEvents","",20,0,20);
    // hNumberOfEvents -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    

    //Number of Kaon finally getting selected with specified cuts event wise
    hNumberOfKaonEtaLess0     = new TH1D ("hNumberOfKaonEtaLess0","",3000,0,3000);
    fOutputList -> Add(hNumberOfKaonEtaLess0);
    

    //Number of Pion finally getting selected with specified cuts event wise
    hNumberOfPionEtaLess0     = new TH1D ("hNumberOfPionEtaLess0","",3000,0,3000);
    fOutputList -> Add(hNumberOfPionEtaLess0);
    
    //Number of Proton finally getting selected with specified cuts event wise
    hNumberOfProtonEtaLess0     = new TH1D ("hNumberOfProtonEtaLess0","",3000,0,3000);
    fOutputList -> Add(hNumberOfProtonEtaLess0);
    

    
    //TTree object to store variables
    fTreeEvent = new TTree(fTreeName,"Event Tree");
    fTreeEvent->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
    fTreeEvent->Branch("fMeanPt_less0",&fMeanPt_less0,"fMeanPt_less0/F");
    fTreeEvent->Branch("fMeanPt_greaterEtaMin",&fMeanPt_greaterEtaMin,"fMeanPt_greaterEtaMin/F");
    fTreeEvent->Branch("fPt_factor",&fPt_factor,"fPt_factor[14]/F");
    fTreeEvent->Branch("fPt_factor_pion",&fPt_factor_pion,"fPt_factor_pion[14]/F");
    fTreeEvent->Branch("fPt_factor_kaon",&fPt_factor_kaon,"fPt_factor_kaon[14]/F");
    fTreeEvent->Branch("fPt_factor_proton",&fPt_factor_proton,"fPt_factor_proton[14]/F");

    
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeEvent);
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskDiffPtFluc_PIDhadrons::UserExec(Option_t *)  {
  
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

   
    //Initialize number 0f K+, K-, pi+, pi-, p and p-bar per event
    
 
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    Int_t ptBinNo = 0;
    Double_t BinCont = 0;
    Double_t EffWgt = 0;

    if(fListTRKCorr)
      {
	cout<<"## Got efficiency histograms..## "<<endl;
	GetMCEffCorrectionHist();
      }
    else
      cout<<"################ No histograms found ############### "<<endl;



    TH1D *fPt_profile = new TH1D("fPt_profile","fPt_profile", 14, 0.2, 3.0);
    Double_t pT_sum_etaLess0 = 0.0;
    Double_t N_sum_etaLess0 = 0.0;
    Double_t pT_sum_etaGreaterEtamin = 0.0;
    Double_t N_sum_etaGreaterEtamin = 0.0;


    TH1D *fPt_profile_pion = new TH1D("fPt_profile_pion","fPt_profile_pion", 14, 0.2, 3.0);
    TH1D *fPt_profile_kaon = new TH1D("fPt_profile_kaon","fPt_profile_kaon", 14, 0.2, 3.0);
    TH1D *fPt_profile_proton = new TH1D("fPt_profile_proton","fPt_profile_proton", 14, 0.2, 3.0);
    Double_t N_sumPion_etaLess0 = 0.0;
    Double_t N_sumKaon_etaLess0 = 0.0;
    Double_t N_sumProton_etaLess0 = 0.0;

    //Loop on reconstructed tracks
    
    for(Int_t itr=0; itr < fAODevent->GetNumberOfTracks(); itr++)
      {
	
	AliVTrack   *track = (AliVTrack*)fAODevent->GetTrack(itr);
	if(!track)      continue;
	AliAODTrack *aodtrack  = dynamic_cast<AliAODTrack*>(track);
	if(!aodtrack)      continue;


	//Track selectionL FilterBit 96
	//if(!aodtrack->TestFilterBit(96))  continue;
	if(!aodtrack->TestFilterBit(fFBNo))  continue;


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

	
	//cuts on TPCchi2perClstr and ITSchi2perClstr and TPC nCrossedRows
	if (trkTPCcrossedrows < fTPCcrossedrows)
	  continue;
	if (trkTPCchi2perNcls > fChi2TPC)
	  continue;
	if (trkITSchi2perNcls > fChi2ITS)
	  continue;


	//Kinematic cuts on pT and Eta
	if (TMath::Abs(trkEta) > 0.8) continue;
	if (trkPt < 0.2) continue;
	if (trkPt > 3.0) continue;

	
	//PID selection
	Bool_t IsKaon = KaonSelector(track, fPIDnSigmaCut);
	Bool_t IsPion = PionSelector(track, fPIDnSigmaCut);
	Bool_t IsProton = ProtonSelector(track, fPIDnSigmaCut);
	if (!IsKaon && !IsPion && !IsProton) continue;

	//Check if a particle is selected as more than one type
	Int_t flag = 0;
	if(IsKaon) flag+=1;
	if(IsPion) flag+=1;
	if(IsProton) flag+=1;
	if(flag>1)
	  {
	    cout<<"Particle identified as more than on PID: flag= "<<flag<<endl;
	    continue;
	  }
	

	if(TMath::Abs(trkCharge) > 0)
	  {
	    if(trkEta < 0.0)
	      {
		fPt_profile->Fill(trkPt);
		pT_sum_etaLess0 += trkPt;
		N_sum_etaLess0 += 1.0;

		if(IsPion)
		  {
		    fPt_profile_pion->Fill(trkPt);
		    N_sumPion_etaLess0 += 1.0;
		  }
		if(IsKaon)
		  {
		    fPt_profile_kaon->Fill(trkPt);
		    N_sumKaon_etaLess0 += 1.0;
		  }
		if(IsProton && trkPt > 0.4)
		  {
		    fPt_profile_proton->Fill(trkPt);
		    N_sumProton_etaLess0 += 1.0;
		  }

	      }
	    
	    if(trkEta > fEtaMin)
	      {
		pT_sum_etaGreaterEtamin += trkPt;
		N_sum_etaGreaterEtamin += 1.0;
	      }
	  }

	
	
      }      
    //end reconstructed track loop
    
    
    
    //Kaon
    hNumberOfKaonEtaLess0->Fill(N_sumKaon_etaLess0);
    //Pion
    hNumberOfPionEtaLess0->Fill(N_sumPion_etaLess0);
    //Proton
    hNumberOfProtonEtaLess0->Fill(N_sumProton_etaLess0);

    //Tree Variables++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    fTreeVariableCentrality=lV0M;
    fMeanPt_less0=pT_sum_etaLess0/N_sum_etaLess0;
    fMeanPt_greaterEtaMin=pT_sum_etaGreaterEtamin/N_sum_etaGreaterEtamin;
    for(int i=0; i<14; i++)
      {
	fPt_factor[i]=fPt_profile->GetBinContent(i+1)/N_sum_etaLess0;

	fPt_factor_pion[i]=fPt_profile_pion->GetBinContent(i+1)/N_sumPion_etaLess0;
	fPt_factor_kaon[i]=fPt_profile_kaon->GetBinContent(i+1)/N_sumKaon_etaLess0;
	fPt_factor_proton[i]=fPt_profile_proton->GetBinContent(i+1)/N_sumProton_etaLess0;
      }
    if (N_sum_etaGreaterEtamin>0 && N_sum_etaLess0>0 && N_sumPion_etaLess0>0 && N_sumKaon_etaLess0>0 && N_sumProton_etaLess0>0)
      fTreeEvent->Fill();

    fPt_profile->Delete();
    fPt_profile_pion->Delete();
    fPt_profile_kaon->Delete();
    fPt_profile_proton->Delete();
    
    

 
    PostData(1, fOutputList);
    PostData(2, fQAList);
    PostData(3, fTreeEvent);
}    
//_____________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons::GetEvent ()  //event cuts copied from my code written earlier 

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
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons::PassedTrackQualityCuts (AliAODTrack *track)  {
    
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
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons::KaonSelector(AliVTrack *track, Double_t nSigmaCut)  {
 
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
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons::PionSelector(AliVTrack *track, Double_t nSigmaCut)  {
  
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
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons::ProtonSelector(AliVTrack *track, Double_t nSigmaCut)  {
  
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
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons::PassedPIDSelection (AliAODTrack *track, AliPID::EParticleType type)  {
    
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
Bool_t AliAnalysisTaskDiffPtFluc_PIDhadrons::PassedSingleParticlePileUpCuts(AliAODTrack *track)
{
  Bool_t passedTrackPileupCut = (kTRUE);
  if (!(track->HasPointOnITSLayer(1)) && !(track->HasPointOnITSLayer(4)) && !(track->HasPointOnITSLayer(5)) && !(track->GetTOFBunchCrossing() == 0))
    passedTrackPileupCut = (kFALSE);
  return passedTrackPileupCut;
}
//_____________________________________________________________________________________________________________________________________
void AliAnalysisTaskDiffPtFluc_PIDhadrons::GetMCEffCorrectionHist()
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
void AliAnalysisTaskDiffPtFluc_PIDhadrons::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
