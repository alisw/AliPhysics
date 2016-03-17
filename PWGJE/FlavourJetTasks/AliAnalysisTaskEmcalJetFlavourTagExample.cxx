/**************************************************************************
* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/*
*   This sample task demostrates the basics of tagging a jet. The PID response
*   is retrived for both the TPC.  A Hadronic tag is implemented for
*   clusters matched to tracks in a jet with a E/P < 0.2
*
*   Author: Andrew Castro (UTK) and Joel Mazer (UTK)
*/

#include "AliAnalysisTaskEmcalJetFlavourTagExample.h"

// general ROOT includes                                                                                                                                                  
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TParameter.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVector3.h>
#include <TObjArray.h>

// AliROOT includes                                                                                                                         
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliEmcalJet.h"
#include "AliAODJet.h"
#include "AliVCluster.h"
#include "AliVTrack.h"
#include <AliVEvent.h>
#include <AliVParticle.h>
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliEmcalParticle.h"
#include "AliESDCaloCluster.h"
#include <AliESDtrackCuts.h>
#include "AliPID.h"

// event handler (and pico's) includes                                                                                                      
#include <AliInputEventHandler.h>
#include <AliVEventHandler.h>
#include "AliESDInputHandler.h"
#include "AliPicoTrack.h"
#include "AliEventPoolManager.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"

// PID includes                                                                                                                             
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliESDpid.h"

// magnetic field includes
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalJetFlavourTagExample)

//________________________________________________________________________
AliAnalysisTaskEmcalJetFlavourTagExample::AliAnalysisTaskEmcalJetFlavourTagExample() : 
  AliAnalysisTaskEmcalJet("heavyF",kFALSE), 
  event(0),
  fCuts(0),
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0),
  fJetHIpt(50.0),
  fTrackPtCut(2.0),
  fTrackEta(0.9),
  fesdTrackCuts(0),
  fPIDResponse(0x0), fTPCResponse(),
  fJetsCont(0), fTracksCont(0), fCaloClustersCont(0), fTracksJetCont(0), fCaloClustersJetCont(0),
  fESD(0), fAOD(0),
  fHistJetPhi(0),
  fHistCorJetPt(0), fHistJetPt(0),
  fHistHighJetPt(0),
  fHistnSigElecPt(0),
  fHistdEdXvsPt(0),
  fHistnJetTrackvnJetClusters(0)
{
  // Default constructor.


  SetMakeGeneralHistograms(kTRUE);

}

//________________________________________________________________________
AliAnalysisTaskEmcalJetFlavourTagExample::AliAnalysisTaskEmcalJetFlavourTagExample(const char *name) :
  AliAnalysisTaskEmcalJet(name,kTRUE),
  event(0),
  fCuts(0),
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0),
  fJetHIpt(50.0),
  fTrackPtCut(2.0),
  fTrackEta(0.9),
  fesdTrackCuts(0),
  fPIDResponse(0x0), fTPCResponse(),
  fJetsCont(0), fTracksCont(0), fCaloClustersCont(0), fTracksJetCont(0), fCaloClustersJetCont(0),
  fESD(0), fAOD(0),
  fHistJetPhi(0),
  fHistCorJetPt(0), fHistJetPt(0),
  fHistHighJetPt(0),
  fHistnSigElecPt(0),
  fHistdEdXvsPt(0),
  fHistnJetTrackvnJetClusters(0)
{ 

   SetMakeGeneralHistograms(kTRUE);
 
   DefineInput(0,TChain::Class());
   DefineOutput(1, TList::Class());
}

//_______________________________________________________________________
AliAnalysisTaskEmcalJetFlavourTagExample::~AliAnalysisTaskEmcalJetFlavourTagExample()
{
  // destructor
  //
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetFlavourTagExample::UserCreateOutputObjects()
{
  if (! fCreateHisto)
    return;
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  //fJetsCont           = GetJetContainer(0);
  if(fJetsCont) { //get particles and clusters connected to jets
    fTracksJetCont       = fJetsCont->GetParticleContainer();
    fCaloClustersJetCont = fJetsCont->GetClusterContainer();
  }
 else {        //no jets, just analysis tracks and clusters
  fTracksCont       = GetParticleContainer(0);
  fCaloClustersCont = GetClusterContainer(0);
}
fTracksCont->SetClassName("AliVTrack");
fCaloClustersCont->SetClassName("AliVCluster");

  fHistJetPhi                = new TH1F("NjetvsPhi", "NjetvsPhi", 288,-2*TMath::Pi(),2*TMath::Pi());
  fHistJetPt                 = new TH1F("NjetvsJetPt", "NjetvsJetPt", 300, 0, 300);
  fOutput->Add(fHistJetPhi);
  fOutput->Add(fHistJetPt);

 
  TString histname;

  fHistCorJetPt		            = new TH1F("CorrJetPt", "CorrJetPt", 300, -100, 200);
  fHistnSigElecPt             = new TH2F("nsigvsPt(TPC)","nSigmaElectronTPC_v_Pt",60,0,60,40,-10,10);
  fHistdEdXvsPt               = new TH2F("dEdXvsTrackPt","dEdXvsTrackPt",60,0,60,80,0,80);
  fHistnJetTrackvnJetClusters = new TH2F("NumbJetTracksvJetClusters","NumbJetTracksvJetClusters",21,0,20,21,0,20);
  fHistHighJetPt              = new TH1F("HighestPtJetPerEvent","HighJetPt",80,0,80);
    
  TString name;
  TString title;

  fOutput->Add(fHistCorJetPt);
  fOutput->Add(fHistnSigElecPt);
  fOutput->Add(fHistdEdXvsPt);
  fOutput->Add(fHistnJetTrackvnJetClusters);
  fOutput->Add(fHistHighJetPt);


  // ****************************** PID *****************************************************                                               
  // set up PID handler                                                                                                                     
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(!inputHandler) {
    AliFatal("Input handler needed");
  }

  // PID response object                                                                                                                    
  //fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();                                                                         
  //  inputHandler->CreatePIDResponse(fIsMC);         // needed to create object, why though?                                                 
  fPIDResponse = inputHandler->GetPIDResponse();
  if (!fPIDResponse) {
    AliError("PIDResponse object was not created");
  }
  // ***************************************************************************************

  PostData(1, fOutput);

}

//________________________________________________________
void AliAnalysisTaskEmcalJetFlavourTagExample::ExecOnce()
{
  //  Initialize the analysis
  AliAnalysisTaskEmcalJet::ExecOnce();
  
  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;


} // end of ExecOnce

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetFlavourTagExample::Run()
{
  // check to see if we have any tracks
  if (!fTracks)  return kTRUE;
  if (!fJets)  return kTRUE;

  // what kind of event do we have: AOD or ESD?
  Bool_t useAOD;
  if (dynamic_cast<AliAODEvent*>(InputEvent())) useAOD = kTRUE;
  else useAOD = kFALSE;
  
  // if we have ESD event, set up ESD object
  if(!useAOD){
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) {
      AliError(Form("ERROR: fESD not available\n"));
      return kTRUE;
    }
  }

  // if we have AOD event, set up AOD object
  if(useAOD){
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) {
      AliError(Form("ERROR: fAOD not available\n"));
      return kTRUE;
    }
  }
  
  // get magnetic field info for DCA
  Double_t  MagF = fESD->GetMagneticField();
  Double_t MagSign = 1.0;
  if(MagF<0)MagSign = -1.0;
  // set magnetic field
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    AliMagF* field = new AliMagF("Maps","Maps", MagSign, MagSign, AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }

  // get vertex information
  Double_t fvertex[3]={0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(fvertex);
  //Double_t zVtx=fvertex[2];

  // create pointer to list of input event                                                                                                  
  TList *list = InputEvent()->GetList();
  if(!list) {
    AliError(Form("ERROR: list not attached\n"));
    return kTRUE;
  }

  // background density                                                                                                                                                                                                                               
  fRhoVal = fRho->GetVal();

  // initialize TClonesArray pointers to jets and tracks                                                                                    
  TClonesArray *jets = 0;
  //TClonesArray *tracks = 0;
  //TClonesArray *clusters = 0;
  
  // get Jets object                                                                                                                        
  jets = dynamic_cast<TClonesArray*>(list->FindObject(fJets));
  if(!jets){
    AliError(Form("Pointer to jets %s == 0", fJets->GetName()));
    return kTRUE;
  } // verify existence of jets
  
  // get number of jets and tracks                                                                                                          
  const Int_t Njets = jets->GetEntries();
  if(Njets<1)     return kTRUE;
  
  if (fTracksCont) {
    fTracksCont->ResetCurrentID();
    AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    while(track) {
      track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    }
  }
  if (fCaloClustersCont) {
    fCaloClustersCont->ResetCurrentID();
    AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster();
    while(cluster) {
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, fVertex);
      cluster = fCaloClustersCont->GetNextAcceptCluster();
    }
  }
  
    //  Start Jet Analysis
    // initialize jet parameters
    Int_t ijethi=-1;
    Double_t highestjetpt=0.0;
  
  // loop over jets in an event - to find highest jet pT and apply some cuts && JetQA Sparse
  for (Int_t ijet = 0; ijet < Njets; ijet++){
    // get our jets
    AliEmcalJet *jet = static_cast<AliEmcalJet*>(jets->At(ijet));
    if (!jet) continue;
    
    // apply jet cuts
    if(!AcceptMyJet(jet)) continue;
    
    if(highestjetpt<jet->Pt()){
      ijethi=ijet;
      highestjetpt=jet->Pt();
    }
  } // end of looping over jets
  
    fHistHighJetPt->Fill(ijethi);
  
 // **********************************************************************
 //                JET LOOP
 // **********************************************************************
  
    // loop over jets in the event and make appropriate cuts
    for (Int_t iJets = 0; iJets < Njets; ++iJets) {
      AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(iJets));
      if (!jet)  // see if we have a jet
        continue;
      
      // apply jet cuts
      if(!AcceptMyJet(jet)) continue;
      
      //AliEmcalJet::EFlavourTag tag=AliEmcalJet::kDStar;
      //jet->AddFlavourTag(tag);
    
 
      Int_t JetClusters = jet->GetNumberOfClusters();
      Int_t JetTracks = jet -> GetNumberOfTracks();
      fHistnJetTrackvnJetClusters->Fill(JetClusters,JetTracks);
      // Initializations and Calculations
      Double_t jetPt = -500;                                     // initialize corr jet pt LOCAL
      Double_t jetarea = -500;                                   // initialize jet area
      jetPt = jet->Pt() - jetarea*fRhoVal;                       // semi-corrected pT of jet from GLOBAL rho value
      fHistCorJetPt->Fill(jetPt);
      
      Bool_t bkgrnd1 = kFALSE;
      
      if(jet->Pt() > fJetHIpt) {
        if(!fTracksCont || !fCaloClustersCont) continue;
          
        //******************************Cluster Matched To Closest Track
        //**************************************************************
        Int_t NumbTrackContainer = -999;
        NumbTrackContainer = fTracksCont->GetNParticles();
        for(int iTracks = 0; iTracks <= NumbTrackContainer; iTracks++){
          AliVTrack *AcceptedTrack =static_cast<AliVTrack*>(fTracksCont->GetParticle(iTracks));
          if(!AcceptedTrack){
            AliError(Form("Couldn't get AliVTrack Container %d\n", iTracks));
            continue;
          }
          if(!IsJetTrack(jet,iTracks,kFALSE))continue;
          //Get matched cluster
          Int_t emc1 = AcceptedTrack->GetEMCALcluster();
        
          Double_t acceptTrackP = AcceptedTrack->P();
          Double_t acceptTrackPt = AcceptedTrack->Pt();
          Double_t nSigmaElectron_TPC = fPIDResponse->NumberOfSigmasTPC(AcceptedTrack,AliPID::kElectron);
          fHistnSigElecPt->Fill(acceptTrackPt,nSigmaElectron_TPC);

          AliESDtrack *ESDacceptedTrack = static_cast<AliESDtrack*>(AcceptedTrack);
          if(!ESDacceptedTrack){
            AliError(Form("Couldn't get AliESDTrack %d\n", iTracks));
            continue;
          }
          Double_t dEdx = AcceptedTrack->GetTPCsignal();
                              fHistdEdXvsPt->Fill(acceptTrackPt,dEdx);
          if(fCaloClustersCont && emc1>=0) {
            AliVCluster *clusMatch = fCaloClustersCont->GetCluster(emc1);
            if(!clusMatch){
              AliError(Form("Couldn't get matched AliVCluster %d\n", emc1));
              continue;
            }
 
            Double_t mClusterE = clusMatch->E();
            Double_t EovP_mc = -999;
            EovP_mc = mClusterE/acceptTrackP;
          
            if(EovP_mc < 0.2) bkgrnd1 = kTRUE;
          }
          
      } //loop over tracks for matching to closest cluster

     

        
    } // highest pt jet cut
    
      
      if(bkgrnd1 == kTRUE) {
        AliEmcalJet::EFlavourTag tag=AliEmcalJet::kBckgrd1;
        jet->AddFlavourTag(tag);
    }
        
  } // LOOP over JETS in event


  return kTRUE;
  
}
//________________________________________________________________________
void AliAnalysisTaskEmcalJetFlavourTagExample::Terminate(Option_t *)
{
  cout<<"###########################"<<endl;
  cout<<"####   Task Finished   ####"<<endl;
  cout<<"###########################"<<endl;
  cout<<"###########################"<<endl;
} // end of terminate

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetFlavourTagExample::AcceptMyJet(AliEmcalJet *jet) {
  //applies all jet cuts except pt
  if ((jet->Phi()<fPhimin)||(jet->Phi()>fPhimax)) return 0;
  if ((jet->Eta()<fEtamin)||(jet->Eta()>fEtamax)) return 0;
  if (jet->Area()<fAreacut) return 0;
  // prevents 0 area jets from sneaking by when area cut == 0
  if (jet->Area()==0) return 0;
  //exclude jets with extremely high pt tracks which are likely misreconstructed
  if(jet->MaxTrackPt()>100) return 0;
  
  //passed all above cuts
  return 1;
}





