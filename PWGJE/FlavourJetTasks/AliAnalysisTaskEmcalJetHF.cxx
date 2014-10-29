//
// Author: A Castro (UTK)
// Last Modified: August 11, 2014

#include "AliAnalysisTaskEmcalJetHF.h"

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
#include "AliTPCdEdxInfo.h"
//#include "AliCaloTrackESDReader.h"
//#include "AliCaloTrackAODReader.h"
//#include "AliCaloTrackReader.h"

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

#include <AliInputEventHandler.h>
#include <AliVEventHandler.h>

// magnetic field includes
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEmcalJetHF)

//________________________________________________________________________
AliAnalysisTaskEmcalJetHF::AliAnalysisTaskEmcalJetHF() : 
  AliAnalysisTaskEmcalJet("heavyF",kFALSE), 
  event(0),
  fillHist(0),
  fEventTrigEMCALL1Gamma1(0),
  fEventTrigEMCALL1Gamma2(0),
  fGlobalQA(0),
  fInputEvent(0x0),
  fCuts(0),
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0),
  fJetHIpt(20.0),
  fTrackPtCut(2.0),
  fTrackEta(0.9),
  fTrkQAcut(0),
  fesdTrackCuts(0),
  fPIDResponse(0x0), fTPCResponse(),
  fEsdtrackCutsITSTPC(),
  fEsdtrackCutsTPC(),
  fEsdtrackCutsITS(),
  fJetsCont(0), fTracksCont(0), fCaloClustersCont(0), fTracksJetCont(0), fCaloClustersJetCont(0),
  fESD(0), fAOD(0),
  fHistRhovsCent(0),
  fHistJetPhi(0),
  fHistCorJetPt(0), fHistJetPt(0),
  fHistHighJetPt(0),
  fHistnSigElecPt(0),
  fHistnJetTrackvnJetClusters(0),
  fHistPtDEtaDPhiTrackClus(0),
  fHistPtDEtaDPhiClusTrack(0),
  fhnPIDHF(0x0), fhnJetQA(0x0), fhnClusterTrackQA(0x0), fhnTrackClusterQA(0x0), fhnPIDHFTtoC(0x0)
{
  // Default constructor.
  for (Int_t i = 0;i<6;++i){
    fHistJetPtvsTrackPt[i]      = 0;
    fHistTrackPt[i]             = 0;
    fHistEP0[i]                 = 0;
    fHistEP0A[i]                = 0;
    fHistEP0C[i]                = 0;
    fHistEPAvsC[i]              = 0;
  }

  SetMakeGeneralHistograms(kTRUE);

}

//________________________________________________________________________
AliAnalysisTaskEmcalJetHF::AliAnalysisTaskEmcalJetHF(const char *name) :
  AliAnalysisTaskEmcalJet(name,kTRUE),
  event(0),
  fillHist(0),
  fEventTrigEMCALL1Gamma1(0),
  fEventTrigEMCALL1Gamma2(0),
  fGlobalQA(0),
  fInputEvent(0x0),
  fCuts(0),
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0),
  fJetHIpt(20.0),
  fTrackPtCut(2.0),
  fTrackEta(0.9),
  fTrkQAcut(0),
  fesdTrackCuts(0),
  fPIDResponse(0x0), fTPCResponse(),
  fEsdtrackCutsITSTPC(),
  fEsdtrackCutsTPC(),
  fEsdtrackCutsITS(),
  fJetsCont(0), fTracksCont(0), fCaloClustersCont(0), fTracksJetCont(0), fCaloClustersJetCont(0),
  fESD(0), fAOD(0),
  fHistRhovsCent(0),
  fHistJetPhi(0),
  fHistCorJetPt(0), fHistJetPt(0),
  fHistHighJetPt(0),
  fHistnSigElecPt(0),
  fHistnJetTrackvnJetClusters(0),
  fHistPtDEtaDPhiTrackClus(0),
  fHistPtDEtaDPhiClusTrack(0),
  fhnPIDHF(0x0), fhnJetQA(0x0), fhnClusterTrackQA(0x0), fhnTrackClusterQA(0x0), fhnPIDHFTtoC(0x0)
{ 
  for (Int_t i = 0;i<6;++i){
    fHistJetPtvsTrackPt[i]      = 0;
    fHistTrackPt[i]             = 0;
    fHistEP0[i]                 = 0;
    fHistEP0A[i]                = 0;
    fHistEP0C[i]                = 0;
    fHistEPAvsC[i]              = 0;
   }
   SetMakeGeneralHistograms(kTRUE);
 
   DefineInput(0,TChain::Class());
   DefineOutput(1, TList::Class());
}

//_______________________________________________________________________
AliAnalysisTaskEmcalJetHF::~AliAnalysisTaskEmcalJetHF()
{
  // destructor
  //
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetHF::UserCreateOutputObjects()
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

  fillHist = 1;
  TString histname;

  if(fillHist>0){
  fHistRhovsCent              = new TH2F("RhovsCent", "RhovsCent", 100, 0.0, 100.0, 400, 0, 400);
  fHistCorJetPt		            = new TH1F("NjetvsCorrJetPt", "NjetvsCorrJetPt", 300, -100, 200);
  fHistnSigElecPt             = new TH2F("nsig_v_pt(TPC)","nsig_v_pt(TPC)",200,0,100,100,-10,10);
  fHistnJetTrackvnJetClusters = new TH2F("NumbJetTracksvJetClusters","NumbJetTracksvJetClusters",21,0,20,21,0,20);
  fHistHighJetPt              = new TH1F("HighestPtJetPerEvent","HighJetPt",300,0,150);
    
  histname = "fHistPtDEtaDPhiTrackClus";
  fHistPtDEtaDPhiTrackClus = new TH3F(histname.Data(),Form("%s;#it{p}_{T}^{track};#Delta#eta;#Delta#varphi",histname.Data()),100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiTrackClus);
    
  histname = "fHistPtDEtaDPhiClusTrack";
  fHistPtDEtaDPhiClusTrack = new TH3F(histname.Data(),Form("%s;#it{p}_{T}^{clus};#Delta#eta;#Delta#varphi",histname.Data()),100,0.,100.,100,-0.1,0.1,100,-0.1,0.1);
  fOutput->Add(fHistPtDEtaDPhiClusTrack);
    
  // PT bins used to be (2000, -100, 300) 
  TString name;
  TString title;

  // creating centrality dependent histos that don't involve Global Rho
  for (Int_t i = 0;i<6;++i){
    name = TString(Form("JetPtvsTrackPt_%i",i));
    title = TString(Form("Jet pT vs Leading Track pT cent bin %i",i));
    fHistJetPtvsTrackPt[i] = new TH2F(name,title, 500, -100, 400, 100,0,100);
    fOutput->Add(fHistJetPtvsTrackPt[i]);

    name = TString(Form("TrackPt_%i",i));
    title = TString(Form("Track pT cent bin %i",i));
    fHistTrackPt[i] = new TH1F(name,title,400,0,200);
    fOutput->Add(fHistTrackPt[i]);   

    name = TString(Form("EP0_%i",i));
    title = TString(Form("EP VZero cent bin %i",i));
    fHistEP0[i] = new TH1F(name,title,144,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEP0[i]);

    name = TString(Form("EP0A_%i",i));
    title = TString(Form("EP VZero cent bin %i",i));
    fHistEP0A[i] = new TH1F(name,title,144,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEP0A[i]);

    name = TString(Form("EP0C_%i",i));
    title = TString(Form("EP VZero cent bin %i",i));
    fHistEP0C[i] = new TH1F(name,title,144,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEP0C[i]);

    name = TString(Form("EPAvsC_%i",i));
    title = TString(Form("EP VZero cent bin %i",i));
    fHistEPAvsC[i] = new TH2F(name,title,144,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
    fOutput->Add(fHistEPAvsC[i]);

  }

  fOutput->Add(fHistRhovsCent);
  fOutput->Add(fHistCorJetPt);
  fOutput->Add(fHistnSigElecPt);
  fOutput->Add(fHistnJetTrackvnJetClusters);
  fOutput->Add(fHistHighJetPt);
  }//Fill Histograms

  // ****************************** PID *****************************************************                                               
  // set up PID handler                                                                                                                     
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if(!inputHandler) {
    AliFatal("Input handler needed");
    return;
  }

  // PID response object                                                                                                                    
  //fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();                                                                         
  //  inputHandler->CreatePIDResponse(fIsMC);         // needed to create object, why though?                                                 
  fPIDResponse = inputHandler->GetPIDResponse();
  if (!fPIDResponse) {
    AliError("PIDResponse object was not created");
    return;
  }
  // ****************************************************************************************
  UInt_t bitcoded = 0;  // bit coded, see GetDimParamsPID() below
  bitcoded = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8 | 1<<9 | 1<<10 | 1<<11 | 1<<12 | 1<<13| 1<<14 | 1<<15 | 1<<16 | 1<<17;
  fhnPIDHF = NewTHnSparseDHF("fhnPIDHFCtoT", bitcoded);
  
  UInt_t bitcoded1 = 0;  // bit coded, see GetDimParamsPID() below
  bitcoded1 = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4;
  fhnJetQA = NewTHnSparseDJetQA("fhnJetQA", bitcoded1);
  
  UInt_t bitcoded2 = 0;  // bit coded, see GetDimParamsPID() below
  bitcoded2 = 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<8 | 1<<9 | 1<<10 | 1<<15 | 1<<16 | 1<<17;
  fhnClusterTrackQA = NewTHnSparseDHF("fhnClusterTrackQA", bitcoded2);
  
  UInt_t bitcoded3 = 0;  // bit coded, see GetDimParamsPID() below
  bitcoded3 = 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<8 | 1<<9 | 1<<10 | 1<<15 | 1<<16 | 1<<17;
  fhnTrackClusterQA = NewTHnSparseDHF("fhnTrackClusterQA", bitcoded3);
  
  UInt_t bitcoded7 = 0;  // bit coded, see GetDimParamsPID() below
  bitcoded7 = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8 | 1<<9 | 1<<10 | 1<<11 | 1<<12 | 1<<13| 1<<14 | 1<<15 | 1<<16 | 1<<17;
  fhnPIDHFTtoC = NewTHnSparseDHF("fhnPIDHFTtoC", bitcoded7);
  
  cout << "_______________Created Sparse__________________" << endl;
  
  fOutput->Add(fhnPIDHF);
  fOutput->Add(fhnJetQA);
  fOutput->Add(fhnClusterTrackQA);
  fOutput->Add(fhnTrackClusterQA);
  fOutput->Add(fhnPIDHFTtoC);

  PostData(1, fOutput);

}

//________________________________________________________
void AliAnalysisTaskEmcalJetHF::ExecOnce()
{
  //  Initialize the analysis
  AliAnalysisTaskEmcalJet::ExecOnce();
  
  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;


} // end of ExecOnce

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetHF::Run()
{
  // check to see if we have any tracks
  if (!fTracks)  return kTRUE;
  if (!fJets)  return kTRUE;

  // what kind of event do we have: AOD or ESD?
  Bool_t useAOD;
  if (dynamic_cast<AliAODEvent*>(InputEvent())) useAOD = kTRUE;
  else useAOD = kFALSE;
  
  fEventTrigEMCALL1Gamma1 = kFALSE;
  fEventTrigEMCALL1Gamma2 = kFALSE;
  
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

  // get centrality bin
  Int_t centbin = GetCentBin(fCent);
  //for pp analyses we will just use the first centrality bin
  if (centbin == -1)  centbin = 0;

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
  //TClonesArray * clusterList = 0;
  
  // get Jets object                                                                                                                        
  jets = dynamic_cast<TClonesArray*>(list->FindObject(fJets));
  if(!jets){
    AliError(Form("Pointer to jets %s == 0", fJets->GetName()));
    return kTRUE;
  } // verify existence of jets
  
  event++;
  //cout<<"Event #: "<<event<<"  Number of Clusters: "<<fCaloClustersCont->GetNClusters()<<"  Number of Tracks: "<<fTracksCont->GetNParticles()<<endl;
  
 // Int_t   nclusters   = fInputEvent->GetNumberOfCaloClusters();
  
  
  // get number of jets and tracks                                                                                                          
  const Int_t Njets = jets->GetEntries();
  if(Njets<1)     return kTRUE;
  
  if (fTracksCont) {
    AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle(0));
    while(track) {
      track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    }
  }
  if (fCaloClustersCont) {
    AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(0);
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
      
      // phi of jet, constrained to 1.6 < Phi < 2.94
      float jetphi = jet->Phi();      // phi of jet
      // apply jet cuts
      if(!AcceptMyJet(jet)) continue;
      
      //AliEmcalJet::EFlavourTag tag=AliEmcalJet::kDStar;
      //jet->AddFlavourTag(tag);
 
      //MV: removed to avoid compiler warnings
      //      Bool_t bkgrnd1  = kFALSE;
      // Bool_t sig1     = kFALSE;
 
      Int_t JetClusters = jet->GetNumberOfClusters();
      Int_t JetTracks = jet -> GetNumberOfTracks();
      fHistnJetTrackvnJetClusters->Fill(JetClusters,JetTracks);
     // Initializations and Calculations
     Double_t jetptraw = jet->Pt();    				// raw pT of jet
     Double_t jetPt = -500;                                     // initialize corr jet pt LOCAL
     Double_t jetarea = -500;					// initialize jet area
     jetarea = jet->Area();		           		// jet area
     jetPt = jet->Pt() - jetarea*fRhoVal;                 // semi-corrected pT of jet from GLOBAL rho value
      fHistCorJetPt->Fill(jetPt);
      
      if(jet->Pt() > fJetHIpt) {
        if(!fTracksCont || !fCaloClustersCont) return kTRUE;
        
        Double_t dEdx = -99;
        Double_t EovP = -99;
        Float_t DCAxy = -999;
        Float_t DCAz = -999;
        
        Double_t deta = 999;
        Double_t dphi = 999;
        Int_t NumbCluster = -999;
        NumbCluster = fCaloClustersCont->GetNClusters();
        Double_t JetQA[5] = {static_cast<Double_t>(Njets), static_cast<Double_t>(jet->GetNumberOfTracks()), static_cast<Double_t>(jet->GetNumberOfClusters()),jet->Eta(), jet->Phi()};
        fhnJetQA->Fill(JetQA);
        
      //***********************************************
      //****************Track Matched to Closest Cluster
        
        
       for(int iCluster = 0; iCluster <= NumbCluster; iCluster++){
          //Get closest track to cluster to track matching!!!!!
          //AliVCluster *cluster = fCaloClustersCont->GetNextAcceptedCluster(iCluster);
          AliVCluster *cluster = fCaloClustersCont->GetCluster(iCluster);
          
          if(! IsJetCluster(jet, iCluster, kFALSE)) continue;
         
          // while(cluster) {
          TLorentzVector nPart;
          cluster->GetMomentum(nPart, fVertex);
          //fHistClustersPt[fCentBin]->Fill(nPart.Pt());
          Double_t fclusE = -999;
          fclusE = cluster->E();
          
          Float_t pos[3];
          cluster->GetPosition(pos);  // Get cluster position
          TVector3 cp(pos);
          
          //Get matched track
          AliVTrack *mt = NULL;
          //AliESDtrack *mt = NULL;
          AliESDtrack *ESDmt = NULL;
          AliAODCaloCluster *acl = dynamic_cast<AliAODCaloCluster*>(cluster);
          if(acl) {
            
            //mt = static_cast<AliVTrack*>(acl->GetTrackMatched(0));
          }
          else {
            
            AliESDCaloCluster *ecl = dynamic_cast<AliESDCaloCluster*>(cluster);
            Int_t im = ecl->GetTrackMatchedIndex();
            
            if(fTracksCont && im>=0) {
              //mt = static_cast<AliVTrack*>(fTracksCont->GetParticle(im));
              mt = static_cast<AliVTrack*>(fTracksCont->GetParticle(im));
              if(!mt) continue;
              
              ESDmt = static_cast<AliESDtrack*>(mt);
              if(!ESDmt) continue;
              
              Double_t pcluster = mt->P();
              //Double_t esdp = ESDmt->P();
              //Int_t LabelNumb, IDNumb;
              //LabelNumb = ESDmt->GetLabel();
              //IDNumb= ESDmt -> GetID();

              dEdx = mt->GetTPCsignal();
              Double_t p = mt->P();
               //TPC nSigma's
              //nsigpion = fPIDResponse->NumberOfSigmasTPC(mt,AliPID::kPion);
              Double_t nSigmaElectron_TPC = fPIDResponse->NumberOfSigmasTPC(mt,AliPID::kElectron);
              Double_t nSigmaElectron_TOF = fPIDResponse->NumberOfSigmasTOF(mt,AliPID::kElectron);
              dEdx = mt->GetTPCsignal();
              ESDmt->GetImpactParameters(DCAxy, DCAz);
              EovP        = fclusE/p;
                Double_t HF_tracks[18] = {fCent, mt->Pt(), pcluster ,mt->Eta(), mt->Phi(), EovP, DCAxy, DCAz, dEdx,nSigmaElectron_TPC, nSigmaElectron_TOF,0 /*nSigmaElectron_EMCAL*/, jetptraw, jetphi, jet->Eta(),fclusE,cp.PseudoRapidity(),cp.Phi()};
                fhnPIDHF->Fill(HF_tracks);    // fill Sparse Histo with trigger entries

       
              
            }
          }
          
          if(mt) {
            AliPicoTrack::GetEtaPhiDiff(mt, cluster, dphi, deta);
            fHistPtDEtaDPhiClusTrack->Fill(nPart.Pt(),deta,dphi);
            
            
            
          }

          //cluster = fCaloClustersCont->GetNextAcceptCluster();
          //if(! IsJetCluster(jet, cluster, kFALSE)) continue;
        //}
          
          /*
           Double_t p = mt->P();
           Float_t pos[3];
           cluster->GetPosition(pos);  // Get cluster position
           //AliESDtrack *trackESD = fESD->GetTrack(Ntracks);
           
           //  nSigmaElectron_TPC = fPIDResponse->NumberOfSigmasTPC(mt,AliPID::kElectron);
           //  nSigmaElectron_TOF= fPIDResponse->NumberOfSigmasTOF(mt,AliPID::kElectron);
           EovP        = fclusE/p;
           TVector3 cp(pos);

           //if(!fesdTrackCuts->AcceptTrack(mt)) continue;
           //dEdx = mt->GetTPCsignal();
           //mt->GetImpactParameters(DCAxy, DCAz);
           //nSigmaElectron_TPC = fPIDResponse->NumberOfSigmasTPC(trackESD,AliPID::kElectron);
           //nSigmaElectron_TOF = fPIDResponse->NumberOfSigmasTOF(trackESD,AliPID::kElectron);
           
           
           //Double_t HF_tracks[18] = {fCent, mt->Pt(), mt->P() ,mt->Eta(), mt->Phi(), EovP, DCAxy, DCAz, dEdx, nSigmaElectron_TPC, nSigmaElectron_TOF, nSigmaElectron_EMCAL, jet->Pt(), jet->Phi(), jet->Eta(),fclusE,cp.PseudoRapidity(),cp.Phi()};
           //fhnPIDHF->Fill(HF_tracks);    // fill Sparse Histo with trigger entries
           */
          //cluster = fCaloClustersCont->GetNextAcceptCluster();
          //if(! IsJetCluster(jet, cluster, kFALSE)) continue;

          
          //AliESDtrack *ESDacceptedTrack = NULL;
        }//loop over cluster

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
        Double_t acceptTrackEta = AcceptedTrack->Eta();
        Double_t acceptTrackPhi = AcceptedTrack->Phi();
        Double_t nSigmaElectron_TPC_at = fPIDResponse->NumberOfSigmasTPC(AcceptedTrack,AliPID::kElectron);
        Double_t nSigmaElectron_TOF_at = fPIDResponse->NumberOfSigmasTOF(AcceptedTrack,AliPID::kElectron);

        
        AliESDtrack *ESDacceptedTrack = static_cast<AliESDtrack*>(AcceptedTrack);
       
        if(!ESDacceptedTrack){
          AliError(Form("Couldn't get AliESDTrack %d\n", iTracks));
          continue;
        }
        //Double_t DCAxy_at, DCAz_at;
        Double_t dEdxat = AcceptedTrack->GetTPCsignal();
        //ESDacceptedTrack->GetImpactParameters(DCAxy_at, DCAz_at);
        
        
        if(fCaloClustersCont && emc1>=0) {
          AliVCluster *clusMatch = fCaloClustersCont->GetCluster(emc1);
          if(!clusMatch){
            AliError(Form("Couldn't get matched AliVCluster %d\n", emc1));
            continue;
          }
 
          
          
          Double_t mClusterE = clusMatch->E();
          Float_t pos_mc[3];
          clusMatch->GetPosition(pos_mc);  // Get cluster position
          TVector3 mcp(pos_mc);
          Double_t EovP_mc = -999;
          EovP_mc = mClusterE/acceptTrackP;
	  //MV: removed to avoid compiler warnings
          // if(EovP_mc < 0.2){
          //   bkgrnd1 = kTRUE;      //Hadron Background
          // }

	  //Code without meaning:
          //if(0.8 < EovP_mc < 1.2){
	  //            if(-1.5<nSigmaElectron_TPC_at<5.0){
	  //              if(4.0<acceptTrackPt<10.0){

	  //Corrected code:
	  // if(EovP_mc >0.8 && EovP_mc<1.2){ 
          //   if(nSigmaElectron_TPC_at>-1.5 && nSigmaElectron_TPC_at<5.0){
          //     if(acceptTrackPt>4.0 && acceptTrackPt<10.0){
          //       sig1 = kTRUE;                    //Electron Candidate
          //     }
          //   }
          // }
          
          Double_t HF_tracks2[18] = {fCent, acceptTrackPt, acceptTrackP ,acceptTrackEta, acceptTrackPhi, EovP_mc, 0, 0, dEdxat,nSigmaElectron_TPC_at, nSigmaElectron_TOF_at,0 , jetPt, jet->Phi(), jet->Eta(),mClusterE,mcp.PseudoRapidity(),mcp.Phi()};
          fhnPIDHFTtoC->Fill(HF_tracks2);    // fill Sparse Histo with trigger entries
        }
        //AcceptedTrack = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
        
      } //loop over tracks for matching to closest cluster

     

        
    } // highest pt jet cut
      /*
      Int_t tag = -999;
      if(bkgrnd1 == kTRUE) {
        AliEmcalJet::EFlavourTag tag=AliEmcalJet::kBckgrd1;
        jet->AddFlavourTag(tag);
      }
      if(sig1 == kTRUE && !bkgrnd1){
        AliEmcalJet::EFlavourTag tag=AliEmcalJet::kSig1;
        jet->AddFlavourTag(tag);
      }
       */
  } // LOOP over JETS in event

  
  
  if(fGlobalQA == 1) CheckClusTrackMatchingQA();

  return kTRUE;
  
}
//________________________________________________________________________
void AliAnalysisTaskEmcalJetHF::Terminate(Option_t *)
{
  cout<<"###########################"<<endl;
  cout<<"####   Task Finished   ####"<<endl;
  cout<<"###########################"<<endl;
  cout<<"###########################"<<endl;
} // end of terminate


//________________________________________________________________________
void AliAnalysisTaskEmcalJetHF::CheckClusTrackMatchingQA()
{
  
  if(!fTracksCont || !fCaloClustersCont) return;
  
  Int_t trkcounter = 0;
  Int_t cluscounter = 0;
  Double_t deta = 999;
  Double_t dphi = 999;
  
  //Get closest cluster to track
  AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle(0));
  while(track) {
    //if(!track) continue;
    AliESDtrack *ESDtrackQA1 = static_cast<AliESDtrack*>(track);
    if(!ESDtrackQA1) continue;
    if(!fPIDResponse) continue;
    Double_t pQA1 = track->P();
    Double_t nSigmaElectron_TPC_QA1 = fPIDResponse->NumberOfSigmasTPC(ESDtrackQA1,AliPID::kElectron);
    Double_t nSigmaElectron_TOF_QA1 = fPIDResponse->NumberOfSigmasTOF(ESDtrackQA1,AliPID::kElectron);
    Double_t dEdxQA1 = ESDtrackQA1->GetTPCsignal();
    //Get matched cluster
    Int_t emc1 = track->GetEMCALcluster();
    if(fCaloClustersCont && emc1>=0) {
      AliVCluster *clusMatch = fCaloClustersCont->GetCluster(emc1);
      //if(!clusMatch) continue;
      if(clusMatch) {
        Double_t ClusterE_QA1 = clusMatch->E();
        Double_t EovPQA1 = ClusterE_QA1/pQA1;
        Float_t pos_mc1[3];
        clusMatch->GetPosition(pos_mc1);  // Get cluster position
        TVector3 mc1(pos_mc1);
        fHistnSigElecPt->Fill(nSigmaElectron_TPC_QA1,track->Pt());
        Double_t HF_tracks3[11] = {track->Pt(), track->P() , track->Eta(), track->Phi(), EovPQA1, dEdxQA1 ,nSigmaElectron_TPC_QA1, nSigmaElectron_TOF_QA1, clusMatch->E(), mc1.PseudoRapidity(),mc1.Phi()};
        fhnTrackClusterQA->Fill(HF_tracks3);
        AliPicoTrack::GetEtaPhiDiff(track, clusMatch, dphi, deta);
        fHistPtDEtaDPhiTrackClus->Fill(track->Pt(),deta,dphi);
        cluscounter++;
      }//clus matching

     }//matched cluster

    track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    trkcounter++;
  }//track loop
 /*
  //Get closest track to cluster
  AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(0);
  while(cluster) {
    if(!cluster) continue;
    
    Double_t ClusterE_QA2 = cluster->E();
    //Double_t EovPQA1 = ClusterE_QA1/pQA1;
    Float_t pos_mc2[3];
    cluster->GetPosition(pos_mc2);  // Get cluster position
    TVector3 mc2(pos_mc2);

    TLorentzVector nPart;
    cluster->GetMomentum(nPart, fVertex);
    
    //Get matched track
    AliVTrack *mt = NULL;
    AliAODCaloCluster *acl = dynamic_cast<AliAODCaloCluster*>(cluster);
    if(!acl) continue;
    if(acl) {
      if(acl->GetNTracksMatched()>1)
      mt = static_cast<AliVTrack*>(acl->GetTrackMatched(0));
      if(!mt) continue;
      AliESDtrack *ESDtrackQA2 = static_cast<AliESDtrack*>(mt);
      if(!ESDtrackQA2) continue;
      Double_t nSigmaElectron_TPC_QA2 = fPIDResponse->NumberOfSigmasTPC(ESDtrackQA1,AliPID::kElectron);
      Double_t nSigmaElectron_TOF_QA2 = fPIDResponse->NumberOfSigmasTOF(ESDtrackQA1,AliPID::kElectron);
      Double_t dEdxQA1 = ESDtrackQA2->GetTPCsignal();
      Double_t EovPQA2 = -999;
      Double_t pQA2 = mt->P();
      EovPQA2 = ClusterE_QA2/pQA2;
  
      //Double_t HF_tracks4[11] = {mt->Pt(), mt->P() , mt->Eta(), mt->Phi(), EovPQA2, dEdxQA2 ,nSigmaElectron_TPC_QA2, nSigmaElectron_TOF_QA2, mc2.PseudoRapidity(),mc2.Phi()};
      //fhnClusterTrackQA->Fill(HF_tracks4);
    }
    else {
      AliESDCaloCluster *ecl = dynamic_cast<AliESDCaloCluster*>(cluster);
      if(!ecl) continue;
      Int_t im = ecl->GetTrackMatchedIndex();
      if(fTracksCont && im>=0) {
        mt = static_cast<AliVTrack*>(fTracksCont->GetParticle(im));
        if(!mt) continue
      }
    }
    AliPicoTrack::GetEtaPhiDiff(mt, cluster, dphi, deta);
    fHistPtDEtaDPhiClusTrack->Fill(nPart.Pt(),deta,dphi);
    cluster = fCaloClustersCont->GetNextAcceptCluster();
  }
  }
*/
/*
  //Get closest track to cluster
  AliVCluster *cluster = fCaloClustersCont->GetNextAcceptCluster(0);
  while(cluster) {
    if(!cluster){
      AliError(Form("Couldn't get CtoT AliCluster Container"));
      continue;
    }
    Double_t ClusterE_QA2 = cluster->E();
    Float_t pos_mc2[3];
    cluster->GetPosition(pos_mc2);  // Get cluster position
    TVector3 mc2(pos_mc2);
    TLorentzVector nPart;
    cluster->GetMomentum(nPart, fVertex);
    //Get matched track
    AliVTrack *mt = NULL;
    AliESDCaloCluster *ecl = dynamic_cast<AliESDCaloCluster*>(cluster);
    if(!ecl){
      AliError(Form("Couldn't get CtoT AliESDCluster Container %d\n"));
      continue;
    }
    Int_t im = ecl->GetTrackMatchedIndex();
    if(fTracksCont && im>=0) {
      //mt = static_cast<AliVTrack*>(fTracksCont->GetParticle(im));
      //if(!mt){
        //AliError(Form("Couldn't get CT AliVTrack Container %d\n",im));
        //continue;
     // }
      //AliPicoTrack::GetEtaPhiDiff(mt, cluster, dphi, deta);
      fHistPtDEtaDPhiClusTrack->Fill(nPart.Pt(),deta,dphi);
      AliESDtrack *ESDtrackQA2 = static_cast<AliESDtrack*>(mt);
      if(!ESDtrackQA2){
        AliError(Form("Couldn't get CT AliESDTrack Container %d\n"));
        continue;
      }
      Double_t nSigmaElectron_TPC_QA2 = fPIDResponse->NumberOfSigmasTPC(ESDtrackQA2,AliPID::kElectron);
      Double_t nSigmaElectron_TOF_QA2 = fPIDResponse->NumberOfSigmasTOF(ESDtrackQA2,AliPID::kElectron);
      Double_t dEdxQA1 = ESDtrackQA2->GetTPCsignal();
      Double_t EovPQA2 = -999;
      Double_t pQA2 = mt->P();
      EovPQA2 = ClusterE_QA2/pQA2;


    }
    //cluster = fCaloClustersCont->GetNextAcceptCluster();
    cluster = static_cast<AliVCluster*>(fCaloClustersCont->GetNextAcceptCluster());
  }
*/
}


//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHF::AcceptMyJet(AliEmcalJet *jet) {
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

//________________________________________________________________________
Int_t AliAnalysisTaskEmcalJetHF::GetCentBin(Double_t cent) const
{  // Get centrality bin.

  Int_t centbin = -1;
  if (cent>=0 && cent<10)
    centbin = 0; 
  else if (cent>=10 && cent<20)
    centbin = 1;
  else if (cent>=20 && cent<30)
    centbin = 2;
  else if (cent>=30 && cent<40)
    centbin = 3;
  else if (cent>=40 && cent<50)
    centbin = 4;
  else if (cent>=50 && cent<90)
    centbin = 5;
  return centbin;
}

//________________________________________________________________________________
/*
void AliAnalysisTaskEmcalJetHF::FlagFlavour(AliEmcalJet *jet){
  
  AliEmcalJet::EFlavourTag tag=AliEmcalJet::kDStar;
  if (fCandidateType==kD0toKpi) tag=AliEmcalJet::kD0;
  if (fIsDInJet) jet->AddFlavourTag(tag);
  
  return;
  
}
*/
//____________________________________________________________________________________________
THnSparse* AliAnalysisTaskEmcalJetHF::NewTHnSparseDHF(const char* name, UInt_t entries)
{
  // generate new THnSparseD PID, axes are defined in GetDimParams()                                                                                                     
  Int_t count = 0;
  UInt_t tmp = entries;
  while(tmp!=0){
    count++;
    tmp = tmp &~ -tmp;  // clear lowest bit                                                                                                                             
  }

  TString hnTitle(name);
  const Int_t dim = count;
  Int_t nbins[dim];
  Double_t xmin[dim];
  Double_t xmax[dim];

  Int_t i=0;
  Int_t c=0;
  while(c<dim && i<32){
    if(entries&(1<<i)){

      TString label("");
      GetDimParamsHF(i, label, nbins[c], xmin[c], xmax[c]);
      hnTitle += Form(";%s",label.Data());
      c++;
    }

    i++;
  }
  hnTitle += ";";

  return new THnSparseD(name, hnTitle.Data(), dim, nbins, xmin, xmax);
} // end of NewTHnSparseF PID

THnSparse* AliAnalysisTaskEmcalJetHF::NewTHnSparseDJetQA(const char* name, UInt_t entries)
{
  // generate new THnSparseD JetQA, axes are defined in GetDimParamsJetQA()
  Int_t count = 0;
  UInt_t tmp = entries;
  while(tmp!=0){
    count++;
    tmp = tmp &~ -tmp;  // clear lowest bit
  }
  
  TString hnTitle(name);
  const Int_t dim = count;
  Int_t nbins[dim];
  Double_t xmin[dim];
  Double_t xmax[dim];
  
  Int_t i=0;
  Int_t c=0;
  while(c<dim && i<32){
    if(entries&(1<<i)){
      
      TString label("");
      GetDimParamsJetQA(i, label, nbins[c], xmin[c], xmax[c]);
      hnTitle += Form(";%s",label.Data());
      c++;
    }
    
    i++;
  }
  hnTitle += ";";
  
  return new THnSparseD(name, hnTitle.Data(), dim, nbins, xmin, xmax);
} // end of NewTHnSparseF JetQA


void AliAnalysisTaskEmcalJetHF::GetDimParamsHF(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
  // stores label and binning of axis for THnSparse                                                                                                                      
  const Double_t pi = TMath::Pi();

  switch(iEntry){

  case 0:
    label = "V0 centrality (%)";
    nbins = 10;
    xmin = 0.;
    xmax = 100.;
    break;

  case 1:
    label = "Track p_{T}";
    nbins = 300;
    xmin = 0.;
    xmax = 75.;
    break;

  case 2:
    label = "Track p";
    nbins = 300;
    xmin = 0.;
    xmax = 75.;
    break;

  case 3:
    label = "Track Eta";
    nbins = 48;
    xmin = -1.2;
    xmax = 1.2;
    break;

  case 4:
    label = "Track Phi";
    nbins = 72;
    xmin = 0;
    xmax = 2*pi;
    break;

  case 5:
    label = "E/p of track";
    nbins = 400;
    xmin = 0;
    xmax = 4.0;
    break;

 case 6:
    label = "DCA xy";
    nbins = 20;
    xmin = -10;
    xmax =  10;
    break;

  case 7:
    label = "DCA z";
    nbins = 20;
    xmin = -10;
    xmax = 10;
    break;

  case 8:                                                                                                                                               
    label = "dEdX of track - TPC";
    nbins = 300;
    xmin = 0;
    xmax = 300;
    break;

  case 9:                                                                                                                                                
    label = "nSigma electron TPC";
    nbins = 50;
    xmin = -5;
    xmax = 5;
    break;

   case 10:
    label = "nSigma electron TOF";
    nbins = 50;
    xmin = -5;
    xmax = 5;
    break;

   case 11:
    label = "nSigma electron Emcal";
    nbins = 50;
    xmin = -5;
    xmax = 5;
    break;
      
  case 12:
    label = "Jet pT";
    nbins = 40;
    xmin  = 0;
    xmax  = 200;
    break;
      
  case 13:
    label = "Jet Phi";
    nbins = 72;
    xmin = 0;
    xmax = 2*pi;
    break;
      
  case 14:
    label = "Jet Eta";
    nbins = 24;
    xmin = -1.2;
    xmax = 1.2;
    break;
      
  case 15:
    label = "Cluster Energy";
    nbins = 150;
    xmin = 0;
    xmax = 15;
    break;
      
  case 16:
    label = "Cluster Eta";
    nbins = 24;
    xmin = -1.2;
    xmax =  1.2;
    break;
      
  case 17:
    label = "Cluster Phi";
    nbins = 72;
    xmin = 0;
    xmax = 2*pi;
    break;



  } // end of switch
} // end of getting dim-params

void AliAnalysisTaskEmcalJetHF::GetDimParamsJetQA(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
  // stores label and binning of axis for THnSparse
  const Double_t pi = TMath::Pi();
  
  switch(iEntry){
      
    case 0:
      label = "number of Jets in Event";
      nbins = 100;
      xmin = 0.;
      xmax = 100.;
      break;
      
    case 1:
      label = "number of Clusters in a Jet";
      nbins = 100;
      xmin = 0.;
      xmax = 100.;
      break;
      
    case 2:
      label = "number of Tracks in a Jet";
      nbins = 100;
      xmin = 0.;
      xmax = 100.;
      break;
      
    case 3:
      label = "Jet Eta";
      nbins = 24;
      xmin = -1.2;
      xmax = 1.2;
      break;
      
    case 4:
      label = "Jet Phi";
      nbins = 72;
      xmin = 0;
      xmax = 2*pi;
      break;
      
    case 5:
      label = "Cluster Energy";
      nbins = 150;
      xmin = 0;
      xmax = 15;
      break;
      
    case 6:
      label = "Cluster Eta";
      nbins = 24;
      xmin = -1.2;
      xmax =  1.2;
      break;
      
    case 7:
      label = "Cluster Phi";
      nbins = 72;
      xmin = 0;
      xmax = 2*pi;
      break;
      
    case 8:
      label = "Is EMCalCluster";
      nbins = 2;
      xmin = 0;
      xmax = 2;
      break;
      
    case 9:
      label = "Number of Matched Tracks to Cluster";
      nbins = 60;
      xmin = 0;
      xmax = 60;
      break;
      
    case 10:
      label = "Track Pt";
      nbins = 300;
      xmin = 0;
      xmax = 75;
      break;
      
    case 11:
      label = "Track Eta";
      nbins = 24;
      xmin = -1.2;
      xmax = 1.2;
      break;
      
    case 12:
      label= "Track Phi";
      nbins = 72;
      xmin = 0;
      xmax = 2*pi;
      break;
      
    case 13:
      label="Is Track EMCal";
      nbins = 2;
      xmin = 0;
      xmax = 2;
      break;
      
    case 14:
      label = "Get Track EMCal Cluster";
      nbins = 100;
      xmin = 0;
      xmax = 100;
      break;
      
    case 15:
      label = "Track Matched Phi";
      nbins = 72;
      xmin = 0;
      xmax = 2*pi;
      
      
      
      
  } // end of switch
} // end of getting dim-params



