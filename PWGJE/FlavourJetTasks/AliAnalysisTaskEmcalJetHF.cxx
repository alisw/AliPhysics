//
// Author: A Castro

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
#include "AliEmcalParticle.h"

// event handler (and pico's) includes                                                                                                      
#include <AliInputEventHandler.h>
#include <AliVEventHandler.h>
#include "AliESDInputHandler.h"
#include "AliPicoTrack.h"
#include "AliEventPoolManager.h"

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
  //isESD(0),
  doGlobalPID(0),
  fCuts(0),
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0),
  fJetHIpt(50.0),
  fTrackPtCut(2.0),
  fTrackEta(0.9),
  fTrkQAcut(0),
  fPIDResponse(0x0), fTPCResponse(),
  fEsdtrackCutsITSTPC(),
  fEsdtrackCutsTPC(),
  fEsdtrackCutsITS(),
  fESD(0), fAOD(0),
  fHistRhovsCent(0),
  fHistJetPhi(0),
  fHistCorJetPt(0), fHistJetPt(0),
  fHistdEdx(0),
  fHistdEdxvPt(0),
  fHistClusE(0),
  fHistEovPTracks(0),
  fHistEovPvsPtTracks(0),
  fHistPID(0), fHistPIDtpc(0), fHistPIDits(0), fHistPIDtof(0),
  fHistnsigelectron(0),
  fHistnSigElecPt(0),
  fHistTrackPhivEta(0),
  fHistClusterPhivEta(0),
  fHistnJetTrackvnJetClusters(0),
  fhnPIDHF(0x0), fhnQA(0x0), fhnJetQA(0x0), fhnClusQA(0x0), fhnTrackQA(0x0), fhnGlobalPID(0x0)
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
  //isESD(0),
  doGlobalPID(0),
  fCuts(0),
  fPhimin(-10), fPhimax(10),
  fEtamin(-0.9), fEtamax(0.9),
  fAreacut(0.0),
  fJetHIpt(50.0),
  fTrackPtCut(2.0),
  fTrackEta(0.9),
  fTrkQAcut(0),
  fPIDResponse(0x0), fTPCResponse(),
  fEsdtrackCutsITSTPC(),
  fEsdtrackCutsTPC(),
  fEsdtrackCutsITS(),
  fESD(0), fAOD(0),
  fHistRhovsCent(0),
  fHistJetPhi(0),
  fHistCorJetPt(0), fHistJetPt(0),
  fHistdEdx(0),
  fHistdEdxvPt(0),
  fHistClusE(0),
  fHistEovPTracks(0),
  fHistEovPvsPtTracks(0),
  fHistPID(0), fHistPIDtpc(0), fHistPIDits(0), fHistPIDtof(0),
  fHistnsigelectron(0),
  fHistnSigElecPt(0),
  fHistTrackPhivEta(0),
  fHistClusterPhivEta(0),
  fHistnJetTrackvnJetClusters(0),
  fhnPIDHF(0x0), fhnQA(0x0), fhnJetQA(0x0), fhnClusQA(0x0), fhnTrackQA(0x0), fhnGlobalPID(0x0)
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

  fOutput = new TList();
  fOutput->SetOwner(kTRUE);

  fHistJetPhi                = new TH1F("NjetvsPhi", "NjetvsPhi", 288,-2*TMath::Pi(),2*TMath::Pi());
  fHistJetPt                 = new TH1F("NjetvsJetPt", "NjetvsJetPt", 300, 0, 300);
  fOutput->Add(fHistJetPhi);
  fOutput->Add(fHistJetPt);

  fillHist = 1;

  fHistClusE = new TH1F("NumberClustersvsEnergy","NumberClustersvsEnergy", 500, 0, 10);

    
  if(fillHist>0){
  fHistRhovsCent              = new TH2F("RhovsCent", "RhovsCent", 100, 0.0, 100.0, 400, 0, 400);
  fHistCorJetPt		            = new TH1F("NjetvsCorrJetPt", "NjetvsCorrJetPt", 300, -100, 200);
  fHistdEdx		                = new TH1F("dEdxSignal", "dEdxSignal", 500, 0, 500);
  fHistdEdxvPt                = new TH2F("dEdxvPt", "dEdxvPt", 200, 0, 100,500, 0 ,500);
  fHistEovPTracks             = new TH1F("EovPTracks","EovPTracks",200,0.0,2.0);
  fHistEovPvsPtTracks         = new TH2F("E/p_vs_Pt","E/p_vs_Pt",200, 0 ,100 ,50, 0, 2);
  fHistPID                    = new TH1F("fHistPID", "Detector PID", 8, 0, 8);
  fHistPIDtpc                 = new TH1F("fHistPIDtpc", "TPC pid", 4, 0, 4);
  fHistPIDits                 = new TH1F("fHistPIDits", "ITS pid", 4, 0, 4);
  fHistPIDtof                 = new TH1F("fHistPIDtof", "TOF pid", 4, 0, 4);
  fHistnsigelectron           = new TH1F("nsig_elec_TPC","nsig_elec_TPC",500,-10,10);
  fHistnSigElecPt             = new TH2F("nsig_v_pt(TPC)","nsig_v_pt(TPC)",200,0,100,100,-10,10);
  fHistTrackPhivEta           = new TH2F("TrackPhi_v_Eta","TrackPhiEta",64,-1.6,1.6,72,1,2*TMath::Pi());
  fHistClusterPhivEta         = new TH2F("ClusterPhi_v_Eta","ClusterPhiEta",64,-1.6,1.6,72,1,2*TMath::Pi());
  fHistnJetTrackvnJetClusters = new TH2F("NumbJetTracksvJetClusters","NumbJetTracksvJetClusters",21,0,20,21,0,20);

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
  fOutput->Add(fHistdEdx);
  fOutput->Add(fHistdEdxvPt);
  fOutput->Add(fHistEovPTracks);
  fOutput->Add(fHistEovPvsPtTracks);
  fOutput->Add(fHistPID);
  fOutput->Add(fHistPIDtpc);
  fOutput->Add(fHistPIDits);
  fOutput->Add(fHistPIDtof);
  fOutput->Add(fHistnsigelectron);
  fOutput->Add(fHistnSigElecPt);
  fOutput->Add(fHistTrackPhivEta);
  fOutput->Add(fHistClusterPhivEta);
  fOutput->Add(fHistnJetTrackvnJetClusters);
  
  }//Fill Histograms

    fOutput->Add(fHistClusE);


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
  // ****************************************************************************************
  UInt_t bitcoded = 0;  // bit coded, see GetDimParamsPID() below
  bitcoded = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4 | 1<<5 | 1<<6 | 1<<7 | 1<<8 | 1<<9 | 1<<10 | 1<<11;
  fhnPIDHF = NewTHnSparseDHF("fhnPIDHF", bitcoded);
  
  UInt_t bitcoded1 = 0;  // bit coded, see GetDimParamsPID() below
  bitcoded1 = 1<<0 | 1<<1 | 1<<2 | 1<<3 | 1<<4;
  fhnJetQA = NewTHnSparseDJetQA("fhnJetQA", bitcoded1);
  
  UInt_t bitcoded2 = 0;  // bit coded, see GetDimParamsPID() below
  bitcoded2 = 1<<5 | 1<<6 | 1<<7 | 1<<8 | 1<<9 | 1<<15;
  fhnClusQA = NewTHnSparseDJetQA("fhnClusQA", bitcoded2);
  
  UInt_t bitcoded3 = 0;  // bit coded, see GetDimParamsPID() below
  bitcoded3 = 1<<10 | 1<<11 | 1<<12 | 1<<13 | 1<<14;
  fhnTrackQA = NewTHnSparseDJetQA("fhnTrackQA", bitcoded3);
  
  cout << "_______________Created Sparse__________________" << endl;
  
  fOutput->Add(fhnPIDHF);
  fOutput->Add(fhnJetQA);
  fOutput->Add(fhnClusQA);
  fOutput->Add(fhnTrackQA);
  fOutput->Add(fhnGlobalPID);


  PostData(1, fOutput);
}

//________________________________________________________
void AliAnalysisTaskEmcalJetHF::ExecOnce()
{
  //  Initialize the analysis 
  AliAnalysisTaskEmcalJet::ExecOnce();

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
  TClonesArray *tracks = 0;
  TClonesArray *clusters = 0;
    
  // get Jets object                                                                                                                        
  jets = dynamic_cast<TClonesArray*>(list->FindObject(fJets));
  if(!jets){
    AliError(Form("Pointer to jets %s == 0", fJets->GetName()));
    return kTRUE;
  } // verify existence of jets                                                                                                             

  // get Tracks object                                                                                                                      
  tracks = dynamic_cast<TClonesArray*>(list->FindObject(fTracks));
  if (!tracks) {
    AliError(Form("Pointer to tracks %s == 0", fTracks->GetName()));
    return kTRUE;
  } // verify existence of tracks  (fTracksName.Data())
    
  //get Clusters object
  clusters = dynamic_cast<TClonesArray*>(list->FindObject(fCaloClusters));
  if (!clusters){
    AliError(Form("Pointer to tracks %s == 0",fCaloClusters->GetName()));
    return kTRUE;
  }  //verify cluster existence
  
  
  //const Int_t nclus = clusters->GetEntries();
  /* Globalcluster
  for(int icluster = 0; icluster < nclus; icluster++) {
    AliVCluster* andyCluster = static_cast<AliVCluster*>(fCaloClusters->At(icluster));
    if (!icluster) continue;
  }
  */
  // get all emcal clusters
  TRefArray* caloClusters = new TRefArray();
  fESD->GetEMCALClusters( caloClusters );
	
  //TObjArray* listcuts = fEsdtrackCutsTPC->GetAcceptedTracks(fESD);
  //Int_t nGoodTracks = list->GetEntries();
  Int_t nCluster = caloClusters->GetEntries();

  // event plane histograms filled
  if(fillHist>0) fHistEP0[centbin]->Fill(fEPV0);
  if(fillHist>0) fHistEP0A[centbin]->Fill(fEPV0A);
  if(fillHist>0) fHistEP0C[centbin]->Fill(fEPV0C);
  if(fillHist>0) fHistEPAvsC[centbin]->Fill(fEPV0A,fEPV0C);
  
  event++;
  cout<<"Event #: "<<event<<endl;
  
  // get number of jets and tracks                                                                                                          
  const Int_t Njets = jets->GetEntries();
  const Int_t Ntracks = tracks->GetEntries();
  //const Int_t Nclusters = clusters->GetEntries();
  if(Ntracks<1)   return kTRUE;
  if(Njets<1)     return kTRUE;
  //if(Nclusters<1)  return kTRUE;
  
  if(doGlobalPID){
  // loop over Globap tracks in the event for PID: Runs ON ESD's
  for (int i = 0;i<Ntracks;i++){
    //AliVParticle *track = static_cast<AliVParticle*>(fTracks->At(i));
    AliVTrack *trackGlobal = static_cast<AliVTrack*>(fTracks->At(i));
    //AliESDtrack* track = fESD->GetTrack(i);
    if (! trackGlobal) {
      AliError(Form("Couldn't get ESD track %d\n", i));
      continue;
    }
    

    // apply track cuts                                                                                                        
    if(TMath::Abs(trackGlobal->Eta())>fTrackEta) continue;
    if (trackGlobal->Pt()<0.15) continue;

    // initialize track info 
    Double_t dEdxGlobal = -99;
    Double_t trackphiGlobal = -99;
    Double_t trackptGlobal = 0;
    Double_t p = trackGlobal->P();
    Double_t fClsEGlobal = -99;
    Double_t EovPGlobal = -99;

    // distance of closest approach
    Float_t DCAxyGlobal = -999;
    Float_t DCAzGlobal = -999;
    //Double_t DCAXYGlobal = -999;
    //Double_t DCAZGlobal = -999;
    
    // track info    
    trackGlobal->GetTPCsignal();
    //trackGlobal->GetImpactParameters(DCAxyGlobal, DCAzGlobal);
    dEdxGlobal = trackGlobal->GetTPCsignal();
    trackptGlobal = trackGlobal->Pt();
    trackphiGlobal = trackGlobal->Phi();

    // fill track histo's
    //if(fillHist>0) fHistdEdx->Fill(dEdx);
    //if(fillHist>0) fHistdEdxvPt->Fill(trackpt,dEdx);
    //if(fillHist>0) fHistTrackPt[centbin]->Fill(track->Pt());

    // clusters                                                                                                        
    Int_t nMatchClusGlobal = -1;
    AliESDCaloCluster *matchedClusGlobal = NULL;
                                                                                                                                   
    //////////////////////////////////////////////////////////////////////////                                                                            
 
    // cut on 1.5 GeV for EMCal Cluster                                                                                                                        //if(track->GetEMCALcluster()<0 || pt<1.5) continue;                                                                                            
    //////////////////////////////////////////////////////////////////////////                                                                                   
    // particles in TOF                                                                                                                                    
    //Double_t nSigmaPion_TOFGlobal, nSigmaProton_TOFGlobal, nSigmaKaon_TOFGlobal = -1.;

    // misc quantities                                                                                                                                     
    //Double_t TOFsigGlobal  = -1.;
    //Double_t nClustersTPCGlobal = -1;
    //Int_t chargeGlobal     = 0;
    //Int_t trackCutsGlobal  = 0;     // set to 0 to get to run
    //Int_t myPIDGlboal      = 0;     // set to 0 because myPID is to compare with MC unless I hard code in cuts
    Double_t nSigmaElectron_TPCGlobal = -999;
    Double_t nSigmaElectron_TOFGlobal = -999;
    //Double_t nSigmaElectron_EMCALGlobal = -999;

    // get clusters                                                                                                                                       
    //Int_t clusiGlobal = trackGlobal->GetEMCALcluster();
    //nClustersTPCGlobal = trackGlobal->GetTPCclusters(0);
    
    AliESDtrack *trackESD = fESD->GetTrack(Ntracks);
 
    nSigmaElectron_TPCGlobal = fPIDResponse->NumberOfSigmasTPC(trackESD,AliPID::kElectron);
    nSigmaElectron_TOFGlobal = fPIDResponse->NumberOfSigmasTOF(trackESD,AliPID::kElectron);
    
    //for EMCAL
    nMatchClusGlobal = trackGlobal->GetEMCALcluster();
    if(nMatchClusGlobal > 0){
      matchedClusGlobal = (AliESDCaloCluster*)fESD->GetCaloCluster(nMatchClusGlobal);
      //AliESDCaloCluster* matchedClus = fESD->GetCaloCluster(clusi);                                                                                                                                                                                    
      //double eop2Global = -1;
      //double ssGl[4]={0.,0.,0.,0.};
      //Double_t nSigmaEop = fPID->GetPIDResponse()-m >NumberOfSigmasEMCAL(track,AliPID::kElectron,eop2,ss);

      fClsEGlobal       = matchedClusGlobal->E();
      EovPGlobal        = fClsEGlobal/p;
      //nSigmaElectron_EMCALGlobal = fPIDResponse->NumberOfSigmasEMCAL(trackESD,AliPID::kElectron,eop2,ss);
      //if(fillHist>0) fHistEovPTracks->Fill(EovP);
      //if(fillHist>0) fHistEovPvsPtTracks->Fill(trackpt,EovP);
    }

    //if(fillHist>0) fHistnsigelectron->Fill(nSigmaElectron_TPC);
    //if(fillHist>0) fHistnSigElecPt->Fill(trackpt,nSigmaElectron_TPC);

    // extra attempt                                                                                                                                      
    AliVEvent *eventQA=InputEvent();
    if (!eventQA||!fPIDResponse) return kTRUE; // just return, maybe put at beginning
    //////////////////////////////////////////////////////////////////////////                                                                             
        // cut on 1.5 GeV for EMCal Cluster
    if(trackGlobal->Pt() < fTrackPtCut) continue;

    Double_t HF_tracks[12] = {fCent, trackGlobal->Pt(), trackGlobal->P() ,trackGlobal->Eta(), trackGlobal->Phi(), EovPGlobal, DCAxyGlobal, DCAzGlobal, dEdxGlobal, nSigmaElectron_TPCGlobal, nSigmaElectron_TOFGlobal, 0};//nSigmaElectron_EMCAL};
    fhnGlobalPID->Fill(HF_tracks);    // fill Sparse Histo with trigger entries

  } // Loop over tracks
  }//PID Switch
  
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
 
        Double_t JetQA[5] = {0, static_cast<Double_t>(jet->GetNumberOfTracks()), static_cast<Double_t>(jet->GetNumberOfClusters()),jet->Eta(), jet->Phi()};
        fhnJetQA->Fill(JetQA);
      
      
        // Loop over clusters for JetQA
        for(int iCluster = 0; iCluster < nCluster; iCluster++) {
          AliVCluster* caloCluster = static_cast<AliVCluster*>(fCaloClusters->At(iCluster));
          //AliVCluster* caloCluster = (AliVCluster* )caloClusters->At(jet->GetNumberOfClusters());
          //AliESDCaloCluster* caloCluster = (AliESDCaloCluster* )caloClusters->At(iCluster);
          //AliESDCaloCluster* clus = fESD->GetCaloCluster(iclus);
          if (!caloCluster){
            AliError(Form("ERROR: Could not get cluster %d", iCluster));
            continue;
          }
          if(!IsJetCluster(jet,iCluster,kFALSE)) continue ;
          Float_t pos[3];
          caloCluster->GetPosition(pos);  // Get cluster position
          TVector3 cp(pos);
          Double_t NtrMatched = -999.0;
          NtrMatched = caloCluster->GetNTracksMatched();
          
          //loop over tracks for Jet QA
          for(int itrack = 0; itrack < NtrMatched; itrack++){
            //AliVParticle *track = static_cast<AliVParticle*>(fTracks->At(i));
            AliVTrack *trackcluster = static_cast<AliVTrack*>(fTracks->At(itrack));
            //AliESDtrack* track = fESD->GetTrack(i);
            if (! trackcluster) {
              AliError(Form("Couldn't get ESD track %d\n", itrack));
              continue;
            }
            if(!IsJetTrack(jet,itrack,kFALSE)) continue;
            Double_t ClusQA[6] = {caloCluster->E(),cp.PseudoRapidity() ,cp.Phi(), static_cast<Double_t>(caloCluster->IsEMCAL()), NtrMatched, trackcluster->Phi()};
            fhnClusQA->Fill(ClusQA); //,1./Njets);    // fill Sparse Histo with trigger entries
          }//loop over tracks for JetQA

        }//loop over clusters for JetQA
      
      // loop over tracks in the event
      for (int iTrack = 0; iTrack<jet->GetNumberOfTracks(); iTrack++){  //loop over tracks in jet for TrackQA
        //AliVParticle *track = static_cast<AliVParticle*>(fTracks->At(i));
        AliVTrack *jetTrack = static_cast<AliVTrack*>(fTracks->At(iTrack));
        //AliESDtrack* track = fESD->GetTrack(i);
        if (! jetTrack) {
          AliError(Form("Couldn't get ESD track %d\n", iTrack));
          continue;
        }
        if(!IsJetTrack(jet,iTrack,kFALSE)) continue;
        Double_t trackQA[5] = {jetTrack->Pt(), jetTrack->Eta(), jetTrack->Phi(), static_cast<Double_t>(jetTrack->IsEMCAL()), static_cast<Double_t>(jetTrack->GetEMCALcluster())};
        fhnTrackQA->Fill(trackQA); //,1./Njets);
        
      }//track loop for TrackQA

        if(highestjetpt<jet->Pt()){
            ijethi=ijet;
            highestjetpt=jet->Pt();
        }
    } // end of looping over jets

  // loop over jets in the event and make appropriate cuts
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {
     AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(iJets));
     if (!jet)  // see if we have a jet
       continue;
      
     // phi of jet, constrained to 1.6 < Phi < 2.94
     float jetphi = jet->Phi();      // phi of jet 
  
    // apply jet cuts
    if(!AcceptMyJet(jet)) continue;
 
    Int_t JetClusters = jet->GetNumberOfClusters();
    Int_t JetTracks = jet -> GetNumberOfTracks();
    fHistnJetTrackvnJetClusters->Fill(JetClusters,JetTracks);
   
     // Initializations and Calculations
     //Double_t jeteta = jet->Eta();    			        // ETA of jet
     Double_t jetptraw = jet->Pt();    				// raw pT of jet
     Double_t jetPt = -500;                                     // initialize corr jet pt LOCAL
     Double_t jetarea = -500;					// initialize jet area
     jetarea = jet->Area();		           		// jet area
     jetPt = jet->Pt() - jetarea*fRhoVal;                 // semi-corrected pT of jet from GLOBAL rho value

     // make histo's
     if(fillHist>0) fHistJetPtvsTrackPt[centbin]->Fill(jetPt,jet->MaxTrackPt());
     fHistJetPhi->Fill(jetphi);
     if(fillHist>0) fHistCorJetPt->Fill(jetPt);
     fHistJetPt->Fill(jetptraw);
    
      if(jet->Pt() > fJetHIpt) {
        
      //loop over clusters
      //for (int i = 0; i < njetclusters; i++){
      for(int iCluster = 0; iCluster < nCluster; iCluster++) {
          AliVCluster* caloCluster = static_cast<AliVCluster*>(fCaloClusters->At(iCluster));
          //AliVCluster* caloCluster = (AliVCluster* )caloClusters->At(iCluster);
          //AliESDCaloCluster* caloCluster = (AliESDCaloCluster* )caloClusters->At(iCluster);
          //AliESDCaloCluster* clus = fESD->GetCaloCluster(iclus);
          if (!caloCluster){
              AliError(Form("ERROR: Could not get cluster %d", iCluster));
              continue;
          }
        
          //if (!caloCluster -> IsEMCAL()) continue; //Check that Cluster is EMCal Cluster
          if(! IsJetCluster(jet, iCluster, kFALSE)) continue;
        
          AliESDtrack *track = 0;
          if (caloCluster->GetTrackMatchedIndex() > 0) // tender's matching
          track = fESD->GetTrack(caloCluster->GetTrackMatchedIndex());
          Double_t fclusE = -999;
          fclusE = caloCluster->E();
        
          if (fclusE<0.) continue;  //Check that cluster has positive energy
        
          fHistClusE->Fill(fclusE);
          //Int_t  cluslabel = caloCluster->GetID();
          //Int_t nmatched = caloCluster->GetNTracksMatched();
        
          Float_t pos[3];
          caloCluster->GetPosition(pos);  // Get cluster position
          TVector3 cp(pos);
          if(fillHist>0) fHistClusterPhivEta->Fill(cp.PseudoRapidity(),cp.Phi());
        
          Int_t trackMatchedindex=caloCluster->GetTrackMatchedIndex();
          if(trackMatchedindex<0)continue;  // Make sure we don't have a bad index
 
          if (caloCluster->GetTrackMatchedIndex() > 0) // tender's matching
          track = fESD->GetTrack(caloCluster->GetTrackMatchedIndex());
          Double_t NtrMatched = -999;
          NtrMatched = caloCluster->GetNTracksMatched();
          for(int itrack = 0; itrack < NtrMatched; itrack++){   // Loop over tracks matched to clusters from jets
            //AliVParticle *trackCluster = static_cast<AliVParticle*>(fTracks->At(i));
            AliVTrack *trackCluster = static_cast<AliVTrack*>(fTracks->At(itrack));
            //AliESDtrack* trackCluster = fESD->GetTrack(itrack);
            if (! trackCluster) {
              AliError(Form("Couldn't get ESD track %d\n", itrack));
              continue;
            }
            if(!IsJetTrack(jet,itrack,kFALSE)) continue;  // Check that track is still part of ith jet
            
              //if (trackCluster->Phi()>3.2 || trackCluster->Phi()<1.4)cout <<"Out of Range Track!    Track Phi:  " << trackCluster->Phi() << endl;
  
              //if(track->Pt() < 4.0) continue;
              // initialize track info
              Double_t dEdx = -99;
              //Double_t trackphi = -99;
              Double_t trackpt = 0;
              Double_t p = trackCluster->P();
              Double_t EovP = -99;
              
              // distance of closest approach
              Float_t DCAxy = -999;
              Float_t DCAz = -999;
              //Double_t DCAXY = -999;
              //Double_t DCAZ = -999;
              
              // track info
              //trackCluster->GetImpactParameters(DCAxy, DCAz);
              trackpt = trackCluster->Pt();
              //trackphi = track->Phi();
              
              // fill track histo's
              if(fillHist>0) fHistTrackPhivEta->Fill(trackCluster->Eta(),trackCluster->Phi());
              if(fillHist>0) fHistdEdx->Fill(dEdx);
              if(fillHist>0) fHistdEdxvPt->Fill(trackpt,dEdx);
              if(fillHist>0) fHistTrackPt[centbin]->Fill(trackCluster->Pt());
            
              Double_t nSigmaElectron_TPC = -999;
              Double_t nSigmaElectron_TOF = -999;
              Double_t nSigmaElectron_EMCAL = -999;
  
            if(!useAOD){
              
              AliESDtrack *trackESD = fESD->GetTrack(Ntracks);
              
              dEdx = trackESD->GetTPCsignal();
              trackESD->GetImpactParameters(DCAxy, DCAz);
              nSigmaElectron_TPC = fPIDResponse->NumberOfSigmasTPC(trackESD,AliPID::kElectron);
              nSigmaElectron_TOF = fPIDResponse->NumberOfSigmasTOF(trackESD,AliPID::kElectron);
            }
           if (useAOD) {
              AliAODTrack *trackAOD = fAOD->GetTrack(Ntracks);
              
              // get detector signals
              dEdx = trackAOD->GetTPCsignal();
             
             //trackAOD->GetImpactParameters(DCAxy,DCAz);
              nSigmaElectron_TPC = fPIDResponse->NumberOfSigmasTPC(trackAOD,AliPID::kElectron);
              nSigmaElectron_TOF = fPIDResponse->NumberOfSigmasTOF(trackAOD,AliPID::kElectron);
              
            } // end of AOD pid
      
              //double eop2 = -1;
              //double ss[4]={0.,0.,0.,0.};
              EovP        = fclusE/p;
              //nSigmaElectron_EMCAL = fPIDResponse->NumberOfSigmasEMCAL(trackESD,AliPID::kElectron,eop2,ss);
              if(fillHist>0) fHistEovPTracks->Fill(EovP);
              if(fillHist>0) fHistEovPvsPtTracks->Fill(trackpt,EovP);
              Double_t HF_tracks[12] = {fCent, trackCluster->Pt(), trackCluster->P() ,trackCluster->Eta(), trackCluster->Phi(), EovP, 0/*DCAxy*/, 0/*DCAz*/, dEdx, nSigmaElectron_TPC, nSigmaElectron_TOF, nSigmaElectron_EMCAL};
              fhnPIDHF->Fill(HF_tracks);    // fill Sparse Histo with trigger entries
          
          //Int_t trackMatchedIndex = caloCluster->GetTrackMatchedIndex();//find the index of the matched track. This by default returns the best match
          //AliESDtrack *track = event->GetTrack(trackMatchedIndex);
          //if this is a good track, accept track will return true. The track matched is good, so not track matched is false
          //Int_t matched = fESD->AcceptTrack(track);//If the track is bad, don't count it.  By default even bad tracks are accepted

        }  //loop over tracks matched to clusters in jet
        
      } // loop over jet clusters
    
          
      } // highest pt jet cut
  } // LOOP over JETS in event

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
    nbins = 750;
    xmin = 0.;
    xmax = 75.;
    break;

  case 2:
    label = "Track p";
    nbins = 750;
    xmin = 0.;
    xmax = 75.;
    break;

  case 3:
    label = "Track Eta";
    nbins = 64;
    xmin = -1.6;
    xmax = 1.6;
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


  } // end of switch                                                                                                                                                     
} // end of getting dim-params

void AliAnalysisTaskEmcalJetHF::GetDimParamsJetQA(Int_t iEntry, TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax)
{
  // stores label and binning of axis for THnSparse
  const Double_t pi = TMath::Pi();
  
  switch(iEntry){
      
    case 0:
      label = "number of Jets in Event";
      nbins = 2;
      xmin = 0.;
      xmax = 1.;
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
      nbins = 1000;
      xmin = 0;
      xmax = 10;
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
      nbins = 750;
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



