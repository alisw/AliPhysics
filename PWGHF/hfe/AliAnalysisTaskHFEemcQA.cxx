#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "THnSparse.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"

#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"

#include "AliAnalysisTaskHFEemcQA.h"

using std::cout;
using std::endl;
//QA task for EMCAL electron analysis 

ClassImp(AliAnalysisTaskHFEemcQA)
  //________________________________________________________________________
  AliAnalysisTaskHFEemcQA::AliAnalysisTaskHFEemcQA(const char *name) 
: AliAnalysisTaskSE(name),
  fVevent(0),
  fESD(0),
  fAOD(0),
  fpidResponse(0),
  fFlagSparse(kFALSE),
  fOutputList(0),
  fNevents(0),
  fVtxZ(0),
  fVtxX(0),
  fVtxY(0),
  fTrigMulti(0),
  fHistClustE(0),
  fEMCClsEtaPhi(0),
  fNegTrkIDPt(0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fdEdx(0),
  fTPCNpts(0),
  fTPCnsig(0),
  fHistPtMatch(0),
  fEMCTrkMatch(0),
  fEMCTrkPt(0),
  fEMCTrketa(0),
  fEMCTrkphi(0),
  fEMCdEdx(0),
  fEMCTPCnsig(0),
  fEMCTPCNpts(0),
  fHistdEdxEop(0),
  fHistNsigEop(0),
  fHistEop(0),
  fEleCanTPCNpts(0),
  fEleCanTPCNCls(0),
  fEleCanITSNCls(0),
  fEleCanITShit(0),
  fEleCanSPD1(0),
  fEleCanSPD2(0),
  fEleCanSPDBoth(0),
  fEleCanSPDOr(0),
  fSparseElectron(0),
  fvalueElectron(0)
{
  // Constructor

  fvalueElectron = new Double_t[6];
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskHFEemcQA::AliAnalysisTaskHFEemcQA() 
  : AliAnalysisTaskSE("DefaultTask_HfeEMCQA"),
  fVevent(0),
  fESD(0),
  fAOD(0),
  fpidResponse(0),
  fFlagSparse(kFALSE),
  fOutputList(0),
  fNevents(0),
  fVtxZ(0),
  fVtxX(0),
  fVtxY(0),
  fTrigMulti(0),
  fHistClustE(0),
  fEMCClsEtaPhi(0),
  fNegTrkIDPt(0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fdEdx(0),
  fTPCNpts(0),
  fTPCnsig(0),
  fHistPtMatch(0),
  fEMCTrkMatch(0),
  fEMCTrkPt(0),
  fEMCTrketa(0),
  fEMCTrkphi(0),
  fEMCdEdx(0),
  fEMCTPCnsig(0),
  fEMCTPCNpts(0),
  fHistdEdxEop(0),
  fHistNsigEop(0),
  fHistEop(0),
  fEleCanTPCNpts(0),
  fEleCanTPCNCls(0),
  fEleCanITSNCls(0),
  fEleCanITShit(0),
  fEleCanSPD1(0),
  fEleCanSPD2(0),
  fEleCanSPDBoth(0),
  fEleCanSPDOr(0),
  fSparseElectron(0),
  fvalueElectron(0)
{
  //Default constructor

  fvalueElectron = new Double_t[6];
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //DefineOutput(3, TTree::Class());
}
//________________________________________________________________________
AliAnalysisTaskHFEemcQA::~AliAnalysisTaskHFEemcQA()
{
  //Destructor 
  delete fOutputList;
  delete fSparseElectron;
  delete []fvalueElectron;
}
//________________________________________________________________________
void AliAnalysisTaskHFEemcQA::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  AliDebug(3, "Creating Output Objects");

  /////////////////////////////////////////////////
  //Automatic determination of the analysis mode//
  ////////////////////////////////////////////////
  AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
    SetAODAnalysis();
  } else {
    SetESDAnalysis();
  }
  printf("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD");

  //////////////// 
  //Output list//
  ///////////////
  fOutputList = new TList();
  fOutputList->SetOwner();  

  fNevents = new TH1F("fNevents","No of events",3,-0.5,2.5);
  fOutputList->Add(fNevents);
  fNevents->GetYaxis()->SetTitle("counts");
  fNevents->GetXaxis()->SetBinLabel(1,"All");
  fNevents->GetXaxis()->SetBinLabel(2,"With >2 Trks");
  fNevents->GetXaxis()->SetBinLabel(3,"Vtx_{z}<10cm");

  fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
  fOutputList->Add(fVtxZ);

  fVtxY = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",1000,-50,50);
  fOutputList->Add(fVtxY);

  fVtxX = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",1000,-50,50);
  fOutputList->Add(fVtxX);

  fTrigMulti = new TH2F("fTrigMulti","Multiplicity distribution for different triggers; Trigger type; multiplicity",11,-1,10,2000,0,2000);
  fOutputList->Add(fTrigMulti);

  fHistClustE = new TH1F("fHistClustE", "EMCAL cluster energy distribution; Cluster E;counts", 500, 0.0, 50.0);
  fOutputList->Add(fHistClustE);

  fEMCClsEtaPhi = new TH2F("fEMCClsEtaPhi","EMCAL cluster #eta and #phi distribution;#eta;#phi",100,-0.9,0.9,200,0,6.3);
  fOutputList->Add(fEMCClsEtaPhi);

  fNegTrkIDPt = new TH1F("fNegTrkIDPt", "p_{T} distribution of tracks with negative track id;p_{T} (GeV/c);counts", 500, 0.0, 50.0); 
  fOutputList->Add(fNegTrkIDPt);

  fTrkPt = new TH1F("fTrkPt","p_{T} distribution of all tracks;p_{T} (GeV/c);counts",1000,0,100);
  fOutputList->Add(fTrkPt);

  fTrketa = new TH1F("fTrketa","All Track #eta distribution;#eta;counts",100,-1.5,1.5);
  fOutputList->Add(fTrketa);

  fTrkphi = new TH1F("fTrkphi","All Track #phi distribution;#phi;counts",100,0,6.3);
  fOutputList->Add(fTrkphi);

  fdEdx = new TH2F("fdEdx","All Track dE/dx distribution;p (GeV/c);dE/dx",200,0,20,500,0,160);
  fOutputList->Add(fdEdx);

  fTPCNpts = new TH2F("fTPCNpts","All track TPC Npoints used for dE/dx calculation;p (GeV/c);N points",200,0,20,200,0.,200.);
  fOutputList->Add(fTPCNpts);

  fTPCnsig = new TH2F("fTPCnsig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
  fOutputList->Add(fTPCnsig);

  fHistPtMatch = new TH1F("fHistPtMatch", "p_{T} distribution of tracks matched to EMCAL;p_{T} (GeV/c);counts",1000, 0.0, 100.0);
  fOutputList->Add(fHistPtMatch);                                      

  fEMCTrkMatch = new TH2F("fEMCTrkMatch","Distance of EMCAL cluster to its closest track;#phi;z",100,-0.3,0.3,100,-0.3,0.3);
  fOutputList->Add(fEMCTrkMatch);

  fEMCTrkPt = new TH1F("fEMCTrkPt","p_{T} distribution of tracks with EMCAL cluster;p_{T} (GeV/c);counts",1000,0,100);
  fOutputList->Add(fEMCTrkPt);

  fEMCTrketa = new TH1F("fEMCTrketa","#eta distribution of tracks matched to EMCAL;#eta;counts",100,-1.5,1.5);
  fOutputList->Add(fEMCTrketa);

  fEMCTrkphi = new TH1F("fEMCTrkphi","#phi distribution of tracks matched to EMCAL;#phi;counts",100,0,6.3);
  fOutputList->Add(fEMCTrkphi);

  fEMCdEdx = new TH2F("fEMCdEdx","dE/dx distribution of tracks matched to EMCAL;p (GeV/c);dE/dx",200,0,20,500,0,160);
  fOutputList->Add(fEMCdEdx);
    
  fEMCTPCnsig = new TH2F("fEMCTPCnsig","TPC Nsigma distribution of tracks matched to EMCAL;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
  fOutputList->Add(fEMCTPCnsig);

  fEMCTPCNpts = new TH2F("fEMCTPCNpts","TPC Npoints used for dE/dx for tracks matched to EMCAL;p (GeV/c);N points",200,0,20,200,0.,200.);
  fOutputList->Add(fEMCTPCNpts);

  fHistEop = new TH2F("fHistEop", "E/p distribution;p_{T} (GeV/c);E/p", 200,0,20,60, 0.0, 3.0);
  fOutputList->Add(fHistEop);

  fHistdEdxEop = new TH2F("fHistdEdxEop", "E/p vs dE/dx;E/p;dE/dx", 60, 0.0, 3.0, 500,0,160);
  fOutputList->Add(fHistdEdxEop);
    
  fHistNsigEop = new TH2F ("fHistNsigEop", "E/p vs TPC nsig",60, 0.0, 3.0, 200, -10,10);
  fOutputList->Add(fHistNsigEop);

  fEleCanTPCNpts = new TH2F("fEleCanTPCNpts","TPC Npoints used for dE/dx for electron candidates;p_{T} (GeV/c);N points",200,0,20,200,0,200);
  fOutputList->Add(fEleCanTPCNpts);

  fEleCanTPCNCls = new TH2F("fEleCanTPCNCls","TPC N clusters for electron candidates;p_{T} (GeV/c);N TPC clusters",200,0,20,171,-0.5,170.5);
  fOutputList->Add(fEleCanTPCNCls);

  fEleCanITSNCls = new TH2F("fEleCanITSNCls","ITS N clusters for electron candidates;p_{T} (GeV/c);N ITS clusters",200,0,20,8,-0.5,7.5);
  fOutputList->Add(fEleCanITSNCls);

  fEleCanITShit = new TH1F("fEleCanITShit","ITS hit map;ITS layer;counts",7,-0.5,6.5);
  fOutputList->Add(fEleCanITShit);

  fEleCanSPD1 = new TH2F("fEleCanSPD1","Hit on SPD layer 1;p_{T} (GeV/c);Hit",200,0,20,1,0,1);
  fOutputList->Add(fEleCanSPD1);

  fEleCanSPD2 = new TH2F("fEleCanSPD2","Hit on SPD layer 2;p_{T} (GeV/c);Hit",200,0,20,1,0,1);
  fOutputList->Add(fEleCanSPD2);

  fEleCanSPDBoth = new TH2F("fEleCanSPDBoth","Tracks with hits on both SPD layer;p_{T} (GeV/c);Hit",200,0,20,1,0,1);
  fOutputList->Add(fEleCanSPDBoth);

  fEleCanSPDOr = new TH2F("fEleCanSPDOr","Tracks with hits on both SPD layer;p_{T} (GeV/c);Hit",200,0,20,1,0,1);
  fOutputList->Add(fEleCanSPDOr);

  Int_t bins[6]={8,500,200,400,400,400}; //trigger, pt, TPCnsig, E/p, M20, M02
  Double_t xmin[6]={-0.5,0,-10,0,0,0};
  Double_t xmax[6]={7.5,25,10,2,2,2};
  fSparseElectron = new THnSparseD ("Electron","Electron",6,bins,xmin,xmax);
  fOutputList->Add(fSparseElectron);

  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHFEemcQA::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  // Post output data.

  UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fVevent) {
    printf("ERROR: fVEvent not available\n");
    return;
  }

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (fESD) {
    //   printf("fESD available\n");
    //return;
  }

  ////////////////////
  //cuts initialised//
  ///////////////////
  AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
  esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
  esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
  esdTrackCutsH->SetDCAToVertex2D(kTRUE);

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (fAOD) {
    // printf("fAOD available\n");
    //return;
  }

  ///////////////////
  //PID initialised//
  //////////////////
  fpidResponse = fInputHandler->GetPIDResponse();

  ////////////////
  //Event vertex//
  ///////////////
  Int_t ntracks;
  ntracks = fVevent->GetNumberOfTracks();
  printf("There are %d tracks in this event\n",ntracks);

  fNevents->Fill(0); //all events
  Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
  Double_t NcontV = pVtx->GetNContributors();
  if(NcontV<2)return;
  fNevents->Fill(1); //events with 2 tracks

  Zvertex = pVtx->GetZ();  
  Yvertex = pVtx->GetY();  
  Xvertex = pVtx->GetX();  
  fVtxZ->Fill(Zvertex);
  fVtxX->Fill(Xvertex);
  fVtxY->Fill(Yvertex);

  /////////////////
  //trigger check//
  /////////////////
  TString firedTrigger;
  TString TriggerEG1("EG1");
  TString TriggerEG2("EG2");
  fVevent->GetFiredTriggerClasses();

  Bool_t EG1tr = kFALSE;
  Bool_t EG2tr = kFALSE;

  if(firedTrigger.Contains(TriggerEG1))EG1tr = kTRUE;
  if(firedTrigger.Contains(TriggerEG2))EG2tr = kTRUE;

  Int_t trigger = -1;
  if (fAOD){
    Double_t multiplicity=fAOD->GetHeader()->GetRefMultiplicity();
    fTrigMulti->Fill(-0.5, multiplicity);
    if(evSelMask & AliVEvent::kAny) fTrigMulti->Fill(0.5, multiplicity);
    if(evSelMask & AliVEvent::kMB) fTrigMulti->Fill(1.5, multiplicity);
    if(evSelMask & AliVEvent::kINT7) fTrigMulti->Fill(2.5, multiplicity);
    if(evSelMask & AliVEvent::kINT8) fTrigMulti->Fill(3.5, multiplicity);
    if(evSelMask & AliVEvent::kEMC1) fTrigMulti->Fill(4.5, multiplicity);
    if(evSelMask & AliVEvent::kEMC7) fTrigMulti->Fill(5.5, multiplicity);
    if(evSelMask & AliVEvent::kEMC8) fTrigMulti->Fill(6.5, multiplicity);
    if(evSelMask & AliVEvent::kEMCEJE) fTrigMulti->Fill(7.5, multiplicity);
    if(evSelMask & AliVEvent::kEMCEGA) fTrigMulti->Fill(8.5, multiplicity);
    if(evSelMask & AliVEvent::kEMCEGA & EG2tr) fTrigMulti->Fill(9.5, multiplicity);

    if(evSelMask & AliVEvent::kMB) trigger =0;
    if(evSelMask & AliVEvent::kINT7) trigger =1;
    if(evSelMask & AliVEvent::kINT8) trigger =2;
    if(evSelMask & AliVEvent::kEMC1) trigger =3;
    if(evSelMask & AliVEvent::kEMC7) trigger =4;
    if(evSelMask & AliVEvent::kEMC8) trigger =5;
    if(evSelMask & AliVEvent::kEMCEJE) trigger =6;
    if(evSelMask & AliVEvent::kEMCEGA) trigger =7;
  }

  ////////////////////
  //event selection//
  ///////////////////
  if(fabs(Zvertex>10.0))return; 
  fNevents->Fill(2); //events after z vtx cut

  /////////////////////////////
  //EMCAL cluster information//
  ////////////////////////////
  Int_t Nclust = 0;
  Nclust = fVevent->GetNumberOfCaloClusters();
  for(Int_t icl=0; icl<Nclust; icl++)
  {
    AliVCluster *clust = 0x0;
    clust = fVevent->GetCaloCluster(icl);
    if(clust && clust->IsEMCAL())
    {
      Double_t clustE = clust->E();
      Float_t  emcx[3]; // cluster pos
      clust->GetPosition(emcx);
      TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
      Double_t emcphi = clustpos.Phi(); 
      Double_t emceta = clustpos.Eta();
      fHistClustE->Fill(clustE);
      fEMCClsEtaPhi->Fill(emceta,emcphi);
    }
  }

  /////////////////////////////////
  //Look for kink mother for AOD//
  /////////////////////////////////
  Int_t numberofvertices = 100;
  if(fAOD) numberofvertices = fAOD->GetNumberOfVertices();
  Double_t listofmotherkink[numberofvertices];
  Int_t numberofmotherkink = 0;
  if(IsAODanalysis())
  {
    for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
      AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
      if(!aodvertex) continue;
      if(aodvertex->GetType()==AliAODVertex::kKink) {
        AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
        if(!mother) continue;
        Int_t idmother = mother->GetID();
        listofmotherkink[numberofmotherkink] = idmother;
        numberofmotherkink++;
      }
    }
  } //+++


  ///////////////
  //Track loop///
  ///////////////
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {

    AliVParticle* Vtrack = fVevent->GetTrack(iTracks);
    if (!Vtrack) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
    AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
    AliESDtrack *etrack = dynamic_cast<AliESDtrack*>(Vtrack);
    AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);

    ////////////////////
    //Apply track cuts//
    ////////////////////
    if(fAOD)
      if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts

    if(fESD)
      if(!esdTrackCutsH->AcceptTrack(etrack))continue;

    //reject kink
    if(IsAODanalysis()){
      Bool_t kinkmotherpass = kTRUE;
      for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
        if(track->GetID() == listofmotherkink[kinkmother]) {
          kinkmotherpass = kFALSE;
          continue;
        }
      }
      if(!kinkmotherpass) continue;
    }
    else{
      if(etrack->GetKinkIndex(0) != 0) continue;
    }


    ////////////////////
    //Track properties//
    ///////////////////
    Double_t dEdx =-999, fTPCnSigma=-999;
    dEdx = track->GetTPCsignal();
    fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);

    if(track->GetID()<0) fNegTrkIDPt->Fill(track->Pt());
    fTrkPt->Fill(track->Pt());
    fTrketa->Fill(track->Eta());
    fTrkphi->Fill(track->Phi());
    fdEdx->Fill(track->P(),dEdx);
    fTPCNpts->Fill(track->P(),track->GetTPCsignalN());
    fTPCnsig->Fill(track->P(),fTPCnSigma);

    ///////////////////////////
    //Track matching to EMCAL//
    //////////////////////////
    Int_t EMCalIndex = -1;
    EMCalIndex = track->GetEMCALcluster();
    if(EMCalIndex < 0) continue;
    fHistPtMatch->Fill(track->Pt());

    AliVCluster *clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex);
    if(clustMatch && clustMatch->IsEMCAL())
    {
      /////////////////////////////////////////////
      //Properties of tracks matched to the EMCAL//
      /////////////////////////////////////////////
      fEMCTrkMatch->Fill(clustMatch->GetTrackDx(),clustMatch->GetTrackDz());
      fEMCTrkPt->Fill(track->Pt());
      fEMCTrketa->Fill(track->Eta());
      fEMCTrkphi->Fill(track->Phi());
      fEMCdEdx->Fill(track->P(),dEdx);
      fEMCTPCnsig->Fill(track->P(),fTPCnSigma);
      fEMCTPCNpts->Fill(track->P(),track->GetTPCsignalN());

      //E/p distribution
      Double_t clustMatchE = clustMatch->E();
      Double_t eop = -1.0;
      if(track->P()>0)eop = clustMatchE/track->P();

      if(track->Pt()>1.0){
        fHistdEdxEop->Fill(eop,dEdx);
        fHistNsigEop->Fill(eop,fTPCnSigma);
      }
      fHistEop->Fill(track->Pt(),eop);

      //EID THnsparse
      fvalueElectron[0] = trigger;
      fvalueElectron[1] = track->Pt();
      fvalueElectron[2] = fTPCnSigma;
      fvalueElectron[3] = eop;
      fvalueElectron[4] = clustMatch->GetM20();
      fvalueElectron[5] = clustMatch->GetM02();

      if(fFlagSparse){
        cout << "filling sparse"<<endl;
        fSparseElectron->Fill(fvalueElectron);
      }

      ////////////////////////////////////////////////
      //Track properties of EMCAL electron cadidates//
      ///////////////////////////////////////////////
      if(fTPCnSigma > -1 && fTPCnSigma < 3 && eop>0.8 && eop<1.2){
        fEleCanTPCNpts->Fill(track->Pt(),track->GetTPCsignalN());
        fEleCanTPCNCls->Fill(track->Pt(),track->GetTPCNcls());

        Int_t fITSncls=0;
        for(Int_t l=0;l<6;l++) {
          if(TESTBIT(track->GetITSClusterMap(),l)) {
            fEleCanITShit->Fill(l);
            if(l==0) fEleCanSPD1->Fill(track->Pt(),0.5);
            if(l==1) fEleCanSPD2->Fill(track->Pt(),0.5);
            if(l==0 && l==1) fEleCanSPDBoth->Fill(track->Pt(),0.5);
            if(l==0 || l==1) fEleCanSPDOr->Fill(track->Pt(),0.5);
            fITSncls++;
          }
        }
        fEleCanITSNCls->Fill(track->Pt(),fITSncls++);
      }
    }
  } //track loop 

  PostData(1, fOutputList);
}      
//________________________________________________________________________
void AliAnalysisTaskHFEemcQA::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }

}
