/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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


//////////////////////////////////////////////
//  B->e & D->e  at high pT with EMC in Run2   
//  Based on AOD
//  Author: Shingo Sakai      //
//////////////////////////////////////////////

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

#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"

#include "AliMultSelection.h"
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

#include "AliCentrality.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliEMCALTriggerPatchInfo.h"

#include "AliAnalysisTaskBeautyCal.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskBeautyCal)
  //________________________________________________________________________
  AliAnalysisTaskBeautyCal::AliAnalysisTaskBeautyCal(const char *name)
: AliAnalysisTaskSE(name),
  fVevent(0),
  fESD(0),
  fAOD(0),
  fpidResponse(0),
  fCFM(0),
  fFlagSparse(kFALSE),
  fUseTender(kTRUE),
  fEMCEG1(kFALSE),
  fEMCEG2(kFALSE),
  fTracks_tender(0),
  fCaloClusters_tender(0),
  fMCparticle(0),
  fMCarray(0),
  fMultSelection(0),
  fTriggersInfo(0),
  fThresholdEG2(89),
  fThresholdEG1(140),
  fFlagClsTypeEMC(kTRUE),
  fFlagClsTypeDCAL(kTRUE),
  fOutputList(0),
  fNevents(0),
  fCent(0),
  fVtxZ(0),
  fHistClustE(0),
  fHistClustEcent(0),
  fEMCClsEtaPhi(0),
  fHistClustEEG1(0),
  fHistClustEEG1cent(0),
  fHistClustEEG2(0),
  fHistClustEEG2cent(0),
  fEMCClsEtaPhiEG1(0),
  fEMCClsEtaPhiEG2(0),
  fNegTrkIDPt(0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fdEdx(0),
  fTPCNpts(0),
  fTPCnsig(0),
  fHistPtMatch(0),
  fClsEAftMatch(0),
  fClsEtaPhiAftMatch(0),
  fClsEtaPhiAftMatchEMCin(0),
  fClsEtaPhiAftMatchEMCout(0),
  fHistdEdxEop(0),
  fHistNsigEop(0),
  fHistEop(0),
  fM20(0),
  fM02(0),
  fM20EovP(0),
  fM02EovP(0),
  fInvmassULS(0),
  fInvmassLS(0),
  fMCcheckMother(0),
  fSparseElectron(0),
  fvalueElectron(0),
  fHistDCAinc(0),
  fHistDCApho(0),
  fHistDCAhfe(0),
  fHistDCAde(0),
  fHistDCAbe(0),
  fHistDCApe(0),
  fHistHFEcorr(0),
  fhfeCuts(0) 
{
  // Constructor

  fvalueElectron = new Double_t[9];
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskBeautyCal::AliAnalysisTaskBeautyCal()
  : AliAnalysisTaskSE("DefaultTask_HfeEMCQA"),
  fVevent(0),
  fESD(0),
  fAOD(0),
  fpidResponse(0),
  fCFM(0),
  fFlagSparse(kFALSE),
  fUseTender(kTRUE),
  fEMCEG1(kFALSE),
  fEMCEG2(kFALSE),
  fTracks_tender(0),
  fCaloClusters_tender(0),
  fMCparticle(0),
  fMCarray(0),
  fMultSelection(0),
  fTriggersInfo(0),
  fThresholdEG2(89),
  fThresholdEG1(140),
  fFlagClsTypeEMC(kTRUE),
  fFlagClsTypeDCAL(kTRUE),
  fOutputList(0),
  fNevents(0),
  fCent(0), 
  fVtxZ(0),
  fHistClustE(0),
  fHistClustEcent(0),
  fEMCClsEtaPhi(0),
  fHistClustEEG1(0),
  fHistClustEEG1cent(0),
  fHistClustEEG2(0),
  fHistClustEEG2cent(0),
  fEMCClsEtaPhiEG1(0),
  fEMCClsEtaPhiEG2(0),
  fNegTrkIDPt(0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fdEdx(0),
  fTPCNpts(0),
  fTPCnsig(0),
  fHistPtMatch(0),
  fClsEAftMatch(0),
  fClsEtaPhiAftMatch(0),
  fClsEtaPhiAftMatchEMCin(0),
  fClsEtaPhiAftMatchEMCout(0),
  fHistdEdxEop(0),
  fHistNsigEop(0),
  fHistEop(0),
  fM20(0),
  fM02(0),
  fM20EovP(0),
  fM02EovP(0),
  fInvmassULS(0),
  fInvmassLS(0),
  fMCcheckMother(0), 
  fSparseElectron(0),
  fvalueElectron(0),
  fHistDCAinc(0),
  fHistDCApho(0),
  fHistDCAhfe(0),
  fHistDCAde(0),
  fHistDCAbe(0),
  fHistDCApe(0),
  fHistHFEcorr(0),
  fhfeCuts(0) 
{
  //Default constructor

  fvalueElectron = new Double_t[10];
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
AliAnalysisTaskBeautyCal::~AliAnalysisTaskBeautyCal()
{
  //Destructor
  delete fOutputList;
  delete fTracks_tender;
  delete fCaloClusters_tender;
  delete fSparseElectron;
  delete []fvalueElectron;
}
//________________________________________________________________________
void AliAnalysisTaskBeautyCal::UserCreateOutputObjects()
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

  /////////////////
  // cut from HFE package
  /////////////////////
	
 fCFM = new AliCFManager;
 const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack;
 fCFM->SetNStepParticle(kNcutSteps);
 for(Int_t istep = 0; istep < kNcutSteps; istep++) fCFM->SetParticleCutsList(istep, NULL);

  if(!fhfeCuts)
	{
	 fhfeCuts = new AliHFEcuts;
	 fhfeCuts->CreateStandardCuts();

	 fhfeCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);	
	 fhfeCuts->SetMinNClustersTPC(100);							                //Minimum number of clusters on TPC
	 fhfeCuts->SetMinNClustersTPCPID(80);										//Minimum number of clusters for dE/dx
	 fhfeCuts->SetMinRatioTPCclusters(0.6);						                    //Number of clusters (Found/Findable)
	 fhfeCuts->SetCutITSpixel(AliHFEextraCuts::kAny);							//Require at least one cluster on SPD
	 fhfeCuts->SetCheckITSLayerStatus(kFALSE); 
	 fhfeCuts->SetMinNClustersITS(3);								            //Minimum number of clusters on ITS
	 fhfeCuts->SetMaxImpactParam(2,3); //changed z to 3
	 fhfeCuts->SetVertexRange(10.);													//
	}
  
  fhfeCuts->Initialize(fCFM);


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

  fCent = new TH1F("fCent","Centrality",100,0,100);
  fOutputList->Add(fCent);

  fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
  fOutputList->Add(fVtxZ);

  fHistClustE = new TH1F("fHistClustE", "EMCAL cluster energy distribution; Cluster E;counts", 5000, 0.0, 50.0);
  fOutputList->Add(fHistClustE);

  fHistClustEcent = new TH2F("fHistClustEcent", "EMCAL cluster energy distribution vs. centrality; Cluster E;counts", 100,0,100,5000, 0.0, 50.0);
  fOutputList->Add(fHistClustEcent);

  fEMCClsEtaPhi = new TH2F("fEMCClsEtaPhi","EMCAL cluster #eta and #phi distribution;#eta;#phi",100,-0.9,0.9,200,0,6.3);
  fOutputList->Add(fEMCClsEtaPhi);

  fHistClustEEG1 = new TH1F("fHistClustEEG1", "EMCAL cluster energy distribution; Cluster E;counts", 5000, 0.0, 50.0);
  fOutputList->Add(fHistClustEEG1);

  fHistClustEEG1cent = new TH2F("fHistClustEEG1cent", "EMCAL cluster energy distribution vs. centrality; Cluster E;counts", 100,0,100,5000, 0.0, 50.0);
  fOutputList->Add(fHistClustEEG1cent);

  fHistClustEEG2 = new TH1F("fHistClustEEG2", "EMCAL cluster energy distribution; Cluster E;counts", 5000, 0.0, 50.0);
  fOutputList->Add(fHistClustEEG2);

  fHistClustEEG2cent = new TH2F("fHistClustEEG2cent", "EMCAL cluster energy distribution vs. centrality; Cluster E;counts", 100,0,100,5000, 0.0, 50.0);
  fOutputList->Add(fHistClustEEG2cent);

  fEMCClsEtaPhiEG1 = new TH2F("fEMCClsEtaPhiEG1","EMCAL cluster #eta and #phi distribution;#eta;#phi",100,-0.9,0.9,200,0,6.3);
  fOutputList->Add(fEMCClsEtaPhiEG1);

  fEMCClsEtaPhiEG2 = new TH2F("fEMCClsEtaPhiEG2","EMCAL cluster #eta and #phi distribution;#eta;#phi",100,-0.9,0.9,200,0,6.3);
  fOutputList->Add(fEMCClsEtaPhiEG2);
  //

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

  fClsEtaPhiAftMatch = new TH2F("fClsEtaPhiAftMatch","EMCAL cluster #eta and #phi distribution after track matching;#eta;#phi",100,-0.9,0.9,200,0,6.3);
  fOutputList->Add(fClsEtaPhiAftMatch);

  fClsEtaPhiAftMatchEMCin = new TH2F("fClsEtaPhiAftMatchEMCin","EMCAL cluster #eta and #phi distribution after track matching inside EMC #phi acceptence;#eta;#phi",100,-0.9,0.9,200,0,6.3);
  fOutputList->Add(fClsEtaPhiAftMatchEMCin);

  fClsEtaPhiAftMatchEMCout = new TH2F("fClsEtaPhiAftMatchEMCout","EMCAL cluster #eta and #phi distribution after track matching outside EMC #phi acceptence;#eta;#phi",100,-0.9,0.9,200,0,6.3);
  fOutputList->Add(fClsEtaPhiAftMatchEMCout);


  fHistEop = new TH2F("fHistEop", "E/p distribution;p_{T} (GeV/c);E/p", 200,0,20,60, 0.0, 3.0);
  fOutputList->Add(fHistEop);

  fHistdEdxEop = new TH2F("fHistdEdxEop", "E/p vs dE/dx;E/p;dE/dx", 60, 0.0, 3.0, 500,0,160);
  fOutputList->Add(fHistdEdxEop);

  fHistNsigEop = new TH2F ("fHistNsigEop", "E/p vs TPC nsig",60, 0.0, 3.0, 200, -10,10);
  fOutputList->Add(fHistNsigEop);

  fM20 = new TH2F ("fM20","M20 vs pt distribution",200,0,20,400,0,2);
  fOutputList->Add(fM20);

  fM02 = new TH2F ("fM02","M02 vs pt distribution",200,0,20,400,0,2);
  fOutputList->Add(fM02);

  fM20EovP = new TH2F ("fM20EovP","M20 vs E/p distribution",400,0,3,400,0,2);
  fOutputList->Add(fM20EovP);

  fM02EovP = new TH2F ("fM02EovP","M02 vs E/p distribution",400,0,3,400,0,2);
  fOutputList->Add(fM02EovP);

  fInvmassLS = new TH1F("fInvmassLS", "Invmass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 200,0,0.4);
  fOutputList->Add(fInvmassLS);

  fInvmassULS = new TH1F("fInvmassULS", "Invmass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 200,0,0.4);
  fOutputList->Add(fInvmassULS);

  fMCcheckMother = new TH1F("fMCcheckMother", "Mother MC PDG", 1000,-0.5,999.5);
  fOutputList->Add(fMCcheckMother);

  //Int_t bins[8]={8,500,200,400,400,400,400,400}; //trigger, pt, TPCnsig, E/p, M20, M02, sqrt(M20),sqrt(M02)
  //Double_t xmin[8]={-0.5,0,-10,0,0,0,0,0};
  //Double_t xmax[8]={7.5,25,10,2,2,2,2,2};
  //fSparseElectron = new THnSparseD ("Electron","Electron;trigger;pT;nSigma;eop;m20;m02;sqrtm20;sqrtm02;",8,bins,xmin,xmax);
  //Int_t bins[10]={8,500,200,400,400,400,400,3,200,10}; //trigger, pt, TPCnsig, E/p, M20, M02, sqrt(M20),sqrt(M02)
  //Double_t xmin[10]={-0.5,   0,  -8,   0,   0,   0,    0, -0.5, -8,   0};
  Int_t bins[10]=      {8, 280, 160, 200, 200, 200,  200,    3, 100,  10}; //trigger, pt, TPCnsig, E/p, M20, M02, sqrt(M20),sqrt(M02)
  Double_t xmin[10]={-0.5,   2,  -8,   0,   0,   0,    0, -0.5,  -5,   0};
  Double_t xmax[10]={ 7.5,  30,   8,   2,   2,   2,    2,  2.5,  15 , 100};
  fSparseElectron = new THnSparseD ("Electron","Electron;trigger;pT;nSigma;eop;m20;m02;sqrtm02m20;eID;nSigma_Pi;cent;",10,bins,xmin,xmax);
  fOutputList->Add(fSparseElectron);

  fHistDCAinc = new TH2D("fHistDCAinc", "DCA of inclusive e; p_{T}(GeV/c);DCAxchargexMag.", 30,0,30,10000,-0.2,0.2);
  fOutputList->Add(fHistDCAinc);
 
  fHistDCApho = new TH2D("fHistDCApho", "DCA of pho e; p_{T}(GeV/c);DCAxchargexMag.", 30,0,30,4000,-0.2,0.2);
  fOutputList->Add(fHistDCApho);

  fHistDCAhfe = new TH2D("fHistDCAhfe", "DCA of hfe e; p_{T}(GeV/c);DCAxchargexMag.", 30,0,30,4000,-0.2,0.2);
  fOutputList->Add(fHistDCAhfe);

  fHistDCAde = new TH2D("fHistDCAde", "DCA of D-> e; p_{T}(GeV/c);DCAxchargexMag.", 30,0,30,10000,-0.2,0.2);
  fOutputList->Add(fHistDCAde);
 
  fHistDCAbe = new TH2D("fHistDCAbe", "DCA of B-> e; p_{T}(GeV/c);DCAxchargexMag.", 30,0,30,4000,-0.2,0.2);
  fOutputList->Add(fHistDCAbe);

  fHistDCApe = new TH2D("fHistDCApe", "DCA of pi0/eta-> e; p_{T}(GeV/c);DCAxchargexMag.", 30,0,30,4000,-0.2,0.2);
  fOutputList->Add(fHistDCApe);

  fHistHFEcorr = new TH1D("fHistHFEcorr", "HFE corr", 720,-3.6,3.6);
  fOutputList->Add(fHistHFEcorr);

  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskBeautyCal::UserExec(Option_t *)
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

  //////////////
  //if Tender //
  //////////////
  if(fUseTender){
    //new branches with calibrated tracks and clusters
    if(IsAODanalysis()) fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
    if(!IsAODanalysis()) fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("ESDFilterTracks"));
    fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
  }

  ////////////////////
  //cuts initialised//
  ////////////////////
  AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
  esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
  esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
  esdTrackCutsH->SetDCAToVertex2D(kTRUE);
  esdTrackCutsH->SetMinNClustersTPC(80);
  esdTrackCutsH->SetMinNClustersITS(3);
  esdTrackCutsH->SetRequireTPCRefit(kTRUE);
  esdTrackCutsH->SetRequireITSRefit(kTRUE);
  esdTrackCutsH->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  esdTrackCutsH->SetAcceptKinkDaughters(kFALSE);
  esdTrackCutsH->SetMaxChi2PerClusterITS(6); //test.....

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (fAOD) {
    // printf("fAOD available\n");
    //return;
  }
  if(fAOD)fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));

  ///////////////////
  //PID initialised//
  ///////////////////
  fpidResponse = fInputHandler->GetPIDResponse();

  ///////////////////
  // centrality
 /////////////////////
  
  Double_t centrality = -1;
  AliCentrality *fCentrality = (AliCentrality*)fAOD->GetCentrality(); 
  //centrality = fCentrality->GetCentralityPercentile("V0M");

  //Double_t centrality = -1;
  if(fAOD)fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
  if( !fMultSelection) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    //AliWarning("AliMultSelection object not found!");
    centrality = fCentrality->GetCentralityPercentile("V0M");
  }else{
   //lPercentile = fMultSelection->GetMultiplicityPercentile("V0M");
   centrality = fMultSelection->GetMultiplicityPercentile("V0M", false); 
 }
  //printf("mim cent selection %d\n",fcentMim);
  //printf("max cent selection %d\n",fcentMax);
  //printf("cent selection %d\n",centrality);

  if(fcentMim>-0.5)
    {
     if(centrality < fcentMim || centrality > fcentMax)return;
    }

  ////////////////
  //Event vertex//
  ////////////////
  Int_t ntracks = -999;
  if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries();
  //if(ntracks < 1) printf("There are %d tracks in this event\n",ntracks);
  //printf("There are %d tracks in this event\n",ntracks);

  fNevents->Fill(0); //all events
  Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
  Double_t NcontV = pVtx->GetNContributors();
  if(NcontV<2)return;
  fNevents->Fill(1); //events with 2 tracks

  Zvertex = pVtx->GetZ();
  fVtxZ->Fill(Zvertex);

  /////////////////
  //trigger check//
  /////////////////
  TString firedTrigger;
  TString TriggerEG1("EG1");
  TString TriggerEG2("EG2");
  fVevent->GetFiredTriggerClasses();
  if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();

  Bool_t EG1tr = kFALSE;
  Bool_t EG2tr = kFALSE;

  if(firedTrigger.Contains(TriggerEG1))EG1tr = kTRUE;
  if(firedTrigger.Contains(TriggerEG2))EG2tr = kTRUE;

  if(fEMCEG1){if(!firedTrigger.Contains(TriggerEG1))return;}
  if(fEMCEG2){if(!firedTrigger.Contains(TriggerEG2))return;}

  Int_t trigger = -1;
  if (fAOD){
    //Double_t multiplicity=fAOD->GetHeader()->GetRefMultiplicity();
    AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
    if(!header) AliFatal("Not a standard AOD");
    Double_t multiplicity = header->GetRefMultiplicity();

    if(evSelMask & AliVEvent::kMB) trigger =0;
    if(evSelMask & AliVEvent::kINT7) trigger =1;
    if(evSelMask & AliVEvent::kINT8) trigger =2;
    if(evSelMask & AliVEvent::kEMC1) trigger =3;
    if(evSelMask & AliVEvent::kEMC7) trigger =4;
    if(evSelMask & AliVEvent::kEMC8) trigger =5;
    if(evSelMask & AliVEvent::kEMCEJE) trigger =6;
    if(evSelMask & AliVEvent::kEMCEGA) trigger =7;
  }


  ////////////////////////
  // Mag. field
  ///////////////////////
 
  Int_t MagSign = 1;
  if(fVevent->GetMagneticField()<0)MagSign = -1;

  ////////////////////
  //event selection///
  ////////////////////
  if(TMath::Abs(Zvertex)>10.0)return;
  fNevents->Fill(2); //events after z vtx cut
  fCent->Fill(centrality); //centrality dist.

  /////////////////////////////
  //EMCAL cluster information//
  /////////////////////////////
  Int_t Nclust = -999;
  if(!fUseTender) Nclust = fVevent->GetNumberOfCaloClusters();
  if(fUseTender) Nclust = fCaloClusters_tender->GetEntries();

  Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;

  for(Int_t icl=0; icl<Nclust; icl++)
  {
    AliVCluster *clust = 0x0;
    if(!fUseTender) clust = fVevent->GetCaloCluster(icl);
    if(fUseTender) clust = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(icl));
    if(!clust)  printf("ERROR: Could not receive cluster matched calibrated from track %d\n", icl);

    fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;

    if(clust && clust->IsEMCAL())
    {
      /////////////////////////////////
      //Select EMCAL or DCAL clusters//
      /////////////////////////////////
      Float_t  emcx[3]; // cluster pos
      clust->GetPosition(emcx);
      TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
      Double_t emcphi = clustpos.Phi();
      Double_t emceta = clustpos.Eta();
      if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.

      if(emcphi > 1.39 && emcphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
      if(emcphi > 4.53 && emcphi < 5.708) fClsTypeDCAL = kTRUE; //DCAL  : 260 < phi < 327

      //----selects EMCAL+DCAL clusters when fFlagClsTypeEMC and fFlagClsTypeDCAL is kTRUE
      if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
        if(!fClsTypeEMC) continue; //selecting only EMCAL clusters

      if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
        if(!fClsTypeDCAL) continue; //selecting only DCAL clusters

      Double_t clustE = clust->E();
      fHistClustE->Fill(clustE);
      if(centrality>-1)fHistClustEcent->Fill(centrality,clustE);
      fEMCClsEtaPhi->Fill(emceta,emcphi);

      //-----Plots for EMC trigger
      Bool_t hasfiredEG1=0;
      Bool_t hasfiredEG2=0;
      FindPatches(hasfiredEG1,hasfiredEG2,emceta,emcphi);
      if(hasfiredEG1){
        fHistClustEEG1->Fill(clustE);
        if(centrality>-1)fHistClustEEG1cent->Fill(centrality,clustE);
        fEMCClsEtaPhiEG1->Fill(emceta,emcphi);
      }
      if(hasfiredEG2){
        fHistClustEEG2->Fill(clustE);
        if(centrality>-1)fHistClustEEG2cent->Fill(centrality,clustE);
        fEMCClsEtaPhiEG2->Fill(emceta,emcphi);
      }

    }
  }

  // cell information
  AliVCaloCells *fCaloCells = fVevent->GetEMCALCells();

  ////////////////////////////////
  //Look for kink mother for AOD//
  ////////////////////////////////
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

    AliVParticle* Vtrack = 0x0;
    if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
    if(fUseTender) Vtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));

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

    //other cuts
    Double_t d0z0[2]={-999,-999}, cov[3];
    Double_t DCAxyCut = 2.4, DCAzCut = 3.2;

    if(fAOD){
      //cout << "AOD track cuts" << endl;
      if(atrack->GetTPCNcls() < 80) continue;
      if(atrack->GetITSNcls() < 3) continue;
      if((!(atrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
      if(!(atrack->HasPointOnITSLayer(0) || atrack->HasPointOnITSLayer(1))) continue;

      if(atrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
        if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;
      //To be done : Add cuts to apply Chi2PerITSCls < 6 and N shared Cls ITS < 4
      if(atrack->P()<2.0)continue;
    }
 
    /* 
     if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
     if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;
     if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
     if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
    */

    Double_t DCAxy = d0z0[0]*atrack->Charge()*MagSign;
    //cout << "DCAxy = " << DCAxy << endl;

    ///////////////////////
    // Get MC information//
    ///////////////////////
    Int_t ilabel = track->GetLabel();
    Int_t pdg = -999;
    Double_t pid_ele = 0.0;
    Int_t pidM = -1;
    if(ilabel>0 && fMCarray)
    {
      fMCparticle = (AliAODMCParticle*) fMCarray->At(ilabel);
      Int_t pdg = fMCparticle->GetPdgCode();
      if(TMath::Abs(pdg)==11)pid_ele = 1.0;
      Int_t ilabelM = -1;
      if(pid_ele==1.0)FindMother(fMCparticle, ilabelM, pidM);
      if(pidM==22) // from pi0 & eta
        {
          AliAODMCParticle* fMCparticleM = (AliAODMCParticle*) fMCarray->At(ilabelM);
          FindMother(fMCparticleM, ilabelM, pidM);
        }
      fMCcheckMother->Fill(abs(pidM));
    }

    if(abs(pdg)==11)cout << " pid_ele = " << pid_ele << " ; pidM = " << pidM << endl;

    Bool_t pid_eleD = IsDdecay(pidM);
    Bool_t pid_eleB = IsBdecay(pidM);
    Bool_t pid_eleP = IsPdecay(pidM);

    ////////////////////
    //Track properties//
    ///////////////////
    Double_t dEdx =-999, fTPCnSigma=-999, fTPCnSigma_Pi=-999;
    Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999;
    dEdx = track->GetTPCsignal();
    fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    TrkPhi = track->Phi();
    TrkPt = track->Pt();
    TrkEta = track->Eta();
    TrkP = track->P();

    if(track->GetID()<0) fNegTrkIDPt->Fill(track->Pt());
    fTrkPt->Fill(TrkPt);
    fTrketa->Fill(TrkEta);
    fTrkphi->Fill(TrkPhi);
    fdEdx->Fill(TrkP,dEdx);
    fTPCNpts->Fill(TrkP,track->GetTPCsignalN());
    fTPCnsig->Fill(TrkP,fTPCnSigma);

    ///////////////////////////
    //Track matching to EMCAL//
    //////////////////////////
    Int_t EMCalIndex = -1;
    EMCalIndex = track->GetEMCALcluster();
    if(EMCalIndex < 0) continue;
    fHistPtMatch->Fill(track->Pt());

    AliVCluster *clustMatch=0x0;
    if(!fUseTender) clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex);
    if(fUseTender) clustMatch = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(EMCalIndex));

    Double_t emcphi = -999, emceta=-999;
    fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;
    if(clustMatch && clustMatch->IsEMCAL())
    {
      if(TMath::Abs(clustMatch->GetTrackDx())>0.05 || TMath::Abs(clustMatch->GetTrackDz())>0.05) continue;

      /////////////////////////////////
      //Select EMCAL or DCAL clusters//
      /////////////////////////////////
      Float_t  emcx[3]; // cluster pos
      clustMatch->GetPosition(emcx);
      TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
      emcphi = clustpos.Phi();
      emceta = clustpos.Eta();
      if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
      if(emcphi > 1.39 && emcphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
      if(emcphi > 4.53 && emcphi < 5.708) fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327

      //----selects EMCAL+DCAL clusters when fFlagClsTypeEMC and fFlagClsTypeDCAL is kTRUE
      if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
        if(!fClsTypeEMC) continue; //selecting only EMCAL clusters

      if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
        if(!fClsTypeDCAL) continue; //selecting only DCAL clusters

      Double_t clustMatchE = clustMatch->E();
      //fClsEAftMatch->Fill(clustMatchE);

      fClsEtaPhiAftMatch->Fill(emceta,emcphi);

      if(TrkPhi > 1.396  && TrkPhi < 3.141) //emc acceptance (80 to 180 degrees)
      {
        fClsEtaPhiAftMatchEMCin->Fill(emceta,emcphi);
      }
      else{
        fClsEtaPhiAftMatchEMCout->Fill(emceta,emcphi);
      }

      Float_t tof = clustMatch->GetTOF()*1e+9; // ns

      //EMCAL EID info
      Double_t eop = -1.0;
      Double_t m02 = -99999,m20 = -99999,sqm02m20=-99999.0;
      if(track->P()>0)eop = clustMatchE/track->P();
      m02 =clustMatch->GetM02();
      m20 =clustMatch->GetM20();
      sqm02m20 = sqrt(pow(m02,2)+pow(m20,2));

      if(track->Pt()>2.0){
        fHistdEdxEop->Fill(eop,dEdx);
        fHistNsigEop->Fill(eop,fTPCnSigma);
        fM20EovP->Fill(eop,clustMatch->GetM20());
        fM02EovP->Fill(eop,clustMatch->GetM02());
      }
      fHistEop->Fill(track->Pt(),eop);
      fM20->Fill(track->Pt(),clustMatch->GetM20());
      fM02->Fill(track->Pt(),clustMatch->GetM02());

      //EID THnsparse
      fvalueElectron[0] = trigger;
      fvalueElectron[1] = track->Pt();
      fvalueElectron[2] = fTPCnSigma;
      fvalueElectron[3] = eop;
      fvalueElectron[4] = clustMatch->GetM20();
      fvalueElectron[5] = clustMatch->GetM02();
      fvalueElectron[6] = sqm02m20;
      fvalueElectron[7] = pid_ele;
      fvalueElectron[8] = fTPCnSigma_Pi;
      fvalueElectron[9] = centrality;

      if(fFlagSparse && track->Pt()>2.0){
        fSparseElectron->Fill(fvalueElectron);
      }

      Bool_t fFlagNonHFE=kFALSE;
      ////////////////////////////////////////////////
      //Track properties of EMCAL electron cadidates//
      ////////////////////////////////////////////////
      if(fTPCnSigma > -1 && fTPCnSigma < 3 && eop>0.9 && eop<1.3 && m02 > 0.006 && m02 < 0.35){ //rough cuts
        //-----Identify Non-HFE
        SelectPhotonicElectron(iTracks,track,fFlagNonHFE);

        fHistDCAinc->Fill(track->Pt(),DCAxy);
        if(fFlagNonHFE)
          { 
           fHistDCApho->Fill(track->Pt(),DCAxy);
          }
        else
         {
          fHistDCAhfe->Fill(track->Pt(),DCAxy);
          ElectronAway(iTracks,track); //e+e-
         }

        if(pid_eleD)fHistDCAde->Fill(track->Pt(),DCAxy);
        if(pid_eleB)fHistDCAbe->Fill(track->Pt(),DCAxy);
        if(pid_eleP)fHistDCApe->Fill(track->Pt(),DCAxy);
   
      }
    }
  } //track loop

  PostData(1, fOutputList);
}

//_________________________________________________________________
Bool_t AliAnalysisTaskBeautyCal::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
		//Check single track cuts for a given cut step
		//Note this function is called inside the UserExec function
	const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
	if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
	return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskBeautyCal::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec)
{
  ///////////////////////////////////////////
  //////Non-HFE - Invariant mass method//////
  ///////////////////////////////////////////

  AliESDtrackCuts* esdTrackCutsAsso = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  esdTrackCutsAsso->SetAcceptKinkDaughters(kFALSE);
  esdTrackCutsAsso->SetRequireTPCRefit(kTRUE);
  esdTrackCutsAsso->SetRequireITSRefit(kTRUE);
  esdTrackCutsAsso->SetEtaRange(-0.9,0.9);
  esdTrackCutsAsso->SetMaxChi2PerClusterTPC(4);
  esdTrackCutsAsso->SetMinNClustersTPC(70);
  esdTrackCutsAsso->SetMaxDCAToVertexZ(3.2);
  esdTrackCutsAsso->SetMaxDCAToVertexXY(2.4);
  esdTrackCutsAsso->SetDCAToVertex2D(kTRUE);

  Bool_t flagPhotonicElec = kFALSE;

  Int_t ntracks = -999;
  if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries();

  for (Int_t jtrack = 0; jtrack < ntracks; jtrack++) {
    AliVParticle* VAssotrack = 0x0;
    if(!fUseTender) VAssotrack  = fVevent->GetTrack(jtrack);
    if(fUseTender) VAssotrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jtrack)); //take tracks from Tender list

    if (!VAssotrack) {
      printf("ERROR: Could not receive track %d\n", jtrack);
      continue;
    }

    AliVTrack *Assotrack = dynamic_cast<AliVTrack*>(VAssotrack);
    AliESDtrack *eAssotrack = dynamic_cast<AliESDtrack*>(VAssotrack);
    AliAODTrack *aAssotrack = dynamic_cast<AliAODTrack*>(VAssotrack);

    //------reject same track
    if(jtrack==itrack) continue;

    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    Double_t ptAsso=-999., nsigma=-999.0, mass=-999., width = -999;
    Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;

    nsigma = fpidResponse->NumberOfSigmasTPC(Assotrack, AliPID::kElectron);
    ptAsso = Assotrack->Pt();
    Int_t chargeAsso = Assotrack->Charge();
    Int_t charge = track->Charge();
    if(charge>0) fPDGe1 = -11;
    if(chargeAsso>0) fPDGe2 = -11;
    if(charge == chargeAsso) fFlagLS = kTRUE;
    if(charge != chargeAsso) fFlagULS = kTRUE;

    //------track cuts applied
    if(fAOD) {
      if(!aAssotrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
      if(aAssotrack->GetTPCNcls() < 70) continue;
      if((!(aAssotrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aAssotrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
    }
    else{
      if(!esdTrackCutsAsso->AcceptTrack(eAssotrack)) continue;
    }

    //-------loose cut on partner electron
    if(ptAsso <0.2) continue;
    if(aAssotrack->Eta()<-0.9 || aAssotrack->Eta()>0.9) continue;
    if(nsigma < -3 || nsigma > 3) continue;

    //-------define KFParticle to get mass
    AliKFParticle::SetField(fVevent->GetMagneticField());
    AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
    AliKFParticle ge2 = AliKFParticle(*Assotrack, fPDGe2);
    AliKFParticle recg(ge1, ge2);

    if(recg.GetNDF()<1) continue;
    Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
    if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;

    //-------Get mass
    Int_t MassCorrect;
    MassCorrect = recg.GetMass(mass,width);

    if(fFlagLS)
      if(track->Pt()>1) fInvmassLS->Fill(mass);
    if(fFlagULS)
      if(track->Pt()>1) fInvmassULS->Fill(mass);

    if(mass<0.1 && fFlagULS && !flagPhotonicElec)
      flagPhotonicElec = kTRUE; //Tag Non-HFE (random mass cut, not optimised)
  }
  fFlagPhotonicElec = flagPhotonicElec;
}


//________________________________________________________________________
void AliAnalysisTaskBeautyCal::ElectronAway(Int_t itrack, AliVTrack *track)
{
  ///////////////////////////////////////////
  //////Non-HFE - Invariant mass method//////
  ///////////////////////////////////////////

  Int_t ntracks = -999;
  if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries();

  for (Int_t jtrack = 0; jtrack < ntracks; jtrack++) {
    AliVParticle* VAssotrack = 0x0;
    if(!fUseTender) VAssotrack  = fVevent->GetTrack(jtrack);
    if(fUseTender) VAssotrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jtrack)); //take tracks from Tender list

    if (!VAssotrack) {
      printf("ERROR: Could not receive track %d\n", jtrack);
      continue;
    }

    AliVTrack *Assotrack = dynamic_cast<AliVTrack*>(VAssotrack);
    AliAODTrack *aAssotrack = dynamic_cast<AliAODTrack*>(VAssotrack);

    //------reject same track
    if(jtrack==itrack) continue;

    Double_t ptAsso=-999., nsigma=-999.0;

    nsigma = fpidResponse->NumberOfSigmasTPC(Assotrack, AliPID::kElectron);
    ptAsso = Assotrack->Pt();

    //------track cuts applied

      if(aAssotrack->GetTPCNcls() < 80) continue;
      if(aAssotrack->GetITSNcls() < 3) continue;
      if((!(aAssotrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aAssotrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
      if(!(aAssotrack->HasPointOnITSLayer(0) || aAssotrack->HasPointOnITSLayer(1))) continue;

    //-------loose cut on partner electron
    if(ptAsso <0.2) continue;
    if(aAssotrack->Eta()<-0.6 || aAssotrack->Eta()>0.6) continue;
    if(nsigma < 0 || nsigma > 3) continue;
    Double_t dphie = aAssotrack->Phi() - track->Phi();
    if(track->Pt()>2.0)fHistHFEcorr->Fill(dphie);
  }
}


void AliAnalysisTaskBeautyCal::FindMother(AliAODMCParticle* part, int &label, int &pid)
{

 if(part->GetMother()>-1)
   {
    label = part->GetMother();
    AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(label);
    pid = partM->GetPdgCode();
   }
   //cout << "Find Mother : label = " << label << " ; pid" << pid << endl;
}

Bool_t AliAnalysisTaskBeautyCal::IsDdecay(int mpid)
{
 int abmpid = fabs(mpid);
 if(abmpid==411 || abmpid==421 || abmpid==413 || abmpid==423 || abmpid==431 || abmpid==433)
   {
    return kTRUE;
   }
 else
   {
    return kFALSE;
   } 
}

Bool_t AliAnalysisTaskBeautyCal::IsBdecay(int mpid)
{
 int abmpid = fabs(mpid);
 if(abmpid==511 || abmpid==521 || abmpid==513 || abmpid==523 || abmpid==431 || abmpid==533)
   {
    return kTRUE;
   }
 else
   {
    return kFALSE;
   } 
}

Bool_t AliAnalysisTaskBeautyCal::IsPdecay(int mpid)
{
 int abmpid = fabs(mpid);
 if(abmpid==22 || abmpid==111 || abmpid==221)
   {
    return kTRUE;
   }
 else
   {
    return kFALSE;
   } 
}
//________________________________________________________________________
void AliAnalysisTaskBeautyCal::FindPatches(Bool_t &hasfiredEG1,Bool_t &hasfiredEG2,Double_t emceta,Double_t emcphi)
{
  //Find trigger patches

  fTriggersInfo = dynamic_cast <TClonesArray*>(InputEvent()->FindListObject("EmcalTriggers"));
  if(!fTriggersInfo) return;
  Int_t nPatch = fTriggersInfo->GetEntries();;
  AliEMCALTriggerPatchInfo* patch=0;
  for( int iPatch = 0; iPatch < nPatch; iPatch++ ){
    patch = (AliEMCALTriggerPatchInfo*)fTriggersInfo->At( iPatch );
    if(patch->GetADCAmp()<fThresholdEG2) continue;
    if(patch->GetEtaMin()>emceta) continue;
    if(patch->GetEtaMax()<emceta) continue;
    if(patch->GetPhiMin()>emcphi) continue;
    if(patch->GetPhiMax()<emcphi) continue;
    if(patch->GetADCAmp()>fThresholdEG2)  hasfiredEG2=1;
    if(patch->GetADCAmp()>fThresholdEG1)  hasfiredEG1=1;
  }
}

//________________________________________________________________________
void AliAnalysisTaskBeautyCal::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}
