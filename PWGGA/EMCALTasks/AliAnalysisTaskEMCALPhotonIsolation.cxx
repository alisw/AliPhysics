// $Id$
//
// Emcal Neutral Cluster analysis base task.
//
// Authors: D.Lodato,L.Ronflette, M.Marquard



#include <TClonesArray.h>
#include <TChain.h>
#include <TList.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THnSparse.h>
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliClusterContainer.h"
#include "AliVTrack.h"
#include "AliEmcalParticle.h"
#include "AliParticleContainer.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVCaloCells.h"
#include "AliPicoTrack.h"

#include "AliAnalysisTaskEMCALPhotonIsolation.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEMCALPhotonIsolation)

//________________________________________________________________________
AliAnalysisTaskEMCALPhotonIsolation::AliAnalysisTaskEMCALPhotonIsolation() :
AliAnalysisTaskEmcal("AliAnalysisTaskEMCALPhotonIsolation",kTRUE),
//fParticleCollArray(),
fNCluster(0),
fWho(-1),
//fOutputList(0),
fTrackMult(0),
fTrackMultEMCAL(0),
fClustMult(0),
fPVZBefore(0),
fEtaPhiCell(0),
fEtaPhiClus(0),
fClusEvsClusT(0),
fGoodEventsOnPVZ(0),
fPT(0),
fM02(0),
fNLM(0),
fDeltaETAClusTrackVSpT(0),
fDeltaPHIClusTrackVSpT(0),
fEtIsoCells(0),
fEtIsoClust(0),
fPtIsoTrack(0),
fPtEtIsoTC(0),
fPhiBandUEClust(0),
fEtaBandUEClust(0),
fPhiBandUECells(0),
fEtaBandUECells(0),
fPhiBandUETracks(0),
fEtaBandUETracks(0),
fPerpConesUETracks(0),
fTPCWithoutIsoConeB2BbandUE(0),
fNTotClus10GeV(0),
fRecoPV(0),
fEtIsolatedCells(0),
fEtIsolatedClust(0),
fEtIsolatedTracks(0),
fTest(0),
fOutputTHnS(0),
fOutPTMAX(0),
fOutputTree(0),
fIsoConeRadius(0.4),
fEtIsoMethod(0),
fEtIsoThreshold(4),
fdetacut(0.025),
fdphicut(0.03),
fM02mincut(0.1),
fM02maxcut(0.3),
fQA(0),
fIsMC(0),
fTPC4Iso(0),
fIsoMethod(0),
fUEMethod(0),
fVertex(0),
fNDimensions(0),
fisLCAnalysis(0),
fevents(0),
flambda0T(0),
fEtT(0),
fPtT(0),
fEtisoT(0),
fPtisoT(0),
fetaT(0),
fphiT(0),
fsumEtisoconeT(0),
fsumEtUE(0)
//tracks(0),
//clusters(0)

{
  // Default constructor.
  
  //fParticleCollArray.SetOwner(kTRUE);
  // for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
  
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEMCALPhotonIsolation::AliAnalysisTaskEMCALPhotonIsolation(const char *name, Bool_t histo) :
AliAnalysisTaskEmcal(name, histo),
//fParticleCollArray(),
fNCluster(0),
fWho(-1),
//fOutputList(0),
fTrackMult(0),
fTrackMultEMCAL(0),
fClustMult(0),
fPVZBefore(0),
fEtaPhiCell(0),
fEtaPhiClus(0),
fClusEvsClusT(0),
fGoodEventsOnPVZ(0),
fPT(0),
fM02(0),
fNLM(0),
fDeltaETAClusTrackVSpT(0),
fDeltaPHIClusTrackVSpT(0),
fEtIsoCells(0),
fEtIsoClust(0),
fPtIsoTrack(0),
fPtEtIsoTC(0),
fPhiBandUEClust(0),
fEtaBandUEClust(0),
fPhiBandUECells(0),
fEtaBandUECells(0),
fPhiBandUETracks(0),
fEtaBandUETracks(0),
fPerpConesUETracks(0),
fTPCWithoutIsoConeB2BbandUE(0),
fNTotClus10GeV(0),
fRecoPV(0),
fEtIsolatedCells(0),
fEtIsolatedClust(0),
fEtIsolatedTracks(0),
fTest(0),
fOutputTHnS(0),
fOutPTMAX(0),
fOutputTree(0),
fIsoConeRadius(0.4),
fEtIsoMethod(0),
fEtIsoThreshold(4),
fdetacut(0.025),
fdphicut(0.03),
fM02mincut(0.1),
fM02maxcut(0.3),
fQA(0),
fIsMC(0),
fTPC4Iso(0),
fIsoMethod(0),
fUEMethod(0),
fVertex(0),
fNDimensions(0),
fisLCAnalysis(0),
fevents(0),
flambda0T(0),
fEtT(0),
fPtT(0),
fEtisoT(0),
fPtisoT(0),
fetaT(0),
fphiT(0),
fsumEtisoconeT(0),
fsumEtUE(0)
//tracks(0),
//clusters(0)

{
  // Standard constructor.
  
  //fParticleCollArray.SetOwner(kTRUE);
  //    for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
  
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEMCALPhotonIsolation::~AliAnalysisTaskEMCALPhotonIsolation(){
  // Destructor
}


//________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::UserCreateOutputObjects(){
  // Create ouput histograms and THnSparse and TTree
  
  AliAnalysisTaskEmcal::UserCreateOutputObjects();
  
  
  if ((fIsoMethod == 0 || fIsoMethod == 1) && fTPC4Iso){
    cout<<"Error: Iso_Methods with CELLS and CLUSTERS work only within EMCAL "<<endl;
    cout<<"Please Set Iso_Method and TPC4Iso Accordingly!!"<<endl;
    return;
  }
  if ((fIsoMethod == 0 || fIsoMethod == 1) && fUEMethod> 1) {
    cout<<"Error: UE_Methods with CELLS and CLUSTERS work only within EMCAL"<<endl;
    cout<<"Please Set Iso_Method and UE_Method Accordingly!!"<<endl;
    return;
  }
  
  
  TString sIsoMethod="\0",sUEMethod="\0",sBoundaries="\0";
  
  if(fIsoMethod==0)
    sIsoMethod = "Cells";
  else if(fIsoMethod==1)
    sIsoMethod = "Clust";
  else if(fIsoMethod==2)
    sIsoMethod = "Tracks";
  
  if(fUEMethod==0)
    sUEMethod = "PhiBand";
  else if(fUEMethod==1)
    sUEMethod = "EtaBand";
  else if(fUEMethod==2)
    sUEMethod = "PerpCones";
  else if(fUEMethod==3)
    sUEMethod = "FullTPC";
  
  if(fTPC4Iso)
    sBoundaries = "TPC Acceptance";
  else
    sBoundaries = "EMCAL Acceptance";
  
  if(fWho>1 || fWho==-1){
    cout<<"Error!!! OutputMode Can Only Be 0: TTree; 1: THnSparse"<<endl;
    return;
  }
  else{
    fOutput = new TList();
    fOutput->SetOwner();
    //Initialize the common Output histograms
    switch (fWho)
    {
      case 0:
        //Initialization by Lucile and Marco
        fOutputTree = new TTree("OutTree","OutTree");
        fOutputTree->Branch("fevents",&fevents);
        fOutputTree->Branch("flambda0T",&flambda0T);
        fOutputTree->Branch("fEtT",&fEtT);
        fOutputTree->Branch("fPtT",&fPtT);
        fOutputTree->Branch("fEtisoT",&fEtisoT);
        fOutputTree->Branch("fPtTiso",&fPtisoT);
        fOutputTree->Branch("fetaT",&fetaT);
        fOutputTree->Branch("fphiT",&fphiT);
        fOutputTree->Branch("fsumEtisoconeT",&fsumEtisoconeT);
        fOutputTree->Branch("fsumEtUE",&fsumEtUE);
        
        fOutput->Add(fOutputTree);
        
        break;
      case 1:
        //Initialization by Davide;
        
        TString sTitle;
        Int_t binTrackMult=100, binPT=70, binM02=100, binETiso=100, binETUE=110, binETisoUE=110, binetacl=140,binphicl=100, binPTMC=70;
        Int_t bins[] = {binTrackMult, binPT, binM02, binETiso, binETUE, binETisoUE, binetacl, binphicl, binPTMC, binPTMC};
        
        fNDimensions = sizeof(bins)/sizeof(Int_t);
        const Int_t ndims =   fNDimensions;
        
        Double_t xmin[]= {  0.,  0., 0., -10., -10., -10.,-1.0, 1. ,-10.,-10.};
        
        Double_t xmax[]= {1000., 70., 2., 100., 100., 100., 1.0, 3.5, 60., 60.};
        
        sTitle = Form("Direct Photons: Track Multiplicity, p_{T} , M02 , E_{T} Iso%s in %s, E_{T} UE %s in %s, E_{T} Iso_%s - E_{T} UE_%s in %s, #eta_{clus} distribution,#phi_{clus} distribution,MC_pT,MC_pT_incone; N_{ch}; p_{T} (GeV/c); M02; E_{T}^{iso%s} (GeV/c) ; E_{T}^{UE%s} (GeV/c); E_{T}^{iso%s}-E_{T}^{UE%s} (GeV/c); #eta_{cl}; #phi_{cl}; MC p_{T}; MC p_{T}^{incone}", sIsoMethod.Data(), sBoundaries.Data(), sUEMethod.Data(), sBoundaries.Data(), sIsoMethod.Data(), sUEMethod.Data(), sBoundaries.Data(), sIsoMethod.Data(), sUEMethod.Data(), sIsoMethod.Data(), sUEMethod.Data());
        
        fOutputTHnS =  new THnSparseF("fHnOutput",sTitle.Data(), ndims, bins, xmin, xmax);
        
        fOutput->Add(fOutputTHnS);
        
        Int_t binsbis[] = {binPT, binM02, binETisoUE};
        Double_t xminbis[] = {  0., 0., -10.};
        Double_t xmaxbis[] = {100., 2., 100.};
        
        fOutPTMAX = new THnSparseF ("fOutPTMAX","3D matrix E_{#gamma} VS M02 VS pT_{max}^{cone}; E_{T}^{#gamma} (GeV/c); M02; p_{T}^{Iso}(GeV/c)",3,binsbis,xminbis,xmaxbis);
        
        fOutput->Add(fOutPTMAX);
        break;
    }
  }
  
  if(fQA){
    //Include QA plots to the OutputList //DEFINE BETTER THE BINNING AND THE AXES LIMITS
    fTrackMult = new TH1D ("hTrackMult","Tracks multiplicity Distribution",250,0.,1000.);
    fTrackMult->Sumw2();
    fOutput->Add(fTrackMult);
    
    fTrackMultEMCAL = new TH1D ("hTrackMultEMCAL","Tracks multiplicity Distribution inside EMCAL acceptance",250,0.,1000.);
    fTrackMultEMCAL->Sumw2();
    fOutput->Add(fTrackMultEMCAL);
    
    fClustMult = new TH1D ("hClustMult","Clusters multiplicity Distribution",250,0.,1000.);
    fClustMult->Sumw2();
    fOutput->Add(fClustMult);
    
    fRecoPV = new TH1D("hRecoPV","Prim. vert. reconstruction (yes or no);reco (0=no, 1=yes) ;",2,-0.5,1.5);
    fRecoPV->Sumw2();
    fOutput->Add(fRecoPV);
    
    fPVZBefore = new TH1D ("hPVZDistr","Z Distribution for the Reconstructed Vertex",200,0.,40.);
    fPVZBefore->Sumw2();
    fOutput->Add(fPVZBefore);
    
    fEtaPhiCell = new TH2D ("hEtaPhiCellActivity","",250,0.,1000., 250, 0., 1000.);
    fEtaPhiCell->Sumw2();
    fOutput->Add(fEtaPhiCell);
    
    fEtaPhiClus = new TH2D ("hEtaPhiClusActivity","",250,0.,1000., 250, 0., 1000.);
    fEtaPhiClus->Sumw2();
    fOutput->Add(fEtaPhiClus);
    
    fClusEvsClusT = new TH2D ("hEnVSTime","",250,0.,1000., 250,0.,1000.);
    fClusEvsClusT->Sumw2();
    fOutput->Add(fClusEvsClusT);
    
    fDeltaETAClusTrackVSpT = new TH2D("hTC_Dz","Track-Cluster Dz vs pT of Cluster",100,0.,100.,100,-0.05,0.05);
    fDeltaETAClusTrackVSpT->Sumw2();
    fOutput->Add(fDeltaETAClusTrackVSpT);
    
    fDeltaPHIClusTrackVSpT = new TH2D("hTC_Dx","Track-Cluster Dx vs pT of Cluster",100,0.,100.,100,-0.05,0.05);
    fDeltaPHIClusTrackVSpT->Sumw2();
    fOutput->Add(fDeltaPHIClusTrackVSpT);
    
  }
  //Initialize only the Common THistos for the Three different output
  
  fGoodEventsOnPVZ = new TH1D ("hGOODwrtPVZ","Number of Selected Events wrt Cut on Primary Vertex Z (0=disregarded,1=selected)",2,0.,2.);
  fGoodEventsOnPVZ->Sumw2();
  fOutput->Add(fGoodEventsOnPVZ);
  
  fPT = new TH1D("hPt_NC","P_{T} distribution for Neutral Clusters",100,0.,100.);
  fPT->Sumw2();
  fOutput->Add(fPT);
  
  fM02 = new TH2D("hM02_NC","M02 distribution for Neutral Clusters vs E",100,0.,100.,500,0.,5.);
  fM02->Sumw2();
  fOutput->Add(fM02);
  
  fNLM = new TH1D("hNLM_NC","NLM distribution for Neutral Clusters",200,0.,4.);
  fNLM->Sumw2();
  fOutput->Add(fNLM);
  
  fEtIsoCells = new TH1D("hEtIsoCell_NC","E_{T}^{iso cone} in iso cone distribution for Neutral Clusters with EMCAL Cells",100,0.,100.);
  fEtIsoCells->SetXTitle("#Sigma E_{T}^{iso cone} (GeV/c)");
  fEtIsoCells->Sumw2();
  fOutput->Add(fEtIsoCells);
  
  fEtIsoClust = new TH1D("hEtIsoClus_NC","#Sigma E_{T}^{iso cone} in iso cone distribution for Neutral Clusters with EMCAL Clusters",100,0.,100.);
  fEtIsoClust->SetXTitle("#Sigma E_{T}^{iso cone} (GeV/c)");
  fEtIsoClust->Sumw2();
  fOutput->Add(fEtIsoClust);
  
  fPtIsoTrack = new TH1D("hPtIsoTrack_NC"," #Sigma P_{T}^{iso cone} in iso cone distribution for Neutral Clusters with Tracks",100,0.,100.);
  fPtIsoTrack->SetXTitle("#Sigma P_{T}^{iso cone} (GeV/c)");
  fPtIsoTrack->Sumw2();
  fOutput->Add(fPtIsoTrack);
  
  fPtEtIsoTC = new TH1D("hPtEtIsoTrackClust_NC","#Sigma P_{T}^{iso cone} + #Sigma E_{T}^{iso cone} in iso cone distribution for Neutral Clusters with Tracks and Clusters",100,0.,100.);
  fPtEtIsoTC->SetXTitle("#Sigma P_{T}^{iso cone} + #Sigma E_{T}^{iso cone} (GeV/c)");
  fPtEtIsoTC->Sumw2();
  fOutput->Add(fPtEtIsoTC);
  
  fPhiBandUEClust = new TH2D(Form("hPhiBandUE_Cluster"),Form("UE Estimation with Phi Band Clusters"),100,0.,100.,100,0.,100.);
  fPhiBandUEClust->SetXTitle("E_{T}");
  fPhiBandUEClust->SetYTitle("#Sigma E_{T}^{UE}");
  fPhiBandUEClust->Sumw2();
  fOutput->Add(fPhiBandUEClust);
  
  fEtaBandUEClust = new TH2D(Form("hEtaBandUE_Cluster"),Form("UE Estimation with Phi Band Clusters"),100,0.,100.,100,0.,100.);
  fEtaBandUEClust->SetXTitle("E_{T}");
  fEtaBandUEClust->SetYTitle("#Sigma E_{T}^{UE}");
  fEtaBandUEClust->Sumw2();
  fOutput->Add(fEtaBandUEClust);
  
  fPhiBandUECells = new TH2D(Form("hPhiBandUE_CELLS"),Form("UE Estimation with Phi Band CELLS"),100,0.,100.,100,0.,100.);
  fPhiBandUECells->SetXTitle("E_{T}");
  fPhiBandUECells->SetYTitle("#Sigma E_{T}^{UE}");
  fPhiBandUECells->Sumw2();
  fOutput->Add(fPhiBandUECells);
  
  fEtaBandUECells = new TH2D(Form("hEtaBandUE_CELLS"),Form("UE Estimation with Phi Band and CELLS"),100,0.,100.,100,0.,100.);
  fEtaBandUECells->SetXTitle("E_{T}");
  fEtaBandUECells->SetYTitle("#Sigma E_{T}^{UE}");
  fEtaBandUECells->Sumw2();
  fOutput->Add(fEtaBandUECells);
  
  fPhiBandUETracks = new TH2D(Form("hPhiBandUE_TPC"),Form("UE Estimation with Phi Band TPC "),100,0.,100.,100,0.,100.);
  fPhiBandUETracks->SetXTitle("E_{T}");
  fPhiBandUETracks->SetYTitle("#Sigma P_{T}^{UE}");
  fPhiBandUETracks->Sumw2();
  fOutput->Add(fPhiBandUETracks);
  
  fEtaBandUETracks = new TH2D(Form("hEtaBandUE_TPC"),Form("UE Estimation with Phi Band and TPC"),100,0.,100.,100,0.,100.);
  fEtaBandUETracks->SetXTitle("E_{T}");
  fEtaBandUETracks->SetYTitle("#Sigma P_{T}^{UE}");
  fEtaBandUETracks->Sumw2();
  fOutput->Add(fEtaBandUETracks);
  
  fPerpConesUETracks = new TH2D("hConesUE","UE Estimation with Perpendicular Cones in TPC",100,0.,100.,100,0.,100.);
  fPerpConesUETracks->SetXTitle("E_{T}");
  fPerpConesUETracks->SetYTitle("#Sigma P_{T}^{UE}");
  fPerpConesUETracks->Sumw2();
  fOutput->Add(fPerpConesUETracks);
  
  fTPCWithoutIsoConeB2BbandUE = new TH2D("hFullTPCUE","UE Estimation with almost Full TPC",100,0.,100.,100,0.,100.);
  fPhiBandUEClust->SetXTitle("E_{T}");
  fPhiBandUEClust->SetYTitle("#Sigma E_{T}^{UE}");
  fTPCWithoutIsoConeB2BbandUE->Sumw2();
  fOutput->Add(fTPCWithoutIsoConeB2BbandUE);
  
  fEtIsolatedClust = new TH1D("hEtIsolatedClust","E_{T} distribution for Isolated Photons with clusters; #Sigma E_{T}^{iso cone}<Ethres",100,0.,100.);
  fEtIsolatedClust->SetXTitle("E_{T}^{iso}");
  fEtIsolatedClust->Sumw2();
  fOutput->Add(fEtIsolatedClust);
  
  fEtIsolatedCells = new TH1D("hEtIsolatedCells","E_{T} distribution for Isolated Photons with cells; #Sigma E_{T}^{iso cone}<Ethres",100,0.,100.);
  fEtIsolatedCells->SetXTitle("E_{T}^{iso}");
  fEtIsolatedCells->Sumw2();
  fOutput->Add(fEtIsolatedCells);
  
  fEtIsolatedTracks = new TH1D("hEtIsolatedTracks","E_{T} distribution for Isolated Photons with tracks; #Sigma P_{T}^{iso cone}<Pthres",100,0.,100.);
  fEtIsolatedTracks->SetXTitle("E_{T}^{iso}");
  fEtIsolatedTracks->Sumw2();
  fOutput->Add(fEtIsolatedTracks);
  
  fTest = new TH1D ("hTest","test cluster collection",100,-2.,6.);
  fTest->Sumw2();
  fOutput->Add(fTest);
  
  
  PostData(1, fOutput);
  //   return;
}

//________________________________________________________________________
Float_t* AliAnalysisTaskEMCALPhotonIsolation::GenerateFixedBinArray(Int_t n, Float_t min, Float_t max) const
{
  // Generate the bin array for the ThnSparse
  
  Float_t *bins = new Float_t[n+1];
  
  Float_t binWidth = (max-min)/n;
  bins[0] = min;
  for (Int_t i = 1; i <= n; i++) {
    bins[i] = bins[i-1]+binWidth;
  }
  
  return bins;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::ExecOnce()
{
  //   Init the analysis.
  
  
  
  if (fParticleCollArray.GetEntriesFast()<2) {
    AliError(Form("Wrong number of particle collections (%d), required 2",fParticleCollArray.GetEntriesFast()));
    return;
  }
  
  
  //  for (Int_t i = 0; i < 2; i++) {
  //    AliParticleContainer *contain = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
  //    // contain->GetClassName("AliEmcalParticle");
  //  }
  
  
  
  AliAnalysisTaskEmcal::ExecOnce();
  if (!fInitialized) {
    
    AliError(Form("AliAnalysisTask not initialized"));
    return;
  }
  
  fevents+=1;
  
}

//______________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::Run()
{
  // Run the analysis
  
  
  AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliEmcalParticle *emccluster = 0;
  
  // delete output USEFUL LATER FOR CONTAINER CREATION !!
  //fOutClusters->Delete();
  
  //Int_t clusCount = 0;    AliError(Form("Should be here each time"));
  // loop over all clusters
  clusters->ResetCurrentID();
  Int_t index=-1;
  
  //Double_t ETleadingclust = 0., M02leadingcluster = 0., lambda0cluster = 0., phileadingclust = 0., etaleadingclust = 0., ptmc = 0.,mcptsum = 0.;
  //Int_t Ntracks;
  //Definition of the Array for Davide's Output
  //const Int_t ndims =   fNDimensions;
  //Double_t outputValues[ndims];
  
  
  if (fisLCAnalysis) {
    //   AliError(Form("Check the loop"));
    // get the leading particle
    
    
    emccluster = static_cast<AliEmcalParticle*>(clusters->GetLeadingParticle());
    
    if(!emccluster){
      
      AliError(Form("no leading one"));
      return kTRUE;
    }
    
    
    
    //     emccluster = clusters->GetLeadingParticle();
    index = emccluster->IdInCollection();
    
    //add a command to get the index of the leading cluster!
    if (!emccluster->IsEMCAL()) return kTRUE;
    
    AliVCluster *coi = emccluster->GetCluster();
    if (!coi) return kTRUE;
    
    TLorentzVector vecCOI;
    coi->GetMomentum(vecCOI,fVertex);
    
    
    if(!CheckBoundaries(vecCOI))
      return kTRUE;
    
    if(ClustTrackMatching(coi))
      return kTRUE;
    
    else
      FillGeneralHistograms(coi,vecCOI, index);
  }
  else{
    //get the entries of the Cluster Container
    //whatever is a RETURN in LCAnalysis here is a CONTINUE,
    //since there are more than 1 Cluster per Event
    
    while((emccluster = static_cast<AliEmcalParticle*>(clusters->GetNextAcceptParticle()))){
      index++;
      if (!emccluster->IsEMCAL()) return kTRUE;
      
      AliVCluster *coi = emccluster->GetCluster();
      if(!coi) return kTRUE;
      
      TLorentzVector vecCOI;
      coi->GetMomentum(vecCOI,fVertex);
      
      if(!CheckBoundaries(vecCOI))
        return kTRUE;
      
      if(ClustTrackMatching(coi))
        continue;
      else
        FillGeneralHistograms(coi,vecCOI, index);
      
    }
  }
  //  PostData(1, fOutput);
  return kTRUE;
}


//__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::ClustTrackMatching(AliVCluster *cluster) {
  // Check if the cluster match to a track
  
  Double_t deta=999;
  Double_t dphi=999;
  
  AliParticleContainer *tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  // for an incoming cluster reject it if there is a track corresponding with a deta and dphi lower than the cuts
  TLorentzVector nPart;
  cluster->GetMomentum(nPart, fVertex);
  
  AliVTrack *mt = NULL;
  AliAODCaloCluster *acl = dynamic_cast<AliAODCaloCluster*>(cluster);
  if(acl) {
    if(acl->GetNTracksMatched()>1)
      mt = static_cast<AliVTrack*>(acl->GetTrackMatched(0));
  }
  else {
    AliESDCaloCluster *ecl = dynamic_cast<AliESDCaloCluster*>(cluster);
    if(!ecl){
      AliError("ClusTrack matching did not work");
      return kTRUE;
    }
    Int_t im = ecl->GetTrackMatchedIndex();
    if(tracks && im>=0) {
      mt = static_cast<AliVTrack*>(tracks->GetParticle(im));
    }
  }
  //  if(mt && mt->TestFilterBit(768)) { //Hybrid tracks with AOD
  AliPicoTrack::GetEtaPhiDiff(mt, cluster, dphi, deta);
  fDeltaETAClusTrackVSpT->Fill(nPart.Pt(), deta);
  fDeltaPHIClusTrackVSpT->Fill(nPart.Pt(), dphi);
  
  
  if(mt) return kTRUE;
  
  else return kFALSE;
  
}



//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoCellPhiBand(TLorentzVector c, Float_t &etIso, Float_t &phiBandcells){
  // Underlying events study with EMCAL cells in phi band
  
  AliEMCALGeometry* emcalGeom = AliEMCALGeometry::GetInstance();
  
  Float_t sumEnergyPhiBandCells=0., sumEnergyConeCells=0.;
  
  
  // check the cell corresponding to the leading cluster
  Int_t absId = 999;
  //maybe best to call it LeadingCellIdinClus or better maxId since it won't be used anywhere else ???
  Bool_t cellLeadingClustID = emcalGeom->GetAbsCellIdFromEtaPhi(c.Eta(),c.Phi(),absId);
  if(!cellLeadingClustID) return;
  
  else{
    Int_t iTower = -1;
    Int_t iModule = -1;
    Int_t imEta=-1, imPhi=-1;
    Int_t iEta =-1, iPhi =-1;
    
    emcalGeom->GetCellIndex(absId,iModule,iTower,imPhi,imEta); // to get the module, the tower, eta and phi for the cell corresponding to the leading cluster
    emcalGeom->GetCellPhiEtaIndexInSModule(iModule,iTower,imPhi,imEta,iPhi,iEta); // to get the cell eta and phi in the super module for the cell coming from the leading cluster
    
    // Get the row and the column of the cell corresponding to the leading cluster in EMCAL
    Int_t colCellLeadingClust = iEta;
    if(iModule % 2) colCellLeadingClust = AliEMCALGeoParams::fgkEMCALCols + iEta ; // if the SM number is even you need to translate to have the one corresponding in EMCAL
    Int_t rowCellLeadingClust = iPhi + AliEMCALGeoParams::fgkEMCALRows*int(iModule/2); // to have the corresponding row in EMCAL
    
    // total number or rows and columns in EMCAL
    Int_t nTotalRows = AliEMCALGeoParams::fgkEMCALRows*16/3 ;  // 5 + 2/3 supermodules in a row
    Int_t nTotalCols = 2*AliEMCALGeoParams::fgkEMCALCols; // 2 supermodules in a column
    
    Int_t nbConeSize = int(fIsoConeRadius/0.0143); //0.0143 tower acceptance
    
    // Get the cells
    AliVCaloCells * cells =InputEvent()->GetEMCALCells();
    
    // define the max and min row and column corresponding to the isolation cone around the seed cell from the leading cluster
    Int_t iRowMinCone = rowCellLeadingClust-nbConeSize;
    if(iRowMinCone<0) iRowMinCone=0;
    
    Int_t iRowMaxCone = rowCellLeadingClust+nbConeSize;
    if(iRowMaxCone>AliEMCALGeoParams::fgkEMCALRows) iRowMaxCone=AliEMCALGeoParams::fgkEMCALRows;  // AliEMCALGeoParams::fgkEMCALRows = 24 in a supermodule
    
    Int_t iColMinCone = colCellLeadingClust - nbConeSize;
    if(iColMinCone<0) iColMinCone = 0;
    
    Int_t iColMaxCone = colCellLeadingClust+nbConeSize;
    if(iColMaxCone>AliEMCALGeoParams::fgkEMCALCols) iColMaxCone=AliEMCALGeoParams::fgkEMCALCols;  // AliEMCALGeoParams::fgkEMCALCols = 48 in a supermodule
    
    // loop on all cells
    for(Int_t iCol=0; iCol<nTotalCols; iCol++){
      for(Int_t iRow=0; iRow<nTotalRows; iRow++){
        // now recover the cell indexes in a supermodule
        Int_t iSector = int(iRow/AliEMCALGeoParams::fgkEMCALRows); // check in which SM is the cell
        if(iSector==5) continue;  //
        Int_t inModule = -1;
        Int_t iColLoc  = -1;
        if(iCol < AliEMCALGeoParams::fgkEMCALCols){ // if the SM number is odd the column is the one corresponding in the supermodule
          inModule = 2*iSector + 1;
          iColLoc  = iCol;
        }
        else if(iCol > AliEMCALGeoParams::fgkEMCALCols - 1){ // if the SM number is even the column isn't the one corresponding in the supermodule
          inModule = 2*iSector;
          iColLoc  = iCol-AliEMCALGeoParams::fgkEMCALCols;
        }
        
        Int_t iRowLoc  = iRow - AliEMCALGeoParams::fgkEMCALRows*iSector ;
        
        if(TMath::Abs(iCol-colCellLeadingClust)<nbConeSize && TMath::Abs(iCol+colCellLeadingClust)>nbConeSize){
          if(iRow<iRowMaxCone && iRow>iRowMinCone){
            Int_t iabsId = emcalGeom->GetAbsCellIdFromCellIndexes(inModule,iRowLoc,iColLoc);  // verifier les iRow et iCol
            sumEnergyPhiBandCells+=cells->GetAmplitude(iabsId); //should be Et
          }
        }
        else if (TMath::Abs(iCol-colCellLeadingClust)>nbConeSize && TMath::Abs(iCol+colCellLeadingClust)<nbConeSize){
          Int_t iabsId = emcalGeom->GetAbsCellIdFromCellIndexes(inModule,iRowLoc,iColLoc);  // verifier les iRow et iCol
          sumEnergyConeCells+=cells->GetAmplitude(iabsId); //should be Et
        }
      }
    }
  }
  etIso = sumEnergyConeCells;
  phiBandcells = sumEnergyPhiBandCells;
}


//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoCellEtaBand(TLorentzVector c, Float_t &etIso, Float_t &etaBandcells){
  // Underlying events study with EMCAL cell in eta band
  
  
  AliEMCALGeometry* emcalGeom = AliEMCALGeometry::GetInstance();
  
  Float_t sumEnergyEtaBandCells=0., sumEnergyConeCells=0.;
  
  
  
  // check the cell corresponding to the leading cluster
  Int_t absId = 999;
  //maybe best to call it LeadingCellIdinClus or better maxId since it won't be used anywhere else ???
  Bool_t cellLeadingClustID = emcalGeom->GetAbsCellIdFromEtaPhi(c.Eta(),c.Phi(),absId);
  if(!cellLeadingClustID) return;
  
  else{
    
    Int_t iTower = -1;
    Int_t iModule = -1;
    Int_t imEta=-1, imPhi=-1;
    Int_t iEta =-1, iPhi =-1;
    
    emcalGeom->GetCellIndex(absId,iModule,iTower,imPhi,imEta); // to get the module, the tower, eta and phi for the cell corresponding to the leading cluster
    emcalGeom->GetCellPhiEtaIndexInSModule(iModule,iTower,imPhi,imEta,iPhi,iEta); // to get the cell eta and phi in the super module for the cell coming from the leading cluster
    
    // Get the row and the column of the cell corresponding to the leading cluster in EMCAL
    Int_t colCellLeadingClust = iEta;
    if(iModule % 2) colCellLeadingClust = AliEMCALGeoParams::fgkEMCALCols + iEta ; // if the SM number is even you need to translate to have the one corresponding in EMCAL
    Int_t rowCellLeadingClust = iPhi + AliEMCALGeoParams::fgkEMCALRows*int(iModule/2); // to have the corresponding row in EMCAL
    
    // total number or rows and columns in EMCAL
    Int_t nTotalRows = AliEMCALGeoParams::fgkEMCALRows*16/3 ;  // 5 + 2/3 supermodules in a row
    Int_t nTotalCols = 2*AliEMCALGeoParams::fgkEMCALCols; // 2 supermodules in a column
    
    Int_t nbConeSize = int(fIsoConeRadius/0.0143); //0.0143 tower acceptance
    
    // Get the cells
    AliVCaloCells * cells =InputEvent()->GetEMCALCells();
    
    // define the max and min row and column corresponding to the isolation cone around the seed cell from the leading cluster
    Int_t iRowMinCone = rowCellLeadingClust-nbConeSize;
    if(iRowMinCone<0) iRowMinCone=0;
    
    Int_t iRowMaxCone = rowCellLeadingClust+nbConeSize;
    if(iRowMaxCone>AliEMCALGeoParams::fgkEMCALRows) iRowMaxCone=AliEMCALGeoParams::fgkEMCALRows;  // AliEMCALGeoParams::fgkEMCALRows = 24 in a supermodule
    
    Int_t iColMinCone = colCellLeadingClust-nbConeSize;
    if(iColMinCone<0) iColMinCone = 0;
    
    Int_t iColMaxCone = colCellLeadingClust+nbConeSize;
    if(iColMaxCone>AliEMCALGeoParams::fgkEMCALCols) iColMaxCone=AliEMCALGeoParams::fgkEMCALCols;  // AliEMCALGeoParams::fgkEMCALCols = 48 in a supermodule
    
    // loop on all cells
    for(Int_t iCol=0; iCol<nTotalCols; iCol++)
    {
      for(Int_t iRow=0; iRow<nTotalRows; iRow++)
      {
        // now recover the cell indexes in a supermodule
        Int_t iSector = int(iRow/AliEMCALGeoParams::fgkEMCALRows); // check in which SM is the cell
        if(iSector==5) continue;  //
        Int_t inModule = -1;
        Int_t iColLoc  = -1;
        if(iCol < AliEMCALGeoParams::fgkEMCALCols){ // if the SM number is odd the column is the one corresponding in the supermodule
          inModule = 2*iSector + 1;
          iColLoc  = iCol;
        }
        else if(iCol > AliEMCALGeoParams::fgkEMCALCols - 1){ // if the SM number is even the column isn't the one corresponding in the supermodule
          inModule = 2*iSector;
          iColLoc  = iCol-AliEMCALGeoParams::fgkEMCALCols;
        }
        
        Int_t iRowLoc  = iRow - AliEMCALGeoParams::fgkEMCALRows*iSector ;
        
        if(TMath::Abs(iCol-colCellLeadingClust)<nbConeSize && TMath::Abs(iCol+colCellLeadingClust)>nbConeSize){
          if(iCol<iColMaxCone && iCol>iColMinCone){
            Int_t iabsId = emcalGeom->GetAbsCellIdFromCellIndexes(inModule,iRowLoc,iColLoc);  // verifier les iRow et iCol
            sumEnergyEtaBandCells+=cells->GetAmplitude(iabsId); //should be Et
          }
        }
        else if (TMath::Abs(iCol-colCellLeadingClust)>nbConeSize && TMath::Abs(iCol+colCellLeadingClust)<nbConeSize){
          Int_t iabsId = emcalGeom->GetAbsCellIdFromCellIndexes(inModule,iRowLoc,iColLoc);  // verifier les iRow et iCol
          sumEnergyConeCells+=cells->GetAmplitude(iabsId); //should be Et
        }
      }
    }
  }
  etIso = sumEnergyConeCells;
  etaBandcells = sumEnergyEtaBandCells;
}


//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoClusPhiBand(TLorentzVector c, Float_t &etIso, Float_t &phiBandclus, Int_t index){
  // Underlying events study with clusters in phi band
  
  Float_t sumEnergyPhiBandClus=0., sumEnergyConeClus=0.;
  
  //needs a check on the same cluster
  AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliEmcalParticle *clust;
  
  Int_t localIndex=-1;
  while((clust = static_cast<AliEmcalParticle*>(clusters->GetNextAcceptParticle()))){ //check the position of other clusters in respect to the leading cluster
    localIndex++;
    if(localIndex==index) continue;
    
    AliVCluster *cluster= clust->GetCluster();
    
    TLorentzVector nClust; //STILL NOT INITIALIZED
    cluster->GetMomentum(nClust,fVertex);
    Float_t phiClust =nClust.Phi();
    Float_t etaClust= nClust.Eta();
    //redefine phi/c.Eta() from the cluster we passed to the function
    
    Float_t  radius = TMath::Sqrt(TMath::Power(phiClust- c.Phi(),2)+TMath::Power(etaClust-c.Eta(),2)); // define the radius between the leading cluster and the considered cluster
    
    if(radius>fIsoConeRadius){ // the cluster is outside the isolation cone in this case study UE
      
      // actually phi band here
      if(TMath::Abs(etaClust - c.Eta()) < fIsoConeRadius){
        sumEnergyPhiBandClus += nClust.Pt();
      }
    }
    else // if the cluster is in the isolation cone, add the cluster pT
      sumEnergyConeClus += nClust.Pt();
    
  }
  etIso = sumEnergyConeClus;
  phiBandclus = sumEnergyPhiBandClus;
}


//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoClusEtaBand(TLorentzVector c, Float_t &etIso, Float_t &etaBandclus, Int_t index){
  // Underlying events study with clusters in eta band
  
  Float_t sumEnergyEtaBandClus =0., sumEnergyConeClus=0.;
  AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  
  AliEmcalParticle *clust;
  
  Int_t localIndex=-1;
  while((clust = static_cast<AliEmcalParticle*>(clusters->GetNextAcceptParticle()))){ //check the position of other clusters in respect to the leading cluster
    localIndex++;
    if(localIndex==index) continue;
    
    AliVCluster *cluster= clust->GetCluster();
    
    TLorentzVector nClust; //STILL NOT INITIALIZED
    cluster->GetMomentum(nClust,fVertex);
    
    Float_t phiClust =nClust.Phi();
    Float_t etaClust= nClust.Eta();
    //redefine phi/c.Eta() from the cluster we passed to the function
    
    // define the radius between the leading cluster and the considered cluster
    Float_t  radius = TMath::Sqrt(TMath::Power(phiClust-c.Phi(),2)+TMath::Power(etaClust-c.Eta(),2));
    
    if(radius>fIsoConeRadius){ // the cluster is outside the isolation cone in this case study UE
      
      // actually eta band here
      if(TMath::Abs(etaClust - c.Eta()) < fIsoConeRadius){
        sumEnergyEtaBandClus += nClust.Pt();
      }
    }
    else  // if the cluster is in the isolation cone, add the cluster pT
      sumEnergyConeClus += nClust.Pt();
    
  }
  etIso = sumEnergyConeClus;
  etaBandclus = sumEnergyEtaBandClus;
}


//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackPhiBand(TLorentzVector c, Float_t &ptIso, Float_t &phiBandtrack){
  // Underlying events study with tracks in phi band in EMCAL acceptance
  
  //INSERT BOUNDARIES ACCORDING TO THE FLAG WE WANT TO USE
  Float_t sumpTConeCharged=0., sumpTPhiBandTrack=0.;
  Float_t minPhi= 0., maxPhi= 2*TMath::Pi(), minEta = -0.9, maxEta= 0.9;
  
  AliParticleContainer *tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  
  if(!fTPC4Iso){
    minEta = -0.7;
    maxEta = 0.7;
    minPhi = 1.4;
    maxPhi = TMath::Pi();
  }
  
  AliVParticle *track = tracks->GetNextAcceptParticle(0);
  
  while(track){
    
    //CHECK IF TRACK IS IN BOUNDARIES
    Float_t phiTrack = track->Phi();
    Float_t etaTrack = track->Eta();
    // define the radius between the leading cluster and the considered cluster
    //redefine phi/c.Eta() from the cluster we passed to the function
    if(phiTrack < maxPhi && phiTrack > minPhi && etaTrack < maxEta && etaTrack > minEta){
      Float_t  radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2));
      
      if(radius>fIsoConeRadius){ // if tracks are outside the isolation cone study
        
        // actually phi band here --- ADD Boundaries conditions
        if(TMath::Abs(etaTrack - c.Eta()) < fIsoConeRadius){
          sumpTPhiBandTrack += track->Pt();
        }
      }
      else
        sumpTConeCharged+=track->Pt(); // should not double count if the track matching is already done
    }
    track=tracks->GetNextAcceptParticle();
  }
  ptIso = sumpTConeCharged;
  phiBandtrack = sumpTPhiBandTrack;
}


//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackEtaBand(TLorentzVector c, Float_t &ptIso, Float_t &etaBandtrack){
  // Underlying events study with tracks in eta band in EMCAL acceptance
  
  //INSERT BOUNDARIES ACCORDING TO THE FLAG WE WANT TO USE
  Float_t sumpTConeCharged=0., sumpTEtaBandTrack=0.;
  Float_t minPhi= 0., maxPhi= 2*TMath::Pi(), minEta = -0.9, maxEta= 0.9;
  
  AliParticleContainer *tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  
  
  if(!fTPC4Iso){
    minEta = -0.7;
    maxEta = 0.7;
    minPhi = 1.4;
    maxPhi = TMath::Pi();
  }
  
  
  AliVParticle *track = tracks->GetNextAcceptParticle(0);
  
  while(track){
    
    Float_t phiTrack = track->Phi();
    Float_t etaTrack = track->Eta();
    //redefine phi/c.Eta() from the cluster we passed to the function
    if(phiTrack < maxPhi && phiTrack > minPhi && etaTrack < maxEta && etaTrack > minEta){
      Float_t  radius = TMath::Sqrt(TMath::Power(phiTrack-c.Phi(),2)+TMath::Power(etaTrack-c.Eta(),2)); // define the radius between the leading cluster and the considered cluster
      if(radius>fIsoConeRadius){ // if tracks are outside the isolation cone study UE
        
        // actually eta band here --- ADD Boundaries conditions
        if(TMath::Abs(phiTrack - c.Phi()) < fIsoConeRadius){
          sumpTEtaBandTrack += track->Pt();
        }
      }
      else
        sumpTConeCharged += track->Pt(); // should not double count if the track matching is already done
    }
    track=tracks->GetNextAcceptParticle();
  }
  ptIso = sumpTConeCharged;
  etaBandtrack = sumpTEtaBandTrack;
}


//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackOrthCones(TLorentzVector c, Float_t &ptIso, Float_t &cones){
  // Underlying events study with tracks in orthogonal cones in TPC
  
  Float_t sumpTConeCharged=0., sumpTPerpConeTrack=0.;
  
  AliParticleContainer *tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  
  Float_t etaClus = c.Eta();
  Float_t phiClus = c.Phi();
  Float_t phiCone1 = phiClus - TMath::PiOver2();
  Float_t phiCone2 = phiClus + TMath::PiOver2();
  
  if (phiCone1 < 0.) phiCone1 += 2*TMath::Pi();
  
  AliVParticle *track = tracks->GetNextAcceptParticle(0);
  
  while(track){
    
    Float_t phiTrack = track->Phi();
    Float_t etaTrack = track->Eta();
    
    Float_t dist2Clust = TMath::Sqrt(TMath::Power(etaTrack-etaClus, 2)+TMath::Power(phiTrack-phiClus, 2));
    if (dist2Clust<fIsoConeRadius)// if tracks are inside the IsoCone
      sumpTConeCharged += track->Pt();
    
    else{//tracks outside the IsoCone
      //Distances from the centres of the two Orthogonal Cones
      Float_t dist2Cone1 = TMath::Sqrt(TMath::Power(etaTrack-etaClus, 2)+TMath::Power(phiTrack-phiCone1, 2));
      Float_t dist2Cone2 = TMath::Sqrt(TMath::Power(etaTrack-etaClus, 2)+TMath::Power(phiTrack-phiCone2, 2));
      
      //Is the Track Inside one of the two Cones ->Add to UE
      if((dist2Cone1 < fIsoConeRadius) || (dist2Cone2 < fIsoConeRadius))
        sumpTPerpConeTrack += track->Pt();
    }
    track=tracks->GetNextAcceptParticle();
  }
  ptIso = sumpTConeCharged;
  cones = sumpTPerpConeTrack;
}

//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackFullTPC(TLorentzVector c, Float_t &ptIso, Float_t &full){
  // Underlying events study with tracks in full TPC except back to back bands
  
  Float_t sumpTConeCharged=0., sumpTTPCexceptB2B=0.;
  
  AliParticleContainer *tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  AliVParticle *track = tracks->GetNextAcceptParticle(0);
  
  while(track){
    
    Float_t phiTrack = track->Phi();
    Float_t etaTrack = track->Eta();
    //redefine phi/c.Eta() from the cluster we passed to the function
    
    Float_t  radius = TMath::Sqrt(TMath::Power(phiTrack-c.Phi(),2)+TMath::Power(etaTrack-c.Eta(),2)); // define the radius between the leading cluster and the considered cluster
    if(radius>fIsoConeRadius){ // if tracks are outside the isolation cone study UE
      Double_t dphiUp = c.Phi() + TMath::Pi() - fIsoConeRadius;
      Double_t dphiDown = c.Phi() + TMath::Pi() + fIsoConeRadius;
      // TPC except B2B
      if(TMath::Power(phiTrack-c.Phi(),2) +TMath::Power(etaTrack-c.Eta(),2) > TMath::Power((fIsoConeRadius+0.1),2) && phiTrack < dphiDown && phiTrack> dphiUp)
        sumpTTPCexceptB2B += track->Pt();
    }
    else
      sumpTConeCharged += track->Pt(); // should not double count if the track matching is already done
    
    track=tracks->GetNextAcceptParticle();
  }
  ptIso = sumpTConeCharged;
  full = sumpTTPCexceptB2B;
}

//__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::CheckBoundaries(TLorentzVector vecCOI){
  // Check if the cone around the considered cluster is in EMCAL acceptance
  
  Float_t minPhiBound= 1.5 , minEtaBound= -0.5, maxPhiBound= TMath::Pi()-0.1, maxEtaBound= 0.5;
  Bool_t isINBoundaries;
  
  if(!fTPC4Iso){
    minPhiBound = 1.4+0.4; //to be changed with fIsoConeR
    maxPhiBound = TMath::Pi()-0.4;
    minEtaBound = -0.7+0.4;
    maxEtaBound = 0.7-0.4;
  }
  
  if(vecCOI.Eta() > maxEtaBound || vecCOI.Eta() < minEtaBound || vecCOI.Phi() > maxPhiBound || vecCOI.Phi() >minPhiBound)
    isINBoundaries=kFALSE;
  else
    isINBoundaries=kTRUE;
  
  return isINBoundaries;
}

//__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::FillGeneralHistograms(AliVCluster *coi, TLorentzVector vecCOI, Int_t index){
  // Fill the histograms for underlying events and isolation studies
  
  AliParticleContainer *tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  
  AliEmcalParticle *emcTrack = 0;
  
  int nTracks=0;
  tracks->ResetCurrentID();
  while ((emcTrack = static_cast<AliEmcalParticle*>(tracks->GetNextAcceptParticle()))) {
    AliVTrack *track = emcTrack->GetTrack();
    if(!track) continue;
    // if(!(track->TestFilterBit("kHybrid"))) continue;
    
    nTracks++;
  }
  
  Double_t eTCOI = 0., m02COI = 0., lambda0cluster = 0., ptmc = 0., mcptsum = 0.;
  //Int_t Ntracks;
  //Definition of the Array for Davide's Output
  const Int_t ndims =   fNDimensions;
  Double_t outputValues[ndims];
  
  eTCOI = vecCOI.Et();
  m02COI = coi->GetM02();

  
  // ******** Isolation and UE calculation with different methods *********
  
  Double_t eTThreshold = 5;
  
  switch(fEtIsoMethod)
  {
    case 0:  // SumEt<EtThr
      if(fM02mincut < m02COI && m02COI < fM02maxcut)  // photon candidate, cuts have to be decided after studies
      {
        eTThreshold = fEtIsoThreshold;
      }
      break;
      
    case 1:  // SumEt<%Ephoton
      if(fM02mincut < m02COI && m02COI < fM02maxcut) // photon candidate, cuts have to be decided after studies
      {
        eTThreshold = fEtIsoThreshold * eTCOI;
      }
      break;
      
    case 2: // Etmax<eTThreshold
      eTThreshold = fEtIsoThreshold;
      if(fM02mincut < m02COI && m02COI < fM02maxcut && eTCOI<eTThreshold) // photon candidate, cuts have to be decided after studies
      {
        fEtIsolatedClust->Fill(eTCOI);
      }
      break;
  }
  
  //DO NOT CHANGE EVER AGAIN THE FOLLOWING DEFINITIONS
  Float_t isoConeArea = TMath::Pi()*fIsoConeRadius*fIsoConeRadius;
  Float_t etaBandArea = (1.4 - 2.*fIsoConeRadius)*2.*fIsoConeRadius;
  Float_t phiBandArea = ((5./9.)*TMath::Pi()- 2.*fIsoConeRadius)*2.*fIsoConeRadius;
  Float_t etaBandAreaTr = (1.8 - 2.*fIsoConeRadius)*2.*fIsoConeRadius;
  Float_t phiBandAreaTr = (2.*TMath::Pi() - 4.*fIsoConeRadius)*2.*fIsoConeRadius;//there is a 2 more because of the JET CONE B2B
  Float_t perpConesArea = 2.*isoConeArea;
  Float_t fullTPCArea = 1.8*2.*TMath::Pi()-fIsoConeRadius*(TMath::Pi()*fIsoConeRadius + 3.6);
  
  Float_t isolation, ue;
  
  if(!fTPC4Iso){ //EMCAL Only for Acceptance of Clusters
    switch(fIsoMethod)
    {
      case 0: //EMCAL CELLS
        
        switch (fUEMethod)
      {
        case 0: //phi band
          EtIsoCellPhiBand(vecCOI, isolation, ue);
          //Normalisation ue wrt Area - TO DO-
          ue = ue * (isoConeArea / phiBandArea);
          fPhiBandUECells->Fill(vecCOI.Pt() , ue);
          fEtIsoCells->Fill(isolation);
          if(isolation<eTThreshold)
          {
            fEtIsolatedCells->Fill(eTCOI);
            fEtisoT=eTCOI;
            fPtisoT=vecCOI.Pt();
          }
          break;
        case 1: //eta band
          EtIsoCellEtaBand(vecCOI, isolation, ue);
          //Normalisation ue wrt Area - TO DO-
          ue = ue * (isoConeArea / etaBandArea);
          fEtaBandUECells->Fill(vecCOI.Pt() , ue);
          fEtIsoCells->Fill(isolation);
          if(isolation<eTThreshold)
          {
            fEtIsolatedCells->Fill(eTCOI);
            fEtisoT=eTCOI;
            fPtisoT=vecCOI.Pt();
          }
          break;
      }
        break;
        
      case 1: //EMCAL CLUSTERS
        
        switch (fUEMethod)
      {
        case 0: //phi band
          EtIsoClusPhiBand(vecCOI, isolation, ue,index);
          //Normalisation ue wrt Area - TO DO-
          ue = ue * (isoConeArea / phiBandArea);
          fPhiBandUEClust->Fill(vecCOI.Pt() , ue);
          fEtIsoClust->Fill(isolation);
          break;
        case 1: //eta band
          EtIsoClusEtaBand(vecCOI, isolation, ue,index);
          //Normalisation ue wrt Area - TO DO-
          ue = ue * (isoConeArea / etaBandArea);
          fEtaBandUEClust->Fill(vecCOI.Pt() , ue);
          fEtIsoClust->Fill(isolation);
          if(isolation<eTThreshold)
          {
            fEtIsolatedClust->Fill(eTCOI);
            fEtisoT=eTCOI;
            fPtisoT=vecCOI.Pt();
          }
          break;
      }
      case 2: //TRACKS (TBD which tracks) //EMCAL tracks
        switch (fUEMethod)
      {
        case 0: //phi band
          PtIsoTrackPhiBand(vecCOI, isolation, ue);
          //Normalisation ue wrt Area - TO DO-
          ue = ue * (isoConeArea / phiBandAreaTr);
          fPhiBandUETracks->Fill(vecCOI.Pt() , ue);
          fPtIsoTrack->Fill(isolation);
          if(isolation<eTThreshold)
          {
            fEtIsolatedTracks->Fill(eTCOI);
            fEtisoT=eTCOI;
            fPtisoT=vecCOI.Pt();
          }
          break;
        case 1: //eta band
          PtIsoTrackEtaBand(vecCOI, isolation, ue);
          //Normalisation ue wrt Area - TO DO-
          ue = ue * (isoConeArea / etaBandAreaTr);
          fEtaBandUETracks->Fill(vecCOI.Pt() , ue);
          fPtIsoTrack->Fill(isolation);
          if(isolation<eTThreshold)
          {
            fEtIsolatedTracks->Fill(eTCOI);
            fEtisoT=eTCOI;
            fPtisoT=vecCOI.Pt();
          }
          break;
          // case 2: //Cones
          // PtIsoTrackOrthCones(vecCOI, absId, isolation, ue);
          // break;
          // case 3: //Full TPC
          // PtIsoTrackFullTPC(vecCOI, absId, isolation, ue);
          // break;
      }
    }
  }
  else{  //EMCAL + TPC (Only tracks for the Isolation since IsoCone Goes Out of EMCAL)
    switch (fUEMethod)
    {
      case 0: //phi band
        PtIsoTrackPhiBand(vecCOI, isolation, ue);
        //Normalisation ue wrt Area - TO DO-
        ue = ue * (isoConeArea / phiBandAreaTr);
        fPhiBandUETracks->Fill(vecCOI.Pt() , ue);
        fPtIsoTrack->Fill(isolation);
        if(isolation<eTThreshold)
        {
          fEtIsolatedTracks->Fill(eTCOI);
          fEtisoT=eTCOI;
          fPtisoT=vecCOI.Pt();
        }
        break;
      case 1: //eta band
        PtIsoTrackEtaBand(vecCOI, isolation, ue);
        //Normalisation ue wrt Area - TO DO-
        ue = ue * (isoConeArea / etaBandAreaTr);
        fEtaBandUETracks->Fill(vecCOI.Pt() , ue);
        fPtIsoTrack->Fill(isolation);
        if(isolation<eTThreshold)
        {
          fEtIsolatedTracks->Fill(eTCOI);
          fEtisoT=eTCOI;
          fPtisoT=vecCOI.Pt();
        }
        break;
      case 2: //Cones
        PtIsoTrackOrthCones(vecCOI, isolation, ue);
        //Normalisation ue wrt Area - TO DO-
        ue = ue * (isoConeArea / perpConesArea);
        fPerpConesUETracks ->Fill(vecCOI.Pt() , ue);
        fPtIsoTrack->Fill(isolation);
        if(isolation<eTThreshold)
        {
          fEtIsolatedTracks->Fill(eTCOI);
          fEtisoT=eTCOI;
          fPtisoT=vecCOI.Pt();
        }
        break;
      case 3: //Full TPC
        PtIsoTrackFullTPC(vecCOI, isolation, ue);
        //Normalisation ue wrt Area - TO DO-
        ue = ue * (isoConeArea / fullTPCArea);
        fTPCWithoutIsoConeB2BbandUE->Fill(vecCOI.Pt() , ue);
        fPtIsoTrack->Fill(isolation);
        if(isolation<eTThreshold)
        {
          fEtIsolatedTracks->Fill(eTCOI);
          fEtisoT=eTCOI;
          fPtisoT=vecCOI.Pt();
        }
        break;
    }
  }
  
  
  //Here we should call something to know the number of tracks...
  //Soon I'll put in this version the "old way"; please let me know if
  //any of you could do the same with the JET framework
  
  switch(fWho) {
    case 0:
      flambda0T=m02COI;
      fEtT=vecCOI.Et();
      fPtT=vecCOI.Pt();
      fetaT=vecCOI.Eta();
      fphiT=vecCOI.Phi();
      fsumEtisoconeT=isolation;
      //	   AliError(Form("lambda 0 %f",flambda0T));
      fsumEtUE=ue;
      
      fOutputTree->Fill();
      break;
      
    case 1:
      outputValues[0] = nTracks;
      outputValues[1] = eTCOI;
      outputValues[2] = lambda0cluster;
      outputValues[3] = isolation;
      outputValues[4] = ue;
      outputValues[5] = isolation-ue;
      outputValues[6] = vecCOI.Eta();
      outputValues[7] = vecCOI.Phi();
      if (fIsMC) {
        outputValues[8] = ptmc;
        outputValues[9] = mcptsum;
      }
      fOutputTHnS -> Fill(outputValues);
      break;
      //                 //            fOutPTMAX -> Fill(outputValues[1],outputValues[2],);
  }
  return kTRUE;
}



