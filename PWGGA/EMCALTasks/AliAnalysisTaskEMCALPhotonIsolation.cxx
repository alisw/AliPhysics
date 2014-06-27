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
//Tracks(0),
//clusters(0)

{
  // Default constructor.
  
	//fParticleCollArray.SetOwner(kTRUE);
  for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
  
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
//Tracks(0),
//clusters(0)

{
  // Standard constructor.
  
	//fParticleCollArray.SetOwner(kTRUE);
  for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
  
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEMCALPhotonIsolation::~AliAnalysisTaskEMCALPhotonIsolation(){
  // Destructor
}


//________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::UserCreateOutputObjects(){
  
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

//     //________________________________________________________________________
Float_t* AliAnalysisTaskEMCALPhotonIsolation::GenerateFixedBinArray(Int_t n, Float_t min, Float_t max) const
{
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
    
    AliVCluster *COI = emccluster->GetCluster();
    if (!COI) return kTRUE;
    
    TLorentzVector VecCOI;
    COI->GetMomentum(VecCOI,fVertex);
    
    
    if(!CheckBoundaries(VecCOI))
      return kTRUE;
    
    if(ClustTrackMatching(COI))
      return kTRUE;
    
    else
      FillGeneralHistograms(COI,VecCOI, index);
  }
  else{
    //get the entries of the Cluster Container
    //whatever is a RETURN in LCAnalysis here is a CONTINUE,
    //since there are more than 1 Cluster per Event
    
    while((emccluster = static_cast<AliEmcalParticle*>(clusters->GetNextAcceptParticle()))){
      index++;
      if (!emccluster->IsEMCAL()) return kTRUE;
      
      
      AliVCluster *COI = emccluster->GetCluster();
      if(!COI) return kTRUE;
      
      TLorentzVector VecCOI;
      COI->GetMomentum(VecCOI,fVertex);
      
      if(!CheckBoundaries(VecCOI))
        return kTRUE;
      
      
      if(ClustTrackMatching(COI))
        continue;
      
      else
        FillGeneralHistograms(COI,VecCOI, index);
      
    }
  }
  //  PostData(1, fOutput);
  return kTRUE;
}


//     //__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::ClustTrackMatching(AliVCluster *cluster) {
  
  Double_t deta=999;
  Double_t dphi=999;
  
  AliParticleContainer *Tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
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
    Int_t im = ecl->GetTrackMatchedIndex();
    if(Tracks && im>=0) {
      mt = static_cast<AliVTrack*>(Tracks->GetParticle(im));
    }
  }
  //  if(mt && mt->TestFilterBit(768)) { //Hybrid Tracks with AOD
  AliPicoTrack::GetEtaPhiDiff(mt, cluster, dphi, deta);
  fDeltaETAClusTrackVSpT->Fill(nPart.Pt(), deta);
  fDeltaPHIClusTrackVSpT->Fill(nPart.Pt(), dphi);
  
  
  
  
  if(mt) return kTRUE;
  
  else return kFALSE;
  
}



//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoCellPhiBand(TLorentzVector c, Float_t &EtIso, Float_t &PhiBandcells){
  
  
  AliEMCALGeometry* EmcalGeom = AliEMCALGeometry::GetInstance();
  
  Float_t sumEnergyPhiBandCells=0., sumEnergyConeCells=0.;
  
  
  // check the cell corresponding to the leading cluster
  Int_t absId = 999;
  //maybe best to call it LeadingCellIdinClus or better maxId since it won't be used anywhere else ???
  Bool_t CellLeadingClustId = EmcalGeom->GetAbsCellIdFromEtaPhi(c.Eta(),c.Phi(),absId);
  if(!CellLeadingClustId) return;
  
  else{
    Int_t iTower = -1;
    Int_t iModule = -1;
    Int_t imEta=-1, imPhi=-1;
    Int_t ieta =-1, iphi =-1;
    
    EmcalGeom->GetCellIndex(absId,iModule,iTower,imPhi,imEta); // to get the module, the tower, eta and phi for the cell corresponding to the leading cluster
    EmcalGeom->GetCellPhiEtaIndexInSModule(iModule,iTower,imPhi,imEta,iphi,ieta); // to get the cell eta and phi in the super module for the cell coming from the leading cluster
    
    // Get the row and the column of the cell corresponding to the leading cluster in EMCAL
    Int_t colcellleadingclust = ieta;
    if(iModule % 2) colcellleadingclust = AliEMCALGeoParams::fgkEMCALCols + ieta ; // if the SM number is even you need to translate to have the one corresponding in EMCAL
    Int_t rowcellleadingclust = iphi + AliEMCALGeoParams::fgkEMCALRows*int(iModule/2); // to have the corresponding row in EMCAL
    
    // total number or rows and columns in EMCAL
    Int_t ntotalrows = AliEMCALGeoParams::fgkEMCALRows*16/3 ;  // 5 + 2/3 supermodules in a row
    Int_t ntotalcols = 2*AliEMCALGeoParams::fgkEMCALCols; // 2 supermodules in a column
    
    Int_t nbconesize = int(fIsoConeRadius/0.0143); //0.0143 tower acceptance
    
    // Get the cells
    AliVCaloCells * cells =InputEvent()->GetEMCALCells();
    
    // define the max and min row and column corresponding to the isolation cone around the seed cell from the leading cluster
    Int_t irowmincone = rowcellleadingclust-nbconesize;
    if(irowmincone<0) irowmincone=0;
    
    Int_t irowmaxcone = rowcellleadingclust+nbconesize;
    if(irowmaxcone>AliEMCALGeoParams::fgkEMCALRows) irowmaxcone=AliEMCALGeoParams::fgkEMCALRows;  // AliEMCALGeoParams::fgkEMCALRows = 24 in a supermodule
    
    Int_t icolmincone = colcellleadingclust - nbconesize;
    if(icolmincone<0) icolmincone = 0;
    
    Int_t icolmaxcone = colcellleadingclust+nbconesize;
    if(icolmaxcone>AliEMCALGeoParams::fgkEMCALCols) icolmaxcone=AliEMCALGeoParams::fgkEMCALCols;  // AliEMCALGeoParams::fgkEMCALCols = 48 in a supermodule
    
    // loop on all cells
    for(Int_t icol=0; icol<ntotalcols; icol++){
      for(Int_t irow=0; irow<ntotalrows; irow++){
        // now recover the cell indexes in a supermodule
        Int_t iSector = int(irow/AliEMCALGeoParams::fgkEMCALRows); // check in which SM is the cell
        if(iSector==5) continue;  //
        Int_t inModule = -1;
        Int_t icolLoc  = -1;
        if(icol < AliEMCALGeoParams::fgkEMCALCols){ // if the SM number is odd the column is the one corresponding in the supermodule
          iModule = 2*iSector + 1;
          icolLoc  = icol;
        }
        else if(icol > AliEMCALGeoParams::fgkEMCALCols - 1){ // if the SM number is even the column isn't the one corresponding in the supermodule
          inModule = 2*iSector;
          icolLoc  = icol-AliEMCALGeoParams::fgkEMCALCols;
        }
        Int_t irowLoc  = irow - AliEMCALGeoParams::fgkEMCALRows*iSector ;
        
        if(TMath::Abs(icol-colcellleadingclust)<nbconesize && TMath::Abs(icol+colcellleadingclust)>nbconesize){
          if(irow<irowmaxcone && irow>irowmincone){
            Int_t iabsId = EmcalGeom->GetAbsCellIdFromCellIndexes(iModule,irow,icol);  // verifier les irow et icol
            sumEnergyPhiBandCells+=cells->GetAmplitude(iabsId); //should be Et
          }
        }
        else if (TMath::Abs(icol-colcellleadingclust)>nbconesize && TMath::Abs(icol+colcellleadingclust)<nbconesize){
          Int_t iabsId = EmcalGeom->GetAbsCellIdFromCellIndexes(iModule,irowLoc,icolLoc);  // verifier les irow et icol
          sumEnergyConeCells+=cells->GetAmplitude(iabsId); //should be Et
        }
      }
    }
  }
  EtIso = sumEnergyConeCells;
  PhiBandcells = sumEnergyPhiBandCells;
}


//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoCellEtaBand(TLorentzVector c, Float_t &EtIso, Float_t &EtaBandcells){
  
  
  
  AliEMCALGeometry* EmcalGeom = AliEMCALGeometry::GetInstance();
  
  Float_t sumEnergyEtaBandCells=0., sumEnergyConeCells=0.;
  
  
  
  // check the cell corresponding to the leading cluster
  Int_t absId = 999;
  //maybe best to call it LeadingCellIdinClus or better maxId since it won't be used anywhere else ???
  Bool_t CellLeadingClustId = EmcalGeom->GetAbsCellIdFromEtaPhi(c.Eta(),c.Phi(),absId);
  if(!CellLeadingClustId) return;
  
  else{
    
    Int_t iTower = -1;
    Int_t iModule = -1;
    Int_t imEta=-1, imPhi=-1;
    Int_t ieta =-1, iphi =-1;
    
    EmcalGeom->GetCellIndex(absId,iModule,iTower,imPhi,imEta); // to get the module, the tower, eta and phi for the cell corresponding to the leading cluster
    EmcalGeom->GetCellPhiEtaIndexInSModule(iModule,iTower,imPhi,imEta,iphi,ieta); // to get the cell eta and phi in the super module for the cell coming from the leading cluster
    
    // Get the row and the column of the cell corresponding to the leading cluster in EMCAL
    Int_t colcellleadingclust = ieta;
    if(iModule % 2) colcellleadingclust = AliEMCALGeoParams::fgkEMCALCols + ieta ; // if the SM number is even you need to translate to have the one corresponding in EMCAL
    Int_t rowcellleadingclust = iphi + AliEMCALGeoParams::fgkEMCALRows*int(iModule/2); // to have the corresponding row in EMCAL
    
    // total number or rows and columns in EMCAL
    Int_t ntotalrows = AliEMCALGeoParams::fgkEMCALRows*16/3 ;  // 5 + 2/3 supermodules in a row
    Int_t ntotalcols = 2*AliEMCALGeoParams::fgkEMCALCols; // 2 supermodules in a column
    
    Int_t nbconesize = int(fIsoConeRadius/0.0143); //0.0143 tower acceptance
    
    // Get the cells
    AliVCaloCells * cells =InputEvent()->GetEMCALCells();
    
    // define the max and min row and column corresponding to the isolation cone around the seed cell from the leading cluster
    Int_t irowmincone = rowcellleadingclust-nbconesize;
    if(irowmincone<0) irowmincone=0;
    
    Int_t irowmaxcone = rowcellleadingclust+nbconesize;
    if(irowmaxcone>AliEMCALGeoParams::fgkEMCALRows) irowmaxcone=AliEMCALGeoParams::fgkEMCALRows;  // AliEMCALGeoParams::fgkEMCALRows = 24 in a supermodule
    
    Int_t icolmincone = colcellleadingclust-nbconesize;
    if(icolmincone<0) icolmincone = 0;
    
    Int_t icolmaxcone = colcellleadingclust+nbconesize;
    if(icolmaxcone>AliEMCALGeoParams::fgkEMCALCols) icolmaxcone=AliEMCALGeoParams::fgkEMCALCols;  // AliEMCALGeoParams::fgkEMCALCols = 48 in a supermodule
    
    // loop on all cells
    for(Int_t icol=0; icol<ntotalcols; icol++)
    {
      for(Int_t irow=0; irow<ntotalrows; irow++)
      {
        // now recover the cell indexes in a supermodule
        Int_t iSector = int(irow/AliEMCALGeoParams::fgkEMCALRows); // check in which SM is the cell
        if(iSector==5) continue;  //
        Int_t inModule = -1;
        Int_t icolLoc  = -1;
        if(icol < AliEMCALGeoParams::fgkEMCALCols){ // if the SM number is odd the column is the one corresponding in the supermodule
          iModule = 2*iSector + 1;
          icolLoc  = icol;
        }
        else if(icol > AliEMCALGeoParams::fgkEMCALCols - 1){ // if the SM number is even the column isn't the one corresponding in the supermodule
          inModule = 2*iSector;
          icolLoc  = icol-AliEMCALGeoParams::fgkEMCALCols;
        }
        
        Int_t irowLoc  = irow - AliEMCALGeoParams::fgkEMCALRows*iSector ;
        
        if(TMath::Abs(icol-colcellleadingclust)<nbconesize && TMath::Abs(icol+colcellleadingclust)>nbconesize){
          if(icol<icolmaxcone && icol>icolmincone){
            Int_t iabsId = EmcalGeom->GetAbsCellIdFromCellIndexes(iModule,irow,icol);  // verifier les irow et icol
            sumEnergyEtaBandCells+=cells->GetAmplitude(iabsId); //should be Et
          }
        }
        else if (TMath::Abs(icol-colcellleadingclust)>nbconesize && TMath::Abs(icol+colcellleadingclust)<nbconesize){
          Int_t iabsId = EmcalGeom->GetAbsCellIdFromCellIndexes(iModule,irowLoc,icolLoc);  // verifier les irow et icol
          sumEnergyConeCells+=cells->GetAmplitude(iabsId); //should be Et
        }
      }
    }
  }
  EtIso = sumEnergyConeCells;
  EtaBandcells = sumEnergyEtaBandCells;
}


//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoClusPhiBand(TLorentzVector c, Float_t &EtIso, Float_t &PhiBandclus, Int_t index){
  
  Float_t sumEnergyPhiBandClus=0., sumEnergyConeClus=0.;
  
  //needs a check on the same cluster
  AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliEmcalParticle *clust;
  
  Int_t localindex=-1;
  while((clust = static_cast<AliEmcalParticle*>(clusters->GetNextAcceptParticle()))){ //check the position of other clusters in respect to the leading cluster
    localindex++;
    if(localindex==index) continue;
    
    AliVCluster *cluster= clust->GetCluster();
    
    TLorentzVector nClust; //STILL NOT INITIALIZED
    cluster->GetMomentum(nClust,fVertex);
    Float_t phiclust =nClust.Phi();
    Float_t etaclust= nClust.Eta();
    //redefine phi/c.Eta() from the cluster we passed to the function
    
    Float_t  radius = TMath::Sqrt(TMath::Power(phiclust- c.Phi(),2)+TMath::Power(etaclust-c.Eta(),2)); // define the radius between the leading cluster and the considered cluster
    
    if(radius>fIsoConeRadius){ // the cluster is outside the isolation cone in this case study UE
      
      // actually phi band here
      if(TMath::Abs(etaclust - c.Eta()) < fIsoConeRadius){
        sumEnergyPhiBandClus += nClust.Pt();
      }
    }
    else // if the cluster is in the isolation cone, add the cluster pT
      sumEnergyConeClus += nClust.Pt();
    
  }
  EtIso = sumEnergyConeClus;
  PhiBandclus = sumEnergyPhiBandClus;
}


//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoClusEtaBand(TLorentzVector c, Float_t &EtIso, Float_t &EtaBandclus, Int_t index){
  
  Float_t sumEnergyEtaBandClus =0., sumEnergyConeClus=0.;
  AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  
  AliEmcalParticle *clust;
  
  Int_t localindex=-1;
  while((clust = static_cast<AliEmcalParticle*>(clusters->GetNextAcceptParticle()))){ //check the position of other clusters in respect to the leading cluster
    localindex++;
    if(localindex==index) continue;
    
    AliVCluster *cluster= clust->GetCluster();
    
    TLorentzVector nClust; //STILL NOT INITIALIZED
    cluster->GetMomentum(nClust,fVertex);
    
    Float_t phiclust =nClust.Phi();
    Float_t etaclust= nClust.Eta();
    //redefine phi/c.Eta() from the cluster we passed to the function
    
    // define the radius between the leading cluster and the considered cluster
    Float_t  radius = TMath::Sqrt(TMath::Power(phiclust-c.Phi(),2)+TMath::Power(etaclust-c.Eta(),2));
    
    if(radius>fIsoConeRadius){ // the cluster is outside the isolation cone in this case study UE
      
      // actually eta band here
      if(TMath::Abs(etaclust - c.Eta()) < fIsoConeRadius){
        sumEnergyEtaBandClus += nClust.Pt();
      }
    }
    else  // if the cluster is in the isolation cone, add the cluster pT
      sumEnergyConeClus += nClust.Pt();
    
  }
  EtIso = sumEnergyConeClus;
  EtaBandclus = sumEnergyEtaBandClus;
}


//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackPhiBand(TLorentzVector c, Float_t &PtIso, Float_t &PhiBandtrack){
  
  //INSERT BOUNDARIES ACCORDING TO THE FLAG WE WANT TO USE
  Float_t sumpTConeCharged=0., sumpTPhiBandTrack=0.;
  Float_t minphi= 0., maxphi= 2*TMath::Pi(), mineta = -0.9, maxeta= 0.9;
  
  AliParticleContainer *Tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  
  if(!fTPC4Iso){
    mineta = -0.7;
    maxeta = 0.7;
    minphi = 1.4;
    maxphi = TMath::Pi();
  }
  
  AliVParticle *track = Tracks->GetNextAcceptParticle(0);
  
  while(track){
    
    //CHECK IF TRACK IS IN BOUNDARIES
    Float_t phitrack = track->Phi();
    Float_t etatrack = track->Eta();
    // define the radius between the leading cluster and the considered cluster
    //redefine phi/c.Eta() from the cluster we passed to the function
    if(phitrack < maxphi && phitrack > minphi && etatrack < maxeta && etatrack > mineta){
      Float_t  radius = TMath::Sqrt(TMath::Power(phitrack - c.Phi(),2)+TMath::Power(etatrack - c.Eta(),2));
      
      if(radius>fIsoConeRadius){ // if tracks are outside the isolation cone study
        
        // actually phi band here --- ADD Boundaries conditions
        if(TMath::Abs(etatrack - c.Eta()) < fIsoConeRadius){
          sumpTPhiBandTrack += track->Pt();
        }
      }
      else
        sumpTConeCharged+=track->Pt(); // should not double count if the track matching is already done
    }
    track=Tracks->GetNextAcceptParticle();
  }
  PtIso = sumpTConeCharged;
  PhiBandtrack = sumpTPhiBandTrack;
}


//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackEtaBand(TLorentzVector c, Float_t &PtIso, Float_t &EtaBandtrack){
  
  //INSERT BOUNDARIES ACCORDING TO THE FLAG WE WANT TO USE
  Float_t sumpTConeCharged=0., sumpTEtaBandTrack=0.;
  Float_t minphi= 0., maxphi= 2*TMath::Pi(), mineta = -0.9, maxeta= 0.9;
  
  AliParticleContainer *Tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  
  
  if(!fTPC4Iso){
    mineta = -0.7;
    maxeta = 0.7;
    minphi = 1.4;
    maxphi = TMath::Pi();
  }
  
  
  AliVParticle *track = Tracks->GetNextAcceptParticle(0);
  
  while(track){
    
    Float_t phitrack = track->Phi();
    Float_t etatrack = track->Eta();
    //redefine phi/c.Eta() from the cluster we passed to the function
    if(phitrack < maxphi && phitrack > minphi && etatrack < maxeta && etatrack > mineta){
      Float_t  radius = TMath::Sqrt(TMath::Power(phitrack-c.Phi(),2)+TMath::Power(etatrack-c.Eta(),2)); // define the radius between the leading cluster and the considered cluster
      if(radius>fIsoConeRadius){ // if tracks are outside the isolation cone study UE
        
        // actually eta band here --- ADD Boundaries conditions
        if(TMath::Abs(phitrack - c.Phi()) < fIsoConeRadius){
          sumpTEtaBandTrack += track->Pt();
        }
      }
      else
        sumpTConeCharged += track->Pt(); // should not double count if the track matching is already done
    }
    track=Tracks->GetNextAcceptParticle();
  }
  PtIso = sumpTConeCharged;
  EtaBandtrack = sumpTEtaBandTrack;
}


//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackOrthCones(TLorentzVector c, Float_t &PtIso, Float_t &Cones){
  
  Float_t sumpTConeCharged=0., sumpTPerpConeTrack=0.;
  
  AliParticleContainer *Tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  
  Float_t EtaClus = c.Eta();
  Float_t PhiClus = c.Phi();
  Float_t PhiCone1 = PhiClus - TMath::PiOver2();
  Float_t PhiCone2 = PhiClus + TMath::PiOver2();
  
  if (PhiCone1 < 0.) PhiCone1 += 2*TMath::Pi();
  
  AliVParticle *track = Tracks->GetNextAcceptParticle(0);
  
  while(track){
    
    Float_t phitrack = track->Phi();
    Float_t etatrack = track->Eta();
    
    Float_t Dist2Clust = TMath::Sqrt(TMath::Power(etatrack-EtaClus, 2)+TMath::Power(phitrack-PhiClus, 2));
    if (Dist2Clust<fIsoConeRadius)// if tracks are inside the IsoCone
      sumpTConeCharged += track->Pt();
    
    else{//tracks outside the IsoCone
      //Distances from the centres of the two Orthogonal Cones
      Float_t Dist2Cone1 = TMath::Sqrt(TMath::Power(etatrack-EtaClus, 2)+TMath::Power(phitrack-PhiCone1, 2));
      Float_t Dist2Cone2 = TMath::Sqrt(TMath::Power(etatrack-EtaClus, 2)+TMath::Power(phitrack-PhiCone2, 2));
      
      //Is the Track Inside one of the two Cones ->Add to UE
      if((Dist2Cone1 < fIsoConeRadius) || (Dist2Cone2 < fIsoConeRadius))
        sumpTPerpConeTrack += track->Pt();
    }
    track=Tracks->GetNextAcceptParticle();
  }
  PtIso = sumpTConeCharged;
  Cones = sumpTPerpConeTrack;
}

//__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackFullTPC(TLorentzVector c, Float_t &PtIso, Float_t &Full){
  
  Float_t sumpTConeCharged=0., sumpTTPCexceptB2B=0.;
  
  AliParticleContainer *Tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  AliVParticle *track = Tracks->GetNextAcceptParticle(0);
  
  while(track){
    
    Float_t phitrack = track->Phi();
    Float_t etatrack = track->Eta();
    //redefine phi/c.Eta() from the cluster we passed to the function
    
    Float_t  radius = TMath::Sqrt(TMath::Power(phitrack-c.Phi(),2)+TMath::Power(etatrack-c.Eta(),2)); // define the radius between the leading cluster and the considered cluster
    if(radius>fIsoConeRadius){ // if tracks are outside the isolation cone study UE
      Double_t dphiup = c.Phi() + TMath::Pi() - fIsoConeRadius;
      Double_t dphidown = c.Phi() + TMath::Pi() + fIsoConeRadius;
      // TPC except B2B
      if(TMath::Power(phitrack-c.Phi(),2) +TMath::Power(etatrack-c.Eta(),2) > TMath::Power((fIsoConeRadius+0.1),2) && phitrack < dphidown && phitrack> dphiup)
        sumpTTPCexceptB2B += track->Pt();
    }
    else
      sumpTConeCharged += track->Pt(); // should not double count if the track matching is already done
    
    track=Tracks->GetNextAcceptParticle();
  }
  PtIso = sumpTConeCharged;
  Full = sumpTTPCexceptB2B;
}

//__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::CheckBoundaries(TLorentzVector VecCOI){
  
  Float_t minphibound= 1.5 , minetabound= -0.5, maxphibound= TMath::Pi()-0.1, maxetabound= 0.5;
  Bool_t isINBoundaries;
  
  if(!fTPC4Iso){
    minphibound = 1.4+0.4; //to be changed with fIsoConeR
    maxphibound = TMath::Pi()-0.4;
    minetabound = -0.7+0.4;
    maxetabound = 0.7-0.4;
  }
  if(VecCOI.Eta() > maxetabound || VecCOI.Eta() < minetabound || VecCOI.Phi() > maxphibound || VecCOI.Phi() >minphibound)
    isINBoundaries=kFALSE;
  else
    isINBoundaries=kTRUE;
  
  return isINBoundaries;
}

//__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::FillGeneralHistograms(AliVCluster *COI, TLorentzVector VecCOI, Int_t index){
  
  
  AliParticleContainer *Tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  
  AliEmcalParticle *emctrack = 0;
  
  int NTracks=0;
  Tracks->ResetCurrentID();
  while ((emctrack = static_cast<AliEmcalParticle*>(Tracks->GetNextAcceptParticle()))) {
    AliVTrack *track = emctrack->GetTrack();
    if(!track) continue;
    // if(!(track->TestFilterBit("kHybrid"))) continue;
    
    NTracks++;
  }
  
  Double_t ETCOI = 0., M02COI = 0., lambda0cluster = 0., phiCOI = 0., etaCOI = 0., ptmc = 0., mcptsum = 0.;
  //Int_t Ntracks;
  //Definition of the Array for Davide's Output
  const Int_t ndims =   fNDimensions;
  Double_t outputValues[ndims];
  
  ETCOI = VecCOI.Et();
  M02COI = COI->GetM02();
  etaCOI = VecCOI.Eta();
  phiCOI = VecCOI.Phi();
  
  // ******** Isolation and UE calculation with different methods *********
  
  Double_t EtThreshold = 5;
  
  switch(fEtIsoMethod)
  {
    case 0:  // SumEt<EtThr
      if(fM02mincut < M02COI && M02COI < fM02maxcut)  // photon candidate, cuts have to be decided after studies
      {
        EtThreshold = fEtIsoThreshold;
      }
      break;
      
    case 1:  // SumEt<%Ephoton
      if(fM02mincut < M02COI && M02COI < fM02maxcut) // photon candidate, cuts have to be decided after studies
      {
        EtThreshold = fEtIsoThreshold * ETCOI;
      }
      break;
      
    case 2: // Etmax<EtThreshold
      EtThreshold = fEtIsoThreshold;
      if(fM02mincut < M02COI && M02COI < fM02maxcut && ETCOI<EtThreshold) // photon candidate, cuts have to be decided after studies
      {
        fEtIsolatedClust->Fill(ETCOI);
      }
      break;
  }
  
  //DO NOT CHANGE EVER AGAIN THE FOLLOWING DEFINITIONS
  Float_t IsoConeArea = TMath::Pi()*fIsoConeRadius*fIsoConeRadius;
  Float_t etabandArea = (1.4 - 2.*fIsoConeRadius)*2.*fIsoConeRadius;
  Float_t phibandArea = ((5./9.)*TMath::Pi()- 2.*fIsoConeRadius)*2.*fIsoConeRadius;
  Float_t etabandAreaTr = (1.8 - 2.*fIsoConeRadius)*2.*fIsoConeRadius;
  Float_t phibandAreaTr = (2.*TMath::Pi() - 4.*fIsoConeRadius)*2.*fIsoConeRadius;//there is a 2 more because of the JET CONE B2B
  Float_t PerpConesArea = 2.*IsoConeArea;
  Float_t FullTPCArea = 1.8*2.*TMath::Pi()-fIsoConeRadius*(TMath::Pi()*fIsoConeRadius + 3.6);
  
  Float_t ISOLATION, UE;
  
  if(!fTPC4Iso){ //EMCAL Only for Acceptance of Clusters
    switch(fIsoMethod)
    {
      case 0: //EMCAL CELLS
        
        switch (fUEMethod)
      {
        case 0: //phi band
          EtIsoCellPhiBand(VecCOI, ISOLATION, UE);
          //Normalisation UE wrt Area - TO DO-
          UE = UE * (IsoConeArea / phibandArea);
          fPhiBandUECells->Fill(VecCOI.Pt() , UE);
          fEtIsoCells->Fill(ISOLATION);
          if(ISOLATION<EtThreshold)
          {
            fEtIsolatedCells->Fill(ETCOI);
            fEtisoT=ETCOI;
            fPtisoT=VecCOI.Pt();
          }
          break;
        case 1: //eta band
          EtIsoCellEtaBand(VecCOI, ISOLATION, UE);
          //Normalisation UE wrt Area - TO DO-
          UE = UE * (IsoConeArea / etabandArea);
          fEtaBandUECells->Fill(VecCOI.Pt() , UE);
          fEtIsoCells->Fill(ISOLATION);
          if(ISOLATION<EtThreshold)
          {
            fEtIsolatedCells->Fill(ETCOI);
            fEtisoT=ETCOI;
            fPtisoT=VecCOI.Pt();
          }
          break;
      }
        break;
        
      case 1: //EMCAL CLUSTERS
        
        switch (fUEMethod)
      {
        case 0: //phi band
          EtIsoClusPhiBand(VecCOI, ISOLATION, UE,index);
          //Normalisation UE wrt Area - TO DO-
          UE = UE * (IsoConeArea / phibandArea);
          fPhiBandUEClust->Fill(VecCOI.Pt() , UE);
          fEtIsoClust->Fill(ISOLATION);
          break;
        case 1: //eta band
          EtIsoClusEtaBand(VecCOI, ISOLATION, UE,index);
          //Normalisation UE wrt Area - TO DO-
          UE = UE * (IsoConeArea / etabandArea);
          fEtaBandUEClust->Fill(VecCOI.Pt() , UE);
          fEtIsoClust->Fill(ISOLATION);
          if(ISOLATION<EtThreshold)
          {
            fEtIsolatedClust->Fill(ETCOI);
            fEtisoT=ETCOI;
            fPtisoT=VecCOI.Pt();
          }
          break;
      }
      case 2: //TRACKS (TBD which tracks) //EMCAL Tracks
        switch (fUEMethod)
      {
        case 0: //phi band
          PtIsoTrackPhiBand(VecCOI, ISOLATION, UE);
          //Normalisation UE wrt Area - TO DO-
          UE = UE * (IsoConeArea / phibandAreaTr);
          fPhiBandUETracks->Fill(VecCOI.Pt() , UE);
          fPtIsoTrack->Fill(ISOLATION);
          if(ISOLATION<EtThreshold)
          {
            fEtIsolatedTracks->Fill(ETCOI);
            fEtisoT=ETCOI;
            fPtisoT=VecCOI.Pt();
          }
          break;
        case 1: //eta band
          PtIsoTrackEtaBand(VecCOI, ISOLATION, UE);
          //Normalisation UE wrt Area - TO DO-
          UE = UE * (IsoConeArea / etabandAreaTr);
          fEtaBandUETracks->Fill(VecCOI.Pt() , UE);
          fPtIsoTrack->Fill(ISOLATION);
          if(ISOLATION<EtThreshold)
          {
            fEtIsolatedTracks->Fill(ETCOI);
            fEtisoT=ETCOI;
            fPtisoT=VecCOI.Pt();
          }
          break;
          // case 2: //Cones
          // PtIsoTrackOrthCones(VecCOI, absId, ISOLATION, UE);
          // break;
          // case 3: //Full TPC
          // PtIsoTrackFullTPC(VecCOI, absId, ISOLATION, UE);
          // break;
      }
    }
  }
  else{  //EMCAL + TPC (Only Tracks for the Isolation since IsoCone Goes Out of EMCAL)
    switch (fUEMethod)
    {
      case 0: //phi band
        PtIsoTrackPhiBand(VecCOI, ISOLATION, UE);
        //Normalisation UE wrt Area - TO DO-
        UE = UE * (IsoConeArea / phibandAreaTr);
        fPhiBandUETracks->Fill(VecCOI.Pt() , UE);
        fPtIsoTrack->Fill(ISOLATION);
        if(ISOLATION<EtThreshold)
        {
          fEtIsolatedTracks->Fill(ETCOI);
          fEtisoT=ETCOI;
          fPtisoT=VecCOI.Pt();
        }
        break;
      case 1: //eta band
        PtIsoTrackEtaBand(VecCOI, ISOLATION, UE);
        //Normalisation UE wrt Area - TO DO-
        UE = UE * (IsoConeArea / etabandAreaTr);
        fEtaBandUETracks->Fill(VecCOI.Pt() , UE);
        fPtIsoTrack->Fill(ISOLATION);
        if(ISOLATION<EtThreshold)
        {
          fEtIsolatedTracks->Fill(ETCOI);
          fEtisoT=ETCOI;
          fPtisoT=VecCOI.Pt();
        }
        break;
      case 2: //Cones
        PtIsoTrackOrthCones(VecCOI, ISOLATION, UE);
        //Normalisation UE wrt Area - TO DO-
        UE = UE * (IsoConeArea / PerpConesArea);
        fPerpConesUETracks ->Fill(VecCOI.Pt() , UE);
        fPtIsoTrack->Fill(ISOLATION);
        if(ISOLATION<EtThreshold)
        {
          fEtIsolatedTracks->Fill(ETCOI);
          fEtisoT=ETCOI;
          fPtisoT=VecCOI.Pt();
        }
        break;
      case 3: //Full TPC
        PtIsoTrackFullTPC(VecCOI, ISOLATION, UE);
        //Normalisation UE wrt Area - TO DO-
        UE = UE * (IsoConeArea / FullTPCArea);
        fTPCWithoutIsoConeB2BbandUE->Fill(VecCOI.Pt() , UE);
        fPtIsoTrack->Fill(ISOLATION);
        if(ISOLATION<EtThreshold)
        {
          fEtIsolatedTracks->Fill(ETCOI);
          fEtisoT=ETCOI;
          fPtisoT=VecCOI.Pt();
        }
        break;
    }
  }
  
  
  //Here we should call something to know the number of tracks...
  //Soon I'll put in this version the "old way"; please let me know if
  //any of you could do the same with the JET framework
  
  switch(fWho) {
    case 0:
      flambda0T=M02COI;
      fEtT=VecCOI.Et();
      fPtT=VecCOI.Pt();
      fetaT=VecCOI.Eta();
      fphiT=VecCOI.Phi();
      fsumEtisoconeT=ISOLATION;
      //	   AliError(Form("lambda 0 %f",flambda0T));
      fsumEtUE=UE;
      
      fOutputTree->Fill();
      break;
      
    case 1:
      outputValues[0] = NTracks;
      outputValues[1] = ETCOI;
      outputValues[2] = lambda0cluster;
      outputValues[3] = ISOLATION;
      outputValues[4] = UE;
      outputValues[5] = ISOLATION-UE;
      outputValues[6] = VecCOI.Eta();
      outputValues[7] = VecCOI.Phi();
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



