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
#include <TH3D.h>
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
#include "AliTrackContainer.h"
#include "AliVTrack.h"
#include "AliEmcalParticle.h"
#include "AliParticleContainer.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVCaloCells.h"
#include "AliPicoTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEMCALRecoUtils.h"
#include "AliLog.h"
#include "TF1.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "AliGenPythiaEventHeader.h"


#include "AliAnalysisTaskEMCALPhotonIsolation.h"

  /// \cond CLASSIMP
ClassImp(AliAnalysisTaskEMCALPhotonIsolation);
  /// \endcond

using std::cout;
using std::endl;
  //________________________________________________________________________
AliAnalysisTaskEMCALPhotonIsolation::AliAnalysisTaskEMCALPhotonIsolation() :
AliAnalysisTaskEmcal("AliAnalysisTaskEMCALPhotonIsolation",kTRUE),
  //fParticleCollArray(),
fAOD(0),
fVevent(0),
fNCluster(0),
fAODMCParticles(0),
fTracksAna(0),
fStack(0),
fEMCALRecoUtils(new AliEMCALRecoUtils),
fWho(-1),
fSSsmearing(0),
fSSsmearwidth(0),
fSSsmear_mean(0),
fWhich(0),
  //fOutputList(0),
fIsoConeRadius(0.4),
fEtIsoMethod(0),
fEtIsoThreshold(2),
fdetacut(0.025),
fdphicut(0.03),
fdetacutIso(0.025),
fdphicutIso(0.03),
fM02mincut(0.1),
fM02maxcut(0.3),
fExtraIsoCuts(kFALSE),
fQA(0),
fIsMC(0),
fTPC4Iso(0),
fIsoMethod(0),
fUEMethod(0),
fNDimensions(0),
fMCDimensions(0),
fMCQAdim(0),
fisLCAnalysis(0),
fIsNLMCut(kFALSE),
fNLMCut(0),
fNLMmin(0),
fTMClusterRejected(kTRUE),
fTMClusterInConeRejected(kTRUE),
fRejectionEventWithoutTracks(kFALSE),
fAnalysispPb(kFALSE),
fTriggerLevel1(0),
fTest1(0),
fTest2(0),
fEClustersT(0),
fPtClustersT(0),
fEtClustersT(0),
fEtaClustersT(0),
fPhiClustersT(0),
fM02ClustersT(0),
fevents(0),
fNClustersT(0),
flambda0T(0),
fM02isoT(0),
fM02noisoT(0),
fPtnoisoT(0),
fEtT(0),
fPtT(0),
fPtisoT(0),
fEtisolatedT(0),
fPtisolatedT(0),
fetaT(0),
fphiT(0),
fsumEtisoconeT(0),
fsumEtUE(0),
fBinsPt(),
fBinsM02(),
fBinsEtiso(),
fBinsEtue(),
fBinsEta(),
fBinsPhi(),
fBinsLabel(),
fBinsPDG(),
fBinsMomPDG(),
fBinsClustPDG(),
fBinsDx(),
fBinsDz(),
fBinsDecay(),
fTrackMult(0),
fEtaPhiClus(0),
fClusEvsClusT(0),
fPT(0),
fE(0),
fNLM(0),
fVz(0),
fEvents(0),
fPtaftTime(0),
fPtaftCell(0),
fPtaftNLM(0),
fPtaftTM(0),
fPtaftDTBC(0),
fPtaftFC(0),
fPtaftM02C(0),
fClusTime(0),
fM02(0),
fDeltaETAClusTrack(0),
fDeltaPHIClusTrack(0),
fDeltaETAClusTrackMatch(0),
fDeltaPHIClusTrackMatch(0),
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
fEtIsolatedCells(0),
fEtIsolatedClust(0),
fPtIsolatedNClust(0),
fPtIsolatedNTracks(0),
fEtIsolatedTracks(0),
fPtvsM02iso(0),
fPtvsM02noiso(0),
fTestIndex(0),
fTestIndexE(0),
fTestLocalIndexE(0),
fTestEnergyCone(0),
fTestEtaPhiCone(0),
fInvMassM02iso(0),
fInvMassM02noiso(0),
fPtvsM02vsSumPi0(0),
fPtvsM02vsSumEta(0),
fPtvsM02vsSum(0),
fPtvsM02vsSumUE (0),
fTrackMultvsSumChargedvsUE(0),
fTrackMultvsPt(0),
fTracksConeEtaPt(0),
fTracksConeEtaM02(0),
fHistXsection(0),
fHistTrials(0),
fPtTracksVSpTNC(0),
fCTdistVSpTNC(0),
fPtTracksVSpTNC_MC(0),
fpi0VSclusterVSIsolation(0),
fpi0VSclusterVSM02(0),
fpi0VSM02VSIsolation(0),
fPhiTracksVSclustPt(0),
fEtaTracksVSclustPt(0),
fOutputTHnS(0),
fOutMCTruth(0),
fOutClustMC(0),
fOutputQATree(0),
fOutputTree(0),
fphietaPhotons(0),
fphietaOthers(0),
fphietaOthersBis(0)
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
fAOD(0),
fVevent(0),
fNCluster(0),
fAODMCParticles(0),
fTracksAna(0),
fStack(0),
fEMCALRecoUtils(new AliEMCALRecoUtils),
fWho(-1),
fSSsmearing(0),
fSSsmearwidth(0),
fSSsmear_mean(0),
fWhich(0),
  //fOutputList(0),
fIsoConeRadius(0.4),
fEtIsoMethod(0),
fEtIsoThreshold(2),
fdetacut(0.025),
fdphicut(0.03),
fdetacutIso(0.025),
fdphicutIso(0.03),
fM02mincut(0.1),
fM02maxcut(0.3),
fExtraIsoCuts(kFALSE),
fQA(0),
fIsMC(0),
fTPC4Iso(0),
fIsoMethod(0),
fUEMethod(0),
fNDimensions(0),
fMCDimensions(0),
fMCQAdim(0),
fisLCAnalysis(0),
fIsNLMCut(kFALSE),
fNLMCut(0),
fNLMmin(0),
fTMClusterRejected(kTRUE),
fTMClusterInConeRejected(kTRUE),
fRejectionEventWithoutTracks(kFALSE),
fAnalysispPb(kFALSE),
fTriggerLevel1(0),
fTest1(0),
fTest2(0),
fEClustersT(0),
fPtClustersT(0),
fEtClustersT(0),
fEtaClustersT(0),
fPhiClustersT(0),
fM02ClustersT(0),
fevents(0),
fNClustersT(0),
flambda0T(0),
fM02isoT(0),
fM02noisoT(0),
fPtnoisoT(0),
fEtT(0),
fPtT(0),
fPtisoT(0),
fEtisolatedT(0),
fPtisolatedT(0),
fetaT(0),
fphiT(0),
fsumEtisoconeT(0),
fsumEtUE(0),
fBinsPt(),
fBinsM02(),
fBinsEtiso(),
fBinsEtue(),
fBinsEta(),
fBinsPhi(),
fBinsLabel(),
fBinsPDG(),
fBinsMomPDG(),
fBinsClustPDG(),
fBinsDx(),
fBinsDz(),
fBinsDecay(),
fTrackMult(0),
fEtaPhiClus(0),
fClusEvsClusT(0),
fPT(0),
fE(0),
fNLM(0),
fVz(0),
fEvents(0),
fPtaftTime(0),
fPtaftCell(0),
fPtaftNLM(0),
fPtaftTM(0),
fPtaftDTBC(0),
fPtaftFC(0),
fPtaftM02C(0),
fClusTime(0),
fM02(0),
fDeltaETAClusTrack(0),
fDeltaPHIClusTrack(0),
fDeltaETAClusTrackMatch(0),
fDeltaPHIClusTrackMatch(0),
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
fEtIsolatedCells(0),
fEtIsolatedClust(0),
fPtIsolatedNClust(0),
fPtIsolatedNTracks(0),
fEtIsolatedTracks(0),
fPtvsM02iso(0),
fPtvsM02noiso(0),
fTestIndex(0),
fTestIndexE(0),
fTestLocalIndexE(0),
fTestEnergyCone(0),
fTestEtaPhiCone(0),
fInvMassM02iso(0),
fInvMassM02noiso(0),
fPtvsM02vsSumPi0(0),
fPtvsM02vsSumEta(0),
fPtvsM02vsSum(0),
fPtvsM02vsSumUE (0),
fTrackMultvsSumChargedvsUE(0),
fTrackMultvsPt(0),
fTracksConeEtaPt(0),
fTracksConeEtaM02(0),
fHistXsection(0),
fHistTrials(0),
fPtTracksVSpTNC(0),
fCTdistVSpTNC(0),
fPtTracksVSpTNC_MC(0),
fpi0VSclusterVSIsolation(0),
fpi0VSclusterVSM02(0),
fpi0VSM02VSIsolation(0),
fPhiTracksVSclustPt(0),
fEtaTracksVSclustPt(0),
fOutputTHnS(0),
fOutMCTruth(0),
fOutClustMC(0),
fOutputQATree(0),
fOutputTree(0),
fphietaPhotons(0),
fphietaOthers(0),
fphietaOthersBis(0)
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
    //printf("Up here all good");
  
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
  
  if(fUEMethod>1 && !fTPC4Iso){
    cout<<"Please set UE Method Accordingly to the Use of the TPC for the Analysis"<<endl;
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
  
  if(fWho>2 || fWho==-1){
    cout<<"Error!!! OutputMode Can Only Be 0: TTree; 1: THnSparse"<<endl;
    return;
  }
  else{
    fOutput = new AliEmcalList(); // RH: Leak? fOutput already exists in base class
    fOutput->SetOwner();
      //Initialize the common Output histograms
    switch (fWho)
    {
      case 0:
      {
          // Tree for QA after cluster selection
        fOutputQATree = new TTree("OutQATree","OutQATree");
        fOutputQATree->Branch("fevents",&fevents);
        fOutputQATree->Branch("fNClustersT",&fNClustersT);
        fOutputQATree->Branch("fEClustersT",&fEClustersT);
        fOutputQATree->Branch("fPtClustersT",&fPtClustersT);
        fOutputQATree->Branch("fEtClustersT",&fEtClustersT);
        fOutputQATree->Branch("fEtaClustersT",&fEtaClustersT);
        fOutputQATree->Branch("fPhiClustersT",&fPhiClustersT);
        fOutputQATree->Branch("fM02ClustersT",&fM02ClustersT);
        
        fOutput->Add(fOutputQATree);
        
        fOutputTree = new TTree("OutTree",Form("OutTree; Iso Method %d, UE Method %d, TPC %d, LC %d, Iso Cone %f, CPV #eta %f #phi %f",fIsoMethod,fUEMethod,fTPC4Iso,fisLCAnalysis,fIsoConeRadius,fdetacut,fdphicut));
        fOutputTree->Branch("flambda0T",&flambda0T);
        fOutputTree->Branch("fEtT",&fEtT);
        fOutputTree->Branch("fPtT",&fPtT);
        fOutputTree->Branch("fEtisolatedT",&fEtisolatedT);
        fOutputTree->Branch("fPtTiso",&fPtisolatedT);
        fOutputTree->Branch("fetaT",&fetaT);
        fOutputTree->Branch("fphiT",&fphiT);
        fOutputTree->Branch("fsumEtisoconeT",&fsumEtisoconeT);
        fOutputTree->Branch("fsumEtUE",&fsumEtUE);
        
        fOutput->Add(fOutputTree);
      }
        
        break;
      case 1:
          //Initialization by Davide;
      {
        TString sTitle;
        


	Int_t binPT = fBinsPt.size()-1;
	Int_t binM02 = fBinsM02.size()-1;
	Int_t binETiso = fBinsEtiso.size()-1;
	Int_t binETUE = fBinsEtue.size()-1;
	Int_t binetacl = fBinsEta.size()-1;
	Int_t binphicl = fBinsPhi.size()-1;
	Int_t binlabel = fBinsLabel.size()-1;
        
        Int_t binMCPDG = fBinsPDG.size()-1;
        Int_t binMCMotherPDG = fBinsMomPDG.size()-1;
        Int_t binMCClustPDG = fBinsClustPDG.size()-1;
	Int_t bindx = fBinsDx.size()-1;
	Int_t bindz = fBinsDz.size()-1 ;
	Int_t binDecayType = fBinsDecay.size()-1;
	//bincells=20;
        
        Int_t bins[] = {binPT, binM02, binETiso, binETUE, binetacl, binphicl};


	        fNDimensions = sizeof(bins)/sizeof(Int_t);
        const Int_t ndims =   fNDimensions;
        
        sTitle = Form("Direct Photons: p_{T} , M02 , E_{T} Iso%s in %s, E_{T} UE %s in %s, #eta_{clus} distr,#phi_{clus} distr; p_{T} (GeV/c); M02; E_{T}^{iso%s} (GeV/c) ; E_{T}^{UE%s} (GeV/c); #eta_{cl}; #phi_{cl}", sIsoMethod.Data(), sBoundaries.Data(), sUEMethod.Data(), sBoundaries.Data(), sIsoMethod.Data(), sUEMethod.Data());
        
        fOutputTHnS =  new THnSparseF("fHnOutput",sTitle.Data(), ndims, bins);
        fOutputTHnS->SetBinEdges(0,fBinsPt.data());
        fOutputTHnS->SetBinEdges(1,fBinsM02.data());
        fOutputTHnS->SetBinEdges(2,fBinsEtiso.data());
        fOutputTHnS->SetBinEdges(3,fBinsEtue.data());
        fOutputTHnS->SetBinEdges(4,fBinsEta.data());
        fOutputTHnS->SetBinEdges(5,fBinsPhi.data());
        fOutputTHnS->Sumw2();
        
        fOutput->Add(fOutputTHnS);
        
        
        if(fIsMC){
          Int_t binsMC[] = {binPT, binETiso, binETUE, binMCPDG ,binetacl,binphicl,binlabel};
          Int_t binsSMC[] = {binPT, binM02, binMCClustPDG, binMCMotherPDG, binPT, bindx, bindz, binETiso,binDecayType};
          
          fMCDimensions = sizeof(binsMC)/sizeof(Int_t);
	  const Int_t ndimsMC = fMCDimensions;
          
          //Double_t xminbis[] = {0., -10., -10., -1000., -1.0,  1.,    0};
          //Double_t xmaxbis[] = {70., 100., 100.,  1000.,  1.0, 3.5, 1500};
          
          //fOutMCTruth = new THnSparseF ("fOutMCTruth","E_{#gamma}, E_{T}^{iso cone}, E_{T}^{UE}, MomPDG, Eta, Phi, Label; E_{T}^{#gamma} (GeV/c); p_{T}^{Iso}(GeV/c);E_{T} ^{UE} (GeV/c); PDG; #eta; #phi; Label",7,binsMC,xminbis,xmaxbis);
          //fOutMCTruth->Sumw2();
          //fOutput->Add(fOutMCTruth);
	  //
          
          fOutMCTruth = new THnSparseF ("fOutMCTruth","E_{#gamma}, E_{T}^{iso cone}, E_{T}^{UE}, MomPDG, Eta, Phi, Label; E_{T}^{#gamma} (GeV/c); p_{T}^{Iso}(GeV/c);E_{T} ^{UE} (GeV/c); PDG; #eta; #phi; Label",ndimsMC,binsMC);
	  fOutMCTruth->SetBinEdges(0,fBinsPt.data());
          fOutMCTruth->SetBinEdges(1,fBinsEtiso.data());
          fOutMCTruth->SetBinEdges(2,fBinsEtue.data());
          fOutMCTruth->SetBinEdges(3,fBinsPDG.data());
          fOutMCTruth->SetBinEdges(4,fBinsEta.data());
          fOutMCTruth->SetBinEdges(5,fBinsPhi.data());
          fOutMCTruth->SetBinEdges(6,fBinsLabel.data());
          fOutMCTruth->Sumw2();
          fOutput->Add(fOutMCTruth);
          
	  
	  fphietaPhotons = new TH3D ("fDphiDeta_Photons","#Delta#phi vs #Delta#eta Clust-MCpart to check why zero M02 clusters; #eta; #phi", 100, -0.5, 0.5, 200, 1.5, 3.5,60,0.,60.);
          fphietaPhotons->Sumw2();
          fOutput->Add(fphietaPhotons);
          
          fphietaOthers = new TH3D ("fDphiDeta_Others","#Delta#phi vs #Delta#eta Clust-MCpart to check why zero M02 clusters; #eta; #phi", 140, -0.7, 0.7, 220, 0.8, 3.5,60,0.,60.);
          fphietaOthers->Sumw2();
          fOutput->Add(fphietaOthers);
          
          fphietaOthersBis = new TH3D ("fDphiDeta_OthersBis","#Delta#phi vs #Delta#eta Clust-MCpart to check why zero M02 clusters; #eta; #phi", 140, -0.7, 0.7, 220, 0.8, 3.5,60,0.,60.);
          fphietaOthersBis->Sumw2();
          fOutput->Add(fphietaOthersBis);
          
          fMCQAdim = sizeof(binsSMC)/sizeof(Int_t);
	  const Int_t ndimsMCQA = fMCQAdim;

          Double_t xminbismix[] = {0.,  0., -3000, -400,  0.,-1., -1., -10,    0.};
          Double_t xmaxbismix[] = {70., 2.,  3000,  400, 70., 1.,  1., 100.,  10.};
          
          //fOutClustMC = new THnSparseF ("fOutClustMC", "E_{T}^{clust}, M02, PDG, MOM PDG, E_{T}^{true}, #Deltax, #Deltaz, E_{T}^{iso},Label;E_{T}^{reco} (GeV/c); M02;PDG Code; Mothers' PDG Code; E_{T}^{MCtrue} (GeV/c); #Delta#phi; #Delta#eta; E_{T}^{iso} (Gev/c);Label",9,binsSMC,xminbismix,xmaxbismix);
          fOutClustMC = new THnSparseF ("fOutClustMC", "E_{T}^{clust}, M02, PDG, MOM PDG, E_{T}^{true}, #Deltax, #Deltaz, E_{T}^{iso},Label;E_{T}^{reco} (GeV/c); M02;PDG Code; Mothers' PDG Code; E_{T}^{MCtrue} (GeV/c); #Delta#phi; #Delta#eta; E_{T}^{iso} (Gev/c);Label",ndimsMCQA,binsSMC);
       	  fOutClustMC->SetBinEdges(0,fBinsPt.data());
          fOutClustMC->SetBinEdges(1,fBinsM02.data());
          fOutClustMC->SetBinEdges(2,fBinsClustPDG.data());
          fOutClustMC->SetBinEdges(3,fBinsMomPDG.data());
          fOutClustMC->SetBinEdges(4,fBinsPt.data());
          fOutClustMC->SetBinEdges(5,fBinsDx.data());
          fOutClustMC->SetBinEdges(6,fBinsDz.data());
          fOutClustMC->SetBinEdges(7,fBinsEtiso.data());
          fOutClustMC->SetBinEdges(8,fBinsDecay.data());
	  fOutClustMC->Sumw2();
          fOutput->Add(fOutClustMC);
        }
      }
        break;
      case 2:
      {
          //initialization Lucile's output
        fPtaftM02C = new TH1D("hPtaftM02C_NC","p_{T} distribution for Clusters after shower shape cut",200,0.,100.);
        fPtaftM02C->Sumw2();
        fOutput->Add(fPtaftM02C);
        
        fM02 = new TH2D("hM02_NC","M02 distribution for Neutral Clusters vs E",100,0.,100.,500,0.,5.);
        fM02->Sumw2();
        fOutput->Add(fM02);
        
        fEtIsoCells = new TH1D("hEtIsoCell_NC","E_{T}^{iso cone} in iso cone distribution for Neutral Clusters with EMCAL Cells",200,-0.25,99.75);
        fEtIsoCells->SetXTitle("#Sigma E_{T}^{iso cone} (GeV/c)");
        fEtIsoCells->Sumw2();
        fOutput->Add(fEtIsoCells);
        
        fEtIsoClust = new TH2D("hEtIsoClus_NC","#Sigma p_{T}^{iso cone} in iso cone distribution for Neutral Clusters with EMCAL Clusters",200,0.,100.,200,0.,100.);
        fEtIsoClust->SetYTitle("#Sigma P_{T}^{iso cone} (GeV/c)");
        fEtIsoClust->SetXTitle("p_{T}^{clust}");
        fEtIsoClust->Sumw2();
        fOutput->Add(fEtIsoClust);
        
        fPtIsoTrack = new TH2D("hPtIsoTrack_NC"," #Sigma p_{T}^{iso cone} in iso cone distribution for Neutral Clusters with Tracks",200,0.,100.,200,0.,100.);
        fPtIsoTrack->SetYTitle("#Sigma p_{T}^{iso cone} (GeV/c)");
        fPtIsoTrack->SetXTitle("p_{T}^{clust}");
        fPtIsoTrack->Sumw2();
        fOutput->Add(fPtIsoTrack);
        
        fPtEtIsoTC = new TH1D("hPtEtIsoTrackClust_NC","#Sigma P_{T}^{iso cone} + #Sigma E_{T}^{iso cone} in iso cone distribution for Neutral Clusters with Tracks and Clusters",200,-0.25,99.75);
        fPtEtIsoTC->SetXTitle("#Sigma P_{T}^{iso cone} + #Sigma E_{T}^{iso cone} (GeV/c)");
        fPtEtIsoTC->Sumw2();
        fOutput->Add(fPtEtIsoTC);
        
        fPhiBandUEClust = new TH2D(Form("hPhiBandUE_Cluster"),Form("UE Estimation with Phi Band Clusters"),200,0.,100.,250,0.,100.);
        fPhiBandUEClust->SetXTitle("E_{T}");
        fPhiBandUEClust->SetYTitle("#Sigma E_{T}^{UE}");
        fPhiBandUEClust->Sumw2();
        fOutput->Add(fPhiBandUEClust);
        
        fEtaBandUEClust = new TH2D(Form("hEtaBandUE_Cluster"),Form("UE Estimation with Eta Band Clusters"),200,0.,100.,250,0.,100.);
        fEtaBandUEClust->SetXTitle("E_{T}");
        fEtaBandUEClust->SetYTitle("#Sigma E_{T}^{UE}");
        fEtaBandUEClust->Sumw2();
        fOutput->Add(fEtaBandUEClust);
        
        fPhiBandUECells = new TH2D(Form("hPhiBandUE_CELLS"),Form("UE Estimation with Phi Band CELLS"),200,0.,100.,250,0.,100.);
        fPhiBandUECells->SetXTitle("E_{T}");
        fPhiBandUECells->SetYTitle("#Sigma E_{T}^{UE}");
        fPhiBandUECells->Sumw2();
        fOutput->Add(fPhiBandUECells);
        
        fEtaBandUECells = new TH2D(Form("hEtaBandUE_CELLS"),Form("UE Estimation with Eta Band CELLS"),200,0.,100.,250,0.,100.);
        fEtaBandUECells->SetXTitle("E_{T}");
        fEtaBandUECells->SetYTitle("#Sigma E_{T}^{UE}");
        fEtaBandUECells->Sumw2();
        fOutput->Add(fEtaBandUECells);
        
        fPhiBandUETracks = new TH2D(Form("hPhiBandUE_TPC"),Form("UE Estimation with Phi Band TPC "),200,0.,100.,250,0.,100.);
        fPhiBandUETracks->SetXTitle("E_{T}");
        fPhiBandUETracks->SetYTitle("#Sigma P_{T}^{UE}");
        fPhiBandUETracks->Sumw2();
        fOutput->Add(fPhiBandUETracks);
        
        fEtaBandUETracks = new TH2D(Form("hEtaBandUE_TPC"),Form("UE Estimation with Eta Band and TPC"),200,0.,100.,250,0.,100.);
        fEtaBandUETracks->SetXTitle("E_{T}");
        fEtaBandUETracks->SetYTitle("#Sigma P_{T}^{UE}");
        fEtaBandUETracks->Sumw2();
        fOutput->Add(fEtaBandUETracks);
        
        fPerpConesUETracks = new TH2D("hConesUE","UE Estimation with Perpendicular Cones in TPC",200,0.,100.,250,0.,100.);
        fPerpConesUETracks->SetXTitle("E_{T}");
        fPerpConesUETracks->SetYTitle("#Sigma P_{T}^{UE}");
        fPerpConesUETracks->Sumw2();
        fOutput->Add(fPerpConesUETracks);
        
        fTPCWithoutIsoConeB2BbandUE = new TH2D("hFullTPCUE","UE Estimation with almost Full TPC",200,0.,100.,250,0.,100.);
        fTPCWithoutIsoConeB2BbandUE->SetXTitle("E_{T}");
        fTPCWithoutIsoConeB2BbandUE->SetYTitle("#Sigma E_{T}^{UE}");
        fTPCWithoutIsoConeB2BbandUE->Sumw2();
        fOutput->Add(fTPCWithoutIsoConeB2BbandUE);
        
        fEtIsolatedClust = new TH1D("hEtIsolatedClust","E_{T} distribution for Isolated Photons with clusters; #Sigma E_{T}^{iso cone}<Ethres",200,0.,100.);
        fEtIsolatedClust->SetXTitle("E_{T}^{iso}");
        fEtIsolatedClust->Sumw2();
        fOutput->Add(fEtIsolatedClust);
        
        fPtIsolatedNClust = new TH1D("hEtIsolatedNClust","p_{T} distribution for neutral clusters; #Sigma p_{T}^{iso cone}<Pthres",200,0.,100.);
        fPtIsolatedNClust->SetXTitle("p_{T}^{iso}");
        fPtIsolatedNClust->Sumw2();
        fOutput->Add(fPtIsolatedNClust);
        
        fPtIsolatedNTracks = new TH1D("hEtIsolatedNTracks","p_{T} distribution for neutral clusters; #Sigma p_{T}^{iso cone}<Pthres",200,0.,100.);
        fPtIsolatedNTracks->SetXTitle("p_{T}^{iso}");
        fPtIsolatedNTracks->Sumw2();
        fOutput->Add(fPtIsolatedNTracks);
        
        fEtIsolatedCells = new TH1D("hEtIsolatedCells","E_{T} distribution for Isolated Photons with cells; #Sigma E_{T}^{iso cone}<Ethres",100,0.,100.);
        fEtIsolatedCells->SetXTitle("E_{T}^{iso}");
        fEtIsolatedCells->Sumw2();
        fOutput->Add(fEtIsolatedCells);
        
        fEtIsolatedTracks = new TH1D("hEtIsolatedTracks","E_{T} distribution for Isolated Photons with tracks; #Sigma P_{T}^{iso cone}<Pthres",100,0.,100.);
        fEtIsolatedTracks->SetXTitle("E_{T}^{iso}");
        fEtIsolatedTracks->Sumw2();
        fOutput->Add(fEtIsolatedTracks);
        
        fPtvsM02iso = new TH2D("hPtvsM02iso","p_{T} vs #lambda_{0}^{2} distribution for isolated clusters",200,0.,100.,500,0.,5.);
        fPtvsM02iso->SetXTitle("p_{T}^{iso}");
        fPtvsM02iso->SetYTitle("#lambda_{0}^{2}");
        fOutput->Add(fPtvsM02iso);
        
        fPtvsM02noiso = new TH2D("hPtvsM02noiso","p_{T} vs #lambda_{0}^{2} distribution for non isolated clusters",200,0.,100.,500,0.,5.);
        fPtvsM02noiso->SetXTitle("p_{T}^{iso}");
        fPtvsM02noiso->SetYTitle("#lambda_{0}^{2}");
        
        fOutput->Add(fPtvsM02noiso);
        fPtvsM02vsSum = new TH3D("hPtvsM02vsSum","p_{T} vs #lambda_{0}^{2} vs  #Sigma E_{T}^{iso cone} distribution for non isolated clusters",200,0.,100.,500,0.,5.,200,0.,100.);
        fPtvsM02vsSum->Sumw2();
        fOutput->Add(fPtvsM02vsSum);
        
        fPtvsM02vsSumUE = new TH3D("hPtvsM02vsSumUE","p_{T} vs #lambda_{0}^{2} vs  #Sigma E_{T}^{iso cone}-UE distribution for clusters",200,0.,100.,500,0.,5.,200,-10.,90.);
        fPtvsM02vsSumUE->Sumw2();
        fOutput->Add(fPtvsM02vsSumUE);
       
        // fPtvsM02vsSumPi0 = new TH3D("hPtvsM02vsSumPi0 when pi0 rejecting","p_{T} vs #lambda_{0}^{2} vs  #Sigma E_{T}^{iso cone}-UE  pi0 rejecting distribution for clusters",200,0.,100.,500,0.,5.,200,-10.,90.);
        // fPtvsM02vsSumPi0->Sumw2();
        // fOutput->Add(fPtvsM02vsSumPi0);
        
        // fPtvsM02vsSumEta = new TH3D("hPtvsM02vsSumEta when pi0+eta rejecting","p_{T} vs #lambda_{0}^{2} vs  #Sigma E_{T}^{iso cone}-UE  pi0+eta rejecting distribution for clusters",200,0.,100.,500,0.,5.,200,-10.,90.);
        // fPtvsM02vsSumPi0->Sumw2();
        // fOutput->Add(fPtvsM02vsSumEta);
          
        fTrackMultvsSumChargedvsUE = new TH3D("hTrackMultvsSumChargedvsUE","Track Multiplicity vs  #Sigma E_{T}^{iso cone} vs UE charged",100,0.,100.,200,-10.,90.,100,0.,100.);
        fTrackMultvsSumChargedvsUE->Sumw2();
        fOutput->Add(fTrackMultvsSumChargedvsUE);
        
        fTrackMultvsPt = new TH2D("hTrackMultvsPt","Track Multiplicity vs  p_{T}-UE distribution for clusters",100,0.,100.,200,0.,100.);
        fTrackMultvsPt->Sumw2();
        fOutput->Add(fTrackMultvsPt);

	if(fIsMC){
            fphietaPhotons = new TH3D ("fDphiDeta_Photons","#Delta#phi vs #Delta#eta Clust-MCpart to check why zero M02 clusters; #eta; #phi", 100, -0.5, 0.5, 200, 1.5, 3.5,60,0.,60.);
           fphietaPhotons->Sumw2();
            fOutput->Add(fphietaPhotons);
            
          fphietaOthers = new TH3D ("fDphiDeta_Others","#Delta#phi vs #Delta#eta Clust-MCpart to check why zero M02 clusters; #eta; #phi", 140, -0.7, 0.7, 220, 0.8, 3.5,60,0.,60.);
            fphietaOthers->Sumw2();
         fOutput->Add(fphietaOthers);
           
           fphietaOthersBis = new TH3D ("fDphiDeta_OthersBis","#Delta#phi vs #Delta#eta Clust-MCpart to check why zero M02 clusters; #eta; #phi", 140, -0.7, 0.7, 220, 0.8, 3.5,60,0.,60.);
          fphietaOthersBis->Sumw2();
           fOutput->Add(fphietaOthersBis);
	}

      }
        break;
        
    }
  }
    //Common histograms QA initialization
  if(fQA){
      //Include QA plots to the OutputList //DEFINE BETTER THE BINNING AND THE AXES LIMITS
    fTrackMult = new TH1D ("hTrackMult","Tracks multiplicity Distribution",100,0.,100.);
    fTrackMult->Sumw2();
    fOutput->Add(fTrackMult);
    
    fClusTime = new TH1D("hClusTime_NC","Time distribution for Clusters",800,-50.,50.);
    fClusTime->Sumw2();
    fOutput->Add(fClusTime);
    
    fEtaPhiClus = new TH2D ("hEtaPhiClusActivity","",250,-0.8,0.8, 250, 1.2, 3.4);
      //  fEtaPhiClus->Sumw2();
    fOutput->Add(fEtaPhiClus);
    
    fDeltaETAClusTrack = new TH1D("h_Dz","Track-Cluster Dz ",1000,-0.5,0.5);
    fDeltaETAClusTrack->Sumw2();
    fOutput->Add(fDeltaETAClusTrack);
    
    fDeltaPHIClusTrack = new TH1D("h_Dx","Track-Cluster Dx",1000,-0.5,0.5);
    fDeltaPHIClusTrack->Sumw2();
    fOutput->Add(fDeltaPHIClusTrack);
    
    fDeltaETAClusTrackMatch = new TH1D("h_DzMatch","Track-Cluster Dz matching ",100,-0.05,0.05);
    fDeltaETAClusTrackMatch ->Sumw2();
    fOutput->Add(fDeltaETAClusTrackMatch);
    
    fDeltaPHIClusTrackMatch = new TH1D("h_DxMatch","Track-Cluster Dx matching",100,-0.05,0.05);
    fDeltaPHIClusTrackMatch->Sumw2();
    fOutput->Add(fDeltaPHIClusTrackMatch);
    
    fPT = new TH1D("hPt_NC","P_{T} distribution for Neutral Clusters",100,0.,100.);
    fPT->Sumw2();
    fOutput->Add(fPT);
    
    fE = new TH1D("hE_NC","E distribution for Clusters",200,0.,100.);
    fE->Sumw2();
    fOutput->Add(fE);
    
    fNLM = new TH2D("hNLM_NC","NLM distribution for Neutral Clusters",10,0.,10.,100,0.,100.);
    fNLM->Sumw2();
    fOutput->Add(fNLM);
  
    fTestIndex= new TH2D("hTestIndex","Test index pour cluster",100,0.,100.,100,0.,100.);
    fTestIndex->SetXTitle("index");
    fTestIndex->SetYTitle("local index");
    fTestIndex->Sumw2();
    fOutput->Add(fTestIndex);
    
    fTestIndexE= new TH2D("hTestIndexE","Test index vs energy pour cluster",200,0.,100.,100,0.,100.);
    fTestIndexE->SetXTitle("cluster energy");
    fTestIndexE->SetYTitle("index");
    fTestIndexE->Sumw2();
    fOutput->Add(fTestIndexE);
    
    fTestLocalIndexE= new TH2D("hTestLocalIndexE","Test local index vs energy for cluster",200,0.,100.,100,0.,100.);
    fTestLocalIndexE->SetXTitle("cluster energy");
    fTestLocalIndexE->SetYTitle("local index");
    fTestLocalIndexE->Sumw2();
    fOutput->Add(fTestLocalIndexE);

    fTestEtaPhiCone= new TH2D("hTestEtatPhiCone","Test eta phi neutral clusters candidates",200,0,TMath::TwoPi(),140,-.7,.7);
    fTestEtaPhiCone->SetXTitle("phi");
    fTestEtaPhiCone->SetYTitle("eta");
    fTestEtaPhiCone->Sumw2();
    fOutput->Add(fTestEtaPhiCone);
 
    // fInvMassM02iso = new TH3D("hInvMassM02iso","Invariant mass vs M02 vs E_{T}^{iso cluster}",100,0.,1.,500,0.,5.,200,0.,100.);
    // fInvMassM02iso->Sumw2();
    // fOutput->Add(fInvMassM02iso);
    
    // fInvMassM02noiso = new TH3D("hInvMassM02noiso","Invariant mass vs M02 vs E_{T}^{no iso cluster}",100,0.,1.,500,0.,5.,200,0.,100.);
    // fInvMassM02noiso->Sumw2();
    // fOutput->Add(fInvMassM02noiso);
   
  }
    //   Initialization of all the Common THistos for the Three different outputs
  fVz = new TH1D("hVz_NC","Vertex Z distribution",100,-50.,50.);
  fVz->Sumw2();
  fOutput->Add(fVz);
  
  fEvents = new TH1D("hEvents_NC","Events",100,0.,100.);
  fEvents->Sumw2();
  fOutput->Add(fEvents);
  
  fPtaftTime = new TH1D("hPtaftTime_NC","p_{T} distribution for Clusters after cluster time cut",200,0.,100.);
  fPtaftTime->Sumw2();
  fOutput->Add(fPtaftTime);
  
  fPtaftCell = new TH1D("hPtaftCell_NC","p_{T} distribution for Clusters after Ncells cut",200,0.,100.);
  fPtaftCell->Sumw2();
  fOutput->Add(fPtaftCell);
  
  fPtaftNLM = new TH1D("hPtaftNLM_NC","p_{T} distribution for Clusters after NLM cut",200,0.,100.);
  fPtaftNLM->Sumw2();
  fOutput->Add(fPtaftNLM);
  
  fPtaftTM = new TH1D("hPtaftTM_NC","p_{T} distribution for Neutral Clusters",200,0.,100.);
  fPtaftTM->Sumw2();
  fOutput->Add(fPtaftTM);
  
  fPtaftDTBC = new TH1D("hPtaftDTBC_NC","p_{T} distribution for Neutral Clusters after DTBC cut",200,0.,100.);
  fPtaftDTBC->Sumw2();
  fOutput->Add(fPtaftDTBC);
  
  fPtaftFC = new TH1D("hPtaftFC_NC","p_{T} distribution for Clusters after fiducial cut",200,0.,100.);
  fPtaftFC->Sumw2();
  fOutput->Add(fPtaftFC);
 
  fTestEnergyCone= new TH3F("hTestEnergyConeVSpT","Test energy clusters and tracks in cone",200,0.,100.,250,0.,100.,250,0.,100.);
  fTestEnergyCone->SetXTitle("p_{T}^{cluster}");
  fTestEnergyCone->SetYTitle("#sum^{cone} p_{T}^{cluster}");
  fTestEnergyCone->SetZTitle("#sum^{cone} p_{T}^{track}");
  fTestEnergyCone->Sumw2();
  fOutput->Add(fTestEnergyCone);
 
    //  fTracksConeEtaPt = new TH3D("hTracksConeEtaPt","#Sigma vs #eta vs E_{T}",200,0.,100.,320,-0.8,0.8,200,0.,100.);
    // fTracksConeEtaPt->Sumw2();
    // fOutput->Add(fTracksConeEtaPt);
  
    // fTracksConeEtaM02 = new TH3D("hTracksConeEtaM02","#Sigma vs #eta vs M02",200,0.,100.,320,-0.8,0.8,500,0.,5.);
    // fTracksConeEtaM02->Sumw2();
    // fOutput->Add(fTracksConeEtaM02);
  
    // fphietaPhotons = new TH3D("hphietaPhotons","Test eta phi photons MC",250,-0.8,0.8, 250, 1.2, 3.4,200,0.,1.);
    // fOutput->Add(fphietaPhotons);
  
    // fphietaOthers = new TH3D("hphietaOthers","Test eta phi others",250,-0.8,0.8, 250, 1.2, 3.4,200,0.,1.);
    // fOutput->Add(fphietaOthers);
  
  fPtTracksVSpTNC = new TH2F ("hTrackPtSpecVSpT","Charged Particle spectrum vs pT Candidate",70,0.,70.,200,0.,20.);
  fPtTracksVSpTNC->Sumw2();
  fOutput->Add(fPtTracksVSpTNC);
  
  fPhiTracksVSclustPt  = new TH2F("hPhiTracks_vs_clustPT","Tracks phi distr vs pT Candidate",70, 0.,70., 200,0.,TMath::TwoPi());
  fPhiTracksVSclustPt->Sumw2();
  fOutput->Add(fPhiTracksVSclustPt);
  
  fEtaTracksVSclustPt  = new TH2F("hEtaTracks_vs_clustPT","Tracks eta distr vs pT Candidate",70, 0.,70., 90,-0.9,0.9);
  fEtaTracksVSclustPt->Sumw2();
  fOutput->Add(fEtaTracksVSclustPt);
  
  fCTdistVSpTNC = new TH2F ("hDistanceC_TrackVSpT","Distance between Neutral Clust and closest Track vs pT Candidate",70,0.,70.,210,-0.1,2.);
  fCTdistVSpTNC->Sumw2();
  fOutput->Add(fCTdistVSpTNC);
  
  
  if(fIsMC){
      //CREATE THE TH2 specific for the MC Analysis Maybe to be added in the THNSparse, or cloning two or three axes and add the specific MC Truth info
    fHistXsection = new TH1F("fHistXsection", "fHistXsection", 1, 0, 1);
    fHistXsection->GetXaxis()->SetBinLabel(1,"<sigma>");
    fOutput->Add(fHistXsection);
    
    fHistTrials = new TH1F("fHistTrials", "fHistTrials", 1, 0, 1);
    fHistTrials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
    fOutput->Add(fHistTrials);
  
    fPtTracksVSpTNC_MC = new TH2F ("hChargedptSpecVSpT_MC","Charged Particle spectrum vs pT Candidate",70,0.,70.,200,0.,50.);
    fPtTracksVSpTNC_MC->Sumw2();
    fOutput->Add(fPtTracksVSpTNC_MC);
    
    fpi0VSclusterVSIsolation = new TH3F ("hpi0VSclusterVSisolation","Energy pi0 vs cluster Energy vs Isolation",95,5.,100.,95,5.,100.,400,0.,100.);
    fpi0VSclusterVSIsolation->SetXTitle("particle (#pi^{0} or #eta) E");
    fpi0VSclusterVSIsolation->SetYTitle("cluster E");
    fpi0VSclusterVSIsolation->SetZTitle("Isolation");
    fpi0VSclusterVSIsolation->Sumw2();
    fOutput->Add(fpi0VSclusterVSIsolation);
    
    fpi0VSM02VSIsolation = new TH3F ("fpi0VSM02VSIsolation","Energy pi0 vs cluster M02 vs Isolation",95,5.,100.,500,0.,5.,400,0.,100.);
    fpi0VSM02VSIsolation->SetXTitle("particle (#pi^{0} or #eta) E");
    fpi0VSM02VSIsolation->SetYTitle("cluster M02");
    fpi0VSM02VSIsolation->SetZTitle("Isolation");
    fpi0VSM02VSIsolation->Sumw2();
    fOutput->Add(fpi0VSM02VSIsolation);
    
    fpi0VSclusterVSM02 = new TH3F ("fpi0VSclusterVSM02","Energy pi0 vs Energy cluster vs cluster M02 ",95,5.,100.,100,0.,100.,500,0.,5.);
    fpi0VSclusterVSM02->SetXTitle("particle (#pi^{0} or #eta) E");
    fpi0VSclusterVSM02->SetYTitle("cluster E");
    fpi0VSclusterVSM02->SetZTitle("cluster M02");
    fpi0VSclusterVSM02->Sumw2();
    fOutput->Add(fpi0VSclusterVSM02);
  }

  
  
  PostData(1, fOutput);
    //     //   return;
}

  //________________________________________________________________________
Double_t* AliAnalysisTaskEMCALPhotonIsolation::GenerateFixedBinArray(Int_t n, Double_t min, Double_t max) const
{
    // Generate the bin array for the ThnSparse
  
  Double_t *bins = new Double_t[n+1];
  
  Double_t binWidth = (max-min)/n;
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
    //tracks for CT Matching
  AliTrackContainer *tracks = GetTrackContainer("tpconlyMatch");
    //  Printf("name of the first track container: %s", tracks->GetClassName().Data());
  if (!tracks) {
    AliError(Form("%s: This task needs a 1particle container!", GetName()));
    return;
  }
    //tracks for Isolation
  AliTrackContainer *tracksANA = GetTrackContainer("filterTracksAna");
    //  Printf("name of the second track container: %s", tracksANA->GetClassName().Data());
  
  if (!tracksANA) {
    AliError(Form("%s: This task needs a 2particle container!", GetName()));
    return;
  }
    //clusters
  AliClusterContainer *clusters = GetClusterContainer(0);
  if (!clusters) {
    AliError(Form("%s: This task needs a cluster container!", GetName()));
    return;
  }
  
    //Init the EMCAL Framework
  AliAnalysisTaskEmcal::ExecOnce();
  if (!fLocalInitialized) {
    
    AliError(Form("AliAnalysisTask not initialized"));
    return;
  }
}

  //______________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::Run()
{
    // Run the analysis.
    //fTest1+=1;
    //vertex cuts
  AliTrackContainer *tracks = GetTrackContainer("tpconlyMatch");
  if(!tracks){
    Printf("Cannot find the tracks for CT Matching");
    return kFALSE;
  }
  
  if(tracks->GetTrackFilterType()!= AliEmcalTrackSelection::kTPCOnlyTracks){
    AliWarning(Form("CT matching NOT performed with TPCOnly Tracks"));
    AliWarning(Form("You better be sure of what you are doing"));
  }
  
  AliTrackContainer *tracksANA = GetTrackContainer("filterTracksAna");
  if(!tracksANA){
    Printf("Cannot find the tracks for Isolation");
    return kFALSE;
  }
  
    //  Printf("FilterType of the tracks for Analysis: %d \t(should be %d)", tracksANA->GetTrackFilterType(),AliEmcalTrackSelection::kHybridTracks);
  
    //    AliError(Form("\n\n\n\nGO CHECK the Settings!!!! Is Isolation calculated with filteredTracks?\n\n\n\n"));
  if(tracksANA->GetTrackFilterType()!= AliEmcalTrackSelection::kHybridTracks){
    AliWarning(Form("Isolation NOT calculated with HybridTracks"));
    AliWarning(Form("You better be sure of what you are doing"));
  }
  
  
  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  
    //  Printf("Vertex Z coordinate for M2: %f", fVertex[2]);
    //  Printf("Vertex Z coordinate for NF: %lf", fVertex[2]);
  
  if (fVertex[2]>10. || fVertex[2]<-10.) return kFALSE;
    //  AliError(Form("La task tourne bien"));
  
  AliClusterContainer* clusters = GetClusterContainer(0);
    ////Printf("Inside Run Method, Number of clusters for this event: %d",clusters->GetNAcceptedClusters());
  Int_t nbTracksEvent;
  nbTracksEvent =InputEvent()->GetNumberOfTracks();
  
  if(fRejectionEventWithoutTracks && nbTracksEvent==0)
    return kFALSE;
  
    // Fill events number histogram
  fEvents->Fill(0);
    //AliError(Form("one event"));
  
  if(fIsMC){
    Float_t fXSection=0;
    Int_t fTrials=0;
    if(fPythiaHeader){
      fXSection = fPythiaHeader->GetXsection();
      fTrials = fPythiaHeader->Trials();
      if(fTrials>0){
        fHistXsection->Fill("<sigma>",fXSection);
        fHistTrials->Fill("#sum{ntrials}",fTrials);
      }
    }
      // AliError(Form("EMCAL L1 trigger for MC simulation anchored to LHC13 data"));
    if(!MCSimTrigger(fVevent,fTriggerLevel1) && fAnalysispPb) return kFALSE;
  }
  
    //Fill Vertex Z histogram
  fVz->Fill(fVertex[2]);
  
    // delete output USEFUL LATER FOR CONTAINER CREATION !!
    //fOutClusters->Delete();
  Int_t index=0;
  
  if(fIsMC){
    AliAODMCHeader *mcHeader;
    
    fAODMCParticles = static_cast <TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
    
    mcHeader = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindObject(AliAODMCHeader::StdBranchName()));
      //AliError(Form("Passe analyze MC"));
    if (!fIsMC)
      return kFALSE;
      //AliInfo(Form("It's a MC analysis %e",fAODMCParticles));
    if(!fStack && !fAODMCParticles){
        AliError("no MC stack saved\n");
        return kFALSE;
    }
    
      //cout<<"there's a List of particles"<<endl;
      //DO THIS ALSO FOR ESDs
    
    if(fAODMCParticles->GetEntries() < 1){
      AliError("number of MC particles insufficient");
      return kFALSE;
    }
    AnalyzeMC();
  }
  
  if (fisLCAnalysis) {
    
      // get the leading particle
      // AliEmcalParticle *emccluster = static_cast<AliEmcalParticle*>(clusters->GetLeadingParticle());
    
    AliVCluster *coi = (clusters->GetLeadingCluster());
    
    if(!coi){
      
      AliError(Form("No leading cluster"));
      return kFALSE;
    }
    
    
      //   index = coi->IdInCollection();
    index = coi->GetID();
      //add a command to get the index of the leading cluster!
    if (!coi->IsEMCAL()) return kFALSE;
    
      // AliVCluster *coi = emccluster->GetCluster();
      // if (!coi) return kFALSE;
    
    TLorentzVector vecCOI;
    coi->GetMomentum(vecCOI,fVertex);
    
    
    Double_t coiTOF = coi->GetTOF()*1e9;
    if(coiTOF<-30. || coiTOF>30.)
      return kFALSE;
    
    
    if(ClustTrackMatching(coi,kTRUE))
      return kFALSE;
    
    if(!CheckBoundaries(vecCOI))
      return kFALSE;
    
    else
      FillGeneralHistograms(coi,vecCOI, index);
  }
  else{
      //get the entries of the Cluster Container
      //whatever is a RETURN in LCAnalysis here is a CONTINUE,
      //since there are more than 1 Cluster per Event
    
      //    Printf("Starting loop on clusters in RUN method");
    for (auto it : clusters->accepted()){
      AliVCluster *coi = static_cast<AliVCluster*>(it);
      if(!coi) {
        AliError("No cluster found");
        return kFALSE;
      }
        //
      
      index=coi->GetID();
      TLorentzVector vecCOI;
      coi->GetMomentum(vecCOI,fVertex);
      Double_t coiTOF = coi->GetTOF()*1e9;
      Double_t coiM02 = coi->GetM02();
      
      if(fQA)  FillQAHistograms(coi,vecCOI);
      
      if(!fIsMC){
        if(coiTOF< -30. || coiTOF > 30.)
          continue;
      }
      fPtaftTime->Fill(vecCOI.Pt());
      
      if((coi->GetNCells() < 2))
        continue;
      
      fPtaftCell->Fill(vecCOI.Pt());
      
      Int_t nlm=0;
      AliVCaloCells * fCaloCells =InputEvent()->GetEMCALCells();
      if(fCaloCells)
      {
        nlm = GetNLM(coi,fCaloCells);
        AliDebug(1,Form("NLM = %d",nlm));
        
          // if a NLM cut is define, this is a loop to reject clusters with more than the defined NLM (should be 1 or 2 (merged photon decay clusters))
        if(coi->E()>=5. && coi->E()<70. && fQA)
          fNLM->Fill(nlm,coi->E());
        if(fIsNLMCut && fNLMCut>0 && fNLMmin>0)
          if(nlm > fNLMCut || nlm < fNLMmin ){
            //AliWarning(Form("NLM = %d --- NLM min = %d --- NLMcut = %d",nlm,fNLMmin,fNLMCut));
            continue;
          }
      }
      else
        AliDebug(1,Form("Can't retrieve EMCAL cells"));
      
      fPtaftNLM->Fill(vecCOI.Pt());
      if(fTMClusterRejected)
      {
        if(ClustTrackMatching(coi,kTRUE)){
          continue;
        }
      }
      fPtaftTM->Fill(vecCOI.Pt());
      
      if((coi->GetDistanceToBadChannel() < 2)) continue;
      
      fPtaftDTBC->Fill(vecCOI.Pt());
      
      
      
      if(!CheckBoundaries(vecCOI))
        continue;
      
      fPtaftFC->Fill(vecCOI.Pt());
      if(fQA){
        fTestIndexE->Fill(vecCOI.Pt(),index);
      }
      
      if(vecCOI.Pt()<5.) continue;
      
        //      Printf("Inside Run: Passing to FillGeneralHistograms for cluster with ID: %d, Pt: %.3f, Eta: %.3f, Phi: %.3f",index,vecCOI.Pt(),vecCOI.Eta(),vecCOI.Phi());
      for (auto it : tracksANA->accepted()){
        AliVTrack *tr = static_cast<AliVTrack*>(it);
        if(!tr) {
          AliError("No track found");
          return kFALSE;
        }
        fPtTracksVSpTNC->Fill(vecCOI.Pt(),tr->Pt());
        fPhiTracksVSclustPt->Fill(vecCOI.Pt(),tr->Phi());
        fEtaTracksVSclustPt->Fill(vecCOI.Pt(),tr->Eta());
      }
      
      FillGeneralHistograms(coi,vecCOI,index);
    }
    
  }
    //  PostData(1, fOutput);
  return kTRUE;
}


  //_________________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::FillQAHistograms(AliVCluster *coi, TLorentzVector vecCOI){
  
  switch(fWho){
      
    case 0:
      fevents=0;
      fEClustersT=vecCOI.E();
      fPtClustersT=vecCOI.Pt();
      fEtClustersT=vecCOI.Et();
      fEtaClustersT=vecCOI.Eta();
      fPhiClustersT=vecCOI.Phi();
      fM02ClustersT=coi->GetM02();
      
      fOutputQATree->Fill();
      
      break;
      
    case 1:
      
      break;

    case 2:
      
      fM02->Fill(vecCOI.E(),coi->GetM02());

      break;

      
  }
  
  fPT->Fill(vecCOI.Pt());
  fE->Fill(vecCOI.E());
  fEtaPhiClus->Fill(vecCOI.Eta(),vecCOI.Phi());
  
  Double_t checktof = coi->GetTOF()*1e9;
  fClusTime->Fill(checktof);

  if(checktof>-30. && checktof<30. && !fIsMC){
      // fPtaftTime->Fill(vecCOI.Pt());
      //  if(!ClustTrackMatching(coi)){
      //  fPtaftTM->Fill(vecCOI.Pt());
    
    if(CheckBoundaries(vecCOI)){
        //  fPtaftFC->Fill(vecCOI.Pt());
      
      Double_t checkM02=coi->GetM02();
      if(fM02mincut < checkM02 && checkM02 < fM02maxcut && fWho==2){
        fPtaftM02C->Fill(vecCOI.Pt());
      }
    }
      // }
  }
}

  //___________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::MCSimTrigger(AliVEvent *eventIn, Int_t triggerLevel){
  
  if(triggerLevel<1 || triggerLevel>2){
    AliError(Form("No trigger level have been choosed for the MC analysis"));
    return kFALSE;
  }
  Double_t threshold=0;
  Double_t spread=0;
  
    //Int_t runNumber = InputEvent()->GetRunNumber();
  
  Int_t runNumber = eventIn->GetRunNumber();
    //AliError(Form("Le numéro de run est le %d",runNumber));
  
  if(!runNumber) return kFALSE;
  
  
  if(runNumber < 195180 || runNumber > 197469) return kFALSE;
  
    //AliError(Form("est bien après le numéro de run %d",runNumber));
  
    // TString fired = InputEvent() ->GetFiredTriggerClasses();
  
    //AliError(Form("Trigger utilisés dans les evts %s",fired));
  
  if(triggerLevel==1){
    threshold = 11.5;
    spread = 0.5;
  }
  
  if(triggerLevel==2){
    threshold = 7.2;
    spread = 0.3;
  }
  
  if (spread != 0.){
    TF1* triggerSmearing =  new TF1("triggerSmearing","[0]*exp(-0.5*((x-[1])/[2])**2)",0,15);
    triggerSmearing->SetParameter(0, 1/(spread*TMath::Sqrt(TMath::Pi()*2)));
    triggerSmearing->SetParameter(1, threshold);
    triggerSmearing->SetParameter(2, spread);
    threshold = triggerSmearing->GetRandom();
    delete triggerSmearing;
  }
  
    //  AliError(Form("passe bien la définition de la fonction du trigger"));
  AliClusterContainer *clusters = GetClusterContainer(0);
  Int_t localIndex=0;
  
  for(auto it : clusters->accepted()){
    AliVCluster* coi = static_cast<AliVCluster*>(it);
      //  AliError(Form("Recupère bien les clusters"));
    
    if(!coi) continue;
      //    AliError(Form("Recupère bien le cluster"));
    if(coi->E() > threshold){
        //	AliError(Form("Un cluster passe le critère en énergie"));
      return kTRUE;
    }
  }
  return kFALSE;
}


  //__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::ClustTrackMatching(AliVCluster *clust, Bool_t candidate) {
    // Check if the cluster match to a track
  
  AliTrackContainer* tracks = GetTrackContainer(0);
  AliVTrack* mt = 0;
  TLorentzVector vecClust;
  clust->GetMomentum(vecClust,fVertex);
  
  Int_t nbMObj = clust -> GetNTracksMatched();
  if(tracks->GetTrackFilterType()!=AliEmcalTrackSelection::kTPCOnlyTracks)  AliError(Form("NO TPC only tracks"));
  
  Double_t distCT=0.;
  
  if (nbMObj == 0) return kFALSE;
  
  for(Int_t i=0;i<nbMObj;i++){
    
    if (fIsEsd) {
      Int_t imt = clust->GetTrackMatchedIndex(0);
      if (imt >= 0) mt = static_cast<AliVTrack*>(tracks->GetAcceptParticle(imt));
    }
    else {
      mt = static_cast<AliVTrack*>(clust->GetTrackMatched(i));
      UInt_t rejectionReason = 0;
      if (!tracks->AcceptParticle(mt, rejectionReason)) mt = 0;
    }
      //  Int_t imt = partC->GetMatchedObjId(i);
    
    if(!mt) continue;
    
      //    printf("Cluster ID %d matched with track ID %d with pT %.3f",clust->GetID(),mt->GetID(),mt->Pt());
    
    Double_t deltaEta,deltaPhi;
    Double_t deta = 999;
    Double_t dphi = 999;
    
    Double_t veta = mt->GetTrackEtaOnEMCal();
    Double_t vphi = mt->GetTrackPhiOnEMCal();
    
    Float_t pos[3] = {0};
    clust->GetPosition(pos);
    TVector3 cpos(pos);
    Double_t ceta     = cpos.Eta();
    Double_t cphi     = cpos.Phi();
    deta=veta-ceta;
    dphi=TVector2::Phi_mpi_pi(vphi-cphi);
      //    printf("distant deta %.3f and dphi %.3f from the cluster",deta, dphi);
    if(fQA && candidate){
      fDeltaETAClusTrack->Fill(deta);
      fDeltaPHIClusTrack->Fill(dphi);
    }
    distCT=TMath::Sqrt(deta*deta+dphi*dphi);
    fCTdistVSpTNC->Fill(vecClust.Pt(),distCT);
    
      //      Printf("dphimin %g dphi %g   i %d  matchidx %d",dphimin,dphi,i, matchidxphi);
    if(candidate){
      deltaEta=fdetacut;
      deltaPhi=fdphicut;
    }
    else{
      deltaEta=fdetacutIso;
      deltaPhi=fdphicutIso;
    }
   
    if(TMath::Abs(dphi)<deltaPhi && TMath::Abs(deta)<deltaEta){
      if(fQA && candidate){
        fDeltaETAClusTrackMatch->Fill(deta);
        fDeltaPHIClusTrackMatch->Fill(dphi);
      }
      return kTRUE;
    }
  }
  return kFALSE;
}

  //_____________________________________________________________________________________________
Int_t AliAnalysisTaskEMCALPhotonIsolation::GetNLM(AliVCluster *coi, AliVCaloCells* cells){
    // find the number of local maxima of a cluster adapted from AliCalorimeterUtils
  
  const Int_t   nc = coi->GetNCells();
  
  Int_t   absIdList[nc];
  Float_t maxEList[nc];
  
  Int_t nMax = GetNLM(coi, cells, absIdList, maxEList);
  
  return nMax;
}

  //_____________________________________________________________________________________________________________________________
Int_t AliAnalysisTaskEMCALPhotonIsolation::GetNLM(AliVCluster* coi, AliVCaloCells* cells, Int_t *absIdList, Float_t *maxEList) {
    // find the cluster number of local maxima adapted from AliCalorimeterUtils
  
  Int_t iDigitN = 0 ;
  Int_t iDigit  = 0 ;
  Int_t absId1 = -1 ;
  Int_t absId2 = -1 ;
  const Int_t nCells = coi->GetNCells();
  
    // Printf("Cluster Energy Before Recalculation: %.4f",coi->E());
  Float_t eCluster = RecalEnClust(coi, cells);// recalculate cluster energy, avoid non lin correction.
  Float_t localMaxCutE = 0.1;
  Float_t locMaxCutEDiff = 0.0;
  
    //  Float_t simuTotWeight = 0;
    //  if(fMCECellClusFracCorrOn)
    //  {
    //    simuTotWeight = RecalibrateClusterEnergyWeightCell(cluster, cells,eCluster);// same but apply a weight
    //    simuTotWeight/= eCluster;
    //  }
  
    //printf("cluster : ncells %d \n",nCells);
  
  Float_t emax  = 0;
  Int_t   idmax =-1;
  for(iDigit = 0; iDigit < nCells ; iDigit++)
  {
    absIdList[iDigit] = coi->GetCellsAbsId()[iDigit]  ;
    Float_t en = cells->GetCellAmplitude(absIdList[iDigit]);
      //    Printf("Cell Energy Before Recalculation: %.4f",en);
    RecalAmpCell(en,absIdList[iDigit]);
    
      //    if(fMCECellClusFracCorrOn)
      //      en*=GetMCECellClusFracCorrection(en,eCluster)/simuTotWeight;
    
    if( en > emax )
    {
      emax  = en ;
      idmax = absIdList[iDigit] ;
    }
  }
  for(iDigit = 0 ; iDigit < nCells; iDigit++)
  {
    if( absIdList[iDigit] >= 0 )
    {
      absId1 = coi->GetCellsAbsId()[iDigit];
      
      Float_t en1 = cells->GetCellAmplitude(absId1);
        //      Printf("Cell-1 Energy Before Recalculation: %.4f",en1);
      RecalAmpCell(en1,absId1);
      
        //      if(fMCECellClusFracCorrOn)
        //        en1*=GetMCECellClusFracCorrection(en1,eCluster)/simuTotWeight;
      
        //printf("%d : absIDi %d, E %f\n",iDigit, absId1,en1);
      
      for(iDigitN = 0; iDigitN < nCells; iDigitN++)
      {
        absId2 = coi->GetCellsAbsId()[iDigitN] ;
        
        if(absId2==-1 || absId2==absId1) continue;
        
          //printf("\t %d : absIDj %d\n",iDigitN, absId2);
        
        Float_t en2 = cells->GetCellAmplitude(absId2);
          //        Printf("Cell-2 Energy Before Recalculation: %.4f",en2);
        RecalAmpCell(en2,absId2);
        
          //        if(fMCECellClusFracCorrOn)
          //          en2*=GetMCECellClusFracCorrection(en2,eCluster)/simuTotWeight;
        
          //printf("\t %d : absIDj %d, E %f\n",iDigitN, absId2,en2);
        
        if ( AreNeighbours(absId1, absId2) )
        {
            // printf("\t \t Neighbours \n");
          if ( en1 > en2 )
          {
            absIdList[iDigitN] = -1 ;
              //printf("\t \t indexN %d not local max\n",iDigitN);
              // but may be digit too is not local max ?
            if(en1 < en2 + locMaxCutEDiff) {
                //printf("\t \t index %d not local max cause locMaxCutEDiff\n",iDigit);
              absIdList[iDigit] = -1 ;
            }
          }
          else
          {
            absIdList[iDigit] = -1 ;
              //printf("\t \t index %d not local max\n",iDigitN);
              // but may be digitN too is not local max ?
            if(en1 > en2 - locMaxCutEDiff)
            {
              absIdList[iDigitN] = -1 ;
                //printf("\t \t indexN %d not local max cause locMaxCutEDiff\n",iDigit);
            }
          }
        } // if Are neighbours
          //else printf("\t \t NOT Neighbours \n");
      } // while digitN
    } // slot not empty
  } // while digit
  
  iDigitN = 0 ;
  for(iDigit = 0; iDigit < nCells; iDigit++)
  {
    if( absIdList[iDigit] >= 0 )
    {
      absIdList[iDigitN] = absIdList[iDigit] ;
      
      Float_t en = cells->GetCellAmplitude(absIdList[iDigit]);
        // Printf("NlocMax Cell Energy Before Recalculation: %.4f",en);
      RecalAmpCell(en,absIdList[iDigit]);
        // Printf("NlocMax Cell Energy After Recalculation: %.4f",en);
      
        //      if(fMCECellClusFracCorrOn)
        //        en*=GetMCECellClusFracCorrection(en,eCluster)/simuTotWeight;
      
      if(en < localMaxCutE) continue; // Maxima only with seed energy at least
      
      maxEList[iDigitN] = en ;
      
        //printf("Local max %d, id %d, en %f\n", iDigit,absIdList[iDigitN],en);
      iDigitN++ ;
    }
  }
  
  if ( iDigitN == 0 ){
    AliDebug(1,Form("No local maxima found, assign highest energy cell as maxima, id %d, en cell %2.2f, en cluster %2.2f",
                    idmax,emax,coi->E()));
    iDigitN      = 1     ;
    maxEList[0]  = emax  ;
    absIdList[0] = idmax ;
  }
  
  AliDebug(1,Form("In coi E %2.2f (wth non lin. %2.2f), M02 %2.2f, M20 %2.2f, N maxima %d",
                  coi->E(),eCluster, coi->GetM02(),coi->GetM20(), iDigitN));
  return iDigitN ;
}

  //__________________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::AreNeighbours(Int_t absId1, Int_t absId2 ) const
{
    // check if two cells are neighbour (adapted from AliCalorimeterUtils)
  
  Bool_t areNeighbours = kFALSE ;
  
  Int_t iSupMod1 = -1, iTower1 = -1, iIphi1 = -1, iIeta1 = -1, iphi1 = -1, ieta1 = -1;
  Int_t iSupMod2 = -1, iTower2 = -1, iIphi2 = -1, iIeta2 = -1, iphi2 = -1, ieta2 = -1;
  
  Int_t phidiff =  0, etadiff =  0;
  
    //first cell
  fGeom->GetCellIndex(absId1,iSupMod1,iTower1,iIphi1,iIeta1);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod1,iTower1,iIphi1, iIeta1,iphi1,ieta1);
  
    // second cell
  fGeom->GetCellIndex(absId2,iSupMod2,iTower2,iIphi2,iIeta2);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod2,iTower2,iIphi2, iIeta2,iphi2,ieta2);
  
  
  if(iSupMod1!=iSupMod2){
      // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2-1
      // C Side impair SM, nSupMod%2=1; A side pair SM nSupMod%2=0
    if(iSupMod1%2) ieta1+=AliEMCALGeoParams::fgkEMCALCols;
    else           ieta2+=AliEMCALGeoParams::fgkEMCALCols;
  }
  
  phidiff = TMath::Abs( iphi1 - iphi2 ) ;
  etadiff = TMath::Abs( ieta1 - ieta2 ) ;
  
    //if (( coldiff <= 1 )  && ( rowdiff <= 1 ) && (coldiff + rowdiff > 0))
  if ((etadiff + phidiff == 1 ))
    areNeighbours = kTRUE ;
  
  return areNeighbours;
}
  //_____________________________________________________________________________________________
  /// Recalculate cell energy if recalibration factor.
  //_____________________________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::RecalAmpCell(Float_t & amp, Int_t id) const
{
  Int_t iSupMod = -1, iTower = -1, iIphi = -1, iIeta = -1, iphi = -1, ieta = -1;
  Float_t ampold=amp;
  fGeom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
  
  amp *= fEMCALRecoUtils->GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
    // if(ampold!=amp)
    //  Printf("Inside cell id %d\t\tOld Amp %.4f and recalculated amplitude: %.4f",id, ampold, amp);
  
}

  //__________________________________________________________________________
  /// Recalibrate the cluster energy, considering the recalibration map and the energy of the cells that compose the cluster.
  //__________________________________________________________________________
Float_t AliAnalysisTaskEMCALPhotonIsolation::RecalEnClust(AliVCluster * coi,
                                                          AliVCaloCells * cells)
{
    //Printf("Inside clust Recal");
    // Initialize some used variables
  Float_t frac  = 0., energy = 0.;
  
  if(cells)
  {
      //Get the cluster number of cells and list of absId, check what kind of cluster do we have.
    
    UShort_t * index    = coi->GetCellsAbsId() ;
    Double_t * fraction = coi->GetCellsAmplitudeFraction() ;
    
    Int_t ncells     = coi->GetNCells();
    
    
      // Loop on the cells, get the cell amplitude and recalibration factor, multiply and and to the new energy
    for(Int_t icell = 0; icell < ncells; icell++)
    {
      Int_t absId = index[icell];
      
      frac =  fraction[icell];
      if(frac < 1e-3) frac = 1; //in case of EMCAL, this is set as 0, not used.
      
      Float_t amp = cells->GetCellAmplitude(absId);
      RecalAmpCell(amp, absId);
      
        //Printf("Recalibrate cell: EMCAL, cell fraction %f, cell energy: before cal %f; after cal %f",frac,cells->GetCellAmplitude(absId),amp);
      
      energy += amp*frac;
    }
    
    AliDebug(1,Form("Energy before %f, after %f",coi->E(),energy));
    
  } // cells available
  else
  {
    AliFatal("Cells pointer does not exist!");
  }
    //Printf("recalculated energy: %.4f",energy);
  return energy;
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoCellPhiBand(TLorentzVector c, Double_t &etIso, Double_t &phiBandcells){
    // Underlying events study with EMCAL cells in phi band // have to be tested
  
  AliEMCALGeometry* emcalGeom = AliEMCALGeometry::GetInstance();
  
  Double_t sumEnergyPhiBandCells=0., sumEnergyConeCells=0.;
  
  
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
            Int_t iabsId = emcalGeom->GetAbsCellIdFromCellIndexes(iModule,iRow,iCol);  // verifier les iRow et iCol
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
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoCellEtaBand(TLorentzVector c, Double_t &etIso, Double_t &etaBandcells){
    // Underlying events study with EMCAL cell in eta band // have to be tested
  
  AliEMCALGeometry* emcalGeom = AliEMCALGeometry::GetInstance();
  
  Double_t sumEnergyEtaBandCells=0., sumEnergyConeCells=0.;
  
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
          if(iCol<iColMaxCone && iCol>iColMinCone){
            Int_t iabsId = emcalGeom->GetAbsCellIdFromCellIndexes(iModule,iRow,iCol);  // verifier les iRow et iCol
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
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoClusPhiBand(TLorentzVector c, Double_t &ptIso, Double_t &phiBandclus, Int_t index){
    // Underlying events study with clusters in phi band
  
  Double_t sumEnergyPhiBandClus=0., sumEnergyConeClus=0., sumpTConeCharged=0., sumpTPhiBandTracks=0.;
  Double_t clustTOF,phiClust,etaClust, radius;
  Double_t minPhi = 1.4, maxPhi = 3.15, minEta=-0.67, maxEta=0.67;
  
    //needs a check on the same cluster
    // AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliClusterContainer *clusters = GetClusterContainer(0);
  Int_t localIndex=0;
  TLorentzVector nClust; //STILL NOT INITIALIZED
  
  for(auto it : clusters->accepted()){ //check the position of other clusters in respect to the trigger cluster
    
    AliVCluster* coi = static_cast<AliVCluster*>(it);
    
    localIndex=coi->GetID();
    
    if(localIndex==index) continue;
    
    phiClust = etaClust = clustTOF = 0.;
    coi->GetMomentum(nClust,fVertex);
    phiClust =nClust.Phi();
    etaClust= nClust.Eta();
    
    if(fExtraIsoCuts){
      if((coi->GetNCells() < 2)) continue;
      if((coi->GetDistanceToBadChannel() < 2)) continue;
    }
    
    clustTOF = coi->GetTOF()*1e9;
    
    if(!fIsMC)
      if(clustTOF<-30. || clustTOF>30.) continue;
    
    if(fTMClusterInConeRejected)
      if(ClustTrackMatching(coi,kFALSE)) continue;
    
    if(nClust.E()<0.3) continue;
      //redefine phi/c.Eta() from the cluster we passed to the function
    
    radius = TMath::Sqrt(TMath::Power(phiClust- c.Phi(),2)+TMath::Power(etaClust-c.Eta(),2)); // define the radius between the leading cluster and the considered cluster
    
    if(radius>fIsoConeRadius){ // the cluster is outside the isolation cone in this case study UE
                               // actually phi band here
      if(TMath::Abs(phiClust - c.Phi()) < fIsoConeRadius)
        sumEnergyPhiBandClus += nClust.Pt();
    }
    else if(radius<fIsoConeRadius && radius!=0.){  // if the cluster is in the isolation cone, add the cluster pT
                                                   //Printf("Cluster Not Matched, Inside the Cone, with Energy %.4lf",eTcluster);
      sumEnergyConeClus += nClust.Pt();
      if(fQA){
        fTestEtaPhiCone->Fill(c.Eta(),c.Phi());
        fTestIndex->Fill(index,localIndex);
      
        fTestLocalIndexE->Fill(nClust.Pt(),localIndex);
      }
    }
  }
  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
    //fTracksAna = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("FilterTracksAna"));
    // name hard coded to use the defined tracks for analysis
  
    //if (!fTracksAna) {
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
    //  const Int_t nbTracks = fTracksAna->GetEntries();
    //  Int_t iTracks = 0;
  
  tracksAna->ResetCurrentID();
  AliVTrack *eTrack = 0x0;
  AliAODTrack *aodEtrack=0x0;
  Int_t iTracksCone = 0;
  Double_t phiTrack, etaTrack;
  
  while((eTrack = static_cast<AliVTrack*>(tracksAna->GetNextAcceptParticle()))){
    
    if(!eTrack){
      AliError(Form("No tracks in collection"));
      continue;
    }
      //    if(!(eTrack->IsHybridGlobalConstrainedGlobal())){Printf("skipping track %d because it's not an hybrid\n",eTrack->GetID()); continue;}
    
    if(!fIsEsd){
      aodEtrack = static_cast<AliAODTrack*>(eTrack);
      if(!(aodEtrack->IsHybridGlobalConstrainedGlobal())) continue;
    }
    
    if((eTrack->Pt())<0.2)
      continue;
    
    if((eTrack->GetStatus() & AliAODTrack::kITSrefit)!=AliAODTrack::kITSrefit) continue;
    
    if(!fIsEsd){
      Double_t frac = 0;
      Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
      Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if (  ncls> 0 )  frac =  nclsS / ncls ;
      
      if(frac > 0.4) continue;
    }
    
      //CHECK IF TRACK IS IN BOUNDARIES
    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();
      //skip track outside EMCAL
    if(etaTrack < minEta || etaTrack > maxEta || phiTrack < minPhi || phiTrack > maxPhi) continue;
    
    radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2));
    
    if(radius<fIsoConeRadius ){ // if tracks are outside the isolation cone study
      sumpTConeCharged+=eTrack->Pt(); // should not double count if the track matching is already done
      iTracksCone++;
    }
    else{
        // actually phi band here for Tracks also
      if(TMath::Abs(etaTrack - c.Eta()) < fIsoConeRadius){
          //        Printf("but phi similar, so UE EtaBand!!");
        sumpTPhiBandTracks += eTrack->Pt();
      }
    }
  } // end of tracks loop
  
  fTestEnergyCone->Fill(nClust.Pt(),sumEnergyConeClus,sumpTConeCharged);
  
  if(fWho==2){
    fTrackMultvsSumChargedvsUE->Fill(iTracksCone,sumpTConeCharged, sumpTPhiBandTracks);
    fTrackMultvsPt->Fill(iTracksCone,c.Pt());
  }
  
  ptIso = sumEnergyConeClus + sumpTConeCharged;
  
  if(fWho==2){
    fPhiBandUEClust->Fill(c.Pt() , sumEnergyPhiBandClus);
    fPhiBandUETracks->Fill(c.Pt() , sumpTPhiBandTracks);
  }
  phiBandclus = sumEnergyPhiBandClus + sumpTPhiBandTracks;
}


  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoClusEtaBand(TLorentzVector c, Double_t &ptIso, Double_t &etaBandclus, Int_t index){
    // Underlying events study with clusters in eta band
  
  Float_t sumEnergyEtaBandClus =0., sumEnergyConeClus=0., sumpTConeCharged=0, sumpTEtaBandTracks=0.;
  Double_t clustTOF=0;
  Double_t minEta = -0.67, maxEta = 0.67, minPhi=1.398,maxPhi=TMath::Pi();
  
    // AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliClusterContainer *clusters = GetClusterContainer(0);
  Int_t localIndex=0;
  
  TLorentzVector nClust; //STILL NOT INITIALIZED
  
  for(auto it : clusters->accepted()){ //check the position of other clusters in respect to the trigger cluster
    
    AliVCluster* coi = static_cast<AliVCluster*>(it);
    
    localIndex=coi->GetID();
      //  AliError(Form("Tourne bien sur les clusters"));
    if(localIndex==index) continue;
    
    coi->GetMomentum(nClust,fVertex);
    
    Double_t phiClust =nClust.Phi();
    Double_t etaClust= nClust.Eta();
    Double_t eTcluster=0, radius;
    
    if(fExtraIsoCuts){
      if((coi->GetNCells() < 2)) continue;
      if((coi->GetDistanceToBadChannel() < 2)) continue;
    }
      //    printf("\nCluster ID %d with pT: %.4f\t\t Eta: %.4f \t Phi: %.4f \t time %.3f", coi->GetID(), nClust.Pt(),etaClust,phiClust,coi->GetTOF()*1e9);
    clustTOF = coi->GetTOF()*1e9;
    if(!fIsMC)
      if(clustTOF<-30. || clustTOF>30.) continue;
    
    if(fTMClusterInConeRejected)
      if(ClustTrackMatching(coi,kFALSE)) continue;
    
    if(nClust.E()<0.3) continue;
      //redefine phi/c.Eta() from the cluster we passed to the function
    
      // define the radius between the leading cluster and the considered cluster
    radius = TMath::Sqrt(TMath::Power(phiClust-c.Phi(),2)+TMath::Power(etaClust-c.Eta(),2));
      //    printf("\tCluster Not Matched, with distance from Candidate %.3f",radius);
    if(radius>fIsoConeRadius){ // the cluster is outside the isolation cone in this case study UE
                               // actually eta band here
      if(TMath::Abs(etaClust - c.Eta()) < fIsoConeRadius)
        sumEnergyEtaBandClus += nClust.Pt();
    }
    else if(radius<fIsoConeRadius && radius!=0.){
        // if the cluster is in the isolation cone, add the cluster pT
        //      printf("\t Inside the Cone, with pT %.4lf",nClust.Pt());
      
      sumEnergyConeClus += nClust.Pt();
      if(fQA){
        fTestEtaPhiCone->Fill(c.Eta(),c.Phi());
        fTestIndex->Fill(index,localIndex);
      
        fTestLocalIndexE->Fill(nClust.Pt(),localIndex);
      }
    }
  } // end of clusters loop
    //  Printf("\ntotal activity in isolation Cone from Clusters %.4lf",sumEnergyConeClus);
  
  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
    //  Printf("ntracks accepted: %d",tracksAna->GetNAcceptedTracks());
  if(tracksAna->GetTrackFilterType()!=AliEmcalTrackSelection::kHybridTracks)
    AliError(Form("NOT Hybrid Tracks"));
  
    //  Printf("Name of the tracks used for Isolation: %s",tracksAna->GetName());
  tracksAna->ResetCurrentID();
  AliVTrack *eTrack = 0x0;
  Int_t iTracksCone = 0;
  AliAODTrack *aodEtrack=0x0;
    //  Printf("Starting loop on Tracks in EtIsoClusEtaBand");
  while((eTrack = static_cast<AliVTrack*>(tracksAna->GetNextAcceptParticle()))){
    if(!eTrack){
      AliError(Form("No tracks in collection"));
      continue;
    }
    Double_t phiTrack, etaTrack,radius;
    
      //AliError("On a bien des traces pour l'analysos");
    if(!fIsEsd){
      aodEtrack = static_cast<AliAODTrack*>(eTrack);
      if(!(aodEtrack->IsHybridGlobalConstrainedGlobal())){
          //Printf("skipping track %d because it's not an hybrid\n",eTrack->GetID());
        continue;
      }
    }
    
    if((eTrack->Pt())<0.2)
      continue;
    
    if((eTrack->GetStatus() & AliAODTrack::kITSrefit)!=AliAODTrack::kITSrefit) continue;
    
    if(!fIsEsd){
      Double_t frac = 0;
      Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
      Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if (  ncls> 0 )  frac =  nclsS / ncls ;
      
      if(frac > 0.4) continue;
    }
      //CHECK IF TRACK IS IN BOUNDARIES
    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();
      //skip track outside EMCAL
    
    if(etaTrack < minEta || etaTrack > maxEta || phiTrack < minPhi || phiTrack > maxPhi) continue;
      //    Printf("Track with pT: %.4lf\t\t Eta: %.4lf \t Phi: %.4lf", eTrack->Pt(),etaTrack,phiTrack);
    radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2));
      //    Printf("Distance C-T: %.4lf",radius);
    if(radius<fIsoConeRadius){ // if tracks are outside the isolation cone study
      sumpTConeCharged += eTrack->Pt(); // should not double count if the track matching is already done
      iTracksCone++;
        //      Printf("Inside with pT %.4f",eTrack->Pt());//	AliError(Form("Tourne bien aussi sur les traces %f",sumpTConeCharged));
    }
    else{
        // actually eta band here for Tracks also
      if(TMath::Abs(phiTrack - c.Phi()) < fIsoConeRadius){
          //        Printf("but phi similar, so UE EtaBand!!");
        sumpTEtaBandTracks += eTrack->Pt();
      }
    }
  } // end of tracks loop
    //  Printf("\ntotal activity in isolation Cone from Tracks %.4lf",sumpTConeCharged);
  
  
  fTestEnergyCone->Fill(c.Pt(),sumEnergyConeClus,sumpTConeCharged);
  
  if(fWho==2){
    fTrackMultvsSumChargedvsUE->Fill(iTracksCone,sumpTConeCharged, sumpTEtaBandTracks);
    fTrackMultvsPt->Fill(iTracksCone,c.Pt());
  }
  
  ptIso = sumEnergyConeClus + sumpTConeCharged;
  if(fWho==2){
    fEtaBandUEClust->Fill(c.Pt() , sumEnergyEtaBandClus);
    fEtaBandUETracks->Fill(c.Pt() , sumpTEtaBandTracks);
  }
  etaBandclus = sumEnergyEtaBandClus +sumpTEtaBandTracks;
}


  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackPhiBand(TLorentzVector c, Double_t &ptIso, Double_t &phiBandtrack){
    // Underlying events study with tracks in phi band in EMCAL acceptance
  
    //INSERT BOUNDARIES ACCORDING TO THE FLAG WE WANT TO USE
  Double_t sumpTConeCharged=0., sumpTPhiBandTrack=0.;
  Double_t minPhi= 0., maxPhi= 2*TMath::Pi(), minEta = -0.87, maxEta= 0.87;
  
  if(!fTPC4Iso){
    minEta = -0.67;
    maxEta = 0.67;
    minPhi = 1.4;
    maxPhi = 3.15;
  }
  
  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  
    //  Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());
  tracksAna->ResetCurrentID();
  AliVTrack *eTrack = 0x0;
  Double_t phiTrack,etaTrack,radius;
  Int_t iTracksCone = 0;
  AliAODTrack *aodEtrack=0x0;
  while((eTrack = static_cast<AliVTrack*>(tracksAna->GetNextAcceptParticle()))){
    if(!eTrack){
      AliError(Form("No tracks in collection"));
      continue;
    }
    if(!fIsEsd){
      aodEtrack = static_cast<AliAODTrack*>(eTrack);
      if(!(aodEtrack->IsHybridGlobalConstrainedGlobal())) continue;
    }
    
    if((eTrack->Pt())<0.2)
      continue;
    
    if((eTrack->GetStatus() & AliAODTrack::kITSrefit)!=AliAODTrack::kITSrefit) continue;
    
    if(!fIsEsd){
      Double_t frac = 0;
      Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
      Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if (  ncls> 0 )  frac =  nclsS / ncls ;
      
      if(frac > 0.4) continue;
    }
    
      //CHECK IF TRACK IS IN BOUNDARIES
    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();
      // define the radius between the leading cluster and the considered cluster
      //redefine phi/c.Eta() from the cluster we passed to the function
    if( (phiTrack < maxPhi) && (phiTrack > minPhi) && (etaTrack < maxEta) && (etaTrack > minEta)){
      radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2));
      
      if(radius>fIsoConeRadius){ // if tracks are outside the isolation cone study
        
          // actually phi band here --- ADD Boundaries conditions
        if(TMath::Abs(etaTrack - c.Eta()) < fIsoConeRadius)
          sumpTPhiBandTrack += eTrack->Pt();
      }
      else{
        sumpTConeCharged += eTrack->Pt(); // should not double count if the track matching is already done
        iTracksCone++;
      }
    }
  }
  if(fWho==2){
    fTrackMultvsSumChargedvsUE->Fill(iTracksCone,sumpTConeCharged,sumpTPhiBandTrack);
    fTrackMultvsPt->Fill(iTracksCone,c.Pt());
  }
  ptIso = sumpTConeCharged;
  phiBandtrack = sumpTPhiBandTrack;
}


  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackEtaBand(TLorentzVector c, Double_t &ptIso, Double_t &etaBandtrack){
    // Underlying events study with tracks in eta band in EMCAL acceptance
  
    //INSERT BOUNDARIES ACCORDING TO THE FLAG WE WANT TO USE
  Double_t sumpTConeCharged=0., sumpTEtaBandTrack=0.;
  Double_t minPhi= 0., maxPhi= 2*TMath::Pi(), minEta = -0.87, maxEta= 0.87;
  
  if(!fTPC4Iso){
    minEta = -0.67;
    maxEta = 0.67;
    minPhi = 1.4;
    maxPhi = 3.15;
  }
  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  
    //  Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());
  tracksAna->ResetCurrentID();
  AliVTrack *eTrack = 0x0;
  Double_t phiTrack,etaTrack,radius;
  Int_t iTracksCone = 0;
  AliAODTrack *aodEtrack=0x0;
  while((eTrack = static_cast<AliVTrack*>(tracksAna->GetNextAcceptParticle()))){
    if(!eTrack){
      AliError(Form("No tracks in collection"));
      continue;
    }
    if(!fIsEsd){
      aodEtrack = static_cast<AliAODTrack*>(eTrack);
      if(!(aodEtrack->IsHybridGlobalConstrainedGlobal())) continue;
    }
    
    if((eTrack->Pt())<0.2)
      continue;
    
    if((eTrack->GetStatus() & AliAODTrack::kITSrefit)!=AliAODTrack::kITSrefit) continue;
    
    if(!fIsEsd){
      Double_t frac = 0;
      Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
      Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if (  ncls> 0 )  frac =  nclsS / ncls ;
      
      if(frac > 0.4) continue;
    }
    
      //CHECK IF TRACK IS IN BOUNDARIES
    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();
      // define the radius between the leading cluster and the considered cluster
      //redefine phi/c.Eta() from the cluster we passed to the function
    if( (phiTrack < maxPhi) && (phiTrack > minPhi) && (etaTrack < maxEta) && (etaTrack > minEta)){
      radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2));
      
      if(radius>fIsoConeRadius){ // if tracks are outside the isolation cone study
        
          // actually phi band here --- ADD Boundaries conditions
        if(TMath::Abs(phiTrack - c.Phi()) < fIsoConeRadius)
          sumpTEtaBandTrack += eTrack->Pt();
      }
      else
        sumpTConeCharged += eTrack->Pt(); // should not double count if the track matching is already done
      iTracksCone++;
    }
  }
  if(fWho==2){
    fTrackMultvsSumChargedvsUE->Fill(iTracksCone,sumpTConeCharged, sumpTEtaBandTrack,sumpTEtaBandTrack);
    fTrackMultvsPt->Fill(iTracksCone,c.Pt());
  }
  ptIso = sumpTConeCharged;
  etaBandtrack = sumpTEtaBandTrack;
}


  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackOrthCones(TLorentzVector c, Double_t &ptIso, Double_t &cones){
    // Underlying events study with tracks in orthogonal cones in TPC
  
  Double_t sumpTConeCharged=0., sumpTPerpConeTrack=0.;
  
  Double_t etaClus = c.Eta();
  Double_t phiClus = c.Phi();
  Double_t phiCone1 = phiClus - TMath::PiOver2();
  phiCone1=fmod(phiCone1,TMath::TwoPi());
  Double_t phiCone2 = phiClus + TMath::PiOver2();
  phiCone2=fmod(phiCone2,TMath::TwoPi());
  
  if (phiCone1 < 0.) phiCone1 += 2*TMath::Pi();
  
  AliTrackContainer *tracks = GetTrackContainer(0);
  
  if(!tracks){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  
    //  Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());
  tracks->ResetCurrentID();
  AliVTrack *eTrack = 0x0;
  AliAODTrack *aodEtrack=0x0;
  Double_t phiTrack,etaTrack,dist2Clust,dist2Cone1,dist2Cone2;
  Int_t iTracksCone = 0;
  
  while((eTrack = static_cast<AliVTrack*>(tracks->GetNextAcceptParticle()))){
    if(!eTrack){
      AliError(Form("No tracks in collection"));
      continue;
    }
    if(!fIsEsd){
      aodEtrack = static_cast<AliAODTrack*>(eTrack);
      if(!(aodEtrack->IsHybridGlobalConstrainedGlobal())) continue;
    }
    
    if((eTrack->Pt())<0.2)
      continue;
    
    if((eTrack->GetStatus() & AliAODTrack::kITSrefit)!=AliAODTrack::kITSrefit) continue;
    
    if(!fIsEsd){
      Double_t frac = 0;
      Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
      Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if (  ncls> 0 )  frac =  nclsS / ncls ;
      
      if(frac > 0.4) continue;
    }
    
    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();
    dist2Clust = TMath::Sqrt(TMath::Power(etaTrack-etaClus, 2)+TMath::Power(phiTrack-phiClus, 2));
    
    if (dist2Clust<fIsoConeRadius){
      sumpTConeCharged += eTrack->Pt(); // tracks are inside the isolation cone
      iTracksCone++;
    }
    
    else{//tracks outside the IsoCone
         //Distances from the centres of the two Orthogonal Cones
      dist2Cone1 = TMath::Sqrt(TMath::Power(etaTrack-etaClus, 2)+TMath::Power(phiTrack-phiCone1, 2));
      dist2Cone2 = TMath::Sqrt(TMath::Power(etaTrack-etaClus, 2)+TMath::Power(phiTrack-phiCone2, 2));
      
        //Is the Track Inside one of the two Cones ->Add to UE
      if((dist2Cone1 < fIsoConeRadius) || (dist2Cone2 < fIsoConeRadius)) sumpTPerpConeTrack += eTrack->Pt();
    }
  }
  if(fWho==2){
    fTrackMultvsSumChargedvsUE->Fill(iTracksCone,sumpTConeCharged,sumpTPerpConeTrack );
    fTrackMultvsPt->Fill(iTracksCone,c.Pt());
  }
  
  ptIso = sumpTConeCharged;
  cones = sumpTPerpConeTrack;
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackFullTPC(TLorentzVector c, Double_t &ptIso, Double_t &full){
    // Underlying events study with tracks in full TPC except back to back bands
  
  Double_t sumpTConeCharged=0., sumpTTPCexceptB2B=0.;
  
  AliTrackContainer *tracks = GetTrackContainer(0);
  
  if(!tracks){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  
    //  Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());
  tracks->ResetCurrentID();
  AliVTrack *eTrack = 0x0;
  Double_t phiTrack,etaTrack,radius, dphiUp, dphiDown;
  Int_t iTracksCone = 0;
  AliAODTrack *aodEtrack=0x0;
  while((eTrack = static_cast<AliVTrack*>(tracks->GetNextAcceptParticle()))){
    if(!eTrack){
      AliError(Form("No tracks in collection"));
      continue;
    }
    if(!fIsEsd){
      aodEtrack = static_cast<AliAODTrack*>(eTrack);
      if(!(aodEtrack->IsHybridGlobalConstrainedGlobal())) continue;
    }
    
    if((eTrack->Pt())<0.2)
      continue;
    
    if((eTrack->GetStatus() & AliAODTrack::kITSrefit)!=AliAODTrack::kITSrefit) continue;
    
    if(!fIsEsd){
      Double_t frac = 0;
      Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
      Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if (  ncls> 0 )  frac =  nclsS / ncls ;
      
      if(frac > 0.4) continue;
    }
    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();
      //redefine phi/c.Eta() from the cluster we passed to the function
    radius = TMath::Sqrt(TMath::Power(phiTrack-c.Phi(),2)+TMath::Power(etaTrack-c.Eta(),2)); // define the radius between the leading cluster and the considered cluster
    
    if(radius>fIsoConeRadius){ // if tracks are outside the isolation cone study UE
      dphiUp = c.Phi() + TMath::Pi() - fIsoConeRadius;
      dphiDown = c.Phi() + TMath::Pi() + fIsoConeRadius;
        // TPC except B2B
      if(phiTrack < dphiDown && phiTrack> dphiUp) sumpTTPCexceptB2B += eTrack->Pt();
    }
    else {
      sumpTConeCharged += eTrack->Pt(); // should not double count if the track matching is already done
      iTracksCone++;
    }
  }
  if(fWho==2){
    fTrackMultvsSumChargedvsUE->Fill(iTracksCone,sumpTConeCharged, sumpTTPCexceptB2B);
    fTrackMultvsPt->Fill(iTracksCone,c.Pt());
  }
  ptIso = sumpTConeCharged;
  full = sumpTTPCexceptB2B;
}

  //__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::CheckBoundaries(TLorentzVector vecCOI){
    // Check if the cone around the considered cluster is in EMCAL acceptance
    //AliInfo("Inside CheckBoundaries\n");
  
  Double_t minPhiBound= 0. , minEtaBound= 0., maxPhiBound= 0., maxEtaBound= 0.;
  Bool_t isINBoundaries;
  
  if(fTPC4Iso) // to avoid to have the isolation partially outside the TPC acceptance in eta
  {
    minPhiBound = 1.4;
    
    if(!fAnalysispPb)
      maxPhiBound = 2.740; // normally 110° but shorter cut to avoid EMCAL border
    else
      maxPhiBound = 3.15; // normally 110° but shorter cut to avoid EMCAL border
    
    minEtaBound = -0.87+fIsoConeRadius; // ""
    maxEtaBound = 0.87-fIsoConeRadius; // ""
  }
  else{
    minEtaBound = -0.67+fIsoConeRadius; // ""
    maxEtaBound = 0.67-fIsoConeRadius; // ""
    
    if(!fAnalysispPb){
      minPhiBound = 1.798;
      maxPhiBound = 2.740; // normally 110° but shorter cut to avoid EMCAL border
    }
    else{
      minPhiBound = 1.4+fIsoConeRadius;
      maxPhiBound = 3.15-fIsoConeRadius; // normally 110° but shorter cut to avoid EMCAL border
    }
  }
  
  if(vecCOI.Eta() > maxEtaBound || vecCOI.Eta() < minEtaBound || vecCOI.Phi() > maxPhiBound || vecCOI.Phi() <minPhiBound)
    isINBoundaries=kFALSE;
  else
    isINBoundaries=kTRUE;
  
  return isINBoundaries;
}
  //_________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::LookforParticle(Int_t clusterlabel, Double_t energyCLS, Double_t phiCLS, Double_t etaCLS, Double_t /*time*/, Double_t ss, Double_t isolation){
  
    //cout<<"\n\n\n\n\n\n\nInside Look4Particle \n For Cluster "<<clusterlabel<<"\t\t"<<energyCLS<<"\t\t"<<etaCLS<<"\t\t"<<phiCLS<<"\t\t"<<ss<<"\t\t"<<isolation<<"\n\n\n\n"<<endl;
  
  if (!fIsMC){
    AliWarning("not a montecarlo run!!!!!!");
    return;
  } //AliInfo(Form("It's a MC analysis %e",fAODMCParticles));
  if(!fStack && !fAODMCParticles){
    AliWarning("No Particle Stack !!!!!");
    //cout<<"No Particle Stack !!!!!"<<endl;
    return;
  }
  if(fAODMCParticles->GetEntries() < 1){
    AliWarning("number of tracks insufficient!!!!");
    return;
  }
  
  Int_t ndimsMCmix = fMCQAdim;
  
  Double_t outputvalueMCmix[ndimsMCmix];
    //cout<<"dimensions of the array: "<<ndimsMCmix<<endl;
  
  Int_t npart=fAODMCParticles->GetEntries();
    //cout<<"Number of particles in the event: "<<npart<<endl;
  
  AliAODMCParticle *particle2Check, *momP2Check;
  
  Int_t clustPDG, p2clabel, p2ccharge;
  Double_t enTrue,phiTrue, etaTrue;
  Double_t dPhi,dEta ;
  Int_t clusterFromPromptPhoton=-1;
  
  particle2Check = static_cast<AliAODMCParticle*>(fAODMCParticles->At(clusterlabel));
  clustPDG=particle2Check->GetPdgCode();
  int mom2checkidx = particle2Check->GetMother();
  momP2Check = static_cast<AliAODMCParticle*>(fAODMCParticles->At(mom2checkidx));
  
  phiTrue = particle2Check->Phi(); //
  etaTrue = particle2Check->Eta(); //Basic quantities from MCtruth.
  enTrue  = particle2Check->E()*TMath::Sin(particle2Check->Theta()); // Now let's check if we need corrections to the Energy.
  dPhi = phiCLS-phiTrue;
  dEta = etaCLS-etaTrue;
  p2ccharge = particle2Check->Charge();
    
  if(clustPDG==22 || (TMath::Abs(clustPDG) == 11 && momP2Check->GetPdgCode() == 22)){
    //Direct Photon or e+/e- coming from Photon Conversion
      //Checking if the Photon is a decay product of pi0 or eta meson. //Maybe include omega?
      //printf("Cluster Label: %d Asso. with a MCpar w/ PDG %d",clusterlabel,clustPDG);
    if(clustPDG==22){
        //printf("\twhose mother is a %d",momP2Check->GetPdgCode());
      if( momP2Check->GetPdgCode()!=22){
          //printf("  which is not a photon so the cluster is from a decay ");
        if (momP2Check->GetPdgCode()==111 || momP2Check->GetPdgCode()==221) {
            //printf(" of a pi0 or a eta mesons");
          clusterFromPromptPhoton=5;
          
          Int_t idxdaug1 = momP2Check->GetFirstDaughter();
          Int_t idxdaug2 = momP2Check->GetLastDaughter();
          if ( idxdaug1 == clusterlabel ){ //Cluster associated with the 1st daughter? Then look if also the 2nd daughter contributes to the cluster Energy
            if ( idxdaug2<npart ){//2nd daughter within List of Particles.
              
              AliAODMCParticle *daug2 = static_cast<AliAODMCParticle*>(fAODMCParticles->At(idxdaug2));
              if ( daug2->GetPdgCode() == 22 && ( daug2->Phi() - phiTrue )< 0.01 && ( daug2->Eta() - etaTrue )< 0.01 ){//very tight cut because they are photons
                enTrue += daug2->E()*TMath::Sin(daug2->Theta());
                clusterFromPromptPhoton=7;//contribution from both daughters
              }
              else
                clusterFromPromptPhoton=6;//contribution from one daughter
            }
          }
          else{//then the Cluster MUST BE associated with the 2nd daughter!! look if also the 1st daughter contributes to the cluster Energy
            
            AliAODMCParticle *daug1 = static_cast<AliAODMCParticle*>(fAODMCParticles->At(idxdaug1));
            if ( daug1->GetPdgCode() == 22 && ( daug1->Phi()-phiTrue )< 0.01 && ( daug1->Eta()-etaTrue )< 0.01 ){//very tight cut because they are photons
              enTrue += daug1->E()*TMath::Sin(daug1->Theta());
              clusterFromPromptPhoton=6;//contribution from both daughters
            }
            else
              clusterFromPromptPhoton=7;//contribution from one daughter
          }
          fpi0VSclusterVSIsolation->Fill(momP2Check->E()*TMath::Sin(momP2Check->Theta()), energyCLS, isolation);
          fpi0VSclusterVSM02->Fill(momP2Check->E()*TMath::Sin(momP2Check->Theta()), energyCLS,ss);
          fpi0VSM02VSIsolation->Fill(momP2Check->E()*TMath::Sin(momP2Check->Theta()), ss, isolation);
        }
        else{
            //printf("  of a non considered meson/baryon");
          clusterFromPromptPhoton=8;//Undefined
        }
      }
      else{
        clusterFromPromptPhoton=1;//TruePromptPhoton
                                  //printf("  so we have a prompt photon\n");
      }
    }
    else{//Cluster created by e+/e- from Photon Conversion
         //printf(" whose mother PDG is %d and occupies the stack position at: %d\n",momP2Check->GetPdgCode(),mom2checkidx);
      Int_t g_momindex = momP2Check->GetMother();
      AliAODMCParticle *gMomP2Check=static_cast<AliAODMCParticle*>(fAODMCParticles->At(g_momindex));
      if ( mom2checkidx == 8 || (g_momindex == 8 && gMomP2Check->GetPdgCode()==22 && momP2Check->GetPdgCode()==22)){
        clusterFromPromptPhoton=1; //e+/e- from Converted Prompt Photon
                                   //printf(" This means it is a e+/e- cluster from a Converted PromptPhoton.\n");
      }
      else{
          //printf(" This means it is a e+/e- cluster from a Converted DECAYPhoton.\n");
        clusterFromPromptPhoton=5; //cluster created by a Photon but not a prompt one
      }
      Int_t firstidx=momP2Check->GetFirstDaughter();
      Int_t lastidx=momP2Check->GetLastDaughter();
      if(clusterFromPromptPhoton==1){
        if( firstidx == clusterlabel ){//Cluster associated with the 1st electron? Then look if also the 2nd electron contributes to the cluster Energy
          if( lastidx < npart ){//2nd daughter within List of Particles.
            
            AliAODMCParticle *last=static_cast<AliAODMCParticle*>(fAODMCParticles->At(lastidx));
            if(( last->Phi() - phiTrue ) < fdphicut && ( last->Eta() - etaTrue ) < fdetacut ){ //same Proximity cut as the CPV
              enTrue += last->E()*TMath::Sin(last->Theta());
              clusterFromPromptPhoton=3; //contribution from both daughters
                                         //printf(" The cluster HAS actually contribution from both daughters but is asso to the second daughter (lower energy)\n");
            }
            else clusterFromPromptPhoton=2; //contribution from one daughter
          }
        }
        else{ //Cluster associated to 2nd daughter!! look if also the 1st daughter contributes to the cluster Energy
          
          AliAODMCParticle *first=static_cast<AliAODMCParticle*>(fAODMCParticles->At(firstidx));
          if(( first->Phi() - phiTrue ) < fdphicut && ( first->Eta() - etaTrue ) < fdetacut ){//same Proximity cut as the CPV
            enTrue += first->E()*TMath::Sin(first->Theta());
            clusterFromPromptPhoton=3;//contribution from both daughters
                                      //printf("cluster HAS actually contribution from both daughters but is asso to the first daughter (higher energy)\n");
          }
          else
            clusterFromPromptPhoton=2;//contribution from one daughter
        }
      }
      if(clusterFromPromptPhoton >= 5){//Check on wheter also the 2nd gamma from pi0/eta decay contributes to the Energy of the cluster
                                       //This Further check is implemented to take care of very asymmetric decays.
        Int_t idxgrandma = momP2Check->GetMother();
        AliAODMCParticle *grandma=static_cast<AliAODMCParticle*>(fAODMCParticles->At(idxgrandma));
        if( grandma->GetPdgCode() == 111 || momP2Check->GetPdgCode() == 221 ){ //Add also omega mesons, Lambda barion, neutral Kaons?
          
          Int_t idxaunt1 = grandma->GetFirstDaughter();
          Int_t idxaunt2 = grandma->GetLastDaughter();
          if( idxaunt1 ==  mom2checkidx ){ //the 1st daughter of the pi0/eta is the mother of the electron that produced the Cluster.
                                           //Check if the 2nd pi0/eta daughter contributes to the cluster energy.
            AliAODMCParticle *aunt=static_cast<AliAODMCParticle*>(fAODMCParticles->At(idxaunt2));
            if(( aunt->Phi() - phiTrue ) < 0.01 && ( aunt->Eta() - etaTrue ) < 0.01 ){
              enTrue += aunt->E()*TMath::Sin(aunt->Theta());
              clusterFromPromptPhoton=7;//contribution from both daughters
            }
            else
              clusterFromPromptPhoton=6;//contribution from one daughter
          }
          else{//the 2nd daughter of the pi0/eta is the mother of the electron that produced the Cluster.
               //Check if the 1st pi0/eta daughter contributes to the cluster energy.
            AliAODMCParticle *aunt=static_cast<AliAODMCParticle*>(fAODMCParticles->At(idxaunt1));
            if(( aunt->Phi() - phiTrue ) < 0.01 && ( aunt->Eta() - etaTrue ) < 0.01 ){
              enTrue += aunt->E()*TMath::Sin(aunt->Theta());
              clusterFromPromptPhoton=6;//contribution from both daughters
            }
            else
              clusterFromPromptPhoton=7;//contribution from one daughter
          }
          fpi0VSclusterVSIsolation->Fill(grandma->E()*TMath::Sin(grandma->Theta()), energyCLS, isolation);
          fpi0VSclusterVSM02->Fill(momP2Check->E()*TMath::Sin(momP2Check->Theta()), energyCLS, ss);
          fpi0VSM02VSIsolation->Fill(grandma->E()*TMath::Sin(grandma->Theta()), ss, isolation);

        }
        else
          clusterFromPromptPhoton=8;//Undefined
      }
    }
    
      //printf("\nCluster %d  PDG: %d  (Mom is a %d) with pT: %f  ",clusterlabel,clustPDG,momP2Check->GetPdgCode(),energyCLS);
      //printf(" with clusterFromPromptPhoton stored: %d for cluster w/label %d\n",clusterFromPromptPhoton,clusterlabel);
   
  }
  else{
      if(p2ccharge==0){
        clusterFromPromptPhoton=10;
      }
      else{
          clusterFromPromptPhoton=11;
      }
      //AliWarning(Form("Hadronic cluster with energy: %f, M02: %f, cluster PDG: %d, mother PDG: %d, truth energy: %f, dphi: %f, deta %f, isolation energy: %f, cluster label: %d, charge: %d",energyCLS, ss, clustPDG, momP2Check->GetPdgCode(), enTrue, dPhi, dEta, isolation, clusterFromPromptPhoton, p2ccharge ));
  }
  outputvalueMCmix[0] = energyCLS;
  outputvalueMCmix[1] = ss;
  outputvalueMCmix[2] = clustPDG;
  outputvalueMCmix[3] = momP2Check->GetPdgCode();
  outputvalueMCmix[4] = enTrue;
  outputvalueMCmix[5] = dPhi;
  outputvalueMCmix[6] = dEta;
  outputvalueMCmix[7] = isolation;
  outputvalueMCmix[8] = clusterFromPromptPhoton;
    //clusterFromPP=1 ->clusterlabel = 8 TruePromptPhoton;
    //clusterFromPP=2 ->clusterlabel = indexe+/e- with 1 contribution to the Energy;
    //clusterFromPP=3 ->clusterlabel = indexe+/e- with 2 contributions to the Energy;
    //clusterFromPP=6 -> clusterlabel= indexgamma1/2 (or e1e2e3e4) with contribution from max 2 electrons to the Energy;
    //clusterFromPP=7 -> clusterlabel= indexgamma1/2 (or e1e2e3e4) with 4 contribution to the energy;
    //clusterFromPP=8 -> clusterlabel= Gamma decay NOT from pi0/eta decay.
    //clusterFromPP=10 -> clusterlabel= Neutral hadronic particle.
    //clusterFromPP=11 -> clusterlabel= Charged hadronic particle.
  
  if(fWho==1) fOutClustMC->Fill(outputvalueMCmix);
  
  return;
}
  //__________________________________________________________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::FillInvMassHistograms(Bool_t iso, Double_t m02COI, TLorentzVector c, Int_t index, Double_t isolation){

  return;
  
  Double_t clustTOF=0, invMass=0;
  
  Double_t invMasspi0Min=0.135-0.035;
  Double_t invMasspi0Max=0.135+0.035;
  
  Double_t invMassetaMin=0.548-0.035;
  Double_t invMassetaMax=0.548+0.035;
  
    //AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliClusterContainer *clusters = GetClusterContainer(0);
  Int_t localIndex=0;
  for(auto it :clusters->accepted()){ //check the position of other clusters in respect to the leading cluster
    AliVCluster* coi = static_cast<AliVCluster*>(it);
    
    localIndex=coi->GetID();
    if(localIndex==index) continue;
    else{
      localIndex++;
      
      TLorentzVector nClust; //STILL NOT INITIALIZED
      coi->GetMomentum(nClust,fVertex);
      
        //   Float_t phiClust =nClust.Phi();
        //   Float_t etaClust= nClust.Eta();
        //    Float_t eTcluster=0;
      
      if((coi->GetNCells() < 2)) continue;
      if((coi->GetDistanceToBadChannel() < 2)) continue;
      
      clustTOF = coi->GetTOF()*1e9;
      if(!fIsMC)
        if(clustTOF<-30. || clustTOF>30.) continue;
      
      
        // the TM is commented out in order to avoid to reject the cluster partner if it is near to a track
        //    if(ClustTrackMatching(coi)) continue;
      
        //redefine phi/c.Eta() from the cluster we passed to the function
      
      invMass=(c+nClust).M();
      
      if(fWho==2){
        // if(iso)  fInvMassM02iso -> Fill(invMass,m02COI,c.Et());
        // if(!iso) fInvMassM02noiso -> Fill(invMass,m02COI,c.Et());
        
	// if(invMass>invMasspi0Min && invMass<invMasspi0Max){
        //   fPtvsM02vsSumPi0->Fill(c.Pt(),m02COI,isolation);
        //     //  AliError(Form("Passe dans la boucle des pi0"));
        //   break;
        // }
        
        
        // if(invMass>invMassetaMin && invMass<invMassetaMax){
        //   fPtvsM02vsSumEta->Fill(c.Pt(),m02COI,isolation);
        //     //   AliError(Form("Passe dans la boucle des pi0 et eta"));
        //   break;
        // }
      }
    }
  } // end of clusters loop
  
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::IsolationAndUEinEMCAL(AliVCluster *coi, Double_t& isolation,Double_t& ue,Double_t eTThreshold, Int_t index){
  
  Double_t isoConeArea = TMath::Pi()*fIsoConeRadius*fIsoConeRadius;
  
  Double_t etaBandArea = 1.4*2.*fIsoConeRadius-isoConeArea;
  Double_t phiBandArea = 1.75*2.*fIsoConeRadius-isoConeArea;
  
  TLorentzVector vecCOI;
  coi->GetMomentum(vecCOI,fVertex);
  Double_t m02COI=coi->GetM02();
  Double_t eTCOI=vecCOI.Et();
  
    //EMCAL Only for Acceptance of Cell/Clusters/or Tracks
  switch(fIsoMethod)
  {
    case 0: //EMCAL CELLS
      
      switch (fUEMethod)
    {
      case 0: //phi band
        EtIsoCellPhiBand(vecCOI, isolation, ue);
          //Normalisation ue wrt Area - TO DO-
        ue = ue * (isoConeArea / phiBandArea);
        
        if(fWho==2){
          fPhiBandUECells->Fill(vecCOI.Pt() , ue);
          fEtIsoCells->Fill(isolation);
        }
        if(isolation<eTThreshold)
        {
          if(fWho==2) fEtIsolatedCells->Fill(eTCOI);
          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
        break;
      case 1: //eta band
        EtIsoCellEtaBand(vecCOI, isolation, ue);
          //Normalisation ue wrt Area - TO DO-
        ue = ue * (isoConeArea / etaBandArea);
        if(fWho==2){
          fEtaBandUECells->Fill(vecCOI.Pt() , ue);
          fEtIsoCells->Fill(isolation);
        }
        if(isolation<eTThreshold)
        {
          if(fWho==2) fEtIsolatedCells->Fill(eTCOI);
          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
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
        break;
      case 1: //eta band
        EtIsoClusEtaBand(vecCOI, isolation, ue,index);
          //Normalisation ue wrt Area - TO DO-
        ue = ue * (isoConeArea / etaBandArea);
        break;
    }
      if(fWho==2){
        fPtvsM02vsSum->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
          //AliError(Form("Passe bien Fill General histograms"));
        isolation=isolation-ue;  // ue subscraction
        fPtvsM02vsSumUE->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        fEtIsoClust->Fill(vecCOI.Pt(),isolation);
      }
      if(isolation<eTThreshold)
      {
        FillInvMassHistograms(kTRUE, m02COI, vecCOI, index, isolation);
          //	      AliError(Form("rempli histo"));
        if(fWho==2){
          fPtvsM02iso->Fill(vecCOI.Pt(),coi->GetM02());
          fPtIsolatedNClust->Fill(vecCOI.Pt());
        }
        fPtisoT=vecCOI.Pt();
        fM02isoT=m02COI;
        
        if(fM02mincut < m02COI && m02COI < fM02maxcut)
        {
          if(fWho==2) fEtIsolatedClust->Fill(eTCOI);
          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
      }
      else
      {
        if (isolation>3.) FillInvMassHistograms(kFALSE, m02COI, vecCOI, index, isolation);
        
        if(fWho==2) fPtvsM02noiso->Fill(vecCOI.Pt(),coi->GetM02());
        fPtnoisoT=vecCOI.Pt();
        fM02noisoT=m02COI;
      }
      break;
      
    case 2: //TRACKS (TBD which tracks) //EMCAL tracks
      switch (fUEMethod)
    {
      case 0: //phi band
        PtIsoTrackPhiBand(vecCOI, isolation, ue);
          //Normalisation ue wrt Area - TO DO-
        if(fWho==2) fPhiBandUETracks->Fill(vecCOI.Pt() , ue);
        ue = ue * (isoConeArea / phiBandArea);
      case 1: //eta band
        PtIsoTrackEtaBand(vecCOI, isolation, ue);
          //Normalisation ue wrt Area - TO DO-
        if(fWho==2) fEtaBandUETracks->Fill(vecCOI.Pt() , ue);
        ue = ue * (isoConeArea / etaBandArea);
        break;
          // case 2: //Cones
          // PtIsoTrackOrthCones(vecCOI, absId, isolation, ue);
          // break;
          // case 3: //Full TPC
          // PtIsoTrackFullTPC(vecCOI, absId, isolation, ue);
          // break;
    }
        // Fill histograms for isolation
      if(fWho==2){
        fPtvsM02vsSum->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        
        isolation=isolation-ue;  // ue subscraction
        fPtvsM02vsSumUE->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        
        fPtIsoTrack->Fill(vecCOI.Pt() , isolation);
      }
      if(isolation<eTThreshold)
      {
        FillInvMassHistograms(kTRUE, m02COI, vecCOI, index, isolation);
        if(fWho==2) {
          fPtvsM02iso->Fill(vecCOI.Pt(),coi->GetM02());
          fPtIsolatedNTracks->Fill(vecCOI.Pt());
        }
        fPtisoT=vecCOI.Pt();
        fM02isoT=m02COI;
        
        if(fM02mincut < m02COI && m02COI < fM02maxcut)
        {
          if(fWho==2) fEtIsolatedTracks->Fill(eTCOI);
          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
      }
      else
      {
        if(isolation>3.)  FillInvMassHistograms(kFALSE, m02COI, vecCOI, index, isolation);
        
        if(fWho==2) fPtvsM02noiso->Fill(vecCOI.Pt(),coi->GetM02());
        fPtnoisoT=vecCOI.Pt();
        fM02noisoT=m02COI;
      }
      break;
  }
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::IsolationAndUEinTPC(AliVCluster *coi, Double_t& isolation,Double_t& ue,Double_t eTThreshold, Int_t index){
  
    //EMCAL + TPC (Only tracks for the Isolation since IsoCone Goes Out of EMCAL)
  
  Double_t isoConeArea = TMath::Pi()*fIsoConeRadius*fIsoConeRadius;
  
  Double_t etaBandAreaTr = 1.74*2.*fIsoConeRadius-isoConeArea;
  Double_t phiBandAreaTr = 2.*TMath::Pi()*2.*fIsoConeRadius-isoConeArea;
  Double_t perpConesArea = 2.*isoConeArea;
  Double_t fullTPCArea = (1.8*2.*TMath::Pi())-(fIsoConeRadius*2*1.8)-isoConeArea;
  
  TLorentzVector vecCOI;
  coi->GetMomentum(vecCOI,fVertex);
  Double_t m02COI=coi->GetM02();
  Double_t eTCOI=vecCOI.Et();
  
  switch (fUEMethod)
  {
    case 0: //phi band
      PtIsoTrackPhiBand(vecCOI, isolation, ue);
        //Normalisation ue wrt Area - TO DO-
        //      ue = ue * (isoConeArea / phiBandAreaTr);
      fPhiBandUETracks->Fill(vecCOI.Pt() , ue);
      ue = ue * (isoConeArea / phiBandAreaTr);
        // fill histograms for isolation
      if(fWho==2) fPtvsM02vsSum->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
      
      if(fAnalysispPb) isolation=isolation-ue;  // ue subscraction
      
      if(fWho==2) {
        fPtvsM02vsSumUE->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        fPtIsoTrack->Fill(vecCOI.Pt(), isolation);
      }
      
      if(isolation<eTThreshold)
      {
        FillInvMassHistograms(kTRUE, m02COI, vecCOI, index,isolation);
        if(fWho==2) {
          fPtvsM02iso->Fill(vecCOI.Pt(),coi->GetM02());
          fPtIsolatedNTracks->Fill(vecCOI.Pt());
        }
        fPtisoT=vecCOI.Pt();
        fM02isoT=m02COI;
        
        if(fM02mincut < m02COI && m02COI < fM02maxcut)
        {
          if(fWho==2) fEtIsolatedTracks->Fill(eTCOI);
          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
      }
      else
      {
        if(isolation>3.)   FillInvMassHistograms(kFALSE, m02COI, vecCOI, index,isolation);
        if(fWho==2) fPtvsM02noiso->Fill(vecCOI.Pt(),coi->GetM02());
        fPtnoisoT=vecCOI.Pt();
        fM02noisoT=m02COI;
      }
      break;
    case 1: //eta band
      PtIsoTrackEtaBand(vecCOI, isolation, ue);
        //Normalisation ue wrt Area - TO DO-
        //      ue = ue * (isoConeArea / etaBandAreaTr);
      if(fWho==2) fEtaBandUETracks->Fill(vecCOI.Pt() , ue);
      ue = ue * (isoConeArea / etaBandAreaTr);
        // fill histograms for isolation
      if(fWho==2) fPtvsM02vsSum->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
      
      if(fAnalysispPb) isolation=isolation-ue;  // ue subscraction
      
      if(fWho==2) {
        fPtvsM02vsSumUE->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        fPtIsoTrack->Fill(vecCOI.Pt(), isolation);
      }
      if(isolation<eTThreshold)
      {
        FillInvMassHistograms(kTRUE, m02COI, vecCOI, index,isolation);
        if(fWho==2) {
          fPtvsM02iso->Fill(vecCOI.Pt(),coi->GetM02());
          fPtIsolatedNTracks->Fill(vecCOI.Pt());
        }
        fPtisoT=vecCOI.Pt();
        fM02isoT=m02COI;
        
        if(fM02mincut < m02COI && m02COI < fM02maxcut)
        {
          if(fWho==2) fEtIsolatedTracks->Fill(eTCOI);
          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
      }
      else
      {
        if(isolation>3.)   FillInvMassHistograms(kFALSE, m02COI, vecCOI, index,isolation);
        if(fWho==2) fPtvsM02noiso->Fill(vecCOI.Pt(),coi->GetM02());
        fPtnoisoT=vecCOI.Pt();
        fM02noisoT=m02COI;
      }
      break;
    case 2: //Cones
      PtIsoTrackOrthCones(vecCOI, isolation, ue);
        //Normalisation ue wrt Area - TO DO-
        //      ue = ue * (isoConeArea / perpConesArea);
      if(fWho==2) fPerpConesUETracks ->Fill(vecCOI.Pt() , ue);
      ue = ue * (isoConeArea / perpConesArea);
        // fill histograms for isolation
      if(fWho==2) fPtvsM02vsSum->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
      
      if(fAnalysispPb) isolation=isolation-ue;  // ue subscraction
      
      if(fWho==2) {
        fPtvsM02vsSumUE->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        fPtIsoTrack->Fill(vecCOI.Pt(), isolation);
      }
      if(isolation<eTThreshold)
      {
        FillInvMassHistograms(kTRUE, m02COI, vecCOI, index, isolation);
        if(fWho==2) {
          fPtvsM02iso->Fill(vecCOI.Pt(),coi->GetM02());
          fPtIsolatedNTracks->Fill(vecCOI.Pt());
        }
        fPtisoT=vecCOI.Pt();
        fM02isoT=m02COI;
        
        if(fM02mincut < m02COI && m02COI < fM02maxcut)
        {
          if(fWho==2) fEtIsolatedTracks->Fill(eTCOI);
          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
      }
      else
      {
        if(isolation>3.)  FillInvMassHistograms(kFALSE, m02COI, vecCOI, index, isolation);
        if(fWho==2) fPtvsM02noiso->Fill(vecCOI.Pt(),coi->GetM02());
        fPtnoisoT=vecCOI.Pt();
        fM02noisoT=m02COI;
      }
      break;
    case 3: //Full TPC
      PtIsoTrackFullTPC(vecCOI, isolation, ue);
        //Normalisation ue wrt Area - TO DO-
        //      ue = ue * (isoConeArea / fullTPCArea);
      if(fWho==2) fTPCWithoutIsoConeB2BbandUE->Fill(vecCOI.Pt() , ue);
      ue = ue * (isoConeArea / fullTPCArea);
        // fill histograms for isolation
      
      if(fWho==2) fPtvsM02vsSum->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
      
      if(fAnalysispPb) isolation=isolation-ue;  // ue subscraction

      if(fWho==2) {
        fPtvsM02vsSumUE->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        fPtIsoTrack->Fill(vecCOI.Pt(), isolation);
      } //	fTracksConeEtaPt->Fill(isolation, vecCOI.Eta(), vecCOI.Pt());
        //	fTracksConeEtaM02->Fill(isolation, vecCOI.Eta(), coi->GetM02());
      if(isolation<eTThreshold)
      {
        FillInvMassHistograms(kTRUE, m02COI, vecCOI, index, isolation);
        if(fWho==2) {
          fPtvsM02iso->Fill(vecCOI.Pt(),coi->GetM02());
          fPtIsolatedNTracks->Fill(vecCOI.Pt());
        }
        fPtisoT=vecCOI.Pt();
        fM02isoT=m02COI;
        
        if(fM02mincut < m02COI && m02COI < fM02maxcut)
        {
          if(fWho==2) fEtIsolatedTracks->Fill(eTCOI);
          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
      }
      else
      {
        if(isolation>3.) FillInvMassHistograms(kFALSE, m02COI, vecCOI, index, isolation);
        if(fWho==2) fPtvsM02noiso->Fill(vecCOI.Pt(),coi->GetM02());
        fPtnoisoT=vecCOI.Pt();
        fM02noisoT=m02COI;
      }
      break;
  }
  
}

  //__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::FillGeneralHistograms(AliVCluster *coi, TLorentzVector vecCOI, Int_t index){
    //AliInfo("Inside FillGeneralHistograms\n");
    //    Printf("Inside FillGeneralHistograms\n");
    // Fill the histograms for underlying events and isolation studies
    // AliError(Form("Arrive bien dans fill general histograms"));
  
    // I would like to remove this part and fill the tracks multiplicity histogram in FillQAHistograms, is that ok for thnSparses? (especially cause here the histogram is filled several times per event)
    // AliParticleContainer *tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));
  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
    //Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());
  tracksAna->ResetCurrentID();
  
  const Int_t nTracks = tracksAna->GetNAcceptedTracks();
    //  Printf("Ntracks for the event with this cluster: %d", nTracks);
  if(fQA) fTrackMult->Fill(nTracks);
  
    //  Printf("After Loop on Tracks");
  Double_t eTCOI = 0., m02COI = 0.;
  
    //Definition of the Array for Davide's Output
  const Int_t ndims =   fNDimensions;
  Double_t outputValues[ndims];
  
  eTCOI = vecCOI.Et();
  
  Int_t nlm=0;
  AliVCaloCells * fCaloCells =InputEvent()->GetEMCALCells();
  if(fCaloCells)
    nlm = GetNLM(coi,fCaloCells);
  
    //  Printf("cluster ID %d with nlm %d . M02 BEFORE possible smearing %.4lf . Do we set the smearing ? %s ",coi->GetID(),nlm, m02COI = coi->GetM02(), fSSsmearing? "Yes":"No");
  
  if(fSSsmearing){
    if(coi->GetM02()>0.1){
        //      Printf("Smearing for only clusters with nlm = %d" ,fWhich);
      if(nlm==1){
        if((fSSsmearwidth != 0.)){
          TRandom3 *ran=new TRandom3(0);
          if(fWhich==0){//Landau Smearing
            Float_t smear = ran->Landau(fSSsmear_mean,fSSsmearwidth);
            if(fSSsmear_mean==0 || (fSSsmear_mean!=0 && coi->GetID()%3==0))
              m02COI = coi->GetM02() + smear;
          }
          else{//Gaussian Smearing
            Float_t smear = ran->Gaus(fSSsmear_mean,fSSsmearwidth);
            if(fSSsmear_mean==0 || (fSSsmear_mean!=0 && coi->GetID()%3==0))
              m02COI = coi->GetM02() + smear;
          }
        }
        else {
          AliWarning("The Smearing is set but the width of the distribution is null!\nNOT DOING ANYTHING for the Shower Shape!");
          m02COI = coi->GetM02();
        }
      }
      else
        m02COI = coi->GetM02();
    }
  }
  else{
    AliWarning("Smearing not SET!");
    m02COI = coi->GetM02();
  }

  
    //  Printf("M02 AFTER possible smearing: %.4lf",m02COI);
  
    //AliInfo(Form("M02 value: %lf\n",m02COI));
  
    // ******** Isolation and UE calculation with different methods *********
  
  Double_t eTThreshold = 5;
  
  switch(fEtIsoMethod)
  {
    case 0:  // SumEt<EtThr
      eTThreshold = fEtIsoThreshold;
      break;
      
    case 1:  // SumEt<%Ephoton
      eTThreshold = fEtIsoThreshold * eTCOI;
      break;
      
    case 2: // Etmax<eTThreshold
      eTThreshold = fEtIsoThreshold;
      if( eTCOI<eTThreshold && fWho==2) // photon candidate, cuts have to be decided after studies
      {
        fEtIsolatedClust->Fill(eTCOI);
      }
      break;
  }
  
  Double_t isolation=0, ue=0;
  
  if(!fTPC4Iso)
    IsolationAndUEinEMCAL(coi,isolation,ue,eTThreshold,index);
  else
    IsolationAndUEinTPC(coi,isolation,ue,eTThreshold,index);
  
  if(fIsMC)
    LookforParticle(coi->GetLabel(),vecCOI.Et(),vecCOI.Phi(),vecCOI.Eta(),coi->GetTOF()*1e9,m02COI,isolation);
  
  /*  Here we should call something to know the number of tracks...
   Soon I'll put in this version the "old way"; please let me know if
   any of you could do the same with the JET framework*/
  
  switch(fWho) {
    case 0:
      flambda0T=m02COI; // for all neutral clusters
      fEtT=vecCOI.Et(); // for all neutral clusters
      fPtT=vecCOI.Pt(); // for all neutral clusters
      fetaT=vecCOI.Eta(); // for all neutral clusters
      fphiT=vecCOI.Phi(); //for all neutral clusters
      fsumEtisoconeT=isolation;
        //	   AliError(Form("lambda 0 %f",flambda0T));
      fsumEtUE=ue;
        // does not fill the tree anyway files too big
        //    fOutputTree->Fill();
      break;
      
    case 1:
      outputValues[0] = eTCOI;
      outputValues[1] = m02COI;
      outputValues[2] = isolation;
      outputValues[3] = ue;
      outputValues[4] = vecCOI.Eta();
      outputValues[5] = vecCOI.Phi();
      
      fOutputTHnS -> Fill(outputValues);
      break;
      
  }
  return kTRUE;
}

  //_________________________________________________________________________

void AliAnalysisTaskEMCALPhotonIsolation::AddParticleToUEMC(Double_t& sumUE,AliAODMCParticle* mcpp, Double_t eta,Double_t phi){
  
  Double_t etap=mcpp->Eta();
  Double_t phip=mcpp->Phi();
  
  if(!fTPC4Iso){
    if(TMath::Abs(etap)>=0.7 || (phip<=1.4 || phip>= TMath::Pi()))
      return;
    else{
      switch(fUEMethod){
        case 0: //Phi band
          if(TMath::Abs(eta-etap)<fIsoConeRadius) //to be changed with fIsoConeRadius
            sumUE += mcpp->E()*TMath::Sin(mcpp->Theta());
          else
            return;
          
          break;
        case 1: //Eta band
          if(TMath::Abs(phi-phip)<fIsoConeRadius)
            sumUE += mcpp->E()*TMath::Sin(mcpp->Theta());
          else
            return;
          
          break;
      }
    }
  }
  else{
    if(TMath::Abs(etap)>=1.0)
      return;
    else{
      switch(fUEMethod){
        case 0: //Phi band
        {if(TMath::Abs(eta-etap)<fIsoConeRadius) //to be changed with fIsoConeRadius
          sumUE += mcpp->E()*TMath::Sin(mcpp->Theta());
        else
          return;
          break;
        }
        case 1: //Eta band
        {  if(TMath::Abs(phi-phip)<fIsoConeRadius)
          sumUE += mcpp->E()*TMath::Sin(mcpp->Theta());
        else
          return;
          
          break;
        }
        case 2: //Orthogonal Cones
        { double etacone1= eta;
          double etacone2= eta;
          double phicone1= phi - TMath::PiOver2();
          double phicone2= phi + TMath::PiOver2();
          
          if (phicone1 < 0.) phicone1 += 2*TMath::Pi();
          
          if(TMath::Sqrt(TMath::Power(etap-etacone1,2)+TMath::Power(phip-phicone1,2))< fIsoConeRadius ||
             TMath::Sqrt(TMath::Power(etap-etacone2,2)+TMath::Power(phip-phicone2,2))< fIsoConeRadius) //to be changed with fIsoConeRadius
          {sumUE += mcpp->Pt();}
          else
            return;
          
          break;
        }
        case 3: //Full TPC
        {    //                  Double_t phiup= phi +TMath::Pi()+fIsoConeRadius;
             //                  Double_t phidown= phi +TMath::Pi()-fIsoConeRadius;
             //
             //                  if(phip < phidown || phip > phiup ) //TO BE CHECKED
             //                    continue;
          break;
        }
      }
    }
  }
}
  //_________________________________________________________________________

void AliAnalysisTaskEMCALPhotonIsolation::CalculateUEDensityMC(Double_t& sumUE){
  
  Double_t isoConeArea = TMath::Pi()*fIsoConeRadius*fIsoConeRadius;
  
  Double_t etaBandArea = (1.4 - 2.*fIsoConeRadius)*2.*fIsoConeRadius;
  Double_t phiBandArea = ((5./9.)*TMath::Pi()- 2.*fIsoConeRadius)*2.*fIsoConeRadius;
  Double_t etaBandAreaTr = (1.8 - 2.*fIsoConeRadius)*2.*fIsoConeRadius;
  Double_t phiBandAreaTr = (2.*TMath::Pi() - 4.*fIsoConeRadius)*2.*fIsoConeRadius;//there is a 2 more because of the JET CONE B2B
  Double_t perpConesArea = 2.*isoConeArea;
  Double_t fullTPCArea = 1.8*2.*TMath::Pi()-fIsoConeRadius*(TMath::Pi()*fIsoConeRadius + 3.6);
  
  if(!fTPC4Iso){
    switch (fUEMethod){
      case 0:
        sumUE = sumUE * (isoConeArea / phiBandArea);
        break;
      case 1:
        sumUE = sumUE * (isoConeArea / etaBandArea);
        break;
    }
  }
  else{
    switch (fUEMethod){
      case 0:
        sumUE = sumUE * (isoConeArea / phiBandAreaTr);
        break;
      case 1:
        sumUE = sumUE * (isoConeArea / etaBandAreaTr);
        break;
      case 2:
        sumUE = sumUE * (isoConeArea / perpConesArea);
        break;
      case 3:
        sumUE = sumUE * (isoConeArea / fullTPCArea);
        break;
    }
  }
}

  //_________________________________________________________________________

void AliAnalysisTaskEMCALPhotonIsolation::AnalyzeMC(){
  
    //printf("New Event...Analysing Stack!");
  if (!fIsMC)
    return;
    //AliInfo(Form("It's a MC analysis %e",fAODMCParticles));
  if(!fStack && !fAODMCParticles)
  {cout<<"no stack saved\n"; return;}
  
    //cout<<"there's a List of particles"<<endl;
    //DO THIS ALSO FOR ESDs
  
  Double_t eT, sumEiso, sumUE,phi, eta, distance, phip, etap, mcfirstEnergy;
  
  if(fAODMCParticles->GetEntries() < 1){
    AliError("number of tracks insufficient");
    return;
  }
  int nDimMC = fMCDimensions;
  Double_t outputValuesMC[nDimMC];
  
  Int_t nTracks = fAODMCParticles->GetEntriesFast();
  Int_t nFSParticles = 0;
  AliAODMCParticle *multTracks;
  
  for(int a=0; a<nTracks; a++){
    
    multTracks = static_cast<AliAODMCParticle*>(fAODMCParticles->At(a));
    
    if(multTracks->IsPrimary() && multTracks->IsPhysicalPrimary() && multTracks->GetStatus()<10){
      if(TMath::Abs(multTracks->Eta())<=0.9 && multTracks->Charge()!= 0)
        nFSParticles++;
      else
        continue;
    }//implement final state particle condition
    else
      continue;
  }
    //AliInfo(Form("number of particles in the array %d",nTracks));
  AliAODMCParticle *mcpart, *mom, *mcpp,*mcsearch, *mcfirst, *mcfirstmom,*matchingtrack, *mum;
  
    // Bool_t prompt=kFALSE;
  Double_t mcEnergy, maxE, energy;
  Int_t pdg, mompdg, photonlabel;
  Double_t mcFirstEta=0., mcFirstPhi=0.;
  
    // AliAODMCParticle *mcfirst = static_cast<AliAODMCParticle*>(fAODMCParticles->At(0));
    //AliAODMCParticle *mcp, *mcpmaxE, *mcpp, *mom;
  if(!fisLCAnalysis){
      //Loop on the event
    for(int iTr=0;iTr<nTracks;iTr++){
      
      mcEnergy=0.;energy =0;
      eT=0.; phi=0.; eta=0.;
      
      mcpart = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTr));
      
      if(mcpart->GetStatus()>10) {continue;}
      if(!mcpart->IsPrimary()) {continue;}
      if(!mcpart->IsPhysicalPrimary()) {continue;}
      
      pdg = mcpart->GetPdgCode();
      
      if(pdg != 22 /*|| mcpart->GetLabel()!=8*/ )
      {continue;}
      
      eta = mcpart->Eta();
      phi = mcpart->Phi();
      
        //check photons in EMCAL //to be redefined with fIsoConeR
      if(!fTPC4Iso){
      if((TMath::Abs(eta)>0.67-fIsoConeRadius ) || (phi < 1.398 + fIsoConeRadius || phi>(TMath::Pi()-fIsoConeRadius-0.03)))
        continue;
      }
      else{
        if((TMath::Abs(eta)>0.87-fIsoConeRadius ) || (phi < 1.398 || phi>(TMath::Pi()-0.03)))
          continue;
      }
        //printf("\nParticle Position %d  and Label: %d  PDG: %d  Pt: %f  Eta: %f  Phi: %f",iTr, mcpart->GetLabel(),pdg,mcpart->Pt(), eta, phi);
      
      photonlabel = iTr;
      int momidx = mcpart->GetMother();
      
      mom = static_cast<AliAODMCParticle*>(fAODMCParticles->At(momidx));
      mompdg= TMath::Abs(mom->GetPdgCode());
      
        //printf("With Mother at %d with label %d which is a %d",momidx, mom->GetLabel(), mompdg);
      
      eT= mcpart->E()*TMath::Sin(mcpart->Theta()); //transform to transverse Energy
      
      fphietaPhotons->Fill(eta,phi,eT);
      
      bool foundmatch=kFALSE;
      for(int m=0;m<nTracks && foundmatch==kFALSE;m++){
        if(m==iTr) continue;
        
        matchingtrack = static_cast<AliAODMCParticle*>(fAODMCParticles->At(m));
        
        if(! matchingtrack->IsPrimary()) continue;
        if(! matchingtrack->IsPhysicalPrimary()) continue;
        if(matchingtrack->GetStatus()> 10 ) continue;
        
        Double_t etamatching = matchingtrack->Eta();
        Double_t phimatching = matchingtrack->Phi();
        
        if(TMath::Abs(eta-etamatching)<=fdetacut && TMath::Abs(phi-phimatching)<=fdphicut){
          foundmatch=kTRUE;
          fphietaOthers->Fill(matchingtrack->Eta(),matchingtrack->Phi(),eT);
          fphietaOthersBis->Fill(matchingtrack->Eta(),matchingtrack->Phi(),matchingtrack->Pt());
        }
      }
      
      if(foundmatch) continue;
      
      distance=0.;
      phip=0., etap=0.;
      sumEiso=0.,sumUE=0.;
      
      for(int iTrack=0;iTrack<nTracks;iTrack++){
        
        if(iTrack==photonlabel)
          continue;
        
        mcpp = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));
        
        if(!mcpp) {continue;}
        
        if(mcpp->Charge()!=0 && mcpp->GetStatus()>10)
          fPtTracksVSpTNC_MC->Fill(eT,mcpp->Pt());
        
        if(fIsoMethod==2)
          if((mcpp->Charge())==0) {continue;}
        
        if(mcpp->GetStatus()>10) {continue;}
        
        int mumidx=mcpp->GetMother();
        if (mumidx<0 || mumidx>nTracks) continue;
        
        mum = static_cast<AliAODMCParticle*>(fAODMCParticles->At(mumidx));
        if(mumidx == photonlabel || mum->GetPdgCode()==22) continue;
        
        phip = mcpp->Phi();
        etap = mcpp->Eta();
        
          //Depending on which Isolation method and UE method is considered.
        distance= TMath::Sqrt((phi-phip)*(phi-phip) + (eta-etap)*(eta-etap));
        
        if(distance <= fIsoConeRadius){
            //cout<<iTrack<<"\t"<<photonlabel<<endl;
            //mcpp->Print();
          sumEiso += mcpp->E()*TMath::Sin(mcpp->Theta());
          
        }
        else{
          AddParticleToUEMC(sumUE,mcpp, eta, phi);}
      }
      CalculateUEDensityMC(sumUE);
      
        //printf("Storing Particle: Label %d  PDG: %d  Eta: %f  Phi: %f",mcpart->GetLabel(),pdg,eta,phi);
        //printf("With Mother at %d with label %d which is a %d",momidx, mom->GetLabel(), mompdg);
      outputValuesMC[0] = eT;
      outputValuesMC[1] = sumEiso;
      outputValuesMC[2] = sumUE;
      outputValuesMC[3] = mompdg;
      outputValuesMC[4] = eta;
      outputValuesMC[5] = phi;
      outputValuesMC[6] = mcpart->GetLabel();
        // EtaPhiMCPhoton
        // EtMC
        // EtIsoCone
        // EtMother
        // UE Et
        // Mother PDG
        //fill some histograms or a THnSparse or a TTree.
        //	AliError(Form("Fill something in Analize MC"));
      if(fWho==1)  fOutMCTruth -> Fill(outputValuesMC);
    }
  }
  else{
    maxE=0.;
    int indexmaxE=0;
      //getting the index of the particle with the maximum energy.
    for(int iTr=0;iTr<nTracks;iTr++){
      mcsearch= static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTr));
      
      if(!mcsearch) continue;
      
      if(mcsearch->GetStatus()>10) continue;
      if(mcsearch->GetPdgCode()!=22) continue;
      if(TMath::Abs(mcsearch->Eta())>0.67-fIsoConeRadius) continue;
      if(mcsearch->Phi()<= 1.798 ||mcsearch->Phi()>= TMath::Pi()) continue;
      
      mcfirstEnergy= mcsearch->E();
      if(mcfirstEnergy>maxE){
        maxE=mcfirstEnergy;
        indexmaxE=iTr;
      }
      else continue;
    }
    mcfirst= static_cast<AliAODMCParticle*>(fAODMCParticles->At(indexmaxE));
    mcfirstEnergy=mcfirst->E()*TMath::Sin(mcfirst->Theta());
    
    int momidx= mcfirst->GetMother();
    mcfirstmom =  static_cast<AliAODMCParticle*>(fAODMCParticles->At(momidx));
    mompdg= TMath::Abs(mcfirstmom->GetPdgCode());
    mcFirstEta = mcfirst->Eta();
    mcFirstPhi = mcfirst->Phi();
    
    phip=0., etap=0.;
    sumEiso=0,sumUE=0;
    
    for(Int_t iTrack=1;iTrack<nTracks ;iTrack++){
      if(iTrack==indexmaxE) continue;
      
      mcpp= static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));
      phip = mcpp->Phi();
      etap = mcpp->Eta();
      if(!mcpp)
        continue;
      
      if(mcpp->GetStatus()>10) continue;
      if(mcpp->GetPdgCode()==22)continue;
      if(TMath::Abs(etap>0.7)) continue;
      if(phip<=1.4 || phip>= TMath::Pi()) continue;
      distance=0.;
      distance= TMath::Sqrt((mcFirstPhi- phip)*(mcFirstPhi- phip) + (mcFirstEta- etap)*(mcFirstEta- etap));
      
      if(distance<=fIsoConeRadius){
        sumEiso += mcpp->E()*TMath::Sin(mcpp->Theta());
      }
      else{
        AddParticleToUEMC(sumUE,mcpp,mcFirstEta,mcFirstPhi);
      }
    }
      //  cout<<"\n\nTotal Energy inside the Isolation Cone : "<<sumEiso<<endl;
    CalculateUEDensityMC(sumUE);
      //cout<<"Total UE Energy : "<<sumUE<<" calculated with method "<<fUEMethod<<endl;
    
      //Fill the Output TTree for MC Truth
  }
  
  return;
}
