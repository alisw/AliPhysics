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
// fParticleCollArray(),
fAOD(0),
fVevent(0),
fNCluster(0),
fAODMCParticles(0),
fmcHeader(0),
fTracksAna(0),
fStack(0),
fEMCALRecoUtils(new AliEMCALRecoUtils),
fWho(-1),
fSSsmearing(0),
fSSsmearwidth(0),
fSSsmear_mean(0),
fWhich(0),
// fOutputList(0),
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
fMCtruth(0),
fPeriod(""),
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
fNLM2_NC_Acc(0),
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
fPtvsM02vsSumUE(0),
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
fEtVSM02VSEisoclust(0),
fEtVSM02VSPisotrack(0),
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
// tracks(0),
// clusters(0)
{
  
    // Default constructor
    
    // fParticleCollArray.SetOwner(kTRUE);
    // for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
  
  SetMakeGeneralHistograms(kTRUE);
}

  //________________________________________________________________________
AliAnalysisTaskEMCALPhotonIsolation::AliAnalysisTaskEMCALPhotonIsolation(const char *name, Bool_t histo) :
AliAnalysisTaskEmcal(name, histo),
// fParticleCollArray(),
fAOD(0),
fVevent(0),
fNCluster(0),
fAODMCParticles(0),
fmcHeader(0),
fTracksAna(0),
fStack(0),
fEMCALRecoUtils(new AliEMCALRecoUtils),
fWho(-1),
fSSsmearing(0),
fSSsmearwidth(0),
fSSsmear_mean(0),
fWhich(0),
// fOutputList(0),
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
fMCtruth(0),
fPeriod(""),
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
fNLM2_NC_Acc(0),
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
fPtvsM02vsSumUE(0),
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
fEtVSM02VSEisoclust(0),
fEtVSM02VSPisotrack(0),
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
// tracks(0),
// clusters(0)
{
  
    // Standard constructor
  
    // fParticleCollArray.SetOwner(kTRUE);
    //  for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
  
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
  // printf("Up here all good");
  
  if((fIsoMethod == 0 || fIsoMethod == 1 || fIsoMethod==3) && fTPC4Iso){
    cout<<"Error: Iso_Methods with CELLS and CLUSTERS work only within EMCal "<<endl;
    cout<<"Please Set Iso_Method and TPC4Iso Accordingly!!"<<endl;
    return;
  }

  if((fIsoMethod == 0 || fIsoMethod == 1 || fIsoMethod==3) && fUEMethod> 1){
    cout<<"Error: UE_Methods with CELLS and CLUSTERS work only within EMCal"<<endl;
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
  else if(fIsoMethod==3)
    sIsoMethod = "ClustOnly";
  
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
    sBoundaries = "EMCal Acceptance";
  
  if(fWho>2 || fWho==-1){
    cout<<"Error!!! OutputMode Can Only Be 0: TTree; 1: THnSparse; 2: TH**"<<endl;
    return;
  }
  else{
    fOutput = new AliEmcalList(); // RH: Leak? fOutput already exists in base class
    fOutput->SetOwner();

    // Initialize the common output histograms
    switch(fWho)
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
	{
	  // Initialization of THnSparse
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
	  // bincells=20;
        
	  Int_t bins[] = {binPT, binM02, binETiso, binETUE, binetacl, binphicl};
        
	  fNDimensions = sizeof(bins)/sizeof(Int_t);
	  const Int_t ndims = fNDimensions;
        
	  sTitle = Form("Direct Photons: p_{T} , M02 , E_{T} Iso%s in %s, E_{T} UE %s in %s, #eta_{clus} distr,#phi_{clus} distr; p_{T} (GeV/c); M02; E_{T}^{iso%s} (GeV/c) ; E_{T}^{UE%s} (GeV/c); #eta_{cl}; #phi_{cl}", sIsoMethod.Data(), sBoundaries.Data(), sUEMethod.Data(), sBoundaries.Data(), sIsoMethod.Data(), sUEMethod.Data());
        
	  fOutputTHnS = new THnSparseF("fHnOutput",sTitle.Data(), ndims, bins);
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
          
	    // Double_t xminbis[] = {0., -10., -10., -1000., -1.0,  1.,    0};
	    // Double_t xmaxbis[] = {70., 100., 100.,  1000.,  1.0, 3.5, 1500};
          
	    // fOutMCTruth = new THnSparseF ("fOutMCTruth","E_{#gamma}, E_{T}^{iso cone}, E_{T}^{UE}, MomPDG, Eta, Phi, Label; E_{T}^{#gamma} (GeV/c); p_{T}^{Iso}(GeV/c);E_{T} ^{UE} (GeV/c); PDG; #eta; #phi; Label",7,binsMC,xminbis,xmaxbis);
	    // fOutMCTruth->Sumw2();
	    // fOutput->Add(fOutMCTruth);
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
          
	    // fOutClustMC = new THnSparseF ("fOutClustMC", "E_{T}^{clust}, M02, PDG, MOM PDG, E_{T}^{true}, #Deltax, #Deltaz, E_{T}^{iso},Label;E_{T}^{reco} (GeV/c); M02;PDG Code; Mothers' PDG Code; E_{T}^{MCtrue} (GeV/c); #Delta#phi; #Delta#eta; E_{T}^{iso} (Gev/c);Label",9,binsSMC,xminbismix,xmaxbismix);
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
	  // Initialization TH*D/TH*F
	  fPtaftM02C = new TH1D("hPtaftM02C_NC","p_{T} distribution for Clusters after shower shape cut",200,0.,100.);
	  fPtaftM02C->Sumw2();
	  fOutput->Add(fPtaftM02C);
        
	  fM02 = new TH2D("hM02_NC","M02 distribution for Neutral Clusters vs E",100,0.,100.,500,0.,5.);
	  fM02->Sumw2();
	  fOutput->Add(fM02);
        
	  fEtIsoCells = new TH1D("hEtIsoCell_NC","E_{T}^{iso cone} in iso cone distribution for Neutral Clusters with EMCal Cells",200,-0.25,99.75);
	  fEtIsoCells->SetXTitle("#Sigma E_{T}^{iso cone} (GeV/c)");
	  fEtIsoCells->Sumw2();
	  fOutput->Add(fEtIsoCells);
        
	  fEtIsoClust = new TH2D("hEtIsoClus_NC","#Sigma p_{T}^{iso cone} in iso cone distribution for Neutral Clusters with EMCal Clusters",200,0.,100.,200,0.,100.);
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
  
  // Initialize the common QA histograms
  if(fQA){
    // Include QA plots to the OutputList // DEFINE BETTER THE BINNING AND THE AXES LIMITS
    fTrackMult = new TH1D ("hTrackMult","Tracks multiplicity Distribution",100,0.,100.);
    fTrackMult->Sumw2();
    fOutput->Add(fTrackMult);
    
    fClusTime = new TH1D("hClusTime_NC","Time distribution for Clusters",800,-50.,50.);
    fClusTime->Sumw2();
    fOutput->Add(fClusTime);
    
    fEtaPhiClus = new TH2D ("hEtaPhiClusActivity","",250,-0.8,0.8, 250, 1.2, 3.4);
    // fEtaPhiClus->Sumw2();
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
    
    fNLM = new TH2D("hNLM_NC","NLM distribution for Clusters",10,0.,10.,100,0.,100.);
    fNLM->Sumw2();
    fOutput->Add(fNLM);
    
    fNLM2_NC_Acc = new TH2D("hNLM2_NC_Acc","NLM distribution for *Neutral* Clusters in acceptance",10,0.,10.,100,0.,100.);
    fNLM2_NC_Acc->Sumw2();
    fOutput->Add(fNLM2_NC_Acc);
    
    fTestIndex= new TH2D("hTestIndex","Test index for cluster",100,0.,100.,100,0.,100.);
    fTestIndex->SetXTitle("index");
    fTestIndex->SetYTitle("local index");
    fTestIndex->Sumw2();
    fOutput->Add(fTestIndex);
    
    fTestIndexE= new TH2D("hTestIndexE","Test index vs energy for cluster",200,0.,100.,100,0.,100.);
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
    
    fEtVSM02VSPisotrack = new TH3F ("hEtVSM02VSPisotrack","Energy clust vs cluster M02 vs Charged iso",70,0.,700.,200,0.,2.,400,0.,100.);
    fEtVSM02VSPisotrack->SetXTitle("Cluster E (GeV/c)");
    fEtVSM02VSPisotrack->SetYTitle("cluster M02");
    fEtVSM02VSPisotrack->SetZTitle("Isolation Charged (GeV/c)");
    fEtVSM02VSPisotrack->Sumw2();
    fOutput->Add(fEtVSM02VSPisotrack);
    
    fEtVSM02VSEisoclust = new TH3F ("hEtVSM02VSEisoclust","Energy clust vs cluster M02 vs Neutral iso",70,0.,700.,200,0.,2.,400,0.,100.);
    fEtVSM02VSEisoclust->SetXTitle("Cluster E (GeV/c)");
    fEtVSM02VSEisoclust->SetYTitle("cluster M02");
    fEtVSM02VSEisoclust->SetZTitle("Isolation Neutrals (GeV/c)");
    fEtVSM02VSEisoclust->Sumw2();
    fOutput->Add(fEtVSM02VSEisoclust);

    // fInvMassM02iso = new TH3D("hInvMassM02iso","Invariant mass vs M02 vs E_{T}^{iso cluster}",100,0.,1.,500,0.,5.,200,0.,100.);
    // fInvMassM02iso->Sumw2();
    // fOutput->Add(fInvMassM02iso);
    
    // fInvMassM02noiso = new TH3D("hInvMassM02noiso","Invariant mass vs M02 vs E_{T}^{no iso cluster}",100,0.,1.,500,0.,5.,200,0.,100.);
    // fInvMassM02noiso->Sumw2();
    // fOutput->Add(fInvMassM02noiso); 
  }

  // Initialization of all the common THistos for the 3 different outputs
  fVz = new TH1D("hVz_NC","Vertex Z distribution",100,-50.,50.);
  fVz->Sumw2();
  fOutput->Add(fVz);
  
  fEvents = new TH1D("hEvents_NC","Events",100,0.,100.);
  fEvents->Sumw2();
  fOutput->Add(fEvents);
  
  fClusEvsClusT = new TH2D("fClustTimeVSClustEn", "Distribution of cluster Time as a function of the cluster Energy", 70, 0., 70., 120, -40., 80.);
  fClusEvsClusT->SetXTitle("E_{T,clus} (GeV/c)    ");
  fClusEvsClusT->SetYTitle("Time_{clus} (ns)    ");
  fClusEvsClusT->Sumw2();
  fOutput->Add(fClusEvsClusT);
  
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
  
  // fTracksConeEtaPt = new TH3D("hTracksConeEtaPt","#Sigma vs #eta vs E_{T}",200,0.,100.,320,-0.8,0.8,200,0.,100.);
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
    // CREATE THE TH2 specific for the MC Analysis
    // Maybe to be added in the THNSparse, or cloning two or three axes and add the specific MC Truth info
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
  // return;
}

  //________________________________________________________________________
Double_t* AliAnalysisTaskEMCALPhotonIsolation::GenerateFixedBinArray(Int_t n, Double_t min, Double_t max) const
{
    // Generate the bin array for the ThnSparse
  
  Double_t *bins = new Double_t[n+1];
  
  Double_t binWidth = (max-min)/n;
  bins[0] = min;
  for(Int_t i = 1; i <= n; i++){
    bins[i] = bins[i-1]+binWidth;
  }
  
  return bins;
}


  //________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::ExecOnce()
{
    // Init the analysis

  // Tracks for CT Matching
  AliTrackContainer *tracks = GetTrackContainer("tpconlyMatch");
  // Printf("name of the first track container: %s", tracks->GetClassName().Data());
  if(!tracks){
    AliError(Form("%s: This task needs a 1particle container!", GetName()));
    return;
  }

  // Tracks for Isolation
  AliTrackContainer *tracksANA = GetTrackContainer("filterTracksAna");
  // Printf("name of the second track container: %s", tracksANA->GetClassName().Data());
  
  if(!tracksANA){
    AliError(Form("%s: This task needs a 2particle container!", GetName()));
    return;
  }

  // Clusters
  AliClusterContainer *clusters = GetClusterContainer(0);
  if(!clusters){
    AliError(Form("%s: This task needs a cluster container!", GetName()));
    return;
  }
  
  // Init the EMCAL Framework
  AliAnalysisTaskEmcal::ExecOnce();
  if(!fLocalInitialized){
    
    AliError(Form("AliAnalysisTask not initialized"));
    return;
  }
}

  //______________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::SelectCandidate(AliVCluster *coi)
{
  Int_t index=0;
  TLorentzVector vecCOI;
  coi->GetMomentum(vecCOI,fVertex);

  Double_t coiTOF = coi->GetTOF()*1e9;
  index=coi->GetID();
  if(coi->GetM02()>=0.1)
    fClusEvsClusT->Fill(vecCOI.Pt(),coiTOF);
  if(!fIsMC){
    if(coiTOF< -30. || coiTOF > 30.)
      return kFALSE;
  }

  fPtaftTime->Fill(vecCOI.Pt());
  
  if((coi->GetNCells() < 2))
    return kFALSE;

  fPtaftCell->Fill(vecCOI.Pt());
  
  Int_t nlm=0;
  AliVCaloCells * fCaloCells =InputEvent()->GetEMCALCells();
  if(fCaloCells){
    nlm = GetNLM(coi,fCaloCells);
    AliDebug(1,Form("NLM = %d",nlm));
    
    if(coi->E()>=5. && coi->E()<70. && fQA)
      fNLM->Fill(nlm,coi->E());

    if(fIsNLMCut && fNLMCut>0 && fNLMmin>0){ // If the NLM cut is enabled, this is a loop to reject clusters with more than the defined NLM (should be 1 or 2 (merged photon decay clusters))
      if(nlm > fNLMCut || nlm < fNLMmin ){
	// AliWarning(Form("NLM = %d --- NLM min = %d --- NLMcut = %d",nlm,fNLMmin,fNLMCut));
	return kFALSE;
      }
    }
  }
  else{
    AliDebug(1,Form("Can't retrieve EMCal cells"));
    return kFALSE;
  }
  
  fPtaftNLM->Fill(vecCOI.Pt());

  if(fTMClusterRejected){
    if(ClustTrackMatching(coi,kTRUE))
      return kFALSE;
  }

  fPtaftTM->Fill(vecCOI.Pt());

  if((coi->GetDistanceToBadChannel() < 2))
    return kFALSE;
  
  fPtaftDTBC->Fill(vecCOI.Pt());
  
  if(!CheckBoundaries(vecCOI))
    return kFALSE;
  
  fPtaftFC->Fill(vecCOI.Pt());
  
  if(fQA)
    fTestIndexE->Fill(vecCOI.Pt(),index);

  if(vecCOI.Pt()<5.)
    return kFALSE;
  
  if(fQA)
    fNLM2_NC_Acc->Fill(nlm,coi->E());
  
  return kTRUE;
}

  //______________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::Run()
{
    // Run the analysis.

  // fTest1+=1;
  // vertex cuts

  AliTrackContainer *tracks = GetTrackContainer("tpconlyMatch");
  if(!tracks){
    Printf("Cannot find the tracks for CT Matching");
    return kFALSE;
  }
  
  if(tracks->GetTrackFilterType() != AliEmcalTrackSelection::kTPCOnlyTracks){
    AliWarning(Form("CT matching NOT performed with TPCOnly Tracks"));
    AliWarning(Form("You better be sure of what you are doing"));
  }
  
  AliTrackContainer *tracksANA = GetTrackContainer("filterTracksAna");
  if(!tracksANA){
    Printf("Cannot find the tracks for Isolation");
    return kFALSE;
  }
  
  // Printf("FilterType of the tracks for Analysis: %d \t(should be %d)", tracksANA->GetTrackFilterType(),AliEmcalTrackSelection::kHybridTracks);
  
  //  AliError(Form("\n\n\n\nGO CHECK the Settings!!!! Is Isolation calculated with filteredTracks?\n\n\n\n"));
  if(tracksANA->GetTrackFilterType() != AliEmcalTrackSelection::kHybridTracks){
    AliWarning(Form("Isolation NOT calculated with HybridTracks"));
    AliWarning(Form("You better be sure of what you are doing"));
  }
  
  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  
  if(fVertex[2]>10. || fVertex[2]<-10.)
    return kFALSE;
  
  AliClusterContainer* clusters = GetClusterContainer(0);
  Int_t nbTracksEvent;
  nbTracksEvent = InputEvent()->GetNumberOfTracks();
  
  if(fRejectionEventWithoutTracks && (nbTracksEvent == 0))
    return kFALSE;
  
  fEvents->Fill(0); // Fill event number histogram
  
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
    if(!MCSimTrigger(fVevent,fTriggerLevel1) && fAnalysispPb)
      return kFALSE;
  }
  
  fVz->Fill(fVertex[2]); // Fill Vertex Z histogram
  
  // fOutClusters->Delete(); // Delete output USEFUL LATER FOR CONTAINER CREATION !!

  Int_t index=0;
  
  if(fIsMC){
    AliAODMCHeader *mcHeader; // Is this object useful compared to fmcHeader?

    fAODMCParticles = static_cast <TClonesArray*>(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
    fmcHeader = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

    if(!fIsMC)
      return kFALSE;
    // AliInfo(Form("It's a MC analysis %e",fAODMCParticles));
    if(!fStack && !fAODMCParticles){
      AliError("No MC stack saved\n");
      return kFALSE;
    }    
    // cout << "There's a list of particles" << endl;

    // DO THIS ALSO FOR ESDs
    if(fAODMCParticles->GetEntries() < 1){
      AliError("Number of MC particles insufficient");
      return kFALSE;
    }
    
    if(fMCtruth || (fmcHeader->GetEventType()==14 || fmcHeader->GetEventType()==29)){
      AnalyzeMC();
    }
  }
  
  if(fisLCAnalysis){
    // Get the leading particle

    // AliEmcalParticle *emccluster = static_cast<AliEmcalParticle*>(clusters->GetLeadingParticle());    
    AliVCluster *coi = (clusters->GetLeadingCluster());
    
    if(!coi){
      AliError(Form("No leading cluster"));
      return kFALSE;
    }
    if(!coi->IsEMCAL())
      return kFALSE;
    
    index = coi->GetID();
    TLorentzVector vecCOI;
    coi->GetMomentum(vecCOI,fVertex);

    if(fQA)  FillQAHistograms(coi,vecCOI);
    
    
    Bool_t isSelected = SelectCandidate(coi);

    if(isSelected){
      for(auto it : tracksANA->accepted()){
        AliVTrack *tr = static_cast<AliVTrack*>(it);
        if(!tr){
          AliError("No track found");
          return kFALSE;
        }
        fPtTracksVSpTNC->Fill(vecCOI.Pt(),tr->Pt());
        fPhiTracksVSclustPt->Fill(vecCOI.Pt(),tr->Phi());
        fEtaTracksVSclustPt->Fill(vecCOI.Pt(),tr->Eta());
      }

      FillGeneralHistograms(coi,vecCOI, index);
    }
  }
  else{
    // Get the entries of the Cluster Container
    // Whatever is a RETURN in LCAnalysis here is a CONTINUE,
    // since there are more than 1 Cluster per Event
    
    for(auto it : clusters->accepted()){
      AliVCluster *coi = static_cast<AliVCluster*>(it);
      if(!coi){
        AliError("No cluster found");
        return kFALSE;
      }
      if(!coi->IsEMCAL())
	return kFALSE;
      
      index=coi->GetID();
      TLorentzVector vecCOI;
      coi->GetMomentum(vecCOI,fVertex);

      if(fQA)  FillQAHistograms(coi,vecCOI);
      
      Bool_t isSelected = SelectCandidate(coi);
      
      if(isSelected){
        for(auto it : tracksANA->accepted()){
          AliVTrack *tr = static_cast<AliVTrack*>(it);
          if(!tr){
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
  }
  return kTRUE;
}

  //_________________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::FillQAHistograms(AliVCluster *coi, TLorentzVector vecCOI){
  
  switch(fWho)
    {
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
    // if(!ClustTrackMatching(coi)){
    // fPtaftTM->Fill(vecCOI.Pt());
    
    if(CheckBoundaries(vecCOI)){
      // fPtaftFC->Fill(vecCOI.Pt());
      Double_t checkM02=coi->GetM02();
      if(fM02mincut < checkM02 && checkM02 < fM02maxcut && fWho==2)
        fPtaftM02C->Fill(vecCOI.Pt());
    }
    // }
  }
}

  //___________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::MCSimTrigger(AliVEvent *eventIn, Int_t triggerLevel){
  
  if(triggerLevel<1 || triggerLevel>2){
    AliError(Form("No trigger level have been chosen for the MC analysis"));
    return kFALSE;
  }
  Double_t threshold=0;
  Double_t spread=0;
  
  // Int_t runNumber = InputEvent()->GetRunNumber();
  Int_t runNumber = eventIn->GetRunNumber();
  // AliError(Form("The run number is %d",runNumber));
  
  if(!runNumber)
    return kFALSE;
    
  if(runNumber < 195180 || runNumber > 197469) // LHC13a to LHC13f (to be replaced by a condition on fAnalysispPb?)
    return kFALSE;
  
  // AliError(Form("The run is after %d",runNumber));
  // TString fired = InputEvent()->GetFiredTriggerClasses();
  // AliError(Form("Trigger used in the events %s",fired));
  
  if(triggerLevel==1){
    threshold = 11.5;
    spread = 0.5;
  }
  
  if(triggerLevel==2){
    threshold = 7.2;
    spread = 0.3;
  }
  
  if(spread != 0.){
    TF1* triggerSmearing =  new TF1("triggerSmearing","[0]*exp(-0.5*((x-[1])/[2])**2)",0,15);
    triggerSmearing->SetParameter(0, 1/(spread*TMath::Sqrt(TMath::Pi()*2)));
    triggerSmearing->SetParameter(1, threshold);
    triggerSmearing->SetParameter(2, spread);
    threshold = triggerSmearing->GetRandom();
    delete triggerSmearing;
  }
  
  // AliError(Form("Pass the trigger function definition"));
  AliClusterContainer *clusters = GetClusterContainer(0);
  Int_t localIndex=0;
  
  for(auto it : clusters->accepted()){
    AliVCluster* coi = static_cast<AliVCluster*>(it);
    // AliError(Form("Retrieve clusters"));
    
    if(!coi)
      continue;
    //  AliError(Form("Retrieve the cluster"));

    if(coi->E() > threshold){
      // AliError(Form("A cluster passes the energy criterion"));
      return kTRUE;
    }
  }
  return kFALSE;
}

  //__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::ClustTrackMatching(AliVCluster *clust, Bool_t candidate){

    // Check if the cluster matches to a track
  
  AliTrackContainer* tracks = GetTrackContainer(0);
  AliVTrack* mt = 0;
  TLorentzVector vecClust;
  clust->GetMomentum(vecClust,fVertex);
  
  Int_t nbMObj = clust -> GetNTracksMatched();
  if(tracks->GetTrackFilterType() != AliEmcalTrackSelection::kTPCOnlyTracks)  AliError(Form("NO TPC only tracks"));
  
  Double_t distCT=0.;
  
  if(nbMObj == 0)
    return kFALSE;
  
  for(Int_t i=0;i<nbMObj;i++){
    if(fIsEsd){
      Int_t imt = clust->GetTrackMatchedIndex(0);
      if(imt >= 0) mt = static_cast<AliVTrack*>(tracks->GetAcceptParticle(imt));
    }
    else{
      mt = static_cast<AliVTrack*>(clust->GetTrackMatched(i));
      UInt_t rejectionReason = 0;
      if(!tracks->AcceptParticle(mt, rejectionReason)) mt = 0;
    }
    // Int_t imt = partC->GetMatchedObjId(i);
    
    if(!mt)
      continue;
    
    Double_t deltaEta,deltaPhi;
    Double_t deta = 999;
    Double_t dphi = 999;
    
    Double_t veta = mt->GetTrackEtaOnEMCal();
    Double_t vphi = mt->GetTrackPhiOnEMCal();
    
    Float_t pos[3] = {0};
    clust->GetPosition(pos);
    TVector3 cpos(pos);
    Double_t ceta = cpos.Eta();
    Double_t cphi = cpos.Phi();

    deta = veta-ceta;
    dphi = TVector2::Phi_mpi_pi(vphi-cphi);

    if(fQA && candidate){
      fDeltaETAClusTrack->Fill(deta);
      fDeltaPHIClusTrack->Fill(dphi);
    }

    distCT = TMath::Sqrt(deta*deta+dphi*dphi);
    fCTdistVSpTNC->Fill(vecClust.Pt(),distCT);
    
    if(candidate){
      deltaEta = fdetacut;
      deltaPhi = fdphicut;
    }
    else{
      deltaEta = fdetacutIso;
      deltaPhi = fdphicutIso;
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

    // Find the number of local maxima of a cluster (adapted from AliCalorimeterUtils)
  
  const Int_t nc = coi->GetNCells();
  
  Int_t   absIdList[nc];
  Float_t maxEList[nc];
  
  Int_t nMax = GetNLM(coi, cells, absIdList, maxEList);
  
  return nMax;
}

  //_____________________________________________________________________________________________________________________________
Int_t AliAnalysisTaskEMCALPhotonIsolation::GetNLM(AliVCluster* coi, AliVCaloCells* cells, Int_t *absIdList, Float_t *maxEList){

    // Find the cluster number of local maxima (adapted from AliCalorimeterUtils)
  
  Int_t iDigitN = 0 ;
  Int_t iDigit  = 0 ;
  Int_t absId1 = -1 ;
  Int_t absId2 = -1 ;
  const Int_t nCells = coi->GetNCells();
  
  Float_t eCluster = RecalEnClust(coi, cells); // Recalculate cluster energy, avoid non lin. correction
  Float_t localMaxCutE = 0.1;
  Float_t locMaxCutEDiff = 0.0;
  
  Float_t emax  = 0;
  Int_t   idmax = -1;
  for(iDigit = 0; iDigit < nCells ; iDigit++){
    absIdList[iDigit] = coi->GetCellsAbsId()[iDigit];
    Float_t en = cells->GetCellAmplitude(absIdList[iDigit]);
    // Printf("Cell Energy Before Recalculation: %.4f",en);
    RecalAmpCell(en,absIdList[iDigit]);
    
    if(en > emax){
      emax  = en ;
      idmax = absIdList[iDigit] ;
    }
  }

  for(iDigit = 0 ; iDigit < nCells; iDigit++){
    if(absIdList[iDigit] >= 0){
      absId1 = coi->GetCellsAbsId()[iDigit];
      Float_t en1 = cells->GetCellAmplitude(absId1);
      RecalAmpCell(en1,absId1);
      
      for(iDigitN = 0; iDigitN < nCells; iDigitN++){
	absId2 = coi->GetCellsAbsId()[iDigitN];
        
	if(absId2==-1 || absId2==absId1)
	  continue;
        
	Float_t en2 = cells->GetCellAmplitude(absId2);
	RecalAmpCell(en2,absId2);
        
	if(AreNeighbours(absId1, absId2)){
	  if( en1 > en2 ){
	    absIdList[iDigitN] = -1;
	    if(en1 < en2 + locMaxCutEDiff){
	      absIdList[iDigit] = -1;
	    }
	  }
	  else{
	    absIdList[iDigit] = -1;
	    if(en1 > en2 - locMaxCutEDiff){
	      absIdList[iDigitN] = -1;
	    }
	  }
	} // if AreNeighbours
      } // for iDigitN
    } // if absIdList
  } // for iDigit
  
  iDigitN = 0;

  for(iDigit = 0; iDigit < nCells; iDigit++){
    if(absIdList[iDigit] >= 0){
      absIdList[iDigitN] = absIdList[iDigit];
      
      Float_t en = cells->GetCellAmplitude(absIdList[iDigit]);
      RecalAmpCell(en,absIdList[iDigit]);
    
      if(en < localMaxCutE)
	continue; // Maxima only with seed energy at least
      
      maxEList[iDigitN] = en;
      iDigitN++;
    }
  }
  
  if(iDigitN == 0){
    AliDebug(1,Form("No local maxima found, assign highest energy cell as maxima, id %d, en cell %2.2f, en cluster %2.2f",idmax,emax,coi->E()));
    iDigitN      = 1;
    maxEList[0]  = emax;
    absIdList[0] = idmax;
  }
  
  AliDebug(1,Form("In coi E %2.2f (wth non lin. %2.2f), M02 %2.2f, M20 %2.2f, N maxima %d",coi->E(),eCluster,coi->GetM02(),coi->GetM20(),iDigitN));

  return iDigitN ;
}

  //__________________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::AreNeighbours(Int_t absId1, Int_t absId2 ) const {

    // Check if two cells are neighbours (adapted from AliCalorimeterUtils)
  
  Bool_t areNeighbours = kFALSE;
  
  Int_t iSupMod1 = -1, iTower1 = -1, iIphi1 = -1, iIeta1 = -1, iphi1 = -1, ieta1 = -1;
  Int_t iSupMod2 = -1, iTower2 = -1, iIphi2 = -1, iIeta2 = -1, iphi2 = -1, ieta2 = -1;
  
  Int_t phidiff =  0, etadiff =  0;
  
  // First cell
  fGeom->GetCellIndex(absId1,iSupMod1,iTower1,iIphi1,iIeta1);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod1,iTower1,iIphi1, iIeta1,iphi1,ieta1);
  
  // Second cell
  fGeom->GetCellIndex(absId2,iSupMod2,iTower2,iIphi2,iIeta2);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod2,iTower2,iIphi2, iIeta2,iphi2,ieta2);
  
  
  if(iSupMod1 != iSupMod2){
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2-1
    // C Side impair SM, nSupMod%2=1; A side pair SM nSupMod%2=0
    if(iSupMod1%2) ieta1+=AliEMCALGeoParams::fgkEMCALCols;
    else            ieta2+=AliEMCALGeoParams::fgkEMCALCols;
  }
  
  phidiff = TMath::Abs( iphi1 - iphi2 );
  etadiff = TMath::Abs( ieta1 - ieta2 );
  
  if((etadiff + phidiff == 1 ))
    areNeighbours = kTRUE;
  
  return areNeighbours;
}
  //_____________________________________________________________________________________________
  /// Recalculate cell energy if recalibration factor.
  //_____________________________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::RecalAmpCell(Float_t & amp, Int_t id) const {

  Int_t iSupMod = -1, iTower = -1, iIphi = -1, iIeta = -1, iphi = -1, ieta = -1;
  Float_t ampold=amp;
  fGeom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
  
  amp *= fEMCALRecoUtils->GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
}

  //__________________________________________________________________________
  /// Recalibrate the cluster energy, considering the recalibration map and the energy of the cells that compose the cluster.
  //__________________________________________________________________________
Float_t AliAnalysisTaskEMCALPhotonIsolation::RecalEnClust(AliVCluster * coi, AliVCaloCells * cells){

  // Printf("Inside clust Recal");

  // Initialize some used variables
  Float_t frac = 0., energy = 0.;
  
  if(cells){
    // Get the cluster number of cells and list of absId, check what kind of cluster we have
    UShort_t * index    = coi->GetCellsAbsId();
    Double_t * fraction = coi->GetCellsAmplitudeFraction();
    Int_t ncells        = coi->GetNCells();
    
    // Loop on the cells, get the cell amplitude and recalibration factor, multiply to get the new energy
    for(Int_t icell = 0; icell < ncells; icell++){
      Int_t absId = index[icell];
      frac = fraction[icell];
      if(frac < 1e-3)
	frac = 1; // in case of EMCal, this is set as 0, not used.
      
      Float_t amp = cells->GetCellAmplitude(absId);
      RecalAmpCell(amp, absId);
      
      // Printf("Recalibrate cell: EMCal, cell fraction %f, cell energy: before cal %f; after cal %f",frac,cells->GetCellAmplitude(absId),amp);
      
      energy += amp*frac;
    }
    
    AliDebug(1,Form("Energy before %f, after %f",coi->E(),energy));
    
  } // Cells available
  else{
    AliFatal("Cells pointer does not exist!");
  }
  // Printf("recalculated energy: %.4f",energy);

  return energy;
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoCellPhiBand(TLorentzVector c, Double_t &etIso, Double_t &phiBandcells){

    // Underlying events study with EMCal cells in phi band
  
  AliEMCALGeometry* emcalGeom = AliEMCALGeometry::GetInstance();
  Double_t sumEnergyPhiBandCells=0., sumEnergyConeCells=0.;
  
  // Check the cell corresponding to the leading cluster
  Int_t absId = 999;
  Bool_t cellLeadingClustID = emcalGeom->GetAbsCellIdFromEtaPhi(c.Eta(),c.Phi(),absId);

  if(!cellLeadingClustID)
    return;  
  else{
    Int_t iTower = -1;
    Int_t iModule = -1;
    Int_t imEta=-1, imPhi=-1;
    Int_t iEta =-1, iPhi =-1;
    
    emcalGeom->GetCellIndex(absId,iModule,iTower,imPhi,imEta);                    // To get the module, the tower, eta and phi for the cell corresponding to the leading cluster
    emcalGeom->GetCellPhiEtaIndexInSModule(iModule,iTower,imPhi,imEta,iPhi,iEta); // To get the cell eta and phi in the super module for the cell coming from the leading cluster
    
    // Get the row and the column of the cell corresponding to the leading cluster in EMCal
    Int_t colCellLeadingClust = iEta;
    if(iModule % 2)
      colCellLeadingClust = AliEMCALGeoParams::fgkEMCALCols + iEta ;                   // If the SM number is even you need to translate to have the one corresponding in EMCal
    Int_t rowCellLeadingClust = iPhi + AliEMCALGeoParams::fgkEMCALRows*int(iModule/2); // To have the corresponding row in EMCal
    
    // Total number or rows and columns in EMCal
    Int_t nTotalRows = AliEMCALGeoParams::fgkEMCALRows*16/3 ; // 5 + 2/3 supermodules in a row
    Int_t nTotalCols = 2*AliEMCALGeoParams::fgkEMCALCols;     // 2 supermodules in a column
    Int_t nbConeSize = int(fIsoConeRadius/0.0143);            // 0.0143 = tower acceptance
    
    // Get the cells
    AliVCaloCells * cells =InputEvent()->GetEMCALCells();
    
    // Define the max and min row and column corresponding to the isolation cone around the seed cell from the leading cluster
    Int_t iRowMinCone = rowCellLeadingClust-nbConeSize;
    if(iRowMinCone<0)
      iRowMinCone=0;
    
    Int_t iRowMaxCone = rowCellLeadingClust+nbConeSize;
    if(iRowMaxCone>AliEMCALGeoParams::fgkEMCALRows)
      iRowMaxCone=AliEMCALGeoParams::fgkEMCALRows;  // AliEMCALGeoParams::fgkEMCALRows = 24 in a supermodule
    
    Int_t iColMinCone = colCellLeadingClust - nbConeSize;
    if(iColMinCone<0)
      iColMinCone = 0;
    
    Int_t iColMaxCone = colCellLeadingClust+nbConeSize;
    if(iColMaxCone>AliEMCALGeoParams::fgkEMCALCols)
      iColMaxCone=AliEMCALGeoParams::fgkEMCALCols;  // AliEMCALGeoParams::fgkEMCALCols = 48 in a supermodule
    
    // Loop on all cells
    for(Int_t iCol=0; iCol<nTotalCols; iCol++){
      for(Int_t iRow=0; iRow<nTotalRows; iRow++){
	// Now recover the cell indexes in a supermodule
        Int_t iSector = int(iRow/AliEMCALGeoParams::fgkEMCALRows); // Check in which SM is the cell
        if(iSector==5)
	  continue;

        Int_t inModule = -1;
        Int_t iColLoc  = -1;
        if(iCol < AliEMCALGeoParams::fgkEMCALCols){ // If the SM number is odd the column is the one corresponding in the supermodule
          inModule = 2*iSector + 1;
          iColLoc  = iCol;
        }
        else if(iCol > AliEMCALGeoParams::fgkEMCALCols - 1){ // If the SM number is even the column isn't the one corresponding in the supermodule
          inModule = 2*iSector;
          iColLoc  = iCol-AliEMCALGeoParams::fgkEMCALCols;
        }
        
        Int_t iRowLoc  = iRow - AliEMCALGeoParams::fgkEMCALRows*iSector ;
        if(TMath::Abs(iCol-colCellLeadingClust)<nbConeSize && TMath::Abs(iCol+colCellLeadingClust)>nbConeSize){
          if(iRow<iRowMaxCone && iRow>iRowMinCone){
            Int_t iabsId = emcalGeom->GetAbsCellIdFromCellIndexes(iModule,iRow,iCol);
            sumEnergyPhiBandCells+=cells->GetAmplitude(iabsId); // Should be Et
          }
        }
        else if(TMath::Abs(iCol-colCellLeadingClust)>nbConeSize && TMath::Abs(iCol+colCellLeadingClust)<nbConeSize){
          Int_t iabsId = emcalGeom->GetAbsCellIdFromCellIndexes(inModule,iRowLoc,iColLoc);
          sumEnergyConeCells+=cells->GetAmplitude(iabsId); // Should be Et
        }
      }
    }
  }
  etIso = sumEnergyConeCells;
  phiBandcells = sumEnergyPhiBandCells;
}


  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoCellEtaBand(TLorentzVector c, Double_t &etIso, Double_t &etaBandcells){

    // Underlying events study with EMCal cell in eta band
  
  AliEMCALGeometry* emcalGeom = AliEMCALGeometry::GetInstance();
  Double_t sumEnergyEtaBandCells=0., sumEnergyConeCells=0.;
  
  // Check the cell corresponding to the leading cluster
  Int_t absId = 999;
  Bool_t cellLeadingClustID = emcalGeom->GetAbsCellIdFromEtaPhi(c.Eta(),c.Phi(),absId);
  if(!cellLeadingClustID)
    return;
  else{
    Int_t iTower = -1;
    Int_t iModule = -1;
    Int_t imEta=-1, imPhi=-1;
    Int_t iEta =-1, iPhi =-1;
    
    emcalGeom->GetCellIndex(absId,iModule,iTower,imPhi,imEta);                    // To get the module, the tower, eta and phi for the cell corresponding to the leading cluster
    emcalGeom->GetCellPhiEtaIndexInSModule(iModule,iTower,imPhi,imEta,iPhi,iEta); // To get the cell eta and phi in the super module for the cell coming from the leading cluster
    
    // Get the row and the column of the cell corresponding to the leading cluster in EMCal
    Int_t colCellLeadingClust = iEta;
    if(iModule % 2)
      colCellLeadingClust = AliEMCALGeoParams::fgkEMCALCols + iEta ;    // If the SM number is even you need to translate to have the one corresponding in EMCal
    Int_t rowCellLeadingClust = iPhi + AliEMCALGeoParams::fgkEMCALRows*int(iModule/2); // To have the corresponding row in EMCal
    
    // total number or rows and columns in EMCal
    Int_t nTotalRows = AliEMCALGeoParams::fgkEMCALRows*16/3 ;  // 5 + 2/3 supermodules in a row
    Int_t nTotalCols = 2*AliEMCALGeoParams::fgkEMCALCols;      // 2 supermodules in a column
    Int_t nbConeSize = int(fIsoConeRadius/0.0143);             // 0.0143 = tower acceptance
    
    // Get the cells
    AliVCaloCells * cells =InputEvent()->GetEMCALCells();
    
    // Define the max and min row and column corresponding to the isolation cone around the seed cell from the leading cluster
    Int_t iRowMinCone = rowCellLeadingClust-nbConeSize;
    if(iRowMinCone<0)
      iRowMinCone=0;
    
    Int_t iRowMaxCone = rowCellLeadingClust+nbConeSize;
    if(iRowMaxCone>AliEMCALGeoParams::fgkEMCALRows)
      iRowMaxCone=AliEMCALGeoParams::fgkEMCALRows;  // AliEMCALGeoParams::fgkEMCALRows = 24 in a supermodule
    
    Int_t iColMinCone = colCellLeadingClust-nbConeSize;
    if(iColMinCone<0)
      iColMinCone = 0;
    
    Int_t iColMaxCone = colCellLeadingClust+nbConeSize;
    if(iColMaxCone>AliEMCALGeoParams::fgkEMCALCols)
      iColMaxCone=AliEMCALGeoParams::fgkEMCALCols;  // AliEMCALGeoParams::fgkEMCALCols = 48 in a supermodule
    
    // loop on all cells
    for(Int_t iCol=0; iCol<nTotalCols; iCol++){
      for(Int_t iRow=0; iRow<nTotalRows; iRow++){
	// now recover the cell indexes in a supermodule
        Int_t iSector = int(iRow/AliEMCALGeoParams::fgkEMCALRows); // Check in which SM is the cell
        if(iSector==5)
	  continue;

        Int_t inModule = -1;
        Int_t iColLoc  = -1;
        if(iCol < AliEMCALGeoParams::fgkEMCALCols){ // If the SM number is odd the column is the one corresponding in the supermodule
          inModule = 2*iSector + 1;
          iColLoc  = iCol;
        }
        else if(iCol > AliEMCALGeoParams::fgkEMCALCols - 1){ // If the SM number is even the column isn't the one corresponding in the supermodule
          inModule = 2*iSector;
          iColLoc  = iCol-AliEMCALGeoParams::fgkEMCALCols;
        }
        
        Int_t iRowLoc  = iRow - AliEMCALGeoParams::fgkEMCALRows*iSector ;
        if(TMath::Abs(iCol-colCellLeadingClust)<nbConeSize && TMath::Abs(iCol+colCellLeadingClust)>nbConeSize){
          if(iCol<iColMaxCone && iCol>iColMinCone){
            Int_t iabsId = emcalGeom->GetAbsCellIdFromCellIndexes(iModule,iRow,iCol);
            sumEnergyEtaBandCells+=cells->GetAmplitude(iabsId); // Should be Et
          }
        }
        else if(TMath::Abs(iCol-colCellLeadingClust)>nbConeSize && TMath::Abs(iCol+colCellLeadingClust)<nbConeSize){
          Int_t iabsId = emcalGeom->GetAbsCellIdFromCellIndexes(inModule,iRowLoc,iColLoc);
          sumEnergyConeCells+=cells->GetAmplitude(iabsId); // Should be Et
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
  Double_t minPhi = 0., maxPhi = 0., minEta = 0., maxEta = 0.;

  if(fPeriod != ""){
    minEta = fGeom->GetArm1EtaMin()+0.03;
    maxEta = fGeom->GetArm1EtaMax()-0.03;
    minPhi = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

    if(fPeriod.Contains("12") || fPeriod.Contains("13"))
      maxPhi = (fGeom->GetArm1PhiMax()-20.)*TMath::DegToRad()-0.03; // fGeom->GetArm1PhiMax()-20. = 180. deg (in order not to take the two disabled SM into account in 2012-2013)
    else
      maxPhi = (fGeom->GetArm1PhiMax())*TMath::DegToRad()-0.03;
  }
  else{
    minPhi = (4./9.)*TMath::Pi()+0.03;
    maxPhi = TMath::Pi()-0.03;
    minEta = -0.67;
    maxEta = 0.67;
  }

  // Needs a check on the same cluster
  // AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliClusterContainer *clusters = GetClusterContainer(0);
  Int_t localIndex=0;
  TLorentzVector nClust; // STILL NOT INITIALIZED
  AliVCluster *candidate;
  
  for(auto it : clusters->accepted()){ // Check the position of other clusters with respect to the trigger cluster
    
    AliVCluster* coi = static_cast<AliVCluster*>(it);
    localIndex=coi->GetID();
    
    if(localIndex==index){
      candidate = static_cast<AliVCluster*>(it);
      continue;
    }
    
    phiClust = etaClust = clustTOF = 0.;
    coi->GetMomentum(nClust,fVertex);
    phiClust = nClust.Phi();
    etaClust = nClust.Eta();
    
    if(fExtraIsoCuts){
      if((coi->GetNCells() < 2))
	continue;
      if((coi->GetDistanceToBadChannel() < 2))
	continue;
    }
    
    clustTOF = coi->GetTOF()*1e9;
    
    if(!fIsMC){
      if(clustTOF<-30. || clustTOF>30.)
	continue;
    }
    
    if(fTMClusterInConeRejected){
      if(ClustTrackMatching(coi,kFALSE))
	continue;
    }
    
    if(nClust.E()<0.3)
      continue;

    // Redefine c.Eta()/c.Phi() from the cluster we passed to the function
    radius = TMath::Sqrt(TMath::Power(phiClust-c.Phi(),2)+TMath::Power(etaClust-c.Eta(),2)); // Define the radius between the leading cluster and the considered cluster
    
    if(radius>fIsoConeRadius){ // The cluster is outside the isolation cone -> add the cluster pT to pT_UE
      if(TMath::Abs(phiClust - c.Phi()) < fIsoConeRadius)
        sumEnergyPhiBandClus += nClust.Pt();
    }
    else if(radius<fIsoConeRadius && radius != 0.){ // The cluster is inside the isolation cone -> add the cluster pT to pT_iso
      // Printf("Cluster Not Matched, Inside the Cone, with Energy %.4lf",eTcluster);
      sumEnergyConeClus += nClust.Pt();
      if(fQA){
        fTestEtaPhiCone->Fill(c.Eta(),c.Phi());
        fTestIndex->Fill(index,localIndex);
        fTestLocalIndexE->Fill(nClust.Pt(),localIndex);
      }
    }
  }

  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  // fTracksAna = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("FilterTracksAna"));
  // Name hard coded to use the defined tracks for analysis
  
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  // const Int_t nbTracks = fTracksAna->GetEntries();
  // Int_t iTracks = 0;
  
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
    // if(!(eTrack->IsHybridGlobalConstrainedGlobal())){Printf("skipping track %d because it's not an hybrid\n",eTrack->GetID()); continue;}
    
    if(!fIsEsd){
      aodEtrack = static_cast<AliAODTrack*>(eTrack);
      if(!(aodEtrack->IsHybridGlobalConstrainedGlobal()))
	continue;
    }
    
    if((eTrack->Pt())<0.2)
      continue;
    
    if((eTrack->GetStatus() & AliAODTrack::kITSrefit) != AliAODTrack::kITSrefit)
      continue;
    
    if(!fIsEsd){
      Double_t frac = 0;
      Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
      Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls> 0)  frac =  nclsS / ncls ;
      
      if(frac > 0.4)
	continue;
    }
    
    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();
    if(etaTrack < minEta || etaTrack > maxEta || phiTrack < minPhi || phiTrack > maxPhi) // Skip tracks outside EMCal
      continue;
    
    radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2));

    if(radius<fIsoConeRadius ){ // The track is inside the isolation cone -> add the track pT to pT_iso
      sumpTConeCharged+=eTrack->Pt();
      iTracksCone++;
    }
    else{
      if(TMath::Abs(phiTrack - c.Phi()) < fIsoConeRadius){ // The track is outside the isolation cone -> add the track pT to pT_UE
        sumpTPhiBandTracks += eTrack->Pt();
      }
    }
  } // End of tracks loop
  
  fTestEnergyCone->Fill(c.Pt(),sumEnergyConeClus,sumpTConeCharged);

  if(fIsoMethod==1 && fQA){
    fEtVSM02VSPisotrack->Fill(c.Pt(),candidate->GetM02(),sumpTConeCharged);
    fEtVSM02VSEisoclust->Fill(c.Pt(),candidate->GetM02(),sumEnergyConeClus);
  }

  if(fWho==2){
    fTrackMultvsSumChargedvsUE->Fill(iTracksCone,sumpTConeCharged, sumpTPhiBandTracks);
    fTrackMultvsPt->Fill(iTracksCone,c.Pt());
  }

  if(fIsoMethod==1)
    ptIso = sumEnergyConeClus + sumpTConeCharged;
  else if(fIsoMethod==3)
    ptIso = sumEnergyConeClus;
  
  if(fWho==2){
    fPhiBandUEClust->Fill(c.Pt() , sumEnergyPhiBandClus);
    fPhiBandUETracks->Fill(c.Pt() , sumpTPhiBandTracks);
  }

  if(fIsoMethod==1)
    phiBandclus = sumEnergyPhiBandClus + sumpTPhiBandTracks;
  else if(fIsoMethod==3)
    phiBandclus = sumEnergyPhiBandClus;  
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoClusEtaBand(TLorentzVector c, Double_t &ptIso, Double_t &etaBandclus, Int_t index){

    // Underlying events study with clusters in eta band
  
  Float_t sumEnergyEtaBandClus =0., sumEnergyConeClus=0., sumpTConeCharged=0, sumpTEtaBandTracks=0.;
  Double_t clustTOF=0;
  Double_t minPhi = 0., maxPhi = 0., minEta = 0., maxEta = 0.;

  if(fPeriod != ""){
    minEta = fGeom->GetArm1EtaMin()+0.03;
    maxEta = fGeom->GetArm1EtaMax()-0.03;
    minPhi = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

    if(fPeriod.Contains("12") || fPeriod.Contains("13"))
      maxPhi = (fGeom->GetArm1PhiMax()-20.)*TMath::DegToRad()-0.03; // fGeom->GetArm1PhiMax()-20. = 180. deg (in order not to take the two disabled SM into account in 2012-2013)
    else
      maxPhi = (fGeom->GetArm1PhiMax())*TMath::DegToRad()-0.03;
  }
  else{
    minPhi = (4./9.)*TMath::Pi()+0.03;
    maxPhi = TMath::Pi()-0.03;
    minEta = -0.67;
    maxEta = 0.67;
  }
  
  // AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliClusterContainer *clusters = GetClusterContainer(0);
  Int_t localIndex=0;
  AliVCluster *candidate;
  TLorentzVector nClust; // STILL NOT INITIALIZED
  
  for(auto it : clusters->accepted()){ // Check the position of other clusters with respect to the trigger cluster
    
    AliVCluster* coi = static_cast<AliVCluster*>(it);
    localIndex=coi->GetID();

    if(localIndex==index){
      candidate = static_cast<AliVCluster*>(it);
      continue;
    }
    
    coi->GetMomentum(nClust,fVertex);
    
    Double_t phiClust =nClust.Phi();
    Double_t etaClust= nClust.Eta();
    Double_t eTcluster=0, radius;
    
    if(fExtraIsoCuts){
      if((coi->GetNCells() < 2))
	continue;
      if((coi->GetDistanceToBadChannel() < 2))
	continue;
    }
    //  printf("\nCluster ID %d with pT: %.4f\t\t Eta: %.4f \t Phi: %.4f \t time %.3f", coi->GetID(), nClust.Pt(),etaClust,phiClust,coi->GetTOF()*1e9);

    clustTOF = coi->GetTOF()*1e9;

    if(!fIsMC){
      if(clustTOF<-30. || clustTOF>30.)
	continue;
    }
    
    if(fTMClusterInConeRejected){
      if(ClustTrackMatching(coi,kFALSE))
	continue;
    }
    
    if(nClust.E()<0.3)
      continue;
    
    // Redefine c.Eta()/c.Phi() from the cluster we passed to the function
    radius = TMath::Sqrt(TMath::Power(phiClust-c.Phi(),2)+TMath::Power(etaClust-c.Eta(),2)); // Define the radius between the leading cluster and the considered cluster

    if(radius>fIsoConeRadius){ // The cluster is outside the isolation cone -> add the cluster pT to pT_UE
      if(TMath::Abs(etaClust - c.Eta()) < fIsoConeRadius)
        sumEnergyEtaBandClus += nClust.Pt();
    }
    else if(radius<fIsoConeRadius && radius != 0.){ // The cluster is inside the isolation cone -> add the cluster pT to pT_iso
      sumEnergyConeClus += nClust.Pt();
      if(fQA){
        fTestEtaPhiCone->Fill(c.Eta(),c.Phi());
        fTestIndex->Fill(index,localIndex);
        fTestLocalIndexE->Fill(nClust.Pt(),localIndex);
      }
    }
  } // End of clusters loop
  
  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  // Printf("ntracks accepted: %d",tracksAna->GetNAcceptedTracks());

  if(tracksAna->GetTrackFilterType() != AliEmcalTrackSelection::kHybridTracks)
    AliError(Form("NOT Hybrid Tracks"));
  // Printf("Name of the tracks used for Isolation: %s",tracksAna->GetName());

  tracksAna->ResetCurrentID();

  AliVTrack *eTrack = 0x0;
  Int_t iTracksCone = 0;
  AliAODTrack *aodEtrack=0x0;

  while((eTrack = static_cast<AliVTrack*>(tracksAna->GetNextAcceptParticle()))){
    if(!eTrack){
      AliError(Form("No tracks in collection"));
      continue;
    }
    Double_t phiTrack, etaTrack,radius;
    
    if(!fIsEsd){
      aodEtrack = static_cast<AliAODTrack*>(eTrack);
      if(!(aodEtrack->IsHybridGlobalConstrainedGlobal())){
	// Printf("skipping track %d because it's not an hybrid\n",eTrack->GetID());
        continue;
      }
    }
    
    if((eTrack->Pt())<0.2)
      continue;
    
    if((eTrack->GetStatus() & AliAODTrack::kITSrefit) != AliAODTrack::kITSrefit)
      continue;
    
    if(!fIsEsd){
      Double_t frac = 0;
      Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
      Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls> 0) frac =  nclsS / ncls ;
      
      if(frac > 0.4)
	continue;
    }

    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();
    
    if(etaTrack < minEta || etaTrack > maxEta || phiTrack < minPhi || phiTrack > maxPhi) // Skip tracks outside EMCal
      continue;

    radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2));

    if(radius<fIsoConeRadius){ // The track is inside the isolation cone -> add the track pT to pT_iso
      sumpTConeCharged += eTrack->Pt();
      iTracksCone++;
    }
    else{ // The track is outside the isolation cone -> add the track pT to pT_UE
      if(TMath::Abs(etaTrack - c.Eta()) < fIsoConeRadius){
	// Printf("but phi similar, so UE EtaBand!!");
        sumpTEtaBandTracks += eTrack->Pt();
      }
    }
  } // End of tracks loop
  
  // Printf("\ntotal activity in isolation Cone from Tracks %.4lf",sumpTConeCharged);
  
  fTestEnergyCone->Fill(c.Pt(),sumEnergyConeClus,sumpTConeCharged);
  
  if(fIsoMethod==1 && fQA){
    fEtVSM02VSPisotrack->Fill(c.Pt(),candidate->GetM02(),sumpTConeCharged);
    fEtVSM02VSEisoclust->Fill(c.Pt(),candidate->GetM02(),sumEnergyConeClus);
  }

  if(fWho==2){
    fTrackMultvsSumChargedvsUE->Fill(iTracksCone,sumpTConeCharged, sumpTEtaBandTracks);
    fTrackMultvsPt->Fill(iTracksCone,c.Pt());
  }
  
  if(fIsoMethod==1)
    ptIso = sumEnergyConeClus + sumpTConeCharged;
  else if(fIsoMethod==3)
    ptIso = sumEnergyConeClus;

  if(fWho==2){
    fEtaBandUEClust->Fill(c.Pt() , sumEnergyEtaBandClus);
    fEtaBandUETracks->Fill(c.Pt() , sumpTEtaBandTracks);
  }

  if(fIsoMethod==1)
    etaBandclus = sumEnergyEtaBandClus + sumpTEtaBandTracks;
  else if(fIsoMethod==3)
    etaBandclus = sumEnergyEtaBandClus;
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackPhiBand(TLorentzVector c, Double_t &ptIso, Double_t &phiBandtrack){

    // Underlying events study with tracks in phi band
  
  Double_t sumpTConeCharged=0., sumpTPhiBandTrack=0.;
  Double_t minPhi = 0., maxPhi = 2.*TMath::Pi(), minEta = -0.87, maxEta = 0.87;
  
  if(!fTPC4Iso){
    if(fPeriod != ""){
      minEta = fGeom->GetArm1EtaMin()+0.03;
      maxEta = fGeom->GetArm1EtaMax()-0.03;
      minPhi = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

      if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	maxPhi = (fGeom->GetArm1PhiMax()-20.)*TMath::DegToRad()-0.03; // fGeom->GetArm1PhiMax()-20. = 180. deg (in order not to take the two disabled SM into account in 2012-2013)
      else
	maxPhi = (fGeom->GetArm1PhiMax())*TMath::DegToRad()-0.03;
    }
    else{
      minEta = -0.67;
      maxEta = 0.67;
      minPhi = (4./9.)*TMath::Pi()+0.03;
      maxPhi = TMath::Pi()-0.03;
    }
  }
  
  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  // Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());

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
      if(!(aodEtrack->IsHybridGlobalConstrainedGlobal()))
	continue;
    }
    
    if((eTrack->Pt())<0.2)
      continue;
    
    if((eTrack->GetStatus() & AliAODTrack::kITSrefit) != AliAODTrack::kITSrefit)
      continue;
    
    if(!fIsEsd){
      Double_t frac = 0;
      Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
      Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls> 0)  frac =  nclsS / ncls ;
      
      if(frac > 0.4) continue;
    }
    
    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();

    if((phiTrack < maxPhi) && (phiTrack > minPhi) && (etaTrack < maxEta) && (etaTrack > minEta)){
      radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2));
      if(radius>fIsoConeRadius){ // The track is outside the isolation cone -> add the track pT to pT_UE
        if(TMath::Abs(etaTrack - c.Eta()) < fIsoConeRadius)
          sumpTPhiBandTrack += eTrack->Pt();
      }
      else{ // The track is inside the isolation cone -> add the track pT to pT_iso
        sumpTConeCharged += eTrack->Pt();
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

    // Underlying events study with tracks in eta band
  
  Double_t sumpTConeCharged=0., sumpTEtaBandTrack=0.;
  Double_t minPhi = 0., maxPhi = 2.*TMath::Pi(), minEta = -0.87, maxEta = 0.87;
  
  if(!fTPC4Iso){
    if(fPeriod != ""){
      minEta = fGeom->GetArm1EtaMin()+0.03;
      maxEta = fGeom->GetArm1EtaMax()-0.03;
      minPhi = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

      if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	maxPhi = (fGeom->GetArm1PhiMax()-20.)*TMath::DegToRad()-0.03; // fGeom->GetArm1PhiMax()-20. = 180. deg (in order not to take the two disabled SM into account in 2012-2013)
      else
	maxPhi = (fGeom->GetArm1PhiMax())*TMath::DegToRad()-0.03;
    }
    else{
      minEta = -0.67;
      maxEta = 0.67;
      minPhi = (4./9.)*TMath::Pi()+0.03;
      maxPhi = TMath::Pi()-0.03;
    }
  }

  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");  
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  // Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());

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
      if(!(aodEtrack->IsHybridGlobalConstrainedGlobal()))
	continue;
    }
    
    if((eTrack->Pt())<0.2)
      continue;
    
    if((eTrack->GetStatus() & AliAODTrack::kITSrefit) != AliAODTrack::kITSrefit)
      continue;
    
    if(!fIsEsd){
      Double_t frac = 0;
      Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
      Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls> 0)  frac =  nclsS / ncls ;
      
      if(frac > 0.4)
	continue;
    }
    
    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();

    if( (phiTrack < maxPhi) && (phiTrack > minPhi) && (etaTrack < maxEta) && (etaTrack > minEta)){
      radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2));
      if(radius>fIsoConeRadius){ // The track is outside the isolation cone -> add the track pT to pT_UE
        if(TMath::Abs(phiTrack - c.Phi()) < fIsoConeRadius)
          sumpTEtaBandTrack += eTrack->Pt();
      }
      else{ // The track is inside the isolation cone -> add the track pT to pT_iso
        sumpTConeCharged += eTrack->Pt();
	iTracksCone++;
      }
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
  
  if(phiCone1 < 0.)
    phiCone1 += 2*TMath::Pi();
  
  AliTrackContainer *tracks = GetTrackContainer(0);
  if(!tracks){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  // Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());

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
      if(!(aodEtrack->IsHybridGlobalConstrainedGlobal()))
	continue;
    }
    
    if((eTrack->Pt())<0.2)
      continue;
    
    if((eTrack->GetStatus() & AliAODTrack::kITSrefit) != AliAODTrack::kITSrefit)
      continue;
    
    if(!fIsEsd){
      Double_t frac = 0;
      Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
      Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls> 0)  frac =  nclsS / ncls ;
      
      if(frac > 0.4) continue;
    }
    
    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();
    dist2Clust = TMath::Sqrt(TMath::Power(etaTrack-etaClus, 2)+TMath::Power(phiTrack-phiClus, 2));
    
    if(dist2Clust<fIsoConeRadius){ // The track is inside the isolation cone -> add the track pT to pT_iso
      sumpTConeCharged += eTrack->Pt();
      iTracksCone++;
    }
    
    else{
      // Distances from the centres of the two Orthogonal Cones
      dist2Cone1 = TMath::Sqrt(TMath::Power(etaTrack-etaClus, 2)+TMath::Power(phiTrack-phiCone1, 2));
      dist2Cone2 = TMath::Sqrt(TMath::Power(etaTrack-etaClus, 2)+TMath::Power(phiTrack-phiCone2, 2));
      
      // The track is inside one of the two orthogonal cones -> add the track pT to pT_UE
      if((dist2Cone1 < fIsoConeRadius) || (dist2Cone2 < fIsoConeRadius))
	sumpTPerpConeTrack += eTrack->Pt();
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

    // Underlying events study with tracks in full TPC except a back to back band
  
  Double_t sumpTConeCharged=0., sumpTTPCexceptB2B=0.;
  
  AliTrackContainer *tracks = GetTrackContainer(0);  
  if(!tracks){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }
  // Printf("Name of the tracks used for Isolation: %s",(tracks->GetClassName()).Data());

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
      if(!(aodEtrack->IsHybridGlobalConstrainedGlobal()))
	continue;
    }
    
    if((eTrack->Pt())<0.2)
      continue;
    
    if((eTrack->GetStatus() & AliAODTrack::kITSrefit) != AliAODTrack::kITSrefit)
      continue;
    
    if(!fIsEsd){
      Double_t frac = 0;
      Float_t ncls  = Float_t(aodEtrack->GetTPCncls ());
      Float_t nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls> 0)  frac =  nclsS / ncls ;
      
      if(frac > 0.4) continue;
    }

    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();

    radius = TMath::Sqrt(TMath::Power(phiTrack-c.Phi(),2)+TMath::Power(etaTrack-c.Eta(),2));
    
    if(radius>fIsoConeRadius){ // The track is outside the isolation cone -> add the track pT to pT_UE
      dphiUp = c.Phi() + TMath::Pi() - fIsoConeRadius;
      dphiDown = c.Phi() + TMath::Pi() + fIsoConeRadius;

      if(phiTrack < dphiDown && phiTrack> dphiUp)
	sumpTTPCexceptB2B += eTrack->Pt();
    }
    else{ // The track is inside the isolation cone -> add the track pT to pT_iso
      sumpTConeCharged += eTrack->Pt();
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

    // Check if the cone around the considered cluster is in EMCal acceptance
  // AliInfo("Inside CheckBoundaries\n");
  
  Double_t minPhiBound = 0. , minEtaBound = 0., maxPhiBound = 0., maxEtaBound = 0.;
  Bool_t isINBoundaries;
  
  if(fTPC4Iso){
    minEtaBound = -0.87+fIsoConeRadius;
    maxEtaBound = 0.87-fIsoConeRadius;

    if(fPeriod != ""){
      minPhiBound = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

      if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	maxPhiBound = (fGeom->GetArm1PhiMax()-20.)*TMath::DegToRad()-0.03;
      else
	maxPhiBound = (fGeom->GetArm1PhiMax())*TMath::DegToRad()-0.03;
    }
    else{
      minPhiBound = (4./9.)*TMath::Pi()+0.03;
      maxPhiBound = TMath::Pi()-0.03;
    }
  }
  else{
    if(fPeriod != ""){
      minEtaBound = fGeom->GetArm1EtaMin()+0.03+fIsoConeRadius;
      maxEtaBound = fGeom->GetArm1EtaMax()-0.03-fIsoConeRadius;
      minPhiBound = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03+fIsoConeRadius;

      if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	maxPhiBound = (fGeom->GetArm1PhiMax()-20.)*TMath::DegToRad()-0.03-fIsoConeRadius;
      else
	maxPhiBound = (fGeom->GetArm1PhiMax())*TMath::DegToRad()-0.03-fIsoConeRadius;
    }
    else{
      minEtaBound = -0.67+fIsoConeRadius;
      maxEtaBound = 0.67-fIsoConeRadius;
      minPhiBound = (4./9.)*TMath::Pi()+0.03+fIsoConeRadius;
      maxPhiBound = TMath::Pi()-0.03-fIsoConeRadius;
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
  
  // cout<<"\n\n\n\n\n\n\nInside Look4Particle \n For Cluster "<<clusterlabel<<"\t\t"<<energyCLS<<"\t\t"<<etaCLS<<"\t\t"<<phiCLS<<"\t\t"<<ss<<"\t\t"<<isolation<<"\n\n\n\n"<<endl;
  
  if(!fIsMC){
    AliWarning("Not a Monte-Carlo run!!");
    return;
  }
  // AliInfo(Form("It's a MC analysis %e",fAODMCParticles));
  
  if(!fStack && !fAODMCParticles){
    AliWarning("No particle stack!!");
    return;
  }

  if(fAODMCParticles->GetEntries() < 1){
    AliWarning("Number of tracks insufficient!!");
    return;
  }
  
  Int_t ndimsMCmix = fMCQAdim;
  Double_t outputvalueMCmix[ndimsMCmix];
  // cout<<"dimensions of the array: "<<ndimsMCmix<<endl;
  
  Int_t npart=fAODMCParticles->GetEntries();
  // cout<<"Number of particles in the event: "<<npart<<endl;
  
  AliAODMCParticle *particle2Check, *momP2Check;
  Int_t clustPDG, p2clabel, p2ccharge;
  Double_t enTrue,phiTrue, etaTrue;
  Double_t dPhi,dEta ;
  Int_t clusterFromPromptPhoton=-1;
  particle2Check = static_cast<AliAODMCParticle*>(fAODMCParticles->At(clusterlabel));
  clustPDG=particle2Check->GetPdgCode();

  int mom2checkidx = particle2Check->GetMother();
  momP2Check = static_cast<AliAODMCParticle*>(fAODMCParticles->At(mom2checkidx));
  
  phiTrue = particle2Check->Phi();
  etaTrue = particle2Check->Eta();
  enTrue  = particle2Check->E()*TMath::Sin(particle2Check->Theta()); // Check if we need corrections to the energy

  dPhi = phiCLS-phiTrue;
  dEta = etaCLS-etaTrue;
  p2ccharge = particle2Check->Charge();
  
  if(clustPDG==22 || (TMath::Abs(clustPDG) == 11 && momP2Check->GetPdgCode() == 22)){
    // Direct Photon or e+/e- coming from Photon Conversion
    // Checking if the Photon is a decay product of pi0 or eta meson // Maybe include omega?

    // printf("Cluster Label: %d Asso. with a MCpar w/ PDG %d",clusterlabel,clustPDG);

    if(clustPDG==22){
      // printf("\twhose mother is a %d",momP2Check->GetPdgCode());
      if(momP2Check->GetPdgCode() != 22){
	// printf("  which is not a photon so the cluster is from a decay ");
        if(momP2Check->GetPdgCode()==111 || momP2Check->GetPdgCode()==221){
	  // printf(" of a pi0 or a eta mesons");
          clusterFromPromptPhoton=5;
          
          Int_t idxdaug1 = momP2Check->GetFirstDaughter();
          Int_t idxdaug2 = momP2Check->GetLastDaughter();

          if(idxdaug1 == clusterlabel){ // Cluster associated with the 1st daughter? Then look if also the 2nd daughter contributes to the cluster energy
            if(idxdaug2<npart){         // 2nd daughter within List of Particles
              AliAODMCParticle *daug2 = static_cast<AliAODMCParticle*>(fAODMCParticles->At(idxdaug2));

              if(daug2->GetPdgCode() == 22 && ( daug2->Phi() - phiTrue )< 0.01 && ( daug2->Eta() - etaTrue )< 0.01){ // Very tight cut because they are photons
                enTrue += daug2->E()*TMath::Sin(daug2->Theta());
                clusterFromPromptPhoton=7; // Contribution from both daughters
              }
              else
                clusterFromPromptPhoton=6; // Contribution from one daughter
            }
          }
          else{ // Then the cluster MUST BE associated with the 2nd daughter!! Look if also the 1st daughter contributes to the cluster energy
            AliAODMCParticle *daug1 = static_cast<AliAODMCParticle*>(fAODMCParticles->At(idxdaug1));

            if(daug1->GetPdgCode() == 22 && ( daug1->Phi()-phiTrue )< 0.01 && ( daug1->Eta()-etaTrue )< 0.01){ // Very tight cut because they are photons
              enTrue += daug1->E()*TMath::Sin(daug1->Theta());
              clusterFromPromptPhoton=6; // Contribution from both daughters
            }
            else
              clusterFromPromptPhoton=7; // Contribution from one daughter
          }
          fpi0VSclusterVSIsolation->Fill(momP2Check->E()*TMath::Sin(momP2Check->Theta()), energyCLS, isolation);
          fpi0VSclusterVSM02->Fill(momP2Check->E()*TMath::Sin(momP2Check->Theta()), energyCLS,ss);
          fpi0VSM02VSIsolation->Fill(momP2Check->E()*TMath::Sin(momP2Check->Theta()), ss, isolation);
        }
        else{
	  // printf("  of a non considered meson/baryon");
          clusterFromPromptPhoton=8; // Undefined
        }
      }
      else{
        clusterFromPromptPhoton=1; // True prompt photon
	// printf("  so we have a prompt photon\n");
      }
    }
    else{ // Cluster created by e+/e- from Photon Conversion
      // printf(" whose mother PDG is %d and occupies the stack position at: %d\n",momP2Check->GetPdgCode(),mom2checkidx);
      Int_t g_momindex = momP2Check->GetMother();
      AliAODMCParticle *gMomP2Check=static_cast<AliAODMCParticle*>(fAODMCParticles->At(g_momindex));

      if( mom2checkidx == 8 || (g_momindex == 8 && gMomP2Check->GetPdgCode()==22 && momP2Check->GetPdgCode()==22)){
        clusterFromPromptPhoton=1; // e+/e- from converted prompt photon
	// printf(" This means it is a e+/e- cluster from a Converted PromptPhoton.\n");
      }
      else{
	// printf(" This means it is a e+/e- cluster from a Converted DECAYPhoton.\n");
        clusterFromPromptPhoton=5; // Cluster created by a photon but not a prompt one
      }

      Int_t firstidx=momP2Check->GetFirstDaughter();
      Int_t lastidx=momP2Check->GetLastDaughter();

      if(clusterFromPromptPhoton==1){
        if( firstidx == clusterlabel ){ // Cluster associated with the 1st electron? Then look if also the 2nd electron contributes to the cluster energy
          if( lastidx < npart ){        // 2nd daughter within List of Particles
            AliAODMCParticle *last=static_cast<AliAODMCParticle*>(fAODMCParticles->At(lastidx));

            if(( last->Phi() - phiTrue ) < fdphicut && ( last->Eta() - etaTrue ) < fdetacut ){ // Same proximity cut as the CPV
              enTrue += last->E()*TMath::Sin(last->Theta());
              clusterFromPromptPhoton=3; // Contribution from both daughters
	      // printf(" The cluster HAS actually contribution from both daughters but is asso to the second daughter (lower energy)\n");
            }
            else
	      clusterFromPromptPhoton=2; // Contribution from one daughter
          }
        }
        else{ // Cluster associated to the 2nd daughter!! Look if also the 1st daughter contributes to the cluster energy
          
          AliAODMCParticle *first=static_cast<AliAODMCParticle*>(fAODMCParticles->At(firstidx));
          if(( first->Phi() - phiTrue ) < fdphicut && ( first->Eta() - etaTrue ) < fdetacut ){ // Same proximity cut as the CPV
            enTrue += first->E()*TMath::Sin(first->Theta());
            clusterFromPromptPhoton=3; // Contribution from both daughters
	    // printf("cluster HAS actually contribution from both daughters but is asso to the first daughter (higher energy)\n");
          }
          else
            clusterFromPromptPhoton=2;// Contribution from one daughter
        }
      }

      if(clusterFromPromptPhoton >= 5){ // Check on wheter also the 2nd gamma from pi0/eta decay contributes to the energy of the cluster
	// This further check is implemented to take care of very asymmetric decays
        Int_t idxgrandma = momP2Check->GetMother();
        AliAODMCParticle *grandma=static_cast<AliAODMCParticle*>(fAODMCParticles->At(idxgrandma));

        if( grandma->GetPdgCode() == 111 || momP2Check->GetPdgCode() == 221 ){ // Add also omega mesons, Lambda barion, neutral Kaons?
          Int_t idxaunt1 = grandma->GetFirstDaughter();
          Int_t idxaunt2 = grandma->GetLastDaughter();

          if( idxaunt1 ==  mom2checkidx ){ // The 1st daughter of the pi0/eta is the mother of the electron that produced the cluster
	    // Check if the 2nd pi0/eta daughter contributes to the cluster energy
            AliAODMCParticle *aunt=static_cast<AliAODMCParticle*>(fAODMCParticles->At(idxaunt2));

            if(( aunt->Phi() - phiTrue ) < 0.01 && ( aunt->Eta() - etaTrue ) < 0.01 ){
              enTrue += aunt->E()*TMath::Sin(aunt->Theta());
              clusterFromPromptPhoton=7; // Contribution from both daughters
            }
            else
              clusterFromPromptPhoton=6; // Contribution from one daughter
          }
          else{ // The 2nd daughter of the pi0/eta is the mother of the electron that produced the cluster
	    // Check if the 1st pi0/eta daughter contributes to the cluster energy
            AliAODMCParticle *aunt=static_cast<AliAODMCParticle*>(fAODMCParticles->At(idxaunt1));

            if(( aunt->Phi() - phiTrue ) < 0.01 && ( aunt->Eta() - etaTrue ) < 0.01 ){
              enTrue += aunt->E()*TMath::Sin(aunt->Theta());
              clusterFromPromptPhoton=6; // Contribution from both daughters
            }
            else
              clusterFromPromptPhoton=7; // Contribution from one daughter
          }
          fpi0VSclusterVSIsolation->Fill(grandma->E()*TMath::Sin(grandma->Theta()), energyCLS, isolation);
          fpi0VSclusterVSM02->Fill(grandma->E()*TMath::Sin(momP2Check->Theta()), energyCLS, ss);
          fpi0VSM02VSIsolation->Fill(grandma->E()*TMath::Sin(grandma->Theta()), ss, isolation);
        }
        else
          clusterFromPromptPhoton=8; // Undefined
      }
    }
    // printf("\nCluster %d  PDG: %d  (Mom is a %d) with pT: %f  ",clusterlabel,clustPDG,momP2Check->GetPdgCode(),energyCLS);
    // printf(" with clusterFromPromptPhoton stored: %d for cluster w/label %d\n",clusterFromPromptPhoton,clusterlabel);
  }
  else{
    if(p2ccharge==0){
      clusterFromPromptPhoton=10;
    }
    else{
      clusterFromPromptPhoton=11;
    }
    // AliWarning(Form("Hadronic cluster with energy: %f, M02: %f, cluster PDG: %d, mother PDG: %d, truth energy: %f, dphi: %f, deta %f, isolation energy: %f, cluster label: %d, charge: %d",energyCLS, ss, clustPDG, momP2Check->GetPdgCode(), enTrue, dPhi, dEta, isolation, clusterFromPromptPhoton, p2ccharge ));
  }
  // printf("filling fOutClustMC\n");
  
  outputvalueMCmix[0] = energyCLS;
  outputvalueMCmix[1] = ss;
  outputvalueMCmix[2] = clustPDG;
  outputvalueMCmix[3] = momP2Check->GetPdgCode();
  outputvalueMCmix[4] = enTrue;
  outputvalueMCmix[5] = dPhi;
  outputvalueMCmix[6] = dEta;
  outputvalueMCmix[7] = isolation;
  outputvalueMCmix[8] = clusterFromPromptPhoton;
  // clusterFromPP=1 ->clusterlabel = 8 TruePromptPhoton;
  // clusterFromPP=2 ->clusterlabel = indexe+/e- with 1 contribution to the Energy;
  // clusterFromPP=3 ->clusterlabel = indexe+/e- with 2 contributions to the Energy;
  // clusterFromPP=6 -> clusterlabel= indexgamma1/2 (or e1e2e3e4) with contribution from max 2 electrons to the Energy;
  // clusterFromPP=7 -> clusterlabel= indexgamma1/2 (or e1e2e3e4) with 4 contribution to the energy;
  // clusterFromPP=8 -> clusterlabel= Gamma decay NOT from pi0/eta decay.
  // clusterFromPP=10 -> clusterlabel= Neutral hadronic particle.
  // clusterFromPP=11 -> clusterlabel= Charged hadronic particle.
  if(fWho==1)
    fOutClustMC->Fill(outputvalueMCmix);
  
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
  
  // AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliClusterContainer *clusters = GetClusterContainer(0);
  Int_t localIndex=0;

  for(auto it :clusters->accepted()){ // Check the position of other clusters with respect to the leading cluster
    AliVCluster* coi = static_cast<AliVCluster*>(it);
    
    localIndex=coi->GetID();
    if(localIndex==index)
      continue;
    else{
      localIndex++;
      
      TLorentzVector nClust; // STILL NOT INITIALIZED
      coi->GetMomentum(nClust,fVertex);
      
      // Float_t phiClust =nClust.Phi();
      // Float_t etaClust= nClust.Eta();
      // Float_t eTcluster=0;
      
      if((coi->GetNCells() < 2))
	continue;

      if((coi->GetDistanceToBadChannel() < 2))
	continue;
      
      clustTOF = coi->GetTOF()*1e9;
      if(!fIsMC){
        if(clustTOF<-30. || clustTOF>30.)
	  continue;
      }
      
      // The TM cut is commented out in order not to reject the cluster partner if it is near to a track
      // if(ClustTrackMatching(coi))
      // 	continue;
      
      invMass = (c+nClust).M();
    }
  } // End of clusters loop
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::IsolationAndUEinEMCAL(AliVCluster *coi, Double_t& isolation,Double_t& ue,Double_t eTThreshold, Int_t index){

  // Printf("Inside IsolationAncUEinEMCal");

  Double_t isoConeArea = TMath::Pi()*fIsoConeRadius*fIsoConeRadius;
  Double_t etaBandArea = ((fGeom->GetArm1EtaMax()-0.03)-(fGeom->GetArm1EtaMin()+0.03))*2.*fIsoConeRadius-isoConeArea;
  Double_t phiBandArea = 0.;

  if(fPeriod != ""){
    if(fPeriod.Contains("12") || fPeriod.Contains("13"))
      phiBandArea = (((fGeom->GetArm1PhiMax()-20.)*TMath::DegToRad()-0.03)-((fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03))*2.*fIsoConeRadius-isoConeArea;
    else
      phiBandArea = (((fGeom->GetArm1PhiMax())*TMath::DegToRad()-0.03)-((fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03))*2.*fIsoConeRadius-isoConeArea;
  }
  else{
    phiBandArea = (5./9.)*TMath::Pi()*2.*fIsoConeRadius-isoConeArea;
  }

  TLorentzVector vecCOI;
  coi->GetMomentum(vecCOI,fVertex);
  Double_t m02COI=coi->GetM02();
  Double_t eTCOI=vecCOI.Et();
  
  // EMCal Only for Acceptance of Cells/Clusters/Tracks
  switch(fIsoMethod)
    {
    case 0: // EMCal CELLS
      switch(fUEMethod)
	{
	case 0: // Phi band
	  EtIsoCellPhiBand(vecCOI, isolation, ue);

	  ue = ue * (isoConeArea / phiBandArea); // Normalisation of UE wrt UE area
        
	  if(fWho==2){
	    fPhiBandUECells->Fill(vecCOI.Pt() , ue);
	    fEtIsoCells->Fill(isolation);
	  }

	  if(isolation<eTThreshold){
	    if(fWho==2) fEtIsolatedCells->Fill(eTCOI);
	    fEtisolatedT=eTCOI;
	    fPtisolatedT=vecCOI.Pt();
	  }
	  break;

	case 1: // Eta band
	  EtIsoCellEtaBand(vecCOI, isolation, ue);

	  ue = ue * (isoConeArea / etaBandArea); // Normalisation of UE wrt UE area

	  if(fWho==2){
	    fEtaBandUECells->Fill(vecCOI.Pt() , ue);
	    fEtIsoCells->Fill(isolation);
	  }

	  if(isolation<eTThreshold){
	    if(fWho==2) fEtIsolatedCells->Fill(eTCOI);
	    fEtisolatedT=eTCOI;
	    fPtisolatedT=vecCOI.Pt();
	  }
	  break;
	}
      break;
      
    case 1: case 3: // EMCal CLUSTERS + TRACKS
      switch(fUEMethod)
	{
	case 0: // Phi band
	  EtIsoClusPhiBand(vecCOI, isolation, ue,index);

	  ue = ue * (isoConeArea / phiBandArea); // Normalisation of UE wrt UE area
	  break;

	case 1: // Eta band
	  EtIsoClusEtaBand(vecCOI, isolation, ue,index);

	  ue = ue * (isoConeArea / etaBandArea); // Normalisation of UE wrt UE area
	  break;
	}

      if(fWho==2){
        fPtvsM02vsSum->Fill(vecCOI.Pt(),coi->GetM02(),isolation);

        isolation=isolation-ue; // UE subtraction

        fPtvsM02vsSumUE->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        fEtIsoClust->Fill(vecCOI.Pt(),isolation);
      }

      if(isolation<eTThreshold){
	FillInvMassHistograms(kTRUE, m02COI, vecCOI, index, isolation);

	if(fWho==2){
	  fPtvsM02iso->Fill(vecCOI.Pt(),coi->GetM02());
	  fPtIsolatedNClust->Fill(vecCOI.Pt());
	}

	fPtisoT=vecCOI.Pt();
	fM02isoT=m02COI;
        
	if(fM02mincut < m02COI && m02COI < fM02maxcut){
	  if(fWho==2)
	    fEtIsolatedClust->Fill(eTCOI);

	  fEtisolatedT=eTCOI;
	  fPtisolatedT=vecCOI.Pt();
	}
      }
      else{
	if(isolation>3.)
	  FillInvMassHistograms(kFALSE, m02COI, vecCOI, index, isolation);
        
	if(fWho==2)
	  fPtvsM02noiso->Fill(vecCOI.Pt(),coi->GetM02());

	fPtnoisoT=vecCOI.Pt();
	fM02noisoT=m02COI;
      }
      break;
      
    case 2: // EMCal TRACKS
      switch(fUEMethod)
	{
	case 0: // Phi band
	  PtIsoTrackPhiBand(vecCOI, isolation, ue);
	  if(fWho==2)
	    fPhiBandUETracks->Fill(vecCOI.Pt() , ue);

	  ue = ue * (isoConeArea / phiBandArea); // Normalisation of UE wrt UE area

	case 1: // Eta band
	  PtIsoTrackEtaBand(vecCOI, isolation, ue);
	  if(fWho==2)
	    fEtaBandUETracks->Fill(vecCOI.Pt() , ue);

	  ue = ue * (isoConeArea / etaBandArea); // Normalisation of UE wrt UE area
	  break;

	// case 2: // Cones
        //   PtIsoTrackOrthCones(vecCOI, absId, isolation, ue);
        //   break;
	// case 3: // Full TPC
        //   PtIsoTrackFullTPC(vecCOI, absId, isolation, ue);
        //   break;
	}

      if(fWho==2){
        fPtvsM02vsSum->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        
        isolation=isolation-ue; // UE subtraction

        fPtvsM02vsSumUE->Fill(vecCOI.Pt(),coi->GetM02(),isolation);        
        fPtIsoTrack->Fill(vecCOI.Pt() , isolation);
      }

      if(isolation<eTThreshold){
	FillInvMassHistograms(kTRUE, m02COI, vecCOI, index, isolation);

	if(fWho==2){
	  fPtvsM02iso->Fill(vecCOI.Pt(),coi->GetM02());
	  fPtIsolatedNTracks->Fill(vecCOI.Pt());
	}

	fPtisoT=vecCOI.Pt();
	fM02isoT=m02COI;
        
	if(fM02mincut < m02COI && m02COI < fM02maxcut){
	  if(fWho==2)
	    fEtIsolatedTracks->Fill(eTCOI);
	  
	  fEtisolatedT=eTCOI;
	  fPtisolatedT=vecCOI.Pt();
	}
      }
      else{
	if(isolation>3.)
	  FillInvMassHistograms(kFALSE, m02COI, vecCOI, index, isolation);
        
	if(fWho==2)
	  fPtvsM02noiso->Fill(vecCOI.Pt(),coi->GetM02());

	fPtnoisoT=vecCOI.Pt();
	fM02noisoT=m02COI;
      }
      break;
    }
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::IsolationAndUEinTPC(AliVCluster *coi, Double_t& isolation,Double_t& ue,Double_t eTThreshold, Int_t index){
  
    // EMCal + TPC (Only tracks for the Isolation since IsoCone Goes Out of EMCal)
  
  Double_t isoConeArea   = TMath::Pi()*fIsoConeRadius*fIsoConeRadius;
  Double_t etaBandAreaTr = 1.74*2.*fIsoConeRadius-isoConeArea; // 1.74 = 2*0.87
  Double_t phiBandAreaTr = 2.*TMath::Pi()*2.*fIsoConeRadius-isoConeArea;
  Double_t perpConesArea = 2.*isoConeArea;
  Double_t fullTPCArea   = 1.74*2.*TMath::Pi()-1.74*2.*fIsoConeRadius-isoConeArea;
  
  TLorentzVector vecCOI;
  coi->GetMomentum(vecCOI,fVertex);
  Double_t m02COI=coi->GetM02();
  Double_t eTCOI=vecCOI.Et();
  
  switch(fUEMethod)
    {
    case 0: // Phi band
      PtIsoTrackPhiBand(vecCOI, isolation, ue);
      fPhiBandUETracks->Fill(vecCOI.Pt() , ue);

      ue = ue * (isoConeArea / phiBandAreaTr); // Normalisation of UE wrt UE area

      if(fWho==2)
	fPtvsM02vsSum->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
      
      if(fAnalysispPb)
	isolation=isolation-ue; // UE subtraction
      
      if(fWho==2){
        fPtvsM02vsSumUE->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        fPtIsoTrack->Fill(vecCOI.Pt(), isolation);
      }
      
      if(isolation<eTThreshold){
	FillInvMassHistograms(kTRUE, m02COI, vecCOI, index,isolation);

	if(fWho==2){
	  fPtvsM02iso->Fill(vecCOI.Pt(),coi->GetM02());
	  fPtIsolatedNTracks->Fill(vecCOI.Pt());
	}

	fPtisoT=vecCOI.Pt();
	fM02isoT=m02COI;
        
	if(fM02mincut < m02COI && m02COI < fM02maxcut){
	  if(fWho==2)
	    fEtIsolatedTracks->Fill(eTCOI);

	  fEtisolatedT=eTCOI;
	  fPtisolatedT=vecCOI.Pt();
	}
      }
      else{
	if(isolation>3.)
	  FillInvMassHistograms(kFALSE, m02COI, vecCOI, index,isolation);

	if(fWho==2)
	  fPtvsM02noiso->Fill(vecCOI.Pt(),coi->GetM02());

	fPtnoisoT=vecCOI.Pt();
	fM02noisoT=m02COI;
      }
      break;

    case 1: // Eta band
      PtIsoTrackEtaBand(vecCOI, isolation, ue);
      if(fWho==2)
	fEtaBandUETracks->Fill(vecCOI.Pt() , ue);

      ue = ue * (isoConeArea / etaBandAreaTr); // Normalisation of UE wrt UE area

      if(fWho==2)
	fPtvsM02vsSum->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
      
      if(fAnalysispPb)
	isolation=isolation-ue; // UE subtraction
      
      if(fWho==2){
        fPtvsM02vsSumUE->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        fPtIsoTrack->Fill(vecCOI.Pt(), isolation);
      }

      if(isolation<eTThreshold){
	FillInvMassHistograms(kTRUE, m02COI, vecCOI, index,isolation);

	if(fWho==2){
	  fPtvsM02iso->Fill(vecCOI.Pt(),coi->GetM02());
	  fPtIsolatedNTracks->Fill(vecCOI.Pt());
	}

	fPtisoT=vecCOI.Pt();
	fM02isoT=m02COI;
        
	if(fM02mincut < m02COI && m02COI < fM02maxcut){
	  if(fWho==2)
	    fEtIsolatedTracks->Fill(eTCOI);
	  
	  fEtisolatedT=eTCOI;
	  fPtisolatedT=vecCOI.Pt();
	}
      }
      else{
	if(isolation>3.)
	  FillInvMassHistograms(kFALSE, m02COI, vecCOI, index,isolation);

	if(fWho==2)
	  fPtvsM02noiso->Fill(vecCOI.Pt(),coi->GetM02());

	fPtnoisoT=vecCOI.Pt();
	fM02noisoT=m02COI;
      }
      break;

    case 2: // Cones
      PtIsoTrackOrthCones(vecCOI, isolation, ue);

      if(fWho==2)
	fPerpConesUETracks->Fill(vecCOI.Pt() , ue);

      ue = ue * (isoConeArea / perpConesArea); // Normalisation of UE wrt UE area

      if(fWho==2)
	fPtvsM02vsSum->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
      
      if(fAnalysispPb)
	isolation=isolation-ue; // UE subtraction
      
      if(fWho==2){
        fPtvsM02vsSumUE->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        fPtIsoTrack->Fill(vecCOI.Pt(), isolation);
      }

      if(isolation<eTThreshold){
	FillInvMassHistograms(kTRUE, m02COI, vecCOI, index, isolation);

	if(fWho==2){
	  fPtvsM02iso->Fill(vecCOI.Pt(),coi->GetM02());
	  fPtIsolatedNTracks->Fill(vecCOI.Pt());
	}

	fPtisoT=vecCOI.Pt();
	fM02isoT=m02COI;
        
	if(fM02mincut < m02COI && m02COI < fM02maxcut){
	  if(fWho==2)
	    fEtIsolatedTracks->Fill(eTCOI);

	  fEtisolatedT=eTCOI;
	  fPtisolatedT=vecCOI.Pt();
	}
      }
      else{
	if(isolation>3.)
	  FillInvMassHistograms(kFALSE, m02COI, vecCOI, index, isolation);

	if(fWho==2)
	  fPtvsM02noiso->Fill(vecCOI.Pt(),coi->GetM02());

	fPtnoisoT=vecCOI.Pt();
	fM02noisoT=m02COI;
      }
      break;

    case 3: // Full TPC
      PtIsoTrackFullTPC(vecCOI, isolation, ue);

      if(fWho==2)
	fTPCWithoutIsoConeB2BbandUE->Fill(vecCOI.Pt() , ue);

      ue = ue * (isoConeArea / fullTPCArea); // Normalisation of UE wrt UE area
      // fill histograms for isolation
      
      if(fWho==2)
	fPtvsM02vsSum->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
      
      if(fAnalysispPb)
	isolation=isolation-ue; // UE subtraction
      
      if(fWho==2){
        fPtvsM02vsSumUE->Fill(vecCOI.Pt(),coi->GetM02(),isolation);
        fPtIsoTrack->Fill(vecCOI.Pt(), isolation);
      }
      // fTracksConeEtaPt->Fill(isolation, vecCOI.Eta(), vecCOI.Pt());
      // fTracksConeEtaM02->Fill(isolation, vecCOI.Eta(), coi->GetM02());

      if(isolation<eTThreshold){
	FillInvMassHistograms(kTRUE, m02COI, vecCOI, index, isolation);

	if(fWho==2){
	  fPtvsM02iso->Fill(vecCOI.Pt(),coi->GetM02());
	  fPtIsolatedNTracks->Fill(vecCOI.Pt());
	}

	fPtisoT=vecCOI.Pt();
	fM02isoT=m02COI;
        
	if(fM02mincut < m02COI && m02COI < fM02maxcut){
	  if(fWho==2)
	    fEtIsolatedTracks->Fill(eTCOI);

	  fEtisolatedT=eTCOI;
	  fPtisolatedT=vecCOI.Pt();
	}
      }
      else{
	if(isolation>3.)
	  FillInvMassHistograms(kFALSE, m02COI, vecCOI, index, isolation);

	if(fWho==2)
	  fPtvsM02noiso->Fill(vecCOI.Pt(),coi->GetM02());

	fPtnoisoT=vecCOI.Pt();
	fM02noisoT=m02COI;
      }
      break;
    }
}

  //__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::FillGeneralHistograms(AliVCluster *coi, TLorentzVector vecCOI, Int_t index){

    // Fill the histograms for underlying event and isolation studies
  // AliError(Form("Arrive bien dans fill general histograms"));
  
  // I would like to remove this part and fill the tracks multiplicity histogram in FillQAHistograms, is that ok for thnSparses? (especially cause here the histogram is filled several times per event)
  // AliParticleContainer *tracks = static_cast<AliParticleContainer*>(fParticleCollArray.At(0));

  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  tracksAna->ResetCurrentID();
  
  const Int_t nTracks = tracksAna->GetNAcceptedTracks();
  // Printf("Ntracks for the event with this cluster: %d", nTracks);

  if(fQA)
    fTrackMult->Fill(nTracks);
  
  // Printf("After Loop on Tracks");
  Double_t eTCOI = 0., m02COI = 0.;
  
  // Definition of the Array for Davide's Output
  const Int_t ndims =   fNDimensions;
  Double_t outputValues[ndims];
  
  eTCOI = vecCOI.Et();
  
  Int_t nlm=0;
  AliVCaloCells * fCaloCells =InputEvent()->GetEMCALCells();
  if(fCaloCells)
    nlm = GetNLM(coi,fCaloCells);
  
  // Printf("cluster ID %d with nlm %d . M02 BEFORE possible smearing %.4lf . Do we set the smearing ? %s ",coi->GetID(),nlm, m02COI = coi->GetM02(), fSSsmearing? "Yes":"No");
  
  if(fSSsmearing){
    if(coi->GetM02()>0.1){
      // Printf("Smearing for only clusters with nlm = %d" ,fWhich);
      if(nlm==1){
        if((fSSsmearwidth != 0.)){
          TRandom3 *ran=new TRandom3(0);

          if(fWhich==0){ // Landau Smearing
            Float_t smear = ran->Landau(fSSsmear_mean,fSSsmearwidth);

            if(fSSsmear_mean==0 || (fSSsmear_mean !=0 && coi->GetID()%3==0))
              m02COI = coi->GetM02() + smear;
          }
          else{ // Gaussian Smearing
            Float_t smear = ran->Gaus(fSSsmear_mean,fSSsmearwidth);

            if(fSSsmear_mean==0 || (fSSsmear_mean !=0 && coi->GetID()%3==0))
              m02COI = coi->GetM02() + smear;
          }

        }
        else{
          AliWarning("The Smearing is set but the width of the distribution is null!\nNOT DOING ANYTHING for the Shower Shape!");
          m02COI = coi->GetM02();
        }
      }
      else
        m02COI = coi->GetM02();
    } // No else for "if(coi->GetM02()>0.1)"!
  }
  else{
    AliWarning("Smearing not SET!");
    m02COI = coi->GetM02();
  }

  // ******** Isolation and UE calculation with different methods *********
  
  Double_t eTThreshold = 5.;
  
  switch(fEtIsoMethod)
    {
    case 0:  // SumEt < EtThr
      eTThreshold = fEtIsoThreshold;
      break;
      
    case 1:  // SumEt < %Ephoton
      eTThreshold = fEtIsoThreshold * eTCOI;
      break;
      
    case 2:  // Etmax < EtThr
      eTThreshold = fEtIsoThreshold;
      if(eTCOI<eTThreshold && fWho==2){ // Photon candidate, cuts have to be decided after studies
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
  
  // Here we should call something to know the number of tracks...
  // Soon I'll put in this version the "old way", please let me know if
  // any of you could do the same with the JET framework

  switch(fWho)
    {
    case 0:
      flambda0T=m02COI;   // for all neutral clusters
      fEtT=vecCOI.Et();   // for all neutral clusters
      fPtT=vecCOI.Pt();   // for all neutral clusters
      fetaT=vecCOI.Eta(); // for all neutral clusters
      fphiT=vecCOI.Phi(); // for all neutral clusters
      fsumEtisoconeT=isolation;
      // AliError(Form("lambda 0 %f",flambda0T));
      fsumEtUE=ue;
      // Does not fill the tree anyway files too big
      // fOutputTree->Fill();
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
    if(TMath::Abs(etap)>=(fGeom->GetArm1EtaMax()-0.03) || (phip<=((fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03) || phip>= ((fGeom->GetArm1PhiMax())*TMath::DegToRad()-0.03)))
      return;
    else{
      switch(fUEMethod)
	{
	case 0: // Phi band
	  if(TMath::Abs(eta-etap)<fIsoConeRadius)
	    sumUE += mcpp->E()*TMath::Sin(mcpp->Theta());
	  else
	    return;

	  break;

	case 1: // Eta band
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
      switch(fUEMethod)
	{
	case 0:{ // Phi band
	  if(TMath::Abs(eta-etap)<fIsoConeRadius)
	    sumUE += mcpp->E()*TMath::Sin(mcpp->Theta());
	  else
	    return;

	  break;
	} 
	case 1:{ // Eta band
	  if(TMath::Abs(phi-phip)<fIsoConeRadius)
	    sumUE += mcpp->E()*TMath::Sin(mcpp->Theta());
	  else
	    return;
          
	  break;
	}
	case 2:{ // Orthogonal Cones
	  double etacone1= eta;
	  double etacone2= eta;
	  double phicone1= phi - TMath::PiOver2();
	  double phicone2= phi + TMath::PiOver2();
          
	  if(phicone1 < 0.)
	    phicone1 += 2*TMath::Pi();
          
	  if(TMath::Sqrt(TMath::Power(etap-etacone1,2)+TMath::Power(phip-phicone1,2))< fIsoConeRadius ||
	     TMath::Sqrt(TMath::Power(etap-etacone2,2)+TMath::Power(phip-phicone2,2))< fIsoConeRadius)
	    sumUE += mcpp->Pt();
	  else
	    return;
          
	  break;
        }
	case 3:{ // Full TPC
        
	  // Double_t phiup= phi +TMath::Pi()+fIsoConeRadius;
	  // Double_t phidown= phi +TMath::Pi()-fIsoConeRadius;
	  //
	  // if(phip < phidown || phip > phiup ) // TO BE CHECKED
	  // continue;
	  break;
        }
	}
    }
  }
}
  //_________________________________________________________________________

void AliAnalysisTaskEMCALPhotonIsolation::CalculateUEDensityMC(Double_t& sumUE){
  
  Double_t isoConeArea   = TMath::Pi()*fIsoConeRadius*fIsoConeRadius;
  Double_t etaBandArea   = ((fGeom->GetArm1EtaMax()-0.03)-(fGeom->GetArm1EtaMin()+0.03))*2.*fIsoConeRadius-isoConeArea;
  Double_t phiBandArea   = (5./9.)*TMath::Pi()*2.*fIsoConeRadius-isoConeArea;
  Double_t etaBandAreaTr = 1.74*2.*fIsoConeRadius-isoConeArea;
  Double_t phiBandAreaTr = 2.*TMath::Pi()*2.*fIsoConeRadius-isoConeArea;
  Double_t perpConesArea = 2.*isoConeArea;
  Double_t fullTPCArea   = 1.74*2.*TMath::Pi()-1.74*2.*fIsoConeRadius-isoConeArea;

  if(!fTPC4Iso){
    switch(fUEMethod)
      {
      case 0:
	sumUE = sumUE * (isoConeArea / phiBandArea);
	break;

      case 1:
	sumUE = sumUE * (isoConeArea / etaBandArea);
	break;
      }
  }
  else{
    switch(fUEMethod)
      {
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
  
  if(!fIsMC)
    return;
  // AliInfo(Form("It's a MC analysis %e",fAODMCParticles));

  if(!fStack && !fAODMCParticles){
    cout<<"No stack saved\n"; return;
  }
  
  // cout<<"there's a List of particles"<<endl;
  // DO THIS ALSO FOR ESDs
  
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
      if(TMath::Abs(multTracks->Eta())<=0.87 && multTracks->Charge() != 0)
        nFSParticles++;
      else
        continue;
    } // implement final state particle condition
    else
      continue;
  }
  // AliInfo(Form("number of particles in the array %d",nTracks));

  AliAODMCParticle *mcpart, *mom, *mcpp,*mcsearch, *mcfirst, *mcfirstmom,*matchingtrack, *mum;
  
  // Bool_t prompt=kFALSE;
  Double_t mcEnergy, maxE, energy;
  Int_t pdg, mompdg, photonlabel;
  Double_t mcFirstEta=0., mcFirstPhi=0.;
  
  // AliAODMCParticle *mcfirst = static_cast<AliAODMCParticle*>(fAODMCParticles->At(0));
  // AliAODMCParticle *mcp, *mcpmaxE, *mcpp, *mom;

  if(!fisLCAnalysis){
    for(int iTr=0;iTr<nTracks;iTr++){
      mcEnergy=0.;energy =0;
      eT=0.; phi=0.; eta=0.;
      
      mcpart = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTr));
      
      if(mcpart->GetStatus()>10)
	continue;

      if(!mcpart->IsPrimary())
	continue;

      if(!mcpart->IsPhysicalPrimary())
	continue;
      
      pdg = mcpart->GetPdgCode();
      if(pdg != 22 /*|| mcpart->GetLabel()!=8*/)
	continue;
      
      eta = mcpart->Eta();
      phi = mcpart->Phi();

      // Check photons in EMCal
      if(!fTPC4Iso){
        if((TMath::Abs(eta)>(fGeom->GetArm1EtaMax()-0.03)-fIsoConeRadius ) || (phi < ((fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03+fIsoConeRadius) || phi>((fGeom->GetArm1PhiMax())*TMath::DegToRad()-0.03-fIsoConeRadius)))
          continue;
      }
      else{
        if((TMath::Abs(eta)>0.87-fIsoConeRadius ) || (phi < ((fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03) || phi>((fGeom->GetArm1PhiMax())*TMath::DegToRad()-0.03)))
          continue;
      }
      
      // printf("\nParticle Position %d  and Label: %d  PDG: %d  Pt: %f  Eta: %f  Phi: %f",iTr, mcpart->GetLabel(),pdg,mcpart->Pt(), eta, phi);
      
      photonlabel = iTr;

      int momidx = mcpart->GetMother();
      if(momidx>0){
        mom = static_cast<AliAODMCParticle*>(fAODMCParticles->At(momidx));
        mompdg = TMath::Abs(mom->GetPdgCode());
      }
      else
        mompdg=mcpart->GetPdgCode();
      
      // printf("With Mother at %d with label %d which is a %d",momidx, mom->GetLabel(), mompdg);
      
      eT = mcpart->E()*TMath::Sin(mcpart->Theta()); // Transform to transverse Energy
      
      fphietaPhotons->Fill(eta,phi,eT);
      
      bool foundmatch=kFALSE;
      for(int m=0;m<nTracks && foundmatch==kFALSE;m++){
        if(m==iTr)
	  continue;
        
        matchingtrack = static_cast<AliAODMCParticle*>(fAODMCParticles->At(m));
        
        if(! matchingtrack->IsPrimary())
	  continue;
        if(! matchingtrack->IsPhysicalPrimary())
	  continue;
        if(matchingtrack->GetStatus()> 10 )
	  continue;
        
        Double_t etamatching = matchingtrack->Eta();
        Double_t phimatching = matchingtrack->Phi();
        
        if(TMath::Abs(eta-etamatching)<=fdetacut && TMath::Abs(phi-phimatching)<=fdphicut){
          foundmatch=kTRUE;
          fphietaOthers->Fill(matchingtrack->Eta(),matchingtrack->Phi(),eT);
          fphietaOthersBis->Fill(matchingtrack->Eta(),matchingtrack->Phi(),matchingtrack->Pt());
        }
      }
      
      if(foundmatch)
	continue;
      
      distance=0.;
      phip=0., etap=0.;
      sumEiso=0.,sumUE=0.;
      
      for(int iTrack=0;iTrack<nTracks;iTrack++){
        if(iTrack==photonlabel)
          continue;
	
        mcpp = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));
        
        if(!mcpp)
	  continue;
        
        if(mcpp->Charge() != 0 && mcpp->GetStatus()>10)
          fPtTracksVSpTNC_MC->Fill(eT,mcpp->Pt());
        
        if(fIsoMethod==2){
          if((mcpp->Charge())==0)
	    continue;
	}
        
        if(mcpp->GetStatus()>10)
	  continue;
        
        int mumidx=mcpp->GetMother();
        if(mumidx<0 || mumidx>nTracks)
	  continue;
        
        mum = static_cast<AliAODMCParticle*>(fAODMCParticles->At(mumidx));
        if(mumidx == photonlabel || mum->GetPdgCode()==22)
	  continue;
        
        phip = mcpp->Phi();
        etap = mcpp->Eta();

	// Depending on which isolation method and UE method is considered
        distance = TMath::Sqrt((phi-phip)*(phi-phip) + (eta-etap)*(eta-etap));
        
        if(distance <= fIsoConeRadius){
	  // cout<<iTrack<<"\t"<<photonlabel<<endl;
	  // mcpp->Print();
          sumEiso += mcpp->E()*TMath::Sin(mcpp->Theta());
        }
        else
          AddParticleToUEMC(sumUE,mcpp, eta, phi);
      }

      CalculateUEDensityMC(sumUE);
      
      // printf("Storing Particle: Label %d  PDG: %d  Eta: %f  Phi: %f",mcpart->GetLabel(),pdg,eta,phi);
      // printf("With Mother at %d with label %d which is a %d",momidx, mom->GetLabel(), mompdg);
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
      // fill some histograms or a THnSparse or a TTree.
      //	AliError(Form("Fill something in Analize MC"));
      if(fWho==1)
	fOutMCTruth -> Fill(outputValuesMC);
    }
  }
  else{
    maxE=0.;
    int indexmaxE=0;

    // Getting the index of the particle with the maximum energy
    for(int iTr=0;iTr<nTracks;iTr++){
      mcsearch = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTr));
      
      if(!mcsearch)
	continue;
      
      if(mcsearch->GetStatus()>10)
	continue;
      if(mcsearch->GetPdgCode() != 22)
	continue;

      if(!fTPC4Iso){
        if((TMath::Abs(mcsearch->Eta())>(fGeom->GetArm1EtaMax()-0.03)-fIsoConeRadius ) || (mcsearch->Phi() < ((fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03+fIsoConeRadius) || mcsearch->Phi()>((fGeom->GetArm1PhiMax())*TMath::DegToRad()-0.03-fIsoConeRadius)))
          continue;
      }
      else{
        if((TMath::Abs(mcsearch->Eta())>0.87-fIsoConeRadius ) || (mcsearch->Phi() < ((fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03) || mcsearch->Phi()>((fGeom->GetArm1PhiMax())*TMath::DegToRad()-0.03)))
          continue;
      }

      mcfirstEnergy = mcsearch->E()*TMath::Sin(mcsearch->Theta());
      if(mcfirstEnergy>maxE){
        maxE=mcfirstEnergy;
        indexmaxE=iTr;
      }
      else
	continue;
    }

    mcfirst= static_cast<AliAODMCParticle*>(fAODMCParticles->At(indexmaxE));
    mcfirstEnergy=mcfirst->E()*TMath::Sin(mcfirst->Theta());

    int momidx= mcfirst->GetMother();
    if(momidx>0){
      mom = static_cast<AliAODMCParticle*>(fAODMCParticles->At(momidx));
      mompdg= TMath::Abs(mom->GetPdgCode());
    }
    else
      mompdg=mcfirst->GetPdgCode();
    
    mcFirstEta = mcfirst->Eta();
    mcFirstPhi = mcfirst->Phi();
    
    phip=0., etap=0.;
    sumEiso=0,sumUE=0;
    
    for(Int_t iTrack=1;iTrack<nTracks ;iTrack++){
      if(iTrack==indexmaxE)
	continue;

      mcpp = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));
      phip = mcpp->Phi();
      etap = mcpp->Eta();

      if(!mcpp)
        continue;
      
      if(mcpp->GetStatus()>10)
	continue;
      if(!mcpp->IsPrimary())
	continue;
      
      distance=0.;
      distance= TMath::Sqrt((mcFirstPhi- phip)*(mcFirstPhi- phip) + (mcFirstEta- etap)*(mcFirstEta- etap));
      
      if(distance<=fIsoConeRadius)
        sumEiso += mcpp->E()*TMath::Sin(mcpp->Theta());
      else
        AddParticleToUEMC(sumUE,mcpp,mcFirstEta,mcFirstPhi);
    }

    // cout<<"\n\nTotal Energy inside the Isolation Cone : "<<sumEiso<<endl;
    CalculateUEDensityMC(sumUE);
    // cout<<"Total UE Energy : "<<sumUE<<" calculated with method "<<fUEMethod<<endl;
    outputValuesMC[0] = mcfirstEnergy;
    outputValuesMC[1] = sumEiso;
    outputValuesMC[2] = sumUE;
    outputValuesMC[3] = mompdg;
    outputValuesMC[4] = mcFirstEta;
    outputValuesMC[5] = mcFirstPhi;
    outputValuesMC[6] = mcfirst->GetLabel();
    // Fill the Output TTree for MC Truth
    if(fWho==1)
      fOutMCTruth->Fill(outputValuesMC);
  }
  
  return;
}
