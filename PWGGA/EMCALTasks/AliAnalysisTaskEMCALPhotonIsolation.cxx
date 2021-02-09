  // $Id$
  //
  // Emcal Neutral Cluster analysis base task.
  //
  // Authors: D. Lodato, L. Ronflette, M. Marquard, E. Masson

// --- ROOT system ---
#include <Riostream.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <THnSparse.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TRandom3.h>

// --- Analysis system ---
#include "AliAnalysisTaskEMCALPhotonIsolation.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
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
#include "AliLog.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAnalysisUtils.h"
#include "AliEMCALTriggerPatchInfo.h"

// --- Detectors ---
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"

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
fDPMjetHeader(0),
fPythiaVersion(""),
fVariableCPV(kFALSE),
fVariableCPVInCone(kFALSE),
fVariableCPVBoth(kFALSE),
fVariableCPVSyst(""),
fEOverPMin(0.),
fEOverPMax(2000.),
fNonLinRecoEnergyScaling(kFALSE),
fExtraPerpConesFactor(1.363),
fTracksAna(0),
fStack(0),
fEMCALRecoUtils(new AliEMCALRecoUtils),
fWho(-1),
fSSsmearing(0),
fSSsmearwidth(0),
fSSsmear_mean(0),
fWhich(0),
fRejectPileUpEvent(kFALSE),
fNContrToPileUp(3),
// fOutputList(0),
fMinClusterEnergy(5.),
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
fMCtruth(0),
fPeriod(""),
fFiducialCut(0.4),
fAreasPerEvent(kFALSE),
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
fANnoSameTcard(kFALSE),
fBinsPt(0.),
fBinsM02(0.),
fBinsEtiso(0.),
fBinsEtue(0.),
fBinsEta(0.),
fBinsPhi(0.),
fBinsLabel(0.),
fBinsPDG(0.),
fBinsMomPDG(0.),
fBinsClustPDG(0.),
fBinsDx(0.),
fBinsDz(0.),
fBinsDecay(0.),
fTrackMult(0),
fTrackMultInCone(0),
fPtvsSum_MC(0),
fPtvsUE_MC(0),
fPtvsSumUE_MC(0),
fGenPromptPhotonSel(0),
fEtaPhiClus(0),
fEtaPhiClusAftSel(0),
fPtvsDetavsDphi(0),
fPtvsTrackPtvsDeta(0),
fPtvsTrackPtvsDphi(0),
fEOverPvsPt(0),
fEOverPvsPtWithCPV(0),
fClusEvsClusT(0),
fNCellsPerCluster(0),
fDTBCperCluster(0),
fPTbeforeNonLinScaling(0),
fPT(0),
fE(0),
fNLM(0),
fNLM2_NC_Acc(0),
fNLM2_NC_Acc_noTcard(0),
fVz(0),
fEvents(0),
fCutFlowEvents(0),
fCutFlowClusters(0),
fPtaftTime(0),
fPtaftCell(0),
fPtaftNLM(0),
fClusEtVsEtaPhiMatched(0),
fClusEtVsEtaPhiUnmatched(0),
fPtaftTM(0),
fPtaftDTBC(0),
fPtaftFC(0),
fPtaftM02C(0),
fClusTime(0),
fM02(0),
fEtaPhiClusVsM02(0),
fEtaPhiClusVsEtIsoClus(0),
fEtaPhiClusVsPtIsoTrack(0),
fEtaPhiClusVsPtUETrackCside(0),
fEtaPhiClusVsPtUETrackAside(0),
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
fEtaBandUENeutral_MC(0),
fEtaBandUECharged_MC(0),
fPerpConesUETracks(0),
fTPCWithoutIsoConeB2BbandUE(0),
// fNTotClus10GeV(0),
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
fPtVsNtrConeVsChgCone_Norm(0),
fPtVsConeVsUE(0),
fPtVsConeVsUE_Norm(0),
fEtaBandVsConeArea(0),
fPtVsConeVsEtaBand(0),
fPtvsM02vsSumUE_Norm(0),
fTestEtaPhiCone(0),
fInvMassM02iso(0),
fInvMassM02noiso(0),
fPtvsM02vsSum(0),
fPtvsM02vsSumUE(0),
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
fTracksPhiVsPt(0),
fTracksEtaVsPt(0),
fEtaPhiSPDRefit(0),
fEtaPhiNoSPDRefit(0), 
fTrackResolutionPtMC(0),
fVzBeforecut(0),
fOutputTHnS(0),
fOutMCTruth(0),
fOutClustMC(0),
fOutputQATree(0),
fOutputTree(0),
fphietaPhotons(0),
fphietaOthers(0),
fphietaOthersBis(0),
fSPDclustVsSPDtracklets(0),
fnPUevents(0),
f2012EGA(0),
fHistoRangeContainer(0x0)
  // tracks(0),
  // clusters(0)
{

    // Default constructor

    // fParticleCollArray.SetOwner(kTRUE);
    // for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;

  SetMakeGeneralHistograms(kTRUE);
  if(f2012EGA) SetCaloTriggerPatchInfoName("EmcalTriggers");
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
fDPMjetHeader(0),
fPythiaVersion(""),
fVariableCPV(kFALSE),
fVariableCPVInCone(kFALSE),
fVariableCPVBoth(kFALSE),
fVariableCPVSyst(""),
fEOverPMin(0.),
fEOverPMax(2000.),
fNonLinRecoEnergyScaling(kFALSE),
fExtraPerpConesFactor(1.363),
fTracksAna(0),
fStack(0),
fEMCALRecoUtils(new AliEMCALRecoUtils),
fWho(-1),
fSSsmearing(0),
fSSsmearwidth(0),
fSSsmear_mean(0),
fWhich(0),
fRejectPileUpEvent(kFALSE),
fNContrToPileUp(3),
// fOutputList(0),
fMinClusterEnergy(5.),
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
fMCtruth(0),
fPeriod(""),
fFiducialCut(0.4),
fAreasPerEvent(kFALSE),
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
fANnoSameTcard(kFALSE),
fBinsPt(0.),
fBinsM02(0.),
fBinsEtiso(0.),
fBinsEtue(0.),
fBinsEta(0.),
fBinsPhi(0.),
fBinsLabel(0.),
fBinsPDG(0.),
fBinsMomPDG(0.),
fBinsClustPDG(0.),
fBinsDx(0.),
fBinsDz(0.),
fBinsDecay(0.),
fTrackMult(0),
fTrackMultInCone(0),
fPtvsSum_MC(0),
fPtvsUE_MC(0),
fPtvsSumUE_MC(0),
fGenPromptPhotonSel(0),
fEtaPhiClus(0),
fEtaPhiClusAftSel(0),
fPtvsDetavsDphi(0),
fPtvsTrackPtvsDeta(0),
fPtvsTrackPtvsDphi(0),
fEOverPvsPt(0),
fEOverPvsPtWithCPV(0),
fClusEvsClusT(0),
fNCellsPerCluster(0),
fDTBCperCluster(0),
fPTbeforeNonLinScaling(0),
fPT(0),
fE(0),
fNLM(0),
fNLM2_NC_Acc(0),
fNLM2_NC_Acc_noTcard(0),
fVz(0),
fEvents(0),
fCutFlowEvents(0),
fCutFlowClusters(0),
fPtaftTime(0),
fPtaftCell(0),
fPtaftNLM(0),
fClusEtVsEtaPhiMatched(0),
fClusEtVsEtaPhiUnmatched(0),
fPtaftTM(0),
fPtaftDTBC(0),
fPtaftFC(0),
fPtaftM02C(0),
fClusTime(0),
fM02(0),
fEtaPhiClusVsM02(0),
fEtaPhiClusVsEtIsoClus(0),
fEtaPhiClusVsPtIsoTrack(0),
fEtaPhiClusVsPtUETrackCside(0),
fEtaPhiClusVsPtUETrackAside(0),
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
fEtaBandUENeutral_MC(0),
fEtaBandUECharged_MC(0),
fPerpConesUETracks(0),
fTPCWithoutIsoConeB2BbandUE(0),
// fNTotClus10GeV(0),
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
fPtVsNtrConeVsChgCone_Norm(0),
fPtVsConeVsUE(0),
fPtVsConeVsUE_Norm(0),
fEtaBandVsConeArea(0),
fPtVsConeVsEtaBand(0),
fPtvsM02vsSumUE_Norm(0),
fTestEtaPhiCone(0),
fInvMassM02iso(0),
fInvMassM02noiso(0),
fPtvsM02vsSum(0),
fPtvsM02vsSumUE(0),
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
fTracksPhiVsPt(0),
fTracksEtaVsPt(0),
fEtaPhiSPDRefit(0),
fEtaPhiNoSPDRefit(0), 
fTrackResolutionPtMC(0),
fVzBeforecut(0),
fOutputTHnS(0),
fOutMCTruth(0),
fOutClustMC(0),
fOutputQATree(0),
fOutputTree(0),
fphietaPhotons(0),
fphietaOthers(0),
fphietaOthersBis(0),
fSPDclustVsSPDtracklets(0),
fnPUevents(0),
f2012EGA(0),
fHistoRangeContainer(0x0)
  // tracks(0),
  // clusters(0)
{

    // Standard constructor

    // fParticleCollArray.SetOwner(kTRUE);
    //  for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;

  SetMakeGeneralHistograms(kTRUE);
  if(f2012EGA) SetCaloTriggerPatchInfoName("EmcalTriggers");
}

  //________________________________________________________________________
AliAnalysisTaskEMCALPhotonIsolation::~AliAnalysisTaskEMCALPhotonIsolation(){

    // Destructor

  if ( AliAnalysisManager::GetAnalysisManager()->IsProofMode() ) return;

  if ( fOutput ) delete fOutput;

}


  //________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::UserCreateOutputObjects(){

    // Create ouput histograms and THnSparse and TTree

  AliAnalysisTaskEmcal::UserCreateOutputObjects();

    // Histogram binning and ranges (to consider geometry-dependent histograms)
  Int_t   netabins = GetHistogramRangesAndBinning()->GetHistoEtaBins();
  Float_t etamax   = GetHistogramRangesAndBinning()->GetHistoEtaMax();
  Float_t etamin   = GetHistogramRangesAndBinning()->GetHistoEtaMin();

  Int_t   nphibins = GetHistogramRangesAndBinning()->GetHistoPhiBins();
  Float_t phimax   = GetHistogramRangesAndBinning()->GetHistoPhiMax();
  Float_t phimin   = GetHistogramRangesAndBinning()->GetHistoPhiMin();

  if((fIsoMethod == 0 || fIsoMethod == 1 || fIsoMethod == 3) && fTPC4Iso){
    cout<<"Error: Iso_Methods with CELLS and CLUSTERS work only within EMCal "<<endl;
    cout<<"Please Set Iso_Method and TPC4Iso Accordingly!!"<<endl;
    return;
  }

  if((fIsoMethod == 0 || fIsoMethod == 1 || fIsoMethod == 3) && (fUEMethod > 1 && fUEMethod < 4)){
    cout<<"Error: UE_Methods with CELLS and CLUSTERS work only within EMCal"<<endl;
    cout<<"Please Set Iso_Method and UE_Method Accordingly!!"<<endl;
    return;
  }

  if((fUEMethod > 1 && fUEMethod < 4) && !fTPC4Iso){
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
  else if(fUEMethod==4)
    sUEMethod = "ExtraCones";

  if(fTPC4Iso)
    sBoundaries = "TPC Acceptance";
  else
    sBoundaries = "EMCal Acceptance";

  if(fWho>2 || fWho==-1){
    cout<<"Error!!! OutputMode Can Only Be 0: TTree; 1: THnSparse; 2: TH**"<<endl;
    return;
  }
  else{
    fOutput = new AliEmcalList();
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

        sTitle = Form("Direct Photons: #it{p}_{T} , #sigma_{long}^{2} , E_{T} Iso%s in %s, E_{T} UE %s in %s, #eta_{clus} distr,#phi_{clus} distr; #it{p}_{T} (GeV/c); #sigma_{long}^{2}; E_{T}^{iso%s} (GeV/c) ; E_{T}^{UE%s} (GeV/c); #eta_{cl}; #phi_{cl}", sIsoMethod.Data(), sBoundaries.Data(), sUEMethod.Data(), sBoundaries.Data(), sIsoMethod.Data(), sUEMethod.Data());

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

            // fOutMCTruth = new THnSparseF ("fOutMCTruth","E_{#gamma}, E_{T}^{iso cone}, E_{T}^{UE}, MomPDG, Eta, Phi, Label; E_{T}^{#gamma} (GeV/c); #it{p}_{T}^{Iso}(GeV/c);E_{T} ^{UE} (GeV/c); PDG; #eta; #phi; Label",7,binsMC,xminbis,xmaxbis);
            // fOutMCTruth->Sumw2();
            // fOutput->Add(fOutMCTruth);
            //

          fOutMCTruth = new THnSparseF ("fOutMCTruth","E_{#gamma}, E_{T}^{iso cone}, E_{T}^{UE}, MomPDG, Eta, Phi, Label; E_{T}^{#gamma} (GeV/c); #it{p}_{T}^{Iso}(GeV/c);E_{T} ^{UE} (GeV/c); PDG; #eta; #phi; Label",ndimsMC,binsMC);
          fOutMCTruth->SetBinEdges(0,fBinsPt.data());
          fOutMCTruth->SetBinEdges(1,fBinsEtiso.data());
          fOutMCTruth->SetBinEdges(2,fBinsEtue.data());
          fOutMCTruth->SetBinEdges(3,fBinsPDG.data());
          fOutMCTruth->SetBinEdges(4,fBinsEta.data());
          fOutMCTruth->SetBinEdges(5,fBinsPhi.data());
          fOutMCTruth->SetBinEdges(6,fBinsLabel.data());
          fOutMCTruth->Sumw2();
          fOutput->Add(fOutMCTruth);


          fMCQAdim = sizeof(binsSMC)/sizeof(Int_t);
          const Int_t ndimsMCQA = fMCQAdim;

//          Double_t xminbismix[] = {0.,  0., -3000, -400,  0.,-1., -1., -10,    0.};
//          Double_t xmaxbismix[] = {70., 2.,  3000,  400, 70., 1.,  1., 100.,  10.};

            // fOutClustMC = new THnSparseF ("fOutClustMC", "E_{T}^{clust}, #sigma_{long}^{2}, PDG, MOM PDG, E_{T}^{true}, #Deltax, #Deltaz, E_{T}^{iso},Label;E_{T}^{reco} (GeV/c); #sigma_{long}^{2};PDG Code; Mothers' PDG Code; E_{T}^{MCtrue} (GeV/c); #Delta#phi; #Delta#eta; E_{T}^{iso} (Gev/c);Label",9,binsSMC,xminbismix,xmaxbismix);
          fOutClustMC = new THnSparseF ("fOutClustMC", "E_{T}^{clust}, #sigma_{long}^{2}, PDG, MOM PDG, E_{T}^{true}, #Deltax, #Deltaz, E_{T}^{iso},Label;E_{T}^{reco} (GeV/c); #sigma_{long}^{2};PDG Code; Mothers' PDG Code; E_{T}^{MCtrue} (GeV/c); #Delta#phi; #Delta#eta; E_{T}^{iso} (Gev/c);Label",ndimsMCQA,binsSMC);
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

          fTrackResolutionPtMC= new TH2F("hsigmaPTvsPT","track resolution",100.,0.,25., 500, -0.5,0.5);
          fTrackResolutionPtMC->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
          fTrackResolutionPtMC->GetYaxis()->SetTitle("#it{p}_{T}^{gen}-#it{p}_{T}^{reco} (Gev/c)");
          fTrackResolutionPtMC->Sumw2();
          fOutput->Add(fTrackResolutionPtMC);
        }

	if(fQA) {
          fphietaPhotons = new TH3F ("fDphiDeta_Photons","#Delta#phi vs #Delta#eta Clust-MCpart to check why zero #sigma_{long}^{2} clusters; #eta; #phi", 100, -0.5, 0.5, 200, 1.5, 3.5,60,0.,60.);
          fphietaPhotons->Sumw2();
          fOutput->Add(fphietaPhotons);
	/* Not filled anymore	
          fphietaOthers = new TH3F ("fDphiDeta_Others","#Delta#phi vs #Delta#eta Clust-MCpart to check why zero #sigma_{long}^{2} clusters; #eta; #phi", 140, -0.7, 0.7, 220, 0.8, 3.5,60,0.,60.);
          fphietaOthers->Sumw2();
          fOutput->Add(fphietaOthers);

          fphietaOthersBis = new TH3F ("fDphiDeta_OthersBis","#Delta#phi vs #Delta#eta Clust-MCpart to check why zero #sigma_{long}^{2} clusters; #eta; #phi", 140, -0.7, 0.7, 220, 0.8, 3.5,60,0.,60.);
          fphietaOthersBis->Sumw2();
          fOutput->Add(fphietaOthersBis);
	  */ 

	  fNLM2_NC_Acc = new TH2F("hNLM2_NC_Acc","NLM distribution for *Neutral* Clusters in acceptance",10,0.,10.,100,0.,100.);
	  fNLM2_NC_Acc->Sumw2();
	  fOutput->Add(fNLM2_NC_Acc);
    
	  if(fANnoSameTcard){
	    fNLM2_NC_Acc_noTcard = new TH2F("hNLM2_NC_Acc_noTcard","NLM distribution for *Neutral* Clusters in acceptance with NLM=2 NOT in same Tcard",10,0.,10.,100,0.,100.);
	    fNLM2_NC_Acc_noTcard->Sumw2();
	    fOutput->Add(fNLM2_NC_Acc_noTcard);
	  }
	}
      }
        break;

      case 2:
      {
          // Initialization TH*D/TH*F

	if(fIsMC && fNonLinRecoEnergyScaling){
	  fPTbeforeNonLinScaling = new TH1F("hPt_BeforeNonLinScaling","#it{p}_{T} distribution for clusters before candidate selection",100,0.,100.);
	  fPTbeforeNonLinScaling->Sumw2();
	  fOutput->Add(fPTbeforeNonLinScaling);
	}

	if(fQA){
	  fPtaftM02C = new TH1F("hPtaftM02C_NC","#it{p}_{T} distribution for Clusters after shower shape cut",200,0.,100.);
	  fPtaftM02C->Sumw2();
	  fOutput->Add(fPtaftM02C);

	  fM02 = new TH2F("hM02_NC","#sigma_{long}^{2} vs. #it{E}_{T} for clusters",100,0.,100.,500,0.,5.);
	  fM02->Sumw2();
	  fOutput->Add(fM02);

	  fEtaPhiClusVsM02 = new TH3F ("hEtaVsPhiVsM02", "#eta vs. #varphi vs. #sigma_{long}^{2} for clusters with 14 < #it{E}_{T} < 16 GeV", netabins, etamin, etamax, nphibins, phimin, phimax, 75, 0., 1.5);
	  fEtaPhiClusVsM02->Sumw2();
	  fOutput->Add(fEtaPhiClusVsM02);

	  fEtaPhiClusVsEtIsoClus = new TH3F ("hEtaVsPhiVsEtIsoClus", "#eta vs. #varphi vs. #Sigma #it{E}_{T}^{clus, cone} for clusters with 14 < #it{E}_{T} < 16 GeV", netabins, etamin, etamax, nphibins, phimin, phimax, 60, 0., 30.);
	  fEtaPhiClusVsEtIsoClus->Sumw2();
	  fOutput->Add(fEtaPhiClusVsEtIsoClus);

	  fEtaPhiClusVsPtIsoTrack = new TH3F ("hEtaVsPhiVsPtIsoTrack", "#eta vs. #varphi vs. #Sigma #it{p}_{T}^{track, cone} for clusters with 14 < #it{E}_{T} < 16 GeV", netabins, etamin, etamax, nphibins, phimin, phimax, 60, 0., 30.);
	  fEtaPhiClusVsPtIsoTrack->Sumw2();
	  fOutput->Add(fEtaPhiClusVsPtIsoTrack);
	}

        // fEtIsoClust = new TH2F("hEtIsoClus_NC","#Sigma #it{p}_{T}^{iso cone} in iso cone distribution for Neutral Clusters with EMCal Clusters",200,0.,100.,200,0.,100.);
        // fEtIsoClust->SetYTitle("#Sigma #it{p}_{T}^{iso cone} (GeV/c)");
        // fEtIsoClust->SetXTitle("#it{p}_{T}^{clust}");
        // fEtIsoClust->Sumw2();
        // fOutput->Add(fEtIsoClust);

        // fPtIsoTrack = new TH2F("hPtIsoTrack_NC"," #Sigma #it{p}_{T}^{iso cone} in iso cone distribution for Neutral Clusters with Tracks",200,0.,100.,200,0.,100.);
        // fPtIsoTrack->SetYTitle("#Sigma #it{p}_{T}^{iso cone} (GeV/c)");
        // fPtIsoTrack->SetXTitle("#it{p}_{T}^{clust}");
        // fPtIsoTrack->Sumw2();
        // fOutput->Add(fPtIsoTrack);

	if(!fAnalysispPb){
	  fPtvsM02vsSum = new TH3F("hPtvsM02vsSum","#it{p}_{T} vs #sigma_{long}^{2} vs  #Sigma E_{T}^{iso cone} distribution for isolated clusters",140,0.,70.,400,0.,4.,220,-10.,100.);
	  fPtvsM02vsSum->Sumw2();
	  fOutput->Add(fPtvsM02vsSum);
	}

	if(!fAreasPerEvent){
	  fPtvsM02vsSumUE = new TH3F("hPtvsM02vsSumUE","#it{p}_{T} vs #sigma_{long}^{2} vs  #Sigma E_{T}^{iso cone}-UE distribution for clusters",140,0.,70.,400,0.,4.,400,-50.,150.);
	  fPtvsM02vsSumUE->Sumw2();
	  fOutput->Add(fPtvsM02vsSumUE);
	}

	if(!fAreasPerEvent && fUEMethod == 1){
	  fPtVsConeVsEtaBand = new TH3F("hPtVsConeVsEtaBand","Cluster energy vs. cone energy vs. eta band energy (not normalised)",140,0.,70.,250,0.,100.,250,0.,100.);
	  fPtVsConeVsEtaBand->SetXTitle("#it{p}_{T}^{cluster}");
	  fPtVsConeVsEtaBand->SetYTitle("#sum^{cone} #it{p}_{T}");
	  fPtVsConeVsEtaBand->SetZTitle("#sum^{eta-band} #it{p}_{T}");
	  fPtVsConeVsEtaBand->Sumw2();
	  fOutput->Add(fPtVsConeVsEtaBand);
	}

	if(fAreasPerEvent){
	  if(fQA){
	    fEtaBandVsConeArea = new TH2F("hEtaBandVsConeArea","Eta-band vs. cone area (depending on the cluster position)", 140, 0.16, 0.51, 200, 0.30, 0.80);
	    fEtaBandVsConeArea->Sumw2();
	    fOutput->Add(fEtaBandVsConeArea);
	  }

	  if(fAnalysispPb){
	    fPtvsM02vsSumUE_Norm = new TH3F("hPtvsM02vsSumUE_Norm","#it{p}_{T} vs. #sigma_{long}^{2} vs. #Sigma E_{T}^{iso cone}-UE distribution (already normalised by the appropriate areas)",140,0.,70.,400,0.,4.,400,-50.,150.);
	    fPtvsM02vsSumUE_Norm->Sumw2();
	    fOutput->Add(fPtvsM02vsSumUE_Norm);
	  }
	}

	if(fIsoMethod==0){
	  fEtIsoCells = new TH1F("hEtIsoCell_NC","E_{T}^{iso cone} in iso cone distribution for Neutral Clusters with EMCal Cells",200,-0.25,99.75);
	  fEtIsoCells->SetXTitle("#Sigma E_{T}^{iso cone} (GeV/c)");
	  fEtIsoCells->Sumw2();
	  fOutput->Add(fEtIsoCells);

	  fEtIsolatedCells = new TH1F("hEtIsolatedCells","E_{T} distribution for Isolated Photons with cells; #Sigma E_{T}^{iso cone}<Ethres",100,0.,100.);
	  fEtIsolatedCells->SetXTitle("E_{T}^{iso}");
	  fEtIsolatedCells->Sumw2();
	  fOutput->Add(fEtIsolatedCells);

	  fPhiBandUECells = new TH2F(Form("hPhiBandUE_CELLS"),Form("UE Estimation with Phi Band CELLS"),140,0.,70.,250,0.,100.);
	  fPhiBandUECells->SetXTitle("E_{T}");
	  fPhiBandUECells->SetYTitle("#Sigma E_{T}^{UE}");
	  fPhiBandUECells->Sumw2();
	  fOutput->Add(fPhiBandUECells);

	  fEtaBandUECells = new TH2F(Form("hEtaBandUE_CELLS"),Form("UE Estimation with Eta Band CELLS"),140,0.,70.,250,0.,100.);
	  fEtaBandUECells->SetXTitle("E_{T}");
	  fEtaBandUECells->SetYTitle("#Sigma E_{T}^{UE}");
	  fEtaBandUECells->Sumw2();
	  fOutput->Add(fEtaBandUECells);
	}

	if(fUEMethod==0 && fAreasPerEvent){
	  fPhiBandUEClust = new TH2F(Form("hPhiBandUE_Cluster"),Form("UE Estimation with Phi Band Clusters (already normalised by the appropriate areas)"),140,0.,70.,250,0.,100.);
	  fPhiBandUEClust->SetXTitle("E_{T}");
	  fPhiBandUEClust->SetYTitle("#Sigma E_{T}^{UE}");
	  fPhiBandUEClust->Sumw2();
	  fOutput->Add(fPhiBandUEClust);

	  fPhiBandUETracks = new TH2F(Form("hPhiBandUE_TPC"),Form("UE Estimation with Phi Band Tracks (already normalised by the appropriate areas)"),140,0.,70.,250,0.,100.);
	  fPhiBandUETracks->SetXTitle("E_{T}");
	  fPhiBandUETracks->SetYTitle("#Sigma #it{p}_{T}^{UE}");
	  fPhiBandUETracks->Sumw2();
	  fOutput->Add(fPhiBandUETracks);
        }

	if(fUEMethod==1){
	  if(fAreasPerEvent){
	    fEtaBandUEClust = new TH2F(Form("hEtaBandUE_Cluster"),Form("UE Estimation with Eta Band Clusters (already normalised by the appropriate areas)"),140,0.,70.,250,0.,100.);
	    fEtaBandUEClust->SetXTitle("E_{T}");
	    fEtaBandUEClust->SetYTitle("#Sigma E_{T}^{UE}");
	    fEtaBandUEClust->Sumw2();
	    fOutput->Add(fEtaBandUEClust);

	    fEtaBandUETracks = new TH2F(Form("hEtaBandUE_TPC"),Form("UE Estimation with Eta Band Tracks (already normalised by the appropriate areas)"),140,0.,70.,250,0.,100.);
	    fEtaBandUETracks->SetXTitle("E_{T}");
	    fEtaBandUETracks->SetYTitle("#Sigma #it{p}_{T}^{UE}");
	    fEtaBandUETracks->Sumw2();
	    fOutput->Add(fEtaBandUETracks);
	  }

	  if(fIsMC){
	    fEtaBandUENeutral_MC = new TH1F("hEtaBandUE_Neutral_MC", "Neutral UE estimation with Eta Band (generated)",250,0.,100.);
	    fEtaBandUENeutral_MC->Sumw2();
	    fOutput->Add(fEtaBandUENeutral_MC);

	    fEtaBandUECharged_MC = new TH1F("hEtaBandUE_Charged_MC", "Charged UE Estimation with Eta Band (generated)",250,0.,100.);
	    fEtaBandUECharged_MC->Sumw2();
	    fOutput->Add(fEtaBandUECharged_MC);
	  }
        }

	if(fUEMethod==2 || fUEMethod==4){
	  fPerpConesUETracks = new TH2F("hConesUE","UE Estimation with Perpendicular Cones in TPC",140,0.,70.,250,0.,100.);
	  fPerpConesUETracks->SetXTitle("E_{T}");
	  fPerpConesUETracks->SetYTitle("#Sigma #it{p}_{T}^{UE}");
	  fPerpConesUETracks->Sumw2();
	  fOutput->Add(fPerpConesUETracks);
	}

	if(fUEMethod==3){
	  fTPCWithoutIsoConeB2BbandUE = new TH2F("hFullTPCUE","UE Estimation with almost Full TPC",140,0.,70.,250,0.,100.);
	  fTPCWithoutIsoConeB2BbandUE->SetXTitle("E_{T}");
	  fTPCWithoutIsoConeB2BbandUE->SetYTitle("#Sigma E_{T}^{UE}");
	  fTPCWithoutIsoConeB2BbandUE->Sumw2();
	  fOutput->Add(fTPCWithoutIsoConeB2BbandUE);
        }

	//   fPtEtIsoTC = new TH1F("hPtEtIsoTrackClust_NC","#Sigma #it{p}_{T}^{iso cone} + #Sigma E_{T}^{iso cone} in iso cone distribution for Neutral Clusters with Tracks and Clusters",200,-0.25,99.75);
	//   fPtEtIsoTC->SetXTitle("#Sigma #it{p}_{T}^{iso cone} + #Sigma E_{T}^{iso cone} (GeV/c)");
	//   fPtEtIsoTC->Sumw2();
	//   fOutput->Add(fPtEtIsoTC);

	//   fEtIsolatedClust = new TH1F("hEtIsolatedClust","E_{T} distribution for Isolated Photons with clusters; #Sigma E_{T}^{iso cone}<Ethres",140,0.,70.);
	//   fEtIsolatedClust->SetXTitle("E_{T}^{iso}");
	//   fEtIsolatedClust->Sumw2();
	//   fOutput->Add(fEtIsolatedClust);

	//   fPtIsolatedNClust = new TH1F("hEtIsolatedNClust","#it{p}_{T} distribution for neutral clusters; #Sigma #it{p}_{T}^{iso cone}<Pthres",140,0.,70.);
	//   fPtIsolatedNClust->SetXTitle("#it{p}_{T}^{iso}");
	//   fPtIsolatedNClust->Sumw2();
	//   fOutput->Add(fPtIsolatedNClust);

	//   fPtIsolatedNTracks = new TH1F("hEtIsolatedNTracks","#it{p}_{T} distribution for neutral clusters; #Sigma #it{p}_{T}^{iso cone}<Pthres",140,0.,70.);
	//   fPtIsolatedNTracks->SetXTitle("#it{p}_{T}^{iso}");
	//   fPtIsolatedNTracks->Sumw2();
	//   fOutput->Add(fPtIsolatedNTracks);

	//   fEtIsolatedTracks = new TH1F("hEtIsolatedTracks","E_{T} distribution for Isolated Photons with tracks; #Sigma #it{p}_{T}^{iso cone}<Pthres",100,0.,100.);
	//   fEtIsolatedTracks->SetXTitle("E_{T}^{iso}");
	//   fEtIsolatedTracks->Sumw2();
	//   fOutput->Add(fEtIsolatedTracks);

	//   fPtvsM02iso = new TH2F("hPtvsM02iso","#it{p}_{T} vs #sigma_{long}^{2} distribution for isolated clusters",140,0.,70.,500,0.,5.);
	//   fPtvsM02iso->SetXTitle("#it{p}_{T}^{iso}");
	//   fPtvsM02iso->SetYTitle("#sigma_{long}^{2}");
	//   fOutput->Add(fPtvsM02iso);

	//   fPtvsM02noiso = new TH2F("hPtvsM02noiso","#it{p}_{T} vs #sigma_{long}^{2} distribution for non isolated clusters",140,0.,70.,500,0.,5.);
	//   fPtvsM02noiso->SetXTitle("#it{p}_{T}^{iso}");
	//   fPtvsM02noiso->SetYTitle("#sigma_{long}^{2}");
	//   fOutput->Add(fPtvsM02noiso);

	//   fEtaPhiClusVsPtUETrackCside = new TH3F ("hEtaVsPhiVsPtUETrack_Cside", "#eta vs. #varphi vs. #Sigma #it{p}_{T}^{track, UE} (tracks in C side) for clusters with 14 < #it{E}_{T} < 16 GeV", netabins, etamin, etamax, nphibins, phimin, phimax, 60, 0., 30.);
	//   fEtaPhiClusVsPtUETrackCside->Sumw2();
	//   fOutput->Add(fEtaPhiClusVsPtUETrackCside);

	//   fEtaPhiClusVsPtUETrackAside = new TH3F ("hEtaVsPhiVsPtUETrack_Aside", "#eta vs. #varphi vs. #Sigma #it{p}_{T}^{track, UE} (tracks in A side) for clusters with 14 < #it{E}_{T} < 16 GeV", netabins, etamin, etamax, nphibins, phimin, phimax, 60, 0., 30.);
	//   fEtaPhiClusVsPtUETrackAside->Sumw2();
	//   fOutput->Add(fEtaPhiClusVsPtUETrackAside);

	if(fIsMC){
	  if(!fAnalysispPb){
	    fPtvsSum_MC = new TH2F("hPtvsSum_MC","#it{p}_{T} vs #Sigma E_{T}^{iso cone} distribution for isolated particles",140,0.,70.,400,-50.,150.);
	    fPtvsSum_MC->Sumw2();
	    fOutput->Add(fPtvsSum_MC);
	  }

	  if(fAnalysispPb){
	    fPtvsUE_MC = new TH2F("hPtvsUE_MC","Total UE energy vs. particle #it{p}_{T}",140,0.,70.,250,0.,100.);
	    fPtvsUE_MC->Sumw2();
	    fOutput->Add(fPtvsUE_MC);

	    fPtvsSumUE_MC = new TH2F("hPtvsSumUE_MC","#it{p}_{T} vs #Sigma E_{T}^{iso cone}-UE distribution for isolated particles (already normalised by the appropriate areas)",140,0.,70.,400,-50.,150.);
	    fPtvsSumUE_MC->Sumw2();
	    fOutput->Add(fPtvsSumUE_MC);
	  }

	  fGenPromptPhotonSel = new TH1F ("hGenPromptPhotonSel","Generated direct photon selection cuts", 7, 0., 7.);
	  fGenPromptPhotonSel->GetXaxis()->SetBinLabel(1, "All part.");
	  fGenPromptPhotonSel->GetXaxis()->SetBinLabel(2, "Final state");
	  fGenPromptPhotonSel->GetXaxis()->SetBinLabel(3, "Phys. prim.");
	  fGenPromptPhotonSel->GetXaxis()->SetBinLabel(4, "Photons");
	  fGenPromptPhotonSel->GetXaxis()->SetBinLabel(5, "Mother photon");
	  fGenPromptPhotonSel->GetXaxis()->SetBinLabel(6, "Prompt photon");
	  fGenPromptPhotonSel->GetXaxis()->SetBinLabel(7, "Acceptance");
	  fGenPromptPhotonSel->Sumw2();
	  fOutput->Add(fGenPromptPhotonSel);
	}

	if(fQA){
	  fPtvsDetavsDphi = new TH3F("hPtvsDetavsDphi","Cluster-track matching vs. cluster energy", 120, 0., 60., 200, -0.05, 0.05, 200, -0.05, 0.05);
	  fPtvsDetavsDphi->Sumw2();
	  fOutput->Add(fPtvsDetavsDphi);
	}

	fPtvsTrackPtvsDeta = new TH3F("hPtvsTrackPtvsDeta","Cluster-track matching #Delta#eta vs. track #it{p}_{T} vs. cluster energy", 120, 0., 60., 100, 0., 100., 200, -0.05, 0.05);
	fPtvsTrackPtvsDeta->Sumw2();
	fOutput->Add(fPtvsTrackPtvsDeta);

	fPtvsTrackPtvsDphi = new TH3F("hPtvsTrackPtvsDphi","Cluster-track matching #Delta#varphi vs. track #it{p}_{T} vs. cluster energy", 120, 0., 60., 100, 0., 100., 200, -0.05, 0.05);
	fPtvsTrackPtvsDphi->Sumw2();
	fOutput->Add(fPtvsTrackPtvsDphi);

	if(fQA){
	  fEOverPvsPt = new TH2F("hEOverPvsPt","Cluster #it{E} over track #it{p} vs. cluster #it{p}_{T} BEFORE CPV", 120, 0., 60., 300, 0., 15.);
	  fEOverPvsPt->Sumw2();
	  fOutput->Add(fEOverPvsPt);

	  fEOverPvsPtWithCPV = new TH2F("hEOverPvsPt_WithCPV","Cluster #it{E} over track #it{p} vs. cluster #it{p}_{T} AFTER CPV", 120, 0., 60., 300, 0., 15.);
	  fEOverPvsPtWithCPV->Sumw2();
	  fOutput->Add(fEOverPvsPtWithCPV);
	}
      }
      break;
    }
  }

    // Initialize the common QA histograms
  if(fQA){
      // Include QA plots to the OutputList
    fClusTime = new TH1F("hClusTime_NC","Time distribution for Clusters",800,-50.,50.);
    fClusTime->Sumw2();
    fOutput->Add(fClusTime);

    if(fWho != 2){
      fDeltaETAClusTrack = new TH1F("h_Dz","Track-Cluster Dz ",1000,-0.5,0.5);
      fDeltaETAClusTrack->Sumw2();
      fOutput->Add(fDeltaETAClusTrack);

      fDeltaPHIClusTrack = new TH1F("h_Dx","Track-Cluster Dx",1000,-0.5,0.5);
      fDeltaPHIClusTrack->Sumw2();
      fOutput->Add(fDeltaPHIClusTrack);

      fDeltaETAClusTrackMatch = new TH1F("h_DzMatch","Track-Cluster Dz matching ",100,-0.05,0.05);
      fDeltaETAClusTrackMatch ->Sumw2();
      fOutput->Add(fDeltaETAClusTrackMatch);

      fDeltaPHIClusTrackMatch = new TH1F("h_DxMatch","Track-Cluster Dx matching",100,-0.05,0.05);
      fDeltaPHIClusTrackMatch->Sumw2();
      fOutput->Add(fDeltaPHIClusTrackMatch);

      fE = new TH1F("hE_NC","E distribution for Clusters",200,0.,100.);
      fE->Sumw2();
      fOutput->Add(fE);
    }

    if(fWho != 2){
      fTestIndex= new TH2F("hTestIndex","Test index for cluster",100,0.,100.,100,0.,100.);
      fTestIndex->SetXTitle("index");
      fTestIndex->SetYTitle("local index");
      fTestIndex->Sumw2();
      fOutput->Add(fTestIndex);

      fTestIndexE= new TH2F("hTestIndexE","Test index vs energy for cluster",200,0.,100.,100,0.,100.);
      fTestIndexE->SetXTitle("cluster energy");
      fTestIndexE->SetYTitle("index");
      fTestIndexE->Sumw2();
      fOutput->Add(fTestIndexE);

      fTestLocalIndexE= new TH2F("hTestLocalIndexE","Test local index vs energy for cluster",200,0.,100.,100,0.,100.);
      fTestLocalIndexE->SetXTitle("cluster energy");
      fTestLocalIndexE->SetYTitle("local index");
      fTestLocalIndexE->Sumw2();
      fOutput->Add(fTestLocalIndexE);

      fTestEtaPhiCone= new TH2F("hTestEtatPhiCone","Test eta phi neutral clusters candidates",100,0,TMath::TwoPi(),netabins, etamin, etamax);
      fTestEtaPhiCone->SetXTitle("phi");
      fTestEtaPhiCone->SetYTitle("eta");
      fTestEtaPhiCone->Sumw2();
      fOutput->Add(fTestEtaPhiCone);

      fEtVSM02VSPisotrack = new TH3F ("hEtVSM02VSPisotrack","Energy clust vs cluster #sigma_{long}^{2} vs Charged iso",70,0.,700.,200,0.,2.,400,0.,100.);
      fEtVSM02VSPisotrack->SetXTitle("Cluster E (GeV/c)");
      fEtVSM02VSPisotrack->SetYTitle("cluster #sigma_{long}^{2}");
      fEtVSM02VSPisotrack->SetZTitle("Isolation Charged (GeV/c)");
      fEtVSM02VSPisotrack->Sumw2();
      fOutput->Add(fEtVSM02VSPisotrack);

      fEtVSM02VSEisoclust = new TH3F ("hEtVSM02VSEisoclust","Energy clust vs cluster #sigma_{long}^{2} vs Neutral iso",70,0.,700.,200,0.,2.,400,0.,100.);
      fEtVSM02VSEisoclust->SetXTitle("Cluster E (GeV/c)");
      fEtVSM02VSEisoclust->SetYTitle("cluster #sigma_{long}^{2}");
      fEtVSM02VSEisoclust->SetZTitle("Isolation Neutrals (GeV/c)");
      fEtVSM02VSEisoclust->Sumw2();
      fOutput->Add(fEtVSM02VSEisoclust);
    }

      // fInvMassM02iso = new TH3F("hInvMassM02iso","Invariant mass vs #sigma_{long}^{2} vs E_{T}^{iso cluster}",100,0.,1.,500,0.,5.,200,0.,100.);
      // fInvMassM02iso->Sumw2();
      // fOutput->Add(fInvMassM02iso);

      // fInvMassM02noiso = new TH3F("hInvMassM02noiso","Invariant mass vs #sigma_{long}^{2} vs E_{T}^{no iso cluster}",100,0.,1.,500,0.,5.,200,0.,100.);
      // fInvMassM02noiso->Sumw2();
      // fOutput->Add(fInvMassM02noiso);

    fClusEtVsEtaPhiMatched = new TH3F ("hEtaVsPhiVsEt_Matched", "#eta vs. #varphi vs. #it{E}_{T} for TRACK-MATCHED clusters", netabins, etamin, etamax, nphibins, phimin, phimax, 60, 0., 60.);
    fClusEtVsEtaPhiMatched->Sumw2();
    fOutput->Add(fClusEtVsEtaPhiMatched);

    fClusEtVsEtaPhiUnmatched = new TH3F ("hEtaVsPhiVsEt_Unmatched", "#eta vs. #varphi vs. #it{E}_{T} for NON TRACK-MATCHED clusters", netabins, etamin, etamax, nphibins, phimin, phimax, 60, 0., 60.);
    fClusEtVsEtaPhiUnmatched->Sumw2();
    fOutput->Add(fClusEtVsEtaPhiUnmatched);
  }

    // Initialization of all the common THistos for the 3 different outputs

  fNLM = new TH2F("hNLM_NC","Number of local maxima per cluster vs. cluster #it{p}_{T}",140, 0., 70., 30, 0., 30.);
  fNLM->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fNLM->SetYTitle("#it{n}_{LM}");
  fNLM->Sumw2();
  fOutput->Add(fNLM);

  fVz = new TH1F("hVz_NC","Vertex Z distribution",100,-50.,50.);
  fVz->Sumw2();
  fOutput->Add(fVz);

  fVzBeforecut = new TH1F("hVz_ALL", "Inclusive Vertex Z distribution",100,-50.,50.);
  fVzBeforecut->Sumw2();
  fOutput->Add(fVzBeforecut);

  fEvents = new TH1F("hEvents_NC","Events",1,0.,1.);
  fEvents->Sumw2();
  fOutput->Add(fEvents);

  fCutFlowEvents = new TH1F("hCutFlowEvents", "Effect of each cut on event number", 4, 0., 4.);
  fCutFlowEvents->GetXaxis()->SetBinLabel(1, "Initial");
  fCutFlowEvents->GetXaxis()->SetBinLabel(2, "Vertex cut");
  fCutFlowEvents->GetXaxis()->SetBinLabel(3, "SPD pile-up cut");
  fCutFlowEvents->GetXaxis()->SetBinLabel(4, "Trackless cut");
  fCutFlowEvents->Sumw2();
  fOutput->Add(fCutFlowEvents);

  fCutFlowClusters = new TH1F("hCutFlowClusters", "Effect of each cut on candidate cluster number", 9, 0., 9.);
  fCutFlowClusters->GetXaxis()->SetBinLabel(1, "Initial");
  fCutFlowClusters->GetXaxis()->SetBinLabel(2, "IsEMCal cut");
  fCutFlowClusters->GetXaxis()->SetBinLabel(3, "Time cut");
  fCutFlowClusters->GetXaxis()->SetBinLabel(4, "Number of cells cut");
  fCutFlowClusters->GetXaxis()->SetBinLabel(5, "NLM cut");
  fCutFlowClusters->GetXaxis()->SetBinLabel(6, "DTBC cut");
  fCutFlowClusters->GetXaxis()->SetBinLabel(7, "CPV cut");
  fCutFlowClusters->GetXaxis()->SetBinLabel(8, "Fiducial cut");
  fCutFlowClusters->GetXaxis()->SetBinLabel(9, "5 GeV lower cut");
  fCutFlowClusters->Sumw2();
  fOutput->Add(fCutFlowClusters);

  fTrackMult = new TH1F ("hTrackMult","Tracks multiplicity Distribution",100,0.,100.);
  fTrackMult->Sumw2();
  fOutput->Add(fTrackMult);

  fTrackMultInCone = new TH1F ("hTrackMultInCone","Tracks multiplicity in isolation cone (per candidate)",100,0.,100.);
  fTrackMultInCone->Sumw2();
  fOutput->Add(fTrackMultInCone);

  fEtaPhiClus = new TH2F ("hEtaPhiClusActivity", "", netabins, etamin, etamax, nphibins, phimin, phimax);
    // fEtaPhiClus->Sumw2();
  fOutput->Add(fEtaPhiClus);

  fEtaPhiClusAftSel = new TH2F ("hEtaPhiClusAfterSelection", "", netabins, etamin, etamax, nphibins, phimin, phimax);
    // fEtaPhiClusAftSel->Sumw2();
  fOutput->Add(fEtaPhiClusAftSel);

  if(fWho != 2){
    fSPDclustVsSPDtracklets = new TH2F("hSPDclustVsSPDtracklets","Number of SPD clusters VS number of SPD tracklets in events with |Zvtx| < 10",100,0,200,250,0,1000);
    fSPDclustVsSPDtracklets->Sumw2();
    fOutput->Add(fSPDclustVsSPDtracklets);

    fnPUevents = new TH1F("hnPUevents","Number of pileUp events rejected", 2,0.,1.);
    fnPUevents->Sumw2();
    fOutput->Add(fnPUevents);
  }

  fClusEvsClusT = new TH2F("fClustTimeVSClustEn", "Cluster time vs. cluster #it{p}_{T}", 200, 0., 100., 800, -200., 200.);
  fClusEvsClusT->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fClusEvsClusT->SetYTitle("#it{t}_{clus} (ns)");
  fClusEvsClusT->Sumw2();
  fOutput->Add(fClusEvsClusT);

  fNCellsPerCluster = new TH2F ("hNCellsPerCluster","Number of cells per cluster vs. cluster #it{p}_{T}", 140, 0., 70., 60, 0., 60.); 
  fNCellsPerCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fNCellsPerCluster->SetYTitle("#it{n}_{cells}");
  fNCellsPerCluster->Sumw2();
  fOutput->Add(fNCellsPerCluster);

  fDTBCperCluster = new TH2F ("hDTBCperCluster","Distance to bad channel per cluster vs. cluster #it{p}_{T}", 140, 0., 70., 60, 0., 60.); 
  fDTBCperCluster->SetXTitle("#it{p}_{T} (GeV/#it{c})");
  fDTBCperCluster->SetYTitle("#it{d}_{BC}");
  fDTBCperCluster->Sumw2();
  fOutput->Add(fDTBCperCluster);

  fPT = new TH1F("hPt_NC","#it{p}_{T} distribution for clusters before candidate selection",200,0.,100.);
  fPT->Sumw2();
  fOutput->Add(fPT);

  fPtaftTime = new TH1F("hPtaftTime_NC","#it{p}_{T} distribution for Clusters after cluster time cut",200,0.,100.);
  fPtaftTime->Sumw2();
  fOutput->Add(fPtaftTime);

  fPtaftCell = new TH1F("hPtaftCell_NC","#it{p}_{T} distribution for Clusters after Ncells cut",200,0.,100.);
  fPtaftCell->Sumw2();
  fOutput->Add(fPtaftCell);

  fPtaftNLM = new TH1F("hPtaftNLM_NC","#it{p}_{T} distribution for Clusters after NLM cut",200,0.,100.);
  fPtaftNLM->Sumw2();
  fOutput->Add(fPtaftNLM);

  fPtaftTM = new TH1F("hPtaftTM_NC","#it{p}_{T} distribution for Neutral Clusters",200,0.,100.);
  fPtaftTM->Sumw2();
  fOutput->Add(fPtaftTM);

  fPtaftDTBC = new TH1F("hPtaftDTBC_NC","#it{p}_{T} distribution for Neutral Clusters after DTBC cut",200,0.,100.);
  fPtaftDTBC->Sumw2();
  fOutput->Add(fPtaftDTBC);

  fPtaftFC = new TH1F("hPtaftFC_NC","#it{p}_{T} distribution for Clusters after fiducial cut",200,0.,100.);
  fPtaftFC->Sumw2();
  fOutput->Add(fPtaftFC);

  if(!fAreasPerEvent){
    fTestEnergyCone = new TH3F("hTestEnergyConeVSpT","Test energy clusters and tracks in cone",140,0.,70.,250,0.,100.,200,0.,100.);
    fTestEnergyCone->SetXTitle("#it{p}_{T}^{cluster}");
    fTestEnergyCone->SetYTitle("#sum^{cone} #it{p}_{T}^{cluster}");
    fTestEnergyCone->SetZTitle("#sum^{cone} #it{p}_{T}^{track}");
    fTestEnergyCone->Sumw2();
    fOutput->Add(fTestEnergyCone);
  }

  fPtVsNtrConeVsChgCone_Norm = new TH3F("hPtVsNtrConeVsChgCone_Norm","Neutral energy in cone vs. charged energy in cone vs. cluster #it{p}_{T} (already normalised by the appropriate areas)",140,0.,70.,250,0.,100.,200,0.,100.);
  fPtVsNtrConeVsChgCone_Norm->SetXTitle("#it{p}_{T}^{cluster}");
  fPtVsNtrConeVsChgCone_Norm->SetYTitle("#sum^{cone} #it{p}_{T}^{cluster}");
  fPtVsNtrConeVsChgCone_Norm->SetZTitle("#sum^{cone} #it{p}_{T}^{track}");
  fPtVsNtrConeVsChgCone_Norm->Sumw2();
  fOutput->Add(fPtVsNtrConeVsChgCone_Norm);

  if(!fAreasPerEvent){
    fPtVsConeVsUE = new TH3F("hPtVsConeVsUE","Total energy in cone (before area norm) vs. total UE energy (before area norm) vs. cluster #it{p}_{T}",140,0.,70.,250,0.,100.,200,0.,100.);
    fPtVsConeVsUE->SetXTitle("#it{p}_{T}^{cluster}");
    fPtVsConeVsUE->SetYTitle("#it{p}_{T}^{iso}");
    fPtVsConeVsUE->SetZTitle("#it{p}_{T}^{UE}");
    fPtVsConeVsUE->Sumw2();
    fOutput->Add(fPtVsConeVsUE);
  }

  fPtVsConeVsUE_Norm = new TH3F("hPtVsConeVsUE_Norm","Total energy in cone (before area norm) vs. total UE energy (before area norm) vs. cluster #it{p}_{T} (already normalised by the appropriate areas)",140,0.,70.,250,0.,100.,200,0.,100.);
  fPtVsConeVsUE_Norm->SetXTitle("#it{p}_{T}^{cluster}");
  fPtVsConeVsUE_Norm->SetYTitle("#it{p}_{T}^{iso}");
  fPtVsConeVsUE_Norm->SetZTitle("#it{p}_{T}^{UE}");
  fPtVsConeVsUE_Norm->Sumw2();
  fOutput->Add(fPtVsConeVsUE_Norm);

  // fTracksConeEtaPt = new TH3F("hTracksConeEtaPt","#Sigma vs #eta vs E_{T}",200,0.,100.,320,-0.8,0.8,200,0.,100.);
  // fTracksConeEtaPt->Sumw2();
  // fOutput->Add(fTracksConeEtaPt);

  // fTracksConeEtaM02 = new TH3F("hTracksConeEtaM02","#Sigma vs #eta vs #sigma_{long}^{2}",200,0.,100.,320,-0.8,0.8,500,0.,5.);
  // fTracksConeEtaM02->Sumw2();
  // fOutput->Add(fTracksConeEtaM02);

  // fphietaPhotons = new TH3F("hphietaPhotons","Test eta phi photons MC",250,-0.8,0.8, 250, 1.2, 3.4,200,0.,1.);
  // fOutput->Add(fphietaPhotons);

  // fphietaOthers = new TH3F("hphietaOthers","Test eta phi others",250,-0.8,0.8, 250, 1.2, 3.4,200,0.,1.);
  // fOutput->Add(fphietaOthers);

  fPtTracksVSpTNC = new TH2F ("hTrackPtSpecVSpT","Charged Particle spectrum vs pT Candidate",70,0.,70.,100,0.,20.);
  fPtTracksVSpTNC->Sumw2();
  fOutput->Add(fPtTracksVSpTNC);

  if(fWho != 2){
    fPhiTracksVSclustPt  = new TH2F("hPhiTracks_vs_clustPT","Tracks phi distr vs pT Candidate",70, 0.,70., 200,0.,TMath::TwoPi());
    fPhiTracksVSclustPt->Sumw2();
    fOutput->Add(fPhiTracksVSclustPt);

    fEtaTracksVSclustPt  = new TH2F("hEtaTracks_vs_clustPT","Tracks eta distr vs pT Candidate",70, 0.,70., 90,-0.9,0.9);
    fEtaTracksVSclustPt->Sumw2();
    fOutput->Add(fEtaTracksVSclustPt);
  }

  if(fQA){
    fTracksPhiVsPt = new TH2F("hTracksPhiVsPt", "Tracks #varphi vs #it{p}_{T}", 40, 0., 40., 200, 0., TMath::TwoPi());
    fTracksPhiVsPt->Sumw2();
    fOutput->Add(fTracksPhiVsPt);

    fTracksEtaVsPt = new TH2F("hTracksEtaVsPt", "Tracks #eta vs #it{p}_{T}", 40, 0., 40., 180, -0.9, 0.9);
    fTracksEtaVsPt->Sumw2();
    fOutput->Add(fTracksEtaVsPt);
  }

  fEtaPhiSPDRefit = new TH2F ("hEtaPhiSPDRefit","#eta vs #varphi of tracks with SPD and ITS refit", 180, -0.9, 0.9, 200, 0., TMath::TwoPi());
  fEtaPhiSPDRefit->SetXTitle("#eta");
  fEtaPhiSPDRefit->SetYTitle("#varphi (rad)");
  fOutput->Add(fEtaPhiSPDRefit);

  fEtaPhiNoSPDRefit = new TH2F ("hEtaPhiNoSPDRefit","#eta vs #varphi of tracks without SPD and with ITS refit", 180, -0.9, 0.9, 200, 0., TMath::TwoPi());
  fEtaPhiNoSPDRefit->SetXTitle("#eta");
  fEtaPhiNoSPDRefit->SetYTitle("#varphi (rad)");
  fOutput->Add(fEtaPhiNoSPDRefit);

  if(fWho != 2){
    fCTdistVSpTNC = new TH2F ("hDistanceC_TrackVSpT","Distance between Neutral Clust and closest Track vs pT Candidate",70,0.,70.,210,-0.1,2.);
    fCTdistVSpTNC->Sumw2();
    fOutput->Add(fCTdistVSpTNC);
  }

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

    if(fWho != 2 && fQA){
      fpi0VSclusterVSIsolation = new TH3F ("hpi0VSclusterVSisolation","Energy pi0 vs cluster Energy vs Isolation",95,5.,100.,95,5.,100.,400,0.,100.);
      fpi0VSclusterVSIsolation->SetXTitle("particle (#pi^{0} or #eta) E");
      fpi0VSclusterVSIsolation->SetYTitle("cluster E");
      fpi0VSclusterVSIsolation->SetZTitle("Isolation");
      fpi0VSclusterVSIsolation->Sumw2();
      fOutput->Add(fpi0VSclusterVSIsolation);

      fpi0VSM02VSIsolation = new TH3F ("fpi0VSM02VSIsolation","Energy pi0 vs cluster #sigma_{long}^{2} vs Isolation",95,5.,100.,500,0.,5.,400,0.,100.);
      fpi0VSM02VSIsolation->SetXTitle("particle (#pi^{0} or #eta) E");
      fpi0VSM02VSIsolation->SetYTitle("cluster #sigma_{long}^{2}");
      fpi0VSM02VSIsolation->SetZTitle("Isolation");
      fpi0VSM02VSIsolation->Sumw2();
      fOutput->Add(fpi0VSM02VSIsolation);

      fpi0VSclusterVSM02 = new TH3F ("fpi0VSclusterVSM02","Energy pi0 vs Energy cluster vs cluster #sigma_{long}^{2} ",95,5.,100.,100,0.,100.,500,0.,5.);
      fpi0VSclusterVSM02->SetXTitle("particle (#pi^{0} or #eta) E");
      fpi0VSclusterVSM02->SetYTitle("cluster E");
      fpi0VSclusterVSM02->SetZTitle("cluster #sigma_{long}^{2}");
      fpi0VSclusterVSM02->Sumw2();
      fOutput->Add(fpi0VSclusterVSM02);
    }
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
  if(!tracks){
    AliError(Form("%s: This task needs a 1particle container!", GetName()));
    return;
  }

    // Tracks for Isolation
  AliTrackContainer *tracksANA = GetTrackContainer("filterTracksAna");
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
  coi->GetMomentum(vecCOI, fVertex, AliVCluster::kNonLinCorr);
  if(fIsMC && fNonLinRecoEnergyScaling)
    NonLinRecoEnergyScaling(vecCOI);

  Double_t coiTOF = coi->GetTOF()*1e9;
  index=coi->GetID();
  if(coi->GetM02()>=0.1)
    fClusEvsClusT->Fill(vecCOI.Pt(), coiTOF);
  if(!fIsMC){
    if(coiTOF< -30. || coiTOF > 30.)
      return kFALSE;
  }

  fPtaftTime->Fill(vecCOI.Pt());
  fCutFlowClusters->Fill(2.5);
  fNCellsPerCluster->Fill(vecCOI.Pt(), coi->GetNCells());

  if((coi->GetNCells() < 2))
    return kFALSE;

  fPtaftCell->Fill(vecCOI.Pt());
  fCutFlowClusters->Fill(3.5);

  Int_t nlm = 0;
  AliVCaloCells * fCaloCells = InputEvent()->GetEMCALCells();
  if(fCaloCells){
    nlm = GetNLM(coi,fCaloCells);
    AliDebug(1,Form("NLM = %d",nlm));

    fNLM->Fill(vecCOI.Pt(), nlm);

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
  fCutFlowClusters->Fill(4.5);
  fDTBCperCluster->Fill(vecCOI.Pt(), coi->GetDistanceToBadChannel());

  if(coi->GetDistanceToBadChannel() < 2)
    return kFALSE;

  fPtaftDTBC->Fill(vecCOI.Pt());
  fCutFlowClusters->Fill(5.5);
  if(fQA && fWho != 2)
    fTestIndexE->Fill(vecCOI.Pt(),index);

  if(fTMClusterRejected){
    if(ClustTrackMatching(coi,kTRUE)){
      if(fQA)
	fClusEtVsEtaPhiMatched->Fill(vecCOI.Eta(), vecCOI.Phi(), vecCOI.Pt());
      return kFALSE;
    }
    else{
      if(fQA)
	fClusEtVsEtaPhiUnmatched->Fill(vecCOI.Eta(), vecCOI.Phi(), vecCOI.Pt());
    }
  }

  fPtaftTM->Fill(vecCOI.Pt());
  fCutFlowClusters->Fill(6.5);

  if(!CheckBoundaries(vecCOI))
    return kFALSE;

  fPtaftFC->Fill(vecCOI.Pt());
  fCutFlowClusters->Fill(7.5);

  if(vecCOI.Pt() < fMinClusterEnergy)
    return kFALSE;

  fCutFlowClusters->Fill(8.5);

  if(fQA && fWho == 1){
    fNLM2_NC_Acc->Fill(nlm, NonLinRecoEnergyScaling(coi, coi->GetNonLinCorrEnergy()));
    
    if(fANnoSameTcard){
      if(nlm==2){
        Int_t rowdiff,coldiff;
        if(!IsAbsIDsFromTCard(fAbsIDNLM[0],fAbsIDNLM[1],rowdiff,coldiff)){
          fNLM2_NC_Acc_noTcard->Fill(nlm, NonLinRecoEnergyScaling(coi, coi->GetNonLinCorrEnergy()));
          return kTRUE;
        }
        else
          return kFALSE;
      }
      else
	fNLM2_NC_Acc_noTcard->Fill(nlm, NonLinRecoEnergyScaling(coi, coi->GetNonLinCorrEnergy()));
    }
  } //the flag fANnoSameTcard could be used also independently of fQA and fWho,
    //depending on the results of the analysis on full dataset.
  
  return kTRUE;
}

  //______________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::Run()
{
    // Run the analysis

  AliTrackContainer *tracks = GetTrackContainer("tpconlyMatch");
  if(!tracks){
    AliWarning(Form("Cannot find the tracks for CT matching"));
    return kFALSE;
  }

  if(tracks->GetTrackFilterType() != AliEmcalTrackSelection::kTPCOnlyTracks){
    AliWarning(Form("CT matching NOT performed with TPCOnly tracks"));
    AliWarning(Form("You better be sure of what you are doing"));
  }

  AliTrackContainer *tracksANA = GetTrackContainer("filterTracksAna");
  if(!tracksANA){
    AliWarning(Form("Cannot find the tracks for isolation"));
    return kFALSE;
  }

  // AliWarning(Form("FilterType of the tracks for Analysis: %d \t(should be %d)", tracksANA->GetTrackFilterType(),AliEmcalTrackSelection::kHybridTracks));

  if(tracksANA->GetTrackFilterType() != AliEmcalTrackSelection::kHybridTracks){
    AliWarning(Form("Isolation NOT calculated with HybridTracks"));
    AliWarning(Form("You better be sure of what you are doing"));
  }

  fVevent = dynamic_cast<AliVEvent*>(InputEvent());

  fCutFlowEvents->Fill(0.5);
  fVzBeforecut->Fill(fVertex[2]);

  if(fVertex[2]>10. || fVertex[2]<-10.)
    return kFALSE;

  fCutFlowEvents->Fill(1.5);

  if(fWho != 2)
    fSPDclustVsSPDtracklets->Fill(fVevent->GetMultiplicity()->GetNumberOfTracklets(),(fVevent->GetNumberOfITSClusters(0)+fVevent->GetNumberOfITSClusters(1)));

  if(fRejectPileUpEvent && fVevent->IsPileupFromSPD(fNContrToPileUp, 0.8, 3., 2., 5.)){
    if(fWho != 2)
      fnPUevents->Fill(0);

    return kFALSE;
  }

  fCutFlowEvents->Fill(2.5);

  AliClusterContainer* clusters = GetClusterContainer(0);
  Int_t nbTracksEvent;
  nbTracksEvent = InputEvent()->GetNumberOfTracks();

  if(fRejectionEventWithoutTracks && (nbTracksEvent == 0))
    return kFALSE;

  fCutFlowEvents->Fill(3.5);

  // Reject events below 2012 EGA threshold
  Bool_t isL1 = kFALSE;
  if(f2012EGA){
    TClonesArray *triPatchInfo = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcalTriggers"));
    if(triPatchInfo){
      Int_t nPatch = triPatchInfo->GetEntries();
      for(Int_t ip = 0;ip<nPatch;ip++){
        AliEMCALTriggerPatchInfo *pti = static_cast<AliEMCALTriggerPatchInfo*>(triPatchInfo->At(ip));
        if(!pti) continue;
        if(!pti->IsEMCal()) continue;
        if(!pti->IsRecalcGamma()) continue;
	// if(pti->GetEtaGeo()==0. && pti->GetPhiGeo()==0.){
	//   continue;
	// }
        if(!fIsMC && pti->GetADCAmp() > 130){
          isL1 = kTRUE;
          break;
        }
        if(fIsMC && pti->GetSmearedEnergy() > 8.4){
          isL1 = kTRUE;
          break;
        }
      }
    }
  }
  if(!isL1 && !fIsMC && f2012EGA){
    return kFALSE;
  }

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
  }

  Int_t index = 0;

  if(fIsMC){
    fAODMCParticles = static_cast <TClonesArray*  >(InputEvent()->FindListObject(AliAODMCParticle::StdBranchName()));
    fmcHeader       = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

    if(!fIsMC)
      return kFALSE;

    if(!fStack && !fAODMCParticles){
      AliError("No MC stack saved");
      return kFALSE;
    }

    if(fAODMCParticles->GetEntries() < 1){ // DO THIS ALSO FOR ESDs?
      AliError("Number of MC particles insufficient");
      return kFALSE;
    }

    if(fPythiaVersion.Contains("8")){

      // Pythia8 // 201 = qg to qgamma // 202 = qqbar to ggamma
      if(fMCtruth || fPythiaHeader->ProcessType() == 201 || fPythiaHeader->ProcessType() == 202) AnalyzeMC_Pythia8();

    }
    else{

      // Pythia6 // 14 = qqbar to ggamma // 29 = qg to qgamma
      if(fMCtruth || fmcHeader->GetEventType() == 14 || fmcHeader->GetEventType() == 29) AnalyzeMC();

    }
  }

  if(fIsMC && (fTriggerLevel1 != 0)){
    if(!MCSimTrigger(fVevent,fTriggerLevel1) && fAnalysispPb)
      return kFALSE;
  }

  fVz->Fill(fVertex[2]); // Fill Vertex Z histogram

  if(!isL1 && fIsMC && f2012EGA){
    return kFALSE;
  }

  if ( fisLCAnalysis ) {
    // Get the event leading cluster
    AliVCluster *coi = (clusters->GetLeadingCluster());

    if(!coi){
      AliError("No cluster found");
      return kFALSE;
    }

    fCutFlowClusters->Fill(0.5);

    if(!coi->IsEMCAL())
      return kFALSE;

    fCutFlowClusters->Fill(1.5);

    index=coi->GetID();
    TLorentzVector vecCOI;
    coi->GetMomentum(vecCOI, fVertex, AliVCluster::kNonLinCorr);

    if(fIsMC && fNonLinRecoEnergyScaling){
        if(fWho==2){
	    fPTbeforeNonLinScaling->Fill(vecCOI.Pt());
	}
        NonLinRecoEnergyScaling(vecCOI);
    }

    fPT->Fill(vecCOI.Pt());

    if(fQA)
      FillQAHistograms(coi,vecCOI);

    fEtaPhiClus->Fill(vecCOI.Eta(),vecCOI.Phi());

    Bool_t isSelected = SelectCandidate(coi);

    if(isSelected){
      fEtaPhiClusAftSel->Fill(vecCOI.Eta(),vecCOI.Phi());

      for(auto it : tracksANA->accepted()){
	AliVTrack *tr = static_cast<AliVTrack*>(it);
	if(!tr){
	  AliError("No track found");
	  return kFALSE;
	}

	fPtTracksVSpTNC->Fill(vecCOI.Pt(),tr->Pt());
	if(fWho != 2){
	  fPhiTracksVSclustPt->Fill(vecCOI.Pt(),tr->Phi());
	  fEtaTracksVSclustPt->Fill(vecCOI.Pt(),tr->Eta());
	}
	if(fQA){
	  fTracksPhiVsPt->Fill(tr->Pt(),tr->Phi());
	  fTracksEtaVsPt->Fill(tr->Pt(),tr->Eta());
	}

	// Extracting the SPD status of hybrid tracks
	Bool_t        bITSRefit    = (tr->GetStatus() & AliVTrack::kITSrefit) == AliVTrack::kITSrefit;
	Bool_t        bConstrained = kFALSE;

	AliAODTrack * aodTrack     = dynamic_cast<AliAODTrack*>(tr);
	AliESDtrack * esdTrack     = dynamic_cast<AliESDtrack*>(tr);
	if      (aodTrack) bConstrained = aodTrack->IsGlobalConstrained();
	else if (esdTrack) bConstrained = (!esdTrack->HasPointOnITSLayer(0) && !esdTrack->HasPointOnITSLayer(1));

	if (bConstrained) {      
	  if (bITSRefit)
	    fEtaPhiNoSPDRefit->Fill(tr->Eta(), tr->Phi());
	}
	else {
	  if (bITSRefit)
	    fEtaPhiSPDRefit->Fill(tr->Eta(), tr->Phi());
	}
	// Extracting the SPD status of hybrid tracks (end)

      }

      FillGeneralHistograms(coi,vecCOI,index);
    }

  }
  else {
    // Get the entries of the Cluster Container
    for(auto it : clusters->accepted()){
      AliVCluster *coi = static_cast<AliVCluster*>(it);

      if(!coi){
	AliError("No cluster found");
	return kFALSE;
      }

      fCutFlowClusters->Fill(0.5);

      if(!coi->IsEMCAL())
	return kFALSE;

      fCutFlowClusters->Fill(1.5);

      index=coi->GetID();
      TLorentzVector vecCOI;
      coi->GetMomentum(vecCOI, fVertex, AliVCluster::kNonLinCorr);

      if(fIsMC && fNonLinRecoEnergyScaling){
          if(fWho==2){
	  fPTbeforeNonLinScaling->Fill(vecCOI.Pt());
	  }

	  NonLinRecoEnergyScaling(vecCOI);
      }

      fPT->Fill(vecCOI.Pt());

      if(fQA)
	FillQAHistograms(coi,vecCOI);

      fEtaPhiClus->Fill(vecCOI.Eta(),vecCOI.Phi());

      Bool_t isSelected = SelectCandidate(coi);

      if(isSelected){
	fEtaPhiClusAftSel->Fill(vecCOI.Eta(),vecCOI.Phi());

	for(auto it : tracksANA->accepted()){
	  AliVTrack *tr = static_cast<AliVTrack*>(it);
	  if(!tr){
	    AliError("No track found");
	    return kFALSE;
	  }

	  fPtTracksVSpTNC->Fill(vecCOI.Pt(),tr->Pt());
	  if(fWho != 2){
	    fPhiTracksVSclustPt->Fill(vecCOI.Pt(),tr->Phi());
	    fEtaTracksVSclustPt->Fill(vecCOI.Pt(),tr->Eta());
	  }
	  if(fQA){
	    fTracksPhiVsPt->Fill(tr->Pt(),tr->Phi());
	    fTracksEtaVsPt->Fill(tr->Pt(),tr->Eta());
	  }

	  // Extracting the SPD status of hybrid tracks
	  Bool_t        bITSRefit    = (tr->GetStatus() & AliVTrack::kITSrefit) == AliVTrack::kITSrefit;
	  Bool_t        bConstrained = kFALSE;

	  AliAODTrack * aodTrack     = dynamic_cast<AliAODTrack*>(tr);
	  AliESDtrack * esdTrack     = dynamic_cast<AliESDtrack*>(tr);
	  if      (aodTrack) bConstrained = aodTrack->IsGlobalConstrained();
	  else if (esdTrack) bConstrained = (!esdTrack->HasPointOnITSLayer(0) && !esdTrack->HasPointOnITSLayer(1));

	  if (bConstrained) {      
	    if (bITSRefit)
	      fEtaPhiNoSPDRefit->Fill(tr->Eta(), tr->Phi());
	  }
	  else {
	    if (bITSRefit)
	      fEtaPhiSPDRefit->Fill(tr->Eta(), tr->Phi());
	  }
	  // Extracting the SPD status of hybrid tracks (end)

	}

	FillGeneralHistograms(coi,vecCOI,index);
      }
    }
  }

  return kTRUE;
}

  //_________________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::FillQAHistograms(AliVCluster *coi, TLorentzVector vecCOI){

  Double_t m02COI = 0.;
  if(fSSsmearing)
    ApplySmearing(coi, m02COI);
  else
    m02COI = coi->GetM02();

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
      if(fQA)
	fM02->Fill(vecCOI.Pt(), m02COI);
      if(fQA && vecCOI.Pt()>14. && vecCOI.Pt()<16.)
	fEtaPhiClusVsM02->Fill(vecCOI.Eta(),vecCOI.Phi(),m02COI);
      break;
  }

  if(fWho != 2)
    fE->Fill(vecCOI.E());

  Double_t checktof = coi->GetTOF()*1e9;
  fClusTime->Fill(checktof);

  if(checktof>-30. && checktof<30. && !fIsMC){
    if(CheckBoundaries(vecCOI)){
      Double_t checkM02 = coi->GetM02();
      if(fM02mincut < checkM02 && checkM02 < fM02maxcut && fWho==2 && fQA)
        fPtaftM02C->Fill(vecCOI.Pt());
    }
  }
}

  //___________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::MCSimTrigger(AliVEvent *eventIn, Int_t triggerLevel){

  if(triggerLevel < 1 || triggerLevel > 2){
    AliError(Form("No trigger level has been chosen for the MC analysis"));
    return kFALSE;
  }

  Double_t threshold = 0.;
  Double_t spread    = 0.;

  Int_t runNumber = eventIn->GetRunNumber();
  if(!runNumber)
    return kFALSE;
  if(runNumber < 195180 || runNumber > 197469) // LHC13a to LHC13f (to be replaced by a condition on fAnalysispPb?)
    return kFALSE;

  // TString fired = InputEvent()->GetFiredTriggerClasses();
  // AliError(Form("Trigger used in the events %s",fired));

  if(triggerLevel == 1){
    threshold = 11.40; // Value obtained by fit on data (13d, 13e and 13f separately and averaged)
    spread    = 0.50;
  }

  if(triggerLevel == 2){
    threshold = 7.03;
    spread    = 0.30;
  }

  if(spread != 0.){
    TF1* triggerSmearing = new TF1("triggerSmearing","[0]*exp(-0.5*((x-[1])/[2])**2)", 0., 15.);
    triggerSmearing->SetParameter(0, 1./(spread*TMath::Sqrt(2.*TMath::Pi())));
    triggerSmearing->SetParameter(1, threshold);
    triggerSmearing->SetParameter(2, spread);
    threshold = triggerSmearing->GetRandom();
    delete triggerSmearing;
  }

  AliClusterContainer *clusters = GetClusterContainer(0);
  AliVCluster         *coi      = 0x0;

  for(auto it : clusters->accepted()){
    coi = static_cast<AliVCluster*>(it);

    if(!coi)
      continue;

    if(NonLinRecoEnergyScaling(coi, coi->GetNonLinCorrEnergy()) > threshold)
      return kTRUE;
  }

  return kFALSE;
}

  //__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::ClustTrackMatching(AliVCluster *clust, Bool_t candidate){

    // Check if the cluster matches a track

  AliTrackContainer* tracks = GetTrackContainer("tpconlyMatch");
  if(tracks->GetTrackFilterType() != AliEmcalTrackSelection::kTPCOnlyTracks)
    AliError(Form("No TPC only tracks"));

  TLorentzVector vecClust;
  clust->GetMomentum(vecClust, fVertex, AliVCluster::kNonLinCorr);
  if(fIsMC && fNonLinRecoEnergyScaling)
    NonLinRecoEnergyScaling(vecClust);

  Int_t nbMObj = clust->GetNTracksMatched();
  if(nbMObj == 0)
    return kFALSE;

  Double_t distCT = 0., maxdist = 10.;
  Bool_t matched = kFALSE;
  AliVTrack* mt = 0;

  for(Int_t i = 0; i < nbMObj && matched == kFALSE; i++){
    if(fIsEsd){
      Int_t imt = clust->GetTrackMatchedIndex(0);
      if(imt >= 0)
	mt = static_cast<AliVTrack*>(tracks->GetAcceptParticle(imt));
    }
    else{
      mt = static_cast<AliVTrack*>(clust->GetTrackMatched(i));
      UInt_t rejectionReason = 0;
      if(!tracks->AcceptParticle(mt, rejectionReason))
	mt = 0;
    }

    if(!mt)
      continue;

    // Variable CPV parameters (inspired by direct photons analysis)
    Float_t eta_param      [3] = {0.010, 4.07, -2.50}; Float_t phi_param      [3] = {0.015, 3.65, -2.00};
    Float_t eta_param_cone [3] = {0.010, 4.07, -2.50}; Float_t phi_param_cone [3] = {0.015, 3.65, -2.00};

    if      ( fVariableCPVSyst == "loose" ) {
      if ( !fVariableCPVInCone ) {
	eta_param[0] = 0.010; eta_param[1] = 4.78; eta_param[2] = -2.50;
	phi_param[0] = 0.015; phi_param[1] = 3.92; phi_param[2] = -2.00;
      }
      eta_param_cone[0] = 0.010; eta_param_cone[1] = 4.78; eta_param_cone[2] = -2.50;
      phi_param_cone[0] = 0.015; phi_param_cone[1] = 3.92; phi_param_cone[2] = -2.00;
    }
    else if ( fVariableCPVSyst == "tight" ) {
      if ( !fVariableCPVInCone ) {
	eta_param[0] = 0.015; eta_param[1] = 3.46; eta_param[2] = -2.50;
	phi_param[0] = 0.020; phi_param[1] = 3.73; phi_param[2] = -1.75;
      }
      eta_param_cone[0] = 0.015; eta_param_cone[1] = 3.46; eta_param_cone[2] = -2.50;
      phi_param_cone[0] = 0.020; phi_param_cone[1] = 3.73; phi_param_cone[2] = -1.75;
    }
    else {
      eta_param[0] = 0.010; eta_param[1] = 4.07; eta_param[2] = -2.50;
      phi_param[0] = 0.015; phi_param[1] = 3.65; phi_param[2] = -2.00;

      eta_param_cone[0] = 0.010; eta_param_cone[1] = 4.07; eta_param_cone[2] = -2.50;
      phi_param_cone[0] = 0.015; phi_param_cone[1] = 3.65; phi_param_cone[2] = -2.00;
    }

    Double_t deltaEta, deltaPhi;
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

    if(fWho != 2 && fQA && candidate){
      fDeltaETAClusTrack->Fill(deta);
      fDeltaPHIClusTrack->Fill(dphi);
    }

    if(candidate && fWho == 2){
      if(fQA){
	fPtvsDetavsDphi->Fill(vecClust.Pt(), deta, dphi);
	fEOverPvsPt->Fill(vecClust.Pt(), clust->GetNonLinCorrEnergy()/mt->P());
      }

      fPtvsTrackPtvsDeta->Fill(vecClust.Pt(), mt->Pt(), deta);
      fPtvsTrackPtvsDphi->Fill(vecClust.Pt(), mt->Pt(), dphi);
    }

    distCT = TMath::Sqrt(deta*deta+dphi*dphi);
    if(distCT < maxdist)
      maxdist = distCT;

    if(candidate){
      if((fVariableCPVBoth || fVariableCPV) || (fVariableCPVBoth && fVariableCPV)){
	deltaEta = eta_param[0] + TMath::Power((mt->Pt() + eta_param[1]), eta_param[2]);
	deltaPhi = phi_param[0] + TMath::Power((mt->Pt() + phi_param[1]), phi_param[2]);
      }
      else{
	deltaEta = fdetacut;
	deltaPhi = fdphicut;
      }
    }
    else{
      if((fVariableCPVBoth || fVariableCPVInCone) || (fVariableCPVBoth && fVariableCPVInCone)){
	deltaEta = eta_param_cone[0] + TMath::Power((mt->Pt() + eta_param_cone[1]), eta_param_cone[2]);
	deltaPhi = phi_param_cone[0] + TMath::Power((mt->Pt() + phi_param_cone[1]), phi_param_cone[2]);
      }
      else{
	deltaEta = fdetacutIso;
	deltaPhi = fdphicutIso;
      }
    }

    // Spatial match (and E/p cut if tight values set by the user)
    if(TMath::Abs(dphi) < deltaPhi && TMath::Abs(deta) < deltaEta){
      if(fQA && candidate && fWho != 2){
	fDeltaETAClusTrackMatch->Fill(deta);
	fDeltaPHIClusTrackMatch->Fill(dphi);
      }

      if(candidate && fWho == 2 && fQA)
	fEOverPvsPtWithCPV->Fill(vecClust.Pt(), clust->GetNonLinCorrEnergy()/mt->P());

      // E/p cut (extremely loose by default, 0-2000)
      if ( (clust->GetNonLinCorrEnergy()/mt->P() >= fEOverPMin) && (clust->GetNonLinCorrEnergy()/mt->P() < fEOverPMax) )
	matched = kTRUE;
    }
  }

  if(fWho != 2)
    fCTdistVSpTNC->Fill(vecClust.Pt(),distCT);

  return matched;
}

  //_____________________________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::NonLinRecoEnergyScaling ( TLorentzVector cluster_vec ) {

  Int_t    cluster_SM = 0;
  fGeom->SuperModuleNumberFromEtaPhi(cluster_vec.Eta(), cluster_vec.Phi(), cluster_SM);  

  Int_t    maxSM_withoutTRD       = 4;
  Double_t scaleFactor_withoutTRD = 1., scaleFactor_withTRD = 1.;

  if      ( fPeriod.Contains("11") ) {
    maxSM_withoutTRD       = 6;
    scaleFactor_withoutTRD = 1./1.01217;
    scaleFactor_withTRD    = 1./0.99994;
  }
  else if ( fPeriod.Contains("12") ) {
    scaleFactor_withoutTRD = 1./1.00142;
    scaleFactor_withTRD    = 1./0.995343;
  }
  else if ( fPeriod.Contains("13") ) {
    scaleFactor_withoutTRD = 1./0.997278;
    scaleFactor_withTRD    = 1./0.993112;
  }

  if ( cluster_SM < maxSM_withoutTRD )
    cluster_vec *= scaleFactor_withoutTRD; // Scale each component of the 4-vector
  else
    cluster_vec *= scaleFactor_withTRD;

}

  //_____________________________________________________________________________________________
Double_t AliAnalysisTaskEMCALPhotonIsolation::NonLinRecoEnergyScaling ( AliVCluster * cluster, Double_t energy ) {

  if ( fIsMC && fNonLinRecoEnergyScaling ) {
    Double_t scaled_energy          = 0.;
    Int_t    cluster_SM             = GetModuleNumber(cluster);
    Int_t    maxSM_withoutTRD       = 4;
    Double_t scaleFactor_withoutTRD = 1., scaleFactor_withTRD = 1.;

    if      ( fPeriod.Contains("11") ) {
      maxSM_withoutTRD       = 6;
      scaleFactor_withoutTRD = 1./1.01217;
      scaleFactor_withTRD    = 1./0.99994;
    }
    else if ( fPeriod.Contains("12") ) {
      scaleFactor_withoutTRD = 1./1.00142;
      scaleFactor_withTRD    = 1./0.995343;
    }
    else if ( fPeriod.Contains("13") ) {
      scaleFactor_withoutTRD = 1./0.997278;
      scaleFactor_withTRD    = 1./0.993112;
    }

    if ( cluster_SM < maxSM_withoutTRD )
      scaled_energy = scaleFactor_withoutTRD * energy;
    else
      scaled_energy = scaleFactor_withTRD * energy;

    return scaled_energy;
  }
  else
    return energy;

}

  //_____________________________________________________________________________________________
Int_t AliAnalysisTaskEMCALPhotonIsolation::GetModuleNumber ( AliVCluster * cluster ) const {

    // Adapted from AliCalorimeterUtils

  if ( !cluster ) {
    AliDebug(1,"AliCalorimeterUtils::GetModuleNumber() - NUL Cluster, please check!!!");
    return -1;
  }
  
  if ( cluster->GetNCells() <= 0 ) return -1;
  
  Int_t absId = cluster->GetCellAbsId(0);
  
  if ( absId < 0 ) return -1;
  
  if( cluster->IsEMCAL() ) {
    AliDebug(2,Form("EMCAL absid %d, SuperModule %d",absId, fGeom->GetSuperModuleNumber(absId)));
    return fGeom->GetSuperModuleNumber(absId) ;
  }
	
  return -1;
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

  Int_t iDigitN = 0 , iDigit = 0 ;
  Int_t absId1  = -1, absId2 = -1;
  const Int_t nCells = coi->GetNCells();

  Float_t eCluster = RecalEnClust(coi, cells); // Recalculate cluster energy, avoid non lin. correction
  Float_t localMaxCutE = 0.1, locMaxCutEDiff = 0.;

  Float_t en    = 0., emax = 0., en1 = 0., en2 = 0.;
  Int_t   idmax = -1;

  for(iDigit = 0; iDigit < nCells ; iDigit++){
    absIdList[iDigit] = coi->GetCellsAbsId()[iDigit];
    en = cells->GetCellAmplitude(absIdList[iDigit]);
    RecalAmpCell(en,absIdList[iDigit]);

    if(en > emax){
      emax  = en ;
      idmax = absIdList[iDigit] ;
    }
  }

  for(iDigit = 0 ; iDigit < nCells; iDigit++){
    if(absIdList[iDigit] >= 0){
      absId1 = coi->GetCellsAbsId()[iDigit];
      en1 = cells->GetCellAmplitude(absId1);
      RecalAmpCell(en1,absId1);

      for(iDigitN = 0; iDigitN < nCells; iDigitN++){
        absId2 = coi->GetCellsAbsId()[iDigitN];

        if(absId2==-1 || absId2==absId1)
          continue;

        en2 = cells->GetCellAmplitude(absId2);
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

      en = cells->GetCellAmplitude(absIdList[iDigit]);
      RecalAmpCell(en,absIdList[iDigit]);

      if(en < localMaxCutE)
        continue; // Maxima only with seed energy at least

      maxEList[iDigitN] = en;
      iDigitN++;
    }
  }

  if(iDigitN == 0){
    AliDebug(1,Form("No local maxima found, assign highest energy cell as maxima, id %d, en cell %2.2f, en cluster %2.2f",idmax,emax,NonLinRecoEnergyScaling(coi, coi->GetNonLinCorrEnergy())));
    iDigitN      = 1;
    maxEList[0]  = emax;
    absIdList[0] = idmax;
  }

  AliDebug(1,Form("In coi E %2.2f (wth non lin. %2.2f), M02 %2.2f, M20 %2.2f, N maxima %d",NonLinRecoEnergyScaling(coi, coi->GetNonLinCorrEnergy()),eCluster,coi->GetM02(),coi->GetM20(),iDigitN));

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

  //________________________________________________________________________________________
  /// Check, in case of clusters with NLM=2, that they are not in the same T-Card.
  //________________________________________________________________________________________
Bool_t  AliAnalysisTaskEMCALPhotonIsolation::IsAbsIDsFromTCard(Int_t absId1, Int_t absId2,
                                                               Int_t & rowDiff, Int_t & colDiff) const
{
  rowDiff = -100;
  colDiff = -100;
  
  if(absId1 == absId2) return kFALSE;
  
    // Check if in same SM, if not for sure not same TCard
  Int_t sm1 = fGeom->GetSuperModuleNumber(absId1);
  Int_t sm2 = fGeom->GetSuperModuleNumber(absId2);
  if ( sm1 != sm2 ) return kFALSE ;
  
    // Get the column and row of each absId
  Int_t iTower = -1, iIphi = -1, iIeta = -1;
  
  Int_t col1, row1;
  fGeom->GetCellIndex(absId1,sm1,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(sm1,iTower,iIphi, iIeta,row1,col1);
  
  Int_t col2, row2;
  fGeom->GetCellIndex(absId2,sm2,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(sm2,iTower,iIphi, iIeta,row2,col2);
  
  Int_t row0 = Int_t(row1-row1%8);
  Int_t col0 = Int_t(col1-col1%2);
  
  Int_t rowDiff0 = row2-row0;
  Int_t colDiff0 = col2-col0;
  
  rowDiff = row1-row2;
  colDiff = col1-col2;
  
    // TCard is made by 2x8 towers
  if ( colDiff0 >=0 && colDiff0 < 2 && rowDiff0 >=0 && rowDiff0 < 8 )
  {
    
      //    printf("\t absId (%d,%d), sm %d; col (%d,%d), colDiff %d; row (%d,%d),rowDiff %d\n",
      //           absId1 , absId2, sm1,
      //           col1, col2, colDiff,
      //           row1, row2, rowDiff);
    return kTRUE ;
  }
  else
    return kFALSE;
}

  //__________________________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::RecalAmpCell(Float_t & amp, Int_t id) const {

    // Recalculate cell energy if recalibration factor

  Int_t iSupMod = -1, iTower = -1, iIphi = -1, iIeta = -1, iphi = -1, ieta = -1;
  fGeom->GetCellIndex(id,iSupMod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);

  amp *= fEMCALRecoUtils->GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
}

  //__________________________________________________________________________________________
Float_t AliAnalysisTaskEMCALPhotonIsolation::RecalEnClust(AliVCluster * coi, AliVCaloCells * cells){

    // Recalibrate the cluster energy, considering the recalibration map and the energy of the cells that compose the cluster
  Int_t   absId = 0;
  Float_t frac  = 0., amp = 0., energy = 0.;

  if(cells){
      // Get the cluster number of cells and list of absId, check what kind of cluster we have
    UShort_t * index    = coi->GetCellsAbsId();
    Double_t * fraction = coi->GetCellsAmplitudeFraction();
    Int_t      ncells   = coi->GetNCells();

      // Loop on the cells, get the cell amplitude and recalibration factor, multiply to get the new energy
    for(Int_t icell = 0; icell < ncells; icell++){
      absId = index[icell];
      frac  = fraction[icell];
      if(frac < 1.e-3)
        frac = 1.; // in case of EMCal, this is set as 0, not used.

      amp = cells->GetCellAmplitude(absId);
      RecalAmpCell(amp, absId);

      energy += amp*frac;
    }

    AliDebug(1,Form("Energy before %f, after %f",NonLinRecoEnergyScaling(coi, coi->GetNonLinCorrEnergy()),energy));

  } // Cells available
  else{
    AliFatal("Cells pointer does not exist!");
  }

  return energy;
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoCellPhiBand(TLorentzVector c, Double_t &etIso, Double_t &phiBandcells){

    // Underlying events study with EMCal cells in phi band

  AliEMCALGeometry* emcalGeom = AliEMCALGeometry::GetInstance();
  Double_t sumEnergyPhiBandCells = 0., sumEnergyConeCells = 0.;

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
  Double_t sumEnergyEtaBandCells = 0., sumEnergyConeCells = 0.;

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
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoClusPhiBand(TLorentzVector c, Double_t m02candidate, Double_t &ptIso, Double_t &phiBandclus, Int_t index){

    // Underlying events study with clusters in phi band

  Double_t sumEnergyPhiBandClus = 0., sumEnergyConeClus = 0., sumpTConeCharged = 0., sumpTPhiBandTracks = 0.;
  Double_t clustTOF = 0., phiClust = 0., etaClust = 0., radius = 0.;
  Double_t phiMin = 0., phiMax = 0., etaMin = 0., etaMax = 0.;

  if(fPeriod != ""){
    etaMin = fGeom->GetArm1EtaMin()+0.03;
    etaMax = fGeom->GetArm1EtaMax()-0.03;
    phiMin = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

    if(fPeriod.Contains("12") || fPeriod.Contains("13"))
      phiMax = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03; // fGeom->GetEMCALPhiMax()-20. = 180. deg (in order not to take the two disabled SM into account in 2012-2013)
    else
      phiMax = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
  }
  else{
    phiMin = (4./9.)*TMath::Pi()+0.03;
    phiMax = TMath::Pi()-0.03;
    etaMin = -0.67;
    etaMax = 0.67;
  }

    // Needs a check on the same cluster
    // AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliClusterContainer *clusters = GetClusterContainer(0);
  Int_t localIndex = 0;
  TLorentzVector nClust(0., 0., 0., 0.);

  for(auto it : clusters->accepted()){ // Check the position of other clusters with respect to the trigger cluster

    AliVCluster* coi = static_cast<AliVCluster*>(it);
    localIndex=coi->GetID();

    if(localIndex==index)
      continue;

    phiClust = etaClust = clustTOF = 0.;
    coi->GetMomentum(nClust, fVertex, AliVCluster::kNonLinCorr);
    if(fIsMC && fNonLinRecoEnergyScaling)
      NonLinRecoEnergyScaling(nClust);

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

    if(nClust.E() < 0.3)
      continue;

    if((phiClust < phiMax) && (phiClust > phiMin) && (etaClust < etaMax) && (etaClust > etaMin)){
      radius = TMath::Sqrt(TMath::Power(phiClust-c.Phi(),2)+TMath::Power(etaClust-c.Eta(),2)); // Define the radius between the leading cluster and the considered cluster
      if(radius > fIsoConeRadius){                                                             // The cluster is outside the isolation cone -> add the cluster pT to pT_UE
	if(TMath::Abs(etaClust - c.Eta()) < fIsoConeRadius)
	  sumEnergyPhiBandClus += nClust.Pt();
      }
      else{                                                                                    // The cluster is inside the isolation cone -> add the cluster pT to pT_iso
	sumEnergyConeClus += nClust.Pt();
	if(fQA && fWho != 2){
	  fTestEtaPhiCone->Fill(c.Eta(),c.Phi());
	  fTestIndex->Fill(index,localIndex);
	  fTestLocalIndexE->Fill(nClust.Pt(),localIndex);
	}
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

  tracksAna->ResetCurrentID();

  AliVTrack        * eTrack      = 0x0;
  AliAODTrack      * aodEtrack   = 0x0;
  AliAODMCParticle * pMC;
  Int_t              iTracksCone = 0  , tracklabel = 0 ;
  Double_t           phiTrack    = 0. , etaTrack   = 0.;

  Float_t frac = 0., ncls = 0., nclsS = 0.;

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
      frac  = 0.;
      ncls  = Float_t(aodEtrack->GetTPCncls ());
      nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls > 0) frac = nclsS / ncls ;

      if(frac > 0.4)
        continue;
    }

    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();

    if((phiTrack < phiMax) && (phiTrack > phiMin) && (etaTrack < etaMax) && (etaTrack > etaMin)){
      radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2)); // Define the radius between the leading cluster and the considered track
      if(radius > fIsoConeRadius){                                                                 // The track is outside the isolation cone -> add the track pT to pT_UE
	if(TMath::Abs(etaTrack - c.Eta()) < fIsoConeRadius)
	  sumpTPhiBandTracks += eTrack->Pt();
      }
      else{                                                                                        // The track is inside the isolation cone -> add the track pT to pT_iso
	sumpTConeCharged += eTrack->Pt();
	if(fIsMC){
	  tracklabel = TMath::Abs(eTrack->GetLabel());
	  pMC        = static_cast<AliAODMCParticle*>(fAODMCParticles->At(tracklabel));
	  if(fWho==1)
	    fTrackResolutionPtMC->Fill(eTrack->Pt(), pMC->Pt() - eTrack->Pt());
	}
	iTracksCone++;
      }
    }
  } // End of tracks loop

  fTrackMultInCone->Fill(iTracksCone);

  if(!fAreasPerEvent)
    fTestEnergyCone->Fill(c.Pt(),sumEnergyConeClus,sumpTConeCharged);

  if(fQA && c.Pt()>14. && c.Pt()<16. && fWho == 2){
    fEtaPhiClusVsEtIsoClus->Fill(c.Eta(),c.Phi(),sumEnergyConeClus);
    fEtaPhiClusVsPtIsoTrack->Fill(c.Eta(),c.Phi(),sumpTConeCharged);
  }

  if(fIsoMethod==1 && fQA && fWho != 2){
    fEtVSM02VSPisotrack->Fill(c.Pt(),m02candidate,sumpTConeCharged);
    fEtVSM02VSEisoclust->Fill(c.Pt(),m02candidate,sumEnergyConeClus);
  }

  if      (fIsoMethod==1) {
    ptIso       = sumEnergyConeClus + sumpTConeCharged;
    phiBandclus = sumEnergyPhiBandClus + sumpTPhiBandTracks;
  }
  else if (fIsoMethod==3) {
    ptIso       = sumEnergyConeClus;
    phiBandclus = sumEnergyPhiBandClus;
  }

  if(!fAreasPerEvent)
    fPtVsConeVsUE->Fill(c.Pt(), ptIso, phiBandclus);

  Double_t stdConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Standard (full) cone area
  Double_t isoConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Cluster (eta, phi)-dependent cone area
  Double_t phiBandArea = 0.; // Cluster phi-dependent eta-band area

  if ( !fTPC4Iso )
    phiBandArea = ((5./9.)*TMath::Pi()-0.06)*2.*fIsoConeRadius-isoConeArea;
  else
    phiBandArea = 2.*TMath::Pi()*2.*fIsoConeRadius-isoConeArea;

  if(fWho == 2 && fAreasPerEvent){
    ComputeConeAreaInEMCal   (c.Eta(), c.Phi()    , isoConeArea);
    ComputePhiBandAreaInEMCal(c.Eta(), isoConeArea, phiBandArea);

    fPtVsNtrConeVsChgCone_Norm->Fill(c.Pt(), sumEnergyConeClus/isoConeArea, sumpTConeCharged/isoConeArea);
    fPhiBandUEClust->Fill(c.Pt(), sumEnergyPhiBandClus/phiBandArea); // Neutral UE energy vs. candidate energy (area normalised)
    fPhiBandUETracks->Fill(c.Pt(), sumpTPhiBandTracks/phiBandArea);  // Charged UE energy vs. candidate energy (area normalised)
    fPtVsConeVsUE_Norm->Fill(c.Pt(), ptIso/isoConeArea, phiBandclus/phiBandArea);

    if(fAnalysispPb){
      fPtvsM02vsSumUE_Norm->Fill(c.Pt(), m02candidate, ((ptIso - phiBandclus * (isoConeArea / phiBandArea))*(stdConeArea / isoConeArea))); // Cone-UE energy vs. shower shape vs. candidate energy (area normalised)
    }
  }
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoClusEtaBand(TLorentzVector c, Double_t m02candidate, Double_t &ptIso, Double_t &etaBandclus, Int_t index){

    // Underlying events study with clusters in eta band

  Float_t sumEnergyEtaBandClus = 0., sumEnergyConeClus = 0., sumpTConeCharged = 0., sumpTEtaBandTracks = 0., sumpTEtaBandTracks_Cside = 0., sumpTEtaBandTracks_Aside = 0.;
  Double_t clustTOF = 0., phiClust = 0., etaClust = 0., radius = 0.;
  Double_t phiMin = 0., phiMax = 0., etaMin = 0., etaMax = 0.;

  if(fPeriod != ""){
    etaMin = fGeom->GetArm1EtaMin()+0.03;
    etaMax = fGeom->GetArm1EtaMax()-0.03;
    phiMin = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

    if(fPeriod.Contains("12") || fPeriod.Contains("13"))
      phiMax = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03; // fGeom->GetEMCALPhiMax()-20. = 180. deg (in order not to take the two disabled SM into account in 2012-2013)
    else
      phiMax = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
  }
  else{
    phiMin = (4./9.)*TMath::Pi()+0.03;
    phiMax = TMath::Pi()-0.03;
    etaMin = -0.67;
    etaMax = 0.67;
  }

    // AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliClusterContainer *clusters = GetClusterContainer(0);
  Int_t localIndex = 0;
  TLorentzVector nClust(0., 0., 0., 0.);

  for(auto it : clusters->accepted()){ // Check the position of other clusters with respect to the trigger cluster

    AliVCluster* coi = static_cast<AliVCluster*>(it);
    localIndex=coi->GetID();

    if(localIndex==index)
      continue;

    phiClust = etaClust = clustTOF = 0.;
    coi->GetMomentum(nClust, fVertex, AliVCluster::kNonLinCorr);
    if(fIsMC && fNonLinRecoEnergyScaling)
      NonLinRecoEnergyScaling(nClust);

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

    if(nClust.E() < 0.3)
      continue;

    if((phiClust < phiMax) && (phiClust > phiMin) && (etaClust < etaMax) && (etaClust > etaMin)){
      radius = TMath::Sqrt(TMath::Power(phiClust-c.Phi(),2)+TMath::Power(etaClust-c.Eta(),2)); // Define the radius between the leading cluster and the considered cluster
      if(radius > fIsoConeRadius){                                                             // The cluster is outside the isolation cone -> add the cluster pT to pT_UE
	if(TMath::Abs(phiClust - c.Phi()) < fIsoConeRadius)
	  sumEnergyEtaBandClus += nClust.Pt();
      }
      // else if(radius<fIsoConeRadius && radius != 0.){
      else{                                                                                         // The cluster is inside the isolation cone -> add the cluster pT to pT_iso
	sumEnergyConeClus += nClust.Pt();
	if(fQA && fWho != 2){
	  fTestEtaPhiCone->Fill(c.Eta(),c.Phi());
	  fTestIndex->Fill(index,localIndex);
	  fTestLocalIndexE->Fill(nClust.Pt(),localIndex);
	}
      }
    }
  } // End of clusters loop

  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }

  if(tracksAna->GetTrackFilterType() != AliEmcalTrackSelection::kHybridTracks)
    AliError(Form("NOT Hybrid Tracks"));

  tracksAna->ResetCurrentID();

  AliVTrack        * eTrack      = 0x0;
  AliAODTrack      * aodEtrack   = 0x0;
  AliAODMCParticle * pMC;
  Int_t              iTracksCone = 0  , tracklabel = 0 ;
  Double_t           phiTrack    = 0. , etaTrack   = 0.;

  Float_t frac = 0., ncls = 0., nclsS = 0.;

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
      frac  = 0.;
      ncls  = Float_t(aodEtrack->GetTPCncls ());
      nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls > 0) frac = nclsS / ncls ;

      if(frac > 0.4)
        continue;
    }

    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();

    if((phiTrack < phiMax) && (phiTrack > phiMin) && (etaTrack < etaMax) && (etaTrack > etaMin)){
      radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2)); // Define the radius between the leading cluster and the considered track
      if(radius > fIsoConeRadius){                                                                 // The track is outside the isolation cone -> add the track pT to pT_UE
	if(TMath::Abs(phiTrack - c.Phi()) < fIsoConeRadius){
	  sumpTEtaBandTracks += eTrack->Pt();

	  if(fWho == 2 && etaTrack < 0.)
	    sumpTEtaBandTracks_Cside += eTrack->Pt();
	  if(fWho == 2 && etaTrack > 0.)
	    sumpTEtaBandTracks_Aside += eTrack->Pt();
	}
      }
      else{                                                                                             // The track is inside the isolation cone -> add the track pT to pT_iso
	sumpTConeCharged += eTrack->Pt();
	if(fIsMC){
	  tracklabel = TMath::Abs(eTrack->GetLabel());
	  pMC        = static_cast<AliAODMCParticle*>(fAODMCParticles->At(tracklabel));
	  if(fWho==1)
	    fTrackResolutionPtMC->Fill(eTrack->Pt(), pMC->Pt() - eTrack->Pt());
	}
	iTracksCone++;
      }
    }
  } // End of tracks loop

  fTrackMultInCone->Fill(iTracksCone);
  
  if(!fAreasPerEvent)
    fTestEnergyCone->Fill(c.Pt(), sumEnergyConeClus, sumpTConeCharged);

  if(fQA && c.Pt()>14. && c.Pt()<16. && fWho == 2){
    fEtaPhiClusVsEtIsoClus->Fill(c.Eta(),c.Phi(),sumEnergyConeClus);
    fEtaPhiClusVsPtIsoTrack->Fill(c.Eta(),c.Phi(),sumpTConeCharged);
  }

  if(fIsoMethod==1 && fQA && fWho != 2){
    fEtVSM02VSPisotrack->Fill(c.Pt(),m02candidate,sumpTConeCharged);
    fEtVSM02VSEisoclust->Fill(c.Pt(),m02candidate,sumEnergyConeClus);
  }

  if      (fIsoMethod==1) {
    ptIso       = sumEnergyConeClus    + sumpTConeCharged;
    etaBandclus = sumEnergyEtaBandClus + sumpTEtaBandTracks;
  }
  else if (fIsoMethod==3) {
    ptIso       = sumEnergyConeClus;
    etaBandclus = sumEnergyEtaBandClus;
  }

  if(!fAreasPerEvent)
    fPtVsConeVsUE->Fill(c.Pt(), ptIso, etaBandclus);

  if(fWho == 2 && !fAreasPerEvent)
    fPtVsConeVsEtaBand->Fill(c.Pt(), ptIso, etaBandclus);

  Double_t stdConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Standard (full) cone area
  Double_t isoConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Cluster (eta, phi)-dependent cone area
  Double_t etaBandArea = 0.; // Cluster phi-dependent eta-band area

  if ( !fTPC4Iso )
    etaBandArea = ((fGeom->GetArm1EtaMax()-0.03)-(fGeom->GetArm1EtaMin()+0.03))*2.*fIsoConeRadius-isoConeArea;
  else
    etaBandArea = 1.74*2.*fIsoConeRadius-isoConeArea;

  if(fWho == 2 && fAreasPerEvent){
    ComputeConeAreaInEMCal   (c.Eta(), c.Phi()    , isoConeArea);
    ComputeEtaBandAreaInEMCal(c.Phi(), isoConeArea, etaBandArea);

    if (fQA) fEtaBandVsConeArea->Fill(isoConeArea, etaBandArea);

    fPtVsNtrConeVsChgCone_Norm->Fill(c.Pt(), sumEnergyConeClus/isoConeArea, sumpTConeCharged/isoConeArea);
    fEtaBandUEClust->Fill(c.Pt(), sumEnergyEtaBandClus/etaBandArea); // Neutral UE energy vs. candidate energy (area normalised)
    fEtaBandUETracks->Fill(c.Pt(), sumpTEtaBandTracks/etaBandArea);  // Charged UE energy vs. candidate energy (area normalised)
    fPtVsConeVsUE_Norm->Fill(c.Pt(), ptIso/isoConeArea, etaBandclus/etaBandArea);

    if(fAnalysispPb){
      fPtvsM02vsSumUE_Norm->Fill(c.Pt(), m02candidate, ((ptIso - etaBandclus * (isoConeArea / etaBandArea))*(stdConeArea / isoConeArea))); // Cone-UE energy vs. shower shape vs. candidate energy (area normalised)
    }
  }
}
  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::EtIsoClusExtraOrthCones(TLorentzVector c, Double_t m02candidate, Double_t &ptIso, Double_t &cones, Int_t index){

    // Underlying events study with clusters in eta band

  Double_t sumEnergyConeClus = 0.     , sumpTConeCharged = 0.     , sumpTPerpConeTrack = 0.;
  Double_t etaCand           = c.Eta(), phiCand          = c.Phi();
  Double_t phiCone1          = phiCand - TMath::PiOver2(); phiCone1 = fmod(phiCone1, TMath::TwoPi());
  Double_t phiCone2          = phiCand + TMath::PiOver2(); phiCone2 = fmod(phiCone2, TMath::TwoPi());

  if(phiCone1 < 0.)
    phiCone1 += 2*TMath::Pi();

  Double_t clustTOF = 0., phiClust = 0., etaClust = 0., radius = 0.;
  Double_t phiMin = 0., phiMax = 0., etaMin = 0., etaMax = 0.;

  if(fPeriod != ""){
    etaMin = fGeom->GetArm1EtaMin()+0.03;
    etaMax = fGeom->GetArm1EtaMax()-0.03;
    phiMin = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

    if(fPeriod.Contains("12") || fPeriod.Contains("13"))
      phiMax = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03; // fGeom->GetEMCALPhiMax()-20. = 180. deg (in order not to take the two disabled SM into account in 2012-2013)
    else
      phiMax = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
  }
  else{
    phiMin = (4./9.)*TMath::Pi()+0.03;
    phiMax = TMath::Pi()-0.03;
    etaMin = -0.67;
    etaMax = 0.67;
  }

    // AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliClusterContainer *clusters = GetClusterContainer(0);
  Int_t localIndex = 0;
  TLorentzVector nClust(0., 0., 0., 0.);

  for(auto it : clusters->accepted()){ // Check the position of other clusters with respect to the trigger cluster

    AliVCluster* coi = static_cast<AliVCluster*>(it);
    localIndex = coi->GetID();

    if(localIndex == index)
      continue;

    phiClust = etaClust = clustTOF = 0.;
    coi->GetMomentum(nClust, fVertex, AliVCluster::kNonLinCorr);
    if(fIsMC && fNonLinRecoEnergyScaling)
      NonLinRecoEnergyScaling(nClust);

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
      if(clustTOF < -30. || clustTOF > 30.)
        continue;
    }

    if(fTMClusterInConeRejected){
      if(ClustTrackMatching(coi, kFALSE))
        continue;
    }

    if(nClust.E() < 0.3)
      continue;

    if((phiClust < phiMax) && (phiClust > phiMin) && (etaClust < etaMax) && (etaClust > etaMin)){
      radius = TMath::Sqrt(TMath::Power(phiClust-phiCand,2)+TMath::Power(etaClust-etaCand,2)); // Define the radius between the leading cluster and the considered cluster
      if(radius <= fIsoConeRadius){                                                            // The cluster is inside the isolation cone -> add the cluster pT to pT_iso
	sumEnergyConeClus += nClust.Pt();
	if(fQA && fWho != 2){
	  fTestEtaPhiCone->Fill(etaCand,phiCand);
	  fTestIndex->Fill(index,localIndex);
	  fTestLocalIndexE->Fill(nClust.Pt(),localIndex);
	}
      }
    }
  } // End of clusters loop

  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }

  if(tracksAna->GetTrackFilterType() != AliEmcalTrackSelection::kHybridTracks)
    AliError(Form("NOT Hybrid Tracks"));

  tracksAna->ResetCurrentID();

  AliVTrack        * eTrack      = 0x0;
  AliAODTrack      * aodEtrack   = 0x0;
  AliAODMCParticle * pMC;
  Int_t              iTracksCone = 0  , tracklabel = 0 ;
  Double_t           phiTrack    = 0. , etaTrack   = 0., distToCone1 = 0., distToCone2 = 0.;

  Float_t frac = 0., ncls = 0., nclsS = 0.;

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

    if((eTrack->Pt()) < 0.2)
      continue;

    if((eTrack->GetStatus() & AliAODTrack::kITSrefit) != AliAODTrack::kITSrefit)
      continue;

    if(!fIsEsd){
      frac  = 0.;
      ncls  = Float_t(aodEtrack->GetTPCncls ());
      nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls > 0) frac = nclsS / ncls ;

      if(frac > 0.4)
        continue;
    }

    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();
    radius   = TMath::Sqrt(TMath::Power(phiTrack - phiCand, 2) + TMath::Power(etaTrack - etaCand, 2));

    if(radius <= fIsoConeRadius){ // The track is inside the isolation cone -> add the track pT to pT_iso
      sumpTConeCharged += eTrack->Pt();
      if(fIsMC){
        tracklabel = TMath::Abs(eTrack->GetLabel());
        pMC        = static_cast<AliAODMCParticle*>(fAODMCParticles->At(tracklabel));
        if(fWho==1)
          fTrackResolutionPtMC->Fill(eTrack->Pt(), pMC->Pt() - eTrack->Pt());
      }
      iTracksCone ++;
    }
    else{
      distToCone1 = TMath::Sqrt(TMath::Power(etaTrack - etaCand, 2) + TMath::Power(phiTrack - phiCone1, 2));
      distToCone2 = TMath::Sqrt(TMath::Power(etaTrack - etaCand, 2) + TMath::Power(phiTrack - phiCone2, 2));

      if((distToCone1 < fIsoConeRadius) || (distToCone2 < fIsoConeRadius)) // The track is inside one of the two orthogonal cones -> add the track pT to pT_UE
        sumpTPerpConeTrack += eTrack->Pt();
    }
  }

  fTrackMultInCone->Fill(iTracksCone);

  if(!fAreasPerEvent)
    fTestEnergyCone->Fill(c.Pt(), sumEnergyConeClus, sumpTPerpConeTrack);

  ptIso = sumEnergyConeClus + sumpTConeCharged;
  cones = fExtraPerpConesFactor*sumpTPerpConeTrack; // Scaling charged-only UE to neutral + charged UE

  if(!fAreasPerEvent)
    fPtVsConeVsUE->Fill(c.Pt(), ptIso, cones);

  Double_t stdConeArea   = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Standard (full) cone area
  Double_t isoConeArea   = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Cluster (eta, phi)-dependent cone area
  Double_t perpConesArea = 2.*isoConeArea; // Cluster (eta, phi)-dependent perpendicular cones area

  if(fWho == 2 && fAreasPerEvent){
    ComputeConeAreaInEMCal(c.Eta(), c.Phi(), isoConeArea);
    ComputeConeAreaInTPC  (c.Eta(), perpConesArea);
    perpConesArea         *= 2.;

    fPtVsNtrConeVsChgCone_Norm->Fill(c.Pt(), sumEnergyConeClus/isoConeArea, sumpTPerpConeTrack/isoConeArea);
    fPerpConesUETracks->Fill(c.Pt(), sumpTPerpConeTrack/perpConesArea);       // Charged UE
    fPtVsConeVsUE_Norm->Fill(c.Pt(), ptIso/isoConeArea, cones/perpConesArea); // Total extrapolated UE

    if(fAnalysispPb){

      fPtvsM02vsSumUE_Norm->Fill(c.Pt(), m02candidate, ((ptIso - cones * (isoConeArea / perpConesArea))*(stdConeArea / isoConeArea))); // Cone-UE energy vs. shower shape vs. candidate energy (area normalised)
    }
  }
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackPhiBand(TLorentzVector c, Double_t m02candidate, Double_t &ptIso, Double_t &phiBandtrack){

    // Underlying events study with tracks in phi band

  Double_t sumpTConeCharged = 0., sumpTPhiBandTrack = 0.;
  Double_t phiMin = 0., phiMax = 2.*TMath::Pi(), etaMin = -0.87, etaMax = 0.87;

  if(!fTPC4Iso){
    if(fPeriod != ""){
      etaMin = fGeom->GetArm1EtaMin()+0.03;
      etaMax = fGeom->GetArm1EtaMax()-0.03;
      phiMin = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

      if(fPeriod.Contains("12") || fPeriod.Contains("13"))
        phiMax = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03; // fGeom->GetEMCALPhiMax()-20. = 180. deg (in order not to take the two disabled SM into account in 2012-2013)
      else
        phiMax = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
    }
    else{
      etaMin = -0.67;
      etaMax = 0.67;
      phiMin = (4./9.)*TMath::Pi()+0.03;
      phiMax = TMath::Pi()-0.03;
    }
  }

  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }

  tracksAna->ResetCurrentID();

  AliVTrack        * eTrack      = 0x0;
  AliAODTrack      * aodEtrack   = 0x0;
  AliAODMCParticle * pMC;
  Int_t              iTracksCone = 0  , tracklabel = 0 ;
  Double_t           phiTrack    = 0. , etaTrack   = 0., radius = 0.;

  Float_t frac = 0., ncls = 0., nclsS = 0.;

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
      frac  = 0.;
      ncls  = Float_t(aodEtrack->GetTPCncls ());
      nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls > 0) frac = nclsS / ncls ;

      if(frac > 0.4) continue;
    }

    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();

    if((phiTrack < phiMax) && (phiTrack > phiMin) && (etaTrack < etaMax) && (etaTrack > etaMin)){
      radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2)); // Define the radius between the leading cluster and the considered track
      if(radius > fIsoConeRadius){                                                                 // The track is outside the isolation cone -> add the track pT to pT_UE
        if(TMath::Abs(etaTrack - c.Eta()) < fIsoConeRadius)
          sumpTPhiBandTrack += eTrack->Pt();
      }
      else{                                                                                             // The track is inside the isolation cone -> add the track pT to pT_iso
        sumpTConeCharged += eTrack->Pt();
        if(fIsMC){
          tracklabel = TMath::Abs(eTrack->GetLabel());
          pMC        = static_cast<AliAODMCParticle*>(fAODMCParticles->At(tracklabel));
          if(fWho==1)
            fTrackResolutionPtMC->Fill(eTrack->Pt(), pMC->Pt() - eTrack->Pt());
        }
        iTracksCone++;
      }
    }
  }

  fTrackMultInCone->Fill(iTracksCone);

  ptIso        = sumpTConeCharged;
  phiBandtrack = sumpTPhiBandTrack;

  Double_t stdConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Standard (full) cone area
  Double_t isoConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.);
  Double_t phiBandArea = 0.;

  if ( !fTPC4Iso )
    phiBandArea = ((5./9.)*TMath::Pi()-0.06)*2.*fIsoConeRadius-isoConeArea;
  else
    phiBandArea = 2.*TMath::Pi()*2.*fIsoConeRadius-isoConeArea;

  if(fWho == 2 && fAreasPerEvent){
    ComputeConeAreaInTPC   (c.Eta(), isoConeArea);
    ComputePhiBandAreaInTPC(c.Eta(), isoConeArea, phiBandArea);

    fPhiBandUETracks->Fill(c.Pt(), phiBandtrack/phiBandArea);

    if(fAnalysispPb){
      fPtvsM02vsSumUE_Norm->Fill(c.Pt(), m02candidate, ((ptIso - phiBandtrack * (isoConeArea / phiBandArea))*(stdConeArea / isoConeArea))); // Cone-UE energy vs. shower shape vs. candidate energy (area normalised)
    }
  }
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackEtaBand(TLorentzVector c, Double_t m02candidate, Double_t &ptIso, Double_t &etaBandtrack){

    // Underlying events study with tracks in eta band

  Double_t sumpTConeCharged = 0., sumpTEtaBandTrack = 0.;
  Double_t phiMin = 0., phiMax = 2.*TMath::Pi(), etaMin = -0.87, etaMax = 0.87;

  if(!fTPC4Iso){
    if(fPeriod != ""){
      etaMin = fGeom->GetArm1EtaMin()+0.03;
      etaMax = fGeom->GetArm1EtaMax()-0.03;
      phiMin = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

      if(fPeriod.Contains("12") || fPeriod.Contains("13"))
        phiMax = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03; // fGeom->GetEMCALPhiMax()-20. = 180. deg (in order not to take the two disabled SM into account in 2012-2013)
      else
        phiMax = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
    }
    else{
      etaMin = -0.67;
      etaMax = 0.67;
      phiMin = (4./9.)*TMath::Pi()+0.03;
      phiMax = TMath::Pi()-0.03;
    }
  }

  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }

  tracksAna->ResetCurrentID();

  AliVTrack        * eTrack      = 0x0;
  AliAODTrack      * aodEtrack   = 0x0;
  AliAODMCParticle * pMC;
  Int_t              iTracksCone = 0  , tracklabel = 0 ;
  Double_t           phiTrack    = 0. , etaTrack   = 0., radius = 0.;

  Float_t frac = 0., ncls = 0., nclsS = 0.;

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
      frac  = 0.;
      ncls  = Float_t(aodEtrack->GetTPCncls ());
      nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls > 0) frac = nclsS / ncls ;

      if(frac > 0.4)
        continue;
    }

    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();

    if((phiTrack < phiMax) && (phiTrack > phiMin) && (etaTrack < etaMax) && (etaTrack > etaMin)){
      radius = TMath::Sqrt(TMath::Power(phiTrack - c.Phi(),2)+TMath::Power(etaTrack - c.Eta(),2)); // Define the radius between the leading cluster and the considered track
      if(radius > fIsoConeRadius){                                                                 // The track is outside the isolation cone -> add the track pT to pT_UE
        if(TMath::Abs(phiTrack - c.Phi()) < fIsoConeRadius)
          sumpTEtaBandTrack += eTrack->Pt();
      }
      else{                                                                                             // The track is inside the isolation cone -> add the track pT to pT_iso
        sumpTConeCharged += eTrack->Pt();
        if(fIsMC){
          tracklabel = TMath::Abs(eTrack->GetLabel());
          pMC        = static_cast<AliAODMCParticle*>(fAODMCParticles->At(tracklabel));
          if(fWho==1)
            fTrackResolutionPtMC->Fill(eTrack->Pt(), pMC->Pt() - eTrack->Pt());
        }
        iTracksCone++;
      }
    }
  }

  fTrackMultInCone->Fill(iTracksCone);

  ptIso        = sumpTConeCharged;
  etaBandtrack = sumpTEtaBandTrack;

  Double_t stdConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Standard (full) cone area
  Double_t isoConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Cluster (eta, phi)-dependent cone area
  Double_t etaBandArea = 0.; // Cluster phi-dependent eta-band area

  if ( !fTPC4Iso )
    etaBandArea = ((fGeom->GetArm1EtaMax()-0.03)-(fGeom->GetArm1EtaMin()+0.03))*2.*fIsoConeRadius-isoConeArea;
  else
    etaBandArea = 1.74*2.*fIsoConeRadius-isoConeArea;

  if(fWho == 2 && fAreasPerEvent){
    ComputeConeAreaInTPC   (c.Eta()    , isoConeArea);
    ComputeEtaBandAreaInTPC(isoConeArea, etaBandArea);

    fEtaBandUETracks->Fill(c.Pt(), etaBandtrack/etaBandArea);

    if(fAnalysispPb){
      fPtvsM02vsSumUE_Norm->Fill(c.Pt(), m02candidate, ((ptIso - etaBandtrack * (isoConeArea / etaBandArea))*(stdConeArea / isoConeArea))); // Cone-UE energy vs. shower shape vs. candidate energy (area normalised)
    }
  }
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackOrthCones(TLorentzVector c, Double_t m02candidate, Double_t &ptIso, Double_t &cones){

    // Underlying events study with tracks in orthogonal cones in TPC

  Double_t sumpTConeCharged = 0., sumpTPerpConeTrack = 0.;
  Double_t etaClus = c.Eta();
  Double_t phiClus = c.Phi();
  Double_t phiCone1 = phiClus - TMath::PiOver2();
  phiCone1=fmod(phiCone1,TMath::TwoPi());
  Double_t phiCone2 = phiClus + TMath::PiOver2();
  phiCone2=fmod(phiCone2,TMath::TwoPi());

  if(phiCone1 < 0.)
    phiCone1 += 2*TMath::Pi();

  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }

  tracksAna->ResetCurrentID();

  AliVTrack        * eTrack      = 0x0;
  AliAODTrack      * aodEtrack   = 0x0;
  AliAODMCParticle * pMC;
  Int_t              iTracksCone = 0  , tracklabel = 0 ;
  Double_t           phiTrack    = 0. , etaTrack   = 0., distToClust = 0., distToCone1 = 0., distToCone2 = 0.;

  Float_t frac = 0., ncls = 0., nclsS = 0.;

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
      frac  = 0.;
      ncls  = Float_t(aodEtrack->GetTPCncls ());
      nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls > 0) frac = nclsS / ncls ;

      if(frac > 0.4) continue;
    }

    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();
    distToClust = TMath::Sqrt(TMath::Power(etaTrack-etaClus, 2)+TMath::Power(phiTrack-phiClus, 2));

    if(distToClust<fIsoConeRadius){ // The track is inside the isolation cone -> add the track pT to pT_iso
      sumpTConeCharged += eTrack->Pt();
      if(fIsMC){
        tracklabel = TMath::Abs(eTrack->GetLabel());
        pMC        = static_cast<AliAODMCParticle*>(fAODMCParticles->At(tracklabel));
        if(fWho==1)
          fTrackResolutionPtMC->Fill(eTrack->Pt(), pMC->Pt() - eTrack->Pt());
      }
      iTracksCone++;
    }

    else{
        // Distances from the centres of the two Orthogonal Cones
      distToCone1 = TMath::Sqrt(TMath::Power(etaTrack-etaClus, 2)+TMath::Power(phiTrack-phiCone1, 2));
      distToCone2 = TMath::Sqrt(TMath::Power(etaTrack-etaClus, 2)+TMath::Power(phiTrack-phiCone2, 2));

        // The track is inside one of the two orthogonal cones -> add the track pT to pT_UE
      if((distToCone1 < fIsoConeRadius) || (distToCone2 < fIsoConeRadius))
        sumpTPerpConeTrack += eTrack->Pt();
    }
  }

  fTrackMultInCone->Fill(iTracksCone);

  ptIso = sumpTConeCharged;
  cones = sumpTPerpConeTrack;

  Double_t stdConeArea   = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Standard (full) cone area
  Double_t isoConeArea   = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Cluster (eta, phi)-dependent cone area
  Double_t perpConesArea = 2.*isoConeArea; // Cluster (eta, phi)-dependent perpendicular cones area

  if(fWho == 2 && fAreasPerEvent){
    ComputeConeAreaInTPC(c.Eta(), isoConeArea);
    perpConesArea       = 2.*isoConeArea;

    fPerpConesUETracks->Fill(c.Pt(), cones/perpConesArea);
    if(fAnalysispPb) fPtvsM02vsSumUE_Norm->Fill(c.Pt(), m02candidate, ((ptIso - cones * (isoConeArea / perpConesArea))*(stdConeArea / isoConeArea))); // Cone-UE energy vs. shower shape vs. candidate energy (area normalised)
  }
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::PtIsoTrackFullTPC(TLorentzVector c, Double_t &ptIso, Double_t &full){

    // Underlying events study with tracks in full TPC except a back to back band

  Double_t sumpTConeCharged = 0., sumpTTPCexceptB2B = 0.;

  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  if(!tracksAna){
    AliError(Form("Could not retrieve tracks !"));
    return;
  }

  tracksAna->ResetCurrentID();

  AliVTrack        * eTrack      = 0x0;
  AliAODTrack      * aodEtrack   = 0x0;
  AliAODMCParticle * pMC;
  Int_t              iTracksCone = 0  , tracklabel = 0 ;
  Double_t           phiTrack    = 0. , etaTrack   = 0., radius = 0., dphiUp = 0., dphiDown = 0.;

  Float_t frac = 0., ncls = 0., nclsS = 0.;

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
      frac  = 0.;
      ncls  = Float_t(aodEtrack->GetTPCncls ());
      nclsS = Float_t(aodEtrack->GetTPCnclsS());
      if(ncls> 0)  frac =  nclsS / ncls ;

      if(frac > 0.4) continue;
    }

    phiTrack = eTrack->Phi();
    etaTrack = eTrack->Eta();

    radius = TMath::Sqrt(TMath::Power(phiTrack-c.Phi(),2)+TMath::Power(etaTrack-c.Eta(),2)); // Define the radius between the leading cluster and the considered track

    if(radius > fIsoConeRadius){                                                             // The track is outside the isolation cone -> add the track pT to pT_UE
      dphiUp = c.Phi() + TMath::Pi() - fIsoConeRadius;
      dphiDown = c.Phi() + TMath::Pi() + fIsoConeRadius;

      if(phiTrack < dphiDown && phiTrack> dphiUp)
        sumpTTPCexceptB2B += eTrack->Pt();
    }
    else{                                                                                         // The track is inside the isolation cone -> add the track pT to pT_iso
      sumpTConeCharged += eTrack->Pt();
      if(fIsMC){
	tracklabel = TMath::Abs(eTrack->GetLabel());
        pMC        = static_cast<AliAODMCParticle*>(fAODMCParticles->At(tracklabel));
        if(fWho==1)
          fTrackResolutionPtMC->Fill(eTrack->Pt(), pMC->Pt() - eTrack->Pt());
      }
      iTracksCone++;
    }
  }

  fTrackMultInCone->Fill(iTracksCone);

  ptIso = sumpTConeCharged;
  full = sumpTTPCexceptB2B;
}

  //__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::CheckBoundaries(TLorentzVector vecCOI){

    // Check if the cone around the considered cluster is in EMCal acceptance

  Double_t etaClust = 0., phiClust = 0.;
  Double_t etaMax = 0., etaMinDCal_InnerEdge = 0., phiMinEMCal = 0., phiMaxEMCal = 0., phiMinDCal = 0. , phiMaxDCal_FullSM = 0., phiMaxDCal = 0.;
  Bool_t   isINBoundaries = kFALSE;

  etaClust = vecCOI.Eta();
  phiClust = vecCOI.Phi();

  if(fTPC4Iso){
    etaMax = 0.87-fFiducialCut;

    if(!fPeriod.IsNull()){
      phiMinEMCal = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

      if(fPeriod.Contains("12") || fPeriod.Contains("13"))
        phiMaxEMCal = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03;
      else
        phiMaxEMCal = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;

      if(fPeriod.Contains("15") || fPeriod.Contains("16") || fPeriod.Contains("17")){
	phiMinDCal        = (fGeom->GetDCALPhiMin())*TMath::DegToRad()+0.03;
	phiMaxDCal_FullSM = (fGeom->GetEMCGeometry()->GetDCALStandardPhiMax())*TMath::DegToRad()-0.03;
	phiMaxDCal        = (fGeom->GetDCALPhiMax())*TMath::DegToRad()-0.03;
      }

      if(fPeriod.Contains("10") || fPeriod.Contains("11") || fPeriod.Contains("12") || fPeriod.Contains("13")){
	if(TMath::Abs(etaClust) > etaMax || phiClust < phiMinEMCal || phiClust > phiMaxEMCal)
	  isINBoundaries = kFALSE;
	else
	  isINBoundaries = kTRUE;
      } // Run 1
      else if(fPeriod.Contains("15") || fPeriod.Contains("16") || fPeriod.Contains("17")){
	if(TMath::Abs(etaClust) > etaMax || phiClust < phiMinEMCal || phiClust > phiMaxDCal || (phiClust > phiMaxEMCal && phiClust < phiMinDCal))
	  isINBoundaries = kFALSE;
	else
	  isINBoundaries = kTRUE;
      } // Run 2
    }
    else{ // If no period set, default for 2011-2013 (2*5 EMCal SM)
      phiMinEMCal = (4./9.)*TMath::Pi()+0.03;
      phiMaxEMCal = TMath::Pi()-0.03;

      if(TMath::Abs(etaClust) > etaMax || phiClust < phiMinEMCal || phiClust > phiMaxEMCal)
	isINBoundaries = kFALSE;
      else
	isINBoundaries = kTRUE;
    }
  }
  else{
    if(!fPeriod.IsNull()){
      etaMax      = fGeom->GetArm1EtaMax()-0.03-fFiducialCut;
      phiMinEMCal = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03+fFiducialCut;

      if(fPeriod.Contains("12") || fPeriod.Contains("13"))
        phiMaxEMCal = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03-fFiducialCut;
      else
        phiMaxEMCal = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03-fFiducialCut;

      if(fPeriod.Contains("15") || fPeriod.Contains("16") || fPeriod.Contains("17")){
	etaMinDCal_InnerEdge = fGeom->GetEMCGeometry()->GetDCALInnerExtandedEta()+0.03+fFiducialCut;
	phiMinDCal           = (fGeom->GetDCALPhiMin())*TMath::DegToRad()+0.03+fFiducialCut;
	phiMaxDCal_FullSM    = (fGeom->GetEMCGeometry()->GetDCALStandardPhiMax())*TMath::DegToRad()-0.03-fFiducialCut;
	phiMaxDCal           = (fGeom->GetDCALPhiMax())*TMath::DegToRad()-0.03-fFiducialCut;
      }

      if(fPeriod.Contains("10") || fPeriod.Contains("11") || fPeriod.Contains("12") || fPeriod.Contains("13")){
	if(TMath::Abs(etaClust) > etaMax || phiClust < phiMinEMCal || phiClust > phiMaxEMCal)
	  isINBoundaries = kFALSE;
	else
	  isINBoundaries = kTRUE;
      } // Run 1
      else if(fPeriod.Contains("15") || fPeriod.Contains("16") || fPeriod.Contains("17")){
	if(TMath::Abs(etaClust) > etaMax || phiClust < phiMinEMCal || phiClust > phiMaxDCal || (phiClust > phiMaxEMCal && phiClust < phiMinDCal) || (phiClust >= phiMinDCal && phiClust < phiMaxDCal_FullSM && TMath::Abs(etaClust) < etaMinDCal_InnerEdge))
	  isINBoundaries = kFALSE;
	else
	  isINBoundaries = kTRUE;
      } // Run 2
    }
    else{ // If no period set, default for 2011-2013 (2*5 EMCal SM)
      etaMax      = 0.67-fFiducialCut;
      phiMinEMCal = (4./9.)*TMath::Pi()+0.03+fFiducialCut;
      phiMaxEMCal = TMath::Pi()-0.03-fFiducialCut;

      if(TMath::Abs(etaClust) > etaMax || phiClust < phiMinEMCal || phiClust > phiMaxEMCal)
	isINBoundaries = kFALSE;
      else
	isINBoundaries = kTRUE;
    }
  }

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
  Int_t clustPDG, p2ccharge;
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

          Int_t idxdaug1 = momP2Check->GetDaughterFirst();
          Int_t idxdaug2 = momP2Check->GetDaughterLast();

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
	  if(fWho != 2 && fQA){
	    fpi0VSclusterVSIsolation->Fill(momP2Check->E()*TMath::Sin(momP2Check->Theta()), energyCLS, isolation);
	    fpi0VSclusterVSM02->Fill(momP2Check->E()*TMath::Sin(momP2Check->Theta()), energyCLS,ss);
	    fpi0VSM02VSIsolation->Fill(momP2Check->E()*TMath::Sin(momP2Check->Theta()), ss, isolation);
	  }
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

      Int_t firstidx=momP2Check->GetDaughterFirst();
      Int_t lastidx=momP2Check->GetDaughterLast();

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
          Int_t idxaunt1 = grandma->GetDaughterFirst();
          Int_t idxaunt2 = grandma->GetDaughterLast();

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
	  if(fQA){
	    fpi0VSclusterVSIsolation->Fill(grandma->E()*TMath::Sin(grandma->Theta()), energyCLS, isolation);
	    fpi0VSclusterVSM02->Fill(grandma->E()*TMath::Sin(momP2Check->Theta()), energyCLS, ss);
	    fpi0VSM02VSIsolation->Fill(grandma->E()*TMath::Sin(grandma->Theta()), ss, isolation);
	  }
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

  Double_t clustTOF=0;
//  Double_t invMass=0;

//  Double_t invMasspi0Min=0.135-0.035;
//  Double_t invMasspi0Max=0.135+0.035;

//  Double_t invMassetaMin=0.548-0.035;
//  Double_t invMassetaMax=0.548+0.035;

    // AliParticleContainer *clusters = static_cast<AliParticleContainer*>(fParticleCollArray.At(1));
  AliClusterContainer *clusters   = GetClusterContainer(0);
  TLorentzVector       nClust     (0., 0., 0., 0.);
  Int_t                localIndex = 0;

  for(auto it :clusters->accepted()){ // Check the position of other clusters with respect to the leading cluster
    AliVCluster* coi = static_cast<AliVCluster*>(it);

    localIndex=coi->GetID();
    if(localIndex==index)
      continue;
    else{
      localIndex++;

      // TLorentzVector nClust(0., 0., 0., 0.);
      coi->GetMomentum(nClust, fVertex, AliVCluster::kNonLinCorr);
      if(fIsMC && fNonLinRecoEnergyScaling)
	NonLinRecoEnergyScaling(nClust);

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

    }
  } // End of clusters loop
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::IsolationAndUEinEMCAL(AliVCluster *coi, Double_t& isolation, Double_t& ue, Double_t eTThreshold, Int_t index){

  Double_t isoConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.);
  Double_t etaBandArea = ((fGeom->GetArm1EtaMax()-0.03)-(fGeom->GetArm1EtaMin()+0.03))*2.*fIsoConeRadius-isoConeArea;
  Double_t phiBandArea = 0.;

  if(fPeriod != ""){
    if(fPeriod.Contains("12") || fPeriod.Contains("13"))
      phiBandArea = (((fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03)-((fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03))*2.*fIsoConeRadius-isoConeArea;
    else
      phiBandArea = (((fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03)-((fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03))*2.*fIsoConeRadius-isoConeArea;
  }
  else{
    phiBandArea = (5./9.)*TMath::Pi()*2.*fIsoConeRadius-isoConeArea;
  }

  TLorentzVector vecCOI;
  coi->GetMomentum(vecCOI, fVertex, AliVCluster::kNonLinCorr);
  if(fIsMC && fNonLinRecoEnergyScaling)
    NonLinRecoEnergyScaling(vecCOI);

  Double_t eTCOI = vecCOI.Et();

  // Set smearing for MC
  Double_t m02COI = 0.;
  if(fSSsmearing)
    ApplySmearing(coi, m02COI);
  else
    m02COI = coi->GetM02();

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
        EtIsoClusPhiBand(vecCOI, m02COI, isolation, ue, index);

        ue = ue * (isoConeArea / phiBandArea); // Normalisation of UE wrt UE area
        break;

      case 1: // Eta band
        EtIsoClusEtaBand(vecCOI, m02COI, isolation, ue, index);

        ue = ue * (isoConeArea / etaBandArea); // Normalisation of UE wrt UE area
        break;

      case 4: // Extrapolated cones
	EtIsoClusExtraOrthCones(vecCOI, m02COI, isolation, ue, index);

        ue = ue * (isoConeArea / etaBandArea); // Normalisation of UE wrt UE area
        break;
    }

      if(fWho==2){
	if(!fAnalysispPb)
	  fPtvsM02vsSum->Fill(vecCOI.Pt(),m02COI,isolation);

	if(fAnalysispPb)
	  isolation=isolation-ue; // UE subtraction

	if(!fAreasPerEvent)
	  fPtvsM02vsSumUE->Fill(vecCOI.Pt(),m02COI,isolation);
        // fEtIsoClust->Fill(vecCOI.Pt(),isolation);
      }

      if(isolation<eTThreshold){
        FillInvMassHistograms(kTRUE, m02COI, vecCOI, index, isolation);

        // if(fWho==2){
        //   fPtvsM02iso->Fill(vecCOI.Pt(),m02COI);
	//   fPtIsolatedNClust->Fill(vecCOI.Pt());
        // }

        fPtisoT=vecCOI.Pt();
        fM02isoT=m02COI;

        if(fM02mincut < m02COI && m02COI < fM02maxcut){
          // if(fWho==2)
          //   fEtIsolatedClust->Fill(eTCOI);

          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
      }
      else{
        if(isolation>3.)
          FillInvMassHistograms(kFALSE, m02COI, vecCOI, index, isolation);

        // if(fWho==2)
        //   fPtvsM02noiso->Fill(vecCOI.Pt(),m02COI);

        fPtnoisoT=vecCOI.Pt();
        fM02noisoT=m02COI;
      }
      break;

    case 2: // EMCal TRACKS
      switch(fUEMethod)
    {
      case 0: // Phi band
        PtIsoTrackPhiBand(vecCOI, m02COI, isolation, ue);

        ue = ue * (isoConeArea / phiBandArea); // Normalisation of UE wrt UE area

      case 1: // Eta band
        PtIsoTrackEtaBand(vecCOI, m02COI, isolation, ue);

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
	if(!fAnalysispPb)
	  fPtvsM02vsSum->Fill(vecCOI.Pt(),m02COI,isolation);

	if(fAnalysispPb)
	  isolation=isolation-ue; // UE subtraction

	if(!fAreasPerEvent)
	  fPtvsM02vsSumUE->Fill(vecCOI.Pt(),m02COI,isolation);
        // fPtIsoTrack->Fill(vecCOI.Pt() , isolation);
      }

      if(isolation<eTThreshold){
        FillInvMassHistograms(kTRUE, m02COI, vecCOI, index, isolation);

        // if(fWho==2){
        //   fPtvsM02iso->Fill(vecCOI.Pt(),m02COI);
	//   fPtIsolatedNTracks->Fill(vecCOI.Pt());
        // }

        fPtisoT=vecCOI.Pt();
        fM02isoT=m02COI;

        if(fM02mincut < m02COI && m02COI < fM02maxcut){
          // if(fWho==2)
          //   fEtIsolatedTracks->Fill(eTCOI);

          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
      }
      else{
        if(isolation>3.)
          FillInvMassHistograms(kFALSE, m02COI, vecCOI, index, isolation);

        // if(fWho==2)
        //   fPtvsM02noiso->Fill(vecCOI.Pt(),m02COI);

        fPtnoisoT=vecCOI.Pt();
        fM02noisoT=m02COI;
      }
      break;
  }
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::IsolationAndUEinTPC(AliVCluster *coi, Double_t& isolation,Double_t& ue,Double_t eTThreshold, Int_t index){

    // EMCal + TPC (only tracks for isolation since IsoCone goes out of EMCal)

  Double_t isoConeArea   = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.);
  Double_t etaBandAreaTr = 1.74*2.*fIsoConeRadius-isoConeArea; // 1.74 = 2*0.87
  Double_t phiBandAreaTr = 2.*TMath::Pi()*2.*fIsoConeRadius-isoConeArea;
  Double_t perpConesArea = 2.*isoConeArea;
  Double_t fullTPCArea   = 1.74*2.*TMath::Pi()-1.74*2.*fIsoConeRadius-isoConeArea;

  TLorentzVector vecCOI;
  coi->GetMomentum(vecCOI, fVertex, AliVCluster::kNonLinCorr);
  if(fIsMC && fNonLinRecoEnergyScaling)
    NonLinRecoEnergyScaling(vecCOI);

  Double_t eTCOI = vecCOI.Et();

  // Set smearing for MC
  Double_t m02COI = 0.;
  if(fSSsmearing)
    ApplySmearing(coi, m02COI);
  else
    m02COI = coi->GetM02();

  switch(fUEMethod)
  {
    case 0: // Phi band
      PtIsoTrackPhiBand(vecCOI, m02COI, isolation, ue);

      ue = ue * (isoConeArea / phiBandAreaTr); // Normalisation of UE wrt UE area

      if(fWho==2)
	if(!fAnalysispPb)
	  fPtvsM02vsSum->Fill(vecCOI.Pt(),m02COI,isolation);

      if(fAnalysispPb)
        isolation=isolation-ue; // UE subtraction

      if(fWho==2 && !fAreasPerEvent){
	fPtvsM02vsSumUE->Fill(vecCOI.Pt(),m02COI,isolation);
        // fPtIsoTrack->Fill(vecCOI.Pt(), isolation);
      }

      if(isolation<eTThreshold){
        FillInvMassHistograms(kTRUE, m02COI, vecCOI, index,isolation);

        // if(fWho==2){
        //   fPtvsM02iso->Fill(vecCOI.Pt(),m02COI);
	//   fPtIsolatedNTracks->Fill(vecCOI.Pt());
        // }

        fPtisoT=vecCOI.Pt();
        fM02isoT=m02COI;

        if(fM02mincut < m02COI && m02COI < fM02maxcut){
          // if(fWho==2)
          //   fEtIsolatedTracks->Fill(eTCOI);

          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
      }
      else{
        if(isolation>3.)
          FillInvMassHistograms(kFALSE, m02COI, vecCOI, index,isolation);

        // if(fWho==2)
        //   fPtvsM02noiso->Fill(vecCOI.Pt(),m02COI);

        fPtnoisoT=vecCOI.Pt();
        fM02noisoT=m02COI;
      }
      break;

    case 1: // Eta band
      PtIsoTrackEtaBand(vecCOI, m02COI, isolation, ue);

      ue = ue * (isoConeArea / etaBandAreaTr); // Normalisation of UE wrt UE area

      if(fWho==2)
	if(!fAnalysispPb)
	  fPtvsM02vsSum->Fill(vecCOI.Pt(),m02COI,isolation);

      if(fAnalysispPb)
        isolation=isolation-ue; // UE subtraction

      if(fWho==2 && !fAreasPerEvent){
	fPtvsM02vsSumUE->Fill(vecCOI.Pt(),m02COI,isolation);
        // fPtIsoTrack->Fill(vecCOI.Pt(), isolation);
      }

      if(isolation<eTThreshold){
        FillInvMassHistograms(kTRUE, m02COI, vecCOI, index,isolation);

        // if(fWho==2){
        //   fPtvsM02iso->Fill(vecCOI.Pt(),m02COI);
	//   fPtIsolatedNTracks->Fill(vecCOI.Pt());
        // }

        fPtisoT=vecCOI.Pt();
        fM02isoT=m02COI;

        if(fM02mincut < m02COI && m02COI < fM02maxcut){
          // if(fWho==2)
          //   fEtIsolatedTracks->Fill(eTCOI);

          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
      }
      else{
        if(isolation>3.)
          FillInvMassHistograms(kFALSE, m02COI, vecCOI, index,isolation);

        // if(fWho==2)
        //   fPtvsM02noiso->Fill(vecCOI.Pt(),m02COI);

        fPtnoisoT=vecCOI.Pt();
        fM02noisoT=m02COI;
      }
      break;

    case 2: // Cones
      PtIsoTrackOrthCones(vecCOI, m02COI, isolation, ue);

      if(fWho==2)
        fPerpConesUETracks->Fill(vecCOI.Pt() , ue);

      ue = ue * (isoConeArea / perpConesArea); // Normalisation of UE wrt UE area

      if(fWho==2)
	if(!fAnalysispPb)
	  fPtvsM02vsSum->Fill(vecCOI.Pt(),m02COI,isolation);

      if(fAnalysispPb)
        isolation=isolation-ue; // UE subtraction

      if(fWho==2 && !fAreasPerEvent){
	fPtvsM02vsSumUE->Fill(vecCOI.Pt(),m02COI,isolation);
        // fPtIsoTrack->Fill(vecCOI.Pt(), isolation);
      }

      if(isolation<eTThreshold){
        FillInvMassHistograms(kTRUE, m02COI, vecCOI, index, isolation);

        // if(fWho==2){
        //   fPtvsM02iso->Fill(vecCOI.Pt(),m02COI);
	//   fPtIsolatedNTracks->Fill(vecCOI.Pt());
        // }

        fPtisoT=vecCOI.Pt();
        fM02isoT=m02COI;

        if(fM02mincut < m02COI && m02COI < fM02maxcut){
          // if(fWho==2)
          //   fEtIsolatedTracks->Fill(eTCOI);

          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
      }
      else{
        if(isolation>3.)
          FillInvMassHistograms(kFALSE, m02COI, vecCOI, index, isolation);

        // if(fWho==2)
        //   fPtvsM02noiso->Fill(vecCOI.Pt(),m02COI);

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
	if(!fAnalysispPb)
	  fPtvsM02vsSum->Fill(vecCOI.Pt(),m02COI,isolation);

      if(fAnalysispPb)
        isolation=isolation-ue; // UE subtraction

      if(fWho==2 && !fAreasPerEvent){
	fPtvsM02vsSumUE->Fill(vecCOI.Pt(),m02COI,isolation);
        // fPtIsoTrack->Fill(vecCOI.Pt(), isolation);
      }
        // fTracksConeEtaPt->Fill(isolation, vecCOI.Eta(), vecCOI.Pt());
        // fTracksConeEtaM02->Fill(isolation, vecCOI.Eta(), m02COI);

      if(isolation<eTThreshold){
        FillInvMassHistograms(kTRUE, m02COI, vecCOI, index, isolation);

        // if(fWho==2){
        //   fPtvsM02iso->Fill(vecCOI.Pt(),m02COI);
	//   fPtIsolatedNTracks->Fill(vecCOI.Pt());
        // }

        fPtisoT=vecCOI.Pt();
        fM02isoT=m02COI;

        if(fM02mincut < m02COI && m02COI < fM02maxcut){
          // if(fWho==2)
          //   fEtIsolatedTracks->Fill(eTCOI);

          fEtisolatedT=eTCOI;
          fPtisolatedT=vecCOI.Pt();
        }
      }
      else{
        if(isolation>3.)
          FillInvMassHistograms(kFALSE, m02COI, vecCOI, index, isolation);

        // if(fWho==2)
        //   fPtvsM02noiso->Fill(vecCOI.Pt(),m02COI);

        fPtnoisoT=vecCOI.Pt();
        fM02noisoT=m02COI;
      }
      break;
  }
}

  //__________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPhotonIsolation::FillGeneralHistograms(AliVCluster *coi, TLorentzVector vecCOI, Int_t index){

    // Fill the histograms for underlying event and isolation studies

    // I would like to remove this part and fill the tracks multiplicity histogram in FillQAHistograms, is that ok for thnSparses? (especially cause here the histogram is filled several times per event)

  AliTrackContainer *tracksAna = GetTrackContainer("filterTracksAna");
  tracksAna->ResetCurrentID();

  const Int_t nTracks = tracksAna->GetNAcceptedTracks();

  fTrackMult->Fill(nTracks);

    // Definition of the Array for Davide's Output
  const Int_t ndims =   fNDimensions;
  Double_t outputValues[ndims];

  Double_t eTCOI = 0.;
  eTCOI = vecCOI.Et();

  // Set smearing for MC
  Double_t m02COI = 0.;
  if(fSSsmearing)
    ApplySmearing(coi, m02COI);
  else
    m02COI = coi->GetM02();

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
      break;
  }

  Double_t isolation=0, ue=0;

  if(!fTPC4Iso)
    IsolationAndUEinEMCAL(coi,isolation,ue,eTThreshold,index);
  else
    IsolationAndUEinTPC(coi,isolation,ue,eTThreshold,index);

  if(fIsMC && fWho != 2)
    LookforParticle(coi->GetLabel(),vecCOI.Et(),vecCOI.Phi(),vecCOI.Eta(),coi->GetTOF()*1e9,m02COI,isolation);

  switch(fWho)
  {
    case 0:
      flambda0T=m02COI;   // for all neutral clusters
      fEtT=vecCOI.Et();   // for all neutral clusters
      fPtT=vecCOI.Pt();   // for all neutral clusters
      fetaT=vecCOI.Eta(); // for all neutral clusters
      fphiT=vecCOI.Phi(); // for all neutral clusters
      fsumEtisoconeT=isolation;
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

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::ComputeConeAreaInEMCal(Double_t etaCand, Double_t phiCand, Double_t &coneArea){

    // Compute the isolation cone area for neutral + charged isolation depending on the cluster position (for fiducial cuts lower than cone radius)

  Double_t phiMin = 0., phiMax = 0., etaMin  = 0., etaMax  = 0.;
  Double_t d_eta  = 0., d_phi  = 0., sqd_eta = 0., sqd_phi = 0., sqRadius = 0., fullConeArea = 0.;

  coneArea     = 0.;
  sqRadius     = TMath::Power(fIsoConeRadius, 2.);
  fullConeArea = TMath::Pi()*sqRadius;

  if(fPeriod != ""){
    etaMin = fGeom->GetArm1EtaMin()+0.03;
    etaMax = fGeom->GetArm1EtaMax()-0.03;
    phiMin = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

    if(fPeriod.Contains("12") || fPeriod.Contains("13"))
      phiMax = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03;
    else
      phiMax = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
  }
  else{
    phiMin = (4./9.)*TMath::Pi()+0.03;
    phiMax = TMath::Pi()-0.03;
    etaMin = -0.67;
    etaMax = 0.67;
  }

  if((etaCand > etaMax-fIsoConeRadius) && (phiCand > phiMax-fIsoConeRadius)){ // Cluster in EMCal top right corner
    d_eta = TMath::Abs(etaMax-etaCand); sqd_eta = TMath::Power(d_eta, 2.);
    d_phi = TMath::Abs(phiMax-phiCand); sqd_phi = TMath::Power(d_phi, 2.);

    if(TMath::Sqrt(sqd_eta+sqd_phi) >= fIsoConeRadius)
      coneArea = fullConeArea - (sqRadius*TMath::ACos(d_phi/fIsoConeRadius)-d_phi*TMath::Sqrt(sqRadius-sqd_phi)) - (sqRadius*TMath::ACos(d_eta/fIsoConeRadius)-d_eta*TMath::Sqrt(sqRadius-sqd_eta));
    else
      coneArea = 0.25*fullConeArea + d_eta*d_phi + d_eta*TMath::Sqrt(sqRadius-sqd_eta) + d_phi*TMath::Sqrt(sqRadius-sqd_phi) + sqRadius*TMath::ACos(0.5*TMath::Sqrt(4*sqRadius-sqd_eta-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta), 2.))/fIsoConeRadius) - 0.25*TMath::Sqrt(4*sqRadius-sqd_eta-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta), 2.))*TMath::Sqrt(sqd_eta+TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta), 2.)) + 0.5*(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta))*d_eta + sqRadius*TMath::ACos(0.5*TMath::Sqrt(4*sqRadius-sqd_phi-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi), 2.))/fIsoConeRadius) - 0.25*TMath::Sqrt(4*sqRadius-sqd_phi-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi), 2.))*TMath::Sqrt(sqd_phi+TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi), 2.)) + 0.5*(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi))*d_phi;
  }
  else if((etaCand > etaMax-fIsoConeRadius) && (phiCand < phiMin+fIsoConeRadius)){ // Cluster in EMCal bottom right corner
    d_eta = TMath::Abs(etaMax-etaCand); sqd_eta = TMath::Power(d_eta, 2.);
    d_phi = TMath::Abs(phiMin-phiCand); sqd_phi = TMath::Power(d_phi, 2.);

    if(TMath::Sqrt(sqd_eta+sqd_phi) >= fIsoConeRadius)
      coneArea = fullConeArea - (sqRadius*TMath::ACos(d_phi/fIsoConeRadius)-d_phi*TMath::Sqrt(sqRadius-sqd_phi)) - (sqRadius*TMath::ACos(d_eta/fIsoConeRadius)-d_eta*TMath::Sqrt(sqRadius-sqd_eta));
    else
      coneArea = 0.25*fullConeArea + d_eta*d_phi + d_eta*TMath::Sqrt(sqRadius-sqd_eta) + d_phi*TMath::Sqrt(sqRadius-sqd_phi) + sqRadius*TMath::ACos(0.5*TMath::Sqrt(4*sqRadius-sqd_eta-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta), 2.))/fIsoConeRadius) - 0.25*TMath::Sqrt(4*sqRadius-sqd_eta-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta), 2.))*TMath::Sqrt(sqd_eta+TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta), 2.)) + 0.5*(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta))*d_eta + sqRadius*TMath::ACos(0.5*TMath::Sqrt(4*sqRadius-sqd_phi-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi), 2.))/fIsoConeRadius) - 0.25*TMath::Sqrt(4*sqRadius-sqd_phi-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi), 2.))*TMath::Sqrt(sqd_phi+TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi), 2.)) + 0.5*(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi))*d_phi;
  }
  else if((etaCand < etaMin+fIsoConeRadius) && (phiCand < phiMin+fIsoConeRadius)){ // Cluster in EMCal bottom left corner
    d_eta = TMath::Abs(etaMin-etaCand); sqd_eta = TMath::Power(d_eta, 2.);
    d_phi = TMath::Abs(phiMin-phiCand); sqd_phi = TMath::Power(d_phi, 2.);

    if(TMath::Sqrt(sqd_eta+sqd_phi) >= fIsoConeRadius)
      coneArea = fullConeArea - (sqRadius*TMath::ACos(d_phi/fIsoConeRadius)-d_phi*TMath::Sqrt(sqRadius-sqd_phi)) - (sqRadius*TMath::ACos(d_eta/fIsoConeRadius)-d_eta*TMath::Sqrt(sqRadius-sqd_eta));
    else
      coneArea = 0.25*fullConeArea + d_eta*d_phi + d_eta*TMath::Sqrt(sqRadius-sqd_eta) + d_phi*TMath::Sqrt(sqRadius-sqd_phi) + sqRadius*TMath::ACos(0.5*TMath::Sqrt(4*sqRadius-sqd_eta-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta), 2.))/fIsoConeRadius) - 0.25*TMath::Sqrt(4*sqRadius-sqd_eta-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta), 2.))*TMath::Sqrt(sqd_eta+TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta), 2.)) + 0.5*(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta))*d_eta + sqRadius*TMath::ACos(0.5*TMath::Sqrt(4*sqRadius-sqd_phi-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi), 2.))/fIsoConeRadius) - 0.25*TMath::Sqrt(4*sqRadius-sqd_phi-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi), 2.))*TMath::Sqrt(sqd_phi+TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi), 2.)) + 0.5*(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi))*d_phi;
  }
  else if((etaCand < etaMin+fIsoConeRadius) && (phiCand > phiMax-fIsoConeRadius)){ // Cluster in EMCal top left corner
    d_eta = TMath::Abs(etaMin-etaCand); sqd_eta = TMath::Power(d_eta, 2.);
    d_phi = TMath::Abs(phiMax-phiCand); sqd_phi = TMath::Power(d_phi, 2.);

    if(TMath::Sqrt(sqd_eta+sqd_phi) >= fIsoConeRadius)
      coneArea = fullConeArea - (sqRadius*TMath::ACos(d_phi/fIsoConeRadius)-d_phi*TMath::Sqrt(sqRadius-sqd_phi)) - (sqRadius*TMath::ACos(d_eta/fIsoConeRadius)-d_eta*TMath::Sqrt(sqRadius-sqd_eta));
    else
      coneArea = 0.25*fullConeArea + d_eta*d_phi + d_eta*TMath::Sqrt(sqRadius-sqd_eta) + d_phi*TMath::Sqrt(sqRadius-sqd_phi) + sqRadius*TMath::ACos(0.5*TMath::Sqrt(4*sqRadius-sqd_eta-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta), 2.))/fIsoConeRadius) - 0.25*TMath::Sqrt(4*sqRadius-sqd_eta-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta), 2.))*TMath::Sqrt(sqd_eta+TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta), 2.)) + 0.5*(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_eta))*d_eta + sqRadius*TMath::ACos(0.5*TMath::Sqrt(4*sqRadius-sqd_phi-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi), 2.))/fIsoConeRadius) - 0.25*TMath::Sqrt(4*sqRadius-sqd_phi-TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi), 2.))*TMath::Sqrt(sqd_phi+TMath::Power(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi), 2.)) + 0.5*(fIsoConeRadius-TMath::Sqrt(sqRadius-sqd_phi))*d_phi;
  }
  else if((etaCand > etaMax-fIsoConeRadius) && (phiCand > phiMin+fIsoConeRadius && phiCand < phiMax-fIsoConeRadius)){ // Cluster on EMCal right border except corners
    d_eta = TMath::Abs(etaMax-etaCand); sqd_eta = TMath::Power(d_eta, 2.);

    coneArea = fullConeArea - (sqRadius*TMath::ACos(d_eta/fIsoConeRadius)-d_eta*TMath::Sqrt(sqRadius-sqd_eta));
  }
  else if((etaCand < etaMin+fIsoConeRadius) && (phiCand > phiMin+fIsoConeRadius && phiCand < phiMax-fIsoConeRadius)){ // Cluster on EMCal left border except corners
    d_eta = TMath::Abs(etaMin-etaCand); sqd_eta = TMath::Power(d_eta, 2.);

    coneArea = fullConeArea - (sqRadius*TMath::ACos(d_eta/fIsoConeRadius)-d_eta*TMath::Sqrt(sqRadius-sqd_eta));
  }
  else if((phiCand > phiMax-fIsoConeRadius) && (etaCand > etaMin+fIsoConeRadius && etaCand < etaMax-fIsoConeRadius)){ // Cluster on EMCal top border except corners
    d_phi = TMath::Abs(phiMax-phiCand); sqd_phi = TMath::Power(d_phi, 2.);

    coneArea = fullConeArea - (sqRadius*TMath::ACos(d_phi/fIsoConeRadius)-d_phi*TMath::Sqrt(sqRadius-sqd_phi));
  }
  else if((phiCand < phiMin+fIsoConeRadius) && (etaCand > etaMin+fIsoConeRadius && etaCand < etaMax-fIsoConeRadius)){ // Cluster on EMCal bottom border except corners
    d_phi = TMath::Abs(phiMin-phiCand); sqd_phi = TMath::Power(d_phi, 2.);

    coneArea = fullConeArea - (sqRadius*TMath::ACos(d_phi/fIsoConeRadius)-d_phi*TMath::Sqrt(sqRadius-sqd_phi));
  }
  else // Full cone area (EMCal centre)
    coneArea = fullConeArea;
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::ComputeConeAreaInTPC(Double_t etaCand, Double_t &coneArea){

    // Compute the isolation cone area for charged-only isolation depending on the cluster position (for fiducial cuts lower than cone radius)

  Double_t etaMin = -0.87, etaMax = 0.87, d_eta = 0., sqd_eta = 0., sqRadius = 0., fullConeArea = 0.;

  coneArea     = 0.;
  sqRadius     = TMath::Power(fIsoConeRadius, 2.);
  fullConeArea = TMath::Pi()*sqRadius;

  if(etaCand > etaMax-fIsoConeRadius){ // Cluster on EMCal right border, cone going outside TPC
    d_eta = TMath::Abs(etaMax-etaCand); sqd_eta = TMath::Power(d_eta, 2.);

    coneArea = fullConeArea - (sqRadius*TMath::ACos(d_eta/fIsoConeRadius)-d_eta*TMath::Sqrt(sqRadius-sqd_eta));
  }
  else if(etaCand < etaMin+fIsoConeRadius){ // Cluster on EMCal left border, cone going outside TPC
    d_eta = TMath::Abs(etaMin-etaCand); sqd_eta = TMath::Power(d_eta, 2.);

    coneArea = fullConeArea - (sqRadius*TMath::ACos(d_eta/fIsoConeRadius)-d_eta*TMath::Sqrt(sqRadius-sqd_eta));
  }
  else // Full cone area (EMCal centre, cone not going outside TPC)
    coneArea = fullConeArea;
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::ComputeEtaBandAreaInEMCal(Double_t phiCand, Double_t coneArea, Double_t &etaBandArea){

    // Compute the eta-band area for neutral + charged isolation depending on the cluster position (for fiducial cuts lower than cone radius)

  Double_t phiMin = 0., phiMax = 0., etaMin = 0., etaMax = 0.;
  Double_t d_phi  = 0.;

  etaBandArea = 0.;

  if(fPeriod != ""){
    etaMin = fGeom->GetArm1EtaMin()+0.03;
    etaMax = fGeom->GetArm1EtaMax()-0.03;
    phiMin = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

    if(fPeriod.Contains("12") || fPeriod.Contains("13"))
      phiMax = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03;
    else
      phiMax = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
  }
  else{
    phiMin = (4./9.)*TMath::Pi()+0.03;
    phiMax = TMath::Pi()-0.03;
    etaMin = -0.67;
    etaMax = 0.67;
  }

  if(phiCand > phiMax-fIsoConeRadius){ // Cluster on EMCal top border
    d_phi = TMath::Abs(phiMax-phiCand);

    etaBandArea = (etaMax-etaMin)*(fIsoConeRadius+d_phi);
  }
  else if(phiCand < phiMin+fIsoConeRadius){ // Cluster on EMCal bottom border
    d_phi = TMath::Abs(phiMin-phiCand);

    etaBandArea = (etaMax-etaMin)*(fIsoConeRadius+d_phi);
  }
  else // Full band area (EMCal centre)
    etaBandArea = (etaMax-etaMin)*2.*fIsoConeRadius;

  // Whatever the case, remove the cone area (computed in ComputeConeAreaInEMCal(), called before ComputeEtaBandAreaInEMCal())
  etaBandArea -= coneArea;
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::ComputeEtaBandAreaInTPC(Double_t coneArea, Double_t &etaBandArea){

    // Compute the eta-band area for charged-only isolation depending on the cluster associated cone area (for fiducial cuts lower than cone radius)

  Double_t etaMin = -0.87, etaMax = 0.87;

  etaBandArea = 0.;
  
  etaBandArea = (etaMax-etaMin)*2.*fIsoConeRadius - coneArea;
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::ComputePhiBandAreaInEMCal(Double_t etaCand, Double_t coneArea, Double_t &phiBandArea){

    // Compute the phi-band area for neutral + charged isolation depending on the cluster position (for fiducial cuts lower than cone radius)

  Double_t phiMin = 0., phiMax = 0., etaMin = 0., etaMax = 0.;
  Double_t d_eta  = 0.;

  phiBandArea = 0.;

  if(fPeriod != ""){
    etaMin = fGeom->GetArm1EtaMin()+0.03;
    etaMax = fGeom->GetArm1EtaMax()-0.03;
    phiMin = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

    if(fPeriod.Contains("12") || fPeriod.Contains("13"))
      phiMax = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03;
    else
      phiMax = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
  }
  else{
    phiMin = (4./9.)*TMath::Pi()+0.03;
    phiMax = TMath::Pi()-0.03;
    etaMin = -0.67;
    etaMax = 0.67;
  }

  if(etaCand > etaMax-fIsoConeRadius){ // Cluster on EMCal right border
    d_eta = TMath::Abs(etaMax-etaCand);

    phiBandArea = (phiMax-phiMin)*(fIsoConeRadius+d_eta);
  }
  else if(etaCand < etaMin+fIsoConeRadius){ // Cluster on EMCal left border
    d_eta = TMath::Abs(etaMin-etaCand);

    phiBandArea = (phiMax-phiMin)*(fIsoConeRadius+d_eta);
  }
  else // Full band area (EMCal centre)
    phiBandArea = (phiMax-phiMin)*2.*fIsoConeRadius;

  // Whatever the case, remove the cone area (computed in ComputeConeAreaInEMCal(), called before ComputePhiBandAreaInEMCal())
  phiBandArea -= coneArea;
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::ComputePhiBandAreaInTPC(Double_t etaCand, Double_t coneArea, Double_t &phiBandArea){

    // Compute the phi-band area for charged-only isolation depending on the cluster position (for fiducial cuts lower than cone radius)

  Double_t phiMin = 0., phiMax = 2.*TMath::Pi(), etaMin = -0.87, etaMax = 0.87, d_eta  = 0.;

  phiBandArea = 0.;

  if(etaCand > etaMax-fIsoConeRadius){ // Cluster on EMCal right border, cone going outside TPC
    d_eta = TMath::Abs(etaMax-etaCand);

    phiBandArea = (phiMax-phiMin)*(fIsoConeRadius+d_eta);
  }
  else if(etaCand < etaMin+fIsoConeRadius){ // Cluster on EMCal left border, cone going outside TPC
    d_eta = TMath::Abs(etaMin-etaCand);

    phiBandArea = (phiMax-phiMin)*(fIsoConeRadius+d_eta);
  }
  else // Full band area (EMCal centre, cone not going outside TPC)
    phiBandArea = (phiMax-phiMin)*2.*fIsoConeRadius;

  // Whatever the case, remove the cone area (computed in ComputeConeAreaInTPC(), called before ComputePhiBandAreaInTPC())
  phiBandArea -= coneArea;
}

  //__________________________________________________________________________
void AliAnalysisTaskEMCALPhotonIsolation::ApplySmearing(AliVCluster *coi, Double_t &m02COI){

  // Compute and apply smearing to shower shape in MC

  Int_t          nlm        = 0;
  AliVCaloCells *fCaloCells = InputEvent()->GetEMCALCells();
  if(fCaloCells)
    nlm = GetNLM(coi, fCaloCells);

  if(coi->GetM02() > 0.1){
    if(nlm == 1){
      if((fSSsmearwidth != 0.)){
	TRandom3 *ran = new TRandom3(0);

	if(fWhich == 0){
	  Float_t smear = ran->Landau(fSSsmear_mean, fSSsmearwidth);
	  if(fSSsmear_mean == 0 || (fSSsmear_mean != 0 && coi->GetID()%3 == 0))
	    m02COI = coi->GetM02() + smear;
	} // Landau Smearing
	else{
	  Float_t smear = ran->Gaus(fSSsmear_mean, fSSsmearwidth);
	  if(fSSsmear_mean == 0 || (fSSsmear_mean != 0 && coi->GetID()%3 == 0))
	    m02COI = coi->GetM02() + smear;
	} // Gaussian Smearing
      }
      else{
	m02COI = coi->GetM02();
      }
    }
    else
      m02COI = coi->GetM02();
  }
}

  //_________________________________________________________________________

void AliAnalysisTaskEMCALPhotonIsolation::AddParticleToUEMC(Double_t& sumUE, AliAODMCParticle* mcpp, Double_t eta, Double_t phi){

  Double_t etaMax = 0./*, etaMinDCal_InnerEdge = 0.*/, phiMinEMCal = 0., phiMaxEMCal = 0./*, phiMinDCal = 0. , phiMaxDCal_FullSM = 0., phiMaxDCal = 0.*/;

  Double_t etap = mcpp->Eta();
  Double_t phip = mcpp->Phi();

  if(!fTPC4Iso){
    etaMax = fGeom->GetArm1EtaMax()-0.03;
    if(!fPeriod.IsNull()){
      phiMinEMCal = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;        //  80     deg (4/9*pi  rad) plus  0.03 rad = 1.426 rad
      if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	phiMaxEMCal = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03; // 180     deg (pi      rad) minus 0.03 rad = 3.112 rad
      else
	phiMaxEMCal = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;     // 188.137 deg (3.284   rad) minus 0.03 rad = 3.254 rad
    }
    else{ // If no period set, default for 2011-2013 (2*5 EMCal SM)
      phiMinEMCal = (4./9.)*TMath::Pi()+0.03;
      phiMaxEMCal = TMath::Pi()-0.03;
    }

    if((TMath::Abs(etap) <= etaMax) && (phip >= phiMinEMCal) && (phip <= phiMaxEMCal)){
      switch(fUEMethod)
      {
        case 0:{ // Phi band

          if(TMath::Abs(etap-eta) < fIsoConeRadius){
	    if(mcpp->Charge() != 0)
	      sumUE += mcpp->Pt();
	    else if(mcpp->GetPdgCode() == 22)
	      sumUE += mcpp->E()*(TMath::Sin(mcpp->Theta()));
	  }

          break;
        }

        case 1:{ // Eta band

          if(TMath::Abs(phip-phi) < fIsoConeRadius){
	    if(mcpp->Charge() != 0){
	      sumUE += mcpp->Pt();
	      if(fWho == 2){
	        fEtaBandUENeutral_MC->Fill(sumUE);
	      }
	    }
	    else if(mcpp->GetPdgCode() == 22){
	      sumUE += mcpp->E()*(TMath::Sin(mcpp->Theta()));
	      if(fWho == 2){
	        fEtaBandUECharged_MC->Fill(sumUE);
	      }
	    }
	  }

          break;
        }

        case 4:{ // Extrapolated cones
          double etacone1 = eta;
          double etacone2 = eta;
          double phicone1 = phi - TMath::PiOver2();
          double phicone2 = phi + TMath::PiOver2();

          if(phicone1 < 0.)
            phicone1 += 2*TMath::Pi();

          if(TMath::Sqrt(TMath::Power(etap-etacone1,2)+TMath::Power(phip-phicone1,2)) < fIsoConeRadius ||
             TMath::Sqrt(TMath::Power(etap-etacone2,2)+TMath::Power(phip-phicone2,2)) < fIsoConeRadius){
            if(mcpp->Charge() != 0)
	      sumUE += mcpp->Pt();
	  }
	  
	  sumUE *= fExtraPerpConesFactor; // Neutral + charged extrapolation

          break;
        }
      }
    }
  }
  else{ // TPC for isolation/UE : charged particles only (implies Pt() instead of E_T)
    etaMax = 0.87;

    if(TMath::Abs(etap) < etaMax){
      switch(fUEMethod)
      {
        case 0:{ // Phi band
	  
          if(TMath::Abs(etap-eta) < fIsoConeRadius)
	    sumUE += mcpp->Pt();
	  
          break;
        }

        case 1:{ // Eta band

          if(TMath::Abs(phip-phi) < fIsoConeRadius)
	    sumUE += mcpp->Pt();

          break;
        }

        case 2:{ // Orthogonal Cones
          double etacone1 = eta;
          double etacone2 = eta;
          double phicone1 = phi - TMath::PiOver2();
          double phicone2 = phi + TMath::PiOver2();

          if(phicone1 < 0.)
            phicone1 += 2*TMath::Pi();

          if(TMath::Sqrt(TMath::Power(etap-etacone1,2)+TMath::Power(phip-phicone1,2)) < fIsoConeRadius ||
             TMath::Sqrt(TMath::Power(etap-etacone2,2)+TMath::Power(phip-phicone2,2)) < fIsoConeRadius)
            sumUE += mcpp->Pt();

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

void AliAnalysisTaskEMCALPhotonIsolation::CalculateUEDensityMC(Double_t etaCand, Double_t phiCand, Double_t &sumUE){

  Double_t isoConeArea   = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.);
  Double_t etaBandArea   = ((fGeom->GetArm1EtaMax()-0.03)-(fGeom->GetArm1EtaMin()+0.03))*2.*fIsoConeRadius-isoConeArea;
  Double_t phiBandArea   = ((5./9.)*TMath::Pi()-0.06)*2.*fIsoConeRadius-isoConeArea;
  Double_t etaBandAreaTr = 1.74*2.*fIsoConeRadius-isoConeArea;
  Double_t phiBandAreaTr = 2.*TMath::Pi()*2.*fIsoConeRadius-isoConeArea;
  Double_t perpConesArea = 2.*isoConeArea;
  Double_t fullTPCArea   = 1.74*2.*TMath::Pi()-1.74*2.*fIsoConeRadius-isoConeArea;

  if(!fTPC4Iso){
    switch(fUEMethod)
    {
      case 0:{
	if(fWho == 2 && fAreasPerEvent){
	  ComputeConeAreaInEMCal   (etaCand, phiCand    , isoConeArea);
	  ComputePhiBandAreaInEMCal(etaCand, isoConeArea, phiBandArea);
	}

        sumUE = sumUE * (isoConeArea / phiBandArea);
        break;
      }

      case 1:{

	if(fWho == 2 && fAreasPerEvent){
	  ComputeConeAreaInEMCal   (etaCand, phiCand    , isoConeArea);
	  ComputeEtaBandAreaInEMCal(phiCand, isoConeArea, etaBandArea);
	}

	sumUE = sumUE * (isoConeArea / etaBandArea);
	break;
      }

      case 4:{
	if(fWho == 2 && fAreasPerEvent){
	  ComputeConeAreaInEMCal(etaCand, phiCand, isoConeArea);
	  ComputeConeAreaInTPC  (etaCand, perpConesArea);
	  perpConesArea         *= 2.;
	}

        sumUE = sumUE * (isoConeArea / perpConesArea);
        break;
      }
    }
  }
  else{
    switch(fUEMethod)
    {
      case 0:{
	if(fWho == 2 && fAreasPerEvent){
	  ComputeConeAreaInTPC   (etaCand, isoConeArea);
	  ComputePhiBandAreaInTPC(etaCand, isoConeArea, phiBandArea);
	}

	sumUE = sumUE * (isoConeArea / phiBandAreaTr);
	break;
      }

      case 1:{
	if(fWho == 2 && fAreasPerEvent){
	  ComputeConeAreaInTPC   (etaCand    , isoConeArea);
	  ComputeEtaBandAreaInTPC(isoConeArea, etaBandArea);
	}

        sumUE = sumUE * (isoConeArea / etaBandAreaTr);
        break;
      }

      case 2:{
	if(fWho == 2 && fAreasPerEvent){
	  ComputeConeAreaInTPC(etaCand, isoConeArea);
	  perpConesArea = 2.*isoConeArea;
	}

        sumUE = sumUE * (isoConeArea / perpConesArea);
        break;
      }

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

  if(!fStack && !fAODMCParticles){
    cout<<"No stack saved\n"; return;
  }

  Double_t stdConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Standard (full) cone area
  Double_t isoConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.);

  Double_t eT = 0., sumEiso = 0., sumUE = 0., phi = 0., eta = 0., radius = 0., phip = 0., etap = 0.;
  Double_t etaMax_fidu = 0./*, etaMinDCal_InnerEdge = 0.*/, phiMinEMCal_fidu = 0., phiMaxEMCal_fidu = 0./*, phiMinDCal = 0. , phiMaxDCal_FullSM = 0., phiMaxDCal = 0.*/;

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

  AliAODMCParticle *mcpart, *mom, *mcpp;
  Int_t pdg = 0, photonlabel = 0, momidx = 0, mcppmomidx = 0, mompdg = 0;

  for(int iTr = 0; iTr < nTracks; iTr ++){
    eT=0.; phi=0.; eta=0.;

    mcpart = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTr));

    if(fWho==2){
    fGenPromptPhotonSel->Fill(0.5);
    }

    if(mcpart->GetStatus()>10)
      continue;

    if(fWho==2){
    fGenPromptPhotonSel->Fill(1.5);
    }

    if(!mcpart->IsPrimary())
      continue;

    if(!mcpart->IsPhysicalPrimary())
      continue;

    if(fWho==2){
    fGenPromptPhotonSel->Fill(2.5);
    }

    pdg = mcpart->GetPdgCode();
    if(pdg != 22 /*|| mcpart->GetLabel()!=8*/)
      continue;

    if(fWho==2){
    fGenPromptPhotonSel->Fill(3.5);
    }

    momidx = mcpart->GetMother();
    if(momidx > 0 && momidx < nTracks){
      mom    = static_cast<AliAODMCParticle*>(fAODMCParticles->At(momidx));
      mompdg = TMath::Abs(mom->GetPdgCode());
    }
    else
      mompdg = pdg;

    if(fWho==2){
    if(mompdg != 22) continue;              // Discard particles whose mother is not a photon

    fGenPromptPhotonSel->Fill(4.5);
    }

    // if(fPythiaHeader->ProcessType() != 201 && fPythiaHeader->ProcessType() != 202) continue; // Discard particles which do not come from prompt photon processes
    // OR
    if(fmcHeader->GetEventType() != 14 && fmcHeader->GetEventType() != 29 && fWho == 2) continue; // Discard particles which do not come from prompt photon processes

    if(fWho==2){
    fGenPromptPhotonSel->Fill(5.5);
    }

    eta = mcpart->Eta();
    phi = mcpart->Phi();

    if(!fTPC4Iso){ // Check photons in EMCal
      etaMax_fidu = fGeom->GetArm1EtaMax()-0.03-fFiducialCut;

      if(!fPeriod.IsNull()){
	phiMinEMCal_fidu = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03+fFiducialCut;
	if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	  phiMaxEMCal_fidu = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03-fFiducialCut;
	else
	  phiMaxEMCal_fidu = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03-fFiducialCut;
      }
      else{ // If no period set, default for 2011-2013 (2*5 EMCal SM)
	phiMinEMCal_fidu = (4./9.)*TMath::Pi()+0.03+fFiducialCut;
	phiMaxEMCal_fidu = TMath::Pi()-0.03-fFiducialCut;
      }

      if((TMath::Abs(eta) > etaMax_fidu) || (phi < phiMinEMCal_fidu || phi > phiMaxEMCal_fidu))
	continue;
    }
    else{ // Check photons in TPC
      etaMax_fidu = 0.87-fFiducialCut;

      if(!fPeriod.IsNull()){
	phiMinEMCal_fidu = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;
	if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	  phiMaxEMCal_fidu = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03;
	else
	  phiMaxEMCal_fidu = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
      }
      else{ // If no period set, default for 2011-2013 (2*5 EMCal SM)
	phiMinEMCal_fidu = (4./9.)*TMath::Pi()+0.03;
	phiMaxEMCal_fidu = TMath::Pi()-0.03;
      }

      if((TMath::Abs(eta) > etaMax_fidu) || (phi < phiMinEMCal_fidu || phi > phiMaxEMCal_fidu))
	continue;
    }

    if(fWho==2){
    fGenPromptPhotonSel->Fill(6.5);
    }

    photonlabel = iTr;
    eT = mcpart->E()*(TMath::Sin(mcpart->Theta())); // Transform to transverse Energy

    if(fWho == 1 && fQA)
      fphietaPhotons->Fill(eta,phi,eT);

    //Taking out this part of code since it introduces a bias in the efficiency
    //computation
    //      bool foundmatch=kFALSE;
    //      for(int m=0;m<nTracks && foundmatch==kFALSE;m++){
    //        if(m==iTr)
    //          continue;
    //
    //        matchingtrack = static_cast<AliAODMCParticle*>(fAODMCParticles->At(m));
    //
    //        if(! matchingtrack->IsPrimary())
    //          continue;
    //        if(! matchingtrack->IsPhysicalPrimary())
    //          continue;
    //        if(matchingtrack->GetStatus()> 10 )
    //          continue;
    //
    //        Double_t etamatching = matchingtrack->Eta();
    //        Double_t phimatching = matchingtrack->Phi();
    //
    //        if(TMath::Abs(eta-etamatching)<=fdetacut && TMath::Abs(phi-phimatching)<=fdphicut){
    //          foundmatch=kTRUE;
    //	  if(fWho == 1){
    //	    fphietaOthers->Fill(matchingtrack->Eta(),matchingtrack->Phi(),eT);
    //	    fphietaOthersBis->Fill(matchingtrack->Eta(),matchingtrack->Phi(),matchingtrack->Pt());
    //	  }
    //        }
    //      }
    //
    //      if(foundmatch)
    //        continue;

    radius=0.;
    phip=0., etap=0.;
    sumEiso=0.,sumUE=0.;

    for(int iTrack = 0; iTrack < nTracks; iTrack ++){

      if(iTrack == photonlabel)
	continue;

      mcpp = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));
      if(!mcpp)
	continue;

      if(mcpp->Charge() != 0 && mcpp->GetStatus()<10)
	fPtTracksVSpTNC_MC->Fill(eT,mcpp->Pt());

      if(fIsoMethod == 2 && mcpp->Charge() == 0)
	continue;
      else if(fIsoMethod == 3 && mcpp->Charge() != 0)
	continue;

      // Only final state, physical primary and primary particles
      if(mcpp->GetStatus()>10 || (!mcpp->IsPhysicalPrimary() || !mcpp->IsPrimary()))
	continue;

      mcppmomidx = mcpp->GetMother();
      if(mcppmomidx < 0 || mcppmomidx > nTracks)
	continue;

      if(mcppmomidx == photonlabel)
	continue;

      if(mcpp->E() < 0.3) // Minimal energy for clusters allowed at reconstructed level
	continue;

      phip = mcpp->Phi();
      etap = mcpp->Eta();

      // Isolation and UE measurement

      Double_t phiMin = 0., phiMax = 0., etaMin = 0., etaMax = 0.;

      if(fPeriod != ""){
      	etaMin = fGeom->GetArm1EtaMin()+0.03;
      	etaMax = fGeom->GetArm1EtaMax()-0.03;
      	phiMin = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

      	if(fPeriod.Contains("12") || fPeriod.Contains("13"))
      	  phiMax = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03;
      	else
      	  phiMax = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
      }
      else{
      	phiMin = (4./9.)*TMath::Pi()+0.03;
      	phiMax = TMath::Pi()-0.03;
      	etaMin = -0.67;
      	etaMax = 0.67;
      }

      if((phip < phiMax) && (phip > phiMin) && (etap < etaMax) && (etap > etaMin)){
	radius = TMath::Sqrt(TMath::Power(phip-phi,2)+TMath::Power(etap-eta,2));

	if(radius > fIsoConeRadius){ // The cluster is outside the isolation cone -> add the particle pT to pT_UE

	  if(mcpp->Charge() == 0 && mcpp->GetPdgCode() != 22) // Skipping neutral hadrons
	    continue;
	  else
	    AddParticleToUEMC(sumUE, mcpp, eta, phi);

	}
	else{ // The cluster is inside the isolation cone -> add the particle pT to pT_iso

	  if(mcpp->Charge() != 0)        // Using Pt() for charged particles
	    sumEiso += mcpp->Pt();
	  else{
	    if(mcpp->GetPdgCode() == 22) // Using E_T for photons
	      sumEiso += mcpp->E()*(TMath::Sin(mcpp->Theta()));
	    else // Skipping neutral hadrons
	      continue;
	  }

	}
      }
    }

    CalculateUEDensityMC(eta, phi, sumUE);
    
    if(fWho == 2 && fAreasPerEvent){
      if(!fTPC4Iso)
	ComputeConeAreaInEMCal(eta, phi, isoConeArea);
      else
	ComputeConeAreaInTPC(eta, isoConeArea);
    }

    outputValuesMC[0] = eT;
    outputValuesMC[1] = sumEiso;
    outputValuesMC[2] = sumUE;
    outputValuesMC[3] = mompdg;
    outputValuesMC[4] = eta;
    outputValuesMC[5] = phi;
    outputValuesMC[6] = mcpart->GetLabel();

    if(fWho==1)
      fOutMCTruth->Fill(outputValuesMC);

    if(fWho == 2){
      if(fAreasPerEvent){
	if(!fTPC4Iso)
	  ComputeConeAreaInEMCal(eta, phi, isoConeArea);
	else
	  ComputeConeAreaInTPC(eta, isoConeArea);
      }

      if(fAnalysispPb){
	fPtvsUE_MC->Fill(eT, sumUE);
	fPtvsSumUE_MC->Fill(eT, ((sumEiso-sumUE)*(stdConeArea / isoConeArea))); // For etaBand method, output 2, and with fAreasPerEvent flag on: cone and band areas computed candidate-by-candidate
      }
      else
	fPtvsSum_MC->Fill(eT, (sumEiso*(stdConeArea / isoConeArea)));
    }
  }

  return;
}

  //_________________________________________________________________________

void AliAnalysisTaskEMCALPhotonIsolation::AnalyzeMC_Pythia8(){

  if(!fIsMC) return;
  if(!fStack && !fAODMCParticles){ cout << "No stack saved\n"; return; }

  Double_t stdConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.); // Standard (full) cone area
  Double_t isoConeArea = TMath::Pi()*TMath::Power(fIsoConeRadius, 2.);

  Double_t candidateEnergy = 0., candidateEnergyMax = 0., E_T = 0., sumEiso = 0., sumUE = 0., candidatePhi = 0., candidateEta = 0., radius = 0., particlePhi = 0., particleEta = 0.;
  Double_t etaMax_fidu = 0., phiMinEMCal_fidu = 0., phiMaxEMCal_fidu = 0., etaMax = 0., phiMin = 0., phiMax = 0./*, etaMinDCal_InnerEdge = 0., phiMinDCal = 0. , phiMaxDCal_FullSM = 0., phiMaxDCal = 0.*/;

  if(fAODMCParticles->GetEntries() < 1){
    AliError("number of tracks insufficient");
    return;
  }

  Int_t nTracks = fAODMCParticles->GetEntriesFast();
  Int_t iTrack;

  // AliAODMCParticle *track;
  // Int_t nFinalStateParticles = 0;
  // for(iTrack = 0; iTrack < nTracks; iTrack ++){
  //   track = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));

  //   if(track->MCStatusCode() != 1 || !track->IsPhysicalPrimary()) // Discard non final-state particles and non physical primary
  //     continue;

  //   if(TMath::Abs(track->Eta()) <= 0.87 && track->Charge() != 0)  // Count charged particles in TPC acceptance
  //     nFinalStateParticles ++;
  // }

  AliAODMCParticle *candidate, *particle, *candidateMother;
  Int_t candidatePDG = 0, candidatePhotonLabel = 0, particleMotherLabel = 0, candidateMotherLabel = 0, candidateMotherPDG = 0;

  if ( fisLCAnalysis ) {
    for(iTrack = 0; iTrack < nTracks; iTrack ++){
      E_T = 0., candidatePhi = 0., candidateEta = 0.;

      candidate = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));

      // fGenPromptPhotonSel->Fill(0.5);

      if(candidate->MCStatusCode() != 1)  continue;       // Discard non final-state particles

      // fGenPromptPhotonSel->Fill(1.5);

      if(!candidate->IsPhysicalPrimary()) continue;       // Discard non physical primary particles

      // fGenPromptPhotonSel->Fill(2.5);

      candidatePDG = candidate->GetPdgCode();
      if(candidatePDG != 22) continue;                    // Discard particles which are not photons

      // fGenPromptPhotonSel->Fill(3.5);

      if(fPythiaHeader->ProcessType() != 201 && fPythiaHeader->ProcessType() != 202) continue; // Discard particles which do not come from prompt photon processes

      // fGenPromptPhotonSel->Fill(5.5);

      candidateEta = candidate->Eta();
      candidatePhi = candidate->Phi();

      if(!fTPC4Iso){
	etaMax_fidu = fGeom->GetArm1EtaMax()-0.03-fFiducialCut;

	if(!fPeriod.IsNull()){
	  phiMinEMCal_fidu = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03+fFiducialCut;
	  if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	    phiMaxEMCal_fidu = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03-fFiducialCut;
	  else
	    phiMaxEMCal_fidu = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03-fFiducialCut;
	}
	else{ // If no period set, default for 2011-2013 (2*5 EMCal SM)
	  phiMinEMCal_fidu = (4./9.)*TMath::Pi()+0.03+fFiducialCut;
	  phiMaxEMCal_fidu = TMath::Pi()-0.03-fFiducialCut;
	}
      }
      else{
	etaMax_fidu = 0.87-fFiducialCut;

	if(!fPeriod.IsNull()){
	  phiMinEMCal_fidu = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;
	  if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	    phiMaxEMCal_fidu = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03;
	  else
	    phiMaxEMCal_fidu = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
	}
	else{ // If no period set, default for 2011-2013 (2*5 EMCal SM)
	  phiMinEMCal_fidu = (4./9.)*TMath::Pi()+0.03;
	  phiMaxEMCal_fidu = TMath::Pi()-0.03;
	}
      }

      if((TMath::Abs(candidateEta) > etaMax_fidu) || (candidatePhi < phiMinEMCal_fidu || candidatePhi > phiMaxEMCal_fidu)) // Discard photons outside acceptance
	continue;

      // fGenPromptPhotonSel->Fill(6.5);

      // Retrieving the most energetic photon
      candidateEnergy = candidate->E()*TMath::Sin(candidate->Theta());
      if(candidateEnergy > candidateEnergyMax){
	candidateEnergyMax   = candidateEnergy;
	candidatePhotonLabel = iTrack;
      }
      else
	continue;
    }

    candidate = static_cast<AliAODMCParticle*>(fAODMCParticles->At(candidatePhotonLabel));
    E_T       = candidate->E()*(TMath::Sin(candidate->Theta())); // Transform to transverse Energy

    candidateMotherLabel = candidate->GetMother();
    if(candidateMotherLabel > 0 && candidateMotherLabel < nTracks){
      candidateMother    = static_cast<AliAODMCParticle*>(fAODMCParticles->At(candidateMotherLabel));
      candidateMotherPDG = TMath::Abs(candidateMother->GetPdgCode());
    }
    else
      candidateMotherPDG = candidatePDG;

    if(candidateMotherPDG == 22){ // Discard particles whose mother is not a photon

      // fGenPromptPhotonSel->Fill(4.5);

      radius = 0., particlePhi = 0., particleEta = 0., sumEiso = 0., sumUE = 0.;

      // Isolation and UE measurement
      for(iTrack = 0; iTrack < nTracks; iTrack ++){

	if(iTrack == candidatePhotonLabel) continue; // Do not count the candidate photon as particle contributing to isolation/UE energy

	particle = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));
	if(!particle) continue;

	if(particle->Charge() != 0 && particle->MCStatusCode() == 1)
	  fPtTracksVSpTNC_MC->Fill(E_T,particle->Pt());

	if     (fIsoMethod == 2 && particle->Charge() == 0) continue;                                           // Discard neutral particles for charged-only isolation
	else if(fIsoMethod == 3 && particle->Charge() != 0) continue;                                           // Discard charged particles for neutral-only isolation

	if(particle->MCStatusCode() != 1 || !particle->IsPhysicalPrimary()) continue; // Discard non primary, non "detected", non final-state particles

	particleMotherLabel = particle->GetMother();
	if(particleMotherLabel < 0 || particleMotherLabel > nTracks) continue;
	if(particleMotherLabel == candidatePhotonLabel)              continue;                                // Discard mother if it is the candidate photon
	if(particle->E() < 0.3)                                      continue;                                // Discard particles with energy lower than minimal for clusters at reco level

	particlePhi = particle->Phi();
	particleEta = particle->Eta();

	etaMax = fGeom->GetArm1EtaMax()-0.03;
	if(fPeriod != ""){
	  phiMin = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

	  if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	    phiMax = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03;
	  else
	    phiMax = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
	}
	else{
	  phiMin = (4./9.)*TMath::Pi()+0.03;
	  phiMax = TMath::Pi()-0.03;
	}

	if((TMath::Abs(candidateEta) <= etaMax) && (particlePhi >= phiMin) && (particlePhi <= phiMax)){
	  radius = TMath::Sqrt(TMath::Power(particlePhi-candidatePhi,2)+TMath::Power(particleEta-candidateEta,2));

	  if(radius > fIsoConeRadius){                                  // UE energy
	    if(particle->Charge() == 0 && particle->GetPdgCode() != 22)
	      continue;                                                 // Skipping neutral hadrons
	    else
	      AddParticleToUEMC(sumUE, particle, candidateEta, candidatePhi);
	  }
	  else{                                                         // Cone energy
	    if(particle->Charge() != 0)
	      sumEiso += particle->Pt();
	    else{
	      if(particle->GetPdgCode() == 22)
		sumEiso += particle->E()*(TMath::Sin(particle->Theta()));
	      else
		continue;                                               // Skipping neutral hadrons
	    }
	  }
	}
      }

      CalculateUEDensityMC(candidateEta, candidatePhi, sumUE);

      if(fWho == 2){
	if(fAreasPerEvent){
	  if(!fTPC4Iso)
	    ComputeConeAreaInEMCal(candidateEta, candidatePhi, isoConeArea);
	  else
	    ComputeConeAreaInTPC(candidateEta, isoConeArea);
	}

	if(fAnalysispPb){
	  fPtvsUE_MC->Fill(E_T, sumUE);
	  fPtvsSumUE_MC->Fill(E_T, ((sumEiso-sumUE)*(stdConeArea / isoConeArea))); // For etaBand method, output 2, and with fAreasPerEvent flag on: cone and band areas computed candidate-by-candidate
	}
	else
	  fPtvsSum_MC->Fill(E_T, (sumEiso*(stdConeArea / isoConeArea)));
      }
    }
  }
  else {
    for(iTrack = 0; iTrack < nTracks; iTrack ++){
      E_T = 0., candidatePhi = 0., candidateEta = 0.;

      candidate = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));

      fGenPromptPhotonSel->Fill(0.5);

      if(candidate->MCStatusCode() != 1)  continue;       // Discard non final-state particles

      fGenPromptPhotonSel->Fill(1.5);

      if(!candidate->IsPhysicalPrimary()) continue;       // Discard non physical primary particles

      fGenPromptPhotonSel->Fill(2.5);

      candidatePDG = candidate->GetPdgCode();
      if(candidatePDG != 22) continue;                    // Discard particles which are not photons

      fGenPromptPhotonSel->Fill(3.5);

      candidateMotherLabel = candidate->GetMother();
      if(candidateMotherLabel > 0 && candidateMotherLabel < nTracks){
	candidateMother    = static_cast<AliAODMCParticle*>(fAODMCParticles->At(candidateMotherLabel));
	candidateMotherPDG = TMath::Abs(candidateMother->GetPdgCode());
      }
      else
	candidateMotherPDG = candidatePDG;

      if(candidateMotherPDG != 22) continue;              // Discard particles whose mother is not a photon

      fGenPromptPhotonSel->Fill(4.5);

      if(fPythiaHeader->ProcessType() != 201 && fPythiaHeader->ProcessType() != 202) continue; // Discard particles which do not come from prompt photon processes

      fGenPromptPhotonSel->Fill(5.5);

      candidateEta = candidate->Eta();
      candidatePhi = candidate->Phi();

      if(!fTPC4Iso){
	etaMax_fidu = fGeom->GetArm1EtaMax()-0.03-fFiducialCut;

	if(!fPeriod.IsNull()){
	  phiMinEMCal_fidu = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03+fFiducialCut;
	  if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	    phiMaxEMCal_fidu = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03-fFiducialCut;
	  else
	    phiMaxEMCal_fidu = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03-fFiducialCut;
	}
	else{ // If no period set, default for 2011-2013 (2*5 EMCal SM)
	  phiMinEMCal_fidu = (4./9.)*TMath::Pi()+0.03+fFiducialCut;
	  phiMaxEMCal_fidu = TMath::Pi()-0.03-fFiducialCut;
	}
      }
      else{
	etaMax_fidu = 0.87-fFiducialCut;

	if(!fPeriod.IsNull()){
	  phiMinEMCal_fidu = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;
	  if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	    phiMaxEMCal_fidu = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03;
	  else
	    phiMaxEMCal_fidu = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
	}
	else{ // If no period set, default for 2011-2013 (2*5 EMCal SM)
	  phiMinEMCal_fidu = (4./9.)*TMath::Pi()+0.03;
	  phiMaxEMCal_fidu = TMath::Pi()-0.03;
	}
      }

      if((TMath::Abs(candidateEta) > etaMax_fidu) || (candidatePhi < phiMinEMCal_fidu || candidatePhi > phiMaxEMCal_fidu)) // Discard photons outside acceptance
	continue;

      fGenPromptPhotonSel->Fill(6.5);

      candidatePhotonLabel = iTrack;
      E_T = candidate->E()*(TMath::Sin(candidate->Theta())); // Transform to transverse Energy

      radius = 0., particlePhi = 0., particleEta = 0., sumEiso = 0., sumUE = 0.;

      // Isolation and UE measurement
      for(iTrack = 0; iTrack < nTracks; iTrack ++){

	if(iTrack == candidatePhotonLabel) continue; // Do not count the candidate photon as particle contributing to isolation/UE energy

	particle = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));
	if(!particle) continue;

	if(particle->Charge() != 0 && particle->MCStatusCode() == 1)
	  fPtTracksVSpTNC_MC->Fill(E_T,particle->Pt());

	if     (fIsoMethod == 2 && particle->Charge() == 0) continue;                                           // Discard neutral particles for charged-only isolation
	else if(fIsoMethod == 3 && particle->Charge() != 0) continue;                                           // Discard charged particles for neutral-only isolation

	if(particle->MCStatusCode() != 1 || !particle->IsPhysicalPrimary()) continue; // Discard non primary, non "detected", non final-state particles

	particleMotherLabel = particle->GetMother();
	if(particleMotherLabel < 0 || particleMotherLabel > nTracks) continue;
	if(particleMotherLabel == candidatePhotonLabel)              continue;                                // Discard mother if it is the candidate photon
	if(particle->E() < 0.3)                                      continue;                                // Discard particles with energy lower than minimal for clusters at reco level

	particlePhi = particle->Phi();
	particleEta = particle->Eta();

	etaMax = fGeom->GetArm1EtaMax()-0.03;
	if(fPeriod != ""){
	  phiMin = (fGeom->GetArm1PhiMin())*TMath::DegToRad()+0.03;

	  if(fPeriod.Contains("12") || fPeriod.Contains("13"))
	    phiMax = (fGeom->GetEMCALPhiMax()-20.)*TMath::DegToRad()-0.03;
	  else
	    phiMax = (fGeom->GetEMCALPhiMax())*TMath::DegToRad()-0.03;
	}
	else{
	  phiMin = (4./9.)*TMath::Pi()+0.03;
	  phiMax = TMath::Pi()-0.03;
	}

	if((TMath::Abs(candidateEta) <= etaMax) && (particlePhi >= phiMin) && (particlePhi <= phiMax)){
	  radius = TMath::Sqrt(TMath::Power(particlePhi-candidatePhi,2)+TMath::Power(particleEta-candidateEta,2));

	  if(radius > fIsoConeRadius){                                  // UE energy
	    if(particle->Charge() == 0 && particle->GetPdgCode() != 22)
	      continue;                                                 // Skipping neutral hadrons
	    else
	      AddParticleToUEMC(sumUE, particle, candidateEta, candidatePhi);
	  }
	  else{                                                         // Cone energy
	    if(particle->Charge() != 0)
	      sumEiso += particle->Pt();
	    else{
	      if(particle->GetPdgCode() == 22)
		sumEiso += particle->E()*(TMath::Sin(particle->Theta()));
	      else
		continue;                                               // Skipping neutral hadrons
	    }
	  }
	}
      }

      CalculateUEDensityMC(candidateEta, candidatePhi, sumUE);

      if(fWho == 2){
	if(fAreasPerEvent){
	  if(!fTPC4Iso)
	    ComputeConeAreaInEMCal(candidateEta, candidatePhi, isoConeArea);
	  else
	    ComputeConeAreaInTPC(candidateEta, isoConeArea);
	}

	if(fAnalysispPb){
	  fPtvsUE_MC->Fill(E_T, sumUE);
	  fPtvsSumUE_MC->Fill(E_T, ((sumEiso-sumUE)*(stdConeArea / isoConeArea))); // For etaBand method, output 2, and with fAreasPerEvent flag on: cone and band areas computed candidate-by-candidate
	}
	else
	  fPtvsSum_MC->Fill(E_T, (sumEiso*(stdConeArea / isoConeArea)));
      }
    }
  }

  return;
}
