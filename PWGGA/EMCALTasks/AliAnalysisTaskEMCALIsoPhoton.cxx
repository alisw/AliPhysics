// $Id$

#include "AliAnalysisTaskEMCALIsoPhoton.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TLorentzVector.h>
#include <TList.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliDataFile.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDUtils.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliV0vertexer.h"
#include "AliVCluster.h"
#include "AliOADBContainer.h"


#include <iostream>
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEMCALIsoPhoton)

//________________________________________________________________________
AliAnalysisTaskEMCALIsoPhoton::AliAnalysisTaskEMCALIsoPhoton() : 
  AliAnalysisTaskSE(), 
  fESDClusters(0),
  fAODClusters(0),
  fSelPrimTracks(0),
  fTracks(0),
  fAODMCParticles(0),
  fESDCells(0),
  fAODCells(0),
  fVCells(0),
  fPrTrCuts(0),
  fCompTrCuts(0),
  fGeom(0x0),
  fGeoName("EMCAL_COMPLETEV1"),
  fOADBContainer(0),
  fVecPv(0.,0.,0.),
  fPeriod("LHC11c"),
  fTrigBit("kEMC7"),
  fIsTrain(0),
  fIsMc(0),
  fDebug(0),
  fPathStrOpt("/"),
  fExoticCut(0.97),
  fIsoConeR(0.4),
  fNDimensions(7),
  fECut(3.),
  fTrackMult(0),        
  fMcIdFamily(""),
  fNClusForDirPho(0),
  fDirPhoPt(0),
  fHigherPtCone(0),
  fImportGeometryFromFile(0),
  fImportGeometryFilePath(""),
  fMaxPtTrack(0),
  fMaxEClus(0),
  fNCells50(0),
  fFilterBit(0),
  fSelHybrid(kFALSE),
  fFillQA(kFALSE),
  fClusIdFromTracks(""),
  fCpvFromTrack(kFALSE),
  fNBinsPt(200),
  fPtBinLowEdge(0),
  fPtBinHighEdge(100),
  fRemMatchClus(kFALSE),
  fMinIsoClusE(0),
  fNCuts(5),
  fTrCoreRem(kFALSE),
  fClusTDiff(30e-9),
  fPileUpRejSPD(kFALSE),
  fDistToBadChan(0),
  fInConeInvMass(""),
  fInConePairClEt(""),
  fNSigNeutMesonCut(2.0),
  fSigmaSmear(0.0),
  fNLMCut(100),
  fESD(0),
  fAOD(0),
  fVEvent(0),
  fMCEvent(0),
  fStack(0),
  fOutputList(0),
  fEvtSel(0),
  fNClusEt10(0),
  fClusArrayNames(0),
  fRecoPV(0),
  fPVtxZ(0),
  fTrMultDist(0),
  fClusEtCPVSBGISO(0),
  fClusEtCPVBGISO(0),
  fMCDirPhotonPtEtaPhi(0),
  fMCIsoDirPhotonPtEtaPhi(0),
  fMCDirPhotonPtEtIso(0),
  fDecayPhotonPtMC(0),
  fCellAbsIdVsAmpl(0),       
  fNClusHighClusE(0),
  fHigherPtConeM02(0),
  fClusEtMcPt(0),
  fClusMcDetaDphi(0),
  fNClusPerPho(0),
  fMcPtInConeBG(0),
  fMcPtInConeSBG(0),
  fMcPtInConeBGnoUE(0),
  fMcPtInConeSBGnoUE(0),
  fMcPtInConeTrBGnoUE(0),
  fMcPtInConeTrSBGnoUE(0),
  fMcPtInConeMcPhoPt(0),
  fAllIsoEtMcGamma(0),
  fAllIsoNoUeEtMcGamma(0),
  fMCDirPhotonPtEtaPhiNoClus(0),
  fEtCandIsoAndIsoWoPairEt(0),
  fInConePairedClusEtVsCandEt(0),
  fHnOutput(0),
  fQAList(0),
  fNTracks(0),     
  fEmcNCells(0),   
  fEmcNClus(0),    
  fEmcNClusCut(0), 
  fNTracksECut(0), 
  fEmcNCellsCut(0),
  fEmcClusETM1(0),
  fEmcClusETM2(0),
  fEmcClusNotExo(0),
  fEmcClusEClusCuts(0),
  fEmcClusEPhi(0),    
  fEmcClusEPhiCut(0), 
  fEmcClusEEta(0),    
  fEmcClusEEtaCut(0), 
  fTrackPtPhi(0),     
  fTrackPtPhiCut(0),   
  fTrackPtEta(0),     
  fTrackPtEtaCut(0),
  fMaxCellEPhi(0),
  fDetaDphiFromTM(0),
  fEoverPvsE(0),
  fTrackDEtaDPhiPho(0),
  fTrackDEtaDPhiPi(0),
  fTrackDEtaDPhiPPho(0),
  fTrackDzDxIM(0),
  fTrackDzDxPhoSS(0),
  fTrackDzDxIM_bg(0),
  fTrackDzDxPhoSS_bg(0),
  fETrigg(0),
  fM02vsESoftPi0Kid(0),
  fM02vsEHardPi0Kid(0),
  fM02vsESoftPi0BGKid(0),
  fM02vsEHardPi0BGKid(0),
  fClusInvMassPairEt(0),
  fPi0KidClusEtPairEt(0),
  fCellsPi0KidE6(0),
  fCellsPi0KidE11(0),
  fCellsPi0KidE13(0),
  fCellsPi0KidM021st(0),
  fCellsPi0KidM022nd(0),
  fCellsPi0KidM023rd(0),
  fCellsPi0KidM024th(0)
{
  // Default constructor.
  for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
}

//________________________________________________________________________
AliAnalysisTaskEMCALIsoPhoton::AliAnalysisTaskEMCALIsoPhoton(const char *name) : 
  AliAnalysisTaskSE(name), 
  fESDClusters(0),
  fAODClusters(0),
  fSelPrimTracks(0),
  fTracks(0),
  fAODMCParticles(0),
  fESDCells(0),
  fAODCells(0),
  fVCells(0),
  fPrTrCuts(0),
  fCompTrCuts(0),
  fGeom(0x0),
  fGeoName("EMCAL_COMPLETEV1"),
  fOADBContainer(0),
  fVecPv(0.,0.,0.),
  fPeriod("LHC11c"),
  fTrigBit("kEMC7"),
  fIsTrain(0),
  fIsMc(0),
  fDebug(0),
  fPathStrOpt("/"),
  fExoticCut(0.97),
  fIsoConeR(0.4),
  fNDimensions(7),
  fECut(3.),
  fTrackMult(0),        
  fMcIdFamily(""),
  fNClusForDirPho(0),
  fDirPhoPt(0),
  fHigherPtCone(0),
  fImportGeometryFromFile(0),
  fImportGeometryFilePath(""),
  fMaxPtTrack(0),
  fMaxEClus(0),
  fNCells50(0),
  fFilterBit(0),
  fSelHybrid(kFALSE),
  fFillQA(kFALSE),
  fClusIdFromTracks(""),
  fCpvFromTrack(kFALSE),
  fNBinsPt(200),
  fPtBinLowEdge(0.),
  fPtBinHighEdge(100),
  fRemMatchClus(kFALSE),
  fMinIsoClusE(0),
  fNCuts(5),
  fTrCoreRem(kFALSE),
  fClusTDiff(30e-9),
  fPileUpRejSPD(kFALSE),
  fDistToBadChan(0),
  fInConeInvMass(""),
  fInConePairClEt(""),
  fNSigNeutMesonCut(2.0),
  fSigmaSmear(0.0),
  fNLMCut(100),
  fESD(0),
  fAOD(0),
  fVEvent(0),
  fMCEvent(0),
  fStack(0),
  fOutputList(0),
  fEvtSel(0),
  fNClusEt10(0),
  fClusArrayNames(0),
  fRecoPV(0),
  fPVtxZ(0),            
  fTrMultDist(0),
  fClusEtCPVSBGISO(0),
  fClusEtCPVBGISO(0),
  fMCDirPhotonPtEtaPhi(0),
  fMCIsoDirPhotonPtEtaPhi(0),
  fMCDirPhotonPtEtIso(0),
  fDecayPhotonPtMC(0),
  fCellAbsIdVsAmpl(0),       
  fNClusHighClusE(0),   
  fHigherPtConeM02(0),
  fClusEtMcPt(0),
  fClusMcDetaDphi(0),
  fNClusPerPho(0),
  fMcPtInConeBG(0),
  fMcPtInConeSBG(0),
  fMcPtInConeBGnoUE(0),
  fMcPtInConeSBGnoUE(0),
  fMcPtInConeTrBGnoUE(0),
  fMcPtInConeTrSBGnoUE(0),
  fMcPtInConeMcPhoPt(0),
  fAllIsoEtMcGamma(0),
  fAllIsoNoUeEtMcGamma(0),
  fMCDirPhotonPtEtaPhiNoClus(0),
  fEtCandIsoAndIsoWoPairEt(0),
  fInConePairedClusEtVsCandEt(0),
  fHnOutput(0),
  fQAList(0),
  fNTracks(0),     
  fEmcNCells(0),   
  fEmcNClus(0),    
  fEmcNClusCut(0), 
  fNTracksECut(0), 
  fEmcNCellsCut(0),
  fEmcClusETM1(0),
  fEmcClusETM2(0),
  fEmcClusNotExo(0),
  fEmcClusEClusCuts(0),
  fEmcClusEPhi(0),    
  fEmcClusEPhiCut(0), 
  fEmcClusEEta(0),    
  fEmcClusEEtaCut(0), 
  fTrackPtPhi(0),     
  fTrackPtPhiCut(0),   
  fTrackPtEta(0),     
  fTrackPtEtaCut(0),   
  fMaxCellEPhi(0),
  fDetaDphiFromTM(0),
  fEoverPvsE(0),
  fTrackDEtaDPhiPho(0),
  fTrackDEtaDPhiPi(0),
  fTrackDEtaDPhiPPho(0),
  fTrackDzDxIM(0),
  fTrackDzDxPhoSS(0),
  fTrackDzDxIM_bg(0),
  fTrackDzDxPhoSS_bg(0),
  fETrigg(0),
  fM02vsESoftPi0Kid(0),
  fM02vsEHardPi0Kid(0),
  fM02vsESoftPi0BGKid(0),
  fM02vsEHardPi0BGKid(0),
  fClusInvMassPairEt(0),
  fPi0KidClusEtPairEt(0),
  fCellsPi0KidE6(0),
  fCellsPi0KidE11(0),
  fCellsPi0KidE13(0),
  fCellsPi0KidM021st(0),
  fCellsPi0KidM022nd(0),
  fCellsPi0KidM023rd(0),
  fCellsPi0KidM024th(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::UserCreateOutputObjects()
{
  // Create histograms, called once.
    
  fESDClusters = new TObjArray();
  fAODClusters = new TObjArray();
  fSelPrimTracks = new TObjArray();

  
  fOutputList = new TList();
  fOutputList->SetOwner();// Container cleans up all histos (avoids leaks in merging) 
  
  fGeom = AliEMCALGeometry::GetInstance(fGeoName.Data());
  fOADBContainer = new AliOADBContainer("AliEMCALgeo");
  fOADBContainer->InitFromFile(AliDataFile::GetFileNameOADB("EMCAL/EMCALlocal2master.root").data(),"AliEMCALgeo");
 
  fEvtSel = new TH1F("hEvtSel","Event selection counter (0=all trg, 1=pvz cut) ;evt cut ;dN/dcut}",2,0,2);
  fOutputList->Add(fEvtSel);
  
  fNClusEt10 = new TH1F("hNClusEt10","# of cluster with E_{T}>10 per event;E;",101,-0.5,100.5);
  fOutputList->Add(fNClusEt10);

  fClusArrayNames = new TH1F("hClusArrayNames","cluster array names (0=CaloClusters,1=EmcCaloClusters,2=Others);option;#events",3,0,3);
  fOutputList->Add(fClusArrayNames);
  
  fRecoPV = new TH1F("hRecoPV","Prim. vert. reconstruction (yes or no);reco (0=no, 1=yes) ;",2,-0.5,1.5);
  fOutputList->Add(fRecoPV);

  fPVtxZ = new TH1F("hPVtxZ","primary vertex Z before cut;prim-vz(cm) ;",200,-100,100);
  fOutputList->Add(fPVtxZ);

  fTrMultDist = new TH1F("fTrMultDist","track multiplicity;tracks/event;#events",200,0.5,200.5);
  fOutputList->Add(fTrMultDist);

  fClusEtCPVSBGISO = new TH2F("hClusEtCPVSBGISO","ISO^{TRK+EMC} vs. E_{T}^{clus} (after CPV and 0.1<#lambda_{0}^{2}<0.3;E_{T}^{clus} [GeV];ISO^{TRK+EMC} [GeV]",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,1000,0,100);
  fClusEtCPVSBGISO->Sumw2();
  fOutputList->Add(fClusEtCPVSBGISO);

  fClusEtCPVBGISO = new TH2F("hClusEtCPVBGISO","ISO^{TRK+EMC} vs. E_{T}^{clus} (after CPV and 0.5<#lambda_{0}^{2}<2.0;E_{T}^{clus} [GeV];ISO^{TRK+EMC} [GeV]",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,1000,0,100);
  fClusEtCPVBGISO->Sumw2();
  fOutputList->Add(fClusEtCPVBGISO);

  fMCDirPhotonPtEtaPhi = new TH3F("hMCDirPhotonPtEtaPhi","photon (gq->#gammaq) p_{T}, #eta, #phi;GeV/c;#eta;#phi",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,154,-0.77,0.77,130,1.38,3.20);
  fMCDirPhotonPtEtaPhi->Sumw2();
  fOutputList->Add(fMCDirPhotonPtEtaPhi);

  fMCIsoDirPhotonPtEtaPhi = new TH3F("hMCIsoDirPhotonPtEtaPhi","photon (gq->#gammaq, isolated@MC) p_{T}, #eta, #phi;GeV/c;#eta;#phi",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,154,-0.77,0.77,130,1.38,3.20);
  fMCIsoDirPhotonPtEtaPhi->Sumw2();
  fOutputList->Add(fMCIsoDirPhotonPtEtaPhi);

  fMCDirPhotonPtEtIso = new TH2F("hMCDirPhotonPtEtIso",Form("photon (gq->#gammaq @MC) p_{T}, E_{T}^{ISO} (R=%1.1f);GeV/c;E_{T}^{ISO} GeV/c",fIsoConeR),fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,20,-0.25,9.75);
  fMCDirPhotonPtEtIso->Sumw2();
  fOutputList->Add(fMCDirPhotonPtEtIso);


  fDecayPhotonPtMC = new TH1F("hDecayPhotonPtMC","decay photon p_{T};GeV/c;dN/dp_{T} (c GeV^{-1})",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge);
  fDecayPhotonPtMC->Sumw2();
  fOutputList->Add(fDecayPhotonPtMC);

  fCellAbsIdVsAmpl = new TH2F("hCellAbsIdVsAmpl","cell abs id vs cell amplitude (energy);E (GeV);ID",200,0,100,24*48*10,-0.5,24*48*10-0.5);
  fOutputList->Add(fCellAbsIdVsAmpl);

  fNClusHighClusE = new TH2F("hNClusHighClusE","total number of clusters vs. highest clus energy in the event;E (GeV);NClus",200,0,100,301,-0.5,300.5);
  fOutputList->Add(fNClusHighClusE);

  fHigherPtConeM02 = new TH2F("hHigherPtConeM02","#lambda_{0}^{2} vs. in-cone-p_{T}^{max};p_{T}^{max} (GeV/c, in the cone);#lambda_{0}^{2}",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,400,0,4);
  fOutputList->Add(fHigherPtConeM02);

  fClusEtMcPt = new TH2F("hClusEtMcPt","E_{T}^{clus} vs. p_{T}^{mc}; p_{T}^{mc};E_{T}^{clus}",500,0,100,500,0,100);
  fOutputList->Add(fClusEtMcPt);

  fClusMcDetaDphi = new TH2F("hClusMcDetaDphi","#Delta#phi vs. #Delta#eta (reco-mc);#Delta#eta;#Delta#phi",100,-.7,.7,100,-.7,.7);
  fOutputList->Add(fClusMcDetaDphi);

  fNClusPerPho = new TH2F("hNClusPerPho","Number of clusters per prompt photon;p_{T}^{MC};N_{clus}",500,0,100,11,-0.5,10.5);
  fOutputList->Add(fNClusPerPho);

  fMcPtInConeBG = new TH2F("hMcPtInConeBG","#sum_{in-cone}p_{T}^{mc-primaries} vs. ISO^{TRK+EMC} (BG template);ISO^{TRK+EMC} (GeV);#sum_{in-cone}p_{T}^{mc-primaries}",600,-10,50,1000,0,100);
  fOutputList->Add(fMcPtInConeBG);

  fMcPtInConeSBG  = new TH2F("hMcPtInConeSBG","#sum_{in-cone}p_{T}^{mc-primaries} vs. ISO^{TRK+EMC} (SBG range);ISO^{TRK+EMC} (GeV);#sum_{in-cone}p_{T}^{mc-primaries}",600,-10,50,1000,0,100);
  fOutputList->Add(fMcPtInConeSBG);

  fMcPtInConeBGnoUE = new TH2F("hMcPtInConeBGnoUE","#sum_{in-cone}p_{T}^{mc-primaries} vs. ISO^{TRK+EMC} (BG template);ISO^{TRK+EMC} (GeV);#sum_{in-cone}p_{T}^{mc-primaries}",600,-10,50,1000,0,100);
  fOutputList->Add(fMcPtInConeBGnoUE);

  fMcPtInConeSBGnoUE  = new TH2F("hMcPtInConeSBGnoUE","#sum_{in-cone}p_{T}^{mc-primaries} vs. ISO^{TRK+EMC} (SBG range);ISO^{TRK+EMC} (GeV);#sum_{in-cone}p_{T}^{mc-primaries}",600,-10,50,1000,0,100);
  fOutputList->Add(fMcPtInConeSBGnoUE);

  fMcPtInConeTrBGnoUE = new TH2F("hMcPtInConeTrBGnoUE","#sum_{in-cone}p_{T}^{mc-primaries} vs. ISO^{TRK} (BG template);ISO^{TRK} (GeV);#sum_{in-cone}p_{T}^{mc-primaries}",600,-10,50,1000,0,100);
  fOutputList->Add(fMcPtInConeTrBGnoUE);

  fMcPtInConeTrSBGnoUE  = new TH2F("hMcPtInConeTrSBGnoUE","#sum_{in-cone}p_{T}^{mc-primaries} vs. ISO^{TRK} (SBG range);ISO^{TRK} (GeV);#sum_{in-cone}p_{T}^{mc-primaries}",600,-10,50,1000,0,100);
  fOutputList->Add(fMcPtInConeTrSBGnoUE);

  fMcPtInConeMcPhoPt  = new TH2F("hMcPtInConeMcPhoPt","#sum_{in-cone}p_{T}^{mc-primaries} vs. prompt photon p_{T};p_{T}^{mc-#gamma} (GeV);#sum_{in-cone}p_{T}^{mc-primaries}",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,200,-0.25,99.75);
  fOutputList->Add(fMcPtInConeMcPhoPt);

  fAllIsoEtMcGamma  = new TH2F("hAllIsoEtMcGammaE","ISO^{TRK+EMC} vs. E_{T}^{clus} for clusters comming from a MC prompt #gamma; E_{T}^{clus} (GeV);ISO^{TRK+EMC} (GeV);",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,600,-10,50);
  fOutputList->Add(fAllIsoEtMcGamma);

  fAllIsoNoUeEtMcGamma  = new TH2F("hAllIsoNoUeEtMcGammaE","ISO^{TRK+EMC}_{noue} vs. E_{T}^{clus} for clusters comming from a MC prompt #gamma; E_{T}^{clus} (GeV);ISO^{TRK+EMC}_{noue} (GeV);",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,600,-10,50);
  fOutputList->Add(fAllIsoNoUeEtMcGamma);


  fMCDirPhotonPtEtaPhiNoClus = new TH3F("hMCDirPhotonPhiEtaNoClus","p_{T}, #eta and  #phi of prompt photons with no reco clusters;p_{T};#eta;#phi",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,154,-0.77,0.77,130,1.38,3.20);
  fOutputList->Add(fMCDirPhotonPtEtaPhiNoClus);

  fEtCandIsoAndIsoWoPairEt = new TH3F("hEtCandIsoAndIsoWoPairEt","E_{T}^{cand} vs. E_{T}^{ISO} (EMC+Trk) (0.1<M02<0.3, 0.110<m_{#gamma#gamma}<0.165 only);E_{T}^{cand}; E_{T}^{ISO}; E_{T}^{ISO} (w/o #pi^{0} pair E_{T})",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,1000,0,200,1000,0,200);
  fOutputList->Add(fEtCandIsoAndIsoWoPairEt);

  fInConePairedClusEtVsCandEt = new TH2F("hInConePairedClusEtVsCandEt","E_{T}^{partner} vs. E_{T}^{cand} (R<0.4, 0.110<m_{#gamma#gamma}<0.165);E_{T}^{cand};E_{T}^{partner}",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,200,0,40);
  fOutputList->Add(fInConePairedClusEtVsCandEt);

  Int_t nEt=fNBinsPt*5, nM02=400, nCeIso=1000, nTrIso=1000,  nAllIso=1000,  nCeIsoNoUE=1000,  nAllIsoNoUE=1000, nTrClDphi=200, nTrClDeta=100, nClEta=140, nClPhi=128, nTime=60, nMult=100, nPhoMcPt=fNBinsPt, nInConeMass=100, nNLM=11;
  Int_t bins[] = {nEt, nM02, nCeIso, nTrIso, nAllIso, nCeIsoNoUE, nAllIsoNoUE, nTrClDphi, nTrClDeta,nClEta,nClPhi,nTime,nMult,nPhoMcPt,nInConeMass, nNLM};
  fNDimensions = sizeof(bins)/sizeof(Int_t);
  const Int_t ndims =   fNDimensions;
  Double_t xmin[] = { fPtBinLowEdge,   0.,  -10.,   -10., -10., -10., -10., -0.1,-0.05, -0.7, 1.4,-0.15e-06,-0.5,fPtBinLowEdge,0.0,-0.5};
  Double_t xmax[] = { fPtBinHighEdge, 4., 190., 190., 190.,  190., 190., 0.1, 0.05, 0.7, 3.192, 0.15e-06,99.5,fPtBinHighEdge, 1.0,10.5};
  if(fPeriod.Contains("11h")){
    xmax[12]=3999.5;
  }
  fHnOutput =  new THnSparseF("fHnOutput","Output matrix: E_{T},M02,CeIso,TrIso,AllIso, CeIsoNoUESub, AllIsoNoUESub, d#phi_{trk},d#eta_{trk},#eta_{clus},#phi_{clus},T_{max},mult,mc-p_{T}^{#gamma},NLM", ndims, bins, xmin, xmax);
  fHnOutput->Sumw2();
  fOutputList->Add(fHnOutput);

  //QA outputs
  fQAList = new TList();
  fQAList->SetOwner();// Container cleans up all histos (avoids leaks in merging) 

  fNTracks = new TH1F("hNTracks","# of selected tracks;n-tracks;counts",120,-0.5,119.5);
  fNTracks->Sumw2();
  fQAList->Add(fNTracks);

  fEmcNCells = new TH1F("fEmcNCells",";n/event;count",120,-0.5,119.5);  
  fEmcNCells->Sumw2();
  fQAList->Add(fEmcNCells);
  fEmcNClus = new TH1F("fEmcNClus",";n/event;count",120,-0.5,119.5);    
  fEmcNClus->Sumw2();
  fQAList->Add(fEmcNClus);
  fEmcNClusCut = new TH1F("fEmcNClusCut",Form("(at least one E_{clus}>%1.1f);n/event;count",fECut),120,-0.5,119.5); 
  fEmcNClusCut->Sumw2();
  fQAList->Add(fEmcNClusCut);
  fNTracksECut = new TH1F("fNTracksECut",Form("(at least one E_{clus}>%1.1f);n/event;count",fECut),120,-0.5,119.5); 
  fNTracksECut->Sumw2();
  fQAList->Add(fNTracksECut);
  fEmcNCellsCut = new TH1F("fEmcNCellsCut",Form("(at least one E_{clus}>%1.1f);n/event;count",fECut),120,-0.5,119.5);
  fEmcNCellsCut->Sumw2();
  fQAList->Add(fEmcNCellsCut);
  fEmcClusETM1 = new TH1F("fEmcClusETM1","(method clus->GetTrackDx,z);GeV;counts",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge);
  fEmcClusETM1->Sumw2();
  fQAList->Add(fEmcClusETM1);
  fEmcClusETM2 = new TH1F("fEmcClusETM2","(method track->GetEMCALcluster());GeV;counts",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge);
  fEmcClusETM2->Sumw2();
  fQAList->Add(fEmcClusETM2);
  fEmcClusNotExo  = new TH1F("fEmcClusNotExo","exotics removed;GeV;counts",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge);
  fEmcClusNotExo->Sumw2();
  fQAList->Add(fEmcClusNotExo);
  fEmcClusEPhi = new TH2F("fEmcClusEPhi",";GeV;#phi",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,63,0,6.3);    
  fEmcClusEPhi->Sumw2();
  fQAList->Add(fEmcClusEPhi);
  fEmcClusEPhiCut = new TH2F("fEmcClusEPhiCut",Form("(at least one E_{clus}>%1.1f);GeV;#phi",fECut),fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,63,0,6.3); 
  fEmcClusEPhiCut->Sumw2();
  fQAList->Add(fEmcClusEPhiCut);
  fEmcClusEEta = new TH2F("fEmcClusEEta",";GeV;#eta",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,19,-0.9,0.9);    
  fEmcClusEEta->Sumw2();
  fQAList->Add(fEmcClusEEta);
  fEmcClusEEtaCut = new TH2F("fEmcClusEEtaCut",Form("(at least one E_{clus}>%1.1f);GeV;#eta",fECut),fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,18,-0.9,0.9); 
  fEmcClusEEtaCut->Sumw2();
  fQAList->Add(fEmcClusEEtaCut);
  fTrackPtPhi = new TH2F("fTrackPtPhi",";p_{T} [GeV/c];#phi",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,63,0,6.3);     
  fTrackPtPhi->Sumw2();
  fQAList->Add(fTrackPtPhi);
  fTrackPtPhiCut = new TH2F("fTrackPtPhiCut",Form("(at least one E_{clus}>%1.1f);p_{T} [GeV/c];#phi",fECut),fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,63,0,6.3);     
  fTrackPtPhiCut->Sumw2();
  fQAList->Add(fTrackPtPhiCut);
  fTrackPtEta = new TH2F("fTrackPtEta",";p_{T} [GeV/c];#eta",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,18,-0.9,0.9);     
  fTrackPtEta->Sumw2();
  fQAList->Add(fTrackPtEta);
  fTrackPtEtaCut = new TH2F("fTrackPtEtaCut",Form("(at least one E_{clus}>%1.1f);p_{T} [GeV/c];#eta",fECut),fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,18,-0.9,0.9);     
  fTrackPtEtaCut->Sumw2();
  fQAList->Add(fTrackPtEtaCut);
  fEmcClusEClusCuts = new TH2F("fEmcClusEClusCuts",";GeV;cut",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,fNCuts,-0.5,fNCuts-0.5);
  fEmcClusEClusCuts->Sumw2();
  fQAList->Add(fEmcClusEClusCuts);

  fMaxCellEPhi = new TH2F("fMaxCellEPhi","Most energetic cell in event; GeV;#phi",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,63,0,6.3); 
  fMaxCellEPhi->Sumw2();
  fQAList->Add(fMaxCellEPhi);

  fDetaDphiFromTM = new TH2F("fDetaDphiFromTM","d#phi vs. d#eta of clusters from track->GetEMCALcluster();d#eta;d#phi",100,-0.05,0.05,200,-0.1,0.1);
  fDetaDphiFromTM->Sumw2();
  fQAList->Add(fDetaDphiFromTM);

  fEoverPvsE = new TH2F("fEoverPvsE","E^{clus}/p^{track} vs E^{clus} (80<TPCsignal<100);E^{clus} [GeV];E^{clus}/p^{track} [c^{-1}]",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,100,0,2);
  fEoverPvsE->Sumw2();
  fQAList->Add(fEoverPvsE);

  fTrackDEtaDPhiPho = new TH2F("fTrackDEtaDPhiPho","clusters from MC-truth #pi^{0}#rightarrow #gamma ;#Delta#eta_{Tr-Cl};#Delta#phi_{Tr-Cl}",100,-0.05,0.05,200,-0.1,0.1);
  fTrackDEtaDPhiPho->Sumw2();
  fQAList->Add(fTrackDEtaDPhiPho);

  fTrackDEtaDPhiPi = new TH2F("fTrackDEtaDPhiPi","clusters from MC-truth #pi^{#pm};#Delta#eta_{Tr-Cl};#Delta#phi_{Tr-Cl}",100,-0.05,0.05,200,-0.1,0.1);
  fTrackDEtaDPhiPi->Sumw2();
  fQAList->Add(fTrackDEtaDPhiPi);

  fTrackDEtaDPhiPPho = new TH2F("fTrackDEtaDPhiPPho","clusters from MC-truth prompt-#gamma ;#Delta#eta_{Tr-Cl};#Delta#phi_{Tr-Cl}",100,-0.05,0.05,200,-0.1,0.1);
  fTrackDEtaDPhiPPho->Sumw2();
  fQAList->Add(fTrackDEtaDPhiPPho);

  fTrackDzDxIM = new TH2F("fTrackDzDxIM","cluster in #pi^{0} invariant mass range (inclusive energy);#Delta#eta_{Tr-Cl};#Delta#phi_{Tr-Cl}",100,-0.05,0.05,200,-0.1,0.1);
  fTrackDzDxIM->Sumw2();
  fQAList->Add(fTrackDzDxIM);

  fTrackDzDxPhoSS = new TH2F("fTrackDzDxPhoSS","cluster in #pi^{0} #lambda_{0}^{2} range (inclusive energy);#Delta#eta_{Tr-Cl};#Delta#phi_{Tr-Cl}",100,-0.05,0.05,200,-0.1,0.1);
  fTrackDzDxPhoSS->Sumw2();
  fQAList->Add(fTrackDzDxPhoSS);

  fTrackDzDxIM_bg = new TH2F("fTrackDzDxIM_bg","cluster in side bands of #pi^{0} mass (inclusive energy);#Delta#eta_{Tr-Cl};#Delta#phi_{Tr-Cl}",100,-0.05,0.05,200,-0.1,0.1);
  fTrackDzDxIM_bg->Sumw2();
  fQAList->Add(fTrackDzDxIM_bg);

  fTrackDzDxPhoSS_bg =  new TH2F("fTrackDzDxPhoSS_bg","cluster in side bands of merged #pi^{0} #lambda_{0}^{2} (inclusive energy);#Delta#eta_{Tr-Cl};#Delta#phi_{Tr-Cl}",100,-0.05,0.05,200,-0.1,0.1);
  fTrackDzDxPhoSS_bg->Sumw2();
  fQAList->Add(fTrackDzDxPhoSS_bg);

  fETrigg = new TH1F("fETrigg","TrigPatchInfo->GetPatchE();E [GeV];entries/GeV",100,0,100);
  fETrigg->Sumw2();
  fQAList->Add(fETrigg);
  fM02vsESoftPi0Kid = new TH2F("fM02vsESoftPi0Kid","#lambda_{0}^{2} vs E_{clus} (softer #pi^{0} daughter);E_{clus};#lambda_{0}^{2}",fNBinsPt,fPtBinLowEdge,fPtBinHighEdge,400,0,4);
  fM02vsESoftPi0Kid->Sumw2();
  fQAList->Add(fM02vsESoftPi0Kid);
  fM02vsEHardPi0Kid = new TH2F("fM02vsEHardPi0Kid","#lambda_{0}^{2} vs E_{clus} (harder #pi^{0} daughter);E_{clus};#lambda_{0}^{2}",fNBinsPt,fPtBinLowEdge,fPtBinHighEdge,400,0,4);
  fM02vsEHardPi0Kid->Sumw2();
  fQAList->Add(fM02vsEHardPi0Kid);
  fM02vsESoftPi0BGKid = new TH2F("fM02vsESoftPi0BGKid","#lambda_{0}^{2} vs E_{clus} (softer #pi^{0}-BG daughter);E_{clus};#lambda_{0}^{2}",fNBinsPt,fPtBinLowEdge,fPtBinHighEdge,400,0,4);
  fM02vsESoftPi0BGKid->Sumw2();
  fQAList->Add(fM02vsESoftPi0BGKid);
  fM02vsEHardPi0BGKid = new TH2F("fM02vsEHardPi0BGKid","#lambda_{0}^{2} vs E_{clus} (harder #pi^{0}-BG daughter);E_{clus};#lambda_{0}^{2}",fNBinsPt,fPtBinLowEdge,fPtBinHighEdge,400,0,4);
  fM02vsEHardPi0BGKid->Sumw2();
  fQAList->Add(fM02vsEHardPi0BGKid);

  fClusInvMassPairEt = new TH2F("fClusInvMassPairEt","neutral cluster inv mass vs. E_{T}^{pair};E_{T}^{pair};mass_{#gamma#gamma}",fNBinsPt,fPtBinLowEdge,fPtBinHighEdge,100,0,1);
  fClusInvMassPairEt->Sumw2();
  fQAList->Add(fClusInvMassPairEt);

  fPi0KidClusEtPairEt = new TH2F("fPi0KidClusEtPairEt","E_{T}^{clus} vs. E_{T}^{pair};E_{T}^{pair};E_{T}^{clus}",fNBinsPt,fPtBinLowEdge,fPtBinHighEdge,fNBinsPt,fPtBinLowEdge,fPtBinHighEdge);
  fPi0KidClusEtPairEt->Sumw2();
  fQAList->Add(fPi0KidClusEtPairEt);

  fCellsPi0KidE6 = new TH2F("fCellsPi0KidE6","Distribution of cluster energy in its cells (5.9<E_{clus}<6.1,#pi^{0} daughters, all #lambda_{0}^{2};col;row",19,-9,9,19,-9,9);
  fCellsPi0KidE6->Sumw2();
  fQAList->Add(fCellsPi0KidE6);

  fCellsPi0KidE11  = new TH2F("fCellsPi0KidE11","Distribution of cluster energy in its cells (10.5<E_{clus}<11.5,#pi^{0} daughters, all #lambda_{0}^{2};col;row",19,-9,9,19,-9,9);
  fCellsPi0KidE11->Sumw2();
  fQAList->Add(fCellsPi0KidE11);

  fCellsPi0KidE13 = new TH2F("fCellsPi0KidE13","Distribution of cluster energy in its cells (12.5<E_{clus}<13.5,#pi^{0} daughters, all #lambda_{0}^{2};col;row",19,-9,9,19,-9,9);
  fCellsPi0KidE13->Sumw2();
  fQAList->Add(fCellsPi0KidE13);

  fCellsPi0KidM021st  =  new TH2F("fCellsPi0KidM021st","Distribution of cluster energy in its cells (0.10<#lambda_{0}^{2}<0.15, all energies;col;row",19,-9,9,19,-9,9);
  fCellsPi0KidM021st->Sumw2();
  fQAList->Add(fCellsPi0KidM021st);

  fCellsPi0KidM022nd  =  new TH2F("fCellsPi0KidM022nd","Distribution of cluster energy in its cells (0.15<#lambda_{0}^{2}<0.20, all energies;col;row",19,-9,9,19,-9,9);
  fCellsPi0KidM022nd->Sumw2();
  fQAList->Add(fCellsPi0KidM022nd);

  fCellsPi0KidM023rd  =  new TH2F("fCellsPi0KidM023rd","Distribution of cluster energy in its cells (0.20<#lambda_{0}^{2}<0.25, all energies;col;row",19,-9,9,19,-9,9);
  fCellsPi0KidM023rd->Sumw2();
  fQAList->Add(fCellsPi0KidM023rd);

  fCellsPi0KidM024th  =  new TH2F("fCellsPi0KidM024th","Distribution of cluster energy in its cells (0.25<#lambda_{0}^{2}<0.30, all energies;col;row",19,-9,9,19,-9,9);
  fCellsPi0KidM024th->Sumw2();
  fQAList->Add(fCellsPi0KidM024th);


  PostData(1, fOutputList);
  PostData(2, fQAList);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::UserExec(Option_t *) 
{
  // Main loop, called for each event.
  fESDClusters = 0;
  fESDCells = 0;
  fAODClusters = 0;
  fAODCells = 0;
  // event trigger selection
  Bool_t isSelected = 0;
  if(fPeriod.Contains("11a")){
    if(fTrigBit.Contains("kEMC"))
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC1);
    else
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
  }
  else{
    if(fTrigBit.Contains("kEMC"))
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kEMC7);
    else
      if(fTrigBit.Contains("kMB"))
	isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
      else
	isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
  }
  if(fPeriod.Contains("11h")){
    if(fTrigBit.Contains("kAny"))
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAny);
    if(fTrigBit.Contains("kAnyINT"))
      isSelected =  (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kAnyINT);
  }
  if(fIsMc)
    isSelected = kTRUE;
  if(fDebug)
    printf("isSelected = %d, fIsMC=%d\n", isSelected, fIsMc);
  if(!isSelected )
        return; 
  if(fIsMc){
    TTree *tree = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetTree();
    TFile *file = (TFile*)tree->GetCurrentFile();
    TString filename = file->GetName();
    if(!filename.Contains(fPathStrOpt.Data()))
      return;
  }
  fVEvent = (AliVEvent*)InputEvent();
  if (!fVEvent) {
    printf("ERROR: event not available\n");
    return;
  }
  Int_t   runnumber = InputEvent()->GetRunNumber() ;
  if(fDebug)
    printf("run number = %d\n",runnumber);

  fESD = dynamic_cast<AliESDEvent*>(fVEvent);
  if(!fESD){
    fAOD = dynamic_cast<AliAODEvent*>(fVEvent);
    if(!fAOD){
      printf("ERROR: Invalid type of event!!!\n");
      return;
    }
    else if(fDebug)
      printf("AOD EVENT!\n");
  }
  
  fEvtSel->Fill(0);
  if(fDebug)
    printf("event is ok,\n run number=%d\n",runnumber);

  
  AliVVertex *pv = (AliVVertex*)fVEvent->GetPrimaryVertex();
  Bool_t pvStatus = kTRUE;
  if(fESD){
    AliESDVertex *esdv = (AliESDVertex*)fESD->GetPrimaryVertex();
    pvStatus = esdv->GetStatus();
  }
  /*if(fAOD){
    AliAODVertex *aodv = (AliAODVertex*)fAOD->GetPrimaryVertex();
    pvStatus = aodv->GetStatus();
    }*/
  if(!pv)
    return;
  if(!pvStatus)
    fRecoPV->Fill(0);
  else
    fRecoPV->Fill(1);
  fPVtxZ->Fill(pv->GetZ());
  fVecPv.SetXYZ(pv->GetX(),pv->GetY(),pv->GetZ());
  if(TMath::Abs(pv->GetZ())>10)
    return;
  if(fDebug)
    printf("passed vertex cut\n");

  fEvtSel->Fill(1);
  if(fVEvent->IsPileupFromSPD(3, 0.8, 3., 2., 5.) && fPileUpRejSPD){
    if(fDebug)
      printf("Event is SPD pile-up!***\n");
    return;
  }
  if(fESD)
    fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("Tracks"));
  if(fAOD)
    fTracks = dynamic_cast<TClonesArray*>(fAOD->GetTracks());

  if(!fTracks){
    AliError("Track array in event is NULL!");
    if(fDebug)
      printf("returning due to not finding tracks!\n");
    return;
  }
  // Track loop to fill a pT spectrum
  const Int_t Ntracks = fTracks->GetEntriesFast();
  for (Int_t iTracks = 0;  iTracks < Ntracks; ++iTracks) {
    //  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    //AliESDtrack* track = (AliESDtrack*)fESD->GetTrack(iTracks);
    AliVTrack *track = (AliVTrack*)fTracks->At(iTracks);
    if (!track)
      continue;
    AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*>(track);
    AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(track);
    if(esdTrack){
      if(esdTrack->GetEMCALcluster()>0)
	fClusIdFromTracks.Append(Form("%d ",esdTrack->GetEMCALcluster()));
      if (fPrTrCuts && fPrTrCuts->IsSelected(track)){
	fSelPrimTracks->Add(track);
      } else if(fCompTrCuts && fCompTrCuts->IsSelected(track)) {
	fSelPrimTracks->Add(track);
      }
    }
    else if(aodTrack){
      if (fSelHybrid && !aodTrack->IsHybridGlobalConstrainedGlobal())       
	continue ;
      if(!fSelHybrid && !aodTrack->TestFilterBit(fFilterBit))
	continue;
      fSelPrimTracks->Add(track);
      /*if(fTrackMaxPt<track->Pt())
	fTrackMaxPt = track->Pt();*/
    }
  }

  TObjArray *matEMCAL=(TObjArray*)fOADBContainer->GetObject(runnumber,"EmcalMatrices");
  for(Int_t mod=0; mod < (fGeom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
    if(fGeoName=="EMCAL_FIRSTYEARV1" && mod>3)
      break;
    /*if(fESD)
      fGeom->SetMisalMatrix(fESD->GetEMCALMatrix(mod), mod);
      else*/
    // if(fVEvent->GetEMCALMatrix(mod))
    fGeomMatrix[mod] = (TGeoHMatrix*) matEMCAL->At(mod);
    fGeom->SetMisalMatrix(fGeomMatrix[mod] , mod);
  }

  if(fESD){
    AliESDtrackCuts *fTrackCuts = new AliESDtrackCuts();
    fTrackMult = fTrackCuts->GetReferenceMultiplicity(fESD);//kTrackletsITSTPC ,0.5); 
    if(fDebug)
      printf("ESD Track mult= %d\n",fTrackMult);
  }
  else if(fAOD){
    fTrackMult = Ntracks;
    if(fDebug)
      printf("AOD Track mult= %d\n",fTrackMult);
  }
  fTrMultDist->Fill(fTrackMult);
  TList *l = 0;
  TString clusArrayName = "";
  if(fESD){
    l = fESD->GetList();
    if(fDebug)
      l->Print();
    for(int nk=0;nk<l->GetEntries();nk++){
      TObject *obj = (TObject*)l->At(nk);
      TString oname = obj->GetName();
      if(oname.Contains("CaloClus"))
	clusArrayName = oname;
      else
	continue;
      if(clusArrayName=="CaloClusters")
	fClusArrayNames->Fill(0);
      else{
	if(clusArrayName=="EmcCaloClusters")
	  fClusArrayNames->Fill(1);
	else
	  fClusArrayNames->Fill(2);
      }
    }
    fESDClusters =  dynamic_cast<TClonesArray*>(l->FindObject(clusArrayName));
    fESDCells = fESD->GetEMCALCells();
    if(fDebug)
      printf("ESD cluster mult= %d\n",fESDClusters->GetEntriesFast());
  }
  else if(fAOD){
    l = fAOD->GetList();
    if(fDebug)
      l->Print();
    //fAODClusters = dynamic_cast<TClonesArray*>(fAOD->GetCaloClusters());
    for(int nk=0;nk<l->GetEntries();nk++){
      TObject *obj = (TObject*)l->At(nk);
      TString oname = obj->GetName();
      if(oname.Contains("aloClus"))
	clusArrayName = oname;
      else
	continue;
      if(clusArrayName=="caloClusters")
	fClusArrayNames->Fill(0);
      else{
	if(clusArrayName=="EmcCaloClusters")
	  fClusArrayNames->Fill(1);
	else
	  fClusArrayNames->Fill(2);
      }
    }
    fAODClusters = dynamic_cast<TClonesArray*>(l->FindObject(clusArrayName));
    fAODCells = fAOD->GetEMCALCells();
    if(fDebug)
      printf("AOD cluster mult= %d\n",fAODClusters->GetEntriesFast());
  }
  if(fDebug){
    printf("clus array is named %s +++++++++\n",clusArrayName.Data());
  }
  
  
  fMCEvent = MCEvent();
  if(fMCEvent){
    if(fDebug)
      std::cout<<"MCevent exists\n";
    fStack = (AliStack*)fMCEvent->Stack();
    if(!fStack)
      fAODMCParticles = (TClonesArray*)fVEvent->FindListObject("mcparticles");  
  }
  else{
    if(fDebug)
      std::cout<<"ERROR: NO MC EVENT!!!!!!\n";
  }
  fVCells = GetVCaloCells();
  LoopOnCells();
  FollowGamma();
  CheckTriggerPatch();
  FillInvMass();
  if(fDebug)
    printf("passed calling of FollowGamma\n");
  FillClusHists(); 
  if(fDebug)
    printf("passed calling of FillClusHists\n");
  FillMcHists();
  if(fDebug)
    printf("passed calling of FillMcHists\n");
  if(fFillQA)
    FillQA();
  if(fDebug)
    printf("passed calling of FillQA\n");
  if(fESD)
    fESDClusters->Clear();
  fSelPrimTracks->Clear();
  fNClusForDirPho = 0;
  fNCells50 = 0;
  fClusIdFromTracks = "";
  fVecPv.Clear();

  PostData(1, fOutputList);
  PostData(2, fQAList);
}      

//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::FillClusHists()
{
  if(fDebug)
    printf("Inside FillClusHists()....\n");
  // Fill cluster histograms.
  TObjArray *clusters = fESDClusters;

  if (!clusters){
    clusters = fAODClusters;
    if(fDebug)
      printf("ESD clusters empty...");
  }
  if (!clusters){
    if(fDebug)
      printf("  and AOD clusters as well!!!\n"); 
    return;
  }
  if(fDebug)
    printf("\n");

  const Int_t nclus = clusters->GetEntries();
  if(nclus==0)
    return;
  if(fDebug)
    printf("Inside FillClusHists and there are %d clusters\n",nclus);
  Double_t maxE = 0;
  Int_t nclus10 = 0;
  Double_t ptmc=-1;
  for(Int_t ic=0;ic<nclus;ic++){
    maxE=0;
    AliVCluster *c = static_cast<AliVCluster*>(clusters->At(ic));
    if(!c){
      if(fDebug)
	printf("cluster pointer does not exist! xxxx\n");
      continue;
    }
    if(!c->IsEMCAL()){
      if(fDebug)
	printf("cluster is not EMCAL! xxxx\n");
      continue;
    }
    if(c->E()<fECut){
      if(fDebug)
	printf("cluster has E<%1.1f! xxxx\n", fECut);
      continue;
    }
    if(fCpvFromTrack && fClusIdFromTracks.Contains(Form("%d",ic))){
      if(fDebug)
	printf("cluster does not pass CPV criterion! xxxx\n");
       continue;
    }
    if(IsExotic(c)){
      if(fDebug)
	printf("cluster is exotic! xxxx\n");
      continue;
    }
    if(c->GetDistanceToBadChannel()<fDistToBadChan){
      if(fDebug)
	printf("cluster distance to bad channel is %1.1f (<%1.1f) xxxx\n",c->GetDistanceToBadChannel(),fDistToBadChan);
      continue;
    }
    Short_t id;
    Double_t Emax = GetMaxCellEnergy( c, id);
    if(fDebug)
      printf("cluster max cell E=%1.1f",Emax);
    Float_t clsPos[3] = {0,0,0};
    c->GetPosition(clsPos);
    TVector3 clsVec(clsPos);
    clsVec -= fVecPv;
    Double_t Et = c->E()*TMath::Sin(clsVec.Theta());
    if(fDebug)
      printf("\tcluster eta=%1.1f,phi=%1.1f,E=%1.1f\n",clsVec.Eta(),clsVec.Phi(),c->E());
    if(Et>10)
      nclus10++;
    Float_t ceiso=0, cephiband=0, cecore=0;
    Float_t triso=0, trphiband=0, trcore=0;
    Float_t alliso=0, allphiband=0;//, allcore;
    Float_t phibandArea = (1.4 - 2*fIsoConeR)*2*fIsoConeR;
    Float_t netConeArea = TMath::Pi()*(fIsoConeR*fIsoConeR - 0.04*0.04);
    Bool_t isCPV = kFALSE;
    if(TMath::Abs(c->GetTrackDx())>0.03 || TMath::Abs(c->GetTrackDz())>0.02)
      isCPV = kTRUE;
    GetCeIso(clsVec, id, ceiso, cephiband, cecore, Et);
    GetTrIso(clsVec, triso, trphiband, trcore);
    Int_t nInConePairs = 0;
    Double_t onePairMass = 0;
    //---
    //if(c->GetM02()>0.1 && c->GetM02()<0.3 && isCPV){
    TObjArray *inConeInvMassArr = (TObjArray*)fInConeInvMass.Tokenize(";");
    TObjArray *inConePairClEt =  (TObjArray*)fInConePairClEt.Tokenize(";");
    nInConePairs = inConeInvMassArr->GetEntriesFast();
    Int_t nInConePi0 = inConePairClEt->GetEntriesFast();
    Double_t pairEt=0;
    if(nInConePairs != nInConePi0)
      printf("Inconsistent number of in cone pairs!!!\n");
    for(int ipair=0;ipair<nInConePairs;ipair++){
      TObjString *obs = (TObjString*)inConeInvMassArr->At(ipair);
      TObjString *obet = (TObjString*)inConePairClEt->At(ipair);
      TString smass = obs->GetString();
      TString spairEt = obet->GetString();
      Double_t pairmass = smass.Atof();
      pairEt = spairEt.Atof();//this must be zero when inv mass outside pi0 range
      if(0==ipair && nInConePairs==1)
	onePairMass = pairmass;
      if(fDebug)
	printf("=================+++++++++++++++Inv mass inside the cone for photon range: %1.1f,%1.1f,%1.1f+-++++-+-+-+-++-+-+-\n",Et,pairmass,ceiso+triso);
      fEtCandIsoAndIsoWoPairEt->Fill(Et,ceiso+triso,ceiso+triso-pairEt);
    }
    //}
    //---
    Double_t dr = TMath::Sqrt(c->GetTrackDx()*c->GetTrackDx() + c->GetTrackDz()*c->GetTrackDz());
    if(Et>10 && Et<15 && dr>0.025){
      fHigherPtConeM02->Fill(fHigherPtCone,c->GetM02());
      if(fDebug)
	printf("\t\tHigher track pt inside the cone: %1.1f GeV/c; M02=%1.2f\n",fHigherPtCone,c->GetM02());
    }
    alliso = ceiso + triso;
    allphiband = cephiband + trphiband;
    //allcore = cecore + trcore;
    Float_t ceisoue =  cephiband/phibandArea*netConeArea;
    Float_t trisoue =  trphiband/phibandArea*netConeArea;
    Float_t allisoue =  allphiband/phibandArea*netConeArea;
    Float_t mcptsum = GetMcPtSumInCone(clsVec.Eta(), clsVec.Phi(),fIsoConeR); 
    if(fDebug && Et>10)
      printf("\t alliso=%1.1f, Et=%1.1f=-=-=-=-=\n",alliso,Et);
    if(c->GetM02()>0.5 && c->GetM02()<2.0){
      fMcPtInConeBG->Fill(alliso-allisoue, mcptsum);
      fMcPtInConeBGnoUE->Fill(alliso, mcptsum);
      fMcPtInConeTrBGnoUE->Fill(triso, mcptsum);
    }
    if(c->GetM02()>0.1 && c->GetM02()<0.3 && dr>0.03){
      fMcPtInConeSBG->Fill(alliso-allisoue, mcptsum);
      fMcPtInConeSBGnoUE->Fill(alliso, mcptsum);
      fMcPtInConeTrSBGnoUE->Fill(triso, mcptsum);
      if(fMcIdFamily.Contains((Form("%d",c->GetLabel())))){
	fAllIsoEtMcGamma->Fill(Et, alliso-cecore-allisoue);
	fAllIsoNoUeEtMcGamma->Fill(Et, alliso-cecore);
      }
    }
    if(c->GetM02()>0.1 && c->GetM02()<0.3 && isCPV)
      fClusEtCPVSBGISO->Fill(Et,alliso - trcore);
    if(c->GetM02()>0.5 && c->GetM02()<2.0 && isCPV)
      fClusEtCPVBGISO->Fill(Et,alliso - trcore);
    Int_t clusNLM = GetNumberOfLocalMaxima(c,fVCells);
    const Int_t ndims =   fNDimensions;
    Double_t outputValues[ndims];
    if(mcptsum<2)
      ptmc = GetClusSource(c);
    else
      ptmc = -0.1;
    outputValues[0] = Et;
    outputValues[1] = SmearM02(c->GetM02());
    outputValues[2] = ceiso/*cecore*/-ceisoue;
    outputValues[3] = triso-trisoue;
    outputValues[4] = alliso/*cecore*/-allisoue - trcore;
    outputValues[5] = ceiso;
    outputValues[6] = alliso - trcore;
    if(fDebug)
      printf("track-cluster dphi=%1.3f, deta=%1.3f\n",c->GetTrackDx(),c->GetTrackDz());
    if(TMath::Abs(c->GetTrackDx())<0.1)
      outputValues[7] = c->GetTrackDx();
    else
      outputValues[7] = 0.099*c->GetTrackDx()/TMath::Abs(c->GetTrackDx());
    if(TMath::Abs(c->GetTrackDz())<0.05)
      outputValues[8] = c->GetTrackDz();
    else
      outputValues[8] = 0.049*c->GetTrackDz()/TMath::Abs(c->GetTrackDz());
    outputValues[9] = clsVec.Eta();
    outputValues[10] = clsVec.Phi();
    if(fESDCells)
      outputValues[11] = fESDCells->GetCellTime(id);
    else if(fAODCells)
      outputValues[11] = fAODCells->GetCellTime(id);
    outputValues[12] = fTrackMult;
    outputValues[13] = ptmc;
    if(nInConePairs == 1)
      outputValues[14] = onePairMass;
    else
      outputValues[14] = -1;
    if(clusNLM<10)
      outputValues[15] = clusNLM;
    else 
      outputValues[15] = 10;
    fHnOutput->Fill(outputValues);
    if(c->E()>maxE)
      maxE = c->E();
  }
  fNClusHighClusE->Fill(maxE,nclus);
  fMaxEClus = maxE;
  fNClusEt10->Fill(nclus10);
  fNClusPerPho->Fill(fDirPhoPt,fNClusForDirPho);
} 

//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::GetCeIso(TVector3 vec, Int_t maxid, Float_t &iso, Float_t &phiband, Float_t &core, Double_t EtCl)
{
  if(fDebug)
    printf("....indside GetCeIso funtcion\n");
  // Get cell isolation.
  AliVCaloCells *cells = fVCells;
  if (!cells){
     if(fDebug)
      printf("cells are empty  as well!!!\n"); 
    return;
  }

  TObjArray *clusters = fESDClusters;
  if (!clusters)
    clusters = fAODClusters;
  if (!clusters)
    return;
  Double_t pi0mean =  0.1375; 
  Double_t pi0sig = 0.0275;
  Double_t etamean = 0.55;
  Double_t etasig = 0.015;
  Double_t lowpi0mass = pi0mean - pi0sig*fNSigNeutMesonCut;
  Double_t highpi0mass = pi0mean + pi0sig*fNSigNeutMesonCut;
  
  fInConeInvMass = "";
  fInConePairClEt="";
  const Int_t nclus = clusters->GetEntries();
  //const Int_t ncells = cells->GetNumberOfCells();
  Float_t totiso=0;
  Float_t totband=0;
  Float_t totcore=0;
  Float_t etacl = vec.Eta();
  Float_t phicl = vec.Phi();
  if(phicl<0)
    phicl+=TMath::TwoPi();
  /*Int_t absid = -1;
  Float_t eta=-1, phi=-1;  
  for(int icell=0;icell<ncells;icell++){
    absid = TMath::Abs(cells->GetCellNumber(icell));
    Float_t celltime = cells->GetCellTime(absid);
    //if(TMath::Abs(celltime)>2e-8 && (!fIsMc))
    if(TMath::Abs(celltime-maxtcl)>2e-8 )
      continue;
    if(!fGeom)
      return;
    fGeom->EtaPhiFromIndex(absid,eta,phi);
    Float_t dphi = TMath::Abs(phi-phicl);
    Float_t deta = TMath::Abs(eta-etacl);
    Float_t R = TMath::Sqrt(deta*deta + dphi*dphi);*/
  for(int ic=0;ic<nclus;ic++){
    AliVCluster *c = static_cast<AliVCluster*>(clusters->At(ic));
    if(!c)
      continue;
    if(!c->IsEMCAL())
      continue;
    if(c->E()<fMinIsoClusE)
      continue;
    Short_t id=-1;
    GetMaxCellEnergy( c, id);
    Double_t maxct = cells->GetCellTime(id);
    if(TMath::Abs(maxct)>fClusTDiff/*2.5e-9*/ && (!fIsMc))
      continue;
    Float_t clsPos[3] = {0,0,0};
    c->GetPosition(clsPos);
    TVector3 cv(clsPos);
    cv -= fVecPv;
    Double_t Et = c->E()*TMath::Sin(cv.Theta());
    Float_t dphi = (cv.Phi()-phicl);
    Float_t deta = (cv.Eta()-etacl);
    Float_t R = TMath::Sqrt(deta*deta + dphi*dphi);
    if(R<0.007)
      continue;
    if(maxid==id)
      continue;
    Double_t matchedpt =  GetTrackMatchedPt(c->GetTrackMatchedIndex());
    if(fCpvFromTrack){
      if(matchedpt>0 && fRemMatchClus )
	continue;
    } else {
      if(TMath::Abs(c->GetTrackDx())<0.03 && TMath::Abs(c->GetTrackDz())<0.02){
	if(fRemMatchClus){
	  if(fDebug)
	    printf("This isolation cluster is matched to a track!++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	  continue;
	}
      }
    }
    Double_t nEt = TMath::Max(Et-matchedpt, 0.0);
    if(nEt<0){
      printf("nEt=%1.1f\n",nEt);
      continue;
    }
    if(R<fIsoConeR){
      if(c->GetM02()>0.1 && c->GetM02()<0.3 && !(matchedpt>0)){
	TLorentzVector lv, lvec;
	lv.SetPtEtaPhiM(Et,cv.Eta(),cv.Phi(),0);
	lvec.SetPtEtaPhiM(EtCl,vec.Eta(),vec.Phi(),0);
	TLorentzVector lpair = lv + lvec;
	fInConeInvMass += Form("%f;",lpair.M());
	if(lpair.M()>lowpi0mass && lpair.M()<highpi0mass){
	  fInConePairedClusEtVsCandEt->Fill(EtCl,Et);
	  fInConePairClEt += Form("%f;",Et);
	  //continue;
	}
	else 
	  fInConePairClEt += Form("%f;",0.0);
	/*if(lpair.M()>0.52 && lpair.M()<0.58)
	  continue;*/
      }
      totiso += nEt;
      if(R<0.04)
	totcore += nEt;
    }
    else{
      if(dphi>fIsoConeR)
	continue;
      if(deta<fIsoConeR)
	continue;
      totband += nEt;
    }
  }
  iso = totiso;
  phiband = totband;
  core = totcore;
} 
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::GetTrIso(TVector3 vec, Float_t &iso, Float_t &phiband, Float_t &core)
{
  // Get track isolation.

  if(!fSelPrimTracks)
    return;
  fHigherPtCone = 0;
  const Int_t ntracks = fSelPrimTracks->GetEntries();
  Double_t totiso=0;
  Double_t totband=0;
  Double_t totcore=0;
  Double_t etacl = vec.Eta();
  Double_t phicl = vec.Phi();
  if(phicl<0)
    phicl+=TMath::TwoPi();
  for(int itrack=0;itrack<ntracks;itrack++){
    AliVTrack *track = static_cast<AliVTrack*> (fSelPrimTracks->At(itrack));
    if(!track)
      continue;
    Double_t dphi = TMath::Abs(track->Phi()-phicl);
    Double_t deta = TMath::Abs(track->Eta()-etacl);
    Double_t R = TMath::Sqrt(deta*deta + dphi*dphi);
    Double_t pt = track->Pt();
    if(pt>fHigherPtCone)
      fHigherPtCone = pt;
    if(R<fIsoConeR){
      totiso += track->Pt();
      if(R<0.04 && this->fTrCoreRem)
	totcore += pt;
    }
    else{
      if(dphi>fIsoConeR)
	continue;
      if(deta<fIsoConeR)
	continue;
      totband += track->Pt();
    }
  }
  iso = totiso;
  phiband = totband;
  core = totcore;
} 

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALIsoPhoton::GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax)
{
  // Calculate the energy of cross cells around the leading cell.

  AliVCaloCells *cells = fVCells;
  if (!cells)
    return 0;

  if (!fGeom)
    return 0;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Int_t iphis   = -1;
  Int_t ietas   = -1;

  Double_t crossEnergy = 0;

  fGeom->GetCellIndex(idmax,iSupMod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphis,ietas);

  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Int_t cellAbsId = cluster->GetCellAbsId(i);
    fGeom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
    Int_t aphidiff = TMath::Abs(iphi-iphis);
    if (aphidiff>1)
      continue;
    Int_t aetadiff = TMath::Abs(ieta-ietas);
    if (aetadiff>1)
      continue;
    if ( (aphidiff==1 && aetadiff==0) ||
	(aphidiff==0 && aetadiff==1) ) {
      crossEnergy += cells->GetCellAmplitude(cellAbsId);
    }
  }

  return crossEnergy;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALIsoPhoton ::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
{
  // Get maximum energy of attached cell.

  id = -1;

  AliVCaloCells *cells = fVCells;
  if(!cells)
    return 0;

  Double_t maxe = 0;
  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Double_t e = cells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
    if (e>maxe) {
      maxe = e;
      id   = cluster->GetCellAbsId(i);
    }
  }
  return maxe;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton ::FillMcHists()
{
  if(!fStack && !fAODMCParticles)
    return;
  //ESD
  if(fStack){
    Int_t nTracks = fStack->GetNtrack();
    if(fDebug)
      printf("Inside FillMcHists and there are %d mcparts\n",nTracks);
    for(Int_t iTrack=0;iTrack<nTracks;iTrack++){
      TParticle *mcp = static_cast<TParticle*>(fStack->Particle(iTrack));  
      if(!mcp)
	continue;  
      Int_t pdg = mcp->GetPdgCode();
      if(pdg!=22)
	continue;
      if(TMath::Abs(mcp->Eta())>0.7 ||mcp->Phi()<1.4 || mcp->Phi()>3.2)
	continue;
      Int_t imom = mcp->GetMother(0);
      if(imom<0 || imom>nTracks)
	continue;
      TParticle *mcmom = static_cast<TParticle*>(fStack->Particle(imom));  
      if(!mcmom)
	continue;
      Int_t pdgMom = mcmom->GetPdgCode();
      Double_t mcphi = mcp->Phi();
      Double_t mceta = mcp->Eta();
      if((imom==6 || imom==7) && pdgMom==22) {
	fMCDirPhotonPtEtaPhi->Fill(mcp->Pt(),mcp->Eta(),mcp->Phi());
	Float_t ptsum = GetMcPtSumInCone(mcp->Eta(), mcp->Phi(), fIsoConeR);
	fMcPtInConeMcPhoPt->Fill(mcp->Pt(),ptsum);
	if(ptsum<2)
	  fMCIsoDirPhotonPtEtaPhi->Fill(mcp->Pt(),mcp->Eta(),mcp->Phi());
	if(mcphi<(3.14-fIsoConeR) && mcphi>(1.4+fIsoConeR) && TMath::Abs(mceta)<(0.7-fIsoConeR))
	  fMCDirPhotonPtEtIso->Fill(mcp->Pt(),ptsum);
	if(fNClusForDirPho==0)
	  fMCDirPhotonPtEtaPhiNoClus->Fill(mcp->Pt(),mcp->Eta(),mcp->Phi());
	if(fDebug){
	  printf("Found \"photonic\" parton at position %d, with pt=%1.1f, eta=%1.1f and phi=%1.1f, and status=%d,\n",imom,mcmom->Pt(), mcmom->Eta(), mcmom->Phi(), mcmom->GetStatusCode());
	  printf("with a final photon at position %d, with pt=%1.1f, eta=%1.1f and phi=%1.1f, and status=%d\n",iTrack,mcp->Pt(), mcp->Eta(), mcp->Phi(),mcp->GetStatusCode());
	}
      }
      else{
	if(TMath::Abs(pdgMom)>100 && TMath::Abs(pdgMom)<1000)
	  fDecayPhotonPtMC->Fill(mcp->Pt());
      }
    }
  }
  //AOD 
  else if(fAODMCParticles){
    Int_t nTracks = fAODMCParticles->GetEntriesFast();
    if(fDebug)
      printf("Inside FillMcHists and there are %d mcparts\n",nTracks);
    for(Int_t iTrack=0;iTrack<nTracks;iTrack++){
      AliAODMCParticle *mcp = dynamic_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));
      if(!mcp)
	continue;
      Int_t pdg = mcp->GetPdgCode();
      if(pdg!=22)
	continue;
      if(TMath::Abs(mcp->Eta())>0.7 ||mcp->Phi()<1.4 || mcp->Phi()>3.2)
	continue;
      Int_t imom = mcp->GetMother();
      if(imom<0 || imom>nTracks)
	continue;
      AliAODMCParticle *mcmom = static_cast<AliAODMCParticle*>(fAODMCParticles->At(imom));
      if(!mcmom)
	continue;
      Int_t pdgMom = mcmom->GetPdgCode();
      Double_t mcphi = mcp->Phi();
      Double_t mceta = mcp->Eta();
      if((imom==6 || imom==7) && pdgMom==22) {
	fMCDirPhotonPtEtaPhi->Fill(mcp->Pt(),mcp->Eta(),mcp->Phi());
	Float_t ptsum = GetMcPtSumInCone(mcp->Eta(), mcp->Phi(), fIsoConeR);
	fMcPtInConeMcPhoPt->Fill(mcp->Pt(),ptsum);
	if(ptsum<2)
	  fMCIsoDirPhotonPtEtaPhi->Fill(mcp->Pt(),mcp->Eta(),mcp->Phi());
	if(mcphi<(3.14-fIsoConeR) && mcphi>(1.4+fIsoConeR) && TMath::Abs(mceta)<(0.7-fIsoConeR))
	  fMCDirPhotonPtEtIso->Fill(mcp->Pt(),ptsum);
	if(fNClusForDirPho==0)
	  fMCDirPhotonPtEtaPhiNoClus->Fill(mcp->Pt(),mcp->Eta(),mcp->Phi());
	if(fDebug){
	  printf("Found \"photonic\" parton at position %d, with pt=%1.1f, eta=%1.1f and phi=%1.1f, and status=%d,\n",imom,mcmom->Pt(), mcmom->Eta(), mcmom->Phi(), mcmom->GetStatus());
	  printf("with a final photon at position %d, with pt=%1.1f, eta=%1.1f and phi=%1.1f, and status=%d\n",iTrack,mcp->Pt(), mcp->Eta(), mcp->Phi(),mcp->GetStatus());
	}
      }
      else{
	if(TMath::Abs(pdgMom)>100 && TMath::Abs(pdgMom)<1000)
	  fDecayPhotonPtMC->Fill(mcp->Pt());
      }
    }
  }
}
//________________________________________________________________________
Float_t AliAnalysisTaskEMCALIsoPhoton::GetClusSource(const AliVCluster *c)
{
  if(!c)
    return -0.1;
  if(!fStack && !fAODMCParticles)
    return -0.1;
  Int_t clabel = c->GetLabel();
  if(fDebug && fMcIdFamily.Contains(Form("%d",clabel)))
    printf("\n\t ==== Label %d is a descendent of the prompt photon ====\n\n",clabel);
  if(!fMcIdFamily.Contains(Form("%d",clabel)))
    return -0.1;
  fNClusForDirPho++;
  TString partonposstr = (TSubString)fMcIdFamily.operator()(0,1);
  Int_t partonpos = partonposstr.Atoi();
  if(fDebug)
    printf("\t^^^^ parton position = %d ^^^^\n",partonpos);
  Float_t clsPos[3] = {0,0,0};
  c->GetPosition(clsPos);
  TVector3 clsVec(clsPos);
  clsVec -= fVecPv;
  Double_t Et = c->E()*TMath::Sin(clsVec.Theta());
  //ESD
  if(fStack){
    Int_t nmcp = fStack->GetNtrack();
    if(clabel<0 || clabel>nmcp)
      return -0.1;
    TParticle *mcp = static_cast<TParticle*>(fStack->Particle(partonpos));
    if(!mcp)
      return -0.1;
    if(fDebug){
      printf("\tclus mc truth eta=%1.1f,phi=%1.1f,E=%1.1f, pdgcode=%d, stackpos=%d\n",mcp->Eta(),mcp->Phi(),mcp->Energy(),mcp->GetPdgCode(),clabel);
    }
    Int_t lab1 =  mcp->GetFirstDaughter();
    if(lab1<0 || lab1>nmcp)
      return -0.1;
    TParticle *mcd = static_cast<TParticle*>(fStack->Particle(lab1));
    if(!mcd)
      return -0.1;
    if(fDebug)
      printf("\t\tmom mc truth eta=%1.1f,phi=%1.1f,E=%1.1f, pdgcode=%d, stackpos=%d\n",mcd->Eta(),mcd->Phi(),mcd->Energy(),mcd->GetPdgCode(),lab1);
    if(mcd->GetPdgCode()==22){
      fClusEtMcPt->Fill(mcd->Pt(), Et);
      fClusMcDetaDphi->Fill(clsVec.Eta() - mcd->Eta(), clsVec.Phi() - mcd->Phi());
    }
    else{
      if(fDebug)
	printf("Warning: daughter of photon parton is not a photon\n");
      return -0.1;
    }
  }
  //AOD
  else if(fAODMCParticles){
    Int_t nmcp = fAODMCParticles->GetEntriesFast();
    if(clabel<0 || clabel>nmcp)
      return -0.1;
    AliAODMCParticle *mcp = static_cast<AliAODMCParticle*>(fAODMCParticles->At(partonpos));
    if(!mcp)
      return -0.1;
    if(fDebug){
      printf("\tclus mc truth eta=%1.1f,phi=%1.1f,E=%1.1f, pdgcode=%d, stackpos=%d\n",mcp->Eta(),mcp->Phi(),mcp->E(),mcp->GetPdgCode(),clabel);
    }
    Int_t lab1 =  mcp->GetDaughterLabel(0);
    if(lab1<0 || lab1>nmcp)
      return -0.1;
    AliAODMCParticle *mcd = static_cast<AliAODMCParticle*>(fAODMCParticles->At(lab1));
    if(!mcd)
      return -0.1;
    if(fDebug)
      printf("\t\tmom mc truth eta=%1.1f,phi=%1.1f,E=%1.1f, pdgcode=%d, stackpos=%d\n",mcd->Eta(),mcd->Phi(),mcd->E(),mcd->GetPdgCode(),lab1);
    if(mcd->GetPdgCode()==22){
      fClusEtMcPt->Fill(mcd->Pt(), Et);
      fClusMcDetaDphi->Fill(clsVec.Eta() - mcd->Eta(), clsVec.Phi() - mcd->Phi());
    }
    else{
      printf("Warning: daughter of photon parton is not a photon\n");
      return -0.1;
    }
  }
  return fDirPhoPt;
}
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::FollowGamma()
{
  if(!fStack && !fAODMCParticles)
    return;
  Int_t selfid = 6;
  Int_t daug0f=-1, daug0l=-1, nd0=0;
  //ESD
  if(fStack){  
    TParticle *mcp = static_cast<TParticle*>(fStack->Particle(selfid));
    if(!mcp)
      return;
    if(mcp->GetPdgCode()!=22){
      mcp = static_cast<TParticle*>(fStack->Particle(++selfid));
      if(!mcp)
	return;
    }  
    daug0f =  mcp->GetFirstDaughter();
    daug0l =  mcp->GetLastDaughter();
    nd0 = daug0l - daug0f;
    if(fDebug)
      printf("\n\tGenerated gamma (%d) eta=%1.1f,phi=%1.1f,E=%1.1f, pdgcode=%d, n-daug=%d\n",selfid,mcp->Eta(),mcp->Phi(),mcp->Energy(),mcp->GetPdgCode(),nd0+1);
    fMcIdFamily = Form("%d,",selfid);
    GetDaughtersInfo(daug0f,daug0l, selfid,"");
    if(fDebug){
      printf("\t---- end of the prompt  gamma's family tree ----\n\n");
      printf("Family id string = %s,\n\n",fMcIdFamily.Data());
    }
    TParticle *md = static_cast<TParticle*>(fStack->Particle(daug0f));
    if(!md)
      return;
    fDirPhoPt = md->Pt();
  }
  //AOD
  else   if(fAODMCParticles){  
    AliAODMCParticle *mcp = static_cast<AliAODMCParticle*>(fAODMCParticles->At(selfid));
    if(!mcp)
      return;
    if(mcp->GetPdgCode()!=22){
      mcp = static_cast<AliAODMCParticle*>(fAODMCParticles->At(++selfid));
      if(!mcp)
	return;
    }  
    daug0f =  mcp->GetDaughterLabel(0);
    daug0l =  mcp->GetDaughterLabel(mcp->GetNDaughters()-1);
    nd0 = daug0l - daug0f;
    if(fDebug)
      printf("\n\tGenerated gamma (%d) eta=%1.1f,phi=%1.1f,E=%1.1f, pdgcode=%d, n-daug=%d\n",selfid,mcp->Eta(),mcp->Phi(),mcp->E(),mcp->GetPdgCode(),nd0+1);
    fMcIdFamily = Form("%d,",selfid);
    GetDaughtersInfo(daug0f,daug0l, selfid,"");
    if(fDebug){
      printf("\t---- end of the prompt  gamma's family tree ----\n\n");
      printf("Family id string = %s,\n\n",fMcIdFamily.Data());
    }
    AliAODMCParticle *md = static_cast<AliAODMCParticle*>(fAODMCParticles->At(daug0f));
    if(!md)
      return;
    fDirPhoPt = md->Pt();
  }

}
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::GetDaughtersInfo(int firstd, int lastd, int selfid, const char *inputind)
{
  if(fStack){
    int nmcp = fStack->GetNtrack();
    if(firstd<0 || firstd>nmcp)
      return;
    if(lastd<0 || lastd>nmcp)
      return;
    if(firstd>lastd){
      printf("WARNING: First daughter > last (%d,%d)\n",firstd,lastd);
    return;
    }
    TString indenter = Form("\t%s",inputind);
    TParticle *dp = 0x0;
    if(fDebug)
      printf("\t%s--- Daughters of particle %d ---\n", indenter.Data(), selfid);
    for(int id=firstd; id<lastd+1; id++){
      dp =  static_cast<TParticle*>(fStack->Particle(id));
      if(!dp)
	continue;
      Int_t dfd = dp->GetFirstDaughter(); 
      Int_t dld = dp->GetLastDaughter();
      Int_t dnd =  dld - dfd + 1;
      Float_t vr = TMath::Sqrt(dp->Vx()*dp->Vx()+dp->Vy()*dp->Vy());
      if(fDebug)
	printf("\t%sParticle daughter(%d) eta=%1.1f,phi=%1.1f,E=%1.1f, vR=%1.1f, pdgcode=%d, n-daug=%d(%d,%d)\n", indenter.Data(),id, dp->Eta(), dp->Phi(), dp->Energy(), vr, dp->GetPdgCode(), dnd, dfd, dld);
      fMcIdFamily += Form("%d,",id);
      GetDaughtersInfo(dfd,dld,id,indenter.Data());
    }
  }
  if(fAODMCParticles){
    int nmcp = fAODMCParticles->GetEntriesFast();
    if(firstd<0 || firstd>nmcp)
      return;
    if(lastd<0 || lastd>nmcp)
      return;
    if(firstd>lastd){
      printf("WARNING: First daughter > last (%d,%d)\n",firstd,lastd);
    return;
    }
    TString indenter = Form("\t%s",inputind);
    AliAODMCParticle *dp = 0x0;
    if(fDebug)
      printf("\t%s--- Daughters of particle %d ---\n", indenter.Data(), selfid);
    for(int id=firstd; id<lastd+1; id++){
      dp =  static_cast<AliAODMCParticle*>(fAODMCParticles->At(id));
      if(!dp)
	continue;
      Int_t dfd = dp->GetDaughterLabel(0); 
      Int_t dld = dp->GetDaughterLabel(dp->GetNDaughters()-1);
      Int_t dnd =  dld - dfd + 1;
      Float_t vr = TMath::Sqrt(dp->Xv()*dp->Xv()+dp->Xv()*dp->Xv());
      if(fDebug)
	printf("\t%sParticle daughter(%d) eta=%1.1f,phi=%1.1f,E=%1.1f, vR=%1.1f, pdgcode=%d, n-daug=%d(%d,%d)\n", indenter.Data(),id, dp->Eta(), dp->Phi(), dp->E(), vr, dp->GetPdgCode(), dnd, dfd, dld);
      fMcIdFamily += Form("%d,",id);
      GetDaughtersInfo(dfd,dld,id,indenter.Data());
    }
  }
}

//________________________________________________________________________
Float_t AliAnalysisTaskEMCALIsoPhoton::GetMcPtSumInCone(Float_t etaclus, Float_t phiclus, Float_t R)
{
  if(!fStack && !fAODMCParticles)
    return 0;
  if(fDebug)
    printf("Inside GetMcPtSumInCone!!\n");
  Float_t ptsum = 0;
  TString addpartlabels = "";
  //ESD
  if(fStack){
    Int_t nTracks = fStack->GetNtrack();
    for(Int_t iTrack=9;iTrack<nTracks;iTrack++){
      TParticle *mcp = static_cast<TParticle*>(fStack->Particle(iTrack));  
      if(!mcp)
	continue;  
      Int_t status = mcp->GetStatusCode();
      if(status!=1)
	continue;
      Float_t mcvr = TMath::Sqrt(mcp->Vx()*mcp->Vx()+ mcp->Vy()* mcp->Vy() + mcp->Vz()*mcp->Vz());
      if(mcvr>1)
	continue;
      /*else {
	if(fDebug)
	printf("      >>>> mcp Rho, Vx, Vy, Vz = %1.1f,%1.1f,%1.1f,%1.1f.......\n",mcp->Rho(),mcp->Vx(), mcp->Vy(),mcp->Vz());
	}*/
      Float_t dphi = mcp->Phi() - phiclus;
      Float_t deta = mcp->Eta() - etaclus;
      if(fDebug && TMath::Abs(dphi)<0.01)
	printf("      >>>> mcphi = %1.1f, mceta = %1.1f\n>>>> dphi = %1.1f, deta = %1.1f\n", mcp->Phi(), mcp->Eta(),dphi,deta);
      
      if(deta>R || dphi>R)
	continue;
      Float_t dR = TMath::Sqrt(dphi*dphi +  deta*deta);
      if(dR>R)
	continue;
      addpartlabels += Form("%d,",iTrack);
      if(mcp->Pt()<0.2)
	continue;
      ptsum += mcp->Pt();
    }
  }
  //AOD
  if(fAODMCParticles){
    Int_t nTracks = fAODMCParticles->GetEntriesFast();
    for(Int_t iTrack=9;iTrack<nTracks;iTrack++){
      AliAODMCParticle *mcp = static_cast<AliAODMCParticle*>(fAODMCParticles->At(iTrack));  
      if(!mcp)
	continue;  
      Int_t status = mcp->GetStatus();
      if(status!=1)
	continue;
      Float_t mcvr = TMath::Sqrt(mcp->Xv()*mcp->Xv()+ mcp->Yv()* mcp->Yv() + mcp->Zv()*mcp->Zv());
      if(mcvr>1)
	continue;
      /*else {
	if(fDebug)
	printf("      >>>> mcp Rho, Vx, Vy, Vz = %1.1f,%1.1f,%1.1f,%1.1f.......\n",mcp->Rho(),mcp->Vx(), mcp->Vy(),mcp->Vz());
	}*/
      Float_t dphi = mcp->Phi() - phiclus;
      Float_t deta = mcp->Eta() - etaclus;
      if(fDebug && TMath::Abs(dphi)<0.01)
	printf("      >>>> mcphi = %1.1f, mceta = %1.1f\n>>>> dphi = %1.1f, deta = %1.1f\n", mcp->Phi(), mcp->Eta(),dphi,deta);
      
      if(deta>R || dphi>R)
	continue;
      Float_t dR = TMath::Sqrt(dphi*dphi +  deta*deta);
      if(dR>R)
	continue;
      addpartlabels += Form("%d,",iTrack);
      if(mcp->Pt()<0.2)
	continue;
      ptsum += mcp->Pt();
    }
  }
  return ptsum;
}
//________________________________________________________________________
bool AliAnalysisTaskEMCALIsoPhoton::IsMcPDG(Int_t label, Int_t PDG)
{
  bool foundpi0=false;
  if(!fStack && !fAODMCParticles)
    return false;
  bool prompho = false;
  if(PDG==22 && fMcIdFamily.Contains(Form("%d",label)))
    prompho = true;
  int imother=label+1;
  int nmcp=0;
  //ESD
  if(fStack){
    nmcp = fStack->GetNtrack();
    if(label<0 || label>nmcp)
      return false;
    TParticle *mcp = static_cast<TParticle*>(fStack->Particle(label));  
    if(!mcp)
      return false;
    if((mcp->GetPdgCode()==22) && prompho)
      return prompho;
    if(mcp->GetPdgCode()==TMath::Abs(PDG))
      return true;
    imother = mcp->GetMother(0);
    foundpi0 = IsMcPDG(imother,PDG);
  }
  //AOD
  if(fAODMCParticles){
    nmcp = fAODMCParticles->GetEntriesFast();
    if(label<0 || label>nmcp)
      return false;
    AliAODMCParticle *mcp = static_cast<AliAODMCParticle*>(fAODMCParticles->At(label));
    if(!mcp)
      return false;
    if((mcp->GetPdgCode()==22) && prompho)
      return prompho;
    if(mcp->GetPdgCode()==TMath::Abs(PDG))
      return true;
    imother = mcp->GetMother();
    foundpi0 = IsMcPDG(imother,PDG);
  }
  return foundpi0;
}
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::FillQA() 
{

  TObjArray *clusters = fESDClusters;
  //"none", "exotic", "exo+cpv1", "exo+cpv1+time", "exo+cpv1+time+m02"),
  if (!clusters){
    clusters = fAODClusters;
    if(fDebug)
      printf("ESD clusters empty...");
  }
  if (!clusters){
    if(fDebug)
      printf("  and AOD clusters as well!!!\n"); 
    return;
  }
  if(!fSelPrimTracks)
    return;
  const int ntracks = fSelPrimTracks->GetEntriesFast();
  const int ncells = fNCells50;//fESDCells->GetNumberOfCells();
  const Int_t nclus = clusters->GetEntries();
  fNTracks->Fill(ntracks);
  fEmcNCells->Fill(ncells);
  fEmcNClus->Fill(nclus);
  if(fMaxEClus>fECut){
    fNTracksECut->Fill(ntracks);
    fEmcNCellsCut->Fill(ncells);
    fEmcNClusCut->Fill(nclus);
  }
  for(int it=0;it<ntracks;it++){
    AliVTrack *t = (AliVTrack*)fSelPrimTracks->At(it);
    if(!t)
      continue;
    fTrackPtPhi->Fill(t->Pt(),t->Phi());
    fTrackPtEta->Fill(t->Pt(),t->Eta());
    if(fMaxEClus>fECut){
      fTrackPtPhiCut->Fill(t->Pt(), t->Phi());
      fTrackPtEtaCut->Fill(t->Pt(), t->Eta());
    }
    if(t->GetTPCsignal()<80 || t->GetTPCsignal()>100)
      continue;
    if(t->GetEMCALcluster()<=0 || t->GetEMCALcluster()>nclus)
      continue;
    AliVCluster *c = dynamic_cast<AliVCluster*>(clusters->At(t->GetEMCALcluster()));
    if(!c)
      continue;
    fEoverPvsE->Fill(c->E(),c->E()/t->P());
  }
  for(int ic=0;ic<nclus;ic++){
    AliVCluster *c = dynamic_cast<AliVCluster*>(clusters->At(ic));
    //AliESDCaloCluster *c = (AliESDCaloCluster*)clusters->At(ic);
    if(!c)
      continue;
    if(!c->IsEMCAL())
      continue;
    Float_t clsPos[3] = {0,0,0};
    c->GetPosition(clsPos);
    TVector3 clsVec(clsPos);
    clsVec -= fVecPv;
    Double_t cphi = clsVec.Phi();
    Double_t ceta = clsVec.Eta();
    Short_t id = -1;
    GetMaxCellEnergy( c, id);
    fEmcClusEClusCuts->Fill(c->E(),0);
    fEmcClusEPhi->Fill(c->E(), cphi);
    fEmcClusEEta->Fill(c->E(), ceta);
    if(fMaxEClus>fECut){
      fEmcClusEPhiCut->Fill(c->E(), cphi);
      fEmcClusEEtaCut->Fill(c->E(), ceta);
    }
    Double_t maxt=0;
    if(fESDCells)
      maxt = fESDCells->GetCellTime(id);
    else if(fAODCells)
      maxt = fAODCells->GetCellTime(id);
    if(IsExotic(c))
      continue;
    fEmcClusNotExo->Fill(c->E());
    fEmcClusEClusCuts->Fill(c->E(),1);
    if(fClusIdFromTracks.Contains(Form("%d",ic))){
      fEmcClusETM2->Fill(c->E());
      fDetaDphiFromTM->Fill(c->GetTrackDz(),c->GetTrackDx());
    }
    if(IsPi0M02(c->GetM02(),c->E()*TMath::Sin(clsVec.Theta())))
      fTrackDzDxPhoSS->Fill(c->GetTrackDz(),c->GetTrackDx());
    else{
    if(IsPi0M02(c->GetM02()-0.1,c->E()*TMath::Sin(clsVec.Theta())))
      fTrackDzDxPhoSS_bg->Fill(c->GetTrackDz(),c->GetTrackDx());
    if(IsPi0M02(c->GetM02()+0.1,c->E()*TMath::Sin(clsVec.Theta())))
      fTrackDzDxPhoSS_bg->Fill(c->GetTrackDz(),c->GetTrackDx());
    }
    if(fIsMc){
      bool ispi0 = IsMcPDG(c->GetLabel(),111);
      bool ispich = IsMcPDG(c->GetLabel(),211);
      bool isprompho = IsMcPDG(c->GetLabel(),22);
      if(ispi0)
	fTrackDEtaDPhiPho->Fill(c->GetTrackDz(),c->GetTrackDx());
      if(ispich)
	fTrackDEtaDPhiPi->Fill(c->GetTrackDz(),c->GetTrackDx());
      if(isprompho)
	fTrackDEtaDPhiPPho->Fill(c->GetTrackDz(),c->GetTrackDx());
    }
    if(TMath::Abs(c->GetTrackDx())<0.03 && TMath::Abs(c->GetTrackDz())<0.02){
      fEmcClusETM1->Fill(c->E());
      continue;
    }
    fEmcClusEClusCuts->Fill(c->E(),2);
    if(TMath::Abs(maxt)>30e-9 && !fIsMc)
      continue;
    fEmcClusEClusCuts->Fill(c->E(),3);
    if(c->GetM02()>0.1)
      fEmcClusEClusCuts->Fill(c->E(),4);
  }
}
//________________________________________________________________________
Double_t AliAnalysisTaskEMCALIsoPhoton::GetTrackMatchedPt(Int_t matchIndex)
{
  Double_t pt = 0;
  if(!fTracks)
    return pt;
  if(matchIndex<0 || matchIndex>fTracks->GetEntries()){
    if(fDebug)
      printf("track-matched index out of track array range!!!\n");
    return pt;
  }
  AliVTrack* track = static_cast<AliVTrack*>(fTracks->At(matchIndex));
  if(!track){
    if(fDebug)
      printf("track-matched pointer does not exist!!!\n");
    return pt;
  }
  if(fESD){
    if(fPrTrCuts && fPrTrCuts->IsSelected(track))
      pt = track->Pt();
    else {
      if(fCompTrCuts && fCompTrCuts->IsSelected(track))
	pt = track->Pt();
      else
	return pt;
    }
  }
  return pt;
}
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::LoopOnCells()
{
  AliVCaloCells *cells = fVCells;
  if(!cells)
    return;
  Double_t maxe = 0;
  Double_t maxphi = -10;
  Int_t ncells = cells->GetNumberOfCells();
  Double_t eta,phi;
  for (Int_t i=0; i<ncells; i++) {
    Short_t absid = TMath::Abs(cells->GetCellNumber(i));
    Double_t e = cells->GetCellAmplitude(absid);
    if(e>0.05)
      fNCells50++;
    else 
      continue;
    fGeom->EtaPhiFromIndex(absid,eta,phi);
    if(maxe<e){
      maxe = e;
      maxphi = phi;
    }
  }
  fMaxCellEPhi->Fill(maxe,maxphi);

}
//________________________________________________________________________
bool AliAnalysisTaskEMCALIsoPhoton::IsExotic(AliVCluster *c)
{
  bool isExo = 0;
  Short_t id = -1;
  Double_t Emax = GetMaxCellEnergy( c, id);
  Double_t Ecross = GetCrossEnergy( c, id);
  if((1-Ecross/Emax)>fExoticCut)
    isExo = 1;
  return isExo;
}
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::CheckTriggerPatch()
{
  //printf("inside CheckTriggerPatch...+++===***===+++...\n");
  AliVEvent *ve;
  if(fESD){
    ve = (AliVEvent*)fESD;
  }
  else{
    if(fAOD)
      ve = (AliVEvent*)fAOD;
    else
      return;
  }
  if(!ve)
    return;
  //printf("\tevent ok\n");
  TString trigname = "EmcalTrigger";
  TList *l = (TList*)ve->GetList();
  if(l){
    Int_t nentl = l->GetEntries();
    //printf("There are %d objects in the event list\n",nentl);
    for(Int_t il=0;il<nentl;il++){
      TObject *obj = (TObject*)l->At(il);
      TString oname = Form("%s",obj->GetName());
      //printf("\t object %d named %s\n",il,oname.Data());
      if((oname.Contains("Emc") || oname.Contains("EMC")) && (oname.Contains("rigg")))
	trigname = oname;
    }
  }
  //printf("\t\t trigger name  =  %s\n",trigname.Data()); 
  TClonesArray *triPatchInfo = dynamic_cast<TClonesArray*>(ve->FindListObject(trigname));
  if(!triPatchInfo){
    //printf("no patch info array\n");
    return;
  }
  //printf("made it!!!*****************************\n");
  Int_t nPatch = triPatchInfo->GetEntries();
  if(nPatch>1)
    //printf("more than one calo trigger patch in this event!\n");
  for(Int_t ip = 0;ip<nPatch;ip++){
    AliEMCALTriggerPatchInfo *pti = dynamic_cast<AliEMCALTriggerPatchInfo*>(triPatchInfo->At(ip));
    if(!pti)
      continue;
    //printf("\ttrigger patch E=%1.1f\n",pti->GetPatchE());
    fETrigg->Fill(pti->GetPatchE());
    if(!pti->IsLevel0())
      continue;
  }
}
//________________________________________________________________________
Double_t AliAnalysisTaskEMCALIsoPhoton::SmearM02(Double_t m02)
{
  if(!fStack && !fAODMCParticles)
    return m02;
  TRandom3 *r = new TRandom3(0);
  return m02+r->Gaus(0,fSigmaSmear);
}
//________________________________________________________________________
Int_t AliAnalysisTaskEMCALIsoPhoton::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells) 
{  
  const Int_t   nc = cluster->GetNCells();
  
  Int_t   absIdList[nc]; 
  Float_t maxEList[nc]; 
  
  Int_t nMax = GetNumberOfLocalMaxima(cluster, cells, absIdList, maxEList);
  
  return nMax;
}
//________________________________________________________________________
Int_t AliAnalysisTaskEMCALIsoPhoton::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells,
                                                  Int_t *absIdList,     Float_t *maxEList) 
{
  Int_t iDigitN = 0 ;
  Int_t iDigit  = 0 ;
  Int_t absId1 = -1 ;
  Int_t absId2 = -1 ;
  const Int_t nCells = cluster->GetNCells();
  
  Float_t eCluster = cluster->E();
  Float_t fLocalMaxCutE = 0.1;
  Float_t fLocMaxCutEDiff = 0.05;

  
  if(!cluster->IsEMCAL())
    return 0;
  
  //printf("cluster : ncells %d \n",nCells);
  
  Float_t emax  = 0;
  Int_t   idmax =-1;
  for(iDigit = 0; iDigit < nCells ; iDigit++)
  {
    absIdList[iDigit] = cluster->GetCellsAbsId()[iDigit]  ; 
    Float_t en = cells->GetCellAmplitude(absIdList[iDigit]);
    
    
    if( en > emax )
    {
      emax  = en ;
      idmax = absIdList[iDigit] ;
    }
    //Int_t icol = -1, irow = -1, iRCU = -1;
    //Int_t sm = GetModuleNumberCellIndexes(absIdList[iDigit], calorimeter, icol, irow, iRCU) ;
    //printf("\t cell %d, id %d, sm %d, col %d, row %d, e %f\n", iDigit, absIdList[iDigit], sm, icol, irow, en );
  }
  
  for(iDigit = 0 ; iDigit < nCells; iDigit++) 
  {   
    if( absIdList[iDigit] >= 0 ) 
    {
      absId1 = cluster->GetCellsAbsId()[iDigit];
      
      Float_t en1 = cells->GetCellAmplitude(absId1);
      
      
      //printf("%d : absIDi %d, E %f\n",iDigit, absId1,en1);
      
      for(iDigitN = 0; iDigitN < nCells; iDigitN++) 
      {	
        absId2 = cluster->GetCellsAbsId()[iDigitN] ;
        
        if(absId2==-1 || absId2==absId1) continue;
        
        //printf("\t %d : absIDj %d\n",iDigitN, absId2);
        
        Float_t en2 = cells->GetCellAmplitude(absId2);
        

        //printf("\t %d : absIDj %d, E %f\n",iDigitN, absId2,en2);
        
        if ( AreNeighbours( absId1, absId2) ) 
        {
          // printf("\t \t Neighbours \n");
          if ( en1 > en2 ) 
          {    
            absIdList[iDigitN] = -1 ;
            //printf("\t \t indexN %d not local max\n",iDigitN);
            // but may be digit too is not local max ?
            if(en1 < en2 + fLocMaxCutEDiff) {
              //printf("\t \t index %d not local max cause locMaxCutEDiff\n",iDigit);
              absIdList[iDigit] = -1 ;
            }
          }
          else 
          {
            absIdList[iDigit] = -1 ;
            //printf("\t \t index %d not local max\n",iDigitN);
            // but may be digitN too is not local max ?
            if(en1 > en2 - fLocMaxCutEDiff) 
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

      
      
      if(en < fLocalMaxCutE ) continue; // Maxima only with seed energy at least
      
      maxEList[iDigitN] = en ;
      
      //printf("Local max %d, id %d, en %f\n", iDigit,absIdList[iDigitN],en);
      iDigitN++ ; 
    }
  }
  
  if ( iDigitN == 0 )
  {
    AliDebug(1,Form("No local maxima found, assign highest energy cell as maxima, id %d, en cell %2.2f, en cluster %2.2f",
                    idmax,emax,cluster->E()));
    iDigitN      = 1     ;
    maxEList[0]  = emax  ;
    absIdList[0] = idmax ; 
  }
  
  
  AliDebug(1,Form("In cluster E %2.2f (wth non lin. %2.2f), M02 %2.2f, M20 %2.2f, N maxima %d",
                  cluster->E(),eCluster, cluster->GetM02(),cluster->GetM20(), iDigitN));
  
//  if(fDebug > 1) for(Int_t imax = 0; imax < iDigitN; imax++)
//  {
//    printf(" \t i %d, absId %d, Ecell %f\n",imax,absIdList[imax],maxEList[imax]);
//  }
  
  return iDigitN ;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEMCALIsoPhoton::AreNeighbours(Short_t absId1, Short_t absId2)
{
  // Calculate the energy of cross cells around the leading cell.
  Bool_t neighbourhood = kFALSE;

  AliVCaloCells *cells = fVCells;
  if (!cells)
    return neighbourhood;

  if (!fGeom)
    return neighbourhood;

  Int_t iSupMod1 = -1;
  Int_t iTower1  = -1;
  Int_t iIphi1   = -1;
  Int_t iIeta1   = -1;
  Int_t iphi1   = -1;
  Int_t ieta1    = -1;
  Int_t iSupMod2 = -1;
  Int_t iTower2  = -1;
  Int_t iIphi2   = -1;
  Int_t iIeta2   = -1;
  Int_t iphi2   = -1;
  Int_t ieta2  = -1;

  //first cell
  fGeom->GetCellIndex(absId1,iSupMod1,iTower1,iIphi1,iIeta1);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod1,iTower1,iIphi1, iIeta1,iphi1,ieta1);
  //second cell
  fGeom->GetCellIndex(absId2,iSupMod2,iTower2,iIphi2,iIeta2);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod2,iTower2,iIphi2, iIeta2,iphi2,ieta2);

  Int_t aphidiff = TMath::Abs(iphi1-iphi1);
    if (aphidiff>1)
      return neighbourhood;
    Int_t aetadiff = TMath::Abs(ieta1-ieta1);
    if (aetadiff>1)
      return neighbourhood;
    if ((aphidiff + aetadiff) == 1) {
      neighbourhood = kTRUE;
    }

  return neighbourhood;
}
//________________________________________________________________________
AliVCaloCells* AliAnalysisTaskEMCALIsoPhoton::GetVCaloCells()
{
  AliVCaloCells *cells = 0;
  cells = fESDCells;
  if (!cells)
    cells = fAODCells;
  if (!cells)
    return 0;
  return cells;
}
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::FillInvMass()
{
  TObjArray *clusters = fESDClusters;

  if (!clusters)
    clusters = fAODClusters;
  if (!clusters)
    return;
  const Int_t nclus = clusters->GetEntries();
  Bool_t isCPV = kFALSE;
  Bool_t isCPV2 = kFALSE;
  for(Int_t ic = 0; ic<nclus-1; ic++){
      AliVCluster *c = static_cast<AliVCluster*>(clusters->At(ic));
      if(!c)
	continue;
      if(TMath::Abs(c->GetTrackDx())>0.03 || TMath::Abs(c->GetTrackDz())>0.02)
	isCPV = kTRUE;
      /*      if(!isCPV)
	      continue;*/
      Short_t idm1;
      Double_t Emax1 = GetMaxCellEnergy( c, idm1);
      Float_t ec1 = GetCrossEnergy(c,idm1);
      if(Emax1 == 0 || ec1/Emax1<0.03)
	continue;
      Float_t clsPos[3] = {0,0,0};
      c->GetPosition(clsPos);
      TVector3 clsVec(clsPos);
      clsVec -= fVecPv;
      for(Int_t jc = ic+1;jc<nclus;jc++){
	AliVCluster *c2 = static_cast<AliVCluster*>(clusters->At(jc));
	if(!c2)
	  continue;
	Short_t idm2;
	Double_t Emax2 = GetMaxCellEnergy( c2, idm2);
	Float_t ec2 = GetCrossEnergy(c2,idm2);
	if(Emax2==0 || ec2/Emax2<0.03)
	  continue;
	Float_t clsPos2[3] = {0,0,0};
	c2->GetPosition(clsPos2);
	TVector3 clsVec2(clsPos2);
	clsVec2 -= fVecPv;
	if(TMath::Abs(c2->GetTrackDx())>0.03 || TMath::Abs(c2->GetTrackDz())>0.02)
	  isCPV2 = kTRUE;
	TLorentzVector lv1,lv2,lvm;
	lv1.SetPtEtaPhiM(c->E()*TMath::Sin(clsVec.Theta()),clsVec.Eta(),clsVec.Phi(),0.0);
	lv2.SetPtEtaPhiM(c2->E()*TMath::Sin(clsVec2.Theta()),clsVec2.Eta(),clsVec2.Phi(),0.0);
	lvm = lv1 + lv2;
	if(TMath::Abs(lvm.M()-0.135)<0.015){
	  fTrackDzDxIM->Fill(c->GetTrackDz(),c->GetTrackDx());
	  fTrackDzDxIM->Fill(c2->GetTrackDz(),c2->GetTrackDx());
	}
	if(TMath::Abs(lvm.M()-0.135)>0.015 && TMath::Abs(lvm.M()-0.135)<0.02){
	  fTrackDzDxIM_bg->Fill(c->GetTrackDz(),c->GetTrackDx());
	  fTrackDzDxIM_bg->Fill(c2->GetTrackDz(),c2->GetTrackDx());
	}
	if(!isCPV2 || !isCPV)
	  continue;
	fClusInvMassPairEt->Fill(lvm.Pt(),lvm.M());
	if(TMath::Abs(lvm.M()-0.135)<0.015){
	  GetEDistInClusCells(c,idm1);
	  GetEDistInClusCells(c2,idm2);

	  if(c->E()>c2->E()){
	    fM02vsESoftPi0Kid->Fill(c2->E(),c2->GetM02());
	    fM02vsEHardPi0Kid->Fill(c->E(),c->GetM02());
	    fPi0KidClusEtPairEt->Fill(lvm.Pt(),c->E()*TMath::Sin(clsVec.Theta()));
	  }
	  else{
	    fM02vsEHardPi0Kid->Fill(c2->E(),c2->GetM02());
	    fM02vsESoftPi0Kid->Fill(c->E(),c->GetM02());
	    fPi0KidClusEtPairEt->Fill(lvm.Pt(),c2->E()*TMath::Sin(clsVec2.Theta()));
	  }
	}
	if((lvm.M()>0.1 && lvm.M()<0.11) || (lvm.M()>0.15 && lvm.M()<0.16)){
	  if(c->E()>c2->E()){
	    fM02vsESoftPi0BGKid->Fill(c2->E(),c2->GetM02());
	    fM02vsEHardPi0BGKid->Fill(c->E(),c->GetM02());
	  }
	  else{
	    fM02vsEHardPi0BGKid->Fill(c2->E(),c2->GetM02());
	    fM02vsESoftPi0BGKid->Fill(c->E(),c->GetM02());
	  }
	}
      }
    }
}
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::GetEDistInClusCells(const AliVCluster *cluster, Short_t idmax)
{
  // Calculate the energy of cross cells around the leading cell.

  AliVCaloCells *cells = fVCells;
  if (!cells)
    return ;

  if (!fGeom)
    return ;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Int_t iphis   = -1;
  Int_t ietas   = -1;

  fGeom->GetCellIndex(idmax,iSupMod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphis,ietas);

  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Int_t cellAbsId = cluster->GetCellAbsId(i);
    fGeom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
    Int_t phidiff = iphi-iphis;
    /*    if (aphidiff>1)
	  continue;*/
    Int_t etadiff = ieta-ietas;
    /*if (aetadiff>1)
      continue;*/
    if(cluster->E()>5.5 && cluster->E()<6.5)
      fCellsPi0KidE6->Fill(etadiff,phidiff,cells->GetCellAmplitude(cellAbsId));
    if(cluster->E()>10.5 && cluster->E()<11.5)
      fCellsPi0KidE11->Fill(etadiff,phidiff,cells->GetCellAmplitude(cellAbsId));
    if(cluster->E()>12.5 && cluster->E()<13.5)
      fCellsPi0KidE13->Fill(etadiff,phidiff,cells->GetCellAmplitude(cellAbsId));
    if(cluster->GetM02()>0.10 && cluster->GetM02()<0.15)
      fCellsPi0KidM021st->Fill(etadiff,phidiff,cells->GetCellAmplitude(cellAbsId));
    if(cluster->GetM02()>0.15 && cluster->GetM02()<0.20)
      fCellsPi0KidM022nd->Fill(etadiff,phidiff,cells->GetCellAmplitude(cellAbsId));
    if(cluster->GetM02()>0.20 && cluster->GetM02()<0.25)
      fCellsPi0KidM023rd->Fill(etadiff,phidiff,cells->GetCellAmplitude(cellAbsId));
    if(cluster->GetM02()>0.25 && cluster->GetM02()<0.30)
      fCellsPi0KidM024th->Fill(etadiff,phidiff,cells->GetCellAmplitude(cellAbsId));
  }
}
//________________________________________________________________________
Bool_t AliAnalysisTaskEMCALIsoPhoton::IsPi0M02(Double_t M02, Double_t Et)
{
    Double_t M02u;
    if(Et<12)
      M02u = 0.02486*Et*Et - 0.7289*Et + 6.266;
    else
      M02u = 14.32/Et - 0.09863;
    if(M02u<0.65)
      M02u = 0.65;
    Double_t M02l = 12.88/Et - 0.3194;
    if(M02l<0.4)
      M02l = 0.4;
    if(M02<M02u && M02>M02l)
      return kTRUE;
    return kFALSE;
}
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::Terminate(Option_t *) 
{
  // Called once at the end of the query.
}
