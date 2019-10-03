// $Id$

#include "AliAnalysisTaskEMCALPhotonTagged.h"

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
#include "AliDataFile.h"


#include <iostream>
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEMCALPhotonTagged)

//________________________________________________________________________
AliAnalysisTaskEMCALPhotonTagged::AliAnalysisTaskEMCALPhotonTagged() : 
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
  fNDimensions(7),
  fECut(3.),
  fTrackMult(0),        
  fMcIdFamily(""),
  fNClusForDirPho(0),
  fDirPhoPt(0),
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
  fNCuts(5),
  fPileUpRejSPD(kFALSE),
  fDistToBadChan(0),
  fSigmaSmear(0.0),
  fNLMCut(100),
  fLoEtTag(10.),
  fHiEtTag(12.),
  fTagType(1),
  fEtMax(0.),
  fMaxClusLV(0,0,0,0),
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
  fMCDirPhotonPtEtaPhi(0),
  fDecayPhotonPtMC(0),
  fCellAbsIdVsAmpl(0),       
  fNClusHighClusE(0),
  fClusEtMcPt(0),
  fClusMcDetaDphi(0),
  fNClusPerPho(0),
  fInvMassPt(0),
  fMaxEtSpec(0),
  fMaxEtSpecCut(0),
  fMaxEtSpecCutTM(0),
  fTagCandEtType(0),
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
  fETrigg(0)
{
  // Default constructor.
  for(Int_t i = 0; i < 12;    i++)  fGeomMatrix[i] =  0;
}

//________________________________________________________________________
AliAnalysisTaskEMCALPhotonTagged::AliAnalysisTaskEMCALPhotonTagged(const char *name) : 
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
  fNDimensions(7),
  fECut(3.),
  fTrackMult(0),        
  fMcIdFamily(""),
  fNClusForDirPho(0),
  fDirPhoPt(0),
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
  fNCuts(5),
  fPileUpRejSPD(kFALSE),
  fDistToBadChan(0),
  fSigmaSmear(0.0),
  fNLMCut(100),
  fLoEtTag(10.),
  fHiEtTag(12.),
  fTagType(1),
  fEtMax(0.),
  fMaxClusLV(0,0,0,0),
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
  fMCDirPhotonPtEtaPhi(0),
  fDecayPhotonPtMC(0),
  fCellAbsIdVsAmpl(0),       
  fNClusHighClusE(0),   
  fClusEtMcPt(0),
  fClusMcDetaDphi(0),
  fNClusPerPho(0),
  fInvMassPt(0),
  fMaxEtSpec(0),
  fMaxEtSpecCut(0),
  fMaxEtSpecCutTM(0),
  fTagCandEtType(0),
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
  fETrigg(0)
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
void AliAnalysisTaskEMCALPhotonTagged::UserCreateOutputObjects()
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

  fMCDirPhotonPtEtaPhi = new TH3F("hMCDirPhotonPtEtaPhi","photon (gq->#gammaq) p_{T}, #eta, #phi;GeV/c;#eta;#phi",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,154,-0.77,0.77,130,1.38,3.20);
  fMCDirPhotonPtEtaPhi->Sumw2();
  fOutputList->Add(fMCDirPhotonPtEtaPhi);


  fDecayPhotonPtMC = new TH1F("hDecayPhotonPtMC","decay photon p_{T};GeV/c;dN/dp_{T} (c GeV^{-1})",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge);
  fDecayPhotonPtMC->Sumw2();
  fOutputList->Add(fDecayPhotonPtMC);

  fCellAbsIdVsAmpl = new TH2F("hCellAbsIdVsAmpl","cell abs id vs cell amplitude (energy);E (GeV);ID",200,0,100,24*48*10,-0.5,24*48*10-0.5);
  fOutputList->Add(fCellAbsIdVsAmpl);

  fNClusHighClusE = new TH2F("hNClusHighClusE","total number of clusters vs. highest clus energy in the event;E (GeV);NClus",200,0,100,301,-0.5,300.5);
  fOutputList->Add(fNClusHighClusE);

  fClusEtMcPt = new TH2F("hClusEtMcPt","E_{T}^{clus} vs. p_{T}^{mc}; p_{T}^{mc};E_{T}^{clus}",500,0,100,500,0,100);
  fOutputList->Add(fClusEtMcPt);

  fClusMcDetaDphi = new TH2F("hClusMcDetaDphi","#Delta#phi vs. #Delta#eta (reco-mc);#Delta#eta;#Delta#phi",100,-.7,.7,100,-.7,.7);
  fOutputList->Add(fClusMcDetaDphi);

  fNClusPerPho = new TH2F("hNClusPerPho","Number of clusters per prompt photon;p_{T}^{MC};N_{clus}",500,0,100,11,-0.5,10.5);
  fOutputList->Add(fNClusPerPho);

  fInvMassPt = new TH2F("hInvMassPt","leading cluster inv mass vs pair p_T;p_T;m_{#gamma#gamma}",300,0,30,1000,0,10);
  fOutputList->Add(fInvMassPt);

  fMaxEtSpec = new TH1F("hMaxEtSpec","E_T spectrum of leading clusters; E_T; dN/dE_T",500,0,100);
  fOutputList->Add(fMaxEtSpec);

  fMaxEtSpecCut = new TH1F("hMaxEtSpecCut","E_T spectrum of leading clusters (after E_T cut); E_T; dN/dE_T",500,0,100);
  fOutputList->Add(fMaxEtSpecCut);

  fMaxEtSpecCutTM = new TH1F("hMaxEtSpecCutTM","E_T spectrum of leading clusters (after E_T cut and TM); E_T; dN/dE_T",500,0,100);
  fOutputList->Add(fMaxEtSpecCutTM);

  fTagCandEtType = new TH2F("hTagCandEtType","Tag candidate type vs E_T;E_T (GeV);type",500,0,100,5,-0.5,4.5);
  fOutputList->Add(fTagCandEtType);

  Int_t nEt=fNBinsPt*5, nM02=400, nTrClDphi=200, nTrClDeta=100, nClEta=140, nClPhi=128, nTime=60, nMult=100, nPhoMcPt=fNBinsPt,  nNLM=11;
  Int_t bins[] = {nEt, nM02, nTrClDphi, nTrClDeta,nClEta,nClPhi,nTime,nMult,nPhoMcPt, nNLM};
  fNDimensions = sizeof(bins)/sizeof(Int_t);
  const Int_t ndims =   fNDimensions;
  Double_t xmin[] = { fPtBinLowEdge,   0., -0.1,-0.05, -0.7, 1.4,-0.15e-06,-0.5,fPtBinLowEdge,-0.5};
  Double_t xmax[] = { fPtBinHighEdge, 4.,  0.1, 0.05, 0.7, 3.192, 0.15e-06,99.5,fPtBinHighEdge,10.5};
  if(fPeriod.Contains("11h")){
    xmax[7]=3999.5;
  }
  fHnOutput =  new THnSparseF("fHnOutput","Output matrix: E_{T},M02, d#phi_{trk},d#eta_{trk},#eta_{clus},#phi_{clus},T_{max},mult,mc-p_{T}^{#gamma},NLM", ndims, bins, xmin, xmax);
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

  fETrigg = new TH1F("fETrigg","TrigPatchInfo->GetPatchE();E [GeV];entries/GeV",100,0,100);
  fETrigg->Sumw2();

  fQAList->Add(fETrigg);


  PostData(1, fOutputList);
  PostData(2, fQAList);
}

//________________________________________________________________________
void AliAnalysisTaskEMCALPhotonTagged::UserExec(Option_t *) 
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
  TList *l = 0;
  TString clusArrayName = "";
  if(fESD){
    l = fESD->GetList();
    /*if(fDebug)
      l->Print();*/
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
  Int_t leadId = GetLeadEtClusId();
  if(fEtMax<fLoEtTag || fEtMax>fHiEtTag){
    ClearAll();
    return;
  }
  Int_t candStatus = TagEvent(leadId);
  if(0==candStatus && candStatus != fTagType){
    ClearAll();
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
  FollowGamma();
  CheckTriggerPatch();
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
  ClearAll();
  PostData(1, fOutputList);
  PostData(2, fQAList);
}      

//________________________________________________________________________
void AliAnalysisTaskEMCALPhotonTagged::FillClusHists()
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
    Bool_t isCPV = kFALSE;
    if(TMath::Abs(c->GetTrackDx())>0.03 || TMath::Abs(c->GetTrackDz())>0.02)
      isCPV = kTRUE;


    Int_t clusNLM = GetNumberOfLocalMaxima(c,fVCells);
    const Int_t ndims =   fNDimensions;
    Double_t outputValues[ndims];
    ptmc = GetClusSource(c);
    outputValues[0] = Et;
    outputValues[1] = SmearM02(c->GetM02());
    if(fDebug)
      printf("track-cluster dphi=%1.3f, deta=%1.3f\n",c->GetTrackDx(),c->GetTrackDz());
    if(TMath::Abs(c->GetTrackDx())<0.1)
      outputValues[2] = c->GetTrackDx();
    else
      outputValues[2] = 0.099*c->GetTrackDx()/TMath::Abs(c->GetTrackDx());
    if(TMath::Abs(c->GetTrackDz())<0.05)
      outputValues[3] = c->GetTrackDz();
    else
      outputValues[3] = 0.049*c->GetTrackDz()/TMath::Abs(c->GetTrackDz());
    outputValues[4] = clsVec.Eta();
    outputValues[5] = clsVec.Phi();
    if(fESDCells)
      outputValues[6] = fESDCells->GetCellTime(id);
    else if(fAODCells)
      outputValues[6] = fAODCells->GetCellTime(id);
    outputValues[7] = fTrackMult;
    outputValues[8] = ptmc;
    if(clusNLM<10)
      outputValues[9] = clusNLM;
    else 
      outputValues[9] = 10;
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
Double_t AliAnalysisTaskEMCALPhotonTagged::GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax)
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
Double_t AliAnalysisTaskEMCALPhotonTagged ::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
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
void AliAnalysisTaskEMCALPhotonTagged ::FillMcHists()
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
Float_t AliAnalysisTaskEMCALPhotonTagged::GetClusSource(const AliVCluster *c)
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
void AliAnalysisTaskEMCALPhotonTagged::FollowGamma()
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
void AliAnalysisTaskEMCALPhotonTagged::GetDaughtersInfo(int firstd, int lastd, int selfid, const char *inputind)
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
bool AliAnalysisTaskEMCALPhotonTagged::IsMcPDG(Int_t label, Int_t PDG)
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
void AliAnalysisTaskEMCALPhotonTagged::FillQA() 
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
Double_t AliAnalysisTaskEMCALPhotonTagged::GetTrackMatchedPt(Int_t matchIndex)
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
bool AliAnalysisTaskEMCALPhotonTagged::IsExotic(AliVCluster *c)
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
void AliAnalysisTaskEMCALPhotonTagged::CheckTriggerPatch()
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
Double_t AliAnalysisTaskEMCALPhotonTagged::SmearM02(Double_t m02)
{
  if(!fStack && !fAODMCParticles)
    return m02;
  TRandom3 *r = new TRandom3(0);
  return m02+r->Gaus(0,fSigmaSmear);
}
//________________________________________________________________________
Int_t AliAnalysisTaskEMCALPhotonTagged::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells) 
{  
  const Int_t   nc = cluster->GetNCells();
  
  Int_t   absIdList[nc]; 
  Float_t maxEList[nc]; 
  
  Int_t nMax = GetNumberOfLocalMaxima(cluster, cells, absIdList, maxEList);
  
  return nMax;
}
//________________________________________________________________________
Int_t AliAnalysisTaskEMCALPhotonTagged::GetNumberOfLocalMaxima(AliVCluster* cluster, AliVCaloCells* cells,
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
Bool_t AliAnalysisTaskEMCALPhotonTagged::AreNeighbours(Short_t absId1, Short_t absId2)
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
AliVCaloCells* AliAnalysisTaskEMCALPhotonTagged::GetVCaloCells()
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
Bool_t AliAnalysisTaskEMCALPhotonTagged::IsPi0M02(Double_t M02, Double_t Et)
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
Bool_t AliAnalysisTaskEMCALPhotonTagged::HasPi0InvMass(Int_t iCand, TObjArray *clusters)
{
  Bool_t hasIt = kFALSE;
  if(!clusters)
    return kFALSE;
  Int_t nclus = clusters->GetEntriesFast();
  AliVCluster *ccand = static_cast<AliVCluster*>(clusters->At(iCand));

  if(!ccand)
    return kTRUE;
  TLorentzVector lv, pairlv; 
  Double_t invMass=0;
  for(Int_t ic=0;ic<nclus;ic++){
    if(ic==iCand)
      continue;
    AliVCluster *c = static_cast<AliVCluster*>(clusters->At(ic));
    if(!c)
      continue;
    if(c->GetTrackDx()<0.03 && c->GetTrackDz()<0.02)
      continue;
    if(c->GetM02()<0.09 && c->GetM02()>0.31)
      continue;
    Float_t clsPos[3] = {0,0,0};
    c->GetPosition(clsPos);
    TVector3 clsVec(clsPos);
    clsVec -= fVecPv;
    Double_t Et = c->E()*TMath::Sin(clsVec.Theta());
    lv.SetPtEtaPhiM(Et,clsVec.Eta(),clsVec.Phi(),0.);
    pairlv = lv + fMaxClusLV;
    invMass = pairlv.M();
    if(invMass>0.115 && invMass<0.155)
      hasIt = kTRUE;
    fInvMassPt->Fill(invMass,pairlv.Pt());
  }
  return hasIt;
}
//________________________________________________________________________
Int_t AliAnalysisTaskEMCALPhotonTagged::GetLeadEtClusId()
{
  if(fDebug)
    printf("Getting Highest Et cluster...");
  TObjArray *clusters = fESDClusters;
  if (!clusters){
    clusters = fAODClusters;
  }
  if (!clusters){
    return kFALSE;
  }
  Int_t nc = clusters->GetEntriesFast();
  Int_t idMax = -1;
  for(Int_t ic=0;ic<nc;ic++){
    AliVCluster *c = static_cast<AliVCluster*>(clusters->At(ic));
    if(!c)
      continue;
    if(!c->IsEMCAL())
      continue;
    Float_t clsPos[3] = {0,0,0};
    c->GetPosition(clsPos);
    TVector3 clsVec(clsPos);
    clsVec -= fVecPv;
    Double_t Et = c->E()*TMath::Sin(clsVec.Theta());
    if(Et<fEtMax)
      continue;
    fEtMax = Et;
    idMax = ic;
    fMaxClusLV.SetPtEtaPhiM(Et,clsVec.Eta(),clsVec.Phi(),0);
  }
  fMaxEtSpec->Fill(fEtMax);
  return idMax;
}
//________________________________________________________________________
Int_t AliAnalysisTaskEMCALPhotonTagged::TagEvent(Int_t idMax)
{
  if(fDebug)
    printf("Getting Tag Candidates...");
  Int_t candStatus = 0; //0=no cand, 1=photon, 2=pi0 (merged)
  Bool_t hasPi0InvMass = kFALSE;
  TObjArray *clusters = fESDClusters;
  if (!clusters){
    clusters = fAODClusters;
  }
  if (!clusters){
    return kFALSE;
  } 
  fMaxEtSpecCut->Fill(fEtMax);
  Int_t nc = clusters->GetEntriesFast();
  AliVCluster *c = static_cast<AliVCluster*>(clusters->At(idMax));
  if(!c)
    return 0;
  if(!c->IsEMCAL())
    return 0;
  if(c->GetTrackDx()<0.03 && c->GetTrackDz()<0.02){
    fMaxEtSpecCutTM->Fill(fEtMax);
    return 0;
  }
  Float_t clsPos[3] = {0,0,0};
  c->GetPosition(clsPos);
  TVector3 clsVec(clsPos);
  clsVec -= fVecPv;
  Double_t Et = c->E()*TMath::Sin(clsVec.Theta());
  if(Et!=fEtMax)
    AliFatal("Et not consistent with max Et in event!!!\n");
  if(c->GetM02()>0.09 && c->GetM02()<0.31){
    candStatus = 1;
    if(HasPi0InvMass(idMax,clusters))
      candStatus = 2;
  }
  if(IsPi0M02(Et,c->GetM02()))
    candStatus = 2;
  fTagCandEtType->Fill(Et,candStatus);
  return candStatus;
}
//________________________________________________________________________
void AliAnalysisTaskEMCALPhotonTagged::ClearAll()
{
  if(fESD)
    fESDClusters->Clear();
  fSelPrimTracks->Clear();
  fNClusForDirPho = 0;
  fNCells50 = 0;
  fClusIdFromTracks = "";
  fVecPv.Clear();
  fEtMax = 0.;
  fMaxClusLV.SetPtEtaPhiM(0,0,0,0);
}
//________________________________________________________________________
void AliAnalysisTaskEMCALPhotonTagged::Terminate(Option_t *) 
{
  // Called once at the end of the query.
}
