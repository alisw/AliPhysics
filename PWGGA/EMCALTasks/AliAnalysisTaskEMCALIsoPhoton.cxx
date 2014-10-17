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
#include "AliEMCALGeometry.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"
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
  fESD(0),
  fAOD(0),
  fVEvent(0),
  fMCEvent(0),
  fStack(0),
  fOutputList(0),
  fEvtSel(0),
  fNClusEt10(0),
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
  fInvMassWithConeVsEtAndIso(0),
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
  fEoverPvsE(0)
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
  fESD(0),
  fAOD(0),
  fVEvent(0),
  fMCEvent(0),
  fStack(0),
  fOutputList(0),
  fEvtSel(0),
  fNClusEt10(0),
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
  fInvMassWithConeVsEtAndIso(0),
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
  fEoverPvsE(0)
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
  fOADBContainer->InitFromFile(Form("$ALICE_ROOT/OADB/EMCAL/EMCALlocal2master.root"),"AliEMCALgeo");
 
  fEvtSel = new TH1F("hEvtSel","Event selection counter (0=all trg, 1=pvz cut) ;evt cut ;dN/dcut}",2,0,2);
  fOutputList->Add(fEvtSel);
  
  fNClusEt10 = new TH1F("hNClusEt10","# of cluster with E_{T}>10 per event;E;",101,-0.5,100.5);
  fOutputList->Add(fNClusEt10);
  
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

  fInvMassWithConeVsEtAndIso = new TH3F("hInvMassWithConeVsEtAndIso","M_{cand+in_cone_clus} vs E_{T}^{cand} vs. E_{T}^{ISO} (EMC+Trk) (0.1<M02<0.3 only)",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,400,0,1,1000,0,200);
  fOutputList->Add(fInvMassWithConeVsEtAndIso);

  fInConePairedClusEtVsCandEt = new TH2F("hInConePairedClusEtVsCandEt","E_{T}^{partner} vs. E_{T}^{cand} (R<0.4, 0.110<m_{#gamma#gamma}<0.165);E_{T}^{cand};E_{T}^{partner}",fNBinsPt, fPtBinLowEdge,fPtBinHighEdge,200,0,40);
  fOutputList->Add(fInConePairedClusEtVsCandEt);

  Int_t nEt=fNBinsPt*5, nM02=400, nCeIso=1000, nTrIso=1000,  nAllIso=1000,  nCeIsoNoUE=1000,  nAllIsoNoUE=1000, nTrClDphi=200, nTrClDeta=100, nClEta=140, nClPhi=128, nTime=60, nMult=100, nPhoMcPt=fNBinsPt;
  Int_t bins[] = {nEt, nM02, nCeIso, nTrIso, nAllIso, nCeIsoNoUE, nAllIsoNoUE, nTrClDphi, nTrClDeta,nClEta,nClPhi,nTime,nMult,nPhoMcPt};
  fNDimensions = sizeof(bins)/sizeof(Int_t);
  const Int_t ndims =   fNDimensions;
  Double_t xmin[] = { fPtBinLowEdge,   0.,  -10.,   -10., -10., -10., -10., -0.1,-0.05, -0.7, 1.4,-0.15e-06,-0.5,fPtBinLowEdge};
  Double_t xmax[] = { fPtBinHighEdge, 4., 190., 190., 190.,  190., 190., 0.1, 0.05, 0.7, 3.192, 0.15e-06,99.5,fPtBinHighEdge};
  if(fPeriod.Contains("11h")){
    xmax[12]=3999.5;
  }
  fHnOutput =  new THnSparseF("fHnOutput","Output matrix: E_{T},M02,CeIso,TrIso,AllIso, CeIsoNoUESub, AllIsoNoUESub, d#phi_{trk},d#eta_{trk},#eta_{clus},#phi_{clus},T_{max},mult,mc-p_{T}^{#gamma}", ndims, bins, xmin, xmax);
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

  if(fESD){
    TList *l = fESD->GetList();
    if(fDebug){
      for(int nk=0;nk<l->GetEntries();nk++){
	  TObject *obj = (TObject*)l->At(nk);
	  TString oname = obj->GetName();
	  //if(oname.Contains("lus"))
	    printf("Object %d has a clus array named %s +++++++++\n",nk,oname.Data());
	}
    }
    fESDClusters =  dynamic_cast<TClonesArray*>(l->FindObject("CaloClusters"));
    fESDCells = fESD->GetEMCALCells();
    if(fDebug)
      printf("ESD cluster mult= %d\n",fESDClusters->GetEntriesFast());
  }
  else if(fAOD){
    fAODClusters = dynamic_cast<TClonesArray*>(fAOD->GetCaloClusters());
    fAODCells = fAOD->GetEMCALCells();
    if(fDebug)
      printf("AOD cluster mult= %d\n",fAODClusters->GetEntriesFast());
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
  LoopOnCells();
  FollowGamma();
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
    if(c->GetM02()>0.1 && c->GetM02()<0.3 && isCPV){
      TObjArray *inConeInvMassArr = (TObjArray*)fInConeInvMass.Tokenize(";");
      Int_t nInConePairs = inConeInvMassArr->GetEntriesFast();
      Double_t *inConeInvMass = new Double_t[nInConePairs];
      for(int ipair=0;ipair<nInConePairs;ipair++){
	TObjString *obs = (TObjString*)inConeInvMassArr->At(ipair);
	TString smass = obs->GetString();
	Double_t pairmass = smass.Atof();
	if(fDebug)
	  printf("=================+++++++++++++++Inv mass inside the cone for photon range: %1.1f,%1.1f,%1.1f+-++++-+-+-+-++-+-+-\n",Et,pairmass,ceiso+triso);
	fInvMassWithConeVsEtAndIso->Fill(Et,pairmass,ceiso+triso);
      }
      delete [] inConeInvMass;
    }
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
    const Int_t ndims =   fNDimensions;
    Double_t outputValues[ndims];
    if(mcptsum<2)
      ptmc = GetClusSource(c);
    else
      ptmc = -0.1;
    outputValues[0] = Et;
    outputValues[1] = c->GetM02();
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
  AliVCaloCells *cells = fESDCells;
  if (!cells){
    cells = fAODCells;
    if(fDebug)
      printf("ESD cells empty...");
  }
  if (!cells){
     if(fDebug)
      printf("  and AOD cells are empty  as well!!!\n"); 
    return;
  }

  TObjArray *clusters = fESDClusters;
  if (!clusters)
    clusters = fAODClusters;
  if (!clusters)
    return;
  
  fInConeInvMass = "";
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
    Double_t matchedpt =  0;
    if(fCpvFromTrack){
      matchedpt = GetTrackMatchedPt(c->GetTrackMatchedIndex());
      if(matchedpt>0 && fRemMatchClus )
	continue;
    } else {
      if(TMath::Abs(c->GetTrackDx())<0.03 && TMath::Abs(c->GetTrackDz())<0.02){
	matchedpt = GetTrackMatchedPt(c->GetTrackMatchedIndex());
	if(fRemMatchClus){
	  if(fDebug)
	    printf("This isolation cluster is matched to a track!++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	  continue;
	}
      }
    }
    Double_t nEt = TMath::Max(Et-matchedpt, 0.0);
    if(nEt<0)
      printf("nEt=%1.1f\n",nEt);
    if(R<fIsoConeR){
      totiso += nEt;
      if(c->GetM02()>0.1 && c->GetM02()<0.3 && !(matchedpt>0)){
	TLorentzVector lv, lvec;
	lv.SetPtEtaPhiM(Et,cv.Eta(),cv.Phi(),0);
	lvec.SetPtEtaPhiM(EtCl,vec.Eta(),vec.Phi(),0);
	TLorentzVector lpair = lv + lvec;
	fInConeInvMass += Form(";%f",lpair.M());
	if(lpair.M()>0.11 && lpair.M()<0.165)
	  fInConePairedClusEtVsCandEt->Fill(EtCl,Et);
      }
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

  AliVCaloCells *cells = 0;
  cells = fESDCells;
  if (!cells)
    cells = fAODCells;
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

  AliVCaloCells *cells = 0;
  cells = fESDCells;
  if (!cells)
    cells = fAODCells;
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
    Int_t lab1 =  mcp->GetDaughter(0);
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
    Int_t daug0f =  mcp->GetFirstDaughter();
    Int_t daug0l =  mcp->GetLastDaughter();
    Int_t nd0 = daug0l - daug0f;
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
    Int_t daug0f =  mcp->GetDaughter(0);
    Int_t daug0l =  mcp->GetDaughter(mcp->GetNDaughters()-1);
    Int_t nd0 = daug0l - daug0f;
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
      Int_t dfd = dp->GetDaughter(0); 
      Int_t dld = dp->GetDaughter(dp->GetNDaughters()-1);
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
    if(TMath::Abs(c->GetTrackDx())<0.03 && TMath::Abs(c->GetTrackDz())<0.02){
      fEmcClusETM1->Fill(c->E());
      continue;
    }
    fEmcClusEClusCuts->Fill(c->E(),2);
    if(TMath::Abs(maxt)>30e-9)
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
    if(!fPrTrCuts && !fCompTrCuts)
      return pt;
    if(!fPrTrCuts->IsSelected(track) && !fCompTrCuts->IsSelected(track))
      return pt;
    pt = track->Pt();
  }
  return pt;
}
//________________________________________________________________________
void AliAnalysisTaskEMCALIsoPhoton::LoopOnCells()
{
  AliVCaloCells *cells = 0;
  cells = fESDCells;
  if (!cells)
    cells = fAODCells;
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
void AliAnalysisTaskEMCALIsoPhoton::Terminate(Option_t *) 
{
  // Called once at the end of the query.
}
