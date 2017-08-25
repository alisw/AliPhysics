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


////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Task for Heavy Flavour Electron-Hadron DeltaPhi Correlation       //
//  Non-Photonic Electron identified with Invariant mass              //
//  analysis methos in function  SelectPhotonicElectron               //
//  DeltaPhi calculated in function  ElectronHadCorrel                //
//                                                                    //
//  Author: Deepa Thomas (University of Texas at Austin)              //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"

#include "AliAnalysisTaskEHCorrel.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"


#include "AliEventPoolManager.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliESDtrackCuts.h"
#include "AliCentralitySelectionTask.h"
#include "AliMultSelection.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDCaloTrigger.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"

#include "AliCentrality.h"
#include "AliMagF.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TProfile.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "TVector3.h"
#include "TRandom2.h"

  ClassImp(AliAnalysisTaskEHCorrel)
ClassImp(AliehDPhiBasicParticle)
  //________________________________________________________________________
AliAnalysisTaskEHCorrel::AliAnalysisTaskEHCorrel(const char *name)
  : AliAnalysisTaskSE(name),
  fVevent(0),
  fAOD(0),
  fpVtx(0),
  fpidResponse(0),
  fMultSelection(0),
  fHistElecEffi(0),
  fCentrality(-1),
  fMultiplicity(-1),
  fTracks_tender(0),
  fCaloClusters_tender(0),
  fApplyAddPileUpCuts(kFALSE),
  fUseTender(kFALSE),
  fCentralityMin(0),
  fCentralityMax(20),
  fEMCEG1(kFALSE),
  fEMCEG2(kFALSE),
  fFlagClsTypeEMC(kTRUE),
  fFlagClsTypeDCAL(kTRUE),
  fTPCNClsElec(90),
  fFlagEleSPDkFirst(kFALSE),
  fTPCnSigma(-999.0),
  fTPCnSigmaMin(-1),
  fTPCnSigmaMax(3),
  fM02Min(0.01),
  fM02Max(0.35),
  fM20Min(0),
  fM20Max(2),
  fEovPMin(0.9),
  fEovPMax(1.2),
  fTPCNClsHad(80),
  fTPCNClsPartnerE(80),
  fPartElePt(0.3),
  fInvmassCut(0.1),
  fFlagHadSPDkAny(kFALSE),
  fFlagHadITSNCls(kFALSE),
  fFlagHadFiducialCut(kFALSE),
  fFlagHadPosEtaOnly(kFALSE),
  fFlagHadNegEtaOnly(kFALSE),
  fTPCnSigmaHadMin(-10),
  fTPCnSigmaHadMax(-3.5),
  fHadCutCase(2),
  fPoolMgr(0x0),
  fTrigElePtCut(kFALSE),
  fNEle(0),
  fVtxZBin(-999),
  fCentBin(-999),
  fFlagMEBinChange(kFALSE),
  fIsPbPb(kTRUE),
  fEMCClsTimeCut(kFALSE),
  fApplyElectronEffi(kFALSE),
  fEffi(1.0),
  fWeight(1.0),
  fOutputList(0),
  fNevents(0),
  fVtxZ(0),
  fVtxX(0),
  fVtxY(0),
  fCentralityNoPass(0),
  fCentralityPass(0),
  fMultiplicityNoPass(0),
  fMultiplicityPass(0),
  fCentMultiplicityNoPass(0),
  fCentMultiplicityPass(0),
  fHistClustE(0),
  fEMCClsEtaPhi(0),
  fHistoNCells(0),
  fHistoTimeEMC(0),
  fNegTrkIDPt(0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fdEdx(0),
  fTPCnsig(0),
  fHistPtMatch(0),
  fEMCTrkMatch(0),
  fEMCTrkPt(0),
  fEMCTrketa(0),
  fEMCTrkphi(0),
  fEMCTPCnsig(0),
  fClsEAftMatch(0),
  fClsEtaPhiAftMatch(0),
  fHistNsigEop(0),
  fM20EovP(0),
  fM02EovP(0),
  fHistEop(0),
  fM20(0),
  fM02(0),
  fHistEop_AftEID(0),
  fInclsElecPt(0),
  fNElecInEvt(0),
  fHadEop(0),
  fHadPt_AftEID(0),
  fULSElecPt(0),
  fLSElecPt(0),
  fTagULSElecPt(0),
  fTagLSElecPt(0),
  fHadronPhiPt(0),
  fHadronPhi(0),
  fHadronPhiTPChalf(0),
  fHadronPt(0),
  fInvmassLS(0),
  fInvmassULS(0),
  fInvmassLSPt(0),
  fInvmassULSPt(0),
  fNoMixedEvents(0),
  fMixStatCent(0),
  fMixStatVtxZ(0),
  fMixStatCentVtxz(0),
  // fHisHadDphi(0),
  // fHisIncEDphi(0),
  // fHisLSDphi(0),
  // fHisULSDphi(0),
  fSprsAllHadHCorrl(0),
  fSprsMixAllHadHCorrl(0),
  fSprsHadHCorrl(0),
  fSprsInclusiveEHCorrl(0),
  fSprsLSEHCorrl(0),
  fSprsULSEHCorrl(0),
  fSprsLSNoPartnerEHCorrl(0),
  fSprsULSNoPartnerEHCorrl(0),
  fSprsTagULSEHCorrl(0),
  fSprsTagLSEHCorrl(0),
  fSprsMixHadHCorrl(0),
  fSprsMixInclusiveEHCorrl(0),
  fSprsMixLSEHCorrl(0),
  fSprsMixULSEHCorrl(0),
  fSprsMixTagULSEHCorrl(0),
  fSprsMixTagLSEHCorrl(0)
{
  //Named constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEHCorrel::AliAnalysisTaskEHCorrel()
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskEHCorrel"),
  fVevent(0),
  fAOD(0),
  fpVtx(0),
  fpidResponse(0),
  fMultSelection(0),
  fHistElecEffi(0),
  fCentrality(-1),
  fCentralityMin(0),
  fCentralityMax(20),
  fMultiplicity(-1),
  fTracks_tender(0),
  fCaloClusters_tender(0),
  fApplyAddPileUpCuts(kFALSE),
  fUseTender(kFALSE),
  fEMCEG1(kFALSE),
  fEMCEG2(kFALSE),
  fFlagClsTypeEMC(kTRUE),
  fFlagClsTypeDCAL(kTRUE),
  fTPCNClsElec(90),
  fFlagEleSPDkFirst(kFALSE),
  fTPCnSigma(-999.0),
  fTPCnSigmaMin(-1),
  fTPCnSigmaMax(3),
  fM02Min(0.01),
  fM02Max(0.35),
  fM20Min(0),
  fM20Max(2),
  fEovPMin(0.9),
  fEovPMax(1.2),
  fTPCNClsHad(80),
  fTPCNClsPartnerE(80),
  fPartElePt(0.3),
  fInvmassCut(0.1),
  fFlagHadSPDkAny(kFALSE),
  fFlagHadITSNCls(kFALSE),
  fFlagHadFiducialCut(kFALSE),
  fFlagHadPosEtaOnly(kFALSE),
  fFlagHadNegEtaOnly(kFALSE),
  fTPCnSigmaHadMin(-10),
  fTPCnSigmaHadMax(-3.5),
  fHadCutCase(1),
  fPoolMgr(0x0),
  fTrigElePtCut(kFALSE),
  fNEle(0),
  fVtxZBin(-999),
  fCentBin(-999),
  fFlagMEBinChange(kFALSE),
  fIsPbPb(kTRUE),
  fEMCClsTimeCut(kFALSE),
  fApplyElectronEffi(kFALSE),
  fEffi(1.0),
  fWeight(1.0),
  fOutputList(0),
  fNevents(0),
  fVtxZ(0),
  fVtxX(0),
  fVtxY(0),
  fCentralityNoPass(0),
  fCentralityPass(0),
  fMultiplicityNoPass(0),
  fMultiplicityPass(0),
  fCentMultiplicityNoPass(0),
  fCentMultiplicityPass(0),
  fHistClustE(0),
  fEMCClsEtaPhi(0),
  fHistoNCells(0),
  fHistoTimeEMC(0),
  fNegTrkIDPt(0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fdEdx(0),
  fTPCnsig(0),
  fHistPtMatch(0),
  fEMCTrkMatch(0),
  fEMCTrkPt(0),
  fEMCTrketa(0),
  fEMCTrkphi(0),
  fEMCTPCnsig(0),
  fClsEAftMatch(0),
  fClsEtaPhiAftMatch(0),
  fHistNsigEop(0),
  fM20EovP(0),
  fM02EovP(0),
  fHistEop(0),
  fM20(0),
  fM02(0),
  fHistEop_AftEID(0),
  fInclsElecPt(0),
  fNElecInEvt(0),
  fHadEop(0),
  fHadPt_AftEID(0),
  fULSElecPt(0),
  fLSElecPt(0),
  fTagULSElecPt(0),
  fTagLSElecPt(0),
  fHadronPhiPt(0),
  fHadronPhi(0),
  fHadronPhiTPChalf(0),
  fHadronPt(0),
  fInvmassLS(0),
  fInvmassULS(0),
  fInvmassLSPt(0),
  fInvmassULSPt(0),
  fNoMixedEvents(0),
  fMixStatCent(0),
  fMixStatVtxZ(0),
  fMixStatCentVtxz(0),
  // fHisHadDphi(0),
  // fHisIncEDphi(0),
  // fHisLSDphi(0),
  // fHisULSDphi(0),
  fSprsAllHadHCorrl(0),
  fSprsMixAllHadHCorrl(0),
  fSprsHadHCorrl(0),
  fSprsInclusiveEHCorrl(0),
  fSprsLSEHCorrl(0),
  fSprsULSEHCorrl(0),
  fSprsLSNoPartnerEHCorrl(0),
  fSprsULSNoPartnerEHCorrl(0),
  fSprsTagULSEHCorrl(0),
  fSprsTagLSEHCorrl(0),
  fSprsMixHadHCorrl(0),
  fSprsMixInclusiveEHCorrl(0),
  fSprsMixLSEHCorrl(0),
  fSprsMixULSEHCorrl(0),
  fSprsMixTagULSEHCorrl(0),
  fSprsMixTagLSEHCorrl(0)
{
  //Default constructor
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //DefineOutput(3, TTree::Class());
}
//_________________________________________
AliAnalysisTaskEHCorrel::~AliAnalysisTaskEHCorrel()
{
  //Destructor

  delete fOutputList;
  delete fTracks_tender;
  delete fCaloClusters_tender;
  delete fSprsHadHCorrl;
  delete fSprsInclusiveEHCorrl;
  delete fSprsLSEHCorrl;
  delete fSprsULSEHCorrl;
  delete fSprsMixInclusiveEHCorrl;
  delete fSprsMixLSEHCorrl;
  delete fSprsMixULSEHCorrl;
  delete fSprsMixHadHCorrl;
  delete fSprsAllHadHCorrl;
  delete fSprsMixAllHadHCorrl;
  delete fSprsLSNoPartnerEHCorrl;
  delete fSprsULSNoPartnerEHCorrl;
  delete fSprsTagULSEHCorrl;
  delete fSprsTagLSEHCorrl;
  delete fSprsMixTagULSEHCorrl;
  delete fSprsMixTagLSEHCorrl;
    if(fHistElecEffi) {delete fHistElecEffi; fHistElecEffi=0;}

}
//_________________________________________
void AliAnalysisTaskEHCorrel::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  AliDebug(3, "Creating Output Objects");
    
    if(fApplyElectronEffi){
        TString elecEffiFileName;
        
        if(!fIsPbPb){
            elecEffiFileName = "alien:///alice/cern.ch/user/d/dthomas/HFElecEffi_pPb/HFElectronTrackEffi.root";
        }
        if(fIsPbPb){
            elecEffiFileName = "alien:///alice/cern.ch/user/d/dthomas/HFElecEffi_PbPb/HFElectronTrackEffi.root";
        }
        
        TFile* f2 = TFile::Open(elecEffiFileName.Data());
        TH1D *h = (TH1D*)f2->Get("ElecEffi");
        SetElectronEffiMap(h);
    }

  Double_t pi = TMath::Pi();

  ////////////////////////
  //Initiale mixed event//
  ////////////////////////
  Int_t trackDepth = 0;
  if(fIsPbPb) trackDepth = 100000;
  if(!fIsPbPb) trackDepth = 100000;

  Int_t poolsize   = 1000;

  Int_t nZvtxBins  = 6;
  Double_t vertexBins[7];

  Int_t nCentralityBinsPbPb = 6;
  Double_t CentralityBinsPbPb[7];
  Int_t nCentralityBinspp = 4;
  Double_t CentralityBinspp[5];

  if(fIsPbPb){
    if(!fFlagMEBinChange){ //mean of VtxZ is at 0.5
      vertexBins[0] = -10;
      vertexBins[1] = -5;
      vertexBins[2] = -2;
      vertexBins[3] = 0.5;
      vertexBins[4] = 3;
      vertexBins[5] = 6;
      vertexBins[6] = 10;
    }
    if(fFlagMEBinChange){
      vertexBins[0] = -10;
      vertexBins[1] = -5;
      vertexBins[2] = -2.5;
      vertexBins[3] = 0;
      vertexBins[4] = 2.5;
      vertexBins[5] = 5;
      vertexBins[6] = 10;
    }

    if(fCentralityMax == 20)
    {
      if(!fFlagMEBinChange){
        CentralityBinsPbPb[0] = 0;
        CentralityBinsPbPb[1] = 1.5;
        CentralityBinsPbPb[2] = 3.5;
        CentralityBinsPbPb[3] = 6;
        CentralityBinsPbPb[4] = 9;
        CentralityBinsPbPb[5] = 14;
        CentralityBinsPbPb[6] = 20;
      }
      if(fFlagMEBinChange){
        CentralityBinsPbPb[0] = 0;
        CentralityBinsPbPb[1] = 2;
        CentralityBinsPbPb[2] = 4;
        CentralityBinsPbPb[3] = 6.5;
        CentralityBinsPbPb[4] = 10;
        CentralityBinsPbPb[5] = 15;
        CentralityBinsPbPb[6] = 20;
      }
    }
    if(fCentralityMax == 50)
    {
      if(!fFlagMEBinChange){
        CentralityBinsPbPb[0] = 20;
        CentralityBinsPbPb[1] = 25;
        CentralityBinsPbPb[2] = 30;
        CentralityBinsPbPb[3] = 35;
        CentralityBinsPbPb[4] = 40;
        CentralityBinsPbPb[5] = 45;
        CentralityBinsPbPb[6] = 50;
      }
      if(fFlagMEBinChange){
        CentralityBinsPbPb[0] = 20;
        CentralityBinsPbPb[1] = 24;
        CentralityBinsPbPb[2] = 29;
        CentralityBinsPbPb[3] = 35;
        CentralityBinsPbPb[4] = 40;
        CentralityBinsPbPb[5] = 45;
        CentralityBinsPbPb[6] = 50;
      }
    }
    if(fCentralityMax > 50)
    {
      if(!fFlagMEBinChange){
        CentralityBinsPbPb[0] = 50;
        CentralityBinsPbPb[1] = 55;
        CentralityBinsPbPb[2] = 60;
        CentralityBinsPbPb[3] = 65;
        CentralityBinsPbPb[4] = 70;
        CentralityBinsPbPb[5] = 75;
        CentralityBinsPbPb[6] = 80;
      }
    }
  }
  if(!fIsPbPb){
    if(!fFlagMEBinChange){
      vertexBins[0] = -10;
      vertexBins[1] = -4.6;
      vertexBins[2] = -1.6;
      vertexBins[3] = 0.9;
      vertexBins[4] = 3.4;
      vertexBins[5] = 6.1;
      vertexBins[6] = 10;

      CentralityBinspp[0] = 0;
      CentralityBinspp[1] = 25;
      CentralityBinspp[2] = 50;
      CentralityBinspp[3] = 75;
      CentralityBinspp[4] = 100;
    }
    if(fFlagMEBinChange){
      vertexBins[0] = -10;
      vertexBins[1] = -5;
      vertexBins[2] = -2.5;
      vertexBins[3] = 0;
      vertexBins[4] = 2.5;
      vertexBins[5] = 5;
      vertexBins[6] = 10;

      CentralityBinspp[0] = 0;
      CentralityBinspp[1] = 20;
      CentralityBinspp[2] = 40;
      CentralityBinspp[3] = 60;
      CentralityBinspp[4] = 100;
    }
  }

  if(fIsPbPb)
    fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBinsPbPb, (Double_t*) CentralityBinsPbPb, nZvtxBins, (Double_t*) vertexBins);

  if(!fIsPbPb)
    fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBinspp, (Double_t*) CentralityBinspp, nZvtxBins, (Double_t*) vertexBins);

  ///////////////
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

  fCentralityPass = new TH1F("fCentralityPass", "Centrality Pass;centrality;counts", 101, -1, 100);
  fOutputList->Add(fCentralityPass);

  fCentralityNoPass = new TH1F("fCentralityNoPass", "Centrality No Pass;centrality;counts", 101, -1, 100);
  fOutputList->Add(fCentralityNoPass);

  fMultiplicityPass = new TH1F("fMultiplicityPass", "Multiplicity Pass;multiplicity;counts", 2500, 0, 25000);
  fOutputList->Add(fMultiplicityPass);

  fMultiplicityNoPass = new TH1F("fMultiplicityNoPass", "Multiplicity No Pass;multiplicity;counts", 2500, 0, 25000);
  fOutputList->Add(fMultiplicityNoPass);

  fCentMultiplicityNoPass = new TH2F("fCentMultiplicityNoPass", "Multiplicity vs Centrality, No Pass;multiplicity;centrality", 2500, 0, 25000, 101, -1, 100);
  fOutputList->Add(fCentMultiplicityNoPass);

  fCentMultiplicityPass = new TH2F("fCentMultiplicityPass", "Multiplicity vs Centrality, Pass;multiplicity;centrality", 2500, 0, 25000, 101, -1, 100);
  fOutputList->Add(fCentMultiplicityPass);

  fHistClustE = new TH1F("fHistClustE", "EMCAL cluster energy distribution; Cluster E;counts", 500, 0.0, 50.0);
  fOutputList->Add(fHistClustE);

  fEMCClsEtaPhi = new TH2F("fEMCClsEtaPhi","EMCAL cluster #eta and #phi distribution;#eta;#phi",100,-0.9,0.9,200,0,6.3);
  fOutputList->Add(fEMCClsEtaPhi);

  fHistoNCells = new TH2F("fHistoNCells","No of EMCAL cells in a cluster;Cluster E;N^{EMC}_{cells}",300,0,30,50,0,50);
  fOutputList->Add(fHistoNCells);

  fHistoTimeEMC = new TH2F("fHistoTimeEMC","EMCAL Time;E (GeV); t(ns)",500,0,50,1800,-900,900);
  fOutputList->Add(fHistoTimeEMC);

  fNegTrkIDPt = new TH1F("fNegTrkIDPt", "p_{T} distribution of tracks with negative track id;p_{T} (GeV/c);counts", 500, 0.0, 50.0);
  fOutputList->Add(fNegTrkIDPt);

  fTrkPt = new TH1F("fTrkPt","p_{T} distribution of all tracks;p_{T} (GeV/c);counts",1000,0,100);
  fOutputList->Add(fTrkPt);

  fTrketa = new TH1F("fTrketa","All Track #eta distribution;#eta;counts",100,-1.5,1.5);
  fOutputList->Add(fTrketa);

  fTrkphi = new TH1F("fTrkphi","All Track #phi distribution;#phi;counts",100,0,2*pi);
  fOutputList->Add(fTrkphi);

  fdEdx = new TH2F("fdEdx","All Track dE/dx distribution;p (GeV/c);dE/dx",200,0,20,500,0,160);
  fOutputList->Add(fdEdx);

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

  fEMCTrkphi = new TH1F("fEMCTrkphi","#phi distribution of tracks matched to EMCAL;#phi;counts",100,0,2*pi);
  fOutputList->Add(fEMCTrkphi);


  fEMCTPCnsig = new TH2F("fEMCTPCnsig","TPC Nsigma distribution of tracks matched to EMCAL;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
  fOutputList->Add(fEMCTPCnsig);

  fClsEAftMatch = new TH1F("fClsEAftMatch", "EMCAL cluster energy distribution after track matching; Cluster E;counts", 500, 0.0, 50.0);
  fOutputList->Add(fClsEAftMatch);

  fClsEtaPhiAftMatch = new TH2F("fClsEtaPhiAftMatch","EMCAL cluster #eta and #phi distribution after track matching;#eta;#phi",100,-0.9,0.9,200,0,6.3);
  fOutputList->Add(fClsEtaPhiAftMatch);

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

  fHistEop = new TH2F("fHistEop", "E/p distribution;p_{T} (GeV/c);E/p", 200,0,20,60, 0.0, 3.0);
    fHistEop->Sumw2();
  fOutputList->Add(fHistEop);

  fHistEop_AftEID = new TH2F("fHistEop_AftEID", "E/p distribution after nsig, SS cuts;p_{T} (GeV/c);E/p", 200,0,20,60, 0.0, 3.0);
    fHistEop_AftEID->Sumw2();
  fOutputList->Add(fHistEop_AftEID);

  fInclsElecPt = new TH1F("fInclsElecPt","p_{T} distribution of inclusive electrons;p_{T} (GeV/c);counts",500,0,50);
  fInclsElecPt->Sumw2();
  fOutputList->Add(fInclsElecPt);

  fNElecInEvt = new TH1F("fNElecInEvt","No of electrons in the event; N^{ele};counts",20,-0.5,19.5);
  fOutputList->Add(fNElecInEvt);

  fHadEop = new TH2F("fHadEop", "E/p distribution for hadrons;p_{T} (GeV/c);E/p", 200,0,20,60, 0.0, 3.0);
    fHadEop->Sumw2();
  fOutputList->Add(fHadEop);
    
  fHadPt_AftEID = new TH1F("fHadPt_AftEID","p_{T} distribution of hadrons after Eid cuts;p_{T} (GeV/c);counts",500,0,50);
  fHadPt_AftEID->Sumw2();
  fOutputList->Add(fHadPt_AftEID);

  fULSElecPt  = new TH1F("fULSElecPt","p_{T} distribution of ULS electrons;p_{T} (GeV/c);counts",500,0,50);
  fULSElecPt->Sumw2();
  fOutputList->Add(fULSElecPt);

  fLSElecPt= new TH1F("fLSElecPt","p_{T} distribution of LS electrons;p_{T} (GeV/c);counts",500,0,50);
  fLSElecPt->Sumw2();
  fOutputList->Add(fLSElecPt);

  fTagULSElecPt  = new TH1F("fTagULSElecPt","p_{T} distribution of tagged ULS electrons;p_{T} (GeV/c);counts",500,0,50);
  fTagULSElecPt->Sumw2();
  fOutputList->Add(fTagULSElecPt);

  fTagLSElecPt= new TH1F("fTagLSElecPt","p_{T} distribution of tagged LS electrons;p_{T} (GeV/c);counts",500,0,50);
  fTagLSElecPt->Sumw2();
  fOutputList->Add(fTagLSElecPt);

  fHadronPhiPt = new TH2F("fHadronPhiPt", "Hadron phi vs pt; hadron phi; pt (GeV/c)",1000,0,2*pi,500,0,50);
  fOutputList->Add(fHadronPhiPt);

  fHadronPhi = new TH1F("fHadronPhi", "Hadron phi",1000,0,2*pi);
  fOutputList->Add(fHadronPhi);

  fHadronPhiTPChalf = new TH1F("fHadronPhiTPChalf", "Hadron phi for 0<eta<0.9",1000,0,2*pi);
  fOutputList->Add(fHadronPhiTPChalf);

  fHadronPt = new TH1F("fHadronPt","hadron pt distribution",500,0,50);
  fOutputList->Add(fHadronPt);

  fInvmassLS = new TH1F("fInvmassLS", "Inv mass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassLS);

  fInvmassULS = new TH1F("fInvmassULS", "Inv mass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassULS);

  fInvmassLSPt = new TH2F("fInvmassLSPt", "Inv mass of LS (e,e) vs pT; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,1000,0,1.0);
    fInvmassLSPt->Sumw2();
  fOutputList->Add(fInvmassLSPt);

  fInvmassULSPt = new TH2F("fInvmassULSPt", "Inv mass of ULS (e,e) vs pT; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,1000,0,1.0);
    fInvmassULSPt->Sumw2();
  fOutputList->Add(fInvmassULSPt);

  fNoMixedEvents = new TH1F("fNoMixedEvents","No of mixing events",1,-0.5,0.5);
  fOutputList->Add(fNoMixedEvents);

  Int_t nEventBins =500;
  Double_t EventBins[nEventBins+1];
  for(int i=0; i < nEventBins+1; i++)
    EventBins[i] = i;

  if(fIsPbPb){
    fMixStatCent = new TH2F("fMixStatCent","Mix event stats for centrality binning;Nevent in pool;Centrality",nEventBins,EventBins,nCentralityBinsPbPb,CentralityBinsPbPb);
    fOutputList->Add(fMixStatCent);

    fMixStatCentVtxz = new TH2F("fMixStatCentVtxz","Mix event stats Cent vs Zvtx binning;Vtx_{z};Centrality",nZvtxBins,vertexBins,nCentralityBinsPbPb,CentralityBinsPbPb);
    fOutputList->Add(fMixStatCentVtxz);
  }

  if(!fIsPbPb){
    fMixStatCent = new TH2F("fMixStatCent","Mix event stats for centrality binning;Nevent in pool;Centrality",nEventBins,EventBins,nCentralityBinspp,CentralityBinspp);
    fOutputList->Add(fMixStatCent);

    fMixStatCentVtxz = new TH2F("fMixStatCentVtxz","Mix event stats Cent vs Zvtx binning;Vtx_{z};Centrality",nZvtxBins,vertexBins,nCentralityBinspp,CentralityBinspp);
    fOutputList->Add(fMixStatCentVtxz);
  }

  fMixStatVtxZ = new TH2F("fMixStatVtxZ","Mix event stats for Zvtx binning;Nevent in pool;Vtx_{z}",nEventBins,EventBins,nZvtxBins,vertexBins);
  fOutputList->Add(fMixStatVtxZ);

  //  fHisHadDphi = new TH2F("fHisHadDphi","Hadron Dphi;p_{T}^{e};#Delta#varphi",50,0,50,64,-TMath::Pi()/2,(3*TMath::Pi())/2);
  //  fOutputList->Add(fHisHadDphi);

  //  fHisIncEDphi = new TH2F("fHisIncEDphi","InclE Dphi;p_{T}^{e};#Delta#varphi",50,0,50,64,-TMath::Pi()/2,(3*TMath::Pi())/2);
  //  fOutputList->Add(fHisIncEDphi);

  //  fHisLSDphi = new TH2F("fHisLSDphi","LS electron Dphi;p_{T}^{e};#Delta#varphi",50,0,50,64,-TMath::Pi()/2,(3*TMath::Pi())/2);
  //  fOutputList->Add(fHisLSDphi);

  //  fHisULSDphi = new TH2F("fHisULSDphi","ULS electron Dphi;p_{T}^{e};#Delta#varphi",50,0,50,64,-TMath::Pi()/2,(3*TMath::Pi())/2);
  //  fOutputList->Add(fHisULSDphi);


  //------THnsparse------
  // Int_t bin[6] = {30,20,32,50,nZvtxBins,nCentralityBinsPbPb}; //ptElec, ptHad,Dphi, Deta
  // Double_t xmin[6] = {0,0,-TMath::Pi()/2,-1.8,0,0};
  // Double_t xmax[6] = {30,20,(3*TMath::Pi())/2,1.8,6,6};
  Int_t bin[4] = {30,20,32,50}; //ptElec, ptHad,Dphi, Deta
  Double_t xmin[4] = {0,0,-TMath::Pi()/2,-1.8};
  Double_t xmax[4] = {30,20,(3*TMath::Pi())/2,1.8};

  fSprsAllHadHCorrl = new THnSparseD("fSprsAllHadHCorrl","Sparse for Dphi and Deta for all hadrons;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsAllHadHCorrl->Sumw2();
  fOutputList->Add(fSprsAllHadHCorrl);

  fSprsMixAllHadHCorrl = new THnSparseD("fSprsMixAllHadHCorrl","Sparse for Dphi and Deta for all hadrons for mixed event;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsMixAllHadHCorrl->Sumw2();
  fOutputList->Add(fSprsMixAllHadHCorrl);

  fSprsHadHCorrl = new THnSparseD("fSprsHadHCorrl","Sparse for Dphi and Deta hadrons;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsHadHCorrl->Sumw2();
  fOutputList->Add(fSprsHadHCorrl);

  fSprsInclusiveEHCorrl = new THnSparseD("fSprsInclusiveEHCorrl","Sparse for Dphi and Deta with Inclusive electron;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsInclusiveEHCorrl->Sumw2();
  fOutputList->Add(fSprsInclusiveEHCorrl);

  fSprsLSEHCorrl = new THnSparseD("fSprsLSEHCorrl","Sparse for Dphi and Deta with LS Non-HFE electron;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsLSEHCorrl->Sumw2();
  fOutputList->Add(fSprsLSEHCorrl);

  fSprsULSEHCorrl = new THnSparseD("fSprsULSEHCorrl","Sparse for Dphi and Deta with ULS Non-HFE electron;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsULSEHCorrl->Sumw2();
  fOutputList->Add(fSprsULSEHCorrl);

  fSprsLSNoPartnerEHCorrl = new THnSparseD("fSprsLSNoPartnerEHCorrl","Sparse for Dphi and Deta with LS Non-HFE electron, No partner;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsLSNoPartnerEHCorrl->Sumw2();
  fOutputList->Add(fSprsLSNoPartnerEHCorrl);

  fSprsULSNoPartnerEHCorrl = new THnSparseD("fSprsULSNoPartnerEHCorrl","Sparse for Dphi and Deta with ULS Non-HFE electron, No partner;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsULSNoPartnerEHCorrl->Sumw2();
  fOutputList->Add(fSprsULSNoPartnerEHCorrl);

  fSprsTagULSEHCorrl = new THnSparseD("fSprsTagULSEHCorrl","Sparse for Dphi and Deta with Tagged ULS Non-HFE electron;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsTagULSEHCorrl->Sumw2();
  fOutputList->Add(fSprsTagULSEHCorrl);

  fSprsTagLSEHCorrl = new THnSparseD("fSprsTagLSEHCorrl","Sparse for Dphi and Deta with Tagged LS Non-HFE electron;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsTagLSEHCorrl->Sumw2();
  fOutputList->Add(fSprsTagLSEHCorrl);

  fSprsMixHadHCorrl = new THnSparseD("fSprsMixHadHCorrl","Sparse for Dphi and Deta hadrons for mixed events;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsMixHadHCorrl->Sumw2();
  fOutputList->Add(fSprsMixHadHCorrl);

  fSprsMixInclusiveEHCorrl = new THnSparseD("fSprsMixInclusiveEHCorrl","Sparse for Dphi and Deta with Inclusive electron for mixed event;p_{T}^{e};p_{T}^{mixTrk};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsMixInclusiveEHCorrl->Sumw2();
  fOutputList->Add(fSprsMixInclusiveEHCorrl);

  fSprsMixLSEHCorrl = new THnSparseD("fSprsMixLSEHCorrl","Sparse for Dphi and Deta with LS Non-HFE electron for mixed events;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsMixLSEHCorrl->Sumw2();
  fOutputList->Add(fSprsMixLSEHCorrl);

  //fSprsMixULSEHCorrl = new THnSparseD("fSprsMixULSEHCorrl","Sparse for Dphi and Deta with ULS Non-HFE electron for mixed events;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;VtxZBin,CentBin;",6,bin,xmin,xmax);
  fSprsMixULSEHCorrl = new THnSparseD("fSprsMixULSEHCorrl","Sparse for Dphi and Deta with ULS Non-HFE electron for mixed events;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsMixULSEHCorrl->Sumw2();
  fOutputList->Add(fSprsMixULSEHCorrl);

  fSprsMixTagULSEHCorrl = new THnSparseD("fSprsMixTagULSEHCorrl","Sparse for Dphi and Deta with tagged ULS Non-HFE electron for mixed events;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsMixTagULSEHCorrl->Sumw2();
  fOutputList->Add(fSprsMixTagULSEHCorrl);

  fSprsMixTagLSEHCorrl = new THnSparseD("fSprsMixTagLSEHCorrl","Sparse for Dphi and Deta with tagged LS Non-HFE electron for mixed events;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fSprsMixTagLSEHCorrl->Sumw2();
  fOutputList->Add(fSprsMixTagLSEHCorrl);

  PostData(1,fOutputList);
}
//_________________________________________
void AliAnalysisTaskEHCorrel::UserExec(Option_t*)
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
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  fpVtx = fVevent->GetPrimaryVertex();

  ///////////////////
  //PID initialised//
  ///////////////////
  fpidResponse = fInputHandler->GetPIDResponse();

  /////////////////
  //trigger check//
  /////////////////
  TString firedTrigger;
  TString TriggerEG1("EG1");
  TString TriggerEG2("EG2");
  if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();

  if(fEMCEG1){if(!firedTrigger.Contains(TriggerEG1))return;}
  if(fEMCEG2){if(!firedTrigger.Contains(TriggerEG2))return;}

  ////////////////////
  //event selection///
  ////////////////////
  if(!PassEventSelect(fVevent)) return;

  if(fApplyAddPileUpCuts){
    if(!PassAddtionalPileUpCuts()) return;
  }

  /////////////////
  // Centrality ///
  /////////////////
  Bool_t pass = kFALSE;
  if(fCentralityMin > -0.5){
    CheckCentrality(fAOD,pass);
    if(!pass)return;
  }
  /////////////////////////
  //Get VtxZ and Cent Bin//
  /////////////////////////
  GetVtxZCentralityBin();

  //////////////
  //if Tender //
  //////////////
  if(fUseTender){
    //new branches with calibrated tracks and clusters
    fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("tracks"));
    fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("caloClusters"));
  }

  //////////////////////
  //EMcal cluster info//
  //////////////////////
  EMCalClusterInfo();

  /////////////////////
  //Hadron track info//
  /////////////////////
  HadronInfo();

  ///////////////
  //Track loop///
  ///////////////
  fNEle = 0;
  Int_t ntracks = -999;
  if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries();
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliVParticle* Vtrack = 0x0;
    if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
    if(fUseTender) Vtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));
    if (!Vtrack) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
    AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
    AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);

    ////////////////////
    //Apply track cuts//
    ////////////////////
    if(!PassTrackCuts(atrack)) continue;

    Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999;
    TrkPhi = track->Phi();
    TrkPt = track->Pt();
    TrkEta = track->Eta();
    TrkP = track->P();

    if(TrkPt > 2){
      ElectronHadCorrel(iTracks, track, fSprsAllHadHCorrl);
      MixedEvent(track, fSprsMixAllHadHCorrl);
    }

    ///////////////////////////
    //Track matching to EMCAL//
    //////////////////////////
    Int_t EMCalIndex = -1;
    EMCalIndex = track->GetEMCALcluster();
    if(EMCalIndex < 0) continue;
    fHistPtMatch->Fill(TrkPt);

    AliVCluster *clustMatch=0x0;
    if(!fUseTender) clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex);
    if(fUseTender) clustMatch = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(EMCalIndex));

    Double_t emcphi = -999, emceta=-999;
    Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
    if(clustMatch && clustMatch->IsEMCAL())
    {
      Double_t fPhiDiff = -999, fEtaDiff = -999;
      GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);
      fEMCTrkMatch->Fill(fPhiDiff,fEtaDiff);

      if(TMath::Abs(fPhiDiff) > 0.05 || TMath::Abs(fEtaDiff)> 0.05) continue;

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

      Double_t clustTime = clustMatch->GetTOF()*1e+9; // ns;

      if(fEMCClsTimeCut)
        if(TMath::Abs(clustTime) > 50) continue;

      /////////////////////////////////////////////
      //Properties of tracks matched to the EMCAL//
      /////////////////////////////////////////////
      fEMCTrkPt->Fill(TrkPt);
      fEMCTrketa->Fill(TrkEta);
      fEMCTrkphi->Fill(TrkPhi);
      fEMCTPCnsig->Fill(TrkP,fTPCnSigma);
      Double_t clustMatchE = clustMatch->E();
      fClsEAftMatch->Fill(clustMatchE);
      fClsEtaPhiAftMatch->Fill(emceta,emcphi);

      //Select pT>2 GeV/c
      if(TrkPt < 2) continue;

      //////////////////
      //Apply EID cuts//
      //////////////////
      Bool_t fHadTrack = kFALSE, fElectTrack = kFALSE;
      fElectTrack = PassEIDCuts(track, clustMatch, fHadTrack);
        
      //---Get electron weight------
      if(fApplyElectronEffi) {
        fEffi = GetElecEffi(track);
        fWeight = 1/fEffi;
      }

      ///////////////////
      //H-H Correlation//
      ///////////////////
      if(fHadTrack)
      {
        fHadPt_AftEID->Fill(TrkPt,fWeight);
        //ElectronHadCorrel(iTracks, track, fSprsHadHCorrl, fHisHadDphi);
        ElectronHadCorrel(iTracks, track, fSprsHadHCorrl);
        MixedEvent(track, fSprsMixHadHCorrl);
      }

      if(!fElectTrack) continue;
        
      fInclsElecPt->Fill(TrkPt,fWeight);
      fNEle++;

      ///////////////////
      //E-H Correlation//
      ///////////////////
      //HadronInfo(iTracks);

      //Inclusive E-H correl
      //ElectronHadCorrel(iTracks, track, fSprsInclusiveEHCorrl, fHisIncEDphi);
      ElectronHadCorrel(iTracks, track, fSprsInclusiveEHCorrl);
      MixedEvent(track, fSprsMixInclusiveEHCorrl);

      ////////////////////
      //NonHFE selection//
      ////////////////////
      Bool_t fFlagPhotonicElec = kFALSE, fFlagElecLS=kFALSE;
      SelectNonHFElectron(iTracks,track,fFlagPhotonicElec,fFlagElecLS);

      if(fFlagPhotonicElec)
      {
        ElectronHadCorrel(iTracks, track, fSprsTagULSEHCorrl);
        MixedEvent(track, fSprsMixTagULSEHCorrl);
        fTagULSElecPt->Fill(TrkPt,fWeight);
      }
      if(fFlagElecLS)
      {
        ElectronHadCorrel(iTracks, track, fSprsTagLSEHCorrl);
        MixedEvent(track, fSprsMixTagLSEHCorrl);
        fTagLSElecPt->Fill(TrkPt,fWeight);
      }

    }//EMCAL track match
  }//track loop

  fNElecInEvt->Fill(fNEle);

  /////////////////////////
  //Fill Mixed event pool//
  /////////////////////////
  Double_t pVtxZ = fpVtx->GetZ();
  AliEventPool *fPool;
  fPool = fPoolMgr->GetEventPool(fCentrality, pVtxZ); // Get the buffer associated with the current centrality and z-vtx
  if (!fPool)
  {
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCentrality, pVtxZ));
    return;
  }
  TObjArray * fArrayTracksMix = CloneAndReduceTrackList();
  fArrayTracksMix->SetOwner(kTRUE);
  fPool->UpdatePool(fArrayTracksMix);

  PostData(1, fOutputList);
}
//_________________________________________
//void AliAnalysisTaskEHCorrel::ElectronHadCorrel(Int_t itrack, AliVTrack *track, THnSparse *SparseEHCorrl, TH2F *HisDphi)
void AliAnalysisTaskEHCorrel::ElectronHadCorrel(Int_t itrack, AliVTrack *track, THnSparse *SparseEHCorrl)
{
  //Construct Deta Phi between electrons and hadrons

  Double_t fvalueDphi[4] = {-999,999,-999,-999}; //ptElec, ptHad,Dphi, Deta
  Double_t pi = TMath::Pi();

  Double_t ptHad= -999;
  Double_t ptEle = -999;
  Double_t phiEle = -999, phiHad = -999, Dphi = -999;
  Double_t etaEle = -999, etaHad = -999, Deta = -999;

  Int_t ntracks = -999;
  if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries();

  for(Int_t ktracks = 0; ktracks < ntracks; ktracks++){
    if(ktracks == itrack) continue; //do not select the same electron

    AliVParticle* VtrackHad = 0x0;
    if(!fUseTender) VtrackHad = fVevent->GetTrack(ktracks);
    if(fUseTender) VtrackHad = dynamic_cast<AliVTrack*>(fTracks_tender->At(ktracks)); //take tracks from Tender list

    if (!VtrackHad) {
      printf("ERROR: Could not receive track %d\n", ktracks);
      continue;
    }
    AliVTrack *trackHad = dynamic_cast<AliVTrack*>(VtrackHad);
    if(!trackHad) continue;
    AliAODTrack *atrackHad = dynamic_cast<AliAODTrack*>(VtrackHad);
    if(!atrackHad) continue;

    ptHad = trackHad->Pt();
    ptEle = track->Pt();
    phiEle = track->Phi();
    phiHad = trackHad->Phi();
    etaEle = track->Eta();
    etaHad = trackHad->Eta();

    if(!PassHadronCuts(atrackHad)) continue; //apply hadron cuts;
    if(fTrigElePtCut && (ptHad > ptEle)) continue; //select only pTe > pTh is requested

    Dphi = phiEle - phiHad;
    if (Dphi > 3*pi/2)
      Dphi = Dphi - 2*pi;
    if (Dphi < -pi/2)
      Dphi = Dphi + 2*pi;

    Deta = etaEle - etaHad;
      
    fvalueDphi[0] = ptEle;
    fvalueDphi[1] = ptHad;
    fvalueDphi[2] = Dphi;
    fvalueDphi[3] = Deta;
    //    fvalueDphi[4] = fVtxZBin;
    //    fvalueDphi[5] = fCentBin;
    SparseEHCorrl->Fill(fvalueDphi,fWeight);

    //if(ptHad > ptEle) HisDphi->Fill(ptEle,Dphi);
  }
}
//___________________________________________
Double_t AliAnalysisTaskEHCorrel::GetElecEffi(AliVTrack *track){
//Get electron efficiency
   
    Int_t bin = fHistElecEffi->FindBin(track->Pt());
    Int_t bin10 = fHistElecEffi->FindBin(10);
    
    if(fHistElecEffi->IsBinUnderflow(bin)||fHistElecEffi->IsBinOverflow(bin)) return 1.0;
    if(track->Pt()>10) return fHistElecEffi->GetBinContent(bin10);
    
    return fHistElecEffi->GetBinContent(bin);
}
//___________________________________________
void AliAnalysisTaskEHCorrel::ElectronHadCorrelNoPartner(Int_t itrack, Int_t jtrack, AliVTrack *track, THnSparse *SparseEHCorrlNoPartner)
{
  //Construct Deta Phi between electrons and hadrons for electrons from invariant mass calculation excluding partner etrack

  Double_t fvalueDphi[4] = {-999,999,-999,-999}; //ptElec, ptHad,Dphi, Deta
  Double_t pi = TMath::Pi();

  Double_t ptHad= -999;
  Double_t ptEle = -999;
  Double_t phiEle = -999, phiHad = -999, Dphi = -999;
  Double_t etaEle = -999, etaHad = -999, Deta = -999;

  Int_t ntracks = -999;
  if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries();

  for(Int_t ktracks = 0; ktracks < ntracks; ktracks++){
    if(ktracks == itrack || ktracks == jtrack) continue; //do not select the same electron and partner etrack from invmass

    AliVParticle* VtrackHad = 0x0;
    if(!fUseTender) VtrackHad = fVevent->GetTrack(ktracks);
    if(fUseTender) VtrackHad = dynamic_cast<AliVTrack*>(fTracks_tender->At(ktracks)); //take tracks from Tender list

    if (!VtrackHad) {
      printf("ERROR: Could not receive track %d\n", ktracks);
      continue;
    }
    AliVTrack *trackHad = dynamic_cast<AliVTrack*>(VtrackHad);
    if(!trackHad) continue;
    AliAODTrack *atrackHad = dynamic_cast<AliAODTrack*>(VtrackHad);
    if(!atrackHad) continue;

    ptHad = trackHad->Pt();
    ptEle = track->Pt();
    phiEle = track->Phi();
    phiHad = trackHad->Phi();
    etaEle = track->Eta();
    etaHad = trackHad->Eta();

    if(!PassHadronCuts(atrackHad)) continue; //apply hadron cuts;
    if(fTrigElePtCut && (ptHad > ptEle)) continue; //select only pTe > pTh is requested

    Dphi = phiEle - phiHad;
    if (Dphi > 3*pi/2)
      Dphi = Dphi - 2*pi;
    if (Dphi < -pi/2)
      Dphi = Dphi + 2*pi;

    Deta = etaEle - etaHad;
      
    fvalueDphi[0] = ptEle;
    fvalueDphi[1] = ptHad;
    fvalueDphi[2] = Dphi;
    fvalueDphi[3] = Deta;
    //    fvalueDphi[4] = fVtxZBin;
    //    fvalueDphi[5] = fCentBin;
    SparseEHCorrlNoPartner->Fill(fvalueDphi,fWeight);

    //if(ptHad > ptEle) HisDphi->Fill(ptEle,Dphi);
  }
}
//___________________________________________
//void AliAnalysisTaskEHCorrel::HadronInfo(Int_t itrack)
void AliAnalysisTaskEHCorrel::HadronInfo()
{
  //Hadron information

  Int_t ntracks = -999;
  if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries();

  for(Int_t ktracks = 0; ktracks < ntracks; ktracks++){

    //  if(ktracks == itrack) continue; //do not select the same electron
    AliVParticle* VtrackHad = 0x0;
    if(!fUseTender) VtrackHad = fVevent->GetTrack(ktracks);
    if(fUseTender) VtrackHad = dynamic_cast<AliVTrack*>(fTracks_tender->At(ktracks)); //take tracks from Tender list

    if (!VtrackHad) {
      printf("ERROR: Could not receive track %d\n", ktracks);
      continue;
    }
    AliVTrack *trackHad = dynamic_cast<AliVTrack*>(VtrackHad);
    if(!trackHad) continue;
    AliAODTrack *atrackHad = dynamic_cast<AliAODTrack*>(VtrackHad);
    if(!atrackHad) continue;

    if(!PassHadronCuts(atrackHad)) continue; //apply hadron cuts;

    Double_t ptHad= -999, phiHad=-999;
    ptHad = trackHad->Pt();
    phiHad = trackHad->Phi();

    fHadronPhiPt->Fill(phiHad,ptHad);
    fHadronPhi->Fill(phiHad);
    if (trackHad->Eta() >0 && trackHad->Eta() <0.9) fHadronPhiTPChalf->Fill(phiHad);
    fHadronPt->Fill(ptHad);
  }
}
//___________________________________________
Bool_t AliAnalysisTaskEHCorrel::PassHadronCuts(AliAODTrack *HadTrack)
{
  //apply hadron cuts

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
  Double_t d0z0[2]={-999,-999}, cov[3];
  Double_t DCAxyCut = 0.25, DCAzCut = 1;

  if(fHadCutCase == 1)
  {
    if(!HadTrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE;
    if((!(HadTrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(HadTrack->GetStatus()&AliESDtrack::kTPCrefit)))) return kFALSE;
  }
  if(fHadCutCase == 2)
  {
    if(!HadTrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) return kFALSE;
    if((!(HadTrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(HadTrack->GetStatus()&AliESDtrack::kTPCrefit)))) return kFALSE;
  }
  if(fHadCutCase == 3)
  {
    if(!HadTrack->IsHybridGlobalConstrainedGlobal()) return kFALSE;
  }
  if(fHadCutCase == 4)
  {
    if(!HadTrack->IsHybridTPCConstrainedGlobal()) return kFALSE;
  }

  if(HadTrack->GetTPCNcls() < fTPCNClsHad) return kFALSE;
  if(HadTrack->Eta()< -0.9 || HadTrack->Eta()>0.9) return kFALSE;
  if(HadTrack->Pt() < 0.3) return kFALSE;
  if(HadTrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
    if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) return kFALSE;

  if(fFlagHadSPDkAny){
    if(!(HadTrack->HasPointOnITSLayer(0) || HadTrack->HasPointOnITSLayer(1))) return kFALSE;
  }

  if(fFlagHadITSNCls){
    if(HadTrack->GetITSNcls() < 3) return kFALSE;
  }

  if(fFlagHadFiducialCut){
    if(HadTrack->Eta()< -0.8 || HadTrack->Eta()>0.8) return kFALSE;
  }

  if(fFlagHadPosEtaOnly){
    if(HadTrack->Eta()< 0 || HadTrack->Eta()>0.9) return kFALSE;
  }

  if(fFlagHadNegEtaOnly){
    if(HadTrack->Eta()< -0.9 || HadTrack->Eta()>0) return kFALSE;
  }

  return kTRUE;
}
//___________________________________________
Bool_t AliAnalysisTaskEHCorrel::PassEIDCuts(AliVTrack *track, AliVCluster *clust, Bool_t &Hadtrack)
{
  //apply electron identification cuts

  Bool_t hadTrk = kFALSE;
  Double_t eop = -1.0;
  Double_t m02 = -999,m20 = -999;
  Double_t clustE = clust->E();
  Double_t TrkPt = track->Pt();
  if(track->P()>0)eop = clustE/track->P();
  m02 =clust->GetM02();
  m20 =clust->GetM20();

  if(track->Pt()>2.0){
    fHistNsigEop->Fill(eop,fTPCnSigma);
    fM20EovP->Fill(eop,m20);
    fM02EovP->Fill(eop,m02);
  }
  fHistEop->Fill(TrkPt,eop);
  fM20->Fill(TrkPt,m20);
  fM02->Fill(TrkPt,m02);

  //Hadron E/p distribution
  if(fTPCnSigma > fTPCnSigmaHadMin && fTPCnSigma < fTPCnSigmaHadMax)
  {
    if((m02 > fM02Min && m02 < fM02Max) && (m20 > fM20Min && m20 < fM20Max))
    {
      fHadEop->Fill(TrkPt,eop,fWeight);
      if(eop > fEovPMin && eop < fEovPMax) hadTrk=kTRUE;
    }
  }
  Hadtrack = hadTrk;

  if(fTPCnSigma < fTPCnSigmaMin || fTPCnSigma > fTPCnSigmaMax) return kFALSE;
  if(m02 < fM02Min || m02 > fM02Max) return kFALSE;
  if(m20 < fM20Min || m20 > fM20Max) return kFALSE;

  fHistEop_AftEID->Fill(TrkPt,eop,fWeight);

  if(eop < fEovPMin || eop > fEovPMax) return kFALSE;

  return kTRUE;
}
//___________________________________________
Bool_t AliAnalysisTaskEHCorrel::PassTrackCuts(AliAODTrack *atrack)
{
  //apply track cuts

  Double_t d0z0[2]={-999,-999}, cov[3];
  Double_t DCAxyCut = 0.25, DCAzCut = 1;
  Double_t dEdx =-999;
  Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999;

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();

  //kink daughters
  Int_t numberofvertices = 100;
  numberofvertices = fAOD->GetNumberOfVertices();
  Double_t listofmotherkink[numberofvertices];
  Int_t numberofmotherkink = 0;
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
  if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE; //mimimum cuts
  //reject kink
  Bool_t kinkmotherpass = kTRUE;
  for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
    if(atrack->GetID() == listofmotherkink[kinkmother]) {
      kinkmotherpass = kFALSE;
      continue;
    }
  }
  if(!kinkmotherpass) return kFALSE;

  //other cuts
  if(atrack->GetTPCNcls() < fTPCNClsElec) return kFALSE;
  if(atrack->GetITSNcls() < 3) return kFALSE;
  if((!(atrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack->GetStatus()&AliESDtrack::kTPCrefit)))) return kFALSE;
  if(!(atrack->HasPointOnITSLayer(0) || atrack->HasPointOnITSLayer(1))) return kFALSE;

  if(fFlagEleSPDkFirst) {
    if(!(atrack->HasPointOnITSLayer(0))) return kFALSE;
  }

  if(atrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
    if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) return kFALSE;

  ////////////////////
  //Track properties//
  ////////////////////
  dEdx = atrack->GetTPCsignal();
  fTPCnSigma = fpidResponse->NumberOfSigmasTPC(atrack, AliPID::kElectron);
  TrkPhi = atrack->Phi();
  TrkPt = atrack->Pt();
  TrkEta = atrack->Eta();
  TrkP = atrack->P();

  if(atrack->GetID()<0) fNegTrkIDPt->Fill(TrkPt);
  fTrkPt->Fill(TrkPt);
  fTrketa->Fill(TrkEta);
  fTrkphi->Fill(TrkPhi);
  fdEdx->Fill(TrkP,dEdx);
  fTPCnsig->Fill(TrkP,fTPCnSigma);

  return kTRUE;
}
//___________________________________________
Bool_t AliAnalysisTaskEHCorrel::PassEventSelect(AliVEvent *fVevent)
{
  //event selection cuts

  Int_t ntracks = -999;
  ntracks = fVevent->GetNumberOfTracks();
  if(ntracks < 1) printf("There are %d tracks in this event\n",ntracks);

  fNevents->Fill(0); //all events
  Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;

  Double_t NcontV = fpVtx->GetNContributors();

  AliAODVertex* vtSPD = fAOD->GetPrimaryVertexSPD();
  Double_t NcontSPD = vtSPD->GetNContributors();
  if(NcontV<2 || NcontSPD<1)return kFALSE;

  if(!fIsPbPb){
    Double_t covTrc[6],covSPD[6];
    fpVtx->GetCovarianceMatrix(covTrc);
    vtSPD->GetCovarianceMatrix(covSPD);
    Double_t dz = fpVtx->GetZ() - vtSPD->GetZ();
    Double_t errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
    Double_t errTrc = TMath::Sqrt(covTrc[5]);
    Double_t nsigTot = TMath::Abs(dz)/errTot;
    Double_t nsigTrc = TMath::Abs(dz)/errTrc;
    if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return kFALSE;// bad vertexing
  }

  fNevents->Fill(1); //events after vettex cuts

  Zvertex = fpVtx->GetZ();
  Yvertex = fpVtx->GetY();
  Xvertex = fpVtx->GetX();
  fVtxZ->Fill(Zvertex);
  fVtxX->Fill(Xvertex);
  fVtxY->Fill(Yvertex);
  if(TMath::Abs(Zvertex)>10.0) return kFALSE;
  fNevents->Fill(2); //events after z vtx cut

  return kTRUE;
}
//___________________________________________
Bool_t AliAnalysisTaskEHCorrel::PassAddtionalPileUpCuts()
{
  //additional cuts to reject pile-up
  Int_t nTPCout=0;
  Float_t mTotV0=0;

  //get multiplicity
  AliAODVZERO* v0data=(AliAODVZERO*) fAOD->GetVZEROData();
  Float_t mTotV0A=v0data->GetMTotV0A();
  Float_t mTotV0C=v0data->GetMTotV0C();
  mTotV0=mTotV0A+mTotV0C;

  //get no of tracks with kTPCout
  Int_t ntracksEv = fAOD->GetNumberOfTracks();
  for(Int_t itrack=0; itrack<ntracksEv; itrack++) { // loop on tacks
    AliAODTrack * track = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(itrack));
    if(!track) {AliFatal("Not a standard AOD");}
    if(track->GetID()<0)continue;
    if((track->GetFlags())&(AliESDtrack::kTPCout)) nTPCout++;
    else continue;
  }
  Double_t mV0Cut=-2200.+(2.5*nTPCout)+(0.000012*nTPCout*nTPCout); //function to apply to pile-up rejection

  if(mTotV0 < mV0Cut) return kFALSE;

  return kTRUE;
}
//___________________________________________
void AliAnalysisTaskEHCorrel::CheckCentrality(AliAODEvent* fAOD, Bool_t &centralitypass)
{
  //check centrality, Run 2

  if(fAOD)fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
  if(!fMultSelection) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  }else{
    fCentrality = fMultSelection->GetMultiplicityPercentile("V0M", false);
  }

  AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
  if(!header) AliFatal("Not a standard AOD");
  fMultiplicity = header->GetRefMultiplicity();

  if ((fCentrality <= fCentralityMin) || (fCentrality > fCentralityMax))
  {
    fCentralityNoPass->Fill(fCentrality);
    fMultiplicityNoPass->Fill(fMultiplicity);
    fCentMultiplicityNoPass->Fill(fMultiplicity,fCentrality);
    //  cout << "--------------Fill no pass-------------------------"<<endl;
    centralitypass = kFALSE;
  }else
  {
    fCentralityPass->Fill(fCentrality);
    fMultiplicityPass->Fill(fMultiplicity);
    fCentMultiplicityPass->Fill(fMultiplicity,fCentrality);
    //  cout << "--------------Fill pass-------------------------"<<endl;
    centralitypass = kTRUE;
  }
}
//________________________________________________________________________
void AliAnalysisTaskEHCorrel::EMCalClusterInfo()
{
  //EMCAL cluster information

  Int_t Nclust = -999;
  TVector3 clustpos;
  Float_t  emcx[3]; // cluster pos
  Double_t clustE=-999, emcphi = -999, emceta=-999;
  Float_t tof=-999;

  if(!fUseTender) Nclust = fVevent->GetNumberOfCaloClusters();
  if(fUseTender) Nclust = fCaloClusters_tender->GetEntries();
  for(Int_t icl=0; icl<Nclust; icl++)
  {
    AliVCluster *clust = 0x0;
    if(!fUseTender) clust = fVevent->GetCaloCluster(icl);
    if(fUseTender) clust = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(icl));
    if(!clust)  printf("ERROR: Could not receive cluster matched calibrated from track %d\n", icl);

    Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;  
    if(clust && clust->IsEMCAL())
    {
      clustE = clust->E();
      tof = clust->GetTOF()*1e+9; // ns
      clust->GetPosition(emcx);
      clustpos.SetXYZ(emcx[0],emcx[1],emcx[2]);
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

      if(fEMCClsTimeCut)
        if(TMath::Abs(tof) > 50) continue;

      fHistClustE->Fill(clustE);
      fEMCClsEtaPhi->Fill(emceta,emcphi);

      fHistoNCells->Fill(clustE,clust->GetNCells());
      fHistoTimeEMC->Fill(clustE,tof);
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskEHCorrel::GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
{
  // Calculate phi and eta difference between a track and a cluster. The position of the track is obtained on the EMCAL surface

  phidiff = 999;
  etadiff = 999;

  if (!t||!v) return;

  Double_t veta = t->GetTrackEtaOnEMCal();
  Double_t vphi = t->GetTrackPhiOnEMCal();

  Float_t pos[3] = {0};
  v->GetPosition(pos);
  TVector3 cpos(pos);
  Double_t ceta     = cpos.Eta();
  Double_t cphi     = cpos.Phi();
  etadiff=veta-ceta;
  phidiff=TVector2::Phi_mpi_pi(vphi-cphi);
}
//________________________________________________________________________
void AliAnalysisTaskEHCorrel::SelectNonHFElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Bool_t &fFlagElecLS)
{
  //Photonic electron selection

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
  Double_t d0z0[2]={-999,-999}, cov[3];
  Double_t DCAxyCut = 0.25, DCAzCut = 1;

  Bool_t flagPhotonicElec = kFALSE, flagLSElec = kFALSE;
  Double_t ptAsso=-999., nsigmaAsso=-999.0;
  Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;


  Int_t ntracks = -999;
  if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries();

  for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
    if(jTracks==itrack) continue;

    AliVParticle* VtrackAsso = 0x0;
    if(!fUseTender) VtrackAsso  = fVevent->GetTrack(jTracks);
    if(fUseTender) VtrackAsso = dynamic_cast<AliVTrack*>(fTracks_tender->At(jTracks)); //take tracks from Tender list
    if (!VtrackAsso) {
      printf("ERROR: Could not receive track %d\n", jTracks);
      continue;
    }

    AliVTrack *trackAsso = dynamic_cast<AliVTrack*>(VtrackAsso);
    if(!trackAsso) continue;

    AliAODTrack *atrackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);
    if(!atrackAsso) continue;
    if(!atrackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
    if(atrackAsso->GetTPCNcls() < fTPCNClsPartnerE) continue;
    if(!(atrackAsso->GetStatus()&AliESDtrack::kTPCrefit)) continue;
    if(!(atrackAsso->GetStatus()&AliESDtrack::kITSrefit)) continue;

    nsigmaAsso = fpidResponse->NumberOfSigmasTPC(trackAsso, AliPID::kElectron);
    ptAsso = trackAsso->Pt();
    Int_t chargeAsso = trackAsso->Charge();
    Int_t charge = track->Charge();

    if(ptAsso < fPartElePt) continue;
    if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
    if(nsigmaAsso < -3 || nsigmaAsso > 3) continue;

    if(trackAsso->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
      if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;

    Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
    if(charge>0) fPDGe1 = -11;
    if(chargeAsso>0) fPDGe2 = -11;

    fFlagLS=kFALSE; fFlagULS=kFALSE;
    if(charge == chargeAsso) fFlagLS = kTRUE;
    if(charge != chargeAsso) fFlagULS = kTRUE;

    AliKFParticle::SetField(fVevent->GetMagneticField());

    AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
    AliKFParticle ge2 = AliKFParticle(*trackAsso, fPDGe2);
    AliKFParticle recg(ge1, ge2);

    if(recg.GetNDF()<1) continue;
    Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
    if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;

    Double_t mass=-999., width = -999;
    Int_t MassCorrect;
    MassCorrect = recg.GetMass(mass,width);

    if(fFlagLS && track->Pt()>1) fInvmassLS->Fill(mass);
    if(fFlagULS && track->Pt()>1) fInvmassULS->Fill(mass);
    if(fFlagLS) fInvmassLSPt->Fill(track->Pt(),mass,fWeight);
    if(fFlagULS) fInvmassULSPt->Fill(track->Pt(),mass,fWeight);

    if(fFlagLS && mass<fInvmassCut) {
      //ElectronHadCorrel(itrack, track, fSprsLSEHCorrl, fHisLSDphi);
      ElectronHadCorrel(itrack, track, fSprsLSEHCorrl);
      ElectronHadCorrelNoPartner(itrack, jTracks, track, fSprsLSNoPartnerEHCorrl);
      MixedEvent(track, fSprsMixLSEHCorrl);
      fLSElecPt->Fill(track->Pt(),fWeight);
    }
    if(fFlagULS && mass<fInvmassCut){
      //ElectronHadCorrel(itrack, track, fSprsULSEHCorrl, fHisULSDphi);
      ElectronHadCorrel(itrack, track, fSprsULSEHCorrl);
      ElectronHadCorrelNoPartner(itrack, jTracks, track, fSprsULSNoPartnerEHCorrl);
      MixedEvent(track, fSprsMixULSEHCorrl);
      fULSElecPt->Fill(track->Pt(),fWeight);
    }

    if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec){
      flagPhotonicElec = kTRUE;
    }
    if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec){
      flagLSElec = kTRUE;
    }
  }
  fFlagPhotonicElec = flagPhotonicElec;
  fFlagElecLS = flagLSElec;

}
//___________________________________________
void  AliAnalysisTaskEHCorrel::MixedEvent(AliVTrack *track, THnSparse *SparseMixEHCorrl)
{
  //Retrive mixed event pool
  Double_t zVtx;
  zVtx = fpVtx->GetZ();

  Double_t fvalueMixDphi[6] = {-999,999,-999,-999,-999,-999}; //ptElec, ptHad,Dphi, Deta
  Double_t pi = TMath::Pi();

  AliEventPool* fPool;
  fPool = fPoolMgr->GetEventPool(fCentrality, zVtx); // Get the buffer associated with the current centrality and z-vtx
  if (!fPool)
  {
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCentrality, zVtx));
    return;
  }
  //fPool->PrintInfo();
  fMixStatCent->Fill(fPool->GetCurrentNEvents(),fCentrality);
  fMixStatVtxZ->Fill(fPool->GetCurrentNEvents(),zVtx);

  if (fPool->GetCurrentNEvents() >= 3) // start mixing when 3 events are in the buffer
  {
    Int_t nMix = fPool->GetCurrentNEvents();
    fNoMixedEvents->Fill(0);

    fMixStatCentVtxz->Fill(zVtx,fCentrality);

    Double_t ptMixTrk= -999;
    Double_t ptEle = -999;
    Double_t phiEle = -999, phiMixTrk = -999, Dphi = -999;
    Double_t etaEle = -999, etaMixTrk = -999, Deta = -999;

    // cout << "nMix = " << nMix << " tracks in pool = " << pool->NTracksInPool() << endl;
    for (Int_t jMix=0; jMix<nMix; jMix++)  // mix with each event in the buffer
    {
      TObjArray* bgTracks = fPool->GetEvent(jMix);
      for (Int_t i=0;i<bgTracks->GetEntriesFast(); i++)
      {
        AliVParticle* mixtrk = (AliVParticle*) bgTracks->At(i);

        if (!mixtrk) {
          printf("ERROR: Could not receive mix pool track %d\n",i);
          continue;
        }

        ptEle = track->Pt();
        phiEle = track->Phi();
        etaEle = track->Eta();
        phiMixTrk = mixtrk->Phi();
        ptMixTrk = mixtrk->Pt();
        etaMixTrk = mixtrk->Eta();

        if(fTrigElePtCut && (ptMixTrk > ptEle)) continue;

        Dphi = phiEle - phiMixTrk;
        if (Dphi > 3*pi/2)
          Dphi = Dphi - 2*pi;
        if (Dphi < -pi/2)
          Dphi = Dphi + 2*pi;

        Deta = etaEle - etaMixTrk;

        fvalueMixDphi[0] = ptEle;
        fvalueMixDphi[1] = ptMixTrk;
        fvalueMixDphi[2] = Dphi;
        fvalueMixDphi[3] = Deta;
        //fvalueMixDphi[4] = fVtxZBin;
        //fvalueMixDphi[5] = fCentBin;
        SparseMixEHCorrl->Fill(fvalueMixDphi,fWeight);
      }
    }
  }
}
//___________________________________________
TObjArray* AliAnalysisTaskEHCorrel::CloneAndReduceTrackList()
{
  // clones a track list by using AliehDPhiBasicParticle which uses much less memory (used for event mixing)
  TObjArray* fArrayTracksMix = new TObjArray;
  fArrayTracksMix->SetOwner(kTRUE);

  Int_t ntracks = -999;
  if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries();

  for(Int_t ktracks = 0; ktracks < ntracks; ktracks++){
    AliVParticle* Vtrack = 0x0;
    if(!fUseTender) Vtrack  =  fVevent->GetTrack(ktracks);
    if(fUseTender) Vtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(ktracks)); //take tracks from Tender list

    if (!Vtrack) {
      printf("ERROR: Could not receive track %d\n", ktracks);
      continue;
    }

    AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
    if(!track) continue;
    AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
    if(!atrack) continue;

    if(!PassHadronCuts(atrack)) continue; //apply hadron cuts;

    AliVParticle* particle = (AliVParticle*) fVevent->GetTrack(ktracks);
    fArrayTracksMix->Add(new AliehDPhiBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge()));
  }
  return fArrayTracksMix;
}
//___________________________________________
void AliAnalysisTaskEHCorrel::GetVtxZCentralityBin()
{
  //Get the VtxZ and centrality bin for THnSparse

  //Vtx Z
  Double_t pVtxZ = fpVtx->GetZ();
  if(pVtxZ > -10 && pVtxZ <= -5) fVtxZBin = 0;
  if(pVtxZ > -5 && pVtxZ <= -2) fVtxZBin = 1;
  if(pVtxZ > -2 && pVtxZ <= 0) fVtxZBin = 2;
  if(pVtxZ > 0 && pVtxZ <= 2) fVtxZBin = 3;
  if(pVtxZ > 2 && pVtxZ <= 5) fVtxZBin = 4;
  if(pVtxZ > 5 && pVtxZ <= 10) fVtxZBin = 5;

  //Centrality
  if(fCentralityMax == 20)
  {
    if(fCentrality > 0 && fCentrality <= 2) fCentBin = 0;
    if(fCentrality > 2 && fCentrality <= 4) fCentBin = 1;
    if(fCentrality > 4 && fCentrality <= 6) fCentBin = 2;
    if(fCentrality > 6 && fCentrality <= 10) fCentBin = 3;
    if(fCentrality > 10 && fCentrality <= 15) fCentBin = 4;
    if(fCentrality > 15 && fCentrality <= 20) fCentBin = 5;
  }
  if(fCentralityMax == 50)
  {
    if(fCentrality > 20 && fCentrality <= 25) fCentBin = 0;
    if(fCentrality > 25 && fCentrality <= 30) fCentBin = 1;
    if(fCentrality > 30 && fCentrality <= 35) fCentBin = 2;
    if(fCentrality > 35 && fCentrality <= 40) fCentBin = 3;
    if(fCentrality > 40 && fCentrality <= 45) fCentBin = 4;
    if(fCentrality > 45 && fCentrality <= 50) fCentBin = 5;
  }
  if(fCentralityMax > 50)
  {
    if(fCentrality > 50 && fCentrality <= 55) fCentBin = 0;
    if(fCentrality > 55 && fCentrality <= 60) fCentBin = 1;
    if(fCentrality > 60 && fCentrality <= 65) fCentBin = 2;
    if(fCentrality > 65 && fCentrality <= 70) fCentBin = 3;
    if(fCentrality > 70 && fCentrality <= 75) fCentBin = 4;
    if(fCentrality > 75 && fCentrality <= 80) fCentBin = 5;
  }
}
//___________________________________________
void AliAnalysisTaskEHCorrel::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}



