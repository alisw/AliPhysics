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
//  Task for Hadron/LeadingParticle-HFE & HFE-HFE Correlation         //
//  Author: Florian Herrmann & Denise Godoy                           //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH3F.h"
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
#include  "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliVEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"

#include "AliAnalysisTaskHaHFECorrel.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"
#include "AliESDtrackCuts.h"
#include "AliPhysicsSelection.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"
#include "AliESDv0KineCuts.h"
#include "AliESDv0.h"
#include "AliEMCALTrack.h"
#include "AliMagF.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"

#include "AliEventplane.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"

#include "AliEventPoolManager.h"


ClassImp(AliAnalysisTaskHaHFECorrel)
ClassImp(AliBasicParticleHaHFE)

using std::cout;
using std::endl;

//________________________________________________________________________
AliAnalysisTaskHaHFECorrel::AliAnalysisTaskHaHFECorrel(const char *name)
: AliAnalysisTaskSE(name)
,fUseTender(kFALSE)
,fWhichPeriod(2016)
,fUseKFforPhotonicPartner(kFALSE)
,fMaxPtEvent(999)
,fMinPtEvent(0)
,fMaxElectronEta(0.8)
,fMinElectronEta(-0.8)
,fMaxHadronEta(0.8)
,fMinHadronEta(-0.8)
,fTPCnCut(100)
,fTPCndEdxCut(80)
,fITSnCut(3)
,fUseITS(0)
,fSigmaITScut(2.)
,fSigmaTOFcut(2.)
,fSigmaTPCcut(-1.)
,fPhotElecPtCut(0.125)
,fPhotElecSigmaTPCcut(3)
,fPhotElecTPCnCut(80)
,fPhotElecITSrefitCut(kTRUE)
,fAssNonEleTPCcut(-4)
,fHTPCnCut(100)
,fHITSrefitCut(kTRUE)
,fHTPCrefitCut(kTRUE)
,fOpeningAngleCut(1000.)
,fInvmassCut(0.14)
,fChi2Cut(3.5)
,fDCAcut(999)
,fTRDQA(kFALSE)
,fCorrHadron(kTRUE)
,fCorrLParticle(kTRUE)
,fMixedEvent(kTRUE)
,fLParticle(kTRUE)
,fESD(0)
,fesdTrackCuts(0)
,fAOD(0)
,fVevent(0)
,fpidResponse(0)
,fMultSelection(0)
,fCentrality(0)
,fPoolMgr(0)
,fMC(0)
,fStack(0)
,fMCparticle(0)
,fMCarray(0)
,fMCheader(0)
,fTracks_tender(0)
,fCaloClusters_tender(0)
,fEventCuts()
,fCuts(0)
,fIsMC(kFALSE)
,fIsAOD(kTRUE)
,fCFM(0)
,fPID(0)
,fPIDqa(0)
,fOutputList(0)
,fNoEvents(0)
,fTrkpt(0)
,fEtaVtxZ(0)
,fMultiplicity(0)
,fRunList(0)
,fElectronTrackCuts(0)
,fElectronTrackTPCNcls(0)
,fElectronTrackTPCNclsdEdx(0)
,fElectronTrackTPCFrac(0)
,fElectronTrackITSNcls(0)
,fElectronTrackDCA(0)
,fHistITSnSig(0)
,fHistTOFnSig(0)
,fHistTPCnSig(0)
,fHistTPCnSigITScut(0)
,fHistTPCnSigTOFcut(0)
,fHistTPCnSigITSTOFcut(0)
,fCheckNHadronScaling(0)
,fCheckNPhotHadScaling(0)
,fHadContPvsPt(0)
,fHadContEtaPhiPt(0)
,fHadContTPCEtaPhiPt(0)
,fHadContPPhiEtaTPC(0)
,fHadContamination(0)
,fHadContaminationPt(0)
,fHadContMC(0)
,fHadContMCPt(0)
,fInclElecPt(0)
,fInclElecP(0)
,fULSElecPt(0)
,fLSElecPt(0)
,fInvmassLS(0)
,fInvmassULS(0)
,fOpeningAngleLS(0)
,fOpeningAngleULS(0)
,fCheckLSULS(0)
,fTagEtaPhiPt(0)
,fNonTagEtaPhiPt(0)
,fTagMotherPt(0)
,fTagEffIncl(0)
,fTagEffLS(0)
,fTagEffULS(0)
,fTagTruePairs(0)
,fTagEffInclWoWeight(0)
,fTagEffLSWoWeight(0)
,fTagEffULSWoWeight(0)
,fTagTruePairsWoWeight(0)
,fCorrectPiontoData()
,fCorrectEtatoData()
,fAssPtHad_Nbins(0)
,fAssPtHad_Xmin()
,fAssPtHad_Xmax()
,fAssPtElec_Nbins(0)
,fAssPtElec_Xmin()
,fAssPtElec_Xmax()
,fElecTrigger(0)
,fInclElecPhi(0)
,fInclElecEta(0)
,fULSElecPhi(0)
,fLSElecPhi(0)
,fElecDphi(0)
,fULSElecDphi(0)
,fLSElecDphi(0)
,fULSElecDphiDiffMethod(0)
,fLSElecDphiDiffMethod(0)
,fNoPartnerNoT(0)
,fTPartnerNoT(0)
,fElecHadTrigger(0)
,fElecHadTriggerLS(0)
,fElecHadTriggerULS(0)
,fHadContTrigger(0)
,fHadElecTrigger(0)
,fNonElecHadTrigger(0)
,fHadNonElecTrigger(0)
,fInclElecHa(0)
,fLSElecHa(0)
,fULSElecHa(0)
,fMCElecHaHadron(0)
,fElecHaHa(0)
,fElecHaLSNoPartner(0)
,fElecHaULSNoPartner(0)
,fMCElecHaTruePartner(0)
,fMCElecHaNoPartner(0)
,fMCElecHaRemovedPartner(0)
,fMCElecHaTruePartnerTrigger(0)
,fMCElecHaNoPartnerTrigger(0)
,fMCElecHaRemovedPartnerTrigger(0)
,fElecHaMixedEvent(0)
,fLSElecHaMixedEvent(0)
,fULSElecHaMixedEvent(0)
,fTagHaMixedEvent(0)
,fNonTagHaMixedEvent(0)
,fElecLPTrigger(0)
,fElecLPTriggerLS(0)
,fElecLPTriggerULS(0)
,fHadContLPTrigger(0)
,fLPElecTrigger(0)
,fLPNonElecTrigger(0)
,fNonElecLPTrigger(0)
,fInclElecLP(0)
,fLSElecLP(0)
,fULSElecLP(0)
,fMCElecLPHadron(0)
,fElecLPHa(0)
,fElecLPLSNoPartner(0)
,fElecLPULSNoPartner(0)
,fMCElecLPTruePartner(0)
,fMCElecLPNoPartner(0)
,fMCElecLPRemovedPartner(0)
,fMCElecLPTruePartnerTrigger(0)
,fMCElecLPNoPartnerTrigger(0)
,fMCElecLPRemovedPartnerTrigger(0)
,fElecLPMixedEvent(0)
,fLSElecLPMixedEvent(0)
,fULSElecLPMixedEvent(0)
,fCheckMCVertex(0)
,fCheckMCPtvsRecPtHad(0)
,fCheckMCEtavsRecEtaHad(0)
,fCheckMCPhivsRecPhiHad(0)
,fMCHadPtEtaPhiVtx(0)
,fRecHadMCPtEtaPhiVtx(0)
,fRecHadPtEtaPhiVtx(0)
,fCheckMCPtvsRecPtEle(0)
,fMCElecPtEtaPhiVtx(0)
,fRecElecPtEtaPhiVtx(0)
,fRecElecMCPtEtaPhiVtx(0)
,fMCElecPDG(0)
,fMCElecPtEtaPhiStrictVtx(0)
,fMCPi0Prod(0)
,fMCEtaProd(0)
,fMCPiPlusProd(0)
,fMCPiPlusProdV2(0)
,fMCLeadingParticle(0)
,fV0cutsESD(0)
,fV0cuts(0)
,fV0electrons(0)
,fV0pions(0)
,fV0protons(0)
,fhArmenteros(0)
,fEventsPerRun(0)
,fTRDnTrackRun(0)
,fV0tags(0)
,fTRDEtaPhi(0)
,fTRDNTracklets(0)
,fTRDV0NTracklets(0)
,fTRDSpectra(0)
,fTRDV0Spectra(0)
,fTRDMCSpectra(0)
{
    //Named constructor
   
    DefineInput(0, TChain::Class());
    //DefineInput(1, TList::Class());
    DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskHaHFECorrel::AliAnalysisTaskHaHFECorrel()
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisElecHadCorrel")
,fUseTender(kFALSE)
,fWhichPeriod(2016)
,fUseKFforPhotonicPartner(kFALSE)
,fMaxPtEvent(999)
,fMinPtEvent(0)
,fMaxElectronEta(0.8)
,fMinElectronEta(-0.8)
,fMaxHadronEta(0.8)
,fMinHadronEta(-0.8)
,fTPCnCut(100)
,fTPCndEdxCut(80)
,fITSnCut(3)
,fUseITS(0)
,fSigmaITScut(2.)
,fSigmaTOFcut(2.)
,fSigmaTPCcut(-1.)
,fPhotElecPtCut(0.125)
,fPhotElecSigmaTPCcut(3)
,fPhotElecTPCnCut(80)
,fPhotElecITSrefitCut(kTRUE)
,fAssNonEleTPCcut(-4)
,fHTPCnCut(100)
,fHITSrefitCut(kTRUE)
,fHTPCrefitCut(kTRUE)
,fOpeningAngleCut(1000.)
,fInvmassCut(0.14)
,fChi2Cut(3.5)
,fDCAcut(999)
,fTRDQA(kFALSE)
,fCorrHadron(kTRUE)
,fCorrLParticle(kTRUE)
,fMixedEvent(kTRUE)
,fLParticle(kTRUE)
,fESD(0)
,fesdTrackCuts(0)
,fAOD(0)
,fVevent(0)
,fpidResponse(0)
,fMultSelection(0)
,fCentrality(0)
,fPoolMgr(0)
,fMC(0)
,fStack(0)
,fMCparticle(0)
,fMCarray(0)
,fMCheader(0)
,fTracks_tender(0)
,fCaloClusters_tender(0)
,fEventCuts()
,fCuts(0)
,fIsMC(kFALSE)
,fIsAOD(kTRUE)
,fCFM(0)
,fPID(0)
,fPIDqa(0)
,fOutputList(0)
,fNoEvents(0)
,fTrkpt(0)
,fEtaVtxZ(0)
,fMultiplicity(0)
,fRunList(0)
,fElectronTrackCuts(0)
,fElectronTrackTPCNcls(0)
,fElectronTrackTPCNclsdEdx(0)
,fElectronTrackTPCFrac(0)
,fElectronTrackITSNcls(0)
,fElectronTrackDCA(0)
,fHistITSnSig(0)
,fHistTOFnSig(0)
,fHistTPCnSig(0)
,fHistTPCnSigITScut(0)
,fHistTPCnSigTOFcut(0)
,fHistTPCnSigITSTOFcut(0)
,fCheckNHadronScaling(0)
,fCheckNPhotHadScaling(0)
,fHadContPvsPt(0)
,fHadContEtaPhiPt(0)
,fHadContTPCEtaPhiPt(0)
,fHadContPPhiEtaTPC(0)
,fHadContamination(0)
,fHadContaminationPt(0)
,fHadContMC(0)
,fHadContMCPt(0)
,fInclElecPt(0)
,fInclElecP(0)
,fULSElecPt(0)
,fLSElecPt(0)
,fInvmassLS(0)
,fInvmassULS(0)
,fOpeningAngleLS(0)
,fOpeningAngleULS(0)
,fCheckLSULS(0)
,fTagEtaPhiPt(0)
,fNonTagEtaPhiPt(0)
,fTagMotherPt(0)
,fTagEffIncl(0)
,fTagEffLS(0)
,fTagEffULS(0)
,fTagTruePairs(0)
,fTagEffInclWoWeight(0)
,fTagEffLSWoWeight(0)
,fTagEffULSWoWeight(0)
,fTagTruePairsWoWeight(0)
,fCorrectPiontoData()
,fCorrectEtatoData()
,fAssPtHad_Nbins(0)
,fAssPtHad_Xmin()
,fAssPtHad_Xmax()
,fAssPtElec_Nbins(0)
,fAssPtElec_Xmin()
,fAssPtElec_Xmax()
,fElecTrigger(0)
,fInclElecPhi(0)
,fInclElecEta(0)
,fULSElecPhi(0)
,fLSElecPhi(0)
,fElecDphi(0)
,fULSElecDphi(0)
,fLSElecDphi(0)
,fULSElecDphiDiffMethod(0)
,fLSElecDphiDiffMethod(0)
,fNoPartnerNoT(0)
,fTPartnerNoT(0)
,fElecHadTrigger(0)
,fElecHadTriggerLS(0)
,fElecHadTriggerULS(0)
,fHadContTrigger(0)
,fHadElecTrigger(0)
,fNonElecHadTrigger(0)
,fHadNonElecTrigger(0)
,fInclElecHa(0)
,fLSElecHa(0)
,fULSElecHa(0)
,fMCElecHaHadron(0)
,fElecHaHa(0)
,fElecHaLSNoPartner(0)
,fElecHaULSNoPartner(0)
,fMCElecHaTruePartner(0)
,fMCElecHaNoPartner(0)
,fMCElecHaRemovedPartner(0)
,fMCElecHaTruePartnerTrigger(0)
,fMCElecHaNoPartnerTrigger(0)
,fMCElecHaRemovedPartnerTrigger(0)
,fElecHaMixedEvent(0)
,fLSElecHaMixedEvent(0)
,fULSElecHaMixedEvent(0)
,fTagHaMixedEvent(0)
,fNonTagHaMixedEvent(0)
,fElecLPTrigger(0)
,fElecLPTriggerLS(0)
,fElecLPTriggerULS(0)
,fHadContLPTrigger(0)
,fLPElecTrigger(0)
,fLPNonElecTrigger(0)
,fNonElecLPTrigger(0)
,fInclElecLP(0)
,fLSElecLP(0)
,fULSElecLP(0)
,fMCElecLPHadron(0)
,fElecLPHa(0)
,fElecLPLSNoPartner(0)
,fElecLPULSNoPartner(0)
,fMCElecLPTruePartner(0)
,fMCElecLPNoPartner(0)
,fMCElecLPRemovedPartner(0)
,fMCElecLPTruePartnerTrigger(0)
,fMCElecLPNoPartnerTrigger(0)
,fMCElecLPRemovedPartnerTrigger(0)
,fElecLPMixedEvent(0)
,fLSElecLPMixedEvent(0)
,fULSElecLPMixedEvent(0)
,fCheckMCVertex(0)
,fCheckMCPtvsRecPtHad(0)
,fCheckMCEtavsRecEtaHad(0)
,fCheckMCPhivsRecPhiHad(0)
,fMCHadPtEtaPhiVtx(0)
,fRecHadMCPtEtaPhiVtx(0)
,fRecHadPtEtaPhiVtx(0)
,fCheckMCPtvsRecPtEle(0)
,fMCElecPtEtaPhiVtx(0)
,fRecElecPtEtaPhiVtx(0)
,fRecElecMCPtEtaPhiVtx(0)
,fMCElecPDG(0)
,fMCElecPtEtaPhiStrictVtx(0)
,fMCPi0Prod(0)
,fMCEtaProd(0)
,fMCPiPlusProd(0)
,fMCPiPlusProdV2(0)
,fMCLeadingParticle(0)
,fV0cutsESD(0)
,fV0cuts(0)
,fV0electrons(0)
,fV0pions(0)
,fV0protons(0)
,fhArmenteros(0)
,fEventsPerRun(0)
,fTRDnTrackRun(0)
,fV0tags(0)
,fTRDEtaPhi(0)
,fTRDNTracklets(0)
,fTRDV0NTracklets(0)
,fTRDSpectra(0)
,fTRDV0Spectra(0)
,fTRDMCSpectra(0)
{

    // Default constructor
    
    DefineInput(0, TChain::Class());
    //DefineInput(1, TList::Class());
    DefineOutput(1, TList::Class());
}
//_________________________________________

AliAnalysisTaskHaHFECorrel::~AliAnalysisTaskHaHFECorrel()
{
  //Destructor
  delete fOutputList;
  delete fPID;
  delete fCFM;
  delete fCuts;
  delete fPIDqa;
  delete fTracks_tender;
  delete fCaloClusters_tender;
  delete fPoolMgr;
  delete fV0cuts;
  delete fV0cutsESD;
  //delete fTrackCuts;
  //delete fAssTrackCuts;
}
//_________________________________________

void AliAnalysisTaskHaHFECorrel::UserExec(Option_t*)
{
  //Main loop
  //Called for each event 



  // create pointer to event
  fVevent = dynamic_cast<AliVEvent*>(InputEvent());
  if(!fVevent){
    printf("ERROR: fVEvent not available\n");
    return;
  }
  
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!(fESD || fAOD))
    {
      printf("ERROR: fESD & fAOD not available\n");
      return;
    }
  
  if(fUseTender){
    fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
    fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
    if (!fTracks_tender || !fCaloClusters_tender) return;
  }
 
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();

  if(fIsMC) {
    if (fIsAOD) {
      AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
      fMC = eventHandler->MCEvent();
      fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
      if (!fMCheader) AliError("fMCheader could not be initialised");
    }
    else {
      AliMCEventHandler *eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      fMC = eventHandler->MCEvent();
      fStack = fMC->Stack();
    }
    if (!fMC) {
	cout << "No fMC" << endl;
    }
  }

  // suggested by DPG to remove outliers
  const AliVVertex *pVtx=0;  
  const AliVVertex *spdVtx=0;
  if (fAOD) {
    pVtx =   fAOD->GetPrimaryVertex();
    spdVtx = fAOD->GetPrimaryVertexSPD();
  }
  else if (fESD) {
    pVtx =   fESD->GetPrimaryVertex();
    spdVtx = fESD->GetPrimaryVertexSPD();
  }
    
  Double_t pVtxZ = -999.;
  pVtxZ = pVtx->GetZ();
  if (TMath::Abs(pVtxZ)>10.) return;
   
  fNoEvents->Fill(0);
  
  Int_t fNOtrks = fVevent->GetNumberOfTracks();
  if(fNOtrks<2) return;
  fNoEvents->Fill(1);
    
  if(!fEventCuts.AcceptEvent(fVevent)) {
    PostData(1, fOutputList);
    return;
  }
  fNoEvents->Fill(2);

  fpidResponse = fInputHandler->GetPIDResponse();
  if(!fpidResponse){
    AliDebug(1, "Using default PID Response");
    fpidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class());
  }  
  

  Int_t *listofmotherkink=0;
  Int_t nMotherKink=0;
  if (fIsAOD) {
    Int_t nVertices = 1;
    nVertices = fAOD->GetNumberOfVertices();
    listofmotherkink = new Int_t[nVertices];
      for(Int_t ivertex=0; ivertex < nVertices; ivertex++) {
      AliAODVertex *vertex = fAOD->GetVertex(ivertex);
      if(!vertex) continue;
      if(vertex->GetType()==AliAODVertex::kKink) {
	AliAODTrack *mother = (AliAODTrack *) vertex->GetParent();
	if(!mother) continue;
	Int_t idmother = mother->GetID();
	listofmotherkink[nMotherKink] = idmother;
	nMotherKink++;
      }
    }
  }

  if (fMinPtEvent > 0.1 || fMaxPtEvent <100) {
    if (!PassEventBias(pVtx,nMotherKink,listofmotherkink)) {
      delete listofmotherkink;
      return;
    }
  }
  
  fNoEvents->Fill(3);


  Double_t fMultV0Per, fMultSPDPer, fMultV0Tot, fMultSPD;
  Float_t mult = 1;
  // if (!fIsMC) {
    fMultV0Per = -1.;
    fMultSPDPer =-1.;
    // fAOD->GetMultiplicity();
    fMultSelection = (AliMultSelection * ) fVevent->FindListObject("MultSelection");
    fMultV0Per = fMultSelection->GetMultiplicityPercentile("V0M", kFALSE); // Method, Embed Event selection (kFALSE)
    fMultSPDPer = fMultSelection->GetMultiplicityPercentile("SPDTracklets", kFALSE); // Method, Embed Event selection (kFALSE)
    
    // TRY Multiplicity estimates
    AliVVZERO* AODV0 = fVevent->GetVZEROData();
    Float_t multV0A=AODV0->GetMTotV0A();
    Float_t multV0C=AODV0->GetMTotV0C();
    fMultV0Tot=multV0A+multV0C;
    
    AliVMultiplicity* AliMult = fVevent->GetMultiplicity();
    fMultSPD=AliMult->GetNumberOfTracklets();
    
    //  AliAODTracklets *AODtracklets = ((AliAODEvent*)fAOD)->GetTracklets();
    //  fMultAOD = AODtracklets->GetNumberOfTracklets();

    mult=fMultSPDPer;

    Double_t fillSparse[4]={fMultV0Per, fMultSPDPer, fMultV0Tot, fMultSPD};
    fMultiplicity->Fill(fillSparse);
    //}
     
  if(fIsMC) {
    if (fIsAOD) {
      MCEfficiencyCorrections(pVtx); //  Electron reconstruction, Hadron reconstruction
      TList *lh=fMCheader->GetCocktailHeaders();
      Int_t nh=lh->GetEntries();  
      for(Int_t i=0;i<nh;i++)
	{
	  //AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
	  // TString genname=gh->GetName();
	  // Int_t npart=gh->NProduced();
	  //	cout << genname << "with " << npart << " Particles" << endl;
	}
    }
    else {
      //  MCEfficiencyCorrections(pVtx); //  Electron reconstruction, Hadron reconstruction
    }
  }


  ///////////////////////
  // Preparational Tasks
  ///////////////////////


  // List of HFE for analysis and mixed event
  TObjArray* RedTracksHFE = new TObjArray;
  RedTracksHFE->SetOwner(kTRUE);

  // FirstLoop: Find leading particle (LP) and HFE (check for LS/ULS partner)
  AliVTrack* LPtrack;
  fLParticle=kFALSE;
  LPtrack=FindLPAndHFE(RedTracksHFE, pVtx,nMotherKink,listofmotherkink);


  if(fTRDQA) {
    Int_t RunNumber = fVevent->GetRunNumber();
    Int_t RunIndex=0;
    while (RunNumber!=fRunList[RunIndex]) {
      RunIndex++;
      if (RunIndex>199) {
	cout << "Run not found in run list" << endl;
	break;
      }
    }
    if (fIsAOD)  FindV0CandidatesAOD(fAOD);
    else FindV0CandidatesESD(fESD);
    fEventsPerRun->Fill(RunIndex);
    TRDQA(RunIndex,pVtx,nMotherKink,listofmotherkink);
  }

  //////////////////////////
  // Main analysis
  ///////////////////////////

  if (fCorrLParticle) {
    // LP - two different functions for Same Event and mixed Event
    if (fLParticle && RedTracksHFE->GetEntriesFast()>0) CorrelateLP(LPtrack,pVtx, nMotherKink, listofmotherkink, RedTracksHFE);
    if (fLParticle && fMixedEvent  && RedTracksHFE->GetEntriesFast()>0) CorrelateLPMixedEvent(LPtrack, mult, pVtx->GetZ(), LPtrack->Pt());
  }

  // Hadron - only one function for both options, as only one loop over Hadron tracks
  // Mixed event is called inside this function
  if (fCorrHadron && fLParticle &&  RedTracksHFE->GetEntriesFast()>0) {
    CorrelateHadron(RedTracksHFE, pVtx, nMotherKink, listofmotherkink, mult, LPtrack->Pt());
  }
  // Electrons - currently no mixed event (>1 as two electrons are required), work in progress
  // if (RedTracksHFE->GetEntriesFast()>1) CorrelateElectron(RedTracksHFE);
 

  /////////////////////////
  //Fill Mixed event pool//
  /////////////////////////
  if (RedTracksHFE->GetEntriesFast()>0 && fLParticle) {
 
    AliEventPool * HFEPool = fPoolMgr->GetEventPool(mult, pVtxZ, LPtrack->Pt()); // Get the buffer associated with the current centrality and z-vtx
   HFEPool->SetDebug(kTRUE);
   cout << Form("Pool found for centrality = %f, zVtx = %f, leading pt = %f", mult, pVtx->GetZ(), LPtrack->Pt()) << endl;

    if (!HFEPool)
      {
	AliFatal(Form("No pool found for centrality = %f, zVtx = %f, leading pt = %f", mult, pVtx->GetZ(), LPtrack->Pt()));
	return;
      }
    if (RedTracksHFE->GetEntriesFast()>0) {
      HFEPool->UpdatePool(RedTracksHFE);
    }
  }

  ClearV0PIDList();
  delete listofmotherkink;

  PostData(1, fOutputList);
}

//_________________________________________
void AliAnalysisTaskHaHFECorrel::UserCreateOutputObjects()
{      
  
  AliLog::SetClassDebugLevel("AliAnalysisTaskHaHFECorrel", 0);
  AliDebug(1, Form("This is AliDebug level %i", 1));
  AliDebug(2, Form("This is AliDebug level %i", 2));
  AliDebug(3, Form("This is AliDebug level %i", 3));
   
  //---------Output Tlist
  fOutputList = new TList();
  fOutputList->SetOwner();

  fEventCuts.AddQAplotsToList(fOutputList);
  fesdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 

  // V0 Kine cuts 
  fV0cuts = new AliAODv0KineCuts();
  fV0cutsESD = new AliESDv0KineCuts();

  Int_t NumberOfRuns = 192;
  Int_t Runs[192] = {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256591, 256589, 256567, 256565, 256564, 256562, 256560, 256557, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504};
  fRunList = new Int_t[NumberOfRuns];
  for (Int_t i=0; i<NumberOfRuns; i++) {
    fRunList[i]=Runs[i];
  }



  // V0 PID Obj arrays
  fV0electrons = new TObjArray;
  fV0pions     = new TObjArray;
  fV0protons   = new TObjArray;

  fNoEvents = new TH1F("fNoEvents","",4,0,4);
  fNoEvents->Sumw2();
  fOutputList->Add(fNoEvents);
    
  fTrkpt = new TH2F("fTrkpt","track pt",200,0,20,3,0,3);
  fTrkpt->Sumw2();
  fOutputList->Add(fTrkpt);

  fEtaVtxZ = new TH2F("fEtaVtxZ", "Eta vs VtxZ", 90, -0.9, 0.9, 100, -10, 10);
  fEtaVtxZ->Sumw2();
  fOutputList->Add(fEtaVtxZ);
    
    
  // Multiplicity Sparse
  Int_t    nBinsMult[4]={20, 20, 100, 100};
  Double_t xminMult[4]={0,0,0,0};
  Double_t xmaxMult[4]={100,100, 3000, 300};
  fMultiplicity = new THnSparseF("fMultiplicity", "Multiplicity: VOPerc, SPDPerc, TotV0, SPD", 4, nBinsMult, xminMult, xmaxMult);
  fMultiplicity->Sumw2();
  fOutputList->Add(fMultiplicity);


  // Track Cuts
  const char*  ElecTrackCutsLabels[15]= {"All", "FilterBit", "PtCut", "EtaCut", "TPCNcls", "TPCNclsdEdx", "TPCrefit", "TPCfrac", "KinkCut", "ITSNcls", "ITSrefit", "SPDAny", "SPDBoth", "DCACut", ""};
  fElectronTrackCuts = new TH2F("fElectronTrackCuts", "fElectronTrackCuts", 20, 0, 10, 15, 0,15);
  for (Int_t i=0; i<15; i++) fElectronTrackCuts->GetYaxis()->SetBinLabel(i+1, ElecTrackCutsLabels[i]);
  fElectronTrackCuts->Sumw2();
  fOutputList->Add(fElectronTrackCuts);

  fElectronTrackTPCNcls = new TH2F("fElectronTrackTPCNcls", "fElectronTrackTPCNcls", 20, 0, 10,170, 0, 170);
  fElectronTrackTPCNcls->Sumw2();
  fOutputList->Add(fElectronTrackTPCNcls);

  fElectronTrackTPCNclsdEdx= new TH2F("fElectronTrackTPCNclsdEdx", "fElectronTrackTPCNclsdEdx",20, 0, 10, 170, 0, 170);
  fElectronTrackTPCNclsdEdx->Sumw2();
  fOutputList->Add(fElectronTrackTPCNclsdEdx);

  fElectronTrackTPCFrac = new TH2F("fElectronTrackTPCFrac", "fElectronTrackTPCFrac", 20, 0, 10,110, 0, 1.1);
  fElectronTrackTPCFrac->Sumw2();
  fOutputList->Add(fElectronTrackTPCFrac);
  
  fElectronTrackITSNcls = new TH2F("fElectronTrackITSNcls", "fElectronTrackITSNcls", 20, 0, 10,10, 0, 10);
  fElectronTrackITSNcls->Sumw2();
  fOutputList->Add(fElectronTrackITSNcls);

  fElectronTrackDCA = new TH2F("fElectronTrackDCA", "fElectronTrackDCA", 30, 0., 3.0, 30, 0.0, 3.0);
  fElectronTrackDCA->Sumw2();
  fOutputList->Add(fElectronTrackDCA);


  fHistITSnSig = new TH2F("fHistITSnSig","fHistITSnSig",50,0.3,5,200,-10,10);
  fHistITSnSig->Sumw2();
  BinLogX(fHistITSnSig->GetXaxis());
  fOutputList->Add(fHistITSnSig);
  
  fHistTOFnSig = new TH2F("fHistTOFnSig","fHistTOFnSig",100,0.3,10,200,-10,10);
  fHistTOFnSig->Sumw2();
  BinLogX(fHistTOFnSig->GetXaxis());
  fOutputList->Add(fHistTOFnSig);
  
  fHistTPCnSig = new TH2F("fHistTPCnSig","fHistTPCnSig",150,0.3,15,200,-10,10);
  fHistTPCnSig->Sumw2();
  BinLogX(fHistTPCnSig->GetXaxis());
  fOutputList->Add(fHistTPCnSig);
    
  fHistTPCnSigITScut = new TH2F("fHistTPCnSigITScut","fHistTPCnSigITScut",150,0.3,15,200,-10,10);
  fHistTPCnSigITScut->Sumw2();
  BinLogX(fHistTPCnSigITScut->GetXaxis());
  fOutputList->Add(fHistTPCnSigITScut);
    
  fHistTPCnSigTOFcut = new TH2F("fHistTPCnSigTOFcut","fHistTPCnSigTOFcut",150,0.3,15,200,-10,10);
  fHistTPCnSigTOFcut->Sumw2();
  BinLogX(fHistTPCnSigTOFcut->GetXaxis());
  fOutputList->Add(fHistTPCnSigTOFcut);
    
  fHistTPCnSigITSTOFcut = new TH2F("fHistTPCnSigITSTOFcut","fHistTPCnSigITSTOFcut",150,0.3,15,200,-10,10);
  fHistTPCnSigITSTOFcut->Sumw2();
  BinLogX(fHistTPCnSigITSTOFcut->GetXaxis());
  fOutputList->Add(fHistTPCnSigITSTOFcut);


  // General  Binning 		    
  // fAssPtHad_Nbins = 15;
  const Int_t TmpHad_Nbins = 15;
  fAssPtHad_Nbins = TmpHad_Nbins;
  Float_t TmpArrayLow[15]= {0,  0.5, 0.5, 0.5, 1., 1., 1., 2,  2, 5, 5,  5, 10, 10, 15};
  Float_t TmpArrayUp[15]={ 999, 2.0, 5.0, 999, 2., 5.,999, 5,999,10,15,999, 15,999,999};
  fAssPtHad_Xmin.Set(TmpHad_Nbins, TmpArrayLow);
  fAssPtHad_Xmax.Set(TmpHad_Nbins, TmpArrayUp);
  
  const Int_t TmpElec_Nbins =12;
  fAssPtElec_Nbins = TmpElec_Nbins;
  Float_t TmpArrayELow[12]={0 ,  0.5, 0.5, 1., 1., 1., 2,  2, 4,  4,  6,  8};
  Float_t TmpArrayEUp[12]={999, 2.0, 999, 2., 4.,999, 4,999, 6,999, 10,999};
  fAssPtElec_Xmin.Set(TmpElec_Nbins, TmpArrayELow);
  fAssPtElec_Xmax.Set(TmpElec_Nbins, TmpArrayEUp);
  
  // HFE Electron
  Int_t    NBinsElectron =46;
  Double_t XminElectron=0.25;
  Double_t XmaxElectron=6.;
  //Int_t    NBinsElectronRed = 7;
  //Double_t XBinsElectronRed[8]={0.25,0.5, 1., 1.5, 3, 4, 6, 10};

  Int_t     NBinsHadron=200 ;
  Double_t  XminHadron=0.0;
  Double_t  XmaxHadron=100;
  //Int_t     NBinsHadRed=13;
  //Double_t  XBinsHaRed[14]={0.5,1.0, 1.5, 2.0, 2.5, 3.,3.5, 5., 7.5, 15., 30., 50., 100, 1000}; //13 bins

  Int_t    NBinsPhi=64;
  Int_t    NBinsEta=90;
  Double_t XminEta=-0.9;
  Double_t XmaxEta=0.9;
  Int_t    NBinsDPhi=64;
  Double_t DEtaMax=fMaxElectronEta+fMaxHadronEta;
  Double_t DEtaMin=fMinElectronEta+fMinHadronEta;
  Int_t    NBinsDEta=(DEtaMax-DEtaMin)/0.02;


  Int_t    binHadScaling[4]={201,201, 201, 201};
  Double_t xminHadScaling[4]={-0.5,-0.5,-0.5,-0.5};
  Double_t xmaxHadScaling[4]={200.5,200.5, 200.5, 200.5};


  fCheckNHadronScaling = new THnSparseF("fCheckNHadronScaling", "NHadScaling: NHadron, NElectron, NNonElectron, HWrongElectron", 4, binHadScaling, xminHadScaling, xmaxHadScaling);
  fCheckNHadronScaling->Sumw2();
  fOutputList->Add(fCheckNHadronScaling);



  fCheckNPhotHadScaling = new THnSparseF("fCheckNPhotHadScaling", "NHadScaling: NHadron, NElectron, NTagged, NNotTagged", 4, binHadScaling, xminHadScaling, xmaxHadScaling);
  fCheckNPhotHadScaling->Sumw2();
  fOutputList->Add(fCheckNPhotHadScaling);


  // HadContSparse
  Int_t    binHC[4] =  {NBinsElectron  ,NBinsPhi      ,NBinsEta    ,200}; //p, Phi, Eta, TPC
  Double_t xminHC[4] = {XminElectron   ,0             ,-0.9  ,-10};
  Double_t xmaxHC[4] = {XmaxElectron   ,TMath::TwoPi(), 0.9  ,10};

  fHadContPvsPt = new TH2F("fHadContPvsPt", "P vs Pt", 100, 0, 10, 100, 0, 10);
  fHadContPvsPt->Sumw2();
  fOutputList->Add(fHadContPvsPt);  
    
  fHadContPPhiEtaTPC = new THnSparseF("fHadContPPhiEtaTPC", "HadCont: P, Phi, Eta, TPC", 4, binHC, xminHC, xmaxHC);
  fHadContPPhiEtaTPC->Sumw2();
  fOutputList->Add(fHadContPPhiEtaTPC);
 
  Int_t    binHC2[4] =  {NBinsElectron  ,20   , 20 ,200}; //p, ITS, TOF, TPC
  Double_t xminHC2[4] = {XminElectron   ,-10  ,-10 ,-10};
  Double_t xmaxHC2[4] = {XmaxElectron   ,10   , 10 ,10};  
  fHadContamination = new THnSparseF("fHadContamination", "HadCont: P, ITS, TOF, TPC", 4, binHC2, xminHC2, xmaxHC2);
  fHadContamination->Sumw2();
  fOutputList->Add(fHadContamination);

  fHadContaminationPt = new THnSparseF("fHadContaminationPt", "HadCont: Pt, ITS, TOF, TPC", 4, binHC2, xminHC2, xmaxHC2);
  fHadContaminationPt->Sumw2();
  fOutputList->Add(fHadContaminationPt);

  fHadContTPCEtaPhiPt = new TH3F("fHadContTPCEtaPhiPt", "MCHadContTPC: Eta, Phi, Pt", 72, -0.9, 0.9, 64, 0, TMath::TwoPi(), NBinsElectron, XminElectron, XmaxElectron);
  fHadContTPCEtaPhiPt->Sumw2();
  fOutputList->Add(fHadContTPCEtaPhiPt);
  
  if (fIsMC) {
    fHadContEtaPhiPt = new TH3F("fHadContEtaPhiPt", "MCHadCont: Eta, Phi, Pt", 72, -0.9, 0.9, 64, 0, TMath::TwoPi(), NBinsElectron, XminElectron, XmaxElectron);
    fHadContEtaPhiPt->Sumw2();
    fOutputList->Add(fHadContEtaPhiPt);


    


    Int_t    binHMC[4] = {NBinsElectron  , 7,  20, 200}; // p PDG ITS, TPC
    Double_t xminHMC[4] = {XminElectron   , 0, -10, -10};
    Double_t xmaxHMC[4] = {XmaxElectron   , 7,  10,  10};
    fHadContMC = new THnSparseF("fHadContMC", "HadCont: P, PDG, ITS, TPC", 4, binHMC, xminHMC, xmaxHMC);
    fHadContMC->Sumw2();
    fOutputList->Add(fHadContMC);
    fHadContMCPt = new THnSparseF("fHadContMCPt", "HadCont: Pt, PDG, ITS, TPC", 4, binHMC, xminHMC, xmaxHMC);
    fHadContMCPt->Sumw2();
    fOutputList->Add(fHadContMCPt);

  }    
  
  
  fInclElecPt = new TH1F("fInclElePt", "fInclElePt",NBinsElectron, XminElectron, XmaxElectron);
  fInclElecPt->Sumw2();
  fOutputList->Add(fInclElecPt);

  fInclElecP = new TH1F("fInclEleP", "fInclEleP",NBinsElectron, XminElectron, XmaxElectron);
  fInclElecP->Sumw2();
  fOutputList->Add(fInclElecP);

  fULSElecPt = new TH1F("fULSElePt", "fULSElePt",NBinsElectron, XminElectron, XmaxElectron);
  fULSElecPt->Sumw2();
  fOutputList->Add(fULSElecPt);
    
  fLSElecPt = new TH1F("fLSElePt", "fLSElePt",NBinsElectron, XminElectron, XmaxElectron);
  fLSElecPt->Sumw2();
  fOutputList->Add(fLSElecPt);
    
  fInvmassULS = new TH2F("fInvmassULS", "fInvmassULS", 500,0,0.5,NBinsElectron,XminElectron,XmaxElectron);
  fInvmassULS->Sumw2();
  fOutputList->Add(fInvmassULS);
    
  fInvmassLS = new TH2F("fInvmassLS", "fInvmassLS", 500,0,0.5,NBinsElectron,XminElectron,XmaxElectron);
  fInvmassLS->Sumw2();
  fOutputList->Add(fInvmassLS);
  
  fOpeningAngleULS = new TH2F("fOpeningAngleULS","fOpeningAngleULS",100,0,1,NBinsElectron,XminElectron,XmaxElectron);
  fOpeningAngleULS->Sumw2();
  fOutputList->Add(fOpeningAngleULS);

  fOpeningAngleLS = new TH2F("fOpeningAngleLS","fOpeningAngleLS",100,0,1,NBinsElectron,XminElectron,XmaxElectron);
  fOpeningAngleLS->Sumw2();
  fOutputList->Add(fOpeningAngleLS);
    
  fCheckLSULS = new TH2F("fCheckLSULS", "LSULS",5,0,5,5,0,5);
  fCheckLSULS->Sumw2();
  fOutputList->Add(fCheckLSULS);
 

  fTagEtaPhiPt = new TH3F("fTagEtaPhiPt", "Tagged: Eta, Phi, Pt", 72, -0.9, 0.9, 64, 0, TMath::TwoPi(), NBinsElectron, XminElectron, XmaxElectron);
  fTagEtaPhiPt->Sumw2();
  fOutputList->Add(fTagEtaPhiPt);

  if (fIsMC) {
    fNonTagEtaPhiPt = new TH3F("fNonTagEtaPhiPt", "Tagged: Eta, Phi, Pt", 72, -0.9, 0.9, 64, 0, TMath::TwoPi(), NBinsElectron, XminElectron, XmaxElectron);
    fNonTagEtaPhiPt->Sumw2();
    fOutputList->Add(fNonTagEtaPhiPt);


    Int_t binMothPt[4]= {NBinsElectron, 50, 10001, 10001};
    Double_t xminMothPt[4]={XminElectron, 0, -0.5, -0.5};
    Double_t xmaxMothPt[4]={XmaxElectron, 25, 9999.5, 9999.5};
    
    fTagMotherPt = new THnSparseF("fTagMotherPt", "Incl: ptElectron, ptMother, Mother, Grandmother",4, binMothPt, xminMothPt, xmaxMothPt);
    fTagMotherPt->Sumw2();
    fOutputList->Add(fTagMotherPt);

    Int_t    binTagEff[5] =  {NBinsElectron   ,10000 ,10001, 10001, 10001}; //p, pdg, pdgmother
    Double_t xminTagEff[5] = {XminElectron  ,-0.5   , -1.5,-1.5, -1.5};
    Double_t xmaxTagEff[5] = {XmaxElectron    ,9999.5, 9999.5,9999.5, 9999.5};  

    fTagEffIncl = new THnSparseF("fTagEffIncl", "Incl tag.eff.: pt, pdg, pdgmother", 5, binTagEff, xminTagEff, xmaxTagEff);
    fTagEffIncl->Sumw2();
    fOutputList->Add(fTagEffIncl);
    
    fTagEffLS = new THnSparseF("fTagEffLS", "LS tag.eff.: pt, pdg, pdgmother", 5, binTagEff, xminTagEff, xmaxTagEff);
    fTagEffLS->Sumw2();
    fOutputList->Add(fTagEffLS);

    fTagEffULS = new THnSparseF("fTagEffULS", "ULS tag.eff.: pt, pdg, pdgmother", 5, binTagEff, xminTagEff, xmaxTagEff);
    fTagEffULS->Sumw2();
    fOutputList->Add(fTagEffULS);

    fTagTruePairs = new THnSparseF("fTagTruePairs", "ULS true pairs tag.eff.: pt, pdg, pdgmother", 5, binTagEff, xminTagEff, xmaxTagEff);
    fTagTruePairs->Sumw2();
    fOutputList->Add(fTagTruePairs);

    fTagEffInclWoWeight = new THnSparseF("fTagEffInclWoWeight", "Incl tag.eff.: pt, pdg, pdgmother", 5, binTagEff, xminTagEff, xmaxTagEff);
    fTagEffInclWoWeight->Sumw2();
    fOutputList->Add(fTagEffInclWoWeight);
    
    fTagEffLSWoWeight = new THnSparseF("fTagEffLSWoWeight", "LS tag.eff.: pt, pdg, pdgmother", 5, binTagEff, xminTagEff, xmaxTagEff);
    fTagEffLSWoWeight->Sumw2();
    fOutputList->Add(fTagEffLSWoWeight);

    fTagEffULSWoWeight = new THnSparseF("fTagEffULSWoWeight", "ULS tag.eff.: pt, pdg, pdgmother", 5, binTagEff, xminTagEff, xmaxTagEff);
    fTagEffULSWoWeight->Sumw2();
    fOutputList->Add(fTagEffULSWoWeight);

    fTagTruePairsWoWeight = new THnSparseF("fTagTruePairsWoWeight", "ULS true pairs tag.eff.: pt, pdg, pdgmother", 5, binTagEff, xminTagEff, xmaxTagEff);
    fTagTruePairsWoWeight->Sumw2();
    fOutputList->Add(fTagTruePairsWoWeight);

  }
  
  // HFE-HFE Correlation

  fElecTrigger = new TH1F("fEleTrigger","fEleTrigger", NBinsElectron, XminElectron, XmaxElectron);
  fOutputList->Add(fElecTrigger);
    
  fInclElecPhi = new TH2F("fInclElePhi", "fInclElePhi",NBinsElectron, XminElectron, XmaxElectron,NBinsPhi,0,TMath::TwoPi());
  fInclElecPhi->Sumw2();
  fOutputList->Add(fInclElecPhi);
  
  fInclElecEta = new TH2F("fInclEleEta", "fInclEleEta",NBinsElectron, XminElectron, XmaxElectron,NBinsEta,XminEta,XmaxEta);
  fInclElecEta->Sumw2();
  fOutputList->Add(fInclElecEta);
  
  fULSElecPhi= new TH2F("fULSElePhi", "fULSElePhi", NBinsElectron, XminElectron, XmaxElectron ,NBinsPhi,0,TMath::TwoPi());
  fULSElecPhi->Sumw2();
  fOutputList->Add(fULSElecPhi);
  
  fLSElecPhi= new TH2F("fLSElePhi", "fLSElePhi",NBinsElectron, XminElectron, XmaxElectron,NBinsPhi,0,TMath::TwoPi());
  fLSElecPhi->Sumw2();
  fOutputList->Add(fLSElecPhi);
  
  fElecDphi = new TH2F("fEleDphi", "fEleDphi",NBinsElectron, XminElectron, XmaxElectron,NBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fElecDphi->Sumw2();
  fOutputList->Add(fElecDphi);
  
  fULSElecDphi = new TH2F("fULSEleDphi", "fULSEleDphi",NBinsElectron, XminElectron, XmaxElectron,NBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fULSElecDphi->Sumw2();
  fOutputList->Add(fULSElecDphi);
  
  fLSElecDphi = new TH2F("fLSEleDphi", "fLSEleDphi",NBinsElectron, XminElectron, XmaxElectron,NBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fLSElecDphi->Sumw2();
  fOutputList->Add(fLSElecDphi);
  
  fULSElecDphiDiffMethod = new TH2F("fULSEleDphiDiffMethod", "fULSEleDphiDiffMethod",NBinsElectron, XminElectron, XmaxElectron,NBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fULSElecDphi->Sumw2();
  fOutputList->Add(fULSElecDphiDiffMethod);
  
  fLSElecDphiDiffMethod = new TH2F("fLSEleDphiDiffMethod", "fLSEleDphiDiffMethod",NBinsElectron, XminElectron, XmaxElectron,NBinsDPhi,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fLSElecDphi->Sumw2();
  fOutputList->Add(fLSElecDphiDiffMethod);
  
  Int_t     bin[4] = {NBinsHadron   ,NBinsElectron, NBinsDPhi, NBinsDEta}; //ptH, ptE, Dphi, Deta
  Double_t  xmin[4] = {XminHadron   ,XminElectron ,-TMath::Pi()/2    ,DEtaMin};
  Double_t  xmax[4] = {XmaxHadron   ,XmaxElectron ,(3*TMath::Pi())/2 ,DEtaMax};


  fNoPartnerNoT = new TH1F("fNoPartnerNoT", "fNoParnterNoT", bin[1], xmin[1], xmax[1]);
  fNoPartnerNoT->Sumw2();
  fOutputList->Add(fNoPartnerNoT);
  
  fTPartnerNoT = new TH1F("fTPartnerNoT", "fTPartnerNoT", bin[1], xmin[1], xmax[1]);
  fTPartnerNoT->Sumw2();
  fOutputList->Add(fTPartnerNoT);

  fElecHadTrigger = new TH2F("fElecHadTrigger", "fElecHadTrigger", bin[1], xmin[1], xmax[1], fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fElecHadTrigger->Sumw2();
  fOutputList->Add(fElecHadTrigger);						      

  fElecHadTriggerLS = new TH2F("fElecHadTriggerLS", "fElecHadTriggerLS", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fElecHadTriggerLS->Sumw2();
  fOutputList->Add(fElecHadTriggerLS);						      

  fElecHadTriggerULS = new TH2F("fElecHadTriggerULS", "fElecHadTriggerULS", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fElecHadTriggerULS->Sumw2();
  fOutputList->Add(fElecHadTriggerULS);						      

  fHadContTrigger = new TH2F("fHadContTrigger", "fHadContTrigger", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fHadContTrigger->Sumw2();
  fOutputList->Add(fHadContTrigger);						      

  fNonElecHadTrigger = new TH2F("fNonElecHadTrigger", "fNonElecHadTrigger", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fNonElecHadTrigger->Sumw2();
  fOutputList->Add(fNonElecHadTrigger);						      

  fMCElecHaTruePartnerTrigger = new TH2F("fMCElecHaTruePartnerTrigger", "fMCElecHaTruePartnerTrigger", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fMCElecHaTruePartnerTrigger->Sumw2();
  fOutputList->Add(fMCElecHaTruePartnerTrigger);	

  fMCElecHaNoPartnerTrigger = new TH2F("fMCElecHaNoPartnerTrigger", "fMCElecHaNoPartnerTrigger", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fMCElecHaNoPartnerTrigger->Sumw2();
  fOutputList->Add(fMCElecHaNoPartnerTrigger);

  fHadElecTrigger = new TH2F("fHadElecTrigger", "fHadElecTrigger", bin[0], xmin[0], xmax[0],  fAssPtElec_Nbins, -0.5, fAssPtElec_Nbins-0.5);
  fHadElecTrigger->Sumw2();
  fOutputList->Add(fHadElecTrigger);						      

  fHadNonElecTrigger = new TH2F("fHadNonElecTrigger", "fHadNonElecTrigger", bin[0], xmin[0], xmax[0],  fAssPtElec_Nbins, -0.5, fAssPtElec_Nbins-0.5);
  fHadNonElecTrigger->Sumw2();
  fOutputList->Add(fHadNonElecTrigger);	

  fInclElecHa = new THnSparseF("fEleHaIncl", "Sparse for Ele-Had : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fInclElecHa->Sumw2();
  fOutputList->Add(fInclElecHa);

  fLSElecHa = new THnSparseF("fEleHaLS", "Sparse for LSEle-Had : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fLSElecHa->Sumw2();
  fOutputList->Add(fLSElecHa);

  fULSElecHa = new THnSparseF("fEleHaULS", "Sparse for ULSEle-Had : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fULSElecHa->Sumw2();
  fOutputList->Add(fULSElecHa);
  
  fMCElecHaHadron = new THnSparseF("fMCEleHaHadron", "Sparse for Ele-Had : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fMCElecHaHadron->Sumw2();
  fOutputList->Add(fMCElecHaHadron);

  fElecHaLSNoPartner = new THnSparseF("fEleHaLSNoPartner", "Sparse for LSEle-Had with no Partner: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecHaLSNoPartner->Sumw2();
  fOutputList->Add(fElecHaLSNoPartner);

  fElecHaULSNoPartner = new THnSparseF("fEleHaULSNoPartner", "Sparse for ULSEle-Had with no Partner: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecHaULSNoPartner->Sumw2();
  fOutputList->Add(fElecHaULSNoPartner);
 
  fMCElecHaTruePartner = new THnSparseF("fMCElecHaTruePartner", "Sparse for MC true photonics with no Partner: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fMCElecHaTruePartner->Sumw2();
  fOutputList->Add(fMCElecHaTruePartner);

  fMCElecHaNoPartner = new THnSparseF("fMCElecHaNoPartner", "Sparse for MC true photonics with no Partner: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fMCElecHaNoPartner->Sumw2();
  fOutputList->Add(fMCElecHaNoPartner);	

  fElecHaHa = new THnSparseF("fEleHaHa", "Sparse for Hadron Contamination: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecHaHa->Sumw2();
  fOutputList->Add(fElecHaHa);

  fElecHaMixedEvent = new THnSparseF("fEleHaMixedEv", "Sparse for Ele-Had MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecHaMixedEvent->Sumw2();
  fOutputList->Add(fElecHaMixedEvent);
  
  fLSElecHaMixedEvent = new THnSparseF("fEleHaLSMixedEv", "Sparse for LSEle-Had MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fLSElecHaMixedEvent->Sumw2();
  fOutputList->Add(fLSElecHaMixedEvent);
  
  fULSElecHaMixedEvent = new THnSparseF("fEleHaULSMixedEv", "Sparse for ULSEle-Had MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fULSElecHaMixedEvent->Sumw2();
  fOutputList->Add(fULSElecHaMixedEvent);

  fTagHaMixedEvent = new THnSparseF("fTagHaMixedEv", "Sparse for ULSEle-Had MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fTagHaMixedEvent->Sumw2();
  fOutputList->Add(fTagHaMixedEvent);

  fNonTagHaMixedEvent = new THnSparseF("fNonTagMixedEv", "Sparse for ULSEle-Had MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fNonTagHaMixedEvent->Sumw2();
  fOutputList->Add(fNonTagHaMixedEvent);

  // LP HFE

					      

 				      




  fElecLPTrigger = new TH2F("fElecLPTrigger", "fElecLPTrigger", bin[1], xmin[1], xmax[1], fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fElecLPTrigger->Sumw2();
  fOutputList->Add(fElecLPTrigger);	
						      

  fElecLPTriggerLS = new TH2F("fElecLPTriggerLS", "fElecLPTriggerLS", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fElecLPTriggerLS->Sumw2();
  fOutputList->Add(fElecLPTriggerLS);						      

  fElecLPTriggerULS = new TH2F("fElecLPTriggerULS", "fElecLPTriggerULS", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fElecLPTriggerULS->Sumw2();
  fOutputList->Add(fElecLPTriggerULS);	

  fHadContLPTrigger = new TH2F("fHadContLPTrigger", "fHadContLPTrigger", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fHadContLPTrigger->Sumw2();
  fOutputList->Add(fHadContLPTrigger);		

  fNonElecLPTrigger = new TH2F("fNonElecLPTrigger", "fNonElecLPTrigger", bin[1], xmin[1], xmax[1], fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fNonElecLPTrigger->Sumw2();
  fOutputList->Add(fNonElecLPTrigger);					      
					      

  fMCElecLPTruePartnerTrigger = new TH2F("fMCElecLPTruePartnerTrigger", "fMCElecLPTruePartnerTrigger", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fMCElecLPTruePartnerTrigger->Sumw2();
  fOutputList->Add(fMCElecLPTruePartnerTrigger);	

  fMCElecLPNoPartnerTrigger = new TH2F("fMCElecLPNoPartnerTrigger", "fMCElecLPNoPartnerTrigger", bin[1], xmin[1], xmax[1],  fAssPtHad_Nbins, -0.5, fAssPtHad_Nbins-0.5);
  fMCElecLPNoPartnerTrigger->Sumw2();
  fOutputList->Add(fMCElecLPNoPartnerTrigger);

  					      

  fLPElecTrigger = new TH2F("fLPElecTrigger", "fLPElecTrigger", bin[0], xmin[0], xmax[0],  fAssPtElec_Nbins, -0.5, fAssPtElec_Nbins-0.5);
  fLPElecTrigger->Sumw2();
  fOutputList->Add(fLPElecTrigger);						      

  fLPNonElecTrigger = new TH2F("fLPNonElecTrigger", "fLPNonElecTrigger", bin[0], xmin[0], xmax[0], fAssPtElec_Nbins, -0.5, fAssPtElec_Nbins-0.5);
  fLPNonElecTrigger->Sumw2();
  fOutputList->Add(fLPNonElecTrigger);						      




  fInclElecLP = new THnSparseF("fEleLPIncl", "Sparse for Ele-LP : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fInclElecLP->Sumw2();
  fOutputList->Add(fInclElecLP);
    
  fLSElecLP = new THnSparseF("fEleLPLS", "Sparse for LSEle-LP : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fLSElecLP->Sumw2();
  fOutputList->Add(fLSElecLP);
    
  fULSElecLP = new THnSparseF("fEleLPULS", "Sparse for ULSEle-LP : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fULSElecLP->Sumw2();
  fOutputList->Add(fULSElecLP);

  fMCElecLPHadron = new THnSparseF("fMCEleLPHadron", "Sparse for Ele-Had : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fMCElecLPHadron->Sumw2();
  fOutputList->Add(fMCElecLPHadron);
    
  fElecLPLSNoPartner = new THnSparseF("fEleLPLSNoPartner", "Sparse for LSEle-Had with no Partner: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecLPLSNoPartner->Sumw2();
  fOutputList->Add(fElecLPLSNoPartner);
    
  fElecLPULSNoPartner = new THnSparseF("fEleLPULSNoPartner", "Sparse for ULSEle-Had with no Partner: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecLPULSNoPartner->Sumw2();
  fOutputList->Add(fElecLPULSNoPartner);
    
  fElecLPHa = new THnSparseF("fEleLPHa", "Sparse for Hadron Contamination: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecLPHa->Sumw2();
  fOutputList->Add(fElecLPHa);

  fMCElecLPTruePartner = new THnSparseF("fMCElecLPTruePartner", "Sparse for MC true photonics with no Partner: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fMCElecLPTruePartner->Sumw2();
  fOutputList->Add(fMCElecLPTruePartner);

  fMCElecLPNoPartner = new THnSparseF("fMCElecLPNoPartner", "Sparse for MC true photonics with no Partner: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fMCElecLPNoPartner->Sumw2();
  fOutputList->Add(fMCElecLPNoPartner);	
    
  fElecLPMixedEvent = new THnSparseF("fEleLPMixedEv", "Sparse for Ele-LP MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecLPMixedEvent->Sumw2();
  fOutputList->Add(fElecLPMixedEvent);
    
  fLSElecLPMixedEvent = new THnSparseF("fEleLPLSMixedEv", "Sparse for LSEle-LP MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fLSElecLPMixedEvent->Sumw2();
  fOutputList->Add(fLSElecLPMixedEvent);
    
  fULSElecLPMixedEvent = new THnSparseF("fEleLPULSMixedEv", "Sparse for ULSEle-LP MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fULSElecLPMixedEvent->Sumw2();
  fOutputList->Add(fULSElecLPMixedEvent);
      
  if (fIsMC) {
    fCheckMCVertex = new TH2F("fCheckMCVertex", "TruthVsRec Vertex", 200, -10, 10, 200, -10, 10);
    fCheckMCVertex->Sumw2();
    fOutputList->Add(fCheckMCVertex);
  }
  // NBinsEta/2.5 makes 0.5 steps
  Int_t  EffHBins[4]={NBinsHadron, 36, NBinsPhi/2, 50};
  Double_t EffHXmin[4]={XminHadron, -0.9, 0, -10};
  Double_t EffHXmax[4]={XmaxHadron, 0.9, TMath::TwoPi(), 10};

  fRecHadPtEtaPhiVtx=new THnSparseF("fRecHadPtEtaPhiVtx", "Rec hadrons w. rec. pt, eta, phi, vtxz", 4, EffHBins, EffHXmin, EffHXmax);
  fRecHadPtEtaPhiVtx->Sumw2();
  fOutputList->Add(fRecHadPtEtaPhiVtx);
  
  if (fIsMC) {
    fMCHadPtEtaPhiVtx=new THnSparseF("fMCHadPtEtaPhiVtx", "MC truth gen. hadrons pt, eta, phi, vtxz", 4, EffHBins, EffHXmin, EffHXmax);
    fMCHadPtEtaPhiVtx->Sumw2();
    fOutputList->Add(fMCHadPtEtaPhiVtx);
       
    fRecHadMCPtEtaPhiVtx=new THnSparseF("fRecHadMCPtEtaPhiVtx", "MC rec hadrons w. MC pt, eta, phi, vtxz", 4, EffHBins, EffHXmin, EffHXmax);
    fRecHadMCPtEtaPhiVtx->Sumw2();
    fOutputList->Add(fRecHadMCPtEtaPhiVtx);    

    fCheckMCPtvsRecPtHad = new TH2F("fCheckMCPtvsRecPtHad", "MCvsRec Pt", NBinsHadron, XminHadron, XmaxHadron, NBinsHadron, XminHadron, XmaxHadron); //200, 0, 20, 200, 0, 20);
    fCheckMCPtvsRecPtHad->Sumw2();
    //    BinLogX(fCheckMCPtvsRecPtHad->GetXaxis());
    //    BinLogX(fCheckMCPtvsRecPtHad->GetYaxis());
    fOutputList->Add(fCheckMCPtvsRecPtHad);

    fCheckMCEtavsRecEtaHad = new TH2F("fCheckMCEtavsRecEtaHad", "MC vs. Rec; rec. #eta; MC #eta", 36, -0.9, 0.9, 36, -0.9, 0.9);
    fCheckMCEtavsRecEtaHad->Sumw2();
    fOutputList->Add(fCheckMCEtavsRecEtaHad);

    fCheckMCPhivsRecPhiHad = new TH2F("fCheckMCPhivsRecPhiHad", "MC vs. Rec; rec. #eta; MC #eta", NBinsPhi/2, 0, TMath::TwoPi(), NBinsPhi/2, 0, TMath::TwoPi());
    fCheckMCPhivsRecPhiHad->Sumw2();
    fOutputList->Add(fCheckMCPhivsRecPhiHad);

  }   
    

  Int_t  EffEBins[4]={NBinsElectron, 36, NBinsPhi/2, 50};
  Double_t EffEXmin[4]={XminElectron, -0.9, 0, -10};
  Double_t EffEXmax[4]={XmaxElectron, 0.9, TMath::TwoPi(), 10};

  fRecElecPtEtaPhiVtx=new THnSparseF("fRecElePtEtaPhi", "Rec electrons w. rec.  pt, eta, phi", 4, EffEBins, EffEXmin, EffEXmax);
  fRecElecPtEtaPhiVtx->Sumw2();
  fOutputList->Add(fRecElecPtEtaPhiVtx);
    
  if (fIsMC) {
  
    fMCElecPtEtaPhiVtx=new THnSparseF("fMCElePtEtaPhi", "MC truth gen. electrons  pt, eta, phi", 4, EffEBins, EffEXmin, EffEXmax);
    fMCElecPtEtaPhiVtx->Sumw2();
    fOutputList->Add(fMCElecPtEtaPhiVtx);

    fRecElecMCPtEtaPhiVtx=new THnSparseF("fRecEleMCPtEtaPhi", "MC rec electrons w. MC  pt, eta, phi", 4, EffEBins, EffEXmin, EffEXmax);
    fRecElecMCPtEtaPhiVtx->Sumw2();
    fOutputList->Add(fRecElecMCPtEtaPhiVtx);
    
    fCheckMCPtvsRecPtEle = new TH2F("fCheckMCPtvsRecPtEle", "MCvsRec Pt", NBinsElectron, XminElectron, XmaxElectron, NBinsElectron, XminElectron, XmaxElectron); //200, 0, 20, 200, 0, 20);
    fCheckMCPtvsRecPtEle->Sumw2();
    //    BinLogX(fCheckMCPtvsRecPtEle->GetXaxis());
    //    BinLogX(fCheckMCPtvsRecPtEle->GetYaxis());
    fOutputList->Add(fCheckMCPtvsRecPtEle);


    fMCElecPDG=new TH1F("fMCElePDG", "MC truth mother of heavy electrons", 10000, -0.5, 9999.5);
    fMCElecPDG->Sumw2();
    fOutputList->Add(fMCElecPDG);
    
    fMCElecPtEtaPhiStrictVtx=new THnSparseF("fMCElePtEtaPhiStrict", "MC truth electrons  pt, eta, phi", 4, EffEBins, EffEXmin, EffEXmax);
    fMCElecPtEtaPhiStrictVtx->Sumw2();
    fOutputList->Add(fMCElecPtEtaPhiStrictVtx);
  }

  Int_t Pi0EtaBins[5] ={100, 80, 80, 10001,5};
  Double_t Pi0EtaXmin[5]={0.1, -2,-2, -1.5, 0};
  Double_t Pi0EtaXmax[5]={20, 2,2,  9999.5,5};

  if (fIsMC) {
    fMCPi0Prod = new THnSparseF("fMCPi0Prod", "fMCPi0Prod: pt eta pdgmother", 5, Pi0EtaBins, Pi0EtaXmin, Pi0EtaXmax);
    BinLogX(fMCPi0Prod->GetAxis(0));
    fMCPi0Prod->Sumw2();
    fOutputList->Add(fMCPi0Prod);

    fMCEtaProd = new THnSparseF("fMCEtaProd", "fMCEtaProd: pt eta pdgmother", 5, Pi0EtaBins, Pi0EtaXmin, Pi0EtaXmax);
    BinLogX(fMCEtaProd->GetAxis(0));
    fMCEtaProd->Sumw2();
    fOutputList->Add(fMCEtaProd);

    fMCPiPlusProd = new THnSparseF("fMCPiPlusProd", "fMCPiPlusProd: pt eta pdgmother", 5, Pi0EtaBins, Pi0EtaXmin, Pi0EtaXmax);
    BinLogX(fMCPiPlusProd->GetAxis(0));
    fMCPiPlusProd->Sumw2();
    fOutputList->Add(fMCPiPlusProd);

    fMCPiPlusProdV2 = new THnSparseF("fMCPiPlusProdV2", "fMCPiPlusProdV2: pt eta pdgmother", 5, Pi0EtaBins, Pi0EtaXmin, Pi0EtaXmax);
    BinLogX(fMCPiPlusProdV2->GetAxis(0));
    fMCPiPlusProdV2->Sumw2();
    fOutputList->Add(fMCPiPlusProdV2);
    
    
    Int_t LPBins[3]={200, 10000, 10001};
    Double_t LPBinsXmin[3]={0, -0.5, -1.5};
    Double_t LPBinsXmax[3]={200, 9999.5, 9999.5};
    fMCLeadingParticle = new THnSparseF("fMCLeadingParticle", "fMCLeadingParticle: pt, pdg, pdgmother", 3, LPBins, LPBinsXmin, LPBinsXmax);
    fMCLeadingParticle->Sumw2();
    fOutputList->Add(fMCLeadingParticle);

  }

  if (fTRDQA) {
    fhArmenteros  = new TH2F("fhArmenteros","Armenteros plot",200,-1.,1.,200,0.,0.4);
    fhArmenteros->Sumw2();
    fOutputList->Add(fhArmenteros);

    fEventsPerRun = new TH1F("fEventsPerRun", "EventsPerRun",500, 0.5, 500.5);
    fEventsPerRun->Sumw2();
    fOutputList->Add(fEventsPerRun);
    
    fTRDnTrackRun = new TH3F("fTRDnTrackRun", "TrackletCharge, TrackletPID, Run", 6, 0.5, 6.5, 6, 0.5, 6.5, 201, -0.5, 200.5);
    fTRDnTrackRun->Sumw2();
    fOutputList->Add(fTRDnTrackRun);
    
    Int_t TRDBins[6]={11, 30, 54, 7, 3, 201};
    Double_t TRDBMin[6]={0.5, XminEta, 0, -0.5, -1.5,-0.5};
    Double_t TRDBMax[6]={XmaxElectron, XmaxEta, TMath::TwoPi(), 6.5, 1.5, 200.5};

    fTRDEtaPhi = new THnSparseF("fTRDEtaPhi", "fTRDEtaPhi: Pt, Eta, Phi, Layer, Charge, run", 6, TRDBins, TRDBMin, TRDBMax);
    fTRDEtaPhi->Sumw2();
    fOutputList->Add(fTRDEtaPhi);

    Int_t TRDBins3[7]=   {11             ,   9  ,      3,    7,     7,    3, 500};
    Double_t TRD3BMin[7]={0.5         , XminEta,   0.5,   -3,  -0.5, -1.5, 0.5};
    Double_t TRD3BMax[7]={XmaxElectron, XmaxEta,   3.5,    4,   6.5,  1.5, 500.5};

    fTRDV0NTracklets = new THnSparseF("fTRDV0NTracklets", "fTRDnTracklets: Pt, Eta, PID, TPC, NTracklets, Charge, run", 7, TRDBins3, TRD3BMin, TRD3BMax);
    fTRDV0NTracklets->Sumw2();
    fOutputList->Add(fTRDV0NTracklets);
  
    fTRDNTracklets = new THnSparseF("fTRDNTracklets", "fTRDnTracklets: Pt, Eta, PID, TPC, NTracklets, Charge, run",  7, TRDBins3, TRD3BMin, TRD3BMax);
    fTRDNTracklets->Sumw2();
    fOutputList->Add(fTRDNTracklets);

    Int_t TRDBins2[5]={NBinsElectron/2,    7,   3,   7, 3};
    Double_t TRDB2Min[5]={XminElectron, -0.5, 0.5,  -3, -1.5};
    Double_t TRDB2Max[5]={XmaxElectron,  6.5, 3.5,   4,  1.5};   

    fTRDV0Spectra = new THnSparseF("fTRDV0Spectra", "fTRDV0Spectra: Pt, NTracklets, PID, TPC, Charge", 5, TRDBins2, TRDB2Min, TRDB2Max);
    fTRDV0Spectra->Sumw2();
    fOutputList->Add(fTRDV0Spectra);

    fTRDSpectra = new THnSparseF("fTRDSpectra", "fTRDSpectra: Pt, NTracklets, PID, TPC, Charge", 5, TRDBins2, TRDB2Min, TRDB2Max); 
    fTRDSpectra->Sumw2();
    fOutputList->Add(fTRDSpectra);

    fTRDMCSpectra= new THnSparseF("fTRDMCSpectra","fTRDMCSpectra: Pt, NTracklets, PID, TPC, Charge", 5, TRDBins2, TRDB2Min, TRDB2Max);  
    fTRDMCSpectra->Sumw2();
    fOutputList->Add(fTRDMCSpectra);
  }

  Int_t    poolSize = 1000;
  Int_t    trackDepth = 2000; 
  Int_t    nMultBins = 3;
  Double_t multBins[]={0,20,50,100};
  Int_t    nZVtxBins=6;
  Double_t zVtxBins[]={-10,-4.41,-1.41,0.59,3.59,10};
  Int_t    nMaxPtBins=7;
  Double_t maxPtBins[]={0, 0.5, 1., 2., 5., 10, 15, 999};
  //Int_t    nDummyBins=1;
  // Double_t dummyBins[]=


  fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, nMultBins, (Double_t*) multBins, nZVtxBins, (Double_t*) zVtxBins, nMaxPtBins,(Double_t*)maxPtBins);
  fPoolMgr->Validate();
    
    
  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHaHFECorrel::Terminate(Option_t *)
{
  // Info("Terminate");
  AliAnalysisTaskSE::Terminate();
}

//_________________________________________
Double_t AliAnalysisTaskHaHFECorrel::GetDeltaPhi(Double_t phiA,Double_t phiB) const
{
  //Get phi_tag-psi_assoc
  Double_t dPhi = phiA - phiB;
  if (dPhi < -0.5*TMath::Pi()) dPhi = dPhi + TMath::TwoPi();
  if (dPhi >  1.5*TMath::Pi()) dPhi = dPhi - TMath::TwoPi();
  return dPhi;
}


//_________________________________________
Double_t AliAnalysisTaskHaHFECorrel::GetDeltaEta(Double_t etaA,Double_t etaB) const
{
  Double_t dEta=etaB-etaA;  
  return dEta;
}

//_________________________________________
AliVTrack*  AliAnalysisTaskHaHFECorrel::FindLPAndHFE( TObjArray* RedTracks, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[])
{
  AliVTrack* LPtrack=0;
  fLParticle=kFALSE;
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();
  
  // Leading Particle
  Double_t ptH=-9;

  // Check Hadron Correlation
  Int_t NHadrons=0, NElectrons=0, NNonElectrons=0, NHadCont=0;
  Int_t NPhotElectronsUntagged=0, NPhotElectronsTagged=0;

  // Loop over all tracks to find LP and HFE
  for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
          
    AliVParticle* VHtrack = 0x0;
    if(!fUseTender) VHtrack  = fVevent->GetTrack(jTracks);
    if(fUseTender) VHtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jTracks));
    if (!VHtrack) {
      printf("ERROR: Could not receive associated track %d\n", jTracks);
      continue;
    }
    AliVTrack *Vtrack = dynamic_cast<AliVTrack*>(VHtrack);
    if(!Vtrack) continue;
    AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(VHtrack);
    if (fIsAOD && AODtrack==0) continue;

    Double_t p=-9.,pt=-9.,eta =-9.,phi=-9.;
    pt = Vtrack->Pt();
    p = Vtrack->P();
    phi = Vtrack->Phi();
    eta = Vtrack->Eta();
    
    fTrkpt->Fill(pt,0);

    // track cuts
    Bool_t passHadTrackCut=kFALSE;
    passHadTrackCut = ChargedHadronTrackCuts(pVtx, Vtrack, nMother, listMother);
   
    Bool_t passHadPIDCut=kFALSE;
    passHadPIDCut = ChargedHadronPIDCuts(Vtrack); // currently empty
    
    if (passHadTrackCut && passHadPIDCut) NHadrons++;

    // find hadron with the largest pT -> leading particle
    if (passHadTrackCut && passHadPIDCut) {

      // for LP
      if(Vtrack->Pt()>ptH) {
	ptH = Vtrack->Pt();
	LPtrack=Vtrack;
	fLParticle=kTRUE;
      }

      fEtaVtxZ->Fill(eta, pVtx->GetZ());

      Double_t fillSparse[4];
      fillSparse[0]=pt;
      fillSparse[1]=eta;
      fillSparse[2]=phi;
      fillSparse[3]=pVtx->GetZ();
      fRecHadPtEtaPhiVtx->Fill(fillSparse);

      if (fIsAOD && fIsMC) {
	Int_t MClabel=AODtrack->GetLabel();
	AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));  
	Double_t mcVtx[3];
	fMCheader->GetVertex(mcVtx);
	fillSparse[0]=MCParticle->Pt();
	fillSparse[1]=MCParticle->Eta();
	fillSparse[2]=MCParticle->Phi();
	fillSparse[3]=mcVtx[2];
	fRecHadMCPtEtaPhiVtx->Fill(fillSparse);


	// controll plot for MC vs Rec pt;
	fCheckMCPtvsRecPtHad->Fill(MCParticle->Pt(), pt);
	fCheckMCEtavsRecEtaHad->Fill(MCParticle->Eta(), eta);
	fCheckMCPhivsRecPhiHad->Fill(MCParticle->Phi(), phi);
      }
    }  
 
    Bool_t passHFETrackCut=kFALSE;
    passHFETrackCut= InclElecTrackCuts(pVtx,Vtrack,nMother,listMother);
    if (passHFETrackCut)  fTrkpt->Fill(pt,1); // after track cuts
 
    Bool_t passHFEPIDCut=kFALSE;
    if (passHFETrackCut) passHFEPIDCut=  InclElecPIDCuts(Vtrack, kTRUE);
    if (passHFETrackCut && passHFEPIDCut)  fTrkpt->Fill(pt,2);
    
    Int_t lsPartner=0, ulsPartner=0;
    Int_t lsPartnerID[20], ulsPartnerID[20];
    Bool_t trueULSPartner = kFALSE;
    Bool_t isPhotonic = kFALSE;
    Bool_t isHadron = kFALSE;
    if (passHFETrackCut && passHFEPIDCut) { // if HFE is found, look for ls and uls partner
      NElectrons++;
      FindPhotonicPartner(jTracks, Vtrack, pVtx, nMother, listMother, lsPartner, ulsPartner, lsPartnerID, ulsPartnerID, trueULSPartner, isPhotonic);
      if (fIsMC) {
	if (isPhotonic) {
	  if (trueULSPartner) {
	    NPhotElectronsTagged++;
	    fTagEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt());
	  }
	  else  {
	    NPhotElectronsUntagged++;
	    fNonTagEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt());
	  }
	}
	EvaluateTaggingEfficiency(Vtrack, lsPartner, ulsPartner, trueULSPartner);
      }
      else {
	fTagEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), ulsPartner);
	fTagEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), -lsPartner);
      }
      fInclElecPt->Fill(pt);
      fInclElecP->Fill(p);
      fInclElecPhi->Fill(pt,phi); // phi of electron candidates
      fInclElecEta->Fill(pt,eta); 

      if (fIsMC & fIsAOD) {
	Int_t MClabel=AODtrack->GetLabel();
	AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));     
	Int_t PDGCode = abs(MCParticle->GetPdgCode());
	Double_t mcVtx[3];
	fMCheader->GetVertex(mcVtx);

	if (PDGCode ==11) {
	  AliAODMCParticle* mcPartMother=(AliAODMCParticle*)fMC->GetTrack(MCParticle->GetMother());
	  Int_t mcMotherPDG = abs(mcPartMother->GetPdgCode()); 
	  Int_t MIsHeavy  = Int_t (mcMotherPDG / TMath::Power(10, Int_t(TMath::Log10(mcMotherPDG))));
	  if (MIsHeavy>3) {

	    // Fill Reconstructed Histogram
	    Double_t fillSparse[4];
	    fillSparse[0]=pt;
	    fillSparse[1]=eta;
	    fillSparse[2]=phi;
	    fillSparse[3]=pVtx->GetZ();
	    fRecElecPtEtaPhiVtx->Fill(fillSparse);

	    fillSparse[0]=MCParticle->Pt();
	    fillSparse[1]=MCParticle->Eta();
	    fillSparse[2]=MCParticle->Phi();
	    fillSparse[3]=mcVtx[2];
	    fRecElecMCPtEtaPhiVtx->Fill(fillSparse);

	    fCheckMCPtvsRecPtEle->Fill(MCParticle->Pt(), pt);
	  }
	}

	if (PDGCode != 11) {
	  NHadCont++;
	  isHadron=kTRUE;
	  fHadContEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt());
	}
      }
        
      // store only essential informations of these electrons for later correlations and mixed event
      CloneAndReduceTrackList(RedTracks, Vtrack, lsPartner, ulsPartner, lsPartnerID, ulsPartnerID, trueULSPartner, isPhotonic, isHadron);
    }

    Bool_t passNonElecPIDCut=kFALSE;
    passNonElecPIDCut=AssoHadronPIDCuts(Vtrack);
    if (passHFETrackCut && passNonElecPIDCut) {
      NNonElectrons++;
      fHadContTPCEtaPhiPt->Fill(Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt());
    }
  }

  // CheckNoOfElectrons/Event
  Double_t fillSparse[4];
  fillSparse[0]=NHadrons;
  fillSparse[1]=NElectrons;
  fillSparse[2]=NNonElectrons;
  if (fIsMC) fillSparse[3] = NHadCont;
  else fillSparse[3]=0;
  fCheckNHadronScaling->Fill(fillSparse);

  if (fIsMC) { 
    fillSparse[2]= NPhotElectronsTagged;
    fillSparse[3]= NPhotElectronsUntagged;
    fCheckNPhotHadScaling->Fill(fillSparse);
  }

  // if no leading particle is found, mark this in flag fLParticle;
  if (!LPtrack) {
    fLParticle=kFALSE;
  }
  else if (fIsMC && fIsAOD && RedTracks->GetEntriesFast()>0) {
    AliAODTrack *LPAODtrack = dynamic_cast<AliAODTrack*>(LPtrack);   
    if (!LPAODtrack) {
      fLParticle=kFALSE;
      return 0;
    }
    Int_t MClabel=LPAODtrack->GetLabel();
    AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));     
    Int_t PDGCode = abs(MCParticle->GetPdgCode());
    Int_t PDGCodeMother=0;
    if (MCParticle->GetMother()>=0) {
      AliAODMCParticle* MCParticleMother = (AliAODMCParticle*) fMC->GetTrack(MCParticle->GetMother());
      PDGCodeMother = abs(MCParticleMother->GetPdgCode());
    }
    else {
      PDGCodeMother=-1;
    }
    Double_t fillSparse[3];
    fillSparse[0]=LPAODtrack->Pt();
    fillSparse[1]=PDGCode;
    fillSparse[2]=PDGCodeMother;
    fMCLeadingParticle->Fill(fillSparse);
  }

  return LPtrack;
}

//_________________________________________
void AliAnalysisTaskHaHFECorrel::FindPhotonicPartner(Int_t iTracks, AliVTrack* Vtrack, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[], Int_t &lsPartner, Int_t &ulsPartner, Int_t *lsPartnerID, Int_t *ulsPartnerID, Bool_t &trueULSPartner, Bool_t &isPhotonic) {
  

  AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);
  if(!AODtrack) return;

  lsPartner=0;
  ulsPartner=0;
  trueULSPartner=kFALSE; // only for MC use
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();
    
  // Loop over all tracks to find photonic partner
  for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
        
    if(jTracks==iTracks) continue;
    
    AliVParticle* Vassotrack = 0x0;
    if(!fUseTender) Vassotrack  = fVevent->GetTrack(jTracks);
    if(fUseTender)  Vassotrack  = dynamic_cast<AliVTrack*>(fTracks_tender->At(jTracks));
    if (!Vassotrack) {
      printf("ERROR: Could not receive associated track %d\n", jTracks);
      continue;
    }
    AliVTrack *VtrackAsso = dynamic_cast<AliVTrack*>(Vassotrack);
    if(!VtrackAsso) continue;
    AliAODTrack *AODtrackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);   
    if(!AODtrackAsso) continue;
  
    // tagged particle
    Double_t pt=-9.,phi=-9.;
    Int_t charge = 0;      
    pt = Vtrack->Pt();
    phi = Vtrack->Phi();
    charge = Vtrack->Charge();
        
    // associated particle variables
    Double_t phiAsso=-9.;
    Int_t chargeAsso = 0;
    phiAsso = VtrackAsso->Phi();
    chargeAsso = VtrackAsso->Charge();

    double dphi = -9;

    // Cuts for associated electrons for HFE-HFE correlations
    Bool_t passAssoTrackCutIncl = kFALSE;
    passAssoTrackCutIncl = InclElecTrackCuts(pVtx,VtrackAsso,nMother,listMother);

    Bool_t passAssoPIDCutIncl = kFALSE;
    passAssoPIDCutIncl = InclElecPIDCuts(VtrackAsso, kFALSE);

    if (passAssoTrackCutIncl && passAssoPIDCutIncl) {
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      fElecDphi->Fill(pt,dphi);        // HFEHFE
    }

    // looser Track cuts for associated photonic eletron
    Bool_t passAssoTrackCutPhot = kFALSE;
    passAssoTrackCutPhot = PhotElecTrackCuts(pVtx,VtrackAsso,nMother,listMother);
    if(!passAssoTrackCutPhot) continue;
   
    Bool_t passAssoPIDCutPhot = kFALSE;
    passAssoPIDCutPhot = PhotElecPIDCuts(VtrackAsso); 
    if (!passAssoPIDCutPhot) continue;         
   
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    Double_t openingAngle = -999., mass=999., width = -999;
        
    Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
    if(charge>0) fPDGe1 = -11;
    if(chargeAsso>0) fPDGe2 = -11;
        
    if(charge == chargeAsso) fFlagLS = kTRUE;
    if(charge != chargeAsso) fFlagULS = kTRUE;
    
    // Two different Alogrithms in use (check difference) A. KalmanFilter (old), B. DCA
    if (fUseKFforPhotonicPartner) {        
      AliKFParticle::SetField(fVevent->GetMagneticField());
      AliKFParticle ge1(*Vtrack, fPDGe1);
      AliKFParticle ge2(*VtrackAsso, fPDGe2);
      AliKFParticle recg(ge1, ge2);
        
      if(recg.GetNDF()<1) continue;
      Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
      if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
      
      openingAngle = ge1.GetAngle(ge2);
      //if(openingAngle > fOpeningAngleCut) continue;
      
      Int_t MassCorrect=-9;
      MassCorrect = recg.GetMass(mass,width);     
    }
    else {
      //Variables
      Double_t p1[3], p2[3];
      Double_t xt1, xt2; //radial position track 1 and 2  at the DCA point
      Double_t dca12; //DCA 1-2
      Bool_t hasdcaT1, hasdcaT2;
      Double_t bfield = fVevent->GetMagneticField();
            
      AliExternalTrackParam extTrackParam1;
      extTrackParam1.CopyFromVTrack(Vtrack);
      AliExternalTrackParam extTrackParam2;
      extTrackParam2.CopyFromVTrack(VtrackAsso);

      //DCA track1-track2
      dca12 = extTrackParam2.GetDCA(&extTrackParam1,bfield,xt2,xt1);
                
      //Momento of the track extrapolated to DCA track-track
      //Track1
      hasdcaT1 = extTrackParam1.GetPxPyPzAt(xt1,bfield,p1);
      //Track2
      hasdcaT2 = extTrackParam2.GetPxPyPzAt(xt2,bfield,p2);
      
      if(!hasdcaT1 || !hasdcaT2) AliWarning("There could be a problem in the extrapolation");
      //track1-track2 Invariant Mass
      Double_t eMass = 0.000510998910; //Electron mass in GeV
      Double_t pP1 = sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2]); //Track 1 momentum
      Double_t pP2 = sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]); //Track 1 momentum

      TLorentzVector v1(p1[0],p1[1],p1[2],sqrt(eMass*eMass+pP1*pP1));
      TLorentzVector v2(p2[0],p2[1],p2[2],sqrt(eMass*eMass+pP2*pP2));
      mass = (v1+v2).M(); //Invariant Mass
      openingAngle = v1.Angle(v2.Vect()); //Opening Angle (Total Angle)
    }

    if(fFlagLS){
      fOpeningAngleLS->Fill(openingAngle,pt);
      fInvmassLS->Fill(mass,pt);
    }
    if(fFlagULS){
      fOpeningAngleULS->Fill(openingAngle,pt);
      fInvmassULS->Fill(mass,pt);
    }
        
    if(mass<fInvmassCut && fFlagULS) {
      if (fIsAOD && fIsMC && !trueULSPartner) {
	trueULSPartner = HaveSameMother(AODtrack->GetLabel(), AODtrackAsso->GetLabel());
      }
      fULSElecPhi->Fill(pt,phi);
      fULSElecPt->Fill(pt);
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      fULSElecDphi->Fill(pt,dphi);
      ulsPartnerID[ulsPartner]=VtrackAsso->GetID();
      ulsPartner++;
    }
        
    if(mass<fInvmassCut && fFlagLS){
      fLSElecPhi->Fill(pt,phi);
      fLSElecPt->Fill(pt);
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      fLSElecDphi->Fill(pt,dphi);
      lsPartnerID[lsPartner]=VtrackAsso->GetID();
      lsPartner++;
    }
  
  }//track loop
  if (fIsMC) { 
    isPhotonic=IsPhotonicElectron(AODtrack->GetLabel());
  } 
}

void AliAnalysisTaskHaHFECorrel::CorrelateElectron(TObjArray* RedTracksHFE) {

  // work in progress - only correlate electrons which fullfill inclusive cuts

  for (Int_t j=0; j<RedTracksHFE->GetEntriesFast(); j++) {

    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(j);

    // tagged particle
    Double_t pt=-9.,eta =-9.,phi=-9.;
    Int_t charge = 0;
    pt = RedTrack->Pt();
    phi = RedTrack->Phi();
    eta = RedTrack->Eta();
    charge = RedTrack->Charge();

    // fill trigger electron
    fElecTrigger->Fill(pt);

    for (Int_t k=0; k<RedTracksHFE->GetEntriesFast(); k++) {
      if (k==j) continue;
      AliBasicParticleHaHFE *trackAsso = (AliBasicParticleHaHFE*) RedTracksHFE->At(k);
      Double_t ptAsso=-9.,etaAsso =-9.,phiAsso=-9.;
      Int_t lsAsso=0, ulsAsso=0;
        
      ptAsso = trackAsso->Pt();
      phiAsso = trackAsso->Phi();
      etaAsso = trackAsso->Eta();
      lsAsso = trackAsso->LS();
      ulsAsso = trackAsso->ULS();
      
      // Fill Incl-Incl (already done above)
      Double_t dphi = -99;
      dphi = GetDeltaPhi(phi, phiAsso);
    
      for (Int_t l=0; l<lsAsso;l++)   fLSElecDphiDiffMethod->Fill(pt, dphi);
      for (Int_t l=0; l<ulsAsso; l++) fULSElecDphiDiffMethod->Fill(pt, dphi);
      
    }
  }
}


/// Maybe include LP case
void AliAnalysisTaskHaHFECorrel::CorrelateHadron(TObjArray* RedTracksHFE,  const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], Float_t mult, Float_t maxPt) {
    
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();

  Bool_t *HadronTrigger;
  Bool_t **ElectronIsTrigger = new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **HadContIsTrigger=new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **PhotElecWPartnerTrigger=new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **PhotElecWoPartnerTrigger=new Bool_t*[RedTracksHFE->GetEntriesFast()];
  for (Int_t j=0; j<RedTracksHFE->GetEntriesFast(); j++)
    {
      ElectronIsTrigger[j]=new Bool_t[fAssPtHad_Nbins];
      HadContIsTrigger[j]=new Bool_t[fAssPtHad_Nbins];
      PhotElecWPartnerTrigger[j] = new Bool_t[fAssPtHad_Nbins];
      PhotElecWoPartnerTrigger[j] = new Bool_t[fAssPtHad_Nbins];
      for (Int_t m=0; m<fAssPtHad_Nbins; m++) {
	ElectronIsTrigger[j][m]=kFALSE;
	HadContIsTrigger[j][m]=kFALSE;
	PhotElecWPartnerTrigger[j][m]=kFALSE;
	PhotElecWoPartnerTrigger[j][m]=kFALSE;
      }
    }
  Double_t *ElectronIsTriggerPt = new Double_t[RedTracksHFE->GetEntriesFast()];
  Bool_t **NonElectronIsTrigger = new Bool_t*[ntracks];
  for (Int_t j=0; j<ntracks; j++) { 
    NonElectronIsTrigger[j]=new Bool_t[fAssPtHad_Nbins]; 
    for (Int_t m=0; m<fAssPtHad_Nbins; m++) {
      NonElectronIsTrigger[j][m]=kFALSE;
    }
  }
  Double_t *NonElectronIsTriggerPt = new Double_t[ntracks];

  // Track loop for hadron correlations
  for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {     
    HadronTrigger = new Bool_t[fAssPtElec_Nbins];
    for (Int_t l=0; l<fAssPtElec_Nbins; l++) HadronTrigger[l]=kFALSE;

    AliVParticle* Vtrack = 0x0;
    if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
    if(fUseTender)  Vtrack  = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));
    if (!Vtrack) {
      printf("ERROR: Could not receive tagged  track %d\n", iTracks);
      continue;
    }
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(Vtrack);
    if(!track) continue;
        
    //    if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;

    // track cuts
    Bool_t passTrackCut=kFALSE;
    passTrackCut = ChargedHadronTrackCuts(pVtx, track, nMother, listMother);
    if (!passTrackCut) continue;

    // pid cuts
    Bool_t passPIDCut=kFALSE;
    passPIDCut = ChargedHadronPIDCuts(track);
    if (!passPIDCut) continue;
   
    // associated particle
    Double_t ptH=-9.,etaH =-9.,phiH=-9.;
    Int_t idH=-9;
    ptH = track->Pt();
    phiH = track->Phi();
    etaH = track->Eta();
    idH = track->GetID();
 

    CorrelateWithHadrons(track, pVtx, nMother, listMother, kTRUE, kFALSE, NonElectronIsTrigger, NonElectronIsTriggerPt, RedTracksHFE->GetEntriesFast()); // Correlate hadrons with Hadron (kTRUE, kFALSE);

    if (fMixedEvent) CorrelateHadronMixedEvent(track, mult, pVtx->GetZ(), maxPt);

    // loop over all electrons


    for (Int_t k=0; k<RedTracksHFE->GetEntriesFast(); k++) {

      AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(k);

      if (RedTrack->ID()==idH) continue;
      // tagged particle
      Double_t pt=-9.,eta =-9.,phi=-9.;
      Int_t charge = 0, ls=0, uls=0;    
      pt = RedTrack->Pt();
      phi = RedTrack->Phi();
      eta = RedTrack->Eta();
      charge = RedTrack->Charge();
      ls = RedTrack->LS();
      uls = RedTrack->ULS();

      Double_t dphi = -99, deta=-99;
      dphi = GetDeltaPhi(phi, phiH);
      deta = GetDeltaEta(eta, etaH);

      Double_t fillSparse[4]={-999,-999,-999,-999};
      fillSparse[0]=ptH;
      fillSparse[1]=pt;
      fillSparse[2]=dphi;
      fillSparse[3]=deta;
      

      CheckElectronIsTrigger(ptH, ElectronIsTrigger[k]);
      ElectronIsTriggerPt[k]=pt;
      CheckHadronIsTrigger(pt, HadronTrigger);  
      fInclElecHa->Fill(fillSparse);
      for (Int_t j=0; j<ls; j++)  fLSElecHa->Fill(fillSparse);
      for (Int_t j=0; j<uls; j++) fULSElecHa->Fill(fillSparse);
     

      // Fill Histograms with missing partner
      Bool_t HadIsULSPartner=kFALSE;
      Bool_t HadIsLSPartner=kFALSE;
      for (int j=0; j<uls; j++) {
	if (idH==RedTrack->ULSPartner(j)) HadIsULSPartner =kTRUE;
      }
      for (int j=0; j<ls; j++) {
	if (idH==RedTrack->LSPartner(j)) HadIsLSPartner=kTRUE;
      }
      if (!HadIsULSPartner && !HadIsLSPartner) {
	for (Int_t j=0; j<ls; j++)  fElecHaLSNoPartner->Fill(fillSparse);
	for (Int_t j=0; j<uls; j++) fElecHaULSNoPartner->Fill(fillSparse);
      }

      if (fIsMC) {
	if (RedTrack->IsHadron()){
	    CheckElectronIsTrigger(ptH, HadContIsTrigger[k]);
	    fMCElecHaHadron->Fill(fillSparse);
	}
	if  (RedTrack->IsPhotonic()) {
	  if (RedTrack->TruePartner()) {
	    CheckElectronIsTrigger(ptH, PhotElecWPartnerTrigger[k]);
	    fMCElecHaTruePartner->Fill(fillSparse);
	  }
	  else {
	    CheckElectronIsTrigger(ptH,  PhotElecWoPartnerTrigger[k]);
	    fMCElecHaNoPartner->Fill(fillSparse);
	  }
	}
      }
    }

    // Fill Trigger Hadrons
    for (Int_t AssPtBin=0; AssPtBin<fAssPtElec_Nbins; AssPtBin++) {
      if (HadronTrigger[AssPtBin]) fHadElecTrigger->Fill(ptH, AssPtBin);
    }

    delete [] HadronTrigger;
  }//end of track loop

  // Fill Trigger Histograms
  for (Int_t l=0; l<RedTracksHFE->GetEntriesFast(); l++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(l);
    for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
      if (ElectronIsTrigger[l][AssPtBin]) {
	fElecHadTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin);
	if( RedTrack->LS()>0) fElecHadTriggerLS->Fill(ElectronIsTriggerPt[l], AssPtBin,   RedTrack->LS());
	if( RedTrack->ULS()>0) fElecHadTriggerULS->Fill(ElectronIsTriggerPt[l], AssPtBin,  RedTrack->ULS());
      }
      if (PhotElecWPartnerTrigger[l][AssPtBin]) fMCElecHaTruePartnerTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin);
      if (PhotElecWoPartnerTrigger[l][AssPtBin]) fMCElecHaNoPartnerTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin);
      if (HadContIsTrigger[l][AssPtBin]) fHadContTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin);
    }
  }


  for (Int_t l=0; l<RedTracksHFE->GetEntriesFast(); l++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(l);
    //  cout << "Electron" << RedTrack->Pt() << "\t IsPhotonic " << RedTrack->IsPhotonic() << "\t Has Partner " << RedTrack->TruePartner() << endl;
    Bool_t NoPNoT = kTRUE;
    Bool_t TPNoT = kTRUE;
    if (RedTrack->IsPhotonic()) {
      for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
	if (PhotElecWPartnerTrigger[l][AssPtBin] && TPNoT) TPNoT=kFALSE;
	if (PhotElecWoPartnerTrigger[l][AssPtBin] && NoPNoT) NoPNoT = kFALSE;
      }
      if (RedTrack->TruePartner() && TPNoT) fTPartnerNoT->Fill(RedTrack->Pt()); //cout << "TruePartnerNoTrigger" << endl;
      if (!RedTrack->TruePartner() && NoPNoT) fNoPartnerNoT->Fill(RedTrack->Pt()); // cout << "NoPartnerNoTrigger" << endl;
    }
  }


  for (int l=0; l<ntracks; l++) {
    for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
      if (NonElectronIsTrigger[l][AssPtBin]) fNonElecHadTrigger->Fill(NonElectronIsTriggerPt[l], AssPtBin);
    }
  }

  // Clear Trigger Arrays
  for (Int_t j=0; j < RedTracksHFE->GetEntriesFast(); j++) {
    delete [] ElectronIsTrigger[j];
    delete [] HadContIsTrigger[j];
    delete [] PhotElecWPartnerTrigger[j];
    delete [] PhotElecWoPartnerTrigger[j];
  }
  delete [] ElectronIsTrigger;
  delete [] HadContIsTrigger;
  delete [] PhotElecWPartnerTrigger;
  delete [] PhotElecWoPartnerTrigger;
  delete [] ElectronIsTriggerPt;
  for (Int_t j=0; j<ntracks; j++) delete[] NonElectronIsTrigger[j];
  delete [] NonElectronIsTrigger;
  delete [] NonElectronIsTriggerPt;
 
}


void AliAnalysisTaskHaHFECorrel::CorrelateHadronMixedEvent(AliVTrack* Vtrack, Float_t mult, Float_t zVtx, Float_t maxPt) {


  Double_t ptH=-9.,etaH =-9.,phiH=-9.;
  ptH = Vtrack->Pt();
  phiH = Vtrack->Phi();
  etaH = Vtrack->Eta();

  AliEventPool* fPool;
  fPool = fPoolMgr->GetEventPool(mult, zVtx, maxPt); // Get the buffer associated with the current centrality and z-vtx
  //  fPool->SetDebug(kTRUE);
  if (!fPool)
  {
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", mult, zVtx));
    return;
  }
  
  if (fPool->GetCurrentNEvents() >= 3) // start mixing when 3 events are in the buffer
  {

    Int_t nMix = fPool->GetCurrentNEvents();   
    for (Int_t jMix=0; jMix<nMix; jMix++)  // mix with each event in the buffer
    {
      TObjArray* mixTracks = fPool->GetEvent(jMix);
      for (Int_t i=0;i<mixTracks->GetEntriesFast(); i++)
      {	
        AliBasicParticleHaHFE* mixtrk = (AliBasicParticleHaHFE*) mixTracks->At(i);
        if (!mixtrk) {
          printf("ERROR: Could not receive mix pool track %d\n",i);
          continue;
        }

	Double_t ptMix=-9., phiMix=-9., etaMix=-9;
	Int_t ls=-9, uls=-9;
	ptMix=mixtrk->Pt();
	phiMix=mixtrk->Phi();
	etaMix=mixtrk->Eta();
	ls=mixtrk->LS();
	uls=mixtrk->ULS();

	Double_t dphi=-99, deta=-99;
	dphi=GetDeltaPhi(phiMix, phiH);
	deta=GetDeltaEta(etaMix, etaH);

	Double_t fillSparse[4]={-999,-999,-999,-999};
	fillSparse[0]=ptH;
	fillSparse[1]=ptMix;
	fillSparse[2]=dphi;
	fillSparse[3]=deta;

	fElecHaMixedEvent->Fill(fillSparse);
	for (Int_t k=0; k<ls; k++)  fLSElecHaMixedEvent->Fill(fillSparse);
	for (Int_t k=0; k<uls; k++) fULSElecHaMixedEvent->Fill(fillSparse);

	if (fIsMC) {
	  if (mixtrk->IsPhotonic()) {
	    if (mixtrk->TruePartner()) {
	      fTagHaMixedEvent->Fill(fillSparse);
	    }
	    else {
	      fNonTagHaMixedEvent->Fill(fillSparse);
	    }
	  }
	}
      }
    }
  }
}


/////////////--------------

void AliAnalysisTaskHaHFECorrel::CorrelateWithHadrons(AliVTrack* Vtrack, const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], Bool_t FillHadron, Bool_t FillLP, Bool_t **NonElecIsTrigger, Double_t *NonElecIsTriggerPt, Int_t NumElectronsInEvent) {
  
  // Trigger Hadron
  Double_t ptH=-9.,etaH =-9.,phiH=-9.;
  Int_t idH=-9;
  ptH = Vtrack->Pt();
  phiH = Vtrack->Phi();
  etaH = Vtrack->Eta();
  idH = Vtrack->GetID();

  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();

  Bool_t Trigger[fAssPtElec_Nbins];
  for (Int_t AssPt=0; AssPt<fAssPtElec_Nbins; AssPt++) Trigger[AssPt]=kFALSE;           // flag to store information, if used as trigger particle
 

  // Track loop for hadron correlations
  for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {     
    
    AliVParticle* Vtrack = 0x0;
    if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
    if(fUseTender)  Vtrack  = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));
    if (!Vtrack) {
      printf("ERROR: Could not receive tagged  track %d\n", iTracks);
      continue;
    }
    AliVTrack *VtrackAssoHad = dynamic_cast<AliVTrack*>(Vtrack);
    if(!VtrackAssoHad) continue;
        
    if(VtrackAssoHad->GetID()==idH) continue;

    // track cuts
    Bool_t passTrackCut=kFALSE;
    passTrackCut = InclElecTrackCuts(pVtx, VtrackAssoHad, nMother, listMother);
    if (!passTrackCut) continue;

    // pid cuts
    Bool_t passPIDCut=kFALSE;
    passPIDCut = AssoHadronPIDCuts(VtrackAssoHad);
    if (!passPIDCut) continue;
   
    // associated particle
    Double_t ptAsso=-9.,etaAsso=-9.,phiAsso=-9.;
    ptAsso = VtrackAssoHad->Pt();
    phiAsso = VtrackAssoHad->Phi();
    etaAsso = VtrackAssoHad->Eta();
  
 
    CheckHadronIsTrigger(ptAsso, Trigger);
    CheckElectronIsTrigger(ptH, NonElecIsTrigger[iTracks]);
    NonElecIsTriggerPt[iTracks]=ptAsso;
    if(FillLP) {
      //  cout << ptH;
      for (Int_t ptAssHad=0; ptAssHad<fAssPtHad_Nbins; ptAssHad++) {
	if (NonElecIsTrigger[iTracks][ptAssHad]) fNonElecLPTrigger->Fill(ptAsso, ptAssHad);
	//	cout << NonElecIsTrigger[iTracks][ptAssHad];
      }
      //   cout << endl;
    }
  

    Double_t dphi = -99, deta=-99;
    dphi = GetDeltaPhi(phiAsso, phiH);
    deta = GetDeltaEta(etaAsso, etaH);

    // Fill Sparse
    Double_t fillSparse[4]={-999,-999,-999,-999};
    fillSparse[0]=ptH;
    fillSparse[1]=ptAsso;
    fillSparse[2]=dphi;
    fillSparse[3]=deta;
    if (FillHadron) fElecHaHa->Fill(fillSparse);
    if (FillLP)     fElecLPHa->Fill(fillSparse);
  }

  if (FillHadron) {
    for (Int_t ptAssElec=0; ptAssElec<fAssPtElec_Nbins; ptAssElec++) {
      if (Trigger[ptAssElec]) fHadNonElecTrigger->Fill(ptH, ptAssElec);
    }
  }
  if (FillLP){
    for (Int_t ptAssElec=0; ptAssElec<fAssPtElec_Nbins; ptAssElec++) {
      if (Trigger[ptAssElec]) fLPNonElecTrigger->Fill(ptH, ptAssElec);
    }
  }
}



//_______________________

void AliAnalysisTaskHaHFECorrel::CorrelateLP(AliVTrack* LPtrack,  const AliVVertex* pVtx, Int_t nMother, Int_t listMother[], TObjArray *RedTracksHFE) 
{ 

  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();

  Bool_t *LPTrigger= new Bool_t[fAssPtElec_Nbins];
  for (Int_t l=0; l<fAssPtElec_Nbins; l++) LPTrigger[l]=kFALSE;
  Bool_t **ElectronIsTrigger = new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **HadContIsTrigger=new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **PhotElecWPartnerTrigger=new Bool_t*[RedTracksHFE->GetEntriesFast()];
  Bool_t **PhotElecWoPartnerTrigger=new Bool_t*[RedTracksHFE->GetEntriesFast()];
  for (Int_t j=0; j<RedTracksHFE->GetEntriesFast(); j++)
    {
      ElectronIsTrigger[j]=new Bool_t[fAssPtHad_Nbins];
      HadContIsTrigger[j]=new Bool_t[fAssPtHad_Nbins];
      PhotElecWPartnerTrigger[j] = new Bool_t[fAssPtHad_Nbins];
      PhotElecWoPartnerTrigger[j] = new Bool_t[fAssPtHad_Nbins];
      for (Int_t m=0; m<fAssPtHad_Nbins; m++) {
	ElectronIsTrigger[j][m]=kFALSE;
	HadContIsTrigger[j][m]=kFALSE;
	PhotElecWPartnerTrigger[j][m]=kFALSE;
	PhotElecWoPartnerTrigger[j][m]=kFALSE;
      }
    }
  Double_t *ElectronIsTriggerPt = new Double_t[RedTracksHFE->GetEntriesFast()];
  Bool_t **NonElectronIsTrigger = new Bool_t*[ntracks];
  for (Int_t j=0; j<ntracks; j++) { 
    NonElectronIsTrigger[j]=new Bool_t[fAssPtHad_Nbins]; 
    for (Int_t m=0; m<fAssPtHad_Nbins; m++) {
      NonElectronIsTrigger[j][m]=kFALSE;
    }
  }
  Double_t *NonElectronIsTriggerPt = new Double_t[ntracks];

  // leading Particle neq 
  Double_t ptH=-9.,etaH =-9.,phiH=-9.;
  Int_t idH=-9;
  ptH = LPtrack->Pt();
  phiH = LPtrack->Phi();
  etaH = LPtrack->Eta();
  idH = LPtrack->GetID();

  CorrelateWithHadrons(LPtrack, pVtx, nMother, listMother, kFALSE, kTRUE, NonElectronIsTrigger, NonElectronIsTriggerPt, RedTracksHFE->GetEntriesFast()); // correlate LPHadron (kFALSE, kTRUE);
		       
  // Only loop over electrons in event
  for (Int_t k=0; k<RedTracksHFE->GetEntriesFast(); k++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(k);
  


    if (RedTrack->IsHadron()) {
      cout << "RedTrackIsHadron" << endl;
      cout << RedTrack->ID() << "\t" << idH << endl;
    }


    if (RedTrack->ID()==idH) continue; // no self-correlation
    Double_t pt=-9.,eta =-9.,phi=-9.;
    Int_t ls=0, uls=0;    
    pt = RedTrack->Pt();
    phi = RedTrack->Phi();
    eta = RedTrack->Eta();
    ls = RedTrack->LS();
    uls = RedTrack->ULS();

    Double_t dphi = -99, deta = -99;
    dphi = GetDeltaPhi(phi, phiH);
    deta = GetDeltaEta(eta, etaH);

    Double_t fillSparse[4]={-999,-999,-999,-999};
    fillSparse[0]=ptH;
    fillSparse[1]=pt;
    fillSparse[2]=dphi;
    fillSparse[3]=deta;

    CheckElectronIsTrigger(ptH, ElectronIsTrigger[k]);
    ElectronIsTriggerPt[k]=pt;
    CheckHadronIsTrigger(pt, LPTrigger);
    fInclElecLP->Fill(fillSparse);
    for (Int_t j=0; j<ls; j++) fLSElecLP->Fill(fillSparse);
    for (Int_t j=0; j<uls; j++) fULSElecLP->Fill(fillSparse);


    Bool_t HadIsULSPartner=kFALSE;
    Bool_t HadIsLSPartner=kFALSE;
    for (int j=0; j<uls; j++) {
      if (idH==RedTrack->ULSPartner(j)) HadIsULSPartner =kTRUE;
    }
    for (int j=0; j<ls; j++) {
      if (idH==RedTrack->LSPartner(j)) HadIsLSPartner=kTRUE;
    }
    if (!HadIsULSPartner && !HadIsLSPartner) {
      for (Int_t j=0; j<ls; j++) fElecLPLSNoPartner->Fill(fillSparse);
      for (Int_t j=0; j<uls; j++) fElecLPULSNoPartner->Fill(fillSparse);
    }

    if (fIsMC) {
      if (RedTrack->IsHadron()){
	CheckElectronIsTrigger(ptH, HadContIsTrigger[k]);
	fMCElecLPHadron->Fill(fillSparse);
      }
      if (RedTrack->IsPhotonic()) {
	if (RedTrack->TruePartner()) {
	  CheckElectronIsTrigger(ptH, PhotElecWPartnerTrigger[k]);
	  fMCElecLPTruePartner->Fill(fillSparse);
	}
	else {
	  CheckElectronIsTrigger(ptH,  PhotElecWoPartnerTrigger[k]);
	  fMCElecLPNoPartner->Fill(fillSparse);
	}
      }
    }
  }
   
  //save LP as Trigger (existens of LP and HFE are satisfied in function call)
  for (Int_t PtAssElec=0; PtAssElec<fAssPtElec_Nbins; PtAssElec++) {
    if (LPTrigger[PtAssElec]) fLPElecTrigger->Fill(ptH, PtAssElec);
  }


  // Fill Trigger Histograms
  for (Int_t l=0; l<RedTracksHFE->GetEntriesFast(); l++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(l);
    for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
      if (ElectronIsTrigger[l][AssPtBin]) {
	fElecLPTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin);
	if( RedTrack->LS()>0) fElecLPTriggerLS->Fill(ElectronIsTriggerPt[l], AssPtBin,   RedTrack->LS());
	if( RedTrack->ULS()>0) fElecLPTriggerULS->Fill(ElectronIsTriggerPt[l], AssPtBin,  RedTrack->ULS());
      }
      if (PhotElecWPartnerTrigger[l][AssPtBin]) fMCElecLPTruePartnerTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin);
      if (PhotElecWoPartnerTrigger[l][AssPtBin]) fMCElecLPNoPartnerTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin);
      if (HadContIsTrigger[l][AssPtBin]) fHadContLPTrigger->Fill(ElectronIsTriggerPt[l], AssPtBin);
    }
  }

  for (int l=0; l<ntracks; l++) {
    for (Int_t AssPtBin=0; AssPtBin<fAssPtHad_Nbins; AssPtBin++) {
      if (NonElectronIsTrigger[l][AssPtBin]) fNonElecLPTrigger->Fill(NonElectronIsTriggerPt[l], AssPtBin);
    }
  }


  // Clear Trigger Arrays
  for (Int_t j=0; j < RedTracksHFE->GetEntriesFast(); j++) {
    delete [] ElectronIsTrigger[j];
    delete [] HadContIsTrigger[j];
    delete [] PhotElecWPartnerTrigger[j];
    delete [] PhotElecWoPartnerTrigger[j];
  }
  delete [] ElectronIsTrigger;
  delete [] HadContIsTrigger;
  delete [] PhotElecWPartnerTrigger;
  delete [] PhotElecWoPartnerTrigger;
  delete [] LPTrigger;
  delete [] ElectronIsTriggerPt;
  for (Int_t j=0; j<ntracks; j++) delete[] NonElectronIsTrigger[j];
  delete [] NonElectronIsTrigger;
  delete [] NonElectronIsTriggerPt;

}



//_______________________

void AliAnalysisTaskHaHFECorrel::CorrelateLPMixedEvent(AliVTrack* LPtrack, Float_t mult, Float_t zVtx, Float_t maxPt ) 
{ 
  // leading Particle 
  Double_t ptH=-9.,etaH =-9.,phiH=-9.;
  ptH = LPtrack->Pt();
  phiH = LPtrack->Phi();
  etaH = LPtrack->Eta();

  AliEventPool* fPool;
  fPool = fPoolMgr->GetEventPool(mult, zVtx, maxPt); // Get the buffer associated with the current centrality and z-vtx
  //  fPool->SetDebug(kTRUE);
  if (!fPool)
  {
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", mult, zVtx));
    return;
  }
  
  if (fPool->GetCurrentNEvents() >= 3) // start mixing when 3 events are in the buffer
  {
    Int_t nMix = fPool->GetCurrentNEvents();   
    for (Int_t jMix=0; jMix<nMix; jMix++)  // mix with each event in the buffer
    {
      TObjArray* mixTracks = fPool->GetEvent(jMix);
      for (Int_t i=0;i<mixTracks->GetEntriesFast(); i++)
      {	
        AliBasicParticleHaHFE* mixtrk = (AliBasicParticleHaHFE*) mixTracks->At(i);
	if (!mixtrk) {
          printf("ERROR: Could not receive mix pool track %d\n",i);
          continue;
        }
	
	Double_t ptMix=-9., phiMix=-9., etaMix=-9;
	Int_t ls=-9, uls=-9;
	ptMix=mixtrk->Pt();
	phiMix=mixtrk->Phi();
	etaMix=mixtrk->Eta();
	ls=mixtrk->LS();
	uls=mixtrk->ULS();

	Double_t dphi=-99, deta=-99;
	dphi=GetDeltaPhi(phiMix, phiH);
	deta=GetDeltaEta(etaMix, etaH);

	Double_t fillSparse[4]={-999,-999,-999,-999};
	fillSparse[0]=ptH;
	fillSparse[1]=ptMix;
	fillSparse[2]=dphi;
	fillSparse[3]=deta;

	fElecLPMixedEvent->Fill(fillSparse);
	for (Int_t k=0; k<ls; k++) fLSElecLPMixedEvent->Fill(fillSparse);
	for (Int_t k=0; k<uls; k++) fULSElecLPMixedEvent->Fill(fillSparse);
      }
    }
  }
}


//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::InclElecTrackCuts(const AliVVertex *pVietx,AliVTrack *Vtrack,Int_t nMother, Int_t listMother[])
{
  AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);   
  // track cuts for inclusive electrons
  Double_t d0z0[2]={-999,-999}, cov[3];
  fElectronTrackCuts->Fill(Vtrack->Pt(), 0);

 
  fElectronTrackTPCNcls->Fill(Vtrack->Pt(), Vtrack->GetTPCNcls());
  fElectronTrackTPCNclsdEdx->Fill(Vtrack->Pt(), Vtrack->GetTPCsignalN());
  if (fIsAOD) { 
    fElectronTrackTPCFrac->Fill(AODtrack->Pt(), AODtrack->GetTPCFoundFraction());
    fElectronTrackITSNcls->Fill(Vtrack->Pt(), AODtrack->GetITSNcls());
  }
  //  if(ietrack->PropagateToDCA(pVietx, fVevent->GetMagneticField(), 20., d0z0, cov))
  //fElectronTrackDCA->Fill(TMath::Abs(d0z0[0]), TMath::Abs(d0z0[1]));

  if (fIsAOD) {
    if(!AODtrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE;
  }
  else if(!ESDkTrkGlobalNoDCA(Vtrack)) return kFALSE;

  fElectronTrackCuts->Fill(Vtrack->Pt(), 1);

  if(Vtrack->Pt()<0.25) return kFALSE;
  fElectronTrackCuts->Fill(Vtrack->Pt(), 2);

  if(Vtrack->Eta()<fMinElectronEta) return kFALSE;
  if(Vtrack->Eta()>fMaxElectronEta) return kFALSE;
  fElectronTrackCuts->Fill(Vtrack->Pt(), 3);

  if(Vtrack->GetTPCNcls() < fTPCnCut) return kFALSE;
  fElectronTrackCuts->Fill(Vtrack->Pt(), 4);
  if(Vtrack->GetTPCsignalN() < fTPCndEdxCut) return kFALSE ;
  fElectronTrackCuts->Fill(Vtrack->Pt(), 5);
						
  if (fIsAOD) if(!AODtrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  fElectronTrackCuts->Fill(Vtrack->Pt(), 6);
  //if(ietrack->GetTPCFoundFraction() < 0.6) return kFALSE;   
  fElectronTrackCuts->Fill(Vtrack->Pt(), 7);

  if (fIsAOD) {
    Bool_t kinkmotherpass = kTRUE;
    for(Int_t kinkmother = 0; kinkmother < nMother; kinkmother++) {
      if(Vtrack->GetID() == listMother[kinkmother]) {
	kinkmotherpass = kFALSE;
	continue;
      }
    }
    if(!kinkmotherpass) return kFALSE;
  }
  else {
    AliESDtrack *Htrack = dynamic_cast<AliESDtrack*>(Vtrack);   
    if(!Htrack) return kFALSE;
    if(Htrack->GetKinkIndex(0) != 0) return kFALSE;
  }
  fElectronTrackCuts->Fill(Vtrack->Pt(), 8);

  if (fIsAOD){
    if(AODtrack->GetITSNcls() < fITSnCut) return kFALSE;
    fElectronTrackCuts->Fill(Vtrack->Pt(), 9);

    if(!AODtrack->IsOn(AliAODTrack::kITSrefit)) return kFALSE;    
    fElectronTrackCuts->Fill(Vtrack->Pt(), 10);
    
    if(!(AODtrack->HasPointOnITSLayer(0) || AODtrack->HasPointOnITSLayer(1))) return kFALSE; // now kBoth
    fElectronTrackCuts->Fill(Vtrack->Pt(), 11);

    if(!(AODtrack->HasPointOnITSLayer(0) && AODtrack->HasPointOnITSLayer(1))) return kFALSE; // now kBoth
    fElectronTrackCuts->Fill(Vtrack->Pt(), 12);
  
    if(AODtrack->PropagateToDCA(pVietx, fVevent->GetMagneticField(), 20., d0z0, cov))
    if(TMath::Abs(d0z0[0]) > 1 || TMath::Abs(d0z0[1]) > 2) return kFALSE;
    fElectronTrackCuts->Fill(Vtrack->Pt(), 13);

  }

 

  return kTRUE;

 }

 //________________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::InclElecPIDCuts(AliVTrack* Vtrack, Bool_t IsPrimary) {

  AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);   
  // PID cuts
  Double_t fITSnSigma=-9.,fTOFnSigma=-9.,fTPCnSigma=-9.;
  Bool_t PassITSCut=kTRUE, PassTPCCut=kTRUE, PassTOFCut=kTRUE;

  fITSnSigma = fpidResponse->NumberOfSigmasITS(Vtrack, AliPID::kElectron);
  fTPCnSigma = fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
  fTOFnSigma = fpidResponse->NumberOfSigmasTOF(Vtrack, AliPID::kElectron);
  
  
  if ((fITSnSigma < (-1.*fSigmaITScut) ) || (fITSnSigma > fSigmaITScut)) PassITSCut=kFALSE;
  if ((fTOFnSigma < (-1.*fSigmaTOFcut) ) || (fTOFnSigma > fSigmaTOFcut)) PassTOFCut=kFALSE;
  if ((fTPCnSigma<fSigmaTPCcut)  || (fTPCnSigma>3)) PassTPCCut=kFALSE; 


  if (IsPrimary) { //Fill PID histograms only in first round
    fHistITSnSig->Fill(Vtrack->P(),fITSnSigma);
    fHistTOFnSig->Fill(Vtrack->P(),fTOFnSigma);
    fHistTPCnSig->Fill(Vtrack->P(),fTPCnSigma);
    
    if (PassITSCut) fHistTPCnSigITScut->Fill(Vtrack->P(),fTPCnSigma);
    if (PassTOFCut) fHistTPCnSigTOFcut->Fill(Vtrack->P(),fTPCnSigma);
    if (PassITSCut && PassTOFCut) fHistTPCnSigITSTOFcut->Fill(Vtrack->P(),fTPCnSigma);

    if (PassTOFCut) fHadContPvsPt->Fill(Vtrack->P(), Vtrack->Pt());

    Double_t fillSparse[4];   
    fillSparse[0]=Vtrack->P();
    fillSparse[1]=Vtrack->Phi();
    fillSparse[2]=Vtrack->Eta();
    fillSparse[3]=fTPCnSigma;

    if (PassTOFCut) fHadContPPhiEtaTPC->Fill(fillSparse);

    fillSparse[1]=fITSnSigma;
    fillSparse[2]=fTOFnSigma;
    fHadContamination->Fill(fillSparse);

    fillSparse[0]=Vtrack->Pt();
    fHadContaminationPt->Fill(fillSparse);

    if (fIsAOD && PassTOFCut && fIsMC) {

      Int_t MClabel=AODtrack->GetLabel();
      AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));     
      Int_t PDGCode = abs(MCParticle->GetPdgCode());
      fillSparse[0]=AODtrack->P();
      if (PDGCode==11) fillSparse[1]=0; // electron
      else if (PDGCode==13) fillSparse[1]=1; // muon
      else if (PDGCode==111 || PDGCode==211) fillSparse[1]=2; // p0 and p+
      else if (PDGCode==130 || PDGCode==310 || PDGCode==311 || PDGCode==321) fillSparse[1]=3; // K0L KOS K0 K0+
      else if (PDGCode==2212) fillSparse[1]=4; // proton
      else if (PDGCode==1000010020) fillSparse[1]=5;
      else fillSparse[1]=6;
      fillSparse[2]=fITSnSigma;
      fHadContMC->Fill(fillSparse);
      fillSparse[0]=AODtrack->Pt();
      fHadContMCPt->Fill(fillSparse);
    }
   }
 
  if (!fUseITS) PassITSCut = kTRUE;
  if (!PassITSCut || !PassTOFCut || !PassTPCCut) return kFALSE;
  return kTRUE;
}




// _________________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::PhotElecPIDCuts(AliVTrack* Vtrack) {
  // associated particle variables
  Double_t fITSnSigmaAsso=-9.,fTPCnSigmaAsso=-9.;
  
  // looser PID cuts
  fITSnSigmaAsso = fpidResponse->NumberOfSigmasITS(Vtrack, AliPID::kElectron);
  fTPCnSigmaAsso = fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);

  if(TMath::Abs(fTPCnSigmaAsso)>fPhotElecSigmaTPCcut) return kFALSE;
  // if(TMath::Abs(fITSnSigmaAsso)>3) return kFALSE;

  return kTRUE;
  
}


//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::PhotElecTrackCuts(const AliVVertex *pVaetx,AliVTrack *Vtrack,Int_t nMother, Int_t listMother[])
{
  if (fIsAOD){
    AliAODTrack *aetrack = dynamic_cast<AliAODTrack*>(Vtrack);   
    // quality track cuts for associate tracks of photonic electrons
    if(!aetrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) return kFALSE;
  }
  if(Vtrack->Pt() < fPhotElecPtCut) return kFALSE;
  if(TMath::Abs(Vtrack->Eta())>0.9) return kFALSE;
  if(Vtrack->GetTPCNcls() < fPhotElecTPCnCut) return kFALSE;
  if (fPhotElecITSrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
  if(!(Vtrack->GetStatus()&AliESDtrack::kTPCrefit)) return kFALSE;

  return kTRUE;
  
}

 //_________________________________________

Bool_t AliAnalysisTaskHaHFECorrel::ChargedHadronTrackCuts(const AliVVertex *pVaetx,AliVTrack *Vtrack,Int_t nMother, Int_t listMother[])
 {
   // quality track cuts for charged hadrons
   if (fIsAOD) {
     AliAODTrack *Htrack = dynamic_cast<AliAODTrack*>(Vtrack);   
     if (!Htrack) return kFALSE;
     if(!Htrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE; 

     Bool_t kinkmotherpass = kTRUE;
     for(Int_t kinkmother = 0; kinkmother < nMother; kinkmother++) {
       if(Htrack->GetID() == listMother[kinkmother]) {
	 kinkmotherpass = kFALSE;
	 continue;
       }
     }
     if(!kinkmotherpass) return kFALSE;
   }
   else {
     AliESDtrack *Htrack = dynamic_cast<AliESDtrack*>(Vtrack);   
     if(!Htrack) return kFALSE;
     if(Htrack->GetKinkIndex(0) != 0) return kFALSE;
     else if(!ESDkTrkGlobalNoDCA(Vtrack)) return kFALSE;
   }

   if(Vtrack->Pt()<0.5) return kFALSE;
   if(Vtrack->Eta()>fMaxHadronEta) return kFALSE;
   if(Vtrack->Eta()<fMinHadronEta) return kFALSE;


  if(fHITSrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
  if(fHTPCrefitCut && !(Vtrack->GetStatus()&AliESDtrack::kTPCrefit)) return kFALSE;

  if(Vtrack->GetTPCNcls() < fHTPCnCut) return kFALSE;

  return kTRUE;   
}


//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::ChargedHadronPIDCuts(AliVTrack *Vtrack)
{   
  return kTRUE;   
}


Bool_t AliAnalysisTaskHaHFECorrel::AssoHadronPIDCuts(AliVTrack *Vtrack) 
{
  // associated particle variables
  Double_t fTPCnSigmaAsso=-9.;
  Double_t fTOFnSigmaAsso=-9.;
  


  fTOFnSigmaAsso = fpidResponse->NumberOfSigmasTOF(Vtrack, AliPID::kElectron);
  // if ((fTOFnSigmaAsso < (-1.*fSigmaTOFcut) ) || (fTOFnSigmaAsso > fSigmaTOFcut)) return kFALSE;
  // looser PID cuts
  fTPCnSigmaAsso = fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
  if(fTPCnSigmaAsso>fAssNonEleTPCcut) return kFALSE;
   
  return kTRUE;
}

void AliAnalysisTaskHaHFECorrel::MCEfficiencyCorrections(const AliVVertex * RecVertex) {
  for (Int_t i=0; i<fMC->GetNumberOfTracks(); i++) {
    AliAODMCParticle *mcPart  = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(i));
    Double_t mcPt, mcPhi, mcEta, mcVtx[3];
    Int_t mcPDG, mcMotherID;
    mcPt=mcPart->Pt();
    mcPhi=mcPart->Phi();
    mcEta=mcPart->Eta();
    mcPDG=abs(mcPart->GetPdgCode());
    mcMotherID=mcPart->GetMother();
    
    fMCheader->GetVertex(mcVtx);
    fCheckMCVertex->Fill(mcVtx[2], RecVertex->GetZ()); // TH2D

    // PI0 and Eta
    if (mcPart->IsPrimary()) {
      Double_t fillSparse[5];
      if (mcPDG==111 || mcPDG==221 || mcPDG== 211) {
	if (mcMotherID<=0) {
	  fillSparse[4]=1;
	  fillSparse[3]=-1;
	}
	else {
	  AliAODMCParticle* mcPartMother=(AliAODMCParticle*)fMC->GetTrack(mcMotherID);
	  Int_t mcMotherPDG = abs(mcPartMother->GetPdgCode());
	  fillSparse[4]=0;
	  fillSparse[3]=mcMotherPDG;
	}
	fillSparse[0]=mcPt;
	fillSparse[1]=mcEta;
	fillSparse[2]=Eta2y(mcPt, mcPart->M(), mcEta);
	
	if (mcPDG ==111) fMCPi0Prod->Fill(fillSparse);
	else if (mcPDG==221) fMCEtaProd->Fill(fillSparse);
	else if (mcPDG==211) fMCPiPlusProd->Fill(fillSparse);
      }
    }


    if (mcPart->IsPhysicalPrimary() && mcPDG==211) {
      if (!mcPart->IsSecondaryFromMaterial() && ! mcPart->IsSecondaryFromWeakDecay()) {
	Double_t fillSparse[5];
	if (mcMotherID<=0) {
	  fillSparse[4]=1;
	  fillSparse[3]=-1;
	}
	else {
	  AliAODMCParticle* mcPartMother=(AliAODMCParticle*)fMC->GetTrack(mcMotherID);
	  Int_t mcMotherPDG = abs(mcPartMother->GetPdgCode());
	  fillSparse[4]=0;
	  fillSparse[3]=mcMotherPDG;
	}
	fillSparse[0]=mcPt;
	fillSparse[1]=mcEta;
	fillSparse[2]=Eta2y(mcPt, mcPart->M(), mcEta);
	fMCPiPlusProdV2->Fill(fillSparse);
      }
    }

    
    // test chec for MC Generator 
    TString genname;
    Bool_t GetCocktailHeader=fMC->GetCocktailGenerator(i,genname);
    if(!GetCocktailHeader) Printf("no cocktail header list was found for this event");
    if(GetCocktailHeader) {
      //      Printf("cocktail header name is %s", genname.Data()); 
      //you may want to check wether it is Hijing, for example.
      //   if(genname.Contains("Hijing")) Printf("this particle comes from HIJING");
    }
    
    
    // HadElectronRecEff
    if (mcPart->Charge()==0) continue;    
    
    // Fill Histogram for Hadron and Electron Reconstruction
    if (mcPart->IsPhysicalPrimary()) {
      Double_t fillSparse[4];
      fillSparse[0]=mcPt;
      fillSparse[1]=mcEta;
      fillSparse[2]=mcPhi;
      fillSparse[3]=mcVtx[2];

      if (mcEta> fMinHadronEta && mcEta < fMaxHadronEta) {
	fMCHadPtEtaPhiVtx->Fill(fillSparse);
      }
      
      if (abs(mcPDG)==11) {
	
	if (mcMotherID>0) {
	  AliAODMCParticle* mcPartMother=(AliAODMCParticle*)fMC->GetTrack(mcMotherID);
	  Int_t mcMotherPDG = abs(mcPartMother->GetPdgCode()); 
	  Int_t MIsHeavy  = Int_t (mcMotherPDG / TMath::Power(10, Int_t(TMath::Log10(mcMotherPDG))));
	  if (MIsHeavy>3) {
	    fillSparse[2]=mcPhi;
	    fMCElecPtEtaPhiVtx->Fill(fillSparse);
	    fMCElecPDG->Fill(mcMotherPDG);
	    if (mcMotherPDG==411 || mcMotherPDG==421 || mcMotherPDG==431 || mcMotherPDG==511 || mcMotherPDG ==521 || mcMotherPDG ==531 || mcMotherPDG==541) {
	      if (mcEta> fMinElectronEta && mcEta < fMaxElectronEta) {
		 fMCElecPtEtaPhiStrictVtx->Fill(fillSparse);
	      }
	    }
	  } 
	}
      }     
    }
  }
}



void AliAnalysisTaskHaHFECorrel::EvaluateTaggingEfficiency(AliVTrack * Vtrack, Int_t LSPartner, Int_t ULSPartner, Bool_t  trueULSPartner) {
  if (fIsAOD) {
    AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);   
    if (!AODtrack) return;
    Int_t MClabel=AODtrack->GetLabel();
    AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));     
    Int_t PDGCode = abs(MCParticle->GetPdgCode());
    Int_t PDGCodeMother=0;
    Int_t PDGCodeGrandMother=0;
    Int_t PDGCodeGGMother=0;
    Double_t PtMother=0; 
    Double_t PtGrandMother=0;

    // get ancestors (eta-pi0-e, eta-gamma-e, pi0-gamma-e, eta-e, pi0-e, ...) - remember pt of pi0/Eta for correction
    if (MCParticle->GetMother()>=0) {
      AliAODMCParticle* MCParticleMother = (AliAODMCParticle*) fMC->GetTrack(MCParticle->GetMother());
      PDGCodeMother = abs(MCParticleMother->GetPdgCode());
      PtMother = MCParticleMother->Pt();
      if (MCParticleMother->GetMother()>=0) {
	AliAODMCParticle* MCParticleGrandMother=(AliAODMCParticle*) fMC->GetTrack(MCParticleMother->GetMother());
	PDGCodeGrandMother=abs(MCParticleGrandMother->GetPdgCode());
	PtGrandMother = MCParticleGrandMother->Pt();
	if (MCParticleGrandMother->GetMother()>=0) {
	  AliAODMCParticle* MCParticleGGMother = (AliAODMCParticle*) fMC->GetTrack(MCParticleGrandMother->GetMother());
	  PDGCodeGGMother=abs(MCParticleGGMother->GetPdgCode());
	}
	else{
	  //	cout << "No Mother? " << MCParticleGrandMother->GetMother() << endl;
	  PDGCodeGGMother=-1;
	}
      }
      else {
	PDGCodeGrandMother=-1;
	PDGCodeGGMother=-1;
      } 
    }
    else{
      PDGCodeMother=-1;
      PDGCodeGrandMother=-1;
      PDGCodeGGMother=-1;
    } 
  
 

    Double_t fillSparse[5]={-999,-999,-999, -999, -999};
    fillSparse[0]=Vtrack->Pt(); // sparse will be filled with rec pt 
    fillSparse[1]=PDGCode;
    fillSparse[2]=PDGCodeMother;
    fillSparse[3]=PDGCodeGrandMother; 
    fillSparse[4]=PDGCodeGGMother;
  
 
    Double_t PtMotherWeight=1.; 
    // Double_t ExcludePion[9] = {130,211, 310, 321, 2212, 3122,3222,3322,5122}; // K0L, PI+, K0S, K+, p, Lambd0, Sigma+, Xci0, lambda b0
    // Double_t ExcludeEta[0] = {};
    if (PDGCode==11) {
      if (PDGCodeMother==22) {
	if (PDGCodeGrandMother == 221) {
	  PtMotherWeight=fCorrectEtatoData.GetBinContent(fCorrectEtatoData.FindBin(PtGrandMother));
	}
	else if (PDGCodeGrandMother == 111) {
	  PtMotherWeight=fCorrectPiontoData.GetBinContent(fCorrectPiontoData.FindBin(PtGrandMother));
	  //	for (Int_t j=0; j<9; j++) if (PDGCodeGGMother==ExcludePion[j]) PtMotherWeight=1.;
	}
      }
      else  {
	if (PDGCodeMother == 221) {
	  PtMotherWeight=fCorrectEtatoData.GetBinContent(fCorrectEtatoData.FindBin(PtMother));
	}
	else if (PDGCodeMother == 111) {
	  PtMotherWeight=fCorrectPiontoData.GetBinContent(fCorrectPiontoData.FindBin(PtMother));
	  //for (Int_t j=0; j<9; j++) if (PDGCodeGGMother==ExcludePion[j]) PtMotherWeight=1.;
	}
      }
    }
  
    // real tag corr sparse
    fTagEffIncl->Fill(fillSparse, PtMotherWeight);
    if (LSPartner>0)  fTagEffLS->Fill(fillSparse, LSPartner*PtMotherWeight);
    if (ULSPartner>0)  fTagEffULS->Fill(fillSparse, ULSPartner*PtMotherWeight);
    if (trueULSPartner) fTagTruePairs->Fill(fillSparse, PtMotherWeight);

    fTagEffInclWoWeight->Fill(fillSparse);
    if (LSPartner>0)  fTagEffLSWoWeight->Fill(fillSparse, LSPartner);
    if (ULSPartner>0)  fTagEffULSWoWeight->Fill(fillSparse, ULSPartner);
    if (trueULSPartner) fTagTruePairsWoWeight->Fill(fillSparse);


    // extra sparse to control pt mother, pt electron, ratio of primary to secondary
    if (PDGCode==11) {
      if (PDGCodeMother==22) {
	fillSparse[1]=PtGrandMother;
	fillSparse[2]=PDGCodeGrandMother;
	fillSparse[3]=PDGCodeGGMother;
      }
      else {
	fillSparse[1]=PtMother;
	fillSparse[2]=PDGCodeMother;
	fillSparse[3]=PDGCodeGrandMother;
      }
      fTagMotherPt->Fill(fillSparse, PtMotherWeight);
    }
  }
}




//_________________________________
Bool_t  AliAnalysisTaskHaHFECorrel::CloneAndReduceTrackList(TObjArray* RedTracks, AliVTrack* Vtrack, Int_t LSPartner, Int_t ULSPartner, Int_t* LSPartnerID, Int_t* ULSPartnerID, Bool_t trueULSPartner, Bool_t isPhotonic, Bool_t isHadron) {
  
  // Copied form AliAnalysisTaksPhiCorrelations and following instructions on
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PAGCorrelationsFAQ#In_Memory_Event_Mixing
  // Main difference is that track selection is in main to reduce computation time


  fCheckLSULS->Fill(LSPartner, ULSPartner);

  //  cout << "ParticleAdded" << endl;

  // clones a track list by using AliBasicParticle which uses much less memory (used for event mixing) -- modified class for my use
  // Clone only a certain pt bin on demand
  AliBasicParticleHaHFE * bparticle  = 0;
  bparticle = new AliBasicParticleHaHFE(Vtrack->GetID(), Vtrack->Eta(), Vtrack->Phi(), Vtrack->Pt(), Vtrack->Charge(), LSPartner, ULSPartner, LSPartnerID, ULSPartnerID, trueULSPartner, isPhotonic, isHadron);
  if (!bparticle) return kFALSE;
  RedTracks->Add(bparticle);
  return kTRUE;

}

Double_t AliAnalysisTaskHaHFECorrel::Eta2y(Double_t pt, Double_t m, Double_t eta) const
{
  // convert eta to y
  Double_t mt = TMath::Sqrt(m * m + pt * pt);
  return TMath::ASinH(pt / mt * TMath::SinH(eta));
}

void AliAnalysisTaskHaHFECorrel::BinLogX(TAxis *axis)
{
  // Method for the correct logarithmic binning of histograms
  const Int_t bins    = axis->GetNbins();
  const Double_t from = axis->GetXmin();
  const Double_t to   = axis->GetXmax();
    
  if (from<10e-9){
    //printf("BinLogX warning xmin < epsilon! nothing done, axis not set. %e\n", from);  exit(1);
    return;
  }
    Double_t * new_bins = new Double_t[bins + 1];
  new_bins[0] = from;
  const Double_t factor = pow(to/from, 1./bins);
  for (int i = 1; i <= bins; i++) {
    new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete [] new_bins;
}

Bool_t AliAnalysisTaskHaHFECorrel::HaveSameMother(Int_t Label1, Int_t Label2) const
{
  AliAODMCParticle* MCParticle1 = (AliAODMCParticle*) fMC->GetTrack(abs(Label1));  
  AliAODMCParticle* MCParticle2 = (AliAODMCParticle*) fMC->GetTrack(abs(Label2));
  if (MCParticle1->GetMother() ==  MCParticle2->GetMother()) return kTRUE;
  else return kFALSE;
}

Bool_t AliAnalysisTaskHaHFECorrel::IsPhotonicElectron(Int_t Label1) const
{
  Int_t PDGCodeMother, PDGCodeGrandMother;
  AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(Label1));
  if (MCParticle->GetMother()>=0) {
     AliAODMCParticle* MCParticleMother = (AliAODMCParticle*) fMC->GetTrack(MCParticle->GetMother());
     PDGCodeMother = abs(MCParticleMother->GetPdgCode());
     if (PDGCodeMother == 111 || PDGCodeMother == 221) return kTRUE;
     else if (PDGCodeMother == 22 && MCParticleMother->GetMother()>=0) {
       AliAODMCParticle* MCParticleGrandMother=(AliAODMCParticle*) fMC->GetTrack(MCParticleMother->GetMother());
       PDGCodeGrandMother=abs(MCParticleGrandMother->GetPdgCode());
        if (PDGCodeGrandMother == 111 || PDGCodeGrandMother == 221) return kTRUE;
	else return kFALSE;
     }
     else return kFALSE;
  }
  else return kFALSE;
}


void AliAnalysisTaskHaHFECorrel::FindV0CandidatesESD(AliESDEvent *event) 
{

  fV0cutsESD->SetMode(AliESDv0KineCuts::kEffGamma,AliESDv0KineCuts::kPP);
  fV0cutsESD->SetEvent(event);

  const Int_t numTracks = event->GetNumberOfTracks();
  fV0tags = new Int_t[numTracks];
  for (Int_t i = 0; i < numTracks; i++) fV0tags[i] = 0;

  // loop over V0 particles
  for(Int_t iv0=0; iv0<event->GetNumberOfV0s();iv0++){
    AliESDv0 *v0 = (AliESDv0 *) event->GetV0(iv0);
 
    if(!v0) continue;
    if(v0->GetOnFlyStatus()) continue; 
    // Get the particle selection 
    Bool_t foundV0 = kFALSE;
    Int_t pdgV0, pdgP, pdgN;
    foundV0 = fV0cutsESD->ProcessV0(v0, pdgV0, pdgP, pdgN);
    if(!foundV0) continue;
    Int_t iTrackP = v0->GetPindex();  // positive track
    Int_t iTrackN = v0->GetNindex();  // negative track

    // v0 Armenteros plot (QA)
    Float_t armVar[2] = {0.0,0.0};
    fV0cutsESD->Armenteros(v0, armVar);
    fhArmenteros->Fill(armVar[0],armVar[1]);
      
    if( pdgP == -11){
      fV0electrons->Add((AliVTrack*)event->GetTrack(iTrackP));
      // fV0tags[iTrackP] = 11;
    }
    else if( pdgP == 211){
      fV0pions->Add((AliVTrack*)event->GetTrack(iTrackP));
      // fV0tags[iTrackP] = 211;
    }
    else if( pdgP == 2212){
      fV0protons->Add((AliVTrack*)event->GetTrack(iTrackP));
      //fV0tags[iTrackP] = 2212;
    }

    // negative particles
    if( pdgN == 11){
      fV0electrons->Add((AliVTrack*)event->GetTrack(iTrackN));
      // fV0tags[iTrackN] = -11;
    }
    else if( pdgN == -211){
      fV0pions->Add((AliVTrack*)event->GetTrack(iTrackN));
      //        fV0tags[iTrackN] = -211;
    }
    else if( pdgN == -2212){
      fV0protons->Add((AliVTrack*)event->GetTrack(iTrackN));
      //        fV0tags[iTrackN] = -2212;
    }
  }
}



void AliAnalysisTaskHaHFECorrel::FindV0CandidatesAOD(AliAODEvent *event) 
{
  
  fV0cuts->SetMode(AliAODv0KineCuts::kEffGamma,AliAODv0KineCuts::kPP); // set mode (maxEfficiency or Purity kPurity);
  fV0cuts->SetEvent(event);

  const Int_t numTracks = event->GetNumberOfTracks();
  fV0tags = new Int_t[numTracks];
  for (Int_t i = 0; i < numTracks; i++) fV0tags[i] = 0;

  // loop over V0 particles
  for(Int_t iv0=0; iv0<event->GetNumberOfV0s();iv0++){

    AliAODv0 *v0 = (AliAODv0 *) event->GetV0(iv0);
      
    if(!v0) continue;
    if(v0->GetOnFlyStatus()) continue; 
    // Get the particle selection 
    Bool_t foundV0 = kFALSE;
    Int_t pdgV0, pdgP, pdgN;
    foundV0 = fV0cuts->ProcessV0(v0, pdgV0, pdgP, pdgN);
    if(!foundV0) continue;
    Int_t iTrackP = v0->GetPosID();  // positive track
    Int_t iTrackN = v0->GetNegID();  // negative track
      
    if (iTrackP<0 || iTrackN <0) {
      cout << "negative ID" << endl;
      continue;
    }
      
    fhArmenteros->Fill(v0->Alpha(),v0->QtProng());

    // AliAODTrack (V0 Daughters)
    AliAODVertex* vertex = v0->GetSecondaryVtx();
    if (!vertex) {
      Printf("ERROR: Could not retrieve vertex");
      continue;
    }
      
    AliAODTrack *pTrack=0, *nTrack=0;
    Int_t NDaughters = vertex->GetNDaughters(); // should be only two from construction
    for (Int_t i = 0; i<NDaughters; i++) {
      if (((AliAODTrack*)vertex->GetDaughter(i))->GetID()==iTrackP) pTrack= (AliAODTrack*)vertex->GetDaughter(i);
      else if  (((AliAODTrack*)vertex->GetDaughter(i))->GetID()==iTrackN) nTrack= (AliAODTrack*)vertex->GetDaughter(i);
    }
    if (pTrack==0 || nTrack ==0) continue;
    // if (iTrackP<0 || iTrackN<0) continue;
      
    // fill the Object arrays
    // positive particles
    if( pdgP == -11){
      fV0electrons->Add(pTrack);
      //     fV0tags[iTrackP] = 11;
    }
    else if( pdgP == 211){
      fV0pions->Add(pTrack);
      // fV0tags[iTrackP] = 211;
    }
    else if( pdgP == 2212){
      fV0protons->Add(pTrack);
      //	fV0tags[iTrackP] = 2212;
    }
      
    // negative particles
    if( pdgN == 11){
      fV0electrons->Add(nTrack);
      // fV0tags[iTrackN] = -11;
    }
    else if( pdgN == -211){
      fV0pions->Add(nTrack);
      // fV0tags[iTrackN] = -211;
    }
    else if( pdgN == -2212){
      fV0protons->Add(nTrack);
      //fV0tags[iTrackN] = -2212;
    }
  }
}


void AliAnalysisTaskHaHFECorrel::ClearV0PIDList(){

  // Clear the PID object arrays
  fV0electrons->Clear();
  fV0pions->Clear();
  fV0protons->Clear();
  if (fV0tags!=0) delete[] fV0tags;
  fV0tags = 0;
}

void AliAnalysisTaskHaHFECorrel::TRDQA(Int_t RunNumber, const AliVVertex *pVtx, Int_t nMother, Int_t listMother[]) {

  // Fill full Histograms
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();
 
  // Loop over all tracks 
  for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
          
    AliVParticle* VHtrack = 0x0;
    if(!fUseTender) VHtrack  = fVevent->GetTrack(jTracks);
    if(fUseTender) VHtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jTracks));
    if (!VHtrack) {
      printf("ERROR: Could not receive associated track %d\n", jTracks);
      continue;
    }
    AliVTrack *Vtrack = dynamic_cast<AliVTrack*>(VHtrack);
    if(!Vtrack) continue;
    AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);   
    if (fIsAOD && !AODtrack) continue;
    AliESDtrack *ESDtrack = dynamic_cast<AliESDtrack*>(Vtrack);   
    if (!fIsAOD && !ESDtrack) continue;


    if (fIsAOD) {
      if (!AODtrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
    }
    else {
      if(!ESDkTrkGlobalNoDCA(Vtrack)) continue;
    }

    Bool_t passHFETrackCut=kFALSE;
    passHFETrackCut= InclElecTrackCuts(pVtx,Vtrack,nMother,listMother);

    Bool_t passHFEPIDCut=kFALSE;
    if (passHFETrackCut) passHFEPIDCut= InclElecPIDCuts(Vtrack, kFALSE);
  
    Double_t fillSparse[7];
    fillSparse[0]=Vtrack->Pt();
    fillSparse[1]=Vtrack->Eta();
    fillSparse[2]=Vtrack->Phi();
    fillSparse[4]=Vtrack->Charge();
    fillSparse[5]=RunNumber;
    
    Int_t nTracklets=0;
    for (Int_t i=0; i<6; i++) {
      if (Vtrack->GetTRDmomentum(i)>0.01) { // TRDMomentumCriterion
	if (!fIsAOD && !ESDtrack->IsTRDtrkltChmbGood(i)) {
	  cout << "Chamber not good" << i << "\t" <<  Vtrack->Eta() << "\t" << Vtrack->Phi() << endl;
	  continue;
	}
	nTracklets++;
	fillSparse[3]=i+1; // Layer
	fTRDEtaPhi->Fill(fillSparse);
      }
    }
       
    if (passHFETrackCut)  fillSparse[2]=1; //PID - pass HFEcut
    else fillSparse[2]=2; // PID
    fillSparse[3]=fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
    fillSparse[4]=nTracklets;
    fillSparse[5]=Vtrack->Charge();
    fillSparse[6]=RunNumber;
    fTRDNTracklets->Fill(fillSparse);
    
    if (Vtrack->Eta()>0.8 || Vtrack->Eta()<-0.8) continue;
    fTRDnTrackRun->Fill(nTracklets,  (Int_t)Vtrack->GetTRDntrackletsPID(), RunNumber);

    fillSparse[1]=nTracklets;
    if (passHFETrackCut)  fillSparse[2]=1;
    else fillSparse[2]=2;
    fillSparse[3]=fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
    fillSparse[4]=Vtrack->Charge();
    fTRDSpectra->Fill(fillSparse);

    if (fIsMC) {
      Int_t MClabel=-999;
      Int_t PDGCode=-990; 
      if (fIsAOD) {
	if (!AODtrack) continue;
	MClabel=AODtrack->GetLabel();
	AliAODMCParticle* MCParticleAOD = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(abs(MClabel)));  
	if (!MCParticleAOD) continue;
	PDGCode = abs(MCParticleAOD->GetPdgCode());
      }
      else {
	if (!ESDtrack) continue;
	MClabel=ESDtrack->GetLabel();
	fMC->GetTrack(abs(MClabel));
	//cout << MClabel << endl;
	//cout << fMC->GetNumberOfTracks() << endl;
	AliMCParticle* MCParticleESD = dynamic_cast<AliMCParticle*>(fMC->GetTrack(abs(MClabel)));  
	if (!MCParticleESD) continue;
	PDGCode = abs(MCParticleESD->PdgCode());  
	//cout << MClabel << " TrackPt " << Vtrack->Pt() << " MCpt " << MCParticleESD->Pt() << endl;
      }

    
      if (PDGCode == 11) {
	fillSparse[2]=1;
	fTRDMCSpectra->Fill(fillSparse);
      }
      else if (PDGCode == 211) {
	fillSparse[2]=3;
	fTRDMCSpectra->Fill(fillSparse);
      }
      else if (PDGCode == 2212) {
	fillSparse[2]=2;
	fTRDMCSpectra->Fill(fillSparse);
      }
    }
  }



  for (Int_t i=0; i<fV0electrons->GetEntriesFast(); i++) {
     AliAODTrack *track = (AliAODTrack*) fV0electrons->At(i);
     FillV0Histograms(track, 1, RunNumber);
  }
  for (Int_t i=0; i<fV0protons->GetEntriesFast(); i++) {
     AliAODTrack *track = (AliAODTrack*) fV0protons->At(i);
     FillV0Histograms(track, 2, RunNumber);
  }
  for (Int_t i=0; i<fV0pions->GetEntriesFast(); i++) {
     AliAODTrack *track = (AliAODTrack*) fV0pions->At(i);
     FillV0Histograms(track, 3, RunNumber);
  }
}



void AliAnalysisTaskHaHFECorrel::FillV0Histograms(AliVTrack* Vtrack, Int_t Species, Int_t RunNumber) {
   
  AliESDtrack *ESDtrack = dynamic_cast<AliESDtrack*>(Vtrack);
  AliAODTrack *AODtrack = dynamic_cast<AliAODTrack*>(Vtrack);
  if (fIsAOD) {
    if (!AODtrack) return;
    if (!AODtrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return;
  }
  else {
    if (!ESDtrack) return;
    if(!ESDkTrkGlobalNoDCA(Vtrack)) return;
  }

  Double_t fillSparse[7];
  fillSparse[0]=Vtrack->Pt();
  fillSparse[1]=Vtrack->Eta();
  fillSparse[2]=Species;
  fillSparse[3]=fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
  fillSparse[5]=Vtrack->Charge();
  fillSparse[6]=RunNumber;
  
  Int_t nTracklets=0;
  for (Int_t i=0; i<7; i++) {
    if (Vtrack->GetTRDmomentum(i)>0.01) {
      if (!fIsAOD && !ESDtrack->IsTRDtrkltChmbGood(i)) {
	continue;
      }
      nTracklets++;
    }
  }
  
  fillSparse[4]=nTracklets;
  fTRDV0NTracklets->Fill(fillSparse);
  
  if (Vtrack->Eta()>0.8 || Vtrack->Eta()<-0.8) return;
  
  fillSparse[1]= nTracklets;
  fillSparse[2]= Species;
  fillSparse[3]= fpidResponse->NumberOfSigmasTPC(Vtrack, AliPID::kElectron);
  fillSparse[4]= Vtrack->Charge();
  fTRDV0Spectra->Fill(fillSparse);
    
}


void AliAnalysisTaskHaHFECorrel::CheckElectronIsTrigger(Double_t ptH, Bool_t *ElectronIsTrigger) {
  for (Int_t i=0; i<fAssPtHad_Nbins; i++) {
       if (ptH>fAssPtHad_Xmin[i] && ptH< fAssPtHad_Xmax[i]) ElectronIsTrigger[i]=kTRUE;
  }
  return;
}
  
void AliAnalysisTaskHaHFECorrel::CheckHadronIsTrigger(Double_t ptE, Bool_t* HadronIsTrigger) {
  for (Int_t i=0; i<fAssPtElec_Nbins; i++) {
    if (ptE>fAssPtElec_Xmin[i] && ptE<fAssPtElec_Xmax[i]) HadronIsTrigger[i]=kTRUE;
  }
  return;
}
 

Bool_t AliAnalysisTaskHaHFECorrel::PassEventBias( const AliVVertex *pVtx, Int_t nMother, Int_t *listMother) {
 
  
  // Leading Particle
  Double_t ptMax=-999;
  Double_t pt =-999;
  
  if (fIsMC) {
    for (Int_t i=0; i<fMC->GetNumberOfTracks(); i++) {
      AliVParticle *MCVParticle  = dynamic_cast<AliVParticle*>(fMC->GetTrack(i));
      AliAODMCParticle *MCAODParticle = dynamic_cast<AliAODMCParticle*>(fMC->GetTrack(i));
      AliMCParticle *MCESDParticle = dynamic_cast<AliMCParticle*>(fMC->GetTrack(i));
      if (fIsAOD)  {
	if (!MCAODParticle) continue;
	if (!MCAODParticle->IsPhysicalPrimary()) continue;
      }
      else {
	if (!fStack || !MCESDParticle) continue;
	if (!fStack->IsPhysicalPrimary(MCESDParticle->Label())) continue;
      }

      if(MCVParticle->Eta()>fMaxHadronEta) continue;
      if(MCVParticle->Eta()<fMinHadronEta) continue;
      if (pt>ptMax) ptMax=pt;
    }
  }
  else {
    Int_t ntracks = -999;
    if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
    if(fUseTender)  ntracks = fTracks_tender->GetEntries();

    // Loop over all tracks to find LP and HFE
    for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
          
      AliVParticle* VHtrack = 0x0;
      if(!fUseTender) VHtrack  = fVevent->GetTrack(jTracks);
      if(fUseTender) VHtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jTracks));
        
      if (!VHtrack) {
	printf("ERROR: Could not receive associated track %d\n", jTracks);
	continue;
      }
      AliVTrack *Vtrack = dynamic_cast<AliVTrack*>(VHtrack);
      if(!Vtrack) continue;
        
    
      pt = Vtrack->Pt();

      // track cuts
      Bool_t passHadTrackCut=kFALSE;
      passHadTrackCut = ChargedHadronTrackCuts(pVtx, Vtrack, nMother, listMother);
   
      Bool_t passHadPIDCut=kFALSE;
      passHadPIDCut = ChargedHadronPIDCuts(Vtrack); // currently empty
    
      // find hadron with the largest pT -> leading particle
      if (passHadTrackCut && passHadPIDCut) {
	if (pt>ptMax) ptMax=pt;
      }
    }
  }
  
    if (ptMax>fMinPtEvent && ptMax<fMaxPtEvent) return kTRUE;
    else return kFALSE;
}

Bool_t AliAnalysisTaskHaHFECorrel::ESDkTrkGlobalNoDCA(AliVTrack* Vtrack) {
  AliESDtrack *ESDtrack = dynamic_cast<AliESDtrack*>(Vtrack);   
  if (!fIsAOD && !ESDtrack) return kFALSE;
  fesdTrackCuts->SetMaxDCAToVertexXY(2.4);
  fesdTrackCuts->SetMaxDCAToVertexZ(3.2);
  fesdTrackCuts->SetDCAToVertex2D(kTRUE);
  if (fesdTrackCuts->AcceptTrack(ESDtrack)) return kTRUE;
  else return kFALSE;
}
