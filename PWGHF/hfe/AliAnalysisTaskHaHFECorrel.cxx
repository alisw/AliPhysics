
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
,fCorrHadron(kTRUE)
,fCorrLParticle(kTRUE)
,fMixedEvent(kTRUE)
,fLParticle(kTRUE)
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
,fIsAOD(kFALSE)
,fCFM(0)
,fPID(0)
,fPIDqa(0)
,fOutputList(0)
,fNoEvents(0)
,fTrkpt(0)
,fEtaVtxZ(0)
,fMultiplicity(0)
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
,fMCElecHaTruePartnerTrigger(0)
,fMCElecHaNoPartnerTrigger(0)
,fElecHaMixedEvent(0)
,fLSElecHaMixedEvent(0)
,fULSElecHaMixedEvent(0)
,fElecLPTrigger(0)
,fLPElecTrigger(0)
,fLPNonElecTrigger(0)
,fNonElecLPTrigger(0)
,fInclElecLP(0)
,fLSElecLP(0)
,fULSElecLP(0)
,fElecLPHa(0)
,fElecLPLSNoPartner(0)
,fElecLPULSNoPartner(0)
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
,fCorrHadron(kTRUE)
,fCorrLParticle(kTRUE)
,fMixedEvent(kTRUE)
,fLParticle(kTRUE)
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
,fIsAOD(kFALSE)
,fCFM(0)
,fPID(0)
,fPIDqa(0)
,fOutputList(0)
,fNoEvents(0)
,fTrkpt(0)
,fEtaVtxZ(0)
,fMultiplicity(0)
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
,fMCElecHaTruePartnerTrigger(0)
,fMCElecHaNoPartnerTrigger(0)
,fElecHaMixedEvent(0)
,fLSElecHaMixedEvent(0)
,fULSElecHaMixedEvent(0)
,fElecLPTrigger(0)
,fLPElecTrigger(0)
,fLPNonElecTrigger(0)
,fNonElecLPTrigger(0)
,fInclElecLP(0)
,fLSElecLP(0)
,fULSElecLP(0)
,fElecLPHa(0)
,fElecLPLSNoPartner(0)
,fElecLPULSNoPartner(0)
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
  
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD){
    printf("ERROR: fAOD not available\n");
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

  // suggested by DPG to remove outliers
  const AliAODVertex *pVtx =   fAOD->GetPrimaryVertex();
  const AliAODVertex *spdVtx = fAOD->GetPrimaryVertexSPD();
    
  Double_t pVtxZ = -999.;
  pVtxZ = pVtx->GetZ();
  if (TMath::Abs(pVtxZ)>10.) return;
      
  if(fIsMC) {
    fMC = MCEvent();
    fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!fMCheader) AliError("fMCheader could not be initialised");
    MCEfficiencyCorrections(pVtx); //  Electron reconstruction, Hadron reconstruction


    TList *lh=fMCheader->GetCocktailHeaders();
    Int_t nh=lh->GetEntries();  
    for(Int_t i=0;i<nh;i++)
    {
        AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(i);
        TString genname=gh->GetName();
        Int_t npart=gh->NProduced();
	//	cout << genname << "with " << npart << " Particles" << endl;
   }
  }


   
  fNoEvents->Fill(0);
  
  Int_t fNOtrks = fAOD->GetNumberOfTracks();
  if(fNOtrks<2) return;
  fNoEvents->Fill(1);
    
  if(!fEventCuts.AcceptEvent(fAOD)) {
    PostData(1, fOutputList);
    return;
  }
  fNoEvents->Fill(2);

  fpidResponse = fInputHandler->GetPIDResponse();
  if(!fpidResponse){
    AliDebug(1, "Using default PID Response");
    fpidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class());
  }  
      
  Int_t nVertices = 1;
  nVertices = fAOD->GetNumberOfVertices();
  Double_t listofmotherkink[nVertices];
  Int_t nMotherKink = 0;
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

  fNoEvents->Fill(3);


  Double_t fMultV0Per, fMultSPDPer, fMultV0Tot, fMultSPD;
  Float_t mult = 1;
  if (!fIsMC) {
    //multiplicity -- commented because problem with LHC13
    fMultV0Per = -1.;
    fMultSPDPer =-1.;
    // fAOD->GetMultiplicity();
    fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
    fMultV0Per = fMultSelection->GetMultiplicityPercentile("V0M", kFALSE); // Method, Embed Event selection (kFALSE)
    fMultSPDPer = fMultSelection->GetMultiplicityPercentile("SPDTracklets", kFALSE); // Method, Embed Event selection (kFALSE)
    
    // TRY Multiplicity estimates
    AliVVZERO* AODV0 = fVevent->GetVZEROData();
    Float_t multV0A=AODV0->GetMTotV0A();
    Float_t multV0C=AODV0->GetMTotV0C();
    fMultV0Tot=multV0A+multV0C;
    
    AliVMultiplicity* AODMult = fVevent->GetMultiplicity();
    fMultSPD=AODMult->GetNumberOfTracklets();
    
    AliAODTracklets *AODtracklets = ((AliAODEvent*)fAOD)->GetTracklets();
    //  fMultAOD = AODtracklets->GetNumberOfTracklets();

    mult=fMultV0Per;

    Double_t fillSparse[4]={fMultV0Per, fMultSPDPer, fMultV0Tot, fMultSPD};
    fMultiplicity->Fill(fillSparse);
  }

  //cout << "RunNumber " << fVevent->GetRunNumber() << endl;



  ///////////////////////
  // Preparational Tasks
  ///////////////////////

  // List of HFE for analysis and mixed event
  TObjArray* RedTracksHFE = new TObjArray;
  RedTracksHFE->SetOwner(kTRUE);

  // FirstLoop: Find leading particle (LP) and HFE (check for LS/ULS partner)
  Int_t LeadingParticle;
  AliAODTrack* LPtrack;
  fLParticle=kFALSE;
  LPtrack=FindLPAndHFE(RedTracksHFE, pVtx,nMotherKink,listofmotherkink);

  FindV0Candidates(fAOD);
  fEventsPerRun->Fill(fVevent->GetRunNumber());
  TRDQA(fVevent->GetRunNumber(), pVtx,nMotherKink,listofmotherkink);

   //////////////////////////
  // Main analysis
  ///////////////////////////

  // LP - two different functions for Same Event and mixed Event
  if (fLParticle && RedTracksHFE->GetEntriesFast()>0) CorrelateLP(LPtrack,pVtx, nMotherKink, listofmotherkink, RedTracksHFE);
  if (fLParticle && fMixedEvent  && RedTracksHFE->GetEntriesFast()>0) CorrelateLPMixedEvent(LPtrack, mult, pVtx->GetZ());
 

  // Hadron - only one function for both options, as only one loop over Hadron tracks
  // Mixed event is called inside this function
  if (RedTracksHFE->GetEntriesFast()>0) CorrelateHadron(RedTracksHFE, pVtx, nMotherKink, listofmotherkink, mult);

  // Electrons - currently no mixed event (>1 as two electrons are required), work in progress
  if (RedTracksHFE->GetEntriesFast()>1) CorrelateElectron(RedTracksHFE);
 

  /////////////////////////
  //Fill Mixed event pool//
  /////////////////////////
  
  AliEventPool * HFEPool = fPoolMgr->GetEventPool(mult, pVtxZ); // Get the buffer associated with the current centrality and z-vtx
  //  HFEPool->SetDebug(kTRUE);
  if (!HFEPool)
  {
    AliFatal(Form("No pool found for centrality = %f, zVtx = %f", mult, pVtx->GetZ()));
    return;
  }
  if (RedTracksHFE->GetEntriesFast()>0) {
    HFEPool->UpdatePool(RedTracksHFE);
  }

  ClearV0PIDList();

  PostData(1, fOutputList);
}

//_________________________________________
void AliAnalysisTaskHaHFECorrel::UserCreateOutputObjects()
{
  Bool_t UseSparse = kTRUE;
  Bool_t UseTree = kFALSE;		       
  
  AliLog::SetClassDebugLevel("AliAnalysisTaskHaHFECorrel", 0);
  AliDebug(1, Form("This is AliDebug level %i", 1));
  AliDebug(2, Form("This is AliDebug level %i", 2));
  AliDebug(3, Form("This is AliDebug level %i", 3));
   
  //---------Output Tlist
  fOutputList = new TList();
  fOutputList->SetOwner();

  fEventCuts.AddQAplotsToList(fOutputList);

  // V0 Kine cuts 
  fV0cuts = new AliAODv0KineCuts();

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
  // HFE Electron
  Int_t    NBinsElectron =46;
  Double_t XminElectron=0.25;
  Double_t XmaxElectron=6.;
  Int_t    NBinsElectronRed = 7;
  Double_t XBinsElectronRed[8]={0.25,0.5, 1., 1.5, 3, 4, 6, 10};


  Int_t     NBinsHadron=200 ;
  Double_t  XminHadron=0.0;
  Double_t  XmaxHadron=100;
  Int_t     NBinsHadRed=13;
  Double_t  XBinsHaRed[14]={0.5,1.0, 1.5, 2.0, 2.5, 3.,3.5, 5., 7.5, 15., 30., 50., 100, 1000}; //13 bins

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
  
  if (fIsMC) {
    Int_t     binHMC[4] = {NBinsElectron  , 7,  20, 200}; // p PDG ITS, TPC
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
 

  if (fIsMC) {
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


  fElecHadTrigger = new TH1F("fElecHadTrigger", "fElecHadTrigger", bin[1], xmin[1], xmax[1]);
  fElecHadTrigger->Sumw2();
  fOutputList->Add(fElecHadTrigger);						      


  fElecHadTriggerLS = new TH1F("fElecHadTriggerLS", "fElecHadTriggerLS", bin[1], xmin[1], xmax[1]);
  fElecHadTriggerLS->Sumw2();
  fOutputList->Add(fElecHadTriggerLS);						      

  fElecHadTriggerULS = new TH1F("fElecHadTriggerULS", "fElecHadTriggerULS", bin[1], xmin[1], xmax[1]);
  fElecHadTriggerULS->Sumw2();
  fOutputList->Add(fElecHadTriggerULS);						      

  fHadContTrigger = new TH1F("fHadContTrigger", "fHadContTrigger", bin[1], xmin[1], xmax[1]);
  fHadContTrigger->Sumw2();
  fOutputList->Add(fHadContTrigger);						      


  fNonElecHadTrigger = new TH1F("fNonElecHadTrigger", "fNonElecHadTrigger", bin[1], xmin[1], xmax[1]);
  fNonElecHadTrigger->Sumw2();
  fOutputList->Add(fNonElecHadTrigger);						      

  fMCElecHaTruePartnerTrigger = new TH1F("fMCElecHaTruePartnerTrigger", "fMCElecHaTruePartnerTrigger", bin[1], xmin[1], xmax[1]);
  fMCElecHaTruePartnerTrigger->Sumw2();
  fOutputList->Add(fMCElecHaTruePartnerTrigger);	

  fMCElecHaNoPartnerTrigger = new TH1F("fMCElecHaNoPartnerTrigger", "fMCElecHaNoPartnerTrigger", bin[1], xmin[1], xmax[1]);
  fMCElecHaNoPartnerTrigger->Sumw2();
  fOutputList->Add(fMCElecHaNoPartnerTrigger);

  fHadElecTrigger = new TH1F("fHadElecTrigger", "fHadElecTrigger", bin[0], xmin[0], xmax[0]);
  fHadElecTrigger->Sumw2();
  fOutputList->Add(fHadElecTrigger);						      

  fHadNonElecTrigger = new TH1F("fHadNonElecTrigger", "fHadNonElecTrigger", bin[0], xmin[0], xmax[0]);
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


  // LP HFE

  fElecLPTrigger = new TH1F("fElecLPTrigger", "fElecLPTrigger", bin[1], xmin[1], xmax[1]);
  fElecLPTrigger->Sumw2();
  fOutputList->Add(fElecLPTrigger);						      

  fNonElecLPTrigger = new TH1F("fNonElecLPTrigger", "fNonElecLPTrigger", bin[1], xmin[1], xmax[1]);
  fNonElecLPTrigger->Sumw2();
  fOutputList->Add(fNonElecLPTrigger);						      

  fLPElecTrigger = new TH1F("fLPElecTrigger", "fLPElecTrigger", bin[0], xmin[0], xmax[0]);
  fLPElecTrigger->Sumw2();
  fOutputList->Add(fLPElecTrigger);						      

  fLPNonElecTrigger = new TH1F("fLPNonElecTrigger", "fLPNonElecTrigger", bin[0], xmin[0], xmax[0]);
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
    
  fElecLPLSNoPartner = new THnSparseF("fEleLPLSNoPartner", "Sparse for LSEle-Had with no Partner: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecLPLSNoPartner->Sumw2();
  fOutputList->Add(fElecLPLSNoPartner);
    
  fElecLPULSNoPartner = new THnSparseF("fEleLPULSNoPartner", "Sparse for ULSEle-Had with no Partner: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecLPULSNoPartner->Sumw2();
  fOutputList->Add(fElecLPULSNoPartner);
    
  fElecLPHa = new THnSparseF("fEleLPHa", "Sparse for Hadron Contamination: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecLPHa->Sumw2();
  fOutputList->Add(fElecLPHa);
    
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

  fhArmenteros  = new TH2F("fhArmenteros","Armenteros plot",200,-1.,1.,200,0.,0.4);
  fhArmenteros->Sumw2();
  fOutputList->Add(fhArmenteros);

  fEventsPerRun = new TH1F("fEventsPerRun", "EventsPerRun",300000, 100000,400000);
  fEventsPerRun->Sumw2();
  fOutputList->Add(fEventsPerRun);

  fTRDnTrackRun = new TH3F("fTRDnTrackRun", "TrackletCharge, TrackletPID, Run", 6, 0.5, 6.5, 6, 0.5, 6.5, 300000, 100000, 400000);
  fTRDnTrackRun->Sumw2();
  fOutputList->Add(fTRDnTrackRun);




  
  Int_t TRDBins[6]={NBinsElectron, NBinsEta, NBinsPhi*2, 7, 3, 300000};
  Double_t TRDBMin[6]={XminElectron, XminEta, 0, -0.5, -1.5, 100000.5};
  Double_t TRDBMax[6]={XmaxElectron, XmaxEta, TMath::TwoPi(), 6.5, 1.5, 400000.5};


  fTRDEtaPhi = new THnSparseF("fTRDEtaPhi", "fTRDEtaPhi: Pt, Eta, Phi, Layer, Charge, run", 6, TRDBins, TRDBMin, TRDBMax);
  fTRDEtaPhi->Sumw2();
  fOutputList->Add(fTRDEtaPhi);

  Int_t TRDBins3[7]={NBinsElectron,   36       , 3,   14,     7, 3, 300000};
  Double_t TRD3BMin[7]={XminElectron, XminEta,   0.5, -3,  -0.5, -1.5, 100000.5};
  Double_t TRD3BMax[7]={XmaxElectron, XmaxEta,    3.5, 4,   6.5, 1.5 , 400000.5};

  fTRDV0NTracklets = new THnSparseF("fTRDV0NTracklets", "fTRDnTracklets: Pt, Eta, PID, TPC, NTracklets, Charge, run", 7, TRDBins3, TRD3BMin, TRD3BMax);
  fTRDV0NTracklets->Sumw2();
  fOutputList->Add(fTRDV0NTracklets);
  
  fTRDNTracklets = new THnSparseF("fTRDNTracklets", "fTRDnTracklets: Pt, Eta, PID, TPC, NTracklets, Charge, run",  7, TRDBins3, TRD3BMin, TRD3BMax);
  fTRDNTracklets->Sumw2();
  fOutputList->Add(fTRDNTracklets);


  Int_t TRDBins2[5]={NBinsElectron, 7, 3, 14, 3};
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
  

  Int_t    poolSize = 1000;
  Int_t    trackDepth = 2000; 
  Int_t    nMultBins = 4;
  Double_t multBins[]={0,20,50,100};
  Int_t    nZVtxBins=9;
  Double_t zVtxBins[]={-10,-6.91,-4.41,-2.41,-0.41,1.59,3.59,5.59,8.09,10};
  fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, nMultBins, (Double_t*) multBins, nZVtxBins, (Double_t*) zVtxBins);
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
AliAODTrack*  AliAnalysisTaskHaHFECorrel::FindLPAndHFE( TObjArray* RedTracks, const AliAODVertex *pVtx, Int_t nMother, Double_t listMother[])
{
  AliAODTrack* LPtrack=0;
  fLParticle=kFALSE;
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();
  
  // Leading Particle
  Double_t pH=-9.,ptH=-9.,etaH =-9.,phiH=-9.;
  Int_t chargeH = 0, pdgH = 0;

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
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(VHtrack);
    if(!track) continue;
        
    Double_t p=-9.,pt=-9.,eta =-9.,phi=-9.;
    Int_t id=-9;  
    pt = track->Pt();
    p = track->P();
    phi = track->Phi();
    eta = track->Eta();
    id = track->GetID();

   
    fTrkpt->Fill(pt,0);

    // track cuts
    Bool_t passHadTrackCut=kFALSE;
    passHadTrackCut = ChargedHadronTrackCuts(pVtx, track, nMother, listMother);
   
    Bool_t passHadPIDCut=kFALSE;
    passHadPIDCut = ChargedHadronPIDCuts(track); // currently empty
    
    NHadrons++;

    // find hadron with the largest pT -> leading particle
    if (passHadTrackCut && passHadPIDCut) {

      // for LP
      if(track->Pt()>ptH) {
	ptH = track->Pt();
	pH = track->P();
	phiH =track->Phi();
	etaH = track->Eta();
	chargeH = track->Charge();
	LPtrack=track;
	fLParticle=kTRUE;
      }

      fEtaVtxZ->Fill(eta, pVtx->GetZ());

      Double_t fillSparse[3];
      fillSparse[0]=pt;
      fillSparse[1]=eta;
      fillSparse[2]=phi;
      fillSparse[3]=pVtx->GetZ();
      fRecHadPtEtaPhiVtx->Fill(fillSparse);

      if (fIsMC) {
	Int_t MClabel=track->GetLabel();
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
    passHFETrackCut= InclElecTrackCuts(pVtx,track,nMother,listMother);
    if (passHFETrackCut)  fTrkpt->Fill(pt,1); // after track cuts
 

    Bool_t passHFEPIDCut=kFALSE;
    if (passHFETrackCut) passHFEPIDCut=  InclElecPIDCuts(track, kTRUE);
    if (passHFETrackCut && passHFEPIDCut)  fTrkpt->Fill(pt,2);
    
    Int_t lsPartner=0, ulsPartner=0;
    Int_t lsPartnerID[20], ulsPartnerID[20];
    Bool_t trueULSPartner = kFALSE;
    Bool_t isPhotonic = kFALSE;
    Bool_t isHadron = kFALSE;
    if (passHFETrackCut && passHFEPIDCut) { // if HFE is found, look for ls and uls partner
      NElectrons++;
      FindPhotonicPartner(jTracks, track, pVtx, nMother, listMother, lsPartner, ulsPartner, lsPartnerID, ulsPartnerID, trueULSPartner, isPhotonic);
      if (fIsMC) {
	if (isPhotonic) {
	  if (trueULSPartner) NPhotElectronsTagged++;
	  else  NPhotElectronsUntagged++;
	}
	EvaluateTaggingEfficiency(track, lsPartner, ulsPartner, trueULSPartner);
      }
      fInclElecPt->Fill(pt);
      fInclElecP->Fill(p);
      fInclElecPhi->Fill(pt,phi); // phi of electron candidates
      fInclElecEta->Fill(pt,eta); 


      if (fIsMC) {
	Int_t MClabel=track->GetLabel();
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
	    fillSparse[3]==pVtx->GetZ();
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
	}
      }
   

      
      // store only essential informations of these electrons for later correlations and mixed event
      CloneAndReduceTrackList(RedTracks, track, lsPartner, ulsPartner, lsPartnerID, ulsPartnerID, trueULSPartner, isPhotonic, isHadron);
    }


    Bool_t passNonElecPIDCut=kFALSE;
    passNonElecPIDCut=AssoHadronPIDCuts(track);
    if (passHFETrackCut && passNonElecPIDCut) NNonElectrons++;

     					   
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
  else if (fIsMC && RedTracks->GetEntriesFast()>0) {
    Int_t MClabel=LPtrack->GetLabel();
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
    fillSparse[0]=LPtrack->Pt();
    fillSparse[1]=PDGCode;
    fillSparse[2]=PDGCodeMother;
    fMCLeadingParticle->Fill(fillSparse);
  }



  return LPtrack;
}

//_________________________________________
void AliAnalysisTaskHaHFECorrel::FindPhotonicPartner(Int_t iTracks, AliAODTrack* track, const AliAODVertex *pVtx, Int_t nMother, Double_t listMother[], Int_t &lsPartner, Int_t &ulsPartner, Int_t *lsPartnerID, Int_t *ulsPartnerID, Bool_t &trueULSPartner, Bool_t &isPhotonic) {
  
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
    AliAODTrack *trackAsso = dynamic_cast<AliAODTrack*>(Vassotrack);
    if(!trackAsso) continue;
        
    // tagged particle
    Double_t p=-9.,pt=-9.,eta =-9.,phi=-9.;
    Int_t charge = 0, pdg = 0;
        
    pt = track->Pt();
    p = track->P();
    phi = track->Phi();
    eta = track->Eta();
    charge = track->Charge();
        
    // associated particle variables
    Double_t pAsso=-9.,ptAsso=-9.,etaAsso =-9.,phiAsso=-9.;
    Int_t chargeAsso = 0, pdgAsso = 0;
        
    ptAsso = trackAsso->Pt();
    pAsso = trackAsso->P();
    phiAsso = trackAsso->Phi();
    etaAsso = trackAsso->Eta();
    chargeAsso = trackAsso->Charge();

    double dphi = -9;

    // Cuts for associated electrons for HFE-HFE correlations
    Bool_t passAssoTrackCutIncl = kFALSE;
    passAssoTrackCutIncl = InclElecTrackCuts(pVtx,trackAsso,nMother,listMother);

    Bool_t passAssoPIDCutIncl = kFALSE;
    passAssoPIDCutIncl = InclElecPIDCuts(trackAsso, kFALSE);

    if (passAssoTrackCutIncl && passAssoPIDCutIncl) {
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      fElecDphi->Fill(pt,dphi);        // HFEHFE
    }

    // looser Track cuts for associated photonic eletron
    Bool_t passAssoTrackCutPhot = kFALSE;
    passAssoTrackCutPhot = PhotElecTrackCuts(pVtx,trackAsso,nMother,listMother);
    if(!passAssoTrackCutPhot) continue;
   
    Bool_t passAssoPIDCutPhot = kFALSE;
    passAssoPIDCutPhot = PhotElecPIDCuts(trackAsso); 
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
      AliKFParticle ge1(*track, fPDGe1);
      AliKFParticle ge2(*trackAsso, fPDGe2);
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
      extTrackParam1.CopyFromVTrack(track);
      AliExternalTrackParam extTrackParam2;
      extTrackParam2.CopyFromVTrack(trackAsso);

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
        
    if(mass<fInvmassCut && fFlagULS)
      {

	if (fIsMC && !trueULSPartner) {
	  trueULSPartner = HaveSameMother(track->GetLabel(), trackAsso->GetLabel());
	}
	fULSElecPhi->Fill(pt,phi);
	fULSElecPt->Fill(pt);
	dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
	fULSElecDphi->Fill(pt,dphi);
	ulsPartnerID[ulsPartner]=trackAsso->GetID();
	ulsPartner++;
      }
        
    if(mass<fInvmassCut && fFlagLS){
      fLSElecPhi->Fill(pt,phi);
      fLSElecPt->Fill(pt);
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      fLSElecDphi->Fill(pt,dphi);
      lsPartnerID[lsPartner]=trackAsso->GetID();
      lsPartner++;
    }

    if (fIsMC) { 
      isPhotonic=IsPhotonicElectron(track->GetLabel());
    }
      
  }//track loop
}

void AliAnalysisTaskHaHFECorrel::CorrelateElectron(TObjArray* RedTracksHFE) {

  // work in progress - only correlate electrons which fullfill inclusive cuts

  for (Int_t j=0; j<RedTracksHFE->GetEntriesFast(); j++) {

    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(j);

    // tagged particle
    Double_t p=-9.,pt=-9.,eta =-9.,phi=-9.;
    Int_t charge = 0, pdg = 0, ls=0, uls=0;    
    pt = RedTrack->Pt();
    phi = RedTrack->Phi();
    eta = RedTrack->Eta();
    charge = RedTrack->Charge();

    // fill trigger electron
    fElecTrigger->Fill(pt);

    for (Int_t k=0; k<RedTracksHFE->GetEntriesFast(); k++) {
      if (k==j) continue;
      AliBasicParticleHaHFE *trackAsso = (AliBasicParticleHaHFE*) RedTracksHFE->At(k);
      Double_t pAsso=-9.,ptAsso=-9.,etaAsso =-9.,phiAsso=-9.;
      Int_t chargeAsso = 0, lsAsso=0, ulsAsso=0;
        
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

void AliAnalysisTaskHaHFECorrel::CorrelateHadron(TObjArray* RedTracksHFE,  const AliAODVertex* pVtx, Int_t nMother, Double_t listMother[], Float_t mult) {

  Int_t HadronTrigger=kFALSE;
    
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();


  Bool_t ElectronIsTrigger[RedTracksHFE->GetEntriesFast()];
  Bool_t HadContIsTrigger[RedTracksHFE->GetEntriesFast()];
  Double_t ElectronIsTriggerPt[RedTracksHFE->GetEntriesFast()];
  for (Int_t l=0; l<RedTracksHFE->GetEntriesFast(); l++) {
    ElectronIsTrigger[l]=kFALSE;
    HadContIsTrigger[l]=kFALSE;
  }


  Bool_t PhotElecWPartnerTrigger[RedTracksHFE->GetEntriesFast()];
  Bool_t PhotElecWoPartnerTrigger[RedTracksHFE->GetEntriesFast()];
   for (Int_t l=0; l<RedTracksHFE->GetEntriesFast(); l++) {
    PhotElecWPartnerTrigger[l]=kFALSE;
    PhotElecWoPartnerTrigger[l]=kFALSE;
  }



  Bool_t NonElectronIsTrigger[ntracks];
  Double_t NonElectronIsTriggerPt[ntracks];
  for (Int_t l=0; l<ntracks; l++) {
    NonElectronIsTrigger[l]=kFALSE;
  }




  // Track loop for hadron correlations
  for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {     

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
    Double_t pH=-9.,ptH=-9.,etaH =-9.,phiH=-9.;
    Int_t chargeH = 0, pdgH = 0, idH=-9;
    ptH = track->Pt();
    pH = track->P();
    phiH = track->Phi();
    etaH = track->Eta();
    chargeH = track->Charge();
    idH = track->GetID();
    HadronTrigger=kFALSE;

    CorrelateWithHadrons(track, pVtx, nMother, listMother, kTRUE, kFALSE, NonElectronIsTrigger, NonElectronIsTriggerPt, RedTracksHFE->GetEntriesFast()); // Correlate hadrons with Hadron (kTRUE, kFALSE);

    if (fMixedEvent) CorrelateHadronMixedEvent(track, mult, pVtx->GetZ());

    // loop over all electrons


    for (Int_t k=0; k<RedTracksHFE->GetEntriesFast(); k++) {

      AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(k);

      if (RedTrack->ID()==idH) continue;
      // tagged particle
      Double_t p=-9.,pt=-9.,eta =-9.,phi=-9.;
      Int_t charge = 0, pdg = 0, ls=0, uls=0;    
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


    
      ElectronIsTrigger[k]=kTRUE;
      ElectronIsTriggerPt[k]=pt;
      HadronTrigger=kTRUE;

      fInclElecHa->Fill(fillSparse);
      for (Int_t j=0; j<ls; j++)  fLSElecHa->Fill(fillSparse);
      for (Int_t j=0; j<uls; j++) fULSElecHa->Fill(fillSparse);
      if (fIsMC && RedTrack->IsHadron()) {
	fMCElecHaHadron->Fill(fillSparse);
	HadContIsTrigger[k]=kTRUE;
      }


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
      if (fIsMC && RedTrack->IsPhotonic()) {
	if (RedTrack->TruePartner()) {
	  PhotElecWPartnerTrigger[k]=kTRUE;
	  fMCElecHaTruePartner->Fill(fillSparse);
	}
	else {
	  fMCElecHaNoPartner->Fill(fillSparse);
	  PhotElecWoPartnerTrigger[k]=kTRUE;
	}
      }
    }

    // Fill Trigger Hadrons
    if (HadronTrigger) fHadElecTrigger->Fill(ptH);

  }//end of track loop
  for (int l=0; l<RedTracksHFE->GetEntriesFast(); l++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(l);
    if (ElectronIsTrigger[l]) {
      fElecHadTrigger->Fill(ElectronIsTriggerPt[l]);
      if( RedTrack->LS()>0) fElecHadTriggerLS->Fill(ElectronIsTriggerPt[l],  RedTrack->LS());
      if( RedTrack->ULS()>0) fElecHadTriggerULS->Fill(ElectronIsTriggerPt[l],  RedTrack->ULS());
    }

    if (HadContIsTrigger[l]) fHadContTrigger->Fill(ElectronIsTriggerPt[l]);

    if (PhotElecWPartnerTrigger[l]) fMCElecHaTruePartnerTrigger->Fill(ElectronIsTriggerPt[l]);
    if (PhotElecWoPartnerTrigger[l]) fMCElecHaNoPartnerTrigger->Fill(ElectronIsTriggerPt[l]);
    // else cout << "ElectronIsNotTrigger by nevents " << ntracks << endl;
  }
  for (int l=0; l<ntracks; l++) {
    if (NonElectronIsTrigger[l]) fNonElecHadTrigger->Fill(NonElectronIsTriggerPt[l]);
  }
}


void AliAnalysisTaskHaHFECorrel::CorrelateHadronMixedEvent(AliAODTrack* Htrack, Float_t mult, Float_t zVtx) {


  Double_t pH=-9.,ptH=-9.,etaH =-9.,phiH=-9.;
  Int_t chargeH = 0, pdgH = 0, idH=-9;
        
  ptH = Htrack->Pt();
  pH = Htrack->P();
  phiH = Htrack->Phi();
  etaH = Htrack->Eta();
  chargeH = Htrack->Charge();
  idH = Htrack->GetID();

  AliEventPool* fPool;
  fPool = fPoolMgr->GetEventPool(mult, zVtx); // Get the buffer associated with the current centrality and z-vtx
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
      }
    }
  }
}


/////////////--------------

void AliAnalysisTaskHaHFECorrel::CorrelateWithHadrons(AliAODTrack* track, const AliAODVertex* pVtx, Int_t nMother, Double_t listMother[], Bool_t FillHadron, Bool_t FillLP, Bool_t* NonElecIsTrigger, Double_t *NonElecIsTriggerPt, Int_t NumElectronsInEvent) {
  // Trigger Hadron
  Double_t pH=-9.,ptH=-9.,etaH =-9.,phiH=-9.;
  Int_t chargeH = 0, pdgH = 0, idH=-9;
  ptH = track->Pt();
  pH = track->P();
  phiH = track->Phi();
  etaH = track->Eta();
  chargeH = track->Charge();
  idH = track->GetID();
  
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();

  Bool_t Trigger=kFALSE;           // flag to store information, if used as trigger particle
  //  Int_t  CountCorrelations=0;      // counter that makes sure, that all events are given the same weight as for h-e case
 

  // Track loop for hadron correlations
  for(Int_t iTracks = 0; iTracks < ntracks; iTracks++) {     
    

    //if (CountCorrelations==3*NumElectronsInEvent) continue;

    AliVParticle* Vtrack = 0x0;
    if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
    if(fUseTender)  Vtrack  = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));
    if (!Vtrack) {
      printf("ERROR: Could not receive tagged  track %d\n", iTracks);
      continue;
    }
    AliAODTrack *trackAssoHad = dynamic_cast<AliAODTrack*>(Vtrack);
    if(!trackAssoHad) continue;
        
    if(trackAssoHad->GetID()==idH) continue;

    // track cuts
    Bool_t passTrackCut=kFALSE;
    passTrackCut = InclElecTrackCuts(pVtx, trackAssoHad, nMother, listMother);
    if (!passTrackCut) continue;

    // pid cuts
    Bool_t passPIDCut=kFALSE;
    passPIDCut = AssoHadronPIDCuts(trackAssoHad);
    if (!passPIDCut) continue;
   
   

    // associated particle
    Double_t pAsso=-9.,ptAsso=-9.,etaAsso=-9.,phiAsso=-9.;
    Int_t chargeAsso = 0, pdgAsso = 0, idAsso=-999;
    ptAsso = trackAssoHad->Pt();
    pAsso = trackAssoHad->P();
    phiAsso = trackAssoHad->Phi();
    etaAsso = trackAssoHad->Eta();
    chargeAsso = trackAssoHad->Charge();
    idAsso = trackAssoHad->GetID();

    // if (FillLP) cout << CountCorrelations << " 3NELP " << 3*NumElectronsInEvent << endl;
    //else if (FillHadron) cout << CountCorrelations << " 3NEHadron " << 3*NumElectronsInEvent << endl;
    //CountCorrelations++;
    Trigger = kTRUE;
    if(FillLP) fNonElecLPTrigger->Fill(ptAsso);
    if(FillHadron){
      NonElecIsTrigger[iTracks]=kTRUE;
      NonElecIsTriggerPt[iTracks]=ptAsso;
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
  //  if (Trigger==kFALSE) cout << "NoPartner found bei tracks " << ntracks <<  endl;
  if (FillHadron && Trigger) fHadNonElecTrigger->Fill(ptH);
  if (FillLP && Trigger)     fLPNonElecTrigger->Fill(ptH);
}



//_______________________

void AliAnalysisTaskHaHFECorrel::CorrelateLP(AliAODTrack* LPtrack,  const AliAODVertex* pVtx, Int_t nMother, Double_t listMother[], TObjArray *RedTracksHFE) 
{ 
  Bool_t LPTrigger=kFALSE;

  Bool_t   NonElectronIsTrigger[1];
  Double_t NonElectronIsTriggerPt[1];


  // leading Particle neq 
  Double_t pH=-9.,ptH=-9.,etaH =-9.,phiH=-9.;
  Int_t chargeH = 0, pdgH = 0, idH=-9;
  ptH = LPtrack->Pt();
  pH = LPtrack->P();
  phiH = LPtrack->Phi();
  etaH = LPtrack->Eta();
  chargeH = LPtrack->Charge();
  idH = LPtrack->GetID();

  CorrelateWithHadrons(LPtrack, pVtx, nMother, listMother, kFALSE, kTRUE, NonElectronIsTrigger, NonElectronIsTriggerPt, RedTracksHFE->GetEntriesFast()); // correlate LPHadron (kFALSE, kTRUE);
		       
  // Only loop over electrons in event
  for (Int_t k=0; k<RedTracksHFE->GetEntriesFast(); k++) {
    AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(k);
  
    if (RedTrack->ID()==idH) continue; // no self-correlation
    Double_t p=-9.,pt=-9.,eta =-9.,phi=-9.;
    Int_t charge = 0, pdg = 0, ls=0, uls=0;    
    pt = RedTrack->Pt();
    // p = RedTrack->P();
    phi = RedTrack->Phi();
    eta = RedTrack->Eta();
    charge = RedTrack->Charge();
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

    fInclElecLP->Fill(fillSparse);
    fElecLPTrigger->Fill(pt);
    LPTrigger=kTRUE;

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
  }
   
  //save LP as Trigger (existens of LP and HFE are satisfied in function call)
  if (LPTrigger) fLPElecTrigger->Fill(ptH);
}



//_______________________

void AliAnalysisTaskHaHFECorrel::CorrelateLPMixedEvent(AliAODTrack* LPtrack, Float_t mult, Float_t zVtx ) 
{ 
  // leading Particle 
  Double_t pH=-9.,ptH=-9.,etaH =-9.,phiH=-9.;
  Int_t chargeH = 0, pdgH = 0, idH=-9;
  ptH = LPtrack->Pt();
  pH = LPtrack->P();
  phiH = LPtrack->Phi();
  etaH = LPtrack->Eta();
  chargeH = LPtrack->Charge();
  idH=LPtrack->GetID();

  AliEventPool* fPool;
  fPool = fPoolMgr->GetEventPool(mult, zVtx); // Get the buffer associated with the current centrality and z-vtx
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
Bool_t AliAnalysisTaskHaHFECorrel::InclElecTrackCuts(const AliAODVertex *pVietx,AliAODTrack *ietrack,Int_t nMother, Double_t listMother[])
{
  // track cuts for inclusive electrons
  Double_t d0z0[2]={-999,-999}, cov[3];
  fElectronTrackCuts->Fill(ietrack->Pt(), 0);
 
  fElectronTrackTPCNcls->Fill(ietrack->Pt(), ietrack->GetTPCNcls());
  fElectronTrackTPCNclsdEdx->Fill(ietrack->Pt(), ietrack->GetTPCsignalN());
  fElectronTrackTPCFrac->Fill(ietrack->Pt(), ietrack->GetTPCFoundFraction());
  fElectronTrackITSNcls->Fill(ietrack->Pt(), ietrack->GetITSNcls());
  //  if(ietrack->PropagateToDCA(pVietx, fVevent->GetMagneticField(), 20., d0z0, cov))
  //fElectronTrackDCA->Fill(TMath::Abs(d0z0[0]), TMath::Abs(d0z0[1]));

  if(!ietrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE;
  fElectronTrackCuts->Fill(ietrack->Pt(), 1);

  if(ietrack->Pt()<0.25) return kFALSE;
  fElectronTrackCuts->Fill(ietrack->Pt(), 2);

  if(ietrack->Eta()<fMinElectronEta) return kFALSE;
  if(ietrack->Eta()>fMaxElectronEta) return kFALSE;
  fElectronTrackCuts->Fill(ietrack->Pt(), 3);

  if(ietrack->GetTPCNcls() < fTPCnCut) return kFALSE;
  fElectronTrackCuts->Fill(ietrack->Pt(), 4);
  if(ietrack->GetTPCsignalN() < fTPCndEdxCut) return kFALSE ;
  fElectronTrackCuts->Fill(ietrack->Pt(), 5);
						
  if(!ietrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  fElectronTrackCuts->Fill(ietrack->Pt(), 6);
  //if(ietrack->GetTPCFoundFraction() < 0.6) return kFALSE;   
  fElectronTrackCuts->Fill(ietrack->Pt(), 7);

  Bool_t kinkmotherpass = kTRUE;
  for(Int_t kinkmother = 0; kinkmother < nMother; kinkmother++) {
    if(ietrack->GetID() == listMother[kinkmother]) {
      kinkmotherpass = kFALSE;
      continue;
    }
  }
  if(!kinkmotherpass) return kFALSE;
  fElectronTrackCuts->Fill(ietrack->Pt(), 8);

  if(ietrack->GetITSNcls() < fITSnCut) return kFALSE;
  fElectronTrackCuts->Fill(ietrack->Pt(), 9);

  if(!ietrack->IsOn(AliAODTrack::kITSrefit)) return kFALSE;    
  fElectronTrackCuts->Fill(ietrack->Pt(), 10);

  if(!(ietrack->HasPointOnITSLayer(0) || ietrack->HasPointOnITSLayer(1))) return kFALSE; // now kBoth
  fElectronTrackCuts->Fill(ietrack->Pt(), 11);

  if(!(ietrack->HasPointOnITSLayer(0) && ietrack->HasPointOnITSLayer(1))) return kFALSE; // now kBoth
  fElectronTrackCuts->Fill(ietrack->Pt(), 12);
  
  if(ietrack->PropagateToDCA(pVietx, fVevent->GetMagneticField(), 20., d0z0, cov))
    if(TMath::Abs(d0z0[0]) > 1 || TMath::Abs(d0z0[1]) > 2) return kFALSE;
  fElectronTrackCuts->Fill(ietrack->Pt(), 13);



 

  return kTRUE;

 }

 //________________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::InclElecPIDCuts(AliAODTrack* track, Bool_t IsPrimary) {
  // PID cuts
  Double_t fITSnSigma=-9.,fTOFnSigma=-9.,fTPCnSigma=-9.;
  Bool_t PassITSCut=kTRUE, PassTPCCut=kTRUE, PassTOFCut=kTRUE;

  fITSnSigma = fpidResponse->NumberOfSigmasITS(track, AliPID::kElectron);
  fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  fTOFnSigma = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
  
  
  if ((fITSnSigma < (-1.*fSigmaITScut) ) || (fITSnSigma > fSigmaITScut)) PassITSCut=kFALSE;
  if ((fTOFnSigma < (-1.*fSigmaTOFcut) ) || (fTOFnSigma > fSigmaTOFcut)) PassTOFCut=kFALSE;
  if ((fTPCnSigma<fSigmaTPCcut)  || (fTPCnSigma>3)) PassTPCCut=kFALSE; 


  if (IsPrimary) { //Fill PID histograms only in first round
    fHistITSnSig->Fill(track->P(),fITSnSigma);
    fHistTOFnSig->Fill(track->P(),fTOFnSigma);
    fHistTPCnSig->Fill(track->P(),fTPCnSigma);
    
    if (PassITSCut) fHistTPCnSigITScut->Fill(track->P(),fTPCnSigma);
    if (PassTOFCut) fHistTPCnSigTOFcut->Fill(track->P(),fTPCnSigma);
    if (PassITSCut && PassTOFCut) fHistTPCnSigITSTOFcut->Fill(track->P(),fTPCnSigma);

    if (PassTOFCut) fHadContPvsPt->Fill(track->P(), track->Pt());

    Double_t fillSparse[4];   
    fillSparse[0]=track->P();
    fillSparse[1]=track->Phi();
    fillSparse[2]=track->Eta();
    fillSparse[3]=fTPCnSigma;

    if (PassTOFCut) fHadContPPhiEtaTPC->Fill(fillSparse);

    fillSparse[1]=fITSnSigma;
    fillSparse[2]=fTOFnSigma;
    fHadContamination->Fill(fillSparse);

    fillSparse[0]=track->Pt();
    fHadContaminationPt->Fill(fillSparse);

    if (PassTOFCut && fIsMC) {
      Int_t MClabel=track->GetLabel();
      AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));     
      Int_t PDGCode = abs(MCParticle->GetPdgCode());
      fillSparse[0]=track->P();
      if (PDGCode==11) fillSparse[1]=0; // electron
      else if (PDGCode==13) fillSparse[1]=1; // muon
      else if (PDGCode==111 || PDGCode==211) fillSparse[1]=2; // p0 and p+
      else if (PDGCode==130 || PDGCode==310 || PDGCode==311 || PDGCode==321) fillSparse[1]=3; // K0L KOS K0 K0+
      else if (PDGCode==2212) fillSparse[1]=4; // proton
      else if (PDGCode==1000010020) fillSparse[1]=5;
      else fillSparse[1]=6;
      fillSparse[2]=fITSnSigma;
      fHadContMC->Fill(fillSparse);
      fillSparse[0]=track->Pt();
      fHadContMCPt->Fill(fillSparse);
    }
   }
 
  if (!fUseITS) PassITSCut = kTRUE;
  if (!PassITSCut || !PassTOFCut || !PassTPCCut) return kFALSE;
  return kTRUE;
}




// _________________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::PhotElecPIDCuts(AliAODTrack* track) {
  // associated particle variables
  Double_t fITSnSigmaAsso=-9.,fTOFnSigmaAsso=-9.,fTPCnSigmaAsso=-9.;
  
  // looser PID cuts
  fITSnSigmaAsso = fpidResponse->NumberOfSigmasITS(track, AliPID::kElectron);
  fTPCnSigmaAsso = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);

  if(TMath::Abs(fTPCnSigmaAsso)>fPhotElecSigmaTPCcut) return kFALSE;
  // if(TMath::Abs(fITSnSigmaAsso)>3) return kFALSE;

  return kTRUE;
  
}


//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::PhotElecTrackCuts(const AliAODVertex *pVaetx,AliAODTrack *aetrack,Int_t nMother, Double_t listMother[])
{
  // quality track cuts for associate tracks of photonic electrons
  if(!aetrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) return kFALSE;
  // if (aetrack->GetID()<0) cout << "ID <0 passes filter mask" << endl;
  if(aetrack->Pt() < fPhotElecPtCut) return kFALSE;
  if(TMath::Abs(aetrack->Eta())>0.9) return kFALSE;
  if(aetrack->GetTPCNcls() < fPhotElecTPCnCut) return kFALSE;
  if (fPhotElecITSrefitCut && !(aetrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
  if(!(aetrack->GetStatus()&AliESDtrack::kTPCrefit)) return kFALSE;

  return kTRUE;
  
}

 //_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::ChargedHadronTrackCuts(const AliAODVertex *pVaetx,AliAODTrack *Htrack,Int_t nMother, Double_t listMother[])
 {
   // quality track cuts for charged hadrons

   if(!Htrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE; 

   Bool_t kinkmotherpass = kTRUE;
   for(Int_t kinkmother = 0; kinkmother < nMother; kinkmother++) {
     if(Htrack->GetID() == listMother[kinkmother]) {
       kinkmotherpass = kFALSE;
       continue;
     }
   }
   if(!kinkmotherpass) return kFALSE;

   if(Htrack->Pt()<0.5) return kFALSE;
   if(Htrack->Eta()>fMaxHadronEta) return kFALSE;
   if(Htrack->Eta()<fMinHadronEta) return kFALSE;


  if(fHITSrefitCut && !(Htrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
  if(fHTPCrefitCut && !(Htrack->GetStatus()&AliESDtrack::kTPCrefit)) return kFALSE;

  if(Htrack->GetTPCNcls() < fHTPCnCut) return kFALSE;


  //if(Htrack->GetID()<0) return kFALSE; // not required due to appropriated filter mask

  
    
  return kTRUE;
    
}


//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::ChargedHadronPIDCuts(AliAODTrack *Htrack)
{   
  return kTRUE;   
}


Bool_t AliAnalysisTaskHaHFECorrel::AssoHadronPIDCuts(AliAODTrack *track) 
{
  // associated particle variables
  Double_t fTPCnSigmaAsso=-9.;
  
  // looser PID cuts
  fTPCnSigmaAsso = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  
  if(fTPCnSigmaAsso>fAssNonEleTPCcut) return kFALSE;
   
  return kTRUE;
}

void AliAnalysisTaskHaHFECorrel::MCEfficiencyCorrections(const AliAODVertex * RecVertex) {
  for (Int_t i=0; i<fMC->GetNumberOfTracks(); i++) {
    
    AliAODMCParticle *mcPart  = (AliAODMCParticle*)fMC->GetTrack(i);
    Double_t mcP, mcPt, mcPhi, mcEta, mcVtx[3];
    Int_t mcPDG, mcMotherID;
    mcP=mcPart->P();
    mcPt=mcPart->Pt();
    mcPhi=mcPart->Phi();
    mcEta=mcPart->Eta();
    mcPDG=abs(mcPart->GetPdgCode());
    mcMotherID=mcPart->GetMother();
    
    fMCheader->GetVertex(mcVtx);
    fCheckMCVertex->Fill(mcVtx[2], RecVertex->GetZ()); // TH2D

    // PI0 and Eta
    if (mcPart->IsPrimary()) {
      Double_t fillSparse[4];
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
	Double_t fillSparse[4];
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
      Double_t fillSparse[3];
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



void AliAnalysisTaskHaHFECorrel::EvaluateTaggingEfficiency(AliAODTrack * track, Int_t LSPartner, Int_t ULSPartner, Bool_t  trueULSPartner) {
  Int_t MClabel=track->GetLabel();
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
  fillSparse[0]=track->Pt(); // sparse will be filled with rec pt 
  fillSparse[1]=PDGCode;
  fillSparse[2]=PDGCodeMother;
  fillSparse[3]=PDGCodeGrandMother; 
  fillSparse[4]=PDGCodeGGMother;
  
 
  Double_t PtMotherWeight=1.; 
  Double_t ExcludePion[9] = {130,211, 310, 321, 2212, 3122,3222,3322,5122}; // K0L, PI+, K0S, K+, p, Lambd0, Sigma+, Xci0, lambda b0
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




//_________________________________
Bool_t  AliAnalysisTaskHaHFECorrel::CloneAndReduceTrackList(TObjArray* RedTracks, AliAODTrack* track, Int_t LSPartner, Int_t ULSPartner, Int_t* LSPartnerID, Int_t* ULSPartnerID, Bool_t trueULSPartner, Bool_t isPhotonic, Bool_t isHadron) {
  
  // Copied form AliAnalysisTaksPhiCorrelations and following instructions on
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PAGCorrelationsFAQ#In_Memory_Event_Mixing
  // Main difference is that track selection is in main to reduce computation time


  fCheckLSULS->Fill(LSPartner, ULSPartner);

  //  cout << "ParticleAdded" << endl;

  // clones a track list by using AliBasicParticle which uses much less memory (used for event mixing) -- modified class for my use
  // Clone only a certain pt bin on demand
  AliBasicParticleHaHFE * bparticle  = 0;
  bparticle = new AliBasicParticleHaHFE(track->GetID(), track->Eta(), track->Phi(), track->Pt(), track->Charge(), LSPartner, ULSPartner, LSPartnerID, ULSPartnerID, trueULSPartner, isPhotonic, isHadron);
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


void AliAnalysisTaskHaHFECorrel::FindV0Candidates(AliAODEvent *event) 
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

    // v0 Armenteros plot (QA)
    Float_t armVar[2] = {0.0,0.0};
    //  fV0cuts->Armenteros(v0, armVar);

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

  //
  // Clear the PID object arrays
  //
  
  fV0electrons->Clear();
  fV0pions->Clear();
  fV0protons->Clear();
  if (fV0tags!=0) delete[] fV0tags;
  fV0tags = 0;
}

void AliAnalysisTaskHaHFECorrel::TRDQA(Int_t RunNumber, const AliAODVertex *pVtx, Int_t nMother, Double_t listMother[]) {

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
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(VHtrack);
    if(!track) continue;

    if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;

    Bool_t passHFETrackCut=kFALSE;
    passHFETrackCut= InclElecTrackCuts(pVtx,track,nMother,listMother);

    Bool_t passHFEPIDCut=kFALSE;
    if (passHFETrackCut) passHFEPIDCut=  InclElecPIDCuts(track, kFALSE);
  
    Double_t fillSparse[7];
    fillSparse[0]=track->Pt();
    fillSparse[1]=track->Eta();
    fillSparse[2]=track->Phi();
    fillSparse[4]=track->Charge();
    fillSparse[5]=RunNumber;
    
    Int_t nTracklets=0;
    for (Int_t i=0; i<6; i++) {
      if (track->GetTRDmomentum(i)>0.01) { // TRDMomentumCriterion
	nTracklets++;
	fillSparse[3]=i+1;
	fTRDEtaPhi->Fill(fillSparse);
      }
    }
    
    if (nTracklets>0) {
      // if (track->GetTRDncls()/nTracklets<15) continue;
    }  
    
    if (passHFETrackCut)  fillSparse[2]=1;
    else fillSparse[2]=2;
    fillSparse[3]=fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    fillSparse[4]=nTracklets;
    fillSparse[5]=track->Charge();
    fillSparse[6]=RunNumber;
    fTRDNTracklets->Fill(fillSparse);
    
    if (track->Eta()>0.8 || track->Eta()<-0.8) continue;
    fTRDnTrackRun->Fill(nTracklets,  (Int_t)track->GetTRDntrackletsPID(), RunNumber);

    fillSparse[1]=nTracklets;
    if (passHFETrackCut)  fillSparse[2]=1;
    else fillSparse[2]=2;
    fillSparse[3]=fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    fillSparse[4]=track->Charge();
    fTRDSpectra->Fill(fillSparse);

    if (fIsMC) {
      Int_t MClabel=track->GetLabel();
      AliAODMCParticle* MCParticle = (AliAODMCParticle*) fMC->GetTrack(abs(MClabel));     
      Int_t PDGCode = abs(MCParticle->GetPdgCode());
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



void AliAnalysisTaskHaHFECorrel::FillV0Histograms(AliAODTrack* track, Int_t Species, Int_t RunNumber) {

  if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return;

  Double_t fillSparse[7];
  fillSparse[0]=track->Pt();
  fillSparse[1]=track->Eta();
  fillSparse[2]=Species;
  fillSparse[3]=fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  fillSparse[5]=track->Charge();
  fillSparse[6]=RunNumber;
  
  Int_t nTracklets=0;
  for (Int_t i=0; i<7; i++) {
    if (track->GetTRDmomentum(i)>0.01) {
      nTracklets++;
    }
  }
  
  fillSparse[4]=nTracklets;
  fTRDV0NTracklets->Fill(fillSparse);
  
  if (track->Eta()>0.8 || track->Eta()<-0.8) return;
  
  fillSparse[1]=nTracklets;
  fillSparse[2]=Species;
  fillSparse[3]= fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  fillSparse[4]=track->Charge();
  fTRDV0Spectra->Fill(fillSparse);
    
}
