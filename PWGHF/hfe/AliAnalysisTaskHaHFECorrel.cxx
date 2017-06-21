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
,fAssPtCut(0.5)
,fITSnCut(3)
,fAssTPCnCut(80)
,fTPCnCut(100)
,fHTPCnCut(100)
,fAssITSrefitCut(kTRUE)
,fHITSrefitCut(kTRUE)
,fHTPCrefitCut(kTRUE)
,fSigmaITScut(2.)
,fSigmaTOFcut(2.)
,fSigmaTPCcut(0.)
,fWeightSyst(kFALSE)
,fSystTOFcut(kFALSE)
,fEnablePileupRejVZEROTPCout(kFALSE)
,fRejectKinkMother(kFALSE)
,fOpeningAngleCut(1000.)
,fInvmassCut(0.14)
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
,fCuts(0)
,fIdentifiedAsOutInz(kFALSE)
,fPassTheEventCut(kFALSE)
,fIsMC(kFALSE)
,fIsAOD(kFALSE)
,fSetMassConstraint(kFALSE)
,fVz(0.0)
,fCFM(0)
,fPID(0)
,fPIDqa(0)
,fChi2Cut(3.5)
,fDCAcut(999)
,fWhichDecay(0)
,fPi0EtaWeight(1.)
,fOutputList(0)
,fNoEvents(0)
,fTrkpt(0)
,fMultiplicity(0)
,fHistITSnSig(0)
,fHistTOFnSig(0)
,fHistTPCnSig(0)
,fHistTPCnSigITScut(0)
,fHistTPCnSigTOFcut(0)
,fHistTPCnSigITSTOFcut(0)
,fElecTrigger(0)
,fInclElecPt(0)
,fULSElecPt(0)
,fLSElecPt(0)
,fInclElecPhi(0)
,fULSElecPhi(0)
,fLSElecPhi(0)
,fElecDphi(0)
,fULSElecDphi(0)
,fLSElecDphi(0)
,fULSElecDphiDiffMethod(0)
,fLSElecDphiDiffMethod(0)
,fElecHaTrigger(0)
,fElecHa(0)
,fLSElecHa(0)
,fULSElecHa(0)
,fElecHaLSNoPartner(0)
,fElecHaULSNoPartner(0)
,fElecHaMixedEvent(0)
,fLSElecHaMixedEvent(0)
,fULSElecHaMixedEvent(0)
,fElecLPTrigger(0)
,fElecLP(0)
,fLSElecLP(0)
,fULSElecLP(0)
,fElecLPMixedEvent(0)
,fLSElecLPMixedEvent(0)
,fULSElecLPMixedEvent(0)
,fInvmassLS(0)
,fInvmassULS(0)
,fOpeningAngleLS(0)
,fOpeningAngleULS(0)
,fCheckLSULS(0)
,fPi0Pt(0)
,fEtaPt(0)
{
    //Named constructor
    
    fPID = new AliHFEpid("hfePid");
    //fTrackCuts = new AliAODTrackCuts();
    //fAssTrackCuts = new AliAODTrackCuts();
    
    DefineInput(0, TChain::Class());
    //DefineInput(1, TList::Class());
    DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskHaHFECorrel::AliAnalysisTaskHaHFECorrel()
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisElecHadCorrel")
,fUseTender(kFALSE)
,fWhichPeriod(2016)
,fAssPtCut(0.5)
,fITSnCut(3)
,fAssTPCnCut(80)
,fTPCnCut(100)
,fHTPCnCut(100)
,fAssITSrefitCut(kTRUE)
,fHITSrefitCut(kTRUE)
,fHTPCrefitCut(kTRUE)
,fSigmaITScut(2.)
,fSigmaTOFcut(2.)
,fSigmaTPCcut(0.)
,fWeightSyst(kFALSE)
,fSystTOFcut(kFALSE)
,fEnablePileupRejVZEROTPCout(kFALSE)
,fRejectKinkMother(kFALSE)
,fOpeningAngleCut(1000.)
,fInvmassCut(0.14)
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
,fCuts(0)
,fIdentifiedAsOutInz(kFALSE)
,fPassTheEventCut(kFALSE)
,fIsMC(kFALSE)
,fIsAOD(kFALSE)
,fSetMassConstraint(kFALSE)
,fVz(0.0)
,fCFM(0)
,fPID(0)
,fPIDqa(0)
,fChi2Cut(3.5)
,fDCAcut(999)
,fWhichDecay(0)
,fPi0EtaWeight(1.)
,fOutputList(0)
,fNoEvents(0)
,fTrkpt(0)
,fMultiplicity(0)
,fHistITSnSig(0)
,fHistTOFnSig(0)
,fHistTPCnSig(0)
,fHistTPCnSigITScut(0)
,fHistTPCnSigTOFcut(0)
,fHistTPCnSigITSTOFcut(0)
,fElecTrigger(0)
,fInclElecPt(0)
,fULSElecPt(0)
,fLSElecPt(0)
,fInclElecPhi(0)
,fULSElecPhi(0)
,fLSElecPhi(0)
,fElecDphi(0)
,fULSElecDphi(0)
,fLSElecDphi(0)
,fULSElecDphiDiffMethod(0)
,fLSElecDphiDiffMethod(0)
,fElecHaTrigger(0)
,fElecHa(0)
,fLSElecHa(0)
,fULSElecHa(0)
,fElecHaLSNoPartner(0)
,fElecHaULSNoPartner(0)
,fElecHaMixedEvent(0)
,fLSElecHaMixedEvent(0)
,fULSElecHaMixedEvent(0)
,fElecLPTrigger(0)
,fElecLP(0)
,fLSElecLP(0)
,fULSElecLP(0)
,fElecLPMixedEvent(0)
,fLSElecLPMixedEvent(0)
,fULSElecLPMixedEvent(0)
,fInvmassLS(0)
,fInvmassULS(0)
,fOpeningAngleLS(0)
,fOpeningAngleULS(0)
,fCheckLSULS(0)
,fPi0Pt(0)
,fEtaPt(0)
{

    // Default constructor
    
    fPID = new AliHFEpid("hfePid");
    //fTrackCuts = new AliAODTrackCuts();
    //fAssTrackCuts = new AliAODTrackCuts();
    
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
      
  if(!fCuts){
    AliError("HFE cuts not available");
    return;
  }
  
  // not used anymore
  if(!fPID->IsInitialized()){
    // Initialize PID with the given run number
    AliWarning("PID not initialised, get from Run no");
    fPID->InitializePID(fAOD->GetRunNumber());
  }
    
  if(fUseTender){
    fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("AODFilterTracks"));
    fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));
    if (!fTracks_tender || !fCaloClusters_tender) return;
  }
    
  fMCarray  = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
  fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
  
 
  Int_t NpureMC = -1;
  if (fMCheader){
    TList *lh=fMCheader->GetCocktailHeaders();
    if(lh){
      AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(0); //  0 for HIJING
      NpureMC = gh->NProduced();
    }
  }
    
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();
      
  //fIsMC=kTRUE;
  if(fIsMC) fMC = MCEvent();
  if(fIsMC && fMC) fStack = fMC->Stack();
    
  fNoEvents->Fill(0);
    
  Int_t fNOtrks = fAOD->GetNumberOfTracks();
  if(fNOtrks<2) return;
  fNoEvents->Fill(1);
    
  // suggested by DPG to remove outliers
  const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
  const AliAODVertex *spdVtx = fAOD->GetPrimaryVertexSPD();
    
  Double_t pVtxZ = -999;
  pVtxZ = pVtx->GetZ();
    
  if (pVtx->GetNContributors()<2 || spdVtx->GetNContributors()<1) return; // one of vertices is missing
  // pVtx -  vertex not reconstructed from tracks/tracklets; vtx position is the mean vertex in x,y and z=0    

  double covTrc[6],covSPD[6];
  pVtx->GetCovarianceMatrix(covTrc);
  spdVtx->GetCovarianceMatrix(covSPD);
  double dz = pVtx->GetZ()-spdVtx->GetZ();
  double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
  double errTrc = TMath::Sqrt(covTrc[5]);
  double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
  if (TMath::Abs(dz)>0.2 || nsigTot>10 || nsigTrc>20) return; // bad vertexing
  
  if(TMath::Abs(pVtxZ)>10) return;
  fNoEvents->Fill(2);
    
  fpidResponse = fInputHandler->GetPIDResponse();
  if(!fpidResponse){
    AliDebug(1, "Using default PID Response");
    fpidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class());
  }
  fPID->SetPIDResponse(fpidResponse);
  fCFM->SetRecEventInfo(fAOD);

    
  //multiplicity -- commented because problem with LHC13
  Float_t mult = -1.;
  // fAOD->GetMultiplicity();
  fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
  mult = fMultSelection->GetMultiplicityPercentile("V0M", kFALSE);
  fMultiplicity->Fill(mult);
    
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

  ///////////////////////
  // Preparational Tasks
  ///////////////////////

  // List of HFE for analysis and mixed event
  TObjArray* RedTracksHFE = new TObjArray;
  RedTracksHFE->SetOwner(kTRUE);

  // FirstLoop: Find leading particle (LP) and HFE (check for LS/ULS partner)
  Int_t LeadingParticle;
  AliAODTrack* LPtrack;
  LPtrack=FindLPAndHFE(RedTracksHFE, pVtx,nMotherKink,listofmotherkink);

  //////////////////////////
  // Main analysis
  ///////////////////////////

  // LP - two different functions for Same Event and mixed Event
  if (fLParticle && RedTracksHFE->GetEntriesFast()>0) CorrelateLP(LPtrack, RedTracksHFE);
  if (fLParticle && fMixedEvent) CorrelateLPMixedEvent(LPtrack, mult, pVtx->GetZ());
 

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


  PostData(1, fOutputList);
}

//_________________________________________
void AliAnalysisTaskHaHFECorrel::UserCreateOutputObjects()
{
  //--- Check MC
  if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){
    fIsMC = kTRUE;
    printf("+++++ MC Data available");
  }
  
  //--------Initialize PID
  fPID->SetHasMCData(fIsMC);
  if(!fPID->GetNumberOfPIDdetectors())
    {
      fPID->AddDetector("ITS", 0);
      fPID->AddDetector("TOF", 1);
      fPID->AddDetector("TPC", 2);
    }
    
  fPID->SortDetectors();
  fPIDqa = new AliHFEpidQAmanager();
  fPIDqa->Initialize(fPID);
    
  //--------Initialize correction Framework and Cuts -- currently not used
  fCFM = new AliCFManager;
  const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack;
  fCFM->SetNStepParticle(kNcutSteps);
  for(Int_t istep = 0; istep < kNcutSteps; istep++)
    fCFM->SetParticleCutsList(istep, NULL);
  
  if(!fCuts){
    AliWarning("Cuts not available. Default cuts will be used");//same as in the config file (to be be removed)
    fCuts = new AliHFEcuts;
    fCuts->CreateStandardCuts();
    fCuts->SetMinNClustersTPC(fTPCnCut);
    fCuts->SetMinRatioTPCclusters(0.6);
    fCuts->SetMaxChi2perClusterTPC(3.5);
    fCuts->SetTPCmodes(AliHFEextraCuts::kFound, AliHFEextraCuts::kFoundOverFindable);
    //hfecuts->SetMinNClustersITS(ITSncut); // it depends on pt
    fCuts->SetCutITSpixel(AliHFEextraCuts::kAny);
    fCuts->SetCheckITSLayerStatus(kFALSE);
    fCuts->SetVertexRange(10.);
    fCuts->SetPtRange(1.5, 50);
    fCuts->SetMaxImpactParam(2.4,3.2); // radial, z
  }
    
  fCuts->SetAOD();
  fCuts->Initialize(fCFM);
    
  //---------Output Tlist
  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputList->Add(fPIDqa->MakeList("PIDQA"));
    
  fNoEvents = new TH1F("fNoEvents","",3,0,3);
  fOutputList->Add(fNoEvents);
    
  fTrkpt = new TH2F("fTrkpt","track pt",200,0,20,3,0,3);
  fOutputList->Add(fTrkpt);
    
  fMultiplicity = new TH1F("fMultiplicity","Centrality",100,0,100) ;
  fMultiplicity->Sumw2();
  fOutputList->Add(fMultiplicity);
    
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
  

  int nbin_Asso = 4;
  double bin_Asso[7] = {0.25,1.5,3,6,10};
    
  fULSElecPt = new TH1F("fULSElecPt", "fULSElecPt",200,0,20);
  fULSElecPt->Sumw2();
  fOutputList->Add(fULSElecPt);
    
  fLSElecPt = new TH1F("fLSElecPt", "fLSElecPt",200,0,20);
  fLSElecPt->Sumw2();
  fOutputList->Add(fLSElecPt);
    
  fInclElecPt = new TH1F("fInclElecPt", "fInclElecPt", 200,0,20);
  fInclElecPt->Sumw2();
  fOutputList->Add(fInclElecPt);

  fElecTrigger = new TH1F("fElecTrigger","fElecTrigger", 30, 0, 30);
  fOutputList->Add(fElecTrigger);
    
  fInclElecPhi = new TH2F("fInclElecPhi", "fInclElecPhi",nbin_Asso,bin_Asso,100,0,TMath::TwoPi());
  fInclElecPhi->Sumw2();
  fOutputList->Add(fInclElecPhi);
  
  fULSElecPhi= new TH2F("fULSElecPhi", "fULSElecPhi",nbin_Asso,bin_Asso,100,0,TMath::TwoPi());
  fULSElecPhi->Sumw2();
  fOutputList->Add(fULSElecPhi);
  
  fLSElecPhi= new TH2F("fLSElecPhi", "fLSElecPhi",nbin_Asso,bin_Asso,100,0,TMath::TwoPi());
  fLSElecPhi->Sumw2();
  fOutputList->Add(fLSElecPhi);
  
  fElecDphi = new TH2F("fElecDphi", "fElecDphi",nbin_Asso,bin_Asso,100,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fElecDphi->Sumw2();
  fOutputList->Add(fElecDphi);
  
  fULSElecDphi = new TH2F("fULSElecDphi", "fULSElecDphi",nbin_Asso,bin_Asso,100,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fULSElecDphi->Sumw2();
  fOutputList->Add(fULSElecDphi);
  
  fLSElecDphi = new TH2F("fLSElecDphi", "fLSElecDphi",nbin_Asso,bin_Asso,100,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fLSElecDphi->Sumw2();
  fOutputList->Add(fLSElecDphi);
  
  fULSElecDphiDiffMethod = new TH2F("fULSElecDphiDiffMethod", "fULSElecDphiDiffMethod",nbin_Asso,bin_Asso,100,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fULSElecDphi->Sumw2();
  fOutputList->Add(fULSElecDphiDiffMethod);
  
  fLSElecDphiDiffMethod = new TH2F("fLSElecDphiDiffMethod", "fLSElecDphiDiffMethod",nbin_Asso,bin_Asso,100,-0.5*TMath::Pi(),1.5*TMath::Pi());
  fLSElecDphi->Sumw2();
  fOutputList->Add(fLSElecDphiDiffMethod);
  
  // THnSparse Binning			    
  Int_t bin[4] = {99,79,32,50}; //ptH, ptE, Dphi, Deta
  Double_t xmin[4] = {0.5,0.5,-TMath::Pi()/2,-1.7};
  Double_t xmax[4] = {50,40,(3*TMath::Pi())/2,1.7};
		    
  fElecHaTrigger = new TH1F("fElecHaTrigger", "fElecHaTrigger", bin[1], xmin[1], xmax[1]);
  fElecHaTrigger->Sumw2();
  fOutputList->Add(fElecHaTrigger);

  fElecHa = new THnSparseD("fEleHa", "Sparse for Elec-Had : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecHa->Sumw2();
  fOutputList->Add(fElecHa);

  fLSElecHa = new THnSparseD("fEleHaLS", "Sparse for LSElec-Had : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fLSElecHa->Sumw2();
  fOutputList->Add(fLSElecHa);

  fULSElecHa = new THnSparseD("fEleHaULS", "Sparse for ULSElec-Had : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fULSElecHa->Sumw2();
  fOutputList->Add(fULSElecHa);
  
  fElecHaLSNoPartner = new THnSparseD("fEleHaLSNoPartner", "Sparse for LSElec-Had with no Partner: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecHaLSNoPartner->Sumw2();
  fOutputList->Add(fElecHaLSNoPartner);

  fElecHaULSNoPartner = new THnSparseD("fEleHaULSNoPartner", "Sparse for ULSElec-Had with no Partner: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecHaULSNoPartner->Sumw2();
  fOutputList->Add(fElecHaULSNoPartner);

  fElecHaMixedEvent = new THnSparseD("fEleHaMixedEv", "Sparse for Elec-Had MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecHaMixedEvent->Sumw2();
  fOutputList->Add(fElecHaMixedEvent);
  
  fLSElecHaMixedEvent = new THnSparseD("fEleHaLSMixedEv", "Sparse for LSElec-Had MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fLSElecHaMixedEvent->Sumw2();
  fOutputList->Add(fLSElecHaMixedEvent);
  
  fULSElecHaMixedEvent = new THnSparseD("fEleHaULSMixedEv", "Sparse for ULSElec-Had MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fULSElecHaMixedEvent->Sumw2();
  fOutputList->Add(fULSElecHaMixedEvent);

  fElecLPTrigger = new TH1F("fElecLPTrigger", "fElecLPTrigger", bin[1], xmin[1], xmax[1]);
  fElecLPTrigger->Sumw2();
  fOutputList->Add(fElecLPTrigger);						      
  fElecLP = new THnSparseD("fEleLP", "Sparse for Elec-LP : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecLP->Sumw2();
  fOutputList->Add(fElecLP);

  fLSElecLP = new THnSparseD("fEleLPLS", "Sparse for LSElec-LP : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fLSElecLP->Sumw2();
  fOutputList->Add(fLSElecLP);
  
  fULSElecLP = new THnSparseD("fEleLPULS", "Sparse for ULSElec-LP : PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fULSElecLP->Sumw2();
  fOutputList->Add(fULSElecLP);
   
  fElecLPMixedEvent = new THnSparseD("fEleLPMixedEv", "Sparse for Elec-LP MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fElecLPMixedEvent->Sumw2();
  fOutputList->Add(fElecLPMixedEvent);
  
  fLSElecLPMixedEvent = new THnSparseD("fEleLPLSMixedEv", "Sparse for LSElec-LP MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fLSElecLPMixedEvent->Sumw2();
  fOutputList->Add(fLSElecLPMixedEvent);

  fULSElecLPMixedEvent = new THnSparseD("fEleLPULSMixedEv", "Sparse for ULSElec-LP MixEvent: PtH, PtE, Dphi, Deta", 4, bin, xmin, xmax);
  fULSElecLPMixedEvent->Sumw2();
  fOutputList->Add(fULSElecLPMixedEvent);
    
  fInvmassLS = new TH2F("fInvmassLS", "fInvmassLS", 500,0,0.5,200,0,20);
  fInvmassLS->Sumw2();
  fOutputList->Add(fInvmassLS);
  
  fInvmassULS = new TH2F("fInvmassULS", "fInvmassULS", 500,0,0.5,200,0,20);
  fInvmassULS->Sumw2();
  fOutputList->Add(fInvmassULS);
    
  fOpeningAngleLS = new TH2F("fOpeningAngleLS","fOpeningAngleLS",100,0,1,200,0,20);
  fOpeningAngleLS->Sumw2();
  fOutputList->Add(fOpeningAngleLS);
    
  fOpeningAngleULS = new TH2F("fOpeningAngleULS","fOpeningAngleULS",100,0,1,200,0,20);
  fOpeningAngleULS->Sumw2();
  fOutputList->Add(fOpeningAngleULS);
    
  fPi0Pt = new TH2F("fPi0Pt", "Pi0 pT",4,0,4,200,0,20);
  fPi0Pt->Sumw2();
  fOutputList->Add(fPi0Pt);
  
  fEtaPt = new TH2F("fEtaPt", "Eta pT",4,0,4,200,0,20);
  fEtaPt->Sumw2();
  fOutputList->Add(fEtaPt);

  Int_t poolSize = 1000;
  Int_t trackDepth = 2000; 
  Int_t nMultBins = 10;
  Double_t multBins[]={0,10,20,30,40,50,60,70,80,90,100};
  Int_t nZVtxBins=9;
  Double_t zVtxBins[]={-10,-7,-5,-3,-1,1,3,5,7,10};
  fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, nMultBins, (Double_t*) multBins, nZVtxBins, (Double_t*) zVtxBins);
  fPoolMgr->Validate();
    
  fCheckLSULS = new TH2F("fCheckLSULS", "LSULS",5,0,5,5,0,5);
  fCheckLSULS->Sumw2();
  fOutputList->Add(fCheckLSULS);
    
  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHaHFECorrel::Terminate(Option_t *)
{
  // Info("Terminate");
  AliAnalysisTaskSE::Terminate();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
  // Currently not used
  // Check single track cuts for a given cut step
  const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
  if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
  return kTRUE;
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
AliAODTrack*  AliAnalysisTaskHaHFECorrel::FindLPAndHFE( TObjArray* RedTracks, const AliVVertex *pVtx, Int_t nMother, Double_t listMother[])
{
  AliAODTrack* LPtrack=0;
  fLParticle=kTRUE;
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();
  
  // associated particle
  Double_t pH=-9.,ptH=-9.,etaH =-9.,phiH=-9.;
  Int_t chargeH = 0, pdgH = 0;

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
    
    // find hadron with the largest pT -> leading particle
    if (passHadTrackCut && passHadPIDCut && track->Pt()>ptH) {
      ptH = track->Pt();
      pH = track->P();
      phiH =track->Phi();
      etaH = track->Eta();
      chargeH = track->Charge();
      LPtrack=track;
    }  

    if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;

    Bool_t passHFETrackCut=kFALSE;
    passHFETrackCut= InclElecTrackCuts(pVtx,track,nMother,listMother);
    if (!passHFETrackCut) continue;
    fTrkpt->Fill(pt,1); // after track cuts

    Bool_t passHFEPIDCut=kFALSE;
    passHFEPIDCut=  InclElecPIDCuts(track, kTRUE);
    if (!passHFEPIDCut) continue;
    fTrkpt->Fill(pt,2);
    
    Int_t lsPartner=0, ulsPartner=0;
    if (passHFETrackCut && passHFEPIDCut) { // if HFE is found, look for ls and uls partner
      FindPhotonicPartner(jTracks, track, pVtx, nMother, listMother, lsPartner, ulsPartner);
      fInclElecPt->Fill(pt);
      fInclElecPhi->Fill(pt,phi); // phi of electron candidates
    }

    // store only essential informations of these electrons for later correlations and mixed event
    CloneAndReduceTrackList(RedTracks, track, lsPartner, ulsPartner);
    					   
  }

  // if no leading particle is found, mark this in flag fLParticle;
  if (!LPtrack) {
    fLParticle=kFALSE;
  }
 
  return LPtrack;
}

//_________________________________________
void AliAnalysisTaskHaHFECorrel::FindPhotonicPartner(Int_t iTracks, AliAODTrack* track, const AliVVertex *pVtx, Int_t nMother, Double_t listMother[], Int_t &lsPartner, Int_t &ulsPartner) {
  
  lsPartner=0;
  ulsPartner=0;
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

    if(fFlagLS){
      fOpeningAngleLS->Fill(openingAngle,pt);
      fInvmassLS->Fill(mass,pt);
    }
    if(fFlagULS){
      fOpeningAngleULS->Fill(openingAngle,pt);
      fInvmassULS->Fill(mass,pt);
    }
        
    if(mass<fInvmassCut && fFlagULS){
      fULSElecPhi->Fill(pt,phi);
      fULSElecPt->Fill(pt);
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      fULSElecDphi->Fill(pt,dphi);
      ulsPartner++;
    }
        
    if(mass<fInvmassCut && fFlagLS){
      fLSElecPhi->Fill(pt,phi);
      fLSElecPt->Fill(pt);
      dphi = GetDeltaPhi(phi,phiAsso); // electron-electron dphi
      fLSElecDphi->Fill(pt,dphi);
      lsPartner++;
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

void AliAnalysisTaskHaHFECorrel::CorrelateHadron(TObjArray* RedTracksHFE,  const AliVVertex* pVtx, Int_t nMother, Double_t listMother[], Float_t mult) {

  Int_t HadronTrigger=kFALSE;
    
  Int_t ntracks = -999;
  if(!fUseTender) ntracks = fVevent->GetNumberOfTracks();
  if(fUseTender)  ntracks = fTracks_tender->GetEntries();

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

    if (fMixedEvent) CorrelateHadronMixedEvent(track, mult, pVtx->GetZ());

    // loop over all electrons
    for (Int_t k=0; k<RedTracksHFE->GetEntriesFast(); k++) {

      AliBasicParticleHaHFE *RedTrack = (AliBasicParticleHaHFE*) RedTracksHFE->At(k);

      if (RedTrack->ID()==idH) continue;

      // tagged particle
      Double_t p=-9.,pt=-9.,eta =-9.,phi=-9.;
      Int_t charge = 0, pdg = 0, ls=0, uls=0;    
      pt = RedTrack->Pt();
      // p = RedTrack->P();
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

      fElecHa->Fill(fillSparse);
      HadronTrigger=kTRUE;
      for (Int_t j=0; j<ls; j++) fLSElecHa->Fill(fillSparse);
      for (Int_t j=0; j<uls; j++) fULSElecHa->Fill(fillSparse);
    }

    // Fill Trigger Hadrons
    if (HadronTrigger) fElecHaTrigger->Fill(ptH);

  }//end of track loop
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



//_______________________

void AliAnalysisTaskHaHFECorrel::CorrelateLP(AliAODTrack* LPtrack, TObjArray *RedTracksHFE) 
{ 
  Bool_t LPTrigger=kFALSE;
  
  // leading Particle neq 
  Double_t pH=-9.,ptH=-9.,etaH =-9.,phiH=-9.;
  Int_t chargeH = 0, pdgH = 0, idH=-9;
  ptH = LPtrack->Pt();
  pH = LPtrack->P();
  phiH = LPtrack->Phi();
  etaH = LPtrack->Eta();
  chargeH = LPtrack->Charge();
  idH = LPtrack->GetID();

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

    fElecLP->Fill(fillSparse);
    LPTrigger=kTRUE;
    for (Int_t j=0; j<ls; j++) fLSElecLP->Fill(fillSparse);
    for (Int_t j=0; j<uls; j++) fULSElecLP->Fill(fillSparse);
  }
   
  //save LP as Trigger (existens of LP and HFE are satisfied in function call)
  if (LPTrigger) fElecLPTrigger->Fill(ptH);
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
Bool_t AliAnalysisTaskHaHFECorrel::InclElecTrackCuts(const AliVVertex *pVietx,AliAODTrack *ietrack,Int_t nMother, Double_t listMother[])
{
    // track cuts for inclusive electrons
    
    if(!ietrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE;

    if(TMath::Abs(ietrack->Eta())>0.8) return kFALSE;
    
    if(ietrack->GetTPCNcls() < fTPCnCut) return kFALSE;
    
    if(ietrack->Pt()<3 && ietrack->GetITSNcls() < fITSnCut) return kFALSE;
    if(ietrack->Pt()>=3 && ietrack->GetITSNcls() < 3) return kFALSE;
    
    if(!ietrack->IsOn(AliAODTrack::kITSrefit)) return kFALSE;
    
    if(!ietrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
    
    if(!(ietrack->HasPointOnITSLayer(0) || ietrack->HasPointOnITSLayer(1))) return kFALSE;
    
    if(ietrack->GetTPCFoundFraction() < 0.6) return kFALSE;
    
    if (ietrack->Pt()<0.25) return kFALSE;
    
    Bool_t kinkmotherpass = kTRUE;
    for(Int_t kinkmother = 0; kinkmother < nMother; kinkmother++) {
        if(ietrack->GetID() == listMother[kinkmother]) {
            kinkmotherpass = kFALSE;
            continue;
        }
    }
    if(!kinkmotherpass) return kFALSE;
    
    Double_t d0z0[2]={-999,-999}, cov[3];
    
    if(ietrack->PropagateToDCA(pVietx, fVevent->GetMagneticField(), 20., d0z0, cov))
        if(TMath::Abs(d0z0[0]) > 2.4 || TMath::Abs(d0z0[1]) > 3.2) return kFALSE;
    
    return kTRUE;
    
}

//________________________________________________
 Bool_t AliAnalysisTaskHaHFECorrel::InclElecPIDCuts(AliAODTrack* track, Bool_t IsPrimary) {
  // PID cuts
  Double_t fITSnSigma=-9.,fTOFnSigma=-9.,fTPCnSigma=-9.;
        
  fITSnSigma = fpidResponse->NumberOfSigmasITS(track, AliPID::kElectron);
  fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  fTOFnSigma = fpidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);

  double p=-9;
  p=track->P();
        
  if (IsPrimary) {
    fHistITSnSig->Fill(p,fITSnSigma);
    fHistTOFnSig->Fill(p,fTOFnSigma);
    fHistTPCnSig->Fill(p,fTPCnSigma);
    
    if (fITSnSigma>-2 && fITSnSigma<2) fHistTPCnSigITScut->Fill(p,fTPCnSigma);
    if (fTOFnSigma>-2 && fTOFnSigma<2) fHistTPCnSigTOFcut->Fill(p,fTPCnSigma);
    if (fITSnSigma>-2 && fITSnSigma<2 && fTOFnSigma>-2 && fTOFnSigma<2) fHistTPCnSigITSTOFcut->Fill(p,fTPCnSigma);
  }
        
  if ((fITSnSigma < (-1.*fSigmaITScut) ) || (fITSnSigma > fSigmaITScut)) return kFALSE;
  if ((fTOFnSigma < (-1.*fSigmaTOFcut) ) || (fTOFnSigma > fSigmaTOFcut)) return kFALSE;
  if ((fTPCnSigma<fSigmaTPCcut)  || (fTPCnSigma>3)) return kFALSE;

}


// _________________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::PhotElecPIDCuts(AliAODTrack* track) {
  // associated particle variables
  Double_t fITSnSigmaAsso=-9.,fTOFnSigmaAsso=-9.,fTPCnSigmaAsso=-9.;
  
  // looser PID cuts
  fITSnSigmaAsso = fpidResponse->NumberOfSigmasITS(track, AliPID::kElectron);
  fTPCnSigmaAsso = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  
  if(TMath::Abs(fTPCnSigmaAsso)>3) return kFALSE;
  if(TMath::Abs(fITSnSigmaAsso)>3) return kFALSE;
  
  return kTRUE;
  
}


//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::PhotElecTrackCuts(const AliVVertex *pVaetx,AliAODTrack *aetrack,Int_t nMother, Double_t listMother[])
{
    // quality track cuts for associate tracks of photonic electrons
    if(!aetrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) return kFALSE;
    if(aetrack->Pt() < fAssPtCut) return kFALSE;
    if(TMath::Abs(aetrack->Eta())>0.9) return kFALSE;
    if(aetrack->GetTPCNcls() < fAssTPCnCut) return kFALSE;
    if (fAssITSrefitCut && !(aetrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
    if(!(aetrack->GetStatus()&AliESDtrack::kTPCrefit)) return kFALSE;
    if (aetrack->Pt()<0.25) return kFALSE;
    
    return kTRUE;
    
}

//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::ChargedHadronTrackCuts(const AliVVertex *pVaetx,AliAODTrack *Htrack,Int_t nMother, Double_t listMother[])
{
  // quality track cuts for charged hadrons
  Bool_t kinkmotherpass = kTRUE;
  for(Int_t kinkmother = 0; kinkmother < nMother; kinkmother++) {
    if(Htrack->GetID() == listMother[kinkmother]) {
      kinkmotherpass = kFALSE;
      continue;
    }
  }
  if(!kinkmotherpass) return kFALSE;
  if(fHITSrefitCut && !(Htrack->GetStatus()&AliESDtrack::kITSrefit)) return kFALSE;
  if(fHTPCrefitCut && !(Htrack->GetStatus()&AliESDtrack::kTPCrefit)) return kFALSE;
  if(TMath::Abs(Htrack->Eta())>0.9) return kFALSE;
  if(Htrack->GetTPCNcls() < fHTPCnCut) return kFALSE;
  if(Htrack->Pt()<0.5) return kFALSE;

  if(Htrack->GetID()<0) return kFALSE;

  //    if(!aetrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) return kFALSE; ??
    
  return kTRUE;
    
}


//_________________________________________
Bool_t AliAnalysisTaskHaHFECorrel::ChargedHadronPIDCuts(AliAODTrack *Htrack)
{   
  return kTRUE;   
}

//_________________________________
Bool_t  AliAnalysisTaskHaHFECorrel::CloneAndReduceTrackList(TObjArray* RedTracks, AliAODTrack* track, Int_t LSPartner, Int_t ULSPartner) {
  
  // Copied form AliAnalysisTaksPhiCorrelations and following instructions on
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PAGCorrelationsFAQ#In_Memory_Event_Mixing
  // Main difference is that track selection is in main to reduce computation time


  fCheckLSULS->Fill(LSPartner, ULSPartner);

  //  cout << "ParticleAdded" << endl;

  // clones a track list by using AliBasicParticle which uses much less memory (used for event mixing) -- modified class for my use
  // Clone only a certain pt bin on demand
  AliBasicParticleHaHFE * bparticle  = 0;
  bparticle = new AliBasicParticleHaHFE(track->GetID(), track->Eta(), track->Phi(), track->Pt(), track->Charge(), LSPartner, ULSPartner);
  if (!bparticle) return kFALSE;
  RedTracks->Add(bparticle);
  return kTRUE;

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
