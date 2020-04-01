#include "AliAnalysisTaskJpsiJetFilter.h"

#include "TChain.h"
#include "THashList.h"

#include <AliLog.h>
#include "AliAODCaloCluster.h"
#include "AliAODBranchReplicator.h"

#include "AliDielectronVarManager.h"
#include "AliDielectronVarCuts.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronPairLegCuts.h"

class AliDielectronVarManager;
class AliDielectronVarCuts;
class AliDielectronTrackCuts;
class AliDielectronPairLegCuts;

ClassImp(AliAnalysisTaskJpsiJetFilter)

AliAnalysisTaskJpsiJetFilter::AliAnalysisTaskJpsiJetFilter() : 
  AliAnalysisTaskSE(),
  fIsToMerge(kFALSE),
  fIsToReplace(kFALSE),
  fOutputFileName("AliAOD.Dielectron.root"),
  fExtAOD(0x0),
  fSPD(0x0),
  fEMCALTrigger(0x0),
  fPHOSTrigger(0x0),
  fEMCalCells(0x0),
  fPHOSCells(0x0),
  fAODZDC(0x0),
  fAODAD(0x0),
  fAODTZERO(0x0),
  fPairs(0x0),
  fDaughters(0x0),
  fJets(0x0),
  fDielectron(0),
  fSelectPhysics(kTRUE),
  fTriggerMask(AliVEvent::kMB),
  fExcludeTriggerMask(0),
  fTriggerOnV0AND(kFALSE),
  fRejectPileup(kFALSE),
  fEventStat(0x0),
  fTriggerLogic(kAny),
  fTriggerAnalysis(0x0),
  fStoreLikeSign(kFALSE),
  fStoreRotatedPairs(kFALSE),
  fStoreEventsWithSingleTracks(kFALSE),
  fCreateNanoAOD(kFALSE),
  fEventFilter(0x0)
{}

AliAnalysisTaskJpsiJetFilter::AliAnalysisTaskJpsiJetFilter(const char* name) : 
  AliAnalysisTaskSE(name),
  fIsToMerge(kFALSE),
  fIsToReplace(kFALSE),
  fOutputFileName("AliAOD.Dielectron.root"),
  fExtAOD(0x0),
  fSPD(0x0),
  fEMCALTrigger(0x0),
  fPHOSTrigger(0x0),
  fEMCalCells(0x0),
  fPHOSCells(0x0),
  fAODZDC(0x0),
  fAODAD(0x0),
  fAODTZERO(0x0),
  fPairs(0x0),
  fDaughters(0x0),
  fJets(0x0),
  fDielectron(0),
  fSelectPhysics(kTRUE),
  fTriggerMask(AliVEvent::kMB),
  fExcludeTriggerMask(0),
  fTriggerOnV0AND(kFALSE),
  fRejectPileup(kFALSE),
  fEventStat(0x0),
  fTriggerLogic(kAny),
  fTriggerAnalysis(0x0),
  fStoreLikeSign(kFALSE),
  fStoreRotatedPairs(kFALSE),
  fStoreEventsWithSingleTracks(kFALSE),
  fCreateNanoAOD(kFALSE),
  fEventFilter(0x0)
{
  DefineInput(0,TChain::Class());
  DefineOutput(1, THashList::Class());
  DefineOutput(2, TH1D::Class());
};

AliAnalysisTaskJpsiJetFilter::~AliAnalysisTaskJpsiJetFilter()
{
  // Destructor
  if(fDielectron)
    delete fDielectron;
  if(fEventStat)
    delete fEventStat;
  if(fTriggerAnalysis)
    delete fTriggerAnalysis;
  if(fEventFilter)
    delete fEventFilter;
  if(fExtAOD) delete fExtAOD;
  if(fSPD) delete fSPD;
  if(fEMCALTrigger) delete fEMCALTrigger;
  if(fPHOSTrigger) delete fPHOSTrigger;
  if(fEMCalCells) delete fEMCalCells;
  if(fPHOSCells) delete fPHOSCells;
  if(fAODZDC) delete fAODZDC;
  if(fAODAD) delete fAODAD;
  if(fAODTZERO) delete fAODTZERO;
  if(fPairs){
    fPairs->Clear("C");
    delete fPairs;
  }
  if(fDaughters){
    fDaughters->Clear("C");
    delete fDaughters;
  }
  if(fJets){
    fJets->Clear("C");
    delete fJets;
  }
}

void AliAnalysisTaskJpsiJetFilter::UserCreateOutputObjects()
{
  if (!fDielectron){
    InitDielectron();
    InitHistogramsForDielectron("DieFilterHistos");
  }
  if(fStoreRotatedPairs) fDielectron->SetStoreRotatedPairs(kTRUE);
  fDielectron->SetDontClearArrays();
  fDielectron->SetHasMC(kFALSE);
  fDielectron->Init();

  Int_t nbins=kNbinsEvent+2;
  if (!fEventStat){
    fEventStat=new TH1D("hEventStat","Event statistics",nbins,0,nbins);
    fEventStat->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");

    fEventStat->GetXaxis()->SetBinLabel(3,"Bin3 not used");
    fEventStat->GetXaxis()->SetBinLabel(4,"Bin4 not used");
    fEventStat->GetXaxis()->SetBinLabel(5,"Bin5 not used");

    if(fTriggerOnV0AND) fEventStat->GetXaxis()->SetBinLabel(3,"V0and triggers");
    if (fEventFilter) fEventStat->GetXaxis()->SetBinLabel(4,"After Event Filter");
    if (fRejectPileup) fEventStat->GetXaxis()->SetBinLabel(5,"After Pileup rejection");

    fEventStat->GetXaxis()->SetBinLabel((kNbinsEvent+1),Form("#splitline{1 candidate}{%s}",fDielectron->GetName()));
    fEventStat->GetXaxis()->SetBinLabel((kNbinsEvent+2),Form("#splitline{With >1 candidate}{%s}",fDielectron->GetName()));
  }

  // Create branch to nano AOD - track/vertex/cluster
  TClonesArray *nanoAODTracks = new TClonesArray("AliAODTrack",500);
  nanoAODTracks->SetName("tracks");
  fExtAOD->AddBranch("TClonesArray", &nanoAODTracks);
  TClonesArray *nanoAODVertices = new TClonesArray("AliAODVertex",500);
  nanoAODVertices->SetName("vertices");
  fExtAOD->AddBranch("TClonesArray", &nanoAODVertices);
  TClonesArray *nanoAODCaloCluster = new TClonesArray("AliAODCaloCluster",500);
  nanoAODCaloCluster->SetName("caloClusters");
  fExtAOD->AddBranch("TClonesArray", &nanoAODCaloCluster);
  // Create branch for other objects (not TClonesArray) 
  fSPD = new AliAODTracklets;
  fSPD->SetName("tracklets");
  fExtAOD->AddBranch("AliAODTracklets",&fSPD);

  fEMCALTrigger = new AliAODCaloTrigger;
  fEMCALTrigger->SetName("emcalTrigger");
  fExtAOD->AddBranch("AliAODCaloTrigger",&fEMCALTrigger);
  fPHOSTrigger = new AliAODCaloTrigger;
  fPHOSTrigger->SetName("phosTrigger");
  fExtAOD->AddBranch("AliAODCaloTrigger",&fPHOSTrigger);

  fEMCalCells = new AliAODCaloCells;
  fEMCalCells->SetName("emcalCells");
  fExtAOD->AddBranch("AliAODCaloCells",&fEMCalCells);
  fPHOSCells = new AliAODCaloCells;
  fPHOSCells->SetName("phosCells");
  
  fExtAOD->AddBranch("AliAODCaloCells",&fPHOSCells);
  fAODZDC = new AliAODZDC;
  fExtAOD->AddBranch("AliAODZDC",&fAODZDC);
  fAODAD = new AliAODAD;
  fExtAOD->AddBranch("AliAODAD",&fAODAD);
  fAODTZERO = new AliAODTZERO;
  fExtAOD->AddBranch("AliAODTZERO",&fAODTZERO);
  
  // Create branch for user-defined objects
  fPairs = new TClonesArray("AliDielectronPair",10);
  fPairs->SetName("dielectrons");
  fExtAOD->AddBranch("TClonesArray", &fPairs);
  fJets = new TClonesArray("AliEmcalJet", 100);
  fJets->SetName("jets");
  fExtAOD->AddBranch("TClonesArray", &fJets);

  fExtAOD->GetAOD()->GetStdContent();

  PostData(1, const_cast<THashList*>(fDielectron->GetHistogramList()));
  PostData(2,fEventStat);

  // Init array for pair daughter
  fDaughters = new TClonesArray("TLorentzVector",20);
  fDaughters->SetName("daughters");
}


void AliAnalysisTaskJpsiJetFilter::Init(){
  // Initialization
  if (fDebug > 1) AliInfo("Init() \n");
  // require AOD handler
  AliAODHandler *aodH = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
  if (!aodH) AliFatal("No AOD handler. Halting.");
  fExtAOD = aodH->AddFilteredAOD(fOutputFileName, "DielectronEvents", fIsToMerge);
  if(!fExtAOD) AliFatal("Fail to add filtered AOD");
  // Dielectron
  if (!fDielectron) InitDielectron();
}

// Update: copy title for primary vertex
void AliAnalysisTaskJpsiJetFilter::UserExec(Option_t*){
  //
  // Main loop. Called for every event
  //

  if (!fDielectron)
    return;

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD = man->GetInputEventHandler()->IsA() == AliESDInputHandler::Class();
  Bool_t isAOD = man->GetInputEventHandler()->IsA() == AliAODInputHandler::Class();

  AliInputEventHandler *inputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
  if (!inputHandler)
    return;

  if (inputHandler->GetPIDResponse())
  {
    AliDielectronVarManager::SetPIDResponse(inputHandler->GetPIDResponse());
  }
  else
  {
    AliFatal("This task needs the PID response attached to the input event handler!");
  }

  // Was event selected ?
  ULong64_t isSelected = AliVEvent::kAny;
  Bool_t isRejected = kFALSE;
  if (fSelectPhysics && inputHandler)
  {
    if ((isESD && inputHandler->GetEventSelection()) || isAOD)
    {
      isSelected = inputHandler->IsEventSelected();
      if (fExcludeTriggerMask && (isSelected & fExcludeTriggerMask))
        isRejected = kTRUE;
      if (fTriggerLogic == kAny)
        isSelected &= fTriggerMask;
      else if (fTriggerLogic == kExact)
        isSelected = ((isSelected & fTriggerMask) == fTriggerMask);
    }
  }

  //before physics selection
  fEventStat->Fill(kAllEvents);
  if (isSelected == 0 || isRejected)
  {
    PostData(2, fEventStat);
    return;
  }
  //after physics selection
  fEventStat->Fill(kSelectedEvents);

  //V0and
  if (fTriggerOnV0AND)
  {
    if (isESD)
    {
      if (!fTriggerAnalysis->IsOfflineTriggerFired(static_cast<AliESDEvent *>(InputEvent()), AliTriggerAnalysis::kV0AND))
        return;
    }
    if (isAOD)
    {
      if (!((static_cast<AliAODEvent *>(InputEvent()))->GetVZEROData()->GetV0ADecision() == AliVVZERO::kV0BB &&
            (static_cast<AliAODEvent *>(InputEvent()))->GetVZEROData()->GetV0CDecision() == AliVVZERO::kV0BB))
        return;
    }
  }

  fEventStat->Fill(kV0andEvents);

  //Fill Event histograms before the event filter
  AliDielectronHistos *h = fDielectron->GetHistoManager();

  Double_t values[AliDielectronVarManager::kNMaxValues] = {0};
  if (h)
    AliDielectronVarManager::SetFillMap(h->GetUsedVars());
  else
    AliDielectronVarManager::SetFillMap(0x0);
  AliDielectronVarManager::SetEvent(InputEvent());
  AliDielectronVarManager::Fill(InputEvent(), values);

  if (h && h->GetHistogramList()->FindObject("Event_noCuts")){
    h->FillClass("Event_noCuts", AliDielectronVarManager::kNMaxValues, values);
  }

  //event filter
  if (fEventFilter)
  {
    if (!fEventFilter->IsSelected(InputEvent()))
      return;
  }
  fEventStat->Fill(kFilteredEvents);

  //pileup
  if (fRejectPileup)
  {
    if (InputEvent()->IsPileupFromSPD(3, 0.8, 3., 2., 5.))
      return;
  }
  fEventStat->Fill(kPileupEvents);

  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField(bz);

  AliDielectronPID::SetCorrVal((Double_t)InputEvent()->GetRunNumber());

  fDielectron->Process(InputEvent());

  Bool_t hasCand = kFALSE;
  if (fStoreLikeSign)
    hasCand = (fDielectron->HasCandidates() || fDielectron->HasCandidatesLikeSign());
  else
    hasCand = (fDielectron->HasCandidates());

  if (fStoreRotatedPairs)
    hasCand = (hasCand || fDielectron->HasCandidatesTR());

  if (fStoreEventsWithSingleTracks)
    hasCand = (hasCand || fDielectron->GetTrackArray(0) || fDielectron->GetTrackArray(1));

  // Fill nano AOD
  if ( hasCand)
  {
    AliAODEvent *aodEv = (static_cast<AliAODEvent *>(InputEvent()));
    // Fill ZDC, AD, TZERO data
    AliAODZDC* zdc = aodEv->GetZDCData();
    fAODZDC->Clear("C");
    *fAODZDC = *zdc;
    AliAODAD* ad = aodEv->GetADData();
    fAODAD->Clear("C");
    *fAODAD = *ad;
    AliAODTZERO* tzero = aodEv->GetTZEROData();
    fAODTZERO->Clear("C");
    *fAODTZERO = *tzero;
    // Fil tracklets, EMCal/PHOS calo trigger
    AliAODTracklets* spd = aodEv->GetMultiplicity();
    fSPD->Clear("C");
    *fSPD = *spd;
    AliAODCaloTrigger* caloTrig = aodEv->GetCaloTrigger("EMCAL");
    fEMCALTrigger->Clear("C");
    *fEMCALTrigger = *caloTrig;
    caloTrig = aodEv->GetCaloTrigger("PHOS");
    fPHOSTrigger->Clear("C");
    *fPHOSTrigger = *caloTrig;
    // Fill EMCal/PHOS cells
    AliAODCaloCells* cells = aodEv->GetEMCALCells();
    fEMCalCells->Clear("C");
    *fEMCalCells = *cells;
    cells = aodEv->GetPHOSCells();
    fPHOSCells->Clear("C");
    *fPHOSCells = *cells;
    // Clear arrays
    AliAODEvent *nanoEv = fExtAOD->GetAOD();
    nanoEv->GetTracks()->Clear("C");
    nanoEv->GetVertices()->Clear("C");
    nanoEv->GetCaloClusters()->Clear("C");
    // Fill primary and SPD vertex
    AliAODVertex *vtxPriv = aodEv->GetPrimaryVertex();
    if(!vtxPriv) AliFatal("No primary vertex");
    AliAODVertex *vtxSPD = aodEv->GetPrimaryVertexSPD();
    AliAODVertex *vtxTPC = aodEv->GetPrimaryVertexTPC();
    AliAODVertex *tmp = vtxPriv->CloneWithoutRefs();
    tmp->SetTitle("VertexerTracksMVWithConstraint");
    nanoEv->AddVertex(tmp);
    AliAODVertex *tmpSPD = vtxSPD->CloneWithoutRefs(); 
    tmpSPD->SetTitle(vtxSPD->GetTitle());
    nanoEv->AddVertex(tmpSPD); 
    AliAODVertex *tmpTPC = vtxTPC->CloneWithoutRefs(); 
    tmpTPC->SetTitle(vtxTPC->GetTitle());
    nanoEv->AddVertex(tmpTPC); 
    
    // Fill pairs
    Int_t ncandidates = fDielectron->GetPairArray(1)->GetEntriesFast();
    if (ncandidates == 1)
      fEventStat->Fill((kNbinsEvent));
    else if (ncandidates > 1)
      fEventStat->Fill((kNbinsEvent + 1));
    fPairs->Clear("C");
    fDaughters->Clear("C");
    Int_t nD = 0;
    const TObjArray* candidates = fDielectron->GetPairArray(1);
    for(Int_t i = 0; i < ncandidates; i++){
      AliDielectronPair* pair = (AliDielectronPair*)(candidates->UncheckedAt(i));
      new((*fPairs)[i]) AliDielectronPair(*pair);
      const AliKFParticle& d1 = pair->GetKFFirstDaughter();
      new ((*fDaughters)[nD++]) TLorentzVector(d1.GetPx(),d1.GetPy(),d1.GetPz(),d1.GetE());
      const AliKFParticle& d2 = pair->GetKFSecondDaughter();
      new ((*fDaughters)[nD++]) TLorentzVector(d2.GetPx(),d2.GetPy(),d2.GetPz(),d2.GetE());
    }
    fPairs->Expand(ncandidates);

    // Fill tracks
    Int_t nTracks = aodEv->GetNumberOfTracks();
    Int_t nTrackMatched = 0;
    AliAODTrack* trkTemplate = NULL; // To insert pair as track
    for (int iTrk = 0; iTrk < nTracks; iTrk++)
    {
      AliAODTrack* oldTrack = (AliAODTrack*)(aodEv->GetTrack(iTrk));
      // To remove tracks used as pair daughters
      if(fIsToReplace && FindDaughters(oldTrack)){
        nTrackMatched++;
        if(!trkTemplate || trkTemplate->Pt() < oldTrack->Pt()) trkTemplate = oldTrack;
        continue;
      }
      // Add track for nano AOD
      Int_t trkID = nanoEv->AddTrack(oldTrack);
      AliAODTrack* trk = (AliAODTrack*)(nanoEv->GetTrack(trkID));
      trk->ResetBit(kIsReferenced);
      trk->SetUniqueID(0);
      // Vertex of origin
      Int_t vtxID = nanoEv->AddVertex(oldTrack->GetProdVertex());
      AliAODVertex* vtx = nanoEv->GetVertex(vtxID);
      vtx->ResetBit(kIsReferenced);
      vtx->SetUniqueID(0);
      vtx->RemoveDaughters();
      trk->SetProdVertex(vtx);
      // Calo cluster
      Int_t oldCaloID = oldTrack->GetEMCALcluster();
      if(oldCaloID >= 0){
        AliAODCaloCluster* oldCalo = aodEv->GetCaloCluster(oldCaloID);
        if (oldCalo)
        {
          Int_t caloID = nanoEv->AddCaloCluster(oldCalo);
          trk->SetEMCALcluster(caloID);
          AliAODCaloCluster *calo = nanoEv->GetCaloCluster(caloID);
          for (int u = 0; u < calo->GetNTracksMatched(); u++)
            calo->RemoveTrackMatched(calo->GetTrackMatched(u));
          calo->AddTrackMatched(trk);
        }
      }
    }
    if(fIsToReplace){
      TIter nextPair(fPairs);
      AliDielectronPair* pair = NULL;
      while(pair = static_cast<AliDielectronPair*>(nextPair())){
        Int_t trkID = nanoEv->AddTrack(trkTemplate);
        AliAODTrack* trk = (AliAODTrack*)(nanoEv->GetTrack(trkID));
        SetTrackFromPair(pair, trk);
        trk->SetProdVertex(nanoEv->GetPrimaryVertex());
      }
      AliDebug(2, Form("Remove Ndaughters : %d, Add Npairs : %d", nTrackMatched, ncandidates));
    }
    nanoEv->GetTracks()->Expand(nanoEv->GetNumberOfTracks());
    nanoEv->GetVertices()->Expand(nanoEv->GetNumberOfVertices());
    nanoEv->GetCaloClusters()->Expand(nanoEv->GetNumberOfCaloClusters());


    // Fill jets
    FillJets(aodEv, fJets, "Jet_AKTChargedR040_tracks_pT0150_pt_scheme");

    // Write output
    fExtAOD->SelectEvent();
    // DEBUG - not use FinishEvent() to avoid auto flush
    fExtAOD->GetTree()->Fill();
    
    delete tmp;
    delete tmpSPD;
    delete tmpTPC;
  }

  PostData(1, const_cast<THashList *>(fDielectron->GetHistogramList()));
  PostData(2, fEventStat);
  return;
} 

void AliAnalysisTaskJpsiJetFilter::FillJets(AliAODEvent* aodEv, TClonesArray* jetArray, TString jetName){
  jetArray->Clear("C");
  TClonesArray* jets = dynamic_cast<TClonesArray*>(aodEv->FindListObject(jetName));
  if(jets){
    for(Int_t i = 0; i < jets->GetEntriesFast(); i++){
      AliEmcalJet* jet = dynamic_cast<AliEmcalJet*>(jets->UncheckedAt(i));
      new((*jetArray)[i]) AliEmcalJet(*jet);
    }
    jetArray->Expand(jets->GetEntriesFast());
  }else{
    AliInfo(Form("Could not found jet container : %s", jetName.Data()));
  }

}

Bool_t AliAnalysisTaskJpsiJetFilter::FindDaughters(AliVTrack* trk){
  static const Double_t ERR_LIMIT = 1e-6;
  TIter nextEle(fDaughters);
  TLorentzVector* vec = NULL;
  while(vec = static_cast<TLorentzVector*>(nextEle())){
    if(TMath::Abs(vec->Pt() - trk->Pt()) < ERR_LIMIT &&
       TMath::Abs(vec->Eta() - trk->Eta()) < ERR_LIMIT &&
       TMath::Abs(TVector2::Phi_0_2pi(vec->Phi()) - trk->Phi()) < ERR_LIMIT )
      return kTRUE;
  }
  return kFALSE;
}

void AliAnalysisTaskJpsiJetFilter::SetTrackFromPair(AliDielectronPair* pair, AliAODTrack* trk){
  
  //trk->SetStatus(AliVTrack::kEmbedded);

  trk->SetPt(pair->Pt());
  trk->SetPhi(TVector2::Phi_0_2pi(pair->Phi()));
  trk->SetTheta(TMath::ACos(pair->Pz() / pair->P()));

  // Remove EMCal
  trk->ResetStatus(AliVTrack::kEMCALmatch);
  trk->SetEMCALcluster(AliVTrack::kEMCALNoMatch);

  // Reset reference
  trk->ResetBit(kIsReferenced);
  trk->SetUniqueID(0);

  // DEBUG - Pseudo-proper decay length
  auto aod = (AliAODEvent*)(InputEvent());
  auto priv = aod->GetPrimaryVertex();
  Double_t errPseudoProperTime2 = 0.;
  AliKFParticle kfPair = pair->GetKFParticle();
  Double_t lxy = kfPair.GetPseudoProperDecayTime(*priv, TDatabasePDG::Instance()->GetParticle(443)->Mass(), &errPseudoProperTime2 );
  trk->SetTrackPhiEtaPtOnEMCal(pair->M(), lxy, 0.);

}

void AliAnalysisTaskJpsiJetFilter::InitDielectron(){

  fDielectron = new AliDielectron("Diele","Dielectron with EMCal triggered");

/**
 *  Track cuts
**/
  // Good track
  AliDielectronTrackCuts *trackCuts = new AliDielectronTrackCuts("trackCuts", "trackCuts");
  trackCuts->SetClusterRequirementITS(AliDielectronTrackCuts::kSPD, AliDielectronTrackCuts::kAny);
  trackCuts->SetRequireTPCRefit(kTRUE);
  trackCuts->SetRequireITSRefit(kTRUE);
  fDielectron->GetTrackFilter().AddCuts(trackCuts);
  // Track cuts for electron PID
  AliDielectronVarCuts *ePID = new AliDielectronVarCuts("ePidCut", "Track cuts for electron PID");
  ePID->AddCut(AliDielectronVarManager::kKinkIndex0, 0.);
  // Track quality cuts
  ePID->AddCut(AliDielectronVarManager::kNclsTPC, 70., 160.); // Loose
  ePID->AddCut(AliDielectronVarManager::kEta, -0.9, 0.9);
  ePID->AddCut(AliDielectronVarManager::kImpactParXY, -1., 1.);
  ePID->AddCut(AliDielectronVarManager::kImpactParZ, -3., 3.);
  ePID->AddCut(AliDielectronVarManager::kITSLayerFirstCls, 0., 4.);
  ePID->AddCut(AliDielectronVarManager::kTPCchi2Cl, 0., 4.);
  ePID->AddCut(AliDielectronVarManager::kPt, 0.7, 1e30); // Loose
  // Electron PID with TPC
  ePID->AddCut(AliDielectronVarManager::kTPCnSigmaEle, -3.0, 3.0); // Loose
  fDielectron->GetTrackFilter().AddCuts(ePID);

/**
 *  Pair cuts
**/
  // EMCal gamma/electron energy threshold
  Double_t EMCal_E_Threshold = 5.0; // GeV, ADC setting

  //Invariant mass and rapidity selection
  AliDielectronVarCuts *pairCut = new AliDielectronVarCuts("jpsiCuts", "1<M<5 + |Y|<.9");
  pairCut->AddCut(AliDielectronVarManager::kM, 1.0, 5.0);
  pairCut->AddCut(AliDielectronVarManager::kY, -0.9, 0.9);
  pairCut->AddCut(AliDielectronVarManager::kPt, 1, 1e30);
  pairCut->AddCut(AliDielectronVarManager::kPt, EMCal_E_Threshold, 1e30);
  fDielectron->GetPairFilter().AddCuts(pairCut);
  // Leg cuts
  AliDielectronVarCuts *emcCut = new AliDielectronVarCuts("CutEMCAL", "Jpsi leg cuts for EMCal");
  emcCut->AddCut(AliDielectronVarManager::kEMCALEoverP, 0.75, 1.35); // Loose
  emcCut->AddCut(AliDielectronVarManager::kEMCALE, EMCal_E_Threshold, 1e30);
  //emcCut->AddCut(AliDielectronVarManager::kPhi, 4.377, 5.7071, kTRUE); // Exclude DCal
  //emcCut->AddCut(AliDielectronVarManager::kPhi, 1.396, 3.2637, kTRUE); // Exclude EMCal
  AliDielectronPairLegCuts *legVar = new AliDielectronPairLegCuts();
  legVar->GetLeg1Filter().AddCuts(emcCut);
  legVar->GetLeg2Filter().AddCuts(emcCut);
  legVar->SetCutType(AliDielectronPairLegCuts::kAnyLeg);
  fDielectron->GetPairFilter().AddCuts(legVar);
}

void AliAnalysisTaskJpsiJetFilter::InitHistogramsForDielectron(const char* histMgrName)
{
  //Setup histogram Manager
  AliDielectronHistos *histos = new AliDielectronHistos(histMgrName, "Histograms for dielectron");
  fDielectron->SetHistogramManager(histos);

  //Initialise histogram classes
  histos->SetReservedWords("Track;Pair");
  //Track classes
  for (Int_t i = 0; i < 2; ++i)
  {
    histos->AddClass(Form("Track_%s", AliDielectron::TrackClassName(i)));
  }
  //Pair classes
  for (Int_t i = 0; i < 3; ++i)
    histos->AddClass(Form("Pair_%s", AliDielectron::PairClassName(i)));
  //Legs from pair
  for (Int_t i = 0; i < 3; ++i)
    histos->AddClass(Form("Track_Legs_%s", AliDielectron::PairClassName(i)));

/**
 *  Event
**/
  TString histClass = "Event";
  histos->AddClass(histClass);
  // Event Primary vertex and diamond (IP) stats.
  histos->UserHistogram(histClass, "VtxZ", "Vertex Z;Z[cm];#events", 1000, -50., 50., AliDielectronVarManager::kZvPrim);
  histos->UserHistogram(histClass, "VtxX", "Vertex X;X[cm];#events", 1000, -1.0, 1.0, AliDielectronVarManager::kXvPrim);
  histos->UserHistogram(histClass, "VtxY", "Vertex Y;Y[cm];#events", 2000, -1.0, 1.0, AliDielectronVarManager::kYvPrim);
  histos->UserHistogram(histClass, "NVContrib", "Number of Vertex contributors", 1001, -0.5, 1000.5, AliDielectronVarManager::kNVtxContrib);
  // Event track and SPD (tracklets) stats.
  histos->UserHistogram(histClass, "kNTrk", "Number of tracks;kNTrk;Entries", 4001, -0.5, 4000.5, AliDielectronVarManager::kNTrk);
  histos->UserHistogram(histClass, "kNaccTrcklts", "Number of accepted SPD tracklets in |eta|<1.6;kNaccTrcklts;Entries", 1001, -0.5, 1000.5, AliDielectronVarManager::kNaccTrcklts);
  histos->UserHistogram(histClass, "kNaccTrcklts10Corr", "kNaccTrcklts10Corr;kNaccTrcklts10Corr;Entries", 501, -0.5, 500.5, AliDielectronVarManager::kNaccTrcklts10Corr);
  histos->UserHistogram(histClass, "VtxZ_kNaccTrcklts10Corr", "VtxZ vs. kNaccTrcklts10Corr;VtxZ;kNaccTrcklts10Corr", 800, -40., 40., 501, -0.5, 500.5, AliDielectronVarManager::kZvPrim, AliDielectronVarManager::kNaccTrcklts10Corr);
  //new multiplicity estimator: V0
  histos->UserHistogram(histClass, "kMultV0", "kMultV0;kMultV0;Entries", 1000, 0., 1000., AliDielectronVarManager::kMultV0);
  histos->UserHistogram(histClass, "kMultV0A", "kMultV0;kMultV0;Entries", 1000, 0., 1000., AliDielectronVarManager::kMultV0A);
  histos->UserHistogram(histClass, "kMultV0C", "kMultV0;kMultV0;Entries", 1000, 0., 1000., AliDielectronVarManager::kMultV0C);
  // 2D
  histos->UserHistogram(histClass, "kMultV0A_kMultV0C", "kMultV0A vs. kMultV0C;kMultV0A;kMultV0C", 1000, 0., 1000., 1000, 0., 1000., AliDielectronVarManager::kMultV0A, AliDielectronVarManager::kMultV0C);
  histos->UserHistogram(histClass, "kMultV0_kMultV0A", "kMultV0 vs. kMultV0A;kMultV0;kMultV0A", 1000, 0., 1000., 1000, 0., 1000., AliDielectronVarManager::kMultV0, AliDielectronVarManager::kMultV0A);
  histos->UserHistogram(histClass, "kMultV0_kMultV0C", "kMultV0 vs. kMultV0C;kMultV0;kMultV0C", 1000, 0., 1000., 1000, 0., 1000., AliDielectronVarManager::kMultV0, AliDielectronVarManager::kMultV0C);
  // vs Vertex Z
  histos->UserHistogram(histClass, "VtxZ_kMultV0A", "VtxZ vs. kMultV0A;VtxZ;kMultV0A", 300, -15., 15., 1000, 0., 1000., AliDielectronVarManager::kZvPrim, AliDielectronVarManager::kMultV0A);
  histos->UserHistogram(histClass, "VtxZ_kMultV0C", "VtxZ vs. kMultV0C;VtxZ;kMultV0C", 300, -15., 15., 1000, 0., 1000., AliDielectronVarManager::kZvPrim, AliDielectronVarManager::kMultV0C);
  histos->UserHistogram(histClass, "VtxZ_kMultV0", "VtxZ vs. kMultV0;VtxZ;kMultV0", 300, -15., 15., 1000, 0., 1000., AliDielectronVarManager::kZvPrim, AliDielectronVarManager::kMultV0);

  // Dielectron info.
  histos->UserHistogram(histClass, "Nelectrons", "Number of tracks/electron selected by AliDielectron after cuts;N_{e};#events", 50, -0.5, 49.5, AliDielectronVarManager::kTracks);
  histos->UserHistogram(histClass, "Npairs", "Number of Ev1PM pair candidates after all cuts;J/#psi candidates;#events", 20, -0.5, 19.5, AliDielectronVarManager::kPairs);
  /*
	  Histogram for Track
	*/
  // Track kinetics parameter
  histos->UserHistogram("Track", "Pt", "Pt;Pt [GeV/c];#tracks", 2000, 0, 100, AliDielectronVarManager::kPt, kTRUE);
  histos->UserHistogram("Track", "Eta_Phi", "Eta Phi Map; Eta; Phi;#tracks",
                        100, -1, 1, 144, 0, TMath::TwoPi(), AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi, kTRUE);
  histos->UserHistogram("Track", "dXY", "dXY;dXY [cm];#tracks", 1000, -50, 50, AliDielectronVarManager::kImpactParXY, kTRUE);
  histos->UserHistogram("Track", "dZ", "dZ;dZ [cm];#tracks", 1000, -50., 50., AliDielectronVarManager::kImpactParZ, kTRUE);
  // Tracking quality
  histos->UserHistogram("Track", "ITS_FirstCls", "ITS First Cluster;Layer No. of ITS 1st cluster;#Entries", 6, 0.5, 6.5, AliDielectronVarManager::kITSLayerFirstCls, kTRUE);
  histos->UserHistogram("Track", "TPCnCls", "Number of Clusters TPC;TPC number clusteres;#tracks", 161, -0.5, 160.5, AliDielectronVarManager::kNclsTPC, kTRUE);
  histos->UserHistogram("Track", "TPCchi2Cl", "Chi-2/Clusters TPC;Chi2/ncls number clusteres;#tracks", 100, 0, 10, AliDielectronVarManager::kTPCchi2Cl, kTRUE);
  // PID - TPC
  histos->UserHistogram("Track", "dEdx_P", "dEdx vs PinTPC;P [GeV];TPC signal (a.u.);#tracks",
                        800, 0., 40., 800, 20., 200., AliDielectronVarManager::kPIn, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "dEdx_Pt", "dEdx vs Pt;Pt [GeV];TPC signal (a.u.);#tracks",
                        800, 0., 40., 800, 20., 200., AliDielectronVarManager::kPt, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "dEdx_Phi", "dEdx vs #phi;#phi [rad];TPC signal (a.u.);#tracks",
                        200, 0., 2 * TMath::Pi(), 800, 20., 200., AliDielectronVarManager::kPhi, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "dEdx_Eta", "dEdx vs #eta;#eta;TPC signal (a.u.);#tracks",
                        200, -1., 1., 800, 20., 200., AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "TPCnSigmaEle_P", "n#sigma_{e}(TPC) vs P_{in} TPC;P_{in} [GeV];n#sigma_{e}(TPC);#tracks",
                        800, 0., 40., 800, -12., 12., AliDielectronVarManager::kPIn, AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "TPCnSigmaEle_Pt", "n#sigma_{e}(TPC) vs Pt;Pt [GeV];n#sigma_{e}(TPC);#tracks",
                        800, 0., 40., 800, -12., 12., AliDielectronVarManager::kPt, AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "TPCnSigmaEle_Phi", "n#sigma_{e}(TPC) vs #phi;#phi [rad];n#sigma_{e}(TPC);#tracks",
                        200, 0., 2 * TMath::Pi(), 800, -12., 12., AliDielectronVarManager::kPhi, AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "TPCnSigmaEle_Eta", "n#sigma_{e}(TPC) vs #eta;#eta;n#sigma_{e}(TPC);#tracks",
                        200, -1., 1., 800, -12., 12., AliDielectronVarManager::kEta, AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "dEdx_nSigmaEMCal", "dEdx vs n#sigma_{e}(EMCAL);n#sigma_{e}(EMCAL);TPC signal (a.u.);#tracks",
                        200, -5., 5., 800, 20., 200., AliDielectronVarManager::kEMCALnSigmaEle, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "dEdx_TPCnSigmaEle", "dEdx vs n#sigma_{e}(TPC);n#sigma_{e}(TPC);TPC signal (a.u.);#tracks",
                        100, -10., 10., 800, 20., 200., AliDielectronVarManager::kTPCnSigmaEle, AliDielectronVarManager::kTPCsignal, kTRUE);
  // Track - EMCal
  histos->UserHistogram("Track", "EMCalE", "EmcalE;Cluster Energy [GeV];#Clusters",
                        200, 0., 40., AliDielectronVarManager::kEMCALE, kTRUE);
  histos->UserHistogram("Track", "EMCalE_P", "Cluster energy vs. pT; EMCal_E;pT;#tracks",
                        800, 0., 40, 200, 0., 40.,  AliDielectronVarManager::kPIn, AliDielectronVarManager::kEMCALE, kTRUE);
  histos->UserHistogram("Track", "EMCalE_Pt", "Cluster energy vs. pT; EMCal_E;pT;#tracks",
                        800, 0., 40, 200, 0., 40.,  AliDielectronVarManager::kPt, AliDielectronVarManager::kEMCALE, kTRUE);
  //Ecluster versus Phi to separate EMCal and DCal
  histos->UserHistogram("Track", "EMCalE_Phi", "Cluster energy vs. #phi; EMCal_E;Phi;#tracks",
                        200, 0., TMath::TwoPi(), 200, 0., 40., AliDielectronVarManager::kPhi, AliDielectronVarManager::kEMCALE, kTRUE);
  histos->UserHistogram("Track", "EMCalE_Eta", "Cluster energy vs. #eta; EMCal_E;Eta;#tracks",
                        200, -1.0, 1.0, 200, 0., 40., AliDielectronVarManager::kEta, AliDielectronVarManager::kEMCALE, kTRUE);
  // PID - EMCal
  // E/p ratio
  histos->UserHistogram("Track", "EoverP", "EMCal E/p ratio;E/p;#Clusters",
                        200, 0., 2., AliDielectronVarManager::kEMCALEoverP, kTRUE);
  histos->UserHistogram("Track", "EoverP_P", "E/p ratio vs P;P_{in} (GeV/c);E/p;#tracks",
                        200, 0., 40., 200, 0., 2., AliDielectronVarManager::kPIn, AliDielectronVarManager::kEMCALEoverP, kTRUE);  
  histos->UserHistogram("Track", "EoverP_Pt", "E/p ratio vs Pt;Pt (GeV/c);E/p;#tracks",
                        200, 0., 40., 200, 0., 2., AliDielectronVarManager::kPt, AliDielectronVarManager::kEMCALEoverP, kTRUE);
  histos->UserHistogram("Track", "EoverP_Phi", "E/p ratio vs #phi;Phi;E/p;#tracks",
                        200, 0., TMath::TwoPi(), 200, 0., 2., AliDielectronVarManager::kPhi, AliDielectronVarManager::kEMCALEoverP, kTRUE);
  histos->UserHistogram("Track", "EoverP_Eta", "E/p ratio vs #eta;Eta;E/p;#tracks",
                        200, -1.0, 1.0, 200, 0., 2., AliDielectronVarManager::kEta, AliDielectronVarManager::kEMCALEoverP, kTRUE);
  // EMCal nSigma electron
  histos->UserHistogram("Track", "EMCALnSigmaE_P", "n#sigma_{e} vs Pt;Pt (GeV/c);n#sigma_{e};#tracks",
                        200, 0., 40., 200, -12, 12, AliDielectronVarManager::kPIn, AliDielectronVarManager::kEMCALnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "EMCALnSigmaE_Phi", "n#sigma_{e} vs #phi;Phi;n#sigma_{e};#tracks",
                        200, 0., TMath::TwoPi(), 200, -12, 12, AliDielectronVarManager::kPhi, AliDielectronVarManager::kEMCALnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "EMCALnSigmaE_Eta", "n#sigma_{e} vs #eta;Eta;n#sigma_{e};#tracks",
                        200, -1.0, 1.0, 200, 0., 2., AliDielectronVarManager::kEta, AliDielectronVarManager::kEMCALnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "EMCALnSigmaEle_EoverP", "n#sigma_{e}(EMCal) vs E/p;E/p;n#sigma_{e}(EMCal);#tracks",
                        200, 0., 2., 200, -12., 12., AliDielectronVarManager::kEMCALEoverP, AliDielectronVarManager::kEMCALnSigmaEle, kTRUE);
  // PID - TPC + EMCal
  histos->UserHistogram("Track", "dEdx_EoverP", "dEdx vs E/p;E/P;TPC signal (a.u.);#tracks",
                        200, 0., 2., 800, 20., 200., AliDielectronVarManager::kEMCALEoverP, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "TPCnSigmaEle_EoverP", "n#sigma_{e}(TPC) vs E/p;E/p;n#sigma_{e}(TPC);#tracks",
                        200, 0., 2., 200, -12., 12., AliDielectronVarManager::kEMCALEoverP, AliDielectronVarManager::kTPCnSigmaEle, kTRUE);
  histos->UserHistogram("Track", "dEdx_EMCALnSigmaE", "dEdx vs n#sigma_{e}(EMCAL);n#sigma_{e}(EMCAL);TPC signal (a.u.);#tracks",
                        200, -12., 12., 800, 20., 200., AliDielectronVarManager::kEMCALnSigmaEle, AliDielectronVarManager::kTPCsignal, kTRUE);
  histos->UserHistogram("Track", "nSigmaTPC_EMCal", "n#sigma_{e}(TPC vs EMCAL);n#sigma_{e}(EMCAL);n#sigma_{e}(TPC);#tracks",
                        200, -5., 5., 200, -12., 12., AliDielectronVarManager::kEMCALnSigmaEle, AliDielectronVarManager::kTPCnSigmaEle, kTRUE);

/**
 *  Histograms for Pair
**/
  histos->UserHistogram("Pair", "InvMass", "Inv.Mass;Inv. Mass (GeV/c^{2});#pairs/(40 MeV/c^{2})",
                        100, 1.0, 5.0, AliDielectronVarManager::kM);
  histos->UserHistogram("Pair", "pT", "Pt;Pt (GeV/c);#pairs",
                        2000, 0., 100.0, AliDielectronVarManager::kPt);
  histos->UserHistogram("Pair", "Eta_Phi", "#eta-#phi map of dielectron pairs;#eta_{ee};#phi_{ee};#pairs",
                        200, -1, 1, 200, 0., 10, AliDielectronVarManager::kEta, AliDielectronVarManager::kPhi);
  histos->UserHistogram("Pair", "Rapidity", "Rapidity;Rapidity;#pairs",
                        200, -1., 1., AliDielectronVarManager::kY);
  histos->UserHistogram("Pair", "OpeningAngle", "Opening angle / rad;#pairs",
                        50, 0., 3.15, AliDielectronVarManager::kOpeningAngle);

  histos->UserHistogram("Pair", "PseudoProperTime", "Pseudoproper decay length; pseudoproper-decay-length[cm];#pairs / 40#mum",
                        600, -0.3, 0.3, AliDielectronVarManager::kPseudoProperTime);
  histos->UserHistogram("Pair", "InvMass_Pt", "Inv. Mass vs Pt;Pt (GeV/c); Inv. Mass (GeV/c^{2})",
                        200, 0., 40., 100, 1.0, 5.0, AliDielectronVarManager::kPt, AliDielectronVarManager::kM);
  histos->UserHistogram("Pair", "OpeningAngle_Pt", "Opening angle vs p_{T} ;p_{T} (GeV/c); angle",
                        200, 0., 40., 200, 0, TMath::Pi(), AliDielectronVarManager::kPt, AliDielectronVarManager::kOpeningAngle);
  //InvMass versus Proper time
  histos->UserHistogram("Pair", "InvMass_ProperTime", "InvMass vs. ProperTime;pseudoproper-decay-length[cm]; Inv. Mass [GeV]",
                        600, -0.3, 0.3, 100, 1.0, 5.0, AliDielectronVarManager::kPseudoProperTime, AliDielectronVarManager::kM);
}
