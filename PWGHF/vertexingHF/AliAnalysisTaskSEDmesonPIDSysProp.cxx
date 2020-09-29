/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. */

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// \class AliAnalysisTaskSEDmesonPIDSysProp                                                              //
// \brief analysis task for PID Systematic uncertainty propagation from the single track to the D mesons //
// \author: A. M. Barbano, anastasia.maria.barbano@cern.ch                                               //
// \author: F. Grosa, fabrizio.grosa@cern.ch                                                             //
///////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTaskSEDmesonPIDSysProp.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliPIDResponse.h"
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>

ClassImp(AliAnalysisTaskSEDmesonPIDSysProp)

//________________________________________________________________________
AliAnalysisTaskSEDmesonPIDSysProp::AliAnalysisTaskSEDmesonPIDSysProp():
AliAnalysisTaskSE("taskPIDSysProp"),
fOutput(nullptr),
fHistNEvents(nullptr),
fHistPtDauVsD(nullptr),
fHistSystPIDEffD(nullptr),
fHistEffPionTOF(nullptr),
fHistEffKaonTOF(nullptr),
fHistSystPionTOF(nullptr),
fHistSystKaonTOF(nullptr),
fPartName(""),
fPIDresp(nullptr),
fPIDstrategy(kConservativePID),
fnSigma(3.),
fDecayChannel(kD0toKpi),
fKaonTPCHistoOpt(kKaonTOFtag),
fKaonTOFHistoOpt(kSamePionV0tag),
fAODProtection(1),
fNPtBins(0),
fPtLimits(nullptr),
fAnalysisCuts(nullptr),
fVarForProp(kPt)
{
  for(int iHist=0; iHist<2; iHist++) {
    fHistEffPionTPC[iHist]=nullptr;
    fHistEffKaonTPC[iHist]=nullptr;
    fHistSystPionTPC[iHist]=nullptr;
    fHistSystKaonTPC[iHist]=nullptr;
  }

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskSEDmesonPIDSysProp::AliAnalysisTaskSEDmesonPIDSysProp(int ch, AliRDHFCuts* cuts):
AliAnalysisTaskSE("taskPIDSysProp"),
fOutput(nullptr),
fHistNEvents(nullptr),
fHistPtDauVsD(nullptr),
fHistSystPIDEffD(nullptr),
fHistEffPionTOF(nullptr),
fHistEffKaonTOF(nullptr),
fHistSystPionTOF(nullptr),
fHistSystKaonTOF(nullptr),
fPartName(""),
fPIDresp(nullptr),
fPIDstrategy(kConservativePID),
fnSigma(3.),
fDecayChannel(ch),
fKaonTPCHistoOpt(kKaonTOFtag),
fKaonTOFHistoOpt(kSamePionV0tag),
fAODProtection(1),
fNPtBins(0),
fPtLimits(nullptr),
fAnalysisCuts(cuts),
fVarForProp(kPt)
{
  for(int iHist=0; iHist<2; iHist++) {
    fHistEffPionTPC[iHist]=nullptr;
    fHistEffKaonTPC[iHist]=nullptr;
    fHistSystPionTPC[iHist]=nullptr;
    fHistSystKaonTPC[iHist]=nullptr;
  }
    
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//___________________________________________________________________________
AliAnalysisTaskSEDmesonPIDSysProp::~AliAnalysisTaskSEDmesonPIDSysProp()
{  
  if(fOutput && !fOutput->IsOwner()){
    delete fHistNEvents;
    for(int iHist=0; iHist<2; iHist++) {
      if(fHistEffPionTPC[iHist]) delete fHistEffPionTPC[iHist];
      if(fHistEffKaonTPC[iHist]) delete fHistEffKaonTPC[iHist];
      if(fHistSystPionTPC[iHist]) delete fHistSystPionTPC[iHist];
      if(fHistSystKaonTPC[iHist]) delete fHistSystKaonTPC[iHist];
    }
    delete fHistSystPionTOF;
    delete fHistSystKaonTOF;
    delete fHistEffPionTOF;
    delete fHistEffKaonTOF;
    
    delete fHistPtDauVsD;
    delete fHistSystPIDEffD;
  }
  if(fPIDresp) delete fPIDresp;
  if(fAnalysisCuts) delete fAnalysisCuts;
  if(fPtLimits) delete[] fPtLimits; 
  delete fOutput;
}

//________________________________________________________________________
void AliAnalysisTaskSEDmesonPIDSysProp::UserCreateOutputObjects()
{  
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("OutputHistos");
  
  fHistNEvents = new TH1F("hNEvents", "number of events ",15,-0.5,13.5);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"nEventsRead");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"nEvents Matched dAOD");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"nEvents Mismatched dAOD");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"nEventsAnal");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"n. passing IsEvSelected");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"n. rejected due to trigger");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"n. rejected due to not reco vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"n. rejected for contr vertex");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"n. rejected for vertex out of accept");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"n. rejected for pileup events");
  fHistNEvents->GetXaxis()->SetBinLabel(11,"no. of out centrality events");
  fHistNEvents->GetXaxis()->SetBinLabel(12,"no. of D candidates");
  fHistNEvents->GetXaxis()->SetBinLabel(13,"no. of D after PID cuts");
  fHistNEvents->GetXaxis()->SetBinLabel(14,"no. of not on-the-fly rec D");
  
  fHistNEvents->GetXaxis()->SetNdivisions(1,false);
  
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fOutput->Add(fHistNEvents);
  
  fNPtBins = fAnalysisCuts->GetNPtBins();
  float* ptlims = fAnalysisCuts->GetPtBinLimits();
  fPtLimits = new double[fNPtBins+1];
  for(int iPt=0; iPt<=fNPtBins; iPt++) {
    fPtLimits[iPt] = static_cast<double>(ptlims[iPt]);
  }
  if(fPtLimits[fNPtBins]>100) fPtLimits[fNPtBins] = 100.;
  
  TString varname = "";
  if(fVarForProp==kPt) 
    varname = "#it{p}_{T}";
  else if(fVarForProp==kP) 
    varname = "#it{p}";

  fHistSystPIDEffD = new TH2F("fHistSystPIDEffD","PID efficiency systematic uncertainty; #it{p}_{T}^{D} (GeV/#it{c}); relative systematic uncertainty",fNPtBins,fPtLimits,500,0.,0.5);
  fHistPtDauVsD = new TH2F("fHistPtDauVsD",Form("%s Dau vs #it{p}_{T} D; #it{p}_{T}^{D} (GeV/#it{c}); %s^{daugh} (GeV/#it{c})",varname.Data(),varname.Data()),static_cast<int>(fPtLimits[fNPtBins] * 10),0.,fPtLimits[fNPtBins],static_cast<int>(fPtLimits[fNPtBins] * 10),0.,fPtLimits[fNPtBins]);
  fOutput->Add(fHistSystPIDEffD);
  fOutput->Add(fHistPtDauVsD);
    
  PostData(1,fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskSEDmesonPIDSysProp::Init()
{    
  if(fDecayChannel == kDplustoKpipi) {
    fAnalysisCuts = new AliRDHFCutsDplustoKpipi(*(dynamic_cast<AliRDHFCutsDplustoKpipi*>(fAnalysisCuts)));
  }
  else if(fDecayChannel == kDstoKKpi) {
    fAnalysisCuts = new AliRDHFCutsDstoKKpi(*(dynamic_cast<AliRDHFCutsDstoKKpi*>(fAnalysisCuts)));
  }
  else if(fDecayChannel == kD0toKpi) {
    fAnalysisCuts = new AliRDHFCutsD0toKpi(*(dynamic_cast<AliRDHFCutsD0toKpi*>(fAnalysisCuts)));
  }
  else if(fDecayChannel == kDstartoKpipi) {
    fAnalysisCuts = new AliRDHFCutsDStartoKpipi(*(dynamic_cast<AliRDHFCutsDStartoKpipi*>(fAnalysisCuts)));
  }
  else {
    AliFatal("The decay channel MUST be defined according to AliCFVertexing::DecayChannel - Exit.");
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEDmesonPIDSysProp::UserExec(Option_t *)
{    
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(fInputEvent);
  fHistNEvents->Fill(0); // all events
  
  if(fAODProtection>=0){
    //   Protection against different number of events in the AOD and deltaAOD
    //   In case of discrepancy the event is rejected.
    int matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
    if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)) {
      // AOD/deltaAOD trees have different number of entries || TProcessID do not match while it was required
      fHistNEvents->Fill(2);
      PostData(1,fOutput);
      return;
    }
  }
    TClonesArray *arrayBranch = nullptr, *arrayD0toKpi = nullptr;
  
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD
    // event in memory rather than the input (ESD) event.
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)
    ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      if(fDecayChannel == kDplustoKpipi || fDecayChannel == kDstoKKpi)
        arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
      else if(fDecayChannel == kD0toKpi)
        arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
      else if(fDecayChannel == kDstartoKpipi) {
        arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
        arrayD0toKpi=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
      }
    }
  }
  else {
    if(fDecayChannel == kDplustoKpipi || fDecayChannel == kDstoKKpi)
      arrayBranch=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    else if(fDecayChannel == kD0toKpi)
      arrayBranch=(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
    else if(fDecayChannel == kDstartoKpipi) {
      arrayBranch=(TClonesArray*)aod->GetList()->FindObject("Dstar");
      arrayD0toKpi=(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
    }
  }
  
  if (!arrayBranch) {
    AliError("Could not find array of HF vertices");
    PostData(1,fOutput);
    return;
  }
  
  if (!fPIDresp) fPIDresp = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
  
  AliAODVertex *aodVtx = (AliAODVertex*)aod->GetPrimaryVertex();
  if (!aodVtx || TMath::Abs(aod->GetMagneticField())<0.001) {
    AliDebug(3, "The event was skipped due to missing vertex or magnetic field issue");
    PostData(1,fOutput);
    return;
  }
  fHistNEvents->Fill(3); // count event
  
  bool isEvSel  = fAnalysisCuts->IsEventSelected(aod);
  
  if(fAnalysisCuts->IsEventRejectedDueToTrigger()) fHistNEvents->Fill(5);
  if(fAnalysisCuts->IsEventRejectedDueToNotRecoVertex()) fHistNEvents->Fill(6);
  if(fAnalysisCuts->IsEventRejectedDueToVertexContributors()) fHistNEvents->Fill(7);
  if(fAnalysisCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()) fHistNEvents->Fill(8);
  if(fAnalysisCuts->IsEventRejectedDueToPileup()) fHistNEvents->Fill(9);
  if(fAnalysisCuts->IsEventRejectedDueToCentrality()) fHistNEvents->Fill(10);
  
  int runNumber = aod->GetRunNumber();
  
  TClonesArray *arrayMC    =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  AliAODMCHeader *mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  
  if(!arrayMC) {
    AliError("AliAnalysisTaskSEDmesonPIDSysProp::UserExec: MC particles branch not found!\n");
    PostData(1,fOutput);
    return;
  }
  if(!mcHeader) {
    AliError("AliAnalysisTaskSEDmesonPIDSysProp::UserExec: MC header branch not found!\n");
    PostData(1,fOutput);
    return;
  }
  
  if(aod->GetTriggerMask()==0 && (runNumber>=195344 && runNumber<=195677)) {
    // protection for events with empty trigger mask in p-Pb (Run1)
    PostData(1,fOutput);
    return;
  }
  if(fAnalysisCuts->GetUseCentrality()>0 && fAnalysisCuts->IsEventSelectedInCentrality(aod)!=0) {
    // events not passing the centrality selection can be removed immediately.
    PostData(1,fOutput);
    return;
  }
  double zMCVertex = mcHeader->GetVtxZ();
  if (TMath::Abs(zMCVertex) > fAnalysisCuts->GetMaxVtxZ()) {
    PostData(1,fOutput);
    return;
  }
  if(!isEvSel){
    PostData(1,fOutput);
    return;
  }
  fHistNEvents->Fill(4);
  
  int nCand = arrayBranch->GetEntriesFast();
  
  int nprongs = -1;
  int pdgcode = -1;
  int pdgDaughter[3];
  int pdg2Daughter[2];
  if(fDecayChannel == kDplustoKpipi){
    pdgcode = 411;
    nprongs  = 3;
    fPartName="Dplus";
    pdgDaughter[0]=321;
    pdgDaughter[1]=211;
    pdgDaughter[2]=211;
  }else if(fDecayChannel == kDstoKKpi){
    pdgcode = 431;
    nprongs  = 3;
    fPartName="Ds";
    pdgDaughter[0]=321;
    pdgDaughter[1]=321;
    pdgDaughter[2]=211;
  }else if(fDecayChannel == kD0toKpi){
    pdgcode = 421;
    nprongs  = 2;
    fPartName="D0";
    pdgDaughter[0]=321;
    pdgDaughter[1]=211;
  }else if(fDecayChannel == kDstartoKpipi){
    pdgcode = 413;
    nprongs  = 2;
    fPartName="Dstar";
    pdgDaughter[0]=421;
    pdgDaughter[1]=211;
    pdg2Daughter[0]=321;
    pdg2Daughter[1]=211;
  }else{
    AliError("Wrong decay setting");
    PostData(1,fOutput);
    return;
  }
  
  // vHF object is needed to call the method that refills the missing info of the candidates
  // if they have been deleted in dAOD reconstruction phase
  // in order to reduce the size of the file
  AliAnalysisVertexingHF *vHF=new AliAnalysisVertexingHF();
  
  for (int iCand = 0; iCand < nCand; iCand++) {
    
    AliAODRecoDecayHF* d = nullptr;
    AliAODRecoDecayHF2Prong* dD0 = nullptr;
    bool isDStarCand = false;
    int nDau = 0;

    if(fDecayChannel == kDplustoKpipi || fDecayChannel == kDstoKKpi) {
      d = dynamic_cast<AliAODRecoDecayHF3Prong*>(arrayBranch->UncheckedAt(iCand));
      nDau = 3;
    }
    else if(fDecayChannel == kD0toKpi) {
      d = dynamic_cast<AliAODRecoDecayHF2Prong*>(arrayBranch->UncheckedAt(iCand));
      nDau = 2;
    }
    else if(fDecayChannel == kDstartoKpipi) {
      d = dynamic_cast<AliAODRecoCascadeHF*>(arrayBranch->UncheckedAt(iCand));
      isDStarCand = true;
      nDau = 3;

      if(d && d->GetIsFilled()<1)
        dD0 = dynamic_cast<AliAODRecoDecayHF2Prong*>(arrayD0toKpi->At(d->GetProngID(1)));
      else
        dD0 = (dynamic_cast<AliAODRecoCascadeHF*>(d))->Get2Prong();
      if(!dD0)
        continue;
    }

    if(!d) continue;

    //Preselection to speed up task
    TObjArray arrDauTracks(nDau);
    AliAODTrack *track = nullptr;
    if(fDecayChannel != kDstartoKpipi) {
        for(int iDau=0; iDau<nDau; iDau++){
            AliAODTrack *track = vHF->GetProng(aod,d,iDau);
            arrDauTracks.AddAt(track,iDau);
        }
    }
    else {
        for(int iDau=0; iDau<nDau; iDau++){
            if(iDau == 0) 
                track=vHF->GetProng(aod,d,iDau); //soft pion
            else
                track=vHF->GetProng(aod,dD0,iDau-1); //D0 daughters
            arrDauTracks.AddAt(track,iDau);
        }
    }
    if(!fAnalysisCuts->PreSelect(arrDauTracks)){
        continue;
    }

    if(fDecayChannel == kDplustoKpipi || fDecayChannel == kDstoKKpi) {
      if(!vHF->FillRecoCand(aod,dynamic_cast<AliAODRecoDecayHF3Prong*>(d))) {
        fHistNEvents->Fill(13);
        continue;
      }
    }
    else if(fDecayChannel == kD0toKpi) {
      if(!vHF->FillRecoCand(aod,dynamic_cast<AliAODRecoDecayHF2Prong*>(d))) {
        fHistNEvents->Fill(13);
        continue;
      }
    }
    else if(fDecayChannel == kDstartoKpipi) {
      if(!vHF->FillRecoCasc(aod,dynamic_cast<AliAODRecoCascadeHF*>(d),isDStarCand)) {
        fHistNEvents->Fill(13);
        continue;
      }
      if(!d->GetSecondaryVtx()) continue;
    }
    
    fHistNEvents->Fill(11);
    bool unsetvtx=false;
    if(!d->GetOwnPrimaryVtx()){
      d->SetOwnPrimaryVtx(aodVtx);
      unsetvtx=true;
    }
    
    bool recVtx=false;
    AliAODVertex *origownvtx = nullptr;
    if(fAnalysisCuts->GetIsPrimaryWithoutDaughters()){
	    if(d->GetOwnPrimaryVtx()) 
        origownvtx = new AliAODVertex(*d->GetOwnPrimaryVtx());
	    if(fAnalysisCuts->RecalcOwnPrimaryVtx(d,aod))
        recVtx = true;
	    else fAnalysisCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
    }

    double ptD = d->Pt();
    double rapid  = d->Y(pdgcode);
    bool isFidAcc = fAnalysisCuts->IsInFiducialAcceptance(ptD,rapid);
    
    if(isFidAcc){
      int retCodeAnalysisCuts = fAnalysisCuts->IsSelectedPID(d);
      int retCodeAnalysisTrackCuts = fAnalysisCuts->IsSelected(d,AliRDHFCuts::kTracks,aod); //reject also 
      if(retCodeAnalysisCuts==0 || retCodeAnalysisTrackCuts==0) {
        if(unsetvtx) d->UnsetOwnPrimaryVtx();
        continue;
      }
      fHistNEvents->Fill(12);

      int mcLabel=-1;
      int orig = 0;
      if(!isDStarCand) mcLabel = d->MatchToMC(pdgcode,arrayMC,nprongs,pdgDaughter);
      else mcLabel = (dynamic_cast<AliAODRecoCascadeHF*>(d))->MatchToMC(pdgcode,421,pdgDaughter,pdg2Daughter,arrayMC);
      if(mcLabel<0) {
        if(unsetvtx) d->UnsetOwnPrimaryVtx();
        continue;
      }
      AliAODMCParticle* partD = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mcLabel));
      if(partD) orig = AliVertexingHFUtils::CheckOrigin(arrayMC,partD,true);
      if(orig<4) {
        if(unsetvtx) d->UnsetOwnPrimaryVtx();
        continue;
      }
      
      const int nDau = d->GetNDaughters();
      AliAODTrack* dautrack[nDau];
      if(fDecayChannel != kDstartoKpipi) {     
        for(int iDau=0; iDau<nDau; iDau++)
          dautrack[iDau] = dynamic_cast<AliAODTrack*>(d->GetDaughter(iDau));
      }
      else {
        AliAODRecoDecayHF2Prong* D0prong = dynamic_cast<AliAODRecoDecayHF2Prong*>((dynamic_cast<AliAODRecoCascadeHF*>(d))->Get2Prong());
        if(!D0prong) {
          if(unsetvtx) d->UnsetOwnPrimaryVtx();
          continue;
        }
        for(int iDau=0; iDau<nDau; iDau++) {
          dautrack[iDau] = dynamic_cast<AliAODTrack*>(D0prong->GetDaughter(iDau));
        }
      }

      double syst = GetDmesonPIDuncertainty(dautrack,nDau,arrayMC,ptD);
      if(syst==-999. || syst==0.) {
        if(unsetvtx) d->UnsetOwnPrimaryVtx();
        continue;
      }
      fHistSystPIDEffD->Fill(ptD,syst);
    }

    if(recVtx) fAnalysisCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
    if(unsetvtx) d->UnsetOwnPrimaryVtx();
  }

  delete vHF;
  vHF = nullptr;
  
  PostData(1,fOutput);
  return;
}

//________________________________________________________________________
bool AliAnalysisTaskSEDmesonPIDSysProp::LoadEffSystFile(TString systFileName)
{

  TFile* systfile = TFile::Open(systFileName.Data());
  if(!systfile || !systfile->IsOpen())
    AliFatal("Impossible to load single-track systematic file, check if it is correct! Exit.\n");

  if(fPIDstrategy==kConservativePID || fPIDstrategy==kStrongPID) {
    fHistEffPionTPC[0] = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hEffPionTPCDataV0tag_etaint_3sigma")->Clone()));
    fHistSystPionTPC[0] = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hRatioEffPionTPCDataV0tag_etaint_3sigma")->Clone()));
    if(!fHistSystPionTPC[0] || !fHistEffPionTPC[0]) return false;
    if(fKaonTPCHistoOpt==kKaonTOFtag) {
      fHistEffKaonTPC[0] = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hEffKaonTPCDataTOFtag_etaint_3sigma")->Clone()));
      fHistSystKaonTPC[0] = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hRatioEffKaonTPCDataTOFtag_etaint_3sigma")->Clone()));
    }
    else if(fKaonTPCHistoOpt==kKaonKinkstag) {
      fHistEffKaonTPC[0] = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hEffKaonTPCDataKinktag_etaint_3sigma")->Clone()));
      fHistSystKaonTPC[0] = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hRatioEffKaonTPCDataKinktag_etaint_3sigma")->Clone()));
    }
    if(!fHistSystKaonTPC[0] || !fHistEffKaonTPC[0]) return false;
    fHistEffPionTOF = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hEffPionTOFDataV0tag_etaint_3sigma")->Clone()));
    fHistSystPionTOF = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hRatioEffPionTOFDataV0tag_etaint_3sigma")->Clone()));
    if(!fHistSystPionTOF || !fHistEffPionTOF) return false;
    if(fKaonTOFHistoOpt==kKaonTPCtag) {
      fHistEffKaonTOF = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hEffKaonTOFDataTPCtag_etaint_3sigma")->Clone()));
      fHistSystKaonTOF = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hRatioEffKaonTOFDataTPCtag_etaint_3sigma")->Clone()));
    }
    else if(fKaonTOFHistoOpt==kSamePionV0tag) {
      fHistEffKaonTOF = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hEffPionTOFDataV0tag_etaint_3sigma")->Clone()));
      fHistSystKaonTOF = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hRatioEffPionTOFDataV0tag_etaint_3sigma")->Clone()));
    }
    if(!fHistSystKaonTOF || !fHistEffKaonTOF) return false;
    if(fPIDstrategy==kStrongPID) {
      fHistEffPionTPC[1] = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hEffPionTPCDataV0tag_etaint_2sigma")->Clone()));
      fHistSystPionTPC[1] = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hRatioEffPionTPCDataV0tag_etaint_2sigma")->Clone()));
      if(!fHistSystPionTPC[1] || !fHistEffPionTPC[1]) return false;
      if(fKaonTPCHistoOpt==kKaonTOFtag) {
        fHistEffKaonTPC[1] = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hEffKaonTPCDataTOFtag_etaint_2sigma")->Clone()));
        fHistSystKaonTPC[1] = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hRatioEffKaonTPCDataTOFtag_etaint_2sigma")->Clone()));
      }
      else if(fKaonTPCHistoOpt==kKaonKinkstag) {
        fHistEffKaonTPC[1] = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hEffKaonTPCDataKinktag_etaint_2sigma")->Clone()));
        fHistSystKaonTPC[1] = new TH1F(*static_cast<TH1F*>(systfile->Get("etaint/hRatioEffKaonTPCDataKinktag_etaint_2sigma")->Clone()));
      }
      if(!fHistSystKaonTPC[1] || !fHistEffKaonTPC[1]) return false;
    }
  }
  else if(fPIDstrategy==knSigmaPID) {
    fHistEffPionTPC[0] = new TH1F(*static_cast<TH1F*>(systfile->Get(Form("etaint/hEffPionTPCDataV0tag_%0.fsigma",fnSigma))->Clone()));
    fHistSystPionTPC[0] = new TH1F(*static_cast<TH1F*>(systfile->Get(Form("etaint/hRatioEffPionTPCDataV0tag_%0.fsigma",fnSigma))->Clone()));
    if(!fHistSystPionTPC[0] || !fHistEffPionTPC[0]) return false;
    if(fKaonTPCHistoOpt==kKaonTOFtag) {
      fHistEffKaonTPC[0] = new TH1F(*static_cast<TH1F*>(systfile->Get(Form("etaint/hEffKaonTPCDataTOFtag_%0.fsigma",fnSigma))->Clone()));
      fHistSystKaonTPC[0] = new TH1F(*static_cast<TH1F*>(systfile->Get(Form("etaint/hRatioEffKaonTPCDataTOFtag_%0.fsigma",fnSigma))->Clone()));
    }
    else if(fKaonTPCHistoOpt==kKaonKinkstag) {
      fHistEffKaonTPC[0] = new TH1F(*static_cast<TH1F*>(systfile->Get(Form("etaint/hEffKaonTPCDataKinktag_%0.fsigma",fnSigma))->Clone()));
      fHistSystKaonTPC[0] = new TH1F(*static_cast<TH1F*>(systfile->Get(Form("etaint/hRatioEffKaonTPCDataKinktag_%0.fsigma",fnSigma))->Clone()));
    }
    if(!fHistSystKaonTPC[0] || !fHistEffKaonTPC[0]) return false;
    fHistEffPionTOF = new TH1F(*static_cast<TH1F*>(systfile->Get(Form("etaint/hEffPionTOFDataV0tag_%0.fsigma",fnSigma))->Clone()));
    fHistSystPionTOF = new TH1F(*static_cast<TH1F*>(systfile->Get(Form("etaint/hRatioEffPionTOFDataV0tag_%0.fsigma",fnSigma))->Clone()));
    if(!fHistSystPionTOF || !fHistEffPionTOF) return false;
    if(fKaonTOFHistoOpt==kKaonTPCtag) {
      fHistEffKaonTOF = new TH1F(*static_cast<TH1F*>(systfile->Get(Form("etaint/hEffKaonTOFDataTPCtag_%0.fsigma",fnSigma))->Clone()));
      fHistSystKaonTOF = new TH1F(*static_cast<TH1F*>(systfile->Get(Form("etaint/hRatioEffKaonTOFDataTPCtag_%0.fsigma",fnSigma))->Clone()));
    }
    else if(fKaonTOFHistoOpt==kSamePionV0tag) {
      fHistEffKaonTOF = new TH1F(*static_cast<TH1F*>(systfile->Get(Form("etaint/hEffPionTOFDataV0tag_%0.fsigma",fnSigma))->Clone()));
      fHistSystKaonTOF = new TH1F(*static_cast<TH1F*>(systfile->Get(Form("etaint/hRatioEffPionTOFDataV0tag_%0.fsigma",fnSigma))->Clone()));
    }
    if(!fHistSystKaonTOF || !fHistEffKaonTOF) return false;
  }
  
    fHistSystPionTPC[0]->SetDirectory(0);
    fHistEffPionTPC[0]->SetDirectory(0);
    fHistSystKaonTPC[0]->SetDirectory(0);
    fHistEffKaonTPC[0]->SetDirectory(0);
    if(fPIDstrategy==kStrongPID) {
        fHistSystPionTPC[1]->SetDirectory(0);
        fHistEffPionTPC[1]->SetDirectory(0);
        fHistSystKaonTPC[1]->SetDirectory(0);
        fHistEffKaonTPC[1]->SetDirectory(0);    
    }
    fHistSystPionTOF->SetDirectory(0);
    fHistEffPionTOF->SetDirectory(0);
    fHistEffKaonTOF->SetDirectory(0);
    fHistSystKaonTOF->SetDirectory(0);

  systfile->Close();

  return true;
}

//________________________________________________________________________
double AliAnalysisTaskSEDmesonPIDSysProp::GetDmesonPIDuncertainty(AliAODTrack *track[], const int nDau, TClonesArray *arrayMC, double ptD)
{
  double syst = 0.;

  for(int iDau=0; iDau<nDau; iDau++){
    if(!track[iDau]){
      AliWarning("Daughter-particle track not found");
      return -999.;
    }
    int labDau = track[iDau]->GetLabel();
    if(labDau<0) {
      AliWarning("Daughter particle not found");
      return -999.;
    }
    AliAODMCParticle* p = dynamic_cast<AliAODMCParticle*>(arrayMC->UncheckedAt(TMath::Abs(labDau)));
    if(!p) {
      AliWarning("Daughter particle not found");
      return -999.;
    }
    
    int daupdgcode = TMath::Abs(p->GetPdgCode());
    double daupt = track[iDau]->Pt();
    double daupTPC = track[iDau]->GetTPCmomentum();
    double dauvar = -1.;
    if(fVarForProp==kPt)
      dauvar = daupt;
    else if(fVarForProp==kP)
      dauvar = daupTPC;

    bool isTPCok = false;
    bool isTOFok = false;
    if(fPIDresp->CheckPIDStatus(AliPIDResponse::kTPC,track[iDau]) == AliPIDResponse::kDetPidOk) isTPCok = true;
    if(fPIDresp->CheckPIDStatus(AliPIDResponse::kTOF,track[iDau]) == AliPIDResponse::kDetPidOk) isTOFok = true;

    int bin = fHistSystPionTPC[0]->GetXaxis()->FindBin(dauvar);
    double systTPC=0.;
    double systTOF=0.;
    double probTPC=0.;
    double probTOF=0.;
    double probTPCandTOF=0.;
    double probTPCorTOF=0.;

    if(fPIDstrategy==kConservativePID) {
      if(daupdgcode==211) {
        if(isTPCok) GetSingleTrackSystAndProb(fHistSystPionTPC[0],fHistEffPionTPC[0],bin,systTPC,probTPC);
        if(isTOFok) GetSingleTrackSystAndProb(fHistSystPionTOF,fHistEffPionTOF,bin,systTOF,probTOF);
      }
      else if(daupdgcode==321) {
        if(isTPCok) GetSingleTrackSystAndProb(fHistSystKaonTPC[0],fHistEffKaonTPC[0],bin,systTPC,probTPC);
        if(isTOFok) GetSingleTrackSystAndProb(fHistSystKaonTOF,fHistEffKaonTOF,bin,systTOF,probTOF);
      }
      probTPCorTOF = probTPC+probTOF-probTPC*probTOF;
      if(probTPCorTOF>1.e-20) syst += TMath::Sqrt((1-probTPC)*(1-probTPC)*systTOF*systTOF+(1-probTOF)*(1-probTOF)*systTPC*systTPC)/probTPCorTOF;
    }
    else if(fPIDstrategy==kStrongPID) {
      if(daupdgcode==211) {
        if(isTOFok) {
          if(isTPCok) GetSingleTrackSystAndProb(fHistSystPionTPC[0],fHistEffPionTPC[0],bin,systTPC,probTPC);
          GetSingleTrackSystAndProb(fHistSystPionTOF,fHistEffPionTOF,bin,systTOF,probTOF);
        }
        else if(isTPCok && !isTOFok) {//2sisgma cut in case of no TOF
          GetSingleTrackSystAndProb(fHistSystPionTPC[1],fHistEffPionTPC[1],bin,systTPC,probTPC);
        }
      }
      else if(daupdgcode==321) {
        if(isTOFok) {
          if(isTPCok) GetSingleTrackSystAndProb(fHistSystKaonTPC[0],fHistEffKaonTPC[0],bin,systTPC,probTPC);
          GetSingleTrackSystAndProb(fHistSystKaonTOF,fHistEffKaonTOF,bin,systTOF,probTOF);
        }
        else if(isTPCok && !isTOFok) {//2sisgma cut in case of no TOF
          GetSingleTrackSystAndProb(fHistSystKaonTPC[1],fHistEffKaonTPC[1],bin,systTPC,probTPC);
        }
      }
      probTPCorTOF = probTPC+probTOF-probTPC*probTOF;
      if(probTPCorTOF>1.e-20) syst += TMath::Sqrt((1-probTPC)*(1-probTPC)*systTOF*systTOF+(1-probTOF)*(1-probTOF)*systTPC*systTPC)/probTPCorTOF;
    }
    else if(fPIDstrategy==knSigmaPID) {
      if(daupdgcode==211) {
        if(isTPCok && isTOFok) {
          GetSingleTrackSystAndProb(fHistSystPionTPC[0],fHistEffPionTPC[0],bin,systTPC,probTPC);
          GetSingleTrackSystAndProb(fHistSystPionTOF,fHistEffPionTOF,bin,systTOF,probTOF);
        }
      }
      else if(daupdgcode==321) {
        if(isTPCok && isTOFok) {
          GetSingleTrackSystAndProb(fHistSystKaonTPC[0],fHistEffKaonTPC[0],bin,systTPC,probTPC);
          GetSingleTrackSystAndProb(fHistSystKaonTOF,fHistEffKaonTOF,bin,systTOF,probTOF);
        }
      }
      probTPCandTOF = probTPC*probTOF;
      if(probTPCandTOF>1.e-20) syst += TMath::Sqrt(probTPC*probTPC*systTOF*systTOF+probTOF*probTOF*systTPC*systTPC)/probTPCandTOF;
    }
    fHistPtDauVsD->Fill(ptD,dauvar);
  }
  
  return TMath::Abs(syst);
}

//________________________________________________________________________
void AliAnalysisTaskSEDmesonPIDSysProp::GetSingleTrackSystAndProb(TH1F* hSingleTrackSyst, TH1F* hSingleTrackEff, int bin, double &syst, double &prob)
{
  int binmax = hSingleTrackSyst->GetNbinsX();
  int binmin = 1;
  if(bin<binmax && bin>binmin) {
    syst = TMath::Abs(1-hSingleTrackSyst->GetBinContent(bin))*hSingleTrackEff->GetBinContent(bin); //absolute syst
    prob = hSingleTrackEff->GetBinContent(bin);
  }
  else if(bin>=binmax){
    syst = TMath::Abs(1-hSingleTrackSyst->GetBinContent(binmax))*hSingleTrackEff->GetBinContent(binmax); //absolute syst
    prob = hSingleTrackEff->GetBinContent(binmax);
  }
  else if(bin<=binmin){
    syst = TMath::Abs(1-hSingleTrackSyst->GetBinContent(binmax))*hSingleTrackEff->GetBinContent(binmin); //absolute syst
    prob = hSingleTrackEff->GetBinContent(binmin);
  }
}
