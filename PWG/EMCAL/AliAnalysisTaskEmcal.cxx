//
// Emcal base analysis task.
//
// Author: S.Aiola, M. Verweij

#include <RVersion.h>
#include "AliAnalysisTaskEmcal.h"

#include <TClonesArray.h>
#include <TList.h>
#include <TObject.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TFile.h>
#include <TChain.h>
#include <TKey.h>

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEMCALGeometry.h"
#include "AliESDEvent.h"
#include "AliEmcalParticle.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliVCaloTrigger.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliAnalysisUtils.h"
#include "AliEmcalTriggerPatchInfoAP.h"

#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

ClassImp(AliAnalysisTaskEmcal)

//________________________________________________________________________
AliAnalysisTaskEmcal::AliAnalysisTaskEmcal() : 
  AliAnalysisTaskSE("AliAnalysisTaskEmcal"),
  fForceBeamType(kNA),
  fGeneralHistograms(kFALSE),
  fInitialized(kFALSE),
  fCreateHisto(kTRUE),
  fCaloCellsName(),
  fCaloTriggersName(),
  fCaloTriggerPatchInfoName(),
  fMinCent(-999),
  fMaxCent(-999),
  fMinVz(-999),
  fMaxVz(-999),
  fTrackPtCut(0),
  fMinNTrack(0),
  fUseAliAnaUtils(kFALSE),
  fRejectPileup(kFALSE),
  fTklVsClusSPDCut(kFALSE),
  fOffTrigger(AliVEvent::kAny),
  fTrigClass(),
  fTriggerTypeSel(kND),
  fNbins(500),
  fMinBinPt(0),
  fMaxBinPt(250),
  fMinPtTrackInEmcal(0),
  fEventPlaneVsEmcal(-1),
  fMinEventPlane(-1e6),
  fMaxEventPlane(1e6),
  fCentEst("V0M"),
  fIsEmbedded(kFALSE),
  fIsPythia(kFALSE),
  fSelectPtHardBin(-999),
  fMinMCLabel(0),
  fMCLabelShift(0),
  fNcentBins(4),
  fNeedEmcalGeom(kTRUE),
  fParticleCollArray(),
  fClusterCollArray(),
  fTriggers(0),
  fAliAnalysisUtils(0x0),
  fIsEsd(kFALSE),
  fGeom(0),
  fTracks(0),
  fCaloClusters(0),
  fCaloCells(0),
  fCaloTriggers(0),
  fTriggerPatchInfo(0),
  fCent(0),
  fCentBin(-1),
  fEPV0(-1.0),
  fEPV0A(-1.0),
  fEPV0C(-1.0),
  fNVertCont(0),
  fBeamType(kNA),
  fPythiaHeader(0),
  fPtHard(0),
  fPtHardBin(0),
  fNTrials(0),
  fXsection(0),
  fOutput(0),
  fHistEventCount(0),
  fHistTrialsAfterSel(0),
  fHistEventsAfterSel(0),
  fHistXsectionAfterSel(0),
  fHistTrials(0),
  fHistEvents(0),
  fHistXsection(0),
  fHistPtHard(0),
  fHistCentrality(0),
  fHistZVertex(0),
  fHistEventPlane(0),
  fHistEventRejection(0),
  fHistTriggerClasses(0)
{
  // Default constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

  fParticleCollArray.SetOwner(kTRUE);
  fClusterCollArray.SetOwner(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcal::AliAnalysisTaskEmcal(const char *name, Bool_t histo) : 
  AliAnalysisTaskSE(name),
  fForceBeamType(kNA),
  fGeneralHistograms(kFALSE),
  fInitialized(kFALSE),
  fCreateHisto(histo),
  fCaloCellsName(),
  fCaloTriggersName(),
  fCaloTriggerPatchInfoName(),
  fMinCent(-999),
  fMaxCent(-999),
  fMinVz(-999),
  fMaxVz(-999),
  fTrackPtCut(0),
  fMinNTrack(0),
  fUseAliAnaUtils(kFALSE),
  fRejectPileup(kFALSE),
  fTklVsClusSPDCut(kFALSE),
  fOffTrigger(AliVEvent::kAny),
  fTrigClass(),
  fTriggerTypeSel(kND),
  fNbins(500),
  fMinBinPt(0),
  fMaxBinPt(250),
  fMinPtTrackInEmcal(0),
  fEventPlaneVsEmcal(-1),
  fMinEventPlane(-1e6),
  fMaxEventPlane(1e6),
  fCentEst("V0M"),
  fIsEmbedded(kFALSE),
  fIsPythia(kFALSE),
  fSelectPtHardBin(-999),
  fMinMCLabel(0),
  fMCLabelShift(0),
  fNcentBins(4),
  fNeedEmcalGeom(kTRUE),
  fParticleCollArray(),
  fClusterCollArray(),
  fTriggers(0),
  fAliAnalysisUtils(0x0),
  fIsEsd(kFALSE),
  fGeom(0),
  fTracks(0),
  fCaloClusters(0),
  fCaloCells(0),
  fCaloTriggers(0),
  fTriggerPatchInfo(0),
  fCent(0),
  fCentBin(-1),
  fEPV0(-1.0),
  fEPV0A(-1.0),
  fEPV0C(-1.0),
  fNVertCont(0),
  fBeamType(kNA),
  fPythiaHeader(0),
  fPtHard(0),
  fPtHardBin(0),
  fNTrials(0),
  fXsection(0),
  fOutput(0),
  fHistEventCount(0),
  fHistTrialsAfterSel(0),
  fHistEventsAfterSel(0),
  fHistXsectionAfterSel(0),
  fHistTrials(0),
  fHistEvents(0),
  fHistXsection(0),
  fHistPtHard(0),
  fHistCentrality(0),
  fHistZVertex(0),
  fHistEventPlane(0),
  fHistEventRejection(0),
  fHistTriggerClasses(0)
{
  // Standard constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

  fParticleCollArray.SetOwner(kTRUE);
  fClusterCollArray.SetOwner(kTRUE);

  if (fCreateHisto) {
    DefineOutput(1, TList::Class()); 
  }
}

//________________________________________________________________________
AliAnalysisTaskEmcal::~AliAnalysisTaskEmcal()
{
  // Destructor
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::SetClusPtCut(Double_t cut, Int_t c)
{
  AliClusterContainer *cont = GetClusterContainer(c);
  if (cont) cont->SetClusPtCut(cut);
  else AliError(Form("%s in SetClusPtCut(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::SetClusTimeCut(Double_t min, Double_t max, Int_t c)
{
  AliClusterContainer *cont = GetClusterContainer(c);
  if (cont) cont->SetClusTimeCut(min,max);
  else AliError(Form("%s in SetClusTimeCut(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::SetTrackPtCut(Double_t cut, Int_t c)
{
  AliParticleContainer *cont = GetParticleContainer(c);
  if (cont) cont->SetParticlePtCut(cut);
  else AliError(Form("%s in SetTrackPtCut(...): container %d not found",GetName(),c));

  fTrackPtCut = cut;
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::SetTrackEtaLimits(Double_t min, Double_t max, Int_t c)
{
  AliParticleContainer *cont = GetParticleContainer(c);
  if (cont) cont->SetParticleEtaLimits(min,max);
  else AliError(Form("%s in SetTrackPtCut(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::SetTrackPhiLimits(Double_t min, Double_t max, Int_t c)
{
  AliParticleContainer *cont = GetParticleContainer(c);
  if (cont) cont->SetParticlePhiLimits(min,max);
  else AliError(Form("%s in SetTrackPhiLimits(...): container %d not found",GetName(),c));
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (mgr) {
    AliVEventHandler *evhand = mgr->GetInputEventHandler();
    if (evhand) {
      if (evhand->InheritsFrom("AliESDInputHandler")) {
        fIsEsd = kTRUE;
      }
      else {
        fIsEsd = kFALSE;        
      }
    }
    else {
      AliError("Event handler not found!");
    }
  }
  else {
    AliError("Analysis manager not found!");
  }  

  if (!fCreateHisto)
    return;

  OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();

  if (fForceBeamType == kpp)
    fNcentBins = 1;

  if (!fGeneralHistograms)
    return;

  if (fIsPythia) {
    fHistTrialsAfterSel = new TH1F("fHistTrialsAfterSel", "fHistTrialsAfterSel", 11, 0, 11);
    fHistTrialsAfterSel->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistTrialsAfterSel->GetYaxis()->SetTitle("trials");
    fOutput->Add(fHistTrialsAfterSel);
    
    fHistEventsAfterSel = new TH1F("fHistEventsAfterSel", "fHistEventsAfterSel", 11, 0, 11);
    fHistEventsAfterSel->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistEventsAfterSel->GetYaxis()->SetTitle("total events");
    fOutput->Add(fHistEventsAfterSel);

    fHistXsectionAfterSel = new TProfile("fHistXsectionAfterSel", "fHistXsectionAfterSel", 11, 0, 11);
    fHistXsectionAfterSel->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistXsectionAfterSel->GetYaxis()->SetTitle("xsection");
    fOutput->Add(fHistXsectionAfterSel);
    
    fHistTrials = new TH1F("fHistTrials", "fHistTrials", 11, 0, 11);
    fHistTrials->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistTrials->GetYaxis()->SetTitle("trials");
    fOutput->Add(fHistTrials);

    fHistEvents = new TH1F("fHistEvents", "fHistEvents", 11, 0, 11);
    fHistEvents->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistEvents->GetYaxis()->SetTitle("total events");
    fOutput->Add(fHistEvents);

    fHistXsection = new TProfile("fHistXsection", "fHistXsection", 11, 0, 11);
    fHistXsection->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistXsection->GetYaxis()->SetTitle("xsection");
    fOutput->Add(fHistXsection);

    const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
    const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};
    
    for (Int_t i = 1; i < 12; i++) {
      fHistTrialsAfterSel->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
      fHistEventsAfterSel->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
      
      fHistTrials->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
      fHistXsection->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
      fHistEvents->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
    }

    fHistPtHard = new TH1F("fHistPtHard", "fHistPtHard", fNbins*2, fMinBinPt, fMaxBinPt*4);
    fHistPtHard->GetXaxis()->SetTitle("p_{T,hard} (GeV/c)");
    fHistPtHard->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistPtHard);
  }

  fHistZVertex = new TH1F("fHistZVertex","Z vertex position", 60, -30, 30);
  fHistZVertex->GetXaxis()->SetTitle("z");
  fHistZVertex->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistZVertex);

  if (fForceBeamType != kpp) {
    fHistCentrality = new TH1F("fHistCentrality","Event centrality distribution", 200, 0, 100);
    fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
    fHistCentrality->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistCentrality);
    
    fHistEventPlane = new TH1F("fHistEventPlane","Event plane", 120, -TMath::Pi(), TMath::Pi());
    fHistEventPlane->GetXaxis()->SetTitle("event plane");
    fHistEventPlane->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistEventPlane);
  }
  
  fHistEventRejection = new TH1F("fHistEventRejection","Reasons to reject event",20,0,20);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
  fHistEventRejection->SetBit(TH1::kCanRebin);
#else
  fHistEventRejection->SetCanExtend(TH1::kAllAxes);
#endif
  fHistEventRejection->GetXaxis()->SetBinLabel(1,"PhysSel");
  fHistEventRejection->GetXaxis()->SetBinLabel(2,"trigger");
  fHistEventRejection->GetXaxis()->SetBinLabel(3,"trigTypeSel");
  fHistEventRejection->GetXaxis()->SetBinLabel(4,"Cent");
  fHistEventRejection->GetXaxis()->SetBinLabel(5,"vertex contr.");
  fHistEventRejection->GetXaxis()->SetBinLabel(6,"Vz");
  fHistEventRejection->GetXaxis()->SetBinLabel(7,"trackInEmcal");
  fHistEventRejection->GetXaxis()->SetBinLabel(8,"minNTrack");
  fHistEventRejection->GetXaxis()->SetBinLabel(9,"VtxSel2013pA");
  fHistEventRejection->GetXaxis()->SetBinLabel(10,"PileUp");
  fHistEventRejection->GetXaxis()->SetBinLabel(11,"EvtPlane");
  fHistEventRejection->GetXaxis()->SetBinLabel(12,"SelPtHardBin");
  fHistEventRejection->GetXaxis()->SetBinLabel(13,"Bkg evt");
  fHistEventRejection->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistEventRejection);

  fHistTriggerClasses = new TH1F("fHistTriggerClasses","fHistTriggerClasses",3,0,3);
#if ROOT_VERSION_CODE < ROOT_VERSION(6,4,2)
  fHistTriggerClasses->SetBit(TH1::kCanRebin);
#else
  fHistTriggerClasses->SetCanExtend(TH1::kAllAxes);
#endif
  fOutput->Add(fHistTriggerClasses);

  fHistEventCount = new TH1F("fHistEventCount","fHistEventCount",2,0,2);
  fHistEventCount->GetXaxis()->SetBinLabel(1,"Accepted");
  fHistEventCount->GetXaxis()->SetBinLabel(2,"Rejected");
  fHistEventCount->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistEventCount);

  PostData(1, fOutput);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::FillGeneralHistograms()
{
  if (fIsPythia) {
    fHistEventsAfterSel->Fill(fPtHardBin, 1);
    fHistTrialsAfterSel->Fill(fPtHardBin, fNTrials);
    fHistXsectionAfterSel->Fill(fPtHardBin, fXsection);
    fHistPtHard->Fill(fPtHard);
  }

  fHistZVertex->Fill(fVertex[2]);

  if (fForceBeamType != kpp) {
    fHistCentrality->Fill(fCent);
    fHistEventPlane->Fill(fEPV0);
  }
  
  TObjArray* triggerClasses = InputEvent()->GetFiredTriggerClasses().Tokenize(" ");
  TIter next(triggerClasses);
  TObjString* triggerClass = 0;
  while ((triggerClass = static_cast<TObjString*>(next()))) {
    fHistTriggerClasses->Fill(triggerClass->GetString(), 1);
  }
  delete triggerClasses;
  triggerClasses = 0;

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  if (!fInitialized)
    ExecOnce();

  if (!fInitialized)
    return;

  if (!RetrieveEventObjects())
    return;

  if (IsEventSelected()) {
    if (fGeneralHistograms) fHistEventCount->Fill("Accepted",1);
  }
  else {
    if (fGeneralHistograms) fHistEventCount->Fill("Rejected",1);
    return;
  }

  if (fGeneralHistograms && fCreateHisto) {
    if (!FillGeneralHistograms())
      return;
  }

  if (!Run())
    return;

  if (fCreateHisto) {
    if (!FillHistograms())
      return;
  }
    
  if (fCreateHisto && fOutput) {
    // information for this iteration of the UserExec in the container
    PostData(1, fOutput);
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::AcceptCluster(AliVCluster *clus, Int_t c) const
{
  // Return true if cluster is accepted.

  if (!clus)
    return kFALSE;

  AliClusterContainer *cont = GetClusterContainer(c);
  if (!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }

  return cont->AcceptCluster(clus);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::AcceptTrack(AliVParticle *track, Int_t c) const
{
  // Return true if track is accepted.

  if (!track)
    return kFALSE;

  AliParticleContainer *cont = GetParticleContainer(c);
  if (!cont) {
    AliError(Form("%s:Container %d not found",GetName(),c));
    return 0;
  }

  return cont->AcceptParticle(track);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard)
{
  //
  // Get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  // Get the pt hard bin from the file path
  // This is to called in Notify and should provide the path to the AOD/ESD file
  // (Partially copied from AliAnalysisHelperJetTasks)

  TString file(currFile);  
  fXsec = 0;
  fTrials = 1;

  if (file.Contains(".zip#")) {
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  AliDebug(1,Form("File name: %s",file.Data()));

  // Get the pt hard bin
  TString strPthard(file);

  strPthard.Remove(strPthard.Last('/'));
  strPthard.Remove(strPthard.Last('/'));
  if (strPthard.Contains("AOD")) strPthard.Remove(strPthard.Last('/'));    
  strPthard.Remove(0,strPthard.Last('/')+1);
  if (strPthard.IsDec()) 
    pthard = strPthard.Atoi();
  else 
    AliWarning(Form("Could not extract file number from path %s", strPthard.Data()));

  // problem that we cannot really test the existance of a file in a archive so we have to live with open error message from root
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); 
  
  if (!fxsec) {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if (!fxsec) {
	// not a severe condition but inciate that we have no information
      return kFALSE;
    } else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
      if (!key) {
	fxsec->Close();
	return kFALSE;
      }
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) {
	fxsec->Close();
	return kFALSE;
      }
      fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } else { // no tree pyxsec.root
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if (!xtree) {
      fxsec->Close();
      return kFALSE;
    }
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    fTrials = ntrials;
    fXsec = xsection;
    fxsec->Close();
  }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::UserNotify()
{
  // Called when file changes.

  if (!fIsPythia || !fGeneralHistograms || !fCreateHisto)
    return kTRUE;

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliError(Form("%s - UserNotify: No current tree!",GetName()));
    return kFALSE;
  }

  Float_t xsection    = 0;
  Float_t trials      = 0;
  Int_t   pthardbin   = 0;

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliError(Form("%s - UserNotify: No current file!",GetName()));
    return kFALSE;
  }

  TChain *chain = dynamic_cast<TChain*>(tree);
  if (chain) tree = chain->GetTree();

  Int_t nevents = tree->GetEntriesFast();

  PythiaInfoFromFile(curfile->GetName(), xsection, trials, pthardbin);

  // TODO: Workaround
  if ((pthardbin < 0) || (pthardbin > 10)) pthardbin = 0;

  fHistTrials->Fill(pthardbin, trials);
  fHistXsection->Fill(pthardbin, xsection);
  fHistEvents->Fill(pthardbin, nevents);

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::ExecOnce()
{
  // Init the analysis.

  if (!InputEvent()) {
    AliError(Form("%s: Could not retrieve event! Returning!", GetName()));
    return;
  }

  if (fNeedEmcalGeom) {
    fGeom = AliEMCALGeometry::GetInstance();
    if (!fGeom) {
      AliError(Form("%s: Can not create geometry", GetName()));
      return;
    }
  }

  if (fEventPlaneVsEmcal >= 0) {
    if (fGeom) {
      Double_t ep = (fGeom->GetArm1PhiMax() + fGeom->GetArm1PhiMin()) / 2 * TMath::DegToRad() + fEventPlaneVsEmcal - TMath::Pi();
      fMinEventPlane = ep - TMath::Pi() / 4;
      fMaxEventPlane = ep + TMath::Pi() / 4;
    }
    else {
      AliWarning("Could not set event plane limits because EMCal geometry was not loaded!");
    }
  }

  //Load all requested track branches - each container knows name already
  for (Int_t i =0; i<fParticleCollArray.GetEntriesFast(); i++) {
    AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
    cont->SetArray(InputEvent());
  }

  if (fParticleCollArray.GetEntriesFast()>0) {
    fTracks = GetParticleArray(0);
    if (!fTracks) {
      AliError(Form("%s: Could not retrieve first track branch!", GetName()));
      return;
    }
  }

  //Load all requested cluster branches - each container knows name already
  for (Int_t i =0; i<fClusterCollArray.GetEntriesFast(); i++) {
    AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.At(i));
    cont->SetArray(InputEvent());
  }

  if (fClusterCollArray.GetEntriesFast()>0) {
    fCaloClusters = GetClusterArray(0);
    if (!fCaloClusters) {
      AliError(Form("%s: Could not retrieve first cluster branch!", GetName()));
      return;
    }
  }

  if (!fCaloCellsName.IsNull() && !fCaloCells) {
    fCaloCells =  dynamic_cast<AliVCaloCells*>(InputEvent()->FindListObject(fCaloCellsName));
    if (!fCaloCells) {
      AliError(Form("%s: Could not retrieve cells %s!", GetName(), fCaloCellsName.Data())); 
      return;
    }
  }

  if (!fCaloTriggersName.IsNull() && !fCaloTriggers) {
    fCaloTriggers =  dynamic_cast<AliVCaloTrigger*>(InputEvent()->FindListObject(fCaloTriggersName));
    if (!fCaloTriggers) {
      AliError(Form("%s: Could not retrieve calo triggers %s!", GetName(), fCaloTriggersName.Data())); 
      return;
    }
  }

  if (!fCaloTriggerPatchInfoName.IsNull() && !fTriggerPatchInfo) {
    fTriggerPatchInfo = GetArrayFromEvent(fCaloTriggerPatchInfoName.Data(),"AliEmcalTriggerPatchInfo");
    if (!fTriggerPatchInfo) {
      AliError(Form("%s: Could not retrieve calo trigger patch info %s!", GetName(), fCaloTriggerPatchInfoName.Data())); 
      return;
    }

  }

  fInitialized = kTRUE;
}

//_____________________________________________________
AliAnalysisTaskEmcal::BeamType AliAnalysisTaskEmcal::GetBeamType()
{
  // Get beam type : pp-AA-pA
  // ESDs have it directly, AODs get it from hardcoded run number ranges

  if (fForceBeamType != kNA)
    return fForceBeamType;

  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(InputEvent());
  if (esd) {
    const AliESDRun *run = esd->GetESDRun();
    TString beamType = run->GetBeamType();
    if (beamType == "p-p")
      return kpp;
    else if (beamType == "A-A")
      return kAA;
    else if (beamType == "p-A")
      return kpA;
    else
      return kNA;
  } else {
    Int_t runNumber = InputEvent()->GetRunNumber();
    if ((runNumber >= 136851 && runNumber <= 139517) ||  // LHC10h
	(runNumber >= 166529 && runNumber <= 170593)) {  // LHC11h
      return kAA;
    } else if ((runNumber>=188365 && runNumber <= 188366) ||   // LHC12g
	       (runNumber >= 195344 && runNumber <= 196608)) { // LHC13b-f
      return kpA;
    } else {
      return kpp;
    }
  }  
}

//________________________________________________________________________
ULong_t AliAnalysisTaskEmcal::GetTriggerList()
{
  if (!fTriggerPatchInfo)
    return 0;

  //number of patches in event
  Int_t nPatch = fTriggerPatchInfo->GetEntries();

  //loop over patches to define trigger type of event
  Int_t nG1 = 0;
  Int_t nG2 = 0;
  Int_t nJ1 = 0;
  Int_t nJ2 = 0;
  Int_t nL0 = 0;
  AliEmcalTriggerPatchInfo *patch;
  for (Int_t iPatch = 0; iPatch < nPatch; iPatch++) {
    patch = (AliEmcalTriggerPatchInfo*)fTriggerPatchInfo->At( iPatch );
    if (patch->IsGammaHigh()) nG1++;
    if (patch->IsGammaLow())  nG2++;
    if (patch->IsJetHigh()) nJ1++;
    if (patch->IsJetLow())  nJ2++;
    if (patch->IsLevel0())  nL0++;
  }

  AliDebug(2, "Patch summary: ");
  AliDebug(2, Form("Number of patches: %d", nPatch));
  AliDebug(2, Form("Jet:   low[%d], high[%d]" ,nJ2, nJ1));
  AliDebug(2, Form("Gamma: low[%d], high[%d]" ,nG2, nG1));

  ULong_t triggers(0);
  if (nL0>0) SETBIT(triggers, kL0);
  if (nG1>0) SETBIT(triggers, kG1);
  if (nG2>0) SETBIT(triggers, kG2);
  if (nJ1>0) SETBIT(triggers, kJ1);
  if (nJ2>0) SETBIT(triggers, kJ2);
  return triggers;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::HasTriggerType(TriggerType trigger)
{
  // Check if event has a given trigger type
  if(trigger==kND) {
    AliWarning(Form("%s: Requesting undefined trigger type!", GetName())); 
    return kFALSE;
  }
  //MV: removing this logic which as far as I can see doesn't make any sense
  // if(trigger & kND){
  //   return fTriggers == 0;
  // }
  return TESTBIT(fTriggers, trigger);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::IsEventSelected()
{
  // Check if event is selected

  if (fOffTrigger != AliVEvent::kAny) {
    UInt_t res = 0;
    const AliESDEvent *eev = dynamic_cast<const AliESDEvent*>(InputEvent());
    if (eev) {
      res = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    } else {
      const AliAODEvent *aev = dynamic_cast<const AliAODEvent*>(InputEvent());
      if (aev) {
        res = ((AliVAODHeader*)aev->GetHeader())->GetOfflineTrigger();
      }
    }
    if ((res & fOffTrigger) == 0) {
      if (fGeneralHistograms) fHistEventRejection->Fill("PhysSel",1);
      return kFALSE;
    }
  }

  if (!fTrigClass.IsNull()) {
    TString fired;
    const AliESDEvent *eev = dynamic_cast<const AliESDEvent*>(InputEvent());
    if (eev) {
      fired = eev->GetFiredTriggerClasses();
    } else {
      const AliAODEvent *aev = dynamic_cast<const AliAODEvent*>(InputEvent());
      if (aev) {
        fired = aev->GetFiredTriggerClasses();
      }
    }
    if (!fired.Contains("-B-")) {
      if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
      return kFALSE;
    }
    
    TObjArray *arr = fTrigClass.Tokenize("|");
    if (!arr) {
      if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
      return kFALSE;
    }
    Bool_t match = 0;
    for (Int_t i=0;i<arr->GetEntriesFast();++i) {
      TObject *obj = arr->At(i);
      if (!obj)
        continue;

      //Check if requested trigger was fired, relevant for data sets with 2 emcal trigger thresholds
      TString objStr = obj->GetName();
      if(objStr.Contains("J1") || objStr.Contains("J2") || objStr.Contains("G1") || objStr.Contains("G2")) {
        TString trigType1 = "J1";
        TString trigType2 = "J2";
        if(objStr.Contains("G")) {
          trigType1 = "G1";
          trigType2 = "G2";
        }
        if(objStr.Contains(trigType2) && fired.Contains(trigType2.Data())) { //requesting low threshold + overlap
          match = 1;
          break;
        } 
        else if(objStr.Contains(trigType1) && fired.Contains(trigType1.Data()) && !fired.Contains(trigType2.Data())) { //high threshold only
          match = 1;
          break;
        }
      } else {
        if (fired.Contains(obj->GetName())) {
          match = 1;
          break;
        }
      }
    }
    delete arr;
    if (!match) {
      if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
      return kFALSE;
    }
  }

  if (fTriggerTypeSel != kND) {
    if (!HasTriggerType(fTriggerTypeSel)) {
      if (fGeneralHistograms) fHistEventRejection->Fill("trigTypeSel",1);
      return kFALSE;
    }
  }

  if ((fMinCent != -999) && (fMaxCent != -999)) {
    if (fCent<fMinCent || fCent>fMaxCent) {
      if (fGeneralHistograms) fHistEventRejection->Fill("Cent",1);
      return kFALSE;
    }
  }

  if (fUseAliAnaUtils) {
    if (!fAliAnalysisUtils)
      fAliAnalysisUtils = new AliAnalysisUtils();
    fAliAnalysisUtils->SetMinVtxContr(2);
    fAliAnalysisUtils->SetMaxVtxZ(999);
    if(fMinVz<-10.) fMinVz = -10.; 
    if(fMinVz>10.)  fMaxVz = 10.;

    if (!fAliAnalysisUtils->IsVertexSelected2013pA(InputEvent())) {
      if (fGeneralHistograms) fHistEventRejection->Fill("VtxSel2013pA",1);
      return kFALSE;
    }

    if (fRejectPileup && fAliAnalysisUtils->IsPileUpEvent(InputEvent())) {
      if (fGeneralHistograms) fHistEventRejection->Fill("PileUp",1);
      return kFALSE;
    }

    if(fTklVsClusSPDCut && fAliAnalysisUtils->IsSPDClusterVsTrackletBG(InputEvent())) {
      if (fGeneralHistograms) fHistEventRejection->Fill("Bkg evt",1);
      return kFALSE;
    }
  }

  if ((fMinVz != -999) && (fMaxVz != -999)) {
    if (fNVertCont == 0 ) {
      if (fGeneralHistograms) fHistEventRejection->Fill("vertex contr.",1);
      return kFALSE;
    }
    Double_t vz = fVertex[2];
    if (vz<fMinVz || vz>fMaxVz) {
      if (fGeneralHistograms) fHistEventRejection->Fill("Vz",1);
      return kFALSE;
    }
  }

  if (fMinPtTrackInEmcal > 0 && fGeom) {
    Bool_t trackInEmcalOk = kFALSE;
    Int_t ntracks = GetNParticles(0);
    for (Int_t i = 0; i < ntracks; i++) {
      AliVParticle *track = GetAcceptParticleFromArray(i,0);
      if (!track)
	continue;

      Double_t phiMin = fGeom->GetArm1PhiMin() * TMath::DegToRad();
      Double_t phiMax = fGeom->GetArm1PhiMax() * TMath::DegToRad();
      Int_t runNumber = InputEvent()->GetRunNumber();
      if (runNumber>=177295 && runNumber<=197470) { //small SM masked in 2012 and 2013
	phiMin = 1.4;   
	phiMax = TMath::Pi();
      }

      if (track->Eta() < fGeom->GetArm1EtaMin() || track->Eta() > fGeom->GetArm1EtaMax() || track->Phi() < phiMin || track->Phi() > phiMax)
	continue;
      if (track->Pt() > fMinPtTrackInEmcal) {
	trackInEmcalOk = kTRUE;
	break;
      }
    }
    if (!trackInEmcalOk) {
      if (fGeneralHistograms) fHistEventRejection->Fill("trackInEmcal",1);
      return kFALSE;
    }
  }

  if (fMinNTrack > 0) {
    Int_t nTracksAcc = 0;
    Int_t ntracks = GetNParticles(0);
    for (Int_t i = 0; i < ntracks; i++) {
      AliVParticle *track = GetAcceptParticleFromArray(i,0);
      if (!track)
	continue;
      if (track->Pt() > fTrackPtCut) {
	nTracksAcc++;
	if (nTracksAcc>=fMinNTrack)
	  break;
      }
    }
    if (nTracksAcc<fMinNTrack) {
      if (fGeneralHistograms) fHistEventRejection->Fill("minNTrack",1);
      return kFALSE;
    }
  }

  if (!(fEPV0 > fMinEventPlane && fEPV0 <= fMaxEventPlane) &&
      !(fEPV0 + TMath::Pi() > fMinEventPlane && fEPV0 + TMath::Pi() <= fMaxEventPlane) &&
      !(fEPV0 - TMath::Pi() > fMinEventPlane && fEPV0 - TMath::Pi() <= fMaxEventPlane)) 
    {
      if (fGeneralHistograms) fHistEventRejection->Fill("EvtPlane",1);
      return kFALSE;
    }

  if (fSelectPtHardBin != -999 && fSelectPtHardBin != fPtHardBin)  {
      if (fGeneralHistograms) fHistEventRejection->Fill("SelPtHardBin",1);
      return kFALSE;
    }
  
  return kTRUE;
}

//________________________________________________________________________
TClonesArray *AliAnalysisTaskEmcal::GetArrayFromEvent(const char *name, const char *clname)
{
  // Get array from event.

  TClonesArray *arr = 0;
  TString sname(name);
  if (!sname.IsNull()) {
    arr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(sname));
    if (!arr) {
      AliWarning(Form("%s: Could not retrieve array with name %s!", GetName(), name)); 
      return 0;
    }
  } else {
    return 0;
  }

  if (!clname)
    return arr;

  TString objname(arr->GetClass()->GetName());
  TClass cls(objname);
  if (!cls.InheritsFrom(clname)) {
    AliWarning(Form("%s: Objects of type %s in %s are not inherited from %s!", 
                    GetName(), cls.GetName(), name, clname)); 
    return 0;
  }
  return arr;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::RetrieveEventObjects()
{
  // Retrieve objects from event.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
  fNVertCont = 0;

  const AliVVertex *vert = InputEvent()->GetPrimaryVertex();
  if (vert) {
    vert->GetXYZ(fVertex);
    fNVertCont = vert->GetNContributors();
  }

  fBeamType = GetBeamType();

  if (fBeamType == kAA || fBeamType == kpA ) {
    AliCentrality *aliCent = InputEvent()->GetCentrality();
    if (aliCent) {
      fCent = aliCent->GetCentralityPercentile(fCentEst.Data());
      if (fNcentBins==4) {
	if      (fCent >=  0 && fCent <   10) fCentBin = 0;
	else if (fCent >= 10 && fCent <   30) fCentBin = 1;
	else if (fCent >= 30 && fCent <   50) fCentBin = 2;
	else if (fCent >= 50 && fCent <= 100) fCentBin = 3; 
	else {
	  AliWarning(Form("%s: Negative centrality: %f. Assuming 99", GetName(), fCent));
	  fCentBin = fNcentBins-1;
	}
      } else {
	Double_t centWidth = (fMaxCent-fMinCent)/(Double_t)fNcentBins;
	if(centWidth>0.)
	  fCentBin = TMath::FloorNint(fCent/centWidth);
	else 
	  fCentBin = 0;
	if (fCentBin>=fNcentBins) {
	  AliWarning(Form("%s: fCentBin too large: cent = %f fCentBin = %d. Assuming 99", GetName(),fCent,fCentBin));
	  fCentBin = fNcentBins-1;
	}
      }
    } else {
      AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
      fCentBin = fNcentBins-1;
    }
    AliEventplane *aliEP = InputEvent()->GetEventplane();
    if (aliEP) {
      fEPV0  = aliEP->GetEventplane("V0" ,InputEvent());
      fEPV0A = aliEP->GetEventplane("V0A",InputEvent());
      fEPV0C = aliEP->GetEventplane("V0C",InputEvent());
    } else {
      AliWarning(Form("%s: Could not retrieve event plane information!", GetName()));
    }
  } else {
    fCent = 99;
    fCentBin = 0;
  }

  if (fIsPythia) {

    if (MCEvent()) {
      fPythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
      if (!fPythiaHeader) {
	// Check if AOD
	AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

	if (aodMCH) {
	  for (UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
	    fPythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
	    if (fPythiaHeader) break;
	  }
	}
      }
    }

    if (fPythiaHeader) {
      fPtHard = fPythiaHeader->GetPtHard();

      const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
      const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};
      for (fPtHardBin = 0; fPtHardBin < 11; fPtHardBin++) {
	if (fPtHard >= ptHardLo[fPtHardBin] && fPtHard < ptHardHi[fPtHardBin])
	  break;
      }

      fXsection = fPythiaHeader->GetXsection();
      fNTrials = fPythiaHeader->Trials();
    }
  }

  fTriggers = GetTriggerList();

  return kTRUE;
}

//________________________________________________________________________
AliParticleContainer* AliAnalysisTaskEmcal::AddParticleContainer(const char *n) 
{
  // Add particle container
  // will be called in AddTask macro

  TString tmp = TString(n);
  if (tmp.IsNull()) return 0;

  AliParticleContainer *cont = 0x0;
  cont = new AliParticleContainer();
  cont->SetArrayName(n);
  TString contName = cont->GetArrayName();
 
  fParticleCollArray.Add(cont);

  return cont;
}

//________________________________________________________________________
AliClusterContainer* AliAnalysisTaskEmcal::AddClusterContainer(const char *n) 
{
  // Add cluster container
  // will be called in AddTask macro

  TString tmp = TString(n);
  if (tmp.IsNull()) return 0;

  AliClusterContainer *cont = 0x0;
  cont = new AliClusterContainer();
  cont->SetArrayName(n);

  fClusterCollArray.Add(cont);

  return cont;
}

//________________________________________________________________________
AliParticleContainer* AliAnalysisTaskEmcal::GetParticleContainer(Int_t i) const 
{
  // Get i^th particle container

  if (i<0 || i>fParticleCollArray.GetEntriesFast()) return 0;
  AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.At(i));
  return cont;
}

//________________________________________________________________________
AliClusterContainer* AliAnalysisTaskEmcal::GetClusterContainer(Int_t i) const 
{
  // Get i^th cluster container

  if (i<0 || i>fClusterCollArray.GetEntriesFast()) return 0;
  AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.At(i));
  return cont;
}

//________________________________________________________________________
AliParticleContainer* AliAnalysisTaskEmcal::GetParticleContainer(const char *name) const 
{
  // Get particle container with name

  AliParticleContainer *cont = static_cast<AliParticleContainer*>(fParticleCollArray.FindObject(name));
  return cont;
}

//________________________________________________________________________
AliClusterContainer* AliAnalysisTaskEmcal::GetClusterContainer(const char *name) const 
{
  // Get cluster container with name

  AliClusterContainer *cont = static_cast<AliClusterContainer*>(fClusterCollArray.FindObject(name));
  return cont;
}

//________________________________________________________________________
TClonesArray* AliAnalysisTaskEmcal::GetParticleArray(Int_t i) const 
{
  // Get i^th TClonesArray with AliVParticle

  AliParticleContainer *cont = GetParticleContainer(i);
  if (!cont) {
    AliError(Form("%s: Particle container %d not found",GetName(),i));
    return 0;
  }
  TString contName = cont->GetArrayName();
  return cont->GetArray();
}

//________________________________________________________________________
TClonesArray* AliAnalysisTaskEmcal::GetClusterArray(Int_t i) const 
{
  // Get i^th TClonesArray with AliVCluster

  AliClusterContainer *cont = GetClusterContainer(i);
  if (!cont) {
    AliError(Form("%s:Cluster container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetArray();
}

//________________________________________________________________________
AliVParticle* AliAnalysisTaskEmcal::GetAcceptParticleFromArray(Int_t p, Int_t c) const 
{
  // Get particle p if accepted from  container c
  // If particle not accepted return 0

  AliParticleContainer *cont = GetParticleContainer(c);
  if (!cont) {
    AliError(Form("%s: Particle container %d not found",GetName(),c));
    return 0;
  }
  AliVParticle *vp = cont->GetAcceptParticle(p);

  return vp;
}

//________________________________________________________________________
AliVCluster* AliAnalysisTaskEmcal::GetAcceptClusterFromArray(Int_t cl, Int_t c) const 
{
  // Get particle p if accepted from  container c
  // If particle not accepted return 0

  AliClusterContainer *cont = GetClusterContainer(c);
  if (!cont) {
    AliError(Form("%s: Cluster container %d not found",GetName(),c));
    return 0;
  }
  AliVCluster *vc = cont->GetAcceptCluster(cl);

  return vc;
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcal::GetNParticles(Int_t i) const 
{
  // Get number of entries in particle array i

  AliParticleContainer *cont = GetParticleContainer(i);
  if (!cont) {
    AliError(Form("%s: Particle container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetNEntries();
}

//________________________________________________________________________
Int_t AliAnalysisTaskEmcal::GetNClusters(Int_t i) const 
{
  // Get number of entries in cluster array i

  AliClusterContainer *cont = GetClusterContainer(i);
  if (!cont) {
    AliError(Form("%s: Cluster container %d not found",GetName(),i));
    return 0;
  }
  return cont->GetNEntries();
}

//________________________________________________________________________
AliEmcalTriggerPatchInfo* AliAnalysisTaskEmcal::GetMainTriggerPatch(TriggerCategory trigger, Bool_t doSimpleOffline)
{
  // Get main trigger match
  //
  // For the selection of the main trigger patch, high and low threshold triggers of a given category are grouped
  // If there are more than 1 main patch of a given trigger category (i.e. different high and low threshold patches),
  // the highest one according to the ADC value is taken. In case doSimpleOffline is true, then only the patches from
  // the simple offline trigger are used.

  if (!fTriggerPatchInfo) {
    AliError(Form("%s: fTriggerPatchInfo not available",GetName()));
    return 0;
  }

  //number of patches in event
  Int_t nPatch = fTriggerPatchInfo->GetEntries();

  //extract main trigger patch(es)
  AliEmcalTriggerPatchInfo *patch(NULL), *selected(NULL);
  for (Int_t iPatch = 0; iPatch < nPatch; iPatch++) {

    patch = (AliEmcalTriggerPatchInfo*)fTriggerPatchInfo->At( iPatch );
    if (patch->IsMainTrigger()) {
      if(doSimpleOffline){
        if(patch->IsOfflineSimple()){
          switch(trigger){
          case kTriggerLevel0:
            // option not yet implemented in the trigger maker
            if(patch->IsLevel0()) selected = patch;
            break;
          case kTriggerLevel1Jet: 
            if(patch->IsJetHighSimple() || patch->IsJetLowSimple()){
              if(!selected) selected = patch;
              else if(patch->GetADCOfflineAmp() > selected->GetADCOfflineAmp()) selected = patch;
            }
            break;
          case kTriggerLevel1Gamma:
	    if(patch->IsGammaHighSimple() || patch->IsGammaLowSimple()){
	      if(!selected) selected = patch;
	      else if(patch->GetADCOfflineAmp() > selected->GetADCOfflineAmp()) selected = patch;
	    }
	    break;
          default:   // Silence compiler warnings
            AliError("Untreated case: Main Patch is recalculated; should be in 'else' branch");
          };
        }
      } else {  // Not OfflineSimple
        switch(trigger){
        case kTriggerLevel0:
            if(patch->IsLevel0()) selected = patch;
            break;
        case kTriggerLevel1Jet:
            if(patch->IsJetHigh() || patch->IsJetLow()){
              if(!selected) selected = patch;
              else if (patch->GetADCAmp() > selected->GetADCAmp()) 
                selected = patch;
            }
            break;
        case kTriggerLevel1Gamma:
            if(patch->IsGammaHigh() || patch->IsGammaLow()){
              if(!selected) selected = patch;
              else if (patch->GetADCAmp() > selected->GetADCAmp()) 
                selected = patch;
            }
            break;
         default:
            AliError("Untreated case: Main Patch is recalculated; should be in 'else' branch");
        };
      }
    }
    else if ((trigger == kTriggerRecalcJet &&  patch->IsRecalcJet()) || 
             (trigger == kTriggerRecalcGamma && patch->IsRecalcGamma())) {  // recalculated patches
      if (doSimpleOffline && patch->IsOfflineSimple()) {
        if(!selected) selected = patch;
        else if (patch->GetADCOfflineAmp() > selected->GetADCOfflineAmp())  // this in fact should not be needed, but we have it in teh other branches as well, so keeping it for compleness
          selected = patch;
      }
      else if (!doSimpleOffline && !patch->IsOfflineSimple()) {
        if(!selected) selected = patch;
        else if (patch->GetADCAmp() > selected->GetADCAmp()) 
          selected = patch;
      }
    }
  }
  return selected;
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::AddObjectToEvent(TObject *obj)
{
  // Add object to event

  if (!(InputEvent()->FindListObject(obj->GetName()))) {
    InputEvent()->AddObject(obj);
  } else {
    AliFatal(Form("%s: Container with name %s already present. Aborting", GetName(), obj->GetName()));
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcal::SetRejectionReasonLabels(TAxis* axis)
{
  axis->SetBinLabel(1,  "NullObject");
  axis->SetBinLabel(2,  "Pt");
  axis->SetBinLabel(3,  "Acceptance");
  axis->SetBinLabel(4,  "BitMap");
  axis->SetBinLabel(5,  "Bit4");
  axis->SetBinLabel(6,  "Bit5");
  axis->SetBinLabel(7,  "Bit6");
  axis->SetBinLabel(8,  "NotHybridTrack");
  axis->SetBinLabel(9,  "MCFlag");
  axis->SetBinLabel(10, "MCGenerator");
  axis->SetBinLabel(11, "ChargeCut");
  axis->SetBinLabel(12, "MinDistanceTPCSectorEdge");
  axis->SetBinLabel(13, "MinMCLabelAccept");
  axis->SetBinLabel(14, "IsEMCal");
  axis->SetBinLabel(15, "Time");
  axis->SetBinLabel(16, "Energy");
  axis->SetBinLabel(17, "Bit16");
  axis->SetBinLabel(18, "Bit17");
  axis->SetBinLabel(19, "Area");
  axis->SetBinLabel(20, "AreaEmc");
  axis->SetBinLabel(21, "ZLeadingCh");
  axis->SetBinLabel(22, "ZLeadingEmc");
  axis->SetBinLabel(23, "NEF");
  axis->SetBinLabel(24, "MinLeadPt");
  axis->SetBinLabel(25, "MaxTrackPt");
  axis->SetBinLabel(26, "MaxClusterPt");
  axis->SetBinLabel(27, "Flavour");
  axis->SetBinLabel(28, "TagStatus");
  axis->SetBinLabel(29, "MinNConstituents");
  axis->SetBinLabel(30, "Bit29");
  axis->SetBinLabel(31, "Bit30");
  axis->SetBinLabel(32, "Bit31");
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcal::GetParallelFraction(AliVParticle* part1, AliVParticle* part2)
{
  // Calculates the fraction of momentum of part 1 w.r.t. part 2 in the direction of part 2.
  
  TVector3 vect1(part1->Px(), part1->Py(), part1->Pz());
  TVector3 vect2(part2->Px(), part2->Py(), part2->Pz());
  Double_t z = (vect1 * vect2) / (vect2 * vect2);
  return z;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcal::GetParallelFraction(const TVector3& vect1, AliVParticle* part2)
{
  // Calculates the fraction of momentum of vect 1 w.r.t. part 2 in the direction of part 2.
  
  TVector3 vect2(part2->Px(), part2->Py(), part2->Pz());
  Double_t z = (vect1 * vect2) / (vect2 * vect2);
  return z;
}
