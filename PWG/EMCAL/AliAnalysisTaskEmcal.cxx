// $Id: AliAnalysisTaskEmcal.cxx 56756 2012-05-30 05:03:02Z loizides $
//
// Emcal base analysis task.
//
// Author: S.Aiola

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

ClassImp(AliAnalysisTaskEmcal)

//________________________________________________________________________
AliAnalysisTaskEmcal::AliAnalysisTaskEmcal() : 
  AliAnalysisTaskSE("AliAnalysisTaskEmcal"),
  fAnaType(kTPC),
  fForceBeamType(kNA),
  fGeneralHistograms(kFALSE),
  fInitialized(kFALSE),
  fCreateHisto(kTRUE),
  fTracksName(),
  fCaloName(),
  fCaloCellsName(),
  fCaloTriggersName(),
  fMinCent(-999),
  fMaxCent(-999),
  fMinVz(-999),
  fMaxVz(-999),
  fOffTrigger(AliVEvent::kAny),
  fTrigClass(),
  fNbins(500),
  fMinBinPt(0),
  fMaxBinPt(250),
  fClusPtCut(0.15),
  fTrackPtCut(0.15),
  fTrackMinEta(-0.9),
  fTrackMaxEta(0.9),
  fTrackMinPhi(-10),
  fTrackMaxPhi(10),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fMinPtTrackInEmcal(0),
  fEventPlaneVsEmcal(-1),
  fMinEventPlane(-10),
  fMaxEventPlane(10),
  fCentEst("V0M"),
  fTrackBitMap(0),
  fClusterBitMap(0),
  fMCTrackBitMap(0),
  fMCClusterBitMap(0),
  fIsEmbedded(kFALSE),
  fIsPythia(kFALSE),
  fSelectPtHardBin(-999),
  fMinMCLabel(0),
  fMCLabelShift(0),
  fNcentBins(4),
  fGeom(0),
  fTracks(0),
  fCaloClusters(0),
  fCaloCells(0),
  fCaloTriggers(0),
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
  fOutput(0),
  fHistTrialsAfterSel(0),
  fHistEventsAfterSel(0),
  fHistTrials(0),
  fHistXsection(0),
  fHistEvents(0),
  fHistPtHard(0),
  fHistCentrality(0),
  fHistZVertex(0),
  fHistEventPlane(0)
{
  // Default constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;
}

//________________________________________________________________________
AliAnalysisTaskEmcal::AliAnalysisTaskEmcal(const char *name, Bool_t histo) : 
  AliAnalysisTaskSE(name),
  fAnaType(kTPC),
  fForceBeamType(kNA),
  fGeneralHistograms(kFALSE),
  fInitialized(kFALSE),
  fCreateHisto(histo),
  fTracksName(),
  fCaloName(),
  fCaloCellsName(),
  fCaloTriggersName(),
  fMinCent(-999),
  fMaxCent(-999),
  fMinVz(-999),
  fMaxVz(-999),
  fOffTrigger(AliVEvent::kAny),
  fTrigClass(),
  fNbins(500),
  fMinBinPt(0),
  fMaxBinPt(250),
  fClusPtCut(0.15),
  fTrackPtCut(0.15),
  fTrackMinEta(-0.9),
  fTrackMaxEta(0.9),
  fTrackMinPhi(-10),
  fTrackMaxPhi(10),
  fClusTimeCutLow(-10),
  fClusTimeCutUp(10),
  fMinPtTrackInEmcal(0),
  fEventPlaneVsEmcal(-1),
  fMinEventPlane(-10),
  fMaxEventPlane(10),
  fCentEst("V0M"),
  fTrackBitMap(0),
  fClusterBitMap(0),
  fMCTrackBitMap(0),
  fMCClusterBitMap(0),
  fIsEmbedded(kFALSE),
  fIsPythia(kFALSE),
  fSelectPtHardBin(-999),
  fMinMCLabel(0),
  fMCLabelShift(0),
  fNcentBins(4),
  fGeom(0),
  fTracks(0),
  fCaloClusters(0),
  fCaloCells(0),
  fCaloTriggers(0),
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
  fOutput(0),
  fHistTrialsAfterSel(0),
  fHistEventsAfterSel(0),
  fHistTrials(0),
  fHistXsection(0),
  fHistEvents(0),
  fHistPtHard(0),
  fHistCentrality(0),
  fHistZVertex(0),
  fHistEventPlane(0)
{
  // Standard constructor.

  fVertex[0] = 0;
  fVertex[1] = 0;
  fVertex[2] = 0;

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
void AliAnalysisTaskEmcal::UserCreateOutputObjects()
{
  // Create user output.
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
    
    fHistTrials = new TH1F("fHistTrials", "fHistTrials", 11, 0, 11);
    fHistTrials->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistTrials->GetYaxis()->SetTitle("trials");
    fOutput->Add(fHistTrials);

    fHistXsection = new TProfile("fHistXsection", "fHistXsection", 11, 0, 11);
    fHistXsection->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistXsection->GetYaxis()->SetTitle("xsection");
    fOutput->Add(fHistXsection);

    fHistEvents = new TH1F("fHistEvents", "fHistEvents", 11, 0, 11);
    fHistEvents->GetXaxis()->SetTitle("p_{T} hard bin");
    fHistEvents->GetYaxis()->SetTitle("total events");
    fOutput->Add(fHistEvents);

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

  fHistCentrality = new TH1F("fHistCentrality","Event centrality distribution", 200, 0, 100);
  fHistCentrality->GetXaxis()->SetTitle("Centrality (%)");
  fHistCentrality->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistCentrality);

  fHistZVertex = new TH1F("fHistZVertex","Z vertex position", 60, -30, 30);
  fHistZVertex->GetXaxis()->SetTitle("z");
  fHistZVertex->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistZVertex);

  fHistEventPlane = new TH1F("fHistEventPlane","Event plane", 120, -TMath::Pi(), TMath::Pi());
  fHistEventPlane->GetXaxis()->SetTitle("event plane");
  fHistEventPlane->GetYaxis()->SetTitle("counts");
  fOutput->Add(fHistEventPlane);

  PostData(1, fOutput);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::FillGeneralHistograms()
{
  if (fIsPythia) {
    fHistEventsAfterSel->SetBinContent(fPtHardBin + 1, fHistEventsAfterSel->GetBinContent(fPtHardBin + 1) + 1);
    fHistTrialsAfterSel->SetBinContent(fPtHardBin + 1, fHistTrialsAfterSel->GetBinContent(fPtHardBin + 1) + fNTrials);
    fHistPtHard->Fill(fPtHard);
  }

  fHistCentrality->Fill(fCent);
  fHistZVertex->Fill(fVertex[2]);
  fHistEventPlane->Fill(fEPV0);

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

  if (!IsEventSelected()) 
    return;

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
Bool_t AliAnalysisTaskEmcal::AcceptCluster(AliVCluster *clus) const
{
  // Return true if cluster is accepted.

  if (!clus)
    return kFALSE;
      
  if (!clus->IsEMCAL())
    return kFALSE;

  if (clus->GetLabel() > fMinMCLabel) {
    if (clus->TestBits(fMCClusterBitMap) != (Int_t)fMCClusterBitMap) {
      AliDebug(2,"MC Cluster not accepted because of bit map.");
      return kFALSE;
    }
  }
  else {
    if (clus->TestBits(fClusterBitMap) != (Int_t)fClusterBitMap) {
      AliDebug(2,"Cluster not accepted because of bit map.");
      return kFALSE;
    }
  }

  if (clus->GetTOF() > fClusTimeCutUp || clus->GetTOF() < fClusTimeCutLow)
    return kFALSE;

  TLorentzVector nPart;
  clus->GetMomentum(nPart, const_cast<Double_t*>(fVertex));

  if (nPart.Et() < fClusPtCut)
    return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::AcceptTrack(AliVParticle *track) const
{
  // Return true if track is accepted.

  if (!track)
    return kFALSE;

  if (TMath::Abs(track->GetLabel()) > fMinMCLabel) {
    if(track->TestBits(fMCTrackBitMap) != (Int_t)fMCTrackBitMap) {
      AliDebug(2,"MC Track not accepted because of bit map.");
      return kFALSE;
    }
  }
  else {
    if(track->TestBits(fTrackBitMap) != (Int_t)fTrackBitMap) {
      AliDebug(2,"Track not accepted because of bit map.");
      return kFALSE;
    }
  }

  if (track->Pt() < fTrackPtCut) {
    AliDebug(2,Form("Track not accepted because of minimum pt cut (%f < %f).", track->Pt(), fTrackPtCut));
    return kFALSE;
  }

  if (track->Eta() < fTrackMinEta || track->Eta() > fTrackMaxEta || 
      track->Phi() < fTrackMinPhi || track->Phi() > fTrackMaxPhi) {
    AliDebug(2,"Track not accepted because of acceptance cut.");
    return kFALSE;
  }
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcal::AcceptEmcalPart(AliEmcalParticle *part) const
{
  // Return true if EMCal particle is accepted.

  if (!part)
    return kFALSE;

  if (part->IsTrack()) {
    if (part->IsMC()) { 
      if (part->TestBits(fMCTrackBitMap) != (Int_t)fMCTrackBitMap)
	return kFALSE;
    }
    else {
      if (part->TestBits(fTrackBitMap) != (Int_t)fTrackBitMap)
	return kFALSE;
    }

    if (part->Pt() < fTrackPtCut)
      return kFALSE;

    if (part->Eta() < fTrackMinEta || part->Eta() > fTrackMaxEta || 
	part->Phi() < fTrackMinPhi || part->Phi() > fTrackMaxPhi)
      return kFALSE;
  }

  if (part->IsCluster()) {
    if (part->IsMC(fMinMCLabel)) { 
      if (part->TestBits(fMCClusterBitMap) != (Int_t)fMCClusterBitMap)
	return kFALSE;
    }
    else {
      if (part->TestBits(fClusterBitMap) != (Int_t)fClusterBitMap)
	return kFALSE;
    }

    if (!part->IsEMCAL())
      return kFALSE;

    if (part->Pt() < fClusPtCut)
      return kFALSE;
  }

  return kTRUE;
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

  if(file.Contains(".zip#")){
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  }
  else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  Printf("%s",file.Data());

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

  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
  if(!fxsec){
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if(!fxsec){
	// not a severe condition but inciate that we have no information
      return kFALSE;
    }
    else{
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
      if(!key){
	fxsec->Close();
	return kFALSE;
      }
      TList *list = dynamic_cast<TList*>(key->ReadObj());
      if(!list){
	fxsec->Close();
	return kFALSE;
      }
      fXsec = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } // no tree pyxsec.root
  else {
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if(!xtree){
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
  if (!fIsPythia || !fGeneralHistograms || !fCreateHisto)
    return kTRUE;

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliError(Form("%s - UserNotify: No current tree!",GetName()));
    return kFALSE;
  }

  Float_t xsection = 0;
  Float_t trials   = 0;
  Int_t   pthard   = 0;

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliError(Form("%s - UserNotify: No current file!",GetName()));
    return kFALSE;
  }

  TChain *chain = dynamic_cast<TChain*>(tree);
  if (chain)
    tree = chain->GetTree();

  Int_t nevents = tree->GetEntriesFast();

  PythiaInfoFromFile(curfile->GetName(), xsection, trials, pthard);

  fHistTrials->Fill(pthard, trials);
  fHistXsection->Fill(pthard, xsection);
  fHistEvents->Fill(pthard, nevents);

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

  fGeom = AliEMCALGeometry::GetInstance();
  if (!fGeom) {
    AliError(Form("%s: Can not create geometry", GetName()));
    return;
  }

  if (fEventPlaneVsEmcal >= 0) {
    Double_t ep = (fGeom->GetArm1PhiMax() + fGeom->GetArm1PhiMin()) / 2 * TMath::DegToRad() + fEventPlaneVsEmcal - TMath::Pi();
    fMinEventPlane = ep - TMath::Pi() / 4;
    fMaxEventPlane = ep + TMath::Pi() / 4;
  }

  if (!fCaloName.IsNull() && !fCaloClusters) {
    fCaloClusters =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCaloName));
    if (!fCaloClusters) {
      AliError(Form("%s: Could not retrieve clusters %s!", GetName(), fCaloName.Data())); 
      return;
    } else {
      TClass *cl = fCaloClusters->GetClass();
      if (!cl->GetBaseClass("AliVCluster") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVCluster nor AliEmcalParticle objects!", GetName(), fCaloName.Data())); 
	fCaloClusters = 0;
	return;
      }
    }
  }

  if (!fTracksName.IsNull() && !fTracks) {
    fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
    if (!fTracks) {
      AliError(Form("%s: Could not retrieve tracks %s!", GetName(), fTracksName.Data())); 
      return;
    } else {
      TClass *cl = fTracks->GetClass();
      if (!cl->GetBaseClass("AliVParticle") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVParticle nor AliEmcalParticle objects!", GetName(), fTracksName.Data())); 
	fTracks = 0;
	return;
      }
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
	(runNumber >= 166529 && runNumber <= 170593))    // LHC11h
    { 
      return kAA;
    } 
    else if ((runNumber>=188365 && runNumber <= 188366) || // LHC12g
	     (runNumber >= 195344 && runNumber <= 196608))  // LHC13b-f
    {
      return kpA;
    } else {
      return kpp;
    }
  }  
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
        res = aev->GetHeader()->GetOfflineTrigger();
      }
    }
    if ((res & fOffTrigger) == 0)
      return kFALSE;
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
    if (!fired.Contains("-B-"))
      return kFALSE;
    TObjArray *arr = fTrigClass.Tokenize("|");
    if (!arr)
      return kFALSE;
    Bool_t match = 0;
    for (Int_t i=0;i<arr->GetEntriesFast();++i) {
      TObject *obj = arr->At(i);
      if (!obj)
        continue;
      if (fired.Contains(obj->GetName())) {
        match = 1;
        break;
      }
    }
    delete arr;
    if (!match)
      return kFALSE;
  }

  if ((fMinCent != -999) && (fMaxCent != -999)) {
    if (fCent<fMinCent)
      return kFALSE;
    if (fCent>fMaxCent)
      return kFALSE;
  }

  if ((fMinVz != -999) && (fMaxVz != -999)) {
    if (fNVertCont == 0 )
      return kFALSE;
    Double_t vz = fVertex[2];
    if (vz<fMinVz)
      return kFALSE;
    if (vz>fMaxVz)
      return kFALSE;
  }

  if (fMinPtTrackInEmcal > 0 && fTracks && fGeom) {
    Bool_t trackInEmcalOk = kFALSE;
    Int_t ntracks = fTracks->GetEntries();
    for (Int_t i = 0; i < ntracks; i++) {
      AliVTrack *track = static_cast<AliVTrack*>(fTracks->At(i));
      if (!AcceptTrack(track))
	continue;
      if (track->Eta() < fGeom->GetArm1EtaMin() || track->Eta() > fGeom->GetArm1EtaMax() ||
	  track->Phi() < fGeom->GetArm1PhiMin() * TMath::DegToRad() || track->Phi() > fGeom->GetArm1PhiMax() * TMath::DegToRad())
	continue;
      if (track->Pt() > fMinPtTrackInEmcal) {
	trackInEmcalOk = kTRUE;
	break;
      }
    }
    if (!trackInEmcalOk)
      return kFALSE;
  }


  if (!(fEPV0 > fMinEventPlane && fEPV0 <= fMaxEventPlane) &&
      !(fEPV0 + TMath::Pi() > fMinEventPlane && fEPV0 + TMath::Pi() <= fMaxEventPlane) &&
      !(fEPV0 - TMath::Pi() > fMinEventPlane && fEPV0 - TMath::Pi() <= fMaxEventPlane)) 
    return kFALSE;

  if (fSelectPtHardBin != -999 && fSelectPtHardBin != fPtHardBin) 
    return kFALSE;

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
      if(fNcentBins==4) {
	if      (fCent >=  0 && fCent <   10) fCentBin = 0;
	else if (fCent >= 10 && fCent <   30) fCentBin = 1;
	else if (fCent >= 30 && fCent <   50) fCentBin = 2;
	else if (fCent >= 50 && fCent <= 100) fCentBin = 3; 
	else {
	  AliWarning(Form("%s: Negative centrality: %f. Assuming 99", GetName(), fCent));
	  fCentBin = fNcentBins-1;
	}
      }
      else {
	Double_t centWidth = (fMaxCent-fMinCent)/(Double_t)fNcentBins;
	if(centWidth>0.)
          fCentBin = TMath::FloorNint(fCent/centWidth);
	else {
	  fCentBin = fNcentBins-1;
	  AliWarning(Form("%s: Centrality range not set. Setting fCentBin to %d", GetName(),fCentBin));
	}
	if(fCentBin>=fNcentBins) {
	  AliWarning(Form("%s: fCentBin too large: cent = %f fCentBin = %d. Assuming 99", GetName(),fCent,fCentBin));
	  fCentBin = fNcentBins-1;
	}
      }
    } else {
      AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
      fCent = 99;
      fCentBin = 3;
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
	  for(UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++) {
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
    
      fNTrials = fPythiaHeader->Trials();
    }
  }

  return kTRUE;
}
