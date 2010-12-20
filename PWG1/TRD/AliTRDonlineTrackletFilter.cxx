#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"

#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODHandler.h"
#include "AliMCEventHandler.h"
#include "AliLog.h"
#include "AliESDTrdTrack.h"

#include "AliTRDtrackletMCM.h"
#include "AliTRDtrackletWord.h"
#include "AliVParticle.h"
#include "AliMCParticle.h" 

#include "AliTRDonlineTrackletFilter.h"

ClassImp(AliTRDonlineTrackletFilter)

AliTRDonlineTrackletFilter::AliTRDonlineTrackletFilter(const char *name) :
  AliAnalysisTask(name, ""),
  fESD(0x0),
  fInputHandler(0x0),
  fInputEvent(0x0),
  fOutputAOD(0x0),
  fMCEvent(0x0),
  fTrackletsRaw(new TClonesArray("AliTRDtrackletWord")),
  fTrackletsSim(new TClonesArray("AliTRDtrackletMCM")),
  fTrackletTree(0x0),
  fGeo(new AliTRDgeometry),
  fNevent(0),
  fPath(""),
  fTrackletFile(0x0),
  fNEventsPerFile(0),
  fEvent(0),
  fFileNumber(0),
  fTrackletTreeSim(0x0),
  fTrackletTreeRaw(0x0)
{
  // ctor

  DefineInput(0, TChain::Class());

  DefineOutput(0, TTree::Class()); 
  DefineOutput(1, TTree::Class());
}

AliTRDonlineTrackletFilter::~AliTRDonlineTrackletFilter()
{
  // dtor

  delete fTrackletsRaw;
  delete fTrackletsSim;
  delete fGeo;
}

void AliTRDonlineTrackletFilter::ConnectInputData(const Option_t */* option */)
{
  fInputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if (fInputHandler)
    fInputEvent = fInputHandler->GetEvent();

  AliMCEventHandler *mcH = (AliMCEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
  if (mcH)
    fMCEvent = mcH->MCEvent();
}

void AliTRDonlineTrackletFilter::CreateOutputObjects()
{
  OpenFile(1); 
  
  fTrackletTree = new TTree("tracklets", "on-line tracklets");
  fTrackletTree->Branch("tracklets_sim", fTrackletsSim);
  fTrackletTree->Branch("tracklets_raw", fTrackletsRaw);
}

Bool_t AliTRDonlineTrackletFilter::Notify()
{

  TString filename(AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile()->GetName());

  AliInfo(Form("Now reading from %s", filename.Data()));

  if (filename.Contains("AliAOD.root")) 
    filename.ReplaceAll("AliAOD.root", "");
  else if (filename.Contains("AliESDs.root"))
    filename.ReplaceAll("AliESDs.root", "");
  else if (filename.Contains("galice.root"))
    filename.ReplaceAll("galice.root", "");
  else if (filename.BeginsWith("root:"))
    filename.Append("?ZIP=");

  fPath = filename;

  fTrackletFile = TFile::Open(Form("%sTRD.Tracklets.root", fPath.Data()));

  if (!fTrackletFile) {
    AliError("No tracklet file");
    return kFALSE;
  }

  fNEventsPerFile = fTrackletFile->GetNkeys() - fTrackletFile->GetNProcessIDs(); //???

  fEvent = -1;
  fFileNumber = 0;

  return kTRUE;
}


void AliTRDonlineTrackletFilter::Exec(const Option_t * /* option */)
{
  // execute this for each event

  if (!LoadEvent())
    return;

  // ----- simulated tracklets -----
  if (fTrackletTreeSim) {
    AliTRDtrackletMCM *trkl = 0x0;

    TBranch *br = fTrackletTreeSim->GetBranch("mcmtrklbranch");
    br->SetAddress(&trkl);
    
    for (Int_t iTracklet = 0; iTracklet < br->GetEntries(); iTracklet++) {
      br->GetEntry(iTracklet);
      new ((*fTrackletsSim)[fTrackletsSim->GetEntriesFast()]) AliTRDtrackletMCM(*trkl);
    }
  }

  // ----- raw tracklets -----
  if (fTrackletTreeRaw) {
    Int_t hc;
    TClonesArray *trklArray = 0x0;
    fTrackletTreeRaw->SetBranchAddress("hc", &hc);
    fTrackletTreeRaw->SetBranchAddress("trkl", &trklArray);
    for (Int_t iHCidx = 0; iHCidx < fTrackletTreeRaw->GetEntries(); iHCidx++) {
      fTrackletTreeRaw->GetEntry(iHCidx);
      for (Int_t iTracklet = 0; iTracklet < trklArray->GetEntries(); iTracklet++) {
        AliTRDtrackletWord *trklWord = (AliTRDtrackletWord*) ((*trklArray)[iTracklet]);
	trklWord->SetDetector(hc/2);
        new ((*fTrackletsRaw)[fTrackletsRaw->GetEntriesFast()]) AliTRDtrackletWord(*trklWord);
      }
    }
  }

  AliInfo(Form("%i tracklets", fTrackletsSim->GetEntriesFast()));
  fTrackletTree->SetBranchAddress("tracklets_sim", &fTrackletsSim);
  fTrackletTree->SetBranchAddress("tracklets_raw", &fTrackletsRaw);
  fTrackletTree->Fill();
  PostData(1, fTrackletTree);  
}

void AliTRDonlineTrackletFilter::LocalInit()
{

}

void AliTRDonlineTrackletFilter::Terminate(const Option_t * /* option */)
{

}

Bool_t AliTRDonlineTrackletFilter::LoadEvent()
{
  // load tracklets for the current event

  // ----- cleaning -----
  fTrackletsSim->Delete();
  fTrackletsRaw->Delete();

  // ----- initialization -----
  if (!fInputEvent) {
    AliError("No event found!");
    return kFALSE;
  }
  fESD = dynamic_cast<AliESDEvent*> (fInputEvent);

  fEvent++;
  Int_t inew = fEvent / fNEventsPerFile;
  if ( inew != fFileNumber) {
    fFileNumber++;

    delete fTrackletFile;
    fTrackletFile = TFile::Open(Form("%sTRD.Tracklets%d.root", fPath.Data(), fFileNumber));

    if (!fTrackletFile) {
      AliError("No tracklet file");
      return kFALSE;
    }
  }

  if (!fTrackletFile) {
    AliError("no tracklet file");
    return kFALSE;
  }

  // tracklets from simulation
  char treename[30];
  snprintf(treename, 30, "Event%d/tracklets", fEvent);

  fTrackletTreeSim = (TTree*) fTrackletFile->Get(treename);

  // tracklets from raw
  char treenameRaw[30];
  snprintf(treenameRaw, 30, "Event%d/tracklets-raw", fEvent);
  fTrackletTreeRaw = (TTree*) fTrackletFile->Get(treenameRaw);

  return kTRUE;
}

