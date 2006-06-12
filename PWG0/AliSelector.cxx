/* $Id$ */

// The class definition in esdV0.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called everytime a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("AliSelector.C")
// Root > T->Process("AliSelector.C","some options")
// Root > T->Process("AliSelector.C+")
//

#include "AliSelector.h"

#include <TStyle.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TRegexp.h>
#include <TTime.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TFriendElement.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>

#include <AliLog.h>
#include <AliESD.h>

ClassImp(AliSelector)

AliSelector::AliSelector() :
  TSelector(),
  fTree(0),
  fESD(0),
  fCountFiles(0),
  fKineFile(0)
{
  //
  // Constructor. Initialization of pointers
  //

  AliLog::SetClassDebugLevel("AliSelector", AliLog::kDebug);
}

AliSelector::~AliSelector()
{
  //
  // Destructor
  //

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

void AliSelector::Begin(TTree*)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
}

void AliSelector::SlaveBegin(TTree* tree)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  AliDebug(AliLog::kDebug, "=======SLAVEBEGIN========");
  AliDebug(AliLog::kDebug, Form("Hostname: %s", gSystem->HostName()));
  AliDebug(AliLog::kDebug, Form("Time: %s", gSystem->Now().AsString()));

  if (tree != 0)
    Init(tree);
}

void AliSelector::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.

  AliDebug(AliLog::kDebug, "=========Init==========");

  fTree = tree;

  if (fTree == 0)
  {
    AliDebug(AliLog::kError, "ERROR: tree argument is 0.");
    return;
  }

  // Set branch address
  fTree->SetBranchAddress("ESD", &fESD);
  if (fESD != 0)
    AliDebug(AliLog::kInfo, "INFO: Found ESD branch in chain.");
}

Bool_t AliSelector::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. Typically here the branch pointers
  // will be retrieved. It is normaly not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed.

  AliDebug(AliLog::kDebug, "=========NOTIFY==========");
  AliDebug(AliLog::kDebug, Form("Hostname: %s", gSystem->HostName()));
  AliDebug(AliLog::kDebug, Form("Time: %s", gSystem->Now().AsString()));

  ++fCountFiles;
  if (fTree)
  {
    TFile *f = fTree->GetCurrentFile();
    AliDebug(AliLog::kInfo, Form("Processing %d. file %s", fCountFiles, f->GetName()));
  }
  else
  {
    AliDebug(AliLog::kError, "fTree not available");
  }

  DeleteKinematicsFile();

  return kTRUE;
}

Bool_t AliSelector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either TTree::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.

  // WARNING when a selector is used with a TChain, you must use
  //  the pointer to the current TTree to call GetEntry(entry).
  //  The entry is always the local entry number in the current tree.
  //  Assuming that fTree is the pointer to the TChain being processed,
  //  use fTree->GetTree()->GetEntry(entry).

  AliDebug(AliLog::kDebug, Form("=========PROCESS========== Entry %lld", entry));

  if (!fTree)
  {
    AliDebug(AliLog::kError, "ERROR: fTree is 0.");
    return kFALSE;
  }

  fTree->GetTree()->GetEntry(entry);

  /*
  // debugging
  if (fESD)
    AliDebug(AliLog::kDebug, Form("ESD: We have %d tracks.", fESD->GetNumberOfTracks()));

  if (fHeader)
    AliDebug(AliLog::kDebug, Form("Header: We have %d primaries.", fHeader->GetNprimary()));

  TTree* kinematics = GetKinematics();
  if (kinematics)
    AliDebug(AliLog::kDebug, Form("Kinematics: We have %lld particles.", kinematics->GetEntries()));
  */

  return kTRUE;
}

void AliSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

  DeleteKinematicsFile();
}

void AliSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  AliDebug(AliLog::kDebug, "=========TERMINATE==========");
}

TTree* AliSelector::GetKinematics()
{
  // Returns kinematics tree corresponding to current ESD active in fTree
  // Loads the kinematics from the kinematics file, the file is identified by replacing "AliESDs" to
  // "Kinematics" in the file path of the ESD file. This is a hack, to be changed!

  if (!fKineFile)
  {
    if (!fTree->GetCurrentFile())
      return 0;

    TString fileName(fTree->GetCurrentFile()->GetName());
    fileName.ReplaceAll("AliESDs", "Kinematics");

    AliDebug(AliLog::kInfo, Form("Opening %s", fileName.Data()));

    fKineFile = TFile::Open(fileName);
    if (!fKineFile)
      return 0;
  }

  return dynamic_cast<TTree*> (fKineFile->Get(Form("Event%d/TreeK", fTree->GetTree()->GetReadEntry())));
}

void AliSelector::DeleteKinematicsFile()
{
  //
  // Closes the kinematics file and deletes the pointer.
  //

  if (fKineFile)
  {
    fKineFile->Close();
    delete fKineFile;
    fKineFile = 0;
  }
}

Bool_t AliSelector::IsPrimaryCharged(TParticle* aParticle, Int_t aTotalPrimaries) const
{
  //
  // Returns if the given particle is a primary particle
  // This function or a equivalent should be available in some common place of AliRoot
  //

  // if the particle has a daughter primary, we do not want to count it
  if (aParticle->GetFirstDaughter() != -1 && aParticle->GetFirstDaughter() < aTotalPrimaries)
  {
    AliDebug(AliLog::kDebug+1, "Dropping particle because it has a daughter among the primaries.");
    return kFALSE;
  }

  Int_t pdgCode = TMath::Abs(aParticle->GetPdgCode());

  // skip quarks and gluon
  if (pdgCode <= 10 || pdgCode == 21)
  {
    AliDebug(AliLog::kDebug+1, "Dropping particle because it is a quark or gluon.");
    return kFALSE;
  }

  if (strcmp(aParticle->GetName(),"XXX") == 0)
  {
    AliDebug(AliLog::kDebug, Form("WARNING: There is a particle named XXX."));
    return kFALSE;
  }

  TParticlePDG* pdgPart = aParticle->GetPDG();

  if (strcmp(pdgPart->ParticleClass(),"Unknown") == 0)
  {
    AliDebug(AliLog::kDebug, Form("WARNING: There is a particle with an unknown particle class (pdg code %d).", pdgCode));
    return kFALSE;
  }

  if (pdgPart->Charge() == 0)
  {
    return kFALSE;
    AliDebug(AliLog::kDebug+1, "Dropping particle because it is not charged.");
  }

  return kTRUE;
}
