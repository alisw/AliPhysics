#include <TSelector.h>
#include "ToyModel.C"
#ifndef __CINT__
#include <TProof.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TProofOutputFile.h>
#include <TFile.h>
#else
class TTree;
class TProof;
class TProofOutputFile;
class TFile;
#endif


struct ToyModelSelector : public TSelector
{

  ToyModelSelector()
    : fModel(0),
      fProofFile(0),
      fFile(0) 
  {}
  Int_t Version() const { return 2; }
  void Init(TTree*) { Printf("Initialize"); }
  void Begin(TTree*) { Printf("Begin it"); }
  void SlaveBegin(TTree*)
  {
    Printf("Slave begining");
    fModel = new ToyModel(fParams);

    TString loc(GetParam("PROOF_OUTPUTFILE_LOCATION"));
    fProofFile = new TProofOutputFile("dist.root",
				      (loc.IsNull() ? "M" : loc.Data()));
    TString out(GetParam("PROOF_OUTPUTFILE"));
    if (!out.IsNull()) fProofFile->SetOutputFileName(out);
    fFile = fProofFile->OpenFile("recreate");
    
  }
  const char* GetParam(const char* key)
  {
    if (!fInput) return 0;
    TObject* o = fInput->FindObject(key);
    if (!o) {
      Warning("GetParam", "Key %s not found in input", key);
      return 0;
    }
    const char* val = static_cast<TNamed*>(o)->GetTitle();
    Info("GetParam", "key=%30s  value=%s", key, val);
    return val;
  }
  Bool_t Process(Long64_t)
  {
    // Printf("Processing");
    if (!fModel) return false;

    // Printf("Make event");
    fModel->Event(fNtracks, fVarTracks, false);

    return true;
  }
  void SlaveTerminate()
  {
    Printf("Slave Terminating");
    if (!fModel) return;
    fModel->Output(fFile);
    fProofFile->Print();
    fFile->Print();
    fOutput->Add(fProofFile);
    fFile->Close();
  }
  
  void Terminate() {}

  static const char* Find(const char* macro, const char* post="+g")
  {
    const char* found = gSystem->Which(gROOT->GetMacroPath(), macro);
    if (!found) return 0;
    
    return Form("%s%s", found, post);
  }
  static void Run(const char* url="workers=10",
		  Int_t       nTracks=500,
		  Int_t       nEvents=1000,
		  Double_t    varTracks=0)
  {
    TString fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/dndeta/tracklets3/toymodel";
    if (gSystem->Getenv("ANA_SRC")) fwd = "$ANA_SRC/dndeta/tracklets3/toymodel";
    gROOT->SetMacroPath(Form("%s:%s",fwd.Data(), gROOT->GetMacroPath()));
    TProof* proof = TProof::Open(url);
    proof->Load(Find("ToyModel.C"));
    proof->Load(Find("ToyModelSelector.C"));
    proof->AddInput(new TNamed("PROOF_OUTPUTFILE",
			       Form("dist_t%06d_e%06d_v%1dd%03d.root",
				    nTracks, nEvents,
				    Int_t(varTracks),
				    Int_t(varTracks*1000) % 1000)));

    ToyModelSelector* selector = new ToyModelSelector;
    selector->fNtracks   = nTracks;
    selector->fVarTracks = varTracks;

    proof->Process(selector, nEvents);
  }
      
  Int_t             fNtracks;
  Double_t          fVarTracks;
  ToyModel::Params  fParams;
  TProofOutputFile* fProofFile; //!
  TFile*            fFile;//!
  ToyModel*         fModel; //!

  ClassDef(ToyModelSelector,1); 
};

//
// EOF
//

