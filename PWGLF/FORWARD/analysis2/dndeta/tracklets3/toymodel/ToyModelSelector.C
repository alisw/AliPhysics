#include <TSelector.h>
#include "ToyModel.C"
#ifndef __CINT__
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TFile.h>
#else
class TTree;
class TFile;
#endif


struct ToyModelSelector : public TSelector
{

  ToyModelSelector()
    : fModel(0),
      fFile(0) 
  {}
  Int_t Version() const { return 2; }
  void Init(TTree*) { Printf("Initialize"); }
  void Begin(TTree*) { Printf("Begin it"); }
  void SlaveBegin(TTree*)
  {
    Printf("Slave begining");
    fModel = new ToyModel(fParams);

    fFile = TFile::Open(fOutputName, "RECREATE");

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
    fFile->Print();
    fFile->Close();
  }
  
  void Terminate() {}

  static const char* Find(const char* macro, const char* post="+g")
  {
    const char* found = gSystem->Which(gROOT->GetMacroPath(), macro);
    if (!found) return 0;
    
    return Form("%s%s", found, post);
  }
  static void Run(const char* output="",
		  Int_t       nTracks=500,
		  Int_t       nEvents=1000,
		  Double_t    varTracks=0)
  {
    TString fwd = "$ALICE_PHYSICS/PWGLF/FORWARD/dndeta/tracklets3/toymodel";
    if (gSystem->Getenv("ANA_SRC")) fwd = "$ANA_SRC/dndeta/tracklets3/toymodel";
    gROOT->SetMacroPath(Form("%s:%s",fwd.Data(), gROOT->GetMacroPath()));
    ToyModelSelector* selector = new ToyModelSelector;
    selector->fNtracks   = nTracks;
    selector->fVarTracks = varTracks;

    selector->fOutputName = output;
    if (selector->fOutputName.IsNull())
      selector->fOutputName.Form("dist_t%06d_e%06d_v%1dd%03d.root",
                                nTracks, nEvents, Int_t(varTracks),
                                Int_t(varTracks*1000) % 1000);
    selector->Begin(0);
    selector->SlaveBegin(0);
    for (Int_t i = 0; i < nEvents; ++i) selector->Process(i);
    selector->SlaveTerminate();
    selector->Terminate();
    delete selector;
  }
      
  TString           fOutputName;
  Int_t             fNtracks;
  Double_t          fVarTracks;
  ToyModel::Params  fParams;
  TFile*            fFile;//!
  ToyModel*         fModel; //!

  ClassDef(ToyModelSelector,1); 
};

//
// EOF
//

