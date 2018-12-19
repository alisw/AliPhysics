#ifndef AliSelectorFindableHyperTriton3Body_h
#define AliSelectorFindableHyperTriton3Body_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "AliESDtrack.h"

class TH1D;

class AliSelectorFindableHyperTriton3Body : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<AliESDtrack> fTreeHyp3BodyVarTracks[3] = {
      {fReader, "fTreeHyp3BodyVarTrack0"},
      {fReader, "fTreeHyp3BodyVarTrack1"},
      {fReader, "fTreeHyp3BodyVarTrack2"}
   };
   TTreeReaderValue<Int_t> fTreeHyp3BodyVarPDGcodes[3] = {
      {fReader, "fTreeHyp3BodyVarPDGcode0"},
      {fReader, "fTreeHyp3BodyVarPDGcode1"},
      {fReader, "fTreeHyp3BodyVarPDGcode2"}
   };
   TTreeReaderValue<ULong64_t> fTreeHyp3BodyVarEventId = {fReader, "fTreeHyp3BodyVarEventId"};
   TTreeReaderValue<Int_t> fTreeHyp3BodyVarMotherId = {fReader, "fTreeHyp3BodyVarMotherId"};
   TTreeReaderValue<Float_t> fTreeHyp3BodyVarTruePx = {fReader, "fTreeHyp3BodyVarTruePx"};
   TTreeReaderValue<Float_t> fTreeHyp3BodyVarTruePy = {fReader, "fTreeHyp3BodyVarTruePy"};
   TTreeReaderValue<Float_t> fTreeHyp3BodyVarTruePz = {fReader, "fTreeHyp3BodyVarTruePz"};
   TTreeReaderValue<Float_t> fTreeHyp3BodyVarDecayVx = {fReader, "fTreeHyp3BodyVarDecayVx"};
   TTreeReaderValue<Float_t> fTreeHyp3BodyVarDecayVy = {fReader, "fTreeHyp3BodyVarDecayVy"};
   TTreeReaderValue<Float_t> fTreeHyp3BodyVarDecayVz = {fReader, "fTreeHyp3BodyVarDecayVz"};
   TTreeReaderValue<Float_t> fTreeHyp3BodyVarDecayT = {fReader, "fTreeHyp3BodyVarDecayT"};
   TTreeReaderValue<Float_t> fTreeHyp3BodyVarPVx = {fReader, "fTreeHyp3BodyVarPVx"};
   TTreeReaderValue<Float_t> fTreeHyp3BodyVarPVy = {fReader, "fTreeHyp3BodyVarPVy"};
   TTreeReaderValue<Float_t> fTreeHyp3BodyVarPVz = {fReader, "fTreeHyp3BodyVarPVz"};
   TTreeReaderValue<Float_t> fTreeHyp3BodyVarPVt = {fReader, "fTreeHyp3BodyVarPVt"};
   TTreeReaderValue<Float_t> fTreeHyp3BodyVarMagneticField = {fReader, "fTreeHyp3BodyVarMagneticField"};


   AliSelectorFindableHyperTriton3Body(TTree * /*tree*/ =0) : fReader{} { }
   virtual ~AliSelectorFindableHyperTriton3Body() { }
   AliSelectorFindableHyperTriton3Body(const AliSelectorFindableHyperTriton3Body&) = delete;
   AliSelectorFindableHyperTriton3Body& operator=(const AliSelectorFindableHyperTriton3Body& other) = delete;
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();


   ULong_t fCurrentEventId = 0ull;
   int fLastMother = -1;
   TH1D* fHistInvMass = nullptr;
   TH1D* fHistClones = nullptr;
   ClassDef(AliSelectorFindableHyperTriton3Body,0);

};

#endif

#ifdef AliSelectorFindableHyperTriton3Body_cxx
void AliSelectorFindableHyperTriton3Body::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // AliSelectorFindableHyperTriton3Body, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t AliSelectorFindableHyperTriton3Body::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated AliSelectorFindableHyperTriton3Body, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef AliSelectorFindableHyperTriton3Body_cxx
