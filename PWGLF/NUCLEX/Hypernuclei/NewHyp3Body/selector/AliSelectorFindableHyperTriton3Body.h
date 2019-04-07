#ifndef AliSelectorFindableHyperTriton3Body_h
#define AliSelectorFindableHyperTriton3Body_h

#include <TChain.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

// Headers needed by this particular selector
#include "AliESDtrack.h"
#include "AliPIDResponse.h"
#include "AliVertexerHyperTriton3Body.h"
#include "AliVertexerTracks.h"

#include <TVector3.h>
#include <vector>

class TH1D;
class TH2D;

class HyperTritonCandidate {
public:
  HyperTritonCandidate(AliESDtrack *track1, AliESDtrack *track2, AliESDtrack *track3);

  AliESDtrack fDaughterTrack[3];
  int fPDG[3];
  int fMId;
  float fTrueP[3];
  float fTrueDecayVtx[4];
  float fPrimaryVtx[4];
  float fMagField;
  bool fIsGoodCandidate;
  bool fIsClone;
};

class AliSelectorFindableHyperTriton3Body : public TSelector {
public:
  TTreeReader fReader; //! the tree reader
  TTree *fChain = 0;   //! pointer to the analyzed TTree or TChain

  // Readers to access the data (delete the ones you do not need).
  TTreeReaderValue<AliESDtrack> fTreeHyp3BodyVarTracks[3] = {
      {fReader, "fTreeHyp3BodyVarTrack0"}, {fReader, "fTreeHyp3BodyVarTrack1"}, {fReader, "fTreeHyp3BodyVarTrack2"}};
  TTreeReaderValue<Int_t> fTreeHyp3BodyVarPDGcodes[3] = {{fReader, "fTreeHyp3BodyVarPDGcode0"},
                                                         {fReader, "fTreeHyp3BodyVarPDGcode1"},
                                                         {fReader, "fTreeHyp3BodyVarPDGcode2"}};
  TTreeReaderValue<ULong64_t> fTreeHyp3BodyVarEventId = {fReader, "fTreeHyp3BodyVarEventId"};
  TTreeReaderValue<Int_t> fTreeHyp3BodyVarMotherId    = {fReader, "fTreeHyp3BodyVarMotherId"};
  TTreeReaderValue<Float_t> fTreeHyp3BodyVarTrueP[3]  = {
      {fReader, "fTreeHyp3BodyVarTruePx"}, {fReader, "fTreeHyp3BodyVarTruePy"}, {fReader, "fTreeHyp3BodyVarTruePz"}};
  TTreeReaderValue<Float_t> fTreeHyp3BodyVarDecayVtx[4]   = {{fReader, "fTreeHyp3BodyVarDecayVx"},
                                                           {fReader, "fTreeHyp3BodyVarDecayVy"},
                                                           {fReader, "fTreeHyp3BodyVarDecayVz"},
                                                           {fReader, "fTreeHyp3BodyVarDecayT"}};
  TTreeReaderValue<Float_t> fTreeHyp3BodyVarPVtx[4]       = {{fReader, "fTreeHyp3BodyVarPVx"},
                                                       {fReader, "fTreeHyp3BodyVarPVy"},
                                                       {fReader, "fTreeHyp3BodyVarPVz"},
                                                       {fReader, "fTreeHyp3BodyVarPVt"}};
  TTreeReaderValue<Float_t> fTreeHyp3BodyVarMagneticField = {fReader, "fTreeHyp3BodyVarMagneticField"};

  AliSelectorFindableHyperTriton3Body(TString outputName = "output.root", TString outputPath = ".",
                                      TTree * /*tree*/ = 0);
  virtual ~AliSelectorFindableHyperTriton3Body() {}
  AliSelectorFindableHyperTriton3Body(const AliSelectorFindableHyperTriton3Body &) = delete;
  AliSelectorFindableHyperTriton3Body &operator=(const AliSelectorFindableHyperTriton3Body &other) = delete;
  virtual Int_t Version() const { return 2; }
  virtual void Begin(TTree *tree);
  virtual void SlaveBegin(TTree *tree);
  virtual void Init(TTree *tree);
  virtual Bool_t Notify();
  virtual Bool_t Process(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) {
    return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
  }
  virtual void SetOption(const char *option) { fOption = option; }
  virtual void SetObject(TObject *obj) { fObject = obj; }
  virtual void SetInputList(TList *input) { fInput = input; }
  virtual TList *GetOutputList() const { return fOutput; }
  virtual void SlaveTerminate();
  virtual void Terminate();

  bool IsTrackCutsGood(AliESDtrack *lTracks);
  bool IsSameTrack(AliESDtrack *track1, AliESDtrack *track2);

  TString fOutputFilePath;
  TString fOutputFileName;

  AliPIDResponse *fPIDResponse = {nullptr};

  std::vector<HyperTritonCandidate> fClonesVector;

  AliVertexerHyperTriton3Body fHypertritonVertexer;
  AliVertexerHyperTriton3Body fHypertritonVertexerHard;

  TVector3 fPrimaryVertexXYZ;

  ULong_t fCurrentEventId = 0ull;
  int fLastMother         = -1;

  TH2D *fHistInvMass[2][4]    = {{nullptr}};
  TH1D *fHistPt[2][4]         = {{nullptr}};
  TH1D *fHistVertexChi2       = {nullptr};
  TH2D *fHistResDecayVtx[3]   = {nullptr};
  TH1D *fHistCosPAngle        = {nullptr};
  TH2D *fHistDCA2pV[3]        = {nullptr};
  TH2D *fHistDCA2dV[3]        = {nullptr};
  TH1D *fHistTrackDistance[3] = {nullptr};

  TH1D *fHistPDGCheck[4]       = {nullptr};
  TH2D *fHistPDG               = {nullptr};
  TH1D *fHistTrackCheck[4]     = {nullptr};
  TH1D *fHistChargeCheck       = {nullptr};
  TH1D *fHistPassCheck         = {nullptr};
  TH1D *fHistClonesCheck       = {nullptr};
  TH1D *fHistSameTrackCheck[4] = {nullptr};

  TH1D *fHistNSigma[3] = {nullptr};
  TH1D *fHistNSigmaCheck[4] = {nullptr};
  TH1D *fHistNSigmaFakeCheck[4] = {nullptr};

  TH1D *fHistCheckPDGTrackPID = {nullptr};

  ClassDef(AliSelectorFindableHyperTriton3Body, 0);
};

#endif

//______________________________________________________________________________
#ifdef AliSelectorFindableHyperTriton3Body_cxx

void AliSelectorFindableHyperTriton3Body::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // AliSelectorFindableHyperTriton3Body, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);
}

//______________________________________________________________________________
Bool_t AliSelectorFindableHyperTriton3Body::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated AliSelectorFindableHyperTriton3Body, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif // #ifdef AliSelectorFindableHyperTriton3Body_cxx
