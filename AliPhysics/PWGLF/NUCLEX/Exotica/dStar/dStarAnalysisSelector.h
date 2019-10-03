//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan 22 15:13:47 2018 by ROOT version 6.10/08
// from TTree AliCODEX/Alice COmpressed Dataset for EXotica
// found on file: AliCODEXTree.root
//
// NEEDS ROOT6!!
//
// Selector used to analyze Tree generated with CODEX and provide
// blinded invariant mass distribution or d* candidates
// ([2.280; 2.480] GeV/c^2 blinded region) with different
// cuts on π+ π- invarant mass
//
// (for more see "/AliPhysics/PWGLF/NUCLEX/Utils/CODEX")
// Species convention in AliAnalysisCODEX
//////////////////////////////////////////////////////////

#ifndef dStarAnalysisSelector_h
#define dStarAnalysisSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "TRandom3.h"
#include <vector>

/// Headers needed by this particular selector
#include "AliAnalysisCODEX.h"
#include "TH2.h"


typedef AliAnalysisCODEX::Track Track;
typedef AliAnalysisCODEX::Header Header;
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> FourVector_t;


class dStarAnalysisSelector : public TSelector {

public:
  TTreeReader     fReader;      //! the tree reader
  TTree          *fChain = 0;   //! pointer to the analyzed TTree or TChain

  /// Readers to access the data
  TTreeReaderValue<Short_t> fZvert = {fReader, "mZvert"};
  TTreeReaderArray<Track> fTracks = {fReader, "Tracks"};

  /// old reader methods needed to connect AliAnalysisCODEX::Header to TOFpidLite correctly
  TBranch *header_b;
  Header  *header;

  /// standard selector methods
  dStarAnalysisSelector(TTree * /*tree*/ =0);
  virtual ~dStarAnalysisSelector() { }
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

  /// user defined methods
  bool PionCutGood(const Track& tr, float nSigmaTPC);
  bool DeuteronCutGood(const Track& tr, float nSigmaTPC);
  FourVector_t TrackToFourVector(const Track& tr, int species);
  void BinLogAxis(TH2 *h);

  /// output file name
  string fOutputFileName;

  /// TOFpidLite object from AliCODEX for particle identification with TOF
  AliAnalysisCODEX::TOFpidLite fTOFpid;

  /// support vector instantiation
  vector<FourVector_t> fDeuteron;
  vector<FourVector_t> fPiPlus;
  vector<FourVector_t> fPiMinus;

  /// nSigma definition for Pions and Deuterons
  float fNSigmaPi;
  float fNSigmaDe;

  ///  control histograms
  TH1F *fZVertex;
  TH1F *fPiPlusPT;
  TH1F *fPiMinusPT;
  TH1F *fDeuPT;
  TH1F *fMultiplicity;

  /// pions energy loss in TPC for check
  TH2D *fDedxPiDeu;
  TH2D *fDedx;
  TH2D *fBeta;

  /// invariant mass distribution with blinded region (2.280 < mInv < 2.480) GeV/c^2
  TH2D *fMInvBlind[11];
  /// event mixing background with cuts (last one without cuts)
  TH2D *fMEMInvBackPPSame[11];
  TH2D *fMEMInvBack3Event[11];
  /// like-sign background
  TH2D *fLSPlusMInvBackground[11];
  TH2D *fLSMinusMInvBackground[11];

  /// event mixing pool
  AliAnalysisCODEX::EventMixingPool<1, 10, 3> *fPool;


  ClassDef(dStarAnalysisSelector,0);

};

#endif
