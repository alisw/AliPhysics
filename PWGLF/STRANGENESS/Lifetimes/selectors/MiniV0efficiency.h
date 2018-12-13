#ifndef MiniV0efficiency_h
#define MiniV0efficiency_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TF1.h>
#include "TRandom3.h"
#include <vector>

/// Headers needed by this particular selector
#include "MiniV0.h"
#include "MCparticle.h"
#include "HyperTriton2Body.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TGraph.h>


class MiniV0efficiency : public TSelector {

public:
  TTreeReader     fReader;      //! the tree reader
  TTree          *fChain = 0;   //! pointer to the analyzed TTree or TChain

  /// Readers to access the data
  TTreeReaderValue<float>  fMultiplicity = {fReader, "fMultiplicity"};
  TTreeReaderArray<Lifetimes::MiniV0> V0s = {fReader, "V0s"};
  TTreeReaderArray<Lifetimes::HyperTriton2Body> V0Hyper = {fReader, "V0Hyper"};  
  TTreeReaderArray<Lifetimes::MCparticle> MCparticles = {fReader, "MCparticles"};    

  /// standard selector methods
  MiniV0efficiency(TTree * /*tree*/ =0);
  virtual ~MiniV0efficiency() { }
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

  void BinLogAxis(TH2 *h);

  /// output file name
  string fOutputFileName;


  /// Control histograms to monitor the filtering

  TH1D* fHistV0ptMC[3];
  TH1D* fHistV0ptData[3];
  TH1D* fHistV0ctMC[3];
  TH1D* fHistV0ctData[3]; 
  TH1D* EffvsPt[3];
  TH1D* Effvsct[3];
  TH2D* ctAnalysis[3];
  TH2D* ptAnalysisH;  


  ClassDef(MiniV0efficiency,2);

};

#endif
