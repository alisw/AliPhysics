#ifndef AliAnalysisTaskSDKL_H
#define AliAnalysisTaskSDKL_H

// $Id$

class TH1;
class TH2;
class TH3;
class THnSparse;
//class TTree;
class TNtuple;
class AliJetContainer;
class AliParticleContainer;
class TRandom;

//#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "FJ_includes.h"

namespace fastjet {
  class PseudoJet;
  class ClusterSequenceArea;
  namespace contrib {
    class ConstituentSubtractor;
  }
}

struct split {
    float z;
    float r;
    float m;
    int sd_step;
};

class AliAnalysisTaskSDKL : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskSDKL(const char *name = "AliAnalysisTaskSDKL") ;
  AliAnalysisTaskSDKL(const char *name, Int_t const backgroption);
  virtual ~AliAnalysisTaskSDKL();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  static AliAnalysisTaskSDKL* AddTaskSoftDrop(
    const char *ntracks            = "usedefault",
    const char *njets              = "Jets",
    const char *nrho               = "Rho",
    Int_t       nCentBins          = 1,
    Double_t    jetradius          = 0.4,
    Double_t    jetptcut           = 1,
    Double_t    jetareacut         = 0.6,
    const char *type               = "EMCAL",
    Int_t       backgroption       = 0,
    Int_t       leadhadtype        = 0,
    const char *taskname           = "AliAnalysisTaskSDKL"
  );

  std::vector<split>   ReclusterFindHardSplits(AliEmcalJet *jet);
  std::vector<split>   ReclusterFindHardSplits(fastjet::PseudoJet const & jet);
  std::vector<split>   ReclusterFindHardSplits(std::vector <fastjet::PseudoJet> const & particles);
  std::vector<split>   FindHardSplits(fastjet::PseudoJet const & jet);

  void FilterJets(std::vector<fastjet::PseudoJet> const & jets, std::vector<fastjet::PseudoJet> & jets_filtered, Float_t pt_cut);

  int InitializeSubtractor(std::vector <fastjet::PseudoJet> const & event_full, Double_t & rho, Double_t & rho_sparse, Int_t opt = 0);

  void SetFullEventConstSubtractionMode() { fCSOption = 0; }
  void SetJetByJetConstSubtractionMode()  { fCSOption = 1; }
  void SetAlphaDeltaRmaxConstSubtraction(Double_t alpha, Double_t drmax) {
    fCSAlpha = alpha;
    fCSDeltaRmax = drmax;
  }

  std::vector<fastjet::PseudoJet> GetBackSubJets(std::vector<fastjet::PseudoJet> const & event_full);

  void                 AddTracksToEvent(AliParticleContainer* cont, std::vector <fastjet::PseudoJet> & event, Double_t const efficiency = 10.0, TRandom* fRandom = nullptr);
  void                 FillTree(std::vector<fastjet::PseudoJet> const & jets, TNtuple* tree);
  void                 FillTree(AliJetContainer *jets, TNtuple* tree);

  void                 FillSparseFromSplits(THnSparse *histo, std::vector<split> const & splits, double const jet_pt);

 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;

  THnSparse                  *fhAll;                   //!<! distribution of all
  THnSparse                  *fhAllBackSub;            //!<! distribution of all

  TH1F                       *fhRho;                   //!<!
  TH1F                       *fhRhoSparse;             //!<!

  TNtuple                    *fTree;                   //!<!
  TNtuple                    *fTreeBackSub;            //!<!

  AliJetContainer            *fJetsCont;               //! Jets
  AliParticleContainer       *fTracksCont;             //! Tracks
  Int_t                       fbcoption;
  Int_t                       fCSOption;
  Double_t                    fCSAlpha;
  Double_t                    fCSDeltaRmax;
  fastjet::contrib::ConstituentSubtractor *fCSubtractor; //!
  fastjet::ClusterSequenceArea            *fCSubtractorCS; //!

 private:
  AliAnalysisTaskSDKL(const AliAnalysisTaskSDKL&);            // not implemented
  AliAnalysisTaskSDKL &operator=(const AliAnalysisTaskSDKL&); // not implemented


  ClassDef(AliAnalysisTaskSDKL, 1) // jet sample analysis task

};
#endif
