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

//#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "FJ_includes.h"

namespace fastjet {
    class PseudoJet;
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

  std::vector<fastjet::PseudoJet> GetBackSubEvent(std::vector <fastjet::PseudoJet> const & full_event, Double_t & rho, Double_t & rho_sparse, Int_t opt = 0);

  void                 FillAllTracks(AliParticleContainer* cont1, AliParticleContainer* cont2, std::vector <fastjet::PseudoJet> & full_event);
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

 private:
  AliAnalysisTaskSDKL(const AliAnalysisTaskSDKL&);            // not implemented
  AliAnalysisTaskSDKL &operator=(const AliAnalysisTaskSDKL&); // not implemented
  ClassDef(AliAnalysisTaskSDKL, 1) // jet sample analysis task

};
#endif
