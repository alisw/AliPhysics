#ifndef AliAnalysisTaskSDKLResponse_H
#define AliAnalysisTaskSDKLResponse_H

// $Id$

class TH1;
class TH2;
class TH3;
class TNtuple;
class THnSparse;
class TRandom;
class AliJetContainer;
class AliParticleContainer;
//class AliClusterContainer;

//#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
//#include "AliFJWrapper.h"
#include "FJ_includes.h"
//#include "AliEmcalJet.h"
//#include "AliJetContainer.h"
#include "AliAnalysisTaskSDKL.h"

struct mjet {
  double pt;
  double pt_scalar;
  double eta;
  double phi;
  double area;
  std::vector<double> track_pts;
  std::vector<int>    track_labels;
  std::vector<split>  splits;
  AliEmcalJet        *pointerAJet;
  fastjet::PseudoJet *pointerPJet;
};

class AliAnalysisTaskSDKLResponse : public AliAnalysisTaskSDKL {
 public:

  AliAnalysisTaskSDKLResponse();
  AliAnalysisTaskSDKLResponse(const char *name, Int_t const backgroption, Double_t const fractioneventsfortree);
  virtual ~AliAnalysisTaskSDKLResponse();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  static AliAnalysisTaskSDKLResponse* AddTaskSoftDropResponse(
    const char *ntracks            = "usedefault",
    const char *njets1             = "Jets1",
    const char *njets2             = "Jets2",
    const char *nrho               = "Rho",
    Int_t       nCentBins          = 1,
    Double_t    jetradius          = 0.4,
    Double_t    jetptcut           = 1,
    Double_t    jetareacut         = 0.6,
    const char *type               = "EMCAL",
    Int_t       backgroption       = 0,
    Int_t       leadhadtype        = 0,
    Double_t    fractioneventsfortree = 1.e-6,
    const char *taskname           = "AliAnalysisTaskSDKLResponse"
  );

  void SetCustomMCEfficiency(Double_t mc_eff) { fMCTrackEfficiency = mc_eff; }

 protected:

  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;

  THnSparse                  *fhResponse[4];                //!<! distribution of all
  TH2F                       *fhPtDeltaPt;                  //!<!
  TH2F                       *fhPtDeltaZg[4];               //!<!
  TH2F                       *fhPtDeltaRg[4];               //!<!
  TH2F                       *fhPtDeltaMg[4];               //!<!

  THnSparse                  *fhResponseBackSub[4];                //!<! distribution of all
  TH2F                       *fhPtDeltaPtBackSub;                  //!<!
  TH2F                       *fhPtDeltaZgBackSub[4];               //!<!
  TH2F                       *fhPtDeltaRgBackSub[4];               //!<!
  TH2F                       *fhPtDeltaMgBackSub[4];               //!<!

  THnSparse                  *fhResponseDet[4];                //!<! distribution of all
  TH2F                       *fhPtDeltaPtDet;                  //!<!
  TH2F                       *fhPtDeltaZgDet[4];               //!<!
  TH2F                       *fhPtDeltaRgDet[4];               //!<!
  TH2F                       *fhPtDeltaMgDet[4];               //!<!

//  TH1F                       *fhRho;                          //!<!
//  TH1F                       *fhRhoSparse;                    //!<!

  TH2F                       *fhPtDeltaPtAreaBackSub;          //!<!
  TH2F                       *fhPtDeltaPtAreaBackSubSparse;    //!<!

  TH2F                       *fhEshareDistPlDl; //!<!
  TH2F                       *fhEshareDistDlDluebs; //!<!
  TH2F                       *fhPtProbePtResp; //!<!

  TNtuple                    *fTreeDL;                      //!<!
  TNtuple                    *fTreeDLUEBS;                  //!<!
  TNtuple                    *fTreePL;                      //!<!

  AliJetContainer            *fJetsCont1;                   //! Jets
  AliJetContainer            *fJetsCont2;                   //! Jets
  AliParticleContainer       *fTracksCont1;                 //! Tracks
  AliParticleContainer       *fTracksCont2;                 //! Tracks

 private:
  AliAnalysisTaskSDKLResponse(const AliAnalysisTaskSDKLResponse&);            // not implemented
  AliAnalysisTaskSDKLResponse &operator=(const AliAnalysisTaskSDKLResponse&); // not implemented

  int FillResponseFromMjets( mjet const & mjet1, mjet const & mjet2, THnSparse* hr[4], Float_t const dist, Float_t const eshare, Bool_t const iscl);
  int FillDeltasFromMjets( mjet const & mjet1, mjet const & mjet2, TH2F *hptdpt, TH2F *h1[4], TH2F *h2[4], TH2F *h3[4] );

  Float_t CalcDist(mjet const & mjet1, mjet const & mjet2);
  Float_t CalcEnergyShare(mjet const & mjet1, mjet const & mjet2);
  Bool_t IsClosestPair(int const i, int const j, std::vector <mjet> const & mjet_container1, std::vector <mjet> const & mjet_container2);

  void FillMjetContainer(AliJetContainer *jet_container, std::vector <mjet> & mjet_container);
  void FillMjetContainer(std::vector <fastjet::PseudoJet> & jet_container, std::vector <mjet> & mjet_container);

  void FillRespTree(std::vector<AliEmcalJet*> const & probe_jets, std::vector<fastjet::PseudoJet> const & resp_jets, TNtuple* tree);

  Double_t                   fMCTrackEfficiency;

  Double_t                   fFractionEventsDumpedToTree;
  TRandom                    *fRandom; //!<!


  ClassDef(AliAnalysisTaskSDKLResponse, 1) // jet sample analysis task
};
#endif
