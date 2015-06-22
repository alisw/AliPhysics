#ifndef ALIANALYSISTASKEMCALJETSAMPLE_H
#define ALIANALYSISTASKEMCALJETSAMPLE_H

// $Id$

class TH1;
class TH2;
class TH3;
class TClonesArray;
class TArrayF;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

#include "TObject.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisHFjetTagHFE : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisHFjetTagHFE();
  AliAnalysisHFjetTagHFE(const char *name);
  virtual ~AliAnalysisHFjetTagHFE();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  void                        CheckClusTrackMatching();


  AliVEvent   *fVevent;  //!event object
  //TClonesArray  *ftrack;
  //TClonesArray  *fCaloClusters;
  AliPIDResponse *fpidResponse; //!pid response
  // General histograms
  TH1                       **fHistTracksPt;            //!Track pt spectrum
  TH1                       **fHistClustersPt;          //!Cluster pt spectrum
  TH1                       **fHistLeadingJetPt;        //!Leading jet pt spectrum
  TH2                       **fHistJetsPhiEta;          //!Phi-Eta distribution of jets
  TH2                       **fHistJetsPtArea;          //!Jet pt vs. area
  TH2                       **fHistJetsPtLeadHad;       //!Jet pt vs. leading hadron
  TH2                       **fHistJetsCorrPtArea;      //!Jet pt - bkg vs. area
  TH3                        *fHistPtDEtaDPhiTrackClus; //!track pt, delta eta, delta phi to matched cluster
  TH3                        *fHistPtDEtaDPhiClusTrack; //!cluster pt, delta eta, delta phi to matched track

  TH1                        *fHistClustDx; //!
  TH1                        *fHistClustDz; //!

  TH2                        *fHistTPCnSigma;
  TH2                        *fHistEop;
  TH1                        *fHistIncEle;
  TH2                        *fHistIncjet;
  TH2                        *fHistHFjet;
  TH2                        *fInvmassULS;
  TH2                        *fInvmassLS;

  AliJetContainer            *fJetsCont;                   //!Jets
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters  
  Bool_t tagHFjet(AliEmcalJet* jet, double *epT, int MCpid, double &maxpT_e);
  //void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec);
  void SelectPhotonicElectron(Int_t itrack, TClonesArray *ftrackCL, AliVTrack *track, Bool_t &fFlagPhotonicElec);

 private:

  //TClonesArray  *ftrack;

  AliAnalysisHFjetTagHFE(const AliAnalysisHFjetTagHFE&);            // not implemented
  AliAnalysisHFjetTagHFE &operator=(const AliAnalysisHFjetTagHFE&); // not implemented

  ClassDef(AliAnalysisHFjetTagHFE, 4) // jet sample analysis task
};
#endif
