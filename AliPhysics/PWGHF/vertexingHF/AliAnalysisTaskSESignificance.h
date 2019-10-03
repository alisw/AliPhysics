#ifndef ALIANALYSISTASKSESIGNIFICANCE_H
#define ALIANALYSISTASKSESIGNIFICANCE_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
/// \class AliAnalysisTaskSESignificane
/// \brief class to calculate effects on  significance of D mesons  cut
/// \author Authors: G. Ortona, ortona@to.infn.it
/// \author F. Prino, prino@to.infn.it
/// \author Renu Bala, bala@to.infn.it
/// \author Chiara Bianchin, cbianchi@pd.infn.it
//*************************************************************************

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"

class TH1F;
class AliMultiDimVector;
class AliRDHFCuts;

class AliAnalysisTaskSESignificance : public AliAnalysisTaskSE
{
 public:

  enum FeedDownEnum {kBoth,kCharmOnly,kBeautyOnly};
  enum ChanDs {kAllReson,kPhi,kK0star};

  AliAnalysisTaskSESignificance();
  AliAnalysisTaskSESignificance(const char *name, TList *listMDV,AliRDHFCuts *RDCuts, Int_t decaychannel,Int_t selectionlevel=AliRDHFCuts::kAll);
 
  virtual ~AliAnalysisTaskSESignificance();

  Bool_t CheckConsistency();
  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetBFeedDown(FeedDownEnum flagB);//see enum
  void SetDFromCharmOnly(){SetBFeedDown(kCharmOnly);}
  void SetDFromBeautyOnly(){SetBFeedDown(kBeautyOnly);}
  void SetMassLimits(Float_t range,Int_t pdg);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetNBins(Int_t nbins){fNBins=nbins;}
  void SetFillWithPartAntiPartBoth(Int_t value){fPartOrAndAntiPart=value;}
  void SetDsChannel(Int_t chan){fDsChannel=chan;}
  void SetUseSelBit(Bool_t selBit=kTRUE){fUseSelBit=selBit;}
  void SetAODMismatchProtection(Int_t opt=1) {fAODProtection=opt;}

  //void SetMultiVector(const AliMultiDimVector *MultiDimVec){fMultiDimVec->CopyStructure(MultiDimVec);}
  Float_t GetUpperMassLimit()const {return fUpmasslimit;}
  Float_t GetLowerMassLimit()const {return fLowmasslimit;}
  Int_t GetNBins()const {return fNBins;}
  Int_t GetFillWithPartAntiPartBoth()const {return fPartOrAndAntiPart;}
  Int_t GetBFeedDown()const {return fBFeedDown;}
  Int_t GetDsChannel()const {return fDsChannel;}
  Bool_t GetUseSelBit()const {return fUseSelBit;}

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void LocalInit();// {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:

  void SetPDGdaughterDstoKKpi(){
    fPDGdaughters[0]=321;//K
    fPDGdaughters[1]=321;//K
    fPDGdaughters[2]=211;//pi
    fPDGdaughters[3]=0; //empty
  }
  void SetPDGdaughterDstopiKK(){
    fPDGdaughters[0]=211;//pi
    fPDGdaughters[1]=321;//K
    fPDGdaughters[2]=321;//K
    fPDGdaughters[3]=0; //empty
  }

  AliAnalysisTaskSESignificance(const AliAnalysisTaskSESignificance &source);
  AliAnalysisTaskSESignificance& operator=(const AliAnalysisTaskSESignificance& source);
  void SetPDGCodes();
  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*3;}
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*3+1;}
  Int_t GetBackgroundHistoIndex(Int_t iPtBin) const { return iPtBin*3+2;}
  Int_t GetLSHistoIndex(Int_t iPtBin)const { return iPtBin*5;}
  Int_t CheckOrigin(const AliAODMCParticle* mcPart, const TClonesArray* mcArray) const;

  void FillDplus(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Int_t isSel);
  void FillD02p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index, Int_t isSel);
  void FillDs(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Int_t isSel,Int_t optDecay);
  void FillDstar(AliAODRecoCascadeHF* dstarD0pi,TClonesArray *arrayMC,Int_t index,Int_t isSel);
  void FillD04p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Int_t isSel);
  void FillLambdac(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index, Int_t isSel);


  enum {kMaxPtBins=8};
  enum {kMaxCutVar=10};
  enum {kMaxSteps=10};
  enum {kMaxNHist=500000};
  enum {kDplustoKpipi,kD0toKpi,kDstartoKpipi,kDstoKKpi,kD0toKpipipi,kLambdactopKpi};

  TList   *fOutput; //!<! list send on output slot 0
  TList *fCutList; /// Multidimvector container
  TH1F *fMassHist[kMaxNHist]; //!<!hist. for inv mass
  TH1F *fSigHist[kMaxNHist]; //!<!hist. for inv mass (sig from MC truth)
  TH1F *fBkgHist[kMaxNHist]; //!<!hist. for inv mass (bkg from MC truth)
  TH1F *fRflHist[kMaxNHist]; //!<!hist. for inv mass (bkg from MC truth)
  TH1F* fHistNEvents; //!<! hist of n of aods
  Float_t fUpmasslimit;  /// upper inv mass limit for histos
  Float_t fLowmasslimit; /// lower inv mass limit for histos
  AliRDHFCuts *fRDCuts;/// prong cut values
  Int_t fNPtBins; /// number of pt bins
  Int_t fAODProtection;  /// flag to activate protection against AOD-dAOD mismatch.
                         /// -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
  Bool_t fReadMC;    /// flag for access to MC
  Bool_t fUseSelBit;    /// flag to use selection bit (speed up candidates selection)
  FeedDownEnum fBFeedDown; /// flag to search for D from B decays
  Int_t fDecChannel; /// decay channel identifier
  Int_t fPDGmother;  /// PDG code of D meson
  Int_t fNProngs;         /// number of prong of the decay channel
  Int_t fPDGdaughters[4]; /// PDG codes of daughters
  TString fBranchName;    /// AOD branch name for channel
  Int_t fSelectionlevel;  /// selection level: kALL,kTracks,kCandidate
  Int_t   fNVars;         /// number of selection variables
  Float_t fVars[kMaxCutVar];       /// array with values of cut variables
  Int_t fNBins;  /// number of bins in the mass histograms
  Int_t fPartOrAndAntiPart;  /// fill histograms with particle only (+1), antiparticle only (-1), both (0)
  Int_t fDsChannel;          /// Ds resonant channel selected
  Int_t fPDGDStarToD0pi[2]; /// PDG codes for the particles in the D* -> pi + D0 decay
  Int_t fPDGD0ToKpi[2];    /// PDG codes for the particles in the D0 -> K + pi decay

  /// \cond CLASSIMP    
  ClassDef(AliAnalysisTaskSESignificance,6); /// AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
  /// \endcond
};

#endif
