#ifndef ALIANALYSISTASKCHECKGENKINE
#define ALIANALYSISTASKCHECKGENKINE

//*************************************************************************
/// \class Class AliAnalysisTaskCheckGenKine
/// \brief AliAnalysisTask to check MC production at ESD+Kine level
///
///
/// \author Author: F. Prino, prino@to.infn.it
///
//*************************************************************************

class TList;
class TNtuple;
class TH1F;
class TH2F;
class TH3F;
class TTree;
class TString;
class AliESDEvent;
class AliESDfriend;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCheckGenKine : public AliAnalysisTaskSE {

 public:
  
  AliAnalysisTaskCheckGenKine();
  virtual ~AliAnalysisTaskCheckGenKine();

  virtual void   UserExec(Option_t *option);
  virtual void   UserCreateOutputObjects();
  virtual void   Terminate(Option_t *option);

  void SetppConfiguration(){fIsAA=kFALSE;}
  void SetAAConfiguration(){fIsAA=kTRUE;}

  void AddParticleToCheck(Int_t pdg);
  void PrintSpeciesToCheck();

 private:

  AliAnalysisTaskCheckGenKine(const AliAnalysisTaskCheckGenKine &source);
  AliAnalysisTaskCheckGenKine& operator=(const AliAnalysisTaskCheckGenKine &source);
  
  Int_t GetSpeciesIndex(Int_t pdgCode) const{
    for(Int_t j=0; j<fNumOfSpeciesToCheck; j++){if(pdgCode==fPdgCodes[j]) return j;}
    return -1;
  }

  enum ESpeciesCheck {kMaxNumOfSpeciesToCheck=50};  // pi+,pi-,pi0 e+,e-, K+,K-,K0l,K0s, p,pbar  D0, B, d

  TList* fOutput;              //!<! list of output histos
  TH1F*  fHistoNEvents;        //!<! histo with N of events
  TH1F*  fHistoGenMult;        //!<! histo of generated multiplicity in |eta|<0.9

  // vertex histos
  TH3F*  fHistoVtxContrib;     //!<! histo of vertex contributors (SPD, Tracks, multiplicity)
  TH1F*  fHistoSPD3DVtxX;      //!<! histo of vertex coordinates
  TH1F*  fHistoSPD3DVtxY;      //!<! histo of vertex coordinates
  TH1F*  fHistoSPD3DVtxZ;      //!<! histo of vertex coordinates
  TH1F*  fHistoSPDZVtxZ;       //!<! histo of vertex coordinates
  TH1F*  fHistoTrkVtxX;        //!<! histo of vertex coordinates
  TH1F*  fHistoTrkVtxY;        //!<! histo of vertex coordinates
  TH1F*  fHistoTrkVtxZ;        //!<! histo of vertex coordinates
  TH3F*  fHistoSPD3DVtxResidX; //!<! histo of vertex coorindate residuals
  TH3F*  fHistoSPD3DVtxResidY; //!<! histo of vertex coorindate residuals
  TH3F*  fHistoSPD3DVtxResidZ; //!<! histo of vertex coorindate residuals
  TH3F*  fHistoSPDZVtxResidZ;  //!<! histo of vertex coorindate residuals
  TH3F*  fHistoTrkVtxResidX;   //!<! histo of vertex coorindate residuals
  TH3F*  fHistoTrkVtxResidY;   //!<! histo of vertex coorindate residuals
  TH3F*  fHistoTrkVtxResidZ;   //!<! histo of vertex coorindate residuals

  // multiplicity histos
  TH2F* fHistoTracklets;       //!<! histo of trackelts vs. gen mult
  TH2F* fHistoSelTracks;       //!<! histo of tracks (TPC+ITS) vs. gen mult
  TH2F* fHistoGenMultVsb;      //!<! histo of generated multiplicity in |eta|<0.9
  TH2F* fHistoTrackletsVsb;    //!<! histo of trackelts vs. imp par
  TH2F* fHistoSelTracksVsb;    //!<! histo of tracks (TPC+ITS) vs. imp par


  // per-particle histos
  TH2F*  fSpeciesAbundance;                      //!<! histo with abundances for species
  TH2F*  fEtaPt[kMaxNumOfSpeciesToCheck];        //!<! histo with eta,pt per species
  TH3F*  fPrimSec[kMaxNumOfSpeciesToCheck];      //!<! histo with prim/sec, prod radius, pt
  TH2F*  fNumOfDau[kMaxNumOfSpeciesToCheck];     //!<! histo with n daughters per species
  TH2F*  fDecLen[kMaxNumOfSpeciesToCheck];       //!<! histo with decay length per species
  TH2F*  fCt[kMaxNumOfSpeciesToCheck];           //!<! histo with ct per species
  TH2F*  fMassDiff[kMaxNumOfSpeciesToCheck];     //!<! histo of mass mother - inv mass daughters
  TH2F*  fMomDiff[kMaxNumOfSpeciesToCheck];      //!<! histo of mass mother - inv mass daughters
  TH3F*  fPrimSecb[kMaxNumOfSpeciesToCheck];     //!<! histo with prim/sec, prod radius, impact parameter (A-A)
  Int_t  fPdgCodes[kMaxNumOfSpeciesToCheck];     // array of pdg codes
  Bool_t fIsAA;                                  // flag for AA config
  Int_t  fNumOfSpeciesToCheck;                   // actual number of species to be checked
    
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCheckGenKine,4);
  /// \endcond
};


#endif
