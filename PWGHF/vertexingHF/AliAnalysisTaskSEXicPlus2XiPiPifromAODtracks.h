#ifndef ALIANALYSISTASKSEXICPLUS2XIPIPIFROMAODTRACKS_H
#define ALIANALYSISTASKSEXICPLUS2XIPIPIFROMAODTRACKS_H

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */ 

#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliPID.h"
#include "AliRDHFCutsXicPlustoXiPiPifromAODtracks.h"

/// \class AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks

class THnSparse;
class TH1F;
class TH2F;
class TH3F;
class TClonesArray;
class AliAODRecoCascadeHF3Prong;
class AliAODPidHF;
class AliESDtrackCuts;
class AliESDVertex;
class AliAODMCParticle;
class AliNormalizationCounter;

class AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks : public AliAnalysisTaskSE 
{
 public:
  enum ECandStatus {kGenLimAcc,kGenAccMother,kGenAccMother08,kGenAcc,kGenAcc08,kReco,kReco08,kRecoCuts,kRecoCuts08,kRecoPID,kRecoPID08};
  
  AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks();
  AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks(const Char_t* name, AliRDHFCutsXicPlustoXiPiPifromAODtracks* cuts, Bool_t writeVariableTree=kTRUE, Bool_t fillSparse=kFALSE);
  virtual ~AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks();

  /// Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void FillROOTObjects(AliAODRecoCascadeHF3Prong *xicobj, AliAODMCParticle *mcpart, AliAODMCParticle *mcdau1, AliAODMCParticle *mcdau2, AliAODMCParticle *mcdauxi, Int_t mcnused, Bool_t isXiC, Int_t checkOrigin);
  void MakeAnalysis(AliAODEvent *aod, TClonesArray *mcArray);

  
  /// set MC usage
  void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t GetMC() const {return fUseMCInfo;}
  void SetFillSignalOnly(Bool_t signalOnly) {fFillSignalOnly = signalOnly;}
  Bool_t GetFillSignalOnly() const {return fFillSignalOnly;}
  void SetFillBkgOnly(Bool_t bkgOnly) {fFillBkgOnly = bkgOnly;}
  Bool_t GetFillBkgOnly() const {return fFillBkgOnly;}
  
  void   SetReconstructPrimVert(Bool_t a) { fReconstructPrimVert=a; }
  void SelectCascade( const AliVEvent *event,Int_t nCascades,Int_t &nSeleCasc, Bool_t *seleCascFlags);
  void SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags);
  Bool_t SelectLikeSign(AliAODTrack *trk1, AliAODTrack *trk2);
  AliAODRecoCascadeHF3Prong* MakeCascadeHF3Prong(AliAODcascade *casc, AliAODTrack *trk1, AliAODTrack *trk2, AliAODEvent *aod, AliAODVertex *secvert, Double_t dispersion);

  void LoopOverGenParticles(TClonesArray *mcArray);
  Int_t CheckXic2XiPiPi(TClonesArray* arrayMC, AliAODMCParticle *mcPart, Int_t* arrayDauLab);

  void SetITS3UpgradeAnalysis(Bool_t isITS3Upgrade) {fIsXicPlusUpgradeITS3=isITS3Upgrade;}
  
 private:
  
  AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks(const AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks &source);
  AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks& operator=(const AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks& source); 
  
  void DefineTreeVariables();
  void DefineGeneralHistograms();
  void DefineAnalysisHistograms();

  AliAODVertex *CallPrimaryVertex(AliAODcascade *casc, AliAODTrack *trk1, AliAODTrack *trk2, AliAODEvent *evt);
  AliAODVertex* PrimaryVertex(const TObjArray *trkArray,AliVEvent *event);
  AliAODVertex* CallReconstructSecondaryVertex(AliAODTrack *trk1, AliAODTrack *trk2,Double_t &disp);
  AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray, Double_t &dispersion,Bool_t useTRefArray=kTRUE);
  
  Bool_t fUseMCInfo;          /// Use MC info
  Bool_t fFillSignalOnly;     /// Fill only the signal in MC
  Bool_t fFillBkgOnly;        /// Fill only the bkg in MC
  TList *fOutput;             //!<! User output slot 1 // general histos
  TList *fOutputAll;          //!<! User output slot 3 // Analysis histos
  TList *fListCuts;           //!<! User output slot 2 // Cuts 
  TH1F *fCEvents;             /// Histogram to check selected events
  TH1F *fHTrigger;            /// Histograms to check trigger
  TH1F *fHCentrality;         /// histogram to check centrality
  AliRDHFCutsXicPlustoXiPiPifromAODtracks *fAnalCuts;      /// Cuts - sent to output slot 2
  Bool_t fIsEventSelected;           /// flag for event selected
  Bool_t    fWriteVariableTree;       /// flag to decide whether to write the candidate variables on a tree variables
  Bool_t fFillSparse;                 /// flag to decide whether fill the THnSparse
  TTree    *fVariablesTree;           //!<! tree of the candidate variables after track selection on output slot 4
  Bool_t fReconstructPrimVert;            /// Reconstruct primary vertex excluding candidate tracks
  Bool_t fIsMB;            /// Is MB event
  Bool_t fIsSemi;          /// is semi-central trigger event
  Bool_t fIsCent;          /// is central trigger event
  Bool_t fIsINT7;          /// is int7 trigger event
  Bool_t fIsEMC7;          /// is emc7 trigger event
  Float_t *fCandidateVariables;     //!<! variables to be written to the tree
  AliAODVertex *fVtx1;              /// primary vertex
  AliESDVertex *fV1;                /// primary vertex
  Double_t fBzkG;                   /// magnetic field value [kG]
  Float_t  fCentrality;             /// centrality
  //Float_t  fTriggerCheck;           /// Trigger information
  
  //--------------------- My histograms ------------------
  THnSparse*  fHistoXicMass;        //!<! xic mass spectra
  THnSparse*  fSparseXicMass;       //!<! xic sparse to study cut variation
  
  TH3F*  fHistoMCSpectrumAccXic;    //!<! Spectrum of generated particles
  TH1F*  fHistoDcaPi1Pi2;                    //!<!  DCA between pions
  TH1F*  fHistoDcaPi1Casc;                    //!<! DCA between pi and cascade
  TH1F*  fHistoDcaPi2Casc;                    //!<! DCA between pi and cascade
  TH1F*  fHistoLikeDecayLength;              //!<! Decay length
  TH1F*  fHistoLikeDecayLengthXY;            //!<! Decay length in XY
  TH1F*  fHistoXicCosPA;                   //!<! Xic cosine pointing angle
  
  TH1F*  fHistoXiMass;                       //!<! mass of xi
  TH1F*  fHistoCascDcaXiDaughters;           //!<! DCA of xi daughgers
  TH1F*  fHistoCascDcaV0Daughters;           //!<! DCA of v0 daughters
  TH1F*  fHistoCascDcaV0ToPrimVertex;        //!<! DCA of v0 to primary vertex 
  TH1F*  fHistoCascDcaPosToPrimVertex;       //!<! DCA of positive track to primary vertex 
  TH1F*  fHistoCascDcaNegToPrimVertex;       //!<! DCA of negative track to primary vertex 
  TH1F*  fHistoCascDcaBachToPrimVertex;      //!<! DCA of bachelor track to primary vertex 
  TH1F*  fHistoCascCosPAXiPrim;              //!<! Cosine pointing angle of Xi to primary vertex
  TH1F*  fHistoXiPt;                         //!<! Xi pt
  
  TH1F*  fHistoPiPt;                         //!<! Pion pT
  TH1F*  fHistoPid0;                         //!<! pion d0
  TH1F*  fHistonSigmaTPCpi;                  //!<! nSigma of TPC pion
  TH1F*  fHistonSigmaTOFpi;                  //!<! nSigma of TOF pion
  TH1F*  fHistoProbPion;                     //!<! Probability to be pion

  TH2F*  fHistoXiMassvsPtRef1;               //!<! Reference Xi mass spectra 
  TH2F*  fHistoXiMassvsPtRef2;               //!<! Reference Xi mass spectra 
  TH2F*  fHistoXiMassvsPtRef3;               //!<! Reference Xi mass spectra 
  TH2F*  fHistoXiMassvsPtRef4;               //!<! Reference Xi mass spectra 
  TH2F*  fHistoXiMassvsPtRef5;               //!<! Reference Xi mass spectra 
  TH2F*  fHistoXiMassvsPtRef6;               //!<! Reference Xi mass spectra 
  TH1F*  fHistoPiPtRef;                      //!<! Reference pi spectra
  TH1F*  fHistoPiEtaRef;                      //!<! Reference eta distribution of pi 
  TH1F*  fQAHistoNSelectedTracks;            //!<! QA histo for number of selected tracks/event
  TH1F*  fQAHistoNSelectedCasc;              //!<! QA histo for number of selected Cascades/event
  TH1F*  fQAHistoDCApi1pi2;                  //!<! QA histo for dca betwen two pions from XiC
  TH1F*  fQAHistoAODPrimVertX;               //!<! Coordinates of the primary vertex
  TH1F*  fQAHistoAODPrimVertY;               //!<! Coordinates of the primary vertex
  TH1F*  fQAHistoAODPrimVertZ;               //!<! Coordinates of the primary vertex
  TH1F*  fQAHistoRecoPrimVertX;              //!<! Coordinates of the reconstructed primary vertex without XiC decay tracks
  TH1F*  fQAHistoRecoPrimVertY;              //!<! Coordinates of the reconstructed primary vertex without XiC decay tracks
  TH1F*  fQAHistoRecoPrimVertZ;              //!<! Coordinates of the reconstructed primary vertex without XiC decay tracks
  TH1F*  fQAHistoSecondaryVertexX;           //!<! Coordinates of the reconstructed secondary vertex
  TH1F*  fQAHistoSecondaryVertexY;           //!<! Coordinates of the reconstructed secondary vertex
  TH1F*  fQAHistoSecondaryVertexZ;           //!<! Coordinates of the reconstructed secondary vertex
  TH1F*  fQAHistoSecondaryVertexXY;          //!<! Coordinates of the reconstructed secondary vertex
  AliNormalizationCounter *fCounter;         //!<!Counter for normalization slot 4
  Bool_t fIsXicPlusUpgradeITS3;              ///flag to identify if the analysis is for the ITS3 upgrade
  
  /// \cond CLASSIMP    
  ClassDef(AliAnalysisTaskSEXicPlus2XiPiPifromAODtracks,10); /// class for Xic->Xipipi
  /// \endcond
};
#endif

