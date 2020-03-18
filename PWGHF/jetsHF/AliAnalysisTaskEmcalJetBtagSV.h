#ifndef AliAnalysisTaskEmcalJetBtagSV_H
#define AliAnalysisTaskEmcalJetBtagSV_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* Class AliAnalysisTaskEmcalJetBtagSV:                                    *
 * AliAnalysisTaskSE for the extraction of the B-jet spectrum        *
 * using the properties of secondary verices as tagging observables. */

/*
 Mailto: andrea.rossi@cern.ch, elena.bruna@to.infn.it,
 svallero@to.infn.it, s.lapointe@cern.ch
 ycorrale@cern.ch
 */

//--Root--
class TH1F;
class TList;
class TRandom3;
class TRandom;
class THn;
class THnSparse;
class TProfile;
class AliAODMCHeader;
class AliMultSelection;
class AliFJWrapper; //EMB_clus

//--AliRoot--

#include "AliLog.h"  
#include "AliAnalysisTaskEmcalJet.h"
//#include "AliFJWrapper.h"  //EMB_clus

class AliAnalysisUtils;
class AliParticleContainer;

//--AliHFJetsClass--
#include "AliHFJetsUtils.h"
#include "AliRDHFJetsCutsVertex.h"
#include "AliHFJetsTaggingVertex.h"

class AliHFJetsContainerVertex;

//-------------------------------------------------------------------------------------

class AliAnalysisTaskEmcalJetBtagSV : public AliAnalysisTaskEmcalJet {

public:

  AliAnalysisTaskEmcalJetBtagSV();

  AliAnalysisTaskEmcalJetBtagSV(const char* name);

  virtual ~AliAnalysisTaskEmcalJetBtagSV();

  // Setters/Getters
  void SetCorrectionMode(Bool_t mode ) {fCorrMode = mode;};

  void SetDoBkgRejection(Bool_t dorej) {fDoBkgRej = dorej;};

  void SetDoRandomCone(Bool_t dorndm) {fDoRndmCone = dorndm;};

  void SetDoFillSecVtxQA(Bool_t doqa ) {fDoQAVtx  = doqa;};

  void SetDoFillV0Trks(Bool_t doV0) {fDoFillV0Trks = doV0;};

  void SetDoDetRespMtx(Bool_t domtx) {fDoDetRespMtx = domtx;};

  void SetDoOnlyMtxAna(Bool_t onlyMtx) {fDoOnlyMtxAna = onlyMtx;};

  void SetEmbeddPerpendicular(Bool_t EmbeddPerpendicular) {fEmbeddPerpendicular = EmbeddPerpendicular;};  //EMB_clus

  

  void SetRecJetsBranch(const char* branch) {fRecJetsBranch = branch;}

  void SetGenJetsBranch(const char* branch) {fGenJetsBranch = branch;}

  void SetGenNamePattern(const char* patt)  {fGenNamePattern = patt; }

  void SetPtHardName(const char* pthard) { fPtHardName = pthard; };

  void SetCheckMCCrossSection(Bool_t check) {fCheckMCCrossSection = check;}

  void SetUseWeightOn() {fUseWeight = kTRUE;}

  void SetJetTaggingRadius(Double_t tagradius) {fTaggingRadius = tagradius;};
  void SetDelPtSigmaCut(Double_t SigmaCut) { fSigmaSVCut = SigmaCut; }            //newDeltaPt//


  void SetGlLogLevel(Int_t level) { fGlLogLevel = level; };
  void SetLcDebLevel(Int_t level) { fLcDebLevel = level; };

  void SetTagger(AliHFJetsTaggingVertex* tagger)
  {
    if (fTagger) delete fTagger;
    fTagger = (AliHFJetsTaggingVertex*)tagger->Clone("fTagger");
  }

  void SetCuts(AliRDHFJetsCuts* cuts)
  {
    if (fCutsHFjets) delete fCutsHFjets;
    fCutsHFjets = (AliRDHFJetsCuts*)cuts->Clone("fCutsHFjets");
  }

  void SetJetContName(const char* name)        { fJetContName = name; }

  void SetTrkContName(const char* name)        { fTrkContName = name; }

  void SetRhoTaskName(const char* name)        { fRhoTaskName = name; }

  void SetMCJetContName(const char* name)      { fMCJetContName = name; }

  void SetMCTrkContName(const char* name)      { fMCTrkContName = name; }

  void SetMCRhoTaskName(const char* name)      { fMCRhoTaskName = name; }

  void SetAnalysisUtils(AliAnalysisUtils* utils) { fAnalysisUtils = utils; }

  virtual void     GetPythiaCrossSection();

  static  Bool_t   GetPythiaInfoFromFile(TString currFile,
                                         Float_t& xsec,
                                         Float_t& trials);

  AliAnalysisUtils* GetAnalysisUtils()           { return fAnalysisUtils; }

  void SetMaxFacPtHard(Float_t maxfacpthard){ fMaxFacPtHard = maxfacpthard;} //FK


  void SetEmbedding(TString inhybrjets){  //EMB
      fDoEmbedding         = kTRUE;                            //EMB
      fHybridJetContName   = inhybrjets;                       //EMB
  }


protected:
  // Implementation of interface methods
  virtual void    UserCreateOutputObjects();

  virtual void    ExecOnce();

  virtual void    UserExec(Option_t* option);

  virtual Bool_t  UserNotify();

  virtual void    Terminate(Option_t* option);


  // copy constructo not implemented yet
  AliAnalysisTaskEmcalJetBtagSV           (const AliAnalysisTaskEmcalJetBtagSV&);

  // assignment operator not implemented yet
  AliAnalysisTaskEmcalJetBtagSV& operator=(const AliAnalysisTaskEmcalJetBtagSV&);

  void            AnalyseCorrectionsMode();

  void            AnalyseDataMode();

  void GetFlavour2Methods(AliEmcalJet* jet,
                          Double_t (&partonnat)[2],
                          Double_t (&ptpart)[2],
                          Double_t radius);
  // temporary
  Bool_t    GetArrays();
  void      CheckTrackQAinJets();

  Bool_t    FillMapOfV0gTrkIDs();
  Bool_t    FillVecOfV0gTrkIDs(std::vector<Int_t>& vctrTrkIDs);
  Bool_t    IsAODtrkBelongToV0(std::vector<Int_t>& vctrTrkIDs, Int_t trkID);

  Double_t  GetExternalRho(Bool_t isMC = kFALSE);
  Double_t  GetDeltaPtRandomCone(Double_t jetradius, Double_t rhovalue);
  Int_t  FillDeltaPt(Double_t , Int_t , AliAODVertex* , vctr_pair_dbl_int ,  Double_t , Double_t ,  Int_t ); 	
  //Int_t  FillDeltaPt(Int_t, AliAODVertex* , vctr_pair_dbl_int, Double_t, Int_t);      //newDeltaPt//
  Double_t GetDeltaPtRandomConeWithoutSignalPt (Double_t, Double_t, Double_t, Double_t);	 
  Bool_t IsOutlier(); //FK// Tests if the event is pthard bin outlier 



private:


 
  Bool_t      fCorrMode;             // enable correction or data modes
  Bool_t      fDoBkgRej;             // enable background rejection
  Bool_t      fDoRndmCone;           // enable random cone method
  Bool_t      fDoQAVtx;              // enable output of qa on secondary vertex
  Bool_t      fDoFillV0Trks;         // enable V0 checks
  Bool_t      fDoDetRespMtx;        // enable detector repsonse matrix output
  Bool_t      fDoOnlyMtxAna;        // enable only det. matrix analysis
  Bool_t      fUseTriggerData;       // use emacal trigger
  Bool_t      fEmbeddPerpendicular;  // EMB_clus use perpendicular track embedding


  const char* fRecJetsBranch;        // name of the AOD REC-jets branch
  const char* fGenJetsBranch;        // name of the AOD GEN-jets branch

  TString     fPtHardName;
  TString     fGenNamePattern;

  TString     fJetContName;          //  Name of the found jet array
  TString     fTrkContName;          //  Name of the found track array
  TString     fRhoTaskName;          //  Name of the rho task
  TString     fMCJetContName;        //  Name of the found jet array
  TString     fMCTrkContName;        //  Name of the found track array
  TString     fMCRhoTaskName;        //  Name of the rho task

  Double_t    fTaggingRadius;        // radius used in tagging the jet flavour
  Double_t    fSigmaSVCut;           // cut for delPt histogram                    //newDeltaPt//
  //
  // MC weights
  //
  Double_t    fMCWeight;             //!<!  pT-hard bin MC weight. It is used only internally.
  Double_t    fMCXsec;               //!<!
  Double_t    fMCAvgTrials;          //!<!

  Float_t     fZNApercentile;        //!<! multiplicity percentile, ZNA estimator

  TString     fCurrFileName;         //!<!  Current file path name.

  Bool_t      fCheckMCCrossSection;  ///<  Retrieve from the pyxsec.root file the cross section, only if requested.
  Bool_t      fSkipWeightInfo;       //!<!
  Bool_t      fUseWeight;
  Bool_t      fInitialized;          //!<! 

  TList*                      fOutputList;       //!<! list of output objects

  // AliHFJetsContainerVertex
  AliHFJetsContainerVertex*   fhJetVtxSim;       //!<! properties of vertices within the jet MC
  AliHFJetsContainerVertex*   fhJetVtxData;      //!<! properties of vertices within the jet Data
  AliHFJetsContainerVertex*   fhQaVtx;           //!<! vertices properties

  TH1F*                       fhEntries;         //!<!
  TH1F*	                      fhZNApercentQa;    //!<! QA histo for ZNA percentile
  TH1F*                       fhEvtRej;          //!<! Rejection criteria.
  TH1F*                       fhEvtRejBitmap;          //!<! Rejection criteria bitmap.	
  TH1F*                       fhHFjetQa;         //!<! Various QA check on Jet.
  TH1F*                       fhRhoQa;           //!<! 
  TH1F*                       fhMCRhoQa;         //!<! 
  TH1F*                       fhDeltaPt;         //!<! delta pt distribution
  TH1F*                       fhDeltaPtLxy5;     //!<! delta pt distribution, Lxy=5 //newDeltaPt//
  TH1F*                       fhDeltaPtLxy6;     //!<! delta pt distribution, Lxy=6 //newDeltaPt//
  TH1F*                       fhDeltaPtLxy7;     //!<! delta pt distribution, Lxy=7 //newDeltaPt//	
  TH1F*                       fhDeltaPtTrack10;  //!<! delta pt distribution, pt track > 10 //newDeltaPt//	
	
  TH1F*                       fZVertex;          //!<! Z vertex distribuition //AID//	
  TH2F*                       fhTrackEta;        //!<! eta inclusive track distribuition //AID//	
  TH2F*                       fhTrackPhi;        //!<! phi inclusive track distribuition //AID//	
  TH2F*                       fhJetEta;          //!<! eta inclusive jet distribuition //AID//	
  TH2F*                       fhJetPhi;          //!<! phi inclusice jet distribuition //AID//		
  TH2F*                       fhOneOverPtVsPhiNeg;//!<! 1/p_T,track  versus phi for negative tracks //AID//		
  TH2F*                       fhOneOverPtVsPhiPos;//!<! 1/p_T,track  versus phi for positive tracks //AID//		
  TH2F*                       fhSigmaPtOverPtVsPt;//!<! resolution of 1/p_T,track  versus p_T,track //AID//		
  TH2F*                       fhDCAinXVsPt;      //!<! X DCA versus pT  //AID// 
  TH2F*                       fhDCAinYVsPt;      //!<! Y DCA versus pT  //AID//
  TH2F*                       fhDCAinXVsPtPhysPrimary; //!<! X DCA versus pT for physical primaries  //AID// 
  TH2F*                       fhDCAinYVsPtPhysPrimary; //!<! Y DCA versus pT for physical primaries  //AID//
  TH2F*                       fhDCAinXVsPtSecondary; //!<! X DCA versus pT for secondaries //AID// 
  TH2F*                       fhDCAinYVsPtSecondary; //!<! Y DCA versus pT for secondaries //AID//
  TH2D*                       fhFractionOfSecInJet; //!<! Fraction of jet pT carried by secondaries //AID//
  TH1D*                       fhPtTrkTruePrimRec; //!<! pt spectrum of true reconstructed primary tracks    
  TH1D*                       fhPtTrkTruePrimGen; //!<! pt spectrum of true generated primary track    
  TH1D*                       fhPtTrkSecOrFakeRec; //!<! pt spectrum of reconstructed fake or secondary tracks     
	
  THnSparse*                  fhnDetRespMtx;     //!<! THnSparse to fill response matrix
  THn*                        fhnGenerated;      //!<! THn to fill MC generated histo

  TProfile*                   fhXsec;            //!<! Cross section in PYTHIA.
  TH1F*                       fhTrials;          //!<! Number of event trials in PYTHIA.

  AliVEvent*                  fEvent;            //! Input event
  AliAODMCHeader*             fMCHeader;         //! Input MC header
  AliMultSelection*           fMultSelection;    //! multiplicity/centrality selector

  AliHFJetsTaggingVertex*     fTagger;           // Jet Tagging object

  AliRDHFJetsCuts*            fCutsHFjets;       //  specific algo jet cut object
  AliAnalysisUtils*           fAnalysisUtils;    //  points to class with common analysis utilities
  AliParticleContainer*       fMCTracksCont;     //! MC tracks

  TClonesArray*               fRecJetArray;      //! Array of the found jets
  TClonesArray*               fRecTrkArray;      //! Array of PicoTracks
  TClonesArray*               fMCJetArray;       //! Array of the found mc jets
  TClonesArray*               fMCPartArray;      //! Array of MC particles for given event

  TClonesArray*               fHFvertexing;      //! Array of reconstructed secondary vertex (b-tagged jets)

  map_int_bool*               fV0gTrkMap;        //!

  TRandom3*                   fRandom;     	 //! used for throwing random cone

  Int_t                       fGlLogLevel;       ///<
  Int_t                       fLcDebLevel;       ///<
  Int_t                       fStartBin;         //!
  Float_t                     fMaxFacPtHard;     //<! Cut on  pthat events. How many times can be jet pT larger than pthat //FK
  Double_t                    fPtCut;            //<! min cut on track pT   //AID  
  Double_t                    fEtaCut;           //<! cut on track eta   //AID  

  Bool_t                      fDoEmbedding;         ///< EMB flag to do embedding from file
  TString                     fHybridJetContName;   ///< EMB Name of the hybrid jet container created from tracks and embedded tracks 
  TClonesArray*               fHybridJetCont;       ///< EMB hybrid jet container
  TH1F*                       fhDeltaPtEmbedd;      //!<! EMB delta pt distribution, based on embedding of tracks to event//  //AID_emb
  TH2F*                       fhDeltaPtEmbeddCorrelation;  //!<! EMB


  TH1F*                       fhDeltaPtEmbeddPerpendicular;      //!<! EMB_clus 
  TH2F*                       fhDeltaPtEmbeddCorrelationPerpendicular;  //!<! EMB
  AliFJWrapper*		      fFastJetWrapper;  	///< EMB_clus wrapper for fast jet finding
  TRandom*		      fTrackGenerator; 		///< EMB_clus generator for track perpendicular to signal jet


  ClassDef(AliAnalysisTaskEmcalJetBtagSV, 12);  // analysis task for MC study //AID//
};

//-------------------------------------------------------------------------------------
inline Bool_t AliAnalysisTaskEmcalJetBtagSV::UserNotify()
{
  if (fCheckMCCrossSection) {
    GetPythiaCrossSection();
    AliDebugF(3, MSGDEBUG("MC pT-hard weight: %e"), fMCWeight);
  }

  fSkipWeightInfo = kFALSE;
  return kTRUE;
}

#endif
