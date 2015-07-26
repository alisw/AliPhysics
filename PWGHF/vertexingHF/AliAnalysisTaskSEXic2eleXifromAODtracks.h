#ifndef ALIANALYSISTASKSEXIC2ELEXIFROMAODTRACKS_H
#define ALIANALYSISTASKSEXIC2ELEXIFROMAODTRACKS_H

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
#include "AliRDHFCutsXictoeleXifromAODtracks.h"

class THnSparse;
class TH1F;
class TClonesArray;
class AliAODRecoCascadeHF;
class AliESDVertex;
class AliAODMCParticle;
class AliEventPoolManager;
class AliNormalizationCounter;

class AliAnalysisTaskSEXic2eleXifromAODtracks : public AliAnalysisTaskSE 
{
 public:
  AliAnalysisTaskSEXic2eleXifromAODtracks();
  AliAnalysisTaskSEXic2eleXifromAODtracks(const Char_t* name, AliRDHFCutsXictoeleXifromAODtracks* cuts, Bool_t writeVariableTree=kTRUE);
  virtual ~AliAnalysisTaskSEXic2eleXifromAODtracks();

  // Implementation of interface methods  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void FillROOTObjects(AliAODRecoCascadeHF *elobj, AliAODcascade *casc, AliAODTrack *trk, TClonesArray *mcArray, Bool_t mixing);
  void FillElectronROOTObjects(AliAODTrack *trk, TClonesArray *mcArray);
  void FillCascROOTObjects(AliAODcascade *casc, TClonesArray *mcArray);
  void FillMCROOTObjects(AliAODMCParticle *part, AliAODMCParticle *mcepart, AliAODMCParticle *mcv0part, Int_t decaytype);
  void MakeMCAnalysis(TClonesArray *mcArray);
  void MakeAnalysis(AliAODEvent *aod, TClonesArray *mcArray);

  void SelectCascade( const AliVEvent *event,Int_t nCasc,Int_t &nSeleCasc, Bool_t *seleCascFlags, TClonesArray *mcArray);
  void SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Bool_t *seleFlags, TClonesArray *mcArray);

  // set MC usage
  void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t GetMC() const {return fUseMCInfo;}
  void SetUseCentralityV0M(Bool_t centon) {fUseCentralityV0M = centon;}
  Bool_t GetUseCentralityV0M() const {return fUseCentralityV0M;}
  void SetWriteEachVariableTree(Bool_t a) {fWriteEachVariableTree = a;}
  Bool_t GetWriteEachVariableTree() const {return fWriteEachVariableTree;}
  void SetWriteMCVariableTree(Bool_t a) {fWriteMCVariableTree = a;}
  Bool_t GetWriteMCVariableTree() const {return fWriteMCVariableTree;}

  void SetReconstructPrimVert(Bool_t a) { fReconstructPrimVert=a; }

  AliAODRecoCascadeHF* MakeCascadeHF(AliAODcascade *casc, AliAODTrack *trk, AliAODEvent *aod, AliAODVertex *vert, Bool_t mixing);
  AliAODVertex* ReconstructSecondaryVertex(AliAODcascade *casc, AliAODTrack *trk, AliAODEvent *aod);
	Int_t MatchToMC(AliAODRecoCascadeHF *elobj, TClonesArray *mcArray, Int_t *pdgarray_ele, Int_t *pdgarray_casc, Int_t *labelarray_ele, Int_t *labelarray_casc,  Int_t &ngen_ele, Int_t &ngen_casc);
	Int_t MatchToMCCascade(AliAODcascade *theCascade, Int_t pdgabscasc, Int_t *pdgDgcasc, Int_t *pdgDgv0, TClonesArray *mcArray) const;


  //mixing
  void SetEventMixingWithPools(){fDoEventMixing=1;}
  void SetEventMixingOff(){fDoEventMixing=0;}
	void SetNumberOfEventsForMixing(Int_t events){fNumberOfEventsForMixing=events;}
	void SetPoolPVzBinLimits(Int_t Nzvtxbins,const Double_t *ZvtxBins){
		fNzVtxBins = Nzvtxbins;
		for(int ix = 0;ix<fNzVtxBins+1;ix++){fZvtxBins[ix] = ZvtxBins[ix];}
	}
	void SetPoolCentBinLimits(Int_t Ncentbins,const Double_t *CentBins){
		fNCentBins = Ncentbins;
		for(int ix = 0;ix<fNCentBins+1;ix++){fCentBins[ix] = CentBins[ix];}
	}
  void DoEventMixingWithPools(Int_t index,AliAODEvent *aodEvent, Bool_t *seleFlags);
  void ResetPool(Int_t poolIndex);
  Int_t GetPoolIndex(Double_t zvert, Double_t mult);


 private:

  AliAnalysisTaskSEXic2eleXifromAODtracks(const AliAnalysisTaskSEXic2eleXifromAODtracks &source);
  AliAnalysisTaskSEXic2eleXifromAODtracks& operator=(const AliAnalysisTaskSEXic2eleXifromAODtracks& source); 

  void DefineTreeVariables();
  void DefineEleTreeVariables();
  void DefineCascTreeVariables();
  void DefineMCTreeVariables();
  void DefineGeneralHistograms();
  void DefineAnalysisHistograms();

  AliAODVertex *CallPrimaryVertex(AliAODcascade *casc, AliAODTrack *trk, AliAODEvent *evt);
  AliAODVertex* PrimaryVertex(const TObjArray *trkArray,AliVEvent *event);

  Bool_t fUseMCInfo;                 // Use MC info
  TList *fOutput;                    //! User output slot 1 // general histos
  TList *fOutputAll;                 //! User Output slot 3  //analysis histograms 
  TList *fListCuts;                  //! User output slot 2 // Cuts 
  TH1F *fCEvents;                    //! Histogram to check selected events
  TH1F *fHTrigger;                   //! Histogram to check Trigger
  TH1F *fHCentrality;                //! Histogram to check Centrality
  AliRDHFCutsXictoeleXifromAODtracks *fAnalCuts;// Cuts - sent to output slot 2
  Bool_t fIsEventSelected;          // flag for event selected
  Bool_t    fWriteVariableTree;     // flag to decide whether to write the candidate variables on a tree variables
  TTree    *fVariablesTree;         //! tree of the candidate variables after track selection on output slot 4
  Bool_t    fWriteEachVariableTree;     // flag to decide whether to write the candidate variables on a tree variables
  Bool_t    fWriteMCVariableTree;     // flag to decide whether to write the candidate variables on a tree variables
  TTree    *fEleVariablesTree;         //! tree of the candidate variables after track selection on output slot 4
  TTree    *fCascVariablesTree;         //! tree of the candidate variables after track selection on output slot 4
  TTree    *fMCVariablesTree;         //! tree of the candidate variables after track selection on output slot 4
  Bool_t fReconstructPrimVert;       //Reconstruct primary vertex excluding candidate tracks
  Bool_t fIsMB;       //MB trigger event
  Bool_t fIsSemi;     //SemiCentral trigger event
  Bool_t fIsCent;     //Central trigger event
  Bool_t fIsINT7;     //INT7 trigger event
  Bool_t fIsEMC7;     //EMC7 trigger event
  Float_t *fCandidateVariables;   //! variables to be written to the tree
  Float_t *fCandidateEleVariables;   //! variables to be written to the tree
  Float_t *fCandidateCascVariables;   //! variables to be written to the tree
  Float_t *fCandidateMCVariables;   //! variables to be written to the tree
  AliAODVertex *fVtx1;            // primary vertex
  AliESDVertex *fV1;              // primary vertex
  Float_t  fVtxZ;         // zVertex
  Double_t fBzkG;                 // magnetic field value [kG]
  Float_t  fCentrality;           //Centrality
  Float_t  fTriggerCheck;         //Stores trigger information
  Bool_t  fUseCentralityV0M;         //Stores trigger information
  Int_t  fEvNumberCounter;         //EvNumber counter

  //--------------------- My histograms ------------------
  THnSparse* fHistoEleXiMass;         //e-Xi mass spectra
  THnSparse* fHistoEleXiMassRS;         //e-Xi mass spectra (right-sign)
  THnSparse* fHistoEleXiMassWS;         //e-Xi mass spectra (wrong-sign)
  THnSparse* fHistoEleXiMassRSMix;         //e-Xi mass spectra (right-sign)
  THnSparse* fHistoEleXiMassWSMix;         //e-Xi mass spectra (wrong-sign)
  THnSparse* fHistoEleXiMassvsElePtRS;         //e-Xi mass spectra (right-sign)
  THnSparse* fHistoEleXiMassvsElePtWS;         //e-Xi mass spectra (wrong-sign)
  THnSparse* fHistoEleXiMassvsElePtRSMix;         //e-Xi mass-ept spectra (right-sign)
  THnSparse* fHistoEleXiMassvsElePtWSMix;         //e-Xi mass-ept spectra (wrong-sign)
  TH2F* fHistoElePtRS;         //e spectra (right-sign)
  TH2F* fHistoElePtWS;         //e spectra (wrong-sign)
  TH2F* fHistoElePtRSMix;         //e spectra (right-sign, mix)
  TH2F* fHistoElePtWSMix;         //e spectra (wrong-sign, mix)
  THnSparse* fHistoEleXiMassMCS;         //e-Xi mass spectra (Efficiency numerator)
  THnSparse* fHistoEleXiMassMCGen;         //e-Xi mass spectra (Efficiency denominator)
  THnSparse* fHistoEleXiMassvsElePtMCS;         //e-Xi mass-ept spectra (Efficiency numerator)
  THnSparse* fHistoEleXiMassvsElePtMCGen;         //e-Xi mass-ept spectra (Efficiency denominator)
  TH2F* fHistoElePtMCS;         //e spectra (Efficiency numerator)
  TH2F* fHistoElePtMCGen;         //e spectra (Efficiency denominator)

  THnSparse* fHistoElePtvsEtaRS;         //e spectra (right-sign)
  THnSparse* fHistoElePtvsEtaWS;         //e spectra (wrong-sign)
  THnSparse* fHistoElePtvsEtaRSMix;         //e spectra (right-sign, mix)
  THnSparse* fHistoElePtvsEtaWSMix;         //e spectra (wrong-sign, mix)
  THnSparse* fHistoElePtvsEtaMCS;         //e spectra (right-sign) efficiency numerator
  THnSparse* fHistoElePtvsEtaMCGen;         //e spectra (wrong-sign) efficiency denominator

  THnSparse* fHistoElePtvsXiPtRS;         //e-Xi spectra (right-sign)
  THnSparse* fHistoElePtvsXiPtWS;         //e-Xi spectra (wrong-sign)
  THnSparse* fHistoElePtvsXiPtRSMix;         //e-Xi spectra (right-sign, mix)
  THnSparse* fHistoElePtvsXiPtWSMix;         //e-Xi spectra (wrong-sign, mix)
  THnSparse* fHistoElePtvsXiPtMCS;         //e-Xi spectra (right-sign) efficiency numerator
  THnSparse* fHistoElePtvsXiPtMCGen;         //e-Xi spectra (wrong-sign) efficiency denominator

  THnSparse* fHistoElePtvsd0RS;         //e pt-d0 spectra (right-sign)
  THnSparse* fHistoElePtvsd0WS;         //e pt-d0 spectra (wrong-sign)
  THnSparse* fHistoElePtvsd0RSMix;         //e pt-d0 spectra (right-sign, mix)
  THnSparse* fHistoElePtvsd0WSMix;         //e pt-d0 spectra (wrong-sign, mix)
  THnSparse* fHistoElePtvsd0MCS;         //e pt-d0 spectra (right-sign) 
  THnSparse* fHistoElePtvsd0PromptMCS;         //e pt-d0 spectra (right-sign) 
  THnSparse* fHistoElePtvsd0BFeeddownMCS;         //e pt-d0 spectra (right-sign) 

  TH1F* fHistoBachPt;      //! Bachelor pT histogram
  TH1F* fHistoBachPtMCS;      //! Bachelor pT histogram (efficiency numerator)
  TH1F* fHistoBachPtMCGen;      //! Bachelor pT histogram (efficiency denominator)
  TH1F* fHistod0Bach;      //! Bachelor d0 histogram
  TH2F* fHistoXiMassvsPt;     //! Xi mass vs pt histogram
  TH2F* fHistoXiMassvsPtMCS;     //! Xi mass vs pt histogram
  TH2F* fHistoXiMassvsPtMCGen;     //! Xi mass vs pt histogram
  TH2F* fHistoOmegaMassvsPt;     //! Omega mass vs pt histogram
  TH2F* fHistoElectronTPCPID;     //! TPC electron PID
  TH2F* fHistoElectronTOFPID;     //! TOF electron PID
  TH2F* fHistoElectronTPCSelPID;     //! TPC electron PID after selection
  TH2F* fHistoElectronTOFSelPID;     //! TOF electron PID after selection
  TH2F* fHistoElectronTPCPIDSelTOF;     //! TPC electron PID after TOF 3 sigma cut
  TH2F* fHistoElectronTPCPIDSelTOFSmallEta;     //! TPC electron PID after TOF 3 sigma cut (|eta|<0.6)
  TH2F* fHistoElectronTPCPIDSelTOFLargeEta;     //! TPC electron PID after TOF 3 sigma cut (0.8>|eta|>0.6)
  THnSparse* fHistoElePtvsCutVarsRS[23];         //e pt- cut variables (Right-sign)
  THnSparse* fHistoElePtvsCutVarsWS[23];         //e pt- cut variables (Wrong-sign)
  THnSparse* fHistoElePtvsCutVarsMCS[23];         //e pt- cut variables (MCS)

  AliNormalizationCounter *fCounter;//!Counter for normalization
	TH1F *fHistonEvtvsRunNumber;
	TH1F *fHistonElevsRunNumber;
	TH1F *fHistonXivsRunNumber;

  //Mixing
  Int_t fDoEventMixing; // flag for event mixing
  Int_t  fNumberOfEventsForMixing; // maximum number of events to be used in event mixing
	Int_t fNzVtxBins;								// number of z vrtx bins
	Double_t fZvtxBins[100];						// [fNzVtxBinsDim]
	Int_t fNCentBins;								// number of centrality bins
	Double_t fCentBins[100];						// [fNCentBinsDim]
  Int_t  fNOfPools; // number of pools
  TTree** fEventBuffer;   //! structure for event mixing
	TObjString *fEventInfo; //unique event id for mixed event check
  TObjArray* fElectronTracks; // array of electron-compatible tracks


  ClassDef(AliAnalysisTaskSEXic2eleXifromAODtracks,5); // class for Xic->e Xi
};
#endif

