#ifndef ALIANALYSISTASKSEF01710FROMAODTRACKS_H
#define ALIANALYSISTASKSEF01710FROMAODTRACKS_H

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
#include "TLorentzVector.h"
#include "TVector.h"

#include "AliAnalysisTaskSE.h"
#include "AliPIDResponse.h"
#include "AliAODEvent.h"
#include "AliPID.h"

/// \class AliAnalysisTaskSEF01710fromAODtracks

class THnSparse;
class TH1F;
class TH2F;
class TClonesArray;
class AliESDtrackCuts;
class AliESDVertex;
class AliAODMCParticle;

class AliAnalysisTaskSEF01710fromAODtracks : public AliAnalysisTaskSE 
{
  public:
    AliAnalysisTaskSEF01710fromAODtracks();
    AliAnalysisTaskSEF01710fromAODtracks(const Char_t* name,  Bool_t writeVariableTree=kTRUE);
    virtual ~AliAnalysisTaskSEF01710fromAODtracks();

    // Implementation of interface methods  
    virtual void UserCreateOutputObjects();
    virtual void Init();
    virtual void LocalInit() {Init();}
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);

    void FillROOTObjects(AliAODEvent *aod,AliAODv0 *v01, AliAODv0 *v02,TClonesArray *mcArray);
    void FillROOTObjects_ChargedKaon(AliAODEvent *aod,AliAODTrack *trk1, AliAODTrack *trk2,TClonesArray *mcArray);
    void FillMixROOTObjects(TLorentzVector *v1, TLorentzVector *v2, TVector *v1v, TVector *v2v);
    void MakeAnalysis(AliAODEvent *aod, TClonesArray *mcArray);
    void SetMC(Bool_t ismc){fUseMCInfo=ismc;}
    void SetAnalysisType(Int_t at){fAnalysisType=at;}
    void MakeMCAnalysis(TClonesArray *mcArray);
    Int_t MatchToMC(AliAODv0 *v0a, AliAODv0 *v0b, TClonesArray *mcArray, Int_t *pdgv0a_array, Int_t *pdgv0b_array, Int_t *labelv0a_array, Int_t *labelv0b_array, Int_t &ngen_v0a, Int_t &ngen_v0b);
    Int_t MatchToMC_ChargedKaon(AliAODTrack *ta, AliAODTrack *tb, TClonesArray *mcArray, Int_t *pdgta_array, Int_t *pdgtb_array, Int_t *labelta_array, Int_t *labeltb_array, Int_t &ngen_ta, Int_t &ngen_tb);

    void SelectTrack( const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks,Int_t *seleFlags, TClonesArray *mcArray);
    void SelectV0( const AliVEvent *event,Int_t nV0s,Int_t &nSelev0, Bool_t *seleV0Flags, TClonesArray *mcArray);
    Bool_t SingleV0Cuts(AliAODv0 *v0, AliAODVertex *primVert);
    Bool_t IsEventSelected(AliVEvent *event);
    Bool_t IsEventAcceptedPhysicsSelection(AliVEvent *event);
    Bool_t IsEventAcceptedTrigger(AliVEvent *event);

    /// mixing
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
    void DoEventMixingWithPools(Int_t index);
    void ResetPool(Int_t poolIndex);
    Int_t GetPoolIndex(Double_t zvert, Double_t mult);

  private:

    AliAnalysisTaskSEF01710fromAODtracks(const AliAnalysisTaskSEF01710fromAODtracks &source);
    AliAnalysisTaskSEF01710fromAODtracks& operator=(const AliAnalysisTaskSEF01710fromAODtracks& source); 

    void DefineTreeVariables();
    void DefineGeneralHistograms();
    void DefineAnalysisHistograms();

    Int_t fAnalysisType; /// Analysis Type 0: pp(lhc10bcde), Type 1(),,,,
    Bool_t fUseMCInfo;          /// Use MC info
    TList *fOutput;             //!<! User output slot 1 // general histos
    TList *fOutputAll;          //!<! User output slot 2 // Analysis histos
    TH1F *fCEvents;             //!<! structure for event mixing
    TH1F *fCEventsNorm;             //!<! structure for event mixing
    TH1F *fHTrigger;            //!<! Histograms to check trigger
    TH1F *fHCentrality;         //!<! histogram to check centrality
    TH1F *fHVtxZ;         //!<! histogram to check z-vertex distribution
    TH1F *fHMultiplicityV0A;         //!<! histogram to check vzero a multiplicity
    TH1F *fHMultiplicityV0C;         //!<! histogram to check vzero c multiplicity
    Bool_t fIsEventSelected;           /// flag for event selected
    Int_t fWhyRejection;           /// Why the event is rejected
    Bool_t    fWriteVariableTree;       /// flag to decide whether to write the candidate variables on a tree variables
    TTree    *fVariablesTree;           //!<! tree of the candidate variables after track selection on output slot 4
    Bool_t fIsMB;            /// Is MB event
    Bool_t fIsSemi;          /// is semi-central trigger event
    Bool_t fIsCent;          /// is central trigger event
    Bool_t fIsINT7;          ///  is int7 trigger event
    Bool_t fIsEMC7;          /// is emc7 trigger event
    Float_t *fCandidateVariables;     //!<! variables to be written to the tree
    AliAODVertex *fVtx1;              /// primary vertex
    AliESDVertex *fV1;                /// primary vertex
    Float_t fVtxZ;              /// primary vertex
    Double_t fBzkG;                   /// magnetic field value [kG]
    Float_t  fCentrality;             ///centrality
    Float_t  fTriggerCheck;           /// Trigger information
    Float_t  fMultiplicityVZEROA;             //!<! Multiplicity in vzero A
    Float_t  fMultiplicityVZEROC;             //!<! Multiplicity in vzero C
    AliPIDResponse *fPIDResponse; /// pid response

    //--------------------- My histograms ------------------
    THnSparse*  fHistoF01710Mass;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF01710MassMix;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF01710MassMCGen;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF01710MassMCS;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF21525MassMCGen;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF21525MassMCS;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoA21320MassMCGen;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoA21320MassMCS;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF21270MassMCGen;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF21270MassMCS;        //!<! F0 1710 mass spectra

    THnSparse*  fHistoF01710ChargedMass;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF01710ChargedLikeMass;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF01710ChargedMassMix;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF01710ChargedLikeMassMix;        //! F0 1710 mass spectra
    THnSparse*  fHistoF01710ChargedMassMCGen;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF01710ChargedMassMCS;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF21525ChargedMassMCGen;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF21525ChargedMassMCS;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoA21320ChargedMassMCGen;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoA21320ChargedMassMCS;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF21270ChargedMassMCGen;        //!<! F0 1710 mass spectra
    THnSparse*  fHistoF21270ChargedMassMCS;        //!<! F0 1710 mass spectra

    TH1D *fHistonEvtvsRunNumber;//!<! nevt vs runnumber
    THnSparse *fHistonKaonvsRunNumber;//!<! nkaon vs runnumber
    THnSparse *fHistonK0vsRunNumber;//!<! nk0 vs runnumber

    //Mixing
    Int_t fDoEventMixing; /// flag for event mixing
    Int_t  fNumberOfEventsForMixing; /// maximum number of events to be used in event mixing
    Int_t fNzVtxBins;								/// number of z vrtx bins
    Double_t fZvtxBins[100];						// [fNzVtxBinsDim]
    Int_t fNCentBins;								/// number of centrality bins
    Double_t fCentBins[100];						// [fNCentBinsDim]
    Int_t  fNOfPools; /// number of pools
    TTree** fEventBuffer;   //!<! structure for event mixing
    TObjString *fEventInfo; ///unique event id for mixed event check
    TObjArray* fK0Tracks; /// array of electron-compatible tracks
    TObjArray* fK0CutVarsArray; /// array of cut variables
    TObjArray* fKaonTracks; /// array of electron-compatible tracks
    TObjArray* fKaonCutVarsArray; /// array of cut variables

    /// \cond CLASSIMP 
    ClassDef(AliAnalysisTaskSEF01710fromAODtracks,1); // class for f0(1710)->KK
    /// \endcond 
};
#endif

