#ifndef ALIANLYSISTASKGAMMACONVFLOW_cxx
#define ALIANLYSISTASKGAMMACONVFLOW_cxx

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliV0ReaderV1.h"
#include "AliConvEventCuts.h"
#include "AliKFConversionPhoton.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliConversionMesonCuts.h"
#include "AliAnalysisManager.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"


class AliFlowTrackCuts;
class AliFlowTrack;
class AliFlowEvent;
class AliFlowCandidateTrack;
class AliFlowEventSimple;


class AliAnalysisTaskGammaConvFlow : public AliAnalysisTaskSE {
    
public:
    AliAnalysisTaskGammaConvFlow();
    AliAnalysisTaskGammaConvFlow(const char *name);
    virtual ~AliAnalysisTaskGammaConvFlow();
    
    virtual void   UserCreateOutputObjects();
    virtual Bool_t Notify();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);
//     void InitBack();
    
    void SetIsHeavyIon(Int_t flag){
        fIsHeavyIon = flag;
        
    }
    
    void SetDoMesonAnalysis(Bool_t flag){fDoMesonAnalysis = flag;}
    void SetDoMesonQA(Int_t flag){fDoMesonQA = flag;}
    void SetDoPhotonQA(Int_t flag){fDoPhotonQA = flag;}
    void ProcessPhotonCandidates();
  //  void CalculatePi0Candidates();
    //void CalculateBackground();
  //  void CalculateBackgroundRP();

    void RelabelAODPhotonCandidates(Bool_t mode);

    //void RotateParticle(AliAODConversionPhoton *gamma);
    //void RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP);
    
    void SetEventCutList(Int_t nCuts, TList *CutArray){
        fnCuts = nCuts;
        fEventCutArray = CutArray;
    }
    void SetConversionCutList(Int_t nCuts, TList *CutArray){
        fnCuts = nCuts;
        fCutArray = CutArray;
    }
 //   void SetMesonCutList(Int_t nCuts, TList *CutArray){
 //       fnCuts = nCuts;
 //       fMesonCutArray = CutArray;
 //   }
    
    // BG HandlerSettings
   // void SetMoveParticleAccordingToVertex(Bool_t flag){fMoveParticleAccordingToVertex = flag;}
    void FillPhotonCombinatorialBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode[]);
   // void MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
    void UpdateEventByEventData();
    void SetLogBinningXTH2(TH2* histoRebin);
  //  Int_t GetSourceClassification(Int_t daughter, Int_t pdgCode);
    
    
    //=========================================Flow Staff=========================================================================================
    template <typename T> void           SetNullCuts(T* esd);
    void                                 PrepareFlowEvent(Int_t iMulti, AliFlowEvent *FlowEv) const;
    void                                 SetRPCuts(AliFlowTrackCuts *cutsRP) { fCutsRP = cutsRP; }
    void                                 ProcessPhotonCandidatesforV2();
    void                                 MakeTrack(Double_t mass, Double_t pt, Double_t phi, Double_t eta/*, Int_t nDau, Int_t iID[]*/) const;
    //=========================================ENd Flow Staff=========================================================================================
    
    
protected:
    AliV0ReaderV1 						*fV0Reader;											//
    AliGammaConversionAODBGHandler 		**fBGHandler;										//
    AliConversionAODBGHandlerRP    		**fBGHandlerRP;										//
    AliVEvent 							*fInputEvent;										//
//    AliMCEvent 							*fMCEvent;											//
//    AliStack 							*fMCStack;											//
    TList 								**fCutFolder;										//
    TList 								**fESDList;											//
    TList 								**fBackList;										//
    TList 								**fMotherList;										//
    TList 								**fPhotonDCAList;									//
    TList 								**fMesonDCAList;		 							//
 //   TList 								**fTrueList;										//
 //   TList 								**fMCList;											//
    TList 								**fHeaderNameList;									//
    TList 								*fOutputContainer;									//
    TClonesArray 						*fReaderGammas;										//
    TList 								*fGammaCandidates;									//
    TList 								*fEventCutArray;									//
    AliConvEventCuts 					*fEventCuts;										//
    TList 								*fCutArray;											//
    AliConversionPhotonCuts 			*fConversionCuts;									//
    TList 								*fMesonCutArray;									//
    AliConversionMesonCuts 				*fMesonCuts;										//
    TH1F 								**hESDConvGammaPt;									//
    TH1F 								**hESDConvGammaR;									//
    TH1F 								**hESDConvGammaEta;									//
  //  TTree 								**tESDConvGammaPtDcazCat;							//
    Float_t 							fPtGamma;											//
    Float_t 							fDCAzPhoton;										//
    Float_t 							fRConvPhoton;										//
    Float_t 							fEtaPhoton;											//
    UChar_t 							iCatPhoton;											//
    UChar_t 							iPhotonMCInfo; 										//
    // 0: garbage,
    // 1: background
    // 2: secondary photon not from eta or k0s,
    // 3: secondary photon from eta,
    // 4: secondary photon from k0s,
    // 5: dalitz
    // 6: primary gamma
    TH2F 								**hESDMotherInvMassPt;								//
  //  THnSparseF 							**sESDMotherInvMassPtZM;							//
    TH2F 								**hESDMotherBackInvMassPt;							//
  //  THnSparseF 							**sESDMotherBackInvMassPtZM;						//
    TH2F 								**hESDMotherInvMassEalpha;							//
    TH2F 								**hESDMotherPi0PtY;									//
    TH2F 								**hESDMotherEtaPtY;									//
    TH2F 								**hESDMotherPi0PtAlpha;								//
    TH2F 								**hESDMotherEtaPtAlpha;								//
    TH2F 								**hESDMotherPi0PtOpenAngle;							//
    TH2F 								**hESDMotherEtaPtOpenAngle;							//

    TH1I 								**hNEvents;											//
    TH1I 								**hNGoodESDTracks;									//
    TH1I 								**hNGammaCandidates;								//
    TH2F 								**hNGoodESDTracksVsNGammaCanditates;				//
    TH1I 								**hNV0Tracks;										//
    TProfile 							**hEtaShift;										//
    TTree 								**tESDMesonsInvMassPtDcazMinDcazMaxFlag;			//
    Float_t 							fInvMass;											//
    Float_t 							fPt;												//
    Float_t 							fDCAzGammaMin;										//
    Float_t 							fDCAzGammaMax;										//
    UChar_t 							iFlag;												//
    UChar_t 							iMesonMCInfo; 										//
    // 0: garbage,
    // 1: background
    // 2: secondary meson not from eta or k0s,
    // 3: secondary meson from eta,
    // 4: secondary meson from k0s,
    // 5: dalitz
    // 6: primary meson gamma-gamma-channel
    Double_t 							fEventPlaneAngle; 									// EventPlaneAngle
    TRandom3 							fRandom;											//
    Int_t 								fnGammaCandidates;									//
    Double_t 							*fUnsmearedPx;										//[fnGammaCandidates]
    Double_t 							*fUnsmearedPy;										//[fnGammaCandidates]
    Double_t 							*fUnsmearedPz;										//[fnGammaCandidates]
    Double_t 							*fUnsmearedE;										//[fnGammaCandidates]
    Int_t 								*fMCStackPos;										//[fnGammaCandidates]
    Int_t								*fMCStackNeg;										//[fnGammaCandidates]
    Int_t 								*fESDArrayPos;										//[fnGammaCandidates]
    Int_t 								*fESDArrayNeg;										//[fnGammaCandidates]
    Int_t 								fnCuts;												//
    Int_t 								fiCut;												//
    Bool_t 								fMoveParticleAccordingToVertex;						//
    Int_t 								fIsHeavyIon;										//
    Bool_t 								fDoMesonAnalysis;									//
    Int_t 								fDoMesonQA;											//
    Int_t 								fDoPhotonQA;										//
    Bool_t 								fIsFromMBHeader;									//
    
private:
    TObjArray            *fCandidates; // candidate array
    Int_t                fDebug; // debug level (0 none, 1 fcn calls, 2 verbose)
    AliFlowTrackCuts     *fCutsRP; // track cuts for reference particles
    AliFlowTrackCuts     *fNullCuts; // dummy cuts for flow event tracks
    AliFlowEvent         *fFlowEvent; //! flow events

    
    
    
    
    AliAnalysisTaskGammaConvFlow(const AliAnalysisTaskGammaConvFlow&); // Prevent copy-construction
    AliAnalysisTaskGammaConvFlow &operator=(const AliAnalysisTaskGammaConvFlow&); // Prevent assignment
    ClassDef(AliAnalysisTaskGammaConvFlow, 2);
};

#endif
