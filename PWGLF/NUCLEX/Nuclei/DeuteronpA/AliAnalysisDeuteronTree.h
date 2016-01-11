#ifndef ALIANALYSISDEUTERONTREE_H
#define ALIANALYSISDEUTERONTREE_H

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// This analysis extracts information into TTree for pT-spectra of deuterons //
// Based on AliAnalysisDeuteronpA task of J. Anielski for deuteron analysis  //
// and AliAnalysisTaskExtractV0 by D. Chinellato for TTree interface         //
// L.Barnby October 2015                                                     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

class TH1;
class TH1F;
class TH2;
class TH2F;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliHeader;
class AliPIDResponse;
class AliAnalysisUtils;

#include "AliAnalysisTaskSE.h"

class AliAnalysisDeuteronTree : public AliAnalysisTaskSE {
  public:
    AliAnalysisDeuteronTree(const char *name);
    AliAnalysisDeuteronTree();
    virtual ~AliAnalysisDeuteronTree();
    //
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    //
    void           SetESDtrackCuts(AliESDtrackCuts * trackCuts){fESDtrackCuts = trackCuts;};
    //void           SetAlephParameters(const Double_t * parameters){for(Int_t j=0;j<5;j++) fAlephParameters[j] = parameters[j]; Initialize();};
    void           SetIsMCtrue(Bool_t isMCdata = kTRUE){fMCtrue = isMCdata;};
    void           SetRapCMSpA(Bool_t isRapCMSpA = kTRUE){fRapCMSpA = isRapCMSpA;};
    void           SetMinPtCut(Float_t minPt = 0.5){fMinPtCut = minPt;};
    void           SetMinPtTOFCut(Float_t minPt = 1.0){fMinPtTOFCut = minPt;};
    void           SetTOFNsigmaCutdMin(Float_t nSigma = -5.0){fNsigmaTOFdMin = nSigma;};
    void           SetTOFNsigmaCutdMax(Float_t nSigma =  5.0){fNsigmaTOFdMax = nSigma;};
    void           SetTPCNsigmaCutdMin(Float_t nSigma = -5.0){fNsigmaTPCdMin = nSigma;};
    void           SetTPCNsigmaCutdMax(Float_t nSigma =  5.0){fNsigmaTPCdMax = nSigma;};
    void           SetTPCNsigmaCutpiAbs(Float_t nSigma = 3.0){fNsigmaTPCpiAbsCut = nSigma;};
    //
    void           Initialize();
    //

private:


    //AliESDEvent *fESD;                   //! ESD object
    AliESDtrackCuts * fESDtrackCuts;     // basic cut variables
    //AliESDtrackCuts * fESDTrackCutsMult; // cuts for the MULTIPLICITY DETERMINATION
    AliPIDResponse       * fPIDResponse;      // official PID response
    AliAnalysisUtils  *fUtils;           // For vertex cut and pileup rejection

    // For controlling Task behaviour
    Bool_t        fMCtrue;               // flag if real data or MC is processed
    Bool_t        fRapCMSpA;             // flag if shift to CMS_NN system for pA
  
    TTree	*fTree;                     //! Output Tree
    TList   *fListHist;                 //! list for histograms

    // Few histograms for monitoring
    TH1F* fhZVertex; //! event Z vertex distribution
    TH1F* fhCentrality; //! centrality distribution
    TH2F* fhNsigmaTPCvsMom; //! track N-sigma for TPC dE/dx
    TH2F* fhNsigmaTOFvsMom; //! track N-sigma for TOF
    
    //Variables for Tree
    Float_t fTimeStamp; //
    Float_t fCentrality; //
    Float_t fPtCor; //
    Float_t fPt; //
    Float_t fMom; // Momentum multiplied by charge (saves using one more variable)
    Float_t fRapd; //
    Float_t fPxd; //
    Float_t fPyd; //
    Float_t fPzd; //
    Float_t fNsigmaTPCd; //
    Float_t fNsigmaTOFd; //
    Float_t fDcaXYd; //
    Float_t fMcCode; //
    Int_t fNpion; //
    // Pion tree variables
//    const Int_t maxPions = 400; // array size with const variable didn't work
    Float_t fPx[400];
    Float_t fPz[400];
    Float_t fPy[400];
    Float_t fNsigmaTPCpi[400]; //
    Int_t fCharge[400]; //!

    
    //Variables for empirical momentum correction
    Float_t fMomCorrConstA;
    Float_t fMomCorrConstB;
    Float_t fMomCorrPower;
    
    //Other, including cuts
    Float_t fMinPtCut;
    Float_t fMinPtTOFCut;
    Float_t fNsigmaTPCdMin;
    Float_t fNsigmaTPCdMax;

    Float_t fNsigmaTOFdMin;
    Float_t fNsigmaTOFdMax;
    Float_t fNsigmaTPCpiAbsCut;
    
    AliAnalysisDeuteronTree(const AliAnalysisDeuteronTree&);            // not implemented
    AliAnalysisDeuteronTree& operator=(const AliAnalysisDeuteronTree&); // not implemented
    
    ClassDef(AliAnalysisDeuteronTree, 1);
};
#endif
