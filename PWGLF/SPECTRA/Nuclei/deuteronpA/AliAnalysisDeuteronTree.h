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

    //Variables for Tree
    Float_t fCentrality; //
    Float_t fPt; //
    Float_t fMom; //
    Float_t fRapd; //
    Float_t fNsigmaTPCd; //
    Float_t fNsigmaTOFd; //
    Float_t fDcaXYd; //
    Float_t fMcCode; //
    
    //Variables for empirical momentum correction
    Float_t fMomCorrConstA;
    Float_t fMomCorrConstB;
    Float_t fMomCorrPower;
    
    AliAnalysisDeuteronTree(const AliAnalysisDeuteronTree&);            // not implemented
    AliAnalysisDeuteronTree& operator=(const AliAnalysisDeuteronTree&); // not implemented
    
    ClassDef(AliAnalysisDeuteronTree, 1);
};
#endif
