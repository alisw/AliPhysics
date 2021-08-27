#ifndef AliAnalysisTaskKaon2PC_H
#define AliAnalysisTaskKaon2PC_H

/*class TH1F;
class TH2F;
class TString;
class AliAODEvent;
class AliAODTrack;
class AliAODMCParticle;
class TList;
class AliAODVertex;
class AliPIDResponse;
class AliAODv0;*/

//class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliPIDResponse.h"
#include "AliAODv0.h"
#include "AliAODVertex.h"

//#include "THnSparse.h"
//#include "TRandom3.h"

class AliAnalysisTaskKaon2PC : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskKaon2PC();
    AliAnalysisTaskKaon2PC(const char *name);
    virtual ~AliAnalysisTaskKaon2PC();
    
    virtual void     UserCreateOutputObjects();
    virtual void     SetTrackCuts(Double_t c1, Double_t c2, Double_t c3, Double_t c4);
    virtual void     SetV0TrackCuts(Double_t c5, Double_t c6, Double_t c7, Double_t c8, Double_t c9, Double_t c10, Double_t c11, Double_t c12, Double_t c13, Double_t c14, Double_t c15);
    Bool_t AcceptTrack(const AliAODTrack* track);
    Bool_t AcceptV0(const AliAODv0 *v0, Double_t *vertex);
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t * option);
 
 private:

    AliAODEvent     *fAOD;                   //! AOD event
	TList           *fOutput;        // Output list
    AliAODVertex    *fPrimaryVtx;            //! AOD vertex
    AliPIDResponse	*fPIDResponse;	 // PID
    AliEventCuts    *fEventCuts;
    TH2F            *fHistChKaons;
    TH2F            *fHistEtaCuts;
    TH2F            *fHistPosKaons;
    TH2F            *fHistNegKaons;
    TH1F            *fHistV0M;
//    TH1F            *fHistV0s;
    TH2F            *fHistK0Pairs;
    
    Double_t        fLpTCut;        //not a pointer???
	Double_t        fUpTCut;
    Double_t        fEtaCut;
    Double_t        fSigCut;
    Double_t        fDecayLv0Cut;
    Double_t        fLpTv0Cut;
    Double_t        fUpTv0Cut;
    Double_t        fEtav0Cut;
    Double_t        fDcaPosToPrimVtxv0Cut;
    Double_t        fDcaNegToPrimVtxv0Cut;
    Double_t        fEtaPosv0Cut;
    Double_t        fEtaNegv0Cut;
    Double_t        fCosPACut;
    Double_t        fSigPosv0Cut;
    Double_t        fSigNegv0Cut;
    
	AliAnalysisTaskKaon2PC(const AliAnalysisTaskKaon2PC&);            // not implemented
	AliAnalysisTaskKaon2PC& operator=(const AliAnalysisTaskKaon2PC&); // not implemented
	
	ClassDef(AliAnalysisTaskKaon2PC, 1);                              // example of analysis
};

#endif

