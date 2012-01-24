#ifndef AliAnalysisTaskPtFluc_cxx
#define AliAnalysisTaskPtFluc_cxx

// Analysis of Pt FLuctuations (pp)
// Author: Stefan Heckel
// Version of pp task:   8.0, 19.04.2011


class TList;
class TH1F;
class TProfile;

class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;


#include "AliAnalysisTaskSE.h"


class AliAnalysisTaskPtFluc : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskPtFluc(const char *name = "AliAnalysisTaskPtFluc");
    virtual ~AliAnalysisTaskPtFluc();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetAliESDtrackCuts(AliESDtrackCuts* esdTrackCuts) {fESDTrackCuts = esdTrackCuts;}
    void SetMaxVertexZ(Float_t vZ) {fMaxVertexZ = vZ;}
    void SetNContributors(Int_t nCont) {fNContributors = nCont;}
    void SetMC(Bool_t bMC) {fMC = bMC;}


  private:
    AliESDEvent	     *fESD;		// ESD object
    TList*	      fOutputList;	// List where all the output files are stored
	TH1F*		fPtSpec;	// Pt spectrum
	TH1F*		fMult;		// Multiplicity distribution
	TH1F*		fEta;		// Eta distribution
	TH1F*		fEtaPhiPlus;	// Phi distribution for positive eta
	TH1F*		fEtaPhiMinus;	// Phi distribution for negative eta
	TH1F*		fVtxZ;		// Vertex Z distribution after physics selection before any further cuts
	TH1F*		fVtxZCut;	// Vertex Z dist. after vertex Z cut
	TH1F*		fVtxZCont;	// Vertex Z dist. after vertex cut on nContributors
	TH1F*		fVtxZTrackCuts;	// Vertex Z dist. after all event and track cuts
	TH1F*		fEventMeanPt;		// Event mean pT distribution
	TH1F*		fEventMeanPtSq;		// Event mean pT squared dist.
	TH1F*		fMultEventMeanPt;	// Event mean pT for multiplicity bins
	TH1F*		fMultEventMeanPtSq;	// Event mean pT squared for mult. bins
	TH1F*		fTwoPartCorrEv;		// Two-particle correlator for multiplicity bins
	TH1F*		fTwoPartCorrEvSq;	// Two-part. corr. squared for mult. bins
	TH1F*		fTwoPartCorrEvSample;	// Two-part. corr. for the whole sample
	TH1F*		fTwoPartCorrEvSampleSq;	// Two-part. corr. squared for the whole sample

    AliESDtrackCuts* fESDTrackCuts;	// Esd track cuts
    Float_t          fMaxVertexZ;	// Maximum value for Vertex Z position
    Int_t            fNContributors;	// Minimum contributors to the vertex
    Bool_t           fMC;		// Check for MC

    AliAnalysisTaskPtFluc(const AliAnalysisTaskPtFluc&); // not implemented
    AliAnalysisTaskPtFluc& operator=(const AliAnalysisTaskPtFluc&); // not implemented

    ClassDef(AliAnalysisTaskPtFluc, 1);

};

#endif
