#ifndef AliAnalysisTaskPtFlucPbPb_cxx
#define AliAnalysisTaskPtFlucPbPb_cxx

// Analysis of Pt FLuctuations (PbPb)
// Author: Stefan Heckel
// Version of PbPb task: 5.0, 18.04.2011


class TList;
class TH1F;
class TProfile;

class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;


#include "AliAnalysisTaskSE.h"


class AliAnalysisTaskPtFlucPbPb : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskPtFlucPbPb(const char *name = "AliAnalysisTaskPtFlucPbPb");
    virtual ~AliAnalysisTaskPtFlucPbPb();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetAliESDtrackCuts(AliESDtrackCuts* esdTrackCuts) {fESDTrackCuts = esdTrackCuts;}
    void SetMaxVertexZ(Float_t vZ) {fMaxVertexZ = vZ;}
    void SetNContributors(Int_t nCont) {fNContributors = nCont;}
    void SetUseCentrality(Int_t cent) {fUseCentrality = cent;}
    void SetMC(Bool_t bMC) {fMC = bMC;}


  private:
    AliESDEvent	     *fESD;		// ESD object
    TList*	      fOutputList;	// List where all the output files are stored
	TH1F*		fPtSpec;	// Pt spectrum
	TH1F*		fMult;		// Multiplicity distribution
	TH1F*		fMultSum;	// Sum of number of tracks of all events for each multiplicity bin
	TH1F*		fMultSumPt;	// Sum of pTs for multiplicity bins
	TH1F*		fMultNrPairs;	// Sum of number of pairs in mult. bins
	TH1F*		fCent;		// Centrality distribution
	TH1F*		fCentSum;	// Sum of number of tracks of all events for each centrality bin
	TH1F*		fCentSumPt;	// Sum of pTs for centrality bins
	TH1F*		fCentNrPairs;	// Sum of number of pairs in cent. bins
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
	TH1F*		fCentEventMeanPt;	// Event mean pT for centrality bins
	TH1F*		fCentEventMeanPtSq;	// Event mean pT squared for cent. bins
	TH1F*		fTwoPartCorrEv;		// Two-particle correlator for multiplicity bins
	TH1F*		fTwoPartCorrEvSq;	// Two-part. corr. squared for mult. bins
	TH1F*		fTwoPartCorrEvCent;	// Two-particle correlator for centrality bins
	TH1F*		fTwoPartCorrEvCentSq;	// Two-part. corr. squared for cent. bins

    AliESDtrackCuts* fESDTrackCuts;	// Esd track cuts
    Float_t          fMaxVertexZ;	// Maximum value for Vertex Z position
    Int_t            fNContributors;	// Minimum contributors to the vertex
    Int_t            fUseCentrality;	// Use centrality (0=off, 1=VZERO, 2=SPD(not yet implemented))
    Bool_t           fMC;		// Check for MC

    AliAnalysisTaskPtFlucPbPb(const AliAnalysisTaskPtFlucPbPb&); // not implemented
    AliAnalysisTaskPtFlucPbPb& operator=(const AliAnalysisTaskPtFlucPbPb&); // not implemented

    ClassDef(AliAnalysisTaskPtFlucPbPb, 1);

};

#endif
