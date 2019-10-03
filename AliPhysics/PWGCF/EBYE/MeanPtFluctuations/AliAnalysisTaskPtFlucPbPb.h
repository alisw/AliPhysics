#ifndef AliAnalysisTaskPtFlucPbPb_cxx
#define AliAnalysisTaskPtFlucPbPb_cxx

// Analysis of Pt Fluctuations (PbPb)
// Author: Stefan Heckel
// Version of PbPb task:  9.2, 03.07.2012


class TList;
class TH1F;
class TH2F;
class TRandom;

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
    void SetMaxVertexZDiff1(Float_t vZDiff1) {fMaxVertexZDiff1 = vZDiff1;}
    void SetNContributors(Int_t nCont) {fNContributors = nCont;}
    void SetUseCentrality(Int_t cent) {fUseCentrality = cent;}
    void SetMC(Bool_t bMC) {fMC = bMC;}
    void SetMCType(Int_t bMCType) {fMCType = bMCType;}
    void SetMCAMPT(Bool_t bMCAMPT) {fMCAMPT = bMCAMPT;}


  private:
    AliESDEvent	     *fESD;		// ESD object
    AliMCEvent	     *fMCev;		// MC object
    TRandom	     *fRandom3;		// Random generator
    TList*	      fOutputList;	// List where all the output files are stored
	TH1F*		fPtSpec;	// Pt spectrum - data or MC truth
	TH1F*		fPtSpec2;	// Pt spectrum 2 - MC ESD
	TH1F*		fMult;		// Multiplicity distribution
	TH1F*		fMultNbins;	// Multiplicity distribution for single mult. bins
	TH1F*		fMultSum;	// Sum of number of tracks of all events for each multiplicity bin
	TH1F*		fMultSumPt;	// Sum of pTs for multiplicity bins
	TH1F*		fMultNrPairs;	// Sum of number of pairs in mult. bins
	TH1F*		fMult1;		// Multiplicity distribution (first bin divided in 4)
	TH1F*		fMultSum1;	// Sum of number of tracks of all events for each multiplicity bin (first bin divided in 4)
	TH1F*		fMultSumPt1;	// Sum of pTs for multiplicity bins (first bin divided in 4)
	TH1F*		fMultNrPairs1;	// Sum of number of pairs in mult. bins (first bin divided in 4)
	TH1F*		fMult10;	// Multiplicity distribution (five bins in 0 < Nacc < 50)
	TH1F*		fMultSum10;	// Sum of number of tracks of all events for each multiplicity bin (five bins in 0 < Nacc < 50)
	TH1F*		fMultSumPt10;	// Sum of pTs for multiplicity bins (five bins in 0 < Nacc < 50)
	TH1F*		fMultNrPairs10;	// Sum of number of pairs in mult. bins (five bins in 0 < Nacc < 50)
	TH1F*		fMult80;	// Multiplicity distribution  -- (only for events with centrality < 80%)
	TH1F*		fMultSum80;	// Sum of number of tracks of all events for each multiplicity bin -- ( " < 80%)
	TH1F*		fMultSumPt80;	// Sum of pTs for multiplicity bins-- ( " < 80%)
	TH1F*		fMultNrPairs80; // Sum of number of pairs in mult. bins -- ( " < 80%)
	TH1F*		fMult801;	// Multiplicity distribution (first bin divided in 4) -- (only for events with centrality < 80%)
	TH1F*		fMultSum801;	// Sum of number of tracks of all events for each multiplicity bin (first bin divided in 4) -- ( " < 80%)
	TH1F*		fMultSumPt801;	// Sum of pTs for multiplicity bins (first bin divided in 4) -- ( " < 80%)
	TH1F*		fMultNrPairs801; // Sum of number of pairs in mult. bins (first bin divided in 4) -- ( " < 80%)
	TH1F*		fMult810;	// Multiplicity distribution (five bins in 0 < Nacc < 50) -- ( " < 80%)
	TH1F*		fMultSum810;	// Sum of number of tracks of all events for each multiplicity bin (five bins in 0 < Nacc < 50)--( " < 80%)
	TH1F*		fMultSumPt810;	// Sum of pTs for multiplicity bins (five bins in 0 < Nacc < 50) -- ( " < 80%)
	TH1F*		fMultNrPairs810; // Sum of number of pairs in mult. bins (five bins in 0 < Nacc < 50) -- ( " < 80%)
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
	TH1F*		fVtxZCutDiff;	// Vertex Z dist. after vertex cut on vtx Z Difference
	TH1F*		fVtxZTrackCuts;	// Vertex Z dist. after all event and track cuts
	TH1F*		fVtxZDiff1;	// Difference 1 between vertex Z distributions
	TH1F*		fVtxZDiff2;	// Difference 2 between vertex Z distributions
	TH1F*		fVtxZDiff3;	// Difference 3 between vertex Z distributions
	TH1F*		fVtxZDiff1b;	// Difference 1 between vertex Z distributions after all cuts
	TH1F*		fVtxZDiff2b;	// Difference 2 between vertex Z distributions after all cuts
	TH1F*		fVtxZDiff3b;	// Difference 3 between vertex Z distributions after all cuts
	TH1F*		fEventMeanPt;		// Event mean pT distribution
	TH1F*		fEventMeanPtSq;		// Event mean pT squared dist.
	TH2F*		fEventMeanPtMult;	// Event mean pT distribution vs. multiplicity (scatter plot)
	TH1F*		fMultEventMeanPt;	// Event mean pT for multiplicity bins
	TH1F*		fMultEventMeanPtSq;	// Event mean pT squared for mult. bins
	TH1F*		fMultEventMeanPtNbins;	// Event mean pT for single mult. bins
	TH1F*		fMultEventMeanPtSqNbins;// Event mean pT squared for single mult. bins
	TH1F*		fCentEventMeanPt;	// Event mean pT for centrality bins
	TH1F*		fCentEventMeanPtSq;	// Event mean pT squared for cent. bins
	TH1F*		fEventMeanPtCent05;	// Event mean pT distribution in cent bin 0-5%
	TH1F*		fEventMeanPtCent2030;	// Event mean pT distribution in cent bin 20-30%
	TH1F*		fEventMeanPtCent7080;	// Event mean pT distribution in cent bin 70-80%
	TH1F*		fTwoPartCorrEv;		// Two-particle correlator for multiplicity bins
	TH1F*		fTwoPartCorrEvSq;	// Two-part. corr. squared for mult. bins
	TH1F*		fTwoPartCorrEv1;	// Two-particle correlator for multiplicity bins (first bin divided in 4)
	TH1F*		fTwoPartCorrEvSq1;	// Two-part. corr. squared for mult. bins (first bin divided in 4)
	TH1F*		fTwoPartCorrEv10;	// Two-particle correlator for multiplicity bins (five bins in 0 < Nacc < 50)
	TH1F*		fTwoPartCorrEvSq10;	// Two-part. corr. squared for mult. bins (five bins in 0 < Nacc < 50)
	TH1F*		fTwoPartCorrEv80;	// Two-particle correlator for multiplicity bins -- ( " < 80%)
	TH1F*		fTwoPartCorrEvSq80;	// Two-part. corr. squared for mult. bins -- ( " < 80%)
	TH1F*		fTwoPartCorrEv801;	// Two-particle correlator for multiplicity bins (first bin divided in 4) -- ( " < 80%)
	TH1F*		fTwoPartCorrEvSq801;	// Two-part. corr. squared for mult. bins (first bin divided in 4) -- ( " < 80%)
	TH1F*		fTwoPartCorrEv810;	// Two-particle correlator for multiplicity bins (five bins in 0 < Nacc < 50) -- ( " < 80%)
	TH1F*		fTwoPartCorrEvSq810;	// Two-part. corr. squared for mult. bins (five bins in 0 < Nacc < 50) -- ( " < 80%)
	TH1F*		fTwoPartCorrEvCent;	// Two-particle correlator for centrality bins
	TH1F*		fTwoPartCorrEvCentSq;	// Two-part. corr. squared for cent. bins

    AliESDtrackCuts* fESDTrackCuts;	// Esd track cuts
    Float_t          fMaxVertexZ;	// Maximum value for Vertex Z position
    Float_t          fMaxVertexZDiff1;	// Maximum value for Vertex Z difference TPC - global
    Int_t            fNContributors;	// Minimum contributors to the vertex
    Int_t            fUseCentrality;	// Use centrality (0=off, 1=VZERO, 2=SPD(not yet implemented))
    Bool_t           fMC;		// Check for MC
    Int_t            fMCType;		// Set MC type: ESD, MC truth (generator level), mod. MC truth
    Bool_t           fMCAMPT;		// Set MC = AMPT or other

    AliAnalysisTaskPtFlucPbPb(const AliAnalysisTaskPtFlucPbPb&); // not implemented
    AliAnalysisTaskPtFlucPbPb& operator=(const AliAnalysisTaskPtFlucPbPb&); // not implemented

    ClassDef(AliAnalysisTaskPtFlucPbPb, 1);

};

#endif
