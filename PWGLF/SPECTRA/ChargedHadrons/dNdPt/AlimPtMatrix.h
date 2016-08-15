#ifndef AlimPtMatrix_cxx
#define AlimPtMatrix_cxx

// Generation of Correlation Matrix
// by Philipp Luettig, 05.08.2011
// modified by Marco Marquard
class TList;
class TH1F;
class TH2F;
class TRandom3;

class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
// class THnSparse;


#include "AliAnalysisTaskSE.h"

class AlimPtMatrix : public AliAnalysisTaskSE {
	public:
		AlimPtMatrix(const char *name = "AlimPtMatrix");
		virtual ~AlimPtMatrix();

		//     virtual void   ConnectInputData(Option_t *);
		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);
		//     virtual Bool_t AcceptParticle(TParticle *p);

		void SetAliESDtrackCuts(AliESDtrackCuts* esdTrackCuts) {fESDTrackCuts = esdTrackCuts;}
		void SetMaxVertexZ(Float_t vZ) {fMaxVertexZ = vZ;}
		void SetUseCentrality(Int_t cent) {fUseCentrality = cent;}
		void SetUseMC(Bool_t bMC) {bUseMC = bMC;}
		void SetCentLimit(Int_t highCent=9) {fCentLimit = highCent;}
		void SetSubsample(Double_t SubS ) {fSubsample = SubS ;}
		void SetSubsampleEvent(Double_t SubEv ) {fSubsampleEvent = SubEv;}
		void SetSubsampleTrack(Double_t SubTr ) {fSubsampleTrack = SubTr;}
		void SetCutsGenParticle(Double_t ptMin, Double_t ptMax, Double_t etaMin, Double_t etaMax) {genMinPt = ptMin; genMaxPt = ptMax; genMinEta = etaMin; genMaxEta = etaMax; }
		void SetTrigger(const AliTriggerAnalysis::Trigger trigger)    { fTrigger = trigger; }
		AliTriggerAnalysis::Trigger GetTrigger() const                { return fTrigger; }


	private:
		AliESDEvent             *fESD;                  //ESD object
		AliMCEvent              *fMC;                   //MC object
		TList                   *fOutputList;           // List where all the output files are stored
		TRandom3		*r3;			//MC generator

		// Histogram
		TH2F                    *fNaccCent;         // Nacc:Nch:Centrality
		TH2F                    *fNchCent;         // Nacc:Nch:Centrality
		TH2F                    *fNaccNch;         // Nacc:Nch:Centrality
		TH1F		       	*fEvtsCent;		   // Number of Events:Centrality

		AliESDtrackCuts         *fESDTrackCuts;         // Esd track cuts

		//variables
		Float_t                 fMaxVertexZ;            // Maximum value for Vertex Z position
		Int_t                   fUseCentrality;	        // Use centrality (0=off, 1=VZERO)
		Int_t                   iNbEvents;              // number of events in total
		Bool_t                  bUseMC;                 // use mc
		Int_t 			fCentLimit;
		Int_t			fSubsample;		// Only takes a fraction of tracks and or events, (0 = no subsamples, 1 = subsamples of events, 2 = subsamples of tracks, 3 = subsamples of events and tracks)
		Double_t			fSubsampleEvent;	// Define fraction for event subsample
		Double_t			fSubsampleTrack;	// Define fraction for track subsample

		AliTriggerAnalysis::Trigger fTrigger;         // trigger definition MB1, MB2 ...

		// cuts for generated particles
		Double_t                    genMinPt;               // minimum Pt for generated particles
		Double_t                    genMaxPt;               // maximum Pt for generated particles
		Double_t                    genMinEta;              // minimum Eta for generated particles
		Double_t                    genMaxEta;              // maximum Eta for generated particles



		AlimPtMatrix(const AlimPtMatrix&); // not implemented
		AlimPtMatrix& operator=(const AlimPtMatrix&); // not implemented

		ClassDef(AlimPtMatrix, 1); // example of analysis
};

#endif
