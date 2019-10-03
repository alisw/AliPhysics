#ifndef AliJXtAnalysis_cxx
#define AliJXtAnalysis_cxx

#include <TString.h>
#include <AliAnalysisTaskSE.h>
#include <TClonesArray.h>

#include "AliJHistManager.h"

class AliJEfficiency;

class AliJXtAnalysis : public AliAnalysisTaskSE{

	public:
		AliJXtAnalysis();
		AliJXtAnalysis(const char *name);
		AliJXtAnalysis(const AliJXtAnalysis& a);
		AliJXtAnalysis& operator=(const AliJXtAnalysis& ap);

		virtual ~AliJXtAnalysis();
		
		virtual void   UserCreateOutputObjects();
		virtual void   Init();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);

		void SetInputList( TClonesArray *inputarray){fInputList = inputarray; };
		void SetEventCentrality( float cent ){fCent = cent; };
		void SetEventVertex( double *vtx ){ fVertex = vtx; };
		void SetEtaRange( double etaRange ){fetaRange = etaRange; };
		void SetZVertexRange( double zvertRange ){ fzvertexRange = zvertRange; };
		void SetSqrts( double sqrts ){ fsqrts = sqrts;} ;
		void SetMinIsolPt( double minPt ){ fMinIsolPt = minPt; };
		void SetIsolationRadius( double R ){ fisolR = R; };
		
		void SetRunNumber(int runN){ frunNumber = runN; };
		void SetEfficiencyConfig(int mode, int filterBit){ fEffMode = mode; fTrackFilterBit = filterBit;};
		
		void SetDebugLevel( int dblv ){fDebugLevel = dblv;};
		inline void DEBUG(int level, TString msg){if(level<fDebugLevel) std::cout<<level<<"\t"<<msg<<std::endl;};
		
		void SetEfficiency( AliJEfficiency *effInput ){ fEfficiency = effInput; };
		void SetEfficiencyIsolated( AliJEfficiency *effInput ){ fEfficiencyIsolated = effInput; };
		void SetTrackFilterBit( int inputBit ){ fTrackFilterBit = inputBit; };
		void SetEfficiencyFilterBit( int inputBit ) { fEffFilterBit = inputBit; };
			
		double GetEtaRange(){ return fetaRange; };
		double GetZVertexRange(){ return fzvertexRange; };
		double GetSqrts(){ return fsqrts; };
		double GetMinIsolPt(){ return fMinIsolPt; };
		double GetIsolationRadius(){ return fisolR; };
		
		// Fill histogram
		void FillEfficiencyCheckHistogram(double pT, double effCorr);
		void FillIsolatedEfficiencyCheckHistogram(double pT, double effCorr);
		void FillControlHistograms(double pT, double effCorr, int centBin);
		void FillInclusiveHistograms(double pT, double xT, double eta, double phi, double effCorr, int centBin);
		void FillIsolatedHistograms(double pT, double xT, double eta, double phi, double effCorr, int centBin);
		void FillInclusiveConeActivities(double pT, double sumPt); // Fill TProfile
		void FillIsolatedConeActivities(double pT, double sumPt);  // Fill TProfile
		
	private:
		Long64_t 	 AnaEntry;	// comment me 
		TClonesArray	*fInputList;	// input tracks for chosen track cuts, given by the task
		double 		*fVertex;	//! 3D-vertex (x,y,z)
		Float_t		fCent;		// centrality of the event
		int		fDebugLevel;	// debugging

		int 	fNCent;     // number of centrality bins
		int 	fCBin;      // centrality bin
		double	*fCentBin;  //! centrality binning

		int	fNJacek;    // number of bins in published data
		double	*fPttJacek; //! bin borders in the published data
		
		int	fnBinsPt;    // number of pT bins
		double	*fLogBinsPt; //! pT binning
		
		int	fnBinsXt;    // number of xT bins
		double	*fLogBinsXt; //! xT binning
		
		int frunNumber;	     // run number for the efficiency class
		
		double fsqrts;	      // collision energy
		double fetaRange;     // task does not cut eta-range, do it here
		double fisolR;	      // isolation radius
		double fzvertexRange; // task has very loose vertex cut, tighter in here 
		
		double fMinIsolPt;    // minimum pT of the track to be accepted as an isolated one
		
		AliJEfficiency * fEfficiency;         // non-isolated efficiency
		AliJEfficiency * fEfficiencyIsolated; // isolated efficiency
		int		 fEffMode;            // efficiency: 0:NoEff 1:Period 2:RunNum 3:Auto
		int 		 fTrackFilterBit;     // track selection for efficiency 
		int		 fEffFilterBit;	      // JEfficiency correspondance to track bit
			
		AliJHistManager * fHMG; // output histograms encapsulated here

		AliJBin fHistCentBin; // centrality binning
		AliJBin fVertexBin;   // x, y, z 
		
		AliJTH1D fh_cent;    // for cent dist
		AliJTH1D fh_ntracks; // for number of tracks dist
		
		AliJTProfile fh_effCorr; // to check the efficiency correction	
		AliJTProfile fh_effCorrIsolated; // to check the isolated efficiency correction
		
		// Following histograms have centrality binning
		
		AliJTH1D fh_vertex;  // for z-vertex
		
		AliJTH1D fh_ptJacek; // to compare with the published data
		
		AliJTH1D fh_charPt;      // charged pT in bins of centrality
		AliJTH1D fh_invCharPt;   // invariant charged pT in bins of centrality
		AliJTH1D fh_charPtNoEff; // charged pT, no efficiency correction (for checking)
		AliJTH1D fh_charEta;	 // charged eta
		AliJTH1D fh_charPhi;	 // charged phi
			
		AliJTH1D fh_charXt;      // charged xT in bins of centrality
		AliJTH1D fh_invCharXt;   // invariant charged xT in bins of centrality
		
		AliJTH1D fh_isolCharPt;      // isolated charged pT in bins of centrality
		AliJTH1D fh_isolInvCharPt;   // isolated invariant charged pT in bins of centrality
		AliJTH1D fh_isolCharPtNoEff; // isolated charged pT, no efficiency correction (for checking)
		AliJTH1D fh_isolCharEta;     // isolated charged eta
		AliJTH1D fh_isolCharPhi;     // isolated charged phi

		AliJTH1D fh_isolCharXt;      // isolated charged xT in bins of centrality
		AliJTH1D fh_isolInvCharXt;   // isolated invariant charged xT in bins of centrality
		
		AliJTProfile fh_coneActivity;		// pT sum in cone, to be compared to the ALICE UE results
		AliJTProfile fh_coneActivityIsolated;	// activity for isolated particles
		
		ClassDef(AliJXtAnalysis, 1); // example of analysis
};
#endif
