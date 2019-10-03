/**************************************************************************
 *     AliXtAnalysis:
 * This class constructs inclusive and isolated (based on charged tracks)
 * invariant spectra. Isolated case uses a different efficiency where the
 * contamination from cases where non-isolated particle appears as an
 * isolated one due to finite detector efficiency is taken into account.
 *
 * contact: Sami Räsänen
 *          University of Jyväskylä, Finland
 *          sami.s.rasanen@jyu.fi
 **************************************************************************/

#ifndef ALIXTANALYSIS_H
#define ALIXTANALYSIS_H

#include <AliAnalysisTaskSE.h>

class AliVEvent;
class TList;
class AliAnalysisUtils;
class AliJCard;
class TDirectory;
class AliJXtHistos;
class TH1D;
class TObjArray;

class AliXtAnalysis : public AliAnalysisTaskSE {

	public:
		AliXtAnalysis();
		AliXtAnalysis(const char *name, const char *cardname);
		AliXtAnalysis(const AliXtAnalysis& a);
		AliXtAnalysis& operator=(const AliXtAnalysis& ap);
		
		virtual ~AliXtAnalysis(); 

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);

		bool IsGoodEvent(AliVEvent *event);
		void SetDebugMode( int debug) { fDebugMode = debug; };

	private:
		TList               *fOutputList;	// Output list
		AliAnalysisUtils	*fAnaUtils;     // Analysis utils
		
        //AliIsolatedEfficiency	*fEfficiency;	// TODO: isolated efficiency, currently eff = 1
        AliJCard            *fCard;		// All parameters
    
        TDirectory          *fHistDir;  // histograms to the directory, structure compatible with JCORRAN
        AliJXtHistos 		*fHistos;	// encapsulate (almost) all histograms
        TH1D                *fhEvents;  // general event statistics
    
		TObjArray           *fChargedList;	// accepted charged tracks in an event
		TObjArray           *fIsolatedChargedList;	// isolated charged tracks in an event
		
        int fCentBin;       // centrality bin
        int fZBin;          // z -vertex bin
        double fZVert;      // z -vertex
    
		Int_t fevt;     // event number
		int fDebugMode; // debug mode
				
		ClassDef(AliXtAnalysis, 1);
};

#endif
