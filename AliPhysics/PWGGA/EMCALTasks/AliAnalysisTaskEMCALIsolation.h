#ifndef AliAnalysisTaskEMCALIsolation_cxx
#define AliAnalysisTaskEMCALIsolation_cxx

// Task for isolating gammas with EMCAL
// Author: Marco Marquard

class TH1F;
class AliESDEvent;
class AliStack;
class AliESDCaloCells;
class AliEMCALGeometry;

#include "TString.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEMCALIsolation : public AliAnalysisTaskSE {
	public:
		AliAnalysisTaskEMCALIsolation();
		AliAnalysisTaskEMCALIsolation(const char *name);
		virtual ~AliAnalysisTaskEMCALIsolation();

		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);

		void	SetVerbose(Bool_t b)			{bVerbose	= b;}
		void	SetMC(Bool_t mc)			{bMC	= mc;}

	protected:

		const char * GetParticleName(Int_t);	

	private:
		Bool_t		bVerbose;			// Verbose option
		Bool_t		bMC;				// MC option
		AliESDEvent	*fESD;				// ESD object
		AliAODEvent	*fAOD;				// AOD object
		AliMCEvent	*fMC;				// MC Event
		AliStack	*fStack;			// Ali stack
		TList		*fOutputList;			// Output list
		TH2F		*fHistGlobalHmap;		// Cell hit map for the complete EMCAL
		TH2F		*fHistGlobalHmap0;		// Cell hit map for the complete EMCAL
		TTree		*fTreeEvent;			//Tree with event informations
		TTree		*fTreeCluster;			//Tree with cluster informations
		Int_t		emclus;				//number of cluster per event
		Double_t	prodrad;			//production radius of V0 vertex
		Int_t		contPID;			//PID of contributor
		Int_t		mothPID;			//PID of contributor mother
		Bool_t		trackmatch;				//Track matching

		AliEMCALGeometry	*fGeom;			// geometry utils
		TString			fGeoName;		// geometry name (def = EMCAL_COMPLET (alternative: EMCAL_FIRSTYEARV1))
		AliESDCaloCells	*fESDCells;			//!pointer to esd cells
		AliAODCaloCells	*fAODCells;			//!pointer to aod cells

		TObjArray	*fEsdClusters;			//!pointer to esd clusters
		TObjArray	*fAodClusters;			//!pointer to aod clusters


		AliAnalysisTaskEMCALIsolation(const AliAnalysisTaskEMCALIsolation&); // not implemented
		AliAnalysisTaskEMCALIsolation& operator=(const AliAnalysisTaskEMCALIsolation&); // not implemented

		ClassDef(AliAnalysisTaskEMCALIsolation, 1); // example of analysis
};

#endif
