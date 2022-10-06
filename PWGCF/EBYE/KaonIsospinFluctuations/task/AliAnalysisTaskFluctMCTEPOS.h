#ifndef ALIANALYSISTASKFLUCTMCTEPOS_H
#define ALIANALYSISTASKFLUCTMCTEPOS_H

class TH1D;
class TH2F;
class TH3F;
class TString;
class AliESDEvent;
class AliESDTrack;
class AliESDMCParticle;
class TList;
class AliESDtrackCuts;




#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"


class AliAnalysisTaskFluctMCTEPOS : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskFluctMCTEPOS();
    AliAnalysisTaskFluctMCTEPOS(const char *name);
    virtual ~AliAnalysisTaskFluctMCTEPOS();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
		
	void	SetIsMonteCarlo			(bool isMonteCarlo = true)			{fIsMonteCarlo = isMonteCarlo;}

			
	
	
 private:
	
	bool	fIsMonteCarlo;

	
	TList           *fOutput;        // Output list
	AliPIDResponse	*fPIDResponse;	 // PID

	
	AliESDEvent     *fESD;                   //! AOD event                        
	AliESDVertex    *fPrimaryVtx;            //! AOD vertex                       

	

       const UShort_t  fTrackBuffSize;  

	TH1F            *fHistGoodEvent;        // Pt spectrum
	TH1F            *fHistPt;        // Pt spectrum
	TH1F            *fHistEta;       // pseudorapidity spectrum

	TH1F			*fHistNV0;
	TH2D			*fHistBBK0Pos;		//PID of the positive daughter of K0 candidates
	
    	TH2D			*fHistBBK0Neg;		//PID of the negative daughter of K0 candidates
	TH2D			*fHistBBPion;		//PID of the negative daughter of K0 candidates
	
	
	TH2F			*fHistZVertexCent;	 //	Z coordinate of primary vertex
	
	
	TH2D			*fHistCosPaMK0;		//	Transverse momentum distribution vs CosPa for K0Short Candidates
	
	
	TH2D			*fHistcTauMK0;		//	Transverse momentum distribution vs cTau for K0Short Candidates
	
	
	TH2D			*fHistDcaMK0;		//	Transverse momentum distribution vs Dca for K0Short Candidates
	
	
	TH2D			*fHistRapMK0;		//	Transverse momentum distribution vs Rap for K0Short Candidates
	
	
	
	TH2D			*fHistArmPodK0;		//Armenteros plot for K0 candidates.


	Double_t          fMCImpactParameter;  //---activate For impact parameter
	THnSparseD *fHistoCorrelation;
	THnSparseD *fKshortSparse;
	Float_t            fCentrality;       //---activate for centrality

	
	// NEW HISTO to be declared here
	
	AliAnalysisTaskFluctMCTEPOS(const AliAnalysisTaskFluctMCTEPOS&); // not implemented
	AliAnalysisTaskFluctMCTEPOS& operator=(const AliAnalysisTaskFluctMCTEPOS&); // not implemented
	
	ClassDef(AliAnalysisTaskFluctMCTEPOS, 1); // example of analysis
};

#endif
