#ifndef ALIANALYSISTASKKAONISOSPINFLUCT_H
#define ALIANALYSISTASKKAONISOSPINFLUCT_H

class TH1D;
class TH2F;
class TH3F;
class TString;
class AliAODEvent;
class AliAODTrack;
class AliAODMCParticle;
class TList;
class AliESDtrackCuts;





#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include "TRandom3.h"


class AliAnalysisTaskKaonIsospinFluct : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskKaonIsospinFluct();
    AliAnalysisTaskKaonIsospinFluct(const char *name);
    virtual ~AliAnalysisTaskKaonIsospinFluct();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
		

	void	SetcutCosPa				(double cutCosPa = 0.999)			{fcutCosPa = cutCosPa;}
	void	SetcutcTauMin			(double cutcTauMin = -999)			{fcutcTauMin = cutcTauMin;}
	void	SetcutNcTauMax			(double cutNcTauMax = 3.0)			{fcutNcTauMax = cutNcTauMax;}
	void	SetcutBetheBloch		(double cutBetheBloch = 3.0)		{fcutNsigma = cutBetheBloch;}

	void	SetcutEta				(double cutEta = 0.8)				{fcutEta = cutEta;}
	void	SetcutRapidity			(double cutRapidity = 0.5)			{fcutRapidity = cutRapidity;}
	void	SetcutArmenteros		(double cutArmenteros = 0.22)		{fcutArmenteros = cutArmenteros;}
	void    SetcutDCA                        (double cutDCA = 1.0)
	{fcutDCA = cutDCA;}

	Int_t ProcessHybrid(AliAODTrack *track);
	Int_t ProcessHybridPro(AliAODTrack *track);
	Int_t ProcessTPC(AliAODTrack *track);
	
	



		
	
	
 private:
	
	
	double	fcutCosPa;
	double	fcutcTauMin;
	double	fcutNcTauMax;
	double	fcutNsigma;
	double	fcutEta;
	double	fcutRapidity;
	double	fcutArmenteros;
        double  fcutDCA;
	
	TList           *fOutput;        // Output list
	AliPIDResponse	*fPIDResponse;	 // PID

	
	AliAODEvent     *fAOD;                   //! AOD event                        
	AliAODVertex    *fPrimaryVtx;            //! AOD vertex                       

	

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


	TH2D    *fdEdX1;
	TH2D    *fdEdX2;
	TH2D    *fdEdX3;
	TH2D    *fdEdXPID;
	TH2D    *fdEdXnoPID;


	TH2D    *fnsigmakaon;
	TH2D    *fnsigmaproton;
	TH2D    *fnsigmapion;
	
	THnSparseD *fHistoCorrelation;
	THnSparseD *fKshortSparse;
	
	// NEW HISTO to be declared here
	
	AliAnalysisTaskKaonIsospinFluct(const AliAnalysisTaskKaonIsospinFluct&); // not implemented
	AliAnalysisTaskKaonIsospinFluct& operator=(const AliAnalysisTaskKaonIsospinFluct&); // not implemented
	
	ClassDef(AliAnalysisTaskKaonIsospinFluct, 1); // example of analysis
};

#endif

