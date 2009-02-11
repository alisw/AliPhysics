#ifndef ALIANALYSISTASKESDCHECKCASCADE_H
#define ALIANALYSISTASKESDCHECKCASCADE_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskESDCheckCascade class
//            (AliAnalysisTaskESDCheckCascade)
//            This task is for QAing the Cascades from the ESD
//              Origin:  AliAnalysisTaskESDCheckV0 by B.H. Nov2007, hippolyt@in2p3.fr
//            Adapted to Cascade : A.M Mar2008, antonin.maire@ires.in2p3.fr
//-----------------------------------------------------------------

class TList;
class TH1F;
class AliESDEvent;

#include "AliAnalysisTask.h"

class AliAnalysisTaskESDCheckCascade : public AliAnalysisTask {
 public:
  AliAnalysisTaskESDCheckCascade();
  AliAnalysisTaskESDCheckCascade(const char *name);
  virtual ~AliAnalysisTaskESDCheckCascade() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESDInputHandler *fesdH;		//! InputHandler object
  AliMCEvent  *fMCevent;		//! MC event object
  AliESDEvent *fESD;			//! ESD object
  
 
   	
	

		TList	*fListHistCascade;		//! List of Cascade histograms
	
	// - General histos (filled for any event)
	TH1F	*fHistTrackMultiplicity;		//! Track multiplicity distribution
	TH1F	*fHistCascadeMultiplicity;		//! Cascade multiplicity distribution


	// - Vertex Positions
	TH1F	*fHistVtxStatus;			//! Is there a tracking vertex in the event ?

		// Vtx coming from the full tracking
	TH1F	*fHistPosTrkgPrimaryVtxX;  		//! primary vertex position distribution in x 
	TH1F	*fHistPosTrkgPrimaryVtxY;		//! primary vertex position distribution in y
	TH1F	*fHistPosTrkgPrimaryVtxZ;  		//! primary vertex position distribution in z
	TH1F	*fHistTrkgPrimaryVtxRadius;		//! primary vertex (3D) radius distribution 

		// Best primary Vtx available for the event
	TH1F	*fHistPosBestPrimaryVtxX;  		//! (best) primary vertex position distribution in x 
	TH1F	*fHistPosBestPrimaryVtxY;		//! (best) primary vertex position distribution in y
	TH1F	*fHistPosBestPrimaryVtxZ;  		//! (best) primary vertex position distribution in z
	TH1F	*fHistBestPrimaryVtxRadius;		//! (best) primary vertex radius distribution 
	
		// Correlation Best Vtx / Full Tracking Vtx
	TH2F	*f2dHistTrkgPrimVtxVsBestPrimVtx;	//!  Radius of prim. Vtx from tracks Vs Radius of best Prim. Vtx
	
	
	
	// - Typical histos on the variables used for the selection of cascades
	TH1F	*fHistEffMassXi;      			//! reconstructed cascade effective mass
	TH1F	*fHistChi2Xi;         			//! chi2 value
	TH1F	*fHistDcaXiDaughters; 			//! dca between Xi's daughters
	TH1F 	*fHistDcaBachToPrimVertex;		//! dca of the bachelor track to primary vertex
	TH1F	*fHistXiCosineOfPointingAngle;		//! cosine of Xi pointing angle in a cascade
	TH1F	*fHistXiRadius;				//! (transverse) radius of the cascade vertex 
		
	// - Histos about ~ the "V0 selection part" of the cascade,  coming by inheritance from AliESDv0
	TH1F	*fHistMassLambdaAsCascDghter;		//! Test Invariant Mass of Lambda coming from Cascade
	TH1F	*fHistV0Chi2Xi;				//! V0 chi2 distribution, for the V0 associated to a cascade
	TH1F	*fHistDcaV0DaughtersXi;			//! Dca between V0 daughters, for the V0 associated to a cascade
	TH1F	*fHistDcaV0ToPrimVertexXi;		//! Dca of V0 to primary vertex, for the V0 associated to a cascade	
	TH1F	*fHistV0CosineOfPointingAngleXi;	//! Cosine of V0 pointing angle, for the V0 associated to a cascade
	TH1F	*fHistV0RadiusXi;			//! V0 (transverse) distance distribution, for the V0 associated to a cascade

	TH1F	*fHistDcaPosToPrimVertexXi;		//! Dca of V0 positive daughter to primary vertex, for the V0 associated to a cascade
	TH1F	*fHistDcaNegToPrimVertexXi;		//! Dca of V0 negative daughter to primary vertex, for the V0 associated to a cascade
	

	// - Effective mass histos for cascades.
	TH1F	*fHistMassXiMinus;			//! reconstructed cascade effective mass, under Xi- hyp.
	TH1F	*fHistMassXiPlus;			//! reconstructed cascade effective mass, under Xi+ hyp.
	TH1F	*fHistMassOmegaMinus;			//! reconstructed cascade effective mass, under Omega- hyp.
	TH1F	*fHistMassOmegaPlus;			//! reconstructed cascade effective mass, under Omega+ hyp.

	// - Complements for QA
	TH1F	*fHistXiTransvMom;     			//! Xi transverse momentum 
	TH1F	*fHistXiTotMom;     			//! Xi momentum norm
	
	TH1F	*fHistBachTransvMom;   			//! bachelor transverse momentum 
	TH1F	*fHistBachTotMom;  			//! bachelor momentum norm
				
	TH1F	*fHistChargeXi;				//! Charge sign of the cascade candidate
	TH1F	*fHistV0toXiCosineOfPointingAngle;	//! Cos. of Pointing angle between the V0 mom and the Xi-V0 vtx line
  
	TH1F	*fHistRapXi;				//! rapidity of Xi candidates
	TH1F	*fHistRapOmega;				//! rapidity of Omega candidates
	TH1F	*fHistEta;				//! eta distrib. of all the cascade candidates
	TH1F	*fHistTheta;				//! theta distrib. of all the cascade candidates
	TH1F	*fHistPhi;				//! phi distrib. of all the cascade candidates
	
	TH2F	*f2dHistArmenteros;			//! alpha(casc. cand.) Vs PtArm(casc. cand.)
	TH2F	*f2dHistEffMassXiVsEffMassLambda;	//! Xi- Eff mass Vs V0 Eff mass, under Xi- hyp.
	TH2F	*f2dHistXiRadiusVsEffMassXi;		//! transv. casc. decay radius Vs Xi- Eff mass, under Xi- hyp.

  AliAnalysisTaskESDCheckCascade(const AliAnalysisTaskESDCheckCascade&);            // not implemented
  AliAnalysisTaskESDCheckCascade& operator=(const AliAnalysisTaskESDCheckCascade&); // not implemented
  
  ClassDef(AliAnalysisTaskESDCheckCascade, 2);
};

#endif
