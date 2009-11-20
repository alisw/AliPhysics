#ifndef ALIANALYSISTASKCHECKCASCADE_H
#define ALIANALYSISTASKCHECKCASCADE_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskCheckCascade class
//            (AliAnalysisTaskCheckCascade)
//            This task has four roles :
//              1. QAing the Cascades from ESD and AOD
//                 Origin:  AliAnalysisTaskESDCheckV0 by B.H. Nov2007, hippolyt@in2p3.fr
//              2. Prepare the plots which stand as raw material for yield extraction (wi/wo PID)
//              3. Supply an AliCFContainer meant to define the optimised topological selections
//              4. Rough azimuthal correlation study (Eta, Phi)
//            Adapted to Cascade : A.Maire Mar2008, antonin.maire@ires.in2p3.fr
//            Modified :           A.Maire Nov2009, antonin.maire@ires.in2p3.fr
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;
 
class AliTPCpidESD;
class AliCFContainer;


#include "TString.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCheckCascade : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskCheckCascade();
  AliAnalysisTaskCheckCascade(const char *name);
  virtual ~AliAnalysisTaskCheckCascade() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
	  void   DoAngularCorrelation(const Char_t   *lCascType, 
				            Double_t  lInvMassCascade, 
				      const Int_t    *lArrTrackID,
				            TVector3 &lTVect3MomXi, 
				            Double_t  lEtaXi);
  virtual void   Terminate(Option_t *);
  
  void SetCollidingSystems(Short_t collidingSystems = 0)     {fCollidingSystems = collidingSystems;}
  void SetAnalysisType    (const char* analysisType = "ESD") {fAnalysisType = analysisType;}
    
 private:
  	TString 	fAnalysisType;			// "ESD" or "AOD" analysis type	
	Short_t 	fCollidingSystems;		// 0 = pp collisions or 1 = AA collisions
	AliTPCpidESD*	fTpcPidManager;			//! Tool data member to manage the TPC Bethe-Bloch info
	
	
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
	
	
	// PART 1 : Adavanced QA
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
	
	TH1F	*fHistMassWithCombPIDXiMinus;		//! reconstructed Xi- effective mass, with bach. comb PID
	TH1F	*fHistMassWithCombPIDXiPlus;		//! reconstructed Xi+ effective mass, with bach. comb PID
	TH1F	*fHistMassWithCombPIDOmegaMinus;	//! reconstructed Omega- effective mass, with bach. comb PID
	TH1F	*fHistMassWithCombPIDOmegaPlus;		//! reconstructed Omega+ effective mass, with bach. comb PID

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
	
	TH2F	*f2dHistEffMassLambdaVsEffMassXiMinus;	//! Xi- Eff mass Vs V0 Eff mass, under Xi- hyp.
	TH2F	*f2dHistEffMassXiVsEffMassOmegaMinus;	//! Xi- Eff mass Vs Omega- Eff mass, for negative cascades
	TH2F	*f2dHistEffMassLambdaVsEffMassXiPlus;	//! Xi+ Eff mass Vs V0 Eff mass, under Xi+ hyp.
	TH2F	*f2dHistEffMassXiVsEffMassOmegaPlus;	//! Xi+ Eff mass Vs Omega+ Eff mass, for positive cascades
	
	TH2F	*f2dHistXiRadiusVsEffMassXiMinus;	//! transv. casc. decay radius Vs Xi- Eff mass, under Xi- hyp.
	TH2F	*f2dHistXiRadiusVsEffMassXiPlus;	//! transv. casc. decay radius Vs Xi+ Eff mass, under Xi+ hyp.
	TH2F	*f2dHistXiRadiusVsEffMassOmegaMinus;	//! transv. casc. decay radius Vs Omega- Eff mass, under Omega- hyp.
	TH2F	*f2dHistXiRadiusVsEffMassOmegaPlus;	//! transv. casc. decay radius Vs Omega+ Eff mass, under Omega+ hyp.
	
	
	// PART 2 : TH3F needed for pt spectrum and yield extraction
	// Without any PID
	TH3F	*f3dHistXiPtVsEffMassVsYXiMinus;        //! casc. transv. momemtum Vs Xi- Eff mass Vs Y
	TH3F	*f3dHistXiPtVsEffMassVsYXiPlus;         //! casc. transv. momemtum Vs Xi+ Eff mass Vs Y
	TH3F	*f3dHistXiPtVsEffMassVsYOmegaMinus;     //! casc. transv. momemtum Vs Omega- Eff mass Vs Y
	TH3F	*f3dHistXiPtVsEffMassVsYOmegaPlus;      //! casc. transv. momemtum Vs Omega+ Eff mass Vs Y
	
	// With single PID : proton PID for Xi (pion = useless) / bachelor PID for Omega
	// = a priori, the most important PID info for each species
	TH3F	*f3dHistXiPtVsEffMassVsYWithCombPIDXiMinus;     //! casc. transv. momemtum Vs Xi- Eff mass Vs Y
	TH3F	*f3dHistXiPtVsEffMassVsYWithCombPIDXiPlus;      //! casc. transv. momemtum Vs Xi+ Eff mass Vs Y
	TH3F	*f3dHistXiPtVsEffMassVsYWithCombPIDOmegaMinus;  //! casc. transv. momemtum Vs Omega- Eff mass Vs Y
	TH3F	*f3dHistXiPtVsEffMassVsYWithCombPIDOmegaPlus;   //! casc. transv. momemtum Vs Omega+ Eff mass Vs Y
	
	// With double PID : proton PID + bachelor PID for Omega
	// = "second order" refinement for omegas
	TH3F	*f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaMinus; //! casc. transv. momemtum Vs Omega- Eff mass Vs Y
	TH3F	*f3dHistXiPtVsEffMassVsYWith2CombPIDOmegaPlus;  //! casc. transv. momemtum Vs Omega+ Eff mass Vs Y
	
	// With TPC PID : 3-sigma band on the Bethe-Bloch curve
	// = directly detector-based
	TH3F	*f3dHistXiPtVsEffMassVsYWithTpcPIDOmegaMinus;   //! casc. transv. momemtum Vs Omega- Eff mass Vs Y
	
	// Compilation of all PID plots (3D = casc. transv. momemtum Vs Casc Eff mass Vs Y), stored into an AliCFContainer
	AliCFContainer  *fCFContCascadePIDXiMinus;      //! for Xi-   : Container to store any 3D histos with the different PID flavours
	AliCFContainer  *fCFContCascadePIDXiPlus;       //! for Xi+   : Container to store any 3D histos with the different PID flavours
	AliCFContainer  *fCFContCascadePIDOmegaMinus;   //! for Omega-: Container to store any 3D histos with the different PID flavours
	AliCFContainer  *fCFContCascadePIDOmegaPlus;    //! for Omega+: Container to store any 3D histos with the different PID flavours
	
	
	
	// PART 3 : Towards the optimisation of topological selections
	AliCFContainer  *fCFContCascadeCuts;		//! Container meant to store all the relevant distributions corresponding to the cut variables
	
	
	// PART 4 :  Azimuthal correlation study
	THnSparseF	*fHnSpAngularCorrXiMinus;	//! Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks Vs Eff Mass
	THnSparseF	*fHnSpAngularCorrXiPlus;	//! Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks Vs Eff Mass
	THnSparseF	*fHnSpAngularCorrOmegaMinus;	//! Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks Vs Eff Mass
	THnSparseF	*fHnSpAngularCorrOmegaPlus;	//! Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks Vs Eff Mass


  AliAnalysisTaskCheckCascade(const AliAnalysisTaskCheckCascade&);            // not implemented
  AliAnalysisTaskCheckCascade& operator=(const AliAnalysisTaskCheckCascade&); // not implemented
  
  ClassDef(AliAnalysisTaskCheckCascade, 7);
};

#endif
