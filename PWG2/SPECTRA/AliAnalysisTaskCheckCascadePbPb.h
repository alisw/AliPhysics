#ifndef ALIANALYSISTASKCHECKCASCADEPBPB_H
#define ALIANALYSISTASKCHECKCASCADEPBPB_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//            AliAnalysisTaskCheckCascadePbPb class
//              Origin AliAnalysisTaskCheckCascade
//              This task has four roles :
//                1. QAing the Cascades from ESD and AOD
//                   Origin:  AliAnalysisTaskESDCheckV0 by Boris Hippolyte Nov2007, hippolyt@in2p3.fr
//                2. Prepare the plots which stand as raw material for yield extraction (wi/wo PID)
//                3. Supply an AliCFContainer meant to define the optimised topological selections
//                4. Rough azimuthal correlation study (Eta, Phi)
//                Adapted to Cascade : A.Maire Mar2008, antonin.maire@ires.in2p3.fr
//                Modified :           A.Maire Mar2010, antonin.maire@ires.in2p3.fr
//                Modified for PbPb analysis: M. Nicassio Feb 2011, maria.nicassio@ba.infn.it
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;
 
class AliESDpid;
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;


#include "TString.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCheckCascadePbPb : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskCheckCascadePbPb();
  AliAnalysisTaskCheckCascadePbPb(const char *name);
  virtual ~AliAnalysisTaskCheckCascadePbPb();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
/*	  void   DoAngularCorrelation(const Char_t   *lCascType, 
				            Double_t  lInvMassCascade, 
				      const Int_t    *lArrTrackID,
				            TVector3 &lTVect3MomXi, 
				            Double_t  lEtaXi);*/
  virtual Int_t  DoESDTrackWithTPCrefitMultiplicity(const AliESDEvent *lESDevent);
  
  virtual void   Terminate(Option_t *);
  
  void SetCollidingSystems           (Short_t collidingSystems          = 0    ) { fCollidingSystems            = collidingSystems;           }
  void SetAnalysisType               (const char* analysisType          = "ESD") { fAnalysisType                = analysisType;               }
  void SetRelaunchV0CascVertexers    (Bool_t rerunV0CascVertexers       = 0    ) { fkRerunV0CascVertexers       = rerunV0CascVertexers;       }
  void SetQualityCutZprimVtxPos      (Bool_t qualityCutZprimVtxPos      = kTRUE) { fkQualityCutZprimVtxPos      = qualityCutZprimVtxPos;      }
  void SetQualityCutNoTPConlyPrimVtx (Bool_t qualityCutNoTPConlyPrimVtx = kTRUE) { fkQualityCutNoTPConlyPrimVtx = qualityCutNoTPConlyPrimVtx; }
  void SetQualityCutTPCrefit         (Bool_t qualityCutTPCrefit         = kTRUE) { fkQualityCutTPCrefit         = qualityCutTPCrefit;         }
  void SetQualityCut80TPCcls         (Bool_t qualityCut80TPCcls         = kTRUE) { fkQualityCut80TPCcls         = qualityCut80TPCcls;         }
  void SetExtraSelections            (Bool_t extraSelections            = 0    ) { fkExtraSelections            = extraSelections;            }
  void SetCentralityLowLim           (Float_t centrlowlim               = 0.   ) { fCentrLowLim                 = centrlowlim;                }  
  void SetCentralityUpLim            (Float_t centruplim                = 100. ) { fCentrUpLim                  = centruplim;                 }
  void SetCentralityEst              (TString   centrest                = "V0M") { fCentrEstimator              = centrest;                   }
  void SetVertexRange                (Float_t vtxrange                  = 0.   ) { fVtxRange                    = vtxrange;                   }
  void SetUseCFContCascadeCuts       (Bool_t usecfcontcascadecuts       = kFALSE){fUseCFContCascadeCuts         = usecfcontcascadecuts;       }

 private:
        // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
        // your data member object is created on the worker nodes and streaming is not needed.
        // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14


        TString         fAnalysisType;                  // "ESD" or "AOD" analysis type	
        Short_t         fCollidingSystems;              // 0 = pp collisions or 1 = AA collisions

        AliESDpid       *fESDpid;                       // Tool data member to manage the TPC Bethe-Bloch info
        //TPaveText       *fPaveTextBookKeeping;          // TString to store all the relevant info necessary for book keeping (v0 cuts, cascade cuts, quality cuts, ...)

        Bool_t          fkRerunV0CascVertexers;         // Boolean : kTRUE = relaunch both V0 + Cascade vertexers
        Bool_t          fkQualityCutZprimVtxPos;        // Boolean : kTRUE = cut on the prim.vtx  z-position
        Bool_t          fkQualityCutNoTPConlyPrimVtx;   // Boolean : kTRUE = prim vtx should be SPD or Tracking vertex
        Bool_t          fkQualityCutTPCrefit;           // Boolean : kTRUE = ask for TPCrefit for the 3 daughter tracks
        Bool_t          fkQualityCut80TPCcls;           // Boolean : kTRUE = ask for 80 TPC clusters for each daughter track
        Bool_t          fkExtraSelections;              // Boolean : kTRUE = apply tighter selections, before starting the analysis
        Float_t         fCentrLowLim;                   // Lower limit for centrality percentile selection
        Float_t         fCentrUpLim;                    // Upper limit for centrality percentile selection
        TString         fCentrEstimator;                // string for the centrality estimator
        Float_t         fVtxRange;                      // to select events with |zvtx|<fVtxRange cm
        Bool_t          fUseCFContCascadeCuts;          // to use CF containers for topological cut optimization

       
        Double_t        fAlephParameters[5];            // Array to store the 5 param values for the TPC Bethe-Bloch parametrisation
        Double_t        fV0Sels[7];                     // Array to store the 7 values for the different selections V0 related (if fkRerunV0CascVertexers)
        Double_t        fCascSels[8];                   // Array to store the 8 values for the different selections Casc. related (if fkRerunV0CascVertexers)


        TList      *fListHistCascade;              //! List of Cascade histograms
        
        // - General histos (filled before the trigger selection)
        TH1F    *fHistCascadeMultiplicityBeforeEvSel; //! Cascade multiplicity distribution
         
        // - General histos (filled for any triggered event)
        TH1F    *fHistCascadeMultiplicityForCentrEvt;              //! Cascade multiplicity distribution
        TH1F    *fHistTrackMultiplicityForCentrEvt;                //! Track multiplicity distribution (without any cut = include ITS stand-alone + global tracks)
        TH1F    *fHistTPCrefitTrackMultiplicityForCentrEvt;        //! Track multiplicity distribution for tracks with TPCrefit
        
        // - General histos (filled for any triggered event + (|z(prim. vtx)| < 10 cm ) )
        TH1F    *fHistCascadeMultiplicityForTrigEvtAndZprimVtx;   //! Cascade multiplicity distribution

        // - General histos (filled for events selected in this analysis (|z(prim. vtx)| < 10 cm + prim vtx = SPD or Tracking) )
        TH1F    *fHistCascadeMultiplicityForSelEvt;     //! Cascade multiplicity distribution
        TH1F    *fHistPosBestPrimaryVtxXForSelEvt;      //! (best) primary vertex position distribution in x 
        TH1F    *fHistPosBestPrimaryVtxYForSelEvt;      //! (best) primary vertex position distribution in y
        TH1F    *fHistPosBestPrimaryVtxZForSelEvt;      //! (best) primary vertex position distribution in z
        
        
        

        // - Characteristics for event with >1 cascade : Track Multiplicity, TPC clusters + Prim. vertex positions
        TH1F	*fHistTPCrefitTrackMultiplicityForCascadeEvt;   //! TPCrefit Track multiplicity distribution for event with >1 cascade candidate (NB: after quality sel. within the task)

        TH1F    *fHistPosV0TPCClusters;                 //! TPC clusters distribution for Positive V0 daughter track
        TH1F    *fHistNegV0TPCClusters;                 //! TPC clusters distribution for Negative V0 daughter track
        TH1F    *fHistBachTPCClusters;                  //! TPC clusters distribution for Bachelor             track

        TH1F    *fHistVtxStatus;                        //! Is there a tracking vertex in the cascade event ?

                // Vtx coming from the full tracking, for events containing at least a cascade
        TH1F    *fHistPosTrkgPrimaryVtxXForCascadeEvt;  //! primary vertex position distribution in x 
        TH1F    *fHistPosTrkgPrimaryVtxYForCascadeEvt;  //! primary vertex position distribution in y
        TH1F    *fHistPosTrkgPrimaryVtxZForCascadeEvt;  //! primary vertex position distribution in z
        TH1F    *fHistTrkgPrimaryVtxRadius;             //! primary vertex (3D) radius distribution 

                // Best primary Vtx available, for events containing at least a cascade
        TH1F    *fHistPosBestPrimaryVtxXForCascadeEvt;  //! (best) primary vertex position distribution in x 
        TH1F    *fHistPosBestPrimaryVtxYForCascadeEvt;  //! (best) primary vertex position distribution in y
        TH1F    *fHistPosBestPrimaryVtxZForCascadeEvt;  //! (best) primary vertex position distribution in z
        TH1F    *fHistBestPrimaryVtxRadius;             //! (best) primary vertex radius distribution 

                // Correlation Best Vtx / Full Tracking Vtx
        TH2F    *f2dHistTrkgPrimVtxVsBestPrimVtx;       //!  Radius of prim. Vtx from tracks Vs Radius of best Prim. Vtx


// PART 1 : Adavanced QA
// - Typical histos on the variables used for the selection of cascades
        TH1F    *fHistEffMassXi;      			//! reconstructed cascade effective mass
        TH1F    *fHistChi2Xi;         			//! chi2 value
        TH1F    *fHistDcaXiDaughters; 			//! dca between Xi's daughters
        TH1F    *fHistDcaBachToPrimVertex;		//! dca of the bachelor track to primary vertex
        TH1F    *fHistXiCosineOfPointingAngle;		//! cosine of Xi pointing angle in a cascade
        TH1F    *fHistXiRadius;                         //! (transverse) radius of the cascade vertex 
		
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
        TH1F	*fHistXiTransvMom;                      //! Xi transverse momentum, around the mass peak of Xi-/+ 
	TH1F	*fHistXiTotMom;                         //! Xi momentum norm, around the mass peak of Xi-/+
	
	TH1F	*fHistBachTransvMomXi;                  //! bachelor transverse momentum, for cand. around the mass peak of Xi-/+
	TH1F	*fHistBachTotMomXi;                     //! bachelor momentum norm, for cand. around the mass peak of Xi-/+

	TH1F	*fHistChargeXi;				//! Charge sign of the cascade candidate
	TH1F	*fHistV0toXiCosineOfPointingAngle;	//! Cos. of Pointing angle between the V0 mom and the Xi-V0 vtx line
  
	TH1F	*fHistRapXi;                            //! rapidity of Xi candidates, around the mass peak of Xi-/+
	TH1F	*fHistRapOmega;                         //! rapidity of Omega candidates, around the mass peak of Omega-/+
	TH1F	*fHistEtaXi;                            //! eta distrib. of all the cascade candidates, around the mass peak of Xi-/+
	TH1F	*fHistThetaXi;                          //! theta distrib. of all the cascade candidates, around the mass peak of Xi-/+
	TH1F	*fHistPhiXi;                            //! phi distrib. of all the cascade candidates, around the mass peak of Xi-/+
	
	TH2F	*f2dHistArmenteros;			//! alpha(casc. cand.) Vs PtArm(casc. cand.)
	
	TH2F	*f2dHistEffMassLambdaVsEffMassXiMinus;	//! Xi- Eff mass Vs V0 Eff mass, under Xi- hyp.
	TH2F	*f2dHistEffMassXiVsEffMassOmegaMinus;	//! Xi- Eff mass Vs Omega- Eff mass, for negative cascades
	TH2F	*f2dHistEffMassLambdaVsEffMassXiPlus;	//! Xi+ Eff mass Vs V0 Eff mass, under Xi+ hyp.
	TH2F	*f2dHistEffMassXiVsEffMassOmegaPlus;	//! Xi+ Eff mass Vs Omega+ Eff mass, for positive cascades
	
	TH2F	*f2dHistXiRadiusVsEffMassXiMinus;	//! transv. casc. decay radius Vs Xi- Eff mass, under Xi- hyp.
	TH2F	*f2dHistXiRadiusVsEffMassXiPlus;	//! transv. casc. decay radius Vs Xi+ Eff mass, under Xi+ hyp.
	TH2F	*f2dHistXiRadiusVsEffMassOmegaMinus;	//! transv. casc. decay radius Vs Omega- Eff mass, under Omega- hyp.
	TH2F	*f2dHistXiRadiusVsEffMassOmegaPlus;	//! transv. casc. decay radius Vs Omega+ Eff mass, under Omega+ hyp.
        
        TH2F    *f2dHistTPCdEdxOfCascDghters;           //! TPC Bethe-Bloch curve, populated with the cascade daughters
	
	
	// PART 2 : TH3F needed for pt spectrum and yield extraction
	// Without any PID
/*	TH3F	*f3dHistXiPtVsEffMassVsYXiMinus;        //! casc. transv. momemtum Vs Xi- Eff mass Vs Y
	TH3F	*f3dHistXiPtVsEffMassVsYXiPlus;         //! casc. transv. momemtum Vs Xi+ Eff mass Vs Y
	TH3F	*f3dHistXiPtVsEffMassVsYOmegaMinus;     //! casc. transv. momemtum Vs Omega- Eff mass Vs Y
	TH3F	*f3dHistXiPtVsEffMassVsYOmegaPlus;      //! casc. transv. momemtum Vs Omega+ Eff mass Vs Y
*/	
	// Compilation of all PID plots (3D = casc. transv. momemtum Vs Casc Eff mass Vs Y), stored into an AliCFContainer
	AliCFContainer  *fCFContCascadePIDXiMinus;      //! for Xi-   : Container to store any 3D histos with the different PID flavours
	AliCFContainer  *fCFContCascadePIDXiPlus;       //! for Xi+   : Container to store any 3D histos with the different PID flavours
	AliCFContainer  *fCFContCascadePIDOmegaMinus;   //! for Omega-: Container to store any 3D histos with the different PID flavours
	AliCFContainer  *fCFContCascadePIDOmegaPlus;    //! for Omega+: Container to store any 3D histos with the different PID flavours
	
	
	
	// PART 3 : Towards the optimisation of topological selections / systematics
	AliCFContainer  *fCFContCascadeCuts;            //! Container meant to store all the relevant distributions corresponding to the cut variables
	
	
	// PART 4 :  Azimuthal correlation study
/*	THnSparseF	*fHnSpAngularCorrXiMinus;	//! Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks Vs Eff Mass
	THnSparseF	*fHnSpAngularCorrXiPlus;	//! Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks Vs Eff Mass
	THnSparseF	*fHnSpAngularCorrOmegaMinus;	//! Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks Vs Eff Mass
	THnSparseF	*fHnSpAngularCorrOmegaPlus;	//! Delta Phi(Casc,any trck) Vs Delta Eta(Casc,any trck) Vs Casc Pt Vs Pt of the tracks Vs Eff Mass
*/
        TH1F *fV0Ampl;                                  //! histo to check the V0 amplitude distribution  

        TH2F    *fHistDcaXiDaughtersvsInvMass;          //! cut variables vs inv. mass
        TH2F    *fHistDcaBachToPrimVertexvsInvMass;     //! cut variables vs inv. mass
        TH2F    *fHistXiCosineOfPointingAnglevsInvMass; //! cut variables vs inv. mass
        TH2F    *fHistMassLambdaAsCascDghtervsInvMass;  //! cut variables vs inv. mass
        TH2F    *fHistDcaV0DaughtersXivsInvMass;        //! cut variables vs inv. mass
        TH2F    *fHistDcaV0ToPrimVertexXivsInvMass;     //! cut variables vs inv. mass



  AliAnalysisTaskCheckCascadePbPb(const AliAnalysisTaskCheckCascadePbPb&);            // not implemented
  AliAnalysisTaskCheckCascadePbPb& operator=(const AliAnalysisTaskCheckCascadePbPb&); // not implemented
  
  ClassDef(AliAnalysisTaskCheckCascadePbPb, 1);
};

#endif
