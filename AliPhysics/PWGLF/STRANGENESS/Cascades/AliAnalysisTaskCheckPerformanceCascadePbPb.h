#ifndef ALIANALYSISTASKCHECKPERFORMANCECASCADEPBPB_H
#define ALIANALYSISTASKCHECKPERFORMANCECASCADEPBPB_H

/*  See cxx source for full Copyright notice */

// //-----------------------------------------------------------------
// //        AliAnalysisTaskCheckPerformanceCascadePbPb class
// //            This task is for a performance study of cascade identification.
// //            It works with MC info and ESD and AOD tree 
// //            Origin   : A.Maire Jan2010, antonin.maire@ires.in2p3.fr
// //            Modified : M.Nicassio Feb2011, maria.nicassio@ba.infn.it
// //-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class TH3F; 
class AliESDEvent;
class AliESDtrackCuts;
class AliCFContainer;
class AliPIDResponse;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCheckPerformanceCascadePbPb : public AliAnalysisTaskSE {
 public:
 
  AliAnalysisTaskCheckPerformanceCascadePbPb();
  AliAnalysisTaskCheckPerformanceCascadePbPb(const char *name );
  virtual ~AliAnalysisTaskCheckPerformanceCascadePbPb();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetAnalysisType     (const char* analysisType    = "ESD") { fAnalysisType     = analysisType;}
  
  void SetRelaunchV0CascVertexers    (Bool_t rerunV0CascVertexers       = 0    ) { fkRerunV0CascVertexers         = rerunV0CascVertexers;      }
  void SetQualityCutZprimVtxPos      (Bool_t qualityCutZprimVtxPos      = kTRUE) { fkQualityCutZprimVtxPos        = qualityCutZprimVtxPos;     }
  void SetRejectEventPileUp          (Bool_t rejectPileUp               = kTRUE) { fkRejectEventPileUp            = rejectPileUp;              }
  void SetQualityCutNoTPConlyPrimVtx (Bool_t qualityCutNoTPConlyPrimVtx = kTRUE) { fkQualityCutNoTPConlyPrimVtx   = qualityCutNoTPConlyPrimVtx;}
  void SetQualityCutTPCrefit         (Bool_t qualityCutTPCrefit         = kTRUE) { fkQualityCutTPCrefit           = qualityCutTPCrefit;        }
  void SetQualityCutnTPCcls          (Bool_t qualityCutnTPCcls          = kTRUE) { fkQualityCutnTPCcls            = qualityCutnTPCcls;         }
  void SetQualityCutMinnTPCcls       (Int_t minnTPCcls                  = 70   ) { fMinnTPCcls                    = minnTPCcls;                }
  void SetExtraSelections            (Bool_t extraSelections            = 0    ) { fkExtraSelections              = extraSelections;           }
  void SetCentralityLowLim           (Float_t centrlowlim               = 0.   ) { fCentrLowLim                   = centrlowlim;               }
  void SetCentralityUpLim            (Float_t centruplim                = 100. ) { fCentrUpLim                    = centruplim;                }
  void SetCentralityEst              (TString centrest                  = "V0M") { fCentrEstimator                = centrest;                  }  
  void SetUseCleaning                (Bool_t   usecleaning              = kTRUE) { fkUseCleaning                  = usecleaning;               }
  void SetVertexRange                (Float_t vtxrange                  = 0.   ) { fVtxRange                      = vtxrange;                  }
  void SetApplyAccCut                (Bool_t acccut                     = kFALSE){ fApplyAccCut                   = acccut;                    }    
  void SetMinptCutOnDaughterTracks   (Float_t minptdaughtrks            = 0.   ) { fMinPtCutOnDaughterTracks      = minptdaughtrks;            }
  void SetEtaCutOnDaughterTracks     (Float_t etadaughtrks              = 0.   ) { fEtaCutOnDaughterTracks        = etadaughtrks;              }
  //Setters for the V0 and cascade Vertexer Parameters
  void SetV0VertexerMaxChisquare           (Double_t lParameter){ fV0VertexerSels[0] = lParameter; }
  void SetV0VertexerDCAFirstToPV           (Double_t lParameter){ fV0VertexerSels[1] = lParameter; }
  void SetV0VertexerDCASecondtoPV          (Double_t lParameter){ fV0VertexerSels[2] = lParameter; }
  void SetV0VertexerDCAV0Daughters         (Double_t lParameter){ fV0VertexerSels[3] = lParameter; }
  void SetV0VertexerCosinePA               (Double_t lParameter){ fV0VertexerSels[4] = lParameter; }
  void SetV0VertexerMinRadius              (Double_t lParameter){ fV0VertexerSels[5] = lParameter; }
  void SetV0VertexerMaxRadius              (Double_t lParameter){ fV0VertexerSels[6] = lParameter; }
  void SetCascVertexerMaxChisquare         (Double_t lParameter){ fCascadeVertexerSels[0] = lParameter; }
  void SetCascVertexerMinV0ImpactParameter (Double_t lParameter){ fCascadeVertexerSels[1] = lParameter; }
  void SetCascVertexerV0MassWindow         (Double_t lParameter){ fCascadeVertexerSels[2] = lParameter; }
  void SetCascVertexerDCABachToPV          (Double_t lParameter){ fCascadeVertexerSels[3] = lParameter; }
  void SetCascVertexerDCACascadeDaughters  (Double_t lParameter){ fCascadeVertexerSels[4] = lParameter; }
  void SetCascVertexerCascadeCosinePA      (Double_t lParameter){ fCascadeVertexerSels[5] = lParameter; }
  void SetCascVertexerCascadeMinRadius     (Double_t lParameter){ fCascadeVertexerSels[6] = lParameter; }
  void SetCascVertexerCascadeMaxRadius     (Double_t lParameter){ fCascadeVertexerSels[7] = lParameter; } 



 private:
        // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
        // your data member object is created on the worker nodes and streaming is not needed.
        // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14

        TString         fAnalysisType;          // "ESD" or "AOD" analysis type	
        AliESDtrackCuts *fESDtrackCuts;         // ESD track cuts used for primary track definition
      //TPaveText       *fPaveTextBookKeeping;  // TString to store all the relevant info necessary for book keeping (v0 cuts, cascade cuts, quality cuts, ...)
        AliPIDResponse *fPIDResponse;           //! PID response object        

        Bool_t          fkRerunV0CascVertexers;         // Boolean : kTRUE = relaunch both V0 + Cascade vertexers
        Bool_t          fkQualityCutZprimVtxPos;        // Boolean : kTRUE = cut on the prim.vtx  z-position
        Bool_t          fkRejectEventPileUp;            // Boolean : kTRUE = enable the rejection of events tagged as pile-up by SPD (AliESDEvent::IsPileupFromSPD)
        Bool_t          fkQualityCutNoTPConlyPrimVtx;   // Boolean : kTRUE = prim vtx should be SPD or Tracking vertex
        Bool_t          fkQualityCutTPCrefit;           // Boolean : kTRUE = ask for TPCrefit for the 3 daughter tracks
        Bool_t          fkQualityCutnTPCcls;            // Boolean : kTRUE = ask forat least n TPC clusters for each daughter track
        Int_t           fMinnTPCcls;                    // Minimum number of TPC clusters for each daughter track
        Bool_t          fkExtraSelections;              // Boolean : kTRUE = apply tighter selections, before starting the analysis
        Float_t         fCentrLowLim;                   // Lower limit for centrality percentile selection
        Float_t         fCentrUpLim;                    // Upper limit for centrality percentile selection
        TString         fCentrEstimator;                // String for the centrality estimator
        Bool_t          fkUseCleaning;                  // Boolean : kTRUE = uses all the cleaning criteria of centrality selections (vertex cut + outliers) otherwise only outliers
        Float_t         fVtxRange;                      // to select events with |zvtx|<fVtxRange cm
        Bool_t          fApplyAccCut;                   // flag to apply acceptance cuts to MC cascades        
        Float_t         fMinPtCutOnDaughterTracks;      // minimum pt to cut daughter tracks   
        Float_t         fEtaCutOnDaughterTracks;        // pseudorapidity cut on daughter tracks
 

        Double_t        fV0VertexerSels[7];                     // Array to store the 7 values for the different selections V0 related (if fkRerunV0CascVertexers)
        Double_t        fCascadeVertexerSels[8];                   // Array to store the 8 values for the different selections Casc. related (if fkRerunV0CascVertexers)
	
 	TList	*fListHistCascade;		//! List of Cascade histograms

                // - Histos
        TH2F    *fHistEvtsInCentralityBinsvsNtracks;    //! Events in centrality bins vs N ESDtracks
        TH1F    *fHistBestVtxX;                //! Vertex distribution
        TH1F    *fHistBestVtxY;                //! Vertex distribution
        TH1F    *fHistBestVtxZ;                //! Vertex distribution
        TH1F    *fHistnXiPlusPerEvTot;         //! Cascade multiplicity histogram
        TH1F    *fHistnXiMinusPerEvTot;        //! Cascade multiplicity histogram
        TH1F    *fHistnOmegaPlusPerEvTot;      //! Cascade multiplicity histogram
        TH1F    *fHistnOmegaMinusPerEvTot;     //! Cascade multiplicity histogram

        TH1F    *fHistnXiPlusPerEv;             //! Cascade multiplicity histograms          
        TH1F    *fHistnXiMinusPerEv;            //! Cascade multiplicity histograms
        TH1F    *fHistnOmegaPlusPerEv;          //! Cascade multiplicity histograms
        TH1F    *fHistnOmegaMinusPerEv;         //! Cascade multiplicity histograms 

        TH1F    *fHistnAssoXiMinus;             //! Cascade multiplicity histograms
        TH1F    *fHistnAssoXiPlus;              //! Cascade multiplicity histograms
        TH1F    *fHistnAssoOmegaMinus;          //! Cascade multiplicity histograms
        TH1F    *fHistnAssoOmegaPlus;           //! Cascade multiplicity histograms

	TH1F	*fHistMCTrackMultiplicity;	//! MC Track multiplicity (gen. primaries)
                // - Resolution of the multiplicity estimator
        TH2F    *f2dHistRecoMultVsMCMult;       //! resolution of the multiplicity estimator (based on primary tracks)
	
	
	// proton
	TH1F	*fHistEtaGenProton;   			//! MC Pseudo-rapidity of any generated p+ (no cuts in acceptance)
	TH1F	*fHistEtaGenAntiProton;   		//! MC Pseudo-rapidity of any generated p+ (no cuts in acceptance)

// Part 1 - Generated cascades
	
	//--------------
	// Xi-
	TH1F	*fHistEtaGenCascXiMinus;   		//! MC Pseudo-rapidity of any generated Xi- (no cuts in acceptance)
        TH3D    *f3dHistGenPtVsGenYvsCentXiMinusNat;
        TH3D    *f3dHistGenPtVsGenYvsNtracksXiMinusNat;
	TH3D    *f3dHistGenPtVsGenYvsCentXiMinusInj;
        TH3D    *f3dHistGenPtVsGenYvsNtracksXiMinusInj;
        TH3D    *f3dHistGenPtVsGenctauvsCentXiMinusNat;
        TH3D    *f3dHistGenPtVsGenctauvsCentXiMinusInj;	

        TH1F    *fHistThetaGenCascXiMinusNat;           //! MC Theta angle of the generated Xi-
        TH1F    *fHistThetaGenCascXiMinusInj;           //! MC Theta angle of the injected Xi-
	// - Histos planned for Xi- emitted within the acceptance (cuts in theta + pt of daughters)
	// 	= findable cascades
	TH2D	*f2dHistGenPtVsGenYFdblXiMinus;		//! MC Pt Vs MC y of the findable Xi-
	
	TH1F	*fHistThetaLambdaXiMinus;		//! MC Theta angle of the Lambda daughter of the generated Xi-
	TH1F	*fHistThetaBachXiMinus;			//! MC Theta angle of the Bachelor (pi-)
	
	TH1F	*fHistThetaMesDghterXiMinus;		//! MC Theta angle of the mesonic  V0 daughter in the generated cascade, pi-
	TH1F	*fHistThetaBarDghterXiMinus;		//! MC Theta angle of the baryonic V0 daughter in the generated cascade, p+
	
	TH1F	*fHistPtBachXiMinus;   			//! MC Pt of the Bachelor (pi-)                         (Control Plot)
	TH1F	*fHistPtMesDghterXiMinus;		//! MC Pt of the meson daughter of the 'Lambda0', pi-   (Control Plot)
	TH1F	*fHistPtBarDghterXiMinus;		//! MC Pt of the baryon daughter of the 'Lambda0', p+   (Control Plot)
	
        TH1F    *fHistPtRecBachXiMinus;                 //! Rec Pt of the Bachelor (for Xi-)                    (Control Plot)
        TH1F    *fHistPtRecMesDghterXiMinus;            //! Rec Pt of the meson daughter of the 'Lambda0', pi-   (Control Plot)
        TH1F    *fHistPtRecBarDghterXiMinus;            //! Rec Pt of the baryon daughter of the 'Lambda0', p+   (Control Plot)
	
	
	//--------------
	// Xi+
	TH1F	*fHistEtaGenCascXiPlus;   		//! MC Pseudo-rapidity of any generated Xi+ (no cuts in acceptance)
        TH3D    *f3dHistGenPtVsGenYvsCentXiPlusNat;
        TH3D    *f3dHistGenPtVsGenYvsNtracksXiPlusNat;
        TH3D    *f3dHistGenPtVsGenYvsCentXiPlusInj;
        TH3D    *f3dHistGenPtVsGenYvsNtracksXiPlusInj;
        TH3D    *f3dHistGenPtVsGenctauvsCentXiPlusNat;
        TH3D    *f3dHistGenPtVsGenctauvsCentXiPlusInj;
	
        TH1F    *fHistThetaGenCascXiPlusNat;            //! MC Theta angle of the generated Xi+
        TH1F    *fHistThetaGenCascXiPlusInj;            //! MC Theta angle of the injected Xi+
	// - Histos planned for Xi+ emitted within the acceptance (cuts in theta + pt of daughters)
	TH2D	*f2dHistGenPtVsGenYFdblXiPlus;		//! MC Pt Vs MC y of the findable Xi+
	
	TH1F	*fHistThetaLambdaXiPlus;		//! MC Theta angle of the anti-Lambda daughter of the generated Xi+
	TH1F	*fHistThetaBachXiPlus;			//! MC Theta angle of the Bachelor (pi+)
	
	TH1F	*fHistThetaMesDghterXiPlus;		//! MC Theta angle of the mesonic  V0 daughter in the generated cascade, pi+
	TH1F	*fHistThetaBarDghterXiPlus;		//! MC Theta angle of the baryonic V0 daughter in the generated cascade, p-
	
	TH1F	*fHistPtBachXiPlus;   			//! MC Pt of the Bachelor (pi+)                       (Control Plot)
	TH1F	*fHistPtMesDghterXiPlus;		//! MC Pt of the meson daughter of the 'Lambda0', pi+ (Control Plot)
	TH1F	*fHistPtBarDghterXiPlus;		//! MC Pt of the baryon daughter of the 'Lambda0', p- (Control Plot)
	
	
	
	//--------------
	// Omega-
	TH1F	*fHistEtaGenCascOmegaMinus;   		//! MC Pseudo-rapidity of any generated Omega- (no cuts in acceptance)
        TH3D    *f3dHistGenPtVsGenYvsCentOmegaMinusNat;
        TH3D    *f3dHistGenPtVsGenYvsNtracksOmegaMinusNat;
        TH3D    *f3dHistGenPtVsGenYvsCentOmegaMinusInj;
        TH3D    *f3dHistGenPtVsGenYvsNtracksOmegaMinusInj;
        TH3D    *f3dHistGenPtVsGenctauvsCentOmegaMinusNat;
        TH3D    *f3dHistGenPtVsGenctauvsCentOmegaMinusInj;

        TH1F    *fHistThetaGenCascOmegaMinusNat;        //! MC Theta angle of the generated Omega-
        TH1F    *fHistThetaGenCascOmegaMinusInj;        //! MC Theta angle of the injected Omega-
	// - Histos planned for Omega- emitted within the acceptance (cuts in theta + pt of daughters)
	TH2D	*f2dHistGenPtVsGenYFdblOmegaMinus;	//! MC Pt Vs MC y of the findable Omega-
	
	TH1F	*fHistThetaLambdaOmegaMinus;		//! MC Theta angle of the Lambda daughter of the generated Omega-
	TH1F	*fHistThetaBachOmegaMinus;		//! MC Theta angle of the Bachelor (K-)
	
	TH1F	*fHistThetaMesDghterOmegaMinus;		//! MC Theta angle of the mesonic  V0 daughter in the generated cascade, pi-
	TH1F	*fHistThetaBarDghterOmegaMinus;		//! MC Theta angle of the baryonic V0 daughter in the generated cascade, p+
	
	TH1F	*fHistPtBachOmegaMinus;   		//! MC Pt of the Bachelor (K-)                   (Control Plot)
	TH1F	*fHistPtMesDghterOmegaMinus;		//! MC Pt of the meson daughter of the 'Lambda0', pi- (Control Plot)
	TH1F	*fHistPtBarDghterOmegaMinus;		//! MC Pt of the baryon daughter of the 'Lambda0', p+ (Control Plot)
	
	
	
	//--------------
	// Omega+
	TH1F	*fHistEtaGenCascOmegaPlus;   		//! MC Pseudo-rapidity of any generated Omega+ (no cuts in acceptance)
        TH3D    *f3dHistGenPtVsGenYvsCentOmegaPlusNat;
        TH3D    *f3dHistGenPtVsGenYvsNtracksOmegaPlusNat;
        TH3D    *f3dHistGenPtVsGenYvsCentOmegaPlusInj;
        TH3D    *f3dHistGenPtVsGenYvsNtracksOmegaPlusInj;
        TH3D    *f3dHistGenPtVsGenctauvsCentOmegaPlusNat;
        TH3D    *f3dHistGenPtVsGenctauvsCentOmegaPlusInj;
	
        TH1F    *fHistThetaGenCascOmegaPlusNat;         //! MC Theta angle of the generated Omega+
        TH1F    *fHistThetaGenCascOmegaPlusInj;         //! MC Theta angle of the injected Omega+
	// - Histos planned for Omega+ emitted within the acceptance (cuts in theta + pt of daughters)
	TH2D	*f2dHistGenPtVsGenYFdblOmegaPlus;	//! MC Pt Vs MC y of the findable Omega+
	
	TH1F	*fHistThetaLambdaOmegaPlus;		//! MC Theta angle of the anti-Lambda daughter of the generated Omega+
	TH1F	*fHistThetaBachOmegaPlus;		//! MC Theta angle of the Bachelor (K+)
	
	TH1F	*fHistThetaMesDghterOmegaPlus;		//! MC Theta angle of the mesonic  V0 daughter in the generated cascade, pi+
	TH1F	*fHistThetaBarDghterOmegaPlus;		//! MC Theta angle of the baryonic V0 daughter in the generated cascade, p-
	
	TH1F	*fHistPtBachOmegaPlus;   		//! MC Pt of the Bachelor (K+) (Control Plot)
	TH1F	*fHistPtMesDghterOmegaPlus;		//! MC Pt of the meson daughter of the 'Lambda0', pi+ (Control Plot)
	TH1F	*fHistPtBarDghterOmegaPlus;		//! MC Pt of the baryon daughter of the 'Lambda0', p- (Control Plot)
	
	
	
// Part 2 - Any reconstructed cascades + reconstructed cascades associated with MC
	// 2.1 - Effective mass and PID
	// - Effective mass histos for all the cascade candidates
	TH1F	*fHistMassXiMinus;			//! reconstructed cascade effective mass, under Xi- hyp.
	TH1F	*fHistMassXiPlus;			//! reconstructed cascade effective mass, under Xi+ hyp.
	TH1F	*fHistMassOmegaMinus;			//! reconstructed cascade effective mass, under Omega- hyp.
	TH1F	*fHistMassOmegaPlus;			//! reconstructed cascade effective mass, under Omega+ hyp.
	
	// - Effective mass histos with reconstruction combined PID
	TH1F	*fHistMassWithCombPIDXiMinus;		//! reconstructed Xi- effective mass, with bach. comb PID
	TH1F	*fHistMassWithCombPIDXiPlus;		//! reconstructed Xi+ effective mass, with bach. comb PID
	TH1F	*fHistMassWithCombPIDOmegaMinus;	//! reconstructed Omega- effective mass, with bach. comb PID
	TH1F	*fHistMassWithCombPIDOmegaPlus;		//! reconstructed Omega+ effective mass, with bach. comb PID
	
	// - PID Probability versus MC Pt(bachelor track)
	TH2F	*f2dHistPIDprobaKaonVsMCPtBach;		//! Comb. PID probability for the bach. to be a Kaon Vs MC pt(bach)
	TH2F	*f2dHistPIDprobaPionVsMCPtBach;		//! Comb. PID probability for the bach. to be a Pion Vs MC pt(bach)	

	// - Effective mass histos with perfect MC PID
	TH1F	*fHistMassWithMcPIDXiMinus;		//! reconstructed Xi- effective mass, with MC bach. PID
	TH1F	*fHistMassWithMcPIDXiPlus;		//! reconstructed Xi+ effective mass, with MC bach. PID
	TH1F	*fHistMassWithMcPIDOmegaMinus;		//! reconstructed Omega- effective mass, with MC bach. PID
	TH1F	*fHistMassWithMcPIDOmegaPlus;		//! reconstructed Omega+ effective mass, with MC bach. PID

	
	// 2.2 - Associated candidates
	// - Effective mass histos for the cascade candidates associated with MC, without PID info
	TH1F	*fHistAsMCMassXiMinus;			//! reconstr. cascade effective mass, under Xi- hyp. for Associated cand.
	TH1F	*fHistAsMCMassXiPlus;			//! reconstr. cascade effective mass, under Xi+ hyp. for Associated cand.
	TH1F	*fHistAsMCMassOmegaMinus;		//! reconstr. cascade effective mass, under Omega- hyp. for Associated cand.
	TH1F	*fHistAsMCMassOmegaPlus;		//! reconstr. cascade effective mass, under Omega+ hyp. for Associated cand.
	
	// -  Generated Pt Vs generated Y of the cascade candidates associated with MC 
	//     + having the proper maximum proba of combined PID for the bachelor
	TH2F	*f2dHistAsMCandCombPIDGenPtVsGenYXiMinus;	//! Pt(gen) Vs Y(gen) from the MC Xi- associated with Reco cand + with PID info
	TH2F	*f2dHistAsMCandCombPIDGenPtVsGenYXiPlus;	//! Pt(gen) Vs Y(gen) from the MC Xi+ associated with Reco cand + with PID info
	TH2F	*f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus;	//! Pt(gen) Vs Y(gen) from the MC Omega- associated with Reco cand + with PID info
	TH2F	*f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus;	//! Pt(gen) Vs Y(gen) from the MC Omega+ associated with Reco cand + with PID info
	
	// - Generated Pt Vs generated Y, for the cascade candidates associated with MC, without PID info
	TH2F	*f2dHistAsMCGenPtVsGenYXiMinus;		//! gen. Pt Vs gen. Rap. from the MC Xi- associated with a reconstr. cascade
	TH2F	*f2dHistAsMCGenPtVsGenYXiPlus;		//! gen. Pt Vs gen. Rap. from the MC Xi+ associated with a reconstr. cascade
	TH2F	*f2dHistAsMCGenPtVsGenYOmegaMinus;	//! gen. Pt Vs gen. Rap. from the MC Omega- associated with a reconstr. cascade
	TH2F	*f2dHistAsMCGenPtVsGenYOmegaPlus;	//! gen. Pt Vs gen. Rap. from the MC Omega+ associated with a reconstr. cascade
	
	// - Generated Eta of the the cascade candidates associated with MC, without PID info
	TH1F	*fHistAsMCGenEtaXiMinus;		//! generated Eta from the MC Xi- associated with a reconstr. cascade
	TH1F	*fHistAsMCGenEtaXiPlus;			//! generated Eta from the MC Xi+ associated with a reconstr. cascade
	TH1F	*fHistAsMCGenEtaOmegaMinus;		//! generated Eta from the MC Omega- associated with a reconstr. cascade
	TH1F	*fHistAsMCGenEtaOmegaPlus;		//! generated Eta from the MC Omega+ associated with a reconstr. cascade
	
	// - Resolution in Pt as function of generated Pt
	TH2F	*f2dHistAsMCResPtXiMinus;		//! resolution in Pt as function of gen. Pt, for Xi-
	TH2F	*f2dHistAsMCResPtXiPlus;		//! resolution in Pt as function of gen. Pt, for Xi-
	TH2F	*f2dHistAsMCResPtOmegaMinus;		//! resolution in Pt as function of gen. Pt, for Omega-
	TH2F	*f2dHistAsMCResPtOmegaPlus;		//! resolution in Pt as function of gen. Pt, for Omega+
	
	// - Resolution in R(2D) as function of generated R
	TH2F	*f2dHistAsMCResRXiMinus;		//! resolution in transv. R = f(transv. gen. R), for Xi-
	TH2F	*f2dHistAsMCResRXiPlus;			//! resolution in transv. R = f(transv. gen. R), for Xi+
	TH2F	*f2dHistAsMCResROmegaMinus;		//! resolution in transv. R = f(transv. gen. R), for Omega-
	TH2F	*f2dHistAsMCResROmegaPlus;		//! resolution in transv. R = f(transv. gen. R), for Omega+
        
        // - Resolution in phi as function of generated Pt
        TH2F    *f2dHistAsMCResPhiXiMinus;              //! resolution in azimuth Phi = f(gen. Pt), for Xi-
        TH2F    *f2dHistAsMCResPhiXiPlus;               //! resolution in azimuth Phi = f(gen. Pt), for Xi+
        TH2F    *f2dHistAsMCResPhiOmegaMinus;           //! resolution in azimuth Phi = f(gen. Pt), for Omega-
        TH2F    *f2dHistAsMCResPhiOmegaPlus;            //! resolution in azimuth Phi = f(gen. Pt), for Omega+

        TH2F    *f2dHistAsMCptProtonMCptXiMinus;        //! MC pt proton vs Mc pt Xi-
        TH2F    *f2dHistAsMCptAntiprotonMCptXiPlus;     //! MC pt antiproton vs Mc pt Xi+
        TH2F    *f2dHistAsMCptProtonMCptOmegaMinus;     //! MC pt proton vs Mc pt Omega-
        TH2F    *f2dHistAsMCptAntiprotonMCptOmegaPlus;  //! MC pt antiproton vs Mc pt Omega+

        TH1F    *fHistV0toXiCosineOfPointingAngle;      //! To check new V0 CosPA cut
        TH2F    *fHistV0CosineOfPointingAnglevsPtXi;    //! To check new V0 CosPA cut 
        TH2F    *fHistV0CosineOfPointingAnglevsPtOmega; //! To check new V0 CosPA cut 
        
        // - Compilation of all PID plots (3D = casc. transv. momemtum Vs Casc Eff mass Vs Y), stored into an AliCFContainer
	AliCFContainer  *fCFContCascadePIDAsXiMinus;      //! for Xi-   : Container to store any 3D histos with the different PID flavours
	AliCFContainer  *fCFContCascadePIDAsXiPlus;       //! for Xi+   : Container to store any 3D histos with the different PID flavours
	AliCFContainer  *fCFContCascadePIDAsOmegaMinus;   //! for Omega-: Container to store any 3D histos with the different PID flavours
	AliCFContainer  *fCFContCascadePIDAsOmegaPlus;    //! for Omega+: Container to store any 3D histos with the different PID flavours
	
	// - Towards the optimisation of topological selections/ systematics (on associated candidates)
	AliCFContainer  *fCFContAsCascadeCuts;            //! Container meant to store all the relevant distributions corresponding to the cut variables

        TH1F *fV0Ampl;                                    //! Histo to check the V0 amplitude distribution (centrality estimator)  

         // Control plots for reco pseudorapidity of daughter tracks (Xi- asso only)

        TH1F    *fHistEtaBachXiM;                       //! bachelor pseudorapidity
        TH1F    *fHistEtaPosXiM;                        //! positive daughter pseudorapidity
        TH1F    *fHistEtaNegXiM;                        //! negative daughter pseudorapidity


  AliAnalysisTaskCheckPerformanceCascadePbPb(const AliAnalysisTaskCheckPerformanceCascadePbPb&);            // not implemented
  AliAnalysisTaskCheckPerformanceCascadePbPb& operator=(const AliAnalysisTaskCheckPerformanceCascadePbPb&); // not implemented
  
  ClassDef(AliAnalysisTaskCheckPerformanceCascadePbPb, 8);
};

#endif
