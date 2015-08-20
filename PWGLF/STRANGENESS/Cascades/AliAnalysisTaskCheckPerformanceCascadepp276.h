#ifndef ALIANALYSISTASKCHECKPERFORMANCECASCADEPP276_H
#define ALIANALYSISTASKCHECKPERFORMANCECASCADEPP276_H

/*  See cxx source for full Copyright notice */

// //-----------------------------------------------------------------
// //        AliAnalysisTaskCheckPerformanceCascadePbPb class
// //            This task is for a performance study of cascade identification.
// //            It works with MC info and ESD and AOD tree 
// //            Origin   : A.Maire    Jan2010, antonin.maire@ires.in2p3.fr
// //            Modified : M.Nicassio Feb2011, maria.nicassio@ba.infn.it
// //            Modified : D.Colella  Feb2012, domenico.colella@ba.infn.it
// //-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class TH3F; 

class AliESDEvent;
class AliESDtrackCuts;
class AliPhysicsSelection;
class AliCFContainer;
class AliPIDResponse;
class AliAnalysisUtils;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCheckPerformanceCascadepp276 : public AliAnalysisTaskSE {
 public:
 
  AliAnalysisTaskCheckPerformanceCascadepp276();
  AliAnalysisTaskCheckPerformanceCascadepp276(const char *name );
  virtual ~AliAnalysisTaskCheckPerformanceCascadepp276();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetAnalysisType               (const char* analysisType           = "ESD") { fAnalysisType                 = analysisType;               }
  void SetCollidingSystem            (const char* collidingSystem        = "pp" ) { fCollidingSystem              = collidingSystem;            }  
  void SetRelaunchV0CascVertexers    (Bool_t  rerunV0CascVertexers       = 0    ) { fkRerunV0CascVertexers        = rerunV0CascVertexers;       }
  void SetSDDSelection               (Bool_t  sddOnSelection             = kTRUE) { fkSDDselectionOn              = sddOnSelection;             }
  void SetQualityCutZprimVtxPos      (Bool_t  qualityCutZprimVtxPos      = kTRUE) { fkQualityCutZprimVtxPos       = qualityCutZprimVtxPos;      }
  void SetRejectEventPileUp          (Bool_t  rejectPileUp               = kTRUE) { fkRejectEventPileUp           = rejectPileUp;               }
  void SetQualityCutNoTPConlyPrimVtx (Bool_t  qualityCutNoTPConlyPrimVtx = kTRUE) { fkQualityCutNoTPConlyPrimVtx  = qualityCutNoTPConlyPrimVtx; }
  void SetQualityCutTPCrefit         (Bool_t  qualityCutTPCrefit         = kTRUE) { fkQualityCutTPCrefit          = qualityCutTPCrefit;         }
  void SetQualityCutnTPCcls          (Bool_t  qualityCutnTPCcls          = kTRUE) { fkQualityCutnTPCcls           = qualityCutnTPCcls;          }
  void SetWithSDDOn                  (Bool_t  withsddOn                  = kTRUE) { fwithSDD                      = withsddOn;                  }
  void SetQualityCutMinnTPCcls       (Int_t   minnTPCcls                 = 70   ) { fMinnTPCcls                   = minnTPCcls;                 }
  void SetExtraSelections            (Bool_t  extraSelections            = 0    ) { fkExtraSelections             = extraSelections;            }
  void SetVertexRange                (Float_t vtxrange                   = 0.   ) { fVtxRange                     = vtxrange;                   }
  void SetVertexRangeMin             (Float_t vtxrangemin                = 0.   ) { fVtxRangeMin                  = vtxrangemin;                }
  void SetApplyAccCut                (Bool_t  acccut                     = kFALSE){ fApplyAccCut                  = acccut;                     }    
  void SetMinptCutOnDaughterTracks   (Float_t minptdaughtrks             = 0.   ) { fMinPtCutOnDaughterTracks     = minptdaughtrks;             }
  void SetEtaCutOnDaughterTracks     (Float_t etadaughtrks               = 0.   ) { fEtaCutOnDaughterTracks       = etadaughtrks;               }

 private:
        // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
        // your data member object is created on the worker nodes and streaming is not needed.
        // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14

        TString         fAnalysisType;          // "ESD" or "AOD" analysis type	
        AliESDtrackCuts *fESDtrackCuts;         // ESD track cuts used for primary track definition
        AliAnalysisUtils *fUtils;
        TString         fCollidingSystem;       // "pPb" or "pp" colliding system
       // AliESDtrackCuts *fESDtrackCuts;         // ESD track cuts used for primary track definition
        AliPIDResponse *fPIDResponse;           //! PID response object        
       // AliAnalysisUtils *fUtils;

        Bool_t          fkRerunV0CascVertexers;         // Boolean : kTRUE = relaunch both V0 + Cascade vertexers
        Bool_t          fkSDDselectionOn;               // Boolean : kTRUE = enable the selection based on the SDD status
        Bool_t          fkQualityCutZprimVtxPos;        // Boolean : kTRUE = cut on the prim.vtx  z-position
        Bool_t          fkRejectEventPileUp;            // Boolean : kTRUE = enable the rejection of events tagged as pile-up by SPD (AliESDEvent::IsPileupFromSPD)
        Bool_t          fkQualityCutNoTPConlyPrimVtx;   // Boolean : kTRUE = prim vtx should be SPD or Tracking vertex
        Bool_t          fkQualityCutTPCrefit;           // Boolean : kTRUE = ask for TPCrefit for the 3 daughter tracks
        Bool_t          fkQualityCutnTPCcls;            // Boolean : kTRUE = ask for n TPC clusters for each daughter track
        Bool_t          fwithSDD;                       // Boolean : kTRUE = select events with SDD reco
        Int_t           fMinnTPCcls;                    // Boolean : set the value for the minimum number of TPC clusters
        Bool_t          fkExtraSelections;              // Boolean : kTRUE = apply tighter selections, before starting the analysis
        Float_t         fVtxRange;                      // to select events with |zvtx|<fVtxRange cm
        Float_t         fVtxRangeMin;                   // to select events with |zvtx|>fVtxRangeMin cm
        Bool_t          fApplyAccCut;                   // flag to apply acceptance cuts to MC cascades       
        Float_t         fMinPtCutOnDaughterTracks;      // minimum pt cut on daughter tracks
        Float_t         fEtaCutOnDaughterTracks;        // pseudorapidity cut on daughter tracks 
        
        Double_t        fV0Sels[7];                     // Array to store the 7 values for the different selections V0 related (if fkRerunV0CascVertexers)
        Double_t        fCascSels[8];                   // Array to store the 8 values for the different selections Casc. related (if fkRerunV0CascVertexers)
	
 	TList	*fListHistCascade;		        //! List of Cascade histograms
         // - General Plots
         // Cascade multiplicity plots
         TH1F *fHistCascadeMultiplicityBeforeAnySel;
         TH1F *fHistCascadeMultiplicityAfterSDDSel;
         TH1F *fHistCascadeMultiplicityAfterPhysicsSel;
         TH1F *fHistCascadeMultiplicityForSelEvtNoTPCOnly;
         TH1F *fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup;
         TH1F *fHistCascadeMultiplicityAfterVertexCutSel;
         TH1F *fHistnXiPlusPerEvTot;                  // After any event selections, in all the eta and pt range
         TH1F *fHistnXiMinusPerEvTot;                 // After any event selections, in all the eta and pt range
         TH1F *fHistnOmegaPlusPerEvTot;               // After any event selections, in all the eta and pt range
         TH1F *fHistnOmegaMinusPerEvTot;              // After any event selections, in all the eta and pt range
         TH1F *fHistnXiPlusPerEv;                     // After any event selections, in the detector acceptance and over a pt minimum
         TH1F *fHistnXiMinusPerEv;                    // After any event selections, in the detector acceptance and over a pt minimum
         TH1F *fHistnOmegaPlusPerEv;                  // After any event selections, in the detector acceptance and over a pt minimum
         TH1F *fHistnOmegaMinusPerEv;                 // After any event selections, in the detector acceptance and over a pt minimum
         TH1F *fHistnAssoXiMinus;                     // For the Reconstructed-Associated cascades 
         TH1F *fHistnAssoXiPlus;                      // For the Reconstructed-Associated cascades 
         TH1F *fHistnAssoOmegaMinus;                  // For the Reconstructed-Associated cascades 
         TH1F *fHistnAssoOmegaPlus;                   // For the Reconstructed-Associated cascades 
         // Tracks multiplicity plots
         TH1F *fHistTrackMultiplicityBeforeAnySel;
         TH1F *fHistTrackMultiplicityAfterSDDSel;
         TH1F *fHistTrackMultiplicityAfterPhysicsSel;
         TH1F *fHistTrackMultiplicityForSelEvtNoTPCOnly;
         TH1F *fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup;
         TH1F *fHistTrackMultiplicityAfterVertexCutSel;
         // Vertex position plots (BestVertex)
         TH1F *fHistPVx;                              // After any selections but before |Z| < 10 cm
         TH1F *fHistPVy;                              // After any selections but before |Z| < 10 cm
         TH1F *fHistPVz;                              // After any selections but before |Z| < 10 cm
         TH1F *fHistPVxAnalysis;                      // After any event selections 
         TH1F *fHistPVyAnalysis;                      // After any event selections
         TH1F *fHistPVzAnalysis;                      // After any event selections
         // - Plots before Physics Selection
         TH3D *f3dHistGenPtVsGenYvsNtracksXiMinus;    // After the SDD event selection (For efficinecy calculation)
         TH3D *f3dHistGenPtVsGenctauvsYXiMinus;       // After the SDD event selection (For efficinecy calculation)
         TH3D *f3dHistGenPtVsGenYvsNtracksXiPlus;     // After the SDD event selection (For efficinecy calculation)
         TH3D *f3dHistGenPtVsGenctauvsYXiPlus;        // After the SDD event selection (For efficinecy calculation)
         TH3D *f3dHistGenPtVsGenYvsNtracksOmegaMinus; // After the SDD event selection (For efficinecy calculation)
         TH3D *f3dHistGenPtVsGenctauvsYOmegaMinus;    // After the SDD event selection (For efficinecy calculation)
         TH3D *f3dHistGenPtVsGenYvsNtracksOmegaPlus;  // After the SDD event selection (For efficinecy calculation)
         TH3D *f3dHistGenPtVsGenctauvsYOmegaPlus;     // After the SDD event selection (For efficinecy calculation)

         // - Generated cascade plots
         // After all the event selections 
         //Xi-
         TH1F *fHistEtaGenCascXiMinus;                // In all the eta and pt range (as they are generated)
         TH1F *fHistThetaGenCascXiMinus;              // In all the eta and pt range (as they are generated)
         TH3D *f3dHistGenPtVsGenYvsNtracksXiMinusPhysEff;    // 
         TH3D *f3dHistGenPtVsGenctauvsYXiMinusPhysEff;       // 
         TH2D *f2dHistGenPtVsGenYFdblXiMinus;         // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaLambdaXiMinus;               // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaBachXiMinus;                 // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaMesDghterXiMinus;            // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaBarDghterXiMinus;            // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistPtBachXiMinus;                    // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistPtMesDghterXiMinus;               // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistPtBarDghterXiMinus;               // In the detector acceptance and over a pt minimum (Findable particle)
         //Xi+
         TH1F *fHistEtaGenCascXiPlus;                 // In all the eta and pt range (as they are generated)
         TH1F *fHistThetaGenCascXiPlus;               // In all the eta and pt range (as they are generated)
         TH3D *f3dHistGenPtVsGenYvsNtracksXiPlusPhysEff;    // 
         TH3D *f3dHistGenPtVsGenctauvsYXiPlusPhysEff;       // 
         TH2D *f2dHistGenPtVsGenYFdblXiPlus;          // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaLambdaXiPlus;                // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaBachXiPlus;                  // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaMesDghterXiPlus;             // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaBarDghterXiPlus;             // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistPtBachXiPlus;                     // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistPtMesDghterXiPlus;                // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistPtBarDghterXiPlus;                // In the detector acceptance and over a pt minimum (Findable particle)
         //Omega-
         TH1F *fHistEtaGenCascOmegaMinus;             // In all the eta and pt range (as they are generated)
         TH1F *fHistThetaGenCascOmegaMinus;           // In all the eta and pt range (as they are generated)
         TH3D *f3dHistGenPtVsGenYvsNtracksOmegaMinusPhysEff;    // 
         TH3D *f3dHistGenPtVsGenctauvsYOmegaMinusPhysEff;       //
         TH2D *f2dHistGenPtVsGenYFdblOmegaMinus;      // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaLambdaOmegaMinus;            // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaBachOmegaMinus;              // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaMesDghterOmegaMinus;         // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaBarDghterOmegaMinus;         // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistPtBachOmegaMinus;                 // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistPtMesDghterOmegaMinus;            // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistPtBarDghterOmegaMinus;            // In the detector acceptance and over a pt minimum (Findable particle)
         //Omega+      
         TH1F *fHistEtaGenCascOmegaPlus;              // In all the eta and pt range (as they are generated)
         TH1F *fHistThetaGenCascOmegaPlus;            // In all the eta and pt range (as they are generated)
         TH3D *f3dHistGenPtVsGenYvsNtracksOmegaPlusPhysEff;    // 
         TH3D *f3dHistGenPtVsGenctauvsYOmegaPlusPhysEff;       //
         TH2D *f2dHistGenPtVsGenYFdblOmegaPlus;       // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaLambdaOmegaPlus;             // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaBachOmegaPlus;               // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaMesDghterOmegaPlus;          // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistThetaBarDghterOmegaPlus;          // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistPtBachOmegaPlus;                  // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistPtMesDghterOmegaPlus;             // In the detector acceptance and over a pt minimum (Findable particle)
         TH1F *fHistPtBarDghterOmegaPlus;             // In the detector acceptance and over a pt minimum (Findable particle)

         // - Associated to MC cascade plots
         TH1F *fHistMassXiMinus;                      // For the Reconstructed-Associated cascades
         TH1F *fHistMassXiPlus;                       // For the Reconstructed-Associated cascades
         TH1F *fHistMassOmegaMinus;                   // For the Reconstructed-Associated cascades
         TH1F *fHistMassOmegaPlus;                    // For the Reconstructed-Associated cascades
         // Effective mass histos with combined PID
         TH1F *fHistMassWithCombPIDXiMinus;
         TH1F *fHistMassWithCombPIDXiPlus;
         TH1F *fHistMassWithCombPIDOmegaMinus;
         TH1F *fHistMassWithCombPIDOmegaPlus;
         // PID Probability versus MC Pt(bachelor track)
         TH2F *f2dHistPIDprobaKaonVsMCPtBach; 
         TH2F *f2dHistPIDprobaPionVsMCPtBach;
         // Effective mass histos with perfect MC PID on the bachelor
         TH1F *fHistMassWithMcPIDXiMinus; 
         TH1F *fHistMassWithMcPIDXiPlus;
         TH1F *fHistMassWithMcPIDOmegaMinus; 
         TH1F *fHistMassWithMcPIDOmegaPlus;
         // Effective mass histos for the cascade candidates associated with MC
         TH1F *fHistAsMCMassXiMinus;
         TH1F *fHistAsMCMassXiPlus;
         TH1F *fHistAsMCMassOmegaMinus;
         TH1F *fHistAsMCMassOmegaPlus;
         // Generated Pt Vs generated y, for the cascade candidates associated with MC + Info Comb. PID
         TH2F *f2dHistAsMCandCombPIDGenPtVsGenYXiMinus;
         TH2F *f2dHistAsMCandCombPIDGenPtVsGenYXiPlus;
         TH2F *f2dHistAsMCandCombPIDGenPtVsGenYOmegaMinus;
         TH2F *f2dHistAsMCandCombPIDGenPtVsGenYOmegaPlus;
         // Generated Pt Vs generated y, for the cascade candidates associated with MC
         TH2F *f2dHistAsMCGenPtVsGenYXiMinus;
         TH2F *f2dHistAsMCGenPtVsGenYXiPlus;
         TH2F *f2dHistAsMCGenPtVsGenYOmegaMinus;
         TH2F *f2dHistAsMCGenPtVsGenYOmegaPlus;
         // Generated Eta of the the cascade candidates associated with MC
         TH1F *fHistAsMCGenEtaXiMinus;
         TH1F *fHistAsMCGenEtaXiPlus;
         TH1F *fHistAsMCGenEtaOmegaMinus;
         TH1F *fHistAsMCGenEtaOmegaPlus;
         // Resolution in Pt as function of generated Pt
         TH2F *f2dHistAsMCResPtXiMinus;
         TH2F *f2dHistAsMCResPtXiPlus;
         TH2F *f2dHistAsMCResPtOmegaMinus;
         TH2F *f2dHistAsMCResPtOmegaPlus;
         // Resolution in R(2D) as function of generated R
         TH2F *f2dHistAsMCResRXiMinus;
         TH2F *f2dHistAsMCResRXiPlus;
         TH2F *f2dHistAsMCResROmegaMinus;
         TH2F *f2dHistAsMCResROmegaPlus;
         // Resolution in phi as function of generated Pt
         TH2F *f2dHistAsMCResPhiXiMinus;
         TH2F *f2dHistAsMCResPhiXiPlus;
         TH2F *f2dHistAsMCResPhiOmegaMinus;
         TH2F *f2dHistAsMCResPhiOmegaPlus;
         // Correlation between proton (antiproton) daughter MC pt and Xi/Omega MC pt (to apply Geat/Fluka correction)
         TH2F *f2dHistAsMCptProtonMCptXiMinus;
         TH2F *f2dHistAsMCptAntiprotonMCptXiPlus;
         TH2F *f2dHistAsMCptProtonMCptOmegaMinus;
         TH2F *f2dHistAsMCptAntiprotonMCptOmegaPlus;
         // QA plots
         TH1F *fHistV0toXiCosineOfPointingAngle;
         TH2F *fHistV0CosineOfPointingAnglevsPtXi;
         TH2F *fHistV0CosineOfPointingAnglevsPtOmega;

         // Containers                       
         AliCFContainer  *fCFContCascadePIDAsXiMinus;
         AliCFContainer  *fCFContCascadePIDAsXiPlus;
         AliCFContainer  *fCFContCascadePIDAsOmegaMinus;
         AliCFContainer  *fCFContCascadePIDAsOmegaPlus;
         AliCFContainer  *fCFContAsCascadeCuts;



  AliAnalysisTaskCheckPerformanceCascadepp276(const AliAnalysisTaskCheckPerformanceCascadepp276&);            // not implemented
  AliAnalysisTaskCheckPerformanceCascadepp276& operator=(const AliAnalysisTaskCheckPerformanceCascadepp276&); // not implemented
  
  ClassDef(AliAnalysisTaskCheckPerformanceCascadepp276, 8);
};

#endif
