#ifndef ALIANALYSISTASKCHECKPERFORMANCECASCADEPP_H
#define ALIANALYSISTASKCHECKPERFORMANCECASCADEPP_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//        AliAnalysisTaskCheckPerformanceCascadePbPb class
//            This task is for a performance study of cascade identification.
//            It works with MC info and ESD and AOD tree 
//            Origin   : A.Maire    Jan2010, antonin.maire@ires.in2p3.fr
//            Modified : M.Nicassio Feb2011, maria.nicassio@ba.infn.it
//            Modified : D.Colella  Feb2012, domenico.colella@ba.infn.it
//-----------------------------------------------------------------

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

class AliAnalysisTaskCheckPerformanceCascadepp : public AliAnalysisTaskSE {
 public:
 
  AliAnalysisTaskCheckPerformanceCascadepp();
  AliAnalysisTaskCheckPerformanceCascadepp(const char *name );
  virtual ~AliAnalysisTaskCheckPerformanceCascadepp();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAnalysisType                 (const char* analysisType            = "ESD"  ) { fAnalysisType                   = analysisType;                 }
  void SetCollidingSystem              (const char* collidingSystem         = "pp"   ) { fCollidingSystem                = collidingSystem;              }
  void SetSelectedTriggerClass         (AliVEvent::EOfflineTriggerTypes trigType     ) { fkTriggerClass                  = trigType;                     }
  void SetEventSelSDDstatus            (Bool_t eventselSDDstatus            = kFALSE ) { fApplyEvSelSDDstatus            = eventselSDDstatus;            }
  void SetEventSelDAQIncomplete        (Bool_t eventselDAQincomplete        = kTRUE  ) { fApplyEvSelDAQincomplete        = eventselDAQincomplete;        }
  void SetEventSelSPDclustervstracklet (Bool_t eventselSPDclustervstracklet = kTRUE  ) { fApplyEvSelSPDclustervstracklet = eventselSPDclustervstracklet; }
  void SetEventSelPileup               (Bool_t eventselPileup               = kTRUE  ) { fApplyEvSelPileup               = eventselPileup;               }
  void SetEventSelPhysicsSel           (Bool_t eventselPhysicsSel           = kTRUE  ) { fApplyEvSelPhysicsSel           = eventselPhysicsSel;           }
  void SetEventSelNoTPConlyPrimVtx     (Bool_t eventselNoTPConlyPrimVtx     = kTRUE  ) { fApplyEvSelNoTPConlyPrimVtx     = eventselNoTPConlyPrimVtx;     }
  void SetEventSelSPDvtxres            (Bool_t eventselSPDvtxres            = kTRUE  ) { fApplyEvSelSPDvtxres            = eventselSPDvtxres;            }
  void SetEventSelVtxProximity         (Bool_t eventselVtxProximity         = kTRUE  ) { fApplyEvSelVtxProximity         = eventselVtxProximity;         }
  void SetEventSelZprimVtxPos          (Bool_t eventselZprimVtxPos          = kTRUE  ) { fApplyEvSelZprimVtxPos          = eventselZprimVtxPos;          }
  void SetRelaunchV0CascVertexers      (Bool_t rerunV0CascVertexers         = kFALSE ) { fRerunV0CascVertexers           = rerunV0CascVertexers;         }
  void SetWithSDDOn                    (Bool_t withsddOn                    = kTRUE  ) { fwithSDD                        = withsddOn;                    }
  void SetTrackQualityCutTPCrefit      (Bool_t trackqualityCutTPCrefit      = kTRUE  ) { fTrackQualityCutTPCrefit        = trackqualityCutTPCrefit;      }
  void SetTrackQualityCutnTPCcls       (Bool_t trackqualityCutnTPCcls       = kTRUE  ) { fTrackQualityCutnTPCcls         = trackqualityCutnTPCcls;       }
  void SetQualityCutMinnTPCcls         (Int_t  minnTPCcls                   = 70     ) { fMinnTPCcls                     = minnTPCcls;                   }
  void SetVertexRange                  (Float_t vtxrangemin, Float_t vtxrangemax     ) { fVtxRangeMax = vtxrangemax;     fVtxRangeMin = vtxrangemin;     }
  void SetApplyAccCut                  (Bool_t  acccut                      = kFALSE ) { fApplyAccCut                    = acccut;                       }
  void SetMinptCutOnDaughterTracks     (Float_t minptdaughtrks              = 0.0    ) { fMinPtCutOnDaughterTracks       = minptdaughtrks;               }
  void SetExtraSelections              (Bool_t  extraSelections             = kFALSE ) { fkExtraSelections               = extraSelections;              }
  void SetEtaCutOnDaughterTracks       (Float_t etadaughtrks                = 0.8    ) { fEtaCutOnDaughterTracks         = etadaughtrks;                 }
  void SetSPDPileUpminContributors     (Int_t   spdpileupmincontributors    = 3      ) { fSPDPileUpminContributors       = spdpileupmincontributors;     }
  //Setters for the V0 and cascade Vertexer Parameters
  void SetV0VertexerMaxChisquare           (Double_t lParameter){ fV0Sels[0] = lParameter; }
  void SetV0VertexerDCAFirstToPV           (Double_t lParameter){ fV0Sels[1] = lParameter; }
  void SetV0VertexerDCASecondtoPV          (Double_t lParameter){ fV0Sels[2] = lParameter; }
  void SetV0VertexerDCAV0Daughters         (Double_t lParameter){ fV0Sels[3] = lParameter; }
  void SetV0VertexerCosinePA               (Double_t lParameter){ fV0Sels[4] = lParameter; }
  void SetV0VertexerMinRadius              (Double_t lParameter){ fV0Sels[5] = lParameter; }
  void SetV0VertexerMaxRadius              (Double_t lParameter){ fV0Sels[6] = lParameter; }
  void SetCascVertexerMaxChisquare         (Double_t lParameter){ fCascSels[0] = lParameter; }
  void SetCascVertexerMinV0ImpactParameter (Double_t lParameter){ fCascSels[1] = lParameter; }
  void SetCascVertexerV0MassWindow         (Double_t lParameter){ fCascSels[2] = lParameter; }
  void SetCascVertexerDCABachToPV          (Double_t lParameter){ fCascSels[3] = lParameter; }
  void SetCascVertexerDCACascadeDaughters  (Double_t lParameter){ fCascSels[4] = lParameter; }
  void SetCascVertexerCascadeCosinePA      (Double_t lParameter){ fCascSels[5] = lParameter; }
  void SetCascVertexerCascadeMinRadius     (Double_t lParameter){ fCascSels[6] = lParameter; }
  void SetCascVertexerCascadeMaxRadius     (Double_t lParameter){ fCascSels[7] = lParameter; }

 private:
        // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
        // your data member object is created on the worker nodes and streaming is not needed.
        // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14

        TString           fAnalysisType;                   // "ESD" or "AOD" analysis type      
        TString           fCollidingSystem;                // "pPb" or "pp" colliding system
        AliVEvent::EOfflineTriggerTypes fkTriggerClass;    // Trigger selection: kMB, kINT7, etc as needed
        Bool_t            fApplyEvSelSDDstatus;            //
        Bool_t            fApplyEvSelDAQincomplete;        //       
        Bool_t            fApplyEvSelSPDclustervstracklet; //
        Bool_t            fApplyEvSelPileup;               //
        Bool_t            fApplyEvSelPhysicsSel;           //
        Bool_t            fApplyEvSelNoTPConlyPrimVtx;     //
        Bool_t            fApplyEvSelSPDvtxres;            //
        Bool_t            fApplyEvSelVtxProximity;         //
        Bool_t            fApplyEvSelZprimVtxPos;          //
        AliESDtrackCuts  *fESDtrackCuts;                   // ESD track cuts used for primary track definition
        AliAnalysisUtils *fUtils;                          // analysis utils (for pA vertex selection)
        AliPIDResponse   *fPIDResponse;                    //! PID response object
        Bool_t            fRerunV0CascVertexers;           // Boolean : kTRUE = relaunch both V0 + Cascade vertexers
        Bool_t            fwithSDD;                        // Boolean : kTRUE = select events with SDD reco
        Bool_t            fExtraSelections;                // Boolean : kTRUE = apply tighter selections, before starting the analysis
        Bool_t            fTrackQualityCutTPCrefit;        //
        Bool_t            fTrackQualityCutnTPCcls;         //
        Int_t             fMinnTPCcls;                     // Minimum number of TPC cluster for daughter tracks
        Bool_t            fkExtraSelections;               // Boolean : kTRUE = apply tighter selections, before starting the analysis
        Float_t           fVtxRangeMax;                    // to select events with |zvtx|<fVtxRange cm
        Float_t           fVtxRangeMin;                    // to select events with |zvtx|>fVtxRangeMin cm
        Bool_t            fApplyAccCut;                    // flag to apply acceptance cuts to MC cascades
        Float_t           fMinPtCutOnDaughterTracks;       // minimum pt cut on daughter tracks
        Float_t           fEtaCutOnDaughterTracks;         // pseudorapidity cut on daughter tracks
        Int_t             fSPDPileUpminContributors;       // 

        Double_t        fV0Sels[7];                     // Array to store the 7 values for the different selections V0 related (if fkRerunV0CascVertexers)
        Double_t        fCascSels[8];                   // Array to store the 8 values for the different selections Casc. related (if fkRerunV0CascVertexers)
	
 	TList	*fListHistCascade;		        //! List of Cascade histograms
         // - General Plots
         // Cascade multiplicity plots
         TH1F *fHistCascadeMultiplicityBeforeAnySel;
         TH1F *fHistCascadeMultiplicityAfterSDDstatusSel;
         TH1F *fHistCascadeMultiplicityAfterDAQincompleteEvRej;
         TH1F *fHistCascadeMultiplicityAfterSPDclustervstrackletSel;
         TH1F *fHistCascadeMultiplicityAfterPileupRej;
         TH1F *fHistCascadeMultiplicityAfterPhysicsSel;
         TH1F *fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel;
         TH1F *fHistCascadeMultiplicityAfterZprimVtxPosSel;
         // Tracks multiplicity plots
         TH1F *fHistTrackMultiplicityBeforeAnySel;
         TH1F *fHistTrackMultiplicityAfterSDDstatusSel;
         TH1F *fHistTrackMultiplicityAfterDAQincompleteEvRej;
         TH1F *fHistTrackMultiplicityAfterSPDclustervstrackletSel;
         TH1F *fHistTrackMultiplicityAfterPileupRej;
         TH1F *fHistTrackMultiplicityAfterPhysicsSel;
         TH1F *fHistTrackMultiplicityAfterNoTPConlyPrimVtxSel;
         TH1F *fHistTrackMultiplicityAfterZprimVtxPosSel;
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
         // Vertex position plots (BestVertex)
         TH1F *fHistPVx;                              // After any selections but before |Z| < 10 cm
         TH1F *fHistPVy;                              // After any selections but before |Z| < 10 cm
         TH1F *fHistPVz;                              // After any selections but before |Z| < 10 cm
         TH1F *fHistPVxAnalysis;                      // After any event selections 
         TH1F *fHistPVyAnalysis;                      // After any event selections
         TH1F *fHistPVzAnalysis;                      // After any event selections
         // - Plots before Physics Selection
         TH3D *f3dHistGenPtVsGenYvsNtracksXiMinus_A;    // Before any event selection 
         TH3D *f3dHistGenPtVsGenctauvsYXiMinus_A;       // Before any event selection 
         TH3D *f3dHistGenPtVsGenYvsNtracksXiPlus_A;     // Before any event selection 
         TH3D *f3dHistGenPtVsGenctauvsYXiPlus_A;        // Before any event selection 
         TH3D *f3dHistGenPtVsGenYvsNtracksOmegaMinus_A; // Before any event selection 
         TH3D *f3dHistGenPtVsGenctauvsYOmegaMinus_A;    // Before any event selection 
         TH3D *f3dHistGenPtVsGenYvsNtracksOmegaPlus_A;  // Before any event selection 
         TH3D *f3dHistGenPtVsGenctauvsYOmegaPlus_A;     // Before any event selection
         TH3D *f3dHistGenPtVsGenYvsNtracksXiMinus_B;    // After preliminary event selections (IsIncompleteDAQ, Tracklet vs Clusters Cut, SPD Pileup 
         TH3D *f3dHistGenPtVsGenctauvsYXiMinus_B;       // After preliminary event selections (IsIncompleteDAQ, Tracklet vs Clusters Cut, SPD Pileup
         TH3D *f3dHistGenPtVsGenYvsNtracksXiPlus_B;     // After preliminary event selections (IsIncompleteDAQ, Tracklet vs Clusters Cut, SPD Pileup
         TH3D *f3dHistGenPtVsGenctauvsYXiPlus_B;        // After preliminary event selections (IsIncompleteDAQ, Tracklet vs Clusters Cut, SPD Pileup
         TH3D *f3dHistGenPtVsGenYvsNtracksOmegaMinus_B; // After preliminary event selections (IsIncompleteDAQ, Tracklet vs Clusters Cut, SPD Pileup
         TH3D *f3dHistGenPtVsGenctauvsYOmegaMinus_B;    // After preliminary event selections (IsIncompleteDAQ, Tracklet vs Clusters Cut, SPD Pileup
         TH3D *f3dHistGenPtVsGenYvsNtracksOmegaPlus_B;  // After preliminary event selections (IsIncompleteDAQ, Tracklet vs Clusters Cut, SPD Pileup
         TH3D *f3dHistGenPtVsGenctauvsYOmegaPlus_B;     // After preliminary event selections (IsIncompleteDAQ, Tracklet vs Clusters Cut, SPD Pileup
         TH3D *f3dHistGenPtVsGenYvsNtracksXiMinus_C;    // After physics selection
         TH3D *f3dHistGenPtVsGenctauvsYXiMinus_C;       // After physics selection
         TH3D *f3dHistGenPtVsGenYvsNtracksXiPlus_C;     // After physics selection
         TH3D *f3dHistGenPtVsGenctauvsYXiPlus_C;        // After physics selection 
         TH3D *f3dHistGenPtVsGenYvsNtracksOmegaMinus_C; // After physics selection 
         TH3D *f3dHistGenPtVsGenctauvsYOmegaMinus_C;    // After physics selection 
         TH3D *f3dHistGenPtVsGenYvsNtracksOmegaPlus_C;  // After physics selection 
         TH3D *f3dHistGenPtVsGenctauvsYOmegaPlus_C;     // After physics selection
         TH3D *f3dHistGenPtVsGenYvsNtracksXiMinus_D;    // After all event selections 
         TH3D *f3dHistGenPtVsGenctauvsYXiMinus_D;       // After all event selections
         TH3D *f3dHistGenPtVsGenYvsNtracksXiPlus_D;     // After all event selections
         TH3D *f3dHistGenPtVsGenctauvsYXiPlus_D;        // After all event selections 
         TH3D *f3dHistGenPtVsGenYvsNtracksOmegaMinus_D; // After all event selections 
         TH3D *f3dHistGenPtVsGenctauvsYOmegaMinus_D;    // After all event selections 
         TH3D *f3dHistGenPtVsGenYvsNtracksOmegaPlus_D;  // After all event selections 
         TH3D *f3dHistGenPtVsGenctauvsYOmegaPlus_D;     // After all event selections 
         // - Generated cascade plots
         // After all the event selections 
         //Xi-
         TH1F *fHistEtaGenCascXiMinus;                // In all the eta and pt range (as they are generated)
         TH1F *fHistThetaGenCascXiMinus;              // In all the eta and pt range (as they are generated)
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



  AliAnalysisTaskCheckPerformanceCascadepp(const AliAnalysisTaskCheckPerformanceCascadepp&);            // not implemented
  AliAnalysisTaskCheckPerformanceCascadepp& operator=(const AliAnalysisTaskCheckPerformanceCascadepp&); // not implemented
  
  ClassDef(AliAnalysisTaskCheckPerformanceCascadepp, 11);
};

#endif
