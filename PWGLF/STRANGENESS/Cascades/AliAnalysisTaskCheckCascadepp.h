#ifndef ALIANALYSISTASKCHECKCASCADEPP_H
#define ALIANALYSISTASKCHECKCASCADEPP_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//            AliAnalysisTaskCheckCascadePbPb class
//              Origin AliAnalysisTaskCheckCascade
//              This task has four roles :
//                1. QAing the Cascades from ESD and AOD
//                   Origin:  AliAnalysisTaskESDCheckV0 by Boris Hippolyte Nov2007, hippolyt@in2p3.fr
//                2. Prepare the plots which stand as raw material for yield extraction (wi/wo PID)
//                3. Supply an AliCFContainer meant to define the optimised topological selections
//                Adapted to Cascade : A.Maire Mar2008, antonin.maire@ires.in2p3.fr
//                Modified :           A.Maire Mar2010, antonin.maire@ires.in2p3.fr
//                Modified for PbPb analysis: M. Nicassio Feb 2011, maria.nicassio@ba.infn.it
//                Modified for pp@2.76 analysis: D. Colella Feb2012, domenico.colella@ba.infn.it
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class TH3F;
class TVector3;
class THnSparse;
 
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;
class AliPIDResponse;
class AliAnalysisUtils;

#include "TString.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCheckCascadepp : public AliAnalysisTaskSE {
 public:
        AliAnalysisTaskCheckCascadepp();
        AliAnalysisTaskCheckCascadepp(const char *name);
        virtual ~AliAnalysisTaskCheckCascadepp();
        virtual void   UserCreateOutputObjects();
        virtual void   UserExec(Option_t *option);
        virtual Int_t  DoESDTrackWithTPCrefitMultiplicity(const AliESDEvent *lESDevent);
        virtual void   Terminate(Option_t *);
        //Setters   
        void SetAnalysisType                  (const char* analysisType                 ) { fAnalysisType                   = analysisType;                 }
        void SetCollidingSystem               (Int_t  collidingSystem                   ) { fCollidingSystem                = collidingSystem;              }
        void SetSelectedTriggerClass          (AliVEvent::EOfflineTriggerTypes trigType ) { fkTriggerClass                  = trigType;                     }
        void SetEventSelSDDstatus             (Bool_t   eventselSDDstatus               ) { fApplyEvSelSDDstatus            = eventselSDDstatus;            }
        void SetEventSelDAQIncomplete         (Bool_t   eventselDAQincomplete           ) { fApplyEvSelDAQincomplete        = eventselDAQincomplete;        }
        void SetEventSelSPDclustervstracklet  (Bool_t   eventselSPDclustervstracklet    ) { fApplyEvSelSPDclustervstracklet = eventselSPDclustervstracklet; }
        void SetEventSelPileup                (Bool_t   eventselPileup                  ) { fApplyEvSelPileup               = eventselPileup;               }
        void SetEventSelPhysicsSel            (Bool_t   eventselPhysicsSel              ) { fApplyEvSelPhysicsSel           = eventselPhysicsSel;           }
        void SetEventSelNoTPConlyPrimVtx      (Bool_t   eventselNoTPConlyPrimVtx        ) { fApplyEvSelNoTPConlyPrimVtx     = eventselNoTPConlyPrimVtx;     }
        void SetEventSelSPDvtxres             (Bool_t   eventselSPDvtxres               ) { fApplyEvSelSPDvtxres            = eventselSPDvtxres;            }
        void SetEventSelVtxProximity          (Bool_t   eventselVtxProximity            ) { fApplyEvSelVtxProximity         = eventselVtxProximity;         }
        void SetEventSelZprimVtxPos           (Bool_t   eventselZprimVtxPos             ) { fApplyEvSelZprimVtxPos          = eventselZprimVtxPos;          }
        void SetRelaunchV0CascVertexers       (Bool_t   rerunV0CascVertexers            ) { fRerunV0CascVertexers           = rerunV0CascVertexers;         }
        void SetWithSDDOn                     (Bool_t   withsddOn                       ) { fwithSDD                        = withsddOn;                    }
        void SetExtraSelections               (Bool_t   extraSelections                 ) { fExtraSelections                = extraSelections;              }
        void SetTrackQualityCutTPCrefit       (Bool_t   trackqualityCutTPCrefit         ) { fTrackQualityCutTPCrefit        = trackqualityCutTPCrefit;      }
        void SetTrackQualityCutnTPCcls        (Bool_t   trackqualityCutnTPCcls          ) { fTrackQualityCutnTPCcls         = trackqualityCutnTPCcls;       }
        void SetQualityCutMinnTPCcls          (Int_t    minnTPCcls                      ) { fMinnTPCcls                     = minnTPCcls;                   }
        void SetQualityCutClusterOverFindable (Float_t minTPCcrossrawoverfindable       ) { fMinTPCcrossrawoverfindable     = minTPCcrossrawoverfindable;   }
        void SetVertexRange                   (Float_t vtxrangemin, Float_t vtxrangemax ) { fVtxRangeMax                    = vtxrangemax;     
                                                                                            fVtxRangeMin                    = vtxrangemin;                  }
        void SetMinptCutOnDaughterTracks      (Float_t minptdaughtrks                   ) { fMinPtCutOnDaughterTracks       = minptdaughtrks;               }
        void SetEtaCutOnDaughterTracks        (Float_t etadaughtrks                     ) { fEtaCutOnDaughterTracks         = etadaughtrks;                 }
        void SetSPDPileUpminContributors      (Int_t   spdpileupmincontributors         ) { fSPDPileUpminContributors       = spdpileupmincontributors;     }
        void SetNumTPCPIDsigma                (Double_t ftpcpidsigma                    ) { fTPCPIDsigma                    = ftpcpidsigma;                 }
        void SetSuffix                        (const char* suffix                       ) { fSuffix                         = suffix;                       }
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
        AliESDtrackCuts  *fESDtrackCuts;                   // ESD track cuts used for primary track definition
        AliAnalysisUtils *fUtils;                          // analysis utils (for pA vertex selection)

        Int_t             fCollidingSystem;                // "pPb" or "pp" colliding system
        AliVEvent::EOfflineTriggerTypes fkTriggerClass;    // Trigger selection: kMB, kINT7, etc as needed
        AliPIDResponse   *fPIDResponse;                    //! PID response object

        Bool_t            fApplyEvSelSDDstatus;            
        Bool_t            fApplyEvSelDAQincomplete;        //       
        Bool_t            fApplyEvSelSPDclustervstracklet; //
        Bool_t            fApplyEvSelPileup;               //
        Bool_t            fApplyEvSelPhysicsSel;           //
        Bool_t            fApplyEvSelNoTPConlyPrimVtx;     //
        Bool_t            fApplyEvSelSPDvtxres;            //
        Bool_t            fApplyEvSelVtxProximity;         //
        Bool_t            fApplyEvSelZprimVtxPos;          //
        Bool_t            fRerunV0CascVertexers;           // Boolean : kTRUE = relaunch both V0 + Cascade vertexers
        Bool_t            fwithSDD;                        // Boolean : kTRUE = select events with SDD reco
        Bool_t            fExtraSelections;                // Boolean : kTRUE = apply tighter selections, before starting the analysis
        Bool_t            fTrackQualityCutTPCrefit;        //
        Bool_t            fTrackQualityCutnTPCcls;         //
        Int_t             fMinnTPCcls;                     // Minimum number of TPC cluster for daughter tracks
        Float_t           fMinTPCcrossrawoverfindable;     // Minimum value for clusters/findable ratio 
        Float_t           fVtxRangeMax;                    // to select events with |zvtx|<fVtxRange cm
        Float_t           fVtxRangeMin;                    // to select events with |zvtx|>fVtxRangeMin cm
        Float_t           fMinPtCutOnDaughterTracks;       // minimum pt cut on daughter tracks
        Float_t           fEtaCutOnDaughterTracks;         // pseudorapidity cut on daughter tracks
        Int_t             fSPDPileUpminContributors;       //
        Double_t          fTPCPIDsigma;                    //
        TString           fSuffix;      
 
        Double_t          fV0Sels[7];                      // Array to store the 7 values for the different selections V0 related (if fkRerunV0CascVertexers)
        Double_t          fCascSels[8];                    // Array to store the 8 values for the different selections Casc. related (if fkRerunV0CascVertexers)

        TList      *fListHistCascade;                   //! List of Cascade histograms
        // Cascades multiplicity plots
        TH1F *fHistCascadeMultiplicityBeforeAnySel;
        TH1F *fHistCascadeMultiplicityAfterSDDstatusSel;
        TH1F *fHistCascadeMultiplicityAfterDAQincompleteEvRej;
        TH1F *fHistCascadeMultiplicityAfterSPDclustervstrackletSel;
        TH1F *fHistCascadeMultiplicityAfterPileupRej;
        TH1F *fHistCascadeMultiplicityAfterPhysicsSel;
        TH1F *fHistCascadeMultiplicityAfterRevertexing;
        TH1F *fHistCascadeMultiplicityAfterNoTPConlyPrimVtxSel;
        TH1F *fHistCascadeMultiplicityAfterSPDresolution;
        TH1F *fHistCascadeMultiplicityAfterVerticesProximity;
        TH1F *fHistCascadeMultiplicityAfterZprimVtxPosSel;
        // Vertex position plots (BestVertex)
        TH1F   *fHistPVx;                                             //! Best primary vertex X position distribution after all evnt selection
        TH1F   *fHistPVy;                                             //! Best primary vertex Y position distribution after all evnt selection
        TH1F   *fHistPVz;                                             //! Best primary vertex Z position distribution after all evnt selection
        TH1F   *fHistPVxAnalysis;                                     //! Best primary vertex X position distribution after all evnt selection and |z|>10cm cut
        TH1F   *fHistPVyAnalysis;                                     //! Best primary vertex Y position distribution after all evnt selection and |z|>10cm cut    
        TH1F   *fHistPVzAnalysis;                                     //! Best primary vertex Z position distribution after all evnt selection and |z|>10cm cut
        // TPC cluster distributions for daughters
        TH1F   *fHistPosV0TPCClusters;                                //! TPC clusters distribution for Positive V0 daughter track
        TH1F   *fHistNegV0TPCClusters;                                //! TPC clusters distribution for Negative V0 daughter track
        TH1F   *fHistBachTPCClusters;                                 //! TPC clusters distribution for Bachelor V0 daughter track
        // Cut's variables distributions
        TH1F   *fHistEffMassXi;                                       //! reconstructed cascade effective mass
        TH1F   *fHistDcaXiDaughters;                                  //! dca between Xi's daughters
        TH1F   *fHistDcaBachToPrimVertex;                             //! dca of the bachelor track to primary vertex
        TH1F   *fHistXiCosineOfPointingAngle;                         //! cosine of Xi pointing angle in a cascade
        TH1F   *fHistXiRadius;                                        //! (transverse) radius of the cascade vertex
        TH1F   *fHistMassLambdaAsCascDghter;                          //! Test Invariant Mass of Lambda coming from Cascade 
        TH1F   *fHistDcaV0DaughtersXi;                                //! Dca between V0 daughters, for the V0 associated to a cascade
        TH1F   *fHistDcaV0ToPrimVertexXi;                             //! Dca of V0 to primary vertex, for the V0 associated to a cascade
        TH1F   *fHistV0CosineOfPointingAngleXi;                       //! Cosine of V0 pointing angle, for the V0 associated to a cascade
        TH1F   *fHistV0RadiusXi;                                      //! V0 (transverse) distance distribution, for the V0 associated to a cascade
        TH1F   *fHistDcaPosToPrimVertexXi;                            //! Dca of V0 positive daughter to primary vertex, for the V0 associated to a cascade
        TH1F   *fHistDcaNegToPrimVertexXi;                            //! Dca of V0 negative daughter to primary vertex, for the V0 associated to a cascade
        // Invariant mass distributions
        TH1F   *fHistMassXiMinus;                                     //! reconstructed cascade effective mass, under Xi- hyp.
        TH1F   *fHistMassXiPlus;                                      //! reconstructed cascade effective mass, under Xi+ hyp.
        TH1F   *fHistMassOmegaMinus;                                  //! reconstructed cascade effective mass, under Omega- hyp.
        TH1F   *fHistMassOmegaPlus;                                   //! reconstructed cascade effective mass, under Omega+ hyp.
        // Transverse and total momentum distributions
        TH1F   *fHistXiTransvMom;                                     //! Xi transverse momentum, around the mass peak of Xi-/+
        TH1F   *fHistXiTotMom;                                        //! Xi momentum norm, around the mass peak of Xi-/+
        TH1F   *fHistBachTransvMomXi;                                 //! bachelor transverse momentum, for cand. around the mass peak of Xi-/+
        TH1F   *fHistBachTotMomXi;                                    //! bachelor momentum norm, for cand. around the mass peak of Xi-/+
        // Others QA plots
        TH1F   *fHistChargeXi;                                        //! Charge sign of the cascade candidate
        TH1F   *fHistV0toXiCosineOfPointingAngle;                     //! Cos. of Pointing angle between the V0 mom and the Xi-V0 vtx line
        TH1F   *fHistRapXi;                                           //! rapidity of Xi candidates, around the mass peak of Xi-/+
        TH1F   *fHistRapOmega;                                        //! rapidity of Omega candidates, around the mass peak of Omega-/+
        TH1F   *fHistEtaXi;                                           //! eta distrib. of all the cascade candidates, around the mass peak of Xi-/+
        TH1F   *fHistEtaBachXi;                                       
        TH1F   *fHistEtaPosXi;                                        
        TH1F   *fHistEtaNegXi;                                        
        TH1F   *fHistThetaXi;                                         //! theta distrib. of all the cascade candidates, around the mass peak of Xi-/+
        TH1F   *fHistPhiXi;                                           //! phi distrib. of all the cascade candidates, around the mass peak of Xi-/+
        TH2F   *f2dHistArmenteros;                                    //! alpha(casc. cand.) Vs PtArm(casc. cand.)
        TH2F   *f2dHistEffMassLambdaVsEffMassXiMinus;                 //! Xi- Eff mass Vs V0 Eff mass, under Xi- hyp.
        TH2F   *f2dHistEffMassXiVsEffMassOmegaMinus;                  //! Xi- Eff mass Vs Omega- Eff mass, for negative cascades
        TH2F   *f2dHistEffMassLambdaVsEffMassXiPlus;                  //! Xi+ Eff mass Vs V0 Eff mass, under Xi+ hyp. 
        TH2F   *f2dHistEffMassXiVsEffMassOmegaPlus;                   //! Xi+ Eff mass Vs Omega+ Eff mass, for positive cascades
        TH2F   *f2dHistXiRadiusVsEffMassXiMinus;                      //! transv. casc. decay radius Vs Xi- Eff mass, under Xi- hyp.
        TH2F   *f2dHistXiRadiusVsEffMassXiPlus;                       //! transv. casc. decay radius Vs Xi+ Eff mass, under Xi+ hyp.
        TH2F   *f2dHistXiRadiusVsEffMassOmegaMinus;                   //! transv. casc. decay radius Vs Omega- Eff mass, under Omega- hyp.
        TH2F   *f2dHistXiRadiusVsEffMassOmegaPlus;                    //! transv. casc. decay radius Vs Omega+ Eff mass, under Omega+ hyp.
        TH2F   *f2dHistTPCdEdxOfCascDghters;                          //! TPC Bethe-Bloch curve, populated with the cascade daughters
        TH2F   *f2dHistDcaXiDaughtersvsInvMass;                       //! cut variables vs inv. mass
        TH2F   *f2dHistDcaBachToPrimVertexvsInvMass;                  //! cut variables vs inv. mass
        TH2F   *f2dHistXiCosineOfPointingAnglevsInvMass;              //! cut variables vs inv. mass
        TH2F   *f2dHistMassLambdaAsCascDghtervsInvMass;               //! cut variables vs inv. mass
        TH2F   *f2dHistDcaV0DaughtersXivsInvMass;                     //! cut variables vs inv. mass 
        TH2F   *f2dHistDcaV0ToPrimVertexXivsInvMass;                  //! cut variables vs inv. mass 
        // Containers for cuts study 
        AliCFContainer  *fCFContCascadePIDXiMinus;                       //! for Xi-   : Container to store any 3D histos with the different PID flavours
        AliCFContainer  *fCFContCascadePIDXiPlus;                        //! for Xi+   : Container to store any 3D histos with the different PID flavours
        AliCFContainer  *fCFContCascadePIDOmegaMinus;                    //! for Omega-: Container to store any 3D histos with the different PID flavours
        AliCFContainer  *fCFContCascadePIDOmegaPlus;                     //! for Omega+: Container to store any 3D histos with the different PID flavours
        AliCFContainer  *fCFContCascadeCuts;                             //! Container meant to store all the relevant distributions corresponding to the cut variables


  AliAnalysisTaskCheckCascadepp(const AliAnalysisTaskCheckCascadepp&);            // not implemented
  AliAnalysisTaskCheckCascadepp& operator=(const AliAnalysisTaskCheckCascadepp&); // not implemented
  
  ClassDef(AliAnalysisTaskCheckCascadepp, 11);
};

#endif
