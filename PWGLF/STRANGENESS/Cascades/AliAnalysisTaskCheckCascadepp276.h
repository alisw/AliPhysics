#ifndef ALIANALYSISTASKCHECKCASCADEPP276_H
#define ALIANALYSISTASKCHECKCASCADEPP276_H

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

class AliAnalysisTaskCheckCascadepp276 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskCheckCascadepp276();
  AliAnalysisTaskCheckCascadepp276(const char *name);
  virtual ~AliAnalysisTaskCheckCascadepp276();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual Int_t  DoESDTrackWithTPCrefitMultiplicity(const AliESDEvent *lESDevent);
         //virtual Int_t  Tracks2V0vertices(AliESDEvent *event);  
         //virtual Int_t  V0sTracks2CascadeVertices(AliESDEvent *event); 
         //virtual Double_t Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const;
         //virtual Double_t Det(Double_t a00,Double_t a01,Double_t a02,
         //      Double_t a10,Double_t a11,Double_t a12,
         //      Double_t a20,Double_t a21,Double_t a22) const;

         //virtual Double_t PropagateToDCA(AliESDv0 *vtx,AliExternalTrackParam *trk,Double_t b);
  virtual void   Terminate(Option_t *);
  
  void SetAnalysisType               (const char* analysisType          = "ESD"  ) { fAnalysisType                = analysisType;               }
  void SetCollidingSystem            (const char* collidingSystem       = "pp"   ) { fCollidingSystem             = collidingSystem;            }
  void SetRelaunchV0CascVertexers    (Bool_t rerunV0CascVertexers       = kFALSE ) { fkRerunV0CascVertexers       = rerunV0CascVertexers;       }
  void SetSDDSelection               (Bool_t sddOnSelection             = kTRUE  ) { fkSDDSelectionOn             = sddOnSelection;             }
  void SetQualityCutZprimVtxPos      (Bool_t qualityCutZprimVtxPos      = kTRUE  ) { fkQualityCutZprimVtxPos      = qualityCutZprimVtxPos;      }
  void SetQualityCutNoTPConlyPrimVtx (Bool_t qualityCutNoTPConlyPrimVtx = kTRUE  ) { fkQualityCutNoTPConlyPrimVtx = qualityCutNoTPConlyPrimVtx; }
  void SetQualityCutTPCrefit         (Bool_t qualityCutTPCrefit         = kTRUE  ) { fkQualityCutTPCrefit         = qualityCutTPCrefit;         }
  void SetQualityCutnTPCcls          (Bool_t qualityCutnTPCcls          = kTRUE  ) { fkQualityCutnTPCcls          = qualityCutnTPCcls;          }
  void SetQualityCutPileup           (Bool_t qualityCutPileup           = kTRUE  ) { fkQualityCutPileup           = qualityCutPileup;           }
  void SetWithSDDOn                  (Bool_t withsddOn                  = kTRUE  ) { fwithSDD                     = withsddOn;                  }
  void SetQualityCutMinnTPCcls       (Int_t  minnTPCcls                 = 70     ) { fMinnTPCcls                  = minnTPCcls;                 }
  void SetExtraSelections            (Bool_t extraSelections            = kFALSE ) { fkExtraSelections            = extraSelections;            }
  void SetVertexRange                (Float_t vtxrange                  = 10.0   ) { fVtxRange                    = vtxrange;                   }
  void SetVertexRangeMin             (Float_t vtxrangemin               = 0.0    ) { fVtxRangeMin                 = vtxrangemin;                }
  void SetMinptCutOnDaughterTracks   (Float_t minptdaughtrks            = 0.0    ) { fMinPtCutOnDaughterTracks    = minptdaughtrks;             }
  void SetEtaCutOnDaughterTracks     (Float_t etadaughtrks              = 0.8    ) { fEtaCutOnDaughterTracks      = etadaughtrks;               }

 private:
        // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
        // your data member object is created on the worker nodes and streaming is not needed.
        // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14


        TString          fAnalysisType;                  // "ESD" or "AOD" analysis type	
       // TString          fCollidingSystem;               // "pPb" or "pp" colliding system
        AliESDtrackCuts  *fESDtrackCuts;                 // ESD track cuts used for primary track definition
       // AliPIDResponse   *fPIDResponse;                  //! PID response object
        AliAnalysisUtils *fUtils;                        // analysis utils (for pA vertex selection)
        TString          fCollidingSystem;               // "pPb" or "pp" colliding system
        AliPIDResponse   *fPIDResponse;                  //! PID response object

        Bool_t          fkRerunV0CascVertexers;         // Boolean : kTRUE = relaunch both V0 + Cascade vertexers
        Bool_t          fkSDDSelectionOn;               // Boolena : kTRUE = select events with SDD on
        Bool_t          fkQualityCutZprimVtxPos;        // Boolean : kTRUE = cut on the prim.vtx  z-position
        Bool_t          fkQualityCutNoTPConlyPrimVtx;   // Boolean : kTRUE = prim vtx should be SPD or Tracking vertex
        Bool_t          fkQualityCutTPCrefit;           // Boolean : kTRUE = ask for TPCrefit for the 3 daughter tracks
        Bool_t          fkQualityCutnTPCcls;            // Boolean : kTRUE = ask for fMinnTPCcls TPC clusters for each daughter track
        Bool_t          fkQualityCutPileup;             // Boolean : kTRUE = ask for No Pileup events
        Bool_t          fwithSDD;                       // Boolean : kTRUE = select events with SDD reco
        Int_t           fMinnTPCcls;                    // Minimum number of TPC cluster for daughter tracks
        Bool_t          fkExtraSelections;              // Boolean : kTRUE = apply tighter selections, before starting the analysis
        Float_t         fVtxRange;                      // to select events with |zvtx|<fVtxRange cm
        Float_t         fVtxRangeMin;                   // to select events with |zvtx|>fVtxRangeMin cm
        Float_t         fMinPtCutOnDaughterTracks;      // minimum pt cut on daughter tracks
        Float_t         fEtaCutOnDaughterTracks;        // pseudorapidity cut on daughter tracks
       
        Double_t        fV0Sels[7];                     // Array to store the 7 values for the different selections V0 related (if fkRerunV0CascVertexers)
        Double_t        fCascSels[8];                   // Array to store the 8 values for the different selections Casc. related (if fkRerunV0CascVertexers)

        TList      *fListHistCascade;                   //! List of Cascade histograms
        
        // Cascades multiplicity plots
        TH1F   *fHistCascadeMultiplicityBeforeAnySel;                 //! Cascade multiplicity distribution before any evnt selection 
        TH1F   *fHistCascadeMultiplicityAfterSDDSel;                  //! Cascade multiplicity distribution after evnt selection on the SDD
        TH1F   *fHistCascadeMultiplicityAfterPhysicsSel;              //! Cascade multiplicity distribution after evnt Physics Selection  
        TH1F   *fHistCascadeMultiplicityForSelEvtNoTPCOnly;           //! Cascade multiplicity distribution after evnt noTPCOnly selection
        TH1F   *fHistCascadeMultiplicityForSelEvtNoTPCOnlyNoPileup;   //! Cascade multiplicity distribution after evnt PileUp selection
        TH1F   *fHistCascadeMultiplicityAfterVertexCutSel;            //! Cascade multiplicity distribution after evnt selection on the Z vertex position cut
        // Tracks multiplicity plots
        TH1F   *fHistTrackMultiplicityBeforeAnySel;                   //! Track multiplicity distribution before any evnt selection  
        TH1F   *fHistTrackMultiplicityAfterSDDSel;                    //! Track multiplicity distribution after evnt selection on the SDD
        TH1F   *fHistTrackMultiplicityAfterPhysicsSel;                //! Track multiplicity distribution after evnt Physics Selection
        TH1F   *fHistTrackMultiplicityForSelEvtNoTPCOnly;             //! Track multiplicity distribution after evnt noTPCOnly selection
        TH1F   *fHistTrackMultiplicityForSelEvtNoTPCOnlyNoPileup;     //! Track multiplicity distributionafter evnt PileUp selection
        TH1F   *fHistTrackMultiplicityAfterVertexCutSel;              //! Track multiplicity distribution after evnt selection on the Z vertex position cut
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


  AliAnalysisTaskCheckCascadepp276(const AliAnalysisTaskCheckCascadepp276&);            // not implemented
  AliAnalysisTaskCheckCascadepp276& operator=(const AliAnalysisTaskCheckCascadepp276&); // not implemented
  
  ClassDef(AliAnalysisTaskCheckCascadepp276, 8);
};

#endif
