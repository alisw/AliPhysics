#ifndef ALIANALYSISTASKCHECKCASCADEPBPB_H
#define ALIANALYSISTASKCHECKCASCADEPBPB_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//            AliQAProdMultistrange class
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
 
class AliESDEvent;
class AliPhysicsSelection;
class AliCFContainer;
class AliPIDResponse;

#include "TString.h"

#include "AliAnalysisTaskSE.h"

class AliQAProdMultistrange : public AliAnalysisTaskSE {
 public:
  AliQAProdMultistrange();
  AliQAProdMultistrange(const char *name);
  virtual ~AliQAProdMultistrange();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetAnalysisType               (const char* analysisType          = "ESD" ) { fAnalysisType                = analysisType;               }
  void SetCollidingSystem            (const char* collidingSystem       = "PbPb") { fCollidingSystem             = collidingSystem;            }
  void SetSDDselection               (Bool_t SDDSelection               = kFALSE) { fkSDDSelectionOn             = SDDSelection;               }
  void SetQualityCutZprimVtxPos      (Bool_t qualityCutZprimVtxPos      = kTRUE ) { fkQualityCutZprimVtxPos      = qualityCutZprimVtxPos;      }
  void SetQualityCutNoTPConlyPrimVtx (Bool_t qualityCutNoTPConlyPrimVtx = kTRUE ) { fkQualityCutNoTPConlyPrimVtx = qualityCutNoTPConlyPrimVtx; }
  void SetQualityCutTPCrefit         (Bool_t qualityCutTPCrefit         = kTRUE ) { fkQualityCutTPCrefit         = qualityCutTPCrefit;         }
  void SetQualityCutnTPCcls          (Bool_t qualityCutnTPCcls          = kTRUE ) { fkQualityCutnTPCcls          = qualityCutnTPCcls;          }
  void SetQualityCutMinnTPCcls       (Int_t  minnTPCcls                 = 70    ) { fMinnTPCcls                  = minnTPCcls;                 }
  void SetQualityCutPileup           (Bool_t qualitycutPileup           = kFALSE) { fkQualityCutPileup           = qualitycutPileup;           }
  void SetwithSDD                    (Bool_t withSDD                    = kTRUE ) { fwithSDD                     = withSDD;                    } 
  void SetCentralityLowLim           (Float_t centrlowlim               = 0.    ) { fCentrLowLim                 = centrlowlim;                }  
  void SetCentralityUpLim            (Float_t centruplim                = 100.  ) { fCentrUpLim                  = centruplim;                 }
  void SetCentralityEst              (TString   centrest                = "V0M" ) { fCentrEstimator              = centrest;                   }
  void SetUseCleaning                (Bool_t   usecleaning              = kTRUE ) { fkUseCleaning                = usecleaning;                }
  void SetVertexRange                (Float_t vtxrange                  = 0.    ) { fVtxRange                    = vtxrange;                   }
  void SetMinptCutOnDaughterTracks   (Float_t minptdaughtrks            = 0.    ) { fMinPtCutOnDaughterTracks    = minptdaughtrks;             }
  void SetEtaCutOnDaughterTracks     (Float_t etadaughtrks              = 0.    ) { fEtaCutOnDaughterTracks      = etadaughtrks;               }

 private:
        // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
        // your data member object is created on the worker nodes and streaming is not needed.
        // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14


        TString         fAnalysisType;                  // "ESD" or "AOD" analysis type	
        AliESDtrackCuts *fESDtrackCuts;                 // ESD track cuts used for primary track definition
        TString         fCollidingSystem;               // "PbPb" or "pp" colliding system
        AliPIDResponse *fPIDResponse;                   //! PID response object
        Bool_t          fkSDDSelectionOn;               // Boolean : kTRUE = apply the selection on SDD status
        Bool_t          fkQualityCutZprimVtxPos;        // Boolean : kTRUE = cut on the prim.vtx  z-position
        Bool_t          fkQualityCutNoTPConlyPrimVtx;   // Boolean : kTRUE = prim vtx should be SPD or Tracking vertex
        Bool_t          fkQualityCutTPCrefit;           // Boolean : kTRUE = ask for TPCrefit for the 3 daughter tracks
        Bool_t          fkQualityCutnTPCcls;            // Boolean : kTRUE = ask for at least n TPC clusters for each daughter track
        Bool_t          fkQualityCutPileup;             // Boolean : kTRUE = ask for no pileup events
        Bool_t          fwithSDD;                       // Boolean : kTRUE = Select the events that has and use the info from the SDD
        Int_t           fMinnTPCcls;                    // minimum number of TPC cluster for daughter tracks
        Float_t         fCentrLowLim;                   // Lower limit for centrality percentile selection
        Float_t         fCentrUpLim;                    // Upper limit for centrality percentile selection
        TString         fCentrEstimator;                // string for the centrality estimator
        Bool_t          fkUseCleaning;                  // Boolean : kTRUE = uses all the cleaning criteria of centrality selections (vertex cut + outliers) otherwise only outliers
        Float_t         fVtxRange;                      // to select events with |zvtx|<fVtxRange cm
        Float_t         fMinPtCutOnDaughterTracks;      // minimum pt cut on daughter tracks
        Float_t         fEtaCutOnDaughterTracks;        // pseudorapidity cut on daughter tracks
       

        
	AliCFContainer  *fCFContCascadeCuts;            //! Container meant to store all the relevant distributions corresponding to the cut variables
	
	

  AliQAProdMultistrange(const AliQAProdMultistrange&);            // not implemented
  AliQAProdMultistrange& operator=(const AliQAProdMultistrange&); // not implemented
  
  ClassDef(AliQAProdMultistrange, 7);
};

#endif
