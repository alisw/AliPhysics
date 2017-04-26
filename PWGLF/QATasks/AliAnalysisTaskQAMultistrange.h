#ifndef ALIANALYSISTASKQAMULTISTRANGE_H
#define ALIANALYSISTASKQAMULTISTRANGE_H

/*  See cxx source for full Copyright notice */

/////////////////////////////////////////////////////////////
//
//            AliAnalysisTaskQAMultistrange class
//              Origin AliAnalysisTaskCheckCascade
//              This task has four roles :
//                1. QAing the Cascades from ESD and AOD
//                   Origin:  AliAnalysisTaskESDCheckV0 by Boris Hippolyte Nov2007, hippolyt@in2p3.fr
//                2. Prepare the plots which stand as raw material for yield extraction (wi/wo PID)
//                3. Supply an AliCFContainer meant to define the optimised topological selections
//                Adapted to Cascade : A.Maire Mar2008, antonin.maire@ires.in2p3.fr
//                Modified :           A.Maire Mar2010, antonin.maire@ires.in2p3.fr
//                Modified for PbPb analysis: M. Nicassio Feb 2011, maria.nicassio@ba.infn.it
//                Modified for QA production: D. Colella 2013, domenico.colella@cern.ch
//
/////////////////////////////////////////////////////////////

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

class AliAnalysisTaskQAMultistrange : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskQAMultistrange();
  AliAnalysisTaskQAMultistrange(const char *name);
  virtual ~AliAnalysisTaskQAMultistrange();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetIsMC                       (Bool_t isMC                       = kFALSE) { fisMC                        = isMC;                       } 
  void SetAnalysisType               (const char* analysisType          = "ESD" ) { fAnalysisType                = analysisType;               }
  void SetQualityCutTPCrefit         (Bool_t qualityCutTPCrefit         = kTRUE ) { fkQualityCutTPCrefit         = qualityCutTPCrefit;         }
  void SetQualityCutnTPCcls          (Bool_t qualityCutnTPCcls          = kTRUE ) { fkQualityCutnTPCcls          = qualityCutnTPCcls;          }
  void SetQualityCutMinnTPCcls       (Int_t  minnTPCcls                 = 70    ) { fMinnTPCcls                  = minnTPCcls;                 }
  void SetMinptCutOnDaughterTracks   (Float_t minptdaughtrks            = 0.    ) { fMinPtCutOnDaughterTracks    = minptdaughtrks;             }

 private:
        // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
        // your data member object is created on the worker nodes and streaming is not needed.
        // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14


        Bool_t          fisMC;                          // Boolean : kTRUE = is a MC production
        TString         fAnalysisType;                  // "ESD" or "AOD" analysis type	
        AliPIDResponse *fPIDResponse;                   // PID response object
        Bool_t          fkQualityCutTPCrefit;           // Boolean : kTRUE = ask for TPCrefit for the 3 daughter tracks
        Bool_t          fkQualityCutnTPCcls;            // Boolean : kTRUE = ask for at least n TPC clusters for each daughter track
        Int_t           fMinnTPCcls;                    // minimum number of TPC cluster for daughter tracks
        Float_t         fMinPtCutOnDaughterTracks;      // minimum pt cut on daughter tracks
        Float_t         fEtaCutOnDaughterTracks;        // pseudorapidity cut on daughter tracks
       
        TList  *fListHistMultistrangeQA;                //! List of Cascade histograms
        TH1F *fHistEventSel;                            //! Gives the number of the events after each event selection
        TH1F *fHistMassXiMinus;                         //! Gives the invariant mass distr. for the Xi minus without PID info on daughter tracks
        TH1F *fHistMassXiPlus;                          //! Gives the invariant mass distr. for the Xi plus without PID info on daughter tracks
        TH1F *fHistMassOmegaMinus;                      //! Gives the invariant mass distr. for the Omega minus without PID info on daughter tracks
        TH1F *fHistMassOmegaPlus;                       //! Gives the invariant mass distr. for the Omega plus without PID info on daughter tracks
        TH1F *fHistCascadeMultiplicityXiMinus;          //! Gives the distribution of the number of Xi minus per event
        TH1F *fHistCascadeMultiplicityXiPlus;           //! Gives the distribution of the number of Xi plus per event
        TH1F *fHistCascadeMultiplicityOmegaMinus;       //! Gives the distribution of the number of Omega minus per event
        TH1F *fHistCascadeMultiplicityOmegaPlus;        //! Gives the distribution of the number of Omega plus per event  
	AliCFContainer  *fCFContCascadeCuts;            //! Container meant to store all the relevant distributions corresponding to the cut variables
        AliCFContainer  *fCFContCascadeMCCuts;          //! Container meant to store all the relevant distributions corresponding to the cut variables with Association for MC productions
        AliCFContainer  *fCFContCascadeMCgen;           //! Container meant to store general variables for MC generated particles 
	

  AliAnalysisTaskQAMultistrange(const AliAnalysisTaskQAMultistrange&);            // not implemented
  AliAnalysisTaskQAMultistrange& operator=(const AliAnalysisTaskQAMultistrange&); // not implemented
  
  ClassDef(AliAnalysisTaskQAMultistrange, 10);
};

#endif
