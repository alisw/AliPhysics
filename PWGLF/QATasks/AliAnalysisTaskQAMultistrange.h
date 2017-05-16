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
        TH1F *fHistCascadeMultiplicityXiMinus;          //! Gives the distribution of the number of Xi minus per event
        TH1F *fHistCascadeMultiplicityXiPlus;           //! Gives the distribution of the number of Xi plus per event
        TH1F *fHistCascadeMultiplicityOmegaMinus;       //! Gives the distribution of the number of Omega minus per event
        TH1F *fHistCascadeMultiplicityOmegaPlus;        //! Gives the distribution of the number of Omega plus per event  

        TH2F *fHistVarDcaCascDaughtXiMinus;               //!
        TH2F *fHistVarDcaBachToPrimVertexXiMinus;         //!
        TH2F *fHistVarCascCosineOfPointingAngleXiMinus;   //!
        TH2F *fHistVarCascRadiusXiMinus;                  //!
        TH2F *fHistVarInvMassLambdaAsCascDghterXiMinus;   //!
        TH2F *fHistVarDcaV0DaughtersXiMinus;              //!
        TH2F *fHistVarV0CosineOfPAToCascVertexXiMinus;    //!
        TH2F *fHistVarV0RadiusXiMinus;                    //!
        TH2F *fHistVarDcaV0ToPrimVertexXiMinus;           //!
        TH2F *fHistVarDcaPosToPrimVertexXiMinus;          //!
        TH2F *fHistVarDcaNegToPrimVertexXiMinus;          //!
        TH1F *fHistMassXiMinus;                           //!
        TH1F *fHistVarTransvMomentumXiMinus;              //!
        TH1F *fHistVarRapidityXiMinus;                    //!
        TH2F *fHistVarCascProperLengthXiMinus;            //!
        TH2F *fHistVarV0ProperLengthXiMinus;              //!
        TH1F *fHistGenVarTotMomXiMinus;                   //!
        TH1F *fHistGenVarTransvMomXiMinus;                //!
        TH1F *fHistGenVarYXiMinus;                        //! 
        TH1F *fHistGenVarEtaXiMinus;                      //!
        TH1F *fHistGenVarThetaXiMinus;                    //!
        TH1F *fHistGenVarPhiXiMinus;                      //!
 
        TH2F *fHistVarDcaCascDaughtXiPlus;                //!
        TH2F *fHistVarDcaBachToPrimVertexXiPlus;          //!
        TH2F *fHistVarCascCosineOfPointingAngleXiPlus;    //!
        TH2F *fHistVarCascRadiusXiPlus;                   //!
        TH2F *fHistVarInvMassLambdaAsCascDghterXiPlus;    //!
        TH2F *fHistVarDcaV0DaughtersXiPlus;               //! 
        TH2F *fHistVarV0CosineOfPAToCascVertexXiPlus;     //!
        TH2F *fHistVarV0RadiusXiPlus;                     //!
        TH2F *fHistVarDcaV0ToPrimVertexXiPlus;            //!
        TH2F *fHistVarDcaPosToPrimVertexXiPlus;           //!
        TH2F *fHistVarDcaNegToPrimVertexXiPlus;           //!
        TH1F *fHistMassXiPlus;                            //! 
        TH1F *fHistVarTransvMomentumXiPlus;               //!
        TH1F *fHistVarRapidityXiPlus;                     //!
        TH2F *fHistVarCascProperLengthXiPlus;             //!
        TH2F *fHistVarV0ProperLengthXiPlus;               //!
        TH1F *fHistGenVarTotMomXiPlus;                    //!
        TH1F *fHistGenVarTransvMomXiPlus;                 //!
        TH1F *fHistGenVarYXiPlus;                         //!
        TH1F *fHistGenVarEtaXiPlus;                       //!
        TH1F *fHistGenVarThetaXiPlus;                     //!
        TH1F *fHistGenVarPhiXiPlus;                       //!

        TH2F *fHistVarDcaCascDaughtOmegaMinus;            //!
        TH2F *fHistVarDcaBachToPrimVertexOmegaMinus;      //!
        TH2F *fHistVarCascCosineOfPointingAngleOmegaMinus;//!  
        TH2F *fHistVarCascRadiusOmegaMinus;               //!
        TH2F *fHistVarInvMassLambdaAsCascDghterOmegaMinus;//! 
        TH2F *fHistVarDcaV0DaughtersOmegaMinus;           //!
        TH2F *fHistVarV0CosineOfPAToCascVertexOmegaMinus; //!
        TH2F *fHistVarV0RadiusOmegaMinus;                 //!
        TH2F *fHistVarDcaV0ToPrimVertexOmegaMinus;        //!
        TH2F *fHistVarDcaPosToPrimVertexOmegaMinus;       //!
        TH2F *fHistVarDcaNegToPrimVertexOmegaMinus;       //!
        TH1F *fHistMassOmegaMinus;                        //!
        TH1F *fHistVarTransvMomentumOmegaMinus;           //!
        TH1F *fHistVarRapidityOmegaMinus;                 //!
        TH2F *fHistVarCascProperLengthOmegaMinus;         //!
        TH2F *fHistVarV0ProperLengthOmegaMinus;           //!
        TH1F *fHistGenVarTotMomOmegaMinus;                //!
        TH1F *fHistGenVarTransvMomOmegaMinus;             //!
        TH1F *fHistGenVarYOmegaMinus;                     //!
        TH1F *fHistGenVarEtaOmegaMinus;                   //!
        TH1F *fHistGenVarThetaOmegaMinus;                 //!
        TH1F *fHistGenVarPhiOmegaMinus;                   //!
 
        TH2F *fHistVarDcaCascDaughtOmegaPlus;             //!
        TH2F *fHistVarDcaBachToPrimVertexOmegaPlus;       //!
        TH2F *fHistVarCascCosineOfPointingAngleOmegaPlus; //! 
        TH2F *fHistVarCascRadiusOmegaPlus;                //!
        TH2F *fHistVarInvMassLambdaAsCascDghterOmegaPlus; //!
        TH2F *fHistVarDcaV0DaughtersOmegaPlus;            //!
        TH2F *fHistVarV0CosineOfPAToCascVertexOmegaPlus;  //! 
        TH2F *fHistVarV0RadiusOmegaPlus;                  //!
        TH2F *fHistVarDcaV0ToPrimVertexOmegaPlus;         //!
        TH2F *fHistVarDcaPosToPrimVertexOmegaPlus;        //!
        TH2F *fHistVarDcaNegToPrimVertexOmegaPlus;        //!
        TH1F *fHistMassOmegaPlus;                         //!
        TH1F *fHistVarTransvMomentumOmegaPlus;            //!
        TH1F *fHistVarRapidityOmegaPlus;                  //!
        TH2F *fHistVarCascProperLengthOmegaPlus;          //!
        TH2F *fHistVarV0ProperLengthOmegaPlus;            //!
        TH1F *fHistGenVarTotMomOmegaPlus;                 //!
        TH1F *fHistGenVarTransvMomOmegaPlus;              //!
        TH1F *fHistGenVarYOmegaPlus;                      //!
        TH1F *fHistGenVarEtaOmegaPlus;                    //!
        TH1F *fHistGenVarThetaOmegaPlus;                  //!
        TH1F *fHistGenVarPhiOmegaPlus;                    //! 


  AliAnalysisTaskQAMultistrange(const AliAnalysisTaskQAMultistrange&);            // not implemented
  AliAnalysisTaskQAMultistrange& operator=(const AliAnalysisTaskQAMultistrange&); // not implemented
  
  ClassDef(AliAnalysisTaskQAMultistrange, 12);
};

#endif
