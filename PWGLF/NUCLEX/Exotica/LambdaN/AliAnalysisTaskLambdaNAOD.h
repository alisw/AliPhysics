#ifndef ALIANALYSISTASKLAMBDANAOD_H
#define ALIANALYSISTASKLAMBDANAOD_H


/**************************************************************************
 * Author : Nicole Alice Martin (nicole.alice.martin@cern.ch)                  *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
//               AliAnalysisTaskLambdaNAOD  class
//  task for the investigation of (anti-)lambda-n bound state
//          uses the V0 finders, based on AODs or ESDS
//-----------------------------------------------------------------

class TF1;
class TH1F;
class TH2F;
class TH3F;
class THnSparse;

class TTree;

class AliAODInputHandler;
class AliAODEvent;
class AliAODVertex;
class AliAODv0;

class AliESDtrackCuts;

class AliPIDResponse;
class AliESDpid;

class AliESDInputHandler;
class AliESDEvent;
class AliESDVertex;
class AliESDv0;

#include <fstream>
#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
#include "AliStack.h"
#include "AliVTrack.h"

#define maxNofTracks 100000

class AliAnalysisTaskLambdaNAOD : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskLambdaNAOD();
  AliAnalysisTaskLambdaNAOD(const char *name);
  virtual ~AliAnalysisTaskLambdaNAOD();// {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);
  Int_t  Initialize();
  Int_t  SetupEvent();
  void   ResetEvent();
  void SetAnalysisType               (const char* analysisType          = "ESD") { fAnalysisType                = analysisType;               }


  /** Types of pdgCode */
  enum PdgCodeType_t { 
    kPDGPionPlus,
    kPDGPionMinus,
    kPDGProton,
    kPDGAntiProton,
    kPDGDeuteron,
    kPDGAntiDeuteron,
    kPDGHelium3,
    kPDGAntiHelium3,
    kPDGLambda,
    kPDGAntiLambda,
    kPDGLambdaNeutron,
    kPDGAntiLambdaNeutron,
    kPDGHypertriton,
    kPDGAntiHypertriton,
    kPdgCodeMax  // Number of enum entries
  };

  /** Types of Mass */
  enum MassType_t { 
    kMassPion,
    kMassProton,
    kMassDeuteron,
    kMassTriton,
    kMassHelium3,
    kMassMax  // Number of enum entries
  };

  /** Array of types of pdgCodes */
  static const Int_t fgkPdgCode[];           //! transient

  /** Array of types of pdgCodes */
  static const Double_t fgkMass[];        //! transient

 private:


  TString         fAnalysisType;                  // "ESD" or "AOD" analysis type	

  AliInputEventHandler *fEventHandler;            //for ESDs or AODs

  AliESDtrackCuts    *fESDtrackCutsV0;                            // basic cut variables for v0's
  AliESDpid          *fESDpid;                                  // basic TPC object for n-sigma cuts
  AliPIDResponse     *fPIDResponse;                   //! PID response object

  //
   
  
  TTree             *fTreeV0;                          //!
 
  TH1F              *fHistCentralityClass10;                   //! histo to look at the centrality distribution
  TH1F              *fHistCentralityPercentile;                //! histo to look at the centrality distribution
  TH1F              *fHistTriggerStat;                         //! Trigger statistics
  TH1F              *fHistTriggerStatAfterEventSelection;      //! Trigger statistics
  TH1F              *fHistLambdaNeutronPtGen;                  //! for MC 
  TH1F              *fHistLambdaNeutronPtGenMinBias;           //! for MC
  TH1F              *fHistLambdaNeutronPtGenCentral;           //! for MC
  TH1F              *fHistLambdaNeutronPtGenSemiCentral;       //! for MC
  TH1F              *fHistAntiLambdaNeutronPtGen;              //! for MC                                                                                                           
  TH1F              *fHistAntiLambdaNeutronPtGenMinBias;       //! for MC
  TH1F              *fHistAntiLambdaNeutronPtGenCentral;       //! for MC
  TH1F              *fHistAntiLambdaNeutronPtGenSemiCentral;   //! for MC                                                                                                           
  TH1F              *fHistLambdaNeutronInvaMassGen;            //! for MC
  TH1F              *fHistAntiLambdaNeutronInvaMassGen;        //! for MC
  TH1F              *fHistLambdaNeutronDecayLengthGen;         //! for MC
  TH1F              *fHistAntiLambdaNeutronDecayLengthGen;     //! for MC
  TH1F              *fHistLambdaNeutronPtAso;                  //! for MC
  TH1F              *fHistLambdaNeutronPtAsoCuts;              //! for MC
  TH1F              *fHistAntiLambdaNeutronPtAso;              //! for MC
  TH1F              *fHistAntiLambdaNeutronPtAsoCuts;          //! for MC
  TH1F              *fHistLambdaNeutronInvaMassAso;            //! for MC
  TH1F              *fHistAntiLambdaNeutronInvaMassAso;        //! for MC
  TH1F              *fHistLambdaNeutronDecayLengthAso;         //! for MC
  TH1F              *fHistAntiLambdaNeutronDecayLengthAso;     //! for MC

  TH2F              *fof;                                      //! debug histo for OnTheFlyStatus
  TH2F              *fHistArmenterosPodolanskiDeuteronPion;    //!
  TH2F              *fHistArmenterosPodolanskiAntiDeuteronPion;//!

  TH3F              *fHistDeDxQA;                              //! histo for a QA dE/dx

  Int_t              fNTriggers;                               //! N Triggers used

  Bool_t             fMCtrue;                                  //! flag if real data or MC is processed
  Bool_t             fOnlyQA;                                  //! flag if only QA histograms should be filled
  Bool_t             fTriggerFired[5];                         //! TriggerFired 0: MB | 1: CE | 2: SC | 3: EJE | 4: EGA



  //Tree variables
  //AliAODv0 *fV0object;                                         //! Tree variable
  Int_t  fItrk;                                                //! Tree variable

  Int_t fV0finder[maxNofTracks];                               //! Tree variable
  Int_t fkMB[maxNofTracks];                                    //! Tree variable
  Int_t fkCentral[maxNofTracks];                               //! Tree variable
  Int_t fkSemiCentral[maxNofTracks];                           //! Tree variable
  Int_t fkEMCEJE[maxNofTracks];                                //! Tree variable
  Int_t fkEMCEGA[maxNofTracks];                                //! Tree variable 

  Int_t fCentralityClass10[maxNofTracks];                      //! Tree variable 
  Int_t fCentralityPercentile[maxNofTracks];                   //! Tree variable
  Int_t fMultiplicity[maxNofTracks];                           //! Tree variable
  Int_t fRefMultiplicity[maxNofTracks];                        //! Tree variable

  Double_t fPtotN[maxNofTracks];                               //! Tree variable
  Double_t fPtotP[maxNofTracks];                               //! Tree variable
  Double_t fMotherPt[maxNofTracks];                            //! Tree variable

  Double_t fdEdxN[maxNofTracks];                               //! Tree variable
  Double_t fdEdxP[maxNofTracks];                               //! Tree variable
  Double_t fSignN[maxNofTracks];                               //! Tree variable
  Double_t fSignP[maxNofTracks];                               //! Tree variable

  Double_t fSigmadEdxPionPos[maxNofTracks];                       //! Tree variable
  Double_t fSigmadEdxPionNeg[maxNofTracks];                       //! Tree variable

  Float_t fDCAv0[maxNofTracks];                                //! Tree variable
  Float_t fCosinePAv0[maxNofTracks];                           //! Tree variable
  Float_t fDecayRadiusTree[maxNofTracks];                      //! Tree variable
  Double_t fInvaMassDeuteronPionTree[maxNofTracks];            //! Tree variable
  Int_t fChargeComboDeuteronPionTree[maxNofTracks];            //! Tree variable
  Bool_t fIsCorrectlyAssociated[maxNofTracks];                 //! Tree variable

  Double_t fAmenterosAlphaTree[maxNofTracks];                  //! Tree variable
  Double_t fAmenterosQtTree[maxNofTracks];                     //! Tree variable
  Int_t fRotationTree[maxNofTracks];                           //! Tree variable

  Double_t fImpactParameterDeuteronPos[maxNofTracks];          //! Tree variablen
  Double_t fImpactParameterDeuteronNeg[maxNofTracks];          //! Tree variablen                                                                                     
  Double_t fImpactParameterPionPos[maxNofTracks];              //! Tree variable
  Double_t fImpactParameterPionNeg[maxNofTracks];              //! Tree variable                                                                                      

  Double_t fImpactParameterDeuteronPosAliKF[maxNofTracks];     //! Tree variablen                                                              
  Double_t fImpactParameterDeuteronNegAliKF[maxNofTracks];    //! Tree variablen                                                                
  Double_t fImpactParameterPionPosAliKF[maxNofTracks];        //! Tree variable                                                                          
  Double_t fImpactParameterPionNegAliKF[maxNofTracks];        //! Tree variable      

  UShort_t fMinNClustersTPCPos[maxNofTracks];                     //! Tree variable
  UShort_t fMinNClustersTPCNeg[maxNofTracks];                     //! Tree variable

  Double_t fMaxChi2PerClusterTPCPos[maxNofTracks];              //! Tree variable
  Double_t fMaxChi2PerClusterTPCNeg[maxNofTracks];              //! Tree variable
  
  TObjArray         *fOutputContainer;         //! output data container for the histogramms
  //
  void  BinLogAxis(const THnSparse *h, Int_t axisNumber); //define function for log axis for search for Anti-Alpha candidates 
  //
  
  /** Check if event is triggred */
  Bool_t   IsTriggered();
  Bool_t   DeuteronPID(AliVTrack *trackP, AliVTrack *trackN, Double_t ptotP, Double_t ptotN, Int_t runNumber, Bool_t isDeuteron[3]);
  Bool_t   PionPID(AliVTrack *trackP, AliVTrack *trackN, Double_t ptotP, Double_t ptotN, Int_t runNumber, Bool_t isPion[2]);
  Bool_t   TrackCuts(AliVTrack *track, Bool_t testTrackCuts);
  //Bool_t   FilterBit(AliVTrack *track, Bool_t testFilterBit);
  Double_t MomentumInnerParam(AliVTrack *track, Double_t ptot);
  UShort_t TPCclusters(AliVTrack *track, UShort_t numberOfTPCclusters);
  Double_t TPCchi2(AliVTrack *track, Double_t numberOfChi2clustersTPC, UShort_t numberOfTPCclusters);
  Double_t ImpactParameter(AliVTrack *track, Double_t dcaToVertex);

 
  void               MCGenerated(AliStack* stack, Int_t multiplicity);             //! function to loop over the genrated particles
  
  void               MCTwoBodyDecay (AliStack* stack, const TParticle *tparticleMother, Long_t PDGMother, Long_t PDGFirstDaughter, Long_t PDGSecondDaughter, Double_t massFirstDaughter, Double_t massSecondDaughter, Int_t multiplicity);                         //! function to calculate the invariant mass of two daughters and the pt of the mother
  //void RotateKFParticle(AliKFParticle * kfParticle,Double_t angle, const AliVEvent * const ev);


  AliAnalysisTaskLambdaNAOD(const AliAnalysisTaskLambdaNAOD&); // not implemented
  AliAnalysisTaskLambdaNAOD& operator=(const AliAnalysisTaskLambdaNAOD&); // not implemented
  //
  ClassDef(AliAnalysisTaskLambdaNAOD, 1); // example of analysis
};

#endif
