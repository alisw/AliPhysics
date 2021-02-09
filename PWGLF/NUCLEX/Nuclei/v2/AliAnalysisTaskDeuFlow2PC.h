#ifndef ALIANALYSISTASKDEUFLOW2PC_H
#define ALIANALYSISTASKDEUFLOW2PC_H
//NUOVO
// Author:Ramona Lea based on M. Nicassio m.nicassio@gsi.de maria.nicassio@cern.ch

class TH3F;

class AliESDEvent;
class AliESDtrackCuts;
class AliPIDResponse;
class THnSparse;
#include "AliAnalysisTaskSE.h"
#include "AliAnalysishDEventCollection.h"
#include "AliEventCuts.h"

/* #include "AliCFContainer.h"*/

class AliAnalysisTaskDeuFlow2PC : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskDeuFlow2PC();
  AliAnalysisTaskDeuFlow2PC(const char *name);
  virtual ~AliAnalysisTaskDeuFlow2PC();// {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);

  enum firstpart_t {kElectron,kMuon,kPion,kKaon,kProton,kDeuteron,kAny};

  void SetAnalysisType (const char* analysisType = "ESD")   { fAnalysisType = analysisType; }
  void SetUseStandardCuts (Bool_t usestandardesdtrackcuts) {fUseStandardCuts = usestandardesdtrackcuts;}
  void SetCentrality   (Float_t lowlimcent = 0., Float_t uplimcent = 90.) { fCentrLowLim = lowlimcent;  fCentrUpLim = uplimcent; }
  void SetReadMCTruth (Bool_t readmctruth) {fReadMCTruth = readmctruth;}
  void SetUseContainer (Bool_t kusecontainer) { fUseContainer = kusecontainer;}
  void SetCollidingSystem (const char* collSystem = "pp") { fCollidingSystem = collSystem ;}
  void SetYear (Int_t year = 2010) { fYear = year;}
  void SetHMTrigger (Bool_t isHMtrigger = kFALSE) { fHMtrigger = isHMtrigger;}
  void SetFilterBit (Int_t filterBit = 128) { fFilterBit = filterBit;}

  // void DoPairshh    (const Float_t centralityBin, int fieldsign);
  // void DoPairsh1h2    (const Float_t centralityBin, int fieldsign,Float_t );
  void DoPairshh      (const Float_t lcentrality, int fieldsign, const Double_t fSphericityvalue);
  void DoPairsh1h2    (const Float_t lcentrality, int fieldsign, const Double_t fSphericityvalue);

  double CalculateKstar(double momentum1[3], double momentum2[3], double mass1, double mass2); 
  double CalculateMass(double momentum1[3], double momentum2[3], double mass1, double mass2); 
  void SetSftPosR125(AliVTrack *track, const Float_t bfield, Double_t priVtx[3], Double_t fXSftR125[3] ); 
  double CalculateDphiSatR12m(Double_t pos1SftR125[3], Double_t pos2SftR125[3]);
  double CalculateDphiSatR12m(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2);
  double CalculateDPhiStar(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2,Double_t rad);
  double CalculateDeltaEta(Double_t eta1, Double_t eta2); 
  double CalculateDeltaTheta(Double_t theta1, Double_t theta2); 
  double CalculateSphericityofEvent(AliAODEvent *aodEvent);
  bool IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP);
  bool IsPionNSigma(double mom, float nsigmaTPCPi, float nsigmaTOFPi);

  void SetFirstParticle(firstpart_t firstpart)   {fFirstpart = firstpart;}
  void SetSecondParticle(firstpart_t secondpart) {fSecondpart = secondpart;}
  void SetnSigmaTPCPIDfirstParticle(Float_t nsigma) {fnSigmaTPCPIDfirstParticle = nsigma;}
  void SetnSigmaTPCTOFPIDfirstParticle (Float_t nsigma) { fnSigmaTPCTOFPIDfirstParticle = nsigma;}
  void SetnSigmaTPCPIDsecondParticle(Float_t nsigma) {fnSigmaTPCPIDsecondParticle = nsigma;}
  void SetnSigmaTPCTOFPIDsecondParticle (Float_t nsigma) { fnSigmaTPCTOFPIDsecondParticle = nsigma;}
  void SetTrackBufferSize(UShort_t trackbuffersize) { fTrackBufferSize = trackbuffersize;} 
  void SetFirstPartMaxMult(Int_t firstpartmult)    { fMaxSecondMult  = firstpartmult;}
  void SetSecondPartMaxMult(Int_t secpartmult)          { fMaxSecondMult  = secpartmult;}
  void SetApplyTtc(Bool_t kapplyttc) { fkApplyTtc = kapplyttc;}
  void SetDphisMin(Float_t dphismin) { fDphisMin = dphismin;}
  void SetDetasMin(Float_t detasmin) { fDetasMin = detasmin;}
  void SetNEventsToMix(short nevmixing) { fnEventsToMix = nevmixing;}
  void SetMomentumLimitForTOFPIDfirst(Float_t momemtumlimitforTOFPID) { fMomemtumLimitForTOFPIDfirst = momemtumlimitforTOFPID;} 
  void SetMomentumLimitForTOFPIDsecond(Float_t momemtumlimitforTOFPID) { fMomemtumLimitForTOFPIDsecond = momemtumlimitforTOFPID;} 
  void SetApplyRatioCrRnFindCut(Bool_t kapplycrrowsnfindcut) { fkApplyRatioCrRnFindCut = kapplycrrowsnfindcut;}
  void SetCutOnTPCIP(Bool_t kcutontpcip) { fkCutOnTPCIP = kcutontpcip;}
  void SetIPCutxyPrim(Float_t ipcutxy) { fIPCutxyPrim = ipcutxy;}
  void SetIPCutzPrim (Float_t ipcutz)  { fIPCutzPrim = ipcutz;}
  void SetIPCutxySec(Float_t ipcutxy) { fIPCutxySec = ipcutxy;}
  void SetIPCutzSec (Float_t ipcutz)  { fIPCutzSec = ipcutz;}
  void SetMinPtPrim(Float_t minptforprim) { fMinPtForPrim  = minptforprim;}  
  void SetMaxPtPrim(Float_t maxptforprim) { fMaxPtForPrim  = maxptforprim;}
  void SetMinPtSec(Float_t minptforsec) { fMinPtForSec  = minptforsec;}  
  void SetMaxPtSec(Float_t maxptforsec) { fMaxPtForSec  = maxptforsec;}
  void SetPropagateGlobal(Bool_t propagateglobal) { fkPropagateGlobal = propagateglobal;}
  void DoSphirocity(Bool_t doSphericity) {fkDoSphericity = doSphericity;}
  void SetRadius (Float_t radius)  { fRadius = radius;}
  Double_t EtaS( Double_t posSftR125[3] ) const; 
  Double_t ThetaS( Double_t posSftR125[3] ) const;
  AliEventCuts fEventCuts;
 private:
 
  AliESDEvent *fESDevent;                         //! 
  AliAODEvent *fAODevent;                         //! 


  TString     fAnalysisType;                  // "ESD" or "AOD" analysis type  
  TString     fCollidingSystem;               // "pp", "pPb", "PbPb" 
  Int_t       fYear;
  Bool_t      fHMtrigger;
  Int_t       fFilterBit;
  firstpart_t fFirstpart;
  firstpart_t fSecondpart;
  
  Float_t     fnSigmaTPCPIDfirstParticle;
  Float_t     fnSigmaTPCTOFPIDfirstParticle;
  Float_t     fnSigmaTPCPIDsecondParticle;
  Float_t     fnSigmaTPCTOFPIDsecondParticle;
  
  Bool_t      fReadMCTruth;                            // if read MC truth
  Bool_t      fUseContainer;
  Bool_t      fUseStandardCuts;                        // if to use standard ESD track cuts or user-defined ones 
  Bool_t      fkApplyTtc;

  Float_t     fDphisMin;
  Float_t     fDetasMin; 

  Float_t     fMomemtumLimitForTOFPIDfirst;
  Float_t     fMomemtumLimitForTOFPIDsecond;
  
  Bool_t      fkApplyRatioCrRnFindCut;
  Bool_t      fkCutOnTPCIP;
  
  Float_t     fIPCutxyPrim;
  Float_t     fIPCutzPrim;
  Float_t     fIPCutxySec;
  Float_t     fIPCutzSec;
  
  Float_t     fMinPtForPrim;
  Float_t     fMaxPtForPrim;
  Float_t     fMinPtForSec;
  Float_t     fMaxPtForSec;
  	      
  Float_t     fRadius;
	      
  Bool_t      fkDoSphericity;
  Bool_t      fkPropagateGlobal;

  AliESDtrackCuts    *fESDtrackCuts;              //! basic cut variables for tracks added ! not sure
  AliPIDResponse     *fPIDResponse;               //! PID response object

  Float_t fCentrLowLim;                           // centrality settings: lower limit
  Float_t fCentrUpLim;                            // upper limit

  // Store pointers to global tracks for pid and dca

  Int_t     *farrGT;
  UShort_t  fTrackBufferSize;          // Size fo the above array, ~12000 for PbPb

  //  AliAODTrack     **farrGT;                  //! Array of pointers, just nicely sorted according to the id
  AliAnalysishDEventCollection ***fEventColl; //!
  AliAnalysishDEventCollection ****fEventCollwSp; //!
  AliAnalysishDEvent *fEvt;                                //!

  Int_t     fMaxFirstMult;
  Int_t     fMaxSecondMult;

  Float_t     fPDGp;
  Float_t     fPDGk;
  Float_t     fPDGd;
  Float_t     fPDGpi;

  Float_t     fPDGfirst;
  Float_t     fPDGsecond;

  Float_t     fPDGCodefirst;
  Float_t     fPDGCodesecond;

  int fzVertexBins;
  int fnCentBins;
  int fnSphericityBins;
  
  short fnEventsToMix;

  TH1F  *fHistEventMultiplicity;                //! Counts Event multiplicity
  TH1I  *hmult;                                 //! Counts tracklets
  TH1F  *fHistCentrality;                       //! Counts events in centrality
  TH1F  *fHistVertexDistribution;               //! Vertex distribution
  TH1F  *fHistSphericity;                       //! Sphericity distibution
  TH1F  *fHistMultiplicityOfMixedEvent;         //!

  TH2F  *fHistTriggptvsCentrality;              //!

  TH2F  *fHistTPCdEdx;                          //!   
  TH2F  *fHistFirstTPCdEdx;                     //!
  TH2F  *fHistSecondTPCdEdx;                    //!
  
  TH1F  *fHistnTPCCrossedRFirst;                //!
  TH1F  *fHistRationTPCCrossedRnFindFirst;      //!
  TH1F  *fHistSharedFrTPCclFirst;               //!
  
  TH1F  *fHistnTPCCrossedRSecond;               //!
  TH1F  *fHistRationTPCCrossedRnFindSecond;     //!
  TH1F  *fHistSharedFrTPCclSecond;              //!

  TH3F  *fHistyptFirst;                         //!
  TH3F  *fHistyptSecond;                        //!
 
  TH2F  *fHistphietaFirst;                      //!
  TH2F  *fHistphietaSecond;                     //!
  
  TH2F  *fHistIPtoPVxyzTPCFirst;                //!
  TH2F  *fHistIPtoPVxyzGlobalFirst;             //!

  TH2F  *fHistIPtoPVxyzTPCSecond;               //!
  TH2F  *fHistIPtoPVxyzGlobalSecond;            //!

  TH2F  *fHistFirstTOFmisvspt;                  //!
  TH2F  *fHistFirstTOFmisvsp;                   //!
  TH2F  *fHistFirstTOFnsigmavspt;               //!
  TH2F  *fHistFirstTOFnsigmavsp;                //!
  TH2F  *fHistFirstTOFsignalvsp;                //!
  TH2F  *fHistFirstTOFsignalvspt;               //!

  TH2F  *fHistFirstTOFTPCsignalvspt;            //!
  TH2F  *fHistFirstMultvsCent;                  //!

  TH2F  *fHistSecondTOFmisvspt;                 //!
  TH2F  *fHistSecondTOFmisvsp;                  //!
  TH2F  *fHistSecondTOFnsigmavspt;              //!
  TH2F  *fHistSecondTOFnsigmavsp;               //!
  TH2F  *fHistSecondTOFsignalvsp;               //!
  TH2F  *fHistSecondTOFsignalvspt;              //!

  TH2F  *fHistSecondTOFTPCsignalvspt;           //!
  TH2F  *fHistSecondMultvsCent;                 //!
 
  TH2F  *fHistFirstTPCdEdxAfter;                //!
  TH2F  *fHistSecondTPCdEdxAfter;               //! 
  TH2F  *fHistFirstTOFTPCsignalvsptAfter;       //! 
  TH2F  *fHistSecondTOFTPCsignalvsptAfter;      //! 
  
  TH2F *fHistFirstMassTOFvsPt3sTPC;             //!  
  TH2F *fHistSecondMassTOFvsPt3sTPC;            //!
  
  TH2F *fHistFirstMassTOFvsPt3sTPC3sTOF;        //!
  TH2F *fHistSecondMassTOFvsPt3sTPC3sTOF;       //! 
   
  //----------------
  
  TTree *fHistSparseSignal ;  
  TTree *fHistSparseBkg;      

  Int_t   tSignP1;
  Int_t   tSignP2;
  Float_t tCentrality;
  Float_t tDCAxyP1;
  Float_t tDCAzP1; 
  Float_t tDCAxyP2; 
  Float_t tDCAzP2; 
  Float_t tKtpair;
  Float_t tkStar;
  Float_t tptP1;
  Float_t tptP2;
  Float_t tDEta;
  Float_t tDPhiStar;
  Float_t tDPhi;
  Float_t tMassPair;
  Float_t tSphericity;
  Float_t tMassS;
  Float_t tDTheta;
  Int_t   tMCtruepair;
  Int_t   tMCSameMother;
  Int_t   tMCMotherP1;
  Int_t   tMCMotherP2;
  Int_t   tMCptcTypeP1;
  Int_t   tMCptcTypeP2;
  Int_t   tMCSameGM;
  Int_t   tMotherPDG;
  Int_t   tpdgcodeP1; 
  Int_t   tpdgcodeP2;
  Float_t tKstarGen;

  //-----------------
  
  TH1F              *fHistTrackBufferOverflow;             //! 
  TList             *fOutputContainer;                     //! output list for the histograms

  //
  AliAnalysisTaskDeuFlow2PC(const AliAnalysisTaskDeuFlow2PC&); // not implemented
  AliAnalysisTaskDeuFlow2PC& operator=(const AliAnalysisTaskDeuFlow2PC&); // not implemented
  //
  ClassDef(AliAnalysisTaskDeuFlow2PC, 2);
};

#endif
                 
