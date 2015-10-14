#ifndef ALIANALYSISTASKHCASCADEFEMTO_H
#define ALIANALYSISTASKHCASCADEFEMTO_H

// Author: M. Nicassio m.nicassio@gsi.de maria.nicassio@cern.ch

class TH3F;

class AliESDEvent;
class AliESDtrackCuts;
class AliPIDResponse;

#include "AliAnalysisTaskSE.h"
#include "AliAnalysishCascadeEventCollection.h"
#include "AliCFContainer.h"

class AliAnalysisTaskhCascadeFemto : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskhCascadeFemto();
  AliAnalysisTaskhCascadeFemto(const char *name);
  virtual ~AliAnalysisTaskhCascadeFemto();// {}

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(const Option_t*);

  enum firstpart_t {kPion, kProton, kXi, kOmega};

  void SetAnalysisType (const char* analysisType = "ESD")   { fAnalysisType = analysisType; }
  void SetMassHist     (Float_t nbins = 100, Float_t lowlim = 2.028, Float_t uplim = 2.228) {fNbinsMass = nbins; fLowMassLim = lowlim; fUpMassLim = uplim; }
  void SetMassHistGrandmother (Float_t nbins = 100, Float_t lowlim = 2.028, Float_t uplim = 2.228) {fNbinsMassGm = nbins; fLowMassLimGm = lowlim; fUpMassLimGm = uplim;}
  void SetUseStandardCuts (Bool_t usestandardesdtrackcuts) {fUseStandardCuts = usestandardesdtrackcuts;}
  void SetCentrality   (Float_t lowlimcent = 0., Float_t uplimcent = 90.) { fCentrLowLim = lowlimcent;  fCentrUpLim = uplimcent; }
  void SetReadMCTruth (Bool_t readmctruth) {fReadMCTruth = readmctruth;}
  void SetUseContainer (Bool_t kusecontainer) { fUseContainer = kusecontainer;}

  void DoPairshCasc (const AliAnalysishCascadeEvent *event, const Float_t centralityBin); 
  void DoPairshh    (const AliAnalysishCascadeEvent *event, const Float_t centralityBin, int fieldsign);

  Int_t CheckDaughterTrack (Int_t xiIndex, Int_t daughterIndex, const AliAnalysishCascadeEvent *event, Int_t*  checkeddaughters);
  void SelectBestCandidate ( const AliAnalysishCascadeEvent *event);
  double CalculateKstar(double momentum1[3], double momentum2[3], double mass1, double mass2); 
  void ProtonOrigin();
  void SetSftPosR125(AliVTrack *track, const Float_t bfield, Double_t priVtx[3], Double_t fXSftR125[3] ); 
  double CalculateDphiSatR12m(Double_t pos1SftR125[3], Double_t pos2SftR125[3]);
  double CalculateDphiSatR12m(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2);
  void SetFirstParticle(firstpart_t firstpart) {fFirstpart = firstpart;}
  void SetSecondParticle(firstpart_t secondpart) {fSecondpart = secondpart;}
  void SetMassWindowCascades(Float_t masswincasc) { fMassWindowCascades = masswincasc;}
  void SetnSigmaTPCPIDfirstParticle(Float_t nsigma) { fnSigmaTPCPIDfirstParticle = nsigma;}
  void SetnSigmaTPCTOFPIDfirstParticle (Float_t nsigma) { fnSigmaTPCTOFPIDfirstParticle = nsigma;}
  void SetnSigmaTPCPIDsecondParticleDau (Float_t nsigma) { fnSigmaTPCPIDsecondParticleDau = nsigma;}
  void SetTrackBufferSize(UShort_t trackbuffersize) { fTrackBufferSize = trackbuffersize;} 
  void SetFirstPartMaxMult(Int_t firstpartmult)    { fMaxPMult = firstpartmult;}
  void SetSecondPartMaxMult(Int_t secpartmult)          { fMaxXiMult = secpartmult;}
  void SetApplyTtc(Bool_t kapplyttc) { fkApplyTtc = kapplyttc;}
  void SetDphisMin(Float_t dphismin) { fDphisMin = dphismin;}
  void SetDetasMin(Float_t detasmin) { fDetasMin = detasmin;}
  void SetCascadeSideBands(Bool_t kcascadesidebands) { fkCascadeSideBands = kcascadesidebands;} 
  void SetCascadeTightCuts(Bool_t kcascadetightcuts) { fkTightCutsForCascades = kcascadetightcuts;}
  void SetNEventsToMix(short nevmixing) { fnEventsToMix = nevmixing;}
  void SetMomentumLimitForTOFPID(Float_t momemtumlimitforTOFPID) { fMomemtumLimitForTOFPID = momemtumlimitforTOFPID;} 
  void SetApplyRatioCrRnFindCut(Bool_t kapplycrrowsnfindcut) { fkApplyRatioCrRnFindCut = kapplycrrowsnfindcut;}
  void SetCutOnTPCIP(Bool_t kcutontpcip) { fkCutOnTPCIP = kcutontpcip;}
  void SetIPCutxy(Float_t ipcutxy) { fIPCutxy = ipcutxy;}
  void SetIPCutz (Float_t ipcutz)  { fIPCutz = ipcutz;}
  void SetMinPtPrim(Float_t minptforprim) { fMinPtForPrim  = minptforprim;}  
  void SetMaxPtPrim(Float_t maxptforprim) { fMaxPtForPrim  = maxptforprim;}
  void SetMinPtCasc(Float_t minptforcasc) { fMinPtForCasc  = minptforcasc;}
  void SetMaxPtCasc(Float_t maxptforcasc) { fMaxPtForCasc  = maxptforcasc;}
  void SetIPCutBac(Float_t ipcutbac) { fIPCutBac = ipcutbac;}
  void SetApplyYcutCasc(Bool_t applyycutcasc) { fkApplyYcutCasc = applyycutcasc;} 
  void SetPropagateGlobal(Bool_t propagateglobal) { fkPropagateGlobal = propagateglobal;}

 
  Double_t EtaS( Double_t posSftR125[3] ) const; 
  Double_t ThetaS( Double_t posSftR125[3] ) const;

 private:
 
  AliESDEvent *fESDevent;                         //! 
  AliAODEvent *fAODevent;                         //! 


  TString         fAnalysisType;                  // "ESD" or "AOD" analysis type       
  firstpart_t fFirstpart;
  firstpart_t fSecondpart;
  Float_t fMassWindowCascades;
  Float_t fnSigmaTPCPIDfirstParticle;
  Float_t fnSigmaTPCTOFPIDfirstParticle;
  Float_t fnSigmaTPCPIDsecondParticleDau;
  Bool_t fkCascadeSideBands;
  Bool_t fkTightCutsForCascades;
  Bool_t fReadMCTruth;                            // if read MC truth
  Bool_t fUseContainer;
  Bool_t fUseStandardCuts;                        // if to use standard ESD track cuts or user-defined ones 
  Bool_t fkApplyTtc;
  Float_t fDphisMin;
  Float_t fDetasMin; 

  Float_t fMomemtumLimitForTOFPID;
  Bool_t fkApplyRatioCrRnFindCut;
  Bool_t fkCutOnTPCIP;
  Float_t fIPCutxy;
  Float_t fIPCutz;
  Float_t fMinPtForPrim;
  Float_t fMaxPtForPrim;
  Float_t fMinPtForCasc;
  Float_t fMaxPtForCasc;
  Float_t fIPCutBac;
  Bool_t fkApplyYcutCasc;
  Bool_t fkPropagateGlobal;

  AliESDtrackCuts    *fESDtrackCuts;              //! basic cut variables for tracks added ! not sure
  AliPIDResponse     *fPIDResponse;               //! PID response object

  Float_t fCentrLowLim;                           // centrality settings: lower limit
  Float_t fCentrUpLim;                            // upper limit

  Float_t fNbinsMass;                             // Inv mass histo settings for a two-body decay
  Float_t fLowMassLim;
  Float_t fUpMassLim;

  Float_t fNbinsMassGm;                           // Inv mass histo settings for a casade-like decay 
  Float_t fLowMassLimGm;
  Float_t fUpMassLimGm;

  // Store pointers to global tracks for pid and dca

  Int_t     *farrGT;
  UShort_t  fTrackBufferSize;          // Size fo the above array, ~12000 for PbPb

//  AliAODTrack     **farrGT;                  //! Array of pointers, just nicely sorted according to the id
  AliAnalysishCascadeEventCollection ***fEventColl; //!
  AliAnalysishCascadeEvent *fEvt;                                //!
  Int_t     fMaxPMult;
  Int_t     fMaxXiMult;

  Float_t     fPDGXi;
  Float_t     fPDGOmega;
  Float_t     fPDGp;
  Float_t     fPDGL;
  Float_t     fPDGpi;
  Float_t     fPDGfirst;
  Float_t     fPDGsecond;

  int fzVertexBins;
  int fnCentBins;
  short fnEventsToMix;


  TH1F              *fHistCentrality;                       //! histo to count the number of events in centrality
  TH1F              *fHistVertexDistribution;               //! Vertex distribution
  TH1F              *fHistMultiplicityOfMixedEvent;         //!
  TH1F              *fHistnTPCCrossedR;                     //!
  TH1F              *fHistRationTPCCrossedRnFind;           //!
  TH1F              *fHistSharedFrTPCcl;                    //!
  TH2F              *fHistprimpTPCdEdx;                     //!

  TH3F              *fHistyptProtons;                       //!
  TH2F              *fHistphietaProtons;                    //!
  TH2F              *fHistIPtoPVxyzTPC;                     //!
  TH2F              *fHistIPtoPVxyzGlobal;                  //!

  TH2F              *fHistpTOFmisvspt;                      //!
  TH2F              *fHistpTOFmisvsp;                       //!
  TH2F              *fHistpTOFnsigmavspt;                   //!
  TH2F              *fHistpTOFnsigmavsp;                    //!
  TH2F              *fHistpTOFsignalvsp;                    //!
  TH2F              *fHistpTOFsignalvspt;                   //!

  TH2F              *fHistpTOFTPCsignalvspt;                //!
  TH2F              *fHistProtonMultvsCent;                 //!

  TH2F              *fHistpTPCdEdx;                         //!
  TH2F              *fHistnTPCdEdx;                         //!
  TH2F              *fHistbTPCdEdx;                         //!

  TH1F              *fHistPosV0TPCClusters;                 //!
  TH1F              *fHistNegV0TPCClusters;                 //!
  TH1F              *fHistBachTPCClusters;                  //!

  TH1F              *fHistSharedFrTPCclPos;                 //!
  TH1F              *fHistSharedFrTPCclNeg;                 //!
  TH1F              *fHistSharedFrTPCclBach;                //!

  TH2F              *fHistyptXi;                            //!
  TH2F              *fHistphietaXi;                         //!
  TH1F              *fHistptL;                              //!
  TH1F              *fHistptpL;                             //!
  TH1F              *fHistptnL;                             //!
  TH1F              *fHistptbac;                            //!

  TH2F              *fHistInvMassXiMinus;                   //!
  TH2F              *fHistInvMassL;                         //!
  TH2F              *fHistInvMassXiPlus;                    //!
  TH2F              *fHistInvMassAntiL;                     //!
  TH2F              *fHistInvMassXiInPairs;                 //!
  TH2F              *fHistXiMultvsCent;                     //!
  AliCFContainer    *fCFContCascadeCuts;                    //!

  TH2F              *fHistpXiSignalRealKstar;               //!
  TH2F              *fHistapXiSignalRealKstar;              //!
  TH2F              *fHistpaXiSignalRealKstar;              //!
  TH2F              *fHistapaXiSignalRealKstar;             //!

  TH2F              *fHistpXiSignalBkgKstar;                //!
  TH2F              *fHistapXiSignalBkgKstar;               //!
  TH2F              *fHistpaXiSignalBkgKstar;               //!
  TH2F              *fHistapaXiSignalBkgKstar;              //!
  TH2F              *fHistFractionOfXiWithSharedDaughters;  //!
  TH2F              *fHistFractionOfaXiWithSharedDaughters; //!

  TH2F              *fHistpXibacDetaSDphiS;      //!
  TH2F              *fHistpXiposDetaSDphiS;      //!
  TH2F              *fHistpXinegDetaSDphiS;      //!
  TH2F              *fHistpaXibacDetaSDphiS;     //!
  TH2F              *fHistpaXiposDetaSDphiS;     //!
  TH2F              *fHistpaXinegDetaSDphiS;     //!
  TH2F              *fHistapXibacDetaSDphiS;     //! 
  TH2F              *fHistapXiposDetaSDphiS;     //!
  TH2F              *fHistapXinegDetaSDphiS;     //!
  TH2F              *fHistapaXibacDetaSDphiS;    //!
  TH2F              *fHistapaXiposDetaSDphiS;    //!
  TH2F              *fHistapaXinegDetaSDphiS;    //!

  TH2F              *fHistpXibacDetaSDphiSBkg;   //!
  TH2F              *fHistpXiposDetaSDphiSBkg;   //!
  TH2F              *fHistpXinegDetaSDphiSBkg;   //!
  TH2F              *fHistpaXibacDetaSDphiSBkg;  //!
  TH2F              *fHistpaXiposDetaSDphiSBkg;  //!
  TH2F              *fHistpaXinegDetaSDphiSBkg;  //!
  TH2F              *fHistapXibacDetaSDphiSBkg;  //! 
  TH2F              *fHistapXiposDetaSDphiSBkg;  //!
  TH2F              *fHistapXinegDetaSDphiSBkg;  //!
  TH2F              *fHistapaXibacDetaSDphiSBkg; //!
  TH2F              *fHistapaXiposDetaSDphiSBkg; //!
  TH2F              *fHistapaXinegDetaSDphiSBkg; //!

  TH1F              *fHistTrackBufferOverflow;   //! 

  TList             *fOutputContainer;           //! output list for the histograms


  //
  AliAnalysisTaskhCascadeFemto(const AliAnalysisTaskhCascadeFemto&); // not implemented
  AliAnalysisTaskhCascadeFemto& operator=(const AliAnalysisTaskhCascadeFemto&); // not implemented
  //
  ClassDef(AliAnalysisTaskhCascadeFemto, 1);
};

#endif
                 
