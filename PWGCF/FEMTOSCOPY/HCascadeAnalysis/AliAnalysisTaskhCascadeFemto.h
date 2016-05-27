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
  enum xicuts_t {kloose, kdefault, ktight};

  void SetAnalysisType (const char* analysisType = "ESD")   { fAnalysisType = analysisType; }
  void SetMassHist     (Float_t nbins = 100, Float_t lowlim = 2.028, Float_t uplim = 2.228) {fNbinsMass = nbins; fLowMassLim = lowlim; fUpMassLim = uplim; }
  void SetMassHistGrandmother (Float_t nbins = 100, Float_t lowlim = 2.028, Float_t uplim = 2.228) {fNbinsMassGm = nbins; fLowMassLimGm = lowlim; fUpMassLimGm = uplim;}
  void SetUseStandardCuts (Bool_t usestandardesdtrackcuts) {fUseStandardCuts = usestandardesdtrackcuts;}
  void SetCentrality   (Float_t lowlimcent = 0., Float_t uplimcent = 90.) { fCentrLowLim = lowlimcent;  fCentrUpLim = uplimcent; }
  void SetReadMCTruth (Bool_t readmctruth) {fReadMCTruth = readmctruth;}
  void SetUseContainer (Bool_t kusecontainer) { fUseContainer = kusecontainer;}

  void DoPairshCasc (const Float_t centralityBin); 
  void DoPairshh    (const Float_t centralityBin, int fieldsign);

  Int_t CheckDaughterTrack (Int_t xiIndex, Int_t daughterIndex, Int_t*  checkeddaughters);
  void SelectBestCandidate ();
  double CalculateKstar(double momentum1[3], double momentum2[3], double mass1, double mass2); 
  AliReconstructedProton::MCProtonOrigin_t ProtonOrigin(int trackLabel, TClonesArray *arrayMC, double* pMomentumTruth);
  AliReconstructedXi::MCXiOrigin_t XiOrigin(int labelB, int labelP, int labelN, TClonesArray *arrayMC, double* xiMomentumTruth);
  void SetSftPosR125(AliVTrack *track, const Float_t bfield, Double_t priVtx[3], Double_t fXSftR125[3] ); 
  double CalculateDphiSatR12m(Double_t pos1SftR125[3], Double_t pos2SftR125[3]);
  double CalculateDphiSatR12mAnal(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2, Double_t* dps2);
  void SetFirstParticle(firstpart_t firstpart) {fFirstpart = firstpart;}
  void SetSecondParticle(firstpart_t secondpart) {fSecondpart = secondpart;}
  void SetXiCuts(xicuts_t xicutset) {fXicutSet = xicutset;}

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
  void SetNEventsToMix(short nevmixing) { fnEventsToMix = nevmixing;}
  void SetMomentumLimitForTOFPID(Float_t momemtumlimitforTOFPID) { fMomemtumLimitForTOFPID = momemtumlimitforTOFPID;} 
  void SetApplyRatioCrRnFindCut(Bool_t kapplycrrowsnfindcut) { fkApplyRatioCrRnFindCut = kapplycrrowsnfindcut;}
  void SetApplyCrossedRowCut(Bool_t kapplycrrowscut) { fkCutOnCrossedRows = kapplycrrowscut;}
  void SetCutOnTPCIP(Bool_t kcutontpcip) { fkCutOnTPCIP = kcutontpcip;}
  void SetIPCutxy(Float_t ipcutxy) { fIPCutxy = ipcutxy;}
  void SetIPCutz (Float_t ipcutz)  { fIPCutz = ipcutz;}
  void SetMinPtPrim(Float_t minptforprim) { fMinPtForPrim  = minptforprim;}  
  void SetMaxPtPrim(Float_t maxptforprim) { fMaxPtForPrim  = maxptforprim;}
  void SetMinPtCasc(Float_t minptforcasc) { fMinPtForCasc  = minptforcasc;}
  void SetMaxPtCasc(Float_t maxptforcasc) { fMaxPtForCasc  = maxptforcasc;}
  void SetIPCutBac(Float_t ipcutbac) { fIPBac = ipcutbac;}
  void SetApplyYcutCasc(Bool_t applyycutcasc) { fkApplyYcutCasc = applyycutcasc;}
  void SetPropagateGlobal(Bool_t propagateglobal) { fkPropagateGlobal = propagateglobal;}
  void SetPropagateAtFixedR(Bool_t propagatefixedr) { fkPropagateAtFixedR = propagatefixedr;}

  void SetCutOnttcProp(Bool_t kcutonttcprop) { fkCutOnTtcProp = kcutonttcprop;}
  void SetApplySharedDaughterCutXi(Bool_t kapplyshdaucutxi) { fkApplySharedDaughterCutXi = kapplyshdaucutxi;}
  void SetPosR125(AliVTrack *track, const Float_t bfield, Double_t posSftR125[3] );  
  Double_t EtaS( Double_t posSftR125[3] ) const; 
  Double_t ThetaS( Double_t posSftR125[3] ) const;

 private:
 
  AliESDEvent *fESDevent;                         //! 
  AliAODEvent *fAODevent;                         //! 


  TString         fAnalysisType;                  // "ESD" or "AOD" analysis type       
  firstpart_t fFirstpart;                         //
  firstpart_t fSecondpart;                        //
  xicuts_t fXicutSet;                             //
  Float_t fMassWindowCascades;                    //
  Float_t fnSigmaTPCPIDfirstParticle;             //     
  Float_t fnSigmaTPCTOFPIDfirstParticle;          //
  Float_t fnSigmaTPCPIDsecondParticleDau;         //
  Bool_t fkCascadeSideBands;                      //
  Bool_t fReadMCTruth;                            // if read MC truth
  Bool_t fUseContainer;                           //
  Bool_t fUseStandardCuts;                        // if to use standard ESD track cuts or user-defined ones 
  Bool_t fkApplyTtc;                              //
  Float_t fDphisMin;                              //
  Float_t fDetasMin;                              // 
                 
  Float_t fMomemtumLimitForTOFPID;                //
  Bool_t fkApplyRatioCrRnFindCut;                 //
  Bool_t fkCutOnCrossedRows;                      //
  Bool_t fkCutOnTPCIP;                            // 
  Float_t fIPCutxy;                               //
  Float_t fIPCutz;                                //
  Float_t fMinPtForPrim;                          //
  Float_t fMaxPtForPrim;                          //
  Float_t fMinPtForCasc;                          // 
  Float_t fMaxPtForCasc;                          //
  Float_t fIPBac;                                 //              
  Float_t fDcaXiDaughters;                        //  
  Float_t fXiCosineOfPointingAngle;               //
  Float_t fXiRadiusMin;                           //
  Float_t fXiRadiusMax;                           //
  Float_t fCtau;                                  //
  Float_t fIPV0ToPrimVertexXi;                    // 
  Float_t fMassWindowL;                           //
  Float_t fV0CosineOfPointingAngle;               //
  Float_t fV0RadiusXiMax;                         //
  Float_t fV0RadiusXiMin;                         //
  Float_t fDcaV0Daughters;                        //
  Float_t fIPPosToPrimVertexXi;                   //
  Float_t fIPNegToPrimVertexXi;                   //

  Bool_t fkApplyYcutCasc;                         // 
  Bool_t fkPropagateGlobal;                       //
  Bool_t fkPropagateAtFixedR;                     //
  Bool_t fkCutOnTtcProp;                          //  
  Bool_t fkApplySharedDaughterCutXi;              //

  AliESDtrackCuts    *fESDtrackCuts;              //! basic cut variables for tracks added ! not sure
  AliPIDResponse     *fPIDResponse;               //! PID response object

  Float_t fCentrLowLim;                           // centrality settings: lower limit
  Float_t fCentrUpLim;                            // upper limit

  Float_t fNbinsMass;                             // Inv mass histo settings for a two-body decay
  Float_t fLowMassLim;                            //
  Float_t fUpMassLim;                             //
                
  Float_t fNbinsMassGm;                           // Inv mass histo settings for a casade-like decay 
  Float_t fLowMassLimGm;                          //
  Float_t fUpMassLimGm;                           //

  // Store pointers to global tracks for pid and dca

  Int_t     *farrGT;                                // !
  UShort_t  fTrackBufferSize;                       // Size fo the above array, ~12000 for PbPb

  AliAnalysishCascadeEventCollection ***fEventColl; //!
  AliAnalysishCascadeEvent *fEvt;                   //!
  Int_t     fMaxPMult;                              //
  Int_t     fMaxXiMult;                             // 

  Double_t     fPDGXi;                              //
  Double_t     fPDGOmega;                           //
  Double_t     fPDGp;                               //
  Double_t     fPDGL;                               //
  Double_t     fPDGpi;                              //                     
  Double_t     fPDGfirst;                           //
  Double_t     fPDGsecond;                          //  

  int fzVertexBins;                                 //
  int fnCentBins;                                   //
  short fnEventsToMix;                              //
     

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

  TH2F              *fHistpTOFmisvsp;                       //!
  TH2F              *fHistpTOFnsigmavsp;                    //!
  TH2F              *fHistpTOFsignalvsp;                    //!
  TH3F              *fHistpnsigTOFnsigTPC;                  //!
  TH3F              *fHistpsignalTOFsignalTPC;              //!
  TH2F              *fHistProtonMultvsCent;                 //!
  TH2F              *fHistMCPrimProtons;                    //!
  TH2F              *fHistMCFromWdecayProtons;              //!
  TH2F              *fHistMCFromMaterialProtons;            //! 
  TH2F              *fHistMCOtherProtons;                   //!
  TH2F              *fHistMCPrimAProtons;                   //!
  TH2F              *fHistMCFromWdecayAProtons;             //!
  TH2F              *fHistMCFromMaterialAProtons;           //! 
  TH2F              *fHistMCOtherAProtons;                  //!
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
  TH2F              *fHistIPtoPVxyGlobalvspt;               //!
//  AliCFContainer    *fCFContCascadeCuts;                    //!

  TH2D              *fHistpXiSignalRealKstar;               //!
  TH2D              *fHistapXiSignalRealKstar;              //!
  TH2D              *fHistpaXiSignalRealKstar;              //!
  TH2D              *fHistapaXiSignalRealKstar;             //!
  TH2D              *fHistpXiSignalMixedKstargenvsrec;      //!
  TH2D              *fHistapXiSignalMixedKstargenvsrec;     //!
  TH2D              *fHistpaXiSignalMixedKstargenvsrec;     //!
  TH2D              *fHistapaXiSignalMixedKstargenvsrec;    //!

  TH2D              *fHistpXiSignalBkgKstar;                //!
  TH2D              *fHistapXiSignalBkgKstar;               //!
  TH2D              *fHistpaXiSignalBkgKstar;               //!
  TH2D              *fHistapaXiSignalBkgKstar;              //!
  TH2F              *fHistFractionOfXiWithSharedDaughters;  //!
  TH2F              *fHistTotMaxFractionOfXiWithSharedDaughters; //!
 
  TH2D              *fHistpXibacDetaSDphiS;      //!
  TH2D              *fHistpXiposDetaSDphiS;      //!
  TH2D              *fHistpXinegDetaSDphiS;      //!
  TH2D              *fHistpaXibacDetaSDphiS;     //!
  TH2D              *fHistpaXiposDetaSDphiS;     //!
  TH2D              *fHistpaXinegDetaSDphiS;     //!
  TH2D              *fHistapXibacDetaSDphiS;     //! 
  TH2D              *fHistapXiposDetaSDphiS;     //!
  TH2D              *fHistapXinegDetaSDphiS;     //!
  TH2D              *fHistapaXibacDetaSDphiS;    //!
  TH2D              *fHistapaXiposDetaSDphiS;    //!
  TH2D              *fHistapaXinegDetaSDphiS;    //!

  TH2D              *fHistpXibacDetaSDphiSBkg;   //!
  TH2D              *fHistpXiposDetaSDphiSBkg;   //!
  TH2D              *fHistpXinegDetaSDphiSBkg;   //!
  TH2D              *fHistpaXibacDetaSDphiSBkg;  //!
  TH2D              *fHistpaXiposDetaSDphiSBkg;  //!
  TH2D              *fHistpaXinegDetaSDphiSBkg;  //!
  TH2D              *fHistapXibacDetaSDphiSBkg;  //! 
  TH2D              *fHistapXiposDetaSDphiSBkg;  //!
  TH2D              *fHistapXinegDetaSDphiSBkg;  //!
  TH2D              *fHistapaXibacDetaSDphiSBkg; //!
  TH2D              *fHistapaXiposDetaSDphiSBkg; //!
  TH2D              *fHistapaXinegDetaSDphiSBkg; //!

  TH1F              *fHistTrackBufferOverflow;   //! 

  TList             *fOutputContainer;           //! output list for the histograms


  //
  AliAnalysisTaskhCascadeFemto(const AliAnalysisTaskhCascadeFemto&); // not implemented
  AliAnalysisTaskhCascadeFemto& operator=(const AliAnalysisTaskhCascadeFemto&); // not implemented
  //
  ClassDef(AliAnalysisTaskhCascadeFemto, 8);
};

#endif
                 
