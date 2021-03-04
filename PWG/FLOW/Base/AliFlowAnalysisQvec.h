/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

/*************************************************
 * Save Qvec                                     *
 *                                               *
 * author: Shi Qiu                               *
 *         (s.qiu@nikhef.nl)                     *
 *************************************************/

#ifndef AliFlowAnalysisQvec_H
#define AliFlowAnalysisQvec_H

#include "TMatrixD.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TRandom3.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowAnalysisQvecEvent.h"
#include "TNamed.h"
#include <complex>
#include <cmath>

class TObjArray;
class TList;
class TFile;
class TGraph;
class TH1;
class TH3;
class TProfile;
class TProfile2D;
class TProfile3D;
class TDirectoryFile;
class TRandom3;
class TNtuple;
class THnSparse;

class AliFlowEventSimple;
class AliFlowTrackSimple;
class AliFlowCommonConstants;
class AliFlowCommonHist;
class AliFlowCommonHistResults;
class AliFlowVector;

//==============================================================================================================

class AliFlowAnalysisQvec : public TNamed {
public:
  AliFlowAnalysisQvec(const char* name="AliFlowAnalysisQvec",
                       Int_t nCen=7,
                       Double_t CenWidth=10.);
  virtual ~AliFlowAnalysisQvec();

  enum DataSet { k2010,
    k2011,
    k2015,
    k2015v6,
    k2015pidfix,
    kAny
  };

  enum InteractionRate {
    kHigh, // >0.9
    kLow,  // <0.6
    kPos, // mag.field pol.+
    kNeg,  // mag.field pol.-
    kAll
  };

  enum SelectCharge {
    kPosCh,
    kNegCh,
    kAllCh
  };

  enum POIExtraWeights {
    kNone,
    kEtaPhi,
    kEtaPhiCh,
    kEtaPhiVtx,
    kEtaPhiChPt,
    kEtaPhiRbR,
    kEtaPhiChRbR,
    kEtaPhiVtxRbR,
  };

  // 0.) methods called in the constructor:
  virtual void InitializeArraysForVarious();
  virtual void InitializeArraysForCRC();
  virtual void InitializeArraysForCRCVZ();
  virtual void InitializeArraysForCRCZDC();
  virtual void InitializeArraysForParticleWeights();
  virtual void InitializeCostantsForCRC();
  virtual void InitializeArraysForCME();
  
  // 1.) method Init() and methods called within Init():
  virtual void Init();
  virtual void InitializeArraysForQVec();
  virtual void CrossCheckSettings();
  virtual void CommonConstants(TString method);
  virtual void BookEverythingForIntegratedFlow();
  virtual void BookEverythingForVarious();
  virtual void BookAndNestAllLists();
  virtual void BookEverythingForQVec();
  virtual void BookEverythingForCME();
  virtual void BookAndFillWeightsHistograms();
  virtual void SetRunList();
  virtual void SetCentralityWeights();
  
  // 2.) method Make() and methods called within Make():
  virtual void Make(AliFlowEventSimple *anEvent);
  // 2a.) Common:
  virtual void CheckPointersUsedInMake();
  virtual void ResetEventByEventQuantities();
  // 2b.) Reference flow:
  
  // 2c.) Cross-checking reference flow correlations with nested loops:
  
  // 2d.) Differential flow:
  
  // 2e.) 2D differential flow:
  
  // 2f.) Other differential correlators (i.e. Teaney-Yan correlator):
  
  // 2g.) Distributions of reference flow correlations:

  // 2h.) Cross-checking differential flow correlations with nested loops:

  // 2i.) Charge-Rapidity Correlations
  virtual void CalculateCRCQVec();
  virtual void CalculateCMESPPP(); //@shi add CalculateCMESPPP() for spectator plane participant plane method
  virtual void RecenterCRCQVecZDC();
  virtual Bool_t PassCutZDCQVecDis(Double_t ZCRe, Double_t ZCIm, Double_t ZARe, Double_t ZAIm);
  virtual void RecenterCRCQVecZDC2(); //@Shi load my recentering file
  virtual void RecenterCRCQVecVZERO();
  virtual void PassQAZDCCuts();
  virtual Bool_t MultCut2015o();
  
  // 2h.) Various
  
  // 3.) method Finish() and methods called within Finish():
  virtual void Finish();
  virtual void CheckPointersUsedInFinish();
  // 3a.) integrated flow:
  
  //  nua:
  
  // 3b.) differential flow:
  
  // 3c.) 2D:

  // 3d.) Other differential correlators:

  // 3e.) Mixed harmonics:

  // 3f.) Bootstrap:

  // 3g.) CRC:
  
  // 3h.) Various:
  virtual Bool_t EvaulateIfSplitMergedTracks(AliFlowEventSimple* anEvent, AliFlowTrackSimple* aftsTrack, Int_t it1);
  virtual Double_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign);
  
  // 4.)  method GetOutputHistograms() and methods called within GetOutputHistograms():
  virtual void GetOutputHistograms(TList *outputListHistos);
  virtual void GetPointersForCommonHistograms();
  virtual void GetPointersForParticleWeightsHistograms();
  virtual void GetPointersForQVec();
  virtual void GetPointersForVarious();
  // 5.) other methods:
  virtual void WriteHistograms(TString outputFileName);
  virtual void WriteHistograms(TDirectoryFile *outputFileName);
  virtual Int_t GetCRCRunBin(Int_t RunNum);
  virtual Int_t GetCRCCenBin(Double_t Centrality);
  
  // **** SETTERS and GETTERS ****

  // 0.) base:
  void SetHistList(TList* const hlist) {this->fHistList = hlist;}
  TList* GetHistList() const {return this->fHistList;}

  // 1.) common:

  // 2a.) particle weights:
  void SetWeightsList(TList* const wlist) {this->fWeightsList = wlist;}
  TList* GetWeightsList() const {return this->fWeightsList;}
  void SetPhiWeights(TH1F* const histPhiWeights) {this->fPhiWeightsRPs = histPhiWeights;}
  TH1F* GetPhiWeights() const {return this->fPhiWeightsRPs;}
  void SetUsePtWeights(Bool_t const uPtW) {this->fUsePtWeights = uPtW;};
  Bool_t GetUsePtWeights() const {return this->fUsePtWeights;};
  void SetUsePhiEtaCuts(Bool_t const uPhiEtaW) {this->fUsePhiEtaCuts = uPhiEtaW;};
  Bool_t GetUsePhiEtaCuts() const {return this->fUsePhiEtaCuts;};
  void SetUseZDCESEMulWeights(Bool_t const uPhiEtaW) {this->fUseZDCESEMulWeights = uPhiEtaW;};
  Bool_t GetUseZDCESEMulWeights() const {return this->fUseZDCESEMulWeights;};
  void SetUseZDCESESpecWeights(Bool_t const uPhiEtaW) {this->fUseZDCESESpecWeights = uPhiEtaW;};
  Bool_t GetUseZDCESESpecWeights() const {return this->fUseZDCESESpecWeights;};
  void SetCutMultiplicityOutliers(Bool_t const uPhiEtaW) {this->fCutMultiplicityOutliers = uPhiEtaW;};
  Bool_t GetCutMultiplicityOutliers() const {return this->fCutMultiplicityOutliers;};
  
  // 2b.) event weights:
  void SetMultiplicityWeight(const char *multiplicityWeight) {*this->fMultiplicityWeight = multiplicityWeight;};
  
  // 3.) Reference flow:
  // Flags:
  void SetExactNoRPs(Int_t const enr) {this->fExactNoRPs = enr;};
  Int_t GetExactNoRPs() const {return this->fExactNoRPs;};
  // Reference flow profiles:
  
  // integrated flow histograms holding all results:
  
  // 4.) Differential flow:
  //  Flags:
  
  //  Profiles:
  //   1D:

  //   2D:

  //   Other differential correlators:

  // histograms:
  
  //  2D:
  
  // 5.) distributions of correlations:
  // profile:

  // flags:

  // # of bins for correlation axis in fDistributions[4], fCorrelation2468VsMult[4] and fCorrelationProduct2468VsMult[1]:

  // histograms:

  // min and max values of correlations (ci is correlations index [0=<2>,1=<4>,2=<6>,3=<8>]):

  // min and max values of correlation products:

  // min and max values of QvectorTerms:

  // x.) debugging and cross-checking:

  // 9.) Mixed harmonics:

  // 10.) Control histograms:
  
  // 11.) Bootstrap:

  // 12.) CRC
  void SetUseVZERO(Bool_t const cCRC) {this->fUseVZERO = cCRC;};
  Bool_t GetUseVZERO() const {return this->fUseVZERO;};
  void SetCRCVZEROCalibList(TList* const wlist) {this->fCRCVZEROCalibList = wlist;}
  TList* GetCRCVZEROCalibList() const {return this->fCRCVZEROCalibList;}
  void SetCalculateCRC(Bool_t const cCRC) {this->fCalculateCRC = cCRC;};
  Bool_t GetCalculateCRC() const {return this->fCalculateCRC;};
  void SetUseZDC(Bool_t const cCRC) {this->fUseZDC = cCRC;};
  Bool_t GetUseZDC() const {return this->fUseZDC;};
  void SetCRCQVecList(TList* const CRCL) {this->fCRCQVecList = CRCL;};
  void SetInvertZDC(Bool_t const cCRC) {this->fInvertZDC = cCRC;};
  Bool_t GetInvertZDC() const {return this->fInvertZDC;};
  void SetQAZDCCuts(Bool_t const cCRC) {this->fQAZDCCuts = cCRC;};
  Bool_t GetQAZDCCuts() const {return this->fQAZDCCuts;};
  void SetRecenterZDC(Bool_t const cCRC) {this->fRecenterZDC = cCRC;};
  Bool_t GetRecenterZDC() const {return this->fRecenterZDC;};
  void SetUseTracklets(Bool_t const cCRC) {this->fUseTracklets = cCRC;};
  void SetRemoveSplitMergedTracks(Bool_t const uPhiEtaW) {this->fRemoveSplitMergedTracks = uPhiEtaW;};
  Bool_t GetRemoveSplitMergedTracks() const {return this->fRemoveSplitMergedTracks;};
  void SetTestSin(Bool_t const cCRC) {this->fCRCTestSin = cCRC;};
  Bool_t GetTestSin() const {return this->fCRCTestSin;};
  void SetCalculateCME(Bool_t const cCRC) {this->fCalculateCME = cCRC;};
  Bool_t GetCalculateCME() const {return this->fCalculateCME;};
  void SetCRCQVecListRun(TList* const CRCL, Int_t r) {this->fCRCQVecListRun[r] = CRCL;};
  void SetDivSigma(Bool_t const cCRC) {this->fDivSigma = cCRC;};
  Bool_t GetDivSigma() const {return this->fDivSigma;};
  void SetStoreZDCQVecVtxPos(Bool_t const cCRC) {this->fStoreZDCQVecVtxPos = cCRC;};
  Bool_t GetStoreZDCQVecVtxPos() const {return this->fStoreZDCQVecVtxPos;};
  void SetCRCZDCCalibList(TList* const wlist) {this->fCRCZDCCalibList = wlist;}
  TList* GetCRCZDCCalibList() const {return this->fCRCZDCCalibList;}
  void SetCRCZDCResList(TList* const wlist) {this->fCRCZDCResList = wlist;}
  TList* GetCRCZDCResList() const {return this->fCRCZDCResList;}
  void SetCRCEtaRange(Double_t const etamin, Double_t const etamax) {this->fCRCEtaMin = etamin; this->fCRCEtaMax = etamax;};
  
  //@Shi my ZDC recenter calib hist
  void SetZDCCalibListFinalCommonPart(TList* const kList) {this->fZDCCalibListFinalCommonPart = (TList*)kList->Clone();};
  TList* GetZDCCalibListFinalCommonPart() const {return this->fZDCCalibListFinalCommonPart;};
  void SetZDCCalibListFinalRunByRun(TList* const kList) {this->fZDCCalibListFinalRunByRun = (TList*)kList->Clone();};
  TList* GetZDCCalibListFinalRunByRun() const {return this->fZDCCalibListFinalRunByRun;};
  
  void SetCRCZDC2DCutList(TList* const wlist) {this->fCRCZDC2DCutList = wlist;}
  void SetZDCESEList(TList* const kList) {this->fZDCESEList = kList;};
  TList* GetZDCESEList() const {return this->fZDCESEList;};
  // 12.a) EbE Corr:

  // 12.b) Final histo:

  // 12.c) Covariances:

  // 12.d) NUA corrections:

  // 12.e) Q Vectors:
  void SetCRCZDCQVecRes(TProfile* const TH, Int_t const r, Int_t const c) {this->fCRCZDCQVecRes[r][c] = TH;};
  TProfile* GetCRCZDCQVecRes(Int_t const r, Int_t const c) const {return this->fCRCZDCQVecRes[r][c];};
  void SetCRCZDCQVecAHist(TProfile* const TH, Int_t const r, Int_t const c) {this->fCRCZDCQVecA[r][c] = TH;};
  TProfile* GetCRCZDCQVecAHist(Int_t const r, Int_t const c) const {return this->fCRCZDCQVecA[r][c];};
  void SetCRCZDCQVecCHist(TProfile* const TH, Int_t const r, Int_t const c) {this->fCRCZDCQVecC[r][c] = TH;};
  TProfile* GetCRCZDCQVecCHist(Int_t const r, Int_t const c) const {return this->fCRCZDCQVecC[r][c];};
  void SetCRCZDCQVecACorrHist(TProfile* const TH, Int_t const r, Int_t const c) {this->fCRCZDCQVecACorr[r][c] = TH;};
  TProfile* GetCRCZDCQVecACorrHist(Int_t const r, Int_t const c) const {return this->fCRCZDCQVecACorr[r][c];};
  void SetCRCZDCQVecCCorrHist(TProfile* const TH, Int_t const r, Int_t const c) {this->fCRCZDCQVecCCorr[r][c] = TH;};
  TProfile* GetCRCZDCQVecCCorrHist(Int_t const r, Int_t const c) const {return this->fCRCZDCQVecCCorr[r][c];};
  
  // CRC VZERO:

  // CRC ZDC:
  
  // CRC2:
  void SetCenWeightsHist(TH1D* const n) {this->fCenWeightsHist = n;};
  TH1D* GetCenWeightsHist() const {return this->fCenWeightsHist;};
  void SetRefMultRbRPro(TProfile2D* const n) {this->fRefMultRbRPro = n;};
  void SetAvEZDCRbRPro(TProfile2D* const A, TProfile2D* const B) {this->fAvEZDCCRbRPro = A; this->fAvEZDCARbRPro = B;};
  void SetPtWeightsHist(TH1D* const n, Int_t c) {this->fPtWeightsHist[c] = n;};
  TH1D* GetPtWeightsHist(Int_t c) const {return this->fPtWeightsHist[c];};
  void SetCenvsMul(TH2F* const n, Int_t const h) {this->fhCenvsMul[h] = n;};
  TH2F* GetCenvsMul(Int_t const h) const {return this->fhCenvsMul[h];};
  void SetZDCESEMultWeightsHist(TH2F* const n, Int_t h) {this->fZDCESEMultWeightsHist[h] = n;};
  TH2F* GetZDCESEMultWeightsHist(Int_t h) const {return this->fZDCESEMultWeightsHist[h];};
  void SetZDCESESpecWeightsHist(TH2F* const n, Int_t h) {this->fZDCESESpecWeightsHist[h] = n;};
  TH2F* GetZDCESESpecWeightsHist(Int_t h) const {return this->fZDCESESpecWeightsHist[h];};
  void SetZNvsCen(TH2F* const n, Int_t const h) {this->fhZNvsCen[h] = n;};
  TH2F* GetZNvsCen(Int_t const h) const {return this->fhZNvsCen[h];};
  void SetZNCvsZNA(TH2F* const n, Int_t const h) {this->fhZNCvsZNA[h] = n;};
  TH2F* GetZNCvsZNA(Int_t const h) const {return this->fhZNCvsZNA[h];};
  void SetZNvsMul(TH2F* const n) {this->fhZNvsMul = n;};
  TH2F* GetZNvsMul() const {return this->fhZNvsMul;};
  // Flow QC

  // Flow Generic Framework

  // sub-sampling

  // in wide pt bins

  // Flow SP ZDC

  // v1

  // Flow SP VZ
  Int_t GetnRun() const {return this->fCRCnRun;};
  // CME:

  // CME TPC only:

  // CME TPC-ZDCs:
  
  // EbE Flow
  
  // 15.) Various
  void SetVariousList(TList* const Various) {this->fVariousList = Various;};
  void SetDataSet(DataSet set) {this->fDataSet = set;};
  DataSet GetDataSet() const {return this->fDataSet;}
  void SetInteractionRate(InteractionRate set) {this->fInteractionRate = set;};
  InteractionRate GetInteractionRate() const {return this->fInteractionRate;}
  void SetSelectCharge(SelectCharge set) {this->fSelectCharge = set;};
  SelectCharge GetSelectCharge() const {return this->fSelectCharge;}
  void SetPOIExtraWeights(POIExtraWeights set) {this->fPOIExtraWeights = set;};
  POIExtraWeights GetPOIExtraWeights() const {return this->fPOIExtraWeights;}
  void SetRunNumber(Int_t const n) {this->fRunNum = n;};
  Int_t GetRunNumber() const {return this->fRunNum;}
  void SetMinMulZN(Int_t weights) {this->fMinMulZN = weights;};
  Int_t GetMinMulZN() const {return this->fMinMulZN;};
  
private:

  AliFlowAnalysisQvec(const AliFlowAnalysisQvec& afawQc);
  AliFlowAnalysisQvec& operator=(const AliFlowAnalysisQvec& afawQc);

  // 0.) base:
  TList* fHistList; //! base list to hold all output object
  
  // 0.1) event list
  TList *fEventList;
  
  // 1.) common:
  
  // 2a.) particle weights:
  TList *fWeightsList; // list to hold all histograms with particle weights: fUseParticleWeights, fPhiWeights, fPtWeights and fEtaWeights
  Bool_t fUsePtWeights; // use pt weights
  Bool_t fUsePhiEtaCuts; // use phi,eta cuts (for NUA)
  Bool_t fUseZDCESEMulWeights;       // use ZDC-ESE mult. weights
  Bool_t fUseZDCESESpecWeights;       // use ZDC-ESE spec. weights
  Bool_t fCutMultiplicityOutliers;  // cut on reference multiplicity
  // TProfile *fUseParticleWeights; //! profile with three bins to hold values of fUsePhiWeights, fUsePtWeights and fUseEtaWeights
  // TH1F *fPhiWeightsPOIs[2]; //! histogram holding phi weights
  // TH1D *fPtWeightsPOIs[2]; //! histogram holding pt weights
  // TH1D *fEtaWeightsPOIs[2]; //! histogram holding eta weights
  // TH2D *fPhiEtaWeightsPOIs[2]; //! histogram holding phi,eta weights
  TH1F *fPhiWeightsRPs; //!
  // TH1D *fPtWeightsRPs; //!
  // TH1D *fEtaWeightsRPs; //!
  // TH1F *fPhiDistrRefPOIs[2]; //! histogram holding phi weights
  // TH1D *fPtDistrRefPOIs[2]; //! histogram holding pt weights
  // TH1D *fEtaDistrRefPOIs[2]; //! histogram holding eta weights
  // TH2D *fPhiEtaDistrRefPOIs[2]; //! histogram holding phi,eta weights
  // TH1F *fPhiDistrRefRPs; //!
  // TH1D *fPtDistrRefRPs; //!
  // TH1D *fEtaDistrRefRPs; //!
  // TH2D *fPhiEtaDistrRefRPs; //!
  // TH1D *fPtWeights[2]; //! histogram holding pt weights

  // 2b.) event weights:
  TString *fMultiplicityWeight; //! event-by-event weights for multiparticle correlations

  // 3.) integrated flow
  //  3a.) lists:

  //  3b.) flags:
  Int_t fExactNoRPs; // when shuffled, select only this number of RPs for the analysis

  //  3c.) event-by-event quantities:
  TMatrixD *fReQGF; //! fReQ[m][k] = sum_{i=1}^{M} w_{i}^{k} cos(m*phi_{i})
  TMatrixD *fImQGF; //! fImQ[m][k] = sum_{i=1}^{M} w_{i}^{k} sin(m*phi_{i})
  const static Int_t fkGFPtB = 8;
  TMatrixD *fReQGFPt[fkGFPtB]; //! fReQ[m][k] = sum_{i=1}^{M} w_{i}^{k} cos(m*phi_{i})
  TMatrixD *fImQGFPt[fkGFPtB]; //! fImQ[m][k] = sum_{i=1}^{M} w_{i}^{k} sin(m*phi_{i})
  
  Double_t fNumberOfRPsEBE; // # of Reference Particles
  Double_t fNumberOfPOIsEBE; // # of Particles of Interest
  Double_t fReferenceMultiplicityEBE; // reference multiplicity
  Double_t fReferenceMultiplicityRecEBE; // reference multiplicity - <reference multiplicity>
  Double_t fCentralityEBE; // centrality percentile
  Double_t fNewCentralityEBE; // new centrality percentile
  Double_t fNewMetricLEBE; // new metric L
  Double_t fNewMetricDEBE; // new metric D
  Double_t fNewMetricL2EBE; // new metric L
  Double_t fNewMetricD2EBE; // new metric D
  Double_t fCentralityCL1EBE; // centrality (CL1) percentile
  Double_t fNITSCL1EBE; // centrality (TRK) percentile
  Double_t fCentralityTRKEBE; // centrality (TRK) percentile
  
  //  3d.) profiles:

  //  3e.) histograms with final results:

  // 4.) differential flow
  //  4a.) lists:

  //    4aa.) nested list in list fDiffFlowProfiles:

  //    4ab.) nested list in list fDiffFlowResults:

  //  4b.) flags:

  //  4c.) event-by-event quantities:
  //   1D:

  //   2D:

  //  4d.) profiles:
  //   1D:

  //   2D:

  //  4e.) histograms holding final results:
  //   1D:

  //   2D:

  // 6.) distributions:

  // 8.) debugging and cross-checking:
  
    // 9.) mixed harmonics:
  //  9a.) lists:

  //  9b.) flags:

  //  9c.) profiles:

  //  9d.) results:

  //  9e.) statistical error propagation:

  // 10.) Control histograms:
  //  10a.) list:

  //  10b.) flags:

  //  10c.) histograms:

  // 11.) Bootstrap:
  //  11a) lists:

  //  11b) flags:
  TRandom3 *fRandom; //! local random generator
  //  11c) profiles:

  //  11d) histograms:

  // 12.) CRC
  TList *fTempList; //! list to hold temp histograms
  Bool_t fCalculateCRC; // calculate CRC
  Int_t fRunNum;
  Int_t fCachedRunNum;
  Int_t fRunBin;
  Int_t fCenBin;
  Bool_t fUseVZERO;
  Bool_t fUseZDC;
  Bool_t fInvertZDC;
  Bool_t fRecenterZDC;
  Int_t fCRCnCen;
  Bool_t fRemoveSplitMergedTracks;
  Bool_t fCRCTestSin;
  Bool_t fCalculateCME;
  Bool_t fDivSigma;
  Double_t fCRCEtaMin;
  Double_t fCRCEtaMax;
  
  const static Int_t fCRCMaxnCen = 10;
  const static Int_t fCRCMaxnRun = 211;
  const static Int_t fCRCnHar = 3;
  
  TH2F *fPtWeightsCent; //!
  TH3D *fPhiEtaWeightsVtx[fCRCMaxnCen]; //!
  TH3D *fPhiEtaWeights; //!
  TH3D *fPhiEtaWeightsCh[2]; //!
  TH3D *fPhiEtaWeightsChPt[2][3]; //!
  TH3D *fPhiEtaRbRWeights; //!
  TH3D *fPhiEtaRbRWeightsCh[2]; //!
  
  // Q vectors
  const static Int_t fCRCQVecnCR = 64;
  Int_t fCRCnRun;
  DataSet fDataSet;
  InteractionRate fInteractionRate;
  SelectCharge fSelectCharge;
  POIExtraWeights fPOIExtraWeights;
  TArrayI fRunList;    // Run list
  TArrayD fAvVtxPosX;    // Run list
  TArrayD fAvVtxPosY;    // Run list
  TArrayD fAvVtxPosZ;    // Run list
  TList *fCRCZDCResList; //! ZDC rescaling list
  TList *fCRCQVecListTPC; //! Q Vectors list TPC
  
  //@Shi add ave vtx IR split
  TArrayD fAvVtxPosX15oIRSplit;    // Run list
  TArrayD fAvVtxPosY15oIRSplit;    // Run list
  TArrayD fAvVtxPosZ15oIRSplit;    // Run list
  
  TList *fCRCZDCCalibList; //! ZDC calibration
  TList *fCRCVZEROCalibList; //! ZDC calibration
  TList *fCRCQVecListRun[fCRCMaxnRun]; //! Q Vectors list per run
  TList *fCRCQVecList; //! Q Vectors list
  TList *fZDCCalibListFinalCommonPart; //! Shi my ZDC calbration split to run independent part and run-by-run calib due to limited size of content that can be held by TList
  TList *fZDCCalibListFinalRunByRun; // Shi run dependent calib part, run-by-run was set in AliAnalysisTaskCRC
  TList *fZDCESEList; //! ZDC ESE
  TList *fCRCZDC2DCutList; //! ZDC 2D cut
  // temp
  TProfile *fZDCQHist[12];  //! Run-by-run ZDCQvecHist
  const static Int_t fkNZDCResHist = 4;
  TH1D *fZDCResHist[fkNZDCResHist]; //!
  TProfile2D *fhAvRefMulRbR; //! Average reference multiplicity vs run vs centrality
  TProfile2D *fhAvQMCRbR; //!
  TProfile2D *fhAvQMARbR; //!
  TProfile *fhAvAbsOrbit; //!
  TProfile3D *fZDCVtxHist[4]; //! Run-by-run vtxZDCQvec
  TProfile2D *fZDCEcomHist[4];//! Run-by-run vtxZDCQvec
  TProfile2D *fZDCEcomTotHist[4];//! Run-by-run vtxZDCQvec
  TProfile3D *fZDCEcomTotvsVtxHist[12];//! Run-by-run vtxZDCQvec
  TProfile3D *fZDCVtxCenHist[10][4]; //! Run-by-run vtxZDCQvec
  TProfile3D *fZDCVtxCenHistMagPol[10][8]; //! Run-by-run vtxZDCQvec

  TProfile2D *fVZEROCenHist[3];//! Run-by-run VZERO Q-vector (harmonics 1-3)
  TH1D *fZDCESEAvHist[2]; //!
  const static Int_t fkNsteps = 14;
  TProfile *fCRCZDCQVecRes[fCRCMaxnRun][4]; //! Q Vectors Resolution Terms
  Bool_t fStoreZDCQVecVtxPos; //
  TH2D *fZDCEPHist[2];  //! Run-by-run ZDCQvecHist
  TH3D *fZDCVtxFitHist[4]; //!
  TH1D *fZDCVtxFitCenProjHist[4][3]; //!
  TH3D *fZDCVtxFitHist2[4]; //!
  TH1D *fZDCVtxFitCenProjHist2[4][3]; //!
  TH3D *fZDCBinsCenRefMult[10]; //!
  
  TProfile2D *fZDCBinsCenRefMultRbR[4]; //!
  TProfile2D *fZDCBinsCenRefMultTot[4]; //!
  TProfile *fZDCBinsCenRefMultRbRProf[10][4]; //!
  TProfile *fZDCBinsCenRefMultTotProf[10][4]; //!
  TH1D *fZDCBinsCenRefMultRbRProj[10][4]; //!
  TH1D *fZDCBinsCenRefMultTotProj[10][4]; //!
  TH3D *fZDCBinsVtxCenEZDC[3][4]; //!
  TH3D *fZDCQVecVtxCenEZDC3D[10][10][4]; //!
  TH3D *fZDCQVecVtxCenEZDCFit0; //!
  TH3D *fZDCQVecVtxCenEZDCFit1; //!
  TH2D *fCRCZDC2DCutZDCC[2][10]; //!
  TH2D *fCRCZDC2DCutZDCA[2][10]; //!
  
  //@Shi temp ZDC calib histograms
  TProfile *fAvr_Run_CentQ[4]; //!
  TProfile3D *fAvr_Cent_VtxXYZQ[20][4]; //!
  TProfile3D *fAvr_Run_VtxXYZQ[4]; //!
  
  TProfile2D *fCRCZDCQVecCorSteps[4]; //!
  

  
  // test
  TH2D* fCRCZDCQVecDis[2][fCRCMaxnCen][2]; //!
  
  // CRCVZERO
  AliFlowVector fVZFlowVect[2][fCRCnHar];
  
  // CRCZDC
  AliFlowVector fZDCFlowVect[2];
  
  // CRC Pt differential

  // CME
  const static Int_t fCMEnEtaBin = 2;
  TH1D *fCMEQRe[4][fCRCnHar]; //! real part [0=pos,1=neg][0=back,1=forw][m]
  TH1D *fCMEQIm[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][m]
  TH1D *fCMEMult[4][fCRCnHar]; //! imaginary part [0=pos,1=neg][0=back,1=forw][p][k]
  Int_t fZDCESEclEbE;
  
  //@shi add Qvector for both charge (begin)
  TH1D *fCMEQReBothCharge[2][fCRCnHar]; //! real part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEQImBothCharge[2][fCRCnHar]; //! imaginary part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  TH1D *fCMEMultBothCharge[2][fCRCnHar]; //! imaginary part [2]: power of weight, [fCRCnHar]: cos((h+1)*phi)
  //@shi add Qvector for both charge (end)
  
  // CRC2
  
  // Flow all
  
  // Flow SP ZDC

  // Flow QC
  
  // flow Generic Framework

  // sub-sampling

  // in wide pt bins

  // SC w ZDC

  // Flow SP VZ
  
  // EbE Flow

  // Various:
  TList *fVariousList;
  const static Int_t fZDCESEnCl = 5;
  TH2F* fhCenvsMul[fZDCESEnCl+1]; //! cen vs mul
  TH1D* fCenWeigCalHist; //! Centrality weights
  TProfile *fCRCZDCQVecA[fCRCMaxnRun][2]; //! Q Vectors ZDCN-A
  TProfile *fCRCZDCQVecC[fCRCMaxnRun][2]; //! Q Vectors ZDCN-C
  TProfile *fCRCZDCQVecACorr[fCRCMaxnRun][2]; //! Q Vectors ZDCN-A
  TProfile *fCRCZDCQVecCCorr[fCRCMaxnRun][2]; //! Q Vectors ZDCN-C

  Double_t fCenWeightEbE;
  TProfile2D* fRefMultRbRPro; //! run-by-run average reference multiplicity
  TProfile2D* fAvEZDCCRbRPro; //! run-by-run average EZDC-C
  TProfile2D* fAvEZDCARbRPro; //! run-by-run average EZDC-A
  TH1D *fEventCounter; //! Event counter for different methods
  TH1D* fPtWeightsHist[10]; //! Pt weights
  TH1F* fMultCutMin; //!
  TH1F* fMultCutMax; //!
  TH1F* fMultCutAv; //!
  TH2F* fZDCESEMultWeightsHist[5]; //! ZDC-ESE mult weights
  TH2F* fZDCESESpecWeightsHist[5]; //! ZDC-ESE mult weights
  TH1D* fCenWeightsHist; //! Centrality weights
  TH3F* fTwoTrackDistanceLS[2]; //!
  TH3F* fTwoTrackDistanceUS[2]; //!
  TH2F* fhZNvsCen[2]; //! cen vs mul
  Double_t *fCorrMap; //!
  TH1F* fEZNCutMin; //!
  TH1F* fEZNCutMax; //!
  TH2F* fhZNCvsZNA[fCRCMaxnCen]; //! ZNA-ZNC correlation
  TH2F* fhZNvsMul; //! cen vs mul
  
  Bool_t fbFlagIsPosMagField;
  Bool_t fbFlagIsBadRunForC34;
  
  Double_t fVtxPos[3]; // primary vertex position (x,y,z)
  Double_t fVtxPosCor[3]; // primary vertex position (x,y,z), re-centered at 0
  Double_t fVtxPosCor15oIRSplit[3];  //@Shi primary vertex position (x,y,z), re-centered at 0

  TF1 *fPolMin[2]; //!
  TF1 *fPolMax[2]; //!
  TF1 *fPolAv[2]; //!
  TF1 *fPolDist[2]; //!
  TF1 *fPolSlope[2]; //!
  
  TGraph *fCenMetric; //!
  
  Double_t fZNCQ0; // common tower energy from ZNC-C
  Double_t fZNAQ0; // common tower energy from ZNC-A
  Double_t fZNCen; // total energy from ZNC-C
  Double_t fZNAen; // total energy from ZNC-A
  Double_t fZPCen; // total energy from ZPC-C
  Double_t fZPAen; // total energy from ZPC-A
  Double_t fEnNucl; // energy per nucleon (GeV)
  Bool_t fQAZDCCuts;
  Bool_t fQAZDCCutsFlag;
  Int_t fMinMulZN;
  Bool_t fUseTracklets;
  
  const static Int_t fZDCESEnPol=4;
  TF1 *fPolCuts[fZDCESEnPol]; //!
  TH1D *fZDCESECutsHist[fZDCESEnPol]; //!
  
  TTree* treeEvent; //!
  AliFlowAnalysisQvecEvent* fpQvecEvent; //!
  
  ClassDef(AliFlowAnalysisQvec, 1);

};

//================================================================================================================

#endif
