#ifndef ALIANALYSISEVENT_H
#define ALIANALYSISEVENT_H

#include "TObject.h"
#include "AliTOFPIDResponse.h"
#include "TF1.h"

class AliAnalysisEvent :
public TObject
{

 public:

  AliAnalysisEvent(); // default constructor
  AliAnalysisEvent(const AliAnalysisEvent &source); // copy constructor
  AliAnalysisEvent &operator=(const AliAnalysisEvent &source); // operator=
  virtual ~AliAnalysisEvent(); // default destructor

  Bool_t IsCollisionCandidate() const {return fIsCollisionCandidate;}; // getter
  Bool_t HasVertex() const {return fHasVertex;}; // getter
  Float_t GetVertexZ() const {return fVertexZ;}; // getter
  UChar_t GetCentralityQuality() const {return fCentralityQuality;}; // getter
  Float_t GetCentralityPercentile(Int_t i) const {return i < kNCentralityEstimators ? fCentralityPercentile[i] : (Float_t)GetReferenceMultiplicity();}; // getter
  Int_t GetReferenceMultiplicity() const {return fReferenceMultiplicity;}; // getter
  Float_t *GetTimeZeroTOF() {return fTimeZeroTOF;}; // getter
  Float_t *GetTimeZeroTOFSigma() {return fTimeZeroTOFSigma;}; // getter
  Float_t *GetTimeZeroT0() {return fTimeZeroT0;}; // getter
  Float_t GetTimeZeroTOF(Int_t i) {return i < 10 ? fTimeZeroTOF[i] : 0.;}; // getter
  Float_t GetTimeZeroTOFSigma(Int_t i) {return i < 10 ? fTimeZeroTOFSigma[i] : 0.;}; // getter
  Float_t GetTimeZeroT0(Int_t i) {return i < 3 ? fTimeZeroT0[i] : 0.;}; // getter
  Float_t GetMCTimeZero() const {return fMCTimeZero;}; // getter

  void SetIsCollisionCandidate(Bool_t value) {fIsCollisionCandidate = value;}; // setter
  void SetIsEventSelected(UInt_t value) {fIsEventSelected = value;}; // setter
  void SetIsPileupFromSPD(Bool_t value) {fIsPileupFromSPD = value;}; // setter
  void SetHasVertex(Bool_t value) {fHasVertex = value;}; // setter
  void SetVertexZ(Float_t value) {fVertexZ = value;}; // setter
  void SetCentralityQuality(UChar_t value) {fCentralityQuality = value;}; // setter
  void SetCentralityPercentile(Int_t icent, Float_t value) {if (icent < kNCentralityEstimators) fCentralityPercentile[icent] = value;}; // setter
  void SetReferenceMultiplicity(Int_t value) {fReferenceMultiplicity = value;}; // setter
  void SetMCMultiplicity(Int_t value) {fMCMultiplicity = value;}; // setter
  void SetTimeZeroTOF(Int_t i, Float_t value) {fTimeZeroTOF[i] = value;}; // setter
  void SetTimeZeroTOFSigma(Int_t i, Float_t value) {fTimeZeroTOFSigma[i] = value;}; // setter
  void SetTimeZeroT0(Int_t i, Float_t value) {fTimeZeroT0[i] = value;}; // setter
  void SetMCTimeZero(Float_t value) {fMCTimeZero = value;}; // setter
  
  void Reset(); // reset
  Bool_t CheckLimits(Float_t value, Float_t *params, Float_t nSigma = 3.) const; // check limits

  Bool_t AcceptEvent(Int_t type = 0) const; // accept event proton-proton
  Bool_t AcceptVertex() const; // accept vertex
  Bool_t HasTimeZeroT0_AND() const; // has time-zero T0-AND
  Bool_t HasTimeZeroT0_XOR() const; // has time-zero T0-XOR
  Bool_t HasTimeZeroT0_OR() const; // has time-zero T0-OR
  Bool_t HasTimeZeroTOF(Float_t momentum) const; // has time-zero TOF
  Bool_t HasTimeZeroBest(Float_t momentum) const; // has time-zero TOF
  Bool_t HasTimeZeroSafe(Float_t momentum) const; // has time-zero safe

  Float_t GetTimeZeroT0_AND() const; // get time-zero T0-AND
  Float_t GetTimeZeroT0_XOR() const; // get time-zero T0-XOR
  Float_t GetTimeZeroT0_OR() const; // get time-zero T0-OR
  Float_t GetTimeZeroTOF(Float_t momentum) const; // get time-zero TOF
  Float_t GetTimeZeroBest(Float_t momentum) const; // get time-zero best
  Float_t GetTimeZeroSafe(Float_t momentum) const; // get time-zero safe

  Float_t GetTimeZeroT0Sigma_AND() const; // get time-zero T0-AND sigma
  Float_t GetTimeZeroT0Sigma_XOR() const; // get time-zero T0-XOR sigma
  Float_t GetTimeZeroT0Sigma_OR() const; // get time-zero T0-OR sigma
  Float_t GetTimeZeroTOFSigma(Float_t momentum) const; // get time-zero TOF sigma
  Float_t GetTimeZeroBestSigma(Float_t momentum) const; // get time-zero best sigma
  Float_t GetTimeZeroSafeSigma(Float_t momentum) const; // get time-zero safe sigma

  enum ECentralityEstimator_t {
    kCentEst_V0M, /* V0 multiplicity */
    kCentEst_V0A, /* V0A multiplicity */
    kCentEst_V0C, /* V0C multiplicity */
    kCentEst_TRK, /* N. of tracks */
    kCentEst_TKL, /* N. of tracklets */
    kCentEst_CL1, /*  N. of clusters in layer 1 */
    kCentEst_ZNA, /* ZNA */
    kNCentralityEstimators
  };
  static const Char_t *fgkCentralityEstimatorName[kNCentralityEstimators]; // centrality estimator name

  static void SetVertexZCuts(Float_t min, Float_t max) {fgVertexZ_cuts[0] = min; fgVertexZ_cuts[1] = max;}; // setter

  void ApplyTimeZeroTOFCorrection();
  Double_t GetTimeZeroTOFCorrection() {return fgTimeZeroTOFCentCorrFunc->Eval(fCentralityPercentile[kCentEst_V0M]);};

 private:

  /*** global event info ***/
  Bool_t fIsCollisionCandidate; // is collision candidate
  UInt_t fIsEventSelected; // is event selected
  Bool_t fIsPileupFromSPD; // is pile-up from SPD
  Bool_t fHasVertex; // has vertex
  Float_t fVertexZ; // vertex z
  UChar_t fCentralityQuality; // centrality quality
  Float_t fCentralityPercentile[kNCentralityEstimators]; // centrality percentile
  Int_t fReferenceMultiplicity; // reference multiplicity
  Int_t fMCMultiplicity; // MC multiplicity
  /*** TPC event info ***/
  /*** TOF event info ***/
  Float_t fTimeZeroTOF[10]; // time-zero TOF
  Float_t fTimeZeroTOFSigma[10]; // time-zero TOF sigma
  /*** T0 event info ***/
  Float_t fTimeZeroT0[3]; // time-zero T0
  /*** MC info ***/
  Float_t fMCTimeZero; // MC time-zero

  /*** tools ***/
  static AliTOFPIDResponse fgTOFResponse; // TOF PID response

  /*** cuts ***/
  static Float_t fgVertexZ_cuts[2]; // min,max
  static Float_t fgTimeZeroT0_AND_params[2]; // mean,sigma
  static Float_t fgTimeZeroT0_A_params[2];
  static Float_t fgTimeZeroT0_C_params[2];
  static Float_t fgTimeZeroT0_ACdiff_params[2];
  static Float_t fgTimeZero_TOFT0diff_params[2];

  /*** other ***/
  static Float_t fgTimeZeroSpread;
  static Float_t fgTimeZeroT0_AND_sigma;
  static Float_t fgTimeZeroT0_A_sigma;
  static Float_t fgTimeZeroT0_C_sigma;

  /*** corrections ***/
  static Bool_t fgInitCorrections;
  static const Char_t *fgTimeZeroTOFCentCorrFormula;
  static Double_t fgTimeZeroTOFCentCorrParams[3];
  static TF1 *fgTimeZeroTOFCentCorrFunc;

  ClassDef(AliAnalysisEvent, 6);
};

#endif /* ALIANALYSISEVENT_H */
