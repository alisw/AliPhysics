#ifndef ALIANALYSISPIDEVENT_H
#define ALIANALYSISPIDEVENT_H

#include "TObject.h"
#include "AliTOFPIDResponse.h"
#include "TF1.h"

class AliAnalysisPIDEvent :
public TObject
{

 public:

  AliAnalysisPIDEvent(); // default constructor
  AliAnalysisPIDEvent(const AliAnalysisPIDEvent &source); // copy constructor
  AliAnalysisPIDEvent &operator=(const AliAnalysisPIDEvent &source); // operator=
  virtual ~AliAnalysisPIDEvent(); // default destructor
  Bool_t IsCollisionCandidate() const {return fIsCollisionCandidate;}; // getter
  Bool_t HasVertex() const {return fHasVertex;}; // getter
  Float_t GetVertexZ() const {return fVertexZ;}; // getter
  UChar_t GetCentralityQuality() const {return fCentralityQuality;}; // getter
  Int_t GetReferenceMultiplicity()  {return fReferenceMultiplicity;}; // getter
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
  void SetReferenceMultiplicity(Int_t value) {fReferenceMultiplicity = value; };// setter
  void SetMCMultiplicity(Int_t value) {fMCMultiplicity = value;}; // setter
  void SetTimeZeroTOF(Int_t i, Float_t value) {fTimeZeroTOF[i] = value;}; // setter
  void SetTimeZeroTOFSigma(Int_t i, Float_t value) {fTimeZeroTOFSigma[i] = value;}; // setter
  void SetTimeZeroT0(Int_t i, Float_t value) {fTimeZeroT0[i] = value;}; // setter
  void SetMCTimeZero(Float_t value) {fMCTimeZero = value;}; // setter
  void SetV0Mmultiplicity(Float_t multi) { fV0Mmultiplicity = multi;}; // setter
  void SetEventFlags(Int_t NewFlag) { fEventFlags = NewFlag; };
  Int_t GetEventFlags() { return fEventFlags; };
  Bool_t HasEventFlag(Int_t CheckFlag) { return (fEventFlags&CheckFlag)==CheckFlag; };
  void SetMagneticField(Double_t MagneticFieldValue) { fMagneticField = MagneticFieldValue; };
  void SetRunNumber(Int_t RunNo) { fRunNo = RunNo; };


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
  Float_t GetV0Mmultiplicity() { return fV0Mmultiplicity;}; //getter
  Int_t GetMCMultiplicity() {return fMCMultiplicity; }; // getter
  Bool_t IsPileup() {return fIsPileupFromSPD; }; // getter
  void SetV0CellAmplitude(Int_t i, Float_t CellAmp) { fV0CellAmplitude[i] = CellAmp; }; // setter
    Float_t *GetV0CellAmplitude() { return fV0CellAmplitude; }; //! Return the signal array in V0 cells. 0-31 for V0A, 32-63 for V0C
  Float_t GetV0CellAmplitude(Int_t CellNo) { return fV0CellAmplitude[CellNo]; };//! Return the signal CellNo V0 cell. 0-31 for V0A, 32-63 for V0C
  Double_t GetMagneticField() { return fMagneticField;};
  Int_t GetRunNumber() { return fRunNo; };
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
  enum EventFlags_t {
    kNotPileupInSPD = 1,
    kNotPileupInMV = 2,
    kNotPileupInMB = 4,
    kINELgtZERO = 8,
    kNoInconsistentVtx = 16,
    kNoV0Asym = 32,
    kAll = 63
  };
  static const Char_t *fgkCentralityEstimatorName[kNCentralityEstimators]; // centrality estimator name
  static void SetVertexZCuts(Float_t min, Float_t max) {fgVertexZ_cuts[0] = min; fgVertexZ_cuts[1] = max;}; // setter
  static void SetCheckFlag(Int_t);
  static void AddCheckFlag(EventFlags_t);
  static void RemoveCheckFlag(EventFlags_t);
  static Int_t GetCheckFlag() { return fgFlagToCheck; }; 
  Bool_t CheckFlag();
  static void PrintEventSelection();

 private:

  /*** global event info ***/
  Bool_t fIsCollisionCandidate; // is collision candidate
  UInt_t fIsEventSelected; // is event selected
  Bool_t fIsPileupFromSPD; // is pile-up from SPD
  Bool_t fHasVertex; // has vertex
  Float_t fVertexZ; // vertex z
  UChar_t fCentralityQuality; // centrality quality
  Int_t fReferenceMultiplicity; // reference multiplicity in eta 0.8
  Float_t fV0Mmultiplicity;
  Float_t fV0CellAmplitude[64];
  Int_t fMCMultiplicity; // MC multiplicity
  /*** TPC event info ***/
  /*** TOF event info ***/
  Float_t fTimeZeroTOF[10]; // time-zero TOF
  Float_t fTimeZeroTOFSigma[10]; // time-zero TOF sigma
  /*** T0 event info ***/
  Float_t fTimeZeroT0[3]; // time-zero T0
  /*** MC info ***/
  Float_t fMCTimeZero; // MC time-zero
  //Some flags from PPVsMultUtils class
  Int_t fEventFlags;
  Double_t fMagneticField;
  Int_t fRunNo;


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
  static Int_t fgFlagToCheck; //! check flag

  ClassDef(AliAnalysisPIDEvent, 5);
};

#endif /* ALIANALYSISPIDEVENT_H */
