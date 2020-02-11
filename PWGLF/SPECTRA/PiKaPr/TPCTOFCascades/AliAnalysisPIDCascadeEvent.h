#ifndef ALIANALYSISPIDCASCADEVENT_H
#define ALIANALYSISPIDCASCADEVENT_H

#include "TObject.h"
#include "AliTOFPIDResponse.h"
#include "TF1.h"

class AliAnalysisPIDCascadeEvent :
public TObject
{

 public:

  AliAnalysisPIDCascadeEvent(); // default constructor
  AliAnalysisPIDCascadeEvent(const AliAnalysisPIDCascadeEvent &source); // copy constructor
  AliAnalysisPIDCascadeEvent &operator=(const AliAnalysisPIDCascadeEvent &source); // operator=
  virtual ~AliAnalysisPIDCascadeEvent(); // default destructor
  Bool_t IsCollisionCandidate() const {return fIsCollisionCandidate;}; // getter
  Bool_t HasVertex() const {return fHasVertex;}; // getter
  Float_t GetVertexZ() const {return fVertexZ;}; // getter
  UChar_t GetCentralityQuality() const {return fCentralityQuality;}; // getter
  Int_t GetReferenceMultiplicity()  {return fReferenceMultiplicity;}; // getter

  void SetIsCollisionCandidate(Bool_t value) {fIsCollisionCandidate = value;}; // setter
  void SetIsEventSelected(UInt_t value) {fIsEventSelected = value;}; // setter
  void SetIsPileupFromSPD(Bool_t value) {fIsPileupFromSPD = value;}; // setter
  void SetHasVertex(Bool_t value) {fHasVertex = value;}; // setter
  void SetVertexZ(Float_t value) {fVertexZ = value;}; // setter
  void SetCentralityQuality(UChar_t value) {fCentralityQuality = value;}; // setter
  void SetReferenceMultiplicity(Int_t value) {fReferenceMultiplicity = value; };// setter
  void SetMCMultiplicity(Int_t value) {fMCMultiplicity = value;}; // setter

  void SetV0Mmultiplicity(Float_t multi) { fV0Mmultiplicity = multi;}; // setter
  void SetRefMult08(Float_t multi) { fRefMult08 = multi;}; // setter
  void SetRefMult05(Float_t multi) { fRefMult05 = multi;}; // setter
  void SetSPDTracklets(Float_t multi) { fSPDTracklets = multi;}; // setter

  void SetEventFlags(Int_t NewFlag) { fEventFlags = NewFlag; };
  Int_t GetEventFlags() { return fEventFlags; };
  Bool_t HasEventFlag(Int_t CheckFlag) { return (fEventFlags&CheckFlag)==CheckFlag; };
  void SetMagneticField(Double_t MagneticFieldValue) { fMagneticField = MagneticFieldValue; };
  void SetRunNumber(Int_t RunNo) { fRunNo = RunNo; };


  void Reset(); // reset

  Bool_t AcceptEvent(Bool_t CheckVertex=kTRUE, Int_t type = 0) const; // accept event proton-proton
  Bool_t AcceptVertex() const; // accept vertex

  Float_t GetV0Mmultiplicity() { return fV0Mmultiplicity;}; //getter
  Float_t GetRefMult08() { return fRefMult08;};
  Float_t GetRefMult05() { return fRefMult05;};
  Float_t GetSPDTracklets() { return fSPDTracklets;};

  Int_t GetMCMultiplicity() {return fMCMultiplicity; }; // getter
  Bool_t IsPileup() {return fIsPileupFromSPD; }; // getter
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
    kVertexSelected2015pp=64,
    kSPDandTrkVtxExists=128,
    kPassProximityCut=256,
    kAll = 511
  };
  static const Char_t *fgkCentralityEstimatorName[kNCentralityEstimators]; // centrality estimator name
  static void SetVertexZCuts(Float_t min, Float_t max) {fgVertexZ_cuts[0] = min; fgVertexZ_cuts[1] = max;}; // setter
  static void SetCheckFlag(Int_t);
  static void AddCheckFlag(EventFlags_t);
  static void RemoveCheckFlag(EventFlags_t);
  static Int_t GetCheckFlag() { return fgFlagToCheck; };
  Bool_t CheckFlag();
  static void PrintEventSelection();
  /*** tools ***/
  static AliTOFPIDResponse fgTOFResponse; // TOF PID response

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
  Float_t fRefMult08;
  Float_t fRefMult05;
  Float_t fSPDTracklets;
  //  Float_t fV0CellAmplitude[64];
  Int_t fMCMultiplicity; // MC multiplicity
  /*** TPC event info ***/
  /*** TOF event info ***/
  /*** T0 event info ***/
  /*** MC info ***/
  //Some flags from PPVsMultUtils class
  Int_t fEventFlags;
  Double_t fMagneticField;
  Int_t fRunNo;



  /*** cuts ***/
  static Float_t fgVertexZ_cuts[2]; // min,max

  /*** other ***/

  /*** corrections ***/
  static Int_t fgFlagToCheck; //! check flag

  ClassDef(AliAnalysisPIDCascadeEvent, 2);
};

#endif /* ALIANALYSISPIDCASCADEEVENT_H */
