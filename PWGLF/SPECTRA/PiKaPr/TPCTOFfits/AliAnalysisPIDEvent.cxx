#include "AliAnalysisPIDEvent.h"
#include "AliVEvent.h"

ClassImp(AliAnalysisPIDEvent)

//___________________________________________________________

Bool_t AliAnalysisPIDEvent::fgInitCorrections = kFALSE;

const Char_t *AliAnalysisPIDEvent::fgTimeZeroTOFCentCorrFormula = "[0]+[1]*x+[2]*x*x";

Double_t AliAnalysisPIDEvent::fgTimeZeroTOFCentCorrParams[3] = {
  -14.9262, 0.0212941, 0.000960528
};

TF1 *AliAnalysisPIDEvent::fgTimeZeroTOFCentCorrFunc = NULL;

//___________________________________________________________

AliTOFPIDResponse AliAnalysisPIDEvent::fgTOFResponse;

Float_t AliAnalysisPIDEvent::fgVertexZ_cuts[2] = {-10., 10.}; /* min,max */
Float_t AliAnalysisPIDEvent::fgTimeZeroT0_AND_params[2] = {0., 1000.}; /* mean,sigma */
Float_t AliAnalysisPIDEvent::fgTimeZeroT0_A_params[2] = {0., 1000.}; /* mean,sigma */
Float_t AliAnalysisPIDEvent::fgTimeZeroT0_C_params[2] = {0., 1000.}; /* mean,sigma */
Float_t AliAnalysisPIDEvent::fgTimeZeroT0_ACdiff_params[2] = {0., 1000.}; /* mean,sigma */
Float_t AliAnalysisPIDEvent::fgTimeZero_TOFT0diff_params[2] = {0., 1000.}; /* mean,sigma */

const Char_t *AliAnalysisPIDEvent::fgkCentralityEstimatorName[AliAnalysisPIDEvent::kNCentralityEstimators] = {
  "V0M",
  "V0A",
  "V0C",
  "TRK",
  "TKL",
  "CL1",
  "ZNA"
};

Float_t AliAnalysisPIDEvent::fgTimeZeroSpread = 196.7;
Float_t AliAnalysisPIDEvent::fgTimeZeroT0_AND_sigma = 3.87264325235363032e+01;
Float_t AliAnalysisPIDEvent::fgTimeZeroT0_A_sigma = 8.27180042372880706e+01;
Float_t AliAnalysisPIDEvent::fgTimeZeroT0_C_sigma = 9.73209262235003933e+01;
Int_t AliAnalysisPIDEvent::fgFlagToCheck = 63;

//___________________________________________________________

AliAnalysisPIDEvent::AliAnalysisPIDEvent() :
  TObject(),
  fIsCollisionCandidate(kFALSE),
  fIsEventSelected(0),
  fIsPileupFromSPD(kFALSE),
  fHasVertex(kFALSE),
  fVertexZ(-999.0),
  fCentralityQuality(0),
  fReferenceMultiplicity(0),
  fV0Mmultiplicity(0.),
  fMCMultiplicity(0),
  fTimeZeroTOF(),
  fTimeZeroTOFSigma(),
  fTimeZeroT0(),
  fMCTimeZero(0.),
  fEventFlags(0),
  fMagneticField(0.),
  fRunNo(0)
{
  /*
   * default constructor
   */

  /* init corrections */
  if (!fgInitCorrections) {
    
    /*** TIME-ZERO TOF CENTRALITY CORRECTION ***/
    fgTimeZeroTOFCentCorrFunc = new TF1("fTimeZeroTOFCentCorrFunc", fgTimeZeroTOFCentCorrFormula, 0., 100.);
    for (Int_t iparam = 0; iparam < 3; iparam++)
      fgTimeZeroTOFCentCorrFunc->SetParameter(iparam, fgTimeZeroTOFCentCorrParams[iparam]);
    
    /* set init flag */
    fgInitCorrections = kTRUE;
  }

  /* reset */
  Reset();
}

//___________________________________________________________

AliAnalysisPIDEvent::AliAnalysisPIDEvent(const AliAnalysisPIDEvent &source) :
  TObject(source),
  fIsCollisionCandidate(source.fIsCollisionCandidate),
  fIsEventSelected(source.fIsEventSelected),
  fIsPileupFromSPD(source.fIsPileupFromSPD),
  fHasVertex(source.fHasVertex),
  fVertexZ(source.fVertexZ),
  fCentralityQuality(source.fCentralityQuality),
  fReferenceMultiplicity(0),
  fV0Mmultiplicity(source.fV0Mmultiplicity),
  fMCMultiplicity(source.fMCMultiplicity),
  fTimeZeroTOF(),
  fTimeZeroTOFSigma(),
  fTimeZeroT0(),
  fMCTimeZero(source.fMCTimeZero),
  fEventFlags(source.fEventFlags),
  fMagneticField(source.fMagneticField),
  fRunNo(source.fRunNo)
{
  /*
   * copy constructor
   */

  for (Int_t i = 0; i < 10; i++) {
    fTimeZeroTOF[i] = source.fTimeZeroTOF[i];
    fTimeZeroTOFSigma[i] = source.fTimeZeroTOFSigma[i];
  }
  for (Int_t i = 0; i < 3; i++) 
    fTimeZeroT0[i] = source.fTimeZeroT0[i];

    fReferenceMultiplicity = source.fReferenceMultiplicity;
    fV0Mmultiplicity = source.fV0Mmultiplicity;
}

//___________________________________________________________

AliAnalysisPIDEvent &
AliAnalysisPIDEvent::operator=(const AliAnalysisPIDEvent &source)
{
  /*
   * operator=
   */

  if (&source == this) return *this;
  TObject::operator=(source);

  fIsCollisionCandidate = source.fIsCollisionCandidate;
  fIsEventSelected = source.fIsEventSelected;
  fIsPileupFromSPD = source.fIsPileupFromSPD;
  fHasVertex = source.fHasVertex;
  fVertexZ = source.fVertexZ;
  fCentralityQuality = source.fCentralityQuality;
  fReferenceMultiplicity = source.fReferenceMultiplicity;
  fV0Mmultiplicity = source.fV0Mmultiplicity;
  fMCMultiplicity = source.fMCMultiplicity;
  for (Int_t i = 0; i < 10; i++) {
    fTimeZeroTOF[i] = source.fTimeZeroTOF[i];
    fTimeZeroTOFSigma[i] = source.fTimeZeroTOFSigma[i];
  }
  for (Int_t i = 0; i < 3; i++) 
    fTimeZeroT0[i] = source.fTimeZeroT0[i];
  fMCTimeZero = source.fMCTimeZero;
  fEventFlags = source.fEventFlags;
  fMagneticField=source.fMagneticField;
  fRunNo=source.fRunNo;
  return *this;
}

//___________________________________________________________

AliAnalysisPIDEvent::~AliAnalysisPIDEvent()
{
  /*
   * default destructor
   */

}

//___________________________________________________________


void
AliAnalysisPIDEvent::Reset()
{
  /*
   * reset
   */

  fIsCollisionCandidate = kFALSE;
  fIsEventSelected = 0;
  fIsPileupFromSPD = kFALSE;
  fHasVertex = 0.;
  fVertexZ = -999.;
  fCentralityQuality = 0;
  fReferenceMultiplicity = 0;
  fMCMultiplicity = 0;
  fV0Mmultiplicity = 0;
  fEventFlags = 0;
  fMagneticField=0;
  fRunNo=0;
  for (Int_t i = 0; i < 10; i++) {
    fTimeZeroTOF[i] = 0.;
    fTimeZeroTOFSigma[i] = 0.;
  }
  for (Int_t i = 0; i < 3; i++) 
    fTimeZeroT0[i] = 0.;  
  fMCTimeZero = 0.;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDEvent::CheckLimits(Float_t value, Float_t *params, Float_t nSigma) const
{
  /*
   * check limits
   */

  Float_t min = params[0] - nSigma * params[1];
  Float_t max = params[0] + nSigma * params[1];
  if (value < min || value > max) return kFALSE;
  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDEvent::AcceptEvent(Int_t type) const
{
  /*
   * accept event proton-proton
   */

  if (!fIsCollisionCandidate) return kFALSE;
  if (!AcceptVertex()) return kFALSE;
  if (fCentralityQuality != 0) return kFALSE;
  if (!(fIsEventSelected & AliVEvent::kINT7)) return kFALSE;
  if((fEventFlags&fgFlagToCheck)!=fgFlagToCheck) return kFALSE;
  if (type > 0) {
    if (fIsPileupFromSPD) return kFALSE;
    if (!(fIsEventSelected & AliVEvent::kMB)) return kFALSE;
  };


  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDEvent::AcceptVertex() const
{
  /*
   * accept vertex
   */

  if (!HasVertex()) return kFALSE;
  if (fVertexZ < fgVertexZ_cuts[0] || fVertexZ > fgVertexZ_cuts[1]) return kFALSE;
  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDEvent::HasTimeZeroT0_AND() const
{
  /*
   * has time-zero T0-AND
   */

  /* check T0-AND */
  if (!CheckLimits(fTimeZeroT0[0], fgTimeZeroT0_AND_params)) return kFALSE;
  /* check T0-A */
  if (!CheckLimits(fTimeZeroT0[1], fgTimeZeroT0_A_params)) return kFALSE;
  /* check T0-C */
  if (!CheckLimits(fTimeZeroT0[2], fgTimeZeroT0_C_params)) return kFALSE;
  /* check A-C difference */
  if (!CheckLimits(fTimeZeroT0[1] - fTimeZeroT0[2], fgTimeZeroT0_ACdiff_params)) return kFALSE;
  /* ok */
  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDEvent::HasTimeZeroT0_XOR() const
{
  /*
   * has time-zero T0-XOR
   */

  /* check T0-A only */
  if (CheckLimits(fTimeZeroT0[1], fgTimeZeroT0_A_params) &&
      !CheckLimits(fTimeZeroT0[2], fgTimeZeroT0_C_params)) return kTRUE;
  /* check T0-C only */
  if (!CheckLimits(fTimeZeroT0[1], fgTimeZeroT0_A_params) &&
      CheckLimits(fTimeZeroT0[2], fgTimeZeroT0_C_params)) return kTRUE;
  return kFALSE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDEvent::HasTimeZeroT0_OR() const
{
  /*
   * has time-zero T0-OR
   */

  if (HasTimeZeroT0_AND() || HasTimeZeroT0_XOR()) return kTRUE;
  return kFALSE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDEvent::HasTimeZeroTOF(Float_t momentum) const
{
  /*
   * has time-zero TOF
   */

  Int_t momBin = fgTOFResponse.GetMomBin(momentum);
  if (fTimeZeroTOF[momBin] == 0.) return kFALSE;
  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDEvent::HasTimeZeroBest(Float_t momentum) const
{
  /*
   * has time-zero best
   */

  if (HasTimeZeroT0_OR() || HasTimeZeroTOF(momentum)) return kTRUE;
  return kFALSE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDEvent::HasTimeZeroSafe(Float_t momentum) const
{
  /*
   * has time-zero safe
   */

  if (!HasTimeZeroT0_AND() || !HasTimeZeroTOF(momentum)) return kFALSE;
  Float_t tzeroTOF = GetTimeZeroTOF(momentum);
  Float_t tzeroT0 = GetTimeZeroT0_AND();
  Float_t diff = tzeroTOF - tzeroT0;
  if (!CheckLimits(diff, fgTimeZero_TOFT0diff_params)) return kFALSE;
  return kTRUE;
}

//___________________________________________________________

Float_t
AliAnalysisPIDEvent::GetTimeZeroT0_AND() const
{
  /*
   * get time-zero T0-AND
   */

  return fTimeZeroT0[0] - fgTimeZeroT0_AND_params[0];
}

//___________________________________________________________

Float_t
AliAnalysisPIDEvent::GetTimeZeroT0Sigma_AND() const
{
  /*
   * get time-zero T0-AND sigma
   */

  return fgTimeZeroT0_AND_sigma;
}

//___________________________________________________________

Float_t
AliAnalysisPIDEvent::GetTimeZeroT0_XOR() const
{
  /*
   * get time-zero T0-XOR
   */
  
  if (CheckLimits(fTimeZeroT0[1], fgTimeZeroT0_A_params) &&
      !CheckLimits(fTimeZeroT0[2], fgTimeZeroT0_C_params)) return fTimeZeroT0[1] - fgTimeZeroT0_A_params[0];
  if (!CheckLimits(fTimeZeroT0[1], fgTimeZeroT0_A_params) &&
      CheckLimits(fTimeZeroT0[2], fgTimeZeroT0_C_params)) return fTimeZeroT0[2] - fgTimeZeroT0_C_params[0];
  return 0.;
}

//___________________________________________________________

Float_t
AliAnalysisPIDEvent::GetTimeZeroT0Sigma_XOR() const
{
  /*
   * get time-zero T0-XOR sigma
   */

  if (CheckLimits(fTimeZeroT0[1], fgTimeZeroT0_A_params) &&
      !CheckLimits(fTimeZeroT0[2], fgTimeZeroT0_C_params)) return fgTimeZeroT0_A_sigma;
  if (!CheckLimits(fTimeZeroT0[1], fgTimeZeroT0_A_params) &&
      CheckLimits(fTimeZeroT0[2], fgTimeZeroT0_C_params)) return fgTimeZeroT0_C_sigma;
  return 0.;
}

//___________________________________________________________

Float_t
AliAnalysisPIDEvent::GetTimeZeroT0_OR() const
{
  /*
   * get time-zero T0-OR
   */
  
  if (HasTimeZeroT0_AND()) return GetTimeZeroT0_AND();
  if (HasTimeZeroT0_XOR()) return GetTimeZeroT0_XOR();
  return 0.;
}

//___________________________________________________________

Float_t
AliAnalysisPIDEvent::GetTimeZeroT0Sigma_OR() const
{
  /*
   * get time-zero T0-OR sigma
   */
  
  if (HasTimeZeroT0_AND()) return GetTimeZeroT0Sigma_AND();
  if (HasTimeZeroT0_XOR()) return GetTimeZeroT0Sigma_XOR();
  return 0.;
}

//___________________________________________________________

Float_t
AliAnalysisPIDEvent::GetTimeZeroTOF(Float_t momentum) const
{
  /*
   * get time-zero TOF
   */

  Int_t momBin = fgTOFResponse.GetMomBin(momentum);
  return fTimeZeroTOF[momBin];
}

//___________________________________________________________

Float_t
AliAnalysisPIDEvent::GetTimeZeroTOFSigma(Float_t momentum) const
{
  /*
   * get time-zero TOF sigma
   */

  Int_t momBin = fgTOFResponse.GetMomBin(momentum);
  return fTimeZeroTOFSigma[momBin];
}

//___________________________________________________________

Float_t
AliAnalysisPIDEvent::GetTimeZeroBest(Float_t momentum) const
{
  /*
   * get time-zero best
   */

  if (HasTimeZeroTOF(momentum)) return GetTimeZeroTOF(momentum);
  if (HasTimeZeroT0_OR()) return GetTimeZeroT0_OR();
  return 0.;
}

//___________________________________________________________

Float_t
AliAnalysisPIDEvent::GetTimeZeroBestSigma(Float_t momentum) const
{
  /*
   * get time-zero best sigma
   */

  if (HasTimeZeroTOF(momentum)) return GetTimeZeroTOFSigma(momentum);
  if (HasTimeZeroT0_OR()) return GetTimeZeroT0Sigma_OR();
  return fgTimeZeroSpread;
}

//___________________________________________________________

Float_t
AliAnalysisPIDEvent::GetTimeZeroSafe(Float_t momentum) const
{
  /*
   * get time-zero safe
   */

  return GetTimeZeroTOF(momentum);
}

//___________________________________________________________

Float_t
AliAnalysisPIDEvent::GetTimeZeroSafeSigma(Float_t momentum) const
{
  /*
   * get time-zero safe sigma
   */

  return GetTimeZeroTOFSigma(momentum);
}

//___________________________________________________________
void AliAnalysisPIDEvent::SetCheckFlag(Int_t newval) {
  if(newval>63||newval<0) {
    printf("Flag value %i not defined!\n",newval);
    return;
  };
  fgFlagToCheck = newval;
};
void AliAnalysisPIDEvent::AddCheckFlag(EventFlags_t av) {
  fgFlagToCheck = fgFlagToCheck|av;
};
void AliAnalysisPIDEvent::RemoveCheckFlag(EventFlags_t av) {
  Int_t flagmask = kAll-av;
  fgFlagToCheck = fgFlagToCheck&flagmask;
};
Bool_t AliAnalysisPIDEvent::CheckFlag() {
  return (fEventFlags&fgFlagToCheck)==fgFlagToCheck;
};
void AliAnalysisPIDEvent::PrintEventSelection() {
  printf("AliAnalysisPIDEvent::AcceptEvent() requires:\n");
  printf("Vertex position: %f..%f\n",fgVertexZ_cuts[0],fgVertexZ_cuts[1]);
  printf("Trigger class: kINT7\n");
  printf("Not pileup in SPD:   %s\n",(fgFlagToCheck&kNotPileupInSPD)?"Yes":"No");
  printf("Not pileup in MV:    %s\n",(fgFlagToCheck&kNotPileupInMV)?"Yes":"No");
  printf("Not pileup in MB:    %s\n",(fgFlagToCheck&kNotPileupInMB)?"Yes":"No");
  printf("INEL > 0:            %s\n",(fgFlagToCheck&kINELgtZERO)?"Yes":"No");
  printf("No inconsistent VTX: %s\n",(fgFlagToCheck&kNoInconsistentVtx)?"Yes":"No");
  printf("No asynn. in V0:    %s\n",(fgFlagToCheck&kNoV0Asym)?"Yes":"No");
};
