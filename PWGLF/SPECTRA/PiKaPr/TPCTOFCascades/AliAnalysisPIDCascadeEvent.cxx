#include "AliAnalysisPIDCascadeEvent.h"
#include "AliVEvent.h"

ClassImp(AliAnalysisPIDCascadeEvent)

//___________________________________________________________

//___________________________________________________________

AliTOFPIDResponse AliAnalysisPIDCascadeEvent::fgTOFResponse;

Float_t AliAnalysisPIDCascadeEvent::fgVertexZ_cuts[2] = {-10., 10.}; /* min,max */

const Char_t *AliAnalysisPIDCascadeEvent::fgkCentralityEstimatorName[AliAnalysisPIDCascadeEvent::kNCentralityEstimators] = {
  "V0M",
  "V0A",
  "V0C",
  "TRK",
  "TKL",
  "CL1",
  "ZNA"
};

Int_t AliAnalysisPIDCascadeEvent::fgFlagToCheck = 205;

//___________________________________________________________

AliAnalysisPIDCascadeEvent::AliAnalysisPIDCascadeEvent() :
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
  fEventFlags(0),
  fMagneticField(0.),
  fRunNo(0)
{
  /*
   * default constructor
   */
}

//___________________________________________________________

AliAnalysisPIDCascadeEvent::AliAnalysisPIDCascadeEvent(const AliAnalysisPIDCascadeEvent &source) :
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
  fEventFlags(source.fEventFlags),
  fMagneticField(source.fMagneticField),
  fRunNo(source.fRunNo)
{
  /*
   * copy constructor
   */

    fReferenceMultiplicity = source.fReferenceMultiplicity;
    fV0Mmultiplicity = source.fV0Mmultiplicity;
}

//___________________________________________________________

AliAnalysisPIDCascadeEvent &
AliAnalysisPIDCascadeEvent::operator=(const AliAnalysisPIDCascadeEvent &source)
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
  fEventFlags = source.fEventFlags;
  fMagneticField=source.fMagneticField;
  fRunNo=source.fRunNo;
  return *this;
}

//___________________________________________________________

AliAnalysisPIDCascadeEvent::~AliAnalysisPIDCascadeEvent()
{
  /*
   * default destructor
   */

}

//___________________________________________________________


void
AliAnalysisPIDCascadeEvent::Reset()
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
}

//___________________________________________________________

//___________________________________________________________

Bool_t
AliAnalysisPIDCascadeEvent::AcceptEvent(Bool_t CheckVertex, Int_t type) const
{
  /*
   * accept event proton-proton
   */

  if (!fIsCollisionCandidate) return kFALSE;
  if (fCentralityQuality != 0) return kFALSE;
  if (!(fIsEventSelected & AliVEvent::kINT7)) return kFALSE;
  if((fEventFlags&fgFlagToCheck)!=fgFlagToCheck) return kFALSE;
  if(CheckVertex)
    if(!AcceptVertex()) return kFALSE;
  if (type > 0) {
    if (fIsPileupFromSPD) return kFALSE;
    if (!(fIsEventSelected & AliVEvent::kMB)) return kFALSE;
  };


  return kTRUE;
}

//___________________________________________________________

Bool_t
AliAnalysisPIDCascadeEvent::AcceptVertex() const
{
  /*
   * accept vertex
   */

  if (!HasVertex()) return kFALSE;
  if (fVertexZ < fgVertexZ_cuts[0] || fVertexZ > fgVertexZ_cuts[1]) return kFALSE;
  return kTRUE;
}

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________

//___________________________________________________________
void AliAnalysisPIDCascadeEvent::SetCheckFlag(Int_t newval) {
  if(newval>kAll||newval<0) {
    printf("Flag value %i not defined!\n",newval);
    return;
  };
  fgFlagToCheck = newval;
};
void AliAnalysisPIDCascadeEvent::AddCheckFlag(EventFlags_t av) {
  fgFlagToCheck = fgFlagToCheck|av;
};
void AliAnalysisPIDCascadeEvent::RemoveCheckFlag(EventFlags_t av) {
  Int_t flagmask = kAll-av;
  fgFlagToCheck = fgFlagToCheck&flagmask;
};
Bool_t AliAnalysisPIDCascadeEvent::CheckFlag() {
  return (fEventFlags&fgFlagToCheck)==fgFlagToCheck;
};
void AliAnalysisPIDCascadeEvent::PrintEventSelection() {
  printf("AliAnalysisPIDCascadeEvent::AcceptEvent() requires:\n");
  printf("Vertex position: %f..%f\n",fgVertexZ_cuts[0],fgVertexZ_cuts[1]);
  printf("Trigger class: kINT7\n");
  printf("Not pileup in SPD:   %s\n",(fgFlagToCheck&kNotPileupInSPD)?"Yes":"No");
  printf("Not pileup in MV:    %s\n",(fgFlagToCheck&kNotPileupInMV)?"Yes":"No");
  printf("Not pileup in MB:    %s\n",(fgFlagToCheck&kNotPileupInMB)?"Yes":"No");
  printf("INEL > 0:            %s\n",(fgFlagToCheck&kINELgtZERO)?"Yes":"No");
  printf("No inconsistent VTX: %s\n",(fgFlagToCheck&kNoInconsistentVtx)?"Yes":"No");
  printf("No assym. in V0:     %s\n",(fgFlagToCheck&kNoV0Asym)?"Yes":"No");
  printf("2015 pp vertex cut:  %s\n",(fgFlagToCheck&kVertexSelected2015pp)?"Yes":"No");
  printf("Req. SPD & TRK vtx.: %s\n",(fgFlagToCheck&kSPDandTrkVtxExists)?"Yes":"No");
  printf("Check proximity cut: %s\n",(fgFlagToCheck&kPassProximityCut)?"Yes":"No");
};
