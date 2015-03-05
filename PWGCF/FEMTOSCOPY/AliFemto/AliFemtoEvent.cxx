///////////////////////////////////////////////////////////////////////////
//
//  AliFemtoEvent: hold the information specific to the event and a
//  track list
//  AliFemtoEvent is the "transient microDST"  Objects of this class are
//   generated from the input data by a Reader, and then presented to
//   the Cuts of the various active Analyses.
//
////////////////
#include "AliFemtoEvent.h"
#include "AliFemtoTrack.h"
#include "AliFemtoV0.h"
#include "AliFemtoXi.h"
#include "AliFemtoKink.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoXiCut.h"
#include "AliFemtoKinkCut.h"
#include "PhysicalConstants.h"
#include "SystemOfUnits.h"
#include "AliEventplane.h"

#include <algorithm>

// Mike removed all of the AliFemtoTTree stuff here 21apr2006 - it was not used for a long time.

//___________________
AliFemtoEvent::AliFemtoEvent():
  fEventNumber(0),
  fRunNumber(0),
  fNumberOfTracks(0),
  fNormalizedMult(-2),
  fSPDMult(0),
  fEstimateITSTPC(0),
  fEstimateTracklets(0),
  fEstimateITSPure(0),
  fCentralityV0(0.0f),
  fCentralityV0A(0.0f),
  fCentralityV0C(0.0f),
  fCentralityZNA(0.0f),
  fCentralityZNC(0.0f),
  fCentralityCL1(0.0f),
  fCentralityCL0(0.0f),
  fCentralityTKL(0.0f),
  fCentralityFMD(0.0f),
  fCentralityTrk(0.0f),
  fCentralityCND(0.0f),
  fCentralityNPA(0.0f),
  fCentralitySPD1(0.0f),
  fMagneticField(0.0),
  fIsCollisionCandidate(kTRUE),
  fPrimVertPos(-999.0,-999.0,-999.0),
  fPrimVertCov(),
  fTrackCollection(NULL),
  fV0Collection(NULL),
  fXiCollection(NULL),
  fKinkCollection(NULL),
  fZDCN1Energy(0.0f),
  fZDCP1Energy(0.0f),
  fZDCN2Energy(0.0f),
  fZDCP2Energy(0.0f),
  fZDCEMEnergy(0.0f),
  fZDCParticipants(0),
  fTriggerMask(0),
  fTriggerCluster(0),
  fReactionPlaneAngle(0.0f),
  fEP(NULL)
{
  // Default constructor

  // fill all 6 entries of fPrimVertCov with 0.000000000001
  std::fill_n(fPrimVertCov, 6, 1.0e-12);
  fTrackCollection = new AliFemtoTrackCollection;
  fV0Collection = new AliFemtoV0Collection;
  fXiCollection = new AliFemtoXiCollection;
  fKinkCollection = new AliFemtoKinkCollection;
}
//___________________
AliFemtoEvent::AliFemtoEvent(const AliFemtoEvent& ev, AliFemtoTrackCut* tCut, AliFemtoV0Cut* vCut, AliFemtoXiCut* xCut, AliFemtoKinkCut* kCut):
  fEventNumber(ev.fEventNumber),
  fRunNumber(ev.fRunNumber),
  fNumberOfTracks(ev.fNumberOfTracks),
  fNormalizedMult(ev.fNormalizedMult),
  fSPDMult(ev.fSPDMult),
  fEstimateITSTPC(ev.fEstimateITSTPC),
  fEstimateTracklets(ev.fEstimateTracklets),
  fEstimateITSPure(ev.fEstimateITSPure),
  fCentralityV0(ev.fCentralityV0),
  fCentralityV0A(ev.fCentralityV0A),
  fCentralityV0C(ev.fCentralityV0C),
  fCentralityZNA(ev.fCentralityZNA),
  fCentralityZNC(ev.fCentralityZNC),
  fCentralityCL1(ev.fCentralityCL1),
  fCentralityCL0(ev.fCentralityCL0),
  fCentralityTKL(ev.fCentralityTKL),
  fCentralityFMD(ev.fCentralityFMD),
  fCentralityTrk(ev.fCentralityTrk),
  fCentralityCND(ev.fCentralityCND),
  fCentralityNPA(ev.fCentralityNPA),
  fCentralitySPD1(ev.fCentralitySPD1),
  fMagneticField(ev.fMagneticField),
  fIsCollisionCandidate(ev.fIsCollisionCandidate),
  fPrimVertPos(ev.fPrimVertPos),
  fPrimVertCov(),
  fTrackCollection(NULL),
  fV0Collection(NULL),
  fXiCollection(NULL),
  fKinkCollection(NULL),
  fZDCN1Energy(ev.fZDCN1Energy),
  fZDCP1Energy(ev.fZDCP1Energy),
  fZDCN2Energy(ev.fZDCN2Energy),
  fZDCP2Energy(ev.fZDCP2Energy),
  fZDCEMEnergy(ev.fZDCEMEnergy),
  fZDCParticipants(ev.fZDCParticipants),
  fTriggerMask(ev.fTriggerMask),
  fTriggerCluster(ev.fTriggerCluster),
  fReactionPlaneAngle(ev.fReactionPlaneAngle),
  fEP(ev.fEP)
{ // copy constructor with track and v0 cuts
  //cout << "AliFemtoEvent::AliFemtoEvent(const AliFemtoEvent& ev, AliFemtoTrackCut* tCut, AliFemtoV0Cut* vCut, AliFemtoV0Cut* kCut)" << endl;

  SetPrimVertCov(ev.PrimVertCov());

  // create collections
  fTrackCollection = new AliFemtoTrackCollection;
  fV0Collection = new AliFemtoV0Collection;
  fXiCollection = new AliFemtoXiCollection;
  fKinkCollection = new AliFemtoKinkCollection;
  // copy track collection
  for ( AliFemtoTrackIterator tIter=ev.fTrackCollection->begin(); tIter!=ev.fTrackCollection->end(); tIter++) {
    if ( !tCut || tCut->Pass(*tIter) ) {
      AliFemtoTrack* trackCopy = new AliFemtoTrack(**tIter);
      fTrackCollection->push_back(trackCopy);
    }
  }
  // copy v0 collection
  for ( AliFemtoV0Iterator vIter=ev.fV0Collection->begin(); vIter!=ev.fV0Collection->end(); vIter++) {
    if ( !vCut || vCut->Pass(*vIter) ) {
      AliFemtoV0* v0Copy = new AliFemtoV0(**vIter);
      fV0Collection->push_back(v0Copy);
    }
  }
  // copy xi collection
  for ( AliFemtoXiIterator xIter=ev.fXiCollection->begin(); xIter!=ev.fXiCollection->end(); xIter++) {
    if ( !xCut || xCut->Pass(*xIter) ) {
      AliFemtoXi* xiCopy = new AliFemtoXi(**xIter);
      fXiCollection->push_back(xiCopy);
    }
  }
  // copy kink collection
  for ( AliFemtoKinkIterator kIter=ev.fKinkCollection->begin(); kIter!=ev.fKinkCollection->end(); kIter++) {
    if ( !kCut || kCut->Pass(*kIter) ) {
      //cout << " kinkCut passed " << endl;
      AliFemtoKink* kinkCopy = new AliFemtoKink(**kIter);
      fKinkCollection->push_back(kinkCopy);
    }
  }
}
//___________________
AliFemtoEvent::AliFemtoEvent(const AliFemtoEvent& ev):
  fEventNumber(ev.fEventNumber),
  fRunNumber(ev.fRunNumber),
  fNumberOfTracks(ev.fNumberOfTracks),
  fNormalizedMult(ev.fNormalizedMult),
  fSPDMult(ev.fSPDMult),
  fEstimateITSTPC(ev.fEstimateITSTPC),
  fEstimateTracklets(ev.fEstimateTracklets),
  fEstimateITSPure(ev.fEstimateITSPure),
  fCentralityV0(ev.fCentralityV0),
  fCentralityV0A(ev.fCentralityV0A),
  fCentralityV0C(ev.fCentralityV0C),
  fCentralityZNA(ev.fCentralityZNA),
  fCentralityZNC(ev.fCentralityZNC),
  fCentralityCL1(ev.fCentralityCL1),
  fCentralityCL0(ev.fCentralityCL0),
  fCentralityTKL(ev.fCentralityTKL),
  fCentralityFMD(ev.fCentralityFMD),
  fCentralityTrk(ev.fCentralityTrk),
  fCentralityCND(ev.fCentralityCND),
  fCentralityNPA(ev.fCentralityNPA),
  fCentralitySPD1(ev.fCentralitySPD1),
  fMagneticField(ev.fMagneticField),
  fIsCollisionCandidate(ev.fIsCollisionCandidate),
  fPrimVertPos(ev.fPrimVertPos),
  fPrimVertCov(),
  fTrackCollection(NULL),
  fV0Collection(NULL),
  fXiCollection(NULL),
  fKinkCollection(NULL),
  fZDCN1Energy(ev.fZDCN1Energy),
  fZDCP1Energy(ev.fZDCP1Energy),
  fZDCN2Energy(ev.fZDCN2Energy),
  fZDCP2Energy(ev.fZDCP2Energy),
  fZDCEMEnergy(ev.fZDCEMEnergy),
  fZDCParticipants(ev.fZDCParticipants),
  fTriggerMask(ev.fTriggerMask),
  fTriggerCluster(ev.fTriggerCluster),
  fReactionPlaneAngle(ev.fReactionPlaneAngle),
  fEP(ev.fEP)
{
  // copy constructor
  SetPrimVertCov(ev.PrimVertCov());

  // create collections
  fTrackCollection = new AliFemtoTrackCollection;
  fV0Collection = new AliFemtoV0Collection;
  fXiCollection = new AliFemtoXiCollection;
  fKinkCollection = new AliFemtoKinkCollection;
  // copy track collection
  for ( AliFemtoTrackIterator tIter=ev.fTrackCollection->begin(); tIter!=ev.fTrackCollection->end(); tIter++) {
    AliFemtoTrack* trackCopy = new AliFemtoTrack(**tIter);
    fTrackCollection->push_back(trackCopy);
  }
  // copy v0 collection
  for ( AliFemtoV0Iterator vIter=ev.fV0Collection->begin(); vIter!=ev.fV0Collection->end(); vIter++) {
    AliFemtoV0* v0Copy = new AliFemtoV0(**vIter);
    fV0Collection->push_back(v0Copy);
  }
  // copy xi collection
  for ( AliFemtoXiIterator xIter=ev.fXiCollection->begin(); xIter!=ev.fXiCollection->end(); xIter++) {
    AliFemtoXi* xiCopy = new AliFemtoXi(**xIter);
    fXiCollection->push_back(xiCopy);
  }
  // copy kink collection
  for ( AliFemtoKinkIterator kIter=ev.fKinkCollection->begin(); kIter!=ev.fKinkCollection->end(); kIter++) {
    //cout << " kinkCut passed " << endl;
    AliFemtoKink* kinkCopy = new AliFemtoKink(**kIter);
    fKinkCollection->push_back(kinkCopy);
  }
}
//______________________________
AliFemtoEvent& AliFemtoEvent::operator=(const AliFemtoEvent& aEvent)
{
  // assignment operator
  if (this == &aEvent)
    return *this;

  fEventNumber = aEvent.fEventNumber;
  fRunNumber = aEvent.fRunNumber;

  fZDCN1Energy = aEvent.fZDCN1Energy;
  fZDCP1Energy = aEvent.fZDCP1Energy;
  fZDCN2Energy = aEvent.fZDCN2Energy;
  fZDCP2Energy = aEvent.fZDCP2Energy;
  fZDCEMEnergy = aEvent.fZDCEMEnergy;
  fZDCParticipants = aEvent.fZDCParticipants;
  fNumberOfTracks = aEvent.fNumberOfTracks;
  fEstimateITSTPC = aEvent.fEstimateITSTPC;
  fEstimateTracklets = aEvent.fEstimateTracklets;
  fEstimateITSPure = aEvent.fEstimateITSPure;
  fCentralityV0 = aEvent.fCentralityV0;
  fCentralityV0A = aEvent.fCentralityV0A;
  fCentralityV0C = aEvent.fCentralityV0C;
  fCentralityZNA = aEvent.fCentralityZNA;
  fCentralityZNC = aEvent.fCentralityZNC;
  fCentralityCL1 = aEvent.fCentralityCL1;
  fCentralityCL0 = aEvent.fCentralityCL0;
  fCentralityTKL = aEvent.fCentralityTKL;
  fCentralityFMD = aEvent.fCentralityFMD;
  fCentralityTrk = aEvent.fCentralityTrk;
  fCentralityCND = aEvent.fCentralityCND;
  fCentralityNPA = aEvent.fCentralityNPA;
  fCentralitySPD1 = aEvent.fCentralitySPD1;
  fNormalizedMult = aEvent.fNormalizedMult;
  fEstimateITSTPC = aEvent.fEstimateITSTPC;
  fEstimateTracklets = aEvent.fEstimateTracklets;
  fEstimateITSPure = aEvent.fEstimateITSPure;
  fMagneticField = aEvent.fMagneticField;
  fIsCollisionCandidate = aEvent.fIsCollisionCandidate;
  fPrimVertPos = aEvent.fPrimVertPos;

  fTriggerMask = aEvent.fTriggerMask;     // Trigger Type (mask)
  fTriggerCluster = aEvent.fTriggerCluster;
  fReactionPlaneAngle = aEvent.fReactionPlaneAngle;
  fEP = aEvent.fEP;

  if (fTrackCollection) {
    for (AliFemtoTrackIterator iter=fTrackCollection->begin();iter!=fTrackCollection->end();iter++){
      delete *iter;
    }
    fTrackCollection->clear();
  } else {
    fTrackCollection = new AliFemtoTrackCollection;
  }

  if (fV0Collection) {
    for (AliFemtoV0Iterator tV0iter=fV0Collection->begin();tV0iter!=fV0Collection->end();tV0iter++){
      delete *tV0iter;
    }//added by M Chojnacki To avodid memory leak
    fV0Collection->clear();
  } else {
    fV0Collection = new AliFemtoV0Collection;
  }

  if (fXiCollection) {
    for (AliFemtoXiIterator tXiIter=fXiCollection->begin();tXiIter!=fXiCollection->end();tXiIter++){
      delete *tXiIter;
    }
    fXiCollection->clear();
  } else {
    fXiCollection = new AliFemtoXiCollection;
  }

  if (fKinkCollection) {
    for (AliFemtoKinkIterator kinkIter=fKinkCollection->begin();kinkIter!=fKinkCollection->end();kinkIter++){
      delete *kinkIter;
    }
    fKinkCollection->clear();
  } else {
    fKinkCollection = new AliFemtoKinkCollection;
  }

  // copy track collection
  for ( AliFemtoTrackIterator tIter=aEvent.fTrackCollection->begin(); tIter!=aEvent.fTrackCollection->end(); tIter++) {
    AliFemtoTrack* trackCopy = new AliFemtoTrack(**tIter);
    fTrackCollection->push_back(trackCopy);
  }
  // copy v0 collection
  for ( AliFemtoV0Iterator vIter=aEvent.fV0Collection->begin(); vIter!=aEvent.fV0Collection->end(); vIter++) {
    AliFemtoV0* v0Copy = new AliFemtoV0(**vIter);
    fV0Collection->push_back(v0Copy);
  }
  // copy xi collection
  for ( AliFemtoXiIterator xIter=aEvent.fXiCollection->begin(); xIter!=aEvent.fXiCollection->end(); xIter++) {
    AliFemtoXi* xiCopy = new AliFemtoXi(**xIter);
    fXiCollection->push_back(xiCopy);
  }
  // copy kink collection
  for ( AliFemtoKinkIterator kIter=aEvent.fKinkCollection->begin(); kIter!=aEvent.fKinkCollection->end(); kIter++) {
    AliFemtoKink* kinkCopy = new AliFemtoKink(**kIter);
    fKinkCollection->push_back(kinkCopy);
  }

  return *this;
}

//___________________
AliFemtoEvent::~AliFemtoEvent(){
  // destructor
#ifdef STHBTDEBUG
  cout << " AliFemtoEvent::~AliFemtoEvent() " << endl;
#endif
  for (AliFemtoTrackIterator iter=fTrackCollection->begin();iter!=fTrackCollection->end();iter++){
    delete *iter;
  }
  fTrackCollection->clear();
  delete fTrackCollection;
  //must do the same for the V0 collection
  for (AliFemtoV0Iterator tV0iter=fV0Collection->begin();tV0iter!=fV0Collection->end();tV0iter++){
    delete *tV0iter;
  }//added by M Chojnacki To avodid memory leak
  fV0Collection->clear();
  delete fV0Collection;
  //must do the same for the Xi collection
  for (AliFemtoXiIterator tXiIter=fXiCollection->begin();tXiIter!=fXiCollection->end();tXiIter++){
    delete *tXiIter;
  }
  fXiCollection->clear();
  delete fXiCollection;
  //must do the same for the Kink collection
  for (AliFemtoKinkIterator kinkIter=fKinkCollection->begin();kinkIter!=fKinkCollection->end();kinkIter++){
    delete *kinkIter;
  }
  fKinkCollection->clear();
  delete fKinkCollection;
}
//___________________



void AliFemtoEvent::SetEventNumber(const unsigned short& event){fEventNumber = event;}
void AliFemtoEvent::SetRunNumber(const int& runNum){fRunNumber = runNum;}


void AliFemtoEvent::SetZDCN1Energy(const float& aZDCN1Energy){fZDCN1Energy=aZDCN1Energy;}
void AliFemtoEvent::SetZDCP1Energy(const float& aZDCP1Energy){fZDCP1Energy=aZDCP1Energy;}
void AliFemtoEvent::SetZDCN2Energy(const float& aZDCN2Energy){fZDCN2Energy=aZDCN2Energy;}
void AliFemtoEvent::SetZDCP2Energy(const float& aZDCP2Energy){fZDCP2Energy=aZDCP2Energy;}
void AliFemtoEvent::SetZDCEMEnergy(const float& aZDCEMEnergy){fZDCEMEnergy=aZDCEMEnergy;}
void AliFemtoEvent::SetZDCParticipants(const unsigned int& aZDCParticipants){fZDCParticipants=aZDCParticipants;}

void AliFemtoEvent::SetNumberOfTracks(const unsigned short& tracks){fNumberOfTracks = tracks;}
void AliFemtoEvent::SetNormalizedMult(const int& i){fNormalizedMult = i;}
void AliFemtoEvent::SetSPDMult(const int& i){fSPDMult = i;}

void AliFemtoEvent::SetPrimVertPos(const AliFemtoThreeVector& vp){fPrimVertPos = vp;}
void AliFemtoEvent::SetPrimVertCov(const double* v){
  fPrimVertCov[0] = v[0];
  fPrimVertCov[1] = v[1];
  fPrimVertCov[2] = v[2];
  fPrimVertCov[3] = v[3];
  fPrimVertCov[4] = v[4];
  fPrimVertCov[5] = v[5];
}
void AliFemtoEvent::SetMagneticField(const double& magF){fMagneticField = magF;}
void AliFemtoEvent::SetIsCollisionCandidate(const bool& is){fIsCollisionCandidate = is;}

void AliFemtoEvent::SetTriggerMask(const unsigned long int& aTriggerMask) {fTriggerMask=aTriggerMask;}
void AliFemtoEvent::SetTriggerCluster(const unsigned char& aTriggerCluster) {fTriggerCluster=aTriggerCluster;}


unsigned short AliFemtoEvent::EventNumber() const {return fEventNumber;}
int            AliFemtoEvent::RunNumber() const {return fRunNumber;}



unsigned short AliFemtoEvent::NumberOfTracks() const {return fNumberOfTracks;}

AliFemtoV0Collection* AliFemtoEvent::V0Collection() const {return fV0Collection;}
AliFemtoXiCollection* AliFemtoEvent::XiCollection() const {return fXiCollection;}
AliFemtoKinkCollection* AliFemtoEvent::KinkCollection() const {return fKinkCollection;}
AliFemtoTrackCollection* AliFemtoEvent::TrackCollection() const {return fTrackCollection;}
AliFemtoThreeVector AliFemtoEvent::PrimVertPos() const {return fPrimVertPos;}
const double* AliFemtoEvent::PrimVertCov() const {return fPrimVertCov;}
double AliFemtoEvent::MagneticField() const {return fMagneticField;}
unsigned long int AliFemtoEvent::TriggerMask() const {return fTriggerMask;}
unsigned char AliFemtoEvent::TriggerCluster() const {return fTriggerCluster;}
bool AliFemtoEvent::IsCollisionCandidate() const {return fIsCollisionCandidate;}


float AliFemtoEvent::ZDCN1Energy() const {return fZDCN1Energy;}
float AliFemtoEvent::ZDCP1Energy() const {return fZDCP1Energy;}
float AliFemtoEvent::ZDCN2Energy() const {return fZDCN2Energy;}
float AliFemtoEvent::ZDCP2Energy() const {return fZDCP2Energy;}
float AliFemtoEvent::ZDCEMEnergy() const {return fZDCEMEnergy;}
unsigned int  AliFemtoEvent::ZDCParticipants() const {return fZDCParticipants;}

void AliFemtoEvent::SetReactionPlaneAngle(const float& a) { fReactionPlaneAngle = a;}
float AliFemtoEvent::ReactionPlaneAngle() const { return fReactionPlaneAngle; }
void AliFemtoEvent::SetEP(AliEventplane* ep) { fEP = ep;}
AliEventplane* AliFemtoEvent::EP() const {return fEP; }
 //fV0perEvent->Sumw2();
//----------------------------- below here is only for star

int AliFemtoEvent::UncorrectedNumberOfNegativePrimaries() const
{
  return NumberOfTracks()/2;
}

int AliFemtoEvent::SPDMultiplicity() const
{
  return fSPDMult;
}

int AliFemtoEvent::NumberOfV0s() const
{
  return V0Collection()->size();
}

int AliFemtoEvent::UncorrectedNumberOfPrimaries() const
{
  if (fNormalizedMult < -1) {
    // Count number of normalized charged tracks
    Int_t tNormTrackCount = 0;
    for (AliFemtoTrackIterator iter=fTrackCollection->begin();iter!=fTrackCollection->end();iter++){
      if (!((*iter)->Flags()&(AliFemtoTrack::kTPCrefit))) continue;
      if ((*iter)->TPCncls() < 50) continue;
      if ((*iter)->TPCchi2()/(*iter)->TPCncls() > 60.0) continue;
      if ((*iter)->ImpactD() > 6.0) continue;
      if ((*iter)->ImpactZ() > 6.0) continue;
      if (fabs((*iter)->P().PseudoRapidity()) > 0.9) continue;

      tNormTrackCount++;
    }
    return tNormTrackCount;
  }

  return fNormalizedMult;
  //  return NumberOfTracks();
}

unsigned short AliFemtoEvent::MultiplicityEstimateITSTPC() const
{
  return fEstimateITSTPC;
}

unsigned short AliFemtoEvent::MultiplicityEstimateTracklets() const
{
  return fEstimateTracklets;
}

unsigned short AliFemtoEvent::MultiplicityEstimateITSPure() const
{
  return fEstimateITSPure;
}

void AliFemtoEvent::SetMultiplicityEstimateITSTPC(const unsigned short &s)
{
  fEstimateITSTPC = s;
}

void AliFemtoEvent::SetMultiplicityEstimateTracklets(const unsigned short &s)
{
  fEstimateTracklets = s;
}

void AliFemtoEvent::SetMultiplicityEstimateITSPure(const unsigned short &s)
{
  fEstimateITSPure = s;
}

void AliFemtoEvent::SetCentralityV0(const float &c)
{
  fCentralityV0 = c;
}
void AliFemtoEvent::SetCentralityV0A(const float &c)
{
  fCentralityV0A = c;
}
void AliFemtoEvent::SetCentralityV0C(const float &c)
{
  fCentralityV0C = c;
}
void AliFemtoEvent::SetCentralityZNA(const float &c)
{
  fCentralityZNA = c;
}
void AliFemtoEvent::SetCentralityZNC(const float &c)
{
  fCentralityZNC = c;
}
void AliFemtoEvent::SetCentralityCL1(const float &c)
{
  fCentralityCL1 = c;
}
void AliFemtoEvent::SetCentralityCL0(const float &c)
{
  fCentralityCL0 = c;
}
void AliFemtoEvent::SetCentralityTKL(const float &c)
{
  fCentralityTKL = c;
}
void AliFemtoEvent::SetCentralityFMD(const float &c)
{
  fCentralityFMD = c;
}
void AliFemtoEvent::SetCentralityTrk(const float &c)
{
  fCentralityTrk = c;
}
void AliFemtoEvent::SetCentralityCND(const float &c)
{
  fCentralityCND = c;
}
void AliFemtoEvent::SetCentralityNPA(const float &c)
{
  fCentralityNPA = c;
}
void AliFemtoEvent::SetCentralitySPD1(const float &c)
{
  fCentralitySPD1 = c;
}

float AliFemtoEvent::CentralityV0() const
{
  return fCentralityV0;
}
float AliFemtoEvent::CentralityV0A() const
{
  return fCentralityV0A;
}
float AliFemtoEvent::CentralityV0C() const
{
  return fCentralityV0C;
}
float AliFemtoEvent::CentralityZNA() const
{
  return fCentralityZNA;
}
float AliFemtoEvent::CentralityZNC() const
{
  return fCentralityZNC;
}
float AliFemtoEvent::CentralityCL1() const
{
  return fCentralityCL1;
}
float AliFemtoEvent::CentralityCL0() const
{
  return fCentralityCL0;
}
float AliFemtoEvent::CentralityTKL() const
{
  return fCentralityTKL;
}
float AliFemtoEvent::CentralityFMD() const
{
  return fCentralityFMD;
}
float AliFemtoEvent::CentralityTrk() const
{
  return fCentralityTrk;
}
float AliFemtoEvent::CentralityCND() const
{
  return fCentralityCND;
}
float AliFemtoEvent::CentralityNPA() const
{
  return fCentralityNPA;
}
float AliFemtoEvent::CentralitySPD1() const
{
  return fCentralitySPD1;
}
