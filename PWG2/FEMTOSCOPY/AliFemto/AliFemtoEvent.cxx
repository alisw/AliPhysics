///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  AliFemtoEvent: hold the information specific to the event and a      //
//  track list                                                           //
//  AliFemtoEvent is the "transient microDST"  Objects of this class are //
//   generated from the input data by a Reader, and then presented to	 //
//   the Cuts of the various active Analyses.				 //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
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

// Mike removed all of the AliFemtoTTree stuff here 21apr2006 - it was not used for a long time.



//___________________
AliFemtoEvent::AliFemtoEvent():
  fEventNumber(0),
  fRunNumber(0),
  fNumberOfTracks(0),
  fMagneticField(0),
  fPrimVertPos(0,0,0),
  fTrackCollection(0),
  fV0Collection(0),
  fXiCollection(0),
  fKinkCollection(0),
  fZDCN1Energy(0),   
  fZDCP1Energy(0),   
  fZDCN2Energy(0),   
  fZDCP2Energy(0),   
  fZDCEMEnergy(0),   
  fZDCParticipants(0),
  fTriggerMask(0),  
  fTriggerCluster(0)
{
  // Default constructor
  fPrimVertPos[0]=-999.0;
  fPrimVertPos[1]=-999.0;
  fPrimVertPos[2]=-999.0;
  fTrackCollection = new AliFemtoTrackCollection;
  fV0Collection = new AliFemtoV0Collection;
  fXiCollection = new AliFemtoXiCollection;
  fKinkCollection = new AliFemtoKinkCollection;
  fMagneticField=0.0;
}
//___________________
AliFemtoEvent::AliFemtoEvent(const AliFemtoEvent& ev, AliFemtoTrackCut* tCut, AliFemtoV0Cut* vCut, AliFemtoXiCut* xCut, AliFemtoKinkCut* kCut):
  fEventNumber(0),
  fRunNumber(0),
  fNumberOfTracks(0),
  fMagneticField(0),
  fPrimVertPos(0,0,0),
  fTrackCollection(0),
  fV0Collection(0),
  fXiCollection(0),
  fKinkCollection(0),
  fZDCN1Energy(0),   
  fZDCP1Energy(0),   
  fZDCN2Energy(0),   
  fZDCP2Energy(0),   
  fZDCEMEnergy(0),   
  fZDCParticipants(0),
  fTriggerMask(0),  
  fTriggerCluster(0)
{ // copy constructor with track and v0 cuts
  //cout << "AliFemtoEvent::AliFemtoEvent(const AliFemtoEvent& ev, AliFemtoTrackCut* tCut, AliFemtoV0Cut* vCut, AliFemtoV0Cut* kCut)" << endl;
  fEventNumber = ev.fEventNumber;
  fRunNumber = ev.fRunNumber;
  
  fZDCN1Energy=ev.fZDCN1Energy;     
  fZDCP1Energy=ev.fZDCP1Energy;      
  fZDCN2Energy=ev.fZDCN2Energy;      
  fZDCP2Energy=ev.fZDCP2Energy;      
  fZDCEMEnergy=ev.fZDCEMEnergy;
  fZDCParticipants=ev.fZDCParticipants;
  fNumberOfTracks = ev.fNumberOfTracks;
  fMagneticField= ev.fMagneticField;
  
  fTriggerMask=ev.fTriggerMask;     // Trigger Type (mask)
  fTriggerCluster=ev.fTriggerCluster;
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
  fEventNumber(0),
  fRunNumber(0),
  fNumberOfTracks(0),
  fMagneticField(0),
  fPrimVertPos(0,0,0),
  fTrackCollection(0),
  fV0Collection(0),
  fXiCollection(0),
  fKinkCollection(0),
  fZDCN1Energy(0),   
  fZDCP1Energy(0),   
  fZDCN2Energy(0),   
  fZDCP2Energy(0),   
  fZDCEMEnergy(0),   
  fZDCParticipants(0),
  fTriggerMask(0),  
  fTriggerCluster(0)
{ 
  // copy constructor 
  fEventNumber = ev.fEventNumber;
  fRunNumber = ev.fRunNumber;
  
  fZDCN1Energy=ev.fZDCN1Energy;     
  fZDCP1Energy=ev.fZDCP1Energy;      
  fZDCN2Energy=ev.fZDCN2Energy;      
  fZDCP2Energy=ev.fZDCP2Energy;      
  fZDCEMEnergy=ev.fZDCEMEnergy;
  fZDCParticipants=ev.fZDCParticipants;
  fNumberOfTracks = ev.fNumberOfTracks;
  fMagneticField= ev.fMagneticField;
  
  fTriggerMask=ev.fTriggerMask;     // Trigger Type (mask)
  fTriggerCluster=ev.fTriggerCluster;
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
  
  fZDCN1Energy=aEvent.fZDCN1Energy;     
  fZDCP1Energy=aEvent.fZDCP1Energy;      
  fZDCN2Energy=aEvent.fZDCN2Energy;      
  fZDCP2Energy=aEvent.fZDCP2Energy;      
  fZDCEMEnergy=aEvent.fZDCEMEnergy;
  fZDCParticipants=aEvent.fZDCParticipants;
  fNumberOfTracks = aEvent.fNumberOfTracks;
  fMagneticField= aEvent.fMagneticField;
  
  fTriggerMask=aEvent.fTriggerMask;     // Trigger Type (mask)
  fTriggerCluster=aEvent.fTriggerCluster;
  // create collections
  fTrackCollection = new AliFemtoTrackCollection;
  fV0Collection = new AliFemtoV0Collection;
  fXiCollection = new AliFemtoXiCollection;
  fKinkCollection = new AliFemtoKinkCollection;
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


void AliFemtoEvent::SetZDCN1Energy(const float& ZDCN1Energy){fZDCN1Energy=ZDCN1Energy;}
void AliFemtoEvent::SetZDCP1Energy(const float& ZDCP1Energy){fZDCP1Energy=ZDCP1Energy;}      
void AliFemtoEvent::SetZDCN2Energy(const float& ZDCN2Energy){fZDCN2Energy=ZDCN2Energy;}      
void AliFemtoEvent::SetZDCP2Energy(const float& ZDCP2Energy){fZDCP2Energy=ZDCP2Energy;}      
void AliFemtoEvent::SetZDCEMEnergy(const float& ZDCEMEnergy){fZDCEMEnergy=ZDCEMEnergy;}    
void AliFemtoEvent::SetZDCParticipants(const unsigned int& ZDCParticipants){fZDCParticipants=ZDCParticipants;}
    

void AliFemtoEvent::SetNumberOfTracks(const unsigned short& tracks){fNumberOfTracks = tracks;}



void AliFemtoEvent::SetPrimVertPos(const AliFemtoThreeVector& vp){fPrimVertPos = vp;}
void AliFemtoEvent::SetMagneticField(const double& magF){fMagneticField = magF;}

void AliFemtoEvent::SetTriggerMask(const unsigned long int& TriggerMask) {fTriggerMask=TriggerMask;}
void AliFemtoEvent::SetTriggerCluster(const unsigned char& TriggerCluster) {fTriggerCluster=TriggerCluster;}


unsigned short AliFemtoEvent::EventNumber() const {return fEventNumber;}
int            AliFemtoEvent::RunNumber() const {return fRunNumber;}



unsigned short AliFemtoEvent::NumberOfTracks() const {return fNumberOfTracks;}

AliFemtoV0Collection* AliFemtoEvent::V0Collection() const {return fV0Collection;}
AliFemtoXiCollection* AliFemtoEvent::XiCollection() const {return fXiCollection;}
AliFemtoKinkCollection* AliFemtoEvent::KinkCollection() const {return fKinkCollection;}
AliFemtoTrackCollection* AliFemtoEvent::TrackCollection() const {return fTrackCollection;}
AliFemtoThreeVector AliFemtoEvent::PrimVertPos() const {return fPrimVertPos;}
double AliFemtoEvent::MagneticField() const {return fMagneticField;}
unsigned long int AliFemtoEvent::TriggerMask() const {return fTriggerMask;}
unsigned char AliFemtoEvent::TriggerCluster() const {return fTriggerCluster;}


float AliFemtoEvent::ZDCN1Energy() const {return fZDCN1Energy;}       
float AliFemtoEvent::ZDCP1Energy() const {return fZDCP1Energy;}       
float AliFemtoEvent::ZDCN2Energy() const {return fZDCN2Energy;}       
float AliFemtoEvent::ZDCP2Energy() const {return fZDCP2Energy;}       
float AliFemtoEvent::ZDCEMEnergy() const {return fZDCEMEnergy;}   
unsigned int  AliFemtoEvent::ZDCParticipants() const {return fZDCParticipants;}

//----------------------------- below here is only for star

double AliFemtoEvent::UncorrectedNumberOfNegativePrimaries() const
{
  return NumberOfTracks()/2;
}

double AliFemtoEvent::UncorrectedNumberOfPrimaries() const
{
  return NumberOfTracks();
}


