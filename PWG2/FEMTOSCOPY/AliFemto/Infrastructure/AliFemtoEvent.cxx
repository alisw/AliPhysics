/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *   HbtEvent is the "transient microDST"  Objects of this class are
 *   generated from the input data by a Reader, and then presented to
 *   the Cuts of the various active Analyses.
 *
 ***************************************************************************
 *
 * $Log$
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.2  2007-04-03 16:00:08  mchojnacki
 * Changes to iprove memory managing
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.23  2005/08/19 21:19:11  chajecki
 * line 326: the same change as in line 235
 *
 * Revision 1.22  2005/08/19 11:33:31  chajecki
 * fix due to the last changes in MuDst
 * line 235: TClonesArray* tracks=0; to TObjArray* tracks=0;
 * see for more details:
 * http://www.star.bnl.gov/HyperNews-star/protected/get/starsoft/5949.html
 *
 * Revision 1.21  2003/09/02 17:58:32  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.20  2003/01/31 19:43:20  magestro
 * several casts added to remove compiler warnings
 *
 * Revision 1.19  2003/01/17 16:46:58  mercedes
 * StMuEvent::refMult() added
 *
 * Revision 1.18  2002/11/19 23:27:37  renault
 * New event constructor to find V0 daughters informations(helix for average
 * separation calculation)
 *
 * Revision 1.17  2002/03/21 18:49:31  laue
 * updated for new MuDst reader
 *
 * Revision 1.16  2001/12/06 16:47:13  laue
 * l3 trigger algorithm added
 *
 * Revision 1.15  2001/11/14 21:07:21  lisa
 * Fixed several small things (mostly discarded const) that caused fatal errors with gcc2.95.3
 *
 * Revision 1.14  2001/09/05 20:41:42  laue
 * Updates of the hbtMuDstTree microDSTs
 *
 * Revision 1.13  2001/07/20 20:03:53  rcwells
 * Added pT weighting and moved event angle cal. to event cut
 *
 * Revision 1.12  2001/06/23 21:55:17  laue
 * AliFemtoCheckPdgIdList can take can not check for mother,particle,daughter
 * Some output turned off
 *
 * Revision 1.11  2001/06/21 19:15:45  laue
 * Modified fiels:
 *   CTH.h : new constructor added
 *   AliFemtoEvent, AliFemtoKink, AliFemtoTrack : constructors from the persistent
 *                                   (TTree) classes added
 *   AliFemtoLikeSignAnalysis : minor changes, for debugging
 *   AliFemtoTypes: split into different files
 * Added files: for the new TTree muDst's
 *   StExceptions.cxx StExceptions.h AliFemtoEnumeration.h
 *   AliFemtoHelix.h AliFemtoHisto.h AliFemtoString.h AliFemtoTFile.h
 *   AliFemtoTTreeEvent.cxx AliFemtoTTreeEvent.h AliFemtoTTreeKink.cxx
 *   AliFemtoTTreeKink.h AliFemtoTTreeTrack.cxx AliFemtoTTreeTrack.h
 *   AliFemtoTTreeV0.cxx AliFemtoTTreeV0.h AliFemtoVector.h
 *
 * Revision 1.10  2001/05/15 15:30:16  rcwells
 * Added magnetic field to AliFemtoEvent
 *
 * Revision 1.9  2000/08/31 22:31:31  laue
 * AliFemtoAnalysis: output changed (a little bit less)
 * AliFemtoEvent: new version, members for reference mult added
 * AliFemtoIOBinary: new IO for new AliFemtoEvent version
 * AliFemtoTypes: TTree typedef to AliFemtoTTree added
 * AliFemtoVertexAnalysis: overflow and underflow added
 *
 * Revision 1.8  2000/07/16 21:38:22  laue
 * AliFemtoCoulomb.cxx AliFemtoSectoredAnalysis.cxx : updated for standalone version
 * AliFemtoV0.cc AliFemtoV0.h : some cast to prevent compiling warnings
 * AliFemtoParticle.cc AliFemtoParticle.h : pointers mTrack,mV0 initialized to 0
 * AliFemtoIOBinary.cc : some printouts in #ifdef STHBTDEBUG
 * AliFemtoEvent.cc : B-Field set to 0.25Tesla, we have to think about a better
 *                 solution
 *
 * Revision 1.7  2000/05/25 21:54:16  laue
 * RotateZ implemented. Rotates momentum and helix around the z axis
 *
 * Revision 1.5  2000/02/18 21:32:23  laue
 * franksTrackCut changed. If mCharge is set to '0' there will be no cut
 * on charge. This is important for front-loaded cuts.
 *
 * copy constructor implemented for AliFemtoEvent, AliFemtoTrack and AliFemtoV0.
 *
 * franks1HistoD.cxx franks1HistoD.h franks2HistoD.cxx franks2HistoD.h
 * removed. We can now (CC5 on Solaris) use the versions (no D)
 *
 * Revision 1.4  1999/09/16 18:47:59  lisa
 * replace placeholder HbtV0Track stuff with Helens AliFemtoV0 classes
 *
 * Revision 1.3  1999/07/27 10:47:04  lisa
 * now works in dev on linux and solaris - mistake in deleting picoEvents fixed
 *
 * Revision 1.2  1999/07/19 14:24:05  hardtke
 * modifications to implement uDST
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#include "Infrastructure/AliFemtoEvent.h"
#include "Infrastructure/AliFemtoTrack.h"
#include "Infrastructure/AliFemtoV0.h"
#include "Infrastructure/AliFemtoXi.h"
#include "Infrastructure/AliFemtoKink.h"
#include "Base/AliFemtoTrackCut.h"
#include "Base/AliFemtoV0Cut.h"
#include "Base/AliFemtoXiCut.h"
#include "Base/AliFemtoKinkCut.h"
#include "Base/PhysicalConstants.h"
#include "Base/SystemOfUnits.h"

// Mike removed all of the AliFemtoTTree stuff here 21apr2006 - it was not used for a long time.



//___________________
AliFemtoEvent::AliFemtoEvent(){
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
AliFemtoEvent::AliFemtoEvent(const AliFemtoEvent& ev, AliFemtoTrackCut* tCut, AliFemtoV0Cut* vCut, AliFemtoXiCut* xCut, AliFemtoKinkCut* kCut){ // copy constructor with track and v0 cuts
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
AliFemtoEvent::~AliFemtoEvent(){
#ifdef STHBTDEBUG
  cout << " AliFemtoEvent::~AliFemtoEvent() " << endl;
#endif
  for (AliFemtoTrackIterator iter=fTrackCollection->begin();iter!=fTrackCollection->end();iter++){
    delete *iter;
  }
  fTrackCollection->clear();
  delete fTrackCollection;
  //must do the same for the V0 collection
  for (AliFemtoV0Iterator V0iter=fV0Collection->begin();V0iter!=fV0Collection->end();V0iter++){
    delete *V0iter;
  }//added by M Chojnacki To avodid memory leak 
  fV0Collection->clear();
  delete fV0Collection;
  //must do the same for the Xi collection
  for (AliFemtoXiIterator XiIter=fXiCollection->begin();XiIter!=fXiCollection->end();XiIter++){
    delete *XiIter;
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


