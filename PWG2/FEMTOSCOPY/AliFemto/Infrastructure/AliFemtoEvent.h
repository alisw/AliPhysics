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
 * Revision 1.21 to use in Alice version 1 Chojnacki  
 *
 *
 * $Log$
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.20  2003/01/17 16:46:22  mercedes
 * StMuEvent::refMult() added
 *
 * Revision 1.19  2002/11/19 23:27:25  renault
 * New event constructor to find V0 daughters informations(helix for average
 * separation calculation)
 *
 * Revision 1.18  2002/03/21 18:49:31  laue
 * updated for new MuDst reader
 *
 * Revision 1.17  2001/12/06 16:47:13  laue
 * l3 trigger algorithm added
 *
 * Revision 1.14  2001/06/21 19:15:45  laue
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
 * Revision 1.13  2001/06/04 19:09:52  rcwells
 * Adding B-field, run number, and improved reaction plane functionality
 *
 * Revision 1.12  2001/05/25 23:23:59  lisa
 * Added in AliFemtoKink stuff
 *
 * Revision 1.11  2001/05/15 15:30:16  rcwells
 * Added magnetic field to AliFemtoEvent
 *
 * Revision 1.10  2000/08/31 22:31:31  laue
 * AliFemtoAnalysis: output changed (a little bit less)
 * AliFemtoEvent: new version, members for reference mult added
 * AliFemtoIOBinary: new IO for new AliFemtoEvent version
 * AliFemtoTypes: TTree typedef to AliFemtoTTree added
 * AliFemtoVertexAnalysis: overflow and underflow added
 *
 * Revision 1.9  2000/05/25 21:54:16  laue
 * RotateZ implemented. Rotates momentum and helix around the z axis
 *
 * Revision 1.7  2000/02/18 21:32:23  laue
 * franksTrackCut changed. If mCharge is set to '0' there will be no cut
 * on charge. This is important for front-loaded cuts.
 *
 * copy constructor implemented for AliFemtoEvent, AliFemtoTrack and AliFemtoV0.
 *
 * franks1HistoD.cxx franks1HistoD.h franks2HistoD.cxx franks2HistoD.h
 * removed. We can now (CC5 on Solaris) use the versions (no D)
 *
 * Revision 1.6  1999/09/16 18:47:59  lisa
 * replace placeholder HbtV0Track stuff with Helens AliFemtoV0 classes
 *
 * Revision 1.5  1999/09/03 22:39:15  lisa
 * Readers now MUST have Report() methods and MAY have WriteHbtEvent() methods
 *
 * Revision 1.4  1999/07/19 14:24:06  hardtke
 * modifications to implement uDST
 *
 * Revision 1.3  1999/07/06 22:33:22  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.2  1999/06/29 17:50:27  fisyak
 * formal changes to account new StEvent, does not complie yet
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#ifndef AliFemtoEvent_hh
#define AliFemtoEvent_hh

#include "Infrastructure/AliFemtoTypes.h"
#include "Infrastructure/AliFemtoTrackCollection.h"
#include "Infrastructure/AliFemtoV0Collection.h"
#include "Infrastructure/AliFemtoXiCollection.h"
#include "Infrastructure/AliFemtoKinkCollection.h"

class AliFemtoTrackCut;
class AliFemtoV0Cut;
class AliFemtoXiCut;
class AliFemtoKinkCut;


#ifdef __ROOT__

// the following encapsulation by malisa 21apr2006
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
class StMuDst;
#endif

#endif

class AliFemtoEvent{
public:
  AliFemtoEvent();
#ifdef __ROOT__
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
//
#endif
#endif
  AliFemtoEvent(const AliFemtoEvent&, AliFemtoTrackCut* =0, AliFemtoV0Cut* =0,  AliFemtoXiCut* =0, AliFemtoKinkCut* =0); // copy constructor with track and v0 cuts
  ~AliFemtoEvent();

  unsigned short EventNumber() const;
  int RunNumber() const;
  unsigned short NumberOfTracks() const;
  AliFemtoThreeVector PrimVertPos() const;
  AliFemtoV0Collection* V0Collection() const;
  AliFemtoXiCollection* XiCollection() const;
  AliFemtoKinkCollection* KinkCollection() const;
  AliFemtoTrackCollection* TrackCollection() const;
  double MagneticField() const;

  //functions for alice variables
  float ZDCN1Energy() const;      
  float ZDCP1Energy() const;      
  float ZDCN2Energy() const;      
  float ZDCP2Energy() const;      
  float ZDCEMEnergy() const;    
  unsigned int ZDCParticipants() const; 
  
  unsigned long int     TriggerMask() const;     
  unsigned char      TriggerCluster() const;  
  
  void SetEventNumber(const unsigned short&);
  void SetRunNumber(const int&);
  void SetNumberOfTracks(const unsigned short&);
  void SetPrimVertPos(const AliFemtoThreeVector&);
  void SetMagneticField(const double&);
  
   //functions for alice variables
  void SetZDCN1Energy(const float&);      
  void SetZDCP1Energy(const float&);      
  void SetZDCN2Energy(const float&);      
  void SetZDCP2Energy(const float&);      
  void SetZDCEMEnergy(const float&);    
  void SetZDCParticipants(const unsigned int&);
  
  void SetTriggerMask(const unsigned long int&);     
  void SetTriggerCluster(const unsigned char&); 
  
  double UncorrectedNumberOfNegativePrimaries() const;
  double UncorrectedNumberOfPrimaries() const;

private:
  unsigned short fEventNumber;           //
  unsigned short fRunNumber;
  unsigned short fNumberOfTracks;     // total number of TPC tracks
  double fMagneticField; // magnetic field in Z direction

  AliFemtoThreeVector fPrimVertPos;
  AliFemtoTrackCollection* fTrackCollection;
  AliFemtoV0Collection* fV0Collection;
  AliFemtoXiCollection* fXiCollection;
  AliFemtoKinkCollection* fKinkCollection;

  //for alice changed by Marek Chojnacki
  float      fZDCN1Energy;      // reconstructed energy in the neutron ZDC
  float      fZDCP1Energy;      // reconstructed energy in the proton ZDC
  float      fZDCN2Energy;      // reconstructed energy in the neutron ZDC
  float      fZDCP2Energy;      // reconstructed energy in the proton ZDC
  float      fZDCEMEnergy;     // reconstructed energy in the electromagnetic ZDC
  unsigned int        fZDCParticipants; // number of participants estimated by the ZDC
  
  unsigned long int     fTriggerMask;     // Trigger Type (mask)
  unsigned char      fTriggerCluster;  // Trigger cluster (mask)
};



#endif 
