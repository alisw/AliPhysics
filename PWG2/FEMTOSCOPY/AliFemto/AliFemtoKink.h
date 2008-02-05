///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoKink: main class holding all the necessary information        //
// about a kink (before the identification) that is required during      //
// femtoscopic analysis. This class is filled with information from the  //
// input stream by the reader. A particle has a link back to the Kink    //
// it was created from, so we do not copy the information.               //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
/***********************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, 23May2001
 *
 ***********************************************************************
 *
 * Description: Kink class with information gotten from the StKinkVertex
 *              of Wenshen Deng and Spiros Margetis
 *
 ***********************************************************************
 *
 * $Log$
 * Revision 1.2  2007/05/03 09:42:29  akisiel
 * Fixing Effective C++ warnings
 *
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.4  2003/09/02 17:58:32  perev
 * gcc 3.2 updates + WarnOff
 *
 * Revision 1.3  2001/11/14 21:07:21  lisa
 * Fixed several small things (mostly discarded const) that caused fatal errors with gcc2.95.3
 *
 * Revision 1.2  2001/06/21 19:15:46  laue
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
 * Revision 1.1  2001/05/25 23:23:59  lisa
 * Added in AliFemtoKink stuff
 *
 * 
 *
 ***********************************************************************/
#ifndef ALIFEMTOKINK_H
#define ALIFEMTOKINK_H

class StKinkVertex;
//#include "StEvent/StKinkVertex.h"  // from StEvent
#include "AliFemtoTrack.h"

#include "AliFemtoTypes.h" //same as in AliFemtoTrack.h

class AliFemtoKink {
public:
  AliFemtoKink();
  AliFemtoKink( const AliFemtoKink& k); // copy constructor
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
#ifdef __ROOT__
  AliFemtoKink( const StKinkVertex&, AliFemtoThreeVector PrimaryVertex); // create a AliFemtoKink from a StKinkVertex
#endif
#endif
  ~AliFemtoKink(){/* no-op */}

  // Get's
  float        DcaParentDaughter() const;
  float        DcaDaughterPrimaryVertex() const;
  float        DcaParentPrimaryVertex() const;
  float        HitDistanceParentDaughter() const;
  float        HitDistanceParentVertex() const;
  float        DeltaEnergy(int i=0) const;
  float        DecayAngle() const;
  float        DecayAngleCM() const;
  AliFemtoTrack   Daughter() const;
  AliFemtoTrack   Parent() const;
  AliFemtoThreeVector Position() const; 

  

protected:

  float        fDcaParentDaughter;           // from StKinkVertex class directly 
  float        fDcaDaughterPrimaryVertex;    // from StKinkVertex class directly 
  float        fDcaParentPrimaryVertex;      // from StKinkVertex class directly 
  float        fHitDistanceParentDaughter;   // from StKinkVertex class directly 
  float        fHitDistanceParentVertex;     // from StKinkVertex class directly 
  float        fDeltaEnergy[3];              // from StKinkVertex class directly 
  float        fDecayAngle;                  // from StKinkVertex class directly 
  float        fDecayAngleCM;                // from StKinkVertex class directly 
  AliFemtoTrack   fDaughter;                    // from StKinkVertex class directly 
  AliFemtoTrack   fParent;                      // from StVertex class (which StKinkVertex inherits from)
  AliFemtoThreeVector fPosition;                // from StMeasuredPoint class (which StVertex inherits from)

};

// Get's
inline float        AliFemtoKink::DcaParentDaughter() const {return fDcaParentDaughter;}
inline float        AliFemtoKink::DcaDaughterPrimaryVertex() const {return fDcaDaughterPrimaryVertex;}
inline float        AliFemtoKink::DcaParentPrimaryVertex() const {return fDcaParentPrimaryVertex;}
inline float        AliFemtoKink::HitDistanceParentDaughter() const {return fHitDistanceParentDaughter;}
inline float        AliFemtoKink::HitDistanceParentVertex() const {return fHitDistanceParentVertex;}
inline float        AliFemtoKink::DeltaEnergy(int i) const {return fDeltaEnergy[i];}
inline float        AliFemtoKink::DecayAngle() const {return fDecayAngle;}
inline float        AliFemtoKink::DecayAngleCM() const {return fDecayAngleCM;}
inline AliFemtoTrack   AliFemtoKink::Daughter() const {return fDaughter;}
inline AliFemtoTrack   AliFemtoKink::Parent() const {return fParent;}
inline AliFemtoThreeVector AliFemtoKink::Position() const {return fPosition;}




#endif


















