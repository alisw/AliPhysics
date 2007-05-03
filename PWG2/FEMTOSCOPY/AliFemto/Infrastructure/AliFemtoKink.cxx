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
 * Revision 1.3  2007/04/27 07:24:34  akisiel
 * Make revisions needed for compilation from the main AliRoot tree
 *
 * Revision 1.1.1.1  2007/04/25 15:38:41  panos
 * Importing the HBT code dir
 *
 * Revision 1.1.1.1  2007/03/07 10:14:49  mchojnacki
 * First version on CVS
 *
 * Revision 1.4  2001/11/14 21:07:21  lisa
 * Fixed several small things (mostly discarded const) that caused fatal errors with gcc2.95.3
 *
 * Revision 1.3  2001/09/05 21:55:23  laue
 * typo fixed
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

#include "AliFemtoKink.h"
#include "Infrastructure/phys_constants.h"
#include "Infrastructure/AliFemtoTrack.h"
// -----------------------------------------------------------------------
AliFemtoKink::AliFemtoKink():
  fDcaParentDaughter(0), fDcaDaughterPrimaryVertex(0), 
  fDcaParentPrimaryVertex(0), fHitDistanceParentDaughter(0),  
  fHitDistanceParentVertex(0),
  fDecayAngle(0), fDecayAngleCM(0),             
  fDaughter(),            
  mParent(),               
  mPosition(0,0,0)        
{/* no-op */}
// -----------------------------------------------------------------------
AliFemtoKink::AliFemtoKink(const AliFemtoKink& k):
  fDcaParentDaughter(0), fDcaDaughterPrimaryVertex(0), 
  fDcaParentPrimaryVertex(0), fHitDistanceParentDaughter(0),  
  fHitDistanceParentVertex(0),
  fDecayAngle(0), fDecayAngleCM(0),             
  fDaughter(),            
  mParent(),               
  mPosition(0,0,0)        
{ // copy constructor

  fDcaParentDaughter          =   k.fDcaParentDaughter;           
  fDcaDaughterPrimaryVertex   =   k.fDcaDaughterPrimaryVertex;    
  fDcaParentPrimaryVertex     =   k.fDcaParentPrimaryVertex;      
  fHitDistanceParentDaughter  =   k.fHitDistanceParentDaughter;   
  fHitDistanceParentVertex    =   k.fHitDistanceParentVertex;     
  fDeltaEnergy[0]             =   k.fDeltaEnergy[0];              
  fDeltaEnergy[1]             =   k.fDeltaEnergy[1];              
  fDeltaEnergy[2]             =   k.fDeltaEnergy[2];              
  fDecayAngle                 =   k.fDecayAngle;                  
  fDecayAngleCM               =   k.fDecayAngleCM;                
  fDaughter                   =   k.fDaughter;                    
  mParent                     =   k.mParent;                      
  mPosition                   =   k.mPosition;                

}
// -----------------------------------------------------------------------


//--------------------- below here is ONLY star ----------------
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
#ifdef __ROOT__
#include "StEvent/StTrack.h"
#include "StEvent/StKinkVertex.h"
AliFemtoKink::AliFemtoKink( const StKinkVertex& SKV, AliFemtoThreeVector PrimaryVertex )
{ 

  fDcaParentDaughter          = SKV.dcaParentDaughter();
  fDcaDaughterPrimaryVertex   = SKV.dcaDaughterPrimaryVertex();
  fDcaParentPrimaryVertex     = SKV.dcaParentPrimaryVertex();
  fHitDistanceParentDaughter  = SKV.hitDistanceParentDaughter();
  fHitDistanceParentVertex    = SKV.hitDistanceParentVertex();
  fDeltaEnergy[0]             = SKV.dE(0);
  fDeltaEnergy[1]             = SKV.dE(1);
  fDeltaEnergy[2]             = SKV.dE(2);
  fDecayAngle                 = SKV.decayAngle();
  fDecayAngleCM               = SKV.decayAngleCM();

  // now fill member AliFemtoTrack data...
  const StTrack* StTrk;
  AliFemtoTrack* HbtTrk;
  // Daughter
  StTrk = SKV.daughter(0);
  HbtTrk = new AliFemtoTrack(StTrk,PrimaryVertex); // generate NEW HbtTrack from StTrack
  fDaughter = *HbtTrk;                         // invoke copy ctr of AliFemtoTrack
  delete HbtTrk;                               // get rid of the NEW HbtTrack - we are done with that
  // Parent
  StTrk = SKV.parent();
  HbtTrk = new AliFemtoTrack(StTrk,PrimaryVertex); // generate NEW HbtTrack from StTrack
  mParent = *HbtTrk;                           // invoke copy ctr of AliFemtoTrack
  delete HbtTrk;                               // get rid of the NEW HbtTrack - we are done with that

  // finally, the kink position
  mPosition.setX(SKV.position().x());
  mPosition.setY(SKV.position().y());
  mPosition.setZ(SKV.position().z());

}

// mike removed all AliFemtoTTree stuff 21apr2006

#endif // __ROOT__
#endif  // __NO_STAR_DEPENDENCE_ALLOWED__
