/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id:  */

//_________________________________________________________________________
//  A singleton that retrieves objets from an array stored in a Tree on a disk file
//    1. AliPHOSDigit from TreeD     
//                  
//*-- Author: Yves Schutz (SUBATECH)
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TTree.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSIndexToObject.h"

ClassImp(AliPHOSIndexToObject)
  
  AliPHOSIndexToObject * AliPHOSIndexToObject::fgObjGetter = 0 ; 

//____________________________________________________________________________ 
AliPHOSIndexToObject::AliPHOSIndexToObject(AliPHOS * det)
{
  // ctor called once to initialize the detector in use

  fDetector = det ; 
}

//____________________________________________________________________________ 
AliPHOSIndexToObject * AliPHOSIndexToObject::GetInstance()
{
  // Returns the pointer of the unique instance already defined
  
  AliPHOSIndexToObject * rv = 0 ;
  if ( fgObjGetter )
    rv = fgObjGetter ;
  else
    cout << "AliPHOSIndexToObject::GetInstance ERROR: not yet initialized" << endl ;
  
  return rv ;
}

//____________________________________________________________________________ 
AliPHOSIndexToObject * AliPHOSIndexToObject::GetInstance(AliPHOS * det)
{
  // Creates and returns the pointer of the unique instance
  // Must be called only when the environment has changed (a new event for exemple)

  if ( fgObjGetter )      // delete it if already exists
    delete fgObjGetter ; 

  fgObjGetter = new AliPHOSIndexToObject(det) ; 
  
  return fgObjGetter ; 
  
}

//____________________________________________________________________________ 
AliPHOSDigit * AliPHOSIndexToObject::GimeDigit(Int_t index)
{
  // returns the object AliPHOSDigit stored at array position index in TreeD

  AliPHOSDigit * rv = 0 ; 

  if ( index >= fDetector->Digits()->GetEntries() ) 
    cout << "AliPHOSIndexToObject::GimeDigit: index " << index << " larger than available entries " 
	 <<  fDetector->Digits()->GetEntries() << endl ; 
  else if ( index != -1) 
    rv =  (AliPHOSDigit *) (fDetector->Digits()->At(index) ) ; 

  return rv ;
  
}
//____________________________________________________________________________ 
TParticle * AliPHOSIndexToObject::GimePrimaryParticle(Int_t index)
{
  // returns the object TParticle stored at array position index in TreeK

  TParticle * rv = 0 ; 
     
  if ( index >= gAlice->Particles()->GetEntries() ) 
    cout << "AliPHOSIndexToObject::GimePrimaryParticles: index " << index << " larger than available entries " 
	 <<  gAlice->Particles()->GetEntries() << endl ; 
  else 
    rv =  (TParticle *) (gAlice->Particles()->At(index) ) ; 

  return rv ;
  
}

//____________________________________________________________________________ 
AliPHOSRecParticle * AliPHOSIndexToObject::GimeRecParticle(Int_t index)
{
  // returns the object AliPHOSRecParticle stored at array position index in TreeR/PHOSRP
  // this one takes more work because the detetor object and the objects in TreeR are not saved at the same time
  // therefore the links are lost

  AliPHOSRecParticle * rv = 0 ; 

  AliPHOSRecParticle::RecParticlesList * rplist = *(fDetector->RecParticles()) ; 

  Int_t rpentries  = 0 ; 

  if (rplist) 
    rpentries = rplist->GetEntries() ;
  
  fReconstruct = gAlice->TreeR() ; 
  
  if (!rpentries) {
    fReconstruct->SetBranchAddress( "PHOSRP", &rplist ) ;
    fReconstruct->GetEvent(0) ;
    rpentries = rplist->GetEntries() ;  
  }     
  
  if ( index >= rpentries )  // ERROR 
    cout << "AliPHOSIndexToObject::GimeRecParticle: index " << index << " larger than available entries " 
	   <<  rpentries << endl ; 
  else 
    rv =  (AliPHOSRecParticle *) (*(fDetector->RecParticles()) )->At(index)  ; 
  
  return rv ;
  
}

//____________________________________________________________________________ 
AliRecPoint * AliPHOSIndexToObject::GimeRecPoint(Int_t index, TString type)
{
  // returns the object AliPHOSRecPoint stored at array position index in TreeR/PHOSEmcRP or TreeR/PHOSPpsdRP
  // this one takes more work because the detetor object and the objects in TreeR are not saved at the same time
  // therefore the links are lost

  AliPHOSRecPoint * rv = 0 ; 
  
  AliPHOSRecPoint::RecPointsList * emclist = *(fDetector->EmcRecPoints() ); 
  AliPHOSRecPoint::RecPointsList * ppsdlist = *(fDetector->PpsdRecPoints() ); 

  Int_t emcentries  = 0 ; 
  Int_t ppsdentries = 0 ; 

  if (emclist) 
    emcentries = emclist->GetEntries() ;

  if (ppsdlist)
    ppsdentries= ppsdlist->GetEntries() ;

  fReconstruct = gAlice->TreeR() ; 
  
  if (!emcentries && !ppsdentries) {
    fReconstruct->SetBranchAddress("PHOSEmcRP",&emclist);
    fReconstruct->SetBranchAddress("PHOSPpsdRP",&ppsdlist);
    fReconstruct->GetEvent(0) ;
    emcentries = emclist->GetEntries() ;
    ppsdentries= ppsdlist->GetEntries() ;
  }     

  if ( type == "emc" ) {
    if ( index >= emcentries ) 
      cout << "AliPHOSIndexToObject::GimeRecPoint emc: index " << index << " larger than available entries " 
	   <<  emcentries << endl ; 
    else 
      rv =  (AliPHOSEmcRecPoint *) ( emclist->At(index) ) ;
  } 
  else if ( type == "ppsd" ) {  
    if ( index >= ppsdentries ) 
      cout << "AliPHOSIndexToObject::GimeRecPoint ppsd: index " << index << " larger than available entries " 
	   <<  ppsdentries << endl ; 
    else if (index != -1) 
      rv =  (AliPHOSPpsdRecPoint *) (ppsdlist->At(index) ) ;
  } else
    cout << "AliPHOSIndexToObject::GimeRecPoint: " << type << " is an unknown type " << endl
	 << " valid types are : emc " << endl 
	 << "                   ppsd " << endl ;
  
  return rv ;
  
}

//____________________________________________________________________________ 
AliPHOSTrackSegment * AliPHOSIndexToObject::GimeTrackSegment(Int_t index)
{
  // returns the object AliPHOSTrackSegment stored at array position index in TreeR/PHOSTS
  // this one takes more work because the detetor object and the objects in TreeR are not saved at the same time
  // therefore the links are lost

  AliPHOSTrackSegment * rv = 0 ; 

  AliPHOSTrackSegment::TrackSegmentsList * tslist = *( fDetector->TrackSegments()) ; 

  Int_t tsentries  = 0 ; 

  if (tslist) 
    tsentries = tslist->GetEntries() ;
  
  fReconstruct = gAlice->TreeR() ; 
  
  if (!tsentries) {
    fReconstruct->SetBranchAddress( "PHOSTS", &tslist ) ;
    fReconstruct->GetEvent(0) ;
    tsentries = tslist->GetEntries() ;  
  }     
  
  if ( index >= tsentries )  // ERROR 
      cout << "AliPHOSIndexToObject::GimeTrackSegment: index " << index << " larger than available entries " 
	   <<  tsentries << endl ; 
  else 
    rv =  (AliPHOSTrackSegment *) (tslist->At(index) ) ; 
  
  return rv ;
  
}
