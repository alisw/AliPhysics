
// Flugg tag 

////////////////////////////////////////////////////////////////////
//
// WrapNorml.hh - Sara Vanini 
//
// Wrapper for computing normal unit-vector in global coordinates.
//
// Fluka requires normal vector exiting from final position (of 
// particle) volume, that is: entering in volume of initial position. 
// Geant4 always computes the normal vector exiting from the volume. 
// In GetLocalExitNormal() call the volume is the pre-step volume 
// (so G4 normal vector sign is opposite of fluka-required normal).
// If IsLocalExitNormalValid=false, normal value is computed from 
// init-step volume (in this case sign must change), or from
// end-step volume (the sign is the same). Normal vector is computed 
// always on boundary between volumes, in global coordinates (to take
// rotation of parameterised volumes in hierarchy in consideration). 
// So: nrmlwr returns inwards pointing unit normal of the shape for 
// surface closest to the point returned by the navigator (last step 
// end-point).
// 
// modified 10/III/99
// modified 25/V/00
// modified 7/VI/00 for boundary-crossing in case of relocation
// modified 5/VII/00 geometry error on boundary fixed
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
////////////////////////////////////////////////////////////////////

#include "Wrappers.hh"
#include "WrapUtils.hh"
#include "FGeometryInit.hh"
#include "G4VPhysicalVolume.hh"
#include "FluggNavigator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"


void nrmlwr(G4double& pSx, G4double& pSy, G4double& pSz,
            G4double& pVx, G4double& pVy, G4double& pVz,
	    G4double* norml, const G4int& oldReg, 
	    const G4int& newReg, G4int& flagErr)
{
  //flag
#ifdef G4GEOMETRY_DEBUG
  G4cout << "============ NRMLWR-DBG =============" << G4endl;
#endif
  
  //dummy variables
  flagErr=0;
  
  //navigator pointer
  static FGeometryInit * ptrGeoInit;
  ptrGeoInit = FGeometryInit::GetInstance();
  FluggNavigator * ptrNavig = ptrGeoInit->getNavigatorForTracking();
  
  //variables
  G4ThreeVector normalLoc;
  G4ThreeVector normalGlob;
  
  //normal computing
  G4bool* valid;
  normalLoc=ptrNavig->GetLocalExitNormal(valid);
  if(valid) {
    normalLoc *= -1;        
    
    //global cooordinates normal
    normalGlob = ptrNavig->GetLocalToGlobalTransform().
      TransformAxis(normalLoc);
  }
  else {
    G4VPhysicalVolume *touchvolume;
    G4ThreeVector theLocalPoint;
    
    if(ptrNavig->EnteredDaughterVolume()) {
      //volume from touchable history
      G4TouchableHistory *ptrTouchableHistory=ptrGeoInit->
	GetTouchableHistory();
      touchvolume = ptrTouchableHistory->GetVolume();
      
      //local point from navigator and normal
      theLocalPoint = ptrNavig->GetCurrentLocalCoordinate(); 
      normalLoc = touchvolume->GetLogicalVolume()->GetSolid()->
	SurfaceNormal(theLocalPoint);
      
      //global cooordinates normal
      normalGlob = ptrNavig->GetLocalToGlobalTransform().
	TransformAxis(normalLoc);
    }
    else {
      //volume from old history
      const G4NavigationHistory * ptrOldNavHist = 
	ptrGeoInit->GetOldNavHist()->GetHistory();
      touchvolume = ptrOldNavHist->GetTopVolume();
      
      if(!touchvolume) {
	// Old history has been reseted by LOOKZ relocation,
	// so is necessary to track back and forward to find
	// the right histories.
	
	
	////////////  COMPUTE STEP BACKWARD //////////////////
	G4ThreeVector theGlobalPoint = ptrNavig->GetLocalToGlobalTransform().
	  TransformPoint(ptrNavig->GetCurrentLocalCoordinate());
	
	//compute step and new location
	G4ThreeVector pVec(pVx,pVy,pVz);
	pVec=-pVec;
	G4double appStep = 0;
	G4double safety = 0;
	G4bool onBoundary = false;
	G4double physStep = 1000000000;
	G4int newRegStep;
	G4ThreeVector partLoc =  theGlobalPoint;
	G4bool fErr=false;
	
#ifdef G4GEOMETRY_DEBUG
	G4cout << "Old history not found" << G4endl;
	G4cout << "* NRML needs boundary-crossing: computing step backward..."
	       << G4endl;
#endif
	
	//compute step and location 
	newRegStep=StepAndLocation(partLoc,pVec,physStep,
				   appStep,safety,onBoundary,
				   fErr,oldReg);
	
	G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::
	  GetInstance();
	
	if(appStep<physStep && newRegStep!=G4int(pVolStore->size())+1) {
	  //end-step point is on boundary between volumes;
	  //==> update step-histories
	  
#ifdef G4GEOMETRY_DEBUG
	  G4cout << "* updating step-histories" << G4endl;
#endif
	  
	  G4TouchableHistory * ptrTouchHist = 
	    ptrGeoInit->GetTouchableHistory();
	  ptrGeoInit->UpdateHistories(ptrTouchHist->GetHistory(),2);
	}
	else {
#ifdef G4GEOMETRY_DEBUG
	  G4cout << "ERROR! Boundary not found" << G4endl;
#endif
	}
	
#ifdef G4GEOMETRY_DEBUG
	G4cout << "* computing step forward..." << G4endl;
#endif
	pVec=-pVec;
	safety = 0;
	onBoundary = false;
	
	//compute step and location for boundary crossing
	newRegStep=StepAndLocation(partLoc,pVec,physStep,
				   appStep,safety,onBoundary,fErr,oldReg);
	if(appStep<physStep) {
	  //end-step point is on boundary between volumes;
	  //==> update step-histories
	  
#ifdef G4GEOMETRY_DEBUG
	  G4cout << "* updating step-histories" << G4endl;
#endif
	  
	  G4TouchableHistory * ptrTouchHist = 
	    ptrGeoInit->GetTouchableHistory();
	  ptrGeoInit->UpdateHistories(ptrTouchHist->GetHistory(),2);
	}
      }
      
      // now touchvolume exist.
      // local point from navigator and global point
      // N.B. if particle has exited world volume, 
      // FluggNavigator doesn't update lastLocatedPoint. 
      // So be carefull in building geometry always to have a 
      // big world volume that fluka won't exit.
      
      touchvolume = ptrOldNavHist->GetTopVolume();
      
      G4ThreeVector theGlobalPoint = ptrNavig->GetLocalToGlobalTransform().
	TransformPoint(ptrNavig->GetCurrentLocalCoordinate());
      theLocalPoint = ptrOldNavHist->GetTopTransform().
	TransformPoint(theGlobalPoint);
      normalLoc = (touchvolume->GetLogicalVolume()->GetSolid()->
		   SurfaceNormal(theLocalPoint));
      normalLoc *= -1; 
      
      //global cooordinates normal
      normalGlob = ptrOldNavHist->GetTopTransform().
	Inverse().TransformAxis(normalLoc);	    
    }
  }
  
  //return normal:
  norml[0]=G4double(normalGlob.x());
  norml[1]=G4double(normalGlob.y());
  norml[2]=G4double(normalGlob.z());
  
  //for debugging
#ifdef G4GEOMETRY_DEBUG
  G4cout << "Normal: " << normalGlob << G4endl;
#endif
}


