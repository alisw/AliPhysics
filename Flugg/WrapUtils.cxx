#include "WrapUtils.hh"
#include "FGeometryInit.hh"

////////////////////////////////////////////////////////////////////////
// StepAndLocation
////////////////////////////////////////////////////////////////////////

G4int StepAndLocation(G4ThreeVector &partLoc, const G4ThreeVector &pVec,
		      const G4double &physStep, G4double &appStep,
		      G4double &safety, G4bool & onBound, 
		      G4bool & flagErr, const G4int & oldReg)
{
  //NB onBound=true in the particular case when particle is in boundary 
  //rounding BUT geometry hasn't crossed it.
  //flagErr=true when particle is on boundary but step=infinity, 
  //newSafety>10.e-10 and lLocate function locates end step point in
  // a new region. This is absurb!
  
  
  //Geoinit, Navigator, TouchableHistory, etc. pointers
  static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
  FluggNavigator * ptrNavig = ptrGeoInit->getNavigatorForTracking();
  G4TouchableHistory * ptrTouchableHistory = ptrGeoInit->GetTouchableHistory();
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  
  //compute step
  G4double step = ptrNavig->ComputeStep(partLoc,pVec,physStep,safety);
  
  //compute approved step
  if(step<physStep) {
    //NB:to check if point is really on boundary: 
    //see G4AuxiliaryNavServices.hh
#ifdef G4GEOMETRY_DEBUG
    G4cout << "* Boundary crossing!" << G4endl; 
#endif
    
    appStep=step;
    ptrNavig->SetGeometricallyLimitedStep();
    
    //if SetGeometricallyLimitedStep() is called, next 
    //LocateGlobalPointAndUpdateTouchable(...,true) will 
    //start searching from the end of the previous step, 
    //otherwise will start from last located point. 
    //(e.g.PropagatorInField and SteppingManger.
  }
  else {
    //step=kInfinity (the nearest boundary is >physStep)
    //Fluka always wants step lenght approved, so:
    appStep=physStep;
  }
  
  //new position
  partLoc+=(pVec*appStep);
  
  //locate point and update touchable history
  ptrNavig->LocateGlobalPointAndUpdateTouchable(partLoc,0,
						ptrTouchableHistory,true);
  
  G4VPhysicalVolume * touchvolume=ptrTouchableHistory->GetVolume();
  
  //if volume not found, out of mother volume: 
  //returns [number of volumes]+1
  G4int newReg;
  if(!touchvolume) {
#ifdef G4GEOMETRY_DEBUG
    G4cout << "Out of mother volume! (G1WR)"  << G4endl;
#endif
    G4int numVol = G4int(pVolStore->size());
    newReg = numVol + 1;
  }
  else {
#ifdef G4GEOMETRY_DEBUG
    G4cout <<"Volume after step: " <<touchvolume->GetName() << G4endl;
#endif
    
    //compute new volume index
    G4int volIndex = ~0;
    for (unsigned int i=0; i<pVolStore->size(); i++)
      if ((*pVolStore)[i] == touchvolume) 
	volIndex = i;
    if (volIndex==(~0)) {
      G4cerr << "ERROR in FLUGG: Problem in routine WrapG1 trying to find volume after step" << G4endl;
      exit(-999);
    }
    newReg=volIndex+1;   
  }
  
  onBound=false;
  flagErr=false;
  //check if position is in a boundary rounding while volume isn't changed
  if(step==kInfinity || step==physStep) {
    G4ThreeVector globalpoint = ptrNavig->GetLocalToGlobalTransform().
      TransformPoint(ptrNavig->GetCurrentLocalCoordinate());
    
    //compute new safety
    G4double newSafetydbg = ptrNavig->ComputeSafety(globalpoint,DBL_MAX );
    
#ifdef G4GEOMETRY_DEBUG
    G4cout << "|From StepAndLocation:"  << G4endl;
    G4cout << "|step=" << step << G4endl;
    G4cout << "|phystep=" << physStep << G4endl;
    G4cout << "|safety=" << safety << G4endl;
    G4cout << "|newSafetydbg =" << newSafetydbg << G4endl;
#endif
    
    if(newSafetydbg<1.0e-10) 
      onBound=true;
    else if(newReg != oldReg) { 
#ifdef G4GEOMETRY_DEBUG
      G4cout << "New Volume but ComputeStep didn't notice!" << G4endl;
#endif
      flagErr=true;
    } 
  }
  
  //return
  return newReg;
}




////////////////////////////////////////////////////////////////////////
// EqualHistories
////////////////////////////////////////////////////////////////////////

bool EqualHistories(const G4NavigationHistory* ptrFirstHist,
		    const G4NavigationHistory* ptrSecHist)
{
#ifdef G4GEOMETRY_DEBUG
  G4cout << "* Testing if histories are equal..."  << G4endl;
#endif
     
  G4int depth1 = ptrFirstHist->GetDepth();
  G4int depth2 = ptrSecHist->GetDepth();
  
  if(depth1!=depth2) 
    return false;
  
  for(G4int w=0;w<=depth1;w++) {
    if (*(ptrFirstHist->GetVolume(w))==*(ptrSecHist->GetVolume(w))) {
      //kNormal volume
      if(ptrFirstHist->GetVolumeType(w) == kNormal) {
	if ( ptrFirstHist->GetVolume(w)->GetCopyNo() !=
	     ptrSecHist->GetVolume(w)->GetCopyNo() )
	  return false;
      }
      
      //Replica or Parametric volume
      else {
	if ( ptrFirstHist->GetReplicaNo(w) !=
	     ptrSecHist->GetReplicaNo(w) )
	  return false;
      }
    }
    else
      return false;	
  }
  
#ifdef G4GEOMETRY_DEBUG
  G4cout << "Histories are equal!" << G4endl;	
#endif
  
  return true;
}



////////////////////////////////////////////////////////////////////////
// PrintHeader
////////////////////////////////////////////////////////////////////////
std::ostream& PrintHeader(std::ostream& os, const char* title) {
  os << "*\n" << "*\n" << "*\n";
  os << "*********************  " << title << " *********************\n"
     << "*\n";
  os << "*...+....1....+....2....+....3....+....4....+....5....+....6....+....7..."
     << G4endl;
  os << "*" << G4endl;

  return os;
}

////////////////////////////////////////////////////////////////////////
// ToFlukaString converts an string to something usefull in FLUKA:
// * Capital letters
// * Only 8 letters
// * Replace ' ' by '_'
////////////////////////////////////////////////////////////////////////
G4String ToFlukaString(G4String str) {
  str = str.substr(0,8);
  str.toUpper();
  for (unsigned int i = 0; i < str.length(); i++) {
    if (str.at(i) == ' ')
      str.replace(i,1,"_",1);
  }
  return str;
}
