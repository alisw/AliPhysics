
// Flugg tag 

///////////////////////////////////////////////////////////////////
//
// WrapLookMG.hh - Sara Vanini 26/X/99
//
// Wrapper for localisation of particle for magnetic field tracking 
//
// modified 13/IV/00: check if point belongs to one of the lattice 
// histories stored in jrLtGeant 
//
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
// modified 17/06/02: by I. Gonzalez. STL migration
//
///////////////////////////////////////////////////////////////////


#include "Wrappers.hh"
#include "FGeometryInit.hh"
#include "NavHistWithCount.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4NormalNavigation.hh"
#include "G4VoxelNavigation.hh"
#include "G4ParameterisedNavigation.hh"
#include "G4ReplicaNavigation.hh"
#include "globals.hh"


//auxiliary function declarations
G4bool PointLocate(const G4NavigationHistory &,
			  const G4ThreeVector &, 
			  const G4ThreeVector *,
			  G4int &);
EVolume CharacteriseDaughters(const G4LogicalVolume *);


void lkmgwr(G4double& pSx, G4double& pSy, G4double& pSz,
            G4double* pV, const G4int& oldReg, const G4int& oldLttc,
	    G4int& flagErr, G4int& newReg, G4int& newLttc)
{
//flag
#ifdef G4GEOMETRY_DEBUG
    G4cout<<"============= LKMGWR =============="<<G4endl;
#endif

    //Geoinit, Navigator, etc. pointers
    static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
    G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
    
    //setting variables (and dimension: Fluka uses cm.!)
    G4ThreeVector globalPointcm(pSx,pSy,pSz);
    G4ThreeVector globalPoint =  globalPointcm * 10.; //in mm.
    G4ThreeVector globalDirection(pV[0],pV[1],pV[2]);
    
    //get jrLtGeant array and initialize variables
    G4int * jrLtGeant = ptrGeoInit->GetJrLtGeantArray();
    G4int LttcFlagGeant = ptrGeoInit->GetLttcFlagGeant();
    G4bool belongToVolume = false;
    G4int i = LttcFlagGeant;
    G4int regionIn = 0;
    newReg = -1;
    newLttc = -1;
    
    while(!belongToVolume && i>=0) {
      //get history from jrLtGeant
      NavHistWithCount * ptrNavHistCount = 
	reinterpret_cast<NavHistWithCount*>(jrLtGeant[i]);
      const G4NavigationHistory * ptrNavHist =
	ptrNavHistCount->GetNavHistPtr();
      
      //check if globalPoint belongs to volume (can't call 
      //LocateGlobalPoint... because of flag settings)
      belongToVolume = PointLocate(*ptrNavHist,globalPoint,
				   &globalDirection,regionIn);
      
      
      //if point belongs to surface, check scalar product direction*normal
      if(regionIn==-100) {
#ifdef G4GEOMETRY_DEBUG
	G4cout<<"On surface!"<<G4endl;
#endif
	
	//if entering, then point belongs to volume, 
	//if exiting, point doesn't belong to volume.
	G4int oldReg,flagErr,reg;
	G4double x,y,z,px,py,pz;
	G4double * norml = new G4double[3];
	x = globalPoint.x();
	y = globalPoint.y();
	z = globalPoint.z();
	px = globalDirection.x();
	py = globalDirection.y();
	pz = globalDirection.z();
	
	nrmlwr(x,y,z,px,py,pz,norml,oldReg,reg,flagErr);
	
	G4ThreeVector normal(norml[0],norml[1],norml[2]);
#ifdef G4GEOMETRY_DEBUG
	// G4cout<<"Scalar product="<<globalDirection.dot(normal)<<G4endl;
#endif
	if(globalDirection.dot(normal)>0) {
#ifdef G4GEOMETRY_DEBUG
	  G4cout<<"entering volume!"<<G4endl;
#endif
	  G4int volIndex = ~0;
	  for (unsigned int i=0; i<pVolStore->size(); i++)
	    if ((*pVolStore)[i] == ptrNavHist->GetTopVolume()) 
	      volIndex = i;
	  if (volIndex==(~0)) {
	    G4cerr << "FLUGG: Problem in routine WrapG1 tryingto find volume after step" << G4endl;
	    exit(-999);
	  }
	  //regionIn = G4int(pVolStore->index(ptrNavHist->GetTopVolume()))+1;
	  regionIn = volIndex+1; 
	  belongToVolume = true;
	}
      }
      
      i -= 1;
    }
    
    //output variables
    if(belongToVolume) {
      newReg = regionIn;
      newLttc = jrLtGeant[i+1];
    }

    flagErr = pVolStore->size() + 2;      //flagErr=fluka region number + 1
    
    
#ifdef G4GEOMETRY_DEBUG
    G4cout << "Global point (cm): " << globalPointcm << G4endl;
    G4cout << "Direction: " << globalDirection << G4endl;
    if(newReg!=-1) 
      G4cout << "Point belongs to region " << newReg
	     << " ptr history=" << newLttc << G4endl;
    else 
      G4cout << "No containing volume found!" << G4endl;
#endif
}




//returns false if the point doesn't belong to history top volume (can either 
//belong to one of its daughter volumes or none), otherwise returns true.

G4bool PointLocate(const G4NavigationHistory & blockedHistConst,
		   const G4ThreeVector & globalPoint, 
		   const G4ThreeVector * pGlobalDirection,
		   G4int & reg)
{
  
  //variables
  G4NavigationHistory * fHistory = new G4NavigationHistory(blockedHistConst);
  G4bool belongVolume;
  
  //G4 flags resetted (see: ResetStackAndState)
  G4bool fExiting = false; 
  G4VPhysicalVolume * fBlockedPhysicalVolume=0;
  G4int fBlockedReplicaNo=-1;
  G4bool fLocatedOnEdge=false;  
  
  G4bool notKnownContained=true,noResult;
  G4VPhysicalVolume *targetPhysical;
  G4LogicalVolume *targetLogical;
  G4VSolid *targetSolid;
  G4ThreeVector localPoint,localDirection;
  EInside insideCode;
  
  
  //Helpers/Utility classes
  G4NormalNavigation  fnormalNav;
  G4VoxelNavigation fvoxelNav;
  G4ParameterisedNavigation fparamNav;
  G4ReplicaNavigation freplicaNav;	  
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();
  
  //local variables
  localPoint=fHistory->GetTopTransform().TransformPoint(globalPoint);
  localDirection=fHistory->GetTopTransform().TransformPoint(*pGlobalDirection);
  
  //search if volume contains point
  if (fHistory->GetTopVolumeType()!=kReplica) {
    targetSolid=fHistory->GetTopVolume()->GetLogicalVolume()->GetSolid();
    insideCode=targetSolid->Inside(localPoint);
  }
  else {
    insideCode=freplicaNav.BackLocate(*fHistory,globalPoint,
				      localPoint,fExiting,notKnownContained);
    // !CARE! if notKnownContained returns false then the point is within
    // the containing placement volume of the replica(s). If insidecode
    // will result in the history being backed up one level, then the
    // local point returned is the point in the system of this new level
#ifdef G4GEOMETRY_DEBUG
    G4cout << "Replica: fExiting=" << fExiting << G4endl;
    G4cout << "         notKnownContained=" << notKnownContained << G4endl;
#endif
  }
  
  if (insideCode==kOutside) 
    belongVolume = false;

  else if (insideCode==kSurface) {
    belongVolume = false;
    reg = -100;
  }
  else {
    // Search downwards in daughter volumes
    // until deepest containing volume found
    //
    // 3 Cases:
    //
    // o Parameterised daughters
    //   =>Must be one G4PVParameterised daughter & voxels
    // o Positioned daughters & voxels
    // o Positioned daughters & no voxels
    
    belongVolume = true;    
    noResult=true;  
    
    do {
      // Determine `type' of current mother volume
      targetPhysical=fHistory->GetTopVolume();
      targetLogical=targetPhysical->GetLogicalVolume();
      switch(CharacteriseDaughters(targetLogical)) {
      case kNormal:
	if (targetLogical->GetVoxelHeader()) {
	  noResult=fvoxelNav.LevelLocate(*fHistory,
					 fBlockedPhysicalVolume,
					 fBlockedReplicaNo,
					 globalPoint,
					 pGlobalDirection,
					 fLocatedOnEdge,
					 localPoint);
	}
	else {
	  noResult=fnormalNav.LevelLocate(*fHistory,
					  fBlockedPhysicalVolume,
					  fBlockedReplicaNo,
					  globalPoint,
					  pGlobalDirection,
					  fLocatedOnEdge,
					  localPoint);
	}
	break;
      case kReplica:
	noResult=freplicaNav.LevelLocate(*fHistory,
					 fBlockedPhysicalVolume,
					 fBlockedReplicaNo,
					 globalPoint,
					 pGlobalDirection,
					 fLocatedOnEdge,
					 localPoint);
	break;
      case kParameterised:
	noResult=fparamNav.LevelLocate(*fHistory,
				       fBlockedPhysicalVolume,
				       fBlockedReplicaNo,
				       globalPoint,
				       pGlobalDirection,
				       fLocatedOnEdge,
				       localPoint);
	break;
      } //switch
      
      // LevelLocate search in the first daughter level. 
      // LevelLocate returns noResult=true if it finds a daughter volume
      // in which globalPoint is inside (or on the surface). So point
      // doesn't belong only to mother volume ==> belongVolume=false.
      
      if (noResult) {
	// The blocked volume no longer valid - it was for another level
	fBlockedPhysicalVolume= 0;
	fBlockedReplicaNo= -1;
	belongVolume=false;
      }
    } while (noResult);
    
    //on exit targetPhysical is the volume globalPoint belongs to;
    G4int volIndex = ~0;
    for (unsigned int i=0; i<pVolStore->size(); i++)
      if ((*pVolStore)[i] == targetPhysical) 
	volIndex = i;     
    if (volIndex==(~0)) {
      G4cerr << "FLUGG: Problem in routine WrapG1 tryingto find volume after step" << G4endl;
      exit(-999);
    }
    //reg = G4int(pVolStore->index(targetPhysical))+1;
    reg = volIndex+1;
    
  }
  
  delete fHistory;
  return belongVolume;
}

EVolume CharacteriseDaughters(const G4LogicalVolume *pLog)
{
  EVolume type;
  EAxis axis;
  G4int nReplicas;
  G4double width,offset;
  G4bool consuming;
  G4VPhysicalVolume *pVol;
  
  if (pLog->GetNoDaughters()==1) {
    pVol=pLog->GetDaughter(0);
    if (pVol->IsReplicated()) {
      pVol->GetReplicationData(axis,nReplicas,width,offset,consuming);
      type=(consuming) ? kReplica : kParameterised;
    }
    else
      type=kNormal;
  }
  else
    type=kNormal;
  return type;
}



