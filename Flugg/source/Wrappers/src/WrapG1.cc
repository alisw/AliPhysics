// $Id$
// Flugg tag $Name$

///////////////////////////////////////////////////////////////////
//
// WrapG1.hh - Sara Vanini
//
// Wrapper for geometry tracking: returns approved step of 
// particle and alla variables that fluka G1 computes.
//
// modified 11/III/99 : included jrLtGeant array for storing 
// lattice histories; fixed jump-on-boundaries and back-scattering
// modified 20/IV/99 : history counters of jrLt array are
// incremented when created and decremented when check is required:
// this for not deleting histories twice.
// modified 18/X/99 : LocateGlobalPointWithinVolume used when position
// is changed from last located point.
// modified 1/III/00 : update utilities histories when crossing 
// identical volume boundaries.
//
// modified 22/III/00 : fixed LttcFlag and jrLt return values.
// modified 12/VI/00 : end-step on Boundary bug fixed.
// modified 5/VII/00 : boundary not seen by G4 geometry bug fixed.
// modified 24.10.01: by I. Hrivnacova
//   functions declarations separated from implementation
//   (moved to Wrappers.hh);
//
/////////////////////////////////////////////////////////////////////

#include "Wrappers.hh"
#include "FGeometryInit.hh"
#include "NavHistWithCount.hh"
#include "G4VPhysicalVolume.hh"
#include "G4NavigationHistory.hh"
#include "G4Navigator.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalVolumeStore.hh"
#include "globals.hh"


//StepAndLocation declaration
extern G4int StepAndLocation(G4ThreeVector &, const G4ThreeVector &, 
                             const G4double &, G4double &, G4double &, G4bool &,
			     G4bool &, const G4int &);

//EqualHistories declaration: true if histories are identical, otherwise false
extern bool EqualHistories(const G4NavigationHistory*, 
                           const G4NavigationHistory*);



void g1wr(G4double& pSx, G4double& pSy, G4double& pSz, G4double* pV,
          G4int& oldReg, const G4int& oldLttc, G4double& propStep,
          G4int& nascFlag, G4double& retStep, G4int& newReg,
          G4double& saf, G4int& newLttc, G4int& LttcFlag,
          G4double* sLt, G4int* jrLt)
{
//flag
#ifdef G4GEOMETRY_DEBUG
  G4cout<<"============= G1WR =============="<<G4endl;    
#endif 

//static G4int count=0;
//count+=1;
//G4cout<<"contatore G1="<<count<<G4endl;

///////////////////////// INPUT ///////////////////////////
//Geoinit, Navigator, TouchableHistory pointers
  static FGeometryInit * ptrGeoInit = FGeometryInit::GetInstance();
  G4Navigator * ptrNavig = ptrGeoInit->getNavigatorForTracking();
  G4TouchableHistory * ptrTouchHist = ptrGeoInit->GetTouchableHistory();
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::GetInstance();

//setting variables (and dimension: Fluka uses cm.!)
  G4ThreeVector partLoc(pSx,pSy,pSz);
  partLoc *= 10.0; // in millimeters!
  static G4ThreeVector partLocOld = partLoc;
  static G4ThreeVector oldLocalPoint = 
      ptrNavig->ComputeLocalPoint(partLocOld);

  G4ThreeVector pVec(pV[0],pV[1],pV[2]);
  const G4double physStep=G4double(propStep*10.);

  G4int LttcFlagGeant = ptrGeoInit->GetLttcFlagGeant();


#ifdef G4GEOMETRY_DEBUG
  G4cout.precision(10);
  G4cout<<"Position (cm):"<<pSx<<","<<pSy<<","<<pSz<<G4endl;
  G4cout<<"Direction: "<<pVec<<G4endl;
  G4cout<<"Proposed step :"<<propStep<<G4endl;
#endif 



///////////////////// FLUKA/G4 REGIONS COMPARISON //////////////////////
//get oldReg pointer and G4volume of previous step end-point pointer
  G4VPhysicalVolume * ptrOldReg = (*pVolStore)[oldReg-1];
  G4VPhysicalVolume * lastVolume = ptrTouchHist->GetVolume();
  G4int fluVolCpNb = ptrOldReg->GetCopyNo();
  G4int touVolCpNb = 0;
  if (lastVolume) touVolCpNb = lastVolume->GetCopyNo();

#ifdef  G4GEOMETRY_DEBUG
  G4cout<<"Fluka volume before step: "<<ptrOldReg->GetName()<<
    ","<<fluVolCpNb<<G4endl;
  if(lastVolume) G4cout<<"G4 Touch Hist volume: "<<
  		   lastVolume->GetName()<<","<<touVolCpNb<<G4endl;
  G4cout <<"------------------------------------------------"<<G4endl;
#endif 


//if volume is changed, this is a new particle tracking, or fluka tries
//to reach a boundary more softly with the same particle. In this case
//fluka restart tracking from old history, in general. For tracking in 
//lattice-volume fluka goes back to one of the previous lattice volumes. 
//Check if ptrOldReg is equal to old history, or to one of the N lattice 
//volumes stored in jrLt. Then reinitialise step-histories! Otherwise 
//must relocate.
//NB. jrLtGeant stores lattice volume histories until LttcFlag==-1, 
//then all histories are checked and deleted, and array is reinitialised 


  G4int haveHistNb = -1;
  G4int newRegErr=0;
  G4int indHist = LttcFlagGeant;
  G4int * jrLtGeant = ptrGeoInit->GetJrLtGeantArray();


  G4NavigationHistory* ptrLttcHist;
  if(oldLttc) ptrLttcHist = reinterpret_cast
    <NavHistWithCount*>(oldLttc)->GetNavHistPtr();


  while(indHist>=0 && haveHistNb==-1)
    {
#ifdef  G4GEOMETRY_DEBUG
      G4cout<<"Searching in jrLtc...."<<G4endl;
#endif 
      if(oldLttc==jrLtGeant[indHist]) haveHistNb=indHist;
      indHist-=1;
    }

   if(haveHistNb!=-1)
    {
      //fluka history found in jrLtGeant
      if(haveHistNb<LttcFlagGeant)
	{
#ifdef  G4GEOMETRY_DEBUG
	  G4cout<<"* Fluka reaches boundary more softly..."<<G4endl;
	  G4cout<<"* Re-initializing G4Navigator history"<<G4endl;
	  G4cout<<"* and updating step-histories"<<G4endl;
#endif 

	  ptrNavig->UpdateNavigatorHistory(ptrLttcHist);
	  ptrGeoInit->UpdateHistories(ptrLttcHist,0);
	}
#ifdef  G4GEOMETRY_DEBUG
      if(haveHistNb==LttcFlagGeant) G4cout<<"Continuing step...."<<G4endl;
#endif 
      jrLt[0]=oldLttc;
    }
  else
    {
      //not found fluka history in jrLttGeant!
      G4cout<<"* ERROR! Geometry not correctly initialised in fluka history!"
	  <<G4endl; 

      //relocation!
      ptrNavig->LocateGlobalPointAndUpdateTouchable
	(partLoc,0,ptrTouchHist,true);
	  
      G4cout<<"* ATTENTION: point relocation in: "<<
	ptrTouchHist->GetVolume()->GetName()<<G4endl;

      ptrGeoInit->UpdateHistories(ptrTouchHist->GetHistory(),1);
      
      if(ptrTouchHist->GetVolume() != ptrOldReg)
	{
	  G4cout<<"* ERROR! Point not in fluka volume!"<<G4endl;
	  newRegErr=-3;
	}

      //save new history in jrLt[0] and increment its counter 
#ifdef G4GEOMETRY_DEBUG
      G4cout<<"* ISVHWR call to store new NavHistWithCount in jrLt[0]"
	  <<G4endl;
#endif 
      jrLt[0]=isvhwr(0,0);

#ifdef G4GEOMETRY_DEBUG
      G4cout<<"* CONHWR call to increment counter"<<G4endl;
#endif 
      G4int incrCount2=1;
      conhwr(jrLt[0],&incrCount2);
    }


//jrLtGeant - history check: decrement counter and delete 
//histories, if LttcFlag=-1, then reinitialise array with -1. 
  if(LttcFlag==-1 && LttcFlagGeant>=0)
    {
#ifdef G4GEOMETRY_DEBUG
      G4cout<<"* CONHWR call to check and delete histories in jrLtGeant[]:"
	  <<G4endl;
#endif 
	    
      for(G4int ind=0; ind<=LttcFlagGeant; ind++)
	{
	  G4int incrCount1=-1;
	  if(jrLtGeant[ind]!=jrLt[0]) conhwr(jrLtGeant[ind],&incrCount1);
	  jrLtGeant[ind]=-1;
	} 
      LttcFlagGeant=-1;
    }


//update jrLt and sLt arrays   
  G4int end;
  if(LttcFlag>=0) end=LttcFlag;
  else end=99;
  // Added by A.Solodkov
  if (end>=100) {
      G4cout<<"Problems in WrapG1 routine"<<G4endl;
      G4cout<<"Index LttcFlag="<<end<<" is outside array bounds"<<G4endl;
      G4cout<<"Better to stop immediately !"<<G4endl;
      exit(1);
  }
  //jrLt re-initialization with -1 (jrLt[0] is already set)
  for(G4int vv=1;vv<=end;vv++) jrLt[vv]=-1;
  //sLt re-initialization
  for(G4int vs=0;vs<=end;vs++) sLt[vs]=0;

  LttcFlag=0;


/////////////////////  COMPUTE STEP  ////////////////////////
//update Navigator private flags and voxel stack if point is 
//changed from last located point (otherwise troubles come 
//when fluka changes location or particle because G4 computes 
//from last located point).
  G4ThreeVector newLocalPoint = ptrNavig->ComputeLocalPoint(partLoc);
  G4double moveLenSq = (newLocalPoint-oldLocalPoint).mag2();
  if(moveLenSq>=kCarTolerance*kCarTolerance) 
	ptrNavig->LocateGlobalPointWithinVolume(partLoc);


//compute step and new location
  newReg = oldReg;
  G4double appStep = 0;
  G4double safety = 0;
  G4bool onBoundaryRound = false;
  G4bool crossBound = false;
  G4double physStepTmp = physStep;
  G4bool flagError = false;

  while( (oldReg==newReg && appStep<physStepTmp) || onBoundaryRound )
    {
#ifdef G4GEOMETRY_DEBUG
      G4cout<<"* Computing step..."<<G4endl;
#endif
      
      //update variables
      oldReg = newReg;

      if(onBoundaryRound) 
	{
	  physStepTmp=10.e-10;
	  //compute step and location: returns newReg
	  newReg=StepAndLocation(partLoc,pVec,physStepTmp,appStep,
			     safety,onBoundaryRound,flagError,oldReg);
	  physStepTmp=0.;
	  crossBound=true;
	}
      else 
	{
	  physStepTmp -= appStep;
	  //compute step and location: returns newReg
	  newReg=StepAndLocation(partLoc,pVec,physStepTmp,appStep,
			     safety,onBoundaryRound,flagError,oldReg);
	}
      
      G4bool EqualHist;
      EqualHist = EqualHistories(ptrTouchHist->GetHistory(),
                        ptrGeoInit->GetTempNavHist()->GetHistory());

      if(!EqualHist && flagError)
	{
	  pVec=-pVec;
          newReg=StepAndLocation(partLoc,pVec,physStepTmp,appStep,
                             safety,onBoundaryRound,flagError,oldReg); 
          pVec=-pVec;
  	  physStepTmp+=1;
          newReg=StepAndLocation(partLoc,pVec,physStepTmp,appStep,
                             safety,onBoundaryRound,flagError,oldReg);

	  EqualHist = EqualHistories(ptrTouchHist->GetHistory(),
                        ptrGeoInit->GetTempNavHist()->GetHistory()); 
	}
         	
      //update sLt
      G4double pas = G4double(appStep);
      sLt[LttcFlag] += pas/cm;
      safety = (oldReg!=newReg)?0.:safety; 

      if(!EqualHist)
	{
	  //end-step point is on boundary between volumes;
	  //==> update step-histories, save new NavHistWithCount 
	  //and save its pointer in jrLt; save step in sLt.

	  //set onBoundaryRound=false to avoid re-compute step!
	  onBoundaryRound=false;

#ifdef G4GEOMETRY_DEBUG
	  G4cout<<"* History is changed!"<<G4endl;
          G4cout<<"* updating step-histories, jrLt, LttcFlag"<<G4endl;
#endif
      
	  ptrGeoInit->UpdateHistories(ptrTouchHist->GetHistory(),2);

	  LttcFlag += 1;

#ifdef G4GEOMETRY_DEBUG
	  G4cout<<"* ISVHWR call to store new NavHistWithCount in jrLt"<<G4endl;
#endif 

	  jrLt[LttcFlag] = isvhwr(0,0);

          #ifdef G4GEOMETRY_DEBUG
	  G4cout<<"* CONHWR call to increment counter"<<G4endl;
          #endif 
	  G4int incrCount3=1;
	  conhwr(jrLt[LttcFlag],&incrCount3);

	  sLt[LttcFlag] = sLt[LttcFlag-1];
	}
    }

//////////////////////   OUTPUT   //////////////////////////
//If back-scattering occured, and fluka is in the wrong region, return -3.
//(N. B. Boundary between replicans are not seen when physStep=distance 
//particle-boundary: in this case step=kInfinity and history is unchanged. 
//Following step is =0 and then history changes.) 
   if(nascFlag<0 && !appStep && physStep && newReg!=oldReg && !crossBound)
	{
	//don't need to compare histories because boundary between
	//identical volumes in different replicans are not seen
#ifdef G4GEOMETRY_DEBUG
   	G4cout <<"* Back-scattering!"<<G4endl;
#endif

	newReg=-3;
	}
    else if(newRegErr<0) newReg=newRegErr;


//compute output variables (in cm.!)
   //final step
   retStep = sLt[LttcFlag];
   
   //safety (Fluka sottracts a bit to safety to be sure 
   //not to jump on a boundary)
   G4double s = G4double(safety);
   s -= s*3.0e-09;
   saf=s/cm; 

//update wrapper utility variables
   //copy jrLt in jrLtGeant
   G4int start=0;
   if(haveHistNb!=-1 && LttcFlagGeant!=-1) start=1;
   for(G4int lt=start;lt<=LttcFlag;lt++)
      jrLtGeant[LttcFlagGeant+1+lt-start]=jrLt[lt];
   LttcFlagGeant+=(1+LttcFlag-start);
   newLttc = jrLt[LttcFlag];
   ptrGeoInit->SetLttcFlagGeant(LttcFlagGeant);

   partLocOld=partLoc;
   oldLocalPoint = ptrNavig->ComputeLocalPoint(partLocOld);

//compute new position
   G4ThreeVector oldPos = G4ThreeVector(pSx,pSy,pSz);
   G4ThreeVector newPos = oldPos + retStep*pVec;

#ifdef G4GEOMETRY_DEBUG
   G4cout<<"New position: "<<newPos<<G4endl;
   G4cout<<"Output region: "<<newReg<<G4endl;
   G4cout<<"G4 safety (cm): "<<(safety*0.1)<<G4endl;
   G4cout<<"Fluka safety (cm): "<<saf<<G4endl;
   G4cout<<"Approved step: "<<retStep<<G4endl;
   G4cout<<"LttcFlag = "<<LttcFlag<<G4endl;
   for(G4int i=0;i<=LttcFlag+1;i++)
     {
       G4cout<<"jrLt["<<i<<"]="<<jrLt[i]<<G4endl;
       G4cout<<"sLt["<<i<<"]="<<sLt[i]<<G4endl;
     }

   G4cout<<"LttcFlagGeant ="<<LttcFlagGeant<<G4endl;
   for(G4int ib=0;ib<=LttcFlagGeant+1;ib++)
     {
       G4cout<<"jrLtGeant["<<ib<<"]="<<jrLtGeant[ib]<<G4endl;
     }
   G4cout<<"newLttc="<<newLttc<<G4endl;
#endif
}




extern G4int StepAndLocation(G4ThreeVector &partLoc, const G4ThreeVector &pVec,
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
  G4Navigator * ptrNavig = ptrGeoInit->getNavigatorForTracking();
  G4TouchableHistory * ptrTouchableHistory = ptrGeoInit->
						GetTouchableHistory();
  G4PhysicalVolumeStore * pVolStore = G4PhysicalVolumeStore::
							GetInstance();
  
  //compute step
  G4double step=ptrNavig->ComputeStep(partLoc,pVec,physStep,safety);

  //compute approved step
  if(step<physStep)
    {
      //NB:to check if point is really on boundary: 
      //see G4AuxiliaryNavServices.hh
#ifdef G4GEOMETRY_DEBUG
      G4cout<<"* Boundary crossing!"<<G4endl; 
#endif

      appStep=step;
      ptrNavig->SetGeometricallyLimitedStep();

      //if SetGeometricallyLimitedStep() is called, next 
      //LocateGlobalPointAndUpdateTouchable(...,true) will 
      //start searching from the end of the previous step, 
      //otherwise will start from last located point. 
      //(e.g.PropagatorInField and SteppingManger.
    }
  else
    {
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
  if(!touchvolume)
    {
#ifdef G4GEOMETRY_DEBUG
      G4cout<<"Out of mother volume! (G1WR)"<<G4endl;
#endif
      G4int numVol = G4int(pVolStore->size());
      newReg = numVol + 1;
    }
  else
    {
#ifdef G4GEOMETRY_DEBUG
      G4cout <<"Volume after step: "<<touchvolume->GetName()<<G4endl;
#endif

      //compute new volume index
      G4int volIndex = 0;
      for (G4int i=0; i<pVolStore->size(); i++)
        if ((*pVolStore)[i] == touchvolume) volIndex = i;
      // G4int volIndex=G4int(pVolStore->index(touchvolume));
      newReg=volIndex+1;   
    }

  onBound=false;
  flagErr=false;
  //check if position is in a boundary rounding while volume isn't changed
  if(step==kInfinity || step==physStep)
    {
      G4ThreeVector globalpoint = ptrNavig->
	GetLocalToGlobalTransform().TransformPoint
	(ptrNavig->GetCurrentLocalCoordinate());

      //compute new safety
      G4double newSafetydbg = ptrNavig->ComputeSafety(globalpoint,DBL_MAX );

#ifdef G4GEOMETRY_DEBUG
      G4cout<<"|From StepAndLocation:"<<G4endl;
      G4cout<<"|step="<<step<<G4endl;
      G4cout<<"|phystep="<<physStep<<G4endl;
      G4cout<<"|safety="<<safety<<G4endl;
      G4cout<<"|newSafetydbg ="<<newSafetydbg<<G4endl;
#endif

      if(newSafetydbg<1.0e-10) onBound=true;
      else if(newReg != oldReg)
	 { 
#ifdef G4GEOMETRY_DEBUG
          G4cout<<"New Volume but ComputeStep didn't notice!"<<G4endl;
#endif
	  flagErr=true;
	 } 
    }

  //return
  return newReg;
}





extern bool EqualHistories(const G4NavigationHistory* ptrFirstHist,
			   const G4NavigationHistory* ptrSecHist)
   {
#ifdef G4GEOMETRY_DEBUG
     G4cout<<"* Testing if histories are equal..."<<G4endl;
#endif
     
     G4int depth1 = ptrFirstHist->GetDepth();
     G4int depth2 = ptrSecHist->GetDepth();

     if(depth1!=depth2) return false;
     
     for(G4int w=0;w<=depth1;w++)
       {
	 if (*(ptrFirstHist->GetVolume(w))==*(ptrSecHist->GetVolume(w)))
	   {
	     //kNormal volume
	     if(ptrFirstHist->GetVolumeType(w) == kNormal)
	       {
		 if ( ptrFirstHist->GetVolume(w)->GetCopyNo() !=
		      ptrSecHist->GetVolume(w)->GetCopyNo() )
		   return false;
	       }
	     
	     //Replica or Parametric volume
	     else 
	       {
		 if ( ptrFirstHist->GetReplicaNo(w) !=
		      ptrSecHist->GetReplicaNo(w) )
		   return false;
	       }
	   }
	 else return false;	
       }
     
#ifdef G4GEOMETRY_DEBUG
     G4cout<<"Histories are equal!"<<G4endl;	
#endif

     return true;
   }




