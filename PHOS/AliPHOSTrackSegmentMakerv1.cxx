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
/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.80  2006/08/28 10:01:56  kharlov
 * Effective C++ warnings fixed (Timur Pocheptsov)
 *
 * Revision 1.79  2006/04/25 12:41:15  hristov
 * Moving non-persistent data to AliESDfriend (Yu.Belikov)
 *
 * Revision 1.78  2005/11/18 13:04:51  hristov
 * Bug fix
 *
 * Revision 1.77  2005/11/17 23:34:36  hristov
 * Corrected logics
 *
 * Revision 1.76  2005/11/17 22:29:12  hristov
 * Faster version, no attempt to match tracks outside the PHOS acceptance
 *
 * Revision 1.75  2005/11/17 12:35:27  hristov
 * Use references instead of objects. Avoid to create objects when they are not really needed
 *
 * Revision 1.74  2005/07/08 14:01:36  hristov
 * Tracking in non-uniform nmagnetic field (Yu.Belikov)
 *
 * Revision 1.73  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Implementation version 1 of algorithm class to construct PHOS track segments
// Track segment for PHOS is list of 
//        EMC RecPoint + (possibly) CPV RecPoint
// To find TrackSegments we do the following: 
//  for each EMC RecPoints we look at
//   CPV RecPoints in the radious fRcpv. 
//  If there is such a CPV RecPoint, 
//   we make "Link" it is just indexes of EMC and CPV RecPoint and distance
//   between them in the PHOS plane. 
//  Then we sort "Links" and starting from the 
//   least "Link" pointing to the unassined EMC and CPV RecPoints assing them to 
//   new TrackSegment. 
// If there is no CPV RecPoint we make TrackSegment 
// consisting from EMC alone. There is no TrackSegments without EMC RecPoint.
//// In principle this class should be called from AliPHOSReconstructor, but 
// one can use it as well in standalone mode.
// Use  case:
//  root [0] AliPHOSTrackSegmentMakerv1 * t = new AliPHOSTrackSegmentMaker("galice.root", "tracksegmentsname", "recpointsname")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//               // reads gAlice from header file "galice.root", uses recpoints stored in the branch names "recpointsname" (default = "Default")
//               // and saves recpoints in branch named "tracksegmentsname" (default = "recpointsname")                       
//  root [1] t->ExecuteTask()
//  root [3] t->SetTrackSegmentsBranch("max distance 5 cm")
//  root [4] t->ExecuteTask("deb all time") 
//                 
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH) & Yves Schutz (SUBATECH) 
//

// --- ROOT system ---
#include "TTree.h"
#include "TBenchmark.h"

// --- Standard library ---
#include "Riostream.h"
// --- AliRoot header files ---
#include "AliPHOSGeometry.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSLink.h"
#include "AliPHOSGetter.h"
#include "AliESD.h"
#include "AliESDtrack.h"

ClassImp( AliPHOSTrackSegmentMakerv1) 


//____________________________________________________________________________
AliPHOSTrackSegmentMakerv1::AliPHOSTrackSegmentMakerv1() :
  AliPHOSTrackSegmentMaker(),
  fDefaultInit(kTRUE),
  fWrite(kFALSE),
  fNTrackSegments(0),
  fRcpv(0.f),
  fRtpc(0.f),
  fLinkUpArray(0),
  fEmcFirst(0),
  fEmcLast(0),
  fCpvFirst(0),
  fCpvLast(0),
  fModule(0),
  fTrackSegmentsInRun(0)
  
{
  // default ctor (to be used mainly by Streamer)
  InitParameters() ; 
}

//____________________________________________________________________________
AliPHOSTrackSegmentMakerv1::AliPHOSTrackSegmentMakerv1(const TString & alirunFileName, const TString & eventFolderName) :
  AliPHOSTrackSegmentMaker(alirunFileName, eventFolderName),
  fDefaultInit(kFALSE),
  fWrite(kFALSE),
  fNTrackSegments(0),
  fRcpv(0.f),
  fRtpc(0.f),
  fLinkUpArray(0),
  fEmcFirst(0),
  fEmcLast(0),
  fCpvFirst(0),
  fCpvLast(0),
  fModule(0),
  fTrackSegmentsInRun(0)
{
  // ctor
  InitParameters() ; 
  Init() ;
  fESD = 0;
}


AliPHOSTrackSegmentMakerv1::AliPHOSTrackSegmentMakerv1(const AliPHOSTrackSegmentMakerv1 & tsm) :
  AliPHOSTrackSegmentMaker(tsm),
  fDefaultInit(kFALSE),
  fWrite(kFALSE),
  fNTrackSegments(0),
  fRcpv(0.f),
  fRtpc(0.f),
  fLinkUpArray(0),
  fEmcFirst(0),
  fEmcLast(0),
  fCpvFirst(0),
  fCpvLast(0),
  fModule(0),
  fTrackSegmentsInRun(0)
{
  // cpy ctor: no implementation yet
  // requested by the Coding Convention
  Fatal("cpy ctor", "not implemented") ;
}


//____________________________________________________________________________
 AliPHOSTrackSegmentMakerv1::~AliPHOSTrackSegmentMakerv1()
{ 
  // dtor
  // fDefaultInit = kTRUE if TrackSegmentMaker created by default ctor (to get just the parameters)
  if (!fDefaultInit)  
    delete fLinkUpArray ;
}


//____________________________________________________________________________
const TString AliPHOSTrackSegmentMakerv1::BranchName() const 
{  
 
  return GetName() ;
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::FillOneModule()
{
  // Finds first and last indexes between which 
  // clusters from one PHOS module are

  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  
  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints() ; 
 
  //First EMC clusters
  Int_t totalEmc = emcRecPoints->GetEntriesFast() ;
  for(fEmcFirst = fEmcLast; (fEmcLast < totalEmc) &&  
	((dynamic_cast<AliPHOSRecPoint *>(emcRecPoints->At(fEmcLast)))->GetPHOSMod() == fModule ); 
      fEmcLast ++)  ;
  
  //Now CPV clusters
  Int_t totalCpv = cpvRecPoints->GetEntriesFast() ;

    for(fCpvFirst = fCpvLast; (fCpvLast < totalCpv) && 
         ((dynamic_cast<AliPHOSRecPoint *>(cpvRecPoints->At(fCpvLast)))->GetPHOSMod() == fModule ); 
       fCpvLast ++) ;
      
}

//____________________________________________________________________________
Float_t  AliPHOSTrackSegmentMakerv1::GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * emcClu,AliPHOSCpvRecPoint * cpvClu, Int_t &trackindex) const
{
  // Calculates the distance between the EMC RecPoint and the CPV RecPoint
  // Clusters are sorted in "rows" and "columns" of width 1 cm

  //Float_t delta = 1 ;  // Width of the rows in sorting of RecPoints (in cm)
                       // if you change this value, change it as well in xxxRecPoint::Compare()
  Float_t distance2Cpv   = fRcpv ;
  Float_t distance2Track = fRtpc ; 

  trackindex = -1 ; // closest track within fRCpv 

  TVector3 vecEmc ;   // Local position of EMC recpoint
  TVector3 vecCpv ;   // Local position of CPV recpoint propagated to EMC
  TVector3 vecDist ;  // Distance between local positions of two points
  
  emcClu->GetLocalPosition(vecEmc) ;
  cpvClu->GetLocalPosition(vecCpv) ;

  //toofar = kTRUE ;
  if(emcClu->GetPHOSMod() == cpvClu->GetPHOSMod()){ 

    // Find EMC-CPV distance
    distance2Cpv = (vecCpv - vecEmc).Mag() ;
    
    if (fESD != 0x0) {
      AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
      const AliPHOSGeometry * geom = gime->PHOSGeometry() ; 

      Double_t rPHOS = geom->GetIPtoCrystalSurface();

      //PH Acceptance boundaries for each PHOS module
      Int_t nModules = geom->GetNModules();
      Double_t * thmin = new Double_t[nModules];// theta min
      Double_t * thmax = new Double_t[nModules];// theta max
      Double_t * phmin = new Double_t[nModules];// phi min
      Double_t * phmax = new Double_t[nModules];// phi max
      
      for (Int_t imod=0; imod<nModules; imod++) {
	// Modules are numbered from 1 to 5 in AliPHOSGeometry
	geom->EmcModuleCoverage(imod+1,
				thmin[imod],thmax[imod],
				phmin[imod],phmax[imod]);
      }

      // Extrapolate the global track direction if any to CPV and find the closest track
      Int_t nTracks = fESD->GetNumberOfTracks();
      Int_t iClosestTrack = -1;
      Double_t minDistance = 1e6;
      Double_t pxyz[3], xyz[3];
      TVector3 inPHOS; //PH Used to calculate theta and phi

      //PH Loop on tracks
      AliESDtrack *track;
      for (Int_t iTrack=0; iTrack<nTracks; iTrack++) {
	track = fESD->GetTrack(iTrack);
	if (!track->GetXYZAt(rPHOS, fESD->GetMagneticField(), xyz))
           continue; //track coord on the cylinder of PHOS radius
	if ((TMath::Abs(xyz[0])+TMath::Abs(xyz[1])+TMath::Abs(xyz[2]))<=0)
	   continue;
	//PH Here one has to cut out the tracks which are not inside the PHOS
	//PH acceptance
	inPHOS.SetXYZ(xyz[0],xyz[1],xyz[2]);
	Double_t inPhi = inPHOS.Phi();
	Double_t inTheta = inPHOS.Theta();

	Bool_t skip = kTRUE;
	for (Int_t imod=0; imod<nModules; imod++) {
	  //PH Loop on modules to check if the track enters in the acceptance 
	  if (thmin[imod] < inTheta && thmax[imod] > inTheta && 
	      phmin[imod] < inPhi   && phmax[imod] > inPhi) {
	    skip = kFALSE;
	    break;
	  }
	}
	if (skip) continue; //PH Skip, if not in the PHOS acceptance

	if (!track->GetPxPyPzAt(rPHOS, fESD->GetMagneticField(), pxyz))
           continue; // track momentum ibid.
	PropagateToPlane(vecDist,xyz,pxyz,"CPV",cpvClu->GetPHOSMod());
	// 	Info("GetDistanceInPHOSPlane","Track %d propagation to CPV = (%f,%f,%f)",
 	//     iTrack,vecDist.X(),vecDist.Y(),vecDist.Z());
	vecDist -= vecCpv;
	distance2Track = TMath::Sqrt(vecDist.X()*vecDist.X() + vecDist.Z()*vecDist.Z());
	// Find the closest track to the EMC recpoint
	if (distance2Track < minDistance) {
	  minDistance = distance2Track;
	  iClosestTrack = iTrack;
	}
      }

      delete [] thmin;
      delete [] thmax;
      delete [] phmin;
      delete [] phmax;

      if (iClosestTrack != -1) {
	track = fESD->GetTrack(iClosestTrack);
	if (track->GetPxPyPzAt(rPHOS, fESD->GetMagneticField(), pxyz)) { // track momentum ibid.
	TVector3 vecCpvGlobal; // Global position of the CPV recpoint
	geom->GetGlobal((AliRecPoint*)cpvClu,vecCpvGlobal);
	for (Int_t ixyz=0; ixyz<3; ixyz++)
	  xyz[ixyz] = vecCpvGlobal[ixyz];
	PropagateToPlane(vecDist,xyz,pxyz,"EMC",cpvClu->GetPHOSMod());
// 	Info("GetDistanceInPHOSPlane","Track %d propagation to EMC = (%f,%f,%f)",
// 	     iClosestTrack,vecDist.X(),vecDist.Y(),vecDist.Z());
	vecDist -= vecEmc;
	distance2Track = TMath::Sqrt(vecDist.X()*vecDist.X() + vecDist.Z()*vecDist.Z());
	}
      }
//     } else {
//       // If no ESD exists, than simply find EMC-CPV distance
//       distance = (vecCpv - vecEmc).Mag() ;
    
      //if(distance2Track < fRcpv + 2*delta )
      if(distance2Track < fRtpc )
	trackindex = iClosestTrack ; 
      //      toofar = kFALSE ;
    }
    //     Info("GetDistanceInPHOSPlane","cpv-emc distance is %f cm",
    // 	 distance);
  }
  
  return distance2Cpv ;
}

//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::PropagateToPlane(TVector3& globalIntersection,
						  Double_t *x,
						  Double_t *p,
						  const char *det,
						  Int_t moduleNumber) const
{
  // Propagate a straight-line track from the origin point x
  // along the direction p to the CPV or EMC module moduleNumber
  // Returns a local position of such a propagation

  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  const AliPHOSGeometry * geom = gime->PHOSGeometry() ; 
  TVector3 moduleCenter;
  geom->GetModuleCenter(moduleCenter,det,moduleNumber);
  TVector3 vertex; vertex.SetXYZ(x[0],x[1],x[2]);
  TVector3 direction; direction.SetXYZ(p[0],p[1],p[2]);

//   Info("PropagateToCPV","Center of the %s module %d is (%f,%f,%f)",
//        det,moduleNumber,moduleCenter[0],moduleCenter[1],moduleCenter[2]);

  Double_t time = (moduleCenter.Mag2() - vertex.Dot(moduleCenter)) /
    (direction.Dot(moduleCenter));
  vertex += direction*time;
  geom->Global2Local(globalIntersection,vertex,moduleNumber);
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::Init()
{
  // Make all memory allocations that are not possible in default constructor
  
  AliPHOSGetter* gime = AliPHOSGetter::Instance();
  if(!gime)
    gime = AliPHOSGetter::Instance(GetTitle(), fEventFolderName.Data());
  
  fLinkUpArray  = new TClonesArray("AliPHOSLink", 1000); 
  if ( !gime->TrackSegmentMaker() ) {
    gime->PostTrackSegmentMaker(this);
  }
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::InitParameters()
{
  //Initializes parameters
  fRcpv      = 10. ;
  fRtpc      = 4. ;
  fEmcFirst  = 0 ;    
  fEmcLast   = 0 ;   
  fCpvFirst  = 0 ;   
  fCpvLast   = 0 ;   
  fLinkUpArray = 0 ;
  fWrite                   = kTRUE ;
  fTrackSegmentsInRun       = 0 ; 
  SetEventRange(0,-1) ;
}


//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakeLinks()const
{ 
  // Finds distances (links) between all EMC and CPV clusters, 
  // which are not further apart from each other than fRcpv 
  // and sort them in accordance with this distance
  
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints() ; 

  fLinkUpArray->Clear() ;    

  AliPHOSCpvRecPoint * cpv ;
  AliPHOSEmcRecPoint * emcclu ;

  Int_t iLinkUp  = 0 ;
  
  Int_t iEmcRP;
  for(iEmcRP = fEmcFirst; iEmcRP < fEmcLast; iEmcRP++ ) {
    emcclu = dynamic_cast<AliPHOSEmcRecPoint *>(emcRecPoints->At(iEmcRP)) ;

    //Bool_t toofar ;        
    Int_t iCpv = 0 ;    
    for(iCpv = fCpvFirst; iCpv < fCpvLast;iCpv++ ) { 
      
      cpv = dynamic_cast<AliPHOSCpvRecPoint *>(cpvRecPoints->At(iCpv)) ;
      Int_t track = -1 ; 
      Float_t r = GetDistanceInPHOSPlane(emcclu, cpv, track) ;     
      //      if(toofar)
      //	continue ;	 
      if(r < fRcpv) { 
        new ((*fLinkUpArray)[iLinkUp++])  AliPHOSLink(r, iEmcRP, iCpv, track) ;
      }      
    }
  } 
  
  fLinkUpArray->Sort() ;  //first links with smallest distances
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakePairs()
{ 
  // Using the previously made list of "links", we found the smallest link - i.e. 
  // link with the least distance between EMC and CPV and pointing to still 
  // unassigned RecParticles. We assign these RecPoints to TrackSegment and 
  // remove them from the list of "unassigned". 

  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 

  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints() ; 
  TClonesArray * trackSegments = gime->TrackSegments();
    
  //Make arrays to mark clusters already chosen
  Int_t * emcExist = 0;
  if(fEmcLast > fEmcFirst)
    emcExist = new Int_t[fEmcLast-fEmcFirst] ;
  
  Int_t index;
  for(index = 0; index <fEmcLast-fEmcFirst; index ++)
    emcExist[index] = 1 ;
  
  Bool_t * cpvExist = 0;
  if(fCpvLast > fCpvFirst)
    cpvExist = new Bool_t[fCpvLast-fCpvFirst] ;
  for(index = 0; index <fCpvLast-fCpvFirst; index ++)
    cpvExist[index] = kTRUE ;
  
  
  // Finds the smallest links and makes pairs of CPV and EMC clusters with smallest distance 
  TIter nextUp(fLinkUpArray) ;
  
  AliPHOSLink * linkUp ;
  
  AliPHOSCpvRecPoint * nullpointer = 0 ;
  
  while ( (linkUp =  static_cast<AliPHOSLink *>(nextUp()) ) ){  

    if(emcExist[linkUp->GetEmc()-fEmcFirst] != -1){

      if(cpvExist[linkUp->GetCpv()-fCpvFirst]){ //CPV still exist
	 new ((* trackSegments)[fNTrackSegments]) 
	   AliPHOSTrackSegment(dynamic_cast<AliPHOSEmcRecPoint *>(emcRecPoints->At(linkUp->GetEmc())) , 
			       dynamic_cast<AliPHOSCpvRecPoint *>(cpvRecPoints->At(linkUp->GetCpv())) , 
			       linkUp->GetTrack()) ;
	 
       (dynamic_cast<AliPHOSTrackSegment *>(trackSegments->At(fNTrackSegments)))->SetIndexInList(fNTrackSegments);
       fNTrackSegments++ ;
       emcExist[linkUp->GetEmc()-fEmcFirst] = -1 ; //Mark emc  that Cpv was found 
       //mark CPV recpoint as already used 
       cpvExist[linkUp->GetCpv()-fCpvFirst] = kFALSE ;
      } //if CpvUp still exist
    } 
  }        

  //look through emc recPoints left without CPV
  if(emcExist){ //if there is emc rec point
    Int_t iEmcRP ;
    for(iEmcRP = 0; iEmcRP < fEmcLast-fEmcFirst  ; iEmcRP++ ){
      if(emcExist[iEmcRP] > 0 ){
       new ((*trackSegments)[fNTrackSegments])  
         AliPHOSTrackSegment(dynamic_cast<AliPHOSEmcRecPoint *>(emcRecPoints->At(iEmcRP+fEmcFirst)), 
                           nullpointer) ;
       (dynamic_cast<AliPHOSTrackSegment *>(trackSegments->At(fNTrackSegments)))->SetIndexInList(fNTrackSegments);
       fNTrackSegments++;    
      } 
    }
  }
  delete [] emcExist ; 
  delete [] cpvExist ; 
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::Exec(Option_t *option)
{
  // Steering method to perform track segment construction for events
  // in the range from fFirstEvent to fLastEvent.
  // This range is optionally set by SetEventRange().
  // if fLastEvent=-1 (by default), then process events until the end.
  
  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSTSMaker");
 
  if(strstr(option,"print")) {
    Print() ; 
    return ; 
  }
  
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ;  
 
  const AliPHOSGeometry * geom = gime->PHOSGeometry() ; 

  if (fLastEvent == -1) 
    fLastEvent = gime->MaxEvent() - 1 ;
  else 
    fLastEvent = TMath::Min(fFirstEvent,gime->MaxEvent());
  Int_t nEvents   = fLastEvent - fFirstEvent + 1;

  Int_t ievent ; 
  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
    gime->Event(ievent,"R") ;
   //Make some initializations 
    fNTrackSegments = 0 ;
    fEmcFirst = 0 ;    
    fEmcLast  = 0 ;   
    fCpvFirst = 0 ;   
    fCpvLast  = 0 ;   
    
    gime->TrackSegments()->Clear();

    //    if(!ReadRecPoints(ievent))   continue; //reads RecPoints for event ievent
    
    for(fModule = 1; fModule <= geom->GetNModules() ; fModule++ ) {
      FillOneModule() ; 
      MakeLinks() ;
      MakePairs() ;
    }

    WriteTrackSegments() ;

    if(strstr(option,"deb"))
      PrintTrackSegments(option);
    
    //increment the total number of track segments per run 
    fTrackSegmentsInRun += gime->TrackSegments()->GetEntriesFast() ; 
  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSTSMaker");
    Info("Exec", "took %f seconds for making TS %f seconds per event", 
          gBenchmark->GetCpuTime("PHOSTSMaker"), 
          gBenchmark->GetCpuTime("PHOSTSMaker")/nEvents) ;
   }
  if(fWrite) //do not unload in "on flight" mode
    Unload();
}

//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::Unload() 
{
  // Unloads the task from the folder
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ;  
  gime->PhosLoader()->UnloadRecPoints() ;
  gime->PhosLoader()->UnloadTracks() ;
}

//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::Print(const Option_t *)const
{
  //  Print TrackSegmentMaker parameters

  TString message("") ;
  if( strcmp(GetName(), "") != 0 ) {
    message = "\n======== AliPHOSTrackSegmentMakerv1 ========\n" ; 
    message += "Making Track segments\n" ;
    message += "with parameters:\n" ; 
    message += "     Maximal EMC - CPV distance (cm) %f\n" ;
    message += "============================================\n" ;
    Info("Print", message.Data(),fRcpv) ;
  }
  else
    Info("Print", "AliPHOSTrackSegmentMakerv1 not initialized ") ;
}

//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::WriteTrackSegments()
{
  // Writes found TrackSegments to TreeR. Creates branches 
  // "PHOSTS" and "AliPHOSTrackSegmentMaker" with the same title.
  // In the former branch found TrackSegments are stored, while 
  // in the latter all parameters, with which TS were made. 
  // ROOT does not allow overwriting existing branches, therefore
  // first we check, if branches with the same title already exist.
  // If yes - exits without writing.

  AliPHOSGetter *gime = AliPHOSGetter::Instance() ; 

  TClonesArray * trackSegments = gime->TrackSegments() ; 
  trackSegments->Expand(trackSegments->GetEntriesFast()) ;

  if(fWrite){ //We write TreeT
    TTree * treeT = gime->TreeT();
    
    //First TS
    Int_t bufferSize = 32000 ; 
    TBranch * tsBranch = treeT->Branch("PHOSTS",&trackSegments,bufferSize);
    tsBranch->Fill() ;  
    
    gime->WriteTracks("OVERWRITE");
    gime->WriteTrackSegmentMaker("OVERWRITE");
  }
}


//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::PrintTrackSegments(Option_t * option)
{
  // option deb - prints # of found TrackSegments
  // option deb all - prints as well indexed of found RecParticles assigned to the TS

  TClonesArray * trackSegments = AliPHOSGetter::Instance()->TrackSegments() ; 

  Info("PrintTrackSegments", "Results from TrackSegmentMaker:") ; 
  printf("nevent: %d\n", gAlice->GetEvNumber()) ; 
  printf("        Found %d TrackSegments\n", trackSegments->GetEntriesFast() ); 
  
  if(strstr(option,"all")) {  // printing found TS
    printf("TrackSegment #  EMC RP#  CPV RP#\n") ; 
    Int_t index;
    for (index = 0 ; index <trackSegments->GetEntriesFast() ; index++) {
      AliPHOSTrackSegment * ts = (AliPHOSTrackSegment * )trackSegments->At(index) ; 
      printf("   %d           %d        %d \n", ts->GetIndexInList(), ts->GetEmcIndex(), ts->GetCpvIndex() ) ; 
    }	
  }
}
