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
//_________________________________________________________________________
// Implementation version 1 of algorithm class to construct PHOS track segments
// Track segment for PHOS is list of 
//        EMC RecPoint + (possibly) CPV RecPoint + (possibly) PPSD RecPoint
// To find TrackSegments we do the following: 
//  for each EMC RecPoints we look at
//   CPV/PPSD RecPoints in the radious fRcpv. 
//  If there is such a CPV RecPoint, 
//   we make "Link" it is just indexes of EMC and CPV/PPSD RecPoint and distance
//   between them in the PHOS plane. 
//  Then we sort "Links" and starting from the 
//   least "Link" pointing to the unassined EMC and CPV RecPoints assing them to 
//   new TrackSegment. 
// If there is no CPV/PPSD RecPoint we make TrackSegment 
// consisting from EMC alone. There is no TrackSegments without EMC RecPoint.
//// In principle this class should be called from AliPHOSReconstructor, but 
// one can use it as well in standalone mode.
// Use  case:
//  root [0] AliPHOSTrackSegmentMakerv1 * t = new AliPHOSTrackSegmentMaker("galice.root", "tracksegmentsname", "recpointsname")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//               // reads gAlice from header file "galice.root", uses recpoints stored in the branch names "recpointsname" (default = "Default")
//               // and saves recpoints in branch named "tracksegmentsname" (default = "recpointsname")                       
//  root [1] t->ExecuteTask()
//  root [2] t->SetMaxEmcPpsdDistance(5)
//  root [3] t->SetTrackSegmentsBranch("max distance 5 cm")
//  root [4] t->ExecuteTask("deb all time") 
//                 
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH) & Yves Schutz (SUBATECH) 
//

// --- ROOT system ---
#include "TTree.h"
#include "TBenchmark.h"

// --- Standard library ---

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
  AliPHOSTrackSegmentMakerv1::AliPHOSTrackSegmentMakerv1() : AliPHOSTrackSegmentMaker()
{
  // default ctor (to be used mainly by Streamer)

  InitParameters() ; 
  fDefaultInit = kTRUE ; 
}

//____________________________________________________________________________
 AliPHOSTrackSegmentMakerv1::AliPHOSTrackSegmentMakerv1(const TString alirunFileName, const TString eventFolderName)
   :AliPHOSTrackSegmentMaker(alirunFileName, eventFolderName)
{
  // ctor

  InitParameters() ; 
  Init() ;
  fDefaultInit = kFALSE ; 
  fESD = 0;
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
Float_t  AliPHOSTrackSegmentMakerv1::GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * emcClu,AliPHOSCpvRecPoint * cpvClu, Bool_t &toofar) const
{
  // Calculates the distance between the EMC RecPoint and the PPSD RecPoint
  // Clusters are sorted in "rows" and "columns" of width 1 cm

  Float_t delta = 1 ;  // Width of the rows in sorting of RecPoints (in cm)
                       // if you change this value, change it as well in xxxRecPoint::Compare()
  Float_t distance = fRcpv ;
 
  TVector3 vecEmc ;   // Local position of EMC recpoint
  TVector3 vecCpv ;   // Local position of CPV recpoint propagated to EMC
  TVector3 vecDist ;  // Distance between local positions of two points
  
  emcClu->GetLocalPosition(vecEmc) ;
  cpvClu->GetLocalPosition(vecCpv) ;

  toofar = kTRUE ;
  if(emcClu->GetPHOSMod() == cpvClu->GetPHOSMod()){ 

    if (fESD != 0x0) {
      // Extrapolate the global track direction if any to CPV and find the closest track
      Int_t nTracks = fESD->GetNumberOfTracks();
      Int_t iClosestTrack = -1;
      Double_t minDistance = 1e6;
      Double_t pxyz[3], xyz[3];
      AliESDtrack *track;
      for (Int_t iTrack=0; iTrack<nTracks; iTrack++) {
	track = fESD->GetTrack(iTrack);
	track->GetOuterXYZ(xyz);     // track coord on the cylinder of PHOS radius
	if ((TMath::Abs(xyz[0])+TMath::Abs(xyz[1])+TMath::Abs(xyz[2]))<=0)
	  continue;
	track->GetOuterPxPyPz(pxyz); // track momentum ibid.
	vecDist = PropagateToPlane(xyz,pxyz,"CPV",cpvClu->GetPHOSMod());
// 	Info("GetDistanceInPHOSPlane","Track %d propagation to CPV = (%f,%f,%f)",
// 	     iTrack,vecDist.X(),vecDist.Y(),vecDist.Z());
	vecDist -= vecCpv;
	distance = TMath::Sqrt(vecDist.X()*vecDist.X() + vecDist.Z()*vecDist.Z());
	// Find the closest track to the EMC recpoint
	if (distance < minDistance) {
	  minDistance = distance;
	  iClosestTrack = iTrack;
	}
      }

      if (iClosestTrack != -1) {
	track = fESD->GetTrack(iClosestTrack);
	track->GetOuterPxPyPz(pxyz); // track momentum ibid.
	TVector3 vecCpvGlobal; // Global position of the CPV recpoint
	AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
	const AliPHOSGeometry * geom = gime->PHOSGeometry() ; 
	geom->GetGlobal((AliRecPoint*)cpvClu,vecCpvGlobal);
	for (Int_t ixyz=0; ixyz<3; ixyz++)
	  xyz[ixyz] = vecCpvGlobal[ixyz];
	vecDist = PropagateToPlane(xyz,pxyz,"EMC",cpvClu->GetPHOSMod());
// 	Info("GetDistanceInPHOSPlane","Track %d propagation to EMC = (%f,%f,%f)",
// 	     iClosestTrack,vecDist.X(),vecDist.Y(),vecDist.Z());
	vecDist -= vecEmc;
	distance = TMath::Sqrt(vecDist.X()*vecDist.X() + vecDist.Z()*vecDist.Z());
      }
    } else {
      // If no ESD exists, than simply find EMC-CPV distance
      distance = (vecCpv - vecEmc).Mag() ;
    }
//     Info("GetDistanceInPHOSPlane","cpv-emc distance is %f cm",
// 	 distance);
    if(distance < fRcpv + 2*delta )
      toofar = kFALSE ;
  }
  
  return distance ;
}

//____________________________________________________________________________
TVector3  AliPHOSTrackSegmentMakerv1::PropagateToPlane(Double_t *x, Double_t *p,
						       char *det, Int_t moduleNumber) const
{
  // Propagate a straight-line track from the origin point x
  // along the direction p to the CPV or EMC module moduleNumber
  // Returns a local position of such a propagation

  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  const AliPHOSGeometry * geom = gime->PHOSGeometry() ; 
  TVector3 moduleCenter = geom->GetModuleCenter(det,moduleNumber);
  TVector3 vertex(x);
  TVector3 direction(p);

//   Info("PropagateToCPV","Center of the %s module %d is (%f,%f,%f)",
//        det,moduleNumber,moduleCenter[0],moduleCenter[1],moduleCenter[2]);

  Double_t time = (moduleCenter.Mag2() - vertex.Dot(moduleCenter)) /
    (direction.Dot(moduleCenter));
  TVector3 globalIntersection = vertex + direction*time;
  return geom->Global2Local(globalIntersection,moduleNumber);
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::Init()
{
  // Make all memory allocations that are not possible in default constructor
  
  AliPHOSGetter* gime = AliPHOSGetter::Instance(GetTitle(), fEventFolderName.Data());
  
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
  fEmcFirst  = 0 ;    
  fEmcLast   = 0 ;   
  fCpvFirst  = 0 ;   
  fCpvLast   = 0 ;   
  fLinkUpArray = 0 ;
  fTrackSegmentsInRun       = 0 ; 
  SetEventRange(0,-1) ;
}


//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakeLinks()const
{ 
  // Finds distances (links) between all EMC and PPSD clusters, 
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

    Bool_t toofar ;        
    Int_t iCpv = 0 ;    
    for(iCpv = fCpvFirst; iCpv < fCpvLast;iCpv++ ) { 
      
      cpv = dynamic_cast<AliPHOSCpvRecPoint *>(cpvRecPoints->At(iCpv)) ;
      Float_t r = GetDistanceInPHOSPlane(emcclu, cpv, toofar) ;
      
      if(toofar)
        break ;	 
      if(r < fRcpv) { 
        new ((*fLinkUpArray)[iLinkUp++])  AliPHOSLink(r, iEmcRP, iCpv) ;
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

    if(emcExist[linkUp->GetEmc()-fEmcFirst] != -1){ //without ppsd Up yet 

      if(cpvExist[linkUp->GetPpsd()-fCpvFirst]){ //CPV still exist
       
       new ((* trackSegments)[fNTrackSegments]) 
         AliPHOSTrackSegment(dynamic_cast<AliPHOSEmcRecPoint *>(emcRecPoints->At(linkUp->GetEmc())) , 
			     dynamic_cast<AliPHOSCpvRecPoint *>(cpvRecPoints->At(linkUp->GetPpsd()))) ;
       (dynamic_cast<AliPHOSTrackSegment *>(trackSegments->At(fNTrackSegments)))->SetIndexInList(fNTrackSegments);
       fNTrackSegments++ ;
       
       emcExist[linkUp->GetEmc()-fEmcFirst] = -1 ; //Mark emc  that Cpv was found 
       //mark CPV recpoint as already used 
       cpvExist[linkUp->GetPpsd()-fCpvFirst] = kFALSE ;
      } //if ppsdUp still exist
    } 
  }        

  //look through emc recPoints left without CPV/PPSD
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
  
  AliPHOSGetter * gime = AliPHOSGetter::Instance(GetTitle()) ;  
 
  const AliPHOSGeometry * geom = gime->PHOSGeometry() ; 

  if (fLastEvent == -1) fLastEvent = gime->MaxEvent() - 1 ;
  else fLastEvent = TMath::Min(fFirstEvent,gime->MaxEvent());
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
void AliPHOSTrackSegmentMakerv1::Print()const
{
  //  Print TrackSegmentMaker parameters

  TString message("") ;
  if( strcmp(GetName(), "") != 0 ) {
    message = "\n======== AliPHOSTrackSegmentMakerv1 ========\n" ; 
    message += "Making Track segments\n" ;
    message += "with parameters:\n" ; 
    message += "     Maximal EMC - CPV (PPSD) distance (cm) %f\n" ;
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

  TTree * treeT = gime->TreeT();
 
  //First TS
  Int_t bufferSize = 32000 ; 
  TBranch * tsBranch = treeT->Branch("PHOSTS",&trackSegments,bufferSize);
  tsBranch->Fill() ;  

  gime->WriteTracks("OVERWRITE");
  gime->WriteTrackSegmentMaker("OVERWRITE");
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
