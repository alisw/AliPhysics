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
 * Revision 1.93  2007/10/10 09:05:10  schutz
 * Changing name QualAss to QA
 *
 * Revision 1.92  2007/08/28 12:55:08  policheh
 * Loaders removed from the reconstruction code (C.Cheshkov)
 *
 * Revision 1.91  2007/08/07 14:12:03  kharlov
 * Quality assurance added (Yves Schutz)
 *
 * Revision 1.90  2007/07/11 13:43:30  hristov
 * New class AliESDEvent, backward compatibility with the old AliESD (Christian)
 *
 * Revision 1.89  2007/07/03 08:13:04  kharlov
 * Bug fix in CPV local coordinates
 *
 * Revision 1.88  2007/06/27 09:11:07  kharlov
 * Bug fix for CPV-EMC distance
 *
 * Revision 1.87  2007/05/04 14:49:29  policheh
 * AliPHOSRecPoint inheritance from AliCluster
 *
 * Revision 1.86  2007/04/02 15:00:16  cvetan
 * No more calls to gAlice in the reconstruction
 *
 * Revision 1.85  2007/03/28 19:18:15  kharlov
 * RecPoints recalculation in TSM removed
 *
 * Revision 1.84  2007/03/07 07:01:21  hristov
 * Fixing copy/paste erro. Additional protections
 *
 * Revision 1.83  2007/03/06 21:07:37  kharlov
 * DP: xz CPV-EMC distance filled to TS
 *
 * Revision 1.82  2007/03/06 06:54:48  kharlov
 * DP:Calculation of cluster properties dep. on vertex added
 *
 * Revision 1.81  2007/02/05 10:02:40  kharlov
 * Module numbering is corrected
 *
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
#include "TVector3.h"
#include "TTree.h"
#include "TBenchmark.h"

// --- Standard library ---
#include "Riostream.h"
// --- AliRoot header files ---
#include "AliPHOSGeometry.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSLink.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSCpvRecPoint.h"
#include "AliLog.h"
#include "AliMagF.h"
#include "AliMagF.h"
#include "AliTracker.h"
#include "AliGeomManager.h"
#include "AliCluster.h"
#include "AliKalmanTrack.h"
#include "AliGlobalQADataMaker.h"
#include "AliVParticle.h"


ClassImp( AliPHOSTrackSegmentMakerv1) 


//____________________________________________________________________________
AliPHOSTrackSegmentMakerv1::AliPHOSTrackSegmentMakerv1() :
  AliPHOSTrackSegmentMaker(),
  fDefaultInit(kTRUE),
  fWrite(kFALSE),
  fNTrackSegments(0),
  fRcpv(0.f),
  fRtpc(0.f),
  fVtx(0.f), 
  fLinkUpArray(0),
  fEmcFirst(0),
  fEmcLast(0),
  fCpvFirst(0),
  fCpvLast(0),
  fModule(0),
  fTrackSegments(NULL)
{
  // default ctor (to be used mainly by Streamer)
  InitParameters() ; 
}

//____________________________________________________________________________
AliPHOSTrackSegmentMakerv1::AliPHOSTrackSegmentMakerv1(AliPHOSGeometry *geom) :
  AliPHOSTrackSegmentMaker(geom),
  fDefaultInit(kFALSE),
  fWrite(kFALSE),
  fNTrackSegments(0),
  fRcpv(0.f),
  fRtpc(0.f),
  fVtx(0.f), 
  fLinkUpArray(0),
  fEmcFirst(0),
  fEmcLast(0),
  fCpvFirst(0),
  fCpvLast(0),
  fModule(0),
  fTrackSegments(NULL)
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
  fVtx(0.f), 
  fLinkUpArray(0),
  fEmcFirst(0),
  fEmcLast(0),
  fCpvFirst(0),
  fCpvLast(0),
  fModule(0),
  fTrackSegments(NULL)
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
  if (fTrackSegments) {
    fTrackSegments->Delete();
    delete fTrackSegments;
  }
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::FillOneModule()
{
  // Finds first and last indexes between which 
  // clusters from one PHOS module are

  //First EMC clusters
  Int_t totalEmc = fEMCRecPoints->GetEntriesFast() ;
  for(fEmcFirst = fEmcLast; (fEmcLast < totalEmc) &&  
	((static_cast<AliPHOSRecPoint *>(fEMCRecPoints->At(fEmcLast)))->GetPHOSMod() == fModule ); 
      fEmcLast ++)  ;
  
  //Now CPV clusters
  Int_t totalCpv = fCPVRecPoints->GetEntriesFast() ;

    for(fCpvFirst = fCpvLast; (fCpvLast < totalCpv) && 
         ((static_cast<AliPHOSRecPoint *>(fCPVRecPoints->At(fCpvLast)))->GetPHOSMod() == fModule ); 
       fCpvLast ++) ;
      
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * emcClu,
                                                         AliPHOSCpvRecPoint * cpvClu, 
                                                         Int_t &trackindex, 
                                                         Float_t &dx, Float_t &dz) const
{
  // Calculates the distance between the EMC RecPoint and the CPV RecPoint
  // If no CPV, calculates the distance between the EMC RecPoint and the track
  // prolongation to the PHOS module plane.
  // Clusters are sorted in "rows" and "columns" of width 1 cm

//  Float_t delta = 1 ;  // Width of the rows in sorting of RecPoints (in cm)
//                       // if you change this value, change it as well in xxxRecPoint::Compare()

  trackindex = -1;
  dx         = 999.;
  dz         = 999.;

  if(!cpvClu) {
    
    if(!emcClu) {
      return;
    }

    // *** Start the matching
    Int_t nt=fESD->GetNumberOfTracks();
    Int_t iPHOSMod = emcClu->GetPHOSMod()  ;
    //Calculate actual distance to PHOS module
    TVector3 globaPos ;
    fGeom->Local2Global(iPHOSMod, 0.,0., globaPos) ;
    const Double_t rPHOS = globaPos.Pt() ; //Distance to center of  PHOS module
    const Double_t kYmax = 72.+10. ; //Size of the module (with some reserve) in phi direction
    const Double_t kZmax = 64.+10. ; //Size of the module (with some reserve) in z direction
    const Double_t kAlpha0=330./180.*TMath::Pi() ; //First PHOS module angular direction
    const Double_t kAlpha= 20./180.*TMath::Pi() ; //PHOS module angular size
    Double_t minDistance = 1.e6;

    TVector3 vecEmc ;   // Local position of EMC recpoint
    emcClu->GetLocalPosition(vecEmc) ;

    Double_t gposTrack[3] ; 
    Double_t bz = AliTracker::GetBz() ; //B-Field for approximate matching
    Double_t b[3]; 
    for (Int_t i=0; i<nt; i++) {
      AliESDtrack *esdTrack=fESD->GetTrack(i);

      // Skip the tracks having "wrong" status (has to be checked/tuned)
      ULong_t status = esdTrack->GetStatus();
      if ((status & AliESDtrack::kTPCout)   == 0) continue;
//     if ((status & AliESDtrack::kTRDout)   == 0) continue;
//     if ((status & AliESDtrack::kTRDrefit) == 1) continue;

      //Continue extrapolation from TPC outer surface
      const AliExternalTrackParam *outerParam=esdTrack->GetOuterParam();
      if (!outerParam) continue;
      AliExternalTrackParam t(*outerParam);

      t.GetBxByBz(b) ;
      //Direction to the current PHOS module
      Double_t phiMod=kAlpha0-kAlpha*iPHOSMod ;
      if(!t.Rotate(phiMod))
        continue ;
 
      Double_t y;                       // Some tracks do not reach the PHOS
      if (!t.GetYAt(rPHOS,bz,y)) continue; //    because of the bending

      Double_t z; 
      if(!t.GetZAt(rPHOS,bz,z))
        continue ;
      if (TMath::Abs(z) > kZmax) 
        continue; // Some tracks miss the PHOS in Z
      if(TMath::Abs(y) < kYmax){
        t.PropagateToBxByBz(rPHOS,b);        // Propagate to the matching module
      //t.CorrectForMaterial(...); // Correct for the TOF material, if needed
        t.GetXYZ(gposTrack) ;
        TVector3 globalPositionTr(gposTrack) ;
        TVector3 localPositionTr ;
        fGeom->Global2Local(localPositionTr,globalPositionTr,iPHOSMod) ;
        Double_t ddx = vecEmc.X()-localPositionTr.X();
        Double_t ddz = vecEmc.Z()-localPositionTr.Z();
        Double_t d2 = ddx*ddx + ddz*ddz;
        if(d2 < minDistance) {
	  dx = ddx ;
	  dz = ddz ;
	  trackindex=i;
	  minDistance=d2 ;
        }
      }
    } //Scanned all tracks
    return ;
  }


  TVector3 emcGlobal;
  fGeom->GetGlobalPHOS((AliPHOSRecPoint*)emcClu,emcGlobal);
    
  // Radius from IP to current point
  Double_t rEMC = TMath::Abs(emcGlobal.Pt());
    
  // Extrapolate the global track direction to EMC
  // and find the closest track
    
  Int_t nTracks = fESD->GetNumberOfTracks();
    
  AliESDtrack *track;
  Double_t xyz[] = {-1,-1,-1};
  Double_t pxyz[3];
  Double_t zEMC,xEMC;
  Int_t module;
  TVector3 vecP;
  TVector3 locClu;

  Float_t minDistance = 1.e6;
  Float_t dr;

  for (Int_t iTrack=0; iTrack<nTracks; iTrack++) {
    track = fESD->GetTrack(iTrack);
    if (!track->GetXYZAt(rEMC, fESD->GetMagneticField(), xyz)) continue;
    
    AliDebug(1,Form("Event %d, iTrack: %d, (%.3f,%.3f,%.3f)",
		    fESD->GetEventNumberInFile(),iTrack,xyz[0],xyz[1],xyz[2]));
      
    if (track->GetPxPyPzAt(rEMC,fESD->GetMagneticField(),pxyz)) { 

      vecP.SetXYZ(pxyz[0],pxyz[1],pxyz[2]);
      fGeom->ImpactOnEmc(xyz,vecP.Theta(),vecP.Phi(),module,zEMC,xEMC) ;

      if(!module) continue;
      AliDebug(1,Form("\t\tTrack hit PHOS! Module: %d, (x,z)=(%.3f,%.3f)",module,xEMC,zEMC));
	
      if(emcClu->GetPHOSMod() != module) continue;
	
      // match track to EMC cluster
      emcClu->GetLocalPosition(locClu);
	
      Float_t delta_x = xEMC - locClu.X();
      Float_t delta_z = zEMC - locClu.Z();
      dr = TMath::Sqrt(delta_x*delta_x + delta_z*delta_z);
      AliDebug(1,Form("\tMatch iTrack=%d: (dx,dz)=(%.3f,%.3f)",iTrack,delta_x,delta_z));
	
      if(dr<minDistance) {
	trackindex = iTrack;
	minDistance = dr;
	dx = delta_x;
	dz = delta_z;
      }
    }
      
  }
    
  if(trackindex>=0) {
    AliDebug(1,Form("\t\tBest match for (xClu,zClu,eClu)=(%.3f,%.3f,%.3f): iTrack=%d, dR=%.3f",
		    locClu.X(),locClu.Z(),emcClu->GetEnergy(),
		    trackindex,TMath::Sqrt(dx*dx+dz*dz)));	  
    return;
  }
  
  Float_t distance2Track = fRtpc ; 
  
  trackindex = -1 ; // closest track within fRCpv 
  
  TVector3 vecEmc ;   // Local position of EMC recpoint
  TVector3 vecPloc ;     // Momentum direction at CPV plain
  
  //toofar = kTRUE ;
  if(emcClu->GetPHOSMod() != cpvClu->GetPHOSMod()){
    dx=999. ;
    dz=999. ;
    return ;
  }
  
  emcClu->GetLocalPosition(vecEmc) ;
  
  Double_t xCPV=0,zCPV=0 ; //EMC-projected coordinates of CPV cluster 
  TVector3 cpvGlobal; // Global position of the CPV recpoint
  fGeom->GetGlobalPHOS((AliPHOSRecPoint*)cpvClu,cpvGlobal);
  Double_t vtxCPV[3]={cpvGlobal.X(),cpvGlobal.Y(),cpvGlobal.Z()} ;
  Int_t dummyMod ;

  if (fESD == 0x0) {
    //if no track information available, assume straight line from IP to emcal
    fGeom->ImpactOnEmc(vtxCPV,cpvGlobal.Theta(),cpvGlobal.Phi(),dummyMod,zCPV,xCPV) ;
    dx=xCPV - vecEmc.X() ;
    dz=zCPV - vecEmc.Z() ;
    return ;
  } 
  
  //if there is ESD try to correct distance using TPC information on particle direct in CPV 
  if (fESD != 0x0) {

    Double_t rCPV = cpvGlobal.Pt() ;// Radius from IP to current point 

    // Extrapolate the global track direction if any to CPV and find the closest track
    Int_t iClosestTrack = -1;
    TVector3 inPHOS ; 

    for (Int_t iTrack=0; iTrack<nTracks; iTrack++) {
      track = fESD->GetTrack(iTrack);
      if (!track->GetXYZAt(rCPV, fESD->GetMagneticField(), xyz))
	continue; //track coord on the cylinder of PHOS radius
      if ((TMath::Abs(xyz[0])+TMath::Abs(xyz[1])+TMath::Abs(xyz[2]))<=0)
	continue;
      //Check if this track hits PHOS
      inPHOS.SetXYZ(xyz[0],xyz[1],xyz[2]);
      distance2Track = inPHOS.Angle(cpvGlobal) ;
      // Find the closest track to the CPV recpoint
      if (distance2Track < minDistance) {
	minDistance = distance2Track;
	iClosestTrack = iTrack;
      }
    }

    if (iClosestTrack != -1) {
      track = fESD->GetTrack(iClosestTrack);
      if (track->GetPxPyPzAt(rCPV, fESD->GetMagneticField(), pxyz)) { // track momentum ibid.
        vecP.SetXYZ(pxyz[0],pxyz[1],pxyz[2]);
	fGeom->ImpactOnEmc(vtxCPV,vecP.Theta(),vecP.Phi(),dummyMod,zCPV,xCPV) ;
      }
    }
    
    if(minDistance < fRtpc ){
      trackindex = iClosestTrack ; 
    }
  }
  if(trackindex!=-1){
    // If the closest global track is found, calculate EMC-CPV distance from it
    dx=xCPV - vecEmc.X() ;
    dz=zCPV - vecEmc.Z() ;
  }
  else{
    // If no global track was found, just take the nearest CPV point
    fGeom->ImpactOnEmc(vtxCPV,cpvGlobal.Theta(),cpvGlobal.Phi(),dummyMod,zCPV,xCPV) ;
    dx=xCPV - vecEmc.X() ;
    dz=zCPV - vecEmc.Z() ;
  }
  return ;
}
//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::Init()
{
  // Make all memory allocations that are not possible in default constructor
  
  fLinkUpArray  = new TClonesArray("AliPHOSLink", 1000); 
  fTrackSegments = new TClonesArray("AliPHOSTrackSegment",100);
  fTrackSegments->SetName("TRACKS");
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
}


//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakeLinks()const
{ 
  // Finds distances (links) between all EMC and CPV clusters, 
  // which are not further apart from each other than fRcpv 
  // and sort them in accordance with this distance
  
  fLinkUpArray->Clear() ;    

  AliPHOSCpvRecPoint * cpv ;
  AliPHOSEmcRecPoint * emcclu ;

  Int_t iLinkUp  = 0 ;
  
  Int_t iEmcRP;
  for(iEmcRP = fEmcFirst; iEmcRP < fEmcLast; iEmcRP++ ) {
    emcclu = dynamic_cast<AliPHOSEmcRecPoint *>(fEMCRecPoints->At(iEmcRP)) ;

    //Bool_t toofar ;        
    Int_t iCpv = 0 ;    
    for(iCpv = fCpvFirst; iCpv < fCpvLast;iCpv++ ) { 
      
      cpv = dynamic_cast<AliPHOSCpvRecPoint *>(fCPVRecPoints->At(iCpv)) ;
      Int_t track = -1 ; 
      Float_t dx,dz ;
      GetDistanceInPHOSPlane(emcclu, cpv, track,dx,dz) ;     
      if(TMath::Sqrt(dx*dx+dz*dz) < fRcpv ){ 
        new ((*fLinkUpArray)[iLinkUp++])  AliPHOSLink(dx, dz, iEmcRP, iCpv, track) ;
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

  //Make arrays to mark clusters already chosen
  Int_t index;

  Int_t * emcExist = 0;
  if(fEmcLast > fEmcFirst) {
    emcExist = new Int_t[fEmcLast-fEmcFirst] ;
    for(index = 0; index <fEmcLast-fEmcFirst; index ++)
      emcExist[index] = 1 ;
  }
  else 
    return;

  Bool_t * cpvExist = 0;
  if(fCpvLast > fCpvFirst) {
    cpvExist = new Bool_t[fCpvLast-fCpvFirst] ;
    for(index = 0; index <fCpvLast-fCpvFirst; index ++)
      cpvExist[index] = kTRUE ;
  }
 else 
    return;
  
  
  // Finds the smallest links and makes pairs of CPV and EMC clusters with smallest distance 
  TIter nextUp(fLinkUpArray) ;
  
  AliPHOSLink * linkUp ;
  
  AliPHOSCpvRecPoint * nullpointer = 0 ;
  
  while ( (linkUp =  static_cast<AliPHOSLink *>(nextUp()) ) ){  

    if(emcExist[linkUp->GetEmc()-fEmcFirst] != -1){

      if(cpvExist[linkUp->GetCpv()-fCpvFirst]){ //CPV still exist
         Float_t dx,dz ;
         linkUp->GetXZ(dx,dz) ;
	 new ((* fTrackSegments)[fNTrackSegments]) 
	   AliPHOSTrackSegment(static_cast<AliPHOSEmcRecPoint *>(fEMCRecPoints->At(linkUp->GetEmc())) , 
			       static_cast<AliPHOSCpvRecPoint *>(fCPVRecPoints->At(linkUp->GetCpv())) , 
			       linkUp->GetTrack(),dx,dz) ;
	 
       (static_cast<AliPHOSTrackSegment *>(fTrackSegments->At(fNTrackSegments)))->SetIndexInList(fNTrackSegments);
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
        Int_t track = -1 ;
        Float_t dx,dz ;
        AliPHOSEmcRecPoint *emcclu = dynamic_cast<AliPHOSEmcRecPoint *>(fEMCRecPoints->At(iEmcRP+fEmcFirst));
        GetDistanceInPHOSPlane(emcclu, 0, track,dx,dz);
        if(track<0)
	  new ((*fTrackSegments)[fNTrackSegments]) AliPHOSTrackSegment(emcclu,nullpointer) ;
        else
          new ((*fTrackSegments)[fNTrackSegments]) AliPHOSTrackSegment(emcclu,0,track,dx,dz);
	(dynamic_cast<AliPHOSTrackSegment *>(fTrackSegments->At(fNTrackSegments)))->SetIndexInList(fNTrackSegments);
	fNTrackSegments++;    
      } 
    }
  }
  delete [] emcExist ; 
  delete [] cpvExist ; 
}

//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::Clusters2TrackSegments(Option_t *option)
{
  // Steering method to perform track segment construction for the current event
  // Returns an array with the found track-segments.
  
  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSTSMaker");
 
  if(strstr(option,"print")) {
    Print() ; 
    return ; 
  }
  
  //Make some initializations 
  fNTrackSegments = 0 ;
  fEmcFirst = 0 ;    
  fEmcLast  = 0 ;   
  fCpvFirst = 0 ;   
  fCpvLast  = 0 ;   

  fTrackSegments->Clear();

  //   if(!ReadRecPoints(ievent))   continue; //reads RecPoints for event ievent

  for(fModule = 1; fModule <= fGeom->GetNModules() ; fModule++ ) {
    FillOneModule() ; 
    MakeLinks() ;
    MakePairs() ;
  }
    
  if(strstr(option,"deb"))
    PrintTrackSegments(option);

  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSTSMaker");
    Info("Exec", "took %f seconds for making TS", 
	 gBenchmark->GetCpuTime("PHOSTSMaker")); 
  }
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
void AliPHOSTrackSegmentMakerv1::PrintTrackSegments(Option_t * option)
{
  // option deb - prints # of found TrackSegments
  // option deb all - prints as well indexed of found RecParticles assigned to the TS

  Info("PrintTrackSegments", "Results from TrackSegmentMaker:") ; 
  printf("        Found %d TrackSegments\n", fTrackSegments->GetEntriesFast() ); 
  
  if(strstr(option,"all")) {  // printing found TS
    printf("TrackSegment #  EMC RP#  CPV RP#\n") ; 
    Int_t index;
    for (index = 0 ; index <fTrackSegments->GetEntriesFast() ; index++) {
      AliPHOSTrackSegment * ts = (AliPHOSTrackSegment * )fTrackSegments->At(index) ; 
      printf("   %d           %d        %d \n", ts->GetIndexInList(), ts->GetEmcIndex(), ts->GetCpvIndex() ) ; 
    }	
  }
}

