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
 * Revision 1.5  2007/07/11 13:43:30  hristov
 * New class AliESDEvent, backward compatibility with the old AliESD (Christian)
 *
 * Revision 1.4  2007/05/04 14:49:29  policheh
 * AliPHOSRecPoint inheritance from AliCluster
 *
 * Revision 1.3  2007/04/25 19:39:42  kharlov
 * Track extracpolation improved
 *
 * Revision 1.2  2007/04/01 19:16:52  kharlov
 * D.P.: Produce EMCTrackSegments using TPC/ITS tracks (no CPV)
 *
 *
 */

//_________________________________________________________________________
// Implementation version 2 of algorithm class to construct PHOS track segments
// Track segment for PHOS is list of 
//        EMC RecPoint + (possibly) projection of TPC track
// To find TrackSegments we do the following: 
//  for each EMC RecPoints we look at
//   TPC projections radius fRtpc. 
//  If there is such a track
//   we make "Link" it is just indexes of EMC and TPC track and distance
//   between them in the PHOS plane. 
//  Then we sort "Links" and starting from the 
//   least "Link" pointing to the unassined EMC and TPC assing them to 
//   new TrackSegment. 
// If there is no TPC track we make TrackSegment 
// consisting from EMC alone. There is no TrackSegments without EMC RecPoint.
//// In principle this class should be called from AliPHOSReconstructor, but 
// one can use it as well in standalone mode.
// Use  case:
//  root [0] AliPHOSTrackSegmentMakerv2 * t = new AliPHOSTrackSegmentMaker("galice.root", "tracksegmentsname", "recpointsname")
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
#include "TFile.h"
#include "TTree.h"
#include "TBenchmark.h"

// --- Standard library ---
#include "Riostream.h"
// --- AliRoot header files ---
#include "AliPHOSGeometry.h"
#include "AliPHOSTrackSegmentMakerv2.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSLink.h"
#include "AliPHOSGetter.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"

ClassImp( AliPHOSTrackSegmentMakerv2) 


//____________________________________________________________________________
AliPHOSTrackSegmentMakerv2::AliPHOSTrackSegmentMakerv2() :
  AliPHOSTrackSegmentMaker(),
  fDefaultInit(kTRUE),
  fWrite(kFALSE),
  fNTrackSegments(0),
  fRtpc(0.f),
  fVtx(0.,0.,0.),
  fLinkUpArray(0),
  fTPCtracks(),
  fEmcFirst(0),
  fEmcLast(0),
  fModule(0),
  fTrackSegmentsInRun(0)
{
  // default ctor (to be used mainly by Streamer)
  InitParameters() ; 
}

//____________________________________________________________________________
AliPHOSTrackSegmentMakerv2::AliPHOSTrackSegmentMakerv2(const TString & alirunFileName, const TString & eventFolderName) :
  AliPHOSTrackSegmentMaker(alirunFileName, eventFolderName),
  fDefaultInit(kFALSE),
  fWrite(kFALSE),
  fNTrackSegments(0),
  fRtpc(0.f),
  fVtx(0.,0.,0.),
  fLinkUpArray(0),
  fTPCtracks(),
  fEmcFirst(0),
  fEmcLast(0),
  fModule(0),
  fTrackSegmentsInRun(0)
{
  // ctor
  InitParameters() ; 
  Init() ;
  fESD = 0;
}

//____________________________________________________________________________
AliPHOSTrackSegmentMakerv2::AliPHOSTrackSegmentMakerv2(const AliPHOSTrackSegmentMakerv2 & tsm) :
  AliPHOSTrackSegmentMaker(tsm),
  fDefaultInit(kFALSE),
  fWrite(kFALSE),
  fNTrackSegments(0),
  fRtpc(0.f),
  fVtx(0.,0.,0.),
  fLinkUpArray(0),
  fTPCtracks(),
  fEmcFirst(0),
  fEmcLast(0),
  fModule(0),
  fTrackSegmentsInRun(0)
{
  // cpy ctor: no implementation yet
  // requested by the Coding Convention
  Fatal("cpy ctor", "not implemented") ;
}


//____________________________________________________________________________
 AliPHOSTrackSegmentMakerv2::~AliPHOSTrackSegmentMakerv2()
{ 
  // dtor
  // fDefaultInit = kTRUE if TrackSegmentMaker created by default ctor (to get just the parameters)
  if (!fDefaultInit)  
    delete fLinkUpArray ;
}

//____________________________________________________________________________
const TString AliPHOSTrackSegmentMakerv2::BranchName() const 
{  
 
  return GetName() ;
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv2::FillOneModule()
{
  // Finds first and last indexes between which 
  // clusters from one PHOS module are

  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  const AliPHOSGeometry * geom = gime->PHOSGeometry() ; 
  
  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
 
  //First EMC clusters
  Int_t totalEmc = emcRecPoints->GetEntriesFast() ;
  for(fEmcFirst = fEmcLast; (fEmcLast < totalEmc) &&  
	((dynamic_cast<AliPHOSRecPoint *>(emcRecPoints->At(fEmcLast)))->GetPHOSMod() == fModule ); 
      fEmcLast ++)  ;
  
  //Now TPC tracks
  if(fESD){
    //Do it ones, only first time
    if(fModule==1){
      Int_t nTracks = fESD->GetNumberOfTracks();

      Int_t nPHOSmod = geom->GetNModules() ;
      if(fTPCtracks[0].size()<(UInt_t)nTracks){
        for(Int_t i=0; i<nPHOSmod; i++) 
          fTPCtracks[i].resize(nTracks) ;
      }
      for(Int_t i=0; i<5; i++)fNtpcTracks[i]=0 ;
      TVector3 inPHOS ;
   
      //In this particular case we use fixed vertex position at zero
      Double_t vtx[3]={0.,0.,0.} ;
      AliESDtrack *track;
      Double_t xyz[3] ;
      Double_t rEMC = geom->GetIPtoCrystalSurface() ; //Use here ideal geometry 
      for (Int_t iTrack=0; iTrack<nTracks; iTrack++) {
        track = fESD->GetTrack(iTrack);
        for(Int_t iTestMod=1; iTestMod<=nPHOSmod; iTestMod++){
           Double_t modAngle=270.+geom->GetPHOSAngle(iTestMod) ;
           modAngle*=TMath::Pi()/180. ;
           track->Rotate(modAngle);  
           if (!track->GetXYZAt(rEMC, fESD->GetMagneticField(), xyz))
             continue; //track coord on the cylinder of PHOS radius
           if ((TMath::Abs(xyz[0])+TMath::Abs(xyz[1])+TMath::Abs(xyz[2]))<=0)
             continue;
           //Check if this track hits PHOS
           inPHOS.SetXYZ(xyz[0],xyz[1],xyz[2]);
           Int_t modNum ; 
           Double_t x,z ;
           geom->ImpactOnEmc(vtx, inPHOS.Theta(), inPHOS.Phi(), modNum, z, x) ;
           if(modNum==iTestMod){
             //Mark this track as one belonging to module
             TrackInPHOS_t &t = fTPCtracks[modNum-1][fNtpcTracks[modNum-1]++] ; 
             t.track = track ;
             t.x = x ;
             t.z = z ;
             break ;
           }
        }
      }
    }
  } 

}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv2::GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * emcClu,
                                                         AliESDtrack *track,
                                                         Float_t &dx, Float_t &dz) const
{
  // Calculates the distance between the EMC RecPoint and the CPV RecPoint
  // Clusters are sorted in "rows" and "columns" of width 1 cm

  const AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;
  TVector3 emcGlobal; // Global position of the EMC recpoint
  geom->GetGlobalPHOS((AliPHOSRecPoint*)emcClu,emcGlobal);

  //Calculate actual distance to the PHOS surface
  Double_t modAngle=270.+geom->GetPHOSAngle(emcClu->GetPHOSMod()) ;
  modAngle*=TMath::Pi()/180. ;
  Double_t rEMC = emcGlobal.X()*TMath::Cos(modAngle)+emcGlobal.Y()*TMath::Sin(modAngle) ;
  track->Rotate(modAngle);
  Double_t xyz[3] ;
  if(!track->GetXYZAt(rEMC, fESD->GetMagneticField(), xyz)){
    dx=999. ;
    dz=999. ;
    return ; //track coord on the cylinder of PHOS radius
  }
  if((TMath::Abs(xyz[0])+TMath::Abs(xyz[1])+TMath::Abs(xyz[2]))<=0){
    dx=999. ; 
    dz=999. ;
    return ;
  } 

  dx=(emcGlobal.X()-xyz[0])*TMath::Sin(modAngle)-(emcGlobal.Y()-xyz[1])*TMath::Cos(modAngle) ;
  dx=TMath::Sign(dx,(Float_t)(emcGlobal.X()-xyz[0])) ; //set direction
  dz=emcGlobal.Z()-xyz[2] ;

  return ;
}
//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv2::Init()
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
void  AliPHOSTrackSegmentMakerv2::InitParameters()
{
  //Initializes parameters
  fRtpc      = 4. ;
  fEmcFirst  = 0 ;    
  fEmcLast   = 0 ;   
  fLinkUpArray = 0 ;
  fWrite                   = kTRUE ;
  fTrackSegmentsInRun       = 0 ; 
  SetEventRange(0,-1) ;
}


//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv2::MakeLinks()
{ 
  // Finds distances (links) between all EMC and CPV clusters, 
  // which are not further apart from each other than fRcpv 
  // and sort them in accordance with this distance
  
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 

  fLinkUpArray->Clear() ;    

  AliPHOSEmcRecPoint * emcclu ;

  Int_t iLinkUp  = 0 ;
  
  Int_t iEmcRP;
  for(iEmcRP = fEmcFirst; iEmcRP < fEmcLast; iEmcRP++ ) {
    emcclu = dynamic_cast<AliPHOSEmcRecPoint *>(emcRecPoints->At(iEmcRP)) ;
    TVector3 vecEmc ;
    emcclu->GetLocalPosition(vecEmc) ;
    Int_t mod=emcclu->GetPHOSMod() ;
    for(Int_t itr=0; itr<fNtpcTracks[mod-1] ; itr++){
      TrackInPHOS_t &t = fTPCtracks[mod-1][itr] ;
      //calculate raw distance
      Double_t rawdx = t.x - vecEmc.X() ;
      Double_t rawdz = t.z - vecEmc.Z() ;
      if(TMath::Sqrt(rawdx*rawdx+rawdz*rawdz)<fRtpc){
        Float_t dx,dz ;
        //calculate presize difference accounting misalignment
        GetDistanceInPHOSPlane(emcclu, t.track, dx,dz) ;     
        Int_t itrack = t.track->GetID() ;
        new ((*fLinkUpArray)[iLinkUp++])  AliPHOSLink(dx, dz, iEmcRP, itrack, -1) ;
      }      
    }
  } 
  fLinkUpArray->Sort() ;  //first links with smallest distances
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv2::MakePairs()
{ 
  // Using the previously made list of "links", we found the smallest link - i.e. 
  // link with the least distance between EMC and CPV and pointing to still 
  // unassigned RecParticles. We assign these RecPoints to TrackSegment and 
  // remove them from the list of "unassigned". 

  AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 

  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
  TClonesArray * trackSegments = gime->TrackSegments();
    
  //Make arrays to mark clusters already chosen
  Int_t * emcExist = 0;
  if(fEmcLast > fEmcFirst)
    emcExist = new Int_t[fEmcLast-fEmcFirst] ;
  
  Int_t index;
  for(index = 0; index <fEmcLast-fEmcFirst; index ++)
    emcExist[index] = 1 ;
  
  
  Bool_t * tpcExist = 0;
  Int_t nTracks = fESD->GetNumberOfTracks();
  if(nTracks>0)
    tpcExist = new Bool_t[nTracks] ;
  
  for(index = 0; index <nTracks; index ++)
    tpcExist[index] = 1 ;
  
  
  // Finds the smallest links and makes pairs of CPV and EMC clusters with smallest distance 
  TIter nextUp(fLinkUpArray) ;
  
  AliPHOSLink * linkUp ;
  
  AliPHOSCpvRecPoint * nullpointer = 0 ;
  
  while ( (linkUp =  static_cast<AliPHOSLink *>(nextUp()) ) ){  

    if(emcExist[linkUp->GetEmc()-fEmcFirst] != -1){

      if(tpcExist[linkUp->GetCpv()]){ //Track still exist
         Float_t dx,dz ; 
         linkUp->GetXZ(dx,dz) ;
	 new ((* trackSegments)[fNTrackSegments]) 
	   AliPHOSTrackSegment(dynamic_cast<AliPHOSEmcRecPoint *>(emcRecPoints->At(linkUp->GetEmc())) , 
			       nullpointer,
			       linkUp->GetTrack(),dx,dz) ;
       (dynamic_cast<AliPHOSTrackSegment *>(trackSegments->At(fNTrackSegments)))->SetIndexInList(fNTrackSegments);
       fNTrackSegments++ ;
       emcExist[linkUp->GetEmc()-fEmcFirst] = -1 ; //Mark emc  that Cpv was found 
       //mark track as already used 
       tpcExist[linkUp->GetCpv()] = kFALSE ;
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
  delete [] tpcExist ; 
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv2::Exec(Option_t *option)
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
    gime->Event(ievent,"DR") ;
   //Make some initializations 
    fNTrackSegments = 0 ;
    fEmcFirst = 0 ;    
    fEmcLast  = 0 ;   
    
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
void AliPHOSTrackSegmentMakerv2::Unload() 
{
  // Unloads the task from the folder
  AliPHOSGetter * gime = AliPHOSGetter::Instance() ;  
  gime->PhosLoader()->UnloadRecPoints() ;
  gime->PhosLoader()->UnloadTracks() ;
}

//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv2::Print(const Option_t *)const
{
  //  Print TrackSegmentMaker parameters

  TString message("") ;
  if( strcmp(GetName(), "") != 0 ) {
    message = "\n======== AliPHOSTrackSegmentMakerv2 ========\n" ; 
    message += "Making Track segments\n" ;
    message += "with parameters:\n" ; 
    message += "     Maximal EMC - TPC distance (cm) %f\n" ;
    message += "============================================\n" ;
    Info("Print", message.Data(),fRtpc) ;
  }
  else
    Info("Print", "AliPHOSTrackSegmentMakerv2 not initialized ") ;
}

//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv2::WriteTrackSegments()
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
void AliPHOSTrackSegmentMakerv2::PrintTrackSegments(Option_t * option)
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
      printf("   %d           %d        %d \n", ts->GetIndexInList(), ts->GetEmcIndex(), ts->GetTrackIndex() ) ; 
    }	
  }
}
