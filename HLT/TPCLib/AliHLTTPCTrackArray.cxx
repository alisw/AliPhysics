// @(#) $Id$
// Original: AliHLTTrackArray.cxx,v 1.21 2005/06/14 10:55:21 cvetan 

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Uli Frankenfeld, maintained by                          *
 *                  Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCTrackArray.cxx
    @author Uli Frankenfeld, maintained by Matthias Richter
    @date   
    @brief  Array of AliHLTTPCTracks */

#include "AliLog.h"
#include "TClass.h"
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCTrackArray.h"
// Matthias 17.10.2007 the hough code has been disabled for the moment
//#define INCLUDE_TPC_HOUGH
#ifdef INCLUDE_TPC_HOUGH
#include "AliHLTTPCHoughTrack.h"
#endif
#include "AliHLTTPCModelTrack.h"
#include "AliHLTTPCConfMapTrack.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCConfMapPoint.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCTrackArray)

AliHLTTPCTrackArray::AliHLTTPCTrackArray()
  :
  fTrackType('t'),
  fSize(0),
  fIsPresent(NULL),
  fNAbsent(0),
  fTrack(NULL),
  fNTracks(0)
{
  //Default constructor
  SetSize();
}


AliHLTTPCTrackArray::AliHLTTPCTrackArray(Int_t ntrack)
  :
  fTrackType('t'),
  fSize(0),
  fIsPresent(NULL),
  fNAbsent(0),
  fTrack(NULL),
  fNTracks(0)
{
  //Constructor.
  SetSize(ntrack);
}

AliHLTTPCTrackArray::AliHLTTPCTrackArray(char* tracktype,Int_t ntrack)
  :
  fTrackType('t'),
  fSize(0),
  fIsPresent(NULL),
  fNAbsent(0),
  fTrack(NULL),
  fNTracks(0)
{
  //Constructor.
  if(strcmp(tracktype,"AliHLTTPCTrack")==0) fTrackType='t';
  if(strcmp(tracktype,"AliHLTTPCConfMapTrack")==0) fTrackType='c';
#ifdef INCLUDE_TPC_HOUGH
  if(strcmp(tracktype,"AliHLTTPCHoughTrack")==0) fTrackType='h';
#endif
  if(strcmp(tracktype,"AliHLTTPCModelTrack")==0) fTrackType='m';
  SetSize(ntrack);
}

AliHLTTPCTrackArray::AliHLTTPCTrackArray(char* tracktype)
  :
  fTrackType('t'),
  fSize(0),
  fIsPresent(NULL),
  fNAbsent(0),
  fTrack(NULL),
  fNTracks(0)
{
  //Constructor.
  if(strcmp(tracktype,"AliHLTTPCTrack")==0) fTrackType='t';
  if(strcmp(tracktype,"AliHLTTPCConfMapTrack")==0) fTrackType='c';
#ifdef INCLUDE_TPC_HOUGH
  if(strcmp(tracktype,"AliHLTTPCHoughTrack")==0) fTrackType='h';
#endif
  if(strcmp(tracktype,"AliHLTTPCModelTrack")==0) fTrackType='m';
  SetSize();
}

AliHLTTPCTrackArray::~AliHLTTPCTrackArray()
{
  //Destructor
  DeleteArray();
}


AliHLTTPCTrack *AliHLTTPCTrackArray::NextTrack()
{
  //next track in array
  if(fNTracks<fSize) return fTrack[fNTracks++];
  SetSize(fSize+100);
   return fTrack[fNTracks++]; 
}

void AliHLTTPCTrackArray::DeleteArray()
{
  //delete array
  for(Int_t i=0; i<fSize;i++)
    delete fTrack[i];
  delete[] fIsPresent;
  delete[] fTrack;
  fSize=0;
}

Bool_t AliHLTTPCTrackArray::SetSize(Int_t newsize)
{
  //set size
  if(newsize<=fSize) return kFALSE; //shrink comes later!! 
  if(!fSize){
    fSize = newsize;
    fTrack = new AliHLTTPCTrack*[fSize];
    fIsPresent = new Bool_t[fSize];
    switch(fTrackType){
      case 't':
        for(Int_t i=0;i<fSize;i++){
          fTrack[i]   = new AliHLTTPCTrack();
          fIsPresent[i] = kTRUE;
        }
        break;
      case 'c':  
        for(Int_t i=0;i<fSize;i++){
          fTrack[i]   = new AliHLTTPCConfMapTrack();
          fIsPresent[i] = kTRUE;
        }
        break;
#ifdef INCLUDE_TPC_HOUGH
      case 'h':
        for(Int_t i=0;i<fSize;i++){
          fTrack[i]   = new AliHLTTPCHoughTrack();
          fIsPresent[i] = kTRUE;
        }
        break;
#endif
       case 'm':
        for(Int_t i=0;i<fSize;i++){
          fTrack[i]   = new AliHLTTPCModelTrack();
          fIsPresent[i] = kTRUE;
        }
        break;
      default: 
        return kFALSE;
    }
    return kTRUE;
  }
  AliHLTTPCTrack **tmp = new AliHLTTPCTrack*[fSize];
  Bool_t *pre = new Bool_t[fSize];
  for(Int_t i=0; i<fSize;i++){
    tmp[i] = fTrack[i];
    pre[i] = fIsPresent[i];
  }
  delete[]  fTrack;
  delete[] fIsPresent;
  fTrack =  new AliHLTTPCTrack*[newsize];
  fIsPresent = new Bool_t[newsize];
  for(Int_t i=0; i<fSize;i++){
    fTrack[i]   = tmp[i];
    fIsPresent[i] = pre[i];
  }
  delete[] tmp;
  delete[] pre;
  switch(fTrackType){
    case 't':
      for(Int_t i=fSize;i<newsize;i++){
        fTrack[i]   = new AliHLTTPCTrack();
        fIsPresent[i] = kTRUE;
      }
      break;
    case 'c':  
      for(Int_t i=fSize;i<newsize;i++){
        fTrack[i]   = new AliHLTTPCConfMapTrack();
        fIsPresent[i] = kTRUE;
      }
      break;
#ifdef INCLUDE_TPC_HOUGH
    case 'h':
      for(Int_t i=fSize;i<newsize;i++){
        fTrack[i]   = new AliHLTTPCHoughTrack();
        fIsPresent[i] = kTRUE;
      }
      break;
#endif
    case 'm':
      for(Int_t i=fSize;i<newsize;i++){
        fTrack[i]   = new AliHLTTPCModelTrack();
        fIsPresent[i] = kTRUE;
      }
      break;
    default: 
      return kFALSE;
  }
  fSize = newsize;
  return kTRUE;
}

void AliHLTTPCTrackArray::Reset()
{
  //reset
  fNTracks=0;
  fNAbsent=0;
  for(Int_t i=0; i<fSize;i++)
    fIsPresent[i] = kTRUE; 
}

void AliHLTTPCTrackArray::Remove(Int_t track)
{
  //remove track
  if(fIsPresent[track]){
    fIsPresent[track]=kFALSE;
    fNAbsent++;
  }
}

void AliHLTTPCTrackArray::FillTracks(Int_t ntracks, AliHLTTPCTrackSegmentData* tr,Int_t slice, Int_t bTransform)
{
  //Read tracks from shared memory (or memory)
  AliHLTTPCTrackSegmentData *trs = tr;
  for(Int_t i=0; i<ntracks; i++){
    AliHLTTPCTrack *track = NextTrack(); 
    track->SetPt(trs->fPt);
    track->SetPterr(trs->fPterr);
    Float_t psi[1];
    psi[0]=trs->fPsi;
    if (slice>=0 && bTransform!=0)  {
      AliHLTTPCTransform::Local2GlobalAngle(psi,slice);
    }
    //cout << "psi " << psi[0] << endl;
    track->SetPsi(psi[0]);
    track->SetTgl(trs->fTgl);
    track->SetPsierr(trs->fPsierr);
    track->SetTglerr(trs->fTglerr);
    track->SetNHits(trs->fNPoints);
    track->SetCharge(trs->fCharge);
    Float_t first[3];
    first[0]=trs->fX;first[1]=trs->fY;first[2]=trs->fZ;
    if (slice>=0 && bTransform!=0)  {
      AliHLTTPCTransform::Local2Global(first,slice);
    }
    //cout << "first point: " << first[0] << " " << first[1] << " " << first[3] << endl;
    track->SetFirstPoint(first[0],first[1],first[2]);
    Float_t last[3];
    last[0]=trs->fLastX;last[1]=trs->fLastY;last[2]=trs->fLastZ;
    if (slice>=0 && bTransform!=0)  {
      AliHLTTPCTransform::Local2Global(last,slice);
    }
    //cout << "last point: " << last[0] << " " << last[1] << " " << last[3] << endl;
    track->SetLastPoint(last[0],last[1],last[2]);
    track->SetHits( trs->fNPoints, trs->fPointIDs );

    //if (slice>=0 && bTransform!=0)  {
      // Matthias Feb07: as everything is now in global coordinates, sector should
      // be set to 0. But as the display does a check on the sector, we have to set
      // it to the slice no. I suspect, that the transformation is done twice.
      //track->SetSector(0);
      track->SetSector(slice);
    //} else {
      // the parameters are in local coordinates, set the sector no
      //#ifndef INCLUDE_TPC_HOUGH
      //if (slice<0) track->SetSector(0);
      //else track->SetSector(slice);
      //#else 
      // Matthias Feb 2007: this is some kind of legacy ...
      // the INCLUDE_TPC_HOUGH has never been switched on in the new TPCLib
      // and this line was below in the corresponding block. As the slice
      // parameter is very useful but not available if the define is off
      // we distinguish the two cases here. Should be cleaned up.
      // Matthias April 2007: update, try to integrate Cvetans Hough tracker
      // so we need the support for the AliHLTTPCHoughTrack. I dont have the
      // full control of this code (should we use slice or trs->fSector?)
      // But the FillTracks method is never called from the hough code, so we
      // take 'slice'
      if (GetTrackType()=='h') {
	AliErrorClassStream() << "FillTracks was never used with AliHLTTPCHoughTrack:" 
			   << " CHECK THIS CODE!!!" << endl;
      }
      //track->SetSector(trs->fSector);
      //#endif // INCLUDE_TPC_HOUGH
      //}

    // this is currently a quick hack for straight lines of the first version 
    // of the CA tracker.
    // we have to think about a more general way of treating straight and curved
    // tracks
    if ( trs->fPt == -9876.0 ||  trs->fPt == -1.0) {
      track->SetPhi0(atan2(first[1],first[0]));
      track->SetKappa(1.0);
      track->SetRadius(999999.0);
    } else {
      // Matthias Feb07: just tried to take this away, but this causes the tracks
      // in the display not to be drawn. But we still have to tink about this.
      track->CalculateHelix();
    }

#ifdef INCLUDE_TPC_HOUGH
#ifdef ROWHOUGHPARAMS
    if(GetTrackType()=='h') {
      ((AliHLTTPCHoughTrack *)track)->SetWeight(trs->fWeight);
      ((AliHLTTPCHoughTrack *)track)->SetBinXY(trs->fBinX,trs->fBinY,trs->fBinXSize,trs->fBinYSize);
    }
    track->SetMCid(trs->fTrackID);
    track->SetRowRange(trs->fRowRange1,trs->fRowRange2);
    track->SetPID(trs->fPID);
#endif
#endif // INCLUDE_TPC_HOUGH
    UChar_t *tmpP = (UChar_t*)trs;
    tmpP += sizeof(AliHLTTPCTrackSegmentData)+trs->fNPoints*sizeof(UInt_t);
    trs = (AliHLTTPCTrackSegmentData*)tmpP;
  }
}

UInt_t AliHLTTPCTrackArray::GetOutSize()
{
  //get size for IO
  UInt_t count = GetOutCount();   //use only present tracks
  UInt_t tHits = 0;
  for(Int_t i=0;i<fNTracks;i++){  //loop over all tracks
    AliHLTTPCTrack *track = GetCheckedTrack(i);  //use only present tracks
    if(track)                                       //use only present tracks
      tHits += track->GetNHits();
  }

  //calculate size of track
  return count*sizeof(AliHLTTPCTrackSegmentData)+sizeof(UInt_t)*tHits;
}

UInt_t AliHLTTPCTrackArray::WriteTracks(UInt_t & ntracks,AliHLTTPCTrackSegmentData* tr)
{
  //write tracks
  ntracks = GetOutCount();
  return WriteTracks(tr);
}

UInt_t AliHLTTPCTrackArray::WriteTracks(AliHLTTPCTrackSegmentData* tr)
{
  //if(GetTrackType()=='c') return WriteConfMapTracks(tr);
  AliHLTTPCTrackSegmentData *tP = tr;
  UInt_t *pP;
  UInt_t size = 0;
  for(Int_t i=0; i<fNTracks; i++){  //loop over all tracks
    AliHLTTPCTrack *track = GetCheckedTrack(i); //use only present tracks
    if(!track) continue;                           //use only present tracks
    tP->fX = track->GetFirstPointX();
    tP->fY = track->GetFirstPointY();
    tP->fZ = track->GetFirstPointZ();
    tP->fPt = track->GetPt();
    tP->fPterr = track->GetPterr();
    tP->fLastX = track->GetLastPointX();
    tP->fLastY = track->GetLastPointY();
    tP->fLastZ = track->GetLastPointZ();
    tP->fPsi = track->GetPsi();
    tP->fTgl = track->GetTgl();
    tP->fPsierr = track->GetPsierr();
    tP->fTglerr = track->GetTglerr();
    tP->fCharge = track->GetCharge();
    tP->fNPoints = track->GetNHits();
#ifdef INCLUDE_TPC_HOUGH
#ifdef ROWHOUGHPARAMS
    if(GetTrackType()=='h') {
      tP->fWeight = ((AliHLTTPCHoughTrack *)track)->GetWeight();
      tP->fBinX = ((AliHLTTPCHoughTrack *)track)->GetBinX();
      tP->fBinY = ((AliHLTTPCHoughTrack *)track)->GetBinY();
      tP->fBinXSize = ((AliHLTTPCHoughTrack *)track)->GetSizeX();
      tP->fBinYSize = ((AliHLTTPCHoughTrack *)track)->GetSizeY();
    }
    tP->fTrackID = track->GetMCid();
    tP->fRowRange1 = track->GetFirstRow();
    tP->fRowRange2 = track->GetLastRow();
    tP->fSector = track->GetSector();
    tP->fPID = track->GetPID();
#endif
#endif // INCLUDE_TPC_HOUGH
    pP = (UInt_t*)track->GetHitNumbers();
    for (UInt_t j=0;j<tP->fNPoints;j++){
      tP->fPointIDs[j] = pP[j];
    }
    Byte_t *tmpP = (Byte_t *)tP;
    tmpP += sizeof(AliHLTTPCTrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    size += sizeof(AliHLTTPCTrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    tP = (AliHLTTPCTrackSegmentData*)tmpP;

//    LOG(AliHLTTPCLog::kError,"AliHLTTPCTrackArray::WriteTracks","TRACKPARAMETER") <<ENDLOG;
//    track->Rotate(0,kFALSE);
//    track->Print();

  }
  return size;
}

UInt_t AliHLTTPCTrackArray::WriteConfMapTracks(AliHLTTPCTrackSegmentData* tr)
{
  // use first and last point objects
  AliHLTTPCTrackSegmentData *tP = tr;
  UInt_t *pP;
  UInt_t size = 0;
  for(Int_t i=0; i<fNTracks; i++){  //loop over all tracks
    AliHLTTPCConfMapTrack *track =(AliHLTTPCConfMapTrack *) GetCheckedTrack(i); //use only present tracks
    if(!track) continue;                           //use only present tracks
    AliHLTTPCConfMapPoint *hit = (AliHLTTPCConfMapPoint*)track->GetLastHit();
    AliHLTTPCConfMapPoint *lastHit = (AliHLTTPCConfMapPoint*)track->GetFirstHit();
    tP->fX = hit->GetX();
    tP->fY = hit->GetY();
    tP->fZ = hit->GetZ();
    tP->fLastX = lastHit->GetX();
    tP->fLastY = lastHit->GetY();
    tP->fLastZ = lastHit->GetZ();
   
//    tP->fX = track->GetFirstPointX();
//    tP->fY = track->GetFirstPointY();
//    tP->fZ = track->GetFirstPointZ();
    tP->fPt = track->GetPt();
    tP->fPterr = track->GetPterr();
//    tP->fLastX = track->GetLastPointX();
//    tP->fLastY = track->GetLastPointY();
//    tP->fLastZ = track->GetLastPointZ();
    tP->fPsi = track->GetPsi();
    tP->fTgl = track->GetTgl();
    tP->fPsierr = track->GetPsierr();
    tP->fTglerr = track->GetTglerr();
    tP->fCharge = track->GetCharge();
#ifdef INCLUDE_TPC_HOUGH
#ifdef ROWHOUGHPARAMS
    tP->fTrackID = track->GetMCid();
    tP->fRowRange1 = track->GetFirstRow();
    tP->fRowRange2 = track->GetLastRow();
    tP->fSector = track->GetSector();
    tP->fPID = track->GetPID();
#endif
#endif // INCLUDE_TPC_HOUGH
    tP->fNPoints = track->GetNHits();
    pP = (UInt_t*)track->GetHitNumbers();
    for (UInt_t j=0;j<tP->fNPoints;j++){
      tP->fPointIDs[j] = pP[j];
    }
    Byte_t *tmpP = (Byte_t *)tP;
    tmpP += sizeof(AliHLTTPCTrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    size +=sizeof(AliHLTTPCTrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    tP = (AliHLTTPCTrackSegmentData*)tmpP;
  }
  return size;
}

void AliHLTTPCTrackArray::AddLast(AliHLTTPCTrack *track)
{
  //add track to last position
  AliHLTTPCTrack *tpt = NextTrack();
  tpt->Copy(track);
  
}

void AliHLTTPCTrackArray::AddTracks(AliHLTTPCTrackArray *newtrack,Bool_t remove_old,Int_t slice)
{
  //add tracks
  if(GetTrackType() != newtrack->GetTrackType() && GetTrackType()!='t')
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCTrackArray::AddTracks","Track types")
	<<"Bad idea to add tracks of different types"<<ENDLOG;
      return;
    }
  if(fSize < fNTracks+newtrack->GetNPresent())
    SetSize(fSize+newtrack->GetSize());
  for(Int_t i =0;i<newtrack->GetNTracks();i++){
    AliHLTTPCTrack *tpt = newtrack->GetCheckedTrack(i);
    if(!tpt) continue;
    if(remove_old)
      newtrack->Remove(i);
    AliHLTTPCTrack *track = NextTrack();
    track->Copy(tpt);
    if(slice>=0)
      track->Rotate(slice); //Rotate track to global coordinates
    /*
      AliHLTTPCTrack *track;
#ifdef INCLUDE_TPC_HOUGH
      if(GetTrackType()=='h')
      track = (AliHLTTPCHoughTrack*)NextTrack();
      else
#endif
      track = NextTrack();
      track->Copy(tpt);
    */
  }
}

void AliHLTTPCTrackArray::Compress()
{
  //compress array
  if(GetNPresent()==GetNTracks()) return;
  AliHLTTPCTrack **tmp =  new AliHLTTPCTrack *[fNTracks];
  Int_t present=0;
  Int_t absent=GetNPresent();
  for(Int_t i=0;i<GetNTracks();i++){
    if(fIsPresent[i]) tmp[present++] = fTrack[i];
    else tmp[absent++] = fTrack[i];
  }
  for(Int_t i=0;i<GetNTracks();i++)
    fIsPresent[i]=kTRUE;

  //Copy pointers back
  for(Int_t i=0; i<GetNTracks();i++){
    fTrack[i]=tmp[i];
  }

  delete[] tmp;

  fNTracks = GetNPresent();
  fNAbsent = 0;
}

void AliHLTTPCTrackArray::QSort()
{
  // compress and sort
  Compress();
  QSort(fTrack,0,fNTracks);
}

void AliHLTTPCTrackArray::QSort( AliHLTTPCTrack **a, Int_t first, Int_t last)
{
   // Sort array of AliHLTTPCTrack pointers using a quicksort algorithm.
   // Uses TrackCompare() to compare objects.
   // Thanks to Root! 

   static AliHLTTPCTrack *tmp;
   static int i;           // "static" to save stack space
   int j;

   while (last - first > 1) {
      i = first;
      j = last;
      for (;;) {
	while (++i < last && TrackCompare(a[i], a[first]) < 0)
	  ;
	while (--j > first && TrackCompare(a[j], a[first]) > 0)
	  ;
         if (i >= j)
            break;

         tmp  = a[i];
         a[i] = a[j];
         a[j] = tmp;
      }
      if (j == first) {
         ++first;
         continue;
      }
      tmp = a[first];
      a[first] = a[j];
      a[j] = tmp;
      if (j - first < last - (j + 1)) {
         QSort(a, first, j);
         first = j + 1;   // QSort(j + 1, last);
      } else {
         QSort(a, j + 1, last);
         last = j;        // QSort(first, j);
      }
   }
}

Int_t AliHLTTPCTrackArray::TrackCompare(AliHLTTPCTrack *a, AliHLTTPCTrack *b) const
{
   // Compare the two tracks.
  
  return b->Compare(a);
  
  /*
#ifdef INCLUDE_TPC_HOUGH
    if(fTrackType=='h')
    {
    AliHLTTPCHoughTrack *tra = (AliHLTTPCHoughTrack*)a;
    AliHLTTPCHoughTrack *trb = (AliHLTTPCHoughTrack*)b;
    if(tra->GetWeight() < trb->GetWeight()) return 1;
    if(tra->GetWeight() > trb->GetWeight()) return -1;
    }
    else
#endif
    {
    if(a->GetNHits() < b->GetNHits()) return 1;
    if(a->GetNHits() > b->GetNHits()) return -1;
    }
    
    return 0;
  */
}

AliHLTTPCTrack* AliHLTTPCTrackArray::operator[](int index)
{
  // access operator
  if (index<fNTracks) return fTrack[index];
  return NULL;
}
