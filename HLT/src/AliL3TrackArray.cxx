//Author:        Uli Frankenfeld
//Last Modified: 06.12.2000

#include <math.h>
#include <string.h>
#include <iostream.h>
#include "AliL3Logging.h"
#include "AliL3TrackArray.h"
#include "AliL3HoughTrack.h"
#include "AliL3ConfMapTrack.h"
#include "AliL3TrackSegmentData.h"
#include "AliL3Transform.h"
#include "AliL3ConfMapPoint.h"
//_____________________________________________________________
//
// The L3 TrackArray 
//

ClassImp(AliL3TrackArray)

AliL3TrackArray::AliL3TrackArray(){
  //Default constructor
  fSize = 0;
  fNTracks=0;
  fNAbsent=0;
  fTrackType='t';
  SetSize();
}


AliL3TrackArray::AliL3TrackArray(Int_t ntrack){
  //Constructor.
  fSize = 0;
  fNTracks=0;
  fNAbsent=0;
  fTrackType='t';
  SetSize(ntrack);
}

AliL3TrackArray::AliL3TrackArray(char* tracktype,Int_t ntrack){
  //Constructor.
  fSize = 0;
  fNTracks=0;
  fNAbsent=0;
  if(strcmp(tracktype,"AliL3Track")==0) fTrackType='t';
  if(strcmp(tracktype,"AliL3ConfMapTrack")==0) fTrackType='c';
  if(strcmp(tracktype,"AliL3HoughTrack")==0) fTrackType='h';
  SetSize(ntrack);
}

AliL3TrackArray::AliL3TrackArray(char* tracktype){
  //Constructor.
  fSize = 0;
  fNTracks=0;
  fNAbsent=0;
  if(strcmp(tracktype,"AliL3Track")==0) fTrackType='t';
  if(strcmp(tracktype,"AliL3ConfMapTrack")==0) fTrackType='c';
  if(strcmp(tracktype,"AliL3HoughTrack")==0) fTrackType='h';
  SetSize();
}


AliL3TrackArray::~AliL3TrackArray(){
  //Destructor
  DeleteArray();
}


AliL3Track *AliL3TrackArray::NextTrack(){
  if(fNTracks<fSize) return fTrack[fNTracks++];
  SetSize(fSize+100);
   return fTrack[fNTracks++]; 
}

void AliL3TrackArray::DeleteArray(){
  for(Int_t i=0; i<fSize;i++)
    delete fTrack[i];
  delete[] fIsPresent;
  delete[] fTrack;
}

Bool_t AliL3TrackArray::SetSize(Int_t newsize){
  if(newsize<=fSize) return kFALSE; //shrink comes later!! 
  if(!fSize){
    fSize = newsize;
    fTrack = new AliL3Track*[fSize];
    fIsPresent = new Bool_t[fSize];
    switch(fTrackType){
      case 't':
        for(Int_t i=0;i<fSize;i++){
          fTrack[i]   = new AliL3Track();
          fIsPresent[i] = kTRUE;
        }
        break;
      case 'c':  
        for(Int_t i=0;i<fSize;i++){
          fTrack[i]   = new AliL3ConfMapTrack();
          fIsPresent[i] = kTRUE;
        }
        break;
      case 'h':
        for(Int_t i=0;i<fSize;i++){
          fTrack[i]   = new AliL3HoughTrack();
          fIsPresent[i] = kTRUE;
        }
        break;
      default: 
        return kFALSE;
    }
    return kTRUE;
  }
  AliL3Track **tmp = new AliL3Track*[fSize];
  Bool_t *pre = new Bool_t[fSize];
  for(Int_t i=0; i<fSize;i++){
    tmp[i] = fTrack[i];
    pre[i] = fIsPresent[i];
  }
  delete[]  fTrack;
  delete[] fIsPresent;
  fTrack =  new AliL3Track*[newsize];
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
        fTrack[i]   = new AliL3Track();
        fIsPresent[i] = kTRUE;
      }
      break;
    case 'c':  
      for(Int_t i=fSize;i<newsize;i++){
        fTrack[i]   = new AliL3ConfMapTrack();
        fIsPresent[i] = kTRUE;
      }
      break;
    case 'h':
      for(Int_t i=fSize;i<newsize;i++){
        fTrack[i]   = new AliL3HoughTrack();
        fIsPresent[i] = kTRUE;
      }
      break;
    default: 
      return kFALSE;
  }
  fSize = newsize;
  return kTRUE;
}

void AliL3TrackArray::Reset(){
  fNTracks=0;
  fNAbsent=0;
  for(Int_t i=0; i<fSize;i++)
    fIsPresent[i] = kTRUE; 
}

void AliL3TrackArray::Remove(Int_t track){
  if(fIsPresent[track]){
    fIsPresent[track]=kFALSE;
    fNAbsent++;
  }
}

void AliL3TrackArray::FillTracks(Int_t ntracks, AliL3TrackSegmentData* tr){
  //Read tracks from shared memory (or memory)
  AliL3TrackSegmentData *trs = tr;
  for(Int_t i=0; i<ntracks; i++){
    AliL3Track *track = NextTrack(); 
    track->SetPt(trs->fPt);
    track->SetPsi(trs->fPsi);
    track->SetTgl(trs->fTgl);
    track->SetNHits(trs->fNPoints);
    track->SetCharge(trs->fCharge);
    track->SetFirstPoint(trs->fX,trs->fY,trs->fZ);
    track->SetLastPoint(trs->fLastX,trs->fLastY,trs->fLastZ);
    track->SetHits( trs->fNPoints, trs->fPointIDs );
    UChar_t *tmpP = (UChar_t*)trs;
    tmpP += sizeof(AliL3TrackSegmentData)+trs->fNPoints*sizeof(UInt_t);
    trs = (AliL3TrackSegmentData*)tmpP;
  }
}

void AliL3TrackArray::FillTracks(Int_t ntracks, AliL3TrackSegmentData* tr,Int_t slice, AliL3Transform* trans){
  //Read tracks from shared memory (or memory)
  AliL3TrackSegmentData *trs = tr;
  for(Int_t i=0; i<ntracks; i++){
    AliL3Track *track = NextTrack(); 
    track->SetPt(trs->fPt);
    Float_t psi[1];
    psi[0]=trs->fPsi;
    trans->Local2GlobalAngle(psi,slice);
    track->SetPsi(psi[0]);
    track->SetTgl(trs->fTgl);
    track->SetNHits(trs->fNPoints);
    track->SetCharge(trs->fCharge);
    Float_t first[3];
    first[0]=trs->fX;first[1]=trs->fY;first[2]=trs->fZ;
    trans->Local2Global(first,slice);
    track->SetFirstPoint(first[0],first[1],first[2]);
    Float_t last[3];
    last[0]=trs->fLastX;last[1]=trs->fLastY;last[2]=trs->fLastZ;
    trans->Local2Global(last,slice);
    track->SetLastPoint(last[0],last[1],last[2]);
    track->SetHits( trs->fNPoints, trs->fPointIDs );
    UChar_t *tmpP = (UChar_t*)trs;
    tmpP += sizeof(AliL3TrackSegmentData)+trs->fNPoints*sizeof(UInt_t);
    trs = (AliL3TrackSegmentData*)tmpP;
  }
}

UInt_t AliL3TrackArray::GetOutSize(){
  UInt_t count = GetOutCount();   //use only present tracks
  UInt_t tHits = 0;
  for(Int_t i=0;i<fNTracks;i++){  //loop over all tracks
    AliL3Track *track = GetCheckedTrack(i);  //use only present tracks
    if(track)                                       //use only present tracks
      tHits += track->GetNHits();
  }

  //calculate size of track
  return count*sizeof(AliL3TrackSegmentData)+sizeof(UInt_t)*tHits;
}

UInt_t AliL3TrackArray::WriteTracks(UInt_t & ntracks,AliL3TrackSegmentData* tr){
  ntracks = GetOutCount();
  return WriteTracks(tr);
}

UInt_t AliL3TrackArray::WriteTracks(AliL3TrackSegmentData* tr){
  if(GetTrackType()=='c') return WriteConfMapTracks(tr);
  AliL3TrackSegmentData *tP = tr;
  UInt_t *pP;
  UInt_t size = 0;
  for(Int_t i=0; i<fNTracks; i++){  //loop over all tracks
    AliL3Track *track = GetCheckedTrack(i); //use only present tracks
    if(!track) continue;                           //use only present tracks
    tP->fX = track->GetFirstPointX();
    tP->fY = track->GetFirstPointY();
    tP->fZ = track->GetFirstPointZ();
    tP->fPt = track->GetPt();
    tP->fLastX = track->GetLastPointX();
    tP->fLastY = track->GetLastPointY();
    tP->fLastZ = track->GetLastPointZ();
    tP->fPsi = track->GetPsi();
    tP->fTgl = track->GetTgl();
    tP->fCharge = track->GetCharge();
    tP->fNPoints = track->GetNHits();
    pP = (UInt_t*)track->GetHitNumbers();
    for (UInt_t j=0;j<tP->fNPoints;j++){
      tP->fPointIDs[j] = pP[j];
    }
    Byte_t *tmpP = (Byte_t *)tP;
    tmpP += sizeof(AliL3TrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    size += sizeof(AliL3TrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    tP = (AliL3TrackSegmentData*)tmpP;
  }
  return size;
}

UInt_t AliL3TrackArray::WriteConfMapTracks(AliL3TrackSegmentData* tr){
  // use first and last point objects
  AliL3TrackSegmentData *tP = tr;
  UInt_t *pP;
  UInt_t size = 0;
  for(Int_t i=0; i<fNTracks; i++){  //loop over all tracks
    AliL3ConfMapTrack *track =(AliL3ConfMapTrack *) GetCheckedTrack(i); //use only present tracks
    if(!track) continue;                           //use only present tracks
    AliL3ConfMapPoint *hit = (AliL3ConfMapPoint*)track->lastHit;
    AliL3ConfMapPoint *lastHit = (AliL3ConfMapPoint*)track->firstHit;
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
//    tP->fLastX = track->GetLastPointX();
//    tP->fLastY = track->GetLastPointY();
//    tP->fLastZ = track->GetLastPointZ();
    tP->fPsi = track->GetPsi();
    tP->fTgl = track->GetTgl();
    tP->fCharge = track->GetCharge();
    tP->fNPoints = track->GetNHits();
    pP = (UInt_t*)track->GetHitNumbers();
    for (UInt_t j=0;j<tP->fNPoints;j++){
      tP->fPointIDs[j] = pP[j];
    }
    Byte_t *tmpP = (Byte_t *)tP;
    tmpP += sizeof(AliL3TrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    size +=sizeof(AliL3TrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    tP = (AliL3TrackSegmentData*)tmpP;
  }
  return size;
}

void AliL3TrackArray::AddTracks(AliL3TrackArray *newtrack){
  if(GetTrackType() != newtrack->GetTrackType())
    return;
  if(fSize < fNTracks+newtrack->GetNPresent())
    SetSize(fSize+newtrack->GetSize());
  for(Int_t i =0;i<newtrack->GetNTracks();i++){
    AliL3Track *tpt = newtrack->GetCheckedTrack(i);
    if(!tpt) continue;
    newtrack->Remove(i);
    AliL3Track *track = NextTrack();
    
    track->Set(tpt);

  }
}


void AliL3TrackArray::Compress(){
  if(GetNPresent()==GetNTracks()) return;
  AliL3Track **tmp =  new AliL3Track *[fNTracks];
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

void AliL3TrackArray::QSort(){
  // compress an sort
  Compress();
  QSort(fTrack,0,fNTracks);
}

 void AliL3TrackArray::QSort( AliL3Track **a, Int_t first, Int_t last){

   // Sort array of AliL3Track pointers using a quicksort algorithm.
   // Uses TrackCompare() to compare objects.
   // Thanks to Root! 

   static AliL3Track *tmp;
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

Int_t AliL3TrackArray::TrackCompare(AliL3Track *a, AliL3Track *b){
   // Compare the two tracks.

  if(a->GetNHits() < b->GetNHits()) return 1;
  if(a->GetNHits() > b->GetNHits()) return -1;
  return 0;
}


