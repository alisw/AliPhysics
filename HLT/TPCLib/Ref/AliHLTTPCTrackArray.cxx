// @(#) $Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCStandardIncludes.h"

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCHoughTrack.h"
#include "AliHLTTPCModelTrack.h"
#include "AliHLTTPCConfMapTrack.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCConfMapPoint.h"

/** \class AliHLTTPCTrackArray
<pre>
//_____________________________________________________________
// AliHLTTPCTrackArray
//
// Track array class 
//
</pre>
*/

ClassImp(AliHLTTPCTrackArray)

AliHLTTPCTrackArray::AliHLTTPCTrackArray(){
  //Default constructor
  fSize = 0;
  fNTracks=0;
  fNAbsent=0;
  fTrackType='t';
  SetSize();
}


AliHLTTPCTrackArray::AliHLTTPCTrackArray(Int_t ntrack){
  //Constructor.
  fSize = 0;
  fNTracks=0;
  fNAbsent=0;
  fTrackType='t';
  SetSize(ntrack);
}

AliHLTTPCTrackArray::AliHLTTPCTrackArray(char* tracktype,Int_t ntrack){
  //Constructor.
  fSize = 0;
  fNTracks=0;
  fNAbsent=0;
  if(strcmp(tracktype,"AliHLTTPCTrack")==0) fTrackType='t';
  if(strcmp(tracktype,"AliHLTTPCConfMapTrack")==0) fTrackType='c';
  if(strcmp(tracktype,"AliHLTTPCHoughTrack")==0) fTrackType='h';
  if(strcmp(tracktype,"AliHLTTPCModelTrack")==0) fTrackType='m';
  SetSize(ntrack);
}

AliHLTTPCTrackArray::AliHLTTPCTrackArray(char* tracktype){
  //Constructor.
  fSize = 0;
  fNTracks=0;
  fNAbsent=0;
  if(strcmp(tracktype,"AliHLTTPCTrack")==0) fTrackType='t';
  if(strcmp(tracktype,"AliHLTTPCConfMapTrack")==0) fTrackType='c';
  if(strcmp(tracktype,"AliHLTTPCHoughTrack")==0) fTrackType='h';
  if(strcmp(tracktype,"AliHLTTPCModelTrack")==0) fTrackType='m';
  SetSize();
}

AliHLTTPCTrackArray::~AliHLTTPCTrackArray(){
  //Destructor
  DeleteArray();
}


AliHLTTPCTrack *AliHLTTPCTrackArray::NextTrack(){
  if(fNTracks<fSize) return fTrack[fNTracks++];
  SetSize(fSize+100);
   return fTrack[fNTracks++]; 
}

void AliHLTTPCTrackArray::DeleteArray(){
  for(Int_t i=0; i<fSize;i++)
    delete fTrack[i];
  delete[] fIsPresent;
  delete[] fTrack;
}

Bool_t AliHLTTPCTrackArray::SetSize(Int_t newsize){
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
      case 'h':
        for(Int_t i=0;i<fSize;i++){
          fTrack[i]   = new AliHLTTPCHoughTrack();
          fIsPresent[i] = kTRUE;
        }
        break;
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
    case 'h':
      for(Int_t i=fSize;i<newsize;i++){
        fTrack[i]   = new AliHLTTPCHoughTrack();
        fIsPresent[i] = kTRUE;
      }
      break;
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

void AliHLTTPCTrackArray::Reset(){
  fNTracks=0;
  fNAbsent=0;
  for(Int_t i=0; i<fSize;i++)
    fIsPresent[i] = kTRUE; 
}

void AliHLTTPCTrackArray::Remove(Int_t track){
  if(fIsPresent[track]){
    fIsPresent[track]=kFALSE;
    fNAbsent++;
  }
}

void AliHLTTPCTrackArray::FillTracks(Int_t ntracks, AliHLTTPCTrackSegmentData* tr){
  //Read tracks from shared memory (or memory)
  AliHLTTPCTrackSegmentData *trs = tr;
   for(Int_t i=0; i<ntracks; i++){
    AliHLTTPCTrack *track = NextTrack(); 
    track->SetPt(trs->fPt);
    track->SetPsi(trs->fPsi);
    track->SetTgl(trs->fTgl);
    track->SetNHits(trs->fNPoints);
    track->SetCharge(trs->fCharge);
    track->SetFirstPoint(trs->fX,trs->fY,trs->fZ);
    track->SetLastPoint(trs->fLastX,trs->fLastY,trs->fLastZ);
    track->SetHits( trs->fNPoints, trs->fPointIDs );
    UChar_t *tmpP = (UChar_t*)trs;
    tmpP += sizeof(AliHLTTPCTrackSegmentData)+trs->fNPoints*sizeof(UInt_t);
    trs = (AliHLTTPCTrackSegmentData*)tmpP;
  }
}

void AliHLTTPCTrackArray::FillTracks(Int_t ntracks, AliHLTTPCTrackSegmentData* tr,Int_t slice){
  //Read tracks from shared memory (or memory)
  AliHLTTPCTrackSegmentData *trs = tr;
  for(Int_t i=0; i<ntracks; i++){
    AliHLTTPCTrack *track = NextTrack(); 
    track->SetPt(trs->fPt);
    Float_t psi[1];
    psi[0]=trs->fPsi;
    AliHLTTPCTransform::Local2GlobalAngle(psi,slice);
    track->SetPsi(psi[0]);
    track->SetTgl(trs->fTgl);
    track->SetNHits(trs->fNPoints);
    track->SetCharge(trs->fCharge);
    Float_t first[3];
    first[0]=trs->fX;first[1]=trs->fY;first[2]=trs->fZ;
    AliHLTTPCTransform::Local2Global(first,slice);
    track->SetFirstPoint(first[0],first[1],first[2]);
    Float_t last[3];
    last[0]=trs->fLastX;last[1]=trs->fLastY;last[2]=trs->fLastZ;
    AliHLTTPCTransform::Local2Global(last,slice);
    track->SetLastPoint(last[0],last[1],last[2]);
    track->SetHits( trs->fNPoints, trs->fPointIDs );
    UChar_t *tmpP = (UChar_t*)trs;
    tmpP += sizeof(AliHLTTPCTrackSegmentData)+trs->fNPoints*sizeof(UInt_t);
    trs = (AliHLTTPCTrackSegmentData*)tmpP;
  }
}

UInt_t AliHLTTPCTrackArray::GetOutSize(){
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

UInt_t AliHLTTPCTrackArray::WriteTracks(UInt_t & ntracks,AliHLTTPCTrackSegmentData* tr){
  ntracks = GetOutCount();
  return WriteTracks(tr);
}

UInt_t AliHLTTPCTrackArray::WriteTracks(AliHLTTPCTrackSegmentData* tr){
  if(GetTrackType()=='c') return WriteConfMapTracks(tr);
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
    tmpP += sizeof(AliHLTTPCTrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    size += sizeof(AliHLTTPCTrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    tP = (AliHLTTPCTrackSegmentData*)tmpP;
  }
  return size;
}

UInt_t AliHLTTPCTrackArray::WriteConfMapTracks(AliHLTTPCTrackSegmentData* tr){
  // use first and last point objects
  AliHLTTPCTrackSegmentData *tP = tr;
  UInt_t *pP;
  UInt_t size = 0;
  for(Int_t i=0; i<fNTracks; i++){  //loop over all tracks
    AliHLTTPCConfMapTrack *track =(AliHLTTPCConfMapTrack *) GetCheckedTrack(i); //use only present tracks
    if(!track) continue;                           //use only present tracks
    AliHLTTPCConfMapPoint *hit = (AliHLTTPCConfMapPoint*)track->lastHit;
    AliHLTTPCConfMapPoint *lastHit = (AliHLTTPCConfMapPoint*)track->firstHit;
/*
    tP->fX = hit->GetX();
    tP->fY = hit->GetY();
    tP->fZ = hit->GetZ();
    tP->fLastX = lastHit->GetX();
    tP->fLastY = lastHit->GetY();
    tP->fLastZ = lastHit->GetZ();
*/   
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
    tmpP += sizeof(AliHLTTPCTrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    size +=sizeof(AliHLTTPCTrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    tP = (AliHLTTPCTrackSegmentData*)tmpP;
  }
  return size;
}

void AliHLTTPCTrackArray::AddLast(AliHLTTPCTrack *track)
{
  AliHLTTPCTrack *tpt = NextTrack();
  tpt->Set(track);
  
}

void AliHLTTPCTrackArray::AddTracks(AliHLTTPCTrackArray *newtrack,Bool_t remove_old,Int_t slice)
{
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
    track->Set(tpt);
    if(slice>=0)
      track->Rotate(slice); //Rotate track to global coordinates
    /*
      AliHLTTPCTrack *track;
      if(GetTrackType()=='h')
      track = (AliHLTTPCHoughTrack*)NextTrack();
      else
      track = NextTrack();
      track->Set(tpt);
    */
  }
}


void AliHLTTPCTrackArray::Compress(){
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

void AliHLTTPCTrackArray::QSort(){
  // compress an sort
  Compress();
  QSort(fTrack,0,fNTracks);
}

void AliHLTTPCTrackArray::QSort( AliHLTTPCTrack **a, Int_t first, Int_t last){

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

Int_t AliHLTTPCTrackArray::TrackCompare(AliHLTTPCTrack *a, AliHLTTPCTrack *b){
   // Compare the two tracks.
  
  return b->Compare(a);
  
  /*
    if(fTrackType=='h')
    {
    AliHLTTPCHoughTrack *tra = (AliHLTTPCHoughTrack*)a;
    AliHLTTPCHoughTrack *trb = (AliHLTTPCHoughTrack*)b;
    if(tra->GetWeight() < trb->GetWeight()) return 1;
    if(tra->GetWeight() > trb->GetWeight()) return -1;
    }
    else
    {
    if(a->GetNHits() < b->GetNHits()) return 1;
    if(a->GetNHits() > b->GetNHits()) return -1;
    }
    
    return 0;
  */
}


