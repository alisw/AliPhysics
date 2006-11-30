// @(#) $Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTTrackArray.h"
#include "AliHLTHoughTrack.h"
#include "AliHLTModelTrack.h"
#include "AliHLTConfMapTrack.h"
#include "AliHLTTrackSegmentData.h"
#include "AliHLTTransform.h"
#include "AliHLTConfMapPoint.h"

/** \class AliHLTTrackArray
<pre>
//_____________________________________________________________
// AliHLTTrackArray
//
// Track array class 
//
</pre>
*/

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTrackArray)

AliHLTTrackArray::AliHLTTrackArray()
{
  //Default constructor
  fSize = 0;
  fNTracks=0;
  fNAbsent=0;
  fTrackType='t';
  SetSize();
}


AliHLTTrackArray::AliHLTTrackArray(Int_t ntrack)
{
  //Constructor.
  fSize = 0;
  fNTracks=0;
  fNAbsent=0;
  fTrackType='t';
  SetSize(ntrack);
}

AliHLTTrackArray::AliHLTTrackArray(char* tracktype,Int_t ntrack)
{
  //Constructor.
  fSize = 0;
  fNTracks=0;
  fNAbsent=0;
  if(strcmp(tracktype,"AliHLTTrack")==0) fTrackType='t';
  if(strcmp(tracktype,"AliHLTConfMapTrack")==0) fTrackType='c';
  if(strcmp(tracktype,"AliHLTHoughTrack")==0) fTrackType='h';
  if(strcmp(tracktype,"AliHLTModelTrack")==0) fTrackType='m';
  SetSize(ntrack);
}

AliHLTTrackArray::AliHLTTrackArray(char* tracktype)
{
  //Constructor.
  fSize = 0;
  fNTracks=0;
  fNAbsent=0;
  if(strcmp(tracktype,"AliHLTTrack")==0) fTrackType='t';
  if(strcmp(tracktype,"AliHLTConfMapTrack")==0) fTrackType='c';
  if(strcmp(tracktype,"AliHLTHoughTrack")==0) fTrackType='h';
  if(strcmp(tracktype,"AliHLTModelTrack")==0) fTrackType='m';
  SetSize();
}

AliHLTTrackArray::~AliHLTTrackArray()
{
  //Destructor
  DeleteArray();
}


AliHLTTrack *AliHLTTrackArray::NextTrack()
{
  //next track in array
  if(fNTracks<fSize) return fTrack[fNTracks++];
  SetSize(fSize+100);
   return fTrack[fNTracks++]; 
}

void AliHLTTrackArray::DeleteArray()
{
  //delete array
  for(Int_t i=0; i<fSize;i++)
    delete fTrack[i];
  delete[] fIsPresent;
  delete[] fTrack;
}

Bool_t AliHLTTrackArray::SetSize(Int_t newsize)
{
  //set size
  if(newsize<=fSize) return kFALSE; //shrink comes later!! 
  if(!fSize){
    fSize = newsize;
    fTrack = new AliHLTTrack*[fSize];
    fIsPresent = new Bool_t[fSize];
    switch(fTrackType){
      case 't':
        for(Int_t i=0;i<fSize;i++){
          fTrack[i]   = new AliHLTTrack();
          fIsPresent[i] = kTRUE;
        }
        break;
      case 'c':  
        for(Int_t i=0;i<fSize;i++){
          fTrack[i]   = new AliHLTConfMapTrack();
          fIsPresent[i] = kTRUE;
        }
        break;
      case 'h':
        for(Int_t i=0;i<fSize;i++){
          fTrack[i]   = new AliHLTHoughTrack();
          fIsPresent[i] = kTRUE;
        }
        break;
       case 'm':
        for(Int_t i=0;i<fSize;i++){
          fTrack[i]   = new AliHLTModelTrack();
          fIsPresent[i] = kTRUE;
        }
        break;
      default: 
        return kFALSE;
    }
    return kTRUE;
  }
  AliHLTTrack **tmp = new AliHLTTrack*[fSize];
  Bool_t *pre = new Bool_t[fSize];
  for(Int_t i=0; i<fSize;i++){
    tmp[i] = fTrack[i];
    pre[i] = fIsPresent[i];
  }
  delete[]  fTrack;
  delete[] fIsPresent;
  fTrack =  new AliHLTTrack*[newsize];
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
        fTrack[i]   = new AliHLTTrack();
        fIsPresent[i] = kTRUE;
      }
      break;
    case 'c':  
      for(Int_t i=fSize;i<newsize;i++){
        fTrack[i]   = new AliHLTConfMapTrack();
        fIsPresent[i] = kTRUE;
      }
      break;
    case 'h':
      for(Int_t i=fSize;i<newsize;i++){
        fTrack[i]   = new AliHLTHoughTrack();
        fIsPresent[i] = kTRUE;
      }
      break;
    case 'm':
      for(Int_t i=fSize;i<newsize;i++){
        fTrack[i]   = new AliHLTModelTrack();
        fIsPresent[i] = kTRUE;
      }
      break;
    default: 
      return kFALSE;
  }
  fSize = newsize;
  return kTRUE;
}

void AliHLTTrackArray::Reset()
{
  //reset
  fNTracks=0;
  fNAbsent=0;
  for(Int_t i=0; i<fSize;i++)
    fIsPresent[i] = kTRUE; 
}

void AliHLTTrackArray::Remove(Int_t track)
{
  //remove track
  if(fIsPresent[track]){
    fIsPresent[track]=kFALSE;
    fNAbsent++;
  }
}

void AliHLTTrackArray::FillTracks(Int_t ntracks, AliHLTTrackSegmentData* tr){
  //Read tracks from shared memory (or memory)
  AliHLTTrackSegmentData *trs = tr;
   for(Int_t i=0; i<ntracks; i++){
    AliHLTTrack *track = NextTrack(); 
    track->SetPt(trs->fPt);
    track->SetPsi(trs->fPsi);
    track->SetTgl(trs->fTgl);
    track->SetPterr(trs->fPterr);
    track->SetPsierr(trs->fPsierr);
    track->SetTglerr(trs->fTglerr);
    track->SetNHits(trs->fNPoints);
    track->SetCharge(trs->fCharge);
    track->SetFirstPoint(trs->fX,trs->fY,trs->fZ);
    track->SetLastPoint(trs->fLastX,trs->fLastY,trs->fLastZ);
    track->SetHits( trs->fNPoints, trs->fPointIDs );
#ifdef ROWHOUGHPARAMS
    if(GetTrackType()=='h') {
      ((AliHLTHoughTrack *)track)->SetWeight(trs->fWeight);
      ((AliHLTHoughTrack *)track)->SetBinXY(trs->fBinX,trs->fBinY,trs->fBinXSize,trs->fBinYSize);
    }
    track->SetMCid(trs->fTrackID);
    track->SetRowRange(trs->fRowRange1,trs->fRowRange2);
    track->SetSector(trs->fSector);
    track->SetPID(trs->fPID);
#endif
    UChar_t *tmpP = (UChar_t*)trs;
    tmpP += sizeof(AliHLTTrackSegmentData)+trs->fNPoints*sizeof(UInt_t);
    trs = (AliHLTTrackSegmentData*)tmpP;
  }
}

void AliHLTTrackArray::FillTracks(Int_t ntracks, AliHLTTrackSegmentData* tr,Int_t slice)
{
  //Read tracks from shared memory (or memory)
  AliHLTTrackSegmentData *trs = tr;
  for(Int_t i=0; i<ntracks; i++){
    AliHLTTrack *track = NextTrack(); 
    track->SetPt(trs->fPt);
    track->SetPterr(trs->fPterr);
    Float_t psi[1];
    psi[0]=trs->fPsi;
    AliHLTTransform::Local2GlobalAngle(psi,slice);
    track->SetPsi(psi[0]);
    track->SetTgl(trs->fTgl);
    track->SetPsierr(trs->fPsierr);
    track->SetTglerr(trs->fTglerr);
    track->SetNHits(trs->fNPoints);
    track->SetCharge(trs->fCharge);
    Float_t first[3];
    first[0]=trs->fX;first[1]=trs->fY;first[2]=trs->fZ;
    AliHLTTransform::Local2Global(first,slice);
    track->SetFirstPoint(first[0],first[1],first[2]);
    Float_t last[3];
    last[0]=trs->fLastX;last[1]=trs->fLastY;last[2]=trs->fLastZ;
    AliHLTTransform::Local2Global(last,slice);
    track->SetLastPoint(last[0],last[1],last[2]);
    track->SetHits( trs->fNPoints, trs->fPointIDs );
#ifdef ROWHOUGHPARAMS
    if(GetTrackType()=='h') {
      ((AliHLTHoughTrack *)track)->SetWeight(trs->fWeight);
      ((AliHLTHoughTrack *)track)->SetBinXY(trs->fBinX,trs->fBinY,trs->fBinXSize,trs->fBinYSize);
    }
    track->SetMCid(trs->fTrackID);
    track->SetRowRange(trs->fRowRange1,trs->fRowRange2);
    track->SetSector(slice);
    track->SetPID(trs->fPID);
#endif
    UChar_t *tmpP = (UChar_t*)trs;
    tmpP += sizeof(AliHLTTrackSegmentData)+trs->fNPoints*sizeof(UInt_t);
    trs = (AliHLTTrackSegmentData*)tmpP;
  }
}

UInt_t AliHLTTrackArray::GetOutSize()
{
  //get size for IO
  UInt_t count = GetOutCount();   //use only present tracks
  UInt_t tHits = 0;
  for(Int_t i=0;i<fNTracks;i++){  //loop over all tracks
    AliHLTTrack *track = GetCheckedTrack(i);  //use only present tracks
    if(track)                                       //use only present tracks
      tHits += track->GetNHits();
  }

  //calculate size of track
  return count*sizeof(AliHLTTrackSegmentData)+sizeof(UInt_t)*tHits;
}

UInt_t AliHLTTrackArray::WriteTracks(UInt_t & ntracks,AliHLTTrackSegmentData* tr)
{
  //write tracks
  ntracks = GetOutCount();
  return WriteTracks(tr);
}

UInt_t AliHLTTrackArray::WriteTracks(AliHLTTrackSegmentData* tr)
{
  //if(GetTrackType()=='c') return WriteConfMapTracks(tr);
  AliHLTTrackSegmentData *tP = tr;
  UInt_t *pP;
  UInt_t size = 0;
  for(Int_t i=0; i<fNTracks; i++){  //loop over all tracks
    AliHLTTrack *track = GetCheckedTrack(i); //use only present tracks
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
#ifdef ROWHOUGHPARAMS
    if(GetTrackType()=='h') {
      tP->fWeight = ((AliHLTHoughTrack *)track)->GetWeight();
      tP->fBinX = ((AliHLTHoughTrack *)track)->GetBinX();
      tP->fBinY = ((AliHLTHoughTrack *)track)->GetBinY();
      tP->fBinXSize = ((AliHLTHoughTrack *)track)->GetSizeX();
      tP->fBinYSize = ((AliHLTHoughTrack *)track)->GetSizeY();
    }
    tP->fTrackID = track->GetMCid();
    tP->fRowRange1 = track->GetFirstRow();
    tP->fRowRange2 = track->GetLastRow();
    tP->fSector = track->GetSector();
    tP->fPID = track->GetPID();
#endif
    pP = (UInt_t*)track->GetHitNumbers();
    for (UInt_t j=0;j<tP->fNPoints;j++){
      tP->fPointIDs[j] = pP[j];
    }
    Byte_t *tmpP = (Byte_t *)tP;
    tmpP += sizeof(AliHLTTrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    size += sizeof(AliHLTTrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    tP = (AliHLTTrackSegmentData*)tmpP;
  }
  return size;
}

UInt_t AliHLTTrackArray::WriteConfMapTracks(AliHLTTrackSegmentData* tr)
{
  // use first and last point objects
  AliHLTTrackSegmentData *tP = tr;
  UInt_t *pP;
  UInt_t size = 0;
  for(Int_t i=0; i<fNTracks; i++){  //loop over all tracks
    AliHLTConfMapTrack *track =(AliHLTConfMapTrack *) GetCheckedTrack(i); //use only present tracks
    if(!track) continue;                           //use only present tracks
    AliHLTConfMapPoint *hit = (AliHLTConfMapPoint*)track->GetLastHit();
    AliHLTConfMapPoint *lastHit = (AliHLTConfMapPoint*)track->GetFirstHit();
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
#ifdef ROWHOUGHPARAMS
    tP->fTrackID = track->GetMCid();
    tP->fRowRange1 = track->GetFirstRow();
    tP->fRowRange2 = track->GetLastRow();
    tP->fSector = track->GetSector();
    tP->fPID = track->GetPID();
#endif
    tP->fNPoints = track->GetNHits();
    pP = (UInt_t*)track->GetHitNumbers();
    for (UInt_t j=0;j<tP->fNPoints;j++){
      tP->fPointIDs[j] = pP[j];
    }
    Byte_t *tmpP = (Byte_t *)tP;
    tmpP += sizeof(AliHLTTrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    size +=sizeof(AliHLTTrackSegmentData)+tP->fNPoints*sizeof(UInt_t);
    tP = (AliHLTTrackSegmentData*)tmpP;
  }
  return size;
}

void AliHLTTrackArray::AddLast(AliHLTTrack *track)
{
  //add track to last position
  AliHLTTrack *tpt = NextTrack();
  tpt->Set(track);
  
}

void AliHLTTrackArray::AddTracks(AliHLTTrackArray *newtrack,Bool_t remove_old,Int_t slice)
{
  //add tracks
  if(GetTrackType() != newtrack->GetTrackType() && GetTrackType()!='t')
    {
      LOG(AliHLTLog::kError,"AliHLTTrackArray::AddTracks","Track types")
	<<"Bad idea to add tracks of different types"<<ENDLOG;
      return;
    }
  if(fSize < fNTracks+newtrack->GetNPresent())
    SetSize(fSize+newtrack->GetSize());
  for(Int_t i =0;i<newtrack->GetNTracks();i++){
    AliHLTTrack *tpt = newtrack->GetCheckedTrack(i);
    if(!tpt) continue;
    if(remove_old)
      newtrack->Remove(i);
    AliHLTTrack *track = NextTrack();
    track->Set(tpt);
    if(slice>=0)
      track->Rotate(slice); //Rotate track to global coordinates
    /*
      AliHLTTrack *track;
      if(GetTrackType()=='h')
      track = (AliHLTHoughTrack*)NextTrack();
      else
      track = NextTrack();
      track->Set(tpt);
    */
  }
}

void AliHLTTrackArray::Compress()
{
  //compress array
  if(GetNPresent()==GetNTracks()) return;
  AliHLTTrack **tmp =  new AliHLTTrack *[fNTracks];
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

void AliHLTTrackArray::QSort()
{
  // compress and sort
  Compress();
  QSort(fTrack,0,fNTracks);
}

void AliHLTTrackArray::QSort( AliHLTTrack **a, Int_t first, Int_t last)
{
   // Sort array of AliHLTTrack pointers using a quicksort algorithm.
   // Uses TrackCompare() to compare objects.
   // Thanks to Root! 

   static AliHLTTrack *tmp;
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

Int_t AliHLTTrackArray::TrackCompare(AliHLTTrack *a, AliHLTTrack *b) const
{
   // Compare the two tracks.
  
  return b->Compare(a);
  
  /*
    if(fTrackType=='h')
    {
    AliHLTHoughTrack *tra = (AliHLTHoughTrack*)a;
    AliHLTHoughTrack *trb = (AliHLTHoughTrack*)b;
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


