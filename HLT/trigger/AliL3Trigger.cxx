
#include "AliL3Trigger.h"
#include "AliL3TrackArray.h"
#include "AliL3Transform.h"
#include "AliL3Vertex.h"
#include "AliL3Defs.h"
#include "AliL3DigitData.h"
#include "AliL3Logging.h"
#include "AliL3Track.h"
#include "AliL3MemHandler.h"

ClassImp(AliL3Trigger)

AliL3Trigger::AliL3Trigger()
{
  fDataSize=0;
  fTracks=0;
  fTransform=0;
  fDigitRowData=0;
  fOutput=0;
  fVertex=0;
}

AliL3Trigger::~AliL3Trigger()
{
  if(fTracks)
    delete fTracks;
  if(fTransform)
    delete fTransform;
}

void AliL3Trigger::InitTrigger()
{
  if(fTracks)
    delete fTracks;
  if(fTransform)
    delete fTransform;
  fTracks = new AliL3TrackArray();
  fTransform = new AliL3Transform;
}

void AliL3Trigger::InitPatch(Int_t slice,Int_t patch)
{
  fSlice=slice;
  fPatch=patch;
  fTracks->Reset();
}

void AliL3Trigger::FillTracks(Int_t ntracks,AliL3TrackSegmentData *tr)
{
  fTracks->FillTracks(ntracks,tr);
}

void AliL3Trigger::FillData(AliL3DigitRowData *data)
{
  fDigitRowData = data;
}

void AliL3Trigger::SetParameters(Float_t zcut,Int_t timematch,Int_t padmatch)
{
  fZcut=zcut;
  fTimeMatch=timematch;
  fPadMatch=padmatch;
}

void AliL3Trigger::SetOutputData(AliL3DigitRowData *ptr)
{
  fOutput=ptr;
}

void AliL3Trigger::RemovePileupTracks()
{
  Double_t xc,yc,zc;
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliL3Track *track = fTracks->GetCheckedTrack(i);
      if(!track) continue;
      track->Rotate(fSlice,kTRUE);
      track->CalculateHelix();
      track->GetClosestPoint(fVertex,xc,yc,zc);
      if(fabs(zc) > fZcut)
	{
	  fTracks->Remove(i);
	  continue;
	}
    }
  fTracks->Compress();
}

void AliL3Trigger::RemovePileupData()
{
  Float_t hit[3];
  Int_t sector,row;
  struct rowhit {Int_t pad; Int_t time;};
  rowhit row_cross[(fTracks->GetNTracks())];
  Int_t digitcount[(NumRows[fPatch])];
  Int_t totalcount=0;
  AliL3DigitRowData *rowPt = fDigitRowData;
  for(Int_t i=NRows[fPatch][0]; i<=NRows[fPatch][1]; i++)
    {
      digitcount[(i-NRows[fPatch][0])]=0;
      for(Int_t j=0; j<fTracks->GetNTracks(); j++)
	{
	  AliL3Track *track = fTracks->GetCheckedTrack(j);
	  if(!track) continue;
	  track->GetCrossingPoint(i,hit);
	  fTransform->Slice2Sector(fSlice,i,sector,row);
	  fTransform->Local2Raw(hit,sector,row);
	  row_cross[j].pad = (Int_t)rint(hit[1]);
	  row_cross[j].time = (Int_t)rint(hit[2]);
	}
      AliL3DigitData *digPt = (AliL3DigitData*)rowPt->fDigitData;
      Bool_t mark;
      for(Int_t k=0; k<rowPt->fNDigit; k++)
	{
	  mark = kFALSE;
	  for(Int_t l=0; l<fTracks->GetNTracks(); l++)
	    {
	      if(abs((Int_t)digPt[k].fPad-row_cross[l].pad) < fPadMatch &&
		 abs((Int_t)digPt[k].fTime-row_cross[l].time) < fTimeMatch)
		{
		  digitcount[(i-NRows[fPatch][0])]++;
		  totalcount++;
		  mark=kTRUE;
		  break;
		}
	    }
	  if(mark==kTRUE)
	    digPt[k].fCharge=0;
	}
      AliL3MemHandler::UpdateRowPointer(rowPt);
    }
  
  Int_t size = totalcount*sizeof(AliL3DigitData) + NumRows[fPatch]*sizeof(AliL3DigitRowData);
  fDataSize = size;
  LOG(AliL3Log::kDebug,"AliL3Trigger::RemovePileupData","Memory")
    <<"Allocating "<<size<<" bytes of data for trigger event"<<ENDLOG;
  Byte_t *data = new Byte_t[size];
  memset(data,0,size);
  AliL3DigitRowData *tempPt = (AliL3DigitRowData*)data;
  rowPt = fDigitRowData;
  
  Int_t localcount;
  for(Int_t i=NRows[fPatch][0]; i<=NRows[fPatch][1]; i++)
    {
      tempPt->fRow = i;
      tempPt->fNDigit = digitcount[(i-NRows[fPatch][0])];
      AliL3DigitData *digPt = (AliL3DigitData*)rowPt->fDigitData;
      localcount=0;
      for(Int_t j=0; j<rowPt->fNDigit; j++)
	{
	  if(digPt[j].fCharge==0) continue;
	  if(localcount >= digitcount[(i-NRows[fPatch][0])])
	    {
	      LOG(AliL3Log::kFatal,"AliL§Trigger::RemovePileupData","Array")
		<<"Mismatch in digitcount: "<<localcount<<" "<<digitcount[(i-NRows[fPatch][0])]<<ENDLOG;
	      return;
	    }
	  tempPt->fDigitData[localcount].fCharge=digPt[j].fCharge;
	  tempPt->fDigitData[localcount].fPad=digPt[j].fPad;
	  tempPt->fDigitData[localcount].fTime=digPt[j].fTime;
	  localcount++;
	}
      if(digitcount[(i-NRows[fPatch][0])]!=localcount)
	{
	  LOG(AliL3Log::kFatal,"AliL§Trigger::RemovePileupData","Array")
	    <<"Mismatch in digitcount: "<<localcount<<" "<<digitcount[(i-NRows[fPatch][0])]<<ENDLOG;
	}
      AliL3MemHandler::UpdateRowPointer(rowPt);
      Byte_t *tmp = (Byte_t*)tempPt;
      Int_t size = sizeof(AliL3DigitRowData)+digitcount[(i-NRows[fPatch][0])]*sizeof(AliL3DigitData);
      tmp += size;
      tempPt = (AliL3DigitRowData*)tmp;
    }
  
  fOutput=(AliL3DigitRowData*)data;
}


