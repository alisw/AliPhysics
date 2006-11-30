
#include "AliHLTTrigger.h"
#include "AliHLTTrackArray.h"
#include "AliHLTTransform.h"
#include "AliHLTVertex.h"
#include "AliHLTDefs.h"
#include "AliHLTDigitData.h"
#include "AliHLTLogging.h"
#include "AliHLTTrack.h"
#include "AliHLTMemHandler.h"

ClassImp(AliHLTTrigger)

AliHLTTrigger::AliHLTTrigger()
{
  fDataSize=0;
  fTracks=0;
  fDigitRowData=0;
  fOutput=0;
  fVertex=0;
}

AliHLTTrigger::~AliHLTTrigger()
{
  if(fTracks)
    delete fTracks;
}

void AliHLTTrigger::InitTrigger()
{
  if(fTracks)
    delete fTracks;
  fTracks = new AliHLTTrackArray();
}

void AliHLTTrigger::InitPatch(Int_t slice,Int_t patch)
{
  fSlice=slice;
  fPatch=patch;
  fTracks->Reset();
}

void AliHLTTrigger::FillTracks(Int_t ntracks,AliHLTTrackSegmentData *tr)
{
  fTracks->FillTracks(ntracks,tr);
}

void AliHLTTrigger::FillData(AliHLTDigitRowData *data)
{
  fDigitRowData = data;
}

void AliHLTTrigger::SetParameters(Float_t zcut,Int_t timematch,Int_t padmatch)
{
  fZcut=zcut;
  fTimeMatch=timematch;
  fPadMatch=padmatch;
}

void AliHLTTrigger::SetOutputData(AliHLTDigitRowData *ptr)
{
  fOutput=ptr;
}

void AliHLTTrigger::RemovePileupTracks()
{
  Double_t xc,yc,zc;
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliHLTTrack *track = fTracks->GetCheckedTrack(i);
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

void AliHLTTrigger::RemovePileupData()
{
  Float_t hit[3];
  Int_t sector,row;
  struct rowhit {Int_t pad; Int_t time;};
  rowhit row_cross[(fTracks->GetNTracks())];
  Int_t digitcount[(NumRows[fPatch])];
  Int_t totalcount=0;
  AliHLTDigitRowData *rowPt = fDigitRowData;
  for(Int_t i=NRows[fPatch][0]; i<=NRows[fPatch][1]; i++)
    {
      digitcount[(i-NRows[fPatch][0])]=0;
      for(Int_t j=0; j<fTracks->GetNTracks(); j++)
	{
	  AliHLTTrack *track = fTracks->GetCheckedTrack(j);
	  if(!track) continue;
	  track->GetCrossingPoint(i,hit);
	  AliHLTTransform::Slice2Sector(fSlice,i,sector,row);
	  AliHLTTransform::Local2Raw(hit,sector,row);
	  row_cross[j].pad = (Int_t)rint(hit[1]);
	  row_cross[j].time = (Int_t)rint(hit[2]);
	}
      AliHLTDigitData *digPt = (AliHLTDigitData*)rowPt->fDigitData;
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
      AliHLTMemHandler::UpdateRowPointer(rowPt);
    }
  
  Int_t size = totalcount*sizeof(AliHLTDigitData) + NumRows[fPatch]*sizeof(AliHLTDigitRowData);
  fDataSize = size;
  LOG(AliHLTLog::kDebug,"AliHLTTrigger::RemovePileupData","Memory")
    <<"Allocating "<<size<<" bytes of data for trigger event"<<ENDLOG;
  Byte_t *data = new Byte_t[size];
  memset(data,0,size);
  AliHLTDigitRowData *tempPt = (AliHLTDigitRowData*)data;
  rowPt = fDigitRowData;
  
  Int_t localcount;
  for(Int_t i=NRows[fPatch][0]; i<=NRows[fPatch][1]; i++)
    {
      tempPt->fRow = i;
      tempPt->fNDigit = digitcount[(i-NRows[fPatch][0])];
      AliHLTDigitData *digPt = (AliHLTDigitData*)rowPt->fDigitData;
      localcount=0;
      for(Int_t j=0; j<rowPt->fNDigit; j++)
	{
	  if(digPt[j].fCharge==0) continue;
	  if(localcount >= digitcount[(i-NRows[fPatch][0])])
	    {
	      LOG(AliHLTLog::kFatal,"AliL§Trigger::RemovePileupData","Array")
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
	  LOG(AliHLTLog::kFatal,"AliL§Trigger::RemovePileupData","Array")
	    <<"Mismatch in digitcount: "<<localcount<<" "<<digitcount[(i-NRows[fPatch][0])]<<ENDLOG;
	}
      AliHLTMemHandler::UpdateRowPointer(rowPt);
      Byte_t *tmp = (Byte_t*)tempPt;
      Int_t size = sizeof(AliHLTDigitRowData)+digitcount[(i-NRows[fPatch][0])]*sizeof(AliHLTDigitData);
      tmp += size;
      tempPt = (AliHLTDigitRowData*)tmp;
    }
  
  fOutput=(AliHLTDigitRowData*)data;
}


