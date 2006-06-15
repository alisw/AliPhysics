// $Header$

#include "TPCLoader.h"
#include "TPCData.h"
#include <Alieve/TPCSector2D.h>
#include <Alieve/TPCSector3D.h>
#include <Reve/RGTopFrame.h>

#include <AliRawReaderRoot.h>
#include <AliTPCRawStream.h>

#include <TSystem.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// TPCLoader
//

ClassImp(TPCLoader)

TPCLoader::TPCLoader(const Text_t* n, const Text_t* t) :
  RenderElementList(n, t),

  fEvent(-1),
  fDoubleSR(kFALSE),

  fReader(0),
  fData(0),

  fSec2Ds(36),
  fSec3Ds(36),

  fSetInitSectorParams(kFALSE)
{
  SetMainColorPtr(0);
}

TPCLoader::~TPCLoader()
{
  delete fReader;
  delete fData;
}

/**************************************************************************/


void TPCLoader::SetData(TPCData* d)
{
  delete fData;
  fData = d;
}

/**************************************************************************/

void TPCLoader::OpenFile()
{
  static const Exc_t eH("TPCLoader::OpenFile ");

  if(gSystem->AccessPathName(fFile, kReadPermission))
      throw(eH + "can not read '" + fFile + "'.");

  if(fData == 0)
    fData = new TPCData;

  delete fReader;
  fReader =  0;
  fEvent  = -1;

  fReader = new AliRawReaderRoot(fFile);
  fReader->LoadEquipmentIdsMap
    (gSystem->ExpandPathName("$(ALICE_ROOT)/TPC/mapping/EquipmentIdMap.data"));

  NextEvent();
  LoadEvent();
  UpdateSectors();
}

void TPCLoader::LoadEvent()
{
  static const Exc_t eH("TPCLoader::LoadEvent ");

  if(fReader == 0)
    throw(eH + "data file not opened.");

  printf("Now loading event %d\n", fEvent);
  fReader->Reset();
  AliTPCRawStream input(fReader);
  input.SetOldRCUFormat(kTRUE);
  fReader->Select(0);

  fData->DropAllSectors();
  fData->LoadRaw(input, kTRUE, kTRUE);  
}

void TPCLoader::NextEvent(Bool_t rewindOnEnd)
{
  static const Exc_t eH("TPCLoader::NextEvent ");

  if(fReader == 0)
    throw(eH + "data file not opened.");

  if(fReader->NextEvent() == kTRUE) {
    ++fEvent;
  } else {
    if(fEvent == -1)
      throw(eH + "no events available.");
    if(rewindOnEnd) {
      printf("Reached end of stream (event=%d), rewinding to first event.\n", fEvent);
      fReader->RewindEvents();
      fReader->NextEvent();
      fEvent = 0;
    } else {
      throw(eH + "last event reached.");
    }
  }
}

void TPCLoader::GotoEvent(Int_t event)
{
  static const Exc_t eH("TPCLoader::GotoEvent ");

  if(fReader == 0)
    throw(eH + "data file not opened.");

  if(event == fEvent)
    return;
  Bool_t checkEnd;
  if(event < fEvent) {
    fReader->RewindEvents();
    fEvent = -1;
    checkEnd = kFALSE;
  } else {
    checkEnd = kTRUE;
  }
  do {
    NextEvent();
  } while(fEvent != event && !(checkEnd == kTRUE && fEvent == 0));
  LoadEvent();
  UpdateSectors();
}

/**************************************************************************/

void TPCLoader::UpdateSectors()
{
  gReve->DisableRedraw();
  for(Int_t i=0; i<=35; ++i) {
    if(fSec2Ds[i] != 0) {
      fSec2Ds[i]->IncRTS();
    } else {
      TPCSectorData* sd = fData->GetSectorData(i);
      if(sd != 0) {
	TPCSector2D* s = new TPCSector2D(Form("Sector2D %d", i));
	fSec2Ds[i] = s;
	s->SetSectorID(i);
	s->SetDataSource(fData);

	if(fDoubleSR)
	  s->SetMaxTime(1023);

        if(fSetInitSectorParams) {
          s->SetMinTime(fInitMinTime);
          s->SetMaxTime(fInitMaxTime);
          s->SetThreshold(fInitThreshold);
        }

	s->SetTrans(kTRUE);
	s->SetFrameColor(36);

	gReve->AddRenderElement(s);
	gReve->DrawRenderElement(s);
      }
    }

    if(fSec3Ds[i] != 0) {
      fSec3Ds[i]->IncRTS();
    }
  }
  gReve->EnableRedraw();
}

void TPCLoader::CreateSectors3D()
{
  gReve->DisableRedraw();
  for(Int_t i=0; i<=35; ++i) {
    TPCSectorData* sd = fData->GetSectorData(i);
    if(sd != 0 && fSec3Ds[i] == 0) {
      TPCSector3D* s = new TPCSector3D(Form("Sector3D %d", i));
      fSec3Ds[i] = s;
      s->SetSectorID(i);
      s->SetDataSource(fData);

      if(fDoubleSR)
	s->SetDriftVel(2.273);
      if(fSec2Ds[i] != 0)
	s->CopyVizParams(*fSec2Ds[i]);

      s->SetTrans(kTRUE);
      s->SetFrameColor(36);

      gReve->AddRenderElement(s);
      gReve->DrawRenderElement(s);
    }
  }
  gReve->EnableRedraw();
}

void TPCLoader::DeleteSectors3D()
{
  gReve->DisableRedraw();
  for(Int_t i=0; i<=35; ++i) {
    RenderElement* re = fSec3Ds[i];
    if(re != 0) {
      gReve->UndrawRenderElement(re);
      delete re;
      fSec3Ds[i] = 0;
    }
  }
  gReve->EnableRedraw();
}

/**************************************************************************/

void TPCLoader::SetInitParams(Int_t mint, Int_t maxt, Int_t thr)
{
  fSetInitSectorParams = kTRUE;
  fInitMinTime   = mint;
  fInitMaxTime   = maxt;
  fInitThreshold = thr;
}
