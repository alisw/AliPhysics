// $Header$

// Functions to read rootified raw-data from TPC sector test.
//
// Use tpc_sector_raw_test("filename.root") for initialization,
// next_event() to advance along the data stream.
// When there is no more data ROOT will crash.

class AliRawReaderRoot;

namespace Alieve {
class TPCData;
class TPCSector2D;
class TPCSector3D;
}

using namespace Alieve;

TPCData*     x = 0;
TPCSector2D* s = 0;
TPCSector3D* t = 0;

AliRawReaderRoot* reader =  0;
Int_t             event  = -1;
Int_t     default_sector =  13;

void tpc_sector_raw_test(const char *file = "", Int_t ievent = 0)
{
  gStyle->SetPalette(1, 0);

  reader = new AliRawReaderRoot(file);
  //reader->LoadEquipmentIdsMap
  //  (gSystem->ExpandPathName("$(ALICE_ROOT)/TPC/mapping/EquipmentIdMap.data"));
  //  (gSystem->ExpandPathName("EquipmentIdMap.data"));

  reader->Reset();
  for(Int_t i=0; i<ievent; ++i, ++event) {
    if(reader->NextEvent() == kFALSE) {
      printf("End of raw stream at event %d (reqired event %d).\n", i, ievent);
      return;
    }
  }

  x = new TPCData;
  // x->SetLoadPedestal(5);
  x->SetLoadThreshold(5);
  x->SetAutoPedestal(kTRUE);

  s = new TPCSector2D();
  s->SetSectorID(default_sector);
  s->SetAutoTrans(kTRUE); // place on proper 3D coordinates
  s->SetDataSource(x);
  s->SetFrameColor(36);
  gReve->AddRenderElement(s);
  gReve->DrawRenderElement(s);

  t = new TPCSector3D();
  t->SetSectorID(default_sector);
  t->SetAutoTrans(kTRUE);
  t->SetDataSource(x);
  t->SetMaxTime(1023);
  t->SetDriftVel(2.273);
  gReve->AddRenderElement(t);
  gReve->DrawRenderElement(t);

  next_event();
}

void next_event()
{
  if(reader->NextEvent() == kTRUE) {
    ++event;
  } else {
    printf("Reached end of stream, rewinding to first event.\n");
    event = 0;
    reader->RewindEvents();
    reader->NextEvent();
  }

  printf("Now loading event %d\n", event);
  reader->Reset();
  AliTPCRawStream input(reader);
  input.SetOldRCUFormat(kTRUE);
  reader->Select("TPC"); // ("TPC", firstRCU, lastRCU);

  x->DropAllSectors();
  x->LoadRaw(input, kTRUE, kTRUE);

  printf("Updating scene\n");
  s->IncRTS();
  t->IncRTS();
  gReve->Redraw3D();
}

void tpc_raw_pad_dump(Int_t s, Int_t r, Int_t p)
{
  if(r >= TPCSectorData::GetInnSeg().GetNRows()) {
    r -=  TPCSectorData::GetInnSeg().GetNRows();
    s += 36;
  }

  reader->Reset();
  AliTPCRawStream input(reader);
  input.SetOldRCUFormat(kTRUE);
  // reader->Select(0, firstRCU, lastRCU);

  Int_t sector = input.GetSector();
  Int_t row    = input.GetRow();

  while (input.Next()) {
    if (input.IsNewRow()) {
      sector = input.GetSector();
      row    = input.GetRow();
    }
    if(sector != s || row != r) continue;

    Int_t signal = input.GetSignal();
    Int_t pad    = input.GetPad();
    Int_t time   = input.GetTime();

    if(pad == p)
      printf("%d %d\n", time, signal);
  }
}
