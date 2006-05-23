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

AliRawReaderRoot* reader = 0;
Int_t             event  = -1;

void tpc_sector_raw_test(const char *file = "", Int_t ievent = 0)
{
  // gROOT->Macro("alieve_loadlibs.C");
  // Only sub-set of ALICE libraries is needed:
  gSystem->Load("libESD");
  gSystem->Load("libSTEER");
  gSystem->Load("libRAWData");
  gSystem->Load("libTPCbase");
  gSystem->Load("libTPCrec");
  gSystem->Load("libTPCsim");
  // ALICE visualization
  gSystem->Load("libAlieve");

  gStyle->SetPalette(1, 0);

  reader = new AliRawReaderRoot(file);
  reader->RequireHeader(kFALSE);
  reader->Reset();
  for(Int_t i=0; i<ievent; ++i, ++event)
    reader->NextEvennt();

  x = new TPCData;
  //x->SetLoadPedestal(5);
  x->SetLoadThreshold(5);
  x->SetAutoPedestal(kTRUE);

  s = new TPCSector2D();
  // s->SetSectorID(0);  // 0 is default
  // s->SetTrans(kTRUE); // place on proper 3D coordinates
  s->SetDataSource(x);
  s->SetFrameColor(36);
  gReve->AddRenderElement(s);
  gReve->DrawRenderElement(s);

  t = new TPCSector3D();
  // t->SetSectorID(0);
  // t->SetTrans(kTRUE);
  t->SetDataSource(x);
  t->SetMaxTime(1023);
  t->SetDriftVel(2.273);
  gReve->AddRenderElement(t);
  gReve->DrawRenderElement(t);

  next_event();
}

void next_event()
{
  reader->NextEvent();
  ++event;

  printf("Now loading event %d\n", event);
  AliTPCRawStreamOld input(reader);
  reader->SelectEquipment(-1);
  x->LoadRaw(input, kTRUE, kTRUE);

  printf("Updating scene\n");
  s->IncRTS();
  t->IncRTS();
  gReve->Redraw3D();
}

void tpc_raw_pad_dump(Int_t s, Int_t r, Int_t p)
{
  reader->Reset();
  reader->NextEvent();

  if(r >= TPCSectorData::GetInnSeg().GetNRows()) {
    r -=  TPCSectorData::GetInnSeg().GetNRows();
    s += 36;
  }

  // AliTPCRawStream input(reader);
  AliTPCRawStreamOld input(reader);
  reader->SelectEquipment(-1);

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
