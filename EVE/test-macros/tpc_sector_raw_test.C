class AliRawReaderRoot;

namespace Alieve {
class TPCData;
}

Alieve::TPCData*  x = 0;
AliRawReaderRoot* reader = 0;

void tpc_sector_raw_test(const char *file = "",Int_t ievent = 0)
{
  gROOT->Macro("alieve_loadlibs.C");
  gSystem->Load("libAlieve");

  reader = new AliRawReaderRoot(file);
  reader->RequireHeader(kFALSE);
  reader->Reset();
  for(Int_t i = 0; i <= ievent; i++)
    reader->NextEvent();
  AliTPCRawStreamOld input(reader);
  reader->SelectEquipment(-1);

  x = new Alieve::TPCData;
  //x->SetSectorBlockSize(8192);
  x->SetLoadThreshold(65);
  x->CreateAllSectors();
  x->LoadRaw(input, kFALSE);

  gStyle->SetPalette(1, 0);

  Alieve::TPCSector2D* s = new Alieve::TPCSector2D();
  s->SetDataSource(x);
  s->SetMainColor(36);
  gReve->AddRenderElement(s);
  gReve->DrawRenderElement(s);

}

void tpc_raw_pad_dump(Int_t s, Int_t r, Int_t p)
{
  reader->Reset();
  reader->NextEvent();

  if(r >= Alieve::TPCSectorData::GetInnSeg().fNRows) {
    r -= Alieve::TPCSectorData::GetInnSeg().fNRows;
    s += 36;
  }

  AliTPCRawStream input(reader);
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
