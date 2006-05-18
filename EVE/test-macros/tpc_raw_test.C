class AliRawReaderFile;

namespace Alieve {
class TPCData;
}

Alieve::TPCData*  x = 0;
AliRawReaderFile* reader = 0;

void tpc_raw_test()
{
  gROOT->Macro("alieve_loadlibs.C");
  gSystem->Load("libAlieve");

  reader = new AliRawReaderFile("raw0");
  reader->Reset();
  reader->NextEvent();
  AliTPCRawStream input(reader);

  x = new Alieve::TPCData;
  //x->SetSectorBlockSize(8192);
  //x->SetLoadThreshold(5);
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
