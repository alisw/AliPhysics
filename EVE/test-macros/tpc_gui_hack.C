// $Header$

// Function to spawn a gui for reading rootified raw-data from TPC sector test.

#ifdef __CINT__

class AliRawReaderRoot;

namespace Alieve {
class TPCData;
class TPCSector2D;
class TPCSector3D;
}

#else

#include <Reve/Reve.h>
#include <Reve/ReveManager.h>
#include <Alieve/TPCData.h>
#include <Alieve/TPCSector2D.h>
#include <Alieve/TPCSector3D.h>

#include <RAW/AliRawReaderRoot.h>
#include <TPC/AliTPCRawStream.h>

#include <TSystem.h>
#include <TStyle.h>

#endif


using namespace Alieve;

TPCSectorData* su = 0;
TPCSectorData* sl = 0;

TPCLoader* loader = 0;

void tpc_gui_hack(const char *file=0, Int_t ievent=0)
{
  gStyle->SetPalette(1, 0);

  TPCLoader* l = new TPCLoader; loader = l;
  TPCData*   d = new TPCData;
  // d->SetLoadPedestal(5);
  d->SetLoadThreshold(5);
  d->SetAutoPedestal(kTRUE);
  l->SetData(d);
  l->SetDoubleSR(kTRUE);
  l->SetInitParams(50, 980, 15);
  // l->SetTPCEquipementMap("EquipmentIdMap.data");

  su = d->GetSectorData( 4, kTRUE);
  sl = d->GetSectorData(13, kTRUE);

  gReve->AddRenderElement(l);
  gReve->NotifyBrowser(l);
  gReve->DrawRenderElement(l);

  if(file != 0) {
    l->SetFile(file);
    l->OpenFile();
    l->GotoEvent(ievent);
  }
}

void disable_pad(Int_t row, Int_t pad, TPCSectorData* sd,
                 Int_t thrExt=10, Float_t thrFac=2)
{
  if(row < 0 || row >= TPCSectorData::GetNAllRows())
    { printf("row off, %d\n", row); return; }

  Int_t npads = TPCSectorData::GetNPadsInRow(row);
  if(pad < 0) pad = npads + pad;
  if(pad < 0 || pad >= npads) { printf("pad off\n"); return; }

  sd->AddPadRowHack(row, pad, thrExt, thrFac);
}

void disable_std()
{
  // The noisy pads in lower inner seg, left edge
  disable_pad(30, 1, sl, 200, 5);
  disable_pad(30, 2, sl, 200, 5);
  disable_pad(31, 1, sl, 200, 5);

  // The noisy pads in lower inner seg, middle
  disable_pad(31, 44, sl);
  disable_pad(31, 43, sl);
  disable_pad(30, 43, sl);

  for(Int_t r=16; r<32; ++r) {
    disable_pad(r, 0, su, 20, 3);
    disable_pad(r, 1, su, 20, 3);
    disable_pad(r, 2, su, 20, 3);
    disable_pad(r, 3, su, 20, 3);
    disable_pad(r, 4, su, 20, 3);
  }

  { // Top 12, 4 pads on the negative side.
    TPCSectorData::SegmentInfo& o1si = TPCSectorData::GetOut1Seg();
    Int_t last = o1si.GetLastRow();
    for(Int_t r=last - 11; r<=last; ++r) {
      disable_pad(r, -1, su, 20, 4);
      disable_pad(r, -2, su, 10, 3);
      disable_pad(r, -3, su, 10, 3);
      disable_pad(r, -4, su,  5, 2.5);
      disable_pad(r, -5, su,  5, 2.5);
      disable_pad(r,  0, su,  5, 2.5);
      disable_pad(r,  1, su, 20, 4);
      disable_pad(r,  2, su, 10, 3);
      disable_pad(r,  3, su, 10, 3);
      disable_pad(r,  4, su,  5, 2.5);

      disable_pad(r, -1, sl, 20, 4);
      disable_pad(r, -2, sl, 10, 3);
      disable_pad(r, -3, sl, 10, 3);
      disable_pad(r, -4, sl,  5, 2.5);
      disable_pad(r, -5, sl,  5, 2.5);
    }

    disable_pad(last-8, -4, sl, 5, 2.5);

    disable_pad(last-9, -2, sl, 5, 2);
    disable_pad(last-9, -3, sl, 5, 2.5);
    disable_pad(last-9, -4, sl, 5, 2.5);
    disable_pad(last-9, -5, sl, 5, 2.5);
    disable_pad(last-8, -5, sl, 5, 2.5);
  }

  {
    TPCSectorData::SegmentInfo& o2si = TPCSectorData::GetOut2Seg();
    Int_t first = o2si.GetFirstRow();
    Int_t last  = o2si.GetLastRow();

    for(Int_t r=first; r<=last; ++r) {
      disable_pad(r, -1, su, 30, 3.5);
      disable_pad(r, -2, su, 25, 3);
      disable_pad(r, -3, su, 20, 3);
      disable_pad(r,  0, su, 30, 3.5);
      disable_pad(r,  1, su, 25, 3);

      disable_pad(r,  0, sl, 30, 3.5);
      disable_pad(r,  1, sl, 25, 3);
      disable_pad(r,  2, sl, 20, 3);
      disable_pad(r,  3, sl, 20, 3);
      disable_pad(r,  4, sl, 20, 3);
      disable_pad(r,  5, sl, 20, 3);
      disable_pad(r, -1, sl, 30, 3.5);
      disable_pad(r, -2, sl, 25, 3);
      disable_pad(r, -3, sl, 25, 3);
      disable_pad(r, -4, sl, 25, 3);
      disable_pad(r, -5, sl, 25, 3);
    }

    for(Int_t pad=3; pad<30; ++pad) {
      disable_pad(last, -pad-1, sl, 20, 3);
      disable_pad(last,    pad, su, 20, 3);
    }
    for(Int_t pad=3; pad<15; ++pad) {
      disable_pad(last-1, -pad-1, sl, 20, 3);
      disable_pad(last-2, -pad-1, sl, 20, 3);
    }
    disable_pad(last, 2, su, 20, 3);
  }

  loader->LoadEvent(); loader->UpdateSectors();
}
