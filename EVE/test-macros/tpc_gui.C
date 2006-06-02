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
#include <Reve/RGTopFrame.h>
#include <Alieve/TPCData.h>
#include <Alieve/TPCSector2D.h>
#include <Alieve/TPCSector3D.h>

#include <RAW/AliRawReaderRoot.h>
#include <TPC/AliTPCRawStream.h>

#include <TSystem.h>
#include <TStyle.h>

#endif


using namespace Alieve;

void tpc_gui(const char *file=0, Int_t ievent=0)
{
  gStyle->SetPalette(1, 0);

  TPCLoader* l = new TPCLoader;
  TPCData*   d = new TPCData;
  // d->SetLoadPedestal(5);
  d->SetLoadThreshold(5);
  d->SetAutoPedestal(kTRUE);
  l->SetData(d);

  TGListTreeItem* loader_item = gReve->AddRenderElement(l);

  if(file != 0) {
    l->SetFile(file);
    l->OpenFile();
    l->GotoEvent(ievent);
  }
}
