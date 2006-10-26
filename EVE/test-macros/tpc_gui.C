// $Header$

// Function to spawn a gui for reading rootified raw-data from TPC sector test.
//
// To use:
// a) select TPCLoader entry in the list-tree view;
//    you'll get a dialog to steer the data-loading process in an adjacent window
// b) to select a ROOT file containing the raw-data double-click on 'File:'
//    text entry to spawn a file-dialog or type in the name
// c) click open to actually open the file and load an event

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

TPCLoader* tpc_loader = 0;

void tpc_gui(const char *file=0, Int_t ievent=0)
{
  gStyle->SetPalette(1, 0);

  TPCLoader* l = tpc_loader = new TPCLoader;
  TPCData*   d = new TPCData;
  // d->SetLoadPedestal(5);
  d->SetLoadThreshold(5);
  d->SetAutoPedestal(kTRUE);
  l->SetData(d);
  l->SetDoubleSR(kTRUE);
  // l->SetInitParams(40, 980, 10); // min-time, max-time, threshold
  // l->SetTPCEquipementMap("EquipmentIdMap.data");

  gReve->AddRenderElement(l);
  gReve->NotifyBrowser(l);
  gReve->DrawRenderElement(l);

  if(file != 0) {
    l->SetFile(file);
    l->OpenFile();
    l->GotoEvent(ievent);
  }
}
