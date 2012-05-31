// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// Function to spawn a gui for reading rootified raw-data from TPC sector test.
//
// To use:
// a) select AliEveTPCLoader entry in the list-tree view;
//    you'll get a dialog to steer the data-loading process in an adjacent window
// b) to select a ROOT file containing the raw-data double-click on 'File:'
//    text entry to spawn a file-dialog or type in the name
// c) click open to actually open the file and load an event

#ifdef __CINT__

class AliEveTPCData;
class AliEveTPCSector2D;
class AliEveTPCSector3D;

#else

#include <TEveManager.h>
#include <EveDet/AliEveTPCData.h>
#include <EveDet/AliEveTPCSector2D.h>
#include <EveDet/AliEveTPCSector3D.h>
#include <EveDet/AliEveTPCLoader.h>

#include <TSystem.h>
#include <TStyle.h>

#endif


AliEveTPCLoader* tpc_loader = 0;

void tpc_gui(const char *file=0, Int_t ievent=0)
{
  gStyle->SetPalette(1, 0);

  AliEveTPCLoader* l = tpc_loader = new AliEveTPCLoader;
  AliEveTPCData*   d = new AliEveTPCData;
  // d->SetLoadPedestal(5);
  d->SetLoadThreshold(5);
  d->SetAutoPedestal(kTRUE);
  l->SetData(d);
  // l->SetDoubleSR(kTRUE);
  l->SetInitParams(40, 980, 10); // min-time, max-time, threshold
  // l->SetTPCEquipementMap("EquipmentIdMap.data");

  gEve->AddElement(l);

  if(file != 0) {
    l->SetFile(file);
    l->OpenFile();
    l->GotoEvent(ievent);
  }
}
