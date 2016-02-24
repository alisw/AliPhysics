// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TTree.h>
#include <TString.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TGeoMatrix.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEveTrans.h>
#include <TEveGeoShape.h>

#include <AliEveEventManager.h>
#include <AliRawReader.h>
#include <AliACORDERawStream.h>
#endif

TString acorde_module_path(Int_t module);

Color_t g_acorde_raw_color_on  = kRed;
Color_t g_acorde_raw_color_off = kBlue;

UChar_t g_acorde_raw_transp_on  = 30;
UChar_t g_acorde_raw_transp_off = 60;

void acorde_raw()
{
    printf("*** RAW ACORDE ***");
    
  // From Mario RC

  AliEveEventManager::Instance()->AssertGeometry();

  AliRawReader       * reader = AliEveEventManager::AssertRawReader();
  AliACORDERawStream * stream = new AliACORDERawStream(reader);

  stream->Reset();
  stream->Next();

  UInt_t dy[4];
  dy[0] = stream->GetWord(0);
  dy[1] = stream->GetWord(1);
  dy[2] = stream->GetWord(2);
  dy[3] = stream->GetWord(3);

  printf ("ACORDE event 0x%08x 0x%08x 0x%08x 0x%08x\n", dy[0], dy[1], dy[2], dy[3]);

  if (acorde_module_path(0).IsNull())
  {
    Warning("acorde_raw", "Missing / wrong ACORDE module geometry.");
    return;
  }

  TEveElementList* acorde = new TEveElementList("ACORDE_Raw");

  gEve->AddElement(acorde);

  for (Int_t module=0; module < 60; module++)
  {
    TString path = acorde_module_path(module);

    if ( ! gGeoManager->cd(path))
    {
      Warning("acorde_raw", "Module id=%d, path='%s' not found.", module, path.Data());
      continue;
    }

    // From Matevz:
    // Here check state and assign color, I do it partially for now.
    Int_t  word_idx = module / 30;
    Int_t  bit_idx  = module % 30;
    Bool_t val      = (dy[word_idx] & (1 << bit_idx)) != 0; 

    TEveGeoShape* eg_shape = new TEveGeoShape(TString::Format("Module %d", module),
                                              TString::Format("Module %d, %s", module, val ? "Fired" : "Not fired"));
    eg_shape->SetMainColor       (val ? g_acorde_raw_color_on  : g_acorde_raw_color_off);
    eg_shape->SetMainTransparency(val ? g_acorde_raw_transp_on : g_acorde_raw_transp_off);
    eg_shape->SetPickable(kTRUE);
    eg_shape->RefMainTrans().SetFrom(*gGeoManager->GetCurrentMatrix());
    eg_shape->SetShape((TGeoShape*) gGeoManager->GetCurrentVolume()->GetShape()->Clone());

    acorde->AddElement(eg_shape);
  }

  delete stream;
  gEve->Redraw3D();
}

//==============================================================================
//==============================================================================

TString acorde_module_path(Int_t module)
{
  if (module < 0 || module > 59)
  {
    Error("acorde_module_path", "module %d out of range.", module);
    return "";
  }

  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(Form("ACORDE/Array%d", module));
  if (!pne) return "";

  return Form("%s/ACORDESCINTILLATORMODULE_6", pne->GetTitle());

}
