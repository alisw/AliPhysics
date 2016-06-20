
/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TString.h>
#include <TGeoManager.h>
#include <TGeoPhysicalNode.h>
#include <TEveTrans.h>
#include <TEveManager.h>
#include <TEveElement.h>
#include <TEveGeoShape.h>

#include <AliEveEventManager.h>
#include <AliESDEvent.h>
#include <AliESDACORDE.h>
#endif

TString acorde_module_path(Int_t module);

Color_t g_acorde_esd_color_on  = kRed;
Color_t g_acorde_esd_color_off = kBlue;

UChar_t g_acorde_esd_transp_on  = 30;
UChar_t g_acorde_esd_transp_off = 60;

void acorde_esd()
{

  // From Mario RC:

  AliEveEventManager::Instance()->AssertGeometry();
  //AliESDEvent* fESD = new AliESDEvent();
  AliESDACORDE *acordeESD = AliEveEventManager::Instance()->AssertESD()->GetACORDEData();


  if (acorde_module_path(0).IsNull())
  {
    Warning("acorde_esd", "Missing / wrong ACORDE module geometry.");
    return;
  }

  TEveElementList* acorde = new TEveElementList("ACORDE_Esd");

  gEve->AddElement(acorde);

  for (Int_t module=0; module < 60; module++)
  {
    TString path = acorde_module_path(module);

    if ( ! gGeoManager->cd(path))
    {
      Warning("acorde_esd", "Module id=%d, path='%s' not found.", module, path.Data());
      continue;
    }

    Bool_t val = kFALSE;
    if (acordeESD->GetHitChannel(module)) val=kTRUE;  
    // From Matevz:
    TEveGeoShape* eg_shape = new TEveGeoShape(TString::Format("Module %d", module),
                                              TString::Format("Module %d, %s", module, val ? "Fired" : "Not fired"));
    eg_shape->SetMainColor       (val ? g_acorde_esd_color_on  : g_acorde_esd_color_off);
    eg_shape->SetMainTransparency(val ? g_acorde_esd_transp_on : g_acorde_esd_transp_off);
    eg_shape->SetPickable(kTRUE);
    eg_shape->RefMainTrans().SetFrom(*gGeoManager->GetCurrentMatrix());
    eg_shape->SetShape((TGeoShape*) gGeoManager->GetCurrentVolume()->GetShape()->Clone());

    acorde->AddElement(eg_shape);
  }

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
