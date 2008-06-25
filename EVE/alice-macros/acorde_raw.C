// $Id: geom_acorde.C 23412 2008-01-18 21:04:54Z mtadel $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

TString acorde_module_path(Int_t module);


void acorde_raw()
{
  AliEveEventManager::AssertGeometry();

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

  TEveElementList* acorde = new TEveElementList("ACORDE Raw");

  gEve->AddElement(acorde);

  Int_t shape_offset = TEveGeoShape::Class()->GetDataMemberOffset("fShape");

  for (Int_t module=0; module < 60; ++module)
  {
    TString path = acorde_module_path(module);
    // printf("%2d - %s\n", i, path.Data());

    if ( ! gGeoManager->cd(path))
    {
      Warning("acorde_raw", "Module id=%d, path='%s' not found.\n", module, path.Data());
      continue;
    }

    TEveGeoShape* eg_shape = new TEveGeoShape(Form("Module %d", module));

    eg_shape->RefMainTrans().SetFrom(* gGeoManager->GetCurrentMatrix());

    // @@NEWROOT@@ Temporary hack.
    // Hack to set shape pointer, no interface available in TEveGeoShape.
    * (TGeoShape**) (((char*)eg_shape) + shape_offset) = gGeoManager->GetCurrentVolume()->GetShape();

    // From Matevz:
    // Here check ctate and assign color, I do it partially for now.
    Int_t  word_idx = module / 30;
    Int_t  bit_idx  = module % 30;
    Bool_t val      = (dy[word_idx] & (1 << bit_idx)) != 0;
    //printf("Module %2d: word_idx = %d, bit_idx = %2d => val = %d\n",
    //       module, word_idx, bit_idx, val);
    if (val)
      eg_shape->SetMainColor(2);
    else
      eg_shape->SetMainColor(4);
    eg_shape->StampColorSelection();

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

  TGeoPNEntry* pne = gGeoManager->GetAlignableEntry(Form("ACORDE/Array%d", module + 1));
  if(!pne) return "missing_pne";

  return Form("%s/ACORDE2_5", pne->GetTitle());
}
