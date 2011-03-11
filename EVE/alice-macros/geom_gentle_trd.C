// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

TEveGeoShape* geom_gentle_trd()
{
  TFile f("$ALICE_ROOT/EVE/alice-data/gentle_geo_trd.root");
  TEveGeoShapeExtract* gse = (TEveGeoShapeExtract*) f.Get("Gentle TRD");
  TEveGeoShape* gsre = TEveGeoShape::ImportShapeExtract(gse);
  gEve->AddGlobalElement(gsre);
  f.Close();

  const Int_t smInstalled[]={0, 1, 7, 8, 9, 10, 11, 15, 16, 17};
  const Int_t nInstalled = static_cast<Int_t>(sizeof(smInstalled)/sizeof(Int_t));
  Int_t sm = 0;
  // Fix visibility, color and transparency
  gsre->SetRnrSelf(kFALSE);
  for (TEveElement::List_i i = gsre->BeginChildren(); i != gsre->EndChildren(); ++i)
  {
    TEveGeoShape* lvl1 = (TEveGeoShape*) *i;
    lvl1->SetRnrSelf(kFALSE);
    for (TEveElement::List_i j = lvl1->BeginChildren(); j != lvl1->EndChildren(); ++j)
    {
      TEveGeoShape* lvl2 = (TEveGeoShape*) *j;
      lvl2->SetRnrSelf(kFALSE);
      for(Int_t ism(nInstalled); ism--;){
        if ( sm == smInstalled[ism] ){
          lvl2->SetRnrSelf(kTRUE);
          break; 
        }
      }
      lvl2->SetMainColor(3);
      lvl2->SetMainTransparency(80);

      ++sm;
    }
  }

  return gsre;
}
