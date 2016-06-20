// $Id: geom_trd_tof.C 23442 2008-01-21 16:02:24Z mtadel $
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

void geom_extr_trd()
{
  // Extract reasonably top-level mother-volumes of TRD and
  // store them as geo-shape-extract.

  gGeoManager = gEve->GetGeometry("$ALICE_ROOT/EVE/resources/geometry/alice_fullgeo.root");

  TGeoNode* tnode = gGeoManager->GetTopVolume()->FindNode("B077_1");

  TEveGeoTopNode* eve_tnode = new TEveGeoTopNode(gGeoManager, tnode);


  for (Int_t i = 0; i < 18; ++i)
  {
    TGeoNode    * node1 = tnode->GetVolume()->FindNode(Form("BSEGMO%d_1", i));
    printf("%2d - node1 = %p\n", i, node1);
    TEveGeoNode * eve_node1 = new TEveGeoNode(node1);
    eve_tnode->AddElement(eve_node1);

    TGeoNode    * node2 = node1->GetVolume()->FindNode(Form("BTRD%d_1", i));
    printf("%2d - node2 = %p\n", i, node2);
    TEveGeoNode * eve_node2 = new TEveGeoNode(node2);
    eve_node1->AddElement(eve_node2);
  }

  eve_tnode->Save("gentle_geo_trd.root", "Gentle TRD");

}
