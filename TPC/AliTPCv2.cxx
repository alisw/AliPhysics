/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber version 2 -- detailed TPC and slow simulation    //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTPCv2Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//#include <stdlib.h>

#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TString.h>
#include "AliLog.h"
#include "AliMathBase.h"
#include "AliTrackReference.h"
#include "AliTPCParam.h"
#include "AliTPCTrackHitsV2.h"
#include "AliTPCv2.h"
#include "AliGeomManager.h"
#include "TGeoVolume.h"
#include "TGeoPcon.h"
#include "TGeoTube.h"
#include "TGeoPgon.h"
#include "TGeoTrd1.h"
#include "TGeoCompositeShape.h"
#include "TGeoPara.h"
#include "TGeoPhysicalNode.h"
#include "TGeoHalfSpace.h"
#include "TTreeStream.h"

ClassImp(AliTPCv2)
 
//_____________________________________________________________________________
AliTPCv2::AliTPCv2(const char *name, const char *title) :
  AliTPC(name, title),
  fIdSens(0),
  fIDrift(0),
  fSecOld(0) 
{
  //
  // Standard constructor for Time Projection Chamber version 2
  //


  SetBufferSize(128000);


//   if (fTPCParam)
//      fTPCParam->Write(fTPCParam->GetTitle());
}
 
//_____________________________________________________________________________
void AliTPCv2::CreateGeometry()
{
  //
  // Create the geometry of Time Projection Chamber version 2
  //
  //Begin_Html
  /*
    <img src="picts/AliTPC.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliTPCv2Tree.gif">
  */
  //End_Html

  //----------------------------------------------------------
  // This geometry is written using TGeo class
  // Firstly the shapes are defined, and only then the volumes
  // What is recognized by the MC are volumes
  //----------------------------------------------------------
  //
  //  tpc - this will be the mother volume
  //

  //
  // here I define a volume TPC
  // retrive the medium name with "TPC_" as a leading string
  //
  TGeoPcon *tpc = new TGeoPcon(0.,360.,30); //30 sections
  //
  tpc->DefineSection(0,-289.6,77.,278.);
  tpc->DefineSection(1,-262.1,77.,278.);
  //
  tpc->DefineSection(2,-262.1,83.1,278.);
  tpc->DefineSection(3,-260.,83.1,278.);
  //
  tpc->DefineSection(4,-260.,70.,278.);    
  tpc->DefineSection(5,-259.6,70.,278.);
  //
  tpc->DefineSection(6,-259.6,68.1,278.);
  tpc->DefineSection(7,-253.6,68.1,278.);
  //
  tpc->DefineSection(8,-253.6,68.,278.);
  tpc->DefineSection(9,-74.0,60.8,278.);
  //
  tpc->DefineSection(10,-74.0,60.1,278.);
  tpc->DefineSection(11,-73.3,60.1,278.);
  //
  tpc->DefineSection(12,-73.3,56.9,278.); 
  tpc->DefineSection(13,-68.5,56.9,278.);
  //
  tpc->DefineSection(14,-68.5,60.,278.);
  tpc->DefineSection(15,-64.7,60.,278.);
  //
  tpc->DefineSection(16,-64.7,56.9,278.); 
  tpc->DefineSection(17,73.3,56.9,278.);
  //
  tpc->DefineSection(18,73.3,60.1,278.);
  tpc->DefineSection(19,74.0,60.1,278.);
  //
  tpc->DefineSection(20,74.0,60.8,278.);
  tpc->DefineSection(21,253.6,65.5,278.);
  //
  tpc->DefineSection(22,253.6,65.6,278.);
  tpc->DefineSection(23,259.6,65.6,278.);
  //
  tpc->DefineSection(24,259.6,70.0,278.);
  tpc->DefineSection(25,260.,70.0,278.);
  //
  tpc->DefineSection(26,260.,83.1,278.);
  tpc->DefineSection(27,262.1,83.1,278.);
  //
  tpc->DefineSection(28,262.1,77.,278);
  tpc->DefineSection(29,289.6,77.,278.);
  //
  TGeoMedium *m1 = gGeoManager->GetMedium("TPC_Air");
  TGeoVolume *v1 = new TGeoVolume("TPC_M",tpc,m1);
  //
  // drift volume - sensitive volume, extended beyond the
  // endcaps, because of the alignment
  //
  TGeoPcon *dvol = new TGeoPcon(0.,360.,6);
  dvol->DefineSection(0,-260.,74.5,264.4);
  dvol->DefineSection(1,-253.6,74.5,264.4);
  //
  dvol->DefineSection(2,-253.6,76.6774,258.);
  dvol->DefineSection(3,253.6,76.6774,258.); 
  //
  dvol->DefineSection(4,253.6,74.5,264.4);
  dvol->DefineSection(5,260.,74.5,264.4);
  //
  TGeoMedium *m5 = gGeoManager->GetMedium("TPC_Ne-CO2-2");
  TGeoVolume *v9 = new TGeoVolume("TPC_Drift",dvol,m5);
  //
  v1->AddNode(v9,1);
  //
  // outer insulator
  //
  TGeoPcon *tpco = new TGeoPcon(0.,360.,6); //insulator
  //
  tpco->DefineSection(0,-256.6,264.8,278.);
  tpco->DefineSection(1,-253.6,264.8,278.);
  //
  tpco->DefineSection(2,-253.6,258.,278.);
  tpco->DefineSection(3,250.6,258.,278.);
  //
  tpco->DefineSection(4,250.6,258.,275.5);
  tpco->DefineSection(5,253.6,258.,275.5);
  //
  TGeoMedium *m2 = gGeoManager->GetMedium("TPC_CO2");
  TGeoVolume *v2 = new TGeoVolume("TPC_OI",tpco,m2);
  //
  TGeoRotation *segrot;//segment rotations
  //
  // outer containment vessel
  //
  TGeoPcon *tocv = new TGeoPcon(0.,360.,6);  // containment vessel
  //
  tocv->DefineSection(0,-256.6,264.8,278.);
  tocv->DefineSection(1,-253.6,264.8,278.);
  //
  tocv->DefineSection(2,-253.6,274.8124,278.);
  tocv->DefineSection(3,247.6,274.8124,278.);  
  //
  tocv->DefineSection(4,247.6,270.4,278.);
  tocv->DefineSection(5,250.6,270.4,278.);
  //
  TGeoMedium *m3 = gGeoManager->GetMedium("TPC_Al");
  TGeoVolume *v3 = new TGeoVolume("TPC_OCV",tocv,m3); 
  //
  TGeoTubeSeg *to1 = new TGeoTubeSeg(274.8174,277.995,252.1,0.,59.9); //epoxy
  TGeoTubeSeg *to2 = new TGeoTubeSeg(274.8274,277.985,252.1,0.,59.9); //tedlar
  TGeoTubeSeg *to3 = new TGeoTubeSeg(274.8312,277.9812,252.1,0.,59.9);//prepreg2
  TGeoTubeSeg *to4 = new TGeoTubeSeg(274.9062,277.9062,252.1,0.,59.9);//nomex
  TGeoTubeSeg *tog5 = new TGeoTubeSeg(274.8174,277.995,252.1,59.9,60.);//epoxy
  //
  TGeoMedium *sm1 = gGeoManager->GetMedium("TPC_Epoxy");
  TGeoMedium *sm2 = gGeoManager->GetMedium("TPC_Tedlar");
  TGeoMedium *sm3 = gGeoManager->GetMedium("TPC_Prepreg2");
  TGeoMedium *sm4 = gGeoManager->GetMedium("TPC_Nomex");
  //
  TGeoVolume *tov1 = new TGeoVolume("TPC_OCV1",to1,sm1);
  TGeoVolume *tov2 = new TGeoVolume("TPC_OCV2",to2,sm2);
  TGeoVolume *tov3 = new TGeoVolume("TPC_OCV3",to3,sm3);
  TGeoVolume *tov4 = new TGeoVolume("TPC_OCV4",to4,sm4);
  TGeoVolume *togv5 = new TGeoVolume("TPC_OCVG5",tog5,sm1);
  //-------------------------------------------------------
  //  Tpc Outer Field Cage
  //  daughters - composite (sandwich)
  //-------------------------------------------------------

  TGeoPcon *tofc = new TGeoPcon(0.,360.,6);
  //
  tofc->DefineSection(0,-253.6,258.,269.6);
  tofc->DefineSection(1,-250.6,258.,269.6);
  //
  tofc->DefineSection(2,-250.6,258.,260.0676); 
  tofc->DefineSection(3,250.6,258.,260.0676);
  //
  tofc->DefineSection(4,250.6,258.,275.5);
  tofc->DefineSection(5,253.6,258.,275.5);
  //
  TGeoVolume *v4 = new TGeoVolume("TPC_TOFC",tofc,m3); 
  //sandwich
  TGeoTubeSeg *tf1 = new TGeoTubeSeg(258.0,260.0676,252.1,0.,59.9); //tedlar
  TGeoTubeSeg *tf2 = new TGeoTubeSeg(258.0038,260.0638,252.1,0.,59.9); //prepreg3
  TGeoTubeSeg *tf3 = new TGeoTubeSeg(258.0338,260.0338,252.1,0.,59.9);//nomex
  TGeoTubeSeg *tfg4 = new TGeoTubeSeg(258.0,260.0676,252.1,59.9,60.); //epoxy glue
  //
  TGeoMedium *sm5 = gGeoManager->GetMedium("TPC_Prepreg3");
  //
  TGeoVolume *tf1v = new TGeoVolume("TPC_OFC1",tf1,sm2);
  TGeoVolume *tf2v = new TGeoVolume("TPC_OFC2",tf2,sm5);
  TGeoVolume *tf3v = new TGeoVolume("TPC_OFC3",tf3,sm4);
  TGeoVolume *tfg4v = new TGeoVolume("TPC_OFCG4",tfg4,sm1);
  //
  // outer part - positioning
  //
  tov1->AddNode(tov2,1); tov2->AddNode(tov3,1); tov3->AddNode(tov4,1);//ocv
  //
  tf1v->AddNode(tf2v,1); tf2v->AddNode(tf3v,1);//ofc
  //
  TGeoVolumeAssembly *t200 = new TGeoVolumeAssembly("TPC_OCVSEG");
  TGeoVolumeAssembly *t300 = new TGeoVolumeAssembly("TPC_OFCSEG");
  //
  // assembly OCV and OFC
  //
  // 1st - no rotation
  t200->AddNode(tov1,1); t200->AddNode(togv5,1);
  t300->AddNode(tf1v,1); t300->AddNode(tfg4v,1);
  // 2nd - rotation 60 deg
  segrot = new TGeoRotation();
  segrot->RotateZ(60.);
  t200->AddNode(tov1,2,segrot); t200->AddNode(togv5,2,segrot);
  t300->AddNode(tf1v,2,segrot); t300->AddNode(tfg4v,2,segrot);
  // 3rd rotation 120 deg
  segrot = new TGeoRotation();
  segrot->RotateZ(120.);
  t200->AddNode(tov1,3,segrot); t200->AddNode(togv5,3,segrot);
  t300->AddNode(tf1v,3,segrot); t300->AddNode(tfg4v,3,segrot);
  //4th rotation 180 deg
  segrot = new TGeoRotation();
  segrot->RotateZ(180.);
  t200->AddNode(tov1,4,segrot); t200->AddNode(togv5,4,segrot);
  t300->AddNode(tf1v,4,segrot); t300->AddNode(tfg4v,4,segrot);
  //5th rotation 240 deg
  segrot = new TGeoRotation();
  segrot->RotateZ(240.);
  t200->AddNode(tov1,5,segrot); t200->AddNode(togv5,5,segrot);
  t300->AddNode(tf1v,5,segrot); t300->AddNode(tfg4v,5,segrot);
  //6th rotation 300 deg
  segrot = new TGeoRotation();
  segrot->RotateZ(300.);
  t200->AddNode(tov1,6,segrot); t200->AddNode(togv5,6,segrot);
  t300->AddNode(tf1v,6,segrot); t300->AddNode(tfg4v,6,segrot);
  //
  v3->AddNode(t200,1,new TGeoTranslation(0.,0.,-1.5)); v4->AddNode(t300,1);
  //
  v2->AddNode(v3,1); v2->AddNode(v4,1); 
  //
  v1->AddNode(v2,1);
  //--------------------------------------------------------------------
  // Tpc Inner INsulator (CO2) 
  // the cones, the central drum and the inner f.c. sandwich with a piece
  // of the flane will be placed in the TPC
  //--------------------------------------------------------------------
  TGeoPcon *tpci = new TGeoPcon(0.,360.,4);
  //
  tpci->DefineSection(0,-253.6,68.4,76.6774);
  tpci->DefineSection(1,-74.0,61.2,76.6774);
  //
  tpci->DefineSection(2,74.0,61.2,76.6774);  
  //
  tpci->DefineSection(3,253.6,65.9,76.6774);
  //
  TGeoVolume *v5 = new TGeoVolume("TPC_INI",tpci,m2);
  //
  // now the inner field cage - only part of flanges (2 copies)
  //
  TGeoTube *tif1 = new TGeoTube(69.9,76.6774,1.5); 
  TGeoVolume *v6 = new TGeoVolume("TPC_IFC1",tif1,m3);
  //
 //---------------------------------------------------------
  // Tpc Inner Containment vessel - Muon side
  //---------------------------------------------------------
  TGeoPcon *tcms = new TGeoPcon(0.,360.,10);
  //
  tcms->DefineSection(0,-259.1,68.1,74.2);
  tcms->DefineSection(1,-253.6,68.1,74.2);
  //
  tcms->DefineSection(2,-253.6,68.1,68.4);
  tcms->DefineSection(3,-74.0,60.9,61.2);
  //
  tcms->DefineSection(4,-74.0,60.1,61.2);
  tcms->DefineSection(5,-73.3,60.1,61.2);
  //
  tcms->DefineSection(6,-73.3,56.9,61.2);
  tcms->DefineSection(7,-73.0,56.9,61.2);
  //
  tcms->DefineSection(8,-73.0,56.9,58.8);
  tcms->DefineSection(9,-71.3,56.9,58.8);
  //
  TGeoVolume *v7 = new TGeoVolume("TPC_ICVM",tcms,m3);
  //-----------------------------------------------
  // inner containment vessel - shaft side
  //-----------------------------------------------
  TGeoPcon *tcss = new TGeoPcon(0.,360.,10);
  //
  tcss->DefineSection(0,71.3,56.9,58.8);
  tcss->DefineSection(1,73.0,56.9,58.8);
  //
  tcss->DefineSection(2,73.0,56.9,61.2);
  tcss->DefineSection(3,73.3,56.9,61.2);
  //  
  tcss->DefineSection(4,73.3,60.1,61.2);
  tcss->DefineSection(5,74.0,60.1,61.2);
  //
  tcss->DefineSection(6,74.0,60.9,61.2);
  tcss->DefineSection(7,253.6,65.6,65.9);
  //
  tcss->DefineSection(8,253.6,65.6,74.2);
  tcss->DefineSection(9,258.1,65.6,74.2);
  //
  TGeoVolume *v8 = new TGeoVolume("TPC_ICVS",tcss,m3);
  //-----------------------------------------------
  //  Inner field cage
  //  define 4 parts and make an assembly
  //-----------------------------------------------
  // part1 - Al - 2 copies
  TGeoTube *t1 = new TGeoTube(76.6774,78.845,0.75);
  TGeoVolume *tv1 = new TGeoVolume("TPC_IFC2",t1,m3);
  // sandwich - outermost parts - 2 copies
  //
  // segment outermost
  //
  TGeoTubeSeg *t2 = new TGeoTubeSeg(76.6774,78.845,74.175,350.,109.4); // tedlar 38 microns
  TGeoTubeSeg *t3 = new TGeoTubeSeg(76.6812,78.8412,74.175,350.,109.4); // prepreg2 500 microns
  TGeoTubeSeg *t4 = new TGeoTubeSeg(76.7312,78.7912,74.175,350.,109.4); // prepreg3 300 microns
  TGeoTubeSeg *t5 = new TGeoTubeSeg(76.7612,78.7612,74.175,350.,109.4); // nomex 2 cm
  TGeoTubeSeg *tepox1 = new TGeoTubeSeg(76.6774,78.845,74.175,109.4,110.);//epoxy
  TGeoTubeSeg *tpr1 = new TGeoTubeSeg(78.845,78.885,74.175,109.,111.);
  
  // volumes for the outer part  
  TGeoVolume *tv2 = new TGeoVolume("TPC_IFC3",t2,sm2);
  TGeoVolume *tv3 = new TGeoVolume("TPC_IFC4",t3,sm3);
  TGeoVolume *tv4 = new TGeoVolume("TPC_IFC5",t4,sm5);
  TGeoVolume *tv5 = new TGeoVolume("TPC_IFC6",t5,sm4);
  TGeoVolume *tvep1 = new TGeoVolume("TPC_IFEPOX1",tepox1,sm1); 
  TGeoVolume *tvpr1 = new TGeoVolume("TPC_PRSTR1",tpr1,sm2); 
  //
  // middle parts - 2 copies
  //
  // segment middle
  //
  TGeoTubeSeg *t6 = new TGeoTubeSeg(76.6774,78.795,5.,350.,109.4); // tedlar 38 microns
  TGeoTubeSeg *t7 = new TGeoTubeSeg(76.6812,78.7912,5.,350.,109.4); // prepreg2 250 microns
  TGeoTubeSeg *t8 = new TGeoTubeSeg(76.7062,78.7662,5.,350.,109.4); // prepreg3 300 microns
  TGeoTubeSeg *t9 = new TGeoTubeSeg(76.7362,78.7362,5.,350.,109.4); // nomex 2 cm
  TGeoTubeSeg *tepox2 = new TGeoTubeSeg(76.6774,78.795,5.,109.4,110.);//epoxy
  TGeoTubeSeg *tpr2 = new TGeoTubeSeg(78.795,78.835,5.,109.,111.);
  // volumes for the middle part
  TGeoVolume *tv6 = new TGeoVolume("TPC_IFC7",t6,sm2);
  TGeoVolume *tv7 = new TGeoVolume("TPC_IFC8",t7,sm3);
  TGeoVolume *tv8 = new TGeoVolume("TPC_IFC9",t8,sm5);
  TGeoVolume *tv9 = new TGeoVolume("TPC_IFC10",t9,sm4);
  TGeoVolume *tvep2 = new TGeoVolume("TPC_IFEPOX2",tepox2,sm1);
  TGeoVolume *tvpr2 = new TGeoVolume("TPC_PRSTR2",tpr2,sm2);
  // central part - 1 copy
  // 
  // segment central part
  //
  TGeoTubeSeg *t10 = new TGeoTubeSeg(76.6774,78.785,93.75,350.,109.4); // tedlar 38 microns 
  TGeoTubeSeg *t11 = new TGeoTubeSeg(76.6812,78.7812,93.75,350.,109.4); // prepreg3 500 microns
  TGeoTubeSeg *t12 = new TGeoTubeSeg(76.7312,78.7312,93.75,350.,109.4); // nomex 2 cm 
  TGeoTubeSeg *tepox3 = new TGeoTubeSeg(76.6774,78.785,93.75,109.4,110.);//epoxy
  TGeoTubeSeg *tpr3 = new TGeoTubeSeg(78.785,78.825,93.75,109.,111.);
  // volumes for the central part
  TGeoVolume *tv10 = new TGeoVolume("TPC_IFC11",t10,sm2);
  TGeoVolume *tv11 = new TGeoVolume("TPC_IFC12",t11,sm5);
  TGeoVolume *tv12 = new TGeoVolume("TPC_IFC13",t12,sm4);
  TGeoVolume *tvep3 = new TGeoVolume("TPC_IFEPOX3",tepox3,sm1);
  TGeoVolume *tvpr3 = new TGeoVolume("TPC_PRSTR3",tpr3,sm2);  
  //
  // creating a sandwich for the outer par,t tv2 is the mother
  //
  tv2->AddNode(tv3,1); tv3->AddNode(tv4,1); tv4->AddNode(tv5,1);
  //
  // creating a sandwich for the middle part, tv6 is the mother
  //
  tv6->AddNode(tv7,1); tv7->AddNode(tv8,1); tv8->AddNode(tv9,1);
  //
  // creating a sandwich for the central part, tv10 is the mother
  //
  tv10->AddNode(tv11,1); tv11->AddNode(tv12,1);
  //
  TGeoVolumeAssembly *tv100 = new TGeoVolumeAssembly("TPC_IFC"); // ifc itself - 3 segments

  //
  // first segment - no rotation
  //
  // central
  tv100->AddNode(tv10,1); //sandwich
  tv100->AddNode(tvep3,1);//epoxy
  tv100->AddNode(tvpr3,1);//prepreg strip				     
  // middle
  tv100->AddNode(tv6,1,new TGeoTranslation(0.,0.,-98.75)); //sandwich1
  tv100->AddNode(tv6,2,new TGeoTranslation(0.,0.,98.75)); // sandwich2
  tv100->AddNode(tvep2,1,new TGeoTranslation(0.,0.,-98.75)); //epoxy
  tv100->AddNode(tvep2,2,new TGeoTranslation(0.,0.,98.75)); //epoxy
  tv100->AddNode(tvpr2,1,new TGeoTranslation(0.,0.,-98.75));//prepreg strip
  tv100->AddNode(tvpr2,2,new TGeoTranslation(0.,0.,98.75));
  // outer
  tv100->AddNode(tv2,1,new TGeoTranslation(0.,0.,-177.925)); //sandwich
  tv100->AddNode(tv2,2,new TGeoTranslation(0.,0.,177.925));
  tv100->AddNode(tvep1,1,new TGeoTranslation(0.,0.,-177.925)); //epoxy
  tv100->AddNode(tvep1,2,new TGeoTranslation(0.,0.,177.925));
  tv100->AddNode(tvpr1,1,new TGeoTranslation(0.,0.,-177.925));//prepreg strip
  tv100->AddNode(tvpr1,2,new TGeoTranslation(0.,0.,-177.925));
  //
  // second segment - rotation 120 deg.
  //
  segrot = new TGeoRotation();
  segrot->RotateZ(120.);
  //
  // central	
  tv100->AddNode(tv10,2,segrot); //sandwich
  tv100->AddNode(tvep3,2,segrot);//epoxy
  tv100->AddNode(tvpr3,2,segrot);//prepreg strip
  // middle
  tv100->AddNode(tv6,3,new TGeoCombiTrans(0.,0.,-98.75,segrot)); //sandwich1
  tv100->AddNode(tv6,4,new TGeoCombiTrans(0.,0.,98.75,segrot)); // sandwich2
  tv100->AddNode(tvep2,3,new TGeoCombiTrans(0.,0.,-98.75,segrot)); //epoxy
  tv100->AddNode(tvep2,4,new TGeoCombiTrans(0.,0.,98.75,segrot)); //epoxy
  tv100->AddNode(tvpr2,3,new TGeoCombiTrans(0.,0.,-98.75,segrot));//prepreg strip
  tv100->AddNode(tvpr2,4,new TGeoCombiTrans(0.,0.,98.75,segrot));
  //outer
  tv100->AddNode(tv2,3,new TGeoCombiTrans(0.,0.,-177.925,segrot));//sandwich
  tv100->AddNode(tv2,4,new TGeoCombiTrans(0.,0.,177.925,segrot));
  tv100->AddNode(tvep1,3,new TGeoCombiTrans(0.,0.,-177.925,segrot));//epoxy
  tv100->AddNode(tvep1,4,new TGeoCombiTrans(0.,0.,177.925,segrot));
  tv100->AddNode(tvpr1,3,new TGeoCombiTrans(0.,0.,-177.925,segrot));//prepreg strip
  tv100->AddNode(tvpr1,4,new TGeoCombiTrans(0.,0.,177.925,segrot));
  //
  //  third segment - rotation 240 deg.
  //
  segrot = new TGeoRotation();
  segrot->RotateZ(240.);
  //
  // central	
  tv100->AddNode(tv10,3,segrot); //sandwich
  tv100->AddNode(tvep3,3,segrot);//epoxy
  tv100->AddNode(tvpr3,3,segrot);//prepreg strip
  // middle
  tv100->AddNode(tv6,5,new TGeoCombiTrans(0.,0.,-98.75,segrot)); //sandwich1
  tv100->AddNode(tv6,6,new TGeoCombiTrans(0.,0.,98.75,segrot)); // sandwich2
  tv100->AddNode(tvep2,5,new TGeoCombiTrans(0.,0.,-98.75,segrot)); //epoxy
  tv100->AddNode(tvep2,6,new TGeoCombiTrans(0.,0.,98.75,segrot)); //epoxy
  tv100->AddNode(tvpr2,5,new TGeoCombiTrans(0.,0.,-98.75,segrot));//prepreg strip
  tv100->AddNode(tvpr2,6,new TGeoCombiTrans(0.,0.,98.75,segrot));
  //outer
  tv100->AddNode(tv2,5,new TGeoCombiTrans(0.,0.,-177.925,segrot));//sandwich
  tv100->AddNode(tv2,6,new TGeoCombiTrans(0.,0.,177.925,segrot));
  tv100->AddNode(tvep1,5,new TGeoCombiTrans(0.,0.,-177.925,segrot));//epoxy
  tv100->AddNode(tvep1,6,new TGeoCombiTrans(0.,0.,177.925,segrot));
  tv100->AddNode(tvpr1,5,new TGeoCombiTrans(0.,0.,-177.925,segrot));//prepreg strip
  tv100->AddNode(tvpr1,6,new TGeoCombiTrans(0.,0.,177.925,segrot));
  // Al parts - rings
  tv100->AddNode(tv1,1,new TGeoTranslation(0.,0.,-252.85));
  tv100->AddNode(tv1,2,new TGeoTranslation(0.,0.,252.85));
  //
  v5->AddNode(v6,1, new TGeoTranslation(0.,0.,-252.1));
  v5->AddNode(v6,2, new TGeoTranslation(0.,0.,252.1));
  v1->AddNode(v5,1); v1->AddNode(v7,1); v1->AddNode(v8,1); 
  v9->AddNode(tv100,1);
  //
  // central drum 
  //
  // flange + sandwich
  //
  TGeoPcon *cfl = new TGeoPcon(0.,360.,6);
  cfl->DefineSection(0,-71.1,59.7,61.2);
  cfl->DefineSection(1,-68.6,59.7,61.2);
  //
  cfl->DefineSection(2,-68.6,60.6124,61.2);
  cfl->DefineSection(3,68.6,60.6124,61.2); 
  //
  cfl->DefineSection(4,68.6,59.7,61.2);
  cfl->DefineSection(5,71.1,59.7,61.2);  
  //
  TGeoVolume *cflv = new TGeoVolume("TPC_CDR",cfl,m3);
  // sandwich
  TGeoTubeSeg *cd1 = new TGeoTubeSeg(60.6224,61.19,71.1,0.1,119.9);
  TGeoTubeSeg *cd2 = new TGeoTubeSeg(60.6262,61.1862,71.1,0.1,119.9);
  TGeoTubeSeg *cd3 = new TGeoTubeSeg(60.6462,61.1662,71.1,0.1,119.9);
  TGeoTubeSeg *cd4 = new TGeoTubeSeg(60.6562,61.1562,71.1,0.1,119.9);
  TGeoTubeSeg *tepox4 = new TGeoTubeSeg(60.6224,61.19,71.1,359.9,0.1);
  //  TGeoTube *cd1 = new TGeoTube(60.6224,61.19,71.1);
  //  TGeoTube *cd2 = new TGeoTube(60.6262,61.1862,71.1);
  //  TGeoTube *cd3 = new TGeoTube(60.6462,61.1662,71.1);  
  //  TGeoTube *cd4 = new TGeoTube(60.6562,61.1562,71.1);  
  //
  TGeoMedium *sm6 = gGeoManager->GetMedium("TPC_Prepreg1");
  TGeoMedium *sm8 = gGeoManager->GetMedium("TPC_Epoxyfm");
  TGeoVolume *cd1v = new TGeoVolume("TPC_CDR1",cd1,sm2); //tedlar
  TGeoVolume *cd2v = new TGeoVolume("TPC_CDR2",cd2,sm6);// prepreg1
  TGeoVolume *cd3v = new TGeoVolume("TPC_CDR3",cd3,sm8); //epoxy film
  TGeoVolume *cd4v = new TGeoVolume("TPC_CDR4",cd4,sm4); //nomex
  TGeoVolume *tvep4 = new TGeoVolume("TPC_IFEPOX4",tepox4,sm1);

  //
  // seals for central drum 2 copies
  //
  TGeoTube *cs = new TGeoTube(56.9,61.2,0.1);
  TGeoMedium *sm7 = gGeoManager->GetMedium("TPC_Mylar");
  TGeoVolume *csv = new TGeoVolume("TPC_CDRS",cs,sm7);
  v1->AddNode(csv,1,new TGeoTranslation(0.,0.,-71.2));
  v1->AddNode(csv,2,new TGeoTranslation(0.,0.,71.2));
  //
  // seal collars 
  TGeoPcon *se = new TGeoPcon(0.,360.,6);
  se->DefineSection(0,-72.8,59.7,61.2);
  se->DefineSection(1,-72.3,59.7,61.2);
  //
  se->DefineSection(2,-72.3,58.85,61.2);
  se->DefineSection(3,-71.6,58.85,61.2); 
  //
  se->DefineSection(4,-71.6,59.7,61.2);
  se->DefineSection(5,-71.3,59.7,61.2);  
  //
  TGeoVolume *sev = new TGeoVolume("TPC_CDCE",se,m3);
  //
  TGeoTube *si = new TGeoTube(56.9,58.8,1.); 
  TGeoVolume *siv = new TGeoVolume("TPC_CDCI",si,m3);
  //
  // define reflection matrix 
  //
  TGeoRotation *ref = new TGeoRotation("ref",90.,0.,90.,90.,180.,0.);
  //
  cd1v->AddNode(cd2v,1); cd2v->AddNode(cd3v,1); cd3v->AddNode(cd4v,1); //sandwich
  // first segment  
  cflv->AddNode(cd1v,1); cflv->AddNode(tvep4,1);
  // second segment
  segrot = new TGeoRotation();
  segrot->RotateZ(120.);
  cflv->AddNode(cd1v,2,segrot); cflv->AddNode(tvep4,2,segrot);
  // third segment
  segrot = new TGeoRotation();
  segrot->RotateZ(240.);
  cflv->AddNode(cd1v,3,segrot); cflv->AddNode(tvep4,3,segrot);
  //
  v1->AddNode(siv,1,new TGeoTranslation(0.,0.,-69.9));
  v1->AddNode(siv,2,new TGeoTranslation(0.,0.,69.9));
  v1->AddNode(sev,1); v1->AddNode(sev,2,ref); v1->AddNode(cflv,1);
  //
  // central membrane - 2 rings and a mylar membrane - assembly
  //
  TGeoTube *ih = new TGeoTube(81.05,84.05,0.3);
  TGeoTube *oh = new TGeoTube(250.,256.,0.5);
  TGeoTube *mem = new TGeoTube(84.05,250.,0.00115);

  //
  TGeoMedium *m4 = gGeoManager->GetMedium("TPC_G10");
  //
  TGeoVolume *ihv = new TGeoVolume("TPC_IHVH",ih,m3);
  TGeoVolume *ohv = new TGeoVolume("TPC_OHVH",oh,m3);
  
  TGeoVolume *memv = new TGeoVolume("TPC_HV",mem,sm7);
  //
  TGeoVolumeAssembly *cm = new TGeoVolumeAssembly("TPC_HVMEM");
  cm->AddNode(ihv,1);
  cm->AddNode(ohv,1);
  cm->AddNode(memv,1);
 
  v9->AddNode(cm,1);
  //
  // end caps - they are make as an assembly of single segments
  // containing both readout chambers
  //
  Double_t openingAngle = 10.*TMath::DegToRad();
  Double_t thick=1.5; // rib
  Double_t shift = thick/TMath::Sin(openingAngle);
  //
  Double_t lowEdge = 86.3; // hole in the wheel
  Double_t upEdge = 240.4; // hole in the wheel
  //
  new TGeoTubeSeg("sec",74.5,264.4,3.,0.,20.);
  //
  TGeoPgon *hole = new TGeoPgon("hole",0.,20.,1,4);
  //
  hole->DefineSection(0,-3.5,lowEdge-shift,upEdge-shift);
  hole->DefineSection(1,-1.5,lowEdge-shift,upEdge-shift);
  //
  hole->DefineSection(2,-1.5,lowEdge-shift,upEdge+3.-shift);
  hole->DefineSection(3,3.5,lowEdge-shift,upEdge+3.-shift);
  //
  Double_t ys = shift*TMath::Sin(openingAngle); 
  Double_t xs = shift*TMath::Cos(openingAngle);
  TGeoTranslation *tr = new TGeoTranslation("tr",xs,ys,0.);  
  tr->RegisterYourself();
  TGeoCompositeShape *chamber = new TGeoCompositeShape("sec-hole:tr");
  TGeoVolume *sv = new TGeoVolume("TPC_WSEG",chamber,m3);
  TGeoPgon *bar = new TGeoPgon("bar",0.,20.,1,2);
  bar->DefineSection(0,-3.,131.5-shift,136.5-shift);
  bar->DefineSection(1,1.5,131.5-shift,136.5-shift);
  TGeoVolume *barv = new TGeoVolume("TPC_WBAR",bar,m3);
  TGeoVolumeAssembly *ch = new TGeoVolumeAssembly("TPC_WCH");//empty segment
  //
  ch->AddNode(sv,1); ch->AddNode(barv,1,tr);
  //
  // readout chambers
  //
  // IROC first
  //
   TGeoTrd1 *ibody = new TGeoTrd1(13.8742,21.3328,4.29,21.15);
   TGeoVolume *ibdv = new TGeoVolume("TPC_IROCB",ibody,m3);
  // empty space
   TGeoTrd1 *emp = new TGeoTrd1(12.3742,19.8328,3.99,19.65);
   TGeoVolume *empv = new TGeoVolume("TPC_IROCE",emp,m1);
   ibdv->AddNode(empv,1,new TGeoTranslation(0.,-0.3,0.));
   //bars
   Double_t tga = (19.8328-12.3742)/39.3;
   Double_t xmin,xmax;
   xmin = 9.55*tga+12.3742;
   xmax = 9.95*tga+12.3742;
   TGeoTrd1 *ib1 = new TGeoTrd1(xmin,xmax,3.29,0.2);
   TGeoVolume *ib1v = new TGeoVolume("TPC_IRB1",ib1,m3);
   empv->AddNode(ib1v,1,new TGeoTranslation("tt1",0.,0.7,-9.9));
   xmin=19.4*tga+12.3742;
   xmax=19.9*tga+12.3742;
   TGeoTrd1 *ib2 = new TGeoTrd1(xmin,xmax,3.29,0.25);
   TGeoVolume *ib2v = new TGeoVolume("TPC_TRB2",ib2,m3);
   empv->AddNode(ib2v,1,new TGeoTranslation(0.,0.7,0.));
   xmin=29.35*tga+12.3742;
   xmax=29.75*tga+12.3742;
   TGeoTrd1 *ib3 = new TGeoTrd1(xmin,xmax,3.29,0.2); 
   TGeoVolume *ib3v = new TGeoVolume("TPC_IRB3",ib3,m3);    
   empv->AddNode(ib3v,1,new TGeoTranslation(0.,0.7,9.9));
   //
   // holes for connectors
   //
   TGeoBBox *conn = new TGeoBBox(0.4,0.3,4.675); // identical for iroc and oroc
   TGeoVolume *connv = new TGeoVolume("TPC_RCCON",conn,m1);
   TString fileName(gSystem->Getenv("ALICE_ROOT"));
   fileName += "/TPC/conn_iroc.dat";
   ifstream in;
   in.open(fileName.Data(), ios_base::in); // asci file
   TGeoRotation *rrr[86];
   for(Int_t i =0;i<86;i++){
      Double_t y = 3.99;
      Double_t x,z,ang;
      in>>x>>z>>ang;
      z-=26.5;
      rrr[i]= new TGeoRotation();
      rrr[i]->RotateY(ang);
      ibdv->AddNode(connv,i+1,new TGeoCombiTrans(x,y,z,rrr[i]));
   }
   in.close();
   // "cap"
   new TGeoTrd1("icap",14.5974,23.3521,1.19,24.825);
   // "hole"
   new TGeoTrd1("ihole",13.8742,21.3328,1.2,21.15);
   TGeoTranslation *tr1 = new TGeoTranslation("tr1",0.,0.,1.725);  
   tr1->RegisterYourself();
   TGeoCompositeShape *ic = new TGeoCompositeShape("icap-ihole:tr1");
   TGeoVolume *icv = new TGeoVolume("TPC_IRCAP",ic,m3);
   //
   // pad plane and wire fixations
   //
   TGeoTrd1 *pp = new TGeoTrd1(14.5974,23.3521,0.3,24.825); //pad+iso
   TGeoVolume *ppv = new TGeoVolume("TPC_IRPP",pp,m4);
   TGeoPara *f1 = new TGeoPara(.6,.5,24.825,0.,-10.,0.);
   TGeoVolume *f1v = new TGeoVolume("TPC_IRF1",f1,m4);
   TGeoPara *f2 = new TGeoPara(.6,.5,24.825,0.,10.,0.);
   TGeoVolume *f2v = new TGeoVolume("TPC_IRF2",f2,m4);
   //
   TGeoVolumeAssembly *iroc = new TGeoVolumeAssembly("TPC_IROC");
   //
   iroc->AddNode(ibdv,1);
   iroc->AddNode(icv,1,new TGeoTranslation(0.,3.1,-1.725));
   iroc->AddNode(ppv,1,new TGeoTranslation(0.,4.59,-1.725));
   tga =(23.3521-14.5974)/49.65; 
   Double_t xx = 24.825*tga+14.5974-0.6;
   iroc->AddNode(f1v,1,new TGeoTranslation(-xx,5.39,-1.725));
   iroc->AddNode(f2v,1,new TGeoTranslation(xx,5.39,-1.725));
   //
   // OROC
   //
   TGeoTrd1 *obody = new TGeoTrd1(22.2938,40.5084,4.19,51.65);
   TGeoVolume *obdv = new TGeoVolume("TPC_OROCB",obody,m3);
   TGeoTrd1 *oemp = new TGeoTrd1(20.2938,38.5084,3.89,49.65);
   TGeoVolume *oempv = new TGeoVolume("TPC_OROCE",oemp,m1);
   obdv->AddNode(oempv,1,new TGeoTranslation(0.,-0.3,0.));
   //horizontal bars
   tga=(38.5084-20.2938)/99.3;
   xmin=tga*10.2+20.2938;
   xmax=tga*10.6+20.2938;
   TGeoTrd1 *ob1 = new TGeoTrd1(xmin,xmax,2.915,0.2);
   TGeoVolume *ob1v = new TGeoVolume("TPC_ORB1",ob1,m3);
   //
   xmin=22.55*tga+20.2938;
   xmax=24.15*tga+20.2938;
   TGeoTrd1 *ob2 = new TGeoTrd1(xmin,xmax,2.915,0.8);
   TGeoVolume *ob2v = new TGeoVolume("TPC_ORB2",ob2,m3);
   //
   xmin=36.1*tga+20.2938;
   xmax=36.5*tga+20.2938;
   TGeoTrd1 *ob3 = new TGeoTrd1(xmin,xmax,2.915,0.2);
   TGeoVolume *ob3v = new TGeoVolume("TPC_ORB3",ob3,m3);
   //
   xmin=49.0*tga+20.2938;
   xmax=50.6*tga+20.2938;   
   TGeoTrd1 *ob4 = new TGeoTrd1(xmin,xmax,2.915,0.8);
   TGeoVolume *ob4v = new TGeoVolume("TPC_ORB4",ob4,m3);
   //
   xmin=63.6*tga+20.2938;
   xmax=64.0*tga+20.2938;
   TGeoTrd1 *ob5 = new TGeoTrd1(xmin,xmax,2.915,0.2);
   TGeoVolume *ob5v = new TGeoVolume("TPC_ORB5",ob5,m3);
   //
   xmin=75.5*tga+20.2938;
   xmax=77.15*tga+20.2938;
   TGeoTrd1 *ob6 = new TGeoTrd1(xmin,xmax,2.915,0.8);
   TGeoVolume *ob6v = new TGeoVolume("TPC_ORB6",ob6,m3);
   //
   xmin=88.7*tga+20.2938;
   xmax=89.1*tga+20.2938;
   TGeoTrd1 *ob7 = new TGeoTrd1(xmin,xmax,2.915,0.2);
   TGeoVolume *ob7v = new TGeoVolume("TPC_ORB7",ob7,m3);
   //
   oempv->AddNode(ob1v,1,new TGeoTranslation(0.,0.975,-39.25));
   oempv->AddNode(ob2v,1,new TGeoTranslation(0.,0.975,-26.3));
   oempv->AddNode(ob3v,1,new TGeoTranslation(0.,0.975,-13.35));
   oempv->AddNode(ob4v,1,new TGeoTranslation(0.,0.975,0.15));
   oempv->AddNode(ob5v,1,new TGeoTranslation(0.,0.975,14.15));
   oempv->AddNode(ob6v,1,new TGeoTranslation(0.,0.975,26.7));
   oempv->AddNode(ob7v,1,new TGeoTranslation(0.,0.975,39.25));
   // vertical bars
   TGeoBBox *ob8 = new TGeoBBox(0.8,2.915,5.1); 
   TGeoBBox *ob9 = new TGeoBBox(0.8,2.915,5.975);
   TGeoBBox *ob10 = new TGeoBBox(0.8,2.915,5.775);
   TGeoBBox *ob11 = new TGeoBBox(0.8,2.915,6.25);
   TGeoBBox *ob12 = new TGeoBBox(0.8,2.915,6.5);
   //
   TGeoVolume *ob8v = new TGeoVolume("TPC_ORB8",ob8,m3);
   TGeoVolume *ob9v = new TGeoVolume("TPC_ORB9",ob9,m3);
   TGeoVolume *ob10v = new TGeoVolume("TPC_ORB10",ob10,m3);
   TGeoVolume *ob11v = new TGeoVolume("TPC_ORB11",ob11,m3);
   TGeoVolume *ob12v = new TGeoVolume("TPC_ORB12",ob12,m3);
   //
   oempv->AddNode(ob8v,1,new TGeoTranslation(0.,0.975,-44.55));
   oempv->AddNode(ob8v,2,new TGeoTranslation(0.,0.975,44.55));
   oempv->AddNode(ob9v,1,new TGeoTranslation(0.,0.975,-33.075));
   oempv->AddNode(ob9v,2,new TGeoTranslation(0.,0.975,-19.525));
   oempv->AddNode(ob10v,1,new TGeoTranslation(0.,0.975,20.125));
   oempv->AddNode(ob10v,2,new TGeoTranslation(0.,0.975,33.275));
   oempv->AddNode(ob11v,1,new TGeoTranslation(0.,0.975,-6.9));
   oempv->AddNode(ob12v,1,new TGeoTranslation(0.,0.975,7.45));
   //
   // holes for connectors
   //
   fileName = gSystem->Getenv("ALICE_ROOT");
   fileName += "/TPC/conn_oroc.dat";
   in.open(fileName.Data(), ios_base::in); // asci file
   TGeoRotation *rr[78];
   for(Int_t i =0;i<78;i++){
      Double_t y =3.89;
      Double_t x,z,ang;
      Double_t x1,z1,x2,z2;
      in>>x>>z>>ang;        
      Double_t xr = 4.7*TMath::Sin(ang*TMath::DegToRad());
      Double_t zr = 4.7*TMath::Cos(ang*TMath::DegToRad());
      //
      x1=xr+x; x2=-xr+x; z1=zr+z; z2 = -zr+z;      
      //
      rr[i]= new TGeoRotation();
      rr[i]->RotateY(ang); 
      z1-=54.95;
      z2-=54.95;
      //
      obdv->AddNode(connv,i+1,new TGeoCombiTrans(x1,y,z1,rr[i]));
      obdv->AddNode(connv,i+79,new TGeoCombiTrans(x2,y,z2,rr[i]));
   }
   in.close();
   // cap
   new TGeoTrd1("ocap",23.3874,43.5239,1.09,57.1);
   new TGeoTrd1("ohole",22.2938,40.5084,1.09,51.65);
   TGeoTranslation *tr5 = new TGeoTranslation("tr5",0.,0.,-2.15);
   tr5->RegisterYourself();
   TGeoCompositeShape *oc = new TGeoCompositeShape("ocap-ohole:tr5");
   TGeoVolume *ocv = new TGeoVolume("TPC_ORCAP",oc,m3);
   //
   // pad plane and wire fixations
   //
   TGeoTrd1 *opp = new TGeoTrd1(23.3874,43.5239,0.3,57.1);
   TGeoVolume *oppv = new TGeoVolume("TPC_ORPP",opp,m4);
   //
   tga=(43.5239-23.3874)/114.2;
   TGeoPara *f3 = new TGeoPara(.7,.6,57.1,0.,-10.,0.);
   TGeoPara *f4 = new TGeoPara(.7,.6,57.1,0.,10.,0.);  
   xx = 57.1*tga+23.3874-0.7;
   TGeoVolume *f3v = new TGeoVolume("TPC_ORF1",f3,m4);
   TGeoVolume *f4v = new TGeoVolume("TPC_ORF2",f4,m4);
   //
   TGeoVolumeAssembly *oroc = new TGeoVolumeAssembly("TPC_OROC");
   //
   oroc->AddNode(obdv,1);
   oroc->AddNode(ocv,1,new TGeoTranslation(0.,3.1,2.15));
   oroc->AddNode(oppv,1,new TGeoTranslation(0.,4.49,2.15));
   oroc->AddNode(f3v,1,new TGeoTranslation(-xx,5.39,2.15));
   oroc->AddNode(f4v,1,new TGeoTranslation(xx,5.39,2.15));
   // 
   // now iroc and oroc are placed into a sector...
   //
   TGeoVolumeAssembly *secta = new TGeoVolumeAssembly("TPC_SECT"); // a-side
   TGeoVolumeAssembly *sectc = new TGeoVolumeAssembly("TPC_SECT"); // c-side
   TGeoRotation rot1("rot1",90.,90.,0.);
   TGeoRotation rot2("rot2");
   rot2.RotateY(10.);
   TGeoRotation *rot = new TGeoRotation("rot");
   *rot=rot1*rot2;
   //
   Double_t x0,y0;
   x0=110.2*TMath::Cos(openingAngle);
   y0=110.2*TMath::Sin(openingAngle);
   TGeoCombiTrans *combi1a = new TGeoCombiTrans("combi1",x0,y0,1.09+0.195,rot); //a-side 
   TGeoCombiTrans *combi1c = new TGeoCombiTrans("combi1",x0,y0,1.09+0.222,rot); //c-side
   x0=188.45*TMath::Cos(openingAngle);
   y0=188.45*TMath::Sin(openingAngle);
   TGeoCombiTrans *combi2a = new TGeoCombiTrans("combi2",x0,y0,0.99+0.195,rot); //a-side
   TGeoCombiTrans *combi2c = new TGeoCombiTrans("combi2",x0,y0,0.99+0.222,rot); //c-side
   //
   //
   // A-side
   //
   secta->AddNode(ch,1);
   secta->AddNode(iroc,1,combi1a);
   secta->AddNode(oroc,1,combi2a);
   //
   // C-side
   //
   sectc->AddNode(ch,1);
   sectc->AddNode(iroc,1,combi1c);
   sectc->AddNode(oroc,1,combi2c);
   //
   // now I try to make  wheels...
   //
   TGeoVolumeAssembly *wheela = new TGeoVolumeAssembly("TPC_ENDCAP");
   TGeoVolumeAssembly *wheelc = new TGeoVolumeAssembly("TPC_ENDCAP");
   //
   TGeoRotation *rwh[18]; 
   for(Int_t i =0;i<18;i++){
     Double_t phi = (20.*i);
     rwh[i]=new TGeoRotation();
     rwh[i]->RotateZ(phi);
     wheela->AddNode(secta,i+1,rwh[i]);
     wheelc->AddNode(sectc,i+1,rwh[i]); 
    
   }
   // wheels in the drift volume!   

   TGeoCombiTrans *combi3 = new TGeoCombiTrans("combi3",0.,0.,256.6,ref);
   v9->AddNode(wheela,1,combi3);
   v9->AddNode(wheelc,2,new TGeoTranslation(0.,0.,-256.6));
   //_____________________________________________________________
   // service support wheel
   //_____________________________________________________________
  TGeoPgon *sw = new TGeoPgon(0.,20.,1,2);
  sw->DefineSection(0,-4.,80.5,251.75);
  sw->DefineSection(1,4.,80.5,251.75); 
  TGeoVolume *swv = new TGeoVolume("TPC_SWSEG",sw,m3); //Al
  //
  thick=1.;
  shift = thick/TMath::Sin(openingAngle);
  TGeoPgon *sh = new TGeoPgon(0.,20.,1,2);
  sh->DefineSection(0,-4.,81.5-shift,250.75-shift);
  sh->DefineSection(1,4.,81.5-shift,250.75-shift);
  TGeoVolume *shv = new TGeoVolume("TPC_SWS1",sh,m1); //Air
  //
  TGeoMedium *m9 =  gGeoManager->GetMedium("TPC_Si"); 
  TGeoPgon *el = new TGeoPgon(0.,20.,1,2);
  el->DefineSection(0,-1.872,81.5-shift,250.75-shift);
  el->DefineSection(1,1.872,81.5-shift,250.75-shift);
  TGeoVolume *elv = new TGeoVolume("TPC_ELEC",el,m9); //Si 
  //
  shv->AddNode(elv,1);
  //
  //
  ys = shift*TMath::Sin(openingAngle);
  xs = shift*TMath::Cos(openingAngle);
  swv->AddNode(shv,1,new TGeoTranslation(xs,ys,0.));
  // cover
  TGeoPgon *co = new TGeoPgon(0.,20.,1,2);
  co->DefineSection(0,-0.5,77.,255.25);
  co->DefineSection(1,0.5,77.,255.25);
  TGeoVolume *cov = new TGeoVolume("TPC_SWC1",co,m3);//Al
  // hole in a cover
  TGeoPgon *coh = new TGeoPgon(0.,20.,1,2);
  shift=4./TMath::Sin(openingAngle);
  coh->DefineSection(0,-0.5,85.-shift,247.25-shift);
  coh->DefineSection(1,0.5,85.-shift,247.25-shift);  
  //
  TGeoVolume *cohv = new TGeoVolume("TPC_SWC2",coh,m1);
  //
  ys = shift*TMath::Sin(openingAngle);
  xs = shift*TMath::Cos(openingAngle);  
  cov->AddNode(cohv,1,new TGeoTranslation(xs,ys,0.));
  //
  // Sector as an Assembly
  //
  TGeoVolumeAssembly *swhs = new TGeoVolumeAssembly("TPC_SSWSEC");
  swhs->AddNode(swv,1);
  swhs->AddNode(cov,1,new TGeoTranslation(0.,0.,-4.5));
  swhs->AddNode(cov,2,new TGeoTranslation(0.,0.,4.5));
  //
  // SSW as an Assembly of sectors
  //
  TGeoRotation *rsw[18];
  TGeoVolumeAssembly *swheel = new TGeoVolumeAssembly("TPC_SSWHEEL");
   for(Int_t i =0;i<18;i++){
     Double_t phi = (20.*i);
     rsw[i] = new TGeoRotation();
     rsw[i]->RotateZ(phi);
     swheel->AddNode(swhs,i+1,rsw[i]);   
   }
   v1->AddNode(swheel,1,new TGeoTranslation(0.,0.,-284.6));
   v1->AddNode(swheel,2,new TGeoTranslation(0.,0.,284.6));

   // sensitive strips - strip "0" is always set
   // conditional
   Int_t totrows;
   totrows = fTPCParam->GetNRowLow() + fTPCParam->GetNRowUp();
   Double_t *upar;
   upar=NULL;
   gGeoManager->Volume("TPC_Strip","PGON",m5->GetId(),upar);
   upar=new Double_t [10];
   upar[0]=0.;
   upar[1]=360.;
   upar[2]=18.;
   upar[3]=2.;
   //
   upar[4]=-124.8;
   upar[7]=124.8;

   Double_t rlow=fTPCParam->GetPadRowRadiiLow(0);

   upar[5]=rlow;
   upar[6]=rlow+.01;
   upar[8]=upar[5];
   upar[9]=upar[6];
   //
   gGeoManager->Node("TPC_Strip",1,"TPC_Drift",0.,0.,124.82,0,kTRUE,upar,10);
   gGeoManager->Node("TPC_Strip",totrows+1,
                     "TPC_Drift",0.,0.,-124.82,0,kTRUE,upar,10);
   //
   // now, strips optionally
   //
   if(fSens){
     //lower sectors
     for(Int_t i=2;i<fTPCParam->GetNRowLow()+1;i++){
       rlow=fTPCParam->GetPadRowRadiiLow(i-1);
       upar[5]=rlow;
       upar[6]=rlow+.01;
       upar[8]=upar[5];
       upar[9]=upar[6];
       gGeoManager->Node("TPC_Strip",i,
                         "TPC_Drift",0.,0.,124.82,0,kTRUE,upar,10);
       gGeoManager->Node("TPC_Strip",totrows+i,
                         "TPC_Drift",0.,0.,-124.82,0,kTRUE,upar,10);       
     }
     //upper sectors
     for(Int_t i=1;i<fTPCParam->GetNRowUp()+1;i++){
       rlow=fTPCParam->GetPadRowRadiiUp(i-1); 
       upar[5]=rlow;
       upar[6]=rlow+.01;
       upar[8]=upar[5];
       upar[9]=upar[6];
       gGeoManager->Node("TPC_Strip",i+fTPCParam->GetNRowLow(),
                         "TPC_Drift",0.,0.,124.82,0,kTRUE,upar,10); 
       gGeoManager->Node("TPC_Strip",totrows+i+fTPCParam->GetNRowLow(),
                         "TPC_Drift",0.,0.,-124.82,0,kTRUE,upar,10); 
     }
   }//strips
  //----------------------------------------------------------
  // TPC Support Rods - MAKROLON
  //----------------------------------------------------------
  TGeoMedium *m6=gGeoManager->GetMedium("TPC_Makrolon");
  TGeoMedium *m7=gGeoManager->GetMedium("TPC_Cu");
  TGeoMedium *m10 =  gGeoManager->GetMedium("TPC_Alumina");
  TGeoMedium *m11 =  gGeoManager->GetMedium("TPC_Peek");
  TGeoMedium *m12 =  gGeoManager->GetMedium("TPC_Water");
  TGeoMedium *m13 = gGeoManager->GetMedium("TPC_Brass");
  // 
  // tpc rod is an assembly of 10 long parts and 2 short parts
  // connected with alu rings and plagged on both sides.
  //
  //
// tpc rod long
//
  TGeoPcon *rod = new TGeoPcon("rod",0.,360.,6);
 rod->DefineSection(0,-10.43,1.92,2.08);
 rod->DefineSection(1,-9.75,1.92,2.08);

 rod->DefineSection(2,-9.75,1.8,2.2);
 rod->DefineSection(3,9.75,1.8,2.2);

 rod->DefineSection(4,9.75,1.92,2.08);
 rod->DefineSection(5,10.43,1.92,2.08);
 //
 TGeoVolume *mrodl = new TGeoVolume("TPC_mrodl",rod,m6);
 //
 // tpc rod short 
 //
 TGeoPcon *rod1 = new TGeoPcon("rod1",0.,360.,6);
 rod1->DefineSection(0,-8.93,1.92,2.08);
 rod1->DefineSection(1,-8.25,1.92,2.08);

 rod1->DefineSection(2,-8.25,1.8,2.2);
 rod1->DefineSection(3,8.25,1.8,2.2);

 rod1->DefineSection(4,8.25,1.92,2.08);
 rod1->DefineSection(5,8.93,1.92,2.08);
 //
 TGeoVolume *mrods = new TGeoVolume("TPC_mrods",rod1,m6);
 //
 // below is for the resistor rod
 //
 // hole for the brass connectors
 //

 new TGeoTube("hhole",0.,0.3,0.3);
 //
 //transformations for holes - initialy they
 // are placed at x=0 and negative y
 //
  TGeoRotation *rhole = new TGeoRotation();
  rhole->RotateX(90.);
  TGeoCombiTrans *transf[13];
  Char_t name[30];
  for(Int_t i=0;i<13;i++){
    //sprintf(name,"transf%d",i);
    snprintf(name,30,"transf%d",i);
    transf[i]= new TGeoCombiTrans(name,0.,-2.,-9.+i*1.5,rhole);
    transf[i]->RegisterYourself();
  }
  // union expression for holes
  TString operl("hhole:transf0");
  for (Int_t i=1;i<13;i++){
    //sprintf(name,"+hhole:transf%d",i);
    snprintf(name,30,"+hhole:transf%d",i);
    operl.Append(name);   
  }
  //
 TString opers("hhole:transf1");
  for (Int_t i=2;i<12;i++){
    //sprintf(name,"+hhole:transf%d",i);
    snprintf(name,30,"+hhole:transf%d",i); 
    opers.Append(name);   
  }
  //union of holes
  new TGeoCompositeShape("hlv",operl.Data());
  new TGeoCompositeShape("hsv",opers.Data());
  //
  TGeoCompositeShape *rodl = new TGeoCompositeShape("rodl","rod-hlv");
  TGeoCompositeShape *rods = new TGeoCompositeShape("rods","rod1-hsv");
 //rods - volumes - makrolon rods with holes
  TGeoVolume *rodlv = new TGeoVolume("TPC_rodl",rodl,m6);
  TGeoVolume *rodsv = new TGeoVolume("TPC_rods",rods,m6);
  //brass connectors
 //connectors
  TGeoTube *bcon = new TGeoTube(0.,0.3,0.3);//connectors
  TGeoVolume *bconv = new TGeoVolume("TPC_bcon",bcon,m13);
 //
 // hooks holding strips
 //
 new TGeoBBox("hk1",0.625,0.015,0.75);
 new TGeoBBox("hk2",0.625,0.015,0.15);
 TGeoTranslation *tr21 = new TGeoTranslation("tr21",0.,-0.03,-0.6);
 TGeoTranslation *tr12 = new TGeoTranslation("tr12",0.,-0.03,0.6);
 tr21->RegisterYourself();
 tr12->RegisterYourself();
 
 TGeoCompositeShape *hook = new TGeoCompositeShape("hook","hk1+hk2:tr21+hk2:tr12");
 TGeoVolume *hookv = new TGeoVolume("TPC_hook",hook,m13);
 //
 // assembly of the short rod with connectors and hooks
 //
 //
 // short rod
 //
  TGeoVolumeAssembly *spart = new TGeoVolumeAssembly("TPC_spart");
  //
  spart->AddNode( rodsv,1);
  for(Int_t i=1;i<12;i++){
  spart->AddNode(bconv,i,transf[i]);
  }
  for(Int_t i =0;i<11;i++){
    spart->AddNode(hookv,i+1,new TGeoTranslation(0.,-2.315,-7.5+i*1.5));
  } 
 //
 // long rod
 //
  TGeoVolumeAssembly *lpart = new TGeoVolumeAssembly("TPC_lpart");
  //
  lpart->AddNode( rodlv,1);
  for(Int_t i=0;i<13;i++){
  lpart->AddNode(bconv,i+12,transf[i]);
  }
  for(Int_t i =0;i<13;i++){
    lpart->AddNode(hookv,i+12,new TGeoTranslation(0.,-2.315,-9.+i*1.5));
  }   
  //
  // alu ring
  //
  new TGeoTube("ring1",2.1075,2.235,0.53);
  new TGeoTube("ring2",1.7925,1.89,0.43);
  new TGeoTube("ring3",1.89,2.1075,0.05);
  TGeoCompositeShape *ring = new TGeoCompositeShape("ring","ring1+ring2+ring3");
  TGeoVolume *ringv = new TGeoVolume("TPC_ring",ring,m3);
  //
  // rod assembly
  //
  TGeoVolumeAssembly *tpcrrod = new TGeoVolumeAssembly("TPC_rrod");//rrod
  TGeoVolumeAssembly *tpcmrod = new TGeoVolumeAssembly("TPC_mrod");//makrolon rod  
  //long pieces
  for(Int_t i=0;i<11;i++){
    tpcrrod->AddNode(ringv,i+1,new TGeoTranslation(0.,0.,-105.+i*21));
    tpcmrod->AddNode(ringv,i+12,new TGeoTranslation(0.,0.,-105.+i*21));
  }
  for(Int_t i=0;i<10;i++){
    tpcrrod->AddNode(lpart,i+1,new TGeoTranslation(0.,0.,-94.5+i*21));
    tpcmrod->AddNode(mrodl,i+1,new TGeoTranslation(0.,0.,-94.5+i*21));	    
  }
  //
  // short pieces
  //
  tpcrrod->AddNode(spart,1,new TGeoTranslation(0.,0.,-114.));
  tpcrrod->AddNode(spart,2,new TGeoTranslation(0.,0.,114.)); 
  tpcrrod->AddNode(ringv,23,new TGeoTranslation(0.,0.,-123.));
  tpcrrod->AddNode(ringv,24,new TGeoTranslation(0.,0.,123.));
  //
  tpcmrod->AddNode(mrods,1,new TGeoTranslation(0.,0.,-114.));
  tpcmrod->AddNode(mrods,2,new TGeoTranslation(0.,0.,114.)); 
  tpcmrod->AddNode(ringv,25,new TGeoTranslation(0.,0.,-123.));
  tpcmrod->AddNode(ringv,26,new TGeoTranslation(0.,0.,123.));
 //
  // left plaque
  //
  TGeoPcon *lp = new TGeoPcon(0.,360.,4);
  //
  lp->DefineSection(0,-125.8,1.92,2.235);
  lp->DefineSection(1,-124.8,1.92,2.235);
  //
  lp->DefineSection(2,-124.8,1.92,2.08);
  lp->DefineSection(3,-123.1,1.92,2.08);
  //
  TGeoVolume *lpv1 = new TGeoVolume("TPC_lpv1",lp,m6);
  TGeoVolume *lpv2 = new TGeoVolume("TPC_lpv2",lp,m6);
  // left ring
  TGeoTube *lr = new TGeoTube(2.1075,2.235,0.5);
  TGeoVolume *lrv = new TGeoVolume("TPC_lrv",lr,m3);
  //
  lpv2->AddNode(lrv,1,new TGeoTranslation(0.,0.,-125.3));
  //
  tpcrrod->AddNode(lpv2,1); 
  tpcmrod->AddNode(lpv1,1);
  //
  //  right plaque
  //
  TGeoTube *rp = new TGeoTube(1.92,2.08,2.025);
  TGeoVolume *rpv = new TGeoVolume("TPC_rpv",rp,m6);
  tpcrrod->AddNode(rpv,1, new TGeoTranslation(0.,0.,125.125)); 
  tpcmrod->AddNode(rpv,2,new TGeoTranslation(0.,0.,125.125));
  //
  //
  //HV rods - makrolon + 0.58cm (diameter) Cu
  TGeoTube *hvr = new TGeoTube(0.,1.465,126.5);
  TGeoTube *hvc = new TGeoTube(0.,0.29,126.5);
  //
  TGeoVolume *hvrv = new TGeoVolume("TPC_HV_Rod",hvr,m6);
  TGeoVolume *hvcv = new TGeoVolume("TPC_HV_Cable",hvc,m7);
  hvrv->AddNode(hvcv,1);
  //
  //resistor rod
  //
  TGeoTube *cr = new TGeoTube(0.,0.45,126.2);   
  TGeoTube *cw = new TGeoTube(0.,0.15,126.2);  
  TGeoVolume *crv = new TGeoVolume("TPC_CR",cr,m10);
  TGeoVolume *cwv = new TGeoVolume("TPC_W",cw,m12);   
  //
  // ceramic rod with water
  //
  crv->AddNode(cwv,1);
  //
  //peek rod
  //
  TGeoTube *pr =new TGeoTube(0.2,0.35,126.2);  
  TGeoVolume *prv = new TGeoVolume("TPC_PR",pr,m11); 
  //
  // copper plates with connectors
  //
  new TGeoTube("tub",0.,1.7,0.025);
  //
  // half space - points on the plane and a normal vector
  //
  Double_t n[3],p[3];
  Double_t slope = TMath::Tan(22.*TMath::DegToRad());
  Double_t intp = 1.245;
  //
  Double_t b = slope*slope+1.;
  p[0]=intp*slope/b;
  p[1]=-intp/b;
  p[2]=0.;
  //
  n[0]=-p[0];
  n[1]=-p[1];
  n[2]=0.;
  Double_t norm;
  norm=TMath::Sqrt(n[0]*n[0]+n[1]*n[1]);
  n[0] /= norm;
  n[1] /=norm;
  //
  new TGeoHalfSpace("sp1",p,n);
  //
  slope = -slope; 
  //
  p[0]=intp*slope/b;
  p[1]=-intp/b;
  //
  n[0]=-p[0];
  n[1]=-p[1];
  norm=TMath::Sqrt(n[0]*n[0]+n[1]*n[1]);
  n[0] /= norm;
  n[1] /=norm;
  //
  new TGeoHalfSpace("sp2",p,n);
  // holes for rods
 //holes
 new TGeoTube("h1",0.,0.5,0.025);
 new TGeoTube("h2",0.,0.35,0.025);
 //translations:
 TGeoTranslation *ttr11 = new TGeoTranslation("ttr11",-0.866,0.5,0.);
 TGeoTranslation *ttr22 = new TGeoTranslation("ttr22",0.866,0.5,0.);
 ttr11->RegisterYourself();
 ttr22->RegisterYourself();
 // elastic connector
 new TGeoBBox("elcon",0.72,0.005,0.3);
 TGeoRotation *crr1 = new TGeoRotation();
 crr1->RotateZ(-22.);
TGeoCombiTrans *ctr1 = new TGeoCombiTrans("ctr1",-0.36011, -1.09951,-0.325,crr1);
ctr1->RegisterYourself();
 TGeoCompositeShape *cs1 = new TGeoCompositeShape("cs1",
"(((((tub-h1:ttr11)-h1:ttr22)-sp1)-sp2)-h2)+elcon:ctr1");
 //
 TGeoVolume *csvv = new TGeoVolume("TPC_RR_CU",cs1,m7);
 //
 // resistor rod assembly 2 ceramic rods, peak rod, Cu plates
 // and resistors
 //
 TGeoVolumeAssembly *rrod = new TGeoVolumeAssembly("TPC_RRIN");
 // rods
 rrod->AddNode(crv,1,ttr11);
 rrod->AddNode(crv,2,ttr22); 
 rrod->AddNode(prv,1);
 //Cu plates
 for(Int_t i=0;i<165;i++){
   rrod->AddNode(csvv,i+1,new TGeoTranslation(0.,0.,-122.675+i*1.5));
 }
 //resistors
 TGeoTube *res = new TGeoTube(0.,0.15,0.5);
 TGeoVolume *resv = new TGeoVolume("TPC_RES",res,m10);
 TGeoVolumeAssembly *ress = new TGeoVolumeAssembly("TPC_RES_CH");
 ress->AddNode(resv,1,new TGeoTranslation(0.2,0.,0.));
 ress->AddNode(resv,2,new TGeoTranslation(-0.2,0.,0.));
 //
 TGeoRotation *crr2 = new TGeoRotation();
 crr2->RotateY(30.);
 TGeoRotation *crr3 = new TGeoRotation();
 crr3->RotateY(-30.); 
 //
 for(Int_t i=0;i<164;i+=2){
   rrod->AddNode(ress,i+1, new TGeoCombiTrans(0.,1.2,-121.925+i*1.5,crr2));
   rrod->AddNode(ress,i+2, new TGeoCombiTrans(0.,1.2,-121.925+(i+1)*1.5,crr3));
 }

 tpcrrod->AddNode(rrod,1,new TGeoCombiTrans(0.,0.,0.5,crr1));
 //
 // guard ring resistor chain
 //

 TGeoTube *gres1 = new TGeoTube(0.,0.375,125.);// inside ifc
 //
 TGeoVolume *vgres1 = new TGeoVolume("TPC_GRES1",gres1,m10);

 //
 Double_t xrc,yrc;
 //
 xrc=79.3*TMath::Cos(350.*TMath::DegToRad());
 yrc=79.3*TMath::Sin(350.*TMath::DegToRad());
 //
 v9->AddNode(vgres1,1,new TGeoTranslation(xrc,yrc,126.9));
 v9->AddNode(vgres1,2,new TGeoTranslation(xrc,yrc,-126.9));
 //
 xrc=79.3*TMath::Cos(190.*TMath::DegToRad());
 yrc=79.3*TMath::Sin(190.*TMath::DegToRad()); 
 //
 v9->AddNode(vgres1,3,new TGeoTranslation(xrc,yrc,126.9));
 v9->AddNode(vgres1,4,new TGeoTranslation(xrc,yrc,-126.9));
 //------------------------------------------------------------------
 TGeoRotation refl("refl",90.,0.,90.,90.,180.,0.);
 TGeoRotation rotrod("rotrod");
 //
 TGeoRotation *rotpos[2]; 
 //
 TGeoRotation *rotrod1[2]; 
 TGeoTubeSeg *irh = new TGeoTubeSeg(78.825,79.25,1.5,358.5,1.5);
 TGeoTubeSeg *orh = new TGeoTubeSeg(256.5,257.95,1.5,359.5,0.5);
 TGeoTubeSeg *ohh = new TGeoTubeSeg(256.5,257.95,1.5,9.5,10.5);
  TGeoVolume *irhv = new TGeoVolume("TPC_IRHH",irh,m4);
  TGeoVolume *orhv = new TGeoVolume("TPC_ORHH",orh,m4);
  TGeoVolume *ohhv = new TGeoVolume("TPC_OHVHH",ohh,m4);
  
 //v9 - drift gas

  for(Int_t i=0;i<18;i++){
    Double_t angle,x,y;
    Double_t z,r; 
    angle=TMath::DegToRad()*20.*(Double_t)i;
    TGeoRotation *roth = new TGeoRotation(); //rotation for rod holders
    roth->RotateZ(angle);
    //inner rods
    r=81.5;
    x=r * TMath::Cos(angle);
    y=r * TMath::Sin(angle);
    z = 126.1;
    //
    v9->AddNode(irhv,i+1,roth);
    v9->AddNode(orhv,i+1,roth);
    v9->AddNode(ohhv,i+1,roth);
    //
    if(i==11){//resistor rod inner
       rotrod.RotateZ(-90.+angle);
       rotrod1[0]= new TGeoRotation();
       rotpos[0]= new TGeoRotation();
       //
       rotrod1[0]->RotateZ(-90.+angle);
       *rotpos[0] = refl*rotrod; //rotation+reflection
	v9->AddNode(tpcrrod,1,new TGeoCombiTrans(x,y, z, rotrod1[0])); //A
	v9->AddNode(tpcrrod,2,new TGeoCombiTrans(x,y,-z, rotpos[0])); //C      
    } 
    else { 
      v9->AddNode(tpcmrod,i+1,new TGeoTranslation(x,y,z));//shaft
      v9->AddNode(tpcmrod,i+19,new TGeoCombiTrans(x,y,-z,ref));//muon
    }
    // outer rods
    r=254.25;
    x=r * TMath::Cos(angle);
    y=r * TMath::Sin(angle);
    z=126.3;
    //
    if(i==3){//resistor rod outer
      rotrod.RotateZ(90.+angle);
      rotrod1[1]= new TGeoRotation();
      rotpos[1]= new TGeoRotation();
      rotrod1[1]->RotateZ(90.+angle);
      *rotpos[1] = refl*rotrod;//rotation+reflection
      v9->AddNode(tpcrrod,3,new TGeoCombiTrans(x,y, z, rotrod1[1])); //A 
      v9->AddNode(tpcrrod,4,new TGeoCombiTrans(x,y, -z, rotpos[1])); //C
    }
    else {
      v9->AddNode(tpcmrod,i+37,new TGeoTranslation(x,y,z));//shaft
      v9->AddNode(tpcmrod,i+55,new TGeoCombiTrans(x,y,-z,ref));//muon      
    }
    if(i==15){
      v9->AddNode(hvrv,1,new TGeoTranslation(x,y,z+0.7)); //hv->A-side only      
    }
  } //end of rods positioning

  TGeoVolume *alice = gGeoManager->GetVolume("ALIC");
  alice->AddNode(v1,1);  
  
} // end of function
 
//_____________________________________________________________________________
void AliTPCv2::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  // 
  SetInnerChambersAlignable();
  SetOuterChambersAlignable();
}
 
//_____________________________________________________________________________
void AliTPCv2::SetInnerChambersAlignable() const
{
  //
  AliGeomManager::ELayerID idTPC1 = AliGeomManager::kTPC1;
  Int_t modUID, modnum = 0;
  TString vpstr1 = "ALIC_1/TPC_M_1/TPC_Drift_1/TPC_ENDCAP_1/TPC_SECT_";
  TString vpstr2 = "ALIC_1/TPC_M_1/TPC_Drift_1/TPC_ENDCAP_2/TPC_SECT_";
  TString vpappend = "/TPC_IROC_1";
  TString snstr1="TPC/EndcapA/Sector";
  TString snstr2="TPC/EndcapC/Sector";
  TString snappend="/InnerChamber";
  TString volpath, symname;
  
  for(Int_t cnt=1; cnt<=18; cnt++){
    modUID = AliGeomManager::LayerToVolUID(idTPC1,modnum++);
    volpath = vpstr1;
    volpath += cnt;
    volpath += vpappend;
    symname = snstr1;
    symname += cnt;
    symname += snappend;
    if(!gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data(),modUID))
      AliFatal(Form("Alignable entry %s not created. Volume path %s not valid", symname.Data(),volpath.Data()));
    TGeoPNEntry *alignableEntry = gGeoManager->GetAlignableEntryByUID(modUID);
    TGeoHMatrix* globMatrix = alignableEntry->GetGlobalOrig();
    TGeoHMatrix* matTtoL = fTPCParam->Tracking2LocalMatrix(globMatrix,cnt-1);
    alignableEntry->SetMatrix(matTtoL);
  }

  for(Int_t cnt=1; cnt<=18; cnt++){
    modUID = AliGeomManager::LayerToVolUID(idTPC1,modnum++);
    volpath = vpstr2;
    volpath += cnt;
    volpath += vpappend;
    symname = snstr2;
    symname += cnt;
    symname += snappend;
    if(!gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data(),modUID))
      AliFatal(Form("Alignable entry %s not created. Volume path %s not valid", symname.Data(),volpath.Data()));
    TGeoPNEntry *alignableEntry = gGeoManager->GetAlignableEntryByUID(modUID);
    TGeoHMatrix* globMatrix = alignableEntry->GetGlobalOrig();
    TGeoHMatrix* matTtoL = fTPCParam->Tracking2LocalMatrix(globMatrix,18+cnt-1);
    alignableEntry->SetMatrix(matTtoL);
  }
}

//_____________________________________________________________________________
void AliTPCv2::SetOuterChambersAlignable() const
{
  //
  AliGeomManager::ELayerID idTPC2 = AliGeomManager::kTPC2;
  Int_t modUID, modnum = 0;
  TString vpstr1 = "ALIC_1/TPC_M_1/TPC_Drift_1/TPC_ENDCAP_1/TPC_SECT_";
  TString vpstr2 = "ALIC_1/TPC_M_1/TPC_Drift_1/TPC_ENDCAP_2/TPC_SECT_";
  TString vpappend = "/TPC_OROC_1";
  TString snstr1="TPC/EndcapA/Sector";
  TString snstr2="TPC/EndcapC/Sector";
  TString snappend="/OuterChamber";
  TString volpath, symname;
  
  for(Int_t cnt=1; cnt<=18; cnt++){
    modUID = AliGeomManager::LayerToVolUID(idTPC2,modnum++);
    volpath = vpstr1;
    volpath += cnt;
    volpath += vpappend;
    symname = snstr1;
    symname += cnt;
    symname += snappend;
    if(!gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data(),modUID))
      AliFatal(Form("Alignable entry %s not created. Volume path %s not valid", symname.Data(),volpath.Data()));
    TGeoPNEntry *alignableEntry = gGeoManager->GetAlignableEntryByUID(modUID);
    TGeoHMatrix* globMatrix = alignableEntry->GetGlobalOrig();
    TGeoHMatrix* matTtoL = fTPCParam->Tracking2LocalMatrix(globMatrix,36+cnt-1);
    alignableEntry->SetMatrix(matTtoL);
  }

  for(Int_t cnt=1; cnt<=18; cnt++){
    modUID = AliGeomManager::LayerToVolUID(idTPC2,modnum++);
    volpath = vpstr2;
    volpath += cnt;
    volpath += vpappend;
    symname = snstr2;
    symname += cnt;
    symname += snappend;
    if(!gGeoManager->SetAlignableEntry(symname.Data(),volpath.Data(),modUID))
      AliFatal(Form("Alignable entry %s not created. Volume path %s not valid", symname.Data(),volpath.Data()));
    TGeoPNEntry *alignableEntry = gGeoManager->GetAlignableEntryByUID(modUID);
    TGeoHMatrix* globMatrix = alignableEntry->GetGlobalOrig();
    TGeoHMatrix* matTtoL = fTPCParam->Tracking2LocalMatrix(globMatrix,36+18+cnt-1);
    alignableEntry->SetMatrix(matTtoL);
  }
}


//_____________________________________________________________________________
void AliTPCv2::CreateMaterials()
{
  //
  // Define materials for version 2 of the Time Projection Chamber
  //
 
  AliTPC::CreateMaterials();
}

//_____________________________________________________________________________
void AliTPCv2::Init()
{
  //
  // Initialises version 2 of the TPC after that it has been built
  //

  Int_t *idtmed = fIdtmed->GetArray();
  
  AliTPC::Init();

 
  fIdSens=gMC->VolId("TPC_Strip");  // one strip is always selected...

  fIDrift=gMC->VolId("TPC_Drift");
  fSecOld=-100; // fake number 

  gMC->SetMaxNStep(-30000); // max. number of steps increased

  if (fPrimaryIonisation) {
    // for FLUKA
      gMC->Gstpar(idtmed[2],"PRIMIO_E", 20.77); // 1st ionisation potential
 
      gMC->Gstpar(idtmed[2],"PRIMIO_N", 14.35);
      gMC->Gstpar(idtmed[2],"LOSS", 14); // specific energy loss
      gMC->Gstpar(idtmed[2],"STRA",4);
  } 
  // specific energy loss for geant3 is now defined in galice.cuts


  AliDebug(1,"*** TPC version 2 initialized ***");
  AliDebug(1,Form("Maximum number of steps = %d",gMC->GetMaxNStep()));

  //
  
}

//_____________________________________________________________________________
void AliTPCv2::StepManager()
{
  //
  // Called for every step in the Time Projection Chamber
  //

  //
  // parameters used for the energy loss calculations
  //
  const Float_t kprim = 14.35; // number of primary collisions per 1 cm
  const Float_t kpoti = 20.77e-9; // first ionization potential for Ne/CO2
  const Float_t kwIon = 35.97e-9; // energy for the ion-electron pair creation 
  const Int_t   kMaxDistRef =15;     // maximal difference between 2 stored references 
 
  const Float_t kbig = 1.e10;

  Int_t id,copy;
  Float_t hits[5];
  Int_t vol[2];  
  TLorentzVector p;
  
  vol[1]=0; // preset row number to 0
  //
  if (!fPrimaryIonisation) gMC->SetMaxStep(kbig);
  
  if(!gMC->IsTrackAlive()) return; // particle has disappeared
  
  Float_t charge = gMC->TrackCharge();
  
  if(TMath::Abs(charge)<=0.) return; // take only charged particles
  
  // check the sensitive volume

  id = gMC->CurrentVolID(copy); // vol ID and copy number (starts from 1!)
  if(id != fIDrift && id != fIdSens) return; // not in the sensitive folume 

  if ( fPrimaryIonisation && id == fIDrift ) {
    Double_t rnd = gMC->GetRandom()->Rndm();
    gMC->SetMaxStep(0.2+(2.*rnd-1.)*0.05);  // 2 mm +- rndm*0.5mm step
  }   

  //if ( fPrimaryIonisation && id == fIDrift && gMC->IsTrackEntering()) {
  //  gMC->SetMaxStep(0.2);  // 2 mm 
  //}   
  
  gMC->TrackPosition(p);
  Double_t r = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
  //
  
  //
  Double_t angle = TMath::ACos(p[0]/r);  
  angle = (p[1]<0.) ? TMath::TwoPi()-angle : angle;
  //
  // angular segment, it is not a real sector number...
  //
  Int_t sector=TMath::Nint((angle-fTPCParam->GetInnerAngleShift())/
               fTPCParam->GetInnerAngle());
  // rotate to segment "0"
  Float_t cos,sin;
  fTPCParam->AdjustCosSin(sector,cos,sin);
  Float_t x1=p[0]*cos + p[1]*sin;
  // check if within sector's limits
  if((x1>=fTPCParam->GetInnerRadiusLow()&&x1<=fTPCParam->GetInnerRadiusUp())
     ||(x1>=fTPCParam->GetOuterRadiusLow()&&x1<=fTPCParam->GetOuterRadiusUp())){
  // calculate real sector number...
  if (x1>fTPCParam->GetOuterRadiusLow()){
    sector = TMath::Nint((angle-fTPCParam->GetOuterAngleShift())/
             fTPCParam->GetOuterAngle())+fTPCParam->GetNInnerSector();
    if (p[2]<0) 	sector+=(fTPCParam->GetNOuterSector()>>1);
  }
    else   
      if (p[2]<0) sector+=(fTPCParam->GetNInnerSector()>>1);  
  //
  // here I have a sector number
  //

  vol[0]=sector;

  static Double_t lastReferenceR=0;
  if (TMath::Abs(lastReferenceR-r)>kMaxDistRef){
    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kTPC);
    lastReferenceR = r;
  }

  // check if change of sector
  if(sector != fSecOld){
    fSecOld=sector;
    // add track reference
    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kTPC);
  }  
  // track is in the sensitive strip
  if(id == fIdSens){
    // track is entering the strip
    if (gMC->IsTrackEntering()){
      Int_t totrows = fTPCParam->GetNRowLow()+fTPCParam->GetNRowUp();
      vol[1] = (copy<=totrows) ? copy-1 : copy-1-totrows;
      // row numbers are autonomous for lower and upper sectors
      if(vol[0] > fTPCParam->GetNInnerSector()) {
        vol[1] -= fTPCParam->GetNRowLow();
      }
    //
      if(vol[0]<fTPCParam->GetNInnerSector()&&vol[1] == 0){
  
        // lower sector, row 0, because Jouri wants to have this

        gMC->TrackMomentum(p);
        hits[0]=p[0];
        hits[1]=p[1];
        hits[2]=p[2];
        hits[3]=0.; // this hit has no energy loss
        // Get also the track time for pileup simulation
        hits[4]=gMC->TrackTime();

        AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol,hits);  
      }
    //

       gMC->TrackPosition(p);
       hits[0]=p[0];
       hits[1]=p[1];
       hits[2]=p[2];
       hits[3]=0.; // this hit has no energy loss
       // Get also the track time for pileup simulation
       hits[4]=gMC->TrackTime();

       AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol,hits);  
    
    }
    else return;
  }  
  //-----------------------------------------------------------------
  //  charged particle is in the sensitive drift volume
  //-----------------------------------------------------------------
  if(gMC->TrackStep() > 0) {
    Int_t nel=0;
    if (!fPrimaryIonisation) {
      nel = (Int_t)(((gMC->Edep())-kpoti)/kwIon) + 1;
    }
    else {
      static Double_t deForNextStep = 0.;
      // Geant4 (the meaning of Edep as in Geant3) - wrong
      //nel = (Int_t)(((gMC->Edep())-kpoti)/kwIon) + 1;

      // Geant4 (the meaning of Edep as in Geant3) - NEW
      Double_t eAvailable = gMC->Edep() + deForNextStep;
      nel = (Int_t)(eAvailable/kwIon);
      deForNextStep = eAvailable - nel*kwIon;
    }
    nel=TMath::Min(nel,300); // 300 electrons corresponds to 10 keV
    //
    gMC->TrackPosition(p);
    hits[0]=p[0];
    hits[1]=p[1];
    hits[2]=p[2];
    hits[3]=(Float_t)nel;

    // Add this hit

    //    if (fHitType&&2){
    if(fHitType){
      gMC->TrackMomentum(p);
      Float_t momentum = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
      Float_t precision =   (momentum>0.1) ? 0.002 :0.01;
      fTrackHits->SetHitPrecision(precision);
    }

    // Get also the track time for pileup simulation
    hits[4]=gMC->TrackTime();
 
    AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(), vol,hits);
    if (fDebugStreamer){   
      // You can dump here what you need
      // function  CreateDebugStremer() to be called in the Config.C  macro
      // if you want to enable it
      // By default debug streaemer is OFF
      Float_t edep = gMC->Edep();
      Float_t tstep = gMC->TrackStep();
      Int_t pid=gMC->TrackPid();
      (*fDebugStreamer)<<"hit"<<      
	"x="<<hits[0]<<  // hit position
 	"y="<<hits[1]<<
 	"z="<<hits[2]<<
 	"nel="<<hits[3]<<  // number of electorns
 	"tof="<<hits[4]<<  // hit TOF
	"edep="<<edep<<    // energy deposit
	"pid="<<pid<<      // pid
	"step="<<tstep<<
  	"p.="<<&p<<
	"\n";
    }
    
  } // step>0 
  } //within sector's limits
  // Stemax calculation for the next step
  
  Float_t pp;
  TLorentzVector mom;
  // below is valid only for Geant3 (fPromaryIonisation not set)
  if(!fPrimaryIonisation){
    gMC->TrackMomentum(mom);
    Float_t ptot=mom.Rho();
    Float_t betaGamma = ptot/gMC->TrackMass();

    Int_t pid=gMC->TrackPid();
    if((pid==kElectron || pid==kPositron) && ptot > 0.002)
      { 
        pp = kprim*1.58; // electrons above 20 MeV/c are on the plateau!
      }
    else
      {

        betaGamma = TMath::Max(betaGamma,(Float_t)7.e-3); // protection against too small bg
        pp=kprim*AliMathBase::BetheBlochAleph(betaGamma); 
   
    }
  
    Double_t rnd = gMC->GetRandom()->Rndm();
  
    gMC->SetMaxStep(-TMath::Log(rnd)/pp);
  }
  
}


