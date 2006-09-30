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
//  Time Projection Chamber version 0 -- "coarse" TPC                        //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTPCv0Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>

#include <TGeometry.h>
#include <TMath.h>
#include <TNode.h>
#include <TTUBE.h>
#include <TVirtualMC.h>
#include <TString.h>
#include <TSystem.h>

#include "AliConst.h"
#include "AliRun.h"
#include "AliTPCDigitsArray.h"
#include "AliTPCParamSR.h"
#include "AliTPCv0.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoPcon.h"
#include "TGeoTube.h"
#include "TGeoPgon.h"
#include "TGeoTrd1.h"
#include "TGeoCompositeShape.h"
#include "TGeoPara.h"
ClassImp(AliTPCv0)
 
//_____________________________________________________________________________
AliTPCv0::AliTPCv0(const char *name, const char *title) 
         :AliTPC(name, title)
{
  //
  // Standard creator for TPC version 0
  //


  if (fTPCParam)
     fTPCParam->Write(fTPCParam->GetTitle());
}

//_____________________________________________________________________________
void AliTPCv0::CreateGeometry()
{

  //
  // Creation of the TPC version 0, i.e. no sensitive volumes,
  // only the material distribution
  //
  //Begin_Html
  /*
    <img src="picts/AliTPC.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliTPCv0Tree.gif">
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
  TGeoPcon *tpc = new TGeoPcon(0.,360.,18); //18 sections
  tpc->DefineSection(0,-290.,77.,278.);
  tpc->DefineSection(1,-259.6,70.,278.);
  //
  tpc->DefineSection(2,-259.6,68.1,278.);
  tpc->DefineSection(3,-253.6,68.1,278.);
  //
  tpc->DefineSection(4,-253.6,68.,278.);
  tpc->DefineSection(5,-74.0,60.8,278.);
  //
  tpc->DefineSection(6,-74.0,60.1,278.);
  tpc->DefineSection(7,-73.3,60.1,278.);
  //
  tpc->DefineSection(8,-73.3,56.9,278.); 
  tpc->DefineSection(9,73.3,56.9,278.);
  //
  tpc->DefineSection(10,73.3,60.1,278.);
  tpc->DefineSection(11,74.0,60.1,278.);
  //
  tpc->DefineSection(12,74.0,60.8,278.);
  tpc->DefineSection(13,253.6,65.5,278.);
  //
  tpc->DefineSection(14,253.6,65.6,278.);
  tpc->DefineSection(15,259.6,65.6,278.);
  //
  tpc->DefineSection(16,259.6,70.0,278.);
  tpc->DefineSection(17,290.,77.,278.);
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
  TGeoMedium *m5 = gGeoManager->GetMedium("TPC_Ne-CO2-N-1");
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
  TGeoTube *to1 = new TGeoTube(274.8174,277.995,252.1); //epoxy
  TGeoTube *to2 = new TGeoTube(274.8274,277.985,252.1); //tedlar
  TGeoTube *to3 = new TGeoTube(274.8312,277.9812,252.1);//prepreg2
  TGeoTube *to4 = new TGeoTube(274.9062,277.9062,252.1);//nomex
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
  TGeoTube *tf1 = new TGeoTube(258.0,260.0676,252.1); //tedlar
  TGeoTube *tf2 = new TGeoTube(258.0038,260.0638,252.1); //prepreg3
  TGeoTube *tf3 = new TGeoTube(258.0338,260.0338,252.1);//nomex
  //
  TGeoMedium *sm5 = gGeoManager->GetMedium("TPC_Prepreg3");
  //
  TGeoVolume *tf1v = new TGeoVolume("TPC_OFC1",tf1,sm2);
  TGeoVolume *tf2v = new TGeoVolume("TPC_OFC2",tf2,sm5);
  TGeoVolume *tf3v = new TGeoVolume("TPC_OFC3",tf3,sm4);
  //
  // outer part - positioning
  //
  tov1->AddNode(tov2,1); tov2->AddNode(tov3,1); tov3->AddNode(tov4,1);
  //
  tf1v->AddNode(tf2v,1); tf2v->AddNode(tf3v,1);
  //
  v3->AddNode(tov1,1,new TGeoTranslation(0.,0.,-1.5)); v4->AddNode(tf1v,1);
  //
  v2->AddNode(v3,1); v2->AddNode(v4,1); 
  //
  v1->AddNode(v2,1);
  //-----------------------------------------------------------
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
  TGeoTube *t2 = new TGeoTube(76.6774,78.845,74.175); // tedlar 38 microns
  TGeoTube *t3 = new TGeoTube(76.6812,78.8412,74.175); // prepreg2 500 microns
  TGeoTube *t4 = new TGeoTube(76.7312,78.7912,74.175); // prepreg3 300 microns
  TGeoTube *t5 = new TGeoTube(76.7612,78.7612,74.175); // nomex 2 cm
  //
  TGeoVolume *tv2 = new TGeoVolume("TPC_IFC3",t2,sm2);
  TGeoVolume *tv3 = new TGeoVolume("TPC_IFC4",t3,sm3);
  TGeoVolume *tv4 = new TGeoVolume("TPC_IFC5",t4,sm5);
  TGeoVolume *tv5 = new TGeoVolume("TPC_IFC6",t5,sm4);
  //
  // middle parts - 2 copies
  TGeoTube *t6 = new TGeoTube(76.6774,78.795,5.); // tedlar 38 microns
  TGeoTube *t7 = new TGeoTube(76.6812,78.7912,5.); // prepreg2 250 microns
  TGeoTube *t8 = new TGeoTube(76.7062,78.7662,5.); // prepreg3 300 microns
  TGeoTube *t9 = new TGeoTube(76.7362,78.7362,5.); // nomex 2 cm
  //
  TGeoVolume *tv6 = new TGeoVolume("TPC_IFC7",t6,sm2);
  TGeoVolume *tv7 = new TGeoVolume("TPC_IFC8",t7,sm3);
  TGeoVolume *tv8 = new TGeoVolume("TPC_IFC9",t8,sm5);
  TGeoVolume *tv9 = new TGeoVolume("TPC_IFC10",t9,sm4);
  // central part - 1 copy
  TGeoTube *t10 = new TGeoTube(76.6774,78.745,93.75); // tedlar 38 microns 
  TGeoTube *t11 = new TGeoTube(76.6812,78.7412,93.75); // prepreg3 300 microns
  TGeoTube *t12 = new TGeoTube(76.7112,78.7112,93.75); // nomex 2 cm
  //
  TGeoVolume *tv10 = new TGeoVolume("TPC_IFC11",t10,sm2);
  TGeoVolume *tv11 = new TGeoVolume("TPC_IFC12",t11,sm5);
  TGeoVolume *tv12 = new TGeoVolume("TPC_IFC13",t12,sm4);
  //
  // inner part - positioning
  //
  // creating a sandwich
  tv2->AddNode(tv3,1); tv3->AddNode(tv4,1); tv4->AddNode(tv5,1);
  //
  tv6->AddNode(tv7,1); tv7->AddNode(tv8,1); tv8->AddNode(tv9,1);
  //
  tv10->AddNode(tv11,1); tv11->AddNode(tv12,1);
  //
  TGeoVolumeAssembly *tv100 = new TGeoVolumeAssembly("TPC_IFC");
  //
  tv100->AddNode(tv10,1);
  tv100->AddNode(tv6,1,new TGeoTranslation(0.,0.,-98.75));
  tv100->AddNode(tv6,2,new TGeoTranslation(0.,0.,98.75));
  tv100->AddNode(tv2,1,new TGeoTranslation(0.,0.,-177.925));
  tv100->AddNode(tv2,2,new TGeoTranslation(0.,0.,177.925));
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
  cfl->DefineSection(2,-68.6,60.6324,61.2);
  cfl->DefineSection(3,68.6,60.6324,61.2); 
  //
  cfl->DefineSection(4,68.6,59.7,61.2);
  cfl->DefineSection(5,71.1,59.7,61.2);  
  //
  TGeoVolume *cflv = new TGeoVolume("TPC_CDR",cfl,m3);
  // sandwich
  TGeoTube *cd1 = new TGeoTube(60.6424,61.19,71.1);
  TGeoTube *cd2 = new TGeoTube(60.6462,61.1862,71.1);
  TGeoTube *cd3 = new TGeoTube(60.6662,61.1662,71.1);  
  //
  TGeoMedium *sm6 = gGeoManager->GetMedium("TPC_Prepreg1");
  TGeoVolume *cd1v = new TGeoVolume("TPC_CDR1",cd1,sm2); //tedlar
  TGeoVolume *cd2v = new TGeoVolume("TPC_CDR2",cd2,sm6);// prepreg1
  TGeoVolume *cd3v = new TGeoVolume("TPC_CDR3",cd3,sm4); //nomex
  //
  // seals for central drum 2 copies
  //
  TGeoTube *cs = new TGeoTube(56.9,61.2,0.1);
  TGeoMedium *sm7 = gGeoManager->GetMedium("TPC_Mylar");
  TGeoVolume *csv = new TGeoVolume("TPC_CDRS",cs,sm7);
  v1->AddNode(csv,1,new TGeoTranslation(0.,0.,-71.));
  v1->AddNode(csv,2,new TGeoTranslation(0.,0.,71.));
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
  cd1v->AddNode(cd2v,1); cd2v->AddNode(cd3v,1); cflv->AddNode(cd1v,1);
  //
  v1->AddNode(siv,1,new TGeoTranslation(0.,0.,-70.1));
  v1->AddNode(siv,2,new TGeoTranslation(0.,0.,70.1));
  v1->AddNode(sev,1); v1->AddNode(sev,2,ref); v1->AddNode(cflv,1);
  //
  // central membrane - 2 rings and a mylar membrane - assembly
  //
  TGeoTube *ih = new TGeoTube(81.05,84.05,0.3);
  TGeoTube *oh = new TGeoTube(250.,256.,.5);
  TGeoTube *mem = new TGeoTube(84.05,250,0.01);
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
   for(Int_t i =0;i<86;i++){
      Double_t y = 3.99;
      Double_t x,z,ang;
      in>>x>>z>>ang;
      z-=26.5;
      TGeoRotation *rrr = new TGeoRotation();
      rrr->RotateY(ang);
      TGeoCombiTrans *trans = new TGeoCombiTrans("trans",x,y,z,rrr);
      ibdv->AddNode(connv,i+1,trans);
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
   TGeoMedium *m4 = gGeoManager->GetMedium("TPC_G10");
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
      TGeoRotation *rr = new TGeoRotation();
      rr->RotateY(ang); 
      z1-=54.95;
      z2-=54.95;
      TGeoCombiTrans *trans1 = new TGeoCombiTrans("trans1",x1,y,z1,rr);
      TGeoCombiTrans *trans2 = new TGeoCombiTrans("trans2",x2,y,z2,rr);
      obdv->AddNode(connv,i+1,trans1);
      obdv->AddNode(connv,i+79,trans2);
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
   TGeoVolumeAssembly *sect = new TGeoVolumeAssembly("TPC_SECT");
   TGeoRotation rot1("rot1",90.,90.,0.);
   TGeoRotation rot2("rot2");
   rot2.RotateY(10.);
   TGeoRotation *rot = new TGeoRotation("rot");
   *rot=rot1*rot2;
   //
   Double_t x0,y0;
   x0=110.2*TMath::Cos(openingAngle);
   y0=110.2*TMath::Sin(openingAngle);
   TGeoCombiTrans *combi1 = new TGeoCombiTrans("combi1",x0,y0,1.09,rot);
   x0=188.45*TMath::Cos(openingAngle);
   y0=188.45*TMath::Sin(openingAngle);
   TGeoCombiTrans *combi2 = new TGeoCombiTrans("combi2",x0,y0,0.99,rot);
   //
   sect->AddNode(ch,1);
   sect->AddNode(iroc,1,combi1);
   sect->AddNode(oroc,1,combi2);
   //
   // segment is ready...
   // now I try to make a wheel...
   //
   TGeoVolumeAssembly *wheel = new TGeoVolumeAssembly("TPC_ENDCAP");
   //
   for(Int_t i =0;i<18;i++){
     Double_t phi = (20.*i);
     TGeoRotation *r = new TGeoRotation();
     r->RotateZ(phi);
     wheel->AddNode(sect,i+1,r);
    
   }
   // wheels in the drift volume!
   TGeoCombiTrans *combi3 = new TGeoCombiTrans("combi3",0.,0.,256.6,ref);
   v9->AddNode(wheel,1,combi3);
   v9->AddNode(wheel,2,new TGeoTranslation(0.,0.,-256.6));
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
 TGeoVolumeAssembly *swheel = new TGeoVolumeAssembly("TPC_SSWHEEL");
   for(Int_t i =0;i<18;i++){
     Double_t phi = (20.*i);
     TGeoRotation *r = new TGeoRotation();
     r->RotateZ(phi);
     swheel->AddNode(swhs,i+1,r);   
   }
   v1->AddNode(swheel,1,new TGeoTranslation(0.,0.,-284.6));
   v1->AddNode(swheel,2,new TGeoTranslation(0.,0.,284.6));
   //
  //----------------------------------------------------------
  // TPc Support Rods - MAKROLON
  //----------------------------------------------------------
  TGeoMedium *m6=gGeoManager->GetMedium("TPC_Makrolon");
  TGeoMedium *m7=gGeoManager->GetMedium("TPC_Cu");
  // upper and lower rods differ in length!
  Double_t *upar;
  upar=NULL;
  gGeoManager->Volume("TPC_Rod","TUBE",m6->GetId(),upar);
  upar=new Double_t [3];
  upar[0]=1.8;
  upar[1]=2.2;
  
  //
  //HV rods - makrolon + 0.58cm (diameter) Cu
  TGeoTube *hvr = new TGeoTube(0.,2.2,126.64);
  TGeoTube *hvc = new TGeoTube(0.,0.29,126.64);
  //
  TGeoVolume *hvrv = new TGeoVolume("TPC_HV_Rod",hvr,m6);
  TGeoVolume *hvcv = new TGeoVolume("TPC_HV_Cable",hvc,m7);
  hvrv->AddNode(hvcv,1);
  for(Int_t i=0;i<18;i++){
    Double_t angle,x,y;
    Double_t z,r; 
    angle=TMath::DegToRad()*20.*(Double_t)i;
    r=81.5;
    x=r * TMath::Cos(angle);
    y=r * TMath::Sin(angle);
    upar[2]=126.64; //lower
    z= 126.96;
    if(i==15){
      v9->AddNode(hvrv,1,new TGeoTranslation(x,y,z));
      v9->AddNode(hvrv,2,new TGeoTranslation(x,y,-z));
    }
    else{
     gGeoManager->Node("TPC_Rod",i+1,"TPC_Drift",x,y,z,0,kTRUE,upar,3);//shaft
     gGeoManager->Node("TPC_Rod",i+18,"TPC_Drift",x,y,-z,0,kTRUE,upar,3);//muon
    }
    r=254.25;
    x=r * TMath::Cos(angle);
    y=r * TMath::Sin(angle);
    upar[2]=126.54; //upper
    z=127.06;
    gGeoManager->Node("TPC_Rod",i+36,"TPC_Drift",x,y,z,0,kTRUE,upar,3);
    gGeoManager->Node("TPC_Rod",i+54,"TPC_Drift",x,y,-z,0,kTRUE,upar,3);
  }

  TGeoVolume *alice = gGeoManager->GetVolume("ALIC");
  alice->AddNode(v1,1);
  

} // end of function


//_____________________________________________________________________________
void AliTPCv0::CreateMaterials()
{
  //
  // Define materials for the TPC
  //
  AliTPC::CreateMaterials();
}

//_____________________________________________________________________________
void AliTPCv0::DrawDetector()
{
  //
  // Draw a shaded view of the Time Projection Chamber version 0
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("TPC ","SEEN",0);
  gMC->Gsatt("TOIN","SEEN",1);
  gMC->Gsatt("TOIN","COLO",7);
  gMC->Gsatt("TOCV","SEEN",1);
  gMC->Gsatt("TOCV","COLO",4);
  gMC->Gsatt("TSA1","SEEN",0);
  gMC->Gsatt("TSA2","SEEN",0);
  gMC->Gsatt("TSA3","SEEN",0);
  gMC->Gsatt("TSA4","SEEN",0);  
  gMC->Gsatt("TSA5","SEEN",0);
  gMC->Gsatt("TOFC","SEEN",1);
  gMC->Gsatt("TOFC","COLO",4);
  gMC->Gsatt("TSA6","SEEN",0);
  gMC->Gsatt("TSA7","SEEN",0);
  gMC->Gsatt("TSA8","SEEN",0);    
  gMC->Gsatt("TIIN","SEEN",1);
  gMC->Gsatt("TIIN","COLO",7);
  gMC->Gsatt("TII1","SEEN",0);
  gMC->Gsatt("TIFC","SEEN",1);
  gMC->Gsatt("TIFC","COLO",4);
  gMC->Gsatt("TSA9","SEEN",0); 
  gMC->Gsatt("TS10","SEEN",0);
  gMC->Gsatt("TS11","SEEN",0);
  gMC->Gsatt("TS12","SEEN",0);
  gMC->Gsatt("TS13","SEEN",0);
  gMC->Gsatt("TS14","SEEN",0);
  gMC->Gsatt("TICC","SEEN",0);
  gMC->Gsatt("TICM","SEEN",0);
  gMC->Gsatt("TS15","SEEN",0);
  gMC->Gsatt("TS16","SEEN",0);
  gMC->Gsatt("TS17","SEEN",0);
  gMC->Gsatt("TS18","SEEN",0);  
  gMC->Gsatt("TS19","SEEN",0); 
  gMC->Gsatt("TPJ1","SEEN",0);
  gMC->Gsatt("TPJ2","SEEN",0);
  gMC->Gsatt("TICS","SEEN",0);
  gMC->Gsatt("TDGN","SEEN",0); 
  gMC->Gsatt("TIRC","SEEN",0);
  gMC->Gsatt("TIC1","SEEN",1);
  gMC->Gsatt("TIPP","SEEN",0);
  gMC->Gsatt("TIC3","SEEN",0);
  gMC->Gsatt("TRCE","SEEN",0);
  gMC->Gsatt("TPSC","SEEN",0);
  gMC->Gsatt("TPCC","SEEN",0); 
  gMC->Gsatt("TORC","SEEN",0);
  gMC->Gsatt("TOPP","SEEN",0);
  gMC->Gsatt("TOC3","SEEN",0);
  gMC->Gsatt("TOC1","SEEN",1);
  gMC->Gsatt("TSSW","SEEN",1);
  gMC->Gsatt("TSWC","SEEN",1);
  gMC->Gsatt("TSSW","COLO",3);
  gMC->Gsatt("TSWC","COLO",3);
  gMC->Gsatt("TSCE","COLO",6);
  gMC->Gsatt("TSCE","SEEN",1);
  gMC->Gsatt("TWES","SEEN",0);
  gMC->Gsatt("TSWB","SEEN",0);
  gMC->Gsatt("TPEL","SEEN",0);
  gMC->Gsatt("TPMW","SEEN",1);
  gMC->Gsatt("TESR","SEEN",1);
  gMC->Gsatt("TPMW","COLO",12);
  gMC->Gsatt("TIC1","COLO",5);
  gMC->Gsatt("TOC1","COLO",5);
  gMC->Gsatt("TESB","SEEN",0);
  gMC->Gsatt("THVM","SEEN",1);
  gMC->Gsatt("THVM","COLO",11);
  gMC->Gsatt("THVH","SEEN",0);
  gMC->Gsatt("TPSR","SEEN",0); 
  gMC->Gsatt("THVL","SEEN",0);
  gMC->Gsatt("THVC","SEEN",0);  
  gMC->Gsatt("THVE","SEEN",0);
  gMC->Gsatt("THVR","SEEN",0);

  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("TPMW",-300,300,-300,300,254.,270.);
  gMC->SetClipBox("TESR",-300,300,-300,300,254.,270.);
  gMC->SetClipBox("TSSW",-300,300,-300,300,283.,284.);
  gMC->SetClipBox("TSWC",-300,300,-300,300,283.,284.);
  gMC->SetClipBox("*", 0, 300, -300, 300, -290, 290);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 9.5, .025, .025);
  gMC->Gdhead(1111, "Time Projection Chamber");
  gMC->Gdman(18, 4, "MAN");
  gMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliTPCv0::Init()
{
  //
  // Initialise Time Projection Chamber version 0
  //
 
  printf("%s: *** TPC version 0 initialized***\n",ClassName()); 

  //  printf("TPC version 0 initialized\n");
}

//_____________________________________________________________________________
void AliTPCv0::StepManager()
{
  //
  // Procedure called at each step in the TPC
  //
}
