/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
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

//////////////////////////////////////////////////////////////////////////////
//
// Utility class to help implement collection of FMD modules into
// rings.  This is used by AliFMDSubDetector and AliFMD.  
//
// The AliFMD object owns the AliFMDRing objects, and the
// AliFMDSubDetector objects reference these.  That is, the AliFMDRing
// objects are share amoung the AliFMDSubDetector objects.
//
// Latest changes by Christian Holm Christensen
//
//////////////////////////////////////////////////////////////////////////////
#ifndef ALIFMDRING_H
# include "AliFMDRing.h"
#endif
#ifndef ALILOG_H
# include "AliLog.h"
#endif
#ifndef ROOT_TMath
# include <TMath.h>
#endif
#ifndef ROOT_TH2
# include <TH2.h>
#endif
#ifndef ROOT_TVirtualMC
# include <TVirtualMC.h>
#endif
#ifndef ROOT_TVector2
# include <TVector2.h>
#endif
#ifndef ROOT_TBrowser
# include <TBrowser.h>
#endif
#ifndef ROOT_TString
# include <TString.h>
#endif
#ifndef ROOT_TArc
# include <TArc.h>
#endif
#ifndef ROOT_TObjArray
# include <TObjArray.h>
#endif
#ifndef ROOT_TXTRU
# include <TXTRU.h>
#endif
#ifndef ROOT_TNode
# include <TNode.h>
#endif
#ifndef ROOT_TRotMatrix
# include <TRotMatrix.h>
#endif
#ifndef ROOT_TList
# include <TList.h>
#endif
#ifndef __IOSTREAM__
# include <iostream>
#endif

ClassImp(AliFMDRing);

//____________________________________________________________________
// Construct a alifmdring. 
// 
//     id		Id of the ring (either 'i' or 'o').
//     lowr		Lower radius of ring (in centimeters).
//     highr		Upper radius of ring (in centimeters).
//     r		Radius of the silicon wafers (in centimeters). 
//     theta		Opening angle of the silicon wafers.
//     strips		Number of strips. 
AliFMDRing::AliFMDRing(Char_t id, Bool_t detailed) 
  : fId(id), 
    fDetailed(detailed),
    fWaferRadius(0), 
    fSiThickness(0),
    fLowR(0), 
    fHighR(0), 
    fTheta(0), 
    fNStrips(0), 
    fShape(0),
    fRotMatricies(0)
{}

//____________________________________________________________________
void 
AliFMDRing::Init() 
{
  // Initialize the ring object.
  // DebugGuard guard("AliFMDRing::Init");
  AliDebug(10, "AliFMDRing::Init");
  fPolygon.Clear();
  SetupCoordinates();  
}

//____________________________________________________________________
AliFMDRing::~AliFMDRing() 
{
  if (fShape) delete fShape;
  if (fRotMatricies) delete fRotMatricies;
}


//____________________________________________________________________
void 
AliFMDRing::Browse(TBrowser* /* b */)
{
  // DebugGuard guard("AliFMDRing::Browse");
  AliDebug(10, "AliFMDRing::Browse");
}

  
//____________________________________________________________________
void 
AliFMDRing::SetupCoordinates() 
{
  // Calculates the parameters of the polygon shape. 
  // 
  // DebugGuard guard("AliFMDRing::SetupCoordinates");
  AliDebug(10, "AliFMDRing::SetupCoordinates");
  // Get out immediately if we have already done all this 
  if (fPolygon.GetNVerticies() > 1) return;

  double tan_theta  = TMath::Tan(fTheta * TMath::Pi() / 180.);
  double tan_theta2 = TMath::Power(tan_theta,2);
  double r2         = TMath::Power(fWaferRadius,2);
  double y_A        = tan_theta * fLowR;
  double lr2        = TMath::Power(fLowR, 2);
  double hr2        = TMath::Power(fHighR,2);
  double x_D        = fLowR + TMath::Sqrt(r2 - tan_theta2 * lr2);
  double x_D2       = TMath::Power(x_D,2);
  //double x_D_2      = fLowR - TMath::Sqrt(r2 - tan_theta2 * lr2);
  double y_B        = sqrt(r2 - hr2 + 2 * fHighR * x_D - x_D2);
  double x_C        = ((x_D + TMath::Sqrt(-tan_theta2 * x_D2 + r2 
					  + r2 * tan_theta2)) 
		       / (1 + tan_theta2));
  double y_C        = tan_theta * x_C;

  fPolygon.AddVertex(fLowR,  -y_A);
  fPolygon.AddVertex(x_C,    -y_C);
  fPolygon.AddVertex(fHighR, -y_B);
  fPolygon.AddVertex(fHighR,  y_B);
  fPolygon.AddVertex(x_C,     y_C);
  fPolygon.AddVertex(fLowR,   y_A);
}

//____________________________________________________________________
bool
AliFMDRing::IsWithin(size_t moduleNo, double x, double y) const
{
  // Checks if a point (x,y) is inside the module with number moduleNo 
  //
  // DebugGuard guard("AliFMDRing::IsWithin");
  AliDebug(10, "AliFMDRing::IsWithin");
  bool   ret            = false;
  double r2             = x * x + y * y;
  if (r2 < fHighR * fHighR && r2 > fLowR * fLowR) {
    // double point_angle    = TMath::ATan2(y, x);
    // int    n_modules      = 360 / Int_t(fTheta * 2);
    double m_angle        = (.5 + moduleNo) * 2 * fTheta;
    double m_radians      = TMath::Pi() * m_angle / 180.;
    
    // Rotate the point.
    double xr = x * TMath::Cos(-m_radians) - y * TMath::Sin(-m_radians);
    double yr = x * TMath::Sin(-m_radians) + y * TMath::Cos(-m_radians);
  
    ret = fPolygon.Contains(xr,yr);
  }
  return ret;
}




//____________________________________________________________________
void
AliFMDRing::Draw(Option_t* option) const
{
  // Draw a the shape of the ring into a 2D histogram.  Useful for
  // superimposing the actual shape of the ring onto a scatter plot of
  // hits in the detector. 
  // 
  // DebugGuard guard("AliFMDRing::Draw");
  AliDebug(10, "AliFMDRing::Draw");
  // The unrotated coordinates of the polygon verticies
  if (fPolygon.GetNVerticies() < 1) return;
  
  TVector2 v[6];
  for (size_t i = 0; i < fPolygon.GetNVerticies(); i++) 
    v[i] = fPolygon.GetVertex(i);
  
  Int_t    nModules  = 360 / Int_t(fTheta * 2);
  Double_t dTheta    = fTheta * 2;
  
  TString opt(option);
  if (opt.Contains("B", TString::kIgnoreCase)) {
    opt.Remove(opt.Index("B", 1, TString::kIgnoreCase),1);
    TH1* null = new TH2F("null", "Null",
			 100, -fHighR * 1.1, fHighR * 1.1, 
			 100, -fHighR * 1.1, fHighR * 1.1);
    null->SetStats(0);
    null->Draw(opt.Data());
  }
   
  for (int i = 0; i < nModules; i++) {
    Double_t theta = (i + .5) * dTheta;
    AliFMDPolygon p;
    for (int j = 0; j < 6; j++) {
      TVector2 vr(v[j].Rotate(TMath::Pi() * theta / 180.));
      if (!p.AddVertex(vr.X(),vr.Y())) {
	// std::cerr << "Draw of polygon " << i << " failed" << std::endl;
	break;
      }
    }
    p.Draw(opt.Data(), Form("MOD%c_%d", fId, i));
  }
  if (opt.Contains("0", TString::kIgnoreCase)) {
    TArc* arcH = new TArc(0,0, fHighR);
    arcH->SetLineStyle(2);
    arcH->SetLineColor(4);
    arcH->Draw();

    TArc* arcL = new TArc(0,0, fLowR);
    arcL->SetLineStyle(2);
    arcL->SetLineColor(4);
    arcL->Draw();
  }
}

//____________________________________________________________________
void 
AliFMDRing::SetupGeometry(Int_t vacuumId, Int_t siId, Int_t pcbId, 
			  Int_t pbRotId, Int_t idRotId)
{
  // Setup the geometry of the ring.  It defines the volumes 
  // RNGI or RNGO which can later be positioned in a sub-detector
  // volume. 
  // 
  // The hieracy of the RNGx volume is 
  // 
  //    RNGx                        // Ring volume
  //      VFx                       // Container of hybrid + legs
  //        ACTx                    // Active volume (si sensor approx)
  //          SECx                  // Section division
  //            STRx                // Strip division 
  //        PBTx                    // Print board (bottom)
  //        PTTx                    // Print board (top)  
  //        LLEG                    // Support leg (long version)
  //      VBx                       // Container of hybrid + legs
  //        ACTx                    // Active volume (si sensor approx)
  //          SECx                  // Section division
  //            STRx                // Strip division 
  //        PBTx                    // Print board (bottom)
  //        PTTx                    // Print board (top)  
  //        SLEG                    // Support leg (long version)
  //        
  // Parameters: 
  //
  //   vacuumId        Medium of inactive virtual volumes 
  //   siId            Medium of Silicon sensor (active)
  //   pcbId           Medium of print boards 
  //   pbRotId         Print board rotation matrix 
  //   idRotId         Identity rotation matrix 
  //
  // DebugGuard guard("AliFMDRing::SetupGeometry");
  AliDebug(10, "AliFMDRing::SetupGeometry");

  const TVector2& bCorner   = fPolygon.GetVertex(3); // Third  corner
  const TVector2& aCorner   = fPolygon.GetVertex(5); // First  corner
  const TVector2& cCorner   = fPolygon.GetVertex(4); // Second corner
  TString name;
  TString name2;
  Double_t dStrip     = (bCorner.Mod() - aCorner.Mod()) / fNStrips;
  Double_t stripOff   = aCorner.Mod();
  Double_t rmin       = fLowR;
  Double_t rmax       = bCorner.Mod();
  Double_t pars[10];
  fRingDepth          = (fSiThickness 
			 + fPrintboardThickness 
			 + fLegLength 
			 + fModuleSpacing);

  // Ring virtual volume 
  pars[0]             = rmin;
  pars[1]             = fHighR;
  pars[2]             = fRingDepth / 2;
  name                = Form("RNG%c", fId);
  fRingId             = gMC->Gsvolu(name.Data(), "TUBE", vacuumId, pars, 3);
  
  // Virtual volume for modules with long legs 
  pars[1]             = rmax;
  pars[3]             = -fTheta;
  pars[4]             =  fTheta;
  name                = Form("VF%c", fId);
  fVirtualFrontId     = gMC->Gsvolu(name.Data(), "TUBS", vacuumId, pars, 5);

  // Virtual volume for modules with long legs 
  pars[2]             =  (fRingDepth - fModuleSpacing) / 2;
  name                =  Form("VB%c", fId);
  fVirtualBackId      =  gMC->Gsvolu(name.Data(), "TUBS", vacuumId, pars, 5);
  
  // Virtual mother volume for silicon
  pars[2]             =  fSiThickness/2;
  name2               =  name;
  name                =  Form("ACT%c",fId);
  fActiveId           =  gMC->Gsvolu(name.Data(), "TUBS", vacuumId , pars, 5);

  if (fDetailed) {
    // Virtual sector volumes 
    name2               = name;
    name                = Form("SEC%c",fId);
    gMC->Gsdvn2(name.Data(), name2.Data(), 2, 2, -fTheta, vacuumId);
    fSectionId          = gMC->VolId(name.Data());
    
    // Active strip volumes 
    name2               = name;
    name                = Form("STR%c", fId);
    gMC->Gsdvt2(name.Data(), name2.Data(), dStrip, 1, 
		stripOff, siId, fNStrips);
    fStripId            = gMC->VolId(name.Data());
  }
  
  // Print-board on back of module 
  pars[4]             = TMath::Tan(TMath::Pi() * fTheta / 180) * fBondingWidth;
  // Top of the print board
  pars[0]             = cCorner.Y() - pars[4];
  pars[1]             = bCorner.Y() - pars[4];
  pars[2]             = fPrintboardThickness / 2; // PCB half thickness
  pars[3]             = (bCorner.X() - cCorner.X()) / 2;
  name                = Form("PBT%c", fId);
  fPrintboardTopId    = gMC->Gsvolu(name.Data(), "TRD1", pcbId, pars, 4);

  // Bottom of the print board
  pars[0]             = aCorner.Y() - pars[4];
  pars[1]             = cCorner.Y() - pars[4];
  pars[3]             = (cCorner.X() - aCorner.X()) / 2;
  name                = Form("PBB%c", fId);
  fPrintboardBottomId = gMC->Gsvolu(name.Data(), "TRD1", pcbId, pars, 4);

  // Define rotation matricies
  Int_t    nModules  = 360 / Int_t(fTheta * 2);
  Double_t dTheta    = fTheta * 2;
  fRotations.Set(nModules);
  for (int i = 0; i < nModules; i++) {
    Double_t theta  = (i + .5) * dTheta;
    Int_t    idrot  = 0;
    // Rotation matrix for virtual module volumes
    gMC->Matrix(idrot, 90, theta, 90, fmod(90 + theta, 360), 0, 0);
    fRotations[i] = idrot;
  }


  // Int_t    nModules  = 360 / Int_t(fTheta * 2);
  // Double_t dTheta    = fTheta * 2;
  Double_t pbTopL    = (bCorner.X() - cCorner.X());
  Double_t pbBotL    = (cCorner.X() - aCorner.X());
  Double_t yoffset   = ((TMath::Tan(TMath::Pi() * fTheta / 180) 
			 * fBondingWidth)); 
  
  for (int i = 0; i < nModules; i++) {
    TString  name2    = Form("RNG%c", fId);

    Int_t     id      = i;
    // Double_t  theta   = (i + .5) * dTheta;
    Bool_t    isFront = (i % 2 == 1);
    Double_t  dz      = 0;
    Double_t  w       = fRingDepth - (isFront ? 0 : fModuleSpacing);

    // Place virtual module volume 
    name = Form("V%c%c", (isFront ? 'F' : 'B'), fId);
    dz   = (w - fRingDepth) / 2;
    gMC->Gspos(name.Data(), id, name2.Data(), 0., 0., dz,fRotations[i]);

    // We only need to place the children once, they are copied when
    // we place the other virtual volumes. 
    if (i > 1) continue;
    name2 = name;

    // Place active silicon wafer - this is put so that the front of
    // the silicon is on the edge of the virtual volume. 
    name  = Form("ACT%c", fId);
    dz    = (w - fSiThickness) / 2;
    gMC->Gspos(name.Data(), id, name2.Data(),0.,0.,dz,idRotId);

    // Place print board.  This is put immediately behind the silicon
    name = Form("PBT%c", fId);
    dz   =  w / 2 - fSiThickness - fPrintboardThickness / 2;
    gMC->Gspos(name.Data(), id, name2.Data(), 
	       fLowR + pbBotL + pbTopL / 2, 0, dz, pbRotId, "ONLY");
    name = Form("PBB%c", fId);
    gMC->Gspos(name.Data(), id, name2.Data(), 
	       fLowR + pbBotL / 2, 0, dz, pbRotId, "ONLY");

    // Support legs 
    // This is put immediately behind the pringboard. 
    dz     = (w / 2 - fSiThickness - fPrintboardThickness 
	     - (fLegLength + (isFront ? fModuleSpacing : 0)) /2);
    name  = (isFront ? "LLEG" : "SLEG");
    gMC->Gspos(name.Data(), id*10 + 1, name2.Data(), 
	       aCorner.X() + fLegOffset + fLegRadius, 0., dz, idRotId, "");
    Double_t y = cCorner.Y() - yoffset - fLegOffset - fLegRadius;
    gMC->Gspos(name.Data(),id*10+2,name2.Data(),cCorner.X(), y,dz,idRotId,"");
    gMC->Gspos(name.Data(),id*10+3,name2.Data(),cCorner.X(), -y,dz ,idRotId,"");
  }
}
//____________________________________________________________________
void 
AliFMDRing::Geometry(const char* mother, Int_t baseId, Double_t z, 
		     Int_t /* pbRotId */, Int_t idRotId)
{
  // Positions a RNGx volume inside a mother. 
  // 
  // Parameters
  //  
  //    mother    Mother volume to position the RNGx volume in 
  //    baseId    Base copy number 
  //    z         Z coordinate where the front of the active silicon
  //              should be in the mother volume, so we need to
  //              subtract half the ring width.  
  //    idRotId   Identity rotation matrix 
  // 
  // DebugGuard guard("AliFMDRing::Geometry");
  AliDebug(10, "AliFMDRing::Geometry");
  TString  name;
  Double_t offsetZ   = (fSiThickness 
			+ fPrintboardThickness 
			+ fLegLength + fModuleSpacing) / 2;
  name = Form("RNG%c", fId);
  gMC->Gspos(name.Data(), baseId, mother, 0., 0., z - offsetZ, idRotId, "");
}

//____________________________________________________________________
void 
AliFMDRing::SimpleGeometry(TList* nodes, 
			   TNode* mother, 
			   Int_t colour, 
			   Double_t z, 
			   Int_t n) 
{
  // Make a simple geometry of the ring for event display. 
  // 
  // The simple geometry is made from ROOT TNode and TShape objects. 
  // Note, that we cache the TShape and TRotMatrix objects used for
  // this. 
  // 
  // Parameters
  // 
  //    nodes     List of nodes to register all create nodes in 
  //    mother    Mother node to put the ring in. 
  //    colour    Colour of the nodes 
  //    z         Z position of the node in the mother volume 
  //    n         Detector number
  //
  // DebugGuard guard("AliFMDRing::SimpleGeometry");
  AliDebug(10, "AliFMDRing::SimpleGeometry");
  SetupCoordinates();

  // If the shape hasn't been defined yet, we define it here. 
  if (!fShape) {
    TString name(Form("ACT%c", fId));
    TString title(Form("Shape of modules in %c Rings", fId));
    Int_t n = fPolygon.GetNVerticies();
    TXTRU* shape = new TXTRU(name.Data(), title.Data(), "void", n, 2);
    for (Int_t i = 0; i < n; i++) {
      const TVector2& v = fPolygon.GetVertex(i);
      shape->DefineVertex(i, v.X(), v.Y());
    }
    shape->DefineSection(0, - fSiThickness / 2, 1, 0, 0);
    shape->DefineSection(1, + fSiThickness / 2, 1, 0, 0);
    fShape = shape;
    fShape->SetLineColor(colour);
  }
  
  Int_t    nModules  = 360 / Int_t(fTheta * 2);
  Double_t dTheta    = fTheta * 2;

  // If the roation matricies hasn't been defined yet, we do so here
  if (!fRotMatricies) {
    fRotMatricies = new TObjArray(nModules);
    for (int i = 0; i < nModules; i++) {
      Double_t theta  = (i + .5) * dTheta;
      TString name(Form("FMD_ring_%c_rot", fId));
      TString title(Form("FMD Ring %c Rotation", fId));
      TRotMatrix* rot = 
	new TRotMatrix(name.Data(), title.Data(), 
		       90, theta, 90, fmod(90 + theta, 360), 0, 0);
      fRotMatricies->AddAt(rot, i);
    }
  }

  Double_t offsetZ   = (fSiThickness 
			+ fPrintboardThickness
			+ fLegLength + fModuleSpacing) / 2;

  // Make all the nodes
  for (int i = 0; i < nModules; i++) {
    Bool_t    isFront = (i % 2 == 1);
    mother->cd();
    TRotMatrix* rot = static_cast<TRotMatrix*>(fRotMatricies->At(i));
    TString name(Form("ACT%c_%d_%d", fId, n, i));
    TString title(Form("Active FMD%d volume in %c Ring", n, fId));
    TNode* node = new TNode(name.Data(), title.Data(), fShape, 
			    0, 0, 
			    z - offsetZ + (isFront ? fModuleSpacing : 0), 
			    rot);
    node->SetLineColor(colour);
    nodes->Add(node);
  }
}

  

//____________________________________________________________________
void 
AliFMDRing::Gsatt() 
{
  // Set drawing attributes for the RING 
  // 
  // DebugGuard guard("AliFMDRing::Gsatt");
  AliDebug(10, "AliFMDRing::Gsatt");
  TString name;
  name = Form("RNG%c",fId);
  gMC->Gsatt(name.Data(), "SEEN", 0);

  name = Form("VF%c",fId);
  gMC->Gsatt(name.Data(), "SEEN", 0);

  name = Form("VB%c",fId);
  gMC->Gsatt(name.Data(), "SEEN", 0);

  name = Form("ACT%c",fId);
  gMC->Gsatt(name.Data(), "SEEN", 1);

  name = Form("PBT%c",fId);
  gMC->Gsatt(name.Data(), "SEEN", 1);

  name = Form("PBB%c",fId);
  gMC->Gsatt(name.Data(), "SEEN", 1);
}

//
// EOF
//
