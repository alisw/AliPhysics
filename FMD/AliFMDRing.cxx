/**************************************************************************
 * Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved.      *
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

//__________________________________________________________________
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

#include <math.h>               // fmod
#include <AliLog.h>		// ALILOG_H
#include "AliFMDRing.h"		// ALIFMDRING_H
#include "AliFMD.h"		// ALIFMD_H
#include <TMath.h>		// ROOT_TMath
#include <TH2.h>		// ROOT_TH2
#include <TVirtualMC.h>		// ROOT_TVirtualMC
#include <TVector2.h>		// ROOT_TVector2
#include <TBrowser.h>		// ROOT_TBrowser
#include <TString.h>		// ROOT_TString
#include <TArc.h>		// ROOT_TArc
#include <TObjArray.h>		// ROOT_TObjArray
#include <TXTRU.h>		// ROOT_TXTRU
#include <TNode.h>		// ROOT_TNode
#include <TRotMatrix.h>		// ROOT_TRotMatrix
#include <TList.h>		// ROOT_TList

const Char_t* AliFMDRing::fgkRingFormat         = "F%cRG";
const Char_t* AliFMDRing::fgkVirtualFormat      = "F%cV%c";
const Char_t* AliFMDRing::fgkActiveFormat       = "F%cAC";
const Char_t* AliFMDRing::fgkSectorFormat       = "F%cAP";
const Char_t* AliFMDRing::fgkStripFormat        = "F%cAR";
const Char_t* AliFMDRing::fgkPrintboardFormat   = "F%cP%c";


//____________________________________________________________________
ClassImp(AliFMDRing)

//____________________________________________________________________
AliFMDRing::AliFMDRing(Char_t id, Bool_t detailed) 
  : fId(id), 
    fDetailed(detailed),
    fActiveId(0),
    fPrintboardBottomId(0),
    fPrintboardTopId(0),
    fRingId(0),
    fSectionId(0),
    fStripId(0),
    fVirtualBackId(0),
    fVirtualFrontId(0),
    fBondingWidth(0),
    fWaferRadius(0), 
    fSiThickness(0),
    fLowR(0), 
    fHighR(0), 
    fTheta(0), 
    fNStrips(0), 
    fRingDepth(0),
    fLegRadius(0),
    fLegLength(0),
    fLegOffset(0),
    fModuleSpacing(0),
    fPrintboardThickness(0),
    fShape(0),
    fRotMatricies(0)
{
  // Construct a alifmdring. 
  // 
  //     id		Id of the ring (either 'i' or 'o').
  //     detailed       Whether the strips are made or not.
  // 
}

//____________________________________________________________________
AliFMDRing::AliFMDRing(const AliFMDRing& other) 
  : TObject(other),
    fId(other.fId), 
    fDetailed(other.fDetailed),
    fActiveId(other.fActiveId),
    fPrintboardBottomId(other.fPrintboardBottomId),
    fPrintboardTopId(other.fPrintboardTopId),
    fRingId(other.fRingId),
    fSectionId(other.fSectionId),
    fStripId(other.fStripId),
    fVirtualBackId(other.fVirtualBackId),
    fVirtualFrontId(other.fVirtualFrontId),
    fBondingWidth(other.fBondingWidth),
    fWaferRadius(other.fWaferRadius), 
    fSiThickness(other.fSiThickness),
    fLowR(other.fLowR), 
    fHighR(other.fHighR), 
    fTheta(other.fTheta), 
    fNStrips(other.fNStrips), 
    fRingDepth(other.fRingDepth),
    fLegRadius(other.fLegRadius),
    fLegLength(other.fLegLength),
    fLegOffset(other.fLegOffset),
    fModuleSpacing(other.fModuleSpacing),
    fPrintboardThickness(other.fPrintboardThickness),
    fRotations(other.fRotations),
    fShape(other.fShape), 
    fRotMatricies(other.fRotMatricies)
{
  // Copy constructor of a AliFMDRing. 
}

//____________________________________________________________________
AliFMDRing&
AliFMDRing::operator=(const AliFMDRing& other) 
{
  // Assignment operator 
  // 
  fId			= other.fId; 
  fDetailed		= other.fDetailed;
  fActiveId		= other.fActiveId;
  fPrintboardBottomId	= other.fPrintboardBottomId;
  fPrintboardTopId	= other.fPrintboardTopId;
  fRingId		= other.fRingId;
  fSectionId		= other.fSectionId;
  fStripId		= other.fStripId;
  fVirtualBackId	= other.fVirtualBackId;
  fVirtualFrontId	= other.fVirtualFrontId;
  fBondingWidth		= other.fBondingWidth;
  fWaferRadius		= other.fWaferRadius; 
  fSiThickness		= other.fSiThickness;
  fLowR			= other.fLowR; 
  fHighR		= other.fHighR; 
  fTheta		= other.fTheta; 
  fNStrips		= other.fNStrips; 
  fRingDepth		= other.fRingDepth;
  fLegRadius		= other.fLegRadius;
  fLegLength		= other.fLegLength;
  fLegOffset		= other.fLegOffset;
  fModuleSpacing	= other.fModuleSpacing;
  fPrintboardThickness	= other.fPrintboardThickness;
  fRotations            = other.fRotations;
  if (other.fShape)  {
    if (other.fShape->IsA() == TXTRU::Class()) 
      ((TXTRU*)other.fShape)->Copy(*fShape);
    else 
      fShape = 0;
  }
  if (other.fRotMatricies) {
    Int_t n = other.fRotMatricies->GetEntries();
    if (!fRotMatricies) fRotMatricies = new TObjArray(n);
    else                fRotMatricies->Expand(n);
    TIter next(other.fRotMatricies);
    TObject* o = 0;
    while ((o = next())) fRotMatricies->Add(o);
  }
  return *this;
}

  

//____________________________________________________________________
void 
AliFMDRing::Init() 
{
  // Initialize the ring object.
  // DebugGuard guard("AliFMDRing::Init");
  AliDebug(30, Form("\tInitializing ring %c", fId));
  fPolygon.Clear();
  SetupCoordinates();  
}

//____________________________________________________________________
AliFMDRing::~AliFMDRing() 
{
  // Destructor - deletes shape and rotation matricies 
  AliDebug(30, Form("\tDestructing ring %c", fId));
  if (fShape) delete fShape;
  if (fRotMatricies) delete fRotMatricies;
}


//____________________________________________________________________
void 
AliFMDRing::Browse(TBrowser* /* b */)
{
  // DebugGuard guard("AliFMDRing::Browse");
  AliDebug(30, Form("\tBrowsing ring %c", fId));
}

  
//____________________________________________________________________
void 
AliFMDRing::SetupCoordinates() 
{
  // Calculates the parameters of the polygon shape. 
  // 
  //

  // Get out immediately if we have already done all this 
  if (fPolygon.GetNVerticies() > 1) return;
  AliDebug(10, Form("\tSetting up the coordinates for ring %c", fId));

  double tanTheta  = TMath::Tan(fTheta * TMath::Pi() / 180.);
  double tanTheta2 = TMath::Power(tanTheta,2);
  double r2         = TMath::Power(fWaferRadius,2);
  double yA        = tanTheta * fLowR;
  double lr2        = TMath::Power(fLowR, 2);
  double hr2        = TMath::Power(fHighR,2);
  double xD        = fLowR + TMath::Sqrt(r2 - tanTheta2 * lr2);
  double xD2       = TMath::Power(xD,2);
  //double xD_2      = fLowR - TMath::Sqrt(r2 - tanTheta2 * lr2);
  double yB        = TMath::Sqrt(r2 - hr2 + 2 * fHighR * xD - xD2);
  double xC        = ((xD + TMath::Sqrt(-tanTheta2 * xD2 + r2 
					+ r2 * tanTheta2)) 
		       / (1 + tanTheta2));
  double yC        = tanTheta * xC;

  fPolygon.AddVertex(fLowR,  -yA);
  fPolygon.AddVertex(xC,     -yC);
  fPolygon.AddVertex(fHighR, -yB);
  fPolygon.AddVertex(fHighR,  yB);
  fPolygon.AddVertex(xC,      yC);
  fPolygon.AddVertex(fLowR,   yA);
}

//____________________________________________________________________
bool
AliFMDRing::IsWithin(size_t moduleNo, double x, double y) const
{
  // Checks if a point (x,y) is inside the module with number moduleNo 
  //
  // DebugGuard guard("AliFMDRing::IsWithin");
  AliDebug(20, Form("\tChecking wether the hit at (%lf,%lf) in module %d "
		    "is within this ring (%c)", x, y, moduleNo, fId));
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
  AliDebug(20, Form("\tDrawing ring %c", fId));
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
  //    FRGx                        // Ring volume
  //      FVFx                      // Container of hybrid + legs
  //        FACx                    // Active volume (si sensor approx)
  //          FSEx                  // Section division
  //            FSTx                // Strip division 
  //        FPTx                    // Print board (bottom)
  //        FPBx                    // Print board (top)  
  //        FLL                     // Support leg (long version)
  //      FVBx                      // Container of hybrid + legs
  //        FACx                    // Active volume (si sensor approx)
  //          FSEx                  // Section division
  //            FSTx                // Strip division 
  //        FPTx                    // Print board (bottom)
  //        FPBx                    // Print board (top)  
  //        FSL                     // Support leg (long version)
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
  AliDebug(10, Form("\tSetting up the geometry for ring %c", fId));

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
  pars[1]             = rmax;
  pars[2]             = fRingDepth / 2;
  name                = Form(fgkRingFormat, fId);
  fRingId             = gMC->Gsvolu(name.Data(), "TUBE", vacuumId, pars, 3);
  
  // Virtual volume for modules with long legs 
  pars[1]             = rmax;
  pars[3]             = -fTheta;
  pars[4]             =  fTheta;
  name                = Form(fgkVirtualFormat, fId, 'F');
  fVirtualFrontId     = gMC->Gsvolu(name.Data(), "TUBS", vacuumId, pars, 5);

  // Virtual volume for modules with long legs 
  pars[2]             =  (fRingDepth - fModuleSpacing) / 2;
  name                =  Form(fgkVirtualFormat, fId, 'B');
  fVirtualBackId      =  gMC->Gsvolu(name.Data(), "TUBS", vacuumId, pars, 5);
  
  // Virtual mother volume for silicon
  pars[2]             =  fSiThickness/2;
  name2               =  name;
  name                =  Form(fgkActiveFormat, fId);
  fActiveId           =  gMC->Gsvolu(name.Data(), "TUBS", vacuumId , pars, 5);

  if (fDetailed) {
    // Virtual sector volumes 
    name2               = name;
    name                = Form(fgkSectorFormat, fId);
    gMC->Gsdvn2(name.Data(), name2.Data(), 2, 2, -fTheta, vacuumId);
    fSectionId          = gMC->VolId(name.Data());
    
    // Active strip volumes 
    name2               = name;
    name                = Form(fgkStripFormat, fId);
    gMC->Gsdvt2(name.Data(), name2.Data(), dStrip, 1,stripOff, siId, fNStrips);
    fStripId            = gMC->VolId(name.Data());
  }
  
  // Print-board on back of module 
  pars[4]             = TMath::Tan(TMath::Pi() * fTheta / 180) * fBondingWidth;
  // Top of the print board
  pars[0]             = cCorner.Y() - pars[4];
  pars[1]             = bCorner.Y() - pars[4];
  pars[2]             = fPrintboardThickness / 2; // PCB half thickness
  pars[3]             = (bCorner.X() - cCorner.X()) / 2;
  name                = Form(fgkPrintboardFormat, fId, 'T');
  fPrintboardTopId    = gMC->Gsvolu(name.Data(), "TRD1", pcbId, pars, 4);

  // Bottom of the print board
  pars[0]             = aCorner.Y() - pars[4];
  pars[1]             = cCorner.Y() - pars[4];
  pars[3]             = (cCorner.X() - aCorner.X()) / 2;
  name                = Form(fgkPrintboardFormat, fId, 'B');
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
    TString  name2    = Form(fgkRingFormat, fId);

    Int_t     id      = i;
    // Double_t  theta   = (i + .5) * dTheta;
    Bool_t    isFront = (i % 2 == 1);
    Double_t  dz      = 0;
    Double_t  w       = fRingDepth - (isFront ? 0 : fModuleSpacing);

    // Place virtual module volume 
    name = Form(fgkVirtualFormat, fId, (isFront ? 'F' : 'B'));
    dz   = (w - fRingDepth) / 2;
    gMC->Gspos(name.Data(), id, name2.Data(), 0., 0., dz,fRotations[i], 
	       "ONLY");

    // We only need to place the children once, they are copied when
    // we place the other virtual volumes. 
    if (i > 1) continue;
    name2 = name;

    // Place active silicon wafer - this is put so that the front of
    // the silicon is on the edge of the virtual volume. 
    name  = Form(fgkActiveFormat, fId);
    dz    = (w - fSiThickness) / 2;
    gMC->Gspos(name.Data(), id, name2.Data(),0.,0.,dz,idRotId, "ONLY");

    // Place print board.  This is put immediately behind the silicon
    name = Form(fgkPrintboardFormat, fId, 'T');
    dz   =  w / 2 - fSiThickness - fPrintboardThickness / 2;
    gMC->Gspos(name.Data(), id, name2.Data(), 
	       fLowR + pbBotL + pbTopL / 2, 0, dz, pbRotId, "ONLY");
    name = Form(fgkPrintboardFormat, fId, 'B');
    gMC->Gspos(name.Data(), id, name2.Data(), 
	       fLowR + pbBotL / 2, 0, dz, pbRotId, "ONLY");

    // Support legs 
    // This is put immediately behind the pringboard. 
    dz     = (w / 2 - fSiThickness - fPrintboardThickness 
	     - (fLegLength + (isFront ? fModuleSpacing : 0)) /2);
    name  = (isFront ? AliFMD::fgkLongLegName : AliFMD::fgkShortLegName);
    gMC->Gspos(name.Data(), id*10 + 1, name2.Data(), 
	       aCorner.X() + fLegOffset + fLegRadius, 0., dz, idRotId, "ONLY");
    Double_t y = cCorner.Y() - yoffset - fLegOffset - fLegRadius;
    gMC->Gspos(name.Data(),id*10+2,name2.Data(),cCorner.X(),y,dz,
	       idRotId,"ONLY");
    gMC->Gspos(name.Data(),id*10+3,name2.Data(),cCorner.X(),-y,dz,
	       idRotId,"ONLY");
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
  TString  name;
  Double_t offsetZ   = (fSiThickness 
			+ fPrintboardThickness 
			+ fLegLength + fModuleSpacing) / 2;
  name = Form(fgkRingFormat, fId);
  AliDebug(10, Form("\tPlacing ring %s in %s at z=%lf-%lf=%lf (base ID: %d)", 
		    name.Data(), mother, z, offsetZ, z-offsetZ, baseId));
  gMC->Gspos(name.Data(), baseId, mother, 0., 0., z - offsetZ, idRotId, 
	     "ONLY");
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
  SetupCoordinates();

  AliDebug(10, Form("\tCreating simple geometry for "
		    "ring %c at z=%lf cm in %s", 
		    fId, z, mother->GetName()));
  // If the shape hasn't been defined yet, we define it here. 
  if (!fShape) {

    TString name(Form(fgkActiveFormat, fId));
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
    TString name(Form("FAC%c_%d_%d", fId, n, i));
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
  AliDebug(10, Form("\tSetting drawing attributes for Ring %c", fId));
  TString name;
  name = Form(fgkRingFormat,fId);
  gMC->Gsatt(name.Data(), "SEEN", 0);

  name = Form(fgkVirtualFormat, fId, 'F');
  gMC->Gsatt(name.Data(), "SEEN", 0);

  name = Form(fgkVirtualFormat, fId, 'B');
  gMC->Gsatt(name.Data(), "SEEN", 0);

  name = Form(fgkActiveFormat,fId);
  gMC->Gsatt(name.Data(), "SEEN", 1);

  name = Form(fgkSectorFormat,fId);
  gMC->Gsatt(name.Data(), "SEEN", 0);

  name = Form(fgkStripFormat,fId);
  gMC->Gsatt(name.Data(), "SEEN", 0);

  name = Form(fgkPrintboardFormat, fId, 'T');
  gMC->Gsatt(name.Data(), "SEEN", 1);

  name = Form(fgkPrintboardFormat, fId, 'B');
  gMC->Gsatt(name.Data(), "SEEN", 1);
}

//
// EOF
//
