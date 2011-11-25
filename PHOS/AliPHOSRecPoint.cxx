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
//_________________________________________________________________________
//  Base Class for PHOS Reconstructed Points  
//  Why should I put meaningless comments
//  just to satisfy
//  the code checker                
//-- Author: Gines Martinez (SUBATECH)

// --- ROOT system ---
#include "TPad.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TGeoMatrix.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliPHOSLoader.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"
#include "AliGeomManager.h"

ClassImp(AliPHOSRecPoint)


//____________________________________________________________________________
AliPHOSRecPoint::AliPHOSRecPoint()
  : AliCluster(),fPHOSMod(0),
    fMulTrack(0),fMaxDigit(100),fMulDigit(0),fMaxTrack(200),
    fDigitsList(0),fTracksList(0),fAmp(0),
    fIndexInList(-1), // to be set when the point is already stored
    fLocPos(0,0,0)
{
  // ctor

}

//____________________________________________________________________________
AliPHOSRecPoint::AliPHOSRecPoint(const char * ) 
  : AliCluster(),fPHOSMod(0),
    fMulTrack(0),fMaxDigit(100),fMulDigit(0),fMaxTrack(200),
    fDigitsList(new Int_t[fMaxDigit]),fTracksList(new Int_t[fMaxTrack]),fAmp(0),
    fIndexInList(-1), // to be set when the point is already stored
    fLocPos(0,0,0)

{
  // ctor
  
}
//_______________________________________________________________________
AliPHOSRecPoint::~AliPHOSRecPoint()
{
  // dtor
  
  delete [] fDigitsList ; 
  delete [] fTracksList ;  
  
}
//____________________________________________________________________________
AliPHOSRecPoint::AliPHOSRecPoint(const AliPHOSRecPoint &rp) : 
  AliCluster(rp),
  fPHOSMod(rp.fPHOSMod),fMulTrack(rp.fMulTrack),fMaxDigit(rp.fMaxDigit),
  fMulDigit(rp.fMulDigit),fMaxTrack(rp.fMaxTrack),fDigitsList(0x0),
  fTracksList(0x0),fAmp(rp.fAmp),fIndexInList(rp.fIndexInList), 
  fLocPos(rp.fLocPos)
{
  //copy ctor

  if (rp.fMulDigit>0) fDigitsList = new Int_t[rp.fMulDigit];
  for(Int_t i=0; i<fMulDigit; i++)
    fDigitsList[i] = rp.fDigitsList[i];

  if (rp.fMulTrack>0) fTracksList = new Int_t[rp.fMulTrack];
  for(Int_t i=0; i<fMulTrack; i++)
    fTracksList[i] = rp.fTracksList[i];
  
}
//____________________________________________________________________________
AliPHOSRecPoint& AliPHOSRecPoint::operator= (const AliPHOSRecPoint &rp)
{
  if(&rp == this) return *this;

  fPHOSMod = rp.fPHOSMod;
  fMulTrack = rp.fMulTrack;
  fMaxDigit = rp.fMaxDigit;
  fMulDigit = rp.fMulDigit;
  fMaxTrack = rp.fMaxTrack;
  fAmp = rp.fAmp;
  fIndexInList = rp.fIndexInList; 
  fLocPos = rp.fLocPos;

  if (rp.fMulDigit>0) fDigitsList = new Int_t[rp.fMulDigit];
  for(Int_t i=0; i<fMaxDigit; i++)
    fDigitsList[i] = rp.fDigitsList[i];

  if (rp.fMulTrack>0) fTracksList = new Int_t[rp.fMulTrack];
  for(Int_t i=0; i<fMaxTrack; i++)
    fTracksList[i] = rp.fTracksList[i];

  return *this;
}
//____________________________________________________________________________
Int_t AliPHOSRecPoint::DistancetoPrimitive(Int_t px, Int_t py)
{
  // Compute distance from point px,py to  a AliPHOSRecPoint considered as a Tmarker
  // Compute the closest distance of approach from point px,py to this marker.
  // The distance is computed in pixels units.

  TVector3 pos(0.,0.,0.) ;
  GetLocalPosition( pos) ;
  Float_t x =  pos.X() ;
  Float_t y =  pos.Z() ;
  const Int_t kMaxDiff = 10;
  Int_t pxm  = gPad->XtoAbsPixel(x);
  Int_t pym  = gPad->YtoAbsPixel(y);
  Int_t dist = (px-pxm)*(px-pxm) + (py-pym)*(py-pym);
  
  if (dist > kMaxDiff) return 9999;
  return dist;
}

//___________________________________________________________________________
 void AliPHOSRecPoint::Draw(Option_t *option)
 {
   // Draw this AliPHOSRecPoint with its current attributes
   
   AppendPad(option);
 }

//______________________________________________________________________________
void AliPHOSRecPoint::ExecuteEvent(Int_t event, Int_t, Int_t)
{
  // Execute action corresponding to one event
  // This member function is called when a AliPHOSRecPoint is clicked with the locator
  //
  // If Left button is clicked on AliPHOSRecPoint, the digits are switched on    
  // and switched off when the mouse button is released.

  //  static Int_t pxold, pyold;

  static TGraph *  digitgraph = 0 ;
  static TPaveText* clustertext = 0 ;
  
  if (!gPad->IsEditable()) return;
  
  switch (event) {
    
    
  case kButton1Down:{
    AliPHOSDigit * digit ;
  
    AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance() ;

    Int_t iDigit;
    Int_t relid[4] ;
  
    const Int_t kMulDigit=AliPHOSRecPoint::GetDigitsMultiplicity() ;
    Float_t * xi = new Float_t [kMulDigit] ; 
    Float_t * zi = new Float_t [kMulDigit] ;
    
    for(iDigit = 0; iDigit < kMulDigit; iDigit++) {
      Fatal("AliPHOSRecPoint::ExecuteEvent", "-> Something wrong with the code"); 
      digit = 0 ; //dynamic_cast<AliPHOSDigit *>((fDigitsList)[iDigit]);
      phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
      phosgeom->RelPosInModule(relid, xi[iDigit], zi[iDigit]) ;
    }
    
    if (!digitgraph) {
      digitgraph = new TGraph(fMulDigit,xi,zi);
      digitgraph-> SetMarkerStyle(5) ; 
      digitgraph-> SetMarkerSize(1.) ;
      digitgraph-> SetMarkerColor(1) ;
      digitgraph-> Draw("P") ;
    }
    if (!clustertext) {
      
      TVector3 pos(0.,0.,0.) ;
      GetLocalPosition(pos) ;
      clustertext = new TPaveText(pos.X()-10,pos.Z()+10,pos.X()+50,pos.Z()+35,"") ;
      Text_t  line1[40] ;
      Text_t  line2[40] ;
      snprintf(line1,40,"Energy=%1.2f GeV",GetEnergy()) ;
      snprintf(line2,40,"%d Digits",GetDigitsMultiplicity()) ;
      clustertext ->AddText(line1) ;
      clustertext ->AddText(line2) ;
      clustertext ->Draw("");
    }
    gPad->Update() ; 
    Print("dummy") ;
    delete[] xi ; 
    delete[] zi ; 
   }
  
break;
  
  case kButton1Up:
    if (digitgraph) {
      delete digitgraph  ;
      digitgraph = 0 ;
    }
    if (clustertext) {
      delete clustertext ;
      clustertext = 0 ;
    }
    
    break;
    
  }
}
//____________________________________________________________________________
void AliPHOSRecPoint::EvalAll(TClonesArray * /* digits */) 
{
  //evaluates (if necessary) all RecPoint data members 

}

//____________________________________________________________________________
void AliPHOSRecPoint::EvalPHOSMod(AliPHOSDigit * digit) 
{
  // Returns the PHOS module in which the RecPoint is found

  if( fPHOSMod == 0){
  Int_t relid[4] ; 
 
  AliPHOSGeometry * phosgeom = (AliPHOSGeometry::GetInstance());

  phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
  fPHOSMod = relid[0];
  }
}

//____________________________________________________________________________
void AliPHOSRecPoint::GetGlobalPosition(TVector3 & gpos, TMatrixF & gmat) const
{
  // returns the position of the cluster in the global reference system of ALICE
  // and the uncertainty on this position

   AliPHOSGeometry * phosgeom = (AliPHOSGeometry::GetInstance());
   phosgeom->GetGlobalPHOS(this, gpos, gmat);
  
}

//______________________________________________________________________________
void AliPHOSRecPoint::Paint(Option_t *)
{
  // Paint this ALiRecPoint as a TMarker  with its current attributes
  
  TVector3 pos(0.,0.,0.)  ;
  GetLocalPosition(pos)   ;
  Coord_t x = pos.X()     ;
  Coord_t y = pos.Z()     ;
  Color_t markercolor = 1 ;
  Size_t  markersize = 1. ;
  Style_t markerstyle = 5 ;
  
  if (!gPad->IsBatch()) {
    gVirtualX->SetMarkerColor(markercolor) ;
    gVirtualX->SetMarkerSize (markersize)  ;
    gVirtualX->SetMarkerStyle(markerstyle) ;
  }
  gPad->SetAttMarkerPS(markercolor,markerstyle,markersize) ;
  gPad->PaintPolyMarker(1,&x,&y,"") ;
}
//______________________________________________________________________________
void AliPHOSRecPoint::GetLocalPosition(TVector3 & pos) const
{
  // returns the position of the cluster in the local reference system 
  // of the sub-detector
  
  pos = fLocPos;
}

//____________________________________________________________________________
void AliPHOSRecPoint::EvalLocal2TrackingCSTransform()
{
  //Evaluates local to "tracking" c.s. transformation (B.P.).
  //All evaluations should be completed before calling for this function.
  //See ALICE PPR Chapter 5 p.18 for "tracking" c.s. definition,
  //or just ask Jouri Belikov. :)

  if(IsEmc()) {
    SetVolumeId(AliGeomManager::LayerToVolUID(AliGeomManager::kPHOS1,GetPHOSMod()-1));
  }
  else
    SetVolumeId(AliGeomManager::LayerToVolUID(AliGeomManager::kPHOS2,GetPHOSMod()-1));

  Double_t lxyz[3] = {fLocPos.X(),fLocPos.Y(),fLocPos.Z()};
  Double_t txyz[3] = {0,0,0};
  Double_t dy;
  Double_t crystalShift;

  AliPHOSGeometry * phosgeom = AliPHOSGeometry::GetInstance();
  AliPHOSEMCAGeometry* geoEMCA = phosgeom->GetEMCAGeometry(); 

  //Calculate offset to crystal surface.
  //See fCrystalShift code in AliPHOSGeometry::Init()).

  const Float_t * inthermo = geoEMCA->GetInnerThermoHalfSize() ;
  const Float_t * strip    = geoEMCA->GetStripHalfSize() ;
  const Float_t * splate   = geoEMCA->GetSupportPlateHalfSize();
  const Float_t * crystal  = geoEMCA->GetCrystalHalfSize() ;
  const Float_t * pin      = geoEMCA->GetAPDHalfSize() ;
  const Float_t * preamp   = geoEMCA->GetPreampHalfSize() ;
  crystalShift = -inthermo[1]+strip[1]+splate[1]+crystal[1]-geoEMCA->GetAirGapLed()/2.+pin[1]+preamp[1] ;

  if(IsEmc()) {
    dy = crystalShift;
    lxyz[2] = -lxyz[2]; //Opposite z directions in EMC matrix and local frame!!!
  }
  else
    dy = phosgeom->GetCPVBoxSize(1)/2.; //center of CPV module

  lxyz[1] = lxyz[1] - dy;

  const TGeoHMatrix* tr2loc = GetTracking2LocalMatrix();
  if(!tr2loc) AliFatal(Form("No Tracking2LocalMatrix found for VolumeID=%d",GetVolumeId()));

  tr2loc->MasterToLocal(lxyz,txyz);
  SetX(txyz[0]); SetY(txyz[1]); SetZ(txyz[2]);

  if(AliLog::GetGlobalDebugLevel()>0) {
    TVector3 gpos; TMatrixF gmat;
    GetGlobalPosition(gpos,gmat);
    Float_t gxyz[3];
    GetGlobalXYZ(gxyz);
    TString emc;
    if(IsEmc()) 
      emc="EMC";
    else
      emc="CPV";
    AliInfo(Form("lCS-->(%.3f,%.3f,%.3f), tCS-->(%.3f,%.3f,%.3f), gCS-->(%.3f,%.3f,%.3f),  gCScalc-->(%.3f,%.3f,%.3f), module %d %s",
		 fLocPos.X(),fLocPos.Y(),fLocPos.Z(),
		 GetX(),GetY(),GetZ(),
		 gpos.X(),gpos.Y(),gpos.Z(),
		 gxyz[0],gxyz[1],gxyz[2],GetPHOSMod(),emc.Data()));
  }

}
