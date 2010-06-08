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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliTPCCorrection class                                                     //
//                                                                            //
// This class provides a general framework to deal with space point           //
// distortions. An correction class which inherits from here is for example   //
// AliTPCExBBShape or AliTPCExBTwist                                          //
//                                                                            //
// General functions are (for example):                                       //
//   CorrectPoint(x,roc) where x is the vector of inital positions in         //
//   cartesian coordinates and roc represents the Read Out chamber number     //
//   according to the offline naming convention. The vector x is overwritten  //
//   with the corrected coordinates.                                          //
//                                                                            //
// An alternative usage would be CorrectPoint(x,roc,dx), which leaves the     //
//   vector x untouched, put returns the distortions via the vector dx        //
//                                                                            //
// The class allows "effective Omega Tau" corrections to be shifted to the    //
// single distortion classes.                                                 //
//                                                                            //
// Note: This class is normally used via the class AliTPCComposedCorrection   //
//                                                                            //
// date: 27/04/2010                                                           //
// Authors: Magnus Mager, Stefan Rossegger, Jim Thomas                        //
////////////////////////////////////////////////////////////////////////////////
#include "Riostream.h"

#include <TH2F.h>
#include <TMath.h>
#include <TROOT.h>
#include <TTreeStream.h>
#include <TTree.h>
#include <TFile.h>
#include <TTimeStamp.h>
#include <AliCDBStorage.h>
#include <AliCDBId.h>
#include <AliCDBMetaData.h>
#include  "TVectorD.h"

#include "TRandom.h"
#include "AliTPCTransform.h"
#include "AliTPCcalibDB.h"
#include "AliTPCExB.h"
#include "AliTPCCorrection.h"
#include "AliTPCRecoParam.h"

#include  "AliExternalTrackParam.h"
#include  "AliTrackPointArray.h"
#include  "TDatabasePDG.h"
#include  "AliTrackerBase.h"
#include  "AliTPCROC.h"
#include  "THnSparse.h"
#include  "AliTPCLaserTrack.h"

#include "AliTPCCorrection.h"

ClassImp(AliTPCCorrection)

// FIXME: the following values should come from the database
const Double_t AliTPCCorrection::fgkTPCZ0    =249.7;     // nominal gating grid position 
const Double_t AliTPCCorrection::fgkIFCRadius= 83.06;    // Mean Radius of the Inner Field Cage ( 82.43 min,  83.70 max) (cm)
const Double_t AliTPCCorrection::fgkOFCRadius=254.5;     // Mean Radius of the Outer Field Cage (252.55 min, 256.45 max) (cm)
const Double_t AliTPCCorrection::fgkZOffSet  = 0.2;      // Offset from CE: calculate all distortions closer to CE as if at this point
const Double_t AliTPCCorrection::fgkCathodeV =-100000.0; // Cathode Voltage (volts)
const Double_t AliTPCCorrection::fgkGG       =-70.0;     // Gating Grid voltage (volts)


// FIXME: List of interpolation points (course grid in the middle, fine grid on the borders)
const Double_t AliTPCCorrection::fgkRList[AliTPCCorrection::kNR] =  {   
84.0,   84.5,   85.0,   85.5,   86.0,  87.0,    88.0,
90.0,   92.0,   94.0,   96.0,   98.0,  100.0,  102.0,  104.0,  106.0,  108.0, 
110.0,  112.0,  114.0,  116.0,  118.0,  120.0,  122.0,  124.0,  126.0,  128.0, 
130.0,  132.0,  134.0,  136.0,  138.0,  140.0,  142.0,  144.0,  146.0,  148.0, 
150.0,  152.0,  154.0,  156.0,  158.0,  160.0,  162.0,  164.0,  166.0,  168.0, 
170.0,  172.0,  174.0,  176.0,  178.0,  180.0,  182.0,  184.0,  186.0,  188.0,
190.0,  192.0,  194.0,  196.0,  198.0,  200.0,  202.0,  204.0,  206.0,  208.0,
210.0,  212.0,  214.0,  216.0,  218.0,  220.0,  222.0,  224.0,  226.0,  228.0,
230.0,  232.0,  234.0,  236.0,  238.0,  240.0,  242.0,  244.0,  246.0,  248.0,
249.0,  249.5,  250.0,  251.5,  252.0  } ;
  
const Double_t AliTPCCorrection::fgkZList[AliTPCCorrection::kNZ]     =   { 
-249.5, -249.0, -248.5, -248.0, -247.0, -246.0, -245.0, -243.0, -242.0, -241.0,
-240.0, -238.0, -236.0, -234.0, -232.0, -230.0, -228.0, -226.0, -224.0, -222.0,
-220.0, -218.0, -216.0, -214.0, -212.0, -210.0, -208.0, -206.0, -204.0, -202.0,
-200.0, -198.0, -196.0, -194.0, -192.0, -190.0, -188.0, -186.0, -184.0, -182.0,
-180.0, -178.0, -176.0, -174.0, -172.0, -170.0, -168.0, -166.0, -164.0, -162.0,
-160.0, -158.0, -156.0, -154.0, -152.0, -150.0, -148.0, -146.0, -144.0, -142.0,
-140.0, -138.0, -136.0, -134.0, -132.0, -130.0, -128.0, -126.0, -124.0, -122.0,
-120.0, -118.0, -116.0, -114.0, -112.0, -110.0, -108.0, -106.0, -104.0, -102.0,
-100.0,  -98.0,  -96.0,  -94.0,  -92.0,  -90.0,  -88.0,  -86.0,  -84.0,  -82.0,
-80.0,  -78.0,  -76.0,  -74.0,  -72.0,  -70.0,  -68.0,  -66.0,  -64.0,  -62.0,
-60.0,  -58.0,  -56.0,  -54.0,  -52.0,  -50.0,  -48.0,  -46.0,  -44.0,  -42.0,
-40.0,  -38.0,  -36.0,  -34.0,  -32.0,  -30.0,  -28.0,  -26.0,  -24.0,  -22.0,
-20.0,  -18.0,  -16.0,  -14.0,  -12.0,  -10.0,   -8.0,   -6.0,   -4.0,   -2.0,
-1.0,   -0.5,   -0.2,   -0.1,  -0.05,   0.05,    0.1,    0.2,    0.5,    1.0, 
 2.0,    4.0,    6.0,    8.0,   10.0,   12.0,   14.0,   16.0,   18.0,   20.0,
 22.0,   24.0,   26.0,   28.0,   30.0,   32.0,   34.0,   36.0,   38.0,   40.0, 
 42.0,   44.0,   46.0,   48.0,   50.0,   52.0,   54.0,   56.0,   58.0,   60.0, 
 62.0,   64.0,   66.0,   68.0,   70.0,   72.0,   74.0,   76.0,   78.0,   80.0, 
 82.0,   84.0,   86.0,   88.0,   90.0,   92.0,   94.0,   96.0,   98.0,  100.0, 
102.0,  104.0,  106.0,  108.0,  110.0,  112.0,  114.0,  116.0,  118.0,  120.0, 
122.0,  124.0,  126.0,  128.0,  130.0,  132.0,  134.0,  136.0,  138.0,  140.0, 
142.0,  144.0,  146.0,  148.0,  150.0,  152.0,  154.0,  156.0,  158.0,  160.0, 
162.0,  164.0,  166.0,  168.0,  170.0,  172.0,  174.0,  176.0,  178.0,  180.0, 
182.0,  184.0,  186.0,  188.0,  190.0,  192.0,  194.0,  196.0,  198.0,  200.0,
202.0,  204.0,  206.0,  208.0,  210.0,  212.0,  214.0,  216.0,  218.0,  220.0,
222.0,  224.0,  226.0,  228.0,  230.0,  232.0,  234.0,  236.0,  238.0,  240.0,
242.0,  243.0,  244.0,  245.0,  246.0,  247.0,  248.0,  248.5,  249.0,  249.5   } ;



AliTPCCorrection::AliTPCCorrection() 
  : TNamed("correction_unity","unity"),fJLow(0),fKLow(0), fT1(1), fT2(1)
{
  //
  // default constructor
  //
}

AliTPCCorrection::AliTPCCorrection(const char *name,const char *title)
: TNamed(name,title),fJLow(0),fKLow(0), fT1(1), fT2(1)
{
  //
  // default constructor, that set the name and title
  //
}

AliTPCCorrection::~AliTPCCorrection() {
  // 
  // virtual destructor
  //
}

void AliTPCCorrection::CorrectPoint(Float_t x[],const Short_t roc) {
  //
  // Corrects the initial coordinates x (cartesian coordinates)
  // according to the given effect (inherited classes)
  // roc represents the TPC read out chamber (offline numbering convention)
  //
  Float_t dx[3];
  GetCorrection(x,roc,dx);
  for (Int_t j=0;j<3;++j) x[j]+=dx[j];
}

void AliTPCCorrection::CorrectPoint(const Float_t x[],const Short_t roc,Float_t xp[]) {
  //
  // Corrects the initial coordinates x (cartesian coordinates) and stores the new 
  // (distorted) coordinates in xp. The distortion is set according to the given effect (inherited classes)
  // roc represents the TPC read out chamber (offline numbering convention)
  //
  Float_t dx[3];
  GetCorrection(x,roc,dx);
  for (Int_t j=0;j<3;++j) xp[j]=x[j]+dx[j];
}

void AliTPCCorrection::DistortPoint(Float_t x[],const Short_t roc) {
  //
  // Distorts the initial coordinates x (cartesian coordinates)
  // according to the given effect (inherited classes)
  // roc represents the TPC read out chamber (offline numbering convention)
  //
  Float_t dx[3];
  GetDistortion(x,roc,dx);
  for (Int_t j=0;j<3;++j) x[j]+=dx[j];
}

void AliTPCCorrection::DistortPoint(const Float_t x[],const Short_t roc,Float_t xp[]) {
  //
  // Distorts the initial coordinates x (cartesian coordinates) and stores the new 
  // (distorted) coordinates in xp. The distortion is set according to the given effect (inherited classes)
  // roc represents the TPC read out chamber (offline numbering convention)
  //
  Float_t dx[3];
  GetDistortion(x,roc,dx);
  for (Int_t j=0;j<3;++j) xp[j]=x[j]+dx[j];
}

void AliTPCCorrection::GetCorrection(const Float_t /*x*/[],const Short_t /*roc*/,Float_t dx[]) {
  //
  // This function delivers the correction values dx in respect to the inital coordinates x
  // roc represents the TPC read out chamber (offline numbering convention)
  // Note: The dx is overwritten by the inherited effectice class ...
  //
  for (Int_t j=0;j<3;++j) { dx[j]=0.; }
}

void AliTPCCorrection::GetDistortion(const Float_t x[],const Short_t roc,Float_t dx[]) {
  //
  // This function delivers the distortion values dx in respect to the inital coordinates x
  // roc represents the TPC read out chamber (offline numbering convention)
  //
  GetCorrection(x,roc,dx);
  for (Int_t j=0;j<3;++j) dx[j]=-dx[j];
}

void AliTPCCorrection::Init() {
  //
  // Initialization funtion (not used at the moment)
  //
}

void AliTPCCorrection::Update(const TTimeStamp &/*timeStamp*/) {
  //
  // Update function 
  //
}

void AliTPCCorrection::Print(Option_t* /*option*/) const {
  //
  // Print function to check which correction classes are used 
  // option=="d" prints details regarding the setted magnitude 
  // option=="a" prints the C0 and C1 coefficents for calibration purposes
  //
  printf("TPC spacepoint correction: \"%s\"\n",GetTitle());
}

void AliTPCCorrection:: SetOmegaTauT1T2(Float_t /*omegaTau*/,Float_t t1,Float_t t2) {
  //
  // Virtual funtion to pass the wt values (might become event dependent) to the inherited classes
  // t1 and t2 represent the "effective omegaTau" corrections and were measured in a dedicated
  // calibration run
  //
  fT1=t1;
  fT2=t2;
  //SetOmegaTauT1T2(omegaTau, t1, t2);
}

TH2F* AliTPCCorrection::CreateHistoDRinXY(Float_t z,Int_t nx,Int_t ny) {
  //
  // Simple plot functionality.
  // Returns a 2d hisogram which represents the corrections in radial direction (dr)
  // in respect to position z within the XY plane.
  // The histogramm has nx times ny entries. 
  //
  
  TH2F *h=CreateTH2F("dr_xy",GetTitle(),"x [cm]","y [cm]","dr [cm]",
		     nx,-250.,250.,ny,-250.,250.);
  Float_t x[3],dx[3];
  x[2]=z;
  Int_t roc=z>0.?0:18; // FIXME
  for (Int_t iy=1;iy<=ny;++iy) {
    x[1]=h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix=1;ix<=nx;++ix) {
      x[0]=h->GetXaxis()->GetBinCenter(ix);
      GetCorrection(x,roc,dx);
      Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));
      if (90.<=r0 && r0<=250.) {
	Float_t r1=TMath::Sqrt((x[0]+dx[0])*(x[0]+dx[0])+(x[1]+dx[1])*(x[1]+dx[1]));
	h->SetBinContent(ix,iy,r1-r0);
      }
      else
	h->SetBinContent(ix,iy,0.);
    }
  }
  return h;
}

TH2F* AliTPCCorrection::CreateHistoDRPhiinXY(Float_t z,Int_t nx,Int_t ny) {
  //
  // Simple plot functionality.
  // Returns a 2d hisogram which represents the corrections in rphi direction (drphi) 
  // in respect to position z within the XY plane.
  // The histogramm has nx times ny entries. 
  //

  TH2F *h=CreateTH2F("drphi_xy",GetTitle(),"x [cm]","y [cm]","drphi [cm]",
		     nx,-250.,250.,ny,-250.,250.);
  Float_t x[3],dx[3];
  x[2]=z;
  Int_t roc=z>0.?0:18; // FIXME
  for (Int_t iy=1;iy<=ny;++iy) {
    x[1]=h->GetYaxis()->GetBinCenter(iy);
    for (Int_t ix=1;ix<=nx;++ix) {
      x[0]=h->GetXaxis()->GetBinCenter(ix);
      GetCorrection(x,roc,dx);
      Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));
      if (90.<=r0 && r0<=250.) {
 	Float_t phi0=TMath::ATan2(x[1]      ,x[0]      );
	Float_t phi1=TMath::ATan2(x[1]+dx[1],x[0]+dx[0]);

	Float_t dphi=phi1-phi0;
	if (dphi<TMath::Pi()) dphi+=TMath::TwoPi();
	if (dphi>TMath::Pi()) dphi-=TMath::TwoPi();
      
	h->SetBinContent(ix,iy,r0*dphi);
      }
      else
	h->SetBinContent(ix,iy,0.);
    }
  }
  return h;
}

TH2F* AliTPCCorrection::CreateHistoDRinZR(Float_t phi,Int_t nz,Int_t nr) {
  //
  // Simple plot functionality.
  // Returns a 2d hisogram which represents the corrections in r direction (dr) 
  // in respect to angle phi within the ZR plane.
  // The histogramm has nx times ny entries. 
  //
  TH2F *h=CreateTH2F("dr_zr",GetTitle(),"z [cm]","r [cm]","dr [cm]",
		     nz,-250.,250.,nr,85.,250.);
  Float_t x[3],dx[3];
  for (Int_t ir=1;ir<=nr;++ir) {
    Float_t radius=h->GetYaxis()->GetBinCenter(ir);
    x[0]=radius*TMath::Cos(phi);
    x[1]=radius*TMath::Sin(phi);
    for (Int_t iz=1;iz<=nz;++iz) {
      x[2]=h->GetXaxis()->GetBinCenter(iz);
      Int_t roc=x[2]>0.?0:18; // FIXME
      GetCorrection(x,roc,dx);
      Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));
      Float_t r1=TMath::Sqrt((x[0]+dx[0])*(x[0]+dx[0])+(x[1]+dx[1])*(x[1]+dx[1]));
      h->SetBinContent(iz,ir,r1-r0);
    }
  }
  printf("SDF\n");
  return h;

}

TH2F* AliTPCCorrection::CreateHistoDRPhiinZR(Float_t phi,Int_t nz,Int_t nr) {
  //
  // Simple plot functionality.
  // Returns a 2d hisogram which represents the corrections in rphi direction (drphi) 
  // in respect to angle phi within the ZR plane.
  // The histogramm has nx times ny entries. 
  //
  TH2F *h=CreateTH2F("drphi_zr",GetTitle(),"z [cm]","r [cm]","drphi [cm]",
		     nz,-250.,250.,nr,85.,250.);
  Float_t x[3],dx[3];
  for (Int_t iz=1;iz<=nz;++iz) {
    x[2]=h->GetXaxis()->GetBinCenter(iz);
    Int_t roc=x[2]>0.?0:18; // FIXME
    for (Int_t ir=1;ir<=nr;++ir) {
      Float_t radius=h->GetYaxis()->GetBinCenter(ir);
      x[0]=radius*TMath::Cos(phi);
      x[1]=radius*TMath::Sin(phi);
      GetCorrection(x,roc,dx);
      Float_t r0=TMath::Sqrt((x[0]      )*(x[0]      )+(x[1]      )*(x[1]      ));
      Float_t phi0=TMath::ATan2(x[1]      ,x[0]      );
      Float_t phi1=TMath::ATan2(x[1]+dx[1],x[0]+dx[0]);
      
      Float_t dphi=phi1-phi0;
      if (dphi<TMath::Pi()) dphi+=TMath::TwoPi();
      if (dphi>TMath::Pi()) dphi-=TMath::TwoPi();
      
      h->SetBinContent(iz,ir,r0*dphi);
    }
  }
  return h;
}

TH2F* AliTPCCorrection::CreateTH2F(const char *name,const char *title,
				   const char *xlabel,const char *ylabel,const char *zlabel,
				  Int_t nbinsx,Double_t xlow,Double_t xup,
				  Int_t nbinsy,Double_t ylow,Double_t yup) {
  //
  // Helper function to create a 2d histogramm of given size
  //
  
  TString hname=name;
  Int_t i=0;
  if (gDirectory) {
    while (gDirectory->FindObject(hname.Data())) {
      hname =name;
      hname+="_";
      hname+=i;
      ++i;
    }
  }
  TH2F *h=new TH2F(hname.Data(),title,
		   nbinsx,xlow,xup,
		   nbinsy,ylow,yup);
  h->GetXaxis()->SetTitle(xlabel);
  h->GetYaxis()->SetTitle(ylabel);
  h->GetZaxis()->SetTitle(zlabel);
  h->SetStats(0);
  return h;
}


// Simple Interpolation functions: e.g. with bi(tri)cubic interpolations (not yet in TH2 and TH3)

void AliTPCCorrection::Interpolate2DEdistortion( const Int_t order, const Double_t r, const Double_t z, 
						  const Double_t er[kNZ][kNR], Double_t &erValue ) {
  //
  // Interpolate table - 2D interpolation
  //
  Double_t saveEr[10] ;

  Search( kNZ,   fgkZList,  z,   fJLow   ) ;
  Search( kNR,   fgkRList,  r,   fKLow   ) ;
  if ( fJLow < 0 ) fJLow = 0 ;   // check if out of range
  if ( fKLow < 0 ) fKLow = 0 ;
  if ( fJLow + order  >=    kNZ - 1 ) fJLow =   kNZ - 1 - order ;
  if ( fKLow + order  >=    kNR - 1 ) fKLow =   kNR - 1 - order ;

  for ( Int_t j = fJLow ; j < fJLow + order + 1 ; j++ ) {
      saveEr[j-fJLow]     = Interpolate( &fgkRList[fKLow], &er[j][fKLow], order, r )   ;
  }
  erValue = Interpolate( &fgkZList[fJLow], saveEr, order, z )   ;

}


Double_t AliTPCCorrection::Interpolate( const Double_t xArray[], const Double_t yArray[], 
				       const Int_t order, const Double_t x ) {
  //
  // Interpolate function Y(x) using linear (order=1) or quadratic (order=2) interpolation.
  //

  Double_t y ;
  if ( order == 2 ) {                // Quadratic Interpolation = 2 
    y  = (x-xArray[1]) * (x-xArray[2]) * yArray[0] / ( (xArray[0]-xArray[1]) * (xArray[0]-xArray[2]) ) ; 
    y += (x-xArray[2]) * (x-xArray[0]) * yArray[1] / ( (xArray[1]-xArray[2]) * (xArray[1]-xArray[0]) ) ; 
    y += (x-xArray[0]) * (x-xArray[1]) * yArray[2] / ( (xArray[2]-xArray[0]) * (xArray[2]-xArray[1]) ) ; 
  } else {                           // Linear Interpolation = 1
    y  = yArray[0] + ( yArray[1]-yArray[0] ) * ( x-xArray[0] ) / ( xArray[1] - xArray[0] ) ;
  }

  return (y);

}


void AliTPCCorrection::Search( const Int_t n, const Double_t xArray[], const Double_t x, Int_t &low ) {
  //
  // Search an ordered table by starting at the most recently used point
  //

  Long_t middle, high ;
  Int_t  ascend = 0, increment = 1 ;

  if ( xArray[n-1] >= xArray[0] ) ascend = 1 ;  // Ascending ordered table if true
  
  if ( low < 0 || low > n-1 ) { 
    low = -1 ; high = n ; 
  } else {                                            // Ordered Search phase
    if ( (Int_t)( x >= xArray[low] ) == ascend )  {
      if ( low == n-1 ) return ;          
      high = low + 1 ;
      while ( (Int_t)( x >= xArray[high] ) == ascend ) {
	low = high ;
	increment *= 2 ;
	high = low + increment ;
	if ( high > n-1 )  {  high = n ; break ;  }
      }
    } else {
      if ( low == 0 )  {  low = -1 ;  return ;  }
      high = low - 1 ;
      while ( (Int_t)( x < xArray[low] ) == ascend ) {
	high = low ;
	increment *= 2 ;
	if ( increment >= high )  {  low = -1 ;  break ;  }
	else  low = high - increment ;
      }
    }
  }
  
  while ( (high-low) != 1 ) {                     // Binary Search Phase
    middle = ( high + low ) / 2 ;
    if ( (Int_t)( x >= xArray[middle] ) == ascend )
      low = middle ;
    else
      high = middle ;
  }
  
  if ( x == xArray[n-1] ) low = n-2 ;
  if ( x == xArray[0]   ) low = 0 ;
  
}


AliExternalTrackParam * AliTPCCorrection::FitDistortedTrack(AliExternalTrackParam & trackIn, Double_t refX, Int_t dir, TTreeSRedirector * const pcstream){
  //
  // Fit the track parameters - without and with distortion
  // 1. Space points in the TPC are simulated along the trajectory  
  // 2. Space points distorted
  // 3. Fits  the non distorted and distroted track to the reference plane at refX
  // 4. For visualization and debugging  purposes the space points and tracks can be stored  in the tree - using the TTreeSRedirector functionality 
  //
  // trackIn   - input track parameters
  // refX     - reference X to fit the track
  // dir      - direction - out=1 or in=-1
  // pcstream -  debug streamer to check the results
  //
  // see AliExternalTrackParam.h documentation:
  // track1.fP[0] - local y (rphi)
  // track1.fP[1] - z 
  // track1.fP[2] - sinus of local inclination angle
  // track1.fP[3] - tangent of deep angle
  // track1.fP[4] - 1/pt
  AliTPCROC * roc = AliTPCROC::Instance();
  const Int_t    npoints0=roc->GetNRows(0)+roc->GetNRows(36);
  const Double_t kRTPC0  =roc->GetPadRowRadii(0,0);
  const Double_t kRTPC1  =roc->GetPadRowRadii(36,roc->GetNRows(36)-1);
    
  const Double_t kMaxSnp = 0.85;  
  const Double_t kSigmaY=0.1;
  const Double_t kSigmaZ=0.1;
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();

  AliExternalTrackParam  track(trackIn); // 
  // generate points
  AliTrackPointArray pointArray0(npoints0);
  AliTrackPointArray pointArray1(npoints0);
  Double_t xyz[3];
  AliTrackerBase::PropagateTrackToBxByBz(&track,kRTPC0,kMass,3,kTRUE,kMaxSnp);
  //
  // simulate the track
  Int_t npoints=0;
  Float_t covPoint[6]={0,0,0, kSigmaY*kSigmaY,0,kSigmaZ*kSigmaZ};  //covariance at the local frame
  for (Double_t radius=kRTPC0; radius<kRTPC1; radius++){
    AliTrackerBase::PropagateTrackToBxByBz(&track,radius,kMass,3,kTRUE,kMaxSnp);
    track.GetXYZ(xyz);
    xyz[0]+=gRandom->Gaus(0,0.005);
    xyz[1]+=gRandom->Gaus(0,0.005);
    xyz[2]+=gRandom->Gaus(0,0.005);
    AliTrackPoint pIn0;                               // space point          
    AliTrackPoint pIn1;
    Int_t sector= (xyz[2]>0)? 0:18;
    pointArray0.GetPoint(pIn0,npoints);
    pointArray1.GetPoint(pIn1,npoints);
    Double_t alpha = TMath::ATan2(xyz[1],xyz[0]);
    Float_t distPoint[3]={xyz[0],xyz[1],xyz[2]};
    DistortPoint(distPoint, sector);
    pIn0.SetXYZ(xyz[0], xyz[1],xyz[2]);
    pIn1.SetXYZ(distPoint[0], distPoint[1],distPoint[2]);
    //
    track.Rotate(alpha);
    AliTrackPoint prot0 = pIn0.Rotate(alpha);   // rotate to the local frame - non distoted  point
    AliTrackPoint prot1 = pIn1.Rotate(alpha);   // rotate to the local frame -     distorted point
    prot0.SetXYZ(prot0.GetX(),prot0.GetY(), prot0.GetZ(),covPoint);
    prot1.SetXYZ(prot1.GetX(),prot1.GetY(), prot1.GetZ(),covPoint);
    pIn0=prot0.Rotate(-alpha);                       // rotate back to global frame
    pIn1=prot1.Rotate(-alpha);                       // rotate back to global frame
    pointArray0.AddPoint(npoints, &pIn0);      
    pointArray1.AddPoint(npoints, &pIn1);
    npoints++;
    if (npoints>=npoints0) break;
  }
  if (npoints<npoints0/2) return 0;
  //
  // refit track
  //
  AliExternalTrackParam *track0=0;
  AliExternalTrackParam *track1=0;
  AliTrackPoint   point1,point2,point3;
  if (dir==1) {  //make seed inner
    pointArray0.GetPoint(point1,1);
    pointArray0.GetPoint(point2,30);
    pointArray0.GetPoint(point3,60);
  }
  if (dir==-1){ //make seed outer
    pointArray0.GetPoint(point1,npoints-60);
    pointArray0.GetPoint(point2,npoints-30);
    pointArray0.GetPoint(point3,npoints-1);
  }  
  track0 = AliTrackerBase::MakeSeed(point1, point2, point3);
  track1 = AliTrackerBase::MakeSeed(point1, point2, point3);

  for (Int_t jpoint=0; jpoint<npoints; jpoint++){
    Int_t ipoint= (dir>0) ? jpoint: npoints-1-jpoint;
    //
    AliTrackPoint pIn0;
    AliTrackPoint pIn1;
    pointArray0.GetPoint(pIn0,ipoint);
    pointArray1.GetPoint(pIn1,ipoint);
    AliTrackPoint prot0 = pIn0.Rotate(track0->GetAlpha());   // rotate to the local frame - non distoted  point
    AliTrackPoint prot1 = pIn1.Rotate(track1->GetAlpha());   // rotate to the local frame -     distorted point
    //
    AliTrackerBase::PropagateTrackToBxByBz(track0,prot0.GetX(),kMass,3,kFALSE,kMaxSnp);
    AliTrackerBase::PropagateTrackToBxByBz(track1,prot0.GetX(),kMass,3,kFALSE,kMaxSnp);
    track.GetXYZ(xyz);  // distorted track also propagated to the same reference radius
    //
    Double_t pointPos[2]={0,0};
    Double_t pointCov[3]={0,0,0};
    pointPos[0]=prot0.GetY();//local y
    pointPos[1]=prot0.GetZ();//local z
    pointCov[0]=prot0.GetCov()[3];//simay^2
    pointCov[1]=prot0.GetCov()[4];//sigmayz
    pointCov[2]=prot0.GetCov()[5];//sigmaz^2
    track0->Update(pointPos,pointCov);
    //
    Double_t deltaX=prot1.GetX()-prot0.GetX();   // delta X 
    Double_t deltaYX=deltaX*TMath::Tan(TMath::ASin(track1->GetSnp()));  // deltaY due  delta X
    Double_t deltaZX=deltaX*track1->GetTgl();                           // deltaZ due  delta X

    pointPos[0]=prot1.GetY()-deltaYX;//local y is sign correct? should be minus
    pointPos[1]=prot1.GetZ()-deltaZX;//local z is sign correct? should be minus
    pointCov[0]=prot1.GetCov()[3];//simay^2
    pointCov[1]=prot1.GetCov()[4];//sigmayz
    pointCov[2]=prot1.GetCov()[5];//sigmaz^2
    track1->Update(pointPos,pointCov);
  }

  AliTrackerBase::PropagateTrackToBxByBz(track0,refX,kMass,2.,kTRUE,kMaxSnp);
  track1->Rotate(track0->GetAlpha());
  track1->PropagateTo(track0->GetX(),AliTrackerBase::GetBz());

  if (pcstream) (*pcstream)<<Form("fitDistort%s",GetName())<<
    "point0.="<<&pointArray0<<   //  points
    "point1.="<<&pointArray1<<   //  distorted points
    "trackIn.="<<&track<<       //  original track
    "track0.="<<track0<<         //  fitted track
    "track1.="<<track1<<         //  fitted distorted track
    "\n";
  new(&trackIn) AliExternalTrackParam(*track0);
  delete track0;
  return track1;
}





TTree* AliTPCCorrection::CreateDistortionTree(Double_t step){
  //
  // create the distortion tree on a mesh with granularity given by step
  // return the tree with distortions at given position 
  // Map is created on the mesh with given step size
  //
  TTreeSRedirector *pcstream = new TTreeSRedirector(Form("correction%s.root",GetName()));
  Float_t xyz[3];
  for (Double_t x= -250; x<250; x+=step){
    for (Double_t y= -250; y<250; y+=step){
      Double_t r    = TMath::Sqrt(x*x+y*y);
      if (r<80) continue;
      if (r>250) continue;      
      for (Double_t z= -250; z<250; z+=step){
	Int_t roc=(z>0)?0:18;
	xyz[0]=x;
	xyz[1]=y;
	xyz[2]=z;
	Double_t phi  = TMath::ATan2(y,x);
	DistortPoint(xyz,roc);
	Double_t r1    = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
	Double_t phi1  = TMath::ATan2(xyz[1],xyz[0]);
	if ((phi1-phi)>TMath::Pi()) phi1-=TMath::Pi();
	if ((phi1-phi)<-TMath::Pi()) phi1+=TMath::Pi();
	Double_t dx = xyz[0]-x;
	Double_t dy = xyz[1]-y;
	Double_t dz = xyz[2]-z;
	Double_t dr=r1-r;
	Double_t drphi=(phi1-phi)*r;
	(*pcstream)<<"distortion"<<
	  "x="<<x<<           // original position        
	  "y="<<y<<
	  "z="<<z<<
	  "r="<<r<<
	  "phi="<<phi<<	  
	  "x1="<<xyz[0]<<      // distorted position
	  "y1="<<xyz[1]<<
	  "z1="<<xyz[2]<<
	  "r1="<<r1<<
	  "phi1="<<phi1<<
	  //
	  "dx="<<dx<<          // delta position
	  "dy="<<dy<<
	  "dz="<<dz<<
	  "dr="<<dr<<
	  "drphi="<<drphi<<
	  "\n";
      }
    }	
  }
  delete pcstream;
  TFile f(Form("correction%s.root",GetName()));
  TTree * tree = (TTree*)f.Get("distortion");
  TTree * tree2= tree->CopyTree("1");
  tree2->SetName(Form("dist%s",GetName()));
  tree2->SetDirectory(0);
  delete tree;
  return tree2;
}




void AliTPCCorrection::MakeTrackDistortionTree(TTree *tinput, Int_t dtype, Int_t ptype, const TObjArray * corrArray, Int_t step, Bool_t debug ){
  //
  // Make a fit tree:
  // For each partial correction (specified in array) and given track topology (phi, theta, snp, refX)
  // calculates partial distortions
  // Partial distortion is stored in the resulting tree
  // Output is storred in the file distortion_<dettype>_<partype>.root
  // Partial  distortion is stored with the name given by correction name
  //
  //
  // Parameters of function:
  // input     - input tree
  // dtype     - distortion type 0 - ITSTPC,  1 -TPCTRD, 2 - TPCvertex 
  // ppype     - parameter type
  // corrArray - array with partial corrections
  // step      - skipe entries  - if 1 all entries processed - it is slow
  // debug     0 if debug on also space points dumped - it is slow
  const Double_t kMaxSnp = 0.85;  
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  //  const Double_t kB2C=-0.299792458e-3;
  const Int_t kMinEntries=50; 
  Double_t phi,theta, snp, mean,rms, entries;
  tinput->SetBranchAddress("theta",&theta);
  tinput->SetBranchAddress("phi", &phi);
  tinput->SetBranchAddress("snp",&snp);
  tinput->SetBranchAddress("mean",&mean);
  tinput->SetBranchAddress("rms",&rms);
  tinput->SetBranchAddress("entries",&entries);
  TTreeSRedirector *pcstream = new TTreeSRedirector(Form("distortion%d_%d.root",dtype,ptype));
  //
  Int_t nentries=tinput->GetEntries();
  Int_t ncorr=corrArray->GetEntries();
  Double_t corrections[100]={0}; //
  Double_t tPar[5];
  Double_t cov[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t refX=0;
  Int_t dir=0;
  if (dtype==0) {refX=85.; dir=-1;}
  if (dtype==1) {refX=275.; dir=1;}
  if (dtype==2) {refX=85.; dir=-1;}
  if (dtype==3) {refX=360.; dir=-1;}
  //
  for (Int_t ientry=0; ientry<nentries; ientry+=step){
    tinput->GetEntry(ientry);
    if (TMath::Abs(snp)>kMaxSnp) continue;
    tPar[0]=0;
    tPar[1]=theta*refX;
    tPar[2]=snp;
    tPar[3]=theta;
    tPar[4]=(gRandom->Rndm()-0.5)*0.02;  // should be calculated - non equal to 0
    Double_t bz=AliTrackerBase::GetBz();
    if (refX>10. && TMath::Abs(bz)>0.1 )  tPar[4]=snp/(refX*bz*kB2C*2);
    tPar[4]+=(gRandom->Rndm()-0.5)*0.02;
    AliExternalTrackParam track(refX,phi,tPar,cov);
    Double_t xyz[3];
    track.GetXYZ(xyz);
    Int_t id=0;
    Double_t dRrec=0; // dummy value - needed for points - e.g for laser
    if (ptype==4 &&bz<0) mean*=-1;  // interpret as curvature
    (*pcstream)<<"fit"<<
      "bz="<<bz<<         // magnetic filed used
      "dtype="<<dtype<<   // detector match type
      "ptype="<<ptype<<   // parameter type
      "theta="<<theta<<   // theta
      "phi="<<phi<<       // phi 
      "snp="<<snp<<       // snp
      "mean="<<mean<<     // mean dist value
      "rms="<<rms<<       // rms
      "gx="<<xyz[0]<<         // global position at reference
      "gy="<<xyz[1]<<         // global position at reference
      "gz="<<xyz[2]<<         // global position at reference	
      "dRrec="<<dRrec<<      // delta Radius in reconstruction
      "id="<<id<<             // track id
      "entries="<<entries;// number of entries in bin
    //
    for (Int_t icorr=0; icorr<ncorr; icorr++) {
      AliTPCCorrection *corr = (AliTPCCorrection*)corrArray->At(icorr);
      corrections[icorr]=0;
      if (entries>kMinEntries){
	AliExternalTrackParam trackIn(refX,phi,tPar,cov);
	AliExternalTrackParam *trackOut = 0;
	if (debug) trackOut=corr->FitDistortedTrack(trackIn, refX, dir,pcstream);
	if (!debug) trackOut=corr->FitDistortedTrack(trackIn, refX, dir,0);
	if (dtype==0) {refX=85.; dir=-1;}
	if (dtype==1) {refX=275.; dir=1;}
	if (dtype==2) {refX=0; dir=-1;}
	if (dtype==3) {refX=360.; dir=-1;}
	//
	if (trackOut){
	  AliTrackerBase::PropagateTrackToBxByBz(&trackIn,refX,kMass,3,kTRUE,kMaxSnp);
	  trackOut->Rotate(trackIn.GetAlpha());
	  trackOut->PropagateTo(trackIn.GetX(),AliTrackerBase::GetBz());
	  //
	  corrections[icorr]= trackOut->GetParameter()[ptype]-trackIn.GetParameter()[ptype];
	  delete trackOut;      
	}else{
	  corrections[icorr]=0;
	}
	if (ptype==4 &&bz<0) corrections[icorr]*=-1;  // interpret as curvature
      }      
      Double_t dRdummy=0;
      (*pcstream)<<"fit"<<
	Form("%s=",corr->GetName())<<corrections[icorr]<<   // dump correction value
	Form("dR%s=",corr->GetName())<<dRdummy;   // dump dummy correction value not needed for tracks 
                                                  // for points it is neccessary
    }
    (*pcstream)<<"fit"<<"\n";
  }
  delete pcstream;
}



void AliTPCCorrection::MakeLaserDistortionTree(TTree* tree, TObjArray *corrArray, Int_t itype){
  //
  // Make a laser fit tree for global minimization
  //
  const Double_t cutErrY=0.1;
  const Double_t cutErrZ=0.1;
  const Double_t kEpsilon=0.00000001;
  TVectorD *vecdY=0;
  TVectorD *vecdZ=0;
  TVectorD *veceY=0;
  TVectorD *veceZ=0;
  AliTPCLaserTrack *ltr=0;
  AliTPCLaserTrack::LoadTracks();
  tree->SetBranchAddress("dY.",&vecdY);
  tree->SetBranchAddress("dZ.",&vecdZ);
  tree->SetBranchAddress("eY.",&veceY);
  tree->SetBranchAddress("eZ.",&veceZ);
  tree->SetBranchAddress("LTr.",&ltr);
  Int_t entries= tree->GetEntries();
  TTreeSRedirector *pcstream= new TTreeSRedirector("distortion4_0.root");
  Double_t bz=AliTrackerBase::GetBz();
  // 

  for (Int_t ientry=0; ientry<entries; ientry++){
    tree->GetEntry(ientry);
    if (!ltr->GetVecGX()){
      ltr->UpdatePoints();
    }
    TVectorD * delta= (itype==0)? vecdY:vecdZ;
    TVectorD * err= (itype==0)? veceY:veceZ;
    
    for (Int_t irow=0; irow<159; irow++){
      Int_t nentries = 1000;
      if (veceY->GetMatrixArray()[irow]>cutErrY||veceZ->GetMatrixArray()[irow]>cutErrZ) nentries=0;
      if (veceY->GetMatrixArray()[irow]<kEpsilon||veceZ->GetMatrixArray()[irow]<kEpsilon) nentries=0;
      Int_t dtype=4;
      Double_t phi   =(*ltr->GetVecPhi())[irow];
      Double_t theta =ltr->GetTgl();
      Double_t mean=delta->GetMatrixArray()[irow];
      Double_t gx=0,gy=0,gz=0;
      Double_t snp = (*ltr->GetVecP2())[irow];
      Double_t rms = 0.1+err->GetMatrixArray()[irow];
      gx = (*ltr->GetVecGX())[irow];
      gy = (*ltr->GetVecGY())[irow];
      gz = (*ltr->GetVecGZ())[irow];
      Int_t bundle= ltr->GetBundle();
      Double_t dRrec=0;
      //
      // get delta R used in reconstruction
      AliTPCcalibDB*  calib=AliTPCcalibDB::Instance();  
      AliTPCCorrection * correction = calib->GetTPCComposedCorrection();
      const AliTPCRecoParam * recoParam = calib->GetTransform()->GetCurrentRecoParam();
      Double_t xyz0[3]={gx,gy,gz};
      Double_t oldR=TMath::Sqrt(gx*gx+gy*gy);
      //
      // old ExB correction 
      //      
      if(recoParam&&recoParam->GetUseExBCorrection()) {	
	Double_t xyz1[3]={gx,gy,gz};
	calib->GetExB()->Correct(xyz0,xyz1);
	Double_t newR=TMath::Sqrt(xyz1[0]*xyz1[0]+xyz1[1]*xyz1[1]);
	dRrec=oldR-newR;
      } 
      if(recoParam&&recoParam->GetUseComposedCorrection()&&correction) {
	Float_t xyz1[3]={gx,gy,gz};
	Int_t sector=(gz>0)?0:18;
	correction->CorrectPoint(xyz1, sector);
	Double_t newR=TMath::Sqrt(xyz1[0]*xyz1[0]+xyz1[1]*xyz1[1]);
	dRrec=oldR-newR;
      } 


      (*pcstream)<<"fit"<<
	"bz="<<bz<<         // magnetic filed used
	"dtype="<<dtype<<   // detector match type
	"ptype="<<itype<<   // parameter type
	"theta="<<theta<<   // theta
	"phi="<<phi<<       // phi 
	"snp="<<snp<<       // snp
	"mean="<<mean<<     // mean dist value
	"rms="<<rms<<       // rms
	"gx="<<gx<<         // global position
	"gy="<<gy<<         // global position
	"gz="<<gz<<         // global position
	"dRrec="<<dRrec<<      // delta Radius in reconstruction
	"id="<<bundle<<     //bundle
	"entries="<<nentries;// number of entries in bin
      //
      //    
      Double_t ky = TMath::Tan(TMath::ASin(snp));
      Int_t ncorr = corrArray->GetEntries();
      Double_t r0   = TMath::Sqrt(gx*gx+gy*gy);
      Double_t phi0 = TMath::ATan2(gy,gx);
      Double_t distortions[1000]={0};
      Double_t distortionsR[1000]={0};
      for (Int_t icorr=0; icorr<ncorr; icorr++) {
	AliTPCCorrection *corr = (AliTPCCorrection*)corrArray->At(icorr);
	Float_t distPoint[3]={gx,gy,gz}; 
	Int_t sector= (gz>0)? 0:18;
	if (r0>80){
	  corr->DistortPoint(distPoint, sector);
	}
	Double_t value=distPoint[2]-gz;
	if (itype==0){
	  Double_t r1   = TMath::Sqrt(distPoint[0]*distPoint[0]+distPoint[1]*distPoint[1]);
	  Double_t phi1 = TMath::ATan2(distPoint[1],distPoint[0]);
	  Double_t drphi= r0*(phi1-phi0);
	  Double_t dr   = r1-r0;
	  distortions[icorr]  = drphi-ky*dr;
	  distortionsR[icorr] = dr;
	}
	(*pcstream)<<"fit"<<
	  Form("%s=",corr->GetName())<<distortions[icorr]<<    // dump correction value
	  Form("dR%s=",corr->GetName())<<distortionsR[icorr];   // dump correction R  value
      }
      (*pcstream)<<"fit"<<"\n";
    }
  }
  delete pcstream;
}



void   AliTPCCorrection::MakeDistortionMap(THnSparse * his0, TTreeSRedirector * const pcstream, const char* hname, Int_t run){
  //
  // make a distortion map out ou fthe residual histogram
  // Results are written to the debug streamer - pcstream
  // Parameters:
  //   his0       - input (4D) residual histogram
  //   pcstream   - file to write the tree
  //   run        - run number
  // marian.ivanov@cern.ch
  const Int_t kMinEntries=50;
  Int_t nbins1=his0->GetAxis(1)->GetNbins();
  Int_t first1=his0->GetAxis(1)->GetFirst();
  Int_t last1 =his0->GetAxis(1)->GetLast();
  //
  Double_t bz=AliTrackerBase::GetBz();
  Int_t idim[4]={0,1,2,3};
  for (Int_t ibin1=first1; ibin1<last1; ibin1++){   //axis 1 - theta
    //
    his0->GetAxis(1)->SetRange(TMath::Max(ibin1,1),TMath::Min(ibin1,nbins1));
    Double_t       x1= his0->GetAxis(1)->GetBinCenter(ibin1);
    THnSparse * his1 = his0->Projection(4,idim);  // projected histogram according range1
    Int_t nbins3     = his1->GetAxis(3)->GetNbins();
    Int_t first3     = his1->GetAxis(3)->GetFirst();
    Int_t last3      = his1->GetAxis(3)->GetLast();
    //

    for (Int_t ibin3=first3-1; ibin3<last3; ibin3+=1){   // axis 3 - local angle
      his1->GetAxis(3)->SetRange(TMath::Max(ibin3-1,1),TMath::Min(ibin3+1,nbins3));
      Double_t      x3= his1->GetAxis(3)->GetBinCenter(ibin3);
      if (ibin3<first3) {
	his1->GetAxis(3)->SetRangeUser(-1,1);
	x3=0;
      }
      THnSparse * his3= his1->Projection(4,idim);         //projected histogram according selection 3
      Int_t nbins2    = his3->GetAxis(2)->GetNbins();
      Int_t first2    = his3->GetAxis(2)->GetFirst();
      Int_t last2     = his3->GetAxis(2)->GetLast();
      //
      for (Int_t ibin2=first2; ibin2<last2; ibin2+=1){
	his3->GetAxis(2)->SetRange(TMath::Max(ibin2-1,1),TMath::Min(ibin2+1,nbins2));
	Double_t x2= his3->GetAxis(2)->GetBinCenter(ibin2);
	TH1 * hisDelta = his3->Projection(0);
	//
	Double_t entries = hisDelta->GetEntries();
	Double_t mean=0, rms=0;
	if (entries>kMinEntries){
	  mean    = hisDelta->GetMean(); 
	  rms = hisDelta->GetRMS(); 
	}
	(*pcstream)<<hname<<
	  "run="<<run<<
	  "bz="<<bz<<
	  "theta="<<x1<<
	  "phi="<<x2<<
	  "snp="<<x3<<
	  "entries="<<entries<<
	  "mean="<<mean<<
	  "rms="<<rms<<
	  "\n";
	delete hisDelta;
	printf("%f\t%f\t%f\t%f\t%f\n",x1,x3,x2, entries,mean);
      }
      delete his3;
    }
    delete his1;
  }
}





void AliTPCCorrection::StoreInOCDB(Int_t startRun, Int_t endRun, const char *comment){
  //
  // Store object in the OCDB
  // By default the object is stored in the current directory 
  // default comment consit of user name and the date
  //
  TString ocdbStorage="";
  ocdbStorage+="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("AliTPCCorrection");
  metaData->SetResponsible("Marian Ivanov");
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion("05-25-01"); //root version
  TString userName=gSystem->GetFromPipe("echo $USER");
  TString date=gSystem->GetFromPipe("date");

  if (!comment) metaData->SetComment(Form("Space point distortion calibration\n User: %s\n Data%s",userName.Data(),date.Data()));
  if (comment) metaData->SetComment(comment);
  AliCDBId* id1=NULL;
  id1=new AliCDBId("TPC/Calib/Correction", startRun, endRun);
  AliCDBStorage* gStorage = AliCDBManager::Instance()->GetStorage(ocdbStorage);
  gStorage->Put(this, (*id1), metaData);
}

