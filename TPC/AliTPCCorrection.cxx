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
#include "TVectorD.h"
#include "AliTPCParamSR.h"

#include "AliTPCCorrection.h"
#include "AliLog.h"

#include "AliExternalTrackParam.h"
#include "AliTrackPointArray.h"
#include "TDatabasePDG.h"
#include "AliTrackerBase.h"
#include "AliTPCROC.h"
#include "THnSparse.h"

#include "AliTPCLaserTrack.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "TDatabasePDG.h"
#include "TF1.h"
#include "TRandom.h"

#include "TDatabasePDG.h"

#include "AliTPCTransform.h"
#include "AliTPCcalibDB.h"
#include "AliTPCExB.h"

#include "AliTPCRecoParam.h"


ClassImp(AliTPCCorrection)


TObjArray *AliTPCCorrection::fgVisualCorrection=0;
// instance of correction for visualization


// FIXME: the following values should come from the database
const Double_t AliTPCCorrection::fgkTPCZ0    = 249.7;     // nominal gating grid position 
const Double_t AliTPCCorrection::fgkIFCRadius=  83.5;     // radius which renders the "18 rod manifold" best -> compare calc. of Jim Thomas
// compare gkIFCRadius=  83.05: Mean Radius of the Inner Field Cage ( 82.43 min,  83.70 max) (cm)
const Double_t AliTPCCorrection::fgkOFCRadius= 254.5;     // Mean Radius of the Outer Field Cage (252.55 min, 256.45 max) (cm)
const Double_t AliTPCCorrection::fgkZOffSet  =   0.2;     // Offset from CE: calculate all distortions closer to CE as if at this point
const Double_t AliTPCCorrection::fgkCathodeV = -100000.0; // Cathode Voltage (volts)
const Double_t AliTPCCorrection::fgkGG       =     -70.0; // Gating Grid voltage (volts)

const Double_t  AliTPCCorrection::fgkdvdE = 0.0024; // [cm/V] drift velocity dependency on the E field (from Magboltz for NeCO2N2 at standard environment)

const Double_t AliTPCCorrection::fgkEM = -1.602176487e-19/9.10938215e-31; // charge/mass in [C/kg]
const Double_t AliTPCCorrection::fgke0 = 8.854187817e-12;                 // vacuum permittivity [A·s/(V·m)]
 

AliTPCCorrection::AliTPCCorrection() 
  : TNamed("correction_unity","unity"),fILow(0),fJLow(0),fKLow(0), fT1(1), fT2(1)
{
  //
  // default constructor
  //
  if (!fgVisualCorrection) fgVisualCorrection= new TObjArray;

  InitLookUpfulcrums();

}

AliTPCCorrection::AliTPCCorrection(const char *name,const char *title)
: TNamed(name,title),fILow(0),fJLow(0),fKLow(0), fT1(1), fT2(1)
{
  //
  // default constructor, that set the name and title
  //
  if (!fgVisualCorrection) fgVisualCorrection= new TObjArray;

  InitLookUpfulcrums();

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
  AliTPCParam* tpcparam = new AliTPCParamSR;

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
      if (tpcparam->GetPadRowRadii(0,0)<=r0 && r0<=tpcparam->GetPadRowRadii(36,95)) {
	Float_t r1=TMath::Sqrt((x[0]+dx[0])*(x[0]+dx[0])+(x[1]+dx[1])*(x[1]+dx[1]));
	h->SetBinContent(ix,iy,r1-r0);
      }
      else
	h->SetBinContent(ix,iy,0.);
    }
  }
  delete tpcparam;
  return h;
}

TH2F* AliTPCCorrection::CreateHistoDRPhiinXY(Float_t z,Int_t nx,Int_t ny) {
  //
  // Simple plot functionality.
  // Returns a 2d hisogram which represents the corrections in rphi direction (drphi) 
  // in respect to position z within the XY plane.
  // The histogramm has nx times ny entries. 
  //

  AliTPCParam* tpcparam = new AliTPCParamSR;

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
      if (tpcparam->GetPadRowRadii(0,0)<=r0 && r0<=tpcparam->GetPadRowRadii(36,95)) {
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
  delete tpcparam;
  return h;
}

TH2F* AliTPCCorrection::CreateHistoDZinXY(Float_t z,Int_t nx,Int_t ny) {
  //
  // Simple plot functionality.
  // Returns a 2d hisogram which represents the corrections in longitudinal direction (dz)
  // in respect to position z within the XY plane.
  // The histogramm has nx times ny entries. 
  //

  AliTPCParam* tpcparam = new AliTPCParamSR;
 
  TH2F *h=CreateTH2F("dz_xy",GetTitle(),"x [cm]","y [cm]","dz [cm]",
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
      if (tpcparam->GetPadRowRadii(0,0)<=r0 && r0<=tpcparam->GetPadRowRadii(36,95)) {
	h->SetBinContent(ix,iy,dx[2]);
      }
      else
	h->SetBinContent(ix,iy,0.);
    }
  }
  delete tpcparam;
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

TH2F* AliTPCCorrection::CreateHistoDZinZR(Float_t phi,Int_t nz,Int_t nr) {
  //
  // Simple plot functionality.
  // Returns a 2d hisogram which represents the corrections in longitudinal direction (dz) 
  // in respect to angle phi within the ZR plane.
  // The histogramm has nx times ny entries. 
  //
  TH2F *h=CreateTH2F("dz_zr",GetTitle(),"z [cm]","r [cm]","dz [cm]",
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
      h->SetBinContent(iz,ir,dx[2]);
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
  Double_t saveEr[5] = {0,0,0,0,0};

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

void AliTPCCorrection::Interpolate3DEdistortion( const Int_t order, const Double_t r, const Float_t phi, const Double_t z, 
						 const Double_t er[kNZ][kNPhi][kNR], const Double_t ephi[kNZ][kNPhi][kNR], const Double_t ez[kNZ][kNPhi][kNR],
						 Double_t &erValue, Double_t &ephiValue, Double_t &ezValue) {
  //
  // Interpolate table - 3D interpolation
  //
  
  Double_t saveEr[5]= {0,0,0,0,0};
  Double_t savedEr[5]= {0,0,0,0,0} ;

  Double_t saveEphi[5]= {0,0,0,0,0};
  Double_t savedEphi[5]= {0,0,0,0,0} ;

  Double_t saveEz[5]= {0,0,0,0,0};
  Double_t savedEz[5]= {0,0,0,0,0} ;

  Search( kNZ,   fgkZList,   z,   fILow   ) ;
  Search( kNPhi, fgkPhiList, z,   fJLow   ) ;
  Search( kNR,   fgkRList,   r,   fKLow   ) ;

  if ( fILow < 0 ) fILow = 0 ;   // check if out of range
  if ( fJLow < 0 ) fJLow = 0 ;
  if ( fKLow < 0 ) fKLow = 0 ;

  if ( fILow + order  >=    kNZ - 1 ) fILow =   kNZ - 1 - order ;
  if ( fJLow + order  >=  kNPhi - 1 ) fJLow = kNPhi - 1 - order ;
  if ( fKLow + order  >=    kNR - 1 ) fKLow =   kNR - 1 - order ;

  for ( Int_t i = fILow ; i < fILow + order + 1 ; i++ ) {
    for ( Int_t j = fJLow ; j < fJLow + order + 1 ; j++ ) {
      saveEr[j-fJLow]     = Interpolate( &fgkRList[fKLow], &er[i][j][fKLow], order, r )   ;
      saveEphi[j-fJLow]   = Interpolate( &fgkRList[fKLow], &ephi[i][j][fKLow], order, r ) ;
      saveEz[j-fJLow]     = Interpolate( &fgkRList[fKLow], &ez[i][j][fKLow], order, r )   ;
    }
    savedEr[i-fILow]     = Interpolate( &fgkPhiList[fJLow], saveEr, order, phi )   ; 
    savedEphi[i-fILow]   = Interpolate( &fgkPhiList[fJLow], saveEphi, order, phi ) ; 
    savedEz[i-fILow]     = Interpolate( &fgkPhiList[fJLow], saveEz, order, phi )   ; 
  }
  erValue     = Interpolate( &fgkZList[fILow], savedEr, order, z )    ;
  ephiValue   = Interpolate( &fgkZList[fILow], savedEphi, order, z )  ;
  ezValue     = Interpolate( &fgkZList[fILow], savedEz, order, z )    ;

}

Double_t AliTPCCorrection::Interpolate2DTable( const Int_t order, const Double_t x, const Double_t y, 
					      const Int_t nx,  const Int_t ny, const Double_t xv[], const Double_t yv[], 
					      const TMatrixD &array ) {
  //
  // Interpolate table (TMatrix format) - 2D interpolation
  //

  static  Int_t jlow = 0, klow = 0 ;
  Double_t saveArray[5] = {0,0,0,0,0} ;

  Search( nx,  xv,  x,   jlow  ) ;
  Search( ny,  yv,  y,   klow  ) ;
  if ( jlow < 0 ) jlow = 0 ;   // check if out of range
  if ( klow < 0 ) klow = 0 ;
  if ( jlow + order  >=    nx - 1 ) jlow =   nx - 1 - order ;
  if ( klow + order  >=    ny - 1 ) klow =   ny - 1 - order ;

  for ( Int_t j = jlow ; j < jlow + order + 1 ; j++ )
    {
      Double_t *ajkl = &((TMatrixD&)array)(j,klow);
      saveArray[j-jlow]  = Interpolate( &yv[klow], ajkl , order, y )   ;
    }

  return( Interpolate( &xv[jlow], saveArray, order, x ) )   ;

}

Double_t AliTPCCorrection::Interpolate3DTable( const Int_t order, const Double_t x,   const Double_t y,   const Double_t z,
					      const Int_t  nx,    const Int_t  ny,    const Int_t  nz,
					      const Double_t xv[], const Double_t yv[], const Double_t zv[],
					      TMatrixD **arrayofArrays ) {
  //
  // Interpolate table (TMatrix format) - 3D interpolation
  //

  static  Int_t ilow = 0, jlow = 0, klow = 0 ;
  Double_t saveArray[5]= {0,0,0,0,0};
  Double_t savedArray[5]= {0,0,0,0,0} ;

  Search( nx, xv, x, ilow   ) ;
  Search( ny, yv, y, jlow   ) ;
  Search( nz, zv, z, klow   ) ;  

  if ( ilow < 0 ) ilow = 0 ;   // check if out of range
  if ( jlow < 0 ) jlow = 0 ;
  if ( klow < 0 ) klow = 0 ;

  if ( ilow + order  >=    nx - 1 ) ilow =   nx - 1 - order ;
  if ( jlow + order  >=    ny - 1 ) jlow =   ny - 1 - order ;
  if ( klow + order  >=    nz - 1 ) klow =   nz - 1 - order ;

  for ( Int_t k = klow ; k < klow + order + 1 ; k++ )
    {
      TMatrixD &table = *arrayofArrays[k] ;
      for ( Int_t i = ilow ; i < ilow + order + 1 ; i++ )
	{
	  saveArray[i-ilow] = Interpolate( &yv[jlow], &table(i,jlow), order, y )   ;
	}
      savedArray[k-klow] = Interpolate( &xv[ilow], saveArray, order, x )   ; 
    }
  return( Interpolate( &zv[klow], savedArray, order, z ) )   ;

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

void AliTPCCorrection::InitLookUpfulcrums() {
  //
  // Initialization of interpolation points - for main look up table
  //   (course grid in the middle, fine grid on the borders)
  //

  AliTPCROC * roc = AliTPCROC::Instance();
  const Double_t rLow =  TMath::Floor(roc->GetPadRowRadii(0,0))-1; // first padRow plus some margin 

  // fulcrums in R
  fgkRList[0] = rLow;
  for (Int_t i = 1; i<kNR; i++) {
    fgkRList[i] = fgkRList[i-1] + 3.5;     // 3.5 cm spacing    
    if (fgkRList[i]<90 ||fgkRList[i]>245) 
       fgkRList[i] = fgkRList[i-1] + 0.5; // 0.5 cm spacing
    else if (fgkRList[i]<100 || fgkRList[i]>235) 
       fgkRList[i] = fgkRList[i-1] + 1.5;  // 1.5 cm spacing
    else if (fgkRList[i]<120 || fgkRList[i]>225) 
       fgkRList[i] = fgkRList[i-1] + 2.5;  // 2.5 cm spacing
  }

  // fulcrums in Z
  fgkZList[0] = -249.5;
  fgkZList[kNZ-1] = 249.5;
  for (Int_t j = 1; j<kNZ/2; j++) {
    fgkZList[j] = fgkZList[j-1];
    if      (TMath::Abs(fgkZList[j])< 0.15)
      fgkZList[j] = fgkZList[j-1] + 0.09; // 0.09 cm spacing
    else if(TMath::Abs(fgkZList[j])< 0.6)
      fgkZList[j] = fgkZList[j-1] + 0.4; // 0.4 cm spacing
    else if      (TMath::Abs(fgkZList[j])< 2.5 || TMath::Abs(fgkZList[j])>248) 
      fgkZList[j] = fgkZList[j-1] + 0.5; // 0.5 cm spacing
    else if (TMath::Abs(fgkZList[j])<10 || TMath::Abs(fgkZList[j])>235) 
      fgkZList[j] = fgkZList[j-1] + 1.5;  // 1.5 cm spacing
    else if (TMath::Abs(fgkZList[j])<25 || TMath::Abs(fgkZList[j])>225) 
      fgkZList[j] = fgkZList[j-1] + 2.5;  // 2.5 cm spacing
    else 
      fgkZList[j] = fgkZList[j-1] + 4;  // 4 cm spacing

    fgkZList[kNZ-j-1] = -fgkZList[j];
  }
  
  // fulcrums in phi
  for (Int_t k = 0; k<kNPhi; k++) 
    fgkPhiList[k] = TMath::TwoPi()*k/(kNPhi-1);    
  
  
}


void AliTPCCorrection::PoissonRelaxation2D(TMatrixD &arrayV, TMatrixD &chargeDensity, 
					   TMatrixD &arrayErOverEz, TMatrixD &arrayDeltaEz, 
					   const Int_t rows, const Int_t columns, const Int_t iterations,
					   const Bool_t rocDisplacement ) {
  //
  // Solve Poisson's Equation by Relaxation Technique in 2D (assuming cylindrical symmetry)
  //
  // Solve Poissons equation in a cylindrical coordinate system. The arrayV matrix must be filled with the 
  // boundary conditions on the first and last rows, and the first and last columns.  The remainder of the 
  // array can be blank or contain a preliminary guess at the solution.  The Charge density matrix contains 
  // the enclosed spacecharge density at each point. The charge density matrix can be full of zero's if 
  // you wish to solve Laplaces equation however it should not contain random numbers or you will get 
  // random numbers back as a solution. 
  // Poisson's equation is solved by iteratively relaxing the matrix to the final solution.  In order to 
  // speed up the convergence to the best solution, this algorithm does a binary expansion of the solution 
  // space.  First it solves the problem on a very sparse grid by skipping rows and columns in the original 
  // matrix.  Then it doubles the number of points and solves the problem again.  Then it doubles the 
  // number of points and solves the problem again.  This happens several times until the maximum number
  // of points has been included in the array.  
  //
  // NOTE: In order for this algorithmto work, the number of rows and columns must be a power of 2 plus one.
  // So rows == 2**M + 1 and columns == 2**N + 1.  The number of rows and columns can be different.
  // 
  // NOTE: rocDisplacement is used to include (or ignore) the ROC misalignment in the dz calculation
  //
  // Original code by Jim Thomas (STAR TPC Collaboration)
  //

  Double_t ezField = (fgkCathodeV-fgkGG)/fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm; 

  const Float_t  gridSizeR   =  (fgkOFCRadius-fgkIFCRadius) / (rows-1) ;
  const Float_t  gridSizeZ   =  fgkTPCZ0 / (columns-1) ;
  const Float_t  ratio       =  gridSizeR*gridSizeR / (gridSizeZ*gridSizeZ) ;

  TMatrixD  arrayEr(rows,columns) ;
  TMatrixD  arrayEz(rows,columns) ;

  //Check that number of rows and columns is suitable for a binary expansion
  
  if ( !IsPowerOfTwo(rows-1) ) {
    AliError("PoissonRelaxation - Error in the number of rows. Must be 2**M - 1");
    return;
  }
  if ( !IsPowerOfTwo(columns-1) ) {
    AliError("PoissonRelaxation - Error in the number of columns. Must be 2**N - 1");
    return;
  }
  
  // Solve Poisson's equation in cylindrical coordinates by relaxation technique
  // Allow for different size grid spacing in R and Z directions
  // Use a binary expansion of the size of the matrix to speed up the solution of the problem
  
  Int_t iOne = (rows-1)/4 ;
  Int_t jOne = (columns-1)/4 ;
  // Solve for N in 2**N, add one.
  Int_t loops = 1 + (int) ( 0.5 + TMath::Log2( (double) TMath::Max(iOne,jOne) ) ) ;  

  for ( Int_t count = 0 ; count < loops ; count++ ) { 
    // Loop while the matrix expands & the resolution increases.

    Float_t tempGridSizeR = gridSizeR * iOne ;
    Float_t tempRatio     = ratio * iOne * iOne / ( jOne * jOne ) ;
    Float_t tempFourth    = 1.0 / (2.0 + 2.0*tempRatio) ;
    
    // Do this the standard C++ way to avoid gcc extensions for Float_t coef1[rows]
    std::vector<float> coef1(rows) ;  
    std::vector<float> coef2(rows) ;  

    for ( Int_t i = iOne ; i < rows-1 ; i+=iOne ) {
       Float_t radius = fgkIFCRadius + i*gridSizeR ;
      coef1[i] = 1.0 + tempGridSizeR/(2*radius);
      coef2[i] = 1.0 - tempGridSizeR/(2*radius);
    }
    
    TMatrixD sumChargeDensity(rows,columns) ;

    for ( Int_t i = iOne ; i < rows-1 ; i += iOne ) {
      Float_t radius = fgkIFCRadius + iOne*gridSizeR ;
      for ( Int_t j = jOne ; j < columns-1 ; j += jOne ) {
	if ( iOne == 1 && jOne == 1 ) sumChargeDensity(i,j) = chargeDensity(i,j) ;
	else {        
	  // Add up all enclosed charge density contributions within 1/2 unit in all directions
	  Float_t weight = 0.0 ;
	  Float_t sum    = 0.0 ;
	  sumChargeDensity(i,j) = 0.0 ;
	  for ( Int_t ii = i-iOne/2 ; ii <= i+iOne/2 ; ii++ ) {
	    for ( Int_t jj = j-jOne/2 ; jj <= j+jOne/2 ; jj++ ) {
	      if ( ii == i-iOne/2 || ii == i+iOne/2 || jj == j-jOne/2 || jj == j+jOne/2 ) weight = 0.5 ;
	      else
		weight = 1.0 ;
	      // Note that this is cylindrical geometry
	      sumChargeDensity(i,j) += chargeDensity(ii,jj)*weight*radius ;  
	      sum += weight*radius ;
	    }
	  }
	  sumChargeDensity(i,j) /= sum ;
	}
        sumChargeDensity(i,j) *= tempGridSizeR*tempGridSizeR; // just saving a step later on
       }
    }

    for ( Int_t k = 1 ; k <= iterations; k++ ) {               
      // Solve Poisson's Equation
      // Over-relaxation index, must be >= 1 but < 2.  Arrange for it to evolve from 2 => 1 
      // as interations increase.
      Float_t overRelax   = 1.0 + TMath::Sqrt( TMath::Cos( (k*TMath::PiOver2())/iterations ) ) ; 
      Float_t overRelaxM1 = overRelax - 1.0 ;
      Float_t overRelaxtempFourth, overRelaxcoef5 ;
      overRelaxtempFourth = overRelax * tempFourth ;
      overRelaxcoef5 = overRelaxM1 / overRelaxtempFourth ; 

      for ( Int_t i = iOne ; i < rows-1 ; i += iOne ) {
	for ( Int_t j = jOne ; j < columns-1 ; j += jOne ) {

	  arrayV(i,j) = (   coef2[i]       *   arrayV(i-iOne,j)
			  + tempRatio      * ( arrayV(i,j-jOne) + arrayV(i,j+jOne) )
			  - overRelaxcoef5 *   arrayV(i,j) 
			  + coef1[i]       *   arrayV(i+iOne,j) 
			  + sumChargeDensity(i,j) 
			) * overRelaxtempFourth;
	}
      }

      if ( k == iterations ) {    
	// After full solution is achieved, copy low resolution solution into higher res array
	for ( Int_t i = iOne ; i < rows-1 ; i += iOne ) {
	  for ( Int_t j = jOne ; j < columns-1 ; j += jOne ) {

	    if ( iOne > 1 ) {              
	      arrayV(i+iOne/2,j)                    =  ( arrayV(i+iOne,j) + arrayV(i,j)     ) / 2 ;
	      if ( i == iOne )  arrayV(i-iOne/2,j) =  ( arrayV(0,j)       + arrayV(iOne,j) ) / 2 ;
	    }
	    if ( jOne > 1 ) {
	      arrayV(i,j+jOne/2)                    =  ( arrayV(i,j+jOne) + arrayV(i,j) )     / 2 ;
	      if ( j == jOne )  arrayV(i,j-jOne/2) =  ( arrayV(i,0)       + arrayV(i,jOne) ) / 2 ;
	    }
	    if ( iOne > 1 && jOne > 1 ) {
	      arrayV(i+iOne/2,j+jOne/2) =  ( arrayV(i+iOne,j+jOne) + arrayV(i,j) ) / 2 ;
	      if ( i == iOne ) arrayV(i-iOne/2,j-jOne/2) =   ( arrayV(0,j-jOne) + arrayV(iOne,j) ) / 2 ;
	      if ( j == jOne ) arrayV(i-iOne/2,j-jOne/2) =   ( arrayV(i-iOne,0) + arrayV(i,jOne) ) / 2 ;
	      // Note that this leaves a point at the upper left and lower right corners uninitialized. 
	      // -> Not a big deal.
	    }

	  }
	}
      }

    }

    iOne = iOne / 2 ; if ( iOne < 1 ) iOne = 1 ;
    jOne = jOne / 2 ; if ( jOne < 1 ) jOne = 1 ;

    sumChargeDensity.Clear();
  }      

  // Differentiate V(r) and solve for E(r) using special equations for the first and last rows
  for ( Int_t j = 0 ; j < columns ; j++ ) {	  
    for ( Int_t i = 1 ; i < rows-1 ; i++ ) arrayEr(i,j) = -1 * ( arrayV(i+1,j) - arrayV(i-1,j) ) / (2*gridSizeR) ;
    arrayEr(0,j)      =  -1 * ( -0.5*arrayV(2,j) + 2.0*arrayV(1,j) - 1.5*arrayV(0,j) ) / gridSizeR ;  
    arrayEr(rows-1,j) =  -1 * ( 1.5*arrayV(rows-1,j) - 2.0*arrayV(rows-2,j) + 0.5*arrayV(rows-3,j) ) / gridSizeR ; 
  }

  // Differentiate V(z) and solve for E(z) using special equations for the first and last columns
  for ( Int_t i = 0 ; i < rows ; i++) {
    for ( Int_t j = 1 ; j < columns-1 ; j++ ) arrayEz(i,j) = -1 * ( arrayV(i,j+1) - arrayV(i,j-1) ) / (2*gridSizeZ) ;
    arrayEz(i,0)         =  -1 * ( -0.5*arrayV(i,2) + 2.0*arrayV(i,1) - 1.5*arrayV(i,0) ) / gridSizeZ ;  
    arrayEz(i,columns-1) =  -1 * ( 1.5*arrayV(i,columns-1) - 2.0*arrayV(i,columns-2) + 0.5*arrayV(i,columns-3) ) / gridSizeZ ; 
  }
  
  for ( Int_t i = 0 ; i < rows ; i++) {
    // Note: go back and compare to old version of this code.  See notes below.
    // JT Test ... attempt to divide by real Ez not Ez to first order
    for ( Int_t j = 0 ; j < columns ; j++ ) {
      arrayEz(i,j) += ezField;
      // This adds back the overall Z gradient of the field (main E field component)
    } 
    // Warning: (-=) assumes you are using an error potetial without the overall Field included
  }                                 
  
  // Integrate Er/Ez from Z to zero
  for ( Int_t j = 0 ; j < columns ; j++ )  {	  
    for ( Int_t i = 0 ; i < rows ; i++ ) {
      
      Int_t index = 1 ;   // Simpsons rule if N=odd.  If N!=odd then add extra point by trapezoidal rule.  
      arrayErOverEz(i,j) = 0.0 ;
      arrayDeltaEz(i,j) = 0.0 ;
      
      for ( Int_t k = j ; k < columns ; k++ ) {
	arrayErOverEz(i,j)  +=  index*(gridSizeZ/3.0)*arrayEr(i,k)/arrayEz(i,k) ;
	arrayDeltaEz(i,j)   +=  index*(gridSizeZ/3.0)*(arrayEz(i,k)-ezField) ;
	if ( index != 4 )  index = 4; else index = 2 ;
      }
      if ( index == 4 ) {
	arrayErOverEz(i,j)  -=  (gridSizeZ/3.0)*arrayEr(i,columns-1)/arrayEz(i,columns-1) ;
	arrayDeltaEz(i,j)   -=  (gridSizeZ/3.0)*(arrayEz(i,columns-1)-ezField) ;
      }
      if ( index == 2 ) {
	arrayErOverEz(i,j)  +=  (gridSizeZ/3.0) * ( 0.5*arrayEr(i,columns-2)/arrayEz(i,columns-2) 
						    -2.5*arrayEr(i,columns-1)/arrayEz(i,columns-1));
	arrayDeltaEz(i,j)   +=  (gridSizeZ/3.0) * ( 0.5*(arrayEz(i,columns-2)-ezField) 
						    -2.5*(arrayEz(i,columns-1)-ezField));
      }
      if ( j == columns-2 ) {
	arrayErOverEz(i,j) =  (gridSizeZ/3.0) * ( 1.5*arrayEr(i,columns-2)/arrayEz(i,columns-2)
						  +1.5*arrayEr(i,columns-1)/arrayEz(i,columns-1) ) ;
	arrayDeltaEz(i,j)  =  (gridSizeZ/3.0) * ( 1.5*(arrayEz(i,columns-2)-ezField)
						  +1.5*(arrayEz(i,columns-1)-ezField) ) ;
      }
      if ( j == columns-1 ) {
	arrayErOverEz(i,j) =  0.0 ;
	arrayDeltaEz(i,j)  =  0.0 ;
      }
    }
  }
  
  // calculate z distortion from the integrated Delta Ez residuals
  // and include the aquivalence (Volt to cm) of the ROC shift !!

  for ( Int_t j = 0 ; j < columns ; j++ )  {	  
    for ( Int_t i = 0 ; i < rows ; i++ ) {

      // Scale the Ez distortions with the drift velocity pertubation -> delivers cm
      arrayDeltaEz(i,j) = arrayDeltaEz(i,j)*fgkdvdE;

      // ROC Potential in cm aquivalent
      Double_t dzROCShift =  arrayV(i, columns -1)/ezField;  
      if ( rocDisplacement ) arrayDeltaEz(i,j) = arrayDeltaEz(i,j) + dzROCShift;  // add the ROC misaligment

    }
  }
 
  arrayEr.Clear();
  arrayEz.Clear();

}

void AliTPCCorrection::PoissonRelaxation3D( TMatrixD**arrayofArrayV, TMatrixD**arrayofChargeDensities, 
		    TMatrixD**arrayofEroverEz, TMatrixD**arrayofEPhioverEz, TMatrixD**arrayofDeltaEz,
		    const Int_t rows, const Int_t columns,  const Int_t phislices, 
		    const Float_t deltaphi, const Int_t iterations, const Int_t symmetry,
		    Bool_t rocDisplacement  ) {
  //
  // 3D - Solve Poisson's Equation in 3D by Relaxation Technique
  //
  //    NOTE: In order for this algorith to work, the number of rows and columns must be a power of 2 plus one.  
  //    The number of rows and COLUMNS can be different.
  //
  //    ROWS       ==  2**M + 1  
  //    COLUMNS    ==  2**N + 1  
  //    PHISLICES  ==  Arbitrary but greater than 3
  //
  //    DeltaPhi in Radians
  //
  //    SYMMETRY = 0 if no phi symmetries, and no phi boundary conditions
  //             = 1 if we have reflection symmetry at the boundaries (eg. sector symmetry or half sector symmetries).
  //
  // NOTE: rocDisplacement is used to include (or ignore) the ROC misalignment in the dz calculation

  const Double_t ezField = (fgkCathodeV-fgkGG)/fgkTPCZ0; // = ALICE Electric Field (V/cm) Magnitude ~ -400 V/cm; 

  const Float_t  gridSizeR   =  (fgkOFCRadius-fgkIFCRadius) / (rows-1) ;
  const Float_t  gridSizePhi =  deltaphi ;
  const Float_t  gridSizeZ   =  fgkTPCZ0 / (columns-1) ;
  const Float_t  ratioPhi    =  gridSizeR*gridSizeR / (gridSizePhi*gridSizePhi) ;
  const Float_t  ratioZ      =  gridSizeR*gridSizeR / (gridSizeZ*gridSizeZ) ;

  TMatrixD arrayE(rows,columns) ;

  // Check that the number of rows and columns is suitable for a binary expansion
  if ( !IsPowerOfTwo((rows-1))    ) {  
    AliError("Poisson3DRelaxation - Error in the number of rows. Must be 2**M - 1"); 
    return; }
  if ( !IsPowerOfTwo((columns-1)) ) { 
    AliError("Poisson3DRelaxation - Error in the number of columns. Must be 2**N - 1");
    return; }
  if ( phislices <= 3   )  { 
    AliError("Poisson3DRelaxation - Error in the number of phislices. Must be larger than 3");
    return; }
  if  ( phislices > 1000 ) { 
    AliError("Poisson3D  phislices > 1000 is not allowed (nor wise) ");  
    return; }  
  
  // Solve Poisson's equation in cylindrical coordinates by relaxation technique
  // Allow for different size grid spacing in R and Z directions
  // Use a binary expansion of the matrix to speed up the solution of the problem

  Int_t loops, mplus, mminus, signplus, signminus  ;
  Int_t ione = (rows-1)/4 ;
  Int_t jone = (columns-1)/4 ;
  loops = TMath::Max(ione, jone) ;      // Calculate the number of loops for the binary expansion
  loops = 1 + (int) ( 0.5 + TMath::Log2((double)loops) ) ;  // Solve for N in 2**N

  TMatrixD* arrayofSumChargeDensities[1000] ;    // Create temporary arrays to store low resolution charge arrays

  for ( Int_t i = 0 ; i < phislices ; i++ ) { arrayofSumChargeDensities[i] = new TMatrixD(rows,columns) ; }

  for ( Int_t count = 0 ; count < loops ; count++ ) {      // START the master loop and do the binary expansion
   
    Float_t  tempgridSizeR   =  gridSizeR  * ione ;
    Float_t  tempratioPhi    =  ratioPhi * ione * ione ; // Used tobe divided by ( m_one * m_one ) when m_one was != 1
    Float_t  tempratioZ      =  ratioZ   * ione * ione / ( jone * jone ) ;

    std::vector<float> coef1(rows) ;  // Do this the standard C++ way to avoid gcc extensions for Float_t coef1[rows]
    std::vector<float> coef2(rows) ;  // Do this the standard C++ way to avoid gcc extensions for Float_t coef1[rows]
    std::vector<float> coef3(rows) ;  // Do this the standard C++ way to avoid gcc extensions for Float_t coef1[rows]
    std::vector<float> coef4(rows) ;  // Do this the standard C++ way to avoid gcc extensions for Float_t coef1[rows]

    for ( Int_t i = ione ; i < rows-1 ; i+=ione )  {
      Float_t radius = fgkIFCRadius + i*gridSizeR ;
      coef1[i] = 1.0 + tempgridSizeR/(2*radius);
      coef2[i] = 1.0 - tempgridSizeR/(2*radius);
      coef3[i] = tempratioPhi/(radius*radius);
      coef4[i] = 0.5 / (1.0 + tempratioZ + coef3[i]);
    }

    for ( Int_t m = 0 ; m < phislices ; m++ ) {
      TMatrixD &chargeDensity    = *arrayofChargeDensities[m] ;
      TMatrixD &sumChargeDensity = *arrayofSumChargeDensities[m] ;
      for ( Int_t i = ione ; i < rows-1 ; i += ione ) {
	Float_t radius = fgkIFCRadius + i*gridSizeR ;
	for ( Int_t j = jone ; j < columns-1 ; j += jone ) {
	  if ( ione == 1 && jone == 1 ) sumChargeDensity(i,j) = chargeDensity(i,j) ;
	  else {           // Add up all enclosed charge density contributions within 1/2 unit in all directions
	    Float_t weight = 0.0 ;
	    Float_t sum    = 0.0 ;
	    sumChargeDensity(i,j) = 0.0 ;
	    for ( Int_t ii = i-ione/2 ; ii <= i+ione/2 ; ii++ ) {
	      for ( Int_t jj = j-jone/2 ; jj <= j+jone/2 ; jj++ ) {
		if ( ii == i-ione/2 || ii == i+ione/2 || jj == j-jone/2 || jj == j+jone/2 ) weight = 0.5 ;
		else
		  weight = 1.0 ; 
		sumChargeDensity(i,j) += chargeDensity(ii,jj)*weight*radius ;  
		sum += weight*radius ;
	      }
	    }
	    sumChargeDensity(i,j) /= sum ;
	  }
          sumChargeDensity(i,j) *= tempgridSizeR*tempgridSizeR; // just saving a step later on
	}
      }
    }

    for ( Int_t k = 1 ; k <= iterations; k++ ) {

      // over-relaxation index, >= 1 but < 2
      Float_t overRelax   = 1.0 + TMath::Sqrt( TMath::Cos( (k*TMath::PiOver2())/iterations ) ) ; 
      Float_t overRelaxM1 = overRelax - 1.0 ;

      std::vector<float> overRelaxcoef4(rows) ;  // Do this the standard C++ way to avoid gcc extensions
      std::vector<float> overRelaxcoef5(rows) ;  // Do this the standard C++ way to avoid gcc extensions

      for ( Int_t i = ione ; i < rows-1 ; i+=ione ) { 
	overRelaxcoef4[i] = overRelax * coef4[i] ;
	overRelaxcoef5[i] = overRelaxM1 / overRelaxcoef4[i] ; 
      }

      for ( Int_t m = 0 ; m < phislices ; m++ ) {

	mplus  = m + 1;   signplus  = 1 ; 
	mminus = m - 1 ;  signminus = 1 ;
	if (symmetry==1) {  // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
	  if ( mplus  > phislices-1 ) mplus  = phislices - 2 ;
	  if ( mminus < 0 )           mminus = 1 ;
	}
	else if (symmetry==-1) {   // Anti-symmetry in phi
	  if ( mplus  > phislices-1 ) { mplus  = phislices - 2 ; signplus  = -1 ; }
	  if ( mminus < 0 )           { mminus = 1 ;	         signminus = -1 ; } 
	}
		else { // No Symmetries in phi, no boundaries, the calculation is continuous across all phi
	  if ( mplus  > phislices-1 ) mplus  = m + 1 - phislices ;
	  if ( mminus < 0 )           mminus = m - 1 + phislices ;
	}
	TMatrixD& arrayV    =  *arrayofArrayV[m] ;
	TMatrixD& arrayVP   =  *arrayofArrayV[mplus] ;
	TMatrixD& arrayVM   =  *arrayofArrayV[mminus] ;
	TMatrixD& sumChargeDensity =  *arrayofSumChargeDensities[m] ;

	for ( Int_t i = ione ; i < rows-1 ; i+=ione )  {
	  for ( Int_t j = jone ; j < columns-1 ; j+=jone ) {

            arrayV(i,j) = (   coef2[i]          *   arrayV(i-ione,j)
			    + tempratioZ        * ( arrayV(i,j-jone)  +  arrayV(i,j+jone) )
			    - overRelaxcoef5[i] *   arrayV(i,j) 
			    + coef1[i]          *   arrayV(i+ione,j)  
			    + coef3[i]          * ( signplus*arrayVP(i,j)       +  signminus*arrayVM(i,j) )
			    + sumChargeDensity(i,j) 
			  ) * overRelaxcoef4[i] ;     
	    // Note: over-relax the solution at each step.  This speeds up the convergance.

	  }
	}

	if ( k == iterations ) {   // After full solution is achieved, copy low resolution solution into higher res array
	  for ( Int_t i = ione ; i < rows-1 ; i+=ione )  {
	    for ( Int_t j = jone ; j < columns-1 ; j+=jone ) {
	      
	      if ( ione > 1 ) {              
		arrayV(i+ione/2,j)                    =  ( arrayV(i+ione,j) + arrayV(i,j)     ) / 2 ;
		if ( i == ione )  arrayV(i-ione/2,j) =  ( arrayV(0,j)       + arrayV(ione,j) ) / 2 ;
	      }
	      if ( jone > 1 ) {
		arrayV(i,j+jone/2)                    =  ( arrayV(i,j+jone) + arrayV(i,j) )     / 2 ;
		if ( j == jone )  arrayV(i,j-jone/2) =  ( arrayV(i,0)       + arrayV(i,jone) ) / 2 ;
	      }
	      if ( ione > 1 && jone > 1 ) {
		arrayV(i+ione/2,j+jone/2) =  ( arrayV(i+ione,j+jone) + arrayV(i,j) ) / 2 ;
		if ( i == ione ) arrayV(i-ione/2,j-jone/2) =   ( arrayV(0,j-jone) + arrayV(ione,j) ) / 2 ;
		if ( j == jone ) arrayV(i-ione/2,j-jone/2) =   ( arrayV(i-ione,0) + arrayV(i,jone) ) / 2 ;
		// Note that this leaves a point at the upper left and lower right corners uninitialized. Not a big deal.
	      }
	    }	    
	  }
	}

      }
    }      

    ione = ione / 2 ; if ( ione < 1 ) ione = 1 ;
    jone = jone / 2 ; if ( jone < 1 ) jone = 1 ;

  }
  
  //Differentiate V(r) and solve for E(r) using special equations for the first and last row
  //Integrate E(r)/E(z) from point of origin to pad plane

  for ( Int_t m = 0 ; m < phislices ; m++ ) {
    TMatrixD& arrayV    =  *arrayofArrayV[m] ;
    TMatrixD& eroverEz  =  *arrayofEroverEz[m] ;
    
    for ( Int_t j = columns-1 ; j >= 0 ; j-- ) {  // Count backwards to facilitate integration over Z
      
      // Differentiate in R
      for ( Int_t i = 1 ; i < rows-1 ; i++ )  arrayE(i,j) = -1 * ( arrayV(i+1,j) - arrayV(i-1,j) ) / (2*gridSizeR) ;
      arrayE(0,j)      =  -1 * ( -0.5*arrayV(2,j) + 2.0*arrayV(1,j) - 1.5*arrayV(0,j) ) / gridSizeR ;  
      arrayE(rows-1,j) =  -1 * ( 1.5*arrayV(rows-1,j) - 2.0*arrayV(rows-2,j) + 0.5*arrayV(rows-3,j) ) / gridSizeR ; 
      // Integrate over Z
      for ( Int_t i = 0 ; i < rows ; i++ ) {
	Int_t index = 1 ;   // Simpsons rule if N=odd.  If N!=odd then add extra point by trapezoidal rule.  
	eroverEz(i,j) = 0.0 ;
	for ( Int_t k = j ; k < columns ; k++ ) {
	  
	  eroverEz(i,j)  +=  index*(gridSizeZ/3.0)*arrayE(i,k)/(-1*ezField) ;
	  if ( index != 4 )  index = 4; else index = 2 ;
	}
	if ( index == 4 ) eroverEz(i,j)  -=  (gridSizeZ/3.0)*arrayE(i,columns-1)/ (-1*ezField) ;
	if ( index == 2 ) eroverEz(i,j)  +=  
	  (gridSizeZ/3.0)*(0.5*arrayE(i,columns-2)-2.5*arrayE(i,columns-1))/(-1*ezField) ;
	if ( j == columns-2 ) eroverEz(i,j) =  
	  (gridSizeZ/3.0)*(1.5*arrayE(i,columns-2)+1.5*arrayE(i,columns-1))/(-1*ezField) ;
	if ( j == columns-1 ) eroverEz(i,j) =  0.0 ;
      }
    }
    // if ( m == 0 ) { TCanvas*  c1 =  new TCanvas("erOverEz","erOverEz",50,50,840,600) ;  c1 -> cd() ;
    // eroverEz.Draw("surf") ; } // JT test
  }
  
  //Differentiate V(r) and solve for E(phi) 
  //Integrate E(phi)/E(z) from point of origin to pad plane

  for ( Int_t m = 0 ; m < phislices ; m++ ) {
    
    mplus  = m + 1;   signplus  = 1 ; 
    mminus = m - 1 ;  signminus = 1 ; 
    if (symmetry==1) { // Reflection symmetry in phi (e.g. symmetry at sector boundaries, or half sectors, etc.)
      if ( mplus  > phislices-1 ) mplus  = phislices - 2 ;
      if ( mminus < 0 )           mminus = 1 ;
    }
    else if (symmetry==-1) {       // Anti-symmetry in phi
      if ( mplus  > phislices-1 ) { mplus  = phislices - 2 ;  signplus  = -1 ; }
      if ( mminus < 0 )           { mminus = 1 ;	            signminus = -1 ; } 
    }
    else { // No Symmetries in phi, no boundaries, the calculations is continuous across all phi
      if ( mplus  > phislices-1 ) mplus  = m + 1 - phislices ;
      if ( mminus < 0 )           mminus = m - 1 + phislices ;
    }
    TMatrixD &arrayVP     =  *arrayofArrayV[mplus] ;
    TMatrixD &arrayVM     =  *arrayofArrayV[mminus] ;
    TMatrixD &ePhioverEz  =  *arrayofEPhioverEz[m] ;
    for ( Int_t j = columns-1 ; j >= 0 ; j-- ) { // Count backwards to facilitate integration over Z
      // Differentiate in Phi
      for ( Int_t i = 0 ; i < rows ; i++ ) {
	Float_t radius = fgkIFCRadius + i*gridSizeR ;
	arrayE(i,j) = -1 * (signplus * arrayVP(i,j) - signminus * arrayVM(i,j) ) / (2*radius*gridSizePhi) ;
      }
      // Integrate over Z
      for ( Int_t i = 0 ; i < rows ; i++ ) {
	Int_t index = 1 ;   // Simpsons rule if N=odd.  If N!=odd then add extra point by trapezoidal rule.  
	ePhioverEz(i,j) = 0.0 ;
	for ( Int_t k = j ; k < columns ; k++ ) {
	  
	  ePhioverEz(i,j)  +=  index*(gridSizeZ/3.0)*arrayE(i,k)/(-1*ezField) ;
	  if ( index != 4 )  index = 4; else index = 2 ;
	}
	if ( index == 4 ) ePhioverEz(i,j)  -=  (gridSizeZ/3.0)*arrayE(i,columns-1)/ (-1*ezField) ;
	if ( index == 2 ) ePhioverEz(i,j)  +=  
	  (gridSizeZ/3.0)*(0.5*arrayE(i,columns-2)-2.5*arrayE(i,columns-1))/(-1*ezField) ;
	if ( j == columns-2 ) ePhioverEz(i,j) =  
	  (gridSizeZ/3.0)*(1.5*arrayE(i,columns-2)+1.5*arrayE(i,columns-1))/(-1*ezField) ;
	if ( j == columns-1 ) ePhioverEz(i,j) =  0.0 ;
      }
    }
    // if ( m == 5 ) { TCanvas* c2 =  new TCanvas("arrayE","arrayE",50,50,840,600) ;  c2 -> cd() ;
    // arrayE.Draw("surf") ; } // JT test
  }
  

  // Differentiate V(r) and solve for E(z) using special equations for the first and last row
  // Integrate (E(z)-Ezstd) from point of origin to pad plane
  
  for ( Int_t m = 0 ; m < phislices ; m++ ) {
    TMatrixD& arrayV   =  *arrayofArrayV[m] ;
    TMatrixD& deltaEz  =  *arrayofDeltaEz[m] ;
    
    // Differentiate V(z) and solve for E(z) using special equations for the first and last columns
    for ( Int_t i = 0 ; i < rows ; i++) {
      for ( Int_t j = 1 ; j < columns-1 ; j++ ) arrayE(i,j) = -1 * ( arrayV(i,j+1) - arrayV(i,j-1) ) / (2*gridSizeZ) ;
      arrayE(i,0)         =  -1 * ( -0.5*arrayV(i,2) + 2.0*arrayV(i,1) - 1.5*arrayV(i,0) ) / gridSizeZ ;  
      arrayE(i,columns-1) =  -1 * ( 1.5*arrayV(i,columns-1) - 2.0*arrayV(i,columns-2) + 0.5*arrayV(i,columns-3) ) / gridSizeZ ; 
    }
    
    for ( Int_t j = columns-1 ; j >= 0 ; j-- ) {  // Count backwards to facilitate integration over Z
      // Integrate over Z
      for ( Int_t i = 0 ; i < rows ; i++ ) {
	Int_t index = 1 ;   // Simpsons rule if N=odd.  If N!=odd then add extra point by trapezoidal rule.  
	deltaEz(i,j) = 0.0 ;
	for ( Int_t k = j ; k < columns ; k++ ) {
	  deltaEz(i,j)  +=  index*(gridSizeZ/3.0)*arrayE(i,k) ;
	  if ( index != 4 )  index = 4; else index = 2 ;
	}
	if ( index == 4 ) deltaEz(i,j)  -=  (gridSizeZ/3.0)*arrayE(i,columns-1) ;
	if ( index == 2 ) deltaEz(i,j)  +=  
	  (gridSizeZ/3.0)*(0.5*arrayE(i,columns-2)-2.5*arrayE(i,columns-1)) ;
	if ( j == columns-2 ) deltaEz(i,j) =  
	  (gridSizeZ/3.0)*(1.5*arrayE(i,columns-2)+1.5*arrayE(i,columns-1)) ;
	if ( j == columns-1 ) deltaEz(i,j) =  0.0 ;
      }
    }
    // if ( m == 0 ) { TCanvas*  c1 =  new TCanvas("erOverEz","erOverEz",50,50,840,600) ;  c1 -> cd() ;
    // eroverEz.Draw("surf") ; } // JT test
    
    // calculate z distortion from the integrated Delta Ez residuals
    // and include the aquivalence (Volt to cm) of the ROC shift !!
    
    for ( Int_t j = 0 ; j < columns ; j++ )  {	  
      for ( Int_t i = 0 ; i < rows ; i++ ) {
	
	// Scale the Ez distortions with the drift velocity pertubation -> delivers cm
	deltaEz(i,j) = deltaEz(i,j)*fgkdvdE;
	
	// ROC Potential in cm aquivalent
	Double_t dzROCShift =  arrayV(i, columns -1)/ezField;  
	if ( rocDisplacement ) deltaEz(i,j) = deltaEz(i,j) + dzROCShift;  // add the ROC misaligment
	
      }
    }

  } // end loop over phi
  
 

  for ( Int_t k = 0 ; k < phislices ; k++ )
    {
      arrayofSumChargeDensities[k]->Delete() ;
    }
  


  arrayE.Clear();
}


Int_t AliTPCCorrection::IsPowerOfTwo(Int_t i) const {
  //
  // Helperfunction: Check if integer is a power of 2
  //
  Int_t j = 0;
  while( i > 0 ) { j += (i&1) ; i = (i>>1) ; }
  if ( j == 1 ) return(1) ;  // True
  return(0) ;                // False
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
  const Double_t kMaxR=500;
  const Double_t kMaxZ=500;
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  Int_t npoints1=0;
  Int_t npoints2=0;

  AliExternalTrackParam  track(trackIn); // 
  // generate points
  AliTrackPointArray pointArray0(npoints0);
  AliTrackPointArray pointArray1(npoints0);
  Double_t xyz[3];
  if (!AliTrackerBase::PropagateTrackToBxByBz(&track,kRTPC0,kMass,3,kTRUE,kMaxSnp)) return 0;
  //
  // simulate the track
  Int_t npoints=0;
  Float_t covPoint[6]={0,0,0, kSigmaY*kSigmaY,0,kSigmaZ*kSigmaZ};  //covariance at the local frame
  for (Double_t radius=kRTPC0; radius<kRTPC1; radius++){
    if (!AliTrackerBase::PropagateTrackToBxByBz(&track,radius,kMass,3,kTRUE,kMaxSnp)) return 0;
    track.GetXYZ(xyz);
    xyz[0]+=gRandom->Gaus(0,0.00005);
    xyz[1]+=gRandom->Gaus(0,0.00005);
    xyz[2]+=gRandom->Gaus(0,0.00005);
    if (TMath::Abs(track.GetZ())>kMaxZ) break;
    if (TMath::Abs(track.GetX())>kMaxR) break;
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
    if (!AliTrackerBase::PropagateTrackToBxByBz(track0,prot0.GetX(),kMass,3,kFALSE,kMaxSnp)) break;
    if (!AliTrackerBase::PropagateTrackToBxByBz(track1,prot0.GetX(),kMass,3,kFALSE,kMaxSnp)) break;
    if (TMath::Abs(track0->GetZ())>kMaxZ) break;
    if (TMath::Abs(track0->GetX())>kMaxR) break;
    if (TMath::Abs(track1->GetZ())>kMaxZ) break;
    if (TMath::Abs(track1->GetX())>kMaxR) break;

    track.GetXYZ(xyz);  // distorted track also propagated to the same reference radius
    //
    Double_t pointPos[2]={0,0};
    Double_t pointCov[3]={0,0,0};
    pointPos[0]=prot0.GetY();//local y
    pointPos[1]=prot0.GetZ();//local z
    pointCov[0]=prot0.GetCov()[3];//simay^2
    pointCov[1]=prot0.GetCov()[4];//sigmayz
    pointCov[2]=prot0.GetCov()[5];//sigmaz^2
    if (!track0->Update(pointPos,pointCov)) break;
    //
    Double_t deltaX=prot1.GetX()-prot0.GetX();   // delta X 
    Double_t deltaYX=deltaX*TMath::Tan(TMath::ASin(track1->GetSnp()));  // deltaY due  delta X
    Double_t deltaZX=deltaX*track1->GetTgl();                           // deltaZ due  delta X

    pointPos[0]=prot1.GetY()-deltaYX;//local y is sign correct? should be minus
    pointPos[1]=prot1.GetZ()-deltaZX;//local z is sign correct? should be minus
    pointCov[0]=prot1.GetCov()[3];//simay^2
    pointCov[1]=prot1.GetCov()[4];//sigmayz
    pointCov[2]=prot1.GetCov()[5];//sigmaz^2
    if (!track1->Update(pointPos,pointCov)) break;
    npoints1++;
    npoints2++;
  }
  if (npoints2<npoints)  return 0;
  AliTrackerBase::PropagateTrackToBxByBz(track0,refX,kMass,2.,kTRUE,kMaxSnp);
  track1->Rotate(track0->GetAlpha());
  AliTrackerBase::PropagateTrackToBxByBz(track1,refX,kMass,2.,kTRUE,kMaxSnp);

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
	// Double_t value=distPoint[2]-gz;
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


void AliTPCCorrection::FastSimDistortedVertex(Double_t orgVertex[3], Int_t nTracks, AliESDVertex &aV, AliESDVertex &avOrg, AliESDVertex &cV, AliESDVertex &cvOrg, TTreeSRedirector * const pcstream, Double_t etaCuts){
  //
  // Fast method to simulate the influence of the given distortion on the vertex reconstruction
  //

  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
  if (!magF) AliError("Magneticd field - not initialized");
  Double_t bz = magF->SolenoidField(); //field in kGauss
  printf("bz: %lf\n",bz);
  AliVertexerTracks *vertexer = new AliVertexerTracks(bz); // bz in kGauss

  TObjArray   aTrk;              // Original Track array of Aside
  TObjArray   daTrk;             // Distorted Track array of A side
  UShort_t    *aId = new UShort_t[nTracks];      // A side Track ID
  TObjArray   cTrk;               
  TObjArray   dcTrk;
  UShort_t    *cId = new UShort_t [nTracks];
  Int_t id=0; 
  Double_t mass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  TF1 fpt("fpt",Form("x*(1+(sqrt(x*x+%f^2)-%f)/([0]*[1]))^(-[0])",mass,mass),0.4,10);
  fpt.SetParameters(7.24,0.120);
  fpt.SetNpx(10000);
  for(Int_t nt=0; nt<nTracks; nt++){
    Double_t phi = gRandom->Uniform(0.0, 2*TMath::Pi());
    Double_t eta = gRandom->Uniform(-etaCuts, etaCuts);
    Double_t pt = fpt.GetRandom(); // momentum for f1
    //   printf("phi %lf  eta %lf pt %lf\n",phi,eta,pt);
    Short_t sign=1;
    if(gRandom->Rndm() < 0.5){
      sign =1;
    }else{
      sign=-1;
    }

    Double_t theta = 2*TMath::ATan(TMath::Exp(-eta))-TMath::Pi()/2.;
    Double_t pxyz[3];
    pxyz[0]=pt*TMath::Cos(phi);
    pxyz[1]=pt*TMath::Sin(phi);
    pxyz[2]=pt*TMath::Tan(theta);
    Double_t cv[21]={0};
    AliExternalTrackParam *t= new AliExternalTrackParam(orgVertex, pxyz, cv, sign);

    Double_t refX=1.;
    Int_t dir=-1;
    AliExternalTrackParam *td = FitDistortedTrack(*t, refX, dir,  NULL);
    if (!td) continue;
    if (pcstream) (*pcstream)<<"track"<<
      "eta="<<eta<<
      "theta="<<theta<<
      "tOrig.="<<t<<
      "td.="<<td<<
      "\n";
    if(( eta>0.07 )&&( eta<etaCuts )) { // - log(tan(0.5*theta)), theta = 0.5*pi - ATan(5.0/80.0)
      if (td){
	daTrk.AddLast(td);
	aTrk.AddLast(t);
	Int_t nn=aTrk.GetEntriesFast();
	aId[nn]=id;
      }
    }else if(( eta<-0.07 )&&( eta>-etaCuts )){
      if (td){
	dcTrk.AddLast(td);
	cTrk.AddLast(t);
	Int_t nn=cTrk.GetEntriesFast();
	cId[nn]=id;
      }
    }
    id++;  
  }// end of track loop

  vertexer->SetTPCMode();
  vertexer->SetConstraintOff();

  aV = *((AliESDVertex*)vertexer->FindPrimaryVertex(&daTrk,aId));  
  avOrg = *((AliESDVertex*)vertexer->FindPrimaryVertex(&aTrk,aId));
  cV = *((AliESDVertex*)vertexer->FindPrimaryVertex(&dcTrk,cId));  
  cvOrg = *((AliESDVertex*)vertexer->FindPrimaryVertex(&cTrk,cId));
  if (pcstream) (*pcstream)<<"vertex"<<
    "x="<<orgVertex[0]<<
    "y="<<orgVertex[1]<<
    "z="<<orgVertex[2]<<
    "av.="<<&aV<<              // distorted vertex A side
    "cv.="<<&cV<<              // distroted vertex C side
    "avO.="<<&avOrg<<         // original vertex A side
    "cvO.="<<&cvOrg<<
    "\n";
  delete []aId;
  delete []cId;
}

void AliTPCCorrection::AddVisualCorrection(AliTPCCorrection* corr, Int_t position){
  //
  // make correction available for visualization using 
  // TFormula, TFX and TTree::Draw 
  // important in order to check corrections and also compute dervied variables 
  // e.g correction partial derivatives
  //
  // NOTE - class is not owner of correction
  //     
  if (!fgVisualCorrection) fgVisualCorrection=new TObjArray;
  if (position!=0&&position>=fgVisualCorrection->GetEntriesFast())
    fgVisualCorrection->Expand(position*2);
  fgVisualCorrection->AddAt(corr, position);
}



Double_t AliTPCCorrection::GetCorrSector(Double_t sector, Double_t r, Double_t kZ, Int_t axisType, Int_t corrType){
  //
  // calculate the correction at given position - check the geffCorr
  //
  if (!fgVisualCorrection) return 0;
  AliTPCCorrection *corr = (AliTPCCorrection*)fgVisualCorrection->At(corrType);
  if (!corr) return 0;

  Double_t phi=sector*TMath::Pi()/9.;
  Double_t gx = r*TMath::Cos(phi);
  Double_t gy = r*TMath::Sin(phi);
  Double_t gz = r*kZ;
  Int_t nsector=(gz>0) ? 0:18; 
  //
  //
  //
  Float_t distPoint[3]={gx,gy,gz};
  corr->DistortPoint(distPoint, nsector);
  Double_t r0=TMath::Sqrt(gx*gx+gy*gy);
  Double_t r1=TMath::Sqrt(distPoint[0]*distPoint[0]+distPoint[1]*distPoint[1]);
  Double_t phi0=TMath::ATan2(gy,gx);
  Double_t phi1=TMath::ATan2(distPoint[1],distPoint[0]);
  if (axisType==0) return r1-r0;
  if (axisType==1) return (phi1-phi0)*r0;
  if (axisType==2) return distPoint[2]-gz;
  return phi1-phi0;
}

Double_t AliTPCCorrection::GetCorrXYZ(Double_t gx, Double_t gy, Double_t gz, Int_t axisType, Int_t corrType){
  //
  // return correction at given x,y,z
  // 
  if (!fgVisualCorrection) return 0;
  AliTPCCorrection *corr = (AliTPCCorrection*)fgVisualCorrection->At(corrType);
  if (!corr) return 0;
  Double_t phi0= TMath::ATan2(gy,gx);
  Int_t nsector=(gz>0) ? 0:18; 
  Float_t distPoint[3]={gx,gy,gz};
  corr->DistortPoint(distPoint, nsector);
  Double_t r0=TMath::Sqrt(gx*gx+gy*gy);
  Double_t r1=TMath::Sqrt(distPoint[0]*distPoint[0]+distPoint[1]*distPoint[1]);
  Double_t phi1=TMath::ATan2(distPoint[1],distPoint[0]);
  if (axisType==0) return r1-r0;
  if (axisType==1) return (phi1-phi0)*r0;
  if (axisType==2) return distPoint[2]-gz;
  return phi1-phi0;
}
