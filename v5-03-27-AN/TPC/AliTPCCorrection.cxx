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

// _________________________________________________________________
//
// Begin_Html
//   <h2>  AliTPCCorrection class   </h2>    
//  
//   The AliTPCCorrection class provides a general framework to deal with space point distortions. 
//   An correction class which inherits from here is for example AliTPCExBBShape or AliTPCExBTwist. <br> 
//   General virtual functions are (for example) CorrectPoint(x,roc) where x is the vector of initial 
//   positions in cartesian coordinates and roc represents the read-out chamber number according to 
//   the offline numbering convention. The vector x is overwritten with the corrected coordinates. <br> 
//   An alternative usage would be CorrectPoint(x,roc,dx), which leaves the vector x untouched, but 
//   returns the distortions via the vector dx. <br>
//   This class is normally used via the general class AliTPCComposedCorrection.   
//   <p>
//   Furthermore, the class contains basic geometrical descriptions like field cage radii 
//   (fgkIFCRadius, fgkOFCRadius) and length (fgkTPCZ0) plus the voltages. Also, the definitions 
//   of size and widths of the fulcrums building the grid of the final look-up table, which is 
//   then interpolated, is defined in kNX and fgkXList).
//   <p>
//   All physics-model classes below are derived from this class in order to not duplicate code 
//   and to allow a uniform treatment of all physics models.
//   <p>
//   <h3> Poisson solver </h3>    
//   A numerical solver of the Poisson equation (relaxation technique) is implemented for 2-dimensional 
//   geometries (r,z) as well as for 3-dimensional problems (r,$\phi$,z). The corresponding function 
//   names are PoissonRelaxation?D. The relevant function arguments are the arrays of the boundary and 
//   initial conditions (ArrayofArrayV, ArrayofChargeDensities) as well as the grid granularity which 
//   is used during the calculation. These inputs can be chosen according to the needs of the physical 
//   effect which is supposed to be simulated. In the 3D version, different symmetry conditions can be set
//   in order to reduce the calculation time (used in AliTPCFCVoltError3D).
//   <p>
//   <h3> Unified plotting functionality  </h3>    
//   Generic plot functions were implemented. They return a histogram pointer in the chosen plane of 
//   the TPC drift volume with a selectable grid granularity and the magnitude of the correction vector.
//   For example, the function CreateHistoDZinXY(z,nx,ny) returns a 2-dimensional histogram which contains 
//   the longitudinal corrections $dz$ in the (x,y)-plane at the given z position with the granularity of 
//   nx and ny. The magnitude of the corrections is defined by the class from which this function is called.
//   In the same manner, standard plots for the (r,$\phi$)-plane and for the other corrections like $dr$ and $rd\phi$ are available  
//   <p>                                                                      
//   Note: This class is normally used via the class AliTPCComposedCorrection
// End_Html
//
// Begin_Macro(source)
//   {
//   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
//   TCanvas *c2 = new TCanvas("cAliTPCCorrection","cAliTPCCorrection",700,1050);  c2->Divide(2,3);
//   AliTPCROCVoltError3D roc; // EXAMPLE PLOTS - SEE BELOW
//   roc.SetOmegaTauT1T2(0,1,1); // B=0
//   Float_t z0 = 1; // at +1 cm -> A side
//   c2->cd(1); roc.CreateHistoDRinXY(1.,300,300)->Draw("cont4z"); 
//   c2->cd(3);roc.CreateHistoDRPhiinXY(1.,300,300)->Draw("cont4z"); 
//   c2->cd(5);roc.CreateHistoDZinXY(1.,300,300)->Draw("cont4z"); 
//   Float_t phi0=0.5;
//   c2->cd(2);roc.CreateHistoDRinZR(phi0)->Draw("surf2"); 
//   c2->cd(4);roc.CreateHistoDRPhiinZR(phi0)->Draw("surf2"); 
//   c2->cd(6);roc.CreateHistoDZinZR(phi0)->Draw("surf2"); 
//   return c2;
//   } 
// End_Macro
//
// Begin_Html
//   <p>
//   Date: 27/04/2010  <br>
//   Authors: Magnus Mager, Stefan Rossegger, Jim Thomas                     
// End_Html 
// _________________________________________________________________


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

void AliTPCCorrection::DistortPointLocal(Float_t x[],const Short_t roc) {
  //
  // Distorts the initial coordinates x (cartesian coordinates)
  // according to the given effect (inherited classes)
  // roc represents the TPC read out chamber (offline numbering convention)
  //
  Float_t gxyz[3]={0,0,0};
  Double_t alpha = TMath::Pi()*(roc%18+0.5)/18;
  Double_t ca=TMath::Cos(alpha), sa= TMath::Sin(alpha);
  gxyz[0]=  ca*x[0]+sa*x[1];
  gxyz[1]= -sa*x[0]+ca*x[1];
  gxyz[2]= x[2];
  DistortPoint(gxyz,roc);
  x[0]=  ca*gxyz[0]-sa*gxyz[1];
  x[1]= +sa*gxyz[0]+ca*gxyz[1];
  x[2]= gxyz[2];
}
void AliTPCCorrection::CorrectPointLocal(Float_t x[],const Short_t roc) {
  //
  // Distorts the initial coordinates x (cartesian coordinates)
  // according to the given effect (inherited classes)
  // roc represents the TPC read out chamber (offline numbering convention)
  //
  Float_t gxyz[3]={0,0,0};
  Double_t alpha = TMath::Pi()*(roc%18+0.5)/18;
  Double_t ca=TMath::Cos(alpha), sa= TMath::Sin(alpha);
  gxyz[0]=  ca*x[0]+sa*x[1];
  gxyz[1]= -sa*x[0]+ca*x[1];
  gxyz[2]= x[2];
  CorrectPoint(gxyz,roc);
  x[0]=  ca*gxyz[0]-sa*gxyz[1];
  x[1]=  sa*gxyz[0]+ca*gxyz[1];
  x[2]=  gxyz[2];
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

Float_t AliTPCCorrection::Interpolate2DTable( const Int_t order, const Double_t x, const Double_t y, 
					      const Int_t nx,  const Int_t ny, const Double_t xv[], const Double_t yv[], 
					      const TMatrixF &array ) {
  //
  // Interpolate table (TMatrix format) - 2D interpolation
  // Float version (in order to decrease the OCDB size)
  //

  static  Int_t jlow = 0, klow = 0 ;
  Float_t saveArray[5] = {0.,0.,0.,0.,0.} ;

  Search( nx,  xv,  x,   jlow  ) ;
  Search( ny,  yv,  y,   klow  ) ;
  if ( jlow < 0 ) jlow = 0 ;   // check if out of range
  if ( klow < 0 ) klow = 0 ;
  if ( jlow + order  >=    nx - 1 ) jlow =   nx - 1 - order ;
  if ( klow + order  >=    ny - 1 ) klow =   ny - 1 - order ;

  for ( Int_t j = jlow ; j < jlow + order + 1 ; j++ )
    {
      Float_t *ajkl = &((TMatrixF&)array)(j,klow);
      saveArray[j-jlow]  = Interpolate( &yv[klow], ajkl , order, y )   ;
    }

  return( Interpolate( &xv[jlow], saveArray, order, x ) )   ;

}

Float_t AliTPCCorrection::Interpolate3DTable( const Int_t order, const Double_t x,   const Double_t y,   const Double_t z,
					      const Int_t  nx,    const Int_t  ny,    const Int_t  nz,
					      const Double_t xv[], const Double_t yv[], const Double_t zv[],
					      TMatrixF **arrayofArrays ) {
  //
  // Interpolate table (TMatrix format) - 3D interpolation 
  // Float version (in order to decrease the OCDB size)
  //

  static  Int_t ilow = 0, jlow = 0, klow = 0 ;
  Float_t saveArray[5]= {0.,0.,0.,0.,0.};
  Float_t savedArray[5]= {0.,0.,0.,0.,0.} ;

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
      TMatrixF &table = *arrayofArrays[k] ;
      for ( Int_t i = ilow ; i < ilow + order + 1 ; i++ )
	{
	  saveArray[i-ilow] = Interpolate( &yv[jlow], &table(i,jlow), order, y )   ;
	}
      savedArray[k-klow] = Interpolate( &xv[ilow], saveArray, order, x )   ; 
    }
  return( Interpolate( &zv[klow], savedArray, order, z ) )   ;

}
Float_t AliTPCCorrection::Interpolate( const Double_t xArray[], const Float_t yArray[], 
				       const Int_t order, const Double_t x ) {
  //
  // Interpolate function Y(x) using linear (order=1) or quadratic (order=2) interpolation.
  // Float version (in order to decrease the OCDB size)
  //

  Float_t y ;
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
  
  const Double_t kMaxZ0=220;
  const Double_t kZcut=3;
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  Int_t npoints1=0;
  Int_t npoints2=0;

  AliExternalTrackParam  track(trackIn); // 
  // generate points
  AliTrackPointArray pointArray0(npoints0);
  AliTrackPointArray pointArray1(npoints0);
  Double_t xyz[3];
  if (!AliTrackerBase::PropagateTrackTo(&track,kRTPC0,kMass,5,kTRUE,kMaxSnp)) return 0;
  //
  // simulate the track
  Int_t npoints=0;
  Float_t covPoint[6]={0,0,0, kSigmaY*kSigmaY,0,kSigmaZ*kSigmaZ};  //covariance at the local frame
  for (Double_t radius=kRTPC0; radius<kRTPC1; radius++){
    if (!AliTrackerBase::PropagateTrackTo(&track,radius,kMass,5,kTRUE,kMaxSnp)) return 0;
    track.GetXYZ(xyz);
    xyz[0]+=gRandom->Gaus(0,0.000005);
    xyz[1]+=gRandom->Gaus(0,0.000005);
    xyz[2]+=gRandom->Gaus(0,0.000005);
    if (TMath::Abs(track.GetZ())>kMaxZ0) continue;
    if (TMath::Abs(track.GetX())<kRTPC0) continue;
    if (TMath::Abs(track.GetX())>kRTPC1) continue;
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
  if (npoints<npoints0/4.) return 0;
  //
  // refit track
  //
  AliExternalTrackParam *track0=0;
  AliExternalTrackParam *track1=0;
  AliTrackPoint   point1,point2,point3;
  if (dir==1) {  //make seed inner
    pointArray0.GetPoint(point1,1);
    pointArray0.GetPoint(point2,11);
    pointArray0.GetPoint(point3,21);
  }
  if (dir==-1){ //make seed outer
    pointArray0.GetPoint(point1,npoints-21);
    pointArray0.GetPoint(point2,npoints-11);
    pointArray0.GetPoint(point3,npoints-1);
  } 
  if ((TMath::Abs(point1.GetX()-point3.GetX())+TMath::Abs(point1.GetY()-point3.GetY()))<10){
    printf("fit points not properly initialized\n");
    return 0;
  }
  track0 = AliTrackerBase::MakeSeed(point1, point2, point3);
  track1 = AliTrackerBase::MakeSeed(point1, point2, point3);
  track0->ResetCovariance(10);
  track1->ResetCovariance(10);
  if (TMath::Abs(AliTrackerBase::GetBz())<0.01){
    ((Double_t*)track0->GetParameter())[4]=  trackIn.GetParameter()[4];    
    ((Double_t*)track1->GetParameter())[4]=  trackIn.GetParameter()[4];
  }
  for (Int_t jpoint=0; jpoint<npoints; jpoint++){
    Int_t ipoint= (dir>0) ? jpoint: npoints-1-jpoint;
    //
    AliTrackPoint pIn0;
    AliTrackPoint pIn1;
    pointArray0.GetPoint(pIn0,ipoint);
    pointArray1.GetPoint(pIn1,ipoint);
    AliTrackPoint prot0 = pIn0.Rotate(track0->GetAlpha());   // rotate to the local frame - non distoted  point
    AliTrackPoint prot1 = pIn1.Rotate(track1->GetAlpha());   // rotate to the local frame -     distorted point
    if (TMath::Abs(prot0.GetX())<kRTPC0) continue;
    if (TMath::Abs(prot0.GetX())>kRTPC1) continue;
    //
    if (!AliTrackerBase::PropagateTrackTo(track0,prot0.GetX(),kMass,5,kFALSE,kMaxSnp)) break;
    if (!AliTrackerBase::PropagateTrackTo(track1,prot0.GetX(),kMass,5,kFALSE,kMaxSnp)) break;
    if (TMath::Abs(track0->GetZ())>kMaxZ) break;
    if (TMath::Abs(track0->GetX())>kMaxR) break;
    if (TMath::Abs(track1->GetZ())>kMaxZ) break;
    if (TMath::Abs(track1->GetX())>kMaxR) break;
    if (dir>0 && track1->GetX()>refX) continue;
    if (dir<0 && track1->GetX()<refX) continue;
    if (TMath::Abs(track1->GetZ())<kZcut)continue;
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
  if (npoints2<npoints/4.)  return 0;
  AliTrackerBase::PropagateTrackTo(track0,refX,kMass,5.,kTRUE,kMaxSnp);
  AliTrackerBase::PropagateTrackTo(track0,refX,kMass,1.,kTRUE,kMaxSnp);
  track1->Rotate(track0->GetAlpha());
  AliTrackerBase::PropagateTrackTo(track1,track0->GetX(),kMass,5.,kFALSE,kMaxSnp);

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




void AliTPCCorrection::MakeTrackDistortionTree(TTree *tinput, Int_t dtype, Int_t ptype, const TObjArray * corrArray, Int_t step, Int_t offset, Bool_t debug ){
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
  // dtype     - distortion type 0 - ITSTPC,  1 -TPCTRD, 2 - TPCvertex , 3 - TPC-TOF,  4 - TPCTPC track crossing 
  // ppype     - parameter type
  // corrArray - array with partial corrections
  // step      - skipe entries  - if 1 all entries processed - it is slow
  // debug     0 if debug on also space points dumped - it is slow

  const Double_t kMaxSnp = 0.85;  
  const Double_t kcutSnp=0.25;
  const Double_t kcutTheta=1.;
  const Double_t kRadiusTPC=85;
  //  AliTPCROC *tpcRoc =AliTPCROC::Instance();  
  //
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  //  const Double_t kB2C=-0.299792458e-3;
  const Int_t kMinEntries=20; 
  Double_t phi,theta, snp, mean,rms, entries,sector,dsec;
  Float_t refX;  
  Int_t run;
  tinput->SetBranchAddress("run",&run);
  tinput->SetBranchAddress("theta",&theta);
  tinput->SetBranchAddress("phi", &phi);
  tinput->SetBranchAddress("snp",&snp);
  tinput->SetBranchAddress("mean",&mean);
  tinput->SetBranchAddress("rms",&rms);
  tinput->SetBranchAddress("entries",&entries);
  tinput->SetBranchAddress("sector",&sector);
  tinput->SetBranchAddress("dsec",&dsec);
  tinput->SetBranchAddress("refX",&refX);
  TTreeSRedirector *pcstream = new TTreeSRedirector(Form("distortion%d_%d_%d.root",dtype,ptype,offset));
  //
  Int_t nentries=tinput->GetEntries();
  Int_t ncorr=corrArray->GetEntries();
  Double_t corrections[100]={0}; //
  Double_t tPar[5];
  Double_t cov[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t dir=0;
  if (dtype==5 || dtype==6) dtype=4;
  if (dtype==0) { dir=-1;}
  if (dtype==1) { dir=1;}
  if (dtype==2) { dir=-1;}
  if (dtype==3) { dir=1;}
  if (dtype==4) { dir=-1;}
  //
  for (Int_t ientry=offset; ientry<nentries; ientry+=step){
    tinput->GetEntry(ientry);
    if (TMath::Abs(snp)>kMaxSnp) continue;
    tPar[0]=0;
    tPar[1]=theta*refX;
    if (dtype==2)  tPar[1]=theta*kRadiusTPC;
    tPar[2]=snp;
    tPar[3]=theta;
    tPar[4]=(gRandom->Rndm()-0.5)*0.02;  // should be calculated - non equal to 0
    if (dtype==4){
      // tracks crossing CE
      tPar[1]=0;   // track at the CE
      //if (TMath::Abs(theta) <0.05) continue;  // deep cross
    }

    if (TMath::Abs(snp) >kcutSnp) continue;
    if (TMath::Abs(theta) >kcutTheta) continue;
    printf("%f\t%f\t%f\t%f\t%f\t%f\n",entries, sector,theta,snp, mean,rms);
    Double_t bz=AliTrackerBase::GetBz();
    if (dtype !=4) { //exclude TPC  - for TPC mainly non primary tracks
      if (dtype!=2 && TMath::Abs(bz)>0.1 )  tPar[4]=snp/(refX*bz*kB2C*2);
      
      if (dtype==2 && TMath::Abs(bz)>0.1 )  {
	tPar[4]=snp/(kRadiusTPC*bz*kB2C*2);//
	// snp at the TPC inner radius in case the vertex match used
      }
    }
    //
    tPar[4]+=(gRandom->Rndm()-0.5)*0.02;
    AliExternalTrackParam track(refX,phi,tPar,cov);
    Double_t xyz[3];
    track.GetXYZ(xyz);
    Int_t id=0;
    Double_t pt=1./tPar[4];
    Double_t dRrec=0; // dummy value - needed for points - e.g for laser
    //if (ptype==4 &&bz<0) mean*=-1;  // interpret as curvature -- COMMENTED out - in lookup signed 1/pt used
    Double_t refXD=refX;
    (*pcstream)<<"fit"<<
      "run="<<run<<       // run number
      "bz="<<bz<<         // magnetic filed used
      "dtype="<<dtype<<   // detector match type
      "ptype="<<ptype<<   // parameter type
      "theta="<<theta<<   // theta
      "phi="<<phi<<       // phi 
      "snp="<<snp<<       // snp
      "mean="<<mean<<     // mean dist value
      "rms="<<rms<<       // rms
      "sector="<<sector<<
      "dsec="<<dsec<<
      "refX="<<refXD<<         // referece X as double
      "gx="<<xyz[0]<<         // global position at reference
      "gy="<<xyz[1]<<         // global position at reference
      "gz="<<xyz[2]<<         // global position at reference	
      "dRrec="<<dRrec<<      // delta Radius in reconstruction
      "pt="<<pt<<            // pt
      "id="<<id<<             // track id
      "entries="<<entries;// number of entries in bin
    //
    Bool_t isOK=kTRUE;
    if (entries<kMinEntries) isOK=kFALSE;
    //
    if (dtype!=4) for (Int_t icorr=0; icorr<ncorr; icorr++) {
      AliTPCCorrection *corr = (AliTPCCorrection*)corrArray->At(icorr);
      corrections[icorr]=0;
      if (entries>kMinEntries){
	AliExternalTrackParam trackIn(refX,phi,tPar,cov);
	AliExternalTrackParam *trackOut = 0;
	if (debug) trackOut=corr->FitDistortedTrack(trackIn, refX, dir,pcstream);
	if (!debug) trackOut=corr->FitDistortedTrack(trackIn, refX, dir,0);
	if (dtype==0) {dir= -1;}
	if (dtype==1) {dir=  1;}
	if (dtype==2) {dir= -1;}
	if (dtype==3) {dir=  1;}
	//
	if (trackOut){
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn,refX,kMass,5,kTRUE,kMaxSnp)) isOK=kFALSE;
	  if (!trackOut->Rotate(trackIn.GetAlpha())) isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(trackOut,trackIn.GetX(),kMass,5,kFALSE,kMaxSnp)) isOK=kFALSE;
	  //	  trackOut->PropagateTo(trackIn.GetX(),AliTrackerBase::GetBz());
	  //	  
	  corrections[icorr]= trackOut->GetParameter()[ptype]-trackIn.GetParameter()[ptype];
	  delete trackOut;      
	}else{
	  corrections[icorr]=0;
	  isOK=kFALSE;
	}
	//if (ptype==4 &&bz<0) corrections[icorr]*=-1;  // interpret as curvature - commented out
      }      
      (*pcstream)<<"fit"<<
	Form("%s=",corr->GetName())<<corrections[icorr];   // dump correction value
    }
  
    if (dtype==4) for (Int_t icorr=0; icorr<ncorr; icorr++) {
      //
      // special case of the TPC tracks crossing the CE
      //
      AliTPCCorrection *corr = (AliTPCCorrection*)corrArray->At(icorr);
      corrections[icorr]=0;
      if (entries>kMinEntries){
	AliExternalTrackParam trackIn0(refX,phi,tPar,cov); //Outer - direction to vertex
	AliExternalTrackParam trackIn1(refX,phi,tPar,cov); //Inner - direction magnet 
	AliExternalTrackParam *trackOut0 = 0;
	AliExternalTrackParam *trackOut1 = 0;
	//
	if (debug)  trackOut0=corr->FitDistortedTrack(trackIn0, refX, dir,pcstream);
	if (!debug) trackOut0=corr->FitDistortedTrack(trackIn0, refX, dir,0);
	if (debug)  trackOut1=corr->FitDistortedTrack(trackIn1, refX, -dir,pcstream);
	if (!debug) trackOut1=corr->FitDistortedTrack(trackIn1, refX, -dir,0);
	//
	if (trackOut0 && trackOut1){
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn0,refX,kMass,5,kTRUE,kMaxSnp))  isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn0,refX,kMass,1,kFALSE,kMaxSnp)) isOK=kFALSE;
	  if (!trackOut0->Rotate(trackIn0.GetAlpha())) isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(trackOut0,trackIn0.GetX(),kMass,5,kFALSE,kMaxSnp)) isOK=kFALSE;
	  //
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn1,refX,kMass,5,kTRUE,kMaxSnp)) isOK=kFALSE;
	  if (!trackIn1.Rotate(trackIn0.GetAlpha()))  isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn1,trackIn0.GetX(),kMass,1,kFALSE,kMaxSnp)) isOK=kFALSE;
	  if (!trackOut1->Rotate(trackIn1.GetAlpha())) isOK=kFALSE;	  
	  if (!AliTrackerBase::PropagateTrackTo(trackOut1,trackIn1.GetX(),kMass,5,kFALSE,kMaxSnp)) isOK=kFALSE;
	  //
	  corrections[icorr] = (trackOut0->GetParameter()[ptype]-trackIn0.GetParameter()[ptype]);
	  corrections[icorr]-= (trackOut1->GetParameter()[ptype]-trackIn1.GetParameter()[ptype]);
	  if (isOK)
	    if ((TMath::Abs(trackOut0->GetX()-trackOut1->GetX())>0.1)||
		(TMath::Abs(trackOut0->GetX()-trackIn1.GetX())>0.1)||
		(TMath::Abs(trackOut0->GetAlpha()-trackOut1->GetAlpha())>0.00001)||
		(TMath::Abs(trackOut0->GetAlpha()-trackIn1.GetAlpha())>0.00001)||
		(TMath::Abs(trackIn0.GetTgl()-trackIn1.GetTgl())>0.0001)||
		(TMath::Abs(trackIn0.GetSnp()-trackIn1.GetSnp())>0.0001)
		){
	      isOK=kFALSE;
	    }	  	  
	  delete trackOut0;      
	  delete trackOut1;    	  
	}else{
	  corrections[icorr]=0;
	  isOK=kFALSE;
	}
	//
	//if (ptype==4 &&bz<0) corrections[icorr]*=-1;  // interpret as curvature - commented out no in lookup
      }      
      (*pcstream)<<"fit"<<
	Form("%s=",corr->GetName())<<corrections[icorr];   // dump correction value
    }
    //
    (*pcstream)<<"fit"<<"isOK="<<isOK<<"\n";
  }


  delete pcstream;
}



void AliTPCCorrection::MakeSectorDistortionTree(TTree *tinput, Int_t dtype, Int_t ptype, const TObjArray * corrArray, Int_t step, Int_t offset, Bool_t debug ){
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
  // dtype     - distortion type 10 - IROC-OROC 
  // ppype     - parameter type
  // corrArray - array with partial corrections
  // step      - skipe entries  - if 1 all entries processed - it is slow
  // debug     0 if debug on also space points dumped - it is slow

  const Double_t kMaxSnp = 0.8;  
  const Int_t kMinEntries=200; 
  //  AliTPCROC *tpcRoc =AliTPCROC::Instance();  
  //
  const Double_t kMass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  //  const Double_t kB2C=-0.299792458e-3;
  Double_t phi,theta, snp, mean,rms, entries,sector,dsec,globalZ;
  Int_t isec1, isec0;
  Double_t refXD;
  Float_t refX;
  Int_t run;
  tinput->SetBranchAddress("run",&run);
  tinput->SetBranchAddress("theta",&theta);
  tinput->SetBranchAddress("phi", &phi);
  tinput->SetBranchAddress("snp",&snp);
  tinput->SetBranchAddress("mean",&mean);
  tinput->SetBranchAddress("rms",&rms);
  tinput->SetBranchAddress("entries",&entries);
  tinput->SetBranchAddress("sector",&sector);
  tinput->SetBranchAddress("dsec",&dsec);
  tinput->SetBranchAddress("refX",&refXD);
  tinput->SetBranchAddress("z",&globalZ);
  tinput->SetBranchAddress("isec0",&isec0);
  tinput->SetBranchAddress("isec1",&isec1);
  TTreeSRedirector *pcstream = new TTreeSRedirector(Form("distortionSector%d_%d_%d.root",dtype,ptype,offset));
  //
  Int_t nentries=tinput->GetEntries();
  Int_t ncorr=corrArray->GetEntries();
  Double_t corrections[100]={0}; //
  Double_t tPar[5];
  Double_t cov[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t dir=0;
  //
  for (Int_t ientry=offset; ientry<nentries; ientry+=step){
    tinput->GetEntry(ientry);
    refX=refXD;
    Int_t id=-1;
    if (TMath::Abs(TMath::Abs(isec0%18)-TMath::Abs(isec1%18))==0) id=1;  // IROC-OROC - opposite side
    if (TMath::Abs(TMath::Abs(isec0%36)-TMath::Abs(isec1%36))==0) id=2;  // IROC-OROC - same side
    if (dtype==10  && id==-1) continue;
    //
    dir=-1;
    tPar[0]=0;
    tPar[1]=globalZ;
    tPar[2]=snp;
    tPar[3]=theta;
    tPar[4]=(gRandom->Rndm()-0.1)*0.2;  //
    Double_t pt=1./tPar[4];
    //
    printf("%f\t%f\t%f\t%f\t%f\t%f\n",entries, sector,theta,snp, mean,rms);
    Double_t bz=AliTrackerBase::GetBz();
    AliExternalTrackParam track(refX,phi,tPar,cov);    
    Double_t xyz[3],xyzIn[3],xyzOut[3];
    track.GetXYZ(xyz);
    track.GetXYZAt(85,bz,xyzIn);    
    track.GetXYZAt(245,bz,xyzOut);    
    Double_t phiIn  = TMath::ATan2(xyzIn[1],xyzIn[0]);
    Double_t phiOut = TMath::ATan2(xyzOut[1],xyzOut[0]);
    Double_t phiRef = TMath::ATan2(xyz[1],xyz[0]);
    Int_t sectorRef = TMath::Nint(9.*phiRef/TMath::Pi()-0.5);
    Int_t sectorIn  = TMath::Nint(9.*phiIn/TMath::Pi()-0.5);
    Int_t sectorOut = TMath::Nint(9.*phiOut/TMath::Pi()-0.5);
    //
    Bool_t isOK=kTRUE; 
    if (sectorIn!=sectorOut) isOK=kFALSE;  // requironment - cluster in the same sector
    if (sectorIn!=sectorRef) isOK=kFALSE;  // requironment - cluster in the same sector
    if (entries<kMinEntries/(1+TMath::Abs(globalZ/100.))) isOK=kFALSE;  // requironment - minimal amount of tracks in bin
    // Do downscale
    if (TMath::Abs(theta)>1) isOK=kFALSE;
    //
    Double_t dRrec=0; // dummy value - needed for points - e.g for laser
    //
    (*pcstream)<<"fit"<<
      "run="<<run<<       //run
      "bz="<<bz<<         // magnetic filed used
      "dtype="<<dtype<<   // detector match type
      "ptype="<<ptype<<   // parameter type
      "theta="<<theta<<   // theta
      "phi="<<phi<<       // phi 
      "snp="<<snp<<       // snp
      "mean="<<mean<<     // mean dist value
      "rms="<<rms<<       // rms
      "sector="<<sector<<
      "dsec="<<dsec<<
      "refX="<<refXD<<         // referece X
      "gx="<<xyz[0]<<         // global position at reference
      "gy="<<xyz[1]<<         // global position at reference
      "gz="<<xyz[2]<<         // global position at reference	
      "dRrec="<<dRrec<<      // delta Radius in reconstruction
      "pt="<<pt<<      //pt
      "id="<<id<<             // track id
      "entries="<<entries;// number of entries in bin
    //
    AliExternalTrackParam *trackOut0 = 0;
    AliExternalTrackParam *trackOut1 = 0;
    AliExternalTrackParam *ptrackIn0 = 0;
    AliExternalTrackParam *ptrackIn1 = 0;

    for (Int_t icorr=0; icorr<ncorr; icorr++) {
      //
      // special case of the TPC tracks crossing the CE
      //
      AliTPCCorrection *corr = (AliTPCCorrection*)corrArray->At(icorr);
      corrections[icorr]=0;
      if (entries>kMinEntries &&isOK){
	AliExternalTrackParam trackIn0(refX,phi,tPar,cov);
	AliExternalTrackParam trackIn1(refX,phi,tPar,cov);
	ptrackIn1=&trackIn0;
	ptrackIn0=&trackIn1;
	//
	if (debug)  trackOut0=corr->FitDistortedTrack(trackIn0, refX, dir,pcstream);
	if (!debug) trackOut0=corr->FitDistortedTrack(trackIn0, refX, dir,0);
	if (debug)  trackOut1=corr->FitDistortedTrack(trackIn1, refX, -dir,pcstream);
	if (!debug) trackOut1=corr->FitDistortedTrack(trackIn1, refX, -dir,0);
	//
	if (trackOut0 && trackOut1){
	  //
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn0,refX,kMass,1,kTRUE,kMaxSnp))  isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn0,refX,kMass,1,kFALSE,kMaxSnp)) isOK=kFALSE;
	  // rotate all tracks to the same frame
	  if (!trackOut0->Rotate(trackIn0.GetAlpha())) isOK=kFALSE;
	  if (!trackIn1.Rotate(trackIn0.GetAlpha()))  isOK=kFALSE;
	  if (!trackOut1->Rotate(trackIn0.GetAlpha())) isOK=kFALSE;	  
	  //
	  if (!AliTrackerBase::PropagateTrackTo(trackOut0,refX,kMass,1,kFALSE,kMaxSnp)) isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(&trackIn1,refX,kMass,1,kFALSE,kMaxSnp)) isOK=kFALSE;
	  if (!AliTrackerBase::PropagateTrackTo(trackOut1,refX,kMass,1,kFALSE,kMaxSnp)) isOK=kFALSE;
	  //
	  corrections[icorr] = (trackOut0->GetParameter()[ptype]-trackIn0.GetParameter()[ptype]);
	  corrections[icorr]-= (trackOut1->GetParameter()[ptype]-trackIn1.GetParameter()[ptype]);
	  (*pcstream)<<"fitDebug"<< // just to debug the correction
	    "mean="<<mean<<
	    "pIn0.="<<ptrackIn0<<
	    "pIn1.="<<ptrackIn1<<
	    "pOut0.="<<trackOut0<<
	    "pOut1.="<<trackOut1<<
	    "refX="<<refXD<<
	    "\n";
	  delete trackOut0;      
	  delete trackOut1;      
	}else{
	  corrections[icorr]=0;
	  isOK=kFALSE;
	}
      }      
      (*pcstream)<<"fit"<<
	Form("%s=",corr->GetName())<<corrections[icorr];   // dump correction value
    }
    //
    (*pcstream)<<"fit"<<"isOK="<<isOK<<"\n";
  }
  delete pcstream;
}



void AliTPCCorrection::MakeLaserDistortionTreeOld(TTree* tree, TObjArray *corrArray, Int_t itype){
  //
  // Make a laser fit tree for global minimization
  //
  const Double_t cutErrY=0.1;
  const Double_t cutErrZ=0.1;
  const Double_t kEpsilon=0.00000001;
  const Double_t kMaxDist=1.;  // max distance - space correction
  const Double_t kMaxRMS=0.05;  // max distance -between point and local mean
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
  TTreeSRedirector *pcstream= new TTreeSRedirector("distortionLaser_0.root");
  Double_t bz=AliTrackerBase::GetBz();
  // 

  for (Int_t ientry=0; ientry<entries; ientry++){
    tree->GetEntry(ientry);
    if (!ltr->GetVecGX()){
      ltr->UpdatePoints();
    }
    TVectorD * delta= (itype==0)? vecdY:vecdZ;
    TVectorD * err= (itype==0)? veceY:veceZ;
    TLinearFitter  fitter(2,"pol1");
    for (Int_t iter=0; iter<2; iter++){
      Double_t kfit0=0, kfit1=0;
      Int_t npoints=fitter.GetNpoints();
      if (npoints>80){
	fitter.Eval();
	kfit0=fitter.GetParameter(0);
	kfit1=fitter.GetParameter(1);
      }
      for (Int_t irow=0; irow<159; irow++){
	Bool_t isOK=kTRUE;
	Int_t isOKF=0;
	Int_t nentries = 1000;
	if (veceY->GetMatrixArray()[irow]>cutErrY||veceZ->GetMatrixArray()[irow]>cutErrZ) nentries=0;
	if (veceY->GetMatrixArray()[irow]<kEpsilon||veceZ->GetMatrixArray()[irow]<kEpsilon) nentries=0;
	Int_t dtype=5;
	Double_t array[10];
	Int_t first3=TMath::Max(irow-3,0);
	Int_t last3 =TMath::Min(irow+3,159);
	Int_t counter=0;
	if ((*ltr->GetVecSec())[irow]>=0 && err) {
	  for (Int_t jrow=first3; jrow<=last3; jrow++){
	    if ((*ltr->GetVecSec())[irow]!= (*ltr->GetVecSec())[jrow]) continue;
	    if ((*err)[jrow]<kEpsilon) continue;
	    array[counter]=(*delta)[jrow];
	    counter++;
	  }
	}    
	Double_t rms3  = 0;
	Double_t mean3 = 0;
	if (counter>2){
	  rms3  = TMath::RMS(counter,array);
	  mean3  = TMath::Mean(counter,array);
	}else{
	  isOK=kFALSE;
	}
	Double_t phi   =(*ltr->GetVecPhi())[irow];
	Double_t theta =ltr->GetTgl();
	Double_t mean=delta->GetMatrixArray()[irow];
	Double_t gx=0,gy=0,gz=0;
	Double_t snp = (*ltr->GetVecP2())[irow];
	Double_t dRrec=0;
	//      Double_t rms = err->GetMatrixArray()[irow];
	//
	gx = (*ltr->GetVecGX())[irow];
	gy = (*ltr->GetVecGY())[irow];
	gz = (*ltr->GetVecGZ())[irow];
	//
	// get delta R used in reconstruction
	AliTPCcalibDB*  calib=AliTPCcalibDB::Instance();  
	AliTPCCorrection * correction = calib->GetTPCComposedCorrection(AliTrackerBase::GetBz());
	//      const AliTPCRecoParam * recoParam = calib->GetTransform()->GetCurrentRecoParam();
	//Double_t xyz0[3]={gx,gy,gz};
	Double_t oldR=TMath::Sqrt(gx*gx+gy*gy);
	Double_t fphi = TMath::ATan2(gy,gx);      
	Double_t fsector = 9.*fphi/TMath::Pi();
	if (fsector<0) fsector+=18;
	Double_t dsec = fsector-Int_t(fsector)-0.5;
	Double_t refX=0;
	Int_t id= ltr->GetId();
	Double_t pt=0;
	//
	if (1 && oldR>1) {
	  Float_t xyz1[3]={gx,gy,gz};
	  Int_t sector=(gz>0)?0:18;
	  correction->CorrectPoint(xyz1, sector);
	  refX=TMath::Sqrt(xyz1[0]*xyz1[0]+xyz1[1]*xyz1[1]);
	  dRrec=oldR-refX;
	} 
	if (TMath::Abs(rms3)>kMaxRMS) isOK=kFALSE;
	if (TMath::Abs(mean-mean3)>kMaxRMS) isOK=kFALSE;
	if (counter<4) isOK=kFALSE;	
	if (npoints<90) isOK=kFALSE;	
	if (isOK){
	  fitter.AddPoint(&refX,mean);
	}
	Double_t deltaF=kfit0+kfit1*refX;
	if (iter==1){
	  (*pcstream)<<"fitFull"<<  // dumpe also intermediate results
	    "bz="<<bz<<         // magnetic filed used
	    "dtype="<<dtype<<   // detector match type
	    "ptype="<<itype<<   // parameter type
	    "theta="<<theta<<   // theta
	    "phi="<<phi<<       // phi 
	    "snp="<<snp<<       // snp
	    "mean="<<mean3<<     // mean dist value
	    "rms="<<rms3<<       // rms
	    "deltaF="<<deltaF<<
	    "npoints="<<npoints<<  //number of points
	    "mean3="<<mean3<<     // mean dist value
	    "rms3="<<rms3<<       // rms
	    "counter="<<counter<<
	    "sector="<<fsector<<
	    "dsec="<<dsec<<
	    //
	    "refX="<<refX<<      // reference radius
	    "gx="<<gx<<         // global position
	    "gy="<<gy<<         // global position
	    "gz="<<gz<<         // global position
	    "dRrec="<<dRrec<<      // delta Radius in reconstruction
	    "id="<<id<<     //bundle	
	    "entries="<<nentries<<// number of entries in bin
	    "\n";
	}
	if (iter==1) (*pcstream)<<"fit"<<  // dump valus for fit
	  "bz="<<bz<<         // magnetic filed used
	  "dtype="<<dtype<<   // detector match type
	  "ptype="<<itype<<   // parameter type
	  "theta="<<theta<<   // theta
	  "phi="<<phi<<       // phi 
	  "snp="<<snp<<       // snp
	  "mean="<<mean3<<     // mean dist value
	  "rms="<<rms3<<       // rms
	  "sector="<<fsector<<
	  "dsec="<<dsec<<
	  //
	  "refX="<<refX<<      // reference radius
	  "gx="<<gx<<         // global position
	  "gy="<<gy<<         // global position
	  "gz="<<gz<<         // global position
	  "dRrec="<<dRrec<<      // delta Radius in reconstruction
	  "pt="<<pt<<           //pt
	  "id="<<id<<     //bundle	
	  "entries="<<nentries;// number of entries in bin
	//
	//    
	Double_t ky = TMath::Tan(TMath::ASin(snp));
	Int_t ncorr = corrArray->GetEntries();
	Double_t r0   = TMath::Sqrt(gx*gx+gy*gy);
	Double_t phi0 = TMath::ATan2(gy,gx);
	Double_t distortions[1000]={0};
	Double_t distortionsR[1000]={0};
	if (iter==1){
	  for (Int_t icorr=0; icorr<ncorr; icorr++) {
	    AliTPCCorrection *corr = (AliTPCCorrection*)corrArray->At(icorr);
	    Float_t distPoint[3]={gx,gy,gz}; 
	    Int_t sector= (gz>0)? 0:18;
	    if (r0>80){
	      corr->DistortPoint(distPoint, sector);
	    }
	    // Double_t value=distPoint[2]-gz;
	    if (itype==0 && r0>1){
	      Double_t r1   = TMath::Sqrt(distPoint[0]*distPoint[0]+distPoint[1]*distPoint[1]);
	      Double_t phi1 = TMath::ATan2(distPoint[1],distPoint[0]);
	      Double_t drphi= r0*(phi1-phi0);
	      Double_t dr   = r1-r0;
	      distortions[icorr]  = drphi-ky*dr;
	      distortionsR[icorr] = dr;
	    }
	    if (TMath::Abs(distortions[icorr])>kMaxDist) {isOKF=icorr+1; isOK=kFALSE; }
	    if (TMath::Abs(distortionsR[icorr])>kMaxDist) {isOKF=icorr+1; isOK=kFALSE;}
	    (*pcstream)<<"fit"<<
	      Form("%s=",corr->GetName())<<distortions[icorr];    // dump correction value
	  }
	  (*pcstream)<<"fit"<<"isOK="<<isOK<<"\n";
	}
      }
    }
  }
  delete pcstream;
}



void   AliTPCCorrection::MakeDistortionMap(THnSparse * his0, TTreeSRedirector * const pcstream, const char* hname, Int_t run, Float_t refX, Int_t type, Int_t integ){
  //
  // make a distortion map out ou fthe residual histogram
  // Results are written to the debug streamer - pcstream
  // Parameters:
  //   his0       - input (4D) residual histogram
  //   pcstream   - file to write the tree
  //   run        - run number
  //   refX       - track matching reference X
  //   type       - 0- y 1-z,2 -snp, 3-theta, 4=1/pt
  // THnSparse axes:
  // OBJ: TAxis     #Delta  #Delta
  // OBJ: TAxis     tanTheta        tan(#Theta)
  // OBJ: TAxis     phi     #phi
  // OBJ: TAxis     snp     snp

  // marian.ivanov@cern.ch
  const Int_t kMinEntries=10;
  Double_t bz=AliTrackerBase::GetBz();
  Int_t idim[4]={0,1,2,3};
  //
  //
  //
  Int_t nbins3=his0->GetAxis(3)->GetNbins();
  Int_t first3=his0->GetAxis(3)->GetFirst();
  Int_t last3 =his0->GetAxis(3)->GetLast();
  //
  for (Int_t ibin3=first3; ibin3<last3; ibin3+=1){   // axis 3 - local angle
    his0->GetAxis(3)->SetRange(TMath::Max(ibin3-integ,1),TMath::Min(ibin3+integ,nbins3));
    Double_t      x3= his0->GetAxis(3)->GetBinCenter(ibin3);
    THnSparse * his3= his0->Projection(3,idim);         //projected histogram according selection 3
    //
    Int_t nbins2    = his3->GetAxis(2)->GetNbins();
    Int_t first2    = his3->GetAxis(2)->GetFirst();
    Int_t last2     = his3->GetAxis(2)->GetLast();
    //
    for (Int_t ibin2=first2; ibin2<last2; ibin2+=1){   // axis 2 - phi
      his3->GetAxis(2)->SetRange(TMath::Max(ibin2-integ,1),TMath::Min(ibin2+integ,nbins2));
      Double_t      x2= his3->GetAxis(2)->GetBinCenter(ibin2);
      THnSparse * his2= his3->Projection(2,idim);         //projected histogram according selection 2
      Int_t nbins1     = his2->GetAxis(1)->GetNbins();
      Int_t first1     = his2->GetAxis(1)->GetFirst();
      Int_t last1      = his2->GetAxis(1)->GetLast();
      for (Int_t ibin1=first1; ibin1<last1; ibin1++){   //axis 1 - theta
	//
	Double_t       x1= his2->GetAxis(1)->GetBinCenter(ibin1);
	his2->GetAxis(1)->SetRange(TMath::Max(ibin1-1,1),TMath::Min(ibin1+1,nbins1));
	if (TMath::Abs(x1)<0.1){
	  if (x1<0) his2->GetAxis(1)->SetRange(TMath::Max(ibin1-1,1),TMath::Min(ibin1,nbins1));
	  if (x1>0) his2->GetAxis(1)->SetRange(TMath::Max(ibin1,1),TMath::Min(ibin1+1,nbins1));
	}
	if (TMath::Abs(x1)<0.06){
	  his2->GetAxis(1)->SetRange(TMath::Max(ibin1,1),TMath::Min(ibin1,nbins1));
	}
	TH1 * hisDelta = his2->Projection(0);
	//
	Double_t entries = hisDelta->GetEntries();
	Double_t mean=0, rms=0;
	if (entries>kMinEntries){
	  mean    = hisDelta->GetMean(); 
	  rms = hisDelta->GetRMS(); 
	}
	Double_t sector = 9.*x2/TMath::Pi();
	if (sector<0) sector+=18;
	Double_t dsec = sector-Int_t(sector)-0.5;
	Double_t z=refX*x1;
	(*pcstream)<<hname<<
	  "run="<<run<<
	  "bz="<<bz<<
	  "theta="<<x1<<
	  "phi="<<x2<<
	  "z="<<z<<            // dummy z
	  "snp="<<x3<<
	  "entries="<<entries<<
	  "mean="<<mean<<
	  "rms="<<rms<<
	  "refX="<<refX<<   // track matching refernce plane
	  "type="<<type<<   //
	  "sector="<<sector<<
	  "dsec="<<dsec<<
	  "\n";
	delete hisDelta;
	//printf("%f\t%f\t%f\t%f\t%f\n",x3,x2,x1, entries,mean);
      }
      delete his2;
    }
    delete his3;
  }
}




void   AliTPCCorrection::MakeDistortionMapCosmic(THnSparse * hisInput, TTreeSRedirector * const pcstream, const char* hname, Int_t run, Float_t refX, Int_t type){
  //
  // make a distortion map out ou fthe residual histogram
  // Results are written to the debug streamer - pcstream
  // Parameters:
  //   his0       - input (4D) residual histogram
  //   pcstream   - file to write the tree
  //   run        - run number
  //   refX       - track matching reference X
  //   type       - 0- y 1-z,2 -snp, 3-theta, 4=1/pt
  // marian.ivanov@cern.ch
  //
  //  Histo axeses
  //   Collection name='TObjArray', class='TObjArray', size=16
  //  0. OBJ: TAxis     #Delta  #Delta
  //  1. OBJ: TAxis     N_{cl}  N_{cl}
  //  2. OBJ: TAxis     dca_{r} (cm)    dca_{r} (cm)
  //  3. OBJ: TAxis     z (cm)  z (cm)
  //  4. OBJ: TAxis     sin(#phi)       sin(#phi)
  //  5. OBJ: TAxis     tan(#theta)     tan(#theta)
  //  6. OBJ: TAxis     1/pt (1/GeV)    1/pt (1/GeV)
  //  7. OBJ: TAxis     pt (GeV)        pt (GeV)
  //  8. OBJ: TAxis     alpha   alpha
  const Int_t kMinEntries=10;
  //
  //  1. make default selections
  //
  TH1 * hisDelta=0;
  Int_t idim0[4]={0 , 5, 8,  3};   // delta, theta, alpha, z
  hisInput->GetAxis(1)->SetRangeUser(110,190);   //long tracks
  hisInput->GetAxis(2)->SetRangeUser(-10,35);    //tracks close to beam pipe
  hisInput->GetAxis(4)->SetRangeUser(-0.3,0.3); //small snp at TPC entrance
  hisInput->GetAxis(7)->SetRangeUser(3,100); //"high pt tracks"
  hisDelta= hisInput->Projection(0);
  hisInput->GetAxis(0)->SetRangeUser(-6.*hisDelta->GetRMS(), +6.*hisDelta->GetRMS());
  delete hisDelta;
  THnSparse *his0=  hisInput->Projection(4,idim0);
  //
  // 2. Get mean in diferent bins
  //
  Int_t nbins1=his0->GetAxis(1)->GetNbins();
  Int_t first1=his0->GetAxis(1)->GetFirst();
  Int_t last1 =his0->GetAxis(1)->GetLast();
  //
  Double_t bz=AliTrackerBase::GetBz();
  Int_t idim[4]={0,1, 2,  3};  // delta, theta,alpha,z
  //
  for (Int_t ibin1=first1; ibin1<=last1; ibin1++){   //axis 1 - theta
    //
    Double_t       x1= his0->GetAxis(1)->GetBinCenter(ibin1);  
    his0->GetAxis(1)->SetRange(TMath::Max(ibin1-1,1),TMath::Min(ibin1+1,nbins1));
    //
    THnSparse * his1 = his0->Projection(4,idim);  // projected histogram according range1
    Int_t nbins3     = his1->GetAxis(3)->GetNbins();
    Int_t first3     = his1->GetAxis(3)->GetFirst();
    Int_t last3      = his1->GetAxis(3)->GetLast();
    //
    for (Int_t ibin3=first3-1; ibin3<=last3; ibin3+=1){   // axis 3 - z at "vertex"
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
      for (Int_t ibin2=first2; ibin2<=last2; ibin2+=1){
	his3->GetAxis(2)->SetRange(TMath::Max(ibin2-1,1),TMath::Min(ibin2+1,nbins2));
	Double_t x2= his3->GetAxis(2)->GetBinCenter(ibin2);
	hisDelta = his3->Projection(0);
	//
	Double_t entries = hisDelta->GetEntries();
	Double_t mean=0, rms=0;
	if (entries>kMinEntries){
	  mean    = hisDelta->GetMean(); 
	  rms = hisDelta->GetRMS(); 
	}
	Double_t sector = 9.*x2/TMath::Pi();
	if (sector<0) sector+=18;
	Double_t dsec = sector-Int_t(sector)-0.5;
	Double_t snp=0;  // dummy snp - equal 0
	(*pcstream)<<hname<<
	  "run="<<run<<
	  "bz="<<bz<<            // magnetic field
	  "theta="<<x1<<         // theta
	  "phi="<<x2<<           // phi (alpha)
	  "z="<<x3<<             // z at "vertex"
	  "snp="<<snp<<          // dummy snp
	  "entries="<<entries<<  // entries in bin
	  "mean="<<mean<<        // mean
	  "rms="<<rms<<
	  "refX="<<refX<<        // track matching refernce plane
	  "type="<<type<<        // parameter type
	  "sector="<<sector<<    // sector
	  "dsec="<<dsec<<        // dummy delta sector
	  "\n";
	delete hisDelta;
	printf("%f\t%f\t%f\t%f\t%f\n",x1,x3,x2, entries,mean);
      }
      delete his3;
    }
    delete his1;
  }
  delete his0;
}



void   AliTPCCorrection::MakeDistortionMapSector(THnSparse * hisInput, TTreeSRedirector * const pcstream, const char* hname, Int_t run, Int_t type){
  //
  // make a distortion map out of the residual histogram
  // Results are written to the debug streamer - pcstream
  // Parameters:
  //   his0       - input (4D) residual histogram
  //   pcstream   - file to write the tree
  //   run        - run number
  //   type       - 0- y 1-z,2 -snp, 3-theta
  // marian.ivanov@cern.ch

  //Collection name='TObjArray', class='TObjArray', size=16
  //0  OBJ: TAxis     delta   delta
  //1  OBJ: TAxis     phi     phi
  //2  OBJ: TAxis     localX  localX
  //3  OBJ: TAxis     kY      kY
  //4  OBJ: TAxis     kZ      kZ
  //5  OBJ: TAxis     is1     is1
  //6  OBJ: TAxis     is0     is0
  //7. OBJ: TAxis     z       z
  //8. OBJ: TAxis     IsPrimary       IsPrimary

  const Int_t kMinEntries=10;
  THnSparse * hisSector0=0;
  TH1 * htemp=0;    // histogram to calculate mean value of parameter
  Double_t bz=AliTrackerBase::GetBz();

  //
  // Loop over pair of sector:
  // isPrim         - 8  ==> 8
  // isec0          - 6  ==> 7
  //   isec1        - 5  ==> 6
  //     refX       - 2  ==> 5
  //
  //     phi        - 1  ==> 4
  //       z        - 7  ==> 3
  //         snp    - 3  ==> 2
  //           theta- 4  ==> 1
  //                  0  ==> 0;           
  for (Int_t isec0=0; isec0<72; isec0++){
    Int_t index0[9]={0, 4, 3, 7, 1, 2, 5, 6,8}; //regroup indeces
    //
    //hisInput->GetAxis(8)->SetRangeUser(-0.1,0.4);  // select secondaries only ? - get out later ?
    hisInput->GetAxis(6)->SetRangeUser(isec0-0.1,isec0+0.1);
    hisSector0=hisInput->Projection(7,index0);
    //
    //
    for (Int_t isec1=isec0+1; isec1<72; isec1++){    
      //if (isec1!=isec0+36) continue;
      if ( TMath::Abs((isec0%18)-(isec1%18))>1.5 && TMath::Abs((isec0%18)-(isec1%18))<16.5) continue;
      printf("Sectors %d\t%d\n",isec1,isec0);
      hisSector0->GetAxis(6)->SetRangeUser(isec1-0.1,isec1+0.1);      
      TH1 * hisX=hisSector0->Projection(5);
      Double_t refX= hisX->GetMean();
      delete hisX;
      TH1 *hisDelta=hisSector0->Projection(0);
      Double_t dmean = hisDelta->GetMean();
      Double_t drms = hisDelta->GetRMS();
      hisSector0->GetAxis(0)->SetRangeUser(dmean-5.*drms, dmean+5.*drms);
      delete hisDelta;
      //
      //  1. make default selections
      //
      Int_t idim0[5]={0 , 1, 2, 3, 4}; // {delta, theta, snp, z, phi }
      THnSparse *hisSector1=  hisSector0->Projection(5,idim0);
      //
      // 2. Get mean in diferent bins
      //
      Int_t idim[5]={0, 1, 2,  3, 4};  // {delta, theta-1,snp-2 ,z-3, phi-4}
      //
      //      Int_t nbinsPhi=hisSector1->GetAxis(4)->GetNbins();
      Int_t firstPhi=hisSector1->GetAxis(4)->GetFirst();
      Int_t lastPhi =hisSector1->GetAxis(4)->GetLast();
      //
      for (Int_t ibinPhi=firstPhi; ibinPhi<=lastPhi; ibinPhi+=1){   //axis 4 - phi
	//
	// Phi loop
	//
	Double_t       xPhi= hisSector1->GetAxis(4)->GetBinCenter(ibinPhi);         
	Double_t psec    = (9*xPhi/TMath::Pi());
	if (psec<0) psec+=18;
	Bool_t isOK0=kFALSE;
	Bool_t isOK1=kFALSE;
	if (TMath::Abs(psec-isec0%18-0.5)<1. || TMath::Abs(psec-isec0%18-17.5)<1.)  isOK0=kTRUE;
	if (TMath::Abs(psec-isec1%18-0.5)<1. || TMath::Abs(psec-isec1%18-17.5)<1.)  isOK1=kTRUE;
	if (!isOK0) continue;
	if (!isOK1) continue;
	//
	hisSector1->GetAxis(4)->SetRange(TMath::Max(ibinPhi-2,firstPhi),TMath::Min(ibinPhi+2,lastPhi));
	if (isec1!=isec0+36) {
	  hisSector1->GetAxis(4)->SetRange(TMath::Max(ibinPhi-3,firstPhi),TMath::Min(ibinPhi+3,lastPhi));
	}
	//
	htemp = hisSector1->Projection(4);
	xPhi=htemp->GetMean();
	delete htemp;
	THnSparse * hisPhi = hisSector1->Projection(4,idim);
	//Int_t nbinsZ     = hisPhi->GetAxis(3)->GetNbins();
	Int_t firstZ     = hisPhi->GetAxis(3)->GetFirst();
	Int_t lastZ      = hisPhi->GetAxis(3)->GetLast();
	//
	for (Int_t ibinZ=firstZ; ibinZ<=lastZ; ibinZ+=1){   // axis 3 - z
	  //
	  // Z loop
	  //
	  hisPhi->GetAxis(3)->SetRange(TMath::Max(ibinZ,firstZ),TMath::Min(ibinZ,lastZ));
	  if (isec1!=isec0+36) {
	    hisPhi->GetAxis(3)->SetRange(TMath::Max(ibinZ-1,firstZ),TMath::Min(ibinZ-1,lastZ));	    
	  }
	  htemp = hisPhi->Projection(3);
	  Double_t      xZ= htemp->GetMean();
	  delete htemp;
	  THnSparse * hisZ= hisPhi->Projection(3,idim);         
	  //projected histogram according selection 3 -z
	  //
	  //
	  //Int_t nbinsSnp    = hisZ->GetAxis(2)->GetNbins();
	  Int_t firstSnp    = hisZ->GetAxis(2)->GetFirst();
	  Int_t lastSnp     = hisZ->GetAxis(2)->GetLast();
	  for (Int_t ibinSnp=firstSnp; ibinSnp<=lastSnp; ibinSnp+=2){   // axis 2 - snp
	    //
	    // Snp loop
	    //
	    hisZ->GetAxis(2)->SetRange(TMath::Max(ibinSnp-1,firstSnp),TMath::Min(ibinSnp+1,lastSnp));
	    if (isec1!=isec0+36) {
	      hisZ->GetAxis(2)->SetRange(TMath::Max(ibinSnp-2,firstSnp),TMath::Min(ibinSnp+2,lastSnp));
	    }
	    htemp = hisZ->Projection(2);
	    Double_t      xSnp= htemp->GetMean();
	    delete htemp;
	    THnSparse * hisSnp= hisZ->Projection(2,idim);         
	    //projected histogram according selection 2 - snp
	    
	    //Int_t nbinsTheta    = hisSnp->GetAxis(1)->GetNbins();
	    Int_t firstTheta    = hisSnp->GetAxis(1)->GetFirst();
	    Int_t lastTheta     = hisSnp->GetAxis(1)->GetLast();
	    //
	    for (Int_t ibinTheta=firstTheta; ibinTheta<=lastTheta; ibinTheta+=2){  // axis1 theta
	      
	      
	      hisSnp->GetAxis(1)->SetRange(TMath::Max(ibinTheta-2,firstTheta),TMath::Min(ibinTheta+2,lastTheta));
	      if (isec1!=isec0+36) {
		 hisSnp->GetAxis(1)->SetRange(TMath::Max(ibinTheta-3,firstTheta),TMath::Min(ibinTheta+3,lastTheta));		 
	      }
	      htemp = hisSnp->Projection(1);	      
	      Double_t xTheta=htemp->GetMean();
	      delete htemp;
	      hisDelta = hisSnp->Projection(0);
	      //
	      Double_t entries = hisDelta->GetEntries();
	      Double_t mean=0, rms=0;
	      if (entries>kMinEntries){
		mean    = hisDelta->GetMean(); 
		rms = hisDelta->GetRMS(); 
	      }
	      Double_t sector = 9.*xPhi/TMath::Pi();
	      if (sector<0) sector+=18;
	      Double_t dsec = sector-Int_t(sector)-0.5;
	      Int_t dtype=1;  // TPC alignment type
	      (*pcstream)<<hname<<
		"run="<<run<<
		"bz="<<bz<<             // magnetic field
		"ptype="<<type<<         // parameter type
		"dtype="<<dtype<<         // parameter type
		"isec0="<<isec0<<       // sector 0 
		"isec1="<<isec1<<       // sector 1		
		"sector="<<sector<<     // sector as float
		"dsec="<<dsec<<         // delta sector
		//
		"theta="<<xTheta<<      // theta
		"phi="<<xPhi<<          // phi (alpha)	      
		"z="<<xZ<<              // z
		"snp="<<xSnp<<          // snp
		//
		"entries="<<entries<<  // entries in bin
		"mean="<<mean<<        // mean
		"rms="<<rms<<          // rms 
		"refX="<<refX<<        // track matching reference plane
		"\n";
	      delete hisDelta;
	      printf("%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n",isec0, isec1, xPhi,xZ,xSnp, xTheta, entries,mean);
	      //
	    }//ibinTheta
	    delete hisSnp;
	  } //ibinSnp
	  delete hisZ;
	}//ibinZ
	delete hisPhi;
      }//ibinPhi
      delete hisSector1;      
    }//isec1
    delete hisSector0;
  }//isec0
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
  printf("bz: %f\n",bz);
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
  if (!fgVisualCorrection) fgVisualCorrection=new TObjArray(10000);
  if (position>=fgVisualCorrection->GetEntriesFast())
    fgVisualCorrection->Expand((position+10)*2);
  fgVisualCorrection->AddAt(corr, position);
}



Double_t AliTPCCorrection::GetCorrSector(Double_t sector, Double_t r, Double_t kZ, Int_t axisType, Int_t corrType){
  //
  // calculate the correction at given position - check the geffCorr
  //
  // corrType return values
  // 0 - delta R
  // 1 - delta RPhi
  // 2 - delta Z
  // 3 - delta RPHI
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
  if (axisType==3) return (TMath::Cos(phi)*(distPoint[0]-gx)+ TMath::Cos(phi)*(distPoint[1]-gy));
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





void AliTPCCorrection::MakeLaserDistortionTree(TTree* tree, TObjArray */*corrArray*/, Int_t /*itype*/){
  //
  // Make a laser fit tree for global minimization
  //  
  AliTPCcalibDB*  calib=AliTPCcalibDB::Instance();  
  AliTPCCorrection * correction = calib->GetTPCComposedCorrection();  
  if (!correction) correction = calib->GetTPCComposedCorrection(AliTrackerBase::GetBz());  
  correction->AddVisualCorrection(correction,0);  //register correction

  //  AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform() ;
  //AliTPCParam     *param     = AliTPCcalibDB::Instance()->GetParameters();
  //
  const Double_t cutErrY=0.05;
  const Double_t kSigmaCut=4;
  //  const Double_t cutErrZ=0.03;
  const Double_t kEpsilon=0.00000001;
  //  const Double_t kMaxDist=1.;  // max distance - space correction
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
  TTreeSRedirector *pcstream= new TTreeSRedirector("distortionLaser_0.root");
  Double_t bz=AliTrackerBase::GetBz();
  // 
  //  Double_t globalXYZ[3];
  //Double_t globalXYZCorr[3];
  for (Int_t ientry=0; ientry<entries; ientry++){
    tree->GetEntry(ientry);
    if (!ltr->GetVecGX()){
      ltr->UpdatePoints();
    }
    //
    TVectorD fit10(5);
    TVectorD fit5(5);
    printf("Entry\t%d\n",ientry);
    for (Int_t irow0=0; irow0<158; irow0+=1){
      //       
      TLinearFitter fitter10(4,"hyp3");
      TLinearFitter fitter5(2,"hyp1");
      Int_t sector= (Int_t)(*ltr->GetVecSec())[irow0];
      if (sector<0) continue;
      //if (TMath::Abs(vecdY->GetMatrixArray()[irow0])<kEpsilon) continue;

      Double_t refX= (*ltr->GetVecLX())[irow0];
      Int_t firstRow1 = TMath::Max(irow0-10,0);
      Int_t lastRow1  = TMath::Min(irow0+10,158);
      Double_t padWidth=(irow0<64)?0.4:0.6;
      // make long range fit
      for (Int_t irow1=firstRow1; irow1<=lastRow1; irow1++){
	if (TMath::Abs((*ltr->GetVecSec())[irow1]-sector)>kEpsilon) continue;
	if (veceY->GetMatrixArray()[irow1]>cutErrY) continue;
	if (TMath::Abs(vecdY->GetMatrixArray()[irow1])<kEpsilon) continue;
	Double_t idealX= (*ltr->GetVecLX())[irow1];
	Double_t idealY= (*ltr->GetVecLY())[irow1];
	//	Double_t idealZ= (*ltr->GetVecLZ())[irow1];
	Double_t gx= (*ltr->GetVecGX())[irow1];
	Double_t gy= (*ltr->GetVecGY())[irow1];
	Double_t gz= (*ltr->GetVecGZ())[irow1];
	Double_t measY=(*vecdY)[irow1]+idealY;
	Double_t deltaR = GetCorrXYZ(gx, gy, gz, 0,0);
	// deltaR = R distorted -R ideal
	Double_t xxx[4]={idealX+deltaR-refX,TMath::Cos(idealY/padWidth), TMath::Sin(idealY/padWidth)};
	fitter10.AddPoint(xxx,measY,1);
      }
      Bool_t isOK=kTRUE;
      Double_t rms10=0;//TMath::Sqrt(fitter10.GetChisquare()/(fitter10.GetNpoints()-4));
      Double_t mean10  =0;//   fitter10.GetParameter(0);
      Double_t slope10  =0;//   fitter10.GetParameter(0);
      Double_t cosPart10  = 0;//  fitter10.GetParameter(2);
      Double_t sinPart10   =0;//  fitter10.GetParameter(3); 

      if (fitter10.GetNpoints()>10){
	fitter10.Eval();
	rms10=TMath::Sqrt(fitter10.GetChisquare()/(fitter10.GetNpoints()-4));
	mean10      =   fitter10.GetParameter(0);
	slope10     =   fitter10.GetParameter(1);
	cosPart10   =   fitter10.GetParameter(2);
	sinPart10   =  fitter10.GetParameter(3); 
	//
	// make short range fit
	//
	for (Int_t irow1=firstRow1+5; irow1<=lastRow1-5; irow1++){
	  if (TMath::Abs((*ltr->GetVecSec())[irow1]-sector)>kEpsilon) continue;
	  if (veceY->GetMatrixArray()[irow1]>cutErrY) continue;
	  if (TMath::Abs(vecdY->GetMatrixArray()[irow1])<kEpsilon) continue;
	  Double_t idealX= (*ltr->GetVecLX())[irow1];
	  Double_t idealY= (*ltr->GetVecLY())[irow1];
	  //	  Double_t idealZ= (*ltr->GetVecLZ())[irow1];
	  Double_t gx= (*ltr->GetVecGX())[irow1];
	  Double_t gy= (*ltr->GetVecGY())[irow1];
	  Double_t gz= (*ltr->GetVecGZ())[irow1];
	  Double_t measY=(*vecdY)[irow1]+idealY;
	  Double_t deltaR = GetCorrXYZ(gx, gy, gz, 0,0);
	  // deltaR = R distorted -R ideal 
	  Double_t expY= mean10+slope10*(idealX+deltaR-refX);
	  if (TMath::Abs(measY-expY)>kSigmaCut*rms10) continue;
	  //
	  Double_t corr=cosPart10*TMath::Cos(idealY/padWidth)+sinPart10*TMath::Sin(idealY/padWidth);
	  Double_t xxx[4]={idealX+deltaR-refX,TMath::Cos(idealY/padWidth), TMath::Sin(idealY/padWidth)};
	  fitter5.AddPoint(xxx,measY-corr,1);
	}     
      }else{
	isOK=kFALSE;
      }
      if (fitter5.GetNpoints()<8) isOK=kFALSE;

      Double_t rms5=0;//TMath::Sqrt(fitter5.GetChisquare()/(fitter5.GetNpoints()-4));
      Double_t offset5  =0;//  fitter5.GetParameter(0);
      Double_t slope5   =0;//  fitter5.GetParameter(0); 
      if (isOK){
	fitter5.Eval();
	rms5=TMath::Sqrt(fitter5.GetChisquare()/(fitter5.GetNpoints()-4));
	offset5  =  fitter5.GetParameter(0);
	slope5   =  fitter5.GetParameter(0); 
      }
      //
      Double_t dtype=5;
      Double_t ptype=0;
      Double_t phi   =(*ltr->GetVecPhi())[irow0];
      Double_t theta =ltr->GetTgl();
      Double_t mean=(vecdY)->GetMatrixArray()[irow0];
      Double_t gx=0,gy=0,gz=0;
      Double_t snp = (*ltr->GetVecP2())[irow0];
      Int_t bundle= ltr->GetBundle();
      Int_t id= ltr->GetId();
      //      Double_t rms = err->GetMatrixArray()[irow];
      //
      gx = (*ltr->GetVecGX())[irow0];
      gy = (*ltr->GetVecGY())[irow0];
      gz = (*ltr->GetVecGZ())[irow0];
      Double_t dRrec = GetCorrXYZ(gx, gy, gz, 0,0);
      fitter10.GetParameters(fit10);
      fitter5.GetParameters(fit5);      
      Double_t idealY= (*ltr->GetVecLY())[irow0];
      Double_t measY=(*vecdY)[irow0]+idealY;
      Double_t corr=cosPart10*TMath::Cos(idealY/padWidth)+sinPart10*TMath::Sin(idealY/padWidth);
      if (TMath::Max(rms5,rms10)>0.06) isOK=kFALSE;
      //
      (*pcstream)<<"fitFull"<<  // dumpe also intermediate results
	"bz="<<bz<<         // magnetic filed used
	"dtype="<<dtype<<   // detector match type
	"ptype="<<ptype<<   // parameter type
	"theta="<<theta<<   // theta
	"phi="<<phi<<       // phi 
	"snp="<<snp<<       // snp
	"sector="<<sector<<
	"bundle="<<bundle<<
// 	//	"dsec="<<dsec<<
	"refX="<<refX<<      // reference radius
	"gx="<<gx<<         // global position
	"gy="<<gy<<         // global position
	"gz="<<gz<<         // global position
	"dRrec="<<dRrec<<      // delta Radius in reconstruction
 	"id="<<id<<     //bundle
	"rms10="<<rms10<<
	"rms5="<<rms5<<
	"fit10.="<<&fit10<<
	"fit5.="<<&fit5<<
	"measY="<<measY<<
	"mean="<<mean<<
	"idealY="<<idealY<<
	"corr="<<corr<<
	"isOK="<<isOK<<
	"\n";
    }
  }
  delete pcstream;
}
