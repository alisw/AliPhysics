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

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.54  2006/08/28 10:01:56  kharlov
 * Effective C++ warnings fixed (Timur Pocheptsov)
 *
 * Revision 1.53  2005/12/20 14:28:47  hristov
 * Additional protection
 *
 * Revision 1.52  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
//  RecPoint implementation for PHOS-EMC 
//  An EmcRecPoint is a cluster of digits   
//*--
//*-- Author: Dmitri Peressounko (RRC KI & SUBATECH)


// --- ROOT system ---
#include "TH2.h"
#include "TMath.h" 
#include "TCanvas.h" 
#include "TGraph.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliPHOSLoader.h"
#include "AliGenerator.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#include "AliPHOSEmcRecPoint.h"
 
ClassImp(AliPHOSEmcRecPoint)

//____________________________________________________________________________
AliPHOSEmcRecPoint::AliPHOSEmcRecPoint() : 
  AliPHOSRecPoint(),
  fCoreEnergy(0.), fDispersion(0.),
  fEnergyList(0), fTime(-1.), fNExMax(0),
  fM2x(0.), fM2z(0.), fM3x(0.), fM4z(0.),
  fPhixe(0.), fDebug(0)
{
  // ctor
  fMulDigit   = 0 ;  
  fAmp   = 0. ;   
  fLocPos.SetX(1000000.)  ;      //Local position should be evaluated
}

//____________________________________________________________________________
AliPHOSEmcRecPoint::AliPHOSEmcRecPoint(const char * opt) : 
  AliPHOSRecPoint(opt),
  fCoreEnergy(0.), fDispersion(0.),
  fEnergyList(0), fTime(-1.), fNExMax(0),
  fM2x(0.), fM2z(0.), fM3x(0.), fM4z(0.),
  fPhixe(0.), fDebug(0)
{
  // ctor
  fMulDigit   = 0 ;  
  fAmp   = 0. ;   
  fLocPos.SetX(1000000.)  ;      //Local position should be evaluated
}

//____________________________________________________________________________
AliPHOSEmcRecPoint::AliPHOSEmcRecPoint(const AliPHOSEmcRecPoint & rp) : 
  AliPHOSRecPoint(rp),
  fCoreEnergy(rp.fCoreEnergy), fDispersion(rp.fDispersion),
  fEnergyList(0), fTime(rp.fTime), fNExMax(rp.fNExMax),
  fM2x(rp.fM2x), fM2z(rp.fM2z), fM3x(rp.fM3x), fM4z(rp.fM4z),
  fPhixe(rp.fPhixe), fDebug(rp.fDebug)
{
  // cpy ctor
  fMulDigit   = rp.fMulDigit ;  
  fAmp        = rp.fAmp ;   
  fEnergyList = new Float_t[rp.fMulDigit] ;
  Int_t index ; 
  for(index = 0 ; index < fMulDigit ; index++) 
    fEnergyList[index] = rp.fEnergyList[index] ; 
}

//____________________________________________________________________________
AliPHOSEmcRecPoint::~AliPHOSEmcRecPoint()
{
  // dtor
  if ( fEnergyList )
    delete[] fEnergyList ; 
}

//____________________________________________________________________________
void AliPHOSEmcRecPoint::AddDigit(AliPHOSDigit & digit, Float_t Energy)
{
  // Adds a digit to the RecPoint
  // and accumulates the total amplitude and the multiplicity 
  
  if(fEnergyList == 0)
    fEnergyList =  new Float_t[fMaxDigit]; 

  if ( fMulDigit >= fMaxDigit ) { // increase the size of the lists 
    fMaxDigit*=2 ; 
    Int_t * tempo = new Int_t[fMaxDigit]; 
    Float_t * tempoE =  new Float_t[fMaxDigit];

    Int_t index ;     
    for ( index = 0 ; index < fMulDigit ; index++ ){
      tempo[index]  = fDigitsList[index] ;
      tempoE[index] = fEnergyList[index] ; 
    }
    
    delete [] fDigitsList ; 
    fDigitsList =  new Int_t[fMaxDigit];
 
    delete [] fEnergyList ;
    fEnergyList =  new Float_t[fMaxDigit];

    for ( index = 0 ; index < fMulDigit ; index++ ){
      fDigitsList[index] = tempo[index] ;
      fEnergyList[index] = tempoE[index] ; 
    }
 
    delete [] tempo ;
    delete [] tempoE ; 
  } // if
  
  fDigitsList[fMulDigit]   = digit.GetIndexInList()  ; 
  fEnergyList[fMulDigit]   = Energy ;
  fMulDigit++ ; 
  fAmp += Energy ; 

  EvalPHOSMod(&digit) ;
}

//____________________________________________________________________________
Bool_t AliPHOSEmcRecPoint::AreNeighbours(AliPHOSDigit * digit1, AliPHOSDigit * digit2 ) const
{
  // Tells if (true) or not (false) two digits are neighbors
  
  Bool_t aren = kFALSE ;
  
  AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance() ;

  Int_t relid1[4] ; 
  phosgeom->AbsToRelNumbering(digit1->GetId(), relid1) ; 

  Int_t relid2[4] ; 
  phosgeom->AbsToRelNumbering(digit2->GetId(), relid2) ; 
  
  Int_t rowdiff = TMath::Abs( relid1[2] - relid2[2] ) ;  
  Int_t coldiff = TMath::Abs( relid1[3] - relid2[3] ) ;  

  if (( coldiff <= 1 )  && ( rowdiff <= 1 ) && (coldiff + rowdiff > 0)) 
    aren = kTRUE ;
  
  return aren ;
}

//____________________________________________________________________________
Int_t AliPHOSEmcRecPoint::Compare(const TObject * obj) const
{
  // Compares two RecPoints according to their position in the PHOS modules
  
  const Float_t delta = 1 ; //Width of "Sorting row". If you changibg this 
                      //value (what is senseless) change as vell delta in
                      //AliPHOSTrackSegmentMakerv* and other RecPoints...
  Int_t rv ; 

  AliPHOSEmcRecPoint * clu = (AliPHOSEmcRecPoint *)obj ; 

 
  Int_t phosmod1 = GetPHOSMod() ;
  Int_t phosmod2 = clu->GetPHOSMod() ;

  TVector3 locpos1; 
  GetLocalPosition(locpos1) ;
  TVector3 locpos2;  
  clu->GetLocalPosition(locpos2) ;  

  if(phosmod1 == phosmod2 ) {
    Int_t rowdif = (Int_t)TMath::Ceil(locpos1.X()/delta)-(Int_t)TMath::Ceil(locpos2.X()/delta) ;
    if (rowdif> 0) 
      rv = 1 ;
     else if(rowdif < 0) 
       rv = -1 ;
    else if(locpos1.Z()>locpos2.Z()) 
      rv = -1 ;
    else 
      rv = 1 ; 
     }

  else {
    if(phosmod1 < phosmod2 ) 
      rv = -1 ;
    else 
      rv = 1 ;
  }

  return rv ; 
}
//______________________________________________________________________________
void AliPHOSEmcRecPoint::ExecuteEvent(Int_t event, Int_t, Int_t) /*const*/
{
  
  // Execute action corresponding to one event
  //  This member function is called when a AliPHOSRecPoint is clicked with the locator
  //
  //  If Left button is clicked on AliPHOSRecPoint, the digits are switched on    
  //  and switched off when the mouse button is released.
  
  
  AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance();
  
  static TGraph *  digitgraph = 0 ;
  
  if (!gPad->IsEditable()) return;
  
  TH2F * histo = 0 ;
  TCanvas * histocanvas ; 

  
  //try to get run loader from default event folder
  AliRunLoader* rn = AliRunLoader::GetRunLoader(AliConfig::GetDefaultEventFolderName());
  if (rn == 0x0) 
    {
      AliError(Form("Cannot find Run Loader in Default Event Folder"));
      return;
    }
  AliPHOSLoader* gime = dynamic_cast<AliPHOSLoader*>(rn->GetLoader("PHOSLoader"));
  if (gime == 0x0) 
    {
      AliError(Form("Cannot find PHOS Loader from Run Loader"));
      return;
    }
  
  
  const TClonesArray * digits = gime->Digits() ;
  
  switch (event) {
    
  case kButton1Down: {
    AliPHOSDigit * digit ;
    Int_t iDigit;
    Int_t relid[4] ;
    
    const Int_t kMulDigit = AliPHOSEmcRecPoint::GetDigitsMultiplicity() ; 
    Float_t * xi = new Float_t[kMulDigit] ; 
    Float_t * zi = new Float_t[kMulDigit] ; 
    
    // create the histogram for the single cluster 
    // 1. gets histogram boundaries
    Float_t ximax = -999. ; 
    Float_t zimax = -999. ; 
    Float_t ximin = 999. ; 
    Float_t zimin = 999. ;
    
    for(iDigit=0; iDigit<kMulDigit; iDigit++) {
      digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit])  ;
      phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
      phosgeom->RelPosInModule(relid, xi[iDigit], zi[iDigit]);
      if ( xi[iDigit] > ximax )
	ximax = xi[iDigit] ; 
      if ( xi[iDigit] < ximin )
	ximin = xi[iDigit] ; 
      if ( zi[iDigit] > zimax )
	zimax = zi[iDigit] ; 
      if ( zi[iDigit] < zimin )
	zimin = zi[iDigit] ;     
    }
    ximax += phosgeom->GetCrystalSize(0) / 2. ;
    zimax += phosgeom->GetCrystalSize(2) / 2. ;
    ximin -= phosgeom->GetCrystalSize(0) / 2. ;
    zimin -= phosgeom->GetCrystalSize(2) / 2. ;
    Int_t xdim = (int)( (ximax - ximin ) / phosgeom->GetCrystalSize(0) + 0.5  ) ; 
    Int_t zdim = (int)( (zimax - zimin ) / phosgeom->GetCrystalSize(2) + 0.5 ) ;
    
    // 2. gets the histogram title
    
    Text_t title[100] ; 
    sprintf(title,"Energy=%1.2f GeV ; Digits ; %d ", GetEnergy(), GetDigitsMultiplicity()) ;
    
    if (!histo) {
      delete histo ; 
      histo = 0 ; 
    }
    histo = new TH2F("cluster3D", title,  xdim, ximin, ximax, zdim, zimin, zimax)  ;
    
    Float_t x, z ; 
    for(iDigit=0; iDigit<kMulDigit; iDigit++) {
      digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit])  ;
      phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
      phosgeom->RelPosInModule(relid, x, z);
      histo->Fill(x, z, fEnergyList[iDigit] ) ;
    }
    
    if (!digitgraph) {
      digitgraph = new TGraph(kMulDigit,xi,zi);
      digitgraph-> SetMarkerStyle(5) ; 
      digitgraph-> SetMarkerSize(1.) ;
      digitgraph-> SetMarkerColor(1) ;
      digitgraph-> Paint("P") ;
    }
    
    //    Print() ;
    histocanvas = new TCanvas("cluster", "a single cluster", 600, 500) ; 
    histocanvas->Draw() ; 
    histo->Draw("lego1") ; 
    
    delete[] xi ; 
    delete[] zi ; 
    
    break;
  }
  
  case kButton1Up: 
    if (digitgraph) {
      delete digitgraph  ;
      digitgraph = 0 ;
    }
    break;
  
   }
}

//____________________________________________________________________________
void  AliPHOSEmcRecPoint::EvalDispersion(Float_t logWeight,TClonesArray * digits)
{
  // Calculates the dispersion of the shower at the origine of the RecPoint

  Float_t d    = 0. ;
  Float_t wtot = 0. ;

  Float_t x = 0.;
  Float_t z = 0.;

  AliPHOSDigit * digit ;
 
  AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance();

 // Calculates the center of gravity in the local PHOS-module coordinates 
  
  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit]) ;
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    if (fAmp>0 && fEnergyList[iDigit]>0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ) ) ;
      x    += xi * w ;
      z    += zi * w ;
      wtot += w ;
    }
    else
      AliError(Form("Wrong energy %f and/or amplitude %f\n", fEnergyList[iDigit], fAmp));
  }
  if (wtot>0) {
    x /= wtot ;
    z /= wtot ;
  }
  else
    AliError(Form("Wrong weight %f\n", wtot));


// Calculates the dispersion in coordinates 
  wtot = 0.;
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit])  ;
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    if (fAmp>0 && fEnergyList[iDigit]>0) {
      Float_t w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
      d += w*((xi-x)*(xi-x) + (zi-z)*(zi-z) ) ; 
      wtot+=w ;
    }
    else
      AliError(Form("Wrong energy %f and/or amplitude %f\n", fEnergyList[iDigit], fAmp));
  }
  

  if (wtot>0) {
    d /= wtot ;
  }
  else
    AliError(Form("Wrong weight %f\n", wtot));

  fDispersion = 0;
  if (d>=0)
    fDispersion = TMath::Sqrt(d) ;
 
}
//______________________________________________________________________________
void AliPHOSEmcRecPoint::EvalCoreEnergy(Float_t logWeight, TClonesArray * digits)
{
  // This function calculates energy in the core, 
  // i.e. within a radius rad = 3cm around the center. Beyond this radius
  // in accordance with shower profile the energy deposition 
  // should be less than 2%

  Float_t coreRadius = 3 ;

  Float_t x = 0 ;
  Float_t z = 0 ;

  AliPHOSDigit * digit ;

  AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance();

  Int_t iDigit;

// Calculates the center of gravity in the local PHOS-module coordinates 
  Float_t wtot = 0;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit]) ;
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    if (fAmp>0 && fEnergyList[iDigit]>0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ) ) ;
      x    += xi * w ;
      z    += zi * w ;
      wtot += w ;
    }
    else
      AliError(Form("Wrong energy %f and/or amplitude %f\n", fEnergyList[iDigit], fAmp));
  }
  if (wtot>0) {
    x /= wtot ;
    z /= wtot ;
  }
  else
    AliError(Form("Wrong weight %f\n", wtot));


  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) ( digits->At(fDigitsList[iDigit]) ) ;
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);    
    Float_t distance = TMath::Sqrt((xi-x)*(xi-x)+(zi-z)*(zi-z)) ;
    if(distance < coreRadius)
      fCoreEnergy += fEnergyList[iDigit] ;
  }

}

//____________________________________________________________________________
void  AliPHOSEmcRecPoint::EvalElipsAxis(Float_t logWeight,TClonesArray * digits)
{
  // Calculates the axis of the shower ellipsoid

  Double_t wtot = 0. ;
  Double_t x    = 0.;
  Double_t z    = 0.;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;

  AliPHOSDigit * digit ;

  AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance();

  Int_t iDigit;


  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit])  ;
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    if (fAmp>0 && fEnergyList[iDigit]>0) {
      Double_t w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
      dxx  += w * xi * xi ;
      x    += w * xi ;
      dzz  += w * zi * zi ;
      z    += w * zi ; 
      dxz  += w * xi * zi ; 
      wtot += w ;
    }
    else
      AliError(Form("Wrong energy %f and/or amplitude %f\n", fEnergyList[iDigit], fAmp));
  }
  if (wtot>0) {
    dxx /= wtot ;
    x   /= wtot ;
    dxx -= x * x ;
    dzz /= wtot ;
    z   /= wtot ;
    dzz -= z * z ;
    dxz /= wtot ;
    dxz -= x * z ;

//   //Apply correction due to non-perpendicular incidence
//   Double_t CosX ;
//   Double_t CosZ ;
//   AliPHOSGetter * gime = AliPHOSGetter::Instance() ; 
//   AliPHOSGeometry * phosgeom =  (AliPHOSGeometry*)gime->PHOSGeometry();
  //   Double_t DistanceToIP= (Double_t ) phosgeom->GetIPtoCrystalSurface() ;
  
//   CosX = DistanceToIP/TMath::Sqrt(DistanceToIP*DistanceToIP+x*x) ;
//   CosZ = DistanceToIP/TMath::Sqrt(DistanceToIP*DistanceToIP+z*z) ;

//   dxx = dxx/(CosX*CosX) ;
//   dzz = dzz/(CosZ*CosZ) ;
//   dxz = dxz/(CosX*CosZ) ;


    fLambda[0] =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    if(fLambda[0] > 0)
      fLambda[0] = TMath::Sqrt(fLambda[0]) ;

    fLambda[1] =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    if(fLambda[1] > 0) //To avoid exception if numerical errors lead to negative lambda.
      fLambda[1] = TMath::Sqrt(fLambda[1]) ;
    else
      fLambda[1]= 0. ;
  }
  else {
    AliError(Form("Wrong weight %f\n", wtot));
    fLambda[0]=fLambda[1]=0.;
  }
}

//____________________________________________________________________________
void  AliPHOSEmcRecPoint::EvalMoments(Float_t logWeight,TClonesArray * digits)
{
  // Calculate the shower moments in the eigen reference system
  // M2x, M2z, M3x, M4z
  // Calculate the angle between the shower position vector and the eigen vector

  Double_t wtot = 0. ;
  Double_t x    = 0.;
  Double_t z    = 0.;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;
  Double_t lambda0=0, lambda1=0;

  AliPHOSDigit * digit ;

  AliPHOSGeometry * phosgeom = AliPHOSGeometry::GetInstance() ;

  Int_t iDigit;

  // 1) Find covariance matrix elements:
  //    || dxx dxz ||
  //    || dxz dzz ||

  Float_t xi ;
  Float_t zi ;
  Int_t relid[4] ;
  Double_t w;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit])  ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    if (fAmp>0 && fEnergyList[iDigit]>0) {
      w     = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
      x    += w * xi ;
      z    += w * zi ; 
      dxx  += w * xi * xi ;
      dzz  += w * zi * zi ;
      dxz  += w * xi * zi ; 
      wtot += w ;
    }
    else
      AliError(Form("Wrong energy %f and/or amplitude %f\n", fEnergyList[iDigit], fAmp));
    
  }
  if (wtot>0) {
    x   /= wtot ;
    z   /= wtot ;
    dxx /= wtot ;
    dzz /= wtot ;
    dxz /= wtot ;
    dxx -= x * x ;
    dzz -= z * z ;
    dxz -= x * z ;

  // 2) Find covariance matrix eigen values lambda0 and lambda1

    lambda0 =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    lambda1 =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
  }
  else {
    AliError(Form("Wrong weight %f\n", wtot));
    lambda0=lambda1=0.;
  }

  // 3) Find covariance matrix eigen vectors e0 and e1

  TVector2 e0,e1;
  if (dxz != 0)
    e0.Set(1.,(lambda0-dxx)/dxz);
  else
    e0.Set(0.,1.);

  e0 = e0.Unit();
  e1.Set(-e0.Y(),e0.X());

  // 4) Rotate cluster tensor from (x,z) to (e0,e1) system
  //    and calculate moments M3x and M4z

  Float_t cosPhi = e0.X();
  Float_t sinPhi = e0.Y();

  Float_t xiPHOS ;
  Float_t ziPHOS ;
  Double_t dx3, dz3, dz4;
  wtot = 0.;
  x    = 0.;
  z    = 0.;
  dxx  = 0.;
  dzz  = 0.;
  dxz  = 0.;
  dx3  = 0.;
  dz3  = 0.;
  dz4  = 0.;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit])  ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xiPHOS, ziPHOS);
    xi    = xiPHOS*cosPhi + ziPHOS*sinPhi;
    zi    = ziPHOS*cosPhi - xiPHOS*sinPhi;
    if (fAmp>0 && fEnergyList[iDigit]>0) {
      w     = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
      x    += w * xi ;
      z    += w * zi ; 
      dxx  += w * xi * xi ;
      dzz  += w * zi * zi ;
      dxz  += w * xi * zi ; 
      dx3  += w * xi * xi * xi;
      dz3  += w * zi * zi * zi ;
      dz4  += w * zi * zi * zi * zi ;
      wtot += w ;
    }
    else
      AliError(Form("Wrong energy %f and/or amplitude %f\n", fEnergyList[iDigit], fAmp));
  }
  if (wtot>0) {
    x   /= wtot ;
    z   /= wtot ;
    dxx /= wtot ;
    dzz /= wtot ;
    dxz /= wtot ;
    dx3 /= wtot ;
    dz3 /= wtot ;
    dz4 /= wtot ;
    dx3 += -3*dxx*x + 2*x*x*x;
    dz4 += -4*dz3*z + 6*dzz*z*z -3*z*z*z*z;
    dxx -= x * x ;
    dzz -= z * z ;
    dxz -= x * z ;
  }
  else
    AliError(Form("Wrong weight %f\n", wtot));

  // 5) Find an angle between cluster center vector and eigen vector e0

  Float_t phi = 0;
  if ( (x*x+z*z) > 0 ) 
    phi = TMath::ACos ((x*e0.X() + z*e0.Y()) / sqrt(x*x + z*z));

  fM2x   = lambda0;
  fM2z   = lambda1;
  fM3x   = dx3;
  fM4z   = dz4;
  fPhixe = phi;

}
//____________________________________________________________________________
void AliPHOSEmcRecPoint::EvalAll(Float_t logWeight, TClonesArray * digits )
{
  // Evaluates all shower parameters
  EvalLocalPosition(logWeight, digits) ;
  EvalElipsAxis(logWeight, digits) ;
  EvalMoments(logWeight, digits) ;
  EvalDispersion(logWeight, digits) ;
  EvalCoreEnergy(logWeight, digits);
  EvalTime(digits) ;
  AliPHOSRecPoint::EvalAll(digits) ;
}
//____________________________________________________________________________
void AliPHOSEmcRecPoint::EvalLocalPosition(Float_t logWeight, TClonesArray * digits)
{
  // Calculates the center of gravity in the local PHOS-module coordinates 
  Float_t wtot = 0. ;
 
  Int_t relid[4] ;

  Float_t x = 0. ;
  Float_t z = 0. ;
  
  AliPHOSDigit * digit ;

  AliPHOSGeometry * phosgeom =  AliPHOSGeometry::GetInstance() ;

  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliPHOSDigit *) digits->At(fDigitsList[iDigit]) ;

    Float_t xi ;
    Float_t zi ;
    phosgeom->AbsToRelNumbering(digit->GetId(), relid) ;
    phosgeom->RelPosInModule(relid, xi, zi);
    if (fAmp>0 && fEnergyList[iDigit]>0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ) ) ;
      x    += xi * w ;
      z    += zi * w ;
      wtot += w ;
    }
    else
      AliError(Form("Wrong energy %f and/or amplitude %f\n", fEnergyList[iDigit], fAmp));
  }
  if (wtot>0) {
    x /= wtot ;
    z /= wtot ;
  }
  else
    AliError(Form("Wrong weight %f\n", wtot));

  // Correction for the depth of the shower starting point (TDR p 127)  
  Float_t para = 0.925 ; 
  Float_t parb = 6.52 ; 

  Float_t xo,yo,zo ; //Coordinates of the origin
  //We should check all 3 possibilities to avoid seg.v.
  if(gAlice && gAlice->GetMCApp() && gAlice->Generator())
    gAlice->Generator()->GetOrigin(xo,yo,zo) ;
  else{
    xo=yo=zo=0.;
  }
  Float_t phi = phosgeom->GetPHOSAngle(relid[0]) ;

  //Transform to the local ref.frame
  Float_t xoL,yoL ;
  xoL = xo*TMath::Cos(phi)-yo*TMath::Sin(phi) ;
  yoL = xo*TMath::Sin(phi)+yo*TMath::Cos(phi) ;
  
  Float_t radius = phosgeom->GetIPtoCrystalSurface()-yoL;
 
  Float_t incidencephi = TMath::ATan((x-xoL ) / radius) ; 
  Float_t incidencetheta = TMath::ATan((z-zo) / radius) ;
 
  Float_t depthx = 0.; 
  Float_t depthz = 0.;
  if (fAmp>0) {
    depthx = ( para * TMath::Log(fAmp) + parb ) * TMath::Sin(incidencephi) ; 
    depthz = ( para * TMath::Log(fAmp) + parb ) * TMath::Sin(incidencetheta) ; 
  }
  else 
    AliError(Form("Wrong amplitude %f\n", fAmp));

  fLocPos.SetX(x - depthx)  ;
  fLocPos.SetY(0.) ;
  fLocPos.SetZ(z - depthz)  ;

  fLocPosM = 0 ;
}

//____________________________________________________________________________
Float_t AliPHOSEmcRecPoint::GetMaximalEnergy(void) const
{
  // Finds the maximum energy in the cluster
  
  Float_t menergy = 0. ;

  Int_t iDigit;

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
 
    if(fEnergyList[iDigit] > menergy) 
      menergy = fEnergyList[iDigit] ;
  }
  return menergy ;
}

//____________________________________________________________________________
Int_t AliPHOSEmcRecPoint::GetMultiplicityAtLevel(Float_t H) const
{
  // Calculates the multiplicity of digits with energy larger than H*energy 
  
  Int_t multipl   = 0 ;
  Int_t iDigit ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {

    if(fEnergyList[iDigit] > H * fAmp) 
      multipl++ ;
  }
  return multipl ;
}

//____________________________________________________________________________
Int_t  AliPHOSEmcRecPoint::GetNumberOfLocalMax( AliPHOSDigit **  maxAt, Float_t * maxAtEnergy,
					       Float_t locMaxCut,TClonesArray * digits) const
{ 
  // Calculates the number of local maxima in the cluster using fLocalMaxCut as the minimum
  // energy difference between two local maxima

  AliPHOSDigit * digit ;
  AliPHOSDigit * digitN ;
  

  Int_t iDigitN ;
  Int_t iDigit ;

  for(iDigit = 0; iDigit < fMulDigit; iDigit++)
    maxAt[iDigit] = (AliPHOSDigit*) digits->At(fDigitsList[iDigit])  ;

  
  for(iDigit = 0 ; iDigit < fMulDigit; iDigit++) {   
    if(maxAt[iDigit]) {
      digit = maxAt[iDigit] ;
          
      for(iDigitN = 0; iDigitN < fMulDigit; iDigitN++) {	
	if(iDigit == iDigitN)
	  continue ;
	
	digitN = (AliPHOSDigit *) digits->At(fDigitsList[iDigitN]) ; 
	
	if ( AreNeighbours(digit, digitN) ) {
	  if (fEnergyList[iDigit] > fEnergyList[iDigitN] ) {    
	    maxAt[iDigitN] = 0 ;
	    // but may be digit too is not local max ?
	    if(fEnergyList[iDigit] < fEnergyList[iDigitN] + locMaxCut) 
	      maxAt[iDigit] = 0 ;
	  }
	  else {
	    maxAt[iDigit] = 0 ;
	    // but may be digitN too is not local max ?
	    if(fEnergyList[iDigit] > fEnergyList[iDigitN] - locMaxCut) 
	      maxAt[iDigitN] = 0 ; 
	  } 
	} // if Areneighbours
      } // while digitN
    } // slot not empty
  } // while digit
  
  iDigitN = 0 ;
  for(iDigit = 0; iDigit < fMulDigit; iDigit++) { 
    if(maxAt[iDigit]){
      maxAt[iDigitN] = maxAt[iDigit] ;
      maxAtEnergy[iDigitN] = fEnergyList[iDigit] ;
      iDigitN++ ; 
    }
  }
  return iDigitN ;
}
//____________________________________________________________________________
void AliPHOSEmcRecPoint::EvalTime(TClonesArray * digits)
{
  // Define a rec.point time as a time in the cell with the maximum energy

  Float_t maxE = 0;
  Int_t maxAt = 0;
  for(Int_t idig=0; idig < fMulDigit; idig++){
    if(fEnergyList[idig] > maxE){
      maxE = fEnergyList[idig] ;
      maxAt = idig;
    }
  }
  fTime = ((AliPHOSDigit*) digits->At(fDigitsList[maxAt]))->GetTime() ;
  
}
//____________________________________________________________________________
void AliPHOSEmcRecPoint::Purify(Float_t threshold){
  //Removes digits below threshold

  Int_t * tempo = new Int_t[fMaxDigit]; 
  Float_t * tempoE =  new Float_t[fMaxDigit];

  Int_t mult = 0 ;
  for(Int_t iDigit=0;iDigit< fMulDigit ;iDigit++){
    if(fEnergyList[iDigit] > threshold){
      tempo[mult] = fDigitsList[iDigit] ;
      tempoE[mult] = fEnergyList[iDigit] ;
      mult++ ;
    }
  }
  
  fMulDigit = mult ;
  delete [] fDigitsList ;
  delete [] fEnergyList ;
  fDigitsList = new Int_t[fMulDigit];
  fEnergyList = new Float_t[fMulDigit];

  fAmp = 0.; //Recalculate total energy
  for(Int_t iDigit=0;iDigit< fMulDigit ;iDigit++){
     fDigitsList[iDigit] = tempo[iDigit];
     fEnergyList[iDigit] = tempoE[iDigit] ;
     fAmp+=tempoE[iDigit];
  }
      
  delete [] tempo ;
  delete [] tempoE ;

}
//____________________________________________________________________________
void AliPHOSEmcRecPoint::Print(Option_t *) const
{
  // Print the list of digits belonging to the cluster
  
  TString message ; 
  message  = "AliPHOSEmcRecPoint:\n" ;
  message += "Digit multiplicity = %d" ;
  message += ", cluster Energy = %f" ; 
  message += ", number of primaries = %d" ; 
  message += ", rec.point index = %d \n" ; 
  printf(message.Data(), fMulDigit, fAmp, fMulTrack,GetIndexInList() ) ;  

  Int_t iDigit;
  printf(" digits ids = ") ; 
  for(iDigit=0; iDigit<fMulDigit; iDigit++)
    printf(" %d ", fDigitsList[iDigit] ) ;  
  
  printf("\n digit energies = ") ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) 
    printf(" %f ", fEnergyList[iDigit] ) ;

  printf("\n digit primaries  ") ;
  if (fMulTrack<1) printf("... no primaries");
  for(iDigit = 0;iDigit < fMulTrack; iDigit++)
    printf(" %d ", fTracksList[iDigit]) ;
  printf("\n") ; 	

}
 
  
