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
//  Reconstructed Points for the EMCAL
//  A RecPoint is a cluster of digits
//*-- Author: Yves Schutz (SUBATECH)
//*-- Author: Dmitri Peressounko (RRC KI & SUBATECH)
//*-- Author: Heather Gray (LBL) merged AliEMCALRecPoint and AliEMCALTowerRecPoint 02/04

// --- ROOT system ---
#include "TPad.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TMath.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliGenerator.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRecPoint.h"
#include "AliEMCALGetter.h"

ClassImp(AliEMCALRecPoint)

//____________________________________________________________________________
AliEMCALRecPoint::AliEMCALRecPoint()
  : AliRecPoint()
{
  // ctor
  fMaxTrack = 0 ;
  fMulDigit   = 0 ;  
  fMaxParent = 0;
  fMulParent = 0;
  fAmp   = 0. ;   
  fCoreEnergy = 0 ; 
  fEnergyList = 0 ;
  fParentsList = 0;
  fTime = 0. ;
  fLocPos.SetX(0.)  ;      //Local position should be evaluated
  fCoreRadius = 10;        //HG Check this
}

//____________________________________________________________________________
AliEMCALRecPoint::AliEMCALRecPoint(const char * opt) : AliRecPoint(opt)
{
  // ctor
  fMaxTrack = 1000 ;
  fMaxParent = 1000;
  fMulDigit   = 0 ; 
  fMulParent = 0; 
  fAmp   = 0. ;   
  fCoreEnergy = 0 ; 
  fEnergyList = 0 ;
  fParentsList = new Int_t[fMaxParent];
  fTime = -1. ;
  fLocPos.SetX(1000000.)  ;      //Local position should be evaluated
  fCoreRadius = 10;        //HG Check this
}
//____________________________________________________________________________
AliEMCALRecPoint::~AliEMCALRecPoint()
{
  // dtor
  if ( fEnergyList )
    delete[] fEnergyList ; 
   if ( fParentsList)
    delete[] fParentsList;
}

//____________________________________________________________________________
void AliEMCALRecPoint::AddDigit(AliEMCALDigit & digit, Float_t Energy)
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

}
//____________________________________________________________________________
Bool_t AliEMCALRecPoint::AreNeighbours(AliEMCALDigit * digit1, AliEMCALDigit * digit2 ) const
{
  // Tells if (true) or not (false) two digits are neighbours
  // A neighbour is defined as being two digits which share a corner
  
  Bool_t areNeighbours = kFALSE ;
  
  AliEMCALGeometry * geom =  (AliEMCALGetter::Instance())->EMCALGeometry();

  Int_t relid1[2] ; 
    //copied for shish-kebab geometry, ieta,iphi is cast as float as eta,phi conversion
    // for this geometry does not exist
    int nSupMod=0, nTower=0, nIphi=0, nIeta=0;
    int iphi=0, ieta=0;
       geom->GetCellIndex(digit1->GetId(), nSupMod,nTower,nIphi,nIeta);
       geom->GetCellPhiEtaIndexInSModule(nSupMod,nTower,nIphi,nIeta, iphi,ieta);
       relid1[0]=ieta;
       relid1[1]=iphi;
//  geom->AbsToRelNumbering(digit1->GetId(), relid1) ; 

  Int_t relid2[2] ; 
    //copied for shish-kebab geometry, ieta,iphi is cast as float as eta,phi conversion
    // for this geometry does not exist
    int nSupMod1=0, nTower1=0, nIphi1=0, nIeta1=0;
    int iphi1=0, ieta1=0;
       geom->GetCellIndex(digit2->GetId(), nSupMod1,nTower1,nIphi1,nIeta1);
       geom->GetCellPhiEtaIndexInSModule(nSupMod1,nTower1,nIphi1,nIeta1, iphi1,ieta1);
       relid2[0]=ieta1;
       relid2[1]=iphi1;
//  geom->AbsToRelNumbering(digit2->GetId(), relid2) ; 
  
  Int_t rowdiff = TMath::Abs( relid1[0] - relid2[0] ) ;  
  Int_t coldiff = TMath::Abs( relid1[1] - relid2[1] ) ;  

  if (( coldiff <= 1 )  && ( rowdiff <= 1 ) && (coldiff + rowdiff > 0)) 
    areNeighbours = kTRUE ;
  
  return areNeighbours;
}

//____________________________________________________________________________
Int_t AliEMCALRecPoint::Compare(const TObject * obj) const
{
  // Compares two RecPoints according to their position in the EMCAL modules

  Float_t delta = 1 ; //Width of "Sorting row". If you change this 
                      //value (what is senseless) change as well delta in
                      //AliEMCALTrackSegmentMakerv* and other RecPoints...
  Int_t rv ; 

  AliEMCALRecPoint * clu = (AliEMCALRecPoint *)obj ; 

  TVector3 locpos1; 
  GetLocalPosition(locpos1);
  TVector3 locpos2;  
  clu->GetLocalPosition(locpos2);  

  Int_t rowdif = (Int_t)TMath::Ceil(locpos1.X()/delta)-(Int_t)TMath::Ceil(locpos2.X()/delta) ;
  if (rowdif> 0) 
    rv = 1 ;
  else if(rowdif < 0) 
    rv = -1 ;
  else if(locpos1.Y()>locpos2.Y()) 
    rv = -1 ;
  else 
    rv = 1 ; 

  return rv ; 
}

//____________________________________________________________________________
Int_t AliEMCALRecPoint::DistancetoPrimitive(Int_t px, Int_t py)
{
  // Compute distance from point px,py to  a AliEMCALRecPoint considered as a Tmarker
  // Compute the closest distance of approach from point px,py to this marker.
  // The distance is computed in pixels units.
  // HG Still need to update -> Not sure what this should achieve

  TVector3 pos(0.,0.,0.) ;
  GetLocalPosition(pos) ;
  Float_t x =  pos.X() ;
  Float_t y =  pos.Y() ;
  const Int_t kMaxDiff = 10;
  Int_t pxm  = gPad->XtoAbsPixel(x);
  Int_t pym  = gPad->YtoAbsPixel(y);
  Int_t dist = (px-pxm)*(px-pxm) + (py-pym)*(py-pym);
  
  if (dist > kMaxDiff) return 9999;
  return dist;
}

//___________________________________________________________________________
 void AliEMCALRecPoint::Draw(Option_t *option)
 {
   // Draw this AliEMCALRecPoint with its current attributes
   
   AppendPad(option);
 }

//______________________________________________________________________________
void AliEMCALRecPoint::ExecuteEvent(Int_t /*event*/, Int_t, Int_t)
{
  // Execute action corresponding to one event
  // This member function is called when a AliEMCALRecPoint is clicked with the locator
  //
  // If Left button is clicked on AliEMCALRecPoint, the digits are switched on    
  // and switched off when the mouse button is released.

  //  static Int_t pxold, pyold;

  /* static TGraph *  digitgraph = 0 ;
  static TPaveText* clustertext = 0 ;
  
  if (!gPad->IsEditable()) return;
  
  switch (event) {
    
    
  case kButton1Down:{
    AliEMCALDigit * digit ;
    AliEMCALGeometry * emcalgeom =  (AliEMCALGetter::Instance())->EMCALGeometry() ;

    Int_t iDigit;
    Int_t relid[2] ;
  
    const Int_t kMulDigit=AliEMCALRecPoint::GetDigitsMultiplicity() ;
    Float_t * xi = new Float_t [kMulDigit] ; 
    Float_t * zi = new Float_t [kMulDigit] ;
    
    for(iDigit = 0; iDigit < kMulDigit; iDigit++) {
      Fatal("AliEMCALRecPoint::ExecuteEvent", " -> Something wrong with the code"); 
      digit = 0 ; //dynamic_cast<AliEMCALDigit *>((fDigitsList)[iDigit]);
      emcalgeom->AbsToRelNumbering(digit->GetId(), relid) ;
      emcalgeom->PosInAlice(relid, xi[iDigit], zi[iDigit]) ;
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
      sprintf(line1,"Energy=%1.2f GeV",GetEnergy()) ;
      sprintf(line2,"%d Digits",GetDigitsMultiplicity()) ;
      clustertext ->AddText(line1) ;
      clustertext ->AddText(line2) ;
      clustertext ->Draw("");
    }
    gPad->Update() ; 
    Print("") ;
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
    
    }*/
}
//____________________________________________________________________________
void AliEMCALRecPoint::EvalAll(Float_t logWeight,TClonesArray * digits) 
{
  // Evaluates all shower parameters

  EvalLocalPosition(logWeight, digits) ;
  //  printf("eval position done\n");
  EvalElipsAxis(logWeight, digits) ;
  //  printf("eval axis done\n");
  EvalDispersion(logWeight, digits) ;
  //  printf("eval dispersion done\n");
 // EvalCoreEnergy(logWeight, digits);
 // printf("eval energy done\n");
  EvalTime(digits) ;
  //  printf("eval time done\n");

  EvalPrimaries(digits) ;
  //  printf("eval pri done\n");
  EvalParents(digits);
  //  printf("eval parent done\n");
}

//____________________________________________________________________________
void  AliEMCALRecPoint::EvalDispersion(Float_t logWeight, TClonesArray * digits)
{
  // Calculates the dispersion of the shower at the origin of the RecPoint

  Float_t d    = 0. ;
  Float_t wtot = 0. ;

  AliEMCALDigit * digit ;
 
  AliEMCALGeometry * geom = (AliEMCALGetter::Instance())->EMCALGeometry();
  
  // Calculates the centre of gravity in the local EMCAL-module coordinates   
  Int_t iDigit;

  if (!fLocPos.X() || !fLocPos.Y()) 
    EvalLocalPosition(logWeight, digits) ;
  
  //Sub const Float_t kDeg2Rad = TMath::DegToRad() ; 
    
  Float_t cluEta =  fLocPos.X() ; 
  Float_t cluPhi =  fLocPos.Y() ; 
  Float_t cluR =  fLocPos.Z() ; 
  
  if (gDebug == 2) 
    printf("EvalDispersion: eta,phi,r = %f,%f,%f", cluEta, cluPhi, cluR) ;
  
  // Calculates the dispersion in coordinates 
  wtot = 0.;
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;
    Float_t etai = 0.;
    Float_t phii = 0.;
    //copied for shish-kebab geometry, ieta,iphi is cast as float as eta,phi conversion
    // for this geometry does not exist
    int nSupMod=0, nTower=0, nIphi=0, nIeta=0;
    int iphi=0, ieta=0;
       geom->GetCellIndex(digit->GetId(), nSupMod,nTower,nIphi,nIeta);
       geom->GetCellPhiEtaIndexInSModule(nSupMod,nTower,nIphi,nIeta, iphi,ieta);
	etai=(Float_t)ieta;
	phii=(Float_t)iphi;
	//        printf("%f,%d,%d \n", fAmp, ieta, iphi) ;

/* Sub
  geom->EtaPhiFromIndex(digit->GetId(), etai, phii);
    phii = phii * kDeg2Rad;
*/
///////////////////////////
//  if(fAmp>0)printf("%f %d %f", fAmp,iDigit,fEnergyList[iDigit]) ;
/////////////////////////

    if (gDebug == 2) 
      printf("EvalDispersion: id = %d, etai,phii = %f,%f", digit->GetId(), etai, phii) ;

    Float_t w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
    d += w * ( (etai-cluEta)*(etai-cluEta) + (phii-cluPhi)*(phii-cluPhi) ) ; 
    wtot+=w ;
  }
  
  if ( wtot > 0 ) 
    d /= wtot ;
  else 
    d = 0. ; 

  fDispersion = TMath::Sqrt(d) ;
  //      printf("Dispersion: = %f", fDispersion) ;
}

//____________________________________________________________________________
void AliEMCALRecPoint::EvalLocalPosition(Float_t logWeight, TClonesArray * digits)
{
  // Calculates the center of gravity in the local EMCAL-module coordinates 
  Float_t wtot = 0. ;
 
  //  Int_t relid[3] ;
  
  AliEMCALDigit * digit ;
  AliEMCALGeometry * geom  =  (AliEMCALGetter::Instance())->EMCALGeometry();
  Int_t iDigit;
  Float_t cluEta = 0;
  Float_t cluPhi = 0;
 //Sub const Float_t kDeg2Rad = TMath::DegToRad();

  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(fDigitsList[iDigit])) ;

    Float_t etai ;
    Float_t phii ;
    //copied for shish-kebab geometry, ieta,iphi is cast as float as eta,phi conversion
    // for this geometry does not exist
    int nSupMod=0, nTower=0, nIphi=0, nIeta=0;
    int iphi=0, ieta=0;
       geom->GetCellIndex(digit->GetId(), nSupMod,nTower,nIphi,nIeta);
       geom->GetCellPhiEtaIndexInSModule(nSupMod, nTower,nIphi,nIeta, iphi,ieta); //19-oct-05
	etai=(Float_t)ieta;
	phii=(Float_t)iphi;
//Sub    geom->EtaPhiFromIndex(digit->GetId(), etai, phii);
//Sub    phii = phii * kDeg2Rad;
    Float_t w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ) ) ;
    cluEta += (etai * w) ;
    cluPhi += (phii * w );
    wtot += w ;
  }

  if ( wtot > 0 ) { 
    cluEta /= wtot ;
    cluPhi /= wtot ;
  } else {
    cluEta = -1 ; 
    cluPhi = -1.; 
  }
  
  fLocPos.SetX(cluEta);
  fLocPos.SetY(cluPhi);
  fLocPos.SetZ(geom->GetIP2ECASection());
    
//  if (gDebug==2)
//    printf("EvalLocalPosition: eta,phi,r = %f,%f,%f", fLocPos.X(), fLocPos.Y(), fLocPos.Z()) ; 
  fLocPosM = 0 ;
}

//______________________________________________________________________________
void AliEMCALRecPoint::EvalCoreEnergy(Float_t logWeight, TClonesArray * digits)
{
  // This function calculates energy in the core, 
  // i.e. within a radius rad = 3cm around the center. Beyond this radius
  // in accordance with shower profile the energy deposition 
  // should be less than 2%

  AliEMCALDigit * digit ;
  const Float_t kDeg2Rad = TMath::DegToRad() ;
  AliEMCALGeometry * geom = (AliEMCALGetter::Instance())->EMCALGeometry();    
  Int_t iDigit;

  if (!fLocPos.X() || !fLocPos.Y() ) {
    EvalLocalPosition(logWeight, digits);
  }
  
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) ( digits->At(fDigitsList[iDigit]) ) ;
    Float_t etai = 0. ;
    Float_t phii = 0. ;
    geom->PosInAlice(digit->GetId(), etai, phii);
    phii = phii * kDeg2Rad;
  
    Float_t distance = TMath::Sqrt((etai-fLocPos.X())*(etai-fLocPos.X())+(phii-fLocPos.Y())*(phii-fLocPos.Y())) ;
    if(distance < fCoreRadius)
      fCoreEnergy += fEnergyList[iDigit] ;
  }
  
}
//____________________________________________________________________________
void  AliEMCALRecPoint::EvalElipsAxis(Float_t logWeight,TClonesArray * digits)
{
  // Calculates the axis of the shower ellipsoid in eta and phi

  Double_t wtot = 0. ;
  Double_t x    = 0.;
  Double_t z    = 0.;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;

  const Float_t kDeg2Rad = TMath::DegToRad();
  AliEMCALDigit * digit ;

  AliEMCALGeometry * geom = (AliEMCALGetter::Instance())->EMCALGeometry();
  TString gn(geom->GetName());

  Int_t iDigit;
  
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;
    Float_t etai = 0. ;
    Float_t phii = 0. ; 
    if(gn.Contains("SHISH")) {
    //copied for shish-kebab geometry, ieta,iphi is cast as float as eta,phi conversion
    // for this geometry does not exist
      int nSupMod=0, nTower=0, nIphi=0, nIeta=0;
      int iphi=0, ieta=0;
      geom->GetCellIndex(digit->GetId(), nSupMod,nTower,nIphi,nIeta);
      geom->GetCellPhiEtaIndexInSModule(nSupMod,nTower,nIphi,nIeta, iphi,ieta);
      etai=(Float_t)ieta;
      phii=(Float_t)iphi;
    } else {
      geom->EtaPhiFromIndex(digit->GetId(), etai, phii);
      phii = phii * kDeg2Rad;
    }

    Double_t w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;

    dxx  += w * etai * etai ;
    x    += w * etai ;
    dzz  += w * phii * phii ;
    z    += w * phii ; 

    dxz  += w * etai * phii ; 

    wtot += w ;
  }

  if ( wtot > 0 ) { 
    dxx /= wtot ;
    x   /= wtot ;
    dxx -= x * x ;
    dzz /= wtot ;
    z   /= wtot ;
    dzz -= z * z ;
    dxz /= wtot ;
    dxz -= x * z ;

    fLambda[0] =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    if(fLambda[0] > 0)
      fLambda[0] = TMath::Sqrt(fLambda[0]) ;
    else
      fLambda[0] = 0;
    
    fLambda[1] =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;

    if(fLambda[1] > 0) //To avoid exception if numerical errors lead to negative lambda.
      fLambda[1] = TMath::Sqrt(fLambda[1]) ;
    else
      fLambda[1]= 0. ;
  } else { 
    fLambda[0]= 0. ;
    fLambda[1]= 0. ;
  }

  //  printf("Evalaxis: lambdas  = %f,%f", fLambda[0],fLambda[1]) ; 

}

//______________________________________________________________________________
void  AliEMCALRecPoint::EvalPrimaries(TClonesArray * digits)
{
  // Constructs the list of primary particles (tracks) which have contributed to this RecPoint
  
  AliEMCALDigit * digit ;
  Int_t * tempo    = new Int_t[fMaxTrack] ;

  Int_t index ;  
  for ( index = 0 ; index < GetDigitsMultiplicity() ; index++ ) { // all digits
    digit = dynamic_cast<AliEMCALDigit *>(digits->At( fDigitsList[index] )) ; 
    Int_t nprimaries = digit->GetNprimary() ;
    if ( nprimaries == 0 ) continue ;
    Int_t * newprimaryarray = new Int_t[nprimaries] ;
    Int_t ii ; 
    for ( ii = 0 ; ii < nprimaries ; ii++)
      newprimaryarray[ii] = digit->GetPrimary(ii+1) ; 

    Int_t jndex ;
    for ( jndex = 0 ; jndex < nprimaries ; jndex++ ) { // all primaries in digit
      if ( fMulTrack > fMaxTrack ) {
	fMulTrack = fMaxTrack ;
	Error("GetNprimaries", "increase fMaxTrack ")  ;
	break ;
      }
      Int_t newprimary = newprimaryarray[jndex] ;
      Int_t kndex ;
      Bool_t already = kFALSE ;
      for ( kndex = 0 ; kndex < fMulTrack ; kndex++ ) { //check if not already stored
	if ( newprimary == tempo[kndex] ){
	  already = kTRUE ;
	  break ;
	}
      } // end of check
      if ( !already && (fMulTrack < fMaxTrack)) { // store it
	tempo[fMulTrack] = newprimary ; 
	fMulTrack++ ;
      } // store it
    } // all primaries in digit
    delete [] newprimaryarray ; 
  } // all digits

  
  fTracksList = new Int_t[fMulTrack] ;
  for(index = 0; index < fMulTrack; index++)
   fTracksList[index] = tempo[index] ;
 
  delete [] tempo ;

}

//______________________________________________________________________________
void  AliEMCALRecPoint::EvalParents(TClonesArray * digits)
{
  // Constructs the list of parent particles (tracks) which have contributed to this RecPoint
 
  AliEMCALDigit * digit ;
  Int_t * tempo    = new Int_t[fMaxParent] ;

  Int_t index ;  
  for ( index = 0 ; index < GetDigitsMultiplicity() ; index++ ) { // all digits
    digit = dynamic_cast<AliEMCALDigit *>(digits->At( fDigitsList[index] )) ; 
    Int_t nparents = digit->GetNiparent() ;
    if ( nparents == 0 ) continue ;
    Int_t * newparentarray = new Int_t[nparents] ;
    Int_t ii ; 
    for ( ii = 0 ; ii < nparents ; ii++)
      newparentarray[ii] = digit->GetIparent(ii+1) ; 

    Int_t jndex ;
    for ( jndex = 0 ; jndex < nparents ; jndex++ ) { // all primaries in digit
      if ( fMulParent > fMaxParent ) {
	fMulTrack = - 1 ;
	Error("GetNiparent", "increase fMaxParent")  ;
	break ;
      }
      Int_t newparent = newparentarray[jndex] ;
      Int_t kndex ;
      Bool_t already = kFALSE ;
      for ( kndex = 0 ; kndex < fMulParent ; kndex++ ) { //check if not already stored
	if ( newparent == tempo[kndex] ){
	  already = kTRUE ;
	  break ;
	}
      } // end of check
      if ( !already && (fMulTrack < fMaxTrack)) { // store it
	tempo[fMulParent] = newparent ; 
	fMulParent++ ;
      } // store it
    } // all parents in digit
    delete [] newparentarray ; 
  } // all digits

  if (fMulParent>0) {
    fParentsList = new Int_t[fMulParent] ;
    for(index = 0; index < fMulParent; index++)
      fParentsList[index] = tempo[index] ;
  }
 
  delete [] tempo ;

}

//____________________________________________________________________________
void AliEMCALRecPoint::GetLocalPosition(TVector3 & lpos) const
{
  // returns the position of the cluster in the local reference system of ALICE
  // X = eta, Y = phi, Z = r (a constant for the EMCAL)
  
  lpos.SetX(fLocPos.X()) ;
  lpos.SetY(fLocPos.Y()) ;
  lpos.SetZ(fLocPos.Z()) ;
}

//____________________________________________________________________________
void AliEMCALRecPoint::GetGlobalPosition(TVector3 & gpos) const
{
  // returns the position of the cluster in the global reference system of ALICE
  // These are now the Cartesian X, Y and Z

  AliEMCALGeometry * geom = (AliEMCALGetter::Instance())->EMCALGeometry();
  Int_t absid = geom->TowerIndexFromEtaPhi(fLocPos.X(), TMath::RadToDeg()*fLocPos.Y());
  geom->XYZFromIndex(absid, gpos);
}

//____________________________________________________________________________
Float_t AliEMCALRecPoint::GetMaximalEnergy(void) const
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
Int_t AliEMCALRecPoint::GetMultiplicityAtLevel(Float_t H) const
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
Int_t  AliEMCALRecPoint::GetNumberOfLocalMax(AliEMCALDigit **  maxAt, Float_t * maxAtEnergy,
					       Float_t locMaxCut,TClonesArray * digits) const
{ 
  // Calculates the number of local maxima in the cluster using fLocalMaxCut as the minimum
  // energy difference between two local maxima

  AliEMCALDigit * digit ;
  AliEMCALDigit * digitN ;
  
  Int_t iDigitN ;
  Int_t iDigit ;

  for(iDigit = 0; iDigit < fMulDigit; iDigit++)
    maxAt[iDigit] = (AliEMCALDigit*) digits->At(fDigitsList[iDigit])  ;
  
  for(iDigit = 0 ; iDigit < fMulDigit; iDigit++) {   
    if(maxAt[iDigit]) {
      digit = maxAt[iDigit] ;
          
      for(iDigitN = 0; iDigitN < fMulDigit; iDigitN++) {	
	digitN = (AliEMCALDigit *) digits->At(fDigitsList[iDigitN]) ; 
	
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
    if(maxAt[iDigit] ){
      maxAt[iDigitN] = maxAt[iDigit] ;
      maxAtEnergy[iDigitN] = fEnergyList[iDigit] ;
      iDigitN++ ; 
    }
  }
  return iDigitN ;
}
//____________________________________________________________________________
void AliEMCALRecPoint::EvalTime(TClonesArray * digits){
  // time is set to the time of the digit with the maximum energy

  Float_t maxE = 0;
  Int_t maxAt = 0;
  for(Int_t idig=0; idig < fMulDigit; idig++){
    if(fEnergyList[idig] > maxE){
      maxE = fEnergyList[idig] ;
      maxAt = idig;
    }
  }
  fTime = ((AliEMCALDigit*) digits->At(fDigitsList[maxAt]))->GetTime() ;
  
}

//______________________________________________________________________________
void AliEMCALRecPoint::Paint(Option_t *)
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
Float_t AliEMCALRecPoint::EtaToTheta(Float_t arg) const
{
  //Converts Theta (Radians) to Eta(Radians)
  return (2.*TMath::ATan(TMath::Exp(-arg)));
}

//______________________________________________________________________________
Float_t AliEMCALRecPoint::ThetaToEta(Float_t arg) const
{
  //Converts Eta (Radians) to Theta(Radians)
  return (-1 * TMath::Log(TMath::Tan(0.5 * arg)));
}

//____________________________________________________________________________
void AliEMCALRecPoint::Print(Option_t *) const
{
  // Print the list of digits belonging to the cluster
  
  TString message ; 
  message  = "AliPHOSEmcRecPoint:\n" ;
  message +=  " digits # = " ; 
  Info("Print", message.Data()) ; 

  Int_t iDigit;
  for(iDigit=0; iDigit<fMulDigit; iDigit++)
    printf(" %d ", fDigitsList[iDigit] ) ;  
  
  Info("Print", " Energies = ") ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) 
    printf(" %f ", fEnergyList[iDigit] ) ;
  printf("\n") ; 
   Info("Print", " Primaries  ") ;
  for(iDigit = 0;iDigit < fMulTrack; iDigit++)
    printf(" %d ", fTracksList[iDigit]) ;
  printf("\n") ; 	
  message  = "       Multiplicity    = %d" ;
  message += "       Cluster Energy  = %f" ; 
  message += "       Core energy     = %f" ; 
  message += "       Core radius     = %f" ; 
  message += "       Number of primaries %d" ; 
  message += "       Stored at position %d" ; 
 
  Info("Print", message.Data(), fMulDigit, fAmp, fCoreEnergy, fCoreRadius, fMulTrack, GetIndexInList() ) ;  
}
