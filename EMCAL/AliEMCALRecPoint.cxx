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
class Riostream;
#include <TPad.h>
class TGraph;
class TPaveText;
#include <TClonesArray.h>
#include <TMath.h>

// --- Standard library ---

// --- AliRoot header files ---
//#include "AliGenerator.h"
class AliGenerator;
#include "AliRunLoader.h"
#include "AliRun.h"
class AliEMCAL;
#include "AliEMCALLoader.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALHit.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRecPoint.h"

ClassImp(AliEMCALRecPoint)

//____________________________________________________________________________
AliEMCALRecPoint::AliEMCALRecPoint()
  : AliRecPoint(),
    fGeomPtr(0),
    fClusterType(-1),
    fCoreEnergy(0),
    fDispersion(0),
    fEnergyList(0),
    fTimeList(0),
    fAbsIdList(0),
    fTime(0.),
    fCoreRadius(10),  //HG check this
    fMulParent(0),
    fMaxParent(0),
    fParentsList(0),
    fSuperModuleNumber(0)
{
  // ctor
  fMaxTrack = 0 ;
  fMulDigit = 0 ;
  fAmp   = 0. ;   
  //  fLocPos.SetX(1.e+6)  ;      //Local position should be evaluated

  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  fGeomPtr = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
  //fGeomPtr = AliEMCALGeometry::GetInstance();
  fGeomPtr->GetTransformationForSM(); // Global <-> Local
}

//____________________________________________________________________________
AliEMCALRecPoint::AliEMCALRecPoint(const char * opt) 
  : AliRecPoint(opt),
    fGeomPtr(0),
    fClusterType(-1),
    fCoreEnergy(0),
    fDispersion(0),
    fEnergyList(0),
    fTimeList(0),
    fAbsIdList(0),
    fTime(-1.),
    fCoreRadius(10),  //HG check this
    fMulParent(0),
    fMaxParent(1000),
    fParentsList(0),
    fSuperModuleNumber(0)
{
  // ctor
  fMaxTrack = 1000 ;
  fMulDigit   = 0 ; 
  fAmp   = 0. ;   
  fParentsList = new Int_t[fMaxParent];

  //fLocPos.SetX(1.e+6)  ;      //Local position should be evaluated
  //fGeomPtr = AliEMCALGeometry::GetInstance();
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  fGeomPtr = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
  fGeomPtr->GetTransformationForSM(); // Global <-> Local
}

//____________________________________________________________________________
AliEMCALRecPoint::AliEMCALRecPoint(const AliEMCALRecPoint & rp) 
  : AliRecPoint(rp),
    fGeomPtr(rp.fGeomPtr),
    fClusterType(rp.fClusterType),
    fCoreEnergy(rp.fCoreEnergy),
    fDispersion(rp.fDispersion),
    fEnergyList(0),
    fTimeList(0),
    fAbsIdList(0),
    fTime(rp.fTime),
    fCoreRadius(rp.fCoreRadius),
    fMulParent(rp.fMulParent),
    fMaxParent(rp.fMaxParent),
    fParentsList(0),
    fSuperModuleNumber(rp.fSuperModuleNumber)
{
  //copy ctor
  fLambda[0] = rp.fLambda[0];
  fLambda[1] = rp.fLambda[1];

  fEnergyList = new Float_t[rp.fMulDigit];
  fTimeList = new Float_t[rp.fMulDigit];
  fAbsIdList = new Int_t[rp.fMulDigit];
  for(Int_t i = 0; i < rp.fMulDigit; i++) {
    fEnergyList[i] = rp.fEnergyList[i];
    fTimeList[i] = rp.fTimeList[i];
    fAbsIdList[i] = rp.fAbsIdList[i];
  }
  fParentsList = new Int_t[rp.fMulParent];
  for(Int_t i = 0; i < rp.fMulParent; i++) fParentsList[i] = rp.fParentsList[i];

}
//____________________________________________________________________________
AliEMCALRecPoint::~AliEMCALRecPoint()
{
  // dtor
  if ( fEnergyList )
    delete[] fEnergyList ; 
  if ( fTimeList )
    delete[] fTimeList ; 
  if ( fAbsIdList )
    delete[] fAbsIdList ; 
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
  if(fTimeList == 0)
    fTimeList =  new Float_t[fMaxDigit]; 
  if(fAbsIdList == 0) {
    fAbsIdList =  new Int_t[fMaxDigit];
    fSuperModuleNumber = fGeomPtr->GetSuperModuleNumber(digit.GetId());
  }

  if ( fMulDigit >= fMaxDigit ) { // increase the size of the lists 
    fMaxDigit*=2 ; 
    Int_t * tempo = new Int_t[fMaxDigit]; 
    Float_t * tempoE =  new Float_t[fMaxDigit];
    Float_t * tempoT =  new Float_t[fMaxDigit];
    Int_t * tempoId = new Int_t[fMaxDigit]; 

    Int_t index ;     
    for ( index = 0 ; index < fMulDigit ; index++ ){
      tempo[index]   = fDigitsList[index] ;
      tempoE[index]  = fEnergyList[index] ; 
      tempoT[index]  = fTimeList[index] ; 
      tempoId[index] = fAbsIdList[index] ; 
    }
    
    delete [] fDigitsList ; 
    fDigitsList =  new Int_t[fMaxDigit];
 
    delete [] fEnergyList ;
    fEnergyList =  new Float_t[fMaxDigit];

    delete [] fTimeList ;
    fTimeList =  new Float_t[fMaxDigit];

    delete [] fAbsIdList ;
    fAbsIdList =  new Int_t[fMaxDigit];

    for ( index = 0 ; index < fMulDigit ; index++ ){
      fDigitsList[index] = tempo[index] ;
      fEnergyList[index] = tempoE[index] ; 
      fTimeList[index] = tempoT[index] ; 
      fAbsIdList[index]  = tempoId[index] ; 
    }
 
    delete [] tempo ;
    delete [] tempoE ; 
    delete [] tempoT ; 
    delete [] tempoId ; 
  } // if
  
  fDigitsList[fMulDigit]   = digit.GetIndexInList()  ; 
  fEnergyList[fMulDigit]   = Energy ;
  fTimeList[fMulDigit]     = digit.GetTime() ;
  fAbsIdList[fMulDigit]    = digit.GetId();
  fMulDigit++ ; 
  fAmp += Energy ; 

}
//____________________________________________________________________________
Bool_t AliEMCALRecPoint::AreNeighbours(AliEMCALDigit * digit1, AliEMCALDigit * digit2 ) const
{
  // Tells if (true) or not (false) two digits are neighbours
  // A neighbour is defined as being two digits which share a corner
  
  static Bool_t areNeighbours = kFALSE ;
  static Int_t nSupMod=0, nTower=0, nIphi=0, nIeta=0;
  static int nSupMod1=0, nTower1=0, nIphi1=0, nIeta1=0;
  static Int_t relid1[2] , relid2[2] ; // ieta, iphi
  static Int_t rowdiff=0, coldiff=0;

  areNeighbours = kFALSE ;

  fGeomPtr->GetCellIndex(digit1->GetId(), nSupMod,nTower,nIphi,nIeta);
  fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nTower,nIphi,nIeta, relid1[0],relid1[1]);

  fGeomPtr->GetCellIndex(digit2->GetId(), nSupMod1,nTower1,nIphi1,nIeta1);
  fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod1,nTower1,nIphi1,nIeta1, relid2[0],relid2[1]);
  
  rowdiff = TMath::Abs( relid1[0] - relid2[0] ) ;  
  coldiff = TMath::Abs( relid1[1] - relid2[1] ) ;  

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

  Int_t rowdif = (Int_t)(TMath::Ceil(locpos1.X()/delta)-TMath::Ceil(locpos2.X()/delta)) ;
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
  //EvalCoreEnergy(logWeight, digits);
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
  // in cell units - Nov 16,2006

  Double_t d = 0., wtot = 0., w = 0.;
  Int_t iDigit=0, nstat=0, i=0;
  AliEMCALDigit * digit ;
  
  // Calculates the dispersion in cell units 
  Double_t etai, phii, etaMean=0.0, phiMean=0.0; 
  int nSupMod=0, nTower=0, nIphi=0, nIeta=0;
  int iphi=0, ieta=0;
  // Calculate mean values
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;

    if (fAmp>0 && fEnergyList[iDigit]>0) {
      fGeomPtr->GetCellIndex(digit->GetId(), nSupMod,nTower,nIphi,nIeta);
      fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nTower,nIphi,nIeta, iphi,ieta);
      etai=(Double_t)ieta;
      phii=(Double_t)iphi;
      w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;

      if(w>0.0) {
        phiMean += phii*w;
        etaMean += etai*w;
        wtot    += w;
      }
    }
  }
  if (wtot>0) {
    phiMean /= wtot ;
    etaMean /= wtot ;
  } else AliError(Form("Wrong weight %f\n", wtot));

  // Calculate dispersion
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;

    if (fAmp>0 && fEnergyList[iDigit]>0) {
      fGeomPtr->GetCellIndex(digit->GetId(), nSupMod,nTower,nIphi,nIeta);
      fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nTower,nIphi,nIeta, iphi,ieta);
      etai=(Double_t)ieta;
      phii=(Double_t)iphi;
      w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;

      if(w>0.0) {
        nstat++;
        d += w*((etai-etaMean)*(etai-etaMean)+(phii-phiMean)*(phii-phiMean));
      }
    }
  }
  
  if ( wtot > 0 && nstat>1) d /= wtot ;
  else                      d = 0. ; 

  fDispersion = TMath::Sqrt(d) ;
}

//____________________________________________________________________________
void AliEMCALRecPoint::EvalLocalPosition(Float_t logWeight, TClonesArray * digits)
{
  // Calculates the center of gravity in the local EMCAL-module coordinates 
  //  Info("Print", " logWeight %f : cluster energy %f ", logWeight, fAmp); // for testing
  
  AliEMCALDigit * digit;
  Int_t i=0, nstat=0;
  Double_t clXYZ[3]={0.,0.,0.}, clRmsXYZ[3]={0.,0.,0.}, xyzi[3], wtot=0., w=0.;

  for(Int_t iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(fDigitsList[iDigit])) ;

    fGeomPtr->RelPosCellInSModule(digit->GetId(), xyzi[0], xyzi[1], xyzi[2]);
    // printf(" Id %i : Local x,y,z %f %f %f \n", digit->GetId(), xyzi[0], xyzi[1], xyzi[2]);

    if(logWeight > 0.0) w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ));
    else                w = fEnergyList[iDigit]; // just energy

    if(w>0.0) {
      wtot += w ;
      nstat++;
      for(i=0; i<3; i++ ) {
        clXYZ[i]    += (w*xyzi[i]);
        clRmsXYZ[i] += (w*xyzi[i]*xyzi[i]);
      }
    }
  }
  //  cout << " wtot " << wtot << endl;
  if ( wtot > 0 ) { 
    //    xRMS   = TMath::Sqrt(x2m - xMean*xMean);
    for(i=0; i<3; i++ ) {
      clXYZ[i] /= wtot;
      if(nstat>1) {
        clRmsXYZ[i] /= (wtot*wtot);  
        clRmsXYZ[i]  =  clRmsXYZ[i] - clXYZ[i]*clXYZ[i];
	if(clRmsXYZ[i] > 0.0) {
          clRmsXYZ[i] =  TMath::Sqrt(clRmsXYZ[i]);
        } else      clRmsXYZ[i] = 0;
      } else        clRmsXYZ[i] = 0;
    }    
  } else {
    for(i=0; i<3; i++ ) {
      clXYZ[i] = clRmsXYZ[i] = -1.;
    }
  }
  // clRmsXYZ[i] ??
  fLocPos.SetX(clXYZ[0]);
  fLocPos.SetY(clXYZ[1]);
  fLocPos.SetZ(clXYZ[2]);
    
//  if (gDebug==2)
//    printf("EvalLocalPosition: eta,phi,r = %f,%f,%f", fLocPos.X(), fLocPos.Y(), fLocPos.Z()) ; 
  fLocPosM = 0 ; // covariance matrix
}

//void AliEMCALRecPoint::EvalLocalPositionSimple()
//{ // Weight is proportional of cell energy 
//}

//______________________________________________________________________________
void AliEMCALRecPoint::EvalCoreEnergy(Float_t logWeight, TClonesArray * digits)
{
  // This function calculates energy in the core, 
  // i.e. within a radius rad = fCoreEnergy around the center. Beyond this radius
  // in accordance with shower profile the energy deposition 
  // should be less than 2%
  // Unfinished - Nov 15,2006
  // Distance is calculate in (phi,eta) units

  AliEMCALDigit * digit ;

  Int_t iDigit;

  if (!fLocPos.Mag()) {
    EvalLocalPosition(logWeight, digits);
  }
  
  Double_t phiPoint = fLocPos.Phi(), etaPoint = fLocPos.Eta();
  Double_t eta, phi, distance;
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) ( digits->At(fDigitsList[iDigit]) ) ;
    
    eta = phi = 0.0;
    fGeomPtr->EtaPhiFromIndex(digit->GetId(),eta, phi) ;
    phi = phi * TMath::DegToRad();
  
    distance = TMath::Sqrt((eta-etaPoint)*(eta-etaPoint)+(phi-phiPoint)*(phi-phiPoint));
    if(distance < fCoreRadius)
      fCoreEnergy += fEnergyList[iDigit] ;
  }
  
}
//____________________________________________________________________________
void  AliEMCALRecPoint::EvalElipsAxis(Float_t logWeight,TClonesArray * digits)
{
  // Calculates the axis of the shower ellipsoid in eta and phi
  // in cell units

  static TString gn(fGeomPtr->GetName());

  Double_t wtot = 0.;
  Double_t x    = 0.;
  Double_t z    = 0.;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;

  AliEMCALDigit * digit = 0;

  Double_t etai , phii, w; 
  int nSupMod=0, nTower=0, nIphi=0, nIeta=0;
  int iphi=0, ieta=0;
  for(Int_t iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;
    etai = phii = 0.; 
    if(gn.Contains("SHISH")) { 
    // Nov 15,2006 - use cell numbers as coordinates
    // Copied for shish-kebab geometry, ieta,iphi is cast as double as eta,phi
    // We can use the eta,phi(or coordinates) of cell
      nSupMod = nTower = nIphi = nIeta = iphi = ieta = 0;

      fGeomPtr->GetCellIndex(digit->GetId(), nSupMod,nTower,nIphi,nIeta);
      fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nTower,nIphi,nIeta, iphi,ieta);
      etai=(Double_t)ieta;
      phii=(Double_t)iphi;
    } else { // 
      fGeomPtr->EtaPhiFromIndex(digit->GetId(), etai, phii);
      phii = phii * TMath::DegToRad();
    }

    w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
    // fAmp summed amplitude of digits, i.e. energy of recpoint 
    // Gives smaller value of lambda than log weight  
    // w = fEnergyList[iDigit] / fAmp; // Nov 16, 2006 - try just energy

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
  
  lpos.SetX(fLocPos.X()) ;
  lpos.SetY(fLocPos.Y()) ;
  lpos.SetZ(fLocPos.Z()) ;
}

//____________________________________________________________________________
void AliEMCALRecPoint::GetGlobalPosition(TVector3 & gpos) const
{
  // returns the position of the cluster in the global reference system of ALICE
  // These are now the Cartesian X, Y and Z
  //  cout<<" geom "<<geom<<endl;
  fGeomPtr->GetGlobal(fLocPos, gpos, fSuperModuleNumber);
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
Int_t AliEMCALRecPoint::GetPrimaryIndex() const  
{
  // Get the primary track index in TreeK which deposits the most energy 
  // in Digits which forms RecPoint. Kinematics, Hits and Digits must be 
  // loaded before the call of the method.

  AliRunLoader *rl = AliRunLoader::GetRunLoader(); 
  if (!rl) 
    AliError(Form(" No Runloader ")) ; 
 
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>
    (rl->GetDetectorLoader("EMCAL"));

  // Get the list of digits forming this RecPoint
  Int_t  nDigits   = fMulDigit   ;
  Int_t *digitList = fDigitsList ;
  
  // Find the digit with maximum amplitude
  AliEMCALDigit *digit = 0;
  TClonesArray *digits = emcalLoader->Digits();
  Int_t maxAmp = 0;
  Int_t bestDigitIndex = -1;
  for (Int_t iDigit=0; iDigit<nDigits; iDigit++) {
    digit = static_cast<AliEMCALDigit *>(digits->At(digitList[iDigit]));
    if (digit->GetAmp() > maxAmp) {
      maxAmp = digit->GetAmp();
      bestDigitIndex = iDigit;
    }
  }

  digit = static_cast<AliEMCALDigit *>(digits->At(digitList[bestDigitIndex]));  

  // Get the list of hits producing this digit,
  // find which hit has deposited more energy 
  // and find the primary track.

  AliEMCALHit *hit = 0;
  TClonesArray *hits = emcalLoader->Hits();

  Double_t maxedep  =  0;
  Int_t    maxtrack = -1;
  Int_t    nHits    = hits ->GetEntries();
  Int_t    id       = digit->GetId();
  for (Int_t iHit=0; iHit<nHits; iHit++) {
    hit = static_cast<AliEMCALHit*> (hits->At(iHit)) ;
    if(hit->GetId() == id){
      Double_t edep  = hit->GetEnergy();
      Int_t    track = hit->GetIparent();//Primary();
      if(edep > maxedep){
	maxedep  = edep;
	maxtrack = track;
      }
    }
  }
  if (maxtrack != -1) return maxtrack; 
  return -12345;                       // no track found :(
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
  return;
  TString message ; 
  message  = "AliEMCALRecPoint:\n" ;
  message +=  " digits # = " ; 
  Info("Print", message.Data()) ; 

  Int_t iDigit;
  for(iDigit=0; iDigit<fMulDigit; iDigit++)
    printf(" %d ", fDigitsList[iDigit] ) ;  
  printf("\n");

  Info("Print", " Energies = ") ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) 
    printf(" %f ", fEnergyList[iDigit] ) ;
  printf("\n");

  Info("Print", "\n Abs Ids = ") ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) 
    printf(" %i ", fAbsIdList[iDigit] ) ;
  printf("\n");

  Info("Print", " Primaries  ") ;
  for(iDigit = 0;iDigit < fMulTrack; iDigit++)
    printf(" %d ", fTracksList[iDigit]) ;

  printf("\n Local x %6.2f y %7.2f z %7.1f \n", fLocPos[0], fLocPos[1], fLocPos[2]);

  message  = "       ClusterType     = %d" ;
  message += "       Multiplicity    = %d" ;
  message += "       Cluster Energy  = %f" ; 
  message += "       Core energy     = %f" ; 
  message += "       Core radius     = %f" ; 
  message += "       Number of primaries %d" ; 
  message += "       Stored at position %d" ; 
  Info("Print", message.Data(), fClusterType, fMulDigit, fAmp, fCoreEnergy, fCoreRadius, fMulTrack, GetIndexInList() ) ;  
}

Double_t  AliEMCALRecPoint::GetPointEnergy() const
{
  static double e;
  e=0.0;
  for(int ic=0; ic<GetMultiplicity(); ic++) e += double(fEnergyList[ic]);
  return e;
}
