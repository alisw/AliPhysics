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

//_________________________________________________________________________
// Algorithm class to construct track segments connection RecPoints in 
// EMCA and Ppsd. Unfolds also the clusters in EMCA. 
//*-- Author : D. Peressounko  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TObjArray.h"
#include "TClonesArray.h"
#include "TObjectTable.h"

// --- Standard library ---

#include <iostream>
#include <cassert>

// --- AliRoot header files ---

#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSLink.h"
#include "AliPHOSv0.h"
#include "AliRun.h"

extern void UnfoldingChiSquare(Int_t &NPar, Double_t *Grad, Double_t & fret, Double_t *x, Int_t iflag) ; 

ClassImp( AliPHOSTrackSegmentMakerv1) 


//____________________________________________________________________________
 AliPHOSTrackSegmentMakerv1::AliPHOSTrackSegmentMakerv1() 
{
  // ctor
  fR0 = 4. ;   
  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;
  //clusters are sorted in "rows" and "columns" of width geom->GetCrystalSize(0),
  fDelta = fR0 + geom->GetCrystalSize(0) ;
  fMinuit = new TMinuit(100) ; 
}

//____________________________________________________________________________
 AliPHOSTrackSegmentMakerv1::~AliPHOSTrackSegmentMakerv1()
{ 
  // dtor
  delete fMinuit ; 
}
//____________________________________________________________________________
Bool_t  AliPHOSTrackSegmentMakerv1::FindFit(AliPHOSEmcRecPoint * emcRP, int * maxAt, Float_t * maxAtEnergy,
				    Int_t NPar, Float_t * FitParameters)
{ 
  // gObjectTable->Print() ; 
  // Calls TMinuit for fitting cluster with several maxima 
  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;
  assert( NPar < 100 ) ; 

  gMinuit->SetPrintLevel(-1) ;           // No Printout
  gMinuit->SetFCN(UnfoldingChiSquare) ;  // To set the address of the minimization function 
  gMinuit->SetObjectFit(emcRP) ;         // To tranfer pointer to UnfoldingChiSquare

  // filling initial values for fit parameters
  AliPHOSDigit * digit ;

  Int_t ierflg  = 0; 
  Int_t index   = 0 ;
  Int_t NDigits = (Int_t) NPar / 3 ;

  Int_t iDigit ;


  for(iDigit = 0 ; iDigit < NDigits ; iDigit++){
    digit = (AliPHOSDigit *) maxAt[iDigit]; 

    Int_t RelId[4] ;
    Float_t x ;
    Float_t z ;
    geom->AbsToRelNumbering(digit->GetId(), RelId) ;
    geom->RelPosInModule(RelId, x, z) ;

    Float_t Energy = maxAtEnergy[iDigit] ;

    gMinuit->mnparm(index, "x",  x, 0.1, 0, 0, ierflg) ;
    index++ ;   
    if(ierflg != 0){ 
      cout << "PHOS Unfolding>  Unable to set initial value for fit procedure : x = " << x << endl ;
      return kFALSE;
    }
    gMinuit->mnparm(index, "z",  z, 0.1, 0, 0, ierflg) ;
    index++ ;   
    if(ierflg != 0){
      cout << "PHOS Unfolding>  Unable to set initial value for fit procedure : z = " << z << endl ;
      return kFALSE;
    }
    gMinuit->mnparm(index, "Energy",  Energy , 0.05*Energy, 0., 4.*Energy, ierflg) ;
    index++ ;   
    if(ierflg != 0){
      cout << "PHOS Unfolding>  Unable to set initial value for fit procedure : Energy = " << Energy << endl ;      
      return kFALSE;
}
  }

  Double_t p0 = 0.1 ; // "Tolerance" Evaluation stops when EDM = 0.0001*p0 ; The number of function call slightly
                      //  depends on it. 
  Double_t p1 = 1.0 ;
  Double_t p2 = 0.0 ;

  gMinuit->mnexcm("SET STR", &p2, 0, ierflg) ;   // force TgMinuit to reduce function calls  
  gMinuit->mnexcm("SET GRA", &p1, 1, ierflg) ;   // force TMinuit to use my gradient  
  gMinuit->SetMaxIterations(5);
  gMinuit->mnexcm("SET NOW", &p2 , 0, ierflg) ;  // No Warnings
  gMinuit->mnexcm("MIGRAD", &p0, 0, ierflg) ;    // minimize 
  if(ierflg == 4){  // Minimum not found   
    cout << "PHOS Unfolding>  Fit not converged, cluster abandoned "<< endl ;      
    return kFALSE ;
  }            
  for(index = 0; index < NPar; index++){
    Double_t err ;
    Double_t val ;
    gMinuit->GetParameter(index, val, err) ;    // Returns value and error of parameter index
    FitParameters[index] = val ;
   }
  return kTRUE;

}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::FillOneModule(DigitsList * Dl, RecPointsList * emcIn, TObjArray * emcOut, 
					RecPointsList * ppsdIn, TObjArray * ppsdOutUp,
					TObjArray * ppsdOutLow, Int_t & PHOSMod, Int_t & emcStopedAt, 
					Int_t & ppsdStopedAt)
{
  // Unfold clusters and fill xxxOut arrays with clusters from one PHOS module
 
  AliPHOSEmcRecPoint *  emcRecPoint  ; 
  AliPHOSPpsdRecPoint * ppsdRecPoint ;
  Int_t index ;
  cout << "Fill 1" << endl ;
  Int_t NemcUnfolded = emcIn->GetEntries() ;
  for(index = emcStopedAt; index < NemcUnfolded; index++){
    emcRecPoint = (AliPHOSEmcRecPoint *) (*emcIn)[index] ;
    cout << "Fill 2" << endl ;
   
    if(emcRecPoint->GetPHOSMod() != PHOSMod )  
       break ;
    
    
    Int_t NMultipl = emcRecPoint->GetMultiplicity() ; 
    int maxAt[NMultipl] ;
    Float_t maxAtEnergy[NMultipl] ;
    Int_t Nmax = emcRecPoint->GetNumberOfLocalMax(maxAt, maxAtEnergy) ;

    

    if(Nmax <= 1)     // if cluster is very flat, so that no prononsed maximum, then Nmax = 0 
      emcOut->Add(emcRecPoint) ;
    else {
      UnfoldClusters(Dl, emcIn, emcRecPoint, Nmax, maxAt, maxAtEnergy, emcOut) ;
      emcIn->Remove(emcRecPoint); 
      cout << "Fill 3" << endl ;
      emcIn->Compress() ;
      NemcUnfolded-- ;
      index-- ;
    }
  }
  emcStopedAt = index ;

  for(index = ppsdStopedAt; index < ppsdIn->GetEntries(); index++){
    ppsdRecPoint = (AliPHOSPpsdRecPoint *) (*ppsdIn)[index] ;
    if(ppsdRecPoint->GetPHOSMod() != PHOSMod )   
      break ;
    if(ppsdRecPoint->GetUp() ) 
      ppsdOutUp->Add(ppsdRecPoint) ;
    else  
      ppsdOutLow->Add(ppsdRecPoint) ;
  }
  ppsdStopedAt = index ;
  
  PHOSMod ++ ; 
  
  emcOut->Sort() ;
  ppsdOutUp->Sort() ;
  ppsdOutLow->Sort() ; 
  
}
//____________________________________________________________________________
Float_t  AliPHOSTrackSegmentMakerv1::GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * EmcClu,AliPHOSPpsdRecPoint * PpsdClu, Bool_t &TooFar)
{
  Float_t R = fR0 ;
 
  TVector3 vecEmc ;
  TVector3 vecPpsd ;
  
  EmcClu->GetLocalPosition(vecEmc) ;
  PpsdClu->GetLocalPosition(vecPpsd)  ; 
  if(EmcClu->GetPHOSMod() == PpsdClu->GetPHOSMod()){ 
    if(vecPpsd.X() >= vecEmc.X() - fDelta ){ 
      if(vecPpsd.Z() >= vecEmc.Z() - fDelta ){
	AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;
	//Correct to difference in CPV and EMC position due to different distance to center.
	//we assume, that particle moves from center
	Float_t DCPV = geom->GetIPtoOuterCoverDistance();
	Float_t DEMC = geom->GetIPtoCrystalSurface() ;
	DEMC         = DEMC / DCPV ;
        vecPpsd = DEMC * vecPpsd  - vecEmc ; 
        R = vecPpsd.Mag() ;
      } // if  zPpsd >= zEmc - fDelta
      TooFar = kFALSE ;
    } // if  xPpsd >= xEmc - fDelta
    else 
      TooFar = kTRUE ;
  } 
  else 
    TooFar = kTRUE ;
  
  return R ;
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakeLinks(TObjArray * EmcRecPoints, TObjArray * PpsdRecPointsUp, 
				     TObjArray * PpsdRecPointsLow, TClonesArray * LinkLowArray, 
				     TClonesArray *LinkUpArray) 
{ 
  //Finds distanses (links) between all EMC and PPSD clusters, which are not further from each other than fR0 

  TIter nextEmc(EmcRecPoints) ;
  Int_t iEmcClu = 0 ; 
  
  AliPHOSPpsdRecPoint * PpsdLow ; 
  AliPHOSPpsdRecPoint * PpsdUp ;
  AliPHOSEmcRecPoint * EmcClu ;
  
  Int_t iLinkLow = 0 ;
  Int_t iLinkUp  = 0 ;
  
  while( (EmcClu = (AliPHOSEmcRecPoint*)nextEmc() ) ) {
    Bool_t TooFar ;
    TIter nextPpsdLow(PpsdRecPointsLow ) ;
    Int_t iPpsdLow = 0 ;
    
    while( (PpsdLow = (AliPHOSPpsdRecPoint*)nextPpsdLow() ) ) { 
      Float_t R = GetDistanceInPHOSPlane(EmcClu, PpsdLow, TooFar) ;
      
      if(TooFar) 
	break ;	 
      if(R < fR0) 
	new( (*LinkLowArray)[iLinkLow++]) AliPHOSLink(R, iEmcClu, iPpsdLow) ;
      
      iPpsdLow++ ;
      
    }
    
    TIter nextPpsdUp(PpsdRecPointsUp ) ;
    Int_t iPpsdUp = 0 ;
    
    while( (PpsdUp = (AliPHOSPpsdRecPoint*)nextPpsdUp() ) ) { 
      Float_t R = GetDistanceInPHOSPlane(EmcClu, PpsdUp, TooFar) ;
      
      if(TooFar)
	break ;	 
      if(R < fR0) 
	new( (*LinkUpArray)[iLinkUp++]) AliPHOSLink(R, iEmcClu, iPpsdUp) ;
      
      iPpsdUp++ ;
      
    }
    
    iEmcClu++ ; 
    
  } // while nextEmC
  
  LinkLowArray->Sort() ; //first links with smallest distances
  LinkUpArray->Sort() ;
}
    
//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakePairs(TObjArray * EmcRecPoints, TObjArray * PpsdRecPointsUp, 
				    TObjArray * PpsdRecPointsLow, TClonesArray * LinkLowArray, 
				    TClonesArray * LinkUpArray, TrackSegmentsList * trsl) 
{ // Finds the smallest links and makes pairs of PPSD and EMC clusters with smallest distance 
  TIter nextLow(LinkLowArray) ;
  TIter nextUp(LinkUpArray) ;
  
  AliPHOSLink * linkLow ;
  AliPHOSLink * linkUp ;

  AliPHOSEmcRecPoint * emc ;
  AliPHOSPpsdRecPoint * ppsdLow ;
  AliPHOSPpsdRecPoint * ppsdUp ;

  AliPHOSRecPoint * NullPointer = 0 ;

  while ( (linkLow =  (AliPHOSLink *)nextLow() ) ){
    emc = (AliPHOSEmcRecPoint *) EmcRecPoints->At(linkLow->GetEmc()) ;
    ppsdLow = (AliPHOSPpsdRecPoint *) PpsdRecPointsLow->At(linkLow->GetPpsd()) ;
    if((emc)&&(ppsdLow)){ // RecPoints not removed yet 
	 ppsdUp = 0 ;
	 
	 while ( (linkUp =  (AliPHOSLink *)nextUp() ) ){  
	   if(linkLow->GetEmc() == linkUp->GetEmc() ){
	     ppsdUp = (AliPHOSPpsdRecPoint *) PpsdRecPointsUp->At(linkUp->GetPpsd()) ;
	     break ;
	   }
	
	 } // while nextUp
	 
	 nextUp.Reset();
         AliPHOSTrackSegment * subtr = new AliPHOSTrackSegment(emc, ppsdUp, ppsdLow ) ;
	 trsl->Add(subtr) ;  
	 EmcRecPoints->AddAt(NullPointer,linkLow->GetEmc()) ;	  
	 PpsdRecPointsLow->AddAt(NullPointer,linkLow->GetPpsd()) ;
	 
	 if(ppsdUp)  
	   PpsdRecPointsUp->AddAt(NullPointer,linkUp->GetPpsd()) ;
	 
    } // if NLocMax
  } 
   
  TIter nextEmc(EmcRecPoints) ;          
  nextEmc.Reset() ;

  while( (emc = (AliPHOSEmcRecPoint*)nextEmc()) ){ //to create pairs if no PpsdLow
    ppsdLow = 0 ; 
    ppsdUp  = 0 ;
    
    while ( (linkUp =  (AliPHOSLink *)nextUp() ) ){
      
      if(EmcRecPoints->IndexOf(emc) == linkUp->GetEmc() ){
	ppsdUp = (AliPHOSPpsdRecPoint *) PpsdRecPointsUp->At(linkUp->GetPpsd()) ;
	break ;
      }
      
    }
    AliPHOSTrackSegment * subtr = new AliPHOSTrackSegment(emc, ppsdUp, ppsdLow ) ;
    trsl->Add(subtr) ;   
 
    if(ppsdUp)  
      PpsdRecPointsUp->AddAt(NullPointer,linkUp->GetPpsd()) ;
  }
     
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakeTrackSegments(DigitsList * DL, RecPointsList * emcl, 
					RecPointsList * ppsdl, TrackSegmentsList * trsl)
{
  // main function, does the job

  Int_t PHOSMod      = 1 ;
  Int_t emcStopedAt  = 0 ; 
  Int_t ppsdStopedAt = 0 ; 
  
  TObjArray * EmcRecPoints     = new TObjArray(100) ;  // these arrays keep pointers 
  TObjArray * PpsdRecPointsUp  = new TObjArray(100) ;  // to RecPoints, which are 
  TObjArray * PpsdRecPointsLow = new TObjArray(100) ;  // kept in TClonesArray's emcl and ppsdl
  
  
  TClonesArray * LinkLowArray = new TClonesArray("AliPHOSLink", 100);
  TClonesArray * LinkUpArray  = new TClonesArray("AliPHOSLink", 100); 
  
  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;
  
  while(PHOSMod <= geom->GetNModules() ){

    cout << PHOSMod << " Track1 " << endl ;
    
    FillOneModule(DL, emcl, EmcRecPoints, ppsdl, PpsdRecPointsUp, PpsdRecPointsLow, PHOSMod , emcStopedAt, ppsdStopedAt) ;
    cout << PHOSMod << " Track2 " << endl ;
    MakeLinks(EmcRecPoints, PpsdRecPointsUp, PpsdRecPointsLow, LinkLowArray, LinkUpArray) ; 

    cout << PHOSMod << " Track3 " << endl ;
    MakePairs(EmcRecPoints, PpsdRecPointsUp, PpsdRecPointsLow, LinkLowArray, LinkUpArray, trsl) ;
 
    EmcRecPoints->Clear() ;
 
    PpsdRecPointsUp->Clear() ;
  
    PpsdRecPointsLow->Clear() ;

    LinkUpArray->Clear();
   
    LinkLowArray->Clear();
   
  }
  delete EmcRecPoints ; 
  EmcRecPoints = 0 ; 

  delete PpsdRecPointsUp ; 
  PpsdRecPointsUp = 0 ; 

  delete PpsdRecPointsLow ; 
  PpsdRecPointsLow = 0 ; 

  delete LinkUpArray ; 
  LinkUpArray = 0  ; 

  delete LinkLowArray ; 
  LinkLowArray = 0 ; 
}

//____________________________________________________________________________
Double_t  AliPHOSTrackSegmentMakerv1::ShowerShape(Double_t r)
{ 
// If you change this function, change also gradiend evaluation  in ChiSquare()
  Double_t r4 = r*r*r*r ;
  Double_t r295 = TMath::Power(r, 2.95) ;
  Double_t shape = TMath::Exp( -r4 * (1. / (2.32 + 0.26 * r4) + 0.0316 / (1 + 0.0652 * r295) ) ) ;
  return shape ;
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::UnfoldClusters(DigitsList * DL, RecPointsList * emcIn,  AliPHOSEmcRecPoint * iniEmc, 
					 Int_t Nmax, int * maxAt, Float_t * maxAtEnergy, TObjArray * emcList)
{ 
  // fits cluster with Nmax overlapping showers 
  
  Int_t NPar = 3 * Nmax ;
  Float_t FitParameters[NPar] ;
  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;

  Bool_t rv = FindFit(iniEmc, maxAt, maxAtEnergy, NPar, FitParameters) ;
  if( !rv )  // Fit failed, return and remove cluster
    return ;
  
  Float_t xDigit ;
  Float_t zDigit ;
  Int_t RelId[4] ;

  Int_t Ndigits = iniEmc->GetMultiplicity() ;  
  Float_t xpar  ;
  Float_t zpar  ;
  Float_t Epar  ;
  Float_t Distance ;
  Float_t Ratio ;
  Float_t Efit[Ndigits] ;
  Int_t iparam ;
  Int_t iDigit ;
  
  AliPHOSDigit * digit ;
  AliPHOSEmcRecPoint * emcRP ;  
  int * emcDigits = iniEmc->GetDigitsList() ;
  Float_t * emcEnergies = iniEmc->GetEnergiesList() ;

  Int_t iRecPoint = emcIn->GetEntries() ;

  for(iDigit = 0 ; iDigit < Ndigits ; iDigit ++){
    digit = (AliPHOSDigit *) emcDigits[iDigit];
    geom->AbsToRelNumbering(digit->GetId(), RelId) ;
    geom->RelPosInModule(RelId, xDigit, zDigit) ;
    Efit[iDigit] = 0;
    iparam = 0 ;
    
    while(iparam < NPar ){
      xpar = FitParameters[iparam] ;
      zpar = FitParameters[iparam+1] ;
      Epar = FitParameters[iparam+2] ;
      iparam += 3 ;
      Distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;
      Distance =  TMath::Sqrt(Distance) ;
      Efit[iDigit] += Epar * ShowerShape(Distance) ;
    }

  }

  iparam = 0 ;
  Float_t eDigit ;

  while(iparam < NPar ){
    xpar = FitParameters[iparam] ;
    zpar = FitParameters[iparam+1] ;
    Epar = FitParameters[iparam+2] ;
    iparam += 3 ;
    new ((*emcIn)[iRecPoint]) AliPHOSEmcRecPoint( iniEmc->GetLogWeightCut(), iniEmc->GetLocMaxCut() ) ;
    emcRP = (AliPHOSEmcRecPoint *) emcIn->At(iRecPoint++);

    for(iDigit = 0 ; iDigit < Ndigits ; iDigit ++){
      digit = (AliPHOSDigit *) emcDigits[iDigit];
      geom->AbsToRelNumbering(digit->GetId(), RelId) ;
      geom->RelPosInModule(RelId, xDigit, zDigit) ;
      Distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;
      Distance =  TMath::Sqrt(Distance) ;
      Ratio = Epar * ShowerShape(Distance) / Efit[iDigit] ; 
      eDigit = emcEnergies[iDigit] * Ratio ;
      emcRP->AddDigit( *digit, eDigit ) ;
    }

    emcList->Add(emcRP) ;

  }
}
//______________________________________________________________________________
void UnfoldingChiSquare(Int_t &NPar, Double_t *Grad, Double_t & fret, Double_t *x, Int_t iflag)
{

// Number of paramters, Gradient , Chi squared, parameters, what to do

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;

  AliPHOSEmcRecPoint * emcRP = (AliPHOSEmcRecPoint *) gMinuit->GetObjectFit() ; // EmcRecPoint to fit
  int * emcDigits = emcRP->GetDigitsList() ;
  Float_t * emcEnergies = emcRP->GetEnergiesList() ;
  fret = 0. ;     
  Int_t iparam ;

  if(iflag == 2)
    for(iparam = 0 ; iparam < NPar ; iparam++)    
      Grad[iparam] = 0 ; // Will evaluate gradient
  
  Double_t Efit ;  
  
  AliPHOSDigit * digit ;
  Int_t iDigit = 0 ;

  while ( (digit = (AliPHOSDigit *)emcDigits[iDigit] )){
    Int_t RelId[4] ;
    Float_t xDigit ;
    Float_t zDigit ;
    geom->AbsToRelNumbering(digit->GetId(), RelId) ;
    geom->RelPosInModule(RelId, xDigit, zDigit) ;
    
     if(iflag == 2){  // calculate gradient
       Int_t iParam = 0 ;
       Efit = 0 ;
       while(iParam < NPar ){
	 Double_t Distance = (xDigit - x[iParam]) * (xDigit - x[iParam]) ;
	 iParam++ ; 
	 Distance += (zDigit - x[iParam]) * (zDigit - x[iParam]) ; 
	 Distance = TMath::Sqrt( Distance ) ; 
	 iParam++ ; 	 
	 Efit += x[iParam] * AliPHOSTrackSegmentMakerv1::ShowerShape(Distance) ;
	 iParam++ ;
       }
       Double_t sum = 2. * (Efit - emcEnergies[iDigit]) / emcEnergies[iDigit] ; // Here we assume, that sigma = sqrt(E) 
       iParam = 0 ;
       while(iParam < NPar ){
	 Double_t xpar = x[iParam] ;
	 Double_t zpar = x[iParam+1] ;
	 Double_t Epar = x[iParam+2] ;
	 Double_t dr = TMath::Sqrt( (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar) );
	 Double_t shape = sum * AliPHOSTrackSegmentMakerv1::ShowerShape(dr) ;
	 Double_t r4 = dr*dr*dr*dr ;
	 Double_t r295 = TMath::Power(dr,2.95) ;
	 Double_t deriv =-4. * dr*dr * ( 2.32 / ( (2.32 + 0.26 * r4) * (2.32 + 0.26 * r4) ) +
					 0.0316 * (1. + 0.0171 * r295) / ( ( 1. + 0.0652 * r295) * (1. + 0.0652 * r295) ) ) ;
	 
	 Grad[iParam] += Epar * shape * deriv * (xpar - xDigit) ;  // Derivative over x    
	 iParam++ ; 
	 Grad[iParam] += Epar * shape * deriv * (zpar - zDigit) ;  // Derivative over z         
	 iParam++ ; 
	 Grad[iParam] += shape ;                                  // Derivative over energy     	
	 iParam++ ; 
       }
     }
     Efit = 0;
     iparam = 0 ;
     while(iparam < NPar ){
       Double_t xpar = x[iparam] ;
       Double_t zpar = x[iparam+1] ;
       Double_t Epar = x[iparam+2] ;
       iparam += 3 ;
       Double_t Distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;
       Distance =  TMath::Sqrt(Distance) ;
       Efit += Epar * AliPHOSTrackSegmentMakerv1::ShowerShape(Distance) ;
     }
     fret += (Efit-emcEnergies[iDigit])*(Efit-emcEnergies[iDigit])/emcEnergies[iDigit] ; 
     // Here we assume, that sigma = sqrt(E)
     iDigit++ ;
  }
}
