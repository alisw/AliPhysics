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
// Implementation version 1 of algorithm class to construct PHOS track segments
// Associates EMC and PPSD clusters
// Unfolds the EMC cluster   
//                  
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH)
//

// --- ROOT system ---

#include "TObjArray.h"
#include "TClonesArray.h"
#include "TObjectTable.h"

// --- Standard library ---

#include <iostream.h>

// --- AliRoot header files ---

#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSLink.h"
#include "AliPHOSv0.h"
#include "AliRun.h"

extern void UnfoldingChiSquare(Int_t &nPar, Double_t *Grad, Double_t & fret, Double_t *x, Int_t iflag) ; 

ClassImp( AliPHOSTrackSegmentMakerv1) 


//____________________________________________________________________________
 AliPHOSTrackSegmentMakerv1::AliPHOSTrackSegmentMakerv1() : AliPHOSTrackSegmentMaker()
{
  // ctor

  fR0 = 4. ;   
  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;
  //clusters are sorted in "rows" and "columns" of width geom->GetCrystalSize(0),
  fDelta = fR0 + geom->GetCrystalSize(0) ;
  fMinuit = new TMinuit(100) ;
  fUnfoldFlag = kTRUE ; 
}

//____________________________________________________________________________
 AliPHOSTrackSegmentMakerv1::~AliPHOSTrackSegmentMakerv1()
{ 
  // dtor
 
   delete fMinuit ; 
}

//____________________________________________________________________________
Bool_t  AliPHOSTrackSegmentMakerv1::FindFit(AliPHOSEmcRecPoint * emcRP, int * maxAt, Float_t * maxAtEnergy,
				    Int_t nPar, Float_t * fitparameters)
{ 
  // Calls TMinuit to fit the energy distribution of a cluster with several maxima 

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;

  gMinuit->SetPrintLevel(-1) ;           // No Printout
  gMinuit->SetFCN(UnfoldingChiSquare) ;  // To set the address of the minimization function 
  gMinuit->SetObjectFit(emcRP) ;         // To tranfer pointer to UnfoldingChiSquare

  // filling initial values for fit parameters
  AliPHOSDigit * digit ;

  Int_t ierflg  = 0; 
  Int_t index   = 0 ;
  Int_t nDigits = (Int_t) nPar / 3 ;

  Int_t iDigit ;


  for(iDigit = 0; iDigit < nDigits; iDigit++){
    digit = (AliPHOSDigit *) maxAt[iDigit]; 

    Int_t relid[4] ;
    Float_t x ;
    Float_t z ;
    geom->AbsToRelNumbering(digit->GetId(), relid) ;
    geom->RelPosInModule(relid, x, z) ;

    Float_t energy = maxAtEnergy[iDigit] ;

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
    gMinuit->mnparm(index, "Energy",  energy , 0.05*energy, 0., 4.*energy, ierflg) ;
    index++ ;   
    if(ierflg != 0){
      cout << "PHOS Unfolding>  Unable to set initial value for fit procedure : energy = " << energy << endl ;      
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
  for(index = 0; index < nPar; index++){
    Double_t err ;
    Double_t val ;
    gMinuit->GetParameter(index, val, err) ;    // Returns value and error of parameter index
    fitparameters[index] = val ;
   }
  return kTRUE;

}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::FillOneModule(DigitsList * dl, RecPointsList * emcIn, TObjArray * emcOut, 
					RecPointsList * ppsdIn, TObjArray * ppsdOutUp,
					TObjArray * ppsdOutLow, Int_t & phosmod, Int_t & emcStopedAt, 
					Int_t & ppsdStopedAt)
{
  // Unfold clusters and fill xxxOut arrays with clusters from one PHOS module
 
  AliPHOSEmcRecPoint *  emcRecPoint  ; 
  AliPHOSPpsdRecPoint * ppsdRecPoint ;
  Int_t index ;
  
  Int_t nEmcUnfolded = emcIn->GetEntries() ;
  for(index = emcStopedAt; index < nEmcUnfolded; index++){
    emcRecPoint = (AliPHOSEmcRecPoint *) (*emcIn)[index] ;

    if(emcRecPoint->GetPHOSMod() != phosmod )  
       break ;
    
    Int_t nMultipl = emcRecPoint->GetMultiplicity() ; 
    Int_t * maxAt = new Int_t[nMultipl] ;
    Float_t * maxAtEnergy = new Float_t[nMultipl] ;
    Int_t nMax = emcRecPoint->GetNumberOfLocalMax(maxAt, maxAtEnergy) ;

    if(nMax <= 1 )     // if cluster is very flat (no pronounced maximum) then nMax = 0 
      emcOut->Add(emcRecPoint) ;
    else if (fUnfoldFlag) {
      UnfoldClusters(dl, emcIn, emcRecPoint, nMax, maxAt, maxAtEnergy, emcOut) ;
      emcIn->Remove(emcRecPoint); 
      emcIn->Compress() ;
      nEmcUnfolded-- ;
      index-- ;
    }
    
    delete[] maxAt ; 
    delete[] maxAtEnergy ; 
  }
  emcStopedAt = index ;

  for(index = ppsdStopedAt; index < ppsdIn->GetEntries(); index++){
    ppsdRecPoint = (AliPHOSPpsdRecPoint *) (*ppsdIn)[index] ;
    if(ppsdRecPoint->GetPHOSMod() != phosmod )   
      break ;
    if(ppsdRecPoint->GetUp() ) 
      ppsdOutUp->Add(ppsdRecPoint) ;
    else  
      ppsdOutLow->Add(ppsdRecPoint) ;
  }
  ppsdStopedAt = index ;
   
  emcOut->Sort() ;
  ppsdOutUp->Sort() ;
  ppsdOutLow->Sort() ;   
}
//____________________________________________________________________________
Float_t  AliPHOSTrackSegmentMakerv1::GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * emcclu,AliPHOSPpsdRecPoint * PpsdClu, Bool_t &toofar)
{
  // Calculates the distance between the EMC RecPoint and the PPSD RecPoint
 
  Float_t r = fR0 ;
 
  TVector3 vecEmc ;
  TVector3 vecPpsd ;
  
  emcclu->GetLocalPosition(vecEmc) ;
  PpsdClu->GetLocalPosition(vecPpsd)  ; 
  if(emcclu->GetPHOSMod() == PpsdClu->GetPHOSMod()){ 
    if(vecPpsd.X() >= vecEmc.X() - fDelta ){ 
      if(vecPpsd.Z() >= vecEmc.Z() - fDelta ){
	AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;
	// Correct to difference in CPV and EMC position due to different distance to center.
	// we assume, that particle moves from center
	Float_t dCPV = geom->GetIPtoOuterCoverDistance();
	Float_t dEMC = geom->GetIPtoCrystalSurface() ;
	dEMC         = dEMC / dCPV ;
        vecPpsd = dEMC * vecPpsd  - vecEmc ; 
        r = vecPpsd.Mag() ;
      } // if  zPpsd >= zEmc - fDelta
      toofar = kFALSE ;
    } // if  xPpsd >= xEmc - fDelta
    else 
      toofar = kTRUE ;
  } 
  else 
    toofar = kTRUE ;
  
  return r ;
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakeLinks(TObjArray * emcRecPoints, TObjArray * ppsdRecPointsUp, 
				     TObjArray * ppsdRecPointsLow, TClonesArray * linklowArray, 
				     TClonesArray *linkupArray) 
{ 
  // Finds distances (links) between all EMC and PPSD clusters, which are not further apart from each other than fR0 
  
  TIter nextEmc(emcRecPoints) ;
  Int_t iEmcClu = 0 ; 
  
  AliPHOSPpsdRecPoint * ppsdlow ; 
  AliPHOSPpsdRecPoint * ppsdup ;
  AliPHOSEmcRecPoint * emcclu ;
  
  Int_t iLinkLow = 0 ;
  Int_t iLinkUp  = 0 ;
  
  while( (emcclu = (AliPHOSEmcRecPoint*)nextEmc() ) ) {
    Bool_t toofar ;
    TIter nextPpsdLow(ppsdRecPointsLow ) ;
    Int_t iPpsdLow = 0 ;
    
    while( (ppsdlow = (AliPHOSPpsdRecPoint*)nextPpsdLow() ) ) { 
      Float_t r = GetDistanceInPHOSPlane(emcclu, ppsdlow, toofar) ;
      
      if(toofar) 
	break ;	 
      if(r < fR0) 
	new( (*linklowArray)[iLinkLow++]) AliPHOSLink(r, iEmcClu, iPpsdLow) ;
      
      iPpsdLow++ ;
      
    }
    
    TIter nextPpsdUp(ppsdRecPointsUp ) ;
    Int_t iPpsdUp = 0 ;
    
    while( (ppsdup = (AliPHOSPpsdRecPoint*)nextPpsdUp() ) ) { 
      Float_t r = GetDistanceInPHOSPlane(emcclu, ppsdup, toofar) ;
      
      if(toofar)
	break ;	 
      if(r < fR0) 
	new( (*linkupArray)[iLinkUp++]) AliPHOSLink(r, iEmcClu, iPpsdUp) ;
      
      iPpsdUp++ ;
      
    }
    
    iEmcClu++ ; 
    
  } // while nextEmC
  
  linklowArray->Sort() ; //first links with smallest distances
  linkupArray->Sort() ;
}
    
//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakePairs(TObjArray * emcRecPoints, TObjArray * ppsdRecPointsUp, 
				    TObjArray * ppsdRecPointsLow, TClonesArray * linklowArray, 
				    TClonesArray * linkupArray, TrackSegmentsList * trsl) 
{ 

  // Finds the smallest links and makes pairs of PPSD and EMC clusters with smallest distance 
  
  TIter nextLow(linklowArray) ;
  TIter nextUp(linkupArray) ;
  
  AliPHOSLink * linkLow ;
  AliPHOSLink * linkUp ;

  AliPHOSEmcRecPoint * emc ;
  AliPHOSPpsdRecPoint * ppsdLow ;
  AliPHOSPpsdRecPoint * ppsdUp ;

  AliPHOSRecPoint * nullpointer = 0 ;

  while ( (linkLow =  (AliPHOSLink *)nextLow() ) ){
    emc = (AliPHOSEmcRecPoint *) emcRecPoints->At(linkLow->GetEmc()) ;
    ppsdLow = (AliPHOSPpsdRecPoint *) ppsdRecPointsLow->At(linkLow->GetPpsd()) ;
    if( (emc) && (ppsdLow) ){ // RecPoints not removed yet 
	 ppsdUp = 0 ;
	 
	 while ( (linkUp =  (AliPHOSLink *)nextUp() ) ){  
	   if(linkLow->GetEmc() == linkUp->GetEmc() ){
	     ppsdUp = (AliPHOSPpsdRecPoint *) ppsdRecPointsUp->At(linkUp->GetPpsd()) ;
	     break ;
	   }
	
	 } // while nextUp
	 
	 nextUp.Reset();
//          AliPHOSTrackSegment * subtr = new AliPHOSTrackSegment(emc, ppsdUp, ppsdLow ) ;
// 	 trsl->Add(subtr) ;  
	 new( (*trsl)[fNTrackSegments] ) AliPHOSTrackSegment(emc, ppsdUp, ppsdLow ) ;
	 fNTrackSegments++ ;
	 emcRecPoints->AddAt(nullpointer,linkLow->GetEmc()) ;	  
	 ppsdRecPointsLow->AddAt(nullpointer,linkLow->GetPpsd()) ;
	 
	 if(ppsdUp)  
	   ppsdRecPointsUp->AddAt(nullpointer,linkUp->GetPpsd()) ;
	 
    } 
  } 
   
  TIter nextEmc(emcRecPoints) ;          
  nextEmc.Reset() ;

  while( (emc = (AliPHOSEmcRecPoint*)nextEmc()) ){ //to create pairs if no ppsdlow
    ppsdLow = 0 ; 
    ppsdUp  = 0 ;
    
    while ( (linkUp =  (AliPHOSLink *)nextUp() ) ){
      
      if(emcRecPoints->IndexOf(emc) == linkUp->GetEmc() ){
	ppsdUp = (AliPHOSPpsdRecPoint *) ppsdRecPointsUp->At(linkUp->GetPpsd()) ;
	break ;
      }
      
    }
//     AliPHOSTrackSegment * subtr = new AliPHOSTrackSegment(emc, ppsdUp, ppsdLow ) ;
//     trsl->Add(subtr) ;   
    new( (*trsl)[fNTrackSegments] ) AliPHOSTrackSegment(emc, ppsdUp, ppsdLow ) ;
    fNTrackSegments++ ;
    

    if(ppsdUp)  
      ppsdRecPointsUp->AddAt(nullpointer,linkUp->GetPpsd()) ;
  }
     
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakeTrackSegments(DigitsList * dl, RecPointsList * emcl, 
					RecPointsList * ppsdl, TrackSegmentsList * trsl)
{
  // Makes the track segments out of the list of EMC and PPSD Recpoints and stores them in a list

  Int_t phosmod      = 1 ;
  Int_t emcStopedAt  = 0 ; 
  Int_t ppsdStopedAt = 0 ; 
  
  TObjArray * emcRecPoints     = new TObjArray(100) ;  // these arrays keep pointers 
  TObjArray * ppsdRecPointsUp  = new TObjArray(100) ;  // to RecPoints, which are 
  TObjArray * ppsdRecPointsLow = new TObjArray(100) ;  // kept in TClonesArray's emcl and ppsdl
  
  
  TClonesArray * linklowArray = new TClonesArray("AliPHOSLink", 100);
  TClonesArray * linkupArray  = new TClonesArray("AliPHOSLink", 100); 
  
  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;
  
  while(phosmod <= geom->GetNModules() ){
    
    FillOneModule(dl, emcl, emcRecPoints, ppsdl, ppsdRecPointsUp, ppsdRecPointsLow, phosmod, emcStopedAt, ppsdStopedAt) ;

    MakeLinks(emcRecPoints, ppsdRecPointsUp, ppsdRecPointsLow, linklowArray, linkupArray) ; 

    MakePairs(emcRecPoints, ppsdRecPointsUp, ppsdRecPointsLow, linklowArray, linkupArray, trsl) ;
 
    emcRecPoints->Clear() ;
 
    ppsdRecPointsUp->Clear() ;
  
    ppsdRecPointsLow->Clear() ;

    linkupArray->Clear() ;
   
    linklowArray->Clear() ;
   
    phosmod++ ; 
  }
  delete emcRecPoints ; 
  emcRecPoints = 0 ; 

  delete ppsdRecPointsUp ; 
  ppsdRecPointsUp = 0 ; 

  delete ppsdRecPointsLow ; 
  ppsdRecPointsLow = 0 ; 

  delete linkupArray ; 
  linkupArray = 0  ; 

  delete linklowArray ; 
  linklowArray = 0 ; 
}

//____________________________________________________________________________
Double_t  AliPHOSTrackSegmentMakerv1::ShowerShape(Double_t r)
{ 
  // Shape of the shower (see PHOS TDR)
  // If you change this function, change also the gradien evaluation  in ChiSquare()

  Double_t r4    = r*r*r*r ;
  Double_t r295  = TMath::Power(r, 2.95) ;
  Double_t shape = TMath::Exp( -r4 * (1. / (2.32 + 0.26 * r4) + 0.0316 / (1 + 0.0652 * r295) ) ) ;
  return shape ;
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::UnfoldClusters(DigitsList * dl, RecPointsList * emcIn,  AliPHOSEmcRecPoint * iniEmc, 
					 Int_t nMax, int * maxAt, Float_t * maxAtEnergy, TObjArray * emcList)
{ 
  // Performs the unfolding of a cluster with nMax overlapping showers 
  // This is time consuming (use the (Un)SetUnfolFlag()  )

  Int_t nPar = 3 * nMax ;
  Float_t * fitparameters = new Float_t[nPar] ;
  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;

  Bool_t rv = FindFit(iniEmc, maxAt, maxAtEnergy, nPar, fitparameters) ;
  if( !rv ) {
    // Fit failed, return and remove cluster
    delete[] fitparameters ; 
    return ;
  }
  
  Float_t xDigit ;
  Float_t zDigit ;
  Int_t relid[4] ;

  Int_t nDigits = iniEmc->GetMultiplicity() ;  
  Float_t xpar  ;
  Float_t zpar  ;
  Float_t epar  ;
  Float_t distance ;
  Float_t ratio ;
  Float_t * efit = new Float_t[nDigits] ;
  Int_t iparam ;
  Int_t iDigit ;
  
  AliPHOSDigit * digit ;
  AliPHOSEmcRecPoint * emcRP ;  
  Int_t * emcDigits = iniEmc->GetDigitsList() ;
  Float_t * emcEnergies = iniEmc->GetEnergiesList() ;

  Int_t iRecPoint = emcIn->GetEntries() ;

  for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
    digit = (AliPHOSDigit *) emcDigits[iDigit];
    geom->AbsToRelNumbering(digit->GetId(), relid) ;
    geom->RelPosInModule(relid, xDigit, zDigit) ;
    efit[iDigit] = 0;
    iparam = 0 ;
    
    while(iparam < nPar ){
      xpar = fitparameters[iparam] ;
      zpar = fitparameters[iparam+1] ;
      epar = fitparameters[iparam+2] ;
      iparam += 3 ;
      distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;
      distance =  TMath::Sqrt(distance) ;
      efit[iDigit] += epar * AliPHOSTrackSegmentMakerv1::ShowerShape(distance) ;
    }
  }

  iparam = 0 ;
  Float_t eDigit ;

  while(iparam < nPar ){
    xpar = fitparameters[iparam] ;
    zpar = fitparameters[iparam+1] ;
    epar = fitparameters[iparam+2] ;
    iparam += 3 ;
    new ((*emcIn)[iRecPoint]) AliPHOSEmcRecPoint( iniEmc->GetLogWeightCut(), iniEmc->GetLocMaxCut() ) ;
    emcRP = (AliPHOSEmcRecPoint *) emcIn->At(iRecPoint++);

    for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
      digit = (AliPHOSDigit *) emcDigits[iDigit];
      geom->AbsToRelNumbering(digit->GetId(), relid) ;
      geom->RelPosInModule(relid, xDigit, zDigit) ;
      distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;
      distance =  TMath::Sqrt(distance) ;
      ratio = epar * AliPHOSTrackSegmentMakerv1::ShowerShape(distance) / efit[iDigit] ; 
      eDigit = emcEnergies[iDigit] * ratio ;
      emcRP->AddDigit( *digit, eDigit ) ;
    }

    emcList->Add(emcRP) ;

  }
  
  delete[] fitparameters ; 
  delete[] efit ; 

}

//______________________________________________________________________________
void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)
{
  // Calculates th Chi square for the cluster unfolding minimization
  // Number of parameters, Gradient, Chi squared, parameters, what to do

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;

  AliPHOSEmcRecPoint * emcRP = (AliPHOSEmcRecPoint *) gMinuit->GetObjectFit() ; // EmcRecPoint to fit
  Int_t * emcDigits     = emcRP->GetDigitsList() ;
  Float_t * emcEnergies = emcRP->GetEnergiesList() ;
  fret = 0. ;     
  Int_t iparam ;

  if(iflag == 2)
    for(iparam = 0 ; iparam < nPar ; iparam++)    
      Grad[iparam] = 0 ; // Will evaluate gradient
  
  Double_t efit ;  
  
  AliPHOSDigit * digit ;
  Int_t iDigit = 0 ;

  while ( (digit = (AliPHOSDigit *)emcDigits[iDigit] )){
    Int_t relid[4] ;
    Float_t xDigit ;
    Float_t zDigit ;
    geom->AbsToRelNumbering(digit->GetId(), relid) ;
    geom->RelPosInModule(relid, xDigit, zDigit) ;
    
     if(iflag == 2){  // calculate gradient
       Int_t iParam = 0 ;
       efit = 0 ;
       while(iParam < nPar ){
	 Double_t distance = (xDigit - x[iParam]) * (xDigit - x[iParam]) ;
	 iParam++ ; 
	 distance += (zDigit - x[iParam]) * (zDigit - x[iParam]) ; 
	 distance = TMath::Sqrt( distance ) ; 
	 iParam++ ; 	 
	 efit += x[iParam] * AliPHOSTrackSegmentMakerv1::ShowerShape(distance) ;
	 iParam++ ;
       }
       Double_t sum = 2. * (efit - emcEnergies[iDigit]) / emcEnergies[iDigit] ; // Here we assume, that sigma = sqrt(E) 
       iParam = 0 ;
       while(iParam < nPar ){
	 Double_t xpar = x[iParam] ;
	 Double_t zpar = x[iParam+1] ;
	 Double_t epar = x[iParam+2] ;
	 Double_t dr = TMath::Sqrt( (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar) );
	 Double_t shape = sum * AliPHOSTrackSegmentMakerv1::ShowerShape(dr) ;
	 Double_t r4 = dr*dr*dr*dr ;
	 Double_t r295 = TMath::Power(dr,2.95) ;
	 Double_t deriv =-4. * dr*dr * ( 2.32 / ( (2.32 + 0.26 * r4) * (2.32 + 0.26 * r4) ) +
					 0.0316 * (1. + 0.0171 * r295) / ( ( 1. + 0.0652 * r295) * (1. + 0.0652 * r295) ) ) ;
	 
	 Grad[iParam] += epar * shape * deriv * (xpar - xDigit) ;  // Derivative over x    
	 iParam++ ; 
	 Grad[iParam] += epar * shape * deriv * (zpar - zDigit) ;  // Derivative over z         
	 iParam++ ; 
	 Grad[iParam] += shape ;                                  // Derivative over energy     	
	 iParam++ ; 
       }
     }
     efit = 0;
     iparam = 0 ;
     while(iparam < nPar ){
       Double_t xpar = x[iparam] ;
       Double_t zpar = x[iparam+1] ;
       Double_t epar = x[iparam+2] ;
       iparam += 3 ;
       Double_t distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;
       distance =  TMath::Sqrt(distance) ;
       efit += epar * AliPHOSTrackSegmentMakerv1::ShowerShape(distance) ;
     }
     fret += (efit-emcEnergies[iDigit])*(efit-emcEnergies[iDigit])/emcEnergies[iDigit] ; 
     // Here we assume, that sigma = sqrt(E)
     iDigit++ ;
  }
}
