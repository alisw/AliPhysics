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
#include "AliPHOSIndexToObject.h"
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

  fR0 = 10. ;   
  //clusters are sorted in "rows" and "columns" of width geom->GetCrystalSize(0),
  fDelta = fR0 + fGeom->GetCrystalSize(0) ;
  if(!gMinuit) gMinuit = new TMinuit(100) ;
  fUnfoldFlag = kTRUE ; 
}

//____________________________________________________________________________
 AliPHOSTrackSegmentMakerv1::~AliPHOSTrackSegmentMakerv1()
{ 
  // dtor

  delete gMinuit ; 
  gMinuit = 0 ;
}

//____________________________________________________________________________
Bool_t  AliPHOSTrackSegmentMakerv1::FindFit(AliPHOSEmcRecPoint * emcRP, int * maxAt, Float_t * maxAtEnergy,
				    Int_t nPar, Float_t * fitparameters)
{ 
  // Calls TMinuit to fit the energy distribution of a cluster with several maxima 

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
    fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
    fGeom->RelPosInModule(relid, x, z) ;

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

  gMinuit->mnexcm("SET STR", &p2, 0, ierflg) ;   // force TMinuit to reduce function calls  
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
void  AliPHOSTrackSegmentMakerv1::FillOneModule(AliPHOSRecPoint::RecPointsList * emcIn, 
						TArrayI * emcOut, 
						AliPHOSRecPoint::RecPointsList * ppsdIn, 
						TArrayI * ppsdOutUp,
						TArrayI * ppsdOutLow, 
						Int_t & phosmod, 
						Int_t & emcStopedAt, 
						Int_t & ppsdStopedAt)
{
  // Fill xxxOut arrays with clusters from one PHOS module
 
  AliPHOSEmcRecPoint *  emcRecPoint  ; 
  AliPHOSPpsdRecPoint * ppsdRecPoint ;
  Int_t index ;
  
  Int_t nEmcUnfolded = emcIn->GetEntries() ;
  emcOut->Set(nEmcUnfolded);
  Int_t inEmcOut = 0 ;
  for(index = emcStopedAt; index < nEmcUnfolded; index++){

    emcRecPoint = (AliPHOSEmcRecPoint *) emcIn->At(index) ;
    
    if(emcRecPoint->GetPHOSMod() != phosmod )  
      break ;
    
    emcOut->AddAt(emcRecPoint->GetIndexInList(),inEmcOut) ;
    inEmcOut++ ; 
  }
  emcOut->Set(inEmcOut) ;

  emcStopedAt = index ;

  ppsdOutLow->Set(ppsdIn->GetEntries()) ;
  ppsdOutUp->Set(ppsdIn->GetEntries()) ;
  Int_t inPpsdLow = 0;
  Int_t inPpsdUp = 0;
  for(index = ppsdStopedAt; index < ppsdIn->GetEntries(); index++){
    ppsdRecPoint = (AliPHOSPpsdRecPoint *) ppsdIn->At(index) ;
    if(ppsdRecPoint->GetPHOSMod() != phosmod )   
      break ;
    if(phosmod <= fGeom->GetNCPVModules())   //in CPV
      ppsdOutUp->AddAt(index,inPpsdUp++) ;
    else{			             //in PPSD
      if(ppsdRecPoint->GetUp() ) 
	ppsdOutUp->AddAt(index,inPpsdUp++) ;
      else  
	ppsdOutLow->AddAt(index,inPpsdLow++) ;
    }
  }
  ppsdOutLow->Set(inPpsdLow);
  ppsdOutUp->Set(inPpsdUp);
  ppsdStopedAt = index ;
   
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
    //    if(vecPpsd.X() >= vecEmc.X() - fDelta ){ 
    //  if(vecPpsd.Z() >= vecEmc.Z() - fDelta ){
	// Correct to difference in CPV and EMC position due to different distance to center.
	// we assume, that particle moves from center
	Float_t dCPV = fGeom->GetIPtoOuterCoverDistance();
	Float_t dEMC = fGeom->GetIPtoCrystalSurface() ;
	dEMC         = dEMC / dCPV ;
        vecPpsd = dEMC * vecPpsd  - vecEmc ; 
        r = vecPpsd.Mag() ;
	//    } // if  zPpsd >= zEmc - fDelta
      toofar = kFALSE ;
      //} // if  xPpsd >= xEmc - fDelta
      // else 
      //toofar = kTRUE ;
  } 
  else 
    toofar = kTRUE ;

  //toofar = kFALSE ;
 
  
  return r ;
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakeLinks(TArrayI * emcRecPoints, TArrayI * ppsdRecPointsUp, 
				     TArrayI * ppsdRecPointsLow, TClonesArray * linklowArray, 
				     TClonesArray *linkupArray) 
{ 
  // Finds distances (links) between all EMC and PPSD clusters, which are not further apart from each other than fR0 
  
  
  AliPHOSPpsdRecPoint * ppsdlow ; 
  AliPHOSPpsdRecPoint * ppsdup ;
  AliPHOSEmcRecPoint * emcclu ;

  Int_t iLinkLow = 0 ;
  Int_t iLinkUp  = 0 ;
  
  Int_t iEmcRP;
  
  for(iEmcRP = 0; iEmcRP < emcRecPoints->GetSize(); iEmcRP++ ) {
    emcclu = (AliPHOSEmcRecPoint *) fPlease->GimeRecPoint(emcRecPoints->At(iEmcRP),"emc") ;
    Bool_t toofar ;
    
    Int_t iPpsdLow ;
    
    for(iPpsdLow = 0; iPpsdLow < ppsdRecPointsLow->GetSize();iPpsdLow++ ) {
      
      ppsdlow = (AliPHOSPpsdRecPoint *) fPlease->GimeRecPoint(ppsdRecPointsLow->At(iPpsdLow),"ppsd") ;
      Float_t r = GetDistanceInPHOSPlane(emcclu, ppsdlow, toofar) ;
      
      if(toofar) 
	break ;	 
      if(r < fR0){
	new( (*linklowArray)[iLinkLow++]) AliPHOSLink(r, iEmcRP, iPpsdLow) ;
      }
    }
    
    Int_t iPpsdUp = 0 ;    
    for(iPpsdUp = 0; iPpsdUp < ppsdRecPointsUp->GetSize();iPpsdUp++ ) { 
      
      ppsdup = (AliPHOSPpsdRecPoint *)fPlease->GimeRecPoint(ppsdRecPointsUp->At(iPpsdUp),"ppsd") ;
      Float_t r = GetDistanceInPHOSPlane(emcclu, ppsdup, toofar) ;
      
      if(toofar)
	break ;	 
      if(r < fR0) { 
	new( (*linkupArray)[iLinkUp++]) AliPHOSLink(r, iEmcRP, iPpsdUp) ;
      }      
    }
  } 
  
  linklowArray->Sort() ; //first links with smallest distances
  linkupArray->Sort() ;
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakePairs(TArrayI * emcRecPoints, 
					    TArrayI * ppsdRecPointsUp, 
					    TArrayI * ppsdRecPointsLow, 
					    TClonesArray * linklowArray, 
					    TClonesArray * linkupArray, 
					    AliPHOSTrackSegment::TrackSegmentsList * trsl) 
{ 

  // Finds the smallest links and makes pairs of PPSD and EMC clusters with smallest distance 
  
  TIter nextLow(linklowArray) ;
  TIter nextUp(linkupArray) ;
  
  AliPHOSLink * linkLow ;
  AliPHOSLink * linkUp ;

  Int_t emc ;
  Int_t ppsdLow ;
  Int_t ppsdUp ;

  AliPHOSPpsdRecPoint * nullpointer = 0 ;
  ppsdUp = 0 ;

  while ( (linkLow =  (AliPHOSLink *)nextLow() ) ){
  
    emc = emcRecPoints->At(linkLow->GetEmc()) ;
    ppsdLow = ppsdRecPointsLow->At(linkLow->GetPpsd()) ;

    if( (emc >= 0) && (ppsdLow >= 0) ){    // RecPoints not removed yet 
      
      new( (*trsl)[fNTrackSegments] ) AliPHOSTrackSegment((AliPHOSEmcRecPoint *)fPlease->GimeRecPoint(emc,"emc"), 
							  nullpointer, 
							  (AliPHOSPpsdRecPoint *)fPlease->GimeRecPoint(ppsdLow,"ppsd") ) ;
      ((AliPHOSTrackSegment* )trsl->At(fNTrackSegments))->SetIndexInList(fNTrackSegments);    
      //replace index of emc to negative and shifted index of TS      
      emcRecPoints->AddAt(-2 - fNTrackSegments,linkLow->GetEmc()) ;  
      //replace index of PPSD Low to negative and shifted index of TS      
      ppsdRecPointsLow->AddAt(-2 - fNTrackSegments,linkLow->GetPpsd()) ; 
      fNTrackSegments++ ;

    } 
  } 
	 
  while ( (linkUp =  (AliPHOSLink *)nextUp() ) ){  
    emc = emcRecPoints->At(linkUp->GetEmc()) ;
    if(emc != -1){ //without ppsd Up yet 

      ppsdUp = ppsdRecPointsUp->At(linkUp->GetPpsd()) ;
      if(ppsdUp >= 0){ //ppsdUp still exist
	
	if(emc >= 0){ //without ppsd Low => create new TS

	  fNTrackSegments = trsl->GetEntries() ; 
	  new( (*trsl)[fNTrackSegments] ) AliPHOSTrackSegment((AliPHOSEmcRecPoint *) fPlease->GimeRecPoint(emc,"emc"), 
							      (AliPHOSPpsdRecPoint *)fPlease->GimeRecPoint(ppsdUp,"ppsd"), 
							      nullpointer) ;
	  ((AliPHOSTrackSegment *) trsl->At(fNTrackSegments))->SetIndexInList(fNTrackSegments);
	  fNTrackSegments++ ;
	}
	else{ // append ppsd Up to existing TS
	  ((AliPHOSTrackSegment *)trsl->At(-2-emc))->SetPpsdUpRecPoint((AliPHOSPpsdRecPoint *)fPlease->GimeRecPoint(ppsdUp,"ppsd"));
	}

	emcRecPoints->AddAt(-1,linkUp->GetEmc()) ; //Mark that PPSD Up found 
	//replace index of PPSD Up to negative and shifted index of TS      
	ppsdRecPointsUp->AddAt(-2 - fNTrackSegments,linkUp->GetPpsd()) ; 
      } //if ppsdUp still exist
    } 
  }	 

  Int_t iEmcRP ;
  for(iEmcRP = 0; iEmcRP <emcRecPoints->GetSize() ; iEmcRP++ ){
    emc = emcRecPoints->At(iEmcRP) ;
    if(emc >=0 ){
      ppsdUp = 0;
      ppsdLow = 0;
      new( (*trsl)[fNTrackSegments] ) AliPHOSTrackSegment((AliPHOSEmcRecPoint *) fPlease->GimeRecPoint(emc,"emc"), 
							  nullpointer, nullpointer ) ;
      ((AliPHOSTrackSegment *) trsl->At(fNTrackSegments))->SetIndexInList(fNTrackSegments);
      fNTrackSegments++;    
    }
    
  }
  
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakeTrackSegments(DigitsList * dl, 
						    AliPHOSRecPoint::RecPointsList * emcl, 
						    AliPHOSRecPoint::RecPointsList * ppsdl, 
						    AliPHOSTrackSegment::TrackSegmentsList * trsl)
{
  // Makes the track segments out of the list of EMC and PPSD Recpoints and stores them in a list
  
  Int_t emcStopedAt  = 0 ; 
  Int_t ppsdStopedAt = 0 ; 

  fNTrackSegments = 0 ; 
  
  TArrayI * emcRecPoints     = new TArrayI(1000) ;  // these arrays keep indexes 
  TArrayI * ppsdRecPointsUp  = new TArrayI(1000) ;  // of RecPoints, which are 
  TArrayI * ppsdRecPointsLow = new TArrayI(1000) ;  // kept in TClonesArray's emcl, ppsdl, cpv
  
  TClonesArray * linklowArray = new TClonesArray("AliPHOSLink", 1000);
  TClonesArray * linkupArray  = new TClonesArray("AliPHOSLink", 1000); 

  if(fUnfoldFlag){
    UnfoldAll(dl, emcl) ; // Unfolds all EMC clusters
    UnfoldAll(dl, ppsdl) ; // Unfolds all CPV clusters
  }

  Int_t phosmod  = 1 ;
  while(phosmod <= fGeom->GetNModules() ){
    
    FillOneModule(emcl, emcRecPoints, ppsdl, ppsdRecPointsUp, ppsdRecPointsLow, phosmod, emcStopedAt, ppsdStopedAt) ;
    
    MakeLinks(emcRecPoints, ppsdRecPointsUp, ppsdRecPointsLow, linklowArray, linkupArray) ; 
    
    MakePairs(emcRecPoints, ppsdRecPointsUp, ppsdRecPointsLow, linklowArray, linkupArray, trsl) ;
    
    emcRecPoints->Reset() ;
    
    ppsdRecPointsUp->Reset() ;
    
    ppsdRecPointsLow->Reset() ;
    
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
void  AliPHOSTrackSegmentMakerv1::MakeTrackSegmentsCPV(DigitsList * dl, 
                                                       AliPHOSRecPoint::RecPointsList * emcl, 
                                                       AliPHOSRecPoint::RecPointsList * cpvl)
{
  // Unfold clusters in EMC and CPV and refill reconstructed point lists emcl and ppsdl
  // Yuri Kharlov. 19 October 2000
  
  fNTrackSegments = 0 ; 

  TArrayI * emcRecPoints     = new TArrayI(1000) ;  // these arrays keep indexes 
  TArrayI * cpvRecPoints     = new TArrayI(1000) ;  // of RecPoints, which are kept in emcl and ppsdl
  
  if(fUnfoldFlag){
    UnfoldAll(dl, emcl) ;   // Unfolds all EMC clusters
    UnfoldAll(dl, cpvl) ;   // Unfolds all CPV clusters
  }

//    Int_t phosmod      = 1 ;
//    Int_t emcStopedAt  = 0 ; 
//    Int_t cpvStopedAt  = 0 ; 
//    while(phosmod <= fGeom->GetNModules() ){
//      FillOneModule(emcl, emcRecPoints, ppsdl, cpvRecPoints, phosmod, emcStopedAt, cpvStopedAt) ;
//      emcRecPoints->Reset() ;
//      cpvRecPoints->Reset() ;
//      phosmod++ ; 
//    }

  delete emcRecPoints ; emcRecPoints = 0 ; 
  delete cpvRecPoints ; cpvRecPoints = 0 ; 
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
void  AliPHOSTrackSegmentMakerv1::UnfoldAll(DigitsList * dl, AliPHOSRecPoint::RecPointsList * emcIn) 
{
  // Performs unfolding of all EMC/CPV but NOT ppsd clusters, sorts them and resets indexes in RecPoints

  AliPHOSEmcRecPoint *  emcRecPoint  ; 
  Int_t index ;
  Int_t nEmcUnfolded = emcIn->GetEntries() ;

  Int_t nModulesToUnfold ;

  if(emcIn->GetEntries() > 0){

    if(((AliPHOSRecPoint *)emcIn->At(0))->IsEmc())
      nModulesToUnfold = fGeom->GetNModules() ; 
    else
      nModulesToUnfold = fGeom->GetNCPVModules() ;
    
    for(index = 0 ; index < nEmcUnfolded; index++){
      
      emcRecPoint = (AliPHOSEmcRecPoint *) emcIn->At(index) ;
      if(emcRecPoint->GetPHOSMod()> nModulesToUnfold)
	break ;
      
      Int_t nMultipl = emcRecPoint->GetMultiplicity() ; 
      Int_t * maxAt = new Int_t[nMultipl] ;
      Float_t * maxAtEnergy = new Float_t[nMultipl] ;
      Int_t nMax = emcRecPoint->GetNumberOfLocalMax(maxAt, maxAtEnergy) ;
      
      if( nMax > 1 ) {     // if cluster is very flat (no pronounced maximum) then nMax = 0       
	UnfoldClusters(dl, emcIn, emcRecPoint, nMax, maxAt, maxAtEnergy) ;
	emcIn->Remove(emcRecPoint); 
	emcIn->Compress() ;
	index-- ;
	nEmcUnfolded-- ;
      }
      
      delete[] maxAt ; 
      delete[] maxAtEnergy ; 
    } //Unfolding finished
    
    emcIn->Sort() ;
    
    // to set index to new and correct index of old RecPoints
    for( index = 0 ; index < emcIn->GetEntries() ; index++){
      
      ((AliPHOSEmcRecPoint *) emcIn->At(index))->SetIndexInList(index) ;   
      
    }
    
    emcIn->Sort() ;
  }

}
//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::UnfoldClusters(DigitsList * dl, 
						 AliPHOSRecPoint::RecPointsList * emcIn, 
						 AliPHOSEmcRecPoint * iniEmc, 
						 Int_t nMax, 
						 int * maxAt, 
						 Float_t * maxAtEnergy)
{ 
  // Performs the unfolding of a cluster with nMax overlapping showers 
  // This is time consuming (use the (Un)SetUnfolFlag()  )

  Int_t nPar = 3 * nMax ;
  Float_t * fitparameters = new Float_t[nPar] ;


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
    digit = fPlease->GimeDigit( emcDigits[iDigit] ) ;   
    fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
    fGeom->RelPosInModule(relid, xDigit, zDigit) ;
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

    if(iRecPoint >= emcIn->GetSize())
      emcIn->Expand(2*iRecPoint) ;
    (*emcIn)[iRecPoint] = new AliPHOSEmcRecPoint( iniEmc->GetLogWeightCut(), iniEmc->GetLocMaxCut() ) ;
    
    emcRP = (AliPHOSEmcRecPoint *) emcIn->At(iRecPoint);
    iRecPoint++ ;

    for(iDigit = 0 ; iDigit < nDigits ; iDigit ++){
      digit = fPlease->GimeDigit( emcDigits[iDigit] ) ; 
      fGeom->AbsToRelNumbering(digit->GetId(), relid) ;
      fGeom->RelPosInModule(relid, xDigit, zDigit) ;
      distance = (xDigit - xpar) * (xDigit - xpar) + (zDigit - zpar) * (zDigit - zpar)  ;
      distance =  TMath::Sqrt(distance) ;
      ratio = epar * AliPHOSTrackSegmentMakerv1::ShowerShape(distance) / efit[iDigit] ; 
      eDigit = emcEnergies[iDigit] * ratio ;
      emcRP->AddDigit( *digit, eDigit ) ;
    }

  }
  
  delete[] fitparameters ; 
  delete[] efit ; 

}

//______________________________________________________________________________
void UnfoldingChiSquare(Int_t & nPar, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)
{
  // Calculates th Chi square for the cluster unfolding minimization
  // Number of parameters, Gradient, Chi squared, parameters, what to do
  
  AliPHOSEmcRecPoint * emcRP = (AliPHOSEmcRecPoint *) gMinuit->GetObjectFit() ; // EmcRecPoint to fit

  Int_t * emcDigits     = emcRP->GetDigitsList() ;

  Int_t nOfDigits = emcRP->GetDigitsMultiplicity() ; 

  Float_t * emcEnergies = emcRP->GetEnergiesList() ;

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;

  AliPHOSIndexToObject * please = AliPHOSIndexToObject::GetInstance() ;

  fret = 0. ;     
  Int_t iparam ;

  if(iflag == 2)
    for(iparam = 0 ; iparam < nPar ; iparam++)    
      Grad[iparam] = 0 ; // Will evaluate gradient
  
  Double_t efit ;    

  AliPHOSDigit * digit ;
  Int_t iDigit ;

  for( iDigit = 0 ; iDigit < nOfDigits ; iDigit++) {

    digit = please->GimeDigit( emcDigits[iDigit] ) ; 

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
  }

}
