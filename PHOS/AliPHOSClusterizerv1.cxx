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
//  Implementation version 1 of the clusterization algorithm 
// 
//*-- Author: Yves Schutz (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TMath.h" 

// --- Standard library ---

#include <iostream.h>

// --- AliRoot header files ---

#include "AliPHOSClusterizerv1.h"
#include "AliPHOSDigit.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSPpsdRecPoint.h"
#include "AliPHOSv0.h" 
#include "AliRun.h" 

ClassImp(AliPHOSClusterizerv1)

//____________________________________________________________________________
AliPHOSClusterizerv1::AliPHOSClusterizerv1()
{
  // default ctor (to be used)

  fA                       = 0.;
  fB                       = 0.01 ;
  fNumberOfEmcClusters     = 0 ; 
  fNumberOfPpsdClusters    = 0 ; 
  fEmcClusteringThreshold  = 0.1;   
  fEmcEnergyThreshold      = 0.01;    
  fPpsdClusteringThreshold = 0.00000015; 
  fPpsdEnergyThreshold     = 0.0000001;  
  fW0                      = 4.5 ;
  fLocMaxCut               = 0.06 ;
}

//____________________________________________________________________________
Int_t AliPHOSClusterizerv1::AreNeighbours(AliPHOSDigit * d1, AliPHOSDigit * d2)
{
  // Gives the neighbourness of two digits = 0 are not neighbour but continue searching 
  //                                       = 1 are neighbour
  //                                       = 2 are not neighbour but do not continue searching
  // neighbours are defined as digits having at least common vertex
  // The order of d1 and d2 is important: first (d1) should be a digit already in a cluster 
  //                                      which is compared to a digit (d2)  not yet in a cluster  

  Int_t rv = 0 ; 

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;

  Int_t relid1[4] ; 
  geom->AbsToRelNumbering(d1->GetId(), relid1) ; 

  Int_t relid2[4] ; 
  geom->AbsToRelNumbering(d2->GetId(), relid2) ; 
 
  if ( (relid1[0] == relid2[0]) && (relid1[1]==relid2[1]) ) { // inside the same PHOS module and the same PPSD Module 
    Int_t rowdiff = TMath::Abs( relid1[2] - relid2[2] ) ;  
    Int_t coldiff = TMath::Abs( relid1[3] - relid2[3] ) ;  
    
    if (( coldiff <= 1 )  && ( rowdiff <= 1 )){
      rv = 1 ; 
    }
    else {
      if((relid2[2] > relid1[2]) && (relid2[3] > relid1[3]+1)) 
	rv = 2; //  Difference in row numbers is too large to look further 
    }

  } 
  else {
    
    if( (relid1[0] < relid2[0]) || (relid1[1] < relid2[1]) )  
      rv=2 ;

  }
  
  return rv ; 
}

//____________________________________________________________________________
void AliPHOSClusterizerv1::FillandSort(const DigitsList * dl, TObjArray * tl) 
{
  // Copies the digits with energy above thershold and sorts the list
  // according to increasing Id number

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;
  Int_t relid[4] ;  
  
  TIter next(dl) ; 
  AliPHOSDigit * digit ;
  
 
 

  while ( (digit = (AliPHOSDigit *)next()) ) { 

//     cout << " clusterizerv1 " << endl ;
//     int nprim = digit->GetNprimary() ;
//     int * aprim = digit->GetPrimary() ;
//     for ( int ii = 0 ; ii < nprim ; ii++)
//       cout << ii << " prim = " << aprim[ii] << endl ;

    Int_t id    = digit->GetId() ; 
    Float_t ene = Calibrate(digit->GetAmp()) ; 
    geom->AbsToRelNumbering(id, relid) ;
    if(relid[1]==0){ // EMC
      if ( ene > fEmcEnergyThreshold )
	tl->Add(digit) ;
    }

    else { //Ppsd
      if ( ene > fPpsdEnergyThreshold )
	tl->Add(digit) ; 
    }

  }
  tl->Sort() ; 
}

//____________________________________________________________________________
void AliPHOSClusterizerv1:: GetNumberOfClustersFound(Int_t * numb) 
{
  // Fills numb with the number of EMC  (numb[0]) clusters found
  //                               PPSD (numb[1]) clusters found

  numb[0] = fNumberOfEmcClusters ; 
  numb[1] = fNumberOfPpsdClusters ; 
}

//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::IsInEmc(AliPHOSDigit * digit) 
{
  // Tells if (true) or not (false) the digit is in a PHOS-EMC module
 
  Bool_t rv = kFALSE ; 

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;  

  Int_t relid[4] ; 
  geom->AbsToRelNumbering(digit->GetId(), relid) ; 

  if ( relid[1] == 0  )
    rv = kTRUE; 

  return rv ; 
}

//____________________________________________________________________________
void AliPHOSClusterizerv1::MakeClusters(const DigitsList * dl, 
					AliPHOSRecPoint::RecPointsList * emcl, 
					AliPHOSRecPoint::RecPointsList * ppsdl)
{
  // Steering method to construct the clusters stored in a list of Reconstructed Points
  // A cluster is defined as a list of neighbour digits
  
  // Fill and sort the working digits list
  TObjArray tempodigitslist( dl->GetEntries() ) ;
  FillandSort(dl, &tempodigitslist) ; 

  // Clusterization starts  
  TIter nextdigit(&tempodigitslist) ; 
  AliPHOSDigit * digit ; 
  Bool_t notremoved = kTRUE ;

  while ( (digit = (AliPHOSDigit *)nextdigit()) ) { // scan over the list of digits
    AliPHOSRecPoint * clu ; 
   
    AliPHOSDigit ** clusterdigitslist = new AliPHOSDigit*[dl->GetEntries()] ;   
    Int_t index ;
    if (( ( IsInEmc(digit) ) && ( Calibrate(digit->GetAmp() ) > fEmcClusteringThreshold ) ) || 
        ( ( !IsInEmc(digit) ) && ( Calibrate(digit->GetAmp() ) > fPpsdClusteringThreshold ) ) ) {
  
      Int_t iDigitInCluster = 0 ; 

      if  ( IsInEmc(digit) ) {   
	// start a new EMC RecPoint
	//        new ((*emcl)[fNumberOfEmcClusters]) AliPHOSEmcRecPoint(fW0, fLocMaxCut) ; if TClonesArray
	(*emcl)[fNumberOfEmcClusters] = new  AliPHOSEmcRecPoint(fW0, fLocMaxCut) ;

	clu = (AliPHOSEmcRecPoint *) (*emcl)[fNumberOfEmcClusters] ; 
	fNumberOfEmcClusters++ ; 
	clu->AddDigit(*digit, Calibrate(digit->GetAmp())) ; 

	clusterdigitslist[iDigitInCluster] = digit ;	
	iDigitInCluster++ ; 
	tempodigitslist.Remove(digit) ; 


      }

      else { 
	
	// start a new PPSD cluster
	// new ((*ppsdl)[fNumberOfPpsdClusters]) AliPHOSPpsdRecPoint() ;  if TClonesArray
	(*ppsdl)[fNumberOfPpsdClusters] = new AliPHOSPpsdRecPoint() ;
	
	clu =  (AliPHOSPpsdRecPoint *) ppsdl->At(fNumberOfPpsdClusters)  ;  
	fNumberOfPpsdClusters++ ; 
	clu->AddDigit(*digit, Calibrate(digit->GetAmp()) ) ;	
	clusterdigitslist[iDigitInCluster] = digit  ;	
	iDigitInCluster++ ; 
	tempodigitslist.Remove(digit) ; 
        nextdigit.Reset() ;
	
	// Here we remove resting EMC digits, which cannot make cluster

        if( notremoved ) { 
	  
	  while( ( digit = (AliPHOSDigit *)nextdigit() ) ) {
	    
            if( IsInEmc(digit) ) 
	      tempodigitslist.Remove(digit) ;
            else 
	      break ;
	  
	  } // while digit  
	  
	} // if notremoved 
	
      } // else        
      
      nextdigit.Reset() ;

      AliPHOSDigit * digitN ; 
      index = 0 ;
      while (index < iDigitInCluster){ // scan over digits already in cluster 
	digit =  clusterdigitslist[index]  ;      
	index++ ; 
        while ( (digitN = (AliPHOSDigit *)nextdigit()) ) { // scan over the reduced list of digits 
	  Int_t ineb = AreNeighbours(digit, digitN);   //  call (digit,digitN) in THAT oder !!!!!
          switch (ineb ) {
          case 0 :   // not a neighbour
	    break ;	 
	  case 1 :   // are neighbours 
	    clu->AddDigit(*digitN, Calibrate( digitN->GetAmp() ) ) ;
	    clusterdigitslist[iDigitInCluster] = digitN ; 
	    iDigitInCluster++ ; 
	    tempodigitslist.Remove(digitN) ;
	    break ;
          case 2 :   // too far from each other
	    goto endofloop;   
	  } // switch
	  
	} // while digitN

      endofloop: ;
	nextdigit.Reset() ; 
	
      } // loop over cluster     
 
   }  //below energy theshold  
  
    delete[] clusterdigitslist ; 
    
  } // while digit

  tempodigitslist.Clear() ; 

}

//____________________________________________________________________________
void AliPHOSClusterizerv1::PrintParameters() 
{
  // Print the energy thresholds 

  cout << "PHOS Clusterizer version 1 :" << endl 
       << "                       EMC  Clustering threshold = " << fEmcClusteringThreshold << endl
       << "                       EMC  Energy threshold     = " << fEmcEnergyThreshold << endl                  
       << "                      PPSD  Clustering threshold = " << fPpsdClusteringThreshold << endl
       << "                      PPSD  Energy threshold     = " << fPpsdEnergyThreshold << endl ;                
}
