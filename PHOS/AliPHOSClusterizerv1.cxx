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
// A brief description of the class
//*-- Author : Yves Schutz  SUBATECH 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TMath.h" 

// --- Standard library ---

#include <iostream>

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
  // neigbours are defined as digits having at least common vertex
  // The order of A and B in AreNeighbours(A,B) is important: first (A) should be digit 
  // in cluster, which compared with digits not clusterized yet  
  Int_t rv = 0 ; 

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;

  Int_t relid1[4] ; 
  geom->AbsToRelNumbering(d1->GetId(), relid1) ; 

  Int_t relid2[4] ; 
  geom->AbsToRelNumbering(d2->GetId(), relid2) ; 
 
  if ( (relid1[0] == relid2[0]) && (relid1[1]==relid2[1]) ) { // inside the same PHOS module and the same PPSD Module 
    Int_t RowDiff = TMath::Abs( relid1[2] - relid2[2] ) ;  
    Int_t ColDiff = TMath::Abs( relid1[3] - relid2[3] ) ;  
    
    if (( ColDiff<=1 )  && ( RowDiff <= 1 )){
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
  // copies the digits with energy above thershold and sorts the list
  // according to increasing Id number

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;
  Int_t relid[4] ;  
  
  TIter next(dl) ; 
  AliPHOSDigit * digit ;
  
  while ( (digit = (AliPHOSDigit *)next()) ) { 
    
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
 numb[0] = fNumberOfEmcClusters ; 
 numb[1] = fNumberOfPpsdClusters ; 
}

//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::IsInEmc(AliPHOSDigit * digit) 
{
  Bool_t rv = kFALSE ; 

  AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ;  

  Int_t relid[4] ; 
  geom->AbsToRelNumbering(digit->GetId(), relid) ; 

  if ( relid[1] == 0  )
    rv = kTRUE; 

  return rv ; 
}

//____________________________________________________________________________
void AliPHOSClusterizerv1::MakeClusters(const DigitsList * dl, RecPointsList * emcl, RecPointsList * ppsdl)
{
  // Fill and sort the working digits list
  TObjArray TempoDigitsList( dl->GetEntries() ) ;
  FillandSort(dl, &TempoDigitsList) ; 

  // Clusterization starts  
  TIter nextdigit(&TempoDigitsList) ; 
  AliPHOSDigit * digit ; 
  Bool_t NotRemoved = kTRUE ;
  cout << "Cluster1" << endl;
  while ( (digit = (AliPHOSDigit *)nextdigit()) ) { // scan over the list of digits
    AliPHOSRecPoint * clu ; 
   
    int * ClusterDigitsList[dl->GetEntries()] ;   
    Int_t index ;
    if (( ( IsInEmc(digit) ) && ( Calibrate(digit->GetAmp() ) > fEmcClusteringThreshold ) ) || 
        ( ( !IsInEmc(digit) ) && ( Calibrate(digit->GetAmp() ) > fPpsdClusteringThreshold ) ) ) {
  
      Int_t iDigitInCluster = 0 ; 

      if  ( IsInEmc(digit) ) { 
        emcl->Print("");
	cout << "Cluster2 EMC" << endl;  
        new ((*emcl)[fNumberOfEmcClusters]) AliPHOSEmcRecPoint(fW0, fLocMaxCut) ; // start a new EMC RecPoint
	clu = (AliPHOSEmcRecPoint *) (*emcl)[fNumberOfEmcClusters] ; 
	fNumberOfEmcClusters++ ; 
	clu->AddDigit(*digit, Calibrate(digit->GetAmp())) ; 

	ClusterDigitsList[iDigitInCluster] = (int* ) digit ;	
	iDigitInCluster++ ; 
	TempoDigitsList.Remove(digit) ; 
      }

      else { 
	cout << "Cluster2 PPSD" << endl;
	new ((*ppsdl)[fNumberOfPpsdClusters]) AliPHOSPpsdRecPoint() ; // start a new PPSD cluster
	clu =  (AliPHOSPpsdRecPoint *) ppsdl->At(fNumberOfPpsdClusters)  ;  
	fNumberOfPpsdClusters++ ; 
	clu->AddDigit(*digit, 0.) ;	
	ClusterDigitsList[iDigitInCluster] = (int* ) digit ;	
	iDigitInCluster++ ; 
	TempoDigitsList.Remove(digit) ; 
        nextdigit.Reset() ;
	
	// Here we remove resting EMC digits, which cannot make cluster

        if( NotRemoved ) { 
	  
	  while( ( digit = (AliPHOSDigit *)nextdigit() ) ) {
	    
            if( IsInEmc(digit) ) 
	      TempoDigitsList.Remove(digit) ;
            else 
	      break ;
	  
	  } // while digit  
	  
	} // if NotRemoved 
	
      } // else        
      
      nextdigit.Reset() ;

      AliPHOSDigit * digitN ; 
      index = 0 ;
      while (index < iDigitInCluster){ // scan over digits already in cluster 
	digit = (AliPHOSDigit *) ClusterDigitsList[index]  ;      
	index++ ; 
        while ( (digitN = (AliPHOSDigit *)nextdigit()) ) { // scan over the reduced list of digits 
	  Int_t ineb = AreNeighbours(digit, digitN);   //  call (digit,digitN) in THAT oder !!!!!
          switch (ineb ) {
          case 0 :   // not a neighbour
	    break ;	 
	  case 1 :   // are neighbours 
	    clu->AddDigit( *digitN, Calibrate( digitN->GetAmp() ) ) ;
	    ClusterDigitsList[iDigitInCluster] =(int*) digitN ; 
	    iDigitInCluster++ ; 
	    TempoDigitsList.Remove(digitN) ;
	    break ;
          case 2 :   // too far from each other
	    goto endofloop;   
	  } // switch
	  
	} // while digitN

      endofloop: ;
	nextdigit.Reset() ; 
	
      } // loop over cluster     
 
   }  //below energy theshold  
  
  } // while digit

  TempoDigitsList.Clear() ; 
}

//____________________________________________________________________________
void AliPHOSClusterizerv1::PrintParameters() 
{
  cout << "PHOS Clusterizer version 1 :" << endl 
       << "                       EMC  Clustering threshold = " << fEmcClusteringThreshold << endl
       << "                       EMC  Energy threshold     = " << fEmcEnergyThreshold << endl                  
       << "                      PPSD  Clustering threshold = " << fPpsdClusteringThreshold << endl
       << "                      PPSD  Energy threshold     = " << fPpsdEnergyThreshold << endl ;                
}
