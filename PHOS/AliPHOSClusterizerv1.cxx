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

/* $Log:
   1 October 2000. Yuri Kharlov:
     AreNeighbours()
     PPSD upper layer is considered if number of layers>1

   18 October 2000. Yuri Kharlov:
     AliPHOSClusterizerv1()
     CPV clusterizing parameters added

     MakeClusters()
     After first PPSD digit remove EMC digits only once
*/

//_________________________________________________________________________
//  Implementation version 1 of the clusterization algorithm 
// 
//*-- Author: Yves Schutz (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

#include <assert.h>

// --- ROOT system ---

#include "TMath.h" 

// --- Standard library ---

#include <iostream.h>

// --- AliRoot header files ---

#include "AliPHOSClusterizerv1.h"
#include "AliPHOSDigit.h"
#include "AliPHOSEmcRecPoint.h"
#include "AliPHOSPpsdRecPoint.h"
#include "AliPHOSCpvRecPoint.h"
#include "AliPHOSv0.h" 
#include "AliRun.h" 

ClassImp(AliPHOSClusterizerv1)

//____________________________________________________________________________
AliPHOSClusterizerv1::AliPHOSClusterizerv1()
{
  // default ctor (to be used)

  fA                       = 0.;
  fB                       = 0.0000001 ;
  fGeom                    = AliPHOSGeometry::GetInstance();
  fNumberOfEmcClusters     = 0 ; 
  fNumberOfPpsdClusters    = 0 ; 
  fEmcClusteringThreshold  = 0.2;   
  fEmcEnergyThreshold      = 0.05;    
  fPpsdClusteringThreshold = 0.0000002 ;
  fPpsdEnergyThreshold     = 0.0000002 ;  
  fCpvClusteringThreshold  = 0.0;
  fCpvEnergyThreshold      = 0.09;  
  fW0                      = 4.5 ;
  fLocMaxCut               = 0.03 ;
  fW0CPV                   = 4.0 ;
  fLocMaxCutCPV            = 0.03 ;

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

  Int_t relid1[4] ; 
  fGeom->AbsToRelNumbering(d1->GetId(), relid1) ; 

  Int_t relid2[4] ; 
  fGeom->AbsToRelNumbering(d2->GetId(), relid2) ; 
 
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

  //Do NOT clusterize upper PPSD  
  if( IsInPpsd(d1) && IsInPpsd(d2) &&
     relid1[1] > 0                 &&
     relid1[1] < fGeom->GetNumberOfPadsPhi()*fGeom->GetNumberOfPadsPhi() ) rv = 2 ;

  return rv ; 
}

//____________________________________________________________________________
void AliPHOSClusterizerv1::FillandSort(const DigitsList * dl, TObjArray * tl) 
{
  // Copies the digits with energy above thershold and sorts the list
  // according to increasing Id number

  TIter next(dl) ; 
  AliPHOSDigit * digit ;
  
  while ( (digit = (AliPHOSDigit *)next()) ) { 
    Float_t ene = Calibrate(digit->GetAmp()) ; 
    if      ( IsInEmc  (digit) && ene > fEmcEnergyThreshold  )
      tl->Add(digit) ;
    else if ( IsInPpsd (digit) && ene > fPpsdEnergyThreshold )
      tl->Add(digit) ;
    else if ( IsInCpv  (digit) && ene > fCpvEnergyThreshold  )
      tl->Add(digit) ;
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

  Int_t relid[4] ; 
  fGeom->AbsToRelNumbering(digit->GetId(), relid) ; 

  if ( relid[1] == 0  ) rv = kTRUE; 

  return rv ; 
}

//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::IsInPpsd(AliPHOSDigit * digit) 
{
  // Tells if (true) or not (false) the digit is in a PHOS-PPSD module
 
  Bool_t rv = kFALSE ; 

  Int_t relid[4] ; 
  fGeom->AbsToRelNumbering(digit->GetId(), relid) ; 

  if ( relid[1] > 0 && relid[0] > fGeom->GetNCPVModules() ) rv = kTRUE; 

  return rv ; 
}

//____________________________________________________________________________
Bool_t AliPHOSClusterizerv1::IsInCpv(AliPHOSDigit * digit) 
{
  // Tells if (true) or not (false) the digit is in a PHOS-CPV module
 
  Bool_t rv = kFALSE ; 

  Int_t relid[4] ; 
  fGeom->AbsToRelNumbering(digit->GetId(), relid) ; 

  if ( relid[1] > 0 && relid[0] <= fGeom->GetNCPVModules() ) rv = kTRUE; 

  return rv ; 
}

//____________________________________________________________________________
void AliPHOSClusterizerv1::MakeClusters(const DigitsList * dl, 
					AliPHOSRecPoint::RecPointsList * emcl, 
					AliPHOSRecPoint::RecPointsList * ppsdl)
{
  // Steering method to construct the clusters stored in a list of Reconstructed Points
  // A cluster is defined as a list of neighbour digits

  fNumberOfEmcClusters  = 0 ;
  fNumberOfPpsdClusters = 0 ;
  fNumberOfCpvClusters  = 0 ;

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

    if (( IsInEmc (digit) && Calibrate(digit->GetAmp()) > fEmcClusteringThreshold  ) || 
        ( IsInPpsd(digit) && Calibrate(digit->GetAmp()) > fPpsdClusteringThreshold ) ||
        ( IsInCpv (digit) && Calibrate(digit->GetAmp()) > fCpvClusteringThreshold  ) ) {
      
      Int_t iDigitInCluster = 0 ; 

      if  ( IsInEmc(digit) ) {   
	// start a new EMC RecPoint
	if(fNumberOfEmcClusters >= emcl->GetSize()) emcl->Expand(2*fNumberOfEmcClusters+1) ;
	(*emcl)[fNumberOfEmcClusters] = new  AliPHOSEmcRecPoint(fW0, fLocMaxCut) ;
	clu = (AliPHOSEmcRecPoint *) emcl->At(fNumberOfEmcClusters) ; 
	fNumberOfEmcClusters++ ; 
	clu->AddDigit(*digit, Calibrate(digit->GetAmp())) ; 
	clusterdigitslist[iDigitInCluster] = digit ;	
	iDigitInCluster++ ; 
	tempodigitslist.Remove(digit) ; 

      } else { 
	
	// start a new PPSD/CPV cluster
	if(fNumberOfPpsdClusters >= ppsdl->GetSize()) ppsdl->Expand(2*fNumberOfPpsdClusters+1);
	if      (IsInPpsd(digit)) 
	  (*ppsdl)[fNumberOfPpsdClusters] = new AliPHOSPpsdRecPoint() ;
	else
	  (*ppsdl)[fNumberOfPpsdClusters] = new AliPHOSCpvRecPoint(fW0CPV, fLocMaxCutCPV) ;
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
	  }
	  notremoved = kFALSE ;
	}
	
      } // else        
      
      nextdigit.Reset() ;
      
      AliPHOSDigit * digitN ; 
      index = 0 ;
      while (index < iDigitInCluster){ // scan over digits already in cluster 
	digit =  clusterdigitslist[index]  ;      
	index++ ; 
        while ( (digitN = (AliPHOSDigit *)nextdigit()) ) { // scan over the reduced list of digits 
	  Int_t ineb = AreNeighbours(digit, digitN);       // call (digit,digitN) in THAT oder !!!!!
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

  ppsdl->Sort() ;
  Int_t index ;
  for(index = 0; index < ppsdl->GetEntries(); index++)
    ((AliPHOSPpsdRecPoint *)ppsdl->At(index))->SetIndexInList(index) ;
}

//____________________________________________________________________________
void AliPHOSClusterizerv1::PrintParameters() 
{
  // Print the energy thresholds 

  cout << "PHOS Clusterizer version 1 :" << endl 
       << "                       EMC  Clustering threshold = " << fEmcClusteringThreshold << endl
       << "                       EMC  Energy threshold     = " << fEmcEnergyThreshold << endl                  
       << "                      PPSD  Clustering threshold = " << fPpsdClusteringThreshold << endl
       << "                      PPSD  Energy threshold     = " << fPpsdEnergyThreshold << endl
       << "                       CPV  Clustering threshold = " << fCpvClusteringThreshold << endl
       << "                       CPV  Energy threshold     = " << fCpvEnergyThreshold << endl ;                
}
