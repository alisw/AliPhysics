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

/*
$Log$
Revision 1.2  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

Revision 1.1.2.1  2000/06/09 21:48:28  morsch
Code from AliMUONSegResTrigger.cxx

*/

#include "AliMUONResponseTrigger.h"
#include "AliSegmentation.h"
#include <TMath.h>
#include <TRandom.h>
#include <iostream.h> 


ClassImp(AliMUONResponseTrigger)

//------------------------------------------------------------------   
Float_t AliMUONResponseTrigger::IntXY(AliSegmentation * segmentation){    
// Returns 1 or 0 if the current strip is fired or not according 
// to the cluster size and the width of the main strip.
// For the time being the probability to fire a neighbour depends
// only on the width of the main strip and is limited to a maximum 
// cluster-size of 2. 
// The corresponding probabilities are given below (O.Roig PhD Thesis)  
// This will be improved in the future by including a parametrization
// of the cluster size as a function of the position of the physical 
// hit with respect to the center of the strip.
//------------------------------------------------------------------
//  clust. size =      1      2      3     4      5     >5
//  strip width = 1 | 54.7 | 44.5 | 0.7 | 0.06 | 0.04 | 0.0 |
//  strip width = 2 | 89.0 | 10.7 | 0.2 | 0.1  | 0.0  | 0.0 |
//  strip width = 4 | 99.0 |  1.0 | 0.0 | 0.0  | 0.0  | 0.0 |
//------------------------------------------------------------------

//  cout << "in AliMUONResponseTrigger::IntXY" << "\n";    

  // get the "parameters" needed to evaluate the strip response
  // x1    : hit x(y) position
  // x2    : x(y) coordinate of the main strip
  // x3    : current strip real x(y) coordinate  
  // width : width of the main strip 
    Float_t x1,x2,x3,width;
    segmentation->IntegrationLimits(x1,x2,x3,width);  
    //  cout << " x or y main & current = " << x2 << " , " << x3 
    //    << " width main = " << width << "\n";

    /*
    if (TMath::Abs(x3-x1)<TMath::Abs(x3-x2)) { // find neighbour candidate
	Int_t iwidth=Int_t(width);
	Float_t rand = gRandom->Rndm()*100.; 
	if (iwidth==1) {
	    if (rand<44.5) { return 1.; } 
	    else           { return 0.; } 
	}  else if (iwidth==2) {
	    if (rand<10.7) { return 1.; } 
	    else           { return 0.; } 
	} else if (iwidth==4)  {
	    if (rand<1.)   { return 1.; } 
	    else           { return 0.; }
	}
    } else { return 0.;} 
    return -1;    
    */
    return 0;
}


//------------------------------------------------------------------   
Int_t  AliMUONResponseTrigger::DigitResponse(Int_t digit)
{
//
//  only digital (0/1) information available
  if (digit) digit=1;
  return digit;
}






