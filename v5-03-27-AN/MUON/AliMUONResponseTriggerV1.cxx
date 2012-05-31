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

//-----------------------------------------------------------------------------
// Class AliMUONTriggerResponseV1
// ------------------
// Trigger chamber response 
// with cluster size activated
//-----------------------------------------------------------------------------

#include "AliMUONResponseTriggerV1.h"
#include "AliMUON.h"
#include "AliMUONDigit.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONConstants.h"

#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpCathodType.h"

#include "AliRun.h"

#include <TMath.h>
#include <TRandom.h>

/// \cond CLASSIMP
ClassImp(AliMUONResponseTriggerV1)
/// \endcond

namespace
{
  AliMUON* muon()
  {
    return static_cast<AliMUON*>(gAlice->GetModule("MUON"));
  }

  void Global2Local(Int_t detElemId, Double_t xg, Double_t yg, Double_t zg,
                  Double_t& xl, Double_t& yl, Double_t& zl)
  {  
  // ideally should be : 
  // Double_t x,y,z;
  // AliMUONGeometry::Global2Local(detElemId,xg,yg,zg,x,y,z);
  // but while waiting for this geometry singleton, let's go through
  // AliMUON still.
  
    const AliMUONGeometryTransformer* transformer = muon()->GetGeometryTransformer();
    transformer->Global2Local(detElemId,xg,yg,zg,xl,yl,zl);
  }

}

//------------------------------------------------------------------   
AliMUONResponseTriggerV1::AliMUONResponseTriggerV1()
    : AliMUONResponseTrigger(),
      fGenerCluster(0),
      fA(0),
      fB(0),       
      fC(0)
{
/// default constructor 
  Float_t hv=9.2;
  SetParameters(hv);
}

//------------------------------------------------------------------   
AliMUONResponseTriggerV1::AliMUONResponseTriggerV1(Float_t hv)
    : AliMUONResponseTrigger(),
      fGenerCluster(0),
      fA(0),
      fB(0),       
      fC(0)
{
/// Constructor 
  SetParameters(hv);
}

//------------------------------------------------------------------   
AliMUONResponseTriggerV1::~AliMUONResponseTriggerV1()
{
/// destructor 
}

//------------------------------------------------------------------   
void AliMUONResponseTriggerV1::SetParameters(Float_t hv)
{
/// initialize parameters accoring to HV
/// (see V.Barret B.Espagnon and P.Rosnet Alice/note xxx)
/// this parametrisation is valid only for the "streamer" mode
  fA = 6.089 * hv - 52.70;
  fB = 2.966;
  fC = 4.3e-4 * hv - 3.5e-3;
}

//------------------------------------------------------------------   
Int_t AliMUONResponseTriggerV1::SetGenerCluster()
{
/// Set the GenerCluster parameter and return 1
  fGenerCluster = gRandom->Rndm();
  return 1;
}

//------------------------------------------------------------------   
Float_t AliMUONResponseTriggerV1::FireStripProb(Float_t x4, Float_t theta)
const
{
/// parametrisation of the probability that a strip neighbour of the main 
/// strip is fired (V.Barret B.Espagnon and P.Rosnet INT/DIM/01-04 (2001)
/// WARNING : need to convert x4 from cm to mm
/// this parametrisation is valid only for the "streamer" mode

 return 
     (TMath::Cos(theta)*fA/(fA+TMath::Cos(theta)*TMath::Power(x4*10.,fB))+fC)/
     (TMath::Cos(theta)+fC);
}

//------------------------------------------------------------------  
void AliMUONResponseTriggerV1::DisIntegrate(const AliMUONHit& hit, TList& digits, Float_t /*timeDif*/)
{
  /// Generate digits (on each cathode) from 1 hit, with cluster-size
  /// generation.
  
  digits.Clear();
  
  Float_t xhit = hit.X();
  Float_t yhit = hit.Y();
  Float_t zhit = hit.Z();
  Int_t detElemId = hit.DetElemId();  
  
  Double_t x,y,z;
  Global2Local(detElemId,xhit,yhit,zhit,x,y,z);
  
  Float_t tof = hit.Age();
  Int_t twentyNano(100);
  if (tof<AliMUONConstants::TriggerTofLimit())
  {
    twentyNano=1;
  }
  
  Int_t nboard = 0;

  for ( Int_t cath = AliMp::kCath0; cath <= AliMp::kCath1; ++cath )
  {
    const AliMpVSegmentation* seg 
      = AliMpSegmentation::Instance()
        ->GetMpSegmentation(detElemId,AliMp::GetCathodType(cath));

    AliMpPad pad = seg->PadByPosition(x,y,kFALSE);
    Int_t ix = pad.GetIx();
    Int_t iy = pad.GetIy();
    
    AliDebug(1,Form("xhit,yhit=%e,%e lx,ly,lz=%e,%e,%e ix,iy=%d,%d",
                    xhit,yhit,x,y,z,ix,iy));
    
    if ( !pad.IsValid() )
    {
      AliWarning(Form("hit w/o strip %d-%d xhit,yhit=%e,%e local x,y,z "
                      "%e,%e,%e ix,iy=%d,%d",detElemId,
                      cath,
                      xhit,yhit,x,y,z,ix,iy));
      continue;
    }
    
    if ( cath == AliMp::kCath0 ) nboard = pad.GetLocalBoardId(0);
        
    AliMUONDigit* d = new AliMUONDigit(detElemId,nboard,
                                       pad.GetLocalBoardChannel(0),
                                       cath);
    d->SetPadXY(ix,iy);

    d->SetCharge(twentyNano);


    digits.Add(d);

    SetGenerCluster(); // 1 randum number per cathode (to be checked)

    Int_t xList[10], yList[10];
    Neighbours(cath,ix,iy,xList,yList);
    
//    cout << " detElemId cath ix iy = " << detElemId << " " << cath 
//	 << " " << ix << " " << iy << "\n";
//    for (Int_t i=0; i<10; i++) cout << " " << xList[i] << " " << yList[i];
//    cout << "\n";

    Int_t qp = 0; // fired/no-fired strip = 1/0
    for (Int_t i=0; i<10; i++) { // loop on neighbors
	if (i==0||i==5||qp!=0) { // built-up cluster
	    
	    // need to iterate in iy/ix for bending/non-bending plane
	    Int_t ixNeigh = ( cath == 0 ) ? ix : xList[i];
	    Int_t iyNeigh = ( cath == 0 ) ? yList[i] : iy;
	    
	    AliMpPad padNeigh = seg->PadByIndices(ixNeigh,iyNeigh,kFALSE);
	    if(padNeigh.IsValid()){ // existing neighbourg		
		
		Int_t dix=-(ixNeigh-ix);
		Int_t diy=-(iyNeigh-iy);
		Float_t xlocalNeigh = padNeigh.GetPositionX();
		Float_t ylocalNeigh = padNeigh.GetPositionY();
		Float_t dpx = padNeigh.GetDimensionX();
		Float_t dpy = padNeigh.GetDimensionY();
		Float_t distX = TMath::Abs((Float_t)dix) * ((Float_t)dix * dpx + xlocalNeigh - x);
		Float_t distY = TMath::Abs((Float_t)diy) * ((Float_t)diy * dpy + ylocalNeigh - y);
		Float_t dist = TMath::Sqrt(distX*distX+distY*distY);
		
//		cout << " here " << dist << " " << fGenerCluster << " " << FireStripProb(dist,0) << "\n";
		
		if (fGenerCluster<FireStripProb(dist,0)) qp = 1;
		else qp = 0;
		
		if (qp == 1) { // this digit is fired    
        Int_t neighBoard = 0;
        if ( cath == AliMp::kCath0 ) neighBoard = padNeigh.GetLocalBoardId(0);
        else {
          const AliMpVSegmentation* seg0 
            = AliMpSegmentation::Instance()
              ->GetMpSegmentation(detElemId,AliMp::GetCathodType(AliMp::kCath0));
          AliMpPad padNeigh0 = seg0->PadByPosition(xlocalNeigh, y, kFALSE);
          if ( ! padNeigh0.IsValid() ) continue; // This can happen only on the cut RPC, at boards 25, 30, 142 and 147
          neighBoard = padNeigh0.GetLocalBoardId(0);
        }
		    AliMUONDigit* dNeigh = new AliMUONDigit(detElemId,neighBoard,
                                                padNeigh.GetLocalBoardChannel(0),
                                                cath);
		    
		    dNeigh->SetPadXY(ixNeigh,iyNeigh);      
		    dNeigh->SetCharge(twentyNano);
		    digits.Add(dNeigh);
		} // digit fired		
	    } // pad is valid
	} // built-up cluster
    } // loop on neighbors
  } // loop on cathode
}

//------------------------------------------------------------------  
void AliMUONResponseTriggerV1::Neighbours(const Int_t cath, 
					  const Int_t ix, const Int_t iy, 
					  Int_t Xlist[10], Int_t Ylist[10]) const
{
    ///-----------------BENDING-----------------------------------------      /n
    /// Returns list of 10 next neighbours for given X strip (ix, iy)         /n
    /// neighbour number 4 in the list -                                      /n    
    /// neighbour number 3 in the list  |                                     /n   
    /// neighbour number 2 in the list  |_ Upper part                         /n         
    /// neighbour number 1 in the list  |                                     /n    
    /// neighbour number 0 in the list -                                      /n   
    ///      X strip (ix, iy)                                                 /n
    /// neighbour number 5 in the list -                                      /n
    /// neighbour number 6 in the list  | _ Lower part                        /n
    /// neighbour number 7 in the list  |                                     /n
    /// neighbour number 8 in the list  |                                     /n
    /// neighbour number 9 in the list -                                      /n
    ///                                                                       /n
    ///-----------------NON-BENDING-------------------------------------      /n
    /// Returns list of 10 next neighbours for given Y strip (ix, iy)         /n 
    /// neighbour number 9 8 7 6 5 (Y strip (ix, iy)) 0 1 2 3 4 in the list   /n 
    ///                  |_______|                    |_______/               /n 

    ///                    left                         right                 /n
    
    for (Int_t i=0; i<10; i++) {
	Xlist[i]=-1;
	Ylist[i]=-1;
    }
    
    Int_t iList[10]={9,8,7,6,5, 0,1,2,3,4};

    // need to iterate in iy/ix for bending/non-bending plane
    Int_t iNeigh = ( cath == 0 ) ? iy : ix;
    
    Int_t i=0;
    for (Int_t j=iNeigh-5; j<=iNeigh+5; j++){
	if (j == iNeigh) continue;
	// need to iterate in iy/ix for bending/non-bending plane
	Int_t ixNeigh = ( cath == 0 ) ? ix : j;
	Int_t iyNeigh = ( cath == 0 ) ? j : iy;
	
//	cout << " " << cath << " " << ix << " " << iy 
//	     << " "  << ixNeigh << " " << iyNeigh << "\n";
	
	Xlist[iList[i]]=ixNeigh;	
	Ylist[iList[i]]=iyNeigh;	
	i++;
    } 
}

