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

#include "AliCDBManager.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"

#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpCathodType.h"

#include "AliRun.h"
#include "AliDCSValue.h"

#include <TMath.h>
#include <TRandom.h>
#include <TMap.h>

/// \cond CLASSIMP
ClassImp(AliMUONResponseTriggerV1)
/// \endcond

namespace
{
  AliMUON* muon()
  {
    return static_cast<AliMUON*>(gAlice->GetModule("MUON"));
  }

  void Global2Local(Int_t detElemId, Double_t xg, Double_t yg, Double_t zg, Double_t& xl, Double_t& yl, Double_t& zl)
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
AliMUONResponseTriggerV1::AliMUONResponseTriggerV1() : AliMUONResponseTrigger(), fGenerCluster(0), fHVvalues(), fBValues(), fWorkCondition(2)
{
}

//------------------------------------------------------------------   
AliMUONResponseTriggerV1::AliMUONResponseTriggerV1(Int_t mode) : AliMUONResponseTrigger(), fGenerCluster(0), fHVvalues(), fBValues(), fWorkCondition(mode)
{
  //mode: 1=streamer - 2=avalanche
  SetHV();
  SetBValues();
}

//------------------------------------------------------------------   
AliMUONResponseTriggerV1::~AliMUONResponseTriggerV1()
{
}

//------------------------------------------------------------------   
Int_t AliMUONResponseTriggerV1::SetGenerCluster()
{
  // Set the GenerCluster parameter and return 1
  fGenerCluster = gRandom->Rndm();
  return 1;
}

//------------------------------------------------------------------   
Float_t AliMUONResponseTriggerV1::FireStripProb(Float_t x4,Float_t theta,Int_t rpc,Int_t plane,Int_t cath) const
{
/// parametrisation of the probability that a strip neighbour of the main 
/// strip is fired
/// WARNING : need to convert x4 from cm to mm

  Float_t hv = fHVvalues.At(18*plane+rpc);
  Float_t parA, parB, parC;
  
  if(fWorkCondition == 2) //avalanche
    parB = fBValues.At(72*cath+18*plane+rpc);
  else //streamer
    parB = 2.966;
  
  
  parA = 6.089 * hv - 52.70;
  parC = 8.3e-4 * hv - 0.5e-3;  

 return (TMath::Cos(theta)*parA/(parA+TMath::Cos(theta)*TMath::Power(x4*10.,parB))+parC)/(TMath::Cos(theta)+parC);
}

//------------------------------------------------------------------
void AliMUONResponseTriggerV1::DisIntegrate(const AliMUONHit& hit, TList& digits, Float_t /*timeDif*/)
{
  /// Generate digits (on each cathode) from 1 hit, with cluster-size generation.
  
  digits.Clear();
  
  Float_t xhit = hit.X();
  Float_t yhit = hit.Y();
  Float_t zhit = hit.Z();
  Int_t detElemId = hit.DetElemId();
  Int_t plane = detElemId/100 - 11; //plane from 0 to 3
  Int_t rpc = detElemId%100; //rpc from 0 to 3
  
  Double_t x,y,z;
  Global2Local(detElemId,xhit,yhit,zhit,x,y,z);
  
  Float_t tof = hit.Age();
  Int_t twentyNano(100);
  if (tof<AliMUONConstants::TriggerTofLimit())
  {
    twentyNano=1;
  }
  
  
  Int_t nboard = 0;

  for(Int_t cath = AliMp::kCath0; cath <= AliMp::kCath1; ++cath)
  {
    const AliMpVSegmentation* seg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::GetCathodType(cath));
    
    AliMpPad pad = seg->PadByPosition(x,y,kFALSE);
    Int_t ix = pad.GetIx();
    Int_t iy = pad.GetIy();
    
    AliDebug(1,Form("xhit,yhit=%e,%e lx,ly,lz=%e,%e,%e ix,iy=%d,%d",xhit,yhit,x,y,z,ix,iy));
    
    if ( !pad.IsValid() )
    {
      AliWarning(Form("hit w/o strip %d-%d xhit,yhit=%e,%e local x,y,z ""%e,%e,%e ix,iy=%d,%d",detElemId,cath,xhit,yhit,x,y,z,ix,iy));
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
    
    Int_t xList[30], yList[30];
    Neighbours(cath,ix,iy,xList,yList);
    
    
    Int_t qp = 0; // fired/no-fired strip = 1/0
    for (Int_t i=0; i<30; i++)  // loop on neighbors
    {
      if (i==0 || i==15 || qp!=0) // built-up cluster // need to iterate in iy/ix for bending/non-bending plane
      {
	Int_t ixNeigh = (cath == 0) ? ix : xList[i];
	Int_t iyNeigh = (cath == 0) ? yList[i] : iy;
	
	AliMpPad padNeigh = seg->PadByIndices(ixNeigh,iyNeigh,kFALSE);
	if(padNeigh.IsValid()) // existing neighbourg
	{
	  Int_t dix=-(ixNeigh-ix);
	  Int_t diy=-(iyNeigh-iy);
	  Float_t xlocalNeigh = padNeigh.GetPositionX();
	  Float_t ylocalNeigh = padNeigh.GetPositionY();
	  Float_t dpx = padNeigh.GetDimensionX();
	  Float_t dpy = padNeigh.GetDimensionY();
	  Float_t distX = TMath::Abs((Float_t)dix) * ((Float_t)dix * dpx + xlocalNeigh - x);
	  Float_t distY = TMath::Abs((Float_t)diy) * ((Float_t)diy * dpy + ylocalNeigh - y);
	  Float_t dist = TMath::Sqrt(distX*distX+distY*distY);
		
	  if(fGenerCluster < FireStripProb(dist,0,rpc,plane,cath))
	    qp = 1;
	  else
	    qp = 0;
		
	  if(qp == 1)
	  { // this digit is fired    
	    AliMUONDigit* dNeigh = new AliMUONDigit(detElemId,padNeigh.GetLocalBoardId(0),padNeigh.GetLocalBoardChannel(0),cath);
	    
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
void AliMUONResponseTriggerV1::SetHV()
{
  //
  /// Set HV values from OCDB
  //
  fHVvalues.Set(72);
  TString side;
  Int_t newRPC=0,newPlane=0;
  
  AliCDBManager *manager = AliCDBManager::Instance();
  AliCDBPath path("MUON/Calib/TriggerDCS");
  AliCDBEntry *entry = manager->Get(path);
  if (entry == NULL) {
    AliWarning("No map found in MUON/Calib/TriggerDCS");
    return;
  }
  TMap *hvMap = dynamic_cast<TMap*>(entry->GetObject());
  TObjArray *objArr = 0x0;
  
  AliDCSValue *dcsValue = 0x0;
  UInt_t time1,time2,timebegin=0,timeend=0;
  Int_t nEntries;
  Float_t voltage = 0;
  
  for(Int_t iPlane=0; iPlane<4; iPlane++) //loop on MT
  {
    for(Int_t iRPC=0; iRPC<18; iRPC++) //loop on RPC
    {
      if(iRPC>=5 && iRPC<=13)
      {
        side = "OUTSIDE";
        newRPC = 14-iRPC;
      }
  
      else
      {
        side = "INSIDE";
    
        if(iRPC>=14)
          newRPC = iRPC-13;
        else
          newRPC = iRPC+5;
      }
  
      switch(iPlane)
      {
        case 0: newPlane = 11; break;
        case 1: newPlane = 12; break;
        case 2: newPlane = 21; break;
        case 3: newPlane = 22; break;
      }
  
      objArr = (TObjArray*)hvMap->GetValue(Form("MTR_%s_MT%d_RPC%d_HV.vEff",side.Data(),newPlane,newRPC));
      nEntries = objArr->GetEntries();
	
      for(Int_t i=0; i<nEntries-1; i++)
      {	  
        dcsValue = (AliDCSValue*)objArr->At(i+1);
        time2 = dcsValue->GetTimeStamp();
    
        if(i==nEntries-2)
          timeend = time2;
    
        dcsValue = (AliDCSValue*)objArr->At(i);
          time1 = dcsValue->GetTimeStamp();
    
        if(i==0)
          timebegin = time1;
    
        voltage += (dcsValue->GetFloat())*(time2-time1);
      }
      
      Double_t deltaTime = timeend - timebegin;
      Double_t meanVoltage = ( deltaTime == 0. ) ? 0. : voltage/deltaTime/1000.;
      fHVvalues.AddAt(meanVoltage,18*iPlane+iRPC); //voltage in kV, not in V
      
      voltage=0;
      AliDebug(1,Form("HV value for MTR_%s_MT%d_RPC%d_HV.vEff = %g (kV)",side.Data(),newPlane,newRPC,meanVoltage));
    }
  }
}

//------------------------------------------------------------------  
void AliMUONResponseTriggerV1::SetBValues()
{
  //
  /// Set B values for cluster size function
  //
  
  fBValues.Set(144);
  
  Float_t bValues[2][4][18] =                												              {{{1.97,2.47,2.47,2.47,2.97,2.97,2.47,2.47,1.97,2.22,1.97,2.47,1.97,2.97,2.97,2.47,2.47,1.97},  //MT11BP
    {2.22,2.22,1.97,2.47,2.97,2.97,1.97,2.47,1.97,1.97,1.97,2.47,1.97,2.97,2.97,1.97,1.97,1.97},  //MT12BP
    {2.22,2.22,2.47,2.47,2.97,2.97,2.47,2.47,2.22,1.97,1.97,2.47,1.97,2.97,2.97,1.97,1.97,1.97},  //MT21BP
    {1.97,1.97,2.97,2.97,2.97,2.97,2.47,1.97,1.97,1.97,1.72,2.47,2.22,2.97,2.97,1.97,1.97,1.97}}, //MT22BP
   {{1.97,2.47,2.47,2.97,2.97,2.97,2.97,2.47,1.97,1.97,2.22,2.47,2.97,2.97,2.97,2.97,1.97,1.72},  //MT11NBP
    {2.47,1.97,2.22,2.97,2.97,2.97,2.47,2.97,1.97,1.97,1.97,2.97,2.97,2.97,2.97,2.97,1.97,1.97},  //MT12NBP
    {1.97,2.47,2.47,2.97,2.97,2.97,2.97,2.47,2.22,1.97,2.22,2.47,2.97,2.97,2.97,2.47,1.97,1.97},  //MT21NBP
    {1.72,1.97,2.97,2.97,2.97,2.97,2.97,1.97,1.72,2.22,1.97,2.47,2.97,2.47,2.97,1.97,1.97,1.97}}};//MT22NBP
                                
  for(Int_t iCath=0; iCath<2; iCath++) //loop on side
  {
    for(Int_t iPlane=0; iPlane<4; iPlane++) //loop on MT
    {
      for(Int_t iRPC=0; iRPC<18; iRPC++) //loop on RPC
      {
	fBValues.AddAt(bValues[iCath][iPlane][iRPC],72*iCath+18*iPlane+iRPC);
      }
    }
  }
}

//------------------------------------------------------------------  
void AliMUONResponseTriggerV1::Neighbours(const Int_t cath, const Int_t ix, const Int_t iy, Int_t Xlist[30], Int_t Ylist[30]) const
{
  ///-----------------BENDING-----------------------------------------      /n
  /// Returns list of 30 next neighbours for given X strip (ix, iy)         /n
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
  /// Returns list of 30 next neighbours for given Y strip (ix, iy)         /n 
  /// neighbour number 9 8 7 6 5 (Y strip (ix, iy)) 0 1 2 3 4 in the list   /n 
  ///                  |_______|                    |_______|               /n 
  ///                    left                         right                 /n
  
  for (Int_t i=0; i<30; i++)
  {
    Xlist[i] = -1;
    Ylist[i] = -1;
  }
  
  Int_t iList[30]={29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
  
  // need to iterate in iy/ix for bending/non-bending plane
  Int_t iNeigh = (cath == 0) ? iy : ix;
    
  Int_t i=0;
  for (Int_t j=iNeigh-15; j<=iNeigh+15; j++)
  {
    if (j == iNeigh)
      continue;
	
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
