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
// Class AliMUONResponseV0
// --------------------------
// Implementation of 
// Mathieson response
//-----------------------------------------------------------------------------

#include "AliMUONResponseV0.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONDigit.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONHit.h"
#include "AliMUONConstants.h"

#include "AliMpArea.h"
#include "AliMpDEManager.h"
#include "AliMpVPadIterator.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "AliMpCathodType.h"

#include "AliRun.h"
#include "AliLog.h"

#include "Riostream.h"
#include "TVector2.h"
#include <TMath.h>
#include <TRandom.h>

/// \cond CLASSIMP
ClassImp(AliMUONResponseV0)
/// \endcond
	
AliMUON* muon()
{
    return static_cast<AliMUON*>(gAlice->GetModule("MUON"));
}

void Global2Local(Int_t detElemId, Double_t xg, Double_t yg, Double_t zg,
                  Double_t& xl, Double_t& yl, Double_t& zl)
{  
  /// ideally should be : 
  /// Double_t x,y,z;
  /// AliMUONGeometry::Global2Local(detElemId,xg,yg,zg,x,y,z);
  /// but while waiting for this geometry singleton, let's go through
  /// AliMUON still.
  
  const AliMUONGeometryTransformer* transformer = muon()->GetGeometryTransformer();
  transformer->Global2Local(detElemId,xg,yg,zg,xl,yl,zl);
}

//__________________________________________________________________________
AliMUONResponseV0::AliMUONResponseV0()
  : AliMUONResponse(),
  fChargeSlope(0.0),
  fChargeSpreadX(0.0),
  fChargeSpreadY(0.0),
  fSigmaIntegration(0.0),
  fMaxAdc(0),
  fSaturation(0),
  fZeroSuppression(0),
  fChargeCorrel(0.0),
  fMathieson(new AliMUONMathieson),
  fChargeThreshold(1e-4),
  fIsTailEffect(kFALSE)
{
    /// Normal constructor
    AliDebug(1,Form("Default ctor"));
}

//__________________________________________________________________________
AliMUONResponseV0::AliMUONResponseV0(const AliMUONResponseV0& other)
: AliMUONResponse(),
fChargeSlope(0.0),
fChargeSpreadX(0.0),
fChargeSpreadY(0.0),
fSigmaIntegration(0.0),
fMaxAdc(0),
fSaturation(0),
fZeroSuppression(0),
fChargeCorrel(0.0),
fMathieson(0),
fChargeThreshold(1e-4),
fIsTailEffect(kFALSE)
{
  /// copy ctor
  other.CopyTo(*this);
}

//__________________________________________________________________________
AliMUONResponseV0& 
AliMUONResponseV0::operator=(const AliMUONResponseV0& other)
{
  /// Assignment operator
  other.CopyTo(*this);
  return *this;
}

//__________________________________________________________________________
AliMUONResponseV0::~AliMUONResponseV0()
{
/// Destructor

  AliDebug(1,"");
  delete fMathieson;
}

//______________________________________________________________________________
void
AliMUONResponseV0::CopyTo(AliMUONResponseV0& other) const
{
  /// Copy *this to other
  other.fChargeSlope=fChargeSlope;
  other.fChargeSpreadX=fChargeSpreadX;
  other.fChargeSpreadY=fChargeSpreadY;
  other.fSigmaIntegration=fSigmaIntegration;
  other.fMaxAdc=fMaxAdc;
  other.fSaturation=fSaturation;
  other.fZeroSuppression=fZeroSuppression;
  other.fChargeCorrel=fChargeCorrel;
  delete other.fMathieson;
  other.fMathieson = new AliMUONMathieson(*fMathieson);
  other.fChargeThreshold=fChargeThreshold;
}

//______________________________________________________________________________
void
AliMUONResponseV0::Print(Option_t*) const
{
/// Printing

  cout << " ChargeSlope=" << fChargeSlope
    << " ChargeSpreadX,Y=" << fChargeSpreadX
    << fChargeSpreadY
    << " ChargeCorrelation=" << fChargeCorrel
    << endl;
}

  //__________________________________________________________________________
void AliMUONResponseV0::SetSqrtKx3AndDeriveKx2Kx4(Float_t SqrtKx3)
{
  /// Set to "SqrtKx3" the Mathieson parameter K3 ("fSqrtKx3")
  /// in the X direction, perpendicular to the wires,
  /// and derive the Mathieson parameters K2 ("fKx2") and K4 ("fKx4")
  /// in the same direction
  fMathieson->SetSqrtKx3AndDeriveKx2Kx4(SqrtKx3);
}
	
  //__________________________________________________________________________
void AliMUONResponseV0::SetSqrtKy3AndDeriveKy2Ky4(Float_t SqrtKy3)
{
  /// Set to "SqrtKy3" the Mathieson parameter K3 ("fSqrtKy3")
  /// in the Y direction, along the wires,
  /// and derive the Mathieson parameters K2 ("fKy2") and K4 ("fKy4")
  /// in the same direction
  fMathieson->SetSqrtKy3AndDeriveKy2Ky4(SqrtKy3);
}
  //__________________________________________________________________________
Float_t AliMUONResponseV0::IntPH(Float_t eloss) const
{
  /// Calculate charge from given ionization energy loss
  Int_t nel;
  nel= Int_t(eloss*1.e9/27.4);
  Float_t charge=0;
  if (nel == 0) nel=1;
  for (Int_t i=1;i<=nel;i++) {
      Float_t arg=0.;
      while(!arg) arg = gRandom->Rndm();
      charge -= fChargeSlope*TMath::Log(arg);    
  }
  return charge;
}

//_____________________________________________________________________________
Float_t
AliMUONResponseV0::GetAnod(Float_t x) const
{
  /// Return wire coordinate closest to x.

  Int_t n = Int_t(x/Pitch());
  Float_t wire = (x>0) ? n+0.5 : n-0.5;
  return Pitch()*wire;
}

//______________________________________________________________________________
void 
AliMUONResponseV0::DisIntegrate(const AliMUONHit& hit, TList& digits, Float_t timeDif)
{
  /// Go from 1 hit to a list of digits.
  /// The energy deposition of that hit is first converted into charge
  /// (in IntPH() method), and then this charge is dispatched on several
  /// pads, according to the Mathieson distribution.
  
  digits.Clear();
  
  Int_t detElemId = hit.DetElemId();
  Double_t hitX = hit.X() ;
  Double_t hitY = hit.Y() ;
  Double_t hitZ = hit.Z() ;

  // Width of the integration area
  Double_t dx = SigmaIntegration()*ChargeSpreadX();
  Double_t dy = SigmaIntegration()*ChargeSpreadY();
  
  //Modify to take the tailing effect.
  if(fIsTailEffect){
    Double_t locX,locY,locZ,globXCenter,globYCenter,globZ;
    Int_t para = 5; // This parameter is a natural number(excluding zero), higher the value less is the tailing effect 
    Double_t termA = 1.0;
    Double_t termB = 1.0;
    if(para>0){
      for ( Int_t cath = AliMp::kCath0; cath <= AliMp::kCath1; ++cath )
	{
	  // Get an iterator to loop over pads, within the given area.
	  const AliMpVSegmentation* seg = 
	    AliMpSegmentation::Instance()
	    ->GetMpSegmentation(detElemId,AliMp::GetCathodType(cath));
	  AliMp::PlaneType plane = seg->PlaneType();
	  
	  if(plane == AliMp::kBendingPlane) {
	    Global2Local(detElemId,hitX,hitY,hitZ,locX,locY,locZ);
	    AliMpPad pad = seg->PadByPosition(locX,locY,kFALSE);
            if(pad.IsValid()){
              Double_t locYCenter = pad.GetPositionY();
              Double_t locXCenter = pad.GetPositionX();
              const AliMUONGeometryTransformer* transformer = muon()->GetGeometryTransformer();
              transformer->Local2Global(detElemId,locXCenter,locYCenter,locZ,globXCenter,globYCenter,globZ);
              for(Int_t itime = 0; itime<para; itime++)
                termA *= 10.0;
            
              for(Int_t itime = 0; itime<Int_t((2*para) + 1); itime++)
                termB *= (hitY - globYCenter) ; 
            
              hitY = hitY + termA*termB;
            }// if the pad is a valid one
	  }// if bending plane
	}// cathode loop
    }// if para > 0 condn
  }// if tail effect

  // Use that (dx,dy) to specify the area upon which
  // we will iterate to spread charge into.
  Double_t x,y,z;
  Global2Local(detElemId,hitX,hitY,hitZ,x,y,z);
  x = GetAnod(x);
  AliMpArea area(x,y,dx,dy);
  
  // Get pulse height from energy loss.
  Float_t qtot = IntPH(hit.Eloss());
  
  // If from a pileup event we apply a reduction factor to the charge
  if (timeDif!=0){
    qtot = AliMUONConstants::ReducedQTot(qtot,timeDif);
  }
  
  // Scale the charge to it'll (roughly) be in fC
  qtot *= AliMUONConstants::DefaultADC2MV()*AliMUONConstants::DefaultA0()*AliMUONConstants::DefaultCapa();
  
  // Get the charge correlation between cathodes.
  Float_t currentCorrel = TMath::Exp(gRandom->Gaus(0.0,ChargeCorrel()/2.0));

  for ( Int_t cath = AliMp::kCath0; cath <= AliMp::kCath1; ++cath )
  {
    Float_t qcath = qtot * ( cath == 0 ? currentCorrel : 1.0/currentCorrel);
    
    // Get an iterator to loop over pads, within the given area.
    const AliMpVSegmentation* seg = 
        AliMpSegmentation::Instance()
          ->GetMpSegmentation(detElemId,AliMp::GetCathodType(cath));
      
    AliMpVPadIterator* it = seg->CreateIterator(area);
      
    if (!it)
    {
      AliError(Form("Could not get iterator for detElemId %d",detElemId));
      return;
    }
    
    // Start loop over pads.
    it->First();
    
    if ( it->IsDone() )
    {
      // Exceptional case : iterator is built, but is invalid from the start.
      AliMpPad pad = seg->PadByPosition(area.GetPositionX(),area.GetPositionY(),
                                        kFALSE);
      if ( pad.IsValid() )
      {
        AliDebug(1, Form("Got an invalid iterator bug (area.Position() is within "
                      " DE but the iterator is void) for detElemId %d cath %d",
                      detElemId,cath));        
      }
      else
      {
        AliDebug(1, Form("Got an invalid iterator bug for detElemId %d cath %d."
                      "Might be a bad hit ? area.Position()=(%e,%e) "
                      "Dimensions()=(%e,%e)",
                      detElemId,cath,area.GetPositionX(),area.GetPositionY(),
                      area.GetDimensionX(),area.GetDimensionY()));
      }
      delete it;
      return;
    }
    
    while ( !it->IsDone() )
    {
      // For each pad given by the iterator, compute the charge of that
      // pad, according to the Mathieson distribution.
      AliMpPad pad = it->CurrentItem();      
      TVector2 lowerLeft(TVector2(x,y)-TVector2(pad.GetPositionX(),pad.GetPositionY())-
                         TVector2(pad.GetDimensionX(),pad.GetDimensionY()));
      TVector2 upperRight(lowerLeft + TVector2(pad.GetDimensionX(),pad.GetDimensionY())*2.0);
      Float_t qp = TMath::Abs(fMathieson->IntXY(lowerLeft.X(),lowerLeft.Y(),
                                                upperRight.X(),upperRight.Y()));
            
      if ( qp > fChargeThreshold && 
           qp*qcath > AliMUONConstants::DefaultADC2MV()*AliMUONConstants::DefaultA0()*AliMUONConstants::DefaultCapa() )
      {
        // If we're above threshold, then we create a digit,
        // and fill it with relevant information, including electronics.
        
        // note that the second condition above is to be backward compatible (when 
        // the sdigitizer was making a cut on Int_t(qp*qcath) > 0 and qcath was in ADC, not in fC)
        
        AliMUONDigit* d = new AliMUONDigit(detElemId,pad.GetManuId(),
                                           pad.GetManuChannel(),cath);
        d->SetPadXY(pad.GetIx(),pad.GetIy());
        d->SetCharge(qp*qcath);
        digits.Add(d);   
      }       
      it->Next();
    }
    delete it;
  }
}



