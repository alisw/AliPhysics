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

// -------------------------------------
// Class AliMUONTriggerSegmentation
// -------------------------------------
// Segmentation for MUON trigger stations using 
// the mapping package

#include "AliMUONTriggerSegmentation.h"

#include "AliMpPCB.h"
#include "AliMpTrigger.h"
#include "AliMpTriggerSegmentation.h"
#include "AliMpSlat.h"

#include "AliLog.h"

#include "Riostream.h"
#include "TClass.h"
#include "TString.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerSegmentation)
/// \endcond

namespace
{
//  Int_t SPECIAL_SECTOR = 8;
  Int_t fgIntOffset(1);
  Float_t FMAX(1E9);
  Int_t CODEMAKER(1000);
  
  Int_t Code(Int_t ixLA, Int_t iyLA)
  {
    return iyLA*CODEMAKER + ixLA;
  }
  
  void Decode(Int_t code, Int_t& ixLA, Int_t& iyLA)
  {
    iyLA = code/CODEMAKER;
    ixLA = code - iyLA*CODEMAKER;
  }
}

//_____________________________________________________________________________
AliMUONTriggerSegmentation::AliMUONTriggerSegmentation() 
: AliMUONVGeometryDESegmentation(),
  fDetElemId(-1),
  fPlaneType(kNonBendingPlane),
  fSlat(0),
  fSlatSegmentation(0),
  fCurrentPad(),
  fXhit(FMAX),
  fYhit(FMAX),
  fLineNumber(-1)
{
/// Default ctor (empty).

  AliDebug(1,Form("this=%p default (empty) ctor",this));
}

//_____________________________________________________________________________
AliMUONTriggerSegmentation::AliMUONTriggerSegmentation(
                                   AliMpVSegmentation* segmentation,
                                   Int_t detElemId, AliMpPlaneType bendingOrNonBending)
    : AliMUONVGeometryDESegmentation(),
      fDetElemId(detElemId),
      fPlaneType(bendingOrNonBending),
      fSlat(0),
      fSlatSegmentation(0),
      fCurrentPad(),
      fXhit(FMAX),
      fYhit(FMAX),
      fLineNumber(-1)
{
/// Normal ctor.

  fSlatSegmentation = dynamic_cast<AliMpTriggerSegmentation*>(segmentation);
  if (fSlatSegmentation)
    fSlat = fSlatSegmentation->Slat();
  else 
    AliFatal("Wrong mapping segmentation type");
		
  TString id(fSlat->GetID());
  Ssiz_t pos = id.Last('L');
  if ( pos <= 0 )
  {
    AliFatal(Form("Cannot infer line number for slat %s",id.Data()));
  }
  fLineNumber = TString(id(pos+1),1).Atoi();
		
  AliDebug(1,Form("this=%p detElemId=%3d %s fSlatSegmentation=%p",this,detElemId,
									( (bendingOrNonBending==kBendingPlane)?"Bending":"NonBending" ),
									fSlatSegmentation));
}

//_____________________________________________________________________________
AliMUONTriggerSegmentation::~AliMUONTriggerSegmentation() 
{ 
/// Destructor

  AliDebug(1,Form("this=%p",this));			
}

//_____________________________________________________________________________
TF1*
AliMUONTriggerSegmentation::CorrFunc(Int_t) const
{
/// Not implemented

  AliFatal("Not implemented");
  return 0x0;
}

//_____________________________________________________________________________
Float_t
AliMUONTriggerSegmentation::Distance2AndOffset(Int_t, Int_t, 
                                                 Float_t, Float_t, Int_t*)
{
/// Not implemented

  AliFatal("Not implemented");
  return 0;
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::Draw(Option_t*)
{
/// Not implemented

  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
Float_t
AliMUONTriggerSegmentation::Dpx() const
{
/// Not implemented

  AliFatal("Not Implemented");
  return 0.0;
}

//_____________________________________________________________________________
Float_t
AliMUONTriggerSegmentation::Dpy() const
{
/// Not implemented

  AliFatal("Not Implemented");
  return 0.0;
}

//_____________________________________________________________________________
Float_t
AliMUONTriggerSegmentation::Dpx(Int_t sectorCode) const
{
/// Get pad size in x

 Int_t ixLA, iyLA;
  Decode(sectorCode,ixLA,iyLA);
  AliMpPad pad = fSlatSegmentation->PadByIndices(AliMpIntPair(ixLA,iyLA),kFALSE);
  if ( !pad.IsValid() ) return 0.0;
  return pad.Dimensions().X()*2.0;
}

//_____________________________________________________________________________
Float_t
AliMUONTriggerSegmentation::Dpy(Int_t sectorCode) const
{
/// Get pad size in y

  Int_t ixLA, iyLA;
  Decode(sectorCode,ixLA,iyLA);
  AliMpPad pad = fSlatSegmentation->PadByIndices(AliMpIntPair(ixLA,iyLA),kFALSE);
  if ( !pad.IsValid() ) return 0.0;
  return pad.Dimensions().Y()*2.0;
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::FirstPad(Float_t /*xhit*/, Float_t /*yhit*/, Float_t /*zhit*/,
                                       Float_t /*dx*/, Float_t /*dy*/)
{
/// Not implemented

  AliFatal("Not implemented");
}

//_____________________________________________________________________________
Float_t
AliMUONTriggerSegmentation::GetAnod(Float_t) const
{
/// Not implemented

  AliFatal("Not implemented");
  return 0.0;
}

//_____________________________________________________________________________
AliMUONGeometryDirection
AliMUONTriggerSegmentation::GetDirection()
{
/// Not implemented

  //AliWarning("Not Implemented");
  return kDirUndefined;
}

//______________________________________________________________________________
const AliMpVSegmentation*  
AliMUONTriggerSegmentation::GetMpSegmentation() const
{
/// Returns the mapping segmentation
/// (provides access to electronics info)

  return fSlatSegmentation;
}  

//_____________________________________________________________________________
void 
AliMUONTriggerSegmentation::GetNParallelAndOffset(Int_t,Int_t,Int_t*,Int_t*)
{
/// Not implemented

  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::GetPadC(Int_t ix, Int_t iy, 
                                      Float_t& x, Float_t& y, Float_t& z)
{
/// Transform from pad to real coordinates

  z = 0;
  GetPadC(ix,iy,x,y);
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::GetPadC(Int_t ixGlo, Int_t iyGlo, 
                                      Float_t& x, Float_t& y)
{
/// Transform from pad to real coordinates

  Int_t ixLA,iyLA;
  IGlo2ILoc(ixGlo,iyGlo,ixLA,iyLA);
  AliMpPad pad = fSlatSegmentation->PadByIndices(AliMpIntPair(ixLA,iyLA),kTRUE);
  x = pad.Position().X();
  y = pad.Position().Y();
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::GetPadI(Float_t x, Float_t y, Float_t,
                                      Int_t& ix, Int_t& iy)
{
///  Returns pad coordinates (ix,iy) for given real coordinates (x,y)

  GetPadI(x,y,ix,iy);
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::GetPadI(Float_t x, Float_t y,
                                      Int_t& ixGlo, Int_t& iyGlo)
{
///  Returns pad coordinates (ix,iy) for given real coordinates (x,y)

  AliDebug(2,Form("%s x=%e y=%e ixGlo,iyGlo=%d,%d\n",
                  fSlatSegmentation->GetName(),
                  x,y,ixGlo,iyGlo));
  
  AliMpPad pad = 
    fSlatSegmentation->PadByPosition(TVector2(x,y), kTRUE);
	
  if ( pad != AliMpPad::Invalid() )
	{
		Int_t ix = pad.GetIndices().GetFirst();
		Int_t iy = pad.GetIndices().GetSecond();
    ILoc2IGlo(ix,iy,ixGlo,iyGlo);
	}
  else
	{
		ixGlo=iyGlo=-1;
	}
  AliDebug(2,Form("ixGlo,iyGlo=%d,%d\n",ixGlo,iyGlo));
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::GetPadLoc2Glo(Int_t ix, Int_t iy,
                                            Int_t& ixGlo, Int_t& iyGlo) const
{
/// Converts from local (in PC convention) to (ix,iy) to global (ix,iy)

  ixGlo=iyGlo=-1; // starts with invalid values
  
  if ( fPlaneType == kBendingPlane )
  {
    ixGlo = 10*LineNumber() + ix;
    iyGlo = iy - fgIntOffset;
  }
  else if ( fPlaneType == kNonBendingPlane )
  {
    Int_t i = fSlat->GetLayer(0)->FindPCBIndex(ix-fgIntOffset);
    if (i<0)
    {
      AliError(Form("Invalid local (ix=%d,iy=%d) ?",ix,iy));
      return ;
    }
    AliMpPCB* pcb = fSlat->GetLayer(0)->FindPCB(ix-fgIntOffset);
    iyGlo = ix - pcb->Ixmin() - fgIntOffset;
    if ( LineNumber() == 5 ) ++i;
    ixGlo = 10*LineNumber() + i + fgIntOffset; 
  }
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::GetPadGlo2Loc(Int_t ixGlo, Int_t iyGlo,
                                            Int_t& ix, Int_t& iy) const
{
/// Converts from global (ix,iy) to local (ix,iy) (in PC convention)
  
  ix=iy=-1; // starts with invalid values

  if ( abs(ixGlo) == 51 ) return;
  
  Int_t column = ModuleColNum(ixGlo);

  if ( fPlaneType == kBendingPlane )
  {
    ix = column + fgIntOffset;
    iy = iyGlo + fgIntOffset;
  }
  else if ( fPlaneType == kNonBendingPlane )
  {
    if ( LineNumber()==5 ) --column;
    AliMpPCB* pcb = fSlat->GetLayer(0)->GetPCB(column);
    if (!pcb)
    {
      AliError(Form("Invalid global (ix=%d,iy=%d)",ixGlo,iyGlo));
      return;
    }
    ix = pcb->Ixmin() + iyGlo + fgIntOffset;
    iy = fgIntOffset;
  }
}

//_____________________________________________________________________________
void 
AliMUONTriggerSegmentation::GiveTestPoints(Int_t& /*n*/, 
                                             Float_t* /*x*/, Float_t*/*y*/) const
{
// Not implemented

  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
Bool_t
AliMUONTriggerSegmentation::HasPad(Float_t x, Float_t y, Float_t /*z*/)
{
/// Returns true if a pad exists in the given position
///
/// Well, 2 implementations are possible here
/// Either reuse HasPad(int,int), or get it from scratch using
/// underlying fSlatSegmentation.
/// Took second option, but w/o checking whether this is the faster.
/// The second option is commented out below, for the record.
  
//  Int_t ix, iy;
//  GetPadI(x,y,z,ix,iy);
//  Int_t ixLA, iyLA;
//  IGlo2ILoc(ix,iy,ixLA,iyLA);
//  Int_t ixPC, iyPC;
//  LA2PC(ixLA,iyLA,ixPC,iyPC);
//  Bool_t ok1 = HasPad(ixPC,iyPC);

  AliMpPad pad = 
  fSlatSegmentation->PadByPosition(TVector2(x,y),kFALSE);
  return pad.IsValid();
}

//_____________________________________________________________________________
Bool_t
AliMUONTriggerSegmentation::HasPad(Int_t ixGlo, Int_t iyGlo)
{
/// Returns true if a pad with given indices exists

  Int_t ixLA, iyLA;
  IGlo2ILoc(ixGlo,iyGlo,ixLA,iyLA);
  return fSlatSegmentation->HasPad(AliMpIntPair(ixLA,iyLA));
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::IGlo2ILoc(Int_t ixGlo, Int_t iyGlo,
                                        Int_t& ixLA, Int_t& iyLA) const
{
/// \todo FIXME: add comment

  Int_t ixPC, iyPC;
  GetPadGlo2Loc(ixGlo,iyGlo,ixPC,iyPC);
  PC2LA(ixPC,iyPC,ixLA,iyLA);
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::ILoc2IGlo(Int_t ixLA, Int_t iyLA,
                                        Int_t& ixGlo, Int_t& iyGlo) const
{
/// \todo FIXME: add comment

  Int_t ixPC, iyPC;
  LA2PC(ixLA,iyLA,ixPC,iyPC);
  GetPadLoc2Glo(ixPC,iyPC,ixGlo,iyGlo);
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentation::ISector()
{
/// \todo FIXME: remove the usage of ISector from all the code.

  return -10;
}

//_____________________________________________________________________________
void AliMUONTriggerSegmentation::IntegrationLimits(Float_t& x1, 
                                                     Float_t& x2,
                                                     Float_t& x3, 
                                                     Float_t& x4) 
{
/// \param x1 : hit x(y) position
/// \param  x2 : x(y) coordinate of the main strip
/// \param  x3 : current strip real x(y) coordinate  
/// \param  x4 : dist. between x(y) hit pos. and the closest border of the current strip
///
/// Note : need to return (only) x4.
  
  AliFatal("Check me before usage. ResponseTrigger does not use me, while"
           "ResponseTriggerV1 does ?");
    
  AliMpPad strip = 
  fSlatSegmentation->PadByPosition(TVector2(fXhit,fYhit),kFALSE);
  if ( !strip.IsValid() )
  {
    AliWarning(Form("%s got invalid fXhit,fYhit=%e,%e\n",
                    fSlatSegmentation->GetName(),fXhit,fYhit));
    x1=x2=x3=x4=0;
  }
  else
  {
    Double_t xstrip = strip.Position().X();
    Double_t ystrip = strip.Position().Y();
    AliDebug(1,Form("fXhit,Yhit=%e,%e xstrip,ystrip=%e,%e\n",
                    fXhit,fYhit,xstrip,ystrip));
    x1 = (fPlaneType==kBendingPlane) ? fYhit : fXhit;
    x2 = (fPlaneType==kBendingPlane) ? ystrip : xstrip;
    x3 = (fPlaneType==kBendingPlane) ? 
      fCurrentPad.Position().Y() : fCurrentPad.Position().X();
    Double_t xmin = 0.0;
    Double_t xmax = 0.0;
    if (fPlaneType==kBendingPlane)
    {
      xmin = x3 - fCurrentPad.Dimensions().X();
      xmax = x3 + fCurrentPad.Dimensions().X();
    }
    else
    {
      xmin = x3 - fCurrentPad.Dimensions().Y();
      xmax = x3 + fCurrentPad.Dimensions().Y();
    }
    // dist. between the hit and the closest border of the current strip
    x4 = (TMath::Abs(xmax-x1) > TMath::Abs(xmin-x1)) ? 
      TMath::Abs(xmin-x1):TMath::Abs(xmax-x1);    
    
    AliDebug(1,Form("Bending %d x1=%e x2=%e x3=%e x4=%e xmin,max=%e,%e\n",
                    fPlaneType,x1,x2,x3,x4,xmin,xmax));

  }
}  

//_____________________________________________________________________________
Int_t 
AliMUONTriggerSegmentation::Ix()
{
/// Current pad cursor during disintegration
/// x, y-coordinate

  if ( fCurrentPad.IsValid() )
  {
    Int_t ixGlo,iyGlo;
    ILoc2IGlo(fCurrentPad.GetIndices().GetFirst(),
              fCurrentPad.GetIndices().GetSecond(),ixGlo,iyGlo);
    return ixGlo;
  }
  return -1;
}

//_____________________________________________________________________________
Int_t 
AliMUONTriggerSegmentation::Iy()
{
/// Current pad cursor during disintegration
/// x, y-coordinate

  if ( fCurrentPad.IsValid() )
  {
    Int_t ixGlo,iyGlo;
    ILoc2IGlo(fCurrentPad.GetIndices().GetFirst(),
              fCurrentPad.GetIndices().GetSecond(),ixGlo,iyGlo);
    return iyGlo;
  }
  return -1;
}


//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::LA2PC(Int_t ixLA, Int_t iyLA,
                                    Int_t& ixPC, Int_t& iyPC) const
{
/// From LA to PC conventions for integers indices.

  ixPC=iyPC=-1;
  
  if ( ixLA<0 || iyLA<0 ) return;
  
  ixPC = ixLA + 1;
  iyPC = iyLA + 1;
  
  if ( fPlaneType == kBendingPlane )
  {
    if ( LineNumber()==5 )
    {
      ++ixPC;
    }
    if ( LineNumber()==4 && ixLA==0 )
    {
      iyPC -= 16;
    }
  }
  
  AliDebug(3,Form("ix,iy LA (%d,%d) -> PC (%d,%d)",ixLA,iyLA,ixPC,iyPC));
  
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentation::LineNumber() const
{
/// \todo FIXME: add comment

 return 10-fLineNumber;
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentation::ModuleColNum(Int_t ixGlo) const
{
/// returns column number (from 0 to 6) in which the (global) module 
/// ix is sitting (could return 7 if ix=isec)

  return TMath::Abs(ixGlo)-Int_t(TMath::Abs(ixGlo)/10)*10-1;
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentation::MorePads()
{
/// Not implemented

  AliFatal("Not implemented");
  return 0;
}

//_____________________________________________________________________________
void AliMUONTriggerSegmentation::Neighbours(Int_t /*iX*/, Int_t /*iY*/, 
                                              Int_t* /*Nlist*/, 
                                              Int_t /*Xlist*/[10], 
                                              Int_t /*Ylist*/[10]) 
{
/// Not implemented

  //-----------------BENDING-----------------------------------------
  // Returns list of 10 next neighbours for given X strip (ix, iy)  
  // neighbour number 4 in the list -                     
  // neighbour number 3 in the list  |                    
  // neighbour number 2 in the list  |_ Upper part             
  // neighbour number 1 in the list  |            
  // neighbour number 0 in the list -           
  //      X strip (ix, iy) 
  // neighbour number 5 in the list -       
  // neighbour number 6 in the list  | _ Lower part
  // neighbour number 7 in the list  |
  // neighbour number 8 in the list  | 
  // neighbour number 9 in the list -
  
  //-----------------NON-BENDING-------------------------------------
  // Returns list of 10 next neighbours for given Y strip (ix, iy)  
  // neighbour number 9 8 7 6 5 (Y strip (ix, iy)) 0 1 2 3 4 in the list
  //                  \_______/                    \_______/
  //                    left                         right
  AliFatal("Please implement me");
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::NextPad()
{
/// Not implemented

  AliFatal("Not implemented");
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentation::Npx() const
{
/// Maximum number of Pads in x
/// hard coded for the time being

  return 124;// FIXME: this should not have to be done, if only we'd stick 
  // to a local (ix,iy) convention !!! 
  // return fSlatSegmentation->MaxPadIndexX()+1;
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentation::Npy() const
{
/// Maximum number of Pads in y
/// hard coded for the time being

  return 64;
//  return fSlatSegmentation->MaxPadIndexY()+1;
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::PC2LA(Int_t ixPC, Int_t iyPC,
                                    Int_t& ixLA, Int_t& iyLA) const
{
/// From PC to LA conventions for integers indices.

  ixLA=iyLA=-1;
  
  if ( ixPC<0 || iyPC<0 ) return;
  
  ixLA = ixPC - 1;
  iyLA = iyPC - 1;
  
  if ( fPlaneType == kBendingPlane )
  {
    if ( LineNumber()==5 )
    {
      --ixLA;
    }
    if ( LineNumber()==4 && ixLA==0 )
    {
      iyLA += 16;
    }
  }
  
  AliDebug(3,Form("ix,iy PC (%d,%d) -> LA (%d,%d)",ixPC,iyPC,ixLA,iyLA));
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::Print(Option_t* opt) const
{
/// Printing

  TString sopt(opt);
  
  cout << "DetElemId=" << fDetElemId << " PlaneType=" 
    << fPlaneType << " Npx=" << Npx() << " Npy=" << Npy() << endl;
  if ( ( sopt.Contains("SEG") || sopt.Contains("ALL") ) && fSlatSegmentation )
  {
    fSlatSegmentation->Print();
  }
  if ( sopt.Contains("ALL") && fSlat )
  {
    fSlat->Print();
  }
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentation::Sector(Int_t ix, Int_t iy)
{
/// Calculate sector from pad coordinates

  Int_t ixLA, iyLA;
  IGlo2ILoc(ix,iy,ixLA,iyLA);
  AliMpPad pad = fSlatSegmentation->PadByIndices(AliMpIntPair(ixLA,iyLA),kFALSE);
  if ( !pad.IsValid() ) return -1;
  return Code(ixLA,iyLA);
  
//  AliMpPCB* pcb = fSlat->GetLayer(0)->FindPCB(ixLA);
//  if (!pcb)
//  {
//    AliError(Form("Could not find a pcb at (%d,%d) for slat %s",
//             ix,iy,fSlat->GetName()));
//    return -1;
//  }
//  if ( pcb->PadSizeX()==-1.0 )
//  {
//    // special case of column 7 non-bending.
//    return SPECIAL_SECTOR;
//  }
//  return fSlat->GetLayer(0)->FindPCBIndex(ixLA);
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentation::Sector(Float_t, Float_t)
{
/// Not implemented

  AliFatal("Not implemented");
  return 0;
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::SetCorrFunc(Int_t,TF1*)
{
/// Not implemented

  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::SetDAnod(Float_t)
{
/// Not implemented

  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::SetHit(Float_t x, Float_t y)
{
/// Set hit position
/// Sets virtual hit position, needed for evaluating pad response 
/// outside the tracking program 

  fXhit = x;
  fYhit = y;
	
  //
  // insure we're within the slat limits. If not, issue an error and sets
  // the current hit to slat center.
  // FIXME: this should probably a) not happen at all b) be a fatal error
  //
  if ( fXhit < -fSlat->DX() || fXhit > fSlat->DX() ||
       fYhit < -fSlat->DY() || fYhit > fSlat->DY() )
	{
		AliError(Form("Hit outside slat %s limits (x,y)hit = (%e,%e)."
                  " Forcing to (0,0)",fSlat->GetID(),fXhit,fYhit));
		fXhit = 0.0;
		fYhit = 0.0;
	}
  
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::SetHit(Float_t x, Float_t y, Float_t)
{
/// Set hit position
/// Sets virtual hit position, needed for evaluating pad response 
/// outside the tracking program 

  SetHit(x,y);
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::SetPad(Int_t ix, Int_t iy)
{
/// Set pad position.
/// Sets virtual pad coordinates, needed for evaluating pad response 
/// outside the tracking program.

  Int_t ixLA, iyLA;
  IGlo2ILoc(ix,iy,ixLA,iyLA);
  fCurrentPad = fSlatSegmentation->PadByIndices(AliMpIntPair(ixLA,iyLA),kTRUE);
  if ( !fCurrentPad.IsValid() )
	{
		AliError(Form("Setting current pad to invalid ! (ix,iy)=(%4d,%4d)",ix,iy));
	}
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentation::SetPadSize(Float_t,Float_t)
{
/// Not implemented

  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
Int_t 
AliMUONTriggerSegmentation::SigGenCond(Float_t,Float_t,Float_t)
{
/// Not implemented

  AliFatal("Not Implemented");
  return 0;
}

//_____________________________________________________________________________
void 
AliMUONTriggerSegmentation::SigGenInit(Float_t,Float_t,Float_t)
{
/// Not implemented

  AliFatal("Not Implemented");
}





