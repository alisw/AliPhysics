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

#include "AliMUONTriggerSegmentationV2.h"

#include "AliMpPCB.h"
#include "AliMpTrigger.h"
#include "AliMpTriggerSegmentation.h"
#include "AliMpSlat.h"

#include "AliLog.h"

#include "Riostream.h"
#include "TClass.h"
#include "TString.h"

ClassImp(AliMUONTriggerSegmentationV2)

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
AliMUONTriggerSegmentationV2::AliMUONTriggerSegmentationV2() 
: AliMUONVGeometryDESegmentation(),
fDetElemId(-1),
fPlaneType(kNonBendingPlane),
fSlat(0),
fSlatSegmentation(0),
fXhit(FMAX),
fYhit(FMAX),
fLineNumber(-1)
{
  //
  // Default ctor (empty).
  //
  AliDebug(1,Form("this=%p default (empty) ctor",this));
}

//_____________________________________________________________________________
AliMUONTriggerSegmentationV2::AliMUONTriggerSegmentationV2(
                                   AliMpVSegmentation* segmentation,
                                   Int_t detElemId, AliMpPlaneType bendingOrNonBending)
: AliMUONVGeometryDESegmentation(),
fDetElemId(detElemId),
fPlaneType(bendingOrNonBending),
fSlat(0),
fSlatSegmentation(0),
fXhit(FMAX),
fYhit(FMAX),
fLineNumber(-1)
{
  //
  // Normal ctor.
  //

  fSlatSegmentation = dynamic_cast<AliMpTriggerSegmentation*>(segmentation);
  if (fSlatSegmentation)
    fSlat = fSlatSegmentation->Slat();
  else 
    AliFatal("Wrong mapping segmentation type");
		
		
  AliDebug(1,Form("this=%p detElemId=%3d %s fSlatSegmentation=%p",this,detElemId,
									( (bendingOrNonBending==kBendingPlane)?"Bending":"NonBending" ),
									fSlatSegmentation));
}

//_____________________________________________________________________________
AliMUONTriggerSegmentationV2::AliMUONTriggerSegmentationV2(const AliMUONTriggerSegmentationV2& rhs) : AliMUONVGeometryDESegmentation(rhs)
{
  AliFatal("Not implemented.");
}

//_____________________________________________________________________________
AliMUONTriggerSegmentationV2::~AliMUONTriggerSegmentationV2() 
{ 
  AliDebug(1,Form("this=%p",this));			
  // Destructor
}

//_____________________________________________________________________________
AliMUONTriggerSegmentationV2& AliMUONTriggerSegmentationV2::operator=(const AliMUONTriggerSegmentationV2& rhs)
{
// Protected assignement operator
  if (this == &rhs) return *this;
  AliFatal("Not implemented.");
  return *this;  
}

//_____________________________________________________________________________
TF1*
AliMUONTriggerSegmentationV2::CorrFunc(Int_t) const
{
  AliFatal("Not implemented");
  return 0x0;
}

//_____________________________________________________________________________
Float_t
AliMUONTriggerSegmentationV2::Distance2AndOffset(Int_t, Int_t, 
                                                 Float_t, Float_t, Int_t*)
{
  AliFatal("Not implemented");
  return 0;
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::Draw(Option_t*)
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
Float_t
AliMUONTriggerSegmentationV2::Dpx() const
{
  AliFatal("Not Implemented");
  return 0.0;
}

//_____________________________________________________________________________
Float_t
AliMUONTriggerSegmentationV2::Dpy() const
{
  AliFatal("Not Implemented");
  return 0.0;
}

//_____________________________________________________________________________
Float_t
AliMUONTriggerSegmentationV2::Dpx(int sectorCode) const
{
  Int_t ixLA, iyLA;
  Decode(sectorCode,ixLA,iyLA);
  AliMpPad pad = fSlatSegmentation->PadByIndices(AliMpIntPair(ixLA,iyLA),kFALSE);
  if ( !pad.IsValid() ) return 0.0;
  return pad.Dimensions().X()*2.0;
}

//_____________________________________________________________________________
Float_t
AliMUONTriggerSegmentationV2::Dpy(int sectorCode) const
{
  Int_t ixLA, iyLA;
  Decode(sectorCode,ixLA,iyLA);
  AliMpPad pad = fSlatSegmentation->PadByIndices(AliMpIntPair(ixLA,iyLA),kFALSE);
  if ( !pad.IsValid() ) return 0.0;
  return pad.Dimensions().Y()*2.0;
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::FirstPad(Float_t /*xhit*/, Float_t /*yhit*/, 
																				 Float_t /*zhit*/, 
																				 Float_t /*dx*/, Float_t /*dy*/)
{
  AliFatal("Not implemented");
}

//_____________________________________________________________________________
Float_t
AliMUONTriggerSegmentationV2::GetAnod(Float_t) const
{
  AliFatal("Not implemented");
  return 0.0;
}

//_____________________________________________________________________________
AliMUONGeometryDirection
AliMUONTriggerSegmentationV2::GetDirection()
{
  //AliWarning("Not Implemented");
  return kDirUndefined;
}

//______________________________________________________________________________
const AliMpVSegmentation*  
AliMUONTriggerSegmentationV2::GetMpSegmentation() const
{
// Returns the mapping segmentation
// (provides access to electronics info)

  return fSlatSegmentation;
}  

//_____________________________________________________________________________
void 
AliMUONTriggerSegmentationV2::GetNParallelAndOffset(Int_t,Int_t,Int_t*,Int_t*)
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::GetPadC(Int_t ix, Int_t iy, 
                                      Float_t& x, Float_t& y, Float_t& z)
{
  z = 0;
  GetPadC(ix,iy,x,y);
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::GetPadC(Int_t ixGlo, Int_t iyGlo, 
                                      Float_t& x, Float_t& y)
{
  Int_t ixLA,iyLA;
  IGlo2ILoc(ixGlo,iyGlo,ixLA,iyLA);
  AliMpPad pad = fSlatSegmentation->PadByIndices(AliMpIntPair(ixLA,iyLA),kTRUE);
  x = pad.Position().X();
  y = pad.Position().Y();
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::GetPadI(Float_t x, Float_t y, Float_t,
                                      Int_t& ix, Int_t& iy)
{
  GetPadI(x,y,ix,iy);
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::GetPadI(Float_t x, Float_t y,
                                      Int_t& ixGlo, Int_t& iyGlo)
{
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
AliMUONTriggerSegmentationV2::GetPadLoc2Glo(Int_t ix, Int_t iy,
                                            Int_t& ixGlo, Int_t& iyGlo) const
{
  //
  // Converts from local (in PC convention) to (ix,iy) to global (ix,iy)
  //

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
AliMUONTriggerSegmentationV2::GetPadGlo2Loc(Int_t ixGlo, Int_t iyGlo,
                                            Int_t& ix, Int_t& iy) const
{
  //
  // Converts from global (ix,iy) to local (ix,iy) (in PC convention)
  //
  
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
AliMUONTriggerSegmentationV2::GiveTestPoints(Int_t&,Float_t*,Float_t*) const
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
Bool_t
AliMUONTriggerSegmentationV2::HasPad(Float_t x, Float_t y, Float_t)
{
  //
  // Well, 2 implementations are possible here
  // Either reuse HasPad(int,int), or get it from scratch using
  // underlying fSlatSegmentation.
  // Took second option, but w/o checking whether this is the faster.
  // The second option is commented out below, for the record.
  
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
AliMUONTriggerSegmentationV2::HasPad(Int_t ixGlo, Int_t iyGlo)
{
  Int_t ixLA, iyLA;
  IGlo2ILoc(ixGlo,iyGlo,ixLA,iyLA);
  return fSlatSegmentation->HasPad(AliMpIntPair(ixLA,iyLA));
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::IGlo2ILoc(Int_t ixGlo, Int_t iyGlo,
                                        Int_t& ixLA, Int_t& iyLA)
{
  Int_t ixPC, iyPC;
  GetPadGlo2Loc(ixGlo,iyGlo,ixPC,iyPC);
  PC2LA(ixPC,iyPC,ixLA,iyLA);
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::ILoc2IGlo(Int_t ixLA, Int_t iyLA,
                                        Int_t& ixGlo, Int_t& iyGlo)
{
  Int_t ixPC, iyPC;
  LA2PC(ixLA,iyLA,ixPC,iyPC);
  GetPadLoc2Glo(ixPC,iyPC,ixGlo,iyGlo);
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::Init(int)
{
  AliWarning("Not Implemented because not needed ;-)");
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentationV2::ISector()
{
  // FIXME: remove the usage of ISector from all the code.
  return -10;
}

//_____________________________________________________________________________
void AliMUONTriggerSegmentationV2::IntegrationLimits(Float_t& x1, 
                                                     Float_t& x2,
                                                     Float_t& x3, 
                                                     Float_t& x4) 
{
  // x1 : hit x(y) position
  // x2 : x(y) coordinate of the main strip
  // x3 : current strip real x(y) coordinate  
  // x4 : dist. between x(y) hit pos. and the closest border of the current strip
  //
  // Note : need to return (only) x4.
  
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
AliMUONTriggerSegmentationV2::Ix()
{
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
AliMUONTriggerSegmentationV2::Iy()
{
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
AliMUONTriggerSegmentationV2::LA2PC(Int_t ixLA, Int_t iyLA,
                                    Int_t& ixPC, Int_t& iyPC)
{
  //
  // From LA to PC conventions for integers indices.
  //
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
AliMUONTriggerSegmentationV2::LineNumber() const
{
  return 10-fLineNumber;
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentationV2::ModuleColNum(Int_t ixGlo) const
{
  // returns column number (from 0 to 6) in which the (global) module 
  // ix is sitting (could return 7 if ix=isec)
  return TMath::Abs(ixGlo)-Int_t(TMath::Abs(ixGlo)/10)*10-1;
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentationV2::MorePads()
{
  AliFatal("Not implemented");
  return 0;
}

//_____________________________________________________________________________
void AliMUONTriggerSegmentationV2::Neighbours(Int_t /*iX*/, Int_t /*iY*/, 
                                              Int_t* /*Nlist*/, 
                                              Int_t /*Xlist*/[10], 
                                              Int_t /*Ylist*/[10]) 
{
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
AliMUONTriggerSegmentationV2::NextPad()
{
  AliFatal("Not implemented");
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentationV2::Npx() const
{
  return 124;// FIXME: this should not have to be done, if only we'd stick 
  // to a local (ix,iy) convention !!! 
  // return fSlatSegmentation->MaxPadIndexX()+1;
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSegmentationV2::Npy() const
{
  return 64;
//  return fSlatSegmentation->MaxPadIndexY()+1;
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::PC2LA(Int_t ixPC, Int_t iyPC,
                                    Int_t& ixLA, Int_t& iyLA)
{
  //
  // From PC to LA conventions for integers indices.
  //
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
AliMUONTriggerSegmentationV2::Print(Option_t* opt) const
{
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
AliMUONTriggerSegmentationV2::Sector(Int_t ix, Int_t iy)
{
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
AliMUONTriggerSegmentationV2::Sector(Float_t, Float_t)
{
  AliFatal("Not implemented");
  return 0;
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::SetCorrFunc(Int_t,TF1*)
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::SetDAnod(float)
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::SetHit(Float_t x, Float_t y)
{
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
AliMUONTriggerSegmentationV2::SetHit(Float_t x, Float_t y, Float_t)
{
  SetHit(x,y);
}

//_____________________________________________________________________________
void
AliMUONTriggerSegmentationV2::SetPad(Int_t ix, Int_t iy)
{
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
AliMUONTriggerSegmentationV2::SetPadSize(float,float)
{
  AliFatal("Not Implemented");
}

//_____________________________________________________________________________
Int_t 
AliMUONTriggerSegmentationV2::SigGenCond(Float_t,Float_t,Float_t)
{
  AliFatal("Not Implemented");
  return 0;
}

//_____________________________________________________________________________
void 
AliMUONTriggerSegmentationV2::SigGenInit(Float_t,Float_t,Float_t)
{
  AliFatal("Not Implemented");
}





