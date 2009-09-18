// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveITSModule.h"

#include <AliITSgeomTGeo.h>
#include <AliITSsegmentationSPD.h>
#include <AliITSsegmentationSDD.h>
#include <AliITSsegmentationSSD.h>

#include <AliITSdigitSPD.h>
#include <AliITSdigitSDD.h>
#include <AliITSdigitSSD.h>

#include <TEveTrans.h>
#include <TClonesArray.h>
#include <TStyle.h>


//______________________________________________________________________________
//
// Visualization of an ITS module.

ClassImp(AliEveITSModule)

Bool_t AliEveITSModule::fgStaticInitDone = kFALSE;

TEveFrameBox*    AliEveITSModule::fgSPDFrameBox = 0;
TEveFrameBox*    AliEveITSModule::fgSDDFrameBox = 0;
TEveFrameBox*    AliEveITSModule::fgSSDFrameBox = 0;

TEveRGBAPalette* AliEveITSModule::fgSPDPalette  = 0;
TEveRGBAPalette* AliEveITSModule::fgSDDPalette  = 0;
TEveRGBAPalette* AliEveITSModule::fgSSDPalette  = 0;

/******************************************************************************/

AliEveITSModule::AliEveITSModule(const Text_t* n, const Text_t* t) :
  TEveQuadSet(n, t),
  fInfo(0),
  fID(-1), fDetID(-1),
  fLayer(-1), fLadder(-1), fDet(-1),
  fDx(0), fDz(0), fDy(0)
{
  // Constructor.
}

AliEveITSModule::AliEveITSModule(Int_t gid, AliEveITSDigitsInfo* info) :
  TEveQuadSet(Form("ITS module %d", gid)),
  fInfo  (0),
  fID(-1), fDetID(-1),
  fLayer(-1), fLadder(-1), fDet(-1),
  fDx(0), fDz(0), fDy(0)
{
  // Constructor with module id and data-source.

  SetDigitsInfo(info);
  SetID(gid);
}

AliEveITSModule::~AliEveITSModule()
{
  // Destructor.

  if (fInfo) fInfo->DecRefCount();
}

/******************************************************************************/

void AliEveITSModule::InitStatics(AliEveITSDigitsInfo* info)
{
  // Initialize static variables.
  //
  // Warning all sensor sizes are in microns, here we transform them
  // to cm. In Eve half-lengths/widths are used, hence another 1/2.

  if (fgStaticInitDone) return;
  fgStaticInitDone = kTRUE;

  {
    Float_t dx = info->fSegSPD->Dx()*0.00005;
    Float_t dz = info->fSegSPD->Dz()*0.00005;

    fgSPDFrameBox = new TEveFrameBox();
    fgSPDFrameBox->SetAAQuadXZ(-dx, 0, -dz, 2*dx, 2*dz);
    fgSPDFrameBox->SetFrameColor(Color_t(31));
    fgSPDFrameBox->SetFrameFill(kTRUE);
    fgSPDFrameBox->IncRefCount();
    fgSPDPalette  = new TEveRGBAPalette(info->fSPDMinVal,info->fSPDMaxVal);
    fgSPDPalette->IncRefCount();
  }

  {
    Float_t dx = info->fSegSDD->Dx()*0.0001;
    Float_t dz = info->fSegSDD->Dz()*0.00005;

    fgSDDFrameBox = new TEveFrameBox();
    fgSDDFrameBox->SetAAQuadXZ(-dx, 0, -dz, 2*dx, 2*dz);
    fgSDDFrameBox->SetFrameColor(Color_t(32));
    fgSDDFrameBox->SetFrameFill(kTRUE);
    fgSDDFrameBox->IncRefCount();
    fgSDDPalette  = new TEveRGBAPalette(info->fSDDMinVal,info->fSDDMaxVal);
    fgSDDPalette->SetLimits(0, info->fSDDHighLim); // Set proper ADC range.
    fgSDDPalette->IncRefCount();
  }

  {
    Float_t dx = info->fSegSSD->Dx()*0.00005;
    Float_t dz = info->fSegSSD->Dz()*0.00005;

    fgSSDFrameBox = new TEveFrameBox();
    fgSSDFrameBox->SetAAQuadXZ(-dx, 0, -dz, 2*dx, 2*dz);
    fgSSDFrameBox->SetFrameColor(Color_t(33));
    fgSSDFrameBox->SetFrameFill(kTRUE);
    fgSSDFrameBox->IncRefCount();
    fgSSDPalette  = new TEveRGBAPalette(info->fSSDMinVal,info->fSSDMaxVal);
    fgSSDPalette->SetLimits(0, info->fSSDHighLim); // Set proper ADC range.
    fgSSDPalette->IncRefCount();
  }

}

/******************************************************************************/

void AliEveITSModule::SetDigitsInfo(AliEveITSDigitsInfo* info)
{
  // Set data and geometry source.

  if (fInfo == info) return;
  if (fInfo) fInfo->DecRefCount();
  fInfo = info;
  if (fInfo) fInfo->IncRefCount();
}

/******************************************************************************/

void AliEveITSModule::SetID(Int_t gid, Bool_t trans)
{
  // Set detector id.

  static const TEveException kEH("AliEveITSModule::SetID ");

  if (fInfo == 0)
    throw(kEH + "AliEveITSDigitsInfo not set.");

  Int_t firstSPD = AliITSgeomTGeo::GetModuleIndex(1,1,1);
  Int_t lastSSD  = AliITSgeomTGeo::GetNModules() - 1;
  if (gid < firstSPD || gid > lastSSD)
  {
    throw(kEH + Form("%d is not valid. ID range from %d to %d", gid,
		    firstSPD, lastSSD ));
  }

  fID = gid;

  if (!fgStaticInitDone)
  {
    InitStatics(fInfo);
  }

  AliITSgeomTGeo::GetModuleId(fID, fLayer, fLadder, fDet);
  TString strLadder = "Ladder";
  TString strSensor = "Sensor";
  TString symname;
  Int_t   id, nsector, nstave, nladder, rest;

  if (fID <= (AliITSgeomTGeo::GetModuleIndex(3,1,1) - 1))
  {
    // SPD

    SetFrame(fgSPDFrameBox);
    SetPalette(fgSPDPalette);

    symname += strLadder;
    if (fID < 80)
    {
      nsector = fID/8;
      rest    = fID - 8*nsector;
      nstave  = 1;
    }
    else
    {
      id      = fID - 80;
      nsector = id/8;
      rest    = id - 8*nsector;
      nstave  = 1;
    }
    if (rest < 4) nstave = 0;
    rest    -= 4*nstave;
    symname += rest;
    SetName(symname);
    fDetID = 0;
    fDx = fInfo->fSegSPD->Dx()*0.00005;
    fDz = 3.50;
    fDy = fInfo->fSegSPD->Dy()*0.00005;
  }
  else if (fID <= (AliITSgeomTGeo::GetModuleIndex(5,1,1) - 1))
  {
    // SDD

    SetFrame(fgSDDFrameBox);
    SetPalette(fgSDDPalette);

    symname += strSensor;
    if (fID < 324)
    {
      id      = fID - 240;
      nladder = id/6;
      rest    = id - 6*nladder;
    }
    else
    {
      id      = fID - 324;
      nladder = id/8;
      rest    = id - 8*nladder;
    }
    symname += rest;
    SetName(symname);
    fDetID = 1;
    fDx = fInfo->fSegSDD->Dx()*0.0001;
    fDz = fInfo->fSegSDD->Dz()*0.00005;
    fDy = fInfo->fSegSDD->Dy()*0.00005;
  }
  else
  {
    // SSD

    SetFrame(fgSSDFrameBox);
    SetPalette(fgSSDPalette);

    symname += strSensor;
    if (fID < 1248)
    {
      id      = fID - 500;
      nladder = id/22;
      rest    = id - nladder*22;
    }
    else
    {
      id      = fID - 1248;
      nladder = id/25;
      rest    = id - nladder*25;
    }
    symname += rest;
    SetName(symname);
    fDetID = 2;
    fInfo->fSegSSD->SetLayer(fLayer);
    fDx = fInfo->fSegSSD->Dx()*0.00005;
    fDz = fInfo->fSegSSD->Dz()*0.00005;
    fDy = fInfo->fSegSSD->Dy()*0.00005;
  }

  LoadQuads();
  ComputeBBox();
  InitMainTrans();
  if (trans)
    SetTrans();
}

void AliEveITSModule::LoadQuads()
{
  // Read module data from source and create low-level objects for
  // visualization - called quads.

  TClonesArray *digits  = fInfo->GetDigits(fID, fDetID);
  if (!digits) return;

  Int_t         ndigits = digits ? digits->GetEntriesFast() : 0;

  Float_t       x, z, dpx, dpz;
  Int_t         i, j;

  switch(fDetID)
  {
    case 0:
    {
      AliITSsegmentationSPD* seg =  fInfo->fSegSPD;

      Reset(kQT_RectangleXZFixedY, kFALSE, 32);

      for (Int_t k=0; k<ndigits; ++k)
      {
	AliITSdigit *d = (AliITSdigit*) digits->UncheckedAt(k);
	j = d->GetCoord1();
	i = d->GetCoord2();
	x  = -0.5*seg->Dx() + i*seg->Dpx(0);
	x *=  0.0001;
	fInfo->GetSPDLocalZ(j, z);
	dpx = seg->Dpx(i)*0.0001;
	dpz = seg->Dpz(j)*0.0001;

	AddQuad(x, z, dpx, dpz);
	QuadValue(1); // In principle could have color based on number of neigbours
	QuadId(d);
      }
      break;
    }

    case 1:
    {
      AliITSsegmentationSDD *seg =  fInfo->fSegSDD;

      Reset(kQT_RectangleXZFixedY, kFALSE, 32);

      for (Int_t k=0; k<ndigits; ++k)
      {
	AliITSdigit* d = (AliITSdigit*) digits->UncheckedAt(k);

	j = d->GetCoord1();
	i = d->GetCoord2();
	seg->DetToLocal(i, j, x, z);
	dpx = seg->Dpx(i)*0.0001;
	dpz = seg->Dpz(j)*0.0001;

	AddQuad(x-2*dpx, z - dpz*0.5, 4*dpx, dpz);
	QuadValue(d->GetSignal());
	QuadId(d);
      }
      break;
    }

    case 2:
    {
      AliITSsegmentationSSD* seg = fInfo->fSegSSD;

      Reset(kQT_LineXZFixedY, kFALSE, 32);

      Float_t ap, an; // positive/negative angles -> offsets
      seg->Angles(ap, an);
      ap =   TMath::Tan(ap) * fDz;
      an = - TMath::Tan(an) * fDz;

      for (Int_t k=0; k<ndigits; ++k)
      {
	AliITSdigit *d = (AliITSdigit*) digits->UncheckedAt(k);

	j = d->GetCoord1();
	i = d->GetCoord2();
	// !!!! The following function complains about not being verified.
	// !!!! Experts should check.
	seg->DetToLocal(i,j,x,z);

	Float_t a = (d->GetCoord1() == 1) ? ap : an;

	AddLine(x-a, -fDz, 2*a, 2*fDz);
	QuadValue(d->GetSignal());
	QuadId(d);
	// printf("%3d -> %3d -> %8x\n", d->GetSignal(), ci, fQuads.back().color);
      }
      break;
    }

  } // end switch

  RefitPlex();
}

/******************************************************************************/

void AliEveITSModule::SetTrans()
{
  // Set transformation matrix based on module id (use geometry to
  // retrieve this information).

  fMainTrans->SetFrom(*AliITSgeomTGeo::GetMatrix(fID));
}

/******************************************************************************/

void AliEveITSModule::DigitSelected(Int_t idx)
{
  // Override secondary select (alt-click) from TEveQuadSet.

  DigitBase_t *qb  = GetDigit(idx);
  TObject     *obj = qb->fId.GetObject();
  AliITSdigit *d   = dynamic_cast<AliITSdigit*>(obj);
  printf("AliEveITSModule::QuadSelected "); Print();
  printf("  idx=%d, value=%d, obj=0x%lx, digit=0x%lx\n",
	 idx, qb->fValue, (ULong_t)obj, (ULong_t)d);
  if (d)
    printf("  coord1=%3d coord2=%3d signal=%d\n",
	   d->GetCoord1(), d->GetCoord2(), d->GetSignal());

}

/******************************************************************************/

void AliEveITSModule::Print(Option_t* ) const
{
  // Print object summary information.

  printf("AliEveITSModule: ID %d, layer %d, ladder %d, det %d\n",
         fID, fLayer, fLadder, fDetID);
}
