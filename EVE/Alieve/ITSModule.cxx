#include "ITSModule.h"

#include <AliITSdigitSPD.h>
#include <AliITSdigitSDD.h>
#include <AliITSdigitSSD.h>

#include <TStyle.h>

using namespace Reve;
using namespace Alieve;

/*
Short_t ITSModule::fgSDDThreshold  = 5;
Short_t ITSModule::fgSDDMaxVal     = 80;
Short_t ITSModule::fgSSDThreshold  = 2;
Short_t ITSModule::fgSSDMaxVal     = 100;
*/

Bool_t       ITSModule::fgStaticInitDone = kFALSE;

FrameBox*    ITSModule::fgSPDFrameBox = 0;
FrameBox*    ITSModule::fgSDDFrameBox = 0;
FrameBox*    ITSModule::fgSSDFrameBox = 0;

RGBAPalette* ITSModule::fgSPDPalette  = 0;
RGBAPalette* ITSModule::fgSDDPalette  = 0;
RGBAPalette* ITSModule::fgSSDPalette  = 0;

//__________________________________________________________________________
// ITSModule
//
//

ClassImp(ITSModule)

/**************************************************************************/

ITSModule::ITSModule(const Text_t* n, const Text_t* t) :
  QuadSet(n, t),
  fInfo(0),
  fID(-1), fDetID(-1),
  fLayer(-1), fLadder(-1), fDet(-1),
  fDx(0), fDz(0), fDy(0)
{}

ITSModule::ITSModule(Int_t id, ITSDigitsInfo* info) :
  QuadSet(Form("ITS module %d", id)),
  fInfo  (0),
  fID(-1), fDetID(-1),
  fLayer(-1), fLadder(-1), fDet(-1),
  fDx(0), fDz(0), fDy(0)
{
  SetDigitsInfo(info);
  SetID(id);
}

ITSModule::~ITSModule()
{
  if(fInfo) fInfo->DecRefCount();
}

/**************************************************************************/

void ITSModule::SetDigitsInfo(ITSDigitsInfo* info)
{
  if(fInfo) fInfo->DecRefCount();
  fInfo = info;
  if(fInfo) fInfo->IncRefCount();
}

void ITSModule::SetID(Int_t id)
{
  static const Exc_t eH("ITSModule::SetID ");

  if(fInfo == 0)
    throw(eH + "ITSDigitsInfo not set.");

  if (id < fInfo->fGeom->GetStartSPD() || id > fInfo->fGeom->GetLastSSD())
    throw(eH + Form("%d is not valid. ID range from %d to %d", id,
		     fInfo->fGeom->GetStartSPD(), fInfo->fGeom->GetLastSSD()));

  fID = id;
  InitModule();
}

/**************************************************************************/

void ITSModule::InitStatics(ITSDigitsInfo* info)
{
  if (fgStaticInitDone) return;

  {
    Float_t dx = info->fSegSPD->Dx()*0.00005;
    Float_t dz = 3.48; 

    fgSPDFrameBox = new FrameBox();
    fgSPDFrameBox->SetAAQuadXZ(-dx, 0, -dz, 2*dx, 2*dz);
    fgSPDFrameBox->SetFrameColor((Color_t) 31);
    fgSPDPalette  = new RGBAPalette(0, 1, kFALSE, kFALSE);
  }

  {
    Float_t dx = info->fSegSDD->Dx()*0.0001;
    Float_t dz = info->fSegSDD->Dz()*0.00005;

    fgSDDFrameBox = new FrameBox();
    fgSDDFrameBox->SetAAQuadXZ(-dx, 0, -dz, 2*dx, 2*dz);
    fgSDDFrameBox->SetFrameColor((Color_t) 32);
    fgSDDPalette  = new RGBAPalette(5, 80, kTRUE, kFALSE);
  }

  {
    Float_t dx = info->fSegSSD->Dx()*0.00005;
    Float_t dz = info->fSegSSD->Dz()*0.00005;

    fgSSDFrameBox = new FrameBox();
    fgSSDFrameBox->SetAAQuadXZ(-dx, 0, -dz, 2*dx, 2*dz);
    fgSSDFrameBox->SetFrameColor((Color_t) 33);
    fgSSDPalette  = new RGBAPalette(2, 100, kTRUE, kFALSE);
  }

}

/**************************************************************************/

void ITSModule::InitModule()
{
  if (!fgStaticInitDone)  InitStatics(fInfo);

  fInfo->fGeom->GetModuleId(fID,fLayer,fLadder,fDet);
  TString strLadder = "Ladder";
  TString strSensor = "Sensor";
  TString symname;
  Int_t id, nsector, nstave, nladder, rest;
  
  if (fID <= fInfo->fGeom->GetLastSPD())
  {
    // SPD

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
    fDz = 3.48; 
    fDy = fInfo->fSegSPD->Dy()*0.00005;

  }
  else if (fID <= fInfo->fGeom->GetLastSDD())
  {
    // SDD

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
  SetTrans();
}

void ITSModule::LoadQuads()
{
  // Here we still use 'z' for the name of axial coordinates.
  // The transforamtion matrix aplied rotates y -> z.
  // We need this as QuadSet offers optimized treatment for
  // quads in the x-y plane.

  // printf("its module load quads \n");

  TClonesArray *digits;
  Float_t       x, z, dpx, dpz; 
  Int_t         i, j, ndigits;
  digits  = fInfo->GetDigits(fID, fDetID);
  ndigits = digits->GetEntriesFast(); 

  switch(fDetID)
  {

    case 0: { // SPD
      AliITSsegmentationSPD* seg =  fInfo->fSegSPD; 

      Reset(QT_AxisAlignedFixedY, kFALSE, 32);
      fFrame   = fgSPDFrameBox;
      fPalette = fgSPDPalette;

      for (Int_t k=0; k<ndigits; ++k)
      {
	AliITSdigitSPD *d = (AliITSdigitSPD*) digits->UncheckedAt(k);
	j = d->GetCoord1();
	i = d->GetCoord2();
	x  = -0.5*seg->Dx() + i*seg->Dpx(0);
	x *=  0.0001;
	fInfo->GetSPDLocalZ(j, z);
	dpx = seg->Dpx(i)*0.0001;
	dpz = seg->Dpz(j)*0.0001;

	AddQuad(x, z, dpx, dpz);
	QuadValue(1); // In principle could have color based on number of neigbours
      }
      break;
    }

    case 1: { // SDD
      AliITSsegmentationSDD *seg =  fInfo->fSegSDD; 

      Reset(QT_AxisAlignedFixedY, kFALSE, 32);
      fFrame   = fgSDDFrameBox;
      fPalette = fgSDDPalette;

      for (Int_t k=0; k<ndigits; ++k)
      {
	AliITSdigitSDD* d = (AliITSdigitSDD*) digits->UncheckedAt(k);

	// if (d->GetSignal() > fgSDDThreshold)
	{
	  j = d->GetCoord1();
	  i = d->GetCoord2();
	  seg->DetToLocal(i, j, x, z);
	  dpx = seg->Dpx(i)*0.0001;
	  dpz = seg->Dpz(j)*0.0001;

	  AddQuad(x, z, dpx, dpz);
	  QuadValue(d->GetSignal());
	}
      }
      break;
    }

    case 2: { // SSD
      AliITSsegmentationSSD* seg = fInfo->fSegSSD; 

      Reset(QT_LineFixedY, kFALSE, 32);
      fFrame   = fgSSDFrameBox;
      fPalette = fgSSDPalette;

      Float_t ap, an; // positive/negative angles -> offsets
      seg->Angles(ap, an);
      ap =   TMath::Tan(ap) * fDz;
      an = - TMath::Tan(an) * fDz;

      for (Int_t k=0; k<ndigits; ++k)
      {
	AliITSdigitSSD *d = (AliITSdigitSSD*) digits->UncheckedAt(k);
	// if(d->GetSignal() > fgSSDThreshold)
	{
	  j = d->GetCoord1();
	  i = d->GetCoord2();
	  seg->DetToLocal(i,j,x,z);

	  Float_t a = ( d->GetCoord1() == 1) ? ap : an;

	  AddQuad(x-a, -fDz, 2*a, 2*fDz);
	  QuadValue(d->GetSignal());
	  // printf("%3d -> %3d -> %8x\n", d->GetSignal(), ci, fQuads.back().color);
	}
      }
      break;
    }

  } // end switch

  RefitPlex();
}

/**************************************************************************/

void ITSModule::SetTrans()
{
  Double_t x[9];
  fHMTrans.UnitTrans();

   // column major
  fInfo->fGeom->GetRotMatrix(fID, x);
  fHMTrans.SetBaseVec(1, x[0], x[3], x[6]);
  fHMTrans.SetBaseVec(2, x[1], x[4], x[7]);
  fHMTrans.SetBaseVec(3, x[2], x[5], x[8]);
  // translation
  fInfo->fGeom->GetTrans(fID, x);  
  fHMTrans.SetBaseVec(4, x);
}

/**************************************************************************/

void ITSModule::Print(Option_t* ) const
{
  printf("ID %d, layer %d, ladder %d, det %d \n", fID, fLayer, fLadder, fDetID);
}
