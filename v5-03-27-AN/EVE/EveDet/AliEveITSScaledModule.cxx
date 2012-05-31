// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveITSScaledModule.h"

#include <AliITSsegmentationSPD.h>
#include <AliITSsegmentationSDD.h>
#include <AliITSsegmentationSSD.h>

#include <AliITSdigitSPD.h>
#include <AliITSdigitSDD.h>
#include <AliITSdigitSSD.h>

#include <TMath.h>
#include <TClonesArray.h>

//==============================================================================
//==============================================================================
// AliEveDigitScaleInfo
//==============================================================================

//______________________________________________________________________________
//
// Encapsulates current state of scaling and agglomeration.

ClassImp(AliEveDigitScaleInfo)

AliEveDigitScaleInfo::AliEveDigitScaleInfo():
  fScale(1),
  fStatType (kSTAverage),
  fSyncPalette(kFALSE)
{
}

void AliEveDigitScaleInfo::ScaleChanged(Int_t s)
{
  fScale = s;

  AliEveITSScaledModule* sm;
  RefMap_i i = fBackRefs.begin();
  // #endif
  while (i != fBackRefs.end())
  {
    sm = dynamic_cast<AliEveITSScaledModule*>(i->first);
    // #endif
    if(sm) sm->LoadQuads();
    ++i;
  }
}

void AliEveDigitScaleInfo::StatTypeChanged(Int_t t)
{
  fStatType = t;
  fSyncPalette = kTRUE;

  AliEveITSScaledModule* sm;
  RefMap_i i = fBackRefs.begin();
  // #endif
  while (i != fBackRefs.end())
  {
    sm = dynamic_cast<AliEveITSScaledModule*>(i->first);
    // #endif
    if (sm) sm->SetQuadValues();
    ++i;
  }
}

//______________________________________________________________________________
// ScaledDigit_t
//

AliEveITSScaledModule::ScaledDigit_t::ScaledDigit_t() :
  TObject(),
  fN(0),
  fSum(0), fSqrSum(0),
  fMinI(-1), fMinJ(-1), fMaxI(-1), fMaxJ(-1)
{
}

AliEveITSScaledModule::ScaledDigit_t::ScaledDigit_t(Int_t di, Int_t dj) :
  TObject(),
  fN(0),
  fSum(0), fSqrSum(0),
  fMinI(di), fMinJ(dj), fMaxI(di), fMaxJ(dj)
{
}

void AliEveITSScaledModule::ScaledDigit_t::Dump() const
{
  printf("N %d, sum %f, sqr_sum %f", fN, fSum, fSqrSum);
}


//==============================================================================
//==============================================================================
// AliEveITSScaledModule
//==============================================================================

//______________________________________________________________________________
//
// Visualization of an ITS module with digits aggregated
// on a grid of pre-defined size.

ClassImp(AliEveITSScaledModule)

AliEveITSScaledModule::AliEveITSScaledModule(Int_t gid, AliEveITSDigitsInfo* info, AliEveDigitScaleInfo* si):
  AliEveITSModule("AliEveITSScaledModule", "AliEveITSScaledModule"),
  fNx(-1),
  fNz(-1),
  fNCx(-1),
  fNCz(-1),
  fScaleInfo(si),
  fDigitsMap()
{
  SetOwnIds(kTRUE);

  SetDigitsInfo(info);
  SetID(gid);
  fScaleInfo->IncRefCount(this);
}

AliEveITSScaledModule::~AliEveITSScaledModule()
{
  fScaleInfo->DecRefCount(this);
}

/******************************************************************************/

void AliEveITSScaledModule::LoadQuads()
{
  // Here we still use 'z' for the name of axial coordinates.
  // The transforamtion matrix aplied rotates y -> z.
  // We need this as TEveQuadSet offers optimized treatment for
  // quads in the x-y plane.

  TClonesArray *digits = fInfo->GetDigits(fID, fDetID);
  if (!digits) return;

  Int_t ndigits = digits->GetEntriesFast();

  Float_t       x, z, zo, dpx, dpz; // orig cells size, pos
  Int_t         i, j;               // orig cells idx
  Int_t         c1, c2;             // original coordinates

  Int_t id;
  std::map<Int_t, Int_t> dmap;
  std::map<Int_t, Int_t>::iterator miter;

  Int_t scale = fScaleInfo->GetScale() -1;
  switch(fDetID)
  {
    case 0:
    {
      // SPD
      Reset(kQT_RectangleXZFixedY, kFALSE, 32);

      fNCz = fInfo->fSPDScaleZ[scale];
      fNCx = fInfo->fSPDScaleX[scale];
      fNz  = Int_t(fInfo->fSegSPD->Npz()/fNCz);
      fNx  = Int_t(fInfo->fSegSPD->Npx()/fNCx);
      dpz = 2*fDz/fNz;
      dpx = 2*fDx/fNx;
      //printf("SPD orig cells (%d, %d) (%d, %d)\n", fInfo->fSegSPD->Npx(), fInfo->fSegSPD->Npz(), Nx, Nz);

      AliITSdigitSPD *od ;
      for (Int_t k=0; k<ndigits; ++k)
      {
	od = (AliITSdigitSPD*) digits->UncheckedAt(k);

	fInfo->GetSPDLocalZ(od->GetCoord1(),zo);
        c1 = od->GetCoord1(); c2 = od->GetCoord2();
	i  = Int_t((zo+fDz)/dpz);
	j  = Int_t((od->GetCoord2()*fNx)/fInfo->fSegSPD->Npx());
	id = j*fNx + i;

        ScaledDigit_t* sd = 0;
        miter = dmap.find(id);
	if(miter == dmap.end())
	{
          dmap[id] = fPlex.Size();
          z = dpz*(i) - fDz;
          x = dpx*(j) - fDx;
          AddQuad(x, z, dpx, dpz);
          sd = new ScaledDigit_t(c1, c2);
	  QuadId(sd);
	}
        else
	{
	  sd = static_cast<ScaledDigit_t*>(GetId(miter->second));
          if(c1 < sd->fMinI)
	    sd->fMinI = c1;
	  else if( c1 > sd->fMaxI)
            sd->fMaxI = c1;

          if(c2 < sd->fMinJ)
	    sd->fMinJ = c2;
	  else if( c2 > sd->fMaxJ)
	    sd->fMaxJ = c2;
	}

	sd->fN++;
	sd->fSum  += od->GetSignal();
	sd->fSqrSum += od->GetSignal()*od->GetSignal();
      }
      break;
    }
    case 1:
    {
      // SDD
      Reset(kQT_RectangleXZFixedY, kFALSE, 32);

      fNCz = fInfo->fSDDScaleZ[scale];
      fNCx = fInfo->fSDDScaleX[scale];
      fNz  = Int_t(fInfo->fSegSDD->Npz()/fNCz);
      fNx  = Int_t(fInfo->fSegSDD->Npx()/fNCx);
      dpz  = 2*fDz/fNz;
      dpx  = 2*fDx/fNx;

      AliITSdigitSDD *od = 0;
      for (Int_t k = 0; k < ndigits; ++k)
      {
	od = (AliITSdigitSDD*)digits->UncheckedAt(k);
	fInfo->fSegSDD->DetToLocal(od->GetCoord2(), od->GetCoord1(),x,z);
	z += fDz;
	x += fDx;
	i  = Int_t(z/dpz);
	j  = Int_t(x/dpx);
	//printf("Mod %d coord %d,%d out of %d,%d :: ORIG coord %d,%d out of %d,%d \n",fID,
	//       i,j,Nz,Nx,od->GetCoord1(),od->GetCoord2(),fInfo->fSegSDD->Npz(),fInfo->fSegSDD->Npx());

	id = j*fNx + i;
	c1 = od->GetCoord1(); c2 = od->GetCoord2();

	ScaledDigit_t* sd = 0;
	miter = dmap.find(id);
	if(miter == dmap.end())
	{
	  dmap[id] = fPlex.Size();
	  z = dpz*(i) - fDz;
	  x = dpx*(j) - fDx;
	  AddQuad(x, z, dpx, dpz);
	  sd = new ScaledDigit_t(od->GetCoord1(),od->GetCoord2());
	  QuadId(sd);
	}
	else
	{
	  sd = static_cast<ScaledDigit_t*>(GetId(miter->second));
	  if(c1 < sd->fMinI)
	    sd->fMinI = c1;
	  else if( c1 > sd->fMaxI)
	    sd->fMaxI = c1;

	  if(c2 < sd->fMinJ)
	    sd->fMinJ = c2;
	  else if( c2 > sd->fMaxJ)
	    sd->fMaxJ = c2;
	}
	sd->fN++;
	sd->fSum  += od->GetSignal();
	sd->fSqrSum += od->GetSignal()*od->GetSignal();
      }
      break;
    }
    case 2:
    {
      // SSD
      Reset(kQT_LineXZFixedY, kFALSE, 32);

      AliITSsegmentationSSD* seg = fInfo->fSegSSD;
      Float_t ap, an; // positive/negative angles -> offsets
      seg->Angles(ap, an);
      ap =   TMath::Tan(ap) * fDz;
      an = - TMath::Tan(an) * fDz;

      fNCx  = fInfo->fSSDScale[scale];
      fNz  = 1;
      fNx  = Int_t(fInfo->fSegSSD->Npx()/fNCx);
      dpz = 2*fDz/fNz;
      dpx = 2*fDx/fNx;

      AliITSdigitSSD *od = 0;
      for (Int_t k=0; k<ndigits; k++) {
	od=(AliITSdigitSSD*)digits->UncheckedAt(k);
	if(od->GetCoord1() == 1)
	  i = 1; // p side
	else
	  i= -1; // n side
	j = Int_t(od->GetCoord2()/fNCx);
	c1 = od->GetCoord1(); c2 = od->GetCoord2();
	id = j*i;

	ScaledDigit_t* sd = 0;
	miter = dmap.find(id);
	if(miter == dmap.end())
	{
	  // printf("orig digit %d,%d scaled %d,%d \n",od->GetCoord1(),od->GetCoord2(),i,j);
	  dmap[id] = fPlex.Size();
	  z = dpz*(i) - fDz;
	  x = dpx*(j) - fDx;
	  Float_t a = ( od->GetCoord1() == 1) ? ap : an;
	  AddLine(x-a, -fDz, 2*a, 2*fDz);

	  sd = new ScaledDigit_t(c1, c2);
	  QuadId(sd);
	}
	else
	{
	  sd = static_cast<ScaledDigit_t*>(GetId(miter->second));
	  if(c1 < sd->fMinI)
	    sd->fMinI = c1;
	  else if( c1 > sd->fMaxI)
	    sd->fMaxI = c1;

	  if(c2 < sd->fMinJ)
	    sd->fMinJ = c2;
	  else if( c2 > sd->fMaxJ)
	    sd->fMaxJ = c2;
	}
	sd->fN++;
	sd->fSum  += od->GetSignal();
	sd->fSqrSum += od->GetSignal()*od->GetSignal();
      } // for digits
      break;
    } // end case 2
  } // end switch

  SetQuadValues();
  RefitPlex();
}

/******************************************************************************/

void AliEveITSScaledModule::SetQuadValues()
{
  if(fScaleInfo->GetSyncPalette()) SyncPalette();

  Int_t num = fPlex.Size();
  for (Int_t i = 0 ; i < num; i++)
  {
    ScaledDigit_t* sd = static_cast<ScaledDigit_t*>(GetId(i));
    Int_t v = 0;
    switch(fScaleInfo->GetStatType())
    {
      using namespace TMath;

      case AliEveDigitScaleInfo::kSTOccup:
	v = Nint((100.0*sd->fN) / (fNCx*fNCz));
	break;
      case AliEveDigitScaleInfo::kSTAverage:
	v = Nint((Double_t) sd->fSum / sd->fN);
	break;
      case AliEveDigitScaleInfo::kSTRms:
	v = Nint(Sqrt(sd->fSqrSum) / sd->fN);
	break;
    }
    DigitBase_t* qb = GetDigit(i);
    qb->fValue = v;
  }
}

/******************************************************************************/

void AliEveITSScaledModule::SyncPalette()
{
  // printf("AliEveITSScaledModule::SyncPalette()\n");
  if(fScaleInfo->GetStatType() == AliEveDigitScaleInfo::kSTOccup)
  {
    // SPD
    AliEveITSModule::fgSPDPalette->SetLimits(0, 100);
    AliEveITSModule::fgSPDPalette->SetMinMax(0, 100);

    // SDD
    AliEveITSModule::fgSDDPalette->SetLimits(0, 100);
    AliEveITSModule::fgSDDPalette->SetMinMax(0, 100);

    // SSD
    AliEveITSModule::fgSSDPalette->SetLimits(0, 100);
    AliEveITSModule::fgSDDPalette->SetMinMax(0, 100);
  }
  else
  {
    AliEveITSDigitsInfo& di = *fInfo;
    // SPD
    AliEveITSModule::fgSPDPalette->SetLimits(0, di.fSPDHighLim);
    AliEveITSModule::fgSPDPalette->SetMinMax(di.fSPDMinVal, di.fSPDMaxVal);

    // SDD
    AliEveITSModule::fgSDDPalette->SetLimits(0, di.fSDDHighLim);
    AliEveITSModule::fgSDDPalette->SetMinMax(di.fSDDMinVal, di.fSDDMaxVal);

    // SSD
    AliEveITSModule::fgSSDPalette->SetLimits(0, di.fSSDHighLim);
    AliEveITSModule::fgSSDPalette->SetMinMax(di.fSSDMinVal, di.fSSDMaxVal);
  }

  fScaleInfo->SetSyncPalette(kFALSE);
}

/******************************************************************************/

void AliEveITSScaledModule::GetScaleData(Int_t& cnx, Int_t& cnz, Int_t& total) const
{
  cnx   = fNx;
  cnz   = fNz;
  total = cnx*cnz;
}

/******************************************************************************/

void  AliEveITSScaledModule::DigitSelected(Int_t idx)
{
  // Override control-click from TEveQuadSet
  printf("AliEveITSScaledModule::DigitSelected "); Print();

  // DigitBase_t *qb  = GetDigit(idx);
  TObject     *obj = GetId(idx);
  ScaledDigit_t *sd = static_cast<ScaledDigit_t*>(obj);
  TClonesArray *digits = fInfo->GetDigits(fID, fDetID);
  Int_t ndigits = digits->GetEntriesFast();

  printf("%d digits in cell scaleX = %d,  scaleZ = %d \n", sd->fN, fNCx, fNCz);

  Int_t il = 0;
  for(Int_t k=0; k<ndigits; k++)
  {
    AliITSdigit *d = (AliITSdigit*) digits->UncheckedAt(k);

    if(d->GetCoord1()>=sd->fMinI && d->GetCoord1()<=sd->fMaxI &&
       d->GetCoord2()>=sd->fMinJ && d->GetCoord2()<=sd->fMaxJ)
    {
      printf("%3d, %3d: %3d", d->GetCoord1(), d->GetCoord2(), d->GetSignal());
      printf(" | ");
      il++;
      if(il>5) {
	printf("\n");
	il = 0;
      }
    }
  }
  if(il) printf("\n");
}
