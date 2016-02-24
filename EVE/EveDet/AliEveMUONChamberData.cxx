// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel & Bogdan Vulpescu: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliMpVSegmentation.h"
#include "AliEveMUONChamberData.h"

#include <AliMUONGeometryTransformer.h>
#include <AliMpDEIterator.h>
#include <AliMpSector.h>
#include <AliMpPad.h>
#include <AliMpSegmentation.h>
#include <AliMpVSegmentation.h>

#include <TVector2.h>

#include <AliEveEventManager.h>

///////////////////////////////////////////////////////////////////////////////
///
/// AliEveMUONChamberData: geometry and digits
///
///////////////////////////////////////////////////////////////////////////////


ClassImp(AliEveMUONChamberData)

AliMUONGeometryTransformer* AliEveMUONChamberData::fgTransformer = 0;

//______________________________________________________________________________
AliEveMUONChamberData::AliEveMUONChamberData(Int_t chamber) :
  TObject(),
  fChamberID(0),
  fNDetElem(0),
  fNDigits(0),
  fNClusters(0),
  fNHits(0)
{
  //
  // constructor
  //

  fChamberID = chamber;
  fNDetElem  = 0;
  fNDigits   = 0;
  fNClusters = 0;
  fNHits     = 0;

  for (Int_t i = 0; i < 26; i++) {
    for (Int_t j = 0; j < 4; j++) {
      fFrameCoord[i][j] = 0.0;
    }
  }
  for (Int_t i = 0; i < 7*4096; i++) {
    fDigitBuffer[i] = 0.0;
  }
  for (Int_t i = 0; i < 5*256; i++) {
    fClusterBuffer[i] = 0.0;
  }
  for (Int_t i = 0; i < 3*256; i++) {
    fHitBuffer[i] = 0.0;
  }

  for (Int_t i = 0; i < 3; i++) {
    fChamberBox[i*2  ] = +9999;
    fChamberBox[i*2+1] = -9999;
  }

  if (fgTransformer == 0) {
    AliEveEventManager::Instance()->AssertGeometry();
    fgTransformer = new AliMUONGeometryTransformer();
    fgTransformer->LoadGeometryData();
  }

  Init(chamber);
}

//______________________________________________________________________________
AliEveMUONChamberData::~AliEveMUONChamberData()
{
  //
  // destructor
  //

}

//______________________________________________________________________________
void AliEveMUONChamberData::DropData()
{
  //
  // release the chamber data
  //

  fNDigits   = 0;
  fNClusters = 0;
  fNHits     = 0;

  return;

}

//______________________________________________________________________________
void AliEveMUONChamberData::Init(Int_t chamber)
{
  //
  // initialize the drawing coordinates of the chamber
  //

  Float_t locP[3], gloP[3], locD[3], gloD[3];
  Float_t deltax, deltay;
  AliMpDEIterator it;
  const AliMpVSegmentation *vseg;
  const AliMpSector *sector;
  TVector2 position;
  TVector2 dimension;

  for ( it.First(chamber); ! it.IsDone(); it.Next() ) {

    Int_t detElemId = it.CurrentDEId();

    if (chamber < 4) {

      sector = AliMpSegmentation::Instance()->GetSector(detElemId,AliMp::kCath0);

      position  = TVector2(sector->GetPositionX(), sector->GetPositionY());
      dimension = TVector2(sector->GetDimensionX(), sector->GetDimensionY()); // half length

      locP[0] =  position.Px();
      locP[1] =  position.Py();
      locD[0] =  dimension.Px() * 2.;
      locD[1] =  dimension.Py() * 2.;

      locP[2] = 0.0;
      locD[2] = 0.0;

      fgTransformer->Local2Global(detElemId,
				  locP[0], locP[1], locP[2],
				  gloP[0], gloP[1], gloP[2]);

      fgTransformer->Local2Global(detElemId,
				  locD[0], locD[1], locD[2],
				  gloD[0], gloD[1], gloD[2]);

      fFrameCoord[fNDetElem][0] = gloP[0];
      fFrameCoord[fNDetElem][1] = gloP[1];
      fFrameCoord[fNDetElem][2] = gloD[0];
      fFrameCoord[fNDetElem][3] = gloD[1];
      fFrameCoord[fNDetElem][4] = gloP[2]; // Z position

      fChamberBox[0] = TMath::Min(fChamberBox[0],gloP[0]-gloD[0]);
      fChamberBox[1] = TMath::Max(fChamberBox[1],gloP[0]+gloD[0]);
      fChamberBox[2] = TMath::Min(fChamberBox[2],gloP[1]-gloD[1]);
      fChamberBox[3] = TMath::Max(fChamberBox[3],gloP[1]+gloD[1]);
      fChamberBox[4] = TMath::Min(fChamberBox[4],gloP[2]);
      fChamberBox[5] = TMath::Max(fChamberBox[5],gloP[2]);

    } else {

//      if (!fgSegmentation->HasDE(detElemId)) {
//	printf("Segmentation has no %d detElemId! \n",detElemId);
//	continue;
//      }

      vseg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::kCath0);

      if (vseg == 0) {
	printf("No MpVSegmentation for %d detElemId! \n",detElemId);
	continue;
      }

      deltax = vseg->GetDimensionX();
      deltay = vseg->GetDimensionY();
      locP[0] =  -deltax;
      locP[1] =  -deltay;
      locD[0] =  +deltax;
      locD[1] =  +deltay;

      locP[2] = 0.0;
      locD[2] = 0.0;

      fgTransformer->Local2Global(detElemId,
				  locP[0], locP[1], locP[2],
				  gloP[0], gloP[1], gloP[2]);

      fgTransformer->Local2Global(detElemId,
				  locD[0], locD[1], locD[2],
				  gloD[0], gloD[1], gloD[2]);

      fFrameCoord[fNDetElem][0] = gloP[0];
      fFrameCoord[fNDetElem][1] = gloP[1];
      fFrameCoord[fNDetElem][2] = gloD[0];
      fFrameCoord[fNDetElem][3] = gloD[1];
      fFrameCoord[fNDetElem][4] = gloP[2]; // Z position

      fChamberBox[0] = TMath::Min(fChamberBox[0],gloP[0]);
      fChamberBox[0] = TMath::Min(fChamberBox[0],gloD[0]);
      fChamberBox[1] = TMath::Max(fChamberBox[1],gloP[0]);
      fChamberBox[1] = TMath::Max(fChamberBox[1],gloD[0]);
      fChamberBox[2] = TMath::Min(fChamberBox[0],gloP[1]);
      fChamberBox[2] = TMath::Min(fChamberBox[0],gloD[1]);
      fChamberBox[3] = TMath::Max(fChamberBox[1],gloP[1]);
      fChamberBox[3] = TMath::Max(fChamberBox[1],gloD[1]);
      fChamberBox[4] = TMath::Min(fChamberBox[4],gloP[2]);
      fChamberBox[5] = TMath::Max(fChamberBox[5],gloP[2]);

    }

    fNDetElem++;

  }  // end detElemId loop

}

//______________________________________________________________________________
void AliEveMUONChamberData::RegisterDigit(Int_t detElemId, Int_t cathode, Int_t ix, Int_t iy, Int_t charge)
{
  //
  // add a digit to this chamber
  //

  if ((fNDigits/7) == (4096-1)) return;

  Float_t locP[3], gloP[3], locD[3], gloD[3];

  const AliMpVSegmentation* vseg = AliMpSegmentation::Instance()
    ->GetMpSegmentation(detElemId,AliMp::GetCathodType(cathode));

  AliMpPad pad = vseg->PadByIndices(ix,iy,kTRUE);

  locP[0] = pad.GetPositionX();
  locP[1] = pad.GetPositionY();
  locD[0] = pad.GetDimensionX();
  locD[1] = pad.GetDimensionY();

  locP[2] = 0.0;
  locD[2] = 0.0;

  fgTransformer->Local2Global(detElemId,
			      locP[0], locP[1], locP[2],
			      gloP[0], gloP[1], gloP[2]);

  gloD[0] = locD[0];
  gloD[1] = locD[1];
  gloD[2] = gloP[2];

  if (cathode == 0) gloP[2] += 0.1;
  if (cathode == 1) gloP[2] -= 0.1;

  fDigitBuffer[fNDigits  ] = gloP[0];
  fDigitBuffer[fNDigits+1] = gloP[1];
  fDigitBuffer[fNDigits+2] = gloD[0];
  fDigitBuffer[fNDigits+3] = gloD[1];
  fDigitBuffer[fNDigits+4] = gloP[2];
  fDigitBuffer[fNDigits+5] = charge;
  fDigitBuffer[fNDigits+6] = cathode;

  fNDigits += 7;

}

//______________________________________________________________________________
void AliEveMUONChamberData::RegisterCluster(Int_t /*detElemId*/, Int_t cathode, Float_t clsX, Float_t clsY, Float_t clsZ, Float_t charge)
{
  //
  // add a reconstructed point (cluster) to this chamber
  //
  // identical clusters are registered for both cathode planes ...
  //

  if ((fNClusters/5) == (256-1)) return;

  fClusterBuffer[fNClusters  ] = clsX;
  fClusterBuffer[fNClusters+1] = clsY;
  fClusterBuffer[fNClusters+2] = clsZ;
  fClusterBuffer[fNClusters+3] = charge;
  fClusterBuffer[fNClusters+4] = cathode;

  fNClusters += 5;

}

//______________________________________________________________________________
void AliEveMUONChamberData::RegisterHit(Int_t /*detElemId*/, Float_t hitX, Float_t hitY, Float_t hitZ)
{
  //
  // add a simulation hit to this chamber
  //

  if ((fNHits/3) == (256-1)) return;

  fHitBuffer[fNHits  ] = hitX;
  fHitBuffer[fNHits+1] = hitY;
  fHitBuffer[fNHits+2] = hitZ;

  fNHits += 3;

}
