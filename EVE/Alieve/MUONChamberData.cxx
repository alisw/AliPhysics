#include "MUONChamberData.h"

#include <AliMUONSegmentation.h>
#include <AliMUONConstants.h>
#include <AliMUONGeometryTransformer.h>
#include <AliMUONSegFactory.h>
#include <mapping/AliMpDEIterator.h>
#include <mapping/AliMpSectorSegmentation.h>
#include <mapping/AliMpSector.h>
#include <mapping/AliMpPad.h>
#include <mapping/AliMpStationType.h>
#include <mapping/AliMpDEManager.h>
#include <mapping/AliMpSegmentation.h>

#include <TMath.h>
#include <TVector2.h>

///////////////////////////////////////////////////////////////////////////////
///
/// MUONChamberData: geometry and digits
///
///////////////////////////////////////////////////////////////////////////////

using namespace Reve;
using namespace Alieve;

ClassImp(MUONChamberData)

AliMUONSegmentation* MUONChamberData::fgSegmentation = 0;
AliMUONGeometryTransformer* MUONChamberData::fgTransformer = 0;

//______________________________________________________________________
MUONChamberData::MUONChamberData(Int_t chamber)
{
  //
  // constructor
  //

  fChamberID = chamber;
  fNDetElem  = 0;
  fNDigits   = 0;
  for (Int_t i = 0; i < 26; i++) {
    for (Int_t j = 0; j < 4; j++) {
      fFrameCoord[i][j] = 0.0;
    }
  }
  for (Int_t i = 0; i < 7*4096; i++) {
    fDigitBuffer[i] = 0.0;
  }

  for (Int_t i = 0; i < 3; i++) {
    fChamberBox[i*2  ] = +9999;
    fChamberBox[i*2+1] = -9999;
  }

  if (fgSegmentation == 0) {
    AliMUONSegFactory segFactory("volpath.dat","transform.dat");
    fgSegmentation = segFactory.CreateSegmentation("FactoryV4");
    fgTransformer = new AliMUONGeometryTransformer(true);
    fgTransformer->ReadGeometryData("volpaths.dat","transform.dat");
  }

  Init(chamber);

}

//______________________________________________________________________
MUONChamberData::~MUONChamberData()
{
  //
  // destructor
  //

}

//______________________________________________________________________
void MUONChamberData::DropData()
{
  //
  // release the chamber data
  //

  return;

}

//______________________________________________________________________
void MUONChamberData::Init(Int_t chamber)
{
  //
  // initialize the drawing coordinates of the chamber
  //

  Float_t locP[3], gloP[3], locD[3], gloD[3];
  Float_t deltax, deltay;
  AliMpDEIterator it;
  const AliMpVSegmentation *vseg;
  const AliMpSectorSegmentation *sseg;
  const AliMpSector *sector;
  TVector2 position;
  TVector2 dimension;

  for ( it.First(chamber); ! it.IsDone(); it.Next() ) {

    Int_t detElemId = it.CurrentDE();

    if (chamber < 4) {

      sseg = (AliMpSectorSegmentation*)AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,0);
      sector = sseg->GetSector();
      
      position  = sector->Position(); 
      dimension = sector->Dimensions(); // half length
      
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

      //AliMpStationType stationType = AliMpDEManager::GetStationType(detElemId);

      if (!fgSegmentation->HasDE(detElemId)) {
	printf("Segmentation has no %d detElemId! \n",detElemId);
	continue;
      }

      vseg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,0);

      if (vseg == 0) {
	printf("No MpVSegmentation for %d detElemId! \n",detElemId);
	continue;
      }

      deltax = vseg->Dimensions().X();
      deltay = vseg->Dimensions().Y();
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

  //printf("ChamberBox %d \n",chamber);
  //printf("%f %f \n",fChamberBox[0],fChamberBox[1]);
  //printf("%f %f \n",fChamberBox[2],fChamberBox[3]);
  //printf("%f %f \n",fChamberBox[4],fChamberBox[5]);

}

//______________________________________________________________________
void MUONChamberData::RegisterDigit(Int_t detElemId, Int_t cathode, Int_t ix, Int_t iy, Int_t charge)
{
  //
  // add a digit to this chamber
  //

  Float_t locP[3], gloP[3], locD[3], gloD[3];

  const AliMpVSegmentation* vseg = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,cathode);

  AliMpPad pad = vseg->PadByIndices(AliMpIntPair(ix,iy),kTRUE);
  
  locP[0] = pad.Position().X();
  locP[1] = pad.Position().Y();
  locD[0] = pad.Dimensions().X();
  locD[1] = pad.Dimensions().Y();
  
  locP[2] = 0.0;
  locD[2] = 0.0;

  fgTransformer->Local2Global(detElemId, 
			      locP[0], locP[1], locP[2], 
			      gloP[0], gloP[1], gloP[2]);
  
  gloD[0] = locD[0];
  gloD[1] = locD[1];
  gloD[2] = gloP[2];

  //printf("DigitP %f %f %f \n",gloP[0],gloP[1],gloP[2]);
  //printf("DigitD %f %f \n",gloD[0],gloD[1]);

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
