// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#include "AliEveTRDData.h"
#include "AliEveTRDModuleImp.h"

#include "AliLog.h"
#include "AliTRDhit.h"
#include "AliTRDcluster.h"
#include "AliTRDcalibDB.h"
#include "AliTRDpadPlane.h"
#include "AliTRDgeometry.h"
#include "AliTRDdigitsManager.h"

using namespace std;

ClassImp(AliEveTRDHits)
ClassImp(AliEveTRDDigits)
ClassImp(AliEveTRDClusters)

///////////////////////////////////////////////////////////
/////////////   AliEveTRDDigits             /////////////////////
///////////////////////////////////////////////////////////

//________________________________________________________
AliEveTRDDigits::AliEveTRDDigits(AliEveTRDChamber *p): TEveQuadSet("digits", ""), fParent(p)
{}

//________________________________________________________
void	AliEveTRDDigits::SetData(AliTRDdigitsManager *digits)
{

	fData.Allocate(fParent->rowMax, fParent->colMax, fParent->timeMax);
//	digits->Expand();
	for (Int_t  row = 0;  row <  fParent->rowMax;  row++)
		for (Int_t  col = 0;  col <  fParent->colMax;  col++)
			for (Int_t time = 0; time < fParent->timeMax; time++) {
				if(digits->GetDigitAmp(row, col, time, fParent->GetID()) < 0) continue;
				fData.SetDataUnchecked(row, col, time, digits->GetDigitAmp(row, col, time, fParent->GetID()));
	}
}

//________________________________________________________
void AliEveTRDDigits::ComputeRepresentation()
{
  // Calculate digits representation according to user settings. The
  // user can set the following parameters:
  // - digits scale (log/lin)
  // - digits threshold
  // - digits apparence (quads/boxes)

  TEveQuadSet::Reset(TEveQuadSet::kQT_FreeQuad, kTRUE, 64);
  // MT fBoxes.fBoxes.clear();

  Double_t colSize, rowSize, scale;
  Double_t x, y, z;

  Int_t charge;
  Float_t t0;
  Float_t timeBinSize;

  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  Double_t cloc[4][3], cglo[3];
  Int_t color, dimension;
  fData.Expand();
  for (Int_t  row = 0;  row <  fParent->rowMax;  row++) {
    rowSize = .5 * fParent->fPadPlane->GetRowSize(row);
    z = fParent->fPadPlane->GetRowPos(row) - rowSize;

    for (Int_t  col = 0;  col <  fParent->colMax;  col++) {
      colSize = .5 * fParent->fPadPlane->GetColSize(col);
      y = fParent->fPadPlane->GetColPos(col) - colSize;
      t0 = calibration->GetT0(fParent->fDet, col, row);
      timeBinSize = calibration->GetVdrift(fParent->fDet, col, row)/fParent->samplingFrequency;

      for (Int_t time = 0; time < fParent->timeMax; time++) {
	charge = fData.GetDataUnchecked(row, col, time);
	if (charge < fParent->GetDigitsThreshold()) continue;

	x = fParent->fX0 - (time+0.5-t0)*timeBinSize;
	scale = fParent->GetDigitsLog() ? TMath::Log(float(charge))/TMath::Log(1024.) : charge/1024.;
	color  = 50+int(scale*50.);

	cloc[0][2] = z - rowSize * scale;
	cloc[0][1] = y - colSize * scale;
	cloc[0][0] = x;

	cloc[1][2] = z - rowSize * scale;
	cloc[1][1] = y + colSize * scale;
	cloc[1][0] = x;

	cloc[2][2] = z + rowSize * scale;
	cloc[2][1] = y + colSize * scale;
	cloc[2][0] = x;

	cloc[3][2] = z + rowSize * scale;
	cloc[3][1] = y - colSize * scale;
	cloc[3][0] = x;

	Float_t* p = 0;
	if( fParent->GetDigitsBox()){
	  // MT fBoxes.fBoxes.push_back(Box());
	  // MT fBoxes.fBoxes.back().color[0] = (UChar_t)color;
	  // MT fBoxes.fBoxes.back().color[1] = (UChar_t)color;
	  // MT fBoxes.fBoxes.back().color[2] = (UChar_t)color;
	  // MT fBoxes.fBoxes.back().color[3] = (UChar_t)color;
	  // MT p = fBoxes.fBoxes.back().vertices;
	  dimension = 2;
	} else {
	  AddQuad((Float_t*)0);
	  QuadColor(color);
	  p = ((QFreeQuad_t*) fLastDigit)->fVertices;
	  dimension = 1;
	}

	for(int id=0; id<dimension; id++)
	  for (Int_t ic = 0; ic < 4; ic++) {
	    cloc[ic][0] -= .5 * id * timeBinSize;
	    fParent->fGeo->RotateBack(fParent->fDet,cloc[ic],cglo);
	    p[0] = cglo[0]; p[1] = cglo[1]; p[2] = cglo[2];
	    p+=3;
	  }
      }  // end time loop
    }  // end col loop
  }  // end row loop
  fData.Compress(1);
}

//________________________________________________________
void AliEveTRDDigits::Paint(Option_t *option)
{
	if(fParent->GetDigitsBox()) fBoxes.Paint(option);
	else TEveQuadSet::Paint(option);
}

//________________________________________________________
void AliEveTRDDigits::Reset()
{
	TEveQuadSet::Reset(TEveQuadSet::kQT_FreeQuad, kTRUE, 64);
	// MT fBoxes.fBoxes.clear();
	fData.Reset();
}

///////////////////////////////////////////////////////////
/////////////   AliEveTRDHits               /////////////////////
///////////////////////////////////////////////////////////

//________________________________________________________
AliEveTRDHits::AliEveTRDHits(AliEveTRDChamber *p):TEvePointSet("hits", 20), fParent(p)
{}

//________________________________________________________
void AliEveTRDHits::PointSelected(Int_t n)
{
	fParent->SpawnEditor();
	AliTRDhit *h = dynamic_cast<AliTRDhit*>(GetPointId(n));
	printf("\nDetector             : %d\n", h->GetDetector());
	printf("Region of production : %c\n", h->FromAmplification() ? 'A' : 'D');
	printf("TR photon            : %s\n", h->FromTRphoton() ? "Yes" : "No");
	printf("Charge               : %d\n", h->GetCharge());
	printf("MC track label       : %d\n", h->GetTrack());
	printf("Time from collision  : %f\n", h->GetTime());
}


///////////////////////////////////////////////////////////
/////////////   AliEveTRDHits               /////////////////////
///////////////////////////////////////////////////////////

//________________________________________________________
AliEveTRDClusters::AliEveTRDClusters(AliEveTRDChamber *p):AliEveTRDHits(p)
{}

//________________________________________________________
void AliEveTRDClusters::PointSelected(Int_t n)
{
	fParent->SpawnEditor();
	AliTRDcluster *c = dynamic_cast<AliTRDcluster*>(GetPointId(n));
	printf("\nDetector             : %d\n", c->GetDetector());
	printf("Charge               : %f\n", c->GetQ());
	printf("Sum S                : %4.0f\n", c->GetSumS());
	printf("Time bin             : %d\n", c->GetLocalTimeBin());
	printf("Signals              : ");
	Short_t *cSignals = c->GetSignals();
	for(Int_t ipad=0; ipad<7; ipad++) printf("%d ", cSignals[ipad]); printf("\n");
	printf("Central pad          : %d\n", c->GetPadCol());
	printf("MC track labels      : ");
	for(Int_t itrk=0; itrk<3; itrk++) printf("%d ", c->GetLabel(itrk)); printf("\n");
// Bool_t	AliCluster::GetGlobalCov(Float_t* cov) const
// Bool_t	AliCluster::GetGlobalXYZ(Float_t* xyz) const
// Float_t	AliCluster::GetSigmaY2() const
// Float_t	AliCluster::GetSigmaYZ() const
// Float_t	AliCluster::GetSigmaZ2() const
}

///////////////////////////////////////////////////////////
/////////////   AliEveTRDHitsEditor         /////////////////////
///////////////////////////////////////////////////////////
AliEveTRDHitsEditor::AliEveTRDHitsEditor(const TGWindow* p, Int_t width, Int_t height, UInt_t options, Pixel_t back) : TGedFrame(p, width, height, options, back)
{
	MakeTitle("TRD Hits");

}

AliEveTRDHitsEditor::~AliEveTRDHitsEditor()
{}

void AliEveTRDHitsEditor::SetModel(TObject* obj)
{
	fM = dynamic_cast<AliEveTRDHits*>(obj);

// 	Float_t x, y, z;
// 	for(int ihit=0; ihit<fM->GetN(); ihit++){
// 		fM->GetPoint(ihit, x, y, z);
// 		printf("%3d : x=%6.3f y=%6.3f z=%6.3f\n", ihit, x, y, z);
// 	}
}

///////////////////////////////////////////////////////////
/////////////   AliEveTRDDigitsEditor       /////////////////////
///////////////////////////////////////////////////////////
AliEveTRDDigitsEditor::AliEveTRDDigitsEditor(const TGWindow* p, Int_t width, Int_t height, UInt_t options, Pixel_t back) : TGedFrame(p, width, height, options, back)
{
	MakeTitle("TRD Digits");

}

AliEveTRDDigitsEditor::~AliEveTRDDigitsEditor()
{}

void AliEveTRDDigitsEditor::SetModel(TObject* obj)
{
	fM = dynamic_cast<AliEveTRDDigits*>(obj);
	fM->fParent->SpawnEditor();

// 	printf("Chamber %d", fM->fParent->GetID());
// 	for (Int_t  row = 0;  row <  fM->fParent->GetRowMax();  row++)
// 		for (Int_t  col = 0;  col <  fM->fParent->GetColMax();  col++)
// 			for (Int_t time = 0; time < fM->fParent->GetTimeMax(); time++) {
// 				printf("\tA(%d %d %d) = %d\n", row, col, time, fM->fData.GetDataUnchecked(row, col, time));
// 			}
}
