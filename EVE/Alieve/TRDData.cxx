#include "TRDData.h"
#include "TRDModuleImp.h"

#include "AliLog.h"
#include "AliTRDhit.h"
#include "AliTRDcluster.h"
#include "AliTRDcalibDB.h"
#include "AliTRDpadPlane.h"
#include "AliTRDgeometry.h"
#include "AliTRDdigitsManager.h"

using namespace Reve;
using namespace Alieve;
using namespace std;


ClassImp(TRDHits)
ClassImp(TRDDigits)
ClassImp(TRDClusters)

///////////////////////////////////////////////////////////
/////////////   TRDDigits             /////////////////////
///////////////////////////////////////////////////////////

//________________________________________________________
TRDDigits::TRDDigits(TRDChamber *p): OldQuadSet("digits", ""), RenderElement(), fParent(p)
{}

//________________________________________________________
void	TRDDigits::SetData(AliTRDdigitsManager *digits)
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
void	TRDDigits::ComputeRepresentation()
{
// Calculate digits representation according to user settings. The
// user can set the following parameters:
// - digits scale (log/lin)
// - digits threshold
// - digits apparence (quads/boxes)

	fQuads.clear();
	fBoxes.fBoxes.clear();
		
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
	
				Float_t* p;
				if( fParent->GetDigitsBox()){
					fBoxes.fBoxes.push_back(Reve::Box());
					fBoxes.fBoxes.back().color[0] = (UChar_t)color;
					fBoxes.fBoxes.back().color[1] = (UChar_t)color;
					fBoxes.fBoxes.back().color[2] = (UChar_t)color;
					fBoxes.fBoxes.back().color[3] = (UChar_t)color;
					p = fBoxes.fBoxes.back().vertices;
					dimension = 2;
				} else {
					fQuads.push_back(Reve::Quad());
					fQuads.back().ColorFromIdx(color);
					p = fQuads.back().vertices;
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
void TRDDigits::Paint(Option_t *option)
{
	if(fParent->GetDigitsBox()) fBoxes.Paint(option);
	else OldQuadSet::Paint(option);
}

//________________________________________________________
void TRDDigits::Reset()
{
	fQuads.clear();
	fBoxes.fBoxes.clear();
	fData.Reset();
}

///////////////////////////////////////////////////////////
/////////////   TRDHits               /////////////////////
///////////////////////////////////////////////////////////

//________________________________________________________
TRDHits::TRDHits(TRDChamber *p):PointSet("hits", 20), fParent(p)
{}

//________________________________________________________
void TRDHits::PointSelected(Int_t n)
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
/////////////   TRDHits               /////////////////////
///////////////////////////////////////////////////////////

//________________________________________________________
TRDClusters::TRDClusters(TRDChamber *p):TRDHits(p)
{}

//________________________________________________________
void TRDClusters::PointSelected(Int_t n)
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
/////////////   TRDHitsEditor         /////////////////////
///////////////////////////////////////////////////////////
TRDHitsEditor::TRDHitsEditor(const TGWindow* p, Int_t width, Int_t height, UInt_t options, Pixel_t back) : TGedFrame(p, width, height, options, back)
{
	MakeTitle("TRD Hits");

}

TRDHitsEditor::~TRDHitsEditor()
{}

void TRDHitsEditor::SetModel(TObject* obj)
{
	fM = dynamic_cast<TRDHits*>(obj);

// 	Float_t x, y, z;
// 	for(int ihit=0; ihit<fM->GetN(); ihit++){
// 		fM->GetPoint(ihit, x, y, z);
// 		printf("%3d : x=%6.3f y=%6.3f z=%6.3f\n", ihit, x, y, z);
// 	}
}

///////////////////////////////////////////////////////////
/////////////   TRDDigitsEditor       /////////////////////
///////////////////////////////////////////////////////////
TRDDigitsEditor::TRDDigitsEditor(const TGWindow* p, Int_t width, Int_t height, UInt_t options, Pixel_t back) : TGedFrame(p, width, height, options, back)
{
	MakeTitle("TRD Digits");

}

TRDDigitsEditor::~TRDDigitsEditor()
{}

void TRDDigitsEditor::SetModel(TObject* obj)
{
	fM = dynamic_cast<TRDDigits*>(obj);
	fM->fParent->SpawnEditor();
	
// 	printf("Chamber %d", fM->fParent->GetID());
// 	for (Int_t  row = 0;  row <  fM->fParent->GetRowMax();  row++)
// 		for (Int_t  col = 0;  col <  fM->fParent->GetColMax();  col++)
// 			for (Int_t time = 0; time < fM->fParent->GetTimeMax(); time++) {
// 				printf("\tA(%d %d %d) = %d\n", row, col, time, fM->fData.GetDataUnchecked(row, col, time));
// 			}
}
