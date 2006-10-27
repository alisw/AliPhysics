#include "TRDData.h"
#include "TRDModuleImp.h"

#include "AliTRDcalibDB.h"
#include "AliTRDpadPlane.h"
#include "AliTRDgeometry.h"

using namespace Reve;
using namespace Alieve;
using namespace std;

ClassImp(TRDDigits)
ClassImp(TRDHits)

///////////////////////////////////////////////////////////
/////////////   TRDDigits             /////////////////////
///////////////////////////////////////////////////////////

//________________________________________________________
TRDDigits::TRDDigits(TRDChamber *p): QuadSet("digits", ""), RenderElement()
{
	fChamber = p;
	
	kLog = kFALSE;
	kBox  = kFALSE;
	fThreshold = 10;
}

//________________________________________________________
void	TRDDigits::SetData(AliTRDdataArrayI *digits)
{
	
	fData.Allocate(fChamber->rowMax, fChamber->colMax, fChamber->timeMax);
	digits->Expand();
	for (Int_t  row = 0;  row <  fChamber->rowMax;  row++)
		for (Int_t  col = 0;  col <  fChamber->colMax;  col++)
			for (Int_t time = 0; time < fChamber->timeMax; time++) {
//		if(digits->GetDataUnchecked(row, col, time) > 20) printf("%d %d %d %d\n", row, col, time, digits->GetDataUnchecked(row, col, time));
		fData.SetDataUnchecked(row, col, time, digits->GetDataUnchecked(row, col, time));
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
	for (Int_t  row = 0;  row <  fChamber->rowMax;  row++) {
		rowSize = .5 * fChamber->fPadPlane->GetRowSize(row);
		z = fChamber->fPadPlane->GetRowPos(row) - rowSize;
		
		for (Int_t  col = 0;  col <  fChamber->colMax;  col++) {
			colSize = .5 * fChamber->fPadPlane->GetColSize(col);
			y = fChamber->fPadPlane->GetColPos(col) - colSize;
			t0 = calibration->GetT0(fChamber->fDet, col, row);
			timeBinSize = calibration->GetVdrift(fChamber->fDet, col, row)/fChamber->samplingFrequency;
			
			for (Int_t time = 0; time < fChamber->timeMax; time++) {
				charge = fData.GetDataUnchecked(row, col, time);
	  		if (charge < fThreshold) continue;
				
				x = fChamber->fX0 - (time+0.5-t0)*timeBinSize;
				scale = kLog ? TMath::Log(float(charge))/TMath::Log(1024.) : charge/1024.;
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
				if(kBox){
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
					fChamber->fGeo->RotateBack(fChamber->fDet,cloc[ic],cglo);
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
	if(kBox) fBoxes.Paint(option);
	else QuadSet::Paint(option);
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
TRDHits::TRDHits(const Text_t* name, Int_t n_points):PointSet(name, n_points)
{

}

//________________________________________________________
void TRDHits::PointSelected(Int_t n)
{
	printf("void TRDHits::PointSelected(%d)\n", n);
//	printf("Detector %d\n", ((TRDChamber*)GetPointId(n))->GetDetector());
}

