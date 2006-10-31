#include "MUONChamber.h"

#include <Alieve/MUONData.h>
#include <Alieve/MUONChamberData.h>

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

#include <TStyle.h>
#include <TColor.h>
#include <TMath.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// MUONChamber
//

ClassImp(MUONChamber)

//______________________________________________________________________
MUONChamber::MUONChamber(const Text_t* n, const Text_t* t) :
Reve::RenderElement(fFrameColor),
TNamed(n,t),
fMUONData(0),
fFrameColor((Color_t)2),
fRTS(1),
fChamberID(0),
fQuadSet1(n,t),
fQuadSet2(n,t),
fThreshold(0),
fMaxVal(1024)
{
  //
  // constructor
  //

  ComputeBBox();

}

//______________________________________________________________________
MUONChamber::~MUONChamber()
{
  //
  // destructor
  //

  if(fMUONData) fMUONData->DecRefCount();

}

//______________________________________________________________________
void MUONChamber::ComputeBBox()
{
  //
  // bounding box
  //

#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  bbox_init();
#else
  BBoxInit();
#endif
  
  fBBox[0] = -300.0;
  fBBox[1] = +300.0;
  fBBox[2] = -300.0;
  fBBox[3] = +300.0;
  fBBox[4] = -1800.0;
  fBBox[5] = -500.0;

  Float_t* b1 = fQuadSet1.AssertBBox();
  for(Int_t i=0; i<6; ++i) { b1[i] = fBBox[i]; }
  Float_t* b2 = fQuadSet2.AssertBBox();
  for(Int_t i=0; i<6; ++i) { b2[i] = fBBox[i]; }
  
}

//______________________________________________________________________
void MUONChamber::Paint(Option_t*)
{
  //
  // draw...
  //

  if(fRnrElement == kFALSE)
    return;

  TBuffer3D buffer(TBuffer3DTypes::kGeneric);

  buffer.fID           = this;
  buffer.fColor        = 2;
  buffer.fTransparency = 0;
  buffer.fLocalFrame   = 0;

  buffer.SetSectionsValid(TBuffer3D::kCore);
  Int_t reqSections = gPad->GetViewer3D()->AddObject(buffer);
  if (reqSections == TBuffer3D::kNone) {
    //printf("MUONChamber::Paint viewer was happy with Core buff3d.\n");
    return;
  }

  printf("MUONChamber::Paint only GL supported.\n");
  return;

}

//______________________________________________________________________
void MUONChamber::SetThreshold(Short_t t)
{
  //
  // digits amplitude threshold
  //

  fThreshold = TMath::Min(t, (Short_t)(fMaxVal - 1));
  ClearColorArray();
  IncRTS();

}

//______________________________________________________________________
void MUONChamber::SetMaxVal(Int_t mv)
{
  //
  // digits amplitude maximum value
  //

  fMaxVal = TMath::Max(mv, (Int_t)(fThreshold + 1));
  ClearColorArray();
  IncRTS();

}

//______________________________________________________________________
void MUONChamber::SetupColor(Int_t val, UChar_t* pixel) const
{
  //
  // RGBA color for amplitude "val"
  //

  Float_t div  = TMath::Max(1, fMaxVal - fThreshold);
  Int_t   nCol = gStyle->GetNumberOfColors();
  Int_t   cBin = (Int_t) TMath::Nint(nCol*(val - fThreshold)/div);

  ColorFromIdx(gStyle->GetColorPalette(TMath::Min(nCol - 1, cBin)), pixel);

}

//______________________________________________________________________
Int_t MUONChamber::ColorIndex(Int_t val) const
{
  //
  // index color 
  //

  if(val < fThreshold) val = fThreshold;
  if(val > fMaxVal)    val = fMaxVal;

  Float_t div  = TMath::Max(1, fMaxVal - fThreshold);
  Int_t   nCol = gStyle->GetNumberOfColors();
  Int_t   cBin = (Int_t) TMath::Nint(nCol*(val - fThreshold)/div);

  return gStyle->GetColorPalette(TMath::Min(nCol - 1, cBin));

}

//______________________________________________________________________
void MUONChamber::SetupColorArray() const
{
  //
  // build array of colors
  //

  if(fColorArray)
    return;

  fColorArray = new UChar_t [4 * (fMaxVal - fThreshold + 1)];
  UChar_t* p = fColorArray;
  for(Int_t v=fThreshold; v<=fMaxVal; ++v, p+=4)
    SetupColor(v, p);

}

//______________________________________________________________________
void MUONChamber::ClearColorArray()
{
  //
  // delete array of colors
  //

  if(fColorArray) {
    delete [] fColorArray;
    fColorArray = 0;
  }
}

//______________________________________________________________________
void MUONChamber::SetDataSource(MUONData* data)
{

  if (data == fMUONData) return;
  if(fMUONData) fMUONData->DecRefCount();
  fMUONData = data;
  if(fMUONData) fMUONData->IncRefCount();
  IncRTS();

}

//______________________________________________________________________
MUONChamberData* MUONChamber::GetChamberData() const
{
  
  return fMUONData ? fMUONData->GetChamberData(fChamberID) : 0;

}

//______________________________________________________________________
void MUONChamber::UpdateQuads()
{

  fQuadSet1.Quads().clear();
  fQuadSet2.Quads().clear();

  MUONChamberData* data = GetChamberData();
  
  Float_t *buffer;
  Float_t x0, y0, x1, y1, z;
  Int_t charge, cathode;
  if (data != 0) {

    SetupColorArray();

    Int_t ndigits = data->GetNDigits();
    
    for (Int_t id = 0; id < ndigits; id++) {

      buffer = data->GetDigitBuffer(id);

      x0 = buffer[0]-buffer[2];
      y0 = buffer[1]-buffer[3];
      x1 = buffer[0]+buffer[2];
      y1 = buffer[1]+buffer[3];
      z  = buffer[4];
      charge = (Int_t)buffer[5];
      cathode = (Int_t)buffer[6];
      
      if (cathode == 0) {

	fQuadSet1.Quads().push_back(Reve::Quad());
	
	fQuadSet1.Quads().back().ColorFromIdx(ColorIndex(charge));
	//ColorFromArray(charge,(UChar_t*)&fQuadSet1.fQuads.back().color);
	
	//UChar_t* c = (UChar_t*)&fQuadSet1.fQuads.back().color; 
	//printf("%d %d %d %d \n",c[0],c[1],c[2],c[3]);
	
	Float_t* p = fQuadSet1.Quads().back().vertices;
	
	p[0] = x0;  p[1] = y0;  p[2] = z;  p += 3;
	p[0] = x1;  p[1] = y0;  p[2] = z;  p += 3;
	p[0] = x1;  p[1] = y1;  p[2] = z;  p += 3;
	p[0] = x0;  p[1] = y1;  p[2] = z;  p += 3;
	
      }

      if (cathode == 1) {

	fQuadSet2.Quads().push_back(Reve::Quad());
	
	fQuadSet2.Quads().back().ColorFromIdx(ColorIndex(charge));
	//ColorFromArray(charge,(UChar_t*)&fQuadSet2.fQuads.back().color);
	
	//UChar_t* c = (UChar_t*)&fQuadSet2.fQuads.back().color; 
	//printf("%d %d %d %d \n",c[0],c[1],c[2],c[3]);
	
	Float_t* p = fQuadSet2.Quads().back().vertices;
	
	p[0] = x0;  p[1] = y0;  p[2] = z;  p += 3;
	p[0] = x1;  p[1] = y0;  p[2] = z;  p += 3;
	p[0] = x1;  p[1] = y1;  p[2] = z;  p += 3;
	p[0] = x0;  p[1] = y1;  p[2] = z;  p += 3;
	
      }

    } // end digits loop

  }

}

//______________________________________________________________________
void MUONChamber::SetChamberID(Int_t id)
{

  if (id <  0) id = 0;
  if (id > 13) id = 13;

  fChamberID = id;
  IncRTS();

}

