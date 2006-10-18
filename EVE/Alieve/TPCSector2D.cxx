// $Header$

#include "TPCSector2D.h"

#include <Alieve/TPCData.h>
#include <Alieve/TPCSectorData.h>

#include <AliTPCParam.h>

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

using namespace Reve;
using namespace Alieve;
using namespace std;

//______________________________________________________________________
// TPCSector2D
//
// Displays TPC raw-data in 2D.
//
// fShowMax: true  - display maximum value for given time interval
//           false - display integral
// fAverage: only available when fShowMax = false; divide by time window width
//
// fUseTexture: use OpenGL textures to display data (fast rendering,
//   updates take the same time)
//

ClassImp(TPCSector2D)

/**************************************************************************/

TPCSector2D::TPCSector2D(const Text_t* n, const Text_t* t) :
  TPCSectorViz(n,t),

  fShowMax (kTRUE),
  fAverage (kFALSE),

  fUseTexture (kTRUE),
  fPickEmpty  (kFALSE),
  fPickMode   (0)
{}

TPCSector2D::~TPCSector2D()
{}

/**************************************************************************/

void TPCSector2D::ComputeBBox()
{
  const TPCSectorData::SegmentInfo&  iSeg = TPCSectorData::GetInnSeg();
  const TPCSectorData::SegmentInfo& o2Seg = TPCSectorData::GetOut2Seg();

#if ROOT_VERSION_CODE <= ROOT_VERSION(5,11,2)
  bbox_init();
#else
  BBoxInit();
#endif
  Float_t w = o2Seg.GetNMaxPads()*o2Seg.GetPadWidth()/2;
  fBBox[0] = -w;
  fBBox[1] =  w;
  fBBox[2] =  iSeg.GetRLow();
  fBBox[3] =  o2Seg.GetRLow() + o2Seg.GetNRows()*o2Seg.GetPadHeight();
  fBBox[4] = -0.5; // Fake z-width to 1 cm.
  fBBox[5] =  0.5;
}

/**************************************************************************/

void TPCSector2D::Paint(Option_t* )
{
  if(fRnrElement == kFALSE)
    return;

  TBuffer3D buffer(TBuffer3DTypes::kGeneric);

  // Section kCore
  buffer.fID           = this;
  buffer.fColor        = 1;
  buffer.fTransparency = 0;
  fHMTrans.SetBuffer3D(buffer);
  buffer.SetSectionsValid(TBuffer3D::kCore);
   
  Int_t reqSections = gPad->GetViewer3D()->AddObject(buffer);
  if (reqSections == TBuffer3D::kNone) {
    return;
  }

  Error("TPCSector2D::Paint", "only direct OpenGL rendering supported.");
}
