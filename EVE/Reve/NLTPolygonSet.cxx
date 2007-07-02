#include "NLTPolygonSet.h"
#include "PODs.h"

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

#include <list>

using namespace Reve;

/**************************************************************************/

ClassImp(NLTPolygonSet)

/**************************************************************************/
NLTPolygonSet::NLTPolygonSet(const Text_t* n, const Text_t* t) :
 RenderElement(),
 TNamed(n,t), 
 TAtt3D(), 
 fNPnts(0),
 fPnts(0),
 fNPols(0),
 fPols(0),
 fFillColor(5),
 fLineColor(3),
 fLineWidth(1),
 fZDepth(0)
{}

/**************************************************************************/

NLTPolygonSet::~NLTPolygonSet()
{
  for(Int_t i=0; i<fNPols; i++)
    delete [] fPols[i].fPnts;

  delete [] fPols;
  delete [] fPnts;
}

/**************************************************************************/
void NLTPolygonSet::ComputeBBox()
{
  if(fNPols == 0) {
    BBoxZero();
    return;
  }
  
  BBoxInit();
  for(Int_t pi = 0; pi<fNPnts; pi++) {
    BBoxCheckPoint(fPnts[pi].x, fPnts[pi].y, fPnts[pi].z );
  }
}

/**************************************************************************/
void NLTPolygonSet::Paint(Option_t* )
{
  TBuffer3D buffer(TBuffer3DTypes::kGeneric);

  // Section kCore
  buffer.fID           = this;
  buffer.fColor        = 2;
  buffer.fTransparency = 0;
  buffer.fLocalFrame   = false; 

  buffer.SetSectionsValid(TBuffer3D::kCore);
   
  // We fill kCore on first pass and try with viewer
  Int_t reqSections = gPad->GetViewer3D()->AddObject(buffer);
  if (reqSections == TBuffer3D::kNone) {
    return;
  }
}

/**************************************************************************/
void NLTPolygonSet::Dump() const
{
  printf("NLTPolygonSet %d polygond\n", fNPols);
  for(Int_t p = 0; p<fNPols; p++) {
    printf("%d polygon %d points :\n", p, fPols[p].fNPnts);    
    for(Int_t i = 0; i<fPols[p].fNPnts; i++) {
      Int_t pi = fPols[p].fPnts[i];
      printf("(%f, %f, %f)", fPnts[pi].x, fPnts[pi].y, fPnts[pi].z);
    }
    printf("\n");
  }
}
