// $Header$

#include <Reve/Pad.h>

#include <THashList.h>

//______________________________________________________________________
// Pad
//
// This was intended as a TPad wrapper to allow smart updates of
// groups of pads. Uses THashList instead of TList for faster removal
// of objects from the pad.
// Currently not used ...

using namespace Reve;

ClassImp(Pad)

Pad::Pad()
{
  fPrimitives = new THashList;
}

Pad::Pad(const char *name, const char *title, Double_t xlow,
	 Double_t ylow, Double_t xup, Double_t yup,
	 Color_t color, Short_t bordersize, Short_t bordermode)
  : TPad(name,title,xlow,ylow,xup,yup,color,bordersize,bordermode)
{
  delete fPrimitives;
  fPrimitives = new THashList;
}
