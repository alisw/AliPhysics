// $Header$

#include "ITSModuleStepper.h"
#include "ITSDigitsInfo.h"
#include "ITSModule.h"

#include "Reve/RGTopFrame.h"
#include "Reve/GridStepper.h"

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// ITSModuleStepper
//

ClassImp(ITSModuleStepper)


ITSModuleStepper::ITSModuleStepper(ITSDigitsInfo* di):
  RenderElementList("ITS 2DStore", "ITSModuleStepper"),
  fDigitsInfo(di),
  fStepper(0)
{
  fStepper = new GridStepper();
}

/**************************************************************************/
ITSModuleStepper::~ITSModuleStepper()
{
  delete fStepper;
}

/**************************************************************************/
void ITSModuleStepper::SetStepper(Int_t nx, Int_t ny, Float_t dx, Float_t dy)
{
  fStepper->SetNs(nx, ny, 1);
  RenderElement::DestroyElements();

  Int_t nmod = nx*ny;
  for(Int_t m = 0; m<nmod; m++) 
  {
    AddElement( new ITSModule(m, fDigitsInfo));
  }

  if(dx > 0 && dy > 0)
    fStepper->SetDs(dx, dy);
}

/**************************************************************************/
void  ITSModuleStepper::Start()
{
  fPosition = fIDs.begin();
  fStepper->Reset(); 
  Apply();
}

void  ITSModuleStepper::Next()
{
  if(fPosition == fIDs.end())
    return;
  fStepper->Reset(); 
  Apply();
}

/**************************************************************************/
void  ITSModuleStepper::Apply()
{
  for(List_i  childit=fChildren.begin();  childit!=fChildren.end(); ++childit)
  {
    if(fPosition != fIDs.end()) 
    {
      ITSModule* mod = dynamic_cast<ITSModule*>(*childit);
      mod->SetID(*fPosition); 
      ZTrans& mx = mod->RefHMTrans();

      Float_t dx, dy, dz;
      mod->GetFrameDimensions(dx, dy, dz);
      Double_t sh = fStepper->Dy;
      Double_t sw = (dx*fStepper->Dy)/dz;
      if(sw > fStepper->Dx)
      {
	sw =  fStepper->Dx;
	sh =  (dz*fStepper->Dx)/dx;
      }
      mx.UnitTrans();
      mx.RotateLF(3,2,TMath::PiOver2());
      mx.Scale(sw/dx, sh/dz,1);
      // mx.Scale(fStepper->Dx/dx, fStepper->Dy/dz,1);
      fStepper->SetTransAdvance(&mx);
      mod->SetRnrSelf(kTRUE);      
      fPosition++;
    }
    else {
      (*childit)->SetRnrSelf(kFALSE);
    }
  }
  gReve->Redraw3D();
}
