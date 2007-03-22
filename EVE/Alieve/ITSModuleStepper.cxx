// $Header$

#include "ITSModuleStepper.h"
#include "ITSDigitsInfo.h"
#include "ITSScaledModule.h"

#include "Reve/RGTopFrame.h"
#include "Reve/RGEditor.h"
#include "Reve/GridStepper.h"

#include <TObject.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// ITSModuleStepper
//

ClassImp(ITSModuleStepper)

ITSModuleStepper::ITSModuleStepper(ITSDigitsInfo* di):
  RenderElementList("ITS 2DStore", "ITSModuleStepper"),
  fDigitsInfo(di),
  fStepper(0),
  fExpand(0.85)
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
    AddElement( new ITSScaledModule(m, fDigitsInfo));
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
  // check editor
  for(List_i  childit=fChildren.begin();  childit!=fChildren.end(); ++childit)
  {
    if(fPosition != fIDs.end()) 
    {
      ITSScaledModule* mod = dynamic_cast<ITSScaledModule*>(*childit);
      mod->SetID(*fPosition, kFALSE); 
      ZTrans& mx = mod->RefHMTrans();

      Float_t dx, dz;
      Float_t* fp = mod->GetFrame()->GetFramePoints();
      // switch x,z it will be rotated afterwards
      dz = -2*fp[0];
      dx = -2*fp[2];

      Double_t sh = fStepper->Dy;
      Double_t sw = (dx*fStepper->Dy)/dz;
      if(sw > fStepper->Dx)
      {
        printf("fit width \n");
	sw =  fStepper->Dx;
	sh =  (dz*fStepper->Dx)/dx;
      }
      mx.UnitTrans();
      mx.RotateLF(3,2,TMath::PiOver2());
      mx.Scale(sw/dx, sh/dz,1);
      fStepper->SetTransAdvance(&mx);
      mx.Scale(fExpand, fExpand,1);
      mx.RotateLF(2,1,TMath::PiOver2());
      mod->SetRnrSelf(kTRUE);
  
      if(mod->GetSubDetID() == 2)
	mod->SetName(Form("SSD %d", *fPosition));
      else if(mod->GetSubDetID() == 1)
	mod->SetName(Form("SDD %d", *fPosition));
      else
	mod->SetName(Form("SPD %d", *fPosition));
      mod->UpdateItems();

      fPosition++;
    }
    else {
      (*childit)->SetRnrSelf(kFALSE);
    }
  }
  // update in case scaled module is a model in the editor
  gReve->GetEditor()->DisplayObject(gReve->GetEditor()->GetModel());
  gReve->Redraw3D();
}
