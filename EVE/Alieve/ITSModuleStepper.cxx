// $Header$

#include "ITSModuleStepper.h"
#include "ITSDigitsInfo.h"
#include "ITSScaledModule.h"

#include "Reve/RGTopFrame.h"
#include "Reve/RGEditor.h"
#include "Reve/GridStepper.h"
#include "Reve/GLTextNS.h"

#include <TObject.h>
#include <TSystem.h>

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

using namespace Reve;
using namespace Alieve;

//______________________________________________________________________
// ITSModuleStepper
//

ClassImp(ITSModuleStepper)

ITSModuleStepper::ITSModuleStepper(ITSDigitsInfo* di):
    RenderElement(fWColor),
    TNamed("ITS 2DStore", "ITSModuleStepper"),   

    fDigitsInfo(di),
    fScaleInfo(0),
    fStepper(0),
    fExpand(0.85),

    fRnrFrame(kTRUE),

    fWCorner(PT_BottomLeft),
    fWWidth(0.05),
    fWHeight(0.2)
{
  fStepper = new GridStepper();
  
  fScaleInfo = new DigitScaleInfo();

   fWColor = 5;

   //  GLTextNS::LoadDefaultFont(Form("%s/icons/fontdefault.txf",gSystem->Getenv("REVESYS")));
  GLTextNS::LoadDefaultFont(Form("%s/icons/fonthelvetica34.txf",gSystem->Getenv("REVESYS")));
}

/**************************************************************************/
ITSModuleStepper::~ITSModuleStepper()
{
  delete fStepper;
  delete fScaleInfo;
}

/**************************************************************************/

void ITSModuleStepper::SetStepper(Int_t nx, Int_t ny, Float_t dx, Float_t dy)
{
  fStepper->SetNs(nx, ny, 1);
  RenderElement::DestroyElements();

  Int_t nmod = nx*ny;
  for(Int_t m = 0; m<nmod; m++) 
  {
    AddElement( new ITSScaledModule(m, fDigitsInfo, fScaleInfo));
  }

  if(dx > 0 && dy > 0)
    fStepper->SetDs(dx, dy);
}

/**************************************************************************/

void  ITSModuleStepper::SetFirst(Int_t first)
{
  Int_t lastpage = fIDs.size()/Nxy();
  if(fIDs.size() % Nxy() ) lastpage++;  
        
  Int_t first_lastpage = (lastpage -1)*Nxy();
  if(first > first_lastpage) first = first_lastpage;
  if(first < 0) first = 0;
  fPosition = first;
  fStepper->Reset(); 
  Apply();
}

void  ITSModuleStepper::Start()
{
  fPosition = 0;
  fStepper->Reset(); 
  Apply();
}

void  ITSModuleStepper::Next()
{
  SetFirst( fPosition + Nxy());
}

void  ITSModuleStepper::Previous()
{
  // move to the top left corner first
  SetFirst( fPosition - Nxy());
}

void  ITSModuleStepper::End()
{ 
  Int_t lastpage = fIDs.size()/Nxy();
  if(fIDs.size() % Nxy() ) lastpage++;  
  fPosition = (lastpage -1)*Nxy();

  fStepper->Reset(); 
  Apply();
}

/**************************************************************************/

void  ITSModuleStepper::Apply()
{
  // printf("ITSModuleStepper::Apply fPosition %d \n", fPosition);
  gReve->DisableRedraw();

  UInt_t idx = fPosition;
  Float_t  p[3];

  for(List_i  childit=fChildren.begin();  childit!=fChildren.end(); ++childit)
  {
    if(idx < fIDs.size()) 
    {
      ITSScaledModule* mod = dynamic_cast<ITSScaledModule*>(*childit);
      mod->SetID(fIDs[idx], kFALSE); 
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
	// printf("fit width \n");
	sw =  fStepper->Dx;
	sh =  (dz*fStepper->Dx)/dx;
      }
      mx.UnitTrans();
      mx.RotateLF(3,2,TMath::PiOver2());
      mx.Scale(sw/dx, sh/dz,1);
     
      fStepper->GetPosition(p);
      mx.SetPos(p[0]+0.5*fStepper->Dx, p[1]+0.5*fStepper->Dy, p[2]+0.5*fStepper->Dz);
     
      mx.Scale(fExpand, fExpand,1);
      mx.RotateLF(2,1,TMath::PiOver2());
      
      mx.MultLeft(fHMTrans);

      mod->SetRnrSelf(kTRUE);
  
      if(mod->GetSubDetID() == 2)
	mod->SetName(Form("SSD %d", idx));
      else if(mod->GetSubDetID() == 1)
	mod->SetName(Form("SDD %d", idx));
      else
	mod->SetName(Form("SPD %d", idx));
      mod->UpdateItems();

      fStepper->Step();
      idx++;
    }
    else {
      (*childit)->SetRnrSelf(kFALSE);
    }
  }

  // update in case scaled module is a model in the editor
  gReve->EnableRedraw();
  gReve->GetEditor()->DisplayObject(gReve->GetEditor()->GetModel());
}

/**************************************************************************/

void  ITSModuleStepper::DisplayDet(Int_t det, Int_t layer)
{
  ITSModuleSelection sel = ITSModuleSelection();
  sel.fType = det; sel.fLayer=layer;
  fDigitsInfo->GetModuleIDs(&sel, fIDs);
  Start();
}

/**************************************************************************/

void  ITSModuleStepper::DisplayTheta(Float_t min, Float_t max)
{
  fIDs.clear();
  ITSModuleSelection sel = ITSModuleSelection();
  sel.fMaxTheta = max; sel.fMinTheta=min; 
  fDigitsInfo->GetModuleIDs(&sel, fIDs);
  Start();
}

/**************************************************************************/

void ITSModuleStepper::Paint(Option_t* /*option*/)
{
  static const Exc_t eH("ITSModuleStepper::Paint ");

  TBuffer3D buff(TBuffer3DTypes::kGeneric);

  // Section kCore
  buff.fID           = this;
  buff.fColor        = fWColor;
  buff.fTransparency = 0;
  fHMTrans.SetBuffer3D(buff);
  buff.SetSectionsValid(TBuffer3D::kCore);

  Int_t reqSections = gPad->GetViewer3D()->AddObject(buff);
  if (reqSections != TBuffer3D::kNone)
    Error(eH, "only direct GL rendering supported.");
}

/**************************************************************************/

void ITSModuleStepper::ComputeBBox()
{
  // printf("ITSModuleStepper::ComputeBBox \n");
  BBoxInit();
 
  Float_t W = fStepper->Dx*fStepper->Nx;
  Float_t H = fStepper->Dy*fStepper->Ny;
  
  BBoxCheckPoint(0, 0, -0.1);
  BBoxCheckPoint(W, H,  0.1);
  
  // Float_t dx = W*fWWidth);
  Float_t dy = H*fWHeight;

  if(fWCorner == PT_BottomLeft || fWCorner == PT_BottomRight) {
    BBoxCheckPoint(0, -dy, 0);
  }
  else 
    BBoxCheckPoint(0,  H+dy, 0);
}


/**************************************************************************/

Int_t ITSModuleStepper::GetCurrentPage()
{
  Int_t idx = fPosition +1; 
  Int_t n = idx/Nxy();
  if(idx % Nxy()) n++;
  return n;
}
/**************************************************************************/

Int_t ITSModuleStepper::GetPages()
{
  Int_t n = fIDs.size()/Nxy(); 
  if(fIDs.size() % Nxy()) n++; 
  return n;
}
  
