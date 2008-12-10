// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveITSModuleStepper.h"
#include "AliEveITSDigitsInfo.h"
#include "AliEveITSScaledModule.h"

#include <TEveManager.h>
#include <TEveGedEditor.h>
#include <TEveGridStepper.h>
#include <TEveTrans.h>

#include <TGLRnrCtx.h>
#include <TGLIncludes.h>
#include <TGLSelectRecord.h>
#include <TGLUtil.h>
#include <TGLViewer.h>
#include <TGLAxis.h>

#include <TMath.h>
#include <THLimitsFinder.h>
#include <TVirtualPad.h>

#include <RVersion.h>

//______________________________________________________________________________
//
// Display scaled ITS modules in a paged layout, also providing
// GL-overaly control GUI.


ClassImp(AliEveITSModuleStepper)

AliEveITSModuleStepper::AliEveITSModuleStepper(AliEveITSDigitsInfo* di) :
  TEveElementList("ITS 2DStore", "AliEveITSModuleStepper", kTRUE),

  fDigitsInfo(di),
  fScaleInfo(0),
  fStepper(0),

  fModuleIDs(),
  fPosition(0),
  fSubDet(-1),

  fModuleFont(), fTextFont(), fSymbolFont(),
  fAxis(0),

  fMenuHeight(0.13),
  fTextSize(64),
  fTextCol(kGray+1),
  fActiveCol(kRed-4),

  fActiveID(-1)
{
  // Constructor.

  SetMainColorPtr(&fTextCol);
  fAxis = new TGLAxis();

  // override member from base TEveElementList
  fChildClass = AliEveITSScaledModule::Class();

  fDigitsInfo->IncRefCount();

  fStepper = new TEveGridStepper();
  fStepper->SetNs(5, 4);

  fScaleInfo = new AliEveDigitScaleInfo();
  fScaleInfo->IncRefCount();

  gEve->GetDefaultGLViewer()->AddOverlayElement(this);
}

AliEveITSModuleStepper::~AliEveITSModuleStepper()
{
  // Destructor.

  gEve->GetDefaultGLViewer()->RemoveOverlayElement(this);

  fScaleInfo->DecRefCount();
  fDigitsInfo->DecRefCount();

  delete fStepper;
  delete fAxis;
}

/******************************************************************************/

void AliEveITSModuleStepper::Capacity()
{
  // Make sure we have just enough children (module representations)
  // to store as many modules as required by the grid-stepper
  // configuration.

  Int_t n = fStepper->GetNx()*fStepper->GetNy();
  if (n != NumChildren())
  {
    DestroyElements();
    for (Int_t m=0; m<n; ++m)
    {
      AddElement(new AliEveITSScaledModule(m, fDigitsInfo, fScaleInfo));
    }
  }
}

/******************************************************************************/

void AliEveITSModuleStepper::SetFirst(Int_t first)
{
  // Se module ID which apply to first item in stepper.

  Int_t lastpage = fModuleIDs.size()/Nxy();
  if (fModuleIDs.size() % Nxy() ) lastpage++;

  Int_t firstLastpage = (lastpage - 1)*Nxy();
  if (first > firstLastpage) first = firstLastpage;
  if (first < 0) first = 0;
  fPosition = first;
  Apply();
}

void AliEveITSModuleStepper::Start()
{
  // Go to first page.

  fPosition = 0;
  Apply();
}

void AliEveITSModuleStepper::Next()
{
  // Go to next page.

  SetFirst(fPosition + Nxy());
}

void AliEveITSModuleStepper::Previous()
{
  // Go to previous page.

  SetFirst(fPosition - Nxy());
}

void AliEveITSModuleStepper::End()
{
  // Go to last page.

  Int_t lastpage = fModuleIDs.size()/Nxy();
  if (fModuleIDs.size() % Nxy()) lastpage++;
  fPosition = (lastpage - 1)*Nxy();

  fStepper->Reset();
  Apply();
}

/******************************************************************************/

void AliEveITSModuleStepper::DisplayDet(Int_t det, Int_t layer)
{
  // Select modules to display by sub-det type / layer.

  fSubDet = det;
  fModuleIDs.clear();
  AliEveITSModuleSelection sel = AliEveITSModuleSelection();
  sel.SetType (det);
  sel.SetLayer(layer);
  fDigitsInfo->GetModuleIDs(&sel, fModuleIDs);
  //in reder menu define a space between left and right pager
  Start();
}

/******************************************************************************/

Int_t AliEveITSModuleStepper::GetCurrentPage() const
{
  // Get number of current page.

  Int_t idx = fPosition + 1;
  Int_t n   = idx/Nxy();
  if (idx % Nxy()) n++;
  return n;
}

/******************************************************************************/

Int_t AliEveITSModuleStepper::GetPages()
{
  // Get number of all pages.

  Int_t n = fModuleIDs.size()/Nxy();
  if (fModuleIDs.size() % Nxy()) n++;
  return n;
}

/******************************************************************************/

void  AliEveITSModuleStepper::Apply()
{
  // Apply current settings to children modules.

  gEve->DisableRedraw();
  Capacity();

  UInt_t idx = fPosition;
  for(List_i childit=fChildren.begin(); childit!=fChildren.end(); ++childit)
  {
    if (idx < fModuleIDs.size())
    {
      AliEveITSScaledModule* mod = dynamic_cast<AliEveITSScaledModule*>(*childit);
      mod->SetID(fModuleIDs[idx], kFALSE);
      TEveTrans& tr = mod->RefMainTrans();
      tr.UnitTrans();
      tr.RotateLF(3,2,TMath::PiOver2());
      tr.RotateLF(1,3,TMath::PiOver2());

      // scaling
      Float_t mz, mx;
      Float_t* fp = mod->GetFrame()->GetFramePoints();
      // switch x,z it will be rotated afterwards
      mx = -2*fp[0];
      mz = -2*fp[2];

      // fit width first
      Double_t sx = fStepper->GetDx();
      Double_t sy = (mx*fStepper->GetDx())/mz;
      if (sy > fStepper->GetDy())
      {
	sy =  fStepper->GetDy();
	sx =  (mz*fStepper->GetDx())/mx;
      }
      Float_t scale = (0.85*sx)/mz;
      tr.Scale(scale, scale, scale);

      Float_t  p[3];
      fStepper->GetPosition(p);
      tr.SetPos(p[0]+0.5*fStepper->GetDx(), p[1]+0.5*fStepper->GetDy(), p[2]+0.5*fStepper->GetDz());

      if (mod->GetSubDetID() == 2)
	mod->SetName(Form("SSD %d", idx));
      else if (mod->GetSubDetID() == 1)
	mod->SetName(Form("SDD %d", idx));
      else
	mod->SetName(Form("SPD %d", idx));
      mod->SetRnrSelf(kTRUE);

      fStepper->Step();
      idx++;
    }
    else {
      (*childit)->SetRnrSelf(kFALSE);
    }
  }

  fStepper->Reset();
  ElementChanged();
  gEve->EnableRedraw();
}


/******************************************************************************/
// Virtual event handlers from TGLOverlayElement
/******************************************************************************/

//______________________________________________________________________________
Bool_t AliEveITSModuleStepper::Handle(TGLRnrCtx          & /*rnrCtx*/,
                                      TGLOvlSelectRecord & rec,
                                      Event_t            * event)
{
  // Handle overlay event.
  // Return TRUE if event was handled.

  switch (event->fType)
  {
    case kMotionNotify:
    {
      Int_t item = rec.GetN() < 2 ? -1 : (Int_t)rec.GetItem(1);
      if (fActiveID != item) {
        fActiveID = item;
        return kTRUE;
      } else {
        return kFALSE;
      }
      break;
    }
    case kButtonPress:
    {
      if (event->fCode != kButton1) {
        return kFALSE;
      }
      switch (rec.GetItem(1))
      {
        case 1:
          Previous();
          break;
        case 2:
          Start();
          break;
        case 3:
          Next();
          break;
        case 4:
          End();
          break;
        case 5:
        {
          AliEveDigitScaleInfo* si = fScaleInfo;
          if (si->GetScale() < 5)
          {
            si->ScaleChanged(si->GetScale() + 1);
            ElementChanged(kTRUE, kTRUE);
          }
          break;
        }
        case 6:
        {
          AliEveDigitScaleInfo* si = fScaleInfo;
          if (si->GetScale() > 1)
          {
            si->ScaleChanged(si->GetScale() - 1);
            ElementChanged(kTRUE, kTRUE);
          }
          break;
        }
        case 7:
          gEve->GetEditor()->DisplayElement(*BeginChildren());
          break;

        case 8:
	  DisplayDet(0, -1);
          break;
        case 9:
	  DisplayDet(1, -1);
          break;
        case 10:
	  DisplayDet(2, -1);
          break;
        default:
          break;
      }
      return kTRUE;
      break;
    }
    default:
      break;
  } // end switch
  return kFALSE;
}

//______________________________________________________________________________
Bool_t AliEveITSModuleStepper::MouseEnter(TGLOvlSelectRecord& /*rec*/)
{
  // Mouse has entered overlay area.

  return kTRUE;
}

//______________________________________________________________________________
void AliEveITSModuleStepper::MouseLeave()
{
  // Mouse has left overlay area.

  fActiveID = -1;
}


/******************************************************************************/
// Protected sub-renderers
/******************************************************************************/

//______________________________________________________________________________
void AliEveITSModuleStepper::RenderText(const char* txt, Int_t id, const TGLFont &font, Float_t step)
{
  // Render text for button id.

  Float_t llx, lly, llz, urx, ury, urz;
  font.BBox(txt, llx, lly, llz, urx, ury, urz);
  (fActiveID == id && id > 0) ? TGLUtil::Color(fActiveCol) :TGLUtil::Color(fTextCol);

  if (step > 0)
  {
    // center text in the step interval
    glPushMatrix();
    if (step>urx)
    glTranslatef((step-urx+llx)*0.5f-llx, 0, 0);
    glLoadName(id);
    font.Render(txt);
    glPopMatrix();
    glTranslatef(step, 0, 0);
  }
  else 
  {
    glLoadName(id);
    glPushMatrix();
    font.Render(txt);
    glPopMatrix();
    glTranslatef(urx, 0, 0);
  }
}

//______________________________________________________________________________
void AliEveITSModuleStepper::RenderPalette(TEveRGBAPalette* p)
{
  // Render color palette with number axis.

  Float_t length = 7*fTextSize;
  Float_t x = 1.5*fTextSize;
  Float_t y = 0.2*fTextSize;

  glTranslatef(x, 0.8*fTextSize, 0);

  TGLCapabilitySwitch lights_off(GL_LIGHTING, kFALSE);

  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDisable(GL_CULL_FACE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

  glBegin(GL_QUAD_STRIP);
  TGLUtil::Color4ubv(p->ColorFromValue(p->GetMinVal()));
  glVertex2f(0, 0);
  glVertex2f(0, y);
  if (p->GetMaxVal() > p->GetMinVal() + 1)
  {
    Float_t xs = length/(p->GetMaxVal() - p->GetMinVal());
    Float_t x0 = xs;
    for(Int_t i=p->GetMinVal() + 1; i<p->GetMaxVal(); i++)
    {
      TGLUtil::Color4ubv(p->ColorFromValue(i));
      glVertex2f(x0, 0);
      glVertex2f(x0, y);
      x0+=xs;
    }
  }
  TGLUtil::Color4ubv(p->ColorFromValue(p->GetMaxVal()));
  glVertex2f(length, 0);
  glVertex2f(length, y);
  glEnd();

  glRotatef(-90, 1, 0, 0 );
  Double_t v1[3] = {0., 0., 0.};
  Double_t v2[3] = {length, 0, 0.};
  fAxis->SetTextColor(kGray+1);
  fAxis->SetLineColor(kGray+1);
  fAxis->PaintGLAxis(v1, v2, p->GetMinVal(), p->GetMaxVal(), 5);
  glPopAttrib();
}

//______________________________________________________________________________
void AliEveITSModuleStepper::RenderMenu(Int_t curP, Int_t maxP, Int_t scaleX, Int_t scaleZ)
{
  // Make UI to set page in stepper and UI to scale in the AliEveITSScaledModule.

  TGLUtil::Color(fTextCol);
  fTextFont.PreRender();
  glTranslatef(0, fTextSize*0.3, 0);
  {
    // pager
    glTranslatef(fTextSize*0.2, 0 , 0);
    RenderText("9", 2, fSymbolFont); // last page
    RenderText("3", 1, fSymbolFont);//last page
    RenderText(Form("%d/%d", curP, maxP),-1, fTextFont, 2.7*fTextSize); //status
    {
      // bugg in webdings font , bbox does not give realistic value
      Float_t llx, lly, llz, urx, ury, urz;
      fSymbolFont.BBox("4", llx, lly, llz, urx, ury, urz);
      glTranslatef(-llx, 0, 0);
    }
    RenderText("4", 3, fSymbolFont); // next page
    RenderText(":",4, fSymbolFont); // last page
  }
  {
    // scale
    glTranslatef(fTextSize,0, 0);
    RenderText(Form("Zoom:"), -1, fTextFont);
    RenderText("6", 6, fSymbolFont);
    RenderText("5", 5, fSymbolFont);
    RenderText(Form("%dx%d", scaleX, scaleZ), -1, fTextFont, 2*fTextSize);
  }
  {
    // detectors
    glTranslatef(fTextSize, 0, 0);
    RenderText("SPD ", 8, fTextFont);
    RenderText("SDD ", 9, fTextFont);
    RenderText("SSD ", 10, fTextFont);
    fTextFont.PostRender();
  }
}

//______________________________________________________________________________
void AliEveITSModuleStepper::RenderModuleIDs()
{
  // Render module-ids.

  Double_t x, y, z;
  UInt_t idx = fPosition;
  Float_t llx, lly, llz, urx, ury, urz;
  fModuleFont.PreRender();
  TGLUtil::Color(kWhite);
  for (List_i childit=fChildren.begin(); childit!=fChildren.end(); ++childit)
  {
    if (idx < fModuleIDs.size())
    {
      AliEveITSScaledModule* mod = dynamic_cast<AliEveITSScaledModule*>(*childit);
      TEveTrans& tr = mod->RefMainTrans();
      tr.GetPos(x,y,z);
      x += fStepper->GetDx()*0.5;
      y -= fStepper->GetDy()*0.5;
      z += 0.4; // !!! MT hack - cross check with overlay rendering.
      const char* txt = Form("%d",mod->GetID());
      fModuleFont.BBox(txt, llx, lly, llz, urx, ury, urz);
      glRasterPos3f(x, y, z);
      glBitmap(0, 0, 0, 0,-urx, 0, 0);
      fModuleFont.Render(txt);
      idx++;
    }
  }
  fModuleFont.PostRender();
}

/******************************************************************************/

void AliEveITSModuleStepper::Render(TGLRnrCtx& rnrCtx)
{
  // Render the overlay elements.

  AliEveITSScaledModule* sm = dynamic_cast<AliEveITSScaledModule*>(*BeginChildren());
  Int_t scaleIdx = fScaleInfo->GetScale() - 1;
  Int_t cnx = 0, cnz = 0;
  switch(sm->GetSubDetID())
  {
    case 0:
      cnx = fDigitsInfo->fSPDScaleX[scaleIdx];
      cnz = fDigitsInfo->fSPDScaleZ[scaleIdx];
      break;
    case 1:
      cnx = fDigitsInfo->fSDDScaleX[scaleIdx];
      cnz = fDigitsInfo->fSDDScaleZ[scaleIdx];
      break;
    case 2:
      cnx = fDigitsInfo->fSSDScale[scaleIdx];
      cnz = 1;
      break;
  }

  // init fonts
  if (fTextFont.GetMode() == TGLFont::kUndef)
  {
#if ROOT_VERSION_CODE >= 332547
    rnrCtx.RegisterFont(fTextSize, 4, TGLFont::kTexture, fTextFont);
    rnrCtx.RegisterFont(72,       31, TGLFont::kTexture, fSymbolFont);
    rnrCtx.RegisterFont(14,        4, TGLFont::kPixmap,  fModuleFont);
#else
    fTextFont = rnrCtx.GetFont(fTextSize, 4, TGLFont::kTexture);
    fSymbolFont =  rnrCtx.GetFont(72, 31, TGLFont::kTexture);
    fModuleFont =  rnrCtx.GetFont(14, 4, TGLFont::kPixmap);
#endif
  }

  {
    // toolbar
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    if (rnrCtx.Selection())
    {
      TGLRect rect(*rnrCtx.GetPickRectangle());
      rnrCtx.GetCamera()->WindowToViewport(rect);
      gluPickMatrix(rect.X(), rect.Y(), rect.Width(), rect.Height(),
                    (Int_t*) rnrCtx.GetCamera()->RefViewport().CArr());
    }
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(-1, -1, 0); // translate to lower left corner
    Float_t txtScale = fMenuHeight/fTextSize*0.5; // scale text
    glScalef(txtScale, txtScale, 1.);

    //menu
    glPushName(0);
    RenderMenu(GetCurrentPage(), GetPages(), cnx, cnz);
    glPopName();
    //palette
    Double_t labelSize = 1.6*txtScale*fTextSize;
    fAxis->SetLabelsSize(labelSize);
    fAxis->SetLabelsOffset(1.2*labelSize);
    RenderPalette(sm->GetPalette());

    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
  }
  RenderModuleIDs();
}
