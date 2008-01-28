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
#include <TEveGLText.h>
#include <TEveTrans.h>

#include <TObject.h>
#include <TMath.h>

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

#include <TGLRnrCtx.h>
#include <TGLSelectRecord.h>
#include <TGLText.h>
// #include <FTFont.h>
#include <TGLAxis.h>
#include <TGLViewer.h>


//______________________________________________________________________________
//
// Display scaled ITS modules in a paged layout, also providing
// GL-overaly control GUI.


ClassImp(AliEveITSModuleStepper)

AliEveITSModuleStepper::AliEveITSModuleStepper(AliEveITSDigitsInfo* di) :
  TEveElementList("ITS 2DStore", "AliEveITSModuleStepper", kTRUE),

  fIDs(),
  fPosition(0),

  fDigitsInfo(di),
  fScaleInfo(0),

  fSubDet(-1),

  fStepper(0),
  fAxis(0),
  fText(0),
  fTextSize(0.05),
  fPagerGap(0.1),
  fRnrFrame(kFALSE),

  fExpandCell(0.85),
  fModuleFrameCol(2),

  fPaletteOffset(0.2),
  fPaletteLength(0.6),

  fWActive(-1),
  fWWidth(0.025),
  fWHeight(0.032),
  fWOff(0.05),
  fWCol(30),
  fWActiveCol(45),
  fFontCol(8)
{
  // Constructor.

  // override member from base TEveElementList
  fChildClass = AliEveITSScaledModule::Class();

  SetMainColorPtr(&fWCol);

  fDigitsInfo->IncRefCount();

  fStepper = new TEveGridStepper();
  fStepper->SetNs(5, 4);

  fScaleInfo = new AliEveDigitScaleInfo();
  fScaleInfo->IncRefCount();

  fAxis = new TGLAxis();
  fAxis->SetLineColor(4);
  fAxis->SetTextColor(fFontCol);

  fText = new TGLText();
  fText->SetTextColor(fFontCol);
  fText->SetGLTextFont(40);
  fText->SetGLTextAngles(0, 0, 0);
  fText->SetTextSize(fTextSize);

  gEve->GetGLViewer()->AddOverlayElement(this);
}

AliEveITSModuleStepper::~AliEveITSModuleStepper()
{
  // Destructor.

  gEve->GetGLViewer()->RemoveOverlayElement(this);

   fScaleInfo->DecRefCount();
  fDigitsInfo->DecRefCount();

  delete fStepper;

  delete fAxis;
  delete fText;
}

/******************************************************************************/

void AliEveITSModuleStepper::Capacity()
{
  // Make sure we have just enough children (module representations)
  // to store as many modules as required by the grid-stepper
  // configuration.

  Int_t N = fStepper->GetNx()*fStepper->GetNy();
  if (N != GetNChildren())
  {
    DestroyElements();
    for (Int_t m=0; m<N; m++)
    {
      AddElement(new AliEveITSScaledModule(m, fDigitsInfo, fScaleInfo));
    }
  }
}

/******************************************************************************/

void AliEveITSModuleStepper::SetFirst(Int_t first)
{
  Int_t lastpage = fIDs.size()/Nxy();
  if(fIDs.size() % Nxy() ) lastpage++;

  Int_t firstLastpage = (lastpage - 1)*Nxy();
  if(first > firstLastpage) first = firstLastpage;
  if(first < 0) first = 0;
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

  Int_t lastpage = fIDs.size()/Nxy();
  if (fIDs.size() % Nxy()) lastpage++;
  fPosition = (lastpage - 1)*Nxy();

  fStepper->Reset();
  Apply();
}

/******************************************************************************/

void AliEveITSModuleStepper::DisplayDet(Int_t det, Int_t layer)
{
  // Select modules to display by sub-det type / layer. 

  fSubDet = det;
  fIDs.clear();
  AliEveITSModuleSelection sel = AliEveITSModuleSelection();
  sel.SetType (det);
  sel.SetLayer(layer);
  fDigitsInfo->GetModuleIDs(&sel, fIDs);
  //in reder menu define a space between left and right pager
  fPagerGap = 1.2*TextLength(Form("%d/%d",GetPages(), GetPages()));
  Start();
}

/******************************************************************************/

void AliEveITSModuleStepper::DisplayTheta(Float_t min, Float_t max)
{
  // Select modules to display by theta range.

  fIDs.clear();
  AliEveITSModuleSelection sel = AliEveITSModuleSelection();
  sel.SetThetaRange(min, max);
  fDigitsInfo->GetModuleIDs(&sel, fIDs);
  Start();
}

/******************************************************************************/

Int_t AliEveITSModuleStepper::GetCurrentPage()
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

  Int_t n = fIDs.size()/Nxy();
  if(fIDs.size() % Nxy()) n++;
  return n;
}

/******************************************************************************/

void  AliEveITSModuleStepper::Apply()
{
  // Apply current settings to children modules.

  // printf("AliEveITSModuleStepper::Apply fPosition %d \n", fPosition);
  gEve->DisableRedraw();
  Capacity();

  UInt_t idx = fPosition;
  for(List_i childit=fChildren.begin(); childit!=fChildren.end(); ++childit)
  {
    if(idx < fIDs.size())
    {
      AliEveITSScaledModule* mod = dynamic_cast<AliEveITSScaledModule*>(*childit);
      mod->SetID(fIDs[idx], kFALSE);
      TEveTrans& tr = mod->RefHMTrans();
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
      if(sy > fStepper->GetDy())
      {
        //	printf("fit width \n");
	sy =  fStepper->GetDy();
	sx =  (mz*fStepper->GetDx())/mx;
      }
      Float_t scale = (fExpandCell*sx)/mz;
      tr.Scale(scale, scale, scale);

      Float_t  p[3];
      fStepper->GetPosition(p);
      tr.SetPos(p[0]+0.5*fStepper->GetDx(), p[1]+0.5*fStepper->GetDy(), p[2]+0.5*fStepper->GetDz());

      if(mod->GetSubDetID() == 2)
	mod->SetName(Form("SSD %d", idx));
      else if(mod->GetSubDetID() == 1)
	mod->SetName(Form("SDD %d", idx));
      else
	mod->SetName(Form("SPD %d", idx));
      mod->SetRnrSelf(kTRUE);
      mod->UpdateItems();

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

void AliEveITSModuleStepper::Render(TGLRnrCtx& rnrCtx)
{
  // Render the overlay elements.

  // render everyting in relative coordinates
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  if (rnrCtx.Selection())
  {
    // Should be
    // glLoadMatrix(rnrCtx.GetCamera()->GetProjMBase());
    TGLRect rect(*rnrCtx.GetPickRectangle());
    rnrCtx.GetCamera()->WindowToViewport(rect);
    gluPickMatrix(rect.X(), rect.Y(), rect.Width(), rect.Height(),
                  (Int_t*) rnrCtx.GetCamera()->RefViewport().CArr());
  }

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  GLboolean lightp;
  glGetBooleanv(GL_LIGHTING, &lightp);
  if (lightp) glDisable(GL_LIGHTING);

  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDisable(GL_CULL_FACE);
  glEnable(GL_BLEND);
  glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
  RenderMenu();
  RenderPalette(fPaletteLength, 1.6*fWWidth, fWHeight*0.6);
  glPopMatrix();
  glPopAttrib();

  if (lightp) glEnable(GL_LIGHTING);

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);
  RenderCellIDs();
}


/******************************************************************************/
// Protected sub-renderers
/******************************************************************************/

//______________________________________________________________________________
Float_t AliEveITSModuleStepper::TextLength(const char* txt)
{
  // Calculate length of text txt.

  Float_t llx, lly, llz, urx, ury, urz;
  fText->BBox(txt, llx, lly, llz, urx, ury, urz);
  return (urx-llx)*fTextSize;
}

//______________________________________________________________________________
void AliEveITSModuleStepper::RenderString(TString string, Int_t id)
{
  // Render text for button id.

  Float_t txtY = fWHeight*0.5;
  Float_t txtl = TextLength(string.Data());

  if(id > 0) glLoadName(id);
  if(id>0 && fWActive == id)
    fText->SetTextColor(fWActiveCol);
  else
    fText->SetTextColor(fFontCol);


  if(id>0)
  {
    if(fWActive == id)
      fText->SetTextColor(fWActiveCol);
    else
      fText->SetTextColor(fFontCol);

    glLoadName(id);
    Float_t ss = fWWidth*0.4;
    fText->PaintGLText(ss, txtY, -0.8, string.Data());
    // box
    Float_t bw =2*ss+txtl;
    RenderFrame(bw,fWHeight*2,id);
    glTranslatef( bw, 0, 0);
  }
  else
  {
    fText->SetTextColor(fFontCol);
    fText->PaintGLText(0, txtY, -0.8, string.Data());
    glTranslatef(txtl, 0, 0);
  }
}

//______________________________________________________________________________
void AliEveITSModuleStepper::RenderFrame(Float_t dx, Float_t dy, Int_t id)
{
  // Render frame for button id, taking into account if it is currently
  // below mouse.

  if (fRnrFrame == kFALSE)return;

  glPushAttrib(GL_ENABLE_BIT | GL_POLYGON_BIT);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  UChar_t color[4];
  if (fWActive == id)
    TEveUtil::TEveUtil::ColorFromIdx(fWActiveCol, color);
  else
    TEveUtil:: TEveUtil::ColorFromIdx(fWCol, color);
  glColor4ubv(color);

  glBegin(GL_QUADS);
  glVertex2f(0, 0);   glVertex2f(dx, 0);
  glVertex2f(dx, dy); glVertex2f(0, dy);
  glEnd();
  glPopAttrib();
}

//______________________________________________________________________________
void AliEveITSModuleStepper::RenderSymbol(Float_t dx, Float_t dy, Int_t id)
{
  // Render an overlay / GUI symbol, based on button id:
  // 1 ~ <, 2 ~ <<, 3 ~ >, 4 ~ >>, 5 ~ ^, 6 ~ v.

  glLoadName(id);

  UChar_t color[4];
  if (fWActive == id)
    TEveUtil::TEveUtil::ColorFromIdx(fWActiveCol, color);
  else
    TEveUtil::TEveUtil::ColorFromIdx(fWCol, color);
  glColor4ubv(color);

  Float_t xs = dx/4, ys = dy/4;
  if(id == 0) {
    glBegin(GL_QUADS);
    glVertex2f(0,ys); glVertex2f(0, ys*3);
    glVertex2f(dx, ys*3); glVertex2f(dx, ys);
    glEnd();
    return;
  }

  glBegin(GL_TRIANGLES);
  switch (id) {
    case 1:
    {
      // left
      //      glVertex2f(xs*2.5, ys*3); glVertex2f(xs*1.5, ys*2); glVertex2f(xs*2.5, ys);
      glVertex2f(xs*3, ys*3); glVertex2f(xs*1, ys*2); glVertex2f(xs*3, ys);
      break;
    }
    case 2:
    {
      //double left
      glVertex2f(xs*2, ys*3); glVertex2f(xs, ys*2);    glVertex2f(xs*2, ys);
      glVertex2f(xs*3, ys*3); glVertex2f(xs*2, ys*2);  glVertex2f(xs*3, ys);
      break;
    }
    case 3:
    {
      // right
      //glVertex2f(xs*1.5, ys); glVertex2f(xs*2.5, ys*2); glVertex2f(xs*1.5, ys*3);
      glVertex2f(xs*1, ys); glVertex2f(xs*3, ys*2); glVertex2f(xs*1, ys*3);
      break;
    }
    case 4:
    {
      // double right
      glVertex2f(xs, ys);     glVertex2f(xs*2, ys*2);   glVertex2f(xs, ys*3);
      glVertex2f(xs*2, ys);   glVertex2f(xs*3, ys*2);   glVertex2f(xs*2, ys*3);
      break;
    }
    case 5:
    {
      // up
      glVertex2f(xs, ys*2.5);  glVertex2f(xs*2, ys*3.5); glVertex2f(xs*3, ys*2.5);
      break;
    }
    case 6:
    {
      // down
      glVertex2f(xs, ys*1.5);  glVertex2f(xs*2, ys*0.5); glVertex2f(xs*3, ys*1.5);
      break;
    }

    default:
      break;
  }
  glEnd();
  glLoadName(0);
}

//______________________________________________________________________________
void AliEveITSModuleStepper::RenderPalette(Float_t dx, Float_t x, Float_t y)
{
  // Render color palette with number axis.

  glPushMatrix();
  glLoadIdentity();
  glTranslatef(1 -x- dx, -1+y*4, 0);
  AliEveITSModule* qs = dynamic_cast<AliEveITSModule*>(*BeginChildren());
  TEveRGBAPalette* p = qs->GetPalette();
  glBegin(GL_QUAD_STRIP);
  glColor4ubv(p->ColorFromValue(p->GetMinVal()));
  glVertex2f(0, 0);
  glVertex2f(0, y);
  if (p->GetMaxVal() > p->GetMinVal() + 1)
  {
    Float_t xs = dx/(p->GetMaxVal() - p->GetMinVal());
    Float_t x0 = xs;
    for(Int_t i=p->GetMinVal() + 1; i<p->GetMaxVal(); i++)
    {
      glColor4ubv(p->ColorFromValue(i));
      glVertex2f(x0, 0);
      glVertex2f(x0, y);
      x0+=xs;
    }
  }
  glColor4ubv(p->ColorFromValue(p->GetMaxVal()));
  glVertex2f(dx, 0);
  glVertex2f(dx, y);
  glEnd();

  if (p->GetMaxVal() > p->GetMinVal())
  {
    glRotatef(-90,1, 0, 0 );
    Double_t v1[3] = {0., 0., 0.};
    Double_t v2[3] = {dx, 0, 0.};
    fAxis->SetLabelsSize(fTextSize/dx);
    fAxis->PaintGLAxis(v1, v2, p->GetMinVal(), p->GetMaxVal(), 206);
  }
  glPopMatrix();
}

//______________________________________________________________________________
void AliEveITSModuleStepper::RenderMenu()
{
  // Render menu: page control, scale control, detector type buttons.

  Float_t ww = 2*fWWidth;
  Float_t wh = 2*fWHeight;

  // transparent bar
  Float_t a=0.3;
  glColor4f(a, a, a, a);
  Float_t H = 1.9*wh*(1+ 2*fWOff);
  if(1) {
    glBegin(GL_QUADS);
    glVertex3f(-1, -1,   0.1); glVertex3f(-1, -1+H, 0.1);
    glVertex3f(1 , -1+H, 0.1); glVertex3f( 1, -1  , 0.1);
    glEnd();
  }

  Float_t y_base = -1 + wh*0.35;
  glTranslatef(-1, y_base, 0.);
  glPushName(0);
  // pager
  glPushMatrix();
  glTranslatef(ww, 0, 0.);
  fText->SetTextSize(fTextSize);
  Float_t soff = ww*1.3;
  glTranslatef(0, fWOff*wh, 0);
  RenderSymbol(ww, wh, 2);
  RenderFrame(ww,wh,2);
  glTranslatef(soff, 0, 0);
  RenderSymbol(ww, wh, 1);
  RenderFrame(ww,wh,1);
  glTranslatef(soff, 0, 0);
  // text info
  {
    const char* txt =  Form("%d/%d ", GetCurrentPage(), GetPages());
    Float_t dx = (fPagerGap - TextLength(txt))*0.5;
    fText->SetTextColor(fFontCol);
    fText->PaintGLText(dx, wh*0.25, -0.8, txt);
  }
  glTranslatef(fPagerGap, 0, 0);

  RenderSymbol(ww, wh, 3);
  RenderFrame(ww,wh,3);
  glTranslatef(soff, 0, 0);
  RenderSymbol(ww, wh, 4);
  RenderFrame(ww,wh,4);
  glTranslatef(2*ww, 0, 0);
  glPopMatrix();

  // scale info
  glPushMatrix();
  AliEveITSDigitsInfo* di = fDigitsInfo;
  Int_t scale = fScaleInfo->GetScale() - 1;
  AliEveITSScaledModule* sm = dynamic_cast<AliEveITSScaledModule*>(*BeginChildren());
  Int_t cnx = 0, cnz = 0;
  switch(sm->GetSubDetID())
  {
    case 0:
      cnx = di->fSPDScaleX[scale], cnz = di->fSPDScaleZ[scale];
      break;
    case 1:
      cnx = di->fSDDScaleX[scale], cnz = di->fSDDScaleZ[scale];
      break;
    case 2:
      cnx = di->fSSDScale[scale], cnz = 1;
      break;
  }
  glTranslatef(10*ww,0, 0);
  RenderString(Form("Zoom: "));
  glPushMatrix();
  glTranslatef(0, 0.2*wh, 0);
  RenderSymbol(ww, wh*0.9, 5);
  glTranslatef(0, 0.4*wh, 0);
  RenderFrame(ww, wh*0.5, 5);
  glPopMatrix();
  RenderSymbol(ww, wh*0.9, 6);
  RenderFrame(ww, wh*0.5, 6);
  glTranslatef(ww, 0, 0);
  RenderString(Form("%dx%d ", cnx, cnz));
  glPopMatrix();

  //choose detector
  glPushMatrix();
  glTranslatef(18*ww, 0, 0);
  Float_t bs = ww*0.2;
  RenderString("SPD", 8);
  glTranslatef(bs, 0, 0);
  RenderString("SDD", 9);
  glTranslatef(bs, 0, 0);
  RenderString("SSD", 10);
  glPopMatrix();

  glPopName();
}

//______________________________________________________________________________
void AliEveITSModuleStepper::RenderCellIDs()
{
  // Render module-ids under their cells.

  fText->SetTextSize(fStepper->GetDy()*0.1);
  fText->SetTextColor(fFontCol);
  Double_t x, y, z;
  Double_t sx, sy, sz;
  UInt_t idx = fPosition;
  for (List_i childit=fChildren.begin(); childit!=fChildren.end(); ++childit)
  {
    if(idx < fIDs.size())
    {
      AliEveITSScaledModule* mod = dynamic_cast<AliEveITSScaledModule*>(*childit);
      TEveTrans& tr = mod->RefHMTrans();
      TString name = Form("%d",mod->GetID());
      tr.GetPos(x,y,z);
      x += fStepper->GetDx()*0.5;
      y -= fStepper->GetDy()*0.5;
      z += 0.4; // !!! MT hack - cross check with overlay rendering.
      Float_t llx, lly, llz, urx, ury, urz;
      fText->BBox(name, llx, lly, llz, urx, ury, urz);
      tr.GetScale(sx, sy, sz);
      fText->PaintGLText(x-(urx-llx)*sx, y, z, name);
      idx++;
    }
  }
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
      if (fWActive != item) {
        fWActive = item;
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
          if(si->GetScale() < 5)
          {
            si->ScaleChanged(si->GetScale() + 1);
            ElementChanged(kTRUE, kTRUE);
          }
          break;
        }
        case 6:
        {
          AliEveDigitScaleInfo* si = fScaleInfo;
          if(si->GetScale() > 1)
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

  fWActive = -1;
}
