// $Header$

#include "NLTProjectorGL.h"
#include <Reve/NLTProjector.h>

#include <TGLRnrCtx.h>
#include <TGLIncludes.h>
#include <TGLText.h>
#include <TMath.h>

#include <list>

using namespace Reve;

//______________________________________________________________________
// NLTProjectorGL
//

ClassImp(NLTProjectorGL)

NLTProjectorGL::NLTProjectorGL() :
  TGLObject(),

  fRange(300),
  fLabelSize(0.02),
  fLabelOff(0.018),
  fTMSize(0.02),

  fM(0),
  fText(0)
{
  fDLCache = kFALSE; // Disable display list.
  fText = new TGLText();
  fText->SetGLTextFont(40);
  fText->SetTextColor(0);
}

NLTProjectorGL::~NLTProjectorGL()
{}
/**************************************************************************/
const char* NLTProjectorGL::GetText(Float_t x) const
{
  Float_t v = 10*TMath::Nint(x/10.0f);
  return Form("%.0f", v);
}

/**************************************************************************/

void NLTProjectorGL::DrawTickMarks(Float_t tm) const
{
  glBegin(GL_LINES);
  for( std::list<Float_t>::iterator pi = fPos.begin(); pi!= fPos.end(); pi++)
  {
    glVertex3f(*pi, 0,   0.);
    glVertex3f(*pi, tm, 0.);
  }
  glEnd();
}

/**************************************************************************/

void NLTProjectorGL::DrawHInfo() const
{
  Float_t tms = fTMSize*fRange;
  DrawTickMarks(-tms);

  glPushMatrix();
  glRotatef(-90, 1, 0, 0);
  glTranslatef(0, 0, -tms -fLabelOff*fRange);
  const char* txt;
  Float_t llx, lly, llz, urx, ury, urz;
  std::list<Float_t>::iterator vi = fVals.begin();
  for( std::list<Float_t>::iterator pi = fPos.begin(); pi!= fPos.end(); pi++)
  {
    txt = GetText(*vi);
    fText->BBox(txt, llx, lly, llz, urx, ury, urz);
    fText->PaintGLText(*pi -(urx-llx)*fText->GetTextSize()*0.5, 0, 0, txt);
    vi++;
  }
  glPopMatrix();

  fPos.clear(); fVals.clear();
}

/**************************************************************************/
void NLTProjectorGL::DrawVInfo() const
{
  Float_t tms = fTMSize*fRange;
  glRotatef(90, 0, 0, 1);
  DrawTickMarks(tms);
  glRotatef(-90, 0, 0, 1);

  glPushMatrix();
  glRotatef(-90, 1, 0, 0);
  glTranslatef(-fLabelOff*fRange -tms, 0, 0);
  const char* txt;
  Float_t llx, lly, llz, urx, ury, urz;
  std::list<Float_t>::iterator vi = fVals.begin();
  for( std::list<Float_t>::iterator pi = fPos.begin(); pi!= fPos.end(); pi++)
  {
    txt= GetText(*vi);
    fText->BBox(txt, llx, lly, llz, urx, ury, urz);
    fText->PaintGLText( -tms, 0, *pi - (ury - lly)*fText->GetTextSize()*0.5, txt);
    vi++;
  }
  glPopMatrix();

  fPos.clear(); fVals.clear();
}

/**************************************************************************/
void NLTProjectorGL::SplitInterval(Int_t ax) const
{
  if (fM->GetSplitInfoLevel())
  {
    if(fM->GetSplitInfoMode())
      SplitIntervalByVal(fVals.front(), fVals.back(), ax, 0);
    else
      SplitIntervalByPos(fPos.front(), fPos.back(), ax, 0);
  }
}

/**************************************************************************/
void NLTProjectorGL::SplitIntervalByPos(Float_t minp, Float_t maxp, Int_t ax, Int_t level) const
{
  Float_t p = (minp+maxp)*0.5;
  fPos.push_back(p);
  Float_t v = fM->GetProjection()->GetValForScreenPos(ax, p);
  fVals.push_back(v);
  // printf("level %d position %f value %f\n", level, p,v);
  level++;
  if(level<fM->GetSplitInfoLevel())
  {
    SplitIntervalByPos(minp, p , ax, level);
    SplitIntervalByPos(p, maxp, ax, level);
  }
}

/**************************************************************************/
void NLTProjectorGL::SplitIntervalByVal(Float_t minv, Float_t maxv, Int_t ax, Int_t level) const
{
  Float_t v = (minv+maxv)*0.5;
  fVals.push_back(v);
  Float_t p = fM->GetProjection()->GetScreenVal(ax, v);
  fPos.push_back(p);
  //printf("level %d position %f value %f MINMAX val(%f, %f)\n", level, p,v, minv, maxv);
  level++;
  if(level<fM->GetSplitInfoLevel())
  {
    SplitIntervalByVal(minv, v , ax, level);
    SplitIntervalByVal(v, maxv, ax, level);
  }
}


/**************************************************************************/
void NLTProjectorGL::DirectDraw(TGLRnrCtx & /*rnrCtx*/) const
{
  // printf("NLTProjectorGL::DirectDraw %d\n.", fM->GetMainColor());
  GLboolean lightp;
  glGetBooleanv(GL_LIGHTING, &lightp);
  if (lightp) glDisable(GL_LIGHTING);

  Float_t* bbox = fM->GetBBox();
  fRange = bbox[1] - bbox[0];

  // value of projected bbox
  Float_t bbv[4];
  bbv[0] = fM->GetProjection()->GetValForScreenPos(0,bbox[0]);
  bbv[1] = fM->GetProjection()->GetValForScreenPos(0,bbox[1]);
  bbv[2] = fM->GetProjection()->GetValForScreenPos(1,bbox[2]);
  bbv[3] = fM->GetProjection()->GetValForScreenPos(1,bbox[3]);

  Vector zeroPos;
  fM->GetProjection()->ProjectVector(zeroPos);
  fText->SetTextSize(fLabelSize*fRange);
  fText->SetTextColor(fM->GetAxisColor());

  { // horizontal
    glPushMatrix();
    glTranslatef(0, bbox[2], 0);
    // left
    fPos.push_back(bbox[0]); fPos.push_back(zeroPos.x);
    fVals.push_back(bbv[0]); fVals.push_back(0);
    SplitInterval(0);
    DrawHInfo();
    // right
    fPos.push_back(zeroPos.x); fPos.push_back(bbox[1]);
    fVals.push_back(0); fVals.push_back(bbv[1]);
    SplitInterval(0); fVals.pop_front(); fPos.pop_front();
    DrawHInfo();
    glPopMatrix();
  }
  { // vertical
    glPushMatrix();
    glTranslatef(bbox[0], 0, 0);
    // bottom
    fPos.push_back(bbox[2]); fPos.push_back(zeroPos.y);
    fVals.push_back(bbv[2]); fVals.push_back(0);
    SplitInterval(1);
    DrawVInfo();
    // top
    fPos.push_back(zeroPos.y); fPos.push_back(bbox[3]);
    fVals.push_back(0); fVals.push_back(bbv[3]);
    SplitInterval(1);fPos.pop_front(); fVals.pop_front();
    DrawVInfo();
    glPopMatrix();
  }

  // body
  glBegin(GL_LINES);
  glVertex3f(bbox[0], bbox[2], 0.);
  glVertex3f(bbox[1], bbox[2], 0.);
  glVertex3f(bbox[0], bbox[2], 0.);
  glVertex3f(bbox[0], bbox[3], 0.);
  glEnd();

  Float_t d = 10;
  if(fM->GetDrawCenter())
  {
    Float_t* c = fM->GetProjection()->GetProjectedCenter();
    glColor3f(1., 0., 0.);
    glBegin(GL_LINES);
    glVertex3f(c[0] +d, c[1],    c[2]);     glVertex3f(c[0] - d, c[1]   , c[2]);
    glVertex3f(c[0] ,   c[1] +d, c[2]);     glVertex3f(c[0]    , c[1] -d, c[2]);
    glVertex3f(c[0] ,   c[1],    c[2] + d); glVertex3f(c[0]    , c[1]   , c[2] - d);
    glEnd();

  }

  if(fM->GetDrawOrigin())
  {
      Vector zero;
      fM->GetProjection()->ProjectVector(zero);
      glColor3f(1., 1., 1.);
      glBegin(GL_LINES);
      glVertex3f(zero[0] +d, zero[1],    zero[2]);     glVertex3f(zero[0] - d, zero[1]   , zero[2]);
      glVertex3f(zero[0] ,   zero[1] +d, zero[2]);     glVertex3f(zero[0]    , zero[1] -d, zero[2]);
      glVertex3f(zero[0] ,   zero[1],    zero[2] + d); glVertex3f(zero[0]    , zero[1]   , zero[2] - d);
      glEnd();
  }
  if (lightp) glEnable(GL_LIGHTING);
}


/**************************************************************************/

Bool_t NLTProjectorGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  if(SetModelCheckClass(obj, NLTProjector::Class())) {
    fM = dynamic_cast<NLTProjector*>(obj);
    return kTRUE;
  }
  return kFALSE;
}

/**************************************************************************/

void NLTProjectorGL::SetBBox()
{
  // !! This ok if master sub-classed from TAttBBox
  SetAxisAlignedBBox(((NLTProjector*)fExternalObj)->AssertBBox());
}
