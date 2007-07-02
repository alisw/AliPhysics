#include "NLTPolygonSetGL.h"
#include "NLTPolygonSet.h"
#include "PODs.h"

#include <TGLRnrCtx.h>
#include <TGLIncludes.h>


using namespace Reve;

/**************************************************************************/

NLTPolygonSetGL::NLTPolygonSetGL() : TGLObject()
{
  // fDLCache = false; // Disable DL.
}

NLTPolygonSetGL::~NLTPolygonSetGL()
{}

/**************************************************************************/
Bool_t NLTPolygonSetGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  return SetModelCheckClass(obj, NLTPolygonSet::Class());
}

/**************************************************************************/

void NLTPolygonSetGL::SetBBox()
{
  SetAxisAlignedBBox(((NLTPolygonSet*)fExternalObj)->AssertBBox());
}

/**************************************************************************/
static GLUtriangulatorObj *GetTesselator()
{
   static struct Init {
      Init()
      {
#if defined(R__WIN32)
         typedef void (CALLBACK *tessfuncptr_t)();
#elif defined(R__AIXGCC)
         typedef void (*tessfuncptr_t)(...);
#else
         typedef void (*tessfuncptr_t)();
#endif
         fTess = gluNewTess();

         if (!fTess) {
            Error("GetTesselator::Init", "could not create tesselation object");
         } else {
            gluTessCallback(fTess, (GLenum)GLU_BEGIN, (tessfuncptr_t)glBegin);
            gluTessCallback(fTess, (GLenum)GLU_END, (tessfuncptr_t)glEnd);
            gluTessCallback(fTess, (GLenum)GLU_VERTEX, (tessfuncptr_t)glVertex3fv);
         }
      }
      ~Init()
      {
         if(fTess)
            gluDeleteTess(fTess);
      }
      GLUtriangulatorObj *fTess;
   }singleton;

   return singleton.fTess;
}

/**************************************************************************/
void NLTPolygonSetGL::DirectDraw(TGLRnrCtx & /*rnrCtx*/) const
{
  //printf("NLTPolygonSetGL::DirectDraw %s \n",fExternalObj->GetName() );
  NLTPolygonSet& PS = * (NLTPolygonSet*) fExternalObj;
  if(PS.fNPols == 0) return;

  glPushAttrib(GL_ENABLE_BIT | GL_LINE_BIT | GL_POLYGON_BIT);

  glDisable(GL_LIGHTING);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  glPolygonMode(GL_FRONT, GL_FILL);
  glPolygonMode(GL_BACK,  GL_FILL);
  glDisable(GL_CULL_FACE);

  // polygons
  UChar_t col[4];
  ColorFromIdx(PS.fFillColor, col);
  glColor4ubv(col);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.,1.);
  GLUtriangulatorObj *tessObj = GetTesselator();
  for(Int_t i = 0; i < PS.fNPols; i++ ) 
  {
    Int_t ppi;
    {
      if(PS.fPols[i].fNPnts < 4) 
      {
	glBegin(GL_POLYGON);
	Int_t ppi;
	for(Int_t k=0; k<PS.fPols[i].fNPnts; k++)
	{
	  ppi = PS.fPols[i].fPnts[k]; 
	  glVertex3fv(PS.fPnts[ppi].c_vec());
	}
	glEnd();
      }
      else {
	gluBeginPolygon(tessObj);
	gluNextContour(tessObj, (GLenum)GLU_UNKNOWN);
	glNormal3f(0., 0., 1.);
	Double_t coords[3];
	coords[2] = 0.;
	for (Int_t k = 0; k <PS.fPols[i].fNPnts; k++)
	{
	  ppi = PS.fPols[i].fPnts[k];
	  coords[0] = PS.fPnts[ppi].x;
	  coords[1] = PS.fPnts[ppi].y;
	  gluTessVertex(tessObj, coords, PS.fPnts[ppi].c_vec());
	}
	gluEndPolygon(tessObj);
      }
    }
  }
  glDisable(GL_POLYGON_OFFSET_FILL);

  // outline 
  UChar_t lcol[4];
  ColorFromIdx(PS.fLineColor, lcol);
  glColor4ubv(lcol);
  glEnable(GL_LINE_SMOOTH);

  glLineWidth(PS.fLineWidth);
  Int_t ppi;
  for(Int_t i = 0; i < PS.fNPols; i++ ) 
  {
    glBegin(GL_LINE_LOOP); 
    for(Int_t k=0; k<PS.fPols[i].fNPnts; k++)
    {
      ppi = PS.fPols[i].fPnts[k];
      glVertex3fv(PS.fPnts[ppi].c_vec());
    }
    glEnd();
  }

  glPopAttrib();
}

