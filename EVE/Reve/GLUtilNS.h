#ifndef Reve_GLUtilNS_H
#define Reve_GLUtilNS_H
#ifndef __CINT__

#include <TGLIncludes.h>
#include <TObject.h>

class TAttMarker;

namespace GLUtilNS
{

class GL_Capability_Switch {
  GLenum    fWhat;
  GLboolean fState;
  bool      fFlip;

  void set_state(GLboolean s)
  { if(s) glEnable(fWhat); else glDisable(fWhat); }

public:
  GL_Capability_Switch(GLenum what, GLboolean state) :
    fWhat(what), fState(kFALSE), fFlip(kFALSE)
  {
    fState = glIsEnabled(fWhat);
    fFlip  = (fState != state);
    if(fFlip)	set_state(state);
  }
  ~GL_Capability_Switch()
  { if(fFlip) set_state(fState); }
};

class GL_Float_Holder
{
  GL_Float_Holder(const GL_Float_Holder&);            // Not implemented
  GL_Float_Holder& operator=(const GL_Float_Holder&); // Not implemented

  GLenum    fWhat;
  GLfloat   fState;
  bool      fFlip;
  void    (*fFoo)(GLfloat);

public:
  GL_Float_Holder(GLenum what, GLfloat state, void (*foo)(GLfloat)) :
    fWhat(what), fState(kFALSE), fFlip(kFALSE), fFoo(foo)
  {
    glGetFloatv(fWhat, &fState);
    fFlip = (fState != state);
    if(fFlip) fFoo(state);
  }
  ~GL_Float_Holder()
  { if(fFlip) fFoo(fState); }
};


void RenderPolyMarkers(TAttMarker& marker, Float_t* p, Int_t n, 
		       Bool_t selection, Bool_t sec_selection);

void RenderPoints(TAttMarker& marker, Float_t* p, Int_t n,
		  Bool_t selection, Bool_t sec_selection);

void RenderCrosses(TAttMarker& marker, Float_t* p, Int_t n, Bool_t sec_selection);

}

#endif
#endif
