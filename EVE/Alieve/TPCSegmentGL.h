// $Header$

#ifndef ALIEVE_TPCSegmentGL_H
#define ALIEVE_TPCSegmentGL_H


#include <TGLObject.h>

#include <Alieve/TPCSegment.h>

#include <GL/gl.h>
#include <GL/glu.h>

namespace Alieve {

  const int ImageWidth  = 256;
  const int ImageHeight = 128;

  class TPCSegmentGL : public TGLObject
  {
  protected:
    virtual       void DirectDraw(const TGLDrawFlags & flags) const;

  private:
    void          LoadPadrow(Int_t row, Int_t off) const;
    void          DisplayTexture(Float_t pw, Float_t pl, Float_t vR, Int_t nMaxPads, 
				  Int_t nRows, Int_t startRow,Int_t startCol) const;
    void          DisplayQuads(Float_t pw, Float_t pl, Float_t vR, Int_t nMaxPads, 
				Int_t nRows, Int_t startRow,Int_t startCol) const;
    void          DisplayFrame(TPCDigitsInfo* info) const;

    GLubyte* GetRow(Int_t row) const;
    GLubyte* GetRowCol(Int_t row, Int_t col) const;

    // protected:
    TPCSegment*	                 fSegment;
    mutable UChar_t*             fImage;
    mutable UInt_t               fTexture;
    mutable UInt_t               fRTS;
 
  public:
    TPCSegmentGL();
    virtual       ~TPCSegmentGL();
    virtual        Bool_t SetModel(TObject* obj); 
    virtual void   SetBBox();

    void           SetCol(Float_t z, UChar_t* pixel) const;
    void           CreateTexture(TPCDigitsInfo* info) const;

  }; // endclass TPCSegment_GL_Rnr
  

  inline  UChar_t* TPCSegmentGL::GetRowCol(Int_t row, Int_t col) const
  {
    if ( row > ImageHeight) printf("ERROR row %d ImageHeight %d\n", row, col);
    return fImage + (ImageWidth*row +col)*4*sizeof(UChar_t); //*sizeof();
  }

  inline void LoopStepsUp(Alieve::TPCSeg* seg) 
  {
    Float_t x = -(seg->fNMaxPads*1.0/2 - seg->fNsteps)*seg->fPadWidth;
    Float_t y  = seg->fRlow;
    glVertex3f(x,y,0.);
    for(int s = 0; s <seg->fNsteps ;s++){
      y = seg->fStepY[s];
      glVertex3f(x,y,0.);
      x -= seg->fPadWidth;
      glVertex3f(x,y,0.);
    }
    y =  seg->fRlow + seg->fNRows*seg->fPadLength;
    glVertex3f(-seg->fNMaxPads*seg->fPadWidth/2,y,0.);
  }

  inline void LoopStepsDown(TPCSeg* seg) 
  {
    Float_t x = seg->fNMaxPads*seg->fPadWidth/2;
    Float_t y = seg->fRlow + seg->fNRows*seg->fPadLength;
    glVertex3f(x,y,0.);
    for (int s = (seg->fNsteps -1); s >= 0 ;s--){
      y =  seg->fStepY[s];
      glVertex3f(x,y,0.);
      x -= seg->fPadWidth;
      glVertex3f(x,y,0.);
    }
    y = seg->fRlow;
    glVertex3f((seg->fNMaxPads*1.0/2 - seg->fNsteps)*seg->fPadWidth,y, 0.);
  }
}

#endif
