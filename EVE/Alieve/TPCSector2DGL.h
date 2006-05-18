// $Header$

#ifndef ALIEVE_TPCSector2DGL_H
#define ALIEVE_TPCSector2DGL_H

#include <TGLObject.h>

#include <Alieve/TPCSector2D.h>
#include <Alieve/TPCSectorData.h>


namespace Alieve {

class TPCSector2DGL : public TGLObject
{
protected:
  virtual void DirectDraw(const TGLDrawFlags & flags) const;

  void LoadPadrow(TPCSectorData::RowIterator& iter, Int_t row, Int_t off) const;
  void CreateTexture() const;

  void DisplayTexture(Float_t padW,     Float_t padH, Float_t startR,
                      Int_t numMaxPads, Int_t numRows,
                      Int_t startCol,   Int_t startRow) const;
  void DisplayQuads  (Float_t padW,     Float_t padH, Float_t startR,
		      Int_t numMaxPads, Int_t numRows,
		      Int_t startCol,   Int_t startRow) const;
  void DisplayFrame  () const;

  UChar_t* GetRowCol(Int_t row, Int_t col) const;

  TPCSector2D*                 fSector;
  mutable TPCSectorData*       fSectorData;

  mutable UChar_t*             fImage;
  mutable UInt_t               fTexture;
  mutable UInt_t               fRTS;
 
public:
  TPCSector2DGL();
  virtual ~TPCSector2DGL();

  virtual Bool_t SetModel(TObject* obj); 
  virtual void   SetBBox();

  void SetCol(Float_t z, UChar_t* pixel) const;

  static void TraceStepsUp  (const TPCSectorData::SegmentInfo& s);
  static void TraceStepsDown(const TPCSectorData::SegmentInfo& s);

  static const Int_t fgkTextureWidth;
  static const Int_t fgkTextureHeight;
  static const Int_t fgkTextureByteSize;

}; // endclass TPCSector2D_GL_Rnr
  

inline UChar_t* TPCSector2DGL::GetRowCol(Int_t row, Int_t col) const
{
  return fImage + 4*(row*fgkTextureWidth + col);
}


}

#endif
