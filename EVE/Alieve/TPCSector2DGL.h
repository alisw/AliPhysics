// $Header$

#ifndef ALIEVE_TPCSector2DGL_H
#define ALIEVE_TPCSector2DGL_H

#include <TGLObject.h>

#include <Alieve/TPCSector2D.h>
#include <Alieve/TPCSectorData.h>

class TGLViewer;
class TGLScene;

namespace Alieve {

class TPCSector2DGL : public TGLObject
{
  TPCSector2DGL(const TPCSector2DGL&);            // Not implemented
  TPCSector2DGL& operator=(const TPCSector2DGL&); // Not implemented

protected:
  virtual void DirectDraw(const TGLDrawFlags & flags) const;

  void LoadPadrow(TPCSectorData::RowIterator& iter, Int_t row, Int_t off) const;
  void CreateTexture() const;

  void DisplayTexture(const TPCSectorData::SegmentInfo& seg,
                      Int_t startCol, Int_t startRow) const;
  void DisplayQuads(const TPCSectorData::SegmentInfo& seg,
		    Int_t startCol, Int_t startRow) const;
  void DisplayNamedQuads(const TPCSectorData::SegmentInfo& seg,
			 Int_t startCol, Int_t startRow) const;
  void DisplayFrame() const;

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
  virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }

  virtual void ProcessSelection(UInt_t* ptr, TGLViewer*, TGLScene*);

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
