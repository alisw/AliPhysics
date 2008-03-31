// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTPCSector2DGL_H
#define AliEveTPCSector2DGL_H

#include <TGLObject.h>
class TGLViewer;
class TGLScene;

#include <EveDet/AliEveTPCSectorData.h>
class AliEveTPCSector2D;

//------------------------------------------------------------------------------
// AliEveTPCSector2DGL
//
// GL renderer for AliEveTPCSector2D.
//

class AliEveTPCSector2DGL : public TGLObject
{
public:
  AliEveTPCSector2DGL();
  virtual ~AliEveTPCSector2DGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();
  virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }

  virtual void ProcessSelection(TGLRnrCtx & rnrCtx, TGLSelectRecord & rec);

  static void TraceStepsUp  (const AliEveTPCSectorData::SegmentInfo& s);
  static void TraceStepsDown(const AliEveTPCSectorData::SegmentInfo& s);

protected:
  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

  void LoadPadrow(AliEveTPCSectorData::RowIterator& iter, Int_t row, Int_t colOff) const;
  void CreateTexture() const;

  void DisplayTexture(const AliEveTPCSectorData::SegmentInfo& seg,
                      Int_t startCol, Int_t startRow) const;
  void DisplayQuads(const AliEveTPCSectorData::SegmentInfo& seg,
		    Int_t startCol, Int_t startRow) const;
  void DisplayNamedQuads(const AliEveTPCSectorData::SegmentInfo& seg,
			 Int_t startCol, Int_t startRow) const;
  void DisplayFrame() const;

  UChar_t* GetRowCol(Int_t row, Int_t col) const;

  AliEveTPCSector2D              *fSector;     // Model object.
  mutable AliEveTPCSectorData    *fSectorData; // Data of displayed sector.

  mutable UChar_t                *fImage;      // Texture data.
  mutable UInt_t                  fTexture;    // Texture id.
  mutable UInt_t                  fRTS;        // Render time-stamp.

  static const Int_t fgkTextureWidth;     // Width of the texture.
  static const Int_t fgkTextureHeight;    // Height of the texture.
  static const Int_t fgkTextureByteSize;  // Size of the texture in bytes.

private:
  AliEveTPCSector2DGL(const AliEveTPCSector2DGL&);            // Not implemented
  AliEveTPCSector2DGL& operator=(const AliEveTPCSector2DGL&); // Not implemented

  ClassDef(AliEveTPCSector2DGL, 0); // GL renderer for AliEveTPCSector2D.
};


inline UChar_t* AliEveTPCSector2DGL::GetRowCol(Int_t row, Int_t col) const
{
  return fImage + 4*(row*fgkTextureWidth + col);
}


#endif
