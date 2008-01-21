// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_TPCSector2DGL_H
#define ALIEVE_TPCSector2DGL_H

#include <TGLObject.h>

#include <EveDet/AliEveTPCSector2D.h>
#include <EveDet/AliEveTPCSectorData.h>

class TGLViewer;
class TGLScene;


class AliEveTPCSector2DGL : public TGLObject
{
  AliEveTPCSector2DGL(const AliEveTPCSector2DGL&);            // Not implemented
  AliEveTPCSector2DGL& operator=(const AliEveTPCSector2DGL&); // Not implemented

protected:
  virtual void DirectDraw(TGLRnrCtx & rnrCtx) const;

  void LoadPadrow(AliEveTPCSectorData::RowIterator& iter, Int_t row, Int_t off) const;
  void CreateTexture() const;

  void DisplayTexture(const AliEveTPCSectorData::SegmentInfo& seg,
                      Int_t startCol, Int_t startRow) const;
  void DisplayQuads(const AliEveTPCSectorData::SegmentInfo& seg,
		    Int_t startCol, Int_t startRow) const;
  void DisplayNamedQuads(const AliEveTPCSectorData::SegmentInfo& seg,
			 Int_t startCol, Int_t startRow) const;
  void DisplayFrame() const;

  UChar_t* GetRowCol(Int_t row, Int_t col) const;

  AliEveTPCSector2D*                 fSector;
  mutable AliEveTPCSectorData*       fSectorData;

  mutable UChar_t*             fImage;
  mutable UInt_t               fTexture;
  mutable UInt_t               fRTS;

public:
  AliEveTPCSector2DGL();
  virtual ~AliEveTPCSector2DGL();

  virtual Bool_t SetModel(TObject* obj, const Option_t* opt=0);
  virtual void   SetBBox();
  virtual Bool_t SupportsSecondarySelect() const { return kTRUE; }

  virtual void ProcessSelection(TGLRnrCtx & rnrCtx, TGLSelectRecord & rec);

  static void TraceStepsUp  (const AliEveTPCSectorData::SegmentInfo& s);
  static void TraceStepsDown(const AliEveTPCSectorData::SegmentInfo& s);

  static const Int_t fgkTextureWidth;
  static const Int_t fgkTextureHeight;
  static const Int_t fgkTextureByteSize;

}; // endclass TPCSector2D_GL_Rnr


inline UChar_t* AliEveTPCSector2DGL::GetRowCol(Int_t row, Int_t col) const
{
  return fImage + 4*(row*fgkTextureWidth + col);
}


#endif
