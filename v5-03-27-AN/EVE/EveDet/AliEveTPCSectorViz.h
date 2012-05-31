// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveTPCSectorViz_H
#define AliEveTPCSectorViz_H

#include <TEveElement.h>

#include <TNamed.h>
#include <TAtt3D.h>
#include <TAttBBox.h>

class AliEveTPCData; class AliEveTPCSectorData;

class AliEveTPCSectorVizEditor;
class AliEveTPCSector2D;  class AliEveTPCSector2DEditor;  class AliEveTPCSector2DGL;
class AliEveTPCSector3D;  class AliEveTPCSector3DEditor;  class AliEveTPCSector3DGL;

//------------------------------------------------------------------------------
// AliEveTPCSectorViz
//
// Base-class for visualization of data for one TPC sector.
//

class AliEveTPCSectorViz : public TEveElement,
			   public TNamed,
			   public TAtt3D,
			   public TAttBBox
{
  friend class AliEveTPCSectorVizEditor;
  friend class AliEveTPCSector2D;
  friend class AliEveTPCSector2DEditor;
  friend class AliEveTPCSector2DGL;
  friend class AliEveTPCSector3D;
  friend class AliEveTPCSector3DEditor;
  friend class AliEveTPCSector3DGL;

public:
  AliEveTPCSectorViz(const Text_t* n="AliEveTPCSectorViz", const Text_t* t=0);
  virtual ~AliEveTPCSectorViz();

  virtual void CopyVizParams(const TEveElement* el);

  virtual UInt_t IncRTS() { return ++fRTS; }
  virtual Bool_t CanEditMainColor() const { return kTRUE; }

  void SetDataSource(AliEveTPCData* data);
  void SetSectorID(Int_t id);

  AliEveTPCData* GetData()     const { return fTPCData; }
  Int_t          GetSectorID() const { return fSectorID; }
  AliEveTPCSectorData* GetSectorData() const;

  Int_t GetMinTime() const { return fMinTime; }
  Int_t GetMaxTime() const { return fMaxTime; }
  void SetMinTime(Int_t mt)    { fMinTime   = mt; IncRTS(); }
  void SetMaxTime(Int_t mt)    { fMaxTime   = mt; IncRTS(); }
  void SetThreshold(Short_t t);
  void SetMaxVal(Int_t mv);

  void SetRnrInn(Bool_t r)     { fRnrInn  = r; IncRTS(); }
  void SetRnrOut1(Bool_t r)    { fRnrOut1 = r; IncRTS(); }
  void SetRnrOut2(Bool_t r)    { fRnrOut2 = r; IncRTS(); }

  void SetFrameColor(Color_t col)     { fFrameColor = col; IncRTS(); }
  virtual void SetRnrFrame(Bool_t rf) { fRnrFrame = rf;  IncRTS(); }
  void SetAutoTrans(Bool_t t);

  void SetUseTrans(Bool_t t);

protected:
  AliEveTPCData    *fTPCData;    //  Source of data.
  Int_t             fSectorID;   //  Id of the displayed sector.

  Int_t             fMinTime;    //  Min time-bin to display.
  Int_t             fMaxTime;    //  Max time-bin to display.
  Short_t           fThreshold;  //  Threshold for display/
  Int_t             fMaxVal;     //  Maximum signal-value, all above is of the same color.

  Bool_t            fRnrInn;     //  Render inner segment.
  Bool_t            fRnrOut1;    //  Render middle segment.
  Bool_t            fRnrOut2;    //  Render outer segment.

  Color_t           fFrameColor; //  Color of the frame, the main color.
  Bool_t            fRnrFrame;   //  Render frame.
  Bool_t            fAutoTrans;  //  Automatically calculate transformation based on sector id.
  UInt_t            fRTS;        //! Rendering TimeStamp

  mutable UChar_t  *fColorArray; //  Color array caching signal to color mapping.

  void SetupColor(Int_t val, UChar_t* pix) const;
  void ClearColorArray();
  void SetupColorArray() const;
  UChar_t* ColorFromArray(Int_t val) const;
  void     ColorFromArray(Int_t val, UChar_t* pix) const;

private:
  AliEveTPCSectorViz(const AliEveTPCSectorViz&);            // Not implemented
  AliEveTPCSectorViz& operator=(const AliEveTPCSectorViz&); // Not implemented

  ClassDef(AliEveTPCSectorViz, 0); // Base-class for visualization of data for one TPC sector.
};


// --- Inlines ---

inline UChar_t* AliEveTPCSectorViz::ColorFromArray(Int_t val) const
{
  if(val < fThreshold) val = fThreshold;
  if(val > fMaxVal)    val = fMaxVal;
  return fColorArray + 4 * (val - fThreshold);
}

inline void AliEveTPCSectorViz::ColorFromArray(Int_t val, UChar_t* pix) const
{
  UChar_t* c = ColorFromArray(val);
  pix[0] = c[0]; pix[1] = c[1]; pix[2] = c[2]; pix[3] = c[3];
}

#endif
