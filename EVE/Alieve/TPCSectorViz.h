// $Header$

#ifndef ALIEVE_TPCSectorViz_H
#define ALIEVE_TPCSectorViz_H

#include <Reve/RenderElement.h>

#include <TNamed.h> 
#include <TAtt3D.h>
#include <TAttBBox.h>


namespace Alieve {

class TPCData; class TPCSectorData;

class TPCSectorVizEditor;
class TPCSector2D;  class TPCSector2DEditor;  class TPCSector2DGL;
class TPCSector3D;  class TPCSector3DEditor;  class TPCSector3DGL;

class TPCSectorViz : public Reve::RenderElement,
                     public TNamed,
                     public TAtt3D,
                     public TAttBBox
{
  friend class TPCSectorVizEditor;
  friend class TPCSector2D;
  friend class TPCSector2DEditor;
  friend class TPCSector2DGL;
  friend class TPCSector3D;
  friend class TPCSector3DEditor;
  friend class TPCSector3DGL;

  TPCSectorViz(const TPCSectorViz&);            // Not implemented
  TPCSectorViz& operator=(const TPCSectorViz&); // Not implemented

protected:
  TPCData*    fTPCData; 
  Int_t       fSectorID;

  Int_t       fMinTime;     
  Int_t       fMaxTime;
  Short_t     fThreshold;
  Int_t       fMaxVal;

  Bool_t      fRnrInn;
  Bool_t      fRnrOut1;
  Bool_t      fRnrOut2;

  Color_t     fFrameColor;
  Bool_t      fRnrFrame;
  Bool_t      fTrans;
  Double_t    fMatrix[16];
  UInt_t      fRTS;       //! Rendering TimeStamp

  void SetupColor(Int_t val, UChar_t* pix) const;

  mutable UChar_t* fColorArray;
  void ClearColorArray();
  void SetupColorArray() const;
  UChar_t* ColorFromArray(Int_t val) const;
  void     ColorFromArray(Int_t val, UChar_t* pix) const;

public:
  TPCSectorViz(const Text_t* n="TPCSectorViz", const Text_t* t=0);
  virtual ~TPCSectorViz();

  virtual void CopyVizParams(const TPCSectorViz& v);

  virtual UInt_t IncRTS()           { return ++fRTS; }
  virtual Bool_t CanEditMainColor() { return true; }

  void SetDataSource(TPCData* data);
  void SetSectorID(Int_t id);

  TPCData*       GetData()     const { return fTPCData; }
  Int_t          GetSectorID() const { return fSectorID; }
  TPCSectorData* GetSectorData() const;

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
  void SetTrans(Bool_t t);

  ClassDef(TPCSectorViz, 1); // Base-class for TPC raw-data visualization
}; // endclass TPCSectorViz


inline UChar_t* TPCSectorViz::ColorFromArray(Int_t val) const
{
  if(val < fThreshold) val = fThreshold;
  if(val > fMaxVal)    val = fMaxVal;
  return fColorArray + 4 * (val - fThreshold);
}

inline void TPCSectorViz::ColorFromArray(Int_t val, UChar_t* pix) const
{
  UChar_t* c = ColorFromArray(val);
  pix[0] = c[0]; pix[1] = c[1]; pix[2] = c[2]; pix[3] = c[3];
}

}

#endif
