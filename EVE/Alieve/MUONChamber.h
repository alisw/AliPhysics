#ifndef ALIEVE_MUONChamber_H
#define ALIEVE_MUONChamber_H

#include <Reve/RenderElement.h>
#include <Reve/QuadSet.h>

#include <TNamed.h>
#include <TAtt3D.h>
#include <TAttBBox.h>

namespace Alieve {

class MUONData;
class MUONChamberData;
class MUONChamberEditor;
class MUONChamberGL;

class MUONChamber : public Reve::RenderElement,
                    public TNamed,
                    public TAtt3D,
                    public TAttBBox

{

  friend class MUONChamberGL;
  friend class MUONChamberEditor;

  MUONChamber(const MUONChamber&);            // Not implemented
  MUONChamber& operator=(const MUONChamber&); // Not implemented

 protected:

  void UpdateQuads();

  MUONData*         fMUONData;      // data for the current event
  Color_t           fFrameColor;    // main coloring
  UInt_t            fRTS;           //! Rendering Time Stamp
  Int_t             fChamberID;     // number of the chamber, 0 to 13
  Reve::OldQuadSet  fQuadSet1;      // 1st cathode plane digits
  Reve::OldQuadSet  fQuadSet2;      // 2nd cathode plane digits
  Short_t           fThreshold;     // digit amplitude threshold
  Int_t             fMaxVal;        // digit amplitude maximum value

  void SetupColor(Int_t val, UChar_t* pix) const;

  mutable UChar_t* fColorArray;
  void ClearColorArray();
  void SetupColorArray() const;
  UChar_t* ColorFromArray(Int_t val) const;
  void     ColorFromArray(Int_t val, UChar_t* pix) const;
  Int_t    ColorIndex(Int_t val) const;

public:

  MUONChamber(const Text_t* n = "MUONChamber", const Text_t* t = 0);
  virtual ~MUONChamber();

  virtual void ComputeBBox();
  virtual void Paint(Option_t* option = "");
  virtual UInt_t IncRTS()     { return ++fRTS; };
  virtual Bool_t CanEditMainColor() { return kTRUE; }

  void SetDataSource(MUONData *data);
  void SetChamberID(Int_t id);
  void SetFrameColor(Color_t col)     { fFrameColor = col; IncRTS(); };
  MUONData* GetData() const { return fMUONData; };
  MUONChamberData* GetChamberData() const;
  Int_t GetID() const { return fChamberID; };
  void SetThreshold(Short_t t);
  void SetMaxVal(Int_t mv);

  ClassDef(MUONChamber,1);  // Visualisation of the MUON chambers

};

inline UChar_t* MUONChamber::ColorFromArray(Int_t val) const
{
  if(val < fThreshold) val = fThreshold;
  if(val > fMaxVal)    val = fMaxVal;
  return fColorArray + 4 * (val - fThreshold);
}

inline void MUONChamber::ColorFromArray(Int_t val, UChar_t* pix) const
{
  UChar_t* c = ColorFromArray(val);
  pix[0] = c[0]; pix[1] = c[1]; pix[2] = c[2]; pix[3] = c[3];
}

}

#endif
