// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel & Bogdan Vulpescu: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveMUONChamber_H
#define AliEveMUONChamber_H

#include <TEveElement.h>
#include <TEveQuadSet.h>
#include <TEvePointSet.h>

#include <TNamed.h>
#include <TAtt3D.h>
#include <TAttBBox.h>


class AliEveMUONData;
class AliEveMUONChamberData;
class AliEveMUONChamberEditor;
class AliEveMUONChamberGL;

class AliEveMUONChamber : public TEveElement,
			  public TNamed,
			  public TAtt3D,
			  public TAttBBox
{
  friend class AliEveMUONChamberGL;
  friend class AliEveMUONChamberEditor;

public:
  AliEveMUONChamber(Int_t id, const Text_t* n = "AliEveMUONChamber", const Text_t* t = 0);
  virtual ~AliEveMUONChamber();

  virtual void   ComputeBBox();
  virtual void   Paint(Option_t* option = "");
  virtual UInt_t IncRTS()     { return ++fRTS; };
  virtual Bool_t CanEditMainColor() const { return kTRUE; }

  void SetDataSource(AliEveMUONData *data);
  void SetChamberID(Int_t id);
  void SetFrameColor(Color_t col)     { fFrameColor = col; IncRTS(); };
  AliEveMUONData* GetData() const { return fMUONData; };
  AliEveMUONChamberData* GetChamberData() const;
  Int_t GetID() const { return fChamberID; };
  void  SetThreshold(Short_t t);
  void  SetMaxVal(Int_t mv);
  void  SetClusterSize(Int_t size);
  void  SetHitSize(Int_t size);

protected:
  void UpdateQuads();

  AliEveMUONData   *fMUONData;      // data for the current event
  Color_t           fFrameColor;    // main coloring
  UInt_t            fRTS;           //! Rendering Time Stamp
  Int_t             fChamberID;     // number of the chamber, 0 to 13
  TEveQuadSet       fQuadSet1;      // 1st cathode plane digits
  TEveQuadSet       fQuadSet2;      // 2nd cathode plane digits
  TEvePointSet      fPointSet1;     // reconstructed points (1st cathode)
  TEvePointSet      fPointSet2;     // simulation hits
  Short_t           fThreshold;     // digit amplitude threshold
  Int_t             fMaxVal;        // digit amplitude maximum value
  Int_t             fClusterSize;   // cluster point size
  Int_t             fHitSize;       // hit point size

  void SetupColor(Int_t val, UChar_t* pix) const;

  mutable UChar_t  *fColorArray;    // color-cache

  void     ClearColorArray();
  void     SetupColorArray() const;
  UChar_t* ColorFromArray(Int_t val) const;
  void     ColorFromArray(Int_t val, UChar_t* pix) const;
  Int_t    ColorIndex(Int_t val) const;

private:
  AliEveMUONChamber(const AliEveMUONChamber&);            // Not implemented
  AliEveMUONChamber& operator=(const AliEveMUONChamber&); // Not implemented

  ClassDef(AliEveMUONChamber, 0);  // Visualisation of the MUON chambers
};


// --- Inlines ---

inline UChar_t* AliEveMUONChamber::ColorFromArray(Int_t val) const
{
  if(val < fThreshold) val = fThreshold;
  if(val > fMaxVal)    val = fMaxVal;
  return fColorArray + 4 * (val - fThreshold);
}

inline void AliEveMUONChamber::ColorFromArray(Int_t val, UChar_t* pix) const
{
  UChar_t* c = ColorFromArray(val);
  pix[0] = c[0]; pix[1] = c[1]; pix[2] = c[2]; pix[3] = c[3];
}

#endif
