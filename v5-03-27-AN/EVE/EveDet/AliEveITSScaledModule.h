// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveITSScaledModule_H
#define AliEveITSScaledModule_H

#include <TEveUtil.h>
#include <AliEveITSModule.h>

#include <TQObject.h>


/******************************************************************************/
// AliEveDigitScaleInfo
/******************************************************************************/

class AliEveDigitScaleInfo : public TQObject, public TEveRefBackPtr
{
public:
  enum StatType_e { kSTOccup, kSTAverage, kSTRms };

  AliEveDigitScaleInfo();
  virtual ~AliEveDigitScaleInfo() {}

  Int_t  GetScale() const { return fScale; }
  void   ScaleChanged(Int_t s);

  Int_t  GetStatType() const { return fStatType; }
  void   StatTypeChanged(Int_t t);

  Bool_t GetSyncPalette() const   { return fSyncPalette; }
  void   SetSyncPalette(Bool_t x) { fSyncPalette = x; }

protected:
  Int_t            fScale;        // Current scale.
  Int_t            fStatType;     // Digit scaling algorithm, see StatType_e.

  Bool_t           fSyncPalette;  // Synchronize palette on next usage.

private:
  AliEveDigitScaleInfo(const AliEveDigitScaleInfo&);            // Not implemented
  AliEveDigitScaleInfo& operator=(const AliEveDigitScaleInfo&); // Not implemented

  ClassDef(AliEveDigitScaleInfo, 0);
};

/******************************************************************************/
// ScaledDigit
/******************************************************************************/

/******************************************************************************/
// AliEveITSScaledModule
/******************************************************************************/

class AliEveITSScaledModule : public AliEveITSModule
{
  friend class ITSSDSubEditor;

public:
  AliEveITSScaledModule(Int_t gid, AliEveITSDigitsInfo* info, AliEveDigitScaleInfo* si );
  virtual ~AliEveITSScaledModule();

  virtual void DigitSelected(Int_t idx);

  virtual void LoadQuads();
  void         SetQuadValues();

  void         SyncPalette();

  void         GetScaleData(Int_t& cnx, Int_t& cnz, Int_t& total) const;
  AliEveDigitScaleInfo*  GetScaleInfo() { return fScaleInfo; }


  // --- Inner structs

  struct ScaledDigit_t : public TObject
  {
  public:
    Int_t   fN;
    Float_t fSum;
    Float_t fSqrSum;
    Int_t   fMinI, fMinJ;
    Int_t   fMaxI, fMaxJ;

    ScaledDigit_t();
    ScaledDigit_t(Int_t di, Int_t dj);

    void Dump() const;
  };

protected:
  Int_t       fNx;   // per module
  Int_t       fNz;

  Int_t       fNCx;  // per cell
  Int_t       fNCz;

  AliEveDigitScaleInfo* fScaleInfo;

private:
  std::map<Int_t, ScaledDigit_t> fDigitsMap;

  AliEveITSScaledModule(const AliEveITSScaledModule&);            // Not implemented
  AliEveITSScaledModule& operator=(const AliEveITSScaledModule&); // Not implemented

  ClassDef(AliEveITSScaledModule, 0); // Visualization of an ITS module with digits aggregated on a grid of pre-defined size.
};

#endif
