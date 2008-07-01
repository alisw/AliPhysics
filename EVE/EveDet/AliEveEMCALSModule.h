//*************************************************************************
// EMCAL event display
// Visualization of an EMCAL super module.
//
//  Author: Magali Estienne (magali.estienne@cern.ch)
//  June 30 2008
//*************************************************************************

#ifndef ALIEVEEMCALSMODULE_H
#define ALIEVEEMCALSMODULE_H

#include <TEveQuadSet.h>
#include <TEveElement.h>
#include <TEveBoxSet.h>
#include <TEveFrameBox.h>
#include <TEvePointSet.h>
#include <TClonesArray.h>
#include <TTree.h>

#include <TGedFrame.h>

class AliEveEMCALData;
class AliEveEMCALSModuleData;


class AliEveEMCALSModule : public TEveElement,
                           public TNamed,
                           public TAtt3D,
			   public TAttBBox
{

 public:
  AliEveEMCALSModule(Int_t smodid, const Text_t* n, const Text_t* t);
  ~AliEveEMCALSModule();

  void DropData();

  virtual void   ComputeBBox();
  virtual UInt_t IncRTS()     { return ++fRTS; };
  virtual Bool_t CanEditMainColor() const { return kTRUE; }

  void  SetDataSource(AliEveEMCALData *data);
  void  SetSModuleID(Int_t id);
  void  SetFrameColor(Color_t col) { fFrameColor = col; IncRTS(); };
  AliEveEMCALData* GetData() const { return fEMCALData; };
  AliEveEMCALSModuleData* GetSModuleData() const;
  Int_t GetID() const { return fSModuleID; };
  void  SetThreshold(Short_t t);
  void  SetMaxVal(Int_t mv);
  void  SetClusterSize(Int_t size);
  void  SetHitSize(Int_t size);

  void UpdateQuads();

 protected:
  AliEveEMCALData         *fEMCALData;      // data for the current event
  AliEveEMCALSModuleData  *fEMCALSModuleData; // data of super module
  Color_t                 fFrameColor;    // main coloring
  UInt_t                  fRTS;           //! Rendering Time Stamp
  Int_t                   fSModuleID;     // Id of super module, 0 to 11
  TEveQuadSet             *fQuadSet;      // 1st cathode plane digits
  TEveQuadSet             *fQuadSet2;      // 2nd cathode plane digits
  TEvePointSet            *fPointSet;     // reconstructed points (1st cathode)
  Short_t                 fThreshold;     // digit amplitude threshold
  Int_t                   fMaxVal;        // digit amplitude maximum value
  Int_t                   fClusterSize;   // cluster point size
  Int_t                   fHitSize;       // hit point size
  mutable UChar_t         *fColorArray;    // color-cache
  Int_t                   fDebug;         // Debug option

  Float_t                 fSMBigBBox[3];
  Float_t                 fSMSmallBBox[3];
  Float_t                 fSMBBoxCenter[3];

  static   TEveFrameBox    *fFrameBigBox;      // Frame box per super module
  static   TEveFrameBox    *fFrameSmallBox;    // Frame box per super module
  static   TEveRGBAPalette *fFrameDigPalette;  // Signal to color mapping for EMCAL
  static   TEveRGBAPalette *fFrameCluPalette;  // Signal to color mapping for EMCAL

  void SetupColor(Int_t val, UChar_t* pix) const;

  //****** Not used yet **************
  void     ClearColorArray();
  void     SetupColorArray() const;
  UChar_t* ColorFromArray(Int_t val) const;
  void     ColorFromArray(Int_t val, UChar_t* pix) const;
  Int_t    ColorIndex(Int_t val) const;
  //**********************************

 private:
  AliEveEMCALSModule(const AliEveEMCALSModule&);            // Not implemented
  AliEveEMCALSModule& operator=(const AliEveEMCALSModule&); // Not implemented

  ClassDef(AliEveEMCALSModule, 0); // Base class for TRD hits visualisation
};

#endif
