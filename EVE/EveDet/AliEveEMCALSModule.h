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
  AliEveEMCALData         *fEMCALData;        //  Data for the current event
  AliEveEMCALSModuleData  *fEMCALSModuleData; //  Data of Super Module (SM)
  Color_t                 fFrameColor;        //  Main coloring
  UInt_t                  fRTS;               //! Rendering Time Stamp
  Int_t                   fSModuleID;         //  Id of super module, 0 to 11
  TEveQuadSet             *fQuadSet;          //  Digit container
  TEveQuadSet             *fQuadSet2;         //  Cluster container
  TEvePointSet            *fPointSet;         //  Hit container
  Short_t                 fThreshold;         //  Digit amplitude threshold
  Int_t                   fMaxVal;            //  Digit amplitude maximum value
  Int_t                   fClusterSize;       //  Cluster point size
  Int_t                   fHitSize;           //  Hit point size
  mutable UChar_t         *fColorArray;       //  Color-cache
  Int_t                   fDebug;             //  Debug option

  Float_t                 fSMBigBBox[3];      //  Bounding Box of full SM
  Float_t                 fSMSmallBBox[3];    //  Bounding Box of half SM
  Float_t                 fSMBBoxCenter[3];   //  Bounding Box Center of full SM

  static   Bool_t          fStaticInit;        // Flag for static variable initialization.
  static   TEveFrameBox    *fFrameBigBox;      // Frame box per full SM
  static   TEveFrameBox    *fFrameSmallBox;    // Frame box per half SM
  static   TEveRGBAPalette *fFrameDigPalette;  // Signal to color mapping for EMCAL digits
  static   TEveRGBAPalette *fFrameCluPalette;  // Signal to color mapping for EMCAL clusters

  void SetupColor(Int_t val, UChar_t* pix) const;

  //****** Not used yet **************
  void     ClearColorArray();
  void     SetupColorArray() const;
  UChar_t* ColorFromArray(Int_t val) const;
  void     ColorFromArray(Int_t val, UChar_t* pix) const;
  Int_t    ColorIndex(Int_t val) const;
  //**********************************

 private:
  AliEveEMCALSModule(const AliEveEMCALSModule&);            
  AliEveEMCALSModule& operator=(const AliEveEMCALSModule&); // Not implemented

  ClassDef(AliEveEMCALSModule, 0); // Base class for TRD hits visualisation
};

#endif
