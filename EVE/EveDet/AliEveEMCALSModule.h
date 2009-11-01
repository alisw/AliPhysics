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
                           public TAtt3D
{

 public:
  AliEveEMCALSModule(Int_t smodid, const Text_t* n, const Text_t* t);
  ~AliEveEMCALSModule();

  void DropData();

  virtual Bool_t CanEditMainColor() const { return kTRUE; }

  void  SetDataSource(AliEveEMCALData *data);
  void  SetSModuleID(Int_t id);
  void  SetFrameColor(Color_t col) { fFrameColor = col; };
  AliEveEMCALData* GetData() const { return fEMCALData; };
  AliEveEMCALSModuleData* GetSModuleData() const;
  Int_t GetID() const { return fSModuleID; };
  void  SetClusterSize(Int_t size);
  void  SetHitSize(Int_t size);

  void UpdateQuads();

 protected:
  AliEveEMCALData         *fEMCALData;        //  Data for the current event
  AliEveEMCALSModuleData  *fEMCALSModuleData; //  Data of Super Module (SM)
  Color_t                 fFrameColor;        //  Main coloring
  Int_t                   fSModuleID;         //  Id of super module, 0 to 11
  TEveQuadSet             *fQuadSet;          //  Digit container
  TEveQuadSet             *fQuadSet2;         //  Cluster container
  TEvePointSet            *fPointSet;         //  Hit container
  Int_t                   fClusterSize;       //  Cluster point size
  Int_t                   fHitSize;           //  Hit point size
  Int_t                   fDebug;             //  Debug option

  static void InitStatics(AliEveEMCALSModuleData* md);

  static   Bool_t           fStaticInit;       // Flag for static variable initialization.
  static   Float_t          fgSMBigBBox[3];    //  Bounding Box of full SM
  static   Float_t          fgSMSmallBBox[3];  //  Bounding Box of half SM
  static   TEveFrameBox    *fgFrameBigBox;     // Frame box per full SM
  static   TEveFrameBox    *fgFrameSmallBox;   // Frame box per half SM
  static   TEveRGBAPalette *fgFrameDigPalette; // Signal to color mapping for EMCAL digits
  static   TEveRGBAPalette *fgFrameCluPalette; // Signal to color mapping for EMCAL clusters

  void SetupColor(Int_t val, UChar_t* pix) const;

 private:
  AliEveEMCALSModule(const AliEveEMCALSModule&);            
  AliEveEMCALSModule& operator=(const AliEveEMCALSModule&); // Not implemented

  ClassDef(AliEveEMCALSModule, 0); // Base class for TRD hits visualisation
};

#endif
