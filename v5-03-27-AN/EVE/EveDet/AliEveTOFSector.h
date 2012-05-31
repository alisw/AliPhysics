#ifndef ALIEVETOFSECTOR_H
#define ALIEVETOFSECTOR_H

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

//
// Class to visualize the TOF digit information
// in TOF sector frame
//

#include <TEveQuadSet.h>

class TTree;

class TClonesArray;
class TGeoManager;

class TEveFrameBox;
class TEveRGBAPalette;
class TEveElement;

class AliTOFGeometry;
 
class AliEveTOFSector : public TEveQuadSet
{
public:
  AliEveTOFSector(const Text_t* n="AliEveTOFSector", const Text_t* t=0);
  AliEveTOFSector(TGeoManager *localGeoManager, Int_t nSector);
    
  AliEveTOFSector(TGeoManager *localGeoManager, Int_t nSector,
                  TClonesArray *tofArray);
  AliEveTOFSector(TGeoManager *localGeoManager,
                  Int_t nSector, TTree *tofTree);
  virtual ~AliEveTOFSector();

  virtual void InitModule();
  virtual void SetTrans(); 

  static void      InitStatics();

  void SetSectorID(Int_t id);
  void SetAutoTrans(Bool_t r){fAutoTrans=r;};
  void SetThreshold(Short_t t);
  void SetMaxVal(Int_t mv);
  Bool_t GetPlate(Int_t nPlate) const {return fPlateFlag[nPlate];};
  Short_t GetThreshold() const {return fThreshold;};
  Int_t GetMaxVal() const {return fMaxVal;};
  Bool_t GetAutoTrans() const {return fAutoTrans;};
  Int_t GetSectorID() const {return fSectorID;};
  virtual void DigitSelected(Int_t idx);
  ///////////////////////////////////////////

  void SetPlate(Int_t nPlate, Bool_t r);

protected:

  static Bool_t           fgStaticInitDone;    // flag to check on/off inizialization
  static TEveFrameBox    *fgTOFsectorFrameBox; // EVE container for TOF sector
  static TEveRGBAPalette *fgTOFsectorPalette;  // EVE container for setting of visualization parameters
 
  AliTOFGeometry *fTOFgeometry; // pointer to the TOF geometry container class

  TClonesArray   *fTOFarray;    // TOF digit array container
  TTree          *fTOFtree;     // TOF digit tree container

  Int_t fSector;     // TOF sector index
  //Int_t fPlate;
  //Int_t fStrip;

  Float_t  fDx;     // x position of TOF digit (in TOF strip RF)
  Float_t  fDy;     // y position of TOF digit (in TOF strip RF)
  Float_t  fDz;     // z position of TOF digit (in TOF strip RF)
  ///////////////////////////////

  Bool_t      fAutoTrans; // to choose if visualize the TOF sector in ALICE RF or in local RF
  //Int_t       fMinTime;     
  //Int_t       fMaxTime;
  Short_t     fThreshold; // threshold to cut on visualization
  Int_t       fMaxVal;    // max value to cut on visualization
  Int_t       fSectorID;  // TOF sector identifier
  Bool_t     *fPlateFlag; // flag to switch on/off the TOF module visualization inside a TOF SM

  //Bool_t      fPlateFlag0;
  //Bool_t      fPlateFlag1;
  //Bool_t      fPlateFlag2;
  //Bool_t      fPlateFlag3;
  //Bool_t      fPlateFlag4;

  //Color_t     fFrameColor;
  //Bool_t      fRnrFrame;

  TGeoManager *fGeoManager; // pointer to the ALICE geometry

private:
  void LoadQuads();

  AliEveTOFSector(const AliEveTOFSector&);            // Not implemented
  AliEveTOFSector& operator=(const AliEveTOFSector&); // Not implemented

  ClassDef(AliEveTOFSector, 0); // Representation of a TOF sector.
};
#endif
