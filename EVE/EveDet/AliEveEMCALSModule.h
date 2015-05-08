#ifndef ALIEVEEMCALSMODULE_H
#define ALIEVEEMCALSMODULE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///
/// \class AliEveEMCALSModule
/// \brief Visualization of an EMCAL super module.
///
///  Visualization of an EMCAL super module for event display.
///
/// \author Magali Estienne <magali.estienne@cern.ch>, SUBATECH. EMCal implementation, June 2008
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-IN2P3-CNRS. DCal implementation + doxygen, May 2015.
///

#include "AliEveEMCALSModuleData.h"

class TStyle;
class TBuffer3DTypes;
class TBuffer3D;
class TVirtualPad;
class TVirtualViewer3D;
class TClonesArray;
class TTree;
class TGedFrame;

class TEveQuadSet;
class TEveBoxSet;
class TEveFrameBox;
class TEvePointSet;
class TEveTrans;

class AliRun;
class AliESDEvent;

class AliEveEMCALData;

class AliEveEMCALSModule : public TEveElement,
                           public TNamed,
                           public TAtt3D
{

 public:
  
  AliEveEMCALSModule(Int_t smodid, const Text_t* n, const Text_t* t);
  
  ~AliEveEMCALSModule();

  virtual Bool_t CanEditMainColor()        const { return kTRUE      ; } // Remove?

  void  SetDataSource(AliEveEMCALData * data);
  
  void  SetSModuleID(Int_t id);
  
  void  SetFrameColor(Color_t col)               { fFrameColor = col ; }
  
  const AliEveEMCALData* GetData()         const { return fEMCALData ; }
  
  void  DropData()                               { fEMCALSModuleData->DropData() ; }
  
  AliEveEMCALSModuleData* GetSModuleData() const ;
  
  Int_t GetID()                            const { return fSModuleID ; }
  
  void  SetClusterSize(Int_t size);
  
  void  SetHitSize(Int_t size);

  void  UpdateQuads(Bool_t iHits, Bool_t iDigits, Bool_t iClusters);

  TEveQuadSet * GetDigitQuadSet()          const { return fQuadSet   ; }
  
  TEveQuadSet * GetClusterQuadSet()        const { return fQuadSet2  ; }
  
  TEvePointSet* GetHitPointSet()           const { return fPointSet  ; }
    
 protected:
  
  AliEveEMCALData          *fEMCALData;        ///<  Data for the current event
  AliEveEMCALSModuleData   *fEMCALSModuleData; ///<  Data of Super Module (SM)
  Color_t                   fFrameColor;       ///<  Main coloring
  Int_t                     fSModuleID;        ///<  Id of super module, 0 to 11
  TEveQuadSet              *fQuadSet;          ///<  Digit container
  TEveQuadSet              *fQuadSet2;         ///<  Cluster container
  TEvePointSet             *fPointSet;         ///<  Hit container
  Int_t                     fClusterSize;      ///<  Cluster point size
  Int_t                     fHitSize;          ///<  Hit point size
  Int_t                     fDebug;            ///<  Debug option

  static void InitStatics(AliEveEMCALSModuleData* md);

  static   Bool_t           fgStaticInit;      ///< Flag for static variable initialization.
  
  static   Float_t          fgSMBigBBox[3];    ///<  Bounding Box of full SM
  static   Float_t          fgSMSmallBBox[3];  ///<  Bounding Box of 1/3 SM  
  static   Float_t          fgSMDCalBBox[3];   ///<  Bounding Box of DCal SM
  static   Float_t          fgSMSmallDBBox[3]; ///<  Bounding Box of 1/3 DCal SM
  
  static   TEveFrameBox    *fgFrameBigBox;     ///< Frame box per full SM
  static   TEveFrameBox    *fgFrameSmallBox;   ///< Frame box per 1/3 SM
  static   TEveFrameBox    *fgFrameDCalBox;    ///< Frame box per DCal SM
  static   TEveFrameBox    *fgFrameSmallDBox;  ///< Frame box per 1/3 DCal SM
  
  static   TEveRGBAPalette *fgFrameDigPalette; ///< Signal to color mapping for EMCAL digits
  static   TEveRGBAPalette *fgFrameCluPalette; ///< Signal to color mapping for EMCAL clusters

  void SetupColor(Int_t val, UChar_t* pix) const;

 private:
  
  AliEveEMCALSModule           (const AliEveEMCALSModule &esm);  
  
  /// Assignment operator not implemented.
  AliEveEMCALSModule& operator=(const AliEveEMCALSModule &esm); 

  /// \cond CLASSIMP
  ClassDef(AliEveEMCALSModule, 1) ;
  /// \endcond

};

#endif //ALIEVEEMCALSMODULE_H
