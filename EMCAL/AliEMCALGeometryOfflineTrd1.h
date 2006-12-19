#ifndef ALIEMCALGEOMETRYOFFLINETRD1_H
#define ALIEMCALGEOMETRYOFFLINETRD1_H
/* Copyright(c) 1998-2004, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// GeometryOfflineTrd1 class  for EMCAL : singleton
//*-- Author: Aleksei Pavlinov (WSU)

#include <TNamed.h>
#include <TArrayD.h>
#include <TRotation.h>

class TBrowser;
class TString;
class TVector2;
class TVector3;
class TObjArray;
class AliEMCALGeometry;
class AliEMCALShishKebabTrd1Module; 
// --- ROOT system ---
// --- AliRoot header files ---

class AliEMCALGeometryOfflineTrd1 : public TNamed {
 public:
  AliEMCALGeometryOfflineTrd1(const AliEMCALGeometryOfflineTrd1& geom);

  //assignment operator for coding convention
  const AliEMCALGeometryOfflineTrd1 & operator = (const AliEMCALGeometryOfflineTrd1 &) {return *this;}

  virtual ~AliEMCALGeometryOfflineTrd1() { /* nothing */ };
  static   AliEMCALGeometryOfflineTrd1* GetInstance();
  // positon in SuperModule
  TVector3&  PosInSuperModule(int nSupMod, int nModule, int nIphi, int nIeta); 

 public:
  // One Super Module
  void PositionInSuperModule(int iphi, int ieta, double &lphi, double &leta);
  void PositionInSuperModule(int nSupMod, int nModule, int nIphi, int nIeta, double &lphi, double &leta);
  // Position towers(cells)
  TVector3* CellPosition(int absId); // from 0 to fGeometry->GetNCells()
  // Global System
  TRotation* Rotation(Int_t module); // module change from 1 to 12;

  // service methods
  void    PrintSuperModule();       // *MENU*
  void    PrintCell(Int_t absId=1); // *MENU*
  virtual void Browse(TBrowser* b);
  virtual Bool_t  IsFolder() const;

 private:
  AliEMCALGeometryOfflineTrd1();
  void Init();

  static  AliEMCALGeometryOfflineTrd1* fgGeomOfflineTrd1; //!
  AliEMCALGeometry *fGeometry;                            //!  

  Int_t fMaxInEta;                                        //! max module in ETA direction 
  AliEMCALShishKebabTrd1Module* fTrd1Modules[26];         //! 
  // positon in SuperModule
  Int_t     fSMMaxEta;                                    // = 2*fMaxInEta(52 now)
  TVector2* fSMPositionEta;                               //! array of TVector2 [fSMMaxEta]
  TArrayD   fSMPositionPhi;                               //! [24] now; start from 0;
  Double_t  fShiftOnPhi;                                  //! 

  Int_t fNPhiSuperModule;                                //! number phi position for super modules  
  TRotation fSuperModuleRotationZ[6];                    //!
  TString   fNameSuperModuleRotationZ[6];                //!
  TRotation fSuperModuleRotationX;                       //!
  TRotation fSuperModuleRotation[12];                    //! 
  // position of cells in global coordinate system
  TObjArray *fXYZofCells;                                //! 

  ClassDef(AliEMCALGeometryOfflineTrd1, 0) // EMCAL geometry for offline 
};

#endif // ALIEMCALGEOMETRYOFFLINETRD1_H
