#ifndef AliEMCALSensorTemp_H
#define AliEMCALSensorTemp_H
/* Copyright(c) 2006-07, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


////////////////////////////////////////////////////////////////////////////
//              Container class for temperature sensor positions            
////////////////////////////////////////////////////////////////////////////


#include "TMath.h"
#include "AliSplineFit.h"
#include "AliDCSSensor.h"
#include "TTree.h"

class TObject;
class TClonesArray;
class TObjArray;
class TGraph;
class TVector3;
class TFile;
class TString;
class TTimeStamp;


////////////////////////////////////////////////////////////////////////
//              Class AliEMCALSensorTempSensors
////////////////////////////////////////////////////////////////////////

const TString kAmandaString = "EMC_PT_%02d.Temperature";

class AliEMCALSensorTemp : public AliDCSSensor {

public:
  AliEMCALSensorTemp();
  AliEMCALSensorTemp(const AliEMCALSensorTemp& source);
  virtual ~AliEMCALSensorTemp(){}
  AliEMCALSensorTemp& operator=(const AliEMCALSensorTemp& source);
  
  Int_t       GetSide()   const {return fSide;	 }
  Int_t       GetSector() const {return fSector;	 }
  Int_t       GetNum()	  const {return fNum;	 }

  void SetSide   (Int_t side)      {fSide   = side;	 }
  void SetSector (Int_t sector)    {fSector = sector;}
  void SetNum    (Int_t num)       {fNum    = num;   }


  static TClonesArray * ReadList(const char *fname,
                                 const TString& amandaString = kAmandaString);
  static TClonesArray * ReadTree(TTree *tree, 
                                 const TString& amandaString = kAmandaString);

protected:
  // A SuperModule is defined in hardware land with a sector and a side index
  Int_t fSide;      // EMCAL side; 0:Shaft Side (A) -- 1:Muon Side (C)
  Int_t fSector;    // Number of sector             (0-5)
  Int_t fNum;       // Number within a SuperModule: 8 sensors => index range 0-7

  ClassDef(AliEMCALSensorTemp,1)
};
#endif
