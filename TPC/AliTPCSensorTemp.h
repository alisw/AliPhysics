#ifndef AliTPCSensorTemp_H
#define AliTPCSensorTemp_H
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
//              Class AliTPCSensorTempSensors
////////////////////////////////////////////////////////////////////////

const TString kAmandaString = "tpc_temp:PT_%d.Temperature";

class AliTPCSensorTemp : public AliDCSSensor {

public:
  AliTPCSensorTemp();
  AliTPCSensorTemp(const AliTPCSensorTemp& source);
  virtual ~AliTPCSensorTemp(){}
  AliTPCSensorTemp& operator=(const AliTPCSensorTemp& source);
  
  Int_t       GetType()   const {return fType;   }
  Int_t       GetSide()   const {return fSide;	 }
  Int_t       GetSector() const {return fSector;	 }
  Int_t       GetNum()	  const {return fNum;	 }

  void SetType   (Int_t type)      {fType   = type;   }
  void SetSide   (Int_t side)      {fSide   = side;	 }
  void SetSector (Int_t sector)    {fSector = sector;}
  void SetNum    (Int_t num)       {fNum    = num;   }


  static TClonesArray * ReadList(const char *fname,
                                 const TString& amandaString = kAmandaString);
  static TClonesArray * ReadTree(TTree *tree, 
                                 const TString& amandaString = kAmandaString);

protected:
  Int_t fType;      // Position of sensors on fieldcage
                    //  (0=ROC,1=OFC,2=IFC,3=TPC,4=ELM,5=TS,6=COOL)
  Int_t fSide;      // TPC side; 0:Shaft Side (A) -- 1:Muon Side (C)
  Int_t fSector;    // Number of sector             (0-17)
  Int_t fNum;       // Position depands from type of sensor.
                    //    fType=0(0-4) from inside to outside
                    //    fType=1(0-5) fom A side to C side
                    //    fType=2(0-5) fom A side to C side
                    //    fType=3(0)   one per sector 
                    //    fType=4()
                    //    fType=5()
                    //	fType=6(0-1) 0:input -- 1:output


  ClassDef(AliTPCSensorTemp,1)
};
#endif
