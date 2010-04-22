// $Id: AliJRunHeader.h,v 1.1 2008/02/04 13:28:47 rak Exp $
#ifndef ALIJRUNHEADER_H
#define ALIJRUNHEADER_H
////////////////////////////////////////////////////
/*!
        \file AliJRunHeader.h
        \brief
        \author J. Rak, D.J.Kim, F.Krizek  (Jyvaskyla || HIP)
        \email: djkim@cc.jyu.fi
        \version $Revision: 1.1 $
        \date $Date: 2008/02/04 13:28:47 $
*/
////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include <TObject.h>
#include <TNamed.h>
#endif

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>

#include "JConst.h"

using namespace std;

//class TObjString;
//class TObject;

class AliJRunHeader : public TNamed {

public:
  AliJRunHeader();//constructor
  AliJRunHeader(const AliJRunHeader& ap);

  virtual ~AliJRunHeader(){;}              //destructor

  virtual Int_t  GetRunNumber()     const {return fRunNumber;}
  virtual void SetRunNumber(Int_t runN) { fRunNumber = runN;}

  //         s e t t e r s   a n d    g e t t e r s
  void SetL3Field(Short_t polarity,Double_t MagnetFieldInL3){
    fL3MagnetPolarity = polarity;
    fMagneticFieldL3  = MagnetFieldInL3;
  }

  Short_t  GetL3MagnetFieldPolarity()  const { return fL3MagnetPolarity;}
  Double_t GetL3MagnetFieldIntensity() const { return fMagneticFieldL3;}

  //-- Alice trigger table --
  void SetActiveTriggersAlice(TString *triggers);

  Int_t GetActiveTriggerBitAlice(TString TriggerName);

  TString GetActiveTriggerAlice(Int_t TriggerBit) const {
    return ((TObjString*) (fActiveTriggersAlice.At(TriggerBit)))->GetString();
  }

  //-- JCorran trigger table -- 
  void SetActiveTriggersJCorran(TString *triggers, Int_t range);

  TString GetActiveTriggerJCorran(Int_t TriggerBit) const {
    return ((TObjString*) (fActiveTriggersJCorran.At(TriggerBit)))->GetString();
  }

  void PrintOut();

  AliJRunHeader& operator=(const  AliJRunHeader& header);

protected:
  Int_t     fRunNumber;        //run number 
  Short_t   fL3MagnetPolarity; //Polarity of magnetic filed in L3 magnet (LHC convention: + -> +Bz)
  Double_t  fMagneticFieldL3;  //Solenoid Magnetic Field in kG   
  TObjArray fActiveTriggersAlice;   //array maping between trigger bit and trigger names

  Int_t     fSizeOfTableJCorran;  //size of jcorran table
  TObjArray fActiveTriggersJCorran;   //array maping between trigger bit and trigger names
                                      //TBit 0 = MB 
  ClassDef(AliJRunHeader,1)

};

#endif
