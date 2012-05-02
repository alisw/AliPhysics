// $Id: AliJRunHeader.h,v 1.1 2008/02/04 13:28:47 rak Exp $
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

#ifndef ALIJRUNHEADER_H
#define ALIJRUNHEADER_H
#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <TNamed.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <vector>
#include <map>


class AliJRunHeader : public TNamed {

//#pragma link C++ class AliJRunHeader+;

 public:
  AliJRunHeader();//constructor
  AliJRunHeader(const AliJRunHeader& ap);

  virtual ~AliJRunHeader(){;}              //destructor

  virtual Int_t  GetRunNumber()     const {return fRunNumber;}
  virtual void SetRunNumber(Int_t runN) { fRunNumber = runN;}


  void SetRunType(const TString info ) { fRunType = info; }
  TString GetRunType() const { return fRunType; }

  void SetESDInfo(const TString info ) { fESDInfo = info; }
  TString GetESDInfo() const { return fESDInfo; }

//  void SetESDInfo(const char* info ) { fESDInfo = info; }
//  char* GetESDInfo() const { return fESDInfo; }

  //         s e t t e r s   a n d    g e t t e r s
  void SetL3Field(Short_t polarity,Double_t MagnetFieldInL3){
    fL3MagnetPolarity = polarity;
    fMagneticFieldL3  = MagnetFieldInL3;
  }

  Short_t  GetL3MagnetFieldPolarity()  const { return fL3MagnetPolarity;}
  Double_t GetL3MagnetFieldIntensity() const { return fMagneticFieldL3;}

  //-- Alice trigger table -- by Filip. "Trigger Class" like "+CMBACS2-B-NOPF-ALL"
  void SetActiveTriggersAlice( const TString *triggers);
  Int_t GetActiveTriggerBitAlice(TString TriggerName);
  TString GetActiveTriggerAlice(Int_t TriggerBit) const {
    return ((TObjString*) (fActiveTriggersAlice.At(TriggerBit)))->GetString();
  }

  //-- JCorran trigger table -- by Filip
  void SetActiveTriggersJCorran(const TString *triggers, Int_t range);
  TString GetActiveTriggerJCorran(Int_t TriggerBit) const {
    return ((TObjString*) (fActiveTriggersJCorran.At(TriggerBit)))->GetString();
  }

  AliJRunHeader& operator=(const  AliJRunHeader& header);

  void PrintOut();

 protected:
  Int_t       fRunNumber;        //run number 
  TString     fRunType;       // ex) LHC10h
  TString     fESDInfo;       // information of aliroot,  root version while esd production
//  Char_t*     fESDInfo;       // information of aliroot,  root version while esd production
  Short_t     fL3MagnetPolarity; //Polarity of magnetic filed in L3 magnet (LHC convention: + -> +Bz)
  Double32_t  fMagneticFieldL3;  //Solenoid Magnetic Field in kG   
  TObjArray   fActiveTriggersAlice;   //array maping between trigger bit and trigger names

  // AOD
  TString     fFiredTriggers;       // String with fired triggers
  ULong64_t   fTriggerMask;         // Trigger Type (mask)
  UChar_t     fTriggerCluster;      // Trigger cluster (mask)


  Int_t       fSizeOfTableJCorran;  //size of jcorran table
  TObjArray   fActiveTriggersJCorran;   //array maping between trigger bit and trigger names
  //TBit 0 = MB 
  ClassDef(AliJRunHeader,2)

};

#endif
