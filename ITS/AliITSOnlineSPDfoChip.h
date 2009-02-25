#ifndef ALIITSONLINESPDFOCHIP_H
#define ALIITSONLINESPDFOCHIP_H  
/* Copyright(c) 2008-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


/////////////////////////////////////////////////////////////////
// Author: A. Mastroserio                                      //
// This class is the container for FastOR online calibration.  //
//                                                             //
/////////////////////////////////////////////////////////////////

#include <TObject.h>

class TObjArray;
class AliITSOnlineSPDfoChipConfig;

class AliITSOnlineSPDfoChip : public TObject{

 public:
  AliITSOnlineSPDfoChip();//ctor
  AliITSOnlineSPDfoChip(Short_t nparams); 
  
  virtual ~AliITSOnlineSPDfoChip(); //dctor

  // SETTERS
  void SetActiveHS(Int_t hs)               {fActiveHS=hs;}
  void SetChipId(Int_t chipId)             {fChipId=chipId;}
  void SetDACParameter(Int_t i, UShort_t par) {fDACparams[i] = par;}
  void AddMeasurement(AliITSOnlineSPDfoChipConfig *ChipConfiginfo);

  // GETTERS
  Short_t GetActiveHS() const                 {return fActiveHS;}
  Short_t GetChipId()  const                  {return fChipId;} 
  Short_t GetNumberOfDACParams() const        {return fNumDACparams;}
  Short_t GetDACParameters(Int_t ipar)  const {return fDACparams[ipar];}
  Int_t GetNumberOfChipConfigs() const        {return fNumChipConfigs;} 
  TObjArray* GetChipConfigInfo() const        {return fChipConfigArray;}
   
  void PrintInfo();              // prints the container content
  
 protected:
  Short_t     fActiveHS;         //number of the activated HS
  Short_t     fChipId;           //id number of the chip
  Int_t       fNumDACparams;     //number of DAC parameters to be considered 
  Short_t     fNumChipConfigs;   //number of ChipConfigs used in the chip     
  Short_t     *fDACparams;       //[fNumDACparams] 
  TObjArray   *fChipConfigArray; // array of ChipConfigs in the chip

  private:
  AliITSOnlineSPDfoChip(const AliITSOnlineSPDfoChip &c);
  AliITSOnlineSPDfoChip& operator= (const AliITSOnlineSPDfoChip& c);

    ClassDef(AliITSOnlineSPDfoChip,1)
  };
    
#endif
