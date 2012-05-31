#ifndef ALIITSONLINESPDFOCHIPCONFIG_H
#define ALIITSONLINESPDFOCHIPCONFIG_H  
/* Copyright(c) 2008-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////////////
// Author: A. Mastroserio                                      //
// This class is the brick of the container for the FastOR     //
// online calibration.                                         //
/////////////////////////////////////////////////////////////////

#include <TObject.h>
class TObjArray;

class AliITSOnlineSPDfoChipConfig : public TObject{
  
 public:
  AliITSOnlineSPDfoChipConfig(); // ctor
  AliITSOnlineSPDfoChipConfig(Short_t measure[4]);
  AliITSOnlineSPDfoChipConfig(const AliITSOnlineSPDfoChipConfig& p);           //  copy constructor
  
  virtual ~AliITSOnlineSPDfoChipConfig(){;} //dtor    
  
  //SETTERS
  void SetChipConfigMatrixId(Short_t idMatrix)    {fMatrixId = idMatrix;}
  void SetChipConfigRow(Short_t pixRow)           {fChipConfigRow = pixRow;}
  void SetChipConfigCol(Short_t pixCol)           {fChipConfigCol = pixCol;}
  void SetChipConfigCounter(UShort_t ncounts)     {fCounter = ncounts;}
  //GETTERS
  void GetChipConfigId(Int_t pixId[2]) const      {pixId[0] = fChipConfigRow; pixId[1]=fChipConfigCol;}
  Short_t  GetChipConfigCounter() const           {return fCounter;}
  Short_t  GetChipConfigMatrixId() const          {return fMatrixId;}
  
  void PrintInfo() {printf(" MatrixId %d -  Row %d   Col %d  -  Counts %d \n",fMatrixId,fChipConfigRow,fChipConfigCol,fCounter);}
  
 private:
  Short_t fChipConfigRow; 
  Short_t fChipConfigCol; 
  Short_t fCounter;
  Short_t fMatrixId;
  AliITSOnlineSPDfoChipConfig &operator=(const AliITSOnlineSPDfoChipConfig& p);
  
  ClassDef(AliITSOnlineSPDfoChipConfig,1)
    
};
    
#endif
