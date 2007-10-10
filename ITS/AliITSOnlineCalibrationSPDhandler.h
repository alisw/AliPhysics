#ifndef ALI_ITS_ONLINECALIBRATIONSPDHANDLER_H
#define ALI_ITS_ONLINECALIBRATIONSPDHANDLER_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                               //
// Class that  simplifies the managing of dead and noisy pixels.        //
// Has interface to the AliITSOnlineCalibrationSPD container objects    //
// through reading and writing to TFile.                                //
//////////////////////////////////////////////////////////////////////////


#include "AliITSRawStreamSPD.h"
#include <TString.h>

class TArrayI;
class AliITSIntMap;

class AliITSOnlineCalibrationSPDhandler {

 public:
  AliITSOnlineCalibrationSPDhandler();
  AliITSOnlineCalibrationSPDhandler(const AliITSOnlineCalibrationSPDhandler& handle);
  virtual ~AliITSOnlineCalibrationSPDhandler();
  AliITSOnlineCalibrationSPDhandler& operator=(const AliITSOnlineCalibrationSPDhandler& handle);

  void    SetFileLocation(const Char_t* loc) {fFileLocation = loc;}
  TString GetFileLocation() const {return fFileLocation;}

  void    ClearMaps();

  Bool_t  ReadFromFiles();
  Bool_t  ReadDeadFromFiles();
  Bool_t  ReadNoisyFromFiles();
  Bool_t  ReadFromFile(UInt_t module);
  Bool_t  ReadDeadFromFile(UInt_t module);
  Bool_t  ReadNoisyFromFile(UInt_t module);
  Bool_t  ReadFromFileName(const char *fileName);
  Bool_t  ReadDeadFromFileName(const char *fileName);
  Bool_t  ReadNoisyFromFileName(const char *fileName);

  void    WriteToFilesAlways();
  void    WriteToFiles();
  void    WriteDeadToFiles();
  void    WriteNoisyToFiles();
  void    WriteToFile(UInt_t module);  
  void    WriteDeadToFile(UInt_t module);
  void    WriteNoisyToFile(UInt_t module);

#ifndef SPD_DA_OFF
  Bool_t  ReadModuleFromDB(UInt_t module, Int_t runNr);
  Bool_t  ReadFromDB(Int_t runNr);
  Bool_t  WriteToDB(Int_t runNrStart, Int_t runNrEnd);
#endif
  UInt_t  ReadNoisyFromText(const char *fileName);
  UInt_t  ReadDeadFromText(const char *fileName);
  void    GenerateDCSConfigFile(const Char_t* fileName);

  TArrayI GetDeadArray(UInt_t module);
  TArrayI GetNoisyArray(UInt_t module);

  void    ResetDead();
  void    ResetDeadForChip(UInt_t eqId, UInt_t hs, UInt_t chip);
  Bool_t  SetDeadPixel(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row);
  Bool_t  SetDeadPixelM(UInt_t module, UInt_t colM, UInt_t row);
  Bool_t  UnSetDeadPixel(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row);
  Bool_t  UnSetDeadPixelM(UInt_t module, UInt_t colM, UInt_t row);
  Bool_t  IsPixelDead(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const;
  Bool_t  IsPixelDeadM(UInt_t module, UInt_t colM, UInt_t row);

  UInt_t  GetNrDead(UInt_t module) const;
  UInt_t  GetDeadEqIdAt(UInt_t module,UInt_t index);
  UInt_t  GetDeadHSAt(UInt_t module,UInt_t index);
  UInt_t  GetDeadChipAt(UInt_t module,UInt_t index);
  UInt_t  GetDeadColAt(UInt_t module,UInt_t index);
  UInt_t  GetDeadRowAt(UInt_t module,UInt_t index);

  void    ResetNoisy();
  void    ResetNoisyForChip(UInt_t eqId, UInt_t hs, UInt_t chip);
  Bool_t  SetNoisyPixel(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row);
  Bool_t  SetNoisyPixelM(UInt_t module, UInt_t colM, UInt_t row);
  Bool_t  UnSetNoisyPixel(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row);
  Bool_t  UnSetNoisyPixelM(UInt_t module, UInt_t colM, UInt_t row);
  Bool_t  IsPixelNoisy(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const;
  Bool_t  IsPixelNoisyM(UInt_t module, UInt_t colM, UInt_t row);

  UInt_t  GetNrNoisy(UInt_t module) const;
  UInt_t  GetNoisyEqIdAt(UInt_t module, UInt_t index);
  UInt_t  GetNoisyHSAt(UInt_t module, UInt_t index);
  UInt_t  GetNoisyChipAt(UInt_t module, UInt_t index);
  UInt_t  GetNoisyColAt(UInt_t module, UInt_t index);
  UInt_t  GetNoisyRowAt(UInt_t module, UInt_t index);

  UInt_t  GetNrDead() const;
  UInt_t  GetNrNoisy() const;
  void    PrintDead() const;
  void    PrintNoisy() const;

  UInt_t  GetNrDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  UInt_t  GetNrDeadDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  UInt_t  GetNrNoisyDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  AliITSOnlineCalibrationSPDhandler* GetDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  AliITSOnlineCalibrationSPDhandler* GetDeadDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  AliITSOnlineCalibrationSPDhandler* GetNoisyDiff(AliITSOnlineCalibrationSPDhandler* other) const;


 private:
  TString fFileLocation;             // location (dir) of files to read and write from
  AliITSIntMap* fDeadPixelMap[240];  // lists of dead pixels for each module
  AliITSIntMap* fNoisyPixelMap[240]; // lists of noisy pixels for each module
  UInt_t fNrDead[240];               // nr of dead pixels for each module
  UInt_t fNrNoisy[240];              // nr of noisy pixels for each module
  				     
  Bool_t fModuleMapInited;           // flag to know if arrays below are filled 
  UInt_t fiDDL[240];                 // iDDL value for each module (inited when used, fModuleMapInited flag)
  UInt_t fiModule[240];              // iModule value for each module (inited when used, fModuleMapInited flag)

  Int_t    GetKey(UInt_t eqId, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const 
    {return eqId*6*10*32*256 + hs*10*32*256 + chip*32*256 + col*256 + row;}

  UInt_t   GetEqIdFromKey(Int_t key) const 
    {return key/(6*10*32*256);}
  UInt_t   GetHSFromKey(Int_t key) const 
    {return (key%(6*10*32*256))/(10*32*256);}
  UInt_t   GetChipFromKey(Int_t key) const 
    {return ((key%(6*10*32*256))%(10*32*256))/(32*256);}
  UInt_t   GetColFromKey(Int_t key) const 
    {return (((key%(6*10*32*256))%(10*32*256))%(32*256))/256;}
  UInt_t   GetRowFromKey(Int_t key) const 
    {return (((key%(6*10*32*256))%(10*32*256))%(32*256))%256;}
  UInt_t   GetModuleFromKey(Int_t key) const
    {return AliITSRawStreamSPD::GetModuleNumber(GetEqIdFromKey(key),GetHSFromKey(key),GetChipFromKey(key));}
  UInt_t   GetColMFromKey(Int_t key) const 
    {return GetColFromKey(key) + 32 * (GetChipFromKey(key) % 5);}

  void     InitModuleMaps();
  UInt_t   GetEqIdFromOffline(UInt_t module);
  UInt_t   GetHSFromOffline(UInt_t module);
  UInt_t   GetChipFromOffline(UInt_t module, UInt_t colM);
  UInt_t   GetColFromOffline(UInt_t colM) const;

  Bool_t   IsPixelDeadKey(Int_t key) const;
  Bool_t   IsPixelDeadMKey(UInt_t module, Int_t key) const;
  Bool_t   IsPixelNoisyKey(Int_t key) const;
  Bool_t   IsPixelNoisyMKey(UInt_t module, Int_t key) const;

};

#endif
