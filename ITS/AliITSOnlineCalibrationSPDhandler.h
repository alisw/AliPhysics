#ifndef ALI_ITS_ONLINECALIBRATIONSPDHANDLER_H
#define ALI_ITS_ONLINECALIBRATIONSPDHANDLER_H

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
class AliITSCalibrationSPD;

class AliITSOnlineCalibrationSPDhandler {

 public:
  AliITSOnlineCalibrationSPDhandler();
  AliITSOnlineCalibrationSPDhandler(const AliITSOnlineCalibrationSPDhandler& handle);
  virtual ~AliITSOnlineCalibrationSPDhandler();
  AliITSOnlineCalibrationSPDhandler& operator=(const AliITSOnlineCalibrationSPDhandler& handle);

  void    SetFileLocation(const Char_t* loc) {fFileLocation = loc;}
  TString GetFileLocation() const {return fFileLocation;}

  void    ClearMaps();
  void    ResetDead();
  void    ResetNoisy();
  void    ResetDeadForChip(UInt_t eq, UInt_t hs, UInt_t chip);
  void    ResetNoisyForChip(UInt_t eq, UInt_t hs, UInt_t chip);
  void    ResetDeadForEq(UInt_t eq);
  void    ResetNoisyForEq(UInt_t eq);

  Bool_t  ReadFromFiles();
  Bool_t  ReadDeadFromFiles();
  Bool_t  ReadNoisyFromFiles();
  Bool_t  ReadDeadFromFile(UInt_t eq);
  Bool_t  ReadNoisyFromFile(UInt_t eq);
  Bool_t  ReadDeadFromFileName(const char *fileName);
  Bool_t  ReadNoisyFromFileName(const char *fileName);
  UInt_t  ReadDeadFromText(const char *fileName, UInt_t module);
  UInt_t  ReadNoisyFromText(const char *fileName, UInt_t module);

  void    WriteToFilesAlways();
  UInt_t  WriteToFiles();
  void    WriteDeadToFilesAlways();
  void    WriteNoisyToFilesAlways();
  UInt_t  WriteDeadToFiles();
  UInt_t  WriteNoisyToFiles();
  void    WriteDeadToFile(UInt_t eq);
  void    WriteNoisyToFile(UInt_t eq);

#ifndef SPD_DA_OFF
  Bool_t  ReadDeadModuleFromDB(UInt_t module, Int_t runNr, Bool_t treeSerial=kFALSE);
  Bool_t  ReadNoisyModuleFromDB(UInt_t module, Int_t runNr, Bool_t treeSerial=kFALSE);
  Bool_t  ReadFromDB(Int_t runNr, Bool_t treeSerial=kFALSE);
  Bool_t  ReadDeadFromDB(Int_t runNr, Bool_t treeSerial=kFALSE);
  Bool_t  ReadNoisyFromDB(Int_t runNr, Bool_t treeSerial=kFALSE);
  Bool_t  ReadDeadFromCalibObj(TObjArray* calObj);
  Bool_t  ReadNoisyFromCalibObj(TObjArray* calObj);
  Bool_t  WriteToDB(Int_t runNrStart, Int_t runNrEnd);
  Bool_t  WriteDeadToDB(Int_t runNrStart, Int_t runNrEnd);
  Bool_t  WriteNoisyToDB(Int_t runNrStart, Int_t runNrEnd);
#endif

  void    GenerateDCSConfigFile(const Char_t* fileName);

  TArrayI GetDeadArray(UInt_t module, Bool_t treeSerial=kFALSE);
  TArrayI GetNoisyArray(UInt_t module, Bool_t treeSerial=kFALSE);
  TArrayI GetDeadArrayOnline(UInt_t eq);
  TArrayI GetNoisyArrayOnline(UInt_t eq);

  void    PrintEqSummary();
  void    PrintDead() const;
  void    PrintNoisy() const;

  Bool_t  SetDeadPixel(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row);
  Bool_t  SetNoisyPixel(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row);
  Bool_t  SetDeadPixelM(UInt_t module, UInt_t colM, UInt_t row);
  Bool_t  SetNoisyPixelM(UInt_t module, UInt_t colM, UInt_t row);
  Bool_t  UnSetDeadPixel(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row);
  Bool_t  UnSetNoisyPixel(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row);
  Bool_t  UnSetDeadPixelM(UInt_t module, UInt_t colM, UInt_t row);
  Bool_t  UnSetNoisyPixelM(UInt_t module, UInt_t colM, UInt_t row);
  UInt_t  SetDeadChip(UInt_t eq, UInt_t hs, UInt_t chip);
  UInt_t  SetNoisyChip(UInt_t eq, UInt_t hs, UInt_t chip);
  Bool_t  UnSetDeadChip(UInt_t eq, UInt_t hs, UInt_t chip);
  Bool_t  UnSetNoisyChip(UInt_t eq, UInt_t hs, UInt_t chip);

  Bool_t  IsPixelDead(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const;
  Bool_t  IsPixelNoisy(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const;
  Bool_t  IsPixelDeadM(UInt_t module, UInt_t colM, UInt_t row) const;
  Bool_t  IsPixelNoisyM(UInt_t module, UInt_t colM, UInt_t row) const;
  Bool_t  IsPixelDeadKey(Int_t key) const;
  Bool_t  IsPixelNoisyKey(Int_t key) const;

  UInt_t  GetNrDead() const;
  UInt_t  GetNrNoisy() const;
  UInt_t  GetDeadEqIdAt(UInt_t index);
  UInt_t  GetNoisyEqIdAt(UInt_t index);
  UInt_t  GetDeadHSAt(UInt_t index);
  UInt_t  GetNoisyHSAt(UInt_t index);
  UInt_t  GetDeadChipAt(UInt_t index);
  UInt_t  GetNoisyChipAt(UInt_t index);
  UInt_t  GetDeadColAt(UInt_t index);
  UInt_t  GetNoisyColAt(UInt_t index);
  UInt_t  GetDeadRowAt(UInt_t index);
  UInt_t  GetNoisyRowAt(UInt_t index);

  UInt_t  GetNrDead(UInt_t module) const;
  UInt_t  GetNrNoisy(UInt_t module) const;
  UInt_t  GetDeadEqIdAt(UInt_t module,UInt_t index);
  UInt_t  GetNoisyEqIdAt(UInt_t module, UInt_t index);
  UInt_t  GetDeadHSAt(UInt_t module,UInt_t index);
  UInt_t  GetNoisyHSAt(UInt_t module, UInt_t index);
  UInt_t  GetDeadChipAt(UInt_t module,UInt_t index);
  UInt_t  GetNoisyChipAt(UInt_t module, UInt_t index);
  UInt_t  GetDeadColAt(UInt_t module,UInt_t index);
  UInt_t  GetNoisyColAt(UInt_t module, UInt_t index);
  UInt_t  GetDeadRowAt(UInt_t module,UInt_t index);
  UInt_t  GetNoisyRowAt(UInt_t module, UInt_t index);

  UInt_t  GetNrDeadEq(UInt_t eq) const;
  UInt_t  GetNrNoisyEq(UInt_t eq) const;
  UInt_t  GetDeadEqIdAtEq(UInt_t eq, UInt_t index) const;
  UInt_t  GetNoisyEqIdAtEq(UInt_t eq, UInt_t index) const;
  UInt_t  GetDeadHSAtEq(UInt_t eq, UInt_t index);
  UInt_t  GetNoisyHSAtEq(UInt_t eq, UInt_t index);
  UInt_t  GetDeadChipAtEq(UInt_t eq, UInt_t index);
  UInt_t  GetNoisyChipAtEq(UInt_t eq, UInt_t index);
  UInt_t  GetDeadColAtEq(UInt_t eq, UInt_t index);
  UInt_t  GetNoisyColAtEq(UInt_t eq, UInt_t index);
  UInt_t  GetDeadRowAtEq(UInt_t eq, UInt_t index);
  UInt_t  GetNoisyRowAtEq(UInt_t eq, UInt_t index);

  UInt_t  GetNrDeadC(UInt_t eq, UInt_t hs, UInt_t chip) const;
  UInt_t  GetNrNoisyC(UInt_t eq, UInt_t hs, UInt_t chip) const;
  UInt_t  GetDeadEqIdAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetNoisyEqIdAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetDeadHSAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetNoisyHSAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetDeadChipAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetNoisyChipAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetDeadColAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetNoisyColAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetDeadRowAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetNoisyRowAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;

  const Char_t* GetDeadPixelAsTextC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  const Char_t* GetNoisyPixelAsTextC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;

  UInt_t  AddDeadFrom(AliITSOnlineCalibrationSPDhandler* other);
  UInt_t  AddNoisyFrom(AliITSOnlineCalibrationSPDhandler* other);

  UInt_t  GetNrDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  UInt_t  GetNrDeadDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  UInt_t  GetNrNoisyDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  AliITSOnlineCalibrationSPDhandler* GetDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  AliITSOnlineCalibrationSPDhandler* GetDeadDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  AliITSOnlineCalibrationSPDhandler* GetNoisyDiff(AliITSOnlineCalibrationSPDhandler* other) const;


 private:
  TString fFileLocation;              // location (dir) of files to read and write from
  AliITSIntMap* fDeadPixelMap[1200];  // lists of dead pixels for each chip
  AliITSIntMap* fNoisyPixelMap[1200]; // lists of noisy pixels for each chip
  UInt_t fNrDead[1200];               // nr of dead pixels for each chip
  UInt_t fNrNoisy[1200];              // nr of noisy pixels for each chip
  				     
  Int_t    GetKey(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const 
    {return eq*6*10*32*256 + hs*10*32*256 + chip*32*256 + col*256 + row;}
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
    {return AliITSRawStreamSPD::GetOfflineModuleFromOnline(GetEqIdFromKey(key),GetHSFromKey(key),GetChipFromKey(key));}
  UInt_t   GetColMFromKey(Int_t key) const 
    {return AliITSRawStreamSPD::GetOfflineColFromOnline(GetEqIdFromKey(key),GetHSFromKey(key),GetChipFromKey(key),GetColFromKey(key));}
  UInt_t   GetRowMFromKey(Int_t key) const 
    {return AliITSRawStreamSPD::GetOfflineRowFromOnline(GetEqIdFromKey(key),GetHSFromKey(key),GetChipFromKey(key),GetRowFromKey(key));}
  
  UInt_t   GetGloChip(UInt_t eq, UInt_t hs, UInt_t chip) const
    {return eq*6*10 + hs*10 + chip;}

  void     GetChipAndIndexDead(UInt_t module, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex);
  void     GetChipAndIndexNoisy(UInt_t module, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex);
  void     GetChipAndIndexEqDead(UInt_t eq, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex);
  void     GetChipAndIndexEqNoisy(UInt_t eq, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex);
  void     GetChipAndIndexTotDead(UInt_t index, UInt_t& gloChip, UInt_t& chipIndex);
  void     GetChipAndIndexTotNoisy(UInt_t index, UInt_t& gloChip, UInt_t& chipIndex);

  UInt_t   GetEqIdFromOffline(UInt_t module) const;
  UInt_t   GetHSFromOffline(UInt_t module) const;
  UInt_t   GetChipFromOffline(UInt_t module, UInt_t colM) const;
  UInt_t   GetColFromOffline(UInt_t module, UInt_t colM) const;
  UInt_t   GetRowFromOffline(UInt_t module, UInt_t rowM) const;

  void     RecursiveInsertDead(AliITSCalibrationSPD* calibSPD, UInt_t module, Int_t lowInd, Int_t highInd);
  void     RecursiveInsertNoisy(AliITSCalibrationSPD* calibSPD, UInt_t module, Int_t lowInd, Int_t highInd);

  UInt_t   GetNrDeadC2(UInt_t gloChip) const;
  UInt_t   GetNrNoisyC2(UInt_t gloChip) const;
  UInt_t   GetDeadEqIdAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetNoisyEqIdAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetDeadHSAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetNoisyHSAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetDeadChipAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetNoisyChipAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetDeadColAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetNoisyColAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetDeadRowAtC2(UInt_t gloChip, UInt_t index) const;	  
  UInt_t   GetNoisyRowAtC2(UInt_t gloChip, UInt_t index) const;

};

#endif
