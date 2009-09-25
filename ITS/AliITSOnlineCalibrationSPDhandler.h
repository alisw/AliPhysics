#ifndef ALI_ITS_ONLINECALIBRATIONSPDHANDLER_H
#define ALI_ITS_ONLINECALIBRATIONSPDHANDLER_H

//////////////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                                   //
// Class that  simplifies the managing of dead, noisy, and inactive pixels. //
// Has interface to the AliITSOnlineCalibrationSPD container objects        //
// through reading and writing to TFile.                                    //
//////////////////////////////////////////////////////////////////////////////

/* $Id$  */

#include "AliITSRawStreamSPD.h"
#include <TString.h>

class TArrayI;
class TArrayS;
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
  Bool_t  ReadSilentFromFiles();
  Bool_t  ReadDeadFromFiles();
  Bool_t  ReadNoisyFromFiles();
  Bool_t  ReadSilentFromFile(UInt_t eq);
  Bool_t  ReadDeadFromFile(UInt_t eq);
  Bool_t  ReadNoisyFromFile(UInt_t eq);
  Bool_t  ReadSilentFromFileName(const char *fileName);
  Bool_t  ReadDeadFromFileName(const char *fileName, Bool_t inactive=kFALSE);
  Bool_t  ReadNoisyFromFileName(const char *fileName);

  UInt_t  ReadDeadFromText(const char *fileName, UInt_t module);
  UInt_t  ReadNoisyFromText(const char *fileName, UInt_t module);

  void    WriteToFilesAlways();
  UInt_t  WriteToFiles();
  void    WriteSilentToFilesAlways();
  void    WriteDeadToFilesAlways();
  void    WriteNoisyToFilesAlways();
  UInt_t  WriteSilentToFiles();
  UInt_t  WriteDeadToFiles();
  UInt_t  WriteNoisyToFiles();
  void    WriteSilentToFile(UInt_t eq);
  void    WriteDeadToFile(UInt_t eq, Bool_t inactive=kFALSE);
  void    WriteNoisyToFile(UInt_t eq);

#ifndef SPD_DA_OFF
  Bool_t  ReadDeadModuleFromDB(UInt_t module, Int_t runNr, const Char_t *storage="default", Bool_t treeSerial=kFALSE);
  Bool_t  ReadNoisyModuleFromDB(UInt_t module, Int_t runNr, const Char_t *storage="default", Bool_t treeSerial=kFALSE);
  Bool_t  ReadFromDB(Int_t runNr, const Char_t *storage="default", Bool_t treeSerial=kFALSE);
  Bool_t  ReadDeadFromDB(Int_t runNr, const Char_t *storage="default", Bool_t treeSerial=kFALSE);
  Bool_t  ReadNoisyFromDB(Int_t runNr, const Char_t *storage="default", Bool_t treeSerial=kFALSE);
  Bool_t  ReadDeadFromDBasNoisy(Int_t runNr, const Char_t *storage="default", Bool_t treeSerial=kFALSE);
  Bool_t  ReadDeadFromCalibObj(TObjArray* calObj);
  Bool_t  ReadNoisyFromCalibObj(TObjArray* calObj);
  Bool_t  WriteToDB(Int_t runNrStart, Int_t runNrEnd, const Char_t *storage="default");
  Bool_t  WriteDeadToDB(Int_t runNrStart, Int_t runNrEnd, const Char_t *storage="default");
  Bool_t  WriteDeadToDBasNoisy(Int_t runNrStart, Int_t runNrEnd, const Char_t *storage="default");
  Bool_t  WriteNoisyToDB(Int_t runNrStart, Int_t runNrEnd, const Char_t *storage="default");
#endif

  void    GenerateDCSConfigFile(const Char_t* fileName);

  TArrayS GetSilentArray(UInt_t module, Bool_t treeSerial=kFALSE); // temporarily needed
  TArrayS GetDeadArray(UInt_t module, Bool_t treeSerial=kFALSE);
  TArrayS GetNoisyArray(UInt_t module, Bool_t treeSerial=kFALSE);

  TArrayI GetDeadArrayOnline(UInt_t eq);
  TArrayI GetNoisyArrayOnline(UInt_t eq);

  void    PrintEqSummary();
  void    PrintSilent() const; // silent = dead or inactive
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

  Bool_t  IsPixelBad(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const;
  Bool_t  IsPixelSilent(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const; // silent = dead or inactive
  Bool_t  IsPixelDead(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const;
  Bool_t  IsPixelNoisy(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) const;
  Bool_t  IsPixelBadM(UInt_t module, UInt_t colM, UInt_t row) const;
  Bool_t  IsPixelSilentM(UInt_t module, UInt_t colM, UInt_t row) const;
  Bool_t  IsPixelDeadM(UInt_t module, UInt_t colM, UInt_t row) const;
  Bool_t  IsPixelNoisyM(UInt_t module, UInt_t colM, UInt_t row) const;
  Bool_t  IsPixelBadKey(Int_t key) const;
  Bool_t  IsPixelSilentKey(Int_t key) const;
  Bool_t  IsPixelDeadKey(Int_t key) const;
  Bool_t  IsPixelNoisyKey(Int_t key) const;

  UInt_t  GetNrBad() const; // bad = silent or noisy
  UInt_t  GetNrSilent() const; // silent = dead or inactive
  UInt_t  GetNrDead() const;
  UInt_t  GetDeadEqIdAt(UInt_t index) const;
  UInt_t  GetDeadHSAt(UInt_t index) const;
  UInt_t  GetDeadChipAt(UInt_t index) const;
  UInt_t  GetDeadColAt(UInt_t index) const;
  UInt_t  GetDeadRowAt(UInt_t index) const;
  UInt_t  GetNrNoisy() const;
  UInt_t  GetNoisyEqIdAt(UInt_t index) const;
  UInt_t  GetNoisyHSAt(UInt_t index) const;
  UInt_t  GetNoisyChipAt(UInt_t index) const;
  UInt_t  GetNoisyColAt(UInt_t index) const;
  UInt_t  GetNoisyRowAt(UInt_t index) const;


  UInt_t  GetNrBad(UInt_t module) const; // bad = silent or noisy
  UInt_t  GetNrSilent(UInt_t module) const; // silent = dead or inactive
  UInt_t  GetNrDead(UInt_t module) const;
  UInt_t  GetNrDeadSingle(UInt_t module) const;
  UInt_t  GetDeadEqIdAt(UInt_t module,UInt_t index) const;
  UInt_t  GetDeadHSAt(UInt_t module,UInt_t index) const;
  UInt_t  GetDeadChipAt(UInt_t module,UInt_t index) const;
  UInt_t  GetDeadColAt(UInt_t module,UInt_t index) const;
  UInt_t  GetDeadRowAt(UInt_t module,UInt_t index) const;
  UInt_t  GetNrNoisy(UInt_t module) const;
  UInt_t  GetNrNoisySingle(UInt_t module) const;
  UInt_t  GetNoisyEqIdAt(UInt_t module, UInt_t index) const;
  UInt_t  GetNoisyHSAt(UInt_t module, UInt_t index) const;
  UInt_t  GetNoisyChipAt(UInt_t module, UInt_t index) const;
  UInt_t  GetNoisyColAt(UInt_t module, UInt_t index) const;
  UInt_t  GetNoisyRowAt(UInt_t module, UInt_t index) const;

  UInt_t  GetNrBadEq(UInt_t eq) const; // bad = silent or noisy
  UInt_t  GetNrSilentEq(UInt_t eq) const; // silent = dead or inactive
  UInt_t  GetNrDeadEq(UInt_t eq) const;
  UInt_t  GetDeadEqIdAtEq(UInt_t eq, UInt_t index) const;
  UInt_t  GetDeadHSAtEq(UInt_t eq, UInt_t index) const;
  UInt_t  GetDeadChipAtEq(UInt_t eq, UInt_t index) const;
  UInt_t  GetDeadColAtEq(UInt_t eq, UInt_t index) const;
  UInt_t  GetDeadRowAtEq(UInt_t eq, UInt_t index) const;
  UInt_t  GetNrNoisyEq(UInt_t eq) const;
  UInt_t  GetNoisyEqIdAtEq(UInt_t eq, UInt_t index) const;
  UInt_t  GetNoisyHSAtEq(UInt_t eq, UInt_t index) const;
  UInt_t  GetNoisyChipAtEq(UInt_t eq, UInt_t index) const;
  UInt_t  GetNoisyColAtEq(UInt_t eq, UInt_t index) const;
  UInt_t  GetNoisyRowAtEq(UInt_t eq, UInt_t index) const;

  UInt_t  GetNrBadC(UInt_t eq, UInt_t hs, UInt_t chip) const; // bad = silent or noisy
  UInt_t  GetNrSilentC(UInt_t eq, UInt_t hs, UInt_t chip) const; // silent = dead or inactive
  UInt_t  GetNrDeadC(UInt_t eq, UInt_t hs, UInt_t chip) const;
  UInt_t  GetDeadEqIdAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetDeadHSAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetDeadChipAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetDeadColAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetDeadRowAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetNrNoisyC(UInt_t eq, UInt_t hs, UInt_t chip) const;
  UInt_t  GetNoisyEqIdAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetNoisyHSAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetNoisyChipAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetNoisyColAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  UInt_t  GetNoisyRowAtC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;

  const Char_t* GetDeadPixelAsTextC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;
  const Char_t* GetNoisyPixelAsTextC(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t index) const;

  void    ActivateALL();
  void    ActivateEq(UInt_t eq, Bool_t setval = kTRUE);
  void    ActivateHS(UInt_t eq, UInt_t hs, Bool_t setval = kTRUE);
  void    ActivateChip(UInt_t eq, UInt_t hs, UInt_t chip, Bool_t setval = kTRUE);

  void    UnSetDeadALL();
  void    SetDeadEq(UInt_t eq, Bool_t setval = kTRUE);
  void    SetDeadHS(UInt_t eq, UInt_t hs, Bool_t setval = kTRUE);
  void    SetDeadChip(UInt_t eq, UInt_t hs, UInt_t chip, Bool_t setval = kTRUE);

  Bool_t  IsActiveEq(UInt_t eq) const;
  Bool_t  IsActiveHS(UInt_t eq, UInt_t hs) const;
  Bool_t  IsActiveChip(UInt_t eq, UInt_t hs, UInt_t chip) const;

  Bool_t  IsDeadEq(UInt_t eq) const;
  Bool_t  IsDeadHS(UInt_t eq, UInt_t hs) const;
  Bool_t  IsDeadChip(UInt_t eq, UInt_t hs, UInt_t chip) const;

  Bool_t  IsSilentEq(UInt_t eq) const;
  Bool_t  IsSilentHS(UInt_t eq, UInt_t hs) const;
  Bool_t  IsSilentChip(UInt_t eq, UInt_t hs, UInt_t chip) const;

  Bool_t  IsNoisyChip(UInt_t eq, UInt_t hs, UInt_t chip) const; 

  UInt_t  AddSilentFrom(AliITSOnlineCalibrationSPDhandler* other);
  UInt_t  AddDeadFrom(AliITSOnlineCalibrationSPDhandler* other);
  UInt_t  AddNoisyFrom(AliITSOnlineCalibrationSPDhandler* other);

  UInt_t  GetNrDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  UInt_t  GetNrSilentDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  UInt_t  GetNrDeadDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  UInt_t  GetNrNoisyDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  AliITSOnlineCalibrationSPDhandler* GetDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  AliITSOnlineCalibrationSPDhandler* GetSilentDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  AliITSOnlineCalibrationSPDhandler* GetDeadDiff(AliITSOnlineCalibrationSPDhandler* other) const;
  AliITSOnlineCalibrationSPDhandler* GetNoisyDiff(AliITSOnlineCalibrationSPDhandler* other) const;


 private:
  TString fFileLocation;              // location (dir) of files to read and write from
  AliITSIntMap* fDeadPixelMap[1200];  // lists of dead pixels for each chip
  AliITSIntMap* fNoisyPixelMap[1200]; // lists of noisy pixels for each chip
  UInt_t fNrDead[1200];               // nr of dead pixels for each chip
  UInt_t fNrNoisy[1200];              // nr of noisy pixels for each chip
  Bool_t fActiveEq[20];               // active bit for each equipment
  Bool_t fActiveHS[20][6];            // active bit for each half-stave
  Bool_t fActiveChip[20][6][10];      // active bit for each chip
  Bool_t fDeadEq[20];                 // dead bit for each equipment
  Bool_t fDeadHS[20][6];              // dead bit for each half-stave
  Bool_t fDeadChip[20][6][10];        // dead bit for each chip

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
  
  UInt_t   GetGloChip(UInt_t eq, UInt_t hs, UInt_t chip) const {return eq*60 + hs*10 + chip;}
  UInt_t   GetEqGlo(UInt_t gloChip) const {return gloChip/60;}
  UInt_t   GetHSGlo(UInt_t gloChip) const {return (gloChip%60)/10;}
  UInt_t   GetChipGlo(UInt_t gloChip) const {return (gloChip%60)%10;}

  void     GetChipAndIndexDead(UInt_t module, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) const;
  void     GetChipAndIndexNoisy(UInt_t module, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) const;
  void     GetChipAndIndexEqDead(UInt_t eq, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) const;
  void     GetChipAndIndexEqNoisy(UInt_t eq, UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) const;
  void     GetChipAndIndexTotDead(UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) const;
  void     GetChipAndIndexTotNoisy(UInt_t index, UInt_t& gloChip, UInt_t& chipIndex) const;

  UInt_t   GetEqIdFromOffline(UInt_t module) const;
  UInt_t   GetHSFromOffline(UInt_t module) const;
  UInt_t   GetChipFromOffline(UInt_t module, UInt_t colM) const;
  UInt_t   GetColFromOffline(UInt_t module, UInt_t colM) const;
  UInt_t   GetRowFromOffline(UInt_t module, UInt_t rowM) const;

  void     RecursiveInsertDead(AliITSCalibrationSPD* calibSPD, UInt_t module, Int_t lowInd, Int_t highInd);
  void     RecursiveInsertNoisy(AliITSCalibrationSPD* calibSPD, UInt_t module, Int_t lowInd, Int_t highInd);

  UInt_t   GetNrDeadC2(UInt_t gloChip) const;
  UInt_t   GetDeadEqIdAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetDeadHSAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetDeadChipAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetDeadColAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetDeadRowAtC2(UInt_t gloChip, UInt_t index) const;	  
  UInt_t   GetNrNoisyC2(UInt_t gloChip) const;
  UInt_t   GetNoisyEqIdAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetNoisyHSAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetNoisyChipAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetNoisyColAtC2(UInt_t gloChip, UInt_t index) const;
  UInt_t   GetNoisyRowAtC2(UInt_t gloChip, UInt_t index) const;

};

#endif
