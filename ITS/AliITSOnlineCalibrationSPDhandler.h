#ifndef ALI_ITS_ONLINECALIBRATIONSPDHANDLER_H
#define ALI_ITS_ONLINECALIBRATIONSPDHANDLER_H

//////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                           //
// Class that  simplifies the managing of dead and noisy pixels.    //
// Interface to the AliITSOnlineCalibrationSPD container object     //
// through reading and writing to TFile.                            //
//////////////////////////////////////////////////////////////////////  


#include "AliITSIntMap.h"

class TArrayI;

class AliITSOnlineCalibrationSPDhandler {

 public:
  AliITSOnlineCalibrationSPDhandler();
  AliITSOnlineCalibrationSPDhandler(UInt_t module);
  AliITSOnlineCalibrationSPDhandler(const AliITSOnlineCalibrationSPDhandler& handle);
  virtual ~AliITSOnlineCalibrationSPDhandler();
  AliITSOnlineCalibrationSPDhandler& operator=(const AliITSOnlineCalibrationSPDhandler& handle);
  void    ClearMaps();

  void    WriteToFile();
  void    WriteToFile(Char_t* fileName);
  void    ReadFromFile();
  void    ReadDeadFromFile();
  void    ReadNoisyFromFile();
  void    ReadFromFile(Char_t* fileName);
  void    ReadDeadFromFile(Char_t* fileName);
  void    ReadNoisyFromFile(Char_t* fileName);
  void    SetModuleNr(UInt_t mod) {fModuleNr=mod;}
  UInt_t  GetModuleNr() const {return fModuleNr;}
  void    SetFileLocation(Char_t* loc) {sprintf(fFileLocation,"%s",loc);}
  UInt_t  GetNrDead() const {return fDeadPixelMap.GetNrEntries();}
  UInt_t  GetNrNoisy() const {return fNoisyPixelMap.GetNrEntries();}
  TArrayI GetDeadArray();
  TArrayI GetNoisyArray();

  void    ResetDead();
  Bool_t  SetDeadPixel(Int_t col, Int_t row);
  Int_t   GetDeadColAt(UInt_t index) const;
  Int_t   GetDeadRowAt(UInt_t index) const;
  Bool_t  IsPixelDead(Int_t col, Int_t row) const;
  void    ResetNoisy();
  Bool_t  SetNoisyPixel(Int_t col, Int_t row);
  Int_t   GetNoisyColAt(UInt_t index) const;
  Int_t   GetNoisyRowAt(UInt_t index) const;
  Bool_t  IsPixelNoisy(Int_t col, Int_t row) const;

  void    PrintDead() const;
  void    PrintNoisy() const;

 private:
  UInt_t fModuleNr;                 // module nr
  AliITSIntMap fDeadPixelMap;       // list of dead pixels
  AliITSIntMap fNoisyPixelMap;      // list of noisy pixels
  Char_t  fFileLocation[200];       // location (dir) of file to read and write from

  Int_t   GetKey(Int_t col, Int_t row) const {return col*256 + row;}
  Int_t   GetColFromKey(Int_t key) const {return key/256;}
  Int_t   GetRowFromKey(Int_t key) const {return key%256;}

};

#endif
