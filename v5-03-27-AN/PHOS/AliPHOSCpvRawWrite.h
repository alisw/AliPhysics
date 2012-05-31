#ifndef ALIPHOSCPVRAWWRITE_H
#define ALIPHOSCPVRAWWRITE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

/* History:
 *
 * $Log$
 */

//_________________________________________________________________________
//  Create a raw data stream for the CPV detector
//  Input:  AliPHOSDigit or TClonesArray of AliPHOSDigits
//  Output: AliFstream, a raw data stream in DDL format
//  Author: Yuri Kharlov
//  14 April 2008
//_________________________________________________________________________

// --- ROOT system ---

class TObject;
class TClonesArray;
class AliFstream;
 
// --- Standard library ---

// --- AliRoot header files ---
class AliPHOSDigit;

class AliPHOSCpvRawWrite : public TObject
{
public:

  AliPHOSCpvRawWrite() ;
  virtual ~AliPHOSCpvRawWrite() ;  
  void WriteRaw(const TObjArray *digits);
  void HWaddress(const AliPHOSDigit *digit, UInt_t &w32, Int_t &ddl, Int_t &row, Int_t &dilogic, Int_t &address);
  void WriteRowMarker (AliFstream *ddl,UInt_t size);
  void WriteSegMarker (AliFstream *ddl,UInt_t row, Int_t nwInSeg); 
  void WriteEoE       (AliFstream *ddl,UInt_t row,UInt_t dil,UInt_t wordCnt); 

 protected:

  Int_t    fNDDL ;          // Number of DDLs
  Int_t    fNRow ;          // Number of row controllers per DDL
  Int_t    fNDilogic ;      // Number of DLOGIC chips per column
  Int_t    fNPad ;          // Number of pads per DLOGIC
  
  ClassDef(AliPHOSCpvRawWrite,1)  // CPV raw data writer

};

#endif // AliPHOSCPVRAWWRITE_H
