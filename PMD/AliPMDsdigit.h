#ifndef ALIPMDSDIGIT_H
#define ALIPMDSDIGIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//                                                     //
//  Date   : August 05 2003                            //
//  used to store the info into TreeS                  //
//                                                     //
//-----------------------------------------------------//
// Author - B.K. Nandi
//
#include "TObject.h"
class TClonesArray;

class AliPMDsdigit : public TObject
{

 public:
  AliPMDsdigit();
  AliPMDsdigit(Int_t trnumber, Int_t trpid, Int_t det, Int_t smn,
	       Int_t irow, Int_t icol, Float_t edep);
  AliPMDsdigit(AliPMDsdigit *pmdsdigit);
  AliPMDsdigit (const AliPMDsdigit &pmdsdigit);  // copy constructor
  AliPMDsdigit &operator=(const AliPMDsdigit &pmdsdigit); // assignment op

  virtual ~AliPMDsdigit();

  Int_t   GetTrackNumber() const;
  Int_t   GetTrackPid() const;
  Int_t   GetDetector() const;
  Int_t   GetSMNumber() const;
  Int_t   GetRow() const;
  Int_t   GetColumn() const;
  Float_t GetCellEdep() const;

  
 protected:
  Int_t   fTrNumber;   // Parent Track Number
  Int_t   fTrPid;      // Parent Track pid
  Int_t   fDet;        // Detector Number (0:PRE, 1:CPV)
  Int_t   fSMN;        // Serial Module Number
  Int_t   fRow;        // Cell Row Number (0-47)
  Int_t   fColumn;     // Cell Column Number (0-95)
  Float_t fEdep;       // Energy deposition in a hexagonal cell
  
  ClassDef(AliPMDsdigit,5) // SDigits object for Detector set:PMD
};

#endif
