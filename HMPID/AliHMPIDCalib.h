#ifndef AliHMPIDCalib_h
#define AliHMPIDCalib_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class of HMPID to manage digits ---> pads
//.
//.
//.

//#include "TTreePlayer.h"
//#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1S.h>
#include <TMath.h>
#include "AliHMPIDParam.h"
#include "AliHMPIDRawStream.h"

class AliHMPIDCalib: public TObject { 


public:
  AliHMPIDCalib();
  virtual ~AliHMPIDCalib();
          void Init();
//          void FillPedestal(Int_t nDDL,Int_t row, Int_t dil,Int_t adr,Int_t q);
          void FillPedestal(Int_t pad,Int_t q);                             //absolute pad number and the charge of the pad
          void FillErrors(Int_t nDDL,Int_t nErrType, Int_t nErr);
        Bool_t CalcPedestal(ULong_t runNum, Int_t nDDL, Char_t* name, Int_t nEv);       //number of the DDL, name of the output file and the number of events processed
        Bool_t WriteErrors(ULong_t runNum, Int_t nDDL, Char_t* name, Int_t nEv);        //number of the DDL, name of the output file and the number of events processed

   enum {
      kNRows       = 24,                                    // Number of rows (starting from 1 !)//was25
      kNDILOGICAdd = 10,                                    // Number of DILOGIC addresses in a row (starting from 1 !) //was11
      kNPadAdd     = 48,                                    // Number of pad row
      kNDDL = 14
    };
        
protected: 
    Bool_t  faddl[11];                                                            //check is ddl is filled
    Float_t fsq[kNDDL][kNRows][kNDILOGICAdd][kNPadAdd];                           //Sum of pad Q
    Float_t fsq2[kNDDL][kNRows][kNDILOGICAdd][kNPadAdd];                          //Sum of pad Q^2
    Int_t   fNumOfErr[kNDDL][AliHMPIDRawStream::kSumErr];                         // Store the numner of errors for a given error type and a given DDL
    ClassDef(AliHMPIDCalib,2)                                                     //HMPID calibration and pedestal class        
};
#endif
