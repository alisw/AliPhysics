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
#include <TMath.h>

class AliHMPIDCalib: public TObject { 


public:
  AliHMPIDCalib();
  virtual ~AliHMPIDCalib();
          void Init();
          void FillPedestal(Int_t nDDL,Int_t row, Int_t dil,Int_t adr,Int_t q);
        Bool_t CalcPedestal(Int_t nDDL, Char_t* name, Int_t nEv);

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
    ClassDef(AliHMPIDCalib,1)                                                     //HMPID calibration and pedestal class        
};
#endif
