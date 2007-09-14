#ifndef AliHMPIDCalib_h
#define AliHMPIDCalib_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Class of HMPID to manage digits ---> pads
//.
//.
//.
#include "TTreePlayer.h"
#include <TTree.h>
#include <TH1.h>

class AliHMPIDCalib: public TObject { 


public:
  AliHMPIDCalib();
  virtual ~AliHMPIDCalib();
          void Init();
          void FillPedestal(Int_t ddl,Int_t row, Int_t dil,Int_t adr,Int_t q);
        Bool_t CalcPedestal(Int_t ddl, Char_t* name);

        
protected: 
    TTree   *fPedTree;                                                            //Pedestal Tree
    Int_t   fa;                                                                   //DILOGIC address
    Int_t   fd;                                                                   //DILOGIC number
    Int_t   fr;                                                                   //DILOGIC row
    Bool_t  faddl[11];                                                            //check is ddl is filled
    Int_t   fq;                                                                   //Qdc value
    TH1F   *fPedHisto;                                                            //temporary histo for mean and sigma calculation                                                      
 
    ClassDef(AliHMPIDCalib,1)                                                     //HMPID calibration and pedestal class        
};
#endif
