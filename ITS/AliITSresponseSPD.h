#ifndef ALIITSRESPONSESPD_H
#define ALIITSRESPONSESPD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliITSresponse.h"
///////////////////////////////////////////
//                                       //
// ITS response class for SPD            //
///////////////////////////////////////////
  
class AliITSresponseSPD :  public AliITSresponse {
 public:
    AliITSresponseSPD(); // default constructor
    virtual ~AliITSresponseSPD() {;} // destructror

    virtual  void   SetCouplingParam(Double_t col, Double_t row)
        {fCouplCol=col; fCouplRow=row;}   
    virtual  void   GetCouplingParam(Double_t &col, Double_t &row) const 
        {col=fCouplCol; row=fCouplRow;}
    virtual void   SetParamOptions(const char* /* a */,const char* /* b */)
      {NotImplemented("SetParamOptions");}
    virtual void   ParamOptions(char *,char*) const
      {NotImplemented("ParamOptions");} 

 protected:

    static const Float_t fgkDiffCoeffDefault;  //default for fDiffCoeff
    static const Float_t fgkCouplingColDefault;//default for fCouplCol 
    static const Float_t fgkCouplingRowDefault;//default for fCouplRow 

    Float_t fCouplCol;        // Coupling probability along a column
    Float_t fCouplRow;        // Coupling probability along a row


    ClassDef(AliITSresponseSPD,4) // SPD base response class
};

#endif
