#ifndef ALITRDDIGITSPARAM_H
#define ALITRDDIGITSPARAM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class containing parameters for digits                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/* $Id: AliTRDdigitsParam.h 34070 2009-08-04 15:34:53Z cblume $ */

#include "TObject.h"

class AliTRDdigitsParam : public TObject
{

  public:
    
    AliTRDdigitsParam();
    AliTRDdigitsParam(const AliTRDdigitsParam &p);   
    AliTRDdigitsParam &operator=(const AliTRDdigitsParam &p); 
    virtual        ~AliTRDdigitsParam();

    virtual void    Copy(TObject &p) const;

            void    SetCheckOCDB(Bool_t check = kTRUE)    { fCheckOCDB = check; }
            Bool_t  SetNTimeBins(Int_t ntb);

            Bool_t  CheckOCDB() const                     { return fCheckOCDB;  }
            Int_t   GetNTimeBins() const                  { return fNTimeBins;  }

  protected:

	    Bool_t  fCheckOCDB;          //  Do a consistency check with the corresponding OCDB entry
            Int_t   fNTimeBins;          //  Number of timebins
  
    ClassDef(AliTRDdigitsParam,1)        //  The parameters for digits

};
#endif
