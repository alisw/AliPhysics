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

            void    SetNTimeBins(Int_t det, Int_t ntb)          { fNTimeBins[det]       = ntb;      }     
            void    SetPretriggerPhase(Int_t det, UInt_t phase) { fPretriggerPhase[det] = phase;    }
            void    SetADCbaseline(Int_t det, Int_t baseline)   { fADCbaseline[det]     = baseline; }

            void    SetNTimeBinsAll(Int_t ntb)                  { for (Int_t i = 0; i < 540; i++) 
                                                                    { fNTimeBins[i]       = ntb;      } }     
            void    SetPretriggerPhaseAll(UInt_t phase)         { for (Int_t i = 0; i < 540; i++)
		                                                    { fPretriggerPhase[i] = phase;    } }
            void    SetADCbaselineAll(Int_t baseline)           { for (Int_t i = 0; i < 540; i++)
		                                                    { fADCbaseline[i]     = baseline; } }

            Int_t   GetNTimeBins(Int_t det) const               { return fNTimeBins[det];           }
            UInt_t  GetPretriggerPhase(Int_t det) const         { return fPretriggerPhase[det];     }
            Int_t   GetADCbaseline(Int_t det) const             { return fADCbaseline[det];         }

  protected:

            Int_t   fNTimeBins[540];       //  Number of timebins for each detector
            UInt_t  fPretriggerPhase[540]; //  Pretrigger phase for each detector
            Int_t   fADCbaseline[540];     //  ADC baseline for each detector, given in ADC channels

    ClassDef(AliTRDdigitsParam,5)          //  The parameters for digits

};
#endif
