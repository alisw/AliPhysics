//
// Class AliRsnCutRange
//
// General implementation of cuts which check a value inside a range.
// This range can be defined by two integers or two doubles.
// A user-friendly enumeration allows to define what is checked.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNTOFCORRECTIONESD_H
#define ALIRSNTOFCORRECTIONESD_H

#include "AliESDpid.h"

class AliTOFT0maker;
class AliTOFcalib;
class AliESDEvent;

class AliRsnTOFCorrectionESD : public TObject {
public:

   AliRsnTOFCorrectionESD(Bool_t isMC, Double_t tofRes = 100.0);
   AliRsnTOFCorrectionESD(const AliRsnTOFCorrectionESD& copy);
   AliRsnTOFCorrectionESD& operator=(const AliRsnTOFCorrectionESD& copy);
   virtual ~AliRsnTOFCorrectionESD() { if (fOwnESDpid) delete fESDpid; }

   AliESDpid* ESDpid()                    {return fESDpid;}
   void       SetMC(Bool_t yn = kTRUE)    {fgTOFtuneMC = yn;}
   void       SetESDpid(AliESDpid *pid)   {fESDpid = pid; fOwnESDpid = kFALSE;}
   void       ProcessEvent(AliESDEvent *esd);

private:

   Bool_t                 fOwnESDpid;        //  to know if the object is owned or passed
   AliESDpid             *fESDpid;           //  PID utility for ESD

   //static Bool_t          fgTOFcalibrateESD; //! TOF settings
   static Bool_t          fgTOFcorrectTExp;  //! TOF settings
   static Bool_t          fgTOFuseT0;        //! TOF settings
   static Bool_t          fgTOFtuneMC;       //! TOF settings
   static Double_t        fgTOFresolution;   //! TOF settings
   static AliTOFT0maker  *fgTOFmaker;        //! TOF time0 computator
   static AliTOFcalib    *fgTOFcalib;        //! TOF calibration
   static Int_t           fgLastRun;         //! last run number

   ClassDef(AliRsnTOFCorrectionESD, 1)
};

#endif
