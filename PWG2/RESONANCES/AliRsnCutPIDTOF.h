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

#ifndef ALIRSNCUTPIDTOF_H
#define ALIRSNCUTPIDTOF_H

#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODpidUtil.h"

class AliTOFT0maker;
class AliTOFcalib;
class AliESDtrack;
class AliAODTrack;

#include "AliRsnCut.h"

class AliRsnCutPIDTOF : public AliRsnCut
{
  public:

    AliRsnCutPIDTOF(const char *name = "cutPIDTOF", Bool_t isMC = kFALSE, Double_t min = -10.0, Double_t max = 10.0);
    AliRsnCutPIDTOF(const AliRsnCutPIDTOF& copy);
    AliRsnCutPIDTOF& operator=(const AliRsnCutPIDTOF& copy);
    virtual ~AliRsnCutPIDTOF() { }
    
    void             SetMC(Bool_t yn = kTRUE) {fIsMC = yn;}
    virtual Bool_t   IsSelected(TObject *object);

  protected:
  
    void    ProcessCurrentEvent();
    Bool_t  CheckESD(AliESDtrack *track);
    Bool_t  CheckAOD(AliAODTrack *track);
  
    Bool_t                 fIsMC;             //  switch for MC analysis    
    AliPID::EParticleType  fPIDtype;          //  particle type for which PID is checked   
    AliESDpid              fESDpid;           //  PID utility for ESD
    AliAODpidUtil          fAODpid;           //  PID utility for AOD
    
  //static Bool_t          fgTOFcalibrateESD; //! TOF settings
    static Bool_t          fgTOFcorrectTExp;  //! TOF settings
    static Bool_t          fgTOFuseT0;        //! TOF settings
    static Bool_t          fgTOFtuneMC;       //! TOF settings
    static Double_t        fgTOFresolution;   //! TOF settings
    static AliTOFT0maker  *fgTOFmaker;        //! TOF time0 computator
    static AliTOFcalib    *fgTOFcalib;        //! TOF calibration
    static Int_t           fgLastRun;         //! last run number
    static Int_t           fgLastEventID;     //! ID of last event processed
    static AliESDEvent    *fgLastEvent;       //! pointer to last processed event

    ClassDef(AliRsnCutPIDTOF, 1)
};

#endif
