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

#ifndef ALIRSNCUTESD2010_H
#define ALIRSNCUTESD2010_H

#include "AliESDtrackCuts.h"
#include "AliRsnCut.h"

class AliESDpid;
class AliTOFT0maker;
class AliTOFcalib;

class AliRsnCutESD2010 : public AliRsnCut
{
  public:

    AliRsnCutESD2010(const char *name = "cutESD2010", Bool_t isMC = kFALSE);
    AliRsnCutESD2010(const AliRsnCutESD2010& copy);
    AliRsnCutESD2010& operator=(const AliRsnCutESD2010& copy);
    virtual ~AliRsnCutESD2010() {;};

    AliESDtrackCuts* GetCutsTPC() {return &fESDtrackCutsTPC;}
    AliESDtrackCuts* GetCutsITS() {return &fESDtrackCutsITS;}
    virtual Bool_t   IsSelected(TObject *obj1, TObject *obj2 = 0x0);
    
    void             SetMC       (Bool_t yn = kTRUE);
    void             SetCheckITS (Bool_t yn = kTRUE) {fCheckITS = yn;}
    void             SetCheckTPC (Bool_t yn = kTRUE) {fCheckTPC = yn;}
    void             SetCheckTOF (Bool_t yn = kTRUE) {fCheckTOF = yn;}
    void             SetUseGlobal(Bool_t yn = kTRUE) {fUseGlobal = yn;}
    void             SetUseITSSA (Bool_t yn = kTRUE) {fUseITSSA = yn;}
    void             SetMaxEta   (Double_t v)        {fMaxEta = v;}
    
    void             SetITSband(Double_t v) {fMaxITSband = v;}
    
    void             SetTPCpLimit(Double_t v) {fTPCpLimit = v;}
    void             SetTPCrange(Double_t min, Double_t max) {fMinTPCband = min; fMaxTPCband = max;}
    void             SetTPCpar(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4)
                       {fTPCpar[0]=p0;fTPCpar[1]=p1;fTPCpar[2]=p2;fTPCpar[3]=p3;fTPCpar[4]=p4;}

    void             SetTOFcalibrateESD(Bool_t yn = kTRUE)  {fTOFcalibrateESD = yn;}
    void             SetTOFcorrectTExp (Bool_t yn = kTRUE)  {fTOFcorrectTExp = yn;}
    void             SetTOFuseT0       (Bool_t yn = kTRUE)  {fTOFuseT0 = yn;}
    void             SetTOFtuneMC      (Bool_t yn = kTRUE)  {fTOFtuneMC = yn;}
    void             SetTOFresolution  (Double_t v = 100.0) {fTOFresolution = v;}
    void             SetTOFrange       (Double_t v1, Double_t v2) {fMinTOF = v1; fMaxTOF = v2;}
    
    virtual void     SetEvent(AliRsnEvent *event);

  protected:
  
    Bool_t           fIsMC;             //  switch for MC analysis
    Bool_t           fCheckITS;         //  switch for ITS dE/dx check
    Bool_t           fCheckTPC;         //  switch for TPC dE/dx check
    Bool_t           fCheckTOF;         //  switch for TOF time check
    Bool_t           fUseGlobal;        //  switch to use TPC global tracks
    Bool_t           fUseITSSA;         //  switch to use ITS standalone tracks
    
    Double_t         fMaxEta;           //  cut in eta

    Double_t         fMaxITSband;       //  range for ITS de/dx band

    Double_t         fTPCpLimit;        //  limit to choose what band to apply
    Double_t         fTPCpar[5];        //  parameters for TPC bethe-Bloch
    Double_t         fMinTPCband;       //  range for TPC de/dx band - min
    Double_t         fMaxTPCband;       //  range for TPC de/dx band - max
    
    AliESDtrackCuts  fESDtrackCutsTPC;  //  ESD standard defined track cuts for TPC tracks
    AliESDtrackCuts  fESDtrackCutsITS;  //  ESD standard defined track cuts for ITS-SA tracks
    AliESDpid       *fESDpid;           //! PID manager
    AliTOFT0maker   *fTOFmaker;         //! TOF time0 computator
    AliTOFcalib     *fTOFcalib;         //! TOF calibration
    Bool_t           fTOFcalibrateESD;  //  TOF settings
    Bool_t           fTOFcorrectTExp;   //  TOF settings
    Bool_t           fTOFuseT0;         //  TOF settings
    Bool_t           fTOFtuneMC;        //  TOF settings
    Double_t         fTOFresolution;    //  TOF settings
    Double_t         fMinTOF;           //  range for TOF PID (min)
    Double_t         fMaxTOF;           //  range for TOF PID (max)
    Int_t            fLastRun;          //  last run number

    ClassDef(AliRsnCutESD2010, 1)
};

#endif
