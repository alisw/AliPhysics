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

#ifndef ALIRSNCUTAOD2010_H
#define ALIRSNCUTAOD2010_H

#include "AliPID.h"
#include "AliAODTrack.h"
#include "AliAODpidUtil.h"

#include "AliRsnCut.h"

class AliRsnCutAOD2010 : public AliRsnCut
{
  public:

    AliRsnCutAOD2010(const char *name = "cutAOD2010", Bool_t isMC = kFALSE);
    AliRsnCutAOD2010(const AliRsnCutAOD2010& copy);
    virtual ~AliRsnCutAOD2010() {;};

    virtual Bool_t   IsSelected(TObject *object);
    
    void             SetMC       (Bool_t yn = kTRUE);
    void             SetCheckITS (Bool_t yn = kTRUE) {fCheckITS = yn;}
    void             SetCheckTPC (Bool_t yn = kTRUE) {fCheckTPC = yn;}
    void             SetCheckTOF (Bool_t yn = kTRUE) {fCheckTOF = yn;}
    void             SetUseGlobal(Bool_t yn = kTRUE) {fUseGlobal = yn;}
    void             SetUseITSSA (Bool_t yn = kTRUE) {fUseITSSA = yn;}
    void             SetPIDtype  (AliPID::EParticleType pid) {fPIDtype = pid;}
    void             SetMaxEta   (Double_t eta) {fMaxEta = eta;}
    
    void             SetTPCminNclusters(Int_t v)         {fTPCminNclusters = v;}
    void             SetTPCmaxChi2(Double_t v)           {fTPCmaxChi2 = v;}     
    void             SetTPCmaxNSigmaDCA(Double_t v)      {fTPCmaxNSigmaDCA = v;}
    void             SetTPCparamDCA(Int_t i, Double_t v) {if (i >= 0 && i < 3) fTPCparamDCA[i] = v;} 
    void             SetTPClowBand(Double_t v)           {fTPClowBand = v;}     
    void             SetTPChighBand(Double_t v)          {fTPChighBand = v;}    
    void             SetTPClowLimit(Double_t v)          {fTPClowLimit = v;}    
    
    void             SetITSminNclusters(Int_t v)         {fITSminNclusters = v;}
    void             SetITSmaxChi2(Double_t v)           {fITSmaxChi2 = v;}     
    void             SetITSmaxNSigmaDCA(Double_t v)      {fITSmaxNSigmaDCA = v;}
    void             SetITSparamDCA(Int_t i, Double_t v) {if (i >= 0 && i < 3) fITSparamDCA[i] = v;} 
    void             SetITSparamBB(Int_t i, Double_t v)  {if (i >= 0 && i < 3) fITSparamBB[i] = v;}  
    void             SetITSband(Double_t v)              {fITSband = v;}    

    void             SetTOFrange(Double_t v1, Double_t v2) {fTOFlowLimit = v1; fTOFhighLimit = v2;}
    
    AliAODpidUtil*   GetPIDUtil() {return &fPID;}

  protected:
  
    AliRsnCutAOD2010& operator=(const AliRsnCutAOD2010& /*copy*/) {return (*this);}
  
    Bool_t                fIsMC;             // switch for MC analysis
    Bool_t                fCheckITS;         // switch for ITS dE/dx check
    Bool_t                fCheckTPC;         // switch for TPC dE/dx check
    Bool_t                fCheckTOF;         // switch for TOF time check
    Bool_t                fUseGlobal;        // switch to use TPC global tracks
    Bool_t                fUseITSSA;         // switch to use ITS standalone tracks
    
    Double_t              fMaxEta;           // eta range for tracks
    
    AliPID::EParticleType fPIDtype;          // particle type for which PID is checked   
    
    Int_t                 fTPCminNclusters;  // minimum number of clusters in TPC
    Double_t              fTPCmaxChi2;       // maximum chi2 / number of clusters in TPC
    Double_t              fTPCmaxNSigmaDCA;  // cut in DCA (transv) in numbers of sigma (pt-dependent)
    Double_t              fTPCparamDCA[3];   // parameters to compute sigma for DCA
    Double_t              fTPClowBand;       // large band for low momentum PID
    Double_t              fTPChighBand;      // large band for low momentum PID
    Double_t              fTPClowLimit;      // limit of low momentum region
    
    Int_t                 fITSminNclusters;  // minimum number of clusters in TPC
    Double_t              fITSmaxChi2;       // maximum chi2 / number of clusters in TPC
    Double_t              fITSmaxNSigmaDCA;  // cut in DCA (transv) in numbers of sigma (pt-dependent)
    Double_t              fITSparamDCA[3];   // parameters to compute sigma for DCA
    Double_t              fITSparamBB[5];    // parameters for TPC Bethe-Bloch parameterization
    Double_t              fITSband;          // fiducial region for selection of PID
    
    Double_t              fTOFlowLimit;      // low limit in asymmetric TOF PID cut
    Double_t              fTOFhighLimit;     // low limit in asymmetric TOF PID cut
    
    AliAODpidUtil         fPID;              // PID utility

    ClassDef(AliRsnCutAOD2010, 1)
};

#endif
