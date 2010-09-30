// 
// Class AliRsnMonitorTrack
//
// Monitor object used
// for storing info ina TTree
// and studying cut values and variables.
//
// author: A. Pulvirenti
//

#ifndef ALIRSNMONITORTRACK_H
#define ALIRSNMONITORTRACK_H

#include "AliPID.h"

class AliESDtrack;
class AliTOFT0maker;
class AliStack;

class AliRsnMonitorTrack : public TObject
{
  public:

    AliRsnMonitorTrack();
    AliRsnMonitorTrack(const AliRsnMonitorTrack& copy);
    //AliRsnMonitorTrack& operator=(const AliRsnMonitorTrack& copy) {MakeCopy(copy); return (*this);}
    virtual ~AliRsnMonitorTrack() { /*nothing*/ }
    
    void       Reset();
    Bool_t     AdoptMC(Int_t label, AliStack *stack);

    Bool_t&    CutsPassed()     {return fCutsPassed;}
    
    Double_t&  PsimX()          {return fPsim[0];}
    Double_t&  PrecX()          {return fPrec[0];}
    Double_t&  PtpcX()          {return fPtpc[0];}
    
    Double_t&  PsimY()          {return fPsim[1];}
    Double_t&  PrecY()          {return fPrec[1];}
    Double_t&  PtpcY()          {return fPtpc[1];}
    
    Double_t&  PsimZ()          {return fPsim[2];}
    Double_t&  PrecZ()          {return fPrec[2];}
    Double_t&  PtpcZ()          {return fPtpc[2];}
    
    Double_t   PsimT()          {return TMath::Sqrt(fPsim[0]*fPsim[0] + fPsim[1]*fPsim[1]);}
    Double_t   PrecT()          {return TMath::Sqrt(fPrec[0]*fPrec[0] + fPrec[1]*fPrec[1]);}
    Double_t   PtpcT()          {return TMath::Sqrt(fPtpc[0]*fPtpc[0] + fPtpc[1]*fPtpc[1]);}
    
    Double_t   Psim()           {return TMath::Sqrt(PsimT()*PsimT() + fPsim[2]*fPsim[2]);}
    Double_t   Prec()           {return TMath::Sqrt(PrecT()*PrecT() + fPrec[2]*fPrec[2]);}
    Double_t   Ptpc()           {return TMath::Sqrt(PtpcT()*PtpcT() + fPtpc[2]*fPtpc[2]);}
    
    Bool_t&    Prim()           {return fPrim;}
    Int_t&     PDG()            {return fPDG;}
    Int_t&     PDGM()           {return fPDGM;}
    Int_t&     Mother()         {return fMother;}

    Bool_t     OkFlag(UInt_t f) {return ((fStatus & f) != 0);}
    UInt_t&    Status()         {return fStatus;}
    Double_t&  Length()         {return fLength;}
    Int_t&     Charge()         {return fCharge;}
    Double_t&  DCAr()           {return fDCA[0];}
    Double_t&  DCAz()           {return fDCA[1];}
    
    Bool_t&    ITSsa()          {return fITSsa;}
    Bool_t&    TOFok()          {return fTOFok;}
    
    Bool_t&    ITSmap(Int_t i)    {if (i>=0 && i<6) return fITSmap[i]; else return fITSmap[0];}
    Int_t      ITScount()         {return (SPDcount() + SDDcount() + SSDcount());}
    Int_t      SPDcount()         {Int_t count=0; if (fITSmap[0]) count++; if (fITSmap[1]) count++; return count;}
    Int_t      SDDcount()         {Int_t count=0; if (fITSmap[2]) count++; if (fITSmap[3]) count++; return count;}
    Int_t      SSDcount()         {Int_t count=0; if (fITSmap[4]) count++; if (fITSmap[5]) count++; return count;}
    Double_t&  ITSchi2()          {return fITSchi2;}
    Double_t&  ITSsignal()        {return fITSsignal;}
    Double_t&  ITSnsigma(Int_t i) {if (i>=0 && i<AliPID::kSPECIES) return fITSnsigma[i]; else return fITSnsigma[0];}
    
    Int_t&     TPCcount()         {return fTPCcount;}
    Double_t&  TPCchi2()          {return fTPCchi2;}
    Double_t&  TPCsignal()        {return fTPCsignal;}
    Double_t&  TPCnsigma(Int_t i) {if (i>=0 && i<AliPID::kSPECIES) return fTPCnsigma[i]; else return fTPCnsigma[0];}
    
    Double_t&  TOFsignal()        {return fTOFsignal;}
    Double_t&  TOFsigma(Int_t i)  {if (i>=0 && i<AliPID::kSPECIES) return fTOFsigma[i]; else return fTOFsigma[0];}
    Double_t&  TOFref(Int_t i)    {if (i>=0 && i<AliPID::kSPECIES) return fTOFref[i]; else return fTOFref[0];}
    Double_t   TOFnsigma(Int_t i) {if (i>=0 && i<AliPID::kSPECIES) return ((fTOFsignal - fTOFref[i]) / fTOFsigma[i]); return 1E6;}

  private:
  
    Bool_t       fCutsPassed;                     // did it pass all defined cuts?
    
    Bool_t       fPrim;                           // is physical primary?
    Int_t        fPDG;                            // true PDG code
    Int_t        fPDGM;                           // PDG code of mother (if any)
    Int_t        fMother;                         // label of mother (if any)

    UInt_t       fStatus;                         // 'status' flag of track in ESD (0 = none)
    Double_t     fLength;                         // integrated length
    Int_t        fCharge;                         // track charge
    Double_t     fDCA[2];                         // DCA ([0] = xy, [1] = z)
        
    Bool_t       fITSsa;                          // to know if its is ITS standalone (otherwise it is TPC)
    Bool_t       fTOFok;                          // to know if track has a TOF match
    
    Bool_t       fITSmap[6];                      // ITS cluster map
    Double_t     fITSchi2;                        // chi2 in ITS
    Double_t     fITSsignal;                      // ITS signal used for PID
    Double_t     fITSnsigma[AliPID::kSPECIES];    // number of sigmas ITS
    
    Int_t        fTPCcount;                       // # TPC clusters
    Double_t     fTPCchi2;                        // TPC chi 2
    Double_t     fTPCsignal;                      // TPC dEdx signal
    Double_t     fTPCnsigma[AliPID::kSPECIES];    // number of sigmas TPC for: e, mu, pi, K, p
    
    Double_t     fTOFsignal;                      // TOF signal
    Double_t     fTOFsigma[AliPID::kSPECIES];     // TOF sigma for: e, mu, pi, K, p
    Double_t     fTOFref[AliPID::kSPECIES];       // expected times for: e, mu, pi, K, p
    
    Double_t     fPsim[3];                        // simulated momentum
    Double_t     fPrec[3];                        // reconstructed momentum
    Double_t     fPtpc[3];                        // reconstructed momentum at the TPC inner wall

    ClassDef(AliRsnMonitorTrack, 1)
};

#endif
