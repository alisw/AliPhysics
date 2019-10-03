#ifndef ALIMESTRACKINFO_H
#define ALIMESTRACKINFO_H

////////////////////////////////////////////////////////////////////////////
//  Track summary data for the Multiplicity and Event Shape group         //
//  Authors:                                                              //
//    Cristi Andrei <Cristian.Andrei@cern.ch>                             //
//    Andrei Herghelegiu <aherghe@niham.nipne.ro>                         //
//    Madalina Tarzila <mtarzila@niham.nipne.ro>                          //
//                                                                        //
////////////////////////////////////////////////////////////////////////////
#ifndef ALIVParticle_H
#include <AliVParticle.h>
#endif
#ifndef ALIPID_H
#include <AliPID.h>
#endif

class AliESDtrack;
class AliPIDResponse;
class AliPIDCombined;
class AliMCParticle;
class AliStack;
class AliMEStrackInfo : public AliVParticle
{
public:
  enum EMESdetectors{
     kITS = 0
    ,kTPC
    ,kTOF
    ,kNdet
  };

  class AliMESpid : public TObject
  {
  friend class AliMEStrackInfo;
  public:
    AliMESpid();
    const Int_t*    GetNsigma(EMESdetectors det) const { return det>=0&&det<kNdet?fNsigma[det]:NULL; }
    const Double_t* GetProb(EMESdetectors det) const   { return det>=0&&det<kNdet?fProb[det]:NULL; }
    Double_t        GetRaw(EMESdetectors det) const    { return det>=0&&det<kNdet?fRaw[det]:-999.; }
//     Float_t         GetTOFmisProb() const              { return fTOFmisProb; }
    Bool_t         GetTOFmisProb() const              { return fTOFmisProb; }
  protected:
    Int_t    fNsigma[kNdet][AliPID::kSPECIES]; // n sigma for each detector
//     Float_t  fTOFmisProb;                      // TOFMismatchProbability
    Bool_t  fTOFmisProb;                      // TOFMismatchProbability
    Double_t fProb[kNdet][AliPID::kSPECIES];   // probabilty per detector & specie
    Double_t fRaw[kNdet];                      // raw probabilty per detector dEdx(ITS & TPC), beta (TOF)
    //static AliPriors *fPrior;                // priors object
    ClassDef(AliMESpid, 1)                     // PID info for MES
  };

  class AliMESfilterParam : public TObject
  {
  friend class AliMEStrackInfo;
  public:
    AliMESfilterParam();
    Int_t     GetNcl() const        { return fNcl; }
    Double_t  GetChi2Cl() const     { return fChi2Cl; }
  protected:
    Int_t     fNcl;                           // no of clusters in TPC
    Double_t  fChi2Cl;                        // chi2 per cluster in TPC
    ClassDef(AliMESfilterParam, 1)            // Filtering parameters for track
  };
  enum EMESorigin{
    kUnknown=0
    ,kPrimary
    ,kMaterial
    ,kSecondary
  };
  enum EMESfilterId{
     kNone = 0       // not filtered out
    ,kPhys           // filtered by physics selection
    ,kNIHAM          // filtered by NIHAM selection
  };
  enum EMESdetStat{
     kNo   = 0       // not filtered out
    ,kTRD            // TRD component
    ,kAny            // any other detector flag
  };
  AliMEStrackInfo();
  AliMEStrackInfo(const AliMEStrackInfo &t);
  AliMEStrackInfo(AliESDtrack *t, AliPIDResponse *rpid, AliPIDCombined *pidComb);
  AliMEStrackInfo(AliMCParticle *t, AliStack *mcStack);
  ~AliMEStrackInfo();

  Bool_t      PxPyPz(Double_t*) const;
  Double_t    Theta() const            { return 0.5*TMath::Pi()-TMath::ATan(Pz()/(Pt()+1.e-13)); }
  Double_t    E() const                { return TMath::Sqrt(M()*M()+P()*P()); }
  Double_t    M() const;
  Int_t       GetLabel() const         { return fTrackId; }
  Int_t       PdgCode() const;
  const       Double_t* PID() const;

  Short_t     Charge() const           { return fPt>0?1:-1; }
  Double_t    Eta() const              { return fEta; }
  AliMESfilterParam*  GetFilterParam() const { return fFilterParam; }
  Double_t    Phi() const              { return fPhi; }
  Double_t    P() const                { return fP; }
  Double_t    Px() const               { return Pt()*TMath::Cos(fPhi); }
  Double_t    Py() const               { return Pt()*TMath::Sin(fPhi); }
//   Double_t    Pz() const               { return TMath::Sqrt(P()*P() - Pt()*Pt());}  // alex
  Double_t    Pz() const               { return fPz;}
  Double_t    Pt() const               { return TMath::Abs(fPt); }
  const AliMESpid*  GetPID() const     { return &fPID; }
  Short_t     GetTOFbc() const         { return fTOFbc; }
  Bool_t      XvYvZv(Double_t* x) const { x[0] = Xv(); x[1] = Yv(); x[2] = Zv(); return kTRUE; }
  inline Double_t    Rv() const;
  Double_t    Xv() const               { return fPosition[0]; }
  Double_t    Yv() const               { return fPosition[1]; }
  Double_t    Zv() const               { return fPosition[2]; }
  Bool_t    GetDCA(Double_t *dca) const		 { dca[0] = fDCA[0]; dca[1] = fDCA[1]; return kTRUE; }
  Double_t    Y() const                { return fY; }
  Double_t    GetdEdx() const         {return fdEdx;}
  Double_t    Getbeta() const          {return fBeta;}
  Double_t    OneOverPt() const        { return Pt()>kAlmost0?1./Pt():0.; }
  Bool_t      IsFilteredOut() const    { return fFilterId>0; }
  Bool_t      IsFilteredOutBy(EMESfilterId id) const     { return TESTBIT(fFilterId, id); }
  Bool_t      HasDetStat(EMESdetStat stat) const         { return TESTBIT(fDetStat, stat); }
  Bool_t      HasOrigin(EMESorigin o) const              { return TESTBIT(fOrigin, o); }

  void        SetPosition(Double_t x, Double_t y, Double_t z) { fPosition[0]=x; fPosition[1]=y; fPosition[2]=z;}
  void        SetPosition(Double_t pos[3])     { memcpy(fPosition, pos, 3*sizeof(Double_t)); }
  void        SetDCA(Double_t xy, Double_t z) { fDCA[0]=xy; fDCA[1]=z;}
  void        SetDCA(Double_t dca[2])     { memcpy(fDCA, dca, 2*sizeof(Double_t)); }
  void        SetDetStat(EMESdetStat s)   { SETBIT(fDetStat, s); }
  void        SetFilteredOutBy(EMESfilterId id)          { SETBIT(fFilterId, id); }
  void        SetLabel(Int_t l)           { fTrackId = l; }
  void        SetEta(Double_t eta)        { fEta = eta; }
  void        SetOrigin(EMESorigin o)     { SETBIT(fOrigin, o); }
  void        SetPhi(Double_t phi)        { fPhi = phi; }
  void        SetP(Double_t p)            { fP = p; }
  void        SetdEdx(Double_t dEdx)            { fdEdx = dEdx; }
  void        Setbeta(Double_t beta)            { fBeta = beta; }
  void        SetPt(Double_t pt, Bool_t chg) { fPt = pt*(chg?1:-1); }
  void        SetTOFmisProb(Double_t tm)  { fPID.fTOFmisProb = tm; }
//   void        SetTOFmisProb(Bool_t tm)  { fPID.fTOFmisProb = tm; }
private:
  AliMEStrackInfo &operator=(const AliMEStrackInfo &track);

  Int_t              fTrackId;          // MC label || id of MC (ESD tracks) || id of ESD (MC tracks)
  UChar_t            fOrigin;           // origin of track : primary, material, secondary etc.
  Short_t            fTOFbc;            // TOF bunch crossing index
  UChar_t            fFilterId;         // filter idx which rejected the track see EMESfilterId
  UChar_t            fDetStat;          // detector contribution to the track see EMESdetStat
  Double_t           fPt;               // q*pt of the track @ vertex
  Double_t           fP;                // momentum of the track @ vertex
  Double_t           fPz;                // momentum along the z axis - use to compute y
  Double_t           fEta;              // eta of the track @ vertex
  Double_t           fPhi;              // phi of the track @ vertex
  Double_t           fY;                // rapidity
  Double_t           fPosition[3];           // Position x,y and z @ vertex
  Float_t           fDCA[2];           // DCAxy and DCAz
  Double_t        fdEdx;              // energy loss in TPC
  Double_t        fBeta;               // computed beta using TOF
  AliMESpid          fPID;              // PID info
  AliMESfilterParam *fFilterParam;      // Filtering parameters used for debugging purposes
  ClassDef(AliMEStrackInfo, 3)          // Track summary data for MultiEvShape
};

Double_t AliMEStrackInfo::Rv() const
{
  Double_t rv2(fPosition[0]*fPosition[0]+fPosition[1]*fPosition[1]);
  return rv2>kAlmost0?TMath::Sqrt(rv2):0.;
}


#endif
