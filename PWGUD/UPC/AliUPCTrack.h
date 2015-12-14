#ifndef ALIUPCTRACK_H
#define ALIUPCTRACK_H

//_____________________________________________________________________________
//    Class for UPC data
//    Author: Jaroslav Adam
//_____________________________________________________________________________

#include "TObject.h"

class TLorentzVector;

class AliUPCTrack : public TObject
{
public:
  AliUPCTrack();

  virtual ~AliUPCTrack();

  void Clear(Option_t * /*option*/ ="");

  //Setters
  void SetPxPyPz(Double_t p[3]) {fP[0]=p[0]; fP[1]=p[1]; fP[2]=p[2];}
  void SetCharge(Short_t charge) {fCharge=charge;}
  void SetMaskMan(UChar_t mask) {fMaskMan=mask;}
  void SetFilterMap(UInt_t map) {fFilterMap=map;}
  void SetFilterBit(UInt_t ibit) {fFilterMap |= (1 << ibit);}
  void SetChi2perNDF(Double_t chi2) {fChi2perNDF=chi2;}

  void SetTPCmomentum(Double_t momentum) {fTPCmomentum=momentum;}
  void SetTPCsignal(Double_t signal) {fTPCsignal=signal;}
  void SetTPCNcls(UShort_t ncls) {fTPCncls=ncls;}
  void SetTPCCrossedRows(Float_t nrows) {fTPCrows=nrows;}
  void SetTPCNclsF(UShort_t ncls) {fTPCnclsF=ncls;}
  void SetTPCNclsS(UShort_t ncls) {fTPCnclsS=ncls;}

  void SetITSchi2perNDF(Double_t chi2) {fITSchi2perNDF=chi2;}
  void SetITSClusterMap(UChar_t cmap) {fITSClusterMap=cmap;}

  void SetTOFsignal(Double_t signal) {fTOFsignal=signal;}

  void SetNSigmasTPC(Int_t i, Float_t sig) {fNSigmasTPC[i]=sig;}
  void SetNSigmasTOF(Int_t i, Float_t sig) {fNSigmasTOF[i]=sig;}

  void SetImpactParameters(Float_t *p, Float_t *cov);
  void SetImpactParametersSPD(Float_t *p, Float_t *cov);
  void SetImpactParametersIP(Float_t *p, Float_t *cov);

  Int_t MakeArrayInt(Int_t size);
  Int_t MakeArrayD(Int_t size);

  //Getters
  void GetPxPyPz(Double_t p[3]) const {p[0]=fP[0]; p[1]=fP[1]; p[2]=fP[2];}
  void GetMomentum(TLorentzVector *v, Double_t mass) const;
  Short_t GetCharge(void) const { return fCharge; }
  UChar_t GetMaskMan(void) const { return fMaskMan; }
  UInt_t GetFilterMap(void) const { return fFilterMap; }
  Bool_t TestFilterBit(UInt_t filterBit) const {return (Bool_t) ((filterBit & fFilterMap) != 0);}
  Double_t GetChi2perNDF(void) const { return fChi2perNDF; }

  Double_t GetTPCmomentum(void) const { return fTPCmomentum; }
  Double_t GetTPCsignal(void) const { return fTPCsignal; }
  UShort_t GetTPCNcls() const { return fTPCncls; }
  Float_t GetTPCCrossedRows() const { return fTPCrows; }
  UShort_t GetTPCNclsF() const { return fTPCnclsF; }
  UShort_t GetTPCNclsS() const { return fTPCnclsS; }

  Double_t GetITSchi2perNDF(void) const { return fITSchi2perNDF; }
  UChar_t GetITSClusterMap(void) const { return fITSClusterMap; }
  Int_t GetITSNcls(void) const;

  Double_t GetTOFsignal(void) const { return fTOFsignal; }

  Float_t GetNSigmasTPC(Int_t i) const { return fNSigmasTPC[i]; }
  Float_t GetNSigmasTOF(Int_t i) const { return fNSigmasTOF[i]; }

  void GetImpactParameters(Double_t &xy, Double_t &z) const {xy = (Double_t) fDZ[0]; z = (Double_t) fDZ[1];}
  void GetImpactParameters(Double_t *p, Double_t *cov) const;

  void GetImpactParametersSPD(Double_t &xy, Double_t &z) const {xy = (Double_t) fdzSPD[0]; z = (Double_t) fdzSPD[1];}
  void GetImpactParametersSPD(Double_t *p, Double_t *cov) const;

  void GetImpactParametersIP(Double_t &xy, Double_t &z) const {xy = (Double_t) fdzIP[0]; z = (Double_t) fdzIP[1];}
  void GetImpactParametersIP(Double_t *p, Double_t *cov) const;

  TArrayI *GetArrayInt(void) const { return fArrayInt; }
  TArrayD *GetArrayD(void) const { return fArrayD; }

protected:
  AliUPCTrack(const AliUPCTrack &o);
  AliUPCTrack &operator=(const AliUPCTrack &o);

  Double_t fP[3]; // momentum px, py, pz
  Short_t fCharge; // track charge
  UChar_t fMaskMan; // 8-bit filter mask, manual fill
  UInt_t fFilterMap; // // filter information, one bit per set of cuts, 32 bit
  Double_t fChi2perNDF; // chi2/NDF of momentum fit
  Double_t fTPCmomentum; // tpc momentum
  Double_t fTPCsignal; // tpc dEdx signal
  UShort_t fTPCncls; // number of clusters assigned in the TPC
  Float_t fTPCrows; // number of of crossed raws in TPC
  UShort_t fTPCnclsF; // number of findable clusters in the TPC
  UShort_t fTPCnclsS; // number of shared clusters in the TPC
  Double_t fITSchi2perNDF; // chi2 in ITS per cluster
  UChar_t fITSClusterMap; // map of clusters, one bit per a layer
  Double_t fTOFsignal; // TOF PID signal
  Float_t fDZ[2]; // impact parameters in XY and Z to default primary vertex
  Float_t fCov[3]; // Covariance matrix of the impact parameters
  Float_t fdzSPD[2]; // SPD impact parameters in XY and Z
  Float_t fCovSPD[3]; // Covariance matrix of the impact parameters to the SPD vertex
  Float_t fdzIP[2]; // impact parameters in XY and Z to nominal interaction point
  Float_t fCovIP[3]; // Covariance matrix of the impact parameters to nominal interaction point
  Float_t fNSigmasTPC[5]; // TPC PID per species, kElectron = 0,  kMuon = 1,  kPion = 2,  kKaon = 3,  kProton = 4
  Float_t fNSigmasTOF[5]; // TOF PID
  TArrayI *fArrayInt; // extension of the central track for other integer parameters
  TArrayD *fArrayD; // extension of the central track for other double parameters

  ClassDef(AliUPCTrack,1);
};

#endif



















