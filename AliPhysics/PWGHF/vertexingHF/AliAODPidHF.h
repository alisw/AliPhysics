#ifndef ALIAODPIDHF_H
#define ALIAODPIDHF_H

/* Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//***********************************************************
//// \class Class AliAODPidHF
//// \brief class for PID with AliAODRecoDecayHF
//// \author Authors: D. Caffarri caffarri@pd.infn.it, A.Dainese andrea.dainese@pd.infn.it, S. Dash dash@to.infn.it, F. Prino prino@to.infn.it, R. Romita r.romita@gsi.de, Y. Wang yifei@pi0.physi.uni-heidelberg.de,  J. van der Maarel j.vandermaarel@cern.ch
////***********************************************************

#include <TString.h>
#include <TH1F.h>
#include <TObject.h>
#include "AliAODTrack.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliPID.h"

#include "vector"
using std::vector;

class AliAODPidHF : public TObject{
  
public:
  
  enum ECombDetectors {
    kTPC,
    kTOF,
    kTPCTOF,
    kTPCITS
  };

  enum SystemForNsigmaDataCorr {
    kNone=-1,
    kPbPb010,
    kPbPb3050,
    kPbPb6080
  };

  static const int kMaxEtaBins=5;
  static const int kMaxPBins=10;

  AliAODPidHF();
  AliAODPidHF(const AliAODPidHF& pid);
  virtual ~AliAODPidHF();
  
  //Setters
  void SetSigma(Double_t *sigma){
    for(Int_t i=0; i<fnNSigma; i++) fnSigma[i]=sigma[i];
  }
  void SetSigma(Int_t idet,Double_t sigma){fnSigma[idet]=sigma;return;}
  void SetSigmaForTPC(Double_t *sigma){for(Int_t i=0;i<3;i++) fnSigma[i]=sigma[i];return;}
  void SetSigmaForTPCCompat(Double_t sigma){fnSigmaCompat[0]=sigma;return;}
  void SetSigmaForTOFCompat(Double_t sigma){fnSigmaCompat[1]=sigma;return;}
  void SetSigmaForTOF(Double_t sigma){fnSigma[3]=sigma;return;}
  void SetSigmaForITS(Double_t sigma){fnSigma[4]=sigma;return;}
  void SetTofSigma(Double_t sigma){fTOFSigma=sigma;return;}
  
  void SetCutOnTOFmismatchProb(Double_t cut=0.01){fCutTOFmismatch=cut;}
  void DisableCutOnTOFmismatchProb(){fCutTOFmismatch=999.;}
  
  void SetMinNClustersTPCPID(Int_t minc) {fMinNClustersTPCPID=minc;}
  
  void SetCombinednSigmaCutForPiKP(Float_t sigpi, Float_t sigk, Float_t sigp){
    fMaxnSigmaCombined[0]=sigpi;
    fMaxnSigmaCombined[1]=sigk;
    fMaxnSigmaCombined[2]=sigp;
  }
  void SetTPCnSigmaRangeForPions(Float_t smin, Float_t smax){
    fMinnSigmaTPC[0]=smin;
    fMaxnSigmaTPC[0]=smax;
  }
  void SetTOFnSigmaRangeForPions(Float_t smin, Float_t smax){
    fMinnSigmaTOF[0]=smin;
    fMaxnSigmaTOF[0]=smax;
  }
  void SetTPCnSigmaRangeForKaons(Float_t smin, Float_t smax){
    fMinnSigmaTPC[1]=smin;
    fMaxnSigmaTPC[1]=smax;
  }
  void SetTOFnSigmaRangeForKaons(Float_t smin, Float_t smax){
    fMinnSigmaTOF[1]=smin;
    fMaxnSigmaTOF[1]=smax;
  }
  void SetTPCnSigmaRangeForProtons(Float_t smin, Float_t smax){
    fMinnSigmaTPC[2]=smin;
    fMaxnSigmaTPC[2]=smax;
  }
  void SetTOFnSigmaRangeForProtons(Float_t smin, Float_t smax){
    fMinnSigmaTOF[2]=smin;
    fMaxnSigmaTOF[2]=smax;
  }
  
  void SetPriors(Double_t *priors, Int_t npriors);
  void SetPLimit(Double_t *plim, Int_t npLim);

  void SetAsym(Bool_t asym){fAsym=asym;return;}
  void SetUseAsymmnSigmaTOF(Double_t nsmin, Double_t nsmax, Double_t nscompmin, Double_t nscompmax){
    fUseAsymTOF=kTRUE;
    fLownSigmaTOF=nsmin; fUpnSigmaTOF=nsmax;
    fLownSigmaCompatTOF=nscompmin; fUpnSigmaCompatTOF=nscompmax;
  }
  void SetTPC(Bool_t tpc){fTPC=tpc;return;}
  void SetTOF(Bool_t tof){fTOF=tof;return;}
  void SetITS(Bool_t its){fITS=its;return;}
  void SetTRD(Bool_t trd){fTRD=trd;return;}
  void SetMatch(Int_t match){fMatch=match;return;}
  void SetForceTOFforKaons(Bool_t forceTOF){fForceTOFforKaons=forceTOF;return;}
  void SetCompat(Bool_t comp){fCompat=comp;return;}
  void SetMC(Bool_t mc){fMC=mc;return;}
  void SetMClowenpp2011(Bool_t mc){fMCLowEn2011=mc;return;}
  void SetOnePad(Bool_t onepad){fOnePad=onepad;return;}
  void SetppLowEn2011(Bool_t opt){fppLowEn2011=opt;return;}
  void SetPbPb(Bool_t pbpb){fPbPb=pbpb;return;}
  void SetPCompatTOF(Double_t pTOF){fPCompatTOF=pTOF;return;}
  void SetTOFdecide(Bool_t tOFdecide){fTOFdecide=tOFdecide;return;}
  void SetOldPid(Bool_t oldPid){fOldPid=oldPid;return;}
  void SetPtThresholdTPC(Double_t ptThresholdTPC){fPtThresholdTPC=ptThresholdTPC;return;}
  void SetMaxTrackMomForCombinedPID(Double_t mom){fMaxTrackMomForCombinedPID=mom;}
  void SetPidResponse(AliPIDResponse *pidResp) {fPidResponse=pidResp;return;}
  void SetCombDetectors(ECombDetectors pidComb) {
    fCombDetectors=pidComb;
  }
  void SetPionPriorHisto(TH1F* histo){
    if(fPriorsH[AliPID::kPion]) delete fPriorsH[AliPID::kPion];
    fPriorsH[AliPID::kPion] = new TH1F(*histo);
  }
  void SetKaonPriorHisto(TH1F* histo){
    if(fPriorsH[AliPID::kKaon]) delete fPriorsH[AliPID::kKaon];
    fPriorsH[AliPID::kKaon] = new TH1F(*histo);
  }
  void SetProtonPriorHisto(TH1F* histo){
    if(fPriorsH[AliPID::kProton]) delete fPriorsH[AliPID::kProton];
    fPriorsH[AliPID::kProton] = new TH1F(*histo);
  }
  
  
  //Getters
  
  Int_t GetnSigmaTPC(AliAODTrack *track, Int_t species, Double_t &sigma) const;
  Int_t GetnSigmaTOF(AliAODTrack *track, Int_t species, Double_t &sigma) const;
  Int_t GetnSigmaITS(AliAODTrack *track, Int_t species, Double_t &sigma) const;
  Double_t GetSigma(Int_t idet) const{return fnSigma[idet];}
  Double_t GetTofSigma() const{return fTOFSigma;}
  //void GetPriors(Double_t *priors) const{priors=fPriors;return;}
  //void GetPLimit(Double_t *plim) const{plim=fPLimit;}
  void GetPriors(Double_t *priors) const{for(Int_t i=0;i<fnPriors;i++){priors[i]=fPriors[i];}return;}
  void GetPLimit(Double_t *plim) const{for(Int_t i=0;i<fnPLimit;i++){plim[i]=fPLimit[i];}return;}
  Bool_t GetAsym() const{return fAsym;}
  Bool_t GetTPC() const{return fTPC;}
  Bool_t GetTOF() const{return fTOF;}
  Bool_t GetITS() const{return fITS;}
  Bool_t GetTRD() const{return fTRD;}
  Int_t GetMatch() const{return fMatch;}
  Bool_t GetForceTOFforKaons() const{return fForceTOFforKaons;}
  Bool_t GetCompat() const{return fCompat;}
  Bool_t GetMC() const{return fMC;}
  Bool_t GetOnePad() const{return fOnePad;}
  Bool_t GetppLowEn2011() const {return fppLowEn2011;}
  Bool_t GetMCLowEn2011() const {return fMCLowEn2011;}
  Bool_t GetPbPb() const{return fPbPb;}
  Bool_t GetTOFdecide() const{return fTOFdecide;}
  Double_t GetPCompatTOF() const{return fPCompatTOF;}
  Double_t GetnSigmaCompatTPC() const{return fnSigmaCompat[0];}
  Double_t GetnSigmaCompatTOF() const{return fnSigmaCompat[1];}
  Bool_t GetOldPid(){return fOldPid;}
  Double_t GetPtThresholdTPC(){return fPtThresholdTPC;}
  Double_t GetMaxTrackMomForCombinedPID(){return fMaxTrackMomForCombinedPID;}
  AliPIDResponse *GetPidResponse() const {return fPidResponse;}
  AliPIDCombined *GetPidCombined() const {return fPidCombined;}
  ECombDetectors GetCombDetectors() const {
    return fCombDetectors;
  }
  Bool_t GetUseCombined() {return fUseCombined;}
  Bool_t GetDefaultPriors() {return fDefaultPriors;}
  
  Int_t RawSignalPID (AliAODTrack *track, TString detector) const;
  Bool_t IsKaonRaw (AliAODTrack *track, TString detector) const;
  Bool_t IsPionRaw (AliAODTrack *track, TString detector) const;
  Bool_t IsProtonRaw (AliAODTrack *track, TString detector) const;
  Bool_t IsElectronRaw (AliAODTrack *track, TString detector) const;
  void CombinedProbability(AliAODTrack *track,Bool_t *type) const; //0 = pion, 1 = kaon, 2 = proton
  Bool_t CheckStatus(AliAODTrack *track,TString detectors) const;
  
  Bool_t CheckITSPIDStatus(AliAODTrack *track) const;
  Bool_t CheckTPCPIDStatus(AliAODTrack *track) const;
  Bool_t CheckTOFPIDStatus(AliAODTrack *track) const;
  Bool_t CheckTRDPIDStatus(AliAODTrack *track) const;
  
  Bool_t TPCRawAsym(AliAODTrack* track,Int_t specie) const;
  Int_t MatchTPCTOF(AliAODTrack *track,Int_t specie);
  ///PID nSigma strategy closer to the Bayesian approach with Max. prob.
  Int_t MatchTPCTOFMin(AliAODTrack *track, Int_t specie);
  
  Int_t MakeRawPid(AliAODTrack *track,Int_t specie); ///general method to perform PID using raw signals
  
  Bool_t IsTOFPiKexcluded(AliAODTrack *track,Double_t nsigmaK);
  
  Bool_t IsExcluded(AliAODTrack *track, Int_t labelTrack, Double_t nsigmaCut, TString detector);
  
  void GetTPCBetheBlochParams(Double_t alephParameters[5]) const;
  void SetBetheBloch();
  /// method for AliPIDCombined object
  void SetSelectedSpecies(Int_t ispecies = AliPID::kSPECIES){GetPidCombined()->SetSelectedSpecies(ispecies);};
  void SetPriorDistribution(AliPID::EParticleType type,TH1F *prior);
  void DrawPrior(AliPID::EParticleType type);
  void SetPriorsHistos(TString priorFileName);
  void SetUpCombinedPID();
  void SetUseCombined(Bool_t useCombined=kTRUE) {fUseCombined=useCombined;}
  void SetUseDefaultPriors(Bool_t defaultP)	    {fDefaultPriors=defaultP;}
  Int_t ApplyPidTPCRaw(AliAODTrack *track,Int_t specie) const;
  Int_t ApplyPidTOFRaw(AliAODTrack *track,Int_t specie) const;
  Int_t ApplyPidITSRaw(AliAODTrack *track,Int_t specie) const;
  Int_t ApplyTOFCompatibilityBand(AliAODTrack *track,Int_t specie) const;
  
  void PrintAll() const;

  ///Assymetric PID using histograms
  void SetIdBand(AliPID::EParticleType specie, AliPIDResponse::EDetector detector, TH1F *min, TH1F *max);
  void SetIdBand(AliPID::EParticleType specie, AliPIDResponse::EDetector detector, TF1 *min, TF1 *max);
  void SetCompBand(AliPID::EParticleType specie, AliPIDResponse::EDetector detector, TH1F *min, TH1F *max);
  void SetCompBand(AliPID::EParticleType specie, AliPIDResponse::EDetector detector, TF1 *min, TF1 *max);
  Bool_t CheckDetectorPIDStatus(AliPIDResponse::EDetector detector, AliAODTrack *track);
  Float_t NumberOfSigmas(AliPID::EParticleType specie, AliPIDResponse::EDetector detector, AliAODTrack *track);
  Int_t CheckBands(AliPID::EParticleType specie, AliPIDResponse::EDetector detector, AliAODTrack *track);
  TF1 *GetIdBandMin(AliPID::EParticleType specie, AliPIDResponse::EDetector detector) { return fIdBandMin[((int) specie)][((int) detector)]; }
  TF1 *GetIdBandMax(AliPID::EParticleType specie, AliPIDResponse::EDetector detector) { return fIdBandMax[((int) specie)][((int) detector)]; }
  TF1 *GetCompBandMin(AliPID::EParticleType specie, AliPIDResponse::EDetector detector) { return fCompBandMin[((int) specie)][((int) detector)]; }
  TF1 *GetCompBandMax(AliPID::EParticleType specie, AliPIDResponse::EDetector detector) { return fCompBandMax[((int) specie)][((int) detector)]; }

  ///Some suggested asymmetric PID
  void SetShiftedAsymmetricPID();
  void SetIdAsymmetricPID();
  void SetIdCompAsymmetricPID();
  
  ///Set Nsigma data-driven correction
  void EnableNsigmaTPCDataCorr(Int_t run, Int_t system);

  //method to get parameters for NsigmaTPC correction
  static void SetNsigmaTPCDataDrivenCorrection(Int_t run, Int_t system, Int_t &nPbins, Float_t Plims[kMaxPBins+1], Int_t &nEtabins, Float_t absEtalims[kMaxEtaBins+1], 
  vector<vector<Float_t> > &meanNsigmaTPCpion, vector<vector<Float_t> > &meanNsigmaTPCkaon, vector<vector<Float_t> > &meanNsigmaTPCproton, 
  vector<vector<Float_t> > &sigmaNsigmaTPCpion, vector<vector<Float_t> > &sigmaNsigmaTPCkaon, vector<vector<Float_t> > &sigmaNsigmaTPCproton);

protected:
  
  
private:

  AliAODPidHF& operator=(const AliAODPidHF& pid);

  void GetNsigmaTPCMeanSigmaData(Float_t &mean, Float_t &sigma, AliPID::EParticleType species, Float_t pTPC, Float_t eta) const;

  Int_t fnNSigma; /// number of sigmas
  /// sigma for the raw signal PID: 0-2 for TPC, 3 for TOF, 4 for ITS
  Double_t *fnSigma; // [fnNSigma], sigma for the raw signal PID: 0-2 for TPC, 3 for TOF, 4 for ITS
  Double_t fTOFSigma; /// TOF precision
  Double_t fCutTOFmismatch; /// Cut of TOF mismatch probability
  UInt_t fMinNClustersTPCPID;       /// Minimum TPC PID clusters cut
  Int_t fnPriors; /// number of priors
  /// set of priors
  Double_t *fPriors; // [fnPriors], set of priors
  Int_t fnPLimit; /// number of Plimit
  /// limit of p intervals for asimmetric PID: fPLimit<p[0], fPLimit[0]<p<fPLimit[1], p>fPLimit[1]
  Double_t *fPLimit; // [fnPLimit], limit of p intervals for asimmetric PID: fPLimit<p[0], fPLimit[0]<p<fPLimit[1], p>fPLimit[1]
  Bool_t fAsym; /// asimmetric PID required
  Bool_t fTPC; /// switch to include or exclude TPC
  Bool_t fTOF; /// switch to include or exclude TOF
  Bool_t fITS; /// switch to include or exclude ITS
  Bool_t fTRD; /// switch to include or exclude TRD
  Int_t fMatch; /// switch to combine the info from more detectors: 1 = || , 2 = &, 3 = p region
  Bool_t fForceTOFforKaons; /// force TOF for kaons in mode fMatch==5
  Bool_t fCompat; /// compatibility region : useful only if fMatch=1
  Double_t fPCompatTOF; ///  compatibility p limit for TOF
  Bool_t fUseAsymTOF; /// flag for using asymmetrig nSigmaCut in TOF for fMatch==1
  Double_t fLownSigmaTOF;  /// lower nsigma TOF (for fUseAsymTOF)
  Double_t fUpnSigmaTOF;  /// upper nsigma TOF (for fUseAsymTOF)
  Double_t fLownSigmaCompatTOF;  /// lower nsigma TOF (for fUseAsymTOF)
  Double_t fUpnSigmaCompatTOF;  /// upper nsigma TOF (for fUseAsymTOF)
  Int_t fnNSigmaCompat; /// number of sigmas
  /// 0: n sigma for TPC compatibility band, 1: for TOF
  Double_t *fnSigmaCompat; //[fnNSigmaCompat]  0: n sigma for TPC compatibility band, 1: for TOF
  Double_t fMaxnSigmaCombined[3];  /// nSigma cut for pi,K,p (TPC^2+TOF^2)
  Double_t fMinnSigmaTPC[3]; /// min. of nSigma range for pi,K,p in TPC (match==5)
  Double_t fMaxnSigmaTPC[3]; /// max. of nSigma range for pi,K,p in TPC (match==5)
  Double_t fMinnSigmaTOF[3]; /// min. of nSigma range for pi,K,p in TOF (match==5)
  Double_t fMaxnSigmaTOF[3]; /// max. of nSigma range for pi,K,p in TOF (match==5)
  Bool_t fMC; /// MC(kTRUE) or real data (kFALSE, default option)
  Bool_t fOnePad; ///  real data with one pad clusters
  Bool_t fMCLowEn2011; ///  MC for low energy MC
  Bool_t fppLowEn2011; ///  Data for low energy pp 2011
  Bool_t fPbPb; ///  real data PbPb
  Bool_t fTOFdecide; ///  real data PbPb
  Bool_t fOldPid; ///  old PID method implemented
  Double_t fPtThresholdTPC; ///  pT threshold to use TPC PID
  Double_t fMaxTrackMomForCombinedPID; /// momentum threshold to use PID
  AliPIDResponse *fPidResponse; //!<! pid response
  AliPIDCombined* fPidCombined; //!<! combined PID object
  
  AliTPCPIDResponse* fTPCResponse; //!<! TPC response
  
  TH1F* fPriorsH[AliPID::kSPECIES]; /// priors histos
  ECombDetectors fCombDetectors; /// detectors to be involved for combined PID
  Bool_t fUseCombined; /// detectors to be involved for combined PID
  Bool_t fDefaultPriors; /// use default priors for combined PID

  /// Storage of identification/compatibility band for different species and detectors:
  TF1 *fIdBandMin[AliPID::kSPECIES][4];
  TF1 *fIdBandMax[AliPID::kSPECIES][4];
  TF1 *fCompBandMin[AliPID::kSPECIES][4];
  TF1 *fCompBandMax[AliPID::kSPECIES][4];

  Bool_t fApplyNsigmaTPCDataCorr; /// flag to enable data-driven NsigmaTPC correction
  vector<vector<Float_t> > fMeanNsigmaTPCPionData; /// array of NsigmaTPC pion mean in data for different eta bins
  vector<vector<Float_t> > fMeanNsigmaTPCKaonData; /// array of NsigmaTPC kaon mean in data for different eta bins
  vector<vector<Float_t> > fMeanNsigmaTPCProtonData; /// array of NsigmaTPC proton mean in data for different eta bins 
  vector<vector<Float_t> > fSigmaNsigmaTPCPionData; /// array of NsigmaTPC pion mean in data for different eta bins 
  vector<vector<Float_t> > fSigmaNsigmaTPCKaonData; /// array of NsigmaTPC kaon mean in data for different eta bins 
  vector<vector<Float_t> > fSigmaNsigmaTPCProtonData; /// array of NsigmaTPC proton mean in data for different eta bins 
  Float_t fPlimitsNsigmaTPCDataCorr[kMaxPBins+1]; /// array of p limits for data-driven NsigmaTPC correction
  Int_t fNPbinsNsigmaTPCDataCorr;/// number of p bins for data-driven NsigmaTPC correction
  Float_t fEtalimitsNsigmaTPCDataCorr[kMaxEtaBins+1]; /// array of eta limits for data-driven NsigmaTPC correction
  Int_t fNEtabinsNsigmaTPCDataCorr;/// number of eta bins for data-driven NsigmaTPC correction

  /// \cond CLASSIMP
  ClassDef(AliAODPidHF,26); /// AliAODPid for heavy flavor PID
  /// \endcond

};

struct HistFunc {
   HistFunc(TH1F *f): fHist(f) {}
   double operator() (double *x, double * ) const {
      TAxis *axis = fHist->GetXaxis();
      Int_t bin = axis->FindBin(x[0]);
      if (x[0] == axis->GetXmax()) {
        bin = axis->GetNbins();
      }
      return fHist->GetBinContent(bin);

   }
   TH1F *fHist;
};

#endif

