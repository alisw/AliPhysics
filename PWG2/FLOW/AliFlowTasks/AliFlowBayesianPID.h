#ifndef AliFlowBayesianPID_h
#define AliFlowBayesianPID_h

#include "TObject.h"
#include "AliESDpid.h"
#include "AliPIDResponse.h"

class TDatabasePDG;
class AliESDEvent;
class AliESDtrack;
class TH2D;
class TSpline3;
class TF1;
class AliTOFGeometry;
class AliTOFT0maker;

/*
HOW TO

1) in your main program:

  AliFlowBayesianPID *mypid = new AliFlowBayesianPID();

or (if you want a global AliESDpid object)

  AliESDpid *PIDesd = new AliESDpid();
  ....
  gROOT->LoadMacro("AliFlowBayesianPID.cxx++g");
  AliFlowBayesianPID *mypid = new AliFlowBayesianPID(PIDesd);

for data starting from the PbPb pass2 reconstruction

  mypid->SetNewTrackParam();


2) pass the pointer to your task

3) In your Task Exec:
     mypid->SetDetResponse(esdEvent, centrality,AliESDpid::kTOF_T0,kFALSE); // centrality = PbPb centrality class (0-100%) or -1 for pp collisions
     for(...){ // track loop
     mypid->ComputeProb(track,centrality);
     Float_t *prob = mypid->GetProb(); // Bayesian Probability (from 0 to 4) (Combined TPC || TOF) including a tuning of priors and TOF mismatch parameterization

     // for the single detector weights (no priors)
     Float_t *tpcWeight = mypid->GetWeights(0); // TPC weights (equal weights in case of no kTPCpid)
     Float_t *tofWeight = mypid->GetWeights(1); // TOF weights (equal weights in case of no kTOFpid)

     // more tools
     TSpline3 *mismatchShape = mypid->GetMismatch(); // time distribution of mismatched tracks - L/c (L is the distance of the pad from the PV)
     TF1 *tpcprob = mypid->GetTPCprob(); // normalized TPC response function (in Nsigmas)  
     TF1 *tofprob = mypid->GetTOFprob(); // normalized TOF response function (in Nsigmas)  

     TH2D *hPr = mypid->GetHistoPriors(isp); // 2D (centrality - pT) histo for the priors of specie-isp (centrality < 0 means pp collisions)
                                             // all the priors are normalized to the pion ones

  }

*/

class AliFlowBayesianPID : public AliPIDResponse {
 public:
  AliFlowBayesianPID(AliESDpid *esdpid=NULL); 
  virtual ~AliFlowBayesianPID();
  
  // virtual method of AliPIDResponse
  virtual Float_t NumberOfSigmasTOF(const AliVParticle *vtrack, AliPID::EParticleType type) const {if(vtrack) printf("Don't call AliFlowBayesianPID::NumberOfSigmasTOF method (%i)\n",type); return 0.0;} // do not use it

  // setter
  void SetMC(Bool_t flag = kTRUE){fIsMC=flag;} // actually do nothing
  void SetDetResponse(AliESDEvent *esd,Float_t centrality=-1.0,EStartTimeType_t flagStart=AliESDpid::kTOF_T0,Bool_t recomputeT0TOF=kFALSE);
  void SetNewTrackParam(Bool_t flag=kTRUE){fNewTrackParam=flag;};
  void SetDetAND(Int_t idet){if(idet < fNdetectors && idet >= 0) fMaskAND[idet] = kTRUE;};
  void SetDetOR(Int_t idet){if(idet < fNdetectors && idet >= 0) fMaskOR[idet] = kTRUE;};
  void ResetDetAND(Int_t idet){if(idet < fNdetectors && idet >= 0) fMaskAND[idet] = kFALSE;};
  void ResetDetOR(Int_t idet){if(idet < fNdetectors && idet >= 0) fMaskOR[idet] = kFALSE;};

  // getter
  AliESDpid* GetESDpid(){return fPIDesd;};
  TH2D *GetHistoPriors(Int_t specie){if(specie >=0 && specie < fNspecies) return hPriors[specie]; else return NULL;};
  TSpline3 *GetMismatch();  
  TF1 *GetTOFprob(){return fTOFResponse;};
  TF1 *GetTPCprob(){return fTPCResponse;};
  Float_t *GetWeights(Int_t det){if(det < fNdetectors && det >= 0){return fWeights[det];} else{return NULL;}};
  Float_t *GetProb(){return fProb;};
  Float_t GetTOFMismWeight(){return fWTofMism;};
  Float_t GetTOFMismProb(){return fProbTofMism;};
  Float_t GetMassOverZ(){return fMassTOF;};
  Float_t GetZ(){return fZ;};
  Bool_t GetDetANDstatus(Int_t idet){if(idet < fNdetectors && idet >= 0){return fMaskAND[idet];} else{return kFALSE;} };
  Bool_t GetDetORstatus(Int_t idet){if(idet < fNdetectors && idet >= 0){return fMaskOR[idet];} else{return kFALSE;} };
  Bool_t GetCurrentMask(Int_t idet){if(idet < fNdetectors && idet >= 0){return fMaskCurrent[idet];} else{return kFALSE;} };

  // methods for Bayesina Combined PID
  void ComputeWeights(const AliESDtrack *t,Float_t centr=-1.0);
  void ComputeProb(const AliESDtrack *t,Float_t centr=-1.0);

  void SetTOFres(Float_t res){fTOFres=res;};

 private: 
  void SetPriors();

  static const Int_t fNdetectors = 2;
  static const Int_t fNspecies = 8;// 0=el, 1=mu, 2=pi, 3=ka, 4=pr, 5=deuteron, 6=triton, 7=He3 
  static TH2D* hPriors[fNspecies]; // histo with priors (hardcoded)
  static TSpline3 *fMism; // function for mismatch
  static AliTOFGeometry *fTofGeo; // TOF geometry needed to reproduce mismatch shape

  AliESDpid *fPIDesd;//ESDpid object
  TDatabasePDG *fDB; // Database pdg
  Double_t fMass[fNspecies]; // mass for el(0),mu(1),pi(2),K(3),p(4)

  Bool_t fNewTrackParam; // switch for new tracking resolution TOF parameterization
  Bool_t fIsMC; // switch if MC data
  Float_t fTOFres; // TOF res needed only if T0-TOF should be recomputed

  AliFlowBayesianPID(const AliFlowBayesianPID&); // not implemented
  AliFlowBayesianPID& operator=(const AliFlowBayesianPID&); // not implemented 

  TF1 *fTOFResponse; // TOF Gaussian+tail response function (tail at 1.1 sigma)
  TF1 *fTPCResponse; // TPC Gaussian+tail response function (tail at 1.8 sigma)

  AliTOFT0maker *fTOFmaker; //TOF-T0 maker object

  Float_t fWeights[fNdetectors][fNspecies]; // weights: 0=tpc,1=tof
  Float_t fProb[fNspecies],fWTofMism,fProbTofMism; // Bayesian Combined PID + mismatch weights and probability 

  Float_t fZ,fMassTOF; //measured charge(Z) and mass/Z
  TF1 *fBBdata; // Bethe Bloch function (needed to compute the charge of the particle)

  Bool_t fMaskAND[fNdetectors],fMaskOR[fNdetectors],fMaskCurrent[fNdetectors]; // mask detector should be used

  ClassDef(AliFlowBayesianPID, 4); // example of analysis
};

#endif


