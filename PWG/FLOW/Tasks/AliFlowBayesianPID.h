/*
-----------------------------------------------------------------------------------
 AliFlowBayesianPID class implemented by F. Noferini (noferini@bo.infn.it)
 needed to use Bayesian probability in the flow package (used in AliFlowTrackCuts)
-----------------------------------------------------------------------------------
*/
#ifndef ALIFLOWBAYESIANPID_H
#define ALIFLOWBAYESIANPID_H

#include "AliESDpid.h"
#include "AliPIDResponse.h"

class TObject;
class TDatabasePDG;
class AliESDEvent;
class AliESDtrack;
class AliAODEvent;
class AliAODTrack;
class TH2D;
class TSpline3;
class TF1;
class TH1D;

/*
HOW TO

1) in your main program:

  AliFlowBayesianPID *mypid = new AliFlowBayesianPID();

or (if you want a global AliESDpid object)

  AliESDpid *PIDesd = new AliESDpid();
  ....
  AliFlowBayesianPID *mypid = new AliFlowBayesianPID(PIDesd);

for data starting from the PbPb pass2 reconstruction

  mypid->SetNewTrackParam();


before to loop on the on ESD tracks 
 mypid->SetDetResponse(esdEvent, centrality,AliESDpid::kTOF_T0,kFALSE); // centrality > 0 = PbPb or -1 for pp collisions

before to loop on the on AOD tracks 
 mypid->SetDetResponse(aodEvent, centrality); // centrality > 0 = PbPb or -1 for pp collisions

if you want to use a dE/dx depdence on DeltaPhi set the EP with its resolution otherwise it will be skipped
 mypid->SetPsiCorrectionDeDx(psiEP,resEP);


     for(...){ // track loop
       mypid->ComputeProb(track,centrality); // both for ESD and AOD tracks
       Float_t *prob = mypid->GetProb(); // Bayesian Probability (from 0 to 4) (Combined TPC || TOF) including a tuning of priors and TOF mismatch parameterization
     }


More details:
     // for the single detector weights (no priors)
     Float_t *tpcWeight = mypid->GetWeights(0); // TPC weights (equal weights in case of no kTPCpid)
     Float_t *tofWeight = mypid->GetWeights(1); // TOF weights (equal weights in case of no kTOFpid)

     // more tools
     TSpline3 *mismatchShape = mypid->GetMismatch(); // time distribution of mismatched tracks - L/c (L is the distance of the pad from the PV)
     TF1 *tpcprob = mypid->GetTPCprob(); // normalized TPC response function (in Nsigmas)  
     TF1 *tofprob = mypid->GetTOFprob(); // normalized TOF response function (in Nsigmas)  

     TH2D *hPr = mypid->GetHistoPriors(isp); // 2D (centrality - pT) histo for the priors of specie-isp (centrality < 0 means pp collisions)
                                             // all the priors are normalized to the pion ones

*/

class AliFlowBayesianPID : public AliPIDResponse{
 public:
  AliFlowBayesianPID(AliESDpid *esdpid=NULL); 
  virtual ~AliFlowBayesianPID();
  
  // virtual method of AliPIDResponse


  // setter
  void SetDetResponse(AliESDEvent *esd,Float_t centrality=-1.0,EStartTimeType_t flagStart=AliESDpid::kTOF_T0,Bool_t /*recomputeT0TOF*/=kFALSE);
  void SetDetResponse(AliAODEvent *aod,Float_t centrality=-1.0,EStartTimeType_t flagStart=AliESDpid::kTOF_T0);
  void SetNewTrackParam(Bool_t flag=kTRUE){fNewTrackParam=flag;};
  void SetDetAND(Int_t idet){if(idet < fgkNdetectors && idet >= 0) fMaskAND[idet] = kTRUE;};
  void SetDetOR(Int_t idet){if(idet < fgkNdetectors && idet >= 0) fMaskOR[idet] = kTRUE;};
  void ResetDetAND(Int_t idet){if(idet < fgkNdetectors && idet >= 0) fMaskAND[idet] = kFALSE;};
  void ResetDetOR(Int_t idet){if(idet < fgkNdetectors && idet >= 0) fMaskOR[idet] = kFALSE;};
  void SetPsiCorrectionDeDx(Float_t psi,Float_t res);
  void SetMC(Bool_t flag){fIsMC=flag;};

  // getter
  AliESDpid* GetESDpid(){return fPIDesd;};
  const TH2D *GetHistoPriors(Int_t specie) const {if(specie >=0 && specie < fgkNspecies) return fghPriors[specie]; else return NULL;};
  TSpline3 *GetMismatch();  
  const TF1 *GetTOFprob() const {return fTOFResponseF;};
  const TF1 *GetTPCprob() const {return fTPCResponseF;};
  const Float_t *GetWeights(Int_t det) const {if(det < fgkNdetectors && det >= 0){return fWeights[det];} else{return NULL;}};
  Float_t *GetProb() {return fProb;};
  Float_t GetTOFMismWeight() const {return fWTofMism;};
  Float_t GetTOFMismProb() const {return fProbTofMism;};
  Float_t GetMassOverZ() const {return fMassTOF;};
  Float_t GetZ() const {return fZ;};
  Bool_t GetDetANDstatus(Int_t idet) const {if(idet < fgkNdetectors && idet >= 0){return fMaskAND[idet];} else{return kFALSE;} };
  Bool_t GetDetORstatus(Int_t idet) const {if(idet < fgkNdetectors && idet >= 0){return fMaskOR[idet];} else{return kFALSE;} };
  Bool_t GetCurrentMask(Int_t idet) const {if(idet < fgkNdetectors && idet >= 0){return fMaskCurrent[idet];} else{return kFALSE;} };

  Float_t GetExpDeDx(const AliVTrack *t,Int_t iS) const;
  Float_t GetExpDeDx(const AliVTrack *t,Float_t m) const;

  // methods for Bayesina Combined PID
  void ComputeWeights(const AliESDtrack *t);
  void ComputeProb(const AliESDtrack *t,Float_t); // obsolete method
  void ComputeProb(const AliESDtrack *t){ComputeProb(t,0.0);}; 
  void ComputeWeights(const AliAODTrack *t,const AliAODEvent *aod=NULL);
  void ComputeProb(const AliAODTrack *t,const AliAODEvent *aod=NULL); // obsolete method

  void SetTOFres(Float_t res){fTOFresolution=res;};

  Float_t GetDeDx() const {return fDedx;};

 private: 
  void SetPriors();

  static const Int_t fgkNdetectors = 2; // Number of detector used for PID
  static const Int_t fgkNspecies = 9;// 0=el, 1=mu, 2=pi, 3=ka, 4=pr, 5=deuteron, 6=triton, 7=He3 
  static TH2D* fghPriors[fgkNspecies]; // histo with priors (hardcoded)
  static TSpline3 *fgMism; // function for mismatch

  AliESDpid *fPIDesd;//ESDpid object
  TDatabasePDG *fDB; // Database pdg
  Double_t fMass[fgkNspecies]; // mass for el(0),mu(1),pi(2),K(3),p(4)

  Bool_t fNewTrackParam; // switch for new tracking resolution TOF parameterization
  Float_t fTOFresolution; // TOF res needed only if T0-TOF should be recomputed

  AliFlowBayesianPID(const AliFlowBayesianPID&); // not implemented
  AliFlowBayesianPID& operator=(const AliFlowBayesianPID&); // not implemented 

  TF1 *fTOFResponseF; // TOF Gaussian+tail response function (tail at 1.1 sigma)
  TF1 *fTPCResponseF; // TPC Gaussian+tail response function (tail at 1.8 sigma)

  Float_t fWeights[fgkNdetectors][fgkNspecies]; // weights: 0=tpc,1=tof
  Float_t fProb[fgkNspecies],fWTofMism,fProbTofMism; // Bayesian Combined PID + mismatch weights and probability 

  Float_t fZ,fMassTOF; //measured charge(Z) and mass/Z
  TF1 *fBBdata; // Bethe Bloch function (needed to compute the charge of the particle)

  Bool_t fMaskAND[fgkNdetectors],fMaskOR[fgkNdetectors],fMaskCurrent[fgkNdetectors]; // mask detector should be used

  Float_t fCurrCentrality; // Centrality in current event

  Float_t fPsi,fPsiRes;    // RP and its resolution for the event (999 if not available) to correct dEdx for p > 1

  Bool_t fIsMC; // switch for MC analysis

  Float_t fDedx; // dE/dx tuned for MC

  Bool_t fIsTOFheaderAOD; // check the TOF header in AOD

  static TH1D *fgHtofChannelDist; // channel distance from IP

  ClassDef(AliFlowBayesianPID, 9); // example of analysis
};

#endif



