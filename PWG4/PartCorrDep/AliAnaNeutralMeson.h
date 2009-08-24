#ifndef ALIANANEUTRALMESON_H
#define ALIANANEUTRALMESON_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
//a general class to fill two-photon and pi0+gamma invariant mass hisograms 
// to be used to extract pi0, eta and omega raw yield.
// also for PHOS acceptance 
//--
//-- by Renzhuo Wan (Iopp-wuhan) May 13,2009
//_________________________________________________________________________

//Root
class TList;
class TH3F ;
class TH2F ;
class TLorentzVector;
//Analysis
class AliAODEvent ;
class AliESDEvent ;
#include "AliAnaPartCorrBaseClass.h"
class TParticle;

#ifdef __PHOSUTIL__
  class AliPHOSGeoUtils;
#endif

#ifdef __EMCALUTIL__
  class AliEMCALGeoUtils;
#endif


class AliAnaNeutralMeson : public AliAnaPartCorrBaseClass {
  
  public: 
  
  AliAnaNeutralMeson() ; // default ctor
  AliAnaNeutralMeson(const char *name) ; // default ctor
  AliAnaNeutralMeson(const AliAnaNeutralMeson & g) ; // cpy ctor
  AliAnaNeutralMeson & operator = (const AliAnaNeutralMeson & api0) ;//cpy assignment
  virtual ~AliAnaNeutralMeson() ;//virtual dtor
  
  //......
  TList * GetCreateOutputObjects(); 
  void Print(const Option_t * opt) const;
  
  void Init();
  void InitParameters();
  void MakeAnalysisFillHistograms();
  //void Terminate(TList * outList);
  //......
  void Get2GammaAsyDistPtM(AliAODPWG4Particle *p1, AliAODPWG4Particle *p2, Double_t &asy, Double_t &dist, 
                         Double_t &pt, Double_t &mass); //get dist,asy, pair pt, pairmass btween two ptc
  void GetAsyDistPtM(TLorentzVector photon1, TLorentzVector photon2, Double_t &asy, Double_t &dist, 
                   Double_t &pt, Double_t &mass); //get dist,asy, pair pt, pairmass btween two ptc
  Bool_t IsPhotonSelected(AliAODPWG4Particle *p1, Int_t ipid, Int_t nmod); //select the photon based on pid and nmod to be analyzed

  void RealPi0Eta(TClonesArray * array) ; //perform the 2gamma IVM analysis for one event to extract pi0 and eta
  void MixPi0Eta(TClonesArray * array1, TClonesArray *array2) ; //perform the 2gamma IVM analysis for mixed eventsto reconstruct the background
  void RealOmega(TClonesArray *array) ; //perform the pi0+gamma->3gamma IVM analysis for one event to extract omega
  void MixOmegaAB(TClonesArray * array1, TClonesArray * array2); ////perform the 2gamma IVM analysis for mixed events to reconstruct the bacground (r1_event1+r2_event1)+r3_event2 and (r1_event1+r2_event2)+r3_event2
  void MixOmegaC(TClonesArray * array1, TClonesArray * array2, TClonesArray * array3); //perform the 2gamma IVM analysis for mixed events to reconstruct the bacground (r1_event1+r2_event2)+r3_event3
      
  void PhotonAcceptance(AliStack *stack); //to get MC photon pt and photon can be hit on dectector
  void Pi0EtaAcceptance(AliStack *stack); //to get MC pi0 and eta pt and photon can be hit on dectector
  void OmegaAcceptance(AliStack * stack); //to get MC omega pt and photon can be hit on dectector
      
  void GetMCDistAsy(TParticle * p1, TParticle * p2, Double_t & dist, Double_t &asy); //to get MC information between two particle
  
  void SetNAsyBinsMinMax(Int_t bins, Double_t min, Double_t max) {fNbinsAsy=bins; fMinAsy=min; fMaxAsy=max; }   //set Asy bins, min and max
  void SetNPtBinsMinMax(Int_t bins, Double_t min, Double_t max) {fNbinsPt=bins; fMinPt=min; fMaxPt=max; } //set pt bins, min and max
  void SetNMassBinsMinMas(Int_t bins, Double_t min, Double_t max) {fNbinsM=bins; fMinM=min; fMaxM=max; } //set mass pt bins, min and max
  void SetNEventsMixed(Int_t nevents) { fNmaxMixEv=nevents;} //events to be mixed 
  
  TString GetCalorimeter()   const { return fCalorimeter ; } //PHOS or EMCAL
  void SetCalorimeter(TString det);// { fCalorimeter = det ; } 
 
  void SetNPID(Int_t pid) {fNPID=pid;}  //number of PID cut
  void SetInvMassCut(Int_t mass) {fInvMassCut = mass;}  // 
  void SetPi0MassPeakWidthCut(Double_t cut) {fPi0MassPeakWidthCut = cut;}
 
  Bool_t SetAnaPi0Eta(Bool_t ana)  {fAnaPi0Eta = ana; return fAnaPi0Eta;} //analysis pi0 and eta ????
  Bool_t SetAnaOmega(Bool_t ana)   {fAnaOmega = ana; return fAnaOmega ;}  //analysis omega ????
 
  Bool_t IsBadRun(Int_t /*iRun*/) const {return kFALSE;} //Tests if this run bad according to private list
  
  private:

  Bool_t fAnaPi0Eta; //whether ana the 2gamma IVM to extract pi0 and eta spectrum
  Bool_t fAnaOmega;  //whether ana the 3gamma IVM to extract omega spectrum

  Int_t fNPID;       //PID bin
  Int_t fNmaxMixEv ; //number events to be mixed
  Int_t fNAsy;       //number of Asy cut
  Int_t fNMod;       //number of PHOS or EMCAL modules to be analyzed
  Int_t fNHistos;    //number of histograms
  Double_t fPerPtBin;//plot the IVM distribution per pt bin
 
  Double_t *fAsyCut; // asy cut(< 0.7, 0.8, 0.9, 1)
  Int_t    *fModCut; // mod cut(1,3,5) for PHOS, and 0 (default) for EMCAL

  Double_t fDist; //distance between two particle in (eta, phi) plane
  Double_t fAsy;  //asysmetry |(e1-e2)|/(e1+e2)
  Double_t fPt;   //pair pt
  Double_t fMass; //pair mass
  Double_t fInvMassCut; // cut of invmass

  Int_t fNbinsAsy; //Asy bin number, min and max
  Double_t fMinAsy;
  Double_t fMaxAsy;
  Int_t fNbinsPt;  //Pt bin number, min and max
  Double_t fMinPt;
  Double_t fMaxPt;
  Int_t fNbinsM;  //mass bin number, min and max
  Double_t fMinM;
  Double_t fMaxM;

  //! containers for photons in stored events
  TList * fEventsListPi0Eta ; //! events with cluster larger than 2 for pi0 and eta analysis
  TList * fEventsListOmega ; //! events with cluster larger than 3 for omega analysis

  Double_t fPi0Mass;    //nominal pi0 mass
  Double_t fPi0MassPeakWidthCut; //pi0 candidate
  
  TH1I *  fhNClusters;  //cluster multiplicity distribution per event
  TH1D ** fhRecPhoton; //reconstructed photon pt distribution
  TH2F ** fhRecPhotonEtaPhi; //photon (eta, phi) distruction

  TH1F ** fReal2Gamma; //real 2gamma IVM in per pt bin (ipid, imod, iasy)
  TH1F ** fMix2Gamma;  //mix 2gamma  IVM in per pt bin (ipid, imod, iasy, pt)
  TH1F ** fRealOmega;  //real omega  IVM in per pt bin (ipid, imod, iasy, pt)
  TH1F ** fMixOmegaA;  //mixA omega  IVM in per pt bin (ipid, imod, iasy, pt)  (r1_event1+r2_event1)+r3_event2
  TH1F ** fMixOmegaB;  //mixB omega  IVM in per pt bin (ipid, imod, iasy, pt)  (r1_event1+r2_event2)+r3_event2
  TH1F ** fMixOmegaC;  //mixC omega  IVM in per pt bin (ipid, imod, iasy, pt)  (r1_event1+r2_event2)+r3_event3 

  TH3F ** fRealTwoGammaAsyPtM; //real 2gamma IVM(asy, pt, m) without PID, mod cut
  TH3F ** fMixTwoGammaAsyPtM;  //mix  2gamma IVM(asy, pt, m) without PID, mod cut
  TH3F ** fRealPi0GammaAsyPtM; //real omega  IVM(asy, pt, m) without PID, mod cut
  TH3F ** fMixAPi0GammaAsyPtM; //mixA omega  IVM(asy, pt, m) without PID, mod cut
  TH3F ** fMixBPi0GammaAsyPtM; //mixB omega  IVM(asy, pt, m) without PID, mod cut
  TH3F ** fMixCPi0GammaAsyPtM; //mixC omega  IVM(asy, pt, m) without PID, mod cut

  //Histograms filled only if MC data is requested       
  //pt distribution
  //currently only PHOS included 
  TH1D * fhPrimPhotonPt ;    //! Spectrum of Primary 
  TH1D * fhPrimPhotonAccPt ; //! Spectrum of primary with accepted daughters 
  TH1D * fhPrimPhotonY ;     //! Rapidity distribution of primary particles
  TH1D * fhPrimPhotonAccY ;  //! Rapidity distribution of primary with accepted daughters
  TH1D * fhPrimPhotonPhi ;   //! Azimutal distribution of primary particles
  TH1D * fhPrimPhotonAccPhi; //! Azimutal distribution of primary with accepted daughters 
 
  TH1D * fhPrimPi0Pt ;    //! Spectrum of Primary 
  TH1D * fhPrimPi0AccPt ; //! Spectrum of primary with accepted daughters 
  TH1D * fhPrimPi0Y ;     //! Rapidity distribution of primary particles
  TH1D * fhPrimPi0AccY ;  //! Rapidity distribution of primary with accepted daughters
  TH1D * fhPrimPi0Phi ;   //! Azimutal distribution of primary particles
  TH1D * fhPrimPi0AccPhi; //! Azimutal distribution of primary with accepted daughters  
 
  TH1D * fhPrimEtaPt ;    //! Spectrum of Primary 
  TH1D * fhPrimEtaAccPt ; //! Spectrum of primary with accepted daughters 
  TH1D * fhPrimEtaY ;     //! Rapidity distribution of primary particles
  TH1D * fhPrimEtaAccY ;  //! Rapidity distribution of primary with accepted daughters
  TH1D * fhPrimEtaPhi ;   //! Azimutal distribution of primary particles
  TH1D * fhPrimEtaAccPhi; //! Azimutal distribution of primary with accepted daughters
 
  TH1D * fhPrimOmegaPt ;    //! Spectrum of Primary 
  TH1D * fhPrimOmegaAccPt ; //! Spectrum of primary with accepted daughters 
  TH1D * fhPrimOmegaY ;     //! Rapidity distribution of primary particles
  TH1D * fhPrimOmegaAccY ;  //! Rapidity distribution of primary with accepted daughters
  TH1D * fhPrimOmegaPhi ;   //! Azimutal distribution of primary particles
  TH1D * fhPrimOmegaAccPhi; //! Azimutal distribution of primary with accepted daughters

  TString fCalorimeter ;     //Select Calorimeter       
    
#ifdef __PHOSUTIL__
  AliPHOSGeoUtils * fPHOSGeo ; //! PHOS geometry pointer
#endif	
#ifdef __EMCALUTIL__
  AliEMCALGeoUtils * fEMCALGeo ; //! EMCAL geometry pointer
#endif	
  
  ClassDef(AliAnaNeutralMeson,1)
} ;


#endif //ALIANANEUTRALMESON



