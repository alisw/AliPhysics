#ifndef AliJetContainer_H
#define AliJetContainer_H

//
// container with name, TClonesArray and cuts for jets
//

class TClonesArray;
class TString;

class AliEMCALGeometry;
class AliEmcalJet;
class AliVEvent;
class AliParticleContainer;
class AliClusterContainer;

#include "AliRhoParameter.h"

#include "AliEmcalContainer.h"

class AliJetContainer : public AliEmcalContainer {
 public:
 
  enum JetAcceptanceType {
    kTPC       = 0,     // TPC acceptance
    kEMCAL     = 1,     // EMCal acceptance
    kUser      = 2      // User defined acceptance
  };


  AliJetContainer();
  AliJetContainer(const char *name); 
  virtual ~AliJetContainer() {;}

  void SetEMCALGeometry();
  void SetEMCALGeometry(AliEMCALGeometry *p) {fGeom = p;}
  void LoadRho(AliVEvent *event);

  void                        SetJetAcceptanceType(JetAcceptanceType type)         { fJetAcceptanceType          = type ; }

  void                        ResetCuts();

  void                        SetJetEtaPhiEMCAL() ;
  void                        SetJetEtaPhiTPC()   ;
  void                        SetRunNumber(Int_t r)                                { fRunNumber = r;                      }

  void                        SetJetEtaLimits(Float_t min, Float_t max)            { fJetMinEta = min, fJetMaxEta = max ; }
  void                        SetJetPhiLimits(Float_t min, Float_t max)            { fJetMinPhi = min, fJetMaxPhi = max ; }
  void                        SetJetAreaCut(Float_t cut)                           { fJetAreaCut     = cut              ; }
  void                        SetPercAreaCut(Float_t p)                            { if(fJetRadius==0.) AliWarning("JetRadius not set. Area cut will be 0"); fJetAreaCut = p*TMath::Pi()*fJetRadius*fJetRadius; }

  void                        SetAreaEmcCut(Double_t a = 0.99)                     { fAreaEmcCut     = a                ; }
  void                        SetJetPtCut(Float_t cut)                             { fJetPtCut       = cut              ; }
  void                        SetJetRadius(Float_t r)                              { fJetRadius      = r                ; } 
  virtual void                SetRhoName(const char *n)                            { fRhoName        = n                ; }
  void                        SetMaxClusterPt(Float_t b)                           { fMaxClusterPt   = b                ; }
  void                        SetMaxTrackPt(Float_t b)                             { fMaxTrackPt     = b                ; }
  void                        SetPtBiasJetClus(Float_t b)                          { fPtBiasJetClus  = b                ; }
  void                        SetNLeadingJets(Int_t t)                             { fNLeadingJets   = t                ; }
  void                        SetPtBiasJetTrack(Float_t b)                         { fPtBiasJetTrack = b                ; }

  void                        SetLeadingHadronType(Int_t t)                        { fLeadingHadronType = t             ; }
  void                        SetJetBitMap(UInt_t m)                               { fJetBitMap      = m                ; }
  void                        SetJetTrigger(UInt_t t=AliVEvent::kEMCEJE)           { fJetTrigger     = t                ; }

  void                        ConnectParticleContainer(AliParticleContainer *c)    { fParticleContainer = c             ; }
  void                        ConnectClusterContainer(AliClusterContainer *c)      { fClusterContainer  = c             ; }

  AliEmcalJet                *GetLeadingJet(const char* opt="")          ;
  AliEmcalJet                *GetJet(Int_t i)                       const;
  AliEmcalJet                *GetAcceptJet(Int_t i)                 const;
  AliEmcalJet                *GetJetWithLabel(Int_t lab)            const;
  AliEmcalJet                *GetAcceptJetWithLabel(Int_t lab)      const;
  AliEmcalJet                *GetNextAcceptJet(Int_t i=-1)               ;
  AliEmcalJet                *GetNextJet(Int_t i=-1)                     ;
  void                        GetMomentum(TLorentzVector &mom, Int_t i) const;
  Bool_t                      AcceptJet(AliEmcalJet* jet)           const;
  Bool_t                      AcceptBiasJet(AliEmcalJet* jet)       const; 
  Int_t                       GetNJets()                            const    {return GetNEntries();}
  Double_t                    GetLeadingHadronPt(AliEmcalJet* jet)  const;
  void                        GetLeadingHadronMomentum(TLorentzVector &mom, AliEmcalJet* jet)  const;
  AliRhoParameter            *GetRhoParameter()                              {return fRho;}
  Double_t                    GetRhoVal()                           const    { if (fRho) return fRho->GetVal(); else return 0;}
  const TString&              GetRhoName()                          const    {return fRhoName;}
  Double_t                    GetJetPtCorr(Int_t i)                 const;
  Float_t                     GetJetRadius()                        const    {return fJetRadius;}
  Float_t                     GetJetEtaMin()                        const    {return fJetMinEta;}
  Float_t                     GetJetEtaMax()                        const    {return fJetMaxEta;}
  Float_t                     GetJetPhiMin()                        const    {return fJetMinPhi;}
  Float_t                     GetJetPhiMax()                        const    {return fJetMaxPhi;}
  void                        SetClassName(const char *clname);
  void                        SetArray(AliVEvent *event);
  AliParticleContainer       *GetParticleContainer()                         {return fParticleContainer;}
  AliClusterContainer        *GetClusterContainer()                          {return fClusterContainer;}

 protected:
  JetAcceptanceType           fJetAcceptanceType;    //  acceptance type
  Float_t                     fJetRadius;            //  jet radius
  TString                     fRhoName;              //  Name of rho object
  Float_t                     fPtBiasJetTrack;       //  select jets with a minimum pt track
  Float_t                     fPtBiasJetClus;        //  select jets with a minimum pt cluster
  Float_t                     fJetPtCut;             //  cut on jet pt
  Float_t                     fJetAreaCut;           //  cut on jet area
  Float_t                     fAreaEmcCut;           //  minimum cut on jet emcal area
  Float_t                     fJetMinEta;            //  minimum eta jet acceptance
  Float_t                     fJetMaxEta;            //  maximum eta jet acceptance
  Float_t                     fJetMinPhi;            //  minimum phi jet acceptance
  Float_t                     fJetMaxPhi;            //  maximum phi jet acceptance  
  Float_t                     fMaxClusterPt;         //  maximum cluster constituent pt to accept the jet
  Float_t                     fMaxTrackPt;           //  maximum track constituent pt to accept the jet
  Int_t                       fLeadingHadronType;    //  0 = charged, 1 = neutral, 2 = both
  Int_t                       fNLeadingJets;         //  how many jets are to be considered the leading jet(s)
  UInt_t                      fJetBitMap;            //  bit map of accepted jets
  UInt_t                      fJetTrigger;           //  jet trigger
  AliParticleContainer       *fParticleContainer;    //  particle container (jet constituents)
  AliClusterContainer        *fClusterContainer;     //  cluster container (jet constituents)

  AliRhoParameter            *fRho;                  //! event rho for these jets
  AliEMCALGeometry           *fGeom;                 //! emcal geometry
  Int_t                       fRunNumber;            //! run number

 private:
  AliJetContainer(const AliJetContainer& obj); // copy constructor
  AliJetContainer& operator=(const AliJetContainer& other); // assignment

  ClassDef(AliJetContainer,3);

};

#endif


