#ifndef AliJetContainer_H
#define AliJetContainer_H

class TClonesArray;
class TString;

class AliEMCALGeometry;
class AliEmcalJet;
class AliVEvent;
class AliParticleContainer;
class AliClusterContainer;
class AliLocalRhoParameter;

#include <TMath.h>
#include <TLorentzVector.h>
#include "AliRhoParameter.h"
#include "AliEmcalContainer.h"
#include "AliLog.h"
#include "AliVEvent.h"

class AliJetContainer : public AliEmcalContainer {
 public:
 
  enum EJetType_t {
    kFullJet,
    kChargedJet,
    kNeutralJet
  };

  enum EJetAlgo_t {
    kt_algorithm            = 0,
    cambridge_algorithm     = 1,
    antikt_algorithm        = 2,
    undefined_jet_algorithm = 999
  };

  enum ERecoScheme_t {
    E_scheme        = 0,
    pt_scheme       = 1,
    pt2_scheme      = 2,
    Et_scheme       = 3,
    Et2_scheme      = 4,
    BIpt_scheme     = 5,
    BIpt2_scheme    = 6,
    external_scheme = 99
  };

  enum JetAcceptanceType {
    kTPC       ,     // TPC acceptance
    kTPCfid    ,     // TPC fiducial acceptance
    kEMCAL     ,     // EMCal acceptance
    kEMCALfid  ,     // EMCal fiducial acceptance
    kDCAL      ,     // DCal acceptance
    kDCALfid   ,     // DCal fiducial acceptance
    kUser            // User defined acceptance
  };

  AliJetContainer();
  AliJetContainer(const char *name); 
  AliJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag);
  virtual ~AliJetContainer() {;}
  
  void LoadRho(AliVEvent *event);
  void LoadLocalRho(AliVEvent *event);
  void LoadRhoMass(AliVEvent *event);

  void                        SetJetAcceptanceType(JetAcceptanceType type)         { fJetAcceptanceType          = type ; }
  void                        PrintCuts();
  void                        ResetCuts();
  void                        SetJetEtaPhiEMCAL(Double_t r=0.) ;
  void                        SetJetEtaPhiDCAL(Double_t r=0.)  ;
  void                        SetJetEtaPhiTPC(Double_t r=0.)   ;
  void                        SetRunNumber(Int_t r)                                { fRunNumber = r;                      }
  void                        SetJetEtaLimits(Float_t min, Float_t max)            { fJetMinEta = min, fJetMaxEta = max ; }
  void                        SetJetPhiLimits(Float_t min, Float_t max, Float_t offset=0.);
  void                        SetJetPtCut(Float_t cut)                             { fJetPtCut       = cut              ; }
  void                        SetJetPtCutMax(Float_t cut)                          { fJetPtCutMax    = cut              ; }
  void                        SetJetRadius(Float_t r)                              { fJetRadius      = r                ; } 
  void                        SetJetAreaCut(Float_t cut)                           { fJetAreaCut     = cut              ; }
  void                        SetPercAreaCut(Float_t p)                            { if(fJetRadius==0.) AliWarning("JetRadius not set. Area cut will be 0"); 
                                                                                     fJetAreaCut = p*TMath::Pi()*fJetRadius*fJetRadius; }
  void                        SetAreaEmcCut(Double_t a = 0.99)                     { fAreaEmcCut     = a                ; }
  void                        SetZLeadingCut(Float_t zemc, Float_t zch)            { fZLeadingEmcCut = zemc; fZLeadingChCut = zch ; }
  void                        SetNEFCut(Float_t min = 0., Float_t max = 1.)        { fNEFMinCut = min; fNEFMaxCut = max;  }
  void                        SetFlavourCut(Int_t myflavour)                       { fFlavourSelection = myflavour;}
  void                        SetMaxClusterPt(Float_t b)                           { fMaxClusterPt   = b                ; }
  void                        SetMaxTrackPt(Float_t b)                             { fMaxTrackPt     = b                ; }
  void                        SetPtBiasJetClus(Float_t b)                          { fPtBiasJetClus  = b                ; }
  void                        SetNLeadingJets(Int_t t)                             { fNLeadingJets   = t                ; }
  void                        SetMinNConstituents(Int_t n)                         { fMinNConstituents = n              ; }
  void                        SetPtBiasJetTrack(Float_t b)                         { fPtBiasJetTrack = b                ; }
  void                        SetLeadingHadronType(Int_t t)                        { fLeadingHadronType = t             ; }
  void                        SetJetTrigger(UInt_t t=AliVEvent::kEMCEJE)           { fJetTrigger     = t                ; }
  void                        SetTagStatus(Int_t i)                                { fTagStatus      = i                ; }

  void                        SetRhoName(const char *n)                            { fRhoName        = n                ; }
  void                        SetLocalRhoName(const char *n)                       { fLocalRhoName   = n                ; }
  void                        SetRhoMassName(const char *n)                        { fRhoMassName    = n                ; }
    
  void                        ConnectParticleContainer(AliParticleContainer *c)    { fParticleContainer = c             ; }
  void                        ConnectClusterContainer(AliClusterContainer *c)      { fClusterContainer  = c             ; }

  AliEmcalJet                *GetLeadingJet(const char* opt="")          ;
  AliEmcalJet                *GetJet(Int_t i)                       const;
  AliEmcalJet                *GetAcceptJet(Int_t i)                      ;
  AliEmcalJet                *GetJetWithLabel(Int_t lab)            const;
  AliEmcalJet                *GetAcceptJetWithLabel(Int_t lab)           ;
  AliEmcalJet                *GetNextAcceptJet(Int_t i=-1)               ;
  AliEmcalJet                *GetNextJet(Int_t i=-1)                     ;
  Bool_t                      GetMomentum(TLorentzVector &mom, Int_t i);
  Bool_t                      GetAcceptMomentum(TLorentzVector &mom, Int_t i);
  Bool_t                      GetNextMomentum(TLorentzVector &mom, Int_t i=-1);
  Bool_t                      GetNextAcceptMomentum(TLorentzVector &mom, Int_t i=-1);
  Bool_t                      AcceptJet(const AliEmcalJet* jet)          ;
  Bool_t                      AcceptBiasJet(const AliEmcalJet* jet)      ;
  Int_t                       GetFlavourCut()                       const    {return fFlavourSelection;}
  Int_t                       GetNJets()                            const    {return GetNEntries();}
  Int_t                       GetNAcceptedJets()                         ;

  Double_t                    GetLeadingHadronPt(const AliEmcalJet* jet)  const;
  void                        GetLeadingHadronMomentum(TLorentzVector &mom, const AliEmcalJet* jet)  const;
  Double_t                    GetZ(const AliEmcalJet *jet, TLorentzVector mom) const;
  Double_t                    GetZLeadingEmc(const AliEmcalJet *jet)      const;
  Double_t                    GetZLeadingCharged(const AliEmcalJet *jet)  const;
  AliRhoParameter            *GetRhoParameter()                              {return fRho;}
  Double_t                    GetRhoVal()                           const    {if (fRho) return fRho->GetVal(); else return 0;}
  const TString&              GetRhoName()                          const    {return fRhoName;}
  AliLocalRhoParameter       *GetLocalRhoParameter()                const    {return fLocalRho;}
  const TString&              GetLocalRhoName()                     const    {return fLocalRhoName;}
  AliRhoParameter            *GetRhoMassParameter()                          {return fRhoMass;}
  Double_t                    GetRhoMassVal()                       const    {if (fRhoMass) return fRhoMass->GetVal(); else return 0;}
  const TString&              GetRhoMassName()                      const    {return fRhoMassName;}
  Double_t                    GetJetPtCorr(Int_t i)                 const;
  Double_t                    GetJetPtCorrLocal(Int_t i)            const;
  Float_t                     GetJetRadius()                        const    {return fJetRadius;}
  Float_t                     GetJetEtaMin()                        const    {return fJetMinEta;}
  Float_t                     GetJetEtaMax()                        const    {return fJetMaxEta;}
  Float_t                     GetJetPhiMin()                        const    {return fJetMinPhi;}
  Float_t                     GetJetPhiMax()                        const    {return fJetMaxPhi;}
  Float_t                     GetJetPtCut()                         const    {return fJetPtCut;}
  Float_t                     GetJetPtCutMax()                      const    {return fJetPtCutMax;}

  void                        SetClassName(const char *clname);
  void                        SetArray(AliVEvent *event);
  AliParticleContainer       *GetParticleContainer() const                   {return fParticleContainer;}
  AliClusterContainer        *GetClusterContainer() const                    {return fClusterContainer;}
  Double_t                    GetFractionSharedPt(const AliEmcalJet *jet, AliParticleContainer *cont2 = 0x0) const;

  const char*                 GetTitle() const;

  static TString              GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag);

 protected:
  void SetEMCALGeometry();
  
  JetAcceptanceType           fJetAcceptanceType;    //  acceptance type
  Float_t                     fJetRadius;            //  jet radius
  TString                     fRhoName;              //  Name of rho object
  TString                     fLocalRhoName;         //  Name of local rho object
  TString                     fRhoMassName;          //  Name of rho mass object
  Int_t                       fFlavourSelection;     //  selection on jet flavour
  Float_t                     fPtBiasJetTrack;       //  select jets with a minimum pt track
  Float_t                     fPtBiasJetClus;        //  select jets with a minimum pt cluster
  Float_t                     fJetPtCut;             //  cut on jet pt
  Float_t                     fJetPtCutMax;          //  cut on jet pt - MAX
  Float_t                     fJetAreaCut;           //  cut on jet area
  Float_t                     fAreaEmcCut;           //  minimum cut on jet emcal area
  Float_t                     fJetMinEta;            //  minimum eta jet acceptance
  Float_t                     fJetMaxEta;            //  maximum eta jet acceptance
  Float_t                     fJetMinPhi;            //  minimum phi jet acceptance
  Float_t                     fJetMaxPhi;            //  maximum phi jet acceptance  
  Float_t                     fPhiOffset;            //  offset to allow cutting
  Float_t                     fMaxClusterPt;         //  maximum cluster constituent pt to accept the jet
  Float_t                     fMaxTrackPt;           //  maximum track constituent pt to accept the jet
  Float_t                     fZLeadingEmcCut;       //  maximum z,leading neutral
  Float_t                     fZLeadingChCut;        //  maximum z,leading charged
  Float_t                     fNEFMinCut;            //  minimum NEF in a jet
  Float_t                     fNEFMaxCut;            //  maximum NEF in a jet
  Int_t                       fLeadingHadronType;    //  0 = charged, 1 = neutral, 2 = both
  Int_t                       fNLeadingJets;         //  how many jets are to be considered the leading jet(s)
  Int_t                       fMinNConstituents;     //  minimum number of constituents in jet
  UInt_t                      fJetTrigger;           //  jet trigger
  Int_t                       fTagStatus;            //  jet tag status
  AliParticleContainer       *fParticleContainer;    //  particle container (jet constituents)
  AliClusterContainer        *fClusterContainer;     //  cluster container (jet constituents)
  AliRhoParameter            *fRho;                  //! event rho for these jets
  AliLocalRhoParameter       *fLocalRho;             //! event local rho for these jets
  AliRhoParameter            *fRhoMass;              //! event rho mass for these jets
  AliEMCALGeometry           *fGeom;                 //! emcal geometry
  Int_t                       fRunNumber;            //! run number

 private:
  AliJetContainer(const AliJetContainer& obj); // copy constructor
  AliJetContainer& operator=(const AliJetContainer& other); // assignment

  ClassDef(AliJetContainer,14);
};

#endif
