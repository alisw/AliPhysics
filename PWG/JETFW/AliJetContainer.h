#ifndef ALIJETCONTAINER_H
#define ALIJETCONTAINER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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
#include "AliParticleContainer.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliEmcalJet.h"

#if !(defined(__CINT__) || defined(__MAKECINT__))
typedef EMCALIterableContainer::AliEmcalIterableContainerT<AliEmcalJet, EMCALIterableContainer::operator_star_object<AliEmcalJet> > AliJetIterableContainer;
typedef EMCALIterableContainer::AliEmcalIterableContainerT<AliEmcalJet, EMCALIterableContainer::operator_star_pair<AliEmcalJet> > AliJetIterableMomentumContainer;
#endif

/**
 * @class AliJetContainer
 * @brief Container for jet within the EMCAL jet framework
 * @author Marta Verweij
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 *
 * Container with name, TClonesArray and cuts for jets
 */
class AliJetContainer : public AliParticleContainer {
 public:
 
  enum EJetType_t {
    kFullJet,
    kChargedJet,
    kNeutralJet
  };

  enum EJetAlgo_t {
    kt_algorithm                    = 0,
    antikt_algorithm                = 1,
    cambridge_algorithm             = 2,
    genkt_algorithm                 = 3,
    cambridge_for_passive_algorithm = 11,
    genkt_for_passive_algorithm     = 13,
    plugin_algorithm                = 99,
    undefined_jet_algorithm         = 999
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

  /**
   * @enum JetAcceptanceType
   * @brief Bit definition for jet geometry acceptance. Defined here for backwards compatibility. This will be 
   * removed. Please use AliEmcalJet::JetAcceptanceType in your code.
   */
  enum JetAcceptanceType {
    kTPC              = AliEmcalJet::kTPC,          ///< TPC acceptance
    kTPCfid           = AliEmcalJet::kTPCfid,       ///< TPC fiducial acceptance (each eta edge narrowed by jet R)
    kEMCAL            = AliEmcalJet::kEMCAL,        ///< EMCal acceptance
    kEMCALfid         = AliEmcalJet::kEMCALfid,     ///< EMCal fiducial acceptance (each eta, phi edge narrowed by jet R)
    kDCAL             = AliEmcalJet::kDCAL,         ///< DCal acceptance -- spans entire rectangular region in eta-phi (including most of PHOS)
    kDCALfid          = AliEmcalJet::kDCALfid,      ///< DCal fiducial acceptance (each eta, phi edge narrowed by jet R)
    kDCALonly         = AliEmcalJet::kDCALonly,     ///< DCal acceptance -- spans ONLY DCal (no PHOS or gap)
    kDCALonlyfid      = AliEmcalJet::kDCALonlyfid,  ///< DCal fiducial acceptance (each eta, phi edge narrowed by jet R)
    kPHOS             = AliEmcalJet::kPHOS,         ///< PHOS acceptance
    kPHOSfid          = AliEmcalJet::kPHOSfid,      ///< PHOS fiducial acceptance (each eta, phi edge narrowed by jet R)
    kUser             = AliEmcalJet::kUser          ///< Full acceptance, i.e. no acceptance cut applied -- left to user
  };

  AliJetContainer();
  AliJetContainer(const char *name);
  AliJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag = "Jet");
  virtual ~AliJetContainer() {;}
  
  void LoadRho(const AliVEvent *event);
  void LoadLocalRho(const AliVEvent *event);
  void LoadRhoMass(const AliVEvent *event);

  void                        SetJetAcceptanceType(UInt_t type)         { fJetAcceptanceType          = type ; }
  void                        PrintCuts();
  void                        ResetCuts();
  void                        SetJetEtaLimits(Float_t min, Float_t max)            { SetEtaLimits(min, max)             ; }
  void                        SetJetPhiLimits(Float_t min, Float_t max)            { SetPhiLimits(min, max)             ; }
  void                        SetJetPtCut(Float_t cut)                             { SetMinPt(cut)                      ; }
  void                        SetJetPtCutMax(Float_t cut)                          { SetMaxPt(cut)                      ; }
  void                        SetRunNumber(Int_t r)                                { fRunNumber = r;                      }
  void                        SetJetRadius(Float_t r)                              { fJetRadius      = r                ; } 
  void                        SetJetAreaCut(Float_t cut)                           { fJetAreaCut     = cut              ; }
  void                        SetPercAreaCut(Float_t p)                            { if(fJetRadius==0.) AliWarning("JetRadius not set. Area cut will be 0"); 
                                                                                     fJetAreaCut = p*TMath::Pi()*fJetRadius*fJetRadius; }
  void                        SetAreaEmcCut(Double_t a = 0.99)                     { fAreaEmcCut     = a                ; }
  void                        SetZLeadingCut(Float_t zemc, Float_t zch)            { fZLeadingEmcCut = zemc; fZLeadingChCut = zch ; }
  void                        SetNEFCut(Float_t min = 0., Float_t max = 1.)        { fNEFMinCut = min; fNEFMaxCut = max;  }
  void                        SetFlavourCut(Int_t myflavour)                       { fFlavourSelection = myflavour;}
  void                        SetMinClusterPt(Float_t b)                           { fMinClusterPt   = b                ; }
  void                        SetMaxClusterPt(Float_t b)                           { fMaxClusterPt   = b                ; }
  void                        SetMinTrackPt(Float_t b)                             { fMinTrackPt     = b                ; }
  void                        SetMaxTrackPt(Float_t b)                             { fMaxTrackPt     = b                ; }
  void                        SetPtBiasJetClus(Float_t b)                          { SetMinClusterPt(b)                 ; }
  void                        SetNLeadingJets(Int_t t)                             { fNLeadingJets   = t                ; }
  void                        SetMinNConstituents(Int_t n)                         { fMinNConstituents = n              ; }
  void                        SetPtBiasJetTrack(Float_t b)                         { SetMinTrackPt(b)                   ; }
  void                        SetLeadingHadronType(Int_t t)                        { fLeadingHadronType = t             ; }
  void                        SetJetTrigger(UInt_t t=AliVEvent::kEMCEJE)           { fJetTrigger     = t                ; }
  void                        SetTagStatus(Int_t i)                                { fTagStatus      = i                ; }

  void                        SetRhoName(const char *n)                            { fRhoName        = n                ; }
  void                        SetLocalRhoName(const char *n)                       { fLocalRhoName   = n                ; }
  void                        SetRhoMassName(const char *n)                        { fRhoMassName    = n                ; }
    
  void                        SetTpcHolePos(Double_t b)                                {fTpcHolePos       =   b     ;}
  void                        SetTpcHoleWidth(Double_t b)                             {fTpcHoleWidth    =   b     ;} 


  void                        ConnectParticleContainer(AliParticleContainer *c)    { fParticleContainer = c             ; }
  void                        ConnectClusterContainer(AliClusterContainer *c)      { fClusterContainer  = c             ; }

  AliEmcalJet                *GetLeadingJet(const char* opt="")          ;
  AliEmcalJet                *GetJet(Int_t i)                       const;
  AliEmcalJet                *GetAcceptJet(Int_t i)                 const;
  AliEmcalJet                *GetNextAcceptJet()                         ;
  AliEmcalJet                *GetNextJet()                               ;
  Bool_t                      GetMomentumFromJet(TLorentzVector &mom, const AliEmcalJet* jet, Double_t mass) const;
  Bool_t                      GetMomentumFromJet(TLorentzVector &mom, const AliEmcalJet* jet) const;
  Bool_t                      GetMomentum(TLorentzVector &mom, Int_t i) const;
  Bool_t                      GetAcceptMomentum(TLorentzVector &mom, Int_t i) const;
  Bool_t                      GetNextMomentum(TLorentzVector &mom);
  Bool_t                      GetNextAcceptMomentum(TLorentzVector &mom);
  virtual Bool_t              AcceptObject(Int_t i, UInt_t &rejectionReason) const { return AcceptJet(i, rejectionReason);}
  virtual Bool_t              AcceptObject(const TObject* obj, UInt_t &rejectionReason) const { return AcceptJet(dynamic_cast<const AliEmcalJet*>(obj), rejectionReason);}
  virtual Bool_t              AcceptJet(Int_t i, UInt_t &rejectionReason) const;
  virtual Bool_t              AcceptJet(const AliEmcalJet* jet, UInt_t &rejectionReason) const;
  virtual Bool_t              ApplyJetCuts(const AliEmcalJet* clus, UInt_t &rejectionReason) const;
  virtual Bool_t              CheckTpcHolesOverlap(const AliEmcalJet* clus,UInt_t &rejectionReason) const;
  Int_t                       GetFlavourCut()                       const    {return fFlavourSelection;}
  Int_t                       GetNJets()                            const    {return GetNEntries();}
  Int_t                       GetNAcceptedJets()                         ;

  Double_t                    GetLeadingHadronPt(const AliEmcalJet* jet)  const;
  void                        GetLeadingHadronMomentum(TLorentzVector &mom, const AliEmcalJet* jet)  const;
  Double_t                    GetZ(const AliEmcalJet *jet, const TLorentzVector& mom) const;
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
  Float_t                     GetJetRadius()                        const    {return fJetRadius ; }
  Double_t                    GetJetEtaMin()                        const    {return GetMinEta(); }
  Double_t                    GetJetEtaMax()                        const    {return GetMaxEta(); }
  Double_t                    GetJetPhiMin()                        const    {return GetMinPhi(); }
  Double_t                    GetJetPhiMax()                        const    {return GetMaxPhi(); }
  Double_t                    GetJetPtCut()                         const    {return GetMinPt() ; }
  Double_t                    GetJetPtCutMax()                      const    {return GetMaxPt() ; }

  void                        SetArray(const AliVEvent *event);
  AliParticleContainer       *GetParticleContainer() const                   {return fParticleContainer;}
  AliClusterContainer        *GetClusterContainer() const                    {return fClusterContainer;}
  Double_t                    GetFractionSharedPt(const AliEmcalJet *jet, AliParticleContainer *cont2 = 0x0) const;

  const char*                 GetTitle() const;

  static TString              GenerateJetName(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag);

#if !(defined(__CINT__) || defined(__MAKECINT__))
  const AliJetIterableContainer      all() const;
  const AliJetIterableContainer      accepted() const;

  const AliJetIterableMomentumContainer      all_momentum() const;
  const AliJetIterableMomentumContainer      accepted_momentum() const;
#endif

 protected:
  UInt_t                      fJetAcceptanceType;    ///  Jet acceptance type cut, see AliEmcalJet::JetAcceptanceType
  Float_t                     fJetRadius;            ///  jet radius
  TString                     fRhoName;              ///  Name of rho object
  TString                     fLocalRhoName;         ///  Name of local rho object
  TString                     fRhoMassName;          ///  Name of rho mass object
  Int_t                       fFlavourSelection;     ///  selection on jet flavour
  Float_t                     fJetAreaCut;           ///  cut on jet area
  Float_t                     fAreaEmcCut;           ///  minimum cut on jet emcal area
  Float_t                     fMinClusterPt;         ///  maximum cluster constituent pt to accept the jet
  Float_t                     fMaxClusterPt;         ///  maximum cluster constituent pt to accept the jet
  Float_t                     fMinTrackPt;           ///  maximum track constituent pt to accept the jet
  Float_t                     fMaxTrackPt;           ///  maximum track constituent pt to accept the jet
  Float_t                     fZLeadingEmcCut;       ///  maximum z,leading neutral
  Float_t                     fZLeadingChCut;        ///  maximum z,leading charged
  Float_t                     fNEFMinCut;            ///  minimum NEF in a jet
  Float_t                     fNEFMaxCut;            ///  maximum NEF in a jet
  Int_t                       fLeadingHadronType;    ///  0 = charged, 1 = neutral, 2 = both
  Int_t                       fNLeadingJets;         ///  how many jets are to be considered the leading jet(s)
  Int_t                       fMinNConstituents;     ///  minimum number of constituents in jet
  UInt_t                      fJetTrigger;           ///  jet trigger
  Int_t                       fTagStatus;            ///  jet tag status
  AliParticleContainer       *fParticleContainer;    ///  particle container (jet constituents)
  AliClusterContainer        *fClusterContainer;     ///  cluster container (jet constituents)
  AliRhoParameter            *fRho;                  //!<! event rho for these jets
  AliLocalRhoParameter       *fLocalRho;             //!<! event local rho for these jets
  AliRhoParameter            *fRhoMass;              //!<! event rho mass for these jets
  AliEMCALGeometry           *fGeom;                 //!<! emcal geometry
  Int_t                       fRunNumber;            //!<! run number
  Double_t                    fTpcHolePos;           ///   position(in radians) of the malfunctioning TPC sector
  Double_t                    fTpcHoleWidth;         ///   width of the malfunctioning TPC area
 private:
  AliJetContainer(const AliJetContainer& obj); // copy constructor
  AliJetContainer& operator=(const AliJetContainer& other); // assignment

  ClassDef(AliJetContainer, 18);
};

#endif
