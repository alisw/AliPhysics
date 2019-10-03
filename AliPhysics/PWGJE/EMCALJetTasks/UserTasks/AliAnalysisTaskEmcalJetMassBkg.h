#ifndef ALIANALYSISTASKEMCALJETMASSBKG_H
#define ALIANALYSISTASKEMCALJETMASSBKG_H

class TH1;
class TH2;
class TH3;
class TH3F;
class THnSparse;
class TClonesArray;
class TArrayI;
class TLorentzVector;
class AliAnalysisManager;
class AliJetContainer;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetMassBkg : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalJetMassBkg();
  AliAnalysisTaskEmcalJetMassBkg(const char *name);
  virtual ~AliAnalysisTaskEmcalJetMassBkg();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  //Setters
  void                        SetJetContainerBase(Int_t c)                         { fContainerBase           = c          ; }

  void                        SetJetMinRC2LJ(Float_t d)                            { fMinRC2LJ                = d          ; }
  void                        SetRCperEvent(Int_t n)                               { fRCperEvent              = n          ; }
  void                        SetConeRadius(Double_t r)                            { fConeRadius              = r          ; }
  void                        SetConeEtaPhiEMCAL() ;
  void                        SetConeEtaPhiTPC()   ;
  void                        SetConeEtaLimits(Float_t min, Float_t max)           { fConeMinEta = min, fConeMaxEta = max  ; }
  void                        SetConePhiLimits(Float_t min, Float_t max)           { fConeMinPhi = min, fConeMaxPhi = max  ; }

 protected:
  void                        ExecOnce();
  Bool_t                      RetrieveEventObjects();
  Bool_t                      Run();
  Bool_t                      FillHistograms();

  void                        GetRandomCone(TLorentzVector& lvRC, Float_t &pt, Float_t &eta, Float_t &phi, AliParticleContainer* tracks, AliClusterContainer* clusters, AliEmcalJet *jet = 0) const;
  void                        GetCone(TLorentzVector& lvRC,Float_t &pt, Float_t eta, Float_t phi, AliParticleContainer* tracks, AliClusterContainer* clusters) const;
  void                        GetPerpCone(TLorentzVector& lvRC, Float_t &pt, Float_t &eta, Float_t &phi, AliParticleContainer* tracks, AliClusterContainer* clusters, AliEmcalJet *jet = 0) const;
  TLorentzVector              GetSubtractedVector(Double_t pt, Double_t eta, Double_t phi, Double_t e);
  TLorentzVector              GetBkgVector(Double_t eta, Double_t phi, AliJetContainer *cont);

  Int_t                       fContainerBase;              // jets to be tagged

  Float_t                     fMinRC2LJ;                   // Minimum distance random cone to leading jet
  Int_t                       fRCperEvent;                 // No. of random cones per event
  Double_t                    fConeRadius;                 // Radius of the random cones
  Float_t                     fConeMinEta;                 // Minimum eta of the random cones
  Float_t                     fConeMaxEta;                 // Maximum eta of the random cones
  Float_t                     fConeMinPhi;                 // Minimum phi of the random cones
  Float_t                     fConeMaxPhi;                 // Maximum phi of the random cones
  
  AliJetContainer            *fJetsCont;                   //!Jets
  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters  

 private:
  TH2F            **fh2PtVsMassRC;                //!pT vs mass of RC
  TH3F            **fh2PtVsMassRCExLJDPhi;        //!pT vs mass of RC excluding area around leading jet
  TH2F            **fh2PtVsMassPerpConeLJ;        //!pT vs mass of cone perpendicular to leading jet
  TH2F            **fh2PtVsMassPerpConeTJ;        //!pT vs mass of cone perpendicular to all tagged jets

  TH2F            **fh2PtVsERC;                   //!E vs mass of RC
  TH3F            **fh2PtVsERCExLJDPhi;           //!E vs mass of RC excluding area around leading jet
  TH2F            **fh2PtVsEPerpConeLJ;           //!E vs mass of cone perpendicular to leading jet
  TH2F            **fh2PtVsEPerpConeTJ;           //!E vs mass of cone perpendicular to all tagged jets

  TProfile        **fpPtVsMassRC;                 //!pT vs Avg mass of RC
  TProfile        **fpPtVsMassRCExLJ;             //!pT vs Avg mass of RC excluding area around leading jet
  TProfile        **fpPtVsMassPerpConeLJ;         //!pT vs Avg mass of cone perpendicular to leading jet
  TProfile        **fpPtVsMassPerpConeTJ;         //!pT vs Avg mass of cone perpendicular to all tagged jets

  TH2F            **fh2EtaVsMassRC;               //!eta vs mass of RC
  TH2F            **fh2EtaVsMassRCExLJ;           //!eta vs mass of RC excluding area around leading jet
  TH2F            **fh2EtaVsMassPerpConeLJ;       //!eta vs mass of cone perpendicular to leading jet
  TH2F            **fh2EtaVsMassPerpConeTJ;       //!eta vs mass of cone perpendicular to all tagged jets

  TH2F             *fh2CentVsMassRC;              //!cent vs mass of RC
  TH2F             *fh2CentVsMassRCExLJ;          //!cent vs mass of RC excluding area around leading jet
  TH2F             *fh2CentVsMassPerpConeLJ;      //!cent vs mass of RC perpendicular to leading jet
  TH2F             *fh2CentVsMassPerpConeTJ;      //!cent vs mass of RC perpendicular to tagged jets

  TH2F             *fh2MultVsMassRC;              //!track multiplicity vs mass of RC
  TH2F             *fh2MultVsMassRCExLJ;          //!track multiplicity vs mass of RC excluding area around leading jet
  TH2F             *fh2MultVsMassPerpConeLJ;      //!track multiplicity vs mass of RC perpendicular to leading jet
  TH2F             *fh2MultVsMassPerpConeTJ;      //!track multiplicity vs mass of RC perpendicular to tagged jets

  TH2F             *fh2CentVsMedianMassRC;        //!cent vs median mass of all RC in event
  TH2F             *fh2CentVsMedianMassRCExLJ;    //!cent vs meidan mass mass of all RC in event excluding area around leading jet

  TH2F             *fh2MultVsMedianMassRC;        //!cent vs median mass of all RC in event
  TH2F             *fh2MultVsMedianMassRCExLJ;    //!cent vs meidan mass mass of all RC in event excluding area around leading jet

  TH2F             *fh2CentVsMeanMassRC;          //!cent vs median mass of all RC in event
  TH2F             *fh2CentVsMeanMassRCExLJ;      //!cent vs meidan mass mass of all RC in event excluding area around leading jet

  TH2F             *fh2MultVsMeanMassRC;          //!cent vs median mass of all RC in event
  TH2F             *fh2MultVsMeanMassRCExLJ;      //!cent vs meidan mass mass of all RC in event excluding area around leading jet

  TH2F             *fh2CentVsMedianMassPerAreaRC;     //!cent vs median mass of all RC in event
  TH2F             *fh2CentVsMedianMassPerAreaRCExLJ; //!cent vs meidan mass mass of all RC in event excluding area around leading jet

  TH2F             *fh2MultVsMedianMassPerAreaRC;     //!cent vs median mass of all RC in event
  TH2F             *fh2MultVsMedianMassPerAreaRCExLJ; //!cent vs meidan mass mass of all RC in event excluding area around leading jet

  AliAnalysisTaskEmcalJetMassBkg(const AliAnalysisTaskEmcalJetMassBkg&);            // not implemented
  AliAnalysisTaskEmcalJetMassBkg &operator=(const AliAnalysisTaskEmcalJetMassBkg&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetMassBkg, 4)
};
#endif

