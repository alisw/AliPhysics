#ifndef ALIANALYSISTASKJETMASSRESPONSEDET_H
#define ALIANALYSISTASKJETMASSRESPONSEDET_H

class TH1;
class TH2;
class TH3;
class TH3F;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAnalysisManager;
class AliJetContainer;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskJetMassResponseDet : public AliAnalysisTaskEmcalJet {
 public:
  enum JetMassType {
    kRaw   = 0,  //mass form anti-kt 4-vector
    kDeriv = 1   //area based subtracted jet mass
  };

  AliAnalysisTaskJetMassResponseDet();
  AliAnalysisTaskJetMassResponseDet(const char *name);
  virtual ~AliAnalysisTaskJetMassResponseDet();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainerPart(Int_t c)                             { fContainerPart     = c   ; }
  void SetJetContainerDet(Int_t c)                              { fContainerDet      = c   ; }
  void SetJetMassType(JetMassType t)                            { fJetMassType       = t   ; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Double_t                            GetJetMass(AliEmcalJet *jet);
 
  Int_t                               fContainerPart;              // particle level jets
  Int_t                               fContainerDet;               // detector level jets
  JetMassType                         fJetMassType;                // jet mass type to be used

  TH2F            *fh2PtVsMassJetPartAll;            //!pT vs mass of all particle level jets
  TH2F            *fh2PtVsMassJetPartMatch;          //!pT vs mass of all particle level jets matched to a detector level jet
  TH2F            *fh2PtVsMassJetPartTagged;         //!pT vs mass of tagged particle level jets
  TH2F            *fh2PtVsMassJetPartTaggedMatch;    //!pT vs mass of tagged particle level jets matched to a detector level jet
  TH2F            *fh2PtVsMassJetDetAll;             //!pT vs mass of all detector level jets
  TH2F            *fh2PtVsMassJetDetTagged;          //!pT vs mass of tagged detector level jets
  TH2F            *fh2EtaPhiMatchedDet;              //!eta,phi of matched detector level jets
  TH2F            *fh2EtaPhiMatchedPart;             //!eta,phi of matched particle level jets
  THnSparse       *fhnMassResponse;                  //!response matrix

  TH1D            *fh1AreaPartAll;                   //!area of all particle level jets
  TH1D            *fh1AreaDetAll;                    //!area of all detector level jets

 private:
  AliAnalysisTaskJetMassResponseDet(const AliAnalysisTaskJetMassResponseDet&);            // not implemented
  AliAnalysisTaskJetMassResponseDet &operator=(const AliAnalysisTaskJetMassResponseDet&); // not implemented

  ClassDef(AliAnalysisTaskJetMassResponseDet, 3)
};
#endif

