#ifndef ALIITSDEDXANALYZER_H
#define ALIITSDEDXANALYZER_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Class for:                                                    //
// - building histograms of dE/dx in bins of pt/layer/particle   //
// - comparing dEdx signal with Bethe-Bloch expected value       //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TParticle.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>

class AliESDEvent;
class AliStack;
class AliLog;

class AliITSdEdxAnalyzer : public TObject {

 public:
  AliITSdEdxAnalyzer();
  AliITSdEdxAnalyzer(const Int_t npbins, const Float_t pmin, const Float_t pmax);
  virtual ~AliITSdEdxAnalyzer();

  void SetMaterial(Float_t dens, Float_t thick){
    fDensity=dens;
    fThickness=thick;
  }
  void SetMomentumBins(const Int_t npbins, const Float_t pmin, const Float_t pmax);
  void SetUseBBFromAliExternalTrackParam() {fBBmodel=0;}
  void SetUseBBFromAliITSpidESD() {fBBmodel=1;}
  void SetMIPdEdx(Float_t val){fMIP=val;}
  void SetTPCMinimumPIDProb(Float_t min){fTPCpidCut=min;}

  void ReadEvent(const AliESDEvent* ev, AliStack* stack=0);


  Double_t BetheBloch(TParticle* part) const{
    return BetheBloch(part->P(),part->GetMass());
  }
  Double_t BetheBloch(const Float_t p, const Float_t m) const;


  TH1F* GetSingleLayerdEdxHisto(const Int_t lay, const Int_t pdgcode, const Int_t pbin) const {
    if(lay<3 || lay>6) AliFatal("Wrong LayerNumber");
    return fHistodEdx[GetIndex(GetSpecieBin(pdgcode),pbin,lay-3)];
  }
  TH1F* GetTruncatedMeandEdxHisto(const Int_t pdgcode, const Int_t pbin) const {
    return fHistodEdx[GetIndex(GetSpecieBin(pdgcode),pbin,4)];
  }

  TH1F* GetSingleLayerDeltadEdxHisto(const Int_t lay, const Int_t pdgcode, const Int_t pbin) const {
    if(lay<3 || lay>6) AliFatal("Wrong LayerNumber");
    return fHistoDeltadEdx[GetIndex(GetSpecieBin(pdgcode),pbin,lay-3)];
  }
  TH1F* GetTruncatedMeanDeltadEdxHisto(const Int_t pdgcode, const Int_t pbin) const {
    return fHistoDeltadEdx[GetIndex(GetSpecieBin(pdgcode),pbin,4)];
  }

  TH2F* GetSingleLayerdEdxVsPHisto(const Int_t lay, const Int_t pdgcode) const {
    if(lay<3 || lay>6) AliFatal("Wrong LayerNumber");
    return fHistodEdxVsP[GetIndex2(GetSpecieBin(pdgcode),lay-3)];
  }
  TH2F* GetTruncatedMeandEdxVsPHisto(const Int_t pdgcode) const {
    return fHistodEdxVsP[GetIndex2(GetSpecieBin(pdgcode),4)];
  }

  TGraph* GetBetheBlochGraph(const Int_t pdgcode) const;

  void WriteHistos(TString filename="ITS.dEdx.root") const;

 protected:

  Int_t GetNHistos() const {return kNParticles*kNValuesdEdx*fNPBins;}
  Int_t GetNHistos2() const {return kNParticles*kNValuesdEdx;}
  Int_t GetIndex(const Int_t iSp, const Int_t iP, const Int_t iVal) const {
    return iVal+kNValuesdEdx*iP+fNPBins*kNValuesdEdx*iSp;
  }
  Int_t GetIndex2(const Int_t iSp, const Int_t iVal) const {
    return iVal+kNValuesdEdx*iSp;
  }
  Int_t GetMomentumBin(Float_t p) const{
    Int_t iBin=(Int_t)((p-fPMin)/(fPMax-fPMin)*fNPBins);
    if(iBin>=0 && iBin<fNPBins)return iBin; 
    return -1;
  }
  Int_t GetSpecieBin(Int_t absPdgCode) const {
    for(Int_t iS=0; iS<kNParticles; iS++) if(absPdgCode==fgkPdgCode[iS]) return iS;
    return -1;
  }

  Int_t GetPaticleIdFromTPC(const AliESDtrack* track) const;
  void BookHistos();
  void DeleteHistos();
  
 private:

  AliITSdEdxAnalyzer(const AliITSdEdxAnalyzer& dum);
  AliITSdEdxAnalyzer& operator=(const AliITSdEdxAnalyzer& dum);

  enum {kNParticles=3};  
  enum {kNValuesdEdx=5};

  static const Int_t fgkPdgCode[kNParticles]; // initialized in the cxx
  static const Int_t fgkLayerCode[kNValuesdEdx]; // initialized in the cxx
  

  Int_t   fNPBins;         // Number of Momentum bins
  Float_t fPMin;           // Minimum P
  Float_t fPMax;           // Maximum P
  TH1F**  fHistodEdx;      // Array of histograms of dEdx_meas in bins of pt
  TH1F**  fHistoDeltadEdx; // Array of histograms of dEdx_meas-dEdx_expected in bins of pt
  TH2F**  fHistodEdxVsP;   // Array of histograms of dEdx_meas vs. pt
  Float_t fThickness;      // detector thickness (cm)
  Float_t fDensity;        // detector density (g/cm3)
  Int_t   fBBmodel;        // is there MC truth info?
  Float_t fMIP;            // normalization for MIP
  Float_t fTPCpidCut;      // minimum probability for PID in TPC   

  ClassDef(AliITSdEdxAnalyzer,0);
};
#endif
