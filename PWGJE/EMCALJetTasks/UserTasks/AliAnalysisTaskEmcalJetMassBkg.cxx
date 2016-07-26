//
// Jet mass background analysis task.
//
// Author: M.Verweij

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TChain.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <TRandom3.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliEmcalParticle.h"
#include "AliMCEvent.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliJetContainer.h"
#include "AliClusterContainer.h"
#include "AliParticleContainer.h"

#include "AliAnalysisTaskEmcalJetMassBkg.h"

ClassImp(AliAnalysisTaskEmcalJetMassBkg)

//________________________________________________________________________
AliAnalysisTaskEmcalJetMassBkg::AliAnalysisTaskEmcalJetMassBkg() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetMassBkg", kTRUE),
  fContainerBase(0),
  fMinRC2LJ(-1),
  fRCperEvent(10),
  fConeRadius(0.2),
  fConeMinEta(-0.9),
  fConeMaxEta(0.9),
  fConeMinPhi(0),
  fConeMaxPhi(TMath::Pi()*2),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0),
  fh2PtVsMassRC(0),
  fh2PtVsMassRCExLJDPhi(0),
  fh2PtVsMassPerpConeLJ(0),
  fh2PtVsMassPerpConeTJ(0),
  fh2PtVsERC(0),
  fh2PtVsERCExLJDPhi(0),
  fh2PtVsEPerpConeLJ(0),
  fh2PtVsEPerpConeTJ(0),
  fpPtVsMassRC(0),
  fpPtVsMassRCExLJ(0),
  fpPtVsMassPerpConeLJ(0),
  fpPtVsMassPerpConeTJ(0),
  fh2EtaVsMassRC(0),
  fh2EtaVsMassRCExLJ(0),
  fh2EtaVsMassPerpConeLJ(0),
  fh2EtaVsMassPerpConeTJ(0),
  fh2CentVsMassRC(0),
  fh2CentVsMassRCExLJ(0),
  fh2CentVsMassPerpConeLJ(0),
  fh2CentVsMassPerpConeTJ(0),
  fh2MultVsMassRC(0),
  fh2MultVsMassRCExLJ(0),
  fh2MultVsMassPerpConeLJ(0),
  fh2MultVsMassPerpConeTJ(0),
  fh2CentVsMedianMassRC(0),
  fh2CentVsMedianMassRCExLJ(0),
  fh2MultVsMedianMassRC(0),
  fh2MultVsMedianMassRCExLJ(0),
  fh2CentVsMeanMassRC(0),
  fh2CentVsMeanMassRCExLJ(0),
  fh2MultVsMeanMassRC(0),
  fh2MultVsMeanMassRCExLJ(0),
  fh2CentVsMedianMassPerAreaRC(0),
  fh2CentVsMedianMassPerAreaRCExLJ(0),
  fh2MultVsMedianMassPerAreaRC(0),
  fh2MultVsMedianMassPerAreaRCExLJ(0)
{
  // Default constructor.

  fh2PtVsMassRC            = new TH2F*[fNcentBins];
  fh2PtVsMassRCExLJDPhi    = new TH3F*[fNcentBins];
  fh2PtVsMassPerpConeLJ    = new TH2F*[fNcentBins];
  fh2PtVsMassPerpConeTJ    = new TH2F*[fNcentBins];

  fh2PtVsERC               = new TH2F*[fNcentBins];
  fh2PtVsERCExLJDPhi       = new TH3F*[fNcentBins];
  fh2PtVsEPerpConeLJ       = new TH2F*[fNcentBins];
  fh2PtVsEPerpConeTJ       = new TH2F*[fNcentBins];

  fpPtVsMassRC             = new TProfile*[fNcentBins];
  fpPtVsMassRCExLJ         = new TProfile*[fNcentBins];
  fpPtVsMassPerpConeLJ     = new TProfile*[fNcentBins];
  fpPtVsMassPerpConeTJ     = new TProfile*[fNcentBins];

  fh2EtaVsMassRC            = new TH2F*[fNcentBins];
  fh2EtaVsMassRCExLJ        = new TH2F*[fNcentBins];
  fh2EtaVsMassPerpConeLJ    = new TH2F*[fNcentBins];
  fh2EtaVsMassPerpConeTJ    = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2PtVsMassRC[i]             = 0;
    fh2PtVsMassRCExLJDPhi[i]     = 0;
    fh2PtVsMassPerpConeLJ[i]     = 0;
    fh2PtVsMassPerpConeTJ[i]     = 0;

    fh2PtVsERC[i]                = 0;
    fh2PtVsERCExLJDPhi[i]        = 0;
    fh2PtVsEPerpConeLJ[i]        = 0;
    fh2PtVsEPerpConeTJ[i]        = 0;

    fpPtVsMassRC[i]              = 0;
    fpPtVsMassRCExLJ[i]          = 0;
    fpPtVsMassPerpConeLJ[i]      = 0;
    fpPtVsMassPerpConeTJ[i]      = 0;

    fh2EtaVsMassRC[i]             = 0;
    fh2EtaVsMassRCExLJ[i]         = 0;
    fh2EtaVsMassPerpConeLJ[i]     = 0;
    fh2EtaVsMassPerpConeTJ[i]     = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
  
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMassBkg::AliAnalysisTaskEmcalJetMassBkg(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),  
  fContainerBase(0),
  fMinRC2LJ(-1),
  fRCperEvent(10),
  fConeRadius(0.2),
  fConeMinEta(-0.9),
  fConeMaxEta(0.9),
  fConeMinPhi(0),
  fConeMaxPhi(TMath::Pi()*2),
  fJetsCont(0),
  fTracksCont(0),
  fCaloClustersCont(0),
  fh2PtVsMassRC(0),
  fh2PtVsMassRCExLJDPhi(0),
  fh2PtVsMassPerpConeLJ(0),
  fh2PtVsMassPerpConeTJ(0),
  fh2PtVsERC(0),
  fh2PtVsERCExLJDPhi(0),
  fh2PtVsEPerpConeLJ(0),
  fh2PtVsEPerpConeTJ(0),
  fpPtVsMassRC(0),
  fpPtVsMassRCExLJ(0),
  fpPtVsMassPerpConeLJ(0),
  fpPtVsMassPerpConeTJ(0),
  fh2EtaVsMassRC(0),
  fh2EtaVsMassRCExLJ(0),
  fh2EtaVsMassPerpConeLJ(0),
  fh2EtaVsMassPerpConeTJ(0),
  fh2CentVsMassRC(0),
  fh2CentVsMassRCExLJ(0),
  fh2CentVsMassPerpConeLJ(0),
  fh2CentVsMassPerpConeTJ(0),
  fh2MultVsMassRC(0),
  fh2MultVsMassRCExLJ(0),
  fh2MultVsMassPerpConeLJ(0),
  fh2MultVsMassPerpConeTJ(0),
  fh2CentVsMedianMassRC(0),
  fh2CentVsMedianMassRCExLJ(0),
  fh2MultVsMedianMassRC(0),
  fh2MultVsMedianMassRCExLJ(0),
  fh2CentVsMeanMassRC(0),
  fh2CentVsMeanMassRCExLJ(0),
  fh2MultVsMeanMassRC(0),
  fh2MultVsMeanMassRCExLJ(0),
  fh2CentVsMedianMassPerAreaRC(0),
  fh2CentVsMedianMassPerAreaRCExLJ(0),
  fh2MultVsMedianMassPerAreaRC(0),
  fh2MultVsMedianMassPerAreaRCExLJ(0)
{
  // Standard constructor.

  fh2PtVsMassRC            = new TH2F*[fNcentBins];
  fh2PtVsMassRCExLJDPhi    = new TH3F*[fNcentBins];
  fh2PtVsMassPerpConeLJ    = new TH2F*[fNcentBins];
  fh2PtVsMassPerpConeTJ    = new TH2F*[fNcentBins];

  fh2PtVsERC               = new TH2F*[fNcentBins];
  fh2PtVsERCExLJDPhi       = new TH3F*[fNcentBins];
  fh2PtVsEPerpConeLJ       = new TH2F*[fNcentBins];
  fh2PtVsEPerpConeTJ       = new TH2F*[fNcentBins];

  fpPtVsMassRC             = new TProfile*[fNcentBins];
  fpPtVsMassRCExLJ         = new TProfile*[fNcentBins];
  fpPtVsMassPerpConeLJ     = new TProfile*[fNcentBins];
  fpPtVsMassPerpConeTJ     = new TProfile*[fNcentBins];

  fh2EtaVsMassRC            = new TH2F*[fNcentBins];
  fh2EtaVsMassRCExLJ        = new TH2F*[fNcentBins];
  fh2EtaVsMassPerpConeLJ    = new TH2F*[fNcentBins];
  fh2EtaVsMassPerpConeTJ    = new TH2F*[fNcentBins];

  for (Int_t i = 0; i < fNcentBins; i++) {
    fh2PtVsMassRC[i]             = 0;
    fh2PtVsMassRCExLJDPhi[i]     = 0;
    fh2PtVsMassPerpConeLJ[i]     = 0;
    fh2PtVsMassPerpConeTJ[i]     = 0;

    fh2PtVsERC[i]                = 0;
    fh2PtVsERCExLJDPhi[i]        = 0;
    fh2PtVsEPerpConeLJ[i]        = 0;
    fh2PtVsEPerpConeTJ[i]        = 0;

    fpPtVsMassRC[i]              = 0;
    fpPtVsMassRCExLJ[i]          = 0;
    fpPtVsMassPerpConeLJ[i]      = 0;
    fpPtVsMassPerpConeTJ[i]      = 0;

    fh2EtaVsMassRC[i]             = 0;
    fh2EtaVsMassRCExLJ[i]         = 0;
    fh2EtaVsMassPerpConeLJ[i]     = 0;
    fh2EtaVsMassPerpConeTJ[i]     = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetMassBkg::~AliAnalysisTaskEmcalJetMassBkg()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetMassBkg::UserCreateOutputObjects()
{
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fJetsCont         = GetJetContainer(fContainerBase);
  fTracksCont       = fJetsCont->GetParticleContainer();
  fCaloClustersCont = fJetsCont->GetClusterContainer();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  const Int_t nBinsPt  = 250;
  const Double_t minPt = -50.;
  const Double_t maxPt = 200.;

  const Int_t nBinsE  = 250;
  const Double_t minE = -50.;
  const Double_t maxE = 200.;

  const Int_t nBinsM  = 150;
  const Double_t minM = -50.;
  const Double_t maxM = 100.;

  const Int_t nBinsEta  = 100;
  const Double_t minEta = -1.;
  const Double_t maxEta =  1.;

  const Int_t nBinsCent  = 100;
  const Double_t minCent = 0.;
  const Double_t maxCent = 100.;

  const Int_t nBinsMult  = 400;
  const Double_t minMult = 0.;
  const Double_t maxMult = 4000.;

  fh2CentVsMassRC = new TH2F("fh2CentVsMassRC","fh2CentVsMassRC;cent;#it{M}_{RC}",nBinsCent,minCent,maxCent,nBinsM,minM,maxM);
  fOutput->Add(fh2CentVsMassRC);  

  fh2CentVsMassRCExLJ = new TH2F("fh2CentVsMassRCExLJ","fh2CentVsMassRCExLJ;cent;#it{M}_{RC}",nBinsCent,minCent,maxCent,nBinsM,minM,maxM);
  fOutput->Add(fh2CentVsMassRCExLJ);  

  fh2CentVsMassPerpConeLJ = new TH2F("fh2CentVsMassPerpConeLJ","fh2CentVsMassPerpConeLJ;cent;#it{M}_{RC}",nBinsCent,minCent,maxCent,nBinsM,minM,maxM);
  fOutput->Add(fh2CentVsMassPerpConeLJ);  

  fh2CentVsMassPerpConeTJ = new TH2F("fh2CentVsMassPerpConeTJ","fh2CentVsMassPerpConeTJ;cent;#it{M}_{RC}",nBinsCent,minCent,maxCent,nBinsM,minM,maxM);
  fOutput->Add(fh2CentVsMassPerpConeTJ); 

  fh2MultVsMassRC = new TH2F("fh2MultVsMassRC","fh2MultVsMassRC;#it{N}_{track};#it{M}_{RC}",nBinsMult,minMult,maxMult,nBinsM,minM,maxM);
  fOutput->Add(fh2MultVsMassRC);  

  fh2MultVsMassRCExLJ = new TH2F("fh2MultVsMassRCExLJ","fh2MultVsMassRCExLJ;#it{N}_{track};#it{M}_{RC}",nBinsMult,minMult,maxMult,nBinsM,minM,maxM);
  fOutput->Add(fh2MultVsMassRCExLJ);  

  fh2MultVsMassPerpConeLJ = new TH2F("fh2MultVsMassPerpConeLJ","fh2MultVsMassPerpConeLJ;#it{N}_{track};#it{M}_{RC}",nBinsMult,minMult,maxMult,nBinsM,minM,maxM);
  fOutput->Add(fh2MultVsMassPerpConeLJ);  

  fh2MultVsMassPerpConeTJ = new TH2F("fh2MultVsMassPerpConeTJ","fh2MultVsMassPerpConeTJ;#it{N}_{track};#it{M}_{RC}",nBinsMult,minMult,maxMult,nBinsM,minM,maxM);
  fOutput->Add(fh2MultVsMassPerpConeTJ); 

  fh2CentVsMedianMassRC = new TH2F("fh2CentVsMedianMassRC","fh2CentVsMedianMassRC;cent;#it{M}_{RC}",nBinsCent,minCent,maxCent,nBinsM,minM,maxM);
  fOutput->Add(fh2CentVsMedianMassRC);  

  fh2CentVsMedianMassRCExLJ = new TH2F("fh2CentVsMedianMassRCExLJ","fh2CentVsMedianMassRCExLJ;cent;#it{M}_{RC}",nBinsCent,minCent,maxCent,nBinsM,minM,maxM);
  fOutput->Add(fh2CentVsMedianMassRCExLJ);  

  fh2MultVsMedianMassRC = new TH2F("fh2MultVsMedianMassRC","fh2MultVsMedianMassRC;#it{N}_{track};#it{M}_{RC}",nBinsMult,minMult,maxMult,nBinsM,minM,maxM);
  fOutput->Add(fh2MultVsMedianMassRC);  

  fh2MultVsMedianMassRCExLJ = new TH2F("fh2MultVsMedianMassRCExLJ","fh2MultVsMedianMassRCExLJ;#it{N}_{track};#it{M}_{RC}",nBinsMult,minMult,maxMult,nBinsM,minM,maxM);
  fOutput->Add(fh2MultVsMedianMassRCExLJ);  

  fh2CentVsMeanMassRC = new TH2F("fh2CentVsMeanMassRC","fh2CentVsMeanMassRC;cent;#it{M}_{RC}",nBinsCent,minCent,maxCent,nBinsM,minM,maxM);
  fOutput->Add(fh2CentVsMeanMassRC);  

  fh2CentVsMeanMassRCExLJ = new TH2F("fh2CentVsMeanMassRCExLJ","fh2CentVsMeanMassRCExLJ;cent;#it{M}_{RC}",nBinsCent,minCent,maxCent,nBinsM,minM,maxM);
  fOutput->Add(fh2CentVsMeanMassRCExLJ);  

  fh2MultVsMeanMassRC = new TH2F("fh2MultVsMeanMassRC","fh2MultVsMeanMassRC;#it{N}_{track};#it{M}_{RC}",nBinsMult,minMult,maxMult,nBinsM,minM,maxM);
  fOutput->Add(fh2MultVsMeanMassRC);  

  fh2MultVsMeanMassRCExLJ = new TH2F("fh2MultVsMeanMassRCExLJ","fh2MultVsMeanMassRCExLJ;#it{N}_{track};#it{M}_{RC}",nBinsMult,minMult,maxMult,nBinsM,minM,maxM);
  fOutput->Add(fh2MultVsMeanMassRCExLJ);  

  fh2CentVsMedianMassPerAreaRC = new TH2F("fh2CentVsMedianMassPerAreaRC","fh2CentVsMedianMassPerAreaRC;cent;#it{M}_{RC}/A",nBinsCent,minCent,maxCent,nBinsM,minM,maxM);
  fOutput->Add(fh2CentVsMedianMassPerAreaRC);  

  fh2CentVsMedianMassPerAreaRCExLJ = new TH2F("fh2CentVsMedianMassPerAreaRCExLJ","fh2CentVsMedianMassPerAreaRCExLJ;cent;#it{M}_{RC}/A",nBinsCent,minCent,maxCent,nBinsM,minM,maxM);
  fOutput->Add(fh2CentVsMedianMassPerAreaRCExLJ);  

  fh2MultVsMedianMassPerAreaRC = new TH2F("fh2MultVsMedianMassPerAreaRC","fh2MultVsMedianMassPerAreaRC;#it{N}_{track};#it{M}_{RC}/A",nBinsMult,minMult,maxMult,nBinsM,minM,maxM);
  fOutput->Add(fh2MultVsMedianMassPerAreaRC);  

  fh2MultVsMedianMassPerAreaRCExLJ = new TH2F("fh2MultVsMedianMassPerAreaRCExLJ","fh2MultVsMedianMassPerAreaRCExLJ;#it{N}_{track};#it{M}_{RC}/A",nBinsMult,minMult,maxMult,nBinsM,minM,maxM);
  fOutput->Add(fh2MultVsMedianMassPerAreaRCExLJ);  

  TString histName = "";
  TString histTitle = "";
  for (Int_t i = 0; i < fNcentBins; i++) {
    histName = TString::Format("fh2PtVsMassRC_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,RC};#it{M}_{RC}",histName.Data());
    fh2PtVsMassRC[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh2PtVsMassRC[i]);

    histName = TString::Format("fh2PtVsMassRCExLJDPhi_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,RC};#it{M}_{RC}",histName.Data());
    fh2PtVsMassRCExLJDPhi[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,72,-0.5*TMath::Pi(),1.5*TMath::Pi());
    fOutput->Add(fh2PtVsMassRCExLJDPhi[i]);

    histName = TString::Format("fh2PtVsMassPerpConeLJ_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,PerpConeLJ};#it{M}_{PerpConeLJ}",histName.Data());
    fh2PtVsMassPerpConeLJ[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh2PtVsMassPerpConeLJ[i]);

    histName = TString::Format("fh2PtVsMassPerpConeTJ_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,PerpConeTJ};#it{M}_{PerpConeTJ}",histName.Data());
    fh2PtVsMassPerpConeTJ[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsM,minM,maxM);
    fOutput->Add(fh2PtVsMassPerpConeTJ[i]);
    //
    histName = TString::Format("fh2PtVsERC_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,RC};#it{M}_{RC}",histName.Data());
    fh2PtVsERC[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsE,minE,maxE);
    fOutput->Add(fh2PtVsERC[i]);

    histName = TString::Format("fh2PtVsERCExLJDPhi_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,RC};#it{M}_{RC}",histName.Data());
    fh2PtVsERCExLJDPhi[i] = new TH3F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsPt,minPt,maxPt,72,-0.5*TMath::Pi(),1.5*TMath::Pi());
    fOutput->Add(fh2PtVsERCExLJDPhi[i]);

    histName = TString::Format("fh2PtVsEPerpConeLJ_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,PerpConeLJ};#it{M}_{PerpConeLJ}",histName.Data());
    fh2PtVsEPerpConeLJ[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsE,minE,maxE);
    fOutput->Add(fh2PtVsEPerpConeLJ[i]);

    histName = TString::Format("fh2PtVsEPerpConeTJ_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,PerpConeTJ};#it{M}_{PerpConeTJ}",histName.Data());
    fh2PtVsEPerpConeTJ[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt,nBinsE,minE,maxE);
    fOutput->Add(fh2PtVsEPerpConeTJ[i]);
    //
    histName = TString::Format("fh2EtaVsMassRC_%d",i);
    histTitle = TString::Format("%s;#eta_{RC};#it{M}_{RC}",histName.Data());
    fh2EtaVsMassRC[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsEta,minEta,maxEta,nBinsM,minM,maxM);
    fOutput->Add(fh2EtaVsMassRC[i]);

    histName = TString::Format("fh2EtaVsMassRCExLJ_%d",i);
    histTitle = TString::Format("%s;#eta_{RC};#it{M}_{RC}",histName.Data());
    fh2EtaVsMassRCExLJ[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsEta,minEta,maxEta,nBinsM,minM,maxM);
    fOutput->Add(fh2EtaVsMassRCExLJ[i]);

    histName = TString::Format("fh2EtaVsMassPerpConeLJ_%d",i);
    histTitle = TString::Format("%s;#eta_{PerpConeLJ};#it{M}_{PerpConeLJ}",histName.Data());
    fh2EtaVsMassPerpConeLJ[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsEta,minEta,maxEta,nBinsM,minM,maxM);
    fOutput->Add(fh2EtaVsMassPerpConeLJ[i]);

    histName = TString::Format("fh2EtaVsMassPerpConeTJ_%d",i);
    histTitle = TString::Format("%s;#eta_{PerpConeTJ};#it{M}_{PerpConeTJ}",histName.Data());
    fh2EtaVsMassPerpConeTJ[i] = new TH2F(histName.Data(),histTitle.Data(),nBinsEta,minEta,maxEta,nBinsM,minM,maxM);
    fOutput->Add(fh2EtaVsMassPerpConeTJ[i]);

    histName = TString::Format("fpPtVsMassRC_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,RC};Avg #it{M}_{RC}",histName.Data());
    fpPtVsMassRC[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fpPtVsMassRC[i]);

    histName = TString::Format("fpPtVsMassRCExLJ_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,RC};Avg #it{M}_{RC}",histName.Data());
    fpPtVsMassRCExLJ[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fpPtVsMassRCExLJ[i]);

    histName = TString::Format("fpPtVsMassPerpConeLJ_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,RC};Avg #it{M}_{RC}",histName.Data());
    fpPtVsMassPerpConeLJ[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fpPtVsMassPerpConeLJ[i]);

    histName = TString::Format("fpPtVsMassPerpConeTJ_%d",i);
    histTitle = TString::Format("%s;#it{p}_{T,RC};Avg #it{M}_{RC}",histName.Data());
    fpPtVsMassPerpConeTJ[i] = new TProfile(histName.Data(),histTitle.Data(),nBinsPt,minPt,maxPt);
    fOutput->Add(fpPtVsMassPerpConeTJ[i]);
  }

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fOutput->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
    if(hn)hn->Sumw2();
  }

  TH1::AddDirectory(oldStatus);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMassBkg::Run()
{
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMassBkg::FillHistograms()
{
  // Fill histograms.

  const Float_t rcArea = fConeRadius * fConeRadius * TMath::Pi();

  //Event properties
  Double_t rho = GetRhoVal(fContainerBase);
  Int_t trackMult = fTracksCont->GetNAcceptedParticles();

  //Leading jet
  AliEmcalJet* jet = NULL;
  if (fJetsCont)
    jet = fJetsCont->GetLeadingJet("rho");
  
  TLorentzVector lvRC(0.,0.,0.,0.);
  TLorentzVector lvRCS(0.,0.,0.,0.);
  Float_t RCpt = 0;
  Float_t RCeta = 0;
  Float_t RCphi = 0;

  static Double_t massvecRC[999];
  static Double_t massPerAreavecRC[999];

  static Double_t massvecRCExLJ[999];
  static Double_t massPerAreavecRCExLJ[999];

  Int_t nRCAcc = 0;
  Int_t nRCExLJAcc = 0;

  for (Int_t i = 0; i < fRCperEvent; i++) {
    // Simple random cones
    lvRC.SetPxPyPzE(0.,0.,0.,0.);
    RCpt = 0;
    RCeta = 0;
    RCphi = 0;
    GetRandomCone(lvRC,RCpt, RCeta, RCphi, fTracksCont, fCaloClustersCont, 0);
    if (RCpt > 0) {
      lvRCS = GetSubtractedVector(lvRC.Pt(),lvRC.Eta(),lvRC.Phi(),lvRC.E());
      fh2PtVsMassRC[fCentBin]->Fill(RCpt-rho*rcArea,lvRCS.M());
      fh2PtVsERC[fCentBin]->Fill(RCpt-rho*rcArea,lvRCS.E());
      fpPtVsMassRC[fCentBin]->Fill(RCpt-rho*rcArea,lvRCS.M());
      fh2EtaVsMassRC[fCentBin]->Fill(lvRCS.Eta(),lvRCS.M());
      fh2CentVsMassRC->Fill(fCent,lvRCS.M());
      fh2MultVsMassRC->Fill(trackMult,lvRCS.M());

      massvecRC[nRCAcc] = lvRCS.M();
      massPerAreavecRC[nRCAcc] = massvecRC[nRCAcc]/rcArea;
      ++nRCAcc;
    }

    if (fJetsCont && jet) {
      // Random cones far away from leading jet(s)
      lvRC.SetPxPyPzE(0.,0.,0.,0.);
      RCpt = 0;
      RCeta = 0;
      RCphi = 0;
      GetRandomCone(lvRC,RCpt, RCeta, RCphi, fTracksCont, fCaloClustersCont, jet);
      if (RCpt > 0 && jet) {
	lvRCS = GetSubtractedVector(lvRC.Pt(),lvRC.Eta(),lvRC.Phi(),lvRC.E());
	Float_t dphi = RCphi - jet->Phi();
	if (dphi > 1.5*TMath::Pi()) dphi -= TMath::Pi() * 2;
	if (dphi < -0.5*TMath::Pi()) dphi += TMath::Pi() * 2;
	fh2PtVsMassRCExLJDPhi[fCentBin]->Fill(RCpt-rho*rcArea,lvRCS.M(),dphi);
	fh2PtVsERCExLJDPhi[fCentBin]->Fill(RCpt-rho*rcArea,lvRCS.E(),dphi);
	fpPtVsMassRCExLJ[fCentBin]->Fill(RCpt-rho*rcArea,lvRCS.M());
	fh2EtaVsMassRCExLJ[fCentBin]->Fill(lvRCS.Eta(),lvRCS.M());
	fh2CentVsMassRCExLJ->Fill(fCent,lvRCS.M());
	fh2MultVsMassRCExLJ->Fill(trackMult,lvRCS.M());

	massvecRCExLJ[nRCExLJAcc] = lvRCS.M();
	massPerAreavecRCExLJ[nRCExLJAcc] = lvRCS.M()/rcArea;
	++nRCExLJAcc;
      }
    }
  }//RC loop

  Double_t medianRC, medianRCExLJ = 0.;
  medianRC = TMath::Median(nRCAcc,massvecRC);
  medianRCExLJ = TMath::Median(nRCExLJAcc,massvecRCExLJ);

  fh2CentVsMedianMassRC->Fill(fCent,medianRC);
  fh2CentVsMedianMassRCExLJ->Fill(fCent,medianRCExLJ);

  fh2MultVsMedianMassRC->Fill(trackMult,medianRC);
  fh2MultVsMedianMassRCExLJ->Fill(trackMult,medianRCExLJ);

  Double_t meanRC = 0.; Double_t meanRCExLJ = 0.;
  if(nRCAcc>0)   meanRC = TMath::Mean(nRCAcc,massvecRC);
  if(nRCExLJAcc) meanRCExLJ = TMath::Mean(nRCExLJAcc,massvecRCExLJ);

  fh2CentVsMeanMassRC->Fill(fCent,meanRC);
  fh2CentVsMeanMassRCExLJ->Fill(fCent,meanRCExLJ);

  fh2MultVsMeanMassRC->Fill(trackMult,meanRC);
  fh2MultVsMeanMassRCExLJ->Fill(trackMult,meanRCExLJ);

  Double_t medianPerAreaRC, medianPerAreaRCExLJ = 0.;
  medianPerAreaRC = TMath::Median(nRCAcc,massPerAreavecRC);
  medianPerAreaRCExLJ = TMath::Median(nRCExLJAcc,massPerAreavecRCExLJ);

  fh2CentVsMedianMassPerAreaRC->Fill(fCent,medianPerAreaRC);
  fh2CentVsMedianMassPerAreaRCExLJ->Fill(fCent,medianPerAreaRCExLJ);

  fh2MultVsMedianMassPerAreaRC->Fill(trackMult,medianPerAreaRC);
  fh2MultVsMedianMassPerAreaRCExLJ->Fill(trackMult,medianPerAreaRCExLJ);


  if(fJetsCont && jet) {
    //cone perpendicular to leading jet
    TLorentzVector lvPC(0.,0.,0.,0.);
    TLorentzVector lvPCS(0.,0.,0.,0.);
    Float_t PCpt = 0;
    Float_t PCeta = 0;
    Float_t PCphi = 0;
    if(jet) {
      GetPerpCone(lvPC,PCpt, PCeta, PCphi, fTracksCont, fCaloClustersCont, jet);
      if(PCpt>0.) {
	lvPCS = GetSubtractedVector(lvPC.Pt(),lvPC.Eta(),lvPC.Phi(),lvPC.E());
	fh2PtVsMassPerpConeLJ[fCentBin]->Fill(PCpt-rho*rcArea,lvPCS.M());
	fh2PtVsEPerpConeLJ[fCentBin]->Fill(PCpt-rho*rcArea,lvPCS.E());
	fpPtVsMassPerpConeLJ[fCentBin]->Fill(PCpt-rho*rcArea,lvPCS.M());
	fh2EtaVsMassPerpConeLJ[fCentBin]->Fill(lvPCS.Eta(),lvPCS.M());
	fh2CentVsMassPerpConeLJ->Fill(fCent,lvPCS.M());
	fh2MultVsMassPerpConeLJ->Fill(trackMult,lvPCS.M());
      }
    }
    //cone perpendicular to all tagged jets
    for(int i = 0; i < fJetsCont->GetNJets();++i) {
      jet = static_cast<AliEmcalJet*>(fJetsCont->GetAcceptJet(i));
      if(!jet) continue;

      if(jet->GetTagStatus()<1)
	continue;

      lvPC.SetPxPyPzE(0.,0.,0.,0.);
      PCpt = 0;
      PCeta = 0;
      PCphi = 0;
      GetPerpCone(lvPC,PCpt, PCeta, PCphi, fTracksCont, fCaloClustersCont, jet);
      if(PCpt>0.) {
	lvPCS = GetSubtractedVector(lvPC.Pt(),lvPC.Eta(),lvPC.Phi(),lvPC.E());
	fh2PtVsMassPerpConeTJ[fCentBin]->Fill(PCpt-rho*rcArea,lvPCS.M());
	fh2PtVsEPerpConeTJ[fCentBin]->Fill(PCpt-rho*rcArea,lvPCS.E());
	fpPtVsMassPerpConeTJ[fCentBin]->Fill(PCpt-rho*rcArea,lvPCS.M());
	fh2EtaVsMassPerpConeTJ[fCentBin]->Fill(lvPCS.Eta(),lvPCS.M());
	fh2CentVsMassPerpConeTJ->Fill(fCent,lvPCS.M());
	fh2MultVsMassPerpConeTJ->Fill(trackMult,lvPCS.M());
      }
    }//jet loop
  }

  return kTRUE;

}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetMassBkg::GetRandomCone(TLorentzVector& lvRC,Float_t &pt, Float_t &eta, Float_t &phi,
						   AliParticleContainer* tracks, AliClusterContainer* clusters,
						   AliEmcalJet *jet) const
{
  // Get rigid cone.
  lvRC.SetPxPyPzE(0.,0.,0.,0.);

  eta = -999;
  phi = -999;
  pt = 0;

  if (!tracks && !clusters)
    return;

  Float_t LJeta = 999;
  Float_t LJphi = 999;

  if (jet) {
    LJeta = jet->Eta();
    LJphi = jet->Phi();
  }

  Float_t maxEta = fConeMaxEta;
  Float_t minEta = fConeMinEta;
  Float_t maxPhi = fConeMaxPhi;
  Float_t minPhi = fConeMinPhi;

  if (maxPhi > TMath::Pi() * 2) maxPhi = TMath::Pi() * 2;
  if (minPhi < 0) minPhi = 0;
  
  Float_t dLJ = 0;
  Int_t repeats = 0;
  Bool_t reject = kTRUE;
  do {
    eta = gRandom->Rndm() * (maxEta - minEta) + minEta;
    phi = gRandom->Rndm() * (maxPhi - minPhi) + minPhi;
    dLJ = TMath::Sqrt((LJeta - eta) * (LJeta - eta) + (LJphi - phi) * (LJphi - phi));

    repeats++;
  } while (dLJ < fMinRC2LJ && repeats < 999 && reject);

  if (repeats == 999) {
    AliWarning(Form("%s: Could not get random cone!", GetName()));
    return;
  }

  GetCone(lvRC,pt,eta,phi,tracks,clusters);


}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetMassBkg::GetCone(TLorentzVector& lvRC,Float_t &pt, Float_t eta, Float_t phi, AliParticleContainer* tracks, AliClusterContainer* clusters) const
{

  pt = 0.;
  lvRC.SetPxPyPzE(0.,0.,0.,0.);

  if (clusters) {
    clusters->ResetCurrentID();
    AliVCluster* cluster = clusters->GetNextAcceptCluster();
    while (cluster) {     
      TLorentzVector nPart;
      cluster->GetMomentum(nPart, const_cast<Double_t*>(fVertex));

      Float_t cluseta = nPart.Eta();
      Float_t clusphi = nPart.Phi();
      
      if (TMath::Abs(clusphi - phi) > TMath::Abs(clusphi - phi + 2 * TMath::Pi()))
	clusphi += 2 * TMath::Pi();
      if (TMath::Abs(clusphi - phi) > TMath::Abs(clusphi - phi - 2 * TMath::Pi()))
	clusphi -= 2 * TMath::Pi();
     
      Float_t d = TMath::Sqrt((cluseta - eta) * (cluseta - eta) + (clusphi - phi) * (clusphi - phi));
      if (d <= fConeRadius) {
	pt += nPart.Pt();
	TLorentzVector lvcl(nPart.Px(),nPart.Py(),nPart.Pz(),nPart.E());
	lvRC += lvcl;
      }

      cluster = clusters->GetNextAcceptCluster();
    }
  }

  if (tracks) {
    tracks->ResetCurrentID();
    AliVParticle* track = tracks->GetNextAcceptParticle(); 
    while(track) { 
      Float_t tracketa = track->Eta();
      Float_t trackphi = track->Phi();
      
      if (TMath::Abs(trackphi - phi) > TMath::Abs(trackphi - phi + 2 * TMath::Pi()))
	trackphi += 2 * TMath::Pi();
      if (TMath::Abs(trackphi - phi) > TMath::Abs(trackphi - phi - 2 * TMath::Pi()))
	trackphi -= 2 * TMath::Pi();
      
      Float_t d = TMath::Sqrt((tracketa - eta) * (tracketa - eta) + (trackphi - phi) * (trackphi - phi));
      if (d <= fConeRadius) {
	pt += track->Pt();
	TLorentzVector lvtr(track->Px(),track->Py(),track->Pz(),track->E());
	lvRC += lvtr;
      }

      track = tracks->GetNextAcceptParticle(); 
    }
  }

}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetMassBkg::GetPerpCone(TLorentzVector& lvRC,Float_t &pt, Float_t &eta, Float_t &phi, AliParticleContainer* tracks, AliClusterContainer* clusters, AliEmcalJet *jet) const
{
  // Get rigid cone.
  lvRC.SetPxPyPzE(0.,0.,0.,0.);

  eta = -999;
  phi = -999;
  pt = 0;

  if (!tracks && !clusters)
    return;

  if(!jet)
    return;

  Float_t LJeta = jet->Eta();
  Float_t LJphi = jet->Phi();

  eta = LJeta;
  phi = LJphi + 0.5*TMath::Pi();
  if(phi>TMath::TwoPi()) phi-=TMath::TwoPi();

  GetCone(lvRC,pt,eta,phi,tracks,clusters);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetMassBkg::SetConeEtaPhiEMCAL()
{
  // Set default cuts for full cones

  SetConeEtaLimits(-0.7+fConeRadius,0.7-fConeRadius);
  SetConePhiLimits(1.405+fConeRadius,3.135-fConeRadius);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetMassBkg::SetConeEtaPhiTPC()
{
  // Set default cuts for charged cones

  SetConeEtaLimits(-0.9+fConeRadius, 0.9-fConeRadius);
  SetConePhiLimits(-10, 10);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetMassBkg::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;
  if (fCaloClustersCont && fCaloClustersCont->GetArray() == 0) fCaloClustersCont = 0;

  if (fRCperEvent < 0) {
    Double_t area = (fConeMaxEta - fConeMinEta) * (fConeMaxPhi - fConeMinPhi);
    Double_t rcArea = TMath::Pi() * fConeRadius * fConeRadius;
    fRCperEvent = TMath::FloorNint(area / rcArea - 0.5);
    if (fRCperEvent == 0)
      fRCperEvent = 1;
  }

  if (fMinRC2LJ < 0)
    fMinRC2LJ = fConeRadius * 1.5;

  const Float_t maxDist = TMath::Max(fConeMaxPhi - fConeMinPhi, fConeMaxEta - fConeMinEta) / 2;
  if (fMinRC2LJ > maxDist) {
    AliWarning(Form("The parameter fMinRC2LJ = %f is too large for the considered acceptance. "
                    "Will use fMinRC2LJ = %f", fMinRC2LJ, maxDist));
    fMinRC2LJ = maxDist;
  }

}

  //________________________________________________________________________
TLorentzVector AliAnalysisTaskEmcalJetMassBkg::GetBkgVector(Double_t eta, Double_t phi, AliJetContainer *cont) {
  //get background vector

  Double_t rho  = cont->GetRhoVal();
  Double_t rhom = cont->GetRhoMassVal();
  TLorentzVector vpB(0.,0.,0.,0.);
  Double_t aRC = TMath::Pi()*fConeRadius*fConeRadius;
  vpB.SetPxPyPzE(rho*TMath::Cos(phi)*aRC,rho*TMath::Sin(phi)*aRC,(rho+rhom)*TMath::SinH(eta)*aRC,(rho+rhom)*TMath::CosH(eta)*aRC);
  return vpB;
}

//________________________________________________________________________
TLorentzVector AliAnalysisTaskEmcalJetMassBkg::GetSubtractedVector(Double_t pt, Double_t eta, Double_t phi, Double_t e) {
  //get subtracted vector
  TLorentzVector vpS;
  AliJetContainer *jetCont = GetJetContainer(fContainerBase);
  TLorentzVector vpBkg = GetBkgVector(eta,phi,jetCont);
  vpS.SetPxPyPzE(pt*TMath::Cos(phi)-vpBkg.Px(),pt*TMath::Sin(phi)-vpBkg.Py(),pt*TMath::SinH(eta)-vpBkg.Pz(),e-vpBkg.E());
  return vpS;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetMassBkg::RetrieveEventObjects() {
  //
  // retrieve event objects
  //

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;

  return kTRUE;

}

//_______________________________________________________________________
void AliAnalysisTaskEmcalJetMassBkg::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

