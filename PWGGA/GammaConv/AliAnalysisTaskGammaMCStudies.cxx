/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                       *
 * Author: Baldo Sahlmueller, Friederike Bock                     *
 * Version 1.0                                 *
 *                                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims    *
 * about the suitability of this software for any purpose. It is      *
 * provided "as is" without express or implied warranty.               *
 **************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis on conversion photons + calo photons
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPDGCode.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAnalysisTaskGammaMCStudies.h"
#include "AliVParticle.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include <algorithm>
#include <array>
#include <vector>
#include <map>

ClassImp(AliAnalysisTaskGammaMCStudies)

//________________________________________________________________________
AliAnalysisTaskGammaMCStudies::AliAnalysisTaskGammaMCStudies(): AliAnalysisTaskSE(),
  fOutputContainer(nullptr),
  fHistNEvents(nullptr),
  fHistXSection(nullptr),
  fHistPtHard(nullptr),
  fHistPtYPi0(nullptr),
  fHistPtYAllPhotons(nullptr),
  fHistPtYDecayPhotons(nullptr),
  fHistPtYDirectPhotons(nullptr),
  fHistPtYComptonPhotons(nullptr),
  fHistPtYAnnihilationPhotons(nullptr),
  fHistPtYgg2ggammaPhotons(nullptr),
  fHistPtYffbar2gammagammaPhotons(nullptr),
  fHistPtYgg2gammagammaPhotons(nullptr),
  fHistPtYFragmentationPhotons(nullptr),

  fHistIsoYPi0(nullptr),
  fHistIsoYAllPhotons(nullptr),
  fHistIsoYDecayPhotons(nullptr),
  fHistIsoYDirectPhotons(nullptr),
  fHistIsoYComptonPhotons(nullptr),
  fHistIsoYAnnihilationPhotons(nullptr),
  fHistIsoYgg2ggammaPhotons(nullptr),
  fHistIsoYffbar2gammagammaPhotons(nullptr),
  fHistIsoYgg2gammagammaPhotons(nullptr),
  fHistIsoYFragmentationPhotons(nullptr),

  fHistXGluonQComptonPhotons(nullptr),
  fHistXGluonQComptonPhotonsFOCAL(nullptr),
  fHistXGluonQComptonPhotonsLHCb(nullptr),
  fHistXGluonQComptonPhotonsMid(nullptr),

  fHistXGluonPtComptonPhotons(nullptr),
  fHistXGluonPtComptonPhotonsFOCAL(nullptr),
  fHistXGluonPtComptonPhotonsLHCb(nullptr),
  fHistXGluonPtComptonPhotonsMid(nullptr),

  fHistXGluonXQuark(nullptr),
  fHistXGluonXQuarkMid(nullptr),
  fHistXGluonXQuarkFocal(nullptr),
  fHistXGluonXQuarkLHCb(nullptr),

  fHistEtaPhotonEtaQuark(nullptr),
  fHistEtaPhotonEtaQuarkMid(nullptr),
  fHistEtaPhotonEtaQuarkFocal(nullptr),
  fHistEtaPhotonEtaQuarkLHCb(nullptr),

  fHistDeltaPhiPhotonQuark(nullptr),
  fHistDeltaPhiPhotonPionMid(nullptr),
  
  fHistX1vsX1Calc(nullptr),
  fHistX2vsX2Calc(nullptr),

  fHistXSecVsEvent(nullptr),

  fOutPart1(nullptr),
  fOutPart2(nullptr),
  
  fIsMC(1),
  fMaxpT(20),
  fMaxIso(50),
  fMinX(10E-9),
  fProcessCode(-1),
  fX1(-1),
  fX2(-1),
  fId1(-1),
  fId2(-1),
  fQFrac(-1),
  fWeight(-1),
  fEventCounter(0),
  fNTotEvents(1),
  fNRejectEvents(10000)
{

}

//________________________________________________________________________
AliAnalysisTaskGammaMCStudies::AliAnalysisTaskGammaMCStudies(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(nullptr),
  fHistNEvents(nullptr),
  fHistXSection(nullptr),
  fHistPtHard(nullptr),
  fHistPtYPi0(nullptr),
  fHistPtYAllPhotons(nullptr),
  fHistPtYDecayPhotons(nullptr),
  fHistPtYDirectPhotons(nullptr),
  fHistPtYComptonPhotons(nullptr),
  fHistPtYAnnihilationPhotons(nullptr),
  fHistPtYgg2ggammaPhotons(nullptr),
  fHistPtYffbar2gammagammaPhotons(nullptr),
  fHistPtYgg2gammagammaPhotons(nullptr),
  fHistPtYFragmentationPhotons(nullptr),
  
  fHistIsoYPi0(nullptr),
  fHistIsoYAllPhotons(nullptr),
  fHistIsoYDecayPhotons(nullptr),
  fHistIsoYDirectPhotons(nullptr),
  fHistIsoYComptonPhotons(nullptr),
  fHistIsoYAnnihilationPhotons(nullptr),
  fHistIsoYgg2ggammaPhotons(nullptr),
  fHistIsoYffbar2gammagammaPhotons(nullptr),
  fHistIsoYgg2gammagammaPhotons(nullptr),
  fHistIsoYFragmentationPhotons(nullptr),

  fHistXGluonQComptonPhotons(nullptr),
  fHistXGluonQComptonPhotonsFOCAL(nullptr),
  fHistXGluonQComptonPhotonsLHCb(nullptr),
  fHistXGluonQComptonPhotonsMid(nullptr),
  
  fHistXGluonPtComptonPhotons(nullptr),
  fHistXGluonPtComptonPhotonsFOCAL(nullptr),
  fHistXGluonPtComptonPhotonsLHCb(nullptr),
  fHistXGluonPtComptonPhotonsMid(nullptr),

  fHistXGluonXQuark(nullptr),
  fHistXGluonXQuarkMid(nullptr),
  fHistXGluonXQuarkFocal(nullptr),
  fHistXGluonXQuarkLHCb(nullptr),

  fHistEtaPhotonEtaQuark(nullptr),
  fHistEtaPhotonEtaQuarkMid(nullptr),
  fHistEtaPhotonEtaQuarkFocal(nullptr),
  fHistEtaPhotonEtaQuarkLHCb(nullptr),

  fHistDeltaPhiPhotonQuark(nullptr),
  fHistDeltaPhiPhotonPionMid(nullptr),

  fHistX1vsX1Calc(nullptr),
  fHistX2vsX2Calc(nullptr),
  
  fHistXSecVsEvent(nullptr),
  
  fOutPart1(nullptr),
  fOutPart2(nullptr),
  
  fIsMC(1),
  fMaxpT(20),
  fMaxIso(50),
  fMinX(10E-9),
  fX1(-1),
  fX2(-1),
  fId1(-1),
  fId2(-1),
  fQFrac(-1),
  fWeight(-1),
  fEventCounter(0),
  fNTotEvents(1),
  fNRejectEvents(10000)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaMCStudies::~AliAnalysisTaskGammaMCStudies()
{

}

//________________________________________________________________________
void AliAnalysisTaskGammaMCStudies::UserCreateOutputObjects(){

  // Create histograms
  if(fOutputContainer != nullptr){
    delete fOutputContainer;
    fOutputContainer          = nullptr;
  }
  if(fOutputContainer == nullptr){
    fOutputContainer          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  fHistNEvents                		= new TH1F("NEvents", "NEvents", 3, -0.5, 2.5);
  fHistNEvents->Sumw2();
  fOutputContainer->Add(fHistNEvents);

  fHistXSection               		= new TH1D("XSection", "XSection", 1000000, 0, 1e4);

  //   SetLogBinningXTH1(fHistXSection);
  fHistXSection->Sumw2();
  fOutputContainer->Add(fHistXSection);

  fHistPtHard                 		= new TH1F("PtHard", "PtHard", fMaxpT*20, 0, fMaxpT);
  fHistPtHard->Sumw2();
  fOutputContainer->Add(fHistPtHard);

// Spectra
  fHistPtYPi0                		= new TH2F("Pt_Y_Pi0","Pt_Y_Pi0", fMaxpT*20, 0, fMaxpT, fMaxpT, -1.0, 1.0);
  fHistPtYPi0->Sumw2();
  fOutputContainer->Add(fHistPtYPi0);

  fHistPtYAllPhotons                		= new TH2F("fHistPtYAllPhotons","all photons", fMaxpT*20, 0, fMaxpT, 200, -7, 7);
  fHistPtYAllPhotons->Sumw2();
  fOutputContainer->Add(fHistPtYAllPhotons);

  fHistPtYDecayPhotons                		= new TH2F("fHistPtYDecayPhotons","decay photons", fMaxpT*20, 0, fMaxpT, 200, -7, 7);
  fHistPtYDecayPhotons->Sumw2();
  fOutputContainer->Add(fHistPtYDecayPhotons);

  fHistPtYDirectPhotons                		= new TH2F("fHistPtYDirectPhotons","direct photons (including frag photons)", fMaxpT*20, 0, fMaxpT, 200, -7, 7);
  fHistPtYDirectPhotons->Sumw2();
  fOutputContainer->Add(fHistPtYDirectPhotons);

  fHistPtYComptonPhotons                		= new TH2F("fHistPtYComptonPhotons","photons form compton process", fMaxpT*20, 0, fMaxpT, 200, -7, 7);
  fHistPtYComptonPhotons->Sumw2();
  fOutputContainer->Add(fHistPtYComptonPhotons);
  
  fHistPtYAnnihilationPhotons                		= new TH2F("fHistPtYAnnihilationPhotons","photons form Annihilation process", fMaxpT*20, 0, fMaxpT, 200, -7, 7);
  fHistPtYAnnihilationPhotons->Sumw2();
  fOutputContainer->Add(fHistPtYAnnihilationPhotons);

  fHistPtYgg2ggammaPhotons                		= new TH2F("fHistPtYgg2ggammaPhotons","photons form gg #rightarrow g #gamma process", fMaxpT*20, 0, fMaxpT, 200, -7, 7);
  fHistPtYgg2ggammaPhotons->Sumw2();
  fOutputContainer->Add(fHistPtYgg2ggammaPhotons);

  fHistPtYffbar2gammagammaPhotons                		= new TH2F("fHistPtYffbar2gammagammaPhotons","photons form q anti q #rightarrow #gamma #gamma process", fMaxpT*20, 0, fMaxpT, 200, -7, 7);
  fHistPtYffbar2gammagammaPhotons->Sumw2();
  fOutputContainer->Add(fHistPtYffbar2gammagammaPhotons);

  fHistPtYgg2gammagammaPhotons                		= new TH2F("fHistPtYgg2gammagammaPhotons","photons form gg #rightarrow #gamma #gamma process", fMaxpT*20, 0, fMaxpT, 200, -7, 7);
  fHistPtYgg2gammagammaPhotons->Sumw2();
  fOutputContainer->Add(fHistPtYgg2gammagammaPhotons);

  fHistPtYFragmentationPhotons                		= new TH2F("fHistPtYFragmentationPhotons","fragmentation photons", fMaxpT*20, 0, fMaxpT, 200, -7, 7);
  fHistPtYFragmentationPhotons->Sumw2();
  fOutputContainer->Add(fHistPtYFragmentationPhotons);

// Isolation
  fHistIsoYPi0                		= new TH2F("Iso_Y_Pi0","Iso_Y_Pi0", fMaxIso*20, 0, fMaxIso, fMaxIso, -1.0, 1.0);
  fHistIsoYPi0->Sumw2();
  fOutputContainer->Add(fHistIsoYPi0);

  fHistIsoYAllPhotons                		= new TH2F("fHistIsoYAllPhotons","all photons", fMaxIso*20, 0, fMaxIso, 200, -7, 7);
  fHistIsoYAllPhotons->Sumw2();
  fOutputContainer->Add(fHistIsoYAllPhotons);

  fHistIsoYDecayPhotons                		= new TH2F("fHistIsoYDecayPhotons","decay photons", fMaxIso*20, 0, fMaxIso, 200, -7, 7);
  fHistIsoYDecayPhotons->Sumw2();
  fOutputContainer->Add(fHistIsoYDecayPhotons);

  fHistIsoYDirectPhotons                		= new TH2F("fHistIsoYDirectPhotons","direct photons (including frag photons)", fMaxIso*20, 0, fMaxIso, 200, -7, 7);
  fHistIsoYDirectPhotons->Sumw2();
  fOutputContainer->Add(fHistIsoYDirectPhotons);

  fHistIsoYComptonPhotons                		= new TH2F("fHistIsoYComptonPhotons","photons form compton process", fMaxIso*20, 0, fMaxIso, 200, -7, 7);
  fHistIsoYComptonPhotons->Sumw2();
  fOutputContainer->Add(fHistIsoYComptonPhotons);
  
  fHistIsoYAnnihilationPhotons                		= new TH2F("fHistIsoYAnnihilationPhotons","photons form Annihilation process", fMaxIso*20, 0, fMaxIso, 200, -7, 7);
  fHistIsoYAnnihilationPhotons->Sumw2();
  fOutputContainer->Add(fHistIsoYAnnihilationPhotons);

  fHistIsoYgg2ggammaPhotons                		= new TH2F("fHistIsoYgg2ggammaPhotons","photons form gg #rightarrow g #gamma process", fMaxIso*20, 0, fMaxIso, 200, -7, 7);
  fHistIsoYgg2ggammaPhotons->Sumw2();
  fOutputContainer->Add(fHistIsoYgg2ggammaPhotons);

  fHistIsoYffbar2gammagammaPhotons                		= new TH2F("fHistIsoYffbar2gammagammaPhotons","photons form q anti q #rightarrow #gamma #gamma process", fMaxIso*20, 0, fMaxIso, 200, -7, 7);
  fHistIsoYffbar2gammagammaPhotons->Sumw2();
  fOutputContainer->Add(fHistIsoYffbar2gammagammaPhotons);

  fHistIsoYgg2gammagammaPhotons                		= new TH2F("fHistIsoYgg2gammagammaPhotons","photons form gg #rightarrow #gamma #gamma process", fMaxIso*20, 0, fMaxIso, 200, -7, 7);
  fHistIsoYgg2gammagammaPhotons->Sumw2();
  fOutputContainer->Add(fHistIsoYgg2gammagammaPhotons);

  fHistIsoYFragmentationPhotons                		= new TH2F("fHistIsoYFragmentationPhotons","fragmentation photons", fMaxIso*20, 0, fMaxIso, 200, -7, 7);
  fHistIsoYFragmentationPhotons->Sumw2();
  fOutputContainer->Add(fHistIsoYFragmentationPhotons);

  // Process studies

  fHistXGluonQComptonPhotons                		= new TH2F("fHistXGluonQComptonPhotons","q vs x probed by compton photons", 500, fMinX,10E-1 , 1000, 0, 150);
  SetLogBinningXTH2(fHistXGluonQComptonPhotons);
  fHistXGluonQComptonPhotons->SetXTitle("x_{gluon}");
  fHistXGluonQComptonPhotons->SetYTitle("Q^{2}");
  fHistXGluonQComptonPhotons->Sumw2();
  fOutputContainer->Add(fHistXGluonQComptonPhotons);
  
  fHistXGluonQComptonPhotonsFOCAL                		= new TH2F("fHistXGluonQComptonPhotonsFOCAL","q vs x probed by compton photons in FOCAL acc.", 500, fMinX,10E-1 , 1000, 0, 150);
  SetLogBinningXTH2(fHistXGluonQComptonPhotonsFOCAL);
  fHistXGluonQComptonPhotonsFOCAL->SetXTitle("x_{gluon}");
  fHistXGluonQComptonPhotonsFOCAL->SetYTitle("Q^{2}");
  fHistXGluonQComptonPhotonsFOCAL->Sumw2();
  fOutputContainer->Add(fHistXGluonQComptonPhotonsFOCAL);

  fHistXGluonQComptonPhotonsLHCb             		= new TH2F("fHistXGluonQComptonPhotonsLHCb","q vs x probed by compton photons in LHCb acc.", 500, fMinX,10E-1 , 1000, 0, 150);   
  SetLogBinningXTH2(  fHistXGluonQComptonPhotonsLHCb );
  fHistXGluonQComptonPhotonsLHCb->SetXTitle("x_{gluon}");
  fHistXGluonQComptonPhotonsLHCb->SetYTitle("Q^{2}");
  fHistXGluonQComptonPhotonsLHCb->Sumw2();
  fOutputContainer->Add(fHistXGluonQComptonPhotonsLHCb);

  fHistXGluonQComptonPhotonsMid             		= new TH2F("fHistXGluonQComptonPhotonsMid","q vs x probed by compton photons in mid acc.", 500, fMinX,10E-1 , 1000, 0, 150);   
  SetLogBinningXTH2(  fHistXGluonQComptonPhotonsMid );
  fHistXGluonQComptonPhotonsMid->SetXTitle("x_{gluon}");
  fHistXGluonQComptonPhotonsMid->SetYTitle("Q^{2}");
  fHistXGluonQComptonPhotonsMid->Sumw2();
  fOutputContainer->Add(fHistXGluonQComptonPhotonsMid);

  fHistXGluonPtComptonPhotons                		= new TH2F("fHistXGluonPtComptonPhotons","q vs Pt probed of compton photons", 500, fMinX,10E-1 , 100, 0, 25);
  SetLogBinningXTH2(fHistXGluonPtComptonPhotons);
  fHistXGluonPtComptonPhotons->SetXTitle("x_{gluon}");
  fHistXGluonPtComptonPhotons->SetYTitle("Q^{2}");
  fHistXGluonPtComptonPhotons->Sumw2();
  fOutputContainer->Add(fHistXGluonPtComptonPhotons);
  
  fHistXGluonPtComptonPhotonsFOCAL                		= new TH2F("fHistXGluonPtComptonPhotonsFOCAL","q vs Pt probed of compton photons in FOCAL acc.", 500, fMinX,10E-1 , 100, 0, 25);
  SetLogBinningXTH2(fHistXGluonPtComptonPhotonsFOCAL);
  fHistXGluonPtComptonPhotonsFOCAL->SetXTitle("x_{gluon}");
  fHistXGluonPtComptonPhotonsFOCAL->SetYTitle("Q^{2}");
  fHistXGluonPtComptonPhotonsFOCAL->Sumw2();
  fOutputContainer->Add(fHistXGluonPtComptonPhotonsFOCAL);

  fHistXGluonPtComptonPhotonsLHCb             		= new TH2F("fHistXGluonPtComptonPhotonsLHCb","q vs Pt probed of compton photons in LHCb acc.", 500, fMinX,10E-1 , 100, 0, 25);   
  SetLogBinningXTH2(  fHistXGluonPtComptonPhotonsLHCb );
  fHistXGluonPtComptonPhotonsLHCb->SetXTitle("x_{gluon}");
  fHistXGluonPtComptonPhotonsLHCb->SetYTitle("Q^{2}");
  fHistXGluonPtComptonPhotonsLHCb->Sumw2();
  fOutputContainer->Add(fHistXGluonPtComptonPhotonsLHCb);

  fHistXGluonPtComptonPhotonsMid             		= new TH2F("fHistXGluonPtComptonPhotonsMid","q vs Pt probed of compton photons in mid acc.", 500, fMinX,10E-1 , 100, 0, 25);   
  SetLogBinningXTH2(  fHistXGluonPtComptonPhotonsMid );
  fHistXGluonPtComptonPhotonsMid->SetXTitle("x_{gluon}");
  fHistXGluonPtComptonPhotonsMid->SetYTitle("Q^{2}");
  fHistXGluonPtComptonPhotonsMid->Sumw2();
  fOutputContainer->Add(fHistXGluonPtComptonPhotonsMid);

  fHistXGluonXQuark             		= new TH2F("fHistXGluonXQuark","x quark vs x gluon probed by compton photons in all acc.", 500, fMinX,10E-1 ,500, fMinX,10E-1);   
  SetLogBinningXTH2(  fHistXGluonXQuark );
  SetLogBinningYTH2(  fHistXGluonXQuark);
  fHistXGluonXQuark->SetXTitle("x_{gluon}");
  fHistXGluonXQuark->SetYTitle("x_{quark}");
  fHistXGluonXQuark->Sumw2();
  fOutputContainer->Add(fHistXGluonXQuark);

  fHistXGluonXQuarkMid             		= new TH2F("fHistXGluonXQuarkMid","x quark vs x gluon probed by compton photons in mid acc.", 500, fMinX,10E-1 ,500, fMinX,10E-1);   
  SetLogBinningXTH2(  fHistXGluonXQuarkMid );
  SetLogBinningYTH2(  fHistXGluonXQuarkMid);
  fHistXGluonXQuarkMid->SetXTitle("x_{gluon}");
  fHistXGluonXQuarkMid->SetYTitle("x_{quark}");
  fHistXGluonXQuarkMid->Sumw2();
  fOutputContainer->Add(fHistXGluonXQuarkMid);

  fHistXGluonXQuarkFocal             		= new TH2F("fHistXGluonXQuarkFocal","x quark vs x gluon probed by compton photons in mid acc.", 500, fMinX,10E-1 ,500, fMinX,10E-1);   
  SetLogBinningXTH2(  fHistXGluonXQuarkFocal );
  SetLogBinningYTH2(  fHistXGluonXQuarkFocal);
  fHistXGluonXQuarkFocal->SetXTitle("x_{gluon}");
  fHistXGluonXQuarkFocal->SetYTitle("x_{quark}");
  fHistXGluonXQuarkFocal->Sumw2();
  fOutputContainer->Add(fHistXGluonXQuarkFocal);

  fHistXGluonXQuarkLHCb             		= new TH2F("fHistXGluonXQuarkLHCb","x quark vs x gluon probed by compton photons in mid acc.", 500, fMinX,10E-1 ,500, fMinX,10E-1);   
  SetLogBinningXTH2(  fHistXGluonXQuarkLHCb );
  SetLogBinningYTH2(  fHistXGluonXQuarkLHCb);
  fHistXGluonXQuarkLHCb->SetXTitle("x_{gluon}");
  fHistXGluonXQuarkLHCb->SetYTitle("x_{quark}");
  fHistXGluonXQuarkLHCb->Sumw2();
  fOutputContainer->Add(fHistXGluonXQuarkLHCb);

  // Outgoing partons kinematics


  fHistEtaPhotonEtaQuark             		= new TH2F("fHistEtaPhotonEtaQuark","eta quark vs eta Photon probed by compton photons in all acc.", 500, -8,8 ,500, -8,8);   
  fHistEtaPhotonEtaQuark->SetXTitle("#eta_{#gamma}");
  fHistEtaPhotonEtaQuark->SetYTitle("#eta_{quark}");
  fHistEtaPhotonEtaQuark->Sumw2();
  fOutputContainer->Add(  fHistEtaPhotonEtaQuark );

  fHistEtaPhotonEtaQuarkMid             		= new TH2F("fHistEtaPhotonEtaQuarkMid","eta quark vs eta Photon probed by compton photons in mid acc.", 500, -1,1 ,500, -1,1);   
  fHistEtaPhotonEtaQuarkMid->SetXTitle("#eta_{#gamma}");
  fHistEtaPhotonEtaQuarkMid->SetYTitle("#eta_{quark}");
  fHistEtaPhotonEtaQuarkMid->Sumw2();
  fOutputContainer->Add(  fHistEtaPhotonEtaQuarkMid );

  fHistEtaPhotonEtaQuarkFocal             		= new TH2F("fHistEtaPhotonEtaQuarkFocal","x quark vs x #gamma probed by compton photons in mid acc.", 500, 3.2,5.3 ,500, 3.2,5.3);   
  fHistEtaPhotonEtaQuarkFocal->SetXTitle("#eta_{#gamma}");
  fHistEtaPhotonEtaQuarkFocal->SetYTitle("#eta_{quark}");
  fHistEtaPhotonEtaQuarkFocal->Sumw2();
  fOutputContainer->Add(  fHistEtaPhotonEtaQuarkFocal );

  fHistEtaPhotonEtaQuarkLHCb             		= new TH2F("fHistEtaPhotonEtaQuarkLHCb","x quark vs x #gamma probed by compton photons in mid acc.", 500, 1.9,5.1 ,500, 1.9,5.1);   
  fHistEtaPhotonEtaQuarkLHCb->SetXTitle("#eta_{#gamma}");
  fHistEtaPhotonEtaQuarkLHCb->SetYTitle("#eta_{quark}");
  fHistEtaPhotonEtaQuarkLHCb->Sumw2();
  fOutputContainer->Add(  fHistEtaPhotonEtaQuarkLHCb );

  fHistDeltaPhiPhotonQuark             		= new TH2F("fHistDeltaPhiPhotonQuark","phi quark photon correlation", 500, 3.11,3.16 ,500, -8,8);   
  fHistDeltaPhiPhotonQuark->SetXTitle("#Delta #Phi}");
  fHistDeltaPhiPhotonQuark->SetYTitle("#eta_{#gamma}");
  fHistDeltaPhiPhotonQuark->Sumw2();
  fOutputContainer->Add(  fHistDeltaPhiPhotonQuark );  

  fHistDeltaPhiPhotonPionMid             		= new TH1F("fHistDeltaPhiPhotonPion","phi pion photon correlation", 500, 0.,3.16);   
  fHistDeltaPhiPhotonPionMid->SetXTitle("#Delta #Phi}");
  fHistDeltaPhiPhotonPionMid->Sumw2();
  fOutputContainer->Add(  fHistDeltaPhiPhotonPionMid );  

  // Compare calculated x to Pythia x
  fHistX1vsX1Calc             		= new TH2F("fHistX1vsX1Calc","x1 calc vs x1 calc  probed by compton photons in FOCAL acc.", 500, fMinX,10E-1 ,500, fMinX,10E-1);   
  SetLogBinningXTH2(  fHistX1vsX1Calc );
  SetLogBinningYTH2(  fHistX1vsX1Calc);
  fHistX1vsX1Calc->SetXTitle("x_{1, Pythia}");
  fHistX1vsX1Calc->SetYTitle("x_{1, Calc.}");
  fHistX1vsX1Calc->Sumw2();
  fOutputContainer->Add(fHistX1vsX1Calc);
  
  fHistX2vsX2Calc             		= new TH2F("fHistX2vsX2Calc","x2 calc vs x2 calc  probed by compton photons in FOCAL acc.", 500, fMinX,10E-1 ,500, fMinX,10E-1);   
  SetLogBinningXTH2(  fHistX2vsX2Calc );
  SetLogBinningYTH2(  fHistX2vsX2Calc);
  fHistX2vsX2Calc->SetXTitle("x_{1, Pythia}");
  fHistX2vsX2Calc->SetYTitle("x_{1, Calc.}");
  fHistX2vsX2Calc->Sumw2();
  fOutputContainer->Add(fHistX2vsX2Calc);


  // debug
  fHistXSecVsEvent             		= new TH1F("fHistXSecVsEvent","x sec vs nmb of events",50000, 0,50000);   
  fHistXSecVsEvent->SetXTitle("Event");
  fHistXSecVsEvent->SetYTitle("#sigma");
  fHistXSecVsEvent->Sumw2();
  fOutputContainer->Add(fHistXSecVsEvent);
  
  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaMCStudies::UserExec(Option_t *)
{
  fEventCounter++;
  // drop first 10000 events to get good accuracy of sec for weighting
  if(fEventCounter<10000) return;

  fInputEvent = InputEvent();
  //   cout << "I found an Event" << endl;

  fMCEvent = MCEvent();
  if(fMCEvent == nullptr) fIsMC = 0;
  if (fIsMC==0) return;
  //   cout << "I found an MC header" << endl;

  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxZ   = primVtxMC->GetZ(); 

  if (TMath::Abs(mcProdVtxZ) < 10 ){
    fHistNEvents->Fill(0);
  } else {
    fHistNEvents->Fill(1);
  }

  AliGenEventHeader* mcEH = fMCEvent->GenEventHeader();
  AliGenPythiaEventHeader *pyH  = dynamic_cast<AliGenPythiaEventHeader*>(mcEH);
  AliGenHijingEventHeader *hiH  = 0;
  AliGenDPMjetEventHeader *dpmH = 0;

  if(pyH) {
    fProcessCode = pyH->ProcessType();
    fX1          = pyH->GetX1();
    fX2          = pyH->GetX2();
    fQFrac       = pyH->GetQFrac();
   // printf("x1 = %f x2 = %f Q = %f \n",fX1,fX2,fQFrac);
  }
  // it can be only one save some casts
  // assuming PYTHIA and HIJING are the most likely ones...
  if(!pyH){
    hiH = dynamic_cast<AliGenHijingEventHeader*>(mcEH);
    if(!hiH){
      dpmH = dynamic_cast<AliGenDPMjetEventHeader*>(mcEH);
    }
  }
   
  // fetch the trials on a event by event basis, not from pyxsec.root otherwise
  // we will get a problem when running on proof since Notify may be called
  // more than once per file
  // consider storing this information in the AOD output via AliAODHandler
  Float_t ntrials = 0;
  if (!pyH || !hiH || dpmH) {
    AliGenCocktailEventHeader *ccEH = dynamic_cast<AliGenCocktailEventHeader *>(mcEH);
    if (ccEH) {
      TList *genHeaders = ccEH->GetHeaders();
      for (int imch=0; imch<genHeaders->GetEntries(); imch++) {
        if(!pyH)pyH = dynamic_cast<AliGenPythiaEventHeader*>(genHeaders->At(imch));
        if(!hiH)hiH = dynamic_cast<AliGenHijingEventHeader*>(genHeaders->At(imch));
        if(!dpmH)dpmH = dynamic_cast<AliGenDPMjetEventHeader*>(genHeaders->At(imch));
      }
    }
  }

  
  // take the trials from the p+p event
  if(hiH)ntrials = hiH->Trials();
  if(dpmH)ntrials = dpmH->Trials();
  if(pyH)ntrials = pyH->Trials();
  if(ntrials)fHistNEvents->Fill(2,ntrials);

   
  Double_t xSection = 0;
  Double_t ptHard = 0;
  if (pyH) xSection = pyH->GetXsection();
  if (pyH) ptHard = pyH->GetPtHard();
  if (xSection) fHistXSection->Fill(xSection);
  if (ptHard) fHistPtHard->Fill(ptHard);

  // calculate weight
  fWeight = xSection/(ntrials*(fNTotEvents-fNRejectEvents));


  fHistXSecVsEvent->SetBinContent(fEventCounter,xSection);
  fEventCounter++;
  //printf("XSec = %f\t NTrials=%f \t  Weight = %f \n",xSection,ntrials,fWeight);
  ProcessMCParticles();


  PostData(1, fOutputContainer);
}


//________________________________________________________________________
void AliAnalysisTaskGammaMCStudies::ProcessMCParticles()
{

  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {

    
    // fill primary histograms
    TParticle* particle         = nullptr;
    particle                    = (TParticle *)fMCEvent->Particle(i);
    if (!particle) continue;
    // Get partons of process that produced the photon
     TParticle* InPart1 = (TParticle *)fMCEvent->Particle(2); // should be incoming parton
     TParticle* InPart2 = (TParticle *)fMCEvent->Particle(3); // should be incoming parton

     fOutPart1 = fMCEvent->Particle(InPart1->GetDaughter(0)); // should be outgoing parton
     fOutPart2 = fMCEvent->Particle(InPart1->GetDaughter(1)); // should be outgoing parton

     fId1 = InPart1->GetPdgCode();
     fId2 = InPart2->GetPdgCode();
     fStatus1 = InPart1->GetStatusCode();
     fStatus2 = InPart2->GetStatusCode();

    // skip virtual particles
    if (particle->GetStatusCode()<0) continue;

    //const std::array<int, 19> kAcceptPdgCodes = {kPdgPi0, kPdgEta, kPdgEtaPrime, kPdgOmega, kPdgPiPlus, kPdgRho0, kPdgPhi, kPdgJPsi, kPdgSigma0, kPdgK0Short, kPdgDeltaPlus, kPdgDeltaPlusPlus, kPdgDeltaMinus, kPdgDelta0, kPdgRhoPlus, kPdgKStar, kPdgK0Long, kPdgLambda, kPdgKPlus};
    //if(std::find(kAcceptPdgCodes.begin(), kAcceptPdgCodes.end(), TMath::Abs(particle->GetPdgCode())) ==  kAcceptPdgCodes.end()) continue;  // species not supported

    switch(particle->GetPdgCode()){
    case kPdgPi0:
      fHistPtYPi0->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgPhoton:
      ProcessPhoton(particle);
       
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskGammaMCStudies::ProcessPhoton(TParticle* part) {

   TLorentzVector* photonMom = new TLorentzVector();
   photonMom->SetPxPyPzE(part->Px(),part->Py(),part->Pz(),part->Energy());
   Float_t photonPt = photonMom->Pt();
   fHistPtYAllPhotons->Fill(photonPt,photonMom->Rapidity(),fWeight);
   

  // Calculate x for comparison with x obtained using Pythia
  Float_t cme = 14000; // TODO: change this from being hard coded
  Float_t x1Calc = CalculateX1LO(photonMom->Rapidity(),fQFrac,cme);
  Float_t x2Calc = CalculateX2LO(photonMom->Rapidity(),fQFrac,cme);

  // calculate Et within R<0.4
  Float_t EtInCone = CalculateEt(part,0.4);
  fHistIsoYAllPhotons->Fill(EtInCone,photonMom->Rapidity(),fWeight);
   
  // check if photon is a direct photon, frag photon or decay photon

  switch(IsDirectPhoton(part)){
    case -1: 
      // virtual (ignore)
      break;
    case 0:
      // decay photon
        fHistPtYDecayPhotons->Fill(photonMom->Pt(),photonMom->Rapidity(),fWeight);
        fHistIsoYDecayPhotons->Fill(EtInCone,photonMom->Rapidity(),fWeight);
        break;
    case 1: 
        // prompt photon
        fHistPtYDirectPhotons->Fill(photonMom->Pt(),photonMom->Rapidity(),fWeight);
        fHistIsoYDirectPhotons->Fill(EtInCone,photonMom->Rapidity(),fWeight);

        // check correlations for all photons, not only compton
        
        //CalculatePhotonCorrelations(part);


      if(fProcessCode == kPromptPhotonCompton){
        fHistPtYComptonPhotons->Fill(photonMom->Pt(),photonMom->Rapidity(),fWeight);
        fHistIsoYComptonPhotons->Fill(EtInCone,photonMom->Rapidity(),fWeight);
      
        // For photons in FOCAL acceptance, what is the QSquared and x of gluon
        Float_t xgluon = -1;
        Float_t xquark = -1;
        Float_t q2     = fQFrac * fQFrac;
        if(fId1==21){
          xgluon = fX1;
          xquark = fX2;
        } else if(fId2==21){
          xgluon = fX2;
          xquark = fX1;
        }

        fHistXGluonQComptonPhotons->Fill(xgluon,q2,fWeight);
        fHistXGluonPtComptonPhotons->Fill(xgluon,photonPt,fWeight);
        fHistXGluonXQuark->Fill(xgluon,xquark,fWeight);

        // do a cross check here: daughter particles for this process should be photon and quark
        Int_t outPartID1 = TMath::Abs(fOutPart1->GetPdgCode());
        Int_t outPartID2 = TMath::Abs(fOutPart2->GetPdgCode());

        Float_t etaPhoton = -9999.; 
        Float_t etaQuark = -9999.; 
        Float_t phiPhoton = -9999.; 
        Float_t phiQuark = -9999.; 

        Float_t deltaPhi = -9999;
        if(outPartID1<9){ // it is a quark
            if(outPartID2==22){ // id2 is photon
               etaPhoton = fOutPart2->Eta();
               phiPhoton = fOutPart2->Phi();
               etaQuark  = fOutPart1->Eta();
               phiQuark  = fOutPart1->Phi();
            } else{
              printf("ERROR: something is not right during extraction of outgoing particles!\n");
              return;
            }
        } else{
          if((outPartID1==22)&&(outPartID2<9)){ // everthing correct
               etaPhoton = fOutPart1->Eta();
               phiPhoton = fOutPart1->Phi();
               etaQuark  = fOutPart2->Eta();
               phiQuark  = fOutPart2->Phi();
          } else{
               printf("ERROR: something is not right during extraction of outgoing particles!\n");
              return;
          }
        }
        deltaPhi = TMath::Abs(phiPhoton-phiQuark);

        fHistEtaPhotonEtaQuark->Fill(etaPhoton,etaQuark,fWeight);
        fHistDeltaPhiPhotonQuark->Fill(deltaPhi,etaPhoton,fWeight);

        

        if(IsInFOCALAcceptance(part)){
            fHistXGluonQComptonPhotonsFOCAL->Fill(xgluon,q2,fWeight);
            fHistXGluonPtComptonPhotonsFOCAL->Fill(xgluon,photonPt,fWeight);
            fHistXGluonXQuarkFocal->Fill(xgluon,xquark);
            fHistEtaPhotonEtaQuarkFocal->Fill(etaPhoton,etaQuark,fWeight);

            fHistX1vsX1Calc->Fill(fX1,x1Calc,fWeight);
            fHistX2vsX2Calc->Fill(fX2,x2Calc,fWeight);
        }
        if(IsInLHCbAcceptance(part)){
          fHistXGluonQComptonPhotonsLHCb->Fill(xgluon,q2,fWeight);
          fHistXGluonPtComptonPhotonsLHCb->Fill(xgluon,photonPt,fWeight);
          fHistXGluonXQuarkLHCb->Fill(xgluon,xquark,fWeight);
          fHistEtaPhotonEtaQuarkLHCb->Fill(etaPhoton,etaQuark,fWeight);
        }

        if(IsInMidAcceptance(part)){
          fHistXGluonQComptonPhotonsMid->Fill(xgluon,q2,fWeight);
          fHistXGluonPtComptonPhotonsMid->Fill(xgluon,photonPt,fWeight);
          fHistXGluonXQuarkMid->Fill(xgluon,xquark,fWeight);
          fHistEtaPhotonEtaQuarkMid->Fill(etaPhoton,etaQuark,fWeight);
        }
      } else if (fProcessCode == kPromptPhotonAnnihilation){
          fHistPtYAnnihilationPhotons->Fill(photonMom->Pt(),photonMom->Rapidity(),fWeight);
          fHistIsoYAnnihilationPhotons->Fill(EtInCone,photonMom->Rapidity(),fWeight);
      } else if (fProcessCode == kPromptPhotongg2ggamma){
          fHistPtYgg2ggammaPhotons->Fill(photonMom->Pt(),photonMom->Rapidity(),fWeight);
          fHistIsoYgg2ggammaPhotons->Fill(EtInCone,photonMom->Rapidity(),fWeight);
      } else if (fProcessCode == kPromptPhotonffbar2gammagamma){
          fHistPtYffbar2gammagammaPhotons->Fill(photonMom->Pt(),photonMom->Rapidity(),fWeight);
          fHistIsoYffbar2gammagammaPhotons->Fill(EtInCone,photonMom->Rapidity(),fWeight);
      } else if (fProcessCode == kPromptPhotongg2gammagamma){
          fHistPtYgg2gammagammaPhotons->Fill(photonMom->Pt(),photonMom->Rapidity(),fWeight);
          fHistIsoYgg2gammagammaPhotons->Fill(EtInCone,photonMom->Rapidity(),fWeight);
      } else{
        // put unknown case here
      }
      break;
    case 2:
        fHistPtYDirectPhotons->Fill(photonMom->Pt(),photonMom->Rapidity(),fWeight);
        fHistIsoYDirectPhotons->Fill(EtInCone,photonMom->Rapidity(),fWeight);
      // frag photon
        fHistPtYFragmentationPhotons->Fill(photonMom->Pt(),photonMom->Rapidity(),fWeight);
        fHistIsoYFragmentationPhotons->Fill(EtInCone,photonMom->Rapidity(),fWeight);
      break;
  }

  delete photonMom;

  return;
}

//________________________________________________________________________
// -1: unknown
//  0: decay photon
//  1: prompt photon
//  2: fragmentation photon
Int_t AliAnalysisTaskGammaMCStudies::IsDirectPhoton(TParticle* photon) const {
  // ask for status code
   Int_t status = TMath::Abs(photon->GetStatusCode());


  if(status>90){ // decay photon
    return 0;

  } else if( status<90){ // direct photon
     // determine wether or not prompt or frag photon is present
     TParticle* initialPhoton = GetInitialPhoton(photon);
     
     // check status code of initial photon from chain
     Int_t initialStatus = TMath::Abs(initialPhoton->GetStatusCode());

     if((initialStatus>21) && (initialStatus<40)){
       // particle from hardest subprocess or subsequent subprocess -> prompt photon
       return 1;
     } else if ((initialStatus>40) && (initialStatus<60)){
       // particle from initial-state shower or final state shower -> fragmentation photon
       return 2;
     } else{
       // untreated case
       return -1;
     }
  } else{
    // untreated case
    return -1;
  }
}


//________________________________________________________________________
// follow the chain of photons to initial photon that has a none photon mother
TParticle* AliAnalysisTaskGammaMCStudies::GetInitialPhoton(TParticle* part) const{
   if(part->GetMother(0)>-1){ // photon has a mother
      // Check if mother is a photon
      TParticle* motherParticle            = (TParticle *) fMCEvent->Particle(part->GetMother(0));
      if(TMath::Abs(motherParticle->GetPdgCode()) == kPdgPhoton){
        // Mother is a photon, check again
        return GetInitialPhoton(motherParticle);
      } else{
        // found a mother that is not a photon again
        return part;
      }

   } else{ // photon has no more mothers
      return part;
   }
}
//________________________________________________________________________
bool AliAnalysisTaskGammaMCStudies::IsInPCMAcceptance(TParticle* part) const {
  const Double_t kBoundaryEta = 0.900001;
  if (//part->Pt() > 0.050 
  //&& 
  TMath::Abs(part->Eta()) < kBoundaryEta) return true;

  return false;
}

//________________________________________________________________________
bool AliAnalysisTaskGammaMCStudies::IsInPHOSAcceptance(TParticle* part) const {
  const Double_t kBoundaryEtaMin = -0.13;
  const Double_t kBoundaryEtaMax = 0.13;
  const Double_t kBoundaryPhiMin = 4.54;
  const Double_t kBoundaryPhiMax = 5.59;
  //if (part->Pt() < 0.300) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  if (part->Phi() > kBoundaryPhiMax || part->Phi() < kBoundaryPhiMin) return false;
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskGammaMCStudies::IsInEMCalAcceptance(TParticle* part) const {
  const Double_t kBoundaryEtaMin = -0.6687;
  const Double_t kBoundaryEtaMax = 0.66465;
  const Double_t kBoundaryPhiMin = 1.39626;
  const Double_t kBoundaryPhiMax = 3.15;
  //if (part->Pt() < 0.400) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  if (part->Phi() > kBoundaryPhiMax || part->Phi() < kBoundaryPhiMin) return false;
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskGammaMCStudies::IsInFOCALAcceptance(TParticle* part) const {
  const Double_t kBoundaryEtaMin = 3.2;
  const Double_t kBoundaryEtaMax = 5.3;
  //if (part->Pt() < 0.400) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskGammaMCStudies::IsInLHCbAcceptance(TParticle* part) const {
  const Double_t kBoundaryEtaMin = 1.9;
  const Double_t kBoundaryEtaMax = 5.1;
  //if (part->Pt() < 0.400) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskGammaMCStudies::IsInMidAcceptance(TParticle* part) const {
  const Double_t kBoundaryEtaMin = -1;
  const Double_t kBoundaryEtaMax = 1.;
  //if (part->Pt() < 0.400) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  return true;
}
//________________________________________________________________________
void AliAnalysisTaskGammaMCStudies::Terminate(const Option_t *)
{

  //fOutputContainer->Print(); // Will crash on GRID
}


//_________________________________________________________________________________
void AliAnalysisTaskGammaMCStudies::SetLogBinningXTH1(TH1* histoRebin){
  TAxis *axisafter = histoRebin->GetXaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaMCStudies::SetLogBinningXTH2(TH2* histoRebin){
  TAxis *axisafter = histoRebin->GetXaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}

//_________________________________________________________________________________
void AliAnalysisTaskGammaMCStudies::SetLogBinningYTH2(TH2* histoRebin){
  TAxis *axisafter = histoRebin->GetYaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}
//________________________________________________________________________
Float_t AliAnalysisTaskGammaMCStudies::CalculateEt(TParticle* part, Float_t rmax) {
  
  Float_t Et = 0;

  TLorentzVector* p1 = new TLorentzVector();
  p1->SetPxPyPzE(part->Px(),part->Py(),part->Pz(),part->Energy());
  
  // loop over all particles
    for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {

      TParticle* otherParticle         = nullptr;
      otherParticle                    = (TParticle *)fMCEvent->Particle(i);

      if (!otherParticle) continue;
      if(otherParticle->GetStatusCode()<0) continue; // virtual -> do nothing
      if(part==otherParticle) continue;
      
      TLorentzVector p2;
      p2.SetPxPyPzE(otherParticle->Px(),otherParticle->Py(),otherParticle->Pz(),otherParticle->Energy());

      // Calculate Distance
      Float_t deltaR = p1->DeltaR(p2);
      if(deltaR<rmax) Et+=p2.Et();

      
      // check that both the photon and the pion are in givin acceptacne
      if(TMath::Abs(otherParticle->GetPdgCode())==kPdgPiPlus){
          // mid rapidity
          if(IsInMidAcceptance(part) && IsInMidAcceptance(otherParticle)){
            // Calculate Correlation
            Float_t deltaPhi = p1->DeltaPhi(p2);
            fHistDeltaPhiPhotonPionMid->Fill(deltaPhi,fWeight);
          }
      }

    }

  delete p1;
  return Et;
}

//________________________________________________________________________
Float_t AliAnalysisTaskGammaMCStudies::CalculateX1LO(Float_t y, Float_t q, Float_t cme = 14000) {
  Float_t x1 = (q*TMath::Exp(TMath::Abs(y)))/(TMath::Sqrt(cme));
  return x1;
}

//________________________________________________________________________
Float_t AliAnalysisTaskGammaMCStudies::CalculateX2LO(Float_t y, Float_t q, Float_t cme = 14000) {
  Float_t x2 = (q*TMath::Exp((-1)*TMath::Abs(y)))/(TMath::Sqrt(cme));
  return x2;
}

//________________________________________________________________________
void AliAnalysisTaskGammaMCStudies::CalculatePhotonCorrelations(TParticle* part) {
  
  TLorentzVector* p1 = new TLorentzVector();
  p1->SetPxPyPzE(part->Px(),part->Py(),part->Pz(),part->Energy());
  
  // loop over all particles
    for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {

      TParticle* otherParticle         = nullptr;
      otherParticle                    = (TParticle *)fMCEvent->Particle(i);

      if (!otherParticle) continue;

      // do not correlate with same particle (they should point to the same adress in memory)
      if(part==otherParticle) continue;

      if(otherParticle->GetStatusCode()<0) continue; // virtual -> do nothing

      TLorentzVector p2;
      p2.SetPxPyPzE(otherParticle->Px(),otherParticle->Py(),otherParticle->Pz(),otherParticle->Energy());

      // Calculate Correlation
      Float_t deltaPhi = p1->DeltaPhi(p2);

      // check that both the photon and the pion are in givin acceptacne
      if(TMath::Abs(otherParticle->GetPdgCode())==kPdgPiPlus){
          // mid rapidity
          if(IsInMidAcceptance(part) && IsInMidAcceptance(otherParticle)){
            fHistDeltaPhiPhotonPionMid->Fill(deltaPhi,fWeight);
          }
      }
    }
  delete p1;
}