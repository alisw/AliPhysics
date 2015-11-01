/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////////////
// Author: M. Nicassio m.nicassio@gsi.de maria.nicassio@cern.ch 
// First version 
// October 2014
//
// 
// Parts of this code are taken and adapted from the code below 
//
// PWGCF/FEMTOSCOPY/AliFemto and AliFemtoUser (AliFemto framework)
// PWGCF/FEMTOSCOPY/V0LamAnalysis (Jai)
// PWGCF/FEMTOSCOPY/Chaoticity (Dhevan)
// PWGCF/FEMTOSCOPY/K0S/PLamAnalysis (Hans) 
//  
// PWGLF/STRANGENESS 
//
// Possibility to run on
// - first particle: pion (0), proton (1)
// - second particle: pion (0) only when first particle pion (0), 
//   proton (1) only when first particle proton (1), Xi (2), Omega (3), 
//   to be fixed --> randomization, THD for pions
// - AODs
// - ESDs (not mantained)
//
// To be done
// - read MC truth
// . ....
/////////////////////////////////////////////////////////////////////////

class AliESDVertex;

class AliAODVertex;

#include "TChain.h"

#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TDatabasePDG.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliInputEventHandler.h"

#include "AliCentrality.h"
#include "AliPIDResponse.h"

#include "AliVTrack.h"

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDcascade.h"
#include "AliESDv0.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODcascade.h"
#include "AliAODv0.h"
#include "AliAODMCParticle.h"


#include "AliStack.h"
#include "AliMCEvent.h"

#include "AliVVertex.h"

#include "AliAnalysisTaskhCascadeFemto.h"

#define PI 3.1415927
 
using namespace std;


ClassImp(AliAnalysisTaskhCascadeFemto)

AliAnalysisTaskhCascadeFemto::AliAnalysisTaskhCascadeFemto():AliAnalysisTaskSE(),
  fESDevent(NULL),
  fAODevent(NULL),
  fAnalysisType("ESD"),
  fFirstpart(kPion),
  fSecondpart(kXi),
  fMassWindowCascades(0.008),
  fnSigmaTPCPIDfirstParticle(3.),
  fnSigmaTPCTOFPIDfirstParticle(3.),
  fnSigmaTPCPIDsecondParticleDau(4.),
  fkCascadeSideBands(kFALSE),
  fkTightCutsForCascades(kTRUE),
  
  fReadMCTruth(kFALSE),
  fUseContainer(kFALSE),
  fUseStandardCuts(0),
  fkApplyTtc(0),
  fDphisMin(0),
  fDetasMin(0),
  fMomemtumLimitForTOFPID(0.8),
  fkApplyRatioCrRnFindCut(kFALSE),
  fkCutOnTPCIP(kFALSE),  
  fIPCutxy(0.1),
  fIPCutz(0.15),
  fMinPtForPrim(0.7),
  fMaxPtForPrim(4.),
  fMinPtForCasc(0.8),
  fMaxPtForCasc(10.),
  fIPCutBac(0.1),
  fkApplyYcutCasc(kFALSE),
  fkPropagateGlobal(kFALSE),

  fESDtrackCuts(0),
  fPIDResponse(0),

  fCentrLowLim(0),
  fCentrUpLim(0),

  fNbinsMass(0),
  fLowMassLim(0),
  fUpMassLim(0),

  fNbinsMassGm(0),
  fLowMassLimGm(0),
  fUpMassLimGm(0),

  farrGT(0),fTrackBufferSize(20200), // was 18k
  fEventColl(0x0),
  fEvt(0x0),
  fMaxPMult(3000), // 1000 for protons 
  fMaxXiMult(20),  // was 100
  fPDGXi(1.32171),
  fPDGOmega(1.6724),
  fPDGp(0.938272046),
  fPDGL(1.115683),
  fPDGpi(0.13957),
  fPDGfirst(0.),
  fPDGsecond(0.),
  fzVertexBins(10),
  fnCentBins(20),
  fnEventsToMix(7),


  fHistCentrality(0),
  fHistVertexDistribution(0),
  fHistMultiplicityOfMixedEvent(0),

  fHistnTPCCrossedR(0),
  fHistRationTPCCrossedRnFind(0),
  fHistSharedFrTPCcl(0),
  fHistprimpTPCdEdx(0),

  fHistyptProtons(0),
  fHistphietaProtons(0),
  fHistIPtoPVxyzTPC(0),
  fHistIPtoPVxyzGlobal(0),

  fHistpTOFmisvspt(0),
  fHistpTOFmisvsp(0),
  fHistpTOFnsigmavspt(0),
  fHistpTOFnsigmavsp(0),
  fHistpTOFsignalvsp(0),
  fHistpTOFsignalvspt(0),

  fHistpTOFTPCsignalvspt(0),
  fHistProtonMultvsCent(0),

  fHistpTPCdEdx(0),
  fHistnTPCdEdx(0),
  fHistbTPCdEdx(0),

  fHistPosV0TPCClusters(0),
  fHistNegV0TPCClusters(0),
  fHistBachTPCClusters(0),

  fHistSharedFrTPCclPos(0),
  fHistSharedFrTPCclNeg(0),
  fHistSharedFrTPCclBach(0),

  fHistyptXi(0),
  fHistphietaXi(0),
  fHistptL(0),
  fHistptpL(0),
  fHistptnL(0),
  fHistptbac(0),

  fHistInvMassXiMinus(0),
  fHistInvMassL(0),
  fHistInvMassXiPlus(0),
  fHistInvMassAntiL(0),
  fHistInvMassXiInPairs(0),
  fHistXiMultvsCent(0),
  fCFContCascadeCuts(0),

  fHistpXiSignalRealKstar(0),
  fHistapXiSignalRealKstar(0),
  fHistpaXiSignalRealKstar(0),
  fHistapaXiSignalRealKstar(0),

  fHistpXiSignalBkgKstar(0),
  fHistapXiSignalBkgKstar(0),
  fHistpaXiSignalBkgKstar(0),
  fHistapaXiSignalBkgKstar(0),
  fHistFractionOfXiWithSharedDaughters(0),
  fHistFractionOfaXiWithSharedDaughters(0),

  fHistpXibacDetaSDphiS(0),
  fHistpXiposDetaSDphiS(0),
  fHistpXinegDetaSDphiS(0),
  fHistpaXibacDetaSDphiS(0),
  fHistpaXiposDetaSDphiS(0),
  fHistpaXinegDetaSDphiS(0),
  fHistapXibacDetaSDphiS(0),
  fHistapXiposDetaSDphiS(0),
  fHistapXinegDetaSDphiS(0),
  fHistapaXibacDetaSDphiS(0),
  fHistapaXiposDetaSDphiS(0),
  fHistapaXinegDetaSDphiS(0),

  fHistpXibacDetaSDphiSBkg(0),
  fHistpXiposDetaSDphiSBkg(0),
  fHistpXinegDetaSDphiSBkg(0),
  fHistpaXibacDetaSDphiSBkg(0),
  fHistpaXiposDetaSDphiSBkg(0),
  fHistpaXinegDetaSDphiSBkg(0),
  fHistapXibacDetaSDphiSBkg(0),
  fHistapXiposDetaSDphiSBkg(0),
  fHistapXinegDetaSDphiSBkg(0),
  fHistapaXibacDetaSDphiSBkg(0),
  fHistapaXiposDetaSDphiSBkg(0),
  fHistapaXinegDetaSDphiSBkg(0),

  fHistTrackBufferOverflow(0),

  fOutputContainer(NULL)

{
  fPDGp = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  fPDGpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();

  fPDGXi = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  fPDGOmega = TDatabasePDG::Instance()->GetParticle(3334)->Mass();
  fPDGL = TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  

  // default Constructor
  // Define input and output slots here
}

//-------------------------------------------------------------------------------

AliAnalysisTaskhCascadeFemto::AliAnalysisTaskhCascadeFemto(const char *name):AliAnalysisTaskSE(name),
  fESDevent(NULL),
  fAODevent(NULL),
  fAnalysisType("ESD"),
  fFirstpart(kPion),
  fSecondpart(kXi),
  fMassWindowCascades(0.008),
  fnSigmaTPCPIDfirstParticle(3.),
  fnSigmaTPCTOFPIDfirstParticle(3.),
  fnSigmaTPCPIDsecondParticleDau(4.),
  fkCascadeSideBands(kFALSE),
  fkTightCutsForCascades(kTRUE),
 
  fReadMCTruth(kFALSE),
  fUseContainer(kFALSE),
  fUseStandardCuts(0),
  fkApplyTtc(0),
  fDphisMin(0),
  fDetasMin(0),
  fMomemtumLimitForTOFPID(0.8),
  fkApplyRatioCrRnFindCut(kFALSE),
  fkCutOnTPCIP(kFALSE),
  fIPCutxy(0.1),
  fIPCutz(0.15),
  fMinPtForPrim(0.7),
  fMaxPtForPrim(4.),
  fMinPtForCasc(0.8),
  fMaxPtForCasc(10.),
  fIPCutBac(0.1),
  fkApplyYcutCasc(kFALSE),
  fkPropagateGlobal(kFALSE),

  fESDtrackCuts(0),
  fPIDResponse(0),

  fCentrLowLim(0),
  fCentrUpLim(0),                        
  
  fNbinsMass(0),                          
  fLowMassLim(0),
  fUpMassLim(0),
  
  fNbinsMassGm(0),                                             
  fLowMassLimGm(0),
  fUpMassLimGm(0),

  farrGT(0),fTrackBufferSize(20200), // was 18k
  fEventColl(0x0),
  fEvt(0x0),
  fMaxPMult(3000), // 1000 for protons 
  fMaxXiMult(20),  // was 100
  fPDGXi(1.32171),
  fPDGOmega(1.6724), 
  fPDGp(0.938272046),
  fPDGL(1.115683),
  fPDGpi(0.13957),
  fPDGfirst(0.),
  fPDGsecond(0.),
  fzVertexBins(10), 
  fnCentBins(20),
  fnEventsToMix(7),


  fHistCentrality(0),
  fHistVertexDistribution(0),
  fHistMultiplicityOfMixedEvent(0),

  fHistnTPCCrossedR(0),
  fHistRationTPCCrossedRnFind(0),
  fHistSharedFrTPCcl(0),
  fHistprimpTPCdEdx(0),
 
  fHistyptProtons(0),
  fHistphietaProtons(0),
  fHistIPtoPVxyzTPC(0),
  fHistIPtoPVxyzGlobal(0),

  fHistpTOFmisvspt(0),
  fHistpTOFmisvsp(0),
  fHistpTOFnsigmavspt(0),
  fHistpTOFnsigmavsp(0),
  fHistpTOFsignalvsp(0),
  fHistpTOFsignalvspt(0),

  fHistpTOFTPCsignalvspt(0),
  fHistProtonMultvsCent(0),

  fHistpTPCdEdx(0),
  fHistnTPCdEdx(0),
  fHistbTPCdEdx(0),

  fHistPosV0TPCClusters(0),
  fHistNegV0TPCClusters(0),
  fHistBachTPCClusters(0),

  fHistSharedFrTPCclPos(0),
  fHistSharedFrTPCclNeg(0),
  fHistSharedFrTPCclBach(0),

  fHistyptXi(0),
  fHistphietaXi(0),
  fHistptL(0),
  fHistptpL(0),
  fHistptnL(0),
  fHistptbac(0),

  fHistInvMassXiMinus(0),
  fHistInvMassL(0),
  fHistInvMassXiPlus(0),
  fHistInvMassAntiL(0),
  fHistInvMassXiInPairs(0),
  fHistXiMultvsCent(0),
  fCFContCascadeCuts(0),

  fHistpXiSignalRealKstar(0),
  fHistapXiSignalRealKstar(0),
  fHistpaXiSignalRealKstar(0),
  fHistapaXiSignalRealKstar(0),

  fHistpXiSignalBkgKstar(0),
  fHistapXiSignalBkgKstar(0),
  fHistpaXiSignalBkgKstar(0),
  fHistapaXiSignalBkgKstar(0),
  fHistFractionOfXiWithSharedDaughters(0),
  fHistFractionOfaXiWithSharedDaughters(0),

  fHistpXibacDetaSDphiS(0),
  fHistpXiposDetaSDphiS(0),
  fHistpXinegDetaSDphiS(0),
  fHistpaXibacDetaSDphiS(0),
  fHistpaXiposDetaSDphiS(0),
  fHistpaXinegDetaSDphiS(0),
  fHistapXibacDetaSDphiS(0),
  fHistapXiposDetaSDphiS(0),
  fHistapXinegDetaSDphiS(0),
  fHistapaXibacDetaSDphiS(0),
  fHistapaXiposDetaSDphiS(0),
  fHistapaXinegDetaSDphiS(0),

  fHistpXibacDetaSDphiSBkg(0),
  fHistpXiposDetaSDphiSBkg(0),
  fHistpXinegDetaSDphiSBkg(0),
  fHistpaXibacDetaSDphiSBkg(0),
  fHistpaXiposDetaSDphiSBkg(0),
  fHistpaXinegDetaSDphiSBkg(0),
  fHistapXibacDetaSDphiSBkg(0),
  fHistapXiposDetaSDphiSBkg(0),
  fHistapXinegDetaSDphiSBkg(0),
  fHistapaXibacDetaSDphiSBkg(0),
  fHistapaXiposDetaSDphiSBkg(0),
  fHistapaXinegDetaSDphiSBkg(0),

  fHistTrackBufferOverflow(0),

  fOutputContainer(NULL)
{
  fPDGp = TDatabasePDG::Instance()->GetParticle(2212)->Mass(); 
  fPDGpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  fPDGXi = TDatabasePDG::Instance()->GetParticle(3312)->Mass();
  fPDGOmega = TDatabasePDG::Instance()->GetParticle(3334)->Mass();
  fPDGL = TDatabasePDG::Instance()->GetParticle(3122)->Mass(); 
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, AliCFContainer::Class());

}

//--------------------------------------------------------------------------------

AliAnalysisTaskhCascadeFemto::~AliAnalysisTaskhCascadeFemto() {
  //
  // Destructor
  //
  if (fOutputContainer && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())   { delete fOutputContainer;     fOutputContainer = 0x0;    }
  if (fCFContCascadeCuts && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) { delete fCFContCascadeCuts; fCFContCascadeCuts = 0x0; }
  if (fESDtrackCuts) delete fESDtrackCuts;


  if (farrGT)
    delete[] farrGT;
  farrGT=0;


  for (unsigned short i=0; i<fzVertexBins; i++) {
    for (unsigned short j=0; j<fnCentBins; j++) {
      delete fEventColl[i][j];
    }
    delete[] fEventColl[i];
  }
  delete[] fEventColl;

}

//--------------------------------------------------------------------------------

void AliAnalysisTaskhCascadeFemto::UserCreateOutputObjects() {

  // PID object
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  // Event mixing part

  fEventColl = new AliAnalysishCascadeEventCollection **[fzVertexBins]; 

  for (unsigned short i=0; i<fzVertexBins; i++) {
    fEventColl[i] = new AliAnalysishCascadeEventCollection *[fnCentBins];
    for (unsigned short j=0; j<fnCentBins; j++) {
      fEventColl[i][j] = new AliAnalysishCascadeEventCollection(fnEventsToMix+1, fMaxXiMult, fMaxPMult);
    }
  }

  // Set PDG for first and second particle
  if (fFirstpart == kPion) fPDGfirst = fPDGpi; 
  else if (fFirstpart == kProton) fPDGfirst = fPDGp;
  else { cout<<" First particle not known or not possible to deal with in this task for the moment! "<<endl; return;}
  if (fSecondpart == kPion) fPDGsecond = fPDGpi;
  else if (fSecondpart == kProton) fPDGsecond = fPDGp;
  else if (fSecondpart == kXi) fPDGsecond = fPDGXi;
  else if (fSecondpart == kOmega) fPDGsecond = fPDGOmega;
  else { cout<<" Second particle not known or not possible to deal with in this task for the moment! "<<endl; return;}


  // Store pointer to global tracks
  farrGT = new Int_t[fTrackBufferSize];
  
  //Define and fill the OutputContainer
  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);


  // Create histograms 
  fHistCentrality = new TH1F("fHistCentrality", "Number of events", 20, 0., 100.);
  fHistCentrality ->GetXaxis()->SetTitle("Centrality");
  fHistCentrality ->GetYaxis()->SetTitle("Entries");
  fOutputContainer->Add(fHistCentrality);

  fHistVertexDistribution = new TH1F("fHistVertexDistribution", "Primary vertex distribution", 40, -20., 20.);
  fHistVertexDistribution ->GetXaxis()->SetTitle("z_{v} (cm)");
  fHistVertexDistribution ->GetYaxis()->SetTitle("Entries");
  fOutputContainer->Add(fHistVertexDistribution); 

  fHistMultiplicityOfMixedEvent = new TH1F("fHistMultiplicityOfMixedEvent","",100, 0.,100.);
  fOutputContainer->Add(fHistMultiplicityOfMixedEvent);

  fHistnTPCCrossedR = new TH1F("fHistnTPCCrossedR","",200, 0.,200.);
  fOutputContainer->Add(fHistnTPCCrossedR);
  fHistRationTPCCrossedRnFind = new TH1F("fHistRationTPCCrossedRnFind","",2000, 0.,170.);
  fOutputContainer->Add(fHistRationTPCCrossedRnFind);
  fHistSharedFrTPCcl = new TH1F("fHistSharedFrTPCcl","",400, 0.,2.);
  fOutputContainer->Add(fHistSharedFrTPCcl);

  fHistprimpTPCdEdx = new TH2F("fHistprimpTPCdEdx", "", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistprimpTPCdEdx);

  fHistyptProtons = new TH3F("fHistyptProtons","",100,0.,10.,40,-2.,2.,10,0.,100.);
  fOutputContainer->Add(fHistyptProtons);
  fHistphietaProtons = new TH2F("fHistphietaProtons","",400,0.,2*TMath::Pi(),100, -2, 2.);
  fOutputContainer->Add(fHistphietaProtons);
  fHistIPtoPVxyzTPC = new TH2F("fHistIPtoPVxyzTPC","",200, -5.,5., 200, -5., 5.);
  fOutputContainer->Add(fHistIPtoPVxyzTPC);
  fHistIPtoPVxyzGlobal = new TH2F("fHistIPtoPVxyzGlobal","",200, -5.,5., 200, -5., 5.);
  fOutputContainer->Add(fHistIPtoPVxyzGlobal);

  fHistpTOFmisvspt = new TH2F("fHistpTOFmisvspt", "", 200, 0., 10., 101, 0.0, 1.01);
  fOutputContainer->Add(fHistpTOFmisvspt);
  fHistpTOFmisvsp = new TH2F("fHistpTOFmisvsp", "", 200, 0., 10., 101, 0.0, 1.01);
  fOutputContainer->Add(fHistpTOFmisvsp);
  fHistpTOFnsigmavspt = new TH2F("fHistpTOFnsigmavspt", "", 200, 0., 10., 100, -5. , 5);
  fOutputContainer->Add(fHistpTOFnsigmavspt);
  fHistpTOFnsigmavsp = new TH2F("fHistpTOFnsigmavsp", "", 200, 0., 10., 100, -5., 5.);
  fOutputContainer->Add(fHistpTOFnsigmavsp);
  fHistpTOFsignalvsp = new TH2F("fHistpTOFsignalvsp", "", 200, 0., 10, 100,-3000,3000);//1000 , 0., 2.);
  fHistpTOFsignalvsp->GetYaxis()->SetTitle("t_{meas}-t_{0}-t_{exoected} (ps)");
  fHistpTOFsignalvsp->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fOutputContainer->Add(fHistpTOFsignalvsp);
  fHistpTOFsignalvspt = new TH2F("fHistpTOFsignalvspt", "", 200, 0., 10., 100,-3000,3000);//1000, 0.0, 2.);
  fHistpTOFsignalvspt->GetYaxis()->SetTitle("t_{meas}-t_{0}-t_{expected} (ps)");
  fHistpTOFsignalvspt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})"); 
  fOutputContainer->Add(fHistpTOFsignalvspt);

  fHistpTOFTPCsignalvspt = new TH2F("fHistpTOFTPCsignalvspt","",200, 0., 10., 100, 0., 20);
  fOutputContainer->Add(fHistpTOFTPCsignalvspt);
  fHistProtonMultvsCent = new TH2F("fHistProtonMultvsCent","",400,0.,2000.,10,0.,100.);
  fOutputContainer->Add(fHistProtonMultvsCent);
 
  fHistpTPCdEdx = new TH2F("fHistpTPCdEdx", "", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistpTPCdEdx);
  fHistnTPCdEdx = new TH2F("fHistnTPCdEdx", "", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistnTPCdEdx);
  fHistbTPCdEdx = new TH2F("fHistbTPCdEdx", "", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistbTPCdEdx);

  fHistPosV0TPCClusters = new TH1F("fHistPosV0TPCClusters", "", 100, 0.,200.);
  fOutputContainer->Add(fHistPosV0TPCClusters);
  fHistNegV0TPCClusters = new TH1F("fHistNegV0TPCClusters", "", 100, 0.,200.);
  fOutputContainer->Add(fHistNegV0TPCClusters);
  fHistBachTPCClusters = new TH1F("fHistBachTPCClusters", "", 100, 0.,200.);
  fOutputContainer->Add(fHistBachTPCClusters);

  fHistSharedFrTPCclPos = new TH1F("fHistSharedFrTPCclPos","",400, 0.,2.);
  fOutputContainer->Add(fHistSharedFrTPCclPos);
  fHistSharedFrTPCclNeg = new TH1F("fHistSharedFrTPCclNeg","",400, 0.,2.); 
  fOutputContainer->Add(fHistSharedFrTPCclNeg);
  fHistSharedFrTPCclBach = new TH1F("fHistSharedFrTPCclBach","",400, 0.,2.);
  fOutputContainer->Add(fHistSharedFrTPCclBach);

  fHistyptXi = new TH2F("fHistyptXi","",100,0.,10.,40,-2.,2.);
  fOutputContainer->Add(fHistyptXi);
  fHistphietaXi = new TH2F("fHistphietaXi","",400,0.,2*TMath::Pi(),100, -2, 2.);
  fOutputContainer->Add(fHistphietaXi);

  fHistptL = new TH1F("fHistptL","",100,0.,10.);
  fOutputContainer->Add(fHistptL);
  fHistptpL = new TH1F("fHistptpL","",100,0.,10.);
  fOutputContainer->Add(fHistptpL);
  fHistptnL = new TH1F("fHistptnL","",100,0.,10.);
  fOutputContainer->Add(fHistptnL);
  fHistptbac = new TH1F("fHistptbac","",100,0.,10.);
  fOutputContainer->Add(fHistptbac);

  fHistInvMassXiMinus = new TH2F("fHistInvMassXiMinus","", fNbinsMassGm, fLowMassLimGm, fUpMassLimGm,200,0.,10.);
  fOutputContainer->Add(fHistInvMassXiMinus);

  fHistInvMassL = new TH2F("fHistInvMassL","", fNbinsMass, fLowMassLim, fUpMassLim,100, 0., 100.);
  fOutputContainer->Add(fHistInvMassL);
  fHistInvMassXiPlus = new TH2F("fHistInvMassXiPlus","", fNbinsMassGm, fLowMassLimGm, fUpMassLimGm,200,0.,10.);
  fOutputContainer->Add(fHistInvMassXiPlus);
  fHistInvMassAntiL = new TH2F("fHistInvMassAntiL","", fNbinsMass, fLowMassLim, fUpMassLim,100, 0., 100.);
  fOutputContainer->Add(fHistInvMassAntiL);
  fHistXiMultvsCent = new TH2F("fHistXiMultvsCent","",200,0.,200.,10,0.,100.);
  fOutputContainer->Add(fHistXiMultvsCent);

  fHistInvMassXiInPairs = new TH2F("fHistInvMassXiInPairs","", fNbinsMassGm, fLowMassLimGm, fUpMassLimGm,200,0.,100.);
  fOutputContainer->Add(fHistInvMassXiInPairs);

  if(! fCFContCascadeCuts && fUseContainer) {

   // Container meant to store all the relevant distributions corresponding to the cut variables.
   // - NB overflow/underflow of variables on which we want to cut later should be 0!!!
   const Int_t lNbSteps = 2 ;
   const Int_t lNbVariables = 22 ;
   //array for the number of bins in each dimension :
   Int_t lNbBinsPerVar[lNbVariables] = {0};
   lNbBinsPerVar[0] = 100;
   lNbBinsPerVar[1] = 126;
   lNbBinsPerVar[2] = 24;
   lNbBinsPerVar[3] = 220;
   lNbBinsPerVar[4] = 30;
   lNbBinsPerVar[5] = 50;
  
   lNbBinsPerVar[6] = 101;
  
   lNbBinsPerVar[7] = 102;
   lNbBinsPerVar[8] = 101;
   lNbBinsPerVar[9] = 26;
   lNbBinsPerVar[10] = 26;

   lNbBinsPerVar[11] = 150; // 75 2-MeV/c2 bins
   lNbBinsPerVar[12] = 120; // 60 2-MeV/c2 bins

   lNbBinsPerVar[13] = 100;

   lNbBinsPerVar[14] = 44; // 0.05 in rapidity units
   lNbBinsPerVar[15] = 100; // 44 0.05 in rapidity units

   lNbBinsPerVar[16] = 20;

   lNbBinsPerVar[17] = 11;
   lNbBinsPerVar[18] = 100;
   lNbBinsPerVar[19] = 112; // Proper time of cascade
   lNbBinsPerVar[20] = 112; // Proper time of V0
   lNbBinsPerVar[21] = 112; // Distance V0-Xi in transverse plane

   fCFContCascadeCuts = new AliCFContainer("fCFContCascadeCuts","Container for Cascade cuts", lNbSteps, lNbVariables, lNbBinsPerVar );


   //setting the bin limits

   //0
   fCFContCascadeCuts->SetBinLimits(0, 0., 2.); // DcaXiDaughters : 0.0 to 2.0
   //1
   Double_t *lBinLim1 = new Double_t[ lNbBinsPerVar[1]+1 ];
   lBinLim1[0] = 0.0;
   lBinLim1[1] = 0.03;
   for(Int_t i=2; i< lNbBinsPerVar[1];i++) lBinLim1[i] = (Double_t)0.03 + (5. - 0.03 )/(lNbBinsPerVar[1]-2) * (Double_t)(i-1) ;
   lBinLim1[ lNbBinsPerVar[1] ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(1, lBinLim1 );
   delete [] lBinLim1; // DcaBachToPrimVertexXi : 0.0 to 0.5
   //2
   fCFContCascadeCuts->SetBinLimits(2, .9988, 1.); // XiCosineOfPointingAngle : 0.99 to 1.0
   //3
   fCFContCascadeCuts -> SetBinLimits(3, 0., 110.); // XiRadius : 0.0 to 110.0
   //4
   fCFContCascadeCuts->SetBinLimits(4, 1.1, 1.13 ); // InvMassLambdaAsCascDghter
   //5
   fCFContCascadeCuts->SetBinLimits(5, 0., 2.); // DcaV0DaughtersXi : 0.0 to 2.0
   //6
   fCFContCascadeCuts->SetBinLimits(6, .98, 1.0002); // V0CosineOfPointingAngleXi : 0.99 to 1.0
   //7
   Double_t *lBinLim7 = new Double_t[ lNbBinsPerVar[7]+1 ];
   for(Int_t i=0; i< lNbBinsPerVar[7]-1;i++) lBinLim7[i] = (Double_t)0.0 + (100. - 0.0 )/(lNbBinsPerVar[7]-2) * (Double_t)i ;
   lBinLim7[ lNbBinsPerVar[7]-1] = 200.0;
   lBinLim7[ lNbBinsPerVar[7]] = 1000.0;
   fCFContCascadeCuts -> SetBinLimits(7, lBinLim7 );
   delete [] lBinLim7; // V0RadiusXi : 0.0 to 100.0
   //8
   Double_t *lBinLim8 = new Double_t[ lNbBinsPerVar[8]+1 ];
   for(Int_t i=0; i< lNbBinsPerVar[8];i++) lBinLim8[i] = (Double_t)0.0 + (0.4 - 0.0 )/(lNbBinsPerVar[8]-1) * (Double_t)i ;
   lBinLim8[ lNbBinsPerVar[8] ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(8, lBinLim8 );
   delete [] lBinLim8; // DcaV0ToPrimVertexXi : 0.0 to 0.4
   //9
   Double_t *lBinLim9 = new Double_t[ lNbBinsPerVar[9]+1 ];
   for(Int_t i=0; i< lNbBinsPerVar[9];i++) lBinLim9[i] = (Double_t)0.0 + (0.25 - 0.0 )/(lNbBinsPerVar[9]-1) * (Double_t)i ;
   lBinLim9[ lNbBinsPerVar[9] ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(9, lBinLim9 );
   delete [] lBinLim9; // DcaPosToPrimVertexXi : 0.0 to 0.25
   //10
   Double_t *lBinLim10 = new Double_t[ lNbBinsPerVar[10]+1 ];
   for(Int_t i=0; i< lNbBinsPerVar[10];i++) lBinLim10[i] = (Double_t)0.0 + (0.25 - 0.0 )/(lNbBinsPerVar[10]-1) * (Double_t)i ;
   lBinLim10[ lNbBinsPerVar[10] ] = 100.0;
   fCFContCascadeCuts -> SetBinLimits(10, lBinLim10 );
   delete [] lBinLim10; // DcaNegToPrimVertexXi : 0.0 to 0.25
   //11
   fCFContCascadeCuts->SetBinLimits(11, 1.25, 1.40); // InvMassXi
   fCFContCascadeCuts->SetBinLimits(12, 1.62, 1.74); // InvMassOmega
   fCFContCascadeCuts->SetBinLimits(13, 0.0, 10.0); // XiTransvMom
   fCFContCascadeCuts->SetBinLimits(14, -1.1, 1.1); // Y(Xi)
   fCFContCascadeCuts->SetBinLimits(15, 0.,10.);//-1.1, 1.1); // Y(Omega)
   fCFContCascadeCuts->SetBinLimits(16, -10.0, 10.0); // BestPrimaryVtxPosZ
   Double_t *lBinLim17 = new Double_t[ lNbBinsPerVar[17]+1 ];
   for(Int_t i=3; i< lNbBinsPerVar[17]+1;i++) lBinLim17[i] = (Double_t)(i-1)*10.;
   lBinLim17[0] = 0.0;
   lBinLim17[1] = 5.0;
   lBinLim17[2] = 10.0;
   fCFContCascadeCuts -> SetBinLimits(17, lBinLim17 ); // Centrality
   delete [] lBinLim17;
   fCFContCascadeCuts->SetBinLimits(18, 0.0, 6000.0); // Track multiplicity not used now
  
   Double_t *lBinLim19 = new Double_t[ lNbBinsPerVar[19]+1 ];
   for(Int_t i=0; i< lNbBinsPerVar[19];i++) lBinLim19[i] = (Double_t)-1. + (110. + 1.0 )/(lNbBinsPerVar[19]-1) * (Double_t)i ;
   lBinLim19[ lNbBinsPerVar[19] ] = 2000.0;
   fCFContCascadeCuts->SetBinLimits(19, lBinLim19); // Proper time cascade
  
   fCFContCascadeCuts->SetBinLimits(20, lBinLim19); // Proper time V0
  
   fCFContCascadeCuts->SetBinLimits(21, lBinLim19); // Distance V0-Xi in tansverse plane
   delete [] lBinLim19;
   // Setting the number of steps : one for each cascade species (Xi-, Xi+ and Omega-, Omega+)
   fCFContCascadeCuts->SetStepTitle(0, "#Xi^{-} candidates");
   fCFContCascadeCuts->SetStepTitle(1, "#bar{#Xi}^{+} candidates");
 //  fCFContCascadeCuts->SetStepTitle(2, "#Omega^{-} candidates");
//   fCFContCascadeCuts->SetStepTitle(3, "#bar{#Omega}^{+} candidates");
  
   // Setting the variable title, per axis
   fCFContCascadeCuts->SetVarTitle(0, "DCA(#Xi daughters) (cm)");
   fCFContCascadeCuts->SetVarTitle(1, "IP(Bachelor to prim. vertex) (cm)");
   fCFContCascadeCuts->SetVarTitle(2, "cos(#Xi pointing angle)");
   fCFContCascadeCuts->SetVarTitle(3, "R_{xy}(#Xi decay) (cm)");
   fCFContCascadeCuts->SetVarTitle(4, "M_{#Lambda}(as casc. daughter) (GeV/#it{c}^{2})");
   fCFContCascadeCuts->SetVarTitle(5, "DCA(#Lambda daughters) (cm)");
  
   fCFContCascadeCuts->SetVarTitle(6, "cos(#Lambda pointing angle)");
   fCFContCascadeCuts->SetVarTitle(7, "R_{xy}(#Lambda decay) (cm)");
   fCFContCascadeCuts->SetVarTitle(8, "IP(#Lambda to prim. vertex) (cm)");
   fCFContCascadeCuts->SetVarTitle(9, "IP(Pos to prim. vertex) (cm)");
   fCFContCascadeCuts->SetVarTitle(10, "IP(Neg to prim. vertex) (cm)");
  
   fCFContCascadeCuts->SetVarTitle(11, "M_{#Lambda#pi} (GeV/#it{c}^{2})");
   fCFContCascadeCuts->SetVarTitle(12, "M_{#LambdaK} (GeV/#it{c}^{2})");
  
   fCFContCascadeCuts->SetVarTitle(13, "p_{T} (GeV/#it{c})");

   fCFContCascadeCuts->SetVarTitle(14, "y_{#Xi}");
   fCFContCascadeCuts->SetVarTitle(15, "#chi^{2} KFparticle fit");//Y(Omega)"); // not implemented 

   fCFContCascadeCuts->SetVarTitle(16, "Z-position of prim. vertex (cm)");

   fCFContCascadeCuts->SetVarTitle(17, "Centrality");

   fCFContCascadeCuts->SetVarTitle(18, "Event track multiplicity");

   fCFContCascadeCuts->SetVarTitle(19, "Proper time cascade");
   fCFContCascadeCuts->SetVarTitle(20, "Proper time #Lambda");
   fCFContCascadeCuts->SetVarTitle(21, "Distance xy #Lambda-cascade ");
  
 }
  

  fHistpXiSignalRealKstar = new TH2F("fHistpXiSignalRealKstar","",500,0.,5.,20,0.,100.);
  fOutputContainer->Add(fHistpXiSignalRealKstar);
  fHistapXiSignalRealKstar = new TH2F("fHistapXiSignalRealKstar","",500,0.,5.,20,0.,100.);
  fOutputContainer->Add(fHistapXiSignalRealKstar);
  fHistpaXiSignalRealKstar = new TH2F("fHistpaXiSignalRealKstar","",500,0.,5.,20,0.,100.);
  fOutputContainer->Add(fHistpaXiSignalRealKstar);
  fHistapaXiSignalRealKstar = new TH2F("fHistapaXiSignalRealKstar","",500,0.,5.,20,0.,100.);
  fOutputContainer->Add(fHistapaXiSignalRealKstar);

  fHistpXiSignalBkgKstar = new TH2F("fHistpXiSignalBkgKstar","",500,0.,5.,20,0.,100.);
  fOutputContainer->Add(fHistpXiSignalBkgKstar);
  fHistapXiSignalBkgKstar = new TH2F("fHistapXiSignalBkgKstar","",500,0.,5.,20,0.,100.);
  fOutputContainer->Add(fHistapXiSignalBkgKstar);
  fHistpaXiSignalBkgKstar = new TH2F("fHistpaXiSignalBkgKstar","",500,0.,5.,20,0.,100.);
  fOutputContainer->Add(fHistpaXiSignalBkgKstar);
  fHistapaXiSignalBkgKstar = new TH2F("fHistapaXiSignalBkgKstar","",500,0.,5.,20,0.,100.);
  fOutputContainer->Add(fHistapaXiSignalBkgKstar);

  fHistFractionOfXiWithSharedDaughters = new TH2F("fHistFractionOfXiWithSharedDaughters","",200,0.,200.,201,0.,1.005);
  fOutputContainer->Add(fHistFractionOfXiWithSharedDaughters);
  fHistFractionOfaXiWithSharedDaughters = new TH2F("fHistFractionOfaXiWithSharedDaughters","",200,0.,200.,201,0.,1.005);
  fOutputContainer->Add(fHistFractionOfaXiWithSharedDaughters);

  fHistpXibacDetaSDphiS = new TH2F("fHistpXibacDetaSDphiS","",1000,-1,1,200,-2.,2.);    
  fOutputContainer->Add(fHistpXibacDetaSDphiS);
  fHistpXiposDetaSDphiS = new TH2F("fHistpXiposDetaSDphiS","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistpXiposDetaSDphiS);
  fHistpXinegDetaSDphiS = new TH2F("fHistpXinegDetaSDphiS","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistpXinegDetaSDphiS);
  fHistpaXibacDetaSDphiS = new TH2F("fHistpaXibacDetaSDphiS","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistpaXibacDetaSDphiS);
  fHistpaXiposDetaSDphiS = new TH2F("fHistpaXiposDetaSDphiS","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistpaXiposDetaSDphiS);
  fHistpaXinegDetaSDphiS = new TH2F("fHistpaXinegDetaSDphiS","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistpaXinegDetaSDphiS);
  fHistapXibacDetaSDphiS = new TH2F("fHistapXibacDetaSDphiS","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistapXibacDetaSDphiS);
  fHistapXiposDetaSDphiS = new TH2F("fHistapXiposDetaSDphiS","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistapXiposDetaSDphiS);
  fHistapXinegDetaSDphiS = new TH2F("fHistapXinegDetaSDphiS","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistapXinegDetaSDphiS);
  fHistapaXibacDetaSDphiS = new TH2F("fHistapaXibacDetaSDphiS","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistapaXibacDetaSDphiS);
  fHistapaXiposDetaSDphiS = new TH2F("fHistapaXiposDetaSDphiS","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistapaXiposDetaSDphiS);
  fHistapaXinegDetaSDphiS = new TH2F("fHistapaXinegDetaSDphiS","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistapaXinegDetaSDphiS);

  fHistpXibacDetaSDphiSBkg = new TH2F("fHistpXibacDetaSDphiSBkg","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistpXibacDetaSDphiSBkg);
  fHistpXiposDetaSDphiSBkg = new TH2F("fHistpXiposDetaSDphiSBkg","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistpXiposDetaSDphiSBkg);
  fHistpXinegDetaSDphiSBkg = new TH2F("fHistpXinegDetaSDphiSBkg","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistpXinegDetaSDphiSBkg);
  fHistpaXibacDetaSDphiSBkg = new TH2F("fHistpaXibacDetaSDphiSBkg","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistpaXibacDetaSDphiSBkg);
  fHistpaXiposDetaSDphiSBkg = new TH2F("fHistpaXiposDetaSDphiSBkg","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistpaXiposDetaSDphiSBkg);
  fHistpaXinegDetaSDphiSBkg = new TH2F("fHistpaXinegDetaSDphiSBkg","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistpaXinegDetaSDphiSBkg);
  fHistapXibacDetaSDphiSBkg = new TH2F("fHistapXibacDetaSDphiSBkg","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistapXibacDetaSDphiSBkg);
  fHistapXiposDetaSDphiSBkg = new TH2F("fHistapXiposDetaSDphiSBkg","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistapXiposDetaSDphiSBkg);
  fHistapXinegDetaSDphiSBkg = new TH2F("fHistapXinegDetaSDphiSBkg","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistapXinegDetaSDphiSBkg);
  fHistapaXibacDetaSDphiSBkg = new TH2F("fHistapaXibacDetaSDphiSBkg","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistapaXibacDetaSDphiSBkg);
  fHistapaXiposDetaSDphiSBkg = new TH2F("fHistapaXiposDetaSDphiSBkg","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistapaXiposDetaSDphiSBkg);
  fHistapaXinegDetaSDphiSBkg = new TH2F("fHistapaXinegDetaSDphiSBkg","",1000,-1,1,200,-2.,2.);
  fOutputContainer->Add(fHistapaXinegDetaSDphiSBkg);

  fHistTrackBufferOverflow = new TH1F("fHistTrackBufferOverflow","",2,0,2);
  fOutputContainer->Add(fHistTrackBufferOverflow);

  PostData(1, fOutputContainer );
  PostData(2, fCFContCascadeCuts);

}

//---------------------------------------------------------------------------------

void AliAnalysisTaskhCascadeFemto::UserExec(Option_t *) {
  // Main loop
  // Called for each event

  if(!fPIDResponse) {
    AliError("Cannot get pid response");
    return;
  }

  AliVVertex *vertexmain =0x0;
  AliCentrality* centrality = 0x0;
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};


  AliMCEvent  *lMCevent  = 0x0;
  AliStack    *lMCstack  = 0x0;
  TClonesArray *arrayMC = 0x0;

  Int_t ncascades = 0;
  Int_t ntracks = 0;

  if (fAnalysisType == "ESD") {
    fESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!fESDevent) {
      AliWarning("ERROR: ESDevent not available \n");
      return;
    }

    if (fReadMCTruth) {
      lMCevent = MCEvent();
      if (!lMCevent) {
        Printf("ERROR: Could not retrieve MC event \n");
        //cout << "Name of the file with pb :" <<  CurrentFileName() << endl;
        return;
      }
      lMCstack = lMCevent->Stack();
      if (!lMCstack) {
        Printf("ERROR: Could not retrieve MC stack \n");
        //cout << "Name of the file with pb :" <<  CurrentFileName() << endl;
        return;
      }
    }

    if((TMath::Abs(lBestPrimaryVtxPos[2])) > 10.) return;
    ncascades = fESDevent->GetNumberOfCascades();
    ntracks = fESDevent->GetNumberOfTracks(); 
    centrality = fESDevent->GetCentrality();

    const AliESDVertex *lPrimaryBestESDVtx = fESDevent->GetPrimaryVertex();
    if (!lPrimaryBestESDVtx){
      AliWarning("No prim. vertex in ESD... return!");
      PostData(1, fOutputContainer);
      PostData(2, fCFContCascadeCuts);

      return;
    }
    vertexmain = (AliVVertex*) lPrimaryBestESDVtx;
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );

  } else if (fAnalysisType == "AOD") {
    fAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
    if (!fAODevent) {
      AliWarning("ERROR: AODevent not available \n");
      return;
    }

    ncascades  = fAODevent->GetNumberOfCascades();
    ntracks = fAODevent->GetNumberOfTracks();
    centrality = fAODevent->GetCentrality();

    const AliAODVertex *lPrimaryBestAODVtx = fAODevent->GetPrimaryVertex();
    if (!lPrimaryBestAODVtx){
      AliWarning("No prim. vertex in AOD... return!");
      PostData(1, fOutputContainer);
      PostData(2, fCFContCascadeCuts);

      return;
    }
    vertexmain = (AliVVertex*) lPrimaryBestAODVtx;
    lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );

    if (fReadMCTruth) {
//      Printf("Reading MC truth!!! \n");
      arrayMC = (TClonesArray*) fAODevent->GetList()->FindObject(AliAODMCParticle::StdBranchName());

      if (!arrayMC) AliFatal("Error: MC particles branch not found!\n");

    }

  } else {

    Printf("Analysis type (ESD or AOD) not specified \n");
    return;

  }

  Float_t lcentrality = centrality->GetCentralityPercentile("V0M");
  if (lcentrality<fCentrLowLim||lcentrality>=fCentrUpLim)  return;
  if((TMath::Abs(lBestPrimaryVtxPos[2])) > 10.) return;

  fHistCentrality->Fill(lcentrality);
  fHistVertexDistribution->Fill(lBestPrimaryVtxPos[2]); 

  const Float_t bfield = (InputEvent())->GetMagneticField();
  int fieldsign;
  if (bfield >=0.) fieldsign = 1;
  else fieldsign = -1;
  // Store the event in the buffer to do mixing
  // ... find vertex... 
  int zBin=0;

  double zStep=2*10/double(fzVertexBins), zStart=-10.;

  for (int i=0; i<fzVertexBins; i++) {
    if ((lBestPrimaryVtxPos[2] > zStart+i*zStep) && (lBestPrimaryVtxPos[2] < zStart+(i+1)*zStep)) {
      zBin=i;
      break;
    }
  }

  // ... and centrality
  int centralityBin=0;

  if(lcentrality < 5.) centralityBin=19;  // changed <= with < to be consistent with histogram binning, except last bin 
  else if(lcentrality < 10.) centralityBin=18;
  else if(lcentrality < 15.) centralityBin=17;
  else if(lcentrality < 20.) centralityBin=16;
  else if(lcentrality < 25.) centralityBin=15;
  else if(lcentrality < 30.) centralityBin=14;
  else if(lcentrality < 35.) centralityBin=13;
  else if(lcentrality < 40.) centralityBin=12;
  else if(lcentrality < 45.) centralityBin=11;
  else if(lcentrality < 50.) centralityBin=10;
  else if(lcentrality < 55.) centralityBin=9; 
  else if(lcentrality < 60.) centralityBin=8;
  else if(lcentrality < 65.) centralityBin=7;
  else if(lcentrality < 70.) centralityBin=6;
  else if(lcentrality < 75.) centralityBin=5;
  else if(lcentrality < 80.) centralityBin=4;
  else if(lcentrality < 85.) centralityBin=3;
  else if(lcentrality < 90.) centralityBin=2;   // from here on, wont be filled because the range selected in AddTask is 0-90
  else if(lcentrality < 95.) centralityBin=1;
  else if(lcentrality <= 100.) centralityBin=0;

  fEventColl[zBin][centralityBin]->FifoShift();

  fEvt = fEventColl[zBin][centralityBin]->fEvt;

//  printf("buffer size: %d\n",fTrackBufferSize);
//  printf("ntracks: %d\n",ntracks);

  for (Int_t igt = 0; igt < fTrackBufferSize; igt++) farrGT[igt] = -1;

  AliAODTrack *globaltrack = 0x0;
 
  // Read and store global tracks to retrieve PID information for TPC only tracks
  for (Int_t igt = 0; igt < ntracks; igt++) {
    globaltrack = (AliAODTrack*) fAODevent->GetTrack(igt);

    if (!globaltrack) continue; 
    if (globaltrack->GetID()<0 ) continue;
    if (!globaltrack->IsOn(AliAODTrack::kTPCrefit)) continue; // there are such tracks with no TPC clusters

    // Check id is not too big for buffer
    if (globaltrack->GetID()>=fTrackBufferSize) {
      printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n",globaltrack->GetID(),fTrackBufferSize);
      fHistTrackBufferOverflow->Fill(1);
      continue;
    }

//    if ( !(globaltrack->GetFilterMap()) ) { 
//        cout<<" No filter map for this global track!!  "<<globaltrack->GetFilterMap()<<endl;
//        continue;
//    }

    if ( globaltrack->GetTPCNcls()<=0   ) { // such tracks do not have the TPC refit either, filter map is 2 --> ITS constrained
//      cout<<" No TPC cl for this global track!!  "<<igt<<endl;
//      if (!globaltrack->IsOn(AliAODTrack::kTPCrefit)) cout<<" ... and also no tpc refit  "<<globaltrack->GetFilterMap()<<endl;
        continue;
    }
    //cout<<" The array will contain "<<igt<<" , it contains "<<farrGT[globaltrack->GetID()]<<endl; 

    // Warn if we overwrite a track
    if (farrGT[globaltrack->GetID()]>=0) { // Two tracks same id  --> checked that it never happens
//      cout<<" Array already filled "<<farrGT[globaltrack->GetID()]<<endl;
    } else { 
      farrGT[globaltrack->GetID()] = igt;           // solution adopted in the femto framework
//    cout<<" Set array, now it contains  "<<farrGT[globaltrack->GetID()]<<endl;
    }
  }
 
  globaltrack = 0x0; 

  // Read, plot and store pions or protons
  
  AliAODTrack *track = 0x0;
  AliVTrack *vtrackg = 0x0;
  AliVTrack *vtrack = 0x0;

  Bool_t isTOFPIDok = kFALSE;
  Float_t nsigmaTOFp = 0.;
  Float_t nsigmaTPCp = 0.;
  Float_t probMis = 0.;
  Double32_t tTOF = 0.;
  Double_t beta = 0.; 
  Double_t length = 0.;

  Float_t nTPCCrossedRows = 0.; 
  Float_t rationCrnFind = 0.;

  UShort_t nTPCclSp = 0;
  UShort_t nTPCclp = 0;

  Float_t sharedFractionTPCcl = 0.;
  AliPIDResponse::EDetPidStatus statusTOF;
  AliPIDResponse::EDetPidStatus statusTPC;
 
  Double_t expectedTimes[AliPID::kSPECIES];

  Double_t dz[2]= {-999.,-999.}; //Double_t covar[3]={-999.,-999.,-999.};
  Double_t dzg[2]= {-999.,-999.}; Double_t covarg[3]={-999.,-999.,-999.};
  AliExternalTrackParam etp1;


  Float_t rapidity = 0.;
  Short_t charge = -2;


  int pCount = 0;

  Float_t pMomentumTruth[3];
  AliReconstructedProton::MCProtonOrigin_t mcProtonOrigin = AliReconstructedProton::kUnassigned;

  Bool_t isP = kFALSE;  // particle
  Bool_t isaP = kFALSE; // anti-particle 

  for (Int_t ip = 0; ip < ntracks; ip++) {

    vtrack = fAODevent->GetTrack(ip);
    if (!vtrack) continue;
    track = dynamic_cast<AliAODTrack*>(vtrack);
    if(!track) AliFatal("Not a standard AOD");

    // Take TPC only tracks constrained to the SPD PV during filtering, Filter Bit 7 (128) but PID is not stored --> find the corresponding global track filter bit 1
    if(!track->TestFilterBit(AliAODTrack::kTrkTPCOnlyConstrained)) continue; // Bit 7 corresponds to the following cuts (applied to original ESD track, see Ruben presentation 13.02.2014) : 
                                                                             // ncl 50  
                                                                             // dcaxy 2.4 dcaz 3.2 // this cut however is on the global track!! PWGCF meeting (17.02.2014)
                                                                             // no kinks 
                                                                             // chi2percl 4 
                                                                             // dcatovtx2D kTRUE

    // nTPCCrossedRows = track->GetTPCClusterInfo(2, 1); // The method below is equivalent
    nTPCCrossedRows = track->GetTPCNCrossedRows();
//    cout<<"Crossed rows dir "<<track->GetTPCNCrossedRows()<<" Crossed rows indir "<<nTPCCrossedRows<<endl;
//    if(nTPCCrossedRows<70) continue;  

    rationCrnFind = 0.;
    if(track->GetTPCNclsF()) rationCrnFind = nTPCCrossedRows/track->GetTPCNclsF();   
    if (fkApplyRatioCrRnFindCut) if(rationCrnFind< .8) continue;  

    nTPCclSp =  track->GetTPCnclsS();
    nTPCclp = track->GetTPCncls();
    if (nTPCclp<70) continue;

    sharedFractionTPCcl = (Float_t) nTPCclSp/nTPCclp;
    if (sharedFractionTPCcl>0.) continue; // for protons it is not more than 3% (without dca cuts) 
//    cout<<"Shared fraction ncl "<<sharedFractionTPCcl<<endl;

    // Get the corresponding global track to use PID --> stored only for global tracks
    // 
    // Check that the array fGTI isn't too small
    // for the track id
    if (-track->GetID()-1 >= fTrackBufferSize) {
      printf ("Exceeding buffer size!!");
      continue;
    }

    vtrackg = fAODevent->GetTrack(farrGT[-track->GetID()-1]);
    if (!vtrackg) {
      printf ("No global info! iTrack %d, ID %d\n",ip,track->GetID());
      continue;
    }
    if (farrGT[-track->GetID()-1]>=ntracks || farrGT[-track->GetID()-1]<0) { /*cout<<"This index is out of range!!"<<farrGT[-track->GetID()-1]<<endl;*/ continue;}
    globaltrack = dynamic_cast<AliAODTrack*>(vtrackg); 
    if(!globaltrack) AliFatal("Not a standard AOD");

//    cout<<" Filter map for the global track "<<globaltrack->GetFilterMap()<<endl;
 

    // IP to PV of tracks
    dz[0]= -999.; dz[1] = -999.;
    dzg[0]= -999.; dzg[1] = -999.;

    dz[0] = track->DCA();    // the TPC one should be applied the other biases the CF --> from Maciejs note --> FIXME to be checked 
    dz[1] = track->ZAtDCA(); // for those two lines check AliAODTrack.h // FIXME these two lines produce shifted distributions, known problem, asked Mac and Marian. 
                                                                        // Btw dont propagate TPC constrained! 
//    Double_t p[3]; track->GetXYZ(p); // first two elements of p give the same as above, ignore the third, those are original DCA of TPC only tracks (with or w/o constraint?)
//    cout<<"d xy "<<dz[0]<<"while value is for other "<<p[0]<<endl;
//    cout<<"d z "<<dz[1]<<"while value is for other "<<p[1]<<endl;

    etp1.CopyFromVTrack(vtrackg); 
    etp1.PropagateToDCA(vertexmain,(InputEvent())->GetMagneticField(),100.,dzg,covarg);    

    isP = kFALSE;
    isaP = kFALSE;

    charge = globaltrack->Charge();
 
    if (TMath::Abs(track->Eta())> 0.8) continue;
    if (track->Chi2perNDF()>4.) continue; // should be redundant already applied in the TPC only track ESD filtering
    if (track->Pt()<fMinPtForPrim || track->Pt()> fMaxPtForPrim) continue;  

    if (fkCutOnTPCIP) {
      if (TMath::Abs(dz[0])>fIPCutxy ) continue;  // 2.4 proton 1. pion 
      if (TMath::Abs(dz[1])>fIPCutz ) continue;  // 3.2  proton 1. pion
    } else {
      if (TMath::Abs(dzg[0])>fIPCutxy ) continue;  // 0.1 proton/pion
      if (TMath::Abs(dzg[1])>fIPCutz ) continue; // 0.15 proton/pion
    }


    if (fFirstpart == kPion) { 
      if (charge>0) isaP = kTRUE;
      else if (charge<0) isP = kTRUE;
    } else if (fFirstpart == kProton) {
      if (charge>0) isP = kTRUE;
      else if (charge<0) isaP = kTRUE;
    }

    rapidity = 0.5*TMath::Log( (track->E(fPDGfirst) + track->Pz()) / (track->E(fPDGfirst) - track->Pz() +1.e-13) );

//    cout<<" Filter map  "<<track->GetFilterMap()<<endl;


    nsigmaTOFp = 0.; // be careful with those initialization
    nsigmaTPCp = 0.;
    probMis = 0.;

    statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,globaltrack);

    if (AliPIDResponse::kDetPidOk != statusTPC) { 
      //printf ("TPC PID not there"); 
      continue; 
    }

    if (fFirstpart == kPion) nsigmaTPCp = fPIDResponse->NumberOfSigmasTPC(globaltrack, AliPID::kPion); 
    else if (fFirstpart == kProton) nsigmaTPCp = fPIDResponse->NumberOfSigmasTPC(globaltrack, AliPID::kProton); 
   
    if (fFirstpart == kProton) { // For protons use TPC and TOF

      statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,globaltrack); // this checks kTOFout and kTIMEi https://twiki.cern.ch/twiki/bin/viewauth/ALICE/TOF

      tTOF = 0.0;
      beta = -1000.;
      length = 0.;
      isTOFPIDok = kFALSE;

      if ( (statusTOF ==  AliPIDResponse::kDetPidOk) ) {  
        isTOFPIDok = kTRUE;
        probMis = fPIDResponse->GetTOFMismatchProbability(globaltrack);
        nsigmaTOFp = fPIDResponse->NumberOfSigmasTOF(globaltrack, AliPID::kProton); 
        tTOF = globaltrack->GetTOFsignal(); 
        tTOF -= fPIDResponse->GetTOFResponse().GetStartTime(globaltrack->P());      
        globaltrack->GetIntegratedTimes(expectedTimes);
        tTOF -= expectedTimes[AliPID::kProton];  

        // length = globaltrack->GetIntegratedLength();  // // FIXME length is zero!! from a mail february 2014: this info is not available for AODs 115, use AODs 145
        // if (tTOF > 0.) {//&&length>350) {
        // beta = length/(tTOF*2.99792457999999984e-02); 
        // cout<<" rack length  "<<length<<" beta "<<beta<<endl; 
        // gamma = 1/TMath::Sqrt(1 - beta*beta);
        // mass = ptot/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx. 
        // }
      } else { 
        probMis = 0.; nsigmaTOFp = 0.; //cout<<"The corresponding global track has no tof pid!"<<endl;
      } 
     
    
      // PID different momentum intervals Hans 0-0.75-1.-5 TPC-TPCTOF-TOF(+TPC to detect mismatch if present)
                                     // Maciej 0.8       TPC-TPC+TOF
      if (track->P()< fMomemtumLimitForTOFPID) {
        if (TMath::Abs(nsigmaTPCp)> fnSigmaTPCPIDfirstParticle) continue;  
      } else {
        if (isTOFPIDok) {
          if (probMis > 0.01) continue;
          if (TMath::Sqrt(nsigmaTOFp*nsigmaTOFp+nsigmaTPCp*nsigmaTPCp)> fnSigmaTPCTOFPIDfirstParticle) continue;   // this cleans the TOF corrected time plot vs p       
        } else { 
          continue; // if no TOF we dont use this track // NB stat is 4 times smaller 
        //if (TMath::Abs(nsigmaTPCp)> 3.) continue; 
        }
      }
    } else if (fFirstpart == kPion) { // for pions only TPC
      if (TMath::Abs(nsigmaTPCp)> fnSigmaTPCPIDfirstParticle) continue; 
    }

    fHistyptProtons->Fill(track->Pt(),rapidity,lcentrality);
    fHistphietaProtons->Fill(track->Phi(),track->Eta());
    if (track->TestBit(AliAODTrack::kIsDCA)) { fHistIPtoPVxyzTPC->Fill(dz[0],dz[1]); } else { /*cout<<"Bit is not DCA!!!!!!!!!!!"<<endl;*/ }
    fHistIPtoPVxyzGlobal->Fill(dzg[0],dzg[1]);
   
//       cout<<"TOF signal "<<globaltrack->GetTOFsignal()<<endl; 
    if (isTOFPIDok) {
      fHistpTOFmisvspt->Fill(track->Pt(),probMis);
      fHistpTOFmisvsp->Fill(track->P(),probMis);
      fHistpTOFnsigmavspt->Fill(track->Pt(),nsigmaTOFp);
      fHistpTOFnsigmavsp->Fill(track->P(),nsigmaTPCp);
      fHistpTOFsignalvsp->Fill(track->P(),tTOF);  
      fHistpTOFsignalvspt->Fill(track->Pt(),tTOF);

      fHistpTOFTPCsignalvspt->Fill(track->Pt(),TMath::Sqrt(nsigmaTOFp*nsigmaTOFp+nsigmaTPCp*nsigmaTPCp));
    }

    fHistprimpTPCdEdx->Fill(globaltrack->GetTPCmomentum()*charge, globaltrack->GetTPCsignal());
    
   
    fHistnTPCCrossedR->Fill(nTPCCrossedRows);
    fHistRationTPCCrossedRnFind->Fill(rationCrnFind);
    fHistSharedFrTPCcl->Fill(sharedFractionTPCcl);

    // if(fReadMCTruth)  ProtonOrigin(Float_t dcagtxy, AliVTack *track, arrayMC, Float_t pt); //  for purity studies  

    //Save first particle information
    fEvt->fReconstructedProton[pCount].pCharge = charge;
    //cout<<"Charge of pion "<<fEvt->fReconstructedProton[pCount].pCharge<<endl;
    
    fEvt->fReconstructedProton[pCount].pMomentum[0]  = track->Px();
    fEvt->fReconstructedProton[pCount].pMomentum[1]  = track->Py();
    fEvt->fReconstructedProton[pCount].pMomentum[2]  = track->Pz();

    for (int pIndex = 0; pIndex <3; pIndex++) {
      fEvt->fReconstructedProton[pCount].pMomentumTruth[pIndex]  = pMomentumTruth[pIndex];
    }
    fEvt->fReconstructedProton[pCount].pPt     = track->Pt();
    fEvt->fReconstructedProton[pCount].pEta    = track->Eta();
    fEvt->fReconstructedProton[pCount].pPhi    = track->Phi();
    fEvt->fReconstructedProton[pCount].pRap    = rapidity; 
    fEvt->fReconstructedProton[pCount].mcProtonOriginType = mcProtonOrigin;

    if (isP) {fEvt->fReconstructedProton[pCount].isP = kTRUE; fEvt->fReconstructedProton[pCount].isaP = kFALSE;}
    else if (isaP) {fEvt->fReconstructedProton[pCount].isaP = kTRUE; fEvt->fReconstructedProton[pCount].isP = kFALSE;}

    fEvt->fReconstructedProton[pCount].index = globaltrack->GetID();
//    cout<<" Pos 125 cm before set "<<fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition[0]<<" "<<fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition[1]<<" "<<fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition[2]<<endl;
    if (fkPropagateGlobal) SetSftPosR125(vtrackg, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition); // Hans popagates the TPC tracks, I try both, flag to choose 
    else SetSftPosR125(vtrack, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition);
//    cout<<" Pos set 125 cm "<<fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition[0]<<" "<<fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition[1]<<" "<<fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition[2]<<endl; 
    fEvt->fReconstructedProton[pCount].pEtaS = EtaS(fEvt->fReconstructedProton[pCount].pShiftedGlobalPosition); 
//    cout<<" Etas "<<fEvt->fReconstructedProton[pCount].pEtaS<<endl;

//    cout<<" This proton in the stack with pt "<<fEvt->fReconstructedProton[pCount].pPt<<endl;

    pCount++;
    if (fMaxPMult <= pCount){
      cerr<<"Proton counts exceeded "<<fMaxPMult<<"!"<<endl;
      break;
    }

  }

  
  fEvt->fNumberCandidateProton = pCount;
  fHistProtonMultvsCent->Fill(pCount,lcentrality);

  if (fSecondpart == kXi || fSecondpart == kOmega) {

  // Read and plot Xi both ESDs and AODs 
    Double_t lInvMassXiMinus = 0.;
    Double_t lInvMassXiPlus = 0.;
    Double_t lInvMassOmegaMinus = 0.;
    Double_t lInvMassOmegaPlus = 0.;

    Double_t lInvMassXi = 0.;
    Double_t lInvMassLambdaAsCascDghter = 0.;
    Double_t lInvMassAntiLambdaAsCascDghter = 0.;
    Double_t lInvMassLambda = 0.;

    Int_t lPosTPCClusters = -1; 
    Int_t lNegTPCClusters = -1;  
    Int_t lBachTPCClusters = -1; 

    Int_t lPosTPCClustersS;
    Int_t lNegTPCClustersS;
    Int_t lBachTPCClustersS;

    Int_t sharedFractionTPCclPos;
    Int_t sharedFractionTPCclNeg;
    Int_t sharedFractionTPCclBach;

    AliVTrack *pTrackXi = 0x0;
    AliVTrack *nTrackXi = 0x0;
    AliVTrack *bachTrackXi = 0x0;

    Bool_t lIsBachelorPionForTPC = kFALSE;
    Bool_t lIsBachelorKaonForTPC = kFALSE;
    Bool_t lIsNegPionForTPC      = kFALSE;
    Bool_t lIsNegProtonForTPC    = kFALSE;
    Bool_t lIsPosPionForTPC      = kFALSE;
    Bool_t lIsPosProtonForTPC    = kFALSE;

    Double_t lXiMomX       = 0., lXiMomY = 0., lXiMomZ = 0.;
    Double_t lXiTransvMom  = 0., lXiMom  = 0. ;
    Double_t lV0PMomX       = 0., lV0PMomY = 0., lV0PMomZ = 0.;
    Double_t lV0NMomX       = 0., lV0NMomY = 0., lV0NMomZ = 0.;
    Double_t lV0TransvMom     = 0., lV0Mom   = 0. ; 

    Double_t lV0MomX = 0;
    Double_t lV0MomY = 0;
    Double_t lV0MomZ = 0;

    Short_t lChargeXi = -2;

    AliESDcascade *xi = 0x0;
    UInt_t lIdxPosXi  = 0;
    UInt_t lIdxNegXi  = 0;
    UInt_t lBachIdx   = 0;
    Double_t lV0quality = 0;

    ULong_t pStatus    = 0;
    ULong_t nStatus    = 0;
    ULong_t bachStatus = 0;

    const AliAODcascade *xiaod = 0x0;
    AliAODTrack* aodpTrackXi = 0x0;
    AliAODTrack* aodnTrackXi = 0x0;
    AliAODTrack* aodbachTrackXi = 0x0;

    Float_t lRapXi;  
    Float_t lEtaXi;
    Float_t lPhiXi;

    Float_t etaPos;
    Float_t etaNeg;
    Float_t etaBach;

    Float_t xiMomentumTruth[3];
    AliReconstructedXi::MCXiOrigin_t mcXiOrigin = AliReconstructedXi::kUnassigned;
    int xiCount = 0;

    Bool_t isXi = kFALSE;
    Bool_t isaXi = kFALSE;
    Bool_t isOmega = kFALSE;
    Bool_t isaOmega = kFALSE;

    Int_t indexB, indexP, indexN;

    // Topological variables
    Double_t lContainerCutVars[22] = {0.};

    Double_t lChi2Xi = 0. ; // not implemented 
    Double_t lDcaXiDaughters = -1. ;  
    Double_t lXiCosineOfPointingAngle = -1. ; 
    Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 }; 
    Double_t lXiRadius = -1000. ; 
    Float_t lctau = -1.;
    Float_t lctauV0 = -1.;
    Float_t distV0Xi = -1.;
    Float_t ldistTV0Xi = -1.;

    // V0 part in cascades
    Double_t lDcaV0DaughtersXi = -1.; 
    Double_t lV0RadiusXi = -1.;  
    Double_t lPosV0Xi[3] = {-1.,-1.,-1.}; 

    Double_t lDcaBachToPrimVertexXi = -1., lDcaV0ToPrimVertexXi = -1.; 
    Double_t lDcaPosToPrimVertexXi = -1.; 
    Double_t lDcaNegToPrimVertexXi = -1.; 
    Double_t lV0CosineOfPointingAngle = -1. ; 
    Double_t lV0toXiCosineOfPointingAngle = 0. ; // not used

    for (Int_t iXi = 0; iXi < ncascades; iXi++) {

      lInvMassXiMinus = 0.;
      lInvMassXiPlus = 0.;
      lInvMassOmegaMinus = 0.;
      lInvMassOmegaPlus = 0.; 
      lInvMassXi = 0.;
      lInvMassLambdaAsCascDghter = 0.;
      lInvMassAntiLambdaAsCascDghter = 0.;
      lPosTPCClusters = 0;
      lNegTPCClusters = 0;
      lBachTPCClusters = 0;

      pTrackXi = 0x0;
      nTrackXi = 0x0;
      bachTrackXi = 0x0;

      lIsBachelorPionForTPC = kFALSE;
      lIsBachelorKaonForTPC = kFALSE;
      lIsNegPionForTPC      = kFALSE;
      lIsNegProtonForTPC    = kFALSE;
      lIsPosPionForTPC      = kFALSE;
      lIsPosProtonForTPC    = kFALSE;

      lXiMomX       = 0.; lXiMomY  = 0.; lXiMomZ  = 0.;
      lXiTransvMom  = 0.;
      lV0PMomX      = 0.; lV0PMomY = 0.; lV0PMomZ = 0.;
      lV0NMomX      = 0.; lV0NMomY = 0.; lV0NMomZ = 0.;
      lV0TransvMom  = 0.; lV0Mom   = 0. ;


      lChargeXi = -2;

      isXi = kFALSE;
      isaXi = kFALSE;
      isOmega = kFALSE;
      isaOmega = kFALSE;

     
      if (fAnalysisType == "ESD") { // %| not mantained 

        xi = fESDevent->GetCascade(iXi);
        if (!xi) continue;

        lIdxPosXi  = (UInt_t) TMath::Abs( xi->GetPindex() );
        lIdxNegXi  = (UInt_t) TMath::Abs( xi->GetNindex() );
        lBachIdx   = (UInt_t) TMath::Abs( xi->GetBindex() );
        if (lBachIdx == lIdxNegXi) { // prevented already in cascadevertexer
          AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
        }
        if (lBachIdx == lIdxPosXi) {
          AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
        }

        pTrackXi             = (AliVTrack*) fESDevent->GetTrack( lIdxPosXi ); // FIXME check if conversion to VTrack works
        nTrackXi             = (AliVTrack*) fESDevent->GetTrack( lIdxNegXi );
        bachTrackXi          = (AliVTrack*) fESDevent->GetTrack( lBachIdx );
        if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
          AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
          continue;
        }

        lV0quality = 0;
        // Charge
        if ( bachTrackXi->Charge() < 0 )  {  
          lV0quality = 0.;
          xi->ChangeMassHypothesis(lV0quality , 3312);   // Xi- 
          lInvMassXiMinus = xi->GetEffMassXi();
          lInvMassLambdaAsCascDghter = xi->GetEffMass(); // This method without parameters returns fEffMass (set in CascadeVertexer)
        } else if ( bachTrackXi->Charge() >  0 ) {
          lV0quality = 0.;
          xi->ChangeMassHypothesis(lV0quality , -3312);  // Xi+, no need to go back to def hyp 
          lInvMassXiPlus = xi->GetEffMassXi();
          lInvMassAntiLambdaAsCascDghter = xi->GetEffMass(); // This method without parameters returns fEffMass (set in CascadeVertexer) 
        }
        lChargeXi = xi->Charge();

//      if (fkQualityCutTPCrefit) {
                 // 1 - Poor quality related to TPCrefit
        pStatus    = pTrackXi->GetStatus();
        nStatus    = nTrackXi->GetStatus();
        bachStatus = bachTrackXi->GetStatus();

//        if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; }i // Redundant required in V0f
//        if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; } // idem
        if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; }
//      }

        xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );

        xi->GetNPxPyPz(lV0NMomX,lV0NMomY,lV0NMomZ);
        xi->GetPPxPyPz(lV0PMomX,lV0PMomY,lV0PMomZ);

        lV0TransvMom = TMath::Sqrt(TMath::Power(lV0NMomX+lV0PMomX,2)+TMath::Power(lV0NMomY+lV0PMomY,2));


      } else if (fAnalysisType == "AOD") {

        xiaod = fAODevent->GetCascade(iXi);
        if (!xiaod) continue;

        pTrackXi    = dynamic_cast<AliVTrack*>( xiaod->GetDaughter(0) );
        nTrackXi    = dynamic_cast<AliVTrack*>( xiaod->GetDaughter(1) );
        bachTrackXi = dynamic_cast<AliVTrack*>( xiaod->GetDecayVertexXi()->GetDaughter(0) );
        if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
          AliWarning("ERROR: Could not retrieve one of the 3 AOD daughter tracks of the cascade ...");
          continue;
        }

//      lIdxPosXi  = (UInt_t) TMath::Abs( aodpTrackXi->GetID() );  
//      lIdxNegXi  = (UInt_t) TMath::Abs( aodnTrackXi->GetID() );
//      lBachIdx   = (UInt_t) TMath::Abs( aodbachTrackXi->GetID() );

//      if (lBachIdx == lIdxNegXi) {  // Redundand, already checked in Cascade vertexer 
//        AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
//      }
//      if (lBachIdx == lIdxPosXi) {
//        AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
//      }

//      if (fkQualityCutTPCrefit) { 

//        if (!(aodpTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); continue; } // already in vertexer
//        if (!(aodnTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); continue; } // idem
//        if (!(aodbachTrackXi->IsOn(AliAODTrack::kTPCrefit))) { AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); continue; } // TPC PID will remove this anyhow

//      }

     

        aodpTrackXi = dynamic_cast<AliAODTrack*>(pTrackXi); // to convert aod in vtracks
        aodnTrackXi = dynamic_cast<AliAODTrack*>(nTrackXi);
        aodbachTrackXi = dynamic_cast<AliAODTrack*>(bachTrackXi);

        //if (aodpTrackXi->Charge()<0) cout<<"Wrong Sign for pos!!!!!!!!!!!!"<<endl;
        //if (aodnTrackXi->Charge()>0) cout<<"Wrong Sign for neg!!!!!!!!!!!!"<<endl;

        lChargeXi = xiaod->ChargeXi();
        if ( lChargeXi < 0 )  { lInvMassXiMinus = xiaod->MassXi(); lInvMassLambdaAsCascDghter = xiaod->MassLambda(); lInvMassOmegaMinus = xiaod->MassOmega();}
        if ( lChargeXi > 0 )  { lInvMassXiPlus = xiaod->MassXi(); lInvMassAntiLambdaAsCascDghter  = xiaod->MassAntiLambda(); lInvMassOmegaPlus = xiaod->MassOmega();}

//      cout<<" ID bac  "<<bachTrackXi->GetID()<<" n ID  "<<nTrackXi->GetID()<<" p ID "<<pTrackXi->GetID()<<endl;      
        indexB = bachTrackXi->GetID();
        indexP = pTrackXi->GetID();
        indexN = nTrackXi->GetID();
        lXiMomX = xiaod->MomXiX();
        lXiMomY = xiaod->MomXiY();
        lXiMomZ = xiaod->MomXiZ();
        lV0MomX = xiaod->MomV0X();
        lV0MomY = xiaod->MomV0Y();
        lV0MomZ = xiaod->MomV0Z();

 
        lDcaXiDaughters = xiaod->DcaXiDaughters();

        lDcaV0DaughtersXi = xiaod->DcaV0Daughters();
        lDcaV0ToPrimVertexXi = xiaod->DcaV0ToPrimVertex();

        lDcaBachToPrimVertexXi = xiaod->DcaBachToPrimVertex();
        lPosV0Xi[0] = xiaod->DecayVertexV0X();
        lPosV0Xi[1] = xiaod->DecayVertexV0Y();
        lPosV0Xi[2] = xiaod->DecayVertexV0Z();
        lV0RadiusXi = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0] + lPosV0Xi[1]*lPosV0Xi[1] );

        lV0CosineOfPointingAngle = xiaod->CosPointingAngle( lBestPrimaryVtxPos );

        lV0toXiCosineOfPointingAngle = xiaod->CosPointingAngle( xiaod->GetDecayVertexXi() );

        lXiCosineOfPointingAngle = xiaod->CosPointingAngleXi( lBestPrimaryVtxPos[0], lBestPrimaryVtxPos[1], lBestPrimaryVtxPos[2] );
 
        lDcaPosToPrimVertexXi = xiaod->DcaPosToPrimVertex();
        lDcaNegToPrimVertexXi = xiaod->DcaNegToPrimVertex();
 
      //xiaod->DcaXiToPrimVertex(); // to be checked

        lPosXi[0] = xiaod->DecayVertexXiX();
        lPosXi[1] = xiaod->DecayVertexXiY();
        lPosXi[2] = xiaod->DecayVertexXiZ();

        if (fSecondpart == kXi) lRapXi = xiaod->RapXi();
        else if (fSecondpart == kOmega) lRapXi = xiaod->RapOmega();
      
        // those are for the V0 
        //lEta      = xiaod->Eta();                             
        //lTheta    = xiaod->Theta() *180.0/TMath::Pi();       
        //lPhi      = xiaod->Phi()   *180.0/TMath::Pi();       

        lPosTPCClustersS   = aodpTrackXi->GetTPCnclsS();
        lNegTPCClustersS   = aodnTrackXi->GetTPCnclsS();
        lBachTPCClustersS  = aodbachTrackXi->GetTPCnclsS();

        etaPos = aodpTrackXi->Eta();
        etaNeg = aodnTrackXi->Eta();
        etaBach = aodbachTrackXi->Eta();

      }

      lXiRadius = TMath::Sqrt( lPosXi[0]*lPosXi[0] + lPosXi[1]*lPosXi[1] );

      lV0TransvMom = TMath::Sqrt(TMath::Power(lV0MomX,2)+TMath::Power(lV0MomY,2));
      lV0Mom = TMath::Sqrt(TMath::Power(lV0MomX,2)+TMath::Power(lV0MomY,2)+TMath::Power(lV0MomZ,2));
      lV0RadiusXi = TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0] + lPosV0Xi[1]*lPosV0Xi[1] );

      lXiTransvMom = TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
      lXiMom = TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY + lXiMomZ*lXiMomZ);

      lEtaXi = 0.5*TMath::Log((lXiMom+lXiMomZ)/(lXiMom-lXiMomZ+1.e-13));
      lPhiXi = TMath::Pi()+TMath::ATan2(-lXiMomY,-lXiMomX); // from AliAODRecoDecay.h

      lctau = -1.;
      lctau = TMath::Sqrt(TMath::Power((lPosXi[0]-lBestPrimaryVtxPos[0]),2)+TMath::Power((lPosXi[1]-lBestPrimaryVtxPos[1]),2)+TMath::Power(( lPosXi[2]-lBestPrimaryVtxPos[2]),2));
      if (fSecondpart == kXi) {
        if (lXiMom!=0) lctau = lctau*fPDGXi/lXiMom;
        else lctau = -1.;
      } else if (fSecondpart == kOmega ) {
        if (lXiMom!=0) lctau = lctau*fPDGOmega/lXiMom;
        else lctau = -1.;
      }

      distV0Xi = TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2)+TMath::Power((lPosV0Xi[2]-lPosXi[2]),2));
      lctauV0 = -1.;
      if (lV0Mom!=0) lctauV0 = distV0Xi*fPDGL/lV0Mom;
      ldistTV0Xi = TMath::Sqrt(TMath::Power((lPosV0Xi[0]-lPosXi[0]),2)+TMath::Power((lPosV0Xi[1]-lPosXi[1]),2));

      lPosTPCClusters   = pTrackXi->GetTPCNcls();
      lNegTPCClusters   = nTrackXi->GetTPCNcls();
      lBachTPCClusters  = bachTrackXi->GetTPCNcls();

      if (lPosTPCClusters  < 70) { //AliWarning("Pb / V0 Pos. track has less than minn TPC clusters ... continue!");
        continue; 
      }
      if (lNegTPCClusters  < 70) { //AliWarning("Pb / V0 Neg. track has less than minn TPC clusters ... continue!");
        continue; 
      }
      if (lBachTPCClusters < 70) { //AliWarning("Pb / Bach.   track has less than minn TPC clusters ... continue!");
        continue; 
      }
 
//    if (lXiTransvMom<0.8||lXiTransvMom>8.) continue;
//    if (TMath::Abs(lRapXi)>0.5) continue;

      if (TMath::Abs(etaBach)>0.8) continue;
      if (TMath::Abs(etaPos)>0.8) continue;
      if (TMath::Abs(etaNeg)>0.8) continue;
      if (fkApplyYcutCasc) if (TMath::Abs(lRapXi)>0.5) continue; // basically contained in the previous cut

      sharedFractionTPCclPos = (Float_t) lPosTPCClustersS/lPosTPCClusters;
      sharedFractionTPCclNeg = (Float_t) lNegTPCClustersS/lNegTPCClusters;
      sharedFractionTPCclBach = (Float_t) lBachTPCClustersS/lBachTPCClusters;

      if (sharedFractionTPCclPos>0||sharedFractionTPCclNeg>0||sharedFractionTPCclBach>0) continue;

      // TPC PID: 4-sigma bands on Bethe-Bloch curve
      // Bachelor
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kKaon)) < fnSigmaTPCPIDsecondParticleDau) lIsBachelorKaonForTPC = kTRUE;
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( bachTrackXi,AliPID::kPion)) < fnSigmaTPCPIDsecondParticleDau) lIsBachelorPionForTPC = kTRUE;
      // Negative V0 daughter
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kPion   )) < fnSigmaTPCPIDsecondParticleDau) lIsNegPionForTPC   = kTRUE;
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( nTrackXi,AliPID::kProton )) < fnSigmaTPCPIDsecondParticleDau) lIsNegProtonForTPC = kTRUE;
      // Positive V0 daughter
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kPion   )) < fnSigmaTPCPIDsecondParticleDau) lIsPosPionForTPC   = kTRUE;
      if (TMath::Abs(fPIDResponse->NumberOfSigmasTPC( pTrackXi,AliPID::kProton )) < fnSigmaTPCPIDsecondParticleDau) lIsPosProtonForTPC = kTRUE;


      if (fSecondpart == kXi) {

        if (lChargeXi < 0 && lIsBachelorPionForTPC && lIsNegPionForTPC && lIsPosProtonForTPC) {

          isXi = kTRUE;
          lInvMassLambda = lInvMassLambdaAsCascDghter;

          lInvMassXi = lInvMassXiMinus;

        } else if (lChargeXi > 0 && lIsBachelorPionForTPC && lIsPosPionForTPC && lIsNegProtonForTPC) {

          isaXi = kTRUE;
          lInvMassLambda = lInvMassAntiLambdaAsCascDghter;

          lInvMassXi = lInvMassXiPlus;

        } else  continue;

      } else if (fSecondpart == kOmega) {

        if (lChargeXi < 0 && lIsBachelorKaonForTPC && lIsNegPionForTPC && lIsPosProtonForTPC)  {
          isXi = kTRUE; // this flag is used to set the reconstructed second paticle collection 
          isOmega = kTRUE;
          lInvMassLambda = lInvMassLambdaAsCascDghter;

          lInvMassXi = lInvMassOmegaMinus;

        } else if (lChargeXi > 0 && lIsBachelorKaonForTPC && lIsPosPionForTPC && lIsNegProtonForTPC) {

          isaXi = kTRUE;
          isaOmega = kTRUE;
          lInvMassLambda = lInvMassAntiLambdaAsCascDghter;

          lInvMassXi = lInvMassOmegaPlus;

        } else continue;
      }

      if (lDcaBachToPrimVertexXi<fIPCutBac) continue;  
      // Apply tighter topological cuts
      if (fkTightCutsForCascades) {
        if (lDcaXiDaughters>0.3) continue;
        if (lXiCosineOfPointingAngle<0.9992) continue;
        if (fSecondpart == kXi) {
           if (lXiRadius<1.5 || lXiRadius>100) continue;
           if (lctau>15.) continue;
        } else if (fSecondpart == kOmega) {
           if (lXiRadius<1. || lXiRadius>100) continue; 
           if (lctau>8.) continue;
           if (TMath::Abs(xiaod->MassXi()-fPDGXi)<0.008) continue; // reject omega candiates around the xi mass value
        }
        if (!fkCascadeSideBands) if ((lInvMassLambda<fPDGL-0.005)||(lInvMassLambda>fPDGL+0.005)) continue; // if I want side bands and tight cuts I should not apply the cut        
        if (lV0CosineOfPointingAngle<0.998) continue;//lV0toXiCosineOfPointingAngle;
        if (lV0RadiusXi<3.||lV0RadiusXi>100.) continue;
        if (lDcaV0ToPrimVertexXi<0.1) continue;
        if (lDcaPosToPrimVertexXi<0.1) continue;
        if (lDcaNegToPrimVertexXi<0.1) continue;
      }

      if (isXi) {

        fHistInvMassXiMinus->Fill(lInvMassXiMinus, lXiTransvMom);
        fHistInvMassL->Fill(lInvMassLambdaAsCascDghter, lV0TransvMom);
        fHistpTPCdEdx->Fill(pTrackXi->GetTPCmomentum()*pTrackXi->Charge(), pTrackXi->GetTPCsignal());
        fHistnTPCdEdx->Fill(nTrackXi->GetTPCmomentum()*nTrackXi->Charge(), nTrackXi->GetTPCsignal());
        fHistbTPCdEdx->Fill(bachTrackXi->GetTPCmomentum()*bachTrackXi->Charge(), bachTrackXi->GetTPCsignal());


      } else if (isaXi) {

        fHistInvMassXiPlus->Fill(lInvMassXiPlus, lXiTransvMom);
        fHistInvMassAntiL->Fill(lInvMassAntiLambdaAsCascDghter, lV0TransvMom);
        fHistpTPCdEdx->Fill(pTrackXi->GetTPCmomentum()*pTrackXi->Charge(), pTrackXi->GetTPCsignal());
        fHistnTPCdEdx->Fill(nTrackXi->GetTPCmomentum()*nTrackXi->Charge(), nTrackXi->GetTPCsignal());
        fHistbTPCdEdx->Fill(bachTrackXi->GetTPCmomentum()*bachTrackXi->Charge(), bachTrackXi->GetTPCsignal());

      } 
 

      fHistPosV0TPCClusters->Fill(lPosTPCClusters);  
      fHistNegV0TPCClusters->Fill(lNegTPCClusters);
      fHistBachTPCClusters->Fill(lBachTPCClusters); 

      fHistSharedFrTPCclPos->Fill(sharedFractionTPCclPos);
      fHistSharedFrTPCclNeg->Fill(sharedFractionTPCclNeg);
      fHistSharedFrTPCclBach->Fill(sharedFractionTPCclBach);

      fHistyptXi->Fill(lXiTransvMom,lRapXi); 
      fHistphietaXi->Fill(lPhiXi,lEtaXi);
      fHistptL->Fill(lV0TransvMom);
      fHistptpL->Fill(pTrackXi->Pt());
      fHistptnL->Fill(nTrackXi->Pt());
      fHistptbac->Fill(bachTrackXi->Pt());

    

      // Fill the container for topological cuts optimization
      if (fUseContainer) {
        lContainerCutVars[0] = lDcaXiDaughters;
        lContainerCutVars[1] = lDcaBachToPrimVertexXi;
        lContainerCutVars[2] = lXiCosineOfPointingAngle;
        lContainerCutVars[3] = lXiRadius;
        lContainerCutVars[4] = lInvMassLambda;
        lContainerCutVars[5] = lDcaV0DaughtersXi;
        lContainerCutVars[6] = lV0CosineOfPointingAngle;//lV0toXiCosineOfPointingAngle;
        lContainerCutVars[7] = lV0RadiusXi;
        lContainerCutVars[8] = lDcaV0ToPrimVertexXi;
        lContainerCutVars[9] = lDcaPosToPrimVertexXi;
        lContainerCutVars[10] = lDcaNegToPrimVertexXi;

        lContainerCutVars[13] = lXiTransvMom;
        lContainerCutVars[16] = lBestPrimaryVtxPos[2];
        lContainerCutVars[17] = lcentrality;
        lContainerCutVars[18] = 0.;
        // lContainerCutVars[19] = lBachTPCClusters;
        lContainerCutVars[19] = lctau;
        lContainerCutVars[20] = lctauV0;
        lContainerCutVars[21] = ldistTV0Xi;
        if ( lChargeXi < 0 ) {
          lContainerCutVars[11] = lInvMassXiMinus;
          lContainerCutVars[12] = lInvMassOmegaMinus;//1.63;
          lContainerCutVars[14] = lRapXi;
          lContainerCutVars[15] = lChi2Xi;
          if ( lIsBachelorPionForTPC && lIsPosProtonForTPC && lIsNegPionForTPC ) 
            fCFContCascadeCuts->Fill(lContainerCutVars,0); // for Xi-
      
        } else {
          lContainerCutVars[11] = lInvMassXiPlus;
          lContainerCutVars[12] = lInvMassOmegaPlus;//1.63;
          lContainerCutVars[14] = lRapXi;
          lContainerCutVars[15] = lChi2Xi;
          if ( lIsBachelorPionForTPC && lIsNegProtonForTPC && lIsPosPionForTPC )
            fCFContCascadeCuts->Fill(lContainerCutVars,1); // for Xi+

        }
      }

    
      if (!fkCascadeSideBands) {  
        if (TMath::Abs(lInvMassXi-fPDGsecond) > fMassWindowCascades ) continue; // save only xis in the selected mass range
      } else {
        if ((TMath::Abs(lInvMassXi-fPDGsecond) < fMassWindowCascades)|| (lInvMassXi >= fPDGsecond+4*fMassWindowCascades) || (lInvMassXi <= fPDGsecond-4*fMassWindowCascades) )
          continue;
        if (TMath::Abs(lInvMassLambda-fPDGL)<0.005) continue;  
      }   
      if (lXiTransvMom<fMinPtForCasc||lXiTransvMom>fMaxPtForCasc) continue;

      fEvt->fReconstructedXi[xiCount].xiMomentum[0]  = lXiMomX;
      fEvt->fReconstructedXi[xiCount].xiMomentum[1]  = lXiMomY;
      fEvt->fReconstructedXi[xiCount].xiMomentum[2]  = lXiMomZ;

      for (int xiIndex = 0; xiIndex <3; xiIndex++) {
        fEvt->fReconstructedXi[xiIndex].xiMomentumTruth[xiIndex]  = xiMomentumTruth[xiIndex];
      }
      fEvt->fReconstructedXi[xiCount].xiPt     = lXiTransvMom;
      fEvt->fReconstructedXi[xiCount].xiEta    = lEtaXi;
      fEvt->fReconstructedXi[xiCount].xiPhi    = lPhiXi;
      fEvt->fReconstructedXi[xiCount].xiRap    = lRapXi;
      fEvt->fReconstructedXi[xiCount].xiMass   = lInvMassXi;


      //fEvt->fReconstructedXi[xiCount].xiDCA   = ->DcaV0ToPrimVertex();
      fEvt->fReconstructedXi[xiCount].mcXiOriginType = mcXiOrigin;
      fEvt->fReconstructedXi[xiCount].isXi = isXi; 
      fEvt->fReconstructedXi[xiCount].isaXi = isaXi; 
      fEvt->fReconstructedXi[xiCount].indexB = indexB;
      fEvt->fReconstructedXi[xiCount].indexP = indexP;
      fEvt->fReconstructedXi[xiCount].indexN = indexN;
      fEvt->fReconstructedXi[xiCount].doSkipOver = kFALSE;
      // Store for each track the position at 1.2 radius
      SetSftPosR125( bachTrackXi, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedXi[xiCount].daughterBacShiftedGlobalPosition);
      SetSftPosR125( pTrackXi, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedXi[xiCount].daughterPosShiftedGlobalPosition);
      SetSftPosR125( nTrackXi, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedXi[xiCount].daughterNegShiftedGlobalPosition);

      fEvt->fReconstructedXi[xiCount].daughterBacEtaS = EtaS(fEvt->fReconstructedXi[xiCount].daughterBacShiftedGlobalPosition); 
      fEvt->fReconstructedXi[xiCount].daughterPosEtaS = EtaS(fEvt->fReconstructedXi[xiCount].daughterPosShiftedGlobalPosition);
      fEvt->fReconstructedXi[xiCount].daughterNegEtaS = EtaS(fEvt->fReconstructedXi[xiCount].daughterNegShiftedGlobalPosition);

//    cout<<"Array is xi or axi "<<fEvt->fReconstructedXi[xiCount].isXi<<" "<<fEvt->fReconstructedXi[xiCount].isaXi<<" from flag s xi or axi "<<isXi<<"  "<<isaXi<<endl;
//    cout<<"PDGXi mass "<<fPDGXi<<endl; 
//    cout<<"mass "<<lInvMassXi<<endl;
      fEvt->fReconstructedXi[xiCount].isXiMass = kTRUE;  // selection done earlier 
//    else fEvt->fReconstructedXi[xiCount].isXiMass = kFALSE; // modified to reduce the size of the array

//    cout<<"This Xi in the stack with pT "<<fEvt->fReconstructedXi[xiCount].xiPt<<endl; 
//    cout<<"and the pt of the xi was "<<lXiTransvMom<<endl;
      xiCount++;
      if (fMaxXiMult <= xiCount) {
        cerr<<"Xi counts exceeded "<<fMaxXiMult<<"!"<<endl;
        break;
      }
    } 

    fEvt->fNumberCandidateXi = xiCount;
    fHistXiMultvsCent->Fill(xiCount,lcentrality);

    // Check if Xi share daughters!!
    Int_t nXiWithSharedDaughtersWithXiMass = 0.;
    Int_t nXiWithXiMass = 0.; 
    Int_t daughterIndex = -1;
    Int_t ncheckeddaughters = 0;  
    Int_t checkeddaughters[(fEvt->fNumberCandidateXi)*3];

    for (int i=0; i < fEvt->fNumberCandidateXi*3; i++) checkeddaughters[i] = -1;
    for (int i=0; i < fEvt->fNumberCandidateXi; i++) {
      if (fEvt->fReconstructedXi[i].isXiMass == kTRUE) nXiWithXiMass++;
    }

    for (int i=0; i < fEvt->fNumberCandidateXi; i++) {
      nXiWithSharedDaughtersWithXiMass = 0;
      if (fEvt->fReconstructedXi[i].isXiMass == kTRUE) {
        daughterIndex = fEvt->fReconstructedXi[i].indexB;
        //cout<<"Checking bac dau ! "<<endl;
        nXiWithSharedDaughtersWithXiMass = CheckDaughterTrack ( i, daughterIndex, fEvt, checkeddaughters);
        if (nXiWithXiMass>0.&& nXiWithSharedDaughtersWithXiMass!=-1) fHistFractionOfXiWithSharedDaughters->Fill(nXiWithXiMass, nXiWithSharedDaughtersWithXiMass/(Float_t)  nXiWithXiMass); // check on nXiWithXiMass is redundant we cannot be here if this number is zero
        //if (nXiWithSharedDaughtersWithXiMass>nXiWithXiMass) cout<<"Problem xitot<xishared anal bac "<<nXiWithXiMass<<" "<<nXiWithSharedDaughtersWithXiMass<<endl;
        if (nXiWithSharedDaughtersWithXiMass!=-1) {
          checkeddaughters[ncheckeddaughters] = daughterIndex; 
          ncheckeddaughters++;
        }
      }   
    }

    for (int i=0; i < fEvt->fNumberCandidateXi; i++) {
      nXiWithSharedDaughtersWithXiMass = 0;
      if (fEvt->fReconstructedXi[i].isXiMass == kTRUE) {

        //cout<<"Checking pos dau ! "<<endl;
        daughterIndex = fEvt->fReconstructedXi[i].indexP;
        nXiWithSharedDaughtersWithXiMass = CheckDaughterTrack ( i, daughterIndex, fEvt,checkeddaughters);

        if (nXiWithXiMass>0. && nXiWithSharedDaughtersWithXiMass!=-1) fHistFractionOfXiWithSharedDaughters->Fill(nXiWithXiMass, nXiWithSharedDaughtersWithXiMass/(Float_t)  nXiWithXiMass);
        //if (nXiWithSharedDaughtersWithXiMass>nXiWithXiMass) cout<<"Problem xitot<xishared anal p "<<nXiWithXiMass<<" "<<nXiWithSharedDaughtersWithXiMass<<endl;
        if (nXiWithSharedDaughtersWithXiMass!=-1) {
          checkeddaughters[ncheckeddaughters] = daughterIndex;
          ncheckeddaughters++;
        }
      }
    }
      //cout<<"Checking neg dau ! "<<endl;
    for (int i=0; i < fEvt->fNumberCandidateXi; i++) {
      nXiWithSharedDaughtersWithXiMass = 0;
      if (fEvt->fReconstructedXi[i].isXiMass == kTRUE) { 
        daughterIndex = fEvt->fReconstructedXi[i].indexN;
        nXiWithSharedDaughtersWithXiMass = CheckDaughterTrack ( i, daughterIndex, fEvt, checkeddaughters);
        if (nXiWithXiMass>0. && nXiWithSharedDaughtersWithXiMass!=-1) fHistFractionOfXiWithSharedDaughters->Fill(nXiWithXiMass, nXiWithSharedDaughtersWithXiMass/(Float_t)  nXiWithXiMass);
        //if (nXiWithSharedDaughtersWithXiMass>nXiWithXiMass) cout<<"Problem xitot<xishared anal n "<<nXiWithXiMass<<" "<<nXiWithSharedDaughtersWithXiMass<<endl;
        if (nXiWithSharedDaughtersWithXiMass!=-1) {
          checkeddaughters[ncheckeddaughters] = daughterIndex;
          ncheckeddaughters++;
        }

      }
    }
    // cout<<"Fraction of xis with shared daughters "<< nXiWithSharedDaughtersWithXiMass/ (Float_t) nXiWithXiMass<<endl;
    // cout<<"Fraction of antixis with shared daughters "<< naXiWithSharedDaughtersWithXiMass/(Float_t)  naXiWithXiMass<<endl;
 
    DoPairshCasc(fEvt, lcentrality); 
  } else DoPairshh(fEvt, lcentrality, fieldsign);  

// Post output data
  PostData(1, fOutputContainer);
  PostData(2, fCFContCascadeCuts);

}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskhCascadeFemto::DoPairshCasc (const AliAnalysishCascadeEvent *event, const Float_t lcentrality) {

 
  bool ismassXi;
  bool isXi;
  bool isaXi;
  bool isP;
  bool isaP;
  double dphisbac;
  double detasbac;
  double dphispos;
  double detaspos;
  double dphisneg;
  double detasneg;

//  cout<<"Number of xi candates in this event "<< event->fNumberCandidateXi<<endl;


  int evmultmixed = 0;
  bool multmixedcounted = kFALSE;
  double pairKstar = 0.;

  for (int i=0; i < event->fNumberCandidateXi; i++) {
    for (int j=0; j<event->fNumberCandidateProton; j++) {
      if (TMath::Abs(event->fReconstructedXi[i].indexB) == TMath::Abs(event->fReconstructedProton[j].index) || 
          TMath::Abs(event->fReconstructedXi[i].indexP) == TMath::Abs(event->fReconstructedProton[j].index) ||
          TMath::Abs(event->fReconstructedXi[i].indexN) == TMath::Abs(event->fReconstructedProton[j].index)   ) {

        event->fReconstructedXi[i].doSkipOver = kTRUE;
        event->fReconstructedProton[j].doSkipOver = kTRUE;    

      }      
 
    }
  }


  for (int i=0; i < event->fNumberCandidateXi; i++) { //Start looping over reconstructed Xis in this event 

    ismassXi  = event->fReconstructedXi[i].isXiMass;  
    

    if (!ismassXi) continue;
    if (event->fReconstructedXi[i].doSkipOver) continue;
//    cout<<"Found Xi with correct mass"<<endl;
    isXi = event->fReconstructedXi[i].isXi;
    isaXi = event->fReconstructedXi[i].isaXi;    

    fHistInvMassXiInPairs->Fill(event->fReconstructedXi[i].xiMass,lcentrality);

    for (int eventNumber=0; eventNumber<fnEventsToMix+1; eventNumber++) { // Event buffer loop: eventNumber=0 is the current event, all other eventNumbers are past events 
                                                                          

      // For same event pairs

      if (!multmixedcounted && eventNumber!=0 && ((event+eventNumber)->fNumberCandidateProton)!=0.) evmultmixed++; 
      for (int j=0; j<(event+eventNumber)->fNumberCandidateProton; j++) { // Second particle loop (from past or current event)
        if ((event+eventNumber)->fReconstructedProton[j].doSkipOver) continue; // remove also from mixing, if they were used in Xi we are not sure about the quality  
        isP = (event+eventNumber)->fReconstructedProton[j].isP;    
        isaP = (event+eventNumber)->fReconstructedProton[j].isaP;  

        //Calculate k* for the pair
        pairKstar = CalculateKstar(event->fReconstructedXi[i].xiMomentum, (event+eventNumber)->fReconstructedProton[j].pMomentum, fPDGsecond,fPDGfirst);
          
        // Pair histogramming
        dphisbac = CalculateDphiSatR12m( (event+eventNumber)->fReconstructedProton[j].pShiftedGlobalPosition , event->fReconstructedXi[i].daughterBacShiftedGlobalPosition );
        detasbac = (event+eventNumber)->fReconstructedProton[j].pEtaS - event->fReconstructedXi[i].daughterBacEtaS;

        dphispos = CalculateDphiSatR12m( (event+eventNumber)->fReconstructedProton[j].pShiftedGlobalPosition , event->fReconstructedXi[i].daughterPosShiftedGlobalPosition );
        detaspos = (event+eventNumber)->fReconstructedProton[j].pEtaS - event->fReconstructedXi[i].daughterPosEtaS;

        dphisneg = CalculateDphiSatR12m( (event+eventNumber)->fReconstructedProton[j].pShiftedGlobalPosition , event->fReconstructedXi[i].daughterNegShiftedGlobalPosition );
        detasneg = (event+eventNumber)->fReconstructedProton[j].pEtaS - event->fReconstructedXi[i].daughterNegEtaS;

        // Apply two-track cuts
        if (fkApplyTtc) {
          if (isP&&isXi) {
            if (fFirstpart==kPion && fSecondpart==kXi) {
              if (TMath::Abs(dphisneg)<fDphisMin && TMath::Abs(detasneg)<fDetasMin) continue;
              if (TMath::Abs(dphisbac)<fDphisMin && TMath::Abs(detasbac)<fDetasMin) continue;
            } else if (fFirstpart==kProton && fSecondpart==kXi) {
              if (TMath::Abs(dphispos)<fDphisMin && TMath::Abs(detaspos)<fDetasMin) continue;
            }
          } else if (isaP&&isaXi) {
            if (fFirstpart==kPion && fSecondpart==kXi) {
              if (TMath::Abs(dphispos)<fDphisMin && TMath::Abs(detaspos)<fDetasMin) continue;
              if (TMath::Abs(dphisbac)<fDphisMin && TMath::Abs(detasbac)<fDetasMin) continue;
            } else if (fFirstpart==kProton && fSecondpart==kXi) {
              if (TMath::Abs(dphisneg)<fDphisMin && TMath::Abs(detasneg)<fDetasMin) continue;
            }
          } else if (isaP&&isXi) {
     
            if (fFirstpart==kPion && fSecondpart==kXi) {
              if (TMath::Abs(dphispos)<fDphisMin && TMath::Abs(detaspos)<fDetasMin) continue;
            } else if (fFirstpart==kProton && fSecondpart==kXi) {
              if (TMath::Abs(dphisbac)<fDphisMin && TMath::Abs(detasbac)<fDetasMin) continue;
              if (TMath::Abs(dphisneg)<fDphisMin && TMath::Abs(detasneg)<fDetasMin) continue;
            }

          } else if (isP&&isaXi) {
            if (fFirstpart==kPion && fSecondpart==kXi) {
              if (TMath::Abs(dphisneg)<fDphisMin && TMath::Abs(detasneg)<fDetasMin) continue;
            } else if (fFirstpart==kProton && fSecondpart==kXi) {
              if (TMath::Abs(dphisbac)<fDphisMin && TMath::Abs(detasbac)<fDetasMin) continue;
              if (TMath::Abs(dphispos)<fDphisMin && TMath::Abs(detaspos)<fDetasMin) continue;
            }

          }
           
        }



        if (eventNumber==0) {//Same event pair histogramming
          if (isXi&&isP) {
  //          if (TMath::Abs(dphispos)>0.015&&TMath::Abs(detaspos)>0.1) {
            fHistpXiSignalRealKstar->Fill(pairKstar,lcentrality);//}
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistpXibacDetaSDphiS->Fill(dphisbac,detasbac); 
            if (dphispos!=0.&& detaspos!=0.) fHistpXiposDetaSDphiS->Fill(dphispos,detaspos);
            if (dphisneg!=0.&& detasneg!=0.) fHistpXinegDetaSDphiS->Fill(dphisneg,detasneg); 
//            }
          } else if (isXi&&isaP) {
            fHistapXiSignalRealKstar->Fill(pairKstar,lcentrality);
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistapXibacDetaSDphiS->Fill(dphisbac,detasbac);
            if (dphispos!=0.&& detaspos!=0.) fHistapXiposDetaSDphiS->Fill(dphispos,detaspos);
            if (dphisneg!=0.&& detasneg!=0.) fHistapXinegDetaSDphiS->Fill(dphisneg,detasneg);
//            }
          } else if (isaXi&&isP) {
            fHistpaXiSignalRealKstar->Fill(pairKstar,lcentrality);
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistpaXibacDetaSDphiS->Fill(dphisbac,detasbac); 
            if (dphispos!=0.&& detaspos!=0.) fHistpaXiposDetaSDphiS->Fill(dphispos,detaspos); 
            if (dphisneg!=0.&& detasneg!=0.) fHistpaXinegDetaSDphiS->Fill(dphisneg,detasneg); 
//            }
          } else if (isaXi&&isaP) {
         //   if (TMath::Abs(dphisneg)>0.015&&TMath::Abs(detasneg)>0.1) { 
            fHistapaXiSignalRealKstar->Fill(pairKstar,lcentrality);//}
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistapaXibacDetaSDphiS->Fill(dphisbac,detasbac); 
            if (dphispos!=0.&& detaspos!=0.) fHistapaXiposDetaSDphiS->Fill(dphispos,detaspos); 
            if (dphisneg!=0.&& detasneg!=0.) fHistapaXinegDetaSDphiS->Fill(dphisneg,detasneg);
//            }
          }


        } else {//Mixed-event pair histogramming
/*         if (dphisbac==0.&& detasbac==0.) cout<<"Mixing and Deta nad phi bac are 0!!!!"<<endl;
          if (dphispos==0.&& detaspos==0.) cout<<"Mixing and Deta nad phi pos are 0!!!!"<<endl;
          if (dphisneg==0.&& detasneg==0.) cout<<"Mixing and Deta nad phi neg are 0!!!!"<<endl;
*/ 

          if (isXi&&isP) {
  //          if (TMath::Abs(dphispos)>0.015&&TMath::Abs(detaspos)>0.1) {
            fHistpXiSignalBkgKstar->Fill(pairKstar,lcentrality);//}
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistpXibacDetaSDphiSBkg->Fill(dphisbac,detasbac);
            if (dphispos!=0.&& detaspos!=0.) fHistpXiposDetaSDphiSBkg->Fill(dphispos,detaspos);
            if (dphisneg!=0.&& detasneg!=0.) fHistpXinegDetaSDphiSBkg->Fill(dphisneg,detasneg);
//            }
          } else if (isXi&&isaP) {
            fHistapXiSignalBkgKstar->Fill(pairKstar,lcentrality);
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistapXibacDetaSDphiSBkg->Fill(dphisbac,detasbac);
            if (dphispos!=0.&& detaspos!=0.) fHistapXiposDetaSDphiSBkg->Fill(dphispos,detaspos);
            if (dphisneg!=0.&& detasneg!=0.) fHistapXinegDetaSDphiSBkg->Fill(dphisneg,detasneg);
//            }
          } else if (isaXi&&isP) {
            fHistpaXiSignalBkgKstar->Fill(pairKstar,lcentrality);
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistpaXibacDetaSDphiSBkg->Fill(dphisbac,detasbac);
            if (dphispos!=0.&& detaspos!=0.) fHistpaXiposDetaSDphiSBkg->Fill(dphispos,detaspos);
            if (dphisneg!=0.&& detasneg!=0.) fHistpaXinegDetaSDphiSBkg->Fill(dphisneg,detasneg);
//            }

          } else if (isaXi&&isaP) {
  //          if (TMath::Abs(dphisneg)>0.015&&TMath::Abs(detasneg)>0.1) {
            fHistapaXiSignalBkgKstar->Fill(pairKstar,lcentrality);//}
//            if (pairKstar<0.1) {
            if (dphisbac!=0.&& detasbac!=0.) fHistapaXibacDetaSDphiSBkg->Fill(dphisbac,detasbac);
            if (dphispos!=0.&& detaspos!=0.) fHistapaXiposDetaSDphiSBkg->Fill(dphispos,detaspos);
            if (dphisneg!=0.&& detasneg!=0.) fHistapaXinegDetaSDphiSBkg->Fill(dphisneg,detasneg);
//            }
          }

        }//end mixed event pair histogramming

      }


     }//end event loop
     if (evmultmixed!=0) multmixedcounted = kTRUE; // we count only the first time we go through all events and check if there are protons to mix, then the ncounts will always be the same

   }

   if  (multmixedcounted) fHistMultiplicityOfMixedEvent->Fill(evmultmixed);


}

//--------------------------------------------------------------------------
void AliAnalysisTaskhCascadeFemto::DoPairshh (const AliAnalysishCascadeEvent *event, const Float_t lcentrality, int fieldsign) {

  // FIXME randomizaion to be implemented for same particle pairs
  bool isP1;
  bool isaP1;
  bool isP2;
  bool isaP2;
  double dphis;
  double detas;
  double deta;


//  cout<<"we have for this event a number of xi candates "<< event->fNumberCandidateXi<<endl;


  int evmultmixed = 0;
  bool multmixedcounted = kFALSE;
  double pairKstar = 0.;

  for (int i=0; i<event->fNumberCandidateProton; i++) {

    isP1  = event->fReconstructedProton[i].isP;
    isaP1 = event->fReconstructedProton[i].isaP;

    for (int eventNumber=0; eventNumber<fnEventsToMix+1; eventNumber++) { 
      // For same event pairs
      if (!multmixedcounted && eventNumber!=0 && ((event+eventNumber)->fNumberCandidateProton)!=0.) evmultmixed++; 

      for (int j=0; j<(event+eventNumber)->fNumberCandidateProton; j++) {

        //cout<<" event number "<<eventNumber<<endl;
        if ( (eventNumber == 0) && (j<=i)) continue; 
        isP2 = (event+eventNumber)->fReconstructedProton[j].isP;
        isaP2 = (event+eventNumber)->fReconstructedProton[j].isaP;



        //Calculate k* for the pair
        pairKstar = CalculateKstar(event->fReconstructedProton[i].pMomentum, (event+eventNumber)->fReconstructedProton[j].pMomentum, fPDGsecond,fPDGfirst);

        // Pair histogramming
        dphis = CalculateDphiSatR12m(event->fReconstructedProton[i].pCharge, (event+eventNumber)->fReconstructedProton[j].pCharge, fieldsign,event->fReconstructedProton[i].pPt , (event+eventNumber)->fReconstructedProton[j].pPt, event->fReconstructedProton[i].pPhi, (event+eventNumber)->fReconstructedProton[j].pPhi);

        //dphis = CalculateDphiSatR12m( (event+eventNumber)->fReconstructedProton[j].pShiftedGlobalPosition , event->fReconstructedProton[i].pShiftedGlobalPosition );
        deta = (event+eventNumber)->fReconstructedProton[j].pEta - event->fReconstructedProton[i].pEta;

        // Apply two-track cuts
        if (fkApplyTtc) {
          if ( (isP1&&isP2) || (isaP1&&isaP2)) {
            if (fSecondpart == kProton) {
              if (TMath::Abs(dphis)<fDphisMin && TMath::Abs(deta)<fDetasMin) continue; 
            } else if (fSecondpart == kPion) {
              if (TMath::Sqrt((dphis*dphis)+(deta*deta))<fDphisMin && TMath::Abs(deta)<fDetasMin) continue;

            }
          } 
        }
        if (eventNumber==0) {//Same event pair histogramming
           // simplest pair cut // FIXME not needed maybe
          if (TMath::Abs((event+eventNumber)->fReconstructedProton[j].index) == TMath::Abs(event->fReconstructedProton[i].index)) { //cout<<"In the same event the two particle  have the same index!"<<endl; 

            continue;
          }

        
/*          if (dphisbac==0.&& detasbac==0.) cout<<"Deta nad phi bac are 0!!!! Gl pos pion 0 "<<(event+eventNumber)->fReconstructedProton[j].pShiftedGlobalPosition[0]<<" Gl pos bac 0 "<<event->fReconstructedXi[i].daughterBacShiftedGlobalPosition[0]<<endl;
*/
          if (isP1&&isP2) {
            fHistpXiSignalRealKstar->Fill(pairKstar,lcentrality);
            if (dphis!=0.&& deta!=0.) fHistpXibacDetaSDphiS->Fill(dphis,deta);

          } else if (isaP1&&isaP2) {
            fHistapaXiSignalRealKstar->Fill(pairKstar,lcentrality);
            if (dphis!=0.&& deta!=0.) fHistapaXibacDetaSDphiS->Fill(dphis,deta);
          } else if ((isP1&&isaP2) || (isaP1&&isP2)) { // only because we are combining same particles
            fHistpaXiSignalRealKstar->Fill(pairKstar,lcentrality);
            if (dphis!=0.&& deta!=0.) fHistpaXibacDetaSDphiS->Fill(dphis,deta);

          } 


        } else {//Mixed-event pair histogramming
          if (isP1&&isP2) {
            fHistpXiSignalBkgKstar->Fill(pairKstar,lcentrality);
            if (dphis!=0.&& deta!=0.) fHistpXibacDetaSDphiSBkg->Fill(dphis,deta);

          } else if (isaP1&&isaP2) {
            fHistapaXiSignalBkgKstar->Fill(pairKstar,lcentrality);
            if (dphis!=0.&& deta!=0.) fHistapaXibacDetaSDphiSBkg->Fill(dphis,deta);

          } else if ((isP1&&isaP2) || (isaP1&&isP2)) {
            fHistpaXiSignalBkgKstar->Fill(pairKstar,lcentrality);
            if (dphis!=0.&& deta!=0.) fHistpaXibacDetaSDphiSBkg->Fill(dphis,deta);

          } 

        }

      } // second part

    }//end event loop

    if (evmultmixed!=0) multmixedcounted = kTRUE;
    
  } // first part

  if  (multmixedcounted) fHistMultiplicityOfMixedEvent->Fill(evmultmixed);

  
}

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskhCascadeFemto::ProtonOrigin() {

// to study proton purity: take DCA from global tracks (because for TPC only prim and sec are too similar)
// plot the dcaxy for primary, sec from weak decays and material in pt bins. from these templates one fits 
// data with the ROOT TFractionFitter. From fits you derive the fractions vs pt for each contribution

}

//-----------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskhCascadeFemto::CheckDaughterTrack (Int_t xiIndex, Int_t daughterIndex, const AliAnalysishCascadeEvent *event, Int_t* checkeddaughters ) {


   Int_t nXiWithSharedDaughtersWithXiMass = 0;
//   for (int j=xiIndex+1; j < event->fNumberCandidateXi; j++) {  
   for (int j=0; j < event->fNumberCandidateXi; j++) { 
     if (j!=xiIndex && event->fReconstructedXi[j].isXiMass == kTRUE) {
       // tag the track and check if already checked return.
       for(int icd = 0; icd < (event->fNumberCandidateXi)*3; icd++) { if (daughterIndex == checkeddaughters[icd]) return -1; } // if the index was already checked skip it... 
     
       // compare with all daughters of all the other Xi
       Int_t indexB2 = event->fReconstructedXi[j].indexB;
       Int_t indexP2 = event->fReconstructedXi[j].indexP;
       Int_t indexN2 = event->fReconstructedXi[j].indexN;
       if ( daughterIndex == indexB2) { //cout<<"multiple use of this now is bac ! "<<endl; 
         event->fReconstructedXi[j].doPickOne = kTRUE; 
         event->fReconstructedXi[xiIndex].doPickOne = kTRUE;
       } else if ( daughterIndex == indexP2) { //cout<<"multiple use of this now is pos ! "<<endl; 
         event->fReconstructedXi[j].doPickOne = kTRUE; 
         event->fReconstructedXi[xiIndex].doPickOne = kTRUE;
       } else if ( daughterIndex == indexN2) { //cout<<"multiple use of this now is neg! "<<endl; 
         event->fReconstructedXi[j].doPickOne = kTRUE; 
         event->fReconstructedXi[xiIndex].doPickOne = kTRUE;
       } // cannot be equal to more than one otherwise it means there are two daughters used twice!!
   
     }
   }
   // here one should decide for each daughter checked! choose only one that has the best parameters

//   SelectBestCandidate (event ); // FIXME should be done iteratively
   for (int j=0; j < event->fNumberCandidateXi; j++) {
     if (fEvt->fReconstructedXi[j].isXiMass == kTRUE) {
       if (fEvt->fReconstructedXi[j].doPickOne == kTRUE) {
         nXiWithSharedDaughtersWithXiMass++;
         fEvt->fReconstructedXi[j].doPickOne = kFALSE; // reset
       }
     }
   }
   

   return nXiWithSharedDaughtersWithXiMass;

}
//-----------------------------------------------------------------------------------
void AliAnalysisTaskhCascadeFemto::SelectBestCandidate ( const AliAnalysishCascadeEvent *event) {

 // select for each track the best candidate and flag the worst cascades (doSkipOver)
 // process to be reiterated until we have a set of cascaded not sharing daughters... 
 // to be further implemented and tested

 Float_t deltamassi = 0.;
 Float_t deltamassj = 0.;
 for (int i=0; i < event->fNumberCandidateXi; i++) { 

   if (event->fReconstructedXi[i].doPickOne) {
   
     deltamassi = TMath::Abs(fEvt->fReconstructedXi[i].xiMass - fPDGsecond);
   // Take one Xi and is bac
     for (int j=i+1; j < event->fNumberCandidateXi; j++) {
       if (event->fReconstructedXi[j].doPickOne) {
         deltamassj = TMath::Abs(fEvt->fReconstructedXi[i].xiMass - fPDGXi);
         if (deltamassj > deltamassi) // FIXME add other criteria
         event->fReconstructedXi[j].doSkipOver = kTRUE;
         else  { event->fReconstructedXi[i].doSkipOver = kTRUE; i = j;} // continue the check form from j 
       }// check which is best and set flag to skip
     } 
   }
 }
 return;

}


//-----------------------------------------------------------------------------------------------

double AliAnalysisTaskhCascadeFemto::CalculateKstar(double momentum1[3], double momentum2[3], double mass1, double mass2) { // Jai S

   // Calculate k* for any pair of particles, regardless of whether the

   // particles have the same mass.

   double kstar = 0.;

   double e1 = 0.;

   double e2 = 0.;

   for(int i = 0; i < 3; i++){

     kstar -= pow(momentum1[i]-momentum2[i],2);

     e1 += pow(momentum1[i],2);

     e2 += pow(momentum2[i],2);

   }

   e1 += pow(mass1,2);

   e1 = sqrt(e1);

   e2 += pow(mass2,2);

   e2 = sqrt(e2);

   

   kstar += pow(e1-e2,2);

 

   double totalMomentumSquared = 0;

   for(int i = 0; i < 3; i++){

     totalMomentumSquared -= pow(momentum1[i]+momentum2[i],2);

   }

   totalMomentumSquared += pow(e1+e2,2);

   kstar -= pow((pow(mass1,2)-pow(mass2,2)),2)/totalMomentumSquared;

 

   kstar *= -1.;

   kstar = sqrt(kstar); //At this point, we've actually calculated Qinv

   kstar *= 0.5; // kstar is 0.5*Qinv

   return kstar;
 }
//-----------------------------------------------------------------------------------------------
double AliAnalysisTaskhCascadeFemto::CalculateDphiSatR12m(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2) { // AliFemto framework AliFemtoUser/AliFemtoPairCutRadialDistance.cxx + Dhevan not consistent?

  double rad = 1.2;
  double afsi0b = 0.07510020733*chg1*magSign*rad/ptv1; // 0.075 = 0.3=e in H-L units*0.5=B/2 calculation on notebook
  double afsi1b = 0.07510020733*chg2*magSign*rad/ptv2;

  if (fabs(afsi0b) >=1.) return 9999.; // angle is pi/2 or not defined --> dont cut
  if (fabs(afsi1b) >=1.) return 9999.; // MN modified these two lines returning 9999 and not kTRUE
//  Double_t dps = phi2 - phi1 -TMath::ASin(afsi1b) + TMath::ASin(afsi0b);
//  dps = TVector2::Phi_mpi_pi(dps);
  
  double phi1bis =0.;
  double phi2bis =0.;
  phi1bis = phi1-TMath::ASin(afsi0b); 
  if(phi1bis > 2*PI) phi1bis -= 2*PI;
  if(phi1bis < 0) phi1bis += 2*PI;
  phi2bis = phi2 - TMath::ASin(afsi1b);
  if(phi2bis > 2*PI) phi2bis -= 2*PI;
  if(phi2bis < 0) phi2bis += 2*PI;
  double deltaphi = phi2bis - phi1bis;
  if(deltaphi > PI) deltaphi -= PI;
  if(deltaphi < -PI) deltaphi += PI;

//  cout<<" Dphi "<<dps<<" Dhevan "<<deltaphi<<endl;


  return deltaphi;//dps;

}
//-----------------------------------------------------------------------------------------------
double AliAnalysisTaskhCascadeFemto::CalculateDphiSatR12m(Double_t pos1SftR125[3], Double_t pos2SftR125[3]) { // Hans B

  // Returns delta phi star at R = 1.2 m

  const Float_t distSft = TMath::Sqrt(TMath::Power(pos1SftR125[0] - pos2SftR125[0],2)
                          + TMath::Power(pos1SftR125[1] - pos2SftR125[1],2));
  return 2.0 * TMath::ATan(distSft/2./(125.)); 

}

//-----------------------------------------------------------------------------------------------

void AliAnalysisTaskhCascadeFemto::SetSftPosR125(AliVTrack *track, const Float_t bfield, Double_t priVtx[3], Double_t posSftR125[3] ) {  // Hans B
  // Sets the spatial position of the track at the radius R=1.25m in the shifted coordinate system
  
  
  // Initialize the array to something indicating there was no propagation
  posSftR125[0]=-9999.;
  posSftR125[1]=-9999.;
  posSftR125[2]=-9999.;
   // Make a copy of the track to not change parameters of the track
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);
  
  // The global position of the track
  Double_t xyz[3]={-9999.,-9999.,-9999.};  

  // The radius we want to propagate to, squared, for faster code
  const Float_t rSquared = 125.*125.;


  // Propagation is done in local x of the track
  for (Float_t x = 58.;x<247.;x+=1.){
    // Starts at 83 / Sqrt(2) and goes outwards. 85/Sqrt(2) is the smallest local x
    // for global radius 85 cm. x = 245 is the outer radial limit of the TPC when
    // the track is straight, i.e. has inifinite pt and doesn't get bent. 
    // If the track's momentum is smaller than infinite, it will develop a y-component,
    // which adds to the global radius
    // We don't change the propagation steps to not mess up things!

    // Stop if the propagation was not succesful. This can happen for low pt tracks
    // that don't reach outer radii
    if (!etp.PropagateTo(x,bfield)) { //cout<<"propagation failed!! and etss is "<<EtaS(posSftR125)<<endl; 
      break;
    }
    etp.GetXYZ(xyz); // GetXYZ returns global coordinates

    // Calculate the shifted radius we are at, squared. 
    // Compare squared radii for faster code
    Float_t shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
                                 + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);

    // Roughly reached the radius we want
    if(shiftedRadiusSquared > rSquared){
      
      // Bigger loop has bad precision, we're nearly one centimeter too far, 
      // go back in small steps.
      while (shiftedRadiusSquared>rSquared) {
        // Propagate a mm inwards
        x-=.1;
        if (!etp.PropagateTo(x,bfield)){
          // Propagation failed but we're already with a
          // cm precision at R=1.25m so we only break the 
          // inner loop
          //cout<<"propagation failed!! and etss is "<<EtaS(posSftR125)<<endl; 

          break;
        }
        // Get the global position
        etp.GetXYZ(xyz);
        // Calculate shifted radius, squared
        shiftedRadiusSquared = (xyz[0]-priVtx[0])*(xyz[0]-priVtx[0])
                           + (xyz[1]-priVtx[1])*(xyz[1]-priVtx[1]);
      }
      // We reached R=1.25m with a precission of a cm to a mm,
      // set the spatial position
      posSftR125[0]=xyz[0]-priVtx[0];
      posSftR125[1]=xyz[1]-priVtx[1];
      posSftR125[2]=xyz[2]-priVtx[2];
      //cout<<" Pos 125 cm in function end "<<posSftR125[0]<<" "<<posSftR125[1]<<" "<<posSftR125[2]<<endl;
      // Done
      return;
    } // End of if roughly reached radius
  } // End of coarse propagation loop
}

//----------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskhCascadeFemto::ThetaS( Double_t posSftR125[3] ) const { // Hans B
  // Returns the longitudinal angle of the particles propagated
  // position at R=1.25m. See
  // https://edms.cern.ch/file/406391/2/ALICE-INT-2003-038.pdf
  // for the ALICE coordinate system. Theta is zero at positive z,
  // pi/2 at z = 0 aka the xy plane and pi at negative z 

  // R^    ^  
  //  |   /
  //  |?'/
  //  | / ?
  //  |/____>z
  // 
  // Let's compute ?' and ? = pi/2 - ?'
  // where ?' can even be and should 
  // sometimes be negative
  // tan(?') = z/R
  // ?' = arctan(z/R)
  // ? = pi/2 - ?'
  //   = pi/2 - arctan(z/R)
  // Note that in the doc above theta
  // is calculated as arccos(z/sqrt(x^2+y^2+z^2))

  // Array of positions is 85,105,125,..cm,
  // we take the z position at R=1.25m
  // return TMath::Pi()/2. - TMath::ATan(fXshifted[2][2]/125.);
  return TMath::Pi()/2. - TMath::ATan(posSftR125[2]/125.); // ok here R is really there --> transverse plane 
}
//_______________________________________________________________
Double_t AliAnalysisTaskhCascadeFemto::EtaS( Double_t posSftR125[3] ) const {  // Hans B
  // Returns the corresponding eta of a pri. part. 
  // with this particles pos at R=1.25m

  // http://en.wikipedia.org/wiki/Pseudorapidity
  // ? = -ln[ tan(?/2)]
  // printf("z: %+04.0f, thetaS %+03.2f etaS %+1.2f\n"
  // ,fXshifted[2][2],ThetaS(),-TMath::Log( TMath::Tan(ThetaS()/2.) ));
  return -TMath::Log( TMath::Tan(ThetaS(posSftR125 )/2.) );
}



//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskhCascadeFemto::Terminate(const Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  if (!GetOutputData(0)) return;
}

