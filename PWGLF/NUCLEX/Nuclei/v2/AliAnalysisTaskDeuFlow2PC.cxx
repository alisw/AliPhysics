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
//
//////////////////////////////////////////////////////////////////////////
/// This analysis task is taken and adapted from  
/// Ramona.Lea@cern.ch
/// m.nicassio@gsi.de maria.nicassio@cern.ch 
/// First version of the committed code 
/// March 2018
///
/// 
/// Parts of this code are taken and adapted from the code below 
///
/// PWGCF/FEMTOSCOPY/AliFemto and AliFemtoUser (AliFemto framework)
/// PWGCF/FEMTOSCOPY/V0LamAnalysis (Jai)
/// PWGCF/FEMTOSCOPY/Chaoticity (Dhevan)
/// PWGCF/FEMTOSCOPY/K0S/PLamAnalysis (Hans) 
/// PWGCF/FEMTOSCOPY/hCascade (Mariella)
///
/// Possibility to run on
/// - first particle: pi,k, p, d
/// - second particle: //
///   to be added other particles
/// - AODs
/// - ESDs (not mantained)
///--> Purity particle --> to be implemented
/// To be done
/// - read MC truth  
///    --> Momentum resolution correction histos
///  Some more improvements will come
/////////////////////////////////////////////////////////////////////////

class AliESDVertex;

class AliAODVertex;

#include "TChain.h"

#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TDatabasePDG.h"
#include "TRandom.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliInputEventHandler.h"

#include "AliCentrality.h"
#include "AliPIDResponse.h"

#include "AliMultSelection.h"

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

#include "AliEventCuts.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"

#include "AliVHeader.h"
#include "AliVVertex.h"
#include "AliAnalysisTaskDeuFlow2PC.h"
#include "AliAnalysishDEventCollection.h"

#define PI 3.1415927
 
using namespace std;


ClassImp(AliAnalysisTaskDeuFlow2PC)

AliAnalysisTaskDeuFlow2PC::AliAnalysisTaskDeuFlow2PC():AliAnalysisTaskSE(),
  fEventCuts(0),
  fESDevent(NULL),
  fAODevent(NULL),
  fAnalysisType("AOD"),
  fCollidingSystem("pp"),
  fYear(2015),
  fHMtrigger(kFALSE),
  fFilterBit(128),
  fFirstpart(kKaon),
  fSecondpart(kProton),

  fnSigmaTPCPIDfirstParticle(3.),
  fnSigmaTPCTOFPIDfirstParticle(3.),
  fnSigmaTPCPIDsecondParticle(3.),
  fnSigmaTPCTOFPIDsecondParticle(3.),

  fReadMCTruth(kFALSE),
  fUseContainer(kFALSE),
  fUseStandardCuts(0),
  fkApplyTtc(0),

  fDphisMin(0),
  fDetasMin(0),
						       
  fMomemtumLimitForTOFPIDfirst(0.4),
  fMomemtumLimitForTOFPIDsecond(0.8),

  fkApplyRatioCrRnFindCut(kFALSE),
  fkCutOnTPCIP(kTRUE),  

  fIPCutxyPrim(0.1),
  fIPCutzPrim(0.15),
  fIPCutxySec(0.1),
  fIPCutzSec(0.15),

  fMinPtForPrim(0.7),
  fMaxPtForPrim(4.),
  fMinPtForSec(0.8),
  fMaxPtForSec(10.),
  
  fRadius(1.2),
  
  fkDoSphericity(kFALSE),
  fkPropagateGlobal(kFALSE),
  
  fESDtrackCuts(0),
  fPIDResponse(0),

  fCentrLowLim(0),
  fCentrUpLim(0),

  farrGT(0),
  fTrackBufferSize(20200), // was 18k
  
  fEventColl(0x0),
  fEventCollwSp(0x0),
  fEvt(0x0),
  
  fMaxFirstMult(3000), // 1000 for protons 
  fMaxSecondMult(20),  // was 100
  
  fPDGp(0.938272046),
  fPDGk(0.493677),
  fPDGd(1.875612762),
  fPDGpi(0.13957),

  fPDGfirst(0.),
  fPDGsecond(0.),
  
  fPDGCodefirst(0.),
  fPDGCodesecond(0.),

  fzVertexBins(10),
  fnCentBins(20),
  fnSphericityBins(3),
  
  fnEventsToMix(7),
  
  fHistEventMultiplicity(0), 
  hmult(0),  
  fHistCentrality(0),                       
  fHistVertexDistribution(0),                         
  fHistSphericity(0),               
  fHistMultiplicityOfMixedEvent(0),     
  
  fHistTriggptvsCentrality(0),
  
  fHistTPCdEdx(0),
  fHistFirstTPCdEdx(0),
  fHistSecondTPCdEdx(0),
  
  fHistnTPCCrossedRFirst(0),                
  fHistRationTPCCrossedRnFindFirst(0),      
  fHistSharedFrTPCclFirst(0),               

  fHistnTPCCrossedRSecond(0),               
  fHistRationTPCCrossedRnFindSecond(0),     
  fHistSharedFrTPCclSecond(0),              

  fHistyptFirst(0),                         
  fHistyptSecond(0),                        
 
  fHistphietaFirst(0),                      
  fHistphietaSecond(0),                     

  fHistIPtoPVxyzTPCFirst(0),                
  fHistIPtoPVxyzGlobalFirst(0),             
  
  fHistIPtoPVxyzTPCSecond(0),               
  fHistIPtoPVxyzGlobalSecond(0),            
  
  fHistFirstTOFmisvspt(0),                  
  fHistFirstTOFmisvsp(0),                   
  fHistFirstTOFnsigmavspt(0),               
  fHistFirstTOFnsigmavsp(0),                
  fHistFirstTOFsignalvsp(0),                
  fHistFirstTOFsignalvspt(0),               
  
  fHistFirstTOFTPCsignalvspt(0),            
  fHistFirstMultvsCent(0),               
  
  fHistSecondTOFmisvspt(0),                 
  fHistSecondTOFmisvsp(0),                  
  fHistSecondTOFnsigmavspt(0),              
  fHistSecondTOFnsigmavsp(0),               
  fHistSecondTOFsignalvsp(0),               
  fHistSecondTOFsignalvspt(0),              
  
  fHistSecondTOFTPCsignalvspt(0),           
  fHistSecondMultvsCent(0),                

  fHistFirstTPCdEdxAfter(0),
  fHistSecondTPCdEdxAfter(0),     
  fHistFirstTOFTPCsignalvsptAfter(0), 
  fHistSecondTOFTPCsignalvsptAfter(0), 
 
  fHistFirstMassTOFvsPt3sTPC(0), 
  fHistSecondMassTOFvsPt3sTPC(0),  

  fHistFirstMassTOFvsPt3sTPC3sTOF(0), 
  fHistSecondMassTOFvsPt3sTPC3sTOF(0),  

  fHistSparseSignal(0),  
  fHistSparseBkg(0),      
  
  tSignP1(0),
  tSignP2(0),
  tCentrality(0),
  tDCAxyP1(0),
  tDCAzP1(0), 
  tDCAxyP2(0), 
  tDCAzP2(0), 
  tKtpair(0),
  tkStar(0),
  tptP1(0),
  tptP2(0),
  tDEta(0),
  tDPhiStar(0),
  tDPhi(0),
  tMassPair(0),
  tSphericity(0),
  tMassS(0),
  tDTheta(0),
  tMCtruepair(0),
  tMCSameMother(0),
  tMCMotherP1(0),
  tMCMotherP2(0),
  tMCptcTypeP1(0),
  tMCptcTypeP2(0),
  tMCSameGM(0),
  tMotherPDG(0),
  tpdgcodeP1(0), 
  tpdgcodeP2(0),
  tKstarGen(0),

  fHistTrackBufferOverflow(0),             
  fOutputContainer(NULL)

{
  fPDGp  = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  fPDGpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  fPDGk  = TDatabasePDG::Instance()->GetParticle(321)->Mass();

  fPDGd  = 1.875612762;
  

  // default Constructor
  // Define input and output slots here
}

//-------------------------------------------------------------------------------

AliAnalysisTaskDeuFlow2PC::AliAnalysisTaskDeuFlow2PC(const char *name):
  AliAnalysisTaskSE(name),
  fEventCuts(0),
  fESDevent(NULL),
  fAODevent(NULL),
  fAnalysisType("AOD"),
  fCollidingSystem("pp"),
  fYear(2015),
  fHMtrigger(kFALSE),
  fFilterBit(128),
  fFirstpart(kKaon),
  fSecondpart(kProton),

  fnSigmaTPCPIDfirstParticle(3.),
  fnSigmaTPCTOFPIDfirstParticle(3.),
  fnSigmaTPCPIDsecondParticle(3.),
  fnSigmaTPCTOFPIDsecondParticle(3.),

  fReadMCTruth(kFALSE),
  fUseContainer(kFALSE),
  fUseStandardCuts(0),
  fkApplyTtc(0),

  fDphisMin(0),
  fDetasMin(0),
						       
  fMomemtumLimitForTOFPIDfirst(0.4),
  fMomemtumLimitForTOFPIDsecond(0.8),

  fkApplyRatioCrRnFindCut(kFALSE),
  fkCutOnTPCIP(kTRUE),  

  fIPCutxyPrim(0.1),
  fIPCutzPrim(0.15),
  fIPCutxySec(0.1),
  fIPCutzSec(0.15),

  fMinPtForPrim(0.7),
  fMaxPtForPrim(4.),
  fMinPtForSec(0.8),
  fMaxPtForSec(10.),
  
  fRadius(1.2),
  
  fkDoSphericity(kFALSE),
  fkPropagateGlobal(kFALSE),
  
  fESDtrackCuts(0),
  fPIDResponse(0),

  fCentrLowLim(0),
  fCentrUpLim(0),

  farrGT(0),
  fTrackBufferSize(20200), // was 18k
  
  fEventColl(0x0),
  fEventCollwSp(0x0),
  fEvt(0x0),
  
  fMaxFirstMult(3000), // 1000 for protons 
  fMaxSecondMult(20),  // was 100
  
  fPDGp(0.938272046),
  fPDGk(0.493677),
  fPDGd(1.875612762),
  fPDGpi(0.13957),

  fPDGfirst(0.),
  fPDGsecond(0.),
  
  fPDGCodefirst(0.),
  fPDGCodesecond(0.),

  fzVertexBins(10),
  fnCentBins(20),
  fnSphericityBins(3),
  
  fnEventsToMix(7),
  
  fHistEventMultiplicity(0), 
  hmult(0),  
  fHistCentrality(0),                       
  fHistVertexDistribution(0),                         
  fHistSphericity(0),               
  fHistMultiplicityOfMixedEvent(0),     
  
  fHistTriggptvsCentrality(0),
  
  fHistTPCdEdx(0),
  fHistFirstTPCdEdx(0),
  fHistSecondTPCdEdx(0),
  
  fHistnTPCCrossedRFirst(0),                
  fHistRationTPCCrossedRnFindFirst(0),      
  fHistSharedFrTPCclFirst(0),               

  fHistnTPCCrossedRSecond(0),               
  fHistRationTPCCrossedRnFindSecond(0),     
  fHistSharedFrTPCclSecond(0),              

  fHistyptFirst(0),                         
  fHistyptSecond(0),                        
 
  fHistphietaFirst(0),                      
  fHistphietaSecond(0),                     

  fHistIPtoPVxyzTPCFirst(0),                
  fHistIPtoPVxyzGlobalFirst(0),             
  
  fHistIPtoPVxyzTPCSecond(0),               
  fHistIPtoPVxyzGlobalSecond(0),            
  
  fHistFirstTOFmisvspt(0),                  
  fHistFirstTOFmisvsp(0),                   
  fHistFirstTOFnsigmavspt(0),               
  fHistFirstTOFnsigmavsp(0),                
  fHistFirstTOFsignalvsp(0),                
  fHistFirstTOFsignalvspt(0),               
  
  fHistFirstTOFTPCsignalvspt(0),            
  fHistFirstMultvsCent(0),               
  
  fHistSecondTOFmisvspt(0),                 
  fHistSecondTOFmisvsp(0),                  
  fHistSecondTOFnsigmavspt(0),              
  fHistSecondTOFnsigmavsp(0),               
  fHistSecondTOFsignalvsp(0),               
  fHistSecondTOFsignalvspt(0),              
  
  fHistSecondTOFTPCsignalvspt(0),           
  fHistSecondMultvsCent(0),                

  fHistFirstTPCdEdxAfter(0),
  fHistSecondTPCdEdxAfter(0),     
  fHistFirstTOFTPCsignalvsptAfter(0), 
  fHistSecondTOFTPCsignalvsptAfter(0), 
 
  fHistFirstMassTOFvsPt3sTPC(0), 
  fHistSecondMassTOFvsPt3sTPC(0),  

  fHistFirstMassTOFvsPt3sTPC3sTOF(0), 
  fHistSecondMassTOFvsPt3sTPC3sTOF(0),  

  fHistSparseSignal(0),  
  fHistSparseBkg(0),      
  
  tSignP1(0),
  tSignP2(0),
  tCentrality(0),
  tDCAxyP1(0),
  tDCAzP1(0), 
  tDCAxyP2(0), 
  tDCAzP2(0), 
  tKtpair(0),
  tkStar(0),
  tptP1(0),
  tptP2(0),
  tDEta(0),
  tDPhiStar(0),
  tDPhi(0),
  tMassPair(0),
  tSphericity(0),
  tMassS(0),
  tDTheta(0),
  tMCtruepair(0),
  tMCSameMother(0),
  tMCMotherP1(0),
  tMCMotherP2(0),
  tMCptcTypeP1(0),
  tMCptcTypeP2(0),
  tMCSameGM(0),
  tMotherPDG(0),
  tpdgcodeP1(0), 
  tpdgcodeP2(0),
  tKstarGen(0),

  fHistTrackBufferOverflow(0),             
  fOutputContainer(NULL)

  
{
  fPDGp  = TDatabasePDG::Instance()->GetParticle(2212)->Mass(); 
  fPDGpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  fPDGk  = TDatabasePDG::Instance()->GetParticle(321)->Mass();

  fPDGd  = 1.875612762;

  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

//--------------------------------------------------------------------------------

AliAnalysisTaskDeuFlow2PC::~AliAnalysisTaskDeuFlow2PC() {
  //
  // Destructor
  //
  if (fOutputContainer){
    delete fOutputContainer;     
    fOutputContainer = 0x0;    
  }
  if (fESDtrackCuts){
    delete fESDtrackCuts;
    fESDtrackCuts = 0x0;
  }
  if(fHistSparseSignal){
    delete fHistSparseSignal;
    fHistSparseSignal = 0x0;
  }
  if(fHistSparseBkg) {
    delete fHistSparseBkg;
    fHistSparseBkg = 0x0;
  }
  if (farrGT)
    delete[] farrGT;
  farrGT=0;
  
  if(fkDoSphericity == kFALSE){
    for (unsigned short i=0; i<fzVertexBins; i++) {
      for (unsigned short j=0; j<fnCentBins; j++) {
	delete fEventColl[i][j];
      }
      delete[] fEventColl[i];
    }
  }

  else{
    for (unsigned short i=0; i<fzVertexBins; i++) {
      for (unsigned short j=0; j<fnCentBins; j++) {
	for (unsigned short z=0; j<fnSphericityBins; z++) {
	  delete fEventCollwSp[i][j][z];
	}
	delete fEventCollwSp[i][j];
      }
      delete[] fEventCollwSp[i];
    }
  }


  delete[] fEventColl;
  delete[] fEventCollwSp;
 
}

//--------------------------------------------------------------------------------

void AliAnalysisTaskDeuFlow2PC::UserCreateOutputObjects() {

  OpenFile(1);
  // // PID object
  // AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  // AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  // fPIDResponse = inputHandler->GetPIDResponse();

  // Event mixing part
  // Standard 2 variables
  if(fkDoSphericity == kFALSE){
    fEventColl = new AliAnalysishDEventCollection **[fzVertexBins]; 
  
    for (unsigned short i=0; i<fzVertexBins; i++) {
      fEventColl[i] = new AliAnalysishDEventCollection *[fnCentBins];
      for (unsigned short j=0; j<fnCentBins; j++) {
	fEventColl[i][j] = new AliAnalysishDEventCollection(fnEventsToMix+1, fMaxFirstMult, fMaxSecondMult);
      }
    }
  }

  else{
    fEventCollwSp = new AliAnalysishDEventCollection ***[fzVertexBins]; 
  
    for (unsigned short i=0; i<fzVertexBins; i++) {
      fEventCollwSp[i] = new AliAnalysishDEventCollection **[fnCentBins];
      for (unsigned short j=0; j<fnCentBins; j++) {
	fEventCollwSp[i][j] = new AliAnalysishDEventCollection*[fnCentBins];
	for(UChar_t k=0;k<fnSphericityBins;k++){
	  fEventCollwSp[i][j][k] = new AliAnalysishDEventCollection(fnEventsToMix+1, fMaxFirstMult, fMaxSecondMult);
	}
      }
    }
  }
  
 

  // Set PDG for first and second particle
  if (fFirstpart == kKaon) fPDGfirst = fPDGk; 
  else if (fFirstpart == kProton) fPDGfirst = fPDGp;
  else if (fFirstpart == kPion) fPDGfirst = fPDGpi;
  else if (fFirstpart == kDeuteron) fPDGfirst = fPDGd;
  else { cout<<" First particle not known or not possible to deal with in this task for the moment! "<<endl; return;}
  
  if (fSecondpart == kPion) fPDGsecond = fPDGpi;
  else if (fSecondpart == kProton) fPDGsecond = fPDGp;
  else if (fSecondpart == kKaon) fPDGsecond = fPDGk;
  else if (fSecondpart == kDeuteron) fPDGsecond = fPDGd;
  else if (fSecondpart == kAny) fPDGsecond = 0;
  else { cout<<" Second particle not known or not possible to deal with in this task for the moment! "<<endl; return;}

  //SetPDGCode for first and second particle
  // fPDGp  = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  // fPDGpi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
  // fPDGk  = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  

  if      (fFirstpart == kPion)     fPDGCodefirst = 211; 
  else if (fFirstpart == kKaon)     fPDGCodefirst = 321;
  else if (fFirstpart == kProton)   fPDGCodefirst = 2212;
  else if (fFirstpart == kDeuteron) fPDGCodefirst = 1000010020;
  else { cout<<" First particle not known or not possible to deal with in this task for the moment! "<<endl; return;}
  
  if      (fSecondpart == kPion)    fPDGCodesecond = 211;       
  else if (fSecondpart == kKaon)    fPDGCodesecond = 321;       
  else if (fSecondpart == kProton)  fPDGCodesecond = 2212;      
  else if (fSecondpart == kDeuteron)fPDGCodesecond = 1000010020;
  else if (fSecondpart == kAny)     fPDGCodesecond = 0;
  else { cout<<" Second particle not known or not possible to deal with in this task for the moment! "<<endl; return;}
  
  
  // Store pointer to global tracks
  farrGT = new Int_t[fTrackBufferSize];
  
  //Define and fill the OutputContainer
  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);
 
  
  // Create histograms 
  
  fHistEventMultiplicity   = new TH1F( "fHistEventMultiplicity" , "Nb of Events" , 13 , 0.5,13.5);
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(1,"All Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(2,"Events w/PV and PID response");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(3,"Events w/|Vz|<10cm");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(4,"Centrality acc");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(5,"w/o PileUp");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(6,"Any Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(7,"Central Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(8,"Semi-Central Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(9,"MB Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(10,"kInt7 Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(11,"HM Events");
  fHistEventMultiplicity->GetXaxis()->SetBinLabel(12,"Is Selected Events");
  fOutputContainer->Add(fHistEventMultiplicity);

  hmult = new TH1I("hmult","Multiplicity distribution (after cuts on event)",30000,-0.5,2999.5);
  hmult->GetXaxis()->SetTitle("Number of tracklets");
  fOutputContainer->Add(hmult);

  fHistCentrality = new TH1F("fHistCentrality", "Number of events", 10001, -0.5, 100.5);
  fHistCentrality ->GetXaxis()->SetTitle("Centrality");
  fHistCentrality ->GetYaxis()->SetTitle("Entries");
  fOutputContainer->Add(fHistCentrality);

  fHistVertexDistribution = new TH1F("fHistVertexDistribution", "Primary vertex distribution", 40, -20., 20.);
  fHistVertexDistribution ->GetXaxis()->SetTitle("z_{v} (cm)");
  fHistVertexDistribution ->GetYaxis()->SetTitle("Entries");
  fOutputContainer->Add(fHistVertexDistribution); 

  fHistSphericity = new TH1F("fHistSphericity", "Sphericity Distribution", 40, 0., 1.);
  fHistSphericity ->GetXaxis()->SetTitle("Sphericity");
  fHistSphericity ->GetYaxis()->SetTitle("Entries");
  fOutputContainer->Add(fHistSphericity); 
   
  fHistMultiplicityOfMixedEvent = new TH1F("fHistMultiplicityOfMixedEvent","fHistMultiplicityOfMixedEvent",100, 0.,100.);
  fOutputContainer->Add(fHistMultiplicityOfMixedEvent);
     
  fHistTriggptvsCentrality = new TH2F("fHistTriggptvsCentrality","fHistTriggptvsCentrality",600,-30,30,100,0,100);
  fOutputContainer->Add(fHistTriggptvsCentrality);

  fHistTPCdEdx = new TH2F("fHistTPCdEdx", "fHistTPCdEdx", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistTPCdEdx);

  fHistFirstTPCdEdx = new TH2F("fHistFirstTPCdEdx", "fHistFirstTPCdEdx", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistFirstTPCdEdx);

  fHistSecondTPCdEdx = new TH2F("fHistSecondTPCdEdx", "fHistSecondTPCdEdx", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistSecondTPCdEdx);

  fHistnTPCCrossedRFirst = new TH1F("fHistnTPCCrossedRFirst","fHistnTPCCrossedRFirst",200, 0.,200.);
  fOutputContainer->Add(fHistnTPCCrossedRFirst);
 
  fHistRationTPCCrossedRnFindFirst = new TH1F("fHistRationTPCCrossedRnFindFirst","fHistRationTPCCrossedRnFindFirst",2000, 0.,170.);
  fOutputContainer->Add(fHistRationTPCCrossedRnFindFirst);
  
  fHistSharedFrTPCclFirst = new TH1F("fHistSharedFrTPCclFirst","fHistSharedFrTPCclFirst",400, 0.,2.);
  fOutputContainer->Add(fHistSharedFrTPCclFirst);

  fHistnTPCCrossedRSecond = new TH1F("fHistnTPCCrossedRSecond","fHistnTPCCrossedRSecond",200, 0.,200.);
  fOutputContainer->Add(fHistnTPCCrossedRSecond);
  
  fHistRationTPCCrossedRnFindSecond = new TH1F("fHistRationTPCCrossedRnFindSecond","fHistRationTPCCrossedRnFindSecond",2000, 0.,170.);
  fOutputContainer->Add(fHistRationTPCCrossedRnFindSecond);
  
  fHistSharedFrTPCclSecond = new TH1F("fHistSharedFrTPCclSecond","fHistSharedFrTPCclSecond",400, 0.,2.);
  fOutputContainer->Add(fHistSharedFrTPCclSecond);

  fHistyptFirst = new TH3F("fHistyptFirst","fHistyptFirst",100,0.,10.,40,-2.,2.,10,0.,100.);
  fOutputContainer->Add(fHistyptFirst);
  
  fHistyptSecond = new TH3F("fHistyptSecond","fHistyptSecond",100,0.,10.,40,-2.,2.,10,0.,100.);
  fOutputContainer->Add(fHistyptSecond);
  
  fHistphietaFirst = new TH2F("fHistphietaFirst","fHistphietaFirst",400,0.,2*TMath::Pi(),100, -2, 2.);
  fOutputContainer->Add(fHistphietaFirst);

  fHistphietaSecond = new TH2F("fHistphietaSecond","fHistphietaSecond",400,0.,2*TMath::Pi(),100, -2, 2.);
  fOutputContainer->Add(fHistphietaSecond);

  fHistIPtoPVxyzTPCFirst = new TH2F("fHistIPtoPVxyzTPCFirst","fHistIPtoPVxyzTPCFirst",200, -5.,5., 200, -5., 5.);
  fOutputContainer->Add(fHistIPtoPVxyzTPCFirst);

  fHistIPtoPVxyzGlobalFirst = new TH2F("fHistIPtoPVxyzGlobalFirst","fHistIPtoPVxyzGlobalFirst",200, -5.,5., 200, -5., 5.);
  fOutputContainer->Add(fHistIPtoPVxyzGlobalFirst);

  fHistIPtoPVxyzTPCSecond = new TH2F("fHistIPtoPVxyzTPCSecond","fHistIPtoPVxyzTPCSecond",200, -5.,5., 200, -5., 5.);
  fOutputContainer->Add(fHistIPtoPVxyzTPCSecond);

  fHistIPtoPVxyzGlobalSecond = new TH2F("fHistIPtoPVxyzGlobalSecond","fHistIPtoPVxyzGlobalSecond",200, -5.,5., 200, -5., 5.);
  fOutputContainer->Add(fHistIPtoPVxyzGlobalSecond);

  fHistFirstTOFmisvspt = new TH2F("fHistFirstTOFmisvspt", "fHistFirstTOFmisvspt", 200, 0., 10., 101, 0.0, 1.01);
  fOutputContainer->Add(fHistFirstTOFmisvspt);
  
  fHistFirstTOFmisvsp = new TH2F("fHistFirstTOFmisvsp", "fHistFirstTOFmisvsp", 200, 0., 10., 101, 0.0, 1.01);
  fOutputContainer->Add(fHistFirstTOFmisvsp);
  
  fHistFirstTOFnsigmavspt = new TH2F("fHistFirstTOFnsigmavspt", "fHistFirstTOFnsigmavspt", 200, 0., 10., 100, -5. , 5);
  fOutputContainer->Add(fHistFirstTOFnsigmavspt);
  
  fHistFirstTOFnsigmavsp = new TH2F("fHistFirstTOFnsigmavsp", "fHistFirstTOFnsigmavsp", 200, 0., 10., 100, -5., 5.);
  fOutputContainer->Add(fHistFirstTOFnsigmavsp);
  
  fHistFirstTOFsignalvsp = new TH2F("fHistFirstTOFsignalvsp", "fHistFirstTOFsignalvsp", 200, 0., 10, 100,-3000,3000);//1000 , 0., 2.);
  fHistFirstTOFsignalvsp->GetYaxis()->SetTitle("t_{meas}-t_{0}-t_{exoected} (ps)");
  fHistFirstTOFsignalvsp->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fOutputContainer->Add(fHistFirstTOFsignalvsp);
  
  fHistFirstTOFsignalvspt = new TH2F("fHistFirstTOFsignalvspt", "fHistFirstTOFsignalvspt", 200, 0., 10., 100,-3000,3000);//1000, 0.0, 2.);
  fHistFirstTOFsignalvspt->GetYaxis()->SetTitle("t_{meas}-t_{0}-t_{expected} (ps)");
  fHistFirstTOFsignalvspt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})"); 
  fOutputContainer->Add(fHistFirstTOFsignalvspt);

  fHistFirstTOFTPCsignalvspt = new TH2F("fHistFirstTOFTPCsignalvspt","fHistFirstTOFTPCsignalvspt",200, 0., 10., 100, 0., 20);
  fOutputContainer->Add(fHistFirstTOFTPCsignalvspt);

  fHistFirstMultvsCent = new TH2F("fHistFirstMultvsCent","fHistFirstMultvsCent",400,0.,2000.,10,0.,100.);
  fOutputContainer->Add(fHistFirstMultvsCent);
 
  fHistSecondTOFmisvspt = new TH2F("fHistSecondTOFmisvspt", "fHistSecondTOFmisvspt", 200, 0., 10., 101, 0.0, 1.01);
  fOutputContainer->Add(fHistSecondTOFmisvspt);
  
  fHistSecondTOFmisvsp = new TH2F("fHistSecondTOFmisvsp", "fHistSecondTOFmisvsp", 200, 0., 10., 101, 0.0, 1.01);
  fOutputContainer->Add(fHistSecondTOFmisvsp);
  
  fHistSecondTOFnsigmavspt = new TH2F("fHistSecondTOFnsigmavspt", "fHistSecondTOFnsigmavspt", 200, 0., 10., 100, -5. , 5);
  fOutputContainer->Add(fHistSecondTOFnsigmavspt);
  
  fHistSecondTOFnsigmavsp = new TH2F("fHistSecondTOFnsigmavsp", "fHistSecondTOFnsigmavsp", 200, 0., 10., 100, -5., 5.);
  fOutputContainer->Add(fHistSecondTOFnsigmavsp);
  
  fHistSecondTOFsignalvsp = new TH2F("fHistSecondTOFsignalvsp", "fHistSecondTOFsignalvsp", 200, 0., 10, 100,-3000,3000);//1000 , 0., 2.);
  fHistSecondTOFsignalvsp->GetYaxis()->SetTitle("t_{meas}-t_{0}-t_{exoected} (ps)");
  fHistSecondTOFsignalvsp->GetXaxis()->SetTitle("#it{p} (GeV/#it{c})");
  fOutputContainer->Add(fHistSecondTOFsignalvsp);
  
  fHistSecondTOFsignalvspt = new TH2F("fHistSecondTOFsignalvspt", "fHistSecondTOFsignalvspt", 200, 0., 10., 100,-3000,3000);//1000, 0.0, 2.);
  fHistSecondTOFsignalvspt->GetYaxis()->SetTitle("t_{meas}-t_{0}-t_{expected} (ps)");
  fHistSecondTOFsignalvspt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})"); 
  fOutputContainer->Add(fHistSecondTOFsignalvspt);

  fHistSecondTOFTPCsignalvspt = new TH2F("fHistSecondTOFTPCsignalvspt","fHistSecondTOFTPCsignalvspt",200, 0., 10., 100, 0., 20);
  fOutputContainer->Add(fHistSecondTOFTPCsignalvspt);
  
  fHistSecondMultvsCent = new TH2F("fHistSecondMultvsCent","fHistSecondMultvsCent",400,0.,2000.,10,0.,100.);
  fOutputContainer->Add(fHistSecondMultvsCent);

  fHistFirstTPCdEdxAfter = new TH2F("fHistFirstTPCdEdxAfter", "fHistFirstTPCdEdxAfter", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistFirstTPCdEdxAfter);
  
  fHistSecondTPCdEdxAfter = new TH2F("fHistSecondTPCdEdxAfter", "fHistSecondTPCdEdxAfter", 400, -6.0, 6.0, 500, 0.0, 2000);
  fOutputContainer->Add(fHistSecondTPCdEdxAfter);
  
  fHistFirstTOFTPCsignalvsptAfter = new TH2F("fHistFirstTOFTPCsignalvsptAfter","fHistFirstTOFTPCsignalvsptAfter",200, 0., 10., 100, 0., 20);
  fOutputContainer->Add(fHistFirstTOFTPCsignalvsptAfter);
  
  fHistSecondTOFTPCsignalvsptAfter = new TH2F("fHistSecondTOFTPCsignalvsptAfter","fHistSecondTOFTPCsignalvsptAfter",200, 0., 10., 100, 0., 20);
  fOutputContainer->Add(fHistSecondTOFTPCsignalvsptAfter);

  fHistFirstMassTOFvsPt3sTPC       = new TH2F("fHistFirstMassTOFvsPt3sTPC", "fHistFirstMassTOFvsPt3sTPC"  , 200, -2.5, 2.5, 200, -10.0, 10); 
  fHistSecondMassTOFvsPt3sTPC      = new TH2F("fHistSecondMassTOFvsPt3sTPC", "fHistSecondMassTOFvsPt3sTPC", 200, -2.5, 2.5, 200, -10.0, 10);  

  fHistFirstMassTOFvsPt3sTPC3sTOF  = new TH2F("fHistFirstMassTOFvsPt3sTPC3sTOF", "fHistFirstMassTOFvsPt3sTPC3sTOF"  , 200, -2.5 , 2.5, 200, -10.0, 10); 
  fHistSecondMassTOFvsPt3sTPC3sTOF = new TH2F("fHistSecondMassTOFvsPt3sTPC3sTOF", "fHistSecondMassTOFvsPt3sTPC3sTOF", 200, -2.5 , 2.5, 200, -10.0, 10);  
 
  fOutputContainer->Add(fHistFirstMassTOFvsPt3sTPC);
  fOutputContainer->Add(fHistSecondMassTOFvsPt3sTPC);
  fOutputContainer->Add(fHistFirstMassTOFvsPt3sTPC3sTOF);
  fOutputContainer->Add(fHistSecondMassTOFvsPt3sTPC3sTOF);

  //THNSparse
  // UInt_t dimsparse;
  // if(!fReadMCTruth) dimsparse=16;
  // else dimsparse=26;

  fHistTrackBufferOverflow = new TH1F("fHistTrackBufferOverflow","",2,0,2);
  fOutputContainer->Add(fHistTrackBufferOverflow);
  
  fEventCuts.AddQAplotsToList(fOutputContainer);
  
  PostData(1, fOutputContainer );

  //TRee

  if(!fReadMCTruth){
         OpenFile(2);
         fHistSparseSignal = new TTree("fHistSparseSignal","fHistSparseSignal");
	 /*1 */   fHistSparseSignal->Branch("tSignP1",             &tSignP1            , "tSignP1/I" );
	 /*2 */   fHistSparseSignal->Branch("tSignP2",             &tSignP2            , "tSignP2/I" );
	 /*3 */   fHistSparseSignal->Branch("tCentrality",         &tCentrality        , "tCentrality/F" );
	 /*4 */   fHistSparseSignal->Branch("tDCAxyP1",            &tDCAxyP1           , "tDCAxyP1/F" );
	 /*5 */   fHistSparseSignal->Branch("tDCAzP1",             &tDCAzP1            , "tDCAzP1/F" ); 
	 /*6 */   fHistSparseSignal->Branch("tDCAxyP2",            &tDCAxyP2           , "tDCAxyP2/F" ); 
	 /*7 */   fHistSparseSignal->Branch("tDCAzP2",             &tDCAzP2            , "tDCAzP2/F" ); 
	 /*8 */   fHistSparseSignal->Branch("tKtpair",             &tKtpair            , "tKtpair/F" );
	 /*9 */   fHistSparseSignal->Branch("tkStar",              &tkStar             , "tkStar/F" );
	 /*10*/   fHistSparseSignal->Branch("tptP1",               &tptP1              , "tptP1/F" );
	 /*11*/   fHistSparseSignal->Branch("tptP2",               &tptP2              , "tptP2/F" );
	 /*12*/   fHistSparseSignal->Branch("tDEta",               &tDEta              , "tDEta/F" );
	 /*13*/   fHistSparseSignal->Branch("tDPhiStar",           &tDPhiStar          , "tDPhiStar/F" );
	 /*14*/   fHistSparseSignal->Branch("tDPhi",               &tDPhi              , "tDPhi/F" );
	 /*15*/   fHistSparseSignal->Branch("tMassPair",           &tMassPair          , "tMassPair/F" );
	 /*16*/   fHistSparseSignal->Branch("tSphericity",         &tSphericity        , "tSphericity/F" );
	 /*17*/   fHistSparseSignal->Branch("tMassS", &tMassS, "tMassS/F" );
	 /*18*/   fHistSparseSignal->Branch("tDTheta",             &tDTheta            , "tDTheta/F" );
	 fHistSparseSignal->SetAutoSave(100000000);
	 
	 PostData(2, fHistSparseSignal );

         OpenFile(3);
	 fHistSparseBkg = new TTree("fHistSparseBkg","fHistSparseBkg");
	 /*1 */   fHistSparseBkg->Branch("tSignP1",             &tSignP1            , "tSignP1/I" );
	 /*2 */   fHistSparseBkg->Branch("tSignP2",             &tSignP2            , "tSignP2/I" );
	 /*3 */   fHistSparseBkg->Branch("tCentrality",         &tCentrality        , "tCentrality/F" );
	 /*4 */   fHistSparseBkg->Branch("tDCAxyP1",            &tDCAxyP1           , "tDCAxyP1/F" );
	 /*5 */   fHistSparseBkg->Branch("tDCAzP1",             &tDCAzP1            , "tDCAzP1/F" ); 
	 /*6 */   fHistSparseBkg->Branch("tDCAxyP2",            &tDCAxyP2           , "tDCAxyP2/F" ); 
	 /*7 */   fHistSparseBkg->Branch("tDCAzP2",             &tDCAzP2            , "tDCAzP2/F" ); 
	 /*8 */   fHistSparseBkg->Branch("tKtpair",             &tKtpair            , "tKtpair/F" );
	 /*9 */   fHistSparseBkg->Branch("tkStar",              &tkStar             , "tkStar/F" );
	 /*10*/   fHistSparseBkg->Branch("tptP1",               &tptP1              , "tptP1/F" );
	 /*11*/   fHistSparseBkg->Branch("tptP2",               &tptP2              , "tptP2/F" );
	 /*12*/   fHistSparseBkg->Branch("tDEta",               &tDEta              , "tDEta/F" );
	 /*13*/   fHistSparseBkg->Branch("tDPhiStar",           &tDPhiStar          , "tDPhiStar/F" );
	 /*14*/   fHistSparseBkg->Branch("tDPhi",               &tDPhi              , "tDPhi/F" );
	 /*15*/   fHistSparseBkg->Branch("tMassPair",           &tMassPair          , "tMassPair/F" );
	 /*16*/   fHistSparseBkg->Branch("tSphericity",         &tSphericity        , "tSphericity/F" );
	 /*17*/   fHistSparseBkg->Branch("tMassS", &tMassS, "tMassS/F" );
	 /*18*/   fHistSparseBkg->Branch("tDTheta",             &tDTheta            , "tDTheta/F" );
	 fHistSparseBkg->SetAutoSave(100000000);
	 PostData(3, fHistSparseBkg );
  }

  

  else{
    OpenFile(2);
    fHistSparseSignal = new TTree("fHistSparseSignal","fHistSparseSignal");
    /*1 */   fHistSparseSignal->Branch("tSignP1",             &tSignP1            , "tSignP1/I" );
    /*2 */   fHistSparseSignal->Branch("tSignP2",             &tSignP2            , "tSignP2/I" );
    /*3 */   fHistSparseSignal->Branch("tCentrality",         &tCentrality        , "tCentrality/F" );
    /*4 */   fHistSparseSignal->Branch("tDCAxyP1",            &tDCAxyP1           , "tDCAxyP1/F" );
    /*5 */   fHistSparseSignal->Branch("tDCAzP1",             &tDCAzP1            , "tDCAzP1/F" ); 
    /*6 */   fHistSparseSignal->Branch("tDCAxyP2",            &tDCAxyP2           , "tDCAxyP2/F" ); 
    /*7 */   fHistSparseSignal->Branch("tDCAzP2",             &tDCAzP2            , "tDCAzP2/F" ); 
    /*8 */   fHistSparseSignal->Branch("tKtpair",             &tKtpair            , "tKtpair/F" );
    /*9 */   fHistSparseSignal->Branch("tkStar",              &tkStar             , "tkStar/F" );
    /*10*/   fHistSparseSignal->Branch("tptP1",               &tptP1              , "tptP1/F" );
    /*11*/   fHistSparseSignal->Branch("tptP2",               &tptP2              , "tptP2/F" );
    /*12*/   fHistSparseSignal->Branch("tDEta",               &tDEta              , "tDEta/F" );
    /*13*/   fHistSparseSignal->Branch("tDPhiStar",           &tDPhiStar          , "tDPhiStar/F" );
    /*14*/   fHistSparseSignal->Branch("tDPhi",               &tDPhi              , "tDPhi/F" );
    /*15*/   fHistSparseSignal->Branch("tMassPair",           &tMassPair          , "tMassPair/F" );
    /*16*/   fHistSparseSignal->Branch("tSphericity",         &tSphericity        , "tSphericity/F" );
    /*17*/   fHistSparseSignal->Branch("tMassS", &tMassS, "tMassS/F" );
    /*18*/   fHistSparseSignal->Branch("tDTheta",             &tDTheta            , "tDTheta/F" );
    /*19*/   fHistSparseSignal->Branch("tMCtruepair",         &tMCtruepair        , "tMCtruepair/I" );
    /*20*/   fHistSparseSignal->Branch("tMCSameMother",       &tMCSameMother      , "tMCSameMother/I" );
    /*21*/   fHistSparseSignal->Branch("tMCMotherP1",         &tMCMotherP1        , "tMCMotherP1/I" );
    /*22*/   fHistSparseSignal->Branch("tMCMotherP2",         &tMCMotherP2        , "tMCMotherP2/I" );
    /*23*/   fHistSparseSignal->Branch("tMCptcTypeP1",        &tMCptcTypeP1       , "tMCptcTypeP1/I" ); 
    /*24*/   fHistSparseSignal->Branch("tMCptcTypeP2",        &tMCptcTypeP2       , "tMCptcTypeP2/I" ); 
    /*25*/   fHistSparseSignal->Branch("tMCSameGM",           &tMCSameGM          , "tMCSameGM/I" );
    /*26*/   fHistSparseSignal->Branch("tMotherPDG",          &tMotherPDG         , "tMotherPDG/I" );
    /*27*/   fHistSparseSignal->Branch("tpdgcodeP1",          &tpdgcodeP1         , "tpdgcodeP1/I" ); 
    /*28*/   fHistSparseSignal->Branch("tpdgcodeP2",          &tpdgcodeP2         , "tpdgcodeP2/I" );
    /*29*/   fHistSparseSignal->Branch("tKstarGen",           &tKstarGen          , "tKstarGen/F" );
    fHistSparseSignal->SetAutoSave(100000000);
    PostData(2, fHistSparseSignal );
    
    OpenFile(3);
    fHistSparseBkg = new TTree("fHistSparseBkg","fHistSparseBkg");
    /*1 */   fHistSparseBkg->Branch("tSignP1",             &tSignP1            , "tSignP1/I" );
    /*2 */   fHistSparseBkg->Branch("tSignP2",             &tSignP2            , "tSignP2/I" );
    /*3 */   fHistSparseBkg->Branch("tCentrality",         &tCentrality        , "tCentrality/F" );
    /*4 */   fHistSparseBkg->Branch("tDCAxyP1",            &tDCAxyP1           , "tDCAxyP1/F" );
    /*5 */   fHistSparseBkg->Branch("tDCAzP1",             &tDCAzP1            , "tDCAzP1/F" ); 
    /*6 */   fHistSparseBkg->Branch("tDCAxyP2",            &tDCAxyP2           , "tDCAxyP2/F" ); 
    /*7 */   fHistSparseBkg->Branch("tDCAzP2",             &tDCAzP2            , "tDCAzP2/F" ); 
    /*8 */   fHistSparseBkg->Branch("tKtpair",             &tKtpair            , "tKtpair/F" );
    /*9 */   fHistSparseBkg->Branch("tkStar",              &tkStar             , "tkStar/F" );
    /*10*/   fHistSparseBkg->Branch("tptP1",               &tptP1              , "tptP1/F" );
    /*11*/   fHistSparseBkg->Branch("tptP2",               &tptP2              , "tptP2/F" );
    /*12*/   fHistSparseBkg->Branch("tDEta",               &tDEta              , "tDEta/F" );
    /*13*/   fHistSparseBkg->Branch("tDPhiStar",           &tDPhiStar          , "tDPhiStar/F" );
    /*14*/   fHistSparseBkg->Branch("tDPhi",               &tDPhi              , "tDPhi/F" );
    /*15*/   fHistSparseBkg->Branch("tMassPair",           &tMassPair          , "tMassPair/F" );
    /*16*/   fHistSparseBkg->Branch("tSphericity",         &tSphericity        , "tSphericity/F" );
    /*17*/   fHistSparseBkg->Branch("tMassS", &tMassS, "tMassS/F" );
    /*18*/   fHistSparseBkg->Branch("tDTheta",             &tDTheta            , "tDTheta/F" );
    /*19*/   fHistSparseBkg->Branch("tMCtruepair",         &tMCtruepair        , "tMCtruepair/I" );
    /*20*/   fHistSparseBkg->Branch("tMCSameMother",       &tMCSameMother      , "tMCSameMother/I" ); //this was missing
    /*21*/   fHistSparseBkg->Branch("tMCMotherP1",         &tMCMotherP1        , "tMCMotherP1/I" );
    /*22*/   fHistSparseBkg->Branch("tMCMotherP2",         &tMCMotherP2        , "tMCMotherP2/I" );
    /*23*/   fHistSparseBkg->Branch("tMCptcTypeP1",        &tMCptcTypeP1       , "tMCptcTypeP1/I" ); 
    /*24*/   fHistSparseBkg->Branch("tMCptcTypeP2",        &tMCptcTypeP2       , "tMCptcTypeP2/I" ); 
    /*25*/   fHistSparseBkg->Branch("tMCSameGM",           &tMCSameGM          , "tMCSameGM/I" );
    /*26*/   fHistSparseBkg->Branch("tMotherPDG",          &tMotherPDG         , "tMotherPDG/I" );
    /*27*/   fHistSparseBkg->Branch("tpdgcodeP1",          &tpdgcodeP1         , "tpdgcodeP1/I" ); 
    /*28*/   fHistSparseBkg->Branch("tpdgcodeP2",          &tpdgcodeP2         , "tpdgcodeP2/I" );
    /*29*/   fHistSparseBkg->Branch("tKstarGen",           &tKstarGen          , "tKstarGen/F" );
    fHistSparseBkg->SetAutoSave(100000000);
    PostData(3, fHistSparseBkg ); 
  }
  //TTree *fHistSparseBkg;      
  
 
  
}

//---------------------------------------------------------------------------------

void AliAnalysisTaskDeuFlow2PC::UserExec(Option_t *) {
  
  //  cout<<"----------------------------->Enter user exec"<<endl;

  // Main loop
  // Called for each event

  
  AliVVertex *vertexmain =0x0;
  // AliCentrality* centrality = 0x0;
  AliMultSelection* centrality = 0x0;
  Double_t lBestPrimaryVtxPos[3] = {-100.0, -100.0, -100.0};

  fHistEventMultiplicity->Fill(1);
  AliMCEvent   *lMCevent  = 0x0;
  AliStack     *lMCstack  = 0x0;
  TClonesArray *arrayMC = 0x0;

  Int_t ntracks = 0;

  if (fAnalysisType == "ESD") {
    fESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!fESDevent) {
      AliWarning("ERROR: ESDevent not available \n");
      PostData(1, fOutputContainer);
      PostData(2, fHistSparseSignal );
      PostData(3, fHistSparseBkg );
      return;
    }
 
    /// Use the event cut class to apply the required selections
    if (!fEventCuts.AcceptEvent(fESDevent)) {
      return;
    }
    
    if (fReadMCTruth) {
      lMCevent = MCEvent();
      if (!lMCevent) {
        Printf("ERROR: Could not retrieve MC event \n");
        //cout << "Name of the file with pb :" <<  CurrentFileName() << endl; 
	PostData(1, fOutputContainer);
	PostData(2, fHistSparseSignal );
	PostData(3, fHistSparseBkg );
        return;
      }
      
      lMCstack = lMCevent->Stack();
      // cout<<"Stack: "<<lMCstack<<endl;
      if (!lMCstack) {
        Printf("ERROR: Could not retrieve MC stack \n");
        //cout << "Name of the file with pb :" <<  CurrentFileName() << endl;
	PostData(1, fOutputContainer);
	PostData(2, fHistSparseSignal );
	PostData(3, fHistSparseBkg );
      return;
      }
    }
   
     
    ntracks = fESDevent->GetNumberOfTracks(); 
    //    centrality = fESDevent->GetCentrality(); //FIXME : Find out centrality object in pp and pPb
    
    const AliESDVertex *lPrimaryBestESDVtx = fESDevent->GetPrimaryVertex();
    if (!lPrimaryBestESDVtx){
      AliWarning("No prim. vertex in ESD... return!");
      PostData(1, fOutputContainer);
      PostData(2, fHistSparseSignal );
      PostData(3, fHistSparseBkg );

      return;
    }
    vertexmain = (AliVVertex*) lPrimaryBestESDVtx;
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
    
  } else if (fAnalysisType == "AOD") {
    
    fAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
    //    cout<<"fAODevent: "<<fAODevent <<endl;
    
    if (!fAODevent) {
      AliWarning("ERROR: AODevent not available \n"); 
      PostData(1, fOutputContainer);
      PostData(2, fHistSparseSignal );
      PostData(3, fHistSparseBkg );
      return;
    }
    /// Use the event cut class to apply the required selections
    if (!fEventCuts.AcceptEvent(fAODevent)) {   
      PostData(1, fOutputContainer);
      PostData(2, fHistSparseSignal );
      PostData(3, fHistSparseBkg );
      return;
    }
    ntracks = fAODevent->GetNumberOfTracks();
    //    centrality = fAODevent->GetCentrality(); //FIXME :  Find out centrality object in pp and pPb
   
  
    const AliAODVertex *lPrimaryBestAODVtx = fAODevent->GetPrimaryVertex();
    if (!lPrimaryBestAODVtx){
      AliWarning("No prim. vertex in AOD... return!");
      PostData(1, fOutputContainer);
      PostData(2, fHistSparseSignal );
      PostData(3, fHistSparseBkg );
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
    PostData(1, fOutputContainer);
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;

  }
  

  // PID object
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  
  UInt_t mask = inputHandler->IsEventSelected();
  /* to see how many events are rejected
     if (!(mask & 0xffffffff)) {
     PostData(1,fOutputContainer );
     return;
     }
  */
  fPIDResponse = inputHandler->GetPIDResponse();
  
  if(!fPIDResponse) {
    PostData(1,fOutputContainer );
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    AliError("Cannot get pid response");
    return;
  }

  fHistEventMultiplicity->Fill(2);

  if((TMath::Abs(lBestPrimaryVtxPos[2])) > 10.) {
    PostData(1, fOutputContainer);
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;
  }
  fHistEventMultiplicity->Fill(3);

  //event must not be tagged as pileup
  Bool_t isPileUpSpd=kFALSE;
  isPileUpSpd=fAODevent->IsPileupFromSPD();
  if(isPileUpSpd){ 
    PostData(1,fOutputContainer );
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;
  }
  
  fHistEventMultiplicity->Fill(4);



  Float_t lcentrality = -99.;
  
  //  AliMultSelection* centrality = 0x0;
  centrality = (AliMultSelection *) fAODevent->FindListObject("MultSelection");
  //  cout<<"centrality: "<<centrality<<endl;

  if(fCollidingSystem == "PbPb" ){
    //    lcentrality = centrality->GetCentralityPercentile("V0M"); //FIXME : Centrality in pp and pPb works the same???
    lcentrality = centrality->GetMultiplicityPercentile("V0M"); //FIXME : Centrality in pp and pPb works the same???
    hmult->Fill(lcentrality);
  }
  
  else if(fCollidingSystem == "pPb"){
    //    lcentrality = ((AliAODHeader * )fAODevent->GetHeader())->GetRefMultiplicityComb08(); //-->TBChecked
    lcentrality = centrality->GetMultiplicityPercentile("V0A"); //FIXME : Now we use Centrality instead of multiplicity? also for p-Pb
    //    hmult->Fill(lcentrality);
  } 
  
  else if(fCollidingSystem == "pp") {//FIXME : I think up AOD have only refmult as mult estimation
    // lcentrality = ((AliAODHeader * )fAODevent->GetHeader())->GetRefMultiplicityComb08(); //-->RIcambiare
    //    lcentrality = ((AliAODHeader * )fAODevent->GetHeader())->GetRefMultiplicity(); //run on phojet
    lcentrality = centrality->GetMultiplicityPercentile("V0M"); //FIXME : Also for pp? Test on kd
    //    cout<<"Centrality: "<<lcentrality<<endl;
    //    hmult->Fill(lcentrality);
    
  }

  //  cout<<"centrality: "<<lcentrality<<endl;

  if ( lcentrality > 199 ){
    //Event didn't pass Event Selections
    PostData(1,fOutputContainer );
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;
  }
  
  if (lcentrality<fCentrLowLim||lcentrality>=fCentrUpLim){   
    PostData(1,fOutputContainer );
    PostData(2, fHistSparseSignal );
    PostData(3, fHistSparseBkg );
    return;
  }
  fHistEventMultiplicity->Fill(5);
  

  Bool_t isSelectedCentral     = kFALSE;
  Bool_t isSelectedSemiCentral = kFALSE;
  Bool_t isSelectedMB          = kFALSE;
  Bool_t isSelectedInt7        = kFALSE;
  Bool_t isSelectedHM          = kFALSE;
  Bool_t isSelectedAny         = kFALSE;
  Bool_t isSelected            = kFALSE;
 
  if(fCollidingSystem == "PbPb"){
    isSelectedCentral     = (mask & AliVEvent::kCentral);
    isSelectedSemiCentral = (mask & AliVEvent::kSemiCentral);
    isSelectedMB          = (mask & AliVEvent::kMB);
    isSelectedAny         = (mask & AliVEvent::kAny);
  }
  
  else if(fCollidingSystem == "pp"){
  
    isSelectedCentral     = (mask & AliVEvent::kCentral);
    isSelectedSemiCentral = (mask & AliVEvent::kSemiCentral);
    isSelectedMB          = (mask & AliVEvent::kMB);
    isSelectedInt7        = (mask & AliVEvent::kINT7);
    isSelectedHM          = (mask & AliVEvent::kHighMultV0);
    isSelectedAny         = (mask & AliVEvent::kAnyINT);
    
    // if(isSelectedHM)
    // cout<<"isSelectedCentral    : "<<isSelectedCentral     <<endl;
    // cout<<"isSelectedSemiCentral: "<<isSelectedSemiCentral <<endl;
    // cout<<"isSelectedMB         : "<<isSelectedMB          <<endl;
    // cout<<"isSelectedInt7       : "<<isSelectedInt7        <<endl;
    // cout<<"isSelectedHM         : "<<isSelectedHM          <<endl;
    // cout<<"isSelectedAny        : "<<isSelectedAny         <<endl;
      
    if(fYear == 2010 && isSelectedMB )
      isSelected = kTRUE;
    else if(fYear != 2010 && fHMtrigger == kFALSE && isSelectedInt7)
      isSelected = kTRUE;
    else if(fYear != 2010 && fHMtrigger == kTRUE && isSelectedHM)
      isSelected = kTRUE;
    else 
      isSelected = kFALSE;
  }
  
  else if(fCollidingSystem == "pPb"){
    isSelectedCentral     = (mask & AliVEvent::kCentral);
    isSelectedSemiCentral = (mask & AliVEvent::kSemiCentral);
    isSelectedMB          = (mask & AliVEvent::kMB);
    isSelectedInt7        = (mask & AliVEvent::kINT7);
    isSelectedAny         = (mask & AliVEvent::kAny);

    if(isSelectedInt7 )
      isSelected = kTRUE;
  }

  if(isSelectedAny)
    fHistEventMultiplicity->Fill(6);
  if(isSelectedCentral)
    fHistEventMultiplicity->Fill(7);
  if(isSelectedSemiCentral)
    fHistEventMultiplicity->Fill(8);
  if(isSelectedMB)
    fHistEventMultiplicity->Fill(9);
  if(isSelectedInt7)
    fHistEventMultiplicity->Fill(10);
  if(isSelectedHM)
    fHistEventMultiplicity->Fill(11);
  
  //  cout<<"Trigger mask: "<<fAODevent->GetTriggerMask()<<" "<<AliVEvent::kHighMultV0<<endl;
  //  cout<<"Trigger mask: "<<mask<<" "<<AliVEvent::kHighMultV0<<endl;
  // cout<<"Event type  : "<<fAODevent->GetEventType()<<endl;
  // FIXME : event selection to be added.. DONE

  
  if(fCollidingSystem == "PbPb"){
    if(!isSelectedCentral && !isSelectedSemiCentral && !isSelectedMB)  {
      PostData(1, fOutputContainer);
      PostData(2, fHistSparseSignal );
      PostData(3, fHistSparseBkg );
      return;
    }
  }
  
  else  if(fCollidingSystem == "pp"){ // check pPb
    
    if(!isSelected){   
      PostData(1, fOutputContainer);
      PostData(2, fHistSparseSignal );
      PostData(3, fHistSparseBkg );
      return;
    }
  }  
  
  else  if(fCollidingSystem == "pPb" ){ // check pPb
    
    if(!isSelectedInt7){   
    //if(!isSelected){   
      PostData(1, fOutputContainer);
      PostData(2, fHistSparseSignal );
      PostData(3, fHistSparseBkg );
      return;
    }
    
  }

  fHistEventMultiplicity->Fill(12); // is event selected for the analysis
  
  if(fCollidingSystem == "pp" || fCollidingSystem == "pPb" ) 
    hmult->Fill(((AliAODHeader * )fAODevent->GetHeader())->GetRefMultiplicityComb08());

  // cout<<"nTracks: "<<ntracks<<" centrality "<<lcentrality<<endl;
  Double_t fSphericityvalue = CalculateSphericityofEvent(fAODevent);
  fHistSphericity->Fill(fSphericityvalue);
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
 
  int spherBin=0;
  double spherStep=1/double(fnSphericityBins ), spherStart=-0.;
  
  for (int i=0; i<fnSphericityBins ; i++) {
    if ((fSphericityvalue > spherStart+i*spherStep) && (fSphericityvalue < spherStart+(i+1)*spherStep)) {
      spherBin=i;
      break;
    }
  }
  //CENTRALITY!!
  // ... and centrality    //FIXME : find out how centrality wokrs in AOD in pp and pPb
  int centralityBin=0;
  if(fCollidingSystem == "PbPb" ){
    if(lcentrality < 5.) centralityBin=19;      // changed <= with < to be consistent with histogram binning, except last bin 
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
  }

  else if(fCollidingSystem == "pp"  || fCollidingSystem == "pPb"){                // FIXME : To be understand if these selections are fine in the current framework
    /* // This is valid for multiplicity
    if(lcentrality < 4.)        centralityBin=19;  
    else if(lcentrality < 7.)   centralityBin=18;
    else if(lcentrality < 10.)  centralityBin=17;
    else if(lcentrality < 15.)  centralityBin=16;
    else if(lcentrality < 20.)  centralityBin=15;
    else if(lcentrality < 25.)  centralityBin=14;
    else if(lcentrality < 30.)  centralityBin=13;
    else if(lcentrality < 40.)  centralityBin=12;
    else if(lcentrality < 50.)  centralityBin=11;
    else if(lcentrality < 60.)  centralityBin=10;
    else if(lcentrality < 70.)  centralityBin=9; 
    else if(lcentrality < 80.)  centralityBin=8;
    else if(lcentrality < 90.)  centralityBin=7;
    else if(lcentrality < 100.) centralityBin=6;
    else if(lcentrality < 150.) centralityBin=5;
    else if(lcentrality < 200.) centralityBin=4;
    else if(lcentrality < 250.) centralityBin=3;
    else if(lcentrality < 300.) centralityBin=2;   
    else if(lcentrality < 350.) centralityBin=1;
    else if(lcentrality <= 350.)centralityBin=0;
    */
    // this should be valid for centrality...
    if(lcentrality < 0.01)       centralityBin=19;  
    else if(lcentrality < 0.1)   centralityBin=18;
    else if(lcentrality < 0.5)   centralityBin=17;
    else if(lcentrality < 1.0)   centralityBin=16;
    else if(lcentrality < 5.0)   centralityBin=15;
    else if(lcentrality < 10.)   centralityBin=14;
    else if(lcentrality < 20.)   centralityBin=13;
    else if(lcentrality < 30.)   centralityBin=12;
    else if(lcentrality < 40.)   centralityBin=11;
    else if(lcentrality < 50.)   centralityBin=10;
    else if(lcentrality < 70.)   centralityBin=9; 
    else if(lcentrality <= 100.) centralityBin=8;
  }
  
  if(fkDoSphericity == kFALSE){
    fEventColl[zBin][centralityBin]->FifoShift();
    fEvt = fEventColl[zBin][centralityBin]->fEvt;
  }
  
  else{
    fEventCollwSp[zBin][centralityBin][spherBin]->FifoShift();
    fEvt = fEventCollwSp[zBin][centralityBin][spherBin]->fEvt;
  }
  
  //  printf("buffer size: %d\n",fTrackBufferSize);
  //  printf("ntracks: %d\n",ntracks);

  ULong_t  status=0;

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

  // Read, plot and store first and second tracks 
  
  AliAODTrack *track = 0x0;
  AliVTrack *vtrackg = 0x0;
  AliVTrack *vtrack = 0x0;

  Bool_t isTOFPIDok = kFALSE;
  //  Float_t nsigmaTOFf = 10.;
  //  Float_t nsigmaTPCf = 10.;
 
  Float_t nsigmaTOFs = 10.;
  Float_t nsigmaTPCs = 10.;
  
  // Float_t pullsTPC[5] = {0.};
  // Float_t pullsTOF[5] = {10.};

  Float_t probMis = 0.;
  Double32_t tTOF = 0.;

  Float_t beta  = 0.;
  Float_t gamma = 0.;
  Float_t mass  = 0.;
  Float_t length = 0;
  Float_t ptot = 0;

  Float_t nTPCCrossedRows = 0.; 
  Float_t rationCrnFind = 0.;

  UShort_t nTPCclSp = 0;
  UShort_t nTPCclp = 0;
 
  Int_t label = 0;
  Int_t PDGcode=0;

  Float_t sharedFractionTPCcl = 0.;
  //  AliPIDResponse::EDetPidStatus statusTOF;
  //  AliPIDResponse::EDetPidStatus statusTPC;
 
  Double_t expectedTimes[AliPID::kSPECIES];

  Double_t dz[2] = {-999.,-999.}; //Double_t covar[3]={-999.,-999.,-999.};
  Double_t dzg[2]= {-999.,-999.}; Double_t covarg[3]={-999.,-999.,-999.};
  AliExternalTrackParam etp1;

  Float_t rapidityFirst = 0.;
  Float_t rapiditySecond = 0.;
  Short_t charge = -2;

  int fCount = 0;
  int sCount = 0;

  // Float_t pMomentumTruth[3];
  AliReconstructed2pcFirst::MCFirstOrigin_t mcFirstOrigin = AliReconstructed2pcFirst::kUnassigned;
  AliReconstructed2pcSecond::MCSecondOrigin_t mcSecondOrigin = AliReconstructed2pcSecond::kUnassigned;
  
  Bool_t isP = kFALSE;  // particle
  Bool_t isaP = kFALSE; // anti-particle 

  Bool_t isMCfirst  = kFALSE; 
  Bool_t isMCsecond = kFALSE; 
  Int_t MCptcCodeP1 = -999; //0 : primary 2: from weak decay 3: from material
  Int_t MCptcCodeP2 = -999;
  Int_t MCmumIDP1   = -999;
  Int_t MCmumPDGP1  = -999;
  Int_t MCmumIDP2   = -999;
  Int_t MCmumPDGP2  = -999;
  Int_t MCgrammaPDGP1  = -999;
  Int_t MCgrammaPDGP2  = -999;
  Int_t MCgrammaIDP1   = -999;
  Int_t MCgrammaIDP2   = -999;

  // cout<<"Min pt prim: "<<fMinPtForPrim<<" max pt: "<<fMaxPtForPrim<<endl;
  // cout<<"Min pt sec : "<<fMinPtForSec<<"  max pt: "<<fMaxPtForSec<<endl;
  
  for (Int_t ip = 0; ip < ntracks; ip++){
    
    vtrack = fAODevent->GetTrack(ip);
    if (!vtrack) continue;
    track = dynamic_cast<AliAODTrack*>(vtrack);
    if(!track) AliFatal("Not a standard AOD");

    // Take TPC only tracks constrained to the SPD PV during filtering, Filter Bit 7 (128) but PID is not stored --> find the corresponding global track filter bit 1
    //if(!track->TestFilterBit(AliAODTrack::kTrkTPCOnlyConstrained)) continue;  
    // Bit 7 corresponds to the following cuts (applied to original ESD track, see Ruben presentation 13.02.2014) : 
    // ncl 50  
    // dcaxy 2.4 dcaz 3.2 // this cut however is on the global track!! PWGCF meeting (17.02.2014)
    // no kinks 
    // chi2percl 4 
    // dcatovtx2D kTRUE

    if(!track->TestFilterBit(fFilterBit))
      continue;

    // nTPCCrossedRows = track->GetTPCClusterInfo(2, 1); // The method below is equivalent
    nTPCCrossedRows = track->GetTPCNCrossedRows();
    //    cout<<"Crossed rows dir "<<track->GetTPCNCrossedRows()<<" Crossed rows indir "<<nTPCCrossedRows<<endl;
    //    if(nTPCCrossedRows<70) continue;  
    
    rationCrnFind = 0.;
    if(track->GetTPCNclsF()) rationCrnFind = nTPCCrossedRows/track->GetTPCNclsF();   
    if(fkApplyRatioCrRnFindCut) if(rationCrnFind< .8) continue;  
    
    nTPCclSp = track->GetTPCnclsS();
    nTPCclp  = track->GetTPCncls();
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
    if(fFilterBit == 128)
      vtrackg = fAODevent->GetTrack(farrGT[-track->GetID()-1]);
    else
      vtrackg = track;

    if (!vtrackg) {
      printf ("No global info! iTrack %d, ID %d\n",ip,track->GetID());
      continue;
    }
    if (farrGT[-track->GetID()-1]>=ntracks || farrGT[-track->GetID()-1]<0) { /*cout<<"This index is out of range!!"<<farrGT[-track->GetID()-1]<<endl;*/ continue;}
    globaltrack = dynamic_cast<AliAODTrack*>(vtrackg); 
    if(!globaltrack) AliFatal("Not a standard AOD");

    //    cout<<" Filter map for the global track "<<globaltrack->GetFilterMap()<<" "<<globaltrack<<endl;
     
    // IP to PV of tracks
    dz[0] = -999.; dz[1]  = -999.;
    dzg[0]= -999.; dzg[1] = -999.;

    dz[0] = track->DCA();    // the TPC one should be applied the other biases the CF --> from Maciejs note --> FIXME to be checked 
    dz[1] = track->ZAtDCA(); // for those two lines check AliAODTrack.h // FIXME these two lines produce shifted distributions, known problem, asked Mac and Marian. 
    
    if(fFilterBit == 4){//rl
      Double_t covd[3];
      AliAODTrack* track_clone=(AliAODTrack*)globaltrack->Clone("track_clone"); // need to clone because PropagateToDCA updates the track parameters
      Bool_t isDCA = track_clone->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),9999.,dz,covd);
      delete track_clone;
      if(!isDCA)dz[0]=-999.;
    }
    // Btw dont propagate TPC constrained! 
    //    Double_t p[3]; track->GetXYZ(p); // first two elements of p give the same as above, ignore the third, those are original DCA of TPC only tracks (with or w/o constraint?)
    //    cout<<"d xy "<<dz[0]<<"while value is for other "<<p[0]<<endl;
    //    cout<<"d z "<<dz[1]<<"while value is for other "<<p[1]<<endl;
    
    etp1.CopyFromVTrack(vtrackg); 
    etp1.PropagateToDCA(vertexmain,(InputEvent())->GetMagneticField(),100.,dzg,covarg);    

    isP  = kFALSE;
    isaP = kFALSE;

    charge = globaltrack->Charge();

    if (charge>0) {
      isP  = kTRUE;
      isaP = kFALSE;
    }
    else{
      if(charge<0){
	isP  = kFALSE;
	isaP = kTRUE;
      }
    }
    
    if (TMath::Abs(track->Eta())> 0.8) continue;
    if (track->Chi2perNDF()>4.) continue; // should be redundant already applied in the TPC only track ESD filtering

    // if (track->Pt()<fMinPtForPrim || track->Pt()> fMaxPtForPrim) continue;  
   
    //    nsigmaTOFf = 10.; // be careful with those initialization
    //    nsigmaTPCf = 10.;
    
    nsigmaTOFs = 10.; // be careful with those initialization
    nsigmaTPCs = 10.;
    
    probMis = 10.;
 
    //    statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,globaltrack);
    //    cout<<"statusTPC: "<<statusTPC <<" "<<AliPIDResponse::kDetPidOk<<endl;
    //  fHistTPCdEdx->Fill(globaltrack->GetTPCmomentum()*charge, globaltrack->GetTPCsignal());
     
    //FIXME : controlla questo if
    // if (AliPIDResponse::kDetPidOk != statusTPC) { 
    //   //printf ("TPC PID not there"); 
    //   continue; 
    // }
    
    rapidityFirst  = 0; 
    rapiditySecond = 0.5*TMath::Log( (track->E(fPDGsecond) + track->Pz()) / (track->E(fPDGsecond) - track->Pz() +1.e-13));
    
    if(fFirstpart  == 5){
      Double_t energyDeuteron = TMath::Sqrt(track->GetP()*track->GetP() + 4*AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton));
      rapidityFirst = 0.5*TMath::Log((energyDeuteron + track->Pz())/(energyDeuteron - track->Pz()));
    }			     
    if(fSecondpart == 5){
      Double_t energyDeuteron = TMath::Sqrt(track->GetP()*track->GetP() + 4*AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton));
      rapiditySecond = 0.5*TMath::Log((energyDeuteron + track->Pz())/(energyDeuteron - track->Pz()));
    }
    
    //cout<<"Firt particle: "<<fFirstpart<<" "<<rapidityFirst<<" second particle: "<<fSecondpart<<" "<<rapiditySecond<<endl;

    //  if(fFirstpart == kKaon){

    //    nsigmaTPCf = fPIDResponse->NumberOfSigmasTPC(globaltrack, (AliPID::EParticleType)fFirstpart); 
    nsigmaTPCs = fPIDResponse->NumberOfSigmasTPC(globaltrack, (AliPID::EParticleType)fSecondpart); 
    
    //cout<<"Firt particle ns: "<<nsigmaTPCf<<" "<<nsigmaTPCs<<endl;
   
    // //pulls TPC
    // pullsTPC[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(globaltrack, (AliPID::EParticleType)kElectron)); 
    // pullsTPC[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(globaltrack, (AliPID::EParticleType)kPion)); 
    // pullsTPC[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(globaltrack, (AliPID::EParticleType)kKaon)); 
    // pullsTPC[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(globaltrack, (AliPID::EParticleType)kProton)); 
    // pullsTPC[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(globaltrack, (AliPID::EParticleType)kDeuteron)); 
    
    //    nsigmaTPCele = fPIDResponse->NumberOfSigmasTPC(globaltrack, (AliPID::EParticleType)kElectron); 
    //Electron veto
    // if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(globaltrack,(AliPID::EParticleType)kElectron)  <1 ) && globaltrack->GetTPCmomentum() < 0.45/*&&
    //    TMath::Abs(fPIDResponse->NumberOfSigmasTPC(globaltrack, (AliPID::EParticleType)kPion)>3)   &&
    //    TMath::Abs(fPIDResponse->NumberOfSigmasTPC(globaltrack, (AliPID::EParticleType)kKaon)>3)   &&
    //    TMath::Abs(fPIDResponse->NumberOfSigmasTPC(globaltrack, (AliPID::EParticleType)kProton)>3) */)continue;

    //TOF
    
    //    statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,globaltrack); // this checks kTOFout and kTIMEi https://twiki.cern.ch/twiki/bin/viewauth/ALICE/TOF

    tTOF = 0.0;
    isTOFPIDok = kFALSE;

    
    // if ( (statusTOF ==  AliPIDResponse::kDetPidOk) ) {  
    //nsigmaTOFf = fPIDResponse->NumberOfSigmasTOF(globaltrack,  (AliPID::EParticleType)fFirstpart); 
      

      // //pulls TOF
      // pullsTOF[0] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(globaltrack, (AliPID::EParticleType)kElectron)); 
      // pullsTOF[1] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(globaltrack, (AliPID::EParticleType)kPion)); 
      // pullsTOF[2] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(globaltrack, (AliPID::EParticleType)kKaon)); 
      // pullsTOF[3] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(globaltrack, (AliPID::EParticleType)kProton)); 
      // pullsTOF[4] = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(globaltrack, (AliPID::EParticleType)kDeuteron)); 
      //------------------------------
      
      //isTOFPIDok = kTRUE;
    status  = (ULong_t)globaltrack->GetStatus();
    Bool_t hasTOFout  = status&AliVTrack::kTOFout; 
    if (hasTOFout) isTOFPIDok = kTRUE;
    length = globaltrack->GetIntegratedLength(); 
    if (length < 350.) isTOFPIDok = kFALSE;
   
    if(isTOFPIDok){
      
      probMis = fPIDResponse->GetTOFMismatchProbability(globaltrack);
      tTOF = globaltrack->GetTOFsignal(); 
      tTOF -= fPIDResponse->GetTOFResponse().GetStartTime(globaltrack->P());      
      globaltrack->GetIntegratedTimes(expectedTimes);
      nsigmaTOFs = fPIDResponse->NumberOfSigmasTOF(globaltrack,  (AliPID::EParticleType)fSecondpart); 
    }
    else { 
      probMis = 1.; 
      //nsigmaTOFf = 10.;
      nsigmaTOFs = 10.; //cout<<"The corresponding global track has no tof pid!"<<endl;
    } 
    
    
    fHistTPCdEdx->Fill(globaltrack->GetTPCmomentum()*charge, globaltrack->GetTPCsignal());
   
    //
    // HERE ARE THE LOOPs!
    //
    
    AliAODMCParticle *tparticle = 0x0;
    
    if(fSecondpart == kAny) fnSigmaTPCPIDsecondParticle = 0;

    //if (TMath::Abs(nsigmaTPCf)> fnSigmaTPCPIDsecondParticle) {
    if (TMath::Abs(nsigmaTPCs)> fnSigmaTPCPIDsecondParticle) { 
      // exclude all POI candidates 

      if(fReadMCTruth == kTRUE){
	label = track->GetLabel();

	//	if(label< 0)continue; //to run pn phojet
	//AliAODMCParticle *tparticle = (AliAODMCParticle*)arrayMC->At(TMath::Abs(label));
	
	tparticle = (AliAODMCParticle*)arrayMC->At(TMath::Abs(label));
	PDGcode = tparticle->GetPdgCode();
	
	//	if(TMath::Abs(tparticle->GetPdgCode())!=fPDGCodefirst)continue;
	//	if(TMath::Abs(tparticle->GetPdgCode())!=11)continue; //electron veto


	if(TMath::Abs(PDGcode)!=fPDGCodesecond){
	  
	  isMCfirst= kTRUE;
	  
	  Int_t mcMotherLabel = tparticle->GetMother();
	  Int_t mcMotherPdg = 0;
	  AliAODMCParticle *mcMother = (AliAODMCParticle*)arrayMC->At(mcMotherLabel);
	  
	  Int_t mcGrandMotherLabel = mcMother->GetMother();
	  Int_t mcGrandMotherPdg = 0; 
	  AliAODMCParticle *mcGrandMother = (AliAODMCParticle*)arrayMC->At(mcGrandMotherLabel);
	  
	  //	  if (mcMotherLabel < -1) {mcMotherPdg = 0;} else {mcMotherPdg = mcMother->GetPdgCode();} //RAMONA : era questo 02/03/16
	  if(mcMotherLabel < 0) {mcMotherPdg = 0;} else {mcMotherPdg = mcMother->GetPdgCode();} 
	  if(mcGrandMotherLabel < 0){mcGrandMotherPdg=0;}else{mcGrandMotherPdg = mcGrandMother->GetPdgCode();}
	  
	  //Mum id
	  MCmumIDP1     = mcMotherLabel;
	  MCmumPDGP1    = mcMotherPdg;
	  MCgrammaIDP1  = mcGrandMotherLabel;
	  MCgrammaPDGP1 = mcGrandMotherPdg;

	  if (tparticle->IsPhysicalPrimary()){
	    MCptcCodeP1 = 1;
	  }
	  else if (tparticle->IsSecondaryFromMaterial()){
	    MCptcCodeP1 = 2;
	  }
	  else if (tparticle->IsSecondaryFromWeakDecay()){
	    MCptcCodeP1 = 3;
	  }
	  else {
	    MCptcCodeP1 = 4;
	  }
	}
      }//end of MC loop
      
      
      if (fkCutOnTPCIP) {
	if (TMath::Abs(dz[0])>fIPCutxyPrim ) continue;  // 2.4 proton 1. pion 
	if (TMath::Abs(dz[1])>fIPCutzPrim ) continue;  // 3.2  proton 1. pion
      } else {
	   dz[0] = dzg[0];
	   dz[1] = dzg[1];
	   if (TMath::Abs(dz[0])>fIPCutxyPrim ) continue;  // 0.1 proton/pion
	   if (TMath::Abs(dz[1])>fIPCutzPrim ) continue;   // 0.15 proton/pion
	  // if (TMath::Abs(dzg[0])>fIPCutxyPrim ) continue;  // 0.1 proton/pion
	  // if (TMath::Abs(dzg[1])>fIPCutzPrim ) continue; // 0.15 proton/pion
      }
      
      
      if (track->Pt()>fMinPtForPrim  && track->Pt()< fMaxPtForPrim){
	
	if(fReadMCTruth == kTRUE)
	  fHistFirstTPCdEdx->Fill(globaltrack->Pt()*charge, globaltrack->GetTPCsignal());
	else if(fReadMCTruth == kFALSE)
	  fHistFirstTPCdEdx->Fill(globaltrack->GetTPCmomentum()*charge, globaltrack->GetTPCsignal());

	//------------------ Save first particle information

	fEvt->fReconstructedFirst[fCount].fCharge = charge;

	if(fReadMCTruth == kTRUE){
	  fEvt->fReconstructedFirst[fCount].fMomentumTruth[0]  = tparticle->Px();
	  fEvt->fReconstructedFirst[fCount].fMomentumTruth[1]  = tparticle->Py();
	  fEvt->fReconstructedFirst[fCount].fMomentumTruth[2]  = tparticle->Pz();
	}
	
	else{
	  fEvt->fReconstructedFirst[fCount].fMomentumTruth[0]  = 0.;
	  fEvt->fReconstructedFirst[fCount].fMomentumTruth[1]  = 0.;
	  fEvt->fReconstructedFirst[fCount].fMomentumTruth[2]  = 0.;
	}
	
	fEvt->fReconstructedFirst[fCount].fMomentum[0]  = track->Px();
	fEvt->fReconstructedFirst[fCount].fMomentum[1]  = track->Py();
	fEvt->fReconstructedFirst[fCount].fMomentum[2]  = track->Pz();
	fEvt->fReconstructedFirst[fCount].fPt           = track->Pt();
	fEvt->fReconstructedFirst[fCount].fEta          = track->Eta();
	fEvt->fReconstructedFirst[fCount].fPhi          = track->Phi();
	fEvt->fReconstructedFirst[fCount].fTheta        = track->Theta();
	
	fEvt->fReconstructedFirst[fCount].nSigmaFirstTPC[0]  = 0;
	fEvt->fReconstructedFirst[fCount].nSigmaFirstTPC[1]  = 0;
	
	fEvt->fReconstructedFirst[fCount].fRap    = rapidityFirst; 
	fEvt->fReconstructedFirst[fCount].mcFirstOriginType = mcFirstOrigin;
	fEvt->fReconstructedFirst[fCount].isMCptc = isMCfirst;
	fEvt->fReconstructedFirst[fCount].fMCcode = MCptcCodeP1;
	fEvt->fReconstructedFirst[fCount].fPDGcode = PDGcode;

	//cout<<"2a:------------------------------MCmumIDP1: "<<MCmumIDP1<<"   ----------------------------------MCmumPDGP1: "<<MCmumPDGP1<<endl;

	fEvt->fReconstructedFirst[fCount].fDCAxy = dz[0];
	fEvt->fReconstructedFirst[fCount].fDCAz  = dz[1];

	fEvt->fReconstructedFirst[fCount].fMCmumIdx       = MCmumIDP1;
	fEvt->fReconstructedFirst[fCount].fMCmumPDG       = MCmumPDGP1;
	fEvt->fReconstructedFirst[fCount].fMCgrandmumIdx  = MCgrammaIDP1;
	fEvt->fReconstructedFirst[fCount].fMCgrandmumPDG  = MCgrammaPDGP1;
		
	if (isP){
	  fEvt->fReconstructedFirst[fCount].isP  = kTRUE; 
	  fEvt->fReconstructedFirst[fCount].isaP = kFALSE;
	}
	else 
	  if (isaP){
	    fEvt->fReconstructedFirst[fCount].isP  = kFALSE;
	    fEvt->fReconstructedFirst[fCount].isaP = kTRUE; 
	  }

	fEvt->fReconstructedFirst[fCount].index = TMath::Abs(globaltrack->GetID());

	//      cout<<" Pos 125 cm before set "<<fEvt->fReconstructedFirst[fCount].fShiftedGlobalPosition[0]<<" "<<fEvt->fReconstructedFirst[fCount].fShiftedGlobalPosition[1]<<" "<<fEvt->fReconstructedFirst[fCount].fShiftedGlobalPosition[2]<<endl;
	if (fkPropagateGlobal) SetSftPosR125(vtrackg, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedFirst[fCount].fShiftedGlobalPosition); // FIXME: Hans popagates the TPC tracks, I try both, flag to choose 
	else SetSftPosR125(vtrack, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedFirst[fCount].fShiftedGlobalPosition);
	//cout<<" Pos set 125 cm "<<fEvt->fReconstructedFirst[fCount].fShiftedGlobalPosition[0]<<" "<<fEvt->fReconstructedFirst[fCount].fShiftedGlobalPosition[1]<<" "<<fEvt->fReconstructedFirst[fCount].fShiftedGlobalPosition[2]<<endl; 
	fEvt->fReconstructedFirst[fCount].fEtaS = EtaS(fEvt->fReconstructedFirst[fCount].fShiftedGlobalPosition); 
	//cout<<" Etas "<<fEvt->fReconstructedFirst[fCount].fEtaS<<endl;
	//cout<<fCount<<" This proton in the stack with pt "<<fEvt->fReconstructedFirst[fCount].isaP<<endl;
	
	fCount++;
	
	//---------------------- End Save first particle
	
	if(fReadMCTruth == kTRUE && isMCfirst == kTRUE)
	  fHistFirstTPCdEdxAfter->Fill(globaltrack->Pt()*charge, globaltrack->GetTPCsignal());
	else if(fReadMCTruth == kFALSE)
	  fHistFirstTPCdEdxAfter->Fill(globaltrack->GetTPCmomentum()*charge, globaltrack->GetTPCsignal());

	
	fHistIPtoPVxyzTPCFirst->Fill(dz[0],dz[1]);
	fHistIPtoPVxyzGlobalFirst->Fill(dzg[0],dzg[1]);
	
	fHistFirstTOFTPCsignalvsptAfter->Fill(track->Pt(),TMath::Sqrt(nsigmaTOFs*nsigmaTOFs+nsigmaTPCs*nsigmaTPCs));
	
	fHistyptFirst->Fill(track->Pt(),rapidityFirst,lcentrality);
	fHistphietaFirst->Fill(track->Phi(),track->Eta());
	
	fHistnTPCCrossedRFirst->Fill(nTPCCrossedRows);
	fHistRationTPCCrossedRnFindFirst->Fill(rationCrnFind);
	fHistSharedFrTPCclFirst->Fill(sharedFractionTPCcl);
      
	fHistTriggptvsCentrality->Fill(track->Pt()*charge,lcentrality);
      }
      
    } //3 sigma
 
    //cout<<"nprot : "<<fCount<<endl;
    
    if (fMaxFirstMult <= fCount){
      cerr<<"Proton counts exceeded "<<fMaxFirstMult<<"!"<<endl;
      break;
    }
  
    //Second particle
    //if(fSecondpart == kAny) fnSigmaTPCPIDsecondParticle = 0;

    if (fFirstpart != fSecondpart && TMath::Abs(nsigmaTPCs) < fnSigmaTPCPIDsecondParticle && fSecondpart != kAny) { 
      
      if(fReadMCTruth == kTRUE){
	label = track->GetLabel();
	//	cout<<"Label second particle: "<<lMCstack<<endl;
	//AliAODMCParticle *tparticle = (AliAODMCParticle*)arrayMC->At(TMath::Abs(label));
	tparticle = (AliAODMCParticle*)arrayMC->At(TMath::Abs(label));
	PDGcode = tparticle->GetPdgCode();
	//	if(TMath::Abs(tparticle->GetPdgCode())!=fPDGCodesecond)continue;
	//	if(TMath::Abs(tparticle->GetPdgCode())!=11)continue; //electron veto

	

	if(TMath::Abs(PDGcode)==fPDGCodesecond){
	  isMCsecond= kTRUE;
	  /*
	   
	    
	    Int_t mcMotherLabel = tparticle->GetMother();
	    Int_t mcMotherPdg = 0;
	    AliAODMCParticle *mcMother = (AliAODMCParticle*)arrayMC->At(mcMotherLabel);
	  
	    Int_t mcGrandMotherLabel = mcMother->GetMother();
	    Int_t mcGrandMotherPdg = 0; 
	    AliAODMCParticle *mcGrandMother = (AliAODMCParticle*)arrayMC->At(mcGrandMotherLabel);
	  
	    //	  if (mcMotherLabel < -1) {mcMotherPdg = 0;} else {mcMotherPdg = mcMother->GetPdgCode();} //RAMONA : era questo 02/03/16
	    if(mcMotherLabel < 0) {mcMotherPdg = 0;} else {mcMotherPdg = mcMother->GetPdgCode();} 
	    if(mcGrandMotherLabel < 0){mcGrandMotherPdg=0;}else{mcGrandMother = mcGrandMother->GetPdgCode();}
	    // cout<<"mcMotherlabel: "<<mcMotherLabel<<endl;
	    //	  cout<<"------------------------------MCmumIDP1: "<<MCmumIDP1<<"   ----------------------------------MCmumPDGP1: "<<MCmumPDGP1<<endl;

	    //Mum id
	    MCmumIDP1     = mcMotherLabel;
	    MCmumPDGP1    = mcMotherPdg;
	    MCgrammaIDP1  = mcGrandMotherLabel 
	    MCgrammaPDGP1 = mcGrandMotherPdg;

	  */
	  
	  Int_t mcMotherLabel = tparticle->GetMother();
	  Int_t mcMotherPdg = 0;
	  AliAODMCParticle *mcMother = (AliAODMCParticle*)arrayMC->At(mcMotherLabel);
	  
	  Int_t mcGrandMotherLabel = mcMother->GetMother();
	  Int_t mcGrandMotherPdg = 0; 
	  AliAODMCParticle *mcGrandMother = (AliAODMCParticle*)arrayMC->At(mcGrandMotherLabel);
	  
	  //	  if (mcMotherLabel < -1) {mcMotherPdg = 0;} else {mcMotherPdg = mcMother->GetPdgCode();} //RAMONA : era questo 02/03/16
	  if(mcMotherLabel < 0) {mcMotherPdg = 0;} else {mcMotherPdg = mcMother->GetPdgCode();} 
	  if(mcGrandMotherLabel < 0){mcGrandMotherPdg=0;}else{mcGrandMotherPdg = mcGrandMother->GetPdgCode();}
	  // cout<<"mcMotherlabel: "<<mcMotherLabel<<endl;

	  //Mum id
	  MCmumIDP2     = mcMotherLabel;
	  MCmumPDGP2    = mcMotherPdg;
	  MCgrammaIDP2  = mcGrandMotherLabel;
	  MCgrammaPDGP2 = mcGrandMotherPdg;

	  //cout<<"1:------------------------------MCmumIDP2: "<<MCmumIDP2<<"   ----------------------------------MCmumPDGP2: "<<MCmumPDGP2<<endl;
		  
	  if (tparticle->IsPhysicalPrimary()){
	    MCptcCodeP2 = 1;
	  }
	  else if (tparticle->IsSecondaryFromMaterial()){
	    MCptcCodeP2 = 2;
	  }
	  else if (tparticle->IsSecondaryFromWeakDecay()){
	    MCptcCodeP2 = 3;
	  }
	  else {
	    MCptcCodeP2 = 4;
	  }
	}
      }
       

      // if( fSecondpart == 4){ //electron veto in p
      // 	if(track->Pt()>0.600 && track->Pt()<0.800)continue;
      // }

      if (fkCutOnTPCIP) {
	if (TMath::Abs(dz[0])>fIPCutxySec ) continue;  // 2.4 proton 1. pion 
	if (TMath::Abs(dz[1])>fIPCutzSec ) continue;  // 3.2  proton 1. pion
      } else {
	if (TMath::Abs(dzg[0])>fIPCutxySec ) continue;  // 0.1 proton/pion
	if (TMath::Abs(dzg[1])>fIPCutzSec ) continue; // 0.15 proton/pion
      }
      
      if (track->Pt()>=fMinPtForSec && track->Pt()< fMomemtumLimitForTOFPIDsecond){	
	
	if(fReadMCTruth == kTRUE)
	  fHistSecondTPCdEdx->Fill(globaltrack->Pt()*charge, globaltrack->GetTPCsignal());
	else if(fReadMCTruth == kFALSE)
	  fHistSecondTPCdEdx->Fill(globaltrack->GetTPCmomentum()*charge, globaltrack->GetTPCsignal());

	
	//	if(nsigmaTOFs > 5)continue;
	//------------------------------ Save second particle information
	fEvt->fReconstructedSecond[sCount].sCharge = charge;
	//cout<<"Charge of pion "<<fEvt->fReconstructedSecond[sCount].fCharge<<endl;

	if(fReadMCTruth == kTRUE){
	  fEvt->fReconstructedSecond[sCount].sMomentumTruth[0]  = tparticle->Px();
	  fEvt->fReconstructedSecond[sCount].sMomentumTruth[1]  = tparticle->Py();
	  fEvt->fReconstructedSecond[sCount].sMomentumTruth[2]  = tparticle->Pz();
	  //   fEvt->fReconstructedSecond[sCount].fPt           = tparticle->Pt();
	  //   fEvt->fReconstructedSecond[sCount].fEta          = tparticle->Eta();
	  //   fEvt->fReconstructedSecond[sCount].fPhi          = tparticle->Phi();
	}
	
	else{
	  fEvt->fReconstructedSecond[sCount].sMomentumTruth[0]  = 0.;
	  fEvt->fReconstructedSecond[sCount].sMomentumTruth[1]  = 0.;
	  fEvt->fReconstructedSecond[sCount].sMomentumTruth[2]  = 0.;
	}
	
	fEvt->fReconstructedSecond[sCount].sMomentum[0]  = track->Px();
	fEvt->fReconstructedSecond[sCount].sMomentum[1]  = track->Py();
	fEvt->fReconstructedSecond[sCount].sMomentum[2]  = track->Pz();
	fEvt->fReconstructedSecond[sCount].sPt           = track->Pt();
	fEvt->fReconstructedSecond[sCount].sEta          = track->Eta();
	fEvt->fReconstructedSecond[sCount].sPhi          = track->Phi();
	fEvt->fReconstructedSecond[sCount].sTheta        = track->Theta();
	
	// for (int pIndex = 0; pIndex <5; pIndex++) {
	//   fEvt->fReconstructedSecond[sCount].nSigmaSecondTPC[pIndex]  = pullsTPC[pIndex];
	//   fEvt->fReconstructedSecond[sCount].nSigmaSecondTOF[pIndex]  = pullsTOF[pIndex];
	// }	
	
	// for (int pIndex = 0; pIndex <3; pIndex++) {
	//   fEvt->fReconstructedSecond[sCount].sMomentumTruth[pIndex]  = pMomentumTruth[pIndex];
	// }
	
	fEvt->fReconstructedSecond[sCount].sRap    = rapiditySecond; 
	fEvt->fReconstructedSecond[sCount].mcSecondOriginType = mcSecondOrigin;
	fEvt->fReconstructedSecond[sCount].isMCptc = isMCsecond;
	fEvt->fReconstructedSecond[sCount].sMCcode = MCptcCodeP2;
	fEvt->fReconstructedSecond[sCount].sPDGcode = PDGcode;
			    	   	   
	fEvt->fReconstructedSecond[sCount].sDCAxy = dz[0];
	fEvt->fReconstructedSecond[sCount].sDCAz  = dz[1];
			    	   	   
	fEvt->fReconstructedSecond[sCount].sMCmumIdx       = MCmumIDP2;
	fEvt->fReconstructedSecond[sCount].sMCmumPDG       = MCmumPDGP2;
	fEvt->fReconstructedSecond[sCount].sMCgrandmumIdx  = MCgrammaIDP2;
	fEvt->fReconstructedSecond[sCount].sMCgrandmumPDG  = MCgrammaPDGP2;

	//	cout<<"2a:------------------------------MCmumIDP2: "<<MCmumIDP2<<"   ----------------------------------MCmumPDGP2: "<<MCmumPDGP2<<endl;

	if (isP){
	  fEvt->fReconstructedSecond[sCount].isP  = kTRUE; 
	  fEvt->fReconstructedSecond[sCount].isaP = kFALSE;
	}
	else if (isaP){
	  fEvt->fReconstructedSecond[sCount].isP  = kFALSE;
	  fEvt->fReconstructedSecond[sCount].isaP = kTRUE; 
	}
	
	
	fEvt->fReconstructedSecond[sCount].index = TMath::Abs(globaltrack->GetID());
	//      cout<<" Pos 125 cm before set "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[0]<<" "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[1]<<" "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[2]<<endl;
	if (fkPropagateGlobal) SetSftPosR125(vtrackg, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition); // FIXME: Hans popagates the TPC tracks, I try both, flag to choose 
	else SetSftPosR125(vtrack, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition);
	//      cout<<" Pos set 125 cm "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[0]<<" "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[1]<<" "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[2]<<endl; 
	fEvt->fReconstructedSecond[sCount].sEtaS = EtaS(fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition); 
	//      cout<<" Etas "<<fEvt->fReconstructedSecond[sCount].sEtaS<<endl;
	
	//    cout<<" This proton in the stack with pt "<<fEvt->fReconstructedProton[sCount].sPt<<endl;
	
	sCount++;
	
	//------------------------------ Save second particle information

	if(fReadMCTruth == kTRUE && isMCsecond == kTRUE)
	  fHistSecondTPCdEdxAfter->Fill(globaltrack->Pt()*charge, globaltrack->GetTPCsignal());
	else if(fReadMCTruth == kFALSE)
	  fHistSecondTPCdEdxAfter->Fill(globaltrack->GetTPCmomentum()*charge, globaltrack->GetTPCsignal());

	fHistIPtoPVxyzTPCSecond->Fill(dz[0],dz[1]);
	fHistIPtoPVxyzGlobalSecond->Fill(dzg[0],dzg[1]);
	
	fHistSecondTOFTPCsignalvsptAfter->Fill(track->Pt(),TMath::Sqrt(nsigmaTOFs*nsigmaTOFs+nsigmaTPCs*nsigmaTPCs));
	
	fHistyptSecond->Fill(track->Pt(),rapiditySecond,lcentrality);
	fHistphietaSecond->Fill(track->Phi(),track->Eta());
	
	fHistnTPCCrossedRSecond->Fill(nTPCCrossedRows);
	fHistRationTPCCrossedRnFindSecond->Fill(rationCrnFind);
	fHistSharedFrTPCclSecond->Fill(sharedFractionTPCcl);
      }
      
      else if (track->Pt() >= fMomemtumLimitForTOFPIDsecond && track->Pt() < fMaxPtForSec && isTOFPIDok) {
	
	//	tTOF -= expectedTimes[(AliPID::EParticleType)fSecondpart];  
	
	length = globaltrack->GetIntegratedLength();  // // FIXME length is zero!! from a mail february 2014: this info is not available for AODs 115, use AODs 145
	
	ptot = globaltrack->GetTPCmomentum();
		
	fHistSecondTOFmisvspt->Fill(track->Pt(),probMis);
	fHistSecondTOFmisvsp->Fill(track->P(),probMis);
	fHistSecondTOFnsigmavspt->Fill(track->Pt(),nsigmaTOFs);
	fHistSecondTOFnsigmavsp->Fill(track->P(),nsigmaTPCs);
	fHistSecondTOFsignalvsp->Fill(track->P(),tTOF- expectedTimes[(AliPID::EParticleType)fSecondpart]);  
	fHistSecondTOFsignalvspt->Fill(track->Pt(),tTOF- expectedTimes[(AliPID::EParticleType)fSecondpart]);
	  
	fHistSecondTOFTPCsignalvspt->Fill(track->Pt(),TMath::Sqrt(nsigmaTOFs*nsigmaTOFs+nsigmaTPCs*nsigmaTPCs));
	
	//if (track->P()< fMomemtumLimitForTOFPIDsecond)continue; // FIXME : to here
	//if (probMis > 0.01) continue; 
	//if (TMath::Sqrt(nsigmaTOFs*nsigmaTOFs+nsigmaTPCs*nsigmaTPCs)> fnSigmaTPCTOFPIDsecondParticle) continue;   
	// this cleans the TOF corrected time plot vs p       //FIXME : Original line from mariella
	
	//	if (TMath::Abs(nsigmaTOFs)> fnSigmaTPCTOFPIDsecondParticle) continue;   // this cleans the TOF corrected time plot vs p     
	
	beta = length/(tTOF*2.99792457999999984e-02); 
	//	cout<<" rack length  "<<length<<" beta "<<beta<<endl; 
	gamma = 1/TMath::Sqrt(1 - beta*beta);
	mass = ptot/TMath::Sqrt(gamma*gamma - 1); // using inner TPC mom. as approx.
	//	  cout<<" TMath::Sqrt(1 - beta*beta) "<< TMath::Sqrt(1 - beta*beta)<<" mass "<<mass<<endl;

	if(fSecondpart == 5){
	  if(TMath::Abs(mass) > 2.65)continue;
	  if(TMath::Abs(mass) < 1.05)continue;
	}

	fHistSecondMassTOFvsPt3sTPC->Fill(mass-fPDGsecond,track->Pt()*charge);
	
	if (TMath::Abs(nsigmaTOFs)< fnSigmaTPCTOFPIDsecondParticle && tTOF > 0. /*&&length>350*/) {
	  
	  if(fReadMCTruth == kTRUE)
	    fHistSecondTPCdEdx->Fill(globaltrack->Pt()*charge, globaltrack->GetTPCsignal());
	  else if(fReadMCTruth == kFALSE)
	    fHistSecondTPCdEdx->Fill(globaltrack->GetTPCmomentum()*charge, globaltrack->GetTPCsignal());

	  fHistSecondMassTOFvsPt3sTPC3sTOF->Fill(mass-fPDGsecond,track->Pt()*charge);
	  
	  //------------------------------ Save second particle information
	  fEvt->fReconstructedSecond[sCount].sCharge = charge;
	  //cout<<"Charge of pion "<<fEvt->fReconstructedSecond[sCount].fCharge<<endl;
	 	
	  if(fReadMCTruth == kTRUE){
	    fEvt->fReconstructedSecond[sCount].sMomentumTruth[0]  = tparticle->Px();
	    fEvt->fReconstructedSecond[sCount].sMomentumTruth[1]  = tparticle->Py();
	    fEvt->fReconstructedSecond[sCount].sMomentumTruth[2]  = tparticle->Pz();
	    //   fEvt->fReconstructedSecond[sCount].fPt           = tparticle->Pt();
	    //   fEvt->fReconstructedSecond[sCount].fEta          = tparticle->Eta();
	    //   fEvt->fReconstructedSecond[sCount].fPhi          = tparticle->Phi();
	  }
	
	  else{
	    fEvt->fReconstructedSecond[sCount].sMomentumTruth[0]  = 0.;
	    fEvt->fReconstructedSecond[sCount].sMomentumTruth[1]  = 0.;
	    fEvt->fReconstructedSecond[sCount].sMomentumTruth[2]  = 0.;
	  }
	

	  fEvt->fReconstructedSecond[sCount].sMomentum[0]  = track->Px();
	  fEvt->fReconstructedSecond[sCount].sMomentum[1]  = track->Py();
	  fEvt->fReconstructedSecond[sCount].sMomentum[2]  = track->Pz();
	  fEvt->fReconstructedSecond[sCount].sPt           = track->Pt();
	  fEvt->fReconstructedSecond[sCount].sEta          = track->Eta();
	  fEvt->fReconstructedSecond[sCount].sPhi          = track->Phi();
	  fEvt->fReconstructedSecond[sCount].sTheta        = track->Theta();
	  fEvt->fReconstructedSecond[sCount].sTOFMass      = mass;

	  //	  }
	 
	  // for (int pIndex = 0; pIndex <5; pIndex++) {
	  //   fEvt->fReconstructedSecond[sCount].nSigmaSecondTPC[pIndex]  = pullsTPC[pIndex];
	  //   fEvt->fReconstructedSecond[sCount].nSigmaSecondTOF[pIndex]  = pullsTOF[pIndex];
	  // } 

	  //   for (int pIndex = 0; pIndex <3; pIndex++) {
	  //   fEvt->fReconstructedSecond[sCount].sMomentumTruth[pIndex]  = pMomentumTruth[pIndex];
	  // }
	
	  fEvt->fReconstructedSecond[sCount].sRap    = rapiditySecond; 
	  fEvt->fReconstructedSecond[sCount].mcSecondOriginType = mcSecondOrigin;
	  fEvt->fReconstructedSecond[sCount].isMCptc = isMCsecond;
	  fEvt->fReconstructedSecond[sCount].sMCcode = MCptcCodeP2;
	  fEvt->fReconstructedSecond[sCount].sPDGcode = PDGcode;
	
	  fEvt->fReconstructedSecond[sCount].sDCAxy = dz[0];
	  fEvt->fReconstructedSecond[sCount].sDCAz  = dz[1];
	  			     
	  fEvt->fReconstructedSecond[sCount].sMCmumIdx = MCmumIDP2;
	  fEvt->fReconstructedSecond[sCount].sMCmumPDG  = MCmumPDGP2;
	  fEvt->fReconstructedSecond[sCount].sMCgrandmumIdx  = MCgrammaIDP2;
	  fEvt->fReconstructedSecond[sCount].sMCgrandmumPDG  = MCgrammaPDGP2;

	  //	  cout<<"2b:------------------------------MCmumIDP2: "<<MCmumIDP2<<"   ----------------------------------MCmumPDGP2: "<<MCmumPDGP2<<endl;

	  if (isP){
	    fEvt->fReconstructedSecond[sCount].isP  = kTRUE; 
	    fEvt->fReconstructedSecond[sCount].isaP = kFALSE;
	  }
	  else if (isaP){
	    fEvt->fReconstructedSecond[sCount].isP  = kFALSE;
	    fEvt->fReconstructedSecond[sCount].isaP = kTRUE; 
	  }
	
	  fEvt->fReconstructedSecond[sCount].index = TMath::Abs(globaltrack->GetID());
	  //      cout<<" Pos 125 cm before set "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[0]<<" "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[1]<<" "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[2]<<endl;
	  if (fkPropagateGlobal) SetSftPosR125(vtrackg, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition); // FIXME: Hans popagates the TPC tracks, I try both, flag to choose 
	  else SetSftPosR125(vtrack, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition);
	  //      cout<<" Pos set 125 cm "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[0]<<" "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[1]<<" "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[2]<<endl; 
	  fEvt->fReconstructedSecond[sCount].sEtaS = EtaS(fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition); 
	  //      cout<<" Etas "<<fEvt->fReconstructedSecond[sCount].sEtaS<<endl;
	  
	  //    cout<<" This proton in the stack with pt "<<fEvt->fReconstructedProton[sCount].sPt<<endl;
	  
	  sCount++;
	  
	  //------------------------------ Save second particle information
	  //}
	  // }
	  
	  // else { 
	  //   continue; 
	  // }
	  
	  
	  //if (track->TestBit(AliAODTrack::kIsDCA)) {fHistIPtoPVxyzTPCSecond->Fill(dz[0],dz[1]); } else { /*cout<<"Bit is not DCA!!!!!!!!!!!"<<endl;*/ }
	  fHistIPtoPVxyzTPCSecond->Fill(dz[0],dz[1]);
	  fHistIPtoPVxyzGlobalSecond->Fill(dzg[0],dzg[1]);

	  if(fReadMCTruth == kTRUE && isMCsecond == kTRUE)
	    fHistSecondTPCdEdxAfter->Fill(globaltrack->Pt()*charge, globaltrack->GetTPCsignal());
	  else if(fReadMCTruth == kFALSE)
	    fHistSecondTPCdEdxAfter->Fill(globaltrack->GetTPCmomentum()*charge, globaltrack->GetTPCsignal());
	  
	  fHistSecondTOFTPCsignalvsptAfter->Fill(track->Pt(),TMath::Sqrt(nsigmaTOFs*nsigmaTOFs+nsigmaTPCs*nsigmaTPCs));
	  
	  fHistyptSecond->Fill(track->Pt(),rapiditySecond,lcentrality);
	  fHistphietaSecond->Fill(track->Phi(),track->Eta());
      
	  fHistnTPCCrossedRSecond->Fill(nTPCCrossedRows);
	  fHistRationTPCCrossedRnFindSecond->Fill(rationCrnFind);
	  fHistSharedFrTPCclSecond->Fill(sharedFractionTPCcl);
	}//nsigma tof
      }//tof
    }//3sigma 2 ptc
   
    else if(fSecondpart == kAny){
      
      // if( TMath::Abs(globaltrack->GetID() )

      // fEvt->fReconstructedSecond[sCount].index = TMath::Abs(globaltrack->GetID());
  
      
      //ramonaANY
      fHistSecondMassTOFvsPt3sTPC->Fill(mass-fPDGsecond,track->Pt()*charge);
      
	  
      if(fReadMCTruth == kTRUE)
	fHistSecondTPCdEdx->Fill(globaltrack->Pt()*charge, globaltrack->GetTPCsignal());
      else if(fReadMCTruth == kFALSE)
	fHistSecondTPCdEdx->Fill(globaltrack->GetTPCmomentum()*charge, globaltrack->GetTPCsignal());
      
      fHistSecondMassTOFvsPt3sTPC3sTOF->Fill(mass-fPDGsecond,track->Pt()*charge);
      
      //------------------------------ Save second particle information
      fEvt->fReconstructedSecond[sCount].sCharge = charge;
      //cout<<"Charge of pion "<<fEvt->fReconstructedSecond[sCount].fCharge<<endl;
      
      if(fReadMCTruth == kTRUE){
	fEvt->fReconstructedSecond[sCount].sMomentumTruth[0]  = tparticle->Px();
	fEvt->fReconstructedSecond[sCount].sMomentumTruth[1]  = tparticle->Py();
	fEvt->fReconstructedSecond[sCount].sMomentumTruth[2]  = tparticle->Pz();
	//   fEvt->fReconstructedSecond[sCount].fPt           = tparticle->Pt();
	//   fEvt->fReconstructedSecond[sCount].fEta          = tparticle->Eta();
	//   fEvt->fReconstructedSecond[sCount].fPhi          = tparticle->Phi();
      }
      
      else{
	fEvt->fReconstructedSecond[sCount].sMomentumTruth[0]  = 0.;
	fEvt->fReconstructedSecond[sCount].sMomentumTruth[1]  = 0.;
	fEvt->fReconstructedSecond[sCount].sMomentumTruth[2]  = 0.;
      }
      
      
      fEvt->fReconstructedSecond[sCount].sMomentum[0]  = track->Px();
      fEvt->fReconstructedSecond[sCount].sMomentum[1]  = track->Py();
      fEvt->fReconstructedSecond[sCount].sMomentum[2]  = track->Pz();
      fEvt->fReconstructedSecond[sCount].sPt           = track->Pt();
      fEvt->fReconstructedSecond[sCount].sEta          = track->Eta();
      fEvt->fReconstructedSecond[sCount].sPhi          = track->Phi();
      fEvt->fReconstructedSecond[sCount].sTheta        = track->Theta();
      fEvt->fReconstructedSecond[sCount].sTOFMass      = mass;
      
      //	  }
      
      // for (int pIndex = 0; pIndex <5; pIndex++) {
      //   fEvt->fReconstructedSecond[sCount].nSigmaSecondTPC[pIndex]  = pullsTPC[pIndex];
      //   fEvt->fReconstructedSecond[sCount].nSigmaSecondTOF[pIndex]  = pullsTOF[pIndex];
      // } 
      
      //   for (int pIndex = 0; pIndex <3; pIndex++) {
      //   fEvt->fReconstructedSecond[sCount].sMomentumTruth[pIndex]  = pMomentumTruth[pIndex];
      // }
      
      fEvt->fReconstructedSecond[sCount].sRap    = rapiditySecond; 
      fEvt->fReconstructedSecond[sCount].mcSecondOriginType = mcSecondOrigin;
      fEvt->fReconstructedSecond[sCount].isMCptc = isMCsecond;
      fEvt->fReconstructedSecond[sCount].sMCcode = MCptcCodeP2;
      fEvt->fReconstructedSecond[sCount].sPDGcode = PDGcode;
      
      fEvt->fReconstructedSecond[sCount].sDCAxy = dz[0];
      fEvt->fReconstructedSecond[sCount].sDCAz  = dz[1];
      
      fEvt->fReconstructedSecond[sCount].sMCmumIdx = MCmumIDP2;
      fEvt->fReconstructedSecond[sCount].sMCmumPDG  = MCmumPDGP2;
      fEvt->fReconstructedSecond[sCount].sMCgrandmumIdx  = MCgrammaIDP2;
      fEvt->fReconstructedSecond[sCount].sMCgrandmumPDG  = MCgrammaPDGP2;
      
      //	  cout<<"2b:------------------------------MCmumIDP2: "<<MCmumIDP2<<"   ----------------------------------MCmumPDGP2: "<<MCmumPDGP2<<endl;
      
      if (isP){
	fEvt->fReconstructedSecond[sCount].isP  = kTRUE; 
	fEvt->fReconstructedSecond[sCount].isaP = kFALSE;
      }
      else if (isaP){
	fEvt->fReconstructedSecond[sCount].isP  = kFALSE;
	fEvt->fReconstructedSecond[sCount].isaP = kTRUE; 
      }
      
      fEvt->fReconstructedSecond[sCount].index = TMath::Abs(globaltrack->GetID());
      
      //      cout<<" Pos 125 cm before set "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[0]<<" "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[1]<<" "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[2]<<endl;
      if (fkPropagateGlobal) SetSftPosR125(vtrackg, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition); // FIXME: Hans popagates the TPC tracks, I try both, flag to choose 
      else SetSftPosR125(vtrack, bfield, lBestPrimaryVtxPos, fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition);
      //      cout<<" Pos set 125 cm "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[0]<<" "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[1]<<" "<<fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition[2]<<endl; 
      fEvt->fReconstructedSecond[sCount].sEtaS = EtaS(fEvt->fReconstructedSecond[sCount].sShiftedGlobalPosition); 
      //      cout<<" Etas "<<fEvt->fReconstructedSecond[sCount].sEtaS<<endl;
      
      //    cout<<" This proton in the stack with pt "<<fEvt->fReconstructedProton[sCount].sPt<<endl;
      
      sCount++;
      
      //------------------------------ Save second particle information
      //}
      // }
      
      // else { 
      //   continue; 
      // }
      
      
      //if (track->TestBit(AliAODTrack::kIsDCA)) {fHistIPtoPVxyzTPCSecond->Fill(dz[0],dz[1]); } else { /*cout<<"Bit is not DCA!!!!!!!!!!!"<<endl;*/ }
      fHistIPtoPVxyzTPCSecond->Fill(dz[0],dz[1]);
      fHistIPtoPVxyzGlobalSecond->Fill(dzg[0],dzg[1]);
      
      if(fReadMCTruth == kTRUE && isMCsecond == kTRUE)
	fHistSecondTPCdEdxAfter->Fill(globaltrack->Pt()*charge, globaltrack->GetTPCsignal());
      else if(fReadMCTruth == kFALSE)
	fHistSecondTPCdEdxAfter->Fill(globaltrack->GetTPCmomentum()*charge, globaltrack->GetTPCsignal());
      
      fHistSecondTOFTPCsignalvsptAfter->Fill(track->Pt(),TMath::Sqrt(nsigmaTOFs*nsigmaTOFs+nsigmaTPCs*nsigmaTPCs));
      
      fHistyptSecond->Fill(track->Pt(),rapiditySecond,lcentrality);
      fHistphietaSecond->Fill(track->Phi(),track->Eta());
      
      fHistnTPCCrossedRSecond->Fill(nTPCCrossedRows);
      fHistRationTPCCrossedRnFindSecond->Fill(rationCrnFind);
      fHistSharedFrTPCclSecond->Fill(sharedFractionTPCcl);
    }
   
  
 
    if (fMaxSecondMult <= sCount){
      cerr<<"Proton counts exceeded "<<fMaxSecondMult<<"!"<<endl;
      break;
    }
    
  
    // //RAMONA
    // fEvt->fNumberCandidateFirst = fCount;
    // fHistFirstMultvsCent->Fill(fCount,lcentrality);
    
    // fEvt->fNumberCandidateSecond = sCount;
    // fHistFirstMultvsCent->Fill(sCount,lcentrality);
  
  } //14/02/17 --> Parentesi inutile?? --> pero e Loop tracce
  
  fEvt->fNumberCandidateFirst = fCount;
  fHistFirstMultvsCent->Fill(fCount,lcentrality);
  
  fEvt->fNumberCandidateSecond = sCount;
  fHistSecondMultvsCent->Fill(sCount,lcentrality);
  
  // cout<<"Number first: "<<fEvt->fNumberCandidateFirst <<endl;
  // cout<<"Number second: "<<fEvt->fNumberCandidateSecond<<endl;
  if(fFirstpart == fSecondpart || fSecondpart == kAny)
    // DoPairshh(lcentrality, fieldsign);  
    //    DoPairshh(lcentrality, fieldsign,fSphericityvalue);  
    DoPairshh(centralityBin, fieldsign,fSphericityvalue);  
  else{
    //Remove candidates that are at the same time a ptc1 and ptc2
    for (int i=0; i < fEvt->fNumberCandidateFirst; i++) {
      for (int j=0; j<fEvt->fNumberCandidateSecond; j++) {
	if (fEvt->fReconstructedFirst[i].index == fEvt->fReconstructedSecond[j].index) {
	  // cout<<"the track can be both tracks!"<<endl;
	  // cout<<"i: "<<i<<" j: "<<j<<" idxi: "<<fEvt->fReconstructedFirst[i].index<<" idxj: "<<fEvt->fReconstructedSecond[j].index<<endl;
	  fEvt->fReconstructedFirst[i].doSkipOver = kTRUE;
	  fEvt->fReconstructedSecond[j].doSkipOver = kTRUE;
	}
      }
    }
    //--------------------------------------------------------------
    DoPairsh1h2(lcentrality, fieldsign,fSphericityvalue);  
    //DoPairsh1h2(centralityBin, fieldsign,fSphericityvalue);  
  }
  // Post output data
  PostData(1, fOutputContainer);
  PostData(2, fHistSparseSignal );
  PostData(3, fHistSparseBkg );
        
}

//----------------------------------------------------------------------------------------------------
  
//void AliAnalysisTaskDeuFlow2PC::DoPairsh1h2 ( const Float_t lcentrality, int fieldsignf, double fSphericityvalue) {
void AliAnalysisTaskDeuFlow2PC::DoPairsh1h2 ( const Float_t lcentrality, int fieldsign, const Double_t fSphericityvalue) {

  //  cout<<"Ramona Enter DoPair"<<endl;

  // UInt_t dimsparse;
  // if(!fReadMCTruth) dimsparse=16;
  // else dimsparse=25;
  // Double_t xsparseSignal[dimsparse];
  // for(UInt_t counter = 0; counter < dimsparse ; counter++) xsparseSignal[counter]=-999;//Default position for THnSparse 
  // Double_t xsparseBkg[dimsparse];
  // for(UInt_t counter = 0; counter < dimsparse ; counter++) xsparseBkg[counter]=-999;//Default position for THnSparse 

  //-----------
  double DCAxyP1  = -999. ;
  double DCAzP1   = -999. ;  
  double DCAxyP2  = -999. ; 
  double DCAzP2   = -999. ;  

  double  ptP1 = -999.;
  double  ptP2 = -999.;

  // Short_t chargeP1 = -999.;
  // Short_t chargeP2 = -999.;

  bool isP1  = kFALSE;
  bool isaP1 = kFALSE;
  bool isP2  = kFALSE;
  bool isaP2 = kFALSE;

  Int_t  SignP1 = -999;
  Int_t  SignP2 = -999;

  double dphi  = -999.;
  double dphis = -999.;
  double deta  = -999.;
  //  double detas = -999.;
  double dtheta = -999.;
  
  bool isMC1 = kFALSE;
  bool isMC2 = kFALSE;

  bool isMCvector = kFALSE;
  bool sameMother = kFALSE;
  bool sameGrandMother = kFALSE;
  
  Int_t mcMotherLabelP1 = -999;
  Int_t mcMotherLabelP2 = -999;

  // Int_t mcGrandMotherLabelP1 = -999;
  // Int_t mcGrandMotherLabelP2 = -999;
 
  Int_t typeP1 = -999;
  Int_t typeP2 = -999;
 
  Int_t mcPDGMotherP1 = 0;
  Int_t mcPDGMotherP2 = 0;
  //  Int_t mcMotherBin = 0;
 
  //  Int_t mcPDGGrandMother = 0;
  Int_t mcGrandMotherBin = 0;

  Int_t mcPDGcodeP1 = 0;
  Int_t mcPDGcodeP2 = 0;
  // Int_t mcPDG1Bin = 0;
  // Int_t mcPDG2Bin = 0;

  int evmultmixed = 0;
  bool multmixedcounted = kFALSE;
 
  double pairKstar   = 0.;
  double pairKstarMC = 0.;
  double pairMass  = 0.;
  //  double pairMassE = 0.;
  double pairKt    = 0.;

  double massS = 0;

  for (int i=0; i<fEvt->fNumberCandidateFirst; i++) {

    if (fEvt->fReconstructedFirst[i].doSkipOver) continue;
    
    DCAxyP1         = fEvt->fReconstructedFirst[i].fDCAxy;
    DCAzP1          = fEvt->fReconstructedFirst[i].fDCAz;
    ptP1            = fEvt->fReconstructedFirst[i].fPt;
    
    if(ptP1 < fMinPtForPrim)continue;
    
    //    chargeP1        = fEvt->fReconstructedFirst[i].fCharge;
    isMC1           = fEvt->fReconstructedFirst[i].isMCptc;
    mcMotherLabelP1 = fEvt->fReconstructedFirst[i].fMCmumIdx;
    typeP1          = fEvt->fReconstructedFirst[i].fMCcode;
    //  if(typeP1 < 0) continue;
    //    mcGrandMotherLabelP1 = fEvt->fReconstructedFirst[i].fMCgrandmumIdx;
    mcPDGcodeP1          = fEvt->fReconstructedFirst[i].fPDGcode;

    
    isP1  = fEvt->fReconstructedFirst[i].isP;
    isaP1 = fEvt->fReconstructedFirst[i].isaP;

    if(isP1) SignP1 = 1;
    else if (isaP1) SignP1 = -1;

    mcPDGMotherP1    = fEvt->fReconstructedFirst[i].fMCmumPDG;
    //    mcPDGGrandMother = fEvt->fReconstructedFirst[i].fMCgrandmumPDG;
    
    for (int eventNumber=0; eventNumber<fnEventsToMix+1; eventNumber++) { 
      // For same event pairs
      // if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateFirst)!=0.) evmultmixed++; 
      if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateSecond)!=0.) evmultmixed++; 

      
      for (int j=0; j<(fEvt+eventNumber)->fNumberCandidateSecond; j++) {

	// cout<<"Ramona event number "<<eventNumber<<endl;
	// cout<<"eventNumber: "<<eventNumber<<" j: "<<j<<" i: "<<i<<" fEvt+eventNumber: "<<(fEvt+eventNumber)<<endl;
	
	if ((fEvt+eventNumber)->fReconstructedSecond[j].doSkipOver) continue;
        
	//	cout<<"Ramona Doing stuff Second Ptc"<<endl;
	
	DCAxyP2         = (fEvt+eventNumber)->fReconstructedSecond[j].sDCAxy;
	DCAzP2          = (fEvt+eventNumber)->fReconstructedSecond[j].sDCAz;
	ptP2            = (fEvt+eventNumber)->fReconstructedSecond[j].sPt;
	if(ptP2< fMinPtForSec)continue;
	//chargeP2        = (fEvt+eventNumber)->fReconstructedSecond[j].sCharge;
	isMC2           = (fEvt+eventNumber)->fReconstructedSecond[j].isMCptc;
	mcMotherLabelP2 = (fEvt+eventNumber)->fReconstructedSecond[j].sMCmumIdx;
	typeP2          = (fEvt+eventNumber)->fReconstructedSecond[j].sMCcode;
	//	if(typeP2 < 0) continue;
	//	mcGrandMotherLabelP2 = (fEvt+eventNumber)->fReconstructedSecond[j].sMCgrandmumIdx;
  	mcPDGcodeP2     = (fEvt+eventNumber)->fReconstructedSecond[j].sPDGcode;
	massS           = (fEvt+eventNumber)->fReconstructedSecond[j].sTOFMass;

	isP2  = (fEvt+eventNumber)->fReconstructedSecond[j].isP;
        isaP2 = (fEvt+eventNumber)->fReconstructedSecond[j].isaP;
	
	mcPDGMotherP2 = (fEvt+eventNumber)->fReconstructedSecond[j].sMCmumPDG;
	
	if(isP2) SignP2 = 1;
	else if (isaP2) SignP2 = -1;

	if(isMC1 && isMC2) isMCvector = kTRUE;
	else isMCvector = kFALSE;
	
	if(mcMotherLabelP1 == mcMotherLabelP2 && mcMotherLabelP1!=-999)sameMother = kTRUE;
	else sameMother = kFALSE;
	
	//Calculate k* for the pair
        pairKstar = CalculateKstar(fEvt->fReconstructedFirst[i].fMomentum, (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum,fPDGfirst, fPDGsecond);
	
	//mc kstar
	pairKstarMC = CalculateKstar(fEvt->fReconstructedFirst[i].fMomentumTruth, (fEvt+eventNumber)->fReconstructedSecond[j].sMomentumTruth,fPDGfirst, fPDGsecond);
	
	//Invariant Mass of the pair
	pairMass  = CalculateMass(fEvt->fReconstructedFirst[i].fMomentum, (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum,fPDGfirst, fPDGsecond);

	//	pairMassE  = CalculateMass(fEvt->fReconstructedFirst[i].fMomentum, (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum,5.11e-4, 5.11e-4);
	
	//Kt
	pairKt = pow(fEvt->fReconstructedFirst[i].fMomentum[0] + (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum[0],2.);
	pairKt+= pow(fEvt->fReconstructedFirst[i].fMomentum[1] + (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum[1],2.);
	pairKt = sqrt(pairKt)/2.;
	
        // Pair histogramming
     
	//dphis = CalculateDphiSatR12m(fEvt->fReconstructedFirst[i].fCharge, (fEvt+eventNumber)->fReconstructedSecond[j].sCharge, fieldsign,fEvt->fReconstructedFirst[i].fPt , (fEvt+eventNumber)->fReconstructedSecond[j].sPt, fEvt->fReconstructedFirst[i].fPhi, (fEvt+eventNumber)->fReconstructedSecond[j].sPhi); // 2 - 1
	//deta = (fEvt+eventNumber)->fReconstructedSecond[j].sEta - fEvt->fReconstructedFirst[i].fEta; // 2 - 1
	
	dphi = CalculateDPhiStar(fEvt->fReconstructedFirst[i].fCharge, (fEvt+eventNumber)->fReconstructedSecond[j].sCharge, fieldsign,fEvt->fReconstructedFirst[i].fPt , (fEvt+eventNumber)->fReconstructedSecond[j].sPt, fEvt->fReconstructedFirst[i].fPhi, (fEvt+eventNumber)->fReconstructedSecond[j].sPhi,0); // 2 - 1
	
	dphis = CalculateDPhiStar(fEvt->fReconstructedFirst[i].fCharge, (fEvt+eventNumber)->fReconstructedSecond[j].sCharge, fieldsign,fEvt->fReconstructedFirst[i].fPt , (fEvt+eventNumber)->fReconstructedSecond[j].sPt, fEvt->fReconstructedFirst[i].fPhi, (fEvt+eventNumber)->fReconstructedSecond[j].sPhi,fRadius); // 2 - 1

	deta   = CalculateDeltaEta(fEvt->fReconstructedFirst[i].fEta, (fEvt+eventNumber)->fReconstructedSecond[j].sEta);
	dtheta = CalculateDeltaTheta(fEvt->fReconstructedFirst[i].fTheta, (fEvt+eventNumber)->fReconstructedSecond[j].sTheta);
	//detas = fEvt->fReconstructedFirst[i].fEtaS-(fEvt+eventNumber)->fReconstructedFirst[j].fEtaS;
	
	// // Apply two-track cuts //RAMONA
	if (fkApplyTtc) {
	  //	  if (TMath::Abs(dphis)<fDphisMin && TMath::Abs(deta)<fDetasMin) continue;  
	  if (TMath::Abs(dphi)<fDphisMin && TMath::Abs(deta)<fDetasMin) continue;  
	}

	if (eventNumber==0) {//Same event pair histogramming

	  tSignP1             = SignP1;
	  tSignP2             = SignP2;                                  
	  tCentrality         = lcentrality;		                      
	  tDCAxyP1            = DCAxyP1;			                   
	  tDCAzP1             = DCAzP1;  			                   
	  tDCAxyP2            = DCAxyP2; 			                    
	  tDCAzP2             = DCAzP2;  			                   
	  tKtpair             = pairKt;			                  
	  tkStar              = pairKstar;		                 
	  tptP1               = ptP1;			                
	  tptP2               = ptP2;			                
	  tDEta               = deta; 			                
	  tDPhiStar           = dphis;			                    
	  tDPhi               = dphi;			                
	  tMassPair           = pairMass;			                    
	  tSphericity         = fSphericityvalue;		                      
	  tMassS              = massS;		                              
	  tDTheta             = dtheta;                                     
  
	  if(fReadMCTruth == kTRUE){
	    tMCtruepair    =  isMCvector;		
	    tMCSameMother  =  sameMother;		
	    tMCMotherP1    =  mcPDGMotherP1;//mcMotherBin;		
	    tMCMotherP2    =  mcPDGMotherP2;//mcMotherBin;		
	    tMCptcTypeP1   =  typeP1 ;		
	    tMCptcTypeP2   =  typeP2 ;		
	    tMCSameGM	   =  sameGrandMother;	
	    tMotherPDG	   =  mcGrandMotherBin;	
	    tpdgcodeP1 	   =  mcPDGcodeP1;//mcPDG1Bin;		
	    tpdgcodeP2	   =  mcPDGcodeP2;//mcPDG2Bin;		 
	    tKstarGen      =  pairKstarMC;            
	  }
	  
	  //	  fHistSparseSignal->Fill();  
	  if(ptP2<1.5)   
	    fHistSparseSignal->Fill();
	  else
	    if(massS>=1 )
	      fHistSparseSignal->Fill();
	}

	else {//Mixed-event pair histogramming
	  
	  tSignP1             = SignP1;
	  tSignP2             = SignP2;                                  
	  tCentrality         = lcentrality;		                      
	  tDCAxyP1            = DCAxyP1;			                   
	  tDCAzP1             = DCAzP1;  			                   
	  tDCAxyP2            = DCAxyP2; 			                    
	  tDCAzP2             = DCAzP2;  			                   
	  tKtpair             = pairKt;			                  
	  tkStar              = pairKstar;		                 
	  tptP1               = ptP1;			                
	  tptP2               = ptP2;			                
	  tDEta               = deta; 			                
	  tDPhiStar           = dphis;			                    
	  tDPhi               = dphi;			                
	  tMassPair           = pairMass;			                    
	  tSphericity         = fSphericityvalue;		                      
	  tMassS              = massS;		    	                     
	  tDTheta             = dtheta;                                     

	  if(fReadMCTruth == kTRUE){
	    tMCtruepair    =  isMCvector;		
	    tMCSameMother  =  sameMother;		
	    tMCMotherP1    =  mcPDGMotherP1;//mcMotherBin;		
	    tMCMotherP2    =  mcPDGMotherP2;//mcMotherBin;		
	    tMCptcTypeP1   =  typeP1 ;		
	    tMCptcTypeP2   =  typeP2 ;		
	    tMCSameGM	   =  sameGrandMother;	
	    tMotherPDG	   =  mcGrandMotherBin;	
	    tpdgcodeP1 	   =  mcPDGcodeP1;//mcPDG1Bin;		
	    tpdgcodeP2	   =  mcPDGcodeP2;//mcPDG2Bin;		 
	    tKstarGen      =  pairKstarMC;            
	  }
	  
	  // fHistSparseBkg->Fill();  
	  if(ptP2<1.5)   
	    fHistSparseBkg->Fill();
	  else
	    if(massS>=1 )
	      fHistSparseBkg->Fill();
        } //mixed
	
      } // second part

    }//end event loop

    if (evmultmixed!=0) multmixedcounted = kTRUE;
    
  } // first part
  
  if  (multmixedcounted) 
    fHistMultiplicityOfMixedEvent->Fill(evmultmixed);
  
}

//----------------------------------------------------------------------------------------------
//void AliAnalysisTaskDeuFlow2PC::DoPairshh (const Float_t lcentrality, int fieldsign) {
void AliAnalysisTaskDeuFlow2PC::DoPairshh (const Float_t lcentrality, int fieldsign, const Double_t fSphericityvalue) {

  //-----------
  double DCAxyP1  = -999. ;
  double DCAzP1   = -999. ;  
  double DCAxyP2  = -999. ; 
  double DCAzP2   = -999. ;  

  double  ptP1 = -999.;
  double  ptP2 = -999.;

  // Short_t chargeP1 = -999.;
  // Short_t chargeP2 = -999.;

  bool isP1  = kFALSE;
  bool isaP1 = kFALSE;
  bool isP2  = kFALSE;
  bool isaP2 = kFALSE;

  Int_t  SignP1 = -999;
  Int_t  SignP2 = -999;

  double dphi  = -999.;
  double dphis = -999.;
  double deta  = -999.;
  //  double detas = -999.;
  double dtheta = -999.;
  
  bool isMC1 = kFALSE;
  bool isMC2 = kFALSE;

  bool isMCvector = kFALSE;
  bool sameMother = kFALSE;
  bool sameGrandMother = kFALSE;
  
  Int_t mcMotherLabelP1 = -999;
  Int_t mcMotherLabelP2 = -999;

  // Int_t mcGrandMotherLabelP1 = -999;
  // Int_t mcGrandMotherLabelP2 = -999;
 
  Int_t typeP1 = -999;
  Int_t typeP2 = -999;
 
  Int_t mcPDGMotherP1 = 0;
  Int_t mcPDGMotherP2 = 0;
  //  Int_t mcMotherBin = 0;
 
  //  Int_t mcPDGGrandMother = 0;
  Int_t mcGrandMotherBin = 0;

  Int_t mcPDGcodeP1 = 0;
  Int_t mcPDGcodeP2 = 0;
  // Int_t mcPDG1Bin = 0;
  // Int_t mcPDG2Bin = 0;

  int evmultmixed = 0;
  bool multmixedcounted = kFALSE;
 
  double pairKstar   = 0.;
  double pairKstarMC = 0.;
  double pairMass  = 0.;
  //  double pairMassE = 0.;
  double pairKt    = 0.;

  double massS = 0;

  int idx1 = 0;
  int idx2 = 0;

  //fEvt->fReconstructedSecond[sCount].index = TMath::Abs(globaltrack->GetID());


  for (int i=0; i<fEvt->fNumberCandidateFirst; i++) {

    //    if (fEvt->fReconstructedFirst[i].doSkipOver) continue;
    
    DCAxyP1         = fEvt->fReconstructedFirst[i].fDCAxy;
    DCAzP1          = fEvt->fReconstructedFirst[i].fDCAz;
    ptP1            = fEvt->fReconstructedFirst[i].fPt;
    
    if(ptP1 < fMinPtForPrim)continue;
    
    //    chargeP1        = fEvt->fReconstructedFirst[i].fCharge;
    isMC1           = fEvt->fReconstructedFirst[i].isMCptc;
    mcMotherLabelP1 = fEvt->fReconstructedFirst[i].fMCmumIdx;
    typeP1          = fEvt->fReconstructedFirst[i].fMCcode;
    //  if(typeP1 < 0) continue;
    //    mcGrandMotherLabelP1 = fEvt->fReconstructedFirst[i].fMCgrandmumIdx;
    mcPDGcodeP1          = fEvt->fReconstructedFirst[i].fPDGcode;

    
    isP1  = fEvt->fReconstructedFirst[i].isP;
    isaP1 = fEvt->fReconstructedFirst[i].isaP;

    if(isP1) SignP1 = 1;
    else if (isaP1) SignP1 = -1;

    mcPDGMotherP1    = fEvt->fReconstructedFirst[i].fMCmumPDG;
    //    mcPDGGrandMother = fEvt->fReconstructedFirst[i].fMCgrandmumPDG;
    
    idx1 = fEvt->fReconstructedFirst[i].index;

    for (int eventNumber=0; eventNumber<fnEventsToMix+1; eventNumber++) { 
      // For same event pairs
      // if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateFirst)!=0.) evmultmixed++; 
      if (!multmixedcounted && eventNumber!=0 && ((fEvt+eventNumber)->fNumberCandidateSecond)!=0.) evmultmixed++; 

      
      for (int j=0; j<(fEvt+eventNumber)->fNumberCandidateSecond; j++) {

	/*// Ramona
	if ( (eventNumber == 0) && (j<=i)) continue; 
	
	// Instead of loop variables use im and jn to do the swapping in a clean way  
        Int_t tmpi = 0;
        // Pair ramdomization 
        if (eventNumber == 0) {
          if (gRandom->Rndm()>=0.5) { 
	    tmpi = i;
	    i = j;
	    j = tmpi;
	    //            im = j; jn = i;     
          }  
        }
	*/

	// cout<<"Ramona event number "<<eventNumber<<endl;
	// cout<<"eventNumber: "<<eventNumber<<" j: "<<j<<" i: "<<i<<" fEvt+eventNumber: "<<(fEvt+eventNumber)<<endl;
	
	//	if ((fEvt+eventNumber)->fReconstructedSecond[j].doSkipOver) continue;
        
	idx2 = (fEvt+eventNumber)->fReconstructedSecond[j].index;
	
	if ( (eventNumber == 0) && (j<=i) && idx1 == idx2) continue; 
	
	//	cout<<"Ramona Doing stuff Second Ptc"<<endl;
	
	DCAxyP2         = (fEvt+eventNumber)->fReconstructedSecond[j].sDCAxy;
	DCAzP2          = (fEvt+eventNumber)->fReconstructedSecond[j].sDCAz;
	ptP2            = (fEvt+eventNumber)->fReconstructedSecond[j].sPt;

	if(ptP1 == ptP2)continue;

	if(ptP2< fMinPtForSec)continue;
	//chargeP2        = (fEvt+eventNumber)->fReconstructedSecond[j].sCharge;
	isMC2           = (fEvt+eventNumber)->fReconstructedSecond[j].isMCptc;
	mcMotherLabelP2 = (fEvt+eventNumber)->fReconstructedSecond[j].sMCmumIdx;
	typeP2          = (fEvt+eventNumber)->fReconstructedSecond[j].sMCcode;
	//	if(typeP2 < 0) continue;
	//	mcGrandMotherLabelP2 = (fEvt+eventNumber)->fReconstructedSecond[j].sMCgrandmumIdx;
  	mcPDGcodeP2     = (fEvt+eventNumber)->fReconstructedSecond[j].sPDGcode;
	massS           = (fEvt+eventNumber)->fReconstructedSecond[j].sTOFMass;

	isP2  = (fEvt+eventNumber)->fReconstructedSecond[j].isP;
        isaP2 = (fEvt+eventNumber)->fReconstructedSecond[j].isaP;
	
	mcPDGMotherP2 = (fEvt+eventNumber)->fReconstructedSecond[j].sMCmumPDG;
	
	if(isP2) SignP2 = 1;
	else if (isaP2) SignP2 = -1;

	if(isMC1 && isMC2) isMCvector = kTRUE;
	else isMCvector = kFALSE;
	
	if(mcMotherLabelP1 == mcMotherLabelP2 && mcMotherLabelP1!=-999)sameMother = kTRUE;
	else sameMother = kFALSE;
	
	//Calculate k* for the pair
        pairKstar = CalculateKstar(fEvt->fReconstructedFirst[i].fMomentum, (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum,fPDGfirst, fPDGsecond);
	
	//mc kstar
	pairKstarMC = CalculateKstar(fEvt->fReconstructedFirst[i].fMomentumTruth, (fEvt+eventNumber)->fReconstructedSecond[j].sMomentumTruth,fPDGfirst, fPDGsecond);
	
	//Invariant Mass of the pair
	pairMass  = CalculateMass(fEvt->fReconstructedFirst[i].fMomentum, (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum,fPDGfirst, fPDGsecond);

	//	pairMassE  = CalculateMass(fEvt->fReconstructedFirst[i].fMomentum, (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum,5.11e-4, 5.11e-4);
	
	//Kt
	pairKt = pow(fEvt->fReconstructedFirst[i].fMomentum[0] + (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum[0],2.);
	pairKt+= pow(fEvt->fReconstructedFirst[i].fMomentum[1] + (fEvt+eventNumber)->fReconstructedSecond[j].sMomentum[1],2.);
	pairKt = sqrt(pairKt)/2.;
	
        // Pair histogramming
     
	//dphis = CalculateDphiSatR12m(fEvt->fReconstructedFirst[i].fCharge, (fEvt+eventNumber)->fReconstructedSecond[j].sCharge, fieldsign,fEvt->fReconstructedFirst[i].fPt , (fEvt+eventNumber)->fReconstructedSecond[j].sPt, fEvt->fReconstructedFirst[i].fPhi, (fEvt+eventNumber)->fReconstructedSecond[j].sPhi); // 2 - 1
	//deta = (fEvt+eventNumber)->fReconstructedSecond[j].sEta - fEvt->fReconstructedFirst[i].fEta; // 2 - 1
	
	dphi = CalculateDPhiStar(fEvt->fReconstructedFirst[i].fCharge, (fEvt+eventNumber)->fReconstructedSecond[j].sCharge, fieldsign,fEvt->fReconstructedFirst[i].fPt , (fEvt+eventNumber)->fReconstructedSecond[j].sPt, fEvt->fReconstructedFirst[i].fPhi, (fEvt+eventNumber)->fReconstructedSecond[j].sPhi,0); // 2 - 1
	
	dphis = CalculateDPhiStar(fEvt->fReconstructedFirst[i].fCharge, (fEvt+eventNumber)->fReconstructedSecond[j].sCharge, fieldsign,fEvt->fReconstructedFirst[i].fPt , (fEvt+eventNumber)->fReconstructedSecond[j].sPt, fEvt->fReconstructedFirst[i].fPhi, (fEvt+eventNumber)->fReconstructedSecond[j].sPhi,fRadius); // 2 - 1

	deta   = CalculateDeltaEta(fEvt->fReconstructedFirst[i].fEta, (fEvt+eventNumber)->fReconstructedSecond[j].sEta);
	dtheta = CalculateDeltaTheta(fEvt->fReconstructedFirst[i].fTheta, (fEvt+eventNumber)->fReconstructedSecond[j].sTheta);
	//detas = fEvt->fReconstructedFirst[i].fEtaS-(fEvt+eventNumber)->fReconstructedFirst[j].fEtaS;
	
	// // Apply two-track cuts //RAMONA
	if (fkApplyTtc) {
	  //	  if (TMath::Abs(dphis)<fDphisMin && TMath::Abs(deta)<fDetasMin) continue;  
	  if (TMath::Abs(dphi)<fDphisMin && TMath::Abs(deta)<fDetasMin) continue;  
	}

	if (eventNumber==0) {//Same event pair histogramming

	  tSignP1             = SignP1;
	  tSignP2             = SignP2;                                  
	  tCentrality         = lcentrality;		                      
	  tDCAxyP1            = DCAxyP1;			                   
	  tDCAzP1             = DCAzP1;  			                   
	  tDCAxyP2            = DCAxyP2; 			                    
	  tDCAzP2             = DCAzP2;  			                   
	  tKtpair             = pairKt;			                  
	  tkStar              = pairKstar;		                 
	  tptP1               = ptP1;			                
	  tptP2               = ptP2;			                
	  tDEta               = deta; 			                
	  tDPhiStar           = dphis;			                    
	  tDPhi               = dphi;			                
	  tMassPair           = pairMass;			                    
	  tSphericity         = fSphericityvalue;		                      
	  tMassS              = massS;		                              
	  tDTheta             = dtheta;                                     
  
	  if(fReadMCTruth == kTRUE){
	    tMCtruepair    =  isMCvector;		
	    tMCSameMother  =  sameMother;		
	    tMCMotherP1    =  mcPDGMotherP1;//mcMotherBin;		
	    tMCMotherP2    =  mcPDGMotherP2;//mcMotherBin;		
	    tMCptcTypeP1   =  typeP1 ;		
	    tMCptcTypeP2   =  typeP2 ;		
	    tMCSameGM	   =  sameGrandMother;	
	    tMotherPDG	   =  mcGrandMotherBin;	
	    tpdgcodeP1 	   =  mcPDGcodeP1;//mcPDG1Bin;		
	    tpdgcodeP2	   =  mcPDGcodeP2;//mcPDG2Bin;		 
	    tKstarGen      =  pairKstarMC;            
	  }
	  
	  //	  fHistSparseSignal->Fill();  
	  //	  if(ptP2<1.5)   
	  fHistSparseSignal->Fill();
	  // else
	  //   if(massS>=1)
	  //     fHistSparseSignal->Fill();
	}

	else {//Mixed-event pair histogramming
	  
	  tSignP1             = SignP1;
	  tSignP2             = SignP2;                                  
	  tCentrality         = lcentrality;		                      
	  tDCAxyP1            = DCAxyP1;			                   
	  tDCAzP1             = DCAzP1;  			                   
	  tDCAxyP2            = DCAxyP2; 			                    
	  tDCAzP2             = DCAzP2;  			                   
	  tKtpair             = pairKt;			                  
	  tkStar              = pairKstar;		                 
	  tptP1               = ptP1;			                
	  tptP2               = ptP2;			                
	  tDEta               = deta; 			                
	  tDPhiStar           = dphis;			                    
	  tDPhi               = dphi;			                
	  tMassPair           = pairMass;			                    
	  tSphericity         = fSphericityvalue;		                      
	  tMassS              = massS;		    	                     
	  tDTheta             = dtheta;                                     

	  if(fReadMCTruth == kTRUE){
	    tMCtruepair    =  isMCvector;		
	    tMCSameMother  =  sameMother;		
	    tMCMotherP1    =  mcPDGMotherP1;//mcMotherBin;		
	    tMCMotherP2    =  mcPDGMotherP2;//mcMotherBin;		
	    tMCptcTypeP1   =  typeP1 ;		
	    tMCptcTypeP2   =  typeP2 ;		
	    tMCSameGM	   =  sameGrandMother;	
	    tMotherPDG	   =  mcGrandMotherBin;	
	    tpdgcodeP1 	   =  mcPDGcodeP1;//mcPDG1Bin;		
	    tpdgcodeP2	   =  mcPDGcodeP2;//mcPDG2Bin;		 
	    tKstarGen      =  pairKstarMC;            
	  }
	  
	  fHistSparseBkg->Fill();  
	  // if(ptP2<1.5)   
	  //   fHistSparseBkg->Fill();
	  // else
	  //   if(massS>=1)
	  //     fHistSparseBkg->Fill();
        } //mixed
	
      } // second part

    }//end event loop

    if (evmultmixed!=0) multmixedcounted = kTRUE;
    
  } // first part
  
  if  (multmixedcounted) 
    fHistMultiplicityOfMixedEvent->Fill(evmultmixed);
  
  //return;
  
}

//-----------------------------------------------------------------------------------------------

double AliAnalysisTaskDeuFlow2PC::CalculateKstar(double momentum1[3], double momentum2[3], double mass1, double mass2) { // Jai S

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
double AliAnalysisTaskDeuFlow2PC::CalculateMass(double momentum1[3], double momentum2[3], double mass1, double mass2) { // Jai S

  // Calculate Invariant Mass 
  TLorentzVector  vP1,vP2,vSum;
  
  vP1.SetXYZM(momentum1[0],momentum1[1],momentum1[2],mass1); 
  vP2.SetXYZM(momentum2[0],momentum2[1],momentum2[2],mass2);       
  vSum=vP1+vP2;

  double mass = vSum.M();
 
  return mass;
}

//-----------------------------------------------------------------------------------------------
double AliAnalysisTaskDeuFlow2PC::CalculateDphiSatR12m(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2) { // AliFemto framework AliFemtoUser/AliFemtoPairCutRadialDistance.cxx + Dhevan not consistent?
  /*
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
  return deltaphi;//dps;
  */
  //  cout<<" Dphi "<<dps<<" Dhevan "<<deltaphi<<endl;


  //from mariella
  //analitical funcion
  
  double rad = 1.2;
  
  double afsi1b = 0.075*chg1*magSign*rad/ptv1; // 0.07510020733 = - 0.3 (= e in Heaviside-Lorentz units) *0.5 (= B in T) /2 (see later for the -), pT in GeV/c
  double afsi2b = 0.075*chg2*magSign*rad/ptv2;
  
  if (fabs(afsi1b) >=1.) return 9999.; // angle is pi/2 or not defined --> dont cut
  if (fabs(afsi2b) >=1.) return 9999.; // MN modified these two lines returning 9999 and not kTRUE
  double dps = phi2 - phi1  + TMath::ASin(afsi1b) -TMath::ASin(afsi2b); // - sign of e is outside Mariella
  
  dps = TVector2::Phi_mpi_pi(dps);
  
  return dps;

}


//-----------------------------------------------------------------------------------------------
double AliAnalysisTaskDeuFlow2PC::CalculateDPhiStar(Short_t chg1, Short_t chg2, Int_t magSign, Double_t ptv1, Double_t ptv2, Double_t phi1, Double_t phi2,Double_t rad) { //AliFemtoUser/AliFemtoPairCutDetaDphi.h

  const Double_t unit_factor = 0.299792458 / 2.0;
  const Double_t b_field = 0.5006670488586 * magSign;

  Double_t  shift1 = TMath::ASin(unit_factor * chg1 * b_field * rad / ptv1);
  Double_t  shift2 = TMath::ASin(unit_factor * chg2 * b_field * rad / ptv2);

  double dps = (phi1 + shift1) - (phi2 + shift2);
  
  //  dps = TVector2::Phi_mpi_pi(dps); //to be checked

  return dps; //deltaphi;
  
}
//_______________________________________________________________

Double_t AliAnalysisTaskDeuFlow2PC::CalculateDeltaEta( Double_t eta1, Double_t eta2 )  {   //AliFemtoUser/AliFemtoPairCutDetaDphi.h
  const double deta = eta2 - eta1;
  return deta;
}
//_______________________________________________________________
Double_t AliAnalysisTaskDeuFlow2PC::CalculateDeltaTheta( Double_t theta1, Double_t theta2 )  {  
  const double dtheta = theta2 - theta1;
  return dtheta;
}

//-----------------------------------------------------------------------------------------------
double AliAnalysisTaskDeuFlow2PC::CalculateDphiSatR12m(Double_t pos1SftR125[3], Double_t pos2SftR125[3]) { // Hans B

  // Returns delta phi star at R = 1.2 m

  const Float_t distSft = TMath::Sqrt(TMath::Power(pos1SftR125[0] - pos2SftR125[0],2)
				      + TMath::Power(pos1SftR125[1] - pos2SftR125[1],2));
  return 2.0 * TMath::ATan(distSft/2./(125.)); 

}

//-----------------------------------------------------------------------------------------------

void AliAnalysisTaskDeuFlow2PC::SetSftPosR125(AliVTrack *track, const Float_t bfield, Double_t priVtx[3], Double_t posSftR125[3] ) {  // Hans B
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
Double_t AliAnalysisTaskDeuFlow2PC::ThetaS( Double_t posSftR125[3] ) const { // Hans B
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
Double_t AliAnalysisTaskDeuFlow2PC::EtaS( Double_t posSftR125[3] ) const {  // Hans B
  // Returns the corresponding eta of a pri. part. 
  // with this particles pos at R=1.25m

  // http://en.wikipedia.org/wiki/Pseudorapidity
  // ? = -ln[ tan(?/2)]
  // printf("z: %+04.0f, thetaS %+03.2f etaS %+1.2f\n"
  // ,fXshifted[2][2],ThetaS(),-TMath::Log( TMath::Tan(ThetaS()/2.) ));
  return -TMath::Log( TMath::Tan(ThetaS(posSftR125 )/2.) );
}



//_________________________________________________________________
Double_t AliAnalysisTaskDeuFlow2PC::CalculateSphericityofEvent(AliAODEvent *aodEvent)
{ //from Oliver
  Double_t Pt_tot = 0.; //total Pt of all protons and v0s in the event

  Double_t S00 = 0.; //Elements of the sphericity matrix
  Double_t S11 = 0.;
  Double_t S10 = 0.;
  
  Int_t NumOfTracks = aodEvent->GetNumberOfTracks();
  if(NumOfTracks<3) return -9999.;//if already at this point not enough tracks are in the event -> return
  

  Int_t NTracks = 0;
  for(Int_t iTrack=0;iTrack<NumOfTracks;iTrack++)
    {
      AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(iTrack));
      
      if(!aodtrack->TestFilterBit(128)) continue;
      
      Double_t Pt = aodtrack->Pt();
      //Double_t Phi = aodtrack->Phi();
      Double_t Px = aodtrack->Px();
      Double_t Py = aodtrack->Py();
      Double_t eta = aodtrack->Eta();
      

      if(!(eta>-0.8 && eta<0.8)) continue;
      if(Pt<0.5) continue;

      Pt_tot += Pt;
      
      S00 += Px*Px/Pt;
      S11 += Py*Py/Pt;
      S10 += Px*Py/Pt;
      NTracks++;
    }
  
  if(NTracks<3) return -9999.;//new flag: check
  
  //normalize to total Pt to obtain a linear form:
  if(Pt_tot == 0.) return -9999.;
  S00 /= Pt_tot;
  S11 /= Pt_tot;
  S10 /= Pt_tot;
  
  //Calculate the trace of the sphericity matrix:
  Double_t T = S00+S11;
  
  //Calculate the determinant of the sphericity matrix:
  Double_t D = S00*S11 - S10*S10;//S10 = S01

  //Calculate the eigenvalues of the sphericity matrix:
  Double_t lambda1 = 0.5*(T + TMath::Sqrt(T*T - 4.*D));
  Double_t lambda2 = 0.5*(T - TMath::Sqrt(T*T - 4.*D));
  
  if((lambda1 + lambda2) == 0.) return -9999.;
  
  Double_t ST = -1.;
  
  if(lambda2>lambda1)
    {
      ST = 2.*lambda1/(lambda1+lambda2);
    }
  else
    {
      ST = 2.*lambda2/(lambda1+lambda2);
    }

  return ST;
}


//--------------------------------------------------- Methods From AliFemtoESDTrackCut.cxx
bool AliAnalysisTaskDeuFlow2PC::IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP)
{
  //  if(TMath::Abs(nsigmaTPCE)<3 && TMath::Abs(nsigmaTPCPi)>3 && TMath::Abs(nsigmaTPCK)>3 && TMath::Abs(nsigmaTPCP)>3)
  if(TMath::Abs(nsigmaTPCE)<3)
     return false;
  else
     return true;
}

//----------------------------------------------------------

bool AliAnalysisTaskDeuFlow2PC::IsPionNSigma(double mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  //sligly changed w.r.t. the original
  
  return false;

  // if(mom<0.65){   
  //   //      if(nsigmaTOFPi<-999.)
  //   if(nsigmaTOFPi==10)
  //     {
  // 	//use TPC only
  // 	if(mom<0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
  // 	else if(mom<0.5 && mom>=0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
  // 	else if(mom>=0.5 && TMath::Abs(nsigmaTPCPi)<2.0) return true;
  // 	else return false;
  //     } 
 
  //   else if(TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<3.0) return true; //TPC+TOF
  // }
  // //else if(nsigmaTOFPi<-10.) //p > 0.65 + no tof == kfalse
  // else if(mom>0.65 && nsigmaTOFPi>3) //p > 0.65 + no tof == kfalse
  //   {
  //     return false;
  //   }
  // else if(mom<1.5 && TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<5.0) return true;
  // else if(mom>=1.5 && TMath::Abs(nsigmaTOFPi)<2.0 && TMath::Abs(nsigmaTPCPi)<5.0) return true;
  
  // else 
  //   return false;
}
/*

//----------------------------------------------------------
bool AliAnalysisTaskDeuFlow2PC::IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{
  if (fNsigmaTPCTOF) {
    if (mom > 0.5) {
      //        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
      if (TMath::Hypot( nsigmaTOFK, nsigmaTPCK ) < fNsigma)
	return true;
    }
    else {
      if (TMath::Abs(nsigmaTPCK) < fNsigma)
	return true;
    }
  }
  else {

    if(mom<0.4)
      {
	if(nsigmaTOFK<-999.)
	  {
	    if(TMath::Abs(nsigmaTPCK)<2.0) return true;
	  }
	else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
      }
    else if(mom>=0.4 && mom<=0.6)
      {
	if(nsigmaTOFK<-999.)
	  {
	    if(TMath::Abs(nsigmaTPCK)<2.0) return true;
	  }
	else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
      }
    else if(nsigmaTOFK<-999.)
      {
	return false;
      }
    else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
  }
  return false;
}
//----------------------------------------------------------

bool AliAnalysisTaskDeuFlow2PC::IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{
  if (fNsigmaTPCTOF) {
    if (mom > 0.5) {
//        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < fNsigma)
            return true;
    } else if (TMath::Abs(nsigmaTPCP) < fNsigma) {
      return true;
    }
  }
  else if (fNsigmaTPConly) {
    if (TMath::Abs(nsigmaTPCP) < fNsigma)
      return true;
  }
  else {
    if (mom > 0.8 && mom < 2.5) {
      if ( TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 3.0)
        return true;
    }
    else if (mom > 2.5) {
      if ( TMath::Abs(nsigmaTPCP) < 3.0 && TMath::Abs(nsigmaTOFP) < 2.0)
        return true;
    }
    else {
      if (TMath::Abs(nsigmaTPCP) < 3.0)
        return true;
    }
  }

  return false;
}
*/

//-----------------------------------------------------------------------------------------------
void AliAnalysisTaskDeuFlow2PC::Terminate(const Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  if (!GetOutputData(0)) return;
}

