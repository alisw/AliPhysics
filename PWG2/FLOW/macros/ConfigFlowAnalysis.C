//
// configure macro for AliAnalysisTaskFlowEvent
//

#include "TList.h"
#include "TObjArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "AliCFManager.h"
#include "AliCFTrackKineCuts.h"
#include "AliCFParticleGenCuts.h"
#include "AliCFAcceptanceCuts.h"
#include "AliCFTrackKineCuts.h"
#include "AliCFTrackQualityCuts.h"
#include "AliCFTrackIsPrimaryCuts.h"
#include "AliCFTrackCutPid.h"
#include "AliAnalysisTaskFlowEvent.h"


// SETTING THE CUTS

// For integrated flow
const Double_t ptmin1 = 0.0;
const Double_t ptmax1 = 10.0;
const Double_t ymin1  = -1.;
const Double_t ymax1  = 1.;
const Int_t mintrackrefsTPC1 = 2;
const Int_t mintrackrefsITS1 = 3;
const Int_t charge1 = 1; //do not use
Bool_t UsePIDIntegratedFlow = kTRUE;
const Int_t PDG1 = 211;
const Int_t minclustersTPC1 = 50;
const Int_t maxnsigmatovertex1 = 3;

// For differential flow
const Double_t ptmin2 = 0.0;
const Double_t ptmax2 = 10.0;
const Double_t ymin2  = -1.;
const Double_t ymax2  = 1.;
const Int_t mintrackrefsTPC2 = 2;
const Int_t mintrackrefsITS2 = 3;
const Int_t charge2 = 1; //do not use
Bool_t UsePIDDifferentialFlow = kTRUE;
const Int_t PDG2 = 321;
const Int_t minclustersTPC2 = 50;
const Int_t maxnsigmatovertex2 = 3;

void ConfigFlowAnalysis(AliAnalysisTaskFlowEvent* fEventTask)
{

  //Create cuts using correction framework

  Bool_t QA = fEventTask->GetQAOn();
  //Set TList for the QA histograms
  TList* qaIntFE  = new TList(); 
  TList* qaDiffFE = new TList();
   

 //############# cuts on MC
 AliCFTrackKineCuts* mcKineCuts1 = new AliCFTrackKineCuts("mcKineCuts1","MC-level kinematic cuts");
 mcKineCuts1->SetPtRange(ptmin1,ptmax1);
 mcKineCuts1->SetRapidityRange(ymin1,ymax1);
 //mcKineCuts1->SetChargeMC(charge1);  //check what this cut does
 if (QA) { 
   mcKineCuts1->SetQAOn(qaIntFE);
 }

 AliCFTrackKineCuts* mcKineCuts2 = new AliCFTrackKineCuts("mcKineCuts2","MC-level kinematic cuts");
 mcKineCuts2->SetPtRange(ptmin2,ptmax2);
 mcKineCuts2->SetRapidityRange(ymin2,ymax2);
 //mcKineCuts2->SetChargeMC(charge2);  //check what this cut does
 if (QA) { 
   mcKineCuts2->SetQAOn(qaDiffFE);
 }

 AliCFParticleGenCuts* mcGenCuts1 = new AliCFParticleGenCuts("mcGenCuts1","MC particle generation cuts for integrated flow");
 mcGenCuts1->SetRequireIsPrimary();
 if (UsePIDIntegratedFlow) {mcGenCuts1->SetRequirePdgCode(PDG1);}
 if (QA) { 
   mcGenCuts1->SetQAOn(qaIntFE);
 }

 AliCFParticleGenCuts* mcGenCuts2 = new AliCFParticleGenCuts("mcGenCuts2","MC particle generation cuts for differential flow");
 mcGenCuts2->SetRequireIsPrimary();
 if (UsePIDDifferentialFlow) {mcGenCuts2->SetRequirePdgCode(PDG2);}
 if (QA) { 
   mcGenCuts2->SetQAOn(qaDiffFE);
 }

 //############# Acceptance Cuts  
 AliCFAcceptanceCuts *mcAccCuts1 = new AliCFAcceptanceCuts("mcAccCuts1","MC acceptance cuts");
 mcAccCuts1->SetMinNHitITS(mintrackrefsITS1);
 mcAccCuts1->SetMinNHitTPC(mintrackrefsTPC1);
 if (QA) { 
   mcAccCuts1->SetQAOn(qaIntFE);
 }

 AliCFAcceptanceCuts *mcAccCuts2 = new AliCFAcceptanceCuts("mcAccCuts2","MC acceptance cuts");
 mcAccCuts2->SetMinNHitITS(mintrackrefsITS2);
 mcAccCuts2->SetMinNHitTPC(mintrackrefsTPC2);
 if (QA) { 
   mcAccCuts2->SetQAOn(qaDiffFE);
 }
//############# Rec-Level kinematic cuts
 AliCFTrackKineCuts *recKineCuts1 = new AliCFTrackKineCuts("recKineCuts1","rec-level kine cuts");
 recKineCuts1->SetPtRange(ptmin1,ptmax1);
 recKineCuts1->SetRapidityRange(ymin1,ymax1);
 //recKineCuts1->SetChargeRec(charge1); //selects only particles with charge == charge1
 if (QA) { 
   recKineCuts1->SetQAOn(qaIntFE);
 }

 AliCFTrackKineCuts *recKineCuts2 = new AliCFTrackKineCuts("recKineCuts2","rec-level kine cuts");
 recKineCuts2->SetPtRange(ptmin2,ptmax2);
 recKineCuts2->SetRapidityRange(ymin2,ymax2);
 //recKineCuts2->SetChargeRec(charge2); //selects only particles with charge == charge2
 if (QA) { 
   recKineCuts2->SetQAOn(qaDiffFE);
 }

 AliCFTrackQualityCuts *recQualityCuts1 = new AliCFTrackQualityCuts("recQualityCuts1","rec-level quality cuts");
 recQualityCuts1->SetMinNClusterTPC(minclustersTPC1);
 recQualityCuts1->SetStatus(AliESDtrack::kITSrefit);
 if (QA) { 
   recQualityCuts1->SetQAOn(qaIntFE);
 }
AliCFTrackQualityCuts *recQualityCuts2 = new AliCFTrackQualityCuts("recQualityCuts2","rec-level quality cuts");
 recQualityCuts2->SetMinNClusterTPC(minclustersTPC2);
 recQualityCuts2->SetStatus(AliESDtrack::kITSrefit);
 if (QA) { 
   recQualityCuts2->SetQAOn(qaDiffFE);
 }

 AliCFTrackIsPrimaryCuts *recIsPrimaryCuts1 = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts1","rec-level isPrimary cuts");
 recIsPrimaryCuts1->SetMaxNSigmaToVertex(maxnsigmatovertex1);
 if (QA) { 
   recIsPrimaryCuts1->SetQAOn(qaIntFE);
 }

 AliCFTrackIsPrimaryCuts *recIsPrimaryCuts2 = new AliCFTrackIsPrimaryCuts("recIsPrimaryCuts2","rec-level isPrimary cuts");
 recIsPrimaryCuts2->SetMaxNSigmaToVertex(maxnsigmatovertex2);
 if (QA) { 
   recIsPrimaryCuts2->SetQAOn(qaDiffFE);
 }

 int n_species = AliPID::kSPECIES ;
 Double_t* prior = new Double_t[n_species];

prior[0] = 0.0244519 ;
 prior[1] = 0.0143988 ;
 prior[2] = 0.805747  ;
 prior[3] = 0.0928785 ;
 prior[4] = 0.0625243 ;

 AliCFTrackCutPid* cutPID1 = NULL;
 if(UsePIDIntegratedFlow) {
   cutPID1 = new AliCFTrackCutPid("cutPID1","ESD_PID for integrated flow") ;
   cutPID1->SetPriors(prior);
   cutPID1->SetProbabilityCut(0.0);
   cutPID1->SetDetectors("TPC TOF");
   switch(TMath::Abs(PDG1)) {
   case 11   : cutPID1->SetParticleType(AliPID::kElectron, kTRUE); break;
   case 13   : cutPID1->SetParticleType(AliPID::kMuon    , kTRUE); break;
   case 211  : cutPID1->SetParticleType(AliPID::kPion    , kTRUE); break;
   case 321  : cutPID1->SetParticleType(AliPID::kKaon    , kTRUE); break;
   case 2212 : cutPID1->SetParticleType(AliPID::kProton  , kTRUE); break;
   default   : printf("UNDEFINED PID\n"); break;
   }
   if (QA) { 
     cutPID1->SetQAOn(qaIntFE); 
   }
 }
		  
AliCFTrackCutPid* cutPID2 = NULL;
 if (UsePIDDifferentialFlow) {
   cutPID2 = new AliCFTrackCutPid("cutPID2","ESD_PID for differential flow") ;
   cutPID2->SetPriors(prior);
   cutPID2->SetProbabilityCut(0.0);
   cutPID2->SetDetectors("TPC TOF");
   switch(TMath::Abs(PDG2)) {
   case 11   : cutPID2->SetParticleType(AliPID::kElectron, kTRUE); break;
   case 13   : cutPID2->SetParticleType(AliPID::kMuon    , kTRUE); break;
   case 211  : cutPID2->SetParticleType(AliPID::kPion    , kTRUE); break;
   case 321  : cutPID2->SetParticleType(AliPID::kKaon    , kTRUE); break;
   case 2212 : cutPID2->SetParticleType(AliPID::kProton  , kTRUE); break;
   default   : printf("UNDEFINED PID\n"); break;
   }
   if (QA) { 
     cutPID2->SetQAOn(qaIntFE);
   }
 }

 printf("CREATE MC KINE CUTS\n");
 TObjArray* mcList1 = new TObjArray(0);
 mcList1->AddLast(mcKineCuts1);
 mcList1->AddLast(mcGenCuts1);

 TObjArray* mcList2 = new TObjArray(0);
 mcList2->AddLast(mcKineCuts2);
 mcList2->AddLast(mcGenCuts2);

printf("CREATE ACCEPTANCE CUTS\n");
 TObjArray* accList1 = new TObjArray(0) ;
 accList1->AddLast(mcAccCuts1);

 TObjArray* accList2 = new TObjArray(0) ;
 accList2->AddLast(mcAccCuts2);

 printf("CREATE RECONSTRUCTION CUTS\n");
 TObjArray* recList1 = new TObjArray(0) ;
 recList1->AddLast(recKineCuts1);
 recList1->AddLast(recQualityCuts1);
 recList1->AddLast(recIsPrimaryCuts1);

 TObjArray* recList2 = new TObjArray(0) ;
 recList2->AddLast(recKineCuts2);
 recList2->AddLast(recQualityCuts2);
 recList2->AddLast(recIsPrimaryCuts2);

 printf("CREATE PID CUTS\n");
 TObjArray* fPIDCutList1 = new TObjArray(0) ;
 if(UsePIDIntegratedFlow) {fPIDCutList1->AddLast(cutPID1);}

 TObjArray* fPIDCutList2 = new TObjArray(0) ;
 if (UsePIDDifferentialFlow)  {fPIDCutList2->AddLast(cutPID2);}

printf("CREATE INTERFACE AND CUTS\n");
 AliCFManager* cfmgr1 = new AliCFManager();
 cfmgr1->SetNStepParticle(4); 
 cfmgr1->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList1);
 cfmgr1->SetParticleCutsList(AliCFManager::kPartAccCuts,accList1);
 cfmgr1->SetParticleCutsList(AliCFManager::kPartRecCuts,recList1);
 cfmgr1->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList1);

 AliCFManager* cfmgr2 = new AliCFManager();
 cfmgr2->SetNStepParticle(4); 
 cfmgr2->SetParticleCutsList(AliCFManager::kPartGenCuts,mcList2);
 cfmgr2->SetParticleCutsList(AliCFManager::kPartAccCuts,accList2);
 cfmgr2->SetParticleCutsList(AliCFManager::kPartRecCuts,recList2);
 cfmgr2->SetParticleCutsList(AliCFManager::kPartSelCuts,fPIDCutList2);
 
 fEventTask->SetQAList1(qaIntFE);
 fEventTask->SetQAList2(qaDiffFE);
 fEventTask->SetCFManager1(cfmgr1);
 fEventTask->SetCFManager2(cfmgr2);


}

