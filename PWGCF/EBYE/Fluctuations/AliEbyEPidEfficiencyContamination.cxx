/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: ALICE Offline.                                                 *
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


//=========================================================================//
//             AliEbyE Analysis for Net-Particle Higher Moment study       //
//                           Nirbhay K. Behera                             //
//                           nbehera@cern.ch                               //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                                                                         //
//                        (Last Modified 2018/08/27)                       //
//                 Dealing with Wide pT Window Modified to ESDs            //
//Some parts of the code are taken from J. Thaeder/ M. Weber NetParticle   //
//analysis task.                                                           //
//=========================================================================//

//ver: 2018/08/27 support only for proton

#include <Riostream.h>
#include "TList.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"

#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliAODHeader.h"
#include "AliAODpidUtil.h"
#include "AliEventCuts.h"
#include "AliMultSelection.h"
//#include "AliHelperPID.h"
#include "AliAnalysisTaskSE.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliMCEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliEbyEPidEfficiencyContamination.h"

using std::endl;
using std::cout;

ClassImp(AliEbyEPidEfficiencyContamination)

//-----------------------------------------------------------------------
AliEbyEPidEfficiencyContamination::AliEbyEPidEfficiencyContamination() 
: AliAnalysisTaskSE(), 
  fThnList(NULL),
  fInputHandler(NULL),
  fMCEventHandler(NULL),
  fVevent(NULL),
  fArrayMC(NULL),
  fESDtrackCuts(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fEventCuts(NULL),
  fRun("LHC10h"),
  fCentralityEstimator("V0M"),
  
  fAODtrackCutBit(128),
  fVxMax(3.), 
  fVyMax(3.), 
  fVzMax(10.), 
  fPtMin(0.15),   
  fPtMax(3.15), 
  fEtaMin(-1.), 
  fEtaMax(1.),
  fNptBins(22),
  fDcaXy(10.),
  fDcaZ(10.),  
  fNcrossRows(80),
  fChi2NDF(4),

  
  fIsMC(kFALSE),
  fIsAOD(kFALSE),
  fIsRapCut(kFALSE),
  fTotP(kFALSE),

  fNTracks(0),
  fCentrality(-1),
  fEventCounter(NULL), 
  fHistCent(0x0), 

  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fPidType(1),
  fMcPid(211),
  fPidStrategy(0),
  fNSigmaMaxITS(3.),
  fNSigmaMaxTPC(3.),
  fNSigmaMaxTOF(3.),
  fParticleSpecies(AliPID::kPion),
  
  fPtBinNplusNminusCh(NULL),
  fPtBinNplusNminusChTruth(NULL)
{ 
  
  for (Int_t i = 0; i < 2; i++) {
    fHitCentRec[i] = NULL;
    fHitCentGen[i] = NULL;
    fCentPtEtaPhiThnGen[i] = NULL;
    fCentPtEtaPhiThnRec[i] = NULL;
    fCentPtEtaPhiThnRecPrim[i] = NULL;
    fCentPtEtaPhiThnSec[i] = NULL;
    fCentPtEtaPhiThnMat[i] = NULL;
    fCentPtEtaPhiThnMisId[i] = NULL;
  }
}

//-----------------------------------------------------------------------
AliEbyEPidEfficiencyContamination::AliEbyEPidEfficiencyContamination( const char *name ) 
  : AliAnalysisTaskSE(name), 
    fThnList(NULL),
    fInputHandler(NULL),
    fMCEventHandler(NULL),
    fVevent(NULL),
    fArrayMC(NULL),
    fESDtrackCuts(NULL),
    fMCEvent(NULL),
    fMCStack(NULL),
    fEventCuts(NULL),
    fRun("LHC10h"),
    fCentralityEstimator("V0M"),
    
    fAODtrackCutBit(128),
    fVxMax(3.), 
    fVyMax(3.), 
    fVzMax(10.), 
    fPtMin(0.15),   
    fPtMax(3.15), 
    fEtaMin(-1.), 
    fEtaMax(1.),
    fNptBins(22),
    fDcaXy(10.),
    fDcaZ(10.),
    fNcrossRows(80),
    fChi2NDF(4),
    
    fIsMC(kFALSE),
    fIsAOD(kFALSE),
    fIsRapCut(kFALSE),
    fTotP(kFALSE),
    
    fNTracks(0),
    fCentrality(-1),
    fEventCounter(NULL), 
    fHistCent(NULL), 
   
    fPIDResponse(0x0),
    fPIDCombined(0x0),
    fPidType(1),
    fMcPid(211),
    fPidStrategy(0),
    fNSigmaMaxITS(4.),
    fNSigmaMaxTPC(4.),
    fNSigmaMaxTOF(4.),
    fParticleSpecies(AliPID::kPion),
   
    fPtBinNplusNminusCh(NULL),
    fPtBinNplusNminusChTruth(NULL)
{  
  for (Int_t i = 0; i < 2; i++) {
    fHitCentRec[i] = NULL;
    fHitCentGen[i] = NULL;
    fCentPtEtaPhiThnGen[i] = NULL;
    fCentPtEtaPhiThnRec[i] = NULL;
    fCentPtEtaPhiThnRecPrim[i] = NULL;
    fCentPtEtaPhiThnSec[i] = NULL;
    fCentPtEtaPhiThnMat[i] = NULL;
    fCentPtEtaPhiThnMisId[i] = NULL;
  }
  DefineOutput(1, TList::Class()); 
}

//---------------------------------------------------------------------------------
AliEbyEPidEfficiencyContamination::~AliEbyEPidEfficiencyContamination() {
  
   // Default destructor
  if( fThnList ) delete fThnList;
  if( fEventCuts ) delete fEventCuts;
  if( fESDtrackCuts ) delete fESDtrackCuts;
  
}

//---------------------------------------------------------------------------------
void AliEbyEPidEfficiencyContamination::UserCreateOutputObjects(){
  
  fInputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!fInputHandler){
    AliError("No PID response task found !!");
  }
  
  fPIDResponse =  dynamic_cast<AliInputEventHandler *>(fInputHandler)->GetPIDResponse();
  if(!fPIDResponse){
    AliError("No PID response task found !!");
  }
  
  if(fRun == "LHC15o"){
    fEventCuts = new AliEventCuts();
  }
  
  fThnList = new TList();
  fThnList->SetOwner(kTRUE);
  
  if (!fIsAOD) {
    
    if(!fESDtrackCuts) {

      fESDtrackCuts = new AliESDtrackCuts;
      fESDtrackCuts->SetName(Form("NetPesdCut_%d",fPidType));
      fESDtrackCuts->SetEtaRange(fEtaMin, fEtaMax);
      fESDtrackCuts->SetPtRange(0.1, 1e10);
      
      // TPC track
      fESDtrackCuts->SetMinNCrossedRowsTPC( fNcrossRows );// default = 80, sys: 60 and 100 
      fESDtrackCuts->SetMaxChi2PerClusterTPC( fChi2NDF );//default = 4, sys: 3.5, 4.5
      fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);//
      fESDtrackCuts->SetRequireTPCRefit(kTRUE);//
      fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
      
      if (fPidType == 0) {
	fESDtrackCuts->SetMaxDCAToVertexXY( fDcaXy ); //default=2. sys: 1.5, 2.5
      } else {
	if( fDcaXy < 2) fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0156+0.03/pt^1.01"); //tight cut
	else if( fDcaXy > 2) fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0208+0.04/pt^1.01");//losse cut
	else if( fDcaXy == 2) fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.035/pt^1.01");//default
      }
      
      
      fESDtrackCuts->SetRequireITSRefit(kTRUE);//
      fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny); // Reason for the structure in phi
      fESDtrackCuts->SetMaxChi2PerClusterITS(36);
      
      fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);//
      fESDtrackCuts->SetMaxFractionSharedTPCClusters(0.4);//    
      
      // Vertex Related
      fESDtrackCuts->SetMaxDCAToVertexZ( fDcaZ );//
      fESDtrackCuts->SetDCAToVertex2D(kFALSE);
      fESDtrackCuts->SetRequireSigmaToVertex(kFALSE);
    } else 
      Printf(" >>>>  User Track Cuts <<<< ");
    
    /*
      fESDtrackCuts->Print();
      
      Printf(" >>> DCAxy in TC [%8.4f:%8.4f]", 
      fESDtrackCuts->GetMinDCAToVertexXY(), fESDtrackCuts->GetMaxDCAToVertexXY());
      Printf(" >>> DCAz in TC  [%8.4f:%8.4f]", 
      fESDtrackCuts->GetMinDCAToVertexZ(), fESDtrackCuts->GetMaxDCAToVertexZ());
      
      Float_t r1,r2;
      fESDtrackCuts->GetPtRange(r1,r2);
      Printf(" >>> Pt in TC  [%10.4f:%10.4f]",r1,r2);
      
      fESDtrackCuts->GetRapRange(r1,r2);
      Printf(" >>> Rap in TC [%10.4f:%10.4f]",r1,r2);
      
      fESDtrackCuts->GetEtaRange(r1,r2);
      Printf(" >>> Eta in TC [%10.4f:%10.4f]",r1,r2);
    */
  }     
  
  CreateEffCont();
  
  fHistCent = new TH1F("fHistCentPid","Centrality", 100, -0.5, 99.5);			 
  fThnList->Add(fHistCent);
  
  fEventCounter = new TH1D("fEventCounter","EventCounter", 20, -0.5, 19.5);
  fThnList->Add(fEventCounter);
  
  PostData(1, fThnList);
  
}

//----------------------------------------------------------------------------------
void AliEbyEPidEfficiencyContamination::CreateEffCont() {
  
  const Char_t *fgkHistName[4] = {"Nch", "Npi", "Nka", "Npr"};
  
  const Char_t *fgkHistCharge[2] = {"Minus", "Plus"};
  
  const Int_t cbin = 81;
  Double_t CentBins[cbin+1];
  for(Int_t ic = 0; ic <= cbin; ic++) CentBins[ic] = ic - 0.5;
  const Int_t ebin = 8;
  Double_t EtaBins[ebin+1];
  for( Int_t ie = 0; ie <= ebin; ie++) EtaBins[ie] = ie - 0.5;

  const Int_t ptBins = 19;
  Double_t pidPtBins[ptBins+1] = { 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.55 };

  //---For----ThnSparse-----
  const Int_t dim = 305; //1 centrality bin + (19 pt bins) * 2 * 8 (eta bins)
  Int_t bin[dim];
  bin[0] = 81;
  
  for(Int_t ibin = 1; ibin < dim; ibin++) bin[ibin] = 6;
  
  Double_t min[dim];
  for(Int_t jbin = 0; jbin < dim; jbin++) min[jbin] =  -0.5;
  
  Double_t max[dim];
  max[0] = 80.5;
  for(Int_t jbin = 1; jbin < dim; jbin++) max[jbin] = 5.5;
  
  
  fPtBinNplusNminusCh = new THnSparseI("fPtBinNplusNminusCh","cent-nplus-nminus", dim, bin, min, max);
  fThnList->Add(fPtBinNplusNminusCh);
  
  if(fIsMC){

    fPtBinNplusNminusChTruth = new THnSparseI("fPtBinNplusNminusChTruth","cent-nplus-nminus", dim, bin, min, max);
    fThnList->Add(fPtBinNplusNminusChTruth);

    //Book the Eff.cont. histos--------
    for(Int_t k = 0; k < 2; k++) {//loop over charge type; k = 0 -ve charge, k = 1, +ve ch.

      Int_t pidtype = fPidType;
      
      fHitCentRec[k] = new TH2F(Form("fHist%s%sRecCent", fgkHistCharge[k], fgkHistName[pidtype]),Form("Rec%s", fgkHistName[pidtype]), 100, -0.5, 99.5, 1900, -0.5, 1899.5);
      fHitCentRec[k]->Sumw2();
      fThnList->Add(fHitCentRec[k]);
      
      fHitCentGen[k] = new TH2F(Form("fHist%s%sGenCent", fgkHistCharge[k], fgkHistName[pidtype]),Form("Gen%s", fgkHistName[pidtype]), 100, -0.5, 99.5, 1900, -0.5, 1899.5);
      fHitCentGen[k]->Sumw2();
      fThnList->Add(fHitCentGen[k]);
      
      //Gen----
      fCentPtEtaPhiThnGen[k]      = new TH3F( Form("fCentPtEtaThnGen%s%s", fgkHistCharge[k], fgkHistName[pidtype]), "Gen:cent-pt-eta", cbin, CentBins, ptBins, pidPtBins, ebin, EtaBins );
      fCentPtEtaPhiThnGen[k]->GetXaxis()->SetTitle("Centrality");
      fCentPtEtaPhiThnGen[k]->GetYaxis()->SetTitle("p_{T}");
      fCentPtEtaPhiThnGen[k]->GetZaxis()->SetTitle("#eta");
      fCentPtEtaPhiThnGen[k]->Sumw2();
      fThnList->Add( fCentPtEtaPhiThnGen[k] );
      
      //Rec----
      fCentPtEtaPhiThnRec[k]      = new TH3F( Form("fCentPtEtaThnRec%s%s", fgkHistCharge[k], fgkHistName[pidtype]), "Rec:cent-pt-eta", cbin, CentBins, ptBins, pidPtBins, ebin, EtaBins );
      fCentPtEtaPhiThnRec[k]->GetXaxis()->SetTitle("Centrality");
      fCentPtEtaPhiThnRec[k]->GetYaxis()->SetTitle("p_{T}");
      fCentPtEtaPhiThnRec[k]->GetZaxis()->SetTitle("#eta");
      fCentPtEtaPhiThnRec[k]->Sumw2();
      fThnList->Add( fCentPtEtaPhiThnRec[k] );
      
      //Rec. Prim---
      fCentPtEtaPhiThnRecPrim[k]  = new TH3F( Form("fCentPtEtaThnRecPrim%s%s", fgkHistCharge[k], fgkHistName[pidtype]), "RecPrim:cent-pt-eta", cbin, CentBins, ptBins, pidPtBins, ebin, EtaBins );
      fCentPtEtaPhiThnRecPrim[k]->GetXaxis()->SetTitle("Centrality");
      fCentPtEtaPhiThnRecPrim[k]->GetYaxis()->SetTitle("p_{T}");
      fCentPtEtaPhiThnRecPrim[k]->GetZaxis()->SetTitle("#eta");
      fCentPtEtaPhiThnRecPrim[k]->Sumw2();
      fThnList->Add( fCentPtEtaPhiThnRecPrim[k] );
      
      //Rec. Sec-------
      fCentPtEtaPhiThnSec[k]      = new TH3F( Form("fCentPtEtaThnSec%s%s", fgkHistCharge[k], fgkHistName[pidtype]), "RecSec:cent-pt-eta", cbin, CentBins, ptBins, pidPtBins, ebin, EtaBins );
      fCentPtEtaPhiThnSec[k]->GetXaxis()->SetTitle("Centrality");
      fCentPtEtaPhiThnSec[k]->GetYaxis()->SetTitle("p_{T}");
      fCentPtEtaPhiThnSec[k]->GetZaxis()->SetTitle("#eta");
      fCentPtEtaPhiThnSec[k]->Sumw2();
      fThnList->Add( fCentPtEtaPhiThnSec[k] );
      
      //Rec. Material---
      fCentPtEtaPhiThnMat[k]      = new TH3F( Form("fCentPtEtaThnMat%s%s", fgkHistCharge[k], fgkHistName[pidtype]), "RecMat:cent-pt-eta", cbin, CentBins, ptBins, pidPtBins, ebin, EtaBins );
      fCentPtEtaPhiThnMat[k]->GetXaxis()->SetTitle("Centrality");
      fCentPtEtaPhiThnMat[k]->GetYaxis()->SetTitle("p_{T}");
      fCentPtEtaPhiThnMat[k]->GetZaxis()->SetTitle("#eta");
      fCentPtEtaPhiThnMat[k]->Sumw2();
      fThnList->Add( fCentPtEtaPhiThnMat[k] );
      
      //Rec. Mis-ID---
      fCentPtEtaPhiThnMisId[k]    = new TH3F( Form("fCentPtEtaThnMisId%s%s", fgkHistCharge[k], fgkHistName[pidtype]), "RecMisID:cent-pt-eta", cbin, CentBins, ptBins, pidPtBins, ebin, EtaBins );
      fCentPtEtaPhiThnMisId[k]->GetXaxis()->SetTitle("Centrality");
      fCentPtEtaPhiThnMisId[k]->GetYaxis()->SetTitle("p_{T}");
      fCentPtEtaPhiThnMisId[k]->GetZaxis()->SetTitle("#eta");
      fCentPtEtaPhiThnMisId[k]->Sumw2();
      fThnList->Add( fCentPtEtaPhiThnMisId[k] );
        
    }//k ---

  }//isMC---
   
  
  //cout <<" Hisotgrams booked " << endl;
  
}

  
//----------------------------------------------------------------------------------
void AliEbyEPidEfficiencyContamination::LocalPost() {
    PostData(1, fThnList);
    
  }
  
//----------------------------------------------------------------------------------
void AliEbyEPidEfficiencyContamination::UserExec( Option_t * ){


  if( fPidType < 1 || fPidType > 3 ){
    AliError("PID type not supported");
    return;
  }
  
  const Int_t dim = 19*2; //number of pT bins
  const Int_t kEta = 8; //Number of Eta bins---
  Int_t ptCh[dim][kEta];
  Int_t ptChMC[dim][kEta];
  
  for(Int_t idx = 0; idx < dim; idx++){
    for( Int_t it = 0; it < kEta; it++){
      ptCh[idx][it] = 0.;
      ptChMC[idx][it] = 0;
    }
  }
  
  fEventCounter->Fill(1);
  if(!fInputHandler)
    fInputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  
  if(!fInputHandler){
    AliError("No InputHandler");
    return;
  }
  
  fVevent = dynamic_cast<AliVEvent*>(fInputHandler->GetEvent());
  if (!fVevent) {
    printf("ERROR: fVEvent not available\n");
    LocalPost();
    return;
  }

  //---------- Initiate MC
  if (fIsMC){
    fMCEvent = NULL;
    fEventCounter->Fill(8);
    if(fIsAOD) {
      fArrayMC = NULL;
      fArrayMC = dynamic_cast<TClonesArray*>(fVevent->FindListObject(AliAODMCParticle::StdBranchName()));
      if (!fArrayMC)
	AliFatal("No array of MC particles found !!!");
    }
    else{

      fMCEventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      if(!fMCEventHandler){
	AliError("No MC Event Handler available");
	LocalPost(); return;
      }
      fMCEvent = fMCEventHandler->MCEvent();
      if (!fMCEvent) {
	Printf("ERROR: Could not retrieve MC event");
	LocalPost(); return; 
      }
      
      fMCStack = fMCEvent->Stack();
      if (!fMCStack) {
	Printf("ERROR: Could not retrieve MC stack");
	LocalPost(); return; 
      }
    }
  }//if---fIsMC-----
  
  //Pile-up cut for Run2------
  if(fRun == "LHC15o"){
    if(!fEventCuts->AcceptEvent(fVevent)) {
      LocalPost();
      return;
    }
  }
  
  const AliVVertex *vertex = fVevent->GetPrimaryVertex();
  if(!vertex) { LocalPost(); return; }

  Bool_t vtest = kFALSE;
  Double32_t fCov[6];
  vertex->GetCovarianceMatrix(fCov);
  if(vertex->GetNContributors() > 0) {
    if(fCov[5] != 0) {
      vtest = kTRUE;
    }
  }
  if(!vtest) { LocalPost(); return; }
  
  if(TMath::Abs(vertex->GetX()) > fVxMax) { LocalPost(); return; }
  if(TMath::Abs(vertex->GetY()) > fVyMax) { LocalPost(); return; }
  if(TMath::Abs(vertex->GetZ()) > fVzMax) { LocalPost(); return; }

  //-----------Centrality task---------------------------------
  if( fRun == "LHC10h" || fRun == "LHC11h" ){
    AliCentrality *centrality = fVevent->GetCentrality();
    if(!centrality) return;
    if (centrality->GetQuality() != 0) { LocalPost(); return; }
    
    fCentrality = centrality->GetCentralityPercentile(fCentralityEstimator.Data());
  }
  else if (fRun == "LHC15o" ){
    AliMultSelection *fMultSelection= (AliMultSelection *) fVevent->FindListObject("MultSelection");
    if(!fMultSelection ){
      cout << "AliMultSelection object not found!" << endl;
      return;
    }
    else fCentrality = fMultSelection->GetMultiplicityPercentile(fCentralityEstimator.Data(), false);

  }
  else if( fRun == "LHC10hAMPT"){
    //----------centrality range from imp. par. 
    AliGenEventHeader* genHeader = fMCEvent->GenEventHeader();
    if(!genHeader) return;
    
    Double_t impactPar = ((AliGenHijingEventHeader*) genHeader)->ImpactParameter();
    
    if( impactPar >= 0. && impactPar < 3.51 ) fCentrality = 1.; //0-5
    if( impactPar >= 3.51 && impactPar < 4.96 ) fCentrality = 6.; //5-10
    if( impactPar >= 4.96 && impactPar < 6.08 ) fCentrality = 11.;//10-15
    if( impactPar >= 6.08 && impactPar < 7.01 ) fCentrality = 16.;//15-20
    if( impactPar >= 7.01 && impactPar < 7.84 ) fCentrality = 21.;//20-25
    if( impactPar >= 7.84 && impactPar < 8.59 ) fCentrality = 26.;//25-30
    if( impactPar >= 8.59 && impactPar < 9.27 ) fCentrality = 31.;//30-35
    if( impactPar >= 9.27 && impactPar < 9.92 ) fCentrality = 36.;//35-40
    if( impactPar >= 9.92 && impactPar < 10.5 ) fCentrality = 41.;//40-45
    if( impactPar >= 10.5 && impactPar < 11.1 ) fCentrality = 46.;//45-50
    if( impactPar >= 11.1 && impactPar < 11.6 ) fCentrality = 51.;//50-55
    if( impactPar >= 11.6 && impactPar < 12.1 ) fCentrality = 56.;//55-60
    if( impactPar >= 12.1 && impactPar < 12.6 ) fCentrality = 61.;//60-65
    if( impactPar >= 12.6 && impactPar < 13.1 ) fCentrality = 66.;//65-70
    if( impactPar >= 13.1 && impactPar < 13.6 ) fCentrality = 71.;//70-75
    if( impactPar >= 13.6 && impactPar < 14.0 ) fCentrality = 76.;//75-80
  
  }
  else{
    cout <<"Wrong run period for centrality" << endl;
    LocalPost();
    return;
  }
  
  if( fCentrality < 0 || fCentrality >=80 ) return;

  fHistCent->Fill(fCentrality);

  fEventCounter->Fill(2);
  
  
  //----------
  
  fEventCounter->Fill(3);
  fNTracks  = fVevent->GetNumberOfTracks();
 
  Double_t nRec[2] = {0., 0.};
  Double_t nGen[2] = {0., 0.};
  
  //track loop
  
  for (Int_t idxTrack = 0; idxTrack < fNTracks; ++idxTrack) {
    
    AliVTrack *track = static_cast<AliVTrack*>(fVevent->GetTrack(idxTrack));
    
    if(!AcceptTrackL(track)) continue;
    
    Int_t icharge = track->Charge() < 0 ? 0 : 1;    
    Float_t lPt  = (Float_t)track->Pt();
    Float_t lPz  = track->Pz();
    Float_t lP = 0.;
    if( fTotP) lP = track->GetInnerParam()->GetP();
    Float_t lEta = (Float_t)track->Eta();
    Float_t lPhi = (Float_t)track->Phi();

    //Get the pt (p) bin--
    Int_t iptbin = -1;
    if( fTotP) iptbin = GetPtBin(lP);  //total momentum 
    else iptbin = GetPtBin(lPt);
    if( iptbin < 0 || iptbin > fNptBins-1 ) continue;

    //Get the Eta bin--
    Int_t etabin = -1;
    etabin = GetEtaBin( TMath::Abs(lEta) );
    if( etabin < 0 || etabin > 7 ) continue;
    
    //-----------------Fill EffHistos--------
    Float_t RecContainer[3];
    RecContainer[0] = fCentrality;
    if( fTotP ) RecContainer[1] = lP;
    else RecContainer[1] = lPt;
    RecContainer[2] = etabin;

    Bool_t isPid = kFALSE;
    
    if (fPidType != 0 ) {
      isPid = IsPidPassed(track);
      if(isPid){

	if(fIsMC){
	  fCentPtEtaPhiThnRec[icharge]->Fill( RecContainer[0], RecContainer[1], RecContainer[2] );
	}

	nRec[icharge] += 1.;
	
	if(icharge == 1){
	  ptCh[iptbin][etabin] += 1;
	}
	if(icharge == 0){
	  ptCh[iptbin+fNptBins][etabin] += 1;
	}
      }
      
    }//PID----
    
    //---------------------------------
    // - MC Loop for Physical Primary -
    //--------------------------------- 
    if (fIsMC) {
      
      Int_t label  = TMath::Abs(track->GetLabel()); 
      
      Bool_t isPhysicalPrimary        = 0;
      Bool_t isSecondaryFromWeakDecay = 0;
      Bool_t isSecondaryFromMaterial  = 0;
      AliVParticle* particle = NULL;
      
      if (track->InheritsFrom("AliESDtrack")) {
	particle = static_cast<AliVParticle*>(fMCEvent->GetTrack(label));
	if (!particle) return;
	isPhysicalPrimary        = fMCStack->IsPhysicalPrimary(label);
	isSecondaryFromWeakDecay = fMCStack->IsSecondaryFromWeakDecay(label);
	isSecondaryFromMaterial  = fMCStack->IsSecondaryFromMaterial(label);
      }
      else{
	particle                 =  static_cast<AliVParticle*>(fArrayMC->At(label));
	isPhysicalPrimary        =  (static_cast<AliAODMCParticle*>(particle))->IsPhysicalPrimary();
	isSecondaryFromWeakDecay =  (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromWeakDecay();
	isSecondaryFromMaterial  =  (static_cast<AliAODMCParticle*>(particle))->IsSecondaryFromMaterial();
      }
      
      
      Float_t fpTRec   = particle->Pt(); //pT
      Float_t fpZRec   = particle->Pz();
      Float_t fPRec    = particle->P(); //total p
      Float_t fEtaRec  = particle->Eta();
      Float_t fPhiRec  = particle->Phi();
      Int_t pdg        = TMath::Abs(particle->PdgCode());

      //Get the Eta (rap) bin--
      Int_t etabinRecMC = -1;
      etabinRecMC =  GetEtaBin( TMath::Abs(fEtaRec) );
      if( etabinRecMC < 0 || etabinRecMC > 7 ) continue;
      
      Double_t RecMCContainer[3];
      RecMCContainer[0] = fCentrality;
      if( fTotP ) RecMCContainer[1] = fPRec;
      else RecMCContainer[1] = fpTRec;
      RecMCContainer[2] = etabinRecMC;
      
      //PID------
      if(isPid){
	if(isPhysicalPrimary){
	  if ( pdg == fMcPid ) {// For PID
	    fCentPtEtaPhiThnRecPrim[icharge]->Fill( RecMCContainer[0], RecMCContainer[1], RecMCContainer[2] );
	  }
	  else {
	    fCentPtEtaPhiThnMisId[icharge]->Fill( RecMCContainer[0], RecMCContainer[1], RecMCContainer[2] );
	  }
	}//------------isPhysics Primary -------------
	else if (isSecondaryFromWeakDecay) {
	  fCentPtEtaPhiThnSec[icharge]->Fill( RecMCContainer[0], RecMCContainer[1], RecMCContainer[2] );

	}//------isSecondaryFromWeakDecay--------------
	else if (isSecondaryFromMaterial) {
	  fCentPtEtaPhiThnMat[icharge]->Fill( RecMCContainer[0], RecMCContainer[1], RecMCContainer[2] );

	}//------isSecondaryFromMaterial--------------------
	
      }// if(IsPid)---PID finish
      
    }//IsMC
    
  }//rec track loop --

  if( fIsMC ){  
    fHitCentRec[0]->Fill(fCentrality,nRec[0]);
    fHitCentRec[1]->Fill(fCentrality,nRec[1]);
  }
  
  const Int_t thndim = dim*kEta;
  Double_t ptContainer[thndim+1];
  
  ptContainer[0] = (Double_t)fCentrality;
  
  for(Int_t ipt = 0; ipt < dim; ipt++){
    for(Int_t jeta = 0; jeta < kEta; jeta++){
      Int_t k = (ipt*kEta) + jeta;
      ptContainer[k+1] = ptCh[ipt][jeta];
      if( ptCh[ipt][jeta] > 4 ) fEventCounter->Fill( 17 );
    }
  }
  
  fPtBinNplusNminusCh->Fill(ptContainer);
  
  /*
  Int_t recplus = 0, recminus = 0;
  Int_t kk = fNptBins*kEta;
  for(Int_t ik = 1; ik < kk; ik++) {
    recplus  += ptContainer[ik+1];
    recminus += ptContainer[kk+ik+1];
  }
  */
  
    //cout << "In centrality " << fCentrality << " total -ve Rec =" << nRec[0] << "  THnMinus = " << recminus << endl;
    //cout << "In centrality " << fCentrality << " total +ve Rec =" << nRec[1] << "  THnPlus= " << recplus << endl;
  
  fEventCounter->Fill(7);
  //---- - -- - - - - -   -  -- - - - ---- - - - ---
  if (fIsMC) {
    fEventCounter->Fill(8);
    if (fIsAOD) {
      for (Int_t idxMC = 0; idxMC < fArrayMC->GetEntries(); idxMC++) {
	AliAODMCParticle *particle = static_cast<AliAODMCParticle*>(fArrayMC->At(idxMC));
	if (!particle) 
	  continue;
	
	if (!particle->IsPhysicalPrimary()) continue;
	if (!AcceptTrackLMC((AliVParticle*)particle)) continue;
	Int_t icharge = (particle->PdgCode() < 0) ? 0 : 1;
	
	Float_t fpTGen   = particle->Pt();
	Float_t fpzGen   = particle->Pz(); 
	Float_t fpGen    = particle->P();
	Float_t fEtaGen  = particle->Eta();
	Float_t fPhiGen  = particle->Phi();
	
	Int_t pdg  = TMath::Abs(particle->PdgCode());       
	if (pdg != fMcPid) continue;
	
	Int_t iptbinMC = -1;
	if( fTotP ) iptbinMC = GetPtBin(fpGen); // Total p bin
	else iptbinMC = GetPtBin(fpTGen); //pT bin	
	if( iptbinMC < 0 || iptbinMC > fNptBins-1 ) continue;
	
	Int_t etabinMC = -1;
	etabinMC = GetEtaBin( TMath::Abs(fEtaGen) );
	if( etabinMC < 0 || etabinMC > 7 ) continue;
	
	Double_t GenContainer[3];
	GenContainer[0] = fCentrality;
	if( fTotP ) GenContainer[1] = fpGen;
	else  GenContainer[1] = fpTGen;
	GenContainer[2] = etabinMC;

	fCentPtEtaPhiThnGen[icharge]->Fill( GenContainer[0], GenContainer[1], GenContainer[2] );
      
	nGen[icharge] += 1.;
	
	if(icharge == 1){
	  ptChMC[iptbinMC][etabinMC] += 1;
	}
	if(icharge == 0){
	  ptChMC[iptbinMC+fNptBins][etabinMC] += 1;
	}
	
      } //AOD track
      fEventCounter->Fill(9);
    }//if(AOD)----
    else{
      fEventCounter->Fill(10);
      for (Int_t idxMC = 0; idxMC < fMCStack->GetNprimary(); ++idxMC) {
	AliVParticle* particle = fMCEvent->GetTrack(idxMC);
	if (!particle) 
	  continue;
	if(!fMCStack->IsPhysicalPrimary(idxMC))  continue;
	
	if (!AcceptTrackLMC(particle)) continue;
	Int_t icharge = (particle->PdgCode() < 0) ? 0 : 1;
	
	Float_t fpTGen   = particle->Pt();
	Float_t fpzGen   = particle->Pz();
	Float_t fpGen    = particle->P();
	Float_t fEtaGen  = particle->Eta();
	Float_t fPhiGen  = particle->Phi();
	Int_t pdg = TMath::Abs(particle->PdgCode());
	if (pdg != fMcPid) continue;
	
	//Pt (p) bin-------------------
	Int_t iptbinMC = -1;
	if( fTotP ) iptbinMC = GetPtBin(fpGen); // Total p bin
	else iptbinMC = GetPtBin(fpTGen); //pT bin
	if( iptbinMC < 0 || iptbinMC > fNptBins-1 ) continue;

	//Eta bin---
	Int_t etabinMC = -1;
	etabinMC = GetEtaBin( TMath::Abs(fEtaGen) );
	if( etabinMC < 0 || etabinMC > 7 ) continue;
	
	Double_t GenContainer[3];
	GenContainer[0] = fCentrality;
	if( fTotP ) GenContainer[1] = fpGen;
	else  GenContainer[1] = fpTGen;
	GenContainer[2] = etabinMC;
	
	fCentPtEtaPhiThnGen[icharge]->Fill( GenContainer[0], GenContainer[1], GenContainer[2] );
	
	nGen[icharge] += 1.;

	//-------------
	if(icharge == 1){
	  ptChMC[iptbinMC][etabinMC] += 1;
	}
	if(icharge == 0){
	  ptChMC[iptbinMC+fNptBins][etabinMC] += 1;
	}	
      }
      
      
      fEventCounter->Fill(11);
    }//else (ESD)
    
    fHitCentGen[0]->Fill(fCentrality, nGen[0]);
    fHitCentGen[1]->Fill(fCentrality, nGen[1]);
      
    Double_t ptContainerMC[thndim+1];
    ptContainerMC[0] = (Double_t) fCentrality;
    
    for(Int_t ipt = 0; ipt < dim; ipt++){
      for(Int_t jeta = 0; jeta < kEta; jeta++){
	Int_t k = (ipt*kEta) + jeta;
	ptContainerMC[k+1] = ptChMC[ipt][jeta];
	if( ptChMC[ipt][jeta] > 4 ) fEventCounter->Fill( 18 );
      }
    }
    
    fPtBinNplusNminusChTruth->Fill(ptContainerMC);

    /*
    //For debuging--
    Int_t kk = fNptBins*kEta;
    for(Int_t ik = 1; ik < kk; ik++) {
    recplusMC  += ptContainerMC[ik+1];
    recminusMC += ptContainerMC[kk+ik+1];
    }
    
    cout << " Gen positve " << nGen[0] <<" and -ve particle " << nGen[1] << endl;
    cout <<" Thn pos = " << recplusMC << "  and Thn neg = " << recminusMC << endl;
    */
    
  }
  
  fEventCounter->Fill(12);
  PostData(1, fThnList);
}

//___________________________________________________________
Bool_t AliEbyEPidEfficiencyContamination::AcceptTrackL(AliVTrack *track) const {
  
  if (!track) return kFALSE; 
  if (track->Charge() == 0) return kFALSE; 
    
  if (fIsAOD) {  // AOD
    AliAODTrack * trackAOD = dynamic_cast<AliAODTrack*>(track);
    if (!trackAOD) {
      AliError("Pointer to dynamic_cast<AliAODTrack*>(track) = ZERO");
      return kFALSE; 
    }
    if (!trackAOD->TestFilterBit(fAODtrackCutBit))
      return kFALSE;
  } else {      // ESDs
    if(!fESDtrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))  return kFALSE;
  }

  Double_t pt = track->Pt();
  Double_t pz = track->Pz();
  
  if( fTotP ){
    if( !track->GetInnerParam() ) return kFALSE;
    Double_t ptotTPC = track->GetInnerParam()->GetP();//total momentum;
    if( ptotTPC < fPtMin || ptotTPC > fPtMax )  return kFALSE; //cut on momentum (to compare with Identity method result)
  }
  else{
    if( pt < fPtMin || pt > fPtMax )  return kFALSE; //pT cut
  }
  
  if( fIsRapCut ){
    Double_t rap = GetRapidity( pt, pz );
    if( TMath::Abs(rap) > 0.5 ) return kFALSE;//rapidity cut
    if( TMath::Abs(track->Eta()) > fEtaMax ) return kFALSE;
  }
  else{
    if( TMath::Abs(track->Eta()) > fEtaMax ) return kFALSE; 
  }
  
  return kTRUE;
  
}


//___________________________________________________________
Bool_t AliEbyEPidEfficiencyContamination::AcceptTrackLMC(AliVParticle *particle) const {

  if(!particle) return kFALSE;
  if(particle->Charge() == 0.0) return kFALSE;

  Double_t ptotMC = particle->P();
  Double_t ptMC = particle->Pt();
  Double_t pzMC = particle->Pz();
  
  if( fTotP ){
    if ( ptotMC < fPtMin || ptotMC > fPtMax )  return kFALSE; //cut on momentum (to compare with Identity method result)
  }
  else{
    if ( ptMC < fPtMin || ptMC > fPtMax) return kFALSE;
  }
  
  //rapidity cut--------
  if( fIsRapCut ){
    Double_t rapMC = GetRapidity( ptMC, pzMC );
    if( TMath::Abs(rapMC) > 0.5 ) return kFALSE; //rapidity cut
    if( TMath::Abs(particle->Eta()) > fEtaMax ) return kFALSE;
  }
  else{
    if( TMath::Abs(particle->Eta()) > fEtaMax ) return kFALSE;
  }
  
  return kTRUE;
  
}
//---------------------------------------------------------------

//___________________________________________________________
Double_t AliEbyEPidEfficiencyContamination::GetRapidity(Float_t pt, Float_t pz) const{
  
  Double_t partMass = AliPID::ParticleMass(fParticleSpecies);
  Double_t en = TMath::Sqrt( pt*pt + pz*pz + partMass*partMass );
  Double_t rap = -999.;
  if( en != TMath::Abs(pz) ){
    rap = 0.5*TMath::Log( (en + pz)/(en - pz) );
  }
  else rap = -999.;
  
  return rap;
  
}
//---------------------------------------------------------------

//___________________________________________________________
Int_t AliEbyEPidEfficiencyContamination::GetPtBin(Double_t pt){
  
  Int_t bin = -1;
  
  Double_t pidPtBins[20] = { 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.55 };
  for(Int_t iBin = 0; iBin < fNptBins; iBin++){
    
    if( iBin == fNptBins-1){
      if( pt >= pidPtBins[iBin] && pt <= pidPtBins[iBin+1]){
	bin = iBin;
	break;
      }
    }
    else{
      if( pt >= pidPtBins[iBin] && pt < pidPtBins[iBin+1]){
	bin = iBin;
	break;
	
      }
    }
  }//for
  
  return bin;
  
}

//------------------------------------------------------------------------

//________________________________________________________________________
Int_t AliEbyEPidEfficiencyContamination::GetEtaBin(Float_t eta){

  Int_t etabin = -1;

  Float_t EtaRange[9] = {  0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
  
  for(Int_t iBin = 0; iBin < 8; iBin++){
    
    if( iBin == 7){
      if( eta >= EtaRange[iBin] && eta <= EtaRange[iBin+1]){
	etabin = iBin;
	break;
      }
    }
    else{
      if( eta >= EtaRange[iBin] && eta < EtaRange[iBin+1]){
	etabin = iBin;
	break;
	
      }
    }
  }//for
  
  return etabin;
  
}
//___________________________________________________________
void AliEbyEPidEfficiencyContamination::Terminate( Option_t * ){
  //
  //
  Printf(" ------------------------------------------\n"
	 "             Teminating the task           \n"
	 "--------------------------------------------\n");
  
  //
}

void AliEbyEPidEfficiencyContamination::SetPidType(Int_t i) {
  fPidType = i; fMcPid = -9999;
  if (fPidType == 1) { fParticleSpecies = AliPID::kPion; fMcPid = 211; }
  if (fPidType == 2) { fParticleSpecies = AliPID::kKaon; fMcPid = 321;}
  if (fPidType == 3) { fParticleSpecies = AliPID::kProton; fMcPid = 2212;}
 }

Double_t AliEbyEPidEfficiencyContamination::TOFBetaCalc(AliVTrack *track) const{
  //TOF beta calculation
  Double_t tofTime=track->GetTOFsignal();
  
  Double_t c=TMath::C()*1.E-9;// m/ns
  Float_t startTime = fPIDResponse->GetTOFResponse().GetStartTime(((AliVTrack*)track)->P());//in ps
  Double_t length= fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kElectron)*1E-3*c;
  tofTime -= startTime;      // subtract startTime to the signal
  Double_t tof= tofTime*1E-3; // ns, average T0 fill subtracted, no info from T0detector 	 
  tof=tof*c;
  return length/tof;
}

//////////////////////////////////////////////////////////////////////////////////////////////////

Double_t AliEbyEPidEfficiencyContamination::GetMass(AliPID::EParticleType id) const{
  //return Mass according to AliHelperParticleSpecies_t. If undefined return -999.
  Double_t mass=-999.;
  if (id == AliPID::kProton) { mass=9.38271999999999995e-01; }
  if (id == AliPID::kKaon)   { mass=4.93676999999999977e-01; }
  if (id == AliPID::kPion)   { mass=1.39570000000000000e-01; }
  return mass;
}


//___________________________________________________________
Bool_t AliEbyEPidEfficiencyContamination::IsPidPassed(AliVTrack * track) {
  //return  kTRUE; - PID strategy is from Jochen, Hans and Deepika

  Bool_t isAcceptedITS    = kFALSE;
  Bool_t isAcceptedTPC    = kFALSE;
  Bool_t isAcceptedTOF    = kFALSE;
  Bool_t isAccepted       = kFALSE;
  Bool_t hasPIDTOF = kFALSE;
  
  Double_t *pid = new Double_t[3];
  pid[0] = 10.;
  pid[1] = 10.;
  pid[2] = 10.;
  
  Double_t pt = track->Pt();
  Double_t ptot = track->GetInnerParam()->GetP();//total momentum
  
  //---------------------------| el, mu,  pi,  k,    p   | Pt cut offs from spectra
  //ITS--------------
  Double_t ptLowITS[5]       = { 0., 0., 0.2,  0.2,  0.3  };
  Double_t ptHighITS[5]      = { 0., 0., 0.6,  0.6,  1.1  };
  //TPC---------------
  Double_t ptLowTPC[5]       = { 0., 0., 0.2,  0.325, 0.3  };
  Double_t ptHighTPC[5]      = { 0., 0., 2.0,  2.0,   2.0  };
  //TOF----
  Double_t ptLowTOF[5]       = { 0., 0., 0.2,  0.625,  1.1  };
  Double_t ptHighTOF[5]      = { 0., 0., 2.0,  2.0,    2.0  };
  //TPCTOF----------
  Double_t ptLowTPCTOF[5]    = { 0., 0., 0.65, 0.69,   0.8  };
  Double_t ptHighTPCTOF[5]   = { 0., 0., 2.0,  2.00,   2.0  };
  

  //--------------------------------ITS PID--------------------------
  if(fPIDResponse->CheckPIDStatus((AliPIDResponse::EDetector)AliPIDResponse::kITS, track) == AliPIDResponse::kDetPidOk){
    pid[0] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)AliPIDResponse::kITS, track, fParticleSpecies);
    
    if(TMath::Abs(pid[0]) < fNSigmaMaxITS) isAcceptedITS = kTRUE;

    Double_t nSigmaPion = TMath::Abs(fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kPion));
    Double_t nSigmaKaon = TMath::Abs(fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kKaon));
    Double_t nSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kElectron));
    
    if (TMath::Abs(pid[0]) > nSigmaEl) isAcceptedITS = kFALSE;
    
    if( fPidStrategy == 2){
      if( TMath::Abs(pid[0]) > nSigmaPion) isAcceptedITS = kFALSE;
      if( TMath::Abs(pid[0]) > nSigmaKaon) isAcceptedITS = kFALSE;
    }
    
  }
  
  //--------------------------TPC PID----------------
  if (fPIDResponse->CheckPIDStatus((AliPIDResponse::EDetector)AliPIDResponse::kTPC, track) == AliPIDResponse::kDetPidOk) {

    pid[1] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)AliPIDResponse::kTPC, track, fParticleSpecies);
    
    if(fParticleSpecies == 3){////kaon------
      if( track->Pt() > 0.525 && track->Pt() < 0.6){
	if (TMath::Abs(pid[1]) < 1.)  // nsigma < 1
	  isAcceptedTPC = kTRUE;
      }
      else if( track->Pt() >= 0.6 && track->Pt() < 0.8){
	if(pid[1] > -0.5 && pid[1] < 1.)  // asymmetry cut on nsigma
	  isAcceptedTPC = kTRUE;
      }
      else 
	if(TMath::Abs(pid[1]) < fNSigmaMaxTPC ) isAcceptedTPC = kTRUE;
      
    }//kaon------
    else{
      
      if (TMath::Abs(pid[1]) < fNSigmaMaxTPC ) isAcceptedTPC = kTRUE;
    }
    
    Double_t nSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron));
    Double_t nSigmaPion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kPion));
    Double_t nSigmaKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kKaon));
    
    if (TMath::Abs(pid[1]) > nSigmaEl) isAcceptedTPC = kFALSE;

    if( fPidStrategy == 2){
      if( pt > 0.9 ){
	if (TMath::Abs(pid[1]) > nSigmaPion) isAcceptedTPC = kFALSE;
	if (TMath::Abs(pid[1]) > nSigmaKaon) isAcceptedTPC = kFALSE;
      }
    }

    
    
  }
  
  
  //-------------------------------TOF--------------------------
  if ( fPIDResponse->CheckPIDStatus((AliPIDResponse::EDetector)AliPIDResponse::kTOF, track) == AliPIDResponse::kDetPidOk) {
    pid[2] = fPIDResponse->NumberOfSigmas((AliPIDResponse::EDetector)AliPIDResponse::kTOF, track, fParticleSpecies);
    hasPIDTOF = kTRUE;
    if (TMath::Abs(pid[2]) < fNSigmaMaxTOF) isAcceptedTOF = kTRUE;
    
    Double_t nSigmaEl = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kElectron));
    Double_t nSigmaPion = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kPion));
    Double_t nSigmaKaon = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kKaon));
    
    if (TMath::Abs(pid[2]) > nSigmaEl) isAcceptedTOF = kFALSE;

    if( fPidStrategy == 2){
      if( pt > 0.9 ){
	if (TMath::Abs(pid[2]) > nSigmaPion) isAcceptedTOF = kFALSE;
	if (TMath::Abs(pid[2]) > nSigmaKaon) isAcceptedTOF = kFALSE;
      }
    }
    
  }
  
  
  if (fIsMC && isAcceptedTOF) {
    Int_t tofLabel[3];  
    if (track->InheritsFrom("AliESDtrack")) {
      (dynamic_cast<AliESDtrack*>(track))->GetTOFLabel(tofLabel);
    } else if (track->InheritsFrom("AliAODTrack")) {
      (dynamic_cast<AliAODTrack*>(track))->GetTOFLabel(tofLabel);
    }
    
    Bool_t hasMatchTOF = kTRUE;
    if (TMath::Abs(track->GetLabel()) != TMath::Abs(tofLabel[0]) || tofLabel[1] > 0) {
      hasMatchTOF = kFALSE;
    }
    
    if(fIsAOD) {
      //------
    } else {
      TParticle *matchedTrack = fMCStack->Particle(TMath::Abs(tofLabel[0]));
      if (TMath::Abs(matchedTrack->GetFirstMother()) == TMath::Abs(track->GetLabel())) 
	hasMatchTOF = kTRUE;
    }
    isAcceptedTOF = hasMatchTOF;
  }

  //--------Combined--PID------------
  
  if (fParticleSpecies == 2){//for Pion: TPC+TOF---
    
    if(fPidStrategy == 0){
      isAccepted = isAcceptedTPC && isAcceptedTOF;
    }
    else if( fPidStrategy == 1){  
      Double_t nsigCombined = TMath::Sqrt( pid[1]*pid[1] +  pid[2]*pid[2] );
      if( nsigCombined < fNSigmaMaxTOF ) isAccepted = kTRUE;
    }
    else if( fPidStrategy == 2){
      isAccepted = isAcceptedTOF;
    } 
  }
  else if( fParticleSpecies == 3){//for kaon: TPC and/or TOF

    if ( pt > ptLowTOF[fParticleSpecies] && pt < ptHighTOF[fParticleSpecies] ) isAccepted = isAcceptedTOF;
    else isAccepted =  isAcceptedTPC; 
  }
  
  else if( fParticleSpecies == 4){//for proton
    
    if(fPidStrategy == 0){
      //ITS+TPC and TPC+TOF
      if( pt >= ptLowITS[fParticleSpecies] && pt <= ptHighITS[fParticleSpecies] ) isAccepted = isAcceptedITS && isAcceptedTPC;
      else isAccepted = isAcceptedTPC && isAcceptedTOF;
    }
    else if( fPidStrategy == 1){
      
      if( pt >= ptLowITS[fParticleSpecies] && pt <= ptHighITS[fParticleSpecies] ) {
	Double_t nsigCompITSTPC = TMath::Sqrt( pid[1]*pid[1] +  pid[0]*pid[0] );
	if( nsigCompITSTPC < fNSigmaMaxTPC ) isAccepted = kTRUE;
      }
      else {
	Double_t nsigCompTPCTOF = TMath::Sqrt( pid[1]*pid[1] +  pid[2]*pid[2] );
	if( nsigCompTPCTOF < fNSigmaMaxTPC ) isAccepted = kTRUE;
      }
    }
    else if( fPidStrategy == 2){     
      isAccepted = isAcceptedTOF && isAcceptedTPC;
    }
    else if( fPidStrategy == 3){
      if( pt >= 0.3 && pt <= 0.575 ) isAccepted = isAcceptedITS && isAcceptedTPC;
      else if( pt >= 0.825 && pt < 2.0 ) isAccepted = isAcceptedTPC && isAcceptedTOF;
      else isAccepted =  isAcceptedTPC;
      
    }
    else if( fPidStrategy == 4){
      if( pt >= 0.825 && pt < 2.0 ) isAccepted = isAcceptedTPC && isAcceptedTOF;
      else isAccepted =  isAcceptedTPC;
    }
    
  }//for proton
  
  
  delete [] pid;
  return isAccepted;
  
  
}//IsPidPassed


//----------------------------------------------------------------------------------

