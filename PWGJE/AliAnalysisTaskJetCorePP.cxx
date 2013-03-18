
// ******************************************
// This task computes several jet observables like 
// the fraction of energy in inner and outer coronnas,
// jet-track correlations,triggered jet shapes and 
// correlation strength distribution of particles inside jets.    
// Author: lcunquei@cern.ch
// *******************************************


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


#include "TChain.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"

#include "AliLog.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliCentrality.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliInputEventHandler.h"
#include "AliAODJetEventBackground.h"
#include "AliAODMCParticle.h"
//#include "AliAnalysisTaskFastEmbedding.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODJet.h"

#include "AliAnalysisTaskJetCorePP.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskJetCorePP)

//Filip Krizek 1st March 2013

//---------------------------------------------------------------------
AliAnalysisTaskJetCorePP::AliAnalysisTaskJetCorePP() :
AliAnalysisTaskSE(),
fESD(0x0),
fAODIn(0x0),
fAODOut(0x0),
fAODExtension(0x0),
fJetBranchName(""),
fListJets(0x0),
fNonStdFile(""),
fSystem(0), //pp=0  pPb=1
fJetParamR(0.4),
fOfflineTrgMask(AliVEvent::kAny),
fMinContribVtx(1),
fVtxZMin(-10.0),
fVtxZMax(10.0),
fFilterMask(0),
fCentMin(0.0),
fCentMax(100.0),
fJetEtaMin(-0.5),
fJetEtaMax(0.5),
fTriggerEtaCut(0.9),
fTrackEtaCut(0.9),
fTrackLowPtCut(0.15),
fOutputList(0x0),
fHistEvtSelection(0x0),
fh2Ntriggers(0x0),
fHJetSpec(0x0),
fHJetDensity(0x0),
fHJetDensityA4(0x0),
fhJetPhi(0x0),
fhTriggerPhi(0x0),
fhJetEta(0x0),
fhTriggerEta(0x0),
fhVertexZ(0x0),
fhVertexZAccept(0x0),
fhContribVtx(0x0),
fhContribVtxAccept(0x0),
fhDphiTriggerJet(0x0),
fhDphiTriggerJetAccept(0x0),
fhCentrality(0x0),
fhCentralityAccept(0x0),
fHJetPtRaw(0x0),
fHLeadingJetPtRaw(0x0), 
fHDphiVsJetPtAll(0x0), 
fHRhoFastJetVsRhoCone(0x0),
fkAcceptance(2.0*TMath::Pi()*1.8),
fConeArea(TMath::Pi()*0.4*0.4)
{
   // default Constructor
}

//---------------------------------------------------------------------

AliAnalysisTaskJetCorePP::AliAnalysisTaskJetCorePP(const char *name) :
AliAnalysisTaskSE(name),
fESD(0x0),
fAODIn(0x0),
fAODOut(0x0),
fAODExtension(0x0),
fJetBranchName(""),
fListJets(0x0),
fNonStdFile(""),
fSystem(0),  //pp=0   pPb=1
fJetParamR(0.4),
fOfflineTrgMask(AliVEvent::kAny),
fMinContribVtx(1),
fVtxZMin(-10.0),
fVtxZMax(10.0),
fFilterMask(0),
fCentMin(0.0),
fCentMax(100.0),
fJetEtaMin(-0.5),
fJetEtaMax(0.5),
fTriggerEtaCut(0.9),
fTrackEtaCut(0.9),
fTrackLowPtCut(0.15),
fOutputList(0x0),
fHistEvtSelection(0x0),
fh2Ntriggers(0x0),
fHJetSpec(0x0),
fHJetDensity(0x0),
fHJetDensityA4(0x0),
fhJetPhi(0x0),
fhTriggerPhi(0x0),
fhJetEta(0x0),
fhTriggerEta(0x0),
fhVertexZ(0x0),
fhVertexZAccept(0x0),
fhContribVtx(0x0),
fhContribVtxAccept(0x0),
fhDphiTriggerJet(0x0),
fhDphiTriggerJetAccept(0x0),
fhCentrality(0x0),
fhCentralityAccept(0x0),
fHJetPtRaw(0x0),
fHLeadingJetPtRaw(0x0), 
fHDphiVsJetPtAll(0x0), 
fHRhoFastJetVsRhoCone(0x0),
fkAcceptance(2.0*TMath::Pi()*1.8),
fConeArea(TMath::Pi()*0.4*0.4)
{
// Constructor

   DefineOutput(1, TList::Class());
}

//--------------------------------------------------------------
AliAnalysisTaskJetCorePP::AliAnalysisTaskJetCorePP(const AliAnalysisTaskJetCorePP& a):
AliAnalysisTaskSE(a.GetName()),
fESD(a.fESD),
fAODIn(a.fAODIn),
fAODOut(a.fAODOut),
fAODExtension(a.fAODExtension),
fJetBranchName(a.fJetBranchName),
fListJets(a.fListJets),
fNonStdFile(a.fNonStdFile),
fSystem(a.fSystem),  
fJetParamR(a.fJetParamR),
fOfflineTrgMask(a.fOfflineTrgMask),
fMinContribVtx(a.fMinContribVtx),
fVtxZMin(a.fVtxZMin),
fVtxZMax(a.fVtxZMax),
fFilterMask(a.fFilterMask),
fCentMin(a.fCentMin),
fCentMax(a.fCentMax),
fJetEtaMin(a.fJetEtaMin),
fJetEtaMax(a.fJetEtaMax),
fTriggerEtaCut(a.fTriggerEtaCut),
fTrackEtaCut(a.fTrackEtaCut),
fTrackLowPtCut(a.fTrackLowPtCut),
fOutputList(a.fOutputList),
fHistEvtSelection(a.fHistEvtSelection),
fh2Ntriggers(a.fh2Ntriggers),
fHJetSpec(a.fHJetSpec),
fHJetDensity(a.fHJetDensity),
fHJetDensityA4(a.fHJetDensityA4),
fhJetPhi(a.fhJetPhi),
fhTriggerPhi(a.fhTriggerPhi),
fhJetEta(a.fhJetEta),
fhTriggerEta(a.fhTriggerEta),
fhVertexZ(a.fhVertexZ),
fhVertexZAccept(a.fhVertexZAccept),
fhContribVtx(a.fhContribVtx),
fhContribVtxAccept(a.fhContribVtxAccept),
fhDphiTriggerJet(a.fhDphiTriggerJet),
fhDphiTriggerJetAccept(a.fhDphiTriggerJetAccept),
fhCentrality(a.fhCentrality),
fhCentralityAccept(a.fhCentralityAccept),
fHJetPtRaw(a.fHJetPtRaw),
fHLeadingJetPtRaw(a.fHLeadingJetPtRaw),
fHDphiVsJetPtAll(a.fHDphiVsJetPtAll),
fHRhoFastJetVsRhoCone(a.fHRhoFastJetVsRhoCone),
fkAcceptance(a.fkAcceptance),
fConeArea(a.fConeArea)
{
   //Copy Constructor
}
//--------------------------------------------------------------

AliAnalysisTaskJetCorePP& AliAnalysisTaskJetCorePP::operator = (const AliAnalysisTaskJetCorePP& a){
  // assignment operator
  this->~AliAnalysisTaskJetCorePP();
  new(this) AliAnalysisTaskJetCorePP(a);
  return *this;
}
//--------------------------------------------------------------

AliAnalysisTaskJetCorePP::~AliAnalysisTaskJetCorePP()
{
   //Destructor 
   delete fListJets;
   delete fOutputList; // ????
}

//--------------------------------------------------------------

void AliAnalysisTaskJetCorePP::Init()
{
   // check for jet branches
   if(!strlen(fJetBranchName.Data())){
      AliError("Jet branch name not set.");
   }

}

//--------------------------------------------------------------

void AliAnalysisTaskJetCorePP::UserCreateOutputObjects()
{


  // Create histograms
   // Called once
   fListJets = new TList();
 
   OpenFile(1);
   if(!fOutputList) fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, -0.5, 5.5);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"event selection (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(6,"multiplicity (rejected)");
   
   fOutputList->Add(fHistEvtSelection);

   Int_t nBinsCentrality = (fSystem==0) ? 1 : 10; // pp=1 else 10
    
   fh2Ntriggers = new TH2F("fh2Ntriggers","# of triggers",
                             nBinsCentrality,0.0,100.0,50,0.0,50.0);
   fOutputList->Add(fh2Ntriggers);

   //Centrality, A, pT - rho*A, pTtrigg, pT, rho*A
   const Int_t dimSpec   = 6;
   const Int_t nBinsSpec[dimSpec]     = {nBinsCentrality, 100,  140,  50, 100, 100};
   const Double_t lowBinSpec[dimSpec] = {0.0,             0.0,-80.0, 0.0, 0.0, 0.0};
   const Double_t hiBinSpec[dimSpec]  = {100.0,           1.0,200.0,50.0, 200, 50.0};
   fHJetSpec = new THnSparseF("fHJetSpec",
                   "Recoil jet spectrum [cent,A,pTjet-rho*A,pTtrig,pTjet, rho*A]",
                   dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
   fOutputList->Add(fHJetSpec);  

   //------------------- HISTOS FOR DIAGNOSTIC ----------------------
   //Jet number density histos [Trk Mult, jet density, pT trigger]
   const Int_t    dimJetDens   = 3;
   const Int_t    nBinsJetDens[dimJetDens]  = {100,   100, 10};
   const Double_t lowBinJetDens[dimJetDens] = {0.0,   0.0,  0.0};
   const Double_t hiBinJetDens[dimJetDens]  = {500.0, 5.0, 50.0 };

   fHJetDensity = new THnSparseF("fHJetDensity","Jet dens vs trk mult A>0.07",
                                   dimJetDens,nBinsJetDens,lowBinJetDens,hiBinJetDens);

   fHJetDensityA4 =new THnSparseF("fHJetDensityA4","Jet dens vs trk mult A>0.4",
                                   dimJetDens,nBinsJetDens,lowBinJetDens,hiBinJetDens);

   fOutputList->Add(fHJetDensity);
   fOutputList->Add(fHJetDensityA4);
         

   //inclusive azimuthal and pseudorapidity histograms
   fhJetPhi = new TH2D("fhJetPhi","Azim dist jets vs pTjet",
                        50, 0, 100, 50,-TMath::Pi(),TMath::Pi());
   fhTriggerPhi= new TH2D("fhTriggerPhi","azim dist trig had vs pTtrigg",
                        25, 0, 50, 50,-TMath::Pi(),TMath::Pi());
   fhJetEta = new TH2D("fhJetEta","Eta dist jets vs pTjet",
                        50,0, 100, 40,-0.9,0.9);
   fhTriggerEta = new TH2D("fhTriggerEta","Eta dist trig had vs pTtrigg",
                        25, 0, 50, 40,-0.9,0.9);

   fhVertexZ = new TH1D("fhVertexZ","z vertex",40,-20,20);  
   fhVertexZAccept = new TH1D("fhVertexZAccept","z vertex after cut",40,-20,20);  
   fhContribVtx = new TH1D("fhContribVtx","contrib to vtx",200,0,200);   
   fhContribVtxAccept = new TH1D("fhContribVtxAccept","contrib to vtx after cut",200,0,200);   
   fhDphiTriggerJet = new TH1D("fhDphiTriggerJet","Deltaphi trig-jet",50, -TMath::Pi(),TMath::Pi()); 
   fhDphiTriggerJetAccept = new TH1D("fhDphiTriggerJetAccept","Deltaphi trig-jet after cut",50, -TMath::Pi(),TMath::Pi()); 
   fhCentrality = new TH1D("fhCentrality","Centrality",20,0,100);
   fhCentralityAccept = new TH1D("fhCentralityAccept","Centrality after cut",20,0,100);

   fOutputList->Add(fhJetPhi);
   fOutputList->Add(fhTriggerPhi);
   fOutputList->Add(fhJetEta);
   fOutputList->Add(fhTriggerEta);
   fOutputList->Add(fhVertexZ);    
   fOutputList->Add(fhVertexZAccept);    
   fOutputList->Add(fhContribVtx); 
   fOutputList->Add(fhContribVtxAccept); 
   fOutputList->Add(fhDphiTriggerJet);
   fOutputList->Add(fhDphiTriggerJetAccept);
   fOutputList->Add(fhCentrality); 
   fOutputList->Add(fhCentralityAccept);

   // raw spectra of INCLUSIVE jets  
   //Centrality, pTjet, pTjet - rho*A,  A
   const Int_t dimRaw   = 4;
   const Int_t nBinsRaw[dimRaw]     = {nBinsCentrality,  50,    75, 100};
   const Double_t lowBinRaw[dimRaw] = {0.0,             0.0, -50.0, 0.0};
   const Double_t hiBinRaw[dimRaw]  = {100.0,           100, 100.0, 1.0};
   fHJetPtRaw = new THnSparseF("fHJetPtRaw",
                                "Incl. jet spectrum [cent,pTjet,pTjet-rho*A,A]",
                                dimRaw,nBinsRaw,lowBinRaw,hiBinRaw);
   fOutputList->Add(fHJetPtRaw);  

   // raw spectra of LEADING jets  
   //Centrality, pTjet, pTjet - rho*A,  A
   fHLeadingJetPtRaw = new THnSparseF("fHLeadingJetPtRaw",
                                "Leading jet spectrum [cent,pTjet,pTjet-rho*A,A]",
                                dimRaw,nBinsRaw,lowBinRaw,hiBinRaw);
   fOutputList->Add(fHLeadingJetPtRaw);  

   // Dphi versus pT jet 
   //Centrality, Dphi=phiTrig-phiJet, pTjet, pTtrigg 
   const Int_t dimDp   = 4;
   const Int_t nBinsDp[dimDp]     = {nBinsCentrality,  50,     50,    50};
   const Double_t lowBinDp[dimDp] = {0.0,       -TMath::Pi(),   0.0,   0.0};
   const Double_t hiBinDp[dimDp]  = {100.0,      TMath::Pi(), 100.0, 100.0};
   fHDphiVsJetPtAll = new THnSparseF("fHDphiVsJetPtAll",
                                "Dphi vs jet pT [cent,Dphi,pTjet,pTtrigg]",
                                dimDp,nBinsDp,lowBinDp,hiBinDp);
   fOutputList->Add(fHDphiVsJetPtAll);  

   // Rho Fast Jet vs Rho Cone 
   //Centrality, Rho Fast Jet, Rho Perp Cone,   pTjet 
   const Int_t dimRo   = 4;
   const Int_t nBinsRo[dimRo]     = {nBinsCentrality, 100,   100,  50};
   const Double_t lowBinRo[dimRo] = {0.0,             0.0,   0.0,   0.0};
   const Double_t hiBinRo[dimRo]  = {100.0,         100.0, 100.0, 100.0};
   fHRhoFastJetVsRhoCone = new THnSparseF("fHRhoFastJetVsRhoCone",
                      "Rho FastJet vs rho PerpCone [cent,rhoFastJet,rhoCone,pTjet]",
                           dimRo,nBinsRo,lowBinRo,hiBinRo);
   fOutputList->Add(fHRhoFastJetVsRhoCone);  




   // =========== Switch on Sumw2 for all histos ===========
   for(Int_t i=0; i<fOutputList->GetEntries(); i++){
      TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
      if(h1){
         h1->Sumw2();
         continue;
      }
      THnSparse *hn = dynamic_cast<THnSparse*>(fOutputList->At(i));
      if(hn){
         hn->Sumw2();
      }	  
   }
   TH1::AddDirectory(oldStatus);

   PostData(1, fOutputList);
}

//--------------------------------------------------------------------

void AliAnalysisTaskJetCorePP::UserExec(Option_t *)
{

   //Event loop

   if(TMath::Abs((Float_t) fJetParamR)<0.00001){
      AliError("Cone radius is set to zero.");
      return;
   }
   if(!strlen(fJetBranchName.Data())){
      AliError("Jet branch name not set.");
      return;
   }

   fESD = dynamic_cast<AliESDEvent*>(InputEvent());
   if(!fESD){
      AliError("ESD not available");
      fAODIn = dynamic_cast<AliAODEvent*>(InputEvent());
   } 
 
   fAODOut = dynamic_cast<AliAODEvent*>(AODEvent());
   AliAODEvent* aod = NULL;
   // take all other information from the aod we take the tracks from
   if(!aod){
      if(!fESD) aod = fAODIn;
      else      aod = fAODOut;
   }

   if(fNonStdFile.Length()!=0){
      // case that we have an AOD extension we can fetch the jets from the extended output
      AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
      fAODExtension = aodH ? aodH->GetExtension(fNonStdFile.Data()) : 0;
      if(!fAODExtension){
         if(fDebug>1) Printf("AODExtension found for %s",fNonStdFile.Data());
      } 
   }
    
   // ----------------- event selection --------------------------
   fHistEvtSelection->Fill(1); // number of events before event selection

   // physics selection
   AliInputEventHandler* inputHandler = (AliInputEventHandler*)
           ((AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler());

   if(!(inputHandler->IsEventSelected() & fOfflineTrgMask)){
      if(fDebug) Printf(" Trigger Selection: event REJECTED ... ");
      fHistEvtSelection->Fill(2);
      PostData(1, fOutputList);
      return;
   }
  
   //check AOD pointer
   if(!aod){
      if(fDebug) Printf("%s:%d No AOD",(char*)__FILE__,__LINE__);
      fHistEvtSelection->Fill(3);
      PostData(1, fOutputList);
      return;
   }
   
   // vertex selection
   AliAODVertex* primVtx = aod->GetPrimaryVertex();

   if(!primVtx){
      if(fDebug) Printf("%s:%d No primVtx",(char*)__FILE__,__LINE__);
      fHistEvtSelection->Fill(3);
      PostData(1, fOutputList);
      return;
   }

   Int_t nTracksPrim = primVtx->GetNContributors();
   Float_t vtxz = primVtx->GetZ();
   //Input events
   fhContribVtx->Fill(nTracksPrim);
   if( nTracksPrim > 0 ) fhVertexZ->Fill(vtxz);

   if((nTracksPrim < fMinContribVtx) ||
      (vtxz < fVtxZMin) ||
      (vtxz > fVtxZMax)){
      if(fDebug) Printf("%s:%d primary vertex z = %f: event REJECTED...",
                         (char*)__FILE__,__LINE__,vtxz);
      fHistEvtSelection->Fill(3);
      PostData(1, fOutputList);
      return;
   }else{
      //Accepted events
      fhContribVtxAccept->Fill(nTracksPrim);
      fhVertexZAccept->Fill(vtxz);
   }
   //FK// No event class selection imposed
   // event class selection (from jet helper task)
   //Int_t eventClass = AliAnalysisHelperJetTasks::EventClass();
   //if(fDebug) Printf("Event class %d", eventClass);
   //if(eventClass < fEvtClassMin || eventClass > fEvtClassMax){
   //   fHistEvtSelection->Fill(4);
   //   PostData(1, fOutputList);
   //   return;
   //}

   // centrality selection
   AliCentrality *cent = 0x0;
   Double_t centValue  = 0.0; 
   if(fSystem){  //fSystem=0 for pp,   fSystem=1 for pPb
      if(fESD){
         cent = fESD->GetCentrality();
         if(cent) centValue = cent->GetCentralityPercentile("V0M");
      }else{
         centValue = aod->GetHeader()->GetCentrality();
      }   
      if(fDebug) printf("centrality: %f\n", centValue);
      //Input events
      fhCentrality->Fill(centValue); 

      if(centValue < fCentMin || centValue > fCentMax){
         fHistEvtSelection->Fill(4);
         PostData(1, fOutputList);
         return;
      }else{
         //Accepted events
         fhCentralityAccept->Fill( centValue );
      }
   }
 
   if(fDebug) std::cout<<" ACCEPTED EVENT "<<endl;
  
   fHistEvtSelection->Fill(0); // accepted events 
   fConeArea = TMath::Pi()*fJetParamR*fJetParamR;
   // ------------------- end event selection --------------------

   // fetch jets
   TClonesArray *aodJets = 0x0;
   
   if(fAODOut && !aodJets){
      aodJets = dynamic_cast<TClonesArray*>(fAODOut->FindListObject(fJetBranchName.Data()));
   }
   if(fAODExtension && !aodJets){ 
      aodJets = dynamic_cast<TClonesArray*>(fAODExtension->GetAOD()->FindListObject(fJetBranchName.Data()));
   } 
   if(fAODIn && !aodJets){
      aodJets = dynamic_cast<TClonesArray*>(fAODIn->FindListObject(fJetBranchName.Data())); 
   }

   // ------------- Hadron trigger --------------
   TList particleList; //list of tracks
   Int_t indexTrigg = GetListOfTracks(&particleList); //index of trigger hadron in Particle list
   
   if(indexTrigg<0) return; // no trigger track found above 150 MeV/c 

   AliVParticle *triggerHadron = (AliVParticle*) particleList.At(indexTrigg);     
   if(!triggerHadron){  
      PostData(1, fOutputList);
      return;
   }


   fh2Ntriggers->Fill(centValue,triggerHadron->Pt()); //trigger pT

      //Trigger Diagnostics---------------------------------
   fhTriggerPhi->Fill(triggerHadron->Pt(),RelativePhi(triggerHadron->Phi(),0.0)); //phi -pi,pi
   fhTriggerEta->Fill(triggerHadron->Pt(),triggerHadron->Eta());

   //--------- Fill list of jets -------------- 
   fListJets->Clear();
   if(aodJets){
      if(fDebug) Printf("########## %s: %d jets",fJetBranchName.Data(),
                                      aodJets->GetEntriesFast());
      for(Int_t iJet = 0; iJet < aodJets->GetEntriesFast(); iJet++) {
         AliAODJet *jet = dynamic_cast<AliAODJet*>((*aodJets)[iJet]);
         if (jet) fListJets->Add(jet);
      }
   }
   
   Double_t etaJet  = 0.0;
   Double_t pTJet   = 0.0;
   Double_t areaJet = 0.0;
   Double_t phiJet  = 0.0;
   Int_t injet4     = 0;
   Int_t injet      = 0; 
   Int_t indexLeadingJet     = -1;
   Double_t pTLeadingJet     = -10.0; 
   Double_t pTcorrLeadingJet = -10.0; 
   Double_t areaLeadingJet   = -10.0;
   Double_t rhoFastJet       = 0.0; 
   //---------- jet loop ---------
   for(Int_t ij=0; ij<fListJets->GetEntries(); ij++){
      AliAODJet* jet = (AliAODJet*)(fListJets->At(ij));
      if(!jet) continue;
      etaJet  = jet->Eta();
      phiJet  = jet->Phi();
      pTJet   = jet->Pt();
      if(pTJet==0) continue; 
     
      if((etaJet<fJetEtaMin) || (etaJet>fJetEtaMax)) continue;
      areaJet = jet->EffectiveAreaCharged();

      Double_t fastJetbgpT = jet->ChargedBgEnergy();  
      if(areaJet>0) rhoFastJet = fastJetbgpT/areaJet;
      else          rhoFastJet = 0.0;

      //Jet Diagnostics---------------------------------
      fhJetPhi->Fill(pTJet,  RelativePhi(phiJet,0.0)); //phi -pi,pi
      fhJetEta->Fill(pTJet,  etaJet);
      if(areaJet >= 0.07) injet++; 
      if(areaJet >= 0.4)  injet4++;
      //--------------------------------------------------

      Double_t dphi = RelativePhi(triggerHadron->Phi(), phiJet); 
   
      fhDphiTriggerJet->Fill(dphi); //Input
      //Background w.r.t. jet axis
      Double_t pTBckWrtJet = 
         GetBackgroundInPerpCone(fJetParamR, phiJet, etaJet, &particleList);

      Double_t ratioOfAreas = areaJet/fConeArea;
      Double_t rhoA = ratioOfAreas*pTBckWrtJet; //bg activity in a cone of similar size
      Double_t ptcorr = pTJet - rhoA; //Pt Jet UE subtr 

      //search for leading jet
      if(pTJet > pTLeadingJet){
         indexLeadingJet  = ij; 
         pTLeadingJet     = pTJet; 
         pTcorrLeadingJet = ptcorr; 
         areaLeadingJet   = areaJet; 
      } 
 
      // raw spectra of INCLUSIVE jets  
      //Centrality, pTjet, pTjet - rho*A,  A
      Double_t fillraw[] = { centValue,
                             pTJet,
                             ptcorr,
                             areaJet
                           };
      fHJetPtRaw->Fill(fillraw);

      //Dphi versus jet pT   
      //Centrality, Dphi=phiTrig-phiJet, pTjet, pTtrigg 
      Double_t filldp[] = { centValue,
                            dphi,
                            pTJet,
                            triggerHadron->Pt()
                          };
      fHDphiVsJetPtAll->Fill(filldp);

      // Rho Fast Jet vs Rho Cone 
      //Centrality, Rho Fast Jet, Rho Perp Cone, pTjet 
       Double_t fillro[] = { centValue,
                             rhoFastJet,
                             pTBckWrtJet/fConeArea,
                             pTJet
                           };
      fHRhoFastJetVsRhoCone->Fill(fillro);

      // Select back to back trigger - jet pairs
      if(TMath::Abs((Double_t) dphi) < TMath::Pi()-0.6) continue;
      fhDphiTriggerJetAccept->Fill(dphi); //Accepted

 
      //Centrality, A, pT - rho*A, pTtrigg, pT, rho*A
      Double_t fillspec[] = { centValue,
                              areaJet,
                              ptcorr,
                              triggerHadron->Pt(),
                              pTJet,
                              rhoA
                            };
      fHJetSpec->Fill(fillspec);
	   
   }//jet loop
 
   if(indexLeadingJet > -1){ 
     // raw spectra of LEADING jets  
     //Centrality, pTjet, pTjet - rho*A,  A
     Double_t fillleading[] = { centValue,
                                pTLeadingJet,
                                pTcorrLeadingJet,
                                areaLeadingJet
                              };
     fHLeadingJetPtRaw->Fill(fillleading);
   } 

  
   //Fill Jet Density In the Event A>0.07
   if(injet>0){
      Double_t filldens[]={ (Double_t) particleList.GetEntries(),
                            injet/fkAcceptance,
                            triggerHadron->Pt()
                          };
      fHJetDensity->Fill(filldens);
   }

   //Fill Jet Density In the Event A>0.4
   if(injet4>0){ 
      Double_t filldens4[]={ (Double_t) particleList.GetEntries(), 
                             injet4/fkAcceptance,
                             triggerHadron->Pt()
                           };
      fHJetDensityA4->Fill(filldens4);
   }

   PostData(1, fOutputList);
}

//----------------------------------------------------------------------------
void AliAnalysisTaskJetCorePP::Terminate(const Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

   if(fDebug) printf("AliAnalysisTaskJetCorePP DONE\n");
   if(!GetOutputData(1)) return;
}


//----------------------------------------------------------------------------
Int_t  AliAnalysisTaskJetCorePP::GetListOfTracks(TList *list){
   //Fill the list of accepted tracks (passed track cut)
   //return consecutive index of the hardest ch hadron in the list
   Int_t iCount        = 0;
   AliAODEvent *aodevt = NULL;

   if(!fESD) aodevt = fAODIn;
   else      aodevt = fAODOut;   

   if(!aodevt) return -1;

   Int_t    index = -1; //index of the highest particle in the list
   Double_t ptmax = -10;

   for(int it = 0; it < aodevt->GetNumberOfTracks(); it++){
      AliAODTrack *tr = aodevt->GetTrack(it);
      
      //if((fFilterMask > 0) && !(tr->TestFilterBit(fFilterMask))) continue;
      if((fFilterMask > 0) && !(tr->IsHybridGlobalConstrainedGlobal())) continue;
      if(TMath::Abs((Float_t) tr->Eta()) > fTrackEtaCut) continue;
      if(tr->Pt() < fTrackLowPtCut) continue;
      list->Add(tr);
      if(tr->Pt()>ptmax){ 
         ptmax = tr->Pt();	
         index = iCount;
      }
      iCount++;
   }

   if(index>-1){ //check pseudorapidity cut on trigger
      AliAODTrack *trigger = (AliAODTrack*) list->At(index);
      if(trigger && TMath::Abs((Float_t) trigger->Eta())< fTriggerEtaCut){ return index;} 
      return -1;
   }else{
      return -1;
   }
}

//----------------------------------------------------------------------------

Double_t AliAnalysisTaskJetCorePP::GetBackgroundInPerpCone(Float_t jetR, Double_t jetPhi, Double_t jetEta, TList* trkList){
   //calculate sum of track pT in the cone perpendicular in phi to the jet 
   //jetR = cone radius
   // jetPhi, jetEta = direction of the jet 
   Int_t numberOfTrks = trkList->GetEntries();
   Double_t pTsum = 0.0;
   Double_t perpConePhi = jetPhi + TMath::Pi()/2;//perp cone w.r.t. jet in phi
   for(Int_t it=0; it<numberOfTrks; it++){
      AliVParticle *trk = (AliVParticle*) trkList->At(it); 
      Double_t dphi = RelativePhi(perpConePhi,trk->Phi());     
      Double_t deta = trk->Eta()-jetEta;     

      if( (dphi*dphi + deta*deta)< (jetR*jetR)){
         pTsum += trk->Pt();
      } 
   }

   return pTsum;
}

//----------------------------------------------------------------------------

Double_t AliAnalysisTaskJetCorePP::RelativePhi(Double_t mphi,Double_t vphi){
   //Get relative azimuthal angle of two particles -pi to pi
   if      (vphi < -TMath::Pi()) vphi += TMath::TwoPi();
   else if (vphi > TMath::Pi())  vphi -= TMath::TwoPi();

   if      (mphi < -TMath::Pi()) mphi += TMath::TwoPi();
   else if (mphi > TMath::Pi())  mphi -= TMath::TwoPi();

   Double_t dphi = mphi - vphi;
   if      (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
   else if (dphi > TMath::Pi())  dphi -= TMath::TwoPi();

   return dphi;//dphi in [-Pi, Pi]
}


