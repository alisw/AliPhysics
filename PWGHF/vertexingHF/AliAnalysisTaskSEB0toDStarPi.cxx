/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appeuear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//
//
//   Base class for (B0 -> DStar pi -> D0 pi pi -> K pi pi pi) Analysis
//
//
//             Cuts are centralized in AliRDHFCutsB0toDStarPi
//             Like sign background is imlemented in the macro
//
//-----------------------------------------------------------------------
//
//                 Author Lennart van Doremalen
//           Utrecht University - l.v.r.vandoremalen@uu.nl
//
//               Based on the DStartoKpipi macro by:
//
//                         Author A.Grelli 
//              ERC-QGP Utrecht University - a.grelli@uu.nl,
//                         Author Y.Wang
//        University of Heidelberg - yifei@physi.uni-heidelberg.de
//                         Author C.Ivan 
//             ERC-QGP Utrecht University - c.ivan@uu.nl,
//
//-----------------------------------------------------------------------

#include <TSystem.h>
#include <TParticle.h>
#include <TH1I.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliRDHFCutsB0toDStarPi.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliVertex.h"
#include "AliVVertex.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliVertexerTracks.h"
#include "AliExternalTrackParam.h"
#include "AliNeutralTrackParam.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliNormalizationCounter.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskSEB0toDStarPi.h"
#include "AliAODInputHandler.h"
#include <vector>
#include <TMatrix.h>
#include <TVector3.h>
#include <TArrayI.h>
#include <bitset>

// #include "TObjectTable.h" 

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEB0toDStarPi);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSEB0toDStarPi::AliAnalysisTaskSEB0toDStarPi():  
  AliAnalysisTaskSE(),
  fEvents(0),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputD0Pion(0),
  fOutputD0Kaon(0),
  fOutputDStarPion(0),
  fOutputB0Pion(0),
  fOutputD0(0),
  fOutputDStar(0),
  fOutputB0(0),
  fOutputD0_D0Pt(0),
  fOutputD0_DStarPt(0),
  fOutputDStar_DStarPt(0),
  fOutputB0MC(0),
  fCuts(0),
  fQuickSignalAnalysis(0),
  fGetCutInfo(0),
  fCEvents(0),     
  fCounter(0),
  fD0PionTracks(0),
  fD0KaonTracks(0),
  fDStarPionTracks(0),
  fB0PionTracks(0),
  fD0Tracks(0),
  fDStarTracks(0),
  fB0Tracks(0),
  fnPtBins(0),
  fnPtBinsD0forD0ptbin(0),
  fnPtBinsD0forDStarptbin(0),
  fnPtBinsDStarforDStarptbin(0),
  fnPtBinLimits(0),
  fnPtBinsD0forD0ptbinLimits(0),
  fnPtBinsD0forDStarptbinLimits(0),
  fnPtBinsDStarforDStarptbinLimits(0),
  fPtBinLimits(0),
  fPtBinLimitsD0forD0ptbin(0),
  fPtBinLimitsD0forDStarptbin(0),
  fPtBinLimitsDStarforDStarptbin(0),
  fDaughterHistogramArray(),
  fDaughterHistogramArray2D(),
  fMotherHistogramArray(),
  fMotherHistogramArray2D()
{
  //
  /// Default ctor
  //
}
//___________________________________________________________________________
AliAnalysisTaskSEB0toDStarPi::AliAnalysisTaskSEB0toDStarPi(const Char_t* name, AliRDHFCutsB0toDStarPi* cuts) :
  AliAnalysisTaskSE(name),
  fEvents(0),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputD0Pion(0),
  fOutputD0Kaon(0),
  fOutputDStarPion(0),
  fOutputB0Pion(0),
  fOutputD0(0),
  fOutputDStar(0),
  fOutputB0(0),
  fOutputD0_D0Pt(0),
  fOutputD0_DStarPt(0),
  fOutputDStar_DStarPt(0),
  fOutputB0MC(0),
  fCuts(0),
  fQuickSignalAnalysis(0),
  fGetCutInfo(0),
  fCEvents(0),     
  fCounter(0),
  fD0PionTracks(0),
  fD0KaonTracks(0),
  fDStarPionTracks(0),
  fB0PionTracks(0),
  fD0Tracks(0),
  fDStarTracks(0),
  fB0Tracks(0),
  fnPtBins(0),
  fnPtBinsD0forD0ptbin(0),
  fnPtBinsD0forDStarptbin(0),
  fnPtBinsDStarforDStarptbin(0),
  fnPtBinLimits(0),
  fnPtBinsD0forD0ptbinLimits(0),
  fnPtBinsD0forDStarptbinLimits(0),
  fnPtBinsDStarforDStarptbinLimits(0),
  fPtBinLimits(0),
  fPtBinLimitsD0forD0ptbin(0),
  fPtBinLimitsD0forDStarptbin(0),
  fPtBinLimitsDStarforDStarptbin(0),
  fDaughterHistogramArray(),
  fDaughterHistogramArray2D(),
  fMotherHistogramArray(),
  fMotherHistogramArray2D()
{
  //
  /// Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskSEB0toDStarPi","Calling Constructor");

  // we prepare vectors and arrays that will save the candidates during the reconstruction
  fD0PionTracks = new std::vector<Int_t>;
  fD0KaonTracks = new std::vector<Int_t>;
  fDStarPionTracks = new std::vector<Int_t>;
  fB0PionTracks = new std::vector<Int_t>;

  fD0Tracks = new TClonesArray("AliAODRecoDecayHF2Prong", 5000);
  fDStarTracks = new TClonesArray("AliAODRecoCascadeHF", 5000);
  fB0Tracks = new TClonesArray("AliAODRecoCascadeHF", 5000);

  // we get the cut file
  fCuts=cuts;

  // we get information on the pt bins
  fnPtBins = fCuts->GetNPtBins();
  fnPtBinsD0forD0ptbin = fCuts->GetNPtBinsD0forD0ptbin();
  fnPtBinsD0forDStarptbin = fCuts->GetNPtBinsD0forDStarptbin();
  fnPtBinsDStarforDStarptbin = fCuts->GetNPtBinsDStarforDStarptbin();

  fnPtBinLimits = fnPtBins + 1;
  fnPtBinsD0forD0ptbinLimits = fnPtBinsD0forD0ptbin + 1;
  fnPtBinsD0forDStarptbinLimits = fnPtBinsD0forDStarptbin + 1;
  fnPtBinsDStarforDStarptbinLimits = fnPtBinsDStarforDStarptbin + 1;

  fPtBinLimits = fCuts->GetPtBinLimits();
  fPtBinLimitsD0forD0ptbin = fCuts->GetPtBinLimitsD0forD0ptbin();
  fPtBinLimitsD0forDStarptbin = fCuts->GetPtBinLimitsD0forDStarptbin();
  fPtBinLimitsDStarforDStarptbin = fCuts->GetPtBinLimitsDStarforDStarptbin();

  // we create an array of pointers for the histograms. This method is more CPU efficient than looking up each histogram by name.
  // this option is now set manualy in the header file
  // const Int_t numberOfDaughters = 4;
  // const Int_t numberOfDaughterHistogramSets = 5;
  // const Int_t numberOfDaughterHistograms = 15;
  // const Int_t numberOfDaughterHistograms2D = 6;

  // Int_t maxHistogramSets = 6 + 2*fnPtBins;
  // if(2*fnPtBinsD0forD0ptbin > maxHistogramSets) maxHistogramSets = 2*fnPtBinsD0forD0ptbin;
  // if(2*fnPtBinsD0forDStarptbin > maxHistogramSets) maxHistogramSets = 2*fnPtBinsD0forDStarptbin;
  // if(2*fnPtBinsDStarforDStarptbin > maxHistogramSets) maxHistogramSets = 2*fnPtBinsDStarforDStarptbin;

  // const Int_t numberOfOutputs = 6;
  // const Int_t numberOfMotherHistogramSets = maxHistogramSets;
  // const Int_t numberOfMotherHistograms = 46;
  // const Int_t numberOfMotherHistograms2D = 7;

  // fDaughterHistogramArray = new Int_t*[numberOfDaughters][numberOfDaughterHistogramSets][numberOfDaughterHistograms];
  // fDaughterHistogramArray2D = new Int_t*[numberOfDaughters][numberOfDaughterHistograms2D];
  // fMotherHistogramArray = new Int_t*[numberOfOutputs][numberOfMotherHistogramSets][numberOfMotherHistograms];
  // fMotherHistogramArray2D = new Int_t*[numberOfOutputs][numberOfMotherHistograms2D];


  DefineOutput(1,TList::Class());  //counters
  DefineOutput(2,TList::Class());  //All Entries output
  DefineOutput(3,TList::Class());  //3sigma PID output
  DefineOutput(4,AliRDHFCutsB0toDStarPi::Class());   //My private output
  DefineOutput(5,AliNormalizationCounter::Class());   // normalization
  DefineOutput(6,TList::Class());  // D0 pion output
  DefineOutput(7,TList::Class());  // D0 kaon output
  DefineOutput(8,TList::Class());  // DStar pion output
  DefineOutput(9,TList::Class());  // B0 pion output
  DefineOutput(10,TList::Class());  // D0 output
  DefineOutput(11,TList::Class());  // DStar output
  DefineOutput(12,TList::Class());   // B0 output
  DefineOutput(13,TList::Class());   // B0 output
  DefineOutput(14,TList::Class());   // B0 output
  DefineOutput(15,TList::Class());   // B0 output
  DefineOutput(16,TList::Class());   // B0MC output
}

//___________________________________________________________________________
AliAnalysisTaskSEB0toDStarPi::~AliAnalysisTaskSEB0toDStarPi() {
  //
  /// destructor
  //
  Info("~AliAnalysisTaskSEB0toDStarPi","Calling Destructor");

  delete fOutput;
  delete fOutputD0Pion;
  delete fOutputD0Kaon;
  delete fOutputDStarPion;
  delete fOutputB0Pion;
  delete fOutputD0;
  delete fOutputDStar;
  delete fOutputB0;
  delete fOutputD0_D0Pt;
  delete fOutputD0_DStarPt;
  delete fOutputDStar_DStarPt;
  delete fOutputB0MC;
  delete fCuts;
  delete fCEvents;
  delete fD0PionTracks;
  delete fD0KaonTracks; 
  delete fD0Tracks;
  delete fDStarPionTracks;
  delete fDStarTracks;
  delete fB0PionTracks;
  delete fB0Tracks;
  // delete [] fDaughterHistogramArray;
  // delete [] fDaughterHistogramArray2D;
  // delete [] fMotherHistogramArray;
  // delete [] fMotherHistogramArray2D;
}
//_________________________________________________
void AliAnalysisTaskSEB0toDStarPi::Init(){
  //
  /// Initialization
  //

  if(fDebug > 1) printf("AliAnalysisTaskSEB0toDStarPi::Init() \n");
   AliRDHFCutsB0toDStarPi* copyfCuts=new AliRDHFCutsB0toDStarPi(*(static_cast<AliRDHFCutsB0toDStarPi*>(fCuts)));
  // Post the data
  PostData(4,copyfCuts);

  return;
}

//_________________________________________________
void AliAnalysisTaskSEB0toDStarPi::UserExec(Option_t *){

  //==================================================================================
  //  USER EXECUTION FUNCTION - start
  //==================================================================================
  //
  // This is the main function for the heavy flavour analysis.
  //
  //==================================================================================

  if (!fInputEvent) {
    Error("UserExec","NO EVENT FOUND!");
    return;
  }

  if(fEvents%100==0){
    std::cout << "\r" << "Analysing event number: " << fEvents << std::endl;
  }

  fEvents++;

  //add option to show trigger mask
  // std::cout << "event test: " <<((AliAODInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() << std::endl;
  // bitset<32> testEV(((AliAODInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
  // std::cout << "event test: " << testEV << std::endl;
  // std::cout << "triger mask: " << bitset<32>(fCuts->GetTriggerMask()) << std::endl;
  // ( & fTriggerMask);

  //==================================================================================
  //  EVENT INITIALIZATION - start
  //==================================================================================


  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  fCEvents->Fill(1);

  if(!aodEvent && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
  } 

  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
  fCEvents->Fill(2);

  fCounter->StoreEvent(aodEvent,fCuts,fUseMCInfo);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD
  TString trigclass=aodEvent->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD")||trigclass.Contains("C0SMH-B-NOPF-ALL")) fCEvents->Fill(5);

  if(!fCuts->IsEventSelected(aodEvent)) {
    // std::cout << "Event rejected by code: " << fCuts->GetWhyRejection() << std::endl;
    if(fCuts->GetWhyRejection()==6) {// rejected for Z vertex
      fCEvents->Fill(6);
      return;
    }
  }

  Bool_t isEvSel=fCuts->IsEventSelected(aodEvent);
  if(!isEvSel) return;
  fCEvents->Fill(3);

  // counters for efficiencies
  Int_t icountReco = 0;
  
  //get the magnetic field
  Double_t bz = (Double_t)aodEvent->GetMagneticField(); 

  // AOD primary vertex
  AliAODVertex *primaryVertex = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if(!primaryVertex) return;
  if(primaryVertex->GetNContributors()<1) return;
  fCEvents->Fill(4);

  Int_t nSelectedAna =0;
  Int_t nSelectedProd =0;

  //==================================================================================
  //  EVENT INITIALIZATION - end
  //==================================================================================
  //  B0 MC SIGNAL IDENTIFICATION - start
  //==================================================================================

  // We create an array that contains all the monte carlo particles in the event
  TClonesArray *mcTrackArray = 0; 
  if(fUseMCInfo) mcTrackArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));

  // We create an array to save the MC labels of true signal tracks
  TMatrix * B0toDStarPiLabelMatrix = new TMatrix(0,7);
   
  // We fill the array with all B0->DStarPi tracks
  if(fUseMCInfo) {
    B0toDStarPiSignalTracksInMC(mcTrackArray,aodEvent,B0toDStarPiLabelMatrix, fOutputB0MC);
  }
  //==================================================================================
  //  B0 MC SIGNAL IDENTIFICATION - end
  //==================================================================================
  //  PARTICLE SELECTION LOOP - start
  //==================================================================================
  //
  // Here we select and reconstruct the particles for the B0->D*Pion decay. 
  //
  //==================================================================================


  D0PionSelection(aodEvent,mcTrackArray,B0toDStarPiLabelMatrix);
  D0KaonSelection(aodEvent,mcTrackArray,B0toDStarPiLabelMatrix);
  DStarPionSelection(aodEvent,mcTrackArray,B0toDStarPiLabelMatrix);
  B0PionSelection(aodEvent,mcTrackArray,B0toDStarPiLabelMatrix);

  D0Selection(aodEvent,primaryVertex,bz,mcTrackArray,B0toDStarPiLabelMatrix);
  DStarAndB0Selection(aodEvent,primaryVertex,bz,mcTrackArray,B0toDStarPiLabelMatrix);


  // Memory management: for all entries of mother getvertex and delete. 
  for (int i = 0; i < fD0Tracks->GetEntriesFast(); ++i){ 
    AliAODVertex * vertex = (AliAODVertex*)((AliAODRecoDecayHF2Prong*)(*fD0Tracks)[i])->GetOwnSecondaryVtx();
    delete vertex; vertex = NULL;
  }
  for (int i = 0; i < fDStarTracks->GetEntriesFast(); ++i){ 
    AliAODVertex * vertex = (AliAODVertex*)((AliAODRecoCascadeHF*)(*fDStarTracks)[i])->GetOwnSecondaryVtx();
    delete vertex; vertex = NULL;
  }
  for (int i = 0; i < fB0Tracks->GetEntriesFast(); ++i){ 
    AliAODVertex * vertex = (AliAODVertex*)((AliAODRecoCascadeHF*)(*fB0Tracks)[i])->GetOwnSecondaryVtx();
    delete vertex; vertex = NULL;
  }

  fB0Tracks->Clear("C");
  fDStarTracks->Clear("C");
  fD0Tracks->Clear("C");

  fD0PionTracks->erase(fD0PionTracks->begin(),fD0PionTracks->end());
  fD0KaonTracks->erase(fD0KaonTracks->begin(),fD0KaonTracks->end());
  fB0PionTracks->erase(fB0PionTracks->begin(),fB0PionTracks->end());
  fDStarPionTracks->erase(fDStarPionTracks->begin(),fDStarPionTracks->end());
  
  delete B0toDStarPiLabelMatrix; B0toDStarPiLabelMatrix = NULL;

  //==================================================================================
  //  PARTICLE SELECTION LOOP - end
  //==================================================================================

  fCounter->StoreCandidates(aodEvent,nSelectedProd,kTRUE);  
  fCounter->StoreCandidates(aodEvent,nSelectedAna,kFALSE); 

  AliDebug(2, Form("Found %i Reco particles that are D*!!",icountReco));

  PostData(1,fOutput);
  PostData(5,fCounter);
  PostData(6,fOutputD0Pion);
  PostData(7,fOutputD0Kaon);
  PostData(8,fOutputDStarPion);
  PostData(9,fOutputB0Pion);
  PostData(10,fOutputD0);
  PostData(11,fOutputDStar);
  PostData(12,fOutputB0);
  PostData(13,fOutputD0_D0Pt);
  PostData(14,fOutputD0_DStarPt);
  PostData(15,fOutputDStar_DStarPt);
  PostData(16,fOutputB0MC);


  //==================================================================================
  //  USER EXECUTION FUNCTION - end
  //==================================================================================

}
//________________________________________ terminate ___________________________
void AliAnalysisTaskSEB0toDStarPi::Terminate(Option_t*){    
  /// The Terminate() function is the last function to be called during
  /// a query. It always runs on the client, it can be used to present
  /// the results graphically or save the results to file.
  
  //Info("Terminate","");
  AliAnalysisTaskSE::Terminate();
  
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  
  fCEvents        = dynamic_cast<TH1F*>(fOutput->FindObject("fCEvents"));

  fOutputD0Pion = dynamic_cast<TList*> (GetOutputData(6));
  if (!fOutputD0Pion) {
    printf("ERROR: fOutputD0Pion not available\n");
    return;
  }
  fOutputD0Kaon = dynamic_cast<TList*> (GetOutputData(7));
  if (!fOutputD0Kaon) {
    printf("ERROR: fOutputD0Kaon not available\n");
    return;
  }
  fOutputDStarPion = dynamic_cast<TList*> (GetOutputData(8));
  if (!fOutputDStarPion) {
    printf("ERROR: fOutputDStarPion not available\n");
    return;
  }
  fOutputB0Pion = dynamic_cast<TList*> (GetOutputData(9));
  if (!fOutputB0Pion) {
    printf("ERROR: fOutputB0Pion not available\n");
    return;
  }
  fOutputD0 = dynamic_cast<TList*> (GetOutputData(10));
  if (!fOutputD0) {
    printf("ERROR: fOutputD0 not available\n");
    return;
  }
  fOutputDStar = dynamic_cast<TList*> (GetOutputData(11));
  if (!fOutputDStar) {
    printf("ERROR: fOutputDStar not available\n");
    return;
  }
  fOutputB0 = dynamic_cast<TList*> (GetOutputData(12));
  if (!fOutputB0) {
    printf("ERROR: fOutputB0 not available\n");
    return;
  }
  fOutputD0_D0Pt = dynamic_cast<TList*> (GetOutputData(13));
  if (!fOutputD0_D0Pt) {
    printf("ERROR: fOutputD0_D0Pt not available\n");
    return;
  }
  fOutputD0_DStarPt = dynamic_cast<TList*> (GetOutputData(14));
  if (!fOutputD0_DStarPt) {
    printf("ERROR: fOutputD0_DStarPt not available\n");
    return;
  }
  fOutputDStar_DStarPt = dynamic_cast<TList*> (GetOutputData(15));
  if (!fOutputDStar_DStarPt) {
    printf("ERROR: fOutputDStar_DStarPt not available\n");
    return;
  }  
  fOutputB0MC = dynamic_cast<TList*> (GetOutputData(16));
  if (!fOutputB0MC) {
    printf("ERROR: fOutputB0MC not available\n");
    return;
  }
  return;
}

//___________________________________________________________________________
void AliAnalysisTaskSEB0toDStarPi::UserCreateOutputObjects() { 
 /// output
  Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
  
  //slot #1  
  //OpenFile(1);
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("chist0");

  fOutputD0Pion = new TList();
  fOutputD0Pion->SetOwner();
  fOutputD0Pion->SetName("listD0Pion");

  fOutputD0Kaon = new TList();
  fOutputD0Kaon->SetOwner();
  fOutputD0Kaon->SetName("listD0Kaon");

  fOutputDStarPion = new TList();
  fOutputDStarPion->SetOwner();
  fOutputDStarPion->SetName("listDStarPion");

  fOutputB0Pion = new TList();
  fOutputB0Pion->SetOwner();
  fOutputB0Pion->SetName("listB0Pion");

  fOutputD0 = new TList();
  fOutputD0->SetOwner();
  fOutputD0->SetName("listD0");

  fOutputDStar = new TList();
  fOutputDStar->SetOwner();
  fOutputDStar->SetName("listDStar");

  fOutputB0 = new TList();
  fOutputB0->SetOwner();
  fOutputB0->SetName("listB0");

  fOutputD0_D0Pt = new TList();
  fOutputD0_D0Pt->SetOwner();
  fOutputD0_D0Pt->SetName("listD0_D0Pt");

  fOutputD0_DStarPt = new TList();
  fOutputD0_DStarPt->SetOwner();
  fOutputD0_DStarPt->SetName("listD0_DStarPt");

  fOutputDStar_DStarPt = new TList();
  fOutputDStar_DStarPt->SetOwner();
  fOutputDStar_DStarPt->SetName("listDStar_DStarPt");  

  fOutputB0MC = new TList();
  fOutputB0MC->SetOwner();
  fOutputB0MC->SetName("listB0MC");
    
  // define histograms
  DefineHistograms();
  
  //Counter for Normalization
  fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(5)->GetContainer()->GetName()));
  fCounter->Init();

  PostData(1,fOutput);
  PostData(6,fOutputD0Pion);
  PostData(7,fOutputD0Kaon);
  PostData(8,fOutputDStarPion);
  PostData(9,fOutputB0Pion);
  PostData(10,fOutputD0);
  PostData(11,fOutputDStar);
  PostData(12,fOutputB0);
  PostData(13,fOutputD0_D0Pt);
  PostData(14,fOutputD0_DStarPt);
  PostData(15,fOutputDStar_DStarPt);
  PostData(16,fOutputB0MC);

  return;
}
//___________________________________ histograms _______________________________________
void  AliAnalysisTaskSEB0toDStarPi::DefineHistograms(){
  /// Create histograms


  fCEvents = new TH1F("fCEvents","conter",13,0,13);
  fCEvents->SetStats(kTRUE);
  fCEvents->GetXaxis()->SetTitle("1");
  fCEvents->GetYaxis()->SetTitle("counts");
  fCEvents->GetXaxis()->SetBinLabel(2,"no. of events");
  fCEvents->GetXaxis()->SetBinLabel(3,"good prim vtx and B field");
  fCEvents->GetXaxis()->SetBinLabel(4,"no event selected");
  fCEvents->GetXaxis()->SetBinLabel(5,"no vtx contributors");
  fCEvents->GetXaxis()->SetBinLabel(6,"trigger for PbPb");
  fCEvents->GetXaxis()->SetBinLabel(7,"no z vtx");
  fCEvents->GetXaxis()->SetBinLabel(12,"no. of D0 fail to be rec");
  fOutput->Add(fCEvents);

 //====================================================

  TString name_mc_B0_pt ="mc_B0_pt";
  TH1F* hist_mc_B0_pt = new TH1F(name_mc_B0_pt.Data(),"Pt monte carlo B0 in B0->D*#pi; p_{T} [GeV/c]; Entries",400,0,20);
  hist_mc_B0_pt->Sumw2();
  hist_mc_B0_pt->SetLineColor(6);
  hist_mc_B0_pt->SetMarkerStyle(20);
  hist_mc_B0_pt->SetMarkerSize(0.6);
  hist_mc_B0_pt->SetMarkerColor(6);
  TH1F* histogram_mc_B0_pt = (TH1F*)hist_mc_B0_pt->Clone();
  fOutputB0MC->Add(histogram_mc_B0_pt);

  TString name_mc_B0_pion_pt ="mc_B0_pion_pt";
  TH1F* hist_mc_B0_pion_pt = new TH1F(name_mc_B0_pion_pt.Data(),"Pt monte carlo pion of B0 in B0->D*#pi; p_{T} [GeV/c]; Entries",400,0,20);
  hist_mc_B0_pion_pt->Sumw2();
  hist_mc_B0_pion_pt->SetLineColor(6);
  hist_mc_B0_pion_pt->SetMarkerStyle(20);
  hist_mc_B0_pion_pt->SetMarkerSize(0.6);
  hist_mc_B0_pion_pt->SetMarkerColor(6);
  TH1F* histogram_mc_B0_pion_pt = (TH1F*)hist_mc_B0_pion_pt->Clone();
  fOutputB0MC->Add(histogram_mc_B0_pion_pt);

  TString name_mc_DStar_pt ="mc_DStar_pt";
  TH1F* hist_mc_DStar_pt = new TH1F(name_mc_DStar_pt.Data(),"Pt monte carlo DStar in B0->D*#pi; p_{T} [GeV/c]; Entries",400,0,20);
  hist_mc_DStar_pt->Sumw2();
  hist_mc_DStar_pt->SetLineColor(6);
  hist_mc_DStar_pt->SetMarkerStyle(20);
  hist_mc_DStar_pt->SetMarkerSize(0.6);
  hist_mc_DStar_pt->SetMarkerColor(6);
  TH1F* histogram_mc_DStar_pt = (TH1F*)hist_mc_DStar_pt->Clone();
  fOutputB0MC->Add(histogram_mc_DStar_pt);

  TString name_mc_DStar_pion_pt ="mc_DStar_pion_pt";
  TH1F* hist_mc_DStar_pion_pt = new TH1F(name_mc_DStar_pion_pt.Data(),"Pt monte carlo pion of DStar in B0->D*#pi; p_{T} [GeV/c]; Entries",400,0,20);
  hist_mc_DStar_pion_pt->Sumw2();
  hist_mc_DStar_pion_pt->SetLineColor(6);
  hist_mc_DStar_pion_pt->SetMarkerStyle(20);
  hist_mc_DStar_pion_pt->SetMarkerSize(0.6);
  hist_mc_DStar_pion_pt->SetMarkerColor(6);
  TH1F* histogram_mc_DStar_pion_pt = (TH1F*)hist_mc_DStar_pion_pt->Clone();
  fOutputB0MC->Add(histogram_mc_DStar_pion_pt);

  TString name_mc_D0_pt ="mc_D0_pt";
  TH1F* hist_mc_D0_pt = new TH1F(name_mc_D0_pt.Data(),"Pt monte carlo D0 in B0->D*#pi; p_{T} [GeV/c]; Entries",400,0,20);
  hist_mc_D0_pt->Sumw2();
  hist_mc_D0_pt->SetLineColor(6);
  hist_mc_D0_pt->SetMarkerStyle(20);
  hist_mc_D0_pt->SetMarkerSize(0.6);
  hist_mc_D0_pt->SetMarkerColor(6);
  TH1F* histogram_mc_D0_pt = (TH1F*)hist_mc_D0_pt->Clone();
  fOutputB0MC->Add(histogram_mc_D0_pt);

  TString name_mc_D0_pion_pt ="mc_D0_pion_pt";
  TH1F* hist_mc_D0_pion_pt = new TH1F(name_mc_D0_pion_pt.Data(),"Pt monte carlo pion of D0 in B0->D*#pi; p_{T} [GeV/c]; Entries",400,0,20);
  hist_mc_D0_pion_pt->Sumw2();
  hist_mc_D0_pion_pt->SetLineColor(6);
  hist_mc_D0_pion_pt->SetMarkerStyle(20);
  hist_mc_D0_pion_pt->SetMarkerSize(0.6);
  hist_mc_D0_pion_pt->SetMarkerColor(6);
  TH1F* histogram_mc_D0_pion_pt = (TH1F*)hist_mc_D0_pion_pt->Clone();
  fOutputB0MC->Add(histogram_mc_D0_pion_pt);

  TString name_mc_D0_kaon_pt ="mc_D0_kaon_pt";
  TH1F* hist_mc_D0_kaon_pt = new TH1F(name_mc_D0_kaon_pt.Data(),"Pt monte carlo kaon of D0 in B0->D*#pi; p_{T} [GeV/c]; Entries",400,0,20);
  hist_mc_D0_kaon_pt->Sumw2();
  hist_mc_D0_kaon_pt->SetLineColor(6);
  hist_mc_D0_kaon_pt->SetMarkerStyle(20);
  hist_mc_D0_kaon_pt->SetMarkerSize(0.6);
  hist_mc_D0_kaon_pt->SetMarkerColor(6);
  TH1F* histogram_mc_D0_kaon_pt = (TH1F*)hist_mc_D0_kaon_pt->Clone();
  fOutputB0MC->Add(histogram_mc_D0_kaon_pt);

 //==================================================

  TString name_B0s_in_analysis ="B0s_in_analysis";
  TH1F* hist_B0s_in_analysis = new TH1F(name_B0s_in_analysis.Data(),"Number of B0 to kpipipi in the Analysis; Entries",10,0,10);
  hist_B0s_in_analysis->Sumw2();
  hist_B0s_in_analysis->SetLineColor(6);
  hist_B0s_in_analysis->SetMarkerStyle(20);
  hist_B0s_in_analysis->SetMarkerSize(0.6);
  hist_B0s_in_analysis->SetMarkerColor(6);
  hist_B0s_in_analysis->SetStats(kTRUE);
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(1,"no. of B0s");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(2,"no. of B0s to kpipipi");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(3,"no. with all tracks in event");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(4,"no. ...");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(5,"no. ...");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(6,"no. ...");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(7,"no. ...");
  hist_B0s_in_analysis->GetXaxis()->SetBinLabel(8,"no. ...");
  TH1F* hist_B0s_in_analysis_mc = (TH1F*)hist_B0s_in_analysis->Clone();
  fOutputB0MC->Add(hist_B0s_in_analysis_mc);

  TString name_B0s_per_bin ="B0s_per_bin";
  TH1F* hist_B0s_per_bin = new TH1F(name_B0s_per_bin.Data(),"Number of B0 to kpipipi in the Analysis per bin; Entries",6,0,6); // add pt bin info
  hist_B0s_per_bin->GetXaxis()->SetBinLabel(1,"0-3");
  hist_B0s_per_bin->GetXaxis()->SetBinLabel(2,"3-6");
  hist_B0s_per_bin->GetXaxis()->SetBinLabel(3,"6-10");
  hist_B0s_per_bin->GetXaxis()->SetBinLabel(4,"10-18");
  hist_B0s_per_bin->GetXaxis()->SetBinLabel(5,"18-30");
  hist_B0s_per_bin->GetXaxis()->SetBinLabel(6,"30-inf");
  TH1F* hist_B0s_per_bin_mc = (TH1F*)hist_B0s_per_bin->Clone();
  fOutputB0MC->Add(hist_B0s_per_bin_mc);

 //======================================================================================================================================================

  //we make the histograms for the Pions and Kaon
  for (Int_t i = 0; i < 4; i++){
    
    TString add_name = "";
    TList * listout;
    if(i==0) listout = fOutputD0Pion;
    if(i==1) listout = fOutputD0Kaon;
    if(i==2) listout = fOutputDStarPion;
    if(i==3) listout = fOutputB0Pion;

    for (Int_t j = 0; j < 5; j++){
      if(j==0) add_name = "";
      if(j==1) add_name = "Cut";
      if(j==2) add_name = "MC";
      if(j==3) add_name = "Result";
      if(j==4) add_name = "MCResult";

      TString name_Histogram = "";
      TString discription_Histogram = "";
      Int_t numberOfBins = 0;
      Double_t lowerBound = 0.0;
      Double_t upperBound = 0.0;

      for (Int_t k = 0; k < 11; ++k)
      {
        if(k==0){name_Histogram = "ptTrack"; discription_Histogram = "pt track; p_{T} [GeV/c]; Entries"; numberOfBins = 600; lowerBound = 0; upperBound = 30;}
        if(k==1){name_Histogram = "momentumTrack"; discription_Histogram = "momentum track; p [GeV/c]; Entries"; numberOfBins = 600; lowerBound = 0; upperBound = 30;}
        if(k==2){name_Histogram = "energyOverMomentumTrack"; discription_Histogram = "E/P; E/P; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 10;}
        if(k==3){name_Histogram = "dcaTrack"; discription_Histogram = "dca track; distance [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 1;}
        if(k==4){name_Histogram = "momentumdcaTrack"; discription_Histogram = "momentum at dca track; p_{T} [GeV/c]; Entries"; numberOfBins = 600; lowerBound = 0; upperBound = 30;}
        if(k==5){name_Histogram = "numberOfITS"; discription_Histogram = "Number of ITS clusters track; [#]; Entries"; numberOfBins = 7; lowerBound = -0.5; upperBound = 6.5;}
        if(k==6){name_Histogram = "numberOfTPC"; discription_Histogram = "Number of TPC clusters track; [#]; Entries"; numberOfBins = 601; lowerBound = -0.5; upperBound = 600.5;}
        if(k==7){name_Histogram = "pointsOnITS"; discription_Histogram = "Number of ITS clusters track per layer; [#]; Entries"; numberOfBins = 6; lowerBound = -0.5; upperBound = 5.5;}
        if(k==8){name_Histogram = "nSigmaTPC"; discription_Histogram = "n sigma TPC for track PID; sigma; Entries"; numberOfBins = 500; lowerBound = -5; upperBound = 5;}
        if(k==9){name_Histogram = "nSigmaTOF"; discription_Histogram = "n sigma TOF for track PID; sigma; Entries"; numberOfBins = 500; lowerBound = -5; upperBound = 5;}
        if(k==10){name_Histogram = "nSigmaTPCandTOF"; discription_Histogram = "n sigma TPC and TOF for track PID; a.u.; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 10;}

        name_Histogram += add_name;
        TH1F* histogram = new TH1F(name_Histogram.Data(),discription_Histogram.Data(),numberOfBins,lowerBound,upperBound);
        histogram->Sumw2();
        if(j%2==0) histogram->SetLineColor(6);
        if(j%2==1) histogram->SetLineColor(4);
        histogram->SetMarkerStyle(20);
        histogram->SetMarkerSize(0.6);
        if(j%2==0) histogram->SetMarkerColor(6);
        if(j%2==1) histogram->SetMarkerColor(4);
        TH1F* histogram_Clone = (TH1F*)histogram->Clone();
        listout->Add(histogram_Clone);
        fDaughterHistogramArray[i][j][k] = (Int_t*)histogram_Clone;
      }

      TString name_ptdcaTrack ="pt_vs_dcaTrack";
      name_ptdcaTrack += add_name;
      TH2F* hist_ptdcaTrack = new TH2F(name_ptdcaTrack.Data(),"Pt vs dca track; p_{T} [GeV/c]; distance [cm]",100,0,30,200,0,10);
      hist_ptdcaTrack->Sumw2();
      hist_ptdcaTrack->SetLineColor(6);
      hist_ptdcaTrack->SetMarkerStyle(20);
      hist_ptdcaTrack->SetMarkerSize(0.6);
      hist_ptdcaTrack->SetMarkerColor(6);
      TH2F* histogram_ptdcaTrack = (TH2F*)hist_ptdcaTrack->Clone();
      listout->Add(histogram_ptdcaTrack);
      fDaughterHistogramArray[i][j][11] = (Int_t*)histogram_ptdcaTrack;

      TString numberofparticlesperevent="numberofparticlesperevent";
      numberofparticlesperevent += add_name;
      TH1F* hist_numberofparticlesperevent = new TH1F(numberofparticlesperevent.Data(),"Number of particles per event; number of particles in one event; Entries",100,0,100);
      hist_numberofparticlesperevent->Sumw2();
      hist_numberofparticlesperevent->SetLineColor(6);
      hist_numberofparticlesperevent->SetMarkerStyle(20);
      hist_numberofparticlesperevent->SetMarkerSize(0.6);
      hist_numberofparticlesperevent->SetMarkerColor(6);
      TH1F* histogram_numberofparticlesperevent = (TH1F*)hist_numberofparticlesperevent->Clone();
      listout->Add(histogram_numberofparticlesperevent);
      fDaughterHistogramArray[i][j][12] = (Int_t*)histogram_numberofparticlesperevent;
    }

    TH1F * effectOfCuts = new TH1F("effectOfCuts","Removal counter",18,0,18);
    effectOfCuts->SetStats(kTRUE);
    effectOfCuts->GetXaxis()->SetTitle("Cut number");
    effectOfCuts->GetYaxis()->SetTitle("Particles cut");
    effectOfCuts->GetXaxis()->SetBinLabel(1,"total");
    effectOfCuts->GetXaxis()->SetBinLabel(2,"1");
    effectOfCuts->GetXaxis()->SetBinLabel(3,"2");
    effectOfCuts->GetXaxis()->SetBinLabel(4,"3");
    effectOfCuts->GetXaxis()->SetBinLabel(5,"4");
    effectOfCuts->GetXaxis()->SetBinLabel(6,"5");
    effectOfCuts->GetXaxis()->SetBinLabel(7,"6");
    effectOfCuts->GetXaxis()->SetBinLabel(8,"7");
    effectOfCuts->GetXaxis()->SetBinLabel(9,"8");
    effectOfCuts->GetXaxis()->SetBinLabel(10,"9");
    effectOfCuts->GetXaxis()->SetBinLabel(11,"10");
    effectOfCuts->GetXaxis()->SetBinLabel(12,"11");
    effectOfCuts->GetXaxis()->SetBinLabel(13,"12");
    effectOfCuts->GetXaxis()->SetBinLabel(14,"13");
    effectOfCuts->GetXaxis()->SetBinLabel(15,"14");
    effectOfCuts->GetXaxis()->SetBinLabel(16,"15");
    effectOfCuts->GetXaxis()->SetBinLabel(17,"16");
    effectOfCuts->GetXaxis()->SetBinLabel(18,"17");
    listout->Add(effectOfCuts);
    fDaughterHistogramArray2D[i][0] = (Int_t*)effectOfCuts;

    TH1F * effectOfCutsMC = new TH1F("effectOfCutsMC","Removal counter",18,0,18);
    effectOfCutsMC->SetStats(kTRUE);
    effectOfCutsMC->GetXaxis()->SetTitle("Cut number");
    effectOfCutsMC->GetYaxis()->SetTitle("Particles cut");
    effectOfCutsMC->GetXaxis()->SetBinLabel(1,"total");
    effectOfCutsMC->GetXaxis()->SetBinLabel(2,"1");
    effectOfCutsMC->GetXaxis()->SetBinLabel(3,"2");
    effectOfCutsMC->GetXaxis()->SetBinLabel(4,"3");
    effectOfCutsMC->GetXaxis()->SetBinLabel(5,"4");
    effectOfCutsMC->GetXaxis()->SetBinLabel(6,"5");
    effectOfCutsMC->GetXaxis()->SetBinLabel(7,"6");
    effectOfCutsMC->GetXaxis()->SetBinLabel(8,"7");
    effectOfCutsMC->GetXaxis()->SetBinLabel(9,"8");
    effectOfCutsMC->GetXaxis()->SetBinLabel(10,"9");
    effectOfCutsMC->GetXaxis()->SetBinLabel(11,"10");
    effectOfCutsMC->GetXaxis()->SetBinLabel(12,"11");
    effectOfCutsMC->GetXaxis()->SetBinLabel(13,"12");
    effectOfCutsMC->GetXaxis()->SetBinLabel(14,"13");
    effectOfCutsMC->GetXaxis()->SetBinLabel(15,"14");
    effectOfCutsMC->GetXaxis()->SetBinLabel(16,"15");
    effectOfCutsMC->GetXaxis()->SetBinLabel(17,"16");
    effectOfCutsMC->GetXaxis()->SetBinLabel(18,"17");
    listout->Add(effectOfCutsMC);
    fDaughterHistogramArray2D[i][1] = (Int_t*)effectOfCutsMC;

    TString name_particle_pdg ="particle_pdg";
    TH1F* hist_particle_pdg = new TH1F(name_particle_pdg.Data(),"Pdg code particle; pdg code; Entries",2000,-0.5,1999.5);
    hist_particle_pdg->Sumw2();
    hist_particle_pdg->SetLineColor(6);
    hist_particle_pdg->SetMarkerStyle(20);
    hist_particle_pdg->SetMarkerSize(0.6);
    hist_particle_pdg->SetMarkerColor(6);
    TH1F* histogram_particle_pdg = (TH1F*)hist_particle_pdg->Clone();
    listout->Add(histogram_particle_pdg);
    fDaughterHistogramArray2D[i][2] = (Int_t*)histogram_particle_pdg;

    TString name_particle_mother_pdg ="particle_mother_pdg";
    TH1F* hist_particle_mother_pdg = new TH1F(name_particle_mother_pdg.Data(),"Pdg code particle mother; pdg code; Entries",2000,-0.5,1999.5);
    hist_particle_mother_pdg->Sumw2();
    hist_particle_mother_pdg->SetLineColor(6);
    hist_particle_mother_pdg->SetMarkerStyle(20);
    hist_particle_mother_pdg->SetMarkerSize(0.6);
    hist_particle_mother_pdg->SetMarkerColor(6);
    TH1F* histogram_particle_mother_pdg = (TH1F*)hist_particle_mother_pdg->Clone();
    listout->Add(histogram_particle_mother_pdg);
    fDaughterHistogramArray2D[i][3] = (Int_t*)histogram_particle_mother_pdg;

    TString name_ptB0_vs_ptTrack ="ptB0_vs_ptTrack";
    TH2F* hist_ptB0_vs_ptTrack = new TH2F(name_ptB0_vs_ptTrack.Data(),"Pt B0 vs Pt ; p_{T} B0 [GeV/c]; p_{T} track [GeV/c]",100,0,30,100,0,30);
    hist_ptB0_vs_ptTrack->Sumw2();
    hist_ptB0_vs_ptTrack->SetLineColor(6);
    hist_ptB0_vs_ptTrack->SetMarkerStyle(20);
    hist_ptB0_vs_ptTrack->SetMarkerSize(0.6);
    hist_ptB0_vs_ptTrack->SetMarkerColor(6);
    TH2F* histogram_ptB0_vs_ptTrack = (TH2F*)hist_ptB0_vs_ptTrack->Clone();
    listout->Add(histogram_ptB0_vs_ptTrack);
    fDaughterHistogramArray2D[i][4] = (Int_t*)histogram_ptB0_vs_ptTrack;

    TString name_ptB0_vs_ptTrackMC ="ptB0_vs_ptTrackMC";
    TH2F* hist_ptB0_vs_ptTrackMC = new TH2F(name_ptB0_vs_ptTrackMC.Data(),"Pt B0 vs Pt ; p_{T} B0 [GeV/c]; p_{T} track [GeV/c]",100,0,30,100,0,30);
    hist_ptB0_vs_ptTrackMC->Sumw2();
    hist_ptB0_vs_ptTrackMC->SetLineColor(4);
    hist_ptB0_vs_ptTrackMC->SetMarkerStyle(20);
    hist_ptB0_vs_ptTrackMC->SetMarkerSize(0.6);
    hist_ptB0_vs_ptTrackMC->SetMarkerColor(6);
    TH2F* histogram_ptB0_vs_ptTrackMC = (TH2F*)hist_ptB0_vs_ptTrackMC->Clone();
    listout->Add(histogram_ptB0_vs_ptTrackMC);  
    fDaughterHistogramArray2D[i][5] = (Int_t*)histogram_ptB0_vs_ptTrackMC;  
  }

  //we make the histograms for the reconstructed particles
  for (Int_t i = 0; i < 6; i++){
    
    TString add_name = "";
    TList * listout;
    Int_t nHistogramSets = 0;
    if(i==0) {listout = fOutputD0; nHistogramSets = 6 + 2*fnPtBins;}
    if(i==1) {listout = fOutputDStar; nHistogramSets = 6 + 2*fnPtBins;}
    if(i==2) {listout = fOutputB0; nHistogramSets = 6 + 2*fnPtBins;}
    if(i==3) {listout = fOutputD0_D0Pt; nHistogramSets = 2*fnPtBinsD0forD0ptbin;}
    if(i==4) {listout = fOutputD0_DStarPt; nHistogramSets = 2*fnPtBinsD0forDStarptbin;}
    if(i==5) {listout = fOutputDStar_DStarPt; nHistogramSets = 2*fnPtBinsDStarforDStarptbin;}

    for (Int_t j = 0; j < nHistogramSets; j++){
      if(i<3)
      {
        if(j==0) add_name = "";
        if(j==1) add_name = "MC";
        if(j==2) add_name = "Cut";
        if(j==3) add_name = "MCCut";
        if(j==4) add_name = "Result";
        if(j==5) add_name = "MCResult";
        if(j%2==0 && j>5) {add_name = "_ptbin_"; add_name += fPtBinLimits[(j-6)/2]; add_name += "_to_"; add_name += fPtBinLimits[(j-6)/2 + 1];}
        if(j%2==1 && j>5) {add_name = "MC_ptbin_"; add_name += fPtBinLimits[(j-7)/2]; add_name += "_to_"; add_name += fPtBinLimits[(j-7)/2 + 1];}
      }
      if(i==3)
      { 
        if(j%2==0) {add_name = "_ptbin_"; add_name += fPtBinLimitsD0forD0ptbin[j/2]; add_name += "_to_"; add_name += fPtBinLimitsD0forD0ptbin[1+j/2];}
        if(j%2==1) {add_name = "MC_ptbin_"; add_name += fPtBinLimitsD0forD0ptbin[(j-1)/2]; add_name += "_to_"; add_name += fPtBinLimitsD0forD0ptbin[1+(j-1)/2];}
      }
      if(i==4)
      {
        if(j%2==0) {add_name = "_ptbin_"; add_name += fPtBinLimitsD0forDStarptbin[j/2]; add_name += "_to_"; add_name += fPtBinLimitsD0forDStarptbin[1+j/2];}
        if(j%2==1) {add_name = "MC_ptbin_"; add_name += fPtBinLimitsD0forDStarptbin[(j-1)/2]; add_name += "_to_"; add_name += fPtBinLimitsD0forDStarptbin[1+(j-1)/2];}
      }
      if(i==5)
      {
        if(j%2==0) {add_name = "_ptbin_"; add_name += fPtBinLimitsDStarforDStarptbin[j/2]; add_name += "_to_"; add_name += fPtBinLimitsDStarforDStarptbin[1+j/2];}
        if(j%2==1) {add_name = "MC_ptbin_"; add_name += fPtBinLimitsDStarforDStarptbin[(j-1)/2]; add_name += "_to_"; add_name += fPtBinLimitsDStarforDStarptbin[1+(j-1)/2];}
      }


      TString name_Histogram = "";
      TString discription_Histogram  = "";
      Int_t numberOfBins = 0;
      Double_t lowerBound = 0.0;
      Double_t upperBound = 0.0;
      Int_t numberOfBinsTwo = 0;
      Double_t lowerBoundTwo = 0.0;
      Double_t upperBoundTwo = 0.0;

      for (Int_t k = 0; k < 46; ++k)
      {
        if(k==0){name_Histogram = "ptMother"; discription_Histogram = "pt mother; p_{T} [GeV/c]; Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
        if(k==1){name_Histogram = "ptFirstDaughter"; discription_Histogram = "pt first daughter; p_{T} [GeV/c]; Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
        if(k==2){name_Histogram = "ptSecondDaughter"; discription_Histogram = "pt second daughter; p_{T} [GeV/c]; Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
        if(k==3){name_Histogram = "etaMother"; discription_Histogram = "eta mother; #eta; Entries"; numberOfBins = 100; lowerBound = -2; upperBound = 2;}
        if(k==4){name_Histogram = "phiMother"; discription_Histogram = "phi mother; #phi; Entries"; numberOfBins = 25; lowerBound = 0; upperBound = 2*TMath::Pi();}
        if(k==5){name_Histogram = "d0Mother"; discription_Histogram = "d0 mother;  [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 1.0;}
        if(k==6){name_Histogram = "d0FirstDaughter"; discription_Histogram = "d0 first daughter;  [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 1;}

        if(k==7){name_Histogram = "d0SecondDaughter"; discription_Histogram = "d0 second daughter;  [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 1;}

        if(k==8){name_Histogram = "pointingAngleMother"; discription_Histogram = "pointing angle;  [Cos(#theta)]; Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
        if(k==9){name_Histogram = "impactProduct"; discription_Histogram = "impact product; [cm^{2}]; Entries"; numberOfBins = 10000; lowerBound = -0.01; upperBound = 0.01;}
        if(k==10){name_Histogram = "impactProductXY"; discription_Histogram = "impact product XY; [cm^{2}]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 0.5;}
        if(k==11){name_Histogram = "invariantMassMother"; discription_Histogram = "mass mother candidate; m [GeV/c^{2}]; Entries"; numberOfBins = 50000; lowerBound = 0; upperBound = 10;}
        if(k==12){name_Histogram = "deltaMassMother"; discription_Histogram = "mass mother candidate; m [GeV/c^{2}]; Entries"; numberOfBins = 50000; lowerBound = 0; upperBound = 10;}
        if(k==13){name_Histogram = "dcaMother"; discription_Histogram = "dca mother; distance [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 0.25;}
        if(k==14){name_Histogram = "vertexDistance"; discription_Histogram = "vertex distance between mother and primary vertex; distance [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 1;}
        if(k==15){name_Histogram = "normDecayLength"; discription_Histogram = "Normalized decay length w.r.t primary vertex; [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 50;}
        if(k==16){name_Histogram = "pseudoProperDecayTime"; discription_Histogram = "Pseudo Proper Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 10000; lowerBound = -10; upperBound = 10;}
        if(k==17){name_Histogram = "DecayTime"; discription_Histogram = "Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 0.00000001;}
        if(k==18){name_Histogram = "normDecayTime"; discription_Histogram = "Normalized Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 0.0000001;}
        if(k==19){name_Histogram = "angleMotherFirstDaughter"; discription_Histogram = "flight angle mother and first daughter; [Cos(#phi)]; Entries"; numberOfBins = 1000; lowerBound = -1; upperBound = 1;}
        if(k==20){name_Histogram = "angleMotherSecondDaughter"; discription_Histogram = "flight angle mother and second daughter; [Cos(#phi)]; Entries"; numberOfBins = 2000; lowerBound = 0.5; upperBound = 1;}
        if(k==21){name_Histogram = "angleBetweenBothDaughters"; discription_Histogram = "angle between both daughters; [Cos(#phi)]; Entries"; numberOfBins = 1000; lowerBound = -1; upperBound = 1;}
        if(k==22){name_Histogram = "cosThetaStar"; discription_Histogram = "cosThetaStar; [Cos(#theta*)]; Entries"; numberOfBins = 1000; lowerBound = -2; upperBound = 2;}


        if(k==23){if(i==0 || i==3 || i==4){name_Histogram = "pointingAngleToDStar"; discription_Histogram = "Pointing angle w.r.t. DStar decay vertex; [Cos(#theta)]; Entries"; numberOfBins = 1000; lowerBound = -1; upperBound = 1;}
        else continue;} 
        if(k==24){if(i==0 || i==3 || i==4){name_Histogram = "d0MotherToDStar"; discription_Histogram = "d0 Mother w.r.t. DStar decay vertex; [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 1;}
        else continue;} 
        if(k==25){if(i==0 || i==3 || i==4){name_Histogram = "d0FirstDaughterToDStar"; discription_Histogram = "d0 first daughter w.r.t. DStar decay vertex; [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 1;}
        else continue;} 
        if(k==26){if(i==0 || i==3 || i==4){name_Histogram = "d0SecondDaughterToDStar"; discription_Histogram = "d0 second daughter w.r.t. DStar decay vertex; [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 1;}
        else continue;}
        if(k==27){if(i==0 || i==3 || i==4){name_Histogram = "impactProductToDStar"; discription_Histogram = "impact product w.r.t. DStar decay vertex; [cm]; Entries"; numberOfBins = 10000; lowerBound = -0.1; upperBound = 0.1;}
        else continue;} 
        if(k==28){if(i==0 || i==3 || i==4){name_Histogram = "impactProductXYToDStar"; discription_Histogram = "impact product XY w.r.t. DStar decay vertex; [cm^{2}]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 0.5;}
        else continue;} 
        if(k==29){if(i==0 || i==3 || i==4){name_Histogram = "normDecayLengthToDStar"; discription_Histogram = "Normalized decay length w.r.t. DStar decay vertex; [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 50;}
          else continue;} 
        if(k==30){if(i==0 || i==3 || i==4){name_Histogram = "pseudoProperDecayTimeToDStar"; discription_Histogram = "Pseudo Proper Decay Time w.r.t DStar vertex; [a.u.]; Entries"; numberOfBins = 10000; lowerBound = -1; upperBound = 1;}
          else continue;}  
        if(k==31){if(i==0 || i==3 || i==4){name_Histogram = "DecayTimeToDStar"; discription_Histogram = "Decay Time w.r.t DStar vertex; [a.u.]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 0.00000001;}
          else continue;}  
        if(k==32){if(i==0 || i==3 || i==4){name_Histogram = "normDecayTimeToDStar"; discription_Histogram = "Normalized Decay Time w.r.t DStar vertex; [a.u.]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 0.00000001;}
          else continue;} 

        if(k==33){name_Histogram = "topomaticFirstDaughter"; discription_Histogram = "topomatic d0 first daughter; [cm]; Entries"; numberOfBins = 10000; lowerBound = -10; upperBound = 10;}
        if(k==34){name_Histogram = "topomaticSecondDaughter"; discription_Histogram = "topomatic d0 second daughter; [cm]; Entries"; numberOfBins = 10000; lowerBound = -10; upperBound = 10;}
        if(k==35){name_Histogram = "topomaticMax"; discription_Histogram = "Max topomatic; [cm]; Entries"; numberOfBins = 10000; lowerBound = -10; upperBound = 10;}
        if(k==36){name_Histogram = "topomaticMin"; discription_Histogram = "Min topomatic; [cm]; Entries"; numberOfBins = 10000; lowerBound = -10; upperBound = 10;}

        if(k==37){name_Histogram = "ptFirstDaughter_vs_ptSecondDaughter"; discription_Histogram = "Pt first daughter vs Pt second daughter; p_{T} first daughter [GeV/c]; p_{T} second daughter [GeV/c]"; numberOfBins = 200; lowerBound = 0; upperBound = 30; numberOfBinsTwo = 200; lowerBoundTwo = 0; upperBoundTwo = 30;}
        if(k==38){name_Histogram = "ptMother_vs_ptFirstDaughter"; discription_Histogram = "Pt mother vs Pt first daughter; p_{T} Mother [GeV/c]; p_{T} Daughter [GeV/c]"; numberOfBins = 200; lowerBound = 0; upperBound = 30; numberOfBinsTwo = 200; lowerBoundTwo = 0; upperBoundTwo = 30;}
        if(k==39){name_Histogram = "ptMother_vs_ptSecondDaughter"; discription_Histogram = "Pt mother vs Pt second daughter; p_{T} Mother [GeV/c]; p_{T} Daughter [GeV/c]"; numberOfBins = 200; lowerBound = 0; upperBound = 30; numberOfBinsTwo = 200; lowerBoundTwo = 0; upperBoundTwo = 30;}
        if(k==40){name_Histogram = "ptMother_vs_d0FirstDaughter"; discription_Histogram = "Pt mother vs d0 first daughter; p_{T} Mother [GeV/c]; distance [cm]"; numberOfBins = 50; lowerBound = 0; upperBound = 30; numberOfBinsTwo = 50; lowerBoundTwo = 0; upperBoundTwo = 0.01;}
        if(k==41){name_Histogram = "ptMother_vs_d0SecondDaughter"; discription_Histogram = "Pt mother vs d0 second daughter; p_{T} Mother [GeV/c]; distance [cm]"; numberOfBins = 50; lowerBound = 0; upperBound = 30; numberOfBinsTwo = 50; lowerBoundTwo = 0; upperBoundTwo = 0.01;}
        if(k==42){name_Histogram = "ptMother_vs_PointingAngle"; discription_Histogram = "Pt mother vs Pointing Angle; p_{T} [GeV/c];  [Cos()]"; numberOfBins = 50; lowerBound = 0; upperBound = 30; numberOfBinsTwo = 50; lowerBoundTwo = 0; upperBoundTwo = 1;}
        if(k==43){name_Histogram = "ptMother_vs_ImpactProduct"; discription_Histogram = "Pt mother vs Impact product; p_{T} [GeV/c]; distance [cm]"; numberOfBins = 50; lowerBound = 0; upperBound = 30; numberOfBinsTwo = 100; lowerBoundTwo = -0.5; upperBoundTwo = 0.5;}
        if(k==44){name_Histogram = "Momentum_Mother_vs_Angle_one"; discription_Histogram = "Momentum Mother vs Angle First Daughter; P [GeV/c]; Angle [Cos(#phi)]"; numberOfBins = 50; lowerBound = 0; upperBound = 30; numberOfBinsTwo = 100; lowerBoundTwo = -1; upperBoundTwo = 1;}
        if(k==45){name_Histogram = "Momentum_Mother_vs_Angle_two"; discription_Histogram = "Momentum Mother vs Angle Second Daughter; P [GeV/c]; Angle [Cos(#phi)]"; numberOfBins = 50; lowerBound = 0; upperBound = 30; numberOfBinsTwo = 100; lowerBoundTwo = -1; upperBoundTwo = 1;}

        if(k<37){ // 1D
          name_Histogram += add_name;
          TH1F* histogram = new TH1F(name_Histogram.Data(),discription_Histogram.Data(),numberOfBins,lowerBound,upperBound);
          histogram->Sumw2();
          if(j%2==0) histogram->SetLineColor(6);
          if(j%2==1) histogram->SetLineColor(4);
          histogram->SetMarkerStyle(20);
          histogram->SetMarkerSize(0.6);
          if(j%2==0) histogram->SetMarkerColor(6);
          if(j%2==1) histogram->SetMarkerColor(4);
          TH1F* histogram_Clone = (TH1F*)histogram->Clone();
          listout->Add(histogram_Clone);
          fMotherHistogramArray[i][j][k] = (Int_t*)histogram_Clone;
        }
        if(k>36){ // 2D
          name_Histogram += add_name;
          TH2F* histogram = new TH2F(name_Histogram.Data(),discription_Histogram.Data(),numberOfBins,lowerBound,upperBound,numberOfBinsTwo,lowerBoundTwo,upperBoundTwo);
          histogram->Sumw2();
          if(j%2==0) histogram->SetLineColor(6);
          if(j%2==1) histogram->SetLineColor(4);
          histogram->SetMarkerStyle(20);
          histogram->SetMarkerSize(0.6);
          histogram->SetMarkerColor(6);
          TH2F* histogram_Clone = (TH2F*)histogram->Clone();
          listout->Add(histogram_Clone);
          fMotherHistogramArray[i][j][k] = (Int_t*)histogram_Clone;
        }
      }
    }

    TH1F * effectOfCuts = new TH1F("effectOfCuts","Removal counter",100,0,100);
    effectOfCuts->SetStats(kTRUE);
    effectOfCuts->GetXaxis()->SetTitle("Cut number");
    effectOfCuts->GetYaxis()->SetTitle("Particles cut");
    effectOfCuts->GetXaxis()->SetBinLabel(1,"total");
    for (Int_t i = 1; i < 100; ++i)
    {
      TString integerText = "";
      integerText += i;
      effectOfCuts->GetXaxis()->SetBinLabel(i+1,integerText);
    }
    listout->Add(effectOfCuts);
    fMotherHistogramArray2D[i][0] = (Int_t*)effectOfCuts;

    TH1F * effectOfCutsMC = new TH1F("effectOfCutsMC","Removal counter",100,0,100);
    effectOfCutsMC->SetStats(kTRUE);
    effectOfCutsMC->GetXaxis()->SetTitle("Cut number");
    effectOfCutsMC->GetYaxis()->SetTitle("Particles cut");
    effectOfCutsMC->GetXaxis()->SetBinLabel(1,"total");
    for (Int_t i = 1; i < 100; ++i)
    {
      TString integerText = "";
      integerText += i;
      effectOfCutsMC->GetXaxis()->SetBinLabel(i+1,integerText);
    }
    listout->Add(effectOfCutsMC);
    fMotherHistogramArray2D[i][1] = (Int_t*)effectOfCutsMC;

    TString name_particle_pdg ="particle_pdg";
    TH1F* hist_particle_pdg = new TH1F(name_particle_pdg.Data(),"Pdg code particle; pdg code; Entries",2000,-0.5,1999.5);
    hist_particle_pdg->Sumw2();
    hist_particle_pdg->SetLineColor(6);
    hist_particle_pdg->SetMarkerStyle(20);
    hist_particle_pdg->SetMarkerSize(0.6);
    hist_particle_pdg->SetMarkerColor(6);
    TH1F* histogram_particle_pdg = (TH1F*)hist_particle_pdg->Clone();
    listout->Add(histogram_particle_pdg);
    fMotherHistogramArray2D[i][2] = (Int_t*)histogram_particle_pdg;

    TString name_particle_mother_pdg ="particle_mother_pdg";
    TH1F* hist_particle_mother_pdg = new TH1F(name_particle_mother_pdg.Data(),"Pdg code particle mother; pdg code; Entries",2000,-0.5,1999.5);
    hist_particle_mother_pdg->Sumw2();
    hist_particle_mother_pdg->SetLineColor(6);
    hist_particle_mother_pdg->SetMarkerStyle(20);
    hist_particle_mother_pdg->SetMarkerSize(0.6);
    hist_particle_mother_pdg->SetMarkerColor(6);
    TH1F* histogram_particle_mother_pdg = (TH1F*)hist_particle_mother_pdg->Clone();
    listout->Add(histogram_particle_mother_pdg);
    fMotherHistogramArray2D[i][3] = (Int_t*)histogram_particle_mother_pdg;   

    TString name_distance_vertex_from_real ="distance_vertex_from_real";
    TH1F* hist_distance_vertex_from_real = new TH1F(name_distance_vertex_from_real.Data(),"Distance reconstructed vertex from real vertex; distance [cm]; Entries",100,0,1);
    hist_distance_vertex_from_real->Sumw2();
    hist_distance_vertex_from_real->SetLineColor(6);
    hist_distance_vertex_from_real->SetMarkerStyle(20);
    hist_distance_vertex_from_real->SetMarkerSize(0.6);
    hist_distance_vertex_from_real->SetMarkerColor(6);
    TH1F* histogram_distance_vertex_from_real = (TH1F*)hist_distance_vertex_from_real->Clone();
    listout->Add(histogram_distance_vertex_from_real);
    fMotherHistogramArray2D[i][4] = (Int_t*)histogram_distance_vertex_from_real;

    TString name_distance_vertex_from_real_new ="distance_vertex_from_real_new";
    TH1F* hist_distance_vertex_from_real_new = new TH1F(name_distance_vertex_from_real_new.Data(),"Distance reconstructed vertex from real vertex; distance [cm]; Entries",100,0,1);
    hist_distance_vertex_from_real_new->Sumw2();
    hist_distance_vertex_from_real_new->SetLineColor(6);
    hist_distance_vertex_from_real_new->SetMarkerStyle(20);
    hist_distance_vertex_from_real_new->SetMarkerSize(0.6);
    hist_distance_vertex_from_real_new->SetMarkerColor(6);
    TH1F* histogram_distance_vertex_from_real_new = (TH1F*)hist_distance_vertex_from_real_new->Clone();
    listout->Add(histogram_distance_vertex_from_real_new);
    fMotherHistogramArray2D[i][5] = (Int_t*)histogram_distance_vertex_from_real_new;

    TString name_momentum_resolution ="momentum_resolution";
    TH1F* hist_momentum_resolution = new TH1F(name_momentum_resolution.Data(),"Momentum resolution; difference between real and reconstructed momentum [GeV/c]; Entries",1000,0,1);
    hist_momentum_resolution->Sumw2();
    hist_momentum_resolution->SetLineColor(6);
    hist_momentum_resolution->SetMarkerStyle(20);
    hist_momentum_resolution->SetMarkerSize(0.6);
    hist_momentum_resolution->SetMarkerColor(6);
    TH1F* histogram_momentum_resolution = (TH1F*)hist_momentum_resolution->Clone();
    listout->Add(histogram_momentum_resolution);
    fMotherHistogramArray2D[i][6] = (Int_t*)histogram_momentum_resolution;    
  }

  //we make the histograms for the same sign method histograms and the pt bins
  for (Int_t k = 0; k < fnPtBins+3; ++k){
    TString ptBinMother = "";
    if(k==0) ptBinMother = "";
    if(k==1) ptBinMother = "_ptbin_6_to_inf";
    if(k==2) ptBinMother = "_ptbin_3_to_inf";
    if(k>2) {ptBinMother += "_ptbin_"; ptBinMother += fPtBinLimits[k-3]; ptBinMother += "_to_"; ptBinMother += fPtBinLimits[k-2];}

    TString name_deltaMassMother ="deltaInvMassB0";
    name_deltaMassMother += ptBinMother;
    TH1F* hist_deltaMassMother = new TH1F(name_deltaMassMother.Data(),"mass mother candidate; m [GeV/c^2]; Entries",2000,0,20);
    hist_deltaMassMother->Sumw2();
    hist_deltaMassMother->SetLineColor(6);
    hist_deltaMassMother->SetMarkerStyle(20);
    hist_deltaMassMother->SetMarkerSize(0.6);
    hist_deltaMassMother->SetMarkerColor(6);
    TH1F* histogram_deltaMassMother = (TH1F*)hist_deltaMassMother->Clone();
    fOutputB0MC->Add(histogram_deltaMassMother);
   
    for (Int_t i = 0; i < 3; ++i){
      TString signName = "";
      if(i==0) signName = "";
      if(i==1) signName = "_SameSign";
      if(i==2) signName = "_SignSum";
      TString name_invariantMassMother ="invariantMassB0";
      name_invariantMassMother += ptBinMother + signName;
      TH1F* hist_invariantMassMother = new TH1F(name_invariantMassMother.Data(),"mass mother candidate; m [GeV/c^2]; Entries",2000,0,20);
      hist_invariantMassMother->Sumw2();
      hist_invariantMassMother->SetLineColor(6);
      hist_invariantMassMother->SetMarkerStyle(20);
      hist_invariantMassMother->SetMarkerSize(0.6);
      hist_invariantMassMother->SetMarkerColor(6);
      TH1F* histogram_invariantMassMother = (TH1F*)hist_invariantMassMother->Clone();
      fOutputB0MC->Add(histogram_invariantMassMother);
    }
  }

  TString name_cutEffectBackground ="cutEffectBackground";
  TH2I* hist_cutEffectBackground = new TH2I(name_cutEffectBackground.Data(),"Effect of Cuts on background; cut number; cut number",99,0,99,99,0,99);
  for (int i = 0; i < 99; ++i)
  {
    TString integerText = "";
    integerText += i;
    hist_cutEffectBackground->GetXaxis()->SetBinLabel(i+1,integerText);
    hist_cutEffectBackground->GetYaxis()->SetBinLabel(i+1,integerText);
  }
  TH2I* histogram_cutEffectBackground = (TH2I*)hist_cutEffectBackground->Clone();
  fOutputB0MC->Add(histogram_cutEffectBackground);

  TString name_cutEffectSignal ="cutEffectSignal";
  TH2I* hist_cutEffectSignal = new TH2I(name_cutEffectSignal.Data(),"Effect of Cuts on Signal; cut number; cut number",99,0,99,99,0,99);
  for (Int_t i = 0; i < 99; ++i)
  {
    TString integerText = "";
    integerText += i;
    hist_cutEffectSignal->GetXaxis()->SetBinLabel(i+1,integerText);
    hist_cutEffectSignal->GetYaxis()->SetBinLabel(i+1,integerText);
  }
  TH2I* histogram_cutEffectSignal = (TH2I*)hist_cutEffectSignal->Clone();
  fOutputB0MC->Add(histogram_cutEffectSignal);

  TString name_cutEffectUniqueBackground ="cutEffectUniqueBackground";
  TH1I* hist_cutEffectUniqueBackground = new TH1I(name_cutEffectUniqueBackground.Data(),"Effect of Cuts on Signal; cut number; cut number",99,0,99);
  for (Int_t i = 0; i < 99; ++i)
  {
    TString integerText = "";
    integerText += i;
    hist_cutEffectUniqueBackground->GetXaxis()->SetBinLabel(i+1,integerText);
  }
  TH1I* histogram_cutEffectUniqueBackground = (TH1I*)hist_cutEffectUniqueBackground->Clone();
  fOutputB0MC->Add(histogram_cutEffectUniqueBackground);

  TString name_cutEffectUniqueSignal ="cutEffectUniqueSignal";
  TH1I* hist_cutEffectUniqueSignal = new TH1I(name_cutEffectUniqueSignal.Data(),"Effect of Cuts on Signal; cut number; cut number",99,0,99);
  for (Int_t i = 0; i < 99; ++i)
  {
    TString integerText = "";
    integerText += i;
    hist_cutEffectUniqueSignal->GetXaxis()->SetBinLabel(i+1,integerText);
  }
  TH1I* histogram_cutEffectUniqueSignal = (TH1I*)hist_cutEffectUniqueSignal->Clone();
  fOutputB0MC->Add(histogram_cutEffectUniqueSignal);

  TString name_totalITSBackground ="totalITSBackground";
  TH1F* hist_totalITSBackground = new TH1F(name_totalITSBackground.Data(),"Total nr. of ITS hits for the daughters; number [#]; Entries",30,0,30);
  hist_totalITSBackground->Sumw2();
  hist_totalITSBackground->SetLineColor(6);
  hist_totalITSBackground->SetMarkerStyle(20);
  hist_totalITSBackground->SetMarkerSize(0.6);
  hist_totalITSBackground->SetMarkerColor(6);
  TH1F* histogram_totalITSBackground = (TH1F*)hist_totalITSBackground->Clone();
  fOutputB0MC->Add(histogram_totalITSBackground);

  TString name_totalITSSignal ="totalITSSignal";
  TH1F* hist_totalITSSignal = new TH1F(name_totalITSSignal.Data(),"Total nr. of ITS hits for the daughters; number [#]; Entries",30,0,30);
  hist_totalITSSignal->Sumw2();
  hist_totalITSSignal->SetLineColor(6);
  hist_totalITSSignal->SetMarkerStyle(20);
  hist_totalITSSignal->SetMarkerSize(0.6);
  hist_totalITSSignal->SetMarkerColor(6);
  TH1F* histogram_totalITSSignal = (TH1F*)hist_totalITSSignal->Clone();
  fOutputB0MC->Add(histogram_totalITSSignal);

  TString name_totalTPCBackground ="totalTPCBackground";
  TH1F* hist_totalTPCBackground = new TH1F(name_totalTPCBackground.Data(),"Total nr. of TPC hits for the daughters; number [#]; Entries",1000,0,1000);
  hist_totalTPCBackground->Sumw2();
  hist_totalTPCBackground->SetLineColor(6);
  hist_totalTPCBackground->SetMarkerStyle(20);
  hist_totalTPCBackground->SetMarkerSize(0.6);
  hist_totalTPCBackground->SetMarkerColor(6);
  TH1F* histogram_totalTPCBackground = (TH1F*)hist_totalTPCBackground->Clone();
  fOutputB0MC->Add(histogram_totalTPCBackground);

  TString name_totalTPCSignal ="totalTPCSignal";
  TH1F* hist_totalTPCSignal = new TH1F(name_totalTPCSignal.Data(),"Total nr. of TPC hits for the daughters; number [#]; Entries",1000,0,1000);
  hist_totalTPCSignal->Sumw2();
  hist_totalTPCSignal->SetLineColor(6);
  hist_totalTPCSignal->SetMarkerStyle(20);
  hist_totalTPCSignal->SetMarkerSize(0.6);
  hist_totalTPCSignal->SetMarkerColor(6);
  TH1F* histogram_totalTPCSignal = (TH1F*)hist_totalTPCSignal->Clone();
  fOutputB0MC->Add(histogram_totalTPCSignal);

  TString name_totalSigmaPIDBackground ="totalSigmaPIDBackground";
  TH1F* hist_totalSigmaPIDBackground = new TH1F(name_totalSigmaPIDBackground.Data(),"Total sigma of TPC and TOF PID for the daughters; number [#]; Entries",1000,0,100);
  hist_totalSigmaPIDBackground->Sumw2();
  hist_totalSigmaPIDBackground->SetLineColor(6);
  hist_totalSigmaPIDBackground->SetMarkerStyle(20);
  hist_totalSigmaPIDBackground->SetMarkerSize(0.6);
  hist_totalSigmaPIDBackground->SetMarkerColor(6);
  TH1F* histogram_totalSigmaPIDBackground = (TH1F*)hist_totalSigmaPIDBackground->Clone();
  fOutputB0MC->Add(histogram_totalSigmaPIDBackground);

  TString name_totalSigmaPIDSignal ="totalSigmaPIDSignal";
  TH1F* hist_totalSigmaPIDSignal = new TH1F(name_totalSigmaPIDSignal.Data(),"Total sigma of TPC and TOF PID for the daughters; number [#]; Entries",1000,0,100);
  hist_totalSigmaPIDSignal->Sumw2();
  hist_totalSigmaPIDSignal->SetLineColor(6);
  hist_totalSigmaPIDSignal->SetMarkerStyle(20);
  hist_totalSigmaPIDSignal->SetMarkerSize(0.6);
  hist_totalSigmaPIDSignal->SetMarkerColor(6);
  TH1F* histogram_totalSigmaPIDSignal = (TH1F*)hist_totalSigmaPIDSignal->Clone();
  fOutputB0MC->Add(histogram_totalSigmaPIDSignal);
  return;
}
//-------------------------------------------------------------------------------------
AliAODVertex* AliAnalysisTaskSEB0toDStarPi::RecalculateVertex(const AliVVertex *primary,TObjArray *tracks,Double_t bField, Int_t finderAlgorithm){
  //
  // Helper function to recalculate a vertex.
  //
  
  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;

  AliVertexerTracks *vertexer = new AliVertexerTracks(bField);
  // vertexer->SetFinderAlgorithm(finderAlgorithm);

  vertexer->SetVtxStart((AliESDVertex*)primary);//primary vertex
  vertexESD = (AliESDVertex*)vertexer->VertexForSelectedESDTracks(tracks);

  delete vertexer; vertexer=NULL;

  if(!vertexESD) return vertexAOD;


  if(vertexESD->GetNContributors()!=tracks->GetEntriesFast()) {

    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }

  // Double_t pos[3],cov[6],chi2perNDF;
  // vertexESD->GetXYZ(pos);
  // convert to AliAODVertex
  //
  Double_t dispersion;
  Double_t pos[3],cov[6],chi2perNDF;
  for(Int_t a=0;a<3;a++)pos[a]=0.;
  for(Int_t b=0;b<6;b++)cov[b]=0.;
  chi2perNDF=0;
  //
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix


  Double_t vertRadius2=pos[0]*pos[0]+pos[1]*pos[1];
  if(vertRadius2>8.){//(2.82)^2 radius beam pipe

    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }
  
  chi2perNDF = vertexESD->GetChi2toNDF();
  dispersion = vertexESD->GetDispersion();
  delete vertexESD; vertexESD=NULL;
  Int_t nprongs= 2;//tracks->GetEntriesFast();
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF,0x0,-1,AliAODVertex::kUndef,nprongs);
  
  return vertexAOD;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::B0toDStarPiSignalTracksInMC(TClonesArray * mcTrackArray,AliAODEvent*  aodevent,TMatrix * B0toDStarPiLabelMatrix, TList *listout){

  TMatrix &particleMatrix = *B0toDStarPiLabelMatrix;
  for (Int_t i=0; i<mcTrackArray->GetEntriesFast(); i++){ 

    Int_t mcLabelPionB0 = 0;
    Int_t mcLabelPionDStar = 0;
    Int_t mcLabelPionD0 = 0;
    Int_t mcLabelKaon = 0;
    Int_t mcLabelD0 = 0;
    Int_t mcLabelDStar = 0;
    Int_t mcLabelB0 = 0;

    Double_t ptMC[7] = {0.0};


    Bool_t mcPionB0Present = kFALSE;
    Bool_t mcPionDStarPresent = kFALSE;
    Bool_t mcPionD0Present = kFALSE;
    Bool_t mcKaonPresent = kFALSE;

    // Below, we find all the MC labels for the true signal tracks
    AliAODMCParticle *mcTrackParticle = dynamic_cast< AliAODMCParticle*>(mcTrackArray->At(i));
    Int_t pdgCodeMC=TMath::Abs(mcTrackParticle->GetPdgCode());
    
    if (pdgCodeMC==511){ //if the track is a B0 we look at its daughters

      mcLabelB0 = mcTrackParticle->GetLabel();
      Int_t nDaughterB0 = mcTrackParticle->GetNDaughters();
      ptMC[0] = mcTrackParticle->Pt();

      TString fillthis= "B0s_in_analysis";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(0);

      if(nDaughterB0==2){
        for(Int_t iDaughterB0=0; iDaughterB0<2; iDaughterB0++){

          AliAODMCParticle* daughterB0 = (AliAODMCParticle*)mcTrackArray->At(mcTrackParticle->GetDaughter(iDaughterB0));
          if(!daughterB0) break;
          Int_t pdgCodeDaughterB0=TMath::Abs(daughterB0->GetPdgCode());

          if (pdgCodeDaughterB0==211){ //if the track is a pion we save its monte carlo label
            mcLabelPionB0 = daughterB0->GetLabel();
            mcPionB0Present = kTRUE;
            ptMC[1] = daughterB0->Pt();

          } else if (pdgCodeDaughterB0==413){ //if the track is a DStar we look at its daughters
            mcLabelDStar = daughterB0->GetLabel();
            Int_t nDaughterDStar = daughterB0->GetNDaughters();
            ptMC[2] = daughterB0->Pt();

            if(nDaughterDStar==2){
              for(Int_t iDaughterDStar=0; iDaughterDStar<2; iDaughterDStar++){

                AliAODMCParticle* daughterDStar = (AliAODMCParticle*)mcTrackArray->At(daughterB0->GetDaughter(iDaughterDStar));
                if(!daughterDStar) break;
                Int_t pdgCodeDaughterDStar=TMath::Abs(daughterDStar->GetPdgCode());

                if (pdgCodeDaughterDStar==211){ //if the track is a pion we save its monte carlo label
                  mcLabelPionDStar = daughterDStar->GetLabel();
                  mcPionDStarPresent = kTRUE;
                  ptMC[3] = daughterDStar->Pt();

                } else if (pdgCodeDaughterDStar==421){ //if the track is a D0 we look at its daughters
                  mcLabelD0 = daughterDStar->GetLabel();
                  Int_t nDaughterD0 = daughterDStar->GetNDaughters();
                  ptMC[4] = daughterDStar->Pt();

                  if(nDaughterD0==2){
                    for(Int_t iDaughterD0=0; iDaughterD0<2; iDaughterD0++){

                      AliAODMCParticle* daughterD0 = (AliAODMCParticle*)mcTrackArray->At(daughterDStar->GetDaughter(iDaughterD0));
                      if(!daughterD0) break;
                      Int_t pdgCodeDaughterD0=TMath::Abs(daughterD0->GetPdgCode());

                      if (pdgCodeDaughterD0==211){ //if the track is a pion we save its monte carlo label
                        mcLabelPionD0 = daughterD0->GetLabel();
                        ptMC[5] = daughterD0->Pt();
                        mcPionD0Present = kTRUE;

                      } else if (pdgCodeDaughterD0==321){ //if the track is a kaon we save its monte carlo label
                        mcLabelKaon = daughterD0->GetLabel();
                        mcKaonPresent = kTRUE;
                        ptMC[6] = daughterD0->Pt();

                      } else break;
                    }
                  } 
                } else break;
              }
            }
          } else break;
        }
      }
    }
    // Next, we save the labels to our array
    if (mcPionB0Present && mcPionDStarPresent && mcPionD0Present && mcKaonPresent){
      Int_t rows = B0toDStarPiLabelMatrix->GetNrows();

      B0toDStarPiLabelMatrix->ResizeTo(rows+1,7);
      particleMatrix(rows,0) = mcLabelPionB0;
      particleMatrix(rows,1) = mcLabelPionDStar;
      particleMatrix(rows,2) = mcLabelPionD0;
      particleMatrix(rows,3) = mcLabelKaon;
      particleMatrix(rows,4) = mcLabelD0;
      particleMatrix(rows,5) = mcLabelDStar;
      particleMatrix(rows,6) = mcLabelB0;

      // We also save information on the amount of signal tracks that exist in the MC dataset
      TString fillthis= "B0s_in_analysis";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(1);

      fillthis= "B0s_per_bin";
      if(0.0 < ptMC[0] && ptMC[0] < 3.0) {((TH1F*)(listout->FindObject(fillthis)))->Fill(0);}
      if(3.0 < ptMC[0] && ptMC[0] < 6.0) {((TH1F*)(listout->FindObject(fillthis)))->Fill(1);}
      if(6.0 < ptMC[0] && ptMC[0] < 10.0) {((TH1F*)(listout->FindObject(fillthis)))->Fill(2);}
      if(10.0 < ptMC[0] && ptMC[0] < 18.0) {((TH1F*)(listout->FindObject(fillthis)))->Fill(3);}
      if(18.0 < ptMC[0] && ptMC[0] < 30.0) {((TH1F*)(listout->FindObject(fillthis)))->Fill(4);}
      if(30.0 < ptMC[0]) {((TH1F*)(listout->FindObject(fillthis)))->Fill(5);}

      fillthis= "mc_B0_pt";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(ptMC[0]);
      fillthis= "mc_B0_pion_pt";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(ptMC[1]);
      fillthis= "mc_DStar_pt";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(ptMC[2]);
      fillthis= "mc_DStar_pion_pt";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(ptMC[3]);
      fillthis= "mc_D0_pt";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(ptMC[4]);
      fillthis= "mc_D0_pion_pt";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(ptMC[5]);
      fillthis= "mc_D0_kaon_pt";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(ptMC[6]);

    }
  }

  // Not all the tracks can be/are detected by the detector. We are only interested in tracks that are detectable and lie within the acceptance of the detector.
  // We remove the undetectable tracks from the array in order to get accurate information on the amount of signal that lies within the acceptance of our detector.
  Int_t numberOfB0s = 0;
  TArrayI correctLabelArray; 
  for (Int_t i = 0; i < B0toDStarPiLabelMatrix->GetNrows(); i++)
  {
    Int_t particleCounter = 0;
    for (Int_t j = 0; j < 4; j++)
    {
      Int_t labelParticleInList = (Int_t)particleMatrix(i,j);
      for (Int_t k=0; k<aodevent->GetNumberOfTracks(); k++)
      { 
        AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(aodevent->GetTrack(k));
        if(!aodTrack) AliFatal("Not a standard AOD");
        if(TMath::Abs(aodTrack->Eta())>0.9) continue;
        if(aodTrack->GetITSNcls() < 1) continue;
        if(aodTrack->GetTPCNcls() < 1) continue;
        if(aodTrack->GetLabel() == labelParticleInList) 
        {
          particleCounter++;
          break;
        }
      }
    }

    if (particleCounter==4)
    {
      TString fillthis= "B0s_in_analysis";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(2);
      numberOfB0s++;
      correctLabelArray.Set(numberOfB0s);
      correctLabelArray.AddAt(i,numberOfB0s-1);
    }
  }

  for (Int_t i = 0; i < correctLabelArray.GetSize(); i++)
  {
    particleMatrix(i,0) = (Int_t)particleMatrix(correctLabelArray[i],0);
    particleMatrix(i,1) = (Int_t)particleMatrix(correctLabelArray[i],1);
    particleMatrix(i,2) = (Int_t)particleMatrix(correctLabelArray[i],2);
    particleMatrix(i,3) = (Int_t)particleMatrix(correctLabelArray[i],3);
    particleMatrix(i,4) = (Int_t)particleMatrix(correctLabelArray[i],4);
    particleMatrix(i,5) = (Int_t)particleMatrix(correctLabelArray[i],5);
    particleMatrix(i,6) = (Int_t)particleMatrix(correctLabelArray[i],6);
  }
  B0toDStarPiLabelMatrix->ResizeTo(correctLabelArray.GetSize(),7);
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::D0PionSelection(AliAODEvent* aodEvent,TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix){
  
  //we keep track of the number of particles we could use and how many we actually use after cuts
  Int_t numberofparticles = 0;
  Int_t numberofparticlesused = 0;
  Int_t iClonesArray = 0;

  TString fillthis = "";

  //we loop over all tracks in the event
  for (Int_t i=0; i<aodEvent->GetNumberOfTracks(); i++){ 
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));
    if(!aodTrack) AliFatal("Not a standard AOD");

    //quick quality cut
    if(aodTrack->GetITSNcls() < 1) continue;
    if(aodTrack->GetTPCNcls() < 1) continue;

    Int_t mcLabelParticle = -1;
    Int_t pdgParticle = -1;
    mcLabelParticle = aodTrack->GetLabel();

    numberofparticles++;

    // we fill histograms with information of the track
    Double_t pt_track = aodTrack->Pt();
    Double_t momentum_track = aodTrack->P();
    Double_t energy_track = aodTrack->E(0.13957);
    Double_t dca_track = aodTrack->DCA();
    Double_t momentum_dca_track = aodTrack->PAtDCA();
    Int_t numberOfITS = aodTrack->GetITSNcls(); 
    Int_t numberOfTPC = aodTrack->GetTPCNcls();
    Double_t nSigmaTPC = 0;
    Double_t nSigmaTOF = 0;
    Int_t pionPIDnumber = 2;
    Int_t kaonPIDnumber = 3;
    Int_t TPCok = 0;
    Int_t TOFok = 0;

    AliAODPidHF* trackPIDHF = (AliAODPidHF*)fCuts->GetPidHF();
    TPCok = trackPIDHF->GetnSigmaTPC(aodTrack, pionPIDnumber, nSigmaTPC);
    TOFok = trackPIDHF->GetnSigmaTOF(aodTrack, pionPIDnumber, nSigmaTOF);

    Int_t daughterType = 0;
    Int_t histType;

    histType = 0;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(aodTrack->Pt(),aodTrack->DCA());

    //we check if the particle is a signal track
    Bool_t isDesiredCandidate = kFALSE;
    if(fUseMCInfo){
      TMatrix &particleMatrix = *B0toDStarPiLabelMatrix;
      for (Int_t k = 0; k < B0toDStarPiLabelMatrix->GetNrows(); ++k){
        if(mcLabelParticle == (Int_t)particleMatrix(k,2)){
          isDesiredCandidate = kTRUE;
          break;
        }
      }
    }

    if(isDesiredCandidate)
    {
      histType = 2;
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
      if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
      if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
      if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
      ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(aodTrack->Pt(),aodTrack->DCA());
    }

    //we apply a number of cuts on the particle
    Bool_t bCut = kFALSE;

    //we apply a PID cut, 2 is used to indicate we look for a pion
    if(!(fCuts->SelectPID(aodTrack,2))){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[0][1])->Fill(2);
      } else ((TH1F*)fDaughterHistogramArray2D[0][0])->Fill(2);
      bCut = kTRUE;
    }

    if(aodTrack->GetITSNcls() < fCuts->GetMinITSNclsD0Pion()){
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[0][1])->Fill(3);
      } else ((TH1F*)fDaughterHistogramArray2D[0][0])->Fill(3);
      bCut = kTRUE;
    }

    if(aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsD0Pion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[0][1])->Fill(4);
      } else ((TH1F*)fDaughterHistogramArray2D[0][0])->Fill(4);
      bCut = kTRUE;
    }

    if(fCuts->UseITSRefitD0Pion()==kTRUE){
      if(!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArray2D[0][1])->Fill(5);
        } else ((TH1F*)fDaughterHistogramArray2D[0][0])->Fill(5);
        bCut = kTRUE;
      }
    }

    if(fCuts->UseTPCRefitD0Pion()==kTRUE){
      if((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArray2D[0][1])->Fill(6);
        } else ((TH1F*)fDaughterHistogramArray2D[0][0])->Fill(6);
        bCut = kTRUE;
      }
    }

    if(fCuts->UseFilterBitD0Pion()==kTRUE){
      if(!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitD0Pion())))) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArray2D[0][1])->Fill(7);
        } else ((TH1F*)fDaughterHistogramArray2D[0][0])->Fill(7);
        bCut = kTRUE;
      }
    }

    if(aodTrack->Pt() < fCuts->GetMinPtD0Pion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[0][1])->Fill(8);
      } else ((TH1F*)fDaughterHistogramArray2D[0][0])->Fill(8);
      bCut = kTRUE;
    }

    if(!isDesiredCandidate && fQuickSignalAnalysis) bCut = kTRUE;

    if(bCut) {
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[0][1])->Fill(0);
      } else ((TH1F*)fDaughterHistogramArray2D[0][0])->Fill(0);
      continue;
    }

    //we fill histograms with track information of the tracks that pass the cuts
    histType = 1;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(aodTrack->Pt(),aodTrack->DCA());

    fD0PionTracks->push_back(i);
    numberofparticlesused++;
  }
  
  ((TH1F*)fDaughterHistogramArray[0][0][12])->Fill(numberofparticles);
  ((TH1F*)fDaughterHistogramArray[0][1][12])->Fill(numberofparticlesused);
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::D0KaonSelection(AliAODEvent* aodEvent,TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix){
  
  //we keep track of the number of particles we could use and how many we actually use after cuts
  Int_t numberofparticles = 0;
  Int_t numberofparticlesused = 0;
  Int_t iClonesArray = 0;

  TString fillthis = "";

  //we loop over all tracks in the event
  for (Int_t i=0; i<aodEvent->GetNumberOfTracks(); i++){ 
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));
    if(!aodTrack) AliFatal("Not a standard AOD");

    //quick quality cut
    if(aodTrack->GetITSNcls() < 1) continue;
    if(aodTrack->GetTPCNcls() < 1) continue;

    Int_t mcLabelParticle = -1;
    Int_t pdgParticle = -1;
    mcLabelParticle = aodTrack->GetLabel();

    numberofparticles++;

    // we fill histograms with information of the track
    Double_t pt_track = aodTrack->Pt();
    Double_t momentum_track = aodTrack->P();
    Double_t energy_track = aodTrack->E(0.4937);
    Double_t dca_track = aodTrack->DCA();
    Double_t momentum_dca_track = aodTrack->PAtDCA();
    Int_t numberOfITS = aodTrack->GetITSNcls(); 
    Int_t numberOfTPC = aodTrack->GetTPCNcls();
    Double_t nSigmaTPC = 0;
    Double_t nSigmaTOF = 0;
    Int_t pionPIDnumber = 2;
    Int_t kaonPIDnumber = 3;
    Int_t TPCok = 0;
    Int_t TOFok = 0;

    AliAODPidHF* trackPIDHF = (AliAODPidHF*)fCuts->GetPidHF();
    TPCok = trackPIDHF->GetnSigmaTPC(aodTrack, kaonPIDnumber, nSigmaTPC);
    TOFok = trackPIDHF->GetnSigmaTOF(aodTrack, kaonPIDnumber, nSigmaTOF);

    Int_t daughterType = 1;
    Int_t histType;

    histType = 0;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(aodTrack->Pt(),aodTrack->DCA());


    //we check if the particle is a signal track
    Bool_t isDesiredCandidate = kFALSE;
    if(fUseMCInfo){
      TMatrix &particleMatrix = *B0toDStarPiLabelMatrix;
      for (Int_t k = 0; k < B0toDStarPiLabelMatrix->GetNrows(); ++k){
        if(mcLabelParticle == (Int_t)particleMatrix(k,3)){
          isDesiredCandidate = kTRUE;
          break;
        }
      }
    }

    if(isDesiredCandidate)
    {
      histType = 2;
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
      if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
      if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
      if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
      ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(aodTrack->Pt(),aodTrack->DCA());
    }

    //we apply a number of cuts on the particle
    Bool_t bCut = kFALSE;

    //we apply a PID cut, 3 is used to indicate we look for a kaon
    if(!(fCuts->SelectPID(aodTrack,3))){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[1][1])->Fill(2);
      } else ((TH1F*)fDaughterHistogramArray2D[1][0])->Fill(2);
      bCut = kTRUE;
    }

    if(aodTrack->GetITSNcls() < fCuts->GetMinITSNclsD0Kaon()){
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[1][1])->Fill(3);
      } else ((TH1F*)fDaughterHistogramArray2D[1][0])->Fill(3);
      bCut = kTRUE;
    }

    if(aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsD0Kaon()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[1][1])->Fill(4);
      } else ((TH1F*)fDaughterHistogramArray2D[1][0])->Fill(4);
      bCut = kTRUE;
    }

    if(fCuts->UseITSRefitD0Kaon()==kTRUE){
      if(!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArray2D[1][1])->Fill(5);
        } else ((TH1F*)fDaughterHistogramArray2D[1][0])->Fill(5);
        bCut = kTRUE;
      }
    }

    if(fCuts->UseTPCRefitD0Kaon()==kTRUE){
      if((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArray2D[1][1])->Fill(6);
        } else ((TH1F*)fDaughterHistogramArray2D[1][0])->Fill(6);
        bCut = kTRUE;
      }
    }

    if(fCuts->UseFilterBitD0Kaon()==kTRUE){
      if(!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitD0Kaon())))) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArray2D[1][1])->Fill(7);
        } else ((TH1F*)fDaughterHistogramArray2D[1][0])->Fill(7);
        bCut = kTRUE;
      }
    }

    if(aodTrack->Pt() < fCuts->GetMinPtD0Kaon()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[1][1])->Fill(8);
      } else ((TH1F*)fDaughterHistogramArray2D[1][0])->Fill(8);
      bCut = kTRUE;
    }

    if(!isDesiredCandidate && fQuickSignalAnalysis) bCut = kTRUE;

    if(bCut) {
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[1][1])->Fill(0);
      } else ((TH1F*)fDaughterHistogramArray2D[1][0])->Fill(0);
      continue;
    }

    //we fill histograms with track information of the tracks that pass the cuts
    histType = 1;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(aodTrack->Pt(),aodTrack->DCA());

    fD0KaonTracks->push_back(i);
    numberofparticlesused++;
  }

  ((TH1F*)fDaughterHistogramArray[1][0][12])->Fill(numberofparticles);
  ((TH1F*)fDaughterHistogramArray[1][1][12])->Fill(numberofparticlesused);
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::DStarPionSelection(AliAODEvent* aodEvent,TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix){
  
  //we keep track of the number of particles we could use and how many we actually use after cuts
  Int_t numberofparticles = 0;
  Int_t numberofparticlesused = 0;
  Int_t iClonesArray = 0;

  TString fillthis = "";

  //we loop over all tracks in the event
  for (Int_t i=0; i<aodEvent->GetNumberOfTracks(); i++){ 
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));
    if(!aodTrack) AliFatal("Not a standard AOD");

    //quick quality cut
    if(aodTrack->GetITSNcls() < 1) continue;
    if(aodTrack->GetTPCNcls() < 1) continue;

    Int_t mcLabelParticle = -1;
    Int_t pdgParticle = -1;
    mcLabelParticle = aodTrack->GetLabel();

    numberofparticles++;

    // we fill histograms with information of the track
    Double_t pt_track = aodTrack->Pt();
    Double_t momentum_track = aodTrack->P();
    Double_t energy_track = aodTrack->E(0.13957);
    Double_t dca_track = aodTrack->DCA();
    Double_t momentum_dca_track = aodTrack->PAtDCA();
    Int_t numberOfITS = aodTrack->GetITSNcls(); 
    Int_t numberOfTPC = aodTrack->GetTPCNcls();
    Double_t nSigmaTPC = 0;
    Double_t nSigmaTOF = 0;
    Int_t pionPIDnumber = 2;
    Int_t kaonPIDnumber = 3;
    Int_t TPCok = 0;
    Int_t TOFok = 0;

    AliAODPidHF* trackPIDHF = (AliAODPidHF*)fCuts->GetPidHF();
    TPCok = trackPIDHF->GetnSigmaTPC(aodTrack, pionPIDnumber, nSigmaTPC);
    TOFok = trackPIDHF->GetnSigmaTOF(aodTrack, pionPIDnumber, nSigmaTOF);

    Int_t daughterType = 2;
    Int_t histType;

    histType = 0;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(aodTrack->Pt(),aodTrack->DCA());


    //we check if the particle is a signal track
    Bool_t isDesiredCandidate = kFALSE;
    if(fUseMCInfo){
      TMatrix &particleMatrix = *B0toDStarPiLabelMatrix;
      for (Int_t k = 0; k < B0toDStarPiLabelMatrix->GetNrows(); ++k){
        if(mcLabelParticle == (Int_t)particleMatrix(k,1)){
          isDesiredCandidate = kTRUE;
          break;
        }
      }
    }

    if(isDesiredCandidate)
    {
      histType = 2;
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
      if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
      if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
      if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
      ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(aodTrack->Pt(),aodTrack->DCA());
    }

    //we apply a number of cuts on the particle
     Bool_t bCut = kFALSE;

    //we apply a PID cut, 2 is used to indicate we look for a pion
    if(!(fCuts->SelectPID(aodTrack,2))){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[2][1])->Fill(2);
      } else ((TH1F*)fDaughterHistogramArray2D[2][0])->Fill(2);
      bCut = kTRUE;
    }

    if(aodTrack->GetITSNcls() < fCuts->GetMinITSNclsDStarPion()){
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[2][1])->Fill(3);
      } else ((TH1F*)fDaughterHistogramArray2D[2][0])->Fill(3);
      bCut = kTRUE;
    }

    if(aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsDStarPion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[2][1])->Fill(4);
      } else ((TH1F*)fDaughterHistogramArray2D[2][0])->Fill(4);
      bCut = kTRUE;
    }

    if(fCuts->UseITSRefitDStarPion()==kTRUE){
      if(!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArray2D[2][1])->Fill(5);
        } else ((TH1F*)fDaughterHistogramArray2D[2][0])->Fill(5);
        bCut = kTRUE;
      }
    }

    if(fCuts->UseTPCRefitDStarPion()==kTRUE){
      if((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArray2D[2][1])->Fill(6);
        } else ((TH1F*)fDaughterHistogramArray2D[2][0])->Fill(6);
        bCut = kTRUE;
      }
    }

    if(fCuts->UseFilterBitDStarPion()==kTRUE){
      if(!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitDStarPion())))) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArray2D[2][1])->Fill(7);
        } else ((TH1F*)fDaughterHistogramArray2D[2][0])->Fill(7);
        bCut = kTRUE;
      }
    }

    if(aodTrack->Pt() < fCuts->GetMinPtDStarPion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[2][1])->Fill(8);
      } else ((TH1F*)fDaughterHistogramArray2D[2][0])->Fill(8);
      bCut = kTRUE;
    }

    if(!isDesiredCandidate && fQuickSignalAnalysis) bCut = kTRUE;

    if(bCut) {
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[2][1])->Fill(0);
      } else ((TH1F*)fDaughterHistogramArray2D[2][0])->Fill(0);
      continue;
    }

    //we fill histograms with track information of the tracks that pass the cuts
    histType = 1;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(aodTrack->Pt(),aodTrack->DCA());

    fDStarPionTracks->push_back(i);
    numberofparticlesused++;
  }
  
  ((TH1F*)fDaughterHistogramArray[2][0][12])->Fill(numberofparticles);
  ((TH1F*)fDaughterHistogramArray[2][1][12])->Fill(numberofparticlesused);
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::B0PionSelection(AliAODEvent* aodEvent,TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix){
  
  //we keep track of the number of particles we could use and how many we actually use after cuts
  Int_t numberofparticles = 0;
  Int_t numberofparticlesused = 0;
  Int_t iClonesArray = 0;

  TString fillthis = "";

  //we loop over all tracks in the event
  for (Int_t i=0; i<aodEvent->GetNumberOfTracks(); i++){ 
    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(i));
    if(!aodTrack) AliFatal("Not a standard AOD");

    //quick quality cut
    if(aodTrack->GetITSNcls() < 1) continue;
    if(aodTrack->GetTPCNcls() < 1) continue;

    Int_t mcLabelParticle = -1;
    Int_t pdgParticle = -1;
    mcLabelParticle = aodTrack->GetLabel();

    numberofparticles++;

    // we fill histograms with information of the track
    Double_t pt_track = aodTrack->Pt();
    Double_t momentum_track = aodTrack->P();
    Double_t energy_track = aodTrack->E(0.13957);
    Double_t dca_track = aodTrack->DCA();
    Double_t momentum_dca_track = aodTrack->PAtDCA();
    Int_t numberOfITS = aodTrack->GetITSNcls(); 
    Int_t numberOfTPC = aodTrack->GetTPCNcls();
    Double_t nSigmaTPC = 0;
    Double_t nSigmaTOF = 0;
    Int_t pionPIDnumber = 2;
    Int_t kaonPIDnumber = 3;
    Int_t TPCok = 0;
    Int_t TOFok = 0;

    AliAODPidHF* trackPIDHF = (AliAODPidHF*)fCuts->GetPidHF();
    TPCok = trackPIDHF->GetnSigmaTPC(aodTrack, pionPIDnumber, nSigmaTPC);
    TOFok = trackPIDHF->GetnSigmaTOF(aodTrack, pionPIDnumber, nSigmaTOF);

    Int_t daughterType = 3;
    Int_t histType;

    histType = 0;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(aodTrack->Pt(),aodTrack->DCA());


    //we check if the particle is a signal track
    Bool_t isDesiredCandidate = kFALSE;
    if(fUseMCInfo){
      TMatrix &particleMatrix = *B0toDStarPiLabelMatrix;
      for (Int_t k = 0; k < B0toDStarPiLabelMatrix->GetNrows(); ++k){
        if(mcLabelParticle == (Int_t)particleMatrix(k,0)){
          isDesiredCandidate = kTRUE;
          break;
        }
      }
    }

    if(isDesiredCandidate)
    {
      histType = 2;
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
      if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
      if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
      if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
      ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(aodTrack->Pt(),aodTrack->DCA());
    }

    //we apply a number of cuts on the particle
    Bool_t bCut = kFALSE;

    //we apply a PID cut, 2 is used to indicate we look for a pion
    if(!(fCuts->SelectPID(aodTrack,2))){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[3][1])->Fill(2);
      } else ((TH1F*)fDaughterHistogramArray2D[3][0])->Fill(2);
      bCut = kTRUE;
    }

    if(aodTrack->GetITSNcls() < fCuts->GetMinITSNclsB0Pion()){
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[3][1])->Fill(3);
      } else ((TH1F*)fDaughterHistogramArray2D[3][0])->Fill(3);
      bCut = kTRUE;
    }

    if(aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsB0Pion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[3][1])->Fill(4);
      } else ((TH1F*)fDaughterHistogramArray2D[3][0])->Fill(4);
      bCut = kTRUE;
    }

    if(fCuts->UseITSRefitB0Pion()==kTRUE){
      if(!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArray2D[3][1])->Fill(5);
        } else ((TH1F*)fDaughterHistogramArray2D[3][0])->Fill(5);
        bCut = kTRUE;
      }
    }

    if(fCuts->UseTPCRefitB0Pion()==kTRUE){
      if((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArray2D[3][1])->Fill(6);
        } else ((TH1F*)fDaughterHistogramArray2D[3][0])->Fill(6);
        bCut = kTRUE;
      }
    }

    if(fCuts->UseFilterBitB0Pion()==kTRUE){
      if(!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitB0Pion())))) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArray2D[3][1])->Fill(7);
        } else ((TH1F*)fDaughterHistogramArray2D[3][0])->Fill(7);
        bCut = kTRUE;
      }
    }

    if(aodTrack->Pt() < fCuts->GetMinPtB0Pion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[3][1])->Fill(8);
      } else ((TH1F*)fDaughterHistogramArray2D[3][0])->Fill(8);
      bCut = kTRUE;
    }

    if(!isDesiredCandidate && fQuickSignalAnalysis) bCut = kTRUE;

    if(bCut) {
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArray2D[3][1])->Fill(0);
      } else ((TH1F*)fDaughterHistogramArray2D[3][0])->Fill(0);
      continue;
    }

    //we fill histograms with track information of the tracks that pass the cuts
    histType = 1;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(aodTrack->Pt(),aodTrack->DCA());

    fB0PionTracks->push_back(i);
    numberofparticlesused++;
  }

  ((TH1F*)fDaughterHistogramArray[3][0][12])->Fill(numberofparticles);
  ((TH1F*)fDaughterHistogramArray[3][1][12])->Fill(numberofparticlesused);
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::D0Selection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz,TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix){
  
  
  Int_t iClonesArray = 0;

  TString fillthis = "";

  //we loop over all the first daughter candidates
  for (Int_t i = 0; i < (Int_t)fD0PionTracks->size(); i++){
    
    //we get the track of the first daughter
    AliAODTrack * trackFirstDaughter = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack((fD0PionTracks->at(i))));
    if(!trackFirstDaughter) continue;
    
    //next we loop over all the second daughter candidates
    for (Int_t j = 0; j < (Int_t)fD0KaonTracks->size(); j++){
        
      //we get the track of the second daughter
      AliAODTrack * trackSecondDaughter = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack((fD0KaonTracks->at(j))));
      if(!trackSecondDaughter) continue;
      
      //we check if the IDs of the tracks are different
      if(trackFirstDaughter->GetID() == trackSecondDaughter->GetID() ) continue;

      //we check if the charges of the tracks are correct
      if(trackFirstDaughter->Charge() == trackSecondDaughter->Charge() || trackFirstDaughter->Charge() + trackSecondDaughter->Charge() != 0) continue;

      //we get the impact parameter, distance of closes approach, and the covariance.
      AliExternalTrackParam firstTrack;
      firstTrack.CopyFromVTrack(trackFirstDaughter);
      AliExternalTrackParam secondTrack;
      secondTrack.CopyFromVTrack(trackSecondDaughter);

      //we calculate the vertex of the mother candidate
      TObjArray daughterTracks;              
      daughterTracks.Add(&firstTrack);
      daughterTracks.Add(&secondTrack); 
      AliAODVertex *vertexMother = RecalculateVertex(primaryVertex,&daughterTracks,bz);
      if(!vertexMother) continue;
      if(vertexMother==0) {delete vertexMother; vertexMother = NULL; continue;}

      Double_t xdummy=0.,ydummy=0.,dca,e[2];
      Double_t d0z0[2],covd0z0[3],d0[2],d0err[2];

      firstTrack.PropagateToDCA(vertexMother,bz,100.,d0z0,covd0z0);
      secondTrack.PropagateToDCA(vertexMother,bz,100.,d0z0,covd0z0);

      //we reconstruct the mother decay prong
      Double_t px[2],py[2],pz[2];
      px[0] = firstTrack.Px();
      py[0] = firstTrack.Py();
      pz[0] = firstTrack.Pz();
      px[1] = secondTrack.Px();
      py[1] = secondTrack.Py();
      pz[1] = secondTrack.Pz();
      
      UInt_t prongs[2];
      prongs[0] = 211;
      prongs[1] = 321;
      
      UShort_t id[2];
      id[0]= firstTrack.GetID(); 
      id[1]= secondTrack.GetID();

      firstTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
      d0[0] = d0z0[0];
      d0err[0] = TMath::Sqrt(covd0z0[0]);
      secondTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
      d0[1] = d0z0[0];
      d0err[1] = TMath::Sqrt(covd0z0[0]);

      dca = secondTrack.GetDCA(&firstTrack,bz,xdummy,ydummy);   

      //Short_t chargeMother = trackFirstDaughter->Charge() + trackSecondDaughter->Charge();
      AliAODRecoDecayHF2Prong * motherRecoDecayHF2Prong = new AliAODRecoDecayHF2Prong(vertexMother,px,py,pz,d0,d0err,dca); 
      if(!motherRecoDecayHF2Prong) {delete vertexMother; vertexMother = NULL; continue;}
      
      motherRecoDecayHF2Prong->GetSecondaryVtx()->AddDaughter(trackFirstDaughter);
      motherRecoDecayHF2Prong->GetSecondaryVtx()->AddDaughter(trackSecondDaughter);
      motherRecoDecayHF2Prong->SetPrimaryVtxRef((AliAODVertex*)aodEvent->GetPrimaryVertex());
      motherRecoDecayHF2Prong->SetProngIDs(2,id);

      //we save the pdgcode of the used particle and its mother and check if it is a desired candidate
      Float_t pdgCodeMother = -1;
      Float_t pdgCodeGrandMother = -1;
      Bool_t isDesiredCandidate = kFALSE;
      Int_t motherType, histType;
      motherType = 0;

      if(fUseMCInfo){
        Int_t mcLabelMother = -1;
        Int_t mcLabelGrandMother = -1;
        Int_t mcLabelFirstMother = -1;
        Int_t mcLabelSecondMother = -1;

        Int_t mcLabelFirstTrack = trackFirstDaughter->GetLabel();
        Int_t mcLabelSecondTrack = trackSecondDaughter->GetLabel();

        if(mcLabelFirstTrack >= 0 && mcLabelSecondTrack >= 0){
          AliAODMCParticle *mcParticleFirstTrack = (AliAODMCParticle*)mcTrackArray->At(mcLabelFirstTrack);
          AliAODMCParticle *mcParticleSecondTrack = (AliAODMCParticle*)mcTrackArray->At(mcLabelSecondTrack);

          if(mcParticleFirstTrack && mcParticleSecondTrack){
            mcLabelFirstMother = mcParticleFirstTrack->GetMother();
            mcLabelSecondMother = mcParticleSecondTrack->GetMother();       
            if(mcLabelFirstMother == mcLabelSecondMother) mcLabelMother = mcLabelFirstMother;
          }

          if(mcLabelMother >= 0){

            AliAODMCParticle *mcMotherParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelMother);
            pdgCodeMother = TMath::Abs(mcMotherParticle->GetPdgCode());

            ((TH1F*)fMotherHistogramArray2D[motherType][2])->Fill(pdgCodeMother);
            mcLabelGrandMother = mcMotherParticle->GetMother();

            if(mcLabelGrandMother >= 0){
              AliAODMCParticle *mcGrandMotherParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelGrandMother);
              pdgCodeGrandMother = TMath::Abs(mcGrandMotherParticle->GetPdgCode());
             ((TH1F*)fMotherHistogramArray2D[motherType][3])->Fill(pdgCodeGrandMother);
            }

            //check if the reconstructed mother is a desired candidate
            TMatrix &particleMatrix = *B0toDStarPiLabelMatrix;
            for (Int_t k = 0; k < B0toDStarPiLabelMatrix->GetNrows(); ++k){
              if(mcLabelMother == (Int_t)particleMatrix(k,4) && mcLabelFirstTrack == (Int_t)particleMatrix(k,2) && mcLabelSecondTrack == (Int_t)particleMatrix(k,3)){
                
                isDesiredCandidate = kTRUE;
                break;
              }
            }
          }
        }
      }

      // We fill the histograms
      histType = 0;
      FillD0Histograms(motherRecoDecayHF2Prong, primaryVertex, bz, motherType, histType);
      if(isDesiredCandidate && fUseMCInfo){
        histType = 1;
        FillD0Histograms(motherRecoDecayHF2Prong, primaryVertex, bz, motherType, histType);
      }

      // Here we apply cuts on the particle 
      Bool_t cutMother = kFALSE;

      Bool_t bCutArray[25] = {0};
      Int_t cutReturnValue = fCuts->IsD0forD0ptbinSelected(motherRecoDecayHF2Prong, 0, aodEvent, bCutArray);
      if(cutReturnValue == -1) cutMother = kTRUE;
      if(cutReturnValue == 0) cutMother = kTRUE;


      if(fGetCutInfo == kTRUE)
      {
        for (Int_t k = 0; k < 25; ++k)
        {
          if (bCutArray[k] == kTRUE){
            if(isDesiredCandidate){
              ((TH1F*)fMotherHistogramArray2D[motherType][1])->Fill(k+1);
            } else ((TH1F*)fMotherHistogramArray2D[motherType][0])->Fill(k+1);
            cutMother = kTRUE;
          }
        } 
      }


      if(!isDesiredCandidate && fQuickSignalAnalysis) cutMother = kTRUE;

      if(cutMother){
        if(isDesiredCandidate){
          ((TH1F*)fMotherHistogramArray2D[motherType][1])->Fill(0);
        } else ((TH1F*)fMotherHistogramArray2D[motherType][0])->Fill(0);
        delete vertexMother; vertexMother = NULL; 
        delete motherRecoDecayHF2Prong; motherRecoDecayHF2Prong = NULL;
        continue;
      }

      // We fill the cut histograms
      histType = 2;
      FillD0Histograms(motherRecoDecayHF2Prong, primaryVertex, bz, motherType, histType);
      if(isDesiredCandidate && fUseMCInfo){
        histType = 3;
        FillD0Histograms(motherRecoDecayHF2Prong, primaryVertex, bz, motherType, histType);
      }

      //we save the mother to a TClonesArray
      new ((*fD0Tracks)[iClonesArray++]) AliAODRecoDecayHF2Prong(*motherRecoDecayHF2Prong);
      delete motherRecoDecayHF2Prong; motherRecoDecayHF2Prong = NULL;
    }
  }
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::DStarAndB0Selection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz,TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix){
  
  Int_t iClonesArray = 0;

  TString fillthis = "";

  //we loop over all the DStar pion candidates
  for (Int_t i = 0; i < (Int_t)fDStarPionTracks->size(); i++)
  {

    //we get the track of the DStar pion
    AliAODTrack * trackFirstDaughter = (AliAODTrack*)(aodEvent->GetTrack(fDStarPionTracks->at(i)));
    if(!trackFirstDaughter) continue;

    //next we loop over all the D0 candidates
    for (Int_t j = 0; j < fD0Tracks->GetEntriesFast(); j++)
    {

      //we get the track of the D0
      AliAODRecoDecayHF2Prong * trackSecondDaughter = (AliAODRecoDecayHF2Prong*)(fD0Tracks->At(j));
      if(!trackSecondDaughter) {std::cout << "found none" << std::endl; continue;}
      if(trackSecondDaughter == NULL) {std::cout << "found NULL" << std::endl; continue;}

      //we check if the IDs of the tracks are different
      if(trackFirstDaughter->GetID() == trackSecondDaughter->GetProngID(0) || trackFirstDaughter->GetID() == trackSecondDaughter->GetProngID(1)) continue;

      //we check if the charges of the tracks are correct
      if(trackFirstDaughter->Charge() == trackSecondDaughter->Charge() || TMath::Abs(trackFirstDaughter->Charge() + trackSecondDaughter->Charge()) != 1) continue;
      
      //we check if the pions have the same charge
      if(trackFirstDaughter->Charge() != ((AliAODTrack*)trackSecondDaughter->GetDaughter(0))->Charge()) continue;

      //we loop over all the B0 pion candidates
      for (Int_t k = 0; k < (Int_t)fB0PionTracks->size(); k++)
      {

        //we get the track of the first daughter
        AliAODTrack * trackB0Pion = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(fB0PionTracks->at(k)));
        if(!trackB0Pion) continue;

        // std::cout << "DStar Pion - " << i << "/" << fDStarPionTracks->size() << ", " << "D0 - " << j << "/" << fD0Tracks->GetEntriesFast() << ", " << "D0 Pion - " << k << "/" << fB0PionTracks->size() << ". " << std::endl;

        //we check if the IDs of the tracks are different
        AliAODTrack* twoProngdaughter0 = (AliAODTrack*)trackSecondDaughter->GetDaughter(0);
        AliAODTrack* twoProngdaughter1 = (AliAODTrack*)trackSecondDaughter->GetDaughter(1);
        UShort_t idProng0 = twoProngdaughter0->GetID();
        UShort_t idProng1 = twoProngdaughter1->GetID();
        
        if(trackB0Pion->GetID() == trackFirstDaughter->GetID() || trackB0Pion->GetID() == idProng0 || trackB0Pion->GetID() == idProng1) continue;
        
        //we check if the charges of the tracks are correct // later change this for like sign analysis.
        Bool_t bSameSign = kFALSE;
        if(trackB0Pion->Charge() == (trackSecondDaughter->Charge() + trackFirstDaughter->Charge()) && trackB0Pion->Charge() + (trackSecondDaughter->Charge() + trackFirstDaughter->Charge()) != 0)  bSameSign = kTRUE;

        //we use the DStar pion, B0 pion, and D0 tracks to reconstruct the vertex for the B0 and DStar decay
        AliExternalTrackParam firstTrack;
        firstTrack.CopyFromVTrack(trackFirstDaughter);
        AliExternalTrackParam secondTrack;
        secondTrack.CopyFromVTrack(trackSecondDaughter);
        AliExternalTrackParam thirdTrack;
        thirdTrack.CopyFromVTrack(trackB0Pion);

        // we calculate the vertex of the mother candidate
        // for the recalculated vertex we use the DStar pion and the D0 to calculate the vertex position instead of the DStar 
        TObjArray daughterTracksWithRecalculation;

        daughterTracksWithRecalculation.Add(&firstTrack);
        daughterTracksWithRecalculation.Add(&secondTrack);
        daughterTracksWithRecalculation.Add(&thirdTrack);

        AliAODVertex *vertexMother = RecalculateVertex(primaryVertex,&daughterTracksWithRecalculation,bz);
        if(!vertexMother) {delete vertexMother; vertexMother = NULL; continue;}

        Double_t xdummyDStar=0.,ydummyDStar=0.,dcaDStar,eDStar[2];
        Double_t d0z0DStar[2],covd0z0DStar[3],d0DStar[2],d0errDStar[2];

        //DStar creation with the new vertex
        firstTrack.PropagateToDCA(vertexMother,bz,100.,d0z0DStar,covd0z0DStar);
        secondTrack.PropagateToDCA(vertexMother,bz,100.,d0z0DStar,covd0z0DStar);

        // dcaDStar = secondTrack.GetDCA(&firstTrack,bz,xdummyDStar,ydummyDStar);

        Double_t pxDStar[2],pyDStar[2],pzDStar[2];
        pxDStar[0] = firstTrack.Px();
        pyDStar[0] = firstTrack.Py();
        pzDStar[0] = firstTrack.Pz();
        pxDStar[1] = secondTrack.Px();
        pyDStar[1] = secondTrack.Py();
        pzDStar[1] = secondTrack.Pz();

        Double_t xyz_track1[3];
        xyz_track1[0] = firstTrack.GetX();
        firstTrack.GetYAt(xyz_track1[0],bz,xyz_track1[1]);
        firstTrack.GetZAt(xyz_track1[0],bz,xyz_track1[2]);

        Double_t xyz_track2[3];
        xyz_track2[0] = secondTrack.GetX();
        secondTrack.GetYAt(xyz_track2[0],bz,xyz_track2[1]);
        secondTrack.GetZAt(xyz_track2[0],bz,xyz_track2[2]);

        Double_t distanceAtVertex = TMath::Sqrt((xyz_track1[0]-xyz_track2[0])*(xyz_track1[0]-xyz_track2[0]) + (xyz_track1[1]-xyz_track2[1])*(xyz_track1[1]-xyz_track2[1]) + (xyz_track1[2]-xyz_track2[2])*(xyz_track1[2]-xyz_track2[2]));

        firstTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0DStar,covd0z0DStar);
        d0DStar[0] = d0z0DStar[0];
        d0errDStar[0] = TMath::Sqrt(covd0z0DStar[0]);
        secondTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0DStar,covd0z0DStar);
        d0DStar[1] = d0z0DStar[0];
        d0errDStar[1] = TMath::Sqrt(covd0z0DStar[0]); 


        dcaDStar = secondTrack.GetDCA(&firstTrack,bz,xdummyDStar,ydummyDStar);   


        Short_t chargeDStar = trackFirstDaughter->Charge() + trackSecondDaughter->Charge();
        AliAODVertex * vertexDStar = new AliAODVertex(*vertexMother);
        AliAODRecoCascadeHF * trackDStar = new AliAODRecoCascadeHF(vertexDStar,chargeDStar,pxDStar,pyDStar,pzDStar,d0DStar,d0errDStar,distanceAtVertex); 
        if(!trackDStar) 
        {
          delete vertexMother; vertexMother = NULL; 
          delete vertexDStar; vertexDStar = NULL; 
          delete trackDStar; trackDStar = NULL; 
          continue;
        }
        
        UShort_t idDStar[2];
        idDStar[0]= trackFirstDaughter->GetID(); 
        idDStar[1]= 0;

        UInt_t prongsDStar[2];
        prongsDStar[0] = 211;
        prongsDStar[1] = 421;

        trackDStar->GetSecondaryVtx()->AddDaughter(trackFirstDaughter);
        trackDStar->GetSecondaryVtx()->AddDaughter(trackSecondDaughter);
        trackDStar->SetPrimaryVtxRef((AliAODVertex*)aodEvent->GetPrimaryVertex());
        trackDStar->SetProngIDs(2,idDStar);
        

        Bool_t isDesiredCandidate = kFALSE;
        Int_t mcLabelB0 = 0;
        fillthis = "";
        Int_t motherType, histType;
        motherType = 1;

        if(fUseMCInfo)
        {
          Int_t pdgDgD0toKpi[2]={321,211};
          Int_t pdgDgDStartoD0Pi[2]={421,211};
          Int_t mcLabelFirstDaughter = trackB0Pion->GetLabel();
          Int_t mcLabelSecondDaughter = trackDStar->MatchToMC(413,421,pdgDgDStartoD0Pi,pdgDgD0toKpi,mcTrackArray);
          Int_t mcLabelDStarPion = ((AliAODTrack*)trackDStar->GetDaughter(0))->GetLabel();
          if(mcLabelFirstDaughter>=0 && mcLabelSecondDaughter>=0)
          {

            Int_t mcLabelMotherFirstTrack = -1;
            Int_t mcLabelMotherSecondTrack = -1;
            AliAODMCParticle *mcTrackFirstDaughter = (AliAODMCParticle*)mcTrackArray->At(mcLabelFirstDaughter);
            AliAODMCParticle *mcTrackSecondDaughter = (AliAODMCParticle*)mcTrackArray->At(mcLabelSecondDaughter);
            AliAODMCParticle *mcTrackDStarPion = (AliAODMCParticle*)mcTrackArray->At(mcLabelDStarPion);

            if(mcTrackFirstDaughter) mcLabelMotherFirstTrack = mcTrackFirstDaughter->GetMother();
            if(mcTrackSecondDaughter) mcLabelMotherSecondTrack = mcTrackSecondDaughter->GetMother();

            if(mcLabelMotherFirstTrack>=0 && mcLabelMotherFirstTrack == mcLabelMotherSecondTrack)
            {
              TMatrix &particleMatrix = *B0toDStarPiLabelMatrix;
              for (Int_t k = 0; k < B0toDStarPiLabelMatrix->GetNrows(); ++k)
              {
                if(mcLabelMotherFirstTrack == (Int_t)particleMatrix(k,6)  && mcLabelSecondDaughter == (Int_t)particleMatrix(k,5) && mcLabelDStarPion == (Int_t)particleMatrix(k,1) && mcLabelFirstDaughter == (Int_t)particleMatrix(k,0))
                {

                  isDesiredCandidate = kTRUE;
                  mcLabelB0 = mcLabelMotherFirstTrack;
                  break;
                }
              }
            }
          }
        }
        

        // We fill the histograms

        histType = 0;
        FillCascadeMotherHistograms(trackDStar, primaryVertex, bz, motherType, histType);
        if(isDesiredCandidate && fUseMCInfo)
        {
          histType = 1;
          FillCascadeMotherHistograms(trackDStar, primaryVertex, bz, motherType, histType);
        }
     

        // We apply cuts 
        Bool_t cutDStar = kFALSE;

        Bool_t bCutArrayDStar[25] = {0};
        Int_t cutReturnValueDStar = fCuts->IsDStarforDStarptbinSelected(trackDStar, 0, aodEvent, bCutArrayDStar);
        if(cutReturnValueDStar == -1) cutDStar = kTRUE;
        if(cutReturnValueDStar == 0) cutDStar = kTRUE;

        Bool_t bCutArrayD0[35] = {0};
        Int_t cutReturnValueD0 = fCuts->IsD0forDStarptbinSelected(trackDStar, 0, aodEvent, bCutArrayD0);
        if(cutReturnValueD0 == -1) cutDStar = kTRUE;
        if(cutReturnValueD0 == 0) cutDStar = kTRUE;


        if(fGetCutInfo == kTRUE)
        {

          for (Int_t n = 0; n < 25; ++n)
          {
            if(bCutArrayDStar[n] == kTRUE){
              if(isDesiredCandidate){
                ((TH1F*)fMotherHistogramArray2D[motherType][1])->Fill(n+1);
              } else ((TH1F*)fMotherHistogramArray2D[motherType][0])->Fill(n+1);
              cutDStar = kTRUE;
            }
          }

          for (Int_t n = 0; n < 35; ++n)
          {
            if(bCutArrayD0[n] == kTRUE){
              if(isDesiredCandidate){
                ((TH1F*)fMotherHistogramArray2D[motherType][1])->Fill(n+1+35);
              } else ((TH1F*)fMotherHistogramArray2D[motherType][0])->Fill(n+1+35);
              cutDStar = kTRUE;
            }
          }
        }

        if(!isDesiredCandidate && fQuickSignalAnalysis) cutDStar = kFALSE;

        if(cutDStar)
        {
          if(isDesiredCandidate)
          {
            ((TH1F*)fMotherHistogramArray2D[motherType][1])->Fill(0);
          } else ((TH1F*)fMotherHistogramArray2D[motherType][0])->Fill(0);
          delete vertexMother; vertexMother = NULL;
          delete vertexDStar; vertexDStar = NULL;  
          delete trackDStar; trackDStar = NULL;
          continue;
        }

        // We fill the cut histograms
        histType = 2;
        FillCascadeMotherHistograms(trackDStar, primaryVertex, bz, motherType, histType);
        if(isDesiredCandidate && fUseMCInfo)
        {
          histType = 3;
          FillCascadeMotherHistograms(trackDStar, primaryVertex, bz, motherType, histType);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////<<
        //
        // BO Reconstruction
        //
        //
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////<<


        //use the new DStar candidate and the new vertex to create the B0 candidate
        Double_t xdummy=0.,ydummy=0.,dca,e[2];
        Double_t d0z0[2],covd0z0[3],d0[2],d0err[2];

        AliExternalTrackParam fourthTrack;
        fourthTrack.CopyFromVTrack(trackDStar);

        thirdTrack.PropagateToDCA(vertexMother,bz,100.,d0z0,covd0z0);
        fourthTrack.PropagateToDCA(vertexMother,bz,100.,d0z0,covd0z0);

        //we reconstruct the mother decay prong
        Double_t px[2],py[2],pz[2];
        px[0] = thirdTrack.Px();
        py[0] = thirdTrack.Py();
        pz[0] = thirdTrack.Pz();
        px[1] = fourthTrack.Px();
        py[1] = fourthTrack.Py();
        pz[1] = fourthTrack.Pz();
        
        UInt_t prongs[2];
        prongs[0] = 211;
        prongs[1] = 413;
        
        UShort_t id[2];
        id[0]= thirdTrack.GetID(); 
        id[1]= 0;

        thirdTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
        d0[0] = d0z0[0];
        d0err[0] = TMath::Sqrt(covd0z0[0]);
        fourthTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
        d0[1] = d0z0[0];
        d0err[1] = TMath::Sqrt(covd0z0[0]);  
        
        dca = fourthTrack.GetDCA(&thirdTrack,bz,xdummy,ydummy);   

        Short_t chargeMother = trackFirstDaughter->Charge() + trackDStar->Charge();
        AliAODRecoCascadeHF * motherCascadeHF = new AliAODRecoCascadeHF(vertexMother,chargeMother,px,py,pz,d0,d0err,dca); 
        if(!motherCascadeHF)
        {
          delete vertexMother; vertexMother = NULL; 
          delete vertexDStar; vertexDStar = NULL; 
          delete trackDStar; trackDStar = NULL; 
          continue;
        }
        
        motherCascadeHF->GetSecondaryVtx()->AddDaughter(trackB0Pion);
        motherCascadeHF->GetSecondaryVtx()->AddDaughter(trackDStar);
        motherCascadeHF->SetPrimaryVtxRef((AliAODVertex*)aodEvent->GetPrimaryVertex());
        motherCascadeHF->SetProngIDs(2,id);

        fillthis = "";

        motherType = 2;
        histType = 0;

        if(!bSameSign)
        {
          // We fill the histograms
          histType = 0;
          FillCascadeMotherHistograms(motherCascadeHF, primaryVertex, bz, motherType, histType);
          
          if(isDesiredCandidate)
          {
            histType = 1;
            FillCascadeMotherHistograms(motherCascadeHF, primaryVertex, bz, motherType, histType);
          }
        }


        // We apply cuts
        Bool_t cutMother = kFALSE;

        Bool_t bCutArray[85] = {0};
        Int_t numberOfCuts = 85;
        Int_t cutReturnValue = fCuts->IsSelected(motherCascadeHF, 0, aodEvent, bCutArray);
        if(cutReturnValue == -1) cutMother = kTRUE;
        if(cutReturnValue == 0) cutMother = kTRUE;


        // We save information about the cuts
        TString histName = "";
        Double_t invariantMassMother = motherCascadeHF->InvMass(2,prongs);
        Double_t pdgMassMother=TDatabasePDG::Instance()->GetParticle(511)->Mass();
        Double_t massWindow = 0.125; //75; //GeV/c^2
        if(fGetCutInfo == kTRUE)
        {
          for (Int_t n = 0; n < 85; ++n)
          {
            if(bCutArray[n] == kTRUE){
              if(isDesiredCandidate){
                ((TH1F*)fMotherHistogramArray2D[motherType][1])->Fill(n+1);
              } else ((TH1F*)fMotherHistogramArray2D[motherType][0])->Fill(n+1);
              cutMother = kTRUE;
            }
          }
       
          if (TMath::Abs(invariantMassMother-pdgMassMother)<massWindow){
            for (Int_t i = 0; i < numberOfCuts; ++i) //total
            {
              if(bCutArray[i] == kFALSE) continue;
              for (Int_t j = 0; j < numberOfCuts; ++j)
              {
                if(bCutArray[j] == kFALSE) continue;
                if(isDesiredCandidate == kFALSE) histName ="cutEffectBackground";
                if(isDesiredCandidate == kTRUE) histName ="cutEffectSignal";
                ((TH2I*)(fOutputB0MC->FindObject(histName)))->Fill(i,j);
              }
            }

            for (Int_t i = 0; i < numberOfCuts; ++i) //unique
            {
              if(bCutArray[i] == kFALSE) continue;
              Bool_t bFill = kTRUE;
              for (Int_t j = 0; j < numberOfCuts; ++j)
              {
                if(i==j) continue;
                if(bCutArray[j] == kTRUE) 
                {
                  bFill = kFALSE;
                  break;
                }

              }
              if(bFill == kTRUE)
              {
                if(isDesiredCandidate == kFALSE) histName ="cutEffectUniqueBackground";
                if(isDesiredCandidate == kTRUE) histName ="cutEffectUniqueSignal";
                ((TH1I*)(fOutputB0MC->FindObject(histName)))->Fill(i);
              }
            }
          }
        }


        if(!isDesiredCandidate && fQuickSignalAnalysis) cutMother = kTRUE;

        if(cutMother)
        {
          if(isDesiredCandidate)
          {
            ((TH1F*)fMotherHistogramArray2D[motherType][1])->Fill(0);
          } else ((TH1F*)fMotherHistogramArray2D[motherType][0])->Fill(0);
          delete vertexMother; vertexMother = NULL; 
          delete motherCascadeHF; motherCascadeHF = NULL;
          delete vertexDStar; vertexDStar = NULL; 
          delete trackDStar; trackDStar = NULL; 
          continue;
        }


        if(!bSameSign)
        {
          histType = 2;
          FillCascadeMotherHistograms(motherCascadeHF, primaryVertex, bz, motherType, histType);
          if(fUseMCInfo && isDesiredCandidate)
          {
            //fill mc histograms
            histType = 3;
            FillCascadeMotherHistograms(motherCascadeHF, primaryVertex, bz, motherType, histType);
          }
        }

        AliAODRecoCascadeHF* selectedDStar = (AliAODRecoCascadeHF*)motherCascadeHF->GetDaughter(1);
        AliAODRecoDecayHF2Prong* selectedD0 = (AliAODRecoDecayHF2Prong*)selectedDStar->GetDaughter(1);

        // We fill the final cut histograms with the candidates that have an invariant mass close to the PDG value. This way the background we use for optimizing the cuts will not be contaminated with candidates that don't affect the signal region.
        // pdgMassMother=TDatabasePDG::Instance()->GetParticle(511)->Mass();
        // massWindow = 0.125; //75; //GeV/c^2
        if (TMath::Abs(invariantMassMother-pdgMassMother)<massWindow)
        {
          if(!bSameSign) 
          {
            FillFinalTrackHistograms(motherCascadeHF, isDesiredCandidate, mcTrackArray);
            if(!isDesiredCandidate)
            {
              motherType = 0; histType = 4; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histType);
              motherType = 1; histType = 4; FillCascadeMotherHistograms(selectedDStar, primaryVertex, bz, motherType, histType);
              motherType = 2; histType = 4; FillCascadeMotherHistograms(motherCascadeHF, primaryVertex, bz, motherType, histType);
            }
            if(isDesiredCandidate)
            {
              motherType = 0; histType = 5; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histType);
              motherType = 1; histType = 5; FillCascadeMotherHistograms(selectedDStar, primaryVertex, bz, motherType, histType);
              motherType = 2; histType = 5; FillCascadeMotherHistograms(motherCascadeHF, primaryVertex, bz, motherType, histType);
            }
          }
        }

        //Here we fill the histograms per pt bin and apply the same sign method
        TString ptBinMother = "";
        Int_t ptBin = fCuts->PtBin(motherCascadeHF->Pt());
        ptBinMother += "_ptbin_"; ptBinMother += fPtBinLimits[ptBin]; ptBinMother += "_to_"; ptBinMother += fPtBinLimits[ptBin+1];
        histType = 6 + 2 * ptBin; 

        Int_t d0PtBin = fCuts->PtBinD0forD0ptbin(selectedD0->Pt());
        Int_t histTypeD0 = 2 * d0PtBin; 

        Int_t d0DStarPtBin = fCuts->PtBinD0forDStarptbin(selectedDStar->Pt());
        Int_t histTypeD0DStar = 2 * d0DStarPtBin;

        Int_t dstarPtBin = fCuts->PtBinDStarforDStarptbin(selectedDStar->Pt());
        Int_t histTypeDStar = 2 * dstarPtBin;


        if (TMath::Abs(invariantMassMother-pdgMassMother)<massWindow)
        {
          if(!bSameSign && histType > 5) 
          {
            if(!isDesiredCandidate)
            {
              motherType = 0; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histType);
              motherType = 1; FillCascadeMotherHistograms(selectedDStar, primaryVertex, bz, motherType, histType);
              motherType = 2; FillCascadeMotherHistograms(motherCascadeHF, primaryVertex, bz, motherType, histType);
              motherType = 3; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histTypeD0);
              motherType = 4; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histTypeD0DStar);
              motherType = 5; FillCascadeMotherHistograms(selectedDStar, primaryVertex, bz, motherType, histTypeDStar);
            }

            if(isDesiredCandidate)
            {
              histType += 1;
              motherType = 0; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histType);
              motherType = 1; FillCascadeMotherHistograms(selectedDStar, primaryVertex, bz, motherType, histType);
              motherType = 2; FillCascadeMotherHistograms(motherCascadeHF, primaryVertex, bz, motherType, histType);
              motherType = 3; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histTypeD0 + 1);
              motherType = 4; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histTypeD0DStar + 1);
              motherType = 5; FillCascadeMotherHistograms(selectedDStar, primaryVertex, bz, motherType, histTypeDStar + 1);
            }
          }
        }

        Double_t invmassDelta = DeltaInvMassB0Kpipipi(motherCascadeHF);
        if(bSameSign)
        {
          fillthis="invariantMassB0";
          fillthis += "_SameSign";
          ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          fillthis="invariantMassB0";
          fillthis += ptBinMother + "_SameSign";
          ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          fillthis="invariantMassB0";
          fillthis += "_SignSum";
          ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,-1);
          fillthis="invariantMassB0";
          fillthis += ptBinMother + "_SignSum";
          ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,-1);
        }
        if(!bSameSign)
        {

          fillthis="deltaInvMassB0";
          ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          fillthis="deltaInvMassB0";
          fillthis += ptBinMother;
          ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);

          fillthis="invariantMassB0";
          ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          fillthis="invariantMassB0";
          fillthis += ptBinMother;
          ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          fillthis="invariantMassB0";
          fillthis += "_SignSum";
          ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,1);
          fillthis="invariantMassB0";
          fillthis += ptBinMother + "_SignSum";
          ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,1);
        }

        if(motherCascadeHF->Pt() > 6.0)
        {
          ptBinMother = "_ptbin_6_to_inf";
          if(bSameSign)
          {
            fillthis="invariantMassB0";
            fillthis += ptBinMother + "_SameSign";
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
            fillthis="invariantMassB0";
            fillthis += ptBinMother + "_SignSum";
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,-1);
          }
          if(!bSameSign)
          {
            fillthis="deltaInvMassB0";
            fillthis += ptBinMother;
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);

            fillthis="invariantMassB0";
            fillthis += ptBinMother;
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
            fillthis="invariantMassB0";
            fillthis += ptBinMother + "_SignSum";
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,1);
          }

        }

        if(motherCascadeHF->Pt() > 3.0)
        {
          ptBinMother = "_ptbin_3_to_inf";
          if(bSameSign)
          {
            fillthis="invariantMassB0";
            fillthis += ptBinMother + "_SameSign";
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
            fillthis="invariantMassB0";
            fillthis += ptBinMother + "_SignSum";
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,-1);
          }
          if(!bSameSign)
          {
            fillthis="deltaInvMassB0";
            fillthis += ptBinMother;
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);

            fillthis="invariantMassB0";
            fillthis += ptBinMother;
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
            fillthis="invariantMassB0";
            fillthis += ptBinMother + "_SignSum";
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,1);
          }
          
        }

        //we save the mother to a TClonesArray
        //if(!bSameSign) new ((*fB0Tracks)[iClonesArray++]) AliAODRecoCascadeHF(*motherCascadeHF); 
        delete vertexMother; vertexMother = NULL; 
        delete motherCascadeHF; motherCascadeHF = NULL;
        delete vertexDStar; vertexDStar = NULL; 
        delete trackDStar; trackDStar = NULL;
      }  
    }
  }
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::FillFinalTrackHistograms(AliAODRecoCascadeHF * selectedB0, Bool_t isDesiredCandidate,TClonesArray * mcTrackArray){

  //In this function we fill histograms with the properties of all the daughters of our selected signal candidate

  AliAODTrack* selectedB0Pion = (AliAODTrack*)selectedB0->GetDaughter(0);
  AliAODRecoCascadeHF* selectedDStar = (AliAODRecoCascadeHF*)selectedB0->GetDaughter(1);

  AliAODTrack* selectedDStarPion = (AliAODTrack*)selectedDStar->GetDaughter(0);
  AliAODRecoDecayHF2Prong* selectedD0 = (AliAODRecoDecayHF2Prong*)selectedDStar->GetDaughter(1);

  AliAODTrack* selectedD0Pion = (AliAODTrack*)selectedD0->GetDaughter(0);
  AliAODTrack* selectedD0Kaon = (AliAODTrack*)selectedD0->GetDaughter(1);

  Double_t pt_track = 0;
  Double_t momentum_track = 0;
  Double_t energy_track = 0;
  Double_t dca_track = 0;
  Double_t momentum_dca_track = 0;
  Int_t numberOfITS = 0;
  Int_t numberOfTPC = 0;
  Int_t daughterType, histType;
  Int_t totalNumberOfITS = 0;
  Int_t totalNumberOfTPC = 0;
  Double_t nSigmaTPC = 0;
  Double_t nSigmaTOF = 0;
  Double_t nSigmaTPCtotal = 0;
  Double_t nSigmaTOFtotal = 0;
  Int_t pionPIDnumber = 2;
  Int_t kaonPIDnumber = 3;
  Int_t TPCok = 0;
  Int_t TOFok = 0;

  AliAODPidHF* trackPIDHF = (AliAODPidHF*)fCuts->GetPidHF();

  //fill the D0 pion info
  pt_track = selectedD0Pion->Pt();
  momentum_track = selectedD0Pion->P();
  energy_track = selectedD0Pion->E(0.13957);
  dca_track = selectedD0Pion->DCA();
  momentum_dca_track = selectedD0Pion->PAtDCA();
  numberOfITS = selectedD0Pion->GetITSNcls(); 
  numberOfTPC = selectedD0Pion->GetTPCNcls();
  totalNumberOfITS += numberOfITS;
  totalNumberOfTPC += numberOfTPC;
  TPCok = trackPIDHF->GetnSigmaTPC(selectedD0Pion, pionPIDnumber, nSigmaTPC);
  TOFok = trackPIDHF->GetnSigmaTOF(selectedD0Pion, pionPIDnumber, nSigmaTOF);
  if(TPCok != -1) nSigmaTPCtotal += nSigmaTPC*nSigmaTPC;
  if(TOFok != -1) nSigmaTOFtotal += nSigmaTOF*nSigmaTOF;

  Double_t ptB0 = selectedB0->Pt();

  daughterType = 0;
  histType = 3;
  if(!isDesiredCandidate)
  {
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    for (Int_t i = 0; i < 6; ++i)
    {
      if(selectedD0Pion->HasPointOnITSLayer(i)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][7])->Fill(i);
    }
    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));

    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(pt_track,dca_track);
    ((TH2F*)fDaughterHistogramArray2D[daughterType][4])->Fill(ptB0,pt_track);
  }


  if(fUseMCInfo && isDesiredCandidate){
    histType = 4;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    for (Int_t i = 0; i < 6; ++i)
    {
      if(selectedD0Pion->HasPointOnITSLayer(i)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][7])->Fill(i);
    }
    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(pt_track,dca_track);
    ((TH2F*)fDaughterHistogramArray2D[daughterType][5])->Fill(ptB0,pt_track);
  }

  //we save the pdgcode of the used particle and its mother to check PID efficiency
  if(fUseMCInfo)
  {
    Float_t pdgCodeParticle = -1;
    Float_t pdgCodeParticleMother = -1;
    Int_t mcLabelParticle = -1;
    Int_t mcLabelParticleMother = -1;
    mcLabelParticle = selectedD0Pion->GetLabel();

    if(mcLabelParticle >= 0){

      AliAODMCParticle *mcTrackParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticle);
      pdgCodeParticle = TMath::Abs(mcTrackParticle->GetPdgCode());
      ((TH1F*)fDaughterHistogramArray2D[0][2])->Fill(pdgCodeParticle);
      mcLabelParticleMother = mcTrackParticle->GetMother();

      if(mcLabelParticleMother >= 0){
        AliAODMCParticle *mcTrackParticleMother = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticleMother);
        pdgCodeParticleMother = TMath::Abs(mcTrackParticleMother->GetPdgCode());
       ((TH1F*)fDaughterHistogramArray2D[0][3])->Fill(pdgCodeParticleMother);
      }
    }
  }

  //fill the D0 kaon info
  pt_track = selectedD0Kaon->Pt();
  momentum_track = selectedD0Kaon->P();
  energy_track = selectedD0Kaon->E(0.4937);
  dca_track = selectedD0Kaon->DCA();
  momentum_dca_track = selectedD0Kaon->PAtDCA();
  numberOfITS = selectedD0Kaon->GetITSNcls(); 
  numberOfTPC = selectedD0Kaon->GetTPCNcls();
  totalNumberOfITS += numberOfITS;
  totalNumberOfTPC += numberOfTPC;
  TPCok = trackPIDHF->GetnSigmaTPC(selectedD0Kaon, kaonPIDnumber, nSigmaTPC);
  TOFok = trackPIDHF->GetnSigmaTOF(selectedD0Kaon, kaonPIDnumber, nSigmaTOF);
  if(TPCok != -1) nSigmaTPCtotal += nSigmaTPC*nSigmaTPC;
  if(TOFok != -1) nSigmaTOFtotal += nSigmaTOF*nSigmaTOF;

  daughterType = 1;
  histType = 3;
  if(!isDesiredCandidate)
  {
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    for (Int_t i = 0; i < 6; ++i)
    {
      if(selectedD0Kaon->HasPointOnITSLayer(i)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][7])->Fill(i);
    }
    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(pt_track,dca_track);
    ((TH2F*)fDaughterHistogramArray2D[daughterType][4])->Fill(ptB0,pt_track);
  }

  if(fUseMCInfo && isDesiredCandidate){
    histType = 4;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    for (Int_t i = 0; i < 6; ++i)
    {
      if(selectedD0Kaon->HasPointOnITSLayer(i)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][7])->Fill(i);
    }
    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(pt_track,dca_track);
    ((TH2F*)fDaughterHistogramArray2D[daughterType][5])->Fill(ptB0,pt_track);
  }

  //we save the pdgcode of the used particle and its mother to check PID efficiency
  if(fUseMCInfo)
  {
    Float_t pdgCodeParticle = -1;
    Float_t pdgCodeParticleMother = -1;
    Int_t mcLabelParticle = -1;
    Int_t mcLabelParticleMother = -1;
    mcLabelParticle = selectedD0Kaon->GetLabel();

    if(mcLabelParticle >= 0){

      AliAODMCParticle *mcTrackParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticle);
      pdgCodeParticle = TMath::Abs(mcTrackParticle->GetPdgCode());
      ((TH1F*)fDaughterHistogramArray2D[1][2])->Fill(pdgCodeParticle);
      mcLabelParticleMother = mcTrackParticle->GetMother();

      if(mcLabelParticleMother >= 0){
        AliAODMCParticle *mcTrackParticleMother = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticleMother);
        pdgCodeParticleMother = TMath::Abs(mcTrackParticleMother->GetPdgCode());
       ((TH1F*)fDaughterHistogramArray2D[1][3])->Fill(pdgCodeParticleMother);
      }
    }
  }

  //fill the DStar pion info
  pt_track = selectedDStarPion->Pt();
  momentum_track = selectedDStarPion->P();
  energy_track = selectedDStarPion->E(0.13957);
  dca_track = selectedDStarPion->DCA();
  momentum_dca_track = selectedDStarPion->PAtDCA();
  numberOfITS = selectedDStarPion->GetITSNcls(); 
  numberOfTPC = selectedDStarPion->GetTPCNcls();
  totalNumberOfITS += numberOfITS;
  totalNumberOfTPC += numberOfTPC;
  TPCok = trackPIDHF->GetnSigmaTPC(selectedDStarPion, pionPIDnumber, nSigmaTPC);
  TOFok = trackPIDHF->GetnSigmaTOF(selectedDStarPion, pionPIDnumber, nSigmaTOF);
  if(TPCok != -1) nSigmaTPCtotal += nSigmaTPC*nSigmaTPC;

  daughterType = 2;
  histType = 3;
  if(!isDesiredCandidate)
  {
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    for (Int_t i = 0; i < 6; ++i)
    {
      if(selectedDStarPion->HasPointOnITSLayer(i)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][7])->Fill(i);
    }
    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(pt_track,dca_track);
    ((TH2F*)fDaughterHistogramArray2D[daughterType][4])->Fill(ptB0,pt_track);
  }

  if(fUseMCInfo && isDesiredCandidate){
  histType = 4;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    for (Int_t i = 0; i < 6; ++i)
    {
      if(selectedDStarPion->HasPointOnITSLayer(i)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][7])->Fill(i);
    }
    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(pt_track,dca_track);
    ((TH2F*)fDaughterHistogramArray2D[daughterType][5])->Fill(ptB0,pt_track);
  }

  //we save the pdgcode of the used particle and its mother to check PID efficiency
  if(fUseMCInfo)
  {
    Float_t pdgCodeParticle = -1;
    Float_t pdgCodeParticleMother = -1;
    Int_t mcLabelParticle = -1;
    Int_t mcLabelParticleMother = -1;
    mcLabelParticle = selectedDStarPion->GetLabel();

    if(mcLabelParticle >= 0){

      AliAODMCParticle *mcTrackParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticle);
      pdgCodeParticle = TMath::Abs(mcTrackParticle->GetPdgCode());
      ((TH1F*)fDaughterHistogramArray2D[2][2])->Fill(pdgCodeParticle);
      mcLabelParticleMother = mcTrackParticle->GetMother();

      if(mcLabelParticleMother >= 0){
        AliAODMCParticle *mcTrackParticleMother = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticleMother);
        pdgCodeParticleMother = TMath::Abs(mcTrackParticleMother->GetPdgCode());
       ((TH1F*)fDaughterHistogramArray2D[2][3])->Fill(pdgCodeParticleMother);
      }
    }
  }

  //fill the B0 pion info
  pt_track = selectedB0Pion->Pt();
  momentum_track = selectedB0Pion->P();
  energy_track = selectedB0Pion->E(0.13957);
  dca_track = selectedB0Pion->DCA();
  momentum_dca_track = selectedB0Pion->PAtDCA();
  numberOfITS = selectedB0Pion->GetITSNcls(); 
  numberOfTPC = selectedB0Pion->GetTPCNcls();
  totalNumberOfITS += numberOfITS;
  totalNumberOfTPC += numberOfTPC;
  TPCok = trackPIDHF->GetnSigmaTPC(selectedB0Pion, pionPIDnumber, nSigmaTPC);
  TOFok = trackPIDHF->GetnSigmaTOF(selectedB0Pion, pionPIDnumber, nSigmaTOF);
  if(TPCok != -1) nSigmaTPCtotal += nSigmaTPC*nSigmaTPC;
  if(TOFok != -1) nSigmaTOFtotal += nSigmaTOF*nSigmaTOF;

  daughterType = 3;
  histType = 3;
  if(!isDesiredCandidate)
  {
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    for (Int_t i = 0; i < 6; ++i)
    {
      if(selectedB0Pion->HasPointOnITSLayer(i)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][7])->Fill(i);
    }
    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(pt_track,dca_track);
    ((TH2F*)fDaughterHistogramArray2D[daughterType][4])->Fill(ptB0,pt_track);
  }

  if(fUseMCInfo && isDesiredCandidate){
    histType = 4;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(energy_track/momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(momentum_dca_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(numberOfTPC);
    for (Int_t i = 0; i < 6; ++i)
    {
      if(selectedB0Pion->HasPointOnITSLayer(i)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][7])->Fill(i);
    }
    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][9])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][10])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH2F*)fDaughterHistogramArray[daughterType][histType][11])->Fill(pt_track,dca_track);
    ((TH2F*)fDaughterHistogramArray2D[daughterType][5])->Fill(ptB0,pt_track);
  }

    //we save the pdgcode of the used particle and its mother to check PID efficiency
  if(fUseMCInfo)
  {
    Float_t pdgCodeParticle = -1;
    Float_t pdgCodeParticleMother = -1;
    Int_t mcLabelParticle = -1;
    Int_t mcLabelParticleMother = -1;
    mcLabelParticle = selectedB0Pion->GetLabel();

    if(mcLabelParticle >= 0){

      AliAODMCParticle *mcTrackParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticle);
      pdgCodeParticle = TMath::Abs(mcTrackParticle->GetPdgCode());
      ((TH1F*)fDaughterHistogramArray2D[3][2])->Fill(pdgCodeParticle);
      mcLabelParticleMother = mcTrackParticle->GetMother();

      if(mcLabelParticleMother >= 0){
        AliAODMCParticle *mcTrackParticleMother = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticleMother);
        pdgCodeParticleMother = TMath::Abs(mcTrackParticleMother->GetPdgCode());
       ((TH1F*)fDaughterHistogramArray2D[3][3])->Fill(pdgCodeParticleMother);
      }
    }
  }

  if(!isDesiredCandidate)
  {
    ((TH1F*)(fOutputB0MC->FindObject("totalITSBackground")))->Fill(totalNumberOfITS);
    ((TH1F*)(fOutputB0MC->FindObject("totalTPCBackground")))->Fill(totalNumberOfTPC);
    ((TH1F*)(fOutputB0MC->FindObject("totalSigmaPIDBackground")))->Fill(sqrt(nSigmaTPCtotal + nSigmaTOFtotal));
  }
  if(isDesiredCandidate)
  {
    ((TH1F*)(fOutputB0MC->FindObject("totalITSSignal")))->Fill(totalNumberOfITS);
    ((TH1F*)(fOutputB0MC->FindObject("totalTPCSignal")))->Fill(totalNumberOfTPC);
    ((TH1F*)(fOutputB0MC->FindObject("totalSigmaPIDSignal")))->Fill(sqrt(nSigmaTPCtotal + nSigmaTOFtotal));
  }


  return;
}
//-------------------------------------------------------------------------------------
Double_t AliAnalysisTaskSEB0toDStarPi::DeltaInvMassDStarKpipi(AliAODRecoCascadeHF * DStar) const 
{
  ///
  /// 3 prong invariant mass of the D0 daughters and the soft pion
  ///
  Double_t e[3];
  e[0]=DStar->Get2Prong()->EProng(0,211);
  e[1]=DStar->Get2Prong()->EProng(1,321);
  e[2]=DStar->EProng(0,211);

  Double_t esum = e[0]+e[1]+e[2];
  Double_t invMassDStar = TMath::Sqrt(esum*esum-DStar->P()*DStar->P());

  Double_t invMassD0 = DStar->Get2Prong()->InvMassD0();

  return invMassDStar - invMassD0; 
}
//-------------------------------------------------------------------------------------
Double_t AliAnalysisTaskSEB0toDStarPi::DeltaInvMassB0Kpipipi(AliAODRecoCascadeHF * B0) const 
{
  ///
  /// 4 prong invariant mass of the D0 daughters, the soft pion, and the B0 pion
  ///

  AliAODRecoCascadeHF * DStar = (AliAODRecoCascadeHF*)B0->GetDaughter(1);

  Double_t e[4];
  e[0]=DStar->EProng(0,211);
  e[1]=DStar->Get2Prong()->EProng(0,211);
  e[2]=DStar->Get2Prong()->EProng(1,321);
  e[3]=B0->EProng(0,211);

  Double_t esum = e[0]+e[1]+e[2]+e[3];
  Double_t invMassB0 = TMath::Sqrt(esum*esum-B0->P()*B0->P());

  UInt_t prongs[2];
  prongs[0] = 211;
  prongs[1] = 421;

  Double_t invMassD0 = DStar->Get2Prong()->InvMassD0();

  return invMassB0 - invMassD0; 
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::FillD0Histograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType){

  if(histType<0) return;

  //In this function we fill the histograms of the reconstructed mothers
  Double_t ptMother = 0.0;
  Double_t momentumMother = 0.0;
  Double_t etaMother = 0.0;
  Double_t phiMother = 0.0;
  Double_t d0Mother = 0.0;
  Double_t d0firstTrack = 0.0;
  Double_t d0secondTrack = 0.0;
  Double_t pointingAngle = 0.0;
  Double_t impactProduct = 0.0;
  Double_t impactProductXY = 0.0;
  Double_t invariantMassMother = 0.0;
  Double_t invmassDelta = 0.0;
  Double_t dcaMother = 0.0;
  AliAODVertex * vertexMother = 0x0;
  AliAODVertex * vertexDaughter = 0x0;
  Double_t vertexDistance = 0.0;
  Double_t decayLengthDaughter = 0.0;
  Double_t decayTime = 0.0;
  Double_t angleMotherFirstDaughter = 0.0;
  Double_t angleMotherSecondDaughter = 0.0;  
  Double_t topomaticFirstDaughter = 0.0;
  Double_t topomaticSecondDaughter = 0.0;  
  Double_t ptFirstDaughter = 0.0;
  Double_t ptSecondDaughter = 0.0;
  UInt_t prongs[2];
  Double_t angleBetweenBothDaughters = 0;
  Double_t cosThetaStar = 0;
  Double_t normDecayLength = 0;
  Double_t pdgMassMother=TDatabasePDG::Instance()->GetParticle(421)->Mass();


  prongs[0] = 211; prongs[1] = 321;
  AliAODTrack * firstDaughter = (AliAODTrack*)selectedMother->GetDaughter(0);
  AliAODTrack * secondDaughter = (AliAODTrack*)selectedMother->GetDaughter(1);
  vertexMother = selectedMother->GetSecondaryVtx();
  ptFirstDaughter = firstDaughter->Pt();
  ptSecondDaughter = secondDaughter->Pt();

  //Topomatic
  Double_t dd0pr1=0.;
  Double_t dd0pr2=0.;
  Double_t dd0max=0.;
  Double_t dd0min=0.;
  for(Int_t ipr=0; ipr<2; ipr++) 
  {
    Double_t diffIP, errdiffIP;
    selectedMother->Getd0MeasMinusExpProng(ipr,bz,diffIP,errdiffIP);
    Double_t normdd0=0.;
    if(errdiffIP>0.) normdd0=diffIP/errdiffIP;
    if(ipr==0) dd0pr1=normdd0;
    if(ipr==1) dd0pr2=normdd0;

    // else if(TMath::Abs(normdd0)>TMath::Abs(dd0max)) dd0max=normdd0;
  }
  if(TMath::Abs(dd0pr1)>TMath::Abs(dd0pr2)) {dd0max=dd0pr1; dd0min=dd0pr2;}
  else {dd0max=dd0pr2; dd0min=dd0pr1;}

  AliExternalTrackParam motherTrack;
  motherTrack.CopyFromVTrack(selectedMother);
  Double_t d0z0[2],covd0z0[3],d0[2];
  motherTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
  d0[0] = d0z0[0];

  ptMother = selectedMother->Pt();
  momentumMother = selectedMother->P();
  etaMother = selectedMother->Eta();
  phiMother = selectedMother->Phi();

  d0Mother = TMath::Abs(d0[0]);
  d0firstTrack = TMath::Abs(selectedMother->Getd0Prong(0));
  d0secondTrack = TMath::Abs(selectedMother->Getd0Prong(1));
  pointingAngle = selectedMother->CosPointingAngle();
  impactProduct = selectedMother->Prodd0d0();
  impactProductXY = TMath::Abs(selectedMother->ImpParXY());
  invariantMassMother = selectedMother->InvMass(2,prongs);
  dcaMother = selectedMother->GetDCA();
  vertexDistance = vertexMother->DistanceToVertex(primaryVertex);
  angleMotherFirstDaughter = (selectedMother->Px() * firstDaughter->Px() + selectedMother->Py() * firstDaughter->Py() + selectedMother->Pz() * firstDaughter->Pz()) /(selectedMother->P() * firstDaughter->P());
  angleMotherSecondDaughter = (selectedMother->Px() * secondDaughter->Px() + selectedMother->Py() * secondDaughter->Py() + selectedMother->Pz() * secondDaughter->Pz()) /(selectedMother->P() * secondDaughter->P());
  cosThetaStar = selectedMother->CosThetaStar(0,421,211,321);
  angleBetweenBothDaughters  = (firstDaughter->Px() * secondDaughter->Px() + firstDaughter->Py() * secondDaughter->Py() + firstDaughter->Pz() * secondDaughter->Pz()) /(firstDaughter->P() * secondDaughter->P());
  normDecayLength = selectedMother->NormalizedDecayLength();

  Double_t pseudoProperDecayLength = ((vertexMother->GetX() - primaryVertex->GetX()) * selectedMother->Px() / TMath::Abs(selectedMother->Pt())) + ((vertexMother->GetY() - primaryVertex->GetY()) * selectedMother->Py() / TMath::Abs(selectedMother->Pt()));
  Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother/ptMother;
  decayTime = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

  Double_t phi = selectedMother->Phi();
  Double_t theta = selectedMother->Theta();
  Double_t covMatrix[21];
  selectedMother->GetCovarianceXYZPxPyPz(covMatrix);

  Double_t cp = TMath::Cos(phi);
  Double_t sp = TMath::Sin(phi);
  Double_t ct = TMath::Cos(theta);
  Double_t st = TMath::Sin(theta);

  Double_t errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                          +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                          +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                          +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                          +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                          +covMatrix[20]*st*st;  // GetCovPzPz
  Double_t normalizedDecayTime = selectedMother->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

  Double_t eKaon = selectedMother->EProng(1,321);
  Double_t invMassKaon = TMath::Sqrt(eKaon*eKaon-secondDaughter->P()*secondDaughter->P());
  Double_t invMassD0 = selectedMother->InvMassD0();
  invmassDelta = invMassD0 - invMassKaon; 



  ((TH1F*)fMotherHistogramArray[motherType][histType][0])->Fill(ptMother); 
  ((TH1F*)fMotherHistogramArray[motherType][histType][1])->Fill(ptFirstDaughter); 
  ((TH1F*)fMotherHistogramArray[motherType][histType][2])->Fill(ptSecondDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][3])->Fill(etaMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][4])->Fill(phiMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][5])->Fill(d0Mother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][6])->Fill(d0firstTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][7])->Fill(d0secondTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][8])->Fill(pointingAngle);
  ((TH1F*)fMotherHistogramArray[motherType][histType][9])->Fill(impactProduct);
  ((TH1F*)fMotherHistogramArray[motherType][histType][10])->Fill(impactProductXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][11])->Fill(invariantMassMother); 
  ((TH1F*)fMotherHistogramArray[motherType][histType][12])->Fill(invmassDelta); 
  ((TH1F*)fMotherHistogramArray[motherType][histType][13])->Fill(dcaMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][14])->Fill(vertexDistance);
  ((TH1F*)fMotherHistogramArray[motherType][histType][15])->Fill(normDecayLength);
  ((TH1F*)fMotherHistogramArray[motherType][histType][16])->Fill(pseudoProperDecayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][17])->Fill(decayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][18])->Fill(normalizedDecayTime); 
  ((TH1F*)fMotherHistogramArray[motherType][histType][19])->Fill(angleMotherFirstDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][20])->Fill(angleMotherSecondDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][21])->Fill(angleBetweenBothDaughters);
  ((TH1F*)fMotherHistogramArray[motherType][histType][22])->Fill(cosThetaStar);

  ((TH1F*)fMotherHistogramArray[motherType][histType][33])->Fill(dd0pr1);
  ((TH1F*)fMotherHistogramArray[motherType][histType][34])->Fill(dd0pr2);
  ((TH1F*)fMotherHistogramArray[motherType][histType][35])->Fill(dd0max);
  ((TH1F*)fMotherHistogramArray[motherType][histType][36])->Fill(dd0min);

  ((TH2F*)fMotherHistogramArray[motherType][histType][37])->Fill(ptFirstDaughter,ptSecondDaughter);
  ((TH2F*)fMotherHistogramArray[motherType][histType][38])->Fill(ptMother,ptFirstDaughter);
  ((TH2F*)fMotherHistogramArray[motherType][histType][39])->Fill(ptMother,ptSecondDaughter);
  ((TH2F*)fMotherHistogramArray[motherType][histType][40])->Fill(ptMother,d0firstTrack);
  ((TH2F*)fMotherHistogramArray[motherType][histType][41])->Fill(ptMother,d0secondTrack);
  ((TH2F*)fMotherHistogramArray[motherType][histType][42])->Fill(ptMother,pointingAngle);
  ((TH2F*)fMotherHistogramArray[motherType][histType][43])->Fill(ptMother,impactProduct);
  ((TH2F*)fMotherHistogramArray[motherType][histType][44])->Fill(momentumMother,angleMotherFirstDaughter);
  ((TH2F*)fMotherHistogramArray[motherType][histType][45])->Fill(momentumMother,angleMotherSecondDaughter);

  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::FillCascadeMotherHistograms(AliAODRecoCascadeHF * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType){

  if(histType<0) return;

  //In this function we fill the histograms of the reconstructed mothers
  Double_t ptMother = 0.0;
  Double_t momentumMother = 0.0;
  Double_t etaMother = 0.0;
  Double_t phiMother = 0.0;
  Double_t d0Mother = 0.0;
  Double_t d0firstTrack = 0.0;
  Double_t d0secondTrack = 0.0;
  Double_t pointingAngle = 0.0;
  Double_t impactProduct = 0.0;
  Double_t impactProductXY = 0.0;
  Double_t invariantMassMother = 0.0;
  Double_t invmassDelta = 0.0;
  Double_t dcaMother = 0.0;
  AliAODVertex * vertexMother = 0x0;
  AliAODVertex * vertexDaughter = 0x0;
  Double_t vertexDistance = 0.0;
  Double_t normDecayLength = 0.0;
  Double_t decayTime = 0.0;
  Double_t angleMotherFirstDaughter = 0.0;
  Double_t angleMotherSecondDaughter = 0.0;  
  Double_t topomaticFirstDaughter = 0.0;
  Double_t topomaticSecondDaughter = 0.0;  
  Double_t ptFirstDaughter = 0.0;
  Double_t ptSecondDaughter = 0.0;
  UInt_t prongs[2];
  Double_t cosThetaStar = 0;
  Double_t angleBetweenBothDaughters = 0;

  Double_t pdgMassMother = 0;

  if(motherType==1 || motherType==5){ //DStar
    prongs[0] = 211; prongs[1] = 421;
    invmassDelta = DeltaInvMassDStarKpipi(selectedMother);
    AliAODTrack * firstDaughter = (AliAODTrack*)selectedMother->GetDaughter(0);
    AliAODRecoDecayHF2Prong * secondDaughter = (AliAODRecoDecayHF2Prong*)selectedMother->GetDaughter(1);

    ptFirstDaughter = firstDaughter->Pt();
    ptSecondDaughter = secondDaughter->Pt();
    vertexMother = selectedMother->GetSecondaryVtx();
    vertexDaughter = secondDaughter->GetSecondaryVtx();
    angleMotherFirstDaughter = (selectedMother->Px() * firstDaughter->Px() + selectedMother->Py() * firstDaughter->Py() + selectedMother->Pz() * firstDaughter->Pz()) /(selectedMother->P() * firstDaughter->P());
    angleMotherSecondDaughter = (selectedMother->Px() * secondDaughter->Px() + selectedMother->Py() * secondDaughter->Py() + selectedMother->Pz() * secondDaughter->Pz()) /(selectedMother->P() * secondDaughter->P());
    angleBetweenBothDaughters  = (firstDaughter->Px() * secondDaughter->Px() + firstDaughter->Py() * secondDaughter->Py() + firstDaughter->Pz() * secondDaughter->Pz()) /(firstDaughter->P() * secondDaughter->P());
    cosThetaStar = selectedMother->CosThetaStar(0,413,211,421);
    pdgMassMother = TDatabasePDG::Instance()->GetParticle(413)->Mass();
  }
  if(motherType==2){ //B0
    prongs[0] = 211; prongs[1] = 413;
    invmassDelta = DeltaInvMassB0Kpipipi(selectedMother);
    AliAODTrack * firstDaughter = (AliAODTrack*)selectedMother->GetDaughter(0);
    AliAODRecoCascadeHF * secondDaughter = (AliAODRecoCascadeHF*)selectedMother->GetDaughter(1);

    ptFirstDaughter = firstDaughter->Pt();
    ptSecondDaughter = secondDaughter->Pt();
    vertexMother = selectedMother->GetSecondaryVtx();
    vertexDaughter = secondDaughter->GetSecondaryVtx();
    angleMotherFirstDaughter = (selectedMother->Px() * firstDaughter->Px() + selectedMother->Py() * firstDaughter->Py() + selectedMother->Pz() * firstDaughter->Pz()) /(selectedMother->P() * firstDaughter->P());
    angleMotherSecondDaughter = (selectedMother->Px() * secondDaughter->Px() + selectedMother->Py() * secondDaughter->Py() + selectedMother->Pz() * secondDaughter->Pz()) /(selectedMother->P() * secondDaughter->P());
    angleBetweenBothDaughters  = (firstDaughter->Px() * secondDaughter->Px() + firstDaughter->Py() * secondDaughter->Py() + firstDaughter->Pz() * secondDaughter->Pz()) /(firstDaughter->P() * secondDaughter->P());
    cosThetaStar = selectedMother->CosThetaStar(0,511,211,413);
    pdgMassMother = TDatabasePDG::Instance()->GetParticle(511)->Mass();
  }


  //Topomatic
  Double_t dd0pr1=0.;
  Double_t dd0pr2=0.;
  Double_t dd0max=0.;
  Double_t dd0min=0.;
  for(Int_t ipr=0; ipr<2; ipr++) 
  {
    Double_t diffIP, errdiffIP;
    selectedMother->Getd0MeasMinusExpProng(ipr,bz,diffIP,errdiffIP);
    Double_t normdd0=0.;
    if(errdiffIP>0.) normdd0=diffIP/errdiffIP;
    if(ipr==0) dd0pr1=normdd0;
    if(ipr==1) dd0pr2=normdd0;

    // else if(TMath::Abs(normdd0)>TMath::Abs(dd0max)) dd0max=normdd0;
  }
  if(TMath::Abs(dd0pr1)>TMath::Abs(dd0pr2)) {dd0max=dd0pr1; dd0min=dd0pr2;}
  else {dd0max=dd0pr2; dd0min=dd0pr1;}



  AliExternalTrackParam motherTrack;
  motherTrack.CopyFromVTrack(selectedMother);
  Double_t d0z0[2],covd0z0[3],d0[2];
  motherTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0,covd0z0);
  d0[0] = d0z0[0];

  ptMother = selectedMother->Pt();
  momentumMother = selectedMother->P();
  etaMother = selectedMother->Eta();
  phiMother = selectedMother->Phi();
  d0Mother = TMath::Abs(d0[0]);

  pointingAngle = selectedMother->CosPointingAngle();
  impactProduct = selectedMother->Prodd0d0();
  impactProductXY = TMath::Abs(selectedMother->ImpParXY());
  invariantMassMother = selectedMother->InvMass(2,prongs);
  dcaMother = selectedMother->GetDCA();
  vertexDistance = vertexMother->DistanceToVertex(primaryVertex);
  d0firstTrack = TMath::Abs(selectedMother->Getd0Prong(0));
  d0secondTrack = TMath::Abs(selectedMother->Getd0Prong(1));
  normDecayLength = selectedMother->NormalizedDecayLength();

  Double_t pseudoProperDecayLength = ((vertexMother->GetX() - primaryVertex->GetX()) * selectedMother->Px() / TMath::Abs(selectedMother->Pt())) + ((vertexMother->GetY() - primaryVertex->GetY()) * selectedMother->Py() / TMath::Abs(selectedMother->Pt()));
  Double_t pseudoProperDecayTime = pseudoProperDecayLength * pdgMassMother/ptMother;
  decayTime = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

  Double_t phi = selectedMother->Phi();
  Double_t theta = selectedMother->Theta();
  Double_t covMatrix[21];
  selectedMother->GetCovarianceXYZPxPyPz(covMatrix);

  Double_t cp = TMath::Cos(phi);
  Double_t sp = TMath::Sin(phi);
  Double_t ct = TMath::Cos(theta);
  Double_t st = TMath::Sin(theta);

  Double_t errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                          +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                          +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                          +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                          +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                          +covMatrix[20]*st*st;  // GetCovPzPz
  Double_t normalizedDecayTime = selectedMother->NormalizedDecayLength() / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

  ((TH1F*)fMotherHistogramArray[motherType][histType][0])->Fill(ptMother); 
  ((TH1F*)fMotherHistogramArray[motherType][histType][1])->Fill(ptFirstDaughter); 
  ((TH1F*)fMotherHistogramArray[motherType][histType][2])->Fill(ptSecondDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][3])->Fill(etaMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][4])->Fill(phiMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][5])->Fill(d0Mother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][6])->Fill(d0firstTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][7])->Fill(d0secondTrack);
  ((TH1F*)fMotherHistogramArray[motherType][histType][8])->Fill(pointingAngle);
  ((TH1F*)fMotherHistogramArray[motherType][histType][9])->Fill(impactProduct);
  ((TH1F*)fMotherHistogramArray[motherType][histType][10])->Fill(impactProductXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][11])->Fill(invariantMassMother); 
  ((TH1F*)fMotherHistogramArray[motherType][histType][12])->Fill(invmassDelta); 
  ((TH1F*)fMotherHistogramArray[motherType][histType][13])->Fill(dcaMother);
  ((TH1F*)fMotherHistogramArray[motherType][histType][14])->Fill(vertexDistance);
  ((TH1F*)fMotherHistogramArray[motherType][histType][15])->Fill(normDecayLength);
  ((TH1F*)fMotherHistogramArray[motherType][histType][16])->Fill(pseudoProperDecayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][17])->Fill(decayTime);
  ((TH1F*)fMotherHistogramArray[motherType][histType][18])->Fill(normalizedDecayTime); 
  ((TH1F*)fMotherHistogramArray[motherType][histType][19])->Fill(angleMotherFirstDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][20])->Fill(angleMotherSecondDaughter);
  ((TH1F*)fMotherHistogramArray[motherType][histType][21])->Fill(angleBetweenBothDaughters);
  ((TH1F*)fMotherHistogramArray[motherType][histType][22])->Fill(cosThetaStar);

  ((TH1F*)fMotherHistogramArray[motherType][histType][33])->Fill(dd0pr1);
  ((TH1F*)fMotherHistogramArray[motherType][histType][34])->Fill(dd0pr2);
  ((TH1F*)fMotherHistogramArray[motherType][histType][35])->Fill(dd0max);
  ((TH1F*)fMotherHistogramArray[motherType][histType][36])->Fill(dd0min);

  ((TH2F*)fMotherHistogramArray[motherType][histType][37])->Fill(ptFirstDaughter,ptSecondDaughter);
  ((TH2F*)fMotherHistogramArray[motherType][histType][38])->Fill(ptMother,ptFirstDaughter);
  ((TH2F*)fMotherHistogramArray[motherType][histType][39])->Fill(ptMother,ptSecondDaughter);
  ((TH2F*)fMotherHistogramArray[motherType][histType][40])->Fill(ptMother,d0firstTrack);
  ((TH2F*)fMotherHistogramArray[motherType][histType][41])->Fill(ptMother,d0secondTrack);
  ((TH2F*)fMotherHistogramArray[motherType][histType][42])->Fill(ptMother,pointingAngle);
  ((TH2F*)fMotherHistogramArray[motherType][histType][43])->Fill(ptMother,impactProduct);
  ((TH2F*)fMotherHistogramArray[motherType][histType][44])->Fill(momentumMother,angleMotherFirstDaughter);
  ((TH2F*)fMotherHistogramArray[motherType][histType][45])->Fill(momentumMother,angleMotherSecondDaughter);

  if(motherType==1 || motherType==5){
    motherType = motherType -1;

    AliAODRecoDecay* secondDaughter = (AliAODRecoDecay*)selectedMother->GetDaughter(1);
    AliAODTrack * pionD0 = (AliAODTrack*)selectedMother->GetDaughter(0);
    AliAODTrack * kaonD0 = (AliAODTrack*)selectedMother->GetDaughter(1);

    AliAODVertex * vertexDStar = vertexMother;
    AliAODVertex * vertexD0 = secondDaughter->GetSecondaryVtx();
    pdgMassMother = TDatabasePDG::Instance()->GetParticle(421)->Mass();

    AliExternalTrackParam pionD0Track;
    AliExternalTrackParam kaonD0Track;

    Double_t d0z0[2],covd0z0[3],d0[2];

    pionD0Track.CopyFromVTrack(pionD0);
    pionD0Track.PropagateToDCA(vertexDStar,bz,100.,d0z0,covd0z0);
    d0[0] = d0z0[0];

    kaonD0Track.CopyFromVTrack(kaonD0);
    kaonD0Track.PropagateToDCA(vertexDStar,bz,100.,d0z0,covd0z0);
    d0[1] = d0z0[0];

    AliExternalTrackParam D0Track;
    D0Track.CopyFromVTrack(secondDaughter);
    Double_t d0z0D0[2],covd0z0D0[3],d0D0;
    motherTrack.PropagateToDCA(vertexDStar,bz,100.,d0z0D0,covd0z0D0);
    d0D0 = d0z0D0[0];

    Double_t impactProductToDStar = d0[0]*d0[1];
    Double_t impactProductXYToDStar = secondDaughter->ImpParXY(vertexDStar);

    Double_t momentumMother = secondDaughter->P();
    Double_t pointingAngleToDStar = secondDaughter->CosPointingAngle(vertexDStar);
    Double_t d0FirstDaughterToDStar = TMath::Abs(d0[0]);
    Double_t d0SecondDaughterToDStar = TMath::Abs(d0[1]);
    Double_t normDecayLengthToDStar = secondDaughter->NormalizedDecayLength(vertexDStar);

    Double_t pseudoProperDecayLength = ((vertexD0->GetX() - vertexDStar->GetX()) * secondDaughter->Px() / TMath::Abs(secondDaughter->Pt())) + ((vertexD0->GetY() - vertexDStar->GetY()) * secondDaughter->Py() / TMath::Abs(secondDaughter->Pt()));
    Double_t pseudoProperDecayTimeToDStar = pseudoProperDecayLength * pdgMassMother/ptMother;
    Double_t DecayTimeToDStar = vertexDistance / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother/(momentumMother*momentumMother)) + 1)));

    Double_t phi = secondDaughter->Phi();
    Double_t theta = secondDaughter->Theta();
    Double_t covMatrix[21];
    secondDaughter->GetCovarianceXYZPxPyPz(covMatrix);

    Double_t cp = TMath::Cos(phi);
    Double_t sp = TMath::Sin(phi);
    Double_t ct = TMath::Cos(theta);
    Double_t st = TMath::Sin(theta);

    Double_t errorMomentum = covMatrix[9]*cp*cp*ct*ct  // GetCovPxPx
                            +covMatrix[13]*2.*cp*sp*ct*ct  // GetCovPxPy
                            +covMatrix[18]*2.*cp*ct*st  // GetCovPxPz
                            +covMatrix[14]*sp*sp*ct*ct  // GetCovPyPy
                            +covMatrix[19]*2.*sp*ct*st  // GetCovPyPz
                            +covMatrix[20]*st*st;  // GetCovPzPz
    Double_t normDecayTimeToDStar = secondDaughter->NormalizedDecayLength(vertexDStar) / (299792458 * TMath::Sqrt(1/((pdgMassMother*pdgMassMother*errorMomentum*errorMomentum/(momentumMother*momentumMother)) + 1)));

    ((TH1F*)fMotherHistogramArray[motherType][histType][23])->Fill(pointingAngleToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][24])->Fill(d0D0);
    ((TH1F*)fMotherHistogramArray[motherType][histType][25])->Fill(d0FirstDaughterToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][26])->Fill(d0SecondDaughterToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][27])->Fill(impactProductToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][28])->Fill(impactProductXYToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][29])->Fill(normDecayLengthToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][30])->Fill(pseudoProperDecayTimeToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][31])->Fill(DecayTimeToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][32])->Fill(normDecayTimeToDStar);

  }
  return;
}
//-------------------------------------------------------------------------------------
// Bool_t AliAnalysisTaskSEB0toDStarPi::MatchToMonteCarlo(AliAODRecoCascadeHF * mother,TClonesArray * mcTrackArray, finalMotherLabelArray){

//   Int_t i = 0;
//   Int_t numberOfLabels = finalMotherLabelArray->GetEntriesFast();
//   Bool_t finalDaughterCheck = kFALSE;
//   Bool_t decayingDaughterCheck = kTRUE;
//   while(mother->GetDaughter(i)){

//     AliAODTrack* finalDaughter = dynamic_cast<AliAODTrack*>(mother->GetDaughter(i));
//     AliAODRecoCascadeHF* decayingDaughter = dynamic_cast<AliAODRecoCascadeHF*>(mother->GetDaughter(i));

//     std::cout << "finalDaughter: " << finalDaughter << std::endl;
//     std::cout << "decayingDaughter: " << decayingDaughter << std::endl;

//     if(finalDaughter){
//       Int_t mcLabelFinalDaughter = finalDaughter->GetLabel();
//       if(mcLabelFinalDaughter < 0) return kFALSE;
//       AliAODMCParticle *finalParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelFinalDaughter);
//       Int_t labelMother = finalParticle->GetMother();
//       if(labelMother < 0) return kFALSE;
//       while(labelMother >= 0){
//         for (Int_t j = 0; j < numberOfLabels; ++j){
//           if(labelMother == finalMotherLabelArray[j]){finalDaughterCheck = kTRUE; break;}
//         }
//         AliAODMCParticle *mcMother = (AliAODMCParticle*)mcTrackArray->At(labelMother);
//         labelMother = mcMother->GetMother();
//         std::cout << "Mother label: " << labelMother << std::endl;
//       }
//     }

//     if(decayingDaughter) {
//       decayingDaughterCheck = MatchToMonteCarlo(decayingDaughter);
//     }

//     i++;
//   }

//   if(finalDaughterCheck && decayingDaughterCheck) return kTRUE;
//   return kFALSE;
// }
