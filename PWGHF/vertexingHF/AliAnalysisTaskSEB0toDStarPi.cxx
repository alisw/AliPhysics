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
//     Several AliPhysics classes have been used as a basis for this code
//
//-----------------------------------------------------------------------

#include <TSystem.h>
#include <TChain.h>
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
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliVertexingHFUtils.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAnalysisTaskSEB0toDStarPi.h"
#include "AliAODInputHandler.h"
#include <vector>
#include <TMatrix.h>
#include <TVector3.h>
#include <TArrayI.h>
#include <bitset>
#include <TH3F.h>

// #include "TObjectTable.h" 

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskSEB0toDStarPi);
/// \endcond

//__________________________________________________________________________
AliAnalysisTaskSEB0toDStarPi::AliAnalysisTaskSEB0toDStarPi():  
  AliAnalysisTaskSE(),
  fListCuts(0),
  fEvents(0),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputD0FirstDaughter(0),
  fOutputD0SecondDaughter(0),
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
  fDStarPionTracks(0x0),
  fB0PionTracks(0x0),
  fD0Tracks(0x0),
  fShowMask(0),
  fShowRejection(0),
  fnPtBins(0),
  fnPtBinsD0forD0ptbin(0),
  fnPtBinsD0forDStarptbin(0),
  fnPtBinsDStarforDStarptbin(0),
  fnPtBinLimits(0),
  fnPtBinsD0forD0ptbinLimits(0),
  fnPtBinsD0forDStarptbinLimits(0),
  fnPtBinsDStarforDStarptbinLimits(0),
  fPtBinLimits(0x0),
  fPtBinLimitsD0forD0ptbin(0x0),
  fPtBinLimitsD0forDStarptbin(0x0),
  fPtBinLimitsDStarforDStarptbin(0x0),
  fDaughterHistogramArray(),
  fDaughterHistogramArray2D(),
  fDaughterHistogramArrayExtra(),
  fMotherHistogramArray(),
  fMotherHistogramArray2D(),
  fMotherHistogramArrayExtra(), 
  fMotherHistogramArray3D(),
  fUse3DHistograms(0),
  fUpgradeSetting(0),
  fHistMassWindow(0.125),
  fDegreePerRotation(0),
  fNumberOfRotations(0),
  fCheckBackground(0),
  fCheckInjected(1),
  fRemoveInjected(0) 
{
  //
  /// Default ctor
  //

}
//___________________________________________________________________________
AliAnalysisTaskSEB0toDStarPi::AliAnalysisTaskSEB0toDStarPi(const Char_t* name, AliRDHFCutsB0toDStarPi* cuts) :
  AliAnalysisTaskSE(name),
  fListCuts(0),
  fEvents(0),
  fUseMCInfo(kFALSE),
  fOutput(0),
  fOutputD0FirstDaughter(0),
  fOutputD0SecondDaughter(0),
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
  fDStarPionTracks(0x0),
  fB0PionTracks(0x0),
  fD0Tracks(0x0),
  fShowMask(0),
  fShowRejection(0),
  fnPtBins(0),
  fnPtBinsD0forD0ptbin(0),
  fnPtBinsD0forDStarptbin(0),
  fnPtBinsDStarforDStarptbin(0),
  fnPtBinLimits(0),
  fnPtBinsD0forD0ptbinLimits(0),
  fnPtBinsD0forDStarptbinLimits(0),
  fnPtBinsDStarforDStarptbinLimits(0),
  fPtBinLimits(0x0),
  fPtBinLimitsD0forD0ptbin(0x0),
  fPtBinLimitsD0forDStarptbin(0x0),
  fPtBinLimitsDStarforDStarptbin(0x0),
  fDaughterHistogramArray(),
  fDaughterHistogramArray2D(),
  fDaughterHistogramArrayExtra(),
  fMotherHistogramArray(),
  fMotherHistogramArray2D(),
  fMotherHistogramArrayExtra(),
  fMotherHistogramArray3D(),
  fUse3DHistograms(0),
  fUpgradeSetting(0),
  fHistMassWindow(0.125),
  fDegreePerRotation(0),
  fNumberOfRotations(0),
  fCheckBackground(0),
  fCheckInjected(1),
  fRemoveInjected(0)
{
  //
  /// Constructor. Initialization of Inputs and Outputs
  //

  Info("AliAnalysisTaskSEB0toDStarPi","Calling Constructor");

  fCuts = cuts;

  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());   // counters
  DefineOutput(2,TList::Class());   // cut file
  DefineOutput(3,TList::Class());   // D0 pion output
  DefineOutput(4,TList::Class());   // D0 kaon output
  DefineOutput(5,TList::Class());   // DStar pion output
  DefineOutput(6,TList::Class());   // B0 pion output
  DefineOutput(7,TList::Class());   // D0 output
  DefineOutput(8,TList::Class());   // DStar output
  DefineOutput(9,TList::Class());   // B0 output
  DefineOutput(10,TList::Class());  // B0 output
  DefineOutput(11,TList::Class());  // B0 output
  DefineOutput(12,TList::Class());  // B0 output
  DefineOutput(13,TList::Class());  // B0MC output

}

//___________________________________________________________________________
AliAnalysisTaskSEB0toDStarPi::~AliAnalysisTaskSEB0toDStarPi() {
  //
  /// destructor
  //
  Info("~AliAnalysisTaskSEB0toDStarPi","Calling Destructor");

  delete fOutput;
  delete fOutputD0FirstDaughter;
  delete fOutputD0SecondDaughter;
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
  delete fD0Tracks;
  delete fDStarPionTracks;
  delete fB0PionTracks;
  delete fListCuts;
}
//_________________________________________________
void AliAnalysisTaskSEB0toDStarPi::Init(){
  //
  /// Initialization
  //

  if(fDebug > 1) printf("AliAnalysisTaskSEB0toDStarPi::Init() \n");

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

  if(fEvents%50==0){
    std::cout << "\r" << "Analysing event number: " << fEvents << std::endl;
  }

  fEvents++;

  // Show trigger mask
  if(fShowMask)
  {
    std::bitset<32> maskEV(((AliAODInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
    std::cout << "Event mask:   " << maskEV << std::endl;
    std::cout << "Trigger mask: " << std::bitset<32>(fCuts->GetTriggerMask()) << std::endl;
  }

  //==================================================================================
  //  EVENT INITIALIZATION - start
  //==================================================================================


  AliAODEvent * aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TClonesArray * D0TracksFromFriendFile = 0;
  fCEvents->Fill(1);

  if(!aodEvent && AODEvent() && IsStandardAOD()) 
  {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) 
    {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      D0TracksFromFriendFile=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
    }
  } 
  else 
  {
    D0TracksFromFriendFile=(TClonesArray*)aodEvent->GetList()->FindObject("D0toKpi");
  }

  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;
  fCEvents->Fill(2);
  

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD
  TString trigclass=aodEvent->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD")||trigclass.Contains("C0SMH-B-NOPF-ALL")) fCEvents->Fill(5);

  if(!fCuts->IsEventSelected(aodEvent)) 
  {
    if(fShowRejection) std::cout << "Event rejected by code: " << fCuts->GetWhyRejection() << std::endl;
    if(fCuts->GetWhyRejection()==6) // rejected for Z vertex 
    {
      fCEvents->Fill(6);
      return;
    }
  }

  Bool_t isEvSel=fCuts->IsEventSelected(aodEvent);
  if(!isEvSel) return;
  fCEvents->Fill(3);

  //get the magnetic field
  Double_t bz = (Double_t)aodEvent->GetMagneticField(); 

  // AOD primary vertex
  AliAODVertex *primaryVertex = (AliAODVertex*)aodEvent->GetPrimaryVertex();
  if(!primaryVertex) return;
  if(primaryVertex->GetNContributors()<1) return;
  fCEvents->Fill(4);

  if(!D0TracksFromFriendFile)
  {
    AliInfo("Could not find array of HF vertices, skipping the event");
    return;
  }
  else AliDebug(2, Form("Found %d vertices",D0TracksFromFriendFile->GetEntriesFast())); 

  AliAODMCHeader *mcHeader = 0;
  if(fUseMCInfo) 
  { 
    mcHeader = (AliAODMCHeader*)aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf(" MC header branch not found!\n");
      return;
    }
  }


  //==================================================================================
  //  EVENT INITIALIZATION - end
  //==================================================================================
  //  B0 MC SIGNAL IDENTIFICATION - start
  //==================================================================================

  // We create an array that contains all the monte carlo particles in the event
  TClonesArray *mcTrackArray = nullptr; 
  if(fUseMCInfo) mcTrackArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if(fUseMCInfo && !mcTrackArray) return;

  // We create an array to save the MC labels of true signal tracks
  TMatrix * B0toDStarPiLabelMatrix = new TMatrix(0,7);
   
  // We fill the array with all B0->DStarPi tracks
  if(fUseMCInfo) {
    B0toDStarPiSignalTracksInMC(mcTrackArray,aodEvent,B0toDStarPiLabelMatrix,fOutputB0MC);
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

  DStarPionSelection(aodEvent,primaryVertex,bz,mcTrackArray,B0toDStarPiLabelMatrix,mcHeader);
  B0PionSelection(aodEvent,primaryVertex,bz,mcTrackArray,B0toDStarPiLabelMatrix,mcHeader);
  D0Selection(aodEvent,primaryVertex,bz,mcTrackArray,B0toDStarPiLabelMatrix,D0TracksFromFriendFile,mcHeader);

  DStarAndB0Selection(aodEvent,primaryVertex,bz,mcTrackArray,B0toDStarPiLabelMatrix,D0TracksFromFriendFile,mcHeader);

  // Clear arrays and memory management:
  fD0Tracks->erase(fD0Tracks->begin(),fD0Tracks->end());
  fB0PionTracks->erase(fB0PionTracks->begin(),fB0PionTracks->end());
  fDStarPionTracks->erase(fDStarPionTracks->begin(),fDStarPionTracks->end());
  
  delete B0toDStarPiLabelMatrix; B0toDStarPiLabelMatrix = nullptr;

  //==================================================================================
  //  PARTICLE SELECTION LOOP - end
  //==================================================================================

  PostData(1,fOutput);
  PostData(3,fOutputD0FirstDaughter);
  PostData(4,fOutputD0SecondDaughter);
  PostData(5,fOutputDStarPion);
  PostData(6,fOutputB0Pion);
  PostData(7,fOutputD0);
  PostData(8,fOutputDStar);
  PostData(9,fOutputB0);
  PostData(10,fOutputD0_D0Pt);
  PostData(11,fOutputD0_DStarPt);
  PostData(12,fOutputDStar_DStarPt);
  PostData(13,fOutputB0MC);


  //==================================================================================
  //  USER EXECUTION FUNCTION - end
  //==================================================================================
  return;
}
//________________________________________ terminate ___________________________
void AliAnalysisTaskSEB0toDStarPi::Terminate(Option_t*){    
  /// The Terminate() function is the last function to be called during
  /// a query. It always runs on the client, it can be used to present
  /// the results graphically or save the results to file.
  
  // Info("Terminate","");
  AliAnalysisTaskSE::Terminate();
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
  }
  
  fCEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fCEvents"));

  fListCuts = dynamic_cast<TList*> (GetOutputData(2));
  if (!fListCuts) {
    printf("ERROR: fListCuts not available\n");
    return;
  }
  fOutputD0FirstDaughter = dynamic_cast<TList*> (GetOutputData(3));
  if (!fOutputD0FirstDaughter) {
    printf("ERROR: fOutputD0FirstDaughter not available\n");
    return;
  }
  fOutputD0SecondDaughter = dynamic_cast<TList*> (GetOutputData(4));
  if (!fOutputD0SecondDaughter) {
    printf("ERROR: fOutputD0SecondDaughter not available\n");
    return;
  }
  fOutputDStarPion = dynamic_cast<TList*> (GetOutputData(5));
  if (!fOutputDStarPion) {
    printf("ERROR: fOutputDStarPion not available\n");
    return;
  }
  fOutputB0Pion = dynamic_cast<TList*> (GetOutputData(6));
  if (!fOutputB0Pion) {
    printf("ERROR: fOutputB0Pion not available\n");
    return;
  }
  fOutputD0 = dynamic_cast<TList*> (GetOutputData(7));
  if (!fOutputD0) {
    printf("ERROR: fOutputD0 not available\n");
    return;
  }
  fOutputDStar = dynamic_cast<TList*> (GetOutputData(8));
  if (!fOutputDStar) {
    printf("ERROR: fOutputDStar not available\n");
    return;
  }
  fOutputB0 = dynamic_cast<TList*> (GetOutputData(9));
  if (!fOutputB0) {
    printf("ERROR: fOutputB0 not available\n");
    return;
  }
  fOutputD0_D0Pt = dynamic_cast<TList*> (GetOutputData(10));
  if (!fOutputD0_D0Pt) {
    printf("ERROR: fOutputD0_D0Pt not available\n");
    return;
  }
  fOutputD0_DStarPt = dynamic_cast<TList*> (GetOutputData(11));
  if (!fOutputD0_DStarPt) {
    printf("ERROR: fOutputD0_DStarPt not available\n");
    return;
  }
  fOutputDStar_DStarPt = dynamic_cast<TList*> (GetOutputData(12));
  if (!fOutputDStar_DStarPt) {
    printf("ERROR: fOutputDStar_DStarPt not available\n");
    return;
  }  
  fOutputB0MC = dynamic_cast<TList*> (GetOutputData(13));
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

  fOutputD0FirstDaughter = new TList();
  fOutputD0FirstDaughter->SetOwner();
  fOutputD0FirstDaughter->SetName("listD0FirstDaughter");

  fOutputD0SecondDaughter = new TList();
  fOutputD0SecondDaughter->SetOwner();
  fOutputD0SecondDaughter->SetName("listD0SecondDaughter");

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
    
  // we prepare vectors that will save the positions of the daughter tracks in the track list during the reconstruction
  fDStarPionTracks = new std::vector<Int_t>;
  fB0PionTracks = new std::vector<Int_t>;
  fD0Tracks = new std::vector<Int_t>; 

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

  std::cout << "Nr. of B0 meson bins: " <<  fCuts->GetNPtBins() << " limits: " << std::endl;
  for (int i = 0; i < fnPtBinLimits; ++i)
  {
    std::cout << fPtBinLimits[i] << " " << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Nr. of D0 meson bins: " <<  fCuts->GetNPtBinsD0forD0ptbin() << " limits: " << std::endl;
  for (int i = 0; i < fnPtBinsD0forD0ptbinLimits; ++i)
  {
    std::cout << fPtBinLimitsD0forD0ptbin[i] << " " << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Nr. of D0-D* meson bins: " <<  fCuts->GetNPtBinsD0forDStarptbin() << " limits: " << std::endl;
  for (int i = 0; i < fnPtBinsD0forDStarptbinLimits; ++i)
  {
    std::cout << fPtBinLimitsD0forDStarptbin[i] << " " << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Nr. of D* meson bins: " <<  fCuts->GetNPtBinsDStarforDStarptbin() << " limits: " << std::endl;
  for (int i = 0; i < fnPtBinsDStarforDStarptbinLimits; ++i)
  {
    std::cout << fPtBinLimitsDStarforDStarptbin[i] << " " << std::endl;
  }  
  std::cout << std::endl;

  fListCuts=new TList();
  fListCuts->SetOwner();
  fListCuts->SetName("Cuts");
  AliRDHFCutsB0toDStarPi* copyfCuts=new AliRDHFCutsB0toDStarPi(*fCuts);
  // Post the data
  fListCuts->Add(copyfCuts);


  // we create an array of pointers for the histograms. This method is more CPU efficient than looking up each histogram by name.
  // Automatic option is not implemented/complete, the arrays are set manualy in the header file. The array is large enough to accommodate 1 GeV/c pt bins.

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
  // fDaughterHistogramArrayExtra = new Int_t*[numberOfDaughters][numberOfDaughterHistograms2D];
  // fMotherHistogramArray = new Int_t*[numberOfOutputs][numberOfMotherHistogramSets][numberOfMotherHistograms];
  // fMotherHistogramArray2D = new Int_t*[numberOfOutputs][numberOfMotherHistograms2D];

  // define histograms
  DefineHistograms();

  PostData(1,fOutput);
  PostData(2,fListCuts);
  PostData(3,fOutputD0FirstDaughter);
  PostData(4,fOutputD0SecondDaughter);
  PostData(5,fOutputDStarPion);
  PostData(6,fOutputB0Pion);
  PostData(7,fOutputD0);
  PostData(8,fOutputDStar);
  PostData(9,fOutputB0);
  PostData(10,fOutputD0_D0Pt);
  PostData(11,fOutputD0_DStarPt);
  PostData(12,fOutputDStar_DStarPt);
  PostData(13,fOutputB0MC);

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

  TString name_mc_B0_pt_bins ="mc_B0_pt_bins";
  TH1F* hist_mc_B0_pt_bins = new TH1F(name_mc_B0_pt_bins.Data(),"Pt monte carlo B0 in B0->D*#pi; p_{T} [GeV/c]; Entries",fnPtBins,0,fnPtBins);
  hist_mc_B0_pt_bins->Sumw2();
  hist_mc_B0_pt_bins->SetLineColor(6);
  hist_mc_B0_pt_bins->SetMarkerStyle(20);
  hist_mc_B0_pt_bins->SetMarkerSize(0.6);
  hist_mc_B0_pt_bins->SetMarkerColor(6);
  for (Int_t i = 0; i < fnPtBins; ++i)
  {
    TString bin_name = "";
    bin_name += fPtBinLimits[i];
    bin_name += "-";
    bin_name += fPtBinLimits[i+1];
    hist_mc_B0_pt_bins->GetXaxis()->SetBinLabel(i+1,bin_name);
  }
  TH1F* histogram_mc_B0_pt_bins = (TH1F*)hist_mc_B0_pt_bins->Clone();
  fOutputB0MC->Add(histogram_mc_B0_pt_bins);

  TString name_mc_B0_pt_bins_acc ="mc_B0_pt_bins_acc";
  TH1F* hist_mc_B0_pt_bins_acc = new TH1F(name_mc_B0_pt_bins_acc.Data(),"Pt monte carlo B0 in B0->D*#pi; p_{T} [GeV/c]; Entries",fnPtBins,0,fnPtBins);
  hist_mc_B0_pt_bins_acc->Sumw2();
  hist_mc_B0_pt_bins_acc->SetLineColor(6);
  hist_mc_B0_pt_bins_acc->SetMarkerStyle(20);
  hist_mc_B0_pt_bins_acc->SetMarkerSize(0.6);
  hist_mc_B0_pt_bins_acc->SetMarkerColor(6);
  for (Int_t i = 0; i < fnPtBins; ++i)
  {
    TString bin_name = "";
    bin_name += fPtBinLimits[i];
    bin_name += "-";
    bin_name += fPtBinLimits[i+1];
    hist_mc_B0_pt_bins_acc->GetXaxis()->SetBinLabel(i+1,bin_name);
  }
  TH1F* histogram_mc_B0_pt_bins_acc = (TH1F*)hist_mc_B0_pt_bins_acc->Clone();
  fOutputB0MC->Add(histogram_mc_B0_pt_bins_acc);

  TString name_mc_B0_pt_bins_lim_acc ="mc_B0_pt_bins_lim_acc";
  TH1F* hist_mc_B0_pt_bins_lim_acc = new TH1F(name_mc_B0_pt_bins_lim_acc.Data(),"Pt monte carlo B0 in B0->D*#pi; p_{T} [GeV/c]; Entries",fnPtBins,0,fnPtBins);
  hist_mc_B0_pt_bins_lim_acc->Sumw2();
  hist_mc_B0_pt_bins_lim_acc->SetLineColor(6);
  hist_mc_B0_pt_bins_lim_acc->SetMarkerStyle(20);
  hist_mc_B0_pt_bins_lim_acc->SetMarkerSize(0.6);
  hist_mc_B0_pt_bins_lim_acc->SetMarkerColor(6);
  for (Int_t i = 0; i < fnPtBins; ++i)
  {
    TString bin_name = "";
    bin_name += fPtBinLimits[i];
    bin_name += "-";
    bin_name += fPtBinLimits[i+1];
    hist_mc_B0_pt_bins_lim_acc->GetXaxis()->SetBinLabel(i+1,bin_name);
  }
  TH1F* histogram_mc_B0_pt_bins_lim_acc = (TH1F*)hist_mc_B0_pt_bins_lim_acc->Clone();
  fOutputB0MC->Add(histogram_mc_B0_pt_bins_lim_acc);

  TString name_mc_B0_pt_bins_acc_lim_acc ="mc_B0_pt_bins_acc_lim_acc";
  TH1F* hist_mc_B0_pt_bins_acc_lim_acc = new TH1F(name_mc_B0_pt_bins_acc_lim_acc.Data(),"Pt monte carlo B0 in B0->D*#pi; p_{T} [GeV/c]; Entries",fnPtBins,0,fnPtBins);
  hist_mc_B0_pt_bins_acc_lim_acc->Sumw2();
  hist_mc_B0_pt_bins_acc_lim_acc->SetLineColor(6);
  hist_mc_B0_pt_bins_acc_lim_acc->SetMarkerStyle(20);
  hist_mc_B0_pt_bins_acc_lim_acc->SetMarkerSize(0.6);
  hist_mc_B0_pt_bins_acc_lim_acc->SetMarkerColor(6);
  for (Int_t i = 0; i < fnPtBins; ++i)
  {
    TString bin_name = "";
    bin_name += fPtBinLimits[i];
    bin_name += "-";
    bin_name += fPtBinLimits[i+1];
    hist_mc_B0_pt_bins_acc_lim_acc->GetXaxis()->SetBinLabel(i+1,bin_name);
  }
  TH1F* histogram_mc_B0_pt_bins_acc_lim_acc = (TH1F*)hist_mc_B0_pt_bins_acc_lim_acc->Clone();
  fOutputB0MC->Add(histogram_mc_B0_pt_bins_acc_lim_acc);

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

  TString name_mc_B0_rapidity_true ="mc_B0_rapidity_true";
  TH1F* hist_mc_B0_rapidity_true = new TH1F(name_mc_B0_rapidity_true.Data(),"rapidity_true monte carlo B0 in B0->D*#pi; Y; Entries",5000,-20,20);
  hist_mc_B0_rapidity_true->Sumw2();
  hist_mc_B0_rapidity_true->SetLineColor(6);
  hist_mc_B0_rapidity_true->SetMarkerStyle(20);
  hist_mc_B0_rapidity_true->SetMarkerSize(0.6);
  hist_mc_B0_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_B0_rapidity_true = (TH1F*)hist_mc_B0_rapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_B0_rapidity_true);

  TString name_mc_B0_pion_rapidity_true ="mc_B0_pion_rapidity_true";
  TH1F* hist_mc_B0_pion_rapidity_true = new TH1F(name_mc_B0_pion_rapidity_true.Data(),"rapidity_true monte carlo pion of B0 in B0->D*#pi; Y; Entries",5000,-20,20);
  hist_mc_B0_pion_rapidity_true->Sumw2();
  hist_mc_B0_pion_rapidity_true->SetLineColor(6);
  hist_mc_B0_pion_rapidity_true->SetMarkerStyle(20);
  hist_mc_B0_pion_rapidity_true->SetMarkerSize(0.6);
  hist_mc_B0_pion_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_B0_pion_rapidity_true = (TH1F*)hist_mc_B0_pion_rapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_B0_pion_rapidity_true);

  TString name_mc_DStar_rapidity_true ="mc_DStar_rapidity_true";
  TH1F* hist_mc_DStar_rapidity_true = new TH1F(name_mc_DStar_rapidity_true.Data(),"rapidity_true monte carlo DStar in B0->D*#pi; Y; Entries",5000,-20,20);
  hist_mc_DStar_rapidity_true->Sumw2();
  hist_mc_DStar_rapidity_true->SetLineColor(6);
  hist_mc_DStar_rapidity_true->SetMarkerStyle(20);
  hist_mc_DStar_rapidity_true->SetMarkerSize(0.6);
  hist_mc_DStar_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_DStar_rapidity_true = (TH1F*)hist_mc_DStar_rapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_DStar_rapidity_true);

  TString name_mc_DStar_pion_rapidity_true ="mc_DStar_pion_rapidity_true";
  TH1F* hist_mc_DStar_pion_rapidity_true = new TH1F(name_mc_DStar_pion_rapidity_true.Data(),"rapidity_true monte carlo pion of DStar in B0->D*#pi; Y; Entries",5000,-20,20);
  hist_mc_DStar_pion_rapidity_true->Sumw2();
  hist_mc_DStar_pion_rapidity_true->SetLineColor(6);
  hist_mc_DStar_pion_rapidity_true->SetMarkerStyle(20);
  hist_mc_DStar_pion_rapidity_true->SetMarkerSize(0.6);
  hist_mc_DStar_pion_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_DStar_pion_rapidity_true = (TH1F*)hist_mc_DStar_pion_rapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_DStar_pion_rapidity_true);

  TString name_mc_D0_rapidity_true ="mc_D0_rapidity_true";
  TH1F* hist_mc_D0_rapidity_true = new TH1F(name_mc_D0_rapidity_true.Data(),"rapidity_true monte carlo D0 in B0->D*#pi; Y; Entries",5000,-20,20);
  hist_mc_D0_rapidity_true->Sumw2();
  hist_mc_D0_rapidity_true->SetLineColor(6);
  hist_mc_D0_rapidity_true->SetMarkerStyle(20);
  hist_mc_D0_rapidity_true->SetMarkerSize(0.6);
  hist_mc_D0_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_D0_rapidity_true = (TH1F*)hist_mc_D0_rapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_D0_rapidity_true);

  TString name_mc_D0_pion_rapidity_true ="mc_D0_pion_rapidity_true";
  TH1F* hist_mc_D0_pion_rapidity_true = new TH1F(name_mc_D0_pion_rapidity_true.Data(),"rapidity_true monte carlo pion of D0 in B0->D*#pi; Y; Entries",5000,-20,20);
  hist_mc_D0_pion_rapidity_true->Sumw2();
  hist_mc_D0_pion_rapidity_true->SetLineColor(6);
  hist_mc_D0_pion_rapidity_true->SetMarkerStyle(20);
  hist_mc_D0_pion_rapidity_true->SetMarkerSize(0.6);
  hist_mc_D0_pion_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_D0_pion_rapidity_true = (TH1F*)hist_mc_D0_pion_rapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_D0_pion_rapidity_true);

  TString name_mc_D0_kaon_rapidity_true ="mc_D0_kaon_rapidity_true";
  TH1F* hist_mc_D0_kaon_rapidity_true = new TH1F(name_mc_D0_kaon_rapidity_true.Data(),"rapidity_true monte carlo kaon of D0 in B0->D*#pi; Y; Entries",5000,-20,20);
  hist_mc_D0_kaon_rapidity_true->Sumw2();
  hist_mc_D0_kaon_rapidity_true->SetLineColor(6);
  hist_mc_D0_kaon_rapidity_true->SetMarkerStyle(20);
  hist_mc_D0_kaon_rapidity_true->SetMarkerSize(0.6);
  hist_mc_D0_kaon_rapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_D0_kaon_rapidity_true = (TH1F*)hist_mc_D0_kaon_rapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_D0_kaon_rapidity_true);

  TString name_mc_B0_pseudorapidity_true ="mc_B0_pseudorapidity_true";
  TH1F* hist_mc_B0_pseudorapidity_true = new TH1F(name_mc_B0_pseudorapidity_true.Data(),"pseudorapidity_true monte carlo B0 in B0->D*#pi; #eta; Entries",5000,-20,20);
  hist_mc_B0_pseudorapidity_true->Sumw2();
  hist_mc_B0_pseudorapidity_true->SetLineColor(6);
  hist_mc_B0_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_B0_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_B0_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_B0_pseudorapidity_true = (TH1F*)hist_mc_B0_pseudorapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_B0_pseudorapidity_true);

  TString name_mc_B0_pion_pseudorapidity_true ="mc_B0_pion_pseudorapidity_true";
  TH1F* hist_mc_B0_pion_pseudorapidity_true = new TH1F(name_mc_B0_pion_pseudorapidity_true.Data(),"pseudorapidity_true monte carlo pion of B0 in B0->D*#pi; #eta; Entries",5000,-20,20);
  hist_mc_B0_pion_pseudorapidity_true->Sumw2();
  hist_mc_B0_pion_pseudorapidity_true->SetLineColor(6);
  hist_mc_B0_pion_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_B0_pion_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_B0_pion_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_B0_pion_pseudorapidity_true = (TH1F*)hist_mc_B0_pion_pseudorapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_B0_pion_pseudorapidity_true);

  TString name_mc_DStar_pseudorapidity_true ="mc_DStar_pseudorapidity_true";
  TH1F* hist_mc_DStar_pseudorapidity_true = new TH1F(name_mc_DStar_pseudorapidity_true.Data(),"pseudorapidity_true monte carlo DStar in B0->D*#pi; #eta; Entries",5000,-20,20);
  hist_mc_DStar_pseudorapidity_true->Sumw2();
  hist_mc_DStar_pseudorapidity_true->SetLineColor(6);
  hist_mc_DStar_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_DStar_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_DStar_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_DStar_pseudorapidity_true = (TH1F*)hist_mc_DStar_pseudorapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_DStar_pseudorapidity_true);

  TString name_mc_DStar_pion_pseudorapidity_true ="mc_DStar_pion_pseudorapidity_true";
  TH1F* hist_mc_DStar_pion_pseudorapidity_true = new TH1F(name_mc_DStar_pion_pseudorapidity_true.Data(),"pseudorapidity_true monte carlo pion of DStar in B0->D*#pi; #eta; Entries",5000,-20,20);
  hist_mc_DStar_pion_pseudorapidity_true->Sumw2();
  hist_mc_DStar_pion_pseudorapidity_true->SetLineColor(6);
  hist_mc_DStar_pion_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_DStar_pion_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_DStar_pion_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_DStar_pion_pseudorapidity_true = (TH1F*)hist_mc_DStar_pion_pseudorapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_DStar_pion_pseudorapidity_true);

  TString name_mc_D0_pseudorapidity_true ="mc_D0_pseudorapidity_true";
  TH1F* hist_mc_D0_pseudorapidity_true = new TH1F(name_mc_D0_pseudorapidity_true.Data(),"pseudorapidity_true monte carlo D0 in B0->D*#pi; #eta; Entries",5000,-20,20);
  hist_mc_D0_pseudorapidity_true->Sumw2();
  hist_mc_D0_pseudorapidity_true->SetLineColor(6);
  hist_mc_D0_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_D0_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_D0_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_D0_pseudorapidity_true = (TH1F*)hist_mc_D0_pseudorapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_D0_pseudorapidity_true);

  TString name_mc_D0_pion_pseudorapidity_true ="mc_D0_pion_pseudorapidity_true";
  TH1F* hist_mc_D0_pion_pseudorapidity_true = new TH1F(name_mc_D0_pion_pseudorapidity_true.Data(),"pseudorapidity_true monte carlo pion of D0 in B0->D*#pi; #eta; Entries",5000,-20,20);
  hist_mc_D0_pion_pseudorapidity_true->Sumw2();
  hist_mc_D0_pion_pseudorapidity_true->SetLineColor(6);
  hist_mc_D0_pion_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_D0_pion_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_D0_pion_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_D0_pion_pseudorapidity_true = (TH1F*)hist_mc_D0_pion_pseudorapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_D0_pion_pseudorapidity_true);

  TString name_mc_D0_kaon_pseudorapidity_true ="mc_D0_kaon_pseudorapidity_true";
  TH1F* hist_mc_D0_kaon_pseudorapidity_true = new TH1F(name_mc_D0_kaon_pseudorapidity_true.Data(),"pseudorapidity_true monte carlo kaon of D0 in B0->D*#pi; #eta; Entries",5000,-20,20);
  hist_mc_D0_kaon_pseudorapidity_true->Sumw2();
  hist_mc_D0_kaon_pseudorapidity_true->SetLineColor(6);
  hist_mc_D0_kaon_pseudorapidity_true->SetMarkerStyle(20);
  hist_mc_D0_kaon_pseudorapidity_true->SetMarkerSize(0.6);
  hist_mc_D0_kaon_pseudorapidity_true->SetMarkerColor(6);
  TH1F* histogram_mc_D0_kaon_pseudorapidity_true = (TH1F*)hist_mc_D0_kaon_pseudorapidity_true->Clone();
  fOutputB0MC->Add(histogram_mc_D0_kaon_pseudorapidity_true);

 //==================================================

  TString name_dca_D0_DStarPion ="dca_D0_DStarPion";
  TH1F* hist_dca_D0_DStarPion = new TH1F(name_dca_D0_DStarPion.Data(),"dca_D0_DStarPion; DCA [cm]; Entries",1000,0,0.2);
  hist_dca_D0_DStarPion->Sumw2();
  hist_dca_D0_DStarPion->SetLineColor(6);
  hist_dca_D0_DStarPion->SetMarkerStyle(20);
  hist_dca_D0_DStarPion->SetMarkerSize(0.6);
  hist_dca_D0_DStarPion->SetMarkerColor(6);
  TH1F* histogram_dca_D0_DStarPion = (TH1F*)hist_dca_D0_DStarPion->Clone();
  fOutputB0MC->Add(histogram_dca_D0_DStarPion);

  TString name_dca_D0_B0Pion ="dca_D0_B0Pion";
  TH1F* hist_dca_D0_B0Pion = new TH1F(name_dca_D0_B0Pion.Data(),"dca_D0_B0Pion; DCA [cm]; Entries",1000,0,0.2);
  hist_dca_D0_B0Pion->Sumw2();
  hist_dca_D0_B0Pion->SetLineColor(6);
  hist_dca_D0_B0Pion->SetMarkerStyle(20);
  hist_dca_D0_B0Pion->SetMarkerSize(0.6);
  hist_dca_D0_B0Pion->SetMarkerColor(6);
  TH1F* histogram_dca_D0_B0Pion = (TH1F*)hist_dca_D0_B0Pion->Clone();
  fOutputB0MC->Add(histogram_dca_D0_B0Pion);

  TString name_dca_DStarPion_B0Pion ="dca_DStarPion_B0Pion";
  TH1F* hist_dca_DStarPion_B0Pion = new TH1F(name_dca_DStarPion_B0Pion.Data(),"dca_DStarPion_B0Pion; DCA [cm]; Entries",1000,0,0.2);
  hist_dca_DStarPion_B0Pion->Sumw2();
  hist_dca_DStarPion_B0Pion->SetLineColor(6);
  hist_dca_DStarPion_B0Pion->SetMarkerStyle(20);
  hist_dca_DStarPion_B0Pion->SetMarkerSize(0.6);
  hist_dca_DStarPion_B0Pion->SetMarkerColor(6);
  TH1F* histogram_dca_DStarPion_B0Pion = (TH1F*)hist_dca_DStarPion_B0Pion->Clone();
  fOutputB0MC->Add(histogram_dca_DStarPion_B0Pion);

  TString name_dca_Combined ="dca_Combined";
  TH1F* hist_dca_Combined = new TH1F(name_dca_Combined.Data(),"dca_Combined; DCA [cm]; Entries",1000,0,1.0);
  hist_dca_Combined->Sumw2();
  hist_dca_Combined->SetLineColor(6);
  hist_dca_Combined->SetMarkerStyle(20);
  hist_dca_Combined->SetMarkerSize(0.6);
  hist_dca_Combined->SetMarkerColor(6);
  TH1F* histogram_dca_Combined = (TH1F*)hist_dca_Combined->Clone();
  fOutputB0MC->Add(histogram_dca_Combined);

  TString name_dca_Signal_D0_DStarPion ="dca_Signal_D0_DStarPion";
  TH1F* hist_dca_Signal_D0_DStarPion = new TH1F(name_dca_Signal_D0_DStarPion.Data(),"dca_Signal_D0_DStarPion; DCA [cm]; Entries",1000,0,0.2);
  hist_dca_Signal_D0_DStarPion->Sumw2();
  hist_dca_Signal_D0_DStarPion->SetLineColor(4);
  hist_dca_Signal_D0_DStarPion->SetMarkerStyle(20);
  hist_dca_Signal_D0_DStarPion->SetMarkerSize(0.6);
  hist_dca_Signal_D0_DStarPion->SetMarkerColor(4);
  TH1F* histogram_dca_Signal_D0_DStarPion = (TH1F*)hist_dca_Signal_D0_DStarPion->Clone();
  fOutputB0MC->Add(histogram_dca_Signal_D0_DStarPion);

  TString name_dca_Signal_D0_B0Pion ="dca_Signal_D0_B0Pion";
  TH1F* hist_dca_Signal_D0_B0Pion = new TH1F(name_dca_Signal_D0_B0Pion.Data(),"dca_Signal_D0_B0Pion; DCA [cm]; Entries",1000,0,0.2);
  hist_dca_Signal_D0_B0Pion->Sumw2();
  hist_dca_Signal_D0_B0Pion->SetLineColor(4);
  hist_dca_Signal_D0_B0Pion->SetMarkerStyle(20);
  hist_dca_Signal_D0_B0Pion->SetMarkerSize(0.6);
  hist_dca_Signal_D0_B0Pion->SetMarkerColor(4);
  TH1F* histogram_dca_Signal_D0_B0Pion = (TH1F*)hist_dca_Signal_D0_B0Pion->Clone();
  fOutputB0MC->Add(histogram_dca_Signal_D0_B0Pion);

  TString name_dca_Signal_DStarPion_B0Pion ="dca_Signal_DStarPion_B0Pion";
  TH1F* hist_dca_Signal_DStarPion_B0Pion = new TH1F(name_dca_Signal_DStarPion_B0Pion.Data(),"dca_Signal_DStarPion_B0Pion; DCA [cm]; Entries",1000,0,0.2);
  hist_dca_Signal_DStarPion_B0Pion->Sumw2();
  hist_dca_Signal_DStarPion_B0Pion->SetLineColor(4);
  hist_dca_Signal_DStarPion_B0Pion->SetMarkerStyle(20);
  hist_dca_Signal_DStarPion_B0Pion->SetMarkerSize(0.6);
  hist_dca_Signal_DStarPion_B0Pion->SetMarkerColor(4);
  TH1F* histogram_dca_Signal_DStarPion_B0Pion = (TH1F*)hist_dca_Signal_DStarPion_B0Pion->Clone();
  fOutputB0MC->Add(histogram_dca_Signal_DStarPion_B0Pion);

  TString name_dca_Signal_Combined ="dca_Signal_Combined";
  TH1F* hist_dca_Signal_Combined = new TH1F(name_dca_Signal_Combined.Data(),"dca_Signal_Combined; DCA [cm]; Entries",1000,0,1.0);
  hist_dca_Signal_Combined->Sumw2();
  hist_dca_Signal_Combined->SetLineColor(4);
  hist_dca_Signal_Combined->SetMarkerStyle(20);
  hist_dca_Signal_Combined->SetMarkerSize(0.6);
  hist_dca_Signal_Combined->SetMarkerColor(4);
  TH1F* histogram_dca_Signal_Combined = (TH1F*)hist_dca_Signal_Combined->Clone();
  fOutputB0MC->Add(histogram_dca_Signal_Combined);

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
  TH1F* hist_B0s_per_bin = new TH1F(name_B0s_per_bin.Data(),"Number of B0 to kpipipi in the Analysis per bin; Entries",fnPtBins,0,fnPtBins); 
  for (Int_t i = 0; i < fnPtBins; ++i)
  {
    TString bin_name = "";
    bin_name += fPtBinLimits[i];
    bin_name += "-";
    bin_name += fPtBinLimits[i+1];
    hist_B0s_per_bin->GetXaxis()->SetBinLabel(i+1,bin_name);
  }
  TH1F* hist_B0s_per_bin_mc = (TH1F*)hist_B0s_per_bin->Clone();
  fOutputB0MC->Add(hist_B0s_per_bin_mc);

  TString name_B0s_per_bin_in_Acc ="B0s_per_bin_in_Acc";
  TH1F* hist_B0s_per_bin_in_Acc = new TH1F(name_B0s_per_bin_in_Acc.Data(),"Number of B0 to kpipipi in the Analysis per bin with all daughters in acceptance; Entries",fnPtBins,0,fnPtBins); 
  for (Int_t i = 0; i < fnPtBins; ++i)
  {
    TString bin_name = "";
    bin_name += fPtBinLimits[i];
    bin_name += "-";
    bin_name += fPtBinLimits[i+1];
    hist_B0s_per_bin_in_Acc->GetXaxis()->SetBinLabel(i+1,bin_name);
  }
  TH1F* hist_B0s_per_bin_in_Acc_mc = (TH1F*)hist_B0s_per_bin_in_Acc->Clone();
  fOutputB0MC->Add(hist_B0s_per_bin_in_Acc_mc);

  TString name_B0s_per_bin_in_Lim_Acc ="B0s_per_bin_in_Lim_Acc";
  TH1F* hist_B0s_per_bin_in_Lim_Acc = new TH1F(name_B0s_per_bin_in_Lim_Acc.Data(),"Number of B0 to kpipipi in the Analysis per bin with all daughters in acceptance; Entries",fnPtBins,0,fnPtBins); 
  for (Int_t i = 0; i < fnPtBins; ++i)
  {
    TString bin_name = "";
    bin_name += fPtBinLimits[i];
    bin_name += "-";
    bin_name += fPtBinLimits[i+1];
    hist_B0s_per_bin_in_Lim_Acc->GetXaxis()->SetBinLabel(i+1,bin_name);
  }
  TH1F* hist_B0s_per_bin_in_Lim_Acc_mc = (TH1F*)hist_B0s_per_bin_in_Lim_Acc->Clone();
  fOutputB0MC->Add(hist_B0s_per_bin_in_Lim_Acc_mc);

 //======================================================================================================================================================

  //we make the histograms for the Pions and Kaon
  for (Int_t i = 0; i < 4; i++){
    
    TString add_name = "";
    TList * listout;
    if(i==0) listout = fOutputD0FirstDaughter;
    if(i==1) listout = fOutputD0SecondDaughter;
    if(i==2) listout = fOutputDStarPion;
    if(i==3) listout = fOutputB0Pion;

    for (Int_t j = 0; j < 6; j++){
      if(j==0) add_name = "";
      if(j==1) add_name = "Signal";
      if(j==2) add_name = "Cut";
      if(j==3) add_name = "SignalCut";
      if(j==4) add_name = "Result";
      if(j==5) add_name = "SignalResult";

      TString name_Histogram = "";
      TString discription_Histogram = "";
      Int_t numberOfBins = 0;
      Double_t lowerBound = 0.0;
      Double_t upperBound = 0.0;

      for (Int_t k = 0; k < 9; ++k)
      {
        if(k==0){name_Histogram = "ptTrack"; discription_Histogram = "pt track; p_{T} [GeV/c]; Entries"; numberOfBins = 600; lowerBound = 0; upperBound = 30;}
        if(k==1){name_Histogram = "momentumTrack"; discription_Histogram = "momentum track; p [GeV/c]; Entries"; numberOfBins = 600; lowerBound = 0; upperBound = 30;}
        if(k==2){name_Histogram = "numberOfITS"; discription_Histogram = "Number of ITS clusters track; [#]; Entries"; numberOfBins = 10; lowerBound = -0.5; upperBound = 9.5;}
        if(k==3){name_Histogram = "numberOfTPC"; discription_Histogram = "Number of TPC clusters track; [#]; Entries"; numberOfBins = 601; lowerBound = -0.5; upperBound = 600.5;}
        if(k==4){name_Histogram = "pointsOnITS"; discription_Histogram = "Number of ITS clusters track per layer; [#]; Entries"; numberOfBins = 10; lowerBound = -0.5; upperBound = 9.5;}
        if(k==5){name_Histogram = "nSigmaTPC"; discription_Histogram = "n sigma TPC for track PID; sigma; Entries"; numberOfBins = 500; lowerBound = -5; upperBound = 5;}
        if(k==6){name_Histogram = "nSigmaTOF"; discription_Histogram = "n sigma TOF for track PID; sigma; Entries"; numberOfBins = 500; lowerBound = -5; upperBound = 5;}
        if(k==7){name_Histogram = "nSigmaTPCandTOF"; discription_Histogram = "n sigma TPC and TOF for track PID; a.u.; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 10;}
        if(k==8){name_Histogram = "impactParameter"; discription_Histogram = "Impact Parameter track;  [cm]; Entries"; numberOfBins = 2000; lowerBound = 0; upperBound = 0.5;}

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
        fDaughterHistogramArray[i][j][k] = histogram_Clone;
      }

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
      fDaughterHistogramArray[i][j][12] = histogram_numberofparticlesperevent;
    }

    TH1F * effectOfCuts = new TH1F("effectOfCutsOnBackground","Removal counter",18,0,18);
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
    fDaughterHistogramArrayExtra[i][0] = effectOfCuts;

    TH1F * effectOfCutsMC = new TH1F("effectOfCutsOnSignal","Removal counter",18,0,18);
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
    fDaughterHistogramArrayExtra[i][1] = effectOfCutsMC;

    TString name_particle_pdg ="particle_pdg";
    TH1F* hist_particle_pdg = new TH1F(name_particle_pdg.Data(),"Pdg code particle; pdg code; Entries",2000,-0.5,1999.5);
    hist_particle_pdg->Sumw2();
    hist_particle_pdg->SetLineColor(6);
    hist_particle_pdg->SetMarkerStyle(20);
    hist_particle_pdg->SetMarkerSize(0.6);
    hist_particle_pdg->SetMarkerColor(6);
    TH1F* histogram_particle_pdg = (TH1F*)hist_particle_pdg->Clone();
    listout->Add(histogram_particle_pdg);
    fDaughterHistogramArrayExtra[i][2] = histogram_particle_pdg;

    TString name_particle_mother_pdg ="particle_mother_pdg";
    TH1F* hist_particle_mother_pdg = new TH1F(name_particle_mother_pdg.Data(),"Pdg code particle mother; pdg code; Entries",2000,-0.5,1999.5);
    hist_particle_mother_pdg->Sumw2();
    hist_particle_mother_pdg->SetLineColor(6);
    hist_particle_mother_pdg->SetMarkerStyle(20);
    hist_particle_mother_pdg->SetMarkerSize(0.6);
    hist_particle_mother_pdg->SetMarkerColor(6);
    TH1F* histogram_particle_mother_pdg = (TH1F*)hist_particle_mother_pdg->Clone();
    listout->Add(histogram_particle_mother_pdg);
    fDaughterHistogramArrayExtra[i][3] = histogram_particle_mother_pdg;

    TString name_ptB0_vs_ptTrack ="ptB0_vs_ptTrackBackground";
    TH2F* hist_ptB0_vs_ptTrack = new TH2F(name_ptB0_vs_ptTrack.Data(),"Pt B0 vs Pt ; p_{T} B0 [GeV/c]; p_{T} track [GeV/c]",100,0,30,100,0,30);
    hist_ptB0_vs_ptTrack->Sumw2();
    hist_ptB0_vs_ptTrack->SetLineColor(6);
    hist_ptB0_vs_ptTrack->SetMarkerStyle(20);
    hist_ptB0_vs_ptTrack->SetMarkerSize(0.6);
    hist_ptB0_vs_ptTrack->SetMarkerColor(6);
    TH2F* histogram_ptB0_vs_ptTrack = (TH2F*)hist_ptB0_vs_ptTrack->Clone();
    listout->Add(histogram_ptB0_vs_ptTrack);
    fDaughterHistogramArray2D[i][4] = histogram_ptB0_vs_ptTrack;

    TString name_ptB0_vs_ptTrackMC ="ptB0_vs_ptTrackSignal";
    TH2F* hist_ptB0_vs_ptTrackMC = new TH2F(name_ptB0_vs_ptTrackMC.Data(),"Pt B0 vs Pt ; p_{T} B0 [GeV/c]; p_{T} track [GeV/c]",100,0,30,100,0,30);
    hist_ptB0_vs_ptTrackMC->Sumw2();
    hist_ptB0_vs_ptTrackMC->SetLineColor(4);
    hist_ptB0_vs_ptTrackMC->SetMarkerStyle(20);
    hist_ptB0_vs_ptTrackMC->SetMarkerSize(0.6);
    hist_ptB0_vs_ptTrackMC->SetMarkerColor(6);
    TH2F* histogram_ptB0_vs_ptTrackMC = (TH2F*)hist_ptB0_vs_ptTrackMC->Clone();
    listout->Add(histogram_ptB0_vs_ptTrackMC);  
    fDaughterHistogramArray2D[i][5] = histogram_ptB0_vs_ptTrackMC;  
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
        if(j==1) add_name = "Signal";
        if(j==2) add_name = "Cut";
        if(j==3) add_name = "SignalCut";
        if(j==4) add_name = "Result";
        if(j==5) add_name = "SignalResult";
        if(j%2==0 && j>5) {add_name = "_ptbin_"; add_name += fPtBinLimits[(j-6)/2]; add_name += "_to_"; add_name += fPtBinLimits[(j-6)/2 + 1];}
        if(j%2==1 && j>5) {add_name = "Signal_ptbin_"; add_name += fPtBinLimits[(j-7)/2]; add_name += "_to_"; add_name += fPtBinLimits[(j-7)/2 + 1];}
      }
      if(i==3)
      { 
        if(j%2==0) {add_name = "_ptbin_"; add_name += fPtBinLimitsD0forD0ptbin[j/2]; add_name += "_to_"; add_name += fPtBinLimitsD0forD0ptbin[1+j/2];}
        if(j%2==1) {add_name = "Signal_ptbin_"; add_name += fPtBinLimitsD0forD0ptbin[(j-1)/2]; add_name += "_to_"; add_name += fPtBinLimitsD0forD0ptbin[1+(j-1)/2];}
      }
      if(i==4)
      {
        if(j%2==0) {add_name = "_ptbin_"; add_name += fPtBinLimitsD0forDStarptbin[j/2]; add_name += "_to_"; add_name += fPtBinLimitsD0forDStarptbin[1+j/2];}
        if(j%2==1) {add_name = "Signal_ptbin_"; add_name += fPtBinLimitsD0forDStarptbin[(j-1)/2]; add_name += "_to_"; add_name += fPtBinLimitsD0forDStarptbin[1+(j-1)/2];}
      }
      if(i==5)
      {
        if(j%2==0) {add_name = "_ptbin_"; add_name += fPtBinLimitsDStarforDStarptbin[j/2]; add_name += "_to_"; add_name += fPtBinLimitsDStarforDStarptbin[1+j/2];}
        if(j%2==1) {add_name = "Signal_ptbin_"; add_name += fPtBinLimitsDStarforDStarptbin[(j-1)/2]; add_name += "_to_"; add_name += fPtBinLimitsDStarforDStarptbin[1+(j-1)/2];}
      }


      TString name_Histogram = "";
      TString discription_Histogram  = "";
      Int_t numberOfBins = 0;
      Double_t lowerBound = 0.0;
      Double_t upperBound = 0.0;
      Int_t numberOfBinsTwo = 0;
      Double_t lowerBoundTwo = 0.0;
      Double_t upperBoundTwo = 0.0;

      for (Int_t k = 0; k < 43; ++k)
      {
        if(k==0){name_Histogram = "ptMother"; discription_Histogram = "pt mother; p_{T} [GeV/c]; Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
        if(k==1){name_Histogram = "ptFirstDaughter"; discription_Histogram = "pt first daughter; p_{T} [GeV/c]; Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
        if(k==2){name_Histogram = "ptSecondDaughter"; discription_Histogram = "pt second daughter; p_{T} [GeV/c]; Entries"; numberOfBins = 300; lowerBound = 0; upperBound = 30;}
        if(k==3){name_Histogram = "etaMother"; discription_Histogram = "eta mother; #eta; Entries"; numberOfBins = 100; lowerBound = -2; upperBound = 2;}
        if(k==4){name_Histogram = "phiMother"; discription_Histogram = "phi mother; #phi; Entries"; numberOfBins = 25; lowerBound = 0; upperBound = 2*TMath::Pi();}
        if(k==5){name_Histogram = "d0Mother"; discription_Histogram = "d0 mother;  [cm]; Entries"; numberOfBins = 2000; lowerBound = 0; upperBound = 0.5;}
        if(k==6){name_Histogram = "d0FirstDaughter"; discription_Histogram = "d0 first daughter;  [cm]; Entries"; numberOfBins = 2000; lowerBound = 0; upperBound = 0.5;}

        if(k==7){name_Histogram = "d0SecondDaughter"; discription_Histogram = "d0 second daughter;  [cm]; Entries"; numberOfBins = 2000; lowerBound = 0; upperBound = 0.5;}

        if(k==8){name_Histogram = "pointingAngleMother"; discription_Histogram = "pointing angle;  [Cos(#theta)]; Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
        if(k==9){name_Histogram = "impactProduct"; discription_Histogram = "impact product; [cm^{2}]; Entries"; numberOfBins = 500; lowerBound = -0.01; upperBound = 0.01;}
        if(k==10){name_Histogram = "impactProductXY"; discription_Histogram = "impact product XY; [cm^{2}]; Entries"; numberOfBins = 400; lowerBound = 0; upperBound = 0.1;}
        if(k==11){name_Histogram = "invariantMassMother"; discription_Histogram = "mass mother candidate; m [GeV/c^{2}]; Entries"; numberOfBins = 20000; lowerBound = 0; upperBound = 10;}
        if(k==12){name_Histogram = "deltaMassMother"; discription_Histogram = "mass mother candidate; m [GeV/c^{2}]; Entries"; numberOfBins = 10000; lowerBound = 0; upperBound = 10;}
        if(k==13){name_Histogram = "dcaMother"; discription_Histogram = "dca mother; distance [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 0.25;}
        if(k==14){name_Histogram = "vertexDistance"; discription_Histogram = "vertex distance between mother and primary vertex; distance [cm]; Entries"; numberOfBins = 500; lowerBound = 0; upperBound = 1;}
        if(k==15){name_Histogram = "normDecayLength"; discription_Histogram = "Normalized decay length w.r.t primary vertex; [cm]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 50;}
        if(k==16){name_Histogram = "pseudoProperDecayTime"; discription_Histogram = "Pseudo Proper Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 1000; lowerBound = -10; upperBound = 10;}
        if(k==17){name_Histogram = "DecayTime"; discription_Histogram = "Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.00000001;}
        if(k==18){name_Histogram = "normDecayTime"; discription_Histogram = "Normalized Decay Time w.r.t primary vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.0000001;}
        if(k==19){name_Histogram = "angleMotherFirstDaughter"; discription_Histogram = "flight angle mother and first daughter; [Cos(#phi)]; Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
        if(k==20){name_Histogram = "angleMotherSecondDaughter"; discription_Histogram = "flight angle mother and second daughter; [Cos(#phi)]; Entries"; numberOfBins = 100; lowerBound = 0.5; upperBound = 1;}
        if(k==21){name_Histogram = "angleBetweenBothDaughters"; discription_Histogram = "angle between both daughters; [Cos(#phi)]; Entries"; numberOfBins = 100; lowerBound = -1; upperBound = 1;}
        if(k==22){name_Histogram = "cosThetaStar"; discription_Histogram = "cosThetaStar; [Cos(#theta*)]; Entries"; numberOfBins = 200; lowerBound = -2; upperBound = 2;}
        if(k==23){name_Histogram = "vertexX"; discription_Histogram = "Vertex position; [cm]; Entries"; numberOfBins = 200; lowerBound = -5; upperBound = 5;}
        if(k==24){name_Histogram = "vertexY"; discription_Histogram = "Vertex position; [cm]; Entries"; numberOfBins = 200; lowerBound = -5; upperBound = 5;}
        if(k==25){name_Histogram = "vertexZ"; discription_Histogram = "Vertex position; [cm]; Entries"; numberOfBins = 200; lowerBound = -20; upperBound = 20;}


        if(k==26){if(i==0 || i==3 || i==4){name_Histogram = "pointingAngleToDStar"; discription_Histogram = "Pointing angle w.r.t. DStar decay vertex; [Cos(#theta)]; Entries"; numberOfBins = 200; lowerBound = -1; upperBound = 1;}
        else continue;} 
        if(k==27){if(i==0 || i==3 || i==4){name_Histogram = "d0MotherToDStar"; discription_Histogram = "d0 Mother w.r.t. DStar decay vertex; [cm]; Entries"; numberOfBins = 200; lowerBound = 0; upperBound = 1;}
        else continue;} 
        if(k==28){if(i==0 || i==3 || i==4){name_Histogram = "d0FirstDaughterToDStar"; discription_Histogram = "d0 first daughter w.r.t. DStar decay vertex; [cm]; Entries"; numberOfBins = 200; lowerBound = 0; upperBound = 1;}
        else continue;} 
        if(k==29){if(i==0 || i==3 || i==4){name_Histogram = "d0SecondDaughterToDStar"; discription_Histogram = "d0 second daughter w.r.t. DStar decay vertex; [cm]; Entries"; numberOfBins = 200; lowerBound = 0; upperBound = 1;}
        else continue;}
        if(k==30){if(i==0 || i==3 || i==4){name_Histogram = "impactProductToDStar"; discription_Histogram = "impact product w.r.t. DStar decay vertex; [cm]; Entries"; numberOfBins = 200; lowerBound = -0.02; upperBound = 0.02;}
        else continue;} 
        if(k==31){if(i==0 || i==3 || i==4){name_Histogram = "impactProductXYToDStar"; discription_Histogram = "impact product XY w.r.t. DStar decay vertex; [cm^{2}]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.5;}
        else continue;} 
        if(k==32){if(i==0 || i==3 || i==4){name_Histogram = "normDecayLengthToDStar"; discription_Histogram = "Normalized decay length w.r.t. DStar decay vertex; [cm]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 50;}
          else continue;} 
        if(k==33){if(i==0 || i==3 || i==4){name_Histogram = "pseudoProperDecayTimeToDStar"; discription_Histogram = "Pseudo Proper Decay Time w.r.t DStar vertex; [a.u.]; Entries"; numberOfBins = 1000; lowerBound = -1; upperBound = 1;}
          else continue;}  
        if(k==34){if(i==0 || i==3 || i==4){name_Histogram = "DecayTimeToDStar"; discription_Histogram = "Decay Time w.r.t DStar vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.00000001;}
          else continue;}  
        if(k==35){if(i==0 || i==3 || i==4){name_Histogram = "normDecayTimeToDStar"; discription_Histogram = "Normalized Decay Time w.r.t DStar vertex; [a.u.]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 0.00001;}
          else continue;} 

        if(k==36){name_Histogram = "topomaticFirstDaughter"; discription_Histogram = "topomatic d0 first daughter; [cm]; Entries"; numberOfBins = 200; lowerBound = 0; upperBound = 20;}
        if(k==37){name_Histogram = "topomaticSecondDaughter"; discription_Histogram = "topomatic d0 second daughter; [cm]; Entries"; numberOfBins = 200; lowerBound = 0; upperBound = 20;}
        if(k==38){name_Histogram = "topomaticMax"; discription_Histogram = "Max topomatic; [cm]; Entries"; numberOfBins = 200; lowerBound = 0; upperBound = 20;}
        if(k==39){name_Histogram = "topomaticMin"; discription_Histogram = "Min topomatic; [cm]; Entries"; numberOfBins = 200; lowerBound = 0; upperBound = 20;}
        if(k==40){name_Histogram = "pointingAngleMotherXY"; discription_Histogram = "pointing angle XY;  [Cos(#theta)]; Entries"; numberOfBins = 1000; lowerBound = -1; upperBound = 1;}
        if(k==41){name_Histogram = "vertexDistanceXY"; discription_Histogram = "vertex distance between mother and primary vertex XY; distance [cm]; Entries"; numberOfBins = 1000; lowerBound = 0; upperBound = 10;}
        if(k==42){name_Histogram = "normDecayLengthXY"; discription_Histogram = "Normalized decay length w.r.t primary vertex XY; [cm]; Entries"; numberOfBins = 100; lowerBound = 0; upperBound = 50;}

       
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
        fMotherHistogramArray[i][j][k] = histogram_Clone;
      }


      name_Histogram = "";
      discription_Histogram  = "";
      numberOfBins = 0;
      lowerBound = 0.0;
      upperBound = 0.0;
      numberOfBinsTwo = 0;
      lowerBoundTwo = 0.0;
      upperBoundTwo = 0.0;

      //we make the 2D histograms for the reconstructed particles
      Int_t nFirst = 0;
      Int_t nSecond = 1;
      Int_t nVariables = 10;
      Int_t nHistograms = nVariables * (nVariables - 1) / 2;

      TList * list2D = new TList();   
      list2D->SetOwner();   
      TString name2D = "2D_Histograms";   
      name2D += add_name;   
      list2D->SetName(name2D.Data());    
      listout->Add(list2D);
      
      for (Int_t k = 0; k < nHistograms; ++k)
      {
        numberOfBins = 50; numberOfBinsTwo = 50;
        if(nFirst==0){name_Histogram = "d0FirstDaughter"; discription_Histogram = "d0 first daughter [cm];"; lowerBound = 0; upperBound = 1;}
        if(nFirst==1){name_Histogram = "d0SecondDaughter"; discription_Histogram = "d0 second daughter [cm];"; lowerBound = 0; upperBound = 1;}
        if(nFirst==2){name_Histogram = "d0Mother"; discription_Histogram = "d0 mother [cm];"; lowerBound = 0; upperBound = 1;}
        if(nFirst==3){name_Histogram = "pointingAngleMother"; discription_Histogram = "pointing angle [Cos(#theta)];"; lowerBound = -1; upperBound = 1;}
        if(nFirst==4){name_Histogram = "impactProduct"; discription_Histogram = "impact product [cm^{2}];"; lowerBound = -0.01; upperBound = 0.01;}
        if(nFirst==5){name_Histogram = "impactProductXY"; discription_Histogram = "impact product XY [cm^{2}];"; lowerBound = 0; upperBound = 0.5;}
        if(nFirst==6){name_Histogram = "vertexDistance"; discription_Histogram = "vertex distance between mother and primary vertex  [cm];"; lowerBound = 0; upperBound = 1;}
        if(nFirst==7){name_Histogram = "normDecayLength"; discription_Histogram = "Normalized decay length w.r.t primary vertex [cm];"; lowerBound = 0; upperBound = 50;}
        if(nFirst==8){name_Histogram = "pointingAngleMotherXY"; discription_Histogram = "pointing angle XY [Cos(#theta)];"; lowerBound = -1; upperBound = 1;}
        if(nFirst==9){name_Histogram = "vertexDistanceXY"; discription_Histogram = "vertex distance between mother and primary vertex XY [cm];"; lowerBound = 0; upperBound = 1;}
        if(nFirst==10){name_Histogram = "normDecayLengthXY"; discription_Histogram = "Normalized decay length w.r.t primary vertex XY [cm];"; lowerBound = 0; upperBound = 50;}

        if(nSecond==0){name_Histogram += "d0FirstDaughter"; discription_Histogram += "d0 first daughter [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
        if(nSecond==1){name_Histogram += "d0SecondDaughter"; discription_Histogram += "d0 second daughter [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
        if(nSecond==2){name_Histogram += "d0Mother"; discription_Histogram += "d0 mother [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
        if(nSecond==3){name_Histogram += "pointingAngleMother"; discription_Histogram += "pointing angle [Cos(#theta)];"; lowerBoundTwo = -1; upperBoundTwo = 1;}
        if(nSecond==4){name_Histogram += "impactProduct"; discription_Histogram += "impact product [cm^{2}];"; lowerBoundTwo = -0.01; upperBoundTwo = 0.01;}
        if(nSecond==5){name_Histogram += "impactProductXY"; discription_Histogram += "impact product XY [cm^{2}];"; lowerBoundTwo = 0; upperBoundTwo = 0.5;}
        if(nSecond==6){name_Histogram += "vertexDistance"; discription_Histogram += "vertex distance between mother and primary vertex  [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
        if(nSecond==7){name_Histogram += "normDecayLength"; discription_Histogram += "Normalized decay length w.r.t primary vertex [cm];"; lowerBoundTwo = 0; upperBoundTwo = 50;}
        if(nSecond==8){name_Histogram += "_pointingAngleMotherXY"; discription_Histogram += "pointing angle XY [Cos(#theta)];"; lowerBoundTwo = -1; upperBoundTwo = 1;}
        if(nSecond==9){name_Histogram += "_vertexDistanceXY"; discription_Histogram += "vertex distance between mother and primary vertex XY [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
        if(nSecond==10){name_Histogram += "_normDecayLengthXY"; discription_Histogram += "Normalized decay length w.r.t primary vertex XY [cm];"; lowerBoundTwo = 0; upperBoundTwo = 50;}
  
        name_Histogram += add_name;
        TH2F* histogram = new TH2F(name_Histogram.Data(),discription_Histogram.Data(),numberOfBins,lowerBound,upperBound,numberOfBinsTwo,lowerBoundTwo,upperBoundTwo);
        histogram->Sumw2();
        if(j%2==0) histogram->SetLineColor(6);
        if(j%2==1) histogram->SetLineColor(4);
        histogram->SetMarkerStyle(20);
        histogram->SetMarkerSize(0.6);
        histogram->SetMarkerColor(6);
        TH2F* histogram_Clone = (TH2F*)histogram->Clone();
        list2D->Add(histogram_Clone);
        fMotherHistogramArray2D[i][j][k] = histogram_Clone;

        nSecond++;
        if(nSecond>nVariables)
        {
          nFirst++;
          nSecond = nFirst + 1;
        }
      }

      if(fUse3DHistograms)
      {
        name_Histogram = "";
        discription_Histogram  = "";
        numberOfBins = 0;
        lowerBound = 0.0;
        upperBound = 0.0;
        numberOfBinsTwo = 0;
        lowerBoundTwo = 0.0;
        upperBoundTwo = 0.0;
        Int_t numberOfBinsThree = 0;
        Int_t lowerBoundThree = 0.0;
        Int_t upperBoundThree = 0.0;

        //we make the 3D histograms for the reconstructed particles
        nFirst = 0;
        nSecond = 1;
        Int_t nThird = 2;
        nVariables = 10;
        nHistograms = nVariables * (nVariables - 1) * (nVariables - 2) / 6;

        TList * list3D = new TList();   
        list3D->SetOwner();   
        TString name3D = "3D_Histograms";   
        name3D += add_name;   
        list3D->SetName(name2D.Data());    
        listout->Add(list3D);
        
        for (Int_t k = 0; k < nHistograms; ++k)
        {
          numberOfBins = 50; numberOfBinsTwo = 50;
          if(nFirst==0){name_Histogram = "d0FirstDaughter"; discription_Histogram = "d0 first daughter [cm];"; lowerBound = 0; upperBound = 1;}
          if(nFirst==1){name_Histogram = "d0SecondDaughter"; discription_Histogram = "d0 second daughter [cm];"; lowerBound = 0; upperBound = 1;}
          if(nFirst==2){name_Histogram = "d0Mother"; discription_Histogram = "d0 mother [cm];"; lowerBound = 0; upperBound = 1;}
          if(nFirst==3){name_Histogram = "pointingAngleMother"; discription_Histogram = "pointing angle [Cos(#theta)];"; lowerBound = -1; upperBound = 1;}
          if(nFirst==4){name_Histogram = "impactProduct"; discription_Histogram = "impact product [cm^{2}];"; lowerBound = -0.01; upperBound = 0.01;}
          if(nFirst==5){name_Histogram = "impactProductXY"; discription_Histogram = "impact product XY [cm^{2}];"; lowerBound = 0; upperBound = 0.5;}
          if(nFirst==6){name_Histogram = "vertexDistance"; discription_Histogram = "vertex distance between mother and primary vertex  [cm];"; lowerBound = 0; upperBound = 1;}
          if(nFirst==7){name_Histogram = "normDecayLength"; discription_Histogram = "Normalized decay length w.r.t primary vertex [cm];"; lowerBound = 0; upperBound = 50;}
          if(nFirst==8){name_Histogram = "pointingAngleMotherXY"; discription_Histogram = "pointing angle XY [Cos(#theta)];"; lowerBound = -1; upperBound = 1;}
          if(nFirst==9){name_Histogram = "vertexDistanceXY"; discription_Histogram = "vertex distance between mother and primary vertex XY [cm];"; lowerBound = 0; upperBound = 1;}
          if(nFirst==10){name_Histogram = "normDecayLengthXY"; discription_Histogram = "Normalized decay length w.r.t primary vertex XY [cm];"; lowerBound = 0; upperBound = 50;}

          if(nSecond==0){name_Histogram += "_d0FirstDaughter"; discription_Histogram += "d0 first daughter [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
          if(nSecond==1){name_Histogram += "_d0SecondDaughter"; discription_Histogram += "d0 second daughter [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
          if(nSecond==2){name_Histogram += "_d0Mother"; discription_Histogram += "d0 mother [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
          if(nSecond==3){name_Histogram += "_pointingAngleMother"; discription_Histogram += "pointing angle [Cos(#theta)];"; lowerBoundTwo = -1; upperBoundTwo = 1;}
          if(nSecond==4){name_Histogram += "_impactProduct"; discription_Histogram += "impact product [cm^{2}];"; lowerBoundTwo = -0.01; upperBoundTwo = 0.01;}
          if(nSecond==5){name_Histogram += "_impactProductXY"; discription_Histogram += "impact product XY [cm^{2}];"; lowerBoundTwo = 0; upperBoundTwo = 0.5;}
          if(nSecond==6){name_Histogram += "_vertexDistance"; discription_Histogram += "vertex distance between mother and primary vertex  [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
          if(nSecond==7){name_Histogram += "_normDecayLength"; discription_Histogram += "Normalized decay length w.r.t primary vertex [cm];"; lowerBoundTwo = 0; upperBoundTwo = 50;}
          if(nSecond==8){name_Histogram += "_pointingAngleMotherXY"; discription_Histogram += "pointing angle XY [Cos(#theta)];"; lowerBoundTwo = -1; upperBoundTwo = 1;}
          if(nSecond==9){name_Histogram += "_vertexDistanceXY"; discription_Histogram += "vertex distance between mother and primary vertex XY [cm];"; lowerBoundTwo = 0; upperBoundTwo = 1;}
          if(nSecond==10){name_Histogram += "_normDecayLengthXY"; discription_Histogram += "Normalized decay length w.r.t primary vertex XY [cm];"; lowerBoundTwo = 0; upperBoundTwo = 50;}

          if(nThird==0){name_Histogram += "_d0FirstDaughter"; discription_Histogram += "d0 first daughter [cm];"; lowerBoundThree = 0; upperBoundThree = 1;}
          if(nThird==1){name_Histogram += "_d0SecondDaughter"; discription_Histogram += "d0 second daughter [cm];"; lowerBoundThree = 0; upperBoundThree = 1;}
          if(nThird==2){name_Histogram += "_d0Mother"; discription_Histogram += "d0 mother [cm];"; lowerBoundThree = 0; upperBoundThree = 1;}
          if(nThird==3){name_Histogram += "_pointingAngleMother"; discription_Histogram += "pointing angle [Cos(#theta)];"; lowerBoundThree = -1; upperBoundThree = 1;}
          if(nThird==4){name_Histogram += "_impactProduct"; discription_Histogram += "impact product [cm^{2}];"; lowerBoundThree = -0.01; upperBoundThree = 0.01;}
          if(nThird==5){name_Histogram += "_impactProductXY"; discription_Histogram += "impact product XY [cm^{2}];"; lowerBoundThree = 0; upperBoundThree = 0.5;}
          if(nThird==6){name_Histogram += "_vertexDistance"; discription_Histogram += "vertex distance between mother and primary vertex  [cm];"; lowerBoundThree = 0; upperBoundThree = 1;}
          if(nThird==7){name_Histogram += "_normDecayLength"; discription_Histogram += "Normalized decay length w.r.t primary vertex [cm];"; lowerBoundThree = 0; upperBoundThree = 50;}
          if(nThird==8){name_Histogram += "_pointingAngleMotherXY"; discription_Histogram += "pointing angle XY [Cos(#theta)];"; lowerBoundThree = -1; upperBoundThree = 1;}
          if(nThird==9){name_Histogram += "_vertexDistanceXY"; discription_Histogram += "vertex distance between mother and primary vertex XY [cm];"; lowerBoundThree = 0; upperBoundThree = 1;}
          if(nThird==10){name_Histogram += "_normDecayLengthXY"; discription_Histogram += "Normalized decay length w.r.t primary vertex XY [cm];"; lowerBoundThree = 0; upperBoundThree = 50;}


          name_Histogram += add_name;
          TH3F* histogram3D = new TH3F(name_Histogram.Data(),discription_Histogram.Data(),numberOfBins,lowerBound,upperBound,numberOfBinsTwo,lowerBoundTwo,upperBoundTwo,numberOfBinsThree,lowerBoundThree,upperBoundThree);
          histogram3D->Sumw2();
          if(j%2==0) histogram3D->SetLineColor(6);
          if(j%2==1) histogram3D->SetLineColor(4);
          histogram3D->SetMarkerStyle(20);
          histogram3D->SetMarkerSize(0.6);
          histogram3D->SetMarkerColor(6);
          TH3F* histogram_Clone3D = (TH3F*)histogram3D->Clone();
          list3D->Add(histogram_Clone3D);
          fMotherHistogramArray3D[i][j][k] = histogram_Clone3D;


          nThird++;
          if(nThird>nVariables)
          {
            nSecond++;
            nThird = nSecond + 1;
            if(nSecond>nVariables)
            {
              nFirst++;
              nSecond = nFirst + 1;
              nThird = nFirst + 2;

            }
          }
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
    fMotherHistogramArrayExtra[i][0] = effectOfCuts;

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
    fMotherHistogramArrayExtra[i][1] = effectOfCutsMC;

    TString name_particle_pdg ="particle_pdg";
    TH1F* hist_particle_pdg = new TH1F(name_particle_pdg.Data(),"Pdg code particle; pdg code; Entries",2000,-0.5,1999.5);
    hist_particle_pdg->Sumw2();
    hist_particle_pdg->SetLineColor(6);
    hist_particle_pdg->SetMarkerStyle(20);
    hist_particle_pdg->SetMarkerSize(0.6);
    hist_particle_pdg->SetMarkerColor(6);
    TH1F* histogram_particle_pdg = (TH1F*)hist_particle_pdg->Clone();
    listout->Add(histogram_particle_pdg);
    fMotherHistogramArrayExtra[i][2] = histogram_particle_pdg;

    TString name_particle_mother_pdg ="particle_mother_pdg";
    TH1F* hist_particle_mother_pdg = new TH1F(name_particle_mother_pdg.Data(),"Pdg code particle mother; pdg code; Entries",2000,-0.5,1999.5);
    hist_particle_mother_pdg->Sumw2();
    hist_particle_mother_pdg->SetLineColor(6);
    hist_particle_mother_pdg->SetMarkerStyle(20);
    hist_particle_mother_pdg->SetMarkerSize(0.6);
    hist_particle_mother_pdg->SetMarkerColor(6);
    TH1F* histogram_particle_mother_pdg = (TH1F*)hist_particle_mother_pdg->Clone();
    listout->Add(histogram_particle_mother_pdg);
    fMotherHistogramArrayExtra[i][3] = histogram_particle_mother_pdg;   

    TString name_distance_vertex_from_real ="distance_vertex_from_real";
    TH1F* hist_distance_vertex_from_real = new TH1F(name_distance_vertex_from_real.Data(),"Distance reconstructed vertex from real vertex; distance [cm]; Entries",100,0,1);
    hist_distance_vertex_from_real->Sumw2();
    hist_distance_vertex_from_real->SetLineColor(6);
    hist_distance_vertex_from_real->SetMarkerStyle(20);
    hist_distance_vertex_from_real->SetMarkerSize(0.6);
    hist_distance_vertex_from_real->SetMarkerColor(6);
    TH1F* histogram_distance_vertex_from_real = (TH1F*)hist_distance_vertex_from_real->Clone();
    listout->Add(histogram_distance_vertex_from_real);
    fMotherHistogramArrayExtra[i][4] = histogram_distance_vertex_from_real;

    TString name_distance_vertex_from_real_new ="distance_vertex_from_real_new";
    TH1F* hist_distance_vertex_from_real_new = new TH1F(name_distance_vertex_from_real_new.Data(),"Distance reconstructed vertex from real vertex; distance [cm]; Entries",100,0,1);
    hist_distance_vertex_from_real_new->Sumw2();
    hist_distance_vertex_from_real_new->SetLineColor(6);
    hist_distance_vertex_from_real_new->SetMarkerStyle(20);
    hist_distance_vertex_from_real_new->SetMarkerSize(0.6);
    hist_distance_vertex_from_real_new->SetMarkerColor(6);
    TH1F* histogram_distance_vertex_from_real_new = (TH1F*)hist_distance_vertex_from_real_new->Clone();
    listout->Add(histogram_distance_vertex_from_real_new);
    fMotherHistogramArrayExtra[i][5] = histogram_distance_vertex_from_real_new;

    TString name_momentum_resolution ="momentum_resolution";
    TH1F* hist_momentum_resolution = new TH1F(name_momentum_resolution.Data(),"Momentum resolution; difference between real and reconstructed momentum [GeV/c]; Entries",1000,0,1);
    hist_momentum_resolution->Sumw2();
    hist_momentum_resolution->SetLineColor(6);
    hist_momentum_resolution->SetMarkerStyle(20);
    hist_momentum_resolution->SetMarkerSize(0.6);
    hist_momentum_resolution->SetMarkerColor(6);
    TH1F* histogram_momentum_resolution = (TH1F*)hist_momentum_resolution->Clone();
    listout->Add(histogram_momentum_resolution);
    fMotherHistogramArrayExtra[i][6] = histogram_momentum_resolution;    
  }

  //we make the histograms for the same sign method histograms and the pt bins
  for (Int_t k = 0; k < fnPtBins+3; ++k){
    TString ptBinMother = "";
    if(k==0) ptBinMother = "";
    if(k==1) ptBinMother = "_ptbin_6_to_inf";
    if(k==2) ptBinMother = "_ptbin_3_to_inf";
    if(k>2) {ptBinMother += "_ptbin_"; ptBinMother += fPtBinLimits[k-3]; ptBinMother += "_to_"; ptBinMother += fPtBinLimits[k-2];}
  
    for (Int_t i = 0; i < 8; ++i){
      TString signName = "";
      if(i==0) signName = "";
      if(i==1) signName = "_SameSign";
      if(i==2) signName = "_SignSum";
      if(i==3) signName = "_HIJING_Background";
      if(i==4) signName = "_HIJING_Signal";
      if(i==5) signName = "_Background_rotation";
      if(i==6) signName = "_HIJING_Background_rotation";
      if(i==7) signName = "_correlated511";
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

      TString name_deltainvariantMassMother ="deltainvariantMassB0";
      name_deltainvariantMassMother += ptBinMother + signName;
      TH1F* hist_deltainvariantMassMother = new TH1F(name_deltainvariantMassMother.Data(),"delta mass mother candidate; m [GeV/c^2]; Entries",2000,0,20);
      hist_deltainvariantMassMother->Sumw2();
      hist_deltainvariantMassMother->SetLineColor(6);
      hist_deltainvariantMassMother->SetMarkerStyle(20);
      hist_deltainvariantMassMother->SetMarkerSize(0.6);
      hist_deltainvariantMassMother->SetMarkerColor(6);
      TH1F* histogram_deltainvariantMassMother = (TH1F*)hist_deltainvariantMassMother->Clone();
      fOutputB0MC->Add(histogram_deltainvariantMassMother);
    }
  }

  // for (Int_t k = 0; k < fnPtBins+3; ++k){
  //   TString ptBinMother = "";
  //   if(k==0) ptBinMother = "";
  //   if(k==1) ptBinMother = "_ptbin_6_to_inf";
  //   if(k==2) ptBinMother = "_ptbin_3_to_inf";
  //   if(k>2) {ptBinMother += "_ptbin_"; ptBinMother += fPtBinLimits[k-3]; ptBinMother += "_to_"; ptBinMother += fPtBinLimits[k-2];}
 
  //   for (Int_t i = 0; i < 7; ++i){
  //     TString signName = "";
  //     if(i==0) signName = "";
  //     if(i==1) signName = "_SameSign";
  //     if(i==2) signName = "_SignSum";
  //     if(i==3) signName = "_HIJING_Background";
  //     if(i==4) signName = "_HIJING_Signal";
  //     if(i==5) signName = "_Background_rotation";
  //     if(i==6) signName = "_HIJING_Background_rotation";

  //     TString name_invariantMassMother ="fineBin_invariantMassB0";
  //     name_invariantMassMother += ptBinMother + signName;
  //     TH1F* hist_invariantMassMother = new TH1F(name_invariantMassMother.Data(),"mass mother candidate; m [GeV/c^2]; Entries",5000,2.5,7.5);
  //     hist_invariantMassMother->Sumw2();
  //     hist_invariantMassMother->SetLineColor(6);
  //     hist_invariantMassMother->SetMarkerStyle(20);
  //     hist_invariantMassMother->SetMarkerSize(0.6);
  //     hist_invariantMassMother->SetMarkerColor(6);
  //     TH1F* histogram_invariantMassMother = (TH1F*)hist_invariantMassMother->Clone();
  //     fOutputB0MC->Add(histogram_invariantMassMother);

  //     TString name_deltainvariantMassMother ="fineBin_deltainvariantMassB0";
  //     name_deltainvariantMassMother += ptBinMother + signName;
  //     TH1F* hist_deltainvariantMassMother = new TH1F(name_deltainvariantMassMother.Data(),"delta mass mother candidate; m [GeV/c^2]; Entries",2000,0,20);
  //     hist_deltainvariantMassMother->Sumw2();
  //     hist_deltainvariantMassMother->SetLineColor(6);
  //     hist_deltainvariantMassMother->SetMarkerStyle(20);
  //     hist_deltainvariantMassMother->SetMarkerSize(0.6);
  //     hist_deltainvariantMassMother->SetMarkerColor(6);
  //     TH1F* histogram_deltainvariantMassMother = (TH1F*)hist_deltainvariantMassMother->Clone();
  //     fOutputB0MC->Add(histogram_deltainvariantMassMother);
  //   }
  // }

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

  TString name_particle_pdgB0Pion ="particle_pdgB0Pion";
  TH2F* hist_particle_pdgB0Pion = new TH2F(name_particle_pdgB0Pion.Data(),"Pdg code particle; pdg code; Entries",5000,-0.5,4999.5,10,0,10);
  hist_particle_pdgB0Pion->Sumw2();
  hist_particle_pdgB0Pion->SetLineColor(6);
  hist_particle_pdgB0Pion->SetMarkerStyle(20);
  hist_particle_pdgB0Pion->SetMarkerSize(0.6);
  hist_particle_pdgB0Pion->SetMarkerColor(6);
  TH2F* histogram_particle_pdgB0Pion = (TH2F*)hist_particle_pdgB0Pion->Clone();
  fOutputB0MC->Add(histogram_particle_pdgB0Pion); 

  TString name_particle_pdgDStarPion ="particle_pdgDStarPion";
  TH2F* hist_particle_pdgDStarPion = new TH2F(name_particle_pdgDStarPion.Data(),"Pdg code particle; pdg code; Entries",5000,-0.5,4999.5,10,0,10);
  hist_particle_pdgDStarPion->Sumw2();
  hist_particle_pdgDStarPion->SetLineColor(6);
  hist_particle_pdgDStarPion->SetMarkerStyle(20);
  hist_particle_pdgDStarPion->SetMarkerSize(0.6);
  hist_particle_pdgDStarPion->SetMarkerColor(6);
  TH2F* histogram_particle_pdgDStarPion = (TH2F*)hist_particle_pdgDStarPion->Clone();
  fOutputB0MC->Add(histogram_particle_pdgDStarPion); 

  TString name_particle_pdgD0First ="particle_pdgD0First";
  TH2F* hist_particle_pdgD0First = new TH2F(name_particle_pdgD0First.Data(),"Pdg code particle; pdg code; Entries",5000,-0.5,4999.5,10,0,10);
  hist_particle_pdgD0First->Sumw2();
  hist_particle_pdgD0First->SetLineColor(6);
  hist_particle_pdgD0First->SetMarkerStyle(20);
  hist_particle_pdgD0First->SetMarkerSize(0.6);
  hist_particle_pdgD0First->SetMarkerColor(6);
  TH2F* histogram_particle_pdgD0First = (TH2F*)hist_particle_pdgD0First->Clone();
  fOutputB0MC->Add(histogram_particle_pdgD0First); 

  TString name_particle_pdgD0Second ="particle_pdgD0Second";
  TH2F* hist_particle_pdgD0Second = new TH2F(name_particle_pdgD0Second.Data(),"Pdg code particle; pdg code; Entries",5000,-0.5,4999.5,10,0,10);
  hist_particle_pdgD0Second->Sumw2();
  hist_particle_pdgD0Second->SetLineColor(6);
  hist_particle_pdgD0Second->SetMarkerStyle(20);
  hist_particle_pdgD0Second->SetMarkerSize(0.6);
  hist_particle_pdgD0Second->SetMarkerColor(6);
  TH2F* histogram_particle_pdgD0Second = (TH2F*)hist_particle_pdgD0Second->Clone();
  fOutputB0MC->Add(histogram_particle_pdgD0Second); 

  TString name_particle_pdgAll ="particle_pdgAll";
  TH1F* hist_particle_pdgAll = new TH1F(name_particle_pdgAll.Data(),"Pdg code particle; pdg code; Entries",5000,-0.5,4999.5);
  hist_particle_pdgAll->Sumw2();
  hist_particle_pdgAll->SetLineColor(6);
  hist_particle_pdgAll->SetMarkerStyle(20);
  hist_particle_pdgAll->SetMarkerSize(0.6);
  hist_particle_pdgAll->SetMarkerColor(6);
  TH1F* histogram_particle_pdgAll = (TH1F*)hist_particle_pdgAll->Clone();
  fOutputB0MC->Add(histogram_particle_pdgAll); 

  TString name_particle_pdgAllSecond ="particle_pdgAllSecond";
  TH1F* hist_particle_pdgAllSecond = new TH1F(name_particle_pdgAllSecond.Data(),"Pdg code particle; pdg code; Entries",5000,-0.5,4999.5);
  hist_particle_pdgAllSecond->Sumw2();
  hist_particle_pdgAllSecond->SetLineColor(6);
  hist_particle_pdgAllSecond->SetMarkerStyle(20);
  hist_particle_pdgAllSecond->SetMarkerSize(0.6);
  hist_particle_pdgAllSecond->SetMarkerColor(6);
  TH1F* histogram_particle_pdgAllSecond = (TH1F*)hist_particle_pdgAllSecond->Clone();
  fOutputB0MC->Add(histogram_particle_pdgAllSecond); 

  TString name_particle_pdgAllThird ="particle_pdgAllThird";
  TH1F* hist_particle_pdgAllThird = new TH1F(name_particle_pdgAllThird.Data(),"Pdg code particle; pdg code; Entries",5000,-0.5,4999.5);
  hist_particle_pdgAllThird->Sumw2();
  hist_particle_pdgAllThird->SetLineColor(6);
  hist_particle_pdgAllThird->SetMarkerStyle(20);
  hist_particle_pdgAllThird->SetMarkerSize(0.6);
  hist_particle_pdgAllThird->SetMarkerColor(6);
  TH1F* histogram_particle_pdgAllThird = (TH1F*)hist_particle_pdgAllThird->Clone();
  fOutputB0MC->Add(histogram_particle_pdgAllThird); 

  TString name_particle_pdgAllFourth ="particle_pdgAllFourth";
  TH1F* hist_particle_pdgAllFourth = new TH1F(name_particle_pdgAllFourth.Data(),"Pdg code particle; pdg code; Entries",5000,-0.5,4999.5);
  hist_particle_pdgAllFourth->Sumw2();
  hist_particle_pdgAllFourth->SetLineColor(6);
  hist_particle_pdgAllFourth->SetMarkerStyle(20);
  hist_particle_pdgAllFourth->SetMarkerSize(0.6);
  hist_particle_pdgAllFourth->SetMarkerColor(6);
  TH1F* histogram_particle_pdgAllFourth = (TH1F*)hist_particle_pdgAllFourth->Clone();
  fOutputB0MC->Add(histogram_particle_pdgAllFourth); 

  TString name_particle_pdgAllFifth ="particle_pdgAllFifth";
  TH1F* hist_particle_pdgAllFifth = new TH1F(name_particle_pdgAllFifth.Data(),"Pdg code particle; pdg code; Entries",5000,-0.5,4999.5);
  hist_particle_pdgAllFifth->Sumw2();
  hist_particle_pdgAllFifth->SetLineColor(6);
  hist_particle_pdgAllFifth->SetMarkerStyle(20);
  hist_particle_pdgAllFifth->SetMarkerSize(0.6);
  hist_particle_pdgAllFifth->SetMarkerColor(6);
  TH1F* histogram_particle_pdgAllFifth = (TH1F*)hist_particle_pdgAllFifth->Clone();
  fOutputB0MC->Add(histogram_particle_pdgAllFifth); 

  TString name_particle_pdgAllInvMass ="particle_pdgAllInvMass";
  TH2F* hist_particle_pdgAllInvMass = new TH2F(name_particle_pdgAllInvMass.Data(),"Pdg code particle; pdg code; B0 candidate Inv. Mass",5000,-0.5,4999.5,500,4.0,6.0);
  hist_particle_pdgAllInvMass->Sumw2();
  hist_particle_pdgAllInvMass->SetLineColor(6);
  hist_particle_pdgAllInvMass->SetMarkerStyle(20);
  hist_particle_pdgAllInvMass->SetMarkerSize(0.6);
  hist_particle_pdgAllInvMass->SetMarkerColor(6);
  TH2F* histogram_particle_pdgAllInvMass = (TH2F*)hist_particle_pdgAllInvMass->Clone();
  fOutputB0MC->Add(histogram_particle_pdgAllInvMass); 

  TString name_particle_pdgAllInvMassSecond ="particle_pdgAllInvMassSecond";
  TH2F* hist_particle_pdgAllInvMassSecond = new TH2F(name_particle_pdgAllInvMassSecond.Data(),"Pdg code particle; pdg code; B0 candidate Inv. Mass",5000,-0.5,4999.5,500,4.0,6.0);
  hist_particle_pdgAllInvMassSecond->Sumw2();
  hist_particle_pdgAllInvMassSecond->SetLineColor(6);
  hist_particle_pdgAllInvMassSecond->SetMarkerStyle(20);
  hist_particle_pdgAllInvMassSecond->SetMarkerSize(0.6);
  hist_particle_pdgAllInvMassSecond->SetMarkerColor(6);
  TH2F* histogram_particle_pdgAllInvMassSecond = (TH2F*)hist_particle_pdgAllInvMassSecond->Clone();
  fOutputB0MC->Add(histogram_particle_pdgAllInvMassSecond);  

  TString name_particle_pdgAllInvMassThird ="particle_pdgAllInvMassThird";
  TH2F* hist_particle_pdgAllInvMassThird = new TH2F(name_particle_pdgAllInvMassThird.Data(),"Pdg code particle; pdg code; B0 candidate Inv. Mass",5000,-0.5,4999.5,500,4.0,6.0);
  hist_particle_pdgAllInvMassThird->Sumw2();
  hist_particle_pdgAllInvMassThird->SetLineColor(6);
  hist_particle_pdgAllInvMassThird->SetMarkerStyle(20);
  hist_particle_pdgAllInvMassThird->SetMarkerSize(0.6);
  hist_particle_pdgAllInvMassThird->SetMarkerColor(6);
  TH2F* histogram_particle_pdgAllInvMassThird = (TH2F*)hist_particle_pdgAllInvMassThird->Clone();
  fOutputB0MC->Add(histogram_particle_pdgAllInvMassThird);  

  TString name_particle_pdgAllInvMassFourth ="particle_pdgAllInvMassFourth";
  TH2F* hist_particle_pdgAllInvMassFourth = new TH2F(name_particle_pdgAllInvMassFourth.Data(),"Pdg code particle; pdg code; B0 candidate Inv. Mass",5000,-0.5,4999.5,500,4.0,6.0);
  hist_particle_pdgAllInvMassFourth->Sumw2();
  hist_particle_pdgAllInvMassFourth->SetLineColor(6);
  hist_particle_pdgAllInvMassFourth->SetMarkerStyle(20);
  hist_particle_pdgAllInvMassFourth->SetMarkerSize(0.6);
  hist_particle_pdgAllInvMassFourth->SetMarkerColor(6);
  TH2F* histogram_particle_pdgAllInvMassFourth = (TH2F*)hist_particle_pdgAllInvMassFourth->Clone();
  fOutputB0MC->Add(histogram_particle_pdgAllInvMassFourth);  

  TString name_particle_pdgAllInvMassFifth ="particle_pdgAllInvMassFifth";
  TH2F* hist_particle_pdgAllInvMassFifth = new TH2F(name_particle_pdgAllInvMassFifth.Data(),"Pdg code particle; pdg code; B0 candidate Inv. Mass",5000,-0.5,4999.5,500,4.0,6.0);
  hist_particle_pdgAllInvMassFifth->Sumw2();
  hist_particle_pdgAllInvMassFifth->SetLineColor(6);
  hist_particle_pdgAllInvMassFifth->SetMarkerStyle(20);
  hist_particle_pdgAllInvMassFifth->SetMarkerSize(0.6);
  hist_particle_pdgAllInvMassFifth->SetMarkerColor(6);
  TH2F* histogram_particle_pdgAllInvMassFifth = (TH2F*)hist_particle_pdgAllInvMassFifth->Clone();
  fOutputB0MC->Add(histogram_particle_pdgAllInvMassFifth);  

  TString name_particle_daughterPdgOneStep511a ="particle_daughterPdgOneStep511a";
  TH2F* hist_particle_daughterPdgOneStep511a = new TH2F(name_particle_daughterPdgOneStep511a.Data(),"Pdg daughters; n daughter; Pdg daughter",50,-0.5,49.5,5000,-0.5,4999.5);
  hist_particle_daughterPdgOneStep511a->Sumw2();
  hist_particle_daughterPdgOneStep511a->SetLineColor(6);
  hist_particle_daughterPdgOneStep511a->SetMarkerStyle(20);
  hist_particle_daughterPdgOneStep511a->SetMarkerSize(0.6);
  hist_particle_daughterPdgOneStep511a->SetMarkerColor(6);
  TH2F* histogram_particle_daughterPdgOneStep511a = (TH2F*)hist_particle_daughterPdgOneStep511a->Clone();
  fOutputB0MC->Add(histogram_particle_daughterPdgOneStep511a); 

  TString name_particle_daughterPdgOneStep521a ="particle_daughterPdgOneStep521a";
  TH2F* hist_particle_daughterPdgOneStep521a = new TH2F(name_particle_daughterPdgOneStep521a.Data(),"Pdg daughters; n daughter; Pdg daughter",50,-0.5,49.5,5000,-0.5,4999.5);
  hist_particle_daughterPdgOneStep521a->Sumw2();
  hist_particle_daughterPdgOneStep521a->SetLineColor(6);
  hist_particle_daughterPdgOneStep521a->SetMarkerStyle(20);
  hist_particle_daughterPdgOneStep521a->SetMarkerSize(0.6);
  hist_particle_daughterPdgOneStep521a->SetMarkerColor(6);
  TH2F* histogram_particle_daughterPdgOneStep521a = (TH2F*)hist_particle_daughterPdgOneStep521a->Clone();
  fOutputB0MC->Add(histogram_particle_daughterPdgOneStep521a); 

  TString name_particle_daughterPdgOneStep511b ="particle_daughterPdgOneStep511b";
  TH2F* hist_particle_daughterPdgOneStep511b = new TH2F(name_particle_daughterPdgOneStep511b.Data(),"Pdg daughters; n daughter; Pdg daughter",50,-0.5,49.5,5000,-0.5,4999.5);
  hist_particle_daughterPdgOneStep511b->Sumw2();
  hist_particle_daughterPdgOneStep511b->SetLineColor(6);
  hist_particle_daughterPdgOneStep511b->SetMarkerStyle(20);
  hist_particle_daughterPdgOneStep511b->SetMarkerSize(0.6);
  hist_particle_daughterPdgOneStep511b->SetMarkerColor(6);
  TH2F* histogram_particle_daughterPdgOneStep511b = (TH2F*)hist_particle_daughterPdgOneStep511b->Clone();
  fOutputB0MC->Add(histogram_particle_daughterPdgOneStep511b); 

  TString name_particle_daughterPdgOneStep521b ="particle_daughterPdgOneStep521b";
  TH2F* hist_particle_daughterPdgOneStep521b = new TH2F(name_particle_daughterPdgOneStep521b.Data(),"Pdg daughters; n daughter; Pdg daughter",50,-0.5,49.5,5000,-0.5,4999.5);
  hist_particle_daughterPdgOneStep521b->Sumw2();
  hist_particle_daughterPdgOneStep521b->SetLineColor(6);
  hist_particle_daughterPdgOneStep521b->SetMarkerStyle(20);
  hist_particle_daughterPdgOneStep521b->SetMarkerSize(0.6);
  hist_particle_daughterPdgOneStep521b->SetMarkerColor(6);
  TH2F* histogram_particle_daughterPdgOneStep521b = (TH2F*)hist_particle_daughterPdgOneStep521b->Clone();
  fOutputB0MC->Add(histogram_particle_daughterPdgOneStep521b); 

  TString name_particle_daughterPdgTwoStep511a ="particle_daughterPdgTwoStep511a";
  TH2F* hist_particle_daughterPdgTwoStep511a = new TH2F(name_particle_daughterPdgTwoStep511a.Data(),"Pdg daughters; n daughter; Pdg daughter",50,-0.5,49.5,5000,-0.5,4999.5);
  hist_particle_daughterPdgTwoStep511a->Sumw2();
  hist_particle_daughterPdgTwoStep511a->SetLineColor(6);
  hist_particle_daughterPdgTwoStep511a->SetMarkerStyle(20);
  hist_particle_daughterPdgTwoStep511a->SetMarkerSize(0.6);
  hist_particle_daughterPdgTwoStep511a->SetMarkerColor(6);
  TH2F* histogram_particle_daughterPdgTwoStep511a = (TH2F*)hist_particle_daughterPdgTwoStep511a->Clone();
  fOutputB0MC->Add(histogram_particle_daughterPdgTwoStep511a); 

  TString name_particle_daughterPdgTwoStep521a ="particle_daughterPdgTwoStep521a";
  TH2F* hist_particle_daughterPdgTwoStep521a = new TH2F(name_particle_daughterPdgTwoStep521a.Data(),"Pdg daughters; n daughter; Pdg daughter",50,-0.5,49.5,5000,-0.5,4999.5);
  hist_particle_daughterPdgTwoStep521a->Sumw2();
  hist_particle_daughterPdgTwoStep521a->SetLineColor(6);
  hist_particle_daughterPdgTwoStep521a->SetMarkerStyle(20);
  hist_particle_daughterPdgTwoStep521a->SetMarkerSize(0.6);
  hist_particle_daughterPdgTwoStep521a->SetMarkerColor(6);
  TH2F* histogram_particle_daughterPdgTwoStep521a = (TH2F*)hist_particle_daughterPdgTwoStep521a->Clone();
  fOutputB0MC->Add(histogram_particle_daughterPdgTwoStep521a); 

  TString name_particle_daughterPdgTwoStep511b ="particle_daughterPdgTwoStep511b";
  TH2F* hist_particle_daughterPdgTwoStep511b = new TH2F(name_particle_daughterPdgTwoStep511b.Data(),"Pdg daughters; n daughter; Pdg daughter",50,-0.5,49.5,5000,-0.5,4999.5);
  hist_particle_daughterPdgTwoStep511b->Sumw2();
  hist_particle_daughterPdgTwoStep511b->SetLineColor(6);
  hist_particle_daughterPdgTwoStep511b->SetMarkerStyle(20);
  hist_particle_daughterPdgTwoStep511b->SetMarkerSize(0.6);
  hist_particle_daughterPdgTwoStep511b->SetMarkerColor(6);
  TH2F* histogram_particle_daughterPdgTwoStep511b = (TH2F*)hist_particle_daughterPdgTwoStep511b->Clone();
  fOutputB0MC->Add(histogram_particle_daughterPdgTwoStep511b); 

  TString name_particle_daughterPdgTwoStep521b ="particle_daughterPdgTwoStep521b";
  TH2F* hist_particle_daughterPdgTwoStep521b = new TH2F(name_particle_daughterPdgTwoStep521b.Data(),"Pdg daughters; n daughter; Pdg daughter",50,-0.5,49.5,5000,-0.5,4999.5);
  hist_particle_daughterPdgTwoStep521b->Sumw2();
  hist_particle_daughterPdgTwoStep521b->SetLineColor(6);
  hist_particle_daughterPdgTwoStep521b->SetMarkerStyle(20);
  hist_particle_daughterPdgTwoStep521b->SetMarkerSize(0.6);
  hist_particle_daughterPdgTwoStep521b->SetMarkerColor(6);
  TH2F* histogram_particle_daughterPdgTwoStep521b = (TH2F*)hist_particle_daughterPdgTwoStep521b->Clone();
  fOutputB0MC->Add(histogram_particle_daughterPdgTwoStep521b); 

  TString name_particle_daughterPdgTwoStep511c ="particle_daughterPdgTwoStep511c";
  TH2F* hist_particle_daughterPdgTwoStep511c = new TH2F(name_particle_daughterPdgTwoStep511c.Data(),"Pdg daughters; n daughter; Pdg daughter",50,-0.5,49.5,5000,-0.5,4999.5);
  hist_particle_daughterPdgTwoStep511c->Sumw2();
  hist_particle_daughterPdgTwoStep511c->SetLineColor(6);
  hist_particle_daughterPdgTwoStep511c->SetMarkerStyle(20);
  hist_particle_daughterPdgTwoStep511c->SetMarkerSize(0.6);
  hist_particle_daughterPdgTwoStep511c->SetMarkerColor(6);
  TH2F* histogram_particle_daughterPdgTwoStep511c = (TH2F*)hist_particle_daughterPdgTwoStep511c->Clone();
  fOutputB0MC->Add(histogram_particle_daughterPdgTwoStep511c); 

  TString name_particle_daughterPdgTwoStep521c ="particle_daughterPdgTwoStep521c";
  TH2F* hist_particle_daughterPdgTwoStep521c = new TH2F(name_particle_daughterPdgTwoStep521c.Data(),"Pdg daughters; n daughter; Pdg daughter",50,-0.5,49.5,5000,-0.5,4999.5);
  hist_particle_daughterPdgTwoStep521c->Sumw2();
  hist_particle_daughterPdgTwoStep521c->SetLineColor(6);
  hist_particle_daughterPdgTwoStep521c->SetMarkerStyle(20);
  hist_particle_daughterPdgTwoStep521c->SetMarkerSize(0.6);
  hist_particle_daughterPdgTwoStep521c->SetMarkerColor(6);
  TH2F* histogram_particle_daughterPdgTwoStep521c = (TH2F*)hist_particle_daughterPdgTwoStep521c->Clone();
  fOutputB0MC->Add(histogram_particle_daughterPdgTwoStep521c); 

  TString name_invariantMassB0Signal_BA ="invariantMassB0Signal_BA";
  TH1F* hist_invariantMassB0Signal_BA = new TH1F(name_invariantMassB0Signal_BA.Data(),"mass mother candidate; m [GeV/c^2]; Entries",2000,0,20);
  hist_invariantMassB0Signal_BA->Sumw2();
  hist_invariantMassB0Signal_BA->SetLineColor(6);
  hist_invariantMassB0Signal_BA->SetMarkerStyle(20);
  hist_invariantMassB0Signal_BA->SetMarkerSize(0.6);
  hist_invariantMassB0Signal_BA->SetMarkerColor(6);
  TH1F* histogram_invariantMassB0Signal_BA = (TH1F*)hist_invariantMassB0Signal_BA->Clone();
  fOutputB0MC->Add(histogram_invariantMassB0Signal_BA);

  TString name_invariantMassB0Correlated_BA ="invariantMassB0Correlated_BA";
  TH1F* hist_invariantMassB0Correlated_BA = new TH1F(name_invariantMassB0Correlated_BA.Data(),"mass mother candidate; m [GeV/c^2]; Entries",2000,0,20);
  hist_invariantMassB0Correlated_BA->Sumw2();
  hist_invariantMassB0Correlated_BA->SetLineColor(6);
  hist_invariantMassB0Correlated_BA->SetMarkerStyle(20);
  hist_invariantMassB0Correlated_BA->SetMarkerSize(0.6);
  hist_invariantMassB0Correlated_BA->SetMarkerColor(6);
  TH1F* histogram_invariantMassB0Correlated_BA = (TH1F*)hist_invariantMassB0Correlated_BA->Clone();
  fOutputB0MC->Add(histogram_invariantMassB0Correlated_BA);

  TString name_invariantMassB0Background_BA ="invariantMassB0Background_BA";
  TH1F* hist_invariantMassB0Background_BA = new TH1F(name_invariantMassB0Background_BA.Data(),"mass mother candidate; m [GeV/c^2]; Entries",2000,0,20);
  hist_invariantMassB0Background_BA->Sumw2();
  hist_invariantMassB0Background_BA->SetLineColor(6);
  hist_invariantMassB0Background_BA->SetMarkerStyle(20);
  hist_invariantMassB0Background_BA->SetMarkerSize(0.6);
  hist_invariantMassB0Background_BA->SetMarkerColor(6);
  TH1F* histogram_invariantMassB0Background_BA = (TH1F*)hist_invariantMassB0Background_BA->Clone();
  fOutputB0MC->Add(histogram_invariantMassB0Background_BA);

  return;
}
//-------------------------------------------------------------------------------------
AliAODVertex* AliAnalysisTaskSEB0toDStarPi::RecalculateVertex(const AliVVertex *primary,TObjArray *tracks,Double_t bField, Double_t dispersion){
  //
  // Helper function to recalculate a vertex.
  //
  
  AliESDVertex *vertexESD = 0;
  AliAODVertex *vertexAOD = 0;

  AliVertexerTracks vertexer;
  vertexer.SetFieldkG(bField);

  vertexer.SetVtxStart((AliESDVertex*)primary); //primary vertex
  vertexESD = (AliESDVertex*)vertexer.VertexForSelectedESDTracks(tracks);

  // delete vertexer; vertexer=NULL;

  if(!vertexESD) return vertexAOD;


  if(vertexESD->GetNContributors()!=tracks->GetEntriesFast()) 
  {
    delete vertexESD; vertexESD=nullptr;
    return vertexAOD;
  }

  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  for(Int_t a=0;a<3;a++)pos[a]=0.;
  for(Int_t b=0;b<6;b++)cov[b]=0.;
  chi2perNDF=0;

  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix


  Double_t vertRadius2=pos[0]*pos[0]+pos[1]*pos[1];
  if(vertRadius2>8.) //(2.82)^2 radius beam pipe
  {
    delete vertexESD; vertexESD=nullptr;
    return vertexAOD;
  }
  
  chi2perNDF = vertexESD->GetChi2toNDF();
  dispersion = vertexESD->GetDispersion();
  delete vertexESD; vertexESD=nullptr;
  Int_t nprongs = 2; //tracks->GetEntriesFast();
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
    Double_t yMC[7] = {0.0};
    Double_t pseudoYMC[7] = {0.0};


    Bool_t mcPionB0Present = kFALSE;
    Bool_t mcPionDStarPresent = kFALSE;
    Bool_t mcPionD0Present = kFALSE;
    Bool_t mcKaonPresent = kFALSE;

    // Below, we find all the MC labels for the true signal tracks
    AliAODMCParticle *mcTrackParticle = dynamic_cast< AliAODMCParticle*>(mcTrackArray->At(i));
    if(!mcTrackParticle) {std::cout << "no particle" << std::endl; continue;}
    Int_t pdgCodeMC=TMath::Abs(mcTrackParticle->GetPdgCode());
    
    if (pdgCodeMC==511){ //if the track is a B0 we look at its daughters

      mcLabelB0 = i;
      Int_t nDaughterB0 = mcTrackParticle->GetNDaughters();
      ptMC[0] = mcTrackParticle->Pt();
      yMC[0] = mcTrackParticle->Y();
      pseudoYMC[0] = mcTrackParticle->Eta();

      TString fillthis= "B0s_in_analysis";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(0);

      if(nDaughterB0==2){
        for(Int_t iDaughterB0=0; iDaughterB0<2; iDaughterB0++){

          AliAODMCParticle* daughterB0 = (AliAODMCParticle*)mcTrackArray->At(mcTrackParticle->GetDaughterLabel(iDaughterB0));
          if(!daughterB0) break;
          Int_t pdgCodeDaughterB0=TMath::Abs(daughterB0->GetPdgCode());

          if (pdgCodeDaughterB0==211){ //if the track is a pion we save its monte carlo label
            mcLabelPionB0 = mcTrackParticle->GetDaughterLabel(iDaughterB0);
            mcPionB0Present = kTRUE;
            ptMC[1] = daughterB0->Pt();
            yMC[1] = daughterB0->Y();
            pseudoYMC[1] = daughterB0->Eta();

          } else if (pdgCodeDaughterB0==413){ //if the track is a DStar we look at its daughters
            mcLabelDStar = mcTrackParticle->GetDaughterLabel(iDaughterB0);
            Int_t nDaughterDStar = daughterB0->GetNDaughters();
            ptMC[2] = daughterB0->Pt();
            yMC[2] = daughterB0->Y();
            pseudoYMC[2] = daughterB0->Eta();

            if(nDaughterDStar==2){
              for(Int_t iDaughterDStar=0; iDaughterDStar<2; iDaughterDStar++){

                AliAODMCParticle* daughterDStar = (AliAODMCParticle*)mcTrackArray->At(daughterB0->GetDaughterLabel(iDaughterDStar));
                if(!daughterDStar) break;
                Int_t pdgCodeDaughterDStar=TMath::Abs(daughterDStar->GetPdgCode());

                if (pdgCodeDaughterDStar==211){ //if the track is a pion we save its monte carlo label
                  mcLabelPionDStar = daughterB0->GetDaughterLabel(iDaughterDStar);
                  mcPionDStarPresent = kTRUE;
                  ptMC[3] = daughterDStar->Pt();
                  yMC[3] = daughterDStar->Y();
                  pseudoYMC[3] = daughterDStar->Eta();

                } else if (pdgCodeDaughterDStar==421){ //if the track is a D0 we look at its daughters
                  mcLabelD0 = daughterB0->GetDaughterLabel(iDaughterDStar);
                  Int_t nDaughterD0 = daughterDStar->GetNDaughters();
                  ptMC[4] = daughterDStar->Pt();
                  yMC[4] = daughterDStar->Y();
                  pseudoYMC[4] = daughterDStar->Eta();

                  if(nDaughterD0==2){
                    for(Int_t iDaughterD0=0; iDaughterD0<2; iDaughterD0++){

                      AliAODMCParticle* daughterD0 = (AliAODMCParticle*)mcTrackArray->At(daughterDStar->GetDaughterLabel(iDaughterD0));
                      if(!daughterD0) break;
                      Int_t pdgCodeDaughterD0=TMath::Abs(daughterD0->GetPdgCode());

                      if (pdgCodeDaughterD0==211){ //if the track is a pion we save its monte carlo label
                        mcLabelPionD0 = daughterDStar->GetDaughterLabel(iDaughterD0);
                        ptMC[5] = daughterD0->Pt();
                        yMC[5] = daughterD0->Y();
                        pseudoYMC[5] = daughterD0->Eta();
                        mcPionD0Present = kTRUE;

                      } else if (pdgCodeDaughterD0==321){ //if the track is a kaon we save its monte carlo label
                        mcLabelKaon = daughterDStar->GetDaughterLabel(iDaughterD0);;
                        mcKaonPresent = kTRUE;
                        ptMC[6] = daughterD0->Pt();
                        yMC[6] = daughterD0->Y();
                        pseudoYMC[6] = daughterD0->Eta();

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
    if(mcPionB0Present && mcPionDStarPresent && mcPionD0Present && mcKaonPresent){

      // We also save information on the amount of signal tracks that exist in the MC dataset
      TString fillthis= "B0s_in_analysis";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(1);

      fillthis= "B0s_per_bin";
      for (Int_t j = 0; j < fnPtBins; ++j)
      {
        if(fPtBinLimits[j] < ptMC[0] && ptMC[0] < fPtBinLimits[j+1]) {((TH1F*)(listout->FindObject(fillthis)))->Fill(j); break;}
      }

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

      fillthis= "mc_B0_rapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(yMC[0]);
      fillthis= "mc_B0_pion_rapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(yMC[1]);
      fillthis= "mc_DStar_rapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(yMC[2]);
      fillthis= "mc_DStar_pion_rapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(yMC[3]);
      fillthis= "mc_D0_rapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(yMC[4]);
      fillthis= "mc_D0_pion_rapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(yMC[5]);
      fillthis= "mc_D0_kaon_rapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(yMC[6]);

      fillthis= "mc_B0_pseudorapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(pseudoYMC[0]);
      fillthis= "mc_B0_pion_pseudorapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(pseudoYMC[1]);
      fillthis= "mc_DStar_pseudorapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(pseudoYMC[2]);
      fillthis= "mc_DStar_pion_pseudorapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(pseudoYMC[3]);
      fillthis= "mc_D0_pseudorapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(pseudoYMC[4]);
      fillthis= "mc_D0_pion_pseudorapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(pseudoYMC[5]);
      fillthis= "mc_D0_kaon_pseudorapidity_true";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(pseudoYMC[6]);

      fillthis= "mc_B0_pt_bins";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(ptMC[0]);
      if(TMath::Abs(yMC[0]) < 0.5) 
      {
        fillthis= "mc_B0_pt_bins_lim_acc";
        ((TH1F*)(listout->FindObject(fillthis)))->Fill(ptMC[0]);
      }

      if(TMath::Abs(yMC[0]) < 0.5) 
      {
        fillthis= "B0s_per_bin_in_Lim_Acc";
        for (Int_t j = 0; j < fnPtBins; ++j)
        {
          if(fPtBinLimits[j] < ptMC[0] && ptMC[0] < fPtBinLimits[j+1]) {((TH1F*)(listout->FindObject(fillthis)))->Fill(j); break;}
        }
      }

      // We check if the tracks are in acceptance
      if(ptMC[1] < 0.1 || TMath::Abs(pseudoYMC[1]) > 0.9 ) continue;
      if(ptMC[3] < 0.1 || TMath::Abs(pseudoYMC[3]) > 0.9 ) continue;
      if(ptMC[5] < 0.1 || TMath::Abs(pseudoYMC[5]) > 0.9 ) continue;
      if(ptMC[6] < 0.1 || TMath::Abs(pseudoYMC[6]) > 0.9 ) continue;

      // We check if the B0 is in the fiducial region
      if(TMath::Abs(yMC[0]) > 0.8) continue;

      Int_t rows = B0toDStarPiLabelMatrix->GetNrows();

      B0toDStarPiLabelMatrix->ResizeTo(rows+1,7);
      particleMatrix(rows,0) = mcLabelPionB0;
      particleMatrix(rows,1) = mcLabelPionDStar;
      particleMatrix(rows,2) = mcLabelPionD0;
      particleMatrix(rows,3) = mcLabelKaon;
      particleMatrix(rows,4) = mcLabelD0;
      particleMatrix(rows,5) = mcLabelDStar;
      particleMatrix(rows,6) = mcLabelB0;

      fillthis= "mc_B0_pt_bins_acc";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(ptMC[0]);

      if(TMath::Abs(yMC[0]) < 0.5) 
      {
        fillthis= "mc_B0_pt_bins_acc_lim_acc";
        ((TH1F*)(listout->FindObject(fillthis)))->Fill(ptMC[0]);
      }

      fillthis= "B0s_in_analysis";
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(2);


      fillthis= "B0s_per_bin_in_Acc";
      for (Int_t j = 0; j < fnPtBins; ++j)
      {
        if(fPtBinLimits[j] < ptMC[0] && ptMC[0] < fPtBinLimits[j+1]) {((TH1F*)(listout->FindObject(fillthis)))->Fill(j); break;}
      }

    }
  }

  //old method

  // Not all the tracks can be/are detected by the detector. We are only interested in tracks that lie within the acceptance of the detector.
  // We remove the undetectable tracks from the array in order to get accurate information on the amount of signal that lies within the acceptance of our detector.
  // Int_t numberOfB0s = 0;
  // TArrayI correctLabelArray; 
  // for (Int_t i = 0; i < B0toDStarPiLabelMatrix->GetNrows(); i++)
  // {
  //   std::cout << "loop at row = " << i << std::endl;
  //   Int_t particleCounter = 0;
  //   for (Int_t j = 0; j < 4; j++)
  //   {
  //     Int_t labelParticleInList = (Int_t)particleMatrix(i,j);
  //     for (Int_t k=0; k<aodevent->GetNumberOfTracks(); k++)
  //     { 
  //       AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(aodevent->GetTrack(k));
  //       if(!aodTrack) AliFatal("Not a standard AOD");
  //       if(TMath::Abs(aodTrack->Eta())>0.8) continue;
  //       if(aodTrack->GetLabel() == labelParticleInList) 
  //       {
  //         particleCounter++;
  //         break;
  //       }
  //     }
  //   }
  //   if(particleCounter==4) std::cout << "found 4" << std::endl;
  //   if (particleCounter==4)
  //   {
  //     TString fillthis= "B0s_in_analysis";
  //     ((TH1F*)(listout->FindObject(fillthis)))->Fill(2);
  //     numberOfB0s++;
  //     correctLabelArray.Set(numberOfB0s);
  //     correctLabelArray.AddAt(i,numberOfB0s-1);

  //     Int_t labelParticle = (Int_t)particleMatrix(i,0);
  //     AliAODMCParticle * B0track = dynamic_cast< AliAODMCParticle*>(mcTrackArray->At(labelParticle));

  //     fillthis= "B0s_per_bin_in_Acc";
  //     for (Int_t j = 0; j < fnPtBins; ++j)
  //     {
  //       if(fPtBinLimits[j] < B0track->Pt() && B0track->Pt() < fPtBinLimits[j+1]) {((TH1F*)(listout->FindObject(fillthis)))->Fill(j); break;} 
  //     }
  //   }
  // }

  // std::cout << "Number of B0 = " << numberOfB0s << std::endl;

  // for (Int_t i = 0; i < correctLabelArray.GetSize(); i++)
  // {
  //   particleMatrix(i,0) = (Int_t)particleMatrix(correctLabelArray[i],0);
  //   particleMatrix(i,1) = (Int_t)particleMatrix(correctLabelArray[i],1);
  //   particleMatrix(i,2) = (Int_t)particleMatrix(correctLabelArray[i],2);
  //   particleMatrix(i,3) = (Int_t)particleMatrix(correctLabelArray[i],3);
  //   particleMatrix(i,4) = (Int_t)particleMatrix(correctLabelArray[i],4);
  //   particleMatrix(i,5) = (Int_t)particleMatrix(correctLabelArray[i],5);
  //   particleMatrix(i,6) = (Int_t)particleMatrix(correctLabelArray[i],6);
  // }
  // B0toDStarPiLabelMatrix->ResizeTo(correctLabelArray.GetSize(),7);
  return;
}
//-------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSEB0toDStarPi::D0FirstDaughterSelection(AliAODTrack* aodTrack, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, AliAODMCHeader * header){
  
  // we select the D0 pion and save its information
  if(!aodTrack) AliFatal("Not a standard AOD");

  //quick quality cut
  if(aodTrack->GetITSNcls() < 1) return kFALSE;
  if(aodTrack->GetTPCNcls() < 1) return kFALSE;
  if(aodTrack->GetStatus()&AliESDtrack::kITSpureSA) return kFALSE;
  if(!(aodTrack->GetStatus()&AliESDtrack::kITSin)) return kFALSE;
  if(aodTrack->GetID() < 0) return kFALSE;
  Double_t covtest[21];
  if(!aodTrack->GetCovarianceXYZPxPyPz(covtest)) return kFALSE;

  Int_t mcLabelParticle = -1;
  Int_t pdgParticle = -1;
  mcLabelParticle = aodTrack->GetLabel();

  // we fill histograms with information of the track
  Double_t pt_track = aodTrack->Pt();
  Double_t momentum_track = aodTrack->P();
  Int_t numberOfITS = aodTrack->GetITSNcls(); 
  Int_t numberOfTPC = aodTrack->GetTPCNcls();

  AliExternalTrackParam particleTrack;
  particleTrack.CopyFromVTrack(aodTrack);
  Double_t d0[2],covd0[3];
  particleTrack.PropagateToDCA(primaryVertex,bz,100.,d0,covd0);

  //we check if the particle is a signal track, we look at both daughter options therefore the signal will be too high in the D0 daughter signal histograms
  Bool_t isDesiredCandidate = kFALSE;
  if(fUseMCInfo){
    TMatrix &particleMatrix = *B0toDStarPiLabelMatrix;
    for (Int_t k = 0; k < B0toDStarPiLabelMatrix->GetNrows(); ++k){
      if(mcLabelParticle == (Int_t)particleMatrix(k,2) || mcLabelParticle == (Int_t)particleMatrix(k,3)){
        isDesiredCandidate = kTRUE;
        break;
      }
    }
  }

  if(fUseMCInfo){
    if(IsTrackInjected(aodTrack,header,mcTrackArray) && !isDesiredCandidate && fQuickSignalAnalysis == 2) return kFALSE;
  }

  Int_t daughterType = 0;


  Int_t histType = 0;
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

  for (Int_t j = 0; j < 10; ++j)
  {
    if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
      
  }
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

  if(isDesiredCandidate)
  {
    histType = 1;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
  }

  //we apply a number of cuts on the particle
   Bool_t bCut = kFALSE;

  //We do not apply a PID cut at this stage since we don't know if we are dealing with a kaon or a pion

  if(aodTrack->GetITSNcls() < fCuts->GetMinITSNclsD0FirstDaughter()){
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(3);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(3);
    bCut = kTRUE;
  }

  if(aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsD0FirstDaughter()){ 
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(4);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(4);
    bCut = kTRUE;
  }

  if(fCuts->UseITSRefitD0FirstDaughter()==kTRUE){
    if(!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(5);
      } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(5);
      bCut = kTRUE;
    }
  }

  if(fCuts->UseTPCRefitD0FirstDaughter()==kTRUE){
    if((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(6);
      } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(6);
      bCut = kTRUE;
    }
  }

  if(fCuts->UseFilterBitD0FirstDaughter()==kTRUE){
    if(!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitD0FirstDaughter())))) {
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(7);
      } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(7);
      bCut = kTRUE;
    }
  }

  if(aodTrack->Pt() < fCuts->GetMinPtD0FirstDaughter()){ 
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(8);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(8);
    bCut = kTRUE;
  }

  if(TMath::Abs(d0[0]) < fCuts->GetMind0D0FirstDaughter()){ 
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(12);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(12);
    bCut = kTRUE;
  }

  if(TMath::Abs(aodTrack->Eta()) > fCuts->GetMaxAbsEtaD0FirstDaughter()){ 
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(9);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(9);
    bCut = kTRUE;
  }

  Bool_t bHardSelectionArrayITS[7] = {kFALSE};
  fCuts->GetHardSelectionArrayITSD0FirstDaughter(bHardSelectionArrayITS);
  Bool_t bSoftSelectionArrayITS[7] = {kFALSE};
  fCuts->GetSoftSelectionArrayITSD0FirstDaughter(bSoftSelectionArrayITS);

  Bool_t bHardITSPass = kTRUE;
  for (Int_t j = 0; j < 7; ++j)
  {
    if(bHardSelectionArrayITS[j]) 
    {
      if(!aodTrack->HasPointOnITSLayer(j)) bHardITSPass = kFALSE;
    }
  }

  Int_t nCounterSoftSelection = 0;
  Bool_t bSoftITSPass = kTRUE;
  for (Int_t j = 0; j < 7; ++j)
  {
    if(bSoftSelectionArrayITS[j]) 
    {
      if(aodTrack->HasPointOnITSLayer(j)) nCounterSoftSelection++;
    }
  }    
  if(nCounterSoftSelection < fCuts->GetNSoftITSCutD0FirstDaughter()) bSoftITSPass = kFALSE;

  if(!bHardITSPass){ 
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(10);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(10);
    bCut = kTRUE;
  }

  if(!bSoftITSPass){ 
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(11);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(11);
    bCut = kTRUE;
  }

  if(!isDesiredCandidate && fQuickSignalAnalysis == 1) bCut = kTRUE;

  if(bCut) {
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[0][1])->Fill(0);
    } else ((TH1F*)fDaughterHistogramArrayExtra[0][0])->Fill(0);
    return kFALSE;
  }

  //we fill histograms with track information of the tracks that pass the cuts
  histType = 2;
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

  for (Int_t j = 0; j < 10; ++j)
  {
    if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
      
  }
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

  if(isDesiredCandidate)
  {
    histType = 3;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
  }

  return kTRUE;
}
//-------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSEB0toDStarPi::D0SecondDaughterSelection(AliAODTrack* aodTrack, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, AliAODMCHeader * header){
  
  // we select the D0 pion and save its information
  if(!aodTrack) AliFatal("Not a standard AOD");

  //quick quality cut
  if(aodTrack->GetITSNcls() < 1) return kFALSE;
  if(aodTrack->GetTPCNcls() < 1) return kFALSE;
  if(aodTrack->GetStatus()&AliESDtrack::kITSpureSA) return kFALSE;
  if(!(aodTrack->GetStatus()&AliESDtrack::kITSin)) return kFALSE;
  if(aodTrack->GetID() < 0) return kFALSE;
  Double_t covtest[21];
  if(!aodTrack->GetCovarianceXYZPxPyPz(covtest)) return kFALSE;

  Int_t mcLabelParticle = -1;
  Int_t pdgParticle = -1;
  mcLabelParticle = aodTrack->GetLabel();

  // we fill histograms with information of the track
  Double_t pt_track = aodTrack->Pt();
  Double_t momentum_track = aodTrack->P();
  Int_t numberOfITS = aodTrack->GetITSNcls(); 
  Int_t numberOfTPC = aodTrack->GetTPCNcls();

  AliExternalTrackParam particleTrack;
  particleTrack.CopyFromVTrack(aodTrack);
  Double_t d0[2],covd0[3];
  particleTrack.PropagateToDCA(primaryVertex,bz,100.,d0,covd0);

  //we check if the particle is a signal track, we look at both daughter options therefore the signal will be too high in the D0 daughter signal histograms
  Bool_t isDesiredCandidate = kFALSE;
  if(fUseMCInfo){
    TMatrix &particleMatrix = *B0toDStarPiLabelMatrix;
    for (Int_t k = 0; k < B0toDStarPiLabelMatrix->GetNrows(); ++k){
      if(mcLabelParticle == (Int_t)particleMatrix(k,2) || mcLabelParticle == (Int_t)particleMatrix(k,3)){
        isDesiredCandidate = kTRUE;
        break;
      }
    }
  }

  if(fUseMCInfo){
    if(IsTrackInjected(aodTrack,header,mcTrackArray) && !isDesiredCandidate && fQuickSignalAnalysis == 2) return kFALSE;
  }

  Int_t daughterType = 1;


  Int_t histType = 0;
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

  for (Int_t j = 0; j < 10; ++j)
  {
    if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
      
  }
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

  if(isDesiredCandidate)
  {
    histType = 1;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
  }

  //we apply a number of cuts on the particle
   Bool_t bCut = kFALSE;

  //We do not apply a PID cut at this stage since we don't know if we are dealing with a kaon or a pion

  if(aodTrack->GetITSNcls() < fCuts->GetMinITSNclsD0SecondDaughter()){
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(3);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(3);
    bCut = kTRUE;
  }

  if(aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsD0SecondDaughter()){ 
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(4);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(4);
    bCut = kTRUE;
  }

  if(fCuts->UseITSRefitD0SecondDaughter()==kTRUE){
    if(!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(5);
      } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(5);
      bCut = kTRUE;
    }
  }

  if(fCuts->UseTPCRefitD0SecondDaughter()==kTRUE){
    if((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(6);
      } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(6);
      bCut = kTRUE;
    }
  }

  if(fCuts->UseFilterBitD0SecondDaughter()==kTRUE){
    if(!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitD0SecondDaughter())))) {
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(7);
      } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(7);
      bCut = kTRUE;
    }
  }

  if(aodTrack->Pt() < fCuts->GetMinPtD0SecondDaughter()){ 
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(8);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(8);
    bCut = kTRUE;
  }

  if(TMath::Abs(d0[0]) < fCuts->GetMind0D0SecondDaughter()){ 
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(12);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(12);
    bCut = kTRUE;
  }

  if(TMath::Abs(aodTrack->Eta()) > fCuts->GetMaxAbsEtaD0SecondDaughter()){ 
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(9);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(9);
    bCut = kTRUE;
  }

  Bool_t bHardSelectionArrayITS[7] = {kFALSE};
  fCuts->GetHardSelectionArrayITSD0SecondDaughter(bHardSelectionArrayITS);
  Bool_t bSoftSelectionArrayITS[7] = {kFALSE};
  fCuts->GetSoftSelectionArrayITSD0SecondDaughter(bSoftSelectionArrayITS);

  Bool_t bHardITSPass = kTRUE;
  for (Int_t j = 0; j < 7; ++j)
  {
    if(bHardSelectionArrayITS[j]) 
    {
      if(!aodTrack->HasPointOnITSLayer(j)) bHardITSPass = kFALSE;
    }
  }

  Int_t nCounterSoftSelection = 0;
  Bool_t bSoftITSPass = kTRUE;
  for (Int_t j = 0; j < 7; ++j)
  {
    if(bSoftSelectionArrayITS[j]) 
    {
      if(aodTrack->HasPointOnITSLayer(j)) nCounterSoftSelection++;
    }
  }    
  if(nCounterSoftSelection < fCuts->GetNSoftITSCutD0SecondDaughter()) bSoftITSPass = kFALSE;

  if(!bHardITSPass){ 
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(10);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(10);
    bCut = kTRUE;
  }

  if(!bSoftITSPass){ 
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(11);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(11);
    bCut = kTRUE;
  }

  if(!isDesiredCandidate && fQuickSignalAnalysis == 1) bCut = kTRUE;

  if(bCut) {
    if(isDesiredCandidate) {
      ((TH1F*)fDaughterHistogramArrayExtra[1][1])->Fill(0);
    } else ((TH1F*)fDaughterHistogramArrayExtra[1][0])->Fill(0);
    return kFALSE;
  }

  //we fill histograms with track information of the tracks that pass the cuts
  histType = 2;
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

  for (Int_t j = 0; j < 10; ++j)
  {
    if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
      
  }
  ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

  if(isDesiredCandidate)
  {
    histType = 3;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
  }

  return kTRUE;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::DStarPionSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, AliAODMCHeader * header){
  
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
    if(aodTrack->GetStatus()&AliESDtrack::kITSpureSA) continue;
    if(!(aodTrack->GetStatus()&AliESDtrack::kITSin)) continue;
    if(aodTrack->GetID() < 0) continue;
    Double_t covtest[21];
    if(!aodTrack->GetCovarianceXYZPxPyPz(covtest)) continue;

    Int_t mcLabelParticle = -1;
    Int_t pdgParticle = -1;
    mcLabelParticle = aodTrack->GetLabel();

    numberofparticles++;

    // we fill histograms with information of the track
    Double_t pt_track = aodTrack->Pt();
    Double_t momentum_track = aodTrack->P();
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

    AliExternalTrackParam particleTrack;
    particleTrack.CopyFromVTrack(aodTrack);
    Double_t d0[2],covd0[3];
    particleTrack.PropagateToDCA(primaryVertex,bz,100.,d0,covd0);

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

    if(fUseMCInfo){
      if(IsTrackInjected(aodTrack,header,mcTrackArray) && !isDesiredCandidate && fQuickSignalAnalysis == 2) continue;
    }

    Int_t daughterType = 2;


    Int_t histType = 0;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }

    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

    if(isDesiredCandidate)
    {
      histType = 1;
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

      for (Int_t j = 0; j < 10; ++j)
      {
        if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
          
      }

      if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
      if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
      if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
    }

    //we apply a number of cuts on the particle
     Bool_t bCut = kFALSE;

    //we apply a PID cut for a pion
    if(!(fCuts->SelectPID(aodTrack,pionPIDnumber))){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(2);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(2);
      bCut = kTRUE;
    }

    if(aodTrack->GetITSNcls() < fCuts->GetMinITSNclsDStarPion()){
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(3);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(3);
      bCut = kTRUE;
    }

    if(aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsDStarPion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(4);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(4);
      bCut = kTRUE;
    }

    if(fCuts->UseITSRefitDStarPion()==kTRUE){
      if(!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(5);
        } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(5);
        bCut = kTRUE;
      }
    }

    if(fCuts->UseTPCRefitDStarPion()==kTRUE){
      if((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(6);
        } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(6);
        bCut = kTRUE;
      }
    }

    if(fCuts->UseFilterBitDStarPion()==kTRUE){
      if(!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitDStarPion())))) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(7);
        } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(7);
        bCut = kTRUE;
      }
    }

    if(aodTrack->Pt() < fCuts->GetMinPtDStarPion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(8);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(8);
      bCut = kTRUE;
    }

    if(aodTrack->Pt() > fCuts->GetMaxPtDStarPion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(13);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(13);
      bCut = kTRUE;
    }


    if(TMath::Abs(d0[0]) < fCuts->GetMind0DStarPion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(12);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(12);
      bCut = kTRUE;
    }


    if(TMath::Abs(aodTrack->Eta()) > fCuts->GetMaxAbsEtaDStarPion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(9);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(9);
      bCut = kTRUE;
    }

    Bool_t bHardSelectionArrayITS[7] = {kFALSE};
    fCuts->GetHardSelectionArrayITSDStarPion(bHardSelectionArrayITS);
    Bool_t bSoftSelectionArrayITS[7] = {kFALSE};
    fCuts->GetSoftSelectionArrayITSDStarPion(bSoftSelectionArrayITS);

    Bool_t bHardITSPass = kTRUE;
    for (Int_t j = 0; j < 7; ++j)
    {
      if(bHardSelectionArrayITS[j]) 
      {
        if(!aodTrack->HasPointOnITSLayer(j)) bHardITSPass = kFALSE;
      }
    }

    Int_t nCounterSoftSelection = 0;
    Bool_t bSoftITSPass = kTRUE;
    for (Int_t j = 0; j < 7; ++j)
    {
      if(bSoftSelectionArrayITS[j]) 
      {
        if(aodTrack->HasPointOnITSLayer(j)) nCounterSoftSelection++;
      }
    }    
    if(nCounterSoftSelection < fCuts->GetNSoftITSCutDStarPion()) bSoftITSPass = kFALSE;

    if(!bHardITSPass){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(10);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(10);
      bCut = kTRUE;
    }

    if(!bSoftITSPass){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(11);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(11);
      bCut = kTRUE;
    }

    if(!isDesiredCandidate && fQuickSignalAnalysis == 1) bCut = kTRUE;

    if(bCut) {
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[2][1])->Fill(0);
      } else ((TH1F*)fDaughterHistogramArrayExtra[2][0])->Fill(0);
      continue;
    }

    //we fill histograms with track information of the tracks that pass the cuts
    histType = 2;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }

    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

    if(isDesiredCandidate)
    {
      histType = 3;
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

      for (Int_t j = 0; j < 10; ++j)
      {
        if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
          
      }

      if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
      if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
      if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
    }

    fDStarPionTracks->push_back(i);
    numberofparticlesused++;
  }
  
  ((TH1F*)fDaughterHistogramArray[2][0][12])->Fill(numberofparticles);
  ((TH1F*)fDaughterHistogramArray[2][1][12])->Fill(numberofparticlesused);
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::B0PionSelection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz, TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, AliAODMCHeader * header){
  
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
    if(aodTrack->GetStatus()&AliESDtrack::kITSpureSA) continue;
    if(!(aodTrack->GetStatus()&AliESDtrack::kITSin)) continue;
    if(aodTrack->GetID() < 0) continue;
    Double_t covtest[21];
    if(!aodTrack->GetCovarianceXYZPxPyPz(covtest)) continue;


    Int_t mcLabelParticle = -1;
    Int_t pdgParticle = -1;
    mcLabelParticle = aodTrack->GetLabel();

    numberofparticles++;

    // we fill histograms with information of the track
    Double_t pt_track = aodTrack->Pt();
    Double_t momentum_track = aodTrack->P();
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

    AliExternalTrackParam particleTrack;
    particleTrack.CopyFromVTrack(aodTrack);
    Double_t d0[2],covd0[3];
    particleTrack.PropagateToDCA(primaryVertex,bz,100.,d0,covd0);


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

    if(fUseMCInfo){
      if(IsTrackInjected(aodTrack,header,mcTrackArray) && !isDesiredCandidate && fQuickSignalAnalysis == 2) continue;
    }


    Int_t daughterType = 3;

    Int_t histType = 0;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }

    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

    if(isDesiredCandidate)
    {
      histType = 1;
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

      for (Int_t j = 0; j < 10; ++j)
      {
        if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
          
      }

      if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
      if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
      if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
    }

    //we apply a number of cuts on the particle
    Bool_t bCut = kFALSE;

    //we apply a PID cut, 2 is used to indicate we look for a pion
    if(!(fCuts->SelectPID(aodTrack,2))){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[3][1])->Fill(2);
      } else ((TH1F*)fDaughterHistogramArrayExtra[3][0])->Fill(2);
      bCut = kTRUE;
    }

    if(aodTrack->GetITSNcls() < fCuts->GetMinITSNclsB0Pion()){
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[3][1])->Fill(3);
      } else ((TH1F*)fDaughterHistogramArrayExtra[3][0])->Fill(3);
      bCut = kTRUE;
    }

    if(aodTrack->GetTPCNcls() < fCuts->GetMinTPCNclsB0Pion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[3][1])->Fill(4);
      } else ((TH1F*)fDaughterHistogramArrayExtra[3][0])->Fill(4);
      bCut = kTRUE;
    }

    if(fCuts->UseITSRefitB0Pion()==kTRUE){
      if(!(aodTrack->GetStatus()&AliESDtrack::kITSrefit)) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArrayExtra[3][1])->Fill(5);
        } else ((TH1F*)fDaughterHistogramArrayExtra[3][0])->Fill(5);
        bCut = kTRUE;
      }
    }

    if(fCuts->UseTPCRefitB0Pion()==kTRUE){
      if((!(aodTrack->GetStatus()&AliESDtrack::kTPCrefit))) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArrayExtra[3][1])->Fill(6);
        } else ((TH1F*)fDaughterHistogramArrayExtra[3][0])->Fill(6);
        bCut = kTRUE;
      }
    }

    if(fCuts->UseFilterBitB0Pion()==kTRUE){
      if(!(aodTrack->TestFilterMask(BIT(fCuts->GetFilterBitB0Pion())))) {
        if(isDesiredCandidate) {
          ((TH1F*)fDaughterHistogramArrayExtra[3][1])->Fill(7);
        } else ((TH1F*)fDaughterHistogramArrayExtra[3][0])->Fill(7);
        bCut = kTRUE;
      }
    }


    if(aodTrack->Pt() < fCuts->GetMinPtB0Pion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[3][1])->Fill(8);
      } else ((TH1F*)fDaughterHistogramArrayExtra[3][0])->Fill(8);
      bCut = kTRUE;
    }


    if(TMath::Abs(d0[0]) < fCuts->GetMind0B0Pion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[3][1])->Fill(12);
      } else ((TH1F*)fDaughterHistogramArrayExtra[3][0])->Fill(12);
      bCut = kTRUE;
    }

    if(TMath::Abs(aodTrack->Eta()) > fCuts->GetMaxAbsEtaB0Pion()){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[3][1])->Fill(9);
      } else ((TH1F*)fDaughterHistogramArrayExtra[3][0])->Fill(9);
      bCut = kTRUE;
    }

    Bool_t bHardSelectionArrayITS[7] = {kFALSE};
    fCuts->GetHardSelectionArrayITSB0Pion(bHardSelectionArrayITS);
    Bool_t bSoftSelectionArrayITS[7] = {kFALSE};
    fCuts->GetSoftSelectionArrayITSB0Pion(bSoftSelectionArrayITS);

    Bool_t bHardITSPass = kTRUE;
    for (Int_t j = 0; j < 7; ++j)
    {
      if(bHardSelectionArrayITS[j]) 
      {
        if(!aodTrack->HasPointOnITSLayer(j)) bHardITSPass = kFALSE;
      }
    }

    Int_t nCounterSoftSelection = 0;
    Bool_t bSoftITSPass = kTRUE;
    for (Int_t j = 0; j < 7; ++j)
    {
      if(bSoftSelectionArrayITS[j]) 
      {
        if(aodTrack->HasPointOnITSLayer(j)) nCounterSoftSelection++;
      }
    }    
    if(nCounterSoftSelection < fCuts->GetNSoftITSCutB0Pion()) bSoftITSPass = kFALSE;

    if(!bHardITSPass){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[3][1])->Fill(10);
      } else ((TH1F*)fDaughterHistogramArrayExtra[3][0])->Fill(10);
      bCut = kTRUE;
    }

    if(!bSoftITSPass){ 
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[3][1])->Fill(11);
      } else ((TH1F*)fDaughterHistogramArrayExtra[3][0])->Fill(11);
      bCut = kTRUE;
    }


    if(!isDesiredCandidate && fQuickSignalAnalysis == 1) bCut = kTRUE;

    if(bCut) {
      if(isDesiredCandidate) {
        ((TH1F*)fDaughterHistogramArrayExtra[3][1])->Fill(0);
      } else ((TH1F*)fDaughterHistogramArrayExtra[3][0])->Fill(0);
      continue;
    }

    //we fill histograms with track information of the tracks that pass the cuts
    histType = 2;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }

    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);

    if(isDesiredCandidate)
    {
      histType = 3;
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

      for (Int_t j = 0; j < 10; ++j)
      {
        if(aodTrack->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
          
      }

      if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
      if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
      if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
      ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0[0]);
    }


    fB0PionTracks->push_back(i);
    numberofparticlesused++;
  }

  ((TH1F*)fDaughterHistogramArray[3][0][12])->Fill(numberofparticles);
  ((TH1F*)fDaughterHistogramArray[3][1][12])->Fill(numberofparticlesused);
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::D0Selection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz,TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, TClonesArray * D0TracksFromFriendFile, AliAODMCHeader * header){

  TString fillthis = "";

  AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();

  //next we loop over all the D0 candidates
  for (Int_t j = 0; j < D0TracksFromFriendFile->GetEntriesFast(); j++)
  {

    //we get the track of the D0
    AliAODRecoDecayHF2Prong * trackD0 = (AliAODRecoDecayHF2Prong*)(D0TracksFromFriendFile->At(j));
    if(!trackD0) {std::cout << "found none" << std::endl; continue;}
    if(trackD0 == nullptr) {std::cout << "found nullptr" << std::endl; continue;}

    if(!(vHF->FillRecoCand(aodEvent,trackD0))) //Fill the data members of the candidate only if they are empty.
    {
      fCEvents->Fill(12); //monitor how often this fails 
      continue;
    }

    AliAODTrack * trackFirstDaughter = (AliAODTrack*)(trackD0->GetDaughter(0));
    AliAODTrack * trackSecondDaughter = (AliAODTrack*)(trackD0->GetDaughter(1));
    if(!D0FirstDaughterSelection(trackFirstDaughter, primaryVertex, bz, mcTrackArray, B0toDStarPiLabelMatrix,header)) continue; 
    if(!D0SecondDaughterSelection(trackSecondDaughter, primaryVertex, bz, mcTrackArray, B0toDStarPiLabelMatrix,header)) continue; 


    AliAODVertex *vertexMother = (AliAODVertex*)trackD0->GetSecondaryVtx();

    //we save the pdgcode of the used particle and its mother and check if it is a desired candidate
    Int_t pdgCodeMother = -1;
    Float_t pdgCodeGrandMother = -1;
    Bool_t isDesiredCandidate = kFALSE;
    Int_t motherType, histType;
    motherType = 0;
    Int_t mcLabelD0 = -1;

    if(fUseMCInfo)
    {
      mcLabelD0 = MatchCandidateToMonteCarlo(421,trackD0,mcTrackArray,B0toDStarPiLabelMatrix);

      if(mcLabelD0 >= 0)
      {
        isDesiredCandidate = kTRUE;

        Int_t mcLabelFirstTrack = -1;
        mcLabelFirstTrack = trackFirstDaughter->GetLabel();

        if(mcLabelFirstTrack >= 0)
        {
          AliAODMCParticle *mcParticleFirstTrack = (AliAODMCParticle*)mcTrackArray->At(mcLabelFirstTrack);
          AliAODMCParticle *mcMotherParticle = (AliAODMCParticle*)mcTrackArray->At(mcLabelD0);

          if(mcParticleFirstTrack && mcMotherParticle)
          {
            pdgCodeMother = mcMotherParticle->GetPdgCode();

            Double_t vertex_distance = TMath::Sqrt((vertexMother->GetX() - mcParticleFirstTrack->Xv())*(vertexMother->GetX() - mcParticleFirstTrack->Xv()) + (vertexMother->GetY() - mcParticleFirstTrack->Yv())*(vertexMother->GetY() - mcParticleFirstTrack->Yv()) + (vertexMother->GetZ() - mcParticleFirstTrack->Zv())*(vertexMother->GetZ() - mcParticleFirstTrack->Zv()));
            ((TH1F*)fMotherHistogramArrayExtra[motherType][4])->Fill(vertex_distance);

            Double_t momentum_resolution = TMath::Sqrt((trackD0->Px() - mcMotherParticle->Px())*(trackD0->Px() - mcMotherParticle->Px()) + (trackD0->Py() - mcMotherParticle->Py())*(trackD0->Py() - mcMotherParticle->Py()) + (trackD0->Pz() - mcMotherParticle->Pz())*(trackD0->Pz() - mcMotherParticle->Pz()));
            ((TH1F*)fMotherHistogramArrayExtra[motherType][6])->Fill(momentum_resolution);      
          }
        }
      }
    }

    // We fill the histograms
    histType = 0;
    FillD0Histograms(trackD0, primaryVertex, bz, motherType, histType);
    if(isDesiredCandidate && fUseMCInfo){
      histType = 1;
      FillD0Histograms(trackD0, primaryVertex, bz, motherType, histType,pdgCodeMother);
    }

      // Here we apply cuts on the particle 
      Bool_t cutMother = kFALSE;

      Bool_t bCutArray[29] = {0};
      Int_t cutReturnValue = fCuts->IsD0forD0ptbinSelected(trackD0, 0, aodEvent, bCutArray);
      if(cutReturnValue == -1) cutMother = kTRUE;
      if(cutReturnValue == 0) cutMother = kTRUE;


      if(fGetCutInfo == kTRUE)
      {
        for (Int_t k = 0; k < 29; ++k)
        {
          if (bCutArray[k] == kTRUE){
            if(isDesiredCandidate){
              ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(k+1);
            } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(k+1);
            cutMother = kTRUE;
          }
        } 
      }


      if(!isDesiredCandidate && fQuickSignalAnalysis == 1) cutMother = kTRUE;

      if(cutMother){
        if(isDesiredCandidate){
          ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(0);
        } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(0);
        // delete vertexMother; vertexMother = nullptr; 
        // delete trackD0; trackD0 = nullptr;
        continue;
      }

    // We fill the cut histograms
    histType = 2;
    FillD0Histograms(trackD0, primaryVertex, bz, motherType, histType);
    if(isDesiredCandidate && fUseMCInfo){
      histType = 3;
      FillD0Histograms(trackD0, primaryVertex, bz, motherType, histType,pdgCodeMother);
    }

    //we save the location of the D0 candidate
    fD0Tracks->push_back(j);
  }

  delete vHF; vHF = nullptr;
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::DStarAndB0Selection(AliAODEvent* aodEvent, AliAODVertex *primaryVertex, Double_t bz,TClonesArray * mcTrackArray, TMatrix * B0toDStarPiLabelMatrix, TClonesArray * D0TracksFromFriendFile, AliAODMCHeader * header){
  
  TString fillthis = "";

  //we loop over all the DStar pion candidates
  for (Int_t i = 0; i < (Int_t)fDStarPionTracks->size(); i++)
  {
    //Save current Object count
    Int_t ObjectNumber = TProcessID::GetObjectCount();

    //we get the track of the DStar pion
    AliAODTrack * trackFirstDaughter = (AliAODTrack*)(aodEvent->GetTrack(fDStarPionTracks->at(i)));
    if(!trackFirstDaughter) continue;

    Int_t pdgD0 = 421;
    if(trackFirstDaughter->Charge() == -1) pdgD0 = -421;

    //next we loop over all the D0 candidates
    for (Int_t j = 0; j < (Int_t)fD0Tracks->size(); j++)
    {
      //we get the track of the D0
      AliAODRecoDecayHF2Prong * trackSecondDaughter = (AliAODRecoDecayHF2Prong*)(D0TracksFromFriendFile->At(fD0Tracks->at(j)));
      if(!trackSecondDaughter) {std::cout << "found none" << std::endl; continue;}
      if(trackSecondDaughter == nullptr) {std::cout << "found nullptr" << std::endl; continue;}

      //we check if the IDs of the tracks are different
      if(trackFirstDaughter->GetID() == trackSecondDaughter->GetProngID(0) || trackFirstDaughter->GetID() == trackSecondDaughter->GetProngID(1)) continue;

      //we check if the charges of the tracks are correct
      if(trackFirstDaughter->Charge() == trackSecondDaughter->Charge() || TMath::Abs(trackFirstDaughter->Charge() + trackSecondDaughter->Charge()) != 1) continue;
      
                  
      //we check if the pions have the same charge          
      if(trackFirstDaughter->Charge() == -1 && ((AliAODTrack*)trackSecondDaughter->GetDaughter(1))->Charge() != -1) continue;   
      if(trackFirstDaughter->Charge() == 1 && ((AliAODTrack*)trackSecondDaughter->GetDaughter(0))->Charge() != 1) continue;

      //we apply a PID cut on the D0 daughters

      if(trackFirstDaughter->Charge()==1)
      {
        if(!(fCuts->SelectPID(((AliAODTrack*)trackSecondDaughter->GetDaughter(0)),2))) continue;
        if(!(fCuts->SelectPID(((AliAODTrack*)trackSecondDaughter->GetDaughter(1)),3))) continue;
      } else if (trackFirstDaughter->Charge()==-1){
        if(!(fCuts->SelectPID(((AliAODTrack*)trackSecondDaughter->GetDaughter(0)),3))) continue;
        if(!(fCuts->SelectPID(((AliAODTrack*)trackSecondDaughter->GetDaughter(1)),2))) continue;        
      }

      //location DStar pion rotation around PV

      //we make an estimate of the DStar vertex and make an initial broad invariant mass window cut
      AliExternalTrackParam DStarPionTrackParam;
      DStarPionTrackParam.CopyFromVTrack(trackFirstDaughter);
      AliExternalTrackParam D0TrackParam;
      D0TrackParam.CopyFromVTrack(trackSecondDaughter);

      // we calculate the vertex of the mother candidate
      TObjArray tracksTestVertex;

      tracksTestVertex.Add(&DStarPionTrackParam);
      tracksTestVertex.Add(&D0TrackParam);

      Double_t dispersionTest = 0;
      AliAODVertex *testVertex = RecalculateVertex(primaryVertex,&tracksTestVertex,bz,dispersionTest);
      if(!testVertex) {delete testVertex; testVertex = nullptr; continue;}

      Double_t d0z0Test[2],covd0z0Test[3];

      //DStar creation with the new vertex
      DStarPionTrackParam.PropagateToDCA(testVertex,bz,100.,d0z0Test,covd0z0Test);
      D0TrackParam.PropagateToDCA(testVertex,bz,100.,d0z0Test,covd0z0Test);
      delete testVertex; testVertex = nullptr;

      Double_t pdgMassPion = TDatabasePDG::Instance()->GetParticle(211)->Mass();
      Double_t pdgMassD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
      Double_t pdgMassDStar = TDatabasePDG::Instance()->GetParticle(413)->Mass();

      Double_t energyDStarPion = pdgMassPion*pdgMassPion + DStarPionTrackParam.Px()*DStarPionTrackParam.Px()+DStarPionTrackParam.Py()*DStarPionTrackParam.Py()+DStarPionTrackParam.Pz()*DStarPionTrackParam.Pz();
      Double_t energyD0 = pdgMassD0*pdgMassD0 + D0TrackParam.Px()*D0TrackParam.Px()+D0TrackParam.Py()*D0TrackParam.Py()+D0TrackParam.Pz()*D0TrackParam.Pz();
      Double_t energySum = TMath::Sqrt(energyDStarPion) + TMath::Sqrt(energyD0);

      Double_t pxDStarTest = DStarPionTrackParam.Px() + D0TrackParam.Px();
      Double_t pyDStarTest = DStarPionTrackParam.Py() + D0TrackParam.Py();
      Double_t pzDStarTest = DStarPionTrackParam.Pz() + D0TrackParam.Pz();
      Double_t p2DStarTest = pxDStarTest*pxDStarTest + pyDStarTest*pyDStarTest + pzDStarTest*pzDStarTest;

      Double_t invMassDStarTest = TMath::Sqrt(energySum*energySum-p2DStarTest);

      //we use a mass window twice the size of the final cut. We cut here to speed up the code.
      Int_t nCutIndex = 0;
      Bool_t bCutArrayTemp[29];
      Double_t cutVariableValue = TMath::Abs(invMassDStarTest-pdgMassDStar)/2.0;
      Bool_t bPassedCut = fCuts->ApplyCutOnVariableDStarforDStarptbin(nCutIndex,0,cutVariableValue,bCutArrayTemp);
      if(!bPassedCut) continue;

      // Apply an impact product cut. We cut here to speed up the code.
      AliExternalTrackParam firstTrack;
      firstTrack.CopyFromVTrack(trackFirstDaughter);
      AliExternalTrackParam secondTrack;
      secondTrack.CopyFromVTrack(trackSecondDaughter);

      Double_t d0z0DStar[2],covd0z0DStar[3],d0DStar[2],d0errDStar[2];

      firstTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0DStar,covd0z0DStar);
      d0DStar[0] = d0z0DStar[0];
      d0errDStar[0] = TMath::Sqrt(covd0z0DStar[0]);
      secondTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0DStar,covd0z0DStar);
      d0DStar[1] = d0z0DStar[0];
      d0errDStar[1] = TMath::Sqrt(covd0z0DStar[0]); 

      nCutIndex = 10;
      cutVariableValue = d0DStar[0] * d0DStar[1];
      bPassedCut = fCuts->ApplyCutOnVariableDStarforDStarptbin(nCutIndex,0,cutVariableValue,bCutArrayTemp);
      if(!bPassedCut) continue;

      //we loop over all the B0 pion candidates
      for (Int_t k = 0; k < (Int_t)fB0PionTracks->size(); k++)
      {
        //we get the track of the first daughter
        AliAODTrack * trackB0Pion = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(fB0PionTracks->at(k)));
        if(!trackB0Pion) continue;

        //we check if the IDs of the tracks are different
        AliAODTrack* twoProngdaughter0 = (AliAODTrack*)trackSecondDaughter->GetDaughter(0);
        AliAODTrack* twoProngdaughter1 = (AliAODTrack*)trackSecondDaughter->GetDaughter(1);
        UShort_t idProng0 = twoProngdaughter0->GetID();
        UShort_t idProng1 = twoProngdaughter1->GetID();
        
        if(trackB0Pion->GetID() == trackFirstDaughter->GetID() || trackB0Pion->GetID() == idProng0 || trackB0Pion->GetID() == idProng1) continue;

        //we check if the charges of the tracks are correct // later change this for like sign analysis.
        Bool_t bSameSign = kFALSE;
        if(trackB0Pion->Charge() == (trackSecondDaughter->Charge() + trackFirstDaughter->Charge()) && trackB0Pion->Charge() + (trackSecondDaughter->Charge() + trackFirstDaughter->Charge()) != 0)  bSameSign = kTRUE;

        //location B0 pion rotation around PV
        for (Int_t iRot = 0; iRot < fNumberOfRotations + 1; ++iRot)
        {
          //we create a copy of the track that we will rotate
          AliAODTrack * trackB0PionRotated = new AliAODTrack(*trackB0Pion);

          //for iRot == 0, we use the original unrotated track. For iRot > 0 we rotate the track and set the label to -1
          if (iRot != 0)
          {
            //should still check if track is already at PV
            Double_t dPhiRotated = trackB0PionRotated->Phi() + TMath::Pi() - (TMath::Pi() * fDegreePerRotation * fNumberOfRotations / (180.0 * 2.0)) + (TMath::Pi() * fDegreePerRotation * iRot / 180.0);
            trackB0PionRotated->SetPhi(dPhiRotated);
          }


          ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////<<
          //
          // DStar Reconstruction
          //
          ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////<<

          //we use the DStar pion, B0 pion, and D0 tracks to reconstruct the vertex for the B0 and DStar decay
          AliExternalTrackParam thirdTrack;
          thirdTrack.CopyFromVTrack(trackB0PionRotated);

          // we calculate the vertex
          TObjArray daughterTracksWithRecalculation;

          daughterTracksWithRecalculation.Add(&firstTrack);
          daughterTracksWithRecalculation.Add(&secondTrack);
          daughterTracksWithRecalculation.Add(&thirdTrack);

          Double_t dispersion = 0;
          AliAODVertex *vertexMother = RecalculateVertex(primaryVertex,&daughterTracksWithRecalculation,bz,dispersion);
          if(!vertexMother) {
            delete vertexMother; vertexMother = nullptr; 
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;}


          if(vertexMother->GetNDaughters()!=2) 
          {
            std::cout << "bad reconstruction - number of daughters for vertex is incorrect" << std::endl;
            delete vertexMother; vertexMother = nullptr; 
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;
          }

          Double_t xdummyDStar=0.,ydummyDStar=0.,eDStar[2];
          // Double_t d0z0DStar[2],covd0z0DStar[3],d0DStar[2],d0errDStar[2];

          //DStar creation with the new vertex
          firstTrack.PropagateToDCA(vertexMother,bz,100.,d0z0DStar,covd0z0DStar);
          secondTrack.PropagateToDCA(vertexMother,bz,100.,d0z0DStar,covd0z0DStar);

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

          // firstTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0DStar,covd0z0DStar);
          // d0DStar[0] = d0z0DStar[0];
          // d0errDStar[0] = TMath::Sqrt(covd0z0DStar[0]);
          // secondTrack.PropagateToDCA(primaryVertex,bz,100.,d0z0DStar,covd0z0DStar);
          // d0DStar[1] = d0z0DStar[0];
          // d0errDStar[1] = TMath::Sqrt(covd0z0DStar[0]); 

          //Apply cuts on DCA
          Double_t dcaDStarPionD0 = secondTrack.GetDCA(&firstTrack,bz,xdummyDStar,ydummyDStar);
          Double_t dcaDStarPionB0Pion = secondTrack.GetDCA(&thirdTrack,bz,xdummyDStar,ydummyDStar);     
          Double_t dcaB0PionD0 = thirdTrack.GetDCA(&firstTrack,bz,xdummyDStar,ydummyDStar);  

          if(dcaDStarPionD0 > fCuts->GetMaxDCADStarPionD0())
          {
            delete vertexMother; vertexMother = nullptr; 
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;
          }
          if(dcaDStarPionB0Pion > fCuts->GetMaxDCADStarPionB0Pion())
          {
            delete vertexMother; vertexMother = nullptr; 
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;
          }
          if(dcaB0PionD0 > fCuts->GetMaxDCAB0PionD0())
          {
            delete vertexMother; vertexMother = nullptr; 
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;
          }

          Double_t dcaCombined = TMath::Sqrt(TMath::Abs(dcaDStarPionD0) + TMath::Abs(dcaDStarPionB0Pion) + TMath::Abs(dcaB0PionD0));
          if(dcaCombined > fCuts->GetMaxDCACombined())
          {
            delete vertexMother; vertexMother = nullptr; 
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;
          }

          Short_t chargeDStar = trackFirstDaughter->Charge() + trackSecondDaughter->Charge();
          AliAODVertex * vertexDStar = new AliAODVertex(*vertexMother);
          if(!vertexDStar)
          {
            std::cout << "no dstar vertex" << std::endl;
            delete vertexMother; vertexMother = nullptr; 
            delete vertexDStar; vertexDStar = nullptr; 
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;

          }

          Int_t nProngsDStar = 2;
          AliAODRecoDecayHF2Prong trackDStar(vertexDStar,pxDStar,pyDStar,pzDStar,d0DStar,d0errDStar,distanceAtVertex); 
          if(!&trackDStar) 
          {
            delete vertexMother; vertexMother = nullptr; 
            delete vertexDStar; vertexDStar = nullptr; 
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;
          }
          
          trackDStar.SetCharge(chargeDStar);

          UShort_t idDStar[2];
          idDStar[0]= trackFirstDaughter->GetID(); 
          idDStar[1]= 0;

          UInt_t prongsDStar[2];
          prongsDStar[0] = 211;
          prongsDStar[1] = 421;


          if(vertexDStar->GetNDaughters()!=2) 
          {
            std::cout << "bad reconstruction 2 - number of daughters for vertex is incorrect" << std::endl;
            delete vertexMother; vertexMother = nullptr; 
            delete vertexDStar; vertexDStar = nullptr; 
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;
          }

          trackDStar.GetSecondaryVtx()->AddDaughter(trackFirstDaughter);
          trackDStar.GetSecondaryVtx()->AddDaughter(trackSecondDaughter);
          trackDStar.SetPrimaryVtxRef((AliAODVertex*)aodEvent->GetPrimaryVertex());
          trackDStar.SetProngIDs(2,idDStar);
          

          ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////<<
          //
          // BO Reconstruction
          //
          ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////<<

          //location B0 pion rotation around SV

          // Use the new DStar candidate and the new vertex to create the B0 candidate
          Double_t xdummy=0.,ydummy=0.,dca,e[2];
          Double_t d0z0[2],covd0z0[3],d0[2],d0err[2];

          AliExternalTrackParam fourthTrack;
          fourthTrack.CopyFromVTrack(&trackDStar);

          thirdTrack.PropagateToDCA(vertexMother,bz,100.,d0z0,covd0z0);
          fourthTrack.PropagateToDCA(vertexMother,bz,100.,d0z0,covd0z0);

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


          Short_t chargeMother = trackFirstDaughter->Charge() + trackDStar.Charge();
          Int_t nProngsB0 = 2;
          AliAODRecoDecayHF2Prong trackB0(vertexMother,px,py,pz,d0,d0err,dca); 
          if(!&trackB0)
          {
            delete vertexMother; vertexMother = nullptr; 
            delete vertexDStar; vertexDStar = nullptr; 
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;
          }

          trackB0.SetCharge(chargeMother);
          
          trackB0.GetSecondaryVtx()->AddDaughter(trackB0PionRotated);
          trackB0.GetSecondaryVtx()->AddDaughter(&trackDStar);
          trackB0.SetPrimaryVtxRef((AliAODVertex*)aodEvent->GetPrimaryVertex());
          trackB0.SetProngIDs(2,id);

          // Fiducial cut
          if(TMath::Abs(trackB0.Y(511)) > 0.8) {
            delete vertexMother; vertexMother = nullptr; 
            delete vertexDStar; vertexDStar = nullptr; 
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;
          }

          ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////<<
          //
          // Cuts
          //
          ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////<<

          // We check if the B0 candidate is a true signal in Monte Carlo
          Bool_t isDesiredCandidate = kFALSE;
          Int_t mcLabelB0 = -1;
          Int_t mcLabelDStar = -1;
          fillthis = "";
          Int_t motherType, histType;
          motherType = 1;

          if(fUseMCInfo)
          {
            mcLabelDStar = MatchCandidateToMonteCarlo(413,&trackDStar,mcTrackArray,B0toDStarPiLabelMatrix);
            mcLabelB0 = MatchCandidateToMonteCarlo(511,&trackB0,mcTrackArray,B0toDStarPiLabelMatrix);

            if (mcLabelB0 >= 0 && mcLabelDStar >= 0 && trackB0PionRotated->GetLabel() >= 0 && iRot == 0)
            {
              AliAODMCParticle *mcTrackDStarPion = (AliAODMCParticle*)mcTrackArray->At(trackB0PionRotated->GetLabel());
              AliAODMCParticle *mcTrackDStar = (AliAODMCParticle*)mcTrackArray->At(mcLabelDStar);
              // DStar
              Double_t vertex_distance = TMath::Sqrt((vertexMother->GetX() - mcTrackDStarPion->Xv())*(vertexMother->GetX() - mcTrackDStarPion->Xv()) + (vertexMother->GetY() - mcTrackDStarPion->Yv())*(vertexMother->GetY() - mcTrackDStarPion->Yv()) + (vertexMother->GetZ() - mcTrackDStarPion->Zv())*(vertexMother->GetZ() - mcTrackDStarPion->Zv()));
              ((TH1F*)fMotherHistogramArrayExtra[motherType][4])->Fill(vertex_distance);

              Double_t momentum_resolution = TMath::Sqrt((trackDStar.Px() - mcTrackDStar->Px())*(trackDStar.Px() - mcTrackDStar->Px()) + (trackDStar.Py() - mcTrackDStar->Py())*(trackDStar.Py() - mcTrackDStar->Py()) + (trackDStar.Pz() - mcTrackDStar->Pz())*(trackDStar.Pz() - mcTrackDStar->Pz()));
              ((TH1F*)fMotherHistogramArrayExtra[motherType][6])->Fill(momentum_resolution);

              isDesiredCandidate = kTRUE;
            }
          }
          
          // We check if the signal is injected, optionally we can reject injected signals
          Bool_t bIsInjected = kFALSE;
          if(fUseMCInfo)
          {
            if(fCheckInjected)
            {
              bIsInjected = IsCandidateInjected(&trackB0, header,mcTrackArray);
              if(fRemoveInjected == 1 && isDesiredCandidate == kFALSE && bIsInjected) {
                delete vertexMother; vertexMother = nullptr; 
                delete vertexDStar; vertexDStar = nullptr; 
                delete trackB0PionRotated; trackB0PionRotated = nullptr; 
                continue;
              }
              if(fRemoveInjected == 2 && bIsInjected) {
                delete vertexMother; vertexMother = nullptr; 
                delete vertexDStar; vertexDStar = nullptr; 
                delete trackB0PionRotated; trackB0PionRotated = nullptr; 
                continue;
              }
            } 
          }

          // We fill the DStar histograms
          histType = 0;
          FillDStarAndB0Histograms(&trackDStar, primaryVertex, bz, motherType, histType);
          if(isDesiredCandidate && fUseMCInfo)
          {
            histType = 1;
            FillDStarAndB0Histograms(&trackDStar, primaryVertex, bz, motherType, histType);
          }
       
          // We apply cuts on the DStar
          Bool_t cutDStar = kFALSE;

          Bool_t bCutArrayDStar[29] = {0};
          Int_t cutReturnValueDStar = fCuts->IsDStarforDStarptbinSelected(&trackDStar, 0, aodEvent, bCutArrayDStar);
          if(cutReturnValueDStar == -1) cutDStar = kTRUE;
          if(cutReturnValueDStar == 0) cutDStar = kTRUE;

          Bool_t bCutArrayD0[39] = {0};
          Int_t cutReturnValueD0 = fCuts->IsD0forDStarptbinSelected(&trackDStar, 0, aodEvent, bCutArrayD0);
          if(cutReturnValueD0 == -1) cutDStar = kTRUE;
          if(cutReturnValueD0 == 0) cutDStar = kTRUE;


          if(fGetCutInfo == kTRUE)
          {

            for (Int_t n = 0; n < 29; ++n)
            {
              if(bCutArrayDStar[n] == kTRUE){
                if(isDesiredCandidate){
                  ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(n+1);
                } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(n+1);
                cutDStar = kTRUE;
              }
            }

            for (Int_t n = 0; n < 39; ++n)
            {
              if(bCutArrayD0[n] == kTRUE){
                if(isDesiredCandidate){
                  ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(n+1+39);
                } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(n+1+39);
                cutDStar = kTRUE;
              }
            }
          }

          if(!isDesiredCandidate && fQuickSignalAnalysis == 1) cutDStar = kFALSE;

          if(cutDStar)
          {
            if(isDesiredCandidate)
            {
              ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(0);
            } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(0);
            delete vertexMother; vertexMother = nullptr;
            delete vertexDStar; vertexDStar = nullptr;
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;
          }

          // We fill the cut histograms
          histType = 2;
          FillDStarAndB0Histograms(&trackDStar, primaryVertex, bz, motherType, histType);
          if(isDesiredCandidate && fUseMCInfo)
          {
            histType = 3;
            FillDStarAndB0Histograms(&trackDStar, primaryVertex, bz, motherType, histType);
          }

          ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////<<
          //
          // BO Reconstruction
          //
          //
          ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////<<


          //we get information about the reconstructed B0
          Double_t ptMother = trackB0.Pt();

          fillthis = "";

          motherType = 2;
          histType = 0;

          if(isDesiredCandidate)
          {
            AliAODMCParticle *mcTrackFirstDaughter = (AliAODMCParticle*)mcTrackArray->At(trackB0PionRotated->GetLabel());
            AliAODMCParticle *mcTrackB0 = (AliAODMCParticle*)mcTrackArray->At(mcLabelB0);

            Double_t vertex_distance = TMath::Sqrt((vertexMother->GetX() - mcTrackFirstDaughter->Xv())*(vertexMother->GetX() - mcTrackFirstDaughter->Xv()) + (vertexMother->GetY() - mcTrackFirstDaughter->Yv())*(vertexMother->GetY() - mcTrackFirstDaughter->Yv()) + (vertexMother->GetZ() - mcTrackFirstDaughter->Zv())*(vertexMother->GetZ() - mcTrackFirstDaughter->Zv()));
            ((TH1F*)fMotherHistogramArrayExtra[motherType][4])->Fill(vertex_distance);

            Double_t momentum_resolution = TMath::Sqrt((trackB0.Px() - mcTrackB0->Px())*(trackB0.Px() - mcTrackB0->Px()) + (trackB0.Py() - mcTrackB0->Py())*(trackB0.Py() - mcTrackB0->Py()) + (trackB0.Pz() - mcTrackB0->Pz())*(trackB0.Pz() - mcTrackB0->Pz()));
            ((TH1F*)fMotherHistogramArrayExtra[motherType][6])->Fill(momentum_resolution);
          }

          if(!bSameSign)
          {
            // We fill the histograms
            histType = 0;
            FillDStarAndB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
            
            if(isDesiredCandidate)
            {
              histType = 1;
              FillDStarAndB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
            }
          }


          // We apply cuts
          Bool_t cutMother = kFALSE;

          Bool_t bCutArray[97] = {0};
          Int_t numberOfCuts = 97;
          Int_t cutReturnValue = fCuts->IsSelected(&trackB0, 0, aodEvent, bCutArray);
          if(cutReturnValue == -1) cutMother = kTRUE;
          if(cutReturnValue == 0) cutMother = kTRUE;


          // We save information about the cuts
          TString histName = "";
          Double_t invariantMassMother = trackB0.InvMass(2,prongs);
          Double_t pdgMassMother=TDatabasePDG::Instance()->GetParticle(511)->Mass();
          Double_t massWindow = fHistMassWindow; //GeV/c^2
          if(fGetCutInfo == kTRUE)
          {
            for (Int_t n = 0; n < 97; ++n)
            {
              if(bCutArray[n] == kTRUE){
                if(isDesiredCandidate){
                  ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(n+1);
                } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(n+1);
                cutMother = kTRUE;
              }
            }
         
            if (TMath::Abs(invariantMassMother-pdgMassMother)<massWindow){
              for (Int_t l = 0; l < numberOfCuts; ++l) //total
              {
                if(bCutArray[l] == kFALSE) continue;
                for (Int_t j = 0; j < numberOfCuts; ++j)
                {
                  if(bCutArray[j] == kFALSE) continue;
                  if(isDesiredCandidate == kFALSE) histName ="cutEffectBackground";
                  if(isDesiredCandidate == kTRUE) histName ="cutEffectSignal";
                  ((TH2I*)(fOutputB0MC->FindObject(histName)))->Fill(l,j);
                }
              }

              for (Int_t l = 0; l < numberOfCuts; ++l) //unique
              {
                if(bCutArray[l] == kFALSE) continue;
                Bool_t bFill = kTRUE;
                for (Int_t j = 0; j < numberOfCuts; ++j)
                {
                  if(l==j) continue;
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
                  ((TH1I*)(fOutputB0MC->FindObject(histName)))->Fill(l);
                }
              }
            }
          }


          if(!isDesiredCandidate && fQuickSignalAnalysis == 1) cutMother = kTRUE;

          if(cutMother)
          {
            if(isDesiredCandidate)
            {
              ((TH1F*)fMotherHistogramArrayExtra[motherType][1])->Fill(0);
            } else ((TH1F*)fMotherHistogramArrayExtra[motherType][0])->Fill(0);
            delete vertexMother; vertexMother = nullptr; 
            delete vertexDStar; vertexDStar = nullptr; 
            delete trackB0PionRotated; trackB0PionRotated = nullptr; 
            continue;
          }

          // We save the DCA information
          TString name_dca_D0_DStarPion ="dca_D0_DStarPion";
          ((TH1F*)(fOutputB0MC->FindObject(name_dca_D0_DStarPion)))->Fill(dcaDStarPionD0);

          TString name_dca_D0_B0Pion ="dca_D0_B0Pion";
          ((TH1F*)(fOutputB0MC->FindObject(name_dca_D0_B0Pion)))->Fill(dcaB0PionD0);

          TString name_dca_DStarPion_B0Pion ="dca_DStarPion_B0Pion";
          ((TH1F*)(fOutputB0MC->FindObject(name_dca_DStarPion_B0Pion)))->Fill(dcaDStarPionB0Pion);

          TString name_dca_Combined ="dca_Combined";
          ((TH1F*)(fOutputB0MC->FindObject(name_dca_Combined)))->Fill(dcaCombined);

          if(isDesiredCandidate && fUseMCInfo)
          {
            TString name_dca_Signal_D0_DStarPion ="dca_Signal_D0_DStarPion";
            ((TH1F*)(fOutputB0MC->FindObject(name_dca_Signal_D0_DStarPion)))->Fill(dcaDStarPionD0);

            TString name_dca_Signal_D0_B0Pion ="dca_Signal_D0_B0Pion";
            ((TH1F*)(fOutputB0MC->FindObject(name_dca_Signal_D0_B0Pion)))->Fill(dcaB0PionD0);

            TString name_dca_Signal_DStarPion_B0Pion ="dca_Signal_DStarPion_B0Pion";
            ((TH1F*)(fOutputB0MC->FindObject(name_dca_Signal_DStarPion_B0Pion)))->Fill(dcaDStarPionB0Pion);

            TString name_dca_Signal_Combined ="dca_Signal_Combined";
            ((TH1F*)(fOutputB0MC->FindObject(name_dca_Signal_Combined)))->Fill(dcaCombined);
          }


          // Background analysis
          Bool_t bIsCorrelatedBackground = kFALSE;
          Bool_t bIsCorrelatedBackground511 = kFALSE;
          if(!bSameSign && fCheckBackground && fUseMCInfo && !isDesiredCandidate && iRot == 0) 
          {
            Int_t mcLabelB0Pion = trackB0PionRotated->GetLabel();
            Int_t mcLabelDStarPion = trackFirstDaughter->GetLabel();
            Int_t mcLabelD0first = ((AliAODTrack*)trackSecondDaughter->GetDaughter(0))->GetLabel();
            Int_t mcLabelD0second = ((AliAODTrack*)trackSecondDaughter->GetDaughter(1))->GetLabel();

            if(mcLabelB0Pion >= 0 && mcLabelDStarPion >= 0 && mcLabelD0first >= 0 && mcLabelD0second >= 0)
            {
              AliAODMCParticle * mcB0Pion = (AliAODMCParticle*)mcTrackArray->At(mcLabelB0Pion);
              AliAODMCParticle * mcDStarPion = (AliAODMCParticle*)mcTrackArray->At(mcLabelDStarPion);
              AliAODMCParticle * mcD0first = (AliAODMCParticle*)mcTrackArray->At(mcLabelD0first);
              AliAODMCParticle * mcD0second = (AliAODMCParticle*)mcTrackArray->At(mcLabelD0second);

              if(mcB0Pion)
              {
                Int_t iterator = 0;
                while(mcB0Pion->GetMother() >= 0)
                {
                  mcB0Pion = (AliAODMCParticle*)mcTrackArray->At(mcB0Pion->GetMother());
                  fillthis="particle_pdgB0Pion";
                  ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(mcB0Pion->GetPdgCode()),iterator++);
                } 
              }

              if(mcDStarPion)
              {
                Int_t iterator = 0;
                while(mcDStarPion->GetMother() >= 0)
                {
                  mcDStarPion = (AliAODMCParticle*)mcTrackArray->At(mcDStarPion->GetMother());
                  fillthis="particle_pdgB0Pion";
                  ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(mcDStarPion->GetPdgCode()),iterator++);
                } 
              }

              if(mcD0first)
              {
                Int_t iterator = 0;
                while(mcD0first->GetMother() >= 0)
                {
                  mcD0first = (AliAODMCParticle*)mcTrackArray->At(mcD0first->GetMother());
                  fillthis="particle_pdgD0First";
                  ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(mcD0first->GetPdgCode()),iterator++);
                }
              }

              if(mcD0second)
              {
                Int_t iterator = 0;
                while(mcD0second->GetMother() >= 0)
                {
                  mcD0second = (AliAODMCParticle*)mcTrackArray->At(mcD0second->GetMother());
                  fillthis="particle_pdgD0Second";
                  ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(mcD0second->GetPdgCode()),iterator++);
                }
              }

              mcB0Pion = (AliAODMCParticle*)mcTrackArray->At(mcLabelB0Pion);
              mcDStarPion = (AliAODMCParticle*)mcTrackArray->At(mcLabelDStarPion);
              mcD0first = (AliAODMCParticle*)mcTrackArray->At(mcLabelD0first);
              mcD0second = (AliAODMCParticle*)mcTrackArray->At(mcLabelD0second);

              if(mcB0Pion && mcDStarPion && mcD0first && mcD0second)
              {
                // B -> DStar + pion -> D0 + pion
                if(mcD0first->GetMother() == mcD0second->GetMother() && mcD0first->GetMother() >= 0)
                {
                  AliAODMCParticle * D0Mother = (AliAODMCParticle*)mcTrackArray->At(mcD0first->GetMother());
                  AliAODMCParticle * D0GrandMother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetMother());
                  if(D0Mother->GetMother() == mcDStarPion->GetMother() && D0Mother->GetMother() >= 0 && D0GrandMother->GetMother() == mcB0Pion->GetMother() && D0GrandMother->GetMother() >= 0)
                  {
                    AliAODMCParticle * finalMother = (AliAODMCParticle*)mcTrackArray->At(D0GrandMother->GetMother());
                    fillthis="particle_pdgAll";
                    ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(finalMother->GetPdgCode()));
                    fillthis="particle_pdgAllInvMass";
                    ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(finalMother->GetPdgCode()),invariantMassMother);
                    if(TMath::Abs(finalMother->GetPdgCode())==511)
                    {
                      bIsCorrelatedBackground = kTRUE;
                      for (Int_t iDaughter = 0; iDaughter < finalMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterfinalMother = (AliAODMCParticle*)mcTrackArray->At(finalMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep511a";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter,TMath::Abs(daughterfinalMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0GrandMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0GrandMother = (AliAODMCParticle*)mcTrackArray->At(D0GrandMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep511a";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+10,TMath::Abs(daughterD0GrandMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0Mother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0Mother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep511a";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+20,TMath::Abs(daughterD0Mother->GetPdgCode()));
                      }
                    }
                    if(TMath::Abs(finalMother->GetPdgCode())==521)
                    {
                      bIsCorrelatedBackground = kTRUE;
                      for (Int_t iDaughter = 0; iDaughter < finalMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterfinalMother = (AliAODMCParticle*)mcTrackArray->At(finalMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep521a";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter,TMath::Abs(daughterfinalMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0GrandMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0GrandMother = (AliAODMCParticle*)mcTrackArray->At(D0GrandMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep521a";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+10,TMath::Abs(daughterD0GrandMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0Mother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0Mother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep521a";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+20,TMath::Abs(daughterD0Mother->GetPdgCode()));
                      }
                    }
                  }
                }

                // B -> D0 + (B0) pion
                if(mcD0first->GetMother() == mcD0second->GetMother() && mcD0first->GetMother() >= 0)
                {
                  AliAODMCParticle * D0Mother = (AliAODMCParticle*)mcTrackArray->At(mcD0first->GetMother());
                  if(D0Mother->GetMother() == mcB0Pion->GetMother() && D0Mother->GetMother() >= 0)
                  {
                    AliAODMCParticle * finalMother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetMother());
                    fillthis="particle_pdgAllSecond";
                    ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(finalMother->GetPdgCode()));
                    fillthis="particle_pdgAllInvMassSecond";
                    ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(finalMother->GetPdgCode()),invariantMassMother);
                    if(TMath::Abs(finalMother->GetPdgCode())==511)
                    {
                      bIsCorrelatedBackground = kTRUE;
                      for (Int_t iDaughter = 0; iDaughter < finalMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterfinalMother = (AliAODMCParticle*)mcTrackArray->At(finalMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgOneStep511a";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter,TMath::Abs(daughterfinalMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0Mother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0Mother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgOneStep511a";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+10,TMath::Abs(daughterD0Mother->GetPdgCode()));
                      }
                    }
                    if(TMath::Abs(finalMother->GetPdgCode())==521)
                    {
                      if(fUpgradeSetting == 2 || fUpgradeSetting == 3) 
                      {
                        delete vertexMother; vertexMother = nullptr; 
                        delete vertexDStar; vertexDStar = nullptr; 
                        delete trackB0PionRotated; trackB0PionRotated = nullptr; 
                        continue;
                      }
                      bIsCorrelatedBackground = kTRUE;
                      for (Int_t iDaughter = 0; iDaughter < finalMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterfinalMother = (AliAODMCParticle*)mcTrackArray->At(finalMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgOneStep521a";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter,TMath::Abs(daughterfinalMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0Mother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0Mother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgOneStep521a";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+10,TMath::Abs(daughterD0Mother->GetPdgCode()));
                      }
                    }
                  }
                }

                // B -> D0 + (DStar) pion
                if(mcD0first->GetMother() == mcD0second->GetMother() && mcD0first->GetMother() >= 0)
                {
                  AliAODMCParticle * D0Mother = (AliAODMCParticle*)mcTrackArray->At(mcD0first->GetMother());
                  if(D0Mother->GetMother() == mcDStarPion->GetMother() && D0Mother->GetMother() >= 0)
                  {
                    AliAODMCParticle * finalMother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetMother());
                    fillthis="particle_pdgAllThird";
                    ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(finalMother->GetPdgCode()));
                    fillthis="particle_pdgAllInvMassThird";
                    ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(finalMother->GetPdgCode()),invariantMassMother);
                    if(TMath::Abs(finalMother->GetPdgCode())==511)
                    {
                      bIsCorrelatedBackground = kTRUE;
                      for (Int_t iDaughter = 0; iDaughter < finalMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterfinalMother = (AliAODMCParticle*)mcTrackArray->At(finalMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgOneStep511b";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter,TMath::Abs(daughterfinalMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0Mother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0Mother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgOneStep511b";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+10,TMath::Abs(daughterD0Mother->GetPdgCode()));
                      }
                    }
                    if(TMath::Abs(finalMother->GetPdgCode())==521)
                    {
                      bIsCorrelatedBackground = kTRUE;
                      for (Int_t iDaughter = 0; iDaughter < finalMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterfinalMother = (AliAODMCParticle*)mcTrackArray->At(finalMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgOneStep521b";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter,TMath::Abs(daughterfinalMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0Mother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0Mother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgOneStep521b";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+10,TMath::Abs(daughterD0Mother->GetPdgCode()));
                      }
                    }
                  }
                }

                // B -> X + (B0) pion -> D0 
                if(mcD0first->GetMother() == mcD0second->GetMother() && mcD0first->GetMother() >= 0)
                {
                  AliAODMCParticle * D0Mother = (AliAODMCParticle*)mcTrackArray->At(mcD0first->GetMother());
                  AliAODMCParticle * D0GrandMother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetMother());
                  if(D0GrandMother->GetMother() == mcB0Pion->GetMother() && D0GrandMother->GetMother() >= 0)
                  {
                    AliAODMCParticle * finalMother = (AliAODMCParticle*)mcTrackArray->At(D0GrandMother->GetMother());
                    fillthis="particle_pdgAllFourth";
                    ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(finalMother->GetPdgCode()));
                    fillthis="particle_pdgAllInvMassFourth";
                    ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(finalMother->GetPdgCode()),invariantMassMother);
                    if(TMath::Abs(finalMother->GetPdgCode())==511)
                    {
                      if(fUpgradeSetting == 1 || fUpgradeSetting == 3)
                      {
                        delete vertexMother; vertexMother = nullptr; 
                        delete vertexDStar; vertexDStar = nullptr; 
                        delete trackB0PionRotated; trackB0PionRotated = nullptr; 
                        continue;
                      }
                      bIsCorrelatedBackground = kTRUE;
                      bIsCorrelatedBackground511 = kTRUE;
                      for (Int_t iDaughter = 0; iDaughter < finalMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterfinalMother = (AliAODMCParticle*)mcTrackArray->At(finalMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep511b";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter,TMath::Abs(daughterfinalMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0GrandMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0GrandMother = (AliAODMCParticle*)mcTrackArray->At(D0GrandMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep511b";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+10,TMath::Abs(daughterD0GrandMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0Mother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0Mother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep511b";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+20,TMath::Abs(daughterD0Mother->GetPdgCode()));
                      }
                    }
                    if(TMath::Abs(finalMother->GetPdgCode())==521)
                    {
                      bIsCorrelatedBackground = kTRUE;
                      for (Int_t iDaughter = 0; iDaughter < finalMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterfinalMother = (AliAODMCParticle*)mcTrackArray->At(finalMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep521b";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter,TMath::Abs(daughterfinalMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0GrandMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0GrandMother = (AliAODMCParticle*)mcTrackArray->At(D0GrandMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep521b";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+10,TMath::Abs(daughterD0GrandMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0Mother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0Mother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep521b";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+20,TMath::Abs(daughterD0Mother->GetPdgCode()));
                      }
                    }
                  }
                }

                // B -> X + (DStar) pion -> D0 
                if(mcD0first->GetMother() == mcD0second->GetMother() && mcD0first->GetMother() >= 0)
                {
                  AliAODMCParticle * D0Mother = (AliAODMCParticle*)mcTrackArray->At(mcD0first->GetMother());
                  AliAODMCParticle * D0GrandMother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetMother());
                  if(D0GrandMother->GetMother() == mcDStarPion->GetMother() && D0GrandMother->GetMother() >= 0)
                  {
                    AliAODMCParticle * finalMother = (AliAODMCParticle*)mcTrackArray->At(D0GrandMother->GetMother());
                    fillthis="particle_pdgAllFifth";
                    ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(finalMother->GetPdgCode()));
                    fillthis="particle_pdgAllInvMassFifth";
                    ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(TMath::Abs(finalMother->GetPdgCode()),invariantMassMother);
                    if(TMath::Abs(finalMother->GetPdgCode())==511)
                    {
                      bIsCorrelatedBackground = kTRUE;
                      for (Int_t iDaughter = 0; iDaughter < finalMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterfinalMother = (AliAODMCParticle*)mcTrackArray->At(finalMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep511c";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter,TMath::Abs(daughterfinalMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0GrandMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0GrandMother = (AliAODMCParticle*)mcTrackArray->At(D0GrandMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep511c";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+10,TMath::Abs(daughterD0GrandMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0Mother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0Mother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep511c";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+20,TMath::Abs(daughterD0Mother->GetPdgCode()));
                      }
                    }
                    if(TMath::Abs(finalMother->GetPdgCode())==521)
                    {
                      bIsCorrelatedBackground = kTRUE;
                      for (Int_t iDaughter = 0; iDaughter < finalMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterfinalMother = (AliAODMCParticle*)mcTrackArray->At(finalMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep521c";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter,TMath::Abs(daughterfinalMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0GrandMother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0GrandMother = (AliAODMCParticle*)mcTrackArray->At(D0GrandMother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep521c";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+10,TMath::Abs(daughterD0GrandMother->GetPdgCode()));
                      }
                      for (Int_t iDaughter = 0; iDaughter < D0Mother->GetNDaughters(); ++iDaughter) //will work up to ten daughters
                      {
                        AliAODMCParticle* daughterD0Mother = (AliAODMCParticle*)mcTrackArray->At(D0Mother->GetDaughterLabel(0)+iDaughter);
                        fillthis="particle_daughterPdgTwoStep521c";
                        ((TH2F*)(fOutputB0MC->FindObject(fillthis)))->Fill(iDaughter+20,TMath::Abs(daughterD0Mother->GetPdgCode()));
                      }
                    }
                  }
                }
              }
            }
          }
        
          if(!bSameSign)
          {
            histType = 2;
            FillDStarAndB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
            if(fUseMCInfo && isDesiredCandidate)
            {
              //fill mc histograms
              histType = 3;
              FillDStarAndB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
            }
          }

          AliAODRecoDecayHF2Prong* selectedDStar = (AliAODRecoDecayHF2Prong*)trackB0.GetDaughter(1);
          AliAODRecoDecayHF2Prong* selectedD0 = (AliAODRecoDecayHF2Prong*)selectedDStar->GetDaughter(1);

          // We fill the final cut histograms with the candidates that have an invariant mass close to the PDG value. This way the background we use for optimizing the cuts will not be contaminated with candidates that don't affect the signal region.
          massWindow =  fHistMassWindow; // GeV/c^2
          if (TMath::Abs(invariantMassMother-pdgMassMother)<massWindow)
          {
            if(!bSameSign) 
            {
              FillFinalTrackHistograms(&trackB0, isDesiredCandidate, mcTrackArray);
              if(!isDesiredCandidate)
              {
                motherType = 0; histType = 4; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histType, pdgD0);
                motherType = 1; histType = 4; FillDStarAndB0Histograms(selectedDStar, primaryVertex, bz, motherType, histType);
                motherType = 2; histType = 4; FillDStarAndB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
              }
              if(isDesiredCandidate)
              {
                motherType = 0; histType = 5; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histType, pdgD0);
                motherType = 1; histType = 5; FillDStarAndB0Histograms(selectedDStar, primaryVertex, bz, motherType, histType);
                motherType = 2; histType = 5; FillDStarAndB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
              }
            }
          }

          // Here we fill the histograms per pt bin and apply the same sign method
          TString ptBinMother = "";
          Int_t ptBin = fCuts->PtBin(trackB0.Pt());
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
                motherType = 0; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histType, pdgD0);
                motherType = 1; FillDStarAndB0Histograms(selectedDStar, primaryVertex, bz, motherType, histType);
                motherType = 2; FillDStarAndB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
                motherType = 3; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histTypeD0, pdgD0);
                motherType = 4; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histTypeD0DStar, pdgD0);
                motherType = 5; FillDStarAndB0Histograms(selectedDStar, primaryVertex, bz, motherType, histTypeDStar);
              }

              if(isDesiredCandidate)
              {
                histType += 1;
                motherType = 0; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histType, pdgD0);
                motherType = 1; FillDStarAndB0Histograms(selectedDStar, primaryVertex, bz, motherType, histType);
                motherType = 2; FillDStarAndB0Histograms(&trackB0, primaryVertex, bz, motherType, histType);
                motherType = 3; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histTypeD0 + 1, pdgD0);
                motherType = 4; FillD0Histograms(selectedD0, primaryVertex, bz, motherType, histTypeD0DStar + 1, pdgD0);
                motherType = 5; FillDStarAndB0Histograms(selectedDStar, primaryVertex, bz, motherType, histTypeDStar + 1);
              }
            }
          }

          Double_t invmassDelta = DeltaInvMassB0Kpipipi(&trackB0);
          if(bSameSign && iRot == 0)
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
            if(iRot == 0)
            {
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

              fillthis = "invariantMassB0Signal_BA";
              if(isDesiredCandidate) ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);

              fillthis = "invariantMassB0Correlated_BA";
              if(bIsCorrelatedBackground) ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);

              fillthis = "invariantMassB0Background_BA";
              if(!isDesiredCandidate && !bIsCorrelatedBackground) ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);

              if(bIsCorrelatedBackground511)
              {
                fillthis="invariantMassB0_correlated511";
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                fillthis="invariantMassB0";
                fillthis += ptBinMother + "_correlated511";
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
              }
              
              if(!isDesiredCandidate && !bIsInjected)
              {
                TString signName = "_HIJING_Background";
                fillthis="invariantMassB0";
                fillthis += signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                fillthis="invariantMassB0";
                fillthis += ptBinMother + signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
              }
              if(isDesiredCandidate && !bIsInjected)
              {
                TString signName = "_HIJING_Signal";
                fillthis="invariantMassB0";
                fillthis += signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                fillthis="invariantMassB0";
                fillthis += ptBinMother + signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
              }
            }
            else
            {
              TString signName = "_Background_rotation";
              fillthis="invariantMassB0";
              fillthis += signName;
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
              fillthis="invariantMassB0";
              fillthis += ptBinMother + signName;
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
              if(!isDesiredCandidate && !bIsInjected)
              {
                signName = "_HIJING_Background_rotation";
                fillthis="invariantMassB0";
                fillthis += signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                fillthis="invariantMassB0";
                fillthis += ptBinMother + signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
              }
            }
          }

          if(trackB0.Pt() > 6.0)
          {
            TString broadptBinMother = "_ptbin_6_to_inf";
            if(bSameSign  && iRot == 0)
            {
              fillthis="invariantMassB0";
              fillthis += broadptBinMother + "_SameSign";
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
              fillthis="invariantMassB0";
              fillthis += broadptBinMother + "_SignSum";
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,-1);
            }
            if(!bSameSign)
            {
              if(iRot == 0)
              {
                fillthis="invariantMassB0";
                fillthis += broadptBinMother;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                fillthis="invariantMassB0";
                fillthis += broadptBinMother + "_SignSum";
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,1);
                if(!isDesiredCandidate && !bIsInjected)
                {
                  TString signName = "_HIJING_Background";
                  fillthis="invariantMassB0";
                  fillthis += broadptBinMother + signName;
                  ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                }
                if(isDesiredCandidate && !bIsInjected)
                {
                  TString signName = "_HIJING_Signal";
                  fillthis="invariantMassB0";
                  fillthis += broadptBinMother + signName;
                  ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                }
              }
              else
              {
                TString signName = "_Background_rotation";
                fillthis="invariantMassB0";
                fillthis += broadptBinMother + signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                if(!isDesiredCandidate && !bIsInjected)
                {
                  signName = "_HIJING_Background_rotation";
                  fillthis="invariantMassB0";
                  fillthis += broadptBinMother + signName;
                  ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                }
              }
            }
          }

          if(trackB0.Pt() > 3.0)
          {
            TString broadptBinMother = "_ptbin_3_to_inf";
            if(bSameSign  && iRot == 0)
            {
              fillthis="invariantMassB0";
              fillthis += broadptBinMother + "_SameSign";
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
              fillthis="invariantMassB0";
              fillthis += broadptBinMother + "_SignSum";
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,-1);
            }
            if(!bSameSign)
            {
              if(iRot == 0)
              {
                fillthis="invariantMassB0";
                fillthis += broadptBinMother;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                fillthis="invariantMassB0";
                fillthis += broadptBinMother + "_SignSum";
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,1);
                if(!isDesiredCandidate && !bIsInjected)
                {
                  TString signName = "_HIJING_Background";
                  fillthis="invariantMassB0";
                  fillthis += broadptBinMother + signName;
                  ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                }
                if(isDesiredCandidate && !bIsInjected)
                {
                  TString signName = "_HIJING_Signal";
                  fillthis="invariantMassB0";
                  fillthis += broadptBinMother + signName;
                  ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                }
              }
              else
              {
                TString signName = "_Background_rotation";
                fillthis="invariantMassB0";
                fillthis += broadptBinMother + signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                if(!isDesiredCandidate && !bIsInjected)
                {
                  signName = "_HIJING_Background_rotation";
                  fillthis="invariantMassB0";
                  fillthis += broadptBinMother + signName;
                  ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
                }
              }
            }
          }

          // //fine binning
          // if(bSameSign  && iRot == 0)
          // {
          //   fillthis="fineBin_invariantMassB0";
          //   fillthis += "_SameSign";
          //   ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //   fillthis="fineBin_invariantMassB0";
          //   fillthis += ptBinMother + "_SameSign";
          //   ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //   fillthis="fineBin_invariantMassB0";
          //   fillthis += "_SignSum";
          //   ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,-1);
          //   fillthis="fineBin_invariantMassB0";
          //   fillthis += ptBinMother + "_SignSum";
          //   ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,-1);
          // }
          // if(!bSameSign)
          // {
          //   if(iRot == 0)
          //   {
          //     fillthis="fineBin_invariantMassB0";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //     fillthis="fineBin_invariantMassB0";
          //     fillthis += ptBinMother;
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //     fillthis="fineBin_invariantMassB0";
          //     fillthis += "_SignSum";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,1);
          //     fillthis="fineBin_invariantMassB0";
          //     fillthis += ptBinMother + "_SignSum";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,1);

          //     if(!isDesiredCandidate && !bIsInjected)
          //     {
          //       TString signName = "_HIJING_Background";
          //       fillthis="fineBin_invariantMassB0";
          //       fillthis += signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       fillthis="fineBin_invariantMassB0";
          //       fillthis += ptBinMother + signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //     }
          //     if(isDesiredCandidate && !bIsInjected)
          //     {
          //       TString signName = "_HIJING_Signal";
          //       fillthis="fineBin_invariantMassB0";
          //       fillthis += signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       fillthis="fineBin_invariantMassB0";
          //       fillthis += ptBinMother + signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //     }
          //   }
          //   else
          //   {
          //     TString signName = "_Background_rotation";
          //     fillthis="fineBin_invariantMassB0";
          //     fillthis += signName;
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //     fillthis="fineBin_invariantMassB0";
          //     fillthis += ptBinMother + signName;
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //     if(!isDesiredCandidate && !bIsInjected)
          //     {
          //       signName = "_HIJING_Background_rotation";
          //       fillthis="fineBin_invariantMassB0";
          //       fillthis += signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       fillthis="fineBin_invariantMassB0";
          //       fillthis += ptBinMother + signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //     }
          //   }
          // }

          // if(trackB0.Pt() > 6.0)
          // {
          //   TString broadptBinMother = "_ptbin_6_to_inf";
          //   if(bSameSign && iRot == 0)
          //   {
          //     fillthis="fineBin_invariantMassB0";
          //     fillthis += broadptBinMother + "_SameSign";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //     fillthis="fineBin_invariantMassB0";
          //     fillthis += broadptBinMother + "_SignSum";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,-1);
          //   }
          //   if(!bSameSign)
          //   {
          //     if(iRot == 0)
          //     {
          //       fillthis="fineBin_invariantMassB0";
          //       fillthis += broadptBinMother;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       fillthis="fineBin_invariantMassB0";
          //       fillthis += broadptBinMother + "_SignSum";
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,1);
          //       if(!isDesiredCandidate && !bIsInjected)
          //       {
          //         TString signName = "_HIJING_Background";
          //         fillthis="fineBin_invariantMassB0";
          //         fillthis += broadptBinMother + signName;
          //         ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       }
          //       if(isDesiredCandidate && !bIsInjected)
          //       {
          //         TString signName = "_HIJING_Signal";
          //         fillthis="fineBin_invariantMassB0";
          //         fillthis += broadptBinMother + signName;
          //         ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       }
          //     }
          //     else
          //     {
          //       TString signName = "_Background_rotation";
          //       fillthis="fineBin_invariantMassB0";
          //       fillthis += broadptBinMother + signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       if(!isDesiredCandidate && !bIsInjected)
          //       {
          //         signName = "_HIJING_Background_rotation";
          //         fillthis="fineBin_invariantMassB0";
          //         fillthis += broadptBinMother + signName;
          //         ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       }
          //     }
          //   }
          // }

          // if(trackB0.Pt() > 3.0)
          // {
          //   TString broadptBinMother = "_ptbin_3_to_inf";
          //   if(bSameSign && iRot == 0)
          //   {
          //     fillthis="fineBin_invariantMassB0";
          //     fillthis += broadptBinMother + "_SameSign";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //     fillthis="fineBin_invariantMassB0";
          //     fillthis += broadptBinMother + "_SignSum";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,-1);
          //   }
          //   if(!bSameSign)
          //   {
          //     if(iRot == 0)
          //     {
          //       fillthis="fineBin_invariantMassB0";
          //       fillthis += broadptBinMother;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       fillthis="fineBin_invariantMassB0";
          //       fillthis += broadptBinMother + "_SignSum";
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother,1);
          //       if(!isDesiredCandidate && !bIsInjected)
          //       {
          //         TString signName = "_HIJING_Background";
          //         fillthis="fineBin_invariantMassB0";
          //         fillthis += broadptBinMother + signName;
          //         ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       }
          //       if(isDesiredCandidate && !bIsInjected)
          //       {
          //         TString signName = "_HIJING_Signal";
          //         fillthis="fineBin_invariantMassB0";
          //         fillthis += broadptBinMother + signName;
          //         ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       }
          //     }
          //     else
          //     {
          //       TString signName = "_Background_rotation";
          //       fillthis="fineBin_invariantMassB0";
          //       fillthis += broadptBinMother + signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       if(!isDesiredCandidate && !bIsInjected)
          //       {
          //         signName = "_HIJING_Background_rotation";
          //         fillthis="fineBin_invariantMassB0";
          //         fillthis += broadptBinMother + signName;
          //         ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invariantMassMother);
          //       }
          //     }
          //   }
          // }

          //deltamass
          if(bSameSign && iRot == 0)
          {
            fillthis="deltainvariantMassB0";
            fillthis += "_SameSign";
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
            fillthis="deltainvariantMassB0";
            fillthis += ptBinMother + "_SameSign";
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
            fillthis="deltainvariantMassB0";
            fillthis += "_SignSum";
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,-1);
            fillthis="deltainvariantMassB0";
            fillthis += ptBinMother + "_SignSum";
            ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,-1);
          }
          if(!bSameSign)
          {
            if(iRot == 0)
            {
              fillthis="deltainvariantMassB0";
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
              fillthis="deltainvariantMassB0";
              fillthis += ptBinMother;
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
              fillthis="deltainvariantMassB0";
              fillthis += "_SignSum";
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,1);
              fillthis="deltainvariantMassB0";
              fillthis += ptBinMother + "_SignSum";
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,1);
              
              if(bIsCorrelatedBackground511)
              {
                fillthis="deltainvariantMassB0_correlated511";
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                fillthis="deltainvariantMassB0";
                fillthis += ptBinMother + "_correlated511" ;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
              }

              if(!isDesiredCandidate && !bIsInjected)
              {
                TString signName = "_HIJING_Background";
                fillthis="deltainvariantMassB0";
                fillthis += signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                fillthis="deltainvariantMassB0";
                fillthis += ptBinMother + signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
              }
              if(isDesiredCandidate && !bIsInjected)
              {
                TString signName = "_HIJING_Signal";
                fillthis="deltainvariantMassB0";
                fillthis += signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                fillthis="deltainvariantMassB0";
                fillthis += ptBinMother + signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
              }
            }
            else
            {
              TString signName = "_Background_rotation";
              fillthis="deltainvariantMassB0";
              fillthis += signName;
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
              fillthis="deltainvariantMassB0";
              fillthis += ptBinMother + signName;
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
              if(!isDesiredCandidate && !bIsInjected)
              {
                signName = "_HIJING_Background_rotation";
                fillthis="deltainvariantMassB0";
                fillthis += signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                fillthis="deltainvariantMassB0";
                fillthis += ptBinMother + signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
              }
            }
          }

          if(trackB0.Pt() > 6.0)
          {
            TString broadptBinMother = "_ptbin_6_to_inf";
            if(bSameSign && iRot == 0)
            {
              fillthis="deltainvariantMassB0";
              fillthis += broadptBinMother + "_SameSign";
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
              fillthis="deltainvariantMassB0";
              fillthis += broadptBinMother + "_SignSum";
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,-1);
            }
            if(!bSameSign)
            {
              if(iRot == 0)
              {
                fillthis="deltainvariantMassB0";
                fillthis += broadptBinMother;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                fillthis="deltainvariantMassB0";
                fillthis += broadptBinMother + "_SignSum";
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,1);
                if(!isDesiredCandidate && !bIsInjected)
                {
                  TString signName = "_HIJING_Background";
                  fillthis="deltainvariantMassB0";
                  fillthis += broadptBinMother + signName;
                  ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                }
                if(isDesiredCandidate && !bIsInjected)
                {
                  TString signName = "_HIJING_Signal";
                  fillthis="deltainvariantMassB0";
                  fillthis += broadptBinMother + signName;
                  ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                }
              }
              else
              {
                TString signName = "_Background_rotation";
                fillthis="deltainvariantMassB0";
                fillthis += broadptBinMother + signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                if(!isDesiredCandidate && !bIsInjected)
                {
                  signName = "_HIJING_Background_rotation";
                  fillthis="deltainvariantMassB0";
                  fillthis += broadptBinMother + signName;
                  ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                }
              }
            }
          }

          if(trackB0.Pt() > 3.0)
          {
            TString broadptBinMother = "_ptbin_3_to_inf";
            if(bSameSign && iRot == 0)
            {
              fillthis="deltainvariantMassB0";
              fillthis += broadptBinMother + "_SameSign";
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
              fillthis="deltainvariantMassB0";
              fillthis += broadptBinMother + "_SignSum";
              ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,-1);
            }
            if(!bSameSign)
            {
              if(iRot == 0)
              {
                fillthis="deltainvariantMassB0";
                fillthis += broadptBinMother;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                fillthis="deltainvariantMassB0";
                fillthis += broadptBinMother + "_SignSum";
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,1);
                if(!isDesiredCandidate && !bIsInjected)
                {
                  TString signName = "_HIJING_Background";
                  fillthis="deltainvariantMassB0";
                  fillthis += broadptBinMother + signName;
                  ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                }
                if(isDesiredCandidate && !bIsInjected)
                {
                  TString signName = "_HIJING_Signal";
                  fillthis="deltainvariantMassB0";
                  fillthis += broadptBinMother + signName;
                  ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                }
              }
              else
              {
                TString signName = "_Background_rotation";
                fillthis="deltainvariantMassB0";
                fillthis += broadptBinMother + signName;
                ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                if(!isDesiredCandidate && !bIsInjected)
                {
                  signName = "_HIJING_Background_rotation";
                  fillthis="deltainvariantMassB0";
                  fillthis += broadptBinMother + signName;
                  ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
                }
              }
            }
          }

          // //fine binning
          // if(bSameSign && iRot == 0)
          // {
          //   fillthis="fineBin_deltainvariantMassB0";
          //   fillthis += "_SameSign";
          //   ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //   fillthis="fineBin_deltainvariantMassB0";
          //   fillthis += ptBinMother + "_SameSign";
          //   ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //   fillthis="fineBin_deltainvariantMassB0";
          //   fillthis += "_SignSum";
          //   ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,-1);
          //   fillthis="fineBin_deltainvariantMassB0";
          //   fillthis += ptBinMother + "_SignSum";
          //   ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,-1);
          // }
          // if(!bSameSign)
          // {
          //   if(iRot == 0)
          //   {
          //     fillthis="fineBin_deltainvariantMassB0";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //     fillthis="fineBin_deltainvariantMassB0";
          //     fillthis += ptBinMother;
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //     fillthis="fineBin_deltainvariantMassB0";
          //     fillthis += "_SignSum";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,1);
          //     fillthis="fineBin_deltainvariantMassB0";
          //     fillthis += ptBinMother + "_SignSum";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,1);

          //     if(!isDesiredCandidate && !bIsInjected)
          //     {
          //       TString signName = "_HIJING_Background";
          //       fillthis="fineBin_deltainvariantMassB0";
          //       fillthis += signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       fillthis="fineBin_deltainvariantMassB0";
          //       fillthis += ptBinMother + signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //     }
          //     if(isDesiredCandidate && !bIsInjected)
          //     {
          //       TString signName = "_HIJING_Signal";
          //       fillthis="fineBin_deltainvariantMassB0";
          //       fillthis += signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       fillthis="fineBin_deltainvariantMassB0";
          //       fillthis += ptBinMother + signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //     }
          //   }
          //   else
          //   {
          //     TString signName = "_Background_rotation";
          //     fillthis="fineBin_deltainvariantMassB0";
          //     fillthis += signName;
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //     fillthis="fineBin_deltainvariantMassB0";
          //     fillthis += ptBinMother + signName;
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //     if(!isDesiredCandidate && !bIsInjected)
          //     {
          //       signName = "_HIJING_Background_rotation";
          //       fillthis="fineBin_deltainvariantMassB0";
          //       fillthis += signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       fillthis="fineBin_deltainvariantMassB0";
          //       fillthis += ptBinMother + signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //     }
          //   }
          // }

          // if(trackB0.Pt() > 6.0)
          // {
          //   TString broadptBinMother = "_ptbin_6_to_inf";
          //   if(bSameSign && iRot == 0)
          //   {
          //     fillthis="fineBin_deltainvariantMassB0";
          //     fillthis += broadptBinMother + "_SameSign";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //     fillthis="fineBin_deltainvariantMassB0";
          //     fillthis += broadptBinMother + "_SignSum";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,-1);
          //   }
          //   if(!bSameSign)
          //   {
          //     if(iRot == 0)
          //     {
          //       fillthis="fineBin_deltainvariantMassB0";
          //       fillthis += broadptBinMother;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       fillthis="fineBin_deltainvariantMassB0";
          //       fillthis += broadptBinMother + "_SignSum";
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,1);
          //       if(!isDesiredCandidate && !bIsInjected)
          //       {
          //         TString signName = "_HIJING_Background";
          //         fillthis="fineBin_deltainvariantMassB0";
          //         fillthis += broadptBinMother + signName;
          //         ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       }
          //       if(isDesiredCandidate && !bIsInjected)
          //       {
          //         TString signName = "_HIJING_Signal";
          //         fillthis="fineBin_deltainvariantMassB0";
          //         fillthis += broadptBinMother + signName;
          //         ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       }
          //     }
          //     else
          //     {
          //       TString signName = "_Background_rotation";
          //       fillthis="fineBin_deltainvariantMassB0";
          //       fillthis += broadptBinMother + signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       if(!isDesiredCandidate && !bIsInjected)
          //       {
          //         signName = "_HIJING_Background_rotation";
          //         fillthis="fineBin_deltainvariantMassB0";
          //         fillthis += broadptBinMother + signName;
          //         ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       }
          //     }
          //   }
          // }

          // if(trackB0.Pt() > 3.0)
          // {
          //   TString broadptBinMother = "_ptbin_3_to_inf";
          //   if(bSameSign && iRot == 0)
          //   {
          //     fillthis="fineBin_deltainvariantMassB0";
          //     fillthis += broadptBinMother + "_SameSign";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //     fillthis="fineBin_deltainvariantMassB0";
          //     fillthis += broadptBinMother + "_SignSum";
          //     ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,-1);
          //   }
          //   if(!bSameSign)
          //   {
          //     if(iRot == 0)
          //     {
          //       fillthis="fineBin_deltainvariantMassB0";
          //       fillthis += broadptBinMother;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       fillthis="fineBin_deltainvariantMassB0";
          //       fillthis += broadptBinMother + "_SignSum";
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta,1);
          //       if(!isDesiredCandidate && !bIsInjected)
          //       {
          //         TString signName = "_HIJING_Background";
          //         fillthis="fineBin_deltainvariantMassB0";
          //         fillthis += broadptBinMother + signName;
          //         ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       }
          //       if(isDesiredCandidate && !bIsInjected)
          //       {
          //         TString signName = "_HIJING_Signal";
          //         fillthis="fineBin_deltainvariantMassB0";
          //         fillthis += broadptBinMother + signName;
          //         ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       }
          //     }
          //     else
          //     {
          //       TString signName = "_Background_rotation";
          //       fillthis="fineBin_deltainvariantMassB0";
          //       fillthis += broadptBinMother + signName;
          //       ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       if(!isDesiredCandidate && !bIsInjected)
          //       {
          //         signName = "_HIJING_Background_rotation";
          //         fillthis="fineBin_deltainvariantMassB0";
          //         fillthis += broadptBinMother + signName;
          //         ((TH1F*)(fOutputB0MC->FindObject(fillthis)))->Fill(invmassDelta);
          //       }
          //     }
          //   }
          // }

          delete vertexMother; vertexMother = nullptr; 
          delete vertexDStar; vertexDStar = nullptr; 
          delete trackB0PionRotated; trackB0PionRotated = nullptr; 
        }
      }  
    }

    //Restore Object count 
    //To save space in the table keeping track of all referenced objects, 
    //we reset the object count to what it was at the beginning of the loop.
    TProcessID::SetObjectCount(ObjectNumber);    
  }
  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::FillFinalTrackHistograms(AliAODRecoDecayHF2Prong * selectedB0, Bool_t isDesiredCandidate,TClonesArray * mcTrackArray){

  //In this function we fill histograms with the properties of all the daughters of our selected signal candidate

  AliAODTrack* selectedB0Pion = (AliAODTrack*)selectedB0->GetDaughter(0);
  AliAODRecoDecayHF2Prong* selectedDStar = (AliAODRecoDecayHF2Prong*)selectedB0->GetDaughter(1);

  AliAODTrack* selectedDStarPion = (AliAODTrack*)selectedDStar->GetDaughter(0);
  AliAODRecoDecayHF2Prong* selectedD0 = (AliAODRecoDecayHF2Prong*)selectedDStar->GetDaughter(1);

  AliAODTrack* selectedD0Pion;
  AliAODTrack* selectedD0Kaon;

  if(selectedDStarPion->Charge() == 1) selectedD0Pion = (AliAODTrack*)selectedD0->GetDaughter(0);
  if(selectedDStarPion->Charge() == -1) selectedD0Pion = (AliAODTrack*)selectedD0->GetDaughter(1);

  if(selectedDStarPion->Charge() == 1) selectedD0Kaon = (AliAODTrack*)selectedD0->GetDaughter(1);
  if(selectedDStarPion->Charge() == -1) selectedD0Kaon = (AliAODTrack*)selectedD0->GetDaughter(0);

  Double_t d0B0pion = TMath::Abs(selectedB0->Getd0Prong(0));
  Double_t d0DStarpion = TMath::Abs(selectedDStar->Getd0Prong(0));
  Double_t d0D0pion = 0;
  Double_t d0D0kaon = 0;

  if(selectedDStarPion->Charge() == 1) d0D0pion = selectedD0->Getd0Prong(0);
  if(selectedDStarPion->Charge() == -1) d0D0pion = selectedD0->Getd0Prong(1);

  if(selectedDStarPion->Charge() == 1) d0D0kaon = selectedD0->Getd0Prong(1);
  if(selectedDStarPion->Charge() == -1) d0D0kaon = selectedD0->Getd0Prong(0);

  Double_t pt_track = 0;
  Double_t momentum_track = 0;
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
  histType = 4;
  if(!isDesiredCandidate)
  {
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(selectedD0Pion->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }

    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0D0pion);
    ((TH2F*)fDaughterHistogramArray2D[daughterType][4])->Fill(ptB0,pt_track);
  }

  if(isDesiredCandidate)
  {
    histType = 5;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(selectedD0Pion->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }

    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0D0pion);
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
      ((TH1F*)fDaughterHistogramArrayExtra[0][2])->Fill(pdgCodeParticle);
      mcLabelParticleMother = mcTrackParticle->GetMother();

      if(mcLabelParticleMother >= 0){
        AliAODMCParticle *mcTrackParticleMother = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticleMother);
        pdgCodeParticleMother = TMath::Abs(mcTrackParticleMother->GetPdgCode());
       ((TH1F*)fDaughterHistogramArrayExtra[0][3])->Fill(pdgCodeParticleMother);
      }
    }
  }



  //fill the D0 kaon info
  pt_track = selectedD0Kaon->Pt();
  momentum_track = selectedD0Kaon->P();
  numberOfITS = selectedD0Kaon->GetITSNcls(); 
  numberOfTPC = selectedD0Kaon->GetTPCNcls();
  totalNumberOfITS += numberOfITS;
  totalNumberOfTPC += numberOfTPC;
  TPCok = trackPIDHF->GetnSigmaTPC(selectedD0Kaon, kaonPIDnumber, nSigmaTPC);
  TOFok = trackPIDHF->GetnSigmaTOF(selectedD0Kaon, kaonPIDnumber, nSigmaTOF);
  if(TPCok != -1) nSigmaTPCtotal += nSigmaTPC*nSigmaTPC;
  if(TOFok != -1) nSigmaTOFtotal += nSigmaTOF*nSigmaTOF;

  daughterType = 1;
  histType = 4;
  if(!isDesiredCandidate)
  {
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(selectedD0Kaon->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }

    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0D0kaon);
    ((TH2F*)fDaughterHistogramArray2D[daughterType][4])->Fill(ptB0,pt_track);
  }

  if(isDesiredCandidate)
  {
    histType = 5;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(selectedD0Kaon->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }

    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0D0kaon);
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
      ((TH1F*)fDaughterHistogramArrayExtra[1][2])->Fill(pdgCodeParticle);
      mcLabelParticleMother = mcTrackParticle->GetMother();

      if(mcLabelParticleMother >= 0){
        AliAODMCParticle *mcTrackParticleMother = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticleMother);
        pdgCodeParticleMother = TMath::Abs(mcTrackParticleMother->GetPdgCode());
       ((TH1F*)fDaughterHistogramArrayExtra[1][3])->Fill(pdgCodeParticleMother);
      }
    }
  }

  //fill the DStar pion info
  pt_track = selectedDStarPion->Pt();
  momentum_track = selectedDStarPion->P();
  numberOfITS = selectedDStarPion->GetITSNcls(); 
  numberOfTPC = selectedDStarPion->GetTPCNcls();
  totalNumberOfITS += numberOfITS;
  totalNumberOfTPC += numberOfTPC;
  TPCok = trackPIDHF->GetnSigmaTPC(selectedDStarPion, pionPIDnumber, nSigmaTPC);
  TOFok = trackPIDHF->GetnSigmaTOF(selectedDStarPion, pionPIDnumber, nSigmaTOF);
  if(TPCok != -1) nSigmaTPCtotal += nSigmaTPC*nSigmaTPC;

  daughterType = 2;
  histType = 4;
  if(!isDesiredCandidate)
  {
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(selectedDStarPion->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }

    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0DStarpion);
    ((TH2F*)fDaughterHistogramArray2D[daughterType][4])->Fill(ptB0,pt_track);
  }

  if(isDesiredCandidate)
  {
    histType = 5;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(selectedDStarPion->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }

    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0DStarpion);
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
      ((TH1F*)fDaughterHistogramArrayExtra[2][2])->Fill(pdgCodeParticle);
      mcLabelParticleMother = mcTrackParticle->GetMother();

      if(mcLabelParticleMother >= 0){
        AliAODMCParticle *mcTrackParticleMother = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticleMother);
        pdgCodeParticleMother = TMath::Abs(mcTrackParticleMother->GetPdgCode());
       ((TH1F*)fDaughterHistogramArrayExtra[2][3])->Fill(pdgCodeParticleMother);
      }
    }
  }

  //fill the B0 pion info
  pt_track = selectedB0Pion->Pt();
  momentum_track = selectedB0Pion->P();
  numberOfITS = selectedB0Pion->GetITSNcls(); 
  numberOfTPC = selectedB0Pion->GetTPCNcls();
  totalNumberOfITS += numberOfITS;
  totalNumberOfTPC += numberOfTPC;
  TPCok = trackPIDHF->GetnSigmaTPC(selectedB0Pion, pionPIDnumber, nSigmaTPC);
  TOFok = trackPIDHF->GetnSigmaTOF(selectedB0Pion, pionPIDnumber, nSigmaTOF);
  if(TPCok != -1) nSigmaTPCtotal += nSigmaTPC*nSigmaTPC;
  if(TOFok != -1) nSigmaTOFtotal += nSigmaTOF*nSigmaTOF;

  daughterType = 3;
  histType = 4;
  if(!isDesiredCandidate)
  {
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(selectedB0Pion->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }

    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0B0pion);
    ((TH2F*)fDaughterHistogramArray2D[daughterType][4])->Fill(ptB0,pt_track);
  }

  if(isDesiredCandidate)
  {
    histType = 5;
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][0])->Fill(pt_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][1])->Fill(momentum_track);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][2])->Fill(numberOfITS);
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][3])->Fill(numberOfTPC);

    for (Int_t j = 0; j < 10; ++j)
    {
      if(selectedB0Pion->HasPointOnITSLayer(j)) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(j);
        
    }

    if(TPCok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][4])->Fill(nSigmaTPC);
    if(TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][5])->Fill(nSigmaTOF);
    if(TPCok != -1 && TOFok != -1) ((TH1F*)fDaughterHistogramArray[daughterType][histType][6])->Fill(sqrt(nSigmaTPC*nSigmaTPC + nSigmaTOF*nSigmaTOF));
    ((TH1F*)fDaughterHistogramArray[daughterType][histType][8])->Fill(d0B0pion);
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
      ((TH1F*)fDaughterHistogramArrayExtra[3][2])->Fill(pdgCodeParticle);
      mcLabelParticleMother = mcTrackParticle->GetMother();

      if(mcLabelParticleMother >= 0){
        AliAODMCParticle *mcTrackParticleMother = (AliAODMCParticle*)mcTrackArray->At(mcLabelParticleMother);
        pdgCodeParticleMother = TMath::Abs(mcTrackParticleMother->GetPdgCode());
       ((TH1F*)fDaughterHistogramArrayExtra[3][3])->Fill(pdgCodeParticleMother);
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
Double_t AliAnalysisTaskSEB0toDStarPi::DeltaInvMassDStarKpipi(AliAODRecoDecayHF2Prong * DStar) const 
{
  ///
  /// 3 prong invariant mass of the D0 daughters and the soft pion
  ///
  Int_t chargeDStar = DStar->Charge();

  Double_t e[3];
  UInt_t prongs[2];
  if(chargeDStar==1){
    e[0]=((AliAODRecoDecayHF2Prong*)DStar->GetDaughter(1))->EProng(0,211);
    e[1]=((AliAODRecoDecayHF2Prong*)DStar->GetDaughter(1))->EProng(1,321);
    prongs[0] = 211;
    prongs[1] = 321;
  } 
  else if (chargeDStar==-1)
  {
    e[0]=((AliAODRecoDecayHF2Prong*)DStar->GetDaughter(1))->EProng(1,211);
    e[1]=((AliAODRecoDecayHF2Prong*)DStar->GetDaughter(1))->EProng(0,321);
    prongs[1] = 211;
    prongs[0] = 321;
  } 
  else 
  {
    std::cout << "Wrong charge DStar." << std::endl;
    return 0;
  }
  e[2]=DStar->EProng(0,211);

  Double_t esum = e[0]+e[1]+e[2];
  Double_t invMassDStar = TMath::Sqrt(esum*esum-DStar->P2());
  Double_t invMassD0 = ((AliAODRecoDecayHF2Prong*)DStar->GetDaughter(1))->InvMass(2,prongs);

  return invMassDStar - invMassD0; 
}
//-------------------------------------------------------------------------------------
Double_t AliAnalysisTaskSEB0toDStarPi::DeltaInvMassB0Kpipipi(AliAODRecoDecayHF2Prong *Bzero) const 
{
  ///
  /// 4 prong invariant mass of the D0 daughters, the soft pion, and the B0 pion
  ///

  AliAODRecoDecayHF2Prong * DStar = (AliAODRecoDecayHF2Prong*)Bzero->GetDaughter(1);
  Int_t chargeDStar = DStar->Charge();

  Double_t e[4];
  UInt_t prongs[2];
  if(chargeDStar==1){
    e[1]=((AliAODRecoDecayHF2Prong*)DStar->GetDaughter(1))->EProng(0,211);
    e[2]=((AliAODRecoDecayHF2Prong*)DStar->GetDaughter(1))->EProng(1,321);
    prongs[0] = 211;
    prongs[1] = 321;
  } 
  else if (chargeDStar==-1)
  {
    e[1]=((AliAODRecoDecayHF2Prong*)DStar->GetDaughter(1))->EProng(1,211);
    e[2]=((AliAODRecoDecayHF2Prong*)DStar->GetDaughter(1))->EProng(0,321);
    prongs[1] = 211;
    prongs[0] = 321;
  } 
  else 
  {
    std::cout << "Wrong charge DStar." << std::endl;
    return 0;
  }
  e[0]=DStar->EProng(0,211);
  e[3]=Bzero->EProng(0,211);

  Double_t esum = e[0]+e[1]+e[2]+e[3];
  Double_t invMassB0 = TMath::Sqrt(esum*esum-Bzero->P2());

  Double_t invMassD0 = ((AliAODRecoDecayHF2Prong*)DStar->GetDaughter(1))->InvMass(2,prongs);

  return invMassB0 - invMassD0; 
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::FillD0Histograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType, Int_t pdgCodeMother){

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
  UInt_t prongs[2], prongs2[2];
  Double_t angleBetweenBothDaughters = 0;
  Double_t cosThetaStar = 0;
  Double_t normDecayLength = 0;
  Double_t pdgMassMother=TDatabasePDG::Instance()->GetParticle(421)->Mass();


  prongs[0] = 211; prongs[1] = 321;
  prongs2[0] = 321; prongs2[1] = 211;
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
  if(pdgCodeMother == -421) invariantMassMother = selectedMother->InvMass(2,prongs2);

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

  Double_t vertexMotherX = vertexMother->GetX();
  Double_t vertexMotherY = vertexMother->GetY();
  Double_t vertexMotherZ = vertexMother->GetZ();

  Double_t cosPointingAngleXY = selectedMother->CosPointingAngleXY();
  Double_t distanceXYToVertex = vertexMother->DistanceXYToVertex(primaryVertex);
  Double_t normalizedDecayLengthXY = selectedMother->NormalizedDecayLengthXY();

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

  ((TH1F*)fMotherHistogramArray[motherType][histType][23])->Fill(vertexMotherX);
  ((TH1F*)fMotherHistogramArray[motherType][histType][24])->Fill(vertexMotherY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][25])->Fill(vertexMotherZ);

  ((TH1F*)fMotherHistogramArray[motherType][histType][36])->Fill(TMath::Abs(dd0pr1));
  ((TH1F*)fMotherHistogramArray[motherType][histType][37])->Fill(TMath::Abs(dd0pr2));
  ((TH1F*)fMotherHistogramArray[motherType][histType][38])->Fill(TMath::Abs(dd0max));
  ((TH1F*)fMotherHistogramArray[motherType][histType][39])->Fill(TMath::Abs(dd0min));

  ((TH1F*)fMotherHistogramArray[motherType][histType][40])->Fill(cosPointingAngleXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][41])->Fill(distanceXYToVertex);
  ((TH1F*)fMotherHistogramArray[motherType][histType][42])->Fill(normalizedDecayLengthXY);

  //we fill the 2D histograms
  Int_t nFirst = 0;
  Int_t nSecond = 1;
  Int_t nVariables = 10;
  Int_t nHistograms = nVariables * (nVariables - 1) / 2;
  for (Int_t k = 0; k < nHistograms; ++k)
  {
    Double_t firstVariable = 0.0;
    Double_t secondVariable = 0.0;

    if(nFirst==0) firstVariable = d0firstTrack;
    if(nFirst==1) firstVariable = d0secondTrack;
    if(nFirst==2) firstVariable = d0Mother;
    if(nFirst==3) firstVariable = pointingAngle;
    if(nFirst==4) firstVariable = impactProduct;
    if(nFirst==5) firstVariable = impactProductXY;
    if(nFirst==6) firstVariable = vertexDistance;
    if(nFirst==7) firstVariable = normDecayLength;
    if(nFirst==8) firstVariable = cosPointingAngleXY;
    if(nFirst==9) firstVariable = distanceXYToVertex;
    if(nFirst==10) firstVariable = normalizedDecayLengthXY;

    if(nSecond==0) secondVariable = d0firstTrack;
    if(nSecond==1) secondVariable = d0secondTrack;
    if(nSecond==2) secondVariable = d0Mother;
    if(nSecond==3) secondVariable = pointingAngle;
    if(nSecond==4) secondVariable = impactProduct;
    if(nSecond==5) secondVariable = impactProductXY;
    if(nSecond==6) secondVariable = vertexDistance;
    if(nSecond==7) secondVariable = normDecayLength;
    if(nSecond==8) secondVariable = cosPointingAngleXY;
    if(nSecond==9) secondVariable = distanceXYToVertex;
    if(nSecond==10) secondVariable = normalizedDecayLengthXY;

    ((TH2F*)fMotherHistogramArray2D[motherType][histType][k])->Fill(firstVariable,secondVariable);

    nSecond++;
    if(nSecond>nVariables)
    {
      nFirst++;
      nSecond = nFirst + 1;
    }
  }

  if(fUse3DHistograms)
  {
    //we fill the 3D histograms
    nFirst = 0;
    nSecond = 1;
    Int_t nThird = 2;
    nVariables = 10;
    nHistograms = nVariables * (nVariables - 1) * (nVariables - 2) / 6;
    for (Int_t k = 0; k < nHistograms; ++k)
    {
      Double_t firstVariable = 0.0;
      Double_t secondVariable = 0.0;
      Double_t thirdVariable = 0.0;

      if(nFirst==0) firstVariable = d0firstTrack;
      if(nFirst==1) firstVariable = d0secondTrack;
      if(nFirst==2) firstVariable = d0Mother;
      if(nFirst==3) firstVariable = pointingAngle;
      if(nFirst==4) firstVariable = impactProduct;
      if(nFirst==5) firstVariable = impactProductXY;
      if(nFirst==6) firstVariable = vertexDistance;
      if(nFirst==7) firstVariable = normDecayLength;
      if(nFirst==8) firstVariable = cosPointingAngleXY;
      if(nFirst==9) firstVariable = distanceXYToVertex;
      if(nFirst==10) firstVariable = normalizedDecayLengthXY;

      if(nSecond==0) secondVariable = d0firstTrack;
      if(nSecond==1) secondVariable = d0secondTrack;
      if(nSecond==2) secondVariable = d0Mother;
      if(nSecond==3) secondVariable = pointingAngle;
      if(nSecond==4) secondVariable = impactProduct;
      if(nSecond==5) secondVariable = impactProductXY;
      if(nSecond==6) secondVariable = vertexDistance;
      if(nSecond==7) secondVariable = normDecayLength;
      if(nSecond==8) secondVariable = cosPointingAngleXY;
      if(nSecond==9) secondVariable = distanceXYToVertex;
      if(nSecond==10) secondVariable = normalizedDecayLengthXY;

      if(nThird==0) thirdVariable = d0firstTrack;
      if(nThird==1) thirdVariable = d0secondTrack;
      if(nThird==2) thirdVariable = d0Mother;
      if(nThird==3) thirdVariable = pointingAngle;
      if(nThird==4) thirdVariable = impactProduct;
      if(nThird==5) thirdVariable = impactProductXY;
      if(nThird==6) thirdVariable = vertexDistance;
      if(nThird==7) thirdVariable = normDecayLength;
      if(nThird==8) thirdVariable = cosPointingAngleXY;
      if(nThird==9) thirdVariable = distanceXYToVertex;
      if(nThird==10) thirdVariable = normalizedDecayLengthXY;

      ((TH3F*)fMotherHistogramArray3D[motherType][histType][k])->Fill(firstVariable,secondVariable,thirdVariable);

      nThird++;
      if(nThird>nVariables)
      {
        nSecond++;
        nThird = nSecond + 1;
        if(nSecond>nVariables)
        {
          nFirst++;
          nSecond = nFirst + 1;
          nThird = nFirst + 2;

        }
      }
    }
  }

  return;
}
//-------------------------------------------------------------------------------------
void AliAnalysisTaskSEB0toDStarPi::FillDStarAndB0Histograms(AliAODRecoDecayHF2Prong * selectedMother, AliAODVertex *primaryVertex, Double_t bz, Int_t motherType, Int_t histType){

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
    AliAODRecoDecayHF2Prong * secondDaughter = (AliAODRecoDecayHF2Prong*)selectedMother->GetDaughter(1);

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

  Double_t vertexMotherX = vertexMother->GetX();
  Double_t vertexMotherY = vertexMother->GetY();
  Double_t vertexMotherZ = vertexMother->GetZ();

  Double_t cosPointingAngleXY = selectedMother->CosPointingAngleXY();
  Double_t distanceXYToVertex = vertexMother->DistanceXYToVertex(primaryVertex);
  Double_t normalizedDecayLengthXY = selectedMother->NormalizedDecayLengthXY();

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

  ((TH1F*)fMotherHistogramArray[motherType][histType][23])->Fill(vertexMotherX);
  ((TH1F*)fMotherHistogramArray[motherType][histType][24])->Fill(vertexMotherY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][25])->Fill(vertexMotherZ);

  ((TH1F*)fMotherHistogramArray[motherType][histType][36])->Fill(TMath::Abs(dd0pr1));
  ((TH1F*)fMotherHistogramArray[motherType][histType][37])->Fill(TMath::Abs(dd0pr2));
  ((TH1F*)fMotherHistogramArray[motherType][histType][38])->Fill(TMath::Abs(dd0max));
  ((TH1F*)fMotherHistogramArray[motherType][histType][39])->Fill(TMath::Abs(dd0min));

  ((TH1F*)fMotherHistogramArray[motherType][histType][40])->Fill(cosPointingAngleXY);
  ((TH1F*)fMotherHistogramArray[motherType][histType][41])->Fill(distanceXYToVertex);
  ((TH1F*)fMotherHistogramArray[motherType][histType][42])->Fill(normalizedDecayLengthXY);

  //we fill the 2D histograms
  Int_t nFirst = 0;
  Int_t nSecond = 1;
  Int_t nVariables = 10;
  Int_t nHistograms = nVariables * (nVariables - 1) / 2;
  for (Int_t k = 0; k < nHistograms; ++k)
  {
    Double_t firstVariable = 0.0;
    Double_t secondVariable = 0.0;

    if(nFirst==0) firstVariable = d0firstTrack;
    if(nFirst==1) firstVariable = d0secondTrack;
    if(nFirst==2) firstVariable = d0Mother;
    if(nFirst==3) firstVariable = pointingAngle;
    if(nFirst==4) firstVariable = impactProduct;
    if(nFirst==5) firstVariable = impactProductXY;
    if(nFirst==6) firstVariable = vertexDistance;
    if(nFirst==7) firstVariable = normDecayLength;
    if(nFirst==8) firstVariable = cosPointingAngleXY;
    if(nFirst==9) firstVariable = distanceXYToVertex;
    if(nFirst==10) firstVariable = normalizedDecayLengthXY;

    if(nSecond==0) secondVariable = d0firstTrack;
    if(nSecond==1) secondVariable = d0secondTrack;
    if(nSecond==2) secondVariable = d0Mother;
    if(nSecond==3) secondVariable = pointingAngle;
    if(nSecond==4) secondVariable = impactProduct;
    if(nSecond==5) secondVariable = impactProductXY;
    if(nSecond==6) secondVariable = vertexDistance;
    if(nSecond==7) secondVariable = normDecayLength;
    if(nSecond==8) secondVariable = cosPointingAngleXY;
    if(nSecond==9) secondVariable = distanceXYToVertex;
    if(nSecond==10) secondVariable = normalizedDecayLengthXY;

    ((TH2F*)fMotherHistogramArray2D[motherType][histType][k])->Fill(firstVariable,secondVariable);

    nSecond++;
    if(nSecond>nVariables)
    {
      nFirst++;
      nSecond = nFirst + 1;
    }
  }

  if(fUse3DHistograms)
  {
    //we fill the 3D histograms
    nFirst = 0;
    nSecond = 1;
    Int_t nThird = 2;
    nVariables = 10;
    nHistograms = nVariables * (nVariables - 1) * (nVariables - 2) / 6;
    for (Int_t k = 0; k < nHistograms; ++k)
    {
      Double_t firstVariable = 0.0;
      Double_t secondVariable = 0.0;
      Double_t thirdVariable = 0.0;

      if(nFirst==0) firstVariable = d0firstTrack;
      if(nFirst==1) firstVariable = d0secondTrack;
      if(nFirst==2) firstVariable = d0Mother;
      if(nFirst==3) firstVariable = pointingAngle;
      if(nFirst==4) firstVariable = impactProduct;
      if(nFirst==5) firstVariable = impactProductXY;
      if(nFirst==6) firstVariable = vertexDistance;
      if(nFirst==7) firstVariable = normDecayLength;
      if(nFirst==8) firstVariable = cosPointingAngleXY;
      if(nFirst==9) firstVariable = distanceXYToVertex;
      if(nFirst==10) firstVariable = normalizedDecayLengthXY;

      if(nSecond==0) secondVariable = d0firstTrack;
      if(nSecond==1) secondVariable = d0secondTrack;
      if(nSecond==2) secondVariable = d0Mother;
      if(nSecond==3) secondVariable = pointingAngle;
      if(nSecond==4) secondVariable = impactProduct;
      if(nSecond==5) secondVariable = impactProductXY;
      if(nSecond==6) secondVariable = vertexDistance;
      if(nSecond==7) secondVariable = normDecayLength;
      if(nSecond==8) secondVariable = cosPointingAngleXY;
      if(nSecond==9) secondVariable = distanceXYToVertex;
      if(nSecond==10) secondVariable = normalizedDecayLengthXY;

      if(nThird==0) thirdVariable = d0firstTrack;
      if(nThird==1) thirdVariable = d0secondTrack;
      if(nThird==2) thirdVariable = d0Mother;
      if(nThird==3) thirdVariable = pointingAngle;
      if(nThird==4) thirdVariable = impactProduct;
      if(nThird==5) thirdVariable = impactProductXY;
      if(nThird==6) thirdVariable = vertexDistance;
      if(nThird==7) thirdVariable = normDecayLength;
      if(nThird==8) thirdVariable = cosPointingAngleXY;
      if(nThird==9) thirdVariable = distanceXYToVertex;
      if(nThird==10) thirdVariable = normalizedDecayLengthXY;

      ((TH3F*)fMotherHistogramArray3D[motherType][histType][k])->Fill(firstVariable,secondVariable,thirdVariable);

      nThird++;
      if(nThird>nVariables)
      {
        nSecond++;
        nThird = nSecond + 1;
        if(nSecond>nVariables)
        {
          nFirst++;
          nSecond = nFirst + 1;
          nThird = nFirst + 2;

        }
      }
    }
  }

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

    ((TH1F*)fMotherHistogramArray[motherType][histType][26])->Fill(pointingAngleToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][27])->Fill(d0D0);
    ((TH1F*)fMotherHistogramArray[motherType][histType][28])->Fill(d0FirstDaughterToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][29])->Fill(d0SecondDaughterToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][30])->Fill(impactProductToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][31])->Fill(impactProductXYToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][32])->Fill(normDecayLengthToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][33])->Fill(pseudoProperDecayTimeToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][34])->Fill(DecayTimeToDStar);
    ((TH1F*)fMotherHistogramArray[motherType][histType][35])->Fill(normDecayTimeToDStar);


  }
  return;
}
//-------------------------------------------------------------------------------------
Int_t AliAnalysisTaskSEB0toDStarPi::MatchCandidateToMonteCarlo(Int_t pdgabs, AliAODRecoDecayHF2Prong * candidate, TClonesArray *mcArray, TMatrix * B0toDStarPiLabelMatrix) const
{
  //
  // Check if this candidate is matched to a MC signal
  // If no, return -1
  // If yes, return label (>=0) of the AliAODMCParticle


  // Check number of daughters
  Int_t ndg = candidate->GetNDaughters();
  if(!ndg) { AliError("No daughters available"); return -1;}
  if(ndg != 2) return -1;

  // loop on daughters and write the labels
  Int_t dgLabels[2] = {-1};
  Int_t pdgDg[2] = {0};
  Int_t signalPosition = -1;
  if(pdgabs==421)
  {
    AliAODTrack *trk0 = (AliAODTrack*)candidate->GetDaughter(0);
    dgLabels[0] = trk0->GetLabel();
    AliAODTrack *trk1 = (AliAODTrack*)candidate->GetDaughter(1);
    dgLabels[1] = trk1->GetLabel();
    pdgDg[0] = 211; pdgDg[1] = 321;
    signalPosition = 4;
  }
  else if(pdgabs==413)
  {
    AliAODTrack *trk0 = (AliAODTrack*)candidate->GetDaughter(0);
    dgLabels[0] = trk0->GetLabel();
    dgLabels[1] = MatchCandidateToMonteCarlo(421,(AliAODRecoDecayHF2Prong*)candidate->GetDaughter(1), mcArray, B0toDStarPiLabelMatrix);
    pdgDg[0] = 211; pdgDg[1] = 421;
    signalPosition = 5;
  }
  else if(pdgabs==511)
  {
    AliAODTrack *trk0 = (AliAODTrack*)candidate->GetDaughter(0);
    dgLabels[0] = trk0->GetLabel();
    dgLabels[1] = MatchCandidateToMonteCarlo(413,(AliAODRecoDecayHF2Prong*)candidate->GetDaughter(1), mcArray, B0toDStarPiLabelMatrix);
    pdgDg[0] = 211; pdgDg[1] = 413;
    signalPosition = 6;
  }
  else
  { 
    std::cout << "Wrong pdg supplied for function to match candidate to monte carlo signal." << std::endl;
    return -1;
  }
  if(dgLabels[0]==-1) return -1;
  if(dgLabels[1]==-1) return -1;


  Int_t labMom[2]={0,0};
  Int_t i,j,lab,labMother,pdgMother,pdgPart;
  AliAODMCParticle *part=0;
  AliAODMCParticle *mother=0;
  Double_t pxSumDgs=0.,pySumDgs=0.,pzSumDgs=0.;
  Bool_t pdgUsed[2]={kFALSE,kFALSE};

  // loop on daughter labels
  for(i=0; i<ndg; i++) 
  {
    labMom[i]=-1;
    lab = TMath::Abs(dgLabels[i]);
    if(lab<0) 
    {
      printf("daughter with negative label %d\n",lab);
      return -1;
    }
    part = (AliAODMCParticle*)mcArray->At(lab);
    if(!part) 
    { 
      printf("no MC particle\n");
      return -1;
    }

    // check the PDG of the daughter
    pdgPart=TMath::Abs(part->GetPdgCode());
    for(j=0; j<ndg; j++) 
    {
      if(!pdgUsed[j] && pdgPart==pdgDg[j]) 
      {
        pdgUsed[j]=kTRUE;
        break;
      }
    }
    

    mother = part;
    while(mother->GetMother()>=0) 
    {
      labMother=mother->GetMother();
      mother = (AliAODMCParticle*)mcArray->At(labMother);
      if(!mother) 
      {
        printf("no MC mother particle\n");
        break;
      }
      pdgMother = TMath::Abs(mother->GetPdgCode());
      if(pdgMother==pdgabs) 
      {
        labMom[i]=labMother;
        // keep sum of daughters' momenta, to check for mom conservation
        pxSumDgs += part->Px();
        pySumDgs += part->Py();
        pzSumDgs += part->Pz();
        break;
      } 
      else break;
    }
    if(labMom[i]==-1) return -1; // mother PDG not ok for this daughter
  } // end loop on daughters

  // check if the candidate is signal
  labMother=labMom[0];
  // all labels have to be the same and !=-1
  for(i=0; i<ndg; i++) 
  {
    if(labMom[i]==-1)        return -1;
    if(labMom[i]!=labMother) return -1;
  }

  // check that all daughter PDGs are matched
  for(i=0; i<ndg; i++) 
  {
    if(pdgUsed[i]==kFALSE) return -1;
  }

  // Check for mom conservation
  mother = (AliAODMCParticle*)mcArray->At(labMother);
  Double_t pxMother = mother->Px();
  Double_t pyMother = mother->Py();
  Double_t pzMother = mother->Pz();
  // within 0.5%
  if((TMath::Abs(pxMother-pxSumDgs)/(TMath::Abs(pxMother)+1.e-13)) > 0.05 &&
     (TMath::Abs(pyMother-pySumDgs)/(TMath::Abs(pyMother)+1.e-13)) > 0.05 &&
     (TMath::Abs(pzMother-pzSumDgs)/(TMath::Abs(pzMother)+1.e-13)) > 0.05) 
  { 
    return -1;
  }

  // Check if label matches a signal label
  // Int_t bIsSignal = kFALSE;
  // TMatrix &particleMatrix = *B0toDStarPiLabelMatrix;
  // for (Int_t k = 0; k < B0toDStarPiLabelMatrix->GetNrows(); ++k)
  // {
  //   if(labMother == (Int_t)particleMatrix(k,signalPosition))
  //   {
  //     bIsSignal = kTRUE;
  //     break;
  //   }
  // }
  // if(!bIsSignal) return -1;

  return labMother;
}
//-------------------------------------------------------------------------------------
Int_t AliAnalysisTaskSEB0toDStarPi::IsTrackInjected(AliAODTrack *part,AliAODMCHeader *header,TClonesArray *arrayMC){

  AliVertexingHFUtils* ggg = new  AliVertexingHFUtils();

  Int_t lab=part->GetLabel();
  if(lab<0) {delete ggg; ggg = nullptr; return 1;} //
  TString nameGen = ggg->GetGenerator(lab,header);
  TString empty="";
  Int_t countControl =0;
  while(nameGen.IsWhitespace()){
    AliAODMCParticle *mcpart= (AliAODMCParticle*)arrayMC->At(lab);
    if(!mcpart){
            printf("AliVertexingHFUtils::IsTrackInjected - BREAK: No valid AliAODMCParticle at label %i\n",lab);
            break;
    }
    Int_t mother = mcpart->GetMother();
    if(mother<0){
            // printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Reached primary particle without valid mother\n");
            break;
    }
    lab=mother;
    nameGen = ggg->GetGenerator(mother,header);
    countControl++;
    if(countControl>=10){ // 10 = arbitrary number; protection from infinite loops
            printf("AliVertexingHFUtils::IsTrackInjected - BREAK: Protection from infinite loop active\n");
            break;
    }
  }
  if(nameGen.IsWhitespace() || nameGen.Contains("ijing")){delete ggg; ggg = nullptr; return 0;}

  delete ggg; ggg = nullptr;
  return 1;
}
//-------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskSEB0toDStarPi::IsCandidateInjected(AliAODRecoDecayHF2Prong *selectedB0, AliAODMCHeader *header,TClonesArray *arrayMC){

  AliAODTrack* selectedB0Pion = (AliAODTrack*)selectedB0->GetDaughter(0);
  AliAODRecoDecayHF2Prong* selectedDStar = (AliAODRecoDecayHF2Prong*)selectedB0->GetDaughter(1);

  AliAODTrack* selectedDStarPion = (AliAODTrack*)selectedDStar->GetDaughter(0);
  AliAODRecoDecayHF2Prong* selectedD0 = (AliAODRecoDecayHF2Prong*)selectedDStar->GetDaughter(1);

  AliAODTrack* selectedD0FirstDaughter = (AliAODTrack*)selectedD0->GetDaughter(0);
  AliAODTrack* selectedD0SecondDaughter = (AliAODTrack*)selectedD0->GetDaughter(1);

  if(IsTrackInjected(selectedB0Pion,header,arrayMC)) return kTRUE;
  if(IsTrackInjected(selectedDStarPion,header,arrayMC)) return kTRUE;
  if(IsTrackInjected(selectedD0FirstDaughter,header,arrayMC)) return kTRUE;
  if(IsTrackInjected(selectedD0SecondDaughter,header,arrayMC)) return kTRUE;

  return kFALSE;
}
