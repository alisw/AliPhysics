#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TProfile.h"

#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>
#include "AliAODMCParticle.h"
#include <TNtuple.h>
#include <AliInputEventHandler.h>
#include <AliAODTrack.h>
#include "AliPIDResponse.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliTrackReference.h"
#include "AliESDtrackCuts.h"

#include "AliAnalysisH3MC.h"

class AliAnalysisH3MC;
using namespace std;

//const char *GetID(int pdgCode);
ClassImp(AliAnalysisH3MC)

    //________________________________________________________________________
    AliAnalysisH3MC::AliAnalysisH3MC()
    : AliAnalysisTaskSE(), fVEvent(0), fAODEvent(0), fOutputList(0), fOutputEvent(0), fOutputProtons(0), fOutputAProtons(0), fOutputDeuterons(0), fOutputADeuterons(0), fOutputH3(0), fOutputAH3(0), fMCOutputList(0), fMCOutputEvent(0), fMCOutputProtons(0), fMCOutputAProtons(0), fMCOutputDeuterons(0), fMCOutputADeuterons(0), fMCOutputH3(0), fMCOutputAH3(0), fHistTrackCuts(0), fEventStat(0), fVertexZ(0), fVtxContrib(0), fSPDVtxResol(0), fVtxDisplacement(0), fMultV0(0), fCentV0M(0), fHistsProton(), fHistsAProton(), fHistsDeuteron(), fHistsADeuteron(), fPIDResponse(0), fFilterBit(256), fLowPCut(0.0), fHighPCut(1e30), fEtaCut(1.0), fMinClIts(0), fMaxDCAxyCut(10.), fMaxDCAxyAna(10.), fMaxDCAz(10.), fMaxTPCnSigma(10.), fMaxTOFnSigma(10.), fUseTOFPidCut(kFALSE), fMomTOFProt(10.), fMomTOFDeut(10.), fTrackCuts(0)
{
  // default constructor, don't allocate memory here!
  // this is used by root for IO purposes, it needs to remain empty
}

//________________________________________________________________________
AliAnalysisH3MC::AliAnalysisH3MC(const char *name)
    : AliAnalysisTaskSE(name), fVEvent(0), mcEvent(0), fAODEvent(0), fOutputList(0), fOutputEvent(0), fOutputProtons(0), fOutputAProtons(0), fOutputDeuterons(0), fOutputADeuterons(0), fHistTrackCuts(0), fEventStat(0), fVertexZ(0), fVtxContrib(0), fSPDVtxResol(0), fVtxDisplacement(0), fMultV0(0), fCentV0M(0), fHistsProton(), fHistsAProton(), fHistsDeuteron(), fHistsADeuteron(), fPIDResponse(0), fFilterBit(256), fLowPCut(0.0), fHighPCut(1e30), fEtaCut(1.0), fMinClIts(0), fMaxDCAxyCut(10.), fMaxDCAxyAna(10.), fMaxDCAz(10.), fMaxTPCnSigma(10.), fMaxTOFnSigma(10.), fUseTOFPidCut(kFALSE), fMomTOFProt(10.), fMomTOFDeut(10.), fTrackCuts(0)
{
  // constructor

  fHistsProton.clear();
  fHistsAProton.clear();
  fHistsDeuteron.clear();
  fHistsADeuteron.clear();
  fHistsH3.clear();
  fHistsAH3.clear();

  DefineInput(0, TChain::Class()); // a 'chain' of events created by the analysis manager automatically
  DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this case it's a list of histograms
  // one can add more output objects by calling DefineOutput(2, classname::Class())
  // make sure to connect them properly in AddTask and to call PostData() for all of them!
}

//________________________________________________________________________
AliAnalysisH3MC::~AliAnalysisH3MC()
{
  // destructor
  if (fOutputList)
  {
    delete fOutputList;
  }
}

//________________________________________________________________________
void AliAnalysisH3MC::UserCreateOutputObjects()
{
  // Called once at the start of the analysis

  // retrieve PID object from the input handler
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler = (AliInputEventHandler *)(man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  // create output lists
  fOutputList = new TList();
  fOutputList->SetName("Output");
  fOutputList->SetOwner(kTRUE); // the list is owner of all objects it contains and will delete them

  fOutputEvent = new TList();
  fOutputEvent->SetName("Event");
  fOutputEvent->SetOwner(kTRUE);

  fOutputProtons = new TList();
  fOutputProtons->SetName("Protons");
  fOutputProtons->SetOwner(kTRUE);

  fOutputAProtons = new TList();
  fOutputAProtons->SetName("AProtons");
  fOutputAProtons->SetOwner(kTRUE);

  fOutputDeuterons = new TList();
  fOutputDeuterons->SetName("Deuterons");
  fOutputDeuterons->SetOwner(kTRUE);

  fOutputADeuterons = new TList();
  fOutputADeuterons->SetName("ADeuterons");
  fOutputADeuterons->SetOwner(kTRUE);

  fOutputH3 = new TList();
  fOutputH3->SetName("H3");
  fOutputH3->SetOwner(kTRUE);

  fOutputAH3 = new TList();
  fOutputAH3->SetName("AH3");
  fOutputAH3->SetOwner(kTRUE);

  //Start of MC Outputs

  fMCOutputList = new TList();
  fMCOutputList->SetName("MCOutput");
  fMCOutputList->SetOwner(kTRUE);

  fMCOutputEvent = new TList();
  fMCOutputEvent->SetName("Event");
  fMCOutputEvent->SetOwner(kTRUE);

  fMCOutputProtons = new TList();
  fMCOutputProtons->SetName("Protons");
  fMCOutputProtons->SetOwner(kTRUE);

  fMCOutputAProtons = new TList();
  fMCOutputAProtons->SetName("AProtons");
  fMCOutputAProtons->SetOwner(kTRUE);

  fMCOutputDeuterons = new TList();
  fMCOutputDeuterons->SetName("Deuterons");
  fMCOutputDeuterons->SetOwner(kTRUE);

  fMCOutputADeuterons = new TList();
  fMCOutputADeuterons->SetName("ADeuterons");
  fMCOutputADeuterons->SetOwner(kTRUE);

  fMCOutputH3 = new TList();
  fMCOutputH3->SetName("H3");
  fMCOutputH3->SetOwner(kTRUE);

  fMCOutputAH3 = new TList();
  fMCOutputAH3->SetName("AH3");
  fMCOutputAH3->SetOwner(kTRUE);

  primaries = new TList();
  primaries->SetName("Primaries");
  primaries->SetOwner(kTRUE);

  secondariesW = new TList();
  secondariesW->SetName("SecondariesWeak");
  secondariesW->SetOwner(kTRUE);

  secondariesM = new TList();
  secondariesM->SetName("SecondariesMaterial");
  secondariesM->SetOwner(kTRUE);

  std::vector<TH1 *> hists;
  CreateHistosTrack(hists);
  CreateHistosMC(hists);

  for (int i = 0; i < hists.size(); i++)
  {
    primaries->Add((TH1 *)hists.at(i)->Clone());
    secondariesW->Add((TH1 *)hists.at(i)->Clone());
    secondariesM->Add((TH1 *)hists.at(i)->Clone());
  }

  fMCOutputADeuterons->Add(primaries->Clone());
  fMCOutputADeuterons->Add(secondariesW->Clone());
  fMCOutputADeuterons->Add(secondariesM->Clone());

  fMCOutputDeuterons->Add(primaries->Clone());
  fMCOutputDeuterons->Add(secondariesW->Clone());
  fMCOutputDeuterons->Add(secondariesM->Clone());

  fMCOutputProtons->Add(primaries->Clone());
  fMCOutputProtons->Add(secondariesW->Clone());
  fMCOutputProtons->Add(secondariesM->Clone());

  fMCOutputAProtons->Add(primaries->Clone());
  fMCOutputAProtons->Add(secondariesW->Clone());
  fMCOutputAProtons->Add(secondariesM->Clone());

  fMCOutputH3->Add(primaries->Clone());
  fMCOutputH3->Add(secondariesW->Clone());
  fMCOutputH3->Add(secondariesM->Clone());

  fMCOutputAH3->Add(primaries->Clone());
  fMCOutputAH3->Add(secondariesW->Clone());
  fMCOutputAH3->Add(secondariesM->Clone());

  fNtupleH3 = new TNtuple("fNtupleH3", "fNtupleH3", "p:pt:TPCSignal:TPCnSigmaH3:DCAxy:DCAz:TOFm2:TPCNClusters:ITSNClusters:TPCClusters4dEdx:Eta:ITSnSigmaH3:Chi2TPC:Chi2ITS:TPCCrossedRows:label");
  ((TList *)(fMCOutputH3->At(2)))->Add(fNtupleH3);

  //Create sublists and add them to the MC output lists.

  fOutputList->Add(fOutputEvent);
  fOutputList->Add(fOutputProtons);
  fOutputList->Add(fOutputAProtons);
  fOutputList->Add(fOutputDeuterons);
  fOutputList->Add(fOutputADeuterons);
  fOutputList->Add(fOutputH3);
  fOutputList->Add(fOutputAH3);

  fMCOutputList->Add(fMCOutputProtons);
  fMCOutputList->Add(fMCOutputAProtons);
  fMCOutputList->Add(fMCOutputDeuterons);
  fMCOutputList->Add(fMCOutputADeuterons);
  fMCOutputList->Add(fMCOutputH3);
  fMCOutputList->Add(fMCOutputAH3);
  fOutputList->Add(fMCOutputList);

  // create event histograms
  fEventStat = new TH1I("fEventStat", "Event statistics", kNbinsEvent, 0, kNbinsEvent);
  fEventStat->GetXaxis()->SetBinLabel(1, "After Phys. Sel. and trigger");
  fEventStat->GetXaxis()->SetBinLabel(2, "kINT7 selected (cross-check)");
  fEventStat->GetXaxis()->SetBinLabel(3, "DAQ imcomplete (cross-check)");
  fEventStat->GetXaxis()->SetBinLabel(4, "V0 timing (cross-check)");
  fEventStat->GetXaxis()->SetBinLabel(5, "Cluster-tracklet cut (cross-check)");
  fEventStat->GetXaxis()->SetBinLabel(6, "Vertex Z position");
  fEventStat->GetXaxis()->SetBinLabel(7, "N of vtx contrib");
  fEventStat->GetXaxis()->SetBinLabel(8, "SPD pile-up in mult. bins");
  fEventStat->GetXaxis()->SetBinLabel(9, "Vertex displacement");
  fEventStat->GetXaxis()->SetBinLabel(10, "Vertex Z resolution");
  fOutputEvent->Add(fEventStat);

  fVertexZ = new TH1F("fVertexZ", "Vertex Z distribution", 300, -15., 15.);
  fVertexZ->GetXaxis()->SetTitle("vertex Z, cm");
  fVertexZ->Sumw2();
  fOutputEvent->Add(fVertexZ);

  fVtxContrib = new TH1F("fVtxContrib", "N vertex contributors", 200, 0., 200.);
  fVtxContrib->GetXaxis()->SetTitle("N vertex contrib.");
  fOutputEvent->Add(fVtxContrib);

  fSPDVtxResol = new TH1F("fSPDVtxResol", "SPD vertex resolution", 100, 0., 1.);
  fSPDVtxResol->GetXaxis()->SetTitle("SPD vertex resolution, cm");
  fOutputEvent->Add(fSPDVtxResol);

  fVtxDisplacement = new TH2F("fVtxDisplacement", "SPD vertex displacement", 300, -15., 15., 300, -15., 15.);
  fVtxDisplacement->GetXaxis()->SetTitle("SPD vertex position, cm");
  fVtxDisplacement->GetYaxis()->SetTitle("Track vertex position, cm");
  fOutputEvent->Add(fVtxDisplacement);

  fMultV0 = new TH1F("fMultV0", "Mult V0 distrobution", 2000, 0., 2000.);
  fMultV0->GetXaxis()->SetTitle("V0 mult.");
  fOutputEvent->Add(fMultV0);

  fCentV0M = new TH1F("fCentV0M", "Centrality V0M", 300, -50., 250.);
  fCentV0M->GetXaxis()->SetTitle("V0M percentile");
  fOutputEvent->Add(fCentV0M);

  // track cuts config
  fHistTrackCuts = new TProfile("fHistTrackCuts", "TrackCuts config", 10, 0, 10);
  fHistTrackCuts->GetXaxis()->SetBinLabel(1, "Filter Bit");
  fHistTrackCuts->GetXaxis()->SetBinLabel(2, "low #it{p} cut");
  fHistTrackCuts->GetXaxis()->SetBinLabel(3, "#eta cut");
  fHistTrackCuts->GetXaxis()->SetBinLabel(4, "Min. N ITS cl");
  fHistTrackCuts->GetXaxis()->SetBinLabel(5, "Max. DCAxy cut");
  fHistTrackCuts->GetXaxis()->SetBinLabel(6, "Max. DCAxy ana");
  fHistTrackCuts->GetXaxis()->SetBinLabel(7, "Max. DCAz");
  fHistTrackCuts->GetXaxis()->SetBinLabel(8, "TPCnSigma cut");
  fHistTrackCuts->GetXaxis()->SetBinLabel(9, "TOFnSigma cut");
  fHistTrackCuts->GetXaxis()->SetBinLabel(10, "TOFmomentum prot");
  fHistTrackCuts->GetXaxis()->SetBinLabel(11, "TOFmomentum deut");
  fOutputList->Add(fHistTrackCuts);

  fHistTrackCuts->Fill(0.0, fFilterBit);
  fHistTrackCuts->Fill(1, fLowPCut);
  fHistTrackCuts->Fill(2, fEtaCut);
  fHistTrackCuts->Fill(3, fMinClIts);
  fHistTrackCuts->Fill(4, fMaxDCAxyCut);
  fHistTrackCuts->Fill(5, fMaxDCAxyAna);
  fHistTrackCuts->Fill(6, fMaxDCAz);
  fHistTrackCuts->Fill(7, fMaxTPCnSigma);
  fHistTrackCuts->Fill(8, fMaxTOFnSigma);
  fHistTrackCuts->Fill(9, fMomTOFProt);
  fHistTrackCuts->Fill(10, fMomTOFDeut);

  // (anti)proton histograms
  // CreateHistosTrack() creates a list of histograms which are added to the vector to be accessed later. Since the addition of more histograms changes the numberings, these have to be handled carefully.
  CreateHistosTrack(fHistsProton);
  for (Int_t i = 0; i < (int)(fHistsProton.size()); i++)
  {
    fOutputProtons->Add(fHistsProton.at(i));
  }

  CreateHistosTrack(fHistsAProton);
  for (Int_t i = 0; i < (int)(fHistsProton.size()); i++)
  {
    fOutputAProtons->Add(fHistsAProton.at(i));
  }

  // (anti)deuteron histograms
  CreateHistosTrack(fHistsDeuteron);
  for (Int_t i = 0; i < (int)(fHistsDeuteron.size()); i++)
  {
    fOutputDeuterons->Add(fHistsDeuteron.at(i));
  }

  CreateHistosTrack(fHistsADeuteron);
  for (Int_t i = 0; i < (int)(fHistsADeuteron.size()); i++)
  {
    fOutputADeuterons->Add(fHistsADeuteron.at(i));
  }

  // anti(triton) histograms 
  CreateHistosTrack(fHistsH3);
  for (Int_t i = 0; i < (int)(fHistsH3.size()); i++)
  {
    fOutputH3->Add(fHistsH3.at(i));
  }

  CreateHistosTrack(fHistsAH3);
  for (Int_t i = 0; i < (int)(fHistsAH3.size()); i++)
  {
    fOutputAH3->Add(fHistsAH3.at(i));
  }

  PostData(1, fOutputList);
  // postdata will notify the analysis manager of changes / updates to the
  // fOutputList object. the manager will in the end take care of writing output to file
  // so it needs to know what's in the output
}

//________________________________________________________________________
void AliAnalysisH3MC::UserExec(Option_t *)
{
  // Main loop
  // Called for each event
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();

  AliMCEventHandler *mcHandler = (AliMCEventHandler *)man->GetMCtruthEventHandler();
  if (!mcHandler)
  {
    cout << "No MC Hanlder!" << endl;
    return;
  }

  AliMCEvent *testMCEvent = mcHandler->MCEvent();

  if (!testMCEvent and kFALSE)
  {
    ::Error("AliAnalysisH3MC::UserExec", "MC handler not initialised.");
    PostData(1, fOutputList);
    return;
  } //else {cout <<"Handler event found"<<endl;}

  fVEvent = dynamic_cast<AliVEvent *>(InputEvent());
  if (!fVEvent)
  {
    ::Error("AliAnalysisH3MC::UserExec", "AliVEvent not found.");
    PostData(1, fOutputList);
    return;
  }

  fAODEvent = dynamic_cast<AliAODEvent *>(InputEvent());
  if (!fAODEvent)
  {
    ::Error("AliAnalysisH3MC::UserExec", "AliAODEvent not found.");
    PostData(1, fOutputList);
    return;
  }
  mcEvent = MCEvent();
  //mcEvent = dynamic_cast<AliMCEvent*>(InputEvent());

  if (!mcEvent)
  {
    ::Error("AliAnalysisH3MC::UserExec", "MCEvent not found.");
    PostData(1, fOutputList);
    return;
  }
  else
  {
    //cout<<"AliMCEvent found!"<<endl;
  }
  Int_t nMCtracks = mcEvent->GetNumberOfTracks();

  fEventStat->Fill(kSelectedEvents); // all events after trigger and physics selection
  // only for cross-check: is kINT7 selected
  UInt_t fSelectMask = ((AliInputEventHandler *)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isINT7selected = fSelectMask; // & AliVEvent::kINT7;
  /*if (!isINT7selected)
  {
    PostData(1, fOutputList);
    return;
  }*/
  fEventStat->Fill(kINT7selected);

  // only for cross-check: Incomplete events from DAQ
  //    if (fVEvent->IsIncompleteDAQ()) { // probably not implemented for AOD?
  //        PostData(1, fOutputList);
  //        return;
  //    }
  //    fEventStat->Fill(kDAQincomplete);
  //
  //    //only for cross-check: VZERO timing desicion
  //    Int_t fV0ADecision;
  //    Int_t fV0CDecision;
  //    AliVVZERO* vzeroData = fVEvent->GetVZEROData();
  //    fV0ADecision = vzeroData->GetV0ADecision();
  //    fV0CDecision = vzeroData->GetV0CDecision();
  //    if (!fV0ADecision || !fV0CDecision) {
  //        PostData(1, fOutputList);
  //        return;
  //    }
  //    fEventStat->Fill(kV0timing);
  //
  //    //only for cross-check: SPD clusters vs tracklets cut
  //    AliAnalysisUtils *utils = new AliAnalysisUtils();
  //    //utils->SetASPDCvsTCut(100.); // 65 by default
  //    //utils->SetBSPDCvsTCut(5.); // 4 by default
  //    Bool_t ClustersVsTrackletBG = utils->IsSPDClusterVsTrackletBG(fVEvent);
  //    if (ClustersVsTrackletBG) {
  //        PostData(1, fOutputList);
  //        return;
  //    }
  fEventStat->Fill(kClusterTrackletCut);

  // vertex filter
  const AliVVertex *vertex = fVEvent->GetPrimaryVertex();
  Bool_t lVertexAcceptable = (TMath::Abs(vertex->GetZ()) <= 10.);
  Bool_t lVertexNcontrib = (vertex->GetNContributors() > 0);

  if (!lVertexAcceptable)
  {
    PostData(1, fOutputList);
    return;
  }
  fEventStat->Fill(kVertexZ);
  /*
    if (!lVertexNcontrib){
        PostData(1, fOutputList);
        return;
    }*/
  fEventStat->Fill(kVtxNcontrib);

  // SPD pile-up in mult. bins
  //    if (fVEvent->IsPileupFromSPDInMultBins()){
  //        PostData(1, fOutputList);
  //        return;
  //    }
  fEventStat->Fill(kSPDPileUp);

  // spd vertex and its resolution
  const AliVVertex *vtxSPD = fVEvent->GetPrimaryVertexSPD();
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  //    if (TMath::Abs(vtxSPD->GetZ() - vertex->GetZ()) > 0.3){
  //        PostData(1, fOutputList);
  //        return;
  //    }
  fEventStat->Fill(kVtxDisplace);

  //    if (vtxSPD->IsFromVertexerZ() && zRes > 0.3 ) {
  //        PostData(1, fOutputList);
  //        return;
  //    }
  fEventStat->Fill(kVtxRes);

  // fill event histograms
  fVertexZ->Fill(vertex->GetZ());
  fVtxContrib->Fill(vertex->GetNContributors());
  fSPDVtxResol->Fill(zRes);
  fVtxDisplacement->Fill(vtxSPD->GetZ(), vertex->GetZ());

  // V0 information
  Int_t MultV0A = 0;
  Int_t MultV0C = 0;
  //    for(Int_t i=0; i<32; ++i) {
  //        MultV0A += vzeroData->GetMultiplicityV0A(i);
  //        MultV0C += vzeroData->GetMultiplicityV0C(i);
  //    }
  //
  //    fMultV0->Fill(MultV0A + MultV0C);
  //
  // centrality distribution
  /*
    Float_t lPercentile = 300;
    AliMultSelection *MultSelection = 0x0;
    //MultSelection = (AliMultSelection * ) fVEvent->FindListObject("MultSelection");
    if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    }else{
        lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
    }
    fCentV0M->Fill(lPercentile);*/

  // track loop
  int nMCprimaries = mcEvent->GetNumberOfPrimaries();
  //===== MC track loop ====
  for (int i = 0; i < nMCtracks; i++)
  {

    AliAODMCParticle *part = (AliAODMCParticle *)mcEvent->GetTrack(i);
    if (!part)
      continue;
    if (part->Pt() < 0.2 or part->Eta() > 0.8 or part->Eta() < -0.8)
      continue; //Gets rid of weird particles
    Int_t pdgCode = part->PdgCode();
    if (not(pdgCode == -2212 or pdgCode == 2212 or pdgCode == -1000010030 or pdgCode == 1000010030 or pdgCode == -1000020030 or pdgCode == 1000020030 or pdgCode == 1000010020 or pdgCode == -1000010020))
      continue;
    // if (not( antiproton 				or proton	 	 or antitriton 		 			or triton				or antihelium3			 or helium3 				or deuteron 			or antideuteron)) continue; Pseudocode for legibility

    Int_t label = part->Label();

    Bool_t isH3;
    Bool_t isPrimary = (label >= 0 and label <= nMCprimaries) and (part->IsPhysicalPrimary());
    Bool_t isSecondaryFromWeak = (part->IsSecondaryFromWeakDecay());
    Bool_t isSecondaryFromMaterial = (part->IsSecondaryFromMaterial());

    if (not(isPrimary or isSecondaryFromWeak or isSecondaryFromMaterial))
      continue;
    AliAODMCParticle *true_mother; //Will be assigned only if mother is found.

    //##### This checks for secondary He3 or He4  from material and reveals their processes ####
    /*
    if (isSecondaryFromMaterial and (pdgCode==1000020030 or pdgCode==-1000020030 or pdgCode==1000020040 or pdgCode==-1000020040)) {
        int motherlabel = part->GetMother();
        true_mother = (AliAODMCParticle*)mcEvent->GetTrack(motherlabel);
        if (!true_mother) {cout<<"GetTrack dont work"<<endl;}
        const char* PIDmother;
        const char* PIDdaughter;
        if (pdgCode==1000020030) PIDdaughter="He3";
        if (pdgCode==-1000020030) PIDdaughter="AHe3";
        if (true_mother->PdgCode()==1000020030) PIDmother="He3";
        if (true_mother->PdgCode()==-1000020030) PIDmother="AHe3";
        if (true_mother->PdgCode()==1000020040) {PIDmother="He4"; continue;}
        if (true_mother->PdgCode()==-1000020040) PIDmother="AHe4";
        cout<<"Particle was:"<<PIDdaughter<<endl;	
        if (not((true_mother->PdgCode()==1000020030) or (true_mother->PdgCode()==-1000020030)or (true_mother->PdgCode()==1000020040) or (true_mother->PdgCode()==-1000020040))) {
            cout<<"Mother particle was: "<<true_mother->PdgCode()<<endl;
        } else {
            cout<<"Mother particle was:"<<PIDmother<<endl;	
        }
        int ndaughters = true_mother->GetNDaughters();
        int firstd = true_mother->GetDaughterFirst();
        for (int nd=0; nd<ndaughters;nd++) {
            AliAODMCParticle* daughter = (AliAODMCParticle*)mcEvent->GetTrack(firstd+nd);
            if (!daughter) {cout<<"daughter "<<nd<<" doesn't exist"<<endl;continue;}
            cout<<"Daughter "<<nd<<" ID: "<<daughter->PdgCode()<<endl;
            if ((PIDmother=="He3") and ((daughter->PdgCode()==2212) or (daughter->PdgCode()==-2212) or (daughter->PdgCode()==2112) or (daughter->PdgCode()==-2112))) {cout<<endl<<"Mother momentum: "<<true_mother->P()<<endl; cout<<"particle momentum: "<<part->P()<<endl<<endl;}
        }
        cout<<endl<<endl;
    }*/

    Int_t pdgAbs = TMath::Abs(pdgCode);
    //cut down to barebones
    Bool_t testbool = kFALSE; //This is to test if the Ntuple gets filled correctly by filling it with protons.
    switch (pdgCode)
    {
    /*
        case 2212://changed to \Dueterons for analysis
            //cout<<"Proton"<<endl;
            
            if (isPrimary) {FillHistosMC((TList*)fMCOutputProtons->At(0), part, fAODEvent, i); if (testbool) {FillNtuple((TNtuple*)((TList*)((TList*)fMCOutputH3->At(2)))->At(44), part, fAODEvent, i);}}
            if (isSecondaryFromWeak) FillHistosMC((TList*)fMCOutputProtons->At(1), part, fAODEvent, i);
            if (isSecondaryFromMaterial) FillHistosMC((TList*)fMCOutputProtons->At(2), part, fAODEvent, i);
            break;
        case -2212:
            //cout<<"AProton"<<endl;
            if (isPrimary) FillHistosMC((TList*)fMCOutputAProtons->At(0), part, fAODEvent, i);
            if (isSecondaryFromWeak) FillHistosMC((TList*)fMCOutputAProtons->At(1), part, fAODEvent, i);
            if (isSecondaryFromMaterial) FillHistosMC((TList*)fMCOutputAProtons->At(2), part, fAODEvent, i);
            break;
        case 1000010030:
            //cout<<"Triton"<<endl;
            if (isPrimary) FillHistosMC((TList*)fMCOutputDeuterons->At(0), part, fAODEvent, i);
            if (isSecondaryFromWeak) FillHistosMC((TList*)fMCOutputDeuterons->At(1), part, fAODEvent, i);
            if (isSecondaryFromMaterial) FillHistosMC((TList*)fMCOutputDeuterons->At(2), part, fAODEvent, i);
            break;
        case -1000010030:
            //cout<<"ATriton"<<endl;
            if (isPrimary) FillHistosMC((TList*)fMCOutputADeuterons->At(0), part, fAODEvent, i);
            if (isSecondaryFromWeak) FillHistosMC((TList*)fMCOutputADeuterons->At(1), part, fAODEvent, i);
            if (isSecondaryFromMaterial) FillHistosMC((TList*)fMCOutputADeuterons->At(2), part, fAODEvent, i);
            break;
            */
    case 1000010030:
      //cout<<"TRITON"<<endl;
      if (isPrimary)
      {
        FillHistosMC((TList *)fMCOutputH3->At(0), part, fAODEvent, i);
      }
      if (isSecondaryFromWeak)
        FillHistosMC((TList *)fMCOutputH3->At(1), part, fAODEvent, i);
      if (isSecondaryFromMaterial)
      {
        FillHistosMC((TList *)fMCOutputH3->At(2), part, fAODEvent, i);
        FillNtuple((TNtuple *)((TList *)((TList *)fMCOutputH3->At(2)))->At(44), part, fAODEvent, i);
      }
      break;
    case -1000010030:
      //cout<<"anti-TRITON"<<endl;
      if (isPrimary)
        FillHistosMC((TList *)fMCOutputAH3->At(0), part, fAODEvent, i);
      if (isSecondaryFromWeak)
        FillHistosMC((TList *)fMCOutputAH3->At(1), part, fAODEvent, i);
      if (isSecondaryFromMaterial)
        FillHistosMC((TList *)fMCOutputAH3->At(2), part, fAODEvent, i);
      break;
      /*
        case 1000010020:
            //cout<<"Deuteron"<<endl;
            if (isPrimary) FillHistosMC((TList*)fMCOutputTriton->At(0), part, fAODEvent, i);
            if (isSecondaryFromWeak) FillHistosMC((TList*)fMCOutputTriton->At(1), part, fAODEvent, i);
            if (isSecondaryFromMaterial) FillHistosMC((TList*)fMCOutputTriton->At(2), part, fAODEvent, i);
            break;
            */
    }

  } // track loop

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisH3MC::Terminate(Option_t *)
{
  cout << "Analysis Task Terminate Called" << endl;
  // terminate
  // called at the END of the analysis (when all events are processed)
}

//________________________________________________________________________
void AliAnalysisH3MC::CreateHistosTrack(vector<TH1 *> &histos)
{

  /* This function creates the hitograms which are pushed into the vector and can then be accessed (i.e. filled) */
  /* ####### histos vector histogram order #######
     *
     * fHistP					0
     * fHistPt					1
     * fHistEtaPhi				2
     * fHistTPCCrossedRows		3
     * fHistTPCClusters			4
     * fHistTPCCRoverFind		5
     * fHistTPCFracShared		6
     * fHistTPCSignalN			7
     * fHistTPCSignalNfrac		8
     * fHistITSnCls				9
     * fHistChi2				10
     * fHistDCAxy				11
     * fHistDCAz				12
     * fHistDCAxyDCAz			13
     * fHistTPCSignal			14
     * fHistTPCnSigmaProt		15
     * fHistTPCnSigmaDeut		16
     * fHistTPCnSigmaH3			17
     * fHistTOFSignal			18
     * fHistTOFnSigmaProt		19
     * fHistTOFnSigmaDeut		20
     * fHistTOFnSigmaPion		21
     * fHistTOFnSigmaH3			22
     * fHistTOFmass2_DCAxy_p	23
     * fHistITSSignal			24
     * fHistITSnSigmaProt		25
     * fHistITSnSigmaDeut		26
     * fHistITSnSigmaH3			27
     * fHistPvPt				28
     *
     */

  TH1F *fHistP = new TH1F("fHistP", "p distribution;#it{p}/Z, GeV/#it{c};counts", 100, 0., 10.);
  fHistP->Sumw2();
  histos.push_back(fHistP);

  TH1F *fHistPt = new TH1F("fHistPt", "pt distribution;#it{p}_{T}, GeV/#it{c};counts", 100, 0., 10.);
  fHistPt->Sumw2();
  histos.push_back(fHistPt);

  TH2F *fHistEtaPhi = new TH2F("fHistEtaPhi", "eta-phi distribution;#eta;#varphi", 100, -1.0, 1.0, 360, 0, 2 * TMath::Pi());
  histos.push_back(fHistEtaPhi);

  // tracking
  TH1F *fHistTPCCrossedRows = new TH1F("fHistTPCCrossedRows", "TPC crossed rows;N TPC crossed rows;counts", 200, 0, 200);
  histos.push_back(fHistTPCCrossedRows);

  TH1F *fHistTPCClusters = new TH1F("fHistTPCClusters", "TPC clusters;N TPC clusters;counts", 200, 0, 200);
  histos.push_back(fHistTPCClusters);

  TH1F *fHistTPCCRoverFind = new TH1F("fHistTPCCRoverFind", "TPC crossed rows / findable;TPC crossed rows / findable;counts", 150, 0, 1.5);
  histos.push_back(fHistTPCCRoverFind);

  TH1F *fHistTPCFracShared = new TH1F("fHistTPCFracShared", "TPC shared clusters;TPC shared clusters;counts", 250, -1.0, 1.5);
  histos.push_back(fHistTPCFracShared);

  TH1F *fHistTPCSignalN = new TH1F("fHistTPCSignalN", "Number of PID clusters TPC;N of TPC PID clusters;counts", 200, 0., 200.);
  histos.push_back(fHistTPCSignalN);

  TH1F *fHistTPCSignalNfrac = new TH1F("fHistTPCSignalNfrac", "Fraction of PID clusters TPC;Fraction of TPC PID clusters;counts", 120, 0., 1.2);
  histos.push_back(fHistTPCSignalNfrac);

  TH1F *fHistITSnCls = new TH1F("fHistITSnCls", "Number of ITS clusters;N of ITS clusters;counts", 10, 0., 10.);
  histos.push_back(fHistITSnCls);

  TH1F *fHistChi2 = new TH1F("fHistChi2", "track chi2;#chi^{2};counts", 100, 0, 10);
  histos.push_back(fHistChi2);

  TH1F *fHistDCAxy = new TH1F("fHistDCAxy", "DCA xy;DCA_{xy};counts", 400, -2.0, 2.0);
  histos.push_back(fHistDCAxy);

  TH1F *fHistDCAz = new TH1F("fHistDCAz", "DCA z;DCA_{z};counts", 400, -2.0, 2.0);
  histos.push_back(fHistDCAz);

  TH2F *fHistDCAxyDCAz = new TH2F("fHistDCAxyDCAz", "DCAxy vs DCAz;DCA_{xy};DCA_{z}", 400, -2.0, 2.0, 400, -2.0, 2.0);
  histos.push_back(fHistDCAxyDCAz);

  // PID histos
  TH2F *fHistTPCSignal = new TH2F("fHistTPCSignal", "TPC dE/dx;#it{p}/Z, GeV/#it{c};TPC dE/dx", 500, 0., 10.0, 1500, 0., 1500.);
  histos.push_back(fHistTPCSignal);

  TH2F *fHistTPCnSigmaProt = new TH2F("fHistTPCnSigmaProt", "TPC nSigma prot;#it{p}/Z, GeV/#it{c};TPCn#sigma_{prot}", 100, 0., 10.0, 100, -5., 5.);
  histos.push_back(fHistTPCnSigmaProt);

  TH2F *fHistTPCnSigmaDeut = new TH2F("fHistTPCnSigmaDeut", "TPC nSigma deut;#it{p}/Z, GeV/#it{c};TPCn#sigma_{deut}", 100, 0., 10.0, 100, -5., 5.);
  histos.push_back(fHistTPCnSigmaDeut);

  TH2F *fHistTPCnSigmaH3 = new TH2F("fHistTPCnSigmaH3", "TPC nSigma H3;#it{p}/Z, GeV/#it{c};TPCn#sigma_{H3}", 100, 0., 10.0, 100, -5., 5.);
  histos.push_back(fHistTPCnSigmaH3);

  TH2F *fHistTOFSignal = new TH2F("fHistTOFSignal", "TOF signal;#it{p}/Z, GeV/#it{c};TOF signal", 500, 0., 10.0, 200, 0., 2.);
  histos.push_back(fHistTOFSignal);

  TH2F *fHistTOFnSigmaProt = new TH2F("fHistTOFnSigmaProt", "TOF nSigma prot;#it{p}/Z, GeV/#it{c};TOFn#sigma_{prot}", 100, 0., 10.0, 100, -5., 5.);
  histos.push_back(fHistTOFnSigmaProt);

  TH2F *fHistTOFnSigmaDeut = new TH2F("fHistTOFnSigmaDeut", "TOF nSigma deut;#it{p}/Z, GeV/#it{c};TOFn#sigma_{deut}", 100, 0., 10.0, 100, -5., 5.);
  histos.push_back(fHistTOFnSigmaDeut);

  TH2F *fHistTOFnSigmaPion = new TH2F("fHistTOFnSigmaPion", "TOF nSigma pion;#it{p}/Z, GeV/#it{c};TOFn#sigma_{pion}", 100, 0., 10.0, 100, -5., 5.);
  histos.push_back(fHistTOFnSigmaPion);

  TH2F *fHistTOFnSigmaH3 = new TH2F("fHistTOFnSigmaH3", "TOF nSigma H3;#it{p}/Z, GeV/#it{c};TOFn#sigma_{H3}", 100, 0., 10.0, 100, -5., 5.);
  histos.push_back(fHistTOFnSigmaH3);

  TH3F *fHistTOFmass2_DCAxy_p = new TH3F("fHistTOFmass2_DCAxy_p", "TOF mass2 vs DCA_{xy} vs #it{p}/Z;#it{p}/Z, GeV/#it{c};TOF mass^{2} / Z^{2}, (GeV/#it{c}^{2})^{2};DCA_{xy}, cm;", 50, 0., 5., 50, 0.0, 10.0, 100, -1.0, 1.0);
  histos.push_back(fHistTOFmass2_DCAxy_p);

  TH2F *fHistITSSignal = new TH2F("fHistITSSignal", "ITS signal;#it{p}/Z, GeV/#it{c};ITS signal", 500, 0., 10.0, 750, 0., 750.);
  histos.push_back(fHistITSSignal);

  TH2F *fHistITSnSigmaProt = new TH2F("fHistITSnSigmaProt", "ITS nSigma prot;#it{p}/Z, GeV/#it{c};ITSn#sigma_{prot}", 100, 0., 10.0, 200, -10., 10.);
  histos.push_back(fHistITSnSigmaProt);

  TH2F *fHistITSnSigmaDeut = new TH2F("fHistITSnSigmaDeut", "ITS nSigma deut;#it{p}/Z, GeV/#it{c};ITSn#sigma_{deut}", 100, 0., 10.0, 200, -10., 10.);
  histos.push_back(fHistITSnSigmaDeut);

  TH2F *fHistITSnSigmaH3 = new TH2F("fHistITSnSigmaH3", "ITS nSigma H3;#it{p}/Z, GeV/#it{c};ITSn#sigma_{H3}", 100, 0., 10.0, 200, -10., 10.);
  histos.push_back(fHistITSnSigmaH3);

  TH2F *fHistDCAxyP = new TH2F("fHistDCAxyP_wDCAzcut", "fHistDCAxyP_wDCAzcut", 2, -1.0, 1.0, 4, 0, 4);
  fHistDCAxyP->GetXaxis()->SetTitle("DCAxy");
  fHistDCAxyP->GetYaxis()->SetTitle("#it{p}/Z GeV/#it{c}");
  histos.push_back(fHistDCAxyP);

  TH2F *fHistDCAzP = new TH2F("fHistDCAxyP", "fHistDCAxyP", 2, -1.0, 1.0, 1, 0, 4);
  fHistDCAxyP->GetXaxis()->SetTitle("DCAxy");
  fHistDCAxyP->GetYaxis()->SetTitle("#it{p}/Z GeV/#it{c}");
  histos.push_back(fHistDCAzP);

  TH2F *fHistChi2TPCConstrainedGlobal = new TH2F("fHistChi2TPCConstrainedGlobal", "fHistChi2TPCConstrainedGlobal", 400, 0, 8, 500, 0, 100);
  fHistChi2TPCConstrainedGlobal->GetXaxis()->SetTitle("p_Vtx GeV/c");
  fHistChi2TPCConstrainedGlobal->GetYaxis()->SetTitle("Chi2 TPC Constrained global");
  histos.push_back(fHistChi2TPCConstrainedGlobal);

  TH2F *fHistTPCSharedClusters = new TH2F("fHistTPCSharedClusters", "fHistTPCSharedClusters", 400, 0, 8, 500, 0, 1);
  fHistTPCSharedClusters->GetXaxis()->SetTitle("p_Vtx GeV/c");
  fHistTPCSharedClusters->GetYaxis()->SetTitle("TPC shared clusters fraction");
  histos.push_back(fHistTPCSharedClusters);

  TH3F *fHistTOFmass2_DCAxy_pt = new TH3F("fHistTOFmass2_DCAxy_pt", "TOF mass2 vs DCA_{xy} vs #it{pt}/Z;#it{pt}/Z, GeV/#it{c};TOF mass^{2} / Z^{2}, (GeV/#it{c}^{2})^{2};DCA_{xy}, cm;", 1, 0., 5., 1, 0.0, 10.0, 1, -1.0, 1.0);
  histos.push_back(fHistTOFmass2_DCAxy_pt);

  TH3F *fHistDCAxyzP = new TH3F("fHistDCAxyzP", "fHistDCAxyzP;DCAz (cm);DCAxy (cm);#{it}/Z GeV/#it{c}", 2, -1, 1, 2, -1, 1, 4, 0, 4);
  histos.push_back(fHistDCAxyzP); //This causes problems
}

//________________________________________________________________________
void AliAnalysisH3MC::FillHistosTrack(vector<TH1 *> &histos, AliAODTrack *track)
{
  Double_t trackP = track->P();
  Double_t trackPCC = track->P() * track->Charge();
  Double_t trackPt = track->Pt();
  Double_t trackEta = track->Eta();
  Double_t trackPhi = track->Phi();

  // TOF information
  Float_t fbetaTOF = GetTOFBeta(track);
  Float_t fmass2TOF = GetMass2TOF(fbetaTOF, track);
  Float_t fmass2TOFCC = GetMass2TOF(fbetaTOF, track) * track->Charge() * track->Charge();

  // PID information
  Float_t nSigmaTPCprot = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  Float_t nSigmaITSprot = fPIDResponse->NumberOfSigmasITS(track, AliPID::kProton);
  Float_t nSigmaTOFprot = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);

  Float_t nSigmaTPCdeut = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kDeuteron);
  Float_t nSigmaITSdeut = fPIDResponse->NumberOfSigmasITS(track, AliPID::kDeuteron);
  Float_t nSigmaTOFdeut = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kDeuteron);

  Float_t nSigmaTPCH3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kTriton);
  Float_t nSigmaITSH3 = fPIDResponse->NumberOfSigmasITS(track, AliPID::kTriton);
  Float_t nSigmaTOFH3 = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kTriton);

  Float_t nSigmaTOFpion = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);

  // DCA information
  double DCAxy = -99.;
  double DCAz = -99.;
  double dcaVals[2] = {-99., -99.};
  double covar[3] = {0., 0., 0.};
  AliAODTrack copy(*track);
  const AliVVertex *AODeventVtx = copy.GetAODEvent()->GetPrimaryVertex();
  Double_t AODeventMagneticF = copy.GetAODEvent()->GetMagneticField();
  if (copy.PropagateToDCA(AODeventVtx, AODeventMagneticF, 10, dcaVals, covar))
  {
    DCAxy = dcaVals[0];
    DCAz = dcaVals[1];
  }
  else
  {
    DCAxy = -99.; // track->DCA();
    DCAz = -99.;  // track->ZAtDCA();
  }

  // if DCAxy > cut for final analysis -> fill only 3d histogram
  if (kFALSE and (TMath::Abs(DCAxy) > fMaxDCAxyAna))
  { //currently disabled to check what is going on with MC

    ((TH3F *)(histos.at(23)))->Fill(trackP, fmass2TOF, DCAxy);
    return;
  }
  else
  {

    histos.at(0)->Fill(trackP);
    /*
        histos.at(1)->Fill(trackPt);
        histos.at(2)->Fill(trackEta,trackPhi);
        histos.at(3)->Fill(track->GetTPCClusterInfo(2, 1));
        histos.at(4)->Fill(track->GetTPCNcls());
        if (!track->GetTPCNclsF()) {
            histos.at(5)->Fill(0);
        } else {
            histos.at(5)->Fill(track->GetTPCClusterInfo(2, 1)/
                               double(track->GetTPCNclsF()));
        }
        if (!track->GetTPCNcls()){
            histos.at(6)->Fill(-1);
        } else {
            histos.at(6)->Fill(track->GetTPCnclsS()/double(track->GetTPCNcls()));
        }
        histos.at(7)->Fill(track->GetTPCsignalN());
        if (!track->GetTPCNcls()){
            histos.at(8)->Fill(-1);
        } else {
            histos.at(8)->Fill(track->GetTPCsignalN()/double(track->GetTPCNcls())); 
        }
        histos.at(9)->Fill(track->GetITSNcls());
        histos.at(10)->Fill(track->Chi2perNDF());
        histos.at(11)->Fill(DCAxy);
        histos.at(12)->Fill(DCAz);
        histos.at(13)->Fill(DCAxy,DCAz);
        
        // fill PID histos
        histos.at(14)->Fill(trackP,track->GetTPCsignal());
        histos.at(15)->Fill(trackP,nSigmaTPCprot);
        histos.at(16)->Fill(trackP,nSigmaTPCdeut);
    histos.at(17)->Fill(trackP,nSigmaTPCHe3);
        histos.at(18)->Fill(trackP,fbetaTOF);
        histos.at(19)->Fill(trackP,nSigmaTOFprot);
        histos.at(20)->Fill(trackP,nSigmaTOFdeut);
        histos.at(21)->Fill(trackP,nSigmaTOFpion);
    histos.at(22)->Fill(trackP,nSigmaTOFHe3);
        ((TH3F*)(histos.at(23)))->Fill(trackP,fmass2TOF,DCAxy);
        histos.at(24)->Fill(trackP,track->GetITSsignal());
    //cout<<nSigmaITSprot<<endl;
        histos.at(25)->Fill(trackP,nSigmaITSprot);
        histos.at(26)->Fill(trackP,nSigmaITSdeut);
    histos.at(27)->Fill(trackP,nSigmaITSHe3);
    if (TMath::Abs(DCAz)<0.1) {histos.at(28)->Fill(DCAxy, trackP);}
    histos.at(29)->Fill(DCAxy, trackP);
    histos.at(30)->Fill(trackP, track->GetChi2TPCConstrainedVsGlobal());
    histos.at(31)->Fill(trackP, track->GetTPCnclsS());
        ((TH3F*)(histos.at(32)))->Fill(trackPt,fmass2TOF,DCAxy);
        ((TH3F*)(histos.at(33)))->Fill(DCAz, DCAxy, trackP);
    */
    return;
  }
}

//________________________________________________________________________
Float_t AliAnalysisH3MC::GetTOFBeta(AliAODTrack *track)
{
  float beta = -999;
  double integratedTimes[9] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};

  track->GetIntegratedTimes(integratedTimes);

  const float c = 2.99792457999999984e-02;
  float p = track->P();
  float l = integratedTimes[0] * c;

  float trackT0 = fPIDResponse->GetTOFResponse().GetStartTime(p);

  float timeTOF = track->GetTOFsignal() - trackT0;
  if (timeTOF > 0)
  {
    beta = l / timeTOF / c;
  }
  return beta;
}

//________________________________________________________________________
Float_t AliAnalysisH3MC::GetMass2TOF(Float_t beta, AliAODTrack *track)
{
  Float_t p = track->P();
  Float_t mass2sq = -999;
  if (!(beta == 0))
  {
    mass2sq = ((1 / (beta * beta)) - 1) * (p * p);
  }
  return mass2sq;
}

void AliAnalysisH3MC::CreateMCLists(TList *list)
{

  vector<TH1 *> primhists;
  vector<TH1 *> secWhists;
  vector<TH1 *> secMhists;

  cout << "size of primhists @ 0: " << primhists.size() << endl;
  CreateHistosTrack(primhists);
  CreateHistosTrack(secWhists);
  CreateHistosTrack(secMhists);
  cout << "size of primhists @ 1: " << primhists.size() << endl;
  CreateHistosMC(primhists);
  CreateHistosMC(secWhists);
  CreateHistosMC(secMhists);
  cout << "size of primhists @ 2: " << primhists.size() << endl;

  cout << endl
       << "HistsVectorSize = " << fMCHistVector.size() << endl;
  fMCHistVector.push_back(primhists);
  fMCHistVector.push_back(secWhists);
  fMCHistVector.push_back(secMhists);

  for (int i = 0; i < primhists.size(); i++)
  {
    primaries->Add(primhists.at(i));
    secondariesW->Add(secWhists.at(i));
    secondariesM->Add(secMhists.at(i));
  }
  cout << "PList size: " << list->GetSize() << endl;
  list->Add(primaries);
  cout << "PList size: " << list->GetSize() << endl;
  list->Add(secondariesW);
  cout << "PList size: " << list->GetSize() << endl;
  list->Add(secondariesM);
  cout << "PList size: " << list->GetSize() << endl;
}
void AliAnalysisH3MC::CreateHistosMC(vector<TH1 *> &hists)
{
  /* MC histograms. 
     * 
     * DCA is covered by Reco Output
     *
     * p, pTPC and pt spectra for crosschecks and ratio
     * pTOF/pTPC
     * p* / pVTX 
     */

  TH1F *pVTX = new TH1F("pVTX", "pVTX", 500, 0, 10);
  pVTX->GetXaxis()->SetTitle("pVTX GeV/c");
  hists.push_back(pVTX);

  TH1F *pT = new TH1F("pT", "pT", 500, 0, 10);
  pT->GetXaxis()->SetTitle("pT GeV/c");
  hists.push_back(pT);

  TH1F *pTPC = new TH1F("pTPC", "pTPC", 500, 0, 10);
  pTPC->GetXaxis()->SetTitle("pTPC GeV/c");
  hists.push_back(pTPC);

  TH2F *ELossTOF = new TH2F("ELossTOF", "ELossTOF", 500, 0, 10, 110, 0, 1.1);
  ELossTOF->GetXaxis()->SetTitle("pTPC GeV/c");
  ELossTOF->GetYaxis()->SetTitle("pTOF/pTPC");
  hists.push_back(ELossTOF);

  TH2F *ELossITS = new TH2F("ELossTPC", "ELossTPC", 500, 0, 10, 110, 0, 1.1);
  ELossITS->GetXaxis()->SetTitle("pVTX");
  ELossITS->GetYaxis()->SetTitle("pTPC/pVTX");
  hists.push_back(ELossITS);

  TH1F *trackcuts = new TH1F("trackcuts", "trackcuts", 5, 0, 5);
  trackcuts->GetXaxis()->SetBinLabel(1, "MC truth");
  trackcuts->GetXaxis()->SetBinLabel(2, "After Filterbit cut");
  trackcuts->GetXaxis()->SetBinLabel(3, "After ITS N clusters cut");
  trackcuts->GetXaxis()->SetBinLabel(4, "Labels not found");
  trackcuts->GetXaxis()->SetBinLabel(5, "More than 1 label found");
  hists.push_back(trackcuts);

  TH2F *energy_loss = new TH2F("ELoss", "ELoss", 500, 0, 10, 400, -2, 2);
  energy_loss->GetXaxis()->SetTitle("p_Vtx_reco");
  energy_loss->GetYaxis()->SetTitle("p_gen - p_Vtx_reco");
  hists.push_back(energy_loss);

  TH2F *energy_loss2 = new TH2F("ELoss2", "ELoss2", 500, 0, 10, 400, -2, 2);
  energy_loss2->GetXaxis()->SetTitle("p_Vtx_gen");
  energy_loss2->GetYaxis()->SetTitle("p_gen - p_Vtx_reco");
  hists.push_back(energy_loss2);

  TH1F *TPC_yields = new TH1F("TPC_yield", "TPC_yield", 400, 0, 4);
  TPC_yields->GetXaxis()->SetTitle("p GeV/c");
  TPC_yields->GetYaxis()->SetTitle("yield in TPC");
  hists.push_back(TPC_yields);

  TH1F *TOF_yields = new TH1F("TOF_yield", "TOF_yield", 400, 0, 4);
  TOF_yields->GetXaxis()->SetTitle("p GeV/c");
  TOF_yields->GetYaxis()->SetTitle("yield in TOF");
  hists.push_back(TOF_yields);
}

void AliAnalysisH3MC::FillHistosMC(TList *list, AliAODMCParticle *part, AliAODEvent *fAODEvent, Int_t iMC)
{
  AliMCParticle *mcpart = (AliMCParticle *)mcEvent->GetTrack(iMC);
  //AliMCParticle* mcpart = (AliMCParticle*)part;
  vector<TH1 *> histos;
  for (int i = 0; i < list->GetSize(); i++)
  {
    histos.push_back((TH1 *)list->At(i));
  }
  int nRefs = mcpart->GetNumberOfTrackReferences();
  float pVTX = part->P();
  float pT = part->Pt();
  float pITS = 0.f;
  float pTPC = 0.f;
  float pTRD = 0.f;
  float pTOF = 0.f;
  float pStarRef = 0.f;
  float rStar = 0.f;
  bool disappeared = kFALSE;
  bool referenceWorked = kFALSE;
  ;

  Int_t label = mcpart->GetLabel();
  //cout<<"Number of track references: "<<nRefs<<endl<<endl;

  /*
    for (int iRef = nRefs+1; iRef<nRefs;iRef++) {
        cout<<"iRef: "<<iRef<<endl;
        if (!mcpart) {
            cout<<"no MC particle"<<endl; 
            continue;}
        AliTrackReference* ref = (AliTrackReference*)mcpart->GetTrackReference(iRef);
        if (!ref) {
            cout<<"no reference"<<endl; 
            continue;
        } else {referenceWorked=kTRUE;}
        if (ref->DetectorId()==AliTrackReference::kITS) pITS = ref->P();
        if (ref->DetectorId()==AliTrackReference::kTPC) pTPC = ref->P();
        if (ref->DetectorId()==AliTrackReference::kTRD) pTRD = ref->P();
        if (ref->DetectorId()==AliTrackReference::kTOF) pTOF = ref->P();
        if (ref->DetectorId()==AliTrackReference::kDisappeared) {
            rStar = ref->R();
            pStarRef = ref->P();
            disappeared = kTRUE;
        }
    }*/
  //if (!referenceWorked) {AliWarning("No references for this track"); cout<<"No reference for this track"<<endl;return;}

  int i = 34;
  histos.at(i)->Fill(pVTX);
  histos.at(i + 1)->Fill(pT);

  histos.at(i + 5)->Fill("MC truth", 1);
  //cout<<endl<<endl<<histos.size()<<endl<<endl<<endl;

  if (kFALSE and part->IsSecondaryFromMaterial())
  {
    cout << part->GetMother() << endl;
  } //need to acess and identify mother particle, that should give more insight.

  int GotSelected = 0;
  if (!fAODEvent)
  {
    cout << "fAODEvent not Found!" << endl;
    return;
  }
  for (int itrack = 0; itrack < fAODEvent->GetNumberOfTracks(); itrack++)
  {
    AliAODTrack *track = (AliAODTrack *)fAODEvent->GetTrack(itrack);

    if (!track)
    {
      cout << "Error: could not find track " << itrack << endl;
      continue;
    }
    Int_t labelReco = TMath::Abs(track->GetLabel());
    if (kFALSE and part->IsSecondaryFromMaterial() and kFALSE)
    {
      cout << endl
           << "MC label: " << iMC << endl;
      cout << "Reco label: " << track->GetLabel() << endl;
    }
    //if (not(labelReco==label)) {continue;} else {cout<<"Labels matched"<<endl;}
    if (not(labelReco == iMC))
      continue;
    pTPC = track->GetTPCmomentum();
    //if (pTPC>0.0001) {histos.at(i+2)->Fill(pTPC);} else {AliWarning("could not read pTPC");}
    //if (pTPC>0.0001 and pTOF>0.0001) {histos.at(i+3)->Fill(pTPC, pTOF/pTPC);cout<<"pTOF done"<<endl;} else if (disappeared) {AliWarning("particle disappeared before TOF");} else {AliFatal("Particle has no TOF momentum but didnt disappear...");}
    GotSelected++;
    if (!(track->TestFilterBit(fFilterBit)))
      continue;
    //if (!(fTrackCuts->IsSelected(track))) continue;
    histos.at(i + 5)->Fill("After Filterbit cut", 1);
    if (track->GetITSNcls() < fMinClIts)
      continue;
    histos.at(i + 5)->Fill("After ITS N clusters cut", 1);
    //if (TMath::Abs(track->Y())>0.5) continue;
    if (TMath::Abs(track->Eta()) > 0.8)
      continue;

    histos.at(i + 4)->Fill(pVTX, track->GetTPCmomentum() / pVTX);
    histos.at(i + 6)->Fill(track->P(), (track->P()) - pVTX);
    histos.at(i + 7)->Fill(pVTX, (track->P()) - pVTX);
    histos.at(i + 8)->Fill(track->P());

    double TOFm2 = GetMass2TOF(GetTOFBeta(track), track);
    Bool_t reachedTOF = (TOFm2 > 7) and (TOFm2 < 11);
    if (reachedTOF)
    {
      histos.at(i + 9)->Fill(track->P());
    }

    double DCAxy = -99.;
    double DCAz = -99.;
    double dcaVals[2] = {-99., -99.};
    double covar[3] = {0., 0., 0.};
    AliAODTrack copy(*track);
    const AliVVertex *AODeventVtx = copy.GetAODEvent()->GetPrimaryVertex();
    Double_t AODeventMagneticF = copy.GetAODEvent()->GetMagneticField();
    if (copy.PropagateToDCA(AODeventVtx, AODeventMagneticF, 10, dcaVals, covar))
    {
      DCAxy = dcaVals[0];
      DCAz = dcaVals[1];
    }
    else
    {
      DCAxy = -99.; // track->DCA();
      DCAz = -99.;  // track->ZAtDCA();
    }

    //if (TMath::Abs(DCAz ) > fMaxDCAz)        continue;
    //if (TMath::Abs(DCAxy) > fMaxDCAxyCut)    continue;

    FillHistosTrack(histos, track);
  }
  if (GotSelected == 0)
  {
    histos.at(i + 5)->Fill("Labels not found", 1);
  }
  if (GotSelected > 1)
  {
    histos.at(i + 5)->Fill("More than 1 label found", 1);
    //cout<<"Number of reconstructed tracks matched: "<<GotSelected<<endl;
  }
}

/*const char *GetID(int pdgCode)
{
  if (pdgCode == 111 or pdgCode == -111)
    return "Pi0";
  if (pdgCode == 211)
    return "Pi+";
  if (pdgCode == -211)
    return "Pi-";

  if (pdgCode == 2212)
    return "Proton";
  if (pdgCode == -2212)
    return "AProton";

  if (pdgCode == 1000010020)
    return "Deuteron";
  if (pdgCode == -1000010020)
    return "ADeuteron";

  if (pdgCode == 1000010030)
    return "H3";
  if (pdgCode == -1000010030)
    return "AH3";

  if (pdgCode == 1000020040)
    return "He4";
  if (pdgCode == 1000020040)
    return "AHe4";

  return (const char *)pdgCode;
}*/

void AliAnalysisH3MC::FillNtuple(TNtuple *nt, AliAODMCParticle *part, AliAODEvent *fAODEvent, Int_t iMC)
{
  AliMCParticle *mcpart = (AliMCParticle *)mcEvent->GetTrack(iMC);

  if (!fAODEvent)
  {
    cout << "fAODEvent not Found!" << endl;
    return;
  }
  for (int itrack = 0; itrack < fAODEvent->GetNumberOfTracks(); itrack++)
  {
    AliAODTrack *track = (AliAODTrack *)fAODEvent->GetTrack(itrack);

    if (!track)
    {
      cout << "Error: could not find track " << itrack << endl;
      continue;
    }
    Int_t labelReco = TMath::Abs(track->GetLabel());
    if (kFALSE and part->IsSecondaryFromMaterial() and kFALSE)
    {
      cout << endl
           << "MC label: " << iMC << endl;
      cout << "Reco label: " << track->GetLabel() << endl;
    }
    //if (not(labelReco==label)) {continue;} else {cout<<"Labels matched"<<endl;}
    if (not(labelReco == iMC))
      continue;
    if (!(track->TestFilterBit(fFilterBit)))
      continue;
    //if (!(fTrackCuts->IsSelected(track))) continue;
    if (track->GetITSNcls() < fMinClIts)
      continue;
    /*if (TMath::Abs(track->Y()) > 0.5)
      continue;*/
    if (TMath::Abs(track->Eta()) > 0.8)
      continue;

    double TOFm2 = GetMass2TOF(GetTOFBeta(track), track);
    Bool_t reachedTOF = (TOFm2 > 7) and (TOFm2 < 11); //mass of He3: 2.808391482 GeV/c2

    double DCAxy = -99.;
    double DCAz = -99.;
    double dcaVals[2] = {-99., -99.};
    double covar[3] = {0., 0., 0.};
    AliAODTrack copy(*track);
    const AliVVertex *AODeventVtx = copy.GetAODEvent()->GetPrimaryVertex();
    Double_t AODeventMagneticF = copy.GetAODEvent()->GetMagneticField();
    if (copy.PropagateToDCA(AODeventVtx, AODeventMagneticF, 10, dcaVals, covar))
    {
      DCAxy = dcaVals[0];
      DCAz = dcaVals[1];
    }
    else
    {
      DCAxy = -99.; // track->DCA();
      DCAz = -99.;  // track->ZAtDCA();
    }

    Float_t nSigmaTPCH3 = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kTriton);
    //if (TMath::Abs(DCAz ) > fMaxDCAz)        continue;
    //if (TMath::Abs(DCAxy) > fMaxDCAxyCut)    continue;

    //For reference //fNtupleH3 = new TNtuple("fNtupleH3", "fNtupleH3", "p:pt:TPCSignal:TPCnSigmaH3:DCAxy:DCAz:TOFm2:TPCNClusters:ITSNClusters:TPCClusters4dEdx:Eta:Chi2TPC:Chi2ITS:TPCCrossedRows");
    float vars[16];
    vars[0] = track->P();
    vars[1] = track->Pt();
    vars[2] = track->GetTPCsignal();
    vars[3] = nSigmaTPCH3;
    vars[4] = DCAxy;
    vars[5] = DCAz;
    vars[6] = TOFm2;
    vars[7] = track->GetTPCnclsS();
    vars[8] = track->GetITSNcls();
    vars[9] = track->GetTPCsignalN();
    vars[10] = track->Eta();
    vars[12] = track->Chi2perNDF();
    vars[13] = track->GetITSchi2();
    vars[14] = track->GetTPCClusterInfo(2, 1);
    vars[15] = track->GetLabel();
    nt->Fill(vars);
  }
}

