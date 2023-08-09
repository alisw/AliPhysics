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

//______________________________________________________________________________
// Analysis task for high pt particle correlations
// author: Jasper Parkila, D.J. Kim
// ALICE Group University of Jyvaskyla
// Finland
//////////////////////////////////////////////////////////////////////////////

#include <TGrid.h>
#include <TRandom.h>
#include <AliAnalysisUtils.h>
#include <AliAnalysisTaskSE.h>
#include <AliAODMCParticle.h>
#include <AliMCEvent.h>
#include <AliGenHijingEventHeader.h>
#include <AliGenDPMjetEventHeader.h>
#include <AliGenHepMCEventHeader.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisDataContainer.h>
#include <AliAODEvent.h>
#include <AliPWG0Helper.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
//#include <AliInputEventHandler.h>
#include <AliMultSelection.h>
#include "AliJCatalystTask.h"
#include "AliJTrack.h"
#include "AliJHistManager.h"
#include "AliJRunTable.h"
#include "AliJFFlucAnalysis.h" // TEMP for getting bins
#include "TDirectoryFile.h"
#include "TList.h"
#include "TSystem.h"

#include "AliVMultiplicity.h"
#include <string>

//#pragma GCC diagnostic warning "-Wall"
//______________________________________________________________________________
AliJCatalystTask::AliJCatalystTask():
  AliAnalysisTaskSE(),
  fInputList(0),
  fInputListALICE(0),
  //fOutput(0),
  fCentDetName("V0M"),
  paodEvent(0),
  fcent(-999),
  fZvert(-999),
  fnoCentBin(false),
  fDebugLevel(0),
  fEvtNum(0),
  fFilterBit(0),
  fNumTPCClusters(70),
  fEffMode(0),
  fEffFilterBit(0),
  fRunNum(-1),
  fPcharge(0),
  fEta_min(-0.8),
  fEta_max(0.8),
  fPt_min(0.2),
  fPt_max(5.0),
  fzvtxCut(10.0),
  fremovebadarea(kFALSE),
  fremovebadarea18q(kFALSE),
  flags(0),
  fJCatalystEntry(0),
  fIsGoodEvent(false),
  fJCorMapTask(NULL),
  fJCorMapTaskName("JCorrectionMapTask"),
  pPhiWeights(0),
  grEffCor(0),
  fCentBinEff(0),
  phiMapIndex(0),
  bUseAlternativeWeights(kFALSE),
  fAliEventCuts(NULL),
// QA part.
  fMainList(NULL),
  bSaveAllQA(kFALSE),
  bSaveHMOhist(kFALSE),
  bSaveQCNUA(kFALSE),
  fCentralityBins(16),
  fcent_0(0.), fcent_1(0.), fcent_2(0.), fcent_3(0.), fcent_4(0.), fcent_5(0.), fcent_6(0.), fcent_7(0.), fcent_8(0.), fcent_9(0.), 
  fcent_10(0.), fcent_11(0.), fcent_12(0.), fcent_13(0.), fcent_14(0.), fcent_15(0.), fcent_16(0.),
  fChi2perNDF_min(0.1), // TBC: Which default value do we choose?
  fChi2perNDF_max(4.0), // TBC: Same here, this is old 2010 value.
  fDCAxy_max(2.4),  // TBC: Shall we keep 2010 default?
  fDCAz_max(3.2), // TBC: Shall we keep 2010 default?
  bDCABaseCut(false), // DCABaseCut for FB 768
  fUseITSMinClusters(false),
  fITSMinClusters(2),
  fUseTightCuts(false),
  bUseDCAbaseCut(false),
  fAddESDpileupCuts(false),
  fESDpileup_slope(3.38), fESDpileup_inter(15000),
  fSaveESDpileupQA(false),
  fAddTPCpileupCuts(false),
  fSaveTPCpileupQA(false),
  fControlProfileList(NULL)
{
  InitializeArrays(); //
}

//______________________________________________________________________________
AliJCatalystTask::AliJCatalystTask(const char *name):
  AliAnalysisTaskSE(name),
  fInputList(0),
  fInputListALICE(0),
  //fOutput(0),
  fTaskName(name),
  fCentDetName("V0M"),
  paodEvent(0),
  fcent(-999),
  fZvert(-999),
  fnoCentBin(false),
  fDebugLevel(0),
  fEvtNum(0),
  fFilterBit(0),
  fNumTPCClusters(70),
  fEffMode(0),
  fEffFilterBit(0),
  fRunNum(-1),
  fPcharge(0),
  fEta_min(-0.8),
  fEta_max(0.8),
  fPt_min(0.2),
  fPt_max(5.0),
  fzvtxCut(10.0),
  fremovebadarea(kFALSE),
  fremovebadarea18q(kFALSE),
  flags(0),
  fJCatalystEntry(0),
  fIsGoodEvent(false),
  fJCorMapTask(NULL),
  fJCorMapTaskName("JCorrectionMapTask"),
  pPhiWeights(0),
  grEffCor(0),
  fCentBinEff(0),
  phiMapIndex(0),
  bUseAlternativeWeights(kFALSE),
  fAliEventCuts(NULL),
// QA part.
  fMainList(NULL),
  bSaveAllQA(kFALSE),
  bSaveHMOhist(kFALSE),
  bSaveQCNUA(kFALSE),
  fCentralityBins(16),
  fcent_0(0.), fcent_1(0.), fcent_2(0.), fcent_3(0.), fcent_4(0.), fcent_5(0.), fcent_6(0.), fcent_7(0.), fcent_8(0.), fcent_9(0.), 
  fcent_10(0.), fcent_11(0.), fcent_12(0.), fcent_13(0.), fcent_14(0.), fcent_15(0.), fcent_16(0.),
  fChi2perNDF_min(0.1), // TBC: Which default value do we choose?
  fChi2perNDF_max(4.0), // TBC: Same here, this is old 2010 value.
  fDCAxy_max(2.4),  // TBC: Shall we keep 2010 default?
  fDCAz_max(3.2), // TBC: Shall we keep 2010 default?
  bDCABaseCut(false),
  fUseITSMinClusters(false),
  fITSMinClusters(2),
  fUseTightCuts(false),
  bUseDCAbaseCut(false),
  fAddESDpileupCuts(false),
  fESDpileup_slope(3.38), fESDpileup_inter(15000),
  fSaveESDpileupQA(false),
  fAddTPCpileupCuts(false),
  fSaveTPCpileupQA(false),
  fControlProfileList(NULL)
{
// Main list to save the output of the QA.
  fMainList = new TList();
  fMainList->SetName("fJCatalystOutput");
  fMainList->SetOwner(kTRUE);

  InitializeArrays();

  //DefineOutput(1, TDirectory::Class()); // Uncommented in the current version on AliPhysics.
  DefineOutput(1, TList::Class());
}

//____________________________________________________________________________
AliJCatalystTask::AliJCatalystTask(const AliJCatalystTask& ap) :
  AliAnalysisTaskSE(ap.GetName()),
  fInputList(ap.fInputList),
  fInputListALICE(ap.fInputListALICE),
  //fOutput(ap.fOutput),
  fcent(ap.fcent),
  fZvert(ap.fZvert),
  fnoCentBin(ap.fnoCentBin),
  fRunNum(ap.fRunNum),
  fJCorMapTask(ap.fJCorMapTask),
  fJCorMapTaskName(ap.fJCorMapTaskName),
  pPhiWeights(ap.pPhiWeights),
  grEffCor(ap.grEffCor),
  fCentBinEff(ap.fCentBinEff),
  phiMapIndex(ap.phiMapIndex),
  bUseAlternativeWeights(ap.bUseAlternativeWeights),
  fAliEventCuts(ap.fAliEventCuts),
// QA part.
  fMainList(ap.fMainList),
  bSaveAllQA(ap.bSaveAllQA),
  bSaveHMOhist(ap.bSaveHMOhist),
  bSaveQCNUA(ap.bSaveQCNUA),
  fCentralityBins(ap.fCentralityBins),
  fcent_0(ap.fcent_0), fcent_1(ap.fcent_1), fcent_2(ap.fcent_2), fcent_3(ap.fcent_3), fcent_4(ap.fcent_4), fcent_5(ap.fcent_5), fcent_6(ap.fcent_6), fcent_7(ap.fcent_7), fcent_8(ap.fcent_8), fcent_9(ap.fcent_9), fcent_10(ap.fcent_10), fcent_11(ap.fcent_11), fcent_12(ap.fcent_12), fcent_13(ap.fcent_13), fcent_14(ap.fcent_14), fcent_15(ap.fcent_15), fcent_16(ap.fcent_16),
  fChi2perNDF_min(ap.fChi2perNDF_min), fChi2perNDF_max(ap.fChi2perNDF_max),
  fDCAxy_max(ap.fDCAxy_max), fDCAz_max(ap.fDCAz_max),
  bDCABaseCut(ap.bDCABaseCut),
  fUseITSMinClusters(ap.fUseITSMinClusters),fITSMinClusters(ap.fITSMinClusters),
  fUseTightCuts(ap.fUseTightCuts), bUseDCAbaseCut(ap.bUseDCAbaseCut),
  fAddESDpileupCuts(ap.fAddESDpileupCuts),
  fESDpileup_slope(ap.fESDpileup_slope), fESDpileup_inter(ap.fESDpileup_inter),
  fSaveESDpileupQA(ap.fSaveESDpileupQA),
  fAddTPCpileupCuts(ap.fAddTPCpileupCuts),
  fSaveTPCpileupQA(ap.fSaveTPCpileupQA),
  fControlProfileList(ap.fControlProfileList)
{
  AliInfo("----DEBUG AliJCatalystTask COPY ----");
}

//_____________________________________________________________________________
AliJCatalystTask& AliJCatalystTask::operator = (const AliJCatalystTask& ap)
{
  // assignment operator
  AliInfo("----DEBUG AliJCatalystTask operator= ----");
  this->~AliJCatalystTask();
  new(this) AliJCatalystTask(ap);
  return *this;
}

//______________________________________________________________________________
AliJCatalystTask::~AliJCatalystTask()
{
  delete fInputList;
  delete fInputListALICE;
  if (fMainList) {delete fMainList;}
}

//________________________________________________________________________
void AliJCatalystTask::UserCreateOutputObjects()
{
  fInputList = new TClonesArray("AliJBaseTrack" , 2500);
  fInputList->SetOwner(kTRUE);
  fInputListALICE = new TClonesArray("AliJBaseTrack" , 2500);
  fInputListALICE->SetOwner(kTRUE);

  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  OpenFile(1);
  
  fJCorMapTask = (AliJCorrectionMapTask*)man->GetTask(fJCorMapTaskName);
  if(!fJCorMapTask ) AliInfo("----CHECK if AliJCorrectionMapTask Missing ----");
  if( fJCorMapTask ) {
    fCentBinEff = fJCorMapTask->GetCentBinEff();
    fCentBinEff->Print();
  }

  gRandom->SetSeed();

  BookControlHistograms();
  PostData(1, fMainList);

  //OpenFile(1);
  //fOutput = gDirectory;
  //fOutput->cd();
  //PostData(1, fOutput);
}

//______________________________________________________________________________
void AliJCatalystTask::UserExec(Option_t* /*option*/)
{
  // Processing of one event
  if(!((Entry()-1)%100))  AliInfo(Form(" Processing event # %lld",  Entry()));

  // initializing variables from last event
  fJCatalystEntry = fEntry;
  fInputList->Clear();
  fInputListALICE->Clear();

  float fImpactParameter = .0; // setting 0 for the generator which doesn't have this info. 
  double fvertex[3];

  fEvtNum++;
  if(fEvtNum % 1000 == 0)
    cout << "evt : " << fEvtNum << ", fEntry: " << fEntry <<  endl;

  // load current event and save track, event info
  if(flags & FLUC_KINEONLY) {
    if(flags & FLUC_EXCLUDE_EPOS){
        Error("UserExec","FLUC_EXCLUDE_EPOS flag is not set up to work with FLUC_KINEONLY. No EPOS filtering is done.");
    }
    AliMCEvent *mcEvent;
    if(flags & FLUC_KINEONLYEXT) {
      AliInputEventHandler*  fMcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
      mcEvent = fMcHandler->MCEvent();
    }
    else {
      mcEvent = MCEvent();
    }

    if (!mcEvent) {
      AliError("ERROR: mcEvent not available");
      return;
    }

    if(!fnoCentBin) {
      AliGenHijingEventHeader* hijingHeader = dynamic_cast<AliGenHijingEventHeader*>(mcEvent->GenEventHeader());
      AliGenDPMjetEventHeader* dpmHeader = dynamic_cast<AliGenDPMjetEventHeader*>(mcEvent->GenEventHeader());
      AliGenHepMCEventHeader* hepHeader = dynamic_cast<AliGenHepMCEventHeader*>(mcEvent->GenEventHeader());
      if (hijingHeader) {
        fImpactParameter = hijingHeader->ImpactParameter();
      }
      else if (dpmHeader) {
        fImpactParameter = dpmHeader->ImpactParameter();
      }
      else if (hepHeader) {
        fImpactParameter = hepHeader->impact_parameter();
      }
      else {
        DEBUG( 4,  "KineOnly no header in event generator" );                   
      }

      fcent = GetCentralityFromImpactPar(fImpactParameter);
    }

    if(flags & FLUC_ALICE_IPINFO){
      //force to use ALICE impact parameter setting
      double ALICE_Cent[8] = {0, 5, 10, 20, 30, 40, 50, 60};
      double ALICE_IPinfo[8] = {0, 3.50, 4.94, 6.98, 8.55, 9.88, 11.04, 12.09};
      for(int icent=0; icent<8; icent++){
        if(fImpactParameter >= ALICE_IPinfo[icent] && fImpactParameter < ALICE_IPinfo[icent+1])
          fcent = 0.5f*(ALICE_Cent[icent]+ALICE_Cent[icent+1]);
      }
    }

    if(fnoCentBin) fcent = 1.0; // forcing no centrality selection
    if(flags & FLUC_KINEONLYEXT) {
      ReadKineTracks( mcEvent->Stack(), fInputList, fInputListALICE, fcent ) ; // read tracklist
    }
    else {
      ReadKineTracks( mcEvent, fInputList, fInputListALICE, fcent ) ; // read tracklist
    }

    AliGenEventHeader *header = mcEvent->GenEventHeader();
    if(!header)
      return;
    TArrayF gVertexArray;
    header->PrimaryVertex(gVertexArray);
    for(int i = 0; i < 3; i++)
      fvertex[i] = gVertexArray.At(i);
  }

  else { // Kine
    paodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
    if(fnoCentBin) {
      fcent = 1.0;
    }
    else {
      fcent = ReadCentrality(paodEvent,fCentDetName);
    }
    fRunNum = paodEvent->GetRunNumber();

    // Get the centrality bin for the current centrality event.
    Int_t centBin = GetCentralityBin(fcent);
    if (centBin == -1) {return;}

    //ReadAODTracks depends on the vertex info for phi and efficiency weight
    ReadVertexInfo(paodEvent, fvertex); // Read vertex info before selection.
    if (bSaveAllQA) {
      fVertexXHistogram[centBin][0]->Fill(fvertex[0]);
      fVertexYHistogram[centBin][0]->Fill(fvertex[1]);
      fVertexZHistogram[centBin][0]->Fill(fvertex[2]);
    }

    // Apply the event quality selection.
    /// Draw the control histogram sbefore the event selection.
    if (bSaveAllQA) {FillEventQA(paodEvent, centBin, 0);}

    /// Apply the AliEventCut (only for Run2 data).
    AliJRunTable *fRunTable = & AliJRunTable::GetSpecialInstance();
    fRunTable->SetRunNumber(fRunNum);
    int fperiod = fRunTable->GetRunNumberToPeriod(fRunNum);
    if (fperiod == AliJRunTable::kLHC15o || fperiod == AliJRunTable::kLHC18q || fperiod == AliJRunTable::kLHC18r) {
      fAliEventCuts = new AliEventCuts();
      if (fCentDetName == "CL0") {  // CL0 systematics: override the order of the centrality estimators without changing the centrality framework or min/max.
        fAliEventCuts->OverrideCentralityFramework(1);
        fAliEventCuts->SetCentralityEstimators("CL0", "V0M");
        fAliEventCuts->SetCentralityRange(0., 90.);
      }

      //fAliEventCuts->SetMaxVertexZposition(fzvtxCut);   // Zvtx set here to 10cm, recut to 8cm afterwards.
      if (fAddTPCpileupCuts) {fAliEventCuts->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE);}
      fAliEventCuts->AcceptEvent(paodEvent);
    }

    /// Apply the custom event selection.
    fIsGoodEvent = IsGoodEvent(paodEvent, centBin);
    if(!fIsGoodEvent) {
      return;
    }

    if (bSaveAllQA) { // Save the vertex and centrality if PV passes the selection.
      fVertexXHistogram[centBin][1]->Fill(fvertex[0]);
      fVertexYHistogram[centBin][1]->Fill(fvertex[1]);
      fVertexZHistogram[centBin][1]->Fill(fvertex[2]);
      fCentralityHistogram[centBin]->Fill(fcent);

      FillEventQA(paodEvent, centBin, 1);
    }

    // Load correction maps in the event loop
    //if( fEvtNum == 1 && fJCorMapTask ) {

    if(fJCorMapTask) {
      grEffCor = fJCorMapTask->GetEffCorrectionMap(fRunNum,fcent,fEffFilterBit); 
      int fcBin = 
        AliJFFlucAnalysis::GetBin(fcent,AliJFFlucAnalysis::BINNING_CENT_PbPb);
      pPhiWeights = fJCorMapTask->GetCorrectionMap(phiMapIndex,fRunNum,fcBin);
    }

    ReadAODTracks( paodEvent, fInputList, fcent ) ; // read tracklist
  } // AOD
  fZvert = fvertex[2];

}

//______________________________________________________________________________
void AliJCatalystTask::ReadVertexInfo ( AliAODEvent *aod, double*  fvtx )
{
  fvtx[0] = aod->GetPrimaryVertex()->GetX();
  fvtx[1] = aod->GetPrimaryVertex()->GetY();
  fvtx[2] = aod->GetPrimaryVertex()->GetZ();
}

//______________________________________________________________________________
float AliJCatalystTask::ReadCentrality( AliAODEvent *aod, TString Trig ){
  AliMultSelection *pms = (AliMultSelection*)aod->FindListObject("MultSelection");
  return pms?pms->GetMultiplicityPercentile(Trig.Data()):aod->GetCentrality()->GetCentralityPercentile(Trig.Data());
}

//______________________________________________________________________________

void AliJCatalystTask::Init()
{
  AliInfo("Doing initialization") ;
}
//______________________________________________________________________________
void AliJCatalystTask::ReadAODTracks(AliAODEvent *aod, TClonesArray *TrackList, float fcent)
{
  //aod->Print();
  if(flags & FLUC_MC) {  // how to get a flag to check  MC or not !
    TClonesArray *mcArray;
    if(flags & FLUC_EXCLUDE_EPOS){
        mcArray = (TClonesArray*)aod->FindListObject("FilteredParticlesTrueMC");
        if(!mcArray) {
            Error("ReadAODTracks","Could not find FilteredParticlesTrueMC, check that you are using AliAnalysisTaskSVtaskMCFilter properly or disable the FLUC_EXCLUDE_EPOS flag");
            return;
        }
    } else {
        mcArray = (TClonesArray*) aod->FindListObject(AliAODMCParticle::StdBranchName());
        if(!mcArray){ Printf("Error not a proper MC event"); };  // check mc array
    }

    Int_t nt = mcArray->GetEntriesFast();
    Int_t ntrack = 0;
    for( int it=0; it < nt ; it++){
      AliAODMCParticle *track = (AliAODMCParticle*)mcArray->At(it);
      if(!track) {
        Error("ReadEventAODMC","Could not read particle %d",it);
        continue;
      }

      if( track->IsPhysicalPrimary() ){
        // insert AMTP weak decay switch here
        if(flags & FLUC_EXCLUDEWDECAY){
          //cout << "finding weak decaying particle ... " << endl;
          Int_t gMotherIndex = track->GetMother(); // check and ask about this to DJ changed to mother from firstmother
          if(gMotherIndex != -1) {
            //cout << "this mother is " << gMotherIndex << endl;
            AliAODMCParticle *motherParticle = (AliAODMCParticle*)mcArray->At(gMotherIndex);
            //cout << "mother pdg code is " << motherParticle->GetPdgCode() << endl;
            if(motherParticle) {
              if(IsThisAWeakDecayingParticle(motherParticle)){
                //Exclude the decay products of weakly decaying particles
                continue;
              }
            }
          }
        } // weak decay particles are excluded.

        if(fPt_min > 0){
          double Pt = track->Pt();
          if( Pt < fPt_min || Pt > fPt_max )
            continue ; // pt cut
        }

        Int_t pdg = track->GetPdgCode();
        Char_t ch = (Char_t) track->Charge();
        if(ch < 0){
          if(fPcharge == 1)
            continue;
        }
        else if(ch > 0){
          if(fPcharge == -1)
            continue;
        }
        else continue;
        
        if(track->Eta() < fEta_min || track->Eta() > fEta_max) continue; // Need to check this here also
        AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
        itrack->SetLabel(track->GetLabel());
        itrack->SetParticleType( pdg);
        itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
        itrack->SetCharge(ch) ;
      }
    } // for( int it=0; it < nt ; it++){
  } // end: if(flags & FLUC_MC)

  else {  // Read AOD reco tracks.
    TClonesArray *aodArray;
    Int_t nt;
    if(flags & FLUC_EXCLUDE_EPOS){
        aodArray = (TClonesArray*)aod->FindListObject("FilteredTracksDetMC");
        if(aodArray) {
            nt = aodArray->GetEntriesFast();
            if (bSaveAllQA) {fMultHistogram[GetCentralityBin(fcent)][0]->Fill(nt);}
        } else {
            Error("ReadAODTracks","Could not find FilteredTracksDetMC, check that you are using AliAnalysisTaskSVtaskMCFilter properly or disable the FLUC_EXCLUDE_EPOS flag");
            return;
        }
    } else {
        nt = aod->GetNumberOfTracks();
        if (bSaveAllQA) {fMultHistogram[GetCentralityBin(fcent)][0]->Fill(nt);}
    }

    Int_t ntrack =0;
    Double_t PV[3] = {0.};
    ReadVertexInfo(aod, PV);

    for( int it=0; it<nt ; it++){
      AliAODTrack *track;
      if(flags & FLUC_EXCLUDE_EPOS) track = dynamic_cast<AliAODTrack*>(aodArray->At(it));
      else track = dynamic_cast<AliAODTrack*>(aod->GetTrack(it));
      if(!track) {
        Error("ReadEventAOD", "Could not read particle %d", (int) it);
        continue;
      }

      if (bSaveAllQA) {FillControlHistograms(track, 0, fcent, PV);} // Fill the QA histograms before the track selection.

      if(track->GetTPCNcls() < fNumTPCClusters)
        continue;

      if(fUseITSMinClusters && (track->GetITSNcls()<fITSMinClusters))
        continue;

      // New: Get the values of the DCA according to the type of tracks.
      Float_t DCAxy = 0.; // DCA in the transverse plane.
      Float_t DCAz = 0.;    // DCA along the beam axis. 

        track->GetImpactParameters(DCAxy,DCAz);

      // Newer: Set the tighter cuts for PbPb Run2: primary cuts along with hybrids.
      if (fUseTightCuts){
        // Set the golden cut on chi2: < 36
        if ((track->GetChi2TPCConstrainedVsGlobal()) > 36) {continue;}
        // Redefine the max for DCAxy as a function of pT.
        fDCAxy_max = 0.0105 + 0.0350/(TMath::Power(track->Pt(), 1.1));
      }
      if (bUseDCAbaseCut) {
        bDCABaseCut = TMath::Abs(DCAxy)<0.0208+0.04/TMath::Power(track->Pt(),1.01);
      }

      // New: Apply the cuts on the DCA values of the track.
      if(TMath::Abs(DCAxy) > fDCAxy_max) {continue;}
      if(TMath::Abs(DCAz) > fDCAz_max) {continue;}
      if((bUseDCAbaseCut) && (!bDCABaseCut)) {continue;}

      // New: Apply the cut on the chi2 per ndf for the TPC tracks.
      Double_t chi2NDF = track->Chi2perNDF(); // TBC: is this 100% the right one?
      if((chi2NDF < fChi2perNDF_min) || (chi2NDF > fChi2perNDF_max)) {continue;}

      if(track->TestFilterBit( fFilterBit )){ //
        if( fPt_min > 0){
          double Pt = track->Pt();
          if( Pt < fPt_min || Pt > fPt_max )
            continue ; // pt cut
        }

        Char_t ch = (Char_t)track->Charge();
        if(ch < 0){
          if(fPcharge == 1)
            continue;
        }
        else if(ch > 0){
          if(fPcharge == -1)
            continue;
        }
        else continue;

        if(track->Eta() < fEta_min || track->Eta() > fEta_max) continue; // Need to check this here also
                // Removal of bad area, now only with eta symmetric
        Bool_t isBadArea = TMath::Abs(track->Eta()) > 0.6;
        Bool_t isBadArea18q = track->Eta() < -(fZvert+8.)/17.5;
        if(fremovebadarea && isBadArea) continue;

        if (fremovebadarea18q && isBadArea18q) continue;
        
        if (bSaveAllQA) {FillControlHistograms(track, 1, fcent, PV);} // Fill the QA histograms after the track selection.

        AliJBaseTrack *itrack = new( (*TrackList)[ntrack++]) AliJBaseTrack;
        itrack->SetID( TrackList->GetEntriesFast() );
        itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
        itrack->SetParticleType(kJHadron);
        itrack->SetCharge(track->Charge() );
        itrack->SetStatus(track->GetStatus() );
        // Apply eff correction !!!
        Double_t effCorr = 1.0;
        if(grEffCor && fJCorMapTask) {
          effCorr = fJCorMapTask->GetEffCorrection(grEffCor,track->Pt());//fEfficiency->GetCorrection(pt,fEffFilterBit,fCent);
        } 
        itrack->SetTrackEff(effCorr);
        
        // Adding phi weight for a track
        AliJRunTable *fRunTable = & AliJRunTable::GetSpecialInstance();
        fRunTable->SetRunNumber(fRunNum);
        int fperiod = fRunTable->GetRunNumberToPeriod(fRunNum); // Needed for the alternative corrections.

        Double_t phi_module_corr = 1.0;
        if(bUseAlternativeWeights && fperiod == AliJRunTable::kLHC10h){ // Alternative NUA for 10h
           Int_t RunBin = GetRunIndex10h((Int_t)aod->GetRunNumber());
           Int_t CentralityBin = GetCentralityBin(fcent);
           Int_t phiBin = fHistoPhiWeight[CentralityBin][RunBin]->FindBin(track->Phi());
           phi_module_corr = 1./((double)(fHistoPhiWeight[CentralityBin][RunBin])->GetBinContent(phiBin));
        } else {
          if(fJCorMapTask){
            Double_t w;
            if(pPhiWeights) {
              Double_t phi = itrack->Phi();
              Double_t eta = itrack->Eta();
              w = pPhiWeights->GetBinContent(pPhiWeights->FindBin(phi,eta,fZvert));
            }
            else {
              w = 1.0;
            }
            if(w > 1e-6) phi_module_corr = w;
          }
        }
        itrack->SetWeight(phi_module_corr);

        // Uncommentable for local test
        // if(phi_module_corr<0.5){printf("Track info: %.4f \t %.4f \t %.4f \t %.4f \t %d \n", itrack->Phi(), itrack->Eta(), fZvert, phi_module_corr, fRunNum);}

        // Adding centrality weight for a LHC15o track, and given centrality.
        float cent_weight = 1.0;
        if (bUseAlternativeWeights && fperiod == AliJRunTable::kLHC15o) {
          Int_t RunBin = GetRunIndex15o((Int_t)aod->GetRunNumber());
          //printf("RunBin = %d and fcent = %.2f\n", RunBin, fcent);
          if (!fHistoCentWeight[RunBin]) {printf("We don't have a valid histo for cent for this run.\n");}
          Int_t centBin = fHistoCentWeight[RunBin]->FindBin(fcent);
          //printf("centBin = %d\n", centBin);
          cent_weight = (float)(fHistoCentWeight[RunBin]->GetBinContent(centBin));
        }
        itrack->SetCentWeight(cent_weight);

        if (bSaveAllQA){
          fProfileWeights[GetCentralityBin(fcent)]->Fill(itrack->Phi(), phi_module_corr);
          fPhiHistogram[GetCentralityBin(fcent)][2]->Fill(track->Phi(), 1./phi_module_corr);
          fPTHistogram[GetCentralityBin(fcent)][2]->Fill(itrack->Pt(), 1./effCorr);
        }

        if (bSaveQCNUA){
          
          f2DEtaPhiHistogram[GetCentralityBin(fcent)][0]->Fill(itrack->Phi(),itrack->Eta());                    //Before NUA correction
          f2DEtaPhiHistogram[GetCentralityBin(fcent)][1]->Fill(itrack->Phi(),itrack->Eta(),1./phi_module_corr); //After NUA correction
          //cout << "I'm in the QCNUA boolean" << endl;
          for( int ih=0; ih<7; ih++){ //loop over 2-8th harmonic
            //cout << "Entering harmonic loop: " << ih << ", with centrality: " << fcent << ", and Cos: " << TMath::Cos(ih+2*itrack->Phi()) << endl;
            Float_t profcent = GetCentralityBin(fcent)+0.5 ;
            //cout <<"fcent: " << fcent << "centrality bin: " << profcent << endl;
            Double_t arg_phi = (ih+2)*itrack->Phi();
            fProfileCosVSCent[ih]->Fill(fcent,TMath::Cos(arg_phi),1./phi_module_corr);
            fProfileSinVSCent[ih]->Fill(fcent,TMath::Sin(arg_phi),1./phi_module_corr);
          }
        }
      } // End: if(track->TestFilterBit( fFilterBit ))
    }

  } //read aod reco track done.
  if(fDebugLevel>1) cout << "Tracks: " << TrackList->GetEntriesFast() << endl;
  if (bSaveAllQA) {fMultHistogram[GetCentralityBin(fcent)][1]->Fill(TrackList->GetEntriesFast());}
  if (bSaveQCNUA){
     for (int ih=0; ih<7; ih++){
      fControlProfileList->Add(fProfileCosVSCent[ih]);
      fControlProfileList->Add(fProfileSinVSCent[ih]);
     }
   } 
}
//______________________________________________________________________________
Bool_t AliJCatalystTask::IsGoodEvent( AliAODEvent *event, Int_t thisCent){
  //event selection here!
  //check vertex

  AliVVertex *vtx = event->GetPrimaryVertex();
  if(!vtx || vtx->GetNContributors() <= 0)
    return kFALSE;
  double zvert = vtx->GetZ();
  if(zvert < -fzvtxCut || zvert > fzvtxCut)
    return kFALSE;

  if(flags & FLUC_CENT_FLATTENING){
    float fCent = ReadCentrality(event,fCentDetName);
    TH1 *pweightMap = fJCorMapTask->GetCentCorrection();
    if(gRandom->Uniform(0,1) > pweightMap->GetBinContent(pweightMap->GetXaxis()->FindBin(fCent)))
      return kFALSE;
  }

  if(flags & FLUC_KINEONLY)
    return kTRUE;

  //int frunNumber = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetEvent()->GetRunNumber();
  AliJRunTable *fRunTable = & AliJRunTable::GetSpecialInstance();
  fRunTable->SetRunNumber(fRunNum);

  int fperiod = fRunTable->GetRunNumberToPeriod(fRunNum);
  if(fperiod == AliJRunTable::kLHC15o || fperiod == AliJRunTable::kLHC18q || fperiod == AliJRunTable::kLHC18r){
    const AliVVertex* vtTrc = event->GetPrimaryVertex();
    const AliVVertex* vtSPD = event->GetPrimaryVertexSPD();
    double covTrc[6],covSPD[6];
    vtTrc->GetCovarianceMatrix(covTrc);
    vtSPD->GetCovarianceMatrix(covSPD);
    double dz = vtTrc->GetZ()-vtSPD->GetZ();
    double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
    double errTrc = TMath::Sqrt(covTrc[5]);
    double nsigTot = TMath::Abs(dz)/errTot, nsigTrc = TMath::Abs(dz)/errTrc;
    if(TMath::Abs(dz) > 0.2 || nsigTot > 10 || nsigTrc > 20)
      return kFALSE;

    AliMultSelection *pms = (AliMultSelection*)event->FindListObject("MultSelection");
    if(!pms){
      AliError("MultSelection unavailable.");
      return kFALSE;
    }

    double v0mcent = pms->GetMultiplicityPercentile("V0M");
    double cl0cent = pms->GetMultiplicityPercentile("CL0");
    double center = 0.973488*cl0cent+0.0157497;
    double sigma = 0.673612+cl0cent*(0.0290718+cl0cent*(-0.000546728+cl0cent*5.82749e-06));
    if(v0mcent < center-5.0*sigma || v0mcent > center+5.5*sigma || v0mcent < 0.0 || v0mcent > 60.0)
      return kFALSE;
  }

  TPCTracks = 0;
  GlobTracks = 0;
  Int_t nTracks = event->GetNumberOfTracks();
  for(int it = 0; it < nTracks; it++){
    AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(it));
    //AliAODTrack* trackAOD = event->GetTrack(itracks);
    if (!trackAOD || !trackAOD->TestFilterBit(1))
      continue;
    if ((trackAOD->Pt() < 0.2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > 0.8) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.2) )
      continue;
    TPCTracks++;
  }

  for(int it = 0; it < nTracks; it++){
    AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(it));
    if (!trackAOD || !trackAOD->TestFilterBit(16))
      continue;
    if ((trackAOD->Pt() < 0.2) || (trackAOD->Pt() > 5.0) || (TMath::Abs(trackAOD->Eta()) > 0.8) || (trackAOD->GetTPCNcls() < 70) || (trackAOD->GetDetPid()->GetTPCsignal() < 10.0) || (trackAOD->Chi2perNDF() < 0.1) )
      continue;
    Double_t b[2];
    Double_t bCov[3];
    if (!(trackAOD->PropagateToDCA(event->GetPrimaryVertex(), event->GetMagneticField(), 100., b, bCov) ))
      continue;
    if ( (TMath::Abs(b[0]) > 0.3) || (TMath::Abs(b[1]) > 0.3) )
      continue;
    GlobTracks++;
  }

  FB32Tracks = 0;
  FB32TOFTracks = 0;
  for (int it = 0; it < nTracks; it++){
    AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(it));
    if (!trackAOD || !trackAOD->TestFilterBit(32))
      continue;
    FB32Tracks++;
    if (TMath::Abs(trackAOD->GetTOFsignalDz()) <= 10 && trackAOD->GetTOFsignal() >= 12000 && trackAOD->GetTOFsignal() <= 25000)
      FB32TOFTracks++;
  }

  // Define the pileup cuts for EDStracks vs TPConlyTracks.
  UInt_t M_ESD = 0; // Multiplicity corresponding ESD tracks.
  UInt_t M_TPC = 0; // Multiplicity TPC only tracks
  if (fAddESDpileupCuts) {
    M_ESD = ((AliAODHeader*)event->GetHeader())->GetNumberOfESDTracks();
    for (int it = 0; it < nTracks; it++) {
      AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(it));
      if (!trackAOD) {continue;}
      if (trackAOD->TestFilterBit(128)) {M_TPC++;}
    }
  }


  if(flags & FLUC_CUT_OUTLIERS){
    if(fperiod == AliJRunTable::kLHC15o || fperiod == AliJRunTable::kLHC18q || fperiod == AliJRunTable::kLHC18r){
      AliMultSelection *pms = (AliMultSelection*)event->FindListObject("MultSelection");
      if(!pms){
        AliError("MultSelection unavailable.");
        return kFALSE;
      }

      double v0mcent = pms->GetMultiplicityPercentile("V0M");
      double tfbtpc = (double)TPCTracks;
      double lcut = 2.31181837e+03+v0mcent*(-7.79946952e+01+v0mcent*(8.45194500e-01+v0mcent*(-1.72787009e-03-1.86192490e-05*v0mcent)));
      if(tfbtpc < lcut)
        return kFALSE;
      double hcut = 3.15901050e+03+v0mcent*(-9.42636072e+01+v0mcent*(8.06432447e-01+v0mcent*(3.37574557e-03-6.14272547e-05*v0mcent)));
      if(tfbtpc > hcut)
        return kFALSE;

      double tfb32 = (double)FB32Tracks;
      double tfb32tof = (double)FB32TOFTracks;
      double mu32tof = -1.0178+tfb32*(0.333132+tfb32*(9.10282e-05-1.61861e-08*tfb32));
      double sigma32tof = 1.47848+tfb32*(0.0385923+tfb32*(-5.06153e-05+tfb32*(4.37641e-08+tfb32*(-1.69082e-11+tfb32*2.35085e-15))));
      double nsigma[] = {4.0,4.0};
      if(tfb32tof < mu32tof-nsigma[0]*sigma32tof || tfb32tof > mu32tof+nsigma[1]*sigma32tof)
        return kFALSE;

      // if true, apply the pileup cut between ESD and TPConly tracks.
      if (fAddESDpileupCuts) {
        //if (bSaveAllQA && fSaveESDpileupQA) {fESDpileupHistogram[thisCent][0]->Fill(M_TPC, M_ESD);}
        if (!( (double)M_ESD < (fESDpileup_slope*(double)M_TPC + fESDpileup_inter) ))  {return kFALSE;}
        //if (bSaveAllQA && fSaveESDpileupQA) {fESDpileupHistogram[thisCent][1]->Fill(M_TPC, M_ESD);}
      }
    }
    else if (fperiod == AliJRunTable::kLHC10h) {  // High multiplicity outlier cuts for LHC10h based on the SCklm analysis.
      UInt_t MTPC = 0;
      UInt_t Mglobal = 0;
      for(int it = 0; it < nTracks; it++){
        AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(it));
        if (!trackAOD) {continue;}
        if (trackAOD->TestFilterBit(128)) {MTPC++;}
        if (trackAOD->TestFilterBit(256)) {Mglobal++;}
      }

      if (bSaveAllQA && bSaveHMOhist) {fHMOsHistogram[thisCent][0]->Fill(Mglobal, MTPC);} // Fill the HMO histogram only if both bSave* are true.
      if(!((double)MTPC > (-65.0 + 1.54*Mglobal) && (double)MTPC < (90.0 + 2.30*Mglobal))) {
        return kFALSE;
      }
      if (fCentDetName == "V0M") {
        float cent_V0M = ReadCentrality(event,"V0M");
        if ((cent_V0M > 20.) && (cent_V0M < 30.) && (MTPC > 1500)) {return kFALSE;}
        if ((cent_V0M > 30.) && (cent_V0M < 40.) && (MTPC > 1050)) {return kFALSE;}
        if ((cent_V0M > 40.) && (cent_V0M < 50.) && (MTPC > 700)) {return kFALSE;}
      }
      if (bSaveAllQA && bSaveHMOhist) {fHMOsHistogram[thisCent][1]->Fill(Mglobal, MTPC);}
    }
    else {  // TBA: QA histo if chosen, possibility to choose between the two versions for LHC10h
      if(!((double)TPCTracks > (-40.3+1.22*GlobTracks) && (double)TPCTracks < (32.1+1.59*GlobTracks)))
        return kFALSE;
    }
      
  }

  return kTRUE;
}
//
//______________________________________________________________________________
void AliJCatalystTask::Terminate(Option_t *)
{
  //fFFlucAna->Terminate("");
  cout<<"AliJCatalystTask Analysis DONE !!"<<endl;
}


//______________________________________________________________________________
Bool_t AliJCatalystTask::IsThisAWeakDecayingParticle(AliMCParticle *thisGuy)
{
  // In order to prevent analyzing daughters from weak decays
  // - AMPT does not only strong decays, so IsPhysicalPrimary does not catch it
  Int_t pdgcode = TMath::Abs( thisGuy->PdgCode() );
  Int_t myWeakParticles[7] = { 3322, 3312, 3222, // Xi0 Xi+- Sigma-+
    3122, 3112, // Lambda0 Sigma+-
    130, 310 // K_L0 K_S0
  };
  for(Int_t i=0; i!=7; ++i)
    if( myWeakParticles[i] == pdgcode ){
      return kTRUE;
    }

  return kFALSE;
}
//===============================================================================
Bool_t AliJCatalystTask::IsThisAWeakDecayingParticle(AliAODMCParticle *thisGuy)
{
  // In order to prevent analyzing daughters from weak decays
  // - AMPT does not only strong decays, so IsPhysicalPrimary does not catch it
  Int_t pdgcode = TMath::Abs( thisGuy->GetPdgCode() );
  Int_t myWeakParticles[7] = { 3322, 3312, 3222, // Xi0 Xi+- Sigma-+
    3122, 3112, // Lambda0 Sigma+-
    130, 310 // K_L0 K_S0
  };
  for(Int_t i=0; i!=7; ++i)
    if( myWeakParticles[i] == pdgcode ) {
      return kTRUE;
    }
  return kFALSE;
}
//______________________________________________________________________________
void AliJCatalystTask::SetEffConfig( int effMode, int FilterBit)
{
        fEffMode = effMode;
        switch(FilterBit){
        case 128:
                fEffFilterBit = 0;
                break;
        case 96:
                fEffFilterBit = 4;
                break;
        case 768:
                fEffFilterBit = 5;
                break;
        default:
                fEffFilterBit = 0;
                break;
        }
        cout << "setting to EffCorr Mode : " << effMode << endl;
        cout << "setting to EffCorr Filter bit : " << FilterBit  << " = " << fEffFilterBit << endl;
}
//______________________________________________________________________________
void AliJCatalystTask::ReadKineTracks( AliMCEvent *mcEvent, TClonesArray *TrackList, TClonesArray *TrackListALICE, float fcent)
{
  Int_t nt = mcEvent->GetNumberOfPrimaries();
  Int_t ntrack = 0;
  for (Int_t it = 0; it < nt; it++) {
    AliMCParticle* track = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(it));
    if(mcEvent->IsPhysicalPrimary(it)) {
      double Pt = track->Pt();
      if( Pt < fPt_min || Pt > fPt_max )
        continue ; // pt cut
      TParticle *particle = track->Particle();
      if(!particle)
        continue;

      if(flags & FLUC_EXCLUDEWDECAY){
        Int_t gMotherIndex = particle->GetFirstMother(); //
        if(gMotherIndex != -1){
          DEBUG( 4,  "this particle has a mother " );
          AliMCParticle* motherParticle= dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(gMotherIndex));
          if(motherParticle){
            if(IsThisAWeakDecayingParticle( motherParticle)){
              DEBUG ( 4, Form("this particle will be removed because it comes from : %d", motherParticle->PdgCode() ));
              continue;
            }
          }
        }
      }

      Int_t pdg = particle->GetPdgCode();
      Char_t ch = (Char_t) track->Charge();
      if(ch < 0){
        if(fPcharge == 1)
          continue;
      }else
      if(ch > 0){
        if(fPcharge == -1)
          continue;
      }else continue;
      if(track->Eta() < fEta_min || track->Eta() > fEta_max) continue;
      AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
      Int_t label = track->GetLabel();
      itrack->SetLabel( label );
      itrack->SetParticleType(pdg);
      itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
      itrack->SetCharge(ch) ;
      if(TMath::Abs(track->Eta()) < 0.8) {
        AliJBaseTrack *jtrack =  new ((*TrackListALICE)[TrackListALICE->GetEntriesFast()])AliJBaseTrack;
        jtrack->SetLabel( label );
        jtrack->SetParticleType(pdg);
        jtrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->E() );
        jtrack->SetCharge(ch) ;
      }
    }
  }
  if(fDebugLevel>1) cout << "Tracks: " << TrackList->GetEntriesFast() << endl;
}
// To read the track generated from a external alievent generators
void AliJCatalystTask::ReadKineTracks( AliStack *stack, TClonesArray *TrackList, TClonesArray *TrackListALICE, float fcent)
{
  Int_t nt = stack->GetNprimary();
  Int_t ntrack = 0;
  for (Int_t it = 0; it < nt; it++) {
    TParticle* track = stack->Particle(it); if(!track) continue;
    if (!AliPWG0Helper::IsPrimaryCharged(track, nt)) continue; 
      double Pt = track->Pt();
      if( Pt < fPt_min || Pt > fPt_max )
        continue ; // pt cut
      Int_t pdg = track->GetPdgCode();
      Char_t ch = (Char_t) TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
      if(ch < 0){
        if(fPcharge == 1)
          continue;
      }else
      if(ch > 0){
        if(fPcharge == -1)
          continue;
      }else continue;
      if(track->Eta() < fEta_min || track->Eta() > fEta_max) continue;
      AliJBaseTrack *itrack = new ((*TrackList)[ntrack++])AliJBaseTrack;
      Int_t label = 100;//track->GetLabel();
      itrack->SetLabel( label );
      itrack->SetParticleType(pdg);
      itrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->Energy() );
      itrack->SetCharge(ch) ;
      if(TMath::Abs(track->Eta()) < 0.8) {
        AliJBaseTrack *jtrack =  new ((*TrackListALICE)[TrackListALICE->GetEntriesFast()])AliJBaseTrack;
        jtrack->SetLabel( label );
        jtrack->SetParticleType(pdg);
        jtrack->SetPxPyPzE( track->Px(), track->Py(), track->Pz(), track->Energy() );
        jtrack->SetCharge(ch) ;
      }
  }
  if(fDebugLevel>1) cout << "Tracks: " << TrackList->GetEntriesFast() << endl;
}

double AliJCatalystTask::GetCentralityFromImpactPar(double ip) {
  //https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies
  static double bmin[12] = {0.0,1.60,2.27,3.72,5.23,7.31,8.88,10.20,11.38,12.47,14.51,100};
  static double centmean[12] = {0.5,1.5,3.5,7.5,15,25,35,45,55,65,75,90};
  for(UInt_t i = 0; i < 11; i++){
    if(bmin[i+1] > ip)
      return centmean[i];
  }
  return 0.0;
}

//______________________________________________________________________________
void AliJCatalystTask::InitializeArrays() {
  for(Int_t icent=0; icent<fCentralityBins; icent++){

    fControlHistogramsList[icent] = NULL;
    for(Int_t i=0; i<2; i++) {
      fVertexXHistogram[icent][i] = NULL;
      fVertexYHistogram[icent][i] = NULL;
      fVertexZHistogram[icent][i] = NULL;
      fTPCClustersHistogram[icent][i] = NULL;
      fITSClustersHistogram[icent][i] = NULL;
      fChiSquareTPCHistogram[icent][i] = NULL;
      fChiSquareITSHistogram[icent][i] = NULL;      
      fDCAzHistogram[icent][i] = NULL;
      fDCAxyHistogram[icent][i] = NULL;
      fChargeHistogram[icent][i] = NULL;
      fPTHistogram[icent][i] = NULL;
      fPhiHistogram[icent][i] = NULL;
      fEtaHistogram[icent][i] = NULL;
      fMultHistogram[icent][i] = NULL;
      fHMOsHistogram[icent][i] = NULL;
      fESDpileupHistogram[icent][i] = NULL;
      fTPCpileupHistogram[icent][i] = NULL;
      f2DEtaPhiHistogram[icent][i] = NULL;
    }
      fPTHistogram[icent][2] = NULL;
      fPhiHistogram[icent][2] = NULL;

    fCentralityHistogram[icent] = NULL;
    fProfileWeights[icent] = NULL;
    for(int iRun = 0; iRun < 90; iRun++) {fHistoPhiWeight[icent][iRun] = NULL;}
  }
  for (int iR = 0; iR < 138; iR++) {fHistoCentWeight[iR] = NULL;}
  for (int ih=0; ih<7; ih++){    //loop over the harmonics 2-8
    fProfileCosVSCent[ih] = NULL;
    fProfileSinVSCent[ih] = NULL;
  }
}

//______________________________________________________________________________
void AliJCatalystTask::BookControlHistograms(){

fControlProfileList = new TList();
fControlProfileList->SetName("ControlProfiles");
fControlProfileList->SetOwner(kTRUE);
if(bSaveQCNUA){ fMainList->Add(fControlProfileList); }

for(Int_t icent=0; icent<fCentralityBins; icent++) //loop over all centrality bins
{

  //Check if value of this centrality bin is negative -> if yes: break. We do not need anymore
  if(fcentralityArray[icent+1] < 0)
  {
     break; //The next edge is a breaking point -> this bin does not exist anymore
  }

  fControlHistogramsList[icent] = new TList();
  fControlHistogramsList[icent]->SetName(Form("ControlHistograms_%.1f-%.1f", fcentralityArray[icent], fcentralityArray[icent+1]));
  fControlHistogramsList[icent]->SetOwner(kTRUE);
  if(bSaveAllQA){ fMainList->Add(fControlHistogramsList[icent]); }

  fControl2DNUAList[icent] = new TList();
  fControl2DNUAList[icent]->SetName(Form("Control2DNUAHistograms_%.1f-%.1f", fcentralityArray[icent], fcentralityArray[icent+1]));
  fControl2DNUAList[icent]->SetOwner(kTRUE);
  if(bSaveQCNUA){ fMainList->Add(fControl2DNUAList[icent]); }


   // a) Book histogram to hold pt spectra:
   fPTHistogram[icent][0] = new TH1F("fPTHist_BeforeTrackSelection","Pt Distribution",1000,0.,10.);
   fPTHistogram[icent][0]->GetXaxis()->SetTitle("P_t");
   fPTHistogram[icent][0]->SetLineColor(4);
   fControlHistogramsList[icent]->Add(fPTHistogram[icent][0]);

   fPTHistogram[icent][1] = new TH1F("fPTHist_AfterTrackSelection","Pt Distribution",1000,0.,10.);
   fPTHistogram[icent][1]->GetXaxis()->SetTitle("P_t");
   fPTHistogram[icent][1]->SetLineColor(4);
   fControlHistogramsList[icent]->Add(fPTHistogram[icent][1]);

   fPTHistogram[icent][2] = new TH1F("fPTHist_AfterTrackSelection_Weighted","Pt Distribution",1000,0.,10.);
   fPTHistogram[icent][2]->GetXaxis()->SetTitle("P_t");
   fPTHistogram[icent][2]->SetLineColor(4);
   fControlHistogramsList[icent]->Add(fPTHistogram[icent][2]);
   
   // b) Book histogram to hold phi spectra
   fPhiHistogram[icent][0] = new TH1F("fPhiHist_BeforeTrackSelection","Phi Distribution",1000,0.,TMath::TwoPi()); 
   fPhiHistogram[icent][0]->GetXaxis()->SetTitle("Phi");
   fPhiHistogram[icent][0]->SetLineColor(4);
   fControlHistogramsList[icent]->Add(fPhiHistogram[icent][0]);

   fPhiHistogram[icent][1] = new TH1F("fPhiHist_AfterTrackSelection","Phi Distribution",1000,0.,TMath::TwoPi()); 
   fPhiHistogram[icent][1]->GetXaxis()->SetTitle("Phi");
   fPhiHistogram[icent][1]->SetLineColor(4);
   fControlHistogramsList[icent]->Add(fPhiHistogram[icent][1]);

   fPhiHistogram[icent][2] = new TH1F("fPhiHist_AfterTrackSelection_Weighted","Phi Distribution",1000,0.,TMath::TwoPi()); 
   fPhiHistogram[icent][2]->GetXaxis()->SetTitle("Phi");
   fPhiHistogram[icent][2]->SetLineColor(4);
   fControlHistogramsList[icent]->Add(fPhiHistogram[icent][2]);

   // c) Book histogram to hold eta distribution before track selection:
   fEtaHistogram[icent][0] = new TH1F("fEtaHist_BeforeTrackSelection","Eta Distribution",1000,-1.,1.); 
   fEtaHistogram[icent][0]->GetXaxis()->SetTitle("Eta");
   fEtaHistogram[icent][0]->SetLineColor(4);
   fControlHistogramsList[icent]->Add(fEtaHistogram[icent][0]);

   fEtaHistogram[icent][1] = new TH1F("fEtaHist_AfterTrackSelection","Eta Distribution",1000,-1.,1.);
   fEtaHistogram[icent][1]->GetXaxis()->SetTitle("Eta");
   fEtaHistogram[icent][1]->SetLineColor(4);
   fControlHistogramsList[icent]->Add(fEtaHistogram[icent][1]);

   // d) Book histogam to hold multiplicty distributions 
   fMultHistogram[icent][0] = new TH1F("fMultiHisto_BeforeTrackSelection","Multiplicity",30000,0.,30000.); 
   fMultHistogram[icent][0]->GetXaxis()->SetTitle("Multiplicity M");
   fControlHistogramsList[icent]->Add(fMultHistogram[icent][0]);
   
   fMultHistogram[icent][1] = new TH1F("fMultiHisto_AfterTrackSelection","Multiplicity",30000,0.,30000.); 
   fMultHistogram[icent][1]->GetXaxis()->SetTitle("Multiplicity M");
   fControlHistogramsList[icent]->Add(fMultHistogram[icent][1]);

   // e) Book histogam for Vertex X 
   fVertexXHistogram[icent][0] = new TH1F("fVertexX_BeforeEventSelection","VertexXBefore",2000,-10.,10.); 
   fVertexXHistogram[icent][0]->GetXaxis()->SetTitle("");
   fControlHistogramsList[icent]->Add(fVertexXHistogram[icent][0]);

   fVertexXHistogram[icent][1] = new TH1F("fVertexX_AfterEventSelection","VertexXAfter",2000,-10.,10.); 
   fVertexXHistogram[icent][1]->GetXaxis()->SetTitle("");
   fControlHistogramsList[icent]->Add(fVertexXHistogram[icent][1]);

   // f) Book histogam for Vertex Y 
   fVertexYHistogram[icent][0] = new TH1F("fVertexY_BeforeEventSelection","VertexYBefore",2000,-10.,10.); 
   fVertexYHistogram[icent][0]->GetXaxis()->SetTitle("");
   fControlHistogramsList[icent]->Add(fVertexYHistogram[icent][0]);

   fVertexYHistogram[icent][1] = new TH1F("fVertexY_AfterEventSelection","VertexYAfter",2000,-10.,10.); 
   fVertexYHistogram[icent][1]->GetXaxis()->SetTitle("");
   fControlHistogramsList[icent]->Add(fVertexYHistogram[icent][1]);

   // g) Book histogam for Vertex Z 
   fVertexZHistogram[icent][0] = new TH1F("fVertexZ_BeforeEventSelection","VertexZBefore",4000,-20.,20.); 
   fVertexZHistogram[icent][0]->GetXaxis()->SetTitle("");
   fControlHistogramsList[icent]->Add(fVertexZHistogram[icent][0]);

   fVertexZHistogram[icent][1] = new TH1F("fVertexZ_AfterEventSelection","VertexZAfter",4000,-20.,20.); 
   fVertexZHistogram[icent][1]->GetXaxis()->SetTitle("");
   fControlHistogramsList[icent]->Add(fVertexZHistogram[icent][1]);

   // i) Book histogram for number of TPC clustes 
   fTPCClustersHistogram[icent][0] = new TH1F("fTPCClusters_BeforeCut","TPCClustersBeforeCut",170,0.,170.); 
   fControlHistogramsList[icent]->Add(fTPCClustersHistogram[icent][0]);

   fTPCClustersHistogram[icent][1] = new TH1F("fTPCClusters_AfterCut","TPCClustersAfterCut",170,0.,170.); 
   fControlHistogramsList[icent]->Add(fTPCClustersHistogram[icent][1]);

   //j) Book histogram for number of ITC clusters 
   fITSClustersHistogram[icent][0] = new TH1F("fITSClusters_BeforeCut","ITSClustersBeforeCut",10,0.,10.); 
   fControlHistogramsList[icent]->Add(fITSClustersHistogram[icent][0]);

   fITSClustersHistogram[icent][1] = new TH1F("fITSClusters_AfterCut","ITSClustersAfterCut",10,0.,10.); 
   fControlHistogramsList[icent]->Add(fITSClustersHistogram[icent][1]);

   // k) Book histogram for chi square TPC and ITS
   fChiSquareTPCHistogram[icent][0] = new TH1F("fChiSquareTPC_BeforeCut","ChiSquareTPCBeforeCut",1000,0.,20.); 
   fControlHistogramsList[icent]->Add(fChiSquareTPCHistogram[icent][0]);

   fChiSquareTPCHistogram[icent][1] = new TH1F("fChiSquareTPC_AfterCut","ChiSquareTPCAfterCut",1000,0.,20.); 
   fControlHistogramsList[icent]->Add(fChiSquareTPCHistogram[icent][1]);

   fChiSquareITSHistogram[icent][0] = new TH1F("fChiSquareITS_BeforeCut","ChiSquareITSBeforeCut",1000,0.,50.); 
   fControlHistogramsList[icent]->Add(fChiSquareITSHistogram[icent][0]);

   fChiSquareITSHistogram[icent][1] = new TH1F("fChiSquareITS_AfterCut","ChiSquareITSAfterCut",1000,0.,50.); 
   fControlHistogramsList[icent]->Add(fChiSquareITSHistogram[icent][1]);

    // l) Book histogram for DCAz
   fDCAzHistogram[icent][0] = new TH1F("fDCAz_BeforeCut","DCAzBeforeCut",1000,-10.,10.);  
   fControlHistogramsList[icent]->Add(fDCAzHistogram[icent][0]);

   fDCAzHistogram[icent][1] = new TH1F("fDCAz_AfterCut","DCAzAfterCut",1000,-10.,10.); 
   fControlHistogramsList[icent]->Add(fDCAzHistogram[icent][1]);
   
   // m) Book histogram for DCAxy
   fDCAxyHistogram[icent][0] = new TH1F("fDCAxy_BeforeCut","DCAxyBeforeCut",1000,-10.,10.); 
   fControlHistogramsList[icent]->Add(fDCAxyHistogram[icent][0]);

   fDCAxyHistogram[icent][1] = new TH1F("fDCAxy_AfterCut","DCAxyAfterCut",1000,-10.,10.); 
   fControlHistogramsList[icent]->Add(fDCAxyHistogram[icent][1]); 

   // n) Book histogram for Charge
   fChargeHistogram[icent][0] = new TH1I("fCharge_BeforeCut","ChargeBeforeCut",11,-5.5,5.5); 
   fControlHistogramsList[icent]->Add(fChargeHistogram[icent][0]);

   fChargeHistogram[icent][1] = new TH1I("fCharge_AfterCut","ChargeAfterCut",11,-5.5,5.5); 
   fControlHistogramsList[icent]->Add(fChargeHistogram[icent][1]);

   // o) Book histogram Centrality 
   fCentralityHistogram[icent]= new TH1F("fCentralityHistogram_After","CentralityHistogramAfter",22,0.,110.);
   fCentralityHistogram[icent]->GetXaxis()->SetTitle("Centrality");
   fCentralityHistogram[icent]->SetLineColor(4);
   fControlHistogramsList[icent]->Add(fCentralityHistogram[icent]);

   // p) Book the TH2D for the HMOs in LHC10h.
   fHMOsHistogram[icent][0] = new TH2D("fHMOsHistogram_Before","Correlations before HMO cuts", 1000, 0.,5000., 1000, 0., 5000.);
   fHMOsHistogram[icent][0]->GetXaxis()->SetTitle("M_{global}");
   fHMOsHistogram[icent][0]->GetYaxis()->SetTitle("M_{TPC}");
   if (bSaveHMOhist) {fControlHistogramsList[icent]->Add(fHMOsHistogram[icent][0]);}

   fHMOsHistogram[icent][1] = new TH2D("fHMOsHistogram_After","Correlations after HMO cuts", 1000, 0.,5000., 1000, 0., 5000.);
   fHMOsHistogram[icent][1]->GetXaxis()->SetTitle("M_{global}");
   fHMOsHistogram[icent][1]->GetYaxis()->SetTitle("M_{TPC}");
   if (bSaveHMOhist) {fControlHistogramsList[icent]->Add(fHMOsHistogram[icent][1]);}

   // q) Book histogram to hold 2D eta-phi spectra
   f2DEtaPhiHistogram[icent][0] = new TH2F("f2DEtaPhiHist_BeforeCorrection","eta-phi; #phi (rad); #eta",50,0.,TMath::TwoPi(),16,-0.8,0.8); 
   f2DEtaPhiHistogram[icent][0]->GetXaxis()->SetTitle("#phi");
   f2DEtaPhiHistogram[icent][0]->GetYaxis()->SetTitle("#eta");
   if (bSaveQCNUA) {fControl2DNUAList[icent]->Add(f2DEtaPhiHistogram[icent][0]);}

   f2DEtaPhiHistogram[icent][1] = new TH2F("f2DEtaPhiHist_AfterCorrection","eta-phi; #phi (rad); #eta",50,0.,TMath::TwoPi(),16,-0.8,0.8); 
   f2DEtaPhiHistogram[icent][1]->GetXaxis()->SetTitle("#phi");
   f2DEtaPhiHistogram[icent][1]->GetYaxis()->SetTitle("#eta");
   if (bSaveQCNUA) {fControl2DNUAList[icent]->Add(f2DEtaPhiHistogram[icent][0]);}


   // TProfile for the weights to apply.
   fProfileWeights[icent] = new TProfile("fProfileWeights","Phi Weights",1000,-TMath::Pi(),TMath::Pi()); //centrality dependent output
   fProfileWeights[icent]->GetXaxis()->SetTitle("#varphi");
   fProfileWeights[icent]->GetYaxis()->SetTitle("weight");
   fControlHistogramsList[icent]->Add(fProfileWeights[icent]);

   //Control TProfiles requested by pag;  <cos(nphi)>  and <sin(nphi)>  vs centrality
   fProfileCosVSCent[0] = new TProfile("fProfileCosVSCent_n2","<cos(2phi)> vs Cent",12,0,60); 
   fProfileCosVSCent[1] = new TProfile("fProfileCosVSCent_n3","<cos(3phi)> vs Cent",12,0,60); 
   fProfileCosVSCent[2] = new TProfile("fProfileCosVSCent_n4","<cos(4phi)> vs Cent",12,0,60); 
   fProfileCosVSCent[3] = new TProfile("fProfileCosVSCent_n5","<cos(5phi)> vs Cent",12,0,60); 
   fProfileCosVSCent[4] = new TProfile("fProfileCosVSCent_n6","<cos(6phi)> vs Cent",12,0,60); 
   fProfileCosVSCent[5] = new TProfile("fProfileCosVSCent_n7","<cos(7phi)> vs Cent",12,0,60); 
   fProfileCosVSCent[6] = new TProfile("fProfileCosVSCent_n8","<cos(8phi)> vs Cent",12,0,60); 

   fProfileSinVSCent[0] = new TProfile("fProfileSinVSCent_n2","<sin(2phi)> vs Cent",12,0,60); 
   fProfileSinVSCent[1] = new TProfile("fProfileSinVSCent_n3","<sin(3phi)> vs Cent",12,0,60); 
   fProfileSinVSCent[2] = new TProfile("fProfileSinVSCent_n4","<sin(4phi)> vs Cent",12,0,60); 
   fProfileSinVSCent[3] = new TProfile("fProfileSinVSCent_n5","<sin(5phi)> vs Cent",12,0,60); 
   fProfileSinVSCent[4] = new TProfile("fProfileSinVSCent_n6","<sin(6phi)> vs Cent",12,0,60); 
   fProfileSinVSCent[5] = new TProfile("fProfileSinVSCent_n7","<sin(7phi)> vs Cent",12,0,60); 
   fProfileSinVSCent[6] = new TProfile("fProfileSinVSCent_n8","<sin(8phi)> vs Cent",12,0,60);
   
   



   // Save the TH2D for the ESD-TPC pileup cut in Run2.
    fESDpileupHistogram[icent][0] = new TH2D("fESDpileupHistogram_Before","Correlations before ESD-TPC cut", 1200, 0., 6000., 7000, 0., 35000.);
    fESDpileupHistogram[icent][0]->GetXaxis()->SetTitle("M_{TPC}");
    fESDpileupHistogram[icent][0]->GetYaxis()->SetTitle("M_{ESD}");
    if (fSaveESDpileupQA) {fControlHistogramsList[icent]->Add(fESDpileupHistogram[icent][0]);}

    fESDpileupHistogram[icent][1] = new TH2D("fESDpileupHistogram_After","Correlations after ESD-TPC cuts", 1200, 0., 6000., 7000, 0., 35000.);
    fESDpileupHistogram[icent][1]->GetXaxis()->SetTitle("M_{TPC}");
    fESDpileupHistogram[icent][1]->GetYaxis()->SetTitle("M_{ESD}");
    if (fSaveESDpileupQA) {fControlHistogramsList[icent]->Add(fESDpileupHistogram[icent][1]);}

    // Save the TH2D for the ITS-TPC cluster pileup cut in Run2.
    fTPCpileupHistogram[icent][0] = new TH2I("fTPCpileupHistogram_Before","Correlations before ITS-TPC clusters cut", 1e3, 1e5, 5e6, 1e3, 1e8, 5e9);
    fTPCpileupHistogram[icent][0]->GetXaxis()->SetTitle("N_{clusters} in TPC");
    fTPCpileupHistogram[icent][0]->GetYaxis()->SetTitle("N_{clusters} in SDD-SSD");
    if (fSaveTPCpileupQA) {fControlHistogramsList[icent]->Add(fTPCpileupHistogram[icent][0]);}

    fTPCpileupHistogram[icent][1] = new TH2I("fTPCpileupHistogram_After","Correlations after ITS-TPC clusters cut", 1e3, 1e5, 5e6, 1e3, 1e8, 5e9);
    fTPCpileupHistogram[icent][1]->GetXaxis()->SetTitle("N_{clusters} in TPC");
    fTPCpileupHistogram[icent][1]->GetYaxis()->SetTitle("N_{clusters} in SDD-SSD");
    if (fSaveTPCpileupQA) {fControlHistogramsList[icent]->Add(fTPCpileupHistogram[icent][1]);}

  }//for(Int_t icent=0; icent<fCentralityBins; icent++)
}

//______________________________________________________________________________
void AliJCatalystTask::FillControlHistograms(AliAODTrack *thisTrack, Int_t whichHisto, Float_t cent, Double_t *v) {

// Get the corresponding centrality bin.
  Int_t CentralityBin = GetCentralityBin(cent);

  Float_t ValueDCAxy = 999.;   // DCA in the xy-plane.
  Float_t ValueDCAz = 999.;    // DCA along z.

  thisTrack->GetImpactParameters(ValueDCAxy,ValueDCAz);

  fPTHistogram[CentralityBin][whichHisto]->Fill(thisTrack->Pt());
  fPhiHistogram[CentralityBin][whichHisto]->Fill(thisTrack->Phi());
  fEtaHistogram[CentralityBin][whichHisto]->Fill(thisTrack->Eta());
  fTPCClustersHistogram[CentralityBin][whichHisto]->Fill(thisTrack->GetTPCNcls());
  fITSClustersHistogram[CentralityBin][whichHisto]->Fill(thisTrack->GetITSNcls());
  fChiSquareTPCHistogram[CentralityBin][whichHisto]->Fill(thisTrack->Chi2perNDF());
  
  if (!thisTrack->GetITSNcls()==0) { // If check needed to avoid floating point exception for some cases.
    fChiSquareITSHistogram[CentralityBin][whichHisto]->Fill(thisTrack->GetITSchi2()/(float)thisTrack->GetITSNcls());
  }
  fDCAzHistogram[CentralityBin][whichHisto]->Fill(ValueDCAz);
  fDCAxyHistogram[CentralityBin][whichHisto]->Fill(ValueDCAxy);
  fChargeHistogram[CentralityBin][whichHisto]->Fill(thisTrack->Charge());
}

//______________________________________________________________________________
Int_t AliJCatalystTask::GetCentralityBin(Float_t cent)
{
  //if this functions returns a negative value -> Error. No Centrality could be selected
  
  //Check for centrality bin
  for(Int_t icent=0; icent<fCentralityBins+1; icent++) //loop over all centrality bins
  {
    if(fcentralityArray[icent]<0) {return -1;}
    if(cent >= fcentralityArray[icent]) { continue; } 
    else { return icent-1; } 
  } 

  //We went through all centrality edges without returning. This means, that the measured value is bigger than the maximum centrality that we want for our analyis
  return -1;

}

//______________________________________________________________________________
void AliJCatalystTask::SetInitializeCentralityArray()
{

 TString sMethodName = "void AliJCatalystTask::SetInitializeCentralityArray()";

  Float_t ListCentralities[17] = { fcent_0, fcent_1, fcent_2, fcent_3, fcent_4, fcent_5, fcent_6, fcent_7, fcent_8, fcent_9, fcent_10, fcent_11, fcent_12, fcent_13, fcent_14, fcent_15, fcent_16 };

  for(Int_t i=0; i<17; i++) { fcentralityArray[i] = ListCentralities[i]; }

  //Protections
  if(fcentralityArray[0] < 0 || fcentralityArray[1] < 0) { Fatal(sMethodName.Data(),"First Centrality bin not defined"); } //They need at least one well defined centrality bin

  for(Int_t icent=0; icent<fCentralityBins; icent++)
  {
    //The next bin should be a valid boundery, i.e. it is > 0. but it is also smaller than the previous boundery -> Wrong ordering
    if( fcentralityArray[icent+1] > 0. && fcentralityArray[icent+1] < fcentralityArray[icent] ) { Fatal(sMethodName.Data(),"Wrong ordering of centrality bounderies"); }
  }
}

//______________________________________________________________________________
//______________________________________________________________________________
//______________________________________________________________________________

Int_t AliJCatalystTask::GetRunIndex10h(Int_t runNumber)
{
// Return for the given run the index in the run-by-run arrays.
  TString sMethod = "Int_t AliJCatalystTask::GetRunIndex10h()";
  Int_t cRun = -1; // Current index in the loop.

  const int NumberRuns = 90;
  Int_t listRuns[90] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};

// Find the position of the given run into the list of runs.
  for (Int_t iRun = 0; iRun < NumberRuns; iRun++)
  {
    if (listRuns[iRun] == runNumber)
    {
      cRun = iRun;
      break;
    } // End: for (Int_t iRun = 0; iRun < fNumberRuns; iRun++).
  } // End: iRun.

 if(cRun == -1){Fatal(sMethod.Data(), "FATAL: Run Number not in List of Runs!");}  

  return cRun;
} // End: Int_t GetRunIndex(Int_t).

//==========================================================================================================================================================================

void AliJCatalystTask::SetInputAlternativeNUAWeights10h(bool UseAltWeight, TString fileWeight)
{
// Setter to open the external file with the particle weights and import them in the task....
// a.)  Open the external file.                                                                  
// b.)  Parse the runs.                                                                          
// b.1) Open the TDirectoryFile for the current run.                                            
// b.2) Open the list for the current centrality range.                                                                                         
// c.)  Close the external file.                                                               

  TString sMethod = "void AliJCatalystTask::SetInputAlternativeNUAWeights10h()"; 

  bUseAlternativeWeights = UseAltWeight;

  const int NumberRuns = 90;
  Int_t listRuns[90] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};

  // a.) Open the external file.
  TFile *weightsFile = TFile::Open(Form("%s", fileWeight.Data()), "READ");
  if (!weightsFile) {Fatal(sMethod.Data(), "ERROR 404: File not found");}

  for (Int_t iRun = 0; iRun < NumberRuns; iRun++)
  {
    // b.1) Open the TDirectoryFile for the current run.
    Int_t runNumber = listRuns[iRun];

    TDirectoryFile *runTDF = dynamic_cast<TDirectoryFile*>(weightsFile->Get(Form("%d", runNumber)));
    if (!runTDF) {Fatal(sMethod.Data(), "ERROR: Directory not found");}

    for(Int_t icent=0; icent<fCentralityBins; icent++) //loop over all centrality bins //GANESHA: Check about fCentralityBins. Hardcode???
    {
        //Check if value of this centrality bin is negativ -> if yes: break. We do not need anymore
        //The next edge is a breaking point -> this bin does not exist anymore
        if(fcentralityArray[icent+1] < 0) {break;}

        // b.2) Open the list for the current centrality range.
        TList *centralityList = dynamic_cast<TList*>(runTDF->Get(Form("Centrality-%.1f-%.1f", fcentralityArray[icent], fcentralityArray[icent+1])));
        if (!centralityList) {Fatal(sMethod.Data(), "ERROR: List not found");}

        fHistoPhiWeight[icent][iRun] = dynamic_cast<TH1F*>(centralityList->FindObject("phi-weight"));
        if (!fHistoPhiWeight[icent][iRun]) {Fatal(sMethod.Data(), "ERROR: phi-weight histogram not found");}
        else {fHistoPhiWeight[icent][iRun]->SetDirectory(0);}  // Kill the default ownership

        delete centralityList;
     } // End: icent
     delete runTDF;
  }//for irun

  // c.) Close the external file.
  weightsFile->Close();
  delete weightsFile;

} // void AliAnalysisTaskStudentsML::SetInputParticleWeights(TString fileWeight)

//==========================================================================================================================================================================

void AliJCatalystTask::FillEventQA(AliAODEvent *event, int centBin, int stepBin)
{
// Fill the control histograms for the event selection.
// stepBin = 0: before the selection, 1: after the selection.

// ESD-TPC 2d histograms.
  UInt_t M_ESD = 0; // Multiplicity corresponding ESD tracks.
  UInt_t M_TPC = 0; // Multiplicity TPC only tracks
  Int_t nTracks = event->GetNumberOfTracks();
  if (fAddESDpileupCuts) {
    M_ESD = ((AliAODHeader*)event->GetHeader())->GetNumberOfESDTracks();

    for (int it = 0; it < nTracks; it++) {
      AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(event->GetTrack(it));
      if (!trackAOD) {continue;}
      if (trackAOD->TestFilterBit(128)) {M_TPC++;}
    }

  if (fSaveESDpileupQA) {fESDpileupHistogram[centBin][stepBin]->Fill(M_TPC, M_ESD);}
  }

// TPC pileup
  if (fAddTPCpileupCuts) {
    int nTPCclusters = event->GetNumberOfTPCClusters();

    AliVMultiplicity *multi = event->GetMultiplicity();
    int nSDclusters = 0;  // Number of clusters in the 4 external layers of the ITS (SDD-SSD).
    for (int iLay = 2; iLay <=6; iLay++) {nSDclusters += multi->GetNumberOfITSClusters(iLay);}
    //printf("nITSclusters: %d nTPCclusters: %d \n", nSDclusters, nTPCclusters);

    if (fSaveTPCpileupQA) {fTPCpileupHistogram[centBin][stepBin]->Fill(nTPCclusters, nSDclusters);}
  }



} // oid AliJCatalystTask::FillEventQA(AliAODEvent *event, int stepBin)

//==========================================================================================================================================================================

Int_t AliJCatalystTask::GetRunIndex15o(Int_t runNumber)
{
// Return for the given run the index in the run-by-run arrays.
  TString sMethod = "Int_t AliJCatalystTask::GetRunIndex15o()";
  Int_t cRun = -1; // Current index in the loop.

  const int NumberRuns = 138;
  Int_t listRuns[NumberRuns] = {244917, 245347, 245507, 245829, 246113, 246424, 246846, 244918,
      245349, 245535, 245831, 246115, 246431, 246847, 244975, 245353, 245540, 245833, 246148,
      246434, 246851, 244980, 245396, 245542, 245923, 246151, 246750, 246858, 244982, 245397,
      245543, 245949, 246152, 246751, 246859, 244983, 245401, 245544, 245952, 246153, 246757,
      246864, 245064, 245407, 245545, 245954, 246178, 246758, 246865, 245066, 245409, 245554,
      245963, 246180, 246759, 246867, 245068, 245410, 245683, 246001, 246181, 246760, 246870,
      245145, 245411, 245692, 246003, 246182, 246763, 246871, 245146, 245441, 245702, 246012,
      246185, 246765, 246928, 245151, 245446, 245705, 246036, 246217, 246766, 246945, 245152,
      245450, 245729, 246037, 246222, 246804, 246948, 245231, 245453, 245731, 246042, 246225,
      246805, 246982, 245232, 245454, 245752, 246048, 246271, 246807, 246984, 245233, 245496,
      245759, 246049, 246272, 246808, 246989, 245259, 245497, 245766, 246052, 246275, 246809,
      246991, 245343, 245501, 245775, 246053, 246276, 246810, 246994, 245345, 245504, 245785,
      246087, 246391, 246844, 245346, 245505, 245793, 246089, 246392, 246845};

// Find the position of the given run into the list of runs.
  for (Int_t iRun = 0; iRun < NumberRuns; iRun++)
  {
    if (listRuns[iRun] == runNumber)
    {
      cRun = iRun;
      break;
    } // End: for (Int_t iRun = 0; iRun < fNumberRuns; iRun++).
  } // End: iRun.

 if(cRun == -1){Fatal(sMethod.Data(), "FATAL: Run Number not in List of Runs!");}  

  return cRun;
}

//==========================================================================================================================================================================
void AliJCatalystTask::SetInputCentralityWeight15o(bool useCentWeight, TString fileCentWeight)
{
// Get the centrality correction for LHC15o.
  TString sMethod = "void AliJCatalystTask::SetInputCentralityWeight15o()";

  // Enable the use of non-unit centrality correction. Recycled from alternative NUA for 10h.
  bUseAlternativeWeights = useCentWeight;

  // Declare the list of runs for LHC15o, pass2, AOD252.
  const int NumberRuns = 138;
  Int_t listRuns[NumberRuns] = {244917, 245347, 245507, 245829, 246113, 246424, 246846, 244918,
      245349, 245535, 245831, 246115, 246431, 246847, 244975, 245353, 245540, 245833, 246148,
      246434, 246851, 244980, 245396, 245542, 245923, 246151, 246750, 246858, 244982, 245397,
      245543, 245949, 246152, 246751, 246859, 244983, 245401, 245544, 245952, 246153, 246757,
      246864, 245064, 245407, 245545, 245954, 246178, 246758, 246865, 245066, 245409, 245554,
      245963, 246180, 246759, 246867, 245068, 245410, 245683, 246001, 246181, 246760, 246870,
      245145, 245411, 245692, 246003, 246182, 246763, 246871, 245146, 245441, 245702, 246012,
      246185, 246765, 246928, 245151, 245446, 245705, 246036, 246217, 246766, 246945, 245152,
      245450, 245729, 246037, 246222, 246804, 246948, 245231, 245453, 245731, 246042, 246225,
      246805, 246982, 245232, 245454, 245752, 246048, 246271, 246807, 246984, 245233, 245496,
      245759, 246049, 246272, 246808, 246989, 245259, 245497, 245766, 246052, 246275, 246809,
      246991, 245343, 245501, 245775, 246053, 246276, 246810, 246994, 245345, 245504, 245785,
      246087, 246391, 246844, 245346, 245505, 245793, 246089, 246392, 246845};

  // Open the external file.
  TFile *weightsFile = TFile::Open(Form("%s", fileCentWeight.Data()), "READ");
  if (!weightsFile) {Fatal(sMethod.Data(), "ERROR 404: File not found");}

  for (Int_t iRun = 0; iRun < NumberRuns; iRun++)
  {
    // Get the histogram for the current run number.
    Int_t runNumber = listRuns[iRun];
    fHistoCentWeight[iRun] = static_cast<TH1F*>(weightsFile->Get(
      Form("histoCentWeight_%d", runNumber) ));
    if (!fHistoCentWeight[iRun]) {Fatal(sMethod.Data(), "ERROR: cent weight histogram not found");}
    else {fHistoCentWeight[iRun]->SetDirectory(0);}  // Kill the default ownership
  } // for irun

  // c.) Close the external file.
  weightsFile->Close();
  delete weightsFile;
}
