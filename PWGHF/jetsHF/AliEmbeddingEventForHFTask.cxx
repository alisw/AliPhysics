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
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
//Embedding task modified to check D-meson candidates in the base sample
//Antônio Silva (antonio.silva@cern.ch)

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TGrid.h>
#include <TH2C.h>
#include <TList.h>
#include <TStreamerInfo.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TRandom2.h>

// AliRoot
#include "AliVEvent.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliCentrality.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliAODVertex.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliNamedString.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliStack.h"
#include "AliEmcalParticle.h"
#include "AliParticleContainer.h"
#include "AliEmbeddingEventForHFTask.h"
#include "AliAnalysisVertexingHF.h"
#include "AliMultSelectionTask.h"
#include "AliOADBMultSelection.h"
#include "AliOADBContainer.h"

ClassImp(AliEmbeddingEventForHFTask)

//_______________________________________________________________________________
AliEmbeddingEventForHFTask::AliEmbeddingEventForHFTask() :
  AliAnalysisTaskEmcal(),
  fAttempts(5),
  fCurrentAODFileID(0),
  fCurrentAODFile(0),
  fCurrentAODEntry(-1),
  fFirstAODEntry(-1),
  fLastAODEntry(-1),
  fEmbeddingCount(0),
  fCurrentAODTree(0),
  fAODTreeName("aodTree"),
  fAODHeaderName("header"),
  fAODVertexName("vertices"),
  fAODTrackName("tracks"),
  fAODVZEROName("AliAODVZERO"),
  fGridDir("alien:///alice/data/2015/LHC15o"),
  fPassAndDataType("pass1/AOD"),
  fDataFile("AliAOD.root"),
  fMaxFileNumber(500),
  fAODHeader(0),
  fAODFilePath(0),
  fZVertexCut(10),
  fMaxVertexDist(999),
  fMinCentrality(-1),
  fMaxCentrality(-1),
  fTriggerMask(AliVEvent::kAny),
  fAODVertex(0),
  fAODTracks(0),
  fAODVZERO(0),
  fEmbeddedTracks(0),
  fEmbeddedTracksName("EmbeddedTracks"),
  fCuts(0),
  fAodEvent(0),
  fArrayDStartoD0pi(0),
  fInhibitTask(kFALSE),
  fCandidateType(0),
  fBranchName(""),
  fMCarray(0),
  fPDGmother(0),
  fNProngs(0),
  fRejectQuarkNotFound(kTRUE),
  fRejectDfromB(kTRUE),
  fKeepOnlyDfromB(kFALSE),
  fCheckDmeson(kTRUE),
  fDSignalOnly(kTRUE),
  fRandomAccess(kTRUE),
  fIsSimulation(kTRUE),
  fEmbTrackEff(kFALSE),
  fEmbTrackEffPath(""),
  fEmbTrackEffHistName(""),
  fBaseTrackEff(kFALSE),
  fBaseTrackEffPath(""),
  fBaseTrackEffHistName(""),
  fRandomGetter(0),
  fTrackEffHist(0),
  fTrackEffON(),
  fTrackEffOFF(),
  fEmbTrackEffHist(0),
  fEmbTrackEffAccepted(0),
  fEmbTrackEffAll(0),
  fBaseTrackEffHist(0),
  fBaseTrackEffAccepted(0),
  fBaseTrackEffAll(0)
  
{
  //
  // Default constructor
  //
  for (Int_t i=4; i--;) fPDGdaughters[i] = 0;

}

//_______________________________________________________________________________
AliEmbeddingEventForHFTask::AliEmbeddingEventForHFTask(const char *name, AliRDHFCuts *cuts, ACandidateType candtype) :
  AliAnalysisTaskEmcal(name, kTRUE),
  fAttempts(5),
  fCurrentAODFileID(0),
  fCurrentAODFile(0),
  fCurrentAODEntry(-1),
  fFirstAODEntry(-1),
  fLastAODEntry(-1),
  fEmbeddingCount(0),
  fCurrentAODTree(0),
  fAODTreeName("aodTree"),
  fAODHeaderName("header"),
  fAODVertexName("vertices"),
  fAODTrackName("tracks"),
  fAODVZEROName("AliAODVZERO"),
  fGridDir("alien:///alice/data/2015/LHC15o"),
  fPassAndDataType("pass1/AOD"),
  fDataFile("AliAOD.root"),
  fMaxFileNumber(500),
  fAODHeader(0),
  fAODFilePath(0),
  fZVertexCut(10),
  fMaxVertexDist(999),
  fMinCentrality(-1),
  fMaxCentrality(-1),
  fTriggerMask(AliVEvent::kAny),
  fAODVertex(0),
  fAODTracks(0),
  fAODVZERO(0),
  fEmbeddedTracks(0),
  fEmbeddedTracksName("EmbeddedTracks"),
  fCuts(cuts),
  fAodEvent(0),
  fArrayDStartoD0pi(0),
  fInhibitTask(kFALSE),
  fCandidateType(candtype),
  fBranchName(""),
  fMCarray(0),
  fPDGmother(0),
  fNProngs(0),
  fRejectQuarkNotFound(kTRUE),
  fRejectDfromB(kTRUE),
  fKeepOnlyDfromB(kFALSE),
  fCheckDmeson(kTRUE),
  fDSignalOnly(kTRUE),
  fRandomAccess(kTRUE),
  fIsSimulation(kTRUE),
  fEmbTrackEff(kFALSE),
  fEmbTrackEffPath(""),
  fEmbTrackEffHistName(""),
  fBaseTrackEff(kFALSE),
  fBaseTrackEffPath(""),
  fBaseTrackEffHistName(""),
  fRandomGetter(0),
  fTrackEffHist(0),
  fTrackEffON(),
  fTrackEffOFF(),
  fEmbTrackEffHist(0),
  fEmbTrackEffAccepted(0),
  fEmbTrackEffAll(0),
  fBaseTrackEffHist(0),
  fBaseTrackEffAccepted(0),
  fBaseTrackEffAll(0)


{
    
   //
   // Constructor. Initialization of Inputs and Outputs
   //

  Info("AliEmbeddingEventForHFTask","Calling Constructor");  
    
  for (Int_t i=4; i--;) fPDGdaughters[i] = 0;
    
  

  if(fCandidateType==0)
  {
      fBranchName = "D0toKpi";
      fPDGmother = 421;
      fNProngs = 2;
      fPDGdaughters[0] = 211;  // pi
      fPDGdaughters[1] = 321;  // K
      fPDGdaughters[2] = 0;    // empty
      fPDGdaughters[3] = 0;    // empty
  }
  else
  {
      fBranchName = "Dstar";
      fPDGmother = 413;
      fNProngs = 3;
      fPDGdaughters[1] = 211; // pi soft
      fPDGdaughters[0] = 421; // D0
      fPDGdaughters[2] = 211; // pi fromD0
      fPDGdaughters[3] = 321; // K from D0
  }



  DefineOutput(1, TList::Class());       // histos
  
}

//_______________________________________________________________________________
AliEmbeddingEventForHFTask::~AliEmbeddingEventForHFTask()
{
  //
  // Destructor
  //

  Info("~AliEmbeddingEventForHFTask","Calling Destructor");

    // Destructor
    
    if (fCurrentAODFile) {
        fCurrentAODFile->Close();
        delete fCurrentAODFile;
    }
    
    if (fCuts) {
        delete fCuts;
        fCuts   = 0;
    }

    if (fEmbeddedTracks) {
        delete fEmbeddedTracks;
        fEmbeddedTracks = 0;
    }
    
}

//_______________________________________________________________________________

void AliEmbeddingEventForHFTask::UserCreateOutputObjects()
{
  //
  // Create output objects
  //

  Info("UserCreateOutputObjects","CreateOutputObjects of task %s", GetName());

   AliAnalysisTaskEmcal::UserCreateOutputObjects();

    
   fEmbeddedTracks = new TClonesArray("AliAODTrack",50);
    
  
  fEmbeddedTracks->SetOwner();
  fEmbeddedTracks->SetName(fEmbeddedTracksName);
  
  fBaseTrackEffAccepted = new TH1F("fBaseTrackEffAccepted","fBaseTrackEffAccepted",100,0,100);
  fOutput->Add(fBaseTrackEffAccepted);
  fBaseTrackEffAll = new TH1F("fBaseTrackEffAll","fBaseTrackEffAll",100,0,100);
  fOutput->Add(fBaseTrackEffAll);
  
  fEmbTrackEffAccepted = new TH1F("fEmbTrackEffAccepted","fEmbTrackEffAccepted",100,0,100);
  fOutput->Add(fEmbTrackEffAccepted);
  fEmbTrackEffAll = new TH1F("fEmbTrackEffAll","fEmbTrackEffAll",100,0,100);
  fOutput->Add(fEmbTrackEffAll);

  PostData(1, fOutput);
  

  Info("UserCreateOutputObjects","Data posted for task %s", GetName());
}

//_______________________________________________________________________________
void AliEmbeddingEventForHFTask::ExecOnce()
{
  //
  // To be executed only once, for the first event
  //
    delete fRandomGetter;
    fRandomGetter = new TRandom3(0);
    
    if (fInhibitTask) return;
    
    // Load the event
    fAodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
    
    if (fAodEvent) {
        fArrayDStartoD0pi = dynamic_cast<TClonesArray*>(fAodEvent->GetList()->FindObject(fBranchName.Data()));
    }
    else {
        if (AODEvent() && IsStandardAOD()) {
            
            // In case there is an AOD handler writing a standard AOD, use the AOD
            // event in memory rather than the input (ESD) event.
            fAodEvent = dynamic_cast<AliAODEvent*>(AODEvent());
            
            // in this case the branches in the deltaAOD (AliAOD.VertexingHF.root)
            // have to taken from the AOD event hold by the AliAODExtension
            AliAODHandler *aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
            if(aodHandler->GetExtensions()) {
                AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
                AliAODEvent *aodFromExt = ext->GetAOD();
                fArrayDStartoD0pi = (TClonesArray*)aodFromExt->GetList()->FindObject(fBranchName.Data());
            }
            else {
                AliError(Form("This task need an AOD event! Task '%s' will be disabled!", GetName()));
                fInhibitTask = kTRUE;
                return;
            }
        }
    }
    if (fArrayDStartoD0pi) {
        TString objname(fArrayDStartoD0pi->GetClass()->GetName());
        TClass cls(objname);
        if (!cls.InheritsFrom("AliAODRecoDecayHF2Prong")) {
            AliError(Form("%s: Objects of type %s in %s are not inherited from AliAODRecoDecayHF2Prong! Task will be disabled!",
                          GetName(), cls.GetName(), fArrayDStartoD0pi->GetName()));
            fInhibitTask = kTRUE;
            fArrayDStartoD0pi = 0;
            return;
        }
    }
    else {
        AliError(Form("Could not find array %s, skipping the event. Task '%s' will be disabled!", fBranchName.Data(), GetName()));
        fInhibitTask = kTRUE;
        return;
    }

    fMCarray = dynamic_cast<TClonesArray*>(fAodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!fMCarray) {
        AliError(Form("MC particles not found! Task '%s' will be disabled!", GetName()));
        fInhibitTask = kTRUE;
        return;
    }
    
  AliDebug(2, "Entering ExecOnce()");

    fAODFilePath = static_cast<AliNamedString*>(InputEvent()->FindListObject("AODEmbeddingFile"));
    if (!fAODFilePath) {
        fAODFilePath = new AliNamedString("AODEmbeddingFile", "");
        AliDebug(3,"Adding AOD embedding file path object to the event list...");
        InputEvent()->AddObject(fAODFilePath);
    }
    
    AddObjectToEvent(fEmbeddedTracks,kFALSE);
    
    AliAnalysisTaskEmcal::ExecOnce();
  
}
Bool_t AliEmbeddingEventForHFTask::OpenNextFile()
{
    if (fCurrentAODFile) {
        fCurrentAODFile->Close();
        delete fCurrentAODFile;
        fCurrentAODFile = 0;
    }
    
    Int_t i = 0;
    
    while ((!fCurrentAODFile || fCurrentAODFile->IsZombie()) && i < fAttempts) {
        fCurrentAODFile = GetNextFile();
        i++;
    }
    
    if (!fCurrentAODFile || fCurrentAODFile->IsZombie()) {
        AliError("Could not open AOD file to embed!");
        return kFALSE;
    }
    
    fCurrentAODTree = static_cast<TTree*>(fCurrentAODFile->Get(fAODTreeName));
    if (!fCurrentAODTree) {
        AliError(Form("Could not get tree %s from file %s", fAODTreeName.Data(), fCurrentAODFile->GetName()));
        return kFALSE;
    }
    
    if (!fAODHeaderName.IsNull())
        fCurrentAODTree->SetBranchAddress(fAODHeaderName, &fAODHeader);
    
    if (!fAODVertexName.IsNull())
        fCurrentAODTree->SetBranchAddress(fAODVertexName, &fAODVertex);
    
    if (!fAODTrackName.IsNull())
        fCurrentAODTree->SetBranchAddress(fAODTrackName, &fAODTracks);
    
    if (!fAODVZEROName.IsNull())
        fCurrentAODTree->SetBranchAddress(fAODVZEROName, &fAODVZERO);

    
    
    if (fRandomAccess) {
        fFirstAODEntry = TMath::Nint(fRandomGetter->Rndm()*fCurrentAODTree->GetEntries())-1;
    }
    else {
        fFirstAODEntry = -1;
    }
    
    
    fLastAODEntry = fCurrentAODTree->GetEntries();
    fCurrentAODEntry = fFirstAODEntry;
    
    AliDebug(3,Form("Will start embedding from entry %d", fCurrentAODEntry+1));
    
    
    fEmbeddingCount = 0;
    
    return kTRUE;
}

TFile* AliEmbeddingEventForHFTask::GetNextFile()
{
    
    TString fileName;
    TString runNumber;
    
    if(fInputEvent)
    {   
        AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
        if(fIsSimulation)
        {
            runNumber = Form("%d", aodEvent->GetRunNumber());
        }
        else
        {
            runNumber = Form("000%d", aodEvent->GetRunNumber());
        }
        
        TString dirNumberString;
        Int_t dirNumber = TMath::Nint(fRandomGetter->Rndm()*fMaxFileNumber);
        
        if(dirNumber < 10) dirNumberString = Form("000%d",dirNumber);
        else if(dirNumber < 100) dirNumberString = Form("00%d",dirNumber);
        else if(dirNumber < 1000) dirNumberString = Form("0%d",dirNumber);
        else dirNumberString = Form("%d",dirNumber);
        
        fileName = Form("%s/%s/%s/%s/%s",fGridDir.Data(),runNumber.Data(),fPassAndDataType.Data(),dirNumberString.Data(),fDataFile.Data());
    }
    
    
    
    if (fileName.BeginsWith("alien://") && !gGrid) {
        AliInfo("Trying to connect to AliEn ...");
        TGrid::Connect("alien://");
    }
    
    TString baseFileName(fileName);
    if (baseFileName.Contains(".zip#")) {
        Ssiz_t pos = baseFileName.Last('#');
        baseFileName.Remove(pos);
    }
    
    if (gSystem->AccessPathName(baseFileName)) {
        AliError(Form("File %s does not exist!", baseFileName.Data()));
        return 0;
    }
    
    AliDebug(3,Form("Trying to open file %s...", fileName.Data()));
    TFile *file = TFile::Open(fileName);
    
    if (!file || file->IsZombie()) {
        AliError(Form("Unable to open file: %s!", fileName.Data()));
        return 0;
    }
    
    return file;
}

Bool_t AliEmbeddingEventForHFTask::GetNextEntry()
{
    Int_t attempts = -1;
    
    do {
        if (fCurrentAODEntry+1 >= fLastAODEntry) // in case it did not start from the first entry, it will go back
        {
            fLastAODEntry = fFirstAODEntry;
            fFirstAODEntry = -1;
            fCurrentAODEntry = -1;
        }
        
        
        
        if (!fCurrentAODFile || !fCurrentAODTree || fCurrentAODEntry+1 >= fLastAODEntry)
        {
            if (!OpenNextFile())
            {
                AliError("Could not open the next file!");
                return kFALSE;
            }
        }
        
        if (!fCurrentAODTree) {
            AliError("Could not get the tree!");
            return kFALSE;
        }
        
        fCurrentAODEntry++;
        fCurrentAODTree->GetEntry(fCurrentAODEntry);
        
        attempts++;
        if (attempts == 1000)
            AliWarning("After 1000 attempts no event has been accepted by the event selection (trigger, centrality...)!");
        
    } while (!IsAODEventSelected());
    
    
    if (!fCurrentAODTree)
        return kFALSE;
    
    if (fAODHeader)
    {
        AliAODHeader *aodHeader = static_cast<AliAODHeader*>(fAODHeader);
        AliAODHeader *evHeader = static_cast<AliAODHeader*>(InputEvent()->GetHeader());
    }
    
    fEmbeddingCount++;
    
    return kTRUE;
}
Bool_t AliEmbeddingEventForHFTask::IsAODEventSelected()
{
    // AOD event selection.
    
    if (fAODHeader)
    {
        AliAODHeader *aodHeader = static_cast<AliAODHeader*>(fAODHeader);
        
        // Trigger selection
        if (fTriggerMask != 0)
        {
            UInt_t offlineTrigger = aodHeader->GetOfflineTrigger();
            if ((offlineTrigger & fTriggerMask) == 0)
            {
                AliDebug(2,Form("Event rejected due to physics selection. Event trigger mask: %d, trigger mask selection: %d.",
                                offlineTrigger, fTriggerMask));
                return kFALSE;
            }
        }
        
        // Centrality selection
        if (fMinCentrality >= 0)
        {
            Float_t centVal = -999;
            
            
            if(aodHeader->GetRunNumber()>200000)
            {
                Float_t fAmplitude_V0A=fAODVZERO->GetMTotV0A();
                Float_t fAmplitude_V0C=fAODVZERO->GetMTotV0C();
                
                Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
                ((AliVVertex*)fAODVertex->At(0))->GetXYZ(lBestPrimaryVtxPos);
                
                Float_t fEvSel_VtxZ = lBestPrimaryVtxPos[2];
                
                Float_t NewCent = ((fAmplitude_V0A)+(fAmplitude_V0C))/(((fEvSel_VtxZ)<-14.5)*((-1.5616500000)+(-0.3058090000)*((fEvSel_VtxZ)-(1.6873000000))+(-0.0094082000)*TMath::Power((fEvSel_VtxZ)-(1.6873000000),2)) + (TMath::Abs((fEvSel_VtxZ))<14.5)*((0.9996770000)+(0.0014770600)*((fEvSel_VtxZ)-(0.0003504880))+(-0.0000044498)*TMath::Power((fEvSel_VtxZ)-(0.0003504880),2)+(0.0000058547)*TMath::Power((fEvSel_VtxZ)-(0.0003504880),3)+(-0.0000006434)*TMath::Power((fEvSel_VtxZ)-(0.0003504880),4))+((fEvSel_VtxZ)>14.5)*((0.4462690000)+(0.1055920000)*((fEvSel_VtxZ)-(4.0898700000))+(-0.0049480200)*TMath::Power((fEvSel_VtxZ)-(4.0898700000),2)));
                
                TFile *Calibfile = TFile::Open("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/data/OADB-LHC15o.root");
                
                AliOADBContainer * MultContainer = (AliOADBContainer*)Calibfile->Get("MultSel");
                
                AliOADBMultSelection *OadbMultSelection = (AliOADBMultSelection*) MultContainer->GetObject(aodHeader->GetRunNumber(), "Default");
                
                TH1F *CalibHisto = OadbMultSelection->GetCalibHisto("hCalib_V0M");
                                
                centVal = CalibHisto->GetBinContent(CalibHisto->FindBin(NewCent));
                
            }
            else
            {
             
                AliCentrality *cent = aodHeader->GetCentralityP();
                centVal = cent->GetCentralityPercentile("V0M");
            }
            
            if (centVal < fMinCentrality || centVal >= fMaxCentrality)
            {
                AliDebug(2,Form("Event rejected due to centrality selection. Event centrality: %f, centrality range selection: %f to %f",
                                centVal, fMinCentrality, fMaxCentrality));
                return kFALSE;
            }
        }
    }
    
    // Vertex selection
    if (fAODVertex)
    {
        Double_t vert[3]={0};
        ((AliVVertex*)fAODVertex->At(0))->GetXYZ(vert);
        if (TMath::Abs(vert[2]) > fZVertexCut)
        {
            AliDebug(2,Form("Event rejected due to Z vertex selection. Event Z vertex: %f, Z vertex cut: %f",
                            vert[2], fZVertexCut));
            return kFALSE;
        }
        Double_t dist = TMath::Sqrt((vert[0]-fVertex[0])*(vert[0]-fVertex[0])+(vert[1]-fVertex[1])*(vert[1]-fVertex[1])+(vert[2]-fVertex[2])*(vert[2]-fVertex[2]));
        if (dist > fMaxVertexDist)
        {
            AliDebug(2,Form("Event rejected because the distance between the current and embedded verteces is > %f. "
                            "Current event vertex (%f, %f, %f), embedded event vertex (%f, %f, %f). Distance = %f",
                            fMaxVertexDist, fVertex[0], fVertex[1], fVertex[2], vert[0], vert[1], vert[2], dist));
            return kFALSE;
        }
        
    }
    
    return kTRUE;
}


//_______________________________________________________________________________
Bool_t AliEmbeddingEventForHFTask::Run()
{
  //
  // Analysis execution
  //
    Bool_t HFEvent = kFALSE;
    
    if(fCheckDmeson) HFEvent = CheckForDmesons();
    
    if(fCheckDmeson==kTRUE && HFEvent==kFALSE) return kFALSE;
    
    if (!GetNextEntry()) {
        AliError("Unable to get the AOD event to embed. Nothing will be embedded.");
        return kFALSE;
    }
    
    //cout << "Base Event passed cuts" << endl;
    
    fAODFilePath->SetString(fCurrentAODFile->GetName());
    /*
    Printf("*** Starting Embedding ***");
    Printf(Form("Current File: %s",fCurrentAODFile->GetName()));
    
    Printf(Form("Number of particles from the event: %d",fEmbeddedTracks->GetEntriesFast()));
    */
    TIter next(fEmbeddedTracks);
    
    TRandom2 trackRes(0);
    
    if(fBaseTrackEff)
    {
        TFile *trackEffFile = TFile::Open(fBaseTrackEffPath);
        fBaseTrackEffHist = (TH1F*)trackEffFile->Get(fBaseTrackEffHistName);
    }
    
    AliParticleContainer *otracks = GetParticleContainer(0);
    
    otracks->ResetCurrentID();
    AliAODTrack* otrack = 0;
    Int_t on = 0;
    while ((otrack = static_cast<AliAODTrack*>(otracks->GetNextAcceptParticle())))
    {
        if(otrack)
        {
            if(fBaseTrackEff)
            {
                Double_t TrackProp = trackRes.Uniform(0,1);
                if(TrackProp<=fBaseTrackEffHist->GetBinContent(fBaseTrackEffHist->FindBin(otrack->Pt())))
                {
                    //Printf(Form("Accepter with prob: %lf and Track Eff: %lf",TrackProp,fBaseTrackEffHist->GetBinContent(fBaseTrackEffHist->FindBin(otrack->Pt()))));
                    new ((*fEmbeddedTracks)[on]) AliAODTrack(*otrack);
                    on++;
                    fBaseTrackEffAccepted->Fill(otrack->Pt());
                    fBaseTrackEffAll->Fill(otrack->Pt());
                    
                }
                else
                {
                    //Printf(Form("Rejected with prob: %lf and Track Eff: %lf",TrackProp,fBaseTrackEffHist->GetBinContent(fBaseTrackEffHist->FindBin(otrack->Pt()))));
                    fBaseTrackEffAll->Fill(otrack->Pt());
                }
            }
            else
            {
                new ((*fEmbeddedTracks)[on]) AliAODTrack(*otrack);
                on++;
            }
        }
    }
    
    
    if(fEmbTrackEff)
    {
        TFile *trackEffFile = TFile::Open(fEmbTrackEffPath);
        fEmbTrackEffHist = (TH1F*)trackEffFile->Get(fEmbTrackEffHistName);
    }
    
    
    if (fAODTracks)
    {
        AliDebug(3, Form("%d tracks will be processed for embedding.", fAODTracks->GetEntriesFast()));
        for (Int_t i = on; i < fAODTracks->GetEntriesFast(); i++)
        {
            AliAODTrack *atrack = static_cast<AliAODTrack*>(fAODTracks->At(i));
            if (atrack)
            {
                if(fEmbTrackEff)
                {
                    Double_t TrackProp = trackRes.Uniform(0,1);
                    if(TrackProp<=fEmbTrackEffHist->GetBinContent(fEmbTrackEffHist->FindBin(atrack->Pt())))
                    {
                        new ((*fEmbeddedTracks)[i]) AliAODTrack(*atrack);
                        fEmbTrackEffAccepted->Fill(atrack->Pt());
                        fEmbTrackEffAll->Fill(atrack->Pt());
                    }
                    else
                    {
                        fEmbTrackEffAll->Fill(atrack->Pt());
                    }
                }
                else
                {
                    //cout << "Embedding Track: " << i << endl;
                    new ((*fEmbeddedTracks)[i]) AliAODTrack(*atrack);
                }
            }
            else
            {
                AliError(Form("Could not find track %d in branch %s of tree %s!", i, fAODTrackName.Data(), fAODTreeName.Data()));
            }
        }
    }

    //Printf(Form("Number of particles from the event: %d",fEmbeddedTracks->GetEntriesFast()));
    PostData(1, fOutput);
    
    return kTRUE;
}
void AliEmbeddingEventForHFTask::AddObjectToEvent(TObject *obj, Bool_t attempt)
{
    // Add object to event
    
    if (!(InputEvent()->FindListObject(obj->GetName())))
    {
        InputEvent()->AddObject(obj);
    }
    else
    {
        if (!attempt)
        {
            AliFatal(Form("%s: Container with name %s already present. Aborting", GetName(), obj->GetName()));
        }
    }
}

Bool_t AliEmbeddingEventForHFTask::CheckForDmesons()
{
    
    // fix for temporary bug in ESDfilter
    // the AODs with null vertex pointer didn't pass the PhysSel
    if (!fAodEvent->GetPrimaryVertex() || TMath::Abs(fAodEvent->GetMagneticField()) < 0.001) return kFALSE;
    
    //Event selection
    Bool_t iseventselected = fCuts->IsEventSelected(fAodEvent);
    if (!iseventselected) return kFALSE;
    
    
    const Int_t nD = fArrayDStartoD0pi->GetEntriesFast();
    Bool_t HFEvent = kFALSE;
    
    //Fill the vertex info of the candidates
    AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
    
    for (Int_t icharm = 0; icharm < nD; icharm++) {   //loop over D candidates
        Int_t isSelected = 0;
        
        AliAODRecoDecayHF2Prong* charmCand = static_cast<AliAODRecoDecayHF2Prong*>(fArrayDStartoD0pi->At(icharm)); // D candidates
        if (!charmCand) continue;
        
        Int_t nprongs = charmCand->GetNProngs();
        AliDebug(2, Form("Candidate is %d, and nprongs = %d", fCandidateType, nprongs));
        
        AliAODRecoCascadeHF* dstar = 0;
        
        if (fCandidateType == kDstartoKpipi) {
            dstar = dynamic_cast<AliAODRecoCascadeHF*>(charmCand);
            if (!dstar) {
                Error("AliAnalysisTaskSEDmesonsFilterCJ1::UserExec","Candidate type is D* but object type is wrong (should be AliAODRecoCascadeHF)");
                continue;
            }
        }
        
        Int_t pdgMeson = 413;
        if (fCandidateType == kD0toKpi) pdgMeson = 421;
        
        
        
        if(!(vHF->FillRecoCand(fAodEvent,charmCand))) continue;
        
        
        // region of interest + cuts
        if (!fCuts->IsInFiducialAcceptance(charmCand->Pt(), charmCand->Y(pdgMeson))) continue;
        
        //candidate selected by cuts and PID
        isSelected = fCuts->IsSelected(charmCand, AliRDHFCuts::kAll, fAodEvent); //selected
        
        if (!isSelected) continue;
        
        
        AliAODMCParticle* charmPart = 0;
        Int_t mcLabel = charmCand->MatchToMC(421, fMCarray, fNProngs, fPDGdaughters);
        
        if (mcLabel >= 0) {
            charmPart = static_cast<AliAODMCParticle*>(fMCarray->At(mcLabel));
        }
        
        
        if (charmPart) {
            Int_t origin = CheckOrigin(charmPart, fMCarray);
            if (origin < 0) continue;
            
            if (fRejectQuarkNotFound && origin == kQuarkNotFound) {
                continue;
            }
            if (fRejectDfromB && origin == kFromBottom) {
                continue;
            }
            if (fKeepOnlyDfromB && origin != kFromBottom) {
                continue;
            }
        }
        else {
            if(fDSignalOnly==kTRUE) continue;
        }
        
        HFEvent = kTRUE;
        break;
        
    }
    
    delete vHF;
    
    return HFEvent;
}

//_________________________________________________________________________________________________
Int_t AliEmbeddingEventForHFTask::CheckOrigin(AliAODRecoDecay* cand, TClonesArray* mcArray)
{
    // Checks whether the mother of the D meson candidate comes from a charm or a bottom quark.
    
    if (!mcArray) return -1;
    
    Int_t labDau0 = static_cast<AliVTrack*>(cand->GetDaughter(0))->GetLabel();
    if (labDau0 < 0) return -1;
    
    AliAODMCParticle* part = static_cast<AliAODMCParticle*>(mcArray->At(labDau0));
    return CheckOrigin(part, mcArray);
}

//_________________________________________________________________________________________________
Int_t AliEmbeddingEventForHFTask::CheckOrigin(AliAODRecoDecay* cand, AliStack* stack)
{
    // Checks whether the mother of the D meson candidate comes from a charm or a bottom quark.
    
    if (!stack) return -1;
    
    Int_t labDau0 = static_cast<AliVTrack*>(cand->GetDaughter(0))->GetLabel();
    if (labDau0 < 0) return -1;
    
    return CheckOrigin(labDau0, stack);
}

//_________________________________________________________________________________________________
Int_t AliEmbeddingEventForHFTask::CheckOrigin(AliAODMCParticle* part, TClonesArray* mcArray)
{
    // Checks whether the mother of the particle comes from a charm or a bottom quark.
    
    if (!part) return -1;
    if (!mcArray) return -1;
    
    Int_t pdgGranma = 0;
    Int_t mother = part->GetMother();
    Int_t istep = 0;
    Int_t abspdgGranma = 0;
    Bool_t isFromB = kFALSE;
    Bool_t isQuarkFound = kFALSE;
    
    while (mother >= 0) {
        istep++;
        AliAODMCParticle* mcGranma = static_cast<AliAODMCParticle*>(mcArray->At(mother));
        if (mcGranma) {
            pdgGranma = mcGranma->GetPdgCode();
            abspdgGranma = TMath::Abs(pdgGranma);
            if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
                isFromB = kTRUE;
            }
            
            if (abspdgGranma == 4 || abspdgGranma == 5) isQuarkFound = kTRUE;
            mother = mcGranma->GetMother();
        }
        else {
            ::Error("AliAnalysisTaskSEDmesonsFilterCJ1::CheckOrigin", "Could not retrieve mother particle %d!", mother);
            break;
        }
    }
    
    if (isQuarkFound) {
        if (isFromB) {
            return kFromBottom;
        }
        else {
            return kFromCharm;
        }
    }
    else {
        return kQuarkNotFound;
    }
}

//_________________________________________________________________________________________________
Int_t AliEmbeddingEventForHFTask::CheckOrigin(Int_t ipart, AliStack* stack)
{
    // Checks whether the mother of the particle comes from a charm or a bottom quark.
    
    if (!stack) return -1;
    
    TParticle* part = stack->Particle(ipart);
    if (!part) return -1;
    
    Int_t pdgGranma = 0;
    Int_t mother = part->GetFirstMother();
    Int_t istep = 0;
    Int_t abspdgGranma = 0;
    Bool_t isFromB = kFALSE;
    Bool_t isQuarkFound = kFALSE;
    
    while (mother >= 0) {
        istep++;
        TParticle* mcGranma = stack->Particle(mother);
        if (mcGranma) {
            pdgGranma = mcGranma->GetPdgCode();
            abspdgGranma = TMath::Abs(pdgGranma);
            if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
                isFromB = kTRUE;
            }
            
            if (abspdgGranma == 4 || abspdgGranma == 5) isQuarkFound = kTRUE;
            mother = mcGranma->GetFirstMother();
        }
        else {
            ::Error("AliAnalysisTaskSEDmesonsFilterCJ1::CheckOrigin", "Could not retrieve mother particle %d!", mother);
            break;
        }
    }
    
    if (isQuarkFound) {
        if (isFromB) {
            return kFromBottom;
        }
        else {
            return kFromCharm;
        }
    }
    else {
        return kQuarkNotFound;
    }
}

