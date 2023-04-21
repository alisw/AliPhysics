/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Archita Rani Dash                                              *
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
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "THnSparse.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include <TObjString.h>
#include <math.h>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <TSystem.h>
#include <TDatabasePDG.h>

#include "AliEmcalJetTask.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliAODTrack.h"
//#include <AliPIDResponse.h>
#include "AliMultSelection.h"
#include "AliEmcalJet.h"
#include "AliEmcalList.h"
#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliParticleContainer.h"
#include "AliMCEvent.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliEmcalTrackSelection.h"
#include "AliVTrack.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAnalysisTaskEmcalJetValidation.h"
#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEmcalESDHybridTrackCuts.h"
#include "AliEmcalTrackSelResultHybrid.h"


#include "AliFJWrapper.h"

ClassImp(AliAnalysisTaskEmcalJetValidation) // classimp: necessary for root

//using namespace std;            // std namespace: so you can do things like 'cout'


AliAnalysisTaskEmcalJetValidation::AliAnalysisTaskEmcalJetValidation() :
//   AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetValidation", kTRUE),
   AliAnalysisTaskSE("AliAnalysisTaskEmcalJetValidation"),
   fESD(0),
   fOutputList(0),
   fHistJetPt(0),
   fHistJetPhi(0),
   fHistJetEta(0),
   fHistNEvents(0),
   fHistTrackPt(0),
   fHistTrackPhi(0),
   fHistTrackEta(0),
   fFastJetWrapper(0),
   fTrackCuts(0),
   fInitializedLocal(0),
   fMinPt(0.15),
   fJetEtaRange(0.5),
   fJetR(0.4),
   fJetAlgo(AliJetContainer::antikt_algorithm),
   fGhostArea(0.005),
   fRecoScheme(AliJetContainer::E_scheme)
{

}

//_____________________________________________________________________________
AliAnalysisTaskEmcalJetValidation::AliAnalysisTaskEmcalJetValidation(const char* name) :
   //AliAnalysisTaskEmcalJet(name, kTRUE),
   AliAnalysisTaskSE("AliAnalysisTaskEmcalJetValidation"),
   fESD(0),
   fOutputList(0),
   fHistJetPt(0),
   fHistJetPhi(0),
   fHistJetEta(0),
   fHistNEvents(0),
   fHistTrackPt(0),
   fHistTrackPhi(0),
   fHistTrackEta(0),
   fFastJetWrapper(0),
   fTrackCuts(0),
   fInitializedLocal(0),
   fMinPt(0.15),
   fJetEtaRange(0.5),
   fJetR(0.4),
   fJetAlgo(AliJetContainer::antikt_algorithm),
   fGhostArea(0.005),
   fRecoScheme(AliJetContainer::E_scheme)


{
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskEmcalJetValidation::~AliAnalysisTaskEmcalJetValidation()
{
    // destructor
    if(fOutputList && !fOutputList->IsOwner()) {
      delete fHistJetPt;
      delete fHistJetPhi;
      delete fHistJetEta;
      delete fHistNEvents;
      delete fHistTrackPt;
      delete fHistTrackPhi;
      delete fHistTrackEta;
    }

    delete fOutputList;
    delete fFastJetWrapper;
    delete fTrackCuts;
}

//_________________________________________________
AliAnalysisTaskEmcalJetValidation* AliAnalysisTaskEmcalJetValidation::AddTask(TString suffix, TString jsonconfigfile, Bool_t readMC)
{

   // #### DEFINE MANAGER AND DATA CONTAINER NAMES
   AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
   if(!manager){
      ::Error("AliAnalysisTaskEmcalJetValidation", "No analysis manager to connect to.");
      return NULL;
   }

   TString type = manager->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  if(type.Contains("AOD")){
    ::Error("AliAnalysisTaskEmcalJetValidation", "This task requires to run on ESD");
    return NULL;
  }

  // Add MC handler (for kinematics)
  if(readMC){
    AliMCEventHandler* handler = (AliMCEventHandler*)manager->GetMCtruthEventHandler();
    if (!handler) {
      ::Error("AliAnalysisTaskEmcalJetValidation","Macro called with readMC=true but MC handler not present");
      return 0;
    }
  }


   // #### DEFINE MY ANALYSIS TASK
   TString myContName(Form("AliAnalysisTaskEmcalJetValidation"));
   myContName.Append(suffix);

   AliAnalysisTaskEmcalJetValidation* task = new  AliAnalysisTaskEmcalJetValidation(myContName.Data());
   if(jsonconfigfile.Contains("alien://")){
      Bool_t ok = TFile::Cp(jsonconfigfile.Data(),"local_json.txt");
      if(!ok){
         ::Error("AliAnalysisTaskEmcalJetValidation","Copy of JSON file from alien failed");
         jsonconfigfile="";
      }else{ jsonconfigfile="local_json.txt";}
   }
  if(jsonconfigfile != "") task->InitFromJson(jsonconfigfile);


   task->SetDebugLevel(0);   //No debug messages 0

    // output container
   AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChJetSpectra%s", AliAnalysisManager::GetCommonFileName(), myContName.Data()));


   // #### ADD ANALYSIS TASK
   manager->AddTask(task);
   manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
   manager->ConnectOutput(task, 1, contHistos);

   return task;

}

//_________________________________________________
void AliAnalysisTaskEmcalJetValidation::InitFromJson(TString filename){
   /// read configuration from json file
   if (filename != "" && gSystem->Exec(Form("ls %s > /dev/null", filename.Data())) == 0) {
      printf("------Read configuration from JSON file------\n");

      fJetAlgo    = AliJetContainer::antikt_algorithm; //not yet in JSON
      fRecoScheme = AliJetContainer::E_scheme;         //not yet in JSON
      fGhostArea = 0.005;                              //not yet in JSON

      int   njetRadii = 0;
      float* jetRadii = GetJsonArray(filename.Data(), "jet-finder-data", "jetR", njetRadii);
      if(njetRadii > 0){
         fJetR = jetRadii[0];              //set jet radius from JSON
      }else{
         AliFatal("Missing Jet R in JSON, please check it");
      }
   }
   return;
}
//_________________________________________________
void AliAnalysisTaskEmcalJetValidation::ExecOnceLocal(){
    //this function will be call just once from the UserExec
    fInitializedLocal = kTRUE;



    fFastJetWrapper = new AliFJWrapper("FJWrapper", "FJWrapper");     // Initialization of my jet finder
    fFastJetWrapper->SetAreaType(fastjet::active_area);
    fFastJetWrapper->SetGhostArea(fGhostArea); //set ghost area
    fFastJetWrapper->SetR(fJetR);              //set jet radius
    fFastJetWrapper->SetAlgorithm(AliEmcalJetTask::ConvertToFJAlgo(fJetAlgo));
    fFastJetWrapper->SetRecombScheme(AliEmcalJetTask::ConvertToFJRecoScheme(fRecoScheme));


    fJetEtaRange = 0.9 - fJetR; //fiducial cut on jets

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "default");
    fTrackCuts->SetName("Global Hybrid tracks, loose DCA");
    fTrackCuts->SetPtRange(0.15, 1.e15);
    fTrackCuts->SetEtaRange(-0.9, +0.9);
    fTrackCuts->SetRequireITSRefit(kTRUE);
    fTrackCuts->SetRequireTPCRefit(kTRUE);
    //fTrackCuts->SetMinNClustersTPC(70);
    fTrackCuts->SetMinNCrossedRowsTPC(70);    // Marta suggested to use this instead of the cut on l.253
                                              //because basically we do track selections on the no. of crossed rows

    fTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);   // similar to SetMinNCrossedRowsOverFindableClustersTPC(0.8) used in O2
    fTrackCuts->SetMaxChi2PerClusterTPC(4.0);
    fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                          AliESDtrackCuts::kAny);    // similar to SetRequireHitsInITSLayers(1, {0, 1}) in O2
    fTrackCuts->SetMaxChi2PerClusterITS(36.0);
  //  fTrackCuts->SetMaxDCAToVertexXYPtDep("(0.0105 + 0.0350 / TMath::Power(pt, 1.1))"); //similar to "SetMaxDcaXYPtDep" implemented in O2 task (?)

    fTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);              // same as SetRequireGoldenChi2(true) in O2
    fTrackCuts->SetMaxFractionSharedTPCClusters(0.4);           //Not yet in O2

    fTrackCuts->SetAcceptKinkDaughters(kFALSE);                 //Not yet in O2 task
    fTrackCuts->SetMaxDCAToVertexXY(2.4);
    fTrackCuts->SetMaxDCAToVertexZ(3.2);                        // Loose DCA cuts: similar to SetMaxDcaZ in O2 task
    //fTrackCuts->SetDCAToVertex2D(kTRUE);
    fTrackCuts->SetDCAToVertex2D(kFALSE);                          //Marta suggested to set the flag as kFalse
                                                                  //because there's no option to use the 2D cut for the DCA in O2

}
//_________________________________________________
void AliAnalysisTaskEmcalJetValidation::UserCreateOutputObjects()
{
	Printf("Check done %i",__LINE__);

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);


   AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
   if (man){
      AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
      //if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
   }

   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);
   fOutputList->SetName("OutputHistos");

   fHistNEvents = new TH1F("hNEvents", "Number of processed events", 1, 0, 1);

   //JET QA
   fHistJetPt = new TH1F("jetPt", "inclusive jetPt ; p_{T} (GeV/#it{c})", 200, 0, 100);
   fHistJetPhi = new TH1F("jetPhi", "inclusive jetPhi; #phi ", 200, 0, 6.4);
   fHistJetEta = new TH1F("jetEta", "inclusive jetEta; #eta ", 200, -0.9, 0.9);

   //TRACK QA
   fHistTrackPt = new TH1F("jetTrackPt","track Pt;p_{T} (GeV/#it{c})", 200, 0, 100);
   fHistTrackPhi = new TH1F("jetTrackPhi", "track #phi; #phi",200, 0, 6.4);
   fHistTrackEta = new TH1F("jetTrackEta", "track #eta; #eta",200, -0.9, 0.9);

   fOutputList->Add(fHistNEvents);
   fOutputList->Add(fHistJetPt);
   fOutputList->Add(fHistJetPhi);
   fOutputList->Add(fHistJetEta);

   fOutputList->Add(fHistTrackPt);
   fOutputList->Add(fHistTrackPhi);
   fOutputList->Add(fHistTrackEta);

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++
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
//_____________________________________________________________________________

void AliAnalysisTaskEmcalJetValidation::UserExec(Option_t *)
{

   if(!fInitializedLocal) ExecOnceLocal();


   AliVEvent* vevt = InputEvent();
   AliESDEvent* fESD = (AliESDEvent*)(vevt);
   if (!fESD) {
     printf("AliAnalysisTaskEmcalJetValidation::UserExec(): bad ESD\n");
     return;
   }

  //DO SOME EVENT SELECTION HERE
  const AliESDVertex* vertex = (AliESDVertex*)fESD->GetPrimaryVertex();
  if(TMath::Abs(vertex->GetZ()) > 10.) return;

  //EVENTS WHICH PASSED
  fHistNEvents->Fill(0.5);

  fFastJetWrapper->Clear();

  Int_t totTracks = fESD->GetNumberOfTracks();
  TLorentzVector lVec;

  for(Int_t itr = 0; itr < totTracks; itr++) {
     AliESDtrack* track = static_cast< AliESDtrack*>(fESD->GetTrack(itr));      //Feeding my jet finder with tracks to produce jets
     if(!track) continue;
     if(!fTrackCuts->AcceptTrack(track)) continue;

    // lVec.SetPtEtaPhiM(track->Pt(), track->Eta(), track->Phi(), 0.13957);   //assume that track is pion
    // fFastJetWrapper->AddInputVector(lVec.Px(), lVec.Py(), lVec.Pz(), lVec.E()); //fill jet constituents

    // Test for 2 hardcoded tracks : E_scheme
    /* lVec.SetPtEtaPhiM(10, 0, 0, 0.13957);
     fFastJetWrapper->AddInputVector(lVec.Px(), lVec.Py(), lVec.Pz(), lVec.E());

     lVec.SetPtEtaPhiM(10, 0.01, 0.01, 0.13957);
     fFastJetWrapper->AddInputVector(lVec.Px(), lVec.Py(), lVec.Pz(), lVec.E());*/
    // Test for 2 harcoded tracks : pt_scheme

     //Filling Track Histograms
     fHistTrackPt->Fill(track->Pt());
     fHistTrackPhi->Fill(track->Phi());
     fHistTrackEta->Fill(track->Eta());
  }

  fFastJetWrapper->Run();

  std::vector<fastjet::PseudoJet> myJets = fFastJetWrapper->GetInclusiveJets();

  for(UInt_t ijet = 0; ijet < myJets.size(); ++ijet) {
     if(myJets.at(ijet).pt() < fMinPt) continue;  //skip ghost jets
     if(TMath::Abs(myJets.at(ijet).eta()) > fJetEtaRange) continue; //skip jets out of the fiducial volume
     fHistJetPt->Fill(myJets.at(ijet).pt());      //fill jet pT to histogram
     fHistJetPhi->Fill(myJets.at(ijet).phi());   //fill jet phi histogram
     fHistJetEta->Fill(myJets.at(ijet).eta());  //fill jet eta histogram
     //fHistJetPt->Fill(myJets.at(ijet).area());  //fill jet area to histogram
  }



  //_________________________________________________________

  PostData(1, fOutputList);

}
//_____________________________________________________________________________
void AliAnalysisTaskEmcalJetValidation::Terminate(Option_t *)
{
  PostData(1, fOutputList);

  fOutputList = dynamic_cast<AliEmcalList*> (GetOutputData(1)); // '1' refers to the output slot
  if(!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
   }
}

//______________________________________________________________________________
std::string AliAnalysisTaskEmcalJetValidation::GetJsonString(const char* jsonFileName, const char* section, const char* key)
{
  FILE* fj = fopen(jsonFileName, "r");
  char line[500];
  char* value = 0x0;
  bool corrSection = false;
  while (!feof(fj)) {
    char* rc = fgets(line, 500, fj);
    if (rc && strstr(line, section) && !(section && !section[0]))
      corrSection = true;
    if (rc && strstr(line, key)) {
      value = strtok(line, ":");
      value = strtok(NULL, ":");
      if(corrSection || (section && !section[0]))
        break;
    }
  }
  fclose(fj);
  TString sValue = value;
  sValue.ReplaceAll("\",", "");
  sValue.ReplaceAll("\"", "");
  sValue.ReplaceAll("\n", "");
  sValue.ReplaceAll("\t", "");
  sValue.ReplaceAll(" ", "");
  return std::string(sValue.Data());
}
//______________________________________________________________________________
int AliAnalysisTaskEmcalJetValidation::GetJsonBool(const char* jsonFileName, const char* section, const char* key)
{
  FILE* fj = fopen(jsonFileName, "r");
  char line[500];
  int value = -1;
  bool corrSection = false;
  while (!feof(fj)) {
    char* rc = fgets(line, 500, fj);
    if (rc && strstr(line, section) && !(section && !section[0]))
      corrSection = true;
    if (rc && strstr(line, key)) {
      char* token = strtok(line, ":");
      token = strtok(NULL, ":");
      TString temp = token;
      temp.ReplaceAll("\"", "");
      temp.ReplaceAll(",", "");
      if (temp.Contains("true"))
        value = 1;
      if (temp.Contains("false"))
        value = 0;
      if(corrSection || (section && !section[0]))
        break;
    }
  }
  fclose(fj);
  return value;
}
//______________________________________________________________________________

int AliAnalysisTaskEmcalJetValidation::GetJsonInteger(const char* jsonFileName, const char* section, const char* key)
{
  FILE* fj = fopen(jsonFileName, "r");
  char line[500];
  int value = -999;
  bool corrSection = false;
  while (!feof(fj)) {
    char* rc = fgets(line, 500, fj);
    if (rc && strstr(line, section) && !(section && !section[0]))
      corrSection = true;
    if (rc && strstr(line, key)) {
      char* token = strtok(line, ":");
      token = strtok(NULL, ":");
      TString temp = token;
      temp.ReplaceAll("\"", "");
      temp.ReplaceAll(",", "");
      value = temp.Atoi();
      if(corrSection || (section && !section[0]))
        break;
    }
  }
  fclose(fj);
  return value;
}

//______________________________________________________________________________
float AliAnalysisTaskEmcalJetValidation::GetJsonFloat(const char* jsonFileName, const char* section, const char* key)
{
  FILE* fj = fopen(jsonFileName, "r");
  char line[500];
  float value = -999.;
  bool corrSection = false;
  while (!feof(fj)) {
    char* rc = fgets(line, 500, fj);
    if (rc && strstr(line, section) && !(section && !section[0]))
      corrSection = true;
    if (rc && strstr(line, key)) {
      char* token = strtok(line, ":");
      token = strtok(NULL, ":");
      TString temp = token;
      temp.ReplaceAll("\"", "");
      temp.ReplaceAll(",", "");
      value = temp.Atof();
      if(corrSection || (section && !section[0]))
        break;
    }
  }
  fclose(fj);
  return value;
}
//______________________________________________________________________________
float* AliAnalysisTaskEmcalJetValidation::GetJsonArray(const char* jsonFileName, const char* section, const char* key, int& size)
{
  FILE* fj = fopen(jsonFileName, "r");
  char line[500];
  float* arrVals = 0x0;
  bool corrSection = false;
  while (!feof(fj)) {
    fgets(line, 500, fj);
    if (strstr(line, section) && !(section && !section[0]))
      corrSection = true;
    if (strstr(line, key)) {
      TString full = "";
      while (!feof(fj)) {
        fgets(line, 500, fj);
        int len = strlen(line);
        if (line[len - 1] == '\n')
          line[len - 1] = 0;
        full.Append(line);
        if (strstr(line, "}"))
          break;
      }
      full.ReplaceAll("\"values\":", "");
      full.ReplaceAll(" ", "");
      full.ReplaceAll("},", "");
      full.ReplaceAll("}", "");
      TObjArray* arrStr = full.Tokenize(",");
      size = arrStr->GetEntriesFast();
      arrVals = new float[size];
      for (int j = 0; j < size; j++) {
        TObjString* sss = (TObjString*)arrStr->At(j);
        TString strval = sss->GetString();
        strval.ReplaceAll("[", "");
        strval.ReplaceAll("]", "");
        strval.ReplaceAll("\"", "");
        arrVals[j] = strval.Atof();
      }
      arrStr->Delete();
      delete arrStr;
      if(corrSection)
        break;
    }
  }
  return arrVals;
}

//______________________________________________________________________________
float** AliAnalysisTaskEmcalJetValidation::GetJsonMatrix(const char* jsonFileName, const char* section, const char* key, int& size1, int& size2)
{
  FILE* fj = fopen(jsonFileName, "r");
  char line[500];
  float** arrVals = 0x0;
  bool corrSection = false;
  bool moveToValues = false;
  while (!feof(fj)) {
    fgets(line, 500, fj);
    if (strstr(line, section) && !(section && !section[0])){
      corrSection = true;
      if (strstr(section, "selector"))
        moveToValues = true;
    }
    if (strstr(line, key)) {
      TString full = "";
      while (!feof(fj)) {
        fgets(line, 500, fj);
        if (moveToValues && !strstr(line, "values"))
          continue;
        moveToValues = false;
        int len = strlen(line);
        if (line[len - 1] == '\n')
          line[len - 1] = 0;
        full.Append(line);
        if (strstr(line, "}"))
          break;
      }
      full.ReplaceAll("\"values\":", "");
      full.ReplaceAll(" ", "");
      full.ReplaceAll("},", "");
      full.ReplaceAll("}", "");
      TObjArray* rowArrStr = full.Tokenize("]");
      size1 = rowArrStr->GetEntriesFast();
      arrVals = new float*[size1];
      for (int j = 0; j < size1; j++) {
        TObjString* rowStr = (TObjString*)rowArrStr->At(j);
        TString rowStrVal = rowStr->GetString();
        TObjArray* arrStr = rowStrVal.Tokenize(",");
        size2 = arrStr->GetEntriesFast();
        arrVals[j] = new float[size2];
        for (int k = 0; k < size2; k++) {
          TObjString* sss = (TObjString*)arrStr->At(k);
          TString strval = sss->GetString();
          strval.ReplaceAll(",[", "");
          strval.ReplaceAll("[", "");
          strval.ReplaceAll("]", "");
          strval.ReplaceAll("\"", "");
          arrVals[j][k] = strval.Atof();
        }
      }
      if(corrSection)
        break;
    }
  }
  return arrVals;
}
//______________________________________________________________________________
