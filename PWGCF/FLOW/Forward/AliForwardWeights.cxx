#include "AliAnalysisManager.h"
#include "TH3F.h"
#include "AliAnalysisDataContainer.h"
#include "AliForwardWeights.h"
#include "TGrid.h"
#include <iostream>

AliForwardWeights::AliForwardWeights():
  fSettings()
  {
  }

void AliForwardWeights::connectNUA(){
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

    std::cout << "AliForwardWeights connecting NUA weights" << std::endl;
    TObjArray* taskContainers = mgr->GetContainers();
    AliAnalysisDataContainer* weights;
    
    // Looking for name of the .root file without path
    TObjArray *tx = fSettings.nua_file.Tokenize("/");
    TObjArray *ty = ((TObjString *)(tx->At(tx->GetEntries()-1)))->String().Tokenize(".");
    TString nuaobject =  ((TObjString *)(ty->At(0)))->String();
    
    std::cout << "AliForwardWeights: NUA weights = " << nuaobject << std::endl;

    weights = (AliAnalysisDataContainer*) taskContainers->FindObject(nuaobject);
    
    if (!weights) {
      std::cout << "I-AliForwardWeights: " << nuaobject << " NUA weights not defined - reading now. " << std::endl;
      weights = makeWeightContainer(fSettings.nua_file,nuaobject);
    }
    connectContainer( weights);
  }

void AliForwardWeights::connectNUE(){
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

    std::cout << "AliForwardWeights connecting NUE weights" << std::endl;

    if (fSettings.nue_file != ""){
          TObjArray* taskContainers = mgr->GetContainers();
      AliAnalysisDataContainer* nue_weights;
      
    TObjArray *tx = fSettings.nue_file.Tokenize("/");
    TObjArray *ty = ((TObjString *)(tx->At(tx->GetEntries()-1)))->String().Tokenize(".");
    TString nueobject =  ((TObjString *)(ty->At(0)))->String();

    std::cout << "AliForwardWeights: NUE weights = " << nueobject << std::endl;

      nue_weights = (AliAnalysisDataContainer*) taskContainers->FindObject(nueobject);

      if (!nue_weights) {
      std::cout << "I-AliForwardWeights: " << nueobject << " NUE weights not defined - reading now. " << std::endl;
        nue_weights = makeWeightContainerNUE(fSettings.nue_file,nueobject);   
      }
      connectNUEContainer(nue_weights);
    }
  }


void AliForwardWeights::connectSec(){
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

  if (fSettings.sec_file != ""){
    TObjArray* taskContainers = mgr->GetContainers();
    AliAnalysisDataContainer* sec_weights;
    
    sec_weights = (AliAnalysisDataContainer*) taskContainers->FindObject("sec");
    std::cout << "AliForwardWeights: Secondary weights = " << sec_weights << std::endl;

    if (!sec_weights) {
      std::cout << "I-AliForwardWeights: " << "secondary weights not defined - reading now. " << std::endl;

      sec_weights = makeWeightContainerSec(fSettings.sec_file,"sec");   
    }
    connectSecContainer(sec_weights);
  }
}

void AliForwardWeights::connectSecCent(){
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

  if (fSettings.sec_cent_file != ""){
    TObjArray* taskContainers = mgr->GetContainers();

    AliAnalysisDataContainer* sec_weights_cent;
    
    sec_weights_cent = (AliAnalysisDataContainer*) taskContainers->FindObject("sec_cent");
    std::cout << "AliForwardWeights: Secondary weights = " << sec_weights_cent << std::endl;

    if (!sec_weights_cent) {
      sec_weights_cent = makeWeightContainerSecCent(fSettings.sec_cent_file,"sec_cent");   
      std::cout << "I-AliForwardWeights: " << "secondary centrality weights not defined - reading now. " << std::endl;
    }      
    connectSecCentContainer(sec_weights_cent);
  }
}


AliAnalysisDataContainer* AliForwardWeights::makeWeightContainer(TString nua_file, TString containerName){
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  
  AliAnalysisDataContainer* weights;
  if (nua_file.Contains("alien:")) TGrid::Connect("alien:");
  TFile* file;
  file = TFile::Open(nua_file.Data(), "READ");

  if(!file) { printf("E-AliForwardWeights: Input file with differential weights not found!\n"); return NULL; }

  TList* weights_list = new TList();
  weights_list->SetName("nuaWeights");
  
  TH3F* nuacentral = new TH3F();
  TH3F* nuaforward = new TH3F();

  file->GetObject("nuacentral", nuacentral);
  nuacentral->SetDirectory(0);
  nuacentral->SetNameTitle("nuacentral","nuacentral");

  file->GetObject("nuaforward", nuaforward);
  nuaforward->SetDirectory(0);
  nuaforward->SetNameTitle("nuaforward","nuaforward");
  file->Close();

  weights_list->Add(nuacentral);
  weights_list->Add(nuaforward);

  weights = mgr->CreateContainer(containerName,TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));
  weights->SetData(weights_list);
  return weights;
}


AliAnalysisDataContainer* AliForwardWeights::makeWeightContainerNUE(TString nue_file, TString containerName){
      AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

  AliAnalysisDataContainer* weights;
  if (nue_file.Contains("alien:") ) TGrid::Connect("alien:");
  TFile* file;
  //if (sec_file.Contains(containerName))
  file = TFile::Open(nue_file.Data(), "READ");

  if(!file) { printf("E-AliForwardWeights: Input file with secondary weights not found!\n"); return NULL; }

  TList* weights_list = new TList();
  weights_list->SetName("nue_corr");
  
  TH3F* nuehist = new TH3F();

  file->GetObject("correction", nuehist);
  nuehist->SetDirectory(0);
  nuehist->SetNameTitle("nue_corr","nue_corr");

  file->Close();

  weights_list->Add(nuehist);

  weights = mgr->CreateContainer(containerName,TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));
  weights->SetData(weights_list);
  return weights;
}

AliAnalysisDataContainer* AliForwardWeights::makeWeightContainerSec(TString sec_file, TString containerName){
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  
  AliAnalysisDataContainer* weights;
  if (sec_file.Contains("alien:") ) TGrid::Connect("alien:");
  TFile* file;
  file = TFile::Open(sec_file.Data(), "READ");

  if(!file) { printf("E-AliForwardWeights: Input file with secondary weights not found!\n"); return NULL; }

  TList* weights_list = new TList();
  weights_list->SetName("sec_corr");
  
  TH3F* sechist = new TH3F();

  file->GetObject("correction", sechist);
  sechist->SetDirectory(0);
  sechist->SetNameTitle("sec_corr","sec_corr");

  file->Close();

  weights_list->Add(sechist);

  weights = mgr->CreateContainer(containerName,TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));
  weights->SetData(weights_list);
  return weights;
}



AliAnalysisDataContainer* AliForwardWeights::makeWeightContainerSecCent(TString sec_file, TString containerName){
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  
  AliAnalysisDataContainer* weights;
  if (sec_file.Contains("alien:") ) TGrid::Connect("alien:");
  TFile* file;
  file = TFile::Open(sec_file.Data(), "READ");

  if(!file) { printf("E-AliForwardWeights: Input file with secondary weights not found!\n"); return NULL; }

  TList* weights_list = new TList();
  weights_list->SetName("sec_corr_cent");
  
  TH3F* sechist = new TH3F();

  file->GetObject("correction", sechist);
  sechist->SetDirectory(0);
  sechist->SetNameTitle("sec_corr_cent","sec_corr_cent");

  file->Close();

  weights_list->Add(sechist);

  weights = mgr->CreateContainer(containerName,TList::Class(), AliAnalysisManager::kInputContainer,Form("%s", mgr->GetCommonFileName()));
  weights->SetData(weights_list);
  return weights;
}


void AliForwardWeights::connectContainer(AliAnalysisDataContainer* container){

  fSettings.nuacentral = static_cast<TH3F*>( static_cast<TList*>(container->GetData())->FindObject("nuacentral") );
  fSettings.nuaforward = static_cast<TH3F*>( static_cast<TList*>(container->GetData())->FindObject("nuaforward") );
  fSettings.nuacentral->SetDirectory(0);
  fSettings.nuaforward->SetDirectory(0);
}

void AliForwardWeights::connectSecContainer(AliAnalysisDataContainer* container){

  fSettings.seccorr_fwd = static_cast<TH3F*>( static_cast<TList*>(container->GetData())->FindObject("sec_corr") );
  fSettings.seccorr_fwd->SetDirectory(0);
}


void AliForwardWeights::connectNUEContainer(AliAnalysisDataContainer* container){

  fSettings.nuehist = static_cast<TH3F*>( static_cast<TList*>(container->GetData())->FindObject("nue_corr") );
  fSettings.nuehist->SetDirectory(0);
}


void AliForwardWeights::connectSecCentContainer(AliAnalysisDataContainer* container){

  fSettings.seccorr_cent = static_cast<TH3F*>( static_cast<TList*>(container->GetData())->FindObject("sec_corr_cent") );
  fSettings.seccorr_cent->SetDirectory(0);
}