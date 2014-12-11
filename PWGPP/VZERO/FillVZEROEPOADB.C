#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TString.h"
#include "TProfile.h"
#include "TList.h"
#include "TSystem.h"

#include "AliOADBContainer.h"
#include "AliAnalysisManager.h"
#endif

void FillVZEROEPOADB()
{
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");   
  gSystem->Load("libOADB");

  AliOADBContainer * oadbCont = new AliOADBContainer("vzeroEP");

  TList *defaultList = new TList;
  defaultList->SetName("Default");
  TList *inputList = NULL;
  TProfile *profHisto = NULL;
  TFile fInputDefault("VZERO.EPFlatenning.PS.LHC11h_000170162_p1_muon_.root");
  inputList = (TList*)fInputDefault.Get("coutput");
  for(Int_t i = 0; i < 8; ++i) {
    profHisto = (TProfile*)inputList->FindObject(Form("fX2_%d",i))->Clone(Form("fX2_%d",i));
    profHisto->SetDirectory(0);
    defaultList->Add(profHisto);
    profHisto = (TProfile*)inputList->FindObject(Form("fY2_%d",i))->Clone(Form("fY2_%d",i));
    profHisto->SetDirectory(0);
    defaultList->Add(profHisto);
    profHisto = (TProfile*)inputList->FindObject(Form("fX2Y2_%d",i))->Clone(Form("fX2Y2_%d",i));
    profHisto->SetDirectory(0);
    defaultList->Add(profHisto);
    profHisto = (TProfile*)inputList->FindObject(Form("fCos8Psi_%d",i))->Clone(Form("fCos8Psi_%d",i));
    profHisto->SetDirectory(0);
    defaultList->Add(profHisto);
  }
  fInputDefault.Close();
  oadbCont->AddDefaultObject(defaultList);


  TList *list1 = new TList;
  TFile fInput1("VZERO.EPFlatenning.PS.LHC11h_000169683_p1_muon_ESDs.root");
  inputList = (TList*)fInput1.Get("coutput");
  for(Int_t i = 0; i < 8; ++i) {
    profHisto = (TProfile*)inputList->FindObject(Form("fX2_%d",i))->Clone(Form("fX2_%d",i));
    profHisto->SetDirectory(0);
    list1->Add(profHisto);
    profHisto = (TProfile*)inputList->FindObject(Form("fY2_%d",i))->Clone(Form("fY2_%d",i));
    profHisto->SetDirectory(0);
    list1->Add(profHisto);
    profHisto = (TProfile*)inputList->FindObject(Form("fX2Y2_%d",i))->Clone(Form("fX2Y2_%d",i));
    profHisto->SetDirectory(0);
    list1->Add(profHisto);
    profHisto = (TProfile*)inputList->FindObject(Form("fCos8Psi_%d",i))->Clone(Form("fCos8Psi_%d",i));
    profHisto->SetDirectory(0);
    list1->Add(profHisto);
  }
  oadbCont->AppendObject(list1, 169683, 169683);


  TString oadbFileName = Form("%s/COMMON/EVENTPLANE/data/vzero.root", AliAnalysisManager::GetOADBPath());
  oadbCont->WriteToFile(oadbFileName.Data());
}
