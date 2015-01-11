#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TString.h"
#include "TProfile.h"
#include "TList.h"
#include "TSystem.h"

#include "AliOADBContainer.h"
#include "AliAnalysisManager.h"
#endif

void FillVZEROEPOADBFull(const char* filename = "AOD083.txt", Bool_t mbOnly = kFALSE)
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

  {
    TList *defaultList = new TList;
    defaultList->SetName("Default");
    TProfile *profHisto = NULL;
    TFile fInputDefault("minbias/VZERO.EPFlatenning.PS.LHC11h_AOD083_000170162.root");
    TList *inputList = (TList*)fInputDefault.Get("coutput");
    for(Int_t i = 0; i < 11; ++i) {
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
    printf("Run 170162 filled\n");
  }

  {
    TList *list1 = new TList;
    TProfile *profHisto = NULL;
    TFile fInput1("minbias/VZERO.EPFlatenning.PS.LHC11h_AOD083_000169683.root");
    TList *inputList = (TList*)fInput1.Get("coutput");
    for(Int_t i = 0; i < 11; ++i) {
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
    printf("Run 169683 filled\n");
  }

  // loop of over all other runs
  Int_t runList[500];
  ifstream *fruns = new ifstream (filename);
  if (!*fruns) return;
  TString strLine;
  Int_t count = 0;
  while (strLine.ReadLine(*fruns)) {
    runList[count++] = strLine.Atoi();
  }
  delete fruns;

  for(Int_t irun = 0; irun < count; ++irun) {
    TList *list2 = new TList;
    TProfile *profHisto = NULL;
    TFile fInput2(Form("csemi/VZERO.EPFlatenning.PS.LHC11h_AOD083_000%d.root",runList[irun]));
    TList *inputList = (TList*)fInput2.Get("coutput");
    TFile fInput3(Form("cpbi2/VZERO.EPFlatenning.PS.LHC11h_AOD083_000%d.root",runList[irun]));
    TList *inputListBis = (TList*)fInput3.Get("coutput");
    for(Int_t i = 0; i < 11; ++i) {
      profHisto = (TProfile*)inputList->FindObject(Form("fX2_%d",i))->Clone(Form("fX2_%d",i));
      profHisto->SetDirectory(0);
      Int_t ibin = profHisto->FindBin(62.5);
      profHisto->SetBinContent(ibin,0);
      profHisto->SetBinError(ibin,0);
      profHisto->SetBinEntries(ibin,0);
      if (mbOnly) {
	profHisto = (TProfile*)inputListBis->FindObject(Form("fX2_%d",i))->Clone(Form("fX2_%d",i));
	profHisto->SetDirectory(0);
      }
      else
	profHisto->Add((TProfile*)inputListBis->FindObject(Form("fX2_%d",i)));
      list2->Add(profHisto);

      profHisto = (TProfile*)inputList->FindObject(Form("fY2_%d",i))->Clone(Form("fY2_%d",i));
      profHisto->SetDirectory(0);
      profHisto->SetBinContent(ibin,0);
      profHisto->SetBinError(ibin,0);
      profHisto->SetBinEntries(ibin,0);
      if (mbOnly) {
	profHisto = (TProfile*)inputListBis->FindObject(Form("fY2_%d",i))->Clone(Form("fY2_%d",i));
	profHisto->SetDirectory(0);
      }
      else
	profHisto->Add((TProfile*)inputListBis->FindObject(Form("fY2_%d",i)));
      list2->Add(profHisto);

      profHisto = (TProfile*)inputList->FindObject(Form("fX2Y2_%d",i))->Clone(Form("fX2Y2_%d",i));
      profHisto->SetDirectory(0);
      profHisto->SetBinContent(ibin,0);
      profHisto->SetBinError(ibin,0);
      profHisto->SetBinEntries(ibin,0);
      if (mbOnly) {
	profHisto = (TProfile*)inputListBis->FindObject(Form("fX2Y2_%d",i))->Clone(Form("fX2Y2_%d",i));
	profHisto->SetDirectory(0);
      }
      else
	profHisto->Add((TProfile*)inputListBis->FindObject(Form("fX2Y2_%d",i)));
      list2->Add(profHisto);

      profHisto = (TProfile*)inputList->FindObject(Form("fCos8Psi_%d",i))->Clone(Form("fCos8Psi_%d",i));
      profHisto->SetDirectory(0);
      profHisto->SetBinContent(ibin,0);
      profHisto->SetBinError(ibin,0);
      profHisto->SetBinEntries(ibin,0);
      if (mbOnly) {
	profHisto = (TProfile*)inputListBis->FindObject(Form("fCos8Psi_%d",i))->Clone(Form("fCos8Psi_%d",i));
	profHisto->SetDirectory(0);
      }
      else
	profHisto->Add((TProfile*)inputListBis->FindObject(Form("fCos8Psi_%d",i)));
      list2->Add(profHisto);
    }
    oadbCont->AppendObject(list2, runList[irun], runList[irun]);
    printf("Run %d filled\n",runList[irun]);
  }

  TString oadbFileName = Form("%s/COMMON/EVENTPLANE/data/vzero.root", AliAnalysisManager::GetOADBPath());
  oadbCont->WriteToFile(oadbFileName.Data());
}
