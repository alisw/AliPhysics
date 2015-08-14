#include "TCanvas.h"
#include "TH2D.h"
#include "TFile.h"
#include "TPRegexp.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "AliOADBContainer.h"

#include <iostream>

//#define TEST

Int_t addMapToFile(TString filePathNameMapToAdd, TString normalisation, TString period, Int_t pass, Bool_t isMC, Bool_t isSigma, 
                   TString pathMapPackage = "etaMaps", TString fileNameMapPackage = "TPCetaMaps.root")
{
  //normalisation = "NoNormalisation", "SmallEtaNormalisation", "LargeEtaNormalisation"
  
  TString mapNameInFile = "";
  TString par0StringName = "";
  TString mapName = "";
  TString mapTitle = "";
  
  Int_t* runLow = 0x0;
  Int_t* runUp = 0x0;
  Int_t nRunRanges = 1;
  
  TString dataType = "DATA";
  if (isMC)
    dataType = "MC";
  
  if (isSigma) {
    mapName = "sigmaPar1Map";
    //mapNameInFile = "hThetaMapSigmaPar1";
    mapNameInFile = Form("hSigmaPar1_%s_extrapolated", normalisation.Data());
    par0StringName = "c0";
  }
  else {
    mapName = Form("TPCetaMaps_%s_pass%d", dataType.Data(), pass);
    //mapNameInFile = Form("hRefined%s", normalisation.Data());
    mapNameInFile = Form("hRes3DprimeFit_%s_extrapolated", normalisation.Data());
  }
  
  
  
  
  //TPRegexp reg(".*(LHC1[1-2][a-z]+[0-9]+[a-z_]*)/.*");
  
  
  // Find the run range from the period
  period.ToUpper();
  
  Bool_t addDefault = kFALSE;
  TString defaultObj = Form("Default_%s_pass%d", dataType.Data(), pass);
  
  if (period.Contains("LHC10B") == 1 && !isMC) {
    mapTitle = Form("LHC10b.pass%d", pass);
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 114650;
    runUp[0] = 117630;
  }
  else if (period.Contains("LHC10C") == 1 && !isMC) {
    mapTitle = Form("LHC10c.pass%d", pass);
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 117631;
    runUp[0] = 121526;
  }
  else if (period.Contains("LHC10D") == 1 && !isMC) {
    mapTitle = Form("LHC10d.pass%d", pass);
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 121527;
    runUp[0] = 126460;
  }
  else if (period.Contains("LHC10E") == 1 && !isMC) {
    mapTitle = Form("LHC10e.pass%d", pass);
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 126461;
    runUp[0] = 130930;
  }
  else if (period.Contains("LHC10F") == 1 && !isMC) {
    mapTitle = Form("LHC10f.pass%d", pass);
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 130931;
    runUp[0] = 135393;
  }
  else if (period.Contains("LHC10G") == 1 && !isMC) {
    mapTitle = Form("LHC10g.pass%d", pass);
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 135394;
    runUp[0] = 136781;
  }
  else if (period.Contains("LHC11A_7TEV") == 1 && !isMC) {
    //LHC11a -> 7TeV
    mapTitle = Form("LHC11a_7TeV.pass%d", pass);
  
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 139847;
    runUp[0] = 146631;
  }
  else if (period.Contains("LHC11A_2.76TEV") == 1 && !isMC) {
    //LHC11a -> 2.76TeV
    mapTitle = Form("LHC11a_2.76TeV.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 146632;
    runUp[0] = 146974;
  }
  else if (period.Contains("LHC10H") == 1 && !isMC) {
    //LHC10h -> 2.76 ATeV (PbPb)
    mapTitle = Form("LHC10h.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 136782;
    runUp[0] = 139846;
  }
  else if (period.Contains("LHC11H") == 1 && !isMC) {
    //LHC11h -> 2.76 ATeV (PbPb)
    mapTitle = Form("LHC11h.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 165772;
    runUp[0] = 170718;
  }
  else if (period.Contains("LHC10D1") == 1 && isMC) {
    period = "LHC10D1";
    mapTitle = Form("LHC10d1.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    // LHC10b - LHC10c
    runLow[0] = 114650;
    runUp[0] = 121526;
  }
  else if (period.Contains("LHC10F6A") == 1 && isMC) {
    period  = "LHC10F6A";
    mapTitle = Form("LHC10f6a.pass%d", pass);
    
    nRunRanges = 3;
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    // LHC10c - LHC10g
    runLow[0] = 121527;
    runUp[0] = 136781;
    
    // LHC11a_7TeV
    runLow[1] = 139847;
    runUp[1] = 146631;
    
    // NO LHC11a_2.76TeV (-> LHC11b10a)
    
    // LHC11b (146975-150721) + LHC11c (150722-155837) + LHC11d(155838-159649)
    runLow[2] = 146975;
    runUp[2] = 159649;
  }
  else if (period.Contains("LHC11B2") == 1 && isMC) {
    period = "LHC11B2";
    mapTitle = Form("LHC11b2.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    // LHC10b - LHC10e
    runLow[0] = 114650;
    runUp[0] = 130930;
  }
  else if (period.Contains("LHC11B10A") == 1 && isMC) {
    period = "LHC11B10A";
    mapTitle = Form("LHC11b10a.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    // LHC11a, 2.76TeV
    runLow[0] = 146632;
    runUp[0] = 146974;
  }
  else if (period.Contains("LHC12F1") == 1 && isMC) {
    period = "LHC12F1";
    mapTitle = Form("LHC12f1.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    // LHC11a, 2.76TeV
    runLow[0] = 146632;
    runUp[0] = 146974;
  }
  //TODO Hope: Same map for 11A10A and 11A10B! Just use LHC11A10
  else if (period.Contains("LHC11A10") == 1 && isMC) {
    //LHC10h -> 2.76 ATeV (PbPb)
    //TODO maybe also valid MC for LHC11h
    period = "LHC11A10";
    mapTitle = Form("LHC11a10.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    // LHC10h
    runLow[0] = 136782;
    runUp[0] = 139846;
  }
  /*
  else if (period.Contains("LHC11A10A") == 1 && isMC) {
    //LHC10h -> 2.76 ATeV (PbPb)
    period = "LHC11A10A";
    mapTitle = Form("LHC11a10a.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    // LHC10h
    runLow[0] = 136782;
    runUp[0] = 139846;
  }
  /*
  else if (period.Contains("LHC11A10B") == 1 && isMC) {
    //LHC10h -> 2.76 ATeV (PbPb)
    period = "LHC11A10B";
    mapTitle = Form("LHC11a10b.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    // LHC10h
    runLow[0] = 136782;
    runUp[0] = 139846;
  }*/
  else if (period.Contains("LHC12A") == 1 && !isMC) {
    //LHC12a -> 8 TeV (pp)
    mapTitle = Form("LHC12a.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 170719;
    runUp[0] = 177311;
  }
  else if (period.Contains("LHC12B") == 1 && !isMC) {
    //LHC12b -> 8 TeV (pp)
    mapTitle = Form("LHC12b.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 177312;
    runUp[0] = 179356;
  }
  else if (period.Contains("LHC12C") == 1 && !isMC) {
    //LHC12c -> 8 TeV (pp)
    mapTitle = Form("LHC12c.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 179357;
    runUp[0] = 183173;
  }
  else if (period.Contains("LHC12D") == 1 && !isMC) {
    //LHC12d -> 8 TeV (pp)
    mapTitle = Form("LHC12d.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 183174;
    runUp[0] = 186345;
  }
  else if (period.Contains("LHC12E") == 1 && !isMC) {
    //LHC12e -> 8 TeV (pp)
    mapTitle = Form("LHC12e.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 186346;
    runUp[0] = 186635;
  }
  else if (period.Contains("LHC12F") == 1 && !isMC) {
    //LHC12f -> 8 TeV (pp)
    mapTitle = Form("LHC12f.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 186636;
    runUp[0] = 188166;
  }
  else if (period.Contains("LHC12G_PPB") == 1 && isMC) {
    //LHC12g (same map for 12g1, 12g4, ...) -> 5.023 ATeV (pPb)
    period = "LHC12G";
    mapTitle = Form("LHC12g_pPb.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 188167;
    runUp[0] = 188418;
  }
  else if (period.Contains("LHC12G_PP") == 1 && !isMC) {
    //LHC12g -> 8 TeV (pp)
    mapTitle = Form("LHC12g_pp.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 188419;
    runUp[0] = 188719;
  }
  else if (period.Contains("LHC12H") == 1 && !isMC) {
    //LHC12h -> 8 TeV (pp)
    mapTitle = Form("LHC12h.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 188720;
    runUp[0] = 192738;
  }
  else if (period.Contains("LHC12I") == 1 && !isMC) {
    //LHC12i and LHC12j -> 8 TeV (pp)
    mapTitle = Form("LHC12i.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 192739;
    //runUp[0] = 193766;
    //runLow[1] = 193767;
    runUp[0] = 194479;
  }
  else if (period.Contains("LHC13B2_FIXN1") == 1 && isMC) {
    //LHC13b2_fixn1 -> 5.023 ATeV (pPb)
    period = "LHC13B2_FIXn1";
    mapTitle = Form("LHC13b2_fixn1.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    /*
     * runLow[0] = 194480;
    runUp[0] = 196345;*/
  }
  else if (period.Contains("LHC13B2_FIX") == 1 && isMC) {
    //LHC13b2_fix -> 5.023 ATeV (pPb)
    period = "LHC13B2_FIX";
    mapTitle = Form("LHC13b2_fix.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    /*
     * runLow[0] = 194480;
    runUp[0] = 196345;*/
  }
  else if (period.Contains("LHC13B") == 1 && !isMC) {
    //LHC13a-d (beam only since 13b) periods -> 5.023 ATeV (pPb)
    // -> High luminosity periods (e and following) require
    // at the moment separate treatmeant
    period = "LHC13B";
    mapTitle = Form("LHC13b.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 194480;
    runUp[0] = 195874;
  }
  else if (period.Contains("LHC13F") == 1 && !isMC) {
    //LHC13e-f periods -> 5.023 ATeV (pPb) high luminosity
    period = "LHC13F";
    mapTitle = Form("LHC13f.pass%d", pass);
    
    runLow = new Int_t[nRunRanges];
    runUp = new Int_t[nRunRanges];
    
    runLow[0] = 195875;
    runUp[0] = 197411;
  }
  else if (period.Contains("DEFAULT")) {
    mapTitle = Form("Default_%s_pass%d", dataType.Data(), pass);
    addDefault = kTRUE;
    printf("Request for default object...\n");
#ifdef TEST
  runLow = new Int_t[nRunRanges];
  runUp = new Int_t[nRunRanges];
  runLow[0] = 333333;
  runUp[0] = 333333;
#endif
  }
  else {
    printf("Unknown period!\n");
    return -1;
  }
  

  
  
#ifdef TEST
  // TEST created file
  printf("\n\n***********TEST MODE ENABLED*********\n\n");
  Int_t fRun = (runUp[0] + runLow[0]) / 2.0;
 
  TH2D* etaMap = 0x0;
  TH2D* etaSigmaPar1Map = 0x0;
  
  Double_t sigmaPar0 = -1.0;
  
  // Load the eta correction maps
  AliOADBContainer etaMapsCont2(Form("TPCetaMaps_%s_pass%d", dataType.Data(), pass)); 
    
  Int_t statusCont2 = etaMapsCont2.InitFromFile("etaMaps/TPCetaMaps.root", Form("TPCetaMaps_%s_pass%d", dataType.Data(), pass));
  if (statusCont2) {
    printf("Failed initializing TPC eta correction maps from OADB");
    return -1;
  }
  else {
    if (isMC) {
      TString searchMap = Form("TPCetaMaps_%s_%s_pass%d", dataType.Data(), period.Data(), pass);
      etaMap = dynamic_cast<TH2D *>(etaMapsCont2.GetDefaultObject(searchMap.Data()));
      if (!etaMap) {
        // Try default object
        etaMap = dynamic_cast<TH2D *>(etaMapsCont2.GetDefaultObject(defaultObj.Data()));
        
        if (!etaMap) {
          printf("TPC eta correction map not found for period %s and also no default map found -> Disabled eta correction!!!",
                 period.Data());
          return -1;
        }
      }
    }
    else {
      etaMap = dynamic_cast<TH2D *>(etaMapsCont2.GetObject(fRun, defaultObj.Data()));
      if (!etaMap) {
        printf("TPC eta correction map not found for run %d -> Disabled eta correction!!!", fRun);
        return -1;
      }
    }
  }
  
  // Load the sigma parametrisation (1/dEdx vs tanTheta_local (~eta))
  AliOADBContainer etaSigmaMapsCont2(Form("TPCetaSigmaMaps_%s_pass%d", dataType.Data(), pass)); 
  
  statusCont2 = etaSigmaMapsCont2.InitFromFile("etaMaps/TPCetaMaps.root", Form("TPCetaSigmaMaps_%s_pass%d", dataType.Data(), pass));
  if (statusCont2) {
    printf("Failed initializing TPC eta sigma maps from OADB");
    return -1;
  }
  else {
    TObjArray* etaSigmaPars = 0x0;
    if (isMC) {
      TString searchMap = Form("TPCetaSigmaMaps_%s_%s_pass%d", dataType.Data(), period.Data(), pass);
      etaSigmaPars = dynamic_cast<TObjArray *>(etaSigmaMapsCont2.GetDefaultObject(searchMap.Data()));
      if (!etaSigmaPars) {
        // Try default object
        etaSigmaPars = dynamic_cast<TObjArray *>(etaSigmaMapsCont2.GetDefaultObject(defaultObj.Data()));
        if (!etaSigmaPars) {
          printf("TPC eta sigma parametrisation not found for period %s and also no default parametrisation found -> Using old sigma parametrisation!!!",
                 period.Data());
          return -1;
        }
      }
    }
    else {
      etaSigmaPars = dynamic_cast<TObjArray *>(etaSigmaMapsCont2.GetObject(fRun, defaultObj.Data()));
      if (!etaSigmaPars) {
        printf("TPC eta sigma parametrisation not found for run %d -> Using old sigma parametrisation!!!", fRun);
        return -1;
      }
    }
    
    etaSigmaPar1Map = dynamic_cast<TH2D *>(etaSigmaPars->FindObject("sigmaPar1Map"));
    TNamed* sigmaPar0Info = dynamic_cast<TNamed *>(etaSigmaPars->FindObject("sigmaPar0"));
    
    if (sigmaPar0Info) {
      TString sigmaPar0String = sigmaPar0Info->GetTitle();
      sigmaPar0 = sigmaPar0String.Atof();
    }
    else {
      // Something is weired because the object for parameter 0 could not be loaded -> New sigma parametrisation can not be used!
      etaSigmaPar1Map = 0x0;
      return -1;
    }
  }
  
  printf("\n\nLoaded c0: %.4f\n", sigmaPar0);
  TCanvas* c = new TCanvas();
  c->Divide(2,1);
  c->cd(1);
  etaMap->Draw("colz");
  c->cd(2);
  etaSigmaPar1Map->Draw("colz");
  
  return 0;
#endif
  
  
  
  
  
  
  
  
  // Open the map package and retrieve the OADBContainers. If the map package does not exist,
  // create it. If there is no OADBContainer yet, create a new one
  TString filePathNameMapPackage = gSystem->ExpandPathName(Form("%s/%s", pathMapPackage.Data(), fileNameMapPackage.Data()));
 
  AliOADBContainer* etaMapsCont = 0x0;
  
  if (!isSigma) {
    etaMapsCont = new AliOADBContainer(Form("TPCetaMaps_%s_pass%d", dataType.Data(), pass)); 
    Int_t statusCont = etaMapsCont->InitFromFile(filePathNameMapPackage.Data(), Form("TPCetaMaps_%s_pass%d", dataType.Data(), pass));
    if (statusCont) 
      printf("No OADBContainer for the current settings found - creating a new one...\n");
  }
  else {
    etaMapsCont = new AliOADBContainer(Form("TPCetaSigmaMaps_%s_pass%d", dataType.Data(), pass)); 
    Int_t statusCont = etaMapsCont->InitFromFile(filePathNameMapPackage.Data(), Form("TPCetaSigmaMaps_%s_pass%d", dataType.Data(), pass));
    if (statusCont) 
      printf("No OADBContainer for the current settings found - creating a new one...\n");
  }


  // Open the map that should be added
  filePathNameMapToAdd = gSystem->ExpandPathName(filePathNameMapToAdd.Data());
  
  TFile* fMapToAdd = TFile::Open(filePathNameMapToAdd.Data());
  
  if (!fMapToAdd) {
    std::cout << "Failed to open map file \"" << filePathNameMapToAdd.Data() << "\"!" << std::endl;
    return -1;
  }
    
  TH2D* hMap = dynamic_cast<TH2D*>(fMapToAdd->Get(mapNameInFile.Data()));
  if (!hMap) {
    std::cout << "Failed to load map!" << std::endl;
  
    return -1;
  }
  
  hMap->SetName(mapName.Data());
  hMap->SetTitle(mapTitle.Data());
  
  TNamed* c0Info = 0x0;
  TObjArray* sigmaPars = 0x0;
  
  if (isSigma) {
    c0Info = dynamic_cast<TNamed*>(fMapToAdd->Get(par0StringName.Data()));
    if (!c0Info) {
      std::cout << "Failed to load parameter 0!" << std::endl;
      return -1;
    }
    
    c0Info->SetName("sigmaPar0");
    
    sigmaPars = new TObjArray(2);
    sigmaPars->SetOwner(kTRUE);
    
    if (isMC) {
      sigmaPars->SetName(Form("TPCetaSigmaMaps_%s_%s_pass%d", dataType.Data(), period.Data(), pass));
    }
    else {
      sigmaPars->SetName(Form("TPCetaSigmaMaps_%s_pass%d", dataType.Data(), pass));
    }
    
    sigmaPars->Add(hMap);
    sigmaPars->Add(c0Info);
    
    if (addDefault) {
      sigmaPars->SetName(defaultObj.Data());
      if (isMC) {
        TObjArray* oldObj = (TObjArray*)etaMapsCont->GetDefaultList()->Remove(etaMapsCont->GetDefaultObject(sigmaPars->GetName()));
        if (oldObj)
          oldObj->Delete();
        delete oldObj;
      }
      else {
        etaMapsCont->CleanDefaultList();
      }
      etaMapsCont->AddDefaultObject(sigmaPars);
    }
    else {
      if (isMC) {
        TObjArray* oldObj = (TObjArray*)etaMapsCont->GetDefaultList()->Remove(etaMapsCont->GetDefaultObject(sigmaPars->GetName()));
        if (oldObj) {
          printf("Updating existing object \"%s\"...\n", mapTitle.Data());
          oldObj->Delete();
          delete oldObj;
        }
        else {
          printf("Creating new object \"%s\"...\n", mapTitle.Data());
        }
        etaMapsCont->AddDefaultObject(sigmaPars);
      }
      else {
        // Check, if the object already exists by taking the center of the run range.
        // If there is a conflict for the range, AliOADBContainer will give an AliFatal.
        for (Int_t range = 0; range < nRunRanges; range++) {
          Int_t index = etaMapsCont->GetIndexForRun((runUp[range] + runLow[range]) / 2.0);
          if (index < 0) {
            printf("Creating new object for run range %d - %d...\n", runLow[range], runUp[range]);
            etaMapsCont->AppendObject((range == 0) ? sigmaPars : sigmaPars->Clone(), runLow[range], runUp[range]);
          }
          else {
            printf("Updating existing object for run range %d - %d...\n", runLow[range], runUp[range]);
            etaMapsCont->UpdateObject(index, (range == 0) ? sigmaPars : sigmaPars->Clone(), runLow[range], runUp[range]);
          }
        }
      }
    }
  }
  else {
    
    if (isMC) {
      hMap->SetName(Form("TPCetaMaps_%s_%s_pass%d", dataType.Data(), period.Data(), pass));
    }
    
    if (addDefault) {
      hMap->SetName(defaultObj.Data());
      
      if (isMC) {
        TH2D* oldMap = (TH2D*)etaMapsCont->GetDefaultList()->Remove(etaMapsCont->GetDefaultObject(hMap->GetName()));
        delete oldMap;
      }
      else {
        etaMapsCont->CleanDefaultList();
      }
      etaMapsCont->AddDefaultObject(hMap);
    }
    else {
      if (isMC) {
        TH2D* oldMap = (TH2D*)etaMapsCont->GetDefaultList()->Remove(etaMapsCont->GetDefaultObject(hMap->GetName()));
        if (oldMap) {
          printf("Updating existing object \"%s\"...\n", mapTitle.Data());
          delete oldMap;
        }
        else {
          printf("Creating new object \"%s\"...\n", mapTitle.Data());
        }
        etaMapsCont->AddDefaultObject(hMap);
      }
      else {
        // Check, if the object already exists by taking the center of the run range.
        // If there is a conflict for the range, AliOADBContainer will give an AliFatal.
        for (Int_t range = 0; range < nRunRanges; range++) {
          Int_t index = etaMapsCont->GetIndexForRun((runUp[range] + runLow[range]) / 2.0);
          if (index < 0) {
            printf("Creating new object for run range %d - %d...\n", runLow[range], runUp[range]);
            etaMapsCont->AppendObject((range == 0) ? hMap : hMap->Clone(), runLow[range], runUp[range]);
          }
          else {
            printf("Updating existing object for run range %d - %d...\n", runLow[range], runUp[range]);
            etaMapsCont->UpdateObject(index, (range == 0) ? hMap : hMap->Clone(), runLow[range], runUp[range]);
          }
        }
      }
    }
  }

  //etaMapsCont->WriteToFile(filePathNameMapPackage.Data());
  // Does the same as the WriteToFile function, but allow for using "kOverwrite" to save memory!
  TFile* f = new TFile(filePathNameMapPackage.Data(), "update");
  if (!f) {
    printf("Error opening file \"%s\"!\n", filePathNameMapPackage.Data());
    return -1;
  }
  f->Delete(etaMapsCont->GetName());
  etaMapsCont->Write(0, TObject::kOverwrite);
  f->Purge();
  f->Close();
  
  
  gSystem->Exec(Form("echo \"%s%s:\t\t%s#%s\n\" >> %s/UsedFilesForMap.txt", mapTitle.Data(), isSigma ? " (sigma)" : " (etaCorr)",
                     filePathNameMapToAdd.Data(), mapNameInFile.Data(), pathMapPackage.Data()));
  
  fMapToAdd->Close();
  
  // I don't understand why, but if the file list is not cleaned up, there are problems when adding further things to the OADB containers
  gROOT->GetListOfFiles()->RemoveAll();
  
  printf("\n****************Successfully added map \"%s\"!\n\n\n\n", mapTitle.Data());
  return 0;
}