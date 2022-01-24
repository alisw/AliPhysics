#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"

#include "TF1.h"
#include "./PieceWisePoly.cxx"

#include "AliCFContainer.h"
#include "AliCFEffGrid.h"

#include "THnSparseDefinitions.h"

#include <iostream>

using namespace std;

TString getStringFromTObjStrArray(TObjArray *arr, Int_t position)
{
  TObjString* objStr = (TObjString*)(arr->At(position));
  if (!objStr)
    return "";

  return objStr->GetString();
}

const Int_t speciescolors[5] = {kMagenta, kYellow, kRed, kGreen, kBlue};

const Int_t smoothType = 2;     // (matching function and first derivative)

Int_t FitCorrFactors(TString effFile, TString outputfile, TString parameters = "") {
  TFile* fileEff = new TFile(effFile.Data());
  if (!fileEff) {
    printf("Failed to open efficiency file \"%s\"\n", effFile.Data());
    return -1;
  }
  
  AliCFContainer* data = (AliCFContainer*)(fileEff->Get("containerEff"));
  if (!data) {
    printf("Failed to load efficiency container!\n");
    return -1;
  }    

  Int_t parts[AliPID::kSPECIES][2];
  Double_t* cuts[AliPID::kSPECIES][2];
  Int_t* nparameters[AliPID::kSPECIES][2];
  
  TObjArray* parameterArray = parameters.Tokenize(";");
  for (Int_t species=0;species<AliPID::kSPECIES;species++) {
    for (Int_t charge=0;charge<2;charge++) {
      TString speciesString = getStringFromTObjStrArray(parameterArray, 2*species + charge);
      
      TObjArray* speciesParameterArray = speciesString.Tokenize("-");
      
      Int_t nSpeciesParameters = (speciesParameterArray->GetEntriesFast() + 1)/2;
      parts[species][charge] = nSpeciesParameters;
      cuts[species][charge] = new Double_t[nSpeciesParameters-1];
      nparameters[species][charge] = new Int_t[nSpeciesParameters];
      nparameters[species][charge][0] = getStringFromTObjStrArray(speciesParameterArray, 0).Atoi();
      
      for (Int_t cut=1;cut<nSpeciesParameters;cut++) {
        cuts[species][charge][cut-1] = getStringFromTObjStrArray(speciesParameterArray, 2*cut-1).Atof();
        nparameters[species][charge][cut] = getStringFromTObjStrArray(speciesParameterArray, 2*cut).Atoi();
      }
      
      delete speciesParameterArray;
      speciesParameterArray = 0x0;
    }
  }
  delete parameterArray;
  parameterArray = 0x0;

  //For negative electrons
  Int_t species = AliPID::kElectron;
  Int_t charge = 0;
  parts[species][charge] = 4;
  
  cuts[species][charge] = new Double_t[parts[species][charge]-1];
  cuts[species][charge][0] = 0.6;
  cuts[species][charge][1] = 3.2;
  cuts[species][charge][2] = 8.0;
  
  nparameters[species][charge] = new Int_t[parts[species][charge]];
  nparameters[species][charge][0] = 7;
  nparameters[species][charge][1] = 5;
  nparameters[species][charge][2] = 3;
  nparameters[species][charge][3] = 2;
  
  //For positive electrons
  species = AliPID::kElectron;
  charge = 1;
  parts[species][charge] = 4;
  
  cuts[species][charge] = new Double_t[parts[species][charge]-1];
  cuts[species][charge][0] = 0.6;
  cuts[species][charge][1] = 3.2;
  cuts[species][charge][2] = 8.0;
  
  nparameters[species][charge] = new Int_t[parts[species][charge]];
  nparameters[species][charge][0] = 7;
  nparameters[species][charge][1] = 5;
  nparameters[species][charge][2] = 3;
  nparameters[species][charge][3] = 2; 
    
  //For negative muons
  species = AliPID::kMuon;
  charge = 0;
  parts[species][charge] = 6;
  
  cuts[species][charge] = new Double_t[parts[species][charge]-1];
  cuts[species][charge][0] = 0.8;
  cuts[species][charge][1] = 1.6;
  cuts[species][charge][2] = 3.0;
  cuts[species][charge][3] = 10.0;
  cuts[species][charge][4] = 12.0;
  
  nparameters[species][charge] = new Int_t[parts[species][charge]];
  nparameters[species][charge][0] = 9;
  nparameters[species][charge][1] = 4;
  nparameters[species][charge][2] = 4;
  nparameters[species][charge][3] = 3;
  nparameters[species][charge][4] = 3;
  nparameters[species][charge][5] = 2;
   
  //For positive muons
  species = AliPID::kMuon;
  charge = 1;
  parts[species][charge] = 6;
  
  cuts[species][charge] = new Double_t[parts[species][charge]-1];
  cuts[species][charge][0] = 0.8;
  cuts[species][charge][1] = 1.6;
  cuts[species][charge][2] = 3.0;
  cuts[species][charge][3] = 10.0;
  cuts[species][charge][4] = 12.0;
  
  nparameters[species][charge] = new Int_t[parts[species][charge]];
  nparameters[species][charge][0] = 9;
  nparameters[species][charge][1] = 4;
  nparameters[species][charge][2] = 4;
  nparameters[species][charge][3] = 3;
  nparameters[species][charge][4] = 3;
  nparameters[species][charge][5] = 2;
   
    //For negative pions
  species = AliPID::kPion;
  charge = 0;
  parts[species][charge] = 6;
  
  cuts[species][charge] = new Double_t[parts[species][charge]-1];
  cuts[species][charge][0] = 0.8;
  cuts[species][charge][1] = 1.6;
  cuts[species][charge][2] = 3.0;
  cuts[species][charge][3] = 10.0;
  cuts[species][charge][4] = 12.0;
  
  nparameters[species][charge] = new Int_t[parts[species][charge]];
  nparameters[species][charge][0] = 9;
  nparameters[species][charge][1] = 4;
  nparameters[species][charge][2] = 4;
  nparameters[species][charge][3] = 3;
  nparameters[species][charge][4] = 3;
  nparameters[species][charge][5] = 2;
   
  //For positive pions
  species = AliPID::kPion;
  charge = 1;
  parts[species][charge] = 6;
  
  cuts[species][charge] = new Double_t[parts[species][charge]-1];
  cuts[species][charge][0] = 0.8;
  cuts[species][charge][1] = 1.6;
  cuts[species][charge][2] = 3.0;
  cuts[species][charge][3] = 10.0;
  cuts[species][charge][4] = 12.0;
  
  nparameters[species][charge] = new Int_t[parts[species][charge]];
  nparameters[species][charge][0] = 9;
  nparameters[species][charge][1] = 4;
  nparameters[species][charge][2] = 4;
  nparameters[species][charge][3] = 3;
  nparameters[species][charge][4] = 3;
  nparameters[species][charge][5] = 2;
    
    //For negative kaons
  species = AliPID::kKaon;
  charge = 0;
  parts[species][charge] = 5;
  
  cuts[species][charge] = new Double_t[parts[species][charge]-1];
  cuts[species][charge][0] = 0.4;
  cuts[species][charge][1] = 1.2;
  cuts[species][charge][2] = 6.0;
  cuts[species][charge][3] = 15.0;
  
  nparameters[species][charge] = new Int_t[parts[species][charge]];
  nparameters[species][charge][0] = 3;
  nparameters[species][charge][1] = 3;
  nparameters[species][charge][2] = 5;
  nparameters[species][charge][3] = 4;
  nparameters[species][charge][4] = 2;

  //For positive kaons
  species = AliPID::kKaon;
  charge = 1;
  parts[species][charge] = 5;
  
  cuts[species][charge] = new Double_t[parts[species][charge]-1];
  cuts[species][charge][0] = 0.4;
  cuts[species][charge][1] = 1.2;
  cuts[species][charge][2] = 6.0;
  cuts[species][charge][3] = 15.0;
  
  nparameters[species][charge] = new Int_t[parts[species][charge]];
  nparameters[species][charge][0] = 3;
  nparameters[species][charge][1] = 3;
  nparameters[species][charge][2] = 5;
  nparameters[species][charge][3] = 4;
  nparameters[species][charge][4] = 2;
  
    //For negative protons
  species = AliPID::kProton;
  charge = 0;
  parts[species][charge] = 6;
  
  cuts[species][charge] = new Double_t[parts[species][charge]-1];
  cuts[species][charge][0] = 0.4;
  cuts[species][charge][1] = 1.6;
  cuts[species][charge][2] = 2.5;
  cuts[species][charge][3] = 8.0;
  cuts[species][charge][4] = 12.0;
  
  nparameters[species][charge] = new Int_t[parts[species][charge]];
  nparameters[species][charge][0] = 6;
  nparameters[species][charge][1] = 4;
  nparameters[species][charge][2] = 4;
  nparameters[species][charge][3] = 2;
  nparameters[species][charge][4] = 5;
  nparameters[species][charge][5] = 2;
 
  //For positive protons
  species = AliPID::kProton;
  charge = 1;
  parts[species][charge] = 6;  
  
  cuts[species][charge] = new Double_t[parts[species][charge]-1];
  cuts[species][charge][0] = 0.4;
  cuts[species][charge][1] = 1.6;
  cuts[species][charge][2] = 2.5;
  cuts[species][charge][3] = 8.0;
  cuts[species][charge][4] = 12.0;
  
  nparameters[species][charge] = new Int_t[parts[species][charge]];
  nparameters[species][charge][0] = 6;
  nparameters[species][charge][1] = 4;
  nparameters[species][charge][2] = 4;
  nparameters[species][charge][3] = 2;
  nparameters[species][charge][4] = 5;
  nparameters[species][charge][5] = 2;
    
  PieceWisePoly* polynoms[AliPID::kSPECIES][2];
  
  for (Int_t species=0;species < AliPID::kSPECIES;species++) {
    for (Int_t charge=0;charge<2;charge++) {
      polynoms[species][charge] = new PieceWisePoly(parts[species][charge], cuts[species][charge], nparameters[species][charge], 0, 50, 0x0, smoothType);
    }
  }                                                          

  // For backward compatibility:
  // Check whether "P_{T}" or "p_{T}" is used
  TString momentumString = "p";
  for (Int_t i = 0; i < data->GetNVar(); i++) {
    TString temp = data->GetVarTitle(i);
    if (temp.Contains("P_{")) {
      momentumString = "P";
      break; 
    }
    else if (temp.Contains("p_{")) {
      momentumString = "p";
      break;
    }
  }
  
  Int_t iPt     = data->GetVar(Form("%s_{T} (GeV/c)", momentumString.Data()));
  Int_t iMCid   = data->GetVar("MC ID");
  Int_t iEta    = data->GetVar("#eta");
  Int_t iCharge = data->GetVar("Charge (e_{0})");
  Int_t iMult   = data->GetVar("Centrality Percentile");
  
  TFile* output = new TFile(outputfile.Data(),"RECREATE");
  TCanvas* c = new TCanvas("cEfficiencies","Efficiency x Acceptance");
  c->SetLogx(kTRUE);
  
  for (Int_t species=0;species<AliPID::kSPECIES;++species) {
    for (Int_t chargeBin=1;chargeBin<=2;chargeBin++) {
  
      data->SetRangeUser(iMCid,species+1,species+1,kTRUE);
      data->SetRangeUser(iCharge,chargeBin,chargeBin,kTRUE);
//       data->SetRangeUser(iCharge,-1,-1,kTRUE);
      AliCFEffGrid* effSingleTrack = new AliCFEffGrid("effSingleTrack", "Efficiency x Acceptance", *data);
      effSingleTrack->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, kStepGenWithGenCuts);      
      TH1* h = (TH1D*)effSingleTrack->Project(iPt);
      h->SetNameTitle(Form("hSingleTrackEfficiency_%s_%s", AliPID::ParticleShortName(species), chargeBin == 2 ? "pos" : "neg"),Form("%s for %s", AliPID::ParticleLatexName(species), chargeBin == 2 ? "pos" : "neg"));
      h->SetLineColor(speciescolors[species] + 2*chargeBin - 2);
      h->SetMarkerColor(speciescolors[species] + 2*chargeBin -2);
      h->GetXaxis()->SetRangeUser(0.15,50);
      for (Int_t i=1;i<=h->GetNbinsX();++i) {
        if (h->GetBinContent(i) == 1 && h->GetBinError(i) == 0) {
          h->SetBinContent(i,0);
        }
      }      
      TF1* func = new TF1(TString::Format("func_%s_%s",AliPID::ParticleShortName(species), chargeBin == 2 ? "pos" : "neg"),polynoms[species][chargeBin-1],0,50,polynoms[species][chargeBin-1]->GetNOfParam()); 
      func->SetLineColor(speciescolors[species] + 2*chargeBin - 2);

      h->Fit(func,"I","same");

      Int_t nOfParts = polynoms[species][chargeBin-1]->GetNParts();
      TString parameters = TString::Format("%d",nOfParts);
      Double_t* cuts = polynoms[species][chargeBin-1]->GetCuts();
      Int_t* nOfParameters = polynoms[species][chargeBin-1]->GetNParameters();
      for (Int_t part = 0;part<nOfParts - 1;++part) {
        parameters += TString::Format(";%f",cuts[part]);
      }
      for (Int_t part = 0;part<nOfParts;++part) {
        parameters += TString::Format(";%d",nOfParameters[part]);
      }      
      for (Int_t i=0;i<func->GetNpar();++i) {
        parameters += TString::Format(";%f",func->GetParameter(i));
      }
      TNamed* paramSave = new TNamed(TString::Format("fastSimulationParameters_%s_%s",AliPID::ParticleShortName(species),chargeBin == 2 ? "pos" : "neg").Data(),parameters.Data());
      paramSave->Write();
      h->Draw("same");
      h->Write();
      h = 0x0;
      delete effSingleTrack;
      effSingleTrack = 0x0;
    }
  }
  c->Write();
  output->Close();
  

return 0;
}

