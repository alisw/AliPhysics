//Creates from the full MC result the correction factors used to compare the fast MC simulation with the full MC simulation

#include "TH3.h"
#include "./calcEfficiency.C"

Int_t createFileForBbBCorrections(TString pathNameEfficiency, TString outfileName, TString trackObservables = "", TString jetPtLimitsString = "") {
  TFile* fileEff = TFile::Open(pathNameEfficiency.Data());
  if (!fileEff) {
    printf("Failed to open efficiency file \"%s\"\n", pathNameEfficiency.Data());
    return -1;
  }
  
  AliCFContainer* data = (AliCFContainer*)(fileEff->Get("containerEff"));
  if (!data) {
    printf("Failed to load efficiency container!\n");
    return -1;
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
  
  //Variables all defined in calcEfficiency.C
  iPt     = data->GetVar(Form("%s_{T} (GeV/c)", momentumString.Data()));
  iMCid   = data->GetVar("MC ID");
  iEta    = data->GetVar("#eta");
  iCharge = data->GetVar("Charge (e_{0})");
  iMult   = data->GetVar("Centrality Percentile");
  iJetPt = data->GetVar(Form("%s_{T}^{jet} (GeV/c)", momentumString.Data()));
  if (iJetPt == -1) {
    std::cout << "Has to be Jet Correction File" << std::endl;
    return -1;
  }
  iZ     = data->GetVar(Form("z = %s_{T}^{track} / %s_{T}^{jet}", momentumString.Data(), momentumString.Data()));
  iXi    = data->GetVar(Form("#xi = ln(%s_{T}^{jet} / %s_{T}^{track})", momentumString.Data(), momentumString.Data()));
  iDistance = data->GetVar("R");
  ijT = data->GetVar("j_{T} (GeV/c)");
  
  const Int_t nOfTrackObservables = 5;
  Int_t trackObservableBins[nOfTrackObservables] = {iPt, iZ, iXi, iDistance, ijT};
  TString observableNames[nOfTrackObservables] = {"TrackPt", "Z", "Xi", "R", "jT"};  
  
  //Setting jet limits and getting the number of rec/gen jets from the file associated with the efficiency file
  const Int_t nOfJetBins = 5;
  Double_t jetPtLimits[2*nOfJetBins] = {5.0,10.0,10.0,15.0,15.0,20.0,20.0,30.0,30.0,80.0};
  Int_t nOfJets[2][nOfJetBins] = {0};
  Int_t jetBinLimits[2*nOfJetBins] = {0};
  
  TString pathNameDataMC = pathNameEfficiency;
  pathNameDataMC.ReplaceAll("_efficiency", "");
  
  TFile* fDataMC = TFile::Open(pathNameDataMC.Data());
  if (!fDataMC)  {
    std::cout << std::endl;
    std::cout << "Failed to open file \"" << pathNameDataMC.Data() << "\" to obtain num of rec/gen jets!" << std::endl;   
    return -1;
  }
  
  TString listName = pathNameDataMC;
  listName.Replace(0, listName.Last('/') + 1, "");
  listName.ReplaceAll(".root", "");
    
  TObjArray* histList = (TObjArray*)(fDataMC->Get(listName.Data()));
  if (!histList) {
    std::cout << std::endl;
    std::cout << "Failed to load list \"" << listName.Data() << "\" to obtain num of rec/gen jets!" << std::endl;
    return -1;
  }    
  
  TH2* hNjetsGen = (TH2D*)histList->FindObject("fh2FFJetPtGen");
  TH2* hNjetsRec = (TH2D*)histList->FindObject("fh2FFJetPtRec");
  
   // If desired, restrict centrality axis
  Int_t lowerCentrality = 0;
  Int_t upperCentrality = 100;
  Int_t lowerCentralityBinLimit = -1;
  Int_t upperCentralityBinLimit = -2; // Integral(lowerCentBinLimit, uppCentBinLimit) will not be restricted if these values are kept. In particular, under- and overflow bin will be used!
  
  if (lowerCentrality >= -1 && upperCentrality >= -1) {
    // Add subtract a very small number to avoid problems with values right on the border between two bins
    lowerCentralityBinLimit = data->GetAxis(iMult, 0)->FindFixBin(lowerCentrality + 0.001);
    upperCentralityBinLimit = data->GetAxis(iMult, 0)->FindFixBin(upperCentrality - 0.001);
  }

  data->SetRangeUser(iMult, lowerCentralityBinLimit, upperCentralityBinLimit, kTRUE);
  
  for (Int_t i=0;i<nOfJetBins;++i) { 
    jetBinLimits[2*i] = hNjetsRec->GetYaxis()->FindFixBin(jetPtLimits[2*i] + 0.001);
    jetBinLimits[2*i+1] = hNjetsRec->GetYaxis()->FindFixBin(jetPtLimits[2*i+1] - 0.001);  

    nOfJets[0][i] = hNjetsGen ? hNjetsGen->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, jetBinLimits[2*i], jetBinLimits[2*i+1]) : 1.;
    nOfJets[1][i] = hNjetsRec ? hNjetsRec->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, jetBinLimits[2*i], jetBinLimits[2*i+1]) : 1.;

    jetBinLimits[2*i] = data->GetAxis(iJetPt, 0)->FindFixBin(jetPtLimits[2*i] + 0.001);
    jetBinLimits[2*i+1] = data->GetAxis(iJetPt, 0)->FindFixBin(jetPtLimits[2*i+1] - 0.001);      
  }
  
//   geantFlukaCorrection(data, kStepGenWithGenCuts, kStepRecWithRecCutsMeasuredObsPrimaries, kFALSE); //No need for geant-Fluka correction in DPMJet - include into calcEfficiency
  
  TFile* outFile = new TFile(outfileName.Data(),"RECREATE");
  
  //Do for all species
  for (Int_t jetPtStep = 0;jetPtStep<nOfJetBins;++jetPtStep) {
    AliCFContainer *dataRebinned = new AliCFContainer(*data);
    dataRebinned->SetRangeUser(iJetPt,jetBinLimits[2*jetPtStep],jetBinLimits[2*jetPtStep+1],kTRUE);
    Double_t factor_Numerator[2] = { nOfJets[1][jetPtStep] > 0 ? 1. / nOfJets[1][jetPtStep] : 0., 0.  };
    Double_t factor_Denominator[2] = { nOfJets[0][jetPtStep] > 0 ? 1. / nOfJets[0][jetPtStep] : 0., 0.  };
    AliCFEffGrid* eff = new AliCFEffGrid("eff", "Efficiency x Acceptance x pT Resolution", *dataRebinned);
    eff->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, kStepGenWithGenCuts);    
    eff->GetNum()->Scale(factor_Numerator);
    eff->GetDen()->Scale(factor_Denominator);
    
    for (Int_t observable = 0;observable<nOfTrackObservables;++observable) {
      TH1* hEff = eff->Project(trackObservableBins[observable]);
      hEff->SetNameTitle(TString::Format("hBbBCorr%s_%02d_%02d",observableNames[observable].Data(),(Int_t)jetPtLimits[jetPtStep*2],(Int_t)jetPtLimits[jetPtStep*2+1]),"");  
      for (Int_t bin=1;bin<=hEff->GetNbinsX();++bin) {
        const Double_t eff = hEff->GetBinContent(bin);
        const Double_t effErr = hEff->GetBinError(bin);
        
        const Double_t effBbB = eff > 0. ? 1. / eff : 0.;
        const Double_t effBbBErr = eff > 0. ? effErr / (eff * eff) : 0.;
        
        hEff->SetBinContent(bin, effBbB);
        hEff->SetBinError(bin, effBbBErr);        
      }
      hEff->Write();
      delete hEff;
    }
    delete eff;
    delete dataRebinned;    
  }  
    
  //Loop through species  
  for (Int_t species=0;species<AliPID::kSPECIES;++species) {
    TString speciesString = TString("_") + TString(AliPID::ParticleShortName(species));
    for (Int_t jetPtStep = 0;jetPtStep<nOfJetBins;++jetPtStep) {
      AliCFContainer *dataRebinned = new AliCFContainer(*data);
      dataRebinned->SetRangeUser(iJetPt,jetBinLimits[2*jetPtStep],jetBinLimits[2*jetPtStep+1],kTRUE);
      dataRebinned->SetRangeUser(iMCid,species+1,species+1,kTRUE);
      Double_t factor_Numerator[2] = { nOfJets[1][jetPtStep] > 0 ? 1. / nOfJets[1][jetPtStep] : 0., 0.  };
      Double_t factor_Denominator[2] = { nOfJets[0][jetPtStep] > 0 ? 1. / nOfJets[0][jetPtStep] : 0., 0.  };
      AliCFEffGrid* eff = new AliCFEffGrid("eff", "Efficiency x Acceptance x pT Resolution", *dataRebinned);
      eff->CalculateEfficiency(kStepRecWithRecCutsMeasuredObsPrimaries, kStepGenWithGenCuts);
      eff->GetNum()->Scale(factor_Numerator);
      eff->GetDen()->Scale(factor_Denominator);
      for (Int_t observable = 0;observable<nOfTrackObservables;++observable) {        
        TH1* hEff = (TH1*)eff->Project(trackObservableBins[observable]);  
        TString hName = TString::Format("hBbBCorr%s%s_%02d_%02d",observableNames[observable].Data(),speciesString.Data(),(Int_t)jetPtLimits[jetPtStep*2], (Int_t)jetPtLimits[jetPtStep*2+1]);   
        hEff->SetNameTitle(hName.Data(),"");  
        for (Int_t bin=1;bin<=hEff->GetNbinsX();++bin) {
          const Double_t eff = hEff->GetBinContent(bin);
          const Double_t effErr = hEff->GetBinError(bin);
          
          const Double_t effBbB = eff > 0. ? 1. / eff : 0.;
          const Double_t effBbBErr = eff > 0. ? effErr / (eff * eff) : 0.;
          
          hEff->SetBinContent(bin, effBbB);
          hEff->SetBinError(bin, effBbBErr);        
        }      
        hEff->Write();
        delete hEff;
      }
      delete eff;
      delete dataRebinned;
    }
  }  
  
  
  outFile->Close();
  return 0;
  
}
