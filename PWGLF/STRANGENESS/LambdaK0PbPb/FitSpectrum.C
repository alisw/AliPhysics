void FitSpectrum(const char* filename, const char * listName = "clambdak0Histo_00", const char* suffix = "test",Int_t ihist = 0, Int_t iparticle) {
  
  // load basic libs, needed to 
  //gROOT->LoadMacro("run.C");
  gROOT->LoadMacro("runQA.C");
  InitAndLoadLibs();

  // Load Lee's Macro
  gROOT->LoadMacro("AliMassFitControl.h+g");
  gROOT->LoadMacro("PtMassAna2.C");
  gROOT->LoadMacro("MultYields2QA.C");

  char* histName = 0;
  switch (iparticle) {
  case 0:
    switch (ihist){
          case 0:
    		histName = "h2PtVsMassK0";
                break;
          case 1:
		histName = "h2DcaPosToPrimVertexK0vsMassK0";
		break;
          case 2:
		histName = "h2DcaNegToPrimVertexK0vsMassK0";
		break;
          case 3:
		histName = "h2RadiusV0K0vsMassK0";
		break;
          case 4:
		histName = "h2DecayLengthV0K0vsMassK0";   
		break;
          case 5:
		histName = "h2DcaV0DaughtersK0vsMassK0";
		break;
          case 6:
		histName = "h2CosPointAngleK0vsMassK0";
		break;
          default:
                cout << "FitSpectrum - unrecognised histogram, default is standard pt" << endl;
                histName = "h2PtVsMassK0";
          }
    break;
  case 1:
    switch (ihist){
          case 0:
		histName = "h2PtVsMassLambda";   
		break;
          case 1:
		histName = "h2DcaPosToPrimVertexLvsMassL";   
		break;
          case 2:
		histName = "h2DcaNegToPrimVertexLvsMassL";   
		break;
          case 3:
		histName = "h2RadiusV0LvsMassL";   
		break;
          case 4:
		histName = "h2DecayLengthV0LvsMassL";   
		break;
          case 5:
		histName = "h2DcaV0DaughtersLvsMassL";
		break;
          case 6:
		histName = "h2CosPointAngleLvsMassL";
		break;
          default:
                cout << "FitSpectrum - unrecognised histogram, default is standard pt" << endl;
                histName = "h2PtVsMassLambda";   
          }
    break;
    case 2:
      histName = "h2PtVsMassAntiLambda";
      break;
    case 3:
      // This is special case because we have to find two histograms and add them.
      // It is dealt with separately below but we still set histname because it is appended
      // to file names for control plots
      histName = "h2PtVsMassLLbarSummed";
      break;
    default:
      cout << "Particle "<< iparticle << " yet to be implemented" << endl;
      return;
  }
  cout << "FitSpectrum - histogram " << histName << " used" << endl;

  
  TFile *file = new TFile(filename);
  TList *list = file->Get(listName); 
  TH2F* h2;
  if (iparticle == 3) { // Special case of combined Lambda + anti-Lambda
    TH2F* hLam = (TH2F*)list->FindObject("h2PtVsMassLambda");
    h2= (TH2F*)list->FindObject("h2PtVsMassAntiLambda");
    h2->Add(hLam);
  } else {
    h2 = (TH2F*) list->FindObject(histName);
  }

  h2->Draw();
  
  ///// iNorm is used by MultYields2QA to normalize the distributions
  TH1F * h1 = (TH1F*) list->FindObject("h1PrimaryVertexZ");
  Int_t iNorm = h1->GetEntries();
  cout << "number of entries with Zvertex < |10cm|: " << iNorm;

  TString suffixFull = histName;
  if(strlen(suffix)) suffixFull = suffixFull + "_" + suffix;
  //MultYields3((TH3F*)h2,iparticle,0,suffixFull); // FIXME: modify MultYields2 to handle 1D histos
  MultYields2QA(h2,iparticle,ihist,iNorm,0,suffixFull); // FIXME: modify MultYields2 to handle 1D histos

  
}
