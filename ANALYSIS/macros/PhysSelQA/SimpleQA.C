#include <Riostream.h>

// ROOT includes
#include "TFile.h"
#include "TString.h"
#include "TLine.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TMath.h"

/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// Simplified QA of PhysicsSelection (based on the previously Existing code).  //
// Please note that the event_stat.root files must be first downloaded locally //
// using the script CopyFilesToLocal.sh. This ensures the proper naming of the //
// files, namely RUNNUMBER_event_stat.root.                                    //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

void      SimpleQA();
void      MakeQAperPeriod(const Char_t * inputList = "inputListLHC13b.txt", const Char_t * outputFileName = "./Plots/2013/13bPass3.pdf", const Char_t * label = "LHC13bPass3");
Float_t   GetFraction(const Char_t * inputFile = "event_stat.root", const  Char_t * columnLabel = "Accepted",  UInt_t triggerBit = 2, Bool_t returnError = kFALSE);
TList   * InitializeHistograms(const Char_t * label);
TString * GetTriggerBitName(Int_t bitNumber = 0);
UInt_t    GetFillFromRunNumber(UInt_t runNumber);
void      AddFillSeparationLines(TH1D * hist);


//______________________________________________________________________________
void SimpleQA() {
  //
  // execution part -- make the QA for all periods
  //
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2012/inputListLHC12aPass1.txt", "./Plots/2012/12aPass1.pdf", "LHC12aPass1");
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2012/inputListLHC12bPass1.txt", "./Plots/2012/12bPass1.pdf", "LHC12bPass1");
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2012/inputListLHC12cPass1.txt", "./Plots/2012/12cPass1.pdf", "LHC12cPass1");
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2012/inputListLHC12dPass1.txt", "./Plots/2012/12dPass1.pdf", "LHC12dPass1");
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2012/inputListLHC12ePass1.txt", "./Plots/2012/12ePass1.pdf", "LHC12ePass1");
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2012/inputListLHC12fPass1.txt", "./Plots/2012/12fPass1.pdf", "LHC12fPass1");
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2012/inputListLHC12gPass1.txt", "./Plots/2012/12gPass1.pdf", "LHC12gPass1");
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2012/inputListLHC12iPass1.txt", "./Plots/2012/12iPass1.pdf", "LHC12iPass1");
  //
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2013/inputListLHC13bPass3.txt", "./Plots/2013/13bPass3.pdf", "LHC13bPass3");  
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2013/inputListLHC13cPass2.txt", "./Plots/2013/13cPass2.pdf", "LHC13cPass2");
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2013/inputListLHC13dPass2.txt", "./Plots/2013/13dPass2.pdf", "LHC13dPass2");
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2013/inputListLHC13ePass2.txt", "./Plots/2013/13ePass2.pdf", "LHC13ePass2");
  MakeQAperPeriod("./InputFilesFromGridPerPeriod/2013/inputListLHC13fPass2.txt", "./Plots/2013/13fPass2.pdf", "LHC13fPass2");
  
  //
  
}


//______________________________________________________________________________
void MakeQAperPeriod(const Char_t * inputList, const Char_t * outputFileName,const Char_t * periodLabel) {
  //
  // main function: read in list of input files and process them
  //
  ifstream in;
  in.open(inputList);
  // Read the input list of files 
  TString objfile;
  Int_t counter = 0;
  //
  // Initialize histograms
  //
  TList * list = InitializeHistograms(periodLabel); // one histogram per trigger type
  //
  while(in.good()) {
    in >> objfile;
    //
    // check if file is okay and protect
    //
    if (!objfile.Contains("root")) continue; 
    TFile inFile(objfile.Data());
    if (inFile.IsZombie() || !inFile.Get("fHistStatistics")) {
      Printf("*** WARNING: File %s corrupted! ***", objfile.Data());
      inFile.Close();
      continue;
    }
    //
    // get run number from file-name
    //
    Ssiz_t finish1 = objfile.Last('_');
    Ssiz_t start1  = objfile.Last('/');
    TString runNumber( objfile(start1 + 1, finish1 - 7 - start1) );
    //
    // fill relevant histograms
    //
    for(Int_t iTrig = 0; iTrig < 29; iTrig++) { // loop over trigger types
      TH1D * histTrig = (TH1D*) list->FindObject(Form("histAccepted_%s_&%i",periodLabel, iTrig));
      Float_t acceptedFraction = GetFraction((Char_t *) objfile.Data(), "Accepted", TMath::Nint( TMath::Power(2,iTrig)));
      Float_t acceptedFractionErr = GetFraction((Char_t *) objfile.Data(), "Accepted", TMath::Nint(TMath::Power(2,iTrig)), kTRUE);
      if (acceptedFraction > -0.5) {
	histTrig->Fill(runNumber.Data(), acceptedFraction);
	histTrig->SetBinError(histTrig->GetXaxis()->FindBin(runNumber.Data()), acceptedFractionErr);
      }
    }
    //
    counter++;
  }
  //
  // Draw histograms
  //  
  TCanvas * canvAcc = new TCanvas(Form("canvAcc_%s", periodLabel),"fraction of accepted events", 2000, 800);
  canvAcc->SetRightMargin(0.01);
  canvAcc->SetBottomMargin(0.15);
  //
  canvAcc->Print(Form("%s%s",outputFileName,"["));
  Int_t countFilledHists = 0;
  for(Int_t iTrig = 0; iTrig < 29; iTrig++) {
    TH1D * histTrig = (TH1D*) list->FindObject(Form("histAccepted_%s_&%i",periodLabel, iTrig));
    //
    canvAcc->cd();
    if (histTrig->GetEntries() > 0) {
      histTrig->LabelsDeflate();
      histTrig->LabelsOption("v","X");
      //
      histTrig->DrawCopy("Ep");
      AddFillSeparationLines(histTrig);
      //
      TString * triggerBitName = GetTriggerBitName(iTrig); 
      canvAcc->Print(outputFileName, Form("Title:%s", triggerBitName->Data()));
      delete triggerBitName;
      countFilledHists = iTrig;
    }
  }
  TString * triggerBitName = GetTriggerBitName(countFilledHists); 
  canvAcc->Print(Form("%s%s",outputFileName,"]"), Form("Title:%s", triggerBitName->Data()));
  delete triggerBitName;
  //
  //
  in.close();
}


//______________________________________________________________________________
Float_t GetFraction(const Char_t * inputFile, const Char_t * columnLabel, UInt_t triggerBit, Bool_t returnError) {
  //
  // fraction of offline trigger (or accepted) w.r.t. to all input triggers
  // -1 in case of file problems, -2 in case of histogram problems
  //
  TFile * inFile = TFile::Open(inputFile);
  if (!inFile) {
    return -1;
  }
  //
  TH2D * hist = (TH2D*) inFile->Get("fHistStatistics");
  if (!hist)   {
    inFile->Close();
    delete inFile;
    return -2;
  }
  //
  // locate relevant cell
  //
  Int_t binAllX = 1;
  Int_t binColumnLabelX  = -1;
  Int_t binTriggerBitY   = -1;
  //
  // find x-axis bin
  //
  for(Int_t iX =0; iX < hist->GetNbinsX() + 1; iX++) {
    TString label = hist->GetXaxis()->GetBinLabel(iX);
    if (label.CompareTo(columnLabel) == 0) binColumnLabelX = iX;
  }
  //
  // find y-axis bins
  //
  for(Int_t iY =0; iY < hist->GetNbinsY() + 1; iY++) {
    TString label = hist->GetYaxis()->GetBinLabel(iY);
    Ssiz_t start  = label.First('&');
    Ssiz_t finish = label.First('*');
    TString maskLabel( label( start+1, finish - start - 1 ) ); // returns the substring between the & and *
    UInt_t   maskInt = maskLabel.Atoi(); // Convert string to int (the int is in decimal form)
    if (maskInt == triggerBit && label.Contains("-B-")) binTriggerBitY = iY;
  }
  //
  Float_t fraction = 0.;
  Float_t all = hist->GetBinContent(binAllX, binTriggerBitY);
  if (all > 0 &&binColumnLabelX != -1 && binTriggerBitY != -1) {
    fraction = hist->GetBinContent(binColumnLabelX,binTriggerBitY)/all;
  } else {
    inFile->Close();
    delete inFile;
    return -1;
  }
  //
  inFile->Close();
  delete inFile;
  //
  if (returnError) return TMath::Sqrt(fraction*(1-fraction)/all);
  return fraction;

}


//______________________________________________________________________________
TList * InitializeHistograms(const Char_t * label) {
  //
  // Initialize all the relevant histograms and store them in a TList...
  // Cosmetic improvements to be put here...
  //
  TList * listOfHists = new TList();
  //
  //
  TH1F *histAccepted[29]; // one hisogram for each trigger type
  for(Int_t iTrig = 0; iTrig < 29; iTrig++) {
    TString * triggerBitName = GetTriggerBitName(iTrig);
    histAccepted[iTrig] = new TH1F(Form("histAccepted_%s_&%i",label,iTrig), 
				   Form("accepted events %s of trigger &%i: %s; ; accepted fraction", label, iTrig, triggerBitName->Data()),
				   1, 0, 1);
    delete triggerBitName;
    //
    listOfHists->Add(histAccepted[iTrig]);
    //
    histAccepted[iTrig]->SetBit(TH1::kCanRebin);
    histAccepted[iTrig]->GetYaxis()->SetRangeUser(-0.0000001,1.2);
    histAccepted[iTrig]->SetStats(0);
    histAccepted[iTrig]->SetMarkerStyle(20);
    histAccepted[iTrig]->SetMarkerColor(TColor::GetColor("#268bd2"));
    histAccepted[iTrig]->SetLineColor(TColor::GetColor("#268bd2"));
    histAccepted[iTrig]->SetMarkerSize(2.0);
  }
  //
  return listOfHists;

}


//______________________________________________________________________________
TString * GetTriggerBitName(Int_t bitNumber) {
  //
  // returns the name of a trigger associated to the bit represented by 2**bitNumber
  //
  //  TString def = "NOT YET IMPLEMENTED";
  //if (bitNumber > 29) return def;
  //
  TString names[29] = {"kMB",//           = BIT(0), // Minimum bias trigger, i.e. interaction trigger, offline SPD or V0 selection
		       "kINT7",//         = BIT(1), // V0AND trigger, offline V0 selection
		       "kMUON",//         = BIT(2), // Muon trigger, offline SPD or V0 selection
		       "kHighMult",//     = BIT(3), // High-multiplicity trigger (threshold defined online), offline SPD or V0 selection
		       "kEMC1",//         = BIT(4), // EMCAL trigger
		       "kCINT5",//        = BIT(5), // Minimum bias trigger without SPD. i.e. interaction trigger, offline V0 selection
		       "kCMUS5",//        = BIT(6), // Muon trigger, offline V0 selection
		       //kMUSPB        = BIT(6), // idem for PbPb
		       "kMUSH7",//        = BIT(7), // Muon trigger: high pt, single muon, offline V0 selection, CINT7 suite
		       //kMUSHPB//       = BIT(7), // idem for PbPb
		       "kMUL7",//         = BIT(8), // Muon trigger: like sign dimuon, offline V0 selection, CINT7 suite
		       //kMuonLikePB   = BIT(8), // idem for PbPb
		       "kMUU7",//         = BIT(9), // Muon trigger, unlike sign dimuon, offline V0 selection, CINT7 suite
		       //      kMuonUnlikePB = BIT(9), // idem for PbPb
		       "kEMC7",//         = BIT(10), // EMCAL trigger, CINT7 suite
		       //kEMC8         = BIT(10), // EMCAL trigger, CINT8 suite
		       "kMUS7",//         = BIT(11), // Muon trigger: low pt, single muon, offline V0 selection, CINT7 suite
		       "kPHI1",//         = BIT(12), // PHOS trigger, CINT1 suite
		       "kPHI",//          = BIT(13), // PHOS trigger, CINT7 suite
		       //kPHI8,//         = BIT(13), // PHOS trigger, CINT8 suite
		       //kPHOSPb       = BIT(13), // idem for PbPb
		       "kEMCEJE",//       = BIT(14), // EMCAL jet patch trigger
		       "kEMCEGA",//       = BIT(15), // EMCAL gamma trigger
		       "kCentral",//      = BIT(16), // PbPb central collision trigger
		       "kSemiCentral",//  = BIT(17), // PbPb semicentral collision trigger
		       "kDG5",//          = BIT(18), // Double gap diffractive
		       "kZED",//          = BIT(19), // ZDC electromagnetic dissociation
		       //kSPI7,//         = BIT(20), // Power interaction trigger
		       "kSPI",//          = BIT(20), // Power interaction trigger
		       "kINT8",//                 = BIT(21), // CINT8 trigger: 0TVX (T0 vertex) triger
		       "kMuonSingleLowPt8",//     = BIT(22), // Muon trigger : single muon, low pt, T0 selection, CINT8 suite
		       "kMuonSingleHighPt8",//    = BIT(23), // Muon trigger : single muon, high pt, T0 selection, CINT8 suite
		       "kMuonLikeLowPt8",//       = BIT(24), // Muon trigger : like sign muon, low pt, T0 selection, CINT8 suite
		       "kMuonUnlikeLowPt8",//     = BIT(25), // Muon trigger : unlike sign muon, low pt, T0 selection, CINT8 suite
		       "kMuonUnlikeLowPt0",//     = BIT(26), // Muon trigger : unlike sign muon, low pt, no additional L0 requirement
		       "kUserDefined",//  = BIT(27), // Set when custom trigger classes are set in AliPhysicsSelection, offline SPD or V0 selection
		       "kTRD"//          = BIT(28), // TRD trigger
  };
  //
  return new TString(names[bitNumber]);


}


//______________________________________________________________________________
UInt_t GetFillFromRunNumber(UInt_t runNumber) {
  //
  // get the fill from the run number in an ugly way based on a .txt file
  // which is extracted from the logbook:
  // https://alice-logbook.cern.ch/logbook/date_online.php?p_cont=sb&p_rsob=l.run&p_rsob_dir=DESC&prsf_rwb=Yes&p_rspn=0&p_output=TXT&p_oh=yes&p_od=2&p_os=1
  //
  UInt_t fillNumber = 0;
  ifstream in;
  in.open("RunFill.list");
  if (!in) {
    //Printf("Cannot display fill, because file RunFill.list is not available in current directory.");
    return 0;
  }
  //
  // Read the input list of files 
  //
  TString objfile;
  //
  while(in.good()) {
    in >> objfile;
    if (!objfile.Contains(";")) continue; // protection
    Ssiz_t delim = objfile.First(';');
    Ssiz_t end   = objfile.Length();
    TString runString(objfile(0, delim  ) );
    Int_t run = runString.Atoi();
    if ((UInt_t) TMath::Abs(run) == runNumber) {
      TString fillString(objfile(delim+1, end - delim  ) );
      return TMath::Abs(fillString.Atoi());
    }
  }
  //
  //
  return fillNumber;

}


//______________________________________________________________________________
void  AddFillSeparationLines(TH1D * hist) {
  //
  // add the fill separation lines
  //
  for(Int_t iBin = 1; iBin < hist->GetXaxis()->GetNbins(); iBin++) {
    UInt_t runNumberOld = atoi(hist->GetXaxis()->GetBinLabel(iBin));    
    UInt_t runNumberNew = atoi(hist->GetXaxis()->GetBinLabel(iBin + 1));    
    if (GetFillFromRunNumber(runNumberOld) != GetFillFromRunNumber(runNumberNew)) {
      TLine * fillSeparationLine = new TLine(iBin,0,iBin,1.2);
      fillSeparationLine->SetLineColor(kRed);
      fillSeparationLine->SetLineWidth(1);
      fillSeparationLine->Draw();
    }
  }

}
