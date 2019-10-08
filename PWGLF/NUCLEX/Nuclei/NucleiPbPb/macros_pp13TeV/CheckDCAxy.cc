#include "src/Common.h"
#include "src/Plotting.h"
#include "src/Utils.h"
using namespace utils;

#include <TFile.h>
#include <TCanvas.h>
#include <TH3F.h>

void CheckDCAxy(bool isTPC = false){
	TFile mc_file(kMCfilename.data());
  TFile data_file(kDataFilename.data());
  std::string output_string;
	if(isTPC){
		output_string = Form("%sDCAxy_check.root",kBaseOutputDir.data());
	} else {
		output_string = Form("%sDCAxy_TPC_check.root",kBaseOutputDir.data());
	}
  TFile output_file(output_string.data(),"recreate");
	TObjArray obj(2);
  for (auto list_key : *data_file.GetListOfKeys()) {
    if (string(list_key->GetName()).find(kFilterListNames.data()) == string::npos) continue;
    if (string(list_key->GetName()).find("_MV") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid2") != string::npos) continue;
    if (string(list_key->GetName()).find("_pid3") != string::npos) continue;
    if (string(list_key->GetName()).find("_p_selection") != string::npos) continue;
		//
		TTList* mcList = (TTList*)mc_file.Get(list_key->GetName());
    TTList* dtList = (TTList*)data_file.Get(list_key->GetName());

		TH3F* primaries = nullptr;
    TH3F* antideuterons = nullptr;

		if(isTPC){
      primaries = (TH3F*)mcList->FindObject("fMDCAPrimaryTPC");
      Requires(primaries,"Missing primaries");
      antideuterons = (TH3F*)dtList->FindObject("fADCAxyTPC");
      Requires(antideuterons,"Missing data");
    } else {
      primaries = (TH3F*)mcList->FindObject("fMDCAPrimaryTOF");
      Requires(primaries,"Missing primaries");
      antideuterons = (TH3F*)dtList->FindObject("fADCAxyTOF");
      Requires(antideuterons,"Missing data");
    }
		TDirectory *base_dir = output_file.mkdir(list_key->GetName());
		TDirectory *primdir = base_dir->mkdir("Primaries");
		TDirectory *antidir = base_dir->mkdir("Antideuterons");
		TDirectory *compdir = base_dir->mkdir("Comparison");
		for(int iC=0; iC<kCentLength; iC++){
			antidir->mkdir(Form("%d",iC));
			compdir->mkdir(Form("%d",iC));
		}
		for (int iB = 1; iB <= kNPtBins; ++iB) {
      float bin_center = (kPtBins[iB]+kPtBins[iB-1])/2;
      if (((bin_center < kPtRangeMatCorrectionTPC[0] || bin_center > kPtRangeMatCorrectionTPC[1]) && isTPC) || 
          ((bin_center < kPtRangeMatCorrectionTOF[0] || bin_center > kPtRangeMatCorrectionTOF[1]) && !isTPC)) continue;
			int iBin = antideuterons->GetYaxis()->FindBin(bin_center);
      int iBin1 = primaries->GetYaxis()->FindBin(antideuterons->GetYaxis()->GetBinLowEdge(iBin)+0.005);
      int iBin2 = primaries->GetYaxis()->FindBin(antideuterons->GetYaxis()->GetBinLowEdge(iBin+1)-0.005);
			//TH1D* pr_tmp = (TH1D*)primaries->ProjectionZ(Form("pr_tmp_%i",iB),kCentBinsArray[kCentLength-1][0],kCentBinsArray[kCentLength-1][1],iBin,iBin);
      TH1D* pr = (TH1D*)primaries->ProjectionZ(Form("pr_tmp_%i",iB),kCentBinsArray[kCentLength-1][0],kCentBinsArray[kCentLength-1][1],iBin1,iBin2);//(TH1D*)pr_tmp->Rebin(kNDCAbins,Form("pr_%i",iB),kDCAbins);
      pr->SetTitle(Form("%4.1f < p_{T} #leq %4.1f (GeV/#it{c});counts; DCA_{xy} (cm)",kPtBins[iB-1],kPtBins[iB]));
			//pr->Scale(1.,"width");
			primdir->cd();
			pr->Write();
			double pr_integral = pr->Integral(pr->FindBin(-0.12),pr->FindBin(0.12));
			TH1D* anti[kCentLength] = {nullptr};
      for(int iC=kCentLength; iC--;){
        //TH1D* anti_tmp = (TH1D*)antideuterons->ProjectionZ(Form("anti_tmp_%i",iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iBin,iBin,"e");
        anti[iC] = (TH1D*)antideuterons->ProjectionZ(Form("anti_tmp_%i",iB),kCentBinsArray[iC][0],kCentBinsArray[iC][1],iBin,iBin,"e");//(TH1D*) anti_tmp->Rebin(kNDCAbins,Form("dt_%i_C%i",iB,iC),kDCAbins);
        anti[iC]->SetTitle(Form("Mult: %1.0f - %1.0f %% , %4.1f < p_{T} #leq %4.1f (GeV/#it{c})",kCentLabels[iC][0],kCentLabels[iC][1],kPtBins[iB-1],kPtBins[iB]));
				plotting::SetHistStyle(anti[iC],kRed,20);
				//anti[iC]->Scale(1.,"width");
				antidir->cd(Form("%d",iC));
				anti[iC]->Write();
				compdir->cd(Form("%d",iC));
				TCanvas* cv = new TCanvas(Form("comp_%d_%d",iC,iB),Form("comp_%d_%d",iC,iB));
				double anti_integral = anti[iC]->Integral(anti[iC]->FindBin(-0.12),anti[iC]->FindBin(0.12));
				cv->cd();
				double scale_factor = anti_integral/pr_integral;
				TH1D* pr_scaled = (TH1D*)pr->Clone("pr_scaled");
				pr_scaled->Scale(scale_factor);
				anti[iC]->GetXaxis()->SetRangeUser(-0.5,0.5);
				anti[iC]->Draw();
				pr_scaled->Draw("histo same");
				cv->SetLogy();
				cv->Write();
				delete cv;
      }
		}
	}
}