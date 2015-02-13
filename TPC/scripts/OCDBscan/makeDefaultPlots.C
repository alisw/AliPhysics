/// \file makeDefaultPlots.C
///
/// ~~~{.cpp}
/// .L $ALICE_ROOT/TPC/scripts/OCDBscan/makeDefaultPlots.C
/// Init();
/// MakeAliases();
/// ~~~

TTree * guiTree=0;
TObjArray * picArray = new TObjArray;
TCut cutUser="run>114500";

void makeDefaultPlots(){
  /// make default trend plots

  Init();
  MakeAliases();
  MakePlotDrift();
  //
  TFile f("calibTrend.root","recreate");
  picArray->Write();
  f.Close();
}

void Init(){
  /// Init trees and style
  /// summary tree supposed to be in the given file

  TFile *f = new TFile("calibTimeSummary.root");
  guiTree=(TTree*)f->Get("dcs");
  guiTree->SetMarkerStyle(25);
  guiTree->SetMarkerSize(0.5);
}

void MakeAliases(){
  /// Make "reasonable" alaiase for variables
  /// Define te cut variables for visualization
  ///
  /// Drift velocity

  guiTree->SetAlias("vitsOK","abs(dits)<3600*0.5");
  //TOF alaieases
  guiTree->SetAlias("vtofOK","abs(ALIGN_TOF_TPC_DELTATHETA-ALIGN_TOF_TPC_DELTATHETA)<0.0003&&abs(ALIGN_TOF_TPC_DELTAPSI-ALIGN_TOF_TPC_DELTAPSI)<0.0003&&abs(ALIGN_TOF_TPC_DELTATHETA)<0.003&&abs(ALIGN_TOF_TPC_DELTAPSI)<0.003");
  guiTree->SetAlias("vdriftTOF","ALIGN_TOFB_TPC_DRIFTVD");
  //  angular alignent within limits
  //
  //
  // distance to measurement 1/2 hour
  // Goofie
  //
  guiTree->SetAlias("vgoofieOK","abs(goofie.fElements[3]-2.7)<0.15&&abs(goofie.fElements[5]/goofie.fElements[6]-0.8)<0.05&&abs(goofie.fElements[2]-9.)<4&&abs(goofie.fElements[9]/goofieMedian.fElements[9]-1)<0.05&&abs(goofie.fElements[3]-goofieMean.fElements[3])<0.02");
  // Cuts:
  // reasonable drift velocity - outlyer removal
  // ration of the peak are near/far
  // reasonable CO2 information
  // Q value/median in range
  // V drift value-mean in range
  //
  // Goofie normalization:
  guiTree->Draw("goofie.fElements[14]","vgoofieOK");  // nominal pressure
  Double_t goofiePNom = guiTree->GetHistogram()->GetMean();
  guiTree->Draw("goofie.fElements[15]","vgoofieOK");  // nominal pressure
  Double_t goofieT1PNom = guiTree->GetHistogram()->GetMean();
  guiTree->Draw("goofie.fElements[3]","vgoofieOK");
  Double_t goofieDriftNom = guiTree->GetHistogram()->GetMean(); // nominal drift velocity
  guiTree->Draw("goofie.fElements[9]","vgoofieOK");
  Double_t goofieAreaFar = guiTree->GetHistogram()->GetMean();
  guiTree->Draw("goofie.fElements[10]","vgoofieOK");
  Double_t goofieAreaNear = guiTree->GetHistogram()->GetMean();
  //
  guiTree->SetAlias("goofieP",Form("goofie.fElements[14]/%f",goofiePNom));
  guiTree->SetAlias("goofieT1P",Form("goofie.fElements[15]/%f",goofieT1PNom));
  guiTree->SetAlias("goofieVd",Form("goofie.fElements[3]/%f",goofieDriftNom));
  guiTree->SetAlias("goofieAreaNear",Form("(goofie.fElements[9]/%f)",goofieAreaFar));
  guiTree->SetAlias("goofieAreaFar",Form("(goofie.fElements[10]/%f)",goofieAreaNear));
  guiTree->Draw("(goofieVd*goofieP-1)-vdriftITS","vitsOK&&vgoofieOK");
  Double_t goofieITSNorm = guiTree->GetHistogram()->GetMean();
  guiTree->SetAlias("vdriftGoofie",Form("(goofieVd*goofieP-1)-%f",goofieITSNorm));
  //

}

void MakePlotDrift(){
  /// Compare drift velocity calibration
  /// with  Goofie calibration

  Int_t entries =0;
  TGraph *graphITS =0;
  TGraph *graphTOF =0;
  TGraph *graphGoofie =0;
  entries = guiTree->Draw("100*vdriftITS:time","vitsOK"+cutUser,"");
  graphITS = new TGraph(entries,guiTree->GetV2(), guiTree->GetV1());
  entries = guiTree->Draw("100*vdriftTOF:time","vitsOK&&vtofOK"+cutUser,"");
  graphTOF = new TGraph(entries,guiTree->GetV2(), guiTree->GetV1());
  guiTree->Draw("100*vdriftGoofie:time","vitsOK&&vgoofieOK"+cutUser,"same");
  graphGoofie = new TGraph(entries,guiTree->GetV2(), guiTree->GetV1());
  //
  graphITS->SetMarkerStyle(25);
  graphITS->SetMarkerSize(0.5);
  graphITS->SetMarkerColor(2);
  graphTOF->SetMarkerStyle(26);
  graphTOF->SetMarkerSize(0.5);
  graphTOF->SetMarkerColor(3);
  graphGoofie->SetMarkerStyle(26);
  graphGoofie->SetMarkerSize(0.5);
  graphGoofie->SetMarkerColor(4);
  //
  TCanvas * cDrift = new TCanvas("DriftTrend","Driftt Trend",900,600);
  cDrift->cd();
  graphITS->Draw("ap");
  graphITS->GetXaxis()->SetTimeDisplay(1);
  graphITS->GetYaxis()->SetTitle("v_{dcorr} (%)");
  graphGoofie->Draw("p");
  graphTOF->Draw("p");
  TLegend *legend = new TLegend(0.7,0.7,1,1, "Drift velocity correction");
  legend->AddEntry(graphITS,"TPC-ITS");
  legend->AddEntry(graphTOF,"TPC-TOF");
  legend->AddEntry(graphGoofie,"Goofie");
  legend->Draw();
  // Add picture to the array
  //
  picArray->AddLast(cDrift->Clone());
  //
  //

}
