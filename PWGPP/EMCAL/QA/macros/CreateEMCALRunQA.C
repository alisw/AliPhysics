#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "AliEMCALGeometry.h"
#include <TColor.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TTree.h>
#include <TPRegexp.h>
#include <TList.h>
#include <TObjString.h>
#include <TDatime.h>
#include <TError.h>
#include <AliLog.h>
#endif

///
/// \file CreateEMCALRunQA.C
/// \brief QA at run level
///
/// This macro produces runLevelQA for EMCal from a QAresults.root file.
/// Re-factored for automatic QA processing A. Shabetai.
///
/// \author Yaxian Mao, <Yaxian.Mao@cern.ch>, CCNU
/// \author Alexis Mas, <Alexis.Mas@subatech.in2p3.fr>, SUBATECH
/// \author Alexandre Shabetai, <Alexandre.Shabetai@subatech.in2p3.fr>  SUBATECH
/// \author Marie Germain, <Marie.Germain@subatech.in2p3.fr>, SUBATECH
/// \author Erwann Masson, <Erwann.Masson@subatech.in2p3.fr>, SUBATECH
///

Int_t DrawOccupancy(Long_t run, TString period, TString pass, TString fTrigger, TString system, TFile* f, TFile* fout, AliEMCALGeometry* geom, Int_t SavePlots);
Int_t DrawRun(Long_t run, TString period, TString pass, TString fTrigger, TFile *f, TFile* fout, Int_t SavePlots, Int_t nSM, Bool_t kFilter);
Int_t TrendingEMCALTree(Long_t RunId, TString fCalorimeter, TString system, TString period, TString pass, const int n, TList* TriggersList, TFile* f, TFile *fout, Int_t SavePlots);

TH2F* FormatRunHisto(TH2F* aHisto, const char* title, const char* YTitle="");
TH2F* HistoPerMod(TH2F* name, const char* title);
TH2*  AutoZoom(TH2* H, Option_t* aType="all", Int_t EntryMin=0);
int   FindNumberOfSM(TFile* f, TString fTrigger, TString period);

TString QAPATH;
TString QAPATHF = "./";

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void set_plot_style(){

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Double_t pi0massP2(Double_t *x, Double_t *par){

  Double_t gaus;

  if(par[2] != 0.)
    gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) / (2*par[2]*par[2]) );
  else
    gaus = 99999999.;

  Double_t back = par[3] + par[4]*x[0] + par[5]*x[0]*x[0];

  return gaus+back;

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Double_t pi0massP1(Double_t *x, Double_t *par){

  Double_t gaus = par[0] * TMath::Exp( -(x[0]-par[1])*(x[0]-par[1]) / (2*par[2]*par[2]) );
  Double_t back = par[3] + par[4]*x[0];

  return gaus+back;

}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Double_t fitE(Double_t *x, Double_t *par){

  Double_t levy;

  levy = par[0] * TMath::Exp( -par[1]/x[0]) * TMath::Power(x[0], -par[2]) ;

  return levy;

}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
int CreateEMCALRunQA(const char* filename, TString RunId, TString period, TString pass, Int_t SavePlots = 0, Bool_t filter=0, TString fTrigger = "", TString system = "", TString fCalorimeter = "EMCAL"){

  QAPATH = TString(gSystem->Getenv("QAPATH"));
  if(QAPATH.IsNull())
    QAPATH = QAPATHF;

  if(! QAPATH.BeginsWith("./"))
    QAPATH = QAPATH + RunId + "/";

  AliLog::SetGlobalLogLevel(AliLog::kError);
  TFile *f = new TFile(filename);
  AliLog::SetGlobalLogLevel(AliLog::kInfo);

  if(f->IsZombie()){
    Error(__FUNCTION__, Form("Error openning the input file %s", filename));
    return -1;
  }

  TList* TriggersList = new TList();

  if(fTrigger==""){
    TPMERegexp name_re("CaloQA_\\w+");
    TObjLink* link = f->GetListOfKeys()->FirstLink();

    while(link){
      TString name = link->GetObject()->GetName();
      if(name_re.Match(name)){
	TriggersList->Add(link->GetObject());
	if(TString(filename).Contains("barrel") && ! name.Contains("default"))  TriggersList->Remove(link->GetObject());
	if(TString(filename).Contains("outer")  && ! name.Contains("EMC"))      TriggersList->Remove(link->GetObject());
      }
      link = link->Next();
    }
  }
  else{
    TriggersList->Add(new TObjString(fTrigger.Data()));
  }

  if(!TriggersList->GetEntries()){
    Error(__FUNCTION__, "No trigger found!");
    return -2;
  }

  int nSM = FindNumberOfSM(f, ((TObjString*)TriggersList->Last())->GetString(), period);

  if(nSM<0){
    Error(__FUNCTION__, "Could not find the number of super modules!");
    return -3;
  }
  Info(__FUNCTION__, Form("%i super modules were discovered", nSM));

  TString GeomName;
  if(nSM <= 4)      { nSM = 4;  GeomName = "EMCAL_FIRSTYEARv1";            }
  else if(nSM <= 10){ nSM = 10; GeomName = "EMCAL_COMPLETEv1";             }
  else if(nSM <= 12){ nSM = 12; GeomName = "EMCAL_COMPLETE12SMv1";         }
  else              { nSM = 20; GeomName = "EMCAL_COMPLETE12SMv1_DCAL_8SM";}

  AliEMCALGeometry *geom = new AliEMCALGeometry(GeomName.Data(), "EMCAL");
  Info(__FUNCTION__, Form("Using %i super modules and the Geometry %s", nSM, GeomName.Data()));

  TFile *fout = new TFile(TString( QAPATH + period+"_"+pass + fTrigger+"_"+ (Long_t)RunId.Atoi() +"_QAplots.root").Data(), "RECREATE");

  if((system.IsNull()) && (period.EndsWith("h")))
    system = "PbPb";

  Int_t ret=0;
  TIter next(TriggersList);
  while(TObject *obj = next()){
    fTrigger= TString(obj->GetName());
    ret -= DrawOccupancy(RunId.Atoi(), period, pass, fTrigger, system, f, fout, geom, SavePlots);
    ret -= DrawRun(RunId.Atoi(), period, pass, fTrigger, f, fout, SavePlots, nSM, filter);
  }
  ret-= TrendingEMCALTree(RunId.Atoi(), fCalorimeter, system, period, pass, nSM, TriggersList, f, fout, SavePlots);

  f->Close();
  fout->Close();
  delete f;
  delete geom;

  return ret;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------
Int_t DrawOccupancy(Long_t  run, TString period, TString pass, TString fTrigger, TString system, TFile* f, TFile* fout, AliEMCALGeometry* geom, Int_t SavePlots){

  set_plot_style();
  gStyle->SetOptStat(0);

  TH1::AddDirectory(kFALSE);
  TH2D *hEnergyMapReal = new TH2D("hEnergyMapReal", "", 96, -48, 48, 120, -0, 120);
  TH2D *hOccupancyMapReal = new TH2D("hOccupancyMapReal", "", 96, -48, 48, 120, -0, 120);

  hEnergyMapReal->SetXTitle("#eta (bin)");
  hEnergyMapReal->SetYTitle("#varphi (bin)");
  hEnergyMapReal->SetZTitle("E (GeV)/event");
  hEnergyMapReal->GetYaxis()->SetTitleOffset(1.2);
  hEnergyMapReal->GetZaxis()->SetLabelSize(0.02);
  hEnergyMapReal->GetZaxis()->SetTitleOffset(1.36);

  hOccupancyMapReal->SetXTitle("#eta (bin)");
  hOccupancyMapReal->SetYTitle("#varphi (bin)");
  hOccupancyMapReal->GetYaxis()->SetTitleOffset(1.2);
  hOccupancyMapReal->GetZaxis()->SetLabelSize(0.02);

  Int_t nSupMod, nModule, nIphi, nIeta;
  Int_t iphi, ieta;
  Int_t realbineta=0;
  Int_t realbinphi=0;

  // NO MASK
  Int_t mask[1] = {2222222};

  TH2F *hCellAmplitude;
  TH1F *hNEvents;
  Int_t Events;
  Int_t n=0;

  TString direct;
  if(!fTrigger.Contains("QA")){
    direct = "CaloQA_";
  }
  direct += fTrigger;

  Bool_t dirok = f->cd(direct);
  if(!dirok){
    Error(__FUNCTION__, Form("No input drectory %s", direct.Data()));
    return -1;
  }

  TList *outputList = (TList*)gDirectory->Get(direct);
  if(!outputList){
    Error(__FUNCTION__, "No input list! ");
    return -1;
  }
  outputList->SetOwner();

  fout->mkdir(Form("%s/%s/%ld/%s/%s", period.Data(), pass.Data(), run, "RunLevelQA", fTrigger.Data()));
  fout->cd();
  fout->Cd(Form("%s/%s/%ld/%s/%s", period.Data(), pass.Data(), run, "RunLevelQA", fTrigger.Data()));

  hNEvents = (TH1F *)outputList->FindObject("hNEvents");
  if(!hNEvents){
    Error(__FUNCTION__, Form("hNEvent histogram not found for trigger %s ", fTrigger.Data()));
    return -2;
  }

  Events = (Int_t)hNEvents->GetEntries();
  if(Events==0){
    Error(__FUNCTION__, Form("No event in trigger %s", fTrigger.Data()));
    return -3;
  }

  TH1F* hE = (TH1F *)outputList->FindObject("EMCAL_hE");
  if(!hE || (!hE->GetEntries())){
    Error(__FUNCTION__, Form("hE Histogram not found or it is empty for trigger %s, EMCal was not in this run?", fTrigger.Data()));
    return -4;
  }

  Double_t Eth=1;
  if(system=="PbPb"){
    Eth = 5.;
    if(fTrigger.Contains("EMC"))
      Eth=20.;
  }
  if(system=="pp"){
    Eth = 1.;
    if(fTrigger.Contains("EMC"))
      Eth=5.;
  }

  hCellAmplitude = (TH2F *)outputList->FindObject("EMCAL_hAmpId");

  for(Int_t i = 0; i < geom->GetNCells() ; i++){
    Double_t Esum = 0;
    Double_t Nsum = 0;

    for(Int_t j = 1; j <= hCellAmplitude->GetNbinsX(); j++){
      Double_t E = hCellAmplitude->GetXaxis()->GetBinCenter(j);
      Double_t N = hCellAmplitude->GetBinContent(j, i+1);

      if(E < 0.3)
	continue;

      if(E <= Eth){
	Esum += E*N;
	Nsum += N;
      }
    }

    Int_t absId = i;
    if(n!=0){
      if(mask[n]<=mask[n-1])
	Warning(__FUNCTION__, "The list of bad cells is not sorted !!");
    }

    if(i==mask[n]){
      n++ ; 
      continue;
    } // skip bad cells

    geom->GetCellIndex(absId, nSupMod, nModule, nIphi, nIeta);
    geom->GetCellPhiEtaIndexInSModule(nSupMod, nModule, nIphi, nIeta, iphi, ieta);

    realbinphi = 120-(nSupMod/2)*24 -iphi -1;
    if(nSupMod%2==0)
      realbineta= 48-ieta -1;
    if(nSupMod%2==1)
      realbineta= -ieta -1;

    hEnergyMapReal->Fill(realbineta, realbinphi, Esum/(Double_t)Events);
    // hOccupancyMapReal->Fill(realbineta, realbinphi, Nsum/(Double_t)Events);
  }

  cout <<" Run: " << run << " trigger: " << fTrigger << " N_events: "<<Events<<endl;

  TPMERegexp r("_\\w+");
  TString Energy;   Energy   = QAPATH + "MapEnergy"  + fTrigger(r) + ".pdf";
  TString Energy2;  Energy2  = QAPATH + "MapEnergy"  + fTrigger(r) + ".png";
  TString Entries;  Entries  = QAPATH + "MapEntries" + fTrigger(r) + ".pdf";
  TString Entries2; Entries2 = QAPATH + "MapEntries" + fTrigger(r) + ".png";

  TCanvas *c1 = new TCanvas("Energymap", "Energy Map", 600, 600);
  c1->SetFillColor(0);
  c1->SetGrid();
  c1->SetRightMargin(0.14);

  TString title = "run ";
  title += run ;
  if(fTrigger.Contains("EMC"))
    title += " EMC ";
  else
    title += " MB ";
  title += " Summed energy map";

  hEnergyMapReal->SetTitle(title);
  AutoZoom(hEnergyMapReal, "miny")->DrawCopy("colz");

  if(nSupMod<=12){
    if(SavePlots==2)
      c1->SaveAs(Energy);
    if(SavePlots)
      c1->SaveAs(Energy2);
    c1->Write();
  }
  delete c1;

  // TCanvas *c2 = new TCanvas("Occupancy", "Occupancy Map", 600, 600);
  // c2->SetFillColor(0);
  // c2->SetGrid();
  // c2->SetRightMargin(0.14);
  // TString title2 = "run ";
  // title2 += run ;
  // if(fTrigger.Contains("EMC")){ title2 += " EMC ";} else{ title2 += " MB ";}
  // title2 += " Occupancy map";
  // hOccupancyMapReal->SetTitle(title2);
  // AutoZoom(hOccupancyMapReal, "miny")->DrawCopy("colz");
  // if(SavePlots==2) c2->SaveAs(Entries);
  // if(SavePlots) c2->SaveAs(Entries2);
  // c2->Write();
  // delete c2;

  if(outputList)
    outputList->Delete();
  delete hEnergyMapReal;
  delete hOccupancyMapReal;

  return 0;
}

//-----------------------------------------------------------------------------------------------------------------------
Int_t DrawRun(const Long_t  run, TString period, TString pass, TString fTrigger, TFile *f, TFile *fout, Int_t SavePlots, Int_t nSM, Bool_t kFilter){

  TString direct;
  if(!fTrigger.Contains("QA")){
    direct = "CaloQA_";
  }
  direct += fTrigger;

  f->cd(direct);
  if(!direct){
    Error(__FUNCTION__, Form("No input directory %s", direct.Data()));
    return -1;
  }

  TList *outputList = (TList*)gDirectory->Get(direct);
  if(!outputList){
    Error(__FUNCTION__, Form("No input list! %s", direct.Data()));
    return -2;
  }
  outputList->SetOwner();

  if(kFilter){
    fout->mkdir(Form("%s/%s/%ld/%s", period.Data(), pass.Data(), run, fTrigger.Data()));
    fout->cd();
    fout->Cd(Form("%s/%s/%ld/%s", period.Data(), pass.Data(), run, fTrigger.Data()));
    outputList->Write();
  }
  fout->cd();
  fout->Cd(Form("%s/%s/%ld/%s/%s", period.Data(), pass.Data(), run, "RunLevelQA", fTrigger.Data()));

  set_plot_style();
  gStyle->SetPalette(1);
  set_plot_style();
  // gStyle->SetOptStat(1);

  TH1::AddDirectory(kFALSE);
  TString outfilename;
  TString outfilename2;
  const char* legend="";
  TPMERegexp r("_\\w+");

  if(fTrigger.Contains("EMC"))
    legend = Form(" Run %d EMC ", (int)run);
  else
    legend = Form(" Run %d MB ", (int)run);

  TH1F* hNEvents =(TH1F *)outputList->FindObject("hNEvents");
  if(!hNEvents){
    Error(__FUNCTION__, Form("hNEvent histogram not found for trigger %s ", fTrigger.Data()));
    return -3;
  }

  Int_t Events = (Int_t)hNEvents->GetEntries();
  if(Events==0){
    Error(__FUNCTION__, Form("No event in trigger %s", fTrigger.Data()));
    return -4 ;
  }

  TH1F* hE =(TH1F *)outputList->FindObject("EMCAL_hE");
  if(!hE || (!hE->GetEntries())){
    Error(__FUNCTION__, Form("hE Histogram not found or it is empty for trigger %s, EMCal was not in this run?", fTrigger.Data()));
    return -4;
  }

  TCanvas* c00 = new TCanvas("Occupancy map", "Occupancy", 600, 600);
  c00->SetLogz();
  c00->SetFillColor(0);
  c00->SetBorderSize(0);
  c00->SetFrameBorderMode(0);

  TH2F*  hOccupancyMapReal =(TH2F *)outputList->FindObject("EMCAL_hEtaPhi");
  if(!hOccupancyMapReal){
    Error(__FUNCTION__, Form("EMCAL_hEtaPhi: Histogram for trigger %s not found!", fTrigger.Data()));
    return -5;
  }
  FormatRunHisto(hOccupancyMapReal, Form("Occupancy%s", legend), "#varphi (bin)");
  hOccupancyMapReal->SetXTitle("#eta (bin)");

  AutoZoom(hOccupancyMapReal, "all")->DrawCopy("colz");
  outfilename =  QAPATH + "OccupancyMap" + fTrigger(r) + ".pdf" ;
  outfilename2 = QAPATH + "OccupancyMap" + fTrigger(r) + ".png" ;

  if(SavePlots==2)
    c00->SaveAs(outfilename);
  if(SavePlots)
    c00->SaveAs(outfilename2);
  c00->Write();
  delete c00;

  TCanvas* c0 = new TCanvas("Grid Cell", "Occupancy from Grid Cells", 600, 600);
  c0->SetLogz();
  c0->SetFillColor(0);
  c0->SetBorderSize(0);
  c0->SetFrameBorderMode(0);

  TH2F*  hGridCells =(TH2F *)outputList->FindObject("EMCAL_hGridCells");
  if(!hGridCells){
    Error(__FUNCTION__, Form("EMCAL_hGridCells: Histogram for trigger %s not found!", fTrigger.Data()));
    return -5;
  }
  FormatRunHisto(hGridCells, Form("Occupancy from Grid Cells%s", legend), "#varphi row columns");
  hGridCells->SetXTitle("#eta row columns");
  AutoZoom(hGridCells, "all")->DrawCopy("colz");

  outfilename =  QAPATH + "GridCells" + fTrigger(r) + ".pdf" ;
  outfilename2 = QAPATH + "GridCells" + fTrigger(r) + ".png" ;

  if(SavePlots==2)
    c0->SaveAs(outfilename);
  if(SavePlots)
    c0->SaveAs(outfilename2);
  c0->Write();
  delete c0;

  TCanvas* c1 = new TCanvas("TimeVsE", "Cluster Time vs. Energy", 600, 600);
  c1->SetLogz();
  c1->SetFillColor(0);
  c1->SetBorderSize(0);
  c1->SetFrameBorderMode(0);

  TH2F*  hClusterTimeEnergy =(TH2F *)outputList->FindObject("EMCAL_hClusterTimeEnergy");
  if(!hClusterTimeEnergy){
    Error(__FUNCTION__, Form("EMCAL_hClusterTimeEnergy: Histogram for trigger %s not found!", fTrigger.Data()));
    return -5;
  }
  FormatRunHisto(hClusterTimeEnergy, Form("Time vs. Energy%s", legend), "EMCal ToF (ns)");
  AutoZoom(hClusterTimeEnergy, "maxx")->DrawCopy("colz");

  outfilename =  QAPATH + "TimeRun" + fTrigger(r) + ".pdf" ;
  outfilename2 = QAPATH + "TimeRun" + fTrigger(r) + ".png" ;

  if(SavePlots==2)
    c1->SaveAs(outfilename);
  if(SavePlots)
    c1->SaveAs(outfilename2);
  c1->Write();
  delete c1;

  if(pass!="simu"){
    TCanvas* c1a = new TCanvas("TimeVsAbsId", "Cell Id vs. time ", 600, 600);
    c1a->SetLogz();
    c1a->SetFillColor(0);
    c1a->SetBorderSize(0);
    c1a->SetFrameBorderMode(0);

    TH2F*  hTimeAbsId =(TH2F *)outputList->FindObject("EMCAL_hTimeId");
    if(!hTimeAbsId){
      Error(__FUNCTION__, Form("EMCAL_hTimeId: Histogram for trigger %s not found!", fTrigger.Data()));
      return -5;
    }
    FormatRunHisto(hTimeAbsId, Form("Cell ID vs. Time%s", legend), "Cell AbsID");
    hTimeAbsId->RebinY(16);
    hTimeAbsId->SetMinimum(32);
    AutoZoom(hTimeAbsId, "maxx")->DrawCopy("colz");

    outfilename =  QAPATH + "TimeAbsIdRun" + fTrigger(r) + ".pdf" ;
    outfilename2 = QAPATH + "TimeAbsRun" + fTrigger(r) + ".png" ;

    if(SavePlots==2)
      c1a->SaveAs(outfilename);
    if(SavePlots)
      c1a->SaveAs(outfilename2);
    c1a->Write();
    delete c1a;
  }

  // TCanvas* c1aLG = new TCanvas("TimeVsAbsIdLG", "Cell Id vs. time LG ", 600, 600);
  // c1aLG->SetLogz();
  // c1aLG->SetFillColor(0);
  // c1aLG->SetBorderSize(0);
  // c1aLG->SetFrameBorderMode(0);

  // TH2F*  hTimeAbsIdLG =(TH2F *)outputList->FindObject("EMCAL_hTimeIdLG");
  // if(!hTimeAbsIdLG){
  //   Error(__FUNCTION__, Form("EMCAL_hTimeIdLG: Histogram for trigger %s not found!", fTrigger.Data()));
  //   return -5;
  // }
  // FormatRunHisto(hTimeAbsIdLG, Form("Cell ID vs. Time LG%s", legend), "cell AbsID");

  // AutoZoom(hTimeAbsIdLG, "maxx")->DrawCopy("colz");
  // outfilename =  QAPATH + "TimeAbsIdRunLG" + fTrigger(r) + ".pdf" ;
  // outfilename2 = QAPATH + "TimeAbsRunLG" + fTrigger(r) + ".png" ;

  // if(SavePlots==2) c1a->SaveAs(outfilename);
  // if(SavePlots) c1a->SaveAs(outfilename2);
  // c1a->Write();
  // delete c1a;

  TCanvas* c1b = new TCanvas("EVsAbsId", "Cell Id vs. Energy", 600, 600);
  c1b->SetLogz();
  c1b->SetFillColor(0);
  c1b->SetBorderSize(0);
  c1b->SetFrameBorderMode(0);

  TH2F*  hEAbsId =(TH2F *)outputList->FindObject("EMCAL_hAmpId");
  if(!hEAbsId){
    Error(__FUNCTION__, Form("EMCAL_hAmpId: Histogram for trigger %s not found!", fTrigger.Data()));
    return -5;
  }
  FormatRunHisto(hEAbsId, Form("Cell ID vs. Energy%s", legend), "Cell AbsID");
  hEAbsId->RebinY(16);
  hEAbsId->SetMinimum(32);
  AutoZoom(hEAbsId, "maxx")->DrawCopy("colz");

  outfilename =  QAPATH + "EAbsIdRun" + fTrigger(r) + ".pdf" ;
  outfilename2 = QAPATH + "EAbsIdRun" + fTrigger(r) + ".png" ;

  if(SavePlots==2)
    c1b->SaveAs(outfilename);
  if(SavePlots) 
    c1b->SaveAs(outfilename2);
  c1b->Write();
  delete c1b;

  // TCanvas* c1bLG = new TCanvas("EVsAbsIdLG", "Cell Id vs. E LG", 600, 600);
  // c1bLG->SetLogz();
  // c1bLG->SetFillColor(0);
  // c1bLG->SetBorderSize(0);
  // c1bLG->SetFrameBorderMode(0);

  // TH2F*  hEAbsIdLG =(TH2F *)outputList->FindObject("EMCAL_hAmpIdLG");
  // if(!hEAbsIdLG){
  //   Error(__FUNCTION__, Form("EMCAL_hAmpIdLG: Histogram for trigger %s not found!", fTrigger.Data()));
  //   return -5;
  // }
  // FormatRunHisto(hEAbsIdLG, Form("Cell ID vs. Energy%s", legend), "Cell AbsID ");

  // AutoZoom(hEAbsIdLG, "maxx")->DrawCopy("colz");
  // outfilename =  QAPATH + "EAbsIdRunLG" + fTrigger(r) + ".pdf" ;
  // outfilename2 = QAPATH + "EAbsIdRunLG" + fTrigger(r) + ".png" ;

  // if(SavePlots==2) c1bLG->SaveAs(outfilename);
  // if(SavePlots) c1bLG->SaveAs(outfilename2);
  // c1bLG->Write();
  // delete c1bLG;

  if(pass!="muon_calo_pass1"){
    TCanvas  * c2 = new TCanvas("ClusterVsTrack ", "Correlation calo Mult vs. Track Multiplicity", 600, 600);
    c2->SetLogz();
    c2->SetFillColor(0);
    c2->SetBorderSize(0);
    c2->SetFrameBorderMode(0);

    TH2F* hClusterVsTrack =(TH2F *)outputList->FindObject("EMCAL_hCaloTrackMNClusters");
    FormatRunHisto(hClusterVsTrack, Form("#it{N} clusters vs. #it{N} tracks%s", legend), "#it{N} clusters");
    AutoZoom(hClusterVsTrack, "maxx, maxy", 1)->DrawCopy("colz");

    outfilename = QAPATH + "CaloTrackMult" + fTrigger(r) + ".pdf";
    outfilename2 = QAPATH + "CaloTrackMult" + fTrigger(r) + ".png";

    // if(pass!="muon_calo_pass1"){
    if(SavePlots==2)
      c2->SaveAs(outfilename);
    if(SavePlots)
      c2->SaveAs(outfilename2);
    c2->Write();
    delete c2;

    TCanvas* c3 = new TCanvas("ClusterEVsTrack ", "Correlation #it{E} calo vs. Track Multiplicity", 600, 600);
    c3->SetLogz();
    c3->SetFillColor(0);
    c3->SetBorderSize(0);
    c3->SetFrameBorderMode(0);

    TH2F* hClusterEVsTrack =(TH2F*)outputList->FindObject("EMCAL_hCaloTrackMEClusters");
    FormatRunHisto(hClusterEVsTrack, Form("Sum #it{E} clusters vs. #it{N} tracks%s", legend), "#Sigma#it{E} clusters");
    AutoZoom(hClusterEVsTrack, "maxx, maxy", 1)->DrawCopy("colz");

    outfilename =  QAPATH + "ETrackMult" + fTrigger(r) + ".pdf";
    outfilename2 = QAPATH + "ETrackMult" + fTrigger(r) + ".png";

    // if(pass!="muon_calo_pass1"){
    if(SavePlots==2)
      c3->SaveAs(outfilename);
    if(SavePlots)
      c3->SaveAs(outfilename2);
    c3->Write();
    // }
    delete c3;
  }

  TCanvas* c4 = new TCanvas("ClusterEVsV0 ", "Correlation #it{E} calo vs. V0 signal", 600, 600);
  c4->SetLogz();
  c4->SetFillColor(0);
  c4->SetBorderSize(0);
  c4->SetFrameBorderMode(0);

  TH2F* hClusterEVsV0S =(TH2F*)outputList->FindObject("EMCAL_hCaloV0SEClusters");
  FormatRunHisto(hClusterEVsV0S, Form("Sum #it{E} clusters vs. V0 signal%s", legend), "#Sigma#it{E} clusters");
  AutoZoom(hClusterEVsV0S, "maxx, maxy", 1)->DrawCopy("colz");

  outfilename = QAPATH +"EVsV0s" + fTrigger(r) + ".pdf";
  outfilename2 = QAPATH +"EVsV0s" + fTrigger(r) + ".png";

  if(SavePlots==2)
    c4->SaveAs(outfilename);
  if(SavePlots)
    c4->SaveAs(outfilename2);
  c4->Write();
  delete c4;

  TCanvas* c6 = new TCanvas("SumCellEvsSM", "sum cell energy vs. SM", 600, 600);
  c6->SetLogz();
  c6->SetFillColor(0);
  c6->SetBorderSize(0);
  c6->SetFrameBorderMode(0);

  TH2F* hSumCellEVsSM =(TH2F*)outputList->FindObject("EMCAL_hSumCellsAmp_Mod");
  FormatRunHisto(hSumCellEVsSM, Form("Sum #it{E} cells vs. SM%s", legend), "SM number");
  AutoZoom(hSumCellEVsSM, "maxx, maxy", 1)->DrawCopy("colz");

  outfilename = QAPATH +"SumCellEVsSM" + fTrigger(r) + ".pdf";
  outfilename2 = QAPATH +"SumCellEVsSM" + fTrigger(r) + ".png";

  if(SavePlots==2)
    c6->SaveAs(outfilename);
  if(SavePlots)
    c6->SaveAs(outfilename2);
  c6->Write();
  delete c6;

  TCanvas* c5 = new TCanvas("CellsperCluster", "Nb of cells per cluster for each SM", 600, 600);
  c5->SetLogz();
  c5->SetFillColor(0);
  c5->SetBorderSize(0);
  c5->SetFrameBorderMode(0);

  Bool_t mod3=0;
  if(nSM%3)
    mod3=1;
  c5->Divide(3, (nSM/3)+mod3);

  for(int ism = 0; ism < nSM; ism++){
    c5->cd(ism+1);
    gPad->SetLogz();

    if(TString(Form("Nb cells per clus%s Mod %d", legend, ism)).Length() > 60){
      Error(__FUNCTION__, "Title too long!");
      return -6;
    }

    AutoZoom(HistoPerMod((TH2F*)outputList->FindObject(Form("EMCAL_hNCellsPerCluster_Mod%i", ism)), Form("Nb of cells per cluster%s Mod %d", legend, ism)), "all", 1)->DrawCopy("colz");
  }

  outfilename =  QAPATH + "CellsperClusterSM" + fTrigger(r) + ".pdf";
  outfilename2 = QAPATH + "CellsperClusterSM" + fTrigger(r) + ".png";

  if(SavePlots==2)
    c5->SaveAs(outfilename);
  if(SavePlots)
    c5->SaveAs(outfilename2);
  c5->Write();
  delete c5;
  if(outputList)
    outputList->Delete();

  return 0;
}

//----------------------------------------------------------------------------------------------------------------------------------
Int_t TrendingEMCALTree(Long_t RunId, TString fCalorimeter, TString system, TString period, TString pass, int n, TList* TriggersList, TFile* f, TFile *fout, Int_t SavePlots){

  TString fTrigger="";
  TString aCalorimeter;
  if(n<=12)
    aCalorimeter = fCalorimeter;
  else
    aCalorimeter = TString("EMCAL_and_DCAL");
  
  TDatime now;

  Double_t Nevent=0 ;
  Double_t xe=0.5;
  Double_t NTClusters=0;    // total number of clusters
  Double_t NTClustersRMS=0;
  Double_t NCClusters=0;    // number of charged clusters
  Double_t NCClustersRMS=0;
  Double_t NMatchClustersP=0;
  Double_t NMatchClustersPRMS=0;

  Double_t CellMean=0;
  Double_t CellRMS=0;
  Double_t ClusterMean=0;
  Double_t ClusterRMS=0;
  Double_t EtotalMean=0;
  Double_t EtotalRMS=0;

  Double_t CellPerClusterMean=0;
  Double_t CellPerClusterRMS=0;

  Double_t mPDG = 134.9766;
  Double_t Npi0EMCAL=0;
  Double_t Npi0EMCALErr=0;
  Double_t MeanPosEMCAL=0;
  Double_t MeanPosEMCALErr=0;
  Double_t WidthEMCAL=0;
  Double_t WidthEMCALErr=0;
  Double_t Chi2NdfPi0EMCAL=0;
  Double_t NggEMCAL=0;
  Double_t NggEMCALErr=0;
  Double_t SignifEMCAL=0;    // !S/(S+B)
  Double_t SignifEMCALErr=0; // !S/(S+B)

  Double_t Npi0DCAL=0;
  Double_t Npi0DCALErr=0;
  Double_t MeanPosDCAL=0;
  Double_t MeanPosDCALErr=0;
  Double_t WidthDCAL=0;
  Double_t WidthDCALErr=0;
  Double_t Chi2NdfPi0DCAL=0;
  Double_t NggDCAL=0;
  Double_t NggDCALErr=0;
  Double_t SignifDCAL=0;    // !S/(S+B)
  Double_t SignifDCALErr=0; // !S/(S+B)

  TFile* ftree = new TFile(Form("%s/trending.root", QAPATH.Data()), "RECREATE");

  TTree *tree = new TTree("trending", "Trending QA Tree");
  tree->Branch("fDate", &now);
  tree->Branch("fCalorimeter", &aCalorimeter);
  tree->Branch("system", &system);
  tree->Branch("period", &period);
  tree->Branch("pass", &pass);
  tree->Branch("fTrigger", &fTrigger);
  tree->Branch("run", &RunId, "run/I");
  tree->Branch("xe", &xe, "xe/D");

  tree->Branch("Nevent", &Nevent, "Nevent/D");
  tree->Branch("NMatchClustersP", &NMatchClustersP, "NMatchClustersP/D");
  tree->Branch("NMatchClustersPRMS", &NMatchClustersPRMS, "NMatchClustersPRMS/D");
  tree->Branch("CellMean", &CellMean, "CellMean/D");
  tree->Branch("CellRMS", &CellRMS, "CellRMS/D");
  tree->Branch("ClusterMean", &ClusterMean, "ClusterMean/D");
  tree->Branch("ClusterRMS", &ClusterRMS, "ClusterRMS/D");
  tree->Branch("EtotalMean", &EtotalMean, "EtotalMean/D");
  tree->Branch("EtotalRMS", &EtotalRMS, "EtotalRMS/D");

  tree->Branch("CellPerClusterMean", &CellPerClusterMean, "CellPerClusterMean/D");
  tree->Branch("CellPerClusterRMS", &CellPerClusterRMS, "CellPerClusterRMS/D");

  tree->Branch("Npi0EMCAL", &Npi0EMCAL, "Npi0EMCAL/D");
  tree->Branch("Npi0EMCALErr", &Npi0EMCALErr, "Npi0EMCALErr/D");
  tree->Branch("Npi0DCAL", &Npi0DCAL, "Npi0DCAL/D");
  tree->Branch("Npi0DCALErr", &Npi0DCALErr, "Npi0DCALErr/D");
  tree->Branch("MeanPosEMCAL", &MeanPosEMCAL, "MeanPosEMCAL/D");
  tree->Branch("MeanPosEMCALErr", &MeanPosEMCALErr, "MeanPosEMCALErr/D");
  tree->Branch("MeanPosDCAL", &MeanPosDCAL, "MeanPosDCAL/D");
  tree->Branch("MeanPosDCALErr", &MeanPosDCALErr, "MeanPosDCALErr/D");
  tree->Branch("WidthEMCAL", &WidthEMCAL, "WidthEMCAL/D");
  tree->Branch("WidthEMCALErr", &WidthEMCALErr, "WidthEMCALErr/D");
  tree->Branch("WidthDCAL", &WidthDCAL, "WidthDCAL/D");
  tree->Branch("WidthDCALErr", &WidthDCALErr, "WidthDCALErr/D");
  tree->Branch("Chi2NdfPi0EMCAL", &Chi2NdfPi0EMCAL, "Chi2NdfPi0EMCAL/D");
  tree->Branch("Chi2NdfPi0DCAL", &Chi2NdfPi0DCAL, "Chi2NdfPi0DCAL/D");
  tree->Branch("NggEMCAL", &NggEMCAL, "NggEMCAL/D");
  tree->Branch("NggEMCALErr", &NggEMCALErr, "NggEMCALErr/D");
  tree->Branch("NggDCAL", &NggDCAL, "NggDCAL/D");
  tree->Branch("NggDCALErr", &NggDCALErr, "NggDCALErr/D");
  tree->Branch("SignifEMCAL", &SignifEMCAL, "SignifEMCAL/D");
  tree->Branch("SignifDCALErr", &SignifDCALErr, "SignifDCALErr/D");
  tree->Branch("SignifEMCAL", &SignifEMCAL, "SignifEMCAL/D");
  tree->Branch("SignifDCALErr", &SignifDCALErr, "SignifDCALErr/D");

  tree->Branch("nSM", &n, "nSM/I");

  int nMax = 22;
  Double_t CellMeanSM[nMax];
  Double_t CellRMSSM[nMax];
  Double_t ClusterMeanSM[nMax];
  Double_t ClusterTotSM[nMax];
  Double_t ClusterRMSSM[nMax];
  Double_t EtotalMeanSM[nMax];         // mean total energy deposited per event
  Double_t EtotalRMSSM[nMax];
  Double_t CellPerClusterMeanSM[nMax];
  Double_t CellPerClusterRMSSM[nMax];
  Double_t ECell1MeanSM[nMax];         // total energy deposited per event without 1-cell clusters
  Double_t ECell1RMSSM[nMax];

  Double_t MeanPosSM[nMax];
  Double_t MeanPosErrSM[nMax];
  Double_t WidthSM[nMax];
  Double_t WidthErrSM[nMax];
  Double_t Npi0SM[nMax];
  Double_t Npi0ErrSM[nMax];

  tree->Branch("CellMeanSM", CellMeanSM, TString::Format("CellMeanSM[%i]/D", nMax));
  tree->Branch("CellRMSSM", CellRMSSM, TString::Format("CellRMSSM[%i]/D", nMax));
  tree->Branch("ClusterMeanSM", ClusterMeanSM, TString::Format("ClusterMeanSM[%i]/D", nMax));
  tree->Branch("ClusterTotSM", ClusterTotSM, TString::Format("ClusterTotSM[%i]/D", nMax));
  tree->Branch("ClusterRMSSM", ClusterRMSSM, TString::Format("ClusterRMSSM[%i]/D", nMax));
  tree->Branch("EtotalMeanSM", EtotalMeanSM, TString::Format("EtotalMeanSM[%i]/D", nMax));
  tree->Branch("EtotalRMSSM", EtotalRMSSM, TString::Format("EtotalRMSSM[%i]/D", nMax));
  tree->Branch("CellPerClusterMeanSM", CellPerClusterMeanSM, TString::Format("CellPerClusterMeanSM[%i]/D", nMax));
  tree->Branch("CellPerClusterRMSSM", CellPerClusterRMSSM, TString::Format("CellPerClusterRMSSM[%i]/D", nMax));
  tree->Branch("ECell1MeanSM", ECell1MeanSM, TString::Format("ECell1MeanSM[%i]/D", nMax));
  tree->Branch("ECell1RMSSM", ECell1RMSSM, TString::Format("ECell1RMSSM[%i]/D", nMax));

  tree->Branch("MeanPosSM", MeanPosSM, TString::Format("MeanPosSM[%i]/D", nMax));
  tree->Branch("MeanPosErrSM", MeanPosErrSM, TString::Format("MeanPosErrSM[%i]/D", nMax));
  tree->Branch("WidthSM", WidthSM, TString::Format("WidthSM[%i]/D", nMax));
  tree->Branch("WidthErrSM", WidthErrSM, TString::Format("WidthErrSM[%i]/D", nMax));
  tree->Branch("Npi0SM", Npi0SM, TString::Format("Npi0SM[%i]/D", nMax));
  tree->Branch("Npi0ErrSM", Npi0ErrSM, TString::Format("Npi0ErrSM[%i]/D", nMax));

  TF1* fitMass = new TF1("fitMass", pi0massP2, 100, 250, 6);
  fitMass->SetParName(0, "A");
  fitMass->SetParName(1, "m_{0}");
  fitMass->SetParName(2, "sigma");
  fitMass->SetParName(3, "a_{0}");
  fitMass->SetParName(4, "a_{1}");
  fitMass->SetParName(5, "a_{2}");
  fitMass->SetParLimits(0, 1.e-5, 1.e5);
  fitMass->SetParLimits(1, 0.11, 0.16);
  fitMass->SetParLimits(2, 0.001, 0.06);

  TList* outputList;

  TH1F* fhNEvents;
  TH1F* fhE;
  TH1F* fhNClusters = 0x0;
  TH1F* fhECharged;
  TH1F* fhNCells = 0x0;

  TH1F* NCells[n];
  TH1F* NClusters[n];
  TH2F* NCellsPerCluster[n];
  TH1F* E[n];

  TH2F* fhIM;
  TH2F* fhIMDCAL;
  TH1F* fhMggEMCAL;
  TH1F* fhMggDCAL;
  TH2F* IM[n];
  TH1F* MggSM[n];

  TPMERegexp r("_\\w+");
  TIter next(TriggersList);
  int ret = 0;

  while(TObject *obj = next()){
    fTrigger= TString(obj->GetName());
    TString namefile = QAPATH + period + "_" +  pass + fTrigger(r).Data() + "_" + RunId + "_data.txt";
    ofstream QAData(namefile, ios::app); // write checks at the end

    NTClusters=0;
    NTClustersRMS=0;
    NCClusters=0;
    NCClustersRMS=0;
    NMatchClustersP=0;
    NMatchClustersPRMS=0;
    CellMean=0;
    CellRMS=0;
    ClusterMean=0;
    ClusterRMS=0;
    EtotalMean=0;
    EtotalRMS=0;
    CellPerClusterMean=0;
    CellPerClusterRMS=0;

    // EMCal
    Npi0EMCAL=0;
    Npi0EMCALErr=0;
    MeanPosEMCAL=0;
    MeanPosEMCALErr=0;
    WidthEMCAL=0;
    WidthEMCALErr=0;
    Chi2NdfPi0EMCAL=0;
    NggEMCAL=0;
    NggEMCALErr=0;
    SignifEMCAL=0;
    SignifEMCALErr=0;

    // DCal
    Npi0DCAL=0;
    Npi0DCALErr=0;
    MeanPosDCAL=0;
    MeanPosDCALErr=0;
    WidthDCAL=0;
    WidthDCALErr=0;
    Chi2NdfPi0DCAL=0;
    NggDCAL=0;
    NggDCALErr=0;
    SignifDCAL=0;
    SignifDCALErr=0;

    memset (CellMeanSM, 0, sizeof (Double_t) * nMax);
    memset (CellRMSSM, 0, sizeof (Double_t) *  nMax);
    memset (ClusterMeanSM, 0, sizeof (Double_t) * nMax);
    memset (ClusterTotSM, 0, sizeof (Double_t) * nMax);
    memset (ClusterRMSSM, 0, sizeof (Double_t) * nMax);
    memset (EtotalMeanSM, 0, sizeof (Double_t) * nMax);
    memset (EtotalRMSSM, 0, sizeof (Double_t) * nMax);
    memset (CellPerClusterMeanSM, 0, sizeof (Double_t) * nMax);
    memset (CellPerClusterRMSSM, 0, sizeof (Double_t) * nMax);
    memset (ECell1MeanSM, 0, sizeof (Double_t) * nMax);
    memset (ECell1RMSSM, 0, sizeof (Double_t) * nMax);

    memset (MeanPosSM, 0, sizeof (Double_t) * nMax);
    memset (MeanPosErrSM, 0, sizeof (Double_t) * nMax);
    memset (WidthSM, 0, sizeof (Double_t) * nMax);
    memset (WidthErrSM, 0, sizeof (Double_t) * nMax);
    memset (Npi0SM, 0, sizeof (Double_t) * nMax);
    memset (Npi0ErrSM, 0, sizeof (Double_t) * nMax);

    TString dirname;
    if(!fTrigger.Contains("QA")){
      dirname = "CaloQA_";
    }
    dirname += fTrigger;

    Bool_t dirok = f->cd(dirname);
    if(!dirok){
      Error(__FUNCTION__, Form("No input directory %s", dirname.Data()));
      tree->Fill(); ftree->cd();
      tree->Write();
      ret= -1;
      continue;
    }

    outputList = (TList*)gDirectory->Get(dirname);
    if(!outputList){
      Error(__FUNCTION__, Form("No input list! %s", dirname.Data()));
      tree->Fill();  ftree->cd();
      tree->Write();
      ret=-2;
      continue;;
    }
    outputList->SetOwner();

    // number of events
    fhNEvents =(TH1F *)outputList->FindObject("hNEvents");
    if(!fhNEvents){
      Error(__FUNCTION__, Form("NEvent histogram not found for trigger %s", fTrigger.Data()));
      tree->Fill();  ftree->cd();
      tree->Write();  ret=-3;
      continue;
    }

    Nevent=fhNEvents->GetEntries();

    if(Nevent==0){
      Error(__FUNCTION__, Form("No event in trigger %s", fTrigger.Data()));
      tree->Fill();  ftree->cd();
      tree->Write();
      ret=-4;
      continue;
    }
    if(Nevent<20){
      Error(__FUNCTION__, Form("Less than 20 events in trigger %s", fTrigger.Data()));
      tree->Fill();  ftree->cd();
      tree->Write();
      ret=-5;
      continue;
    }

    TH1F* hE =(TH1F *)outputList->FindObject("EMCAL_hE");
    if(!hE || (!hE->GetEntries())){
      Error(__FUNCTION__, Form("hE Histogram not found or it is empty for trigger %s, EMCal was not in this run?", fTrigger.Data())); EtotalMean=0;
      tree->Fill();  ftree->cd();
      tree->Write();
      ret=-4;
      continue;
    }

    // first do clusters trending
    fhE = (TH1F *)outputList->FindObject(fCalorimeter+"_hE");
    fhECharged = (TH1F *)outputList->FindObject(fCalorimeter+"_hECharged");

    Double_t energy = 0. ;

    for(Int_t ibin = fhE->FindBin(0.6) ; ibin <fhE->FindBin(50.) ; ibin++){
      energy+=fhE->GetBinCenter(ibin)*fhE->GetBinContent(ibin);
    }

    if(fhE->Integral(fhE->FindBin(0.6), fhE->FindBin(50.))==0){
      Error(__FUNCTION__, Form("Not enough events"));
      tree->Fill(); ftree->cd();
      tree->Write();
      ret=-6;
      continue;
    }

    EtotalMean=energy/fhE->Integral(fhE->FindBin(0.6), fhE->FindBin(50.)) ;
    EtotalRMS=fhE->GetMeanError();

    TString nameNCell = Form("%s_hNCells_Mod", fCalorimeter.Data());
    TString nameNCluster = Form("%s_hNClusters_Mod", fCalorimeter.Data());
    TString nameE = Form("%s_hE_Mod", fCalorimeter.Data());
    TH2F* hNCellsMod= (TH2F*)outputList->FindObject(nameNCell.Data());
    TH2F* hNClusterMod=(TH2F*)outputList->FindObject(nameNCluster.Data());
    TH2F* hEMod=(TH2F*)outputList->FindObject(nameE.Data());

    if(!hNCellsMod || !hNClusterMod || !hEMod){
      Error(__FUNCTION__, "A requiered histogram was not found (the imput QAresult.root might be too old)!");
      tree->Fill(); ftree->cd();
      tree->Write();
      ret=-7;
      continue;
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // EM addition (Nov. 2016) to observe the PAR Time issue
    TString trig_title = "";
    if(fTrigger.Contains("EMC"))
      trig_title += "EMC";
    else
      trig_title += "MB";

    // BC/4 = 0
    TCanvas* c_TimeSM_BC0 = new TCanvas("TimeSM_BC0", "Cell Time (for BC/4=0) for each SM", 1200, 1200);
    c_TimeSM_BC0->SetFillColor(0);
    c_TimeSM_BC0->SetBorderSize(0);
    c_TimeSM_BC0->SetFrameBorderMode(0);
    Bool_t mod3_BC0=0;
    if(n%3)
      mod3_BC0=1;
    c_TimeSM_BC0->Divide(3, n/3+mod3_BC0);

    TH2F* hTimeSM_BC0 = (TH2F*)outputList->FindObject("EMCAL_hTimePerSM_BC0");
    TH1F* hTimeSM_BC0_proj[n];
    TString TimeSM_BC0_projName;

    if(hTimeSM_BC0){
      cout << "---------------------" << endl;
      cout << " BC/4=0 | Run " << RunId << endl;
      for(Int_t ism = 0 ; ism < n ; ism++){
	c_TimeSM_BC0->cd(ism+1);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.15);
	// gPad->SetLogy();

	TimeSM_BC0_projName = Form("SM_%d", ism);
	hTimeSM_BC0_proj[ism] = (TH1F *)hTimeSM_BC0->ProjectionX(TimeSM_BC0_projName.Data(), ism+1, ism+1, "") ;
	hTimeSM_BC0_proj[ism]->Draw("hist");

	hTimeSM_BC0_proj[ism]->SetTitle(Form("#it{t}_{cell} for super module %i (BC/4=0) %s", ism, trig_title.Data()));
	hTimeSM_BC0_proj[ism]->SetTitleSize(0.1);
	hTimeSM_BC0_proj[ism]->SetXTitle("#it{t}_{cell} (ns)");
	hTimeSM_BC0_proj[ism]->SetYTitle("Nb of entries");
	hTimeSM_BC0_proj[ism]->SetStats(kTRUE);

	hTimeSM_BC0_proj[ism]->GetXaxis()->SetRangeUser(400, 900);
	hTimeSM_BC0_proj[ism]->GetXaxis()->SetLabelSize(0.05);
	hTimeSM_BC0_proj[ism]->GetXaxis()->SetLabelOffset(0.015);
	hTimeSM_BC0_proj[ism]->GetXaxis()->SetTitleSize(0.075);
	hTimeSM_BC0_proj[ism]->GetXaxis()->SetTitleOffset(0.9);

	hTimeSM_BC0_proj[ism]->GetYaxis()->SetLabelSize(0.05);
	hTimeSM_BC0_proj[ism]->GetYaxis()->SetTitleSize(0.06);
	hTimeSM_BC0_proj[ism]->GetYaxis()->SetTitleOffset(0.6);
      }

      TString outfilename_BC0 = QAPATH + "TimeSM_BC0" + fTrigger(r) + ".pdf";
      TString outfilename2_BC0 = QAPATH + "TimeSM_BC0" + fTrigger(r) + ".png";
      if(SavePlots==2)
	c_TimeSM_BC0->SaveAs(outfilename_BC0);
      if(SavePlots)
	c_TimeSM_BC0->SaveAs(outfilename2_BC0);
    }

    // BC/4 = 1
    TCanvas* c_TimeSM_BC1 = new TCanvas("TimeSM_BC1", "Cell Time (for BC/4=1) for each SM", 1200, 1200);
    c_TimeSM_BC1->SetFillColor(0);
    c_TimeSM_BC1->SetBorderSize(0);
    c_TimeSM_BC1->SetFrameBorderMode(0);
    Bool_t mod3_BC1=0;
    if(n%3)
      mod3_BC1=1;
    c_TimeSM_BC1->Divide(3, n/3+mod3_BC1);

    TH2F* hTimeSM_BC1 = (TH2F*)outputList->FindObject("EMCAL_hTimePerSM_BC1");
    TH1F* hTimeSM_BC1_proj[n];
    TString TimeSM_BC1_projName;

    if(hTimeSM_BC1){
      cout << "---------------------" << endl;
      cout << " BC/4=1 | Run " << RunId << endl;
      for(Int_t ism = 0 ; ism < n ; ism++){
	c_TimeSM_BC1->cd(ism+1);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.15);
	// gPad->SetLogy();

	TimeSM_BC1_projName = Form("SM_%d", ism);
	hTimeSM_BC1_proj[ism] = (TH1F *)hTimeSM_BC1->ProjectionX(TimeSM_BC1_projName.Data(), ism+1, ism+1, "") ;
	hTimeSM_BC1_proj[ism]->Draw("hist");

	hTimeSM_BC1_proj[ism]->SetTitle(Form("#it{t}_{cell} for super module %i (BC/4=1) %s", ism, trig_title.Data()));
	hTimeSM_BC1_proj[ism]->SetTitleSize(0.1);
	hTimeSM_BC1_proj[ism]->SetXTitle("#it{t}_{cell} (ns)");
	hTimeSM_BC1_proj[ism]->SetYTitle("Nb of entries");
	hTimeSM_BC1_proj[ism]->SetStats(kTRUE);

	hTimeSM_BC1_proj[ism]->GetXaxis()->SetRangeUser(400, 900);
	hTimeSM_BC1_proj[ism]->GetXaxis()->SetLabelSize(0.05);
	hTimeSM_BC1_proj[ism]->GetXaxis()->SetLabelOffset(0.015);
	hTimeSM_BC1_proj[ism]->GetXaxis()->SetTitleSize(0.075);
	hTimeSM_BC1_proj[ism]->GetXaxis()->SetTitleOffset(0.9);

	hTimeSM_BC1_proj[ism]->GetYaxis()->SetLabelSize(0.05);
	hTimeSM_BC1_proj[ism]->GetYaxis()->SetTitleSize(0.06);
	hTimeSM_BC1_proj[ism]->GetYaxis()->SetTitleOffset(0.6);
      }

      TString outfilename_BC1 = QAPATH + "TimeSM_BC1" + fTrigger(r) + ".pdf";
      TString outfilename2_BC1 = QAPATH + "TimeSM_BC1" + fTrigger(r) + ".png";
      if(SavePlots==2)
	c_TimeSM_BC1->SaveAs(outfilename_BC1);
      if(SavePlots)
	c_TimeSM_BC1->SaveAs(outfilename2_BC1);
    }

    // BC/4 = 2
    TCanvas* c_TimeSM_BC2 = new TCanvas("TimeSM_BC2", "Cell Time (for BC/4=2) for each SM", 1200, 1200);
    c_TimeSM_BC2->SetFillColor(0);
    c_TimeSM_BC2->SetBorderSize(0);
    c_TimeSM_BC2->SetFrameBorderMode(0);
    Bool_t mod3_BC2=0;
    if(n%3)
      mod3_BC2=1;
    c_TimeSM_BC2->Divide(3, n/3+mod3_BC2);

    TH2F* hTimeSM_BC2 = (TH2F*)outputList->FindObject("EMCAL_hTimePerSM_BC2");
    TH1F* hTimeSM_BC2_proj[n];
    TString TimeSM_BC2_projName;

    if(hTimeSM_BC2){
      cout << "---------------------" << endl;
      cout << " BC/4=2 | Run " << RunId << endl;
      for(Int_t ism = 0 ; ism < n ; ism++){
	c_TimeSM_BC2->cd(ism+1);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.15);
	// gPad->SetLogy();

	TimeSM_BC2_projName = Form("SM_%d", ism);
	hTimeSM_BC2_proj[ism] = (TH1F *)hTimeSM_BC2->ProjectionX(TimeSM_BC2_projName.Data(), ism+1, ism+1, "") ;
	hTimeSM_BC2_proj[ism]->Draw("hist");

	hTimeSM_BC2_proj[ism]->SetTitle(Form("#it{t}_{cell} for super module %i (BC/4=2) %s", ism, trig_title.Data()));
	hTimeSM_BC2_proj[ism]->SetTitleSize(0.1);
	hTimeSM_BC2_proj[ism]->SetXTitle("#it{t}_{cell} (ns)");
	hTimeSM_BC2_proj[ism]->SetYTitle("Nb of entries");
	hTimeSM_BC2_proj[ism]->SetStats(kTRUE);

	hTimeSM_BC2_proj[ism]->GetXaxis()->SetRangeUser(400, 900);
	hTimeSM_BC2_proj[ism]->GetXaxis()->SetLabelSize(0.05);
	hTimeSM_BC2_proj[ism]->GetXaxis()->SetLabelOffset(0.015);
	hTimeSM_BC2_proj[ism]->GetXaxis()->SetTitleSize(0.075);
	hTimeSM_BC2_proj[ism]->GetXaxis()->SetTitleOffset(0.9);

	hTimeSM_BC2_proj[ism]->GetYaxis()->SetLabelSize(0.05);
	hTimeSM_BC2_proj[ism]->GetYaxis()->SetTitleSize(0.06);
	hTimeSM_BC2_proj[ism]->GetYaxis()->SetTitleOffset(0.6);
      }

      TString outfilename_BC2 = QAPATH + "TimeSM_BC2" + fTrigger(r) + ".pdf";
      TString outfilename2_BC2 = QAPATH + "TimeSM_BC2" + fTrigger(r) + ".png";
      if(SavePlots==2)
	c_TimeSM_BC2->SaveAs(outfilename_BC2);
      if(SavePlots)
	c_TimeSM_BC2->SaveAs(outfilename2_BC2);
    }

    // BC/4 = 3
    TCanvas* c_TimeSM_BC3 = new TCanvas("TimeSM_BC3", "Cell Time (for BC/4=3) for each SM", 1200, 1200);
    c_TimeSM_BC3->SetFillColor(0);
    c_TimeSM_BC3->SetBorderSize(0);
    c_TimeSM_BC3->SetFrameBorderMode(0);
    Bool_t mod3_BC3=0;
    if(n%3)
      mod3_BC3=1;
    c_TimeSM_BC3->Divide(3, n/3+mod3_BC3);

    TH2F* hTimeSM_BC3 = (TH2F*)outputList->FindObject("EMCAL_hTimePerSM_BC3");
    TH1F* hTimeSM_BC3_proj[n];
    TString TimeSM_BC3_projName;

    if(hTimeSM_BC3){
      cout << "---------------------" << endl;
      cout << " BC/4=3 | Run " << RunId << endl;
      for(Int_t ism = 0 ; ism < n ; ism++){
	c_TimeSM_BC3->cd(ism+1);
	gPad->SetTopMargin(0.1);
	gPad->SetBottomMargin(0.15);
	// gPad->SetLogy();

	TimeSM_BC3_projName = Form("SM_%d", ism);
	hTimeSM_BC3_proj[ism] = (TH1F *)hTimeSM_BC3->ProjectionX(TimeSM_BC3_projName.Data(), ism+1, ism+1, "") ;
	hTimeSM_BC3_proj[ism]->Draw("hist");

	hTimeSM_BC3_proj[ism]->SetTitle(Form("#it{t}_{cell} for super module %i (BC/4=3) %s", ism, trig_title.Data()));
	hTimeSM_BC3_proj[ism]->SetTitleSize(0.1);
	hTimeSM_BC3_proj[ism]->SetXTitle("#it{t}_{cell} (ns)");
	hTimeSM_BC3_proj[ism]->SetYTitle("Nb of entries");
	hTimeSM_BC3_proj[ism]->SetStats(kTRUE);

	hTimeSM_BC3_proj[ism]->GetXaxis()->SetRangeUser(400, 900);
	hTimeSM_BC3_proj[ism]->GetXaxis()->SetLabelSize(0.05);
	hTimeSM_BC3_proj[ism]->GetXaxis()->SetLabelOffset(0.015);
	hTimeSM_BC3_proj[ism]->GetXaxis()->SetTitleSize(0.075);
	hTimeSM_BC3_proj[ism]->GetXaxis()->SetTitleOffset(0.9);

	hTimeSM_BC3_proj[ism]->GetYaxis()->SetLabelSize(0.05);
	hTimeSM_BC3_proj[ism]->GetYaxis()->SetTitleSize(0.06);
	hTimeSM_BC3_proj[ism]->GetYaxis()->SetTitleOffset(0.6);
      }

      TString outfilename_BC3 = QAPATH + "TimeSM_BC3" + fTrigger(r) + ".pdf";
      TString outfilename2_BC3 = QAPATH + "TimeSM_BC3" + fTrigger(r) + ".png";
      if(SavePlots==2)
	c_TimeSM_BC3->SaveAs(outfilename_BC3);
      if(SavePlots)
	c_TimeSM_BC3->SaveAs(outfilename2_BC3);
    }

    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    TCanvas* c1 = new TCanvas("Pi0InvMassSM", "Pi0 Invariant Mass for each SM", 600, 600);
    gStyle->SetOptStat(1); // MG modif
    gStyle->SetOptFit(1); // MG modif
    c1->SetFillColor(0);
    c1->SetBorderSize(0);
    c1->SetFrameBorderMode(0);
    // c1->SetOptStat(1); // MG modif

    Bool_t mod3=0;
    if(n%3)
      mod3=1;
    c1->Divide(3, n/3+mod3);

    // per sm trending
    TString nameNCellPerCluster;
    for(Int_t ism = 0 ; ism < n ; ism++){
      cout << "#########################"<< endl;
      cout      << " Super Module " << ism << " Run " << RunId << endl;

      // first do clusters trending
      nameNCellPerCluster = Form("%s_hNCellsPerCluster_Mod%d", fCalorimeter.Data(), ism);
      NCellsPerCluster[ism] = (TH2F*)outputList->FindObject(nameNCellPerCluster.Data());

      if(!  (TH2F*)outputList->FindObject(nameNCellPerCluster.Data()) ){
	Error(__FUNCTION__, Form("NCellsPerCluster histogram not found for super module %i", ism));
	ret=-8;
	continue;
      }

      NCellsPerCluster[ism] = (TH2F*)outputList->FindObject(nameNCellPerCluster.Data());
      NCells[ism] = (TH1F*)hNCellsMod->ProjectionX(Form("NCells%d", ism), ism+1, ism+1, "");
      NClusters[ism] = (TH1F*)hNClusterMod->ProjectionX(Form("NClusters%d", ism), ism+1, ism+1, "");
      E[ism] = (TH1F*)hEMod->ProjectionX(Form("E%d", ism), ism+1, ism+1, "");

      CellMeanSM[ism]=NCells[ism]->GetMean();
      CellRMSSM[ism]=NCells[ism]->GetMeanError();

      Int_t binmax;
      binmax =  NClusters[ism]->GetNbinsX();
      for(Int_t ibin = 1 ;ibin <=binmax ;ibin++){
	ClusterTotSM[ism]+=(NClusters[ism]->GetBinCenter(ibin))*(NClusters[ism]->GetBinContent(ibin));
      }
      if(NClusters[ism]->GetEntries()>0)
	ClusterTotSM[ism]= ClusterTotSM[ism]/NClusters[ism]->GetEntries();

      // ClusterTotSM[ism]=NClusters[ism]->GetMean();
      ClusterMeanSM[ism]=NClusters[ism]->GetMean();
      ClusterRMSSM[ism]=NClusters[ism]->GetMeanError();
      CellPerClusterMeanSM[ism]=NCellsPerCluster[ism]->GetMean(2);
      CellPerClusterRMSSM[ism]=NCellsPerCluster[ism]->GetMeanError(2);

      ECell1MeanSM[ism] =NCellsPerCluster[ism]->ProjectionX("", 2, 50, "")->Integral(5, 50)/(Nevent);
      ECell1RMSSM[ism] =NCellsPerCluster[ism]->ProjectionX("", 2, 50, "")->GetMeanError();

      Double_t energySM = 0. ;
      for(Int_t ibin = E[ism]->FindBin(0.6) ; ibin <E[ism]->FindBin(50.) ; ibin++){
	energySM+=E[ism]->GetBinCenter(ibin)*(E[ism]->GetBinContent(ibin));
      }

      if(E[ism]->Integral(E[ism]->FindBin(0.6), E[ism]->FindBin(50.))==0){
	Error(__FUNCTION__, Form("Energy: Not enough events/SM"));
	continue;
      }

      EtotalMeanSM[ism]=energySM/(E[ism]->Integral(E[ism]->FindBin(0.6), E[ism]->FindBin(50.)));
      EtotalRMSSM[ism]=E[ism]->GetMeanError();

      if(ism==0){
	fhNCells = (TH1F*)NCells[ism]->Clone("NCells");
	fhNClusters = (TH1F*)NClusters[ism]->Clone("NClusters");
      }
      else{
	fhNCells->Add(NCells[ism], 1);
	fhNClusters->Add(NClusters[ism], 1);
      }

      // Pi0
      c1->cd(ism+1);
      TString namePair = Form("%s_hIM_Mod%d", fCalorimeter.Data(), ism);
      IM[ism] = (TH2F*)outputList->FindObject(namePair.Data());
      IM[ism]->Sumw2();

      TString projname = Form("SM_%d", ism);
      MggSM[ism] = (TH1F *)IM[ism]->ProjectionY(projname.Data(), 2, 20, "") ; // MG modif

      if(MggSM[ism]->GetEntries()>100){
	fitMass->SetParameter(0, MggSM[ism]->GetBinContent(MggSM[ism]->GetMaximumBin()));
	fitMass->SetParameter(1, mPDG);
	fitMass->SetParameter(2, 15.);
	fitMass->SetParameter(3, 0.);
	fitMass->SetParameter(4, MggSM[ism]->GetBinContent(MggSM[ism]->FindBin(0.11)));
	fitMass->SetParameter(5, MggSM[ism]->GetBinContent(MggSM[ism]->FindBin(0.20)));

	if(MggSM[ism]->GetEntries()<1000)
	  MggSM[ism]->Rebin(5);
	else
	  MggSM[ism]->Rebin();
	
	MggSM[ism]->Fit("fitMass", "WL R +", "", 0.1, 0.30);
	MggSM[ism]->SetTitle(Form("Pi0 Mass for super module %i", ism));
	MggSM[ism]->SetTitleSize(0.1);
	MggSM[ism]->SetXTitle("Pi0 Mass");
	MggSM[ism]->SetYTitle("Nb of entries");
	MggSM[ism]->GetXaxis()->SetLabelSize(0.05);
	MggSM[ism]->GetXaxis()->SetTitleSize(0.07);
	MggSM[ism]->GetXaxis()->SetTitleOffset(0.68);
	MggSM[ism]->GetYaxis()->SetLabelSize(0.05);
	MggSM[ism]->GetYaxis()->SetTitleSize(0.06);
	MggSM[ism]->GetYaxis()->SetTitleOffset(0.78);

	MeanPosSM[ism] = MggSM[ism]->GetFunction("fitMass")->GetParameter(1)*1000;
	MeanPosErrSM[ism] = MggSM[ism]->GetFunction("fitMass")->GetParError(1)*1000;
	WidthSM[ism] = MggSM[ism]->GetFunction("fitMass")->GetParameter(2)*1000;
	WidthErrSM[ism] = MggSM[ism]->GetFunction("fitMass")->GetParError(2)*1000;
	Npi0SM[ism] = MggSM[ism]->GetFunction("fitMass")->GetParameter(0)*(MggSM[ism]->GetFunction("fitMass")->GetParameter(2))*TMath::Sqrt(2*TMath::Pi())/(Nevent*MggSM[ism]->GetBinWidth(1));
	Npi0ErrSM[ism] = TMath::Sqrt((MggSM[ism]->GetFunction("fitMass")->GetParError(0)/MggSM[ism]->GetFunction("fitMass")->GetParameter(0))*(MggSM[ism]->GetFunction("fitMass")->GetParError(0)/MggSM[ism]->GetFunction("fitMass")->GetParameter(0))
				     +(MggSM[ism]->GetFunction("fitMass")->GetParError(2)/MggSM[ism]->GetFunction("fitMass")->GetParameter(2))*(MggSM[ism]->GetFunction("fitMass")->GetParError(2)/MggSM[ism]->GetFunction("fitMass")->GetParameter(2)));
	Npi0ErrSM[ism] = 0.;

      } // end if enough events for Pi0 fit and trending
      else
	Info(__FUNCTION__, Form("Not enough events for Pi0 fit and trending for super module %i", ism));
    } // per SM loop

      // Now Pi0 global trending in EMCAL
    TCanvas* c2 = new TCanvas("Pi0InvMassEMCAL", "Pi0 Invariant Mass in EMCal", 600, 600);
    c2->SetFillColor(0);
    c2->SetBorderSize(0);
    c2->SetFrameBorderMode(0);
    // c2->SetOptStat(1); // MG modif

    fhIM = (TH2F *)outputList->FindObject(fCalorimeter+"_hIM");
    fhIM->Sumw2();
    fhMggEMCAL = (TH1F *)fhIM->ProjectionY("MggEMCAL", 2, 20, "") ; // to modify projection range

    if(fhMggEMCAL->GetEntries()==0){
      Error(__FUNCTION__, "The Pi0 histogram in EMCal is empty !");  tree->Fill();
      ret=-8;
      continue;
    }

    fitMass->SetParameter(0, 4500);
    fitMass->SetParameter(1, mPDG);
    fitMass->SetParameter(2, 0.01);
    fitMass->SetParameter(3, 0.);
    fitMass->SetParameter(4, fhMggEMCAL->GetBinContent(fhMggEMCAL->FindBin(0.11)));
    fitMass->SetParameter(5, fhMggEMCAL->GetBinContent(fhMggEMCAL->FindBin(0.20)));

    if(fhMggEMCAL->GetEntries()<5000)
      fhMggEMCAL->Rebin(5); // MG modif
    else
      fhMggEMCAL->Rebin();

    fhMggEMCAL->Fit("fitMass", "L R +", "", 0.05, 0.20);

    fhMggEMCAL->SetTitle("Pi0 Mass in EMCal");
    fhMggEMCAL->SetTitleSize(0.1);
    fhMggEMCAL->SetXTitle("Pi0 Mass");
    fhMggEMCAL->SetYTitle("Nb of entries");
    fhMggEMCAL->GetXaxis()->SetLabelSize(0.03);
    fhMggEMCAL->GetXaxis()->SetTitleSize(0.03);
    fhMggEMCAL->GetXaxis()->SetTitleOffset(1.3);
    fhMggEMCAL->GetYaxis()->SetLabelSize(0.03);
    fhMggEMCAL->GetYaxis()->SetTitleSize(0.03);
    fhMggEMCAL->GetYaxis()->SetTitleOffset(1.3);

    MeanPosEMCAL = fhMggEMCAL->GetFunction("fitMass")->GetParameter(1)*1000;
    MeanPosEMCALErr = fhMggEMCAL->GetFunction("fitMass")->GetParError(1)*1000;
    WidthEMCAL = fhMggEMCAL->GetFunction("fitMass")->GetParameter(2)*1000;
    WidthEMCALErr = fhMggEMCAL->GetFunction("fitMass")->GetParError(2)*1000;
    Chi2NdfPi0EMCAL = fhMggEMCAL->GetFunction("fitMass")->GetChisquare()/fhMggEMCAL->GetFunction("fitMass")->GetNDF();
    Npi0EMCAL = fhMggEMCAL->GetFunction("fitMass")->GetParameter(0)*fhMggEMCAL->GetFunction("fitMass")->GetParameter(2)*TMath::Sqrt(2*TMath::Pi())/(Nevent*fhMggEMCAL->GetBinWidth(10));
    Npi0EMCALErr = TMath::Sqrt((fhMggEMCAL->GetFunction("fitMass")->GetParError(0)/fhMggEMCAL->GetFunction("fitMass")->GetParameter(0))*(fhMggEMCAL->GetFunction("fitMass")->GetParError(0)/fhMggEMCAL->GetFunction("fitMass")->GetParameter(0))+(WidthEMCALErr/WidthEMCAL)*(WidthEMCALErr/WidthEMCAL));
    Npi0EMCALErr = 0.;
    NggEMCAL = fhMggEMCAL->GetFunction("fitMass")->Integral(0.11, 0.16)/(Nevent*fhMggEMCAL->GetBinWidth(10));
    NggEMCALErr = fhMggEMCAL->GetFunction("fitMass")->IntegralError(0.11, 0.16)/(fhMggEMCAL->Integral()*fhMggEMCAL->GetBinWidth(10));
    SignifEMCAL = Npi0EMCAL/NggEMCAL;
    SignifEMCALErr = TMath::Sqrt((Npi0EMCALErr/Npi0EMCAL*(Npi0EMCALErr/Npi0EMCAL)+(NggEMCALErr/NggEMCAL*(NggEMCALErr/NggEMCAL))));
    SignifEMCALErr = SignifEMCAL*SignifEMCALErr;

    TCanvas* c3 = new TCanvas("Pi0InvMassDCAL", "Pi0 Invariant Mass in DCal", 600, 600);
    c3->SetFillColor(0);
    c3->SetBorderSize(0);
    c3->SetFrameBorderMode(0);

    fhIMDCAL = (TH2F *)outputList->FindObject(fCalorimeter+"_hIMDCAL");
    if(fhIMDCAL){
      // Now Pi0 global trending in DCal
      // TCanvas* c3 = new TCanvas("Pi0InvMassDCAL", "Pi0 Invariant Mass in DCal", 600, 600);
      // c3->SetFillColor(0);
      // c3->SetBorderSize(0);
      // c3->SetFrameBorderMode(0);
      // c2->SetOptStat(1); // MG modif

      // fhIMDCAL = (TH2F *)outputList->FindObject(fCalorimeter+"_hIMDCAL");
      // if(!fhIMDCAL) continue;
      fhIMDCAL->Sumw2();

      fhMggDCAL = (TH1F *)fhIMDCAL->ProjectionY("MggDCAL", 2, 20, "") ; // to modify projection range
      if(fhMggDCAL->GetEntries()==0){
	Error(__FUNCTION__, "The Pi0 histogram in DCal is empty !");  tree->Fill();
	ret=-8;
	continue;
      }

      fitMass->SetParameter(0, 4500);
      fitMass->SetParameter(1, mPDG);
      fitMass->SetParameter(2, 0.01);
      fitMass->SetParameter(3, 0.);
      fitMass->SetParameter(4, fhMggDCAL->GetBinContent(fhMggDCAL->FindBin(0.11)));
      fitMass->SetParameter(5, fhMggDCAL->GetBinContent(fhMggDCAL->FindBin(0.20)));

      if(fhMggDCAL->GetEntries()<5000)
	fhMggDCAL->Rebin(5); // MG modif
      else
	fhMggDCAL->Rebin();

      fhMggDCAL->Fit("fitMass", "L R +", "", 0.05, 0.20);

      fhMggDCAL->SetTitle("Pi0 Mass in DCal");
      fhMggDCAL->SetTitleSize(0.1);
      fhMggDCAL->SetXTitle("Pi0 Mass");
      fhMggDCAL->SetYTitle("Nb of entries");
      fhMggDCAL->GetXaxis()->SetLabelSize(0.03);
      fhMggDCAL->GetXaxis()->SetTitleSize(0.03);
      fhMggDCAL->GetXaxis()->SetTitleOffset(1.3);
      fhMggDCAL->GetYaxis()->SetLabelSize(0.03);
      fhMggDCAL->GetYaxis()->SetTitleSize(0.03);
      fhMggDCAL->GetYaxis()->SetTitleOffset(1.3);

      MeanPosDCAL = fhMggDCAL->GetFunction("fitMass")->GetParameter(1)*1000;
      MeanPosDCALErr = fhMggDCAL->GetFunction("fitMass")->GetParError(1)*1000;
      WidthDCAL = fhMggDCAL->GetFunction("fitMass")->GetParameter(2)*1000;
      WidthDCALErr = fhMggDCAL->GetFunction("fitMass")->GetParError(2)*1000;
      Chi2NdfPi0DCAL = fhMggDCAL->GetFunction("fitMass")->GetChisquare()/fhMggDCAL->GetFunction("fitMass")->GetNDF();
      Npi0DCAL = fhMggDCAL->GetFunction("fitMass")->GetParameter(0)*fhMggDCAL->GetFunction("fitMass")->GetParameter(2)*TMath::Sqrt(2*TMath::Pi())/(Nevent*fhMggDCAL->GetBinWidth(10));
      Npi0DCALErr = TMath::Sqrt((fhMggDCAL->GetFunction("fitMass")->GetParError(0)/fhMggDCAL->GetFunction("fitMass")->GetParameter(0))*(fhMggDCAL->GetFunction("fitMass")->GetParError(0)/fhMggDCAL->GetFunction("fitMass")->GetParameter(0))+(WidthDCALErr/WidthDCAL)*(WidthDCALErr/WidthDCAL));
      Npi0DCALErr = 0.;
      NggDCAL = fhMggDCAL->GetFunction("fitMass")->Integral(0.11, 0.16)/(Nevent*fhMggDCAL->GetBinWidth(10));
      NggDCALErr = fhMggDCAL->GetFunction("fitMass")->IntegralError(0.11, 0.16)/(fhMggDCAL->Integral()*fhMggDCAL->GetBinWidth(10));
      SignifDCAL = Npi0DCAL/NggDCAL;
      SignifDCALErr = TMath::Sqrt((Npi0DCALErr/Npi0DCAL*(Npi0DCALErr/Npi0DCAL)+(NggDCALErr/NggDCAL*(NggDCALErr/NggDCAL))));
      SignifDCALErr = SignifDCAL*SignifDCALErr;
    }

    cout<<"******************"<<endl;
    // end of global trending

    // cout<<"********************************** Total number of clusters **********************************"<<endl;
    NTClusters=fhE->IntegralAndError(fhE->GetXaxis()->FindBin(0.), fhE->GetXaxis()->FindBin(200.), NTClustersRMS);
    NCClusters=fhECharged->IntegralAndError(fhECharged->GetXaxis()->FindBin(0.), fhECharged->GetXaxis()->FindBin(50.), NCClustersRMS);

    if(NTClusters!=0)
      NMatchClustersP = NCClusters/NTClusters;
    else
      NMatchClustersP=0;
    // cout<<"   here is searched value "<<NMatchClustersP<<" "<<endl;

    if(NCClusters!=0)
      NMatchClustersPRMS = NMatchClustersP* TMath::Sqrt(TMath::Power(NCClustersRMS/NCClusters, 2)+TMath::Power(NTClustersRMS/NTClusters, 2));
    else
      NMatchClustersPRMS=1;
    
    ClusterMean=fhNClusters->GetMean();
    ClusterRMS=fhNClusters->GetMeanError();
    CellMean=fhNCells->GetMean();
    CellRMS=fhNCells->GetMeanError();

    tree->Fill();

    TString outfilename = QAPATH +  "Pi0InvMassEMCAL" + fTrigger(r) + ".pdf";
    TString outfilename2 = QAPATH + "Pi0InvMassEMCAL" + fTrigger(r) + ".png";
    if(SavePlots==2)
      c2->SaveAs(outfilename);
    if(SavePlots)
      c2->SaveAs(outfilename2);

    outfilename = QAPATH + "Pi0InvMassSM" + fTrigger(r) + ".pdf";
    outfilename2 = QAPATH + "Pi0InvMassSM" + fTrigger(r) + ".png";
    if(SavePlots==2)
      c1->SaveAs(outfilename);
    if(SavePlots)
      c1->SaveAs(outfilename2);
    // cout<<"**********************************  "<<n<<" **********************************"<<endl; // n = number of SMs

    if(fhIMDCAL){
      TString outfilename4 = QAPATH +  "Pi0InvMassDCAL" + fTrigger(r) + ".pdf";
      TString outfilename5= QAPATH + "Pi0InvMassDCAL" + fTrigger(r) + ".png";
      if(SavePlots==2)
	c3->SaveAs(outfilename4);
      if(SavePlots)
	c3->SaveAs(outfilename5);
    }

    fout->cd();
    fout->Cd(Form("%s/%s/%ld/%s/%s", period.Data(), pass.Data(), RunId, "RunLevelQA", fTrigger.Data()));

    // Cell time distribution for each BC and SM
    if(hTimeSM_BC0) c_TimeSM_BC0->Write();
    if(hTimeSM_BC1) c_TimeSM_BC1->Write();
    if(hTimeSM_BC2) c_TimeSM_BC2->Write();
    if(hTimeSM_BC3) c_TimeSM_BC3->Write();

    c2->Write();
    c1->Write();
    if(fhIMDCAL)
      c3->Write();

    delete c1;
    delete c2;
    if(fhIMDCAL)    delete c3;
    if(hTimeSM_BC0) delete c_TimeSM_BC0;
    if(hTimeSM_BC1) delete c_TimeSM_BC1;
    if(hTimeSM_BC2) delete c_TimeSM_BC2;
    if(hTimeSM_BC3) delete c_TimeSM_BC3;
    if(outputList)  outputList->Delete();

    QAData << RunId << "    " << Nevent
	   << "\n";

    QAData.close();
  }

  ftree->cd();
  tree->Write();
  ftree->Close();

  return ret;

}

//-------------------------------------------------------------------------
TH2F* FormatRunHisto(TH2F* aHisto, const char* title, const char* YTitle){

  if(!aHisto){
    Error(__FUNCTION__, Form("The histogram with title \"%s\" was not found!", title));
    return new TH2F();
  }

  aHisto->SetStats(kFALSE);
  aHisto->SetTitle(title);
  aHisto->SetStats(kFALSE);
  aHisto->SetYTitle(YTitle);
  aHisto->GetYaxis()->SetTitleOffset(1.2);
  aHisto->GetYaxis()->SetLabelSize(0.03);
  aHisto->GetZaxis()->SetLabelSize(0.02);

  return aHisto;

}

//--------------------------------------------------------------------------------------------------------------
TH2F* HistoPerMod(TH2F* hTmpPerMod, const char* title){

  if(!hTmpPerMod){
    Error(__FUNCTION__, Form("The histogram with title \"%s\" was not found!", title));
    return new TH2F();
  }

  hTmpPerMod->SetStats(kFALSE);
  hTmpPerMod->SetTitle(title);
  hTmpPerMod->SetTitleSize(0.1);
  hTmpPerMod->GetXaxis()->SetTitleOffset(1.1);
  hTmpPerMod->GetXaxis()->SetTitleSize(0.05);
  hTmpPerMod->GetXaxis()->SetLabelSize(0.06);
  hTmpPerMod->GetYaxis()->SetTitleOffset(1.1);
  hTmpPerMod->GetYaxis()->SetTitleSize(0.05);
  hTmpPerMod->GetYaxis()->SetLabelSize(0.06);
  hTmpPerMod->GetZaxis()->SetLabelSize(0.04);

  return hTmpPerMod;

}

//---------------------------------------------------------------------------------------------------
TH2* AutoZoom(TH2* H, Option_t* aType, Int_t EntryMin){

  Int_t shiftX = (Int_t)(H->GetNbinsX()/30.);
  Int_t shiftY = (Int_t)(H->GetNbinsY()/30.);

  TString opt = aType;
  opt.ToLower();

  int minX = 0;
  int maxX = H->GetNbinsX();
  int New_minX = minX;
  int New_maxX = maxX;

  int minY = 0;
  int maxY = H->GetNbinsY();
  int New_minY = minY;
  int New_maxY = maxY;

  if(opt.Contains("all")) opt = TString("minx, maxx, miny, maxy");

  if(opt.Contains("maxx")){
    for(New_maxX = maxX;New_maxX >=minX; New_maxX--){
      Stat_t c = 0;
      for(int i_y = maxY; i_y >= minY;i_y--){
	c = H->GetBinContent(New_maxX, i_y);
	if(c>EntryMin)
	  break;
      }
      if(c>EntryMin)
	break;
    }
  }

  if(opt.Contains("maxy")){
    for(New_maxY = maxY;New_maxY >=minY;New_maxY--){
      Stat_t c = 0;
      for(int i_x=maxX; i_x>=minX;i_x--){
	c = H->GetBinContent(i_x, New_maxY );
	if(c>EntryMin)
	  break;
      }
      if(c>EntryMin)
	break;
    }
  }

  if(opt.Contains("minx")){
    for(New_minX = minX;New_minX <=maxX; New_minX++){
      Stat_t c = 0;
      for(int i_y = minY; i_y <= maxY;i_y++){
	c = H->GetBinContent(New_minX, i_y);
	if(c>EntryMin)
	  break;
      }
      if(c>EntryMin)
	break;
    }
  }

  if(opt.Contains("miny")){
    for(New_minY = minY;New_minY <=maxY;New_minY++){
      Stat_t c = 0;
      for(int i_x=minX; i_x<=maxX;i_x++){
	c = H->GetBinContent(i_x, New_minY );
	if(c>EntryMin)
	  break;
      }
      if(c>EntryMin)
	break;
    }
  }

  if(New_maxX!=-1 && New_maxY!=-1)
    H->GetXaxis()->SetRange(New_minX - shiftX, New_maxX + shiftX);
  if(New_maxX!=-1 && New_maxY!=-1)
    H->GetYaxis()->SetRange(New_minY - shiftY, New_maxY + shiftY);

  return H;

}

//----------------------------------------------------------------------------------------------------
int FindNumberOfSM(TFile* f, TString fTrigger, TString period){

  TString direct;
  if(!fTrigger.Contains("QA")){
    direct = "CaloQA_";
  }
  direct += fTrigger;

  Int_t nSMt=-1;
  Int_t year = 2000 + TString(period(3, 2)).Atoi();
  if     ( year == 2010 )                { nSMt=4;  }
  else if( year == 2011 || year == 2012 ){ nSMt=10; }
  else if( year == 2013 || year == 2014 ){ nSMt=12; }
  else                                   { nSMt=20; }

  TList *outputList;
  Bool_t dirok = f->cd(direct);

  // ////
  // cout << "==============================================================" << endl;
  // cout << "direct = " << direct << endl;
  // cout << "dirok = " << dirok << endl;
  // cout << "==============================================================" << endl;
  // ////

  if(!dirok)
    Error(__FUNCTION__, Form("No input directory %s, the number SMs will be returned based on the year!", direct.Data()));
  else
    outputList = (TList*)gDirectory->Get(direct);
  
  if(!outputList)
    Error(__FUNCTION__, "No input list, the number SMs will be returned based on the year! ");
  else{
    outputList->SetOwner();
    TH2F* hNSM =(TH2F *)outputList->FindObject("EMCAL_hE_Mod");
    if(!hNSM || (!hNSM->GetEntries()))
      Error(__FUNCTION__, "hNSM Histogram not found or it is empty, the number SMs will be returned based on the year!");
    else
      nSMt = hNSM->GetYaxis()->GetBinUpEdge(hNSM->FindLastBinAbove(0, 2));
  }

  if(outputList)
    outputList->Delete();

  return nSMt;

}
