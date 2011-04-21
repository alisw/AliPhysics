void TriggerStudyResults(const char * filename = "outTrigger/collection_136854.xml/trigger_study.root", TString trigger = "C0SM1-A-NOPF-ALL") {

  TFile * f = new TFile (filename);

  
  LoadLibs();

  // Draw trigger venn diagram
  //    hTrigStat_All_C0OM2-A-NOPF-ALL
  TString vennName = "hTrigStat_All";
  if (trigger != "") vennName = vennName+"_"+trigger;
  cout << "Venn name: " << vennName.Data() << endl;
  
  TH1 * hVenn = (TH1*) f->Get(vennName);
  hVenn->Draw();
  Int_t colors[] = {2,3,4,5,6,7,8,9,10,11,12,13,14,15};
  TPie * pie = new TPie(hVenn);
  Int_t nbin = hVenn->GetNbinsX();
  for(Int_t ibin = 0; ibin < nbin; ibin++){
    pie->SetEntryLabel(ibin, Form("%s, %d (%1.2f%%)", pie->GetEntryLabel(ibin), (int) pie->GetEntryVal(ibin), pie->GetEntryVal(ibin)/hVenn->GetEntries()*100));
  }
  
  pie->SortSlices(1,1); // Sort slices and merge the empty ones
  pie->SortSlices(0,0); // Sort in increasing order. Sorting in increasing order directly leads to segfault
  pie->SetRadius(.31);
  pie->SetX(.53);
  pie->SetY(.34);
  pie->SetLabelsOffset(.01);
  pie->SetLabelFormat("");//#splitline{%val (%perc)}{%txt}");
  pie->SetFillColors(colors);
  //  pie->SetTextSize(0.01);
  pie->Draw("");
  
  pie->MakeLegend(0.224832, 0.66958, 0.833893, 0.97028, trigger.Data());
  cout << pie << endl;
  
  // Make a table of trigger efficiencies for all histos results
  AliLatexTable table(2,"cc");
  //  table.InsertCustomRow("\\multicolumn{c}{2}{Integrated efficiency}");
  table.InsertCustomRow(Form("Trigger Name & Efficiency (%s)\\\\",trigger.Data()));
  table.InsertHline();
  TList * l = gDirectory->GetListOfKeys();
  TIterator * iter = l->MakeIterator();
  TKey * key = 0;
  TH1F * hall = (TH1F*) gDirectory->Get("hTracklets_all"); // FIXME: get the normalization for a given trigger?
  while (key = (TKey*) iter->Next()){
    TString name = key->GetName();
    if(!name.Contains("Tracklets")) continue;
    if(!name.Contains(trigger)) continue;
    if(name.Contains("all")) continue;
    TH1F * h = (TH1F*) gDirectory->Get(name);
    TString label = name(name.Index("_")+1, name.Index("_",name.Index("_")+1)-name.Index("_")-1);
    table.SetNextCol(label);
    table.SetNextCol(h->GetEntries()/hall->GetEntries());
    table.InsertRow();
  }
  cout << "Integrated trigger efficiency" << endl;
  table.PrintTable("ASCII");
}

LoadLibs() {

  gSystem->Load("libCore.so");  
  gSystem->Load("libGeom.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libProof");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libMinuit");
  gSystem->Load("libPWG2Spectra");
  gSystem->Load("libPWG0base"); 
   
  gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG0/multPbPb"));
  gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG1/background"));
  // Load helper classes
  // TODO: replace this by a list of TOBJStrings
  TString taskName("$ALICE_ROOT/PWG0/multPbPb/AliAnalysisTaskMultPbTracks.cxx+");
  TString histoManName("$ALICE_ROOT/PWG0/multPbPb/AliAnalysisMultPbTrackHistoManager.cxx+");
  TString centrName("$ALICE_ROOT/PWG0/multPbPb/AliAnalysisMultPbCentralitySelector.cxx+");
  TString listName("$ALICE_ROOT/PWG1/background/AliHistoListWrapper.cxx+");

  gSystem->ExpandPathName(taskName);
  gSystem->ExpandPathName(histoManName);
  gSystem->ExpandPathName(centrName);
  gSystem->ExpandPathName(listName);

  Bool_t debug=0;
  gROOT->LoadMacro(listName    +(debug?"+g":""));   
  gROOT->LoadMacro(histoManName+(debug?"+g":""));
  gROOT->LoadMacro(centrName   +(debug?"+g":""));   
  gROOT->LoadMacro(taskName    +(debug?"+g":""));   

  // Histo fitter
  // gROOT->LoadMacro("/Users/mfloris/Work/ALICE/ANALYSIS/HistoFitter/fcn.cxx+g");
  // gROOT->LoadMacro("/Users/mfloris/Work/ALICE/ANALYSIS/HistoFitter/AliHistoFitter.cxx+g");


}
