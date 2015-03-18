// #include "AliAnalysisTaskCorrelation3p.h"
#include "TFile.h"
#include "TObject.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TNamed.h"
// #include <ZTrees.h>
//Run with aliroot  -b -l ../run-single-task.C'("compile", "AliAnalysisTaskCorrelation3p.cxx")' DrawThreeParticleCorrelations.C
void DrawThreeParticleCorrelations(const char* options="")
{
  bool bdrawbins=true;
  bool bdrawbinstats=false;
  bool makefile =false;
  bool makefileall =false;
  bool makefilemixed =false;
  bool makefiletriggermixed =false;
  bool makefilegen =false;
  bool ispp = false;
  bool isPbPb = false;
  TString underdirectory("");
  TString delimiter(" ");
  TStringToken token(options, delimiter);
  while (token.NextToken()) {
    const char* key=0;
    TString argument=token;
    if (argument.CompareTo("-h")==0 ||argument.CompareTo("--help")==0) {
      AliInfo(Form("Draw options:"
		    "\n\t  mergefirst   - merges the Multiplicity bins and Vz bins before correcting with mixed events (default merges after mixing)."
		    "\n\t  drawbins     - each multiplicity-Vz bin is drawn alone."
		    "\n\t binstats      - draw bin statistics."
		    "\n\t mixed		- draw mixed."
		    "\n\t triggermixed  - draw trigger mixed."
		    ));
      return;
      }
    key="binstats";
    if (argument.CompareTo(key)==0) {
      bdrawbins=false;
      bdrawbinstats=true;
      underdirectory.Clear();
      underdirectory.Append("/binstats");
      continue;
    }
    key="makefile";
    if (argument.CompareTo(key)==0) {
      bdrawbins=false;
      makefile=true;
      continue;
    }    
    key="makefilemixedall";
    if (argument.CompareTo(key)==0) {
      bdrawbins=false;
      makefileall=true;
      continue;
    }    
    key="makefilemixed";
    if (argument.CompareTo(key)==0) {
      bdrawbins=false;
      makefilemixed=true;
      continue;
    }    
    key="makefiletriggermixed";
    if (argument.CompareTo(key)==0) {
      bdrawbins=false;
      makefiletriggermixed=true;
      continue;
    }    
    key="makefilegen";
    if (argument.CompareTo(key)==0) {
      bdrawbins=false;
      makefilegen=true;
      continue;
    }    
    key="pp";
    if(argument.CompareTo(key)==0){
      ispp=true;
      continue;
    }
    key="PbPb";
    if(argument.CompareTo(key)==0){
      isPbPb=true;
      continue;
    }
  }
  TFile * ffile = TFile::Open("AnalysisResults.root", "READ");
  if(!ffile) return;  

  if(makefilegen){
    TDirectory * dir = ffile->GetDirectory("ThreePartGenerator");
    TList* fList;
    dir->GetObject("ThreePartGeneratorCoutput1",fList);
    AliCorrelation3p* generated;
    if(!(makefilemixed||makefiletriggermixed)) {generated = dynamic_cast<AliCorrelation3p*>fList->FindObject("tracktrigger_correlation_4_8");generated->MakeResultsFile("");}
    if(makefilemixed){generated = dynamic_cast<AliCorrelation3p*>fList->FindObject("tracktrigger_correlation_4_8META");generated->MakeResultsFile("META");}
    if(makefiletriggermixed)  {generated = dynamic_cast<AliCorrelation3p*>fList->FindObject("tracktrigger_correlation_4_8METrigger");generated->MakeResultsFile("METrigger");}
    return;
  }
  
  TList * folderlist = ffile->GetListOfKeys();
  
  
  for(int i = 0;i<folderlist->GetEntries();i++){
    if(TString(folderlist->At(i)->GetName()).BeginsWith("ThreePart")){
      const char* foldername = folderlist->At(i)->GetName();
      TDirectory * folderdir = ffile->GetDirectory(foldername);
      TList * Containerlist = folderdir->GetListOfKeys();
      for(int j=0;j<Containerlist->GetEntries();j++){
	if(TString(Containerlist->At(j)->GetName()).BeginsWith("ThreePart")){
	  const char* containername = Containerlist->At(j)->GetName();
	  cout << containername<<endl;
	  TList * containterdir = dynamic_cast<TList*>(folderdir->Get(containername));
	  for(int k = 0;k<containterdir->GetEntries();k++){
	    if(TString(containterdir->At(k)->GetName()).Contains("trigger_correlation")){
	      if(!TString(containterdir->At(k)->GetName()).Contains("ME")&&makefile){
		AliCorrelation3p* signaltrack;
// 		if(!isPbPb){
// 		  		cout << "2"<<endl;
// 
// 		  signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k));
// 		  		cout << "3"<<endl;
// 
// 		}
// 		if(isPbPb){
		  signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k)->Clone(Form("%sclone",containterdir->At(k)->GetName())));
		  for(int l =0;l<containterdir->GetEntries();l++){
		     if(TString(containterdir->At(l)->GetName()).Contains("trigger_correlation")&&TString(containterdir->At(l)->GetName()).Contains("MEAllfull")){
		       AliCorrelation3p * mixed = dynamic_cast<AliCorrelation3p*>(containterdir->At(l));
		       mixed->SetMixedEvent(NULL);
		       signaltrack->SetMixedEvent(mixed);
		    }}
// 		}		       
		signaltrack->MakeResultsFile(Form("%s",containername));
		delete signaltrack;
	      }
	      
	      if(TString(containterdir->At(k)->GetName()).Contains("META")&&makefilemixed){
		AliCorrelation3p* signaltrack;
// 		if(!isPbPb){signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k));}
// 		if(isPbPb){
		  signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k)->Clone(Form("%sclone",containterdir->At(k)->GetName())));
		  for(int l =0;l<containterdir->GetEntries();l++){
		     if(TString(containterdir->At(l)->GetName()).Contains("trigger_correlation")&&TString(containterdir->At(l)->GetName()).Contains("MEAllfull")){
		       signaltrack->SetMixedEvent(dynamic_cast<AliCorrelation3p*>(containterdir->At(l)));}}
// 	      }
		signaltrack->MakeResultsFile(Form("%s/META",containername));
		delete signaltrack;
	      }
	      
	      if(TString(containterdir->At(k)->GetName()).Contains("METrigger")&&makefiletriggermixed){
		AliCorrelation3p* signaltrack;
// 		if(!isPbPb){signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k));}
// 		if(isPbPb){
		  signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k)->Clone(Form("%sclone",containterdir->At(k)->GetName())));
		  for(int l =0;l<containterdir->GetEntries();l++){
		     if(TString(containterdir->At(l)->GetName()).Contains("trigger_correlation")&&TString(containterdir->At(l)->GetName()).Contains("MEAllfull")){
		       signaltrack->SetMixedEvent(dynamic_cast<AliCorrelation3p*>(containterdir->At(l)));}}
		
// 	      }		
		  signaltrack->MakeResultsFile(Form("%s/METrigger",containername));
		  delete signaltrack;
	      }	      
// 	      if(TString(containterdir->At(k)->GetName()).Contains("MEAllfull")&&makefiletriggermixed){
// 		AliCorrelation3p* signaltrack;
// 		if(!isPbPb){
// 		  signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k));
// 		}
// 		if(isPbPb){continue;}		
// 		  signaltrack->MakeResultsFile(Form("%s/MEALL",containername));
// 	      }
	    }
	  }
	}
      }
      
//       delete Containerlist;
    }
  }
//   delete folderlist;
  ffile->Close();
  /*
  
  TStringToken triggers("Tracks pi0s"," ");
  while(triggers.NextToken()) {
    for (int i =0;i<2;i++){
      TString Directory(Form("ThreePart%s",triggers.Data()));
      if(i == 1) Directory.Append("PbPb");
      TDirectory * dir = ffile->GetDirectory(Directory.Data());
      if(!dir) continue;//directory does not exist.
      TList* fList;
      dir->GetObject(Form("%s%s",Directory.Data(),"Coutput1"),fList);
      if(!fList) continue;//list does not exist.
      if(bdrawbinstats){
	TH1D* nass = (TH1D*)fList->FindObject("NAssociated");
	TH1D* ntrig = (TH1D*)fList->FindObject("Ntriggers");
	TH1D* nasstrig =(TH1D*)fList->FindObject("NAssociatedETriggered");
	TCanvas* c = new TCanvas("numberoftriggers");
	c->Divide(2,2);
	c->cd(1);
	gPad->SetLogy();
	ntrig->GetXaxis()->SetTitle("N_{triggers}");
	ntrig->SetStats(0);
	ntrig->Draw("E");
	double eventwithoutatrigger = ntrig->GetBinContent(1);
	int eventwithoutatriggerint = ntrig->GetBinContent(1);
	double eventswithatrigger = ntrig->Integral(2,ntrig->GetNbinsX());
	int eventswithatriggerint = ntrig->Integral(2,ntrig->GetNbinsX());
	double ratio = eventswithatrigger/(eventwithoutatrigger+eventwithoutatrigger)*100;
	c->cd(2);
	gPad->SetLogy();
	nass->GetXaxis()->SetTitle("N_{associated}");
	nass->SetStats(0);
	nass->Draw("E");
	c->cd(3);
	gPad->SetLogy();
	nasstrig->GetXaxis()->SetTitle("N_{associated}");
	nasstrig->SetTitle("Number of Associated in Events that contain at least one trigger.");
	nasstrig->SetStats(0);
	nasstrig->Draw("E");
	double eventswithatriggerandtwoassociated = nasstrig->Integral(3,nasstrig->GetNbinsX());
	c->cd(4);
	TPaveText* text = new TPaveText(.05,.1,.95,.8,"nb");
	text->SetLineColor(0);
	text->SetFillStyle(0);
	text->AddText(Form("#bullet %4.1f %% of all events contain at least one trigger.",ratio));
	text->AddText(Form("#bullet %.1e events with at least one trigger.",eventswithatrigger));
	text->AddText(Form("#bullet %.1e events with\n at least one trigger and at least two associated.",eventswithatriggerandtwoassociated));
	text->Draw();
	c->Print("Triggersperevent.eps");
	TH3D* eventsbeforeselection = (TH3D*)fList->FindObject("Eventbeforeselection");
	eventsbeforeselection->GetYaxis()->SetTitle("Multiplicity");
	eventsbeforeselection->GetXaxis()->SetTitle("ZVertex");
	eventsbeforeselection->GetZaxis()->SetTitle("Centrality");
	eventsbeforeselection->SetTitleSize(0.045,"xyz");
	eventsbeforeselection->SetTitleOffset(1.3,"xyz");
	TH3D* eventsafterselection  = (TH3D*)fList->FindObject("Eventafterselection");
	eventsafterselection->GetYaxis()->SetTitle("Multiplicity");
	eventsafterselection->GetXaxis()->SetTitle("ZVertex");
	eventsafterselection->GetZaxis()->SetTitle("Centrality");
	eventsafterselection->SetTitleSize(0.045,"xyz");
	eventsafterselection->SetTitleOffset(1.3,"xyz");
	TH2D* multvertexbeforeselection = (TH2D*) eventsbeforeselection->Project3D("yx");
	multvertexbeforeselection->SetTitle("multiplicity vs z vertex before event selection.");
	TH2D* centvertexbeforeselection = (TH2D*) eventsbeforeselection->Project3D("zx");
	centvertexbeforeselection->SetTitle("centrality vs z vertex before event selection.");
	TH1D* multbeforeselection = (TH1D*) eventsbeforeselection->Project3D("x");
	multbeforeselection->SetTitle("multiplicity before event selection");
	TH1D* vertexbeforeselection = (TH1D*)  eventsbeforeselection->Project3D("y");
	vertexbeforeselection->SetTitle("vertex before event selection");
	TH1D* centralitybeforeselection = (TH1D*)  eventsbeforeselection->Project3D("z");
	centralitybeforeselection->SetTitle("centrality before event selection");
	TH2D* multvertexafterselection = (TH2D*) eventsafterselection->Project3D("yx");
	multvertexafterselection->SetTitle("multiplicity vs z vertex after event selection.");
	TH2D* centvertexafterselection = (TH2D*) eventsafterselection->Project3D("zx");
	centvertexafterselection->SetTitle("centrality vs z vertex after event selection.");
	TH1D* multafterselection = (TH1D*) eventsafterselection->Project3D("y");
	multafterselection->SetTitle("multiplicity after event selection");
	if(i==0)multafterselection->GetXaxis()->SetRange(1,multafterselection->GetXaxis()->FindBin(200));
	TH1D* vertexafterselection = (TH1D*)  eventsafterselection->Project3D("x");
	vertexafterselection->SetTitle("vertex after event selection");
	TH1D* centralityafterselection = (TH1D*)  eventsafterselection->Project3D("z");
	centralityafterselection->SetTitle("centrality after event selection");
	c->Clear();
	c->Divide(2,2);
	c->cd(1);
	if(i==0)multvertexbeforeselection->GetYaxis()->SetRange(1,multvertexafterselection->GetYaxis()->FindBin(200));
	multvertexbeforeselection->Draw("surf3");
	multvertexbeforeselection->SetStats(0);
	gPad->SetTheta(30); // default is 30
	gPad->SetPhi(150); // default is 30
	gPad->Update();
	c->cd(2);
	if(i==0)multvertexafterselection->GetYaxis()->SetRange(1,multvertexafterselection->GetYaxis()->FindBin(200));
	multvertexafterselection->Draw("surf3");
	multvertexafterselection->SetStats(0);
	gPad->SetTheta(30); // default is 30
	gPad->SetPhi(150); // default is 30
	gPad->Update();
	c->cd(3);
	multafterselection->SetStats(0);
	multafterselection->SetTitleOffset(1.1,"xy");
	multafterselection->Draw("E");
	gPad->SetLogy();
	c->cd(4);
	vertexafterselection->SetStats(0);
	vertexafterselection->SetTitleOffset(1.1,"xy");
	vertexafterselection->Draw("E");
	c->Print("MultVertex.eps");
	if(i==1){
	  c->Clear();
	  c->Divide(2,2);
	  c->cd(1);
	  centvertexbeforeselection->Draw("surf3");
	  centvertexbeforeselection->SetStats(0);
	  centvertexbeforeselection->GetYaxis()->SetRange(1,centvertexbeforeselection->GetNbinsX());

	  gPad->SetTheta(30); // default is 30
	  gPad->SetPhi(150); // default is 30
	  gPad->Update();
	  c->cd(2);
	  centvertexafterselection->Draw("surf3");
	  centvertexafterselection->SetStats(0);
	  gPad->SetTheta(30); // default is 30
	  gPad->SetPhi(150); // default is 30	  
	  gPad->Update();
	  c->cd(3);
	  centralityafterselection->SetStats(0);
	  centralityafterselection->SetTitleOffset(1.1,"xy");
	  centralityafterselection->Draw("E");
// 	  gPad->SetLogy();
	  c->cd(4);
	  vertexafterselection->SetStats(0);
	  vertexafterselection->SetTitleOffset(1.1,"xy");
	  vertexafterselection->Draw("E");
	  c->Print("CentVertex.eps");
	}
	delete c;
      }
      //Extract the appropriate Correlation object
      TString corobj(Form("%strigger_correlation_4_8",triggers.Data()));
      TString corobj2(Form("%strigger_correlation_8_16",triggers.Data()));

      if(triggers.CompareTo("Tracks")==0){
	corobj.ReplaceAll("Tracks","track");
	corobj2.ReplaceAll("Tracks","track");

	AliCorrelation3p* signaltrack = dynamic_cast<AliCorrelation3p*>(fList->FindObject(corobj.Data()));
	TString save = TString(corobj.Data());
	corobj.Append("META");
	AliCorrelation3p* mixedtrack = dynamic_cast<AliCorrelation3p*>(fList->FindObject(corobj.Data()));
	corobj = TString(save.Data());
	corobj.Append("METrigger");	
	AliCorrelation3p* mixedtriggertrack = dynamic_cast<AliCorrelation3p*>(fList->FindObject(corobj.Data()));
	corobj = TString(save.Data());
	corobj.Append("MEAllfull");	
	AliCorrelation3p* mixedtrackall = dynamic_cast<AliCorrelation3p*>(fList->FindObject(corobj.Data()));
	if(!signaltrack){
	  signaltrack = dynamic_cast<AliCorrelation3p*>(fList->FindObject(corobj2.Data()));
	  if(!signaltrack) continue;
	  mixedtrack= dynamic_cast<AliCorrelation3p*>(fList->FindObject(Form("%s%s",corobj2.Data(),"META")));
	  mixedtrackall= dynamic_cast<AliCorrelation3p*>(fList->FindObject(Form("%s%s",corobj2.Data(),"MEAllfull")));

	  mixedtriggertrack= dynamic_cast<AliCorrelation3p*>(fList->FindObject(Form("%s%s",corobj2.Data(),"METrigger")));
	}

	if(makefile&&0) signaltrack->MakeResultsFile("");

	if(makefile&&1) signaltrack->MakeResultsFile("");
	if(makefilemixed&&0) mixedtrack->MakeResultsFile("META");
	if(makefileall) mixedtrackall->MakeResultsFile("MEALL");
	if(makefilemixed&&1) mixedtrack->MakeResultsFile("META");
	if(makefiletriggermixed&&0) mixedtriggertrack->MakeResultsFile("METrigger");
	if(makefiletriggermixed&&1) mixedtriggertrack->MakeResultsFile("METrigger");	
      }
    }
  }*/
}