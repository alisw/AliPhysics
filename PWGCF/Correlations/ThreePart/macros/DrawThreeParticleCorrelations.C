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
  bool kin48 = false;
  bool kin816 = false;
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
    key="48";
    if(argument.CompareTo(key)==0){
      kin48=true;
      continue;
    }
    key="816";
    if(argument.CompareTo(key)==0){
      kin816=true;
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
  
  int remade = 0;
  for(int i = 0;i<folderlist->GetEntries();i++){
    cout << folderlist->At(i)->GetName()<<endl;
    if(TString(folderlist->At(i)->GetName()).BeginsWith("ThreePart")&&((TString(folderlist->At(i)->GetName()).Contains("4_8")&&kin48))||(TString(folderlist->At(i)->GetName()).Contains("8_16")&&kin816)){
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
		signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k));//->Clone(Form("%sclone",containterdir->At(k)->GetName())));
// 		  for(int l =0;l<containterdir->GetEntries();l++){
// 		     if(TString(containterdir->At(l)->GetName()).Contains("trigger_correlation")&&TString(containterdir->At(l)->GetName()).Contains("MEAllfull")){
// 		       AliCorrelation3p * mixed = dynamic_cast<AliCorrelation3p*>(containterdir->At(l));
// 		       mixed->SetMixedEvent(NULL);
// 		       signaltrack->SetMixedEvent(mixed);
// 		    }}
		if(j==0&&i==0&&remade ==0){signaltrack->MakeResultsFile(Form("%s",containername),true);remade =1;}
		else signaltrack->MakeResultsFile(Form("%s",containername));
// 		delete signaltrack;
	      }
// 	      
	      if(TString(containterdir->At(k)->GetName()).Contains("META")&&makefilemixed){
		AliCorrelation3p* signaltrack;
		  signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k));//->Clone(Form("%sclone",containterdir->At(k)->GetName())));
// 		  for(int l =0;l<containterdir->GetEntries();l++){
// 		     if(TString(containterdir->At(l)->GetName()).Contains("trigger_correlation")&&TString(containterdir->At(l)->GetName()).Contains("MEAllfull")){
// 		       signaltrack->SetMixedEvent(dynamic_cast<AliCorrelation3p*>(containterdir->At(l)));}}
		signaltrack->MakeResultsFile(Form("%s/META",containername));
// 		delete signaltrack;
	      }
	      if(TString(containterdir->At(k)->GetName()).Contains("METrigger")&&makefiletriggermixed){
		AliCorrelation3p* signaltrack;
		  signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k));//->Clone(Form("%sclone",containterdir->At(k)->GetName())));
// 		  for(int l =0;l<containterdir->GetEntries();l++){
// 		     if(TString(containterdir->At(l)->GetName()).Contains("trigger_correlation")&&TString(containterdir->At(l)->GetName()).Contains("MEAllfull")){
// 		       signaltrack->SetMixedEvent(dynamic_cast<AliCorrelation3p*>(containterdir->At(l)));}}
		  signaltrack->MakeResultsFile(Form("%s/METrigger",containername));
// 		  delete signaltrack;
	      }	      
	    }
	  }
	}
      }
    }
  }
  ffile->Close();
}