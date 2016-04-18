#ifndef __CINT__
# include <TObject.h>
# include <TSystem.h>
# include <TROOT.h>
# include <TFile.h>
void SimpleCorrect(UShort_t,UShort_t,const char*,const char*,Int_t,const char*);
#else
class TCanvas;
#endif
const Bool_t kCorrectLoaded = true;
void
AddPath(const TString& dir, Bool_t prepend=true)
{
  TString d(gSystem->ExpandPathName(dir.Data()));
  gSystem->AddIncludePath("-I%s", d.Data());
  const char* oldPath = gROOT->GetMacroPath();
  gROOT->SetMacroPath(Form(".:%s:%s",
			   prepend ? d.Data() : oldPath,
			   prepend ? oldPath  : d.Data()));
}

void
Correct(UShort_t    flags=0x3,
	const char* side="middle",
	const char* var="none",
	Bool_t      forceK=false)
{
  const char* fwd = "$ALICE_ROOT/PWGLF/FORWARD/analysis2";
  AddPath(TString::Format("%s/dndeta/tracklets", fwd));
  if (!gROOT->GetListOfGlobals()->FindObject("kSimpleCorrectLoaded"))
    gROOT->LoadMacro("SimpleCorrect.C");
  SimpleCorrect(flags, TString(var).EqualTo("none") || forceK ? 2 : 3,
		Form("dt_%s_%s/trdt.root", side, "none"),
		Form("mc_%s_%s/trmc.root", side, var),
		9,
		"");
  TString resFile;
  resFile.Form("result_0x%x.root", flags & 0x3);
  TObject* resObj = gROOT->GetListOfFiles()->FindObject(resFile);
  if (resObj) {
    TFile* tmp = static_cast<TFile*>(resObj);
    Printf("Closing %s", tmp->GetName());
    tmp->Close();
  }
  gSystem->mkdir("partial", 1);
  TString dest;
  dest.Form("partial/%s_%s_0x%x.root", side, var, flags&0x3);
  Printf("%s -> %s", resFile.Data(), dest.Data());
  gSystem->Rename(resFile, dest);
}


		
