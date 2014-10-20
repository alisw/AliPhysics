void
CreateIndex(const TString& dir, const TString& tree="esdTree")
{
  gROOT->SetMacroPath(Form("$ALICE_ROOT/PWGLF/FORWARD/trains:%s", 
			   gROOT->GetMacroPath()));
  gROOT->LoadMacro("ChainBuilder.C+");
  gROOT->Macro("$ALICE_ROOT/PWGLF/FORWARD/analysis2/scripts/LoadLibs.C");
  
  Bool_t mc = false;
  if (tree.BeginsWith("mc")) { 
    mc = true;
    tree.Remove(0,2);
  }
  
  TString pat("*.root");
  if      (tree.EqualTo("esdTree", TString::kIgnoreCase)) pat="AliESDs*";
  else if (tree.EqualTo("aodTree", TString::kIgnoreCase)) pat="AliAOD*";
  else    Warning("", "Unknown tree: %s, pattern set to *.root", tree.Data());
  if (mc) {
    pat.Prepend("root_archive.zip@");
    pat.ReplaceAll("*", ".root");
  }


  TString opts;
  opts.Append(Form("pattern=%s", pat.Data()));
  opts.Append("&check");
  opts.Append("&clean");
  opts.Append("&recursive");
  // opts.Append("&verbose");
  if (mc) opts.Append("&mc");

  TString realDir(gSystem->ExpandPathName(dir.Data()));
  if (realDir.EqualTo(".")) realDir = gSystem->WorkingDirectory();

  TUrl url;
  url.SetProtocol("local");
  url.SetPort(0);
  url.SetFile(realDir);
  url.SetAnchor(tree);
  url.SetOptions(opts);
  
  Printf("Running ChainBuilder::CreateCollection(\"%s/index.root\",\"%s\")",
	 realDir.Data(), url.GetUrl());
  ChainBuilder::CreateCollection(Form("%s/index.root", realDir.Data()), url);
}

				 
  
  
