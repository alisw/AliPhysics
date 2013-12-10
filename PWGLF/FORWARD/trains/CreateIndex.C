void
CreateIndex(const TString& dir, const TString& tree="esdTree")
{
  gROOT->SetMacroPath(Form("$ALICE_ROOT/PWGLF/FORWARD/trains:%s", 
			   gROOT->GetMacroPath()));
  gROOT->LoadMacro("ChainBuilder.C+");
  
  Bool_t mc = false;
  if (tree.BeginsWith("mc")) { 
    mc = true;
    tree.Remove(0,2);
  }
  
  TString pat("*.root");
  if      (tree.EqualTo("esdTree", TString::kIgnoreCase)) pat="AliESDs*";
  else if (tree.EqualTo("aodTree", TString::kIgnoreCase)) pat="AliAOD*";
  else    Warning("", "Unknown tree: %s, pattern set to *.root", tree.Data());
  if (mc) pat.Prepend("root_archive.zip@");


  TString opts;
  opts.Append(Form("pattern=%s", pat.Data()));
  opts.Append("&check");
  opts.Append("&clean");
  opts.Append("&recursive");
  if (mc) opts.Append("&mc");

  TUrl url;
  url.SetProtocol("local");
  url.SetPort(0);
  url.SetFile(dir);
  url.SetAnchor(tree);
  url.SetOptions(opts);
  
  Printf("Running ChainBuilder::CreateCollection(\"%s/index.root\",\"%s\")",
	 dir.Data(), url.GetUrl());
  ChainBuilder::CreateCollection(Form("%s/index.root", dir.Data()), url);
}

				 
  
  
