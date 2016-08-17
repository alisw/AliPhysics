void
CreateIndex(const TString& dir, const TString& tree="esdTree",
	    const char* remote=0)
{
  gROOT->SetMacroPath(Form("$ALICE_PHYSICS/PWGLF/FORWARD/trains:%s", 
			   gROOT->GetMacroPath()));
  gROOT->LoadMacro("ChainBuilder.C+");
  gROOT->Macro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/scripts/LoadLibs.C");
  
  Bool_t mc = false;
  Bool_t zip = false;
  if (tree.BeginsWith("mc")) { 
    mc = true;
    zip = true;
    tree.Remove(0,2);
  }
  if (tree.BeginsWith("zip")) { 
    zip = true;
    tree.Remove(0,3);
  }
  
  TString pat("*.root");
  if      (tree.EqualTo("esdTree",  TString::kIgnoreCase)) pat="AliESDs*";
  else if (tree.EqualTo("aodTree",  TString::kIgnoreCase)) pat="AliAOD*";
  else    Warning("", "Unknown tree: %s, pattern set to *.root", tree.Data());
  if (zip) {
    pat.Prepend("root_archive.zip@");
    pat.ReplaceAll("*", ".root");
  }


  TString opts;
  opts.Append(Form("pattern=%s", pat.Data()));
  opts.Append("&check");
  opts.Append("&clean");
  opts.Append("&recursive");
  opts.Append("&verbose");
  if (mc) opts.Append("&mc");

  TString realDir(dir);
  if (!remote)              realDir = gSystem->ExpandPathName(dir.Data());
  if (realDir.EqualTo(".")) realDir = gSystem->WorkingDirectory();

  TUrl url;
  url.SetProtocol("local");
  url.SetPort(0);
  url.SetFile(realDir);
  url.SetAnchor(tree);
  url.SetOptions(opts);
  
  Printf("Running ChainBuilder::CreateCollection(\"%s/index.root\",\"%s\")",
	 realDir.Data(), url.GetUrl());
  TString out(Form("%s/%s.root", realDir.Data(),
		   !remote ? "index" : "remote"));
  ChainBuilder::CreateCollection(out, url, remote);
}

				 
  
  
