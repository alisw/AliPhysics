Bool_t
GetInfoFromFileName(const char* path, 
		    Int_t& firstRun, Int_t& lastRun, 
		    Int_t& version, Int_t& subVersion)
{
  TString fn(gSystem->BaseName(path));
  int ret = sscanf(fn.Data(), "Run%d_%d_v%d_s%d.root", 
		   &firstRun, &lastRun, &version, &subVersion);
  return ret == 4;
}

void 
FixOne(const char* path, const char* id, const char* fn) 
{
  TFile* file = TFile::Open(Form("%s/%s/%s", path, id, fn));
  if (!file) { 
    Warning("FixOne", "%s/%s/%s not found", path, id, fn);
    return;
  }

  AliCDBEntry* entry = static_cast<AliCDBEntry*>(file->Get("AliCDBEntry"));
  if (!entry) { 
    Warning("FixOne", "Did not find an entry in the file");
    return;
  }

  Int_t firstRun, lastRun, version, subVersion;
  if (!GetInfoFromFileName(fn, firstRun, lastRun, version, subVersion)) { 
    Warning("FixOne", "Could not get info from file name");
    return;
  }

  if (firstRun    == entry->GetId().GetFirstRun() && 
      lastRun     == entry->GetId().GetLastRun()  && 
      version     == entry->GetId().GetVersion()  && 
      subVersion  == entry->GetId().GetSubVersion()) { 
    Info("FixOne", "Entry run range and version %d->%d %d.%d match "
	 "filename run range and version %d->%d %d.%d, not creating %s in %s", 
	 entry->GetId().GetFirstRun(), entry->GetId().GetLastRun(), 
	 entry->GetId().GetVersion(),  entry->GetId().GetSubVersion(), 
	 firstRun, lastRun, version, subVersion, fn, 
	 gSystem->WorkingDirectory());
    return;
  }

  AliCDBId cid(id, firstRun, lastRun, version, subVersion);
  entry->SetId(cid);

  Info("FixOne", "Creating %s", fn);
  TFile* out = TFile::Open(fn, "RECREATE");
  out->cd();
  entry->Write();
  out->ls();
  out->Close();
  file->Close();
}

FixCategory(const char* path, const char* id)
{
  TSystemDirectory dir(Form("%s/%s", path, id), Form("%s/%s", path, id));
  TList*       fl = dir.GetListOfFiles();
  TIter        next(fl);
  TSystemFile* file = 0;
  while ((file = static_cast<TSystemFile*>(next()))) { 
    if (file->IsDirectory()) { 
      Info("FixCategory", " skipping %s", file->GetName());
      continue;
    }

    TString fn(gSystem->BaseName(file->GetName()));
    if (!fn.EndsWith(".root")) { 
      Info("FixCategory", " Skipping %s", fn.Data());
      continue; 
    }

    Info("FixCategory", " Now processing %s/%s @ %s", id, fn.Data(), path);
    FixOne(path, id, fn.Data());
  }
}

FixOCDBEntries()
{
  const char* path = gSystem->ExpandPathName("$ALICE_ROOT/OCDB");

  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(Form("local://%s", path));
  cdb->SetRun(0);

  AliFMDParameters* fmdp = AliFMDParameters::Instance();
  fmdp->Init(kTRUE);
  
  const char* ids[] = { AliFMDParameters::PulseGainPath(), 
			AliFMDParameters::PedestalPath(),
			AliFMDParameters::DeadPath(),
			AliFMDParameters::ZeroSuppressionPath(),
			AliFMDParameters::SampleRatePath(),
			AliFMDParameters::StripRangePath(),    
			AliFMDParameters::AltroMapPath(),
			0 };
  const char** ptr = ids;
  while (*ptr) { 
    Info("FixOCDBEntries", "Now processing %s", *ptr);
    FixCategory(path, *ptr);
    ptr++;
  }
}


    
