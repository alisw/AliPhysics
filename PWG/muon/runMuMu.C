///
/// Example macro to run the AliAnalysisTaskMuMu task
///
/// Typical usage is :
///
/// > root
/// root[] .x runMuMu.C(dataset,issimulation,"username@aaf.domain")
///
/// or
///
/// root[] .x runMuMu.C(dataset,issimulation,"pod://") on a VAF
///
/// Note that you must take care of the default Root version your packages use.
/// If you need to change it you must use (from a separate session) :
///
/// TProof::Mgr(where)->SetROOTVersion("VO_ALICE@ROOT::[version]");
///
/// the list of available versions can be found with :
///
/// TProof::Mgr(where)->ShowROOTVersions();
///
/// You have to do this only once, as this Root version will be
/// remembered by the AF until you change it
///
/// This example assumes that you are working with AliPhysics and 
/// a custom PAR file for the PWGmuon library. If you just want to use
/// PWGmuon from the AliPhysics version you chose, just comment out
/// line 45 in SetupLibraries method.
///
/// You also have to change (or at least x-check) the list of triggers 
/// used (see GetTriggerList method below)
///
/// @author L. Aphecetche
///

//______________________________________________________________________________
Bool_t SetupLibraries(Bool_t local, Bool_t debug)
{
  std::vector<std::string> packages;
  
  if (!local)
  {
    packages.push_back("VO_ALICE@AliPhysics::v5-06-39-01");
  }
  packages.push_back("PWGmuon"); // comment this line out to use the AliPhysics version of PWGmuon
  
  Bool_t ok(kTRUE);
  
  for ( std::vector<std::string>::size_type i = 0; i < packages.size(); ++i )
  {
    const std::string& package = packages[i];
    
    if (local)
    {
      if ( debug )
      {
        std::cout << "Loading local library lib" << package << std::endl;
      }
      if (gSystem->Load(Form("lib%s",package.c_str()))<0)
      {
        ok = kFALSE;
      }
    }
    else
    {
      if ( debug )
      {
        std::cout << "Uploading/enabling PAR file " << package << std::endl;
      }

      if ( gProof->UploadPackage(Form("$ALICE_PHYSICS/PARfiles/%s",package.c_str())) )
      {
        ok = kFALSE;
      }
      
      if ( gProof->EnablePackage(package.c_str(),"",kFALSE) )
      {
        ok = kFALSE;
      }
    }
    
    if (!ok)
    {
      std::cerr << "Problem with package " << package.c_str() << std::endl;
      return kFALSE;
    }
  }
  
  return kTRUE;
}

//______________________________________________________________________________
TChain* CreateLocalChain(const char* filelist, const char* treeName)
{
	TChain* c = new TChain(treeName);
	
	char line[1024];
	
	ifstream in(filelist);
	while ( in.getline(line,1024,'\n') )
	{
		c->Add(line);
	}
	return c;
}

//______________________________________________________________________________
TString GetInputType(const TString& sds, TProof* p)
{
  //
  // Get the input type (AOD or ESD) from the dataset type
  // Feel free to update this method to fits your needs !
  //
  
   if (sds.Length()==0 ) 
  {
     if ( gSystem->AccessPathName("list.esd.txt") == kFALSE )
     {
       return "ESD";
     }
     else if ( gSystem->AccessPathName("list.aod.txt") == kFALSE )
     {
       return "AOD";
     }
     else 
     {
        std::cout << "Cannot work in local mode without list.esd.txt or list.aod.txt file... Aborting" << std::endl;
        exit(1);
     }
  }

   if (sds.Contains("SIM_JPSI")) return "AOD";
   
   if (sds.Contains("AOD")) return "AOD";
   if (sds.Contains("ESD")) return "ESD";
   
   if ( gSystem->AccessPathName(gSystem->ExpandPathName(sds.Data())) )
   {
   	// dataset is not a local file so it must be a dataset name
   	if (!p) return "NOPROOF";
   	
   	TFileCollection* fc = p->GetDataSet(sds.Data());
    if (!fc) return "NODATASET";
    
    TIter next(fc->GetList());
    TFileInfo* fi;
    while ( ( fi = static_cast<TFileInfo*>(next()) ) )
    {
      TUrl url(*(fi->GetFirstUrl()));
	  TString surl(url.GetUrl());
	  
	  if (surl.Contains("AOD")) return "AOD";
	  if (surl.Contains("AliESD")) return "ESD";      
    }      

   }
   else
   {
   std::cout << "Will use datasets from file " << sds.Data() << std::endl;
   
   // dataset is a local text file containing a list of dataset names
   	std::ifstream in(sds.Data());
   	char line[1014];
   
   	while (in.getline(line,1023,'\n'))
   	{
    	  TString sline(line);
      	sline.ToUpper();
    	if ( sline.Contains("SIM_JPSI") ) return "AOD";
      	if ( sline.Contains("AOD") ) return "AOD";
      	if ( sline.Contains("ESD") ) return "ESD";
   	}
   }
   
   return "BUG";   
}

//______________________________________________________________________________
AliInputEventHandler* GetInput(const TString& sds, TProof* p)
{
  // Get the input event handler depending on the type of dataset
  
  TString inputType = GetInputType(sds,p);
  
  AliInputEventHandler* input(0x0);
  
  if ( inputType == "AOD" )
  {
    input = new AliAODInputHandler;
  }
  else if ( inputType == "ESD" )
  {
    input = new AliESDInputHandler;
  }
  return input;
}

//______________________________________________________________________________
TString GetOutputName(const TString& sds)
{
  TString outputname = "test.MuMu.AOD.1.root";
  
  if ( sds.Length()>0 && sds.BeginsWith("Find") )
  {
    outputname = sds;
    outputname.ReplaceAll(";Tree=/aodTree","");
    outputname.ReplaceAll(";Tree=/esdTree","");
    outputname.ReplaceAll("Find;","");
    outputname.ReplaceAll("FileName=","");
    outputname.ReplaceAll("BasePath=","");
    outputname.ReplaceAll("/alice","alice");
    outputname.ReplaceAll("/","_");
    outputname.ReplaceAll(";","_");
    outputname.ReplaceAll("*.*","");
    outputname.ReplaceAll("Mode=cache","");
    outputname.ReplaceAll("Mode=remote","");
  }
  else if ( sds.Length() > 0 )
  {
    TString af("local");
    
    if ( gProof )
    {
      af="unknown";
      TString master(gProof->GetSessionTag());
      if (master.Contains("nansaf")) af = "saf";
      if (master.Contains("nansafmaster")) af = "saf3";
      if (master.Contains("skaf")) af = "skaf";
    }
    outputname = Form("%s.%s.root",gSystem->BaseName(sds.Data()),af.Data());
    outputname.ReplaceAll("|","-");
  }
  else if ( gSystem->AccessPathName("list.esd.txt") == kFALSE )
  {
       outputname.ReplaceAll("AOD","ESD");
  }

  cout << "outputname will be " << outputname << endl;
  return outputname;
}

//______________________________________________________________________________
void GetTriggerList(TList& triggers, Bool_t simulations)
{
  // Fill the list of triggers to be analyzed
  //
  // Here you can use class names and/or combinations of classes and input, e.g.
  // triggers->Add(new TObjString("CINT7-B-NOPF-ALLNOTRD & 0MUL"));

  if (!simulations)
  {
    triggers.Add(new TObjString("CMUL7-B-NOPF-MUFAST"));
    triggers.Add(new TObjString("C0MUL-B-NOPF-MUFAST"));
    triggers.Add(new TObjString("CINT7-B-NOPF-MUFAST"));
  }
  else
  {
    // add here the MC-specific trigger you want to analyze (if any)

    triggers.Add(new TObjString("CMULLO-B-NOPF-MUON"));
    triggers.Add(new TObjString("CMSNGL-B-NOPF-MUON"));

    triggers.Add(new TObjString("ANY"));
    
    // e.g. for dpmjet simulations (at least) we have the following "triggers" :
    // C0T0A,C0T0C,MB1,MBBG1,V0L,V0R,MULow,EMPTY,MBBG3,MULL,MULU,MUHigh
  }
}

//______________________________________________________________________________
AliAnalysisTask* runMuMu(const char* dataset="ds.txt",
                         Bool_t simulations=kFALSE,
                         const char* addtask="AddTaskMuMuMinv.C",
                         const char* where="laphecet@nansafmaster2.in2p3.fr/?N")
{
  ///
  /// @param dataset the name of either a dataset or a text file containing dataset names
  /// @param simulations set it to true if the dataset contains MC information (and you want to analyze it)
  /// @param where the connection string for an AF or 0 for local analysis
  ///
  
  // below a few parameters that should be changed less often
  TString workers("workers=8x");
  TString sds(dataset);
  Bool_t baseline(kFALSE); // set to kTRUE in order to debug the AF and/or package list
  // without your analysis but with a baseline one (from the lego train)
  Bool_t debug(kFALSE);
  
  Bool_t prooflite = (strlen(where)==0) || TString(where).Contains("workers");
  Bool_t local = (sds.Length()==0);

  TProof* p(0x0);

  if (prooflite)
  {
    std::cout << "************* Will work in LITE mode *************" << std::endl;
  }
  
  if ( local )
  {
    std::cout << "************* Will work in LOCAL mode *************" << std::endl;
  }
  
  if ( !local )
  {
    p = TProof::Open(where,workers.Data());
    if (!p)
    {
      std::cerr << "Cannot connect to Proof : " << where << std::endl;
      return 0x0;
    }
  }
   
	if (!SetupLibraries(local,debug))
  {
    std::cerr << "Cannot setup libraries. Aborting" << std::endl;
    return 0x0;
  }
  
  AliAnalysisManager* mgr = new AliAnalysisManager("MuMu");
  
  AliInputEventHandler* input = GetInput(sds,p);

  if (!input)
  {
  	std::cerr << "Cannot get input type !" << std::endl;
  	return 0x0;
  }

  mgr->SetInputEventHandler(input);
  
  TList* triggers = new TList;
  triggers->SetOwner(kTRUE);
  
  GetTriggerList(*triggers,simulations);
  
  TString outputname = GetOutputName(sds);
  
  AliAnalysisTask* task(0x0);

  if (!baseline)
  {  
  	gROOT->LoadMacro(Form("$ALICE_PHYSICS/PWG/muon/%s",addtask));
  	task = AddTaskMuMu(outputname.Data(),triggers,"pp",simulations);
  }
  else
  {
  	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddTaskBaseLine.C");
  	task = AddTaskBaseLine();
  }
  
  if (!mgr->InitAnalysis()) 
  {
    std::cerr << "Could not InitAnalysis" << std::endl;
    return 0x0;
  }
  
  mgr->PrintStatus();
  task->Print();

  if ( !local )
  {
    mgr->StartAnalysis("proof",sds.Data());
  }
  else
  {
    TChain* c = 0x0;
    if ( gSystem->AccessPathName("list.esd.txt") == kFALSE )
    {
      c = CreateLocalChain("list.esd.txt","esdTree");      
    }
    else if ( gSystem->AccessPathName("list.aod.txt") == kFALSE )
    {
      c = CreateLocalChain("list.aod.txt","aodTree");
    }
    if (!c)
    {
      std::cerr << "Cannot create input chain !" << std::endl;
      return 0x0;
    }
    if (debug) mgr->SetNSysInfo(10);
    if( debug) mgr->SetDebugLevel(10);
    mgr->StartAnalysis("local",c);
    if (debug)
    {
      mgr->ProfileTask("AliAnalysisTaskMuMu");
      if (baseline) mgr->ProfileTask("baseline");
      AliCodeTimer::Instance()->Print();
    }
  }
    
  delete triggers;
  
  return task;
}

