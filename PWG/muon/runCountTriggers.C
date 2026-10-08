

//______________________________________________________________________________
void AddTask(const char* outputname)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTask", "No analysis manager to connect to.");
    return NULL;
  }
  
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }
  
  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  
  AliAnalysisCountTriggers* task = new AliAnalysisCountTriggers;
  // add the configured task to the analysis manager
  mgr->AddTask(task);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  
  AliAnalysisDataContainer *coutputCC =
  mgr->CreateContainer("CC",AliCounterCollection::Class(),AliAnalysisManager::kOutputContainer,outputname);

  AliAnalysisDataContainer *coutputTM =
  mgr->CreateContainer("TM",TH1I::Class(),AliAnalysisManager::kOutputContainer,outputname);

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 1, coutputCC);
  mgr->ConnectOutput(task, 2, coutputTM);
  
}

//______________________________________________________________________________
Bool_t SetupLibraries(Bool_t local, Bool_t debug)
{
  std::vector<std::string> packages;
  

  packages.push_back("PWGmuon");
  
  Bool_t ok(kTRUE);
  
  for ( std::vector<std::string>::size_type i = 0; i < packages.size(); ++i )
  {
    const std::string& package = packages[i];
    
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
TString GetInputType(const TString& sds)
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
    return "NOINPUT";
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
AliInputEventHandler* GetInput(const TString& sds)
{
  // Get the input event handler depending on the type of dataset
  
  TString inputType = GetInputType(sds);
  
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
  TString outputname = "test.AOD.root";
  
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
  }
  else if ( sds.Length() > 0 )
  {
    TString af("local");
    

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
void runCountTriggers(const char* dataset="",
                                  const char* /*legacyLocation*/="")
{
  // below a few parameters that should be changed less often
  TString sds(dataset);
  Bool_t baseline(kFALSE); // set to kTRUE in order to debug the AF and/or package list
  // without your analysis but with a baseline one (from the lego train)
  Bool_t debug(kFALSE);
  
  Bool_t local = (sds.Length()==0);


  

  {

    std::cout << "************* Will work in LOCAL mode *************" << std::endl;
  
}


   
	if (!SetupLibraries(local,debug))
  {
    std::cerr << "Cannot setup libraries. Aborting" << std::endl;
    return 0x0;
  }
  
  AliAnalysisManager* mgr = new AliAnalysisManager("MuMu");
  
  AliInputEventHandler* input = GetInput(sds);

  if (!input)
  {
  	std::cerr << "Cannot get input type !" << std::endl;
  	return 0x0;
  }

  mgr->SetInputEventHandler(input);
  
  AddTask(GetOutputName(sds));
  
  if (!mgr->InitAnalysis())
  {
    std::cerr << "Could not InitAnalysis" << std::endl;
    return 0x0;
  }
  
  mgr->PrintStatus();

  {
    TChain* c = 0x0;
    if (!sds.IsNull())
    {
      c = CreateLocalChain(sds.Data(), TString(input->GetDataType()) == "AOD" ? "aodTree" : "esdTree");
    }
    else if ( gSystem->AccessPathName("list.esd.txt") == kFALSE )
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
    mgr->StartAnalysis("local",c);
  }
    
}

