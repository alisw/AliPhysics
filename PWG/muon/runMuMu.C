///
/// Example macro to run the AliAnalysisTaskMuMu task.
/// That's just an example. Feel free to use your own version of it...
///
/// Typical usage is :
///
/// > root
/// root[] .x runMuMu.C(dataset,issimulation,"username@aaf.domain")
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
/// Before running it, please add/modify the triggers list in the GetTriggerList function,
/// and, if running with your own version of some packages, the packages vector
/// in the SetupLibraries function.
///
/// @author L. Aphecetche
///

//______________________________________________________________________________
Bool_t SetupLibraries(Bool_t local, Bool_t debug, const char* mainPackage)
{
  std::vector<std::string> packages;
  
 // packages.push_back("PWGmuon"); // uncomment this if you want to use your version of that package
  
  if (!local)
  {
    TList list;
    
    list.Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));
    list.Add(new TNamed("ALIROOT_MODE","base"));
    
    TString extraLibs;
    
    if ( !packages.empty() )
    {
      for ( std::vector<std::string>::size_type i = 0; i < packages.size(); ++i )
      {
        extraLibs += packages[i];
        extraLibs += ":";
      }
    }
    
    extraLibs.Remove(TString::kBoth,':');
    
    list.Add(new TNamed("ALIROOT_EXTRA_LIBS",extraLibs.Data()));
    
    if (gProof->EnablePackage(mainPackage,&list))
    {
      std::cout << "Enable of main package failed" << std::endl;
      return kFALSE;
    }
  }
  
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
        return kFALSE;
      }
    }
    else
    {
      if ( debug )
      {
        std::cout << "Uploading PAR file " << package << std::endl;
      }
      
      if (gProof->UploadPackage(package.c_str()))
      {
        std::cout << "Upload failed" << std::endl;
        return kFALSE;
      }
      
      if ( debug )
      {
        std::cout << "Enabling PAR file " << package << std::endl;
      }
      
      Int_t rv = gProof->EnablePackage(package.c_str(),(TList*)(0x0),kFALSE);
      
      if (rv)
      {
        std::cout << "Enable failed" << std::endl;
        return kFALSE;
      }
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
    outputname.ReplaceAll("*","");
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
      if (master.Contains("nansafmaster2")) af = "saf2";
      if (master.Contains("nansafmaster3")) af = "saf3";
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
void GetTriggerList(Bool_t simulations, TString& triggers, TString& inputs)
{
  // Get the list of triggers to be analyzed
  //
  // Here you can use class names and/or combinations of classes and input, e.g.
  // triggers->Add(new TObjString("CINT7-B-NOPF-ALLNOTRD & 0MUL"));
  
  if (!simulations)
  {
    triggers =  "CMUL7-B-NOPF-MUFAST,"
    "CINT7-B-NOPF-MUFAST,"
    "CMSL7-B-NOPF-MUFAST,"
    "CMUL7-B-NOPF-MUFAST&0MUL,"
    "CMSL7-B-NOPF-MUFAST&0MUL,"
    "CINT7-B-NOPF-MUFAST&0MUL";
    inputs = "0MSL:17,0MSH:18,0MLL:19,0MUL:20";
  }
  else
  {
    // add here the MC-specific trigger you want to analyze (if any)

    triggers = "CMULLO-B-NOPF-MUON,"
    "CMSNGL-B-NOPF-MUON,"
    "ANY";    
    // e.g. for dpmjet simulations (at least) we have the following "triggers" :
    // C0T0A,C0T0C,MB1,MBBG1,V0L,V0R,MULow,EMPTY,MBBG3,MULL,MULU,MUHigh
    inputs = "";
  }
}

//______________________________________________________________________________
AliAnalysisTask* runMuMu(const char* dataset="Find;BasePath=/alice/data/2015/LHC15i/000235839/muon_calo_pass1/AOD/;FileName=AliAOD.Muons.root;Tree=/aodTree",
                         Bool_t simulations=kFALSE,
                         const char* addtask="AddTaskMuMuMinv",
                         const char* where="laphecet@nansafmaster2.in2p3.fr/?N")
{
  ///
  /// @param dataset the name of either a dataset or a text file containing dataset names
  /// @param simulations set it to true if the dataset contains MC information (and you want to analyze it)
  /// @param where the connection string for an AF or 0 for local analysis
  ///
  
  // below a few parameters that should be changed less often
  TString workers;
  
  if ( !TString(where).Contains("pod") )
  {
    workers = "workers=8x";
  }
  
  TString saf2Package("VO_ALICE@AliPhysics::vAN-20150908"); // only care about this if you run on SAF2
  
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
  
  TString mainPackage = saf2Package;
  
  if ( !local )
  {
    p = TProof::Open(where,workers.Data());
    if (!p)
    {
      std::cout << "Cannot connect to Proof : " << where << std::endl;
      return 0x0;
    }
    
    TString master(gProof->GetSessionTag());
    if (master.Contains("nansafmaster3") )
    {
      // dealing with a VAF, main package is the same regardless of the aliphysics version you asked
      mainPackage = "/opt/SAF3/etc/vaf/AliceVaf.par";
    }
    if (master.Contains("alivaf") )
    {
      mainPackage = "/afs/cern.ch/alice/offline/vaf/AliceVaf.par";
    }
  }
  
  if (!SetupLibraries(local,debug,mainPackage.Data()))
  {
    std::cout << "Cannot setup libraries. Aborting" << std::endl;
    return 0x0;
  }
  
  AliAnalysisManager* mgr = new AliAnalysisManager("MuMu");
  
  AliInputEventHandler* input = GetInput(sds,p);
  
  if (!input)
  {
    std::cout << "Cannot get input type !" << std::endl;
    return 0x0;
  }
  
  mgr->SetInputEventHandler(input);
  
  TString triggers, inputs;

  GetTriggerList(simulations,triggers,inputs);
  
  TString outputname = addtask;

  outputname.ReplaceAll("AddTask","");
  
  AliAnalysisTask* task(0x0);
  
  if (!baseline)
  {
    gROOT->LoadMacro(Form("%s.C",addtask));
    task = static_cast<AliAnalysisTaskMuMu*>(gROOT->ProcessLine(Form(".x %s.C(\"%s\",\"%s\",\"%s\",\"%s\",%d)",addtask,outputname.Data(),triggers.Data(),inputs.Data(),"pp",simulations)));
  }
  else
  {
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddTaskBaseLine.C");
    task = AddTaskBaseLine();
  }
  
  if (!mgr->InitAnalysis())
  {
    std::cout << "Could not InitAnalysis" << std::endl;
    return 0x0;
  }
  
  mgr->PrintStatus();
  task->Print();
  
  //  return task;
  
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
      std::cout << "Cannot create input chain !" << std::endl;
      return 0x0;
    }
    if (debug) mgr->SetNSysInfo(10);
    //    if( debug) mgr->SetDebugLevel(10);
    TStopwatch timer;
    mgr->StartAnalysis("local",c);
    timer.Print();
    if (debug)
    {
      mgr->ProfileTask("AliAnalysisTaskMuMu");
      if (baseline) mgr->ProfileTask("baseline");
    }
  }
  
  AliCodeTimer::Instance()->Print();
  
  return task;
}

