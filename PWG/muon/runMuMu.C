///
/// Example macro to run the AliAnalysisTaskMuMu task
///
/// \author L. Aphecetche
///

//______________________________________________________________________________
void LoadLocalLibs(Bool_t localAnalysis=kTRUE)
{
  gSystem->Load("libVMC");
  gSystem->Load("libMinuit");
  gSystem->Load("libTree");
  gSystem->Load("libProofPlayer");
  gSystem->Load("libXMLParser");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");

  if (!localAnalysis)
  {
    gSystem->Load("libCORRFW");
  }
  else
  {
  	gROOT->LoadMacro("AliOADBMuonTrackCutsParam.cxx+g");
  	gROOT->LoadMacro("AliAnalysisMuonUtility.cxx+g");
  	gROOT->LoadMacro("AliMuonTrackCuts.cxx+g");
  	gROOT->LoadMacro("AliMergeableCollection.cxx+g");
  	gROOT->LoadMacro("AliAnalysisMuMuBinning.cxx+g");
  	gROOT->LoadMacro("AliMuonEventCuts.cxx+g");
    
    gROOT->LoadMacro("AliAnalysisMuMuCutElement.cxx+g");
    gROOT->LoadMacro("AliAnalysisMuMuCutCombination.cxx+g");
    gROOT->LoadMacro("AliAnalysisMuMuCutRegistry.cxx+g");
    gROOT->LoadMacro("AliAnalysisMuMuEventCutter.cxx+g");
    gROOT->LoadMacro("AliAnalysisMuMuBase.cxx+g");

    gROOT->LoadMacro("AliAnalysisTaskMuMu.cxx+g");

    gROOT->LoadMacro("AliAnalysisMuMuGlobal.cxx+g");

    gROOT->LoadMacro("AliAnalysisMuMuMinv.cxx+g");
    gROOT->LoadMacro("AliAnalysisMuMuSingle.cxx+g");
    gROOT->LoadMacro("AliAnalysisMuMuNch.cxx+g");

  }
}

//______________________________________________________________________________
TChain* CreateLocalChain(const char* filelist)
{
	TChain* c = new TChain("aodTree");
	
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
   if (sds.Length()==0 ) return "AOD";

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
AliAnalysisTask* runMuMu(const char* dataset="SIM_JPSI_LHC13f_CynthiaTuneWithRejectList_000197388",
                         Bool_t simulations=kTRUE,
                         Bool_t baseline=kFALSE,
                         const char* where="laphecet@nansafmaster.in2p3.fr/?N")
{
  // Create the analysis manager
  
  Bool_t prooflite = (strlen(where)==0) || TString(where).Contains("workers");

  TString sds(dataset);
  
//  if (!prooflite && sds.Length()>0) TProof::Mgr(where)->SetROOTVersion("VO_ALICE@ROOT::v5-34-05");
  
  TProof* p(0x0);
  TString alirootMode("");
  TString workers("workers=8x");

  if (TString(where).Contains("alice-caf"))
  {
    workers="workers=1x";
  }
  if (TString(where).Contains("localhost:2093"))
  {
    workers="workers=8x";
  }
  
  if (prooflite)
  {
    cout << "Will work in LITE mode" << endl;
  }
  
  if ( sds.Length()>0 )
  {
    p = TProof::Open(where,workers.Data());
    
    if (!p)
    {
      cout << "Cannot connect to Proof : " << where << endl;
      return 0;
    }
    
    alirootMode.ToUpper();
    
    if ( alirootMode == "PAR" ) 
    {
      cout << "Will work with PAR files" << endl;
      
      std::vector<std::string> pars;
      
      pars.push_back("STEERBase");
      pars.push_back("ESD");
      pars.push_back("AOD");
      pars.push_back("ANALYSIS");
      pars.push_back("OADB");
      pars.push_back("ANALYSISalice");
      //       pars.push_back("CORRFW");
      //       pars.push_back("PWGmuon");
      
      Bool_t ok(kTRUE);
      
      for ( std::vector<std::string>::size_type i = 0; i < pars.size(); ++i )
      {
        std::string package = pars[i];
        
        if ( gProof->UploadPackage(package.c_str()) ) 
        {
          ok = kFALSE;
        }
        
        if ( gProof->EnablePackage(package.c_str(),"",kTRUE) ) 
        {
          ok = kFALSE;
        }
        
        if (!ok) 
        {
          cout << "Problem with PAR " << package.c_str() << endl;
          return 0;           
        }
      }
    }
    else 
    {
      TList* list = new TList();
      
      //       list->Add(new TNamed("ALIROOT_EXTRA_LIBS", "PWG3base"));
      //       list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES", "PWG3"));
      //       list->Add(new TNamed("ALIROOT_EXTRA_LIBS", "PWG3base"));//:CORRFW:PWG3muon"));
      //       list->Add(new TNamed("ALIROOT_EXTRA_INCLUDES", "PWG3/base"));//:PWG3/muon"));
      
      //    list->Add(new TNamed("ALIROOT_ENABLE_ALIEN", "1"));
      
      if (!alirootMode.IsNull())
      {
        list->Add(new TNamed("ALIROOT_MODE", alirootMode.Data()));  
      }
      else
      {
        list->Add(new TNamed("ALIROOT_MODE",""));
      }
      
      if (!prooflite)
      {
//        p->SetParameter("PROOF_UseTreeCache", 0);
        p->EnablePackage("VO_ALICE@AliRoot::v5-04-65-AN", list, kTRUE);
      }
      else
      {
        //      list->Add(new TNamed("ALIROOT_LOCAL_PATH",gSystem->Getenv("ALICE_ROOT")));       
        p->UploadPackage("$ALICE_ROOT/ANALYSIS/macros/AliRootProofLite.par");
        if (p->EnablePackage("AliRootProofLite",list)) return 0;
      }
    }
    
    // compile task on workers
    if ( alirootMode != "PAR" )
    {
      p->Load("AliOADBMuonTrackCutsParam.cxx+");
      p->Load("AliAnalysisMuonUtility.cxx+");
      p->Load("AliMuonTrackCuts.cxx+");
      p->Load("AliMergeableCollection.cxx+");
      p->Load("AliAnalysisMuMuBinning.cxx+");
      p->Load("AliMuonEventCuts.cxx+");
       p->Load("AliAnalysisMuMuCutElement.cxx+");
       p->Load("AliAnalysisMuMuCutCombination.cxx+");       
       p->Load("AliAnalysisMuMuCutRegistry.cxx+");
      p->Load("AliAnalysisMuMuBase.cxx+");      
      p->Load("AliAnalysisTaskMuMu.cxx+");
       p->Load("AliAnalysisMuMuEventCutter.cxx+");       
     p->Load("AliAnalysisMuMuGlobal.cxx+");
     p->Load("AliAnalysisMuMuNch.cxx+");
     p->Load("AliAnalysisMuMuSingle.cxx+");
     p->Load("AliAnalysisMuMuMinv.cxx+");
	}
  }
  
  LoadLocalLibs(kTRUE);
  
  AliAnalysisManager *mgr = new AliAnalysisManager("MuMu");
  
  AliInputEventHandler* input(0x0);

  TString inputType = GetInputType(sds,p);
  
  if ( inputType == "AOD" ) 
  {
    input = new AliAODInputHandler;
  }
  else if ( inputType == "ESD" ) 
  {
    input = new AliESDInputHandler;
  }
  else
  {
  	std::cout << "Cannot get input type !" << std::endl;
  	return 0;
  }

  mgr->SetInputEventHandler(input);
  
  TList* triggers = new TList;
  triggers->SetOwner(kTRUE);

  if (!simulations)
  {
    triggers->Add(new TObjString("CINT7-B-NOPF-ALLNOTRD"));

//   triggers->Add(new TObjString("CINT7-B-NOPF-ALLNOTRD & 0MUL"));
//   triggers->Add(new TObjString("CINT7-B-NOPF-ALLNOTRD & 0MSL"));
//   triggers->Add(new TObjString("CINT7-B-NOPF-ALLNOTRD & 0MSH"));
//   triggers->Add(new TObjString("CMSL7-B-NOPF-MUON & 0MUL"));
//   triggers->Add(new TObjString("CMSL7-B-NOPF-MUON & 0MSH"));
// 
//   triggers->Add(new TObjString("CMSL7-B-NOPF-MUON"));
//   triggers->Add(new TObjString("CMSH7-B-NOPF-MUON"));
//    triggers->Add(new TObjString("CMUL7-B-NOPF-MUON"));
 
  // below for MB periods only
//  triggers->Add(new TObjString("CMSL7-B-NOPF-ALLNOTRD"));
//   triggers->Add(new TObjString("CMSH7-B-NOPF-ALLNOTRD"));
  triggers->Add(new TObjString("CMUL7-B-NOPF-ALLNOTRD"));
  triggers->Add(new TObjString("CMUL7-B-NOPF-MUON"));
//   triggers->Add(new TObjString("CMSL7-B-NOPF-ALLNOTRD & 0MUL"));
//   triggers->Add(new TObjString("CMSL7-B-NOPF-ALLNOTRD & 0MSH"));
  }

  TString outputname("test.MuMu.AOD.1.root");
  
  if ( sds.Length()>0 ) 
  {
  	TString af("local");
  	
  	if ( gProof )
  	{
  	  af="unknown";
  	  TString master(gProof->GetSessionTag());
      if (master.Contains("lx")) af = "caf";
      if (master.Contains("nansaf")) af = "saf";
      if (master.Contains("skaf")) af = "skaf";
      if (master.Contains("localhost:2093")) af="laf";
  	}
  	outputname = Form("%s.%s.root",gSystem->BaseName(sds.Data()),af.Data());
    outputname.ReplaceAll("|","-");
  	cout << outputname << endl;
  }

  AliAnalysisTask* task(0x0);

  if (!baseline)
  {  
  	gROOT->LoadMacro("AddTaskMuMu.C");

  	task = AddTaskMuMu(outputname.Data(),triggers,"pA",simulations);
  }
  else
  {
  	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddTaskBaseLine.C");
  	task = AddTaskBaseLine();
  }
  
  if (!mgr->InitAnalysis()) 
  {
  	cout << "Could not InitAnalysis" << endl;
    return 0;
  }
  
  if ( sds.Length()>0 )
  {
    TStopwatch timer;
    
    mgr->StartAnalysis("proof",sds.Data());
    
    timer.Print();
  }
  else
  {
    mgr->PrintStatus();
  
   	task->Print();
	
//	return task;
	
    TChain* c = CreateLocalChain("list.aod.txt");
//   	mgr->SetNSysInfo(10);
    TStopwatch timer;
//    mgr->SetDebugLevel(10);
    mgr->StartAnalysis("local",c);
    timer.Print();
//    mgr->ProfileTask("AliAnalysisTaskMuMu");
//    if (baseline) mgr->ProfileTask("baseline");
  }
  
  AliCodeTimer::Instance()->Print();
  
  if (alirootMode=="PAR")
  {
    TProofLog *pl = TProof::Mgr(where)->GetSessionLogs(); pl->Save("*","aod.log");
  }
  
  delete triggers;
  
  return task;
}

