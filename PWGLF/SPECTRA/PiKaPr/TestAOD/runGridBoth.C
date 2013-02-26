
class  AliAnalysisManager;
class  AliAnalysisAlien;
class AliSpectraBothPID;
Int_t mc=-1;
Int_t aod=-1;
TString griddatadir="";
TString datapattern="";
TString gridworkingdir="";
TString aliroot="";
TString root="";


void runGridBoth(TString mode="test",TString configfile="config.txt",Int_t day=15,Int_t month=6, Int_t year=2012) 
{


   	ifstream infile(configfile);
	char buffer[256];
	while (infile.eof()==false)
	{
		buffer[0]='#'; 
		while (buffer[0]=='#'&&infile.eof()==false)
			infile.getline(buffer,256);
		TString tmpstring(buffer);
		cout<<buffer<<endl;
		if(tmpstring.Contains("MC"))
		{
			tmpstring.Remove(0,3);
			cout<<tmpstring<<endl;
		        if (!tmpstring.IsDigit())
			{
				cout<<"Wrong format of MC flag in config file "<<endl;
				return; 	
			}
			mc=tmpstring.Atoi();
		}
		else if(tmpstring.Contains("AOD")&&tmpstring.Length()<6)
		{
			tmpstring.Remove(0,4);
		        if (!tmpstring.IsDigit())
			{
				cout<<"Wrong format of AOD flag in config file "<<endl;
				return; 	
			}
			aod=tmpstring.Atoi();
		}
		else if(tmpstring.Contains("GridDataDir"))
		{
			griddatadir=tmpstring.Remove(0,12);		
		}
		else if(tmpstring.Contains("DataPattern"))
		{	
			datapattern=tmpstring.Remove(0,12);
		}
		else if (tmpstring.Contains("GridWorkingDir"))
		{
			gridworkingdir=tmpstring.Remove(0,15);
		}
		else if(tmpstring.Contains("root")&&!tmpstring.Contains("aliroot"))
		{	
			root=tmpstring.Remove(0,5);
		}
		else if (tmpstring.Contains("aliroot"))
		{
			aliroot=tmpstring.Remove(0,8);
		}

		else
			continue;

	}
	 
	if(mc<0||aod<0||griddatadir.Length()<1||root.Length()<1||datapattern.Length()<1||gridworkingdir.Length()<1||aliroot.Length()<1)
	{
		cout<<"lack of config info"<<endl;
 		return;	
	}	
cout<<mc<<" "<<aod<<" "<<griddatadir.Data()<<" "<<root.Data()<<" "<<datapattern.Data()<<" "<<gridworkingdir.Data()<<" "<<aliroot.Data()<<endl;
	

  TString daystring=Form("%d%d%d",year,month,day);	
  //to be used with Aliroot > v5-03-32-AN
  AliLog::SetGlobalDebugLevel(100);
  // Load common libraries
  gEnv->SetValue("XSec.GSI.DelegProxy", "2");
  
  // Load common libraries
  gSystem->Load("libCore.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so"); 
  gSystem->Load("libGui.so");
  gSystem->Load("libXMLParser.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libCDB.so");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libProof.so");
  gSystem->Load("libRAWDatabase.so");
  gSystem->Load("libSTEER.so");
  gSystem->Load("libTENDER.so");
    gSystem->Load("libTENDERSupplies.so");
gSystem->Load("libPWGLFspectra.so");
  gROOT->ProcessLine(".include $ALICE_ROOT/include ");
  gSystem->SetIncludePath("-I.");
  gROOT->ProcessLine(".include $ALICE_ROOT/TOF ");

  //gSystem->Load("libPWGLFspectra.so");
 // gROOT->LoadMacro("AliSpectraBothTrackCuts.cxx+g");
 // gROOT->LoadMacro("AliSpectraBothEventCuts.cxx+g");
  // gROOT->LoadMacro("HistogramNames.cxx");
 // gROOT->LoadMacro("AliSpectraBothHistoManager.cxx+g");
 // gROOT->LoadMacro("AliSpectraBothPID.cxx+g");
 // gROOT->LoadMacro("AliAnalysisTaskSpectraBoth.cxx+g");
  //__________________________________________________________________________
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include  ");
  gSystem->SetIncludePath("-I.");
  gROOT->ProcessLine(".include $ALICE_ROOT/TOF ");
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD/AddTaskSpectraBoth.C");
	
  // Create and configure the alien handler plugin
  AliAnalysisGrid *alienHandler = CreateAlienHandler(mode,daystring);  
  if (!alienHandler) return;
  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);
  if(aod)
  {	
  	AliAODInputHandler* aodH = new AliAODInputHandler();
  	mgr->SetInputEventHandler(aodH);
  }
  else
  {
	 AliESDInputHandler* esdH = new AliESDInputHandler();	
	 mgr->SetInputEventHandler(esdH);
	 if(mc)
        {
             AliMCEventHandler *mch = new AliMCEventHandler();
             mgr->SetMCtruthEventHandler(mch);
	}
  }		
  
  // Add PID task
 // Bool_t isMC=kFALSE, Bool_t autoMCesd=kTRUE,
   //    Bool_t tuneOnData=kFALSE, Int_t recoPass=2, Bool_t cachePID=kFALSE, TString detResponse="",
//Bool_t useTPCEtaCorrection = kFALSE  gROOT->LoadMacro("./AddTaskSpectraBoth.C");

  Double_t Nsigmapid=3.;
  Double_t pt=5.;
  Double_t p=5.;
  Double_t y=.2;
  Double_t ptTofMatch=.6;
  UInt_t trkbit=1;
  UInt_t trkbit2=16;

  UInt_t trkbitQVector=1;
  Bool_t UseCentPatchAOD049=kFALSE;
  Double_t DCA=100000;
  UInt_t minNclsTPC=70;
  UInt_t minNclsTPC2=0;	
  Int_t rebinfactor=1;		
 Bool_t mcfactor=0;
 if(mc)
	mcfactor=1;	

 gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse(mcfactor,kTRUE,kFALSE,2,kFALSE,"",kTRUE);
taskPID->SetUseTPCEtaCorrection(kTRUE);
 if(!aod)
  {	
  	gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
        AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection();
        if(mc)
        {
             physSelTask->GetPhysicsSelection()->SetAnalyzeMC();
         } 
  }
	

    	AliAnalysisTaskSpectraBoth* task1=0x0;  
	AliAnalysisTaskSpectraBoth* task2=0x0;
	AliAnalysisTaskSpectraBoth* task3=0x0;
	AliAnalysisTaskSpectraBoth* task4=0x0;  



	task1=AddTaskSpectraBoth(mcfactor,-1,-1,-1,-1,-0.8,0.8,Nsigmapid,pt,p,-0.5,0.5,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,rebinfactor,"V0M",2);

	task2=AddTaskSpectraBoth(mcfactor,-1,-1,-1,-1,-0.8,0.8,Nsigmapid,pt,p,-0.465,0.035,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,rebinfactor,"V0M",2);
	
	task3=AddTaskSpectraBoth(mcfactor,-1,-1,-1,-1,-0.8,0.8,2.0,pt,p,-0.465,0.035,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,rebinfactor,"V0M",2);

	task4=AddTaskSpectraBoth(mcfactor,-1,-1,-1,-1,-0.8,0.8,Nsigmapid,pt,p,-0.465,0.035,ptTofMatch,trkbit2,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC2,rebinfactor,"V0M",2);

	


  if(mc)
  {
  	task1->SetdotheMCLoopAfterEventCuts(kTRUE);
	task2->SetdotheMCLoopAfterEventCuts(kTRUE);
	task3->SetdotheMCLoopAfterEventCuts(kTRUE);
	task4->SetdotheMCLoopAfterEventCuts(kTRUE);


  }
  else
 {
	task1->GetEventCuts()->SetTriggerSettings(AliVEvent::kINT7);
	task2->GetEventCuts()->SetTriggerSettings(AliVEvent::kINT7);
	task3->GetEventCuts()->SetTriggerSettings(AliVEvent::kINT7);
	task4->GetEventCuts()->SetTriggerSettings(AliVEvent::kINT7);


}	
  if(!aod)
  {
	AliESDtrackCuts* cuttpc1=AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
	task1->SetAliESDtrackCuts(cuttpc1);
	AliESDtrackCuts* cuttpc2=AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

	task2->SetAliESDtrackCuts(cuttpc2);
  }
	
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  AddTaskPIDqa(); 
   if(mode.Contains("test"))
	mgr->SetDebugLevel(1);
  //mgr->Init();
  if (!mgr->InitAnalysis())return;
  mgr->PrintStatus();

  mgr->StartAnalysis("grid");
}


AliAnalysisGrid* CreateAlienHandler(TString mode="test",TString daystring)
{
  TString filename="*AliAOD.root";
 
		
  if(!aod)
  {
	filename="*AliESDs.root";
  }



 
  TString jobsettings=griddatadir;
  jobsettings+=datapattern;
  if(aod)
	jobsettings+="AOD";
  else
	jobsettings+="ESD";

  if(mc)
	jobsettings+="MC";
  else
	jobsettings+="DATA";
		
		
  jobsettings.ReplaceAll("/","");	

 jobsettings.ReplaceAll("*","");
	
	

  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  //plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TOF  -I$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD");
  //plugin->SetAdditionalLibs("libSTEERBase.so libESD.so libAOD.so libANALYSISalice.so libPWGLFspectra.so libTENDER.so libTENDERSupplies.so");
  //plugin->SetAnalysisSource("AliSpectraBothHistoManager.cxx AliSpectraBothTrackCuts.cxx AliSpectraBothEventCuts.cxx AliSpectraBothPID.cxx AliAnalysisTaskSpectraBoth.cxx");
    plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TOF ");
  plugin->SetAdditionalLibs("libTENDER.so libTENDERSupplies.so libPWGLFspectra.so");
 // plugin->SetAnalysisSource("AliSpectraBothHistoManager.cxx AliSpectraBothTrackCuts.cxx AliSpectraBothEventCuts.cxx AliSpectraBothPID.cxx AliAnalysisTaskSpectraBoth.cxx");
  plugin->SetOverwriteMode();
  plugin->SetExecutableCommand("aliroot -q -b");  
  plugin->SetRunMode(mode.Data());
  plugin->SetNtestFiles(1);
  //Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion(root.Data());
  plugin->SetAliROOTVersion(aliroot.Data());
  // Declare input data to be processed
  plugin->SetGridDataDir(griddatadir.Data());
  plugin->SetDataPattern(Form("%s/%s",datapattern.Data(),filename.Data()));
  plugin->SetGridWorkingDir(Form("%s%s",gridworkingdir.Data(),daystring.Data()));
  plugin->SetAnalysisMacro(Form("%s.C",jobsettings.Data()));
  plugin->SetExecutable(Form("%s.sh",jobsettings.Data()));
  plugin->SetJDLName(Form("Task%s.jdl",jobsettings.Data()));

  if(mc)
  {
      plugin->SetRunPrefix(""); 
    }  
  else
    {
      plugin->SetRunPrefix("000"); 
    }   
  FILE* listruns=fopen("runs.txt","r");
  Int_t irun;
  while(!feof(listruns))
  {
      fscanf(listruns,"%d\n",&irun);
      plugin->AddRunNumber(irun);
  }
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output 
  plugin->SetMaxMergeFiles(25);
  plugin->SetMergeExcludes("EventStat_temp.root"
                            "event_stat.root");
  plugin->SetMergeViaJDL(true);
  plugin->SetTTL(15*3600);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  //plugin->SetSplitMaxInputFileNumber();
  plugin->SetSplitMode("se");
  return plugin;
}
