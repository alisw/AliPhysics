
class  AliAnalysisManager;
class  AliAnalysisAlien;
class AliSpectraBothPID;
void runGridBoth(TString mode="test",Int_t mc=0,Bool_t aod=kTRUE,Int_t day=15,Int_t month=6, Int_t year=2012, Bool_t sddin=kFALSE,Int_t pass=4) 
{
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

  //__________________________________________________________________________
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include  ");
  gSystem->SetIncludePath("-I.");
  gROOT->ProcessLine(".include $ALICE_ROOT/TOF ");
gROOT->ProcessLine(".include $ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD");
  //gSystem->Load("libPWGLFspectra.so");
//  gROOT->LoadMacro("AliSpectraBothTrackCuts.cxx+g");
 // gROOT->LoadMacro("AliSpectraBothEventCuts.cxx+g");
  // gROOT->LoadMacro("HistogramNames.cxx");
 // gROOT->LoadMacro("AliSpectraBothHistoManager.cxx+g");
 // gROOT->LoadMacro("AliSpectraBothPID.cxx+g");
 // gROOT->LoadMacro("AliAnalysisTaskSpectraBoth.cxx+g");
  gROOT->LoadMacro("$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD/AddTaskSpectraBoth.C");
	
  // Create and configure the alien handler plugin
  AliAnalysisGrid *alienHandler = CreateAlienHandler(mode,mc,aod,daystring,sddin,pass);  
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
  Double_t y=.5;
  Double_t ptTofMatch=.6;
  UInt_t trkbit=1;
  UInt_t trkbitQVector=1;
  Bool_t UseCentPatchAOD049=kFALSE;
  Double_t DCA=100000;
  UInt_t minNclsTPC=50;
  Int_t rebinfactor=1;	
 Bool_t mcfactor=0;
 if(mc)
	mcfactor=1;	

 gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
cout<<"test MC "<<mc<<endl;
  AliAnalysisTaskPIDResponse *taskPID=AddTaskPIDResponse(mcfactor,kTRUE,kFALSE,2,kFALSE,"",kTRUE);


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
    	AliAnalysisTaskSpectraBoth*  task2=0x0;
    	AliAnalysisTaskSpectraBoth*  task3=0x0;
    	AliAnalysisTaskSpectraBoth*  task4=0x0;

	task1=AddTaskSpectraBoth(mcfactor,-1,-1,-1,-1,-0.8,0.8,Nsigmapid,pt,p,y,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,rebinfactor);
  	task1->SetUseMinSigma(false);
	task1->GetPID()->SetPIDtype(AliSpectraBothPID::kNSigmaTPCTOF);

	if(!aod)
	{	
		task2=AddTaskSpectraBoth(mcfactor,-2,-2,-2,-2,-0.8,0.8,Nsigmapid,pt,p,y,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,rebinfactor);
  		task2->SetUseMinSigma(false);
  		task2->GetPID()->SetPIDtype(AliSpectraBothPID::kNSigmaTOF);	
		task3=AddTaskSpectraBoth(mcfactor,-3,-3,-3,-3,-0.8,0.8,Nsigmapid,pt,p,y,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,rebinfactor);
  		task3->SetUseMinSigma(false);
		task3->GetPID()->SetPIDtype(AliSpectraBothPID::kNSigmaTOF);
	}
	else
	{
	 task2=AddTaskSpectraBoth(mcfactor,-2,-2,-2,-2,-0.8,0.8,Nsigmapid,pt,p,y,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,rebinfactor);
                task2->SetUseMinSigma(true);
                task2->GetPID()->SetPIDtype(AliSpectraBothPID::kNSigmaTPCTOF);
                task3=AddTaskSpectraBoth(mcfactor,-3,-3,-3,-3,-0.8,0.8,2.0,pt,p,y,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,rebinfactor);
                task3->SetUseMinSigma(false);
                task3->GetPID()->SetPIDtype(AliSpectraBothPID::kNSigmaTPCTOF);
	}
	task4=AddTaskSpectraBoth(mcfactor,-4,-4,-4,-4,-0.8,0.8,Nsigmapid,pt,p,y,ptTofMatch,trkbit,trkbitQVector,UseCentPatchAOD049,DCA,minNclsTPC,rebinfactor);
  	task4->SetUseMinSigma(false);
	task4->GetPID()->SetPIDtype(AliSpectraBothPID::kNSigmaTOF);
 
  if(mc)
  {
  	task1->SetdotheMCLoopAfterEventCuts(kFALSE);
	if(!aod)
	{
		task2->SetdotheMCLoopAfterEventCuts(kFALSE);
		task3->SetdotheMCLoopAfterEventCuts(kFALSE);
	}
	task4->SetdotheMCLoopAfterEventCuts(kFALSE);


  }
  if(!aod)
  {
	AliESDtrackCuts* cuttpc1=AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
	task1->SetAliESDtrackCuts(cuttpc1);
	AliESDtrackCuts* cuttpc2=AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(0);
	task2->SetAliESDtrackCuts(cuttpc2);
	AliESDtrackCuts* cuttpc3=AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(0);
	cuttpc3->SetMaxChi2TPCConstrainedGlobal(36);
	task3->SetAliESDtrackCuts(cuttpc3);
	AliESDtrackCuts* cuttpc4=AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();

	task4->SetAliESDtrackCuts(cuttpc4);
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


AliAnalysisGrid* CreateAlienHandler(TString mode="test",Int_t mc=0,Bool_t aod=kTRUE,TString daystring,Bool_t sddin,Int_t pass)
{
  TString filename="*AliAOD.root";
  TString datapattern="ESDs";
  if(sddin)
	datapattern+=Form("/pass%d_with_SDD",pass);
  else
	datapattern+=Form("/pass%d_without_SDD",pass);
		
  if(!aod)
  {
	filename="*AliESDs.root";
	if(sddin&&mc)
		filename="*AliESDs_wSDD.root";
  }
  else
  {
	if(pass==4)
		datapattern+="/AOD113";
	if(pass==2)
		datapattern+="/AOD052";		
  } 
 TString datadirectory="/alice/data/2011/";		 	  
  TString datasettmp="LHC11a";	

cout<<"mc "<<mc<<" aod "<<aod<<endl;
    		
	if(mc==1)
	{
		datadirectory="/alice/sim/2011";
		datasettmp="LHC11e3a_plus";
		if(aod)
			datapattern="AOD079";
		else
			datapattern="";

		
	}
	else if(mc==2)
	{
		datadirectory="/alice/sim/2012";
		datasettmp="LHC12f1a";
		datapattern="";
	}
	else if(mc==3)	
	{
		datadirectory="/alice/sim/2012";
		datasettmp="LHC12f1b";
		datapattern="";

	}
	else if(mc==4)	
	{
		datadirectory="/alice/sim/";
		datasettmp="LHC11b10a";
		if(aod)
			datapattern="AOD116";
		else
			datapattern="";



	}
	

 
  TString jobsettings=Form("dataset%sAOD%dSDD%d",datasettmp.Data(),aod,sddin);	
  if(!mc)
  {
	jobsettings+="pass";
	jobsettings+=pass;

   }			

  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TOF  -I$ALICE_ROOT/PWGLF/SPECTRA/PiKaPr/TestAOD");
  plugin->SetAdditionalLibs("libSTEERBase.so libESD.so libAOD.so libANALYSISalice.so libPWGLFspectra.so libTENDER.so libTENDERSupplies.so");
  //plugin->SetAnalysisSource("AliSpectraBothHistoManager.cxx AliSpectraBothTrackCuts.cxx AliSpectraBothEventCuts.cxx AliSpectraBothPID.cxx AliAnalysisTaskSpectraBoth.cxx");
  
  plugin->SetOverwriteMode();
  plugin->SetExecutableCommand("aliroot -q -b");  
  plugin->SetRunMode(mode.Data());
  plugin->SetNtestFiles(1);
  //Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  plugin->SetROOTVersion("v5-34-02");
  plugin->SetAliROOTVersion("v5-04-09-AN");
  // Declare input data to be processed.
  plugin->SetGridDataDir(Form("%s/%s/",datadirectory.Data(),datasettmp.Data()));
  plugin->SetDataPattern(Form("%s/%s",datapattern.Data(),filename.Data()));
  plugin->SetGridWorkingDir(Form("/pp2.76TeV%s/%s/",daystring.Data(),jobsettings.Data()));
  plugin->SetAnalysisMacro(Form("pp276%s.C",jobsettings.Data()));
  plugin->SetExecutable(Form("pp276%s.sh",jobsettings.Data()));
  plugin->SetJDLName(Form("Taskpp276%s.jdl",jobsettings.Data()));

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
  plugin->SetTTL(10*3600);
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  //plugin->SetSplitMaxInputFileNumber();
  plugin->SetSplitMode("se");
  return plugin;
}
