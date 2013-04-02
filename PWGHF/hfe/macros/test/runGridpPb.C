#include <fstream>
#include <iostream>
#include <string>

TString g_aliroot_version;
TString g_root_version;
TString g_sample;
TString g_plugin_mode;
TString g_train_dir;
TArrayI g_runlist;

AliAnalysisAlien *CreateGridHandler(){
        //
        // Setup main settings of the Alien plugin
        //
        AliAnalysisAlien *plugin = new AliAnalysisAlien();
        plugin->SetRunMode(g_plugin_mode.Data());
        if(!g_plugin_mode.CompareTo("Terminate"))
                plugin->SetMergeViaJDL(kFALSE);
        else
                plugin->SetMergeViaJDL(kTRUE);
        plugin->SetOverwriteMode();
        plugin->SetNtestFiles(1);

        plugin->SetAPIVersion("V1.1x");
        plugin->SetROOTVersion(g_root_version.Data());
        plugin->SetAliROOTVersion(g_aliroot_version.Data());

        plugin->SetOutputToRunNo();
        plugin->AddIncludePath("-I. .I$ALIEN_ROOT/api/lib -I$ROOTSYS/lib -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/PWGHF/ -I$ALICE_ROOT/PWGHF/hfe/macros -I$ALICE_ROOT/PWGHF/hfe -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/ANALYSIS/Tender -I$ALICE_ROOT/ANALYSIS/TenderSupplies -I$ALICE_ROOT/PWG/ -I$ALICE_ROOT/PWG/FLOW -I$ALICE_ROOT/PWG/Base -I$ALICE_ROOT/PWG/Tasks");
        plugin->SetAdditionalLibs("libGui.so  libXMLParser.so libSTEERBase.so libESD.so libAOD.so libCDB.so libANALYSIS.so libANALYSISalice.so libCORRFW.so  libPWGflowBase.so libPWGflowTasks.so libPWGHFhfe.so libTENDER.so libProof.so libRAWDatabase.so libSTEER.so libTOFbase.so");
   
        plugin->SetDefaultOutputs(kFALSE);
        plugin->SetOutputFiles("AnalysisResults.root"); 
        plugin->SetExecutableCommand("aliroot -b -q");
        plugin->SetTTL(30000);
        plugin->SetInputFormat("xml-single");
        plugin->SetPrice(1);      
        plugin->SetSplitMode("se");
        return plugin;
}

void SplitConfigEntry(const std::string &input, TString &key, TString &value){
        //
        // Decode key and value of a config entry
        //
        std::istringstream stream(input, istringstream::in);
        std::string tmp;
        stream >> tmp;
        key = tmp.c_str();
        stream >> tmp;
        value = tmp.c_str();
}

void DecodeRunlist(const TString &val){
        //
        // Tokenize run list
        //
        TObjArray *runstrings = val.Tokenize(",");
        TObjString *os;
        TString runstr;
        TIter runIter(runstrings);
        g_runlist.Set(runstrings->GetEntries());
        int nruns(0);
        while((os = dynamic_cast<TObjString *>(runIter()))){
                runstr = os->String();
                g_runlist[nruns++] = runstr.Atoi();
        }
        delete runstrings;
}

bool IsMC(const TString &val){
        // 
        // Determine whether sample is MC or Data
        //
        if(!val.CompareTo("MC")) return true;
        return false;
}

bool FindDataSample(const TMap &lookup, TObjArray &sampleinfis){
        //
        // Find Data sample in the list of samples
        //
        TObjArray *entry = dynamic_cast<TObjArray *>(lookup.GetValue(g_sample.Data()));
        if(!entry){
                printf("Sample %s not found in the list of samples", g_sample.Data());
                return false;
        }
        // Copy to output container
        sampleinfis.SetOwner(kFALSE);
        for(int ival = 0; ival < 4; ival++) sampleinfis.AddAt(entry->At(ival), ival);
        return true;
}

bool GetData(TObjArray &in, TString &out, int pos){
        //
        // Helper function reading data string
        //
        TObjString *entry = dynamic_cast<TObjString *>(in.At(pos));
        if(!entry){
                printf("Entry at pos %d not a string\n", pos);
                return false;
        }
        out = entry->String();
        return true;
}

void AddSample(TMap &lookup,
                const char *key, const char* datadir, const char * pattern, const char *sampletype, const char *dataformat){
        //
        // Add sample entry to the lookup table
        //
        TObjArray *infos = new TObjArray(); 
        infos->AddAt(new TObjString(datadir), 0);
        infos->AddAt(new TObjString(pattern), 1);
        infos->AddAt(new TObjString(sampletype), 2);
        infos->AddAt(new TObjString(dataformat), 3);
        lookup.Add(new TObjString(key), infos);
}

void Generate_Sample_Lookup(TMap &lookup){
        // 
        // Create Lookup table for each period
        // Vector contains 
        //   - path
        //   - pattern
        //   - MC/Data 
        //   - ESD/AOD
        //
        AddSample(lookup, "LHC13b.pass2", "/alice/data/2013/LHC13b", "*/pass2/*/AliESDs.root", "Data", "ESD");
        AddSample(lookup, "LHC13b.pass2.AOD", "/alice/data/2013/LHC13b", "*/pass2/AOD/*/AliAOD.root", "Data", "AOD"); 
        AddSample(lookup, "LHC13b.pass2.AOD126", "/alice/data/2013/LHC13b", "*/pass2/AOD126/*/AliAOD.root", "Data", "AOD"); 
        AddSample(lookup, "LHC13c.pass1", "/alice/data/2013/LHC13c/", "*/pass1/*/AliESDs.root", "Data", "ESD");
        AddSample(lookup, "LHC13c.pass1.AOD", "/alice/data/2013/LHC13c/", "*/pass1/AOD/*/AliAOD.root", "Data", "AOD");
        AddSample(lookup, "LHC13c.pass1.AOD126", "/alice/data/2013/LHC13c/", "*/pass1/AOD126/*/AliAOD.root", "Data", "AOD");
        AddSample(lookup, "LHC13b2", "/alice/sim/2013/LHC13b2", "*/*/AliESDs.root", "MC", "ESD");
        AddSample(lookup, "LHC13b2.AOD", "/alice/sim/2013/LHC13b2", "*/AliAOD.root", "MC", "AOD");
        AddSample(lookup, "LHC13b2plus", "/alice/sim/2013/LHC13b2_plus", "*/*/AliESDs.root", "MC", "ESD");
        AddSample(lookup, "LHC13b2plus.AOD", "/alice/sim/2013/LHC13b2_plus", "*/AliAOD.root", "MC", "AOD");
        AddSample(lookup, "LHC13b3", "/alice/sim/2013/LHC13b3", "*/*/AliESDs.root", "MC", "ESD");
        AddSample(lookup, "LHC13b3.AOD", "/alice/sim/2013/LHC13b3", "*/AliAOD.root", "MC", "AOD");
        printf("Lookup table with sample information generated\n");
}

void ConfigParser(const char *configname){
        //
        // Parse configuration
        //
        std::ifstream in(configname);
        std::string config;
        TString key, value;
        while(getline(in, config)){
                SplitConfigEntry(config, key, value);
                key.ToLower();
                if(!key.CompareTo("aliroot")){
                        // Aliroot version
                        g_aliroot_version = value;
                        continue;
                }
                if(!key.CompareTo("root")){
                        // root version
                        g_root_version = value;
                        continue;
                }
                if(!key.CompareTo("sample")){
                        // sample name
                        g_sample = value; 
                        continue;
                }
                if(!key.CompareTo("runlist")){
                        // Runlist
                        DecodeRunlist(value); 
                        continue;
                }
                if(!key.CompareTo("mode")){
                        g_plugin_mode = value;
                        continue;
                }
                if(!key.CompareTo("traindir")){
                        g_train_dir = value;
                        continue;
                }
                printf("Unknown key: %s\n", key.Data());
        }
}

bool MakeSample(AliAnalysisAlien *plugin, TMap &lookup){
        //
        // Fill Sample information (Data dir, pattern, run list) to the Alien plugin
        //
        TObjArray infos;
        bool found = FindDataSample(lookup, infos);
        if(!found){
                printf("sample %s not found\n", g_sample.Data());
                return false;
        }
        TString datadir, pattern, type;
        GetData(infos, datadir, 0);
        GetData(infos, pattern, 1);
        GetData(infos, type, 2);
        plugin->SetGridDataDir(datadir.Data());
        plugin->SetDataPattern(pattern.Data());
        if(!IsMC(type)) plugin->SetRunPrefix("000");
        // Add runs to the sample
        for(int irun = 0; irun < g_runlist.GetSize(); irun++){
                plugin->AddRunNumber(g_runlist[irun]);
        }
        return true;
}

bool CreateTrainDir(AliAnalysisAlien *plugin, const TMap &lookup){
        //
        // Make train data dir name and JDL, C and sh file names
        //
        TObjArray infos;
        bool found = FindDataSample(lookup, infos);
        if(!found){
                printf("sample %s not found\n", g_sample.Data());
                return false;
        }
        TString type; GetData(infos, type, 2);

        // check whether the train dir is already provided or needs to be specified
        if(!g_train_dir.Length()){
                // Query number of train runs before
                const char *gridhome = gGrid->GetHomeDirectory();
                const char gridoutdir[256];
                sprintf(gridoutdir, "%sAnalysis_pPb/%s", gridhome, type.Data());
                TGridResult *trainruns = gGrid->Ls(gridoutdir);
                int nruns = trainruns->GetEntries();
                // Get Date and time
                TDatime time;
                g_train_dir = Form("%d_%d%02d%02d_%02d%02d", nruns, time.GetYear(), time.GetMonth(), time.GetDay(), time.GetHour(), time.GetMinute());
        }
        
        plugin->SetGridWorkingDir(Form("Analysis_pPb/%s/%s", type.Data(), g_train_dir.Data()));
        plugin->SetJDLName(Form("TPCTOFanalysispPb_%s_%s.jdl", type.Data(), g_train_dir.Data()));
        plugin->SetExecutable(Form("TPCTOFanalysispPb_%s_%s.sh", type.Data(), g_train_dir.Data()));
        plugin->SetAnalysisMacro(Form("TPCTOFanalysispPb_%s_%s.C", type.Data(), g_train_dir.Data()));
        return true;
}

void SetupUtil(bool IsMC, bool isAOD){
        //
        // Setup utility packages
        //
        // 1. Physics Selection (ESD only)
        // 2. Tender (ESD only)
        // 3. PID Response (always)
        // 4. Centrality Task (ESD only)
        //

        //==== Physics Selection ====
        if(!isAOD){
                gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPhysicsSelection.C");
                AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);
        }

        //==== Add tender ====
        if(!isAOD){
                gROOT->LoadMacro("AddTaskTender.C");
                AddTaskTender();
        }

        //===== ADD PID RESPONSE: ===
        gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
        AddTaskPIDResponse(IsMC);

        //===== ADD CENTRALITY: ===
        if(!isAOD){
                gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskCentrality.C");
                AddTaskCentrality();
        }
}

void SetupHandlers(bool isMC, bool isAOD){
        //
        // Setup Handlers
        //
        TString macrobase = "$ALICE_ROOT/ANALYSIS/macros/train/";
        TString macroname = macrobase;
        if(isAOD)
                macroname += "AddAODHandler.C";
        else
                macroname += "AddESDHandler.C";
        gROOT->Macro(macroname.Data());

        if(isMC && !isAOD){
                // Add MC truth event handler, only in case of ESDs
                gROOT->LoadMacro(Form("%s/AddMCHandler.C", macrobase.Data()));
                AddMCHandler();
        }
}

void SetupHFEtask(bool isMC, bool isAOD){
        gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/AddTaskHFEpPb.C");
        AddTaskHFEpPb(isMC, isAOD);
        if(!isAOD){
                gROOT->LoadMacro("$ALICE_ROOT/PWGHF/hfe/macros/AddTaskHFEnpepPb.C");
                AddTaskHFEnpepPb();
        }
}

void SetupTrain(const TMap &lookup){
        //
        // Setup train:
        //   Determine whether the trains run on MC or Data 
        //   and ESDs or AODs and Configure Handlers, utils
        //   and HFE task according to this
        //
        bool isMC(false), isAOD(false);
        TObjArray infos;
        bool found = FindDataSample(lookup, infos);
        if(!found) return;
        TString type, mode;
        GetData(infos, type, 2);
        GetData(infos, mode, 3);
        isMC = IsMC(type);
        if(!mode.CompareTo("AOD")) isAOD = true;
      
        SetupHandlers(isMC, isAOD);
        SetupUtil(isMC, isAOD);
        SetupHFEtask(isMC, isAOD);
}

void GenerateMergeConfigs(){
        //
        // Generate configurations for merging 
        // (MergeViaJDL and Terminate)
        //

        // Write config for MergeViaJDL
        std::ofstream outMerge("configMerge.txt");
        outMerge << "aliroot " << g_aliroot_version.Data() << std::endl;
        outMerge << "root " << g_root_version.Data() << std::endl;
        outMerge << "sample " << g_sample.Data() << std::endl;
        outMerge << "mode MergeViaJDL\n";
        outMerge << "traindir " << g_train_dir.Data() << std::endl; 
        outMerge << "runlist ";
        for(int irun = 0; irun < g_runlist.GetSize()-1; irun++) outMerge << g_runlist[irun] << ",";
        outMerge << g_runlist[g_runlist.GetSize()-1] << std::endl;
        outMerge.close();
        // Write config for Terminate
        std::ofstream outTerminate("configTerminate.txt");
        outTerminate << "aliroot " << g_aliroot_version.Data() << std::endl;
        outTerminate << "root " << g_root_version.Data() << std::endl;
        outTerminate << "sample " << g_sample.Data() << std::endl;
        outTerminate << "mode Terminate\n";
        outTerminate << "traindir " << g_train_dir.Data() << std::endl; 
        outTerminate << "runlist ";
        for(int irun = 0; irun < g_runlist.GetSize()-1; irun++) Terminate << g_runlist[irun] << ",";
        outTerminate << g_runlist[g_runlist.GetSize()-1] << std::endl;
        outTerminate.close();

        printf("Configurations for MergeViaJDL and terminate generated\n");
}

void runGridpPb(const char *config = "config.txt"){
        //
        // run analysis 
        //

        TGrid::Connect("alien://");

        // Create Lookup with sample information
        TMap sampleinfos;
        Generate_Sample_Lookup(sampleinfos);

        ConfigParser(config);

        // Configure alien plugin
        AliAnalysisAlien *plugin = CreateGridHandler();
        if(!CreateTrainDir(plugin, sampleinfos)){
                printf("Cannot setup output directory\n");
                return;
        }
        if(!MakeSample(plugin, sampleinfos)){
                printf("Cannot create data sample\n");
                return;
        }
        if(!g_plugin_mode.CompareTo("full")){
                // full mode, creating config files for the merging stage
                GenerateMergeConfigs();
        }

        AliAnalysisManager *mgr = new AliAnalysisManager("tpctofanalysis");
        mgr->SetGridHandler(plugin);
        
        SetupTrain(sampleinfos);

        // Run train
        if (!mgr->InitAnalysis()) return;
        mgr->PrintStatus();
        // Start analysis in grid.
        mgr->StartAnalysis("grid");
} 

