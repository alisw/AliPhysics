// Grid running parameters
TString gGridRunMode = "full";
TString gRootVersion = "v5-34-30-alice7-6";
TString gAlirootVersion = "v5-09-03a-1";
TString gAliphysicsVersion = "vAN-20170804-1";
//TString gGridDataDir = "/alice/data/2015/LHC15o";
//TString gGridDataDir = "/alice/data/2016/LHC16q";
TString gGridDataDir = "/alice/data/2016/LHC16s";
//TString gGridDataDir = "/alice/sim/2017/LHC17d13_cent/";
//TString gGridDataDir = "/alice/cern.ch/user/i/iarsene/work/outputDst";
//TString gGridDataPattern = "*/pass1/*/AliESDs.root";
//TString gGridDataPattern = "*/pass1/PWGDQ/DQ_PbPb/231_20161009-2048/*/dstTree.root";
TString gGridDataPattern = "*/pass1_FAST/*/AliESDs.root";
//TString gGridDataPattern = "*/AOD/*/AliAOD.root";
//TString gGridDataPattern = "*/pass1_CENT_wSDD/AOD/*/AliAOD.root";
TString gGridWorkingDir = "work";
TString gGridOutputDir = "20170805_dstTrees_LHC16rs";
Int_t gGridMaxInputFileNumber = 50;

TChain* makeChain(const Char_t* filename, const Char_t* inputType);

//______________________________________________________________________________________________________________________________________
void runAnalysisTrain_EMCalCorrFramework(const Char_t* infile, const Char_t* runmode = "local", const Char_t* inputType="ESD", Bool_t hasMC = kFALSE,
                                         Int_t reducedEventType = -1, Bool_t writeTree = kFALSE, TString tasks="dst", TString prod = "LHC10h",
                                         Int_t nEntries=-1, Int_t firstEntry=0, TString pathForMacros="$ALICE_PHYSICS/PWGDQ/reducedTree/macros")
{
   //
   // infile: list of input files if mode is local, or list of runs for job submission if mode is grid
   //
   // Load common libraries
   gSystem->Load("libPWGDQdielectron.so"); 
   gSystem->Load("libPWGDQreducedTree.so"); 

   // Use AliRoot includes to compile our task
   gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
   gROOT->ProcessLine(".include $ALICE_ROOT/include");
   gROOT->ProcessLine(".include $ROOTSYS/include");
   
   // Create the analysis manager
   AliAnalysisManager *mgr = new AliAnalysisManager("ReducedTreeAnalysis");
   
   // Create and configure the alien handler plugin
   TString runmodestr(runmode); runmodestr.ToLower();
   AliAnalysisGrid *alienHandler = NULL;
   if(runmodestr.Contains("grid")) {
     // Detect which tasks are being run and add the needed grid output files
     TObjArray* arr = tasks.Tokenize(";");
     TString outputFiles = "";
     if(writeTree) outputFiles = "dstTree.root";
     if((!writeTree && arr->GetEntries()>0) ||
        (writeTree && arr->GetEntries()>1)) outputFiles += " dstAnalysisHistograms.root";
#ifdef __CLING__
     // ROOT6 version
     std::stringstream creatalienhandleradd;
     creatalienhandleradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWGDQ/reducedTree/macros/CreateAlienHandler.C(";
     creatalienhandleradd << "\"" << infile << "\", ";
     creatalienhandleradd << "\"" << gGridRunMode.Data() << "\", ";
     creatalienhandleradd << "\"" << gGridDataDir.Data() << "\", ";
     creatalienhandleradd << "\"" << gGridDataPattern.Data() << "\", ";
     creatalienhandleradd         << gGridMaxInputFileNumber << ", ";
     creatalienhandleradd << "\"" << gGridWorkingDir.Data() << "\", ";
     creatalienhandleradd << "\"" << gGridOutputDir.Data() << "\", ";
     creatalienhandleradd << "\"" << outputFiles.Data() << "\", ";
     creatalienhandleradd << "\"" << gRootVersion.Data() << "\", ";
     creatalienhandleradd << "\"" << gAlirootVersion.Data() << "\", ";
     creatalienhandleradd << "\"" << gAliphysicsVersion.Data() << "\")";
     std::string creatalienhandleraddstr = creatalienhandleradd.str();
     std::cout << "Calling Add macro using command string: " << creatalienhandleraddstr << std::endl;
     alienHandler = (AliAnalysisGrid*)gROOT->ProcessLine(creatalienhandleraddstr.c_str());
#else
     // ROOT5 version
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/reducedTree/macros/CreateAlienHandler.C");
     alienHandler = CreateAlienHandlerPbPb(infile, gGridRunMode, gGridDataDir, gGridDataPattern, gGridMaxInputFileNumber,
                                           gGridWorkingDir, gGridOutputDir, outputFiles,
                                           gRootVersion, gAlirootVersion, gAliphysicsVersion);
#endif
     if (!alienHandler) {
        cout << "runAnalysisTrain.C ::      Could not create the alien handler. Check it out!" << endl;
        return;
     }
     mgr->SetGridHandler(alienHandler);
   }

   // Create the input handler
   TString inputTypeStr(inputType); inputTypeStr.ToLower();
   AliInputEventHandler* inputHandler = NULL; 
   if(inputTypeStr.Contains("esd")) {                               // ESDs
     inputHandler = new AliESDInputHandler();
     ((AliESDInputHandler*)inputHandler)->SetReadFriends(kFALSE);
   }
   if(inputTypeStr.Contains("aod")) {                               // AODs
      inputHandler = new AliAODInputHandler();
   }
   if(inputTypeStr.Contains("reducedevent")) {              // AliReducedEventInfo
      inputHandler = new AliReducedEventInputHandler();
      ((AliReducedEventInputHandler*)inputHandler)->SetInputEventType(AliReducedEventInputHandler::kReducedEventInfo);
   }
   if(inputTypeStr.Contains("baseevent")) {                    // AliReducedBaseEvent
      inputHandler = new AliReducedEventInputHandler();
      ((AliReducedEventInputHandler*)inputHandler)->SetInputEventType(AliReducedEventInputHandler::kReducedBaseEvent);
   }
   mgr->SetInputEventHandler(inputHandler);
   
   // Add the MC handler if needed
   AliMCEventHandler *mc = NULL;
   if(hasMC && (inputTypeStr.Contains("esd") || inputTypeStr.Contains("aod"))) {
      mc = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mc);
   }
   
   //==== Add tender ====
   //   gROOT->LoadMacro("AddTaskTender.C");
   //   AddTaskTender();

   if(inputTypeStr.Contains("esd") || inputTypeStr.Contains("aod")) {         // no need if we run over reduced events
     //==== Physics Selection ====
     Bool_t applyPileupCuts = kTRUE;
#ifdef __CLING__
     // ROOT6 version
     std::stringstream physseladd;
     physseladd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/OADB/macros/AddTaskPhysicsSelection.C(";
     physseladd << hasMC << ", ";
     physseladd << applyPileupCuts << ")";
     std::string physseladdstr = physseladd.str();
     std::cout << "Calling Add macro using command string: " << physseladdstr << std::endl;
     AliPhysicsSelectionTask* physSelTask = (AliPhysicsSelectionTask*)gROOT->ProcessLine(physseladdstr.c_str());
#else
     // ROOT5 version
     gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
     AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(hasMC, applyPileupCuts);
#endif
      //===== ADD CENTRALITY: ===
      if(!prod.CompareTo("LHC10h") || !prod.CompareTo("LHC11h")) {         // Run-1 Pb-Pb
#ifdef __CLING__
        // ROOT6 version
        std::stringstream centradd;
        centradd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/OADB/macros/AddTaskCentrality.C";
        std::string centraddstr = centradd.str();
        std::cout << "Calling Add macro using command string: " << centraddstr << std::endl;
        gROOT->ProcessLine(centraddstr.c_str());
#else
        // ROOT5 version
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
        AddTaskCentrality();
#endif
      }
      else  { // Run-2
      //if(!prod.CompareTo("LHC15o") || !prod.CompareTo("LHC16l") || !prod.CompareTo("LHC16q") || !prod.CompareTo("LHC16t")) {         // Run-2 Pb-Pb
#ifdef __CLING__
        // ROOT6 version
        std::stringstream multsel;
        multsel << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C";
        std::string multselstr = multsel.str();
        std::cout << "Calling Add macro using command string: " << multselstr << std::endl;
        AliMultSelectionTask* multTask = (AliMultSelectionTask*)gROOT->ProcessLine(multselstr.c_str());
#else
        // ROOT5 version
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
        AliMultSelectionTask* multTask = AddTaskMultSelection();
#endif
         if(hasMC && !prod.CompareTo("LHC15o"))
           multTask->SetAlternateOADBforEstimators("LHC15o-DefaultMC-HIJING");
         if(hasMC && (!prod.CompareTo("LHC16q") || !prod.CompareTo("LHC16t")))
            multTask->SetAlternateOADBforEstimators("LHC16q-DefaultMC-HIJING");
         if(hasMC && (!prod.CompareTo("LHC16r")))
            multTask->SetAlternateOADBforEstimators("LHC16r-DefaultMC-EPOSLHC");
         if(hasMC && (!prod.CompareTo("LHC16s")))
            multTask->SetAlternateOADBforEstimators("LHC16s-DefaultMC-EPOSLHC");
      }      

     //===== ADD PID RESPONSE: ===
     Bool_t tuneOnData = kTRUE;
     TString recoPass = "1";
#ifdef __CLING__
     // ROOT6 version
     std::stringstream pidresp;
     if (hasMC) {
       pidresp << ".x " << gSystem->Getenv("ALICE_ROOT") << "/ANALYSIS/macros/AddTaskPIDResponse.C(";
       pidresp << hasMC << ", ";
       pidresp << kTRUE << ", ";
       pidresp << tuneOnData << ", ";
       pidresp << recoPass << ")";
     } else {
       pidresp << ".x " << gSystem->Getenv("ALICE_ROOT") << "/ANALYSIS/macros/AddTaskPIDResponse.C";
     }
     std::string pidrespstr = pidresp.str();
     std::cout << "Calling Add macro using command string: " << pidrespstr << std::endl;
     gROOT->ProcessLine(pidrespstr.c_str());
#else
     // ROOT5 version
     gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
     if(hasMC) AddTaskPIDResponse(hasMC, kTRUE, tuneOnData, recoPass);
     else AddTaskPIDResponse();
#endif

     //===== ADD CDB CONNECT: ===
     // NOTE: required for EMCal correction framework called below
#ifdef __CLING__
     // ROOT6 version
     std::stringstream cdbconnect;
     cdbconnect << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWGPP/PilotTrain/AddTaskCDBconnect.C";
     std::string cdbconnectstr = cdbconnect.str();
     std::cout << "Calling Add macro using command string: " << cdbconnectstr << std::endl;
     AliTaskCDBconnect *taskCDB = (AliTaskCDBconnect*)gROOT->ProcessLine(cdbconnectstr.c_str());
#else
     // ROOT5 version
     gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
     AliTaskCDBconnect *taskCDB = AddTaskCDBconnect();
#endif
     taskCDB->SetFallBackToRaw(kTRUE);

     //===== ADD EMCAL CORRECTION: ===
#ifdef __CLING__
     // ROOT6 version
     std::stringstream emcalcorr;
     emcalcorr << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWG/EMCAL/macros/AddTaskEmcalCorrectionTask.C";
     std::string emcalcorrstr = emcalcorr.str();
     std::cout << "Calling Add macro using command string: " << emcalcorrstr << std::endl;
     AliEmcalCorrectionTask *correctionTask = (AliEmcalCorrectionTask*)gROOT->ProcessLine(emcalcorrstr.c_str());
#else
     // ROOT5 version
     gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalCorrectionTask.C");
     AliEmcalCorrectionTask *correctionTask = AddTaskEmcalCorrectionTask();
#endif
     UInt_t kPhysSel = AliVEvent::kAny;
     correctionTask->SelectCollisionCandidates(kPhysSel);
     correctionTask->SetUserConfigurationFilename("$ALICE_PHYSICS/PWGHF/hfe/macros/configs/pp/userConfigurationEMCele_pp_pPb.yaml");
     correctionTask->Initialize();
   }

#ifdef __CLING__
  // ROOT6 version
  std::stringstream trainadd;
  trainadd << ".x " << gSystem->Getenv("ALICE_PHYSICS") << "/PWGDQ/reducedTree/macros/AddTask_TrainTreeAnalysis.C(";
  trainadd << (!runmodestr.CompareTo("grid") ? kTRUE : kFALSE) << ", ";
  trainadd << "\"" << prod.Data() << "\", ";
  trainadd << reducedEventType << ", ";
  trainadd << writeTree << ", ";
  trainadd << "\"" << tasks.Data() << "\", ";
  trainadd << "\"" << pathForMacros.Data() << "\")";
  std::string trainaddstr = trainadd.str();
  std::cout << "Calling Add macro using command string: " << trainaddstr << std::endl;
  gROOT->ProcessLine(trainaddstr.c_str());
#else
  // ROOT5 version
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/reducedTree/macros/AddTask_TrainTreeAnalysis.C");
  AddTask_TrainTreeAnalysis((!runmodestr.CompareTo("grid") ? kTRUE : kFALSE), prod, reducedEventType, writeTree, tasks, pathForMacros);
#endif

   // Enable debug printouts
   //mgr->SetDebugLevel(10);

   if (!mgr->InitAnalysis()) return;

   TChain* chain = NULL;
   if(!runmodestr.Contains("grid"))
      chain = makeChain(infile, inputType);

   TProof* proof=0x0;
   if(runmodestr.Contains("proof")) {
      proof = TProof::Open("");
      chain->SetProof();
   }

   mgr->PrintStatus();
   // Start analysis
   if(nEntries==-1) nEntries=1234567890;
   if(runmodestr.Contains("local"))
      mgr->StartAnalysis("local", chain, nEntries, firstEntry);
   if(runmodestr.Contains("proof"))
      mgr->StartAnalysis("proof", chain, nEntries, firstEntry);
   if(runmodestr.Contains("grid"))
      mgr->StartAnalysis("grid", nEntries, firstEntry);


	//--------------------------------------------
	// piping event stat histograms from tree file(s) into analysis output(s)
	Bool_t verbose = kFALSE;
	cout << endl;
	TObjArray* taskArr = tasks.Tokenize(";");
	TList* listEventStat = new TList();
	TObjArray* histArr = new TObjArray();

	if (!writeTree && taskArr->GetEntries()>0) {
		// on the fly mode
		if(((TString)taskArr->At(0)->GetName()).Contains("dst"))
			cout << "----------------------" << endl << "On the fly mode: no piping of the TList of event stat histograms into analysis output" << endl << "----------------------" << endl;
		// Analysis running on ReducedTree mode
		else {
			cout << "----------------------" << endl << "Piping the TList of event stat histograms into analysis output(s) from the list of ReducedTree files:" << endl << "----------------------" << endl;
			ifstream in;
			in.open(infile);

			TString line;
			Int_t nFile = 0;
			// loop over input files
			while(in.good()) {
				in >> line;
				if (!line.IsNull()) {
					nFile++;
					TFile* tmpFile = TFile::Open(Form("%s",line.Data()));
					cout << endl << "Looking for list of histograms in " << tmpFile->GetName() << ":" << endl;
					Int_t isOK = 1;
					Bool_t islist = kFALSE;
					Bool_t isOld = kFALSE;
					TIter next(tmpFile->GetListOfKeys());
					TKey* key;
					// loop over keys in file
					while ((key = (TKey*)next())) {
						TClass *cl = gROOT->GetClass(key->GetClassName());
						// Backward compatibilty with former input-file structure
						if (cl->InheritsFrom("TH1")) {
							if (!isOld)
								cout << "INFO: You are using the old input-file structure (ie w/o TList) ==> Please consider to update." << endl;
							isOld = kTRUE;
							TH1 *h = (TH1*)key->ReadObj();
							if(verbose) cout << " " << h->GetName();
							if (nFile==1) {
								if(verbose) cout << " -> added" << endl;
								histArr->Add(h);
							} else {
								if (histArr->FindObject(h->GetName())) {
									if(verbose) cout << " -> merged" << endl;
									((TH1*)histArr->FindObject(h->GetName()))->Add(h);
								} else {
									if(verbose) cout << " -> added" << endl;
									histArr->Add(h);
								}
							}
						}
						// looking for TList "EventStats"
						else if (!cl->InheritsFrom("TList")) continue;
						TList *list = (TList*)key->ReadObj();
						TString lname = list->GetName();
						if(!lname.Contains("EventStats")) continue;
						islist = kTRUE;
						if(verbose) cout << "Input list \"" << lname << "\" found:" << endl;
						// first file --> clone TList
						if (nFile==1) {
							listEventStat = (TList*)list->Clone(list->GetName());
							// verbose
							if(verbose) cout << " -> added histogram(s): " << endl;
							TIter lnext(list);
							while (TObject *lobj = lnext()) {
								if(((TString)lobj->ClassName()).Contains("TList")) {
									if(verbose) cout << "    Input sublist \"" << lobj->GetName() << "\" found:" << endl;
									TIter sublnext((TList*)lobj);
									while (TObject *sublobj = sublnext()) {
										if (!((TString)sublobj->ClassName()).Contains("TH2")) {
											cout << "WARNING: Object \"" << sublobj->GetName() <<"\" of type \""<< sublobj->ClassName() << "\" unexpected. Ignored!!" << endl;
											isOK = -1;
											continue;
										}
										if(verbose) cout << "    " << sublobj->GetName() << endl;
									}
								}
								else if(!((TString)lobj->ClassName()).Contains("TH2")) {
									cout << "WARNING: Object \"" << lobj->GetName() <<"\" of type \""<< lobj->ClassName() << "\" unexpected. Ignored!!" << endl;
									isOK = -1;
									continue;
								}
								else if(verbose) cout << " " << lobj->GetName() << endl;
							}
						}
						// other files --> Add histograms
						else {
							if(!listEventStat->GetEntries()) listEventStat->SetName(list->GetName());
							Bool_t isHere = kTRUE;
							TIter lnext(list);
							// loop over the objects of the input TList to Add objects
							while (TObject *lobj = lnext()) {
								// excisting object
								if (listEventStat->FindObject(lobj->GetName())) {
									TString clname = lobj->ClassName();
									// TH2
									if (clname.Contains("TH2")) {
										if(verbose) {
											if(isHere) cout << " -> merged histogram(s): " << endl;
											cout << " " << lobj->GetName() << endl;
										}
										isHere = kFALSE;
										((TH2*)listEventStat->FindObject(lobj->GetName()))->Add((TH2*)lobj);
									}
									// sub TList
									else if(clname.Contains("TList")) {
										if(verbose) cout << "    Input sublist \"" << lobj->GetName() << "\" found:" << endl;
										TIter sublnext((TList*)lobj);
										// loop over the objects of the sub TList
										while (TObject *sublobj = sublnext()) {
											// looking for TH2 only
											if (!((TString)sublobj->ClassName()).Contains("TH2")) {
												cout << "    WARNING: Unexpected type of object: " << clname << ". Ignored!!" << endl;
												isOK = -1;
												continue;
											}
											// excisting object
											if (listEventStat->FindObject(lobj->GetName())) {
												if(verbose) {
													if(isHere) cout << "     -> merged histogram(s): " << endl;
													cout << "    " << sublobj->GetName() << endl;
												}
												isHere = kFALSE;
												((TH2*)((TList*)listEventStat->FindObject(lobj->GetName()))->FindObject(sublobj->GetName()))->Add((TH2*)sublobj);
											}
											// adding new object
											else {
												isHere = kTRUE;
												isOK = 0;
												cout << "    ERROR: additional object found wrt to previous input files. Check carrefully the consistency of your inputs." << endl;
												cout << "     -> added histogram: " << sublobj->GetName() << endl;
												((TList*)listEventStat->FindObject(lobj->GetName()))->Add((TH2*)sublobj);
											}
										}
									}
									else {
										isOK = -1;
										cout << "WARNING: Unexpected type of object: " << clname << ". Ignored!!" << endl;
									}
								}
								// adding new object
								else {
									isHere = kTRUE;
									isOK = 0;
									cout << "ERROR: additional object found wrt to previous input files. Check carrefully the consistency of your inputs." << endl;
									// TH2
									if(((TString)lobj->ClassName()).Contains("TH2")) {
										cout << " -> added histogram: " << lobj->GetName() << endl;
										listEventStat->Add((TH2*)lobj);
									}
									// sub TList
									else if(((TString)lobj->ClassName()).Contains("TList")) {
										cout << "    Input sublist \"" << lobj->GetName() << "\" found:" << endl;
										TList* sublistEventStat = new TList();
										sublistEventStat->SetName(lobj->GetName());
										TIter sublnext((TList*)lobj);
										// loop over the object of the sub TList
										while (TObject *sublobj = sublnext()) {
											// looking for TH2 only
											if (!((TString)sublobj->ClassName()).Contains("TH2")) {
												cout << "    WARNING: Unexpected type of object: " << sublobj->ClassName() << ". Ignored!!" << endl;
												continue;
											}
											cout << "     -> added histogram: " << sublobj->GetName() << endl;
											// add object to the sub TLIst
											sublistEventStat->Add((TH2*)sublobj);
										}
										// add sub TList to the TList
										listEventStat->Add(sublistEventStat);
									}
									else cout << "WARNING: Unexpected type of object. Ignored!!" << endl;
								}
							}
							// loop over the new TList objects to check if some are missing in the input TList wrt previous files
							TIter newlnext(listEventStat);
							// loop over the objects of the new TList
							while (TObject *lobj = newlnext()) {
								if (!list->FindObject(lobj->GetName())) {
									isOK = 0;
									TString clname = lobj->ClassName();
									// TH2
									if (clname.Contains("TH2")) cout << "ERROR: missing histogram \"" << lobj->GetName() << "\" in the input list wrt previous files!! Check carrefully the consistency of your inputs." << endl;
									// sub TList
									else if(clname.Contains("TList")) {
										TIter sublnext((TList*)lobj);
										// loop over the objects of the sub TList
										while (TObject *sublobj = sublnext()) {
											// looking for TH2 only
											if (!((TString)sublobj->ClassName()).Contains("TH2")) continue;
											if (!list->FindObject(sublobj->GetName())) cout << "ERROR: missing histogram \"" << lobj->GetName() << "\" (sublist \"" << lobj->GetName() << "\") in the input list wrt previous files!! Check carrefully the consistency of your inputs." << endl;
										}
									}
									else {
										isOK = -1;
										cout << "WARNING: missing unexpected object of type \"" << clname << "\". Ignored!!" << endl;
									}
								}
							}
						}
					}
					// no TList
					if(!islist && !isOld) {
						isOK = 0;
						cout << "ERROR: no input TList of histograms found!! Check carrefully the consistency of your inputs." << endl;
					}
					if (isOK==1) cout << " ---> OK" << endl;
					else if(isOK==-1)cout << " ---> OK, with unexpected types of object" << endl;
				}
			}

			// write list of histograms to output file
			TObjArray* outputs = mgr->GetOutputs();
			for (Int_t i=0; i<outputs->GetEntries(); i++) {
				TFile* tmpOut = (TFile*)mgr->OpenFile((AliAnalysisDataContainer*)outputs->At(i), "UPDATE", kTRUE);
				if (!tmpOut) continue;
				TString tmpOutName = tmpOut->GetName();
				if(tmpOutName.Contains("dstTree")) continue;
				cout << endl;
				if (histArr->GetEntries() && listEventStat->GetEntries()) cout << "ERROR: Two different input structures found!! Check carrefully the consistency of your inputs." << endl;
				if (histArr->GetEntries()) {
					cout << "==> Histograms written into " << tmpOut->GetName() << endl;
					for (Int_t j=0; j<histArr->GetEntries(); j++) histArr->At(j)->Write();
				}
				if(listEventStat->GetEntries()) {
					cout << "==> TList of histograms \"" << listEventStat->GetName() << "\" written into " << tmpOut->GetName() << endl;
					listEventStat->Write(listEventStat->GetName(),1);
				}
				tmpOut->Close();
			}
		}
	}

	if(writeTree && taskArr->GetEntries()>1) { // TreeMaker + AnalysisTask mode
		cout << "----------------------" << endl << "Piping the TList of event stat histograms into analysis output(s) from written tree file:" << endl << "----------------------" << endl;
		// Copy list of histograms
		TObjArray* outputs = mgr->GetOutputs();
		for (Int_t i=0; i<outputs->GetEntries(); i++) {
			TFile* tmpOut = (TFile*)mgr->OpenFile((AliAnalysisDataContainer*)outputs->At(i), "UPDATE", kTRUE);
			if (!tmpOut) continue;
			TString tmpOutName = tmpOut->GetName();
			if(tmpOutName.Contains("dstTree")) {
				TIter next(tmpOut->GetListOfKeys());
				// loop over keys in file dstTree
				while (TKey* key = (TKey*)next()) {
					TClass *cl = gROOT->GetClass(key->GetClassName());
					if(!cl->InheritsFrom("TList")) continue;
					TList *list = (TList*)key->ReadObj();
					if(!((TString)list->GetName()).Contains("EventStats")) continue;
					listEventStat = (TList*)list->Clone();
				}
				continue;
			}
			// Write list of histograms
			if(listEventStat->GetEntries()) {
				cout << "==> TList of histograms \"" << listEventStat->GetName() << "\" written into " << tmpOut->GetName() << endl;
				listEventStat->Write(listEventStat->GetName(),1);
			}
			tmpOut->Close();
		}
	}

};

//_______________________________________________________________________________
TChain* makeChain(const Char_t* filename, const Char_t* inputType) {
  //
  // make a chain using the trees from the list of filenames in "filename"
  //
  TString itStr(inputType); itStr.ToLower();
  TChain *chain = NULL;
  if(itStr.Contains("reducedevent") || itStr.Contains("baseevent"))
    chain=new TChain("DstTree");
  if(itStr.Contains("esd"))
    chain=new TChain("esdTree");
  if(itStr.Contains("aod"))
     chain=new TChain("aodTree");
  
  ifstream in;
  in.open(filename);
                                                                                                                                                               
  TString line;
  //loop over file                                                                                                                                                             
  while(in.good()) {
    in >> line;
    if (!line.IsNull()) {
      cout << "Adding file: " << line <<endl;
      chain->AddFile(line);
    }
  }
  cout << "Number of events in chain: " << chain->GetEntries() << endl;
  return chain;
}
