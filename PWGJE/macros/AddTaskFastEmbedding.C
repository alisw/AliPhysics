
AliAnalysisTaskFastEmbedding* AddTaskFastEmbedding(){

   AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
   if(!mgr){
      ::Error("AddTaskFastEmbedding", "No analysis manager to connect ot.");
      return NULL;
   }
   if(!mgr->GetInputEventHandler()){
      ::Error("AddTaskFastEmbedding", "This task requires an input event handler.");
      return NULL;
   }


   AliAnalysisTaskFastEmbedding *task = new AliAnalysisTaskFastEmbedding("FastEmbedding");
   // ## set embedding mode ##
   // kAODFull=0, kAODJetTracks, kAODJet4Mom, kToySingle4Mom
   //task->SetEmbedMode(AliAnalysisTaskFastEmbedding::kToyTracks);
   task->SetEmbedMode(AliAnalysisTaskFastEmbedding::kAODFull);      // embed full event: all tracks in PYTHIA event are added to PbPb data event 
   task->SetJetFriends("AliAOD.Jets.root");

   // ## set ranges for toy ##
   //SetToyTrackRanges(
   Double_t minPt = 50.;   Double_t maxPt = 100.;
   Double_t minEta = -0.5; Double_t maxEta = 0.5;
   //Double_t minEta = -0.4; Double_t maxEta = 0.4;  // for LHC10h pass1
   Double_t minPhi = 0.;   Double_t maxPhi = 2*TMath::Pi();
   //fToyDistributionTrackPt: 0 = uniform distribution
   //                         else = exponential / power law (not implemented yet)
   //task->SetToyNumberOfTrackRange(4,4);
   //task->SetToyTrackRanges(0.15, 300., 5,-.9, .9, 0., 2*TMath::Pi());
   task->SetToyTrackRanges(minPt,maxPt,0.,minEta,maxEta,minPhi,maxPhi);
   task->SetToyFilterMap((1<<32)-1);

   // ## set event selection for events of the addition AOD ##
   // kEventsAll=0; kEventsJetPt
   task->SetEvtSelecMode(AliAnalysisTaskFastEmbedding::kEventsJetPt);

   // event selection
   task->SetOfflineTrgMask(AliTrigger::kMB);
   task->SetCentMin(0.);
   task->SetCentMax(10.);

   //task->SetVtxMin(-10.);
   //task->SetVtxMax(10.);

   // ## set jet pT range for event selection ##
   // SetEvtSelJetPtRange(Float_t minPt, Float_t maxPt)
   task->SetEvtSelJetPtRange(5.,110.);
   //task->SetEvtSelJetEtaRange(-0.4, 0.4); // smaller eta window for LHC10h pass1
   task->SetEvtSelJetEtaRange(-0.5, 0.5);

   task->SetTrackFilterMap(272);

   task->SetJetBranch("clustersAOD_ANTIKT04_B0_Filter00272_Cut00150_Skip00");
   task->SetEvtSelJetPtRange(5.,10000.);  //jet min pt cut
   //task->SetEvtSelJetPtRange(10.,10000.);  //jet min pt cut
   task->SetEvtSelJetEtaRange(-0.5,0.5); 
   task->SetEvtSelJetMinLConstPt(1);
   task->SetFFRadius(0.2); //jet cone size

 //v0 cuts for embedded candidates


  Int_t K0type = AliAnalysisTaskFastEmbedding::kOffl;
  Int_t Latype = AliAnalysisTaskFastEmbedding::kOffl;
  Int_t ALatype = AliAnalysisTaskFastEmbedding::kOffl;

  TString strK0type;
  if(K0type ==  AliAnalysisTaskFastEmbedding::kOnFly)     strK0type = "OnFly";
  if(K0type ==  AliAnalysisTaskFastEmbedding::kOffl)      strK0type = "Offl";
  TString strLatype;
  if(Latype ==  AliAnalysisTaskFastEmbedding::kOnFly)     strLatype = "OnFly";
  if(Latype ==  AliAnalysisTaskFastEmbedding::kOffl)      strLatype = "Offl";
  
  TString strALatype;
  if(ALatype ==  AliAnalysisTaskFastEmbedding::kOnFly)     strALatype = "OnFly";
  if(ALatype ==  AliAnalysisTaskFastEmbedding::kOffl)      strALatype = "Offl";


  //pp V0 cut selection 

  task->SetK0Type(K0type);
  task->SetLaType(Latype); 
  task->SetALaType(ALatype); 

  task->SetCuttrackPosEta(0.8);
  task->SetCuttrackNegEta(0.8);
  task->SetCutV0Eta(0.7); //pseudorapidity cut, dont use 0.5, because too many tracks would fall out of the acceptance; recommended cut for jet analysis of strange particles: 0.75
  task->SetCosOfPointingAngle(0.97);//David:0.97//K0s, Lambda cut is stricter (0.995) and hardcoded inside task
  task->SetAcceptKinkDaughters(kFALSE);//accept kink daughters -> dont use this cut anymore
  task->SetRequireTPCRefit(kTRUE); 
  task->SetCutV0DecayMin(0.);//multiples of ctau, cut on 2D decay distance over transverse mom. (for K0s, Lambda, Antilambda)
  task->SetCutV0DecayMax(30.);//David: 20.(K0s), 30.(La)  //multiples of ctau (for K0s, Lambda, Antilambda) Lee Barnby uses 3.0, use 5.0!!!!!
  task->SetCutDcaV0Daughters(1.);//cut value in multiples of sigma default: 1.
  task->SetCutDcaPosToPrimVertex(0.06); //cut value in cm 
  task->SetCutDcaNegToPrimVertex(0.06); //cut value in cm
  task->SetCutV0RadiusMin(0.5);//0.5//in cm previous value was 0.9 cm
  task->SetCutV0RadiusMax(1000.);//in cm
 
  //task->SetCutArmenteros(0.2);//not needed in pp data/MC








   mgr->AddTask(task);

   // ## create the output containers ##
   AliAnalysisDataContainer *coutputFastEmbedding = mgr->CreateContainer(
   "fastembedding", TList::Class(), AliAnalysisManager::kOutputContainer,
   Form("%s:PWGJE_FastEmbedding", AliAnalysisManager::GetCommonFileName()));

   mgr->ConnectInput (task, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(task, 0, mgr->GetCommonOutputContainer());
   mgr->ConnectOutput(task, 1, coutputFastEmbedding);


   return task;

}

AliAnalysisTaskFastEmbedding* AddTaskFastEmbedding(TObjArray* aodarray){

   AliAnalysisTaskFastEmbedding *task = AddTaskFastEmbedding();
   if(aodarray){
      task->SetArrayOfAODPaths(aodarray);
      task->SetEmbedMode(AliAnalysisTaskFastEmbedding::kAODFull);
   }

   return task;
}


AliAnalysisTaskFastEmbedding* AddTaskFastEmbedding(const char* filepath, Int_t mode = 0){

   AliAnalysisTaskFastEmbedding *task = AddTaskFastEmbedding();
   if(strlen(filepath)){

      if(mode==0){ // path to single AOD
         task->SetAODPath(filepath);
      }
      if(mode==1){ // path to text file with list of paths of multiple AODs
         Printf("Read aod paths from file %s", filepath);
         TObjArray* array = new TObjArray();
         TObjString* ostr = 0;
         TString line;
         ifstream in;
         in.open(filepath);
         while(in.good()){
            in >> line;
            if(line.Length() == 0)  continue;
	    
            Printf("found aod path %s", line.Data());
            ostr = new TObjString(line.Data());
            array->Add(ostr);
         }
         Printf("-> %d aod paths found", array->GetEntries());
         
         task->SetArrayOfAODPaths(array);
      }
      if(mode==2) { //read root file which contains object array
	TFile *f = TFile::Open(filepath);
	TObjArray *objarray;
	f->GetObject("array",objarray);

        Printf("-> %d aod paths found", objarray->GetEntries());
	
	task->SetArrayOfAODPaths(objarray);


	Int_t count = 0;
	Int_t iEntry = -1;
	Int_t iEntrySum = 0;
	Int_t iEntryMax = 0;

	TArrayI* array = new TArrayI();

	for(int i=0; i<objarray->GetEntriesFast(); i++) {
	  TObjString *objStr = (TObjString*) objarray->At(i);
	  TString str = objStr->GetString();
	  iEntry = str.Atoi();
	  array->Set(count+1);
	  array->AddAt(iEntry,count);
	  count++;
	  iEntrySum += iEntry;
	  if(iEntry>iEntryMax) iEntryMax = iEntry;
	}
	
	task->SetArrayOfAODEntries(array);
	task->SetAODEntriesSum(iEntrySum);
	task->SetAODEntriesMax(iEntryMax);

      }

      task->SetEmbedMode(AliAnalysisTaskFastEmbedding::kAODFull);
   }

   return task;
}

AliAnalysisTaskFastEmbedding* AddTaskFastEmbedding(const char* aodpath, const char* entriespath){

   AliAnalysisTaskFastEmbedding *task = AddTaskFastEmbedding(aodpath, 1);

   Printf("Read entries of aod files from %s", entriespath);
   TArrayI* array = new TArrayI();
   Int_t count = 0;
   Int_t iEntry = -1;
   Int_t iEntrySum = 0;
   Int_t iEntryMax = 0;
   TString line;
   ifstream in;
   in.open(entriespath);
   while(in.good()){
      in >> line;
      iEntry = line.Atoi();

      array->Set(count+1);
      array->AddAt(iEntry,count);
      count++;
      iEntrySum += iEntry;
      if(iEntry>iEntryMax) iEntryMax = iEntry;
   }

   task->SetArrayOfAODEntries(array);
   task->SetAODEntriesSum(iEntrySum);
   task->SetAODEntriesMax(iEntryMax);

   return task;
}


