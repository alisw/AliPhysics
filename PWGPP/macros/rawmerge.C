//
// Macro to create the "raw" data file with selected events

void rawmerge( TString inputFileName="wn.xml",
               TString fullEventListFileName="event.list",
               TString outputFileNameTMP="filtered.root")
{
   // Create the filtered raw data files using a file with list of raw chunks with event numbers.
   // inputFileName         - either the text file with chunkname+event number+triggerType or 
   //                         xml collection
   // fullEventListFileName - if 1st arg is an xml collection, this is the full list of
   //                         chunks and events to be filtered acoording the xml contents
   //
   // format of the event list: /path/to/chunk eventNumber offlineTriggerType
   //
   //                     e.g.  /alice/data/2013/LHC13b/000195390/raw/13000195390000.10.root 218 V0s
   //                           /alice/data/2013/LHC13b/000195390/raw/13000195390000.10.root 219 highPt
   //                           /alice/data/2013/LHC13b/000195390/raw/13000195390000.10.root 246 V0s
   //
   // triggerType can be: highPt, V0s, CosmicPairs, Laser
   //
   //
   // if the file list is an xml collection (for running on alien),
   // first extract the available chunks and event numbers from the
   // reference file, and use that as input
   if (inputFileName.Contains(".xml"))
   {
      TString tmp=inputFileName+".list";
      makeAlienInputEventList(tmp,fullEventListFileName,inputFileName);
      inputFileName=tmp;
   }

   TGrid * alien = TGrid::Connect("alien://",0,0,"t");

   Int_t eventNumber;
   ifstream files;
   files.open(inputFileName.Data());
   if (!files.is_open())
   {
      fprintf(stderr,"error: could not read event list file \"%s\". Exiting.\n",inputFileName.Data());
      return;
   }
   
   //TSeqCollection* listOfFiles = (TSeqCollection*)gROOT->GetListOfFiles();
   TList* listOfFiles = new TList();
   TString outputFileName;
   TString line;
   TString iURI;
   TString triggerType;
   TString  iURIold;
   TFile *ifile=0;
   TFile *ofile=0;
   TTree *itree=0;
   TTree *otree=0;
   Long64_t ievent=0;
   Int_t ofilenumber=0;
   Int_t lineNumber=0;
   Int_t eventold=0;

   while (files.good())
   {
      ++lineNumber;
      //read the line, do some checks
      line.ReadLine(files);
      TObjArray* entries = line.Tokenize(" ");
      TObjString* iURIobjstr = (TObjString*)entries->At(0);
      iURI.Clear();
      if (iURIobjstr) iURI=iURIobjstr->String();
      TObjString* ieventobjstr = (TObjString*)entries->At(1);
      if (ieventobjstr) ievent = ieventobjstr->String().Atoi();
      if (iURI.IsNull() || !ieventobjstr)
      {
        printf("malformed line: %s, skipping...\n",line.Data());
        continue;
      }
      TObjString* triggerTypeobjstr = (TObjString*)entries->At(2);
      triggerType.Clear();
      if (triggerTypeobjstr) triggerType=triggerTypeobjstr->String();

      printf("> processing \"%s\" event %d trigger \"%s\"...\n",iURI.Data(),ievent,triggerType.Data());
      if (ievent==eventold && iURI==iURIold)
      {
         printf("duplicated continue\n");
         continue;
      }
      //
      if (!iURIold.Contains(iURI.Data()))
      {
         //if a new file - open it
         printf("new file: %s\n",iURI.Data());
         delete ifile;
         ifile=0;
         ifile=TFile::Open(iURI.Data());
         if (!ifile)
         {
            fprintf(stderr,"warning: could not open file for event \"%s\", skipping it...\n",iURI.Data());
            continue;
         }
      }
      else
      {
         //if same file, reuse it
         printf("using already open file: %s\n",iURI.Data());
      }
      iURIold=iURI;
      //
      TTree *itree=dynamic_cast<TTree*>(ifile->Get("RAW"));
      if (!itree)
      {
         fprintf(stderr,"warning: could not find RAW tree for event \"%s\", skipping it...\n",iURI.Data());
         continue;
      }

      // manage output files
      if (!triggerType.IsNull()) triggerType.Prepend("_");
      outputFileName=outputFileNameTMP;
      outputFileName.ReplaceAll(".root","");
      outputFileName+=triggerType;
      outputFileName+=".root";
      
      ofile=dynamic_cast<TFile*>(listOfFiles->FindObject(outputFileName.Data()));
      if (!ofile) 
      {
        printf("< creating output file \"%s\"\n",outputFileName.Data());
        ofile=TFile::Open(outputFileName,"recreate");
        if (!ofile)
        {
          fprintf(stderr,"error: could not create output file: \"%s\" Exiting.\n",outputFileName.Data());
          return;
        }
        listOfFiles->Add(ofile);
      }
      ofile->cd();
      otree=dynamic_cast<TTree*>(ofile->Get("RAW"));
      if (!otree) 
      {
        otree=itree->CloneTree(0);
      }

      otree->CopyAddresses(itree);
      itree->GetEntry(ievent);
      eventold=ievent;
      printf("filling event %i in file %s\n",ievent,ofile->GetName());
      otree->Fill();
      //otree->CopyEntries(itree,Form("Entry$==%d",ievent),1);

      // reset input
      itree->ResetBranchAddresses();
   }

   //close the files
   for (Int_t i=0; i<listOfFiles->GetEntries(); i++)
   {
     ofile=dynamic_cast<TFile*>(listOfFiles->At(i));
     if (!ofile) {continue;}
     otree=dynamic_cast<TTree*>(ofile->Get("RAW"));
     Long64_t nEntries=0;
     if (otree) { nEntries = otree->GetEntries(); }

     //write the file and close
     ofile->Write();
     printf("closing file: %s with %i entries\n",ofile->GetName(),nEntries);
     delete ofile;

     //remove empty files
     if (nEntries==0)
     {
       gSystem->Unlink(ofile->GetName());
     }
   }
}

Bool_t makeAlienInputEventList( TString outputFileName="wn.list",
                                TString referenceFileName="filteredEvents.list",
                                TString inputCollectionName="wn.xml" )
{
   TGridCollection *coll = gGrid->OpenCollection(inputCollectionName);
   if (!coll)
   {
      ::Error("makeAlienInputEventList", "Cannot open collection from %s", xmlfile);
      return NULL;
   }
   TString configLine;
   TString chunkPath;
   TString chunkName;
   ifstream referenceFile;
   referenceFile.open(referenceFileName.Data());
   if (!referenceFile.is_open())
   {
      printf("could not open %s\n",referenceFileName.Data());
      return kFALSE;
   }
   gSystem->Unlink(outputFileName);
   ofstream outputFile;
   outputFile.open(outputFileName.Data());
   if (!outputFile.is_open())
   {
      printf("could not open %s\n",referenceFileName.Data());
      return kFALSE;
   }

   while (coll->Next())
   {
      chunkPath = coll->GetTURL();
      TObjArray* a = chunkPath.Tokenize("/");
      TObjString* objstr = dynamic_cast<TObjString*>(a->At(a->GetEntries()-1));
      if (!objstr)
      {
         printf("empty chunkPath from collection!\n");
         return kFALSE;
      }
      TString chunkName = objstr->String();
      chunkName.ReplaceAll("alien://","");
      while (referenceFile.good())
      {
         configLine.ReadLine(referenceFile);
         if (configLine.Contains(chunkName))
         {
            configLine=configLine.Strip(TString::kBoth);
            if (!configLine.BeginsWith("alien://"))
            {
               configLine.Prepend("alien://");
            }
            outputFile << configLine << endl;
            cout << configLine << endl;
         }
      }
      //jump to beginning and clear flags
      referenceFile.clear();
      referenceFile.seekg(0,ios::beg);
   }
   outputFile.close();
   referenceFile.close();
   return kTRUE;
}

