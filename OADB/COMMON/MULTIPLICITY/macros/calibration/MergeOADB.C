//____________________________________________________________________________________________
void MergeOADB(const Char_t* fileList, const Char_t* outputFile, const Char_t* objName="MultSel") {
   //
   // Merge the OADB found in fileList
   //
   AliOADBContainer* oadb[1000] = {0x0};

   std::ifstream infile;
   infile.open(fileList);
   Int_t counter=0;
   while(infile.good()) {
      Char_t filename[1000];
      infile.getline(filename, 1000);
      if(filename[0]=='\0') continue;

      TFile* file=TFile::Open(filename);
      oadb[counter] = (AliOADBContainer*)file->Get(objName);
      counter++;
   }

   // Merge the oadb objects
   AliOADBContainer* mergedOADB = 0x0;
   for(Int_t i=0; i<counter; i++) {
      if(i==0) {
         mergedOADB = (AliOADBContainer*)oadb[i]->Clone(objName);
         continue;
      }
      for(Int_t k=0;k<oadb[i]->GetNumberOfEntries(); k++)
         mergedOADB->AppendObject(oadb[i]->GetObjectByIndex(k), oadb[i]->LowerLimit(k), oadb[i]->UpperLimit(k));

      TList* defaultList = oadb[i]->GetDefaultList();
      for(Int_t k=0;k<defaultList->GetEntries();++k)
         mergedOADB->AddDefaultObject(defaultList->At(k));
   }

   TFile* outFile = new TFile(outputFile, "RECREATE");
   mergedOADB->Write();
   outFile->Close();

}
