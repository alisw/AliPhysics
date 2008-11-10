/*
  Macro for merging of the calibration classes:
  marain.ivanov@cern.ch
  //
  //
  Example usage:
  1. Make a list of the calibration files
  2. Load libraries
  3. Merge - the output is stored in the current directory
  //
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  .L $ALICE_ROOT/TPC/macros/CalibFileMerger.C
  CalibFileMerger();
*/

void CalibFileMerger(const char* filename = "mergelist.txt") {


  TObjArray * mergeArray= 0;

  // Open the input stream
  ifstream in;
  in.open(filename);
  
  // Read the input list of files 
  TString objfile;
  Int_t counter=0;
  while(in.good()) {
    in >> objfile;
    if (!objfile.Contains("root")) continue; // protection    
    TFile currentFile(objfile.Data());
    TObjArray * carray =  (TObjArray*)currentFile.Get("TPCCalib");
    if (!carray) continue;
    if (!mergeArray) { mergeArray = carray; continue;}
    printf("Counter\t%d\tMerging file %s\n",counter,objfile.Data());
    for (Int_t icalib=0; icalib<mergeArray->GetEntries(); icalib++){
      AliTPCcalibBase *calib = ( AliTPCcalibBase *)mergeArray->At(icalib);
      //disable calib align
      //if (calib->InheritsFrom("AliTPCcalibAlign")) continue;
      if (calib) calib->Merge(carray);
    }
    AliSysInfo::AddStamp(objfile->Data(),counter,counter,counter);
    //currentFile.Close();
    carray->SetOwner(kTRUE);
    carray->Delete();
    delete carray;
    counter++;
    //if (counter>10) break;
  }

  in.close();
  TFile * output = new TFile("CalibObjects.root", "RECREATE");
  mergeArray->Write("TPCCalib",TObject::kSingleKey);
  output->Close();

  //

}

/*
//
TFile * output = new TFile("CalibObjects0.root", "RECREATE");
mergeArray->At(0)->Write();
output->Close();

TFile * output = new TFile("CalibObjects1.root", "RECREATE");
mergeArray->At(1)->Write();
output->Close();

TFile * output = new TFile("CalibObjects2.root", "RECREATE");
mergeArray->At(2)->Write();
output->Close();

TFile * output = new TFile("CalibObjects3.root", "RECREATE");
mergeArray->At(3)->Write();
output->Close();

TFile * output = new TFile("CalibObjects4.root", "RECREATE");
mergeArray->At(4)->Write();
output->Close();

TFile * output = new TFile("CalibObjects5.root", "RECREATE");
mergeArray->At(5)->Write();
output->Close();

TFile * output = new TFile("CalibObjects6.root", "RECREATE");
mergeArray->At(6)->Write();
output->Close();
*/


