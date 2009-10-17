
//
// Author: G. Balbastre
// Modified by K. Read
//
// Usage: 
// After previously using mymerger.C to merge each PYHIA hard bin,
// do the following to download one file per pt bin:
// aliensh; cd <production name>/output/merged; cp histoss* file:       
// Then, adjust binlist below to reflect which histograms are non-empty.
// Then, locally:
// root -b -x MergeFileInBins.C
//
// Note: This macro overwrites the downloaded histograms with the same names.
// It may be appropriate to make a backup of the histograms before running
// this macro.
//

void MergeFileInBins()
{
  char name[128] ;
  TFile * mfile;
  TFile * sfile;
  TList* list;
  TList* newlist;
  TString binlist = "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15";

  if (binlist.Length()) {
    arr = binlist.Tokenize(" ");
    TObjString *os;
    TIter next(arr);
    while ((os=(TObjString*)next())) {
      printf("Rescaling histosscaled-merged%s.root\n",os->GetString().Data());
	
      //Rescale histograms
      sprintf(name,"histosscaled-merged%s.root",os->GetString().Data());
      mfile = new TFile(name,"read");
      list = (TList*) mfile->Get("histosscaled");
      mfile->Close();
	
      TObject * h ; 
      Int_t split =  ((TH1F*) list->FindObject("hCount"))->GetEntries();
      cout<<"scale with factor "<<split<< " histograms "<<list->GetEntries()<<endl;
      newlist = new TList();
      sprintf(name,"%s2",list->GetName());
      newlist->SetName(name);

      for(Int_t iter = 0; iter < list->GetEntries(); iter++){
        h = list->At(iter);
        if(h && (h->GetName()!="hCount")){
	  if ( !strncmp(h->ClassName(),"TH",2) ) {
	      char name[128] ; 
	      sprintf(name, "%s", h->GetName()) ; 
	      //cout<<iter<<" histo scaled : "<<name<<endl;
	      TH1 * hout = dynamic_cast<TH1*> (h->Clone(name)) ; 
              //if(fSumw2) hout->Sumw2();
	      hout->Scale(1./split) ;  
              newlist->Add(hout) ; 
	  }
        } 
      }

      sprintf(name,"histosscaled-merged%s.root",os->GetString().Data());
      sfile =  new TFile(name,"recreate");
      newlist->Write();
      sfile->Close();	
      cout<<name<<"  has been recreated "<<endl;
    }//while
  }


  //Merge all histos

  cout<<"Merge all bins "<<endl;

  TFileMerger m;
  sprintf(name,"TOTALhistosscaled.root");

  if (binlist.Length()) {
    m.OutputFile(name);
    arr = binlist.Tokenize(" ");
    TObjString *os;
    TIter next(arr);
    while ((os=(TObjString*)next())) {
      printf("Adding bin histosscaled-merged%s.root\n",os->GetString().Data());
      sprintf(name,"histosscaled-merged%s.root",os->GetString().Data());
      m.AddFile(name);
    }
    m.Merge();
  }

}
