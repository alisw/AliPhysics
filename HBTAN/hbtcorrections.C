void hbtcorrections(Int_t first = -1,Int_t last = -1)
{  
  
  const char* basedir=".";
  const char* serie="";
  const char* field = "";

  AliHBTReader* reader = new AliHBTReaderInternal("ESD.old.root");
  TObjArray* dirs=0;
  if ( (first >= 0) && (last>=0) && ( (last-first)>=0 ) )
   {//read from many dirs dirs
     char buff[50];
     dirs = new TObjArray(last-first+1);
     for (Int_t i = first; i<=last; i++)
      {
        sprintf(buff,"%s/%s/%s/%d",basedir,field,serie,i);
        TObjString *odir= new TObjString(buff);
        dirs->Add(odir);
      }
    }
  reader->SetDirs(dirs);

  AliHBTParticleCut* readerpartcut= new AliHBTParticleCut();
  readerpartcut->SetPtRange(0.0,1.0);
  readerpartcut->SetPID(kPiPlus);
  reader->AddParticleCut(readerpartcut);//read this particle type with this cut
  
  AliHBTAnalysis* analysis = new AliHBTAnalysis();
  analysis->SetReader(reader);
  analysis->SetDisplayInfo(100000);
  AliHBTPairCut *paircut = new AliHBTPairCut();
  paircut->SetQInvRange(0.0,0.15);
  analysis->SetGlobalPairCut(paircut);
  


  AliHBTCorrectQInvCorrelFctn* correctqinvCF = new AliHBTCorrectQInvCorrelFctn();
  
  analysis->AddTrackFunction(correctqinvCF);

  TFile* outf = TFile::Open("hbtcorrected.root","recreate");
  if (outf == 0x0)
   {
     cout<<"\n\nERROR: can not open file"<<endl;
     return;
   }
  Int_t iter = 0;
  TString dirname;
  TString ietration("Iteration");
  correctqinvCF->SetInitialValues(0.8,12.0);
  correctqinvCF->SetRadiusConvergenceTreshold(0.005);
  correctqinvCF->SetLambdaConvergenceTreshold(0.01);
  while (1)
   {
    iter++;
    dirname = ietration+iter;
    TDirectory* dir = outf->mkdir(dirname,"results after "+dirname);
    outf->Write();
    if (dir) dir->cd();
    else 
     {
       cout<<"\n\nERROR: can not make an directory in file named"<<dirname<<endl;
       return;
     }

    analysis->Process("Tracks");

/*****************************************/
    dir->cd();
    correctqinvCF->WriteAll();
/*****************************************/
    if (correctqinvCF->IsConverged()) break;
    else
     {
       correctqinvCF->SetInitialValues(correctqinvCF->GetFittedLambda(),correctqinvCF->GetFittedRadius());
     }
   }
  outf->cd();
  analysis->Write();
  delete outf;
}
