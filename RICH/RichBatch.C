void RichBatch(const Int_t iNevents,const Bool_t isDebug,const char *sConfigFileName)
{
  if(isDebug) gAlice->SetDebug(1);

  Info("my/RichBatch.C","%i event(s) requested, debug %i,config file %s",iNevents,isDebug,sConfigFileName);  
  TStopwatch sw;TDatime time;	

  AliSimulation a;  a.Run(iNevents);
   
  cout<<"\nInfo in <my/RichBatch.C>: Start time: ";time.Print();
    cout<<"Info in <my/RichBatch.C>: Stop  time: ";time.Set();  time.Print();
    cout<<"Info in <my/RichBatch.C>: Time  used: ";sw.Print();
  gSystem->Exec("touch ZZZ______finished_______ZZZ");
}
