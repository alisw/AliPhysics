void RichBatch(const Int_t iNevents,const Bool_t isDebug,const char *sConfigFileName)
{
  if(isDebug)
    gAlice->SetDebug(1);

  Info("my/AliceBatch.C","%i event(s) requested, debug %i",iNevents,isDebug);  
  TStopwatch sw;TDatime time;	

  gAlice->Run(iNevents,sConfigFileName);
   
  cout<<"\nInfo in <my/AliceBatch.C>: Start time: ";time.Print();
  cout<<"Info in <my/AliceBatch.C>: Stop  time: ";time.Set();  time.Print();
  cout<<"Info in <my/AliceBatch.C>: Time  used: ";sw.Print();
  gSystem->Exec("touch ZZZfinishedZZZ");
}
