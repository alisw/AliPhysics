void RichBatch(const Int_t iNevents,const Bool_t isDebug,const char *sConfigFileName)
{
  if(isDebug) gAlice->SetDebug(1);

  Info("my/RichBatch.C","%i event(s) requested, debug %i,config file %s",iNevents,isDebug,sConfigFileName);  
  TStopwatch sw;TDatime time;	

  AliSimulation     *pSim=new AliSimulation;     pSim->Run(iNevents); delete pSim;
  AliReconstruction *pRec=new AliReconstruction; pRec->SetRunLocalReconstruction("RICH"); pRec->SetFillESD("RICH");  pRec->Run();         delete pRec;
  
  AliRunLoader* pAL = AliRunLoader::Open();
  pAL->LoadgAlice();
  AliRICH *pRICH=(AliRICH*)pAL->GetAliRun()->GetDetector("RICH");
  if(pRICH) pRICH->ControlPlots();
 
  cout<<"\nInfo in <my/RichBatch.C>: Start time: ";time.Print();
    cout<<"Info in <my/RichBatch.C>: Stop  time: ";time.Set();  time.Print();
    cout<<"Info in <my/RichBatch.C>: Time  used: ";sw.Print();
  gSystem->Exec("touch ZZZ______finished_______ZZZ");
}
