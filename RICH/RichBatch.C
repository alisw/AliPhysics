void RichBatch(const Int_t iNevents,const Bool_t isDebug,const char *sConfigFileName)
{
  gSystem->Exec("rm -rf *.root hlt hough ZZZ*");
  if(isDebug)   AliLog::SetGlobalDebugLevel(AliLog::kDebug);

  TStopwatch sw;TDatime time;	
  
  
  AliSimulation     *pSim=new AliSimulation;     pSim->Run(iNevents); delete pSim;
  AliReconstruction *pRec=new AliReconstruction; 
  pRec->SetRunLocalReconstruction("ITS TPC TRD TOF RICH");
  pRec->SetFillESD("ITS TPC TRD TOF RICH");
  pRec->Run();         delete pRec;
  
  //pRec->SetRunLocalReconstruction("RICH"); pRec->SetFillESD("RICH");    
  
  cout<<"\n!!!!!!!!!!!!Info in <my/RichBatch.C>: Start creating Control Plots \n";
  
  AliRunLoader* pAL = AliRunLoader::Open();  pAL->LoadgAlice();
  AliRICH *pRICH=(AliRICH*)pAL->GetAliRun()->GetDetector("RICH");
  if(pRICH) pRICH->ControlPlots();
 
  cout<<"\n!!!!!!!!!!!!Info in <my/RichBatch.C>: Start time: ";time.Print();
    cout<<"!!!!!!!!!!!!Info in <my/RichBatch.C>: Stop  time: ";time.Set();  time.Print();
    cout<<"!!!!!!!!!!!!Info in <my/RichBatch.C>: Time  used: ";sw.Print();
  gSystem->Exec("touch ZZZ______finished_______ZZZ");
//  gSystem->Exec("playwave .kde/my/end.wav");
}
