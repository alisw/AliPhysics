AliRunLoader *gAL=0; 
Int_t gEvt=0; Int_t gMaxEvt=0;
TObjArray *pNmean;
TTree *gEsdTr;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HESDfromKin()
{//simulate ESD from kinematics

  if(gSystem->IsFileInIncludePath("galice.root")){// tries to open session
    if(gAlice) delete gAlice;                                               //in case we execute this in aliroot delete default AliRun object 
    gAL=AliRunLoader::Open();                                                                    //try to open galice.root from current dir 
    gAL->LoadgAlice();                                                                           //take new AliRun object from galice.root   
    if(gAL->LoadHeader()) return;
    if(gAL->LoadKinematics()) return;

    AliLoader *pHL=gAL->GetDetectorLoader("HMPID");
    pHL->LoadRecPoints();
    AliESD *pEsd = new AliESD();   
    TFile *pEsdFl=TFile::Open("AliESDs.root","recreate"); 
    gEsdTr=new TTree("esdTree","Sim ESD from kinematics"); 
    gEsdTr->Branch("ESD", &pEsd);
         
  }  else return;  

  OpenCalib();
    
  SimEsd(pHL,pEsd);
//  SimEsdHidden(pHL,pEsd);
  
  pEsdFl->Write();pEsdFl->Close();        
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimEsd(AliLoader *pHL,AliESD *pEsd)
{
  Printf("---------okokokok");
  AliHMPIDTracker::SetFieldMap(gAL->GetAliRun()->Field(),kTRUE);
  AliHMPID *pH=(AliHMPID*)gAL->GetAliRun()->GetDetector("HMPID");
  Int_t mtid=-1;
  Int_t iNevt=gAL->GetNumberOfEvents();
  Printf("Number of events to process: %i",iNevt);
  for(Int_t iEvt=0;iEvt<iNevt;iEvt++){//events loop
    if(!(iEvt%50)) Printf("Events processed %i",iEvt);
    gAL->GetEvent(iEvt);    
    pHL->TreeR()->GetEntry(0);
    AliStack *pStack=gAL->Stack();
    for(Int_t i=0;i<pStack->GetNtrack();i++){
      TParticle *pTrack=pStack->Particle(i); 
      mtid=pTrack->GetFirstMother();
      if(mtid>=0) continue; // only primaries
      AliESDtrack trk(pTrack);
      pEsd->AddTrack(&trk);
      AliHMPIDTracker::Recon(pEsd,pH->CluLst(),pNmean);
    }// track loop
    pEsd->SetMagneticField(AliHMPIDTracker::GetBz());
    gEsdTr->Fill();
    pEsd->Reset();
  }// event loop
  Printf("Events processed %i",iEvt);
  gAL->UnloadHeader();  gAL->UnloadKinematics();
}//EsdFromStack()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimEsdHidden(AliLoader *pHL,AliESD *pEsd)
{  
  AliHMPIDTracker::SetFieldMap(gAL->GetAliRun()->Field(),kTRUE);
  AliHMPID *pH=(AliHMPID*)gAL->GetAliRun()->GetDetector("HMPID");
  Int_t mtid=-1;
  Int_t iNevt=gAL->GetNumberOfEvents();
  Printf("Number of events to process: %i",iNevt);
  for(Int_t iEvt=0;iEvt<iNevt;iEvt++){//events loop
    if(!(iEvt%50)) Printf("Events processed %i",iEvt);
    gAL->GetEvent(iEvt);
    pHL->TreeR()->GetEntry(0);
    AliESDtrack trk;
// Hidden Track
    for(int iCh=AliHMPIDDigit::kMinCh;iCh<=AliHMPIDDigit::kMaxCh;iCh++){
      AliESDtrack trk;
       if(AliHMPIDTracker::ReconHiddenTrk(iCh,&trk,pH->CluLst(),pNmean);) continue;
       pEsd->AddTrack(&trk);
       Printf(" theta Cerenkov reconstructed %f",trk.GetHMPIDsignal());
    }
//    AliHMPIDTracker::Recon(pEsd,pH->CluLst(),pNmean);
//      
    pEsd->SetMagneticField(AliHMPIDTracker::GetBz());
    gEsdTr->Fill();
    pEsd->Reset();
  }// event loop
  Printf("Events processed %i",iEvt);
  gAL->UnloadHeader();  gAL->UnloadKinematics();
}//EsdFromStack()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void OpenCalib()
{
  AliCDBManager* pCDB = AliCDBManager::Instance();
  pCDB->SetDefaultStorage("local://$HOME");
  AliCDBEntry *pQthreEnt=pCDB->Get("HMPID/Calib/Qthre",0);
  AliCDBEntry *pNmeanEnt=pCDB->Get("HMPID/Calib/Nmean",0);
  
  if(!pQthreEnt || ! pNmeanEnt) return;
  
  pNmean=(TObjArray*)pNmeanEnt->GetObject(); 
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
