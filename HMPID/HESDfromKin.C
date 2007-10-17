AliRunLoader *gAL=0; 
Int_t gEvt=0; Int_t gMaxEvt=0;
TObjArray *pNmean,*pQthre;
TTree *gEsdTr;
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void HESDfromKin(const char *name="default")
{//simulate ESD from kinematics

  if(gSystem->IsFileInIncludePath("galice.root")){// tries to open session
    if(gAlice) delete gAlice;                                               //in case we execute this in aliroot delete default AliRun object 
    gAL=AliRunLoader::Open();                                                                    //try to open galice.root from current dir 
    gAL->LoadgAlice();                                                                           //take new AliRun object from galice.root   
    if(gAL->LoadHeader()) return;
    if(gAL->LoadKinematics()) return;

    AliLoader *pHL=gAL->GetDetectorLoader("HMPID");
    pHL->LoadRecPoints();
    AliESDEvent *pEsd = new AliESDEvent();   
    TFile *pEsdFl=TFile::Open("AliESDs.root","recreate"); 
    gEsdTr=new TTree("esdTree","Sim ESD from kinematics"); 
    pEsd->CreateStdContent();    pEsd->WriteToTree(gEsdTr);  //clm: new ESD write schema: see Task Force meeting 20th June, 2007
    gEsdTr->GetUserInfo()->Add(pEsd);                        //clm: TList has to be created for ReadFromTree method -- this was not needed by the old ESD
 
       
  }  else return;  

  if(!OpenCalib()) {Printf("Problems in OpenCalib!Bye.");return;}
    
  TString ttl=name;
  Bool_t htaCheck=ttl.Contains("HTA");
  if(!htaCheck) SimEsd(pHL,pEsd); else SimEsdHidden(pHL,pEsd);
  
  pEsdFl->cd();
  pEsdFl->Write();pEsdFl->Close();        
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimEsd(AliLoader *pHL,AliESDEvent *pEsd)
{
  Printf("-----------------------------------------------");
  Printf("| SimESD: Utility to embed ESD from kinematics|");
  Printf("-----------------------------------------------");
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
      AliHMPIDTracker::Recon(pEsd,pH->CluLst(),pNmean,pQthre);
    }// track loop
    pEsd->SetMagneticField(AliHMPIDTracker::GetBz());
    gEsdTr->Fill();
    pEsd->Reset();
  }// event loop
  Printf("Events processed %i",iEvt);
  gAL->UnloadHeader();  gAL->UnloadKinematics();
}//Esd()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimEsdHidden(AliLoader *pHL,AliESDEvent *pEsd)
{
  TH1F *hdC = new TH1F("dC","dC",100,-0.2,0.2);
  TH1F *hCer = new TH1F("Cer","Cer",250,0.,0.75);
  TH2F *htvsp = new TH2F("tvsp","tvsp",100,0.,5.,1000,0.,0.75);
  Double_t rd=TMath::RadToDeg();
  Printf("----------------------------------------------");
  Printf("| SimHTA:Utility to embed ESD from kinematics|");
  Printf("|     with  Hidden Track Algorithm (HTA)     |");
  Printf("----------------------------------------------");
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
      //find the chamber that intersects HMPID
      AliESDtrack trk(pTrack);
      Float_t xPc,yPc;
      Int_t iCh=AliHMPIDTracker::IntTrkCha(&trk,xPc,yPc);                           //get chamber intersected by this track 
      if(iCh<0) continue;                                                           //no intersection at all, go after next track
      Float_t radX,radY,thetaTrk,phiTrk;
      trk.GetHMPIDtrk(radX,radY,thetaTrk,phiTrk);
//      Printf("simulated track theta %f phi %f",thetaTrk*rd,phiTrk*rd);
      TObjArray *pClus = pH->CluLst();
      if(AliHMPIDTracker::ReconHiddenTrk(iCh,&trk,(TClonesArray *)pClus->At(iCh),pNmean,pQthre)!=0) continue;
      trk.GetHMPIDtrk(radX,radY,thetaTrk,phiTrk);
//      Printf("reconstr. track theta %f phi %f",thetaTrk*rd,phiTrk*rd);
      pEsd->AddTrack(&trk);
      Double_t thetaCerSim = TMath::ACos(pTrack->Energy()/(1.292*pTrack->P()));
//      Printf(" theta Cerenkov simulated     %f",thetaCerSim);
//      Printf(" theta Cerenkov reconstructed %f",trk.GetHMPIDsignal());
      hdC->Fill(trk.GetHMPIDsignal()-thetaCerSim);
      hCer->Fill(trk.GetHMPIDsignal());
      htvsp->Fill(pTrack->P(),trk.GetHMPIDsignal());
    }// track loop
    pEsd->SetMagneticField(AliHMPIDTracker::GetBz());
    gEsdTr->Fill();
    pEsd->Reset();
  }// event loop
  Printf("Events processed %i",iEvt);
  gAL->UnloadHeader();  gAL->UnloadKinematics();
}//EsdHidden()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t OpenCalib()
{
  AliCDBManager* pCDB = AliCDBManager::Instance();
  pCDB->SetDefaultStorage("local://$HOME");
  AliCDBEntry *pQthreEnt=pCDB->Get("HMPID/Calib/Qthre",0);
  AliCDBEntry *pNmeanEnt=pCDB->Get("HMPID/Calib/Nmean",0);
  
  if(!pQthreEnt || !pNmeanEnt) return kFALSE;
  
  pNmean=(TObjArray*)pNmeanEnt->GetObject(); 
  pQthre=(TObjArray*)pQthreEnt->GetObject(); 

  if(!pQthre || !pNmean) return kFALSE;  
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
