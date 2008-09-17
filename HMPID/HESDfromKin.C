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
//  if(!htaCheck) SimEsd(pHL,pEsd); else SimEsdHidden(pHL,pEsd);
  SimEsd(pHL,pEsd,htaCheck);
  
  pEsdFl->cd();
  pEsdFl->Write();pEsdFl->Close();        
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimEsd(AliLoader *pHL,AliESDEvent *pEsd,Bool_t htaCheck)
{
  if(htaCheck) {
    TFile *fout = new TFile("HTA.root","recreate");
    TH1F *hdC   = new TH1F("dC"  ,";delta Cerenkov (rad)",100,-0.2,0.2);
    TH1F *hCer  = new TH1F("Cer" ,"Theta Cerenkov (rad)",250,0.,0.75);
    TH2F *htvsp = new TH2F("tvsp",";momentum (GeV/c);theta Cerenkov (rad)",100,0.,5.,1000,0.,0.75);
    TH1F *hdth  = new TH1F("dth" ,";Delta theta Trk (mrad)",100,-250,250);
    TH1F *hdph  = new TH1F("dph" ,";Delta phi Trk (mrad)",100,-500,500);
    Double_t rd=TMath::RadToDeg();
    Printf("----------------------------------------------");
    Printf("| SimHTA:Utility to embed ESD from kinematics|");
    Printf("|     with  Hidden Track Algorithm (HTA)     |");
    Printf("----------------------------------------------");
  } else {
    Printf("-----------------------------------------------");
    Printf("| SimESD: Utility to embed ESD from kinematics|");
    Printf("-----------------------------------------------");
}
  AliHMPIDTracker::SetFieldMap(gAL->GetAliRun()->Field(),kTRUE);
  AliHMPID *pH=(AliHMPID*)gAL->GetAliRun()->GetDetector("HMPID");
  Int_t iNevt=gAL->GetNumberOfEvents();
  pEsd->SetMagneticField(AliHMPIDTracker::GetBz());
  for(Int_t iEvt=0;iEvt<iNevt;iEvt++){//events loop
    gAL->GetEvent(iEvt);    
    pHL->TreeR()->GetEntry(0);
    AliStack *pStack=gAL->Stack();
    Int_t nTrkHMPID=0;
    for(Int_t i=0;i<pStack->GetNtrack();i++){
      if(!pStack->IsPhysicalPrimary(i)) continue;
      TParticle *pTrack=pStack->Particle(i); 
      if(pTrack->GetPDG()->Charge()==0) continue;
      AliESDtrack trk(pTrack); 
      Float_t xPc,yPc,xRa,yRa,thRa,phRa;
      Int_t iCh=AliHMPIDTracker::IntTrkCha(&trk,xPc,yPc,xRa,yRa,thRa,phRa);         //get chamber intersected by this track 
      if(iCh<0) {
        trk.SetHMPIDtrk(0,0,0,0);                                                                //no intersection found
        trk.SetHMPIDcluIdx   (99,99999);                                                         //chamber not found, mip not yet considered
        trk.SetHMPIDsignal(AliHMPIDRecon::kNotPerformed);                                        //ring reconstruction not yet performed
        continue;                                                           //no intersection at all, go after next track
      }
      nTrkHMPID++;
      trk.SetHMPIDcluIdx   (iCh,99999);                                                          //chamber not found, mip not yet considered
      trk.SetHMPIDtrk(xRa,yRa,thRa,phRa);                                                        //store initial infos
      pEsd->AddTrack(&trk);
    }// track loop
    
    if(!(iEvt%50)) Printf("Number of events processed: %i with tracks %i in HMPID",iEvt,nTrkHMPID);
    if(!htaCheck) AliHMPIDTracker::Recon(pEsd,pH->CluLst(),pNmean,pQthre);
// Perform PID
    
    
    for(Int_t iTrk=0;iTrk<pEsd->GetNumberOfTracks();iTrk++){                                       //ESD tracks loop
      
      AliESDtrack *pTrk = pEsd->GetTrack(iTrk);                                                    //get reconstructed track    
      
      Float_t radX,radY,thetaTrk,phiTrk;
      pTrk->GetHMPIDtrk(radX,radY,thetaTrk,phiTrk);
      
      if(htaCheck) {
          Int_t iCh = pTrk->GetHMPIDcluIdx();                                                                  //chamber not found, mip not yet considered
          iCh/=1000000;
        Printf("simulated track theta %f phi %f",thetaTrk*rd,phiTrk*rd);
        TObjArray *pClus = pH->CluLst();

        if(AliHMPIDTracker::ReconHiddenTrk(pEsd,pClus,pNmean,pQthre)!=0) continue;

        Double_t thetaCerSim = TMath::ACos(pTrack->Energy()/(1.292*pTrack->P()));
        Printf("-----------------------------------------------------------");
        Printf(" theta Cerenkov simulated     %f",thetaCerSim);
        Printf(" theta Cerenkov reconstructed %f",pTrk->GetHMPIDsignal());
        Float_t thetaFit,phiFit;
        pTrk->GetHMPIDtrk(radX,radY,thetaFit,phiFit);
//        Printf("reconstr. track theta %f phi %f",thetaFit*rd,phiFit*rd);

        // fill useful histos
        hdC->Fill(pTrk->GetHMPIDsignal()-thetaCerSim);
        hCer->Fill(pTrk->GetHMPIDsignal());
        htvsp->Fill(pTrk->P(),pTrk->GetHMPIDsignal());
        hdth->Fill((thetaFit-thetaTrk)*1000);
        hdph->Fill((phiFit-phiTrk)*1000);
      }
      
      AliHMPIDPid pID;
      Double_t prob[5];
      pID.FindPid(pTrk,5,prob);
      pTrk->SetHMPIDpid(prob);
//      Printf(" Prob e- %6.2f mu %6.2f pi %6.2f k %6.2f p %6.2f",prob[0]*100,prob[1]*100,prob[2]*100,prob[3]*100,prob[4]*100);
    }
    
    gEsdTr->Fill();
    pEsd->Reset();
    
  }// event loop
  Printf("Events processed %i",iEvt);
  gAL->UnloadHeader();  gAL->UnloadKinematics();
}//Esd()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void SimEsdHidden(AliLoader *pHL,AliESDEvent *pEsd)
{
  TFile *fout = new TFile("HTA.root","recreate");
  TH1F *hdC   = new TH1F("dC"  ,";delta Cerenkov (rad)",100,-0.2,0.2);
  TH1F *hCer  = new TH1F("Cer" ,"Theta Cerenkov (rad)",250,0.,0.75);
  TH2F *htvsp = new TH2F("tvsp",";momentum (GeV/c);theta Cerenkov (rad)",100,0.,5.,1000,0.,0.75);
  TH1F *hdth  = new TH1F("dth" ,";Delta theta Trk (mrad)",100,-250,250);
  TH1F *hdph  = new TH1F("dph" ,";Delta phi Trk (mrad)",100,-500,500);
  Double_t rd=TMath::RadToDeg();
  Printf("----------------------------------------------");
  Printf("| SimHTA:Utility to embed ESD from kinematics|");
  Printf("|     with  Hidden Track Algorithm (HTA)     |");
  Printf("----------------------------------------------");
  AliHMPIDTracker::SetFieldMap(gAL->GetAliRun()->Field(),kTRUE);
  AliHMPID *pH=(AliHMPID*)gAL->GetAliRun()->GetDetector("HMPID");
  Int_t iNevt=gAL->GetNumberOfEvents();
  pEsd->SetMagneticField(AliHMPIDTracker::GetBz());
  
  Printf("Number of events to process: %i",iNevt);
  
  for(Int_t iEvt=0;iEvt<iNevt;iEvt++){//events loop
    if(!(iEvt%1)) Printf("Events processed %i",iEvt);
    gAL->GetEvent(iEvt);    
    pHL->TreeR()->GetEntry(0);
    AliStack *pStack=gAL->Stack();
    
    for(Int_t i=0;i<pStack->GetNtrack();i++){
      
      (!pStack->IsPhysicalPrimary(i)) continue;
      TParticle *pTrack=pStack->Particle(i); 
      
      if(pTrack->GetPDG()->Charge()==0) continue;
      //find the chamber that intersects HMPID
      pTrack->Print();
      AliESDtrack trk(pTrack);
      Float_t xPc,yPc,xRa,yRa,thRa,phRa;
      Int_t iCh=AliHMPIDTracker::IntTrkCha(&trk,xPc,yPc,xRa,yRa,thRa,phRa);         //get chamber intersected by this track
      if(iCh<0) continue;                                                           //no intersection at all, go after next track
      trk.SetHMPIDtrk(xRa,yRa,thRa,phRa);                                                        //store initial infos
      Float_t radX,radY,thetaTrk,phiTrk;
      trk.GetHMPIDtrk(radX,radY,thetaTrk,phiTrk);
//      Printf("simulated track theta %f phi %f",thetaTrk*rd,phiTrk*rd);
      TObjArray *pClus = pH->CluLst();
      
      if(AliHMPIDTracker::ReconHiddenTrk(iCh,0,&trk,(TClonesArray *)pClus->At(iCh),pNmean,pQthre)!=0) continue;
      
      pEsd->AddTrack(&trk);
      Double_t thetaCerSim = TMath::ACos(pTrack->Energy()/(1.292*pTrack->P()));
//      Printf(" theta Cerenkov simulated     %f",thetaCerSim);
//      Printf(" theta Cerenkov reconstructed %f",trk.GetHMPIDsignal());
      Float_t thetaFit,phiFit;
      trk.GetHMPIDtrk(radX,radY,thetaFit,phiFit);
//      Printf("reconstr. track theta %f phi %f",thetaFit*rd,phiFit*rd);
      
      // fill useful histos
      hdC->Fill(trk.GetHMPIDsignal()-thetaCerSim);
      hCer->Fill(trk.GetHMPIDsignal());
      htvsp->Fill(pTrack->P(),trk.GetHMPIDsignal());
      hdth->Fill((thetaFit-thetaTrk)*1000);
      hdph->Fill((phiFit-phiTrk)*1000);
    }// track loop
    pEsd->SetMagneticField(AliHMPIDTracker::GetBz());
    gEsdTr->Fill();
    pEsd->Reset();
  }// event loop
  Printf("Events processed %i",iEvt);
  fout->Write();
  fout->Close();
  delete fout;
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
