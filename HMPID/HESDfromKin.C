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

  InitGRP();
//  AliMagF *magFieldMap = (AliMagF*)TGeoGlobalMagField::Instance()->GetField();

//  AliHMPIDTracker::SetFieldMap(gAL->GetAliRun()->Field(),kTRUE);
  AliHMPIDTracker pTracker;
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
      Printf("track n. %i",i);
      AliESDtrack trk(pTrack); 
      Float_t xPc,yPc,xRa,yRa,thRa,phRa;
      Int_t iCh=pTracker.IntTrkCha(&trk,xPc,yPc,xRa,yRa,thRa,phRa);         //get chamber intersected by this track 
      if(iCh<0) {
        trk.SetHMPIDtrk(0,0,0,0);                                                                //no intersection found
        trk.SetHMPIDcluIdx   (99,99999);                                                         //chamber not found, mip not yet considered
        trk.SetHMPIDsignal(AliHMPIDRecon::kNotPerformed);                                        //ring reconstruction not yet performed
        continue;                                                           //no intersection at all, go after next track
      }
      nTrkHMPID++;
      trk.SetHMPIDcluIdx   (iCh,99999);                                                          //chamber not found, mip not yet considered
      
      if(phRa<0) phRa += TMath::TwoPi(); // to be verified
      
      trk.SetHMPIDtrk(xPc,yPc,thRa,phRa);                                                        //store initial infos
      pEsd->AddTrack(&trk);
    
      Int_t status;
      if(!htaCheck) status = pTracker.Recon         (pEsd,pH->CluLst(),pNmean,pQthre);
      else          status = pTracker.ReconHiddenTrk(pEsd,pH->CluLst(),pNmean,pQthre);

//      Printf("status %i",status);
      if(status !=0) continue;

    
    }// track loop
    
    if(!(iEvt%50)) Printf("Number of events processed: %i with tracks %i in HMPID",iEvt,nTrkHMPID);
//      Printf("Number of events processed: %i with tracks %i in HMPID",iEvt,nTrkHMPID);
    
    gEsdTr->Fill();
    pEsd->Reset();
  }// event loop
  Printf("Events processed %i",iEvt);
  if(htaCheck) {
    fout->Write();
    fout->Close();
    delete fout;
  }
  gAL->UnloadHeader();  gAL->UnloadKinematics();
}//SimEsd()
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t OpenCalib()
{
  AliCDBManager* pCDB = AliCDBManager::Instance();
  pCDB->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  pCDB->SetRun(0);
  AliCDBEntry *pQthreEnt=pCDB->Get("HMPID/Calib/Qthre",0);
  AliCDBEntry *pNmeanEnt=pCDB->Get("HMPID/Calib/Nmean",0);
  
  if(!pQthreEnt || !pNmeanEnt) return kFALSE;
  
  pNmean=(TObjArray*)pNmeanEnt->GetObject(); 
  pQthre=(TObjArray*)pQthreEnt->GetObject(); 

  if(!pQthre || !pNmean) return kFALSE;  
  return kTRUE;
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bool_t InitGRP() {
  //------------------------------------
  // Initialization of the GRP entry 
  //------------------------------------
  
  AliGRPObject *fGRPData;
  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");

  if (entry) {

    TMap* m = dynamic_cast<TMap*>(entry->GetObject());  // old GRP entry

    if (!m) {
       Printf("Found a TMap in GRP/GRP/Data, converting it into an AliGRPObject");
       m->Print();
       fGRPData = new AliGRPObject();
       fGRPData->ReadValuesFromMap(m);
    }

    else {
       Printf("Found an AliGRPObject in GRP/GRP/Data, reading it");
       fGRPData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
       entry->SetOwner(0);
    }

    AliCDBManager::Instance()->UnloadFromCache("GRP/GRP/Data");
  }

  if (!fGRPData) {
     Printf("No GRP entry found in OCDB!");
     return kFALSE;
  }

  TString lhcState = fGRPData->GetLHCState();
  if (lhcState==AliGRPObject::GetInvalidString()) {
    Printf("GRP/GRP/Data entry:  missing value for the LHC state ! Using UNKNOWN");
    lhcState = "UNKNOWN";
  }

  TString beamType = fGRPData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    Printf("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  Float_t beamEnergy = fGRPData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    Printf("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }
  // energy is provided in MeV*120
  beamEnergy /= 120E3;

  TString runType = fGRPData->GetRunType();
  if (runType==AliGRPObject::GetInvalidString()) {
    Printf("GRP/GRP/Data entry:  missing value for the run type ! Using UNKNOWN");
    runType = "UNKNOWN";
  }

  Int_t activeDetectors = fGRPData->GetDetectorMask();
  if (activeDetectors==AliGRPObject::GetInvalidUInt()) {
    Printf("GRP/GRP/Data entry:  missing value for the detector mask ! Using 1074790399");
    activeDetectors = 1074790399;
  }

  fRunInfo = new AliRunInfo(lhcState, beamType, beamEnergy, runType, activeDetectors);
  printf("qqqqqqqqqqqqqqqqqqqqqqq %s %s %f %s %d\n", lhcState.Data(), beamType.Data(), beamEnergy, runType.Data(), activeDetectors);
  fRunInfo->Dump();

  //*** Dealing with the magnetic field map
  if ( TGeoGlobalMagField::Instance()->IsLocked() ) {Printf("Running with the externally locked B field !");}
  else {
    // Construct the field map out of the information retrieved from GRP.
    Bool_t ok = kTRUE;
    // L3
    Float_t l3Current = fGRPData->GetL3Current((AliGRPObject::Stats)0);
    if (l3Current == AliGRPObject::GetInvalidFloat()) {
      Prtinf("GRP/GRP/Data entry:  missing value for the L3 current !");
      ok = kFALSE;
    }
    
    Char_t l3Polarity = fGRPData->GetL3Polarity();
    if (l3Polarity == AliGRPObject::GetInvalidChar()) {
      Printf("GRP/GRP/Data entry:  missing value for the L3 polarity !");
      ok = kFALSE;
    }

    // Dipole
    Float_t diCurrent = fGRPData->GetDipoleCurrent((AliGRPObject::Stats)0);
    if (diCurrent == AliGRPObject::GetInvalidFloat()) {
      Printf("GRP/GRP/Data entry:  missing value for the dipole current !");
      ok = kFALSE;
    }

    Char_t diPolarity = fGRPData->GetDipolePolarity();
    if (diPolarity == AliGRPObject::GetInvalidChar()) {
      Printf("GRP/GRP/Data entry:  missing value for the dipole polarity !");
      ok = kFALSE;
    }

    if (ok) { 
      if ( !SetFieldMap(l3Current, diCurrent, l3Polarity ? -1:1, diPolarity ? -1:1) )
	AliFatal("Failed to creat a B field map ! Exiting...");
      Printf("Running with the B field constructed out of GRP !");
    }
    else AliFatal("B field is neither set nor constructed from GRP ! Exitig...");
    
  }
  
  return kTRUE;
} 

//_____________________________________________________________________________
//_____________________________________________________________________________
  Bool_t SetFieldMap(Float_t l3Cur=30000., Float_t diCur=6000., 
		     Float_t l3Pol=1., Float_t diPol=1., Float_t beamenergy=7000., 
		     const Char_t* beamtype="pp",  
		     const Char_t* path="$(ALICE_ROOT)/data/maps/mfchebKGI_sym.root")
{
  //------------------------------------------------
  // The magnetic field map, defined externally...
  // L3 current 30000 A  -> 0.5 T
  // L3 current 12000 A  -> 0.2 T
  // dipole current 6000 A
  // The polarities must be the same
  //------------------------------------------------
  const Float_t l3NominalCurrent1=30000.; // (A)
  const Float_t l3NominalCurrent2=12000.; // (A)
  const Float_t diNominalCurrent =6000. ; // (A)

  const Float_t tolerance=0.03; // relative current tolerance
  const Float_t zero=77.;       // "zero" current (A)
  //
  TString s=(l3Pol < 0) ? "L3: -" : "L3: +";
  //
  AliMagF::BMap_t map = AliMagF::k5kG;
  //
  double fcL3,fcDip;
  //
  l3Cur = TMath::Abs(l3Cur);
  if (TMath::Abs(l3Cur-l3NominalCurrent1)/l3NominalCurrent1 < tolerance) {
    fcL3 = l3Cur/l3NominalCurrent1;
    map  = AliMagF::k5kG;
    s   += "0.5 T;  ";
  } else if (TMath::Abs(l3Cur-l3NominalCurrent2)/l3NominalCurrent2 < tolerance) {
    fcL3 = l3Cur/l3NominalCurrent2;
    map  = AliMagF::k2kG;
    s   += "0.2 T;  ";
  } else if (l3Cur <= zero) {
    fcL3 = 0;
    map  = AliMagF::k5kGUniform;
    s   += "0.0 T;  ";
    fUniformField=kTRUE;        // track with the uniform (zero) B field
  } else {
    AliError(Form("Wrong L3 current (%f A)!",l3Cur));
    return kFALSE;
  }
  //
  diCur = TMath::Abs(diCur);
  if (TMath::Abs(diCur-diNominalCurrent)/diNominalCurrent < tolerance) {
    // 3% current tolerance...
    fcDip = diCur/diNominalCurrent;
    s    += "Dipole ON";
  } else if (diCur <= zero) { // some small current..
    fcDip = 0.;
    s    += "Dipole OFF";
  } else {
    AliError(Form("Wrong dipole current (%f A)!",diCur));
    return kFALSE;
  }
  //
  if (l3Pol!=diPol && (map==AliMagF::k5kG || map==AliMagF::k2kG) && fcDip!=0) {
    AliError("L3 and Dipole polarities must be the same");
    return kFALSE;
  }
  //
  if (l3Pol<0) fcL3  = -fcL3;
  if (diPol<0) fcDip = -fcDip;
  //
  AliMagF::BeamType_t btype = AliMagF::kNoBeamField;
  TString btypestr = beamtype;
  btypestr.ToLower();
  TPRegexp protonBeam("(proton|p)\\s*-?\\s*\\1");
  TPRegexp ionBeam("(lead|pb|ion|a)\\s*-?\\s*\\1");
  if (btypestr.Contains(ionBeam)) btype = AliMagF::kBeamTypeAA;
  else if (btypestr.Contains(protonBeam)) btype = AliMagF::kBeamTypepp;
  else {
    Printf(Form("Cannot determine the beam type from %s, assume no LHC magnet field",beamtype));
  }
  Printf("------------------------------");
  Printf(" Summary for B: %s",s.Data());
  Printf("------------------------------");
  AliMagF* fld = new AliMagF("MagneticFieldMap", s.Data(), 2, fcL3, fcDip, 10., map, path, 
			     btype,beamenergy);
  TGeoGlobalMagField::Instance()->SetField( fld );
  TGeoGlobalMagField::Instance()->Lock();
  //
  return kTRUE;
}
