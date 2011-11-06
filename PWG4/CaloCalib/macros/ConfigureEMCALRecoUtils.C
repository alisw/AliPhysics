
void ConfigureEMCALRecoUtils(
                             AliEMCALRecoUtils* reco,
                             Bool_t  bMC   = kFALSE,
                             TGeoHMatrix* matrix[10],
                             TString path  = "",
                             Int_t   run   = 0, 
                             TString pass  = "pass2"
                             )
{  

  // Configure RecoUtils with OADB objects
  
  printf("**** Configure AliEMCALRecoUtils, LOAD AODB ***\n");
  printf("\t run %d, pass %s\n",run,pass.Data());
  
  // Exotic cells removal
  reco->SwitchOnRejectExoticCell() ;
  reco->SetExoticCellDiffTimeCut(10000); // Open  
  reco->SetExoticCellFractionCut(0.95);  // 1-Ecross/Ecell > 0.95 -> out
  reco->SetExoticCellMinAmplitudeCut(2); // 2 GeV  
  
  gSystem->Load("libOADB");
  
  // Instantiate EMCAL geometry for the first time
  
  AliEMCALGeometry*   geom = 0; 
  if   (run < 140000) geom = AliEMCALGeometry::GetInstance("EMCAL_FIRSTYEARV1");
  else                geom = AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");

  Int_t nSM = geom->GetNumberOfSuperModules();

  // Geometry settings and field
    
  TGeoManager::Import("geometry.root") ; //need file "geometry.root" in local dir!!!!
  
  // Put the automatic initialization like in tender somewhere in the clusterizer or reader!!!
  if(run > 140000 && run < 160000) 
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG)); // for (--)
  else if(run > 160000)             
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps",  1.,  1., AliMagF::k5kG)); // for (++)
  else if (run < 140000 && run > 138280)
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps",  1.,  1., AliMagF::k5kG)); // for (++)
  else if (run < 138280 && run > 122000)
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG)); // for (--)
  else if(run > 119100 && run < 122000)
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps",  1.,  1., AliMagF::k5kG)); // for (++)
  else if(run < 119100 && run > 11900)
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps",  0.,  0., AliMagF::k5kG)); // for (0)
  else if(run < 119000 && run > 115800)
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps",  1.,  1., AliMagF::k5kG)); // for (++)
  else if(run < 115800)
    TGeoGlobalMagField::Instance()->SetField(new AliMagF("Maps","Maps", -1., -1., AliMagF::k5kG)); // for (--)

  // Alignment matrices

  TString fileName="$ALICE_ROOT/OADB/EMCAL/EMCALlocal2master.root";
  if(path!="") fileName=path+"EMCALlocal2master.root";

  AliOADBContainer EMCALgeoCont("AliEMCALgeo");
  EMCALgeoCont.InitFromFile((char*)fileName.Data(),"AliEMCALgeo");
  TObjArray *mobj=(TObjArray*)EMCALgeoCont.GetObject(run,"EmcalMatrices");
  for (Int_t mod=0;mod<nSM;mod++)
    {
      matrix[mod] = (TGeoHMatrix*) mobj->At(mod);
      //matrix[mod]->Print();
    }
      
  reco->SetPositionAlgorithm(AliEMCALRecoUtils::kPosTowerGlobal);   
  
  
  //----------------------------------------------
  
  Bool_t bRecal  = kTRUE;
  Bool_t bBad    = kTRUE;
  if(pass == "pass3"){
    bRecal = kFALSE;
    bBad   = kFALSE;
  }
  
  if(bMC){
    bRecal = kFALSE;
  }
  
  //Recalibration factors
  if(bRecal){
    
    fileName="$ALICE_ROOT/OADB/EMCAL/EMCALRecalib.root";
    if(path!="") fileName=path+"EMCALRecalib.root";
    
    AliOADBContainer *contRF=new AliOADBContainer("");
    contRF->InitFromFile((char*)fileName.Data(),"AliEMCALRecalib");
    
    TObjArray *recal=(TObjArray*)contRF->GetObject(run); 
    if(recal){
      TObjArray *recalpass=(TObjArray*)recal->FindObject(pass);
      if(recalpass){
        TObjArray *recalib=(TObjArray*)recalpass->FindObject("Recalib");
        if(recalib){
          reco->SwitchOnRecalibration();
          printf("AliEMCALRecoUtils - RECALIBRATE \n");
          for (Int_t i=0; i<nSM; ++i) {
            TH2F *h = reco->GetEMCALChannelRecalibrationFactors(i);
            if (h)
              delete h;
            h = (TH2F*)recalib->FindObject(Form("EMCALRecalFactors_SM%d",i));
            if (!h) {
              AliError(Form("Could not load EMCALRecalFactors_SM%d",i));
              continue;
            }
            h->SetDirectory(0);
            reco->SetEMCALChannelRecalibrationFactors(i,h);
          } 
        }else printf("AliEMCALRecoUtils ---Do NOT recalibrate 1\n");
      }else printf("AliEMCALRecoUtils ---Do NOT recalibrate 2\n");
    }else printf("AliEMCALRecoUtils ---Do NOT recalibrate 3\n");
    
    //TFile * f = new TFile("RecalibrationFactors.root","read");
    //for(Int_t i =0; i< 10; i++)  reco->SetEMCALChannelRecalibrationFactors( i, (TH2F*) f->Get(Form("EMCALRecalFactors_SM%d",i)));									 
    // //  reco->SwitchOnTimeDepCorrection();
    // //  //char cmd[200] ; 
    // //  //sprintf(cmd, ".!tar xvfz CorrectionFiles.tgz") ; 
    // //  //gROOT->ProcessLine(cmd) ; 
    // //  	
    
  } else printf("AliEMCALRecoUtils ---Do NOT recalibrate\n");

  // Remove EMCAL hot channels 
  
  if(bBad){
    
    fileName="$ALICE_ROOT/OADB/EMCAL/EMCALBadChannels.root";
    if(path!="") fileName=path+"EMCALBadChannels.root";
    
    AliOADBContainer *contBC=new AliOADBContainer("");
    contBC->InitFromFile((char*)fileName.Data(),"AliEMCALBadChannels"); 
    TObjArray *arrayBC=(TObjArray*)contBC->GetObject(run);
    if(arrayBC){
      TObjArray *arrayBCpass=(TObjArray*)arrayBC->FindObject(pass);
      if(arrayBCpass){
        
        reco->SwitchOnBadChannelsRemoval();
        reco->SwitchOnDistToBadChannelRecalculation();
        printf("AliEMCALRecoUtils - REMOVE bad cells \n");

        for (Int_t i=0; i<nSM; ++i) {
          TH2I *hbm = reco->GetEMCALChannelStatusMap(i);
          if (hbm)
            delete hbm;
          hbm=(TH2I*)arrayBCpass->FindObject(Form("EMCALBadChannelMap_Mod%d",i));
          
          if (!hbm) {
            AliError(Form("Can not get EMCALBadChannelMap_Mod%d",i));
            continue;
          }
          
          hbm->SetDirectory(0);
          reco->SetEMCALChannelStatusMap(i,hbm);
        }
      } else printf("AliEMCALRecoUtils ---Do NOT remove bad channels 1\n");
    }  else printf("AliEMCALRecoUtils ---Do NOT remove bad channels 2\n");
  } else printf("AliEMCALRecoUtils ---Do NOT remove bad channels 3 \n");

  /*  
    Int_t iCol = -1, iRow = -1, iSM =-1, iMod = -1,iIphi =-1,iIeta = -1;
    Int_t badAbsID[]={74, 103};
    
    for(Int_t i=0;i < sizeof(badAbsID)/sizeof(Int_t); i++){
      geom->GetCellIndex(badAbsID[i],iSM,iMod,iIphi,iIeta); 
      // Gives SuperModule and Tower numbers
      geom->GetCellPhiEtaIndexInSModule(iSM,iMod,
					iIphi, iIeta,iRow,iCol);
      //printf("bad ID %d, col %d, row %d, sm %d\n",badAbsID[i],iCol,iRow,iSM);
      reco->SetEMCALChannelStatus(iSM , iCol, iRow,1);
    }

    }
   */
 
  
  reco->SwitchOnTimeRecalibration();
  
  //Waiting for OADB, meanwhile

  TFile * ftime = 0;
  //TString path ="./"; 
  TString path = "alien:///alice/cern.ch/user/g/germain/RecalDB/Time";
  //TGrid::Connect("alien://");
  
  if     (run > 140000 && run < 146500 )
    ftime = TFile::Open(Form("%s/RefLHC11apass1-7TeV.root",path.Data()));
  else if(run > 146500 && run <= 146860 )
    ftime = TFile::Open(Form("%s/RefLHC11apass3-2.76TeV.root",path.Data()));
  else if(run > 146860 && run < 156477 )
    ftime = TFile::Open(Form("%s/RefLHC11cpass1-7TeV.root",path.Data()));
  else if(run >= 156477)
    ftime = TFile::Open(Form("%s/RefLHC11cpass1-7TeV.root",path.Data()));
  else if(run <  140000 && run > 136850)
    ftime = TFile::Open(Form("%s/RefLHC10hpass2PbPb2.76TeV.root",path.Data()));
  else if(run < 136850)
    ftime = TFile::Open(Form("%s/RefLHC10dpass2-7TeV.root",path.Data()));
  else printf("Run %d, not considered for time calibration\n",run);
  
  if(ftime){
    printf("AliEMCALRecoUtils - Time recalibration ON\n");
    
    for(Int_t i =0; i< 4; i++)  reco->SetEMCALChannelTimeRecalibrationFactors( i, (TH1F*) ftime->Get(Form("hAllTimeAvBC%d",i)));	
  }
  else printf("AliEMCALRecoUtils --- Time recalibration OFF\n");
  
}


