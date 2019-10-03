

 char SMP2Name[][100]={
   "SMA0","SMC0","SMA1","SMC1","SMA2","SMC2","SMA3","SMC3",
   "SMA4","SMC4","SMA5","SMC5","SMA9","SMC9","SMA10","SMC10",
   "SMA11","SMC11","SMA12","SMC12"};
 char SMnumber[][100]={
   "0","1","2","3","4","5","6","7","8",
   "9","10","11","12","13","14","15",
   "16","17","18","19"};

 enum detType {kEMCAL,kEMCALthird,kDCAL,kDCALthird};
 int detTypeType[]={kEMCAL,kEMCALthird,kDCAL,kDCALthird};
 char detTypeString[][100]={"EMCAL","EMCALthird","DCAL","DCALthird"};
 int SMdetType[]={
   kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,kEMCAL,
   kEMCAL,kEMCALthird,kEMCALthird,
   kDCAL,kDCAL,kDCAL,kDCAL,kDCAL,kDCAL,kDCALthird,kDCALthird};
 const int kNbColEMCAL=48;
 const int kNbRowEMCAL=24;
 const int kNbSMEMCAL=10;
 const int kNbColEMCALthird=kNbColEMCAL;
 const int kNbRowEMCALthird=(int)(kNbRowEMCAL/3);
 const int kNbSMEMCALthird=2;
 const int kNbColDCAL=32;
 const int kNbRowDCAL=kNbRowEMCAL;
 const int kNbSMDCAL=6;
 const int kNbColDCALthird=kNbColEMCALthird;
 const int kNbRowDCALthird=kNbRowEMCALthird;
 const int kNbSMDCALthird=2;
 const int kNbSMtot=kNbSMEMCAL+kNbSMEMCALthird+kNbSMDCAL+kNbSMDCALthird;
 const int kTabNbCol[4]={kNbColEMCAL,kNbColEMCALthird,kNbColDCAL,kNbColDCALthird};
 const int kTabNbRow[4]={kNbRowEMCAL,kNbRowEMCALthird,kNbRowDCAL,kNbRowDCALthird};
 const int kTabNbSM[4]={kNbSMEMCAL,kNbSMEMCALthird,kNbSMDCAL,kNbSMDCALthird};
 const int kNbColMax=kNbColEMCAL;
 const int kNbRowMax=kNbRowEMCAL;
 const int kNbColOffsetDCAL=kNbColEMCAL-kNbColDCAL;





///
/// \file SetCDB.C
/// \ingroup EMCALOfflineMacrosCalibPi0
/// \brief Put parameters in OCDB file
///
/// Put parameters in OCDB file.
///
/// How to run :
///
///(CAUTION : customize what needs to be)
///
/// aliroot -b
/// .x macros/SetCDB.C
///
///
/// \author Julien Faivre, <Julien.Faivre@cern.ch>, (LPSC-CNRS)
///




///
/// Main method
///
/// Take OCDB used to calibrate existing data, that needs to be updated in
///  alien:/alice/data/2011/OCDB/EMCAL/Calib/Data/Run144484_999999999_v5_s0.root, 
/// take the latest version
///
void SetCDB(void)
{  
   // Julien : this is the OCDB file with the coeffs that have been used to process the data used for the calibration :
   TFile * forg = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/createOCDB_4_with2015data/Run0_999999999_v1_s0_OCDBcoeffsLHC15ij_EMCALnotValid.root","READ");
  
  AliCDBEntry * cdborg = (AliCDBEntry*) forg->Get("AliCDBEntry");    
  AliEMCALCalibData * cparamorg =  cdborg->GetObject();

  //New OCDB container
  AliEMCALCalibData *cparamnew=new AliEMCALCalibData("EMCAL");

  //Get the file with the recalibration factors
  //CUSTOMIZE customize :
  TFile * frecal = new TFile("/cebaf/cebaf/EMCAL/calibPi0_run2/createOCDB_4_with2015data/multiplyPi0CalibrationFactors_TextToHisto_Final.root","READ");
  //File with all coeffs at 1.0, for tests :
  //TFile * frecal = new TFile("/cebaf/cebaf/EMCAL/calibPi0/RecalibrationFactors_allOneForTests_final.root","READ");
  
  //cout<<"============== "<<cparamorg->GetName()<<" ==============="<<endl;
  
  //cparam->Print("gain");
  //cparam->Print("ped");
  // Go tower by tower, get the old calibration value and multiply by the recalibration factor
  // Add the new calibration factor to new OCDB
  for(Int_t m = 0; m < kNbSMtot; m++)
  {
    TH2F* h = frecal->Get(Form("EMCALRecalFactors_SM%d",m));
    for (Int_t c = 0; c < kTabNbCol[SMdetType[m]]; c++){
      for (Int_t r = 0; r < kTabNbRow[SMdetType[m]]; r++){
        Float_t adcorg=cparamorg->GetADCchannel(m,c,r);
        Float_t recal = h->GetBinContent(c,r);
        Float_t newadc = adcorg*recal;
        //printf("m %d, c %d, r %d: org %1.4f, recal %1.4f, new %1.4f\n",m,c,r,adcorg, recal, newadc);
        //printf("m %d, c %d, r %d: org %1.4f\n",m,c,r,newadc);
        //if(recal ==1) printf("m %d, c %d, r %d: org %1.4f, recal %1.4f, new %1.4f\n",m,c,r,adcorg, recal, newadc);
        cparamnew->SetADCchannel(m,c,r,newadc);
      }
    }
    
  }


  //Create OCDB File
  AliCDBMetaData md;
  //CUSTOMIZE customize :
  md.SetComment("bbbb");
  md.SetBeamPeriod(0);
  md.SetResponsible("aaaa");
  md.SetAliRootVersion(gSystem->Getenv("ARVERSION"));

  // Careful, select here the first run where this calibration is valid
  //CUSTOMIZE customize :
  // First run in 2012 : 176326, last run in 2011 : 170593
  // for 2012 calibration, from LHC12a 172439
  AliCDBId id("EMCAL/Calib/Data",172439,AliCDBRunRange::Infinity()); // create in EMCAL/Calib/Data DBFolder 

  AliCDBManager* man = AliCDBManager::Instance();  
  AliCDBStorage* loc = man->GetStorage("local://$ALICE_ROOT/OCDB/");
  loc->Put(cparamnew, id, &md);

 printf("Done.\n");
 return;
 }

