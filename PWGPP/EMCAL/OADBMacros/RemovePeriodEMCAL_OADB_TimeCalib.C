///
/// \file RemovePeriodEMCAL_OADB_TimeCalib.C
/// \ingroup EMCAL_OADB
/// \brief Update OADB file by removing period
///
/// This RemovePeriodEMCAL_OADB_TimeCalib.C removes a given period from an original OADB file EMCALTimeCalib.root
/// and writes the output to EMCALTimeCalib_update1.root
/// Remember to change EMCALTimeCalib_updateX.root to EMCALTimeCalib.root when update done
///
/// ******* Example of Remove periodfrom OADB Container 
/// ******* for EMCal Time calibration
///
///// remove LHC16qrst
///  .x RemovePeriodEMCAL_OADB_TimeCalib.C("OADBfiles_20180120_orginal/EMCALTimeCalib.root",265015,"update21")
///  .x RemovePeriodEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update21.root",265630,"update22")
///  .x RemovePeriodEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update22.root",266405,"update23")
///  .x RemovePeriodEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update23.root",267161,"update24")
///
///  and then update with a new set of calibration constants
///  .x UpdateEMCAL_OADB_TimeCalib.C("EMCALTimeCalib_update24.root","CorrectionFilesTime/Reference_LHC16qrst_merged.root",265015,267166,"TimeCalib16qrst","pass1","update25")
///
/// \author Adam Matyja, <adam.tomasz.matyja@cern.ch>, IFJ PAN Cracow
///
/// update the root file
///
/// \param fileNameOADB: TString, input OADB file
/// \param runMin: Int_t, a run within the period you want to remove from OADB - importent
/// \param outputPrefix: TString, output prefix
///
void RemovePeriodEMCAL_OADB_TimeCalib(TString fileNameOADB="$ALICE_PHYSICS/OADB/EMCAL/EMCALTimeCalib.root",
				      Int_t runMin, 
				      TString outputSufix="update1"){
  AliOADBContainer *con	= new AliOADBContainer("");
  con->InitFromFile(fileNameOADB.Data(), "AliEMCALTimeCalib");

  std::cout<<"Number of periods in OADB: "<<con->GetNumberOfEntries()<<std::endl;
  Int_t index = con->GetIndexForRun(runMin);
  std::cout<<"delete index "<<index<<" period "<<((TObjArray *)con->GetObjectByIndex(index)) ->GetName() <<std::endl;
  con->RemoveObject(index);
  std::cout<<"New number of periods in OADB: "<<con->GetNumberOfEntries()<<std::endl;

  //recreate OADB content (to avoid remnants stored in the file like destructor after removing - BUG??)
  AliOADBContainer* conNew = new AliOADBContainer("AliEMCALTimeCalib");
  for(Int_t index=0;index<con->GetNumberOfEntries();index++){
    //add arrayperiod to main oadb
    conNew->AddDefaultObject((TObject*)con->GetObjectByIndex(index));
    //Establishing run number with the correct objects
    conNew->AppendObject((TObject*)con->GetObjectByIndex(index) ,con->LowerLimit(index),con->UpperLimit(index));
  }
  std::cout<<"Number of entries in recreated OADB: "<<conNew->GetNumberOfEntries()<<std::endl;
  conNew->WriteToFile(Form("EMCALTimeCalib_%s.root",outputSufix.Data()));  

  printf(" written to file\n");
  con->Delete();
}
