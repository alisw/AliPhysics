///
/// \file UpdateEMCAL_OADB_Recalib11hMC.C
/// \ingroup EMCAL_OADB
/// \brief Update OADB file with energy recalibration factors for MC LHC14a1a
///
/// MC production LHC14a1a was screwed up and decalibrated. Here it is the macro 
/// used to put the corresponding factors in OADB.
///
/// \author Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-CNRS ???
/// \author Marcel Figueredo, <marcel.figueredo@cern.ch>, Sao Paulo
///

#if !defined(__CINT__)
#include "CreateEMCAL_OADB_CalibReference.C"
#endif

///
/// Main method
///
void UpdateEMCAL_OADB_Recalib11hMC(const char *fileNameOADB="$ALICE_ROOT/OADB/EMCAL/EMCALRecalib.root")
{
  gSystem->Load("libOADB");  
  
  AliOADBContainer *con	= new AliOADBContainer("");
  con->InitFromFile(fileNameOADB, "AliEMCALRecalib"); //Updating the original OADB file, output will be written into BetaRecalib.root 
  
  TObjArray *array11hMC = GetArrayEnergy(148000,"Run144484_999999999_v5_s0.root","Recalib",0.0162); //last is factor for division.
  array11hMC->SetName("Recalib");
  
  TObjArray *array11hPass = new TObjArray(0);
  array11hPass->SetName("LHC14a1a");
  array11hPass->Add(*&array11hMC);
  
  TObjArray *array11h = new TObjArray(0);
  array11h->SetName("LHC11h");
  array11h->Add(*&array11hPass);
  
  con->AddDefaultObject(*&array11h);
  con->AppendObject(*&array11h,167693,170593);
  con->WriteToFile("BetaRecalib.root");
}
