///
/// \file UpdateEMCAL_OADB_Recalib_Run2_2018Final_PbPb.C
/// \ingroup EMCAL_OADB
/// \brief Update OADB file with energy recalibration for last Run2 Pb-Pb.
///
/// Update Run2 calibration for last Pb-Pb periods, take previously set last 
/// run range [244340,999999] and reset the run range to [244340, 295274]. 
/// Create a new run range [295275,999999] since the SM18-19 online calibration 
/// was updated with offline parameters to allow a uniform triggering in those SM, 
/// so offline parameters are set to 1 for all the channels in SM18-19, and the 
/// other SMs parameters are copied and the same as in the previous run range.
///
/// \author Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-CNRS 
///

#if !defined(__CINT__)
#include <TH2F.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TSystem.h>

#include <Riostream.h>

#include "AliOADBContainer.h"
#endif

 const int kNbSMEMCAL=10;
 const int kNbSMEMCALthird=2;
 const int kNbSMDCAL=6;
 const int kNbSMDCALthird=2;
 const int kNbSMtot=kNbSMEMCAL+kNbSMEMCALthird+kNbSMDCAL+kNbSMDCALthird;

///
/// Main method
///
/// \param fileNameOADB: original OADB file name and path
///
void UpdateEMCAL_OADB_Recalib_Run2_2018Final_PbPb
(
 const char *fileNameOADB="EMCALRecalib.root" 
)
{
  gSystem->Load("libOADB");  
  
  //
  // Open old parameters OADB file
  // Updating the original OADB file, output will be written into BetaRecalib.root 
  //
  AliOADBContainer *con	= new AliOADBContainer("");
  con->InitFromFile(fileNameOADB, "AliEMCALRecalib"); 
  
  //
  // Get the last run range, input for the new one
  //
  
  TObjArray *array15n18p = (TObjArray*)con->GetObject(244340,"LHC15");
  
  //
  // Set new run range for periods before LHC18q (Pb-Pb)
  //
  Int_t lastRun = 295274; // new last run
  con->UpdateObject(con->GetIndexForRun(244340), array15n18p, 244340,lastRun);
  
  //
  // Recover calibrations, needed later for first period with new bias in SM18-19
  //
  TObjArray *pass1  = (TObjArray*)array15n18p->FindObject("pass1");
  TObjArray *recal1 = (TObjArray*)pass1      ->FindObject("Recalib");

  //
  // Create new run range objects from existing ones except last 2 SM
  //
  
  TObjArray *array18q = new TObjArray(kNbSMtot);
  array18q->SetName("LHC18q");
  
  TObjArray *array18q_pass1 = new TObjArray(1);
  TObjArray *array18q_pass2 = new TObjArray(1);
  TObjArray *array18q_pass3 = new TObjArray(1);
  TObjArray *array18q_pass4 = new TObjArray(1);
  
  array18q_pass1->SetName("pass1");
  array18q_pass2->SetName("pass2");
  array18q_pass3->SetName("pass3");
  array18q_pass4->SetName("pass4");
  
  TObjArray *histoArrNew= new TObjArray(kNbSMtot);
  histoArrNew->SetName("Recalib");
  
  //printf("Entries %d \n",recal1->GetEntries());

  char name[30];

  //
  // SM 0/17: copy histograms  to new run range
  //
  for (Int_t iSM=0;iSM<18;iSM++)
  {
    sprintf(name,"EMCALRecalFactors_SM%d",iSM);
//    printf("SM %d, histo name %s, entries %2.4e\n",iSM, name,
//           ((TH2F*) recal1->FindObject(name))->GetEntries());
    histoArrNew->Add(recal1->FindObject(name));
  } // SM 0/17
 
  //
  // SM 18/19: unity histograms to new run range
  //
  for (Int_t iSM=18;iSM<=19;iSM++)
  {
    sprintf(name,"EMCALRecalFactors_SM%d",iSM);
    TH2F* h = (TH2F*) ((TH2F*)recal1->FindObject(name))->Clone();
    h->Reset();
    for(Int_t icol = 0; icol < 48; icol++)
    {
      for(Int_t irow = 0; irow < 8; irow++)
      {    
        h->SetBinContent(icol,irow,1);
      }
    }

    histoArrNew->Add(h);
  } // SM 18/19
  
  array18q_pass1->Add(histoArrNew);
  array18q_pass2->Add(histoArrNew);
  array18q_pass3->Add(histoArrNew);
  array18q_pass4->Add(histoArrNew);
  
  array18q->Add(*&array18q_pass1);
  array18q->Add(*&array18q_pass2);
  array18q->Add(*&array18q_pass3);
  array18q->Add(*&array18q_pass4);
  
  con->AppendObject((TObject*)array18q,lastRun+1,999999);
    
  con->WriteToFile("BetaRecalib.root");

}

