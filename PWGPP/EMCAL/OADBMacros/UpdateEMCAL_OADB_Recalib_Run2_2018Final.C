///
/// \file UpdateEMCAL_OADB_Recalib_Run2_2018Final.C
/// \ingroup EMCAL_OADB
/// \brief Update OADB file with EMCal/DCal energy recalibration for Run2 with 2018 campaign.
///
/// The histograms with energy recalibraton Factors are loaded and some TObjarrays 
/// are filled with these histograms. At the end, a OADB container is created 
/// receiving these arrays.
/// This UpdateEMCAL_OADB_Recalib updates the information of a original OADB file 
/// and writes the output to BetaRecalib.root///
///
/// Update Run2 calibration with the 2018 calibration campaign output. 
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
/// Test what was updated, let's read back the file
///
/// \param runnumber: reference run number
///
void test(int runnumber=195345)
{
  AliOADBContainer *cont=new AliOADBContainer("");
  cont->InitFromFile("BetaRecalib.root", "AliEMCALRecalib");
  // 
  cout<<"_________--------------- dump ---------------------___________"<<endl;
  cont->Dump();
  cout<<"_________--------------- list ---------------------___________"<<endl;
  
  cout<<"cont->GetDefaultList()->Print()"<<endl;
  cont->GetDefaultList()->Print();
  
  TObjArray *recal = (TObjArray*) cont->GetObject(runnumber); //GetObject(int runnumber)
  recal->ls();
  
  TObjArray *recalpass = (TObjArray*) recal->FindObject("pass1");
  
  if(!recalpass)
  {
    cout<<" norecalpass"<<endl;
    return;
  }
  
  TObjArray *recalib = (TObjArray*) recalpass->FindObject("Recalib");
  
  if(!recalib)
  {
    cout<<" no recalib found"<<endl;
    return;
  }
  
  TH2F *h2=(TH2F*)recalib->FindObject("EMCALRecalFactors_SM0");
  if(!h2)
  {
    return;
    cout<<" no histo found"<<endl;
  }
  
  h2->DrawCopy("colz");
  cout<<"That's all folks!"<<endl;
}

///
/// Update OADB Container for EMCal energy recalibration factors
/// from external file. 
///
/// \param fileNameOADB: OADB file name and path
/// \param fileNameRecalibFactors: name and path of input file with new factors
///
void UpdateEMCAL_OADB_Recalib_Run2_2018Final
(
 const char *fileNameOADB="EMCALRecalib.root", 
 // old OADB
 const char *fileNameRecalibFactors="multiplyPi0CalibrationFactors_TextToHisto_Final5.root" 
// Output file of calibration effort with paramters
)
{
  gSystem->Load("libOADB");  

  //  
  // New parameters file from pi0 calibration task on LHC18d calibration campaign
  //
  TFile* fRecalibFactors = new TFile(fileNameRecalibFactors,"read");
  
  // 
  // Open old parameters OADB file
  //
  AliOADBContainer *con	= new AliOADBContainer("");
  con->InitFromFile(fileNameOADB, "AliEMCALRecalib"); //Updating the original OADB file, output will be written into BetaRecalib.root 
  
  //----------------------------------------------------------------------------
  // Update Run2 calibration for EMCal 2012-13 and DCal for periods >= LHC15n
  // duplicate pass1 to pass2, pass3 and pass4
  //----------------------------------------------------------------------------
  
  TObjArray *array15n18 = (TObjArray*)con->GetObject(244340,"LHC15");
  
  con->UpdateObject(con->GetIndexForRun(244340), array15n18, 244340,999999);
  
  TObjArray *pass1  = (TObjArray*)array15n18->FindObject("pass1");
  TObjArray *recal1 = (TObjArray*)pass1     ->FindObject("Recalib"); 
  
  TObjArray *pass2  = (TObjArray*)array15n18->FindObject("pass2");
  TObjArray *recal2 = (TObjArray*)pass2     ->FindObject("Recalib"); // needed later for first Run2 period

  TObjArray *pass3  = (TObjArray*)array15n18->FindObject("pass3");
  TObjArray *recal3 = (TObjArray*)pass3     ->FindObject("Recalib"); // needed later for first Run2 period
   
  TObjArray *pass4  = (TObjArray*)array15n18->FindObject("pass4");
  TObjArray *recal4 = (TObjArray*)pass4     ->FindObject("Recalib"); // needed later for first Run2 period
  
  // Get the histograms per SM in the OADB and the new, update the new multiplying by OADB factors

  char name[30];  
  for (Int_t iSM=0; iSM<kNbSMtot; iSM++)
  {
    sprintf(name,"EMCALRecalFactors_SM%d",iSM);
    printf("SM %d, histo name %s, entries %4.0f\n",
           iSM, name,
           ((TH2F*) recal1->FindObject(name))->GetEntries());

//    for(Int_t icol = 0; icol < 48; icol++)
//    {
//      for(Int_t irow = 0; irow < 24; irow++)
//      {    
//        Float_t org = ( (TH2F*) recal1->FindObject  (name))->GetBinContent(icol,irow);
//        Float_t upd = ( (TH2F*) fRecalibFactors->Get(name))->GetBinContent(icol,irow);
//        ( (TH2F*) recal1->FindObject(name))->SetBinContent(icol,irow,upd*org);
//        ( (TH2F*) recal2->FindObject(name))->SetBinContent(icol,irow,upd*org);
//        ( (TH2F*) recal3->FindObject(name))->SetBinContent(icol,irow,upd*org);
//        ( (TH2F*) recal4->FindObject(name))->SetBinContent(icol,irow,upd*org);
//        
////        printf("\t (%d, %d) : new %1.3f, org %1.3f, new*org %1.3f (%1.3f)\n",
////               icol,irow,upd,org,upd*org,( (TH2F*) recal1->FindObject(name))->GetBinContent(icol,irow));
//      }
//    }
 
    // First, multiply the current calibration factor to the new calibration factor
    // and put it in the histograms containing the new factors
    for(Int_t icol = 0; icol < 48; icol++)
    {
      for(Int_t irow = 0; irow < 24; irow++)
      {    
        Float_t org = ( (TH2F*) recal1->FindObject  (name))->GetBinContent(icol,irow);
        Float_t upd = ( (TH2F*) fRecalibFactors->Get(name))->GetBinContent(icol,irow);
        ( (TH2F*) fRecalibFactors->Get(name))->SetBinContent(icol,irow,upd*org);
        //        printf("\t (%d, %d) : new %1.3f, org %1.3f, new*org %1.3f (%1.3f)\n",
        //               icol,irow,upd,org,upd*org,( (TH2F*) recal1->FindObject(name))->GetBinContent(icol,irow));
      }
    }
 
    // Reset the histograms stored for the different passes
    // since changing the content changes the number of entries, which is misleading
    ((TH2F*) recal1->FindObject(name))->Reset();
    ((TH2F*) recal2->FindObject(name))->Reset();
    ((TH2F*) recal3->FindObject(name))->Reset();
    ((TH2F*) recal4->FindObject(name))->Reset();
    
    // Put the multiplied correction factors in the SM histograms for the different passes
    for(Int_t icol = 0; icol < 48; icol++)
    {
      for(Int_t irow = 0; irow < 24; irow++)
      {    
        Float_t upd = ( (TH2F*) fRecalibFactors->Get(name))->GetBinContent(icol,irow);
        ( (TH2F*) recal1->FindObject(name) )->SetBinContent(icol,irow,upd);
        ( (TH2F*) recal2->FindObject(name) )->SetBinContent(icol,irow,upd);
        ( (TH2F*) recal3->FindObject(name) )->SetBinContent(icol,irow,upd);
        ( (TH2F*) recal4->FindObject(name) )->SetBinContent(icol,irow,upd);
        
        //        printf("\t (%d, %d) : new %1.3f, org %1.3f, new*org %1.3f (%1.3f)\n",
        //               icol,irow,upd,org,upd*org,( (TH2F*) recal1->FindObject(name))->GetBinContent(icol,irow));
      }
    }
    
    printf("SM %d, histo name %s, entries %4.0f, after\n",
           iSM, name,
           ((TH2F*) recal1->FindObject(name))->GetEntries());
  } //iSM
  
  //----------------------------------------------------------------------------
  // Update only EMCal calibration, for periods LHC15a to LHC15m
  // Create pass4 just in case ...
  //----------------------------------------------------------------------------
  
  TObjArray *array15fm = (TObjArray*)con->GetObject(209122,"LHC15fm");
  
  con->UpdateObject(con->GetIndexForRun(209122), array15fm, 209122,244284);
  
  TObjArray *pass1_15fm  = (TObjArray*)array15fm ->FindObject("pass1");
  TObjArray *recal1_15fm = (TObjArray*)pass1_15fm->FindObject("Recalib"); 
  
  TObjArray *pass2_15fm  = (TObjArray*)array15fm ->FindObject("pass2");
  TObjArray *recal2_15fm = (TObjArray*)pass2_15fm->FindObject("Recalib"); // needed later for first Run2 period
  
  TObjArray *pass3_15fm  = (TObjArray*)array15fm ->FindObject("pass3");
  TObjArray *recal3_15fm = (TObjArray*)pass3_15fm->FindObject("Recalib"); // needed later for first Run2 period

  // Update only the Run1 SM from 0 to 10, leave the rest untouched
  for (Int_t iSM=0;iSM<10;iSM++)
  {
    sprintf(name,"EMCALRecalFactors_SM%d",iSM);
    printf("SM %d, histo name %s, entries %4.0f\n",
           iSM, name,
           ((TH2F*) recal1_15fm->FindObject(name))->GetEntries());
    
//    for(Int_t icol = 0; icol < 48; icol++)
//    {
//      for(Int_t irow = 0; irow < 24; irow++)
//      {    
//        Float_t org = ( (TH2F*) recal1_15fm->FindObject(name))->GetBinContent(icol,irow);
//        Float_t upd = ( (TH2F*) fRecalibFactors->Get   (name))->GetBinContent(icol,irow);
//        ( (TH2F*) recal1_15fm->FindObject(name))->SetBinContent(icol,irow,upd*org);
//        ( (TH2F*) recal2_15fm->FindObject(name))->SetBinContent(icol,irow,upd*org);
//        ( (TH2F*) recal3_15fm->FindObject(name))->SetBinContent(icol,irow,upd*org);
//        
////        printf("\t (%d, %d) : new %1.3f, org %1.3f, new*org %1.3f (%1.3f)\n",
////               icol,irow,upd,org,upd*org,( (TH2F*) recal1_15fm->FindObject(name))->GetBinContent(icol,irow));
//      }
//    }
    
    // Reset the histograms stored for the different passes
    // since changing the content changes the number of entries, which is misleading
    ((TH2F*) recal1_15fm->FindObject(name))->Reset();
    ((TH2F*) recal2_15fm->FindObject(name))->Reset();
    ((TH2F*) recal3_15fm->FindObject(name))->Reset();
    
    // Put the multiplied correction factors in the SM histograms for the different passes
    for(Int_t icol = 0; icol < 48; icol++)
    {
      for(Int_t irow = 0; irow < 24; irow++)
      {    
        Float_t upd = ( (TH2F*) fRecalibFactors->Get(name))->GetBinContent(icol,irow); // Modified already for periods >=LHC15n
        ( (TH2F*) recal1_15fm->FindObject(name) )->SetBinContent(icol,irow,upd);
        ( (TH2F*) recal2_15fm->FindObject(name) )->SetBinContent(icol,irow,upd);
        ( (TH2F*) recal3_15fm->FindObject(name) )->SetBinContent(icol,irow,upd);
        
        //        printf("\t (%d, %d) : new %1.3f, org %1.3f, new*org %1.3f (%1.3f)\n",
        //               icol,irow,upd,org,upd*org,( (TH2F*) recal1->FindObject(name))->GetBinContent(icol,irow));
      }
    }

    printf("SM %d, histo name %s, entries %4.0f, after\n",
           iSM, name,
           ((TH2F*) recal1_15fm->FindObject(name))->GetEntries());
  } //iSM
  
  // Create pass4
  TObjArray *array15fm_pass4 = new TObjArray(kNbSMtot);
  
  array15fm_pass4->SetName("pass4");
  
  array15fm_pass4->Add(*&recal1_15fm);
  
  array15fm->Add(*&array15fm_pass4);

  //----------------------------------------------------------------------------
  
  con->WriteToFile("BetaRecalib.root");
  
  //test(195935); // If someone wants to test container
}






