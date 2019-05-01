///
/// \file PrintScaleFactors.C
/// \ingroup CaloTrackCorrMacros
/// \brief Print pT-hard cross section scale factor
///
/// Macro that prints the cross section scaling factors for MC pT-hard binned productions
/// from the output of CaloTrackCorr analysis histgrams.
/// The output histograms must be placed in different directories with ProductionName/PtHardBinNumber/FileName.root
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TH1F.h>
#include <TList.h>
#include <TString.h>
#include <TFile.h>

#endif

///
/// Get the cross section value for a given pT hard bin
///
/// \param bin      : pT hard bin number
/// \param prodName : Name of directory with production factors to be printed
/// \param fileName : Name of analysis output 
/// \param listName : Name of analysis output inside the output file
///
void ForBin
(
 Int_t   bin      = -1, 
 TString prodName = "pp_7TeV_GJ",
 TString fileName = "Merged",
 TString listName = "CTC_EMCAL_Trig_default_ClV1_Ecell100_Eseed500_SPDPileUp"
 )
{ 
  //printf("Scale bin %d\n",bin);
  
  TFile * f= TFile::Open(Form("%s/%d/%s.root",prodName.Data(),bin,fileName.Data()),"read");

  if ( !f )
  {
    //printf("No file available for bin %d: \n",bin);
    return;
  }
  
  TList* list = (TList*) f->Get(listName);

  if ( !list )
  {
    //printf("No list available\n");
    return;
  }
  
  // Extract the cross section from the corresponding histograms in the file

  Float_t nEventsIn = ((TH1F*)list->FindObject("hNEventsIn"))->GetEntries();
  Float_t nEvents   = ((TH1F*)list->FindObject("hNEvents"))  ->GetEntries();
  Float_t entries   = ((TH1F*)list->FindObject("hXsec"))     ->GetEntries();
  Float_t xsec      = ((TH1F*)list->FindObject("hXsec"))     ->GetBinContent(1)   ;
  Float_t trials    = ((TH1F*)list->FindObject("hTrials"))   ->GetBinContent(1) ;
  
  
  Float_t xsecNorm    = xsec  / entries; // average per chunk
  Float_t trialsNorm  = trials/ entries; // trials per chunk
  
  Double_t scale  = xsecNorm/trialsNorm/nEventsIn ;

  printf("Bin %d, nEvents Input %2.4e, nEvents selected %2.4e, Events accepted %2.2f, nXsec files %2.4e,"
         "xsec %2.4e, xSecNorm=xsec/nFiles %2.4e trials %2.4e, trialsNorm=trials/nFiles %2.4e, "
         "factor=xSecNorm/trialsNorm  %2.4e, factor/nEventsInput = %2.4e \n", 
         bin, nEventsIn, nEvents, nEvents/nEventsIn, entries,
         xsec, xsecNorm, trials, trialsNorm, 
         xsecNorm/trialsNorm,scale);
  
  f->Close();
  delete f;
}

///
/// Main method, get the cross section value for all or a given pT hard bin
///
/// \param bin      : pT hard bin number, if "-1", bins from 0 to 20 are checked and printed
/// \param prodName : Name of directory with production factors to be printed
/// \param fileName : Name of analysis output 
/// \param listName : Name of analysis output inside the output file
///
void PrintScaleFactors
(
 TString prodName = "pp_7TeV_GJ",
 TString fileName = "Merged",
 TString listName = "CTC_EMCAL_Trig_default_ClV1_Ecell100_Eseed500_SPDPileUp",
 Int_t   bin      = -1
 )
{ 
  printf("production %s, file %s, list %s\n",prodName.Data(),fileName.Data(),listName.Data());
  
  if ( bin < 0 ) 
    for(Int_t i = 0; i < 20; i++) ForBin(i,prodName,fileName,listName);
  else 
    ForBin(bin,prodName,fileName,listName);
}
