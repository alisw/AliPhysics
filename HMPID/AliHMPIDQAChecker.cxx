/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


/* $Id$ */

//...
//  Checks the quality assurance. 
//  By comparing with reference data
//  Skeleton for HMPID
//...

// --- ROOT system ---
#include <TClass.h>
#include <TH1F.h> 
#include <TH1I.h>
#include <TH2.h> 
#include <TF1.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h>
#include <TLine.h>
#include <TParameter.h> 
#include <TPaveText.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliHMPIDQAChecker.h"
#include "AliCDBEntry.h"
#include "AliQAManager.h"
#include "AliQAThresholds.h"

ClassImp(AliHMPIDQAChecker)
 //_________________________________________________________________
AliHMPIDQAChecker::AliHMPIDQAChecker() : 
AliQACheckerBase("HMPID","HMPID Quality Assurance Data Checker"), 
fNoReference(kTRUE),
fQARefRec(NULL),

fHmpQaThr_NumberOfExcludedDDL(0),
fHmpQaThr_DataSizeLowerThreshold(900),
fHmpQaThr_DataSizeUpperThreshold(1500),
fHmpQaThr_PadOccupancyLowerThreshold(0.005),
fHmpQaThr_PadOccupancyUpperThreshold(0.8),
fHmpQaThr_SectorGainLossWarningThreshold(3),
fHmpQaThr_SectorGainLossErrorThreshold(6),
fHmpQaThr_MissingPadFractionWarningThreshold(0.3),
fHmpQaThr_MissingPadFractionErrorThreshold(0.5),
fIsOnlineThr(kFALSE)                                 
 


{
    //ctor, fetches the reference data from OCDB 
  char * detOCDBDir = Form("HMPID/%s/%s", AliQAv1::GetRefOCDBDirName(), AliQAv1::GetRefDataDirName()) ; 
  AliCDBEntry * QARefRec = AliQAManager::QAManager()->Get(detOCDBDir);
  if(QARefRec) {
    fQARefRec = dynamic_cast<TObjArray*> (QARefRec->GetObject()) ; 
    if (fQARefRec)
      if (fQARefRec->GetEntries()) 
        fNoReference = kFALSE ;            
    if (fNoReference) 
      AliInfo("QA reference data NOT retrieved for Reconstruction check. No HMPIDChecker!");
  }
    
}

//_________________________________________________________________
AliHMPIDQAChecker::~AliHMPIDQAChecker() 
{
  if(fQARefRec) { fQARefRec->Delete() ;   delete fQARefRec ; }
}
//_________________________________________________________________
void AliHMPIDQAChecker::Check(Double_t *  check, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * /*recoParam*/) 
{
//
// Main check function: Depending on the TASK, different checks are applied
// At the moment:       check for empty histograms and checks for RecPoints

  InitOnlineThresholds();    
  
  if(fNoReference)  

  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    check[specie] = 1.0;    
    //printf("+++++++++++++++++++++ specie %d name: %s \n",specie,AliRecoParam::GetEventSpecieName(specie));
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) continue ;     
    // checking for empy histograms
    if(CheckEntries(list[specie]) == 0)  {
      AliWarning("histograms are empty");
      check[specie] = 0.4;//-> Corresponds to kWARNING see AliQACheckerBase::Run
    }
  
    check[specie] = AliQAv1::kINFO ;
    
    
    //check sim
    if(index == AliQAv1::kSIM) check[specie] = CheckSim(list[specie], fQARefRec);

    // checking rec points
    if(index == AliQAv1::kREC) check[specie] = CheckRec(list[specie], fQARefRec);
   
    //checking raw data
    if(index == AliQAv1::kRAW) {       check[specie] = CheckRaw(specie,list[specie]);       }
                               
    
  } // species loop
}
//_________________________________________________________________
Double_t AliHMPIDQAChecker::CheckEntries(TObjArray * list) const
{
  //
  //  check on the QA histograms on the input list: 
  // 

  Double_t test = 0.0  ;
  Int_t count = 0 ; 
  
  if (list->GetEntries() == 0){  
    test = 1. ; // nothing to check
  }
  else {
    TIter next(list) ; 
    TH1 * hdata ;
    count = 0 ; 
    while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
      if (hdata) { 
	Double_t rv = 0.;
	if(hdata->GetEntries()>0)rv=1; 
	count++ ; 
	test += rv ; 
      }
      else{
	AliError("Data type cannot be processed") ;
      }
      
    }
    if (count != 0) { 
      if (test==0) {
	AliWarning("Histograms are booked for THIS specific Task, but they are all empty: setting flag to kWARNING");
	test = 0.;  //upper limit value to set kWARNING flag for a task
      }
      else test = 1 ;
    }
  }

  return test ; 
}  
//_________________________________________________________________
Double_t AliHMPIDQAChecker::CheckSim(TObjArray *listsim, TObjArray *listref) const
{
  //
  //  check on the HMPID RecPoints by using expo fit and Kolmogorov Test:
  //

   Float_t checkresponse = 0;

   Float_t counter = 0 ;
   TIter next(listsim) ;
   TH1* histo;
   while ( (histo = dynamic_cast<TH1 *>(next())) ) {
     //PH The histogram should have at least 10 bins with at least 5 entries
     Int_t nbinsabove = 0;
     for (Int_t ibin=histo->FindBin(1); ibin<=histo->FindBin(50); ibin++) { 
       if (histo->GetBinContent(ibin)>5) nbinsabove++;
     }

   if( nbinsabove < 10 ) counter++;
   else {
    TString h = histo->GetTitle();
    if(h.Contains("Zoom")){
    histo->Fit("expo","LQ0","",5,50);
    if(histo->GetFunction("expo")->GetParameter(1) !=0 ) if(TMath::Abs((-1./(histo->GetFunction("expo"))->GetParameter(1)) - 35 ) > 5) counter++;
    }
    if(h.Contains("size  MIP"))   if(TMath::Abs(histo->GetMean()-5) > 2) counter++;
    if(h.Contains("size  Phots")) if(TMath::Abs(histo->GetMean()-2) > 2) counter++;
    if(h.Contains("distribution")) if(histo->KolmogorovTest((TH1F *)listref->At(0))<0.8) counter++;
    AliDebug(AliQAv1::GetQADebugLevel(),Form(" Kolm. test : %f ",histo->KolmogorovTest((TH1F *)listref->At(0))));  
   }
  }
 Float_t response = counter/(7.+7.+42.+42.); // 7.+7.+42 +42 = N checked histograms (-> To be replaced by listsim->GetEntries())
 
 if(response < 0.1) checkresponse = 0.7;      // <10% of the check histograms show a failing check -> Corresponds to kINFO see AliQACheckerBase::Run
 else if(response < 0.5) checkresponse = 0.4; //  50%  of the check histograms show a failing check -> Corresponds to kWARNING see AliQACheckerBase::Run
 else checkresponse = 0.001;                  // > 50% of the check histograms show a failing check -> Corresponds to kERROR see AliQACheckerBase::Run
 return checkresponse;
}

//___________________________________________________________________________________________________
Double_t AliHMPIDQAChecker::CheckRec(TObjArray *listrec, TObjArray *listref) const
{
  //
  //  check on the HMPID RecPoints by using expo fit and Kolmogorov Test:
  //

   Float_t checkresponse = 0;

   Float_t counter = 0 ;
   TIter next(listrec) ;
   TH1* histo;
   while ( (histo = dynamic_cast<TH1 *>(next())) ) {
     //PH The histogram should have at least 10 bins with at least 5 entries
     Int_t nbinsabove = 0;
     for (Int_t ibin=histo->FindBin(1); ibin<=histo->FindBin(50); ibin++) { 
       if (histo->GetBinContent(ibin)>5) nbinsabove++;
     }

   if( nbinsabove < 10 ) counter++;
   else {
    TString h = histo->GetTitle();
    if(h.Contains("Zoom")){
    histo->Fit("expo","LQ0","",5,50);
    if(histo->GetFunction("expo")->GetParameter(1) !=0 ) if(TMath::Abs((-1./(histo->GetFunction("expo"))->GetParameter(1)) - 35 ) > 5) counter++;
    }
    if(h.Contains("size  MIP"))   if(TMath::Abs(histo->GetMean()-5) > 2) counter++;
    if(h.Contains("size  Phots")) if(TMath::Abs(histo->GetMean()-2) > 2) counter++;
    if(h.Contains("distribution")) if(histo->KolmogorovTest((TH1F *)listref->At(0))<0.8) counter++;
    AliDebug(AliQAv1::GetQADebugLevel(),Form(" Kolm. test : %f ",histo->KolmogorovTest((TH1F *)listref->At(0))));  
   }
  }
 Float_t response = counter/(7.+7.+42.+42.); // 7.+7.+42 +42 = N checked histograms (-> To be replaced by listrec->GetEntries())
 
 if(response < 0.1) checkresponse = 0.7;      // <10% of the check histograms show a failing check -> Corresponds to kINFO see AliQACheckerBase::Run
 else if(response < 0.5) checkresponse = 0.4; //  50%  of the check histograms show a failing check -> Corresponds to kWARNING see AliQACheckerBase::Run
 else checkresponse = 0.001;                  // > 50% of the check histograms show a failing check -> Corresponds to kERROR see AliQACheckerBase::Run
 return checkresponse;
}
//___________________________________________________________________________________________________
Double_t AliHMPIDQAChecker::CheckRaw(Int_t specie, TObjArray* list)
{
  //
  // Check the raw data for offline / online using default or updated thresholds from AMORE
  // As of now (06/07/2012) the quality flag of all histos in raw will be est to the result of teh CheckRaw and displayed
  // in AMORE. But we can pu undividual labels.
  //
  
  Int_t raqQualFlag = AliQAv1::kNULLBit;
  
  Int_t hmpQaFlags[4]={-1}; //init for the 4 shifter histos
  
  TString histname ="null";
  TPaveText text(0.65,0.8,0.9,0.99,"NDC"); 
  TPaveText text1(0.65,0.8,0.9,0.99,"NDC"); 
   
 
  Int_t entries = list->GetEntriesFast();
  
  if ( entries == 0 ) {
     AliWarning(Form("HMPID QA Checker RAWS: no object to analyse! Exiting..."));
     return AliQAv1::kFATAL;
  }
  
  TLine* lineDdlDataSizeMin = new TLine(1536,fHmpQaThr_DataSizeUpperThreshold,1548,fHmpQaThr_DataSizeUpperThreshold);
  TLine* lineDdlDataSizeMax = new TLine(1536,fHmpQaThr_DataSizeLowerThreshold,1548,fHmpQaThr_DataSizeLowerThreshold);
  
  TLine* linePadOccupMin = new TLine(1536,fHmpQaThr_PadOccupancyUpperThreshold,1548,fHmpQaThr_PadOccupancyUpperThreshold);
  TLine* linePadOccupMax = new TLine(1536,fHmpQaThr_PadOccupancyLowerThreshold,1548,fHmpQaThr_PadOccupancyLowerThreshold);
  
  
  Int_t badDdlCnt = 0, badOccCnt = 0;
  
  //___ check data size per ddl
  if( list->FindObject(Form("%s_hHmpDdlDataSize",AliRecoParam::GetEventSpecieName(specie)) )) 
    {
      TH1* h1 = dynamic_cast<TH1*>( list->FindObject(Form("%s_hHmpDdlDataSize",AliRecoParam::GetEventSpecieName(specie)))); 
      if(h1) {
      if( h1->Integral() > 1 ) {
      // no entres -> fatal
      if( h1->GetEntries() == 0) {raqQualFlag = AliQAv1::kFATAL;}
      h1->SetStats(0);
      
      
      // clean up the text, stat, lines, ...
      if( h1->GetListOfFunctions() ) h1->GetListOfFunctions()->Clear();
      
      for ( Int_t iddl = 1 ; iddl <= 14; iddl++)
	{
	  if( h1->GetBinContent(iddl) < fHmpQaThr_DataSizeLowerThreshold || h1->GetBinContent(iddl) > fHmpQaThr_DataSizeUpperThreshold ) badDdlCnt++; 
	}
      //___ check if one or more DDLs are excluded
        
      
      badDdlCnt -= fHmpQaThr_NumberOfExcludedDDL;
      
      
      
      //___ check how many are bad
      if ( badDdlCnt == 0 )  
	{         
	  hmpQaFlags[0] = AliQAv1::kINFO;
	  text1.Clear();
	  text1.AddText(Form("OK (%d)",fIsOnlineThr));
	  text1.SetFillColor(kGreen);
	  h1->GetListOfFunctions()->Add((TPaveText*)text1.Clone());	 
	  lineDdlDataSizeMin->SetLineColor(kGreen);
	  lineDdlDataSizeMax->SetLineColor(kGreen);
	  h1->GetListOfFunctions()->Add((TLine*)lineDdlDataSizeMin->Clone());
	  h1->GetListOfFunctions()->Add((TLine*)lineDdlDataSizeMax->Clone());
	}
      else if ( badDdlCnt == 1 )  
	{         
	  hmpQaFlags[0]  = AliQAv1::kWARNING;
	  text1.Clear();
	  text1.AddText(Form("WARNING CHECK TWIKI (%d)",fIsOnlineThr));
	  text1.SetFillColor(kOrange);
	  h1->GetListOfFunctions()->Add((TPaveText*)text1.Clone());	
	  lineDdlDataSizeMin->SetLineColor(kOrange);
	  lineDdlDataSizeMax->SetLineColor(kOrange);
	  h1->GetListOfFunctions()->Add((TLine*)lineDdlDataSizeMin->Clone());
	  h1->GetListOfFunctions()->Add((TLine*)lineDdlDataSizeMax->Clone());   
	}
      else if (  badDdlCnt >= 2 )  
	{
	  hmpQaFlags[0]  = AliQAv1::kERROR;         
	  text1.Clear();
	  text1.AddText(Form("ERROR CALL ONCALL (%d)",fIsOnlineThr));
	  text1.SetFillColor(kRed);
	  h1->GetListOfFunctions()->Add((TPaveText*)text1.Clone());	      
	  lineDdlDataSizeMin->SetLineColor(kRed);
	  lineDdlDataSizeMax->SetLineColor(kRed);
	  h1->GetListOfFunctions()->Add((TLine*)lineDdlDataSizeMin->Clone());
	  h1->GetListOfFunctions()->Add((TLine*)lineDdlDataSizeMax->Clone());   
	}
      else 
	{
	  hmpQaFlags[0]  = AliQAv1::kFATAL;
	  text1.Clear();
	  text1.AddText(Form("FATAL CALL ONCALL (%d)",fIsOnlineThr));
	  text1.SetFillColor(kRed);
	  h1->GetListOfFunctions()->Add((TPaveText*)text1.Clone());	      
	  lineDdlDataSizeMin->SetLineColor(kRed);
	  lineDdlDataSizeMax->SetLineColor(kRed);
	  h1->GetListOfFunctions()->Add((TLine*)lineDdlDataSizeMin->Clone());
	  h1->GetListOfFunctions()->Add((TLine*)lineDdlDataSizeMax->Clone());   
	}
      }
      }//the histo is filled
    }//___hHmpDdlDataSize
  
  
  
  if( list->FindObject(Form("%s_fHmpPadOcc",AliRecoParam::GetEventSpecieName(specie)) )) 
    {
     
      
      TH1* h1 = dynamic_cast<TH1*>( list->FindObject(Form("%s_fHmpPadOcc",AliRecoParam::GetEventSpecieName(specie)))); 
       if(h1) {
       if( h1->Integral() > 1 ) {
      // no entres -> fatal
      if( h1->GetEntries() == 0) {raqQualFlag = AliQAv1::kFATAL;}
       h1->SetStats(0);
      // clean up the text, stat, lines, ...
      if( h1->GetListOfFunctions() ) h1->GetListOfFunctions()->Clear();
      
      for ( Int_t iddl = 1 ; iddl <= 14; iddl++)
	{
	  if( h1->GetBinContent(iddl) < fHmpQaThr_PadOccupancyLowerThreshold || h1->GetBinContent(iddl) > fHmpQaThr_PadOccupancyUpperThreshold ) badOccCnt++; 
	}
       
      badOccCnt -= fHmpQaThr_NumberOfExcludedDDL;
         
      //___ check how many are bad
      if ( badOccCnt == 0 )  
	{
	  hmpQaFlags[1]  = AliQAv1::kINFO;
	  text.Clear();
	  text.AddText(Form("OK (%d)",fIsOnlineThr));
	  text.SetFillColor(kGreen);
	  h1->GetListOfFunctions()->Add((TPaveText*)text.Clone());	
	  linePadOccupMin->SetLineColor(kGreen);
	  linePadOccupMax->SetLineColor(kGreen);
	  h1->GetListOfFunctions()->Add((TLine*)linePadOccupMin->Clone());
	  h1->GetListOfFunctions()->Add((TLine*)linePadOccupMax->Clone());         
	}
      else if ( badOccCnt == 1 )  
	{
	  hmpQaFlags[1]  = AliQAv1::kWARNING;
	  text.Clear();
	  text.AddText(Form("WARNING CHECK TWIKI (%d)",fIsOnlineThr));
	  text.SetFillColor(kOrange);
	  h1->GetListOfFunctions()->Add((TPaveText*)text.Clone());	
	  linePadOccupMin->SetLineColor(kGreen);
	  linePadOccupMax->SetLineColor(kGreen);
	  h1->GetListOfFunctions()->Add((TLine*)linePadOccupMin->Clone());
	  h1->GetListOfFunctions()->Add((TLine*)linePadOccupMax->Clone());
	}
      else if (  badOccCnt == 2 )  
	{
	  hmpQaFlags[1]  = AliQAv1::kERROR;
	  text.Clear();
	  text.AddText(Form("ERROR CALL ONCALL (%d)",fIsOnlineThr));
	  text.SetFillColor(kRed);
	  h1->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
	  linePadOccupMin->SetLineColor(kGreen);
	  linePadOccupMax->SetLineColor(kGreen);
	  h1->GetListOfFunctions()->Add((TLine*)linePadOccupMin->Clone());
	  h1->GetListOfFunctions()->Add((TLine*)linePadOccupMax->Clone());
	}
      else 
	{
	  hmpQaFlags[1] = AliQAv1::kFATAL;
	  text.Clear();
	  text.AddText(Form("FATAL CALL ONCALL (%d)",fIsOnlineThr));
	  text.SetFillColor(kRed);
	  h1->GetListOfFunctions()->Add((TPaveText*)text.Clone());	      
	  linePadOccupMin->SetLineColor(kGreen);
	  linePadOccupMax->SetLineColor(kGreen);
	  h1->GetListOfFunctions()->Add((TLine*)linePadOccupMin->Clone());
	  h1->GetListOfFunctions()->Add((TLine*)linePadOccupMax->Clone());   
	}
      }
    }
    }//___HmpPadOcc
  
  Int_t sumPadMapChROR[7]={0};
  Int_t sumPadMapChROL[7]={0};
  Int_t bigMapFlag =  AliQAv1::kINFO;
  Int_t errCntBigMap = 0;

  if( list->FindObject(Form("%s_hHmpBigMap",AliRecoParam::GetEventSpecieName(specie)) )) 
    {
      TH2* h2 = dynamic_cast<TH2*>( list->FindObject(Form("%s_hHmpBigMap",AliRecoParam::GetEventSpecieName(specie)))); 
      if(h2) {
        if( h2->Integral() > 1 ) {
      // no entres -> fatal
      if( h2->GetEntries() == 0) {raqQualFlag = AliQAv1::kFATAL;}
       h2->SetStats(0);
      // clean up the text, stat, lines, ...
      if( h2->GetListOfFunctions() ) h2->GetListOfFunctions()->Clear();

      //calculate missing pad fraction
      for(Int_t ich = 0; ich < 7; ich++)
	{
	  for(Int_t iy=1+ich*144;iy<=144+ich*144;iy++) {
	    for(Int_t ix=1;ix<=80;ix++)    if(h2->GetBinContent(ix,iy) > 0) sumPadMapChROL[ich]++;
	    for(Int_t ix=81;ix<=160;ix++)  if(h2->GetBinContent(ix,iy) > 0) sumPadMapChROR[ich]++;        
	  }
	}//ch loop
       //check the calculated missing pad fraction
      for(Int_t ich = 0; ich < 7; ich++)
	{
          
          bigMapFlag =  AliQAv1::kINFO; 
	  if( (1-sumPadMapChROL[ich]/1.0/11520) > fHmpQaThr_MissingPadFractionWarningThreshold ||
	      (1-sumPadMapChROR[ich]/1.0/11520) > fHmpQaThr_MissingPadFractionWarningThreshold ) 
          {
            bigMapFlag =  AliQAv1::kWARNING;
            }
	  if( (1-sumPadMapChROL[ich]/1.0/11520) > fHmpQaThr_MissingPadFractionErrorThreshold ||
	      (1-sumPadMapChROR[ich]/1.0/11520) > fHmpQaThr_MissingPadFractionErrorThreshold ) {
             bigMapFlag =  AliQAv1::kERROR; 
             errCntBigMap++;
           }	  
          
     	}//ch loop
      if( errCntBigMap > 0 ) bigMapFlag =  AliQAv1::kERROR;

      //update labels
       if (  bigMapFlag == AliQAv1::kINFO )  
	{
	  hmpQaFlags[2]  = AliQAv1::kINFO;
	  text.Clear();
	  text.AddText(Form("OK (%d)",fIsOnlineThr));
	  text.SetFillColor(kGreen);
	  h2->GetListOfFunctions()->Add((TPaveText*)text.Clone());	
	}
       else if (  bigMapFlag == AliQAv1::kWARNING )  
	{
	  hmpQaFlags[2]  = AliQAv1::kWARNING;
	  text.Clear();
	  text.AddText(Form("WARNING CHECK TWIKI (%d)",fIsOnlineThr));
	  text.SetFillColor(kOrange);
	  h2->GetListOfFunctions()->Add((TPaveText*)text.Clone());	
	}
       else if (  bigMapFlag == AliQAv1::kERROR )  
	{
	  hmpQaFlags[2]  = AliQAv1::kERROR;
	  text.Clear();
	  text.AddText(Form("ERROR CALL ONCALL (%d)",fIsOnlineThr));
	  text.SetFillColor(kRed);
	  h2->GetListOfFunctions()->Add((TPaveText*)text.Clone());	
	}
       else
	 {
	   hmpQaFlags[2]  = AliQAv1::kFATAL;
	   text.Clear();
	   text.AddText(Form("FATAL CALL ONCALL (%d)",fIsOnlineThr));
	   text.SetFillColor(kRed);
	   h2->GetListOfFunctions()->Add((TPaveText*)text.Clone());	
	 }      
        }       
     }
    }//___HmpBigMap
  

  Int_t numSectorsMissing = 0, numSectorsGainLoss = 0;
  Double_t  hSumQPerSector[42]={0};

  if( list->FindObject(Form("%s_fHmpHvSectorQ",AliRecoParam::GetEventSpecieName(specie))))
    {
      
      TH2* h2 = dynamic_cast<TH2*>( list->FindObject(Form("%s_fHmpHvSectorQ",AliRecoParam::GetEventSpecieName(specie)))); 
      if(h2) {
        if(h2->Integral() > 0 ) {
      // no entres -> fatal
      if( h2->GetEntries() == 0) {raqQualFlag = AliQAv1::kFATAL;}
             h2->SetStats(0);
      // clean up the text, stat, lines, ...
      if( h2->GetListOfFunctions() ) h2->GetListOfFunctions()->Clear();
     
      //___ check sectors 
      for(Int_t isec = 1 ; isec <= 42; isec++)
	{
	  for(Int_t ibiny=100;ibiny<410;ibiny++) {hSumQPerSector[isec-1] += h2->GetBinContent(ibiny,isec); }
	  if(hSumQPerSector[isec-1]  < 0.001) {numSectorsGainLoss++;}  // there is no photon and mip peak, gain loss
	  if(h2->GetBinContent(1,isec) < 0.01 ) {numSectorsMissing++; } //practically there is no charge , the sector is missing
	}
      Int_t  sectorErrors = numSectorsGainLoss+numSectorsMissing;
     
      
	if ( sectorErrors <=  3)
	  {
	    hmpQaFlags[3]  = AliQAv1::kINFO;
	    text.Clear();
	    text.AddText(Form("OK (%d)",fIsOnlineThr));
	    text.SetFillColor(kGreen);
	    h2->GetListOfFunctions()->Add((TPaveText*)text.Clone());	   
	  }
	else if ( sectorErrors > fHmpQaThr_SectorGainLossWarningThreshold)
	  {
	    hmpQaFlags[3]  = AliQAv1::kWARNING;
	    text.Clear();
	    text.AddText(Form("WARNING CHECK TWIKI (%d)",fIsOnlineThr));
            if(numSectorsMissing > 0 ) text.AddText(Form("MISSING SECTOR?"));              
	    text.SetFillColor(kOrange);
	    h2->GetListOfFunctions()->Add((TPaveText*)text.Clone());
	  }
	else if ( sectorErrors > fHmpQaThr_SectorGainLossErrorThreshold)
	  {
	    hmpQaFlags[3]  = AliQAv1::kERROR;
	    text.Clear();
	    text.AddText(Form("ERROR CALL ONCALL (%d)",fIsOnlineThr));
            if(numSectorsMissing > 0 ) text.AddText(Form("MISSING SECTOR?"));
	    text.SetFillColor(kRed);
	    h2->GetListOfFunctions()->Add((TPaveText*)text.Clone());	
	  }
	else
	  {
	    hmpQaFlags[3]  = AliQAv1::kFATAL;
	    text.Clear();
	    text.AddText(Form("FATAL CALL ONCALL (%d)",fIsOnlineThr));
	    text.SetFillColor(kRed);
	    h2->GetListOfFunctions()->Add((TPaveText*)text.Clone());	
	  }
        }
    }
  }
  
  
  //del lines, ...
  lineDdlDataSizeMin->Delete();
  lineDdlDataSizeMax->Delete();
  linePadOccupMin->Delete();
  linePadOccupMax->Delete();
  
  Double_t dflag = -1;
  switch ( TMath::MaxElement(4,hmpQaFlags))
    {
    case  AliQAv1::kINFO:
      dflag = 1.0;
      
      break;
    case AliQAv1::kWARNING:
      dflag = 0.75;
      
      break;
    case AliQAv1::kERROR:
      dflag = 0.25;
      
      break;
    case AliQAv1::kFATAL:
      dflag = -1.0;
     
      break;
    default:
      dflag = AliQAv1::kNULLBit;
     
      break;
    }	  
     
  return dflag;
  
  
  
}
//___________________________________________________________________________________________________
void AliHMPIDQAChecker::InitOnlineThresholds()
{
  //
  // Init the online thresholds from GRP generated by AMORE
  //
  
  AliCDBManager* man = AliCDBManager::Instance(); 
  if(!man)  {     fIsOnlineThr = kFALSE;     return;   }
  

  if(!man->Get("GRP/Calib/QAThresholds"))  {     fIsOnlineThr = kFALSE;    return;    }
 
  AliCDBEntry* entry = man->Get("GRP/Calib/QAThresholds");
  if(!entry)    {     fIsOnlineThr = kFALSE;    return;    }
  
  TObjArray* branch = (TObjArray*) entry->GetObject();
  if(!branch ) {     fIsOnlineThr = kFALSE;     return;    }

  AliQAThresholds* thresholds = (AliQAThresholds*) branch->FindObject("HMP");  
  if(!thresholds) { fIsOnlineThr = kFALSE; return;   }
  else 
     fIsOnlineThr = kTRUE; 
  
  
  Int_t teCnt = 0;
  TString parName = "zero";
  while ( thresholds->GetThreshold(teCnt)) 
  {
     if(!((thresholds->GetThreshold(teCnt))->GetName())) return;
         
     parName = (thresholds->GetThreshold(teCnt))->GetName();
     
     if( parName.Contains("HmpNumberOfExcludedDDLthreshold") )
     {
       TParameter<int>* myParam = (TParameter<int>*) thresholds->GetThreshold(teCnt); 
       fHmpQaThr_NumberOfExcludedDDL = myParam->GetVal();
     }
   
     
     if( parName.Contains("HmpDataSizeLowerThreshold") ) 
     {
       TParameter<int>* myParam = (TParameter<int>*) thresholds->GetThreshold(teCnt); 
       fHmpQaThr_DataSizeLowerThreshold = myParam->GetVal();
     }
             
     if( parName.Contains("HmpDataSizeUpperThreshold") ) 
     {
       TParameter<int>* myParam = (TParameter<int>*) thresholds->GetThreshold(teCnt); 
       fHmpQaThr_DataSizeUpperThreshold = myParam->GetVal();
     }
          
     if( parName.Contains("HmpPadOccupancyLowerThreshold") ) 
     {
       TParameter<float>* myParam = (TParameter<float>*) thresholds->GetThreshold(teCnt); 
       fHmpQaThr_PadOccupancyLowerThreshold = myParam->GetVal();
     }
           
     if( parName.Contains("HmpPadOccupancyUpperThreshold") ) 
     {
       TParameter<float>* myParam = (TParameter<float>*) thresholds->GetThreshold(teCnt); 
       fHmpQaThr_PadOccupancyUpperThreshold = myParam->GetVal();
     }
                    
     if( parName.Contains("HmpSectorGainLossWarningThreshold") ) 
      {
       TParameter<int>* myParam = (TParameter<int>*) thresholds->GetThreshold(teCnt); 
       fHmpQaThr_SectorGainLossWarningThreshold = myParam->GetVal();
     } 
       
     if( parName.Contains("HmpSectorGainLossErrorThreshold") ) 
       {
       TParameter<int>* myParam = (TParameter<int>*) thresholds->GetThreshold(teCnt); 
       fHmpQaThr_SectorGainLossErrorThreshold = myParam->GetVal();
     } 
       
     if( parName.Contains("HmpMissingPadFractionWarningThreshold") ) 
       {
       TParameter<float>* myParam = (TParameter<float>*) thresholds->GetThreshold(teCnt); 
       fHmpQaThr_MissingPadFractionWarningThreshold = myParam->GetVal();
     }
     
     if( parName.Contains("HmpMissingPadFractionErrorThreshold") ) 
         {
       TParameter<float>* myParam = (TParameter<float>*) thresholds->GetThreshold(teCnt); 
       fHmpQaThr_MissingPadFractionErrorThreshold = myParam->GetVal();
     }
      
     teCnt++;    
  }//while
  
  
  // PrintThresholds();
}
//___________________________________________________________________________________________________
void AliHMPIDQAChecker::PrintThresholds()
{
  Printf("--- Printing thresholds ---");  
  Printf(Form("--- Default or online: %d ---",fIsOnlineThr));  
  Printf(Form("--- fHmpQaThr_NumberOfExcludedDDL %d ---",fHmpQaThr_NumberOfExcludedDDL));
  Printf(Form("--- fHmpQaThr_DataSizeLowerThreshold %d ---",fHmpQaThr_DataSizeLowerThreshold));
  Printf(Form("--- fHmpQaThr_DataSizeUpperThreshold %d ---",fHmpQaThr_DataSizeUpperThreshold));
  Printf(Form("--- fHmpQaThr_PadOccupancyLowerThreshold %f ---",fHmpQaThr_PadOccupancyLowerThreshold));
  Printf(Form("--- fHmpQaThr_PadOccupancyUpperThreshold %f ---",fHmpQaThr_PadOccupancyUpperThreshold));
  Printf(Form("--- fHmpQaThr_SectorGainLossWarningThreshold %d ---",fHmpQaThr_SectorGainLossWarningThreshold));
  Printf(Form("--- fHmpQaThr_SectorGainLossErrorThreshold %d ---",fHmpQaThr_SectorGainLossErrorThreshold));
  Printf(Form("--- fHmpQaThr_MissingPadFractionWarningThreshold %f ---",fHmpQaThr_MissingPadFractionWarningThreshold));
  Printf(Form("--- fHmpQaThr_MissingPadFractionErrorThreshold %f ---",fHmpQaThr_MissingPadFractionErrorThreshold));
  Printf("--- Printing thresholds done ---");  

  
  
  
}
//___________________________________________________________________________________________________


