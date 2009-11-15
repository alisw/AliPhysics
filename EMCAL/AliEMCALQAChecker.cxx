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

/*
  Checks the quality assurance. 
  By comparing with reference data

  Based on PHOS code written by
  Y. Schutz CERN July 2007

 For the moment we only implement the checking of raw data QA.
 The checked for ESD and RecPoints will be implemented later.
 

 */


// --- ROOT system ---
#include <TClass.h>
#include <TH1.h> 
#include <TF1.h> 
#include <TH1I.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <TLine.h>
#include <TPaveText.h>
#include <TMath.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliEMCALQAChecker.h"
#include "AliEMCALQADataMakerRec.h"

ClassImp(AliEMCALQAChecker)

//__________________________________________________________________
AliEMCALQAChecker::AliEMCALQAChecker() : 
AliQACheckerBase("EMCAL","EMCAL Quality Assurance Data Maker"),
fLine(new TLine*[fknSM]),
fHref(new TLine*[fknSM]),
fText(NULL)
{
	/// ctor
  for (Int_t sm = 0 ; sm < fknSM ; sm++){
    fLine[sm] = NULL ; 
    fHref[sm] = NULL ; 
  }
}          

//__________________________________________________________________
AliEMCALQAChecker::~AliEMCALQAChecker() 
{
	/// dtor
  delete [] fLine ; 
  delete [] fHref ;
  if (fText) 
    delete fText ; 
}

//__________________________________________________________________
AliEMCALQAChecker::AliEMCALQAChecker(const AliEMCALQAChecker& qac) : 
AliQACheckerBase(qac.GetName(), qac.GetTitle()), 
fLine(new TLine*[fknSM]),
fHref(new TLine*[fknSM]),
fText(qac.fText)
{
	/// copy ctor 
  for (Int_t sm = 0 ; sm < fknSM ; sm++){
    fLine[sm] = qac.fLine[sm] ; 
    fHref[sm] = qac.fHref[sm] ; 
  }
}   
//__________________________________________________________________
AliEMCALQAChecker& AliEMCALQAChecker::operator = (const AliEMCALQAChecker &qac)
{
  fText = qac.fText;

  for (Int_t sm = 0 ; sm < fknSM ; sm++){
    fLine[sm] = qac.fLine[sm] ; 
    fHref[sm] = qac.fHref[sm] ; 
  }
  return *this ;
}

//______________________________________________________________________________
Double_t *
AliEMCALQAChecker::Check(AliQAv1::ALITASK_t index, TObjArray ** list, AliDetectorRecoParam * /*recoParam*/)
{
	/// Check objects in list
	
	if ( index == AliQAv1::kRAW ) 
	{
		return CheckRaws(list);
		printf ("checkers for task %d \n", index) ;		
	}
	
	if ( index == AliQAv1::kREC)
	{
		return CheckRecPoints(list);
	}
	
	if ( index == AliQAv1::kESD )
	{
		return CheckESD(list);
	}
	AliWarning(Form("Checker for task %d not implement for the moment",index));
	return NULL;
}

//______________________________________________________________________________
TH1* 
AliEMCALQAChecker::GetHisto(TObjArray* list, const char* hname, Int_t specie) const
{
	/// Get a given histo from the list
	TH1* h = static_cast<TH1*>(list->FindObject(Form("%s_%s",AliRecoParam::GetEventSpecieName(specie),hname)));
	if (!h)
	{
		AliError(Form("Did not find expected histo %s",hname));
	}
	return h;
}

//______________________________________________________________________________
Double_t 
AliEMCALQAChecker::MarkHisto(TH1& histo, Double_t value) const
{
	/// Mark histo as originator of some QA error/warning
	
	if ( value != 1.0 )
	{
		histo.SetBit(AliQAv1::GetQABit());
	}
	
	return value;
}


//______________________________________________________________________________
Double_t *
AliEMCALQAChecker::CheckRaws(TObjArray ** list)
{
//  Check RAW QA histograms	
//  We count the times of the response for each tower, the propability for all towers should be the same (average is given value).
//  We skip the first few cycles since the statistics is not enough, the average should be always larger than 1 at least.
//  By calculating the difference between the counts for each tower and the average, we decide whether we should recalculate 
//  the average depending the the gaus fitting on the divation distribution. 
//  During the recalutation of the average, we count how many towers are used for the calculation.
//  From the fraction of towers used, we decide whether each SM works fine or not
//  From the divation of average, we set the QA flag for the full detetcor as INFO, WARNING, ERROR or FATAL.
  
//  -- Yaxian Mao, CCNU/CERN/LPSC
					
  Float_t kThreshold = 80. ; 
  
  Double_t * test = new Double_t[AliRecoParam::kNSpecies] ;
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    test[specie] = 1.0 ; 
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie)) 
      continue ; 
    if (list[specie]->GetEntries() == 0)  
      test[specie] = 0. ; // nothing to check
    else {
      TH1 * hdata = (TH1*)list[specie]->At(kTowerHG) ; 
      if(hdata->GetEntries()==0)
	continue;
      Int_t  nbin   = hdata->GetNbinsX() ;
      Int_t  nTower = nbin/4 ;  // ! number of channels for each SM
      if(fText) {
        fText->DeleteText() ; 
        fText->Clear() ;
      }
      else {
        fText = new TPaveText(0.3,0.6,0.7,0.9,"NDC") ;
      }
      fText->AddText(Form("OK if more than %2.2f %% inside aver-sigma < HG counts < aver+sigma", kThreshold));
      //TPaveText * fText[fknSM] = {0};
      Double_t rv    = 0. ;   //value to define the flag for the full detector
      for(Int_t iSM = 0 ; iSM < fknSM ; iSM++){  //number of SMs loop start
        TString projname = Form("proj_%d",iSM);
        Double_t aver  = 0. ;  //average calculated by the input histogram
        Double_t recal = 0 ;   // recalculate the average if the divation is large
        Double_t init  = 0. ;   // value to decide whether we should recalculate the average
        Double_t mean  = 0. ;  //the divation from the average, from the gaus fitting peak 
        Double_t width = 0. ;   // the sigma of the devation, from the gaus fitting
        Int_t flag     = 0 ;  //counts for channels used for recalculation of average
        Double_t ratio = 0. ;  //value to decide whether each SM works
        TH1F * proj    = NULL  ;  //a temp histogram store the difference from the average for the fitting procedure       
       
       for(Int_t iTower = iSM*nTower ; iTower<(iSM+1)*nTower ; iTower++)  //channel loop to calculate the average
          aver += hdata->GetBinContent(iTower);
        
        aver /=nTower;
        AliInfo(Form("SM: %d has average = %f\n",iSM, aver));
       Double_t ymin = hdata->GetBinContent(hdata->GetMinimumBin());
       Double_t ymax = hdata->GetBinContent(hdata->GetMaximumBin());
        proj = new TH1F(projname,projname,nbin,-aver,aver);
       fLine[iSM] = dynamic_cast<TLine*>(hdata->GetListOfFunctions()->FindObject(fLine[iSM]));
        //initialize the lines separate different SM
       if (!fLine[iSM]) {
          fLine[iSM] = new TLine((iSM+1)*nTower,ymin,(iSM+1)*nTower,ymax);
          fLine[iSM]->SetLineColor(1);
          fLine[iSM]->SetLineWidth(2);
          hdata->GetListOfFunctions()->Add(fLine[iSM]);  
          list[specie]->AddAt(hdata, kTowerHG) ; 
        }          
        else {
          fLine[iSM]->SetX1((iSM+1)*nTower) ; 
          fLine[iSM]->SetY1(ymin) ; 
          fLine[iSM]->SetX2((iSM+1)*nTower) ; 
          fLine[iSM]->SetY2(ymax) ; 
        }
        //filling the histogram with the difference and fitting with gaus function to obtain mean and sigma
       for(Int_t iTower = iSM*nTower ; iTower<(iSM+1)*nTower ; iTower++){
	 proj->Fill(hdata->GetBinContent(iTower)-aver);
       }	
       proj->Fit("gaus","","QNO");
       mean=proj->GetFunction("gaus")->GetParameter(1); // ! mean should be peaked at 0 in principal
       width=proj->GetFunction("gaus")->GetParameter(2);
       AliInfo(Form("aver = %f, mean = %f, sigma =%f\n",aver,mean, width));
        
       if(aver) init=TMath::Abs(mean/aver);  //calculate the divation from the average 
       
       //if mean or sigma is too huge or the sigma is too small, the fitting failed 
        if((aver+mean) < 1. || width/aver >2. || width < 1.e-3) 
          flag = -1 ;
       else {   //recalculate the average if the fitting didn't give the initial average  
	       aver+=mean;  //average corrected after fitting is fitting works
	     //if(aver <= 1. ) break ;
	    // if divation is too large, we conside to recalate the average by excluding channels too far from the average
         while(init>0.01 && aver >= 1.){
          if(flag) flag = 0 ;    //for each time recalculation, reset flag 
           for(Int_t iTower = iSM*nTower ; iTower < (iSM+1)*nTower ; iTower++){
             if(hdata->GetBinContent(iTower)>(aver-width) && hdata->GetBinContent(iTower)<(aver+width)){
               flag++ ;
	             recal += hdata->GetBinContent(iTower);
             }  
           }  //end of channels loop 
           
           if(flag == 0 || recal < 1.) {
             flag = -1 ;
             break ;
           }
           else {
             recal/=flag ;
	           init=(aver-recal)/(aver+recal) ; //difference between the new average and the pervious one
            aver =(aver+recal)/2 ; //recalculate the average by the last two values  
           }  
         } //out of while condition
          
         ratio=100.0*flag/nTower ;  //counting how many towers used for average recalculation
	       AliInfo(Form("SM %d has %2.2f %% chanhel works \n", iSM, ratio));
       }  // end of recalculation
        
        rv+=init ; //define the deviation from the final average
    
        //define the average line on each SM
       fHref[iSM] = dynamic_cast<TLine*>(hdata->GetListOfFunctions()->FindObject(fHref[iSM]));
       if(!fHref[iSM]) {
         fHref[iSM] = new TLine(iSM*nTower,aver,(iSM+1)*nTower,aver);
	       hdata->GetListOfFunctions()->Add(fHref[iSM]);  
	       list[specie]->AddAt(hdata, kTowerHG) ; 
       }          
       else {
         fHref[iSM]->Clear() ; 
	       fHref[iSM]->SetX1(iSM*nTower) ; 
         fHref[iSM]->SetY1(aver) ; 
         fHref[iSM]->SetX2((iSM+1)*nTower) ; 
         fHref[iSM]->SetY2(aver) ; 
       }
        
        hdata->Paint() ; 
 
        // if channels used for average recalculation smaller than the threshold,
        // then too much noise channels or channels didn't work 
        if(ratio<= kThreshold && flag >0){                             
          //fText->SetFillColor(2);
          fHref[iSM]->SetLineColor(2);
	        fHref[iSM]->SetLineWidth(2);
	        fText->SetFillColor(2);
	        fText->AddText(Form("SM %d: NOK = %2.0f %% channels OK!!!",iSM,ratio));
	        //fText[iSM]->AddText(Form("SM %d NOT OK, only %5.2f %% works!!!",iSM,flag/nTower));
        }         
        else if (flag == -1 || flag == 0 ) {  //fitting failed or recalculation didn't call
          fText->SetFillColor(2);
          fHref[iSM]->SetLineColor(2);
          fHref[iSM]->SetLineWidth(2);
          fText->SetFillColor(2);
          fText->AddText(Form("SM %d: NOK Fitting failed",iSM,ratio));
	        //fText[iSM]->AddText(Form("SM %d has %5.2f %% towers wok normally",iSM,flag/nTower));
        }  
        else {
          fText->SetFillColor(3);
	        fHref[iSM]->SetLineColor(3);
	        fHref[iSM]->SetLineWidth(2);
	        fText->AddText(Form("SM %d: OK = %2.0f %% channels OK",iSM,ratio));
	        //fText[iSM]->AddText(Form("SM %d has %5.2f %% towers wok normally",iSM,flag/nTower));
        }
        hdata->Paint() ; 
       //hdata->GetListOfFunctions()->Add(fText[iSM]);
       delete proj ; 
      }  // end of SM loop
      rv/=fknSM;
      hdata->GetListOfFunctions()->Add(fText);
      hdata->Paint() ;  
      if ( rv <=0.1 ) {
        AliInfo(Form("The deviation rv = %2.2f is small compared to average, detector works fine",rv));
        test[specie] =  0.05;
      }
      
      if ( rv <=0.5 && rv >0.1 )  {
        AliWarning(Form("The deviation rv = %2.2f is acceptable from average,  the detector works not perfect!",rv));
        test[specie] =  0.3;
      }
      
      if ( rv <=0.8 && rv >0.5 ) {
        AliError(Form("Got a large deviation of %2.2f from average, some error on the detector!!!",rv));
        test[specie] =  0.7;
      }
      
      if ( rv >0.8 ) {
        AliError(Form("Got too large deviation of %2.2f from average, detector out of control???!!!",rv));
        test[specie] =  0.9; 
      }
    
    } //end of checking raw
    
  }  //species loop 
  return test ;
}

//______________________________________________________________________________
void AliEMCALQAChecker::Init(const AliQAv1::DETECTORINDEX_t det) 
{
	/// intialises QA and QA checker settings
	AliQAv1::Instance(det) ; 
	Float_t hiValue[AliQAv1::kNBIT] ; 
	Float_t lowValue[AliQAv1::kNBIT] ;
	lowValue[AliQAv1::kINFO]      = 0.0   ; 
	hiValue[AliQAv1::kINFO]       = 0.1 ; 
	lowValue[AliQAv1::kWARNING]   = 0.1 ; 
  hiValue[AliQAv1::kWARNING]    = 0.5 ; 
	lowValue[AliQAv1::kERROR]     = 0.5   ; 
	hiValue[AliQAv1::kERROR]      = 0.8 ; 
	lowValue[AliQAv1::kFATAL]     = 0.8   ; 
	hiValue[AliQAv1::kFATAL]      = 1.0 ; 
	SetHiLo(&hiValue[0], &lowValue[0]) ; 
}

