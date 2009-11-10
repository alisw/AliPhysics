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

 For the moment we only implement the checking of raw data QA 
 by looking at the average at the each tower.
 We count the times of the response for each tower, the propability for all towers should be the same (given value) 
 depending on that value, the resulting QA flag is info, warning, error or fatal. 
 Also we check the percentage of towers inside the average for each SM
 

 -- Yaxian Mao, CCNU/CERN/LPSC
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
	/// Check raws
  
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
      Int_t  nTower = nbin/4 ;
      if(fText) {
        fText->DeleteText() ; 
        fText->Clear() ;
      }
      else {
        fText = new TPaveText(0.3,0.6,0.7,0.9,"NDC") ;
      }
      fText->AddText(Form("OK if more than %2.2f %% inside aver-sigma < HG counts < aver+sigma", kThreshold));
			//TPaveText * fText[fknSM] = {0};
      Double_t rv    = 0. ;
			for(Int_t iSM = 0 ; iSM < fknSM ; iSM++){
        TString projname = Form("proj_%d",iSM);
        Double_t aver  = 0. ;
        Double_t recal = 0 ;
        Double_t init  = 0. ;
        Double_t mean  = 0. ;
        Double_t width = 0. ;
        Int_t flag     = 0 ;
        Double_t ratio = 0. ;
        TH1F * proj    = NULL  ;         
				for(Int_t iTower = iSM*nTower ; iTower<(iSM+1)*nTower ; iTower++){
					aver += hdata->GetBinContent(iTower);
					//printf("Tower: %d has counts = %f\n",iTower, hdata->GetBinContent(iTower));
				}
				aver /=nTower;
				AliInfo(Form("SM: %d has average = %f\n",iSM, aver));
				Double_t ymin = hdata->GetBinContent(hdata->GetMinimumBin());
				Double_t ymax = hdata->GetBinContent(hdata->GetMaximumBin());
        proj = new TH1F(projname,projname,nbin,-aver,aver);
        
				fLine[iSM] = dynamic_cast<TLine*>(hdata->GetListOfFunctions()->FindObject(fLine[iSM]));
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
        	
				//Double_t size = 0.24 ;
				//fText[iSM] = new TPaveText(0.1+size*iSM,0.7,size*(iSM+1),0.9,"NDC");
				for(Int_t iTower = iSM*nTower ; iTower<(iSM+1)*nTower ; iTower++){
					proj->Fill(hdata->GetBinContent(iTower)-aver);
				}	
				proj->Fit("gaus","","QNO");
				mean=proj->GetFunction("gaus")->GetParameter(1);
				width=proj->GetFunction("gaus")->GetParameter(2);
        AliInfo(Form("aver = %f, mean = %f, sigma =%f\n",aver,mean, width));
        //if mean or sigma is too huge from the fitting, fitting failed 
        if(aver < 0 || width/aver >2.) flag = -1 ;
        else {   //recalculate the average if the fitting didn't give the initial average          
          aver+=mean;
				  init=TMath::Abs(mean/aver); 
           while(init>0.01){
            if(flag) flag = 0 ;
            for(Int_t iTower = iSM*nTower ; iTower < (iSM+1)*nTower ; iTower++){
				        if(hdata->GetBinContent(iTower)>(aver-width) && hdata->GetBinContent(iTower)<(aver+width)){
					      flag++ ;
					      recal += hdata->GetBinContent(iTower);
                }
				    }
            recal/=flag;
            aver =(aver+recal)/2 ;  
            init=(aver-recal)/(aver+recal) ;	
          } //recalculate the average 
          if(flag)ratio=100.0*flag/nTower ;  //counting how many towers used for average recalculation
          rv+=TMath::Abs((aver-recal)/(aver+recal)) ; 
          AliInfo(Form("SM %d has %2.2f %% chanhel works \n", iSM, ratio));
        }
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
				if(ratio<= kThreshold && flag >0){ //if towers used for average recalculation smaller than 10%, then the problem on this SM
                            //fText->SetFillColor(2);
          fHref[iSM]->SetLineColor(2);
          fHref[iSM]->SetLineWidth(2);
          fText->SetFillColor(2);
					fText->AddText(Form("SM %d: NOK = %2.0f %% channels OK!!!",iSM,ratio));
					//fText[iSM]->AddText(Form("SM %d NOT OK, only %5.2f %% works!!!",iSM,flag/nTower));
				}
				else if (flag == -1 || flag == 0 || ratio == 0.) {
          fText->SetFillColor(2);
          fHref[iSM]->SetLineColor(2);
          fHref[iSM]->SetLineWidth(2);
          fText->SetFillColor(2);
					fText->AddText(Form("SM %d: NOK Fitting failed",iSM,ratio));
				    //fText[iSM]->AddText(Form("SM %d has %5.2f %% towers wok normally",iSM,flag/nTower));
				}	else {
          fText->SetFillColor(3);
          fHref[iSM]->SetLineColor(3);
          fHref[iSM]->SetLineWidth(2);
					fText->AddText(Form("SM %d: OK = %2.0f %% channels OK",iSM,ratio));
				    //fText[iSM]->AddText(Form("SM %d has %5.2f %% towers wok normally",iSM,flag/nTower));
				}
        hdata->Paint() ; 
			//hdata->GetListOfFunctions()->Add(fText[iSM]);
        delete proj ; 
      }
			rv/=fknSM;
			hdata->GetListOfFunctions()->Add(fText);
      hdata->Paint() ; 
			if ( rv <=0.1 ) 
			{
				AliInfo(Form("Got a small deviation rv = %f from average, SM works fine",rv));
				test[specie] =  0.05;
			}
			
			if ( rv <=0.5 && rv >0.1 ) 
			{
				AliWarning(Form("Got a deviation rv = %f from average, be careful!",rv));
				test[specie] =  0.3;
			}
			
			if ( rv <=0.5 && rv >0.8 ) 
			{
				AliError(Form("Got a large deviation of %f from average for SMs!!!",rv));
				test[specie] =  0.7;
			}
			if ( rv >0.8 ) 
			{
				AliError(Form("Got too large deviation of %f from average !!!???",rv));
				test[specie] =  0.9;
			}
		}
				
	}
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
	hiValue[AliQAv1::kWARNING]    = 0.5; 
	lowValue[AliQAv1::kWARNING]   = 0.1 ; 
	lowValue[AliQAv1::kERROR]     = 0.5   ; 
	hiValue[AliQAv1::kERROR]      = 0.8 ; 
	lowValue[AliQAv1::kFATAL]     = 0.8   ; 
	hiValue[AliQAv1::kFATAL]      = 1.0 ; 
	SetHiLo(&hiValue[0], &lowValue[0]) ; 
}

