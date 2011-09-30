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


// --- ROOT system ---
#include <TH1F.h> 
#include <TIterator.h> 
#include <TString.h> 
#include <TPaveText.h> 
#include <TObjArray.h>
#include "TList.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliZDCQAChecker.h"

ClassImp(AliZDCQAChecker)

//____________________________________________________________________________
void AliZDCQAChecker::Check(Double_t *  test, AliQAv1::ALITASK_t index, TObjArray ** list,
      const AliDetectorRecoParam * /*recoParam*/) 
{
  // Checks the QA histograms on the input list: 
  //
  const char* taskName = AliQAv1::GetAliTaskName(index);
  //printf("\n\tAliZDCQAChecker -> checking QA histos for task %s\n",taskName);
  //
  Int_t ihitHisto=0, idigHisto=0;
  Int_t irecHisto=0, irawHisto=0, esdInd=0;

  TObjArray messages;
  messages.SetOwner(kTRUE);

  for(Int_t specie = 0; specie<AliRecoParam::kNSpecies; specie++){
    Int_t count = 0; 
    if(!AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(specie))) 
      continue ;
    //printf("\tAliZDCQAChecker -> specie %d, AliRecoParam::ConvertIndex(specie) %d, AliRecoParam::kLowMult %d, IsEventSpecieSet(specie) %d\n",
    //  specie, AliRecoParam::ConvertIndex(specie) ,AliRecoParam::kLowMult,
    //  AliQAv1::Instance(AliQAv1::GetDetIndex(GetName()))->IsEventSpecieSet(AliRecoParam::ConvertIndex(specie)));
    
    // ====================================================================
    // 	Checks for p-p events
    // ====================================================================
    if(AliRecoParam::ConvertIndex(specie) == AliRecoParam::kLowMult){
      if(list[specie]->GetEntries()==0){  
        AliWarning("\t The list to be checked is empty!"); // nothing to check
        return;
      }
      //AliDebug(AliQAv1::GetQADebugLevel(), Form("\n\tAliZDCQAChecker-> checking QA histograms for task %s\n\n",taskName));
      TIter next(list[specie]); 
      count = 0; 
      TH1 * hdata;	  
      
      Float_t res=0., percentageDiff=0.30;
      Float_t meanZNA=0., meanZNC=0., meanZPA=0., meanZPC=0.;
      Float_t pmCZNA=0., pmCZNC=0., pmCZPA=0., pmCZPC=0.;
      Float_t pmQZNA=0., pmQZNC=0., pmQZPA=0., pmQZPC=0.;
      Float_t sumADCZNA=0., sumADCZNC=0., sumADCZPA=0., sumADCZPC=0.;
      Float_t adcCZNA=0., adcCZNC=0., adcCZPA=0., adcCZPC=0.;
      Float_t adcQZNA=0., adcQZNC=0., adcQZPA=0., adcQZPC=0.;
      Int_t nentries = -99;
      
      while((hdata = dynamic_cast<TH1 *>(next()))){
        if(hdata){ 
          // -------------------------------------------------------------------
          if(index == AliQAv1::kSIM){
            //AliDebug(AliQAv1::GetQADebugLevel(), Form("\tAliZDCQAChecker-> checking histo %s",hdata->GetName()));
            // Check HITS histos
            //
	    if(!(strncmp(hdata->GetName(),"hHits",5))){
              if(hdata->GetEntries()>0){
	        if(ihitHisto==0)      meanZNC = hdata->GetMean();
		else if(ihitHisto==1) meanZNA = hdata->GetMean();
		else if(ihitHisto==2) meanZPC = hdata->GetMean();
		else if(ihitHisto==3) meanZPA = hdata->GetMean();
	        else if(ihitHisto==4) pmQZNC = hdata->GetMean();
	        else if(ihitHisto==5) pmQZNA = hdata->GetMean();
	        else if(ihitHisto==6) pmQZPC = hdata->GetMean();
	        else if(ihitHisto==7) pmQZPA = hdata->GetMean();
	        else if(ihitHisto==8)  pmCZNC = hdata->GetMean();
	        else if(ihitHisto==9)  pmCZNA = hdata->GetMean();
	        else if(ihitHisto==10) pmCZPC = hdata->GetMean();
	        else if(ihitHisto==11) pmCZPA = hdata->GetMean();
	      }
	      //
	      // --- Check whether (sum PMQi - PMC)/PMC < percentageDiff
	      if(ihitHisto==11){
	        if(TMath::Abs(meanZNC)>1.e-10){
                  if((TMath::Abs(pmQZNC-pmCZNC)/pmCZNC)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(meanZNA)>1.e-10){
                  if((TMath::Abs(pmQZNA-pmCZNA)/pmCZNA)<percentageDiff) 
                    res=1.;
                  else percentageDiff=
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(meanZPC)>1.e-10){
                  if((TMath::Abs(pmQZPC-pmCZPC)/pmCZPC)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(meanZPA)>1.e-10){
                  if((TMath::Abs(pmQZPA-pmCZPA)/pmCZPA)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	      }
	      ihitHisto++;
            }
            // Check DIGIT HIGH GAIN CHAIN histos
            else if(!(strncmp(hdata->GetName(),"hDig",4))){ 
              if(hdata->GetEntries()>0){
	        if(idigHisto==0)      sumADCZNC = hdata->GetMean();
		else if(idigHisto==1) sumADCZNA = hdata->GetMean();
		else if(idigHisto==2) sumADCZPC = hdata->GetMean();
		else if(idigHisto==3) sumADCZPA = hdata->GetMean();
	        else if(idigHisto==4) pmQZNC = hdata->GetMean();
	        else if(idigHisto==5) pmQZNA = hdata->GetMean();
	        else if(idigHisto==6) pmQZPC = hdata->GetMean();
	        else if(idigHisto==7) pmQZPA = hdata->GetMean();
	        else if(idigHisto==8)  pmCZNC = hdata->GetMean();
	        else if(idigHisto==9)  pmCZNA = hdata->GetMean();
	        else if(idigHisto==10) pmCZPC = hdata->GetMean();
	        else if(idigHisto==11) pmCZPA = hdata->GetMean();
	      }
	      //
	      // --- Check whether (sum PMQi - PMC)/PMC < percentageDiff
	      if(idigHisto==11){
	        if(TMath::Abs(sumADCZNC)>1.e-10){
                  if((TMath::Abs(adcQZNC-adcCZNC)/adcCZNC)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(sumADCZNA)>1.e-10){
                  if((TMath::Abs(adcQZNA-adcCZNA)/adcCZNA)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(sumADCZPC)>1.e-10){
                  if((TMath::Abs(adcQZPC-adcCZPC)/adcCZPC)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(sumADCZPA)>1.e-10){
                  if((TMath::Abs(adcQZPA-adcCZPA)/adcCZPA)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	      }
	      idigHisto++;	      
            }
          } 
          // -------------------------------------------------------------------
	  else if(index == AliQAv1::kRAW) {
	    //
            // Check RAW HIGH GAIN CHAIN histos
            
      	    messages.Clear();
	    
	    Bool_t iDetPM = kTRUE;
	    // --- Checks
	    if(irawHisto==16) nentries = Int_t (hdata->GetEntries());
	    if(irawHisto==18){ 
	      Float_t resADC=0.;
	      for(int ibin=1; ibin<=hdata->GetNbinsX(); ibin++){
		 if((hdata->GetBinContent(ibin))>10.){
		   res=1.;
		 }
		 else if((hdata->GetBinContent(ibin))<10.){
		   res=0.5;
		   if(ibin==1 || ibin==6 || ibin==11 || ibin==12 || ibin==13 || ibin==18){
		     iDetPM = kFALSE;
		   }
		 }
		 //
		 resADC += res;
            	 test[specie] += res;
            	 count++;
	      }
	      if(nentries != -99) messages.Add(new TObjString(Form("#entries %d",nentries)));
	      else messages.Add(new TObjString("#entries not known"));
	      //
	      Float_t rv=1.;
	      if(hdata->GetNbinsX() != 0) rv = resADC/hdata->GetNbinsX();
	      if(rv == 1.) messages.Add(new TObjString("ADCs are OK!")); 
	      else if(iDetPM==kFALSE){
	        messages.Add(new TObjString("Problem with some ADC!"));
                messages.Add(new TObjString("IF THIS IS A PHYSICS RUN"));
	      }
	      else if(iDetPM==kTRUE) messages.Add(new TObjString("Minor problem with ADCs"));
	      SetupHisto(messages, *hdata, rv);
	    }
	    else if(irawHisto==19){
	      // Reference values from RUN 145456
	      Double_t refTDCs[6] = {-322.5,-319.1,-320.9,-319.2,-319.7,-319.2};
	      Float_t resTDC=0.;
	      for(int ibin=1; ibin<=hdata->GetNbinsX(); ibin++){
		 if(TMath::Abs((hdata->GetBinContent(ibin))-refTDCs[ibin-1])<3.){
		   res=1.;
		 }
		 else if(TMath::Abs((hdata->GetBinContent(ibin))-refTDCs[ibin-1])<5.){
		   res=0.8;
            	 }
		 else res=0.5;
		 //
		 resTDC += res;
		 test[specie] += res;
            	 count++;
	      }
	      Float_t rv=1.;
	      if(hdata->GetNbinsX() != 0) rv = resTDC/hdata->GetNbinsX();
	      if(rv == 1.) messages.Add(new TObjString("TDCs are OK!")); 
	      else if(rv<1. && rv>0.75) messages.Add(new TObjString("Minor problem with TDCs"));
	      else{
	        messages.Add(new TObjString("Serious problem in ZDC timing"));
                messages.Add(new TObjString("IF THIS IS A PHYSICS RUN"));
	      }
	      SetupHisto(messages, *hdata, rv);
	    }
	    irawHisto++;
	    
          } 
          // -------------------------------------------------------------------
	  else if(index == AliQAv1::kREC) {
	    //
            // Check REC HIGH GAIN CHAIN histos
            if(hdata->GetEntries()>0){
	      if(irecHisto==0)       meanZNC = hdata->GetMean();
	      else if(irecHisto==1)  meanZNA = hdata->GetMean();
	      else if(irecHisto==2)  meanZPC = hdata->GetMean();
	      else if(irecHisto==3)  meanZPA = hdata->GetMean();
	      else if(irecHisto==4)  pmQZNC = hdata->GetMean();
	      else if(irecHisto==5)  pmQZNA = hdata->GetMean();
	      else if(irecHisto==6)  pmQZPC = hdata->GetMean();
	      else if(irecHisto==7)  pmQZPA = hdata->GetMean();
	      else if(irecHisto==8)  pmCZNC = hdata->GetMean();
	      else if(irecHisto==9)  pmCZNA = hdata->GetMean();
	      else if(irecHisto==10) pmCZPC = hdata->GetMean();
	      else if(irecHisto==11) pmCZPA = hdata->GetMean();
	    }
	    //
	    // --- Check whether (sum PMQi - PMC)/PMC < percentageDiff
	    if(irecHisto==11){
	      if(TMath::Abs(meanZNC)>1.e-10){
            	if((TMath::Abs(pmQZNC-pmCZNC)/pmCZNC)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if(TMath::Abs(meanZNA)>1.e-10){
            	if((TMath::Abs(pmQZNA-pmCZNA)/pmCZNA)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if(TMath::Abs(meanZPC)>1.e-10){
            	if((TMath::Abs(pmQZPC-pmCZPC)/pmCZPC)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if(TMath::Abs(meanZPA)>1.e-10){
            	if((TMath::Abs(pmQZPA-pmCZPA)/pmCZPA)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	    }
	    irecHisto++;	    
          } 
          // -------------------------------------------------------------------
	  else if(index == AliQAv1::kESD) {
	    //
            // Check ESD HIGH GAIN CHAIN histos
            if(hdata->GetEntries()>0){
	      if(esdInd==0)      sumADCZNC = hdata->GetMean();
	      else if(esdInd==1) sumADCZNA = hdata->GetMean();
	      else if(esdInd==2) sumADCZPC = hdata->GetMean();
	      else if(esdInd==3) sumADCZPA = hdata->GetMean();
	      else if(esdInd==6)  pmQZNC = hdata->GetMean();
	      else if(esdInd==7)  pmQZNA = hdata->GetMean();
	      else if(esdInd==8)  pmQZPC = hdata->GetMean();
	      else if(esdInd==9) pmQZPA = hdata->GetMean();
	      else if(esdInd==10) pmCZNC = hdata->GetMean();
	      else if(esdInd==11) pmCZNA = hdata->GetMean();
	      else if(esdInd==12) pmCZPC = hdata->GetMean();
	      else if(esdInd==13) pmCZPA = hdata->GetMean();
	    }
	    //
	    // --- Check whether (sum PMQi - PMC)/PMC < percentageDiff
	    if(esdInd==15){
	      if(TMath::Abs(sumADCZNC)>1.e-10){
            	if((TMath::Abs(pmQZNC-pmCZNC)/pmCZNC)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if(TMath::Abs(sumADCZNA)>1.e-10){
            	if((TMath::Abs(pmQZNA-pmCZNA)/pmCZNA)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if(TMath::Abs(sumADCZPC)>1.e-10){
            	if((TMath::Abs(pmQZPC-pmCZPC)/pmCZPC)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if(TMath::Abs(sumADCZPA)>1.e-10){
            	if((TMath::Abs(pmQZPA-pmCZPA)/pmCZPA)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
            }
            esdInd++;
         }  
	 else {
           AliWarning(Form("\n\t No ZDC QA for %s task\n",taskName)); 
           return ;
         }
        }//if(hdata) 
	else AliError("AliZDCQAChecker-> No histos!!!\n");
      }
    } // LowMult (p-p)
    // ====================================================================
    // 	Checks for A-A events
    // ====================================================================
    if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kHighMult) {
      if(list[specie]->GetEntries()==0){  
        AliWarning("\t The list to be checked is empty!");
        return ;
      }
      //AliDebug(AliQAv1::GetQADebugLevel(), Form("\n\tAliZDCQAChecker-> checking QA histograms for task %s\n\n",taskName));
      //
      TIter next(list[specie]); 
      count = 0; 
      TH1 * hdata;	  
      
      Float_t res=0., percentageDiff=0.10;
      Float_t meanZNA=0., meanZNC=0., meanZPA=0., meanZPC=0.;
      Float_t pmCZNA=0., pmCZNC=0., pmCZPA=0., pmCZPC=0.;
      Float_t pmQZNA=0., pmQZNC=0., pmQZPA=0., pmQZPC=0.;
      Float_t sumADCZNA=0., sumADCZNC=0., sumADCZPA=0., sumADCZPC=0.;
      Float_t adcCZNA=0., adcCZNC=0., adcCZPA=0., adcCZPC=0.;
      Float_t adcQZNA=0., adcQZNC=0., adcQZPA=0., adcQZPC=0.;
      Int_t nentries=-99;
      
      while((hdata = dynamic_cast<TH1 *>(next()))){
        if(hdata){ 
          //AliDebug(AliQAv1::GetQADebugLevel(), Form("\tAliZDCQAChecker-> checking histo %s",hdata->GetName()));
          // -------------------------------------------------------------------
          if(index == AliQAv1::kSIM){
            // Check HITS histos
            if (!(strncmp(hdata->GetName(),"hHits",5))){
              if(hdata->GetEntries()>0){
	        if(ihitHisto==0)      meanZNC = hdata->GetMean();
		else if(ihitHisto==1) meanZNA = hdata->GetMean();
		else if(ihitHisto==2) meanZPC = hdata->GetMean();
		else if(ihitHisto==3) meanZPA = hdata->GetMean();
	        else if(ihitHisto==4) pmQZNC = hdata->GetMean();
	        else if(ihitHisto==5) pmQZNA = hdata->GetMean();
	        else if(ihitHisto==6) pmQZPC = hdata->GetMean();
	        else if(ihitHisto==7) pmQZPA = hdata->GetMean();
	        else if(ihitHisto==8)  pmCZNC = hdata->GetMean();
	        else if(ihitHisto==9)  pmCZNA = hdata->GetMean();
	        else if(ihitHisto==10) pmCZPC = hdata->GetMean();
	        else if(ihitHisto==11) pmCZPA = hdata->GetMean();
	      }
	      //
	      // --- Check whether 2*|Mean ZNA - Mean ZNC|/(Mean ZNA + Mean ZNC) < percentageDiff
	      // --- and 2*|Mean ZPA - Mean ZPC|/(Mean ZPA + Mean ZPC) < 2*percentageDiff
	      if(ihitHisto==3){
	        if(TMath::Abs(meanZNC)>1.e-10 && TMath::Abs(meanZNA)>1.e-10){
                  if((2*TMath::Abs(meanZNC-meanZNA)/(meanZNA+meanZNC))<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(meanZPC)>1.e-10 && TMath::Abs(meanZPA)>1.e-10){
                  if((TMath::Abs(meanZPC-meanZPA)/(meanZPA+meanZPC))<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
              }
	      // --- Check whether (mean PMQi - PMC)/PMC < percentageDiff
	      if(ihitHisto==11){
	        if(TMath::Abs(meanZNC)>1.e-10){
                  if((TMath::Abs(pmQZNC-pmCZNC)/pmCZNC)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(meanZNA)>1.e-10){
                  if((TMath::Abs(pmQZNA-pmCZNA)/pmCZNA)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(meanZPC)>1.e-10){
                  if((TMath::Abs(pmQZPC-pmCZPC)/pmCZPC)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(meanZPA)>1.e-10){
                  if((TMath::Abs(pmQZPA-pmCZPA)/pmCZPA)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	      }
	      ihitHisto++;
            }
            // Check DIGITS histos
            else if (!(strncmp(hdata->GetName(),"hDig",4))){
              if(hdata->GetEntries()>0){
	        if(idigHisto==0)      sumADCZNC = hdata->GetMean();
		else if(idigHisto==1) sumADCZNA = hdata->GetMean();
		else if(idigHisto==2) sumADCZPC = hdata->GetMean();
		else if(idigHisto==3) sumADCZPA = hdata->GetMean();
	        else if(idigHisto==4) adcQZNC = hdata->GetMean();
	        else if(idigHisto==5) adcQZNA = hdata->GetMean();
	        else if(idigHisto==6) adcQZPC = hdata->GetMean();
	        else if(idigHisto==7) adcQZPA = hdata->GetMean();
	        else if(idigHisto==8)  adcCZNC = hdata->GetMean();
	        else if(idigHisto==9)  adcCZNA = hdata->GetMean();
	        else if(idigHisto==10) adcCZPC = hdata->GetMean();
	        else if(idigHisto==11) adcCZPA = hdata->GetMean();
	      }
	      //
	      // --- Check whether 2*|Mean ZNA - Mean ZNC|/(Mean ZNA + Mean ZNC) < percentageDiff
	      // --- and 2*|Mean ZPA - Mean ZPC|/(Mean ZPA + Mean ZPC) < 2*percentageDiff
	      if(idigHisto==3){
	        if(TMath::Abs(sumADCZNC)>1.e-10 && TMath::Abs(sumADCZNA)>1.e-10){
                  if((2*TMath::Abs(sumADCZNC-sumADCZNA)/(sumADCZNA+sumADCZNC))<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(sumADCZPC)>1.e-10 && TMath::Abs(sumADCZPA)>1.e-10){
                  if((TMath::Abs(sumADCZPC-sumADCZPA)/(sumADCZPA+sumADCZPC))<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
              }
	      // --- Check whether (sumADC PMQi - PMC)/PMC < percentageDiff
	      if(idigHisto==11){
	        if(TMath::Abs(sumADCZNC)>1.e-10){
                  if((TMath::Abs(adcQZNC-adcCZNC)/adcCZNC)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(sumADCZNA)>1.e-10){
                  if((TMath::Abs(adcQZNA-adcCZNA)/adcCZNA)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(sumADCZPC)>1.e-10){
                  if((TMath::Abs(adcQZPC-adcCZPC)/adcCZPC)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	        if(TMath::Abs(sumADCZPA)>1.e-10){
                  if((TMath::Abs(adcQZPA-adcCZPA)/adcCZPA)<percentageDiff) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  count++;
		}
	      }
              idigHisto++;
            }
          }
          // -------------------------------------------------------------------
          else if(index == AliQAv1::kRAW){
	    //
            // Check RAW HIGH GAIN CHAIN histos
            
      	    messages.Clear();
	    
	    Bool_t iDetPM = kTRUE;
	    // --- Checks
	    if(irawHisto==16) nentries = Int_t (hdata->GetEntries());
	    if(irawHisto==18){ 
	      Float_t resADC=0.;
	      for(int ibin=1; ibin<=hdata->GetNbinsX(); ibin++){
		 if((hdata->GetBinContent(ibin))>10.){
		   res=1.;
		 }
		 else if((hdata->GetBinContent(ibin))<10. && 
		   (ibin==1 || ibin==6 || ibin==11 || ibin==12 || ibin==13 || ibin==18)){
		   res=0.5;
		   iDetPM = kFALSE;
		 }
		 //
		 resADC += res;
            	 test[specie] += res;
            	 count++;
	      }
	      if(nentries != -99) messages.Add(new TObjString(Form("#entries %d",nentries)));
	      else messages.Add(new TObjString("#entries not known"));
	      //
	      Float_t rv=1.;
	      if(hdata->GetNbinsX() != 0) rv = resADC/hdata->GetNbinsX();
	      if(rv > 0.98) messages.Add(new TObjString("ADCs are OK!")); 
	      else if(iDetPM==kFALSE){
	        messages.Add(new TObjString("Problem with ADCs!"));
                messages.Add(new TObjString("IF THIS IS A PHYSICS RUN"));
	      }
	      else if(iDetPM==kTRUE) messages.Add(new TObjString("Minor problem with ADCs"));
	      SetupHisto(messages, *hdata, rv);
	    }
	    else if(irawHisto==19){
	      // Reference values from RUN 137161
	      Double_t refTDCs[6] = {-320.7,-319.0,-318.6,-319.9,-321.3,-320.8};
	      Float_t resTDC=0.;
	      for(int ibin=1; ibin<=hdata->GetNbinsX(); ibin++){
		 if(TMath::Abs((hdata->GetBinContent(ibin))-refTDCs[ibin-1])<3.){
		   res=1.;
		 }
		 else if(TMath::Abs((hdata->GetBinContent(ibin))-refTDCs[ibin-1])<4.){
		   res=0.8;
            	 }
		 else res=0.5;
		 //
		 resTDC += res;
		 test[specie] += res;
            	 count++;
	      }
	      Float_t rv=1.;
	      if(hdata->GetNbinsX() != 0) rv = resTDC/hdata->GetNbinsX();
	      if(rv == 1.) messages.Add(new TObjString("TDCs are OK!")); 
	      else if(rv<1 && rv>0.75) messages.Add(new TObjString("Minor problem with TDCs"));
	      else{
	        messages.Add(new TObjString("Serious problem in ZDC timing"));
                messages.Add(new TObjString("IF THIS IS A PHYSICS RUN"));
	      }
	      SetupHisto(messages, *hdata, rv);
	    }
	    irawHisto++;
	    
	  }   
          // -------------------------------------------------------------------
          else if(index == AliQAv1::kREC){
	    //
            // Check RAW HIGH GAIN CHAIN histos
            if(hdata->GetEntries()>0){
	      if(irecHisto==0)       meanZNC = hdata->GetMean();
	      else if(irecHisto==1)  meanZNA = hdata->GetMean();
	      else if(irecHisto==2)  meanZPC = hdata->GetMean();
	      else if(irecHisto==3)  meanZPA = hdata->GetMean();
	      else if(irecHisto==4)  pmQZNC = hdata->GetMean();
	      else if(irecHisto==5)  pmQZNA = hdata->GetMean();
	      else if(irecHisto==6)  pmQZPC = hdata->GetMean();
	      else if(irecHisto==7)  pmQZPA = hdata->GetMean();
	      else if(irecHisto==8)  pmCZNC = hdata->GetMean();
	      else if(irecHisto==9)  pmCZNA = hdata->GetMean();
	      else if(irecHisto==10) pmCZPC = hdata->GetMean();
	      else if(irecHisto==11) pmCZPA = hdata->GetMean();
	    }
            //
	    // --- Check whether 2*|Mean ZNA - Mean ZNC|/(Mean ZNA + Mean ZNC) < percentageDiff
	    // --- and 2*|Mean ZPA - Mean ZPC|/(Mean ZPA + Mean ZPC) < 2*percentageDiff
	    if(irecHisto==3){
	      if(TMath::Abs(meanZNC)>1.e-10 && TMath::Abs(meanZNA)>1.e-10){
            	if((2*TMath::Abs(meanZNC-meanZNA)/(meanZNA+meanZNC))<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if(TMath::Abs(meanZPC)>1.e-10 && TMath::Abs(meanZPA)>1.e-10){
            	if((TMath::Abs(meanZPC-meanZPA)/(meanZPA+meanZPC))<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
            }
	    // --- Check whether (sum PMQi - PMC)/PMC < percentageDiff
	    if(irecHisto==11){
	      if((TMath::Abs(meanZNC)>1.e-10) && (pmCZNC>1.e-10)){
            	if((TMath::Abs(pmQZNC-pmCZNC)/pmCZNC)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if((TMath::Abs(meanZNA)>1.e-10) && (pmCZNA>1.e-10)){
            	if((TMath::Abs(pmQZNA-pmCZNA)/pmCZNA)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if((TMath::Abs(meanZPC)>1.e-10) && (pmCZPC>1.e-10)){
            	if((TMath::Abs(pmQZPC-pmCZPC)/pmCZPC)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if((TMath::Abs(meanZPA)>1.e-10) && (pmCZPA>1.e-10)){
            	if((TMath::Abs(pmQZPA-pmCZPA)/pmCZPA)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	    }
	    irecHisto++;	 
	  }   
          // -------------------------------------------------------------------
          else if(index == AliQAv1::kESD){
	    //
            // Check ESD HIGH GAIN CHAIN histos
            if(hdata->GetEntries()>0){
	      if(esdInd==0)      sumADCZNC = hdata->GetMean();
	      else if(esdInd==1) sumADCZNA = hdata->GetMean();
	      else if(esdInd==2) sumADCZPC = hdata->GetMean();
	      else if(esdInd==3) sumADCZPA = hdata->GetMean();
	      else if(esdInd==6)  pmQZNC = hdata->GetMean();
	      else if(esdInd==7)  pmQZNA = hdata->GetMean();
	      else if(esdInd==8) pmQZPC = hdata->GetMean();
	      else if(esdInd==9) pmQZPA = hdata->GetMean();
	      else if(esdInd==10) pmCZNC = hdata->GetMean();
	      else if(esdInd==11) pmCZNA = hdata->GetMean();
	      else if(esdInd==12) pmCZPC = hdata->GetMean();
	      else if(esdInd==13) pmCZPA = hdata->GetMean();
	    }
	    //
	    // --- Check whether 2*|Mean ZNA - Mean ZNC|/(Mean ZNA + Mean ZNC) < percentageDiff
	    // --- and 2*|Mean ZPA - Mean ZPC|/(Mean ZPA + Mean ZPC) < 2*percentageDiff
	    if(esdInd==5){
	      if(TMath::Abs(sumADCZNC)>1.e-10 && TMath::Abs(sumADCZNA)>1.e-10){
            	if((2*TMath::Abs(sumADCZNC-sumADCZNA)/(sumADCZNA+sumADCZNC))<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if(TMath::Abs(sumADCZPC)>1.e-10 && TMath::Abs(sumADCZPA)>1.e-10){
            	if((TMath::Abs(sumADCZPC-sumADCZPA)/(sumADCZPA+sumADCZPC))<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
            }
	    // --- Check whether (sum PMQi - PMC)/PMC < percentageDiff
	    if(esdInd==15){
	      if(TMath::Abs(sumADCZNC)>1.e-10){
            	if((TMath::Abs(pmQZNC-pmCZNC)/pmCZNC)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if(TMath::Abs(sumADCZNA)>1.e-10){
            	if((TMath::Abs(pmQZNA-pmCZNA)/pmCZNA)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
/*	      if(TMath::Abs(sumADCZPC)>1.e-10){
            	if((TMath::Abs(pmQZPC-pmCZPC)/pmCZPC)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
	      if(TMath::Abs(sumADCZPA)>1.e-10){
            	if((TMath::Abs(pmQZPA-pmCZPA)/pmCZPA)<percentageDiff) 
            	  res=1.;
            	else 
            	  res=.5;
            	test[specie] += res;
            	count++;
	      }
*/
            }
            esdInd++;
          }  
	  else {
            AliWarning(Form("\n\t No ZDC QA for %s task\n",taskName)); 
            return ;
          } 
        }//if(hdata) 
	else AliError("\t  No histos found for ZDC!!!\n");
      }
    } // HighMult (Pb-Pb) 
    // ====================================================================
    // 	Checks for Calibration events
    // ====================================================================
    if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCalib) {
      AliWarning(Form("\n\t No check implemented in ZDC QA for %s task in CALIBRATION events\n",taskName)); 
      return ;
    } // Calibration
    // ====================================================================
    // 	Checks for cosmic events
    // ====================================================================
    else if (AliRecoParam::ConvertIndex(specie) == AliRecoParam::kCosmic) {
      AliWarning(Form("\n\t No check needed in ZDC QA for %s task in COSMIC events\n",taskName)); 
      return ; 
    } // Cosmic
    if(TMath::Abs(count)>1.e-10) test[specie] = test[specie]/count;
    AliDebug(AliQAv1::GetQADebugLevel(), Form("\n\t ZDC QA check result = %1.2f\n",test[specie]));
  } // Loop on species
}  

//___________________________________________________________________ 
void AliZDCQAChecker::SetupHisto(const TObjArray& messages, TH1& histo, Float_t& code)
{
  //
  /// Add text to histos
  //

  TPaveText* text = new TPaveText(0.70,0.70,0.99,0.99,"NDC");
    
  TIter next(&messages);
  TObjString* str;
    
  while ( ( str = static_cast<TObjString*>(next()) ) ){
    text->AddText(str->String());
  }

  TString defaultText = "";

  Int_t color = 0;
  if(code==1.){
    color = kGreen;
    defaultText = "Everything is fine!";
  }  
  else if(code<1. && code>=0.85){  
    color = kYellow;
    defaultText = "To be monitored in next runs";
  }
  else if(code<0.85 && code>=0.6){
    color = kOrange;
    defaultText = "notify the expert DURING THE DAY!";
  }
  else if(code<0.6){
    color = kRed;
    defaultText = "PLEASE CALL THE EXPERT!!!!";
  }


  text->AddText(defaultText.Data());
  text->SetFillColor(color);
                      
  //histo.SetFillStyle(1001);
  //histo.SetFillColor(color);

  histo.SetStats(kFALSE);
    
  histo.GetListOfFunctions()->Add(text);
}
