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

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliZDCQAChecker.h"

ClassImp(AliZDCQAChecker)

//____________________________________________________________________________
Double_t * AliZDCQAChecker::Check(AliQAv1::ALITASK_t index, TObjArray ** list) 
{
  // Checks the QA histograms on the input list: 
  //
  Double_t * test   = new Double_t[AliRecoParam::kNSpecies] ;
  Int_t *    ntests = new Int_t[AliRecoParam::kNSpecies]  ; 
  const char* taskName = AliQAv1::GetAliTaskName(index);
  //
  
  //YS Int_t beamType=0; // 0 -> protons, 1 -> ions
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    test[specie] = 1.0 ; 
    ntests[specie] = 0 ; 
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
      continue ; 
    // ====================================================================
    // 	Checks for p-p events
    // ====================================================================
    //if (beamType==0){
    if ( specie == AliRecoParam::kLowMult) {
      if(list[specie]->GetEntries()==0){  
        AliWarning("\tAliZDCQAChecker->The list to be checked is empty!"); // nothing to check
        return test;
      }
      //AliDebug(AliQAv1::GetQADebugLevel(), Form("\n\tAliZDCQAChecker-> checking QA histograms for task %s\n\n",taskName));
      TIter next(list[specie]); 
      ntests[specie] = 0; 
      TH1 * hdata;	  
      while((hdata = dynamic_cast<TH1 *>(next()))){
        if(hdata){ 
          // -------------------------------------------------------------------
          if(index == AliQAv1::kSIM){
            //AliDebug(AliQAv1::GetQADebugLevel(), Form("\tAliZDCQAChecker-> checking histo %s",hdata->GetName()));
            // Check DIGITS histos
            if(!(strncmp(hdata->GetName(),"hDig",4))){
              if(hdata->GetEntries()>0){
                test[specie] += 1.; 
                ntests[specie]++;
              }
            }
            // Check HITS histos
            else{ 
              if(hdata->GetEntries()>0){
                test[specie] += 1.; 
                ntests[specie]++;
              }
            }
            // -------------------------------------------------------------------
          } else if(index == AliQAv1::kRAW) {
            if(hdata->GetEntries()!=0){
              if(hdata->GetMean()>10.) 
                test[specie] += 1.; 
              else 
                test[specie] = 0.5; 
              ntests[specie]++;
            }
            // -------------------------------------------------------------------
          } else if(index == AliQAv1::kESD) {
            Int_t    esdInd=0;
            if(hdata->GetEntries()!=0){
              if(esdInd>1){
                if(hdata->GetMean()>10.) 
                  test[specie] += 1.; 
                else 
                  test[specie] = 0.5; 
                ntests[specie]++;
              }
            }
            esdInd++;
            // -------------------------------------------------------------------
         } else {
           AliWarning(Form("\n\t No ZDC QA for %s task\n",taskName)); 
           return NULL;
         }
        } else
            AliError("AliZDCQAChecker-> No histos!!!\n");
      }
    }
    // ====================================================================
    // 	Checks for A-A events
    // ====================================================================
    //else if (beamType==1){ 
    else if ( specie == AliRecoParam::kHighMult) {
      if(list[specie]->GetEntries()==0){  
        AliWarning("\tAliZDCQAChecker->The list to be checked is empty!");
        return test;
      }
      //AliDebug(AliQAv1::GetQADebugLevel(), Form("\n\tAliZDCQAChecker-> checking QA histograms for task %s\n\n",taskName));
      
      TIter next(list[specie]); 
      TH1 * hdata;	  
      while((hdata = dynamic_cast<TH1 *>(next()))){
        if(hdata){ 
          //AliDebug(AliQAv1::GetQADebugLevel(), Form("\tAliZDCQAChecker-> checking histo %s",hdata->GetName()));
          ntests[specie] = 0; 
          Double_t meanX=0., meanY=0.;
          Double_t meanZNA=0., rmsZNA=0., meanZNC=0.;
          Double_t meanZPA=0., rmsZPA=0., meanZPC=0.;
          Double_t pmCZNA=0., pmCZNC=0., pmCZPA=0., pmCZPC=0.;
          Double_t pmQZNA=0., pmQZNC=0., pmQZPA=0., pmQZPC=0.;
          Float_t  res=0.;
          Int_t    testgood=0;
          // -------------------------------------------------------------------
          if(index == AliQAv1::kSIM){
            Int_t    digInd=0;
            // Check DIGITS histos
            if (!(strncmp(hdata->GetName(),"hDig",4))){
              // [1] check response of ZNC vs. ZNA
              if(digInd==0 || digInd==1){
                if(digInd==0){
                  if(hdata->GetEntries() != 0.){
                    testgood=1;
                    meanZNA = hdata->GetMean();
                    rmsZNA = hdata->GetRMS();
                  }
                }
                else{
                  if(hdata->GetEntries() != 0.){
                    testgood=1;
                    meanZNC = hdata->GetMean();
                  }
                  else testgood=0;
                  // check if the response m.v. of ZNA and ZNC are equal (@ 1sigma level)
                  if(testgood==1){
                    if(TMath::Abs(meanZNA-meanZNC)<rmsZNA) 
                      res=1.;
                    else 
                      res=.5;
                    testgood=0;
                    test[specie] += res;
                    ntests[specie]++;
                    //AliDebug(AliQAv1::GetQADebugLevel(), Form("\t %d performed tests, results %1.2f\n",ntests[specie],test[specie]/ntests[specie]));
                  }
                  else res=0.;
                }
              }
              // [2] check response of ZPC vs. ZPA
              else if(digInd==2 || digInd==3){
                if(digInd==2){
                  if(hdata->GetEntries() != 0.){
                    testgood=1;
                    meanZPA = hdata->GetMean();
                    rmsZPA = hdata->GetRMS();
                  }
                }
                else{
                  if(hdata->GetEntries() != 0.){
                    testgood=1;
                    meanZPC = hdata->GetMean();
                  }
               // check if the response m.v. of ZPA and ZPC are equal (@ 3sigma level)
               if(testgood==1){
                 if(TMath::Abs(meanZPA-meanZPC)<(3.*rmsZPA)) 
                   res=1.;
                 else 
                   res=.5;
                 test[specie] += res;
                 ntests[specie]++;
                 testgood=0;
                 //AliDebug(AliQAv1::GetQADebugLevel(), Form("\t %d performed tests, results %1.2f\n",ntests[specie],test[specie]/ntests[specie]);
               }
               else res=0.;
                }
              }
              // [2] check PMC responses vs. summed PMQ responses
              else if(digInd>3 && digInd<12){
                if(digInd==4) pmQZNC = hdata->GetMean();
                else if(digInd==5) pmQZNA = hdata->GetMean();
                else if(digInd==6) pmQZPC = hdata->GetMean();
                else if(digInd==7) pmQZPA = hdata->GetMean();
                else if(digInd==8){
                  pmCZNC = hdata->GetMean();
                  if(TMath::Abs(pmQZNC-pmCZNC)<(0.1*(pmQZNC+pmCZNC)/2)) 
                    res=1.;
                  else 
                    res=0.5;
                  test[specie] += res;
                  ntests[specie]++;
                  //AliDebug(AliQAv1::GetQADebugLevel(), Form("\t %d performed tests, results %1.2f\n",ntests[specie],test[specie]/ntests[specie]));
                }
                else if(digInd==9){
                  pmCZNA = hdata->GetMean();
                  if(TMath::Abs(pmQZNA-pmCZNA)<(0.1*(pmQZNA+pmCZNA)/2)) 
                    res=1.;
                  else 
                    res=0.5;
                  test[specie] += res;
                  ntests[specie]++;
                  //AliDebug(AliQAv1::GetQADebugLevel(), Form("\t %d performed tests, results %1.2f\n",ntests[specie],test[specie]/ntests[specie]));
                }
                else if(digInd==10){
                  pmCZPC = hdata->GetMean();
                  if(TMath::Abs(pmQZPC-pmCZPC)<(0.1*(pmQZPC+pmCZPC)/2)) 
                    res=1.;
                  else 
                    res=0.5;
                  test[specie] += res;
                  ntests[specie]++;
                  //AliDebug(AliQAv1::GetQADebugLevel(), Form("\t %d performed tests, results %1.2f\n",ntests[specie],test[specie]/ntests[specie]);
                }
                else if(digInd==11){
                  pmCZPA = hdata->GetMean();
                  if(TMath::Abs(pmQZPA-pmCZPA)<(0.1*(pmQZPA+pmCZPA)/2)) 
                    res=1.;
                  else 
                    res=0.5;
                  test[specie] += res;
                  ntests[specie]++;
                  //AliDebug(AliQAv1::GetQADebugLevel(), Form("\t %d performed tests, results %1.2f\n",ntests[specie],test[specie]/ntests[specie]));
                }
              }
              digInd++;
            }
            // Check HITS histos
            else{ 
              // hits histos
              meanX = hdata->GetMean(1);
              meanY = hdata->GetMean(2);
              // check if the spot is centered
              if((TMath::Abs(meanX)<0.2) && (TMath::Abs(meanY)<0.2)) 
                res=1.;
              else 
                res=0.5;
              test[specie] += res;
              ntests[specie]++;
              //AliDebug(AliQAv1::GetQADebugLevel(), Form("\t %d performed tests, results %1.2f\n",ntests[specie],test[specie]/ntests[specie]);
            }
          }
          // -------------------------------------------------------------------
          else if(index == AliQAv1::kRAW){
            Int_t    rawInd=0;
            // [1] check response of ZNC vs. ZNA
            if(rawInd==0 || rawInd==1){
              if(rawInd==0){
                if(hdata->GetEntries() != 0.){
                  testgood=1;
                  meanZNA = hdata->GetMean();
                  rmsZNA = hdata->GetRMS();
                }
              }
              else{
                if(hdata->GetEntries() != 0.){
                  testgood=1;
                  meanZNC = hdata->GetMean();
                }
                else testgood=0;
                // check if the response m.v. of ZNA and ZNC are equal (@ 1sigma level)
                if(testgood==1){
                  if(TMath::Abs(meanZNA-meanZNC)<rmsZNA) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  ntests[specie]++;
                  testgood=0;
                }
                else 
                  res=0.;
              }
            }
            // [2] check response of ZPC vs. ZPA
            else if(rawInd==2 || rawInd==3){
              if(rawInd==2){
                if(hdata->GetEntries() != 0.){
                  testgood=1;
                  meanZPA = hdata->GetMean();
                  rmsZPA = hdata->GetRMS();
                }
              }
              else{
                if(hdata->GetEntries() != 0.){
                  testgood=1;
                  meanZPC = hdata->GetMean();
                }
                // check if the response m.v. of ZPA and ZPC are equal (@ 3sigma level)
                if(testgood==1){
                  if(TMath::Abs(meanZPA-meanZPC)<(3.*rmsZPA)) 
                    res=1.;
                  else 
                    res=.5;
                  test[specie] += res;
                  ntests[specie]++;
                  testgood=0;
                }
                else 
                  res=0.;
              }
            }
            // [2] check PMC responses vs. summed PMQ responses
            else if(rawInd>3 && rawInd<12){
              if(rawInd==4) pmQZNC = hdata->GetMean();
              else if(rawInd==5) pmQZNA = hdata->GetMean();
              else if(rawInd==6) pmQZPC = hdata->GetMean();
              else if(rawInd==7) pmQZPA = hdata->GetMean();
              else if(rawInd==8){
                pmCZNC = hdata->GetMean();
                if(TMath::Abs(pmQZNC-pmCZNC)<(0.1*(pmQZNC+pmCZNC)/2)) 
                  res=1.;
                else 
                  res=0.5;
                test[specie] += res;
                ntests[specie]++;
              }
              else if(rawInd==9){
                pmCZNA = hdata->GetMean();
                if(TMath::Abs(pmQZNA-pmCZNA)<(0.1*(pmQZNA+pmCZNA)/2)) 
                  res=1.;
                else 
                  res=0.5;
                test[specie] += res;
                ntests[specie]++;
              }
              else if(rawInd==10){
                pmCZPC = hdata->GetMean();
                if(TMath::Abs(pmQZPC-pmCZPC)<(0.1*(pmQZPC+pmCZPC)/2)) 
                  res=1.;
                else 
                  res=0.5;
                test[specie] += res;
                ntests[specie]++;
              }
              else if(rawInd==11){
                pmCZPA = hdata->GetMean();
                if(TMath::Abs(pmQZPA-pmCZPA)<(0.1*(pmQZPA+pmCZPA)/2)) 
                  res=1.;
                else 
                  res=0.5;
                test[specie] += res;
                ntests[specie]++;
              }
            }
            rawInd++;
          }
          // -------------------------------------------------------------------
          else if(index == AliQAv1::kESD){
            Double_t eneQZNC, eneQZNA ;  
            Double_t eneQZPC, eneQZPA ;  
            Double_t eneCZNC, eneCZNA ;  
            Double_t eneCZPC, eneCZPA ;            
            Int_t esdInd=0;
            if(esdInd<2){
              // hits histos
              meanX = hdata->GetMean(1);
              meanY = hdata->GetMean(2);
              // check if the spot is centered
              if((TMath::Abs(meanX)<0.2) && (TMath::Abs(meanY)<0.2)) 
                res=1.;
              else 
                res=0.5;
              test[specie] += res;
              ntests[specie]++;
            }
            //
            else{
              // [1] check response of ZNC vs. ZNA
              if(esdInd==0 || esdInd==1){
                if(esdInd==0){
                  if(hdata->GetEntries() != 0.){
                    testgood=1;
                    meanZNA = hdata->GetMean();
                    rmsZNA = hdata->GetRMS();
                  }
                }
                else{
                  if(hdata->GetEntries() != 0.){
                    testgood=1;
                    meanZNC = hdata->GetMean();
                  }
                  else testgood=0;
                  // check if the response m.v. of ZNA and ZNC are equal (@ 1sigma level)
                  if(testgood==1){
                    if(TMath::Abs(meanZNA-meanZNC)<rmsZNA) 
                      res=1.;
                    else 
                      res=.5;
                    testgood=0;
                    test[specie] += res;
                    ntests[specie]++;
                  }
                  else res=0.;
                }
              }
              // [2] check response of ZPC vs. ZPA
              else if(esdInd==2 || esdInd==3){
                if(esdInd==2){
                  if(hdata->GetEntries() != 0.){
                    testgood=1;
                    meanZPA = hdata->GetMean();
                    rmsZPA = hdata->GetRMS();
                  }
                }
                else{
                  if(hdata->GetEntries() != 0.){
                    testgood=1;
                    meanZPC = hdata->GetMean();
                  }
                  // check if the response m.v. of ZPA and ZPC are equal (@ 3sigma level)
                  if(testgood==1){
                    if(TMath::Abs(meanZPA-meanZPC)<(3.*rmsZPA)) 
                      res=1.;
                    else 
                      res=.5;
                    test[specie] += res;
                    ntests[specie]++;
                    testgood=0;
                  }
                  else res=0.;
                }
              }
              // [2] check eneC responses vs. summed eneQ responses
              else if(esdInd>3 && esdInd<12){
                if(esdInd==4) eneQZNC = hdata->GetMean();
                else if(esdInd==5) eneQZNA = hdata->GetMean();
                else if(esdInd==6) eneQZPC = hdata->GetMean();
                else if(esdInd==7) eneQZPA = hdata->GetMean();
                else if(esdInd==8){
                  eneCZNC = hdata->GetMean();
                  if(TMath::Abs(eneQZNC-eneCZNC)<(0.1*(eneQZNC+eneCZNC)/2)) 
                    res=1.;
                  else 
                    res=0.5;
                  test[specie] += res;
                  ntests[specie]++;
                }
                else if(esdInd==9){
                  eneCZNA = hdata->GetMean();
                  if(TMath::Abs(eneQZNA-eneCZNA)<(0.1*(eneQZNA+eneCZNA)/2)) 
                    res=1.;
                  else 
                    res=0.5;
                  test[specie] += res;
                  ntests[specie]++;
                }
                else if(esdInd==10){
                  eneCZPC = hdata->GetMean();
                  if(TMath::Abs(eneQZPC-eneCZPC)<(0.1*(eneQZPC+eneCZPC)/2)) 
                    res=1.;
                  else 
                    res=0.5;
                  test[specie] += res;
                  ntests[specie]++;
                }
                else if(esdInd==11){
                  eneCZPA = hdata->GetMean();
                  if(TMath::Abs(eneQZPA-eneCZPA)<(0.1*(eneQZPA+eneCZPA)/2)) 
                    res=1.;
                  else 
                    res=0.5;
                  test[specie] += res;
                  ntests[specie]++;
                }
              }
              esdInd++;
            }
            //AliDebug(AliQAv1::GetQADebugLevel(), Form("\t %d performed tests, results %1.2f\n",ntests[specie],test[specie]/ntests[specie]));
          }
          else {
          AliWarning(Form("\n\t No ZDC QA for %s task\n",taskName)); 
          return NULL;
          }
        } else
            AliError("\t AliZDCQAChecker->No histos!!!\n");
        }
    } else {
      AliError(Form("Checking not implemented for %s, %s", 
                    AliRecoParam::GetEventSpecieName(AliRecoParam::kCosmic), 
                    AliRecoParam::GetEventSpecieName(AliRecoParam::kCalib))) ; 
    }
  }
   
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    if(ntests[specie]!=0) test[specie] = test[specie]/ntests[specie];
        AliDebug(AliQAv1::GetQADebugLevel(), Form("\n\tAliZDCQAChecker-> QA check result = %1.2f\n",test[specie]));
  }
  delete [] ntests ; 
  return test; 
}  

