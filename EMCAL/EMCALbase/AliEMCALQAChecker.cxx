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
#include <TClass.h>
#include <TH1.h> 
#include <TF1.h> 
#include <TH1I.h> 
#include <TH2F.h>
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <TLine.h>
#include <TText.h>
#include <TPaveText.h>
#include <TMath.h>

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliQAThresholds.h"
#include "AliEMCALQAChecker.h"

ClassImp(AliEMCALQAChecker)

//__________________________________________________________________
/// Constructor.
///
AliEMCALQAChecker::AliEMCALQAChecker() : 
AliQACheckerBase("EMCAL","EMCAL Quality Assurance Data Maker"),
fTextSM(new TText*[fgknSM]),
fLineCol(new TLine(47.5,-0.5,47.5,207.5)),
fText(new TPaveText(0.2,0.7,0.8,0.9,"NDC"))
{
  fLineCol->SetLineColor(1);
  fLineCol->SetLineWidth(2);
  
  for(int iSM = 0; iSM < fgknSM; iSM++) 
  {
    int iside = iSM % 2;
    int isect = iSM / 2;
    int isectNum = isect;
    
    if (isectNum > 5) isectNum += 3; // DCal
    
    if (isectNum < 5)
    {
      if (iside == 0) 
      { // A side
        fTextSM[iSM]= new TText(20, 8+24*isect, Form("SM A%d",isectNum) );
      }
      else 
      { // C side
        fTextSM[iSM]= new TText(64, 8+24*isect, Form("SM C%d",isectNum) );
      }
    }
    else if ( isectNum>4 && isectNum<12)
    {
      if (iside == 0) 
      { // A side
        if(isectNum ==5)
          fTextSM[iSM]= new TText(20, 8+24*(isect-1)+8+6, Form("SM A%d",isectNum) );
        else
          fTextSM[iSM]= new TText(20, 8+24*(isect-1)+8, Form("SM A%d",isectNum) );
      }
      else 
      { // C side
        if(isectNum ==5)
          fTextSM[iSM]= new TText(64, 8+24*(isect-1)+8+6, Form("SM C%d",isectNum) );
        else
          fTextSM[iSM]= new TText(64, 8+24*(isect-1)+8, Form("SM C%d",isectNum) );
      }
    }
    else
    {
      if (iside == 0) { // A side
        fTextSM[iSM]= new TText(20, 8+24*(isect-2)+16+6, Form("SM A%d",isectNum) );
      }
      else { // C side
        fTextSM[iSM]= new TText(64, 8+24*(isect-2)+16+6, Form("SM C%d",isectNum) );
      }
    }
  }
  
  for(int i = 0; i < fgknSectLines; i++) 
  {
    if(i<5)
    {
      fLineRow[i] = new TLine(-0.5,23.5+(24*i),95.5,23.5+(24*i));
      fLineRow[i]->SetLineColor(1);
      fLineRow[i]->SetLineWidth(2);
    }
    else if(i>4 && i<9){
      fLineRow[i] = new TLine(-0.5,23.5+(24*(i-1))+8,95.5,23.5+(24*(i-1))+8);
      fLineRow[i]->SetLineColor(1);
      fLineRow[i]->SetLineWidth(2);
    }
    else{
      fLineRow[i] = new TLine(-0.5,23.5+(24*(i-2))+16,95.5,23.5+(24*(i-2))+16);
      fLineRow[i]->SetLineColor(1);
      fLineRow[i]->SetLineWidth(2);
    }
  }
  
  for(int i = 0; i < 3; i++) {
    fTextL1[i] = new TPaveText(0.2,0.8,0.8,0.9,"NDC");
  }
  
}          

//__________________________________________________________________
/// Destructor.
///
AliEMCALQAChecker::~AliEMCALQAChecker() 
{
  delete [] fTextSM ;
  delete fLineCol ;
  for (Int_t i=0; i<5; ++i) delete fLineRow[i] ;
  delete fText  ; 
}

//______________________________________________________________________________
/// Check objects in list
///
void AliEMCALQAChecker::Check(Double_t * test, AliQAv1::ALITASK_t index, 
                              TObjArray ** list, const AliDetectorRecoParam * /*recoParam*/)
{
  if ( index == AliQAv1::kRAW ) 
  {
    CheckRaws(test, list);
    //printf ("checkers for task %d \n", index) ;		
  }
  
  if ( index == AliQAv1::kREC)
  {
    CheckRecPoints(test, list);
  }
  
  if ( index == AliQAv1::kESD )
  {
    CheckESD(test, list);
  }
  //AliWarning(Form("Checker for task %d not implement for the moment",index));
}

//______________________________________________________________________________
/// Get a given histo from the list.
///
TH1* AliEMCALQAChecker::GetHisto(TObjArray* list, const char* hname, Int_t specie) const
{
  TH1* h = static_cast<TH1*>(list->FindObject(Form("%s_%s",AliRecoParam::GetEventSpecieName(specie),hname)));
  
  if (!h)
  {
    AliError(Form("Did not find expected histo %s",hname));
  }
  
  return h;
}

//______________________________________________________________________________
/// Mark histo as originator of some QA error/warning.
///
Double_t AliEMCALQAChecker::MarkHisto(TH1& histo, Double_t value) const
{
  if ( value != 1.0 )
  {
    histo.SetBit(AliQAv1::GetQABit());
  }
  
  return value;
}

//______________________________________________________________________________
///  Check RAW QA histograms
///  adding new checking method: 25/04/2010, Yaxian Mao
///
///  Comparing the amplitude from current run to the reference run, if the ratio in the range [0.8, .12], count as a good tower.
///  If more than 90% towers are good, EMCAL works fine, otherwise experts should be contacted. 
///
void AliEMCALQAChecker::CheckRaws(Double_t * test, TObjArray ** list)
{
  // Setting the thresholds
  Float_t ratioThresh = 0.9;   // threshold for calibration ratio = good towers/all towers (default 0.9)
  Float_t threshG     = 0.5;   // threshold for L1 Gamma triggers (default 0.5)
  Float_t threshJ     = 0.5;   // threshold for L1 Jet triggers (default 0.5)
  Int_t badLinkThresh = 1;     // threshold for bad links (default 1)
  
  AliCDBManager* man = AliCDBManager::Instance();
  if(man)
  {
    AliCDBEntry* entry = man->Get("GRP/Calib/QAThresholds");
    if(entry)
    {
      TObjArray* branch = (TObjArray*) entry->GetObject();
      if(branch)
      {
        AliQAThresholds* thresholds = (AliQAThresholds*) branch->FindObject("EMC");
        if(thresholds)
        {
          TParameter<float>* paramR  = (TParameter<float>*) thresholds->GetThreshold(0);
          TParameter<float>* paramG  = (TParameter<float>*) thresholds->GetThreshold(1);
          TParameter<float>* paramJ  = (TParameter<float>*) thresholds->GetThreshold(2);
          TParameter<int>*   paramL  = (TParameter<int>*)   thresholds->GetThreshold(3);
          
          if(paramR)
            ratioThresh = paramR->GetVal();
          if(paramG)
            threshG = paramG->GetVal();
          if(paramJ)
            threshJ = paramJ->GetVal();
          if(paramL)
            badLinkThresh = paramL->GetVal();
        }
      }
    }
  }
  
  //cols*rows (in module units) * 4 (each module is 2x2 towers)
  Double_t nTot = AliEMCALTriggerMappingV2::fSTURegionNEta*AliEMCALTriggerMappingV2::fSTURegionNPhi*4;

  //subtracting towers from 6 TRUs (missing in DCAL) * 96 (modules in TRU) * 4 (each module is 2x2 towers) 
  nTot -= 6*AliEMCALTriggerMappingV2::fNModulesInTRU*4;
  
  TList *lstF = 0;
  Int_t calibSpecieId = (Int_t)TMath::Log2( AliRecoParam::kCalib );
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
  {
    test[specie] = 0.0 ; 
    
    if ( !AliQAv1::Instance()->IsEventSpecieSet(specie)) continue ; 
    
    if (list[specie]->GetEntries() == 0) 
    {
      test[specie] = 0. ; // nothing to check
    }
    else 
    {
      //get calib histos
      TH2F * hdata  = (TH2F*)list[specie]->At(k2DRatioAmp) ; 
      TH1F * ratio = (TH1F*)list[specie]->At(kRatioDist) ;
      
      //get L1 histos
      TH2F *hL1GammaPatch = (TH2F*)list[specie]->At(kGL1);
      TH2F *hL1JetPatch = (TH2F*)list[specie]->At(kJL1);
      TH1I *hFrameR = (TH1I*)list[specie]->At(kSTUTRU);
      
      //TH1I *hNTimeSamplesTRU = (TH1I*)list[specie]->At(kNL0TRUSamples);
      //TH1I *hRMSTimeforTRU = (TH1I*)list[specie]->At(kNL0TRURMS);
      //  =======================================================================================
      // calib histo checker first:
      if( hdata && ratio )
      {
        //  first clean lines, text (functions)
        lstF = hdata->GetListOfFunctions();
        CleanListOfFunctions(lstF);
        lstF = ratio->GetListOfFunctions();
        CleanListOfFunctions(lstF);
        
        if(hdata->GetEntries()!=0 && ratio->GetEntries()!=0) 
        {
          lstF = hdata->GetListOfFunctions();
          
          // adding the lines to distinguish different SMs
          lstF->Add(fLineCol->Clone()); 
          for(Int_t iLine = 0; iLine < fgknSectLines; iLine++) 
          {
            lstF->Add(fLineRow[iLine]->Clone());
          } 
          
          // Now adding the text to for each SM
          for(Int_t iSM = 0 ; iSM < fgknSM ; iSM++)
          {  // number of SMs loop start
            lstF->Add(fTextSM[iSM]->Clone()); 
          }			
                    
          // now check the ratio histogram
          lstF = ratio->GetListOfFunctions();
          
          Double_t binContent = 0. ;  
          Int_t nGoodTower = 0 ;
          Double_t rv = 0. ;
          for(Int_t ix = 1; ix <= hdata->GetNbinsX(); ix++) 
          {
            for(Int_t iy = 1; iy <= hdata->GetNbinsY(); iy++) 
            {
              binContent = hdata->GetBinContent(ix, iy) ; 
              //if (binContent < 1.2 && binContent > 0.8) nGoodTower++ ;
            }
          }
          
          //rv = nGoodTower/nTot ; 
          // printf("%2.2f %% towers out of range [0.8, 1.2]\n", (1-rv)*100);
          /*if(fText){
           lstF->Add(fText->Clone()) ;
           fText->Clear() ; 
           
           fText->AddText(Form("%2.2f %% towers out of range [0.8, 1.2]", (1-rv)*100));     
           if (rv < ratioThresh) {
           test[specie] = ratioThresh;
           //  2 lines text info for quality         
           fText->SetFillColor(2) ;
           fText->AddText(Form("EMCAL = NOK, CALL EXPERTS!!!")); 
           }
           else {
           test[specie] = 1 - ratioThresh;
           fText->SetFillColor(3) ;
           fText->AddText(Form("EMCAL = OK, ENJOY...")); 
           }
           }// fText*/
        }// calib histo checking done
      }// histograms NOT NULL
      
       //  ========================================================================================
       // now check L0 (NEW!!!)
      
      //lstF = hNTimeSamplesTRU->GetListOfFunctions();
      //CleanListOfFunctions(lstF);
      //lstF = hRMSTimeforTRU->GetListOfFunctions();
      //CleanListOfFunctions(lstF);
      
      //  ========================================================================================
      //  now L1 checks:
      
      if( hL1GammaPatch )
      {
        //  first clean lines, text (functions)
        lstF = hL1GammaPatch->GetListOfFunctions();
        CleanListOfFunctions(lstF);
        
        if (specie != calibSpecieId) 
        {
          // if(hL1GammaPatch->GetEntries() !=0 ) {
          if(hL1GammaPatch->GetEntries() > 10) 
          { //  need some statistics for hot spot calculation
            lstF = hL1GammaPatch->GetListOfFunctions();
            
            //  Checker for L1GammaPatch 
            // Double_t dL1GmeanTrig    = 1./2961.; 
            // Double_t dL1GmeanTrigTRU = 1./32.; 
            // Int_t sigmaG    = 100; //  deviation from mean value (increased to 100)
            // Int_t sigmaGTRU = 5; //  deviation from mean value for TRUs
            Double_t dL1GEntries = hL1GammaPatch->GetEntries();
            Int_t badL1G[2*AliEMCALTriggerMappingV2::fNPhi][2*AliEMCALTriggerMappingV2::fNEta]   = {{0}} ;
            Int_t badL1GTRU[2][AliEMCALTriggerMappingV2::fNTotalTRU/2] = {{0}} ;
            Int_t nBadL1G    = 0;
            Int_t nBadL1GTRU = 0;
            Double_t binContentTRU[2][AliEMCALTriggerMappingV2::fNTotalTRU/2] = {{0.}};
            for(Int_t ix = 1; ix <=  hL1GammaPatch->GetNbinsX(); ix++) 
            {
              for(Int_t iy = 1; iy <=  hL1GammaPatch->GetNbinsY(); iy++) 
              {
                Double_t binContent = hL1GammaPatch->GetBinContent(ix, iy) ; 
                if (binContent != 0) 
                {
                  //  fill counter for TRUs
                  // binContentTRU[(Int_t)((ix-1)/24)][(Int_t)((iy-1)/4)] += binContent;// OLD TRU SCHEME
                  binContentTRU[(Int_t)((ix-1)/8)][(Int_t)((iy-1)/12)] += binContent;  // NEW TRU SCHEME
          
                  //OLD METHOD (if a patch triggers > sigmaG * mean value (1/#patch positions total) says "hot spot !")
                  //  if ((double)binContent/(double)dL1GEntries > sigmaG*dL1GmeanTrig) {
                  //    badL1G[ix-1][iy-1] += 1;
                  //    nBadL1G += 1;
                  //  }
                  
                  //  NEW METHOD (if Rate > Threshold * ( (Number of towers or TRUs * Average rate) - Rate ) --> "hot spot !")
                  //  Thresold = how much does the noisy tower/TRU contribute to the rate
                  //             1.0 --> Rate of noisy tower/TRU = Rate of all other towers/TRUs  
                  if (binContent/dL1GEntries > threshG / ( 1 + threshG )) 
                  {
                    badL1G[ix-1][iy-1] += 1;
                    nBadL1G += 1;
                  }
                }
              }
            }
            
            //  check TRUs
            for(Int_t ix = 1; ix <=  2; ix++)
            {
              for(Int_t iy = 1; iy <= AliEMCALTriggerMappingV2::fNTotalTRU/2; iy++) 
              {
                if(binContentTRU[ix-1][iy-1]/dL1GEntries >  threshG / ( 1 + threshG )) 
                { 
                  badL1GTRU[ix-1][iy-1] += 1;
                  nBadL1GTRU += 1;
                }
              }
            }
            
            /* if(fTextL1[0]){
             lstF->Add(fTextL1[0]->Clone()) ;
             fTextL1[0]->Clear() ; 
             
             if (nBadL1G == 0 && nBadL1GTRU == 0 ) {
             fTextL1[0]->SetFillColor(3) ;
             fTextL1[0]->AddText(Form("L1 GAMMA TRIGGER = OK, ENJOY...")); 
             }
             else if (nBadL1G == 0){
             fTextL1[0]->SetFillColor(2) ;
             fTextL1[0]->AddText(Form("HOT SPOT IN L1 GAMMA TRIGGER (TRU) = CALL EXPERT!!"));
             
             }
             else{
             fTextL1[0]->SetFillColor(2) ;
             fTextL1[0]->AddText(Form("HOT SPOT IN L1 GAMMA TRIGGER = CALL EXPERT!!"));
             
             for(Int_t ix = 1; ix <=  hL1GammaPatch->GetNbinsX(); ix++) {
             for(Int_t iy = 1; iy <=  hL1GammaPatch->GetNbinsY(); iy++) {
             if(badL1G[ix-1][iy-1] != 0) printf("L1 Gamma patch with position x = %d, y = %d is out of range\n",ix,iy);
             }
             }
             */
            //}
            //}// fTextL1[0]
          }//  L1 gamma patch checking done
        }//  if (specie != calibSpecieId) ..
      }// hL1GammaPatch NOT NULL
      
      if( hL1JetPatch )
      {
        lstF = hL1JetPatch->GetListOfFunctions();
        CleanListOfFunctions(lstF);
        
        if (specie != calibSpecieId) 
        {
          // if(hL1JetPatch->GetEntries() !=0) {
          if(hL1JetPatch->GetEntries() > 10) 
          { //  need some statistics for hot spot calculation
            lstF = hL1JetPatch->GetListOfFunctions();
            
            //  Checker for L1JetPatch
            // Double_t dL1JmeanTrig = 1/126.;
            // Int_t sigmaJ = 5; //  deviation from  mean value
            Double_t dL1JEntries = hL1JetPatch->GetEntries();
            Int_t badL1J[12][16] = {{0}} ;// NEED TO CHECK THIS FOR JETs !!!!!!!
            Int_t nBadL1J = 0;
            for(Int_t ix = 1; ix <=  hL1JetPatch->GetNbinsX(); ix++) 
            {
              for(Int_t iy = 1; iy <=  hL1JetPatch->GetNbinsY(); iy++) 
              {
                Double_t binContent = hL1JetPatch->GetBinContent(ix, iy) ; 
                if (binContent != 0) 
                {
                  //  OLD METHOD  (if a patch triggers > sigmaJ * mean value (1/#patch positions total) says "hot spot !")
                  //  if ((double)binContent/(double)dL1JEntries > sigmaJ*dL1JmeanTrig) {
                  //  	badL1J[ix-1][iy-1] += 1 ;
                  //  	nBadL1J += 1;
                  //  }
                  
                  //  NEW METHOD (if Rate > Threshold * ( (Number of towers or TRUs * Average rate) - Rate ) --> "hot spot !")
                  //  Threshold: same definitionas for Gamma
                  if ((double)binContent/(double)dL1JEntries > threshJ / ( 1 + threshJ )) 
                  {
                    badL1J[ix-1][iy-1] += 1 ;
                    nBadL1J += 1;
                  }
                }
              }
            }
            
            /*
             if(fTextL1[1]){
             lstF->Add(fTextL1[1]->Clone()) ;
             fTextL1[1]->Clear() ; 
             
             if (nBadL1J == 0) {
             fTextL1[1]->SetFillColor(3) ;
             fTextL1[1]->AddText(Form("L1 JET TRIGGER = OK, ENJOY...")); 
             }
             else {
             fTextL1[1]->SetFillColor(2) ;
             fTextL1[1]->AddText(Form("HOT SPOT IN L1 JET TRIGGER = CALL EXPERT!!")); 
             
             for(Int_t ix = 1; ix <=  hL1JetPatch->GetNbinsX(); ix++) {
             for(Int_t iy = 1; iy <=  hL1JetPatch->GetNbinsY(); iy++) {
             if(badL1J[ix-1][iy-1] != 0) printf("L1 Jet patch with position x = %d, y = %d is out of range\n",(4*ix-4),(4*iy-4));
             }
             }
             */
            
            //}
            //}// fTextL1[1]
          } //  L1 Jet patch checking done
        } //  if (specie != calibSpecieId) ..
      }//  hL1JetPatch NOT NULL
      
      if(hFrameR)
      {
        lstF = hFrameR->GetListOfFunctions();
        CleanListOfFunctions(lstF);
        
        if(hFrameR->GetEntries() !=0) 
        {
          lstF = hFrameR->GetListOfFunctions();
          
          Int_t badLink[AliEMCALTriggerMappingV2::fNTotalTRU] = {0};
          Int_t nBadLink = 0;
          for(Int_t ix = 1; ix <= hFrameR->GetNbinsX(); ix++) 
          {
            Double_t binContent = hFrameR->GetBinContent(ix) ; 
            if (binContent == 0) {
              badLink[ix-1] += 1;
              nBadLink += 1;
            }
          }
          
          if(fTextL1[2])
          {
            lstF->Add(fTextL1[2]->Clone()) ;
            fTextL1[2]->Clear() ; 
            
            if (nBadLink < badLinkThresh)
            {
              fTextL1[2]->SetFillColor(3) ;
              fTextL1[2]->AddText(Form("LINK TRU-STU = OK, ENJOY...")); 
            }
            else 
            {
              fTextL1[2]->SetFillColor(2) ;
              fTextL1[2]->AddText(Form("PROBLEM WITH TRU-STU LINK = CALL EXPERT!!"));
              /*
               for(Int_t ix = 0; ix <= hFrameR->GetNbinsX(); ix++) {
               if(badLink[ix] != 0) printf("STU link with TRU %d is out\n",ix);
               }
               */
            }   
          }// fTextL1[2]
        } //  Checker for link TRU-STU done
      } // hFrameR NOT NULL
    } //  species processed		
  } //  specie
}

//______________________________________________________________________________
///
/// Initialises QA and QA checker settings.
///
void AliEMCALQAChecker::Init(const AliQAv1::DETECTORINDEX_t det) 
{
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

//______________________________________________________________________________
///
///  Clean up.
///
void AliEMCALQAChecker::CleanListOfFunctions(TList *list)
{ 
  if (list) 
  {
    TObject *stats = list->FindObject("stats"); list->Remove(stats);
    
    TObject *obj;
    while ( (obj = list->First()) ) 
    {
      while(list->Remove(obj)) { } delete obj; 
    }
    
    if (stats) list->Add(stats);
  }
  else 
  {
    AliWarning(Form("Checker : empty list of data functions; returning"));
    return;
  }
}


