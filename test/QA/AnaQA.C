/*
 *  AnaQA.C
 *  
 *
 *  Created by schutz on 29/09/08.
 *  Copyright 2008 CERN. All rights reserved.
 *
 */
#include <TFile.h>
#include <TList.h>
#include <TNamed.h>
#include "AliQA.h"

void AnaQA(Int_t run) 
{
// Macro to analyse the output of the QAChecker
  
  // Open the file that holds the AliQA object
  const char * resultFileName = "QA.root" ; //AliQA::GetQAResultFileName() ; 
  const char * msgE = "QA ERROR: " ; 
  const char * msgS = "QA SIGNAL: " ; 

  TFile * inputQAFile = TFile::Open(resultFileName) ; 
  if ( ! inputQAFile ) {
    printf("QA ERROR: File %s not found\n", AliQA::GetQAResultFileName()) ;
    exit(1) ; 
  }
  // Get the AliQA object from the file 
 // inputQAFile.ls() ; 
  AliQA * qa = dynamic_cast<AliQA*>(inputQAFile->Get(AliQA::GetQAName())) ; 
  // Show the status of all Detectors
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
    for (Int_t task = 0 ; task < AliQA::kNTASK ; task++) {
      if (qa->IsSetAny(AliQA::DETECTORINDEX_t(det), AliQA::ALITASK_t(task))) {
        qa->ShowStatus(AliQA::DETECTORINDEX_t(det), AliQA::ALITASK_t(task)) ;         
        // found a bad detector, open the QA data file and search for the faulty histogram
        TFile * dataQAFile = TFile::Open(AliQA::GetQADataFileName(AliQA::GetDetName(det), run)) ; 
        if ( ! dataQAFile ) {
          printf("%s File %s not found\n", msgE, AliQA::GetQADataFileName(AliQA::GetDetName(det), run)) ;
          exit(1) ; 
        }    
        dataQAFile->cd(AliQA::GetDetName(AliQA::DETECTORINDEX_t(det))) ; 
        TDirectory * saveDir = gDirectory ; 
        switch (task) {
          case AliQA::kNULLTASK:
            break ; 
          case AliQA::kRAW:
            Bool_t dir = saveDir->cd(AliQA::GetTaskName(AliQA::kRAWS)) ; 
            if ( ! dir ) {
              printf("%s Directory %s not found in %s\n", msgE, AliQA::GetTaskName(AliQA::kRAWS).Data(), AliQA::GetQADataFileName(AliQA::GetDetName(det), run)) ;
            } else {
              TList * listofkeys = gDirectory->GetListOfKeys() ; 
              for (Int_t key = 0 ; key < listofkeys->GetEntries() ; key++) {
                TNamed * obj = dynamic_cast<TNamed*>(listofkeys->At(key)) ; 
                if (obj) {
                  Bool_t rv = obj->TestBit(AliQA::GetQABit()) ;
                  if (rv)
                    printf("%s QA bit set in %s/%s/%s\n", 
                           msgS, 
                           AliQA::GetDetName(det), 
                           AliQA::GetTaskName(AliQA::kRAWS).Data(), 
                           obj->GetName()); 
                }
              }
            }
            break ;  
          case AliQA::kSIM:
            dir = saveDir->cd(AliQA::GetTaskName(AliQA::kHITS)) ; 
            if ( ! dir ) {
              printf("%s Directory %s not found in %s\n", msgE, AliQA::GetTaskName(AliQA::kHITS).Data(), AliQA::GetQADataFileName(AliQA::GetDetName(det), run)) ;
            } else {
              TList * listofkeys = gDirectory->GetListOfKeys() ; 
              for (Int_t key = 0 ; key < listofkeys->GetEntries() ; key++) {
                TNamed * obj = dynamic_cast<TNamed*>(listofkeys->At(key)) ; 
                if (obj) {
                  Bool_t rv = obj->TestBit(AliQA::GetQABit()) ;
                  if (rv)
                    printf("%s QA bit set in %s/%s/%s\n", 
                           msgS, 
                           AliQA::GetDetName(det), 
                           AliQA::GetTaskName(AliQA::kHITS).Data(), 
                           obj->GetName()); 
                }
              }
            }
            dir = saveDir->cd(AliQA::GetTaskName(AliQA::kSDIGITS)) ; 
            if ( ! dir ) {
              printf("%s Directory %s not found in %s\n", msgE, AliQA::GetTaskName(AliQA::kSDIGITS).Data(), AliQA::GetQADataFileName(AliQA::GetDetName(det), run)) ;
            } else {
              TList * listofkeys = gDirectory->GetListOfKeys() ; 
              for (Int_t key = 0 ; key < listofkeys->GetEntries() ; key++) {
                TNamed * obj = dynamic_cast<TNamed*>(listofkeys->At(key)) ; 
                if (obj) {
                  Bool_t rv = obj->TestBit(AliQA::GetQABit()) ;
                  if (rv)
                    printf("%s QA bit set in %s/%s/%s\n", 
                           msgS, 
                           AliQA::GetDetName(det), 
                           AliQA::GetTaskName(AliQA::kSDIGITS).Data(), 
                           obj->GetName()); 
                }
              }
            }
            
            dir = saveDir->cd(AliQA::GetTaskName(AliQA::kDIGITS)) ; 
            if ( ! dir ) {
              printf("%s Directory %s not found in %s\n", msgE, AliQA::GetTaskName(AliQA::kDIGITS).Data(), AliQA::GetQADataFileName(AliQA::GetDetName(det), run)) ;
            } else {
              TList * listofkeys = gDirectory->GetListOfKeys() ; 
              for (Int_t key = 0 ; key < listofkeys->GetEntries() ; key++) {
                TNamed * obj = dynamic_cast<TNamed*>(listofkeys->At(key)) ; 
                if (obj) {
                  Bool_t rv = obj->TestBit(AliQA::GetQABit()) ;
                  if (rv)
                    printf("%s QA bit set in %s/%s/%s\n", 
                           msgS, 
                           AliQA::GetDetName(det), 
                           AliQA::GetTaskName(AliQA::kDIGITS).Data(), 
                           obj->GetName()); 
                }
              }
            }
            break ;
          case AliQA::kREC:
            dir = saveDir->cd(AliQA::GetTaskName(AliQA::kRECPOINTS)) ; 
            if ( ! dir ) {
              printf("%s Directory %s not found in %s\n", msgE, AliQA::GetTaskName(AliQA::kRECPOINTS).Data(), AliQA::GetQADataFileName(AliQA::GetDetName(det), run)) ;
            } else {
              TList * listofkeys = gDirectory->GetListOfKeys() ; 
              for (Int_t key = 0 ; key < listofkeys->GetEntries() ; key++) {
                TNamed * obj = dynamic_cast<TNamed*>(listofkeys->At(key)) ; 
                if (obj) {
                  Bool_t rv = obj->TestBit(AliQA::GetQABit()) ;
                  if (rv)
                    printf("%s QA bit set in %s/%s/%s\n", 
                           msgS, 
                           AliQA::GetDetName(det), 
                           AliQA::GetTaskName(AliQA::kRECPOINTS).Data(), 
                           obj->GetName()) ; 
                }
              }
            }
            break ;
            case AliQA::kESD:
            dir = saveDir->cd(AliQA::GetTaskName(AliQA::kESDS)) ; 
            if ( ! dir ) {
              printf("%s Directory %s not found in %s\n", msgE, AliQA::GetTaskName(AliQA::kESDS).Data(), AliQA::GetQADataFileName(AliQA::GetDetName(det), run)) ;
            } else {
              TList * listofkeys = gDirectory->GetListOfKeys() ; 
              for (Int_t key = 0 ; key < listofkeys->GetEntries() ; key++) {
                TNamed * obj = dynamic_cast<TNamed*>(listofkeys->At(key)) ; 
                if (obj) {
                  Bool_t rv = obj->TestBit(AliQA::GetQABit()) ;
                  if (rv)
                    printf("%s QA bit set in %s/%s/%s\n", 
                           msgS, 
                           AliQA::GetDetName(det), 
                           AliQA::GetTaskName(AliQA::kESDS).Data(), 
                           obj->GetName()); 
                }
              }
            }
            break ;
            case AliQA::kANA:
            break ;
          default:
             break ;
        }
        dataQAFile->Close() ; 
      }
    }
  }
  inputQAFile->Close() ; 
}
