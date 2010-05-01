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
//...
//  Checks the quality assurance. 
//  By comparing with reference data
//  Skeleton for T0
//---------------------------------------------
//checkig without reference data:
//for RAW QA all histograms should have approximatly the same 
//number of entries as RefPoint
//for Rec Points checks 
//  - amplitude measured by 2 methos
//  - online and offline T0 measurements
// for ESD quality of reconstruction ( and measurements):
// RMS of vertex and T0 less than 75ps
//
// Alla.Maevskaya@cern.ch   
//...

// --- ROOT system ---
#include <Riostream.h>
#include <TClass.h>
#include <TH1F.h> 
#include <TH2.h> 
#include <TIterator.h> 
#include <TKey.h> 
#include <TFile.h> 
#include <TMath.h>
#include <TString.h>

// --- Standard library ---

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliQAv1.h"
#include "AliQAChecker.h"
#include "AliCDBEntry.h"
#include "AliQAManager.h"
#include "AliT0QAChecker.h"

ClassImp(AliT0QAChecker)

//__________________________________________________________________
void AliT0QAChecker::Check(Double_t *  test, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * /*recoParam*/)
{

  // Super-basic check on the QA histograms on the input list:
  // look whether they are empty!
    
  char * detOCDBDir = Form("T0/%s/%s", AliQAv1::GetRefOCDBDirName(), AliQAv1::GetRefDataDirName()) ; 

  AliCDBEntry *QARefRec = AliQAManager::QAManager()->Get(detOCDBDir);
  //  QARefRec->Dump();
  if( !QARefRec){
    AliInfo("QA reference data NOT retrieved for Reconstruction check. No T0 reference distribution");
  }
  
   for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) 
    test[specie]    = 10.0 ; 

  Double_t nent[500];
  TString hname[500];
  const char *cname;
  memset(nent,0,500*sizeof(Double_t));
  Double_t w[500];
  memset(w,1,500*sizeof(Double_t));
  TH2 *fhRecDiff[3];  
  TH1 *fhRawEff[500];
  TH2 *fhRawTime[500];
  TH1 *fhESD[2];
  TH2 *fhHits[10];
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    //  TString dataType = AliQAv1::GetAliTaskName(index);
    if (list[specie]->GetEntries() == 0){
      test[specie] = 1. ; // nothing to check
    }
    else {
      TH1 * hdata ;
      TH2 * h ;
      for (Int_t ir=0; ir<list[specie]->GetEntries(); ir++) {
      
	//hits

	if(index == AliQAv1::kSIM && AliQAv1::GetTaskName(AliQAv1::kHITS)){
	  h =  (TH2*) list[specie]->UncheckedAt(ir);
	  cname = h->GetName();
	  hname[ir] = cname;
	  fhHits[ir] = h;
	}

	//raw
      if(index == AliQAv1::kRAW ){
        if(ir > 204 && ir<208 ) {
          hdata = (TH1*) list[specie]->UncheckedAt(ir);
          if(hdata) {
            cname = hdata->GetName();
            hname[ir] = cname;
            fhRawEff[ir] = hdata;
          }
        }
	if((ir>207 && ir<210)  && specie == AliRecoParam::kCalib) {
	  h =  (TH2*) list[specie]->UncheckedAt(ir);
	  if(h) {
	    cname = h->GetName();
	    hname[ir] = cname;
	    fhRawTime[ir] = h;
          AliDebug(AliQAv1::GetQADebugLevel(), Form("count %i %s ",ir, hname[ir].Data()) );
	  }
	}
      }
      //rec
      if(index == AliQAv1::kREC){
        h = (TH2*) list[specie]->UncheckedAt(ir);
        if(h) {
          cname = h->GetName();
          hname[ir] = cname;
           fhRecDiff[ir] = h;
        }
      }
      //esd
      if(index ==  AliQAv1::kESD){
        hdata = (TH1*) list[specie]->UncheckedAt(ir);
        if(hdata){
          fhESD[ir] = hdata;
          AliDebug(AliQAv1::GetQADebugLevel(), Form("count %i %s ",ir, hname[ir].Data()) );
        }
      }
    }



    //raw data

      if (index == AliQAv1::kRAW && specie == AliRecoParam::kCalib){
	test[specie] = CheckRaw(list[specie],dynamic_cast<TObjArray*>(dynamic_cast<TList *>(QARefRec->GetObject())->First()));
      }
      
      if(index == AliQAv1::kREC){
        //rec points
        for (Int_t icase=0; icase<2; icase++) {
          for (Int_t idet=0; idet<24; idet++) {
            Double_t mean = fhRecDiff[icase]->
            ProjectionY(Form("%s_py", fhRecDiff[icase]->GetName()),
                        idet,idet+1)->GetMean();
            Double_t rms= fhRecDiff[icase]->
            ProjectionY(Form("%s_py", fhRecDiff[icase]->GetName()),
                        idet,idet+1)->GetRMS();
            AliDebug(AliQAv1::GetQADebugLevel(), Form("name %s icase %i idet %i mean %f, rms %f\n",
                         fhRecDiff[icase]->GetName(), icase, idet, mean,rms)); 
	  	  
            if(TMath::Abs(mean) >1.5 || rms >1){
              AliDebug(AliQAv1::GetQADebugLevel(), Form(" calibration is nor perfect; test=%f", test)) ;
              test[specie]=0.25;
            }
            if(mean>3 || rms >5) {
              test[specie] = 0.1;
              AliDebug(AliQAv1::GetQADebugLevel(), Form(" wrong calibration test=%f", test[specie])) ;
            } 
          }
        }	 
      }
      if (index == AliQAv1::kESD) {
        //ESD
        for (Int_t icase=0; icase<2; icase++) {
          Double_t rmsVertex = fhESD[icase]->GetRMS();
          Double_t meanVertex = fhESD[icase]->GetMean();
          test[specie]=1;
          AliDebug(AliQAv1::GetQADebugLevel(), Form("numentries %d meanVertex %f rmsVertex %f", fhESD[icase]->GetEntries(), meanVertex, rmsVertex));
          if (TMath::Abs(rmsVertex)>3) {
            test[specie]=0.25;
            AliDebug(AliQAv1::GetQADebugLevel(), Form("Vertex position resolution not good  , rms= %f test=%f",
                            rmsVertex, test[specie])) ; 
          }
          if (TMath::Abs(meanVertex)>3) {
            test[specie]=0.25;
            AliDebug(AliQAv1::GetQADebugLevel(), Form("Vertex position bad calibrated  , Mean= %f test=%f",
                            meanVertex, test[specie])) ; 
          }
        }
      }
    } //  if (list->GetEntries() != 0
    AliDebug(AliQAv1::GetQADebugLevel(), Form("Test Result = %f", test[specie])) ;
  } 
}

//--------------------------------------------------------------------------
Double_t AliT0QAChecker::CheckRaw(TObjArray *listrec , TObjArray *listref) const
{
  
    //TH1 *fhRawEff;
    //TH1 *fhRawRef;
  TH2 *fhRawRec2d;
  TH2 *fhTime;

  TIter next(listref) ;
    //Int_t counter=0;
  Float_t refmean[50][25]; 
  Float_t refrms[50][25]; 
  Float_t checkr = 0;
  
    //Int_t nref = listref->GetEntries(); 
    //Int_t nrec = listrec->GetEntries(); 
  
  for (Int_t iii=4; iii<6; iii++){
    fhRawRec2d =(TH2*) listref->At(iii); 
    for (Int_t idet=1; idet<25; idet++) {
      
      refmean[iii-4][idet] = fhRawRec2d->	
	ProjectionY(Form("%s_py_%i_%i",                                                              fhRawRec2d ->GetName(), idet,iii-4),
		    idet,idet+1)->GetMean();
      
      refrms[iii-4][idet] = fhRawRec2d->
	ProjectionY(Form("%s_py%i_%i", 
			 fhRawRec2d ->GetName(), idet,iii-4),
		    idet,idet+1)->GetRMS();
      
      }
  }

  
  TString nameDev[2] = {"CDF", "LED"};
  for (Int_t icase=208; icase<210; icase++) {
    fhTime = (TH2*) listrec->At(icase);
    for (Int_t idet=1; idet<25; idet++) {
      Double_t binmean = fhTime->
      ProjectionY(Form("%s_py_%i_%i",                                                              fhTime ->GetName(), idet,icase),
		    idet,idet+1)->GetMean();
//        Double_t rms= fhTime ->ProjectionY(Form("%s_py%i_%i", 
//                                              fhTime ->GetName(), idet,icase),
//                                         idet,idet+1)->GetRMS();
      Double_t diffmean = binmean-refmean[icase-208][idet];
      
      if (TMath::Abs(diffmean) < 2 ) {
	checkr = 1;
	//	printf(" Laser calibration signal sits on its place %f for PMT %s %i : check = %f\n",  diffmean, nameDev[icase-208].Data() ,idet, checkr);
	AliDebug(AliQAv1::GetQADebugLevel(),
		 Form(" Laser calibration signal sits on its place %f for PMT %s %i : check = %f\n",  diffmean, nameDev[icase-208].Data(),idet, checkr)) ; 
     }
     if (TMath::Abs(diffmean) <= 5 && TMath::Abs(diffmean) >= 2 ){
       checkr = 0.5;
       // printf(" Laser calibration signal shifted by  %f ps for PMT %s %i : check = %f\n",  diffmean*24.4, nameDev[icase-208].Data(),idet, checkr);
       AliDebug(AliQAv1::GetQADebugLevel(),
		Form(" Laser calibration signal shifted by  %f ps (%f channels) for PMT %s %i : check = %f\n",  diffmean*24.4 ,diffmean , nameDev[icase-208].Data(),idet, checkr)) ; 
      }
     if (TMath::Abs(diffmean) > 5) { 
       checkr = 0.25;
       //   printf(" Big problems :laser calibration signal shifted by  %f ps (%f channels) for PMT %s %i : check = %f\n",  diffmean*24.4, diffmean, nameDev[icase-208].Data(),idet, checkr);
       AliDebug(AliQAv1::GetQADebugLevel(),
		Form(" Big problems :laser calibration signal shifted by  %f ps (%f channels) for PMT %s %i : check = %i\n",  diffmean*24.4, diffmean, nameDev[icase-208].Data(),idet, checkr)) ; 
       
     }
     
   }
 }

  return checkr;
}

