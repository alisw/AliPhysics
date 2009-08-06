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
Double_t * AliT0QAChecker::Check(AliQAv1::ALITASK_t index,TObjArray ** list)
{

  // Super-basic check on the QA histograms on the input list:
  // look whether they are empty!
  
  Double_t * test = new Double_t[AliRecoParam::kNSpecies] ; 
  /*
  char * detOCDBDir = Form("T0/%s/%s", AliQAv1::GetRefOCDBDirName(), AliQAv1::GetRefDataDirName()) ; 
  AliCDBEntry *QARefRec = AliQAManager::QAManager()->Get(detOCDBDir);
  //  QARefRec->Dump();
  if( !QARefRec){
    AliInfo("QA reference data NOT retrieved for Reconstruction check. No T0 reference distribution");
  }
  */
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
  TH1 *fhESD[2];
  
  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
    //  TString dataType = AliQAv1::GetAliTaskName(index);
    if (list[specie]->GetEntries() == 0){
      test[specie] = 1. ; // nothing to check
    }
    else {
    TIter next(list[specie]) ;
    TH1 * hdata ;
    TH2 * h ;
    //  AliDebug(AliQAv1::GetQADebugLevel(), Form(" data type %i %s nentries %i\n",
    //	   index,dataType.Data(),list->GetEntries()));
    for (Int_t ir=0; ir<list[specie]->GetEntries(); ir++) {
      //raw
      if(index==0 ){
        if(ir > 204 && ir<208 ) {
          //	  else{
          hdata = (TH1*) list[specie]->UncheckedAt(ir);
          AliDebug(AliQAv1::GetQADebugLevel(), Form(" index %i ir %i \n", index,ir));
          if(hdata) {
            cname = hdata->GetName();
            hname[ir] = cname;
            AliDebug(AliQAv1::GetQADebugLevel(),Form("count %i %s \n",ir, hname[ir].Data())) ;
            fhRawEff[ir] = hdata;
          }
        }
     }
      
      //rec
      if(index==2){
        h = (TH2*) list[specie]->UncheckedAt(ir);
        if(h) {
          cname = h->GetName();
          hname[ir] = cname;
	  //          AliDebug(AliQAv1::GetQADebugLevel(), Form("count %i %s \n",ir, hname[ir].Data())) ;
           fhRecDiff[ir] = h;
        }
      }
      //esd
      if(index==3){
        hdata = (TH1*) list[specie]->UncheckedAt(ir);
        if(hdata){
          fhESD[ir] = hdata;
          AliDebug(AliQAv1::GetQADebugLevel(), Form("count %i %s ",ir, hname[ir].Data()) );
        }
      }
    }
    //raw data
    if (index == 0) test[specie] = CheckRaw(list[specie]/*,(TObjArray *)QARefRec->GetObject()*/);


      
      if(index == 2){
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
      if (index == 3) {
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
  
  return test ;
}

//--------------------------------------------------------------------------
Double_t AliT0QAChecker::CheckRaw(TObjArray *listrec /*, TObjArray *listref*/) const
{
  
  TH1 *fhRawEff;
  // TH2 *fhRawRef;
  // TIter next(listref) ;
  // Int_t counter=0;
  // Float_t refmean[50]; 
  // Float_t refrms[50]; 
  Float_t checkr = 0;
  /*
  //  Int_t nref = listref->GetEntries(); 
  // Int_t nrec = listrec->GetEntries(); 
  
  cout<<" entries in ref "<<nref<<" in rec "<<nrec<<endl;
  //  for (Int_t iref=0; iref<nref; iref++){
  while  (fhRawRef = dynamic_cast<TH2 *>(next()))  {
  //  fhRawRef->Dump();
    // fhRawRef =(TH2*) listref->At(iref); 
    cout<<counter<<" hist "<<fhRawRef<<endl;
    fhRawRef->Print();
    fhRawRef->Dump();
    refmean[counter] = fhRawRef->GetMean(2);
    cout<<counter<<" mean "<<refmean[counter]<<endl;
    refrms[counter] = fhRawRef->GetRMS(2);
    cout<<counter<<" rms "<<refrms[counter]<<endl;
    counter++;
     cout<<" !!!!! reference "<<counter<<" "<<refmean[counter]<<" "<<refrms[counter]<<endl;
  }
  */
  
  for (Int_t icase=205; icase<208; icase++) {
    fhRawEff = (TH1*) listrec->At(icase);
    for (Int_t idet=0; idet<24; idet++) {
      Double_t mean = fhRawEff->GetBinContent(idet);
      AliDebug(AliQAv1::GetQADebugLevel(), 
	       Form("name %s icase %i idet %i mean %f \n",
		    fhRawEff->GetName(), icase, idet, mean));
       
      if (mean<1.2 && mean> 0.8 ) {
	checkr = 1;
	AliDebug(AliQAv1::GetQADebugLevel(), Form("All channels works meane efficieny %f with  test %f",  mean, checkr)) ; 
      }
      if (mean<=0.8 && mean>= 0.5 ){
	checkr = 0.5;
	AliDebug(AliQAv1::GetQADebugLevel(), 
		 Form("%s problem in channel %i  efficieny %f test %f",
		      fhRawEff->GetName(), idet,  mean,checkr )) ; 
      }
      if (mean<0.5 ) { 
	checkr = 0.25;
	AliDebug(AliQAv1::GetQADebugLevel(),
		 Form("%s big problem in channel %i  efficieny %f test %f",
		      fhRawEff->GetName(), idet,  mean,checkr )) ; 
     }
    }
  }
  /*  
  for (Int_t icase=208; icase<210; icase++) {
   fhRawEff = (TH2*) listrec->At(icase);
   
   for (Int_t idet=0; idet<24; idet++) {
     Double_t mean = fhRawEff->
       ProjectionY(Form("%s_py_%i_%i",                                                              fhRawEff ->GetName(), idet,icase),
		   idet,idet+1)->GetMean();
     Double_t rms= fhRawEff ->
       ProjectionY(Form("%s_py%i_%i", 
			fhRawEff ->GetName(), idet,icase),
		   idet,idet+1)->GetRMS();
     
    
     
     AliDebug(AliQAv1::GetQADebugLevel(), 
	      Form("name %s icase %i idet %i mean %f, rms %f\n",
		   fhRawEff->GetName(), icase, idet, mean,rms));
  
   }

  }
  */
  return checkr;
}

