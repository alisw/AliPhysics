/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

// *****************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  W.Ferrarese  P.Cerello  Mag 2008
//  INFN Torino

// --- ROOT system ---
#include "TH1.h"
#include <Riostream.h>

// --- AliRoot header files ---
#include "AliITSQAChecker.h"
#include "AliITSQASPDChecker.h"
#include "AliITSQASDDChecker.h"
#include "AliITSQASSDChecker.h"
#include "AliITSQADataMakerRec.h"

ClassImp(AliITSQAChecker)

//____________________________________________________________________________
AliITSQAChecker::AliITSQAChecker(Bool_t kMode, Short_t subDet, Short_t ldc) :
AliQACheckerBase("ITS","SDD Quality Assurance Checker"),
fkOnline(0),
fDet(0),  
fLDC(0),
fSPDOffset(0), 
fSDDOffset(0), 
fSSDOffset(0),
fSPDHisto(0),
fSDDHisto(0),
fSSDHisto(0),
fSPDChecker(0),  // SPD Checker
fSDDChecker(0),  // SDD Checker
fSSDChecker(0)  // SSD Checker

{
  // Standard constructor
  fkOnline = kMode; fDet = subDet; fLDC = ldc;
  if(fDet == 0 || fDet == 1) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQAChecker::Create SPD Checker\n");
    fSPDChecker = new AliITSQASPDChecker();
  }
  if(fDet == 0 || fDet == 2) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQAChecker::Create SDD Checker\n");
    fSDDChecker = new AliITSQASDDChecker();
  }
  if(fDet == 0 || fDet == 3) {
    AliDebug(AliQAv1::GetQADebugLevel(),"AliITSQAChecker::Create SSD Checker\n");
    fSSDChecker = new AliITSQASSDChecker();
  }
  InitQACheckerLimits();
}

//____________________________________________________________________________
AliITSQAChecker::AliITSQAChecker(const AliITSQAChecker& qac):
AliQACheckerBase(qac.GetName(), qac.GetTitle()), 
fkOnline(qac.fkOnline), 
fDet(qac.fDet), 
fLDC(qac.fLDC), 
fSPDOffset(qac.fSPDOffset), 
fSDDOffset(qac.fSDDOffset), 
fSSDOffset(qac.fSSDOffset), 
fSPDHisto(qac.fSPDHisto),
fSDDHisto(qac.fSDDHisto),
fSSDHisto(qac.fSSDHisto),
fSPDChecker(qac.fSPDChecker), 
fSDDChecker(qac.fSDDChecker), 
fSSDChecker(qac.fSSDChecker)
{
  // copy constructor
  AliError("Copy should not be used with this class\n");
}
//____________________________________________________________________________
AliITSQAChecker& AliITSQAChecker::operator=(const AliITSQAChecker& qac){
  // assignment operator
  this->~AliITSQAChecker();
  new(this)AliITSQAChecker(qac);
  return *this;
}


//____________________________________________________________________________
AliITSQAChecker::~AliITSQAChecker(){
  // destructor
  if(fSPDChecker)delete fSPDChecker;
  if(fSDDChecker)delete fSDDChecker;
  if(fSSDChecker)delete fSSDChecker;

}
//____________________________________________________________________________
void AliITSQAChecker::Check(Double_t * rv, AliQAv1::ALITASK_t index, TObjArray ** list, const AliDetectorRecoParam * recoParam)
{


  // basic checks on the QA histograms on the input list
  //for the ITS subdetectorQA (Raws Digits Hits RecPoints SDigits) return the worst value of the three result
  if(index == AliQAv1::kESD){
    
    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      rv[specie] = 0.0 ; 
      if ( !AliQAv1::Instance()->IsEventSpecieSet(specie) ) 
        continue ; 
      AliDebug(AliQAv1::GetQADebugLevel(),"Checker for ESD");
      Int_t tested = 0;
      Int_t empty = 0;
      // The following flags are set to kTRUE if the corresponding
      // QA histograms exceed a given quality threshold
      Bool_t cluMapSA = kFALSE;
      Bool_t cluMapMI = kFALSE;
      Bool_t cluMI = kFALSE;
      Bool_t cluSA = kFALSE;
      Bool_t verSPDZ = kFALSE;
      if (list[specie]->GetEntries() == 0) {
        rv[specie] = 0.; // nothing to check
      }
      else {
	Double_t *stepbit=new Double_t[AliQAv1::kNBIT];
	Double_t histonumb= list[specie]->GetEntries();
	CreateStepForBit(histonumb,stepbit); 
        TIter next1(list[specie]);
        TH1 * hdata;
        Int_t nskipped=0;
        Bool_t skipped[6]={kFALSE,kFALSE,kFALSE,kFALSE,kFALSE,kFALSE};
        // look for layers that we wanted to skip
        while ( (hdata = dynamic_cast<TH1 *>(next1())) ) {
          if(!hdata) continue;
          TString hname = hdata->GetName();
          if(!hname.Contains("hESDSkippedLayers")) continue;
          for(Int_t k=1; k<7; k++) {
            if(hdata->GetBinContent(k)>0) { 
              nskipped++; 
              skipped[k-1]=kTRUE; 
            } 
          } 
        }
        TIter next(list[specie]);
        while ( (hdata = dynamic_cast<TH1 *>(next())) ) {
          if(hdata){
            TString hname = hdata->GetName();
            Double_t entries = hdata->GetEntries();
            ++tested;
            if(!(entries>0.))++empty;
            AliDebug(AliQAv1::GetQADebugLevel(),Form("ESD hist name %s - entries %12.1g",hname.Data(),entries));
            if(hname.Contains("hESDClusterMapSA") && entries>0.){
              cluMapSA = kTRUE;
              AliDebug(AliQAv1::GetQADebugLevel(),Form("Processing histogram %s",hname.Data()));
              // Check if there are layers with anomalously low 
              // contributing points to SA reconstructed tracks
              for(Int_t k=1;k<7;k++){
                // check if the layer was skipped
                if(skipped[k-1]) continue;
                if(hdata->GetBinContent(k)<0.5*(entries/6.)){
                  cluMapSA = kFALSE;
                  AliDebug(AliQAv1::GetQADebugLevel(),Form("SA tracks have few points on layer %d - look at histogram hESDClustersSA",k));
                }
              }  
            }

            else if(hname.Contains("hESDClusterMapMI") && entries>0.){
              // Check if there are layers with anomalously low 
              // contributing points to MI reconstructed tracks
              AliDebug(AliQAv1::GetQADebugLevel(),Form("Processing histogram %s",hname.Data()));
              cluMapMI = kTRUE;
              for(Int_t k=1;k<7;k++){
                // check if the layer was skipped
                if(skipped[k-1]) continue;
                if(hdata->GetBinContent(k)<0.5*(entries/6.)){
                  cluMapMI = kFALSE;
                  AliDebug(AliQAv1::GetQADebugLevel(),Form("MI tracks have few points on layer %d - look at histogram hESDClustersMI",k));
                }
              }  
            }

            else if(hname.Contains("hESDClustersMI") && entries>0.){
              // Check if 6 clusters MI tracks are the majority
              AliDebug(AliQAv1::GetQADebugLevel(),Form("Processing histogram %s",hname.Data()));
              cluMI = kTRUE;
              Double_t maxlaytracks = hdata->GetBinContent(7-nskipped);
              for(Int_t k=2; k<7-nskipped; k++){
                if(hdata->GetBinContent(k)>maxlaytracks){
                  cluMI = kFALSE;
                  AliDebug(AliQAv1::GetQADebugLevel(),Form("MI Tracks with %d clusters are more than tracks with %d clusters. Look at histogram hESDClustersMI",k-1,6-nskipped));
                }
              }
            }

            else if(hname.Contains("hESDClustersSA") && entries>0.){
              // Check if 6 clusters SA tracks are the majority
              AliDebug(AliQAv1::GetQADebugLevel(),Form("Processing histogram %s",hname.Data()));
              cluSA = kTRUE;
              Double_t maxlaytracks = hdata->GetBinContent(7-nskipped);
              for(Int_t k=2; k<7-nskipped; k++){
                if(hdata->GetBinContent(k)>maxlaytracks){
                  cluSA = kFALSE;
                  AliDebug(AliQAv1::GetQADebugLevel(), Form("SA Tracks with %d clusters are more than tracks with %d clusters. Look at histogram hESDClustersSA",k-1,6-nskipped));
                }
              }
            }

            else if(hname.Contains("hSPDVertexZ") && entries>0.){
              // Check if average Z vertex coordinate is -5 < z < 5 cm
              AliDebug(AliQAv1::GetQADebugLevel(),Form("Processing histogram %s",hname.Data()));
              verSPDZ = kTRUE;
              if(hdata->GetMean()<-5. && hdata->GetMean()>5.){
                verSPDZ = kFALSE;
                AliDebug(AliQAv1::GetQADebugLevel(),Form("Average z vertex coordinate is at z= %10.4g cm",hdata->GetMean()));
              }
            }
          }
          else{
            AliError("ESD Checker - invalid data type");
          }
	}
	rv[specie] = 0.;
	if(tested>0){
	  if(tested == empty){
	    rv[specie] = 2500.; // set to error
	    AliWarning(Form("All ESD histograms are empty - specie=%d",specie));
	  }
	  else {
	    rv[specie] = 2500.-1500.*(static_cast<Double_t>(tested-empty)/static_cast<Double_t>(tested)); // INFO if all histos are filled
	    if(cluMapSA)rv[specie]-=200.;
	    if(cluMapMI)rv[specie]-=200.;
	    if(cluMI)rv[specie]-=200.;
	    if(cluSA)rv[specie]-=200.;
	    if(verSPDZ)rv[specie]-=199.;  // down to 1 if everything is OK
	  }
	}
      }
      //     AliDebug(AliQAv1::GetQADebugLevel(), Form("ESD - Tested %d histograms, Return value %f \n",tested,rv[specie]));
      AliInfo(Form("ESD - Tested %d histograms, Return value %f \n",tested,rv[specie]));
    }
  }  // end of ESD QA
  else{
    
    //____________________________________________________________________________

    Double_t spdCheck[AliRecoParam::kNSpecies] ;
    Double_t sddCheck[AliRecoParam::kNSpecies] ;
    Double_t ssdCheck[AliRecoParam::kNSpecies] ;



    for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {
      if ( !AliQAv1::Instance()->IsEventSpecieSet(specie)) continue; 
      if ( AliQAv1::Instance()->IsEventSpecieSet(specie) ) {
	Double_t histotot=list[specie]->GetEntries();
	if(histotot!=0)
	  {
	    spdCheck[specie]=0.;
	    sddCheck[specie]=0.;
	    ssdCheck[specie]=0.;
	    rv[specie] = 0.0 ;// 
	    //pixel
	    if(fDet == 0 || fDet == 1) {
	      fSPDChecker->SetTaskOffset(fSPDOffset);
	      //printf("spdoffset = %i \n",fSPDOffset );
	      Double_t histoSPD=double(GetSPDHisto());
	      if(AliITSQADataMakerRec::AreEqual(histoSPD,0)==kFALSE){
		Double_t *stepSPD=new Double_t[AliQAv1::kNBIT];
		CreateStepForBit(histoSPD,stepSPD);
		fSPDChecker->SetStepBit(stepSPD);
		spdCheck[specie] = fSPDChecker->Check(index, list[specie], recoParam);
		if(spdCheck[specie]>fUpTestValue[AliQAv1::kFATAL]||spdCheck[specie]<0.)
		  {
		    AliInfo(Form("SPD check result for %s  is out of range (%f)!!! Retval of specie %s is sit to -1\n ",AliQAv1::GetAliTaskName(index),spdCheck[specie],AliRecoParam::GetEventSpecieName(specie)));
		    spdCheck[specie]=fUpTestValue[AliQAv1::kFATAL];
		  }
		delete []stepSPD;
	      }//end check SPD entries
	      else{spdCheck[specie]=fUpTestValue[AliQAv1::kFATAL];}
	      rv[specie]=spdCheck[specie];
	    }//end SPD check
	    //drift
	    if(fDet == 0 || fDet == 2) {
	      fSDDChecker->SetTaskOffset(fSDDOffset);
	      Double_t histoSDD=double(GetSDDHisto());
	      if(AliITSQADataMakerRec::AreEqual(histoSDD,0)==kFALSE){
		Double_t *stepSDD=new Double_t[AliQAv1::kNBIT];
		CreateStepForBit(histoSDD,stepSDD);
		fSDDChecker->SetStepBit(stepSDD);
		sddCheck[specie] = fSDDChecker->Check(index, list[specie], recoParam);	
		if(sddCheck[specie]>fUpTestValue[AliQAv1::kFATAL]||sddCheck[specie]<0.)
		  {
		    AliInfo(Form("SDD check result for %s  is out of range (%f)!!! Retval of specie %s is sit to -1\n ",AliQAv1::GetAliTaskName(index),sddCheck[specie],AliRecoParam::GetEventSpecieName(specie)));
		    sddCheck[specie]=fUpTestValue[AliQAv1::kFATAL];
		  }
		delete []stepSDD;
	      }//end check SDD entries
	      else{ssdCheck[specie]=fUpTestValue[AliQAv1::kFATAL];}
	      if(sddCheck[specie]>rv[specie])rv[specie]=sddCheck[specie];  
	    }//end SDD
	    //strip
	    if(fDet == 0 || fDet == 3) {
	      fSSDChecker->SetTaskOffset(fSSDOffset);
	      Double_t histoSSD=double(GetSSDHisto());
	      if(AliITSQADataMakerRec::AreEqual(histoSSD,0)==kFALSE){
	      Double_t *stepSSD=new Double_t[AliQAv1::kNBIT];
	      CreateStepForBit(histoSSD,stepSSD);
	      fSSDChecker->SetStepBit(stepSSD);
	      ssdCheck[specie] = fSSDChecker->Check(index, list[specie], recoParam);
	      if(ssdCheck[specie]>fUpTestValue[AliQAv1::kFATAL]||ssdCheck[specie]<0.)
		{
		  AliInfo(Form("SSD check result for %s is out of range (%f)!!! Retval of specie %s is sit to -1\n ",AliQAv1::GetAliTaskName(index),ssdCheck[specie],AliRecoParam::GetEventSpecieName(specie)));
		  ssdCheck[specie]=fUpTestValue[AliQAv1::kFATAL];
		}
	      delete [] stepSSD;
	      }//end check SSD entries
	      else{ssdCheck[specie]=fUpTestValue[AliQAv1::kFATAL];}
	      if(ssdCheck[specie]>rv[specie])rv[specie]=ssdCheck[specie];
	    }//end SSD
	    
	    AliInfo(Form("Check result for %s: \n\t  SPD %f \n\t  SDD %f \n\t  SSD %f \n Check result %f \n ",AliQAv1::GetAliTaskName(index),spdCheck[specie],sddCheck[specie],ssdCheck[specie],rv[specie]));
	    // here merging part for common ITS QA result
	    // 
	  }//end entries
      }//end if event specie
    }//end for
  }
}

//____________________________________________________________________________
void AliITSQAChecker::SetTaskOffset(Int_t SPDOffset, Int_t SDDOffset, Int_t SSDOffset)
{
  //Setting the 3 offsets for each task called
  fSPDOffset = SPDOffset;
  fSDDOffset = SDDOffset;
  fSSDOffset = SSDOffset;
}

//____________________________________________________________________________
void AliITSQAChecker::SetHisto(Int_t SPDhisto, Int_t SDDhisto, Int_t SSDhisto)
{
  //Setting the 3 offsets for each task called
  fSPDHisto = SPDhisto;
  fSDDHisto = SDDhisto;
  fSSDHisto = SSDhisto;
}

 //____________________________________________________________________________
 void AliITSQAChecker::SetDetTaskOffset(Int_t subdet,Int_t offset)
 {
   switch(subdet){
   case 1:
     SetSPDTaskOffset(offset);
     break;
   case 2:
     SetSDDTaskOffset(offset);
     break;
   case 3:
     SetSSDTaskOffset(offset);
     break;
   default:
     AliWarning("No specific (SPD,SDD or SSD) subdetector correspond to to this number!!! all offsets set to zero for all the detectors\n");
     SetTaskOffset(0, 0, 0);
     break;
   }
 }

 //____________________________________________________________________________
 void AliITSQAChecker::SetDetHisto(Int_t subdet,Int_t histo)
 {
   switch(subdet){
   case 1:
     SetSPDHisto(histo);
     break;
   case 2:
     SetSDDHisto(histo);
     break;
   case 3:
     SetSSDHisto(histo);
     break;
   default:
     AliWarning("No specific (SPD,SDD or SSD) subdetector correspond to to this number!!! all offsets set to zero for all the detectors\n");
     SetHisto(0, 0, 0);
     break;
   }
 }

//_____________________________________________________________________________

void AliITSQAChecker::InitQACheckerLimits()
{
  
  AliInfo("Setting of tolerance values\n");

  Float_t lowtolerancevalue[AliQAv1::kNBIT];

  Float_t hightolerancevalue[AliQAv1::kNBIT];
  for(Int_t bit=0;bit<AliQAv1::kNBIT;bit++)
    {
      lowtolerancevalue[bit]=(bit*1000.);
      hightolerancevalue[bit]=((bit+1.)*1000.);
    }
  SetHiLo(hightolerancevalue,lowtolerancevalue);
  //  AliInfo(Form("Range Value  \n INFO    -> %f <  value <  %f \n WARNING -> %f <  value <= %f \n ERROR   -> %f <  value <= %f \n FATAL   -> %f <= value <  %f \n", fLowTestValue[AliQAv1::kINFO], fUpTestValue[AliQAv1::kINFO], fLowTestValue[AliQAv1::kWARNING], fUpTestValue[AliQAv1::kWARNING], fLowTestValue[AliQAv1::kERROR], fUpTestValue[AliQAv1::kERROR], fLowTestValue[AliQAv1::kFATAL], fUpTestValue[AliQAv1::kFATAL]  ));

  if(fDet == 0 || fDet == 1) {
    fSPDChecker->SetSPDLimits( lowtolerancevalue,hightolerancevalue );
  }
  if(fDet == 0 || fDet == 2) {
    fSDDChecker->SetSDDLimits( lowtolerancevalue,hightolerancevalue );
  }
  if(fDet == 0 || fDet == 3) {
    fSSDChecker->SetSSDLimits( lowtolerancevalue,hightolerancevalue );
  }


  
}


//_____________________________________________________________________________

void AliITSQAChecker::CreateStepForBit(Double_t histonumb,Double_t *steprange)
{
  for(Int_t bit=0;bit < AliQAv1::kNBIT; bit++)
    {       
      //printf("%i\t %f \t %f \t %f \n",bit, fUpTestValue[bit],fLowTestValue[AliQAv1::kINFO],histonumb);
      steprange[bit]=double((fUpTestValue[bit] - fLowTestValue[AliQAv1::kINFO])/histonumb);
      //printf("%i\t %f \t %f \t %f \t %f\n",bit, fUpTestValue[bit],fLowTestValue[AliQAv1::kINFO],histonumb,steprange[bit] );
    }
  //AliInfo(Form("StepBitValue:numner of histo %f\n\t INFO %f \t WARNING %f \t ERROR %f \t FATAL %f \n",histonumb, steprange[AliQAv1::kINFO],steprange[AliQAv1::kWARNING],steprange[AliQAv1::kERROR],steprange[AliQAv1::kFATAL]));
}


//_____________________________________________________________________________
void AliITSQAChecker::SetQA(AliQAv1::ALITASK_t index, Double_t * value) const
{

  AliQAv1 * qa = AliQAv1::Instance(index) ;


  for (Int_t specie = 0 ; specie < AliRecoParam::kNSpecies ; specie++) {

    if (! qa->IsEventSpecieSet(AliRecoParam::ConvertIndex(specie)))
      continue ;
    if (  value == NULL ) { // No checker is implemented, set all QA to Fatal
      qa->Set(AliQAv1::kFATAL, specie) ; 
    } else {
      if ( value[specie] > fLowTestValue[AliQAv1::kFATAL] && value[specie] <= fUpTestValue[AliQAv1::kFATAL] ) 
        qa->Set(AliQAv1::kFATAL, AliRecoParam::ConvertIndex(specie)) ; 
      else if ( value[specie] > fLowTestValue[AliQAv1::kERROR] && value[specie] <= fUpTestValue[AliQAv1::kERROR]  )
        qa->Set(AliQAv1::kERROR, AliRecoParam::ConvertIndex(specie)) ; 
      else if ( value[specie] > fLowTestValue[AliQAv1::kWARNING] && value[specie] <= fUpTestValue[AliQAv1::kWARNING]  )
        qa->Set(AliQAv1::kWARNING, AliRecoParam::ConvertIndex(specie)) ;
      else if ( value[specie] > fLowTestValue[AliQAv1::kINFO] && value[specie] <= fUpTestValue[AliQAv1::kINFO] ) 
        qa->Set(AliQAv1::kINFO, AliRecoParam::ConvertIndex(specie)) ; 	
      //else if(value[specie]==0) qa->Set(AliQAv1::kFATAL, AliRecoParam::ConvertIndex(specie)) ; //no ckeck has been done
    }
    qa->ShowStatus(AliQAv1::kITS,index,AliRecoParam::ConvertIndex(specie));
  }//end for

}


