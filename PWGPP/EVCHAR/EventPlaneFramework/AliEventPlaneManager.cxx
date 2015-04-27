/*
***********************************************************
  Manager for event plane corrections framework
  Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
  2014/12/10
  *********************************************************
*/

#ifndef ALIEVENTPLANEMANAGER_H
#include "AliEventPlaneManager.h"
#endif

#include "AliEventPlaneHistos.h"
#include "AliEventPlaneConfiguration.h"
#include "AliEventPlaneDetector.h"
#include "AliEventPlaneVarManager.h"
#include <TMath.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TArrayS.h>
#include <iostream>

ClassImp(AliEventPlaneManager)

TClonesArray* AliEventPlaneManager::fgEventPlaneConfigurations[AliEventPlaneManager::kNdetectors] = {0};
TClonesArray* AliEventPlaneManager::fgDetector[AliEventPlaneManager::kNdetectors] = {0};
TList* AliEventPlaneManager::fgQvectors = 0;


//_______________________________________________________________________________
AliEventPlaneManager::AliEventPlaneManager() :
 fDetector(),
 fNdetectors(),
 fNEventPlaneDetectors(),
 //fNMasterDetectors(0),
 fQvectors(0x0),
 fRunLightWeight(kFALSE),
 fCorrelationHists(),
 fEqualizationHists(),
 fVZEROminMult(0.5),
 fTZEROminMult(1.e-2),
 fZDCminMult(0.01),
 fFMDminMult(1.e-6),
 fInputFriendFileName(""),
 fOutputFriendFileName("")
{
  //
  // default constructor
  //

for(Int_t i=0; i<kNdetectors; ++i) fNEventPlaneDetectors[i]=0;
   


  if(!fgDetector[0])
   for(Int_t i=0; i<kNdetectors; ++i) fgDetector[i] = new TClonesArray("AliEventPlaneDetector", 100000);
  for(Int_t i=0; i<kNdetectors; ++i) fDetector[i] = fgDetector[i];
  if(!fgEventPlaneConfigurations[0])
   for(Int_t i=0; i<kNdetectors; ++i) fgEventPlaneConfigurations[i] = new TClonesArray("AliEventPlaneConfiguration", fgkEPMaxDetectors);
  for(Int_t i=0; i<kNdetectors; ++i) fEventPlaneConfigurations[i] = fgEventPlaneConfigurations[i];

  if(!fgQvectors){
    fgQvectors = new TList();
    for(Int_t idet=0; idet<fgkEPMaxDetectors; ++idet){ 
      TClonesArray * arrayQvectors = new TClonesArray("AliEventPlaneQvector", 100000);
      arrayQvectors->SetName(Form("AliEventPlaneQvector%d", idet));
      fgQvectors->Add(arrayQvectors);
    }
  }
  fQvectors = fgQvectors;

  for(Int_t idet=0; idet<10000; idet++) {fCorrelationHists[idet][0]=-1; fCorrelationHists[idet][1]=-1;}
  for(Int_t idet=0; idet<1000; idet++) fEqualizationHists[idet]=-1;
  
  //if(!fgEventPlaneConfigurations) fgEventPlaneConfigurations = new TClonesArray("AliEventPlaneConfiguration", fgkEPMaxDetectors);
  //fEventPlaneConfigurations = fgEventPlaneConfigurations;

//  if(!fgEventPlaneMasters) fgEventPlaneMasters = new TClonesArray("AliEventPlaneMasterDetector", fgkEPMaxDetectors);
//  fEventPlaneMasters = fgEventPlaneMasters;
}


//____________________________________________________________________________
AliEventPlaneManager::~AliEventPlaneManager()
{
  //
  // De-Constructor
  //
  ClearEvent();
}


//_____________________________________________________________________________
void AliEventPlaneManager::ClearEvent() {
  //
  // clear the event
  //
  for(Int_t i=0; i<kNdetectors; ++i) if(fDetector[i]) fDetector[i]->Clear("C");
  for(Int_t idet=0; idet<NEventPlaneDetectors(); ++idet) {
    if(fQvectors) ((TClonesArray*)fQvectors->At(idet))->Clear("C");
  }
}


//_____________________________________________________________________________
void AliEventPlaneManager::SetEmptyQvectors() {
  //
  // clear the event
  //
  for(Int_t idet=0; idet<NEventPlaneDetectors(); ++idet) {
      TClonesArray& qvec = *(GetQvectors(idet));
      AliEventPlaneQvector *reducedQvector=new(qvec[0]) AliEventPlaneQvector();
  }
}


////____________________________________________________________________________
//void AliEventPlaneManager::CopyEvent(const AliEventPlaneManager* event) {
//  //
//  // copy information from another event to this one
//  //
//  for(Int_t i=0; i<kNdetectors; ++i) {if(event->GetReducedDetector(i)) fDetector[i]=event->GetReducedDetector(i);}
//  fQvectors = event->GetQvectors();
//
//
//
//}


//_________________________________________________________________
void AliEventPlaneManager::GetQvector(Bool_t useFriendChain, Bool_t useEqualizedWeights, Float_t * values){
  //
  // Construct the event plane for a specified detector
  //

  const Int_t nbins=1;


  Int_t eqmethod, epdet;
  for(Int_t idet=0; idet<kNdetectors; idet++){
    if(!GetReducedDetector(idet)) continue;

    if(useEqualizedWeights&&idet==kTPC) continue;  // TPC has no equalizing weights option

    Double_t Qvec[50][nbins][6][2];
    Double_t fillValues[20];
    Double_t weight=0.0, absWeight = 0.0, x=0.0, y=0.0, phi; 
    Int_t * var;
    Int_t dim;
    Double_t mult[50][nbins];
    for(Int_t jdet=0; jdet<fNEventPlaneDetectors[idet]; jdet++) for(Int_t ibin=0; ibin<nbins; ibin++) mult[jdet][ibin] = 0.0;
    Int_t bin=0;
   
    TIter nextEntry(GetReducedDetector(idet));
    
    AliEventPlaneQvector* Qvector = 0x0;
    AliEventPlaneDetector* reducedDetector=0x0;
    AliEventPlaneConfiguration* EPconf = 0x0;
    //TIter nextEPconf(fEventPlaneConfigurations[idet]);


    ///if(useEqualizedWeights){   // Figure out if we need to run equalized qvector calculation
    ///  Bool_t calibrate = kFALSE;
    ///  for(Int_t iconf=0; iconf<fEventPlaneConfigurations[idet]->GetEntriesFast(); iconf++){
    ///    EPconf=EventPlaneConfiguration(idet, iconf);
    ///    //if(!EPconf) break;
    ///    if(EPconf->CalibrationStep()>1&&useFriendChain) continue;
    ///    if(EPconf->CalibrationMethod()!=-1&&EPconf->CalibrationStep()!=0) calibrate=kTRUE;
    ///  }
    ///  if(!calibrate) continue;
    ///}

    ///if(useFriendChain) if(!useEqualizedWeights){  // Figure out if we need to run raw qvector calculation
    ///  Bool_t calculateRaw = kFALSE;
    ///  for(Int_t iconf=0; iconf<1000; iconf++){
    ///    EPconf=EventPlaneConfiguration(idet, iconf);
    ///    if(!EPconf) break;
    ///    if(EPconf->CalibrationStep()>0) continue;
    ///    calculateRaw = kTRUE;}
    ///    if(!calculateRaw) continue;
    ///  }



    Int_t itrack=-1;
    while((reducedDetector=static_cast<AliEventPlaneDetector*>(nextEntry()))) {
      if(!reducedDetector) continue;

      itrack++;
      bin = reducedDetector->Bin();
      weight = reducedDetector->Weight();
      x = reducedDetector->X();
      y = reducedDetector->Y();

      //cout<<itrack<<"  "<<weight<<"  "<<x<<"  "<<y<<endl;
      if(bin<0) continue;

      for(Int_t iconf=0; iconf<fEventPlaneConfigurations[idet]->GetEntriesFast(); iconf++){
        EPconf=EventPlaneConfiguration(idet, iconf);
        //cout<<EPconf->EventPlaneDetectorName()<<"  "<<weight<<"  "<<useEqualizedWeights<<endl;
        //if(useEqualizedWeights) if(!(EPconf->doChannelEqualization()&&EPconf->CalibrationStep()>0)) continue;
        //if(!EPconf) break;


        epdet = EPconf->LocalIndex();
        if(!reducedDetector->CheckEventPlaneDetector(epdet)) continue;
        if(useEqualizedWeights){ 
          if(EPconf->CalibrationStep()<1) continue;
          if(!EPconf->doChannelEqualization()) continue;
        }

        Short_t equalizationMethod=EPconf->EqualizationMethod();
        
        var = EPconf->CalibrationBinning()->Var();
        dim = EPconf->CalibrationBinning()->Dim();

        for(Int_t iv=0; iv<dim; iv++) {fillValues[iv] = values[var[iv]];}
        
        
        if(mult[epdet][bin]==0.0) for(Int_t ih=0; ih<6; ih++) for(Int_t ic=0; ic<2; ic++) Qvec[epdet][bin][ih][ic] = 0.0;  // only reset to 0 if bin is used

        if(useEqualizedWeights) weight = reducedDetector->Weight(epdet,equalizationMethod);

        mult[epdet][bin]+=weight;


        //if(EPconf->EventPlaneDetectorName()=="VZEROA") cout<<EPconf->EventPlaneDetectorName()<<"  "<<weight<<"  "<<useEqualizedWeights<<"   "<<itrack<<endl;

        //  1st harmonic  
        Qvec[epdet][bin][0][0] += weight*x;
        Qvec[epdet][bin][0][1] += weight*y;
        //  2nd harmonic
        Qvec[epdet][bin][1][0] += weight*(2.0*x*x-1);
        Qvec[epdet][bin][1][1] += weight*(2.0*x*y);
        //  3rd harmonic
        Qvec[epdet][bin][2][0] += weight*(4.0*x*x*x-3.0*x);
        Qvec[epdet][bin][2][1] += weight*(3.0*y-4.0*y*y*y);
        //  4th harmonic
        Qvec[epdet][bin][3][0] += weight*(1.0-8.0*x*x*y*y);
        Qvec[epdet][bin][3][1] += weight*(4.0*x*y-8.0*x*y*y*y);
        //  5th harmonic
        Qvec[epdet][bin][4][0] += weight*(16.0*x*x*x*x*x-20.0*x*x*x+5.0*x);
        Qvec[epdet][bin][4][1] += weight*(16.0*y*y*y*y*y-20.0*y*y*y+5.0*y);
        //  6th harmonic
        Qvec[epdet][bin][5][0] += weight*(32.0*x*x*x*x*x*x-48.0*x*x*x*x+18.0*x*x-1.0);
        Qvec[epdet][bin][5][1] += weight*(x*y*(32.0*y*y*y*y-32.0*y*y+6.0)); 
        

        if(EPconf->TwistAndScalingMethod()==0){
          for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) {
            EPconf->U2nHistogram(ih,0)->Fill(fillValues,TMath::Cos(reducedDetector->Phi()*ih*2));
            EPconf->U2nHistogram(ih,1)->Fill(fillValues,TMath::Sin(reducedDetector->Phi()*ih*2));
          }
          EPconf->U2nHistogramE()->Fill(fillValues);
        }
      }
    }
           

  
  
    

 
    Int_t calmethod;
    Double_t div;
    for(Int_t iconf=0; iconf<fEventPlaneConfigurations[idet]->GetEntriesFast(); iconf++){
      EPconf=EventPlaneConfiguration(idet, iconf);
      if(useEqualizedWeights){
        if(EPconf->CalibrationStep()<1) continue;
        if(!EPconf->doChannelEqualization()) continue;
      }
      for(Int_t iv=0; iv<EPconf->CalibrationHistogramQ(0,fgkEPMinHarmonics,0)->GetNdimensions(); iv++) {fillValues[iv] = values[EPconf->CalibrationBinning()->Var(iv)];}
      
      Bool_t once=kTRUE;
      Short_t equalizationMethod=EPconf->EqualizationMethod();
      Int_t filledBin=-1;
      calmethod = EPconf->CalibrationMethod();
      epdet = EPconf->LocalIndex();
      //std::cout<<EPconf->EventPlaneDetectorName()<<"  "<<mult[epdet][0]<<"  "<<epdet<<std::endl;
      
      for(Int_t iBin=0; iBin<nbins; iBin++){
        if(mult[epdet][iBin]<1e-6){
          filledBin++;
          TClonesArray& qvec = *EPconf->Qvectors();
          Qvector=new(qvec[filledBin]) AliEventPlaneQvector();
          Qvector->SetBin(iBin);
          for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) Qvector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kUndefined);
          continue;
        }
        filledBin++;
        if(useEqualizedWeights) {
          Qvector = static_cast<AliEventPlaneQvector*>(EPconf->Qvectors()->At(iBin));
        }
        else{
          TClonesArray& qvect = *EPconf->Qvectors();
          Qvector=new(qvect[filledBin]) AliEventPlaneQvector();
          Qvector->SetBin(iBin);
        }
        Qvector->SetMultiplicity(mult[epdet][iBin]);
        for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) {
          if(Qvector->CheckEventPlaneStatus(ih,AliEventPlaneQvector::kUndefined)) continue;
          if(calmethod==0) div=TMath::Sqrt(mult[epdet][iBin]);
          else if(calmethod==1) div=mult[epdet][iBin];
          else if(calmethod==2) div=TMath::Sqrt(Qvec[epdet][iBin][ih-1][0]*Qvec[epdet][iBin][ih-1][0]+Qvec[epdet][iBin][ih-1][1]*Qvec[epdet][iBin][ih-1][1]);
          else if(calmethod==3) div=1.;
          Qvector->SetQx(ih, Qvec[epdet][iBin][ih-1][0]/div);
          Qvector->SetQy(ih, Qvec[epdet][iBin][ih-1][1]/div);
          Qvector->SetQx(ih, Qvec[epdet][iBin][ih-1][0]/div);
          Qvector->SetQy(ih, Qvec[epdet][iBin][ih-1][1]/div);
          if(useEqualizedWeights) {
            Qvector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kEqualized);
            EPconf->CalibrationHistogramQ(1,ih,0)->Fill(fillValues,Qvector->Qx(ih));
            EPconf->CalibrationHistogramQ(1,ih,1)->Fill(fillValues,Qvector->Qy(ih));
            if(once) EPconf->CalibrationHistogramE(AliEventPlaneQvector::kEqualized)->Fill(fillValues); once=kFALSE;
          }
          else { 
            Qvector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kRaw);
            EPconf->CalibrationHistogramQ(AliEventPlaneQvector::kRaw,ih,0)->Fill(fillValues,Qvector->Qx(ih));
            EPconf->CalibrationHistogramQ(AliEventPlaneQvector::kRaw,ih,1)->Fill(fillValues,Qvector->Qy(ih));
            if(once) EPconf->CalibrationHistogramE(AliEventPlaneQvector::kRaw)->Fill(fillValues); once=kFALSE;
            //std::cout<<EPconf->EventPlaneDetectorName()<<"   "<<Qvector<<"  "<<Qvector->Qx(ih)<<endl;
          }
        }

          
        } 
      }
    }

  return;

}

 
 //_____________________________________________________________________
void AliEventPlaneManager::CalibrateChannels(Float_t* values, Bool_t useFriendChain) {
   //
   // Calibrate channel multiplicities
   //
     

   TString det;
   Int_t* var;
   Int_t dim; 
   Double_t fillValues[20];

  for(Int_t idet=0; idet<kNdetectors; idet++){
   AliEventPlaneConfiguration* EPconf = 0x0;
   TIter nextEPconf(fEventPlaneConfigurations[idet]);
   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
    if(!EPconf) continue;
    if(!EPconf->doChannelEqualization()) continue;

     

     var = EPconf->EqualizationBinning()->Var();
     dim = EPconf->EqualizationBinning()->Dim();



     TIter nextEntry(GetReducedDetector(EPconf->DetectorType()));
     AliEventPlaneDetector* reducedDetector= 0x0;


     while((reducedDetector=static_cast<AliEventPlaneDetector*>(nextEntry()))) {
       if(!reducedDetector) continue;
       
       


       // Gain calibration of channels to average multiplicity per channel
       Double_t average, axentries, axentriesEvent, summ2, summ, n, mult, chanMult;
       Double_t width = 1.0;

       for(Int_t iav=0; iav<dim; iav++) fillValues[iav] = values[var[iav]];
       fillValues[dim-1] = reducedDetector->Id();    

       //cout<<fillValues[dim]<<"  "<<dim<<endl;
       mult = reducedDetector->Weight();
       if(mult>1E-2){
         EPconf->EqualizationHistogramM(0)->Fill(fillValues,mult);
         EPconf->EqualizationHistogramE(0)->Fill(fillValues);
       }
       if(EPconf->CalibrationStep()==0) continue;


       summ2 = EPconf->InputEqualizationHistogramM()->GetBinError(EPconf->InputEqualizationHistogramM()->GetBin(fillValues));
       summ2 = summ2*summ2;
       summ  = EPconf->InputEqualizationHistogramM()->GetBinContent(EPconf->InputEqualizationHistogramM()->GetBin(fillValues));
       n     = EPconf->InputEqualizationHistogramE()->GetBinContent(EPconf->InputEqualizationHistogramE()->GetBin(fillValues));

    //cout<<EPconf->EventPlaneDetectorName()<<endl;
       //cout<<summ<<"  "<<n<<endl;
       //cout<<n<<endl;

       if(n<=1){
         reducedDetector->SetAverageEqualizedWeight(EPconf->LocalIndex(), mult);
         reducedDetector->SetWidthEqualizedWeight(EPconf->LocalIndex(), mult);
         continue;
       }

       if((summ2/n-summ/n*summ/n)<=0.0) width = TMath::Sqrt(summ/n);
       else if((summ2/n-summ/n*summ/n)>0.0) width  = TMath::Sqrt(summ2/n-summ/n*summ/n);
       else width=1.0;

       average = summ/n; 

       chanMult=0.0;
       if(average > 1.0e-6){
         chanMult = mult/average;
       }
       reducedDetector->SetAverageEqualizedWeight(EPconf->LocalIndex(), chanMult);
       EPconf->EqualizationHistogramM(1)->Fill(fillValues,chanMult);
       EPconf->EqualizationHistogramE(1)->Fill(fillValues);

       if(average > 1.0e-6){
         chanMult = 1.+(mult-average)/width*0.1;
       }
       //cout<<EPconf->EventPlaneDetectorName()<<"  "<<mult<<"  "<<chanMult<<endl;
       reducedDetector->SetWidthEqualizedWeight(EPconf->LocalIndex(), chanMult);
       EPconf->EqualizationHistogramM(2)->Fill(fillValues,chanMult);
       EPconf->EqualizationHistogramE(2)->Fill(fillValues);

     }


    }
  }


   return;
 }


//______________________________________
void AliEventPlaneManager::FillCorrelationHistograms(Float_t value, Int_t step){

  TString detA,detB,detC;

  AliEventPlaneQvector* Qvectors[3];
  AliEventPlaneConfiguration* EPconf=0x0,*EPconfB=0x0,*EPconfC=0x0;
  for(Int_t idet=0; idet<AliEventPlaneManager::kNdetectors; idet++){ 
    TClonesArray* epConfList=GetEventPlaneConfigurations(idet);
    for(Int_t iconf=0; iconf<epConfList->GetEntriesFast(); iconf++){
      EPconf = (AliEventPlaneConfiguration*) epConfList->At(iconf);
      if(!EPconf) continue;
      if(EPconf->CalibrationStep()<step) continue;
      detB = EPconf->CorrelationDetectorName(0);
      detC = EPconf->CorrelationDetectorName(1);
      EPconfB = EventPlaneConfiguration(detB);
      EPconfC = EventPlaneConfiguration(detC);

      Qvectors[0] = static_cast<AliEventPlaneQvector*>(EPconf->Qvectors()->At(0));
      Qvectors[1] = static_cast<AliEventPlaneQvector*>(EPconfB->Qvectors()->At(0));
      Qvectors[2] = static_cast<AliEventPlaneQvector*>(EPconfC->Qvectors()->At(0));
          
      for(Int_t icomb=0; icomb<3; ++icomb){ 
        for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){ 
          if(Qvectors[icomb]->CheckEventPlaneStatus(ih,AliEventPlaneQvector::kUndefined)||Qvectors[(icomb+1)%3]->CheckEventPlaneStatus(ih,AliEventPlaneQvector::kUndefined)) continue;
          //cout<<step<<"  "<<EPconf->EventPlaneDetectorName()<<"   "<<EPconfB->EventPlaneDetectorName()<<"   "<<EPconfC->EventPlaneDetectorName()<<"  "<<Qvectors[icomb]->Qx(ih)<<"  "<<Qvectors[(icomb+1)%3]->Qx(ih)<<endl;
          EPconf->CorrelationProf(step,icomb,ih,0)->Fill(value,Qvectors[icomb]->Qx(ih)*Qvectors[(icomb+1)%3]->Qx(ih));
          EPconf->CorrelationProf(step,icomb,ih,1)->Fill(value,Qvectors[icomb]->Qy(ih)*Qvectors[(icomb+1)%3]->Qx(ih));
          EPconf->CorrelationProf(step,icomb,ih,2)->Fill(value,Qvectors[icomb]->Qx(ih)*Qvectors[(icomb+1)%3]->Qy(ih));
          EPconf->CorrelationProf(step,icomb,ih,3)->Fill(value,Qvectors[icomb]->Qy(ih)*Qvectors[(icomb+1)%3]->Qy(ih));
          EPconf->CorrelationEpProf(step,icomb,ih,0)->Fill(value,Qvectors[icomb]->QxNorm(ih)*Qvectors[(icomb+1)%3]->QxNorm(ih));
          EPconf->CorrelationEpProf(step,icomb,ih,1)->Fill(value,Qvectors[icomb]->QyNorm(ih)*Qvectors[(icomb+1)%3]->QxNorm(ih));
          EPconf->CorrelationEpProf(step,icomb,ih,2)->Fill(value,Qvectors[icomb]->QxNorm(ih)*Qvectors[(icomb+1)%3]->QyNorm(ih));
          EPconf->CorrelationEpProf(step,icomb,ih,3)->Fill(value,Qvectors[icomb]->QyNorm(ih)*Qvectors[(icomb+1)%3]->QyNorm(ih));
        }
      }

    }
  }
}


//_____________________________________________________________________
void AliEventPlaneManager::RecenterQvec(Float_t* values, Bool_t useFriendChain) {

  //
  // Recenter the detector event plane
  //


  Int_t bin=0;
  Int_t* var;
  Int_t dim; 
  Double_t fillValues[20];

  Int_t useStep=0;

  for(Int_t idet=0; idet<kNdetectors; idet++){
   AliEventPlaneConfiguration* EPconf = 0x0;
   TIter nextEPconf(fEventPlaneConfigurations[idet]);
   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
    if(!EPconf) continue;



    //if(!EPconf->doRecentering()) continue;
    if( EPconf->CalibrationStep()<2||!EPconf->doRecentering()) continue;

    //cout<<"!!! recenter"<<endl;


    useStep=0;
    if(EPconf->doChannelEqualization()) useStep=1;

    AliEventPlaneBinning* EPbinning =  EPconf->CalibrationBinning();
        
    AliEventPlaneQvector* qVector=0x0;
    TClonesArray* qvecList = EPconf->Qvectors();
    //TClonesArray* qvecList = GetQvectors(EPconf->GlobalIndex());
    for(Int_t ibin=0; ibin<qvecList->GetEntriesFast(); ibin++){
      qVector=static_cast<AliEventPlaneQvector*>(qvecList->At(ibin));
      if(!qVector) break;

      bin=qVector->Bin();

      var = EPbinning->Var();
      dim = EPbinning->Dim();

    //cout<<EPconf->EventPlaneDetectorName()<<"  "<<qVector->Qx(1)<<"   "<<qVector<<endl;
  


      for(Int_t iav=0; iav<dim; iav++) fillValues[iav] = values[var[iav]]; 





      Double_t axentries = EPconf->InputCalibrationHistogramE(useStep)->GetBinContent(EPconf->InputCalibrationHistogramE(useStep)->GetBin(fillValues));

	    for(Int_t ih=fgkEPMinHarmonics; ih<= fgkEPMaxHarmonics; ++ih) {
      
        if(qVector->CheckEventPlaneStatus(ih,AliEventPlaneQvector::kUndefined)) continue;

        if(axentries==1) {
	        qVector->SetQx(ih, qVector->Qx(ih));
   	      qVector->SetQy(ih, qVector->Qy(ih));
        }

        else{

	        qVector->SetQx( ih, 
	                           qVector->Qx(ih) - EPconf->InputCalibrationHistogramQ(useStep,ih,0)->GetBinContent(EPconf->InputCalibrationHistogramQ(useStep, ih,0)->GetBin(fillValues))/axentries);
	        qVector->SetQy( ih,  
	                           qVector->Qy(ih) - EPconf->InputCalibrationHistogramQ(useStep,ih,1)->GetBinContent(EPconf->InputCalibrationHistogramQ(useStep, ih,1)->GetBin(fillValues))/axentries);

	        qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kRecentered);

        }
        EPconf->CalibrationHistogramQ(2,ih,0)->Fill(fillValues,qVector->Qx(ih));
        EPconf->CalibrationHistogramQ(2,ih,1)->Fill(fillValues,qVector->Qy(ih));
        if(ih==fgkEPMinHarmonics) EPconf->CalibrationHistogramE(2)->Fill(fillValues);
      }

      // Perform u2n twist correction (TPC)
       //cout<<"recentering fail  "<<EPconf->CalibrationDetectorName(detector)<<"  "<<bin<<"   "<<qVector->Qx(1)<<"  "<<axentries<<"  "<<values[var[0]]<<"  "<<values[var[1]]<<"   "<<values[var[2]]<<endl;
      //if(!(TMath::Abs(qVector->Qx(1))<1.)&&axentries==0) cout<<"recentering fail  "<<EPconf->EventPlaneDetectorName()<<"  "<<qVector->Qx(1)<<"  "<<axentries<<"  "<<values[var[0]]<<"  "<<values[var[1]]<<"   "<<values[var[2]]<<endl;

    }
  }
  }
  //if(EPconf->TwistAndScalingMethod()==0&&EPconf->CalibrationStep()>2){
  //  U2nTwistQvec(values, AliEventPlaneVarManager::kCentVZERO);
  //}
  return;
}




////_____________________________________________________________________
//void AliEventPlaneManager::TwoDetectorCorrelationTwistQvec(Float_t* values, Int_t corpar) {
//
//  //
//  // Recenter the detector event plane
//  //
//
//
//  Int_t bin=0;
//  Int_t* var;
//  Int_t dim; 
//
//
//  for(Int_t idet=0; idet<kNdetectors; idet++){
//   AliEventPlaneConfiguration* EPconf = 0x0;
//   TIter nextEPconf(fEventPlaneConfigurations[idet]);
//   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
//    if(!EPconf) continue;
//
//    if( EPconf->CalibrationStep()<4) continue;
//    if( EPconf->TwistAndScalingMethod()!=2) continue;
//
//    TProfile * correlationProfiles[3][fgkEPMaxHarmonics][4];
//
//    for(Int_t ic=0; ic<3; ++ic) 
//	    for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  EPconf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist;
//
//	  for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) {
//
//      xx = correlationProfiles[0][ih-1][0]->GetBinContent(correlationProfiles[0][ih-1][0]->FindBin(values[corpar]));
//      yy = correlationProfiles[0][ih-1][3]->GetBinContent(correlationProfiles[0][ih-1][3]->FindBin(values[corpar]));
//      xy = correlationProfiles[0][ih-1][1]->GetBinContent(correlationProfiles[0][ih-1][1]->FindBin(values[corpar]));
//      yx = correlationProfiles[0][ih-1][2]->GetBinContent(correlationProfiles[0][ih-1][2]->FindBin(values[corpar]));
//
//
//      Lp = xy/(xx+TMath::Sqrt(xx*yy-xy*xy));
//      Ln = xy/(yy+TMath::Sqrt(xx*yy-xy*xy));
//
//
//      if(!(Lp>-99999999.&&Lp<99999999.)) continue;
//      if(!(Ln>-99999999.&&Ln<99999999.)) continue;
//
//
//      AliEventPlaneBinning* EPbinning =  EPconf->CalibrationBinning();
//          
//      AliEventPlaneQvector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(EPconf->GlobalIndex());
//      TIter nextQvector(qvecList);
//      while((qVector=static_cast<AliEventPlaneQvector*>(nextQvector()))) {
//        if(!qVector) continue;
//
//	      if(!qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kRecentered)) continue;
//      
//        Qx = qVector->Qx(ih);
//        Qy = qVector->Qy(ih);
//
//        QxTwist = (Qx-Ln*Qy)/(1-Ln*Lp);
//        QyTwist = (Qy-Lp*Qx)/(1-Ln*Lp);
//        
//	      qVector->SetQx( ih, QxTwist);
//     	  qVector->SetQy( ih, QyTwist);
//
//	      qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kDiagonalized);
//      }
//    }
//  }
//  }
//
//  return;
//}
//
//
//
////_____________________________________________________________________
//void AliEventPlaneManager::CorrelationRescalingQvec(Float_t* values, Int_t corpar) {
//
//  //
//  // Recenter the detector event plane
//  //
//
//
//  Int_t bin=0;
//  Int_t* var;
//  Int_t dim; 
//
//
//
//  for(Int_t idet=0; idet<kNdetectors; idet++){
//   AliEventPlaneConfiguration* EPconf = 0x0;
//   TIter nextEPconf(fEventPlaneConfigurations[idet]);
//   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
//    if(!EPconf) continue;
//
//
//    if( EPconf->CalibrationStep()<4) continue;
//    if( EPconf->TwistAndScalingMethod()!=1) continue;
//
//    TProfile * correlationProfiles[3][fgkEPMaxHarmonics][4];
//
//    for(Int_t ic=0; ic<3; ++ic) 
//	    for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  EPconf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxScaled, QyScaled;
//
//	  for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) {
//
//      xx = correlationProfiles[0][ih-1][0]->GetBinContent(correlationProfiles[0][ih-1][0]->FindBin(values[corpar]));
//      yy = correlationProfiles[0][ih-1][3]->GetBinContent(correlationProfiles[0][ih-1][3]->FindBin(values[corpar]));
//      xy = correlationProfiles[0][ih-1][1]->GetBinContent(correlationProfiles[0][ih-1][1]->FindBin(values[corpar]));
//      yx = correlationProfiles[0][ih-1][2]->GetBinContent(correlationProfiles[0][ih-1][2]->FindBin(values[corpar]));
//
//
//      Lp = xy/(xx+TMath::Sqrt(xx*yy-xy*xy));
//      Ln = xy/(yy+TMath::Sqrt(xx*yy-xy*xy));
//      Ap = TMath::Sqrt(2.*xx/(1+Lp*Lp));
//      An = TMath::Sqrt(2.*yy/(1+Ln*Ln));
//
//
//      if(!(Lp>-99999999.&&Lp<99999999.)) continue;
//      if(!(Ln>-99999999.&&Ln<99999999.)) continue;
//      if(!(Ap>-99999999.&&Ap<99999999.)) continue;
//      if(!(An>-99999999.&&An<99999999.)) continue;
//
//
//      AliEventPlaneBinning* EPbinning =  EPconf->CalibrationBinning();
//          
//      AliEventPlaneQvector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(EPconf->GlobalIndex());
//      TIter nextQvector(qvecList);
//      while((qVector=static_cast<AliEventPlaneQvector*>(nextQvector()))) {
//        if(!qVector) continue;
//
//	      if(!qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kDiagonalized)) continue;
//      
//        Qx = qVector->Qx(ih);
//        Qy = qVector->Qy(ih);
//
//        QxScaled = Qx / Ap;
//        QyScaled = Qy / An;
//        
//	      qVector->SetQx( ih, QxScaled);
//     	  qVector->SetQy( ih, QyScaled);
//
//	      qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kRescaled);
//      }
//    }
//  }
//  }
//
//  return;
//}
//
//
//
//
////_____________________________________________________________________
//void AliEventPlaneManager::CorrelationTwistAndRescalingQvec(Float_t* values, Int_t corpar) {
//
//  //
//  // Recenter the detector event plane
//  //
//
//
//  Int_t bin=0;
//  Int_t* var;
//  Int_t dim; 
//
//
// for(Int_t idet=0; idet<kNdetectors; idet++){
//   AliEventPlaneConfiguration* EPconf = 0x0;
//   TIter nextEPconf(fEventPlaneConfigurations[idet]);
//   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
//    if(!EPconf) continue;
//
//    if( EPconf->CalibrationStep()<4) continue;
//    if( EPconf->TwistAndScalingMethod()!=1) continue;
//
//    TProfile * correlationProfiles[3][fgkEPMaxHarmonics][4];
//
//    for(Int_t ic=1; ic<3; ++ic) 
//      for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  EPconf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//
//	  for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) {
//
//      xx = correlationProfiles[0][ih-1][0]->GetBinContent(correlationProfiles[0][ih-1][0]->FindBin(values[corpar]));
//      yy = correlationProfiles[0][ih-1][3]->GetBinContent(correlationProfiles[0][ih-1][3]->FindBin(values[corpar]));
//      xy = correlationProfiles[0][ih-1][1]->GetBinContent(correlationProfiles[0][ih-1][1]->FindBin(values[corpar]));
//      yx = correlationProfiles[0][ih-1][2]->GetBinContent(correlationProfiles[0][ih-1][2]->FindBin(values[corpar]));
//
//
//      Lp = xy/(xx+TMath::Sqrt(xx*yy-xy*xy));
//      Ln = xy/(yy+TMath::Sqrt(xx*yy-xy*xy));
//      Ap = TMath::Sqrt(2.*xx/(1+Lp*Lp));
//      An = TMath::Sqrt(2.*yy/(1+Ln*Ln));
//
//
//
//
//      if(!(Lp>-99999999.&&Lp<99999999.)) continue;
//      if(!(Ln>-99999999.&&Ln<99999999.)) continue;
//      if(!(Ap>-99999999.&&Ap<99999999.)) continue;
//      if(!(An>-99999999.&&An<99999999.)) continue;
//
//
//      AliEventPlaneBinning* EPbinning =  EPconf->CalibrationBinning();
//          
//      AliEventPlaneQvector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(EPconf->GlobalIndex());
//      TIter nextQvector(qvecList);
//      while((qVector=static_cast<AliEventPlaneQvector*>(nextQvector()))) {
//        if(!qVector) continue;
//
//        if(!qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kRecentered)) continue;
//      
//        Qx = qVector->Qx(ih);
//        Qy = qVector->Qy(ih);
//
//        QxTwist = (Qx-Ln*Qy)/(1-Ln*Lp);
//        QyTwist = (Qy-Lp*Qx)/(1-Ln*Lp);
//        
//        QxRescaled = QxTwist / Ap;
//        QyRescaled = QyTwist / An;
//        
//        qVector->SetQx( ih, QxRescaled);
//        qVector->SetQy( ih, QyRescaled);
//
//        qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kDiagonalized);
//        qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kRescaled);
//      }
//    }
//  }
// }
//
//  return;
//}
//
//
//
////_____________________________________________________________________
//void AliEventPlaneManager::ThreeDetectorCorrelationTPCTwistQvec(Float_t* values, Int_t corpar) {
//
//  //
//  // Recenter the detector event plane
//  //
//
//
//  Int_t bin=0;
//  Int_t* var;
//  Int_t dim; 
//
//
// for(Int_t idet=0; idet<kNdetectors; idet++){
//   if(idet==kTPC) continue;
//   AliEventPlaneConfiguration* EPconf = 0x0;
//   TIter nextEPconf(fEventPlaneConfigurations[idet]);
//   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
//    if(!EPconf) continue;
//
//    if( EPconf->CalibrationStep()<4) continue;
//    if(!EPconf->doTwist()) continue;
//    if( EPconf->TwistAndScalingMethod()!=2) continue;
//
//    TProfile * correlationProfiles[3][fgkNHarmonics][4];
//
//    for(Int_t ic=0; ic<3; ++ic) 
//      for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  EPconf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//    Double_t x1xt,y1yt,x1yt,x1x2,x2xt,x2yt;
//
//    for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) {
//
//      //x1xt = correlationProfiles[0][ih-1][0]->GetBinContent(correlationProfiles[0][ih-1][0]->FindBin(corpar));
//      //y1yt = correlationProfiles[0][ih-1][3]->GetBinContent(correlationProfiles[0][ih-1][3]->FindBin(corpar));
//      //x1yt = correlationProfiles[0][ih-1][1]->GetBinContent(correlationProfiles[0][ih-1][1]->FindBin(corpar));
//      x1xt = EPconf->CorrelationProfile(0, ih, 0)->GetBinContent(EPconf->CorrelationProfile(0, ih, 0)->FindBin(corpar));
//      y1yt = EPconf->CorrelationProfile(0, ih, 3)->GetBinContent(EPconf->CorrelationProfile(0, ih, 3)->FindBin(corpar));
//      x1yt = EPconf->CorrelationProfile(0, ih, 1)->GetBinContent(EPconf->CorrelationProfile(0, ih, 1)->FindBin(corpar));
//
//      Lp = x1yt/x1xt;
//      Ln = x1yt/y1yt;
//      //if(EPconf->EventPlaneDetectorName().EqualTo("VZEROA")) if(ih==2) cout<<corpar<<"  "<<x1yt<<"  "<<x1xt<<"   "<<y1yt<<endl;
//      //if(EPConf->EventPlaneDetectorName().EqualTo("VZEROA")) cout<<x1yt<<"  "<<x1xt<<"   "<<y1yt<<endl;
//      //if(ih==2) cout<<EPconf->EventPlaneDetectorName()<<"  "<<x1yt<<"  "<<x1xt<<"   "<<y1yt<<endl;
//
//      //cout<<EPconf->EventPlaneDetectorName()<<"  "<<Lp<<"  "<<Ln<<"  "<<x1xt<<"  "<<y1yt<<"  "<<x1yt<<"  "<<ih<<endl;
//      if(!(Lp>-99999999.&&Lp<99999999.)) continue;
//      if(!(Ln>-99999999.&&Ln<99999999.)) continue;
//
//      AliEventPlaneBinning* EPbinning =  EPconf->CalibrationBinning();
//          
//      AliEventPlaneQvector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(EPconf->GlobalIndex());
//      for(Int_t ibin=0; ibin<qvecList->GetEntriesFast(); ibin++){
//        qVector=static_cast<AliEventPlaneQvector*>(qvecList->At(ibin));
//
//        if(EPconf->EventPlaneDetectorName().Contains("NoRec")){
//          if(!qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kEqualized)) continue;}
//        else if(!qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kRecentered)) continue;
//
//        //if(!qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kRecentered)) continue;
//      
//        Qx = qVector->Qx(ih);
//        Qy = qVector->Qy(ih);
//
//        QxTwist = (Qx-Ln*Qy)/(1-Ln*Lp);
//        QyTwist = (Qy-Lp*Qx)/(1-Ln*Lp);
//        
//        qVector->SetQx( ih, QxTwist);
//        qVector->SetQy( ih, QyTwist);
//
//        qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kDiagonalized);
//      }
//    }
//  }
// }
//
//  return;
//}
//
//
////_____________________________________________________________________
//void AliEventPlaneManager::ThreeDetectorCorrelationTPCRescalingQvec(Float_t* values, Int_t corpar) {
//
//  //
//  // Recenter the detector event plane
//  //
//
//
//  Int_t bin=0;
//  Int_t* var;
//  Int_t dim; 
//
//
// for(Int_t idet=0; idet<kNdetectors; idet++){
//   if(idet==kTPC) continue;
//   AliEventPlaneConfiguration* EPconf = 0x0;
//   TIter nextEPconf(fEventPlaneConfigurations[idet]);
//   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
//    if(!EPconf) continue;
//
//    if( EPconf->CalibrationStep()<5) continue;
//    if(!EPconf->doScaling()) continue;
//    if( EPconf->TwistAndScalingMethod()!=2) continue;
//
//    //TProfile * correlationProfiles[3][fgkEPMaxHarmonics][4];
//
//    //for(Int_t ic=0; ic<3; ++ic) 
//    //  for(Int_t ih=1; ih<=fgkEPMaxHarmonics; ++ih) 
//    //   for(Int_t ip=0; ip<4; ++ip) {
//    //  correlationProfiles[ic][ih-1][ip] =  EPconf->CorrelationProfile(ic, ih, ip);
//    //}
//
//    Double_t Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//    Double_t x1xt,y1yt,x1yt,x1x2,x2xt,x2yt,x1y2,y2yt;
//
//    for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) {
//
//
//
//      x1xt = EPconf->CorrelationProfile(0, ih, 0)->GetBinContent(EPconf->CorrelationProfile(0, ih, 0)->FindBin(corpar));
//      y1yt = EPconf->CorrelationProfile(0, ih, 3)->GetBinContent(EPconf->CorrelationProfile(0, ih, 3)->FindBin(corpar));
//      x1yt = EPconf->CorrelationProfile(0, ih, 1)->GetBinContent(EPconf->CorrelationProfile(0, ih, 1)->FindBin(corpar));
//
//      x1x2 = EPconf->CorrelationProfile(2, ih, 0)->GetBinContent(EPconf->CorrelationProfile(2, ih, 0)->FindBin(corpar));
//
//      x2xt = EPconf->CorrelationProfile(1, ih, 0)->GetBinContent(EPconf->CorrelationProfile(1, ih, 0)->FindBin(corpar));
//      x2yt = EPconf->CorrelationProfile(1, ih, 1)->GetBinContent(EPconf->CorrelationProfile(1, ih, 1)->FindBin(corpar));
//
//
//
//
//
//      //x1xt = correlationProfiles[0][ih-1][0]->GetBinContent(correlationProfiles[0][ih-1][0]->FindBin(corpar));
//      //y1yt = correlationProfiles[0][ih-1][3]->GetBinContent(correlationProfiles[0][ih-1][3]->FindBin(corpar));
//      //x1yt = correlationProfiles[0][ih-1][1]->GetBinContent(correlationProfiles[0][ih-1][1]->FindBin(corpar));
//      ////y1xt = correlationProfiles[0][ih-1][2]->GetBinContent(correlationProfiles[0][ih-1][2]->FindBin(corpar));
//
//      //x1x2 = correlationProfiles[2][ih-1][0]->GetBinContent(correlationProfiles[2][ih-1][0]->FindBin(corpar));
//      ////y1y2 = correlationProfiles[2][ih-1][3]->GetBinContent(correlationProfiles[1][ih-1][3]->FindBin(corpar));
//      ////x1y2 = correlationProfiles[2][ih-1][1]->GetBinContent(correlationProfiles[1][ih-1][1]->FindBin(corpar));
//      ////y1x2 = correlationProfiles[2][ih-1][2]->GetBinContent(correlationProfiles[1][ih-1][2]->FindBin(corpar));
//
//      //x2xt = correlationProfiles[1][ih-1][0]->GetBinContent(correlationProfiles[1][ih-1][0]->FindBin(corpar));
//      ////y2yt = correlationProfiles[1][ih-1][3]->GetBinContent(correlationProfiles[2][ih-1][3]->FindBin(corpar));
//      //x2yt = correlationProfiles[1][ih-1][1]->GetBinContent(correlationProfiles[1][ih-1][1]->FindBin(corpar));
//      ////y2xt = correlationProfiles[1][ih-1][2]->GetBinContent(correlationProfiles[2][ih-1][2]->FindBin(corpar));
//
//
//      //Ap = TMath::Sqrt(2.*x1y2)*x1xt/TMath::Sqrt(x1xt*x2yt+x1yt*y2yt);
//      //An = TMath::Sqrt(2.*x1y2)*y1yt/TMath::Sqrt(x1xt*x2yt+x1yt*y2yt);
//
//      Ap = TMath::Sqrt(2.*x1x2)*x1xt/TMath::Sqrt(x1xt*x2xt+x1yt*x2yt);
//      An = TMath::Sqrt(2.*x1x2)*y1yt/TMath::Sqrt(x1xt*x2xt+x1yt*x2yt);
//
//      //Ap = TMath::Sqrt(2.*x1x2)*x1xt/TMath::Sqrt(x1xt*x2xt + x1yt*x2yt);
//      //An = TMath::Sqrt(2.*x1x2)*y1yt/TMath::Sqrt(x1xt*x2xt + x1yt*x2yt);
//
//      //if(EPconf->EventPlaneDetectorName().Contains("VZEROA")&&ih==2) cout<<"  "<<values[corpar]<<"  "<<ih<<"  "<<Ap<<"  "<<An<<endl;
//
//      if(!(Ap>-99999999.&&Ap<99999999.)) continue;
//      if(!(An>-99999999.&&An<99999999.)) continue;
//
//
//
//      AliEventPlaneBinning* EPbinning =  EPconf->CalibrationBinning();
//          
//      AliEventPlaneQvector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(EPconf->GlobalIndex());
//      for(Int_t ibin=0; ibin<qvecList->GetEntriesFast(); ibin++){
//        qVector=static_cast<AliEventPlaneQvector*>(qvecList->At(ibin));
//
//        //cout<<"hello  "<<qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kEqualized)<<"  "<<qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kUndefined)<<endl;
//	      if(!qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kDiagonalized)) continue;
//      
//        //cout<<"hey  "<<qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kEqualized)<<"  "<<qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kUndefined)<<endl;
//        Qx = qVector->Qx(ih);
//        Qy = qVector->Qy(ih);
//
//        QxRescaled = Qx / Ap;
//        QyRescaled = Qy / An;
//        
//	      qVector->SetQx( ih, QxRescaled);
//     	  qVector->SetQy( ih, QyRescaled);
//
//        //cout<<"hey 2  "<<qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kRescaled)<<endl;
//	      qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kRescaled);
//        //cout<<"hey 3  "<<qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kRescaled)<<endl;
//      }
//    }
//  }
// }
//
//  return;
//}
//
//
//
//
//
//
////_____________________________________________________________________
//void AliEventPlaneManager::ThreeDetectorCorrelationTPCTwistAndRescalingQvec(Float_t* values, Int_t corpar) {
//
//  //
//  // Recenter the detector event plane
//  //
//
//
//  Int_t bin=0;
//  Int_t* var;
//  Int_t dim; 
//
//
// for(Int_t idet=0; idet<kNdetectors; idet++){
//   if(idet==kTPC) continue;
//   AliEventPlaneConfiguration* EPconf = 0x0;
//   TIter nextEPconf(fEventPlaneConfigurations[idet]);
//   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
//    if(!EPconf) continue;
//
//    if( EPconf->CalibrationStep()<3) continue;
//    if( EPconf->TwistAndScalingMethod()!=2) continue;
//
//    TProfile * correlationProfiles[3][fgkEPMaxHarmonics][4];
//
//    for(Int_t ic=1; ic<3; ++ic) 
//	    for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  EPconf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//    Double_t x1xt,y1yt,x1yt,x1x2,x2xt,x2yt;
//
//	  for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) {
//
//      x1xt = correlationProfiles[0][ih-1][0]->GetBinContent(correlationProfiles[0][ih-1][0]->FindBin(values[corpar]));
//      y1yt = correlationProfiles[0][ih-1][3]->GetBinContent(correlationProfiles[0][ih-1][3]->FindBin(values[corpar]));
//      x1yt = correlationProfiles[0][ih-1][1]->GetBinContent(correlationProfiles[0][ih-1][1]->FindBin(values[corpar]));
//      //y1xt = correlationProfiles[0][ih-1][2]->GetBinContent(correlationProfiles[0][ih-1][2]->FindBin(values[corpar]));
//
//      x1x2 = correlationProfiles[1][ih-1][0]->GetBinContent(correlationProfiles[1][ih-1][0]->FindBin(values[corpar]));
//      //y1y2 = correlationProfiles[1][ih-1][3]->GetBinContent(correlationProfiles[1][ih-1][3]->FindBin(values[corpar]));
//      //x1y2 = correlationProfiles[1][ih-1][1]->GetBinContent(correlationProfiles[1][ih-1][1]->FindBin(values[corpar]));
//      //y1x2 = correlationProfiles[1][ih-1][2]->GetBinContent(correlationProfiles[1][ih-1][2]->FindBin(values[corpar]));
//
//      x2xt = correlationProfiles[2][ih-1][0]->GetBinContent(correlationProfiles[2][ih-1][0]->FindBin(values[corpar]));
//      //y2yt = correlationProfiles[2][ih-1][3]->GetBinContent(correlationProfiles[2][ih-1][3]->FindBin(values[corpar]));
//      x2yt = correlationProfiles[2][ih-1][1]->GetBinContent(correlationProfiles[2][ih-1][1]->FindBin(values[corpar]));
//      //y2xt = correlationProfiles[2][ih-1][2]->GetBinContent(correlationProfiles[2][ih-1][2]->FindBin(values[corpar]));
//
//
//      Lp = x1yt/x1xt;
//      Ln = x1yt/y1yt;
//      Ap = TMath::Sqrt(2.*x1x2)*x1xt/TMath::Sqrt(x1xt*x2xt+x1yt*x2yt);
//      An = TMath::Sqrt(2.*x1x2)*y1yt/TMath::Sqrt(x1xt*x2xt+x1yt*x2yt);
//
//
//      if(!(Lp>-99999999.&&Lp<99999999.)) continue;
//      if(!(Ln>-99999999.&&Ln<99999999.)) continue;
//      if(!(Ap>-99999999.&&Ap<99999999.)) continue;
//      if(!(An>-99999999.&&An<99999999.)) continue;
//
//
//      AliEventPlaneBinning* EPbinning =  EPconf->CalibrationBinning();
//          
//      AliEventPlaneQvector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(EPconf->GlobalIndex());
//      TIter nextQvector(qvecList);
//      while((qVector=static_cast<AliEventPlaneQvector*>(nextQvector()))) {
//        if(!qVector) continue;
//
//	      if(!qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kRecentered)) continue;
//      
//        Qx = qVector->Qx(ih);
//        Qy = qVector->Qy(ih);
//
//        QxTwist = (Qx-Ln*Qy)/(1-Ln*Lp);
//        QyTwist = (Qy-Lp*Qx)/(1-Ln*Lp);
//        
//        QxRescaled = QxTwist / Ap;
//        QyRescaled = QyTwist / An;
//        
//	      qVector->SetQx( ih, QxRescaled);
//     	  qVector->SetQy( ih, QyRescaled);
//
//	      qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kDiagonalized);
//	      qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kRescaled);
//      }
//    }
//  }
// }
//
//  return;
//}
//
//
//
//

//_____________________________________________________________________
void AliEventPlaneManager::U2nTwistQvec(Float_t* values, Int_t u2npar) {

  //
  // Recenter the detector event plane
  //


  Int_t bin=0;
  Int_t* var;
  Int_t maxHarmonic;
  Int_t dim; 
  Double_t fillValues[20];


 for(Int_t idet=0; idet<kNdetectors; idet++){
   AliEventPlaneConfiguration* EPconf = 0x0;
   TIter nextEPconf(fEventPlaneConfigurations[idet]);
   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
    if(!EPconf) break;

    if( EPconf->CalibrationStep()<3) continue;
    if(EPconf->TwistAndScalingMethod()!=0&&EPconf->TwistAndScalingMethod()!=1) continue;
    if( !EPconf->doTwist()) continue;
    maxHarmonic = fgkEPMaxHarmonics;



    AliEventPlaneBinning* EPbinning =  EPconf->CalibrationBinning();
        
    var = EPbinning->Var();
    dim = EPbinning->Dim();
    for(Int_t iav=0; iav<dim; iav++) fillValues[iav] = values[var[iav]]; 









    //if(EPconf->TwistAndScalingMethod()==1) maxHarmonic = fgkEPMaxHarmonics;
    //else maxHarmonic = 2*fgkEPMaxHarmonics;

    //TProfile * U2nProfiles[maxHarmonic][2];

	  //  for(Int_t ih=1; ih<=maxHarmonic; ++ih) 
    //   for(Int_t ip=0; ip<2; ++ip) {
    //  U2nProfiles[ih-1][ip] =  EPconf->U2nProfile(ih, ip);
    //}

    Double_t cos2n, sin2n, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, nentries;

	  for(Int_t ih=fgkEPMinHarmonics; ih<=(Int_t) fgkEPMaxHarmonics; ++ih) {


      sin2n = EPconf->InputU2nHistogram(ih,0)->GetBinContent(EPconf->InputU2nHistogram(ih,0)->GetBin(fillValues));
      cos2n = EPconf->InputU2nHistogram(ih,1)->GetBinContent(EPconf->InputU2nHistogram(ih,1)->GetBin(fillValues));
         nentries = EPconf->InputU2nHistogramE()->GetBinContent(EPconf->InputU2nHistogramE()->GetBin(fillValues));
      sin2n/=nentries;
      cos2n/=nentries;

      Ap = 1+cos2n;
      An = 1-cos2n;
      Lp = sin2n/Ap;
      Ln = sin2n/An;

      if(!(Lp>-99999999.&&Lp<99999999.)) continue;
      if(!(Ln>-99999999.&&Ln<99999999.)) continue;
      if(!(Ap>-99999999.&&Ap<99999999.)) continue;
      if(!(An>-99999999.&&An<99999999.)) continue;


          
      AliEventPlaneQvector* qVector=0x0;
      //TClonesArray* qvecList = GetQvectors(EPconf->GlobalIndex());
      TClonesArray* qvecList = EPconf->Qvectors();
      for(Int_t ibin=0; ibin<qvecList->GetEntriesFast(); ibin++){
        qVector=static_cast<AliEventPlaneQvector*>(qvecList->At(ibin));
        //if(!qVector) break;

	      if(!qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kRecentered)) continue;
      
        Qx = qVector->Qx(ih);
        Qy = qVector->Qy(ih);

        QxTwist = (Qx-Ln*Qy)/(1-Ln*Lp);
        QyTwist = (Qy-Lp*Qx)/(1-Ln*Lp);
        
	      qVector->SetQx( ih, QxTwist);
     	  qVector->SetQy( ih, QyTwist);

        //cout<<"Twist "<<EPconf->EventPlaneDetectorName()<<endl;
        //cout<<Qx<<"  "<<Qy<<"   "<<Lp<<"  "<<Ln<<"  "<<Ap<<"  "<<An<<endl;
	      qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kDiagonalized);
        EPconf->CalibrationHistogramQ(4,ih,0)->Fill(fillValues,qVector->Qx(ih));
        EPconf->CalibrationHistogramQ(4,ih,1)->Fill(fillValues,qVector->Qy(ih));
        if(ih==fgkEPMinHarmonics) EPconf->CalibrationHistogramE(4)->Fill(fillValues);
      }
    }
  }
 }

  return;
}



//_____________________________________________________________________
void AliEventPlaneManager::U2nRescalingQvec(Float_t* values, Int_t u2npar) {

  //
  // Recenter the detector event plane
  //


  Int_t bin=0;
  Int_t* var;
  Int_t maxHarmonic;
  Int_t dim; 


 for(Int_t idet=0; idet<kNdetectors; idet++){
   AliEventPlaneConfiguration* EPconf = 0x0;
   TIter nextEPconf(fEventPlaneConfigurations[idet]);
   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
    if(!EPconf) break;

    if( EPconf->CalibrationStep()<4) continue;
    //if( !EPconf->doScaling()) continue;
    if(EPconf->TwistAndScalingMethod()!=0&&EPconf->TwistAndScalingMethod()!=1) continue;
    if( !EPconf->doScaling()) continue;
    //if(EPconf->TwistAndScalingMethod()==1) maxHarmonic = fgkEPMaxHarmonics;
    //else maxHarmonic = 2*fgkEPMaxHarmonics;

    TProfile * U2nProfiles[maxHarmonic][2];

	    for(Int_t ih=1; ih<=maxHarmonic; ++ih) 
       for(Int_t ip=0; ip<2; ++ip) {
      U2nProfiles[ih-1][ip] =  EPconf->U2nProfile(ih, ip);
    }

    Double_t cos2n, sin2n, Lp, Ln, Ap, An, Qx, Qy, QxRescaling, QyRescaling;

	  for(Int_t ih=1; ih<=(Int_t) maxHarmonic/2.; ++ih) {

      cos2n = U2nProfiles[2*ih-1][1]->GetBinContent(U2nProfiles[2*ih-1][1]->FindBin(values[u2npar]));

      Ap = 1+cos2n;
      An = 1-cos2n;

      if(!(Ap>-99999999.&&Ap<99999999.)) continue;
      if(!(An>-99999999.&&An<99999999.)) continue;


      AliEventPlaneBinning* EPbinning =  EPconf->CalibrationBinning();
          
      AliEventPlaneQvector* qVector=0x0;
      //TClonesArray* qvecList = GetQvectors(EPconf->GlobalIndex());
      TClonesArray* qvecList = EPconf->Qvectors();
      for(Int_t ibin=0; ibin<qvecList->GetEntriesFast(); ibin++){
        qVector=static_cast<AliEventPlaneQvector*>(qvecList->At(ibin));
        //if(!qVector) break;

	      if(!qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kRecentered)) continue;
      
        Qx = qVector->Qx(ih);
        Qy = qVector->Qy(ih);

        QxRescaling = Qx / Ap;
        QyRescaling = Qy / An;
        
	      qVector->SetQx( ih, QxRescaling);
     	  qVector->SetQy( ih, QyRescaling);

	      qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kRescaled);
      }
    }
  }
 }

  return;
}


//
//
//
////_____________________________________________________________________
//void AliEventPlaneManager::U2nTwistAndRescalingQvec(Float_t* values, Int_t u2npar) {
//
//  //
//  // Recenter the detector event plane
//  //
//
//
//  Int_t bin=0;
//  Int_t* var;
//  Int_t dim; 
//
//
//
// for(Int_t idet=0; idet<kNdetectors; idet++){
//   AliEventPlaneConfiguration* EPconf = 0x0;
//   TIter nextEPconf(fEventPlaneConfigurations[idet]);
//   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
//    if(!EPconf) continue;
//
//
//    if( EPconf->CalibrationStep()<4) continue;
//    if( EPconf->TwistAndScalingMethod()!=0) continue;
//
//    TProfile * U2nProfiles[fgkNHarmonics*2][2];
//
//	    for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics*2; ++ih) 
//       for(Int_t ip=0; ip<2; ++ip) {
//      U2nProfiles[ih-1][ip] =  EPconf->U2nProfile(ih, ip);
//    }
//
//    Double_t cos2n, sin2n, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//
//	  for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) {
//
//      sin2n = U2nProfiles[ih-1][0]->GetBinContent(U2nProfiles[ih-1][0]->FindBin(values[u2npar]));
//      cos2n = U2nProfiles[ih-1][1]->GetBinContent(U2nProfiles[ih-1][1]->FindBin(values[u2npar]));
//
//      Ap = 1+cos2n;
//      An = 1-cos2n;
//      Lp = sin2n/Ap;
//      Ln = sin2n/An;
//
//      if(!(Lp>-99999999.&&Lp<99999999.)) continue;
//      if(!(Ln>-99999999.&&Ln<99999999.)) continue;
//      if(!(Ap>-99999999.&&Ap<99999999.)) continue;
//      if(!(An>-99999999.&&An<99999999.)) continue;
//
//
//      AliEventPlaneBinning* EPbinning =  EPconf->CalibrationBinning();
//          
//      AliEventPlaneQvector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(EPconf->GlobalIndex());
//      TIter nextQvector(qvecList);
//      while((qVector=static_cast<AliEventPlaneQvector*>(nextQvector()))) {
//        if(!qVector) continue;
//
//	      if(!qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kRecentered)) continue;
//      
//        Qx = qVector->Qx(ih);
//        Qy = qVector->Qy(ih);
//
//        QxTwist = (Qx-Ln*Qy)/(1-Ln*Lp);
//        QyTwist = (Qy-Lp*Qx)/(1-Ln*Lp);
//        
//        QxRescaled = QxTwist / Ap;
//        QyRescaled = QyTwist / An;
//        
//	      qVector->SetQx( ih, QxRescaled);
//     	  qVector->SetQy( ih, QyRescaled);
//
//	      qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kDiagonalized);
//	      qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kRescaled);
//      }
//    }
//  }
// }
//
//  return;
//}
//
//
//
//
//
////_____________________________________________________________________
//void AliEventPlaneManager::RotateQvec() {
//
//  //
//  // Recenter the detector event plane
//  //
//
//
//  Int_t bin=0;
//  Int_t* var;
//  Int_t dim; 
//
//  Double_t min=10.;
//  Double_t max=60.;
//
// for(Int_t idet=0; idet<kNdetectors; idet++){
//   if(idet==kTPC) continue;
//   AliEventPlaneConfiguration* EPconf = 0x0;
//   TIter nextEPconf(fEventPlaneConfigurations[idet]);
//   while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
//    if(!EPconf) continue;
//
//    if(!EPconf->doRotation()) continue;
//    if(EPconf->CalibrationStep()<3) continue;
//    //if( EPconf->CalibrationStep()>3) continue;
//
//    //cout<<EPconf->EventPlaneDetectorName()<<endl;
//    //cout<<EPconf->CorrelationProfile(0,2,0)<<endl;
//    //cout<<EPconf->CorrelationProfile(0,2,0)->GetEntries()<<endl;
//    TProfile * hXX = EPconf->CorrelationProfile(0,2,0);
//    TProfile * hXY = EPconf->CorrelationProfile(0,2,1);
//    TProfile * hYX = EPconf->CorrelationProfile(0,2,2);
//    TProfile * hYY = EPconf->CorrelationProfile(0,2,3);
//    hXX->GetXaxis()->SetRangeUser(min,max);
//    hXY->GetXaxis()->SetRangeUser(min,max);
//    hYX->GetXaxis()->SetRangeUser(min,max);
//    hYY->GetXaxis()->SetRangeUser(min,max);
//    Double_t XX = hXX->GetMean(2);
//    Double_t XY = hXY->GetMean(2);
//    Double_t YX = hYX->GetMean(2);
//    Double_t YY = hYY->GetMean(2);
//    Double_t eXX = hXX->GetMeanError(2);
//    Double_t eXY = hXY->GetMeanError(2);
//    Double_t eYX = hYX->GetMeanError(2);
//    Double_t eYY = hYY->GetMeanError(2);
//
//    Double_t dphi = TMath::ATan((XY-YX)/(XX+YY))/2.;
//    Double_t edenom2 = eXY*eXY+eYX*eYX;
//    Double_t enumer2 = eXX*eXX+eYY*eYY;
//
//    if(TMath::Sqrt((XY-YX)*(XY-YX)/edenom2)<2.) continue;
//    if(!(dphi<1000.)) continue;
//
//    Double_t equot  = TMath::Sqrt(enumer2/(XX+YY)/(XX+YY)+edenom2/(XY-YX)/(XY-YX))*((XY-YX)/(XX+YY));
//    Double_t edphi = equot/(1.+(XY-YX)/(XX+YY)*(XY-YX)/(XX+YY));
//    Double_t sigphi = TMath::Abs(dphi/edphi);
//
//    //cout<<EPconf->EventPlaneDetectorName()<<"  "<<dphi<<"  "<<edphi<<"  "<<sigphi<<"   "<<XX<<"  "<<XY<<"  "<<YX<<"  "<<YY<<endl;
//
//    Double_t Qx, Qy, Qmag, QxRotated, QyRotated;
//    Double_t x1yt,y1x2, x1x2, y1y2;
//
//    for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) {
//
//      AliEventPlaneBinning* EPbinning =  EPconf->CalibrationBinning();
//          
//      AliEventPlaneQvector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(EPconf->GlobalIndex());
//      for(Int_t ibin=0; ibin<qvecList->GetEntriesFast(); ibin++){
//        qVector=static_cast<AliEventPlaneQvector*>(qvecList->At(ibin));
//
//        if(!qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kRecentered)) continue;
//        if(qVector->CheckEventPlaneStatus(ih, AliEventPlaneQvector::kAligned)) continue;
//
//
//        Qx = qVector->Qx(ih);
//        Qy = qVector->Qy(ih);
//        Qmag = TMath::Sqrt(Qx*Qx+Qy*Qy);
//
//        QxRotated = Qx*TMath::Cos(((Double_t) ih)*dphi)+Qy*TMath::Sin(((Double_t) ih)*dphi);
//        QyRotated = Qy*TMath::Cos(((Double_t) ih)*dphi)-Qx*TMath::Sin(((Double_t) ih)*dphi);
//
//        
//        //cout<<EPconf->EventPlaneDetectorName()<<"  "<<Qx<<"   "<<QxRotated<<endl;
//        qVector->SetQx( ih, QxRotated);
//        qVector->SetQy( ih, QyRotated);
//
//        qVector->SetEventPlaneStatus(ih, AliEventPlaneQvector::kAligned);
//      }
//    }
//  }
// }
//
//  return;
//}
//



//
////_________________________________
//void AliEventPlaneManager::FillTPC(const AliESDEvent& esd, Float_t* values){
//
//  Int_t nTrack=-1;
//  AliEventPlaneConfiguration* EPconf = 0x0;
//  TClonesArray* epConfList = GetEventPlaneConfigurations(AliEventPlaneManager::kTPC);
//  TClonesArray& detector = *(GetReducedDetector(AliEventPlaneManager::kTPC));
//
//
//  for (Int_t iTrack = 0; iTrack < esd.GetNumberOfTracks(); ++iTrack)
//  {
//    AliESDtrack* esdTrack = esd.GetTrack(iTrack); //carefull do not modify it othwise  need to work with a copy 
//    AliESDtrack* track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(&esd),esdTrack->GetID());
//    if (!track) continue;
//
//    nTrack++;
//    
//    AliEventPlaneDetector *reducedDetector=new(detector[nTrack]) AliEventPlaneDetector();
//    AliEventPlaneVarManager::FillTrackInfo(track, values);
//
//    reducedDetector->SetPhi(track->Phi());
//    reducedDetector->SetX(TMath::Cos(track->Phi()));
//    reducedDetector->SetY(TMath::Sin(track->Phi()));
//    reducedDetector->SetWeight(1.);
//
//
//
//    Bool_t once = kTRUE;
//    Bool_t trackUsed = kFALSE;
//
//    for(Int_t iconf=0; iconf<epConfList->GetEntriesFast(); iconf++){
//      EPconf = (AliEventPlaneConfiguration*) epConfList->At(iconf);
//      if(EPconf->IsTrackSelected(values)){
//          reducedDetector->SetEventPlaneDetector( EPconf->LocalIndex() );
//          AliEventPlaneHistos::Instance()->FillHistClass("Tracks_"+EPconf->EventPlaneDetectorName(), values);
//          AliEventPlaneHistos::Instance()->FillHistClass("TrackQA_"+EPconf->EventPlaneDetectorName(), values);
//          //if(fRunLightWeight) {if(EPconf->CalibrationStep()==0) AliEventPlaneHistos::FillHistClass("TrackQA_"+EPconf->EventPlaneDetectorName(), values);}
//          //else AliEventPlaneHistos::FillHistClass("TrackQA_"+EPconf->EventPlaneDetectorName(), values, AliReducedVarManager::GetUsedVars());
//          trackUsed = kTRUE;
//      }
//    }
//
//    if(!trackUsed) nTrack--;
//
//  }
//
//
//
//}
//
//
//
////_________________________________
//void AliEventPlaneManager::FillTPC(AliReducedEvent* event, Float_t* values){
//
//  Int_t nTrack=-1;
//  AliReducedTrack* track=0x0;
//  AliEventPlaneConfiguration* EPconf = 0x0;
//  TClonesArray* trackList = event->GetTracks();
//  TClonesArray* epConfList = GetEventPlaneConfigurations(AliEventPlaneManager::kTPC);
//  TClonesArray& detector = *(GetReducedDetector(AliEventPlaneManager::kTPC));
//  //TIter nextTrack(trackList);
//  //while((track=static_cast<AliReducedTrack*>(nextTrack()))) {
//  //  if(!track) continue;
//
//  for(Int_t itrack=0; itrack<trackList->GetEntriesFast(); itrack++){
//  track = (AliReducedTrack*) trackList->At(itrack);
//   if(!track) continue;
//
//    nTrack++;
//    
//    AliEventPlaneDetector *reducedDetector=new(detector[nTrack]) AliEventPlaneDetector();
//
//    AliEventPlaneVarManager::FillTrackInfo(event, track, values);
//
//    reducedDetector->SetPhi(track->Phi());
//    reducedDetector->SetX(TMath::Cos(track->Phi()));
//    reducedDetector->SetY(TMath::Sin(track->Phi()));
//    reducedDetector->SetWeight(1.);
//
//
//
//    Bool_t once = kTRUE;
//    Bool_t trackUsed = kFALSE;
//
//    for(Int_t iconf=0; iconf<epConfList->GetEntriesFast(); iconf++){
//      EPconf = (AliEventPlaneConfiguration*) epConfList->At(iconf);
//      if(EPconf->IsTrackSelected(values)){
//          reducedDetector->SetEventPlaneDetector( EPconf->LocalIndex() );
//          AliEventPlaneHistos::Instance()->FillHistClass("Tracks_"+EPconf->EventPlaneDetectorName(), values);
//          AliEventPlaneHistos::Instance()->FillHistClass("TrackQA_"+EPconf->EventPlaneDetectorName(), values);
//          //if(fRunLightWeight) {if(EPconf->CalibrationStep()==0) AliEventPlaneHistos::FillHistClass("TrackQA_"+EPconf->EventPlaneDetectorName(), values);}
//          //else AliEventPlaneHistos::FillHistClass("TrackQA_"+EPconf->EventPlaneDetectorName(), values, AliReducedVarManager::GetUsedVars());
//          trackUsed = kTRUE;
//      }
//    }
//
//    if(!trackUsed) nTrack--;
//
//  }
//
//
//
//}
//
//
//
////_________________________________
//void AliEventPlaneManager::FillVZERO(TObject* event){
//
//  Double_t weight=0.;
//  const Double_t kX[8] = {0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268, 0.38268, 0.92388};    // cosines of the angles of the VZERO sectors (8 per ring)
//  const Double_t kY[8] = {0.38268, 0.92388, 0.92388, 0.38268, -0.38268, -0.92388, -0.92388, -0.38268};    // sines     -- " --
//  AliEventPlaneConfiguration* EPconf = 0x0;
//  AliVVZERO* vzero = 0x0;
//  AliReducedEvent* rvzero = 0x0;
//  Bool_t isVE=kTRUE;
//  if(event->IsA()==AliReducedEvent::Class()) isVE=kFALSE;
//  if(isVE) vzero= ((AliVEvent*)event)->GetVZEROData();
//  else rvzero=(AliReducedEvent*)event;
//
//  for(Int_t ich=0; ich<64; ich++){
//    weight=0.;
//    if(vzero) weight=vzero->GetMultiplicity(ich);
//    if(rvzero) weight=rvzero->MultChannelVZERO(ich);
//    if(weight<fVZEROminMult) weight=0.;
//
//    TClonesArray& detector = *(GetReducedDetector(AliEventPlaneManager::kVZERO));
//    AliEventPlaneDetector *reducedDetector=new(detector[ich]) AliEventPlaneDetector();
//    // copy vzero data and set respective coordinates
//    reducedDetector->SetX(kX[ich%8]);
//    reducedDetector->SetY(kY[ich%8]);
//    reducedDetector->SetWeight(weight);
//    reducedDetector->SetId(ich);
//    reducedDetector->SetBin(0);
//
//    //// set event plane subdetectors
//    TClonesArray* epConfList = GetEventPlaneConfigurations(AliEventPlaneManager::kVZERO);
//    for(Int_t iconf=0; iconf<epConfList->GetEntriesFast(); iconf++){
//    EPconf = (AliEventPlaneConfiguration*) epConfList->At(iconf);
//     if(!EPconf) continue;
//
//    //TIter nextEPconf(GetEventPlaneConfigurations(AliEventPlaneManager::kVZERO));
//    //while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
//    //  if(!EPconf) continue;
//      if(EPconf->UseChannel(ich)) reducedDetector->SetEventPlaneDetector(EPconf->LocalIndex());
//    }
//  }
//}
//
//
//
////_________________________________
//void AliEventPlaneManager::FillTZERO(TObject* event){
//
//  Double_t weight=0.0;
//  const Double_t kX[24] = {/* Cside */ 0.905348,0.571718,0.0848977,-0.424671,-0.82045,-0.99639,-0.905348,-0.571718,-0.0848977,0.424671,0.82045,0.99639, /* Aside */ 0.99995,0.870982,0.508635,0.00999978,-0.491315,-0.860982,-0.99995,-0.870982,-0.508635,-0.0100001,0.491315,0.860982};
//  const Double_t kY[24] = {/* Cside */ 0.424671,0.82045,0.99639,0.905348,0.571718,0.0848976,-0.424671,-0.82045,-0.99639,-0.905348,-0.571719,-0.0848975, /* Aside */ -0.00999983,0.491315,0.860982,0.99995,0.870982,0.508635,0.00999974,-0.491315,-0.860982,-0.99995,-0.870982,-0.508635};
//
//  AliEventPlaneConfiguration* EPconf = 0x0;
//
//  const AliESDTZERO* tzero = 0x0;
//  AliReducedEvent* rtzero = 0x0;
//  Bool_t isVE=kTRUE;
//  if(event->IsA()==AliReducedEvent::Class()) isVE=kFALSE;
//  if(isVE) tzero= ((AliESDEvent*)event)->GetESDTZERO();
//  else rtzero=(AliReducedEvent*)event;
//
//  for(Int_t ich=0; ich<24; ich++){
//
//    weight=0.;
//    if(tzero) weight=tzero->GetT0amplitude()[ich];
//    if(rtzero) weight=rtzero->AmplitudeTZEROch(ich);
//    if(weight<fTZEROminMult) weight=0.;
//
//    TClonesArray& detector = *(GetReducedDetector(AliEventPlaneManager::kTZERO));
//    AliEventPlaneDetector *reducedDetector=new(detector[ich]) AliEventPlaneDetector();
//    // copy tzero data and set respective coordinates
//    reducedDetector->SetX(kX[ich]);
//    reducedDetector->SetY(kY[ich]);
//    reducedDetector->SetWeight(weight);
//    reducedDetector->SetId(ich);
//    reducedDetector->SetBin(0);
//
//
//
//    // set event plane subdetectors
//    TIter nextEPconf(GetEventPlaneConfigurations(AliEventPlaneManager::kTZERO));
//    while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
//      if(!EPconf) continue;
//      if(EPconf->UseChannel(ich)) reducedDetector->SetEventPlaneDetector(EPconf->LocalIndex());
//    }
//  }
//}
//
//
////_________________________________
//void AliEventPlaneManager::FillZDC(AliReducedEvent* event){
//
//  const Double_t kX[10] = { /* Cside */ 0.0,  -1.75,  1.75, -1.75, 1.75, /* Aside */  0.0,  1.75, -1.75, 1.75, -1.75  };
//  const Double_t kY[10] = { /* Cside */ 0.0,  -1.75, -1.75,  1.75, 1.75, /* Aside */  0.0, -1.75, -1.75, 1.75,  1.75  };
//
//
//  AliEventPlaneConfiguration* EPconf = 0x0;
//
//  Int_t tower=-1;
//  for(Int_t ich=1; ich<10; ich++){
//    if(ich==5) continue;
//    tower++;
//
//    TClonesArray& detector = *(GetReducedDetector(AliEventPlaneManager::kZDC));
//    AliEventPlaneDetector *reducedDetector=new(detector[tower]) AliEventPlaneDetector();
//    // copy tzero data and set respective coordinates
//    reducedDetector->SetX(kX[ich]);
//    reducedDetector->SetY(kY[ich]);
//    reducedDetector->SetWeight(((event->EnergyZDCn(ich)>fZDCminMult) ? event->EnergyZDCn(ich) : 0.0));
//    reducedDetector->SetId(tower);
//    reducedDetector->SetBin(0);
//
//    // set event plane subdetectors
//    TIter nextEPconf(GetEventPlaneConfigurations(AliEventPlaneManager::kZDC));
//    while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
//      if(!EPconf) continue;
//      if(EPconf->UseChannel(ich)) reducedDetector->SetEventPlaneDetector(EPconf->LocalIndex());
//      }
//    }
//}
//
//
//////_________________________________
////void AliEventPlaneManager::FillFMD(AliReducedEvent* event){
////
////  AliEventPlaneConfiguration* EPconf = 0x0;
////  for(Int_t ifmd=0; ifmd<5; ifmd++){
////    AliReducedFMD* fmd = 0x0;
////    TClonesArray* fmdList = event->GetFMD(ifmd);
////    TIter nextChannel(fmdList);
////    Int_t ich=0;
////    Int_t fmddet;
////    if(ifmd==0) fmddet=AliEventPlaneManager::kFMD1;
////    if(ifmd==1) fmddet=AliEventPlaneManager::kFMD2I;
////    if(ifmd==2) fmddet=AliEventPlaneManager::kFMD2O;
////    if(ifmd==3) fmddet=AliEventPlaneManager::kFMD3I;
////    if(ifmd==4) fmddet=AliEventPlaneManager::kFMD3O;
////
////    while((fmd=static_cast<AliReducedFMD*>(nextChannel()))) {
////
////      if(fmd->Multiplicity()<fFMDminMult) continue;
////      TClonesArray& detector = *(GetReducedDetector(fmddet));
////      AliEventPlaneDetector *reducedDetector=new(detector[ich++]) AliEventPlaneDetector();
////      reducedDetector->SetX(TMath::Cos(fmd->Phi(ifmd, (Int_t) event->FMDnEtaSlices())));
////      reducedDetector->SetY(TMath::Sin(fmd->Phi(ifmd, (Int_t) event->FMDnEtaSlices())));
////      reducedDetector->SetWeight(fmd->Multiplicity());
////      reducedDetector->SetId(fmd->Id());
////      reducedDetector->SetBin(0);
////
////        // set event plane subdetectors
////      TIter nextEPconf(GetEventPlaneConfigurations(fmddet));
////      while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
////        if(!EPconf) continue;
////        if(EPconf->UseChannel(fmd->Id())) reducedDetector->SetEventPlaneDetector(EPconf->LocalIndex());
////      }
////    }
////  }
////
////}
//
//
////_________________________________
//void AliEventPlaneManager::FillFMD(AliReducedEvent* event){
//
//
//      AliEventPlaneConfiguration* EPconf = 0x0;
//      AliReducedFMD* fmd = 0x0;
//      TClonesArray* fmdList = event->GetFMD();
//      TIter nextChannel(fmdList);
//      Int_t ich=-1;
//      Int_t fmddet;
//      Int_t nEta=200;
//      Int_t nPhi=20;
//      Int_t iphi,ieta;
//      Double_t multtot=0;
//      Double_t mult;
//
//      fmd=static_cast<AliReducedFMD*>(nextChannel());
//
//      for(Int_t chfmd=0; chfmd<4000; chfmd++){
//          if(fmd&&fmd->Id()<chfmd) fmd=static_cast<AliReducedFMD*>(nextChannel());
//          TClonesArray& detector = *(GetReducedDetector(AliEventPlaneManager::kFMD));
//          AliEventPlaneDetector *reducedDetector=new(detector[chfmd]) AliEventPlaneDetector();
//          iphi=(chfmd-1)%nPhi;
//          ieta=(chfmd-iphi)/nPhi;
//          Double_t eta=-3.975+0.05*(ieta-1);
//          Double_t phi=TMath::Pi()/20.+TMath::Pi()/10.*iphi;
//
//          reducedDetector->SetPhi(phi);
//          reducedDetector->SetX(TMath::Cos(phi));
//          reducedDetector->SetY(TMath::Sin(phi));
//          if(fmd) {
//            if(fmd->Id()!=chfmd) mult=0;
//            else mult = fmd->Multiplicity();
//          }
//          else mult=0;
//          reducedDetector->SetWeight(mult);
//          reducedDetector->SetId(chfmd);
//          reducedDetector->SetBin(0);
//
//        // set event plane subdetectors
//        TIter nextEPconf(GetEventPlaneConfigurations(AliEventPlaneManager::kFMD));
//        while((EPconf=static_cast<AliEventPlaneConfiguration*>(nextEPconf()))) {
//          if(!EPconf) continue;
//              //if(EPconf->UseChannel(ich)) {
//                if(EPconf->EventPlaneDetectorName().Contains("FMDA")) {if(eta>1.0) reducedDetector->SetEventPlaneDetector(EPconf->LocalIndex());}
//                else if(EPconf->EventPlaneDetectorName().Contains("FMDC")) {if(eta<-1.0) reducedDetector->SetEventPlaneDetector(EPconf->LocalIndex());}
//                else reducedDetector->SetEventPlaneDetector(EPconf->LocalIndex());
//              //}
//            }
//         }
//}
//


////_______________________________________________________________________________
void AliEventPlaneManager::AddQvectors(Int_t q1, Int_t q2) {
//
//  Int_t filledBin=-1;
//
//  AliEventPlaneQvector* Qvector=0x0;
//  TClonesArray* QvecList = GetQvectors(q2);
//  TIter nextQvector(QvecList);
//  Int_t bin;
//  while(Qvector=static_cast<AliEventPlaneQvector*>(nextQvector())) {
//    bin=Qvector->Bin();
//
//    AliEventPlaneQvector* MasterQvector=0x0;
//    TClonesArray* MasterQvecList = GetQvectors(q1);
//    TIter nextMasterQvector(MasterQvecList);
//    Int_t Masterbin=-1;
//    while(MasterQvector=static_cast<AliEventPlaneQvector*>(nextMasterQvector())){
//      Masterbin = MasterQvector->Bin();
//      if(bin==Masterbin) break;
//    }
//
//
//    if(bin==Masterbin){
//      MasterQvector->SetMultiplicity(MasterQvector->Multiplicity()+Qvector->Multiplicity());
//      for(Int_t ih=1; ih<=fgkEPMaxHarmonics; ++ih){
//        MasterQvector->SetQx(ih, MasterQvector->Qx(ih)+Qvector->Qx(ih));
//        MasterQvector->SetQy(ih, MasterQvector->Qy(ih)+Qvector->Qy(ih));
//        for(Int_t iflag=0; iflag<AliEventPlaneQvector::kNMaxFlowFlags; iflag++){
//          if(!Qvector->CheckEventPlaneStatus(ih, (AliEventPlaneQvector::EventPlaneStatus) iflag)) continue;
//  	      MasterQvector->SetEventPlaneStatus(ih, (AliEventPlaneQvector::EventPlaneStatus) iflag);
//        }
//      }
//    }
//    else {
//
//      filledBin++;
//
//      TClonesArray& qvec = *(GetQvectors(fNEventPlaneDetectors+q1+1));
//      AliEventPlaneQvector *newQvector=new(qvec[filledBin]) AliEventPlaneQvector();
//
//      newQvector->SetMultiplicity(Qvector->Multiplicity());
//      newQvector->SetBin(Qvector->Bin());
//      for(Int_t ih=1; ih<=fgkEPMaxHarmonics; ++ih){
//        newQvector->SetQx(ih, Qvector->Qx(ih));
//        newQvector->SetQy(ih, Qvector->Qy(ih));
//        for(Int_t iflag=0; iflag<AliEventPlaneQvector::kNMaxFlowFlags; iflag++){
//          if(!Qvector->CheckEventPlaneStatus(ih, (AliEventPlaneQvector::EventPlaneStatus)iflag)) continue;
//  	      newQvector->SetEventPlaneStatus(ih, (AliEventPlaneQvector::EventPlaneStatus) iflag);
//        }
//      }
//    }
//  }
//
//  return;
}


