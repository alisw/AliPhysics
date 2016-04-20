/**************************************************************************************************
 *                                                                                                *
 * Package:       FlowVectorCorrections                                                           *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch                              *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com                             *
 *                Contributors are mentioned in the code where appropriate.                       *
 * Development:   2012-2015                                                                       *
 *                                                                                                *
 * This file is part of FlowVectorCorrections, a software package that corrects Q-vector          *
 * measurements for effects of nonuniform detector acceptance. The corrections in this package    *
 * are based on publication:                                                                      *
 *                                                                                                *
 *  [1] "Effects of non-uniform acceptance in anisotropic flow measurements"                      *
 *  Ilya Selyuzhenkov and Sergei Voloshin                                                         *
 *  Phys. Rev. C 77, 034904 (2008)                                                                *
 *                                                                                                *
 * The procedure proposed in [1] is extended with the following steps:                            *
 * (*) alignment correction between subevents                                                     *
 * (*) possibility to extract the twist and rescaling corrections                                 *
 *      for the case of three detector subevents                                                  *
 *      (currently limited to the case of two “hit-only” and one “tracking” detectors)            *
 * (*) (optional) channel equalization                                                            *
 * (*) flow vector width equalization                                                             *
 *                                                                                                *
 * FlowVectorCorrections is distributed under the terms of the GNU General Public License (GPL)   *
 * (https://en.wikipedia.org/wiki/GNU_General_Public_License)                                     *
 * either version 3 of the License, or (at your option) any later version.                        *
 *                                                                                                *
 **************************************************************************************************/
 
 
  
 

#include "AliQnCorrectionsConstants.h"
#include "AliQnCorrectionsManager.h"
#include "AliQnCorrectionsConfiguration.h"
#include "AliQnCorrectionsDataVector.h"
#include "AliQnCorrectionsHistograms.h"
#include "AliQnCorrectionsSteps.h"
#include <TMath.h>
#include <THashList.h>
#include <THashList.h>
#include <TClonesArray.h>
#include <TRandom3.h>
//#include <TArrayS.h>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

ClassImp(AliQnCorrectionsManager)


//_______________________________________________________________________________
AliQnCorrectionsManager::AliQnCorrectionsManager() :
  TObject(),
  fDataVectors(),
  fConfDataVectors(),
  fAliQnCorrectionsConfigurations(),
  fCorrectedQvectors(),
  fLastStep(),
  fInputHistograms(),
  fOutputHistograms(),
  //fDetectorIdMap(0x0),
  fNumberOfQnConfigurations(0),
  fNumberOfDetectors(0),
  fNumberOfQnConfigurationsForDetector(),
  fNdetectors(0),
  fIndex(),
  fCorrectionStep(0),
  fPassesRequired(0),
  fCalibrateByRun(kFALSE),
  fListInputHistogramsQnCorrections(0x0),
  fDataContainer(),
  fLabel(""),
  fProcessedFirstEvent(kFALSE),
  fQvecOutputList(),
  fSetFillTreeQnVectors(kFALSE),
  fSetFillHistogramsQA(kFALSE),
  fSetFillHistogramsQnCorrections(kFALSE),
  fListQnVectors(0x0),
  fTreeQnVectors(0x0),
  fListOutputHistogramsQnCorrections(0x0),
  fListHistogramsQA(0x0),
  fOutputHistogramsQnCorrectionsFile(0x0),
  fHistogramsQAFile(0x0),
  fTreeQnVectorsFile(0x0),
  fRunList("")
{
  //
  // default constructor
  //

  for(Int_t i=0; i<AliQnCorrectionsConstants::nDataVectors; ++i) fNumberOfQnConfigurationsForDetector[i]=0;
  for(Int_t i=0; i<AliQnCorrectionsConstants::nDataVectors; ++i) fDataVectors[i] = new TClonesArray("AliQnCorrectionsDataVector", 100000);
  for(Int_t i=0; i<AliQnCorrectionsConstants::nQnConfigurations; ++i) fConfDataVectors[i] = new TClonesArray("AliQnCorrectionsDataVector", 100000);
  for(Int_t i=0; i<AliQnCorrectionsConstants::nDataVectors; ++i) fAliQnCorrectionsConfigurations[i] = new TClonesArray("AliQnCorrectionsConfiguration", AliQnCorrectionsConstants::nQnConfigurations);


  for(Int_t iconf=0; iconf<AliQnCorrectionsConstants::nQnConfigurations; ++iconf){ 
    fInputHistograms[iconf]=0x0;
    fOutputHistograms[iconf]=0x0;
    fLastStep[iconf]=0;
    fIndex[iconf][0]=0x0;
    fIndex[iconf][1]=0x0;
  }

  for(Int_t i=0; i<AliQnCorrectionsConstants::nDataContainerVariables; i++) fDataContainer[i]=-999.;

  fListOutputHistogramsQnCorrections   = new THashList();
  fListOutputHistogramsQnCorrections->SetName("CalibrationHistograms");
  fListHistogramsQA = new THashList();
  fListHistogramsQA->SetName("CalibrationHistogramsQA");

  fListOutputHistogramsQnCorrections->SetOwner();
  fListHistogramsQA->SetOwner();

  fListQnVectors = new THashList();

}



//_______________________________________________________________________________
AliQnCorrectionsManager::AliQnCorrectionsManager(const AliQnCorrectionsManager &c) :
  TObject(),
  fDetectorIdMap(c.fDetectorIdMap),
  fNumberOfQnConfigurations(c.fNumberOfQnConfigurations),
  fNumberOfDetectors(c.fNumberOfDetectors),
  fNdetectors(c.fNdetectors),
  fCorrectionStep(c.fCorrectionStep),
  fPassesRequired(c.fPassesRequired),
  fCalibrateByRun(c.fCalibrateByRun),
  fListInputHistogramsQnCorrections(c.fListInputHistogramsQnCorrections),
  fLabel(c.fLabel),
  fProcessedFirstEvent(c.fProcessedFirstEvent),
  fSetFillTreeQnVectors(c.fSetFillTreeQnVectors),
  fSetFillHistogramsQA(c.fSetFillHistogramsQA),
  fSetFillHistogramsQnCorrections(c.fSetFillHistogramsQnCorrections),
  fListQnVectors(c.fListQnVectors),
  fTreeQnVectors(c.fTreeQnVectors),
  fListOutputHistogramsQnCorrections(c.fListOutputHistogramsQnCorrections),
  fListHistogramsQA(c.fListHistogramsQA),
  fOutputHistogramsQnCorrectionsFile(c.fOutputHistogramsQnCorrectionsFile),
  fHistogramsQAFile(c.fHistogramsQAFile),
  fTreeQnVectorsFile(c.fTreeQnVectorsFile),
  fRunList(c.fRunList)
{
  //
  // copy constructor
  //

  for(Int_t i=0; i<AliQnCorrectionsConstants::nDataVectors; ++i) fNumberOfQnConfigurationsForDetector[i]=c.fNumberOfQnConfigurationsForDetector[i];
  for(Int_t i=0; i<AliQnCorrectionsConstants::nDataVectors; ++i) fDataVectors[i] = c.fDataVectors[i];
  for(Int_t i=0; i<AliQnCorrectionsConstants::nQnConfigurations; ++i) fConfDataVectors[i] = c.fConfDataVectors[i] ;
  for(Int_t i=0; i<AliQnCorrectionsConstants::nDataVectors; ++i) fAliQnCorrectionsConfigurations[i] = c.fAliQnCorrectionsConfigurations[i];


  for(Int_t iconf=0; iconf<AliQnCorrectionsConstants::nQnConfigurations; ++iconf){ 
    fInputHistograms[iconf]=c.fInputHistograms[iconf];
    fOutputHistograms[iconf]=c.fOutputHistograms[iconf];
    fLastStep[iconf]=c.fLastStep[iconf];
    fIndex[iconf][0]=c.fIndex[iconf][0];
    fIndex[iconf][1]=c.fIndex[iconf][1];
  }

  for(Int_t i=0; i<AliQnCorrectionsConstants::nDataContainerVariables; i++){
    fDataContainer[i]= c.fDataContainer[i];
  }

}




//_______________________________________________________________________________
AliQnCorrectionsManager & AliQnCorrectionsManager::operator=(const AliQnCorrectionsManager &c) {
  if (this == &c) return *this;
  else {
    fDetectorIdMap=c.fDetectorIdMap;
    fNumberOfQnConfigurations=c.fNumberOfQnConfigurations;
    fNumberOfDetectors=c.fNumberOfDetectors;
    fNdetectors=c.fNdetectors;
    fCorrectionStep=c.fCorrectionStep;
    fPassesRequired=c.fPassesRequired;
    fCalibrateByRun=c.fCalibrateByRun;
    fListInputHistogramsQnCorrections=c.fListInputHistogramsQnCorrections;
    fLabel=c.fLabel;
    fProcessedFirstEvent=c.fProcessedFirstEvent;
    fSetFillTreeQnVectors=c.fSetFillTreeQnVectors;
    fSetFillHistogramsQA=c.fSetFillHistogramsQA;
    fSetFillHistogramsQnCorrections=c.fSetFillHistogramsQnCorrections;
    fListQnVectors=c.fListQnVectors;
    fTreeQnVectors=c.fTreeQnVectors;
    fListOutputHistogramsQnCorrections=c.fListOutputHistogramsQnCorrections;
    fListHistogramsQA=c.fListHistogramsQA;
    fOutputHistogramsQnCorrectionsFile=c.fOutputHistogramsQnCorrectionsFile;
    fHistogramsQAFile=c.fHistogramsQAFile;
    fTreeQnVectorsFile=c.fTreeQnVectorsFile;
    fRunList=c.fRunList;

    for(Int_t i=0; i<AliQnCorrectionsConstants::nDataVectors; ++i) fNumberOfQnConfigurationsForDetector[i]=c.fNumberOfQnConfigurationsForDetector[i];
    for(Int_t i=0; i<AliQnCorrectionsConstants::nDataVectors; ++i) fDataVectors[i] = c.fDataVectors[i];
    for(Int_t i=0; i<AliQnCorrectionsConstants::nQnConfigurations; ++i) fConfDataVectors[i] = c.fConfDataVectors[i] ;
    for(Int_t i=0; i<AliQnCorrectionsConstants::nDataVectors; ++i) fAliQnCorrectionsConfigurations[i] = c.fAliQnCorrectionsConfigurations[i];


    for(Int_t iconf=0; iconf<AliQnCorrectionsConstants::nQnConfigurations; ++iconf){ 
      fInputHistograms[iconf]=c.fInputHistograms[iconf];
      fOutputHistograms[iconf]=c.fOutputHistograms[iconf];
      fLastStep[iconf]=c.fLastStep[iconf];
      fIndex[iconf][0]=c.fIndex[iconf][0];
      fIndex[iconf][1]=c.fIndex[iconf][1];
    }

    for(Int_t i=0; i<AliQnCorrectionsConstants::nDataContainerVariables; i++){
      fDataContainer[i]= c.fDataContainer[i];
    }

      return *this;
  }
    
}










//____________________________________________________________________________
AliQnCorrectionsManager::~AliQnCorrectionsManager()
{
  //
  // De-Constructor
  //
  ClearEvent();
  for(Int_t i=0; i<fNumberOfDetectors; ++i) if(fDataVectors[i]) fDataVectors[i]->Clear("C");
  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; ++iconf) {
    for(Int_t istep=0; istep<AliQnCorrectionsConstants::nCorrectionSteps; istep++) if(fCorrectedQvectors[iconf][istep]) fCorrectedQvectors[iconf][istep]->Clear("C");
    if(fConfDataVectors[iconf]) fConfDataVectors[iconf]->Clear("C");
  }
  for(Int_t i=0; i<AliQnCorrectionsConstants::nDataVectors; ++i) if(fAliQnCorrectionsConfigurations[i]) fAliQnCorrectionsConfigurations[i]->Clear("C");
}


//_____________________________________________________________________________
void AliQnCorrectionsManager::ClearEvent() {
  //
  // clear the event
  //
  for(Int_t i=0; i<fNumberOfDetectors; ++i) if(fDataVectors[i]) fDataVectors[i]->Clear("C");
  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; ++iconf) {
    for(Int_t istep=0; istep<AliQnCorrectionsConstants::nCorrectionSteps; istep++) if(fCorrectedQvectors[iconf][istep]) if(fCorrectedQvectors[iconf][istep]->At(0)) (static_cast<AliQnCorrectionsQnVector*>(fCorrectedQvectors[iconf][istep]->At(0)))->Reset();
    if(fConfDataVectors[iconf]) fConfDataVectors[iconf]->Clear("C");
  }
  for(Int_t i=0; i<AliQnCorrectionsConstants::nDataContainerVariables; i++) fDataContainer[i]=-999.;
  fListQnVectors->Clear();
}

//_____________________________________________________________________________
//Bool_t AliQnCorrectionsManager::SetCalibrationFile(TFile* inputFile){
//    if(inputFile){
//      if(inputFile->GetListOfKeys()->GetEntries()>0){
//        //fListInputHistogramsQnCorrections = (THashList*)((TKey*)inputFile->GetListOfKeys()->At(0))->ReadObj(); return (fListInputHistogramsQnCorrections ? kTRUE : kFALSE);} //TODO make "CalibrationHistograms" name customizable
//        THashList* list = (THashList*)((TKey*)inputFile->GetListOfKeys()->At(0))->ReadObj();
//        if(list){
//          for(Int_t i=0; i<list->GetEntries(); i++){
//            fListInputHistogramsQnCorrections->Add(new TObject(*(list->At(i)));
//          }
//        }
//      
//      return (fListInputHistogramsQnCorrections ? kTRUE : kFALSE);} //TODO make "CalibrationHistograms" name customizable
//    }
//}

//_____________________________________________________________________________
void AliQnCorrectionsManager::Process() {

  if(!fProcessedFirstEvent){
    //Initialize();
    InitializeCalibrationHistograms();
    PrintFrameworkInformation();
  }

  ApplyQnCorrections();

  FillHistograms();

  if(fSetFillTreeQnVectors){
    for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
      fQvecOutputList[iconf]=((AliQnCorrectionsQnVector*) CorrectedQnVector(iconf)->At(0));
    }
    fTreeQnVectors->Fill();
  }

//   ClearEvent();

  fProcessedFirstEvent=kTRUE;

}



////_____________________________________________________________________________
//void AliQnCorrectionsManager::ApplyQnCorrections() {
//
//
//  AliQnCorrectionsConfiguration* QnConf = 0x0;
//
//  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
//    AliQnCorrectionsConfiguration* QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
//    if(!QnConf) continue;
//
//    CallStepEqualizeDataVectorWeights(QnConf);
//
//    BuildQnVectors(QnConf, kFALSE);   // Get Q vectors with raw detector information
//
//    BuildQnVectors(QnConf, kTRUE);    // Get Q vectors with equalized detector information
//
//    CallStepRecenterQnVector(QnConf);
//
//
//  }
//
//}



//_____________________________________________________________________________
void AliQnCorrectionsManager::ApplyQnCorrections() {


  AliQnCorrectionsConfiguration* QnConf = 0x0;


  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
    if(!QnConf) continue;


    BuildQnVectors(QnConf, kFALSE);
    if(QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)) {CallStepEqualizeDataVectorWeights(QnConf);BuildQnVectors(QnConf, kTRUE);}
    
    if(QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kRecentering))             CallStepRecenterQnVector(QnConf);
    
    if(QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kAlignment))               RotateQvec(QnConf);

    if(QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kTwist))                   CallStepTwistAndRescaleQnVector(QnConf);



  }
}




//_____________________________________________________________________
void AliQnCorrectionsManager::FillHistograms() {

  AliQnCorrectionsConfiguration* QnConf;
  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);

    if(fSetFillHistogramsQnCorrections){
      if(QnConf->IsFillHistogram(AliQnCorrectionsConstants::kDataVectorEqualization)||QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)){
        FillHistogramsWeights(QnConf);
        //cout<<"go into fill "<<QnConf->QnConfigurationName()<<"  "<<QnConf->IsFillHistogram(AliQnCorrectionsConstants::kRecentering)<<endl;
        if(QnConf->IsFillHistogram(AliQnCorrectionsConstants::kRecentering))    {FillHistogramsMeanQ(QnConf, (Int_t) AliQnCorrectionsConstants::kDataVectorEqualization);}
      }
      if(QnConf->IsFillHistogram(AliQnCorrectionsConstants::kRecentering)&&!QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization))   
                                                                            {FillHistogramsMeanQ(QnConf, (Int_t) AliQnCorrectionsConstants::kPass0);}
      if(QnConf->IsFillHistogram(AliQnCorrectionsConstants::kAlignment))       {FillHistogramsQnAlignment(QnConf);} 
      if(QnConf->IsFillHistogram(AliQnCorrectionsConstants::kTwist))           {FillHistogramsU2n(QnConf); FillHistogramsQnCorrelations(QnConf);}
    }

    // Fill QA correction histograms

    if(fSetFillHistogramsQA){
                                                                                              //FillHistogramsQnCorrelationsQA(QnConf);
                                                                                           FillHistogramsMeanQ_QA(QnConf);
       if(QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kAlignment))                   FillHistogramsQnAlignmentQA(QnConf); 
       if((QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kTwist)||QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kRescaling))&&QnConf->GetTwistAndRescalingMethod()==2)
         FillHistogramsQnCorrelationsQA(QnConf);
                                                                                              
    }

  }

}

//_____________________________________________________________________________
void AliQnCorrectionsManager::WriteQnVectorsToList() {


  AliQnCorrectionsConfiguration* QnConf=0x0;
  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
    if(!QnConf) continue;
    TClonesArray* qvec = CorrectedQnVector(iconf);
    if(qvec) {
      qvec->SetName(QnConf->QnConfigurationName());
      fListQnVectors->Add(qvec);
    }
    for(Int_t istep=0; istep<AliQnCorrectionsConstants::nCorrectionSteps; istep++) {
      TClonesArray* qvecAll = CorrectedQnVector(iconf,istep);
      if(qvecAll) {
        qvecAll->SetName(Form("%s_%d",QnConf->QnConfigurationName().Data(),istep));
        fListQnVectors->Add(qvecAll);
      }

    }
  }
}

//_____________________________________________________________________________
void AliQnCorrectionsManager::Initialize() {
  //
  // Initialize before looping over data
  // Sets the correlation detector indices (this prevents having to do many string operations), creates the necessary histograms, free space reserved for unused detectors
  //



  AliQnCorrectionsConfiguration* QnConf = 0x0;
  AliQnCorrectionsConfiguration* QnConf2 = 0x0;
  TString detStrCor1 = "";
  TString detStrCor2 = "";
  TString detStrCor3 = "";
  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
    if(!QnConf) continue;
    detStrCor1 = QnConf->QnConfigurationCorrelationName(0);
    detStrCor2 = QnConf->QnConfigurationCorrelationName(1);
    detStrCor3 = QnConf->QnConfigurationCorrelationName(2);
    QnConf2 = 0x0;
    for(Int_t iconf2=0; iconf2<fNumberOfQnConfigurations; iconf2++){
      QnConf2 = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf2);
      if(!QnConf2) continue;
      if(detStrCor1.EqualTo(QnConf2->QnConfigurationName())) QnConf->SetQnConfigurationCorrelationIndex(0,QnConf2->GlobalIndex());
      if(detStrCor2.EqualTo(QnConf2->QnConfigurationName())) QnConf->SetQnConfigurationCorrelationIndex(1,QnConf2->GlobalIndex());
      if(detStrCor3.EqualTo(QnConf2->QnConfigurationName())) QnConf->SetQnConfigurationCorrelationIndex(2,QnConf2->GlobalIndex());
    }
  }

  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
    if(!QnConf) continue;
    fInputHistograms[iconf]  = new AliQnCorrectionsHistograms();
    fOutputHistograms[iconf] = new AliQnCorrectionsHistograms();
    fOutputHistograms[iconf]->CreateCalibrationHistograms(QnConf);
    for(Int_t istep=0; istep<AliQnCorrectionsConstants::kNcorrectionSteps; istep++)
      {fCorrectedQvectors[iconf][istep] = new TClonesArray("AliQnCorrectionsQnVector", 10);
       TClonesArray& arr = *(fCorrectedQvectors[iconf][istep]);
       if(istep>1) continue;
       AliQnCorrectionsQnVector* q = new(arr[0]) AliQnCorrectionsQnVector();
       if(QnConf->GetTwistAndRescalingMethod()==0) q->SetMaximumHarmonic(QnConf->MaximumHarmonic()*2);
       else q->SetMaximumHarmonic(QnConf->MaximumHarmonic()*2);
      }
  }

  fNumberOfDetectors=fDetectorIdMap.size();
  



  if(fSetFillTreeQnVectors){

    fTreeQnVectors = new TTree("Qvectors","Qvector values");
    fTreeQnVectors->SetDirectory(0);

    for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
      QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
      if(!QnConf) continue;
      //fQvecOutputList[QnConf->GlobalIndex()]->At(0)  = new AliQnCorrectionsQnVector(); //CorrectedQnVector(iconf,istep);
      //fQvecOutputList[QnConf->GlobalIndex()]->SetName(QnConf->QnConfigurationName());
      //fTreeQnVectors->Branch(Form("%s_%d",QnConf->QnConfigurationName().Data(),istep),&fQvecOutputList[QnConf->GlobalIndex()],256000,1);
      fQvecOutputList[iconf]=new AliQnCorrectionsQnVector();
      fTreeQnVectors->Branch(QnConf->QnConfigurationName(),&fQvecOutputList[iconf],256000,1);
    }

  }

  //InitializeCalibrationHistograms();
}


//_______________________________________________________________________________
void AliQnCorrectionsManager::SetFileRunLabels(TString path){
  //
  // Add runlist to the framework
  //

  ifstream inBuf;
  inBuf.open(path);
  Int_t index = 0;
  while(inBuf.good()) {
    Char_t str[512];
    inBuf.getline(str,512,'\n');
    AddRunLabel(str);
  }

}



//_______________________________________________________________________________
void AliQnCorrectionsManager::AddDataVector( Int_t detectorId, Double_t phi, Double_t weight, Int_t dataVectorId) {
  //
  // Add datavector to TClonesArray of DataVectors and set its bits
  //

  detectorId=fDetectorIdMap[detectorId]-1;
  Bool_t useForQn;

  AliQnCorrectionsConfiguration* QnConf;
  for(Int_t iconf=0; iconf<fNumberOfQnConfigurationsForDetector[detectorId]; iconf++){
    QnConf=static_cast<AliQnCorrectionsConfiguration*>(fAliQnCorrectionsConfigurations[detectorId]->At(iconf));
    if(!QnConf) continue;

    useForQn=kTRUE;

    if(QnConf->Cuts())         useForQn=QnConf->PassCuts(fDataContainer); 
    if(QnConf->ChannelList())  if(useForQn) useForQn=QnConf->UseChannel(dataVectorId);

    if(useForQn) {
      TClonesArray& dataVectors = *(fConfDataVectors[QnConf->GlobalIndex()]);
      AliQnCorrectionsDataVector *dataVector=new(dataVectors[dataVectors.GetEntriesFast()]) AliQnCorrectionsDataVector();
      dataVector->SetPhi(phi);
      dataVector->SetWeight(weight);
      dataVector->SetId(dataVectorId);
    }
  } 
}







//_________________________________________________________________
Bool_t AliQnCorrectionsManager::BuildQnVectors(AliQnCorrectionsConfiguration* QnConf, Bool_t useEqualizedWeights){
  //
  // Construct the event plane for a specified detector
  //

  //if(useEqualizedWeights){ 
  //  if(QnConf->GetCorrectionStep(fCorrectionStep)<1) return kFALSE;
  //  if(!QnConf->doChannelEqualization()) return kFALSE;
  //}
    
  AliQnCorrectionsQnVector* QnVector = static_cast<AliQnCorrectionsQnVector*>(fCorrectedQvectors[QnConf->GlobalIndex()][(Int_t) useEqualizedWeights]->At(0));
  fLastStep[QnConf->GlobalIndex()]=(Int_t) useEqualizedWeights;


  if(!useEqualizedWeights)   AliQnCorrectionsDataVector::FillQvector(fConfDataVectors[QnConf->GlobalIndex()], QnVector);
  else                       AliQnCorrectionsDataVector::FillQvector(fConfDataVectors[QnConf->GlobalIndex()], QnVector, QnConf->GetDataVectorEqualizationMethod());

  if(QnVector->N()==0) {               // If detector is empty
    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) QnVector->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kUndefined);
    return kTRUE;
  }
  else{
    Int_t maxHar = QnConf->MaximumHarmonic();
    if(QnConf->GetTwistAndRescalingMethod()==0) maxHar*=2;
    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=maxHar; ++ih){
      if(!useEqualizedWeights)   QnVector->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kPass0);
      else                       QnVector->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kDataVectorEqualization);
    }
  }


  if(QnConf->GetQnNormalizationMethod()==0)        {QnVector->SetQoverSquareRootOfM();QnVector->SetQvectorNormalization(0);}
  else if(QnConf->GetQnNormalizationMethod()==1)   {QnVector->SetQoverM();QnVector->SetQvectorNormalization(1);}
  else if(QnConf->GetQnNormalizationMethod()==2)   {QnVector->Normalize();QnVector->SetQvectorNormalization(2);}

  return kTRUE;


}

//_____________________________________________________________________
Bool_t AliQnCorrectionsManager::CallStepEqualizeDataVectorWeights(AliQnCorrectionsConfiguration* QnConf) {
  //
  // Calibrate channel multiplicities
  //

  //if(!QnConf->doChannelEqualization()||QnConf->GetCorrectionStep(fCorrectionStep)<1) return kFALSE;

  Double_t fillValues[20];

  const Int_t* var = QnConf->EqualizationBinning()->Var();
  const Int_t  dim = QnConf->EqualizationBinning()->Dim();

  for(Int_t iav=0; iav<(dim-1); iav++) fillValues[iav] = fDataContainer[var[iav]];
  AliQnCorrectionsSteps::CalibrateDataVector(GetConfDataVectors(QnConf->GlobalIndex()), QnConf,  fInputHistograms[QnConf->GlobalIndex()], fillValues);

  return kTRUE;
}


//_____________________________________________________________________
Bool_t AliQnCorrectionsManager::CallStepRecenterQnVector(AliQnCorrectionsConfiguration* QnConf) {

  //
  // Recenter the detector event plane
  //


  Int_t bin=-1;
  Double_t fillValues[20];

  Int_t useStep=0;


  //if( QnConf->GetCorrectionStep(fCorrectionStep)<2||!QnConf->doRecentering()) return kFALSE;

  useStep=0;
  if(QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)) useStep=1;

  AliQnCorrectionsAxes* EPbinning =  QnConf->GetRecenteringAxes();


  const Int_t* var = EPbinning->Var();
  const Int_t  dim = EPbinning->Dim();

  for(Int_t iav=0; iav<dim; iav++) fillValues[iav] = fDataContainer[var[iav]]; 

  bin=fInputHistograms[QnConf->GlobalIndex()]->CalibrationHistogramE(useStep)->GetBin(fillValues);

  AliQnCorrectionsQnVector* QvectorIn=0x0;
  if(useStep==0) QvectorIn=static_cast<AliQnCorrectionsQnVector*>(fCorrectedQvectors[QnConf->GlobalIndex()][0]->At(0));
  else           QvectorIn=static_cast<AliQnCorrectionsQnVector*>(fCorrectedQvectors[QnConf->GlobalIndex()][1]->At(0));

  TClonesArray& arr = *(fCorrectedQvectors[QnConf->GlobalIndex()][AliQnCorrectionsConstants::kRecentering]);
  AliQnCorrectionsQnVector* QvectorOut = new(arr[0]) AliQnCorrectionsQnVector(*QvectorIn);

  //AliQnCorrectionsQnVector* QvectorOut=static_cast<AliQnCorrectionsQnVector*>(fCorrectedQvectors[QnConf->GlobalIndex()][2]->At(0));
  fLastStep[QnConf->GlobalIndex()]=AliQnCorrectionsConstants::kRecentering;

  AliQnCorrectionsSteps::RecenterQvec( QvectorIn, QvectorOut, fInputHistograms[QnConf->GlobalIndex()], bin, useStep, QnConf->MinimumHarmonic(), QnConf->MaximumHarmonic());


  return kTRUE;

}






//_____________________________________________________________________
void AliQnCorrectionsManager::RotateQvec(AliQnCorrectionsConfiguration* QnConf) {

  //
  // Align Q-vectors
  //


  Double_t fillValues[20];

  const Int_t* var = QnConf->GetAlignmentAxes()->Var();
  const Int_t  dim = QnConf->GetAlignmentAxes()->Dim();
  for(Int_t iav=0; iav<dim; iav++) fillValues[iav] = fDataContainer[var[iav]];

  Int_t iconf=QnConf->GlobalIndex();
  Int_t bin = fInputHistograms[iconf]->GetRotationHistogramE(0)->GetBin(fillValues);



    Double_t XX = fInputHistograms[iconf]->GetRotationHistogram(0,0)->GetBinContent(bin);
    Double_t YY = fInputHistograms[iconf]->GetRotationHistogram(0,1)->GetBinContent(bin);
    Double_t XY = fInputHistograms[iconf]->GetRotationHistogram(0,2)->GetBinContent(bin);
    Double_t YX = fInputHistograms[iconf]->GetRotationHistogram(0,3)->GetBinContent(bin);
    Double_t eXX = fInputHistograms[iconf]->GetRotationHistogram(0,0)->GetBinError(bin);
    Double_t eYY = fInputHistograms[iconf]->GetRotationHistogram(0,1)->GetBinError(bin);
    Double_t eXY = fInputHistograms[iconf]->GetRotationHistogram(0,2)->GetBinError(bin);
    Double_t eYX = fInputHistograms[iconf]->GetRotationHistogram(0,3)->GetBinError(bin);
    Double_t n = fInputHistograms[iconf]->GetRotationHistogramE(0)->GetBinContent(bin);

    eXX =  TMath::Sqrt((eXX*eXX/n-(XX/n*XX/n))/n);
    eYY =  TMath::Sqrt((eYY*eYY/n-(YY/n*YY/n))/n);
    eXY =  TMath::Sqrt((eXY*eXY/n-(XY/n*XY/n))/n);
    eYX =  TMath::Sqrt((eYX*eYX/n-(YX/n*YX/n))/n);

    Double_t dphi = -TMath::ATan((XY-YX)/(XX+YY))/QnConf->AlignmentHarmonic();
    Double_t edenom2 = eXY*eXY+eYX*eYX;
    Double_t enumer2 = eXX*eXX+eYY*eYY;

    AliQnCorrectionsQnVector* QvectorIn= static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(iconf)->At(0));
    TClonesArray &arr = *(fCorrectedQvectors[QnConf->GlobalIndex()][AliQnCorrectionsConstants::kAlignment]);
    AliQnCorrectionsQnVector* QvectorRotated = new(arr[0]) AliQnCorrectionsQnVector(*QvectorIn);

    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) {QvectorRotated->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kAlignment);}

    if(TMath::Sqrt((XY-YX)*(XY-YX)/edenom2)<2.) return;
    //if(!(dphi<1.)) return;

    Double_t equot  = TMath::Sqrt(enumer2/(XX+YY)/(XX+YY)+edenom2/(XY-YX)/(XY-YX))*((XY-YX)/(XX+YY));
    Double_t edphi = equot/(1.+(XY-YX)/(XX+YY)*(XY-YX)/(XX+YY));
    Double_t sigphi = TMath::Abs(dphi/edphi);
    //dphi=-dphi;

    //cout<<QnConf->QnConfigurationName()<<"  "<<dphi<<"  "<<edphi<<"  "<<sigphi<<"   "<<XX<<"  "<<XY<<"  "<<YX<<"  "<<YY<<endl;

    Double_t Qx, Qy, Qmag, QxRotated, QyRotated;
    Double_t x1yt,y1x2, x1x2, y1y2;

    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) {

        Qx = QvectorIn->Qx(ih);
        Qy = QvectorIn->Qy(ih);

        QvectorRotated->SetQx(ih,Qx*TMath::Cos(((Double_t) ih)*dphi)+Qy*TMath::Sin(((Double_t) ih)*dphi));
        QvectorRotated->SetQy(ih,Qy*TMath::Cos(((Double_t) ih)*dphi)-Qx*TMath::Sin(((Double_t) ih)*dphi));


        QvectorRotated->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kAlignment);
    }

  fLastStep[QnConf->GlobalIndex()]=AliQnCorrectionsConstants::kAlignment;
  return;
}



//_____________________________________________________________________
void AliQnCorrectionsManager::FillHistogramsWeights(AliQnCorrectionsConfiguration* QnConf) {


  Double_t fillValues[20];


  TClonesArray* dataVectorArray = GetConfDataVectors(QnConf->GlobalIndex());

  const Int_t* var = QnConf->EqualizationBinning()->Var();
  const Int_t  dim = QnConf->EqualizationBinning()->Dim();
  for(Int_t iav=0; iav<(dim-1); iav++) fillValues[iav] = fDataContainer[var[iav]];

  Int_t globalIndex=QnConf->GlobalIndex();
  AliQnCorrectionsDataVector* dataVector;

  for(Int_t idata=0; idata<dataVectorArray->GetEntriesFast(); idata++){
    dataVector = static_cast<AliQnCorrectionsDataVector*>(dataVectorArray->At(idata));

    //TIter nextEntry(dataVectorArray);
    //while((dataVector=static_cast<AliQnCorrectionsDataVector*>(nextEntry()))) {
    //if(!dataVector) continue;

    //if(dataVector->Weight()<1E-6) continue;
    //if(!dataVector->CheckEventPlaneDetector(localIndex)) continue;

    fillValues[dim-1] = dataVector->Id();
    //bin=fOutputHistograms[globalIndex]->EqualizationHistogramM(0)->GetBin(fillValues);

    fOutputHistograms[globalIndex]->EqualizationHistogramM(0)->Fill(fillValues,dataVector->Weight());
    fOutputHistograms[globalIndex]->EqualizationHistogramE(0)->Fill(fillValues,1.);
    if(!fSetFillHistogramsQA) continue;
    if(!QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)) continue;
    fOutputHistograms[globalIndex]->EqualizationHistogramM(1)->Fill(fillValues,dataVector->Weight(0));
    fOutputHistograms[globalIndex]->EqualizationHistogramE(1)->Fill(fillValues,1.);
    fOutputHistograms[globalIndex]->EqualizationHistogramM(2)->Fill(fillValues,dataVector->Weight(1));
    fOutputHistograms[globalIndex]->EqualizationHistogramE(2)->Fill(fillValues,1.);

  }
}

//______________________________________
void AliQnCorrectionsManager::FillHistogramsQnAlignment(AliQnCorrectionsConfiguration* QnConf){

  Int_t iconf=QnConf->GlobalIndex();

  AliQnCorrectionsQnVector* Qvecs[2];
  AliQnCorrectionsConfiguration* QnConfRef=GetQnConfiguration(QnConf->QnConfigurationCorrelationIndex(2));

  Double_t fillValues[20];

  const Int_t* var = QnConf->GetAlignmentAxes()->Var();
  const Int_t  dim = QnConf->GetAlignmentAxes()->Dim();
  for(Int_t iav=0; iav<dim; iav++) fillValues[iav] = fDataContainer[var[iav]];

  if(QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kPass0))                  Qvecs[0] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(iconf,AliQnCorrectionsConstants::kPass0)->At(0));
  if(QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)) Qvecs[0] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(iconf,AliQnCorrectionsConstants::kDataVectorEqualization)->At(0));
  if(QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kRecentering))            Qvecs[0] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(iconf,AliQnCorrectionsConstants::kRecentering)->At(0));
  
  Qvecs[1] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(QnConf->QnConfigurationCorrelationIndex(2))->At(0));

  //cout<<QnConf->QnConfigurationName()<<"  "<<QnConf->QnConfigurationCorrelationIndex(2)<<"  "<<GetQnConfiguration(QnConf->QnConfigurationCorrelationIndex(2))->QnConfigurationName()<<endl;

  THnF* entr = fOutputHistograms[iconf]->GetRotationHistogramE(0);
  Int_t bin=entr->GetBin(fillValues);

  Double_t qx1,qx2,qy1,qy2;

  Int_t har=QnConf->AlignmentHarmonic();
  if(Qvecs[0]->CheckQnVectorStatus(har,AliQnCorrectionsConstants::kUndefined)||Qvecs[1]->CheckQnVectorStatus(har,AliQnCorrectionsConstants::kUndefined)) return;
  qx1=Qvecs[0]->Qx(har); qx2=Qvecs[1]->Qx(har); qy1=Qvecs[0]->Qy(har); qy2=Qvecs[1]->Qy(har);
  fOutputHistograms[iconf]->GetRotationHistogram(0,0)->AddBinContent(bin,qx1*qx2);
  fOutputHistograms[iconf]->GetRotationHistogram(0,1)->AddBinContent(bin,qy1*qy2);
  fOutputHistograms[iconf]->GetRotationHistogram(0,2)->AddBinContent(bin,qx1*qy2);
  fOutputHistograms[iconf]->GetRotationHistogram(0,3)->AddBinContent(bin,qy1*qx2);
  fOutputHistograms[iconf]->GetRotationHistogram(0,0)->AddBinError2(bin,qx1*qx2*qx1*qx2);
  fOutputHistograms[iconf]->GetRotationHistogram(0,1)->AddBinError2(bin,qy1*qy2*qy1*qy2);
  fOutputHistograms[iconf]->GetRotationHistogram(0,2)->AddBinError2(bin,qx1*qy2*qx1*qy2);
  fOutputHistograms[iconf]->GetRotationHistogram(0,3)->AddBinError2(bin,qy1*qx2*qy1*qx2);
  entr->AddBinContent(bin,1.0);
  Int_t entries = entr->GetEntries();
  entr->SetEntries(++entries);

}


//______________________________________
void AliQnCorrectionsManager::FillHistogramsQnAlignmentQA(AliQnCorrectionsConfiguration* QnConf){

  Int_t iconf=QnConf->GlobalIndex();

  AliQnCorrectionsQnVector* Qvecs[2];
  AliQnCorrectionsConfiguration* QnConfRef=GetQnConfiguration(QnConf->QnConfigurationCorrelationIndex(2));

  Double_t fillValues[20];

  const Int_t* var = QnConf->GetAlignmentAxes()->Var();
  const Int_t  dim = QnConf->GetAlignmentAxes()->Dim();
  for(Int_t iav=0; iav<dim; iav++) fillValues[iav] = fDataContainer[var[iav]];


  Qvecs[0] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(iconf,AliQnCorrectionsConstants::kAlignment)->At(0));
  Qvecs[1] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(QnConf->QnConfigurationCorrelationIndex(2))->At(0));


  //cout<<QnConf->QnConfigurationName()<<"  "<<QnConf->QnConfigurationCorrelationIndex(2)<<"  "<<GetQnConfiguration(QnConf->QnConfigurationCorrelationIndex(2))->QnConfigurationName()<<endl;

  Int_t bin=fOutputHistograms[iconf]->GetRotationHistogramE(1)->GetBin(fillValues);

  Double_t qx1,qx2,qy1,qy2;

  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){ 
    if(ih>QnConfRef->MaximumHarmonic()) continue;
    if(Qvecs[0]->CheckQnVectorStatus(ih,AliQnCorrectionsConstants::kUndefined)||Qvecs[1]->CheckQnVectorStatus(ih,AliQnCorrectionsConstants::kUndefined)) continue;
    qx1=Qvecs[0]->Qx(ih); qx2=Qvecs[1]->Qx(ih); qy1=Qvecs[0]->Qy(ih); qy2=Qvecs[1]->Qy(ih);
    fOutputHistograms[iconf]->GetRotationHistogram(1,0)->AddBinContent(bin,qx1*qx2);
    fOutputHistograms[iconf]->GetRotationHistogram(1,1)->AddBinContent(bin,qy1*qy2);
    fOutputHistograms[iconf]->GetRotationHistogram(1,2)->AddBinContent(bin,qx1*qy2);
    fOutputHistograms[iconf]->GetRotationHistogram(1,3)->AddBinContent(bin,qy1*qx2);
    fOutputHistograms[iconf]->GetRotationHistogram(1,0)->AddBinError2(bin,qx1*qx2*qx1*qx2);
    fOutputHistograms[iconf]->GetRotationHistogram(1,1)->AddBinError2(bin,qy1*qy2*qy1*qy2);
    fOutputHistograms[iconf]->GetRotationHistogram(1,2)->AddBinError2(bin,qx1*qy2*qx1*qy2);
    fOutputHistograms[iconf]->GetRotationHistogram(1,3)->AddBinError2(bin,qy1*qx2*qy1*qx2);
  }
  fOutputHistograms[iconf]->GetRotationHistogramE(1)->AddBinContent(bin,1.0);

}



//______________________________________
void AliQnCorrectionsManager::FillHistogramsQnCorrelations(AliQnCorrectionsConfiguration* QnConfig){

  Float_t value;
  Int_t index1 = QnConfig->QnConfigurationCorrelationIndex(0);
  Int_t index2 = QnConfig->QnConfigurationCorrelationIndex(1);
  if(index1==-1||index2==-1) return;
  Int_t iconf=QnConfig->GlobalIndex();

  AliQnCorrectionsQnVector* Qvecs[3];
  AliQnCorrectionsConfiguration* QnConf[3]={0x0};
  QnConf[0]=QnConfig;

  Int_t dim=QnConf[0]->GetAlignmentAxes()->Dim();

  QnConf[1]=GetQnConfiguration(index1);
  QnConf[2]=GetQnConfiguration(index2);

  Qvecs[1] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(index1)->At(0));
  
  if(QnConf[0]->IsRequestedCorrection(AliQnCorrectionsConstants::kAlignment))
           Qvecs[0] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(iconf,AliQnCorrectionsConstants::kAlignment)->At(0));
  else if(QnConf[0]->IsRequestedCorrection(AliQnCorrectionsConstants::kRecentering))
           Qvecs[0] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(iconf,AliQnCorrectionsConstants::kRecentering)->At(0));
  else if(QnConf[0]->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization))
           Qvecs[0] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(iconf,AliQnCorrectionsConstants::kDataVectorEqualization)->At(0));
  else     Qvecs[0] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(iconf,AliQnCorrectionsConstants::kPass0)->At(0));

  if(QnConf[2]->IsRequestedCorrection(AliQnCorrectionsConstants::kAlignment))
           Qvecs[2] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(index2,AliQnCorrectionsConstants::kAlignment)->At(0));
  else if(QnConf[2]->IsRequestedCorrection(AliQnCorrectionsConstants::kRecentering))
           Qvecs[2] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(index2,AliQnCorrectionsConstants::kRecentering)->At(0));
  else if(QnConf[2]->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization))
           Qvecs[2] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(index2,AliQnCorrectionsConstants::kDataVectorEqualization)->At(0));
  else     Qvecs[2] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(index2,AliQnCorrectionsConstants::kPass0)->At(0));



  for(Int_t icomb=0; icomb<3; ++icomb){ 
    for(Int_t iaxis=0; iaxis<dim; iaxis++){
    value=fDataContainer[QnConf[0]->GetAlignmentAxes()->Var(iaxis)];
      for(Int_t ih=QnConf[icomb]->MinimumHarmonic(); ih<=QnConf[icomb]->MaximumHarmonic(); ++ih){ 
        if(ih>QnConf[(icomb+1)%3]->MaximumHarmonic()) continue;
        if(Qvecs[icomb]->CheckQnVectorStatus(ih,AliQnCorrectionsConstants::kUndefined)||Qvecs[(icomb+1)%3]->CheckQnVectorStatus(ih,AliQnCorrectionsConstants::kUndefined)) continue;
        fOutputHistograms[iconf]->CorrelationProf(0,icomb,ih,0,iaxis)->Fill(value,Qvecs[icomb]->Qx(ih)*Qvecs[(icomb+1)%3]->Qx(ih));
        fOutputHistograms[iconf]->CorrelationProf(0,icomb,ih,1,iaxis)->Fill(value,Qvecs[icomb]->Qy(ih)*Qvecs[(icomb+1)%3]->Qx(ih));
        fOutputHistograms[iconf]->CorrelationProf(0,icomb,ih,2,iaxis)->Fill(value,Qvecs[icomb]->Qx(ih)*Qvecs[(icomb+1)%3]->Qy(ih));
        fOutputHistograms[iconf]->CorrelationProf(0,icomb,ih,3,iaxis)->Fill(value,Qvecs[icomb]->Qy(ih)*Qvecs[(icomb+1)%3]->Qy(ih));
      }
    }
  }

}


//______________________________________
void AliQnCorrectionsManager::FillHistogramsQnCorrelationsQA(AliQnCorrectionsConfiguration* QnConfig){

  Float_t value;
  Int_t index1 = QnConfig->QnConfigurationCorrelationIndex(0);
  Int_t index2 = QnConfig->QnConfigurationCorrelationIndex(1);
  if(index1==-1||index2==-1) return;
  Int_t iconf=QnConfig->GlobalIndex();

  AliQnCorrectionsQnVector* Qvecs[3];
  AliQnCorrectionsConfiguration* QnConf[3]={0x0};
  QnConf[0]=QnConfig;

  Int_t dim=QnConf[0]->GetAlignmentAxes()->Dim();

  QnConf[1]=GetQnConfiguration(index1);
  QnConf[2]=GetQnConfiguration(index2);

  for(Int_t icor=0; icor<2; icor++){

    Int_t correction;
    if(icor==0) correction=AliQnCorrectionsConstants::kTwist;
    if(icor==1) correction=AliQnCorrectionsConstants::kRescaling;
  
  if(QnConf[0]->IsRequestedCorrection(correction))
           Qvecs[0] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(iconf,correction)->At(0));
  if(QnConf[2]->IsRequestedCorrection(correction))
           Qvecs[2] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(index2,correction)->At(0));

  Qvecs[1] = static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(index1)->At(0));
  



  for(Int_t icomb=0; icomb<3; ++icomb){ 
    for(Int_t iaxis=0; iaxis<dim; iaxis++){
    value=fDataContainer[QnConf[0]->GetAlignmentAxes()->Var(iaxis)];
      for(Int_t ih=QnConf[icomb]->MinimumHarmonic(); ih<=QnConf[icomb]->MaximumHarmonic(); ++ih){ 
        if(ih>QnConf[(icomb+1)%3]->MaximumHarmonic()) continue;
        if(Qvecs[icomb]->CheckQnVectorStatus(ih,AliQnCorrectionsConstants::kUndefined)||Qvecs[(icomb+1)%3]->CheckQnVectorStatus(ih,AliQnCorrectionsConstants::kUndefined)) continue;
        fOutputHistograms[iconf]->CorrelationProf(icor+1,icomb,ih,0,iaxis)->Fill(value,Qvecs[icomb]->Qx(ih)*Qvecs[(icomb+1)%3]->Qx(ih));
        fOutputHistograms[iconf]->CorrelationProf(icor+1,icomb,ih,1,iaxis)->Fill(value,Qvecs[icomb]->Qy(ih)*Qvecs[(icomb+1)%3]->Qx(ih));
        fOutputHistograms[iconf]->CorrelationProf(icor+1,icomb,ih,2,iaxis)->Fill(value,Qvecs[icomb]->Qx(ih)*Qvecs[(icomb+1)%3]->Qy(ih));
        fOutputHistograms[iconf]->CorrelationProf(icor+1,icomb,ih,3,iaxis)->Fill(value,Qvecs[icomb]->Qy(ih)*Qvecs[(icomb+1)%3]->Qy(ih));
      }
    }
  }
  }

}


//______________________________________
void AliQnCorrectionsManager::FillHistogramsMeanQ(AliQnCorrectionsConfiguration* QnConf, Int_t step){
        //cout<<"requested fill "<<QnConf->QnConfigurationName()<<"  "<<endl;

  Double_t fillValues[20];
  Int_t iconf=QnConf->GlobalIndex();

  const Int_t* var = QnConf->GetRecenteringAxes()->Var();
  const Int_t  dim = QnConf->GetRecenteringAxes()->Dim();
  
  
  for(Int_t iav=0; iav<dim; iav++) fillValues[iav] = fDataContainer[var[iav]];
  Int_t bin=fOutputHistograms[iconf]->CalibrationHistogramE(0)->GetBin(fillValues);

  AliQnCorrectionsQnVector* Qvec=static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(iconf,step)->At(0));


  Bool_t once=kTRUE;

  Int_t maxHar = QnConf->MaximumHarmonic();
  if(QnConf->GetTwistAndRescalingMethod()==0) maxHar*=2;
  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=maxHar; ++ih){

  Double_t weight=1.;
  //if(QnConf->GetQnNormalizationMethod()==1) weight = Qvec->SumOfWeights();

    if(Qvec->CheckQnVectorStatus(ih, step)){
      fOutputHistograms[iconf]->CalibrationHistogramQ(step,ih,0)->AddBinContent(bin,Qvec->Qx(ih)*weight);
      fOutputHistograms[iconf]->CalibrationHistogramQ(step,ih,1)->AddBinContent(bin,Qvec->Qy(ih)*weight);
      fOutputHistograms[iconf]->CalibrationHistogramQ(step,ih,0)->AddBinError2(bin,Qvec->Qx(ih)*weight*Qvec->Qx(ih)*weight);
      fOutputHistograms[iconf]->CalibrationHistogramQ(step,ih,1)->AddBinError2(bin,Qvec->Qy(ih)*weight*Qvec->Qy(ih)*weight);
      if(once) {
        THnF* entr = fOutputHistograms[iconf]->CalibrationHistogramE(step);
        entr->AddBinContent(bin,weight);
        Int_t entries = entr->GetEntries();
        entr->SetEntries(++entries);
        once=kFALSE;
      }

    }
  }

  return;

}


//______________________________________
void AliQnCorrectionsManager::FillHistogramsMeanQ_QA(AliQnCorrectionsConfiguration* QnConf){


  Double_t fillValues[20];
  Int_t iconf=QnConf->GlobalIndex();


  const Int_t* var = QnConf->GetRecenteringAxes()->Var();
  const Int_t  dim = QnConf->GetRecenteringAxes()->Dim();
  for(Int_t iav=0; iav<dim; iav++) fillValues[iav] = fDataContainer[var[iav]];

  Int_t bin=fOutputHistograms[iconf]->CalibrationHistogramE(0)->GetBin(fillValues);

  for(Int_t istep=0; istep<=fCorrectionStep; istep++){
    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)) {if(istep==1) continue;}
    else if(istep==0) continue;

    AliQnCorrectionsQnVector* Qvec=static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(iconf,istep)->At(0));
    if(!Qvec) continue;
    if(Qvec->N()==0) continue;


    Bool_t once=kTRUE;

    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){

      if(Qvec->CheckQnVectorStatus(ih, istep)){
        fOutputHistograms[iconf]->CalibrationHistogramQ(istep,ih,0)->AddBinContent(bin,Qvec->Qx(ih));
        fOutputHistograms[iconf]->CalibrationHistogramQ(istep,ih,1)->AddBinContent(bin,Qvec->Qy(ih));
        if(once) {
          THnF* entr = fOutputHistograms[iconf]->CalibrationHistogramE(istep);
          entr->AddBinContent(bin);
          Int_t entries = entr->GetEntries();
          entr->SetEntries(++entries);
          once=kFALSE;
        }
      }
    }
  }

  return;

}



//_________________________________________________________________
void AliQnCorrectionsManager::FillHistogramsU2n(AliQnCorrectionsConfiguration* QnConf){

  if(QnConf->GetTwistAndRescalingMethod()!=0) return;
  Double_t fillValues[20];
  const Int_t* var = QnConf->GetTwistAndRescalingAxes()->Var();
  const Int_t  dim = QnConf->GetTwistAndRescalingAxes()->Dim();
  for(Int_t iav=0; iav<dim; iav++) fillValues[iav] = fDataContainer[var[iav]];

  Int_t bin=fOutputHistograms[QnConf->GlobalIndex()]->U2nHistogramE()->GetBin(fillValues);

  TClonesArray* dataVectorArray = GetConfDataVectors(QnConf->GlobalIndex());
  //TIter nextEntry(dataVectorArray);
  AliQnCorrectionsDataVector* dataVector=0x0;
  //while((dataVector=static_cast<AliQnCorrectionsDataVector*>(nextEntry()))) {
  //  if(!dataVector) continue;

  THnF* entr = fOutputHistograms[QnConf->GlobalIndex()]->U2nHistogramE();
  for(Int_t idata=0; idata<dataVectorArray->GetEntriesFast(); idata++){
    dataVector = static_cast<AliQnCorrectionsDataVector*>(dataVectorArray->At(idata));
    //if(!dataVector->CheckEventPlaneDetector(QnConf->LocalIndex())) continue;
    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) {
      fOutputHistograms[QnConf->GlobalIndex()]->U2nHistogram(ih,0)->AddBinContent(bin,TMath::Cos(dataVector->Phi()*ih*2));
      fOutputHistograms[QnConf->GlobalIndex()]->U2nHistogram(ih,1)->AddBinContent(bin,TMath::Sin(dataVector->Phi()*ih*2));
    }
    entr->AddBinContent(bin);
    Int_t entries = entr->GetEntries();
    entr->SetEntries(++entries);

  }
}




////_______________________________________________________________________________
//Bool_t AliQnCorrectionsManager::AddDataVector( Int_t detectorId, Double_t phi, Double_t weight, Int_t dataVectorId) {
//  //
//  // Add datavector to TClonesArray of DataVectors and set its bits
//  //
//
//
//
//  //cout<<fDetectorIdMap[detectorId]<<"  "<<detectorId<<endl;
//  //detectorId=fDetectorIdMap[detectorId]-1;
//  Bool_t keepData=kTRUE;
//  Bool_t useForQn;
//
//  ULong64_t EventPlaneDetectorMask=0;
//
//  // check if data passes selection for Qn of interest and set the bit mask for data used by which Qn configuration
//  AliQnCorrectionsConfiguration* QnConf;
//  for(Int_t iconf=0; iconf<fNumberOfQnConfigurationsForDetector[detectorId]; iconf++){
//    QnConf=static_cast<AliQnCorrectionsConfiguration*>(fAliQnCorrectionsConfigurations[detectorId]->At(iconf));
//    if(!QnConf) continue;
//
//    useForQn=kTRUE;
//
//    if(QnConf->Cuts())         useForQn=QnConf->PassCuts(fDataContainer); 
//    if(QnConf->ChannelList())  if(useForQn) useForQn=QnConf->UseChannel(dataVectorId);
//
//    if(useForQn) EventPlaneDetectorMask |= (1<<QnConf->LocalIndex());
//    //if(useForQn) cout<<QnConf->QnConfigurationName()<<endl;
//
//  }
//
//
//
//  // fill data if it passed selection for Qn of interest
//  if(EventPlaneDetectorMask!=0){
//
//    TClonesArray& dataVectors = *(fDataVectors[detectorId]);
//    Int_t dataVectorPosition=dataVectors.GetEntriesFast();
//    AliQnCorrectionsDataVector *dataVector=new(dataVectors[dataVectorPosition]) AliQnCorrectionsDataVector();
//    dataVector->SetPhi(phi);
//    dataVector->SetWeight(weight);
//    dataVector->SetId(dataVectorId);
//    dataVector->SetEventPlaneDetectorMask(EventPlaneDetectorMask);
//    return kTRUE;
//  } 
//  else return kFALSE;
//
//}



////_______________________________________________________________________________
//Bool_t AliQnCorrectionsManager::AddDataVector( Int_t detectorId, Double_t phi, Double_t weight, Int_t dataVectorId) {
//  //
//  // Add datavector to TClonesArray of DataVectors and set its bits
//  //
//
//  TClonesArray& dataVectors = *(fDataVectors[detectorId]);
//  Int_t dataVectorPosition=dataVectors.GetEntriesFast();
//  AliQnCorrectionsDataVector *dataVector=new(dataVectors[dataVectorPosition]) AliQnCorrectionsDataVector();
//  dataVector->SetPhi(phi);
//  dataVector->SetWeight(weight);
//  dataVector->SetId(dataVectorId);
//
//
//  AliQnCorrectionsConfiguration* QnConf;
//
//  Bool_t keepData=kTRUE;
//  Bool_t useForQn;
//
//
//  for(Int_t iconf=0; iconf<fNumberOfQnConfigurationsForDetector[detectorId]; iconf++){
//
//    QnConf=static_cast<AliQnCorrectionsConfiguration*>(fAliQnCorrectionsConfigurations[detectorId]->At(iconf));
//    if(!QnConf) continue;
//
//    useForQn=kTRUE;
//
//    if(QnConf->Cuts())         useForQn=QnConf->PassCuts(fDataContainer); 
//    if(QnConf->ChannelList())  if(useForQn) useForQn=QnConf->UseChannel(dataVectorId);
//
//    if(useForQn) {dataVector->SetEventPlaneDetector(QnConf->LocalIndex());keepData=kTRUE;}
//    //if(useForQn) {cout<<dataVectorPosition<<"  "<<QnConf->LocalIndex()<<endl;}
//
//  }
//  if(!keepData) dataVectors.RemoveAt(dataVectorPosition);
//  return keepData;
//
//}




////_____________________________________________________________________
//void AliQnCorrectionsManager::TwoDetectorCorrelationTwistQvec(Float_t* values, Int_t corpar) {
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
//  for(Int_t idet=0; idet<fNumberOfDetectors; idet++){
//   AliQnCorrectionsConfiguration* QnConf = 0x0;
//   TIter nextQnConf(fAliQnCorrectionsConfigurations[idet]);
//   while((QnConf=static_cast<AliQnCorrectionsConfiguration*>(nextQnConf()))) {
//    if(!QnConf) continue;
//
//    if( QnConf->CalibrationStep()<4) continue;
//    if( QnConf->GetTwistAndRescalingMethod()!=2) continue;
//
//    TProfile * correlationProfiles[3][QnConf->MaximumHarmonic()][4];
//
//    for(Int_t ic=0; ic<3; ++ic) 
//	    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  QnConf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist;
//
//	  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) {
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
//      AliQnCorrectionsAxes* EPbinning =  QnConf->CalibrationBinning();
//          
//      AliQnCorrectionsQnVector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(QnConf->GlobalIndex());
//      TIter nextQvector(qvecList);
//      while((qVector=static_cast<AliQnCorrectionsQnVector*>(nextQvector()))) {
//        if(!qVector) continue;
//
//	      if(!qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kRecentered)) continue;
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
//	      qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kDiagonalized);
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
//void AliQnCorrectionsManager::CorrelationRescalingQvec(Float_t* values, Int_t corpar) {
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
//  for(Int_t idet=0; idet<fNumberOfDetectors; idet++){
//   AliQnCorrectionsConfiguration* QnConf = 0x0;
//   TIter nextQnConf(fAliQnCorrectionsConfigurations[idet]);
//   while((QnConf=static_cast<AliQnCorrectionsConfiguration*>(nextQnConf()))) {
//    if(!QnConf) continue;
//
//
//    if( QnConf->CalibrationStep()<4) continue;
//    if( QnConf->GetTwistAndRescalingMethod()!=1) continue;
//
//    TProfile * correlationProfiles[3][QnConf->MaximumHarmonic()][4];
//
//    for(Int_t ic=0; ic<3; ++ic) 
//	    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  QnConf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxScaled, QyScaled;
//
//	  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) {
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
//      AliQnCorrectionsAxes* EPbinning =  QnConf->CalibrationBinning();
//          
//      AliQnCorrectionsQnVector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(QnConf->GlobalIndex());
//      TIter nextQvector(qvecList);
//      while((qVector=static_cast<AliQnCorrectionsQnVector*>(nextQvector()))) {
//        if(!qVector) continue;
//
//	      if(!qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kDiagonalized)) continue;
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
//	      qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kRescaled);
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
//void AliQnCorrectionsManager::CorrelationTwistAndRescalingQvec(Float_t* values, Int_t corpar) {
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
// for(Int_t idet=0; idet<fNumberOfDetectors; idet++){
//   AliQnCorrectionsConfiguration* QnConf = 0x0;
//   TIter nextQnConf(fAliQnCorrectionsConfigurations[idet]);
//   while((QnConf=static_cast<AliQnCorrectionsConfiguration*>(nextQnConf()))) {
//    if(!QnConf) continue;
//
//    if( QnConf->CalibrationStep()<4) continue;
//    if( QnConf->GetTwistAndRescalingMethod()!=1) continue;
//
//    TProfile * correlationProfiles[3][QnConf->MaximumHarmonic()][4];
//
//    for(Int_t ic=1; ic<3; ++ic) 
//      for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  QnConf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//
//	  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) {
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
//      AliQnCorrectionsAxes* EPbinning =  QnConf->CalibrationBinning();
//          
//      AliQnCorrectionsQnVector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(QnConf->GlobalIndex());
//      TIter nextQvector(qvecList);
//      while((qVector=static_cast<AliQnCorrectionsQnVector*>(nextQvector()))) {
//        if(!qVector) continue;
//
//        if(!qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kRecentered)) continue;
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
//        qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kDiagonalized);
//        qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kRescaled);
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
//void AliQnCorrectionsManager::ThreeDetectorCorrelationTPCTwistQvec(Float_t* values, Int_t corpar) {
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
// for(Int_t idet=0; idet<fNumberOfDetectors; idet++){
//   if(idet==kTPC) continue;
//   AliQnCorrectionsConfiguration* QnConf = 0x0;
//   TIter nextQnConf(fAliQnCorrectionsConfigurations[idet]);
//   while((QnConf=static_cast<AliQnCorrectionsConfiguration*>(nextQnConf()))) {
//    if(!QnConf) continue;
//
//    if( QnConf->CalibrationStep()<4) continue;
//    if(!QnConf->doTwist()) continue;
//    if( QnConf->GetTwistAndRescalingMethod()!=2) continue;
//
//    TProfile * correlationProfiles[3][AliQnCorrectionsConstants::nHarmonics][4];
//
//    for(Int_t ic=0; ic<3; ++ic) 
//      for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  QnConf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//    Double_t x1xt,y1yt,x1yt,x1x2,x2xt,x2yt;
//
//    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) {
//
//      //x1xt = correlationProfiles[0][ih-1][0]->GetBinContent(correlationProfiles[0][ih-1][0]->FindBin(corpar));
//      //y1yt = correlationProfiles[0][ih-1][3]->GetBinContent(correlationProfiles[0][ih-1][3]->FindBin(corpar));
//      //x1yt = correlationProfiles[0][ih-1][1]->GetBinContent(correlationProfiles[0][ih-1][1]->FindBin(corpar));
//      x1xt = QnConf->CorrelationProfile(0, ih, 0)->GetBinContent(QnConf->CorrelationProfile(0, ih, 0)->FindBin(corpar));
//      y1yt = QnConf->CorrelationProfile(0, ih, 3)->GetBinContent(QnConf->CorrelationProfile(0, ih, 3)->FindBin(corpar));
//      x1yt = QnConf->CorrelationProfile(0, ih, 1)->GetBinContent(QnConf->CorrelationProfile(0, ih, 1)->FindBin(corpar));
//
//      Lp = x1yt/x1xt;
//      Ln = x1yt/y1yt;
//      //if(QnConf->QnConfigurationName().EqualTo("VZEROA")) if(ih==2) cout<<corpar<<"  "<<x1yt<<"  "<<x1xt<<"   "<<y1yt<<endl;
//      //if(QnConf->QnConfigurationName().EqualTo("VZEROA")) cout<<x1yt<<"  "<<x1xt<<"   "<<y1yt<<endl;
//      //if(ih==2) cout<<QnConf->QnConfigurationName()<<"  "<<x1yt<<"  "<<x1xt<<"   "<<y1yt<<endl;
//
//      //cout<<QnConf->QnConfigurationName()<<"  "<<Lp<<"  "<<Ln<<"  "<<x1xt<<"  "<<y1yt<<"  "<<x1yt<<"  "<<ih<<endl;
//      if(!(Lp>-99999999.&&Lp<99999999.)) continue;
//      if(!(Ln>-99999999.&&Ln<99999999.)) continue;
//
//      AliQnCorrectionsAxes* EPbinning =  QnConf->CalibrationBinning();
//          
//      AliQnCorrectionsQnVector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(QnConf->GlobalIndex());
//      for(Int_t ibin=0; ibin<qvecList->GetEntriesFast(); ibin++){
//        qVector=static_cast<AliQnCorrectionsQnVector*>(qvecList->At(ibin));
//
//        if(QnConf->QnConfigurationName().Contains("NoRec")){
//          if(!qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kEqualized)) continue;}
//        else if(!qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kRecentered)) continue;
//
//        //if(!qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kRecentered)) continue;
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
//        qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kDiagonalized);
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
//void AliQnCorrectionsManager::ThreeDetectorCorrelationTPCRescalingQvec(Float_t* values, Int_t corpar) {
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
// for(Int_t idet=0; idet<fNumberOfDetectors; idet++){
//   if(idet==kTPC) continue;
//   AliQnCorrectionsConfiguration* QnConf = 0x0;
//   TIter nextQnConf(fAliQnCorrectionsConfigurations[idet]);
//   while((QnConf=static_cast<AliQnCorrectionsConfiguration*>(nextQnConf()))) {
//    if(!QnConf) continue;
//
//    if( QnConf->CalibrationStep()<5) continue;
//    if(!QnConf->doScaling()) continue;
//    if( QnConf->GetTwistAndRescalingMethod()!=2) continue;
//
//    //TProfile * correlationProfiles[3][QnConf->MaximumHarmonic()][4];
//
//    //for(Int_t ic=0; ic<3; ++ic) 
//    //  for(Int_t ih=1; ih<=QnConf->MaximumHarmonic(); ++ih) 
//    //   for(Int_t ip=0; ip<4; ++ip) {
//    //  correlationProfiles[ic][ih-1][ip] =  QnConf->CorrelationProfile(ic, ih, ip);
//    //}
//
//    Double_t Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//    Double_t x1xt,y1yt,x1yt,x1x2,x2xt,x2yt,x1y2,y2yt;
//
//    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) {
//
//
//
//      x1xt = QnConf->CorrelationProfile(0, ih, 0)->GetBinContent(QnConf->CorrelationProfile(0, ih, 0)->FindBin(corpar));
//      y1yt = QnConf->CorrelationProfile(0, ih, 3)->GetBinContent(QnConf->CorrelationProfile(0, ih, 3)->FindBin(corpar));
//      x1yt = QnConf->CorrelationProfile(0, ih, 1)->GetBinContent(QnConf->CorrelationProfile(0, ih, 1)->FindBin(corpar));
//
//      x1x2 = QnConf->CorrelationProfile(2, ih, 0)->GetBinContent(QnConf->CorrelationProfile(2, ih, 0)->FindBin(corpar));
//
//      x2xt = QnConf->CorrelationProfile(1, ih, 0)->GetBinContent(QnConf->CorrelationProfile(1, ih, 0)->FindBin(corpar));
//      x2yt = QnConf->CorrelationProfile(1, ih, 1)->GetBinContent(QnConf->CorrelationProfile(1, ih, 1)->FindBin(corpar));
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
//      //if(QnConf->QnConfigurationName().Contains("VZEROA")&&ih==2) cout<<"  "<<values[corpar]<<"  "<<ih<<"  "<<Ap<<"  "<<An<<endl;
//
//      if(!(Ap>-99999999.&&Ap<99999999.)) continue;
//      if(!(An>-99999999.&&An<99999999.)) continue;
//
//
//
//      AliQnCorrectionsAxes* EPbinning =  QnConf->CalibrationBinning();
//          
//      AliQnCorrectionsQnVector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(QnConf->GlobalIndex());
//      for(Int_t ibin=0; ibin<qvecList->GetEntriesFast(); ibin++){
//        qVector=static_cast<AliQnCorrectionsQnVector*>(qvecList->At(ibin));
//
//        //cout<<"hello  "<<qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kEqualized)<<"  "<<qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kUndefined)<<endl;
//	      if(!qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kDiagonalized)) continue;
//      
//        //cout<<"hey  "<<qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kEqualized)<<"  "<<qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kUndefined)<<endl;
//        Qx = qVector->Qx(ih);
//        Qy = qVector->Qy(ih);
//
//        QxRescaled = Qx / Ap;
//        QyRescaled = Qy / An;
//        
//	      qVector->SetQx( ih, QxRescaled);
//     	  qVector->SetQy( ih, QyRescaled);
//
//        //cout<<"hey 2  "<<qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kRescaled)<<endl;
//	      qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kRescaled);
//        //cout<<"hey 3  "<<qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kRescaled)<<endl;
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
//void AliQnCorrectionsManager::ThreeDetectorCorrelationTPCTwistAndRescalingQvec(Float_t* values, Int_t corpar) {
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
// for(Int_t idet=0; idet<fNumberOfDetectors; idet++){
//   if(idet==kTPC) continue;
//   AliQnCorrectionsConfiguration* QnConf = 0x0;
//   TIter nextQnConf(fAliQnCorrectionsConfigurations[idet]);
//   while((QnConf=static_cast<AliQnCorrectionsConfiguration*>(nextQnConf()))) {
//    if(!QnConf) continue;
//
//    if( QnConf->CalibrationStep()<3) continue;
//    if( QnConf->GetTwistAndRescalingMethod()!=2) continue;
//
//    TProfile * correlationProfiles[3][QnConf->MaximumHarmonic()][4];
//
//    for(Int_t ic=1; ic<3; ++ic) 
//	    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  QnConf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//    Double_t x1xt,y1yt,x1yt,x1x2,x2xt,x2yt;
//
//	  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) {
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
//      AliQnCorrectionsAxes* EPbinning =  QnConf->CalibrationBinning();
//          
//      AliQnCorrectionsQnVector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(QnConf->GlobalIndex());
//      TIter nextQvector(qvecList);
//      while((qVector=static_cast<AliQnCorrectionsQnVector*>(nextQvector()))) {
//        if(!qVector) continue;
//
//	      if(!qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kRecentered)) continue;
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
//	      qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kDiagonalized);
//	      qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kRescaled);
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
void AliQnCorrectionsManager::CallStepTwistAndRescaleQnVector(AliQnCorrectionsConfiguration* QnConf) {

  //
  // Recenter the detector event plane
  //
  
  
  Int_t bin=0;
  Int_t maxHarmonic;
  Double_t fillValues[20];
  const Int_t* var = QnConf->GetTwistAndRescalingAxes()->Var();
  const Int_t  dim = QnConf->GetTwistAndRescalingAxes()->Dim();
  for(Int_t iav=0; iav<dim; iav++) fillValues[iav] = fDataContainer[var[iav]];

  AliQnCorrectionsQnVector* QvectorIn= static_cast<AliQnCorrectionsQnVector*>(CorrectedQnVector(QnConf->GlobalIndex())->At(0));
  TClonesArray& arr = *(fCorrectedQvectors[QnConf->GlobalIndex()][AliQnCorrectionsConstants::kTwist]);
  AliQnCorrectionsQnVector* QvectorTwist   = new(arr[0]) AliQnCorrectionsQnVector(*QvectorIn);
  TClonesArray& arr2 = *(fCorrectedQvectors[QnConf->GlobalIndex()][AliQnCorrectionsConstants::kRescaling]);
  AliQnCorrectionsQnVector* QvectorRescale = new(arr2[0]) AliQnCorrectionsQnVector(*QvectorIn);

  if(QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kTwist))     fLastStep[QnConf->GlobalIndex()]=AliQnCorrectionsConstants::kTwist;
  if(QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kRescaling)) fLastStep[QnConf->GlobalIndex()]=AliQnCorrectionsConstants::kRescaling;

  //Double_t cos2n, sin2n, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaling, QyRescaling, nentries;



  if(QnConf->GetTwistAndRescalingMethod()==0){

  bin=fInputHistograms[QnConf->GlobalIndex()]->CalibrationHistogramE(0)->GetBin(fillValues);

  AliQnCorrectionsSteps::TwistAndRescale2nQn(QvectorIn, QvectorTwist, QvectorRescale, fInputHistograms[QnConf->GlobalIndex()], bin, QnConf->MinimumHarmonic(), QnConf->MaximumHarmonic(),QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kTwist), QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kRescaling));

  }

  if(QnConf->GetTwistAndRescalingMethod()==2){

    Int_t eventClassParameter=0; // for this method we use profiles and only do twist and rescaling corrections against one user defined axis (currently hardcoded as the first axis)

  bin=fInputHistograms[QnConf->GlobalIndex()]->CorrelationProf(0, 0, QnConf->MinimumHarmonic(), 0, eventClassParameter)->FindBin(fDataContainer[var[eventClassParameter]]);

  AliQnCorrectionsSteps::TwistAndRescale3DetectorCorrelation(QvectorIn, QvectorTwist, QvectorRescale, fInputHistograms[QnConf->GlobalIndex()], bin, QnConf->MinimumHarmonic(), QnConf->MaximumHarmonic(),QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kTwist), QnConf->IsApplyCorrection(AliQnCorrectionsConstants::kRescaling),eventClassParameter);

  }


  return;
}



//_____________________________________________________________________
//void AliQnCorrectionsManager::CallStepRescaleQnVector(Int_t u2npar) {

  //
  // Recenter the detector event plane
  //


  // Int_t bin=0;
  // Int_t* var;
  // Int_t maxHarmonic;
  // Int_t dim;


  // Double_t fillValues[20];

  //AliQnCorrectionsConfiguration* QnConf = 0x0;
  //for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
  //  AliQnCorrectionsConfiguration* QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
  //  if(!QnConf) continue;

  //  if( QnConf->CalibrationStep()<4) continue;
  //  //if( !QnConf->doScaling()) continue;
  //  if(QnConf->GetTwistAndRescalingMethod()!=0&&QnConf->GetTwistAndRescalingMethod()!=1) continue;
  //  if( !QnConf->doScaling()) continue;
  //  //if(QnConf->GetTwistAndRescalingMethod()==1) maxHarmonic = QnConf->MaximumHarmonic();
  //  //else maxHarmonic = 2*QnConf->MaximumHarmonic();

  //  //TProfile * U2nProfiles[maxHarmonic][2];

  //  Double_t cos2n, Scos2n, Ncos2n, sin2n, Lp, Ln, Ap, An, Qx, Qy, QxRescaling, QyRescaling;


  //  const Int_t* var = QnConf->CalibrationBinning().Var();
  //  const Int_t  dim = QnConf->CalibrationBinning().Dim();
  //  for(Int_t iv=0; iv<dim; iv++) {fillValues[iv] = fDataContainer[var[iv]];}

  //  for(Int_t ih=1; ih<=(Int_t) QnConf->MaximumHarmonic()/2.; ++ih) {

  //    Scos2n  = fInputHistograms[iconf]->U2nHistogram(ih,1)->GetBinContent(fInputHistograms[iconf]->U2nHistogram(ih,1)->GetBin(fillValues));
  //    Ncos2n  = fInputHistograms[iconf]->U2nHistogramE()->GetBinContent(fInputHistograms[iconf]->U2nHistogramE()->GetBin(fillValues));

  //    if(Ncos2n>1) cos2n=Scos2n/Ncos2n;
  //    else cos2n=10E7;

  //    Ap = 1+cos2n;
  //    An = 1-cos2n;

  //    if(!(Ap>-10E6&&Ap<10E6)) continue;
  //    if(!(An>-10E6&&An<10E6)) continue;



  //    AliQnCorrectionsQnVector* qVector=0x0;
  //    //TClonesArray* qvecList = GetQvectors(QnConf->GlobalIndex());
  //    TClonesArray* qvecList = Qvectors(iconf);
  //    for(Int_t ibin=0; ibin<qvecList->GetEntriesFast(); ibin++){
  //      qVector=static_cast<AliQnCorrectionsQnVector*>(qvecList->At(ibin));
  //      //if(!qVector) break;

  //      if(!qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kRecentered)) continue;

  //      Qx = qVector->Qx(ih);
  //      Qy = qVector->Qy(ih);

  //      QxRescaling = Qx / Ap;
  //      QyRescaling = Qy / An;

  //      qVector->SetQx( ih, QxRescaling);
  //      qVector->SetQy( ih, QyRescaling);

  //      qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kRescaled);
  //    }
  //  }
  //}

  //return;
//}


//
//
//
////_____________________________________________________________________
//void AliQnCorrectionsManager::U2nTwistAndRescalingQvec(Float_t* values, Int_t u2npar) {
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
// for(Int_t idet=0; idet<fNumberOfDetectors; idet++){
//   AliQnCorrectionsConfiguration* QnConf = 0x0;
//   TIter nextQnConf(fAliQnCorrectionsConfigurations[idet]);
//   while((QnConf=static_cast<AliQnCorrectionsConfiguration*>(nextQnConf()))) {
//    if(!QnConf) continue;
//
//
//    if( QnConf->CalibrationStep()<4) continue;
//    if( QnConf->GetTwistAndRescalingMethod()!=0) continue;
//
//    TProfile * U2nProfiles[AliQnCorrectionsConstants::nHarmonics*2][2];
//
//	    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic()*2; ++ih) 
//       for(Int_t ip=0; ip<2; ++ip) {
//      U2nProfiles[ih-1][ip] =  QnConf->U2nProfile(ih, ip);
//    }
//
//    Double_t cos2n, sin2n, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//
//	  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) {
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
//      AliQnCorrectionsAxes* EPbinning =  QnConf->CalibrationBinning();
//          
//      AliQnCorrectionsQnVector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(QnConf->GlobalIndex());
//      TIter nextQvector(qvecList);
//      while((qVector=static_cast<AliQnCorrectionsQnVector*>(nextQvector()))) {
//        if(!qVector) continue;
//
//	      if(!qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kRecentered)) continue;
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
//	      qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kDiagonalized);
//	      qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kRescaled);
//      }
//    }
//  }
// }
//
//  return;
//}
//
//

////__________________________________________________________________
//void AliQnCorrectionsManager::InitializeCalibrationHistograms(){
//
//  AliQnCorrectionsConfiguration* QnConf = 0x0;
//  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
//    AliQnCorrectionsConfiguration* QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
//    if(!QnConf) continue;
//
//  }
//
//}

//__________________________________________________________________
void AliQnCorrectionsManager::DisableFillHistograms(AliQnCorrectionsConstants::CorrectionSteps step){

    for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    AliQnCorrectionsConfiguration* QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);

      for(Int_t istep=(Int_t) step; istep<AliQnCorrectionsConstants::kNcorrectionSteps; istep++){
        QnConf->SetFillHistogram(istep,kFALSE);
      }

    }

}

//__________________________________________________________________
void AliQnCorrectionsManager::DisableCorrections(AliQnCorrectionsConstants::CorrectionSteps step){

    for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    AliQnCorrectionsConfiguration* QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);

      QnConf->SetApplyCorrection(step,kFALSE);

    }

}


//__________________________________________________________________
void AliQnCorrectionsManager::InitializeCalibrationHistograms(){


  AliQnCorrectionsConfiguration* QnConf = 0x0;
  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
    //if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)) if(fCorrectionStep<AliQnCorrectionsConstants::kDataVectorEqualization) fCorrectionStep=(Int_t) AliQnCorrectionsConstants::kDataVectorEqualization;
    //if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRecentering))            if(fCorrectionStep<AliQnCorrectionsConstants::kRecentering           ) fCorrectionStep=(Int_t) AliQnCorrectionsConstants::kRecentering;
    //if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kAlignment))              if(fCorrectionStep<AliQnCorrectionsConstants::kAlignment             ) fCorrectionStep=(Int_t) AliQnCorrectionsConstants::kAlignment;
    //if(QnConf->GetTwistAndRescalingMethod()==2){                                                                                                   
    //  if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist))                if(fCorrectionStep<AliQnCorrectionsConstants::kTwist                 ) fCorrectionStep=(Int_t) AliQnCorrectionsConstants::kTwist;
    //  if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRescaling))            if(fCorrectionStep<AliQnCorrectionsConstants::kRescaling             ) fCorrectionStep=(Int_t) AliQnCorrectionsConstants::kRescaling;
    //}
    //fCorrectionStep=0;
    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)) if(fCorrectionStep<AliQnCorrectionsConstants::kDataVectorEqualization) {fCorrectionStep++;}
    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRecentering))            if(fCorrectionStep<AliQnCorrectionsConstants::kRecentering           ) {fCorrectionStep++;}
    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kAlignment))              if(fCorrectionStep<AliQnCorrectionsConstants::kAlignment             ) {fCorrectionStep++;}
    if(QnConf->GetTwistAndRescalingMethod()==2){                                                                                                   
      if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist))                if(fCorrectionStep<AliQnCorrectionsConstants::kTwist                 ) fCorrectionStep++;
      if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRescaling))            if(fCorrectionStep<AliQnCorrectionsConstants::kTwist                 ) fCorrectionStep++;
    }
  }

  fPassesRequired=fCorrectionStep;


  //TString label="allData";
  //if(QnConf->CorrectWithEventLabel()) label=fLabel;

  THashList* list = GetInputListWithLabel(fLabel);
  if(!list)  list = GetInputListWithLabel("allData");

  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);

    //cout<<"LIST "<<fListInputHistogramsQnCorrections->GetEntries()<<endl;

    //if(!list) return;


    Int_t step;

    step=(Int_t)AliQnCorrectionsConstants::kDataVectorEqualization;
    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)){
      if(!fInputHistograms[iconf]->ConnectMultiplicityHistograms(list,QnConf)){
        DisableCorrections(AliQnCorrectionsConstants::kDataVectorEqualization);
        fCorrectionStep=step-1;
        if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRecentering)) DisableFillHistograms(AliQnCorrectionsConstants::kRecentering);
      }
    };

    step=(Int_t)AliQnCorrectionsConstants::kRecentering;
    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRecentering)){
      if(!fInputHistograms[iconf]->ConnectMeanQnCalibrationHistograms(list,QnConf)){
         DisableCorrections(AliQnCorrectionsConstants::kRecentering);
         DisableFillHistograms(AliQnCorrectionsConstants::kAlignment);
         DisableFillHistograms(AliQnCorrectionsConstants::kTwist);
         DisableFillHistograms(AliQnCorrectionsConstants::kRescaling);
         if(fCorrectionStep>=step) fCorrectionStep=step-1;
      }
    };

    step=(Int_t)AliQnCorrectionsConstants::kTwist;
    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist)){
      if(QnConf->GetTwistAndRescalingMethod()==0||QnConf->GetTwistAndRescalingMethod()==1){
        //if(!fInputHistograms[iconf]->ConnectU2nQnCalibrationHistograms(list,QnConf)){
        if(!fInputHistograms[iconf]->ConnectMeanQnCalibrationHistograms(list,QnConf)){
          DisableCorrections(AliQnCorrectionsConstants::kTwist);
          DisableCorrections(AliQnCorrectionsConstants::kRescaling);
          if(fCorrectionStep>=step) fCorrectionStep=step-1;
        }
      }
      if(QnConf->GetTwistAndRescalingMethod()==2){
        if(!fInputHistograms[iconf]->ConnectCorrelationQnCalibrationHistograms(list,QnConf)){
          QnConf->SetApplyCorrection(AliQnCorrectionsConstants::kTwist,kFALSE);
          QnConf->SetApplyCorrection(AliQnCorrectionsConstants::kRescaling,kFALSE);
          if(fCorrectionStep>=step) fCorrectionStep=step-1;
      }}
    };

    
    step=(Int_t)AliQnCorrectionsConstants::kAlignment;
    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kAlignment)){
        if(!fInputHistograms[iconf]->ConnectRotationQnCalibrationHistograms(list,QnConf)){
          DisableCorrections(AliQnCorrectionsConstants::kAlignment);
          if(fCorrectionStep>=step) fCorrectionStep=step-1;
          if(QnConf->GetTwistAndRescalingMethod()==2){
            QnConf->SetFillHistogram(AliQnCorrectionsConstants::kTwist,kFALSE);
            QnConf->SetFillHistogram(AliQnCorrectionsConstants::kRescaling,kFALSE);
          }
        }
    };
  }


}

  

//__________________________________________________________________
void AliQnCorrectionsManager::PrintFrameworkInformation(){



	const int widthEntry = 18;
	const int widthBar = 3;

  Bool_t forCorrection=kTRUE;
  PrintFrameworkInformationLine(!forCorrection, "----------------", AliQnCorrectionsConstants::kNothing, widthEntry, widthBar);
  PrintFrameworkInformationLine(forCorrection, Form("FLOW VECTOR FRAMEWORK - PASS %d/%d", fCorrectionStep+1, fPassesRequired+1), AliQnCorrectionsConstants::kNothing, widthEntry, widthBar);
  PrintFrameworkInformationLine(!forCorrection, " ", AliQnCorrectionsConstants::kNothing, widthEntry, widthBar);
  PrintFrameworkInformationLine(forCorrection, "--FLOW VECTORS--", AliQnCorrectionsConstants::kNothing, widthEntry, widthBar);
    
  cout<<setw(widthEntry)<<" "<<setw(widthBar)<<"|";
  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    AliQnCorrectionsConfiguration* QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
    cout<<setw(widthEntry)<<QnConf->QnConfigurationName()<<setw(widthBar)<<"  |";
  }
  cout<<endl;

  PrintFrameworkInformationLine(!forCorrection, "CORRECTIONS", AliQnCorrectionsConstants::kNothing, widthEntry, widthBar);
  PrintFrameworkInformationLine(forCorrection, "Equalization", AliQnCorrectionsConstants::kDataVectorEqualization, widthEntry, widthBar);
  PrintFrameworkInformationLine(forCorrection, "Recentering", AliQnCorrectionsConstants::kRecentering, widthEntry, widthBar);
  PrintFrameworkInformationLine(forCorrection, "Alignment", AliQnCorrectionsConstants::kAlignment, widthEntry, widthBar);
  PrintFrameworkInformationLine(forCorrection, "Twist", AliQnCorrectionsConstants::kTwist, widthEntry, widthBar);
  PrintFrameworkInformationLine(forCorrection, "Rescaling", AliQnCorrectionsConstants::kRescaling, widthEntry, widthBar);
  PrintFrameworkInformationLine(!forCorrection, "FILL HISTS", AliQnCorrectionsConstants::kNothing, widthEntry, widthBar);
  PrintFrameworkInformationLine(!forCorrection, "Equalization", AliQnCorrectionsConstants::kDataVectorEqualization, widthEntry, widthBar);
  PrintFrameworkInformationLine(!forCorrection, "Recentering", AliQnCorrectionsConstants::kRecentering, widthEntry, widthBar);
  PrintFrameworkInformationLine(!forCorrection, "Alignment", AliQnCorrectionsConstants::kAlignment, widthEntry, widthBar);
  PrintFrameworkInformationLine(!forCorrection, "Twist", AliQnCorrectionsConstants::kTwist, widthEntry, widthBar);
  PrintFrameworkInformationLine(!forCorrection, "Rescaling", AliQnCorrectionsConstants::kRescaling, widthEntry, widthBar);
  PrintFrameworkInformationLine(!forCorrection, "----------------", AliQnCorrectionsConstants::kNothing, widthEntry, widthBar);
      
  cout<<setw(widthEntry)<<"Legend "<<setw(widthBar)<<"|";
  cout<<setw(widthEntry)<<"x: this pass"<<setw(widthBar)<<" ";
  cout<<setw(widthEntry)<<"0: future pass"<<setw(widthBar)<<" ";
  cout<<setw(widthEntry)<<"-: N/A"<<setw(widthBar)<<"|"<<endl;

  PrintFrameworkInformationLine(!forCorrection, "----------------", AliQnCorrectionsConstants::kNothing, widthEntry, widthBar);


}

//__________________________________________________________________
void AliQnCorrectionsManager::PrintFrameworkInformationLine(Bool_t HistOrCor,TString correctionname, AliQnCorrectionsConstants::CorrectionSteps stepflag, Int_t widthEntry, Int_t widthBar){





  if(HistOrCor&&stepflag!=AliQnCorrectionsConstants::kNothing){
    cout<<setw(widthEntry)<<correctionname<<setw(widthBar)<<"|";
    for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
      AliQnCorrectionsConfiguration* QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
      if(QnConf->IsApplyCorrection(stepflag)) cout<<setw(widthEntry)<<"x"<<setw(widthBar)<<"|";
      else if(QnConf->IsRequestedCorrection(stepflag)) cout<<setw(widthEntry)<<"0"<<setw(widthBar)<<"|";
      else cout<<setw(widthEntry)<<"-"<<setw(widthBar)<<"|";
    }
  cout<<endl;
  }


  
  if(!HistOrCor&&stepflag!=AliQnCorrectionsConstants::kNothing){
    cout<<setw(widthEntry)<<correctionname<<setw(widthBar)<<"|";
    for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
      AliQnCorrectionsConfiguration* QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
      if(QnConf->IsFillHistogram(stepflag)) cout<<setw(widthEntry)<<"x"<<setw(widthBar)<<"|";
      else if(QnConf->IsRequestedFillHistogram(stepflag)) cout<<setw(widthEntry)<<"0"<<setw(widthBar)<<"|";
      else cout<<setw(widthEntry)<<"-"<<setw(widthBar)<<"|";
    }
  cout<<endl;
  }

  
  
  if(HistOrCor&&stepflag==AliQnCorrectionsConstants::kNothing){
    cout<<setw(widthEntry)<<" "<<setw(widthBar)<<"|";
    for(Int_t iconf=0; iconf<fNumberOfQnConfigurations/2; iconf++){
      cout<<setw(widthEntry)<<" ";
    }
    cout<<correctionname;
    for(Int_t iconf=fNumberOfQnConfigurations/2; iconf<fNumberOfQnConfigurations; iconf++){
      cout<<setw(widthEntry)<<" ";
    }
    cout<<endl;
  }


    
  if(!HistOrCor&&stepflag==AliQnCorrectionsConstants::kNothing){
    cout<<setw(widthEntry)<<correctionname<<setw(widthBar)<<"-";
    for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
        for(Int_t i=0; i<widthEntry; i++) cout<<"-";
        for(Int_t i=0; i<widthBar; i++) cout<<"-";
    }
  cout<<endl;
  }




}



//__________________________________________________________________
void AliQnCorrectionsManager::WriteCalibrationHistogramsToList()
{
  //
  // Finish Task 
  //
  THashList* listall = new THashList();
  TObjArray* runarr = fRunList.Tokenize(";");
  Int_t nruns = runarr->GetEntries();
  if(nruns==0){
    if(fLabel.EqualTo("")) {
      fRunList="allData;";
      fLabel="allData";
    } 
    else {
      fRunList=fLabel+";";
    }
    runarr=fRunList.Tokenize(";");
  };
  for(Int_t ir=0; ir<runarr->GetEntries(); ir++){
    TString runlabel = runarr->At(ir)->GetName();

    AliQnCorrectionsHistograms* histos;

  THashList* dataByRun = new THashList();
  //THashList* dataByRun = fListOutputHistogramsQnCorrections;//new THashList();
  dataByRun->SetName(runlabel);
  dataByRun->SetOwner(kTRUE);

  AliQnCorrectionsConfiguration* QnConf = 0x0;
  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
    if(!QnConf) continue;

    if(runlabel==fLabel) histos = fOutputHistograms[iconf];
    else {histos = new AliQnCorrectionsHistograms();histos->CreateCalibrationHistograms(QnConf);}

    THashList* detector = new THashList();
    THashList* detectorM = new THashList();
    THashList* detectorC = new THashList();
    detector->SetOwner();
    detectorM->SetOwner();
    detectorC->SetOwner();

    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRecentering)||(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist)&&QnConf->GetTwistAndRescalingMethod())==0){
      Int_t maxHar = QnConf->MaximumHarmonic();
      if(QnConf->GetTwistAndRescalingMethod()==0) maxHar*=2;
        for(Int_t ih=QnConf->MinimumHarmonic(); ih<=maxHar; ++ih) {
          for(Int_t ic=0; ic<2; ++ic){
            detector->Add(histos->CalibrationHistogramQ( (QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization) ? 1 : 0) ,ih,ic));
          }
        }
        detector->Add(histos->CalibrationHistogramE((QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization) ? 1 : 0) ));
    }

    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist)&&QnConf->GetTwistAndRescalingMethod()==1){
    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){
      if(histos->U2nHistogram(ih,0)){
        detector->Add(histos->U2nHistogram(ih,0));
        detector->Add(histos->U2nHistogram(ih,1));
      }
    }
    detector->Add(histos->U2nHistogramE());
    }

    detector->SetName(Form("Qvec%s",QnConf->QnConfigurationName().Data()));
    if(detector->GetEntries()!=0){
      dataByRun->Add(detector);
    }


    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)){
      //for(Int_t idim=0; idim<10; idim++) if(histos->EventHistogram(idim)) detectorM->Add(histos->EventHistogram(idim));
      for(Int_t is=0; is<1; ++is){
        if(histos->EqualizationHistogramM(is)) detectorM->Add(histos->EqualizationHistogramM(is));
        if(histos->EqualizationHistogramE(is)) detectorM->Add(histos->EqualizationHistogramE(is));
        if(histos->GroupEqualizationHistogramM(is)) detectorM->Add(histos->GroupEqualizationHistogramM(is));
        if(histos->GroupEqualizationHistogramE(is)) detectorM->Add(histos->GroupEqualizationHistogramE(is));
      }
      detectorM->SetName(Form("Mult%s",QnConf->QnConfigurationName().Data()));
      if(detectorM->GetEntries()!=0){
        dataByRun->Add(detectorM);
      }
    }

  //  

    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kAlignment)){
      for(Int_t is=0; is<4; ++is) if(histos->GetRotationHistogram(0,is)) detectorC->Add(histos->GetRotationHistogram(0,is));
      if(histos->GetRotationHistogramE(0)) detectorC->Add(histos->GetRotationHistogramE(0));
    }




    if((QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist)||QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRescaling))&&QnConf->GetTwistAndRescalingMethod()==2){
    for(Int_t icomb=0; icomb<3; ++icomb){ 
      for(Int_t icomp=0; icomp<4; ++icomp){ 
        for(Int_t iaxis=0; iaxis<=QnConf->GetAlignmentAxes()->Dim(); ++iaxis){ 
          for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){ 
            if(histos->CorrelationProf(0,icomb,ih,icomp,iaxis)) detectorC->Add(histos->CorrelationProf(0,icomb,ih,icomp,iaxis));
      }}}}
      detectorC->SetName(Form("Correlations%s",QnConf->QnConfigurationName().Data()));
    }
    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kAlignment)||((QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist)||QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRescaling))&&QnConf->GetTwistAndRescalingMethod()==2)){
      detectorC->SetName(Form("Correlations%s",QnConf->QnConfigurationName().Data()));
      dataByRun->Add(detectorC);
    }
  }

  
  if(runlabel==fLabel){
    listall->SetName("allData");
    for(Int_t i=0; i<dataByRun->GetEntries(); i++){
      THashList* l = (THashList*) dataByRun->At(i);
      l->SetOwner();
      listall->Add(l->Clone());
    }
  }

  if(!(runlabel.EqualTo("allData"))) fListOutputHistogramsQnCorrections->Add(dataByRun->Clone());
  
  }

  fListOutputHistogramsQnCorrections->Add(listall);
  //fListOutputHistogramsQnCorrections->Add(dataAllRuns);

}

//__________________________________________________________________
//void AliQnCorrectionsManager::WriteQaHistogramsToList()
//{
//  //
//  // Finish Task 
//  //
//  THashList* listall = new THashList();
//  TObjArray* runarr = fRunList.Tokenize(";");
//  for(Int_t ir=0; ir<runarr->GetEntries(); ir++){
//    TString runlabel = runarr->At(ir)->GetName();
//
//    AliQnCorrectionsHistograms* histos;
//    histos = new AliQnCorrectionsHistograms();
//
//  THashList* dataByRun = new THashList();
//  //THashList* dataByRun = fListOutputHistogramsQnCorrections;//new THashList();
//  dataByRun->SetName(runlabel);
//  dataByRun->SetOwner(kTRUE);
//
//  AliQnCorrectionsConfiguration* QnConf = 0x0;
//  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
//    QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
//    if(!QnConf) continue;
//
//    if(runlabel==fLabel) histos = fOutputHistograms[iconf];
//    else histos->CreateCalibrationHistograms(QnConf);
//
//    THashList* detector = new THashList();
//    THashList* detectorM = new THashList();
//    THashList* detectorC = new THashList();
//    detector->SetOwner();
//    detectorM->SetOwner();
//    detectorC->SetOwner();
//
//    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRecentering)||(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist)&&QnConf->GetTwistAndRescalingMethod())==0){
//      Int_t maxHar = QnConf->MaximumHarmonic();
//      if(QnConf->GetTwistAndRescalingMethod()==0) maxHar*=2;
//        for(Int_t ih=QnConf->MinimumHarmonic(); ih<=maxHar; ++ih) {
//          for(Int_t ic=0; ic<2; ++ic){
//            detector->Add(histos->CalibrationHistogramQ( (QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization) ? 1 : 0) ,ih,ic));
//          }
//        }
//        detector->Add(histos->CalibrationHistogramE((QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization) ? 1 : 0) ));
//    }
//
//    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist)&&QnConf->GetTwistAndRescalingMethod()==1){
//    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){
//      if(histos->U2nHistogram(ih,0)){
//        detector->Add(histos->U2nHistogram(ih,0));
//        detector->Add(histos->U2nHistogram(ih,1));
//      }
//    }
//    detector->Add(histos->U2nHistogramE());
//    }
//
//    detector->SetName(Form("Qvec%s",QnConf->QnConfigurationName().Data()));
//    if(detector->GetEntries()!=0){
//      dataByRun->Add(detector);
//    }
//
//
//    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)){
//      //for(Int_t idim=0; idim<10; idim++) if(histos->EventHistogram(idim)) detectorM->Add(histos->EventHistogram(idim));
//      for(Int_t is=0; is<1; ++is){
//        if(histos->EqualizationHistogramM(is)) detectorM->Add(histos->EqualizationHistogramM(is));
//        if(histos->EqualizationHistogramE(is)) detectorM->Add(histos->EqualizationHistogramE(is));
//        if(histos->GroupEqualizationHistogramM(is)) detectorM->Add(histos->GroupEqualizationHistogramM(is));
//        if(histos->GroupEqualizationHistogramE(is)) detectorM->Add(histos->GroupEqualizationHistogramE(is));
//      }
//      detectorM->SetName(Form("Mult%s",QnConf->QnConfigurationName().Data()));
//      if(detectorM->GetEntries()!=0){
//        dataByRun->Add(detectorM);
//      }
//    }
//
//  //  
//
//    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kAlignment)){
//      for(Int_t is=0; is<4; ++is) if(histos->GetRotationHistogram(0,is)) detectorC->Add(histos->GetRotationHistogram(0,is));
//      if(histos->GetRotationHistogramE(0)) detectorC->Add(histos->GetRotationHistogramE(0));
//    }
//
//
//
//
//    if((QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist)||QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRescaling))&&QnConf->GetTwistAndRescalingMethod()==2){
//    for(Int_t icomb=0; icomb<3; ++icomb){ 
//      for(Int_t icomp=0; icomp<4; ++icomp){ 
//        for(Int_t iaxis=0; iaxis<=QnConf->GetAlignmentAxes()->Dim(); ++iaxis){ 
//          for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){ 
//            if(histos->CorrelationProf(0,icomb,ih,icomp,iaxis)) detectorC->Add(histos->CorrelationProf(0,icomb,ih,icomp,iaxis));
//      }}}}
//      detectorC->SetName(Form("Correlations%s",QnConf->QnConfigurationName().Data()));
//    }
//    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kAlignment)||((QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist)||QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRescaling))&&QnConf->GetTwistAndRescalingMethod()==2)){
//      detectorC->SetName(Form("Correlations%s",QnConf->QnConfigurationName().Data()));
//      dataByRun->Add(detectorC);
//    }
//  }
//  
//  if(!(fLabel.EqualTo("allData"))&&runlabel==fLabel){
//    listall->SetName("allData");
//    for(Int_t i=0; i<dataByRun->GetEntries(); i++){
//      THashList* l = (THashList*) dataByRun->At(i);
//      l->SetOwner();
//      listall->Add(l->Clone());
//    }
//  }
//
//  fListHistogramsQA->Add(dataByRun->Clone());
//  
//  }
//
//  fListHistogramsQA->Add(listall);
//  //fListOutputHistogramsQnCorrections->Add(dataAllRuns);
//
//}






//__________________________________________________________________
void AliQnCorrectionsManager::WriteQaHistogramsToList()
{
  //
  // Finish Task 
  //

  THashList* listall = new THashList();


  //detector->SetOwner();
  //detectorMqa->SetOwner();
  //detectorC->SetOwner();
  //dataByRunQA->SetOwner();


  AliQnCorrectionsHistograms* histos;

  THashList* dataByRunQA = new THashList();
  dataByRunQA->SetOwner();
    //dataByRunQA->Clear();


  dataByRunQA->SetName(fLabel);


  AliQnCorrectionsConfiguration* QnConf = 0x0;
  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
    if(!QnConf) continue;
    THashList* detector    = new THashList();
    THashList* detectorMqa = new THashList();
    THashList* detectorC   = new THashList();
    detector->SetOwner();
    detectorMqa->SetOwner();
    detectorC->SetOwner();
    //detector->Clear();
    //detectorMqa->Clear();
    //detectorC->Clear();

    histos = fOutputHistograms[iconf];

     Int_t istep=0;
    for(Int_t istep=0; istep<AliQnCorrectionsConstants::nCorrectionSteps; istep++){
      if(istep==0) if(!QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)) continue;
      if(istep==1) continue;
      if(istep>1)  if(!QnConf->IsRequestedCorrection(istep)) continue;
      

    Int_t maxHar = QnConf->MaximumHarmonic();
      for(Int_t ih=QnConf->MinimumHarmonic(); ih<=maxHar; ++ih) {
        for(Int_t ic=0; ic<2; ++ic){
          detector->Add(histos->CalibrationHistogramQ(istep,ih,ic));
        }
      }
      detector->Add(histos->CalibrationHistogramE(istep));
    }

    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){
      if(histos->U2nHistogram(ih,0)){
        detector->Add(histos->U2nHistogram(ih,0));
        detector->Add(histos->U2nHistogram(ih,1));
      }
    }
    if(histos->U2nHistogramE()) detector->Add(histos->U2nHistogramE());

    detector->SetName(Form("Qvec%s",QnConf->QnConfigurationName().Data()));
    dataByRunQA->Add(detector);


    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kDataVectorEqualization)){
      for(Int_t is=1; is<3; ++is){
        if(histos->EqualizationHistogramM(is)) detectorMqa->Add(histos->EqualizationHistogramM(is));
        if(histos->EqualizationHistogramE(is)) detectorMqa->Add(histos->EqualizationHistogramE(is));
      }
      detectorMqa->SetName(Form("Mult%s",QnConf->QnConfigurationName().Data()));
      dataByRunQA->Add(detectorMqa);
    }

    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kAlignment)){
      for(Int_t is=0; is<4; ++is) if(histos->GetRotationHistogram(1,is)) detectorC->Add(histos->GetRotationHistogram(1,is));
      if(histos->GetRotationHistogramE(1)) detectorC->Add(histos->GetRotationHistogramE(1));
    }

    
    if((QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist)||QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRescaling))&&QnConf->GetTwistAndRescalingMethod()==2){
    for(Int_t icomb=0; icomb<3; ++icomb){ 
      for(Int_t icomp=0; icomp<4; ++icomp){ 
        for(Int_t iaxis=0; iaxis<=QnConf->GetAlignmentAxes()->Dim(); ++iaxis){ 
          for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){ 
            if(histos->CorrelationProf(1,icomb,ih,icomp,iaxis)) detectorC->Add(histos->CorrelationProf(1,icomb,ih,icomp,iaxis));
            if(histos->CorrelationProf(2,icomb,ih,icomp,iaxis)) detectorC->Add(histos->CorrelationProf(2,icomb,ih,icomp,iaxis));
      }}}}
    }

    if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kAlignment)||((QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist)||QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRescaling))&&QnConf->GetTwistAndRescalingMethod()==2)){
      detectorC->SetName(Form("Correlations%s",QnConf->QnConfigurationName().Data()));
      dataByRunQA->Add(detectorC);
    }


    //for(Int_t istep=0; istep<=AliQnCorrectionsConstants::kNcorrectionSteps; istep++){
    //  for(Int_t icomb=0; icomb<3; ++icomb){ 
    //    for(Int_t icomp=0; icomp<4; ++icomp){ 
    //      for(Int_t iaxis=0; iaxis<=QnConf->CalibrationBinning()->Dim(); ++iaxis){ 
    //      for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){ 
    //        //if(histos->CorrelationProf(QnConf->GetCorrectionStep(istep),icomb,ih,icomp,iaxis)&&QnConf->GetCorrectionStep(istep)>2) detectorCqa->Add(histos->CorrelationProf(QnConf->GetCorrectionStep(istep),icomb,ih,icomp,iaxis));
    //        //if(histos->CorrelationEpProf(QnConf->GetCorrectionStep(istep),icomb,ih,icomp,iaxis)) detectorCqa->Add(histos->CorrelationEpProf(QnConf->GetCorrectionStep(istep),icomb,ih,icomp,iaxis));
    //      }}}}}
    //detectorCqa->SetName(Form("Correlations%s",QnConf->QnConfigurationName().Data()));
    //dataByRunQA->Add(detectorCqa);
    //dataAllRunsQA->Add(detectorCqa);

  }
  //if(!(fLabel.EqualTo("allData"))){
  //  listall = (THashList*) dataByRunQA->Clone("allData");
  //  //listall->SetName("allData");
  //  //for(Int_t i=0; i<dataByRunQA->GetEntries(); i++){
  //  //  THashList* l = (THashList*) dataByRunQA->At(i);
  //  //  l->SetOwner();
  //  //  listall->Add(l->Clone());
  //  //}
  //}

  fListHistogramsQA->Add(dataByRunQA->Clone());


  //fListHistogramsQA->Add(listall);


}

//_______________________________________________________________________________
void AliQnCorrectionsManager::AddQnConfiguration(AliQnCorrectionsConfiguration* QnConf, Int_t type)
{
    
  AddDetectorType(type);
  Int_t detId=fDetectorIdMap[type]-1;
  
  QnConf->SetDetectorId( (UShort_t) detId);
  TClonesArray& eparr = *(fAliQnCorrectionsConfigurations[ detId]);
  eparr[fNumberOfQnConfigurationsForDetector[ detId]]=QnConf;
  fIndex[fNumberOfQnConfigurations][0] = detId;
  fIndex[fNumberOfQnConfigurations][1] = fNumberOfQnConfigurationsForDetector[ detId];
  QnConf->SetLocalIndex(fNumberOfQnConfigurationsForDetector[detId]);
  QnConf->SetGlobalIndex(fNumberOfQnConfigurations);
  if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kAlignment)) if(QnConf->MinimumHarmonic()>QnConf->AlignmentHarmonic()) QnConf->SetMinimumHarmonic(QnConf->AlignmentHarmonic());
  //QnConf->SetPassNumbers();
  fNumberOfQnConfigurationsForDetector[detId]++;
  fNumberOfQnConfigurations++;
}


//_______________________________________________________________________________
AliQnCorrectionsConfiguration* AliQnCorrectionsManager::GetQnConfiguration(TString name)
{

  AliQnCorrectionsConfiguration* QnConf = 0x0;
  for(Int_t iconf=0; iconf<fNumberOfQnConfigurations; iconf++){
    QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
    if(!QnConf) continue;
    if(name.EqualTo(QnConf->QnConfigurationName())) return QnConf;
  }
  return 0x0;
}

//_______________________________________________________________________________
void AliQnCorrectionsManager::Finalize()
{
  if(fSetFillHistogramsQnCorrections)  WriteCalibrationHistogramsToList();
  if(fSetFillHistogramsQA)             WriteQaHistogramsToList();
  WriteOutputToFile();
}

//_______________________________________________________________________________
void AliQnCorrectionsManager::WriteOutputToFile()
{
  if(fOutputHistogramsQnCorrectionsFile&&fSetFillHistogramsQnCorrections){
    fOutputHistogramsQnCorrectionsFile->cd();
    fListOutputHistogramsQnCorrections->Write(fListOutputHistogramsQnCorrections->GetName(),TObject::kSingleKey);
    //fListOutputHistogramsQnCorrections->Write(fLabel,TObject::kSingleKey);
    //if(!(fLabel.EqualTo("allData"))) fListOutputHistogramsQnCorrections->Write("allData",TObject::kSingleKey);
  }
  if(fHistogramsQAFile&&fSetFillHistogramsQA){
    fHistogramsQAFile->cd();
    fListHistogramsQA->Write(fListHistogramsQA->GetName(),TObject::kSingleKey);
  }
  if(fTreeQnVectorsFile&&fSetFillTreeQnVectors){
    fTreeQnVectorsFile->cd();
    fTreeQnVectors->Write();
  }
}
