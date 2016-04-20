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
 

#include "AliQnCorrectionsSteps.h"
#include "AliQnCorrectionsConfiguration.h"
#include "AliQnCorrectionsDataVector.h"
#include "AliQnCorrectionsHistograms.h"
#include <TMath.h>
#include <TList.h>
#include <THashList.h>
#include <TClonesArray.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TArrayS.h>

ClassImp(AliQnCorrectionsSteps)


//_______________________________________________________________________________
AliQnCorrectionsSteps::AliQnCorrectionsSteps() 
{
  //
  // default constructor
  //



}


//____________________________________________________________________________
AliQnCorrectionsSteps::~AliQnCorrectionsSteps()
{
  //
  // De-Constructor
  //
}



////_________________________________________________________________
//void AliQnCorrectionsSteps::BuildQnVectors(AliQnCorrectionsQnVector* QvectorOut, TClonesArray* dataVectorArray, Int_t QnConfIndex, Int_t minHar, Int_t maxHar, Int_t EqualizationMethod){
//  //
//  // Construct the event plane for a specified detector
//  //
//
//  if(EqualizationMethod==-1) AliQnCorrectionsDataVector::FillQvector(dataVectorArray,QnConfIndex, QvectorOut);
//  else                       AliQnCorrectionsDataVector::FillQvector(dataVectorArray,QnConfIndex, QvectorOut, EqualizationMethod);
//
//
//  if(QvectorOut->N()==0) {               // If detector is empty
//    for(Int_t ih=minHar; ih<=maxHar; ++ih) QvectorOut->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kUndefined);
//    return;
//  }
//  else{
//    for(Int_t ih=minHar; ih<=maxHar; ++ih){
//      if(EqualizationMethod==-1) QvectorOut->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kRaw);
//      else                       QvectorOut->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kEqualized);
//    }
//  }
//
//  return;
//
//}


//_____________________________________________________________________
void AliQnCorrectionsSteps::CalibrateDataVector(TClonesArray* dataVectorArray, AliQnCorrectionsConfiguration* QnConf,  AliQnCorrectionsHistograms* inputHistos, Double_t* fillValues) {
  //
  // Calibrate channel multiplicities
  //


  TString det;
  Double_t average, width, mult, chanAveMult, chanWidthMult;
  //Int_t* var;
  //Int_t dim; 
  //Double_t fillValues[20];
  Double_t groupWeight=1.0;
  Int_t bin,binGroup;

  //const Int_t* var = QnConf->EqualizationBinning()->Var();
  const Int_t  dim = QnConf->EqualizationBinning()->Dim();

  
  //for(Int_t iav=0; iav<(dim-1); iav++) fillValues[iav] = fDataContainer[var[iav]];

  AliQnCorrectionsDataVector* dataVector= 0x0;
  //TIter nextEntry(dataVectorArray);

  //while((dataVector=static_cast<AliQnCorrectionsDataVector*>(nextEntry()))) {
  //  if(!dataVector) continue;

  for(Int_t idata=0; idata<dataVectorArray->GetEntriesFast(); idata++){
    dataVector = static_cast<AliQnCorrectionsDataVector*>(dataVectorArray->At(idata));

    width = 1.0;

    mult=dataVector->Weight();

    fillValues[dim-1] = dataVector->Id();    
    bin=inputHistos->EqualizationHistogramM(0)->GetBin(fillValues);

    //cout<<"check"<<endl;
    if(QnConf->ChannelGroups()){
      fillValues[dim-1] = QnConf->ChannelGroup(dataVector->Id());
      binGroup=inputHistos->GroupEqualizationHistogramM(0)->GetBin(fillValues);
      groupWeight  = inputHistos->GroupEqualizationHistogramM(0)->GetBinContent(binGroup);
    }
    //cout<<dataVector->Id()<<"  "<<fillValues[dim-1]<<"  "<<groupWeight<<endl;
    //groupWeight=1.0;
    average  = inputHistos->EqualizationHistogramM(0)->GetBinContent(bin);
    width    = inputHistos->EqualizationHistogramM(0)->GetBinError(bin);


    chanAveMult=0.0;
    chanWidthMult=0.0;
    if(average > 1.0e-6){
      chanAveMult   = mult/average*groupWeight;
      chanWidthMult = (1.+(mult-average)/width*0.1)*groupWeight;
    }
    dataVector->SetAverageEqualizedWeight(chanAveMult );
    dataVector->SetWidthEqualizedWeight(  chanWidthMult);
  }

  return;
}



//_____________________________________________________________________
void AliQnCorrectionsSteps::RecenterQvec(AliQnCorrectionsQnVector* QvectorIn, AliQnCorrectionsQnVector* QvectorOut, AliQnCorrectionsHistograms* inputHistos, Int_t bin, Int_t useStep, Int_t minHar, Int_t maxHar) {

  //
  // Recenter the detector event plane
  //


      Double_t axentries = inputHistos->CalibrationHistogramE(useStep)->GetBinContent(bin);



      for(Int_t ih=minHar; ih<= maxHar; ++ih) {

        if(QvectorIn->CheckQnVectorStatus(ih,AliQnCorrectionsConstants::kUndefined)) continue;

        QvectorOut->SetQx( ih, 
            QvectorIn->Qx(ih) - inputHistos->CalibrationHistogramQ(useStep,ih,0)->GetBinContent(bin)/axentries);
        QvectorOut->SetQy( ih,  
            QvectorIn->Qy(ih) - inputHistos->CalibrationHistogramQ(useStep,ih,1)->GetBinContent(bin)/axentries);

        QvectorOut->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kRecentering);
        if(axentries==1) QvectorOut->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kUndefined);   // With one event in the event class, Qvector becomes undefined after recentering

      }

      // Perform u2n twist correction (TPC)
      //cout<<"recentering fail  "<<QnConf->CalibrationDetectorName(detector)<<"  "<<bin<<"   "<<qVector->Qx(1)<<"  "<<axentries<<"  "<<values[var[0]]<<"  "<<values[var[1]]<<"   "<<values[var[2]]<<endl;
      //if(!(TMath::Abs(qVector->Qx(1))<1.)&&axentries==0) cout<<"recentering fail  "<<QnConf->QnConfigurationName()<<"  "<<qVector->Qx(1)<<"  "<<axentries<<"  "<<values[var[0]]<<"  "<<values[var[1]]<<"   "<<values[var[2]]<<endl;

    //}
  //}
//if(QnConf->GetTwistAndRescalingMethod()==0&&QnConf->CalibrationStep()>2){
//  TwistQnVector(values, QnCorrectionsVarManager::kCentVZERO);
//}



return;
}


//_____________________________________________________________________
void AliQnCorrectionsSteps::TwistAndRescale2nQn(AliQnCorrectionsQnVector* QvectorIn, AliQnCorrectionsQnVector* QvectorTwist, AliQnCorrectionsQnVector* QvectorRescale, AliQnCorrectionsHistograms* inputHistos, Int_t bin, Int_t minHar, Int_t maxHar, Bool_t doTwist, Bool_t doRescaling){

  //
  // Twist and rescaling correction with <Q2n>
  //


  Double_t cos2n, sin2n, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaling, QyRescaling;

  Double_t axentries =inputHistos->CalibrationHistogramE(0)->GetBinContent(bin);
  //Double_t axentries =inputHistos->U2nHistogramE()->GetBinContent(bin);

  for(Int_t ih=minHar; ih<=maxHar; ++ih) {
    cos2n = inputHistos->CalibrationHistogramQ(0,ih*2,0)->GetBinContent(bin)/axentries;
    sin2n = inputHistos->CalibrationHistogramQ(0,ih*2,1)->GetBinContent(bin)/axentries;
    //cos2n = inputHistos->U2nHistogram(ih,0)->GetBinContent(bin)/axentries;
    //sin2n = inputHistos->U2nHistogram(ih,1)->GetBinContent(bin)/axentries;

    
     Ap = 1+cos2n;
     An = 1-cos2n;
     Lp = sin2n/Ap;
     Ln = sin2n/An;


    // if(gPrintOnce<3){
    // cout<<inputHistos->CalibrationHistogramQ(0,ih*2,0)->GetName()<<"  "<<minHar<<" Ap,An,Lp,Ln: "<<Ap<<"  "<<An<<"  "<<Lp<<"  "<<Ln<<" , <Qx>,<Qy>:  "<<inputHistos->CalibrationHistogramQ(0,ih,0)->GetBinContent(bin)/axentries<<"  "<<inputHistos->CalibrationHistogramQ(0,ih,1)->GetBinContent(bin)/axentries<<"  , <Q2x>,<Q2y>:  "<<cos2n<<"  "<<sin2n<<endl;
    // gPrintOnce++;
    // }

    if(!(Lp>-99999999.&&Lp<99999999.)) continue;
    if(!(Ln>-99999999.&&Ln<99999999.)) continue;
    if(!(Ap>-99999999.&&Ap<99999999.)) continue;
    if(!(An>-99999999.&&An<99999999.)) continue;


    if(QvectorIn->CheckQnVectorStatus(ih, AliQnCorrectionsConstants::kUndefined)) continue;

      Qx = QvectorIn->Qx(ih);
      Qy = QvectorIn->Qy(ih);
      QxTwist = Qx;
      QyTwist = Qy;
      if(doTwist){
        QxTwist = (Qx-Ln*Qy)/(1-Ln*Lp);
        QyTwist = (Qy-Lp*Qx)/(1-Ln*Lp);
        QvectorTwist->SetQx( ih, QxTwist);
        QvectorTwist->SetQy( ih, QyTwist);
        //cout<<Qx<<"  "<<QxTwist<<"  "<<Lp<<"  "<<Ln<<"  "<<Ap<<"  "<<An<<"  "<<cos2n<<"  "<<sin2n<<endl;
        QvectorTwist->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kTwist);
        QvectorRescale->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kTwist);
      }
      

      if(doRescaling){
        QxRescaling = QxTwist / Ap;
        QyRescaling = QyTwist / An;
        QvectorRescale->SetQx( ih, QxRescaling);
        QvectorRescale->SetQy( ih, QyRescaling);
        QvectorRescale->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kRescaling);
      }

   
  }
  return;
}



//_____________________________________________________________________
void AliQnCorrectionsSteps::TwistAndRescale3DetectorCorrelation(AliQnCorrectionsQnVector* QvectorIn, AliQnCorrectionsQnVector* QvectorTwist, AliQnCorrectionsQnVector* QvectorRescale, AliQnCorrectionsHistograms* inputHistos, Int_t bin, Int_t minHar, Int_t maxHar, Bool_t doTwist, Bool_t doRescaling, Int_t eventClassParameter){

  //
  // Twist and rescaling correction with <Q2n>
  //



    Double_t Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaling, QyRescaling;
    Double_t x1xt,y1yt,x1yt,x1x2,x2xt,x2yt;


	    for(Int_t ih=minHar; ih<=maxHar; ++ih) {

      //x1xt = inputHistos->CorrelationProf(0, 0, ih, 0, eventClassParameter)->GetBinContent(bin);
      //y1yt = inputHistos->CorrelationProf(0, 0, ih, 3, eventClassParameter)->GetBinContent(bin);

      //x1yt = inputHistos->CorrelationProf(0, 0, ih, 1, eventClassParameter)->GetBinContent(bin);

      //x1x2 = inputHistos->CorrelationProf(0, 1, ih, 0, eventClassParameter)->GetBinContent(bin);
      //x2xt = inputHistos->CorrelationProf(0, 2, ih, 0, eventClassParameter)->GetBinContent(bin);
      //x2yt = inputHistos->CorrelationProf(0, 2, ih, 1, eventClassParameter)->GetBinContent(bin);

      x1xt = inputHistos->CorrelationProf(0, 0, ih, 0, eventClassParameter)->GetBinContent(bin);
      y1yt = inputHistos->CorrelationProf(0, 0, ih, 3, eventClassParameter)->GetBinContent(bin);

      x1yt = inputHistos->CorrelationProf(0, 0, ih, 1, eventClassParameter)->GetBinContent(bin);

      x1x2 = inputHistos->CorrelationProf(0, 2, ih, 0, eventClassParameter)->GetBinContent(bin);
      x2xt = inputHistos->CorrelationProf(0, 1, ih, 0, eventClassParameter)->GetBinContent(bin);
      x2yt = inputHistos->CorrelationProf(0, 1, ih, 1, eventClassParameter)->GetBinContent(bin);



      Lp = x1yt/x1xt;
      Ln = x1yt/y1yt;
      Ap = TMath::Sqrt(2.*x1x2)*x1xt/TMath::Sqrt(x1xt*x2xt+x1yt*x2yt);
      An = TMath::Sqrt(2.*x1x2)*y1yt/TMath::Sqrt(x1xt*x2xt+x1yt*x2yt);

      //cout<<Lp<<"  "<<Ln<<"  "<<x1yt<<"  "<<x1xt<<"  "<<y1yt<<endl;

    if(!(Lp>-99999999.&&Lp<99999999.)) continue;
    if(!(Ln>-99999999.&&Ln<99999999.)) continue;
    if(!(Ap>-99999999.&&Ap<99999999.)) continue;
    if(!(An>-99999999.&&An<99999999.)) continue;


    if(QvectorIn->CheckQnVectorStatus(ih, AliQnCorrectionsConstants::kUndefined)) continue;

      Qx = QvectorIn->Qx(ih);
      Qy = QvectorIn->Qy(ih);
      QxTwist = Qx;
      QyTwist = Qy;
      if(doTwist){
        QxTwist = (Qx-Ln*Qy)/(1-Ln*Lp);
        QyTwist = (Qy-Lp*Qx)/(1-Ln*Lp);
        QvectorTwist->SetQx( ih, QxTwist);
        QvectorTwist->SetQy( ih, QyTwist);
        //cout<<Qx<<"  "<<QxTwist<<"  "<<Lp<<"  "<<Ln<<"  "<<Ap<<"  "<<An<<"  "<<cos2n<<"  "<<sin2n<<endl;
        QvectorTwist->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kTwist);
        QvectorRescale->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kTwist);
      }
      

      if(doRescaling){
        QxRescaling = QxTwist / Ap;
        QyRescaling = QyTwist / An;
        QvectorRescale->SetQx( ih, QxRescaling);
        QvectorRescale->SetQy( ih, QyRescaling);
        QvectorRescale->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kRescaling);
      }

   
  }
  return;
}



////_____________________________________________________________________
//void AliQnCorrectionsSteps::TwoDetectorCorrelationTwistQvec(Float_t* values, Int_t corpar) {
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
//    TProfile * correlationProfiles[3][fgkEPMaxHarmonics][4];
//
//    for(Int_t ic=0; ic<3; ++ic) 
//	    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  QnConf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist;
//
//	  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) {
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
//void AliQnCorrectionsSteps::CorrelationRescalingQvec(Float_t* values, Int_t corpar) {
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
//    TProfile * correlationProfiles[3][fgkEPMaxHarmonics][4];
//
//    for(Int_t ic=0; ic<3; ++ic) 
//	    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  QnConf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxScaled, QyScaled;
//
//	  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) {
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
//void AliQnCorrectionsSteps::CorrelationTwistAndRescalingQvec(Float_t* values, Int_t corpar) {
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
//    TProfile * correlationProfiles[3][fgkEPMaxHarmonics][4];
//
//    for(Int_t ic=1; ic<3; ++ic) 
//      for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  QnConf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//
//	  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) {
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
//void AliQnCorrectionsSteps::ThreeDetectorCorrelationTPCTwistQvec(Float_t* values, Int_t corpar) {
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
//      for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  QnConf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//    Double_t x1xt,y1yt,x1yt,x1x2,x2xt,x2yt;
//
//    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) {
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
//void AliQnCorrectionsSteps::ThreeDetectorCorrelationTPCRescalingQvec(Float_t* values, Int_t corpar) {
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
//    //TProfile * correlationProfiles[3][fgkEPMaxHarmonics][4];
//
//    //for(Int_t ic=0; ic<3; ++ic) 
//    //  for(Int_t ih=1; ih<=fgkEPMaxHarmonics; ++ih) 
//    //   for(Int_t ip=0; ip<4; ++ip) {
//    //  correlationProfiles[ic][ih-1][ip] =  QnConf->CorrelationProfile(ic, ih, ip);
//    //}
//
//    Double_t Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//    Double_t x1xt,y1yt,x1yt,x1x2,x2xt,x2yt,x1y2,y2yt;
//
//    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) {
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
//void AliQnCorrectionsSteps::ThreeDetectorCorrelationTPCTwistAndRescalingQvec(Float_t* values, Int_t corpar) {
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
//    TProfile * correlationProfiles[3][fgkEPMaxHarmonics][4];
//
//    for(Int_t ic=1; ic<3; ++ic) 
//	    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) 
//       for(Int_t ip=0; ip<4; ++ip) {
//      correlationProfiles[ic][ih-1][ip] =  QnConf->CorrelationProfile(ic, ih, ip);
//    }
//
//    Double_t xx, yy, xy, yx, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//    Double_t x1xt,y1yt,x1yt,x1x2,x2xt,x2yt;
//
//	  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) {
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
//void AliQnCorrectionsSteps::TwistQnVector() {

  //
  // Recenter the detector event plane
  //


  //Int_t bin=0;
  //Int_t* var;
  //Int_t maxHarmonic;
  //Int_t dim; 
  //Double_t fillValues[20];


 // AliQnCorrectionsConfiguration* QnConf = 0x0;
 // for(Int_t iconf=0; iconf<NumberOfQnConfigurations(); iconf++){
 //   AliQnCorrectionsConfiguration* QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
 //   if(!QnConf) continue;

 //   if( QnConf->CalibrationStep()<3) continue;
 //   if(QnConf->GetTwistAndRescalingMethod()!=0&&QnConf->GetTwistAndRescalingMethod()!=1) continue;
 //   if( !QnConf->doTwist()) continue;
 //   maxHarmonic = fgkEPMaxHarmonics;



 //   AliQnCorrectionsAxes EPbinning =  QnConf->CalibrationBinning();

 //   const Int_t* var = EPbinning.Var();
 //   const Int_t  dim = EPbinning.Dim();
 //   for(Int_t iav=0; iav<dim; iav++) fillValues[iav] = fDataContainer[var[iav]]; 









 //   //if(QnConf->GetTwistAndRescalingMethod()==1) maxHarmonic = fgkEPMaxHarmonics;
 //   //else maxHarmonic = 2*fgkEPMaxHarmonics;

 //   //TProfile * U2nProfiles[maxHarmonic][2];

 //   //  for(Int_t ih=1; ih<=maxHarmonic; ++ih) 
 //   //   for(Int_t ip=0; ip<2; ++ip) {
 //   //  U2nProfiles[ih-1][ip] =  QnConf->U2nProfile(ih, ip);
 //   //}

 //   Double_t cos2n, sin2n, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, nentries;

 //   for(Int_t ih=QnConf->MinimumHarmonic(); ih<=(Int_t) fgkEPMaxHarmonics; ++ih) {


 //     sin2n = fInputHistograms[iconf]->U2nHistogram(ih,0)->GetBinContent(fInputHistograms[iconf]->U2nHistogram(ih,0)->GetBin(fillValues));
 //     cos2n = fInputHistograms[iconf]->U2nHistogram(ih,1)->GetBinContent(fInputHistograms[iconf]->U2nHistogram(ih,1)->GetBin(fillValues));
 //     nentries = fInputHistograms[iconf]->U2nHistogramE()->GetBinContent(fInputHistograms[iconf]->U2nHistogramE()->GetBin(fillValues));
 //     sin2n/=nentries;
 //     cos2n/=nentries;

 //     Ap = 1+cos2n;
 //     An = 1-cos2n;
 //     Lp = sin2n/Ap;
 //     Ln = sin2n/An;

 //     if(!(Lp>-99999999.&&Lp<99999999.)) continue;
 //     if(!(Ln>-99999999.&&Ln<99999999.)) continue;
 //     if(!(Ap>-99999999.&&Ap<99999999.)) continue;
 //     if(!(An>-99999999.&&An<99999999.)) continue;



 //     AliQnCorrectionsQnVector* qVector=0x0;
 //     //TClonesArray* qvecList = GetQvectors(QnConf->GlobalIndex());
 //     TClonesArray* qvecList = Qvectors(iconf);
 //     for(Int_t ibin=0; ibin<qvecList->GetEntriesFast(); ibin++){
 //       qVector=static_cast<AliQnCorrectionsQnVector*>(qvecList->At(ibin));
 //       //if(!qVector) break;

 //       if(!qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kRecentered)) continue;

 //       Qx = qVector->Qx(ih);
 //       Qy = qVector->Qy(ih);

 //       QxTwist = (Qx-Ln*Qy)/(1-Ln*Lp);
 //       QyTwist = (Qy-Lp*Qx)/(1-Ln*Lp);

 //       qVector->SetQx( ih, QxTwist);
 //       qVector->SetQy( ih, QyTwist);

 //       //cout<<"Twist "<<QnConf->QnConfigurationName()<<endl;
 //       //cout<<Qx<<"  "<<Qy<<"   "<<Lp<<"  "<<Ln<<"  "<<Ap<<"  "<<An<<endl;
 //       qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kDiagonalized);
 //       if(fUseEvent){
 //         fOutputHistograms[iconf]->CalibrationHistogramQ(4,ih,0)->Fill(fillValues,qVector->Qx(ih));
 //         fOutputHistograms[iconf]->CalibrationHistogramQ(4,ih,1)->Fill(fillValues,qVector->Qy(ih));
 //         if(ih==QnConf->MinimumHarmonic()) fOutputHistograms[iconf]->CalibrationHistogramE(4)->Fill(fillValues);
 //       }
 //     }
 //   }
 // }

//return;
//}



//_____________________________________________________________________
//void AliQnCorrectionsSteps::RescaleQnVector(Int_t u2npar) {

  //
  // Recenter the detector event plane
  //


  //Int_t bin=0;
  //Int_t* var;
  //Int_t maxHarmonic;
  //Int_t dim; 


  //Double_t fillValues[20];

  //AliQnCorrectionsConfiguration* QnConf = 0x0;
  //for(Int_t iconf=0; iconf<NumberOfQnConfigurations(); iconf++){
  //  AliQnCorrectionsConfiguration* QnConf = (AliQnCorrectionsConfiguration*) GetQnConfiguration(iconf);
  //  if(!QnConf) continue;

  //  if( QnConf->CalibrationStep()<4) continue;
  //  //if( !QnConf->doScaling()) continue;
  //  if(QnConf->GetTwistAndRescalingMethod()!=0&&QnConf->GetTwistAndRescalingMethod()!=1) continue;
  //  if( !QnConf->doScaling()) continue;
  //  //if(QnConf->GetTwistAndRescalingMethod()==1) maxHarmonic = fgkEPMaxHarmonics;
  //  //else maxHarmonic = 2*fgkEPMaxHarmonics;

  //  //TProfile * U2nProfiles[maxHarmonic][2];

  //  Double_t cos2n, Scos2n, Ncos2n, sin2n, Lp, Ln, Ap, An, Qx, Qy, QxRescaling, QyRescaling;


  //  const Int_t* var = QnConf->CalibrationBinning().Var();
  //  const Int_t  dim = QnConf->CalibrationBinning().Dim();
  //  for(Int_t iv=0; iv<dim; iv++) {fillValues[iv] = fDataContainer[var[iv]];}

  //  for(Int_t ih=1; ih<=(Int_t) fgkEPMaxHarmonics/2.; ++ih) {

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
//void AliQnCorrectionsSteps::U2nTwistAndRescalingQvec(Float_t* values, Int_t u2npar) {
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
//	    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics*2; ++ih) 
//       for(Int_t ip=0; ip<2; ++ip) {
//      U2nProfiles[ih-1][ip] =  QnConf->U2nProfile(ih, ip);
//    }
//
//    Double_t cos2n, sin2n, Lp, Ln, Ap, An, Qx, Qy, QxTwist, QyTwist, QxRescaled, QyRescaled;
//
//	  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) {
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
//
//
//
////_____________________________________________________________________
//void AliQnCorrectionsSteps::RotateQvec() {
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
// for(Int_t idet=0; idet<fNumberOfDetectors; idet++){
//   if(idet==kTPC) continue;
//   AliQnCorrectionsConfiguration* QnConf = 0x0;
//   TIter nextQnConf(fAliQnCorrectionsConfigurations[idet]);
//   while((QnConf=static_cast<AliQnCorrectionsConfiguration*>(nextQnConf()))) {
//    if(!QnConf) continue;
//
//    if(!QnConf->doRotation()) continue;
//    if(QnConf->CalibrationStep()<3) continue;
//    //if( QnConf->CalibrationStep()>3) continue;
//
//    //cout<<QnConf->QnConfigurationName()<<endl;
//    //cout<<QnConf->CorrelationProfile(0,2,0)<<endl;
//    //cout<<QnConf->CorrelationProfile(0,2,0)->GetEntries()<<endl;
//    TProfile * hXX = QnConf->CorrelationProfile(0,2,0);
//    TProfile * hXY = QnConf->CorrelationProfile(0,2,1);
//    TProfile * hYX = QnConf->CorrelationProfile(0,2,2);
//    TProfile * hYY = QnConf->CorrelationProfile(0,2,3);
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
//    //cout<<QnConf->QnConfigurationName()<<"  "<<dphi<<"  "<<edphi<<"  "<<sigphi<<"   "<<XX<<"  "<<XY<<"  "<<YX<<"  "<<YY<<endl;
//
//    Double_t Qx, Qy, Qmag, QxRotated, QyRotated;
//    Double_t x1yt,y1x2, x1x2, y1y2;
//
//    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=fgkEPMaxHarmonics; ++ih) {
//
//      AliQnCorrectionsAxes* EPbinning =  QnConf->CalibrationBinning();
//          
//      AliQnCorrectionsQnVector* qVector=0x0;
//      TClonesArray* qvecList = GetQvectors(QnConf->GlobalIndex());
//      for(Int_t ibin=0; ibin<qvecList->GetEntriesFast(); ibin++){
//        qVector=static_cast<AliQnCorrectionsQnVector*>(qvecList->At(ibin));
//
//        if(!qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kRecentered)) continue;
//        if(qVector->CheckQnVectorStatus(ih, AliQnCorrectionsQnVector::kAligned)) continue;
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
//        //cout<<QnConf->QnConfigurationName()<<"  "<<Qx<<"   "<<QxRotated<<endl;
//        qVector->SetQx( ih, QxRotated);
//        qVector->SetQy( ih, QyRotated);
//
//        qVector->SetQnVectorStatus(ih, AliQnCorrectionsQnVector::kAligned);
//      }
//    }
//  }
// }
//
//  return;
//}
//

