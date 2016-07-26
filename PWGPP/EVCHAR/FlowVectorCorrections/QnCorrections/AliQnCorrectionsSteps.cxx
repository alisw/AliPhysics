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




//_____________________________________________________________________
void AliQnCorrectionsSteps::RotateQvec(AliQnCorrectionsQnVector* QvectorIn, AliQnCorrectionsQnVector* QvectorOut, AliQnCorrectionsHistograms* inputHistos, Int_t bin, Int_t minHar, Int_t maxHar, Int_t alignmentHarmonic) {

  //
  // Align Q-vectors
  //


    Double_t XX  = inputHistos->GetRotationHistogram(0,0)->GetBinContent(bin);
    Double_t YY  = inputHistos->GetRotationHistogram(0,1)->GetBinContent(bin);
    Double_t XY  = inputHistos->GetRotationHistogram(0,2)->GetBinContent(bin);
    Double_t YX  = inputHistos->GetRotationHistogram(0,3)->GetBinContent(bin);
    Double_t eXX = inputHistos->GetRotationHistogram(0,0)->GetBinError(bin);
    Double_t eYY = inputHistos->GetRotationHistogram(0,1)->GetBinError(bin);
    Double_t eXY = inputHistos->GetRotationHistogram(0,2)->GetBinError(bin);
    Double_t eYX = inputHistos->GetRotationHistogram(0,3)->GetBinError(bin);
    Double_t n   = inputHistos->GetRotationHistogramE(0)->GetBinContent(bin);

    eXX =  TMath::Sqrt((eXX*eXX/n-(XX/n*XX/n))/n);
    eYY =  TMath::Sqrt((eYY*eYY/n-(YY/n*YY/n))/n);
    eXY =  TMath::Sqrt((eXY*eXY/n-(XY/n*XY/n))/n);
    eYX =  TMath::Sqrt((eYX*eYX/n-(YX/n*YX/n))/n);

    Double_t dphi = -TMath::ATan((XY-YX)/(XX+YY))/alignmentHarmonic;
    Double_t edenom2 = eXY*eXY+eYX*eYX;
    Double_t enumer2 = eXX*eXX+eYY*eYY;

    for(Int_t ih=minHar; ih<=maxHar; ++ih) {QvectorOut->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kAlignment);}

    if(TMath::Sqrt((XY-YX)*(XY-YX)/edenom2)<2.) return;
    //if(!(dphi<1.)) return;

    Double_t equot  = TMath::Sqrt(enumer2/(XX+YY)/(XX+YY)+edenom2/(XY-YX)/(XY-YX))*((XY-YX)/(XX+YY));
    Double_t edphi = equot/(1.+(XY-YX)/(XX+YY)*(XY-YX)/(XX+YY));
    Double_t sigphi = TMath::Abs(dphi/edphi);
    //dphi=-dphi;

    //cout<<QnConf->QnConfigurationName()<<"  "<<dphi<<"  "<<edphi<<"  "<<sigphi<<"   "<<XX<<"  "<<XY<<"  "<<YX<<"  "<<YY<<endl;

    Double_t Qx, Qy, Qmag, QxRotated, QyRotated;
    Double_t x1yt,y1x2, x1x2, y1y2;

    for(Int_t ih=minHar; ih<=maxHar; ++ih) {

        Qx = QvectorIn->Qx(ih);
        Qy = QvectorIn->Qy(ih);

        QvectorOut->SetQx(ih,Qx*TMath::Cos(((Double_t) ih)*dphi)+Qy*TMath::Sin(((Double_t) ih)*dphi));
        QvectorOut->SetQy(ih,Qy*TMath::Cos(((Double_t) ih)*dphi)-Qx*TMath::Sin(((Double_t) ih)*dphi));


        QvectorOut->SetQnVectorStatus(ih, AliQnCorrectionsConstants::kAlignment);
    }


  return;
}


