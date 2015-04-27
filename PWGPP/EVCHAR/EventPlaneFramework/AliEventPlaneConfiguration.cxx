/*
***********************************************************
    Event plane class that contains correction setting for specific event plane
    Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
***********************************************************
*/


#ifndef ALIEVENTPLANECONFIGURATION_H
#include "AliEventPlaneConfiguration.h"
#endif

#include "AliEventPlaneHelper.h"
#include "AliEventPlaneBinning.h"
#include "AliEventPlaneQvector.h"
#include "AliEventPlaneCuts.h"
#include "AliEventPlaneVarManager.h"
//#include <AliReducedHistos.h>
#include <TMath.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TRandom3.h>
#include <TArrayS.h>
#include <iostream>


#ifdef ALIEVENTPLANEVARMANAGER_H
#define VAR AliEventPlaneVarManager
#endif

ClassImp(AliEventPlaneConfiguration)



//_______________________________________________________________________________
AliEventPlaneConfiguration::AliEventPlaneConfiguration() :
  TObject(),
  fCalibrationMethod(-1),
  fEqualizationMethod(-1),
  fTwistAndScalingMethod(-1),
  fCalibrationStep(-1),
  fLocalIndex(-1),
  fGlobalIndex(-1),
  fDetectorType(-1),
  fChannelList(0x0),
  fCalibrationDetectorNames(0x0),
  fEqualizationDetectorNames(0x0),
  fEventPlaneDetectorNames(0x0),
  fCorrelationDetectorNames(),
  fCorrelationDetectorIndices(),
  fTrackCuts(0x0),
  fEqualizationHistM(0x0),
  fEqualizationHistE(0x0),
  fEqualizationHistogramsM(),
  fEqualizationHistogramsE(),
  fInputEqualizationHistogramM(0x0),
  fInputEqualizationHistogramE(0x0),
  fCalibrationHistQ(),
  fCalibrationHistE(0x0),
  fCalibrationHistogramsQ(),
  fCalibrationHistogramsE(),
  fInputCalibrationHistogramsQ(),
  fInputCalibrationHistogramsE(),
  fCorrelationProfiles(),
  fU2nProfiles(),
  fEqualizationBinning(0x0),
  fCalibrationBinning(0x0),
  fEqualizationHistPath(0x0),
  fRecenteringHistPath(0x0),
  fCorrelationHistPath(0x0),
  fChannelEqualization(kFALSE),
  fRecenterQvec(kFALSE),
  fRotateQvec(kFALSE),
  fTwistQvec(kFALSE),
  fScaleQvec(kFALSE),
  //fCalibrationBit(0),
  //fCalibrationIntentBit(0),
  fQvectors(0x0)

  //fTrackCuts()
  //fEventPlaneDetectorId(),
{   
  //
  // Constructor
  //

    fQvectors = new TClonesArray("AliEventPlaneQvector", 100000);
    fQvectors->SetName("AliEventPlaneQvector");

    for(Int_t ic=0; ic<2; ++ic){
      fCorrelationDetectorNames[ic]="";
      fCorrelationDetectorIndices[ic]=-1;
    }

    for(Int_t is=0; is<3; ++is){
      fEqualizationHistogramsM[is] = 0x0;
      fEqualizationHistogramsE[is] = 0x0;
    }
    
    for(Int_t is=0; is<6; ++is){
      for(Int_t ih=0; ih<fgkNHarmonics; ++ih){
        for(Int_t ic=0; ic<2; ++ic){
          fCalibrationHistogramsQ[is][ih][ic]=0x0;
        }
      }
    fCalibrationHistogramsE[is]=0x0;
    }
    for(Int_t ih=0; ih<fgkNHarmonics; ++ih){
      for(Int_t ic=0; ic<2; ++ic){
        fU2nHistograms[ih][ic] = 0x0;
    }}
    fU2nHistogramsE = 0x0;
    for(Int_t is=0; is<6; ++is){
      for(Int_t idet=0; idet<3; ++idet){ 
        for(Int_t icomp=0; icomp<4; ++icomp){ 
          for(Int_t ih=0; ih<fgkNHarmonics; ++ih){ 
            fCorrelationProfs[is][idet][ih][icomp] = 0x0;
    }}}}


    fStages[0]="raw";
    fStages[1]="eq";
    fStages[2]="rec";
    fStages[3]="rot";
    fStages[4]="twist";
    fStages[5]="scal";
}



//_______________________________________________________________________________
AliEventPlaneConfiguration::~AliEventPlaneConfiguration()
{
  //
  // De-Constructor
  //
}


//_______________________________________________________________________________
AliEventPlaneConfiguration::AliEventPlaneConfiguration(const AliEventPlaneConfiguration &c) :
  TObject(),
  fCalibrationMethod(c.CalibrationMethod()),
  fEqualizationMethod(c.EqualizationMethod()),
  fTwistAndScalingMethod(c.TwistAndScalingMethod()),
  fCalibrationStep(c.CalibrationStep()),
  fLocalIndex(c.LocalIndex()),
  fGlobalIndex(c.GlobalIndex()),
  fDetectorType(c.DetectorType()),
  fCalibrationDetectorNames(c.CalibrationDetectorName()),
  fEqualizationDetectorNames(c.EqualizationDetectorName()),
  fEventPlaneDetectorNames(c.EventPlaneDetectorName()),
  fChannelEqualization(c.doChannelEqualization()),

  fRecenterQvec(c.doRecentering()),
  fRotateQvec(c.doRotation()),
  fTwistQvec(c.doTwist()),
  fScaleQvec(c.doScaling()),
  fChannelList(c.ChannelList()),
  fTrackCuts(c.TrackCuts()),
  fEqualizationBinning(c.EqualizationBinning()),
  fCalibrationBinning(c.CalibrationBinning()),
  fEqualizationHistPath(c.EqualizationHistPath()),
  fRecenteringHistPath(c.RecenteringHistPath()),
  fCorrelationHistPath(c.CorrelationHistPath())
  //fQvectors(c.Qvectors())
  //fCalibrationBit(0),
  //fCalibrationIntentBit(0),
  //fTrackCuts()
  //fEventPlaneDetectorId(),
{   
  //
  // Constructor
  //
  fCorrelationDetectorNames[0] = c.CorrelationDetectorName(0);
  fCorrelationDetectorNames[1] = c.CorrelationDetectorName(1);
  fCorrelationDetectorIndices[0] = c.CorrelationDetectorIndex(0);
  fCorrelationDetectorIndices[1] = c.CorrelationDetectorIndex(1);
  fEqualizationHistM  = c.EqualizationHistM();
  fEqualizationHistE = c.EqualizationHistE();
  fCalibrationHistE = c.CalibrationHistE();



    for(Int_t is=0; is<6; ++is){
      for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){
        for(Int_t ic=0; ic<2; ++ic){
          fCalibrationHistogramsQ[is][ih-fgkEPMinHarmonics][ic] = c.CalibrationHistogramQ(is,ih, ic);
        }
      }
      fCalibrationHistogramsE[is] = c.CalibrationHistogramE(is);
    }
    for(Int_t is=0; is<6; ++is){
      for(Int_t idet=0; idet<3; ++idet){ 
        for(Int_t icomp=0; icomp<4; ++icomp){ 
          for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){ 
          fCorrelationProfs[is][idet][ih-fgkEPMinHarmonics][icomp] = c.CorrelationProf(is,idet,ih, icomp);
          fCorrelationEpProfs[is][idet][ih-fgkEPMinHarmonics][icomp] = c.CorrelationEpProf(is,idet,ih, icomp);
    }}}}
    for(Int_t is=0; is<3; ++is){
      fEqualizationHistogramsM[is] = c.EqualizationHistogramM(is);
      fEqualizationHistogramsE[is] = c.EqualizationHistogramE(is);
    }
    
    for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){
      for(Int_t ic=0; ic<2; ++ic){
        fU2nHistograms[ih-fgkEPMinHarmonics][ic] = c.U2nHistogram(ih*2,ic);
    }}
    fU2nHistogramsE = c.U2nHistogramE();

}




//_______________________________________________________________________________
void AliEventPlaneConfiguration::Copy(AliEventPlaneConfiguration* epConf, Int_t localIndex, Int_t globalIndex){   
  //
  // Copy AliEventPlaneConfiguration object
  //

  fCalibrationMethod = epConf->CalibrationMethod();
  fEqualizationMethod = epConf->EqualizationMethod();
  fTwistAndScalingMethod = epConf->TwistAndScalingMethod();
  fCalibrationStep = epConf->CalibrationStep();
  fLocalIndex = localIndex;
  fGlobalIndex = globalIndex;
  fDetectorType = epConf->DetectorType();
  fChannelList = epConf->ChannelList();
  fTrackCuts = epConf->TrackCuts();
  fCalibrationDetectorNames = epConf->CalibrationDetectorName();
  fEqualizationDetectorNames = epConf->EqualizationDetectorName();
  fEventPlaneDetectorNames = epConf->EventPlaneDetectorName();
  fCorrelationDetectorNames[0] = epConf->CorrelationDetectorName(0);
  fCorrelationDetectorNames[1] = epConf->CorrelationDetectorName(1);
  fCorrelationDetectorIndices[0] = epConf->CorrelationDetectorIndex(0);
  fCorrelationDetectorIndices[1] = epConf->CorrelationDetectorIndex(1);
  fEqualizationHistM  = epConf->EqualizationHistM();
  fEqualizationHistE = epConf->EqualizationHistE();
  fCalibrationHistE = epConf->CalibrationHistE();
  fEqualizationBinning = epConf->EqualizationBinning();
  fCalibrationBinning = epConf->CalibrationBinning();
  fEqualizationHistPath = epConf->EqualizationHistPath();
  fRecenteringHistPath = epConf->RecenteringHistPath();
  fCorrelationHistPath = epConf->CorrelationHistPath();

  fChannelEqualization = epConf->doChannelEqualization();
  fRecenterQvec = epConf->doRecentering();
  fRotateQvec = epConf->doRotation();
  fTwistQvec = epConf->doTwist();
  fScaleQvec = epConf->doScaling();
  ////fCalibrationBit = epConf->CalibrationBit();
  //fCalibrationIntentBit = epConf->CalibrationIntentBit();

    for(Int_t is=0; is<6; ++is){
      for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){
        for(Int_t ic=0; ic<2; ++ic){
          fCalibrationHistogramsQ[is][ih-fgkEPMinHarmonics][ic] = epConf->CalibrationHistogramQ(is,ih, ic);
        }
      }
      fCalibrationHistogramsE[is] = epConf->CalibrationHistogramE(is);
    }
    for(Int_t is=0; is<6; ++is){
      for(Int_t idet=0; idet<3; ++idet){ 
        for(Int_t icomp=0; icomp<4; ++icomp){ 
          for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){ 
          fCorrelationProfs[is][idet][ih-fgkEPMinHarmonics][icomp] = epConf->CorrelationProf(is,idet,ih, icomp);
          fCorrelationEpProfs[is][idet][ih-fgkEPMinHarmonics][icomp] = epConf->CorrelationEpProf(is,idet,ih, icomp);
    }}}}
    for(Int_t is=0; is<3; ++is){
      fEqualizationHistogramsM[is] = epConf->EqualizationHistogramM(is);
      fEqualizationHistogramsE[is] = epConf->EqualizationHistogramE(is);
    }
    
    for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){
      for(Int_t ic=0; ic<2; ++ic){
        fU2nHistograms[ih-fgkEPMinHarmonics][ic] = epConf->U2nHistogram(ih*2,ic);
    }}
    fU2nHistogramsE = epConf->U2nHistogramE();

    return;
}


//_______________________________________________________________________________
void AliEventPlaneConfiguration::CreateQvectorHistograms(){

  TString title;

  AliEventPlaneBinning* binning = CalibrationBinning();
  Int_t dim = binning->Dim();
  Int_t* var = binning->Var();
  //TArrayD * binLimits =  binning->Array();
  TAxis * binLimits =  binning->Axes();


  for(Int_t idim=0; idim<dim; idim++){title+=Form("%s;",VAR::VarName(var[idim]));}

  for(Int_t is=0; is<(CalibrationStep()+1); ++is){
    for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){ 
       fCalibrationHistogramsQ[is][ih-fgkEPMinHarmonics][0] = AliEventPlaneHelper::AddHistogram( Form("QvecX_%s_h%d_%s", EventPlaneDetectorName().Data(), ih, fStages[is].Data()), Form("QvecX h%d ;%s ", ih, title.Data()), dim, binLimits);
       fCalibrationHistogramsQ[is][ih-fgkEPMinHarmonics][1] = AliEventPlaneHelper::AddHistogram( Form("QvecY_%s_h%d_%s", EventPlaneDetectorName().Data(), ih, fStages[is].Data()), Form("QvecY h%d ;%s ", ih, title.Data()), dim, binLimits);
    }

    fCalibrationHistogramsE[is] = AliEventPlaneHelper::AddHistogram( Form("Qvec_%s_%s_entries", EventPlaneDetectorName().Data(), fStages[is].Data()), Form("Entries %s; %s", title.Data(), fStages[is].Data()), dim, binLimits);
  }



  if(TwistAndScalingMethod()==0){
    for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){ 
       fU2nHistograms[ih-fgkEPMinHarmonics][0] = AliEventPlaneHelper::AddHistogram( Form("cosnphi_%s_h%d", EventPlaneDetectorName().Data(), ih*2), Form("#LTcos(%d#phi)#GT ;%s ", ih*2, title.Data()), dim, binLimits);
       fU2nHistograms[ih-fgkEPMinHarmonics][1] = AliEventPlaneHelper::AddHistogram( Form("sinnphi_%s_h%d", EventPlaneDetectorName().Data(), ih*2), Form("#LTsin(%d#phi)#GT ;%s ", ih*2, title.Data()), dim, binLimits);

    }
    fU2nHistogramsE = AliEventPlaneHelper::AddHistogram( Form("nphi_%s_entries", EventPlaneDetectorName().Data()), Form("Entries; %s", title.Data()), dim, binLimits);
  }





  return;
}

//_______________________________________________________________________________
Bool_t AliEventPlaneConfiguration::IsTrackSelected(Float_t* values) {
  return fTrackCuts->IsSelected(values);
}


//_______________________________________________________________________________
void AliEventPlaneConfiguration::CreateCorrelationHistograms(){

  TString components[4] = {"XX","XY","YX","YY"};
  TString detectors[3] = {"","",""};

  for(Int_t is=0; is<(CalibrationStep()+1); ++is){
    for(Int_t idet=0; idet<3; ++idet){ 
      for(Int_t icomp=0; icomp<4; ++icomp){ 
        for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih){ 
          detectors[0]=EventPlaneDetectorName();
          detectors[1]=CorrelationDetectorName(0);
          detectors[2]=CorrelationDetectorName(1);

          fCorrelationProfs[is][idet][ih-fgkEPMinHarmonics][icomp] = new TProfile(Form("%s_%sx%s_h%d_%s",components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih, fStages[is].Data()), Form("%s %sx%s h%d ;centrality ", components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih), 100,0,100);
          fCorrelationProfs[is][idet][ih-fgkEPMinHarmonics][icomp]->SetDirectory(0);

          fCorrelationEpProfs[is][idet][ih-fgkEPMinHarmonics][icomp] = new TProfile(Form("%s_EP%sx%s_h%d_%s",components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih, fStages[is].Data()), Form("%s EP %sx%s h%d ;centrality ", components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih), 100,0,100);
          fCorrelationEpProfs[is][idet][ih-fgkEPMinHarmonics][icomp]->SetDirectory(0);
        }
      }
    }
  }

  return;
}



//_______________________________________________________________________________
void AliEventPlaneConfiguration::CreateMultiplicityHistograms(){

  TString title;

  if(EqualizationMethod()<0) return;

  AliEventPlaneBinning* binning = EqualizationBinning();
  Int_t dim = binning->Dim();
  Int_t* var = binning->Var();
  //TArrayD * binLimits =  binning->Array();
  TAxis * binLimits =  binning->Axes();


  for(Int_t idim=0; idim<(dim-1); idim++){title+=Form("%s;",VAR::VarName(var[idim]));}
    title+=";channel number";

    fEqualizationHistogramsM[0] = AliEventPlaneHelper::AddHistogram( Form("Mult_%s_raw", EventPlaneDetectorName().Data()), Form("Multiplicity raw ;%s ", title.Data()), dim, binLimits);
    fEqualizationHistogramsE[0] = AliEventPlaneHelper::AddHistogram( Form("Mult_%s_raw_entries", EventPlaneDetectorName().Data()), Form("Entries raw; %s", title.Data()), dim, binLimits);
    if(CalibrationStep()<1) return;
    fEqualizationHistogramsE[1] = AliEventPlaneHelper::AddHistogram( Form("Mult_%s_mean_entries", EventPlaneDetectorName().Data()), Form("Entries raw; %s", title.Data()), dim, binLimits);
    fEqualizationHistogramsE[2] = AliEventPlaneHelper::AddHistogram( Form("Mult_%s_width_entries", EventPlaneDetectorName().Data()), Form("Entries raw; %s", title.Data()), dim, binLimits);
    fEqualizationHistogramsM[1] = AliEventPlaneHelper::AddHistogram( Form("Mult_%s_mean", EventPlaneDetectorName().Data()), Form("Multiplicity divided by mean ;%s ", title.Data()), dim, binLimits);
    fEqualizationHistogramsM[2] = AliEventPlaneHelper::AddHistogram( Form("Mult_%s_width", EventPlaneDetectorName().Data()), Form("Multiplicity divided by mean and corrected for width ;%s ", title.Data()), dim, binLimits);

  return;
}

//__________________________________________________________________
void AliEventPlaneConfiguration::ConnectInputMultiplicityHistograms(TString file) {
  //
  //  Initialize the calibration file and set the calibration histograms
  //

  AliEventPlaneHelper::InitFile(file);

  TString epdet = EventPlaneDetectorName();

  if(CalibrationStep()<1) return;
  SetEqualizationHistogramM( (THnF*)((THn*)AliEventPlaneHelper::GetHistogram("Mult"+epdet, "Mult_"+epdet+"_raw")));//->Clone(Form("%scalib_total", epdet.Data())));
  SetEqualizationHistogramE( (THnF*)((THn*)AliEventPlaneHelper::GetHistogram("Mult"+epdet, "Mult_"+epdet+"_raw_entries")));//->Clone(Form("%scalib_Etotal", epdet.Data())));
  

}


//_____________________________________________________________________
void AliEventPlaneConfiguration::ConnectInputCalibrationHistograms(TString file) {

  AliEventPlaneHelper::InitFile(file);

  TString det = EventPlaneDetectorName();

  for(Int_t is=0; is<2; ++is) {
    if(CalibrationStep()<=is) continue;
  SetCalibrationHistogramE( is, (THnF*)((THn*)AliEventPlaneHelper::GetHistogram("Qvec"+det, "Qvec_"+det+"_"+fStages[is]+"_entries")));//->Clone(Form("_hCentering%s_Qvec%s_entries", fStages[is].Data(),det.Data())));

  for(Int_t ih=fgkEPMinHarmonics; ih<=fgkEPMaxHarmonics; ++ih) 
    for(Int_t iComponent=0; iComponent<2; ++iComponent) {
      SetCalibrationHistogramQ( ih, is, iComponent, (THnF*)((THn*)AliEventPlaneHelper::GetHistogram("Qvec"+det, Form("Qvec%c_%s_h%d_"+fStages[is],(iComponent==0 ? 'X' : 'Y'), det.Data(), ih))));//->Clone(Form("_hCentering_%s_Qvec%c_%s_h%d",fStages[is].Data(), (iComponent==0 ? 'X' : 'Y'), det.Data(), ih)));
      
    }
  }


  if(TwistAndScalingMethod()==0&&CalibrationStep()>2){
    SetU2nHistogramE( (THnF*)((THn*)AliEventPlaneHelper::GetHistogram("Qvec"+det, Form("nphi_%s_entries", det.Data()))));//->Clone(Form("_hCor_%s_u2nE_%s_h%d", det.Data())));
    for(Int_t ih=fgkEPMaxHarmonics; ih<=fgkEPMaxHarmonics; ++ih) {
      SetU2nHistogram( ih, 0, (THnF*)((THn*)AliEventPlaneHelper::GetHistogram("Qvec"+det, Form("cosnphi_%s_h%d", det.Data(), ih*2))));//->Clone(Form("_hCor_%s_u2nX_%s_h%d", det.Data(), ih*2)));
      SetU2nHistogram( ih, 1, (THnF*)((THn*)AliEventPlaneHelper::GetHistogram("Qvec"+det, Form("sinnphi_%s_h%d", det.Data(), ih*2))));//->Clone(Form("_hCor_%s_u2nY_%s_h%d", det.Data(), ih*2)));
      //std::cout<<"load "<<InputU2nHistogram(ih,0)<<std::endl;
      //std::cout<<"load "<<InputU2nHistogram(ih,0)->GetEntries()<<std::endl;
  }}

}
