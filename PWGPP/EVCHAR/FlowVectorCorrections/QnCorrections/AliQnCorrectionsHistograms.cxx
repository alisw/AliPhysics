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
#include "AliQnCorrectionsHistograms.h"
#include "AliQnCorrectionsConfiguration.h"
#include "AliQnCorrectionsAxes.h"
#include <TArrayS.h>
#include <THn.h>
#include <TList.h>
#include <TProfile.h>

ClassImp(AliQnCorrectionsHistograms)


//_______________________________________________________________________________
AliQnCorrectionsHistograms::AliQnCorrectionsHistograms() :
  TObject(),
  fGroupEqualizationHistogramsM(),
  fGroupEqualizationHistogramsE(),
  fEqualizationHistogramsM(),
  fEqualizationHistogramsE(),
  fRotationHistogram(),
  fRotationHistogramE(),
  fU2nHistogramsQA(),
  fU2nHistogramsEQA(),
  fCalibrationHistogramsQ(),
  fCalibrationHistogramsE(),
  fCorrelationProfs(),
  fCorrelationEpProfs(),
  fU2nHistograms(),
  fU2nHistogramsE(0x0),
  fEventHistograms(),
  fStages()
{   
  //
  // Constructor
  //


  for(Int_t is=0; is<3; ++is){
    fEqualizationHistogramsM[is] = 0x0;
    fEqualizationHistogramsE[is] = 0x0;
    fGroupEqualizationHistogramsM[is] = 0x0;
    fGroupEqualizationHistogramsE[is] = 0x0;
  }

  //fInputEqualizationHistogramM=new THnF();

  for(Int_t is=0; is<AliQnCorrectionsConstants::nCorrectionSteps; ++is){
    for(Int_t ih=0; ih<AliQnCorrectionsConstants::nHarmonics; ++ih){
      for(Int_t ic=0; ic<2; ++ic){
        fCalibrationHistogramsQ[is][ih][ic]=0x0;
      }
    }
    fCalibrationHistogramsE[is]=0x0;
  }
  for(Int_t ih=0; ih<AliQnCorrectionsConstants::nHarmonics; ++ih){
    for(Int_t ic=0; ic<2; ++ic){
      fU2nHistograms[ih][ic] = 0x0;
    }}
  for(Int_t idet=0; idet<3; ++idet){ 
    for(Int_t icomp=0; icomp<4; ++icomp){ 
      for(Int_t iaxis=0; iaxis<AliQnCorrectionsConstants::nHistogramDimensions; ++iaxis){ 
        for(Int_t ih=0; ih<AliQnCorrectionsConstants::nHarmonics; ++ih){ 
          for(Int_t is=0; is<6; ++is){
            fCorrelationProfs[is][idet][ih][icomp][iaxis] = 0x0;
            fCorrelationEpProfs[is][idet][ih][icomp][iaxis] = 0x0;
          }}}}}

  for(Int_t ih=0; ih<AliQnCorrectionsConstants::nHistogramDimensions; ++ih){
    fEventHistograms[ih]=0x0;
  }


  fStages[0]="raw";
  fStages[1]="eq";
  fStages[2]="rec";
  fStages[3]="rot";
  fStages[4]="twist";
  fStages[5]="scal";
}



//_______________________________________________________________________________
AliQnCorrectionsHistograms::~AliQnCorrectionsHistograms()
{
  //
  // De-Constructor
  //
}



//_________________________________________________________________
THnF* AliQnCorrectionsHistograms::AddHistogram( const Char_t* name, const Char_t* title, Int_t nDimensions,TAxis* binLimits){
  //
  // create a multi-dimensional histogram THnF with equal or variable bin widths
  //
  if(!binLimits) return 0x0;
  TString hname = name;

  TString titleStr(title);
  TObjArray* arr=titleStr.Tokenize(";");


  Double_t* xmin = new Double_t[nDimensions];
  Double_t* xmax = new Double_t[nDimensions];
  Int_t* nBins = new Int_t[nDimensions];
  for(Int_t idim=0;idim<nDimensions;++idim) {
    nBins[idim] = binLimits[idim].GetNbins();
    xmin[idim] = binLimits[idim].GetBinLowEdge(1);
    xmax[idim] = binLimits[idim].GetBinUpEdge(nBins[idim]);
  }

  THnF* h=new THnF(hname.Data(),title,nDimensions,nBins,xmin,xmax);
  for(Int_t idim=0;idim<nDimensions;++idim) {
    TAxis* axis=h->GetAxis(idim);
    *axis=TAxis(binLimits[idim]);
    axis->SetTitle(arr->At(idim+1)->GetName());
  }

  h->Sumw2();

  delete [] xmin;
  delete [] xmax;
  delete [] nBins;

  return h;
}



//_________________________________________________________________
THnF* AliQnCorrectionsHistograms::CreateHistogram( const Char_t* name, const Char_t* title, AliQnCorrectionsAxes* axes){
  //
  // create a multi-dimensional histogram THnF with equal or variable bin widths
  //


  const Int_t nDimensions = axes->Dim();
  TAxis * binLimits =  axes->Axes();


  TString hname = name;

  TString titleStr(title);
  TObjArray* arr=titleStr.Tokenize(";");

  Double_t* xmin = new Double_t[nDimensions];
  Double_t* xmax = new Double_t[nDimensions];
  Int_t* nBins = new Int_t[nDimensions];
  for(Int_t idim=0;idim<nDimensions;++idim) {
    nBins[idim] = binLimits[idim].GetNbins();
    xmin[idim] = binLimits[idim].GetBinLowEdge(1);
    xmax[idim] = binLimits[idim].GetBinUpEdge(nBins[idim]);
  }

  THnF* h=new THnF(hname.Data(),title,nDimensions,nBins,xmin,xmax);
  for(Int_t idim=0;idim<nDimensions;++idim) {
    TAxis* axis=h->GetAxis(idim);
    *axis=TAxis(binLimits[idim]);
    axis->SetTitle(arr->At(idim+1)->GetName());
  }

  h->Sumw2();

  delete [] xmin;
  delete [] xmax;
  delete [] nBins;

  return h;
}



//_________________________________________________________________
THnF* AliQnCorrectionsHistograms::DivideTHnF(THnF* t1, THnF* t2){

  THnF* h =  (THnF*) t1->CreateHn(t1->GetName(), t1->GetTitle(), t1);

  
  Double_t n,f,f2,e;
  
  for(Int_t i=0; i<h->GetNbins(); i++){
    f=h->GetBinContent(i);
    f2 = h->GetBinError2(i);
    n=t2->GetBinContent(i);

    if(n>1) {f=f/n;f2=f2/n;}                    // these thn's are divided to calculate averages for calibrations, and is designed to fail for n=1
    else {f=0;f2=0;}

    if((f2-f*f)<=0.0) e = TMath::Sqrt(f*f-f2);
    else if((f2-f*f)>0.0) e  = TMath::Sqrt(f2-f*f);
    else e=1.0;

    h->SetBinContent(i,f);
    h->SetBinError(i,e);


  }


return h;


}


//____________________________________________________________________________________
TObject* AliQnCorrectionsHistograms::GetHistogram(TList* list, const Char_t* listname, const Char_t* hname) {
  //
  // Retrieve a histogram from the list hlist
  //
  if(list->FindObject(listname)) return list->FindObject(listname)->FindObject(hname);
}



//_______________________________________________________________________________
void AliQnCorrectionsHistograms::CreateQvectorHistograms(AliQnCorrectionsConfiguration* QnConf){
  
  if(!QnConf->IsRequestedFillHistogram(AliQnCorrectionsConstants::kRecentering)) return;

  TString title="";

  AliQnCorrectionsAxes* binning = QnConf->GetRecenteringAxes();

  for(Int_t idim=0; idim<binning->Dim(); idim++) title+=+";"+binning->AxisLabel(idim);

  for(Int_t is=0; is<AliQnCorrectionsConstants::kNcorrectionSteps; ++is){
    Int_t maxHar = QnConf->MaximumHarmonic();
    if(QnConf->GetTwistAndRescalingMethod()==0) maxHar*=2;
    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=maxHar; ++ih){
       fCalibrationHistogramsQ[is][ih-1][0] =   CreateHistogram( Form("QvecX_%s_h%d_%s", QnConf->QnConfigurationName().Data(), ih, fStages[is].Data()), Form("QvecX h%d %s ", ih, title.Data()), binning);
       fCalibrationHistogramsQ[is][ih-1][1] =   CreateHistogram( Form("QvecY_%s_h%d_%s", QnConf->QnConfigurationName().Data(), ih, fStages[is].Data()), Form("QvecY h%d %s ", ih, title.Data()), binning);
    }

    fCalibrationHistogramsE[is] =   CreateHistogram( Form("Qvec_%s_%s_entries", QnConf->QnConfigurationName().Data(), fStages[is].Data()), Form("Entries %s %s", title.Data(), fStages[is].Data()), binning);
  }



  if(QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kTwist)||QnConf->IsRequestedCorrection(AliQnCorrectionsConstants::kRescaling)){
  binning = QnConf->GetTwistAndRescalingAxes();
    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){ 
       fU2nHistograms[ih-1][0] =   CreateHistogram( Form("cosnphi_%s_h%d", QnConf->QnConfigurationName().Data(), ih*2), Form("#LTcos(%d#phi)#GT %s ", ih*2, title.Data()), binning);
       fU2nHistograms[ih-1][1] =   CreateHistogram( Form("sinnphi_%s_h%d", QnConf->QnConfigurationName().Data(), ih*2), Form("#LTsin(%d#phi)#GT %s ", ih*2, title.Data()), binning);

    }
    fU2nHistogramsE =   CreateHistogram( Form("nphi_%s_entries", QnConf->QnConfigurationName().Data()), Form("Entries %s", title.Data()), binning);
  }





  return;
}


//_______________________________________________________________________________
void AliQnCorrectionsHistograms::CreateCorrelationHistograms(AliQnCorrectionsConfiguration* QnConf){

  CreateCorrelationHistogramsQA(QnConf);

  if(QnConf->QnConfigurationCorrelationIndex(0)==-1||QnConf->QnConfigurationCorrelationIndex(0)==-1) return;
  if(QnConf->GetTwistAndRescalingMethod()!=2) return;
  if(!QnConf->IsRequestedFillHistogram(AliQnCorrectionsConstants::kTwist)) return;
  //if(QnConf->GetCorrectionStep(currentManagerStep+1)<((Int_t)AliQnCorrectionsConstants::kTwist)) return;

  TString components[4] = {"XX","XY","YX","YY"};
  TString detectors[3] = {"","",""};

  TAxis ax;

    for(Int_t idet=0; idet<3; ++idet){ 
      for(Int_t icomp=0; icomp<4; ++icomp){ 
        for(Int_t iaxis=0; iaxis<QnConf->GetAlignmentAxes()->Dim(); ++iaxis){ 
          for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){ 
            ax= TAxis(QnConf->GetAlignmentAxes()->Axis(iaxis));
            detectors[0]=QnConf->QnConfigurationName();
            detectors[1]=QnConf->QnConfigurationCorrelationName(0);
            detectors[2]=QnConf->QnConfigurationCorrelationName(1);

            fCorrelationProfs[0][idet][ih-1][icomp][iaxis] = new TProfile(Form("Correlation_%s_%sx%s_h%d_%s",components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih, QnConf->GetAlignmentAxes()->AxisLabel(iaxis).Data()), Form("Correlation_%s_%sx%s_h%d_%s",components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih, QnConf->GetAlignmentAxes()->AxisLabel(iaxis).Data()), 1,0.,1.);
            fCorrelationProfs[0][idet][ih-1][icomp][iaxis]->SetBins(ax.GetNbins(),ax.GetXbins()->GetArray());
            fCorrelationProfs[0][idet][ih-1][icomp][iaxis]->SetDirectory(0);

            fCorrelationProfs[1][idet][ih-1][icomp][iaxis] = new TProfile(Form("Correlation_%s_%sx%s_h%d_%s_Twist",components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih, QnConf->GetAlignmentAxes()->AxisLabel(iaxis).Data()), Form("Correlation_%s_%sx%s_h%d_%s",components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih, QnConf->GetAlignmentAxes()->AxisLabel(iaxis).Data()), 1,0.,1.);
            fCorrelationProfs[1][idet][ih-1][icomp][iaxis]->SetBins(ax.GetNbins(),ax.GetXbins()->GetArray());
            fCorrelationProfs[1][idet][ih-1][icomp][iaxis]->SetDirectory(0);

            fCorrelationProfs[2][idet][ih-1][icomp][iaxis] = new TProfile(Form("Correlation_%s_%sx%s_h%d_%s_Rescaled",components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih, QnConf->GetAlignmentAxes()->AxisLabel(iaxis).Data()), Form("Correlation_%s_%sx%s_h%d_%s",components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih, QnConf->GetAlignmentAxes()->AxisLabel(iaxis).Data()), 1,0.,1.);
            fCorrelationProfs[2][idet][ih-1][icomp][iaxis]->SetBins(ax.GetNbins(),ax.GetXbins()->GetArray());
            fCorrelationProfs[2][idet][ih-1][icomp][iaxis]->SetDirectory(0);

          }

        }
      }
    }


  return;
}



//_______________________________________________________________________________
void AliQnCorrectionsHistograms::CreateCorrelationHistogramsQA(AliQnCorrectionsConfiguration* QnConf){

  if(QnConf->QnConfigurationCorrelationIndex(0)==-1||QnConf->QnConfigurationCorrelationIndex(0)==-1) return;

  TString components[4] = {"XX","XY","YX","YY"};
  TString detectors[3] = {"","",""};

  TAxis ax;

  for(Int_t is=0; is<=QnConf->CalibrationStep(); ++is){
    for(Int_t idet=0; idet<3; ++idet){ 
      for(Int_t icomp=0; icomp<4; ++icomp){ 
        for(Int_t iaxis=0; iaxis<QnConf->CalibrationBinning()->Dim(); ++iaxis){ 
          for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){ 
            ax= TAxis(QnConf->CalibrationBinning()->Axis(iaxis));
            detectors[0]=QnConf->QnConfigurationName();
            detectors[1]=QnConf->QnConfigurationCorrelationName(0);
            detectors[2]=QnConf->QnConfigurationCorrelationName(1);


            fCorrelationEpProfs[is][idet][ih-1][icomp][iaxis] = new TProfile(Form("%s_%sx%s_h%d_%s_%s",components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih, QnConf->CalibrationBinning()->AxisLabel(iaxis).Data(), fStages[is].Data()), Form("%s %sx%s h%d ;%s; ", components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih, QnConf->CalibrationBinning()->AxisLabel(iaxis).Data()), 1,0.,1.);
            fCorrelationEpProfs[is][idet][ih-1][icomp][iaxis]->SetBins(ax.GetNbins(),ax.GetXbins()->GetArray());
            fCorrelationEpProfs[is][idet][ih-1][icomp][iaxis]->SetDirectory(0);
          }

        }
      }
    }
  }


  return;
}







//_______________________________________________________________________________
void AliQnCorrectionsHistograms::CreateMultiplicityHistograms(AliQnCorrectionsConfiguration* QnConf){

  TString title="";
  TString det=QnConf->QnConfigurationName();

  if(!QnConf->IsRequestedFillHistogram(AliQnCorrectionsConstants::kDataVectorEqualization)) return;

  AliQnCorrectionsAxes* binning = QnConf->EqualizationBinning();
  Int_t dim = binning->Dim();
  TAxis * binLimits =  binning->Axes();

  for(Int_t idim=0; idim<(dim-1); idim++) title+=";"+binning->AxisLabel(idim);
    title+=";channel number;";

    fEqualizationHistogramsM[0] =   AddHistogram( Form("Mult_%s_raw", det.Data()), Form("Multiplicity raw %s ", title.Data()), dim, binLimits);
    fEqualizationHistogramsE[0] =   AddHistogram( Form("Mult_%s_raw_entries", det.Data()), Form("Entries raw %s", title.Data()), dim, binLimits);
    //if(QnConf->CalibrationStep()<1) return;
    fEqualizationHistogramsE[1] =   AddHistogram( Form("Mult_%s_mean_entries", det.Data()), Form("Entries raw %s", title.Data()), dim, binLimits);
    fEqualizationHistogramsE[2] =   AddHistogram( Form("Mult_%s_width_entries", det.Data()), Form("Entries raw %s", title.Data()), dim, binLimits);
    fEqualizationHistogramsM[1] =   AddHistogram( Form("Mult_%s_mean", det.Data()), Form("Multiplicity divided by mean %s ", title.Data()), dim, binLimits);
    fEqualizationHistogramsM[2] =   AddHistogram( Form("Mult_%s_width", det.Data()), Form("Multiplicity divided by mean and corrected for width %s ", title.Data()), dim, binLimits);

  return;
}




//_______________________________________________________________________________
void AliQnCorrectionsHistograms::CreateRotationHistograms(AliQnCorrectionsConfiguration* QnConf){


  if(!QnConf->IsRequestedFillHistogram(AliQnCorrectionsConstants::kAlignment)) return;

  TString title="";
  TString det=QnConf->QnConfigurationName();

  AliQnCorrectionsAxes* binning = QnConf->GetAlignmentAxes();


  for(Int_t idim=0; idim<binning->Dim(); idim++) title+=";"+binning->AxisLabel(idim);

  fRotationHistogram[0][0]  =   CreateHistogram( Form("Correlation_%sx%s_XX", det.Data(), QnConf->GetReferenceQnForAlignment().Data()), Form("Correlation XX %s", title.Data()), binning);
  fRotationHistogram[0][1]  =   CreateHistogram( Form("Correlation_%sx%s_YY", det.Data(), QnConf->GetReferenceQnForAlignment().Data()), Form("Correlation YY %s", title.Data()), binning);
  fRotationHistogram[0][2]  =   CreateHistogram( Form("Correlation_%sx%s_XY", det.Data(), QnConf->GetReferenceQnForAlignment().Data()), Form("Correlation XY %s", title.Data()), binning);
  fRotationHistogram[0][3]  =   CreateHistogram( Form("Correlation_%sx%s_YX", det.Data(), QnConf->GetReferenceQnForAlignment().Data()), Form("Correlation YX %s", title.Data()), binning);
  fRotationHistogramE[0] =   CreateHistogram( Form("Correlation_%sx%s_entries", det.Data(), QnConf->GetReferenceQnForAlignment().Data()), Form("Correlation entries %s", title.Data()), binning);

  fRotationHistogram[1][0]  =   CreateHistogram( Form("Correlation_%sx%s_XX_rot", det.Data(), QnConf->GetReferenceQnForAlignment().Data()), Form("Correlation XX %s", title.Data()), binning);
  fRotationHistogram[1][1]  =   CreateHistogram( Form("Correlation_%sx%s_YY_rot", det.Data(), QnConf->GetReferenceQnForAlignment().Data()), Form("Correlation YY %s", title.Data()), binning);
  fRotationHistogram[1][2]  =   CreateHistogram( Form("Correlation_%sx%s_XY_rot", det.Data(), QnConf->GetReferenceQnForAlignment().Data()), Form("Correlation XY %s", title.Data()), binning);
  fRotationHistogram[1][3]  =   CreateHistogram( Form("Correlation_%sx%s_YX_rot", det.Data(), QnConf->GetReferenceQnForAlignment().Data()), Form("Correlation YX %s", title.Data()), binning);
  fRotationHistogramE[1] =   CreateHistogram( Form("Correlation_%sx%s_entries_rot", det.Data(), QnConf->GetReferenceQnForAlignment().Data()), Form("Correlation entries %s", title.Data()), binning);

  return;
}



//_______________________________________________________________________________
void AliQnCorrectionsHistograms::CreateEventHistograms(AliQnCorrectionsConfiguration* QnConf){

  const Int_t dim  = QnConf->EqualizationBinning()->Dim();
 
  TString title;
  for(Int_t idim=0; idim<(dim-1); idim++){
    title=QnConf->EqualizationBinning()->AxisLabel(0);
    fEventHistograms[idim]=new TH1F(Form("EventVar%dDistrubution",idim),title,1000,QnConf->EqualizationBinning()->GetLowEdge(idim), QnConf->EqualizationBinning()->GetUpEdge(idim));
    fEventHistograms[idim]->SetDirectory(0);
  }
  

  return;
}


//_______________________________________________________________________________
void AliQnCorrectionsHistograms::FillEventHistograms(Int_t dim, Double_t* fillValues){

 
  for(Int_t idim=0; idim<(dim-1); idim++){
    fEventHistograms[idim]->Fill(fillValues[idim]);
  }
  

  return;
}


//__________________________________________________________________
void AliQnCorrectionsHistograms::CreateGroupMultiplicityHistograms(TList* list, AliQnCorrectionsConfiguration* QnConf) {
  //
  // Creates group multiplicity histograms by summing channels from the multiplicity histograms
  //

  AliQnCorrectionsAxes* binning = QnConf->EqualizationBinning();
  const Int_t dim = binning->Dim();

  Int_t nGroups=1;
  const Int_t nChannelBins=binning->Nbins(dim-1);
  for(Int_t ich=0; ich<nChannelBins; ich++) if(QnConf->ChannelGroup(ich)>nGroups) nGroups=QnConf->ChannelGroup(ich);
  TAxis ax = TAxis(nGroups+1, -0.5, nGroups+0.5);

  TArrayI nChannelsPerGroup = TArrayI(nGroups+1);
  for(Int_t ich=0; ich<nChannelBins; ich++) nChannelsPerGroup.SetAt(nChannelsPerGroup.At(QnConf->ChannelGroup(ich))+1, QnConf->ChannelGroup(ich));
  
  AliQnCorrectionsAxes groupBinning = AliQnCorrectionsAxes(*binning);
  groupBinning.SetAxis(dim-1, 0, ax, "Multiplicity group");
  TAxis* groupLimits = groupBinning.Axes();


  TString title="";
  TString det=QnConf->QnConfigurationName();
  for(Int_t idim=0; idim<(dim-1); idim++) title+=";"+binning->AxisLabel(idim);
  title+=";channel group";

  fGroupEqualizationHistogramsM[0] =   AddHistogram( Form("MultGroup_%s_raw", det.Data()), Form("Multiplicity group raw %s ", title.Data()), dim, groupLimits);
  fGroupEqualizationHistogramsE[0] =   AddHistogram( Form("MultGroup_%s_raw_entries", det.Data()), Form("Entries group raw %s", title.Data()), dim, groupLimits);


  Int_t ndim = fEqualizationHistogramsM[0]->GetNdimensions();
  TAxis* axes[AliQnCorrectionsConstants::nHistogramDimensions]={0x0};
  Int_t nbin[AliQnCorrectionsConstants::nHistogramDimensions]={0};
  Int_t bin[AliQnCorrectionsConstants::nHistogramDimensions];

  for(Int_t id=0; id<ndim; id++){

    axes[id] = fEqualizationHistogramsM[0]->GetAxis(id);
    nbin[id] = axes[id]->GetNbins()+1;
    bin[id]=1;

  }

  Double_t values[AliQnCorrectionsConstants::nHistogramDimensions]={0.0};
  Double_t m,nentries;


  // the following code loops over the bins in the multiplicity histograms and adds the bin content to the associated group multiplicity bin in the group multiplicity histograms
  Int_t id=ndim-1;
  do {
    for(Int_t ibin=0; ibin<nbin[id]-1; ibin++) {

      for(Int_t id2=0; id2<ndim; id2++)  values[id2]=axes[id2]->GetBinCenter(bin[id2]);

      nentries=fEqualizationHistogramsE[0]->GetBinContent(fEqualizationHistogramsE[0]->GetBin(values));
      m=fEqualizationHistogramsM[0]->GetBinContent(fEqualizationHistogramsM[0]->GetBin(values));

      values[ndim-1]=QnConf->ChannelGroup((Int_t) values[ndim-1]);

      fGroupEqualizationHistogramsE[0]->Fill(values,nentries);
      fGroupEqualizationHistogramsM[0]->Fill(values,m);

      bin[id]++;
    }

    for(Int_t id2=0; id2<ndim; id2++) if(id2!=0&&bin[id2]>=(nbin[id2]-1)) {
      bin[id2-1]++;
      if(id2!=1&&bin[id2-1]==nbin[id2-1]) bin[id2-1]=1;
      for(Int_t id3=id2; id3<ndim; id3++) bin[id3]=1;
    }
  } while (bin[0]!=(nbin[0]));






}


//__________________________________________________________________
Bool_t AliQnCorrectionsHistograms::ConnectMultiplicityHistograms(TList* list, AliQnCorrectionsConfiguration* QnConf) {
  //
  //  Initialize the calibration file and set the calibration histograms
  //
  if(!list) return kFALSE;

  TString epdet = QnConf->QnConfigurationName();

  fEqualizationHistogramsM[0]  = (THnF*)((THn*)  GetHistogram(list, "Mult"+epdet, "Mult_"+epdet+"_raw"));
  fEqualizationHistogramsE[0]  = (THnF*)((THn*)  GetHistogram(list, "Mult"+epdet, "Mult_"+epdet+"_raw_entries"));

  if(!fEqualizationHistogramsM[0]) return kFALSE;
  if(fEqualizationHistogramsM[0]->GetEntries()==0) return kFALSE;

  if(QnConf->ChannelGroups()){
    CreateGroupMultiplicityHistograms(list,QnConf);
    SetGroupEqualizationHistogramM(  DivideTHnF(fGroupEqualizationHistogramsM[0],fGroupEqualizationHistogramsE[0]),0);
  }

  SetEqualizationHistogramM(  DivideTHnF(fEqualizationHistogramsM[0],fEqualizationHistogramsE[0]),0);
 
  return kTRUE;

}





//_____________________________________________________________________
Bool_t AliQnCorrectionsHistograms::ConnectMeanQnCalibrationHistograms(TList* list, AliQnCorrectionsConfiguration* QnConf) {

  if(!list) return kFALSE;
  TString det = QnConf->QnConfigurationName();

  Int_t is=-1;

    //if(QnConf->CalibrationStep()<=is) continue;
    //if(QnConf->GetDataVectorEqualizationMethod()<0&&is==1) continue;


  if(QnConf->GetDataVectorEqualizationMethod()<0) is=0;
  else is=1;
  
  SetCalibrationHistogramE( is, (THnF*)((THn*)  GetHistogram(list, "Qvec"+det, "Qvec_"+det+"_"+fStages[is]+"_entries")));//->Clone(Form("_hCentering%s_Qvec%s_entries", fStages[is].Data(),det.Data())));

  if(!fCalibrationHistogramsE[is]) return kFALSE;
  if(fCalibrationHistogramsE[is]->GetEntries()==0) return kFALSE;

  Int_t maxHar = QnConf->MaximumHarmonic();
  if(QnConf->GetTwistAndRescalingMethod()==0) maxHar*=2;
  for(Int_t ih=QnConf->MinimumHarmonic(); ih<=maxHar; ++ih) {
    for(Int_t iComponent=0; iComponent<2; ++iComponent) {
      SetCalibrationHistogramQ( ih, is, iComponent, (THnF*)((THn*)  GetHistogram(list, "Qvec"+det, Form("Qvec%c_%s_h%d_"+fStages[is],(iComponent==0 ? 'X' : 'Y'), det.Data(), ih))));//->Clone(Form("_hCentering_%s_Qvec%c_%s_h%d",fStages[is].Data(), (iComponent==0 ? 'X' : 'Y'), det.Data(), ih)));
      }
      
    }

  return kTRUE;
}



//_____________________________________________________________________
Bool_t AliQnCorrectionsHistograms::ConnectU2nQnCalibrationHistograms(TList* list, AliQnCorrectionsConfiguration* QnConf) {

  if(!list) return kFALSE;

  TString det = QnConf->QnConfigurationName();

    SetU2nHistogramE( (THnF*)((THn*)  GetHistogram(list, "Qvec"+det, Form("nphi_%s_entries", det.Data()))));//->Clone(Form("_hCor_%s_u2nE_%s_h%d", det.Data())));

    if(!fU2nHistogramsE) return kFALSE;
    if(fU2nHistogramsE->GetEntries()==0) return kFALSE;

    for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih) {
      SetU2nHistogram( ih, 0, (THnF*)((THn*)  GetHistogram(list, "Qvec"+det, Form("cosnphi_%s_h%d", det.Data(), ih*2))));//->Clone(Form("_hCor_%s_u2nX_%s_h%d", det.Data(), ih*2)));
      SetU2nHistogram( ih, 1, (THnF*)((THn*)  GetHistogram(list, "Qvec"+det, Form("sinnphi_%s_h%d", det.Data(), ih*2))));//->Clone(Form("_hCor_%s_u2nY_%s_h%d", det.Data(), ih*2)));
  }

  

  return kTRUE;
}



//_____________________________________________________________________
Bool_t AliQnCorrectionsHistograms::ConnectCorrelationQnCalibrationHistograms(TList* list, AliQnCorrectionsConfiguration* QnConf) {

              
  TString components[4] = {"XX","XY","YX","YY"};
  TString detectors[3] = {"","",""};

  TAxis ax;
              

  if(!list) return kFALSE;

  TString det = QnConf->QnConfigurationName();

    for(Int_t idet=0; idet<3; ++idet){ 
      for(Int_t icomp=0; icomp<4; ++icomp){ 
        for(Int_t iaxis=0; iaxis<QnConf->GetAlignmentAxes()->Dim(); ++iaxis){ 
          for(Int_t ih=QnConf->MinimumHarmonic(); ih<=QnConf->MaximumHarmonic(); ++ih){ 
            ax= TAxis(QnConf->GetAlignmentAxes()->Axis(iaxis));
            detectors[0]=QnConf->QnConfigurationName();
            detectors[1]=QnConf->QnConfigurationCorrelationName(0);
            detectors[2]=QnConf->QnConfigurationCorrelationName(1);


            fCorrelationProfs[0][idet][ih-1][icomp][iaxis] = (TProfile*) GetHistogram(list, "Correlations"+det,Form("Correlation_%s_%sx%s_h%d_%s",components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih, QnConf->GetAlignmentAxes()->AxisLabel(iaxis).Data()));

            //cout<<det<<"  "<<fCorrelationProfs[0][idet][ih-1][icomp][iaxis]<<"  "<<Form("Correlation_%s_%sx%s_h%d_%s",components[icomp].Data(), detectors[idet].Data(), detectors[(idet+1)%3].Data(), ih, QnConf->GetAlignmentAxes()->AxisLabel(iaxis).Data())<<endl;
            if(!fCorrelationProfs[0][idet][ih-1][icomp][iaxis]) return kFALSE;
            if(fCorrelationProfs[0][idet][ih-1][icomp][iaxis]->GetEntries()==0) return kFALSE;
  
     }}}}

  return kTRUE;
}



//_____________________________________________________________________
Bool_t AliQnCorrectionsHistograms::ConnectRotationQnCalibrationHistograms(TList* list, AliQnCorrectionsConfiguration* QnConf) {

  if(!list) return kFALSE;

  TString det = QnConf->QnConfigurationName();

    SetRotationHistogramE( (THnF*)((THn*)  GetHistogram(list, "Correlations"+det, Form("Correlation_%sx%s_entries", det.Data(), QnConf->GetReferenceQnForAlignment().Data()))));//->Clone(Form("_hCor_%s_u2nE_%s_h%d", det.Data())));


    if(!fRotationHistogramE[0]) return kFALSE;
    if(fRotationHistogramE[0]->GetEntries()==0) return kFALSE;

    SetRotationHistogram( 0, (THnF*)((THn*)  GetHistogram(list,"Correlations"+det, Form("Correlation_%sx%s_XX", det.Data(), QnConf->GetReferenceQnForAlignment().Data()) )));
    SetRotationHistogram( 1, (THnF*)((THn*)  GetHistogram(list,"Correlations"+det, Form("Correlation_%sx%s_YY", det.Data(), QnConf->GetReferenceQnForAlignment().Data()) )));
    SetRotationHistogram( 2, (THnF*)((THn*)  GetHistogram(list,"Correlations"+det, Form("Correlation_%sx%s_XY", det.Data(), QnConf->GetReferenceQnForAlignment().Data()) )));
    SetRotationHistogram( 3, (THnF*)((THn*)  GetHistogram(list,"Correlations"+det, Form("Correlation_%sx%s_YX", det.Data(), QnConf->GetReferenceQnForAlignment().Data()) )));

  

  return kTRUE;
}
