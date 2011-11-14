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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Class providing the calculation of derived quantities (mean,rms,fits,...) //
//       of calibration entries                                              //
/*


*/
////////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TVectorT.h>
#include <TObjArray.h>
#include <TGraph.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TMap.h>
#include <TGraphErrors.h>
#include <AliCDBStorage.h>
#include <AliDCSSensorArray.h>
#include <AliTPCSensorTempArray.h>
#include <AliDCSSensor.h>
#include <AliLog.h>
#include <AliCDBEntry.h>
#include <AliCDBManager.h>
#include <AliCDBId.h>
#include <AliSplineFit.h>
#include "AliTPCcalibDB.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "AliTPCROC.h"
#include "AliTPCmapper.h"
#include "AliTPCParam.h"
#include "AliTPCCalibRaw.h"
#include "AliTPCPreprocessorOnline.h"
#include "AliTPCdataQA.h"
#include "AliLog.h"
#include "AliTPCcalibDButil.h"
#include "AliTPCCalibVdrift.h"
#include "AliMathBase.h"
#include "AliRelAlignerKalman.h"

const Float_t kAlmost0=1.e-30;

ClassImp(AliTPCcalibDButil)
AliTPCcalibDButil::AliTPCcalibDButil() :
  TObject(),
  fCalibDB(0),
  fPadNoise(0x0),
  fPedestals(0x0),
  fPulserTmean(0x0),
  fPulserTrms(0x0),
  fPulserQmean(0x0),
  fPulserOutlier(new AliTPCCalPad("PulserOutliers","PulserOutliers")),
  fCETmean(0x0),
  fCETrms(0x0),
  fCEQmean(0x0),
  fALTROMasked(0x0),
  fCalibRaw(0x0),
  fDataQA(0x0),
  fRefMap(0x0),
  fCurrentRefMap(0x0),
  fRefValidity(""),
  fRefPadNoise(0x0),
  fRefPedestals(0x0),
  fRefPedestalMasked(0x0),
  fRefPulserTmean(0x0),
  fRefPulserTrms(0x0),
  fRefPulserQmean(0x0),
  fRefPulserOutlier(new AliTPCCalPad("RefPulserOutliers","RefPulserOutliers")),
  fRefPulserMasked(0x0),
  fRefCETmean(0x0),
  fRefCETrms(0x0),
  fRefCEQmean(0x0),
  fRefCEMasked(0x0),
  fRefALTROFPED(0x0),
  fRefALTROZsThr(0x0),
  fRefALTROAcqStart(0x0),
  fRefALTROAcqStop(0x0),
  fRefALTROMasked(0x0),
  fRefCalibRaw(0x0),
  fRefDataQA(0x0),
  fGoofieArray(0x0),
  fMapper(new AliTPCmapper(0x0)),
  fNpulserOutliers(-1),
  fIrocTimeOffset(0),
  fCETmaxLimitAbs(1.5),
  fPulTmaxLimitAbs(1.5),
  fPulQmaxLimitAbs(5),
  fPulQminLimit(11),
  fRuns(0),                         // run list with OCDB info
  fRunsStart(0),                    // start time for given run
  fRunsStop(0)                     // stop time for given run
{
  //
  // Default ctor
  //
}
//_____________________________________________________________________________________
AliTPCcalibDButil::~AliTPCcalibDButil()
{
  //
  // dtor
  //
  delete fPulserOutlier;
  delete fRefPulserOutlier;
  delete fMapper;
  if (fRefPadNoise) delete fRefPadNoise;
  if (fRefPedestals) delete fRefPedestals;
  if (fRefPedestalMasked) delete fRefPedestalMasked;
  if (fRefPulserTmean) delete fRefPulserTmean;
  if (fRefPulserTrms) delete fRefPulserTrms;
  if (fRefPulserQmean) delete fRefPulserQmean;
  if (fRefPulserMasked) delete fRefPulserMasked;
  if (fRefCETmean) delete fRefCETmean;
  if (fRefCETrms) delete fRefCETrms;
  if (fRefCEQmean) delete fRefCEQmean;
  if (fRefCEMasked) delete fRefCEMasked;
  if (fRefALTROFPED) delete fRefALTROFPED;
  if (fRefALTROZsThr) delete fRefALTROZsThr;
  if (fRefALTROAcqStart) delete fRefALTROAcqStart;
  if (fRefALTROAcqStop) delete fRefALTROAcqStop;
  if (fRefALTROMasked) delete fRefALTROMasked;
  if (fRefCalibRaw) delete fRefCalibRaw;
  if (fCurrentRefMap) delete fCurrentRefMap;    
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::UpdateFromCalibDB()
{
  //
  // Update pointers from calibDB
  //
  if (!fCalibDB) fCalibDB=AliTPCcalibDB::Instance();
  fCalibDB->UpdateNonRec();  // load all infromation now
  fPadNoise=fCalibDB->GetPadNoise();
  fPedestals=fCalibDB->GetPedestals();
  fPulserTmean=fCalibDB->GetPulserTmean();
  fPulserTrms=fCalibDB->GetPulserTrms();
  fPulserQmean=fCalibDB->GetPulserQmean();
  fCETmean=fCalibDB->GetCETmean();
  fCETrms=fCalibDB->GetCETrms();
  fCEQmean=fCalibDB->GetCEQmean();
  fALTROMasked=fCalibDB->GetALTROMasked();
  fGoofieArray=fCalibDB->GetGoofieSensors(fCalibDB->GetRun());
  fCalibRaw=fCalibDB->GetCalibRaw();
  fDataQA=fCalibDB->GetDataQA();
  UpdatePulserOutlierMap();
//   SetReferenceRun();
//   UpdateRefDataFromOCDB();
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::ProcessCEdata(const char* fitFormula, TVectorD &fitResultsA, TVectorD &fitResultsC,
                                      Int_t &noutliersCE, Double_t & chi2A, Double_t &chi2C, AliTPCCalPad * const outCE)
{
  //
  // Process the CE data for this run
  // the return TVectorD arrays contian the results of the fit
  // noutliersCE contains the number of pads marked as outliers,
  //   not including masked and edge pads
  //
  
  //retrieve CE and ALTRO data
  if (!fCETmean){
    TString fitString(fitFormula);
    fitString.ReplaceAll("++","#");
    Int_t ndim=fitString.CountChar('#')+2;
    fitResultsA.ResizeTo(ndim);
    fitResultsC.ResizeTo(ndim);
    fitResultsA.Zero();
    fitResultsC.Zero();
    noutliersCE=-1;
    return;
  }
  noutliersCE=0;
  //create outlier map
  AliTPCCalPad *out=0;
  if (outCE) out=outCE;
  else out=new AliTPCCalPad("outCE","outCE");
  AliTPCCalROC *rocMasked=0x0;
  //loop over all channels
  for (UInt_t iroc=0;iroc<fCETmean->kNsec;++iroc){
    AliTPCCalROC *rocData=fCETmean->GetCalROC(iroc);
    if (fALTROMasked) rocMasked=fALTROMasked->GetCalROC(iroc);
    AliTPCCalROC *rocOut=out->GetCalROC(iroc);
    if (!rocData) {
      noutliersCE+=AliTPCROC::Instance()->GetNChannels(iroc);
      rocOut->Add(1.);
      continue;
    }
    //add time offset to IROCs
    if (iroc<AliTPCROC::Instance()->GetNInnerSector())
      rocData->Add(fIrocTimeOffset);
    //select outliers
    UInt_t nrows=rocData->GetNrows();
    for (UInt_t irow=0;irow<nrows;++irow){
      UInt_t npads=rocData->GetNPads(irow);
      for (UInt_t ipad=0;ipad<npads;++ipad){
        rocOut->SetValue(irow,ipad,0);
        //exclude masked pads
        if (rocMasked && rocMasked->GetValue(irow,ipad)) {
          rocOut->SetValue(irow,ipad,1);
          continue;
        }
        //exclude first two rows in IROC and last two rows in OROC
        if (iroc<36){
          if (irow<2) rocOut->SetValue(irow,ipad,1);
        } else {
          if (irow>nrows-3) rocOut->SetValue(irow,ipad,1);
        }
        //exclude edge pads
        if (ipad==0||ipad==npads-1) rocOut->SetValue(irow,ipad,1);
        Float_t valTmean=rocData->GetValue(irow,ipad);
        //exclude values that are exactly 0
        if ( !(TMath::Abs(valTmean)>kAlmost0) ) {
          rocOut->SetValue(irow,ipad,1);
          ++noutliersCE;
        }
        // exclude channels with too large variations
        if (TMath::Abs(valTmean)>fCETmaxLimitAbs) {
          rocOut->SetValue(irow,ipad,1);
          ++noutliersCE;
        }
      }
    }
  }
  //perform fit
  TMatrixD dummy;
  Float_t chi2Af,chi2Cf;
  fCETmean->GlobalSidesFit(out,fitFormula,fitResultsA,fitResultsC,dummy,dummy,chi2Af,chi2Cf);
  chi2A=chi2Af;
  chi2C=chi2Cf;
  if (!outCE) delete out;
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::ProcessCEgraphs(TVectorD &vecTEntries, TVectorD &vecTMean, TVectorD &vecTRMS, TVectorD &vecTMedian,
                     TVectorD &vecQEntries, TVectorD &vecQMean, TVectorD &vecQRMS, TVectorD &vecQMedian,
                     Float_t &driftTimeA, Float_t &driftTimeC )
{
  //
  // Calculate statistical information from the CE graphs for drift time and charge
  //
  
  //reset arrays
  vecTEntries.ResizeTo(72);
  vecTMean.ResizeTo(72);
  vecTRMS.ResizeTo(72);
  vecTMedian.ResizeTo(72);
  vecQEntries.ResizeTo(72);
  vecQMean.ResizeTo(72);
  vecQRMS.ResizeTo(72);
  vecQMedian.ResizeTo(72);
  vecTEntries.Zero();
  vecTMean.Zero();
  vecTRMS.Zero();
  vecTMedian.Zero();
  vecQEntries.Zero();
  vecQMean.Zero();
  vecQRMS.Zero();
  vecQMedian.Zero();
  driftTimeA=0;
  driftTimeC=0;
  TObjArray *arrT=fCalibDB->GetCErocTtime();
  TObjArray *arrQ=fCalibDB->GetCErocQtime();
  if (arrT){
    for (Int_t isec=0;isec<74;++isec){
      TGraph *gr=(TGraph*)arrT->At(isec);
      if (!gr) continue;
      TVectorD values;
      Int_t npoints = gr->GetN();
      values.ResizeTo(npoints);
      Int_t nused =0;
      //skip first points, theres always some problems with finding the CE position
      for (Int_t ipoint=4; ipoint<npoints; ipoint++){
        if (gr->GetY()[ipoint]>500 && gr->GetY()[ipoint]<1020 ){
          values[nused]=gr->GetY()[ipoint];
          nused++;
        }
      }
      //
      if (isec<72) vecTEntries[isec]= nused;
      if (nused>1){
        if (isec<72){
          vecTMedian[isec] = TMath::Median(nused,values.GetMatrixArray());
          vecTMean[isec]   = TMath::Mean(nused,values.GetMatrixArray());
          vecTRMS[isec]    = TMath::RMS(nused,values.GetMatrixArray());
        } else if (isec==72){
          driftTimeA=TMath::Median(nused,values.GetMatrixArray());
        } else if (isec==73){
          driftTimeC=TMath::Median(nused,values.GetMatrixArray());
        }
      }
    }
  }
  if (arrQ){
    for (Int_t isec=0;isec<arrQ->GetEntriesFast();++isec){
      TGraph *gr=(TGraph*)arrQ->At(isec);
      if (!gr) continue;
      TVectorD values;
      Int_t npoints = gr->GetN();
      values.ResizeTo(npoints);
      Int_t nused =0;
      for (Int_t ipoint=0; ipoint<npoints; ipoint++){
        if (gr->GetY()[ipoint]>10 && gr->GetY()[ipoint]<500 ){
          values[nused]=gr->GetY()[ipoint];
          nused++;
        }
      }
      //
      vecQEntries[isec]= nused;
      if (nused>1){
        vecQMedian[isec] = TMath::Median(nused,values.GetMatrixArray());
        vecQMean[isec]   = TMath::Mean(nused,values.GetMatrixArray());
        vecQRMS[isec]    = TMath::RMS(nused,values.GetMatrixArray());
      }
    }
  }
}

//_____________________________________________________________________________________
void AliTPCcalibDButil::ProcessNoiseData(TVectorD &vNoiseMean, TVectorD &vNoiseMeanSenRegions,
                      TVectorD &vNoiseRMS, TVectorD &vNoiseRMSSenRegions,
                      Int_t &nonMaskedZero, Int_t &nNaN)
{
  //
  // process noise data
  // vNoiseMean/RMS contains the Mean/RMS noise of the complete TPC [0], IROCs only [1],
  //    OROCs small pads [2] and OROCs large pads [3]
  // vNoiseMean/RMSsenRegions constains the same information, but only for the sensitive regions (edge pads, corners, IROC spot)
  // nonMaskedZero contains the number of pads which show zero noise and were not masked. This might indicate an error
  //
  
  //set proper size and reset
  const UInt_t infoSize=4;
  vNoiseMean.ResizeTo(infoSize);
  vNoiseMeanSenRegions.ResizeTo(infoSize);
  vNoiseRMS.ResizeTo(infoSize);
  vNoiseRMSSenRegions.ResizeTo(infoSize);
  vNoiseMean.Zero();
  vNoiseMeanSenRegions.Zero();
  vNoiseRMS.Zero();
  vNoiseRMSSenRegions.Zero();
  nonMaskedZero=0;
  nNaN=0;
  //counters
  TVectorD c(infoSize);
  TVectorD cs(infoSize);
  //tpc parameters
  AliTPCParam par;
  par.Update();
  //retrieve noise and ALTRO data
  if (!fPadNoise) return;
  AliTPCCalROC *rocMasked=0x0;
  //create IROC, OROC1, OROC2 and sensitive region masks
  for (UInt_t isec=0;isec<AliTPCCalPad::kNsec;++isec){
    AliTPCCalROC *noiseROC=fPadNoise->GetCalROC(isec);
    if (fALTROMasked) rocMasked=fALTROMasked->GetCalROC(isec);
    UInt_t nrows=noiseROC->GetNrows();
    for (UInt_t irow=0;irow<nrows;++irow){
      UInt_t npads=noiseROC->GetNPads(irow);
      for (UInt_t ipad=0;ipad<npads;++ipad){
        //don't use masked channels;
        if (rocMasked && rocMasked->GetValue(irow,ipad)) continue;
        Float_t noiseVal=noiseROC->GetValue(irow,ipad);
        //check if noise==0
        if (noiseVal<kAlmost0) {
          ++nonMaskedZero;
          continue;
        }
        //check for nan
        if ( !(noiseVal<10000000) ){
	  AliInfo(Form("Warning: nan detected in (sec,row,pad - val): %02d,%02d,%03d - %.1f\n",isec,irow,ipad,noiseVal));
          ++nNaN;
          continue;
        }
        Int_t cpad=(Int_t)ipad-(Int_t)npads/2;
        Int_t masksen=1; // sensitive pards are not masked (0)
        if (ipad<2||npads-ipad-1<2) masksen=0; //don't mask edge pads (sensitive)
        if (isec<AliTPCROC::Instance()->GetNInnerSector()){
          //IROCs
          if (irow>19&&irow<46){
            if (TMath::Abs(cpad)<7) masksen=0; //IROC spot
          }
          Int_t type=1;
          vNoiseMean[type]+=noiseVal;
          vNoiseRMS[type]+=noiseVal*noiseVal;
          ++c[type];
          if (!masksen){
            vNoiseMeanSenRegions[type]+=noiseVal;
            vNoiseRMSSenRegions[type]+=noiseVal*noiseVal;
            ++cs[type];
          }
        } else {
          //OROCs
          //define sensive regions
          if ((nrows-irow-1)<3) masksen=0; //last three rows in OROCs are sensitive
          if ( irow>75 ){
            Int_t padEdge=(Int_t)TMath::Min(ipad,npads-ipad);
            if (padEdge<((((Int_t)irow-76)/4+1))*2) masksen=0; //OROC outer corners are sensitive
          }
          if ((Int_t)irow<par.GetNRowUp1()){
            //OROC1
            Int_t type=2;
            vNoiseMean[type]+=noiseVal;
            vNoiseRMS[type]+=noiseVal*noiseVal;
            ++c[type];
            if (!masksen){
              vNoiseMeanSenRegions[type]+=noiseVal;
              vNoiseRMSSenRegions[type]+=noiseVal*noiseVal;
              ++cs[type];
            }
          }else{
            //OROC2
            Int_t type=3;
            vNoiseMean[type]+=noiseVal;
            vNoiseRMS[type]+=noiseVal*noiseVal;
            ++c[type];
            if (!masksen){
              vNoiseMeanSenRegions[type]+=noiseVal;
              vNoiseRMSSenRegions[type]+=noiseVal*noiseVal;
              ++cs[type];
            }
          }
        }
        //whole tpc
        Int_t type=0;
        vNoiseMean[type]+=noiseVal;
        vNoiseRMS[type]+=noiseVal*noiseVal;
        ++c[type];
        if (!masksen){
          vNoiseMeanSenRegions[type]+=noiseVal;
          vNoiseRMSSenRegions[type]+=noiseVal*noiseVal;
          ++cs[type];
        }
      }//end loop pads
    }//end loop rows
  }//end loop sectors (rocs)
  
  //calculate mean and RMS
  const Double_t verySmall=0.0000000001;
  for (UInt_t i=0;i<infoSize;++i){
    Double_t mean=0;
    Double_t rms=0;
    Double_t meanSen=0;
    Double_t rmsSen=0;
    
    if (c[i]>verySmall){
      AliInfo(Form("i: %d - m: %.3f, c: %.0f, r: %.3f\n",i,vNoiseMean[i],c[i],vNoiseRMS[i]));
      mean=vNoiseMean[i]/c[i];
      rms=vNoiseRMS[i];
      rms=TMath::Sqrt(TMath::Abs(rms/c[i]-mean*mean));
    }
    vNoiseMean[i]=mean;
    vNoiseRMS[i]=rms;
    
    if (cs[i]>verySmall){
      meanSen=vNoiseMeanSenRegions[i]/cs[i];
      rmsSen=vNoiseRMSSenRegions[i];
      rmsSen=TMath::Sqrt(TMath::Abs(rmsSen/cs[i]-meanSen*meanSen));
    }
    vNoiseMeanSenRegions[i]=meanSen;
    vNoiseRMSSenRegions[i]=rmsSen;
  }
}

//_____________________________________________________________________________________
void AliTPCcalibDButil::ProcessQAData(TVectorD &vQaOcc, TVectorD &vQaQtot, 
				      TVectorD &vQaQmax)
{
  //
  // process QA data
  //
  // vQaOcc/Qtot/Qmax contains the Mean occupancy/Qtot/Qmax for each sector
  //


  const UInt_t infoSize = 72;
  //reset counters to error number
  vQaOcc.ResizeTo(infoSize);
  vQaOcc.Zero();
  vQaQtot.ResizeTo(infoSize);
  vQaQtot.Zero();
  vQaQmax.ResizeTo(infoSize);
  vQaQmax.Zero();
  //counter
  //retrieve pulser and ALTRO data
  
  if (!fDataQA) {
    
    AliInfo("No QA data");
    return;
  }
  if (fDataQA->GetEventCounter()<=0) {

    AliInfo("No QA data");
    return; // no data processed
  }
  //
  fDataQA->Analyse();

  TVectorD normOcc(infoSize);
  TVectorD normQ(infoSize);

  for (UInt_t isec=0;isec<AliTPCCalPad::kNsec;++isec){

    AliInfo(Form("Sector %d\n", isec));
    AliTPCCalROC* occupancyROC = fDataQA->GetNoThreshold()->GetCalROC(isec); 
    AliTPCCalROC* nclusterROC = fDataQA->GetNLocalMaxima()->GetCalROC(isec); 
    AliTPCCalROC* qROC = fDataQA->GetMeanCharge()->GetCalROC(isec); 
    AliTPCCalROC* qmaxROC = fDataQA->GetMaxCharge()->GetCalROC(isec); 
    if (!occupancyROC) continue;
    if (!nclusterROC) continue;
    if (!qROC) continue;
    if (!qmaxROC) continue;
    
    const UInt_t nchannels=occupancyROC->GetNchannels();

    AliInfo(Form("Nchannels %d\n", nchannels));

    for (UInt_t ichannel=0;ichannel<nchannels;++ichannel){

      vQaOcc[isec] += occupancyROC->GetValue(ichannel);
      ++normOcc[isec];

      Float_t nClusters = nclusterROC->GetValue(ichannel);
      normQ[isec] += nClusters;
      vQaQtot[isec]+=nClusters*qROC->GetValue(ichannel);
      vQaQmax[isec]+=nClusters*qmaxROC->GetValue(ichannel);
    }
  }

  //calculate mean values
  for (UInt_t isec=0;isec<AliTPCCalPad::kNsec;++isec){
    
    if (normOcc[isec]>0) vQaOcc[isec] /= normOcc[isec];
    else vQaOcc[isec] = 0;

    if (normQ[isec]>0) {
      vQaQtot[isec] /= normQ[isec];
      vQaQmax[isec] /= normQ[isec];
    }else {

      vQaQtot[isec] = 0;
      vQaQmax[isec] = 0;
    }
  }
}

//_____________________________________________________________________________________
void AliTPCcalibDButil::ProcessPulser(TVectorD &vMeanTime)
{
  //
  // Process the Pulser information
  // vMeanTime:     pulser mean time position in IROC-A, IROC-C, OROC-A, OROC-C
  //

  const UInt_t infoSize=4;
  //reset counters to error number
  vMeanTime.ResizeTo(infoSize);
  vMeanTime.Zero();
  //counter
  TVectorD c(infoSize);
  //retrieve pulser and ALTRO data
  if (!fPulserTmean) return;
  //
  //get Outliers
  AliTPCCalROC *rocOut=0x0;
  for (UInt_t isec=0;isec<AliTPCCalPad::kNsec;++isec){
    AliTPCCalROC *tmeanROC=fPulserTmean->GetCalROC(isec);
    if (!tmeanROC) continue;
    rocOut=fPulserOutlier->GetCalROC(isec);
    UInt_t nchannels=tmeanROC->GetNchannels();
    for (UInt_t ichannel=0;ichannel<nchannels;++ichannel){
      if (rocOut && rocOut->GetValue(ichannel)) continue;
      Float_t val=tmeanROC->GetValue(ichannel);
      Int_t type=isec/18;
      vMeanTime[type]+=val;
      ++c[type];
    }
  }
  //calculate mean
  for (UInt_t itype=0; itype<infoSize; ++itype){
    if (c[itype]>0) vMeanTime[itype]/=c[itype];
    else vMeanTime[itype]=0;
  }
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::ProcessALTROConfig(Int_t &nMasked)
{
  //
  // Get Values from ALTRO configuration data
  //
  nMasked=-1;
  if (!fALTROMasked) return;
  nMasked=0;
  for (Int_t isec=0;isec<fALTROMasked->kNsec; ++isec){
    AliTPCCalROC *rocMasked=fALTROMasked->GetCalROC(isec);
    for (UInt_t ichannel=0; ichannel<rocMasked->GetNchannels();++ichannel){
      if (rocMasked->GetValue(ichannel)) ++nMasked;
    }
  }
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::ProcessGoofie(TVectorD & vecEntries, TVectorD & vecMedian, TVectorD &vecMean, TVectorD &vecRMS)
{
  //
  // Proces Goofie values, return statistical information of the currently set goofieArray
  // The meaning of the entries are given below
  /*
  1       TPC_ANODE_I_A00_STAT
  2       TPC_DVM_CO2
  3       TPC_DVM_DriftVelocity
  4       TPC_DVM_FCageHV
  5       TPC_DVM_GainFar
  6       TPC_DVM_GainNear
  7       TPC_DVM_N2
  8       TPC_DVM_NumberOfSparks
  9       TPC_DVM_PeakAreaFar
  10      TPC_DVM_PeakAreaNear
  11      TPC_DVM_PeakPosFar
  12      TPC_DVM_PeakPosNear
  13      TPC_DVM_PickupHV
  14      TPC_DVM_Pressure
  15      TPC_DVM_T1_Over_P
  16      TPC_DVM_T2_Over_P
  17      TPC_DVM_T_Over_P
  18      TPC_DVM_TemperatureS1
   */
  if (!fGoofieArray){
    Int_t nsensors=19;
    vecEntries.ResizeTo(nsensors);
    vecMedian.ResizeTo(nsensors);
    vecMean.ResizeTo(nsensors);
    vecRMS.ResizeTo(nsensors);
    vecEntries.Zero();
    vecMedian.Zero();
    vecMean.Zero();
    vecRMS.Zero();
    return;
  }
  Double_t kEpsilon=0.0000000001;
  Double_t kBig=100000000000.;
  Int_t nsensors = fGoofieArray->NumSensors();
  vecEntries.ResizeTo(nsensors);
  vecMedian.ResizeTo(nsensors);
  vecMean.ResizeTo(nsensors);
  vecRMS.ResizeTo(nsensors);
  TVectorF values;
  for (Int_t isensor=0; isensor<fGoofieArray->NumSensors();isensor++){
    AliDCSSensor *gsensor = fGoofieArray->GetSensor(isensor);
    if (gsensor &&  gsensor->GetGraph()){
      Int_t npoints = gsensor->GetGraph()->GetN();
      // filter zeroes
      values.ResizeTo(npoints);
      Int_t nused =0;
      for (Int_t ipoint=0; ipoint<npoints; ipoint++){
        if (TMath::Abs(gsensor->GetGraph()->GetY()[ipoint])>kEpsilon &&
            TMath::Abs(gsensor->GetGraph()->GetY()[ipoint])<kBig ){
              values[nused]=gsensor->GetGraph()->GetY()[ipoint];
              nused++;
            }
      }
      //
      vecEntries[isensor]= nused;
      if (nused>1){
        vecMedian[isensor] = TMath::Median(nused,values.GetMatrixArray());
        vecMean[isensor]   = TMath::Mean(nused,values.GetMatrixArray());
        vecRMS[isensor]    = TMath::RMS(nused,values.GetMatrixArray());
      }
    }
  }
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::ProcessPedestalVariations(TVectorF &pedestalDeviations)
{
  //
  // check the variations of the pedestal data to the reference pedestal data
  // thresholds are 0.5, 1.0, 1.5 and 2 timebins respectively.
  //
  const Int_t npar=4;
  TVectorF vThres(npar); //thresholds
  Int_t nActive=0;       //number of active channels
  
  //reset and set thresholds
  pedestalDeviations.ResizeTo(npar);
  for (Int_t i=0;i<npar;++i){
    pedestalDeviations.GetMatrixArray()[i]=0;
    vThres.GetMatrixArray()[i]=(i+1)*.5;
  }
  //check all needed data is available
  if (!fRefPedestals || !fPedestals || !fALTROMasked || !fRefALTROMasked) return;
  //loop over all channels
  for (UInt_t isec=0;isec<AliTPCCalPad::kNsec;++isec){
    AliTPCCalROC *pROC=fPedestals->GetCalROC(isec);
    AliTPCCalROC *pRefROC=fRefPedestals->GetCalROC(isec);
    AliTPCCalROC *mROC=fALTROMasked->GetCalROC(isec);
    AliTPCCalROC *mRefROC=fRefALTROMasked->GetCalROC(isec);
    UInt_t nrows=mROC->GetNrows();
    for (UInt_t irow=0;irow<nrows;++irow){
      UInt_t npads=mROC->GetNPads(irow);
      for (UInt_t ipad=0;ipad<npads;++ipad){
        //don't use masked channels;
        if (mROC   ->GetValue(irow,ipad)) continue;
        if (mRefROC->GetValue(irow,ipad)) continue;
        Float_t deviation=TMath::Abs(pROC->GetValue(irow,ipad)-pRefROC->GetValue(irow,ipad));
        for (Int_t i=0;i<npar;++i){
          if (deviation>vThres[i])
            ++pedestalDeviations.GetMatrixArray()[i];
        }
        ++nActive;
      }//end ipad
    }//ind irow
  }//end isec
  if (nActive>0){
    for (Int_t i=0;i<npar;++i){
      pedestalDeviations.GetMatrixArray()[i]/=nActive;
    }
  }
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::ProcessNoiseVariations(TVectorF &noiseDeviations)
{
  //
  // check the variations of the noise data to the reference noise data
  // thresholds are 5, 10, 15 and 20 percent respectively.
  //
  const Int_t npar=4;
  TVectorF vThres(npar); //thresholds
  Int_t nActive=0;       //number of active channels
  
  //reset and set thresholds
  noiseDeviations.ResizeTo(npar);
  for (Int_t i=0;i<npar;++i){
    noiseDeviations.GetMatrixArray()[i]=0;
    vThres.GetMatrixArray()[i]=(i+1)*.05;
  }
  //check all needed data is available
  if (!fRefPadNoise || !fPadNoise || !fALTROMasked || !fRefALTROMasked) return;
  //loop over all channels
  for (UInt_t isec=0;isec<AliTPCCalPad::kNsec;++isec){
    AliTPCCalROC *nROC=fPadNoise->GetCalROC(isec);
    AliTPCCalROC *nRefROC=fRefPadNoise->GetCalROC(isec);
    AliTPCCalROC *mROC=fALTROMasked->GetCalROC(isec);
    AliTPCCalROC *mRefROC=fRefALTROMasked->GetCalROC(isec);
    UInt_t nrows=mROC->GetNrows();
    for (UInt_t irow=0;irow<nrows;++irow){
      UInt_t npads=mROC->GetNPads(irow);
      for (UInt_t ipad=0;ipad<npads;++ipad){
        //don't use masked channels;
        if (mROC   ->GetValue(irow,ipad)) continue;
        if (mRefROC->GetValue(irow,ipad)) continue;
        if (nRefROC->GetValue(irow,ipad)==0) continue;
        Float_t deviation=(nROC->GetValue(irow,ipad)/nRefROC->GetValue(irow,ipad))-1;
        for (Int_t i=0;i<npar;++i){
          if (deviation>vThres[i])
            ++noiseDeviations.GetMatrixArray()[i];
        }
        ++nActive;
      }//end ipad
    }//ind irow
  }//end isec
  if (nActive>0){
    for (Int_t i=0;i<npar;++i){
      noiseDeviations.GetMatrixArray()[i]/=nActive;
    }
  }
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::ProcessPulserVariations(TVectorF &pulserQdeviations, Float_t &varQMean,
                                                Int_t &npadsOutOneTB, Int_t &npadsOffAdd)
{
  //
  // check the variations of the pulserQmean data to the reference pulserQmean data: pulserQdeviations
  // thresholds are .5, 1, 5 and 10 percent respectively.
  // 
  //
  const Int_t npar=4;
  TVectorF vThres(npar); //thresholds
  Int_t nActive=0;       //number of active channels
  
  //reset and set thresholds
  pulserQdeviations.ResizeTo(npar);
  for (Int_t i=0;i<npar;++i){
    pulserQdeviations.GetMatrixArray()[i]=0;
  }
  npadsOutOneTB=0;
  npadsOffAdd=0;
  varQMean=0;
  vThres.GetMatrixArray()[0]=.005;
  vThres.GetMatrixArray()[1]=.01;
  vThres.GetMatrixArray()[2]=.05;
  vThres.GetMatrixArray()[3]=.1;
  //check all needed data is available
  if (!fRefPulserTmean || !fPulserTmean || !fPulserQmean || !fRefPulserQmean || !fALTROMasked || !fRefALTROMasked) return;
  //
  UpdateRefPulserOutlierMap();
  //loop over all channels
  for (UInt_t isec=0;isec<AliTPCCalPad::kNsec;++isec){
    AliTPCCalROC *pqROC=fPulserQmean->GetCalROC(isec);
    AliTPCCalROC *pqRefROC=fRefPulserQmean->GetCalROC(isec);
    AliTPCCalROC *ptROC=fPulserTmean->GetCalROC(isec);
//     AliTPCCalROC *ptRefROC=fRefPulserTmean->GetCalROC(isec);
    AliTPCCalROC *mROC=fALTROMasked->GetCalROC(isec);
    AliTPCCalROC *mRefROC=fRefALTROMasked->GetCalROC(isec);
    AliTPCCalROC *oROC=fPulserOutlier->GetCalROC(isec);
    Float_t ptmean=ptROC->GetMean(oROC);
    UInt_t nrows=mROC->GetNrows();
    for (UInt_t irow=0;irow<nrows;++irow){
      UInt_t npads=mROC->GetNPads(irow);
      for (UInt_t ipad=0;ipad<npads;++ipad){
        //don't use masked channels;
        if (mROC   ->GetValue(irow,ipad)) continue;
        if (mRefROC->GetValue(irow,ipad)) continue;
        //don't user edge pads
        if (ipad==0||ipad==npads-1) continue;
        //data
        Float_t pq=pqROC->GetValue(irow,ipad);
        Float_t pqRef=pqRefROC->GetValue(irow,ipad);
        Float_t pt=ptROC->GetValue(irow,ipad);
//         Float_t ptRef=ptRefROC->GetValue(irow,ipad);
        //comparisons q
        Float_t deviation=TMath::Abs(pq/pqRef-1);
        for (Int_t i=0;i<npar;++i){
          if (deviation>vThres[i])
            ++pulserQdeviations.GetMatrixArray()[i];
        }
        if (pqRef>11&&pq<11) ++npadsOffAdd;
        varQMean+=pq-pqRef;
        //comparisons t
        if (TMath::Abs(pt-ptmean)>1) ++npadsOutOneTB;
        ++nActive;
      }//end ipad
    }//ind irow
  }//end isec
  if (nActive>0){
    for (Int_t i=0;i<npar;++i){
      pulserQdeviations.GetMatrixArray()[i]/=nActive;
      varQMean/=nActive;
    }
  }
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::UpdatePulserOutlierMap()
{
  //
  // Update the outlier map of the pulser data
  //
  PulserOutlierMap(fPulserOutlier,fPulserTmean, fPulserQmean);
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::UpdateRefPulserOutlierMap()
{
  //
  // Update the outlier map of the pulser reference data
  //
  PulserOutlierMap(fRefPulserOutlier,fRefPulserTmean, fRefPulserQmean);
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::PulserOutlierMap(AliTPCCalPad *pulOut, const AliTPCCalPad *pulT, const AliTPCCalPad *pulQ)
{
  //
  // Create a map that contains outliers from the Pulser calibration data.
  // The outliers include masked channels, edge pads and pads with
  //   too large timing and charge variations.
  // fNpulserOutliers is the number of outliers in the Pulser calibration data.
  //   those do not contain masked and edge pads
  //
  if (!pulT||!pulQ) {
    //reset map
    pulOut->Multiply(0.);
    fNpulserOutliers=-1;
    return;
  }
  AliTPCCalROC *rocMasked=0x0;
  fNpulserOutliers=0;
  
  //Create Outlier Map
  for (UInt_t isec=0;isec<AliTPCCalPad::kNsec;++isec){
    AliTPCCalROC *tmeanROC=pulT->GetCalROC(isec);
    AliTPCCalROC *qmeanROC=pulQ->GetCalROC(isec);
    AliTPCCalROC *outROC=pulOut->GetCalROC(isec);
    if (!tmeanROC||!qmeanROC) {
      //reset outliers in this ROC
      outROC->Multiply(0.);
      continue;
    }
    if (fALTROMasked) rocMasked=fALTROMasked->GetCalROC(isec);
//     Double_t dummy=0;
//     Float_t qmedian=qmeanROC->GetLTM(&dummy,.5);
//     Float_t tmedian=tmeanROC->GetLTM(&dummy,.5);
    UInt_t nrows=tmeanROC->GetNrows();
    for (UInt_t irow=0;irow<nrows;++irow){
      UInt_t npads=tmeanROC->GetNPads(irow);
      for (UInt_t ipad=0;ipad<npads;++ipad){
        Int_t outlier=0,masked=0;
        Float_t q=qmeanROC->GetValue(irow,ipad);
        Float_t t=tmeanROC->GetValue(irow,ipad);
        //masked channels are outliers
        if (rocMasked && rocMasked->GetValue(irow,ipad)) masked=1;
        //edge pads are outliers
        if (ipad==0||ipad==npads-1) masked=1;
        //channels with too large charge or timing deviation from the meadian are outliers
//         if (TMath::Abs(q-qmedian)>fPulQmaxLimitAbs || TMath::Abs(t-tmedian)>fPulTmaxLimitAbs) outlier=1;
        if (q<fPulQminLimit && !masked) outlier=1;
        //check for nan
        if ( !(q<10000000) || !(t<10000000)) outlier=1;
        outROC->SetValue(irow,ipad,outlier+masked);
        fNpulserOutliers+=outlier;
      }
    }
  }
}
//_____________________________________________________________________________________
AliTPCCalPad* AliTPCcalibDButil::CreatePadTime0(Int_t model, Double_t &gyA, Double_t &gyC, Double_t &chi2A, Double_t &chi2C )
{
  //
  // Create pad time0 object from pulser and/or CE data, depending on the selected model
  // Model 0: normalise each readout chamber to its mean, outlier cutted, only Pulser
  // Model 1: normalise IROCs/OROCs of each readout side to its mean, only Pulser
  // Model 2: use CE data and a combination CE fit + pulser in the outlier regions.
  //
  // In case model 2 is invoked - gy arival time gradient is also returned
  //
  gyA=0;
  gyC=0;
  AliTPCCalPad *padTime0=new AliTPCCalPad("PadTime0",Form("PadTime0-Model_%d",model));
  // decide between different models
  if (model==0||model==1){
    TVectorD vMean;
    if (model==1) ProcessPulser(vMean);
    for (UInt_t isec=0;isec<AliTPCCalPad::kNsec;++isec){
      AliTPCCalROC *rocPulTmean=fPulserTmean->GetCalROC(isec);
      if (!rocPulTmean) continue;
      AliTPCCalROC *rocTime0=padTime0->GetCalROC(isec);
      AliTPCCalROC *rocOut=fPulserOutlier->GetCalROC(isec);
      Float_t mean=rocPulTmean->GetMean(rocOut);
      //treat case where a whole partition is masked
      if ( TMath::Abs(mean)<kAlmost0 ) mean=rocPulTmean->GetMean();
      if (model==1) {
        Int_t type=isec/18;
        mean=vMean[type];
      }
      UInt_t nrows=rocTime0->GetNrows();
      for (UInt_t irow=0;irow<nrows;++irow){
        UInt_t npads=rocTime0->GetNPads(irow);
        for (UInt_t ipad=0;ipad<npads;++ipad){
          Float_t time=rocPulTmean->GetValue(irow,ipad);
          //in case of an outlier pad use the mean of the altro values.
          //This should be the most precise guess in that case.
          if (rocOut->GetValue(irow,ipad)) {
            time=GetMeanAltro(rocPulTmean,irow,ipad,rocOut);
            if ( TMath::Abs(time)<kAlmost0 ) time=mean;
          }
          Float_t val=time-mean;
          rocTime0->SetValue(irow,ipad,val);
        }
      }
    }
  } else if (model==2){  
    Double_t pgya,pgyc,pchi2a,pchi2c;
    AliTPCCalPad * padPulser = CreatePadTime0(1,pgya,pgyc,pchi2a,pchi2c);
    fCETmean->Add(padPulser,-1.);
    TVectorD vA,vC;
    AliTPCCalPad outCE("outCE","outCE");
    Int_t nOut;
    ProcessCEdata("(sector<36)++gy++gx++(lx-134)++(sector<36)*(lx-134)++(ly/lx)^2",vA,vC,nOut,chi2A, chi2C,&outCE);
    AliTPCCalPad *padFit=AliTPCCalPad::CreateCalPadFit("1++0++gy++0++(lx-134)++0++0",vA,vC);
//     AliTPCCalPad *padFit=AliTPCCalPad::CreateCalPadFit("1++(sector<36)++gy++gx++(lx-134)++(sector<36)*(lx-134)",vA,vC);
    if (!padFit) { delete padPulser; return 0;}
    gyA=vA[2];
    gyC=vC[2];
    fCETmean->Add(padPulser,1.);
    padTime0->Add(fCETmean);
    padTime0->Add(padFit,-1);  
    delete padPulser;
    TVectorD vFitROC;
    TMatrixD mFitROC;
    Float_t chi2;
    for (UInt_t isec=0;isec<AliTPCCalPad::kNsec;++isec){
      AliTPCCalROC *rocPulTmean=fPulserTmean->GetCalROC(isec);
      AliTPCCalROC *rocTime0=padTime0->GetCalROC(isec);
      AliTPCCalROC *rocOutPul=fPulserOutlier->GetCalROC(isec);
      AliTPCCalROC *rocOutCE=outCE.GetCalROC(isec);
      rocTime0->GlobalFit(rocOutCE,kFALSE,vFitROC,mFitROC,chi2);
      AliTPCCalROC *rocCEfit=AliTPCCalROC::CreateGlobalFitCalROC(vFitROC, isec);
      Float_t mean=rocPulTmean->GetMean(rocOutPul);
      if ( TMath::Abs(mean)<kAlmost0 ) mean=rocPulTmean->GetMean();
      UInt_t nrows=rocTime0->GetNrows();
      for (UInt_t irow=0;irow<nrows;++irow){
        UInt_t npads=rocTime0->GetNPads(irow);
        for (UInt_t ipad=0;ipad<npads;++ipad){
          Float_t timePulser=rocPulTmean->GetValue(irow,ipad)-mean;
          if (rocOutCE->GetValue(irow,ipad)){
            Float_t valOut=rocCEfit->GetValue(irow,ipad);
            if (!rocOutPul->GetValue(irow,ipad)) valOut+=timePulser;
            rocTime0->SetValue(irow,ipad,valOut);
          }
        }
      }
      delete rocCEfit;
    }
    delete padFit;
  }
  Double_t median = padTime0->GetMedian();
  padTime0->Add(-median);  // normalize to median
  return padTime0;
}
//_____________________________________________________________________________________
Float_t AliTPCcalibDButil::GetMeanAltro(const AliTPCCalROC *roc, const Int_t row, const Int_t pad, AliTPCCalROC *const rocOut)
{
  //
  // GetMeanAlto information
  //
  if (roc==0) return 0.;
  const Int_t sector=roc->GetSector();
  AliTPCROC *tpcRoc=AliTPCROC::Instance();
  const UInt_t altroRoc=fMapper->GetFEC(sector,row,pad)*8+fMapper->GetChip(sector,row,pad);
  Float_t mean=0;
  Int_t   n=0;
  
  //loop over a small range around the requested pad (+-10 rows/pads)
  for (Int_t irow=row-10;irow<row+10;++irow){
    if (irow<0||irow>(Int_t)tpcRoc->GetNRows(sector)-1) continue;
    for (Int_t ipad=pad-10; ipad<pad+10;++ipad){
      if (ipad<0||ipad>(Int_t)tpcRoc->GetNPads(sector,irow)-1) continue;
      const UInt_t altroCurr=fMapper->GetFEC(sector,irow,ipad)*8+fMapper->GetChip(sector,irow,ipad);
      if (altroRoc!=altroCurr) continue;
      if ( rocOut && rocOut->GetValue(irow,ipad) ) continue;
      Float_t val=roc->GetValue(irow,ipad);
      mean+=val;
      ++n;
    }
  }
  if (n>0) mean/=n;
  return mean;
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::SetRefFile(const char* filename)
{
  //
  // load cal pad objects form the reference file
  //
  TDirectory *currDir=gDirectory;
  TFile f(filename);
  fRefPedestals=(AliTPCCalPad*)f.Get("Pedestals");
  fRefPadNoise=(AliTPCCalPad*)f.Get("PadNoise");
  //pulser data
  fRefPulserTmean=(AliTPCCalPad*)f.Get("PulserTmean");
  fRefPulserTrms=(AliTPCCalPad*)f.Get("PulserTrms");
  fRefPulserQmean=(AliTPCCalPad*)f.Get("PulserQmean");
  //CE data
  fRefCETmean=(AliTPCCalPad*)f.Get("CETmean");
  fRefCETrms=(AliTPCCalPad*)f.Get("CETrms");
  fRefCEQmean=(AliTPCCalPad*)f.Get("CEQmean");
  //Altro data
//   fRefALTROAcqStart=(AliTPCCalPad*)f.Get("ALTROAcqStart");
//   fRefALTROZsThr=(AliTPCCalPad*)f.Get("ALTROZsThr");
//   fRefALTROFPED=(AliTPCCalPad*)f.Get("ALTROFPED");
//   fRefALTROAcqStop=(AliTPCCalPad*)f.Get("ALTROAcqStop");
  fRefALTROMasked=(AliTPCCalPad*)f.Get("ALTROMasked");
  f.Close();
  currDir->cd();
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::UpdateRefDataFromOCDB()
{
  //
  // set reference data from OCDB Reference map
  //
  if (!fRefMap) {
    AliWarning("Referenc map not set!");
    return;
  }
  
  TString cdbPath;
  AliCDBEntry* entry = 0x0;
  Bool_t hasAnyChanged=kFALSE;

  //pedestals
  cdbPath="TPC/Calib/Pedestals";
  if (HasRefChanged(cdbPath.Data())){
    hasAnyChanged=kTRUE;
    //delete old entries
    if (fRefPedestals) delete fRefPedestals;
    if (fRefPedestalMasked) delete fRefPedestalMasked;
    fRefPedestals=fRefPedestalMasked=0x0;
    //get new entries
    entry=GetRefEntry(cdbPath.Data());
    if (entry){
      entry->SetOwner(kTRUE);
      fRefPedestals=GetRefCalPad(entry);
      delete entry;
      fRefPedestalMasked=GetAltroMasked(cdbPath, "MaskedPedestals");
    }
  }

  //noise
  cdbPath="TPC/Calib/PadNoise";
  if (HasRefChanged(cdbPath.Data())){
    hasAnyChanged=kTRUE;
    //delete old entry
    if (fRefPadNoise) delete fRefPadNoise;
    fRefPadNoise=0x0;
    //get new entry
    entry=GetRefEntry(cdbPath.Data());
    if (entry){
      entry->SetOwner(kTRUE);
      fRefPadNoise=GetRefCalPad(entry);
      delete entry;
    }
  }
  
  //pulser
  cdbPath="TPC/Calib/Pulser";
  if (HasRefChanged(cdbPath.Data())){
    hasAnyChanged=kTRUE;
    //delete old entries
    if (fRefPulserTmean) delete fRefPulserTmean;
    if (fRefPulserTrms) delete fRefPulserTrms;
    if (fRefPulserQmean) delete fRefPulserQmean;
    if (fRefPulserMasked) delete fRefPulserMasked;
    fRefPulserTmean=fRefPulserTrms=fRefPulserQmean=fRefPulserMasked=0x0;
    //get new entries
    entry=GetRefEntry(cdbPath.Data());
    if (entry){
      entry->SetOwner(kTRUE);
      fRefPulserTmean=GetRefCalPad(entry,"PulserTmean");
      fRefPulserTrms=GetRefCalPad(entry,"PulserTrms");
      fRefPulserQmean=GetRefCalPad(entry,"PulserQmean");
      delete entry;
      fRefPulserMasked=GetAltroMasked(cdbPath, "MaskedPulser");
    }
  }

  //ce
  cdbPath="TPC/Calib/CE";
  if (HasRefChanged(cdbPath.Data())){
    hasAnyChanged=kTRUE;
    //delete old entries
    if (fRefCETmean) delete fRefCETmean;
    if (fRefCETrms) delete fRefCETrms;
    if (fRefCEQmean) delete fRefCEQmean;
    if (fRefCEMasked) delete fRefCEMasked;
    fRefCETmean=fRefCETrms=fRefCEQmean=fRefCEMasked=0x0;
    //get new entries
    entry=GetRefEntry(cdbPath.Data());
    if (entry){
      entry->SetOwner(kTRUE);
      fRefCETmean=GetRefCalPad(entry,"CETmean");
      fRefCETrms=GetRefCalPad(entry,"CETrms");
      fRefCEQmean=GetRefCalPad(entry,"CEQmean");
      delete entry;
      fRefCEMasked=GetAltroMasked(cdbPath, "MaskedCE");
    }
  }
  
  //altro data
  cdbPath="TPC/Calib/AltroConfig";
  if (HasRefChanged(cdbPath.Data())){
    hasAnyChanged=kTRUE;
    //delete old entries
    if (fRefALTROFPED) delete fRefALTROFPED;
    if (fRefALTROZsThr) delete fRefALTROZsThr;
    if (fRefALTROAcqStart) delete fRefALTROAcqStart;
    if (fRefALTROAcqStop) delete fRefALTROAcqStop;
    if (fRefALTROMasked) delete fRefALTROMasked;
    fRefALTROFPED=fRefALTROZsThr=fRefALTROAcqStart=fRefALTROAcqStop=fRefALTROMasked=0x0;
    //get new entries
    entry=GetRefEntry(cdbPath.Data());
    if (entry){
      entry->SetOwner(kTRUE);
      fRefALTROFPED=GetRefCalPad(entry,"FPED");
      fRefALTROZsThr=GetRefCalPad(entry,"ZsThr");
      fRefALTROAcqStart=GetRefCalPad(entry,"AcqStart");
      fRefALTROAcqStop=GetRefCalPad(entry,"AcqStop");
      fRefALTROMasked=GetRefCalPad(entry,"Masked");
      delete entry;
    }
  }
  
  //raw data
  /*
  cdbPath="TPC/Calib/Raw";
  if (HasRefChanged(cdbPath.Data())){
    hasAnyChanged=kTRUE;
    //delete old entry
    if (fRefCalibRaw) delete fRefCalibRaw;
    //get new entry
    entry=GetRefEntry(cdbPath.Data());
    if (entry){
      entry->SetOwner(kTRUE);
      TObjArray *arr=(TObjArray*)entry->GetObject();
      if (!arr){
        AliError(Form("Could not get object from entry '%s'\nPlease check!!!",entry->GetId().GetPath().Data()));
      } else {
        fRefCalibRaw=(AliTPCCalibRaw*)arr->At(0)->Clone();
      }
    }
  }
  */

  //data qa
  cdbPath="TPC/Calib/QA";
  if (HasRefChanged(cdbPath.Data())){
    hasAnyChanged=kTRUE;
    //delete old entry
    if (fRefDataQA) delete fRefDataQA;
    //get new entry
    entry=GetRefEntry(cdbPath.Data());
    if (entry){
      entry->SetOwner(kTRUE);
      fRefDataQA=dynamic_cast<AliTPCdataQA*>(entry->GetObject());
      if (!fRefDataQA){
        AliError(Form("Could not get object from entry '%s'\nPlease check!!!",entry->GetId().GetPath().Data()));
      } else {
        fRefDataQA=(AliTPCdataQA*)fRefDataQA->Clone();
      }
      delete entry;
    }
  }
  
  
//update current reference maps
  if (hasAnyChanged){
    if (fCurrentRefMap) delete fCurrentRefMap;
    fCurrentRefMap=(TMap*)fRefMap->Clone();
  }
}
//_____________________________________________________________________________________
AliTPCCalPad* AliTPCcalibDButil::GetRefCalPad(AliCDBEntry *entry, const char* objName)
{
  //
  // TObjArray object type case
  // find 'objName' in 'arr' cast is to a calPad and store it in 'pad'
  //
  AliTPCCalPad *pad=0x0;
  TObjArray *arr=(TObjArray*)entry->GetObject();
  if (!arr){
    AliError(Form("Could not get object from entry '%s'\nPlease check!!!",entry->GetId().GetPath().Data()));
    return pad;
  }
  pad=(AliTPCCalPad*)arr->FindObject(objName);
  if (!pad) {
    AliError(Form("Could not get '%s' from TObjArray in entry '%s'\nPlease check!!!",objName,entry->GetId().GetPath().Data()));
    return pad;
  }
  return (AliTPCCalPad*)pad->Clone();
}
//_____________________________________________________________________________________
AliTPCCalPad* AliTPCcalibDButil::GetRefCalPad(AliCDBEntry *entry)
{
  //
  // AliTPCCalPad object type case
  // cast object to a calPad and store it in 'pad'
  //
  AliTPCCalPad *pad=(AliTPCCalPad*)entry->GetObject();
  if (!pad) {
    AliError(Form("Could not get object from entry '%s'\nPlease check!!!",entry->GetId().GetPath().Data()));
    return 0x0;
  }
  pad=(AliTPCCalPad*)pad->Clone();
  return pad;
}
//_____________________________________________________________________________________
AliTPCCalPad* AliTPCcalibDButil::GetAltroMasked(const char* cdbPath, const char* name)
{
  //
  // set altro masked channel map for 'cdbPath'
  //
  AliTPCCalPad* pad=0x0;
  const Int_t run=GetReferenceRun(cdbPath);
  if (run<0) {
    AliError(Form("Could not get reference run number for object '%s'\nPlease check availability!!!",cdbPath));
    return pad;
  }
  AliCDBEntry *entry=AliCDBManager::Instance()->Get("TPC/Calib/AltroConfig", run);
  if (!entry) {
    AliError(Form("Could not get reference object '%s'\nPlease check availability!!!",cdbPath));
    return pad;
  }
  pad=GetRefCalPad(entry,"Masked");
  if (pad) pad->SetNameTitle(name,name);
  entry->SetOwner(kTRUE);
  delete entry;
  return pad;
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::SetReferenceRun(Int_t run){
  //
  // Get Reference map
  //
  if (run<0) run=fCalibDB->GetRun();
  TString cdbPath="TPC/Calib/Ref";
  AliCDBEntry *entry=AliCDBManager::Instance()->Get(cdbPath.Data(), run);
  if (!entry) {
    AliError(Form("Could not get reference object '%s'\nPlease check availability!!!",cdbPath.Data()));
    fRefMap=0;
    return;
  }  
  entry->SetOwner(kTRUE);
  fRefMap=(TMap*)(entry->GetObject());
  AliCDBId &id=entry->GetId();
  fRefValidity.Form("%d_%d_v%d_s%d",id.GetFirstRun(),id.GetLastRun(),id.GetVersion(),id.GetSubVersion());
}
//_____________________________________________________________________________________
Bool_t AliTPCcalibDButil::HasRefChanged(const char *cdbPath)
{
  //
  // check whether a reference cdb entry has changed
  //
  if (!fCurrentRefMap) return kTRUE;
  if (GetReferenceRun(cdbPath)!=GetCurrentReferenceRun(cdbPath)) return kTRUE;
  return kFALSE;
}
//_____________________________________________________________________________________
AliCDBEntry* AliTPCcalibDButil::GetRefEntry(const char* cdbPath)
{
  //
  // get the reference AliCDBEntry for 'cdbPath'
  //
  const Int_t run=GetReferenceRun(cdbPath);
  if (run<0) {
    AliError(Form("Could not get reference run number for object '%s'\nPlease check availability!!!",cdbPath));
    return 0;
  }
  AliCDBEntry *entry=AliCDBManager::Instance()->Get(cdbPath, run);
  if (!entry) {
    AliError(Form("Could not get reference object '%s'\nPlease check availability!!!",cdbPath));
    return 0;
  }
  return entry;
}
//_____________________________________________________________________________________
Int_t AliTPCcalibDButil::GetCurrentReferenceRun(const char* type) const {
  //
  // Get reference run number for the specified OCDB path
  //
  if (!fCurrentRefMap) return -2;
  TObjString *str=dynamic_cast<TObjString*>(fCurrentRefMap->GetValue(type));
  if (!str) return -2;
  return (Int_t)str->GetString().Atoi();
}
//_____________________________________________________________________________________
Int_t AliTPCcalibDButil::GetReferenceRun(const char* type) const{
  //
  // Get reference run number for the specified OCDB path
  //
  if (!fRefMap) return -1;
  TObjString *str=dynamic_cast<TObjString*>(fRefMap->GetValue(type));
  if (!str) return -1;
  return (Int_t)str->GetString().Atoi();
}
//_____________________________________________________________________________________
AliTPCCalPad *AliTPCcalibDButil::CreateCEOutlyerMap( Int_t & noutliersCE, AliTPCCalPad * const ceOut, Float_t minSignal, Float_t cutTrmsMin,  Float_t cutTrmsMax, Float_t cutMaxDistT){
  //
  // Author:  marian.ivanov@cern.ch
  //
  // Create outlier map for CE study
  // Parameters:
  //  Return value - outlyer map
  //  noutlyersCE  - number of outlyers
  //  minSignal    - minimal total Q signal
  //  cutRMSMin    - minimal width of the signal in respect to the median 
  //  cutRMSMax    - maximal width of the signal in respect to the median 
  //  cutMaxDistT  - maximal deviation from time median per chamber
  //
  // Outlyers criteria:
  // 0. Exclude masked pads
  // 1. Exclude first two rows in IROC and last two rows in OROC
  // 2. Exclude edge pads
  // 3. Exclude channels with too large variations
  // 4. Exclude pads with too small signal
  // 5. Exclude signal with outlyers RMS
  // 6. Exclude channels to far from the chamber median	
  noutliersCE=0;
  //create outlier map
  AliTPCCalPad *out=ceOut;
  if (!out)     out= new AliTPCCalPad("outCE","outCE");
  AliTPCCalROC *rocMasked=0x0; 
  if (!fCETmean) return 0;
  if (!fCETrms) return 0;
  if (!fCEQmean) return 0;
  //
  //loop over all channels
  //
  Double_t rmsMedian         = fCETrms->GetMedian();
  for (UInt_t iroc=0;iroc<fCETmean->kNsec;++iroc){
    AliTPCCalROC *rocData=fCETmean->GetCalROC(iroc);
    if (!rocData) continue;
    if (fALTROMasked) rocMasked= fALTROMasked->GetCalROC(iroc);
    AliTPCCalROC *rocOut       = out->GetCalROC(iroc);
    AliTPCCalROC *rocCEQ       = fCEQmean->GetCalROC(iroc);
    AliTPCCalROC *rocCETrms    = fCETrms->GetCalROC(iroc);
    Double_t trocMedian        = rocData->GetMedian();
    //
    if (!rocData || !rocCEQ || !rocCETrms || !rocData) {
      noutliersCE+=AliTPCROC::Instance()->GetNChannels(iroc);
      rocOut->Add(1.);
      continue;
    }
    //
    //select outliers
    UInt_t nrows=rocData->GetNrows();
    for (UInt_t irow=0;irow<nrows;++irow){
      UInt_t npads=rocData->GetNPads(irow);
      for (UInt_t ipad=0;ipad<npads;++ipad){
        rocOut->SetValue(irow,ipad,0);
        Float_t valTmean=rocData->GetValue(irow,ipad);
        Float_t valQmean=rocCEQ->GetValue(irow,ipad);
        Float_t valTrms =rocCETrms->GetValue(irow,ipad);
        //0. exclude masked pads
        if (rocMasked && rocMasked->GetValue(irow,ipad)) {
          rocOut->SetValue(irow,ipad,1);
          continue;
        }
        //1. exclude first two rows in IROC and last two rows in OROC
        if (iroc<36){
          if (irow<2) rocOut->SetValue(irow,ipad,1);
        } else {
          if (irow>nrows-3) rocOut->SetValue(irow,ipad,1);
        }
        //2. exclude edge pads
        if (ipad==0||ipad==npads-1) rocOut->SetValue(irow,ipad,1);
        //exclude values that are exactly 0
        if ( TMath::Abs(valTmean)<kAlmost0) {
          rocOut->SetValue(irow,ipad,1);
          ++noutliersCE;
        }
        //3.  exclude channels with too large variations
        if (TMath::Abs(valTmean)>fCETmaxLimitAbs) {
          rocOut->SetValue(irow,ipad,1);
          ++noutliersCE;
        }
	//
        //4.  exclude channels with too small signal
        if (valQmean<minSignal) {
          rocOut->SetValue(irow,ipad,1);
          ++noutliersCE;
        }
        //
	//5. exclude channels with too small rms
	if (valTrms<cutTrmsMin*rmsMedian || valTrms>cutTrmsMax*rmsMedian){
	  rocOut->SetValue(irow,ipad,1);
          ++noutliersCE;
	}
        //
	//6. exclude channels to far from the chamber median	
	if (TMath::Abs(valTmean-trocMedian)>cutMaxDistT){
	  rocOut->SetValue(irow,ipad,1);
          ++noutliersCE;
	}
      }
    }
  }
  //
  return out;
}


AliTPCCalPad *AliTPCcalibDButil::CreatePulserOutlyerMap(Int_t &noutliersPulser, AliTPCCalPad * const pulserOut,Float_t cutTime, Float_t cutnRMSQ, Float_t cutnRMSrms){
  //
  // Author: marian.ivanov@cern.ch
  //
  // Create outlier map for Pulser
  // Parameters:
  //  Return value     - outlyer map
  //  noutlyersPulser  - number of outlyers
  //  cutTime          - absolute cut - distance to the median of chamber
  //  cutnRMSQ         - nsigma cut from median  q distribution per chamber
  //  cutnRMSrms       - nsigma cut from median  rms distribution 
  // Outlyers criteria:
  // 0. Exclude masked pads
  // 1. Exclude time outlyers (default 3 time bins)
  // 2. Exclude q outlyers    (default 5 sigma)
  // 3. Exclude rms outlyers  (default 5 sigma)
  noutliersPulser=0;
  AliTPCCalPad *out=pulserOut;
  if (!out)     out= new AliTPCCalPad("outPulser","outPulser");
  AliTPCCalROC *rocMasked=0x0; 
  if (!fPulserTmean) return 0;
  if (!fPulserTrms) return 0;
  if (!fPulserQmean) return 0;
  //
  //loop over all channels
  //
  for (UInt_t iroc=0;iroc<fCETmean->kNsec;++iroc){
    if (fALTROMasked)   rocMasked= fALTROMasked->GetCalROC(iroc);
    AliTPCCalROC *rocData       = fPulserTmean->GetCalROC(iroc);
    AliTPCCalROC *rocOut        = out->GetCalROC(iroc);
    AliTPCCalROC *rocPulserQ    = fPulserQmean->GetCalROC(iroc);
    AliTPCCalROC *rocPulserTrms = fPulserTrms->GetCalROC(iroc);
    //
    Double_t rocMedianT         = rocData->GetMedian();
    Double_t rocMedianQ         = rocPulserQ->GetMedian();
    Double_t rocRMSQ            = rocPulserQ->GetRMS();
    Double_t rocMedianTrms      = rocPulserTrms->GetMedian();
    Double_t rocRMSTrms         = rocPulserTrms->GetRMS();
    for (UInt_t ichannel=0;ichannel<rocData->GetNchannels();++ichannel){
      rocOut->SetValue(ichannel,0);
      Float_t valTmean=rocData->GetValue(ichannel);
      Float_t valQmean=rocPulserQ->GetValue(ichannel);
      Float_t valTrms =rocPulserTrms->GetValue(ichannel);
      Float_t valMasked =0;
      if (rocMasked) valMasked = rocMasked->GetValue(ichannel);
      Int_t isOut=0;
      if (valMasked>0.5) isOut=1;
      if (TMath::Abs(valTmean-rocMedianT)>cutTime) isOut=1;
      if (TMath::Abs(valQmean-rocMedianQ)>cutnRMSQ*rocRMSQ) isOut=1;
      if (TMath::Abs(valTrms-rocMedianTrms)>cutnRMSrms*rocRMSTrms) isOut=1;
      rocOut->SetValue(ichannel,isOut);
      if (isOut) noutliersPulser++;
    }
  }
  return out;
}


AliTPCCalPad *AliTPCcalibDButil::CreatePadTime0CE(TVectorD &fitResultsA, TVectorD&fitResultsC, Int_t &nOut, Double_t &chi2A, Double_t &chi2C, const char *dumpfile){
  //
  // Author : Marian Ivanov
  // Create pad time0 correction map using information from the CE and from pulser
  //
  //
  // Return PadTime0 to be used for time0 relative alignment
  // if dump file specified intermediat results are dumped to the fiel and can be visualized 
  // using $ALICE_ROOT/TPC/script/gui application
  //
  // fitResultsA - fitParameters A side
  // fitResultsC - fitParameters C side
  // chi2A       - chi2/ndf for A side (assuming error 1 time bin)
  // chi2C       - chi2/ndf for C side (assuming error 1 time bin)
  //
  //
  // Algorithm:
  // 1. Find outlier map for CE
  // 2. Find outlier map for Pulser
  // 3. Replace outlier by median at given sector  (median without outliers)
  // 4. Substract from the CE data pulser
  // 5. Fit the CE with formula
  //    5.1) (IROC-OROC) offset
  //    5.2) gx
  //    5.3) gy
  //    5.4) (lx-xmid)
  //    5.5) (IROC-OROC)*(lx-xmid)
  //    5.6) (ly/lx)^2
  // 6. Substract gy fit dependence from the CE data
  // 7. Add pulser back to CE data  
  // 8. Replace outliers by fit value - median of diff per given chamber -GY fit
  // 9. return CE data
  //
  // Time0 <= padCE = padCEin  -padCEfitGy  - if not outlier
  // Time0 <= padCE = padFitAll-padCEfitGy  - if outlier 

  // fit formula
  const char *formulaIn="(-1.+2.*(sector<36))*0.5++gx++gy++(lx-134.)++(-1.+2.*(sector<36))*0.5*(lx-134)++((ly/lx)^2/(0.1763)^2)";
  // output for fit formula
  const char *formulaAll="1++(-1.+2.*(sector<36))*0.5++gx++gy++(lx-134.)++(-1.+2.*(sector<36))*0.5*(lx-134)++((ly/lx)^2/(0.1763)^2)";
  // gy part of formula
  const char *formulaOut="0++0*(-1.+2.*(sector<36))*0.5++0*gx++gy++0*(lx-134.)++0*(-1.+2.*(sector<36))*0.5*(lx-134)++0*((ly/lx)^2/(0.1763)^2)";
  //
  //
  if (!fCETmean) return 0;
  Double_t pgya,pgyc,pchi2a,pchi2c;
  AliTPCCalPad * padPulserOut = CreatePulserOutlyerMap(nOut);
  AliTPCCalPad * padCEOut     = CreateCEOutlyerMap(nOut);

  AliTPCCalPad * padPulser    = CreatePadTime0(1,pgya,pgyc,pchi2a,pchi2c);
  AliTPCCalPad * padCE        = new AliTPCCalPad(*fCETmean);
  AliTPCCalPad * padCEIn      = new AliTPCCalPad(*fCETmean);
  AliTPCCalPad * padOut       = new AliTPCCalPad("padOut","padOut");   
  padPulser->SetName("padPulser");
  padPulserOut->SetName("padPulserOut");
  padCE->SetName("padCE");
  padCEIn->SetName("padCEIn");
  padCEOut->SetName("padCEOut");
  padOut->SetName("padOut");

  //
  // make combined outlyers map
  // and replace outlyers in maps with median for chamber
  //
  for (UInt_t iroc=0;iroc<fCETmean->kNsec;++iroc){  
    AliTPCCalROC * rocOut       = padOut->GetCalROC(iroc);
    AliTPCCalROC * rocPulser    = padPulser->GetCalROC(iroc);
    AliTPCCalROC * rocPulserOut = padPulserOut->GetCalROC(iroc);
    AliTPCCalROC * rocCEOut     = padCEOut->GetCalROC(iroc);
    AliTPCCalROC * rocCE        = padCE->GetCalROC(iroc);
    Double_t ceMedian           = rocCE->GetMedian(rocCEOut);
    Double_t pulserMedian       = rocPulser->GetMedian(rocCEOut);
    for (UInt_t ichannel=0;ichannel<rocOut->GetNchannels();++ichannel){
      if (rocPulserOut->GetValue(ichannel)>0) {
	rocPulser->SetValue(ichannel,pulserMedian);  
	rocOut->SetValue(ichannel,1);
      }
      if (rocCEOut->GetValue(ichannel)>0) {
	rocCE->SetValue(ichannel,ceMedian);
	rocOut->SetValue(ichannel,1);
      }
    }
  }
  //
  // remove pulser time 0
  //
  padCE->Add(padPulser,-1);
  //
  // Make fits
  //
  TMatrixD dummy;
  Float_t chi2Af,chi2Cf;  
  padCE->GlobalSidesFit(padOut,formulaIn,fitResultsA,fitResultsC,dummy,dummy,chi2Af,chi2Cf);
  chi2A=chi2Af;
  chi2C=chi2Cf;
  //
  AliTPCCalPad *padCEFitGY=AliTPCCalPad::CreateCalPadFit(formulaOut,fitResultsA,fitResultsC);
  padCEFitGY->SetName("padCEFitGy");
  //
  AliTPCCalPad *padCEFit  =AliTPCCalPad::CreateCalPadFit(formulaAll,fitResultsA,fitResultsC);
  padCEFit->SetName("padCEFit");
  //
  AliTPCCalPad* padCEDiff  = new AliTPCCalPad(*padCE);
  padCEDiff->SetName("padCEDiff");
  padCEDiff->Add(padCEFit,-1.);
  //
  // 
  padCE->Add(padCEFitGY,-1.);

  padCE->Add(padPulser,1.);  
  Double_t padmedian = padCE->GetMedian();
  padCE->Add(-padmedian);  // normalize to median
  //
  // Replace outliers by fit value - median of diff per given chamber -GY fit
  //
  for (UInt_t iroc=0;iroc<fCETmean->kNsec;++iroc){  
    AliTPCCalROC * rocOut       = padOut->GetCalROC(iroc);
    AliTPCCalROC * rocCE        = padCE->GetCalROC(iroc);
    AliTPCCalROC * rocCEFit     = padCEFit->GetCalROC(iroc);
    AliTPCCalROC * rocCEFitGY   = padCEFitGY->GetCalROC(iroc);
    AliTPCCalROC * rocCEDiff    = padCEDiff->GetCalROC(iroc);
    //
    Double_t diffMedian         = rocCEDiff->GetMedian(rocOut);
    for (UInt_t ichannel=0;ichannel<rocOut->GetNchannels();++ichannel){
      if (rocOut->GetValue(ichannel)==0) continue;
      Float_t value=rocCEFit->GetValue(ichannel)-rocCEFitGY->GetValue(ichannel)-diffMedian-padmedian;
      rocCE->SetValue(ichannel,value);
    }    
  }
  //
  //
  if (dumpfile){
    //dump to the file - result can be visualized
    AliTPCPreprocessorOnline preprocesor;
    preprocesor.AddComponent(new AliTPCCalPad(*padCE));
    preprocesor.AddComponent(new AliTPCCalPad(*padCEIn));
    preprocesor.AddComponent(new AliTPCCalPad(*padCEFit));
    preprocesor.AddComponent(new AliTPCCalPad(*padOut));
    //
    preprocesor.AddComponent(new AliTPCCalPad(*padCEFitGY));
    preprocesor.AddComponent(new AliTPCCalPad(*padCEDiff));
    //
    preprocesor.AddComponent(new AliTPCCalPad(*padCEOut));
    preprocesor.AddComponent(new AliTPCCalPad(*padPulser));
    preprocesor.AddComponent(new AliTPCCalPad(*padPulserOut));
    preprocesor.DumpToFile(dumpfile);
  } 
  delete padPulser;
  delete padPulserOut;
  delete padCEIn;
  delete padCEOut;
  delete padOut;
  delete padCEDiff;
  delete padCEFitGY;
  return padCE;
}





Int_t AliTPCcalibDButil::GetNearest(TGraph *graph, Double_t xref, Double_t &dx, Double_t &y){
  //
  // find the closest point to xref  in x  direction
  // return dx and value 
  dx = 0;
  y = 0;

  if(!graph) return 0;
  if(graph->GetN() < 1) return 0;

  Int_t index=0;
  index = TMath::BinarySearch(graph->GetN(), graph->GetX(),xref);
  if (index<0) index=0;
  if(graph->GetN()==1) {
    dx = xref-graph->GetX()[index];
  }
  else {
    if (index>=graph->GetN()-1) index=graph->GetN()-2;
    if (xref-graph->GetX()[index]>graph->GetX()[index]-xref) index++;
    dx = xref-graph->GetX()[index];
  }
  y  = graph->GetY()[index];
  return index;
}

Double_t  AliTPCcalibDButil::GetTriggerOffsetTPC(Int_t run, Int_t timeStamp, Double_t deltaT, Double_t deltaTLaser, Int_t valType){
  //
  // Get the correction of the trigger offset
  // combining information from the laser track calibration 
  // and from cosmic calibration
  //
  // run       - run number
  // timeStamp - tim stamp in seconds
  // deltaT    - integration period to calculate offset 
  // deltaTLaser -max validity of laser data
  // valType   - 0 - median, 1- mean
  // 
  // Integration vaues are just recomendation - if not possible to get points
  // automatically increase the validity by factor 2  
  // (recursive algorithm until one month of data taking)
  //
  //
  const Float_t kLaserCut=0.0005;
  const Int_t   kMaxPeriod=3600*24*30*12; // one year max
  const Int_t   kMinPoints=20;
  //
  TObjArray *array =AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);
  if (!array) {
    AliTPCcalibDB::Instance()->UpdateRunInformations(run,kFALSE); 
  }
  array =AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);
  if (!array) return 0;
  //
  TGraphErrors *laserA[3]={0,0,0};
  TGraphErrors *laserC[3]={0,0,0};
  TGraphErrors *cosmicAll=0;
  laserA[1]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_A");
  laserC[1]=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_C");
  cosmicAll =(TGraphErrors*)array->FindObject("TGRAPHERRORS_MEAN_VDRIFT_COSMICS_ALL");
  //
  //
  if (!cosmicAll) return 0;
  Int_t nmeasC=cosmicAll->GetN();
  Float_t *tdelta = new Float_t[nmeasC];
  Int_t nused=0;
  for (Int_t i=0;i<nmeasC;i++){
    if (TMath::Abs(cosmicAll->GetX()[i]-timeStamp)>deltaT) continue;
    Float_t ccosmic=cosmicAll->GetY()[i];
    Double_t yA=0,yC=0,dA=0,dC=0;
    if (laserA[1]) GetNearest(laserA[1], cosmicAll->GetX()[i],dA,yA);
    if (laserC[1]) GetNearest(laserC[1], cosmicAll->GetX()[i],dC,yC);
    //yA=laserA[1]->Eval(cosmicAll->GetX()[i]);
    //yC=laserC[1]->Eval(cosmicAll->GetX()[i]);
    //
    if (TMath::Sqrt(dA*dA+dC*dC)>deltaTLaser) continue;
    Float_t claser=0;
    if (TMath::Abs(yA-yC)<kLaserCut) {
      claser=(yA-yC)*0.5;
    }else{
      if (i%2==0)  claser=yA;
      if (i%2==1)  claser=yC;
    }
    tdelta[nused]=ccosmic-claser;
    nused++;
  }
  if (nused<kMinPoints &&deltaT<kMaxPeriod) {
    delete [] tdelta;
    return  AliTPCcalibDButil::GetTriggerOffsetTPC(run, timeStamp, deltaT*2,deltaTLaser);
  }
  if (nused<kMinPoints) {
    delete [] tdelta;
    //AliWarning("AliFatal: No time offset calibration available\n");
    return 0;
  }
  Double_t median = TMath::Median(nused,tdelta);
  Double_t mean  = TMath::Mean(nused,tdelta);
  delete [] tdelta;
  return (valType==0) ? median:mean;
}

Double_t  AliTPCcalibDButil::GetVDriftTPC(Double_t &dist, Int_t run, Int_t timeStamp, Double_t deltaT, Double_t deltaTLaser, Int_t valType){
  //
  // Get the correction of the drift velocity
  // combining information from the laser track calibration 
  // and from cosmic calibration
  //
  // dist      - return value - distance to closest point in graph
  // run       - run number
  // timeStamp - tim stamp in seconds
  // deltaT    - integration period to calculate time0 offset 
  // deltaTLaser -max validity of laser data
  // valType   - 0 - median, 1- mean
  // 
  // Integration vaues are just recomendation - if not possible to get points
  // automatically increase the validity by factor 2  
  // (recursive algorithm until one month of data taking)
  //
  //
  //
  TObjArray *array =AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);
  if (!array) {
    AliTPCcalibDB::Instance()->UpdateRunInformations(run,kFALSE); 
  }
  array =AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);
  if (!array) return 0;
  TGraphErrors *cosmicAll=0;
  cosmicAll =(TGraphErrors*)array->FindObject("TGRAPHERRORS_MEAN_VDRIFT_COSMICS_ALL");
  if (!cosmicAll) return 0;
  Double_t grY=0;
  AliTPCcalibDButil::GetNearest(cosmicAll,timeStamp,dist,grY);

  Double_t t0= AliTPCcalibDButil::GetTriggerOffsetTPC(run,timeStamp, deltaT, deltaTLaser,valType);
  Double_t vcosmic =  AliTPCcalibDButil::EvalGraphConst(cosmicAll, timeStamp);
  if (timeStamp>cosmicAll->GetX()[cosmicAll->GetN()-1])  vcosmic=cosmicAll->GetY()[cosmicAll->GetN()-1];
  if (timeStamp<cosmicAll->GetX()[0])  vcosmic=cosmicAll->GetY()[0];
  return  vcosmic-t0;

  /*
    Example usage:
    
    Int_t run=89000
    TObjArray *array =AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);
    cosmicAll =(TGraphErrors*)array->FindObject("TGRAPHERRORS_MEAN_VDRIFT_COSMICS_ALL"); 
    laserA=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_A");
    //
    Double_t *yvd= new Double_t[cosmicAll->GetN()];
    Double_t *yt0= new Double_t[cosmicAll->GetN()];
    for (Int_t i=0; i<cosmicAll->GetN();i++) yvd[i]=AliTPCcalibDButil::GetVDriftTPC(run,cosmicAll->GetX()[i]);
    for (Int_t i=0; i<cosmicAll->GetN();i++) yt0[i]=AliTPCcalibDButil::GetTriggerOffsetTPC(run,cosmicAll->GetX()[i]);

    TGraph *pcosmicVd=new TGraph(cosmicAll->GetN(), cosmicAll->GetX(), yvd);
    TGraph *pcosmicT0=new TGraph(cosmicAll->GetN(), cosmicAll->GetX(), yt0);

  */
  
}

const char* AliTPCcalibDButil::GetGUIRefTreeDefaultName()
{
  //
  // Create a default name for the gui file
  //
  
  return Form("guiRefTreeRun%s.root",GetRefValidity());
}

Bool_t AliTPCcalibDButil::CreateGUIRefTree(const char* filename)
{
  //
  // Create a gui reference tree
  // if dirname and filename are empty default values will be used
  // this is the recommended way of using this function
  // it allows to check whether a file with the given run validity alredy exists
  //
  if (!AliCDBManager::Instance()->GetDefaultStorage()){
    AliError("Default Storage not set. Cannot create reference calibration Tree!");
    return kFALSE;
  }
  
  TString file=filename;
  if (file.IsNull()) file=GetGUIRefTreeDefaultName();
  
  AliTPCPreprocessorOnline prep;
  //noise and pedestals
  if (fRefPedestals) prep.AddComponent(new AliTPCCalPad(*(fRefPedestals)));
  if (fRefPadNoise ) prep.AddComponent(new AliTPCCalPad(*(fRefPadNoise)));
  if (fRefPedestalMasked) prep.AddComponent(new AliTPCCalPad(*fRefPedestalMasked));
  //pulser data
  if (fRefPulserTmean) prep.AddComponent(new AliTPCCalPad(*(fRefPulserTmean)));
  if (fRefPulserTrms ) prep.AddComponent(new AliTPCCalPad(*(fRefPulserTrms)));
  if (fRefPulserQmean) prep.AddComponent(new AliTPCCalPad(*(fRefPulserQmean)));
  if (fRefPulserMasked) prep.AddComponent(new AliTPCCalPad(*fRefPulserMasked));
  //CE data
  if (fRefCETmean) prep.AddComponent(new AliTPCCalPad(*(fRefCETmean)));
  if (fRefCETrms ) prep.AddComponent(new AliTPCCalPad(*(fRefCETrms)));
  if (fRefCEQmean) prep.AddComponent(new AliTPCCalPad(*(fRefCEQmean)));
  if (fRefCEMasked) prep.AddComponent(new AliTPCCalPad(*fRefCEMasked));
  //Altro data
  if (fRefALTROAcqStart ) prep.AddComponent(new AliTPCCalPad(*(fRefALTROAcqStart )));
  if (fRefALTROZsThr    ) prep.AddComponent(new AliTPCCalPad(*(fRefALTROZsThr    )));
  if (fRefALTROFPED     ) prep.AddComponent(new AliTPCCalPad(*(fRefALTROFPED     )));
  if (fRefALTROAcqStop  ) prep.AddComponent(new AliTPCCalPad(*(fRefALTROAcqStop  )));
  if (fRefALTROMasked   ) prep.AddComponent(new AliTPCCalPad(*(fRefALTROMasked   )));
  //QA
  AliTPCdataQA *dataQA=fRefDataQA;
  if (dataQA) {
    if (dataQA->GetNLocalMaxima())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetNLocalMaxima())));
    if (dataQA->GetMaxCharge())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetMaxCharge())));
    if (dataQA->GetMeanCharge())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetMeanCharge())));
    if (dataQA->GetNoThreshold())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetNoThreshold())));
    if (dataQA->GetNTimeBins())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetNTimeBins())));
    if (dataQA->GetNPads())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetNPads())));
    if (dataQA->GetTimePosition())
      prep.AddComponent(new AliTPCCalPad(*(dataQA->GetTimePosition())));
  }
  prep.DumpToFile(file.Data());
  return kTRUE;
}

Double_t  AliTPCcalibDButil::GetVDriftTPCLaserTracks(Double_t &dist, Int_t run, Int_t timeStamp, Double_t deltaT, Int_t side){
  //
  // Get the correction of the drift velocity using the offline laser tracks calbration
  //
  // run       - run number
  // timeStamp - tim stamp in seconds
  // deltaT    - integration period to calculate time0 offset 
  // side      - 0 - A side,  1 - C side, 2 - mean from both sides
  // Note in case no data form both A and C side - the value from active side used
  TObjArray *array =AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);

  return GetVDriftTPCLaserTracksCommon(dist, timeStamp, deltaT, side, array);
}

Double_t  AliTPCcalibDButil::GetVDriftTPCLaserTracksOnline(Double_t &dist, Int_t /*run*/, Int_t timeStamp, Double_t deltaT, Int_t side){
  //
  // Get the correction of the drift velocity using the online laser tracks calbration
  //
  // run       - run number
  // timeStamp - tim stamp in seconds
  // deltaT    - integration period to calculate time0 offset
  // side      - 0 - A side,  1 - C side, 2 - mean from both sides
  // Note in case no data form both A and C side - the value from active side used
  TObjArray *array =AliTPCcalibDB::Instance()->GetCEfitsDrift();

  Double_t dv = GetVDriftTPCLaserTracksCommon(dist, timeStamp, deltaT, side, array);
  AliTPCParam *param  =AliTPCcalibDB::Instance()->GetParameters();
  if (!param) return 0;

  //the drift velocity is hard wired in the AliTPCCalibCE class, since online there is no access to OCDB
  dv*=param->GetDriftV()/2.61301900000000000e+06;
  if (dv>1e-20) dv=1/dv-1;
  else return 0;
  // T/P correction
  TObjArray*  cearray =AliTPCcalibDB::Instance()->GetCEData();
  
  AliTPCSensorTempArray *temp = (AliTPCSensorTempArray*)cearray->FindObject("TempMap");
  AliDCSSensor *press         = (AliDCSSensor*)cearray->FindObject("CavernAtmosPressure");
  
  Double_t corrPTA=0;
  Double_t corrPTC=0;
  
  if (temp&&press) {
    AliTPCCalibVdrift corr(temp,press,0);
    corrPTA=corr.GetPTRelative(timeStamp,0);
    corrPTC=corr.GetPTRelative(timeStamp,1);
  }
  
  if (side==0) dv -=  corrPTA;
  if (side==1) dv -=  corrPTC;
  if (side==2) dv -=  (corrPTA+corrPTC)/2;
  
  return dv;
}

Double_t  AliTPCcalibDButil::GetVDriftTPCLaserTracksCommon(Double_t &dist, Int_t timeStamp, Double_t deltaT,
  Int_t side, TObjArray * const array){
  //
  // common drift velocity retrieval for online and offline method
  //
  TGraphErrors *grlaserA=0;
  TGraphErrors *grlaserC=0;
  Double_t vlaserA=0, vlaserC=0;
  if (!array) return 0;
  grlaserA=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_A");
  grlaserC=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DRIFT_LASER_ALL_C");
  Double_t deltaY;
  if (grlaserA && grlaserA->GetN()>0) {
    AliTPCcalibDButil::GetNearest(grlaserA,timeStamp,dist,deltaY);
    if (TMath::Abs(dist)>deltaT)  vlaserA= deltaY;
    else  vlaserA = AliTPCcalibDButil::EvalGraphConst(grlaserA,timeStamp);
  }
  if (grlaserC && grlaserC->GetN()>0) {
    AliTPCcalibDButil::GetNearest(grlaserC,timeStamp,dist,deltaY);
    if (TMath::Abs(dist)>deltaT)  vlaserC= deltaY;
    else  vlaserC = AliTPCcalibDButil::EvalGraphConst(grlaserC,timeStamp);
  }
  if (side==0) return vlaserA;
  if (side==1) return vlaserC;
  Double_t mdrift=(vlaserA+vlaserC)*0.5;
  if (!grlaserA) return vlaserC;
  if (!grlaserC) return vlaserA;
  return mdrift;
}


Double_t  AliTPCcalibDButil::GetVDriftTPCCE(Double_t &dist,Int_t run, Int_t timeStamp, Double_t deltaT, Int_t side){
  //
  // Get the correction of the drift velocity using the CE laser data
  // combining information from the CE,  laser track calibration
  // and P/T calibration 
  //
  // run       - run number
  // timeStamp - tim stamp in seconds
  // deltaT    - integration period to calculate time0 offset 
  // side      - 0 - A side,  1 - C side, 2 - mean from both sides
  TObjArray *arrT     =AliTPCcalibDB::Instance()->GetCErocTtime();
  if (!arrT) return 0;
  AliTPCParam *param  =AliTPCcalibDB::Instance()->GetParameters();
  TObjArray*  cearray =AliTPCcalibDB::Instance()->GetCEData(); 
  AliTPCCalibVdrift * driftCalib = (AliTPCCalibVdrift *)cearray->FindObject("driftPTCE");
  //
  //
  Double_t corrPTA = 0, corrPTC=0;
  Double_t ltime0A = 0, ltime0C=0;
  Double_t gry=0;
  Double_t corrA=0, corrC=0;
  Double_t timeA=0, timeC=0;
  const Double_t kEpsilon = 0.00001;
  TGraph *graphA = (TGraph*)arrT->At(72);
  TGraph *graphC = (TGraph*)arrT->At(73);
  if (!graphA && !graphC) return 0.;
  if (graphA &&graphA->GetN()>0) {
    AliTPCcalibDButil::GetNearest(graphA,timeStamp,dist,gry);
    timeA   = AliTPCcalibDButil::EvalGraphConst(graphA,timeStamp);
    Int_t mtime   =TMath::Nint((graphA->GetX()[0]+graphA->GetX()[graphA->GetN()-1])*0.5);
    ltime0A       = GetLaserTime0(run,mtime,TMath::Nint(deltaT),0);
    if(ltime0A < kEpsilon) return 0;
    if (driftCalib) corrPTA =  driftCalib->GetPTRelative(timeStamp,0);
    corrA = (param->GetZLength(36)/(timeA*param->GetTSample()*(1.-ltime0A)-param->GetL1Delay()-0*param->GetZSigma()/param->GetDriftV()))/param->GetDriftV()-1;
    corrA-=corrPTA;
  }
  if (graphC&&graphC->GetN()>0){
    AliTPCcalibDButil::GetNearest(graphC,timeStamp,dist,gry);
    timeC=AliTPCcalibDButil::EvalGraphConst(graphC,timeStamp);
    Int_t mtime=TMath::Nint((graphC->GetX()[0]+graphC->GetX()[graphC->GetN()-1])*0.5);
    ltime0C       = GetLaserTime0(run,mtime,TMath::Nint(deltaT),0);
    if(ltime0C < kEpsilon) return 0;   
if (driftCalib) corrPTC =  driftCalib->GetPTRelative(timeStamp,0);
    corrC = (param->GetZLength(54)/(timeC*param->GetTSample()*(1.-ltime0C)-param->GetL1Delay()-0*param->GetZSigma()/param->GetDriftV()))/param->GetDriftV()-1;
    corrC-=corrPTC;
  }
  
  if (side ==0 ) return corrA;
  if (side ==1 ) return corrC;
  Double_t corrM= (corrA+corrC)*0.5;
  if (!graphA) corrM=corrC; 
  if (!graphC) corrM=corrA; 
  return corrM;
}

Double_t  AliTPCcalibDButil::GetVDriftTPCITS(Double_t &dist, Int_t run, Int_t timeStamp){
  //
  // return drift velocity using the TPC-ITS matchin method
  // return also distance to the closest point
  //
  TObjArray *array =AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);
  TGraphErrors *graph=0;
  dist=0;
  if (!array) return 0;
  //array->ls();
  graph = (TGraphErrors*)array->FindObject("ALIGN_ITSB_TPC_DRIFTVD");
  if (!graph) return 0;
  Double_t deltaY;
  AliTPCcalibDButil::GetNearest(graph,timeStamp,dist,deltaY); 
  Double_t value = AliTPCcalibDButil::EvalGraphConst(graph,timeStamp);
  return value;
}

Double_t AliTPCcalibDButil::GetTime0TPCITS(Double_t &dist, Int_t run, Int_t timeStamp){
  //
  // Get time dependent time 0 (trigger delay in cm) correction
  // Arguments:
  // timestamp - timestamp
  // run       - run number
  //
  // Notice - Extrapolation outside of calibration range  - using constant function
  //
  TObjArray *array =AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);
  TGraphErrors *graph=0;
  dist=0;
  if (!array) return 0;
  graph = (TGraphErrors*)array->FindObject("ALIGN_ITSM_TPC_T0");
  if (!graph) return 0;
  Double_t deltaY;
  AliTPCcalibDButil::GetNearest(graph,timeStamp,dist,deltaY); 
  Double_t value = AliTPCcalibDButil::EvalGraphConst(graph,timeStamp);
  return value;
}





Int_t  AliTPCcalibDButil::MakeRunList(Int_t startRun, Int_t stopRun){
  //
  // VERY obscure method - we need something in framework
  // Find the TPC runs with temperature OCDB entry
  // cache the start and end of the run
  //
  AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("TPC/Calib/Temperature");
  if (!storage) storage = AliCDBManager::Instance()->GetDefaultStorage();
  if (!storage) return 0;
  TString path=storage->GetURI(); 
  TString runsT;
  {    
    TString command;
    if (path.Contains("local")){  // find the list if local system
      path.ReplaceAll("local://","");
      path+="TPC/Calib/Temperature";
      command=Form("ls %s  |  sed s/_/\\ /g | awk '{print \"r\"$2}'  ",path.Data());
    }
    runsT=gSystem->GetFromPipe(command);
  }
  TObjArray *arr= runsT.Tokenize("r");
  if (!arr) return 0;
  //
  TArrayI indexes(arr->GetEntries());
  TArrayI runs(arr->GetEntries());
  Int_t naccept=0;
  {for (Int_t irun=0;irun<arr->GetEntries();irun++){
      Int_t irunN = atoi(arr->At(irun)->GetName());
      if (irunN<startRun) continue;
      if (irunN>stopRun) continue;
      runs[naccept]=irunN;
      naccept++;
    }}
  fRuns.Set(naccept);
  fRunsStart.Set(fRuns.fN);
  fRunsStop.Set(fRuns.fN);
  TMath::Sort(fRuns.fN, runs.fArray, indexes.fArray,kFALSE);
  for (Int_t irun=0; irun<fRuns.fN; irun++)  fRuns[irun]=runs[indexes[irun]];
  
  //
  AliCDBEntry * entry = 0;
  {for (Int_t irun=0;irun<fRuns.fN; irun++){
      entry = AliCDBManager::Instance()->Get("TPC/Calib/Temperature",fRuns[irun]);
      if (!entry) continue;
      AliTPCSensorTempArray *  tmpRun = dynamic_cast<AliTPCSensorTempArray*>(entry->GetObject());
      if (!tmpRun) continue;
      fRunsStart[irun]=tmpRun->GetStartTime().GetSec();
      fRunsStop[irun]=tmpRun->GetEndTime().GetSec();
      //AliInfo(Form("irun\t%d\tRun\t%d\t%d\t%d\n",irun,fRuns[irun],tmpRun->GetStartTime().GetSec(),tmpRun->GetEndTime().GetSec()));
    }}
  return fRuns.fN;
}


Int_t AliTPCcalibDButil::FindRunTPC(Int_t    itime, Bool_t debug){
  //
  // binary search - find the run for given time stamp
  //
  Int_t index0  = TMath::BinarySearch(fRuns.fN, fRunsStop.fArray,itime);
  Int_t index1  = TMath::BinarySearch(fRuns.fN, fRunsStart.fArray,itime);
  Int_t cindex  = -1;
  for (Int_t index=index0; index<=index1; index++){
    if (fRunsStart[index]<=itime && fRunsStop[index]>=itime) cindex=index;
    if (debug) {
      AliInfo(Form("%d\t%d\t%d\n",fRuns[index], fRunsStart[index]-itime, fRunsStop[index]-itime));
    }
  }
  if (cindex<0) cindex =(index0+index1)/2;
  if (cindex<0) {
    return 0; 
  }
  return fRuns[cindex];
}





TGraph* AliTPCcalibDButil::FilterGraphMedian(TGraph * graph, Float_t sigmaCut,Double_t &medianY){
  //
  // filter outlyer measurement
  // Only points around median +- sigmaCut filtered 
  if (!graph) return  0;
  Int_t kMinPoints=2;
  Int_t npoints0 = graph->GetN();
  Int_t npoints=0;
  Float_t  rmsY=0;
  //
  //
  if (npoints0<kMinPoints) return 0;

  Double_t *outx=new Double_t[npoints0];
  Double_t *outy=new Double_t[npoints0];
  for (Int_t iter=0; iter<3; iter++){
    npoints=0;
    for (Int_t ipoint=0; ipoint<npoints0; ipoint++){
      if (graph->GetY()[ipoint]==0) continue;
      if (iter>0 &&TMath::Abs(graph->GetY()[ipoint]-medianY)>sigmaCut*rmsY) continue;  
      outx[npoints]  = graph->GetX()[ipoint];
      outy[npoints]  = graph->GetY()[ipoint];
      npoints++;
    }
    if (npoints<=1) break;
    medianY  =TMath::Median(npoints,outy);
    rmsY   =TMath::RMS(npoints,outy);
  }
  TGraph *graphOut=0;
  if (npoints>1) graphOut= new TGraph(npoints,outx,outy); 
  delete [] outx;
  delete [] outy;
  return graphOut;
}


TGraph* AliTPCcalibDButil::FilterGraphMedianAbs(TGraph * graph, Float_t cut,Double_t &medianY){
  //
  // filter outlyer measurement
  // Only points around median +- cut filtered 
  if (!graph) return  0;
  Int_t kMinPoints=2;
  Int_t npoints0 = graph->GetN();
  Int_t npoints=0;
  Float_t  rmsY=0;
  //
  //
  if (npoints0<kMinPoints) return 0;

  Double_t *outx=new Double_t[npoints0];
  Double_t *outy=new Double_t[npoints0];
  for (Int_t iter=0; iter<3; iter++){
    npoints=0;
    for (Int_t ipoint=0; ipoint<npoints0; ipoint++){
      if (graph->GetY()[ipoint]==0) continue;
      if (iter>0 &&TMath::Abs(graph->GetY()[ipoint]-medianY)>cut) continue;  
      outx[npoints]  = graph->GetX()[ipoint];
      outy[npoints]  = graph->GetY()[ipoint];
      npoints++;
    }
    if (npoints<=1) break;
    medianY  =TMath::Median(npoints,outy);
    rmsY   =TMath::RMS(npoints,outy);
  }
  TGraph *graphOut=0;
  if (npoints>1) graphOut= new TGraph(npoints,outx,outy); 
  delete [] outx;
  delete [] outy;
  return graphOut;
}



TGraphErrors* AliTPCcalibDButil::FilterGraphMedianErr(TGraphErrors * const graph, Float_t sigmaCut,Double_t &medianY){
  //
  // filter outlyer measurement
  // Only points with normalized errors median +- sigmaCut filtered
  //
  Int_t kMinPoints=10;
  Int_t npoints0 = graph->GetN();
  Int_t npoints=0;
  Float_t  medianErr=0, rmsErr=0;
  //
  //
  if (npoints0<kMinPoints) return 0;

  Double_t *outx=new Double_t[npoints0];
  Double_t *outy=new Double_t[npoints0];
  Double_t *erry=new Double_t[npoints0];
  Double_t *nerry=new Double_t[npoints0];
  Double_t *errx=new Double_t[npoints0];

  for (Int_t iter=0; iter<3; iter++){
    npoints=0;
    for (Int_t ipoint=0; ipoint<npoints0; ipoint++){
      nerry[npoints]  = graph->GetErrorY(ipoint);
      if (iter>0 &&TMath::Abs(nerry[npoints]-medianErr)>sigmaCut*rmsErr) continue;  
      erry[npoints]  = graph->GetErrorY(ipoint);
      outx[npoints]  = graph->GetX()[ipoint];
      outy[npoints]  = graph->GetY()[ipoint];
      errx[npoints]  = graph->GetErrorY(ipoint);
      npoints++;
    }
    if (npoints==0) break;
    medianErr=TMath::Median(npoints,erry);
    medianY  =TMath::Median(npoints,outy);
    rmsErr   =TMath::RMS(npoints,erry);
  }
  TGraphErrors *graphOut=0;
  if (npoints>1) graphOut= new TGraphErrors(npoints,outx,outy,errx,erry); 
  delete []outx;
  delete []outy;
  delete []erry;
  delete []nerry;
  delete []errx;
  return graphOut;
}


void AliTPCcalibDButil::Sort(TGraph *graph){
  //
  // sort array - neccessay for approx
  //
  Int_t npoints = graph->GetN();
  Int_t *indexes=new Int_t[npoints];
  Double_t *outx=new Double_t[npoints];
  Double_t *outy=new Double_t[npoints];
  TMath::Sort(npoints, graph->GetX(),indexes,kFALSE);
  for (Int_t i=0;i<npoints;i++) outx[i]=graph->GetX()[indexes[i]];
  for (Int_t i=0;i<npoints;i++) outy[i]=graph->GetY()[indexes[i]];
  for (Int_t i=0;i<npoints;i++) graph->GetX()[i]=outx[i];
  for (Int_t i=0;i<npoints;i++) graph->GetY()[i]=outy[i];
 
  delete [] indexes;
  delete [] outx;
  delete [] outy;
}
void AliTPCcalibDButil::SmoothGraph(TGraph *graph, Double_t delta){
  //
  // smmoth graph - mean on the interval
  //
  Sort(graph);
  Int_t npoints = graph->GetN();
  Double_t *outy=new Double_t[npoints];
  
  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    Double_t lx=graph->GetX()[ipoint];
    Int_t index0=TMath::BinarySearch(npoints, graph->GetX(),lx-delta);
    Int_t index1=TMath::BinarySearch(npoints, graph->GetX(),lx+delta);
    if (index0<0) index0=0;
    if (index1>=npoints-1) index1=npoints-1;
    if ((index1-index0)>1){
      outy[ipoint]  = TMath::Mean(index1-index0, &(graph->GetY()[index0]));
    }else{
      outy[ipoint]=graph->GetY()[ipoint];
    }
  }
 //  TLinearFitter  fitter(3,"pol2");
//   for (Int_t ipoint=0; ipoint<npoints; ipoint++){
//     Double_t lx=graph->GetX()[ipoint];
//     Int_t index0=TMath::BinarySearch(npoints, graph->GetX(),lx-delta);
//     Int_t index1=TMath::BinarySearch(npoints, graph->GetX(),lx+delta);
//     if (index0<0) index0=0;
//     if (index1>=npoints-1) index1=npoints-1;
//     fitter.ClearPoints();
//     for (Int_t jpoint=0;jpoint<index1-index0; jpoint++)
//     if ((index1-index0)>1){
//       outy[ipoint]  = TMath::Mean(index1-index0, &(graph->GetY()[index0]));
//     }else{
//       outy[ipoint]=graph->GetY()[ipoint];
//     }
//   }



  for (Int_t ipoint=0; ipoint<npoints; ipoint++){
    graph->GetY()[ipoint] = outy[ipoint];
  }
  delete[] outy;
}

Double_t AliTPCcalibDButil::EvalGraphConst(TGraph * const graph, Double_t xref){
  //
  // Use constant interpolation outside of range 
  //
  if (!graph) {
    AliInfoGeneral("AliTPCcalibDButil","AliTPCcalibDButil::EvalGraphConst: 0 pointer\n");
    return 0;
  }

  if (graph->GetN()<1){
    AliInfoGeneral("AliTPCcalibDButil","AliTPCcalibDButil::EvalGraphConst: Empty graph \n");
    return 0;
  }
 

  if (xref<graph->GetX()[0]) return graph->GetY()[0];
  if (xref>graph->GetX()[graph->GetN()-1]) return graph->GetY()[graph->GetN()-1]; 

  //  AliInfo(Form("graph->Eval(graph->GetX()[0]) %f, graph->Eval(xref) %f \n",graph->Eval(graph->GetX()[0]), graph->Eval(xref)));

  if(graph->GetN()==1)
    return graph->Eval(graph->GetX()[0]);


  return graph->Eval(xref);
}

Double_t AliTPCcalibDButil::EvalGraphConst(AliSplineFit *graph, Double_t xref){
  //
  // Use constant interpolation outside of range also for spline fits
  //
  if (!graph) {
    AliInfoGeneral("AliTPCcalibDButil","AliTPCcalibDButil::EvalGraphConst: 0 pointer\n");
    return 0;
  }
  if (graph->GetKnots()<1){
    AliInfoGeneral("AliTPCcalibDButil","AliTPCcalibDButil::EvalGraphConst: Empty graph");
    return 0;
  }
  if (xref<graph->GetX()[0]) return graph->GetY0()[0];
  if (xref>graph->GetX()[graph->GetKnots()-1]) return graph->GetY0()[graph->GetKnots()-1]; 
  return graph->Eval( xref);
}

Float_t AliTPCcalibDButil::FilterSensor(AliDCSSensor * sensor, Double_t ymin, Double_t ymax, Double_t maxdy,  Double_t sigmaCut){
  //
  // Filter DCS sensor information
  //   ymin     - minimal value
  //   ymax     - max value
  //   maxdy    - maximal deirivative
  //   sigmaCut - cut on values and derivative in terms of RMS distribution
  // Return value - accepted fraction
  // 
  // Algorithm:
  //
  // 0. Calculate median and rms of values in specified range
  // 1. Filter out outliers - median+-sigmaCut*rms
  //    values replaced by median
  //
  AliSplineFit * fit    = sensor->GetFit();
  if (!fit) return 0.;
  Int_t          nknots = fit->GetKnots();
  if (nknots==0) {
    delete fit;
    sensor->SetFit(0);
    return 0;
  }
  //
  Double_t *yin0  = new Double_t[nknots];
  Double_t *yin1  = new Double_t[nknots];
  Int_t naccept=0;
  
  for (Int_t iknot=0; iknot< nknots; iknot++){
    if (fit->GetY0()[iknot]>ymin && fit->GetY0()[iknot]<ymax){
      yin0[naccept]  = fit->GetY0()[iknot];
      yin1[naccept]  = fit->GetY1()[iknot];
      if (TMath::Abs(fit->GetY1()[iknot])>maxdy) yin1[naccept]=0;
      naccept++;
    }
  }
  if (naccept<1) {
    delete fit;
    sensor->SetFit(0);
    delete [] yin0;
    delete [] yin1;
    return 0.;
  }

  Double_t medianY0=0, medianY1=0;
  Double_t rmsY0   =0, rmsY1=0;
  medianY0 = TMath::Median(naccept, yin0);
  medianY1 = TMath::Median(naccept, yin1);
  rmsY0    = TMath::RMS(naccept, yin0);
  rmsY1    = TMath::RMS(naccept, yin1);
  naccept=0;
  //
  // 1. Filter out outliers - median+-sigmaCut*rms
  //    values replaced by median
  //    if replaced the derivative set to 0
  //
  for (Int_t iknot=0; iknot< nknots; iknot++){
    Bool_t isOK=kTRUE;
    if (TMath::Abs(fit->GetY0()[iknot]-medianY0)>sigmaCut*rmsY0) isOK=kFALSE;
    if (TMath::Abs(fit->GetY1()[iknot]-medianY1)>sigmaCut*rmsY1) isOK=kFALSE;
    if (nknots<2) fit->GetY1()[iknot]=0;
    if (TMath::Abs(fit->GetY1()[iknot])>maxdy) fit->GetY1()[iknot]=0;
    if (!isOK){
      fit->GetY0()[iknot]=medianY0;
      fit->GetY1()[iknot]=0;
    }else{
      naccept++;
    }
  }
  delete [] yin0;
  delete [] yin1;
  return Float_t(naccept)/Float_t(nknots);
}

Float_t  AliTPCcalibDButil::FilterTemperature(AliTPCSensorTempArray *tempArray, Double_t ymin, Double_t ymax, Double_t sigmaCut){
  //
  // Filter temperature array
  // tempArray    - array of temperatures         -
  // ymin         - minimal accepted temperature  - default 15
  // ymax         - maximal accepted temperature  - default 22
  // sigmaCut     - values filtered on interval median+-sigmaCut*rms - defaut 5
  // return value - fraction of filtered sensors
  const Double_t kMaxDy=0.1;
  Int_t nsensors=tempArray->NumSensors();
  if (nsensors==0) return 0.;
  Int_t naccept=0;
  for (Int_t isensor=0; isensor<nsensors; isensor++){
    AliDCSSensor *sensor = tempArray->GetSensorNum(isensor);
    if (!sensor) continue;
    FilterSensor(sensor,ymin,ymax,kMaxDy, sigmaCut);
    if (sensor->GetFit()==0){
      //delete sensor;
      tempArray->RemoveSensorNum(isensor);
    }else{
      naccept++;
    }
  }
  return Float_t(naccept)/Float_t(nsensors);
}


void AliTPCcalibDButil::FilterCE(Double_t deltaT, Double_t cutAbs, Double_t cutSigma, TTreeSRedirector * const pcstream){
  //
  // Filter CE data
  // Input parameters:
  //    deltaT   - smoothing window (in seconds)
  //    cutAbs   - max distance of the time info to the median (in time bins)
  //    cutSigma - max distance (in the RMS)
  //    pcstream - optional debug streamer to store original and filtered info
  // Hardwired parameters:
  //    kMinPoints =10;       // minimal number of points to define the CE
  //    kMinSectors=12;       // minimal number of sectors to define sideCE
  // Algorithm:
  // 0. Filter almost emty graphs (kMinPoints=10)
  // 1. calculate median and RMS per side
  // 2. Filter graphs - in respect with side medians 
  //                  - cutAbs and cutDelta used
  // 3. Cut in respect wit the graph median - cutAbs and cutRMS used
  // 4. Calculate mean for A side and C side
  //
  const Int_t kMinPoints =10;       // minimal number of points to define the CE
  const Int_t kMinSectors=12;       // minimal number of sectors to define sideCE
  const Int_t kMinTime   =400;     // minimal arrival time of CE
  TObjArray *arrT=AliTPCcalibDB::Instance()->GetCErocTtime();
  Double_t medianY=0;
  TObjArray*  cearray =AliTPCcalibDB::Instance()->GetCEData(); 
  if (!cearray) return;
  Double_t tmin=-1;
  Double_t tmax=-1;
  //
  //
  AliTPCSensorTempArray *tempMapCE = (AliTPCSensorTempArray *)cearray->FindObject("TempMap");
  AliDCSSensor * cavernPressureCE  = (AliDCSSensor *) cearray->FindObject("CavernAtmosPressure");
  if ( tempMapCE && cavernPressureCE){
    //
    //     Bool_t isOK = FilterTemperature(tempMapCE)>0.1;
    //     FilterSensor(cavernPressureCE,960,1050,10, 5.);
    //     if (cavernPressureCE->GetFit()==0) isOK=kFALSE;
    Bool_t isOK=kTRUE;
    if (isOK)  {      
      // recalculate P/T correction map for time of the CE
      AliTPCCalibVdrift * driftCalib = new AliTPCCalibVdrift(tempMapCE,cavernPressureCE ,0);
      driftCalib->SetName("driftPTCE");
      driftCalib->SetTitle("driftPTCE");
      cearray->AddLast(driftCalib);
    }
  }
  //
  // 0. Filter almost emty graphs
  //

  for (Int_t i=0; i<72;i++){
    TGraph *graph= (TGraph*)arrT->At(i);
    if (!graph) continue; 
    graph->Sort();
    if (graph->GetN()<kMinPoints){
      arrT->AddAt(0,i);
      delete graph;  // delete empty graph
      continue;
    }
    if (tmin<0) tmin = graph->GetX()[0];
    if (tmax<0) tmax = graph->GetX()[graph->GetN()-1];
    //
    if (tmin>graph->GetX()[0]) tmin=graph->GetX()[0];
    if (tmax<graph->GetX()[graph->GetN()-1]) tmax=graph->GetX()[graph->GetN()-1];
  }
  //
  // 1. calculate median and RMS per side
  //
  TArrayF arrA(100000), arrC(100000);
  Int_t nA=0, nC=0;
  Double_t medianA=0, medianC=0;
  Double_t rmsA=0, rmsC=0;
  for (Int_t isec=0; isec<72;isec++){
    TGraph *graph= (TGraph*)arrT->At(isec);
    if (!graph) continue;
    for (Int_t ipoint=kMinPoints-1; ipoint<graph->GetN();ipoint++){
      if (graph->GetY()[ipoint]<kMinTime) continue;
      if (nA>=arrA.fN) arrA.Set(nA*2);
      if (nC>=arrC.fN) arrC.Set(nC*2);
      if (isec%36<18)  arrA[nA++]= graph->GetY()[ipoint];
      if (isec%36>=18) arrC[nC++]= graph->GetY()[ipoint];
    }
  }
  if (nA>0){
    medianA=TMath::Median(nA,arrA.fArray);
    rmsA   =TMath::RMS(nA,arrA.fArray);
  }
  if (nC>0){
    medianC=TMath::Median(nC,arrC.fArray);
    rmsC   =TMath::RMS(nC,arrC.fArray);
  }
  //
  // 2. Filter graphs - in respect with side medians
  //  
  TArrayD vecX(100000), vecY(100000);
  for (Int_t isec=0; isec<72;isec++){
    TGraph *graph= (TGraph*)arrT->At(isec);
    if (!graph) continue;
    Double_t median = (isec%36<18) ? medianA: medianC;
    Double_t rms    = (isec%36<18) ? rmsA:    rmsC;
    Int_t naccept=0;
    //    for (Int_t ipoint=kMinPoints-1; ipoint<graph->GetN();ipoint++){ //not neccessary to remove first points
    for (Int_t ipoint=0; ipoint<graph->GetN();ipoint++){
      if (TMath::Abs(graph->GetY()[ipoint]-median)>cutAbs) continue;
      if (TMath::Abs(graph->GetY()[ipoint]-median)>cutSigma*rms) continue;
      vecX[naccept]= graph->GetX()[ipoint];
      vecY[naccept]= graph->GetY()[ipoint];
      naccept++;
    }
    if (naccept<kMinPoints){
      arrT->AddAt(0,isec);
      delete graph;  // delete empty graph
      continue;
    }
    TGraph *graph2 = new TGraph(naccept, vecX.fArray, vecY.fArray);
    delete graph;
    arrT->AddAt(graph2,isec);
  }
  //
  // 3. Cut in respect wit the graph median
  //
  for (Int_t i=0; i<72;i++){
    TGraph *graph= (TGraph*)arrT->At(i);
    if (!graph) continue;
    //
    // filter in range
    //
    TGraph* graphTS0= FilterGraphMedianAbs(graph,cutAbs,medianY);
    if (!graphTS0) continue;
    if (graphTS0->GetN()<kMinPoints) {
      delete graphTS0;  
      delete graph;
      arrT->AddAt(0,i);
      continue;
    }
    TGraph* graphTS= FilterGraphMedian(graphTS0,cutSigma,medianY);    
    if (!graphTS) continue;
    graphTS->Sort();
    AliTPCcalibDButil::SmoothGraph(graphTS,deltaT);      
    if (pcstream){
      Int_t run = AliTPCcalibDB::Instance()->GetRun();
      (*pcstream)<<"filterCE"<<
	"run="<<run<<
	"isec="<<i<<
	"mY="<<medianY<<
	"graph.="<<graph<<
	"graphTS0.="<<graphTS0<<
	"graphTS.="<<graphTS<<
	"\n";
    }
    delete graphTS0;
    arrT->AddAt(graphTS,i);
    delete graph;
  }
  //
  // Recalculate the mean time A side C side
  //
  TArrayF xA(200), yA(200), eA(200), xC(200),yC(200), eC(200);
  Int_t meanPoints=(nA+nC)/72;  // mean number of points
  for (Int_t itime=0; itime<200; itime++){
    nA=0, nC=0;
    Double_t time=tmin+(tmax-tmin)*Float_t(itime)/200.;
    for (Int_t i=0; i<72;i++){
      TGraph *graph= (TGraph*)arrT->At(i);
      if (!graph) continue;
      if (graph->GetN()<(meanPoints/4)) continue;
      if ( (i%36)<18 )  arrA[nA++]=graph->Eval(time);
      if ( (i%36)>=18 ) arrC[nC++]=graph->Eval(time);
    }
    xA[itime]=time;
    xC[itime]=time;
    yA[itime]=(nA>0)? TMath::Mean(nA,arrA.fArray):0;
    yC[itime]=(nC>0)? TMath::Mean(nC,arrC.fArray):0;
    eA[itime]=(nA>0)? TMath::RMS(nA,arrA.fArray):0;
    eC[itime]=(nC>0)? TMath::RMS(nC,arrC.fArray):0;
  }
  //
  Double_t rmsTA = TMath::RMS(200,yA.fArray)+TMath::Mean(200,eA.fArray);
  Double_t rmsTC = TMath::RMS(200,yC.fArray)+TMath::Mean(200,eC.fArray);
  if (pcstream){
    Int_t run = AliTPCcalibDB::Instance()->GetRun();
    (*pcstream)<<"filterAC"<<
      "run="<<run<<
      "nA="<<nA<<
      "nC="<<nC<<
      "rmsTA="<<rmsTA<<
      "rmsTC="<<rmsTC<<
      "\n";
  }
  //
  TGraphErrors *grA = new TGraphErrors(200,xA.fArray,yA.fArray,0, eA.fArray);
  TGraphErrors *grC = new TGraphErrors(200,xC.fArray,yC.fArray,0, eC.fArray);
  TGraph* graphTSA= FilterGraphMedian(grA,cutSigma,medianY);
  if (graphTSA&&graphTSA->GetN()) SmoothGraph(graphTSA,deltaT);   
  TGraph* graphTSC= FilterGraphMedian(grC,cutSigma,medianY);
  if (graphTSC&&graphTSC->GetN()>0) SmoothGraph(graphTSC,deltaT);   
  delete grA; 
  delete grC;
  if (nA<kMinSectors) arrT->AddAt(0,72);
  else arrT->AddAt(graphTSA,72);
  if (nC<kMinSectors) arrT->AddAt(0,73);
  else arrT->AddAt(graphTSC,73);
}


void AliTPCcalibDButil::FilterTracks(Int_t run, Double_t cutSigma, TTreeSRedirector * const pcstream){
  //
  // Filter Drift velocity measurement using the tracks
  // 0.  remove outlyers - error based
  //     cutSigma      
  //
  //
  const Int_t kMinPoints=1;  // minimal number of points to define value
  TObjArray *arrT=AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);
  Double_t medianY=0;
  if (!arrT) return;
  for (Int_t i=0; i<arrT->GetEntries();i++){
    TGraphErrors *graph= dynamic_cast<TGraphErrors*>(arrT->At(i));
    if (!graph) continue;
    if (graph->GetN()<kMinPoints){
      delete graph;
      arrT->AddAt(0,i);
      continue;
    }
    TGraphErrors *graph2 = NULL;
    if(graph->GetN()<10) {
      graph2 = new TGraphErrors(graph->GetN(),graph->GetX(),graph->GetY(),graph->GetEX(),graph->GetEY()); 
      if (!graph2) {
        delete graph; arrT->AddAt(0,i); continue;
      }
    } 
    else {
      graph2= FilterGraphMedianErr(graph,cutSigma,medianY);
      if (!graph2) {
        delete graph; arrT->AddAt(0,i); continue;
      }
    }
    if (graph2->GetN()<1) {
      delete graph; arrT->AddAt(0,i); continue;
    }
    graph2->SetName(graph->GetName());
    graph2->SetTitle(graph->GetTitle());
    arrT->AddAt(graph2,i);
    if (pcstream){
      (*pcstream)<<"filterTracks"<<
	"run="<<run<<
	"isec="<<i<<
	"mY="<<medianY<<
	"graph.="<<graph<<
	"graph2.="<<graph2<<
	"\n";
    }
    delete graph;
  }
}





Double_t AliTPCcalibDButil::GetLaserTime0(Int_t run, Int_t timeStamp, Int_t deltaT, Int_t side){
  //
  //
  // get laser time offset 
  // median around timeStamp+-deltaT   
  // QA - chi2 needed for later usage - to be added
  //    - currently cut on error
  //
  Int_t kMinPoints=1;
  Double_t kMinDelay=0.01;
  Double_t kMinDelayErr=0.0001;

  TObjArray *array =AliTPCcalibDB::Instance()->GetTimeVdriftSplineRun(run);
  if (!array) return 0;
  TGraphErrors *tlaser=0;
  if (array){
    if (side==0) tlaser=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DELAY_LASER_ALL_A");
    if (side==1) tlaser=(TGraphErrors*)array->FindObject("GRAPH_MEAN_DELAY_LASER_ALL_C");
  }
  if (!tlaser) return 0;
  Int_t npoints0= tlaser->GetN();
  if (npoints0==0) return 0;
  Double_t *xlaser = new Double_t[npoints0];
  Double_t *ylaser = new Double_t[npoints0];
  Int_t npoints=0;
  for (Int_t i=0;i<npoints0;i++){
    //printf("%d\n",i);
    if (tlaser->GetY()[i]<=kMinDelay) continue; // filter zeros  
    if (tlaser->GetErrorY(i)>TMath::Abs(kMinDelayErr)) continue;
    xlaser[npoints]=tlaser->GetX()[npoints];
    ylaser[npoints]=tlaser->GetY()[npoints];
    npoints++;
  }
  //
  //
  Int_t index0=TMath::BinarySearch(npoints, xlaser, Double_t(timeStamp-deltaT))-1;
  Int_t index1=TMath::BinarySearch(npoints, xlaser, Double_t(timeStamp+deltaT))+1;
  //if (index1-index0 <kMinPoints) { index1+=kMinPoints; index0-=kMinPoints;}
  if (index0<0) index0=0;
  if (index1>=npoints-1) index1=npoints-1;
  if (index1-index0<kMinPoints) {
    delete [] ylaser;
    delete [] xlaser;
    return 0;
  }
  //
  //Double_t median = TMath::Median(index1-index0, &(ylaser[index0]));
    Double_t mean = TMath::Mean(index1-index0, &(ylaser[index0]));
  delete [] ylaser;
  delete [] xlaser;
  return mean;
}




void AliTPCcalibDButil::FilterGoofie(AliDCSSensorArray * goofieArray, Double_t deltaT, Double_t cutSigma, Double_t minVd, Double_t maxVd, TTreeSRedirector * const pcstream){
  //
  // Filter Goofie data
  // goofieArray - points will be filtered
  // deltaT      - smmothing time window 
  // cutSigma    - outler sigma cut in rms
  // minVn, maxVd- range absolute cut for variable vd/pt
  //             - to be tuned
  //
  // Ignore goofie if not enough points
  //
  const Int_t kMinPoints = 3;
  //

  TGraph *graphvd = goofieArray->GetSensorNum(2)->GetGraph();
  TGraph *graphan = goofieArray->GetSensorNum(8)->GetGraph();
  TGraph *graphaf = goofieArray->GetSensorNum(9)->GetGraph();
  TGraph *graphpt = goofieArray->GetSensorNum(15)->GetGraph();
  if (!graphvd) return;
  if (graphvd->GetN()<kMinPoints){
    delete graphvd;
    goofieArray->GetSensorNum(2)->SetGraph(0);
    return;
  }
  //
  // 1. Caluclate medians of critical variables
  //    drift velcocity
  //    P/T
  //    area near peak
  //    area far  peak
  //
  Double_t medianpt=0;
  Double_t medianvd=0, sigmavd=0;
  Double_t medianan=0;
  Double_t medianaf=0;
  Int_t    entries=graphvd->GetN();
  Double_t yvdn[10000];
  Int_t nvd=0;
  //
  for (Int_t ipoint=0; ipoint<entries; ipoint++){
    if (graphpt->GetY()[ipoint]<=0.0000001) continue;
    if (graphvd->GetY()[ipoint]/graphpt->GetY()[ipoint]<minVd) continue;
    if (graphvd->GetY()[ipoint]/graphpt->GetY()[ipoint]>maxVd) continue;
    yvdn[nvd++]=graphvd->GetY()[ipoint];
  }
  if (nvd<kMinPoints){
    delete graphvd;
    goofieArray->GetSensorNum(2)->SetGraph(0);
    return;
  }
  //
  Int_t nuni = TMath::Min(TMath::Nint(nvd*0.4+2), nvd-1);
  if (nuni>=kMinPoints){
    AliMathBase::EvaluateUni(nvd, yvdn, medianvd,sigmavd,nuni); 
  }else{
    medianvd = TMath::Median(nvd, yvdn);
  }
  
  TGraph * graphpt0 = AliTPCcalibDButil::FilterGraphMedianAbs(graphpt,10,medianpt);
  TGraph * graphpt1 = AliTPCcalibDButil::FilterGraphMedian(graphpt0,2,medianpt);
  TGraph * graphan0 = AliTPCcalibDButil::FilterGraphMedianAbs(graphan,10,medianan);
  TGraph * graphan1 = AliTPCcalibDButil::FilterGraphMedian(graphan0,2,medianan);
  TGraph * graphaf0 = AliTPCcalibDButil::FilterGraphMedianAbs(graphaf,10,medianaf);
  TGraph * graphaf1 = AliTPCcalibDButil::FilterGraphMedian(graphaf0,2,medianaf);
  delete graphpt0;
  delete graphpt1;
  delete graphan0;
  delete graphan1;
  delete graphaf0;
  delete graphaf1;
  //
  // 2. Make outlyer graph
  //
  Int_t nOK=0;
  TGraph graphOut(*graphvd);
  for (Int_t i=0; i<entries;i++){
    //
    Bool_t isOut=kFALSE;
    if (graphpt->GetY()[i]<=0.0000001) {  graphOut.GetY()[i]=1; continue;}
    if (graphvd->GetY()[i]/graphpt->GetY()[i]<minVd || graphvd->GetY()[i]/graphpt->GetY()[i]>maxVd) {  graphOut.GetY()[i]=1; continue;}
 
    if (TMath::Abs((graphvd->GetY()[i]/graphpt->GetY()[i])/medianvd-1.)<0.05) 
      isOut|=kTRUE;
    if (TMath::Abs(graphpt->GetY()[i]/medianpt-1.)>0.02) isOut|=kTRUE;
    if (TMath::Abs(graphan->GetY()[i]/medianan-1.)>0.2) isOut|=kTRUE;
    if (TMath::Abs(graphaf->GetY()[i]/medianaf-1.)>0.2) isOut|=kTRUE;
    graphOut.GetY()[i]= (isOut)?1:0;
    if (!isOut) nOK++;
  }
  if (nOK<kMinPoints) {
    delete graphvd;
    goofieArray->GetSensorNum(2)->SetGraph(0);
    return;
  } 
  //
  // 3. Filter out outlyers - and smooth 
  //
  TVectorF vmedianArray(goofieArray->NumSensors());
  TVectorF vrmsArray(goofieArray->NumSensors());
  Double_t xnew[10000];
  Double_t ynew[10000]; 
  TObjArray junk;
  junk.SetOwner(kTRUE);
  Bool_t isOK=kTRUE;
  //
  //
  for (Int_t isensor=0; isensor<goofieArray->NumSensors();isensor++){
    isOK=kTRUE;
    AliDCSSensor *sensor = goofieArray->GetSensorNum(isensor); 
    TGraph *graphOld=0, *graphNew=0, * graphNew0=0,*graphNew1=0,*graphNew2=0;
    //
    if (!sensor) continue;
    graphOld = sensor->GetGraph();
    if (graphOld) {
      sensor->SetGraph(0);
      Int_t nused=0;
      for (Int_t i=0;i<entries;i++){
	if (graphOut.GetY()[i]>0.5) continue;
	xnew[nused]=graphOld->GetX()[i];
	ynew[nused]=graphOld->GetY()[i];
	nused++;
      }
      graphNew = new TGraph(nused,xnew,ynew);
      junk.AddLast(graphNew);
      junk.AddLast(graphOld);      
      Double_t median=0;
      graphNew0  = AliTPCcalibDButil::FilterGraphMedian(graphNew,cutSigma,median);
      if (graphNew0!=0){
	junk.AddLast(graphNew0);
	graphNew1  = AliTPCcalibDButil::FilterGraphMedian(graphNew0,cutSigma,median);
	if (graphNew1!=0){
	  junk.AddLast(graphNew1);	  
	  graphNew2  = AliTPCcalibDButil::FilterGraphMedian(graphNew1,cutSigma,median);
	  if (graphNew2!=0) {
	    vrmsArray[isensor]   =TMath::RMS(graphNew2->GetN(),graphNew2->GetY());
	    AliTPCcalibDButil::SmoothGraph(graphNew2,deltaT);
	    AliTPCcalibDButil::SmoothGraph(graphNew2,deltaT);
	    AliTPCcalibDButil::SmoothGraph(graphNew2,deltaT);
	    //	    AliInfo(Form("%d\t%f\t%f\n",isensor, median,vrmsArray[isensor]));
	    vmedianArray[isensor]=median;
	    //
	  }
	}
      }
    }
    if (!graphOld) {  isOK=kFALSE; graphOld =&graphOut;}
    if (!graphNew0) { isOK=kFALSE; graphNew0=graphOld;}
    if (!graphNew1) { isOK=kFALSE; graphNew1=graphOld;}
    if (!graphNew2) { isOK=kFALSE; graphNew2=graphOld;}
    (*pcstream)<<"goofieA"<<
      Form("isOK_%d.=",isensor)<<isOK<<      
      Form("s_%d.=",isensor)<<sensor<<
      Form("gr_%d.=",isensor)<<graphOld<<
      Form("gr0_%d.=",isensor)<<graphNew0<<
      Form("gr1_%d.=",isensor)<<graphNew1<<
      Form("gr2_%d.=",isensor)<<graphNew2;
    if (isOK) sensor->SetGraph(graphNew2);
  }
  (*pcstream)<<"goofieA"<<
    "vmed.="<<&vmedianArray<<
    "vrms.="<<&vrmsArray<<
    "\n";
  junk.Delete();   // delete temoprary graphs
  
}





TMatrixD* AliTPCcalibDButil::MakeStatRelKalman(TObjArray * const array, Float_t minFraction, Int_t minStat, Float_t maxvd){
  //
  // Make a statistic matrix
  // Input parameters:
  //   array        - TObjArray of AliRelKalmanAlign 
  //   minFraction  - minimal ration of accepted tracks
  //   minStat      - minimal statistic (number of accepted tracks)
  //   maxvd        - maximal deviation for the 1
  // Output matrix:
  //    columns    - Mean, Median, RMS
  //    row        - parameter type (rotation[3], translation[3], drift[3])
  if (!array) return 0;
  if (array->GetEntries()<=0) return 0;
  //  Int_t entries = array->GetEntries();
  Int_t entriesFast = array->GetEntriesFast();
  TVectorD state(9);
  TVectorD *valArray[9];
  for (Int_t i=0; i<9; i++){
    valArray[i] = new TVectorD(entriesFast);
  }
  Int_t naccept=0;
  for (Int_t ikalman=0; ikalman<entriesFast; ikalman++){
    AliRelAlignerKalman * kalman = (AliRelAlignerKalman *) array->UncheckedAt(ikalman);
    if (!kalman) continue;
    if (TMath::Abs(kalman->GetTPCvdCorr()-1)>maxvd) continue;
    if (kalman->GetNUpdates()<minStat) continue;
    if (Float_t(kalman->GetNUpdates())/Float_t(kalman->GetNTracks())<minFraction) continue;
    kalman->GetState(state);
    for (Int_t ipar=0; ipar<9; ipar++)
      (*valArray[ipar])[naccept]=state[ipar];
    naccept++;
  }
  //if (naccept<2) return 0;
  if (naccept<1) return 0;
  TMatrixD *pstat=new TMatrixD(9,3);
  TMatrixD &stat=*pstat;
  for (Int_t ipar=0; ipar<9; ipar++){
    stat(ipar,0)=TMath::Mean(naccept, valArray[ipar]->GetMatrixArray());
    stat(ipar,1)=TMath::Median(naccept, valArray[ipar]->GetMatrixArray());
    stat(ipar,2)=TMath::RMS(naccept, valArray[ipar]->GetMatrixArray());
  }
  return pstat;
}


TObjArray *AliTPCcalibDButil::SmoothRelKalman(TObjArray * const array, const TMatrixD & stat, Bool_t direction, Float_t sigmaCut){
  //
  // Smooth the array of AliRelKalmanAlign - detector alignment and drift calibration)
  // Input:
  //   array     - input array
  //   stat      - mean parameters statistic
  //   direction - 
  //   sigmaCut  - maximal allowed deviation from mean in terms of RMS 
  if (!array) return 0;
  if (array->GetEntries()<=0) return 0;
  if (!(&stat)) return 0;
  // error increase in 1 hour
  const Double_t kerrsTime[9]={
    0.00001,  0.00001, 0.00001,
    0.001,    0.001,   0.001,
    0.002,  0.01,   0.001};
  //
  //
  Int_t entries = array->GetEntriesFast();
  TObjArray *sArray= new TObjArray(entries);
  AliRelAlignerKalman * sKalman =0;
  TVectorD state(9);
  for (Int_t i=0; i<entries; i++){
    Int_t index=(direction)? entries-i-1:i;
    AliRelAlignerKalman * kalman = (AliRelAlignerKalman *) array->UncheckedAt(index);
    if (!kalman) continue;
    Bool_t isOK=kTRUE;
    kalman->GetState(state);
    for (Int_t ipar=0; ipar<9; ipar++){
      if (TMath::Abs(state[ipar]-stat(ipar,1))>sigmaCut*stat(ipar,2)) isOK=kFALSE;
    }
    if (!sKalman &&isOK) {
      sKalman=new AliRelAlignerKalman(*kalman);
      sKalman->SetRejectOutliers(kFALSE);
      sKalman->SetRunNumber(kalman->GetRunNumber());
      sKalman->SetTimeStamp(kalman->GetTimeStamp());      
    }
    if (!sKalman) continue;
    Double_t deltaT=TMath::Abs(Int_t(kalman->GetTimeStamp())-Int_t(sKalman->GetTimeStamp()))/3600.;
    for (Int_t ipar=0; ipar<9; ipar++){
//       (*(sKalman->GetStateCov()))(6,6)+=deltaT*errvd*errvd;
//       (*(sKalman->GetStateCov()))(7,7)+=deltaT*errt0*errt0;
//       (*(sKalman->GetStateCov()))(8,8)+=deltaT*errvy*errvy; 
      (*(sKalman->GetStateCov()))(ipar,ipar)+=deltaT*kerrsTime[ipar]*kerrsTime[ipar];
    }
    sKalman->SetRunNumber(kalman->GetRunNumber());
    if (!isOK) sKalman->SetRunNumber(0);
    sArray->AddAt(new AliRelAlignerKalman(*sKalman),index);
    if (!isOK) continue;
    sKalman->SetRejectOutliers(kFALSE);
    sKalman->SetRunNumber(kalman->GetRunNumber());
    sKalman->SetTimeStamp(kalman->GetTimeStamp()); 
    sKalman->Merge(kalman);
    sArray->AddAt(new AliRelAlignerKalman(*sKalman),index);
    //sKalman->Print();
  }
  return sArray;
}

TObjArray *AliTPCcalibDButil::SmoothRelKalman(TObjArray * const arrayP, TObjArray * const arrayM){
  //
  // Merge 2 RelKalman arrays
  // Input:
  //   arrayP    - rel kalman in direction plus
  //   arrayM    - rel kalman in direction minus
  if (!arrayP) return 0;
  if (arrayP->GetEntries()<=0) return 0;
  if (!arrayM) return 0;
  if (arrayM->GetEntries()<=0) return 0;

  Int_t entries = arrayP->GetEntriesFast();
  TObjArray *array = new TObjArray(arrayP->GetEntriesFast()); 

  for (Int_t i=0; i<entries; i++){
     AliRelAlignerKalman * kalmanP = (AliRelAlignerKalman *) arrayP->UncheckedAt(i);
     AliRelAlignerKalman * kalmanM = (AliRelAlignerKalman *) arrayM->UncheckedAt(i);
     if (!kalmanP) continue;
     if (!kalmanM) continue;
  
     AliRelAlignerKalman *kalman = NULL;
     if(kalmanP->GetRunNumber() != 0 && kalmanM->GetRunNumber() != 0) {
       kalman = new AliRelAlignerKalman(*kalmanP);
       kalman->Merge(kalmanM);
     }
     else if (kalmanP->GetRunNumber() == 0) {
       kalman = new AliRelAlignerKalman(*kalmanM);
     }
     else if (kalmanM->GetRunNumber() == 0) {
       kalman = new AliRelAlignerKalman(*kalmanP);
     }
     else 
       continue;

     array->AddAt(kalman,i);
  }

  return array;
}
