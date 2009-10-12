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

#include <AliDCSSensorArray.h>
#include <AliDCSSensor.h>
#include "AliTPCcalibDB.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "AliTPCROC.h"
#include "AliTPCmapper.h"
#include "AliTPCParam.h"
#include "AliTPCCalibRaw.h"

#include "AliTPCcalibDButil.h"
#include "AliTPCPreprocessorOnline.h"

ClassImp(AliTPCcalibDButil)
AliTPCcalibDButil::AliTPCcalibDButil() :
  TObject(),
  fCalibDB(AliTPCcalibDB::Instance()),
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
  fRefPadNoise(0x0),
  fRefPedestals(0x0),
  fRefPulserTmean(0x0),
  fRefPulserTrms(0x0),
  fRefPulserQmean(0x0),
  fRefPulserOutlier(new AliTPCCalPad("RefPulserOutliers","RefPulserOutliers")),
  fRefCETmean(0x0),
  fRefCETrms(0x0),
  fRefCEQmean(0x0),
  fRefALTROMasked(0x0),
  fRefCalibRaw(0x0),
  fGoofieArray(0x0),
  fMapper(new AliTPCmapper(0x0)),
  fNpulserOutliers(-1),
  fIrocTimeOffset(0),
  fCETmaxLimitAbs(1.5),
  fPulTmaxLimitAbs(1.5),
  fPulQmaxLimitAbs(5),
  fPulQminLimit(11)
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
  if (fRefPulserTmean) delete fRefPulserTmean;
  if (fRefPulserTrms) delete fRefPulserTrms;
  if (fRefPulserQmean) delete fRefPulserQmean;
  if (fRefCETmean) delete fRefCETmean;
  if (fRefCETrms) delete fRefCETrms;
  if (fRefCEQmean) delete fRefCEQmean;
  if (fRefALTROMasked) delete fRefALTROMasked;
  if (fRefCalibRaw) delete fRefCalibRaw;
    
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::UpdateFromCalibDB()
{
  //
  // Update pointers from calibDB
  //
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
  UpdatePulserOutlierMap();
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::ProcessCEdata(const char* fitFormula, TVectorD &fitResultsA, TVectorD &fitResultsC,
                                      Int_t &noutliersCE, Double_t & chi2A, Double_t &chi2C, AliTPCCalPad *outCE)
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
        if (valTmean==0) {
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
        if (gr->GetY()[ipoint]>500 && gr->GetY()[ipoint]<1000 ){
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
                      Int_t &nonMaskedZero)
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
        if (noiseVal==0) {
          ++nonMaskedZero;
          continue;
        }
        //check for nan
        if ( !(noiseVal<10000000) ){
          printf ("Warning: nan detected in (sec,row,pad - val): %02d,%02d,%03d - %.1f\n",isec,irow,ipad,noiseVal);
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
//       printf ("i: %d - m: %.3f, c: %.0f, r: %.3f\n",i,vNoiseMean[i],c[i],vNoiseRMS[i]);
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
    Float_t pt_mean=ptROC->GetMean(oROC);
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
        if (TMath::Abs(pt-pt_mean)>1) ++npadsOutOneTB;
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
  //
  //
  PulserOutlierMap(fPulserOutlier,fPulserTmean, fPulserQmean);
}
//_____________________________________________________________________________________
void AliTPCcalibDButil::UpdateRefPulserOutlierMap()
{
  //
  //
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
      if (mean==0) mean=rocPulTmean->GetMean();
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
            if (time==0) time=mean;
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
      if (mean==0) mean=rocPulTmean->GetMean();
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
Float_t AliTPCcalibDButil::GetMeanAltro(const AliTPCCalROC *roc, const Int_t row, const Int_t pad, AliTPCCalROC *rocOut)
{
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




AliTPCCalPad *AliTPCcalibDButil::CreateCEOutlyerMap( Int_t & noutliersCE, AliTPCCalPad *ceOut, Float_t minSignal, Float_t cutTrmsMin,  Float_t cutTrmsMax, Float_t cutMaxDistT){
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
    if (fALTROMasked) rocMasked= fALTROMasked->GetCalROC(iroc);
    AliTPCCalROC *rocOut       = out->GetCalROC(iroc);
    AliTPCCalROC *rocCEQ       = fCEQmean->GetCalROC(iroc);
    AliTPCCalROC *rocCETrms    = fCETrms->GetCalROC(iroc);
    Double_t trocMedian        = rocData->GetMedian();
    //
    if (!rocData) {
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
        if (valTmean==0) {
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


AliTPCCalPad *AliTPCcalibDButil::CreatePulserOutlyerMap(Int_t &noutliersPulser, AliTPCCalPad *pulserOut,Float_t cutTime, Float_t cutnRMSQ, Float_t cutnRMSrms){
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
      Int_t isOut=0;
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

