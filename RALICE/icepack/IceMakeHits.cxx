/*******************************************************************************
 * Copyright(c) 2003, IceCube Experiment at the South Pole. All rights reserved.
 *
 * Author: The IceCube RALICE-based Offline Project.
 * Contributors are mentioned in the code where appropriate.
 *
 * Permission to use, copy, modify and distribute this software and its
 * documentation strictly for non-commercial purposes is hereby granted
 * without fee, provided that the above copyright notice appears in all
 * copies and that both the copyright notice and this permission notice
 * appear in the supporting documentation.
 * The authors make no claims about the suitability of this software for
 * any purpose. It is provided "as is" without express or implied warranty.
 *******************************************************************************/

// $Id$

///////////////////////////////////////////////////////////////////////////
// Class IceMakeHits
// TTask derived class to perform hit extraction from waveforms.
//
// In case an event has been rejected by an AliEventSelector (based) processor,
// this task (and its sub-tasks) is not executed.
//
// The code in this processor is based on the algorithms as developed by
// Nick van Eijndhoven and Garmt de Vries-Uiterweerd (Utrecht University, The Netherlands).
//
// Procedure applied for Amanda TWR data :
// ---------------------------------------
//
// 1) The waveform is fed to a TSpectrum object, and the peak locations 
//    are determined with the TSpectrum::Search() function.
//
// 2) The waveform is divided into regions corresponding to the peaks found by 
//    TSpectrum. The region boundary between two peaks is at the location of 
//    the minimum between the two peaks. 
//
// 3) For each region the "effective baseline" (used in the
//    evaluation of the leading edge value) is determined as :
//    effective baseline = fBasefracXXX * value at lower region boundary.
//    This takes into account the effect from the previous pulse.
//    For the first pulse, the effective baseline is equal to the overall 
//    baseline.
//
// 4) For each region, the point of steepest rise between the lower region
//    boundary and the peak location is determined. The tangent at this point
//    is extrapolated to the effective baseline. The point of intersection yields the
//    leading edge.
//
// 5) For each region the range of charge integration is determined as :
//    - Start of integration at the lower region boundary or at the leading edge,
//      whichever comes last;
//    - End of integration at the upper region boundary or at the point where the
//      signal drops below the overall baseline, whichever comes first.
//
// 6) For each region the integrated charge is determined as :
//    Sum over bins in integration range of (value in bin - overall baseline).
//
// 7) For each pulse the quality is evaluated by requiring that :
//    peak location - lower region boundary > lower region boundary - leading edge.
//    For a too shallow steepest rise, the leading edge value is unreliable, in 
//    which case the pulse is merged with the previous pulse. 
// 
// 8) Each pulse is checked for saturation and discarded if necessary.
//
// 9) TSpectrum needs a minimum number of bins for its Search function, otherwise
//    the clipping window is too large, which causes an error. If a waveform does
//    not contain enough bins, the following alternative approach is used : 
//    - A loop over all bins is performed.
//    - As soon as the signal exceeds a given threshold, a pulse is started.
//    - The pulse ends when the signal drops below the threshold again.
//    - While looping, the charge is integrated for each pulse.
//
// The defaults of the various parameters can be changed by the corresponding
// Set memberfunctions.
//
// Procedure applied for IceCube ATWD data :
// ---------------------------------------
//
// The procedure for ATWD data is essentially the same as for TWR data. The 
// parameters can be set independently.
// 
// Information about the actual parameter settings can be found in the event
// structure itself via the device named "IceMakeHits".
//
//--- Author: Nick van Eijndhoven and Garmt de Vries-Uiterweerd 15-jan-2007 Utrecht University
//- Modified: GdV $Date$ Utrecht University
///////////////////////////////////////////////////////////////////////////
 
#include "IceMakeHits.h"
#include "Riostream.h"

#include "TCanvas.h"

ClassImp(IceMakeHits) // Class implementation to enable ROOT I/O

IceMakeHits::IceMakeHits(const char* name,const char* title) : TTask(name,title)
{
// Default constructor.
 fEvt=0;
 // Parameters for Amanda TWR hit extraction
 fBasefracA=0.5;
 fSigmaA=1.5;
 fMaxPeaksA=10;
 fMinPulseHeightA=50;
 fThresholdA=0.2;
 // Parameters for InIce ATWD hit extraction
 fBasefracI=0.5;
 fSigmaI=2;
 fMaxPeaksI=10;
 fPeakAcceptanceLevelI=2;
 fStopIntegrationI=0.1;
 fMinPulseHeightI=1e-12;
 fThresholdI=0.2;
}
///////////////////////////////////////////////////////////////////////////
IceMakeHits::~IceMakeHits()
{
// Default destructor.
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::SetBasefracA(Float_t val)
{
// Set baseline fractional update for Amanda TWR extraction.
// The default as set in the constructor of this class is 0.5.
 fBasefracA=val;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::SetSigmaA(Float_t val)
{
// Set clipping window width for Amanda TWR extraction.
// The default as set in the constructor of this class is 1.5.
 fSigmaA=val;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::SetMaxPeaksA(Int_t val)
{
// Set maximum number of peaks in a waveform for Amanda TWR extraction.
// The default as set in the constructor of this class is 10.
 fMaxPeaksA=val;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::SetMinPulseHeightA(Float_t val)
{
// Set minimum required pulse height for Amanda TWR extraction.
// This is used only for narrow pulses that cannot be handled with TSpectrum.
// The default as set in the constructor of this class is 50.
 fMinPulseHeightA=val;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::SetThresholdA(Float_t val)
{
// Set threshold for use in analysis of narrow pulses for Amanda TWR extraction.
// A peak is assumed to start when the signal rises above threshold*maxval,
// where maxval is the maximum value found in the waveform.
// The default as set in the constructor of this class is 0.2.
 fThresholdA=val;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::SetBasefracI(Float_t val)
{
// Set baseline fractional update for InIce ATWD hit extraction.
// The default as set in the constructor of this class is 0.5.
 fBasefracI=val;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::SetSigmaI(Float_t val)
{
// Set clipping window width for InIce ATWD hit extraction.
// The default as set in the constructor of this class is 2.
 fSigmaI=val;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::SetMaxPeaksI(Int_t val)
{
// Set maximum number of peaks in a waveform for InIce ATWD hit extraction.
// The default as set in the constructor of this class is 10.
 fMaxPeaksI=val;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::SetPeakAcceptanceLevelI(Float_t val)
{
// Set the InIce minimum height of a peak above the baseline, in terms of baseline spread.
// The default as set in the constructor of this class is 2.
 fPeakAcceptanceLevelI=val;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::SetStopIntegrationI(Int_t val)
{
// Set fraction of the peak height at which charge integration ends in InIce ATWD hit extraction.
// The default as set in the constructor of this class is 0.1.
 fStopIntegrationI=val;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::SetMinPulseHeightI(Float_t val)
{
// Set minimum required pulse height for InIce ATWD hit extraction.
// This is used only for narrow pulses that cannot be handled with TSpectrum.
// The default as set in the constructor of this class is 50.
 fMinPulseHeightI=val;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::SetThresholdI(Float_t val)
{
// Set threshold for use in analysis of narrow pulses for InIce ATWD hit extraction.
// A peak is assumed to start when the signal rises above threshold*maxval,
// where maxval is the maximum value found in the waveform.
// The default as set in the constructor of this class is 0.2.
 fThresholdI=val;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::Exec(Option_t* opt)
{
// Implementation of the hit cleaning procedures.

 TString name=opt;
 AliJob* parent=(AliJob*)(gROOT->GetListOfTasks()->FindObject(name.Data()));

 if (!parent) return;

 fEvt=(IceEvent*)parent->GetObject("IceEvent");
 if (!fEvt) return;

 // Only process accepted events
 AliDevice* seldev=(AliDevice*)fEvt->GetDevice("AliEventSelector");
 if (seldev)
 {
  if (seldev->GetSignal("Select") < 0.1) return;
 }

 // Storage of the used parameters in the IceMakeHits device
 AliSignal params;
 params.SetNameTitle("IceMakeHits","IceMakeHits processor parameters");

 // Amanda hit extraction
 params.SetSlotName("BasefracA",1);
 params.SetSlotName("SigmaA",2);
 params.SetSlotName("MaxPeaksA",3);
 params.SetSlotName("MinPulseHeightA",4);
 params.SetSlotName("ThresholdA",5);
 params.SetSignal(fBasefracA,1);
 params.SetSignal(fSigmaA,2);
 params.SetSignal(fMaxPeaksA,3);
 params.SetSignal(fMinPulseHeightA,4);
 params.SetSignal(fThresholdA,5);

 // InIce hit extraction
 params.SetSlotName("BasefracI",1);
 params.SetSlotName("SigmaI",2);
 params.SetSlotName("MaxPeaksI",3);
 params.SetSlotName("MinPulseHeightI",4);
 params.SetSlotName("ThresholdI",5);
 params.SetSignal(fBasefracI,1);
 params.SetSignal(fSigmaI,2);
 params.SetSignal(fMaxPeaksI,3);
 params.SetSignal(fMinPulseHeightI,4);
 params.SetSignal(fThresholdI,5);

 fEvt->AddDevice(params);

 // Suppress (TSpectrum) warning messages
 gErrorIgnoreLevel=kError; // Only provide error messages

 Amanda();
 IceCube();
// InIce();
// IceTop();

 // Activate all messages
 gErrorIgnoreLevel=kWarning; // Error and warning messages

}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::Amanda()
{
// Hit extraction from the Amanda TWR data.

 // All Amanda OMs with a signal
 TObjArray* aoms=fEvt->GetDevices("IceAOM");
 if (!aoms) return;

 // Arrays for storing info
 Float_t* baseline=new Float_t[fMaxPeaksA];
 Int_t* lowend=new Int_t[fMaxPeaksA];
 Int_t* upend=new Int_t[fMaxPeaksA];
 Int_t* startcharge=new Int_t[fMaxPeaksA];
 Int_t* stopcharge=new Int_t[fMaxPeaksA];
 Int_t* status=new Int_t[fMaxPeaksA]; // 0=OK, 1=rejected, 2=saturation
 Float_t* leadingedge=new Float_t[fMaxPeaksA];
 Float_t* charge=new Float_t[fMaxPeaksA];
 Float_t* tot=new Float_t[fMaxPeaksA];

 // Some objects and variables we will need
 TH1F foo, diff;
 TSpectrum spec(fMaxPeaksA);
 Int_t nrIterations=(Int_t)(7*fSigmaA+0.5); // Number of iterations used in TSpectrum::SearchHighRes()
 Int_t npeaks=0, ibin=0, lookforsteepestuntilbin=0, steep=0;
 Float_t maxval=0, rise=0, rc=0, yyy=0;
 Bool_t pulsegoingon=false;
 Int_t* index=new Int_t[fMaxPeaksA];

 // OM, waveform and hit
 IceAOM* omx=0;
 TH1F* wf=0;
 AliSignal hit;
 hit.SetSlotName("ADC",1);
 hit.SetSlotName("LE",2);
 hit.SetSlotName("TOT",3);

 // Loop over all fired OMs and extract the hit info
 for (Int_t iom=0; iom<aoms->GetEntries(); iom++)
 {
  omx=(IceAOM*)aoms->At(iom);
  if (!omx) continue;
  // Remove all existing hits of this OM 
  omx->RemoveHits();
  // Reset (de)calibration functions to indicate uncalibrated data
  omx->SetCalFunction(0,"ADC");
  omx->SetCalFunction(0,"LE");
  omx->SetCalFunction(0,"TOT");
  omx->SetDecalFunction(0,"ADC");
  omx->SetDecalFunction(0,"LE");
  omx->SetDecalFunction(0,"TOT");
  // Should we skip OMs that we know from the dbase to have problems ?
////  if (omx->GetDeadValue("ADC") || omx->GetDeadValue("LE") || omx->GetDeadValue("TOT")) continue;

  // Investigate all waveforms for this OM
  for (Int_t iwf=1; iwf<=omx->GetNwaveforms(); iwf++)
  {
   wf=omx->GetWaveform(iwf);
   if (!wf) continue; 
   maxval=wf->GetMaximum();
   // Check if clipping window is not too large
   if(wf->GetNbinsX() > 2*nrIterations+1)
   {
    // Find peaks with TSpectrum
    npeaks=spec.Search(wf,fSigmaA,"goff");
    // Discard waveform if no or too many peaks found
    if(npeaks<1 || npeaks>fMaxPeaksA) continue;

    // Get differential of WF
    diff=*wf;
    for(ibin=2;ibin<diff.GetNbinsX();ibin++)
    {
     diff.SetBinContent(ibin,wf->GetBinContent(ibin)-wf->GetBinContent(ibin-1));
    }
    diff.SetBinContent(1,0);
    // Set baseline and lower end for first peak,
    baseline[0]=0;
    lowend[0]=1;

    // Sort peaks in time
    TMath::Sort(npeaks,spec.GetPositionX(),index,false);
    // For each of the peaks,
    for(Int_t ipeak=0; ipeak<npeaks; ipeak++)
    {
     // Find baseline and region around peak
     foo=*wf;
     // (Second and later peaks: lower edge = upper edge previous peak,
     // baseline is average of previous baseline and minimum value between two
     // peaks)
     if(ipeak>0)
     {
      lowend[ipeak]=upend[ipeak-1]+1;
      baseline[ipeak]=fBasefracA*foo.GetBinContent(lowend[ipeak]);
     }
     // (Upper edge range is minimum between this and next peak)
     if(ipeak<npeaks-1)
     {
      foo.SetAxisRange(spec.GetPositionX()[index[ipeak]],spec.GetPositionX()[index[ipeak+1]]);
      upend[ipeak]=foo.GetMinimumBin()-1;
     }
     // (Last peak: upper edge is end of histo)
     else
     {
      upend[ipeak]=wf->GetNbinsX();
     }
     // Find steepest rise
     lookforsteepestuntilbin=wf->FindBin(spec.GetPositionX()[index[ipeak]]);
     foo=diff;
     // Look for steepest rise between lower edge and peak position
     foo.SetAxisRange(wf->GetBinCenter(lowend[ipeak]),wf->GetBinCenter(lookforsteepestuntilbin));
     steep=foo.GetMaximumBin();
     rise=foo.GetBinContent(steep);

     // Extrapolate tangent to find leading edge
     yyy=wf->GetBinContent(steep)-baseline[ipeak];
     rc=rise/foo.GetBinWidth(steep);
     if(rc>0) leadingedge[ipeak]=wf->GetBinCenter(steep)-yyy/rc; else leadingedge[ipeak]=0;

     // Determine peak status
     status[ipeak]=0;
     // Check for saturation
     if(rc<0.1 && wf->GetBinContent(wf->FindBin(spec.GetPositionX()[index[ipeak]])) == maxval)
     {
      status[ipeak]=2;
     }
     // Check quality: LE should not be too far below lower edge
     // Otherwise, ignore this peak and set baseline back to what it was
     else if(wf->GetBinLowEdge(lowend[ipeak]) - leadingedge[ipeak] > spec.GetPositionX()[index[ipeak]] - wf->GetBinLowEdge(lowend[ipeak]))
     {
      status[ipeak]=1;
      if(ipeak>0) baseline[ipeak]=baseline[ipeak-1];
     }
 
     // Start charge integration at LE, or at lower edge of range
     startcharge[ipeak]=wf->FindBin(leadingedge[ipeak]);
     if(lowend[ipeak]>startcharge[ipeak]) startcharge[ipeak]=lowend[ipeak];
 
     // Integrate charge until pulse drop below baseline, or else until edge of range
     stopcharge[ipeak]=upend[ipeak];
     for(ibin=wf->FindBin(spec.GetPositionX()[index[ipeak]]); ibin<=upend[ipeak]; ibin++)
     {
      if(wf->GetBinContent(ibin)<0)
      {
       stopcharge[ipeak]=ibin-1;
       break;
      }
     }

     // Determine time over threshold
     tot[ipeak]=wf->GetBinLowEdge(stopcharge[ipeak]+1)-wf->GetBinLowEdge(startcharge[ipeak]);
 
     // Determine charge
     charge[ipeak]=0;
     for(ibin=startcharge[ipeak]; ibin<=stopcharge[ipeak]; ibin++)
     {
      charge[ipeak]+=wf->GetBinContent(ibin);
     }

    } // end loop over peaks
 
    // Check all peaks, from latest to earliest
    for(int ipeak=npeaks-1; ipeak>=0; ipeak--)
    {

     // If this peak was rejected, add charge and TOT to previous peak (if there is one)
     if(status[ipeak]==1 && ipeak>0)
     {
      charge[ipeak-1]+=charge[ipeak];
      charge[ipeak]=0;
      tot[ipeak-1]+=tot[ipeak];
      tot[ipeak]=0;
     }

     // If this peak is OK, add hit info
     if(status[ipeak]==0)
     {
      hit.Reset();
      hit.SetSignal(charge[ipeak],"ADC");
      hit.SetSignal(leadingedge[ipeak],"LE");
      hit.SetSignal(tot[ipeak],"TOT");
      omx->AddHit(hit);
     }

    } // end loop over peaks
   // If number of bins too small, use different method
   }
   else
   {
    // If maximum value high enough to suspect presence of peak,
    if(maxval>fMinPulseHeightA)
    {
     // Loop over bins
     pulsegoingon=false;
     npeaks=0;
     for(ibin=1; ibin<=wf->GetNbinsX(); ibin++)
     {
      // If bin content above threshold, start pulse
      if(wf->GetBinContent(ibin)>fThresholdA*maxval){
       if(!pulsegoingon)
       {
	// Pulse starts here
        pulsegoingon=true;
        leadingedge[npeaks]=wf->GetBinLowEdge(ibin);
        charge[npeaks]=wf->GetBinContent(ibin);
       }
       else
       {
	// Pulse continues       
        charge[npeaks]+=wf->GetBinContent(ibin);
       }
      }
      else
      {
       if(pulsegoingon)
       {
	// Pulse ends here
        tot[npeaks]=wf->GetBinLowEdge(ibin)-leadingedge[npeaks];

	// Store pulse information
        hit.Reset();
        hit.SetSignal(charge[npeaks],"ADC");
        hit.SetSignal(leadingedge[npeaks],"LE");
        hit.SetSignal(tot[npeaks],"TOT");
        omx->AddHit(hit);

        // Get ready for next pulse
	pulsegoingon=false;
        npeaks++;
       }
      }

     } // End of loop over bins
    } 

   } // End of alternative method for narrow pulses

  } // End of WF loop
 } // End of OM loop

 // Clean up
 delete[] baseline;
 delete[] lowend;
 delete[] upend;
 delete[] startcharge;
 delete[] stopcharge;
 delete[] status;
 delete[] leadingedge;
 delete[] charge;
 delete[] tot;
 delete[] index;
}
///////////////////////////////////////////////////////////////////////////
void IceMakeHits::IceCube()
{
// Hit extraction from IceCube ATWD data.

 // All InIce DOMs with a signal
 TObjArray* idoms=fEvt->GetDevices("IceDOM");
 if (!idoms) return;

 // Arrays for storing info
 Float_t* baseline=new Float_t[fMaxPeaksI];
 Int_t* lowend=new Int_t[fMaxPeaksI];
 Int_t* upend=new Int_t[fMaxPeaksI];
 Int_t* startcharge=new Int_t[fMaxPeaksI];
 Int_t* stopcharge=new Int_t[fMaxPeaksI];
 Int_t* status=new Int_t[fMaxPeaksI]; // 0=OK, 1=rejected, 2=saturation
 Float_t* leadingedge=new Float_t[fMaxPeaksI];
 Float_t* charge=new Float_t[fMaxPeaksI];
 Float_t* tot=new Float_t[fMaxPeaksI];

 // Some objects and variables we will need
 TH1F foo, diff;
 TSpectrum spec(fMaxPeaksI);
 Int_t nrIterations=(Int_t)(7*fSigmaI+0.5); // Number of iterations used in TSpectrum::SearchHighRes()
 Int_t npeaks=0, ibin=0, lookforsteepestuntilbin=0, steep=0;
 Float_t maxval=0, rise=0, rc=0, yyy=0, median=0, spread=0;
 Bool_t pulsegoingon=false;
 Int_t* index=new Int_t[fMaxPeaksI];
 TString slotname;

 // OM, waveform and hit
 IceIDOM* omx=0;
 TH1F* wf=0;
 AliSignal hit;
 hit.SetSlotName("ADC",1);
 hit.SetSlotName("LE",2);
 hit.SetSlotName("TOT",3);

 // Loop over all fired OMs and extract the hit info
 for (Int_t iom=0; iom<idoms->GetEntries(); iom++)
 {
  omx=(IceIDOM*)idoms->At(iom);
  if (!omx) continue;
  // Remove all existing hits of this OM 
  omx->RemoveHits();
  // Reset (de)calibration functions to indicate uncalibrated data
  omx->SetCalFunction(0,"ADC");
  omx->SetCalFunction(0,"LE");
  omx->SetCalFunction(0,"TOT");
  omx->SetDecalFunction(0,"ADC");
  omx->SetDecalFunction(0,"LE");
  omx->SetDecalFunction(0,"TOT");
  // Should we skip OMs that we know from the dbase to have problems ?
////  if (omx->GetDeadValue("ADC") || omx->GetDeadValue("LE") || omx->GetDeadValue("TOT")) continue;

  // Investigate all waveforms for this OM
  for (Int_t iwf=1; iwf<=omx->GetNwaveforms(); iwf++)
  {
   wf=omx->GetWaveform(iwf);
   if (!wf) continue; 
   maxval=wf->GetMaximum();
   // Check if clipping window is not too large
   if(wf->GetNbinsX() > 2*nrIterations+1)
   {
    // Find peaks with TSpectrum
    npeaks=spec.Search(wf,fSigmaI,"goff");
    // Discard waveform if no or too many peaks found
    if(npeaks<1 || npeaks>fMaxPeaksI) continue;

    // Get differential of WF
    diff=*wf;
    for(ibin=2;ibin<diff.GetNbinsX();ibin++)
    {
     diff.SetBinContent(ibin,wf->GetBinContent(ibin)-wf->GetBinContent(ibin-1));
    }
    diff.SetBinContent(1,0);
    // Set baseline and lower end for first peak,
    baseline[0]=0;
    lowend[0]=1;

    // Sort peaks in time
    TMath::Sort(npeaks,spec.GetPositionX(),index,false);
    // For each of the peaks,
    for(Int_t ipeak=0; ipeak<npeaks; ipeak++)
    {
     // Find baseline and region around peak
     foo=*wf;
     // (Second and later peaks: lower edge = upper edge previous peak,
     // baseline is average of previous baseline and minimum value between two
     // peaks)
     if(ipeak>0)
     {
      lowend[ipeak]=upend[ipeak-1]+1;
      baseline[ipeak]=fBasefracI*foo.GetBinContent(lowend[ipeak]);
     }
     // (Upper edge range is minimum between this and next peak)
     if(ipeak<npeaks-1)
     {
      foo.SetAxisRange(spec.GetPositionX()[index[ipeak]],spec.GetPositionX()[index[ipeak+1]]);
      upend[ipeak]=foo.GetMinimumBin()-1;
     }
     // (Last peak: upper edge is end of histo)
     else
     {
      upend[ipeak]=wf->GetNbinsX();
     }
     // Find steepest rise
     lookforsteepestuntilbin=wf->FindBin(spec.GetPositionX()[index[ipeak]]);
     foo=diff;
     // Look for steepest rise between lower edge and peak position
     foo.SetAxisRange(wf->GetBinCenter(lowend[ipeak]),wf->GetBinCenter(lookforsteepestuntilbin));
     steep=foo.GetMaximumBin();
     rise=foo.GetBinContent(steep);

     // Extrapolate tangent to find leading edge
     yyy=wf->GetBinContent(steep)-baseline[ipeak];
     rc=rise/foo.GetBinWidth(steep);
     if(rc>0) leadingedge[ipeak]=wf->GetBinCenter(steep)-yyy/rc; else leadingedge[ipeak]=0;

     // Determine peak status
     status[ipeak]=0;
     // Check for saturation
     if(rc<0.1 && wf->GetBinContent(wf->FindBin(spec.GetPositionX()[index[ipeak]])) == maxval)
     {
      status[ipeak]=2;
     }
     // Check quality: LE should not be too far below lower edge
     // Otherwise, ignore this peak and set baseline back to what it was
     else if(wf->GetBinLowEdge(lowend[ipeak]) - leadingedge[ipeak] > spec.GetPositionX()[index[ipeak]] - wf->GetBinLowEdge(lowend[ipeak]))
     {
      status[ipeak]=1;
      if(ipeak>0) baseline[ipeak]=baseline[ipeak-1];
     }
     // Check if peak is high enough
     slotname="BASELINE-WF";
     slotname+=iwf;
     median=omx->GetSignal(slotname);
     spread=omx->GetSignalError(slotname);
     if(wf->GetBinContent(lookforsteepestuntilbin) < median+fPeakAcceptanceLevelI*spread){
      status[ipeak]=3;
     }
 
     // Start charge integration at LE, or at lower edge of range
     startcharge[ipeak]=wf->FindBin(leadingedge[ipeak]);
     if(lowend[ipeak]>startcharge[ipeak]) startcharge[ipeak]=lowend[ipeak];
 
     // Integrate charge until pulse drop below threshold, or else until edge of range
     stopcharge[ipeak]=upend[ipeak];
     for(ibin=wf->FindBin(spec.GetPositionX()[index[ipeak]]); ibin<=upend[ipeak]; ibin++)
     {
      if(wf->GetBinContent(ibin)<fStopIntegrationI*wf->GetBinContent(lookforsteepestuntilbin))
      {
       stopcharge[ipeak]=ibin-1;
       break;
      }
     }

     // Determine time over threshold
     tot[ipeak]=wf->GetBinLowEdge(stopcharge[ipeak]+1)-wf->GetBinLowEdge(startcharge[ipeak]);
 
     // Determine charge
     charge[ipeak]=0;
     for(ibin=startcharge[ipeak]; ibin<=stopcharge[ipeak]; ibin++)
     {
      charge[ipeak]+=wf->GetBinContent(ibin);
     }

    } // end loop over peaks
 
    // Check all peaks, from latest to earliest
    for(int ipeak=npeaks-1; ipeak>=0; ipeak--)
    {
     // If this peak was rejected, add charge and TOT to previous peak (if there is one)
     if(status[ipeak]==1 && ipeak>0)
     {
      charge[ipeak-1]+=charge[ipeak];
      charge[ipeak]=0;
      tot[ipeak-1]+=tot[ipeak];
      tot[ipeak]=0;
     }

     // If this peak is OK, add hit info
     if(status[ipeak]==0)
     {
      hit.Reset();
      hit.SetSignal(charge[ipeak],"ADC");
      hit.SetSignal(leadingedge[ipeak],"LE");
      hit.SetSignal(tot[ipeak],"TOT");
      omx->AddHit(hit);
     }

    } // end loop over peaks
   // If number of bins too small, use different method
   }
   else
   {
    // If maximum value high enough to suspect presence of peak,
    if(maxval>fMinPulseHeightI)
    {
     // Loop over bins
     pulsegoingon=false;
     npeaks=0;
     for(ibin=1; ibin<=wf->GetNbinsX(); ibin++)
     {
      // If bin content above threshold, start pulse
      if(wf->GetBinContent(ibin)>fThresholdI*maxval){
       if(!pulsegoingon)
       {
	// Pulse starts here
        pulsegoingon=true;
        leadingedge[npeaks]=wf->GetBinLowEdge(ibin);
        charge[npeaks]=wf->GetBinContent(ibin);
       }
       else
       {
	// Pulse continues       
        charge[npeaks]+=wf->GetBinContent(ibin);
       }
      }
      else
      {
       if(pulsegoingon)
       {
	// Pulse ends here
        tot[npeaks]=wf->GetBinLowEdge(ibin)-leadingedge[npeaks];

	// Store pulse information
        hit.Reset();
        hit.SetSignal(charge[npeaks],"ADC");
        hit.SetSignal(leadingedge[npeaks],"LE");
        hit.SetSignal(tot[npeaks],"TOT");
        omx->AddHit(hit);

        // Get ready for next pulse
	pulsegoingon=false;
        npeaks++;
       }
      }

     } // End of loop over bins
    } 

   } // End of alternative method for narrow pulses

  } // End of WF loop
 } // End of OM loop

 // Clean up
 delete[] baseline;
 delete[] lowend;
 delete[] upend;
 delete[] startcharge;
 delete[] stopcharge;
 delete[] status;
 delete[] leadingedge;
 delete[] charge;
 delete[] tot;
 delete[] index;
}
///////////////////////////////////////////////////////////////////////////
