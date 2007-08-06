#include <TPDGCode.h>
#include <TStyle.h>
#include "TMatrixD.h"
#include "TTreeStream.h"
#include "TF1.h"
#include "AliTPCParamSR.h"
#include "AliTPCClusterParam.h"
#include "AliTrackPointArray.h"
#include "TCint.h"
#include "AliTPCcalibTracksGain.h"
#include <TH1.h>
#include <TH3F.h>
#include <TLinearFitter.h>

//
// AliRoot includes
//
#include "AliMagF.h"
//
#include "AliTPCROC.h"
#include "AliTPCParamSR.h"
#include "AliTPCCalROC.h"
//
#include "AliTracker.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h" 
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"



ClassImp(AliTPCcalibTracksGain);

const Double_t AliTPCcalibTracksGain::fgkM = 25.;

AliTPCcalibTracksGain::AliTPCcalibTracksGain(const char* name, const char* title) :
   TNamed(name, title),
   fShortFitter(0),
   fMediumFitter(0),
   fLongFitter(0),
   fSqrtShortFitter(0),
   fSqrtMediumFitter(0),
   fSqrtLongFitter(0),
   fLogShortFitter(0),
   fLogMediumFitter(0),
   fLogLongFitter(0),
   fNShortClusters(0),
   fNMediumClusters(0),
   fNLongClusters(0),
   fTPCparam(0)
 {
   //
   // constructor
   //
   
   //TH1::AddDirectory(kFALSE);
   G__SetCatchException(0);
   //fDebugStream = new TTreeSRedirector("TPCSelectorDebug.root");

   fShortFitter = new TLinearFitter(6, "hyp5");
   fMediumFitter = new TLinearFitter(6, "hyp5");
   fLongFitter = new TLinearFitter(6, "hyp5");
   
   fSqrtShortFitter = new TLinearFitter(6, "hyp5");
   fSqrtMediumFitter = new TLinearFitter(6, "hyp5");
   fSqrtLongFitter = new TLinearFitter(6, "hyp5");
   
   fLogShortFitter = new TLinearFitter(6, "hyp5");
   fLogMediumFitter = new TLinearFitter(6, "hyp5");
   fLogLongFitter = new TLinearFitter(6, "hyp5");

   fShortFitter->StoreData(kFALSE);
   fMediumFitter->StoreData(kFALSE);
   fLongFitter->StoreData(kFALSE);
   fSqrtShortFitter->StoreData(kFALSE);
   fSqrtMediumFitter->StoreData(kFALSE);
   fSqrtLongFitter->StoreData(kFALSE);
   fLogShortFitter->StoreData(kFALSE);
   fLogMediumFitter->StoreData(kFALSE);
   fLogLongFitter->StoreData(kFALSE);
   
   fTPCparam = new AliTPCParamSR();
 }   

AliTPCcalibTracksGain::~AliTPCcalibTracksGain() {
   //
   // destructor
   //
   
   //if (fDebugStream) delete fDebugStream;
   
   if (fShortFitter) delete fShortFitter;
   if (fMediumFitter) delete fMediumFitter;
   if (fLongFitter) delete fLongFitter;

   if (fSqrtShortFitter) delete fSqrtShortFitter;
   if (fSqrtMediumFitter) delete fSqrtMediumFitter;
   if (fSqrtLongFitter) delete fSqrtLongFitter;

   if (fLogShortFitter) delete fLogShortFitter;
   if (fLogMediumFitter) delete fLogMediumFitter;
   if (fLogLongFitter) delete fLogLongFitter;
   
   if (fTPCparam) delete fTPCparam;
}


Bool_t AliTPCcalibTracksGain::AcceptTrack(AliTPCseed * track){
  //
  // Decides whether to accept a track or not.
  // Tracks are discarded, if due to edge effects, the number of clusters
  // is too low, the ratio of the number of clusters and the findable clusters is too low,
  // ...
  //
  
  const Int_t   kMinClusters  = 20;
  const Float_t kMinRatio     = 0.4;
  const Float_t kMax1pt       = 0.5;
  const Float_t kEdgeYXCutNoise    = 0.13;
  const Float_t kEdgeThetaCutNoise = 0.018;
  //
  // edge induced noise tracks - NEXT RELEASE will be removed during tracking
  if (TMath::Abs(track->GetY()/track->GetX())> kEdgeYXCutNoise)
    if (TMath::Abs(track->GetTgl())<kEdgeThetaCutNoise) { /*cerr << "[edge induced] " << flush;*/ return kFALSE; }
  
  //
  if (track->GetNumberOfClusters()<kMinClusters) { /*cerr << "[only " << track->GetNumberOfClusters() << " clusters] " << flush;*/ return kFALSE; }
  Float_t ratio = track->GetNumberOfClusters()/(track->GetNFoundable()+1.);
  if (ratio<kMinRatio) {/*cerr << "[ratio " << ratio << "] " << flush;*/ return kFALSE; }
  Float_t mpt = track->Get1Pt();
  if (TMath::Abs(mpt)>kMax1pt) { /*cerr << "[mpt " << mpt << "] " << flush;*/ return kFALSE; }
  //if (TMath::Abs(track->GetZ())>240.) return kFALSE;
  //if (TMath::Abs(track->GetZ())<10.) return kFALSE;
  //if (TMath::Abs(track->GetTgl())>0.03) return kFALSE;
  
  return kTRUE;
}

void AliTPCcalibTracksGain::AddTrack(AliTPCseed* seed) {
   //
   // The clusters making up the track (seed) are added to various fit functions.
   // See AddCluster(...) for more detail.
   //
   
   for (Int_t iCluster = 0; iCluster < 159; iCluster++) {
      AliTPCclusterMI* cluster = seed->GetClusterPointer(iCluster);
      if (cluster) AddCluster(cluster);
   }
}

void AliTPCcalibTracksGain::AddCluster(AliTPCclusterMI* cluster) {
   //
   // Adds cluster to the appropriate fitter for later analysis.
   // The charge used for the fit is the maximum charge for this specific cluster.
   // It is planned to add a switch to use the accumulated charge per cluster instead.
   // Depending on the pad size where the cluster is registered, the value will be put in
   // the appropriate fitter. Furthermore, for each pad size three different types of fitters
   // are used. The fit functions are the same for all fitters (parabolic functions), but the value
   // added to each fitter is different. The simple fitter gets the charge plugged in as is, the sqrt fitter
   // gets the square root of the charge, and the log fitter gets fgkM*(1+q/fgkM), where q is the original charge
   // and fgkM==25.
   //
   
   Double_t xx[5];
   
   // this block is for using the center of the region with the same pad size as origin for the fit function,
   // comment it out if you want the origin at lx=ly=0.
   Float_t centerPad[3] = {0};
   AliTPCROC* tpcROC = AliTPCROC::Instance();
   Int_t padType = GetPadType(cluster->GetX());
   Int_t IOROC = (padType == 0) ? 0 : tpcROC->GetNInnerSector();
   //tpcROC->GetPositionLocal(IOROC, tpcROC->GetNRows(IOROC)/2, tpcROC->GetNPads(IOROC, tpcROC->GetNRows(IOROC)/2)/2, centerPad);  // use this instead of the switch statement if you want to calculate the center of the ROC and not the center of the regions with the same pad size
   switch (padType) {
      case 0:           // short pads
         tpcROC->GetPositionLocal(IOROC, tpcROC->GetNRows(IOROC)/2, tpcROC->GetNPads(IOROC, tpcROC->GetNRows(IOROC)/2)/2, centerPad);
         break;
      case 1:           // medium pads
         tpcROC->GetPositionLocal(IOROC, 64/2, tpcROC->GetNPads(IOROC, 64/2)/2, centerPad);
         break;
      case 2:           // long pads
         tpcROC->GetPositionLocal(IOROC, 64+32/2, tpcROC->GetNPads(IOROC, 64+32/2)/2, centerPad);
         break;
   }
   xx[0] = cluster->GetX() - centerPad[0];
   xx[1] = cluster->GetY() - centerPad[1];
   
   //xx[0] = cluster->GetX();   // if you want the origin of the fit func at lx=ly=0
   //xx[1] = cluster->GetY();   // if you want the origin of the fit func at lx=ly=0
   xx[2] = xx[0] * xx[0];
   xx[3] = xx[1] * xx[1];
   xx[4] = xx[0] * xx[1];

   Double_t q = ((Double_t)(cluster->GetMax()));  // note: no normalization to pad size!
   Double_t sqrtQ = TMath::Sqrt(q);
   Double_t logQ = fgkM * TMath::Log(1 + q / fgkM);
   if (padType == 0) {
      fShortFitter->AddPoint(xx, q);
      fSqrtShortFitter->AddPoint(xx, sqrtQ);
      fLogShortFitter->AddPoint(xx, logQ);
      fNShortClusters++;
   } else if (padType == 1) {
      fMediumFitter->AddPoint(xx, q);
      fSqrtMediumFitter->AddPoint(xx, sqrtQ);
      fLogMediumFitter->AddPoint(xx, logQ);
      fNMediumClusters++;
   } else if (padType == 2) {
      fLongFitter->AddPoint(xx, q);
      fSqrtLongFitter->AddPoint(xx, sqrtQ);
      fLogLongFitter->AddPoint(xx, logQ);
      fNLongClusters++;
   }
}

Int_t AliTPCcalibTracksGain::Evaluate(UInt_t padType, TObjArray* fitParam, TObjArray* fitError, Double_t* redChi2, Bool_t robust) {
   //
   // Evaluate the tracks for obtaining the calibration information.
   // padType is 0 for evaluating the fitters for the short pads,
   // 1 for the medium pads and 2 for the long pads.
   // The elements of fitParam will be one TVectorD object for each fitter type (simple, sqrt, log fitter),
   // containing the fitted parameter values, the same is valid for the fitError, which
   // will contain the errors of course. redChi2 is an array with the entries being the
   // reduced chi^2 of each fitter type.
   // robust specifies wether the fitter's robust fitting mode shall be used (use with caution, it takes looooong!)
   //

   TObjArray fitters;
   Int_t NClusters = 0;
   switch (padType) {
      case 0:           // short pads
         fitters.Add(fShortFitter);
         fitters.Add(fSqrtShortFitter);
         fitters.Add(fLogShortFitter);
         NClusters = fNShortClusters;
         break;
      case 1:           // medium pads
         fitters.Add(fMediumFitter);
         fitters.Add(fSqrtMediumFitter);
         fitters.Add(fLogMediumFitter);
         NClusters = fNMediumClusters;
         break;
      case 2:           // long pads
         fitters.Add(fLongFitter);
         fitters.Add(fSqrtLongFitter);
         fitters.Add(fLogLongFitter);
         NClusters = fNLongClusters;
         break;
   }

   for (Int_t iFitter = 0; iFitter < fitters.GetEntries(); iFitter++) {
      TLinearFitter* fitter = (TLinearFitter*)fitters[iFitter];
      if (robust) fitter->EvalRobust();
      else fitter->Eval();
      
      if (redChi2) redChi2[iFitter] = fitter->GetChisquare()/(NClusters - 6);
      
      if (fitParam) {
         TVectorD* fitPar = new TVectorD(6);
         fitter->GetParameters(*fitPar);
         fitParam->Add(fitPar);
      }
      
      if (fitError) {
         TVectorD* fitErr = new TVectorD(6);
         fitter->GetErrors(*fitErr);
         *fitErr *= (redChi2) ? (TMath::Sqrt(redChi2[iFitter])) : (TMath::Sqrt(fitter->GetChisquare()/(NClusters - 6)));
         fitError->Add(fitErr);
      }
   }

   return NClusters;
}

AliTPCCalROC* AliTPCcalibTracksGain::CreateFitCalROC(UInt_t sector, UInt_t padType, TVectorD &fitParam, Int_t undoTransformation, Bool_t normalizeToPadSize) {
   //
   // This function is essentially a copy of AliTPCCalROC::CreateGlobalFitCalROC(...), with the
   // modifications, that the center of the region of same pad size is used as the origin
   // of the fit function instead of the center of the ROC.
   // The possibility of a linear fit is removed as well because it is not needed.
   // Only values for pads with the given pad size are calculated, the rest is 0.
   // Set undoTransformation to 0, 1 or 2 for undoing the transformation that was applied to the
   // charge values before they were put into the fitter (thus allowing comparison to the original
   // charge values). Use 0 for the simple fitter, 1 for the sqrt fitter, 2 for the log fitter.
   // Set it to -1 (or any other value) if you don't want any transformations undone (at the moment this is equivalent
   // with the transformation for the simple fitter (because no transformation is applied to the simple
   // fitter)).
   // If normalizeToPadSize is true, the values are normalized to the pad size.
   //
   
   Float_t dlx, dly;
   Float_t centerPad[3] = {0};
   Float_t localXY[3] = {0};
   AliTPCROC* tpcROC = AliTPCROC::Instance();
   if ((padType == 0 && sector >= tpcROC->GetNInnerSector()) || (padType > 0 && sector < tpcROC->GetNInnerSector()) || sector >= tpcROC->GetNSector())
      return 0;
   AliTPCCalROC* ROCfitted = new AliTPCCalROC(sector);
   //tpcROC->GetPositionLocal(sector, ROCfitted->GetNrows()/2, ROCfitted->GetNPads(ROCfitted->GetNrows()/2)/2, centerPad);  // use this instead of the switch statement if you want to calculate the center of the ROC and not the center of the regions with the same pad size
   UInt_t startRow;
   UInt_t endRow;
   switch (padType) {
      case 0:           // short pads
         startRow = 0;
         endRow = ROCfitted->GetNrows();
         tpcROC->GetPositionLocal(sector, endRow/2, ROCfitted->GetNPads(endRow/2)/2, centerPad);
         break;
      case 1:           // medium pads
         startRow = 0;
         endRow = 64;
         tpcROC->GetPositionLocal(sector, endRow/2, ROCfitted->GetNPads(endRow/2)/2, centerPad);
         break;
      case 2:           // long pads
         startRow = 64;
         endRow = ROCfitted->GetNrows();
         tpcROC->GetPositionLocal(sector, (endRow+startRow)/2, ROCfitted->GetNPads((endRow+startRow)/2)/2, centerPad);
         break;
   }
   
   Double_t value = 0;
   for (UInt_t irow = startRow; irow < endRow; irow++) {
      for (UInt_t ipad = 0; ipad < ROCfitted->GetNPads(irow); ipad++) {
         tpcROC->GetPositionLocal(sector, irow, ipad, localXY);   // calculate position localXY by pad and row number
         dlx = localXY[0] - centerPad[0];
         dly = localXY[1] - centerPad[1];
         value = fitParam[0] + fitParam[1]*dlx + fitParam[2]*dly + fitParam[3]*dlx*dlx + fitParam[4]*dly*dly + fitParam[5]*dlx*dly;
         switch (undoTransformation) {
            case  1: value = value * value; break;
            case  2: value = (TMath::Exp(value / fgkM) - 1) * fgkM; break;
            default: break;
         }
         if (normalizeToPadSize) value /= GetPadLength(localXY[0]);
         ROCfitted->SetValue(irow, ipad, value);
      }
   }
   return ROCfitted;
}

AliTPCCalROC* AliTPCcalibTracksGain::CreateCombinedCalROC(const AliTPCCalROC* roc1, const AliTPCCalROC* roc2) {
   //
   // Combines the medium pad size values of roc1 with the long pad size values of roc2 into a new
   // AliTPCCalROC. Returns a null pointer if any one of the ROCs is an IROC; issues a warning message
   // if the sectors of roc1 and roc2 don't match, but still continue and use the sector of roc1 as the
   // sector of the new ROC.
   //

   if (!roc1 || !roc2) return 0;
   if (roc1->GetSector() < fTPCparam->GetNInnerSector()) return 0;
   if (roc2->GetSector() < fTPCparam->GetNInnerSector()) return 0;
   if (roc1->GetSector() != roc2->GetSector()) Warning("CreateCombinedCalROC", "Sector number mismatch.");
   AliTPCCalROC* roc = new AliTPCCalROC(roc1->GetSector());
   
   for (UInt_t iRow = 0; iRow < 64; iRow++) {
      for (UInt_t iPad = 0; iPad < roc->GetNPads(iRow); iPad++)
         roc->SetValue(iRow, iPad, roc1->GetValue(iRow, iPad));
   }
   for (UInt_t iRow = 64; iRow < roc->GetNrows(); iRow++) {
      for (UInt_t iPad = 0; iPad < roc->GetNPads(iRow); iPad++)
         roc->SetValue(iRow, iPad, roc2->GetValue(iRow, iPad));
   }
   return roc;
}

void AliTPCcalibTracksGain::GetParameters(TVectorD &fitParam, UInt_t padType, UInt_t fitType) {
   //
   // Puts the fit parameters for the specified padType and fitType into the fitParam TVectorD
   // (which should contain 6 elements).
   // Note: The fitters have to be evaluated first!
   //    padType: 0 - short pads
   //             1 - medium pads
   //             2 - long pads
   //    fitType: 0 - simple fitter
   //             1 - sqrt fitter
   //             2 - log fitter
   //

   TLinearFitter* fitter = 0;
   switch (padType) {
      case 0: fitter = (fitType == 0) ? fShortFitter : ((fitType == 1) ? fSqrtShortFitter : fLogShortFitter);    break;
      case 1: fitter = (fitType == 0) ? fMediumFitter : ((fitType == 1) ? fSqrtMediumFitter : fLogMediumFitter); break;
      case 2: fitter = (fitType == 0) ? fLongFitter : ((fitType == 1) ? fSqrtLongFitter : fLogLongFitter);       break;
   }
   fitter->GetParameters(fitParam);
}

void AliTPCcalibTracksGain::GetErrors(TVectorD &fitError, UInt_t padType, UInt_t fitType) {
   //
   // Puts the fit parameter errors for the specified padType and fitType into the fitParam TVectorD
   // (which should contain 6 elements).
   // Note: The fitters have to be evaluated first!
   //

   TLinearFitter* fitter = 0;
   switch (padType) {
      case 0: fitter = (fitType == 0) ? fShortFitter : ((fitType == 1) ? fSqrtShortFitter : fLogShortFitter);    break;
      case 1: fitter = (fitType == 0) ? fMediumFitter : ((fitType == 1) ? fSqrtMediumFitter : fLogMediumFitter); break;
      case 2: fitter = (fitType == 0) ? fLongFitter : ((fitType == 1) ? fSqrtLongFitter : fLogLongFitter);       break;
   }
   fitter->GetErrors(fitError);
   fitError *= TMath::Sqrt(GetRedChi2(padType, fitType));
}

Double_t AliTPCcalibTracksGain::GetRedChi2(UInt_t padType, UInt_t fitType) {
   //
   // Returns the reduced chi^2 value for the specified padType and fitType.
   //

   TLinearFitter* fitter = 0;
   Int_t NClusters = 0;
   switch (padType) {
      case 0:
         fitter = (fitType == 0) ? fShortFitter : ((fitType == 1) ? fSqrtShortFitter : fLogShortFitter);
         NClusters = fNShortClusters;
         break;
      case 1:
         fitter = (fitType == 0) ? fMediumFitter : ((fitType == 1) ? fSqrtMediumFitter : fLogMediumFitter);
         NClusters = fNMediumClusters;
         break;
      case 2:
         fitter = (fitType == 0) ? fLongFitter : ((fitType == 1) ? fSqrtLongFitter : fLogLongFitter);
         NClusters = fNLongClusters;
         break;
   }
   return fitter->GetChisquare()/(NClusters - 6);
}

Double_t AliTPCcalibTracksGain::GetPadLength(Double_t lx) {
   //
   // The function returns 0 for an IROC, 1 for an OROC at medium pad size position,
   // 2 for an OROC at long pad size position, -1 if out of bounds.
   //

   Double_t irocLow = fTPCparam->GetPadRowRadiiLow(0) - fTPCparam->GetInnerPadPitchLength()/2;
   Double_t irocUp = fTPCparam->GetPadRowRadiiLow(fTPCparam->GetNRowLow()-1) + fTPCparam->GetInnerPadPitchLength()/2;
   Double_t orocLow1 = fTPCparam->GetPadRowRadiiUp(0) - fTPCparam->GetOuter1PadPitchLength()/2;
   Double_t orocUp1 = fTPCparam->GetPadRowRadiiUp(fTPCparam->GetNRowUp1()-1) + fTPCparam->GetOuter1PadPitchLength()/2;
   Double_t orocLow2 = fTPCparam->GetPadRowRadiiUp(fTPCparam->GetNRowUp1()) - fTPCparam->GetOuter2PadPitchLength()/2;
   Double_t orocUp2 = fTPCparam->GetPadRowRadiiUp(fTPCparam->GetNRowUp()-1) + fTPCparam->GetOuter2PadPitchLength()/2;
   
   // if IROC
   if (lx >= irocLow && lx <= irocUp) return 0.75;
   // if OROC medium pads
   if (lx >= orocLow1 && lx <= orocUp1) return 1.;
   // if OROC long pads
   if (lx >= orocLow2 && lx <= orocUp2) return 1.5;
   // if out of bounds
   return -1;
}

Int_t AliTPCcalibTracksGain::GetPadType(Double_t lx) {
   //
   // The function returns 0 for an IROC, 1 for an OROC at medium pad size position,
   // 2 for an OROC at long pad size position, -1 if out of bounds.
   //
   
   if (GetPadLength(lx) == 0.75) return 0;
   else if (GetPadLength(lx) == 1.) return 1;
   else if (GetPadLength(lx) == 1.5) return 2;
   return -1;
}
