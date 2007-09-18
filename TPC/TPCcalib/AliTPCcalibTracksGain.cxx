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
#include <TTreeStream.h>
#include <TFile.h>

//
// AliRoot includes
//
#include "AliMagF.h"
#include "AliMathBase.h"
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


/*
  
TFile f("TPCCalibTracksGain.root")

gSystem->Load("libPWG1.so")
AliTreeDraw comp
comp.SetTree(dEdx)
Double_t chi2;
TVectorD vec(3)
TMatrixD mat(3,3)
TString * str = comp.FitPlane("Cl.fQ/dedxQ.fElements[0]","Cl.fY++Cl.fX","Cl.fDetector<36",chi2,vec,mat)


*/


ClassImp(AliTPCcalibTracksGain);

const Double_t AliTPCcalibTracksGain::fgkM = 25.;

AliTPCcalibTracksGain::AliTPCcalibTracksGain(const char* name, const char* title) :
   TNamed(name, title),
   fDebugStream(0),
   fShortFitter(0),
   fMediumFitter(0),
   fLongFitter(0),
   fSqrtShortFitter(0),
   fSqrtMediumFitter(0),
   fSqrtLongFitter(0),
   fLogShortFitter(0),
   fLogMediumFitter(0),
   fLogLongFitter(0),
   fTPCparam(0)
 {
   //
   // constructor
   //
   
   //TH1::AddDirectory(kFALSE);
   G__SetCatchException(0);

   for (UInt_t iSegment = 0; iSegment < 36; iSegment++) {
      fNShortClusters[iSegment] = 0;
      fNMediumClusters[iSegment] = 0;
      fNLongClusters[iSegment] = 0;
   }

   fShortFitter = new TObjArray(36);
   fMediumFitter = new TObjArray(36);
   fLongFitter = new TObjArray(36);
   
   fSqrtShortFitter = new TObjArray(36);
   fSqrtMediumFitter = new TObjArray(36);
   fSqrtLongFitter = new TObjArray(36);
   
   fLogShortFitter = new TObjArray(36);
   fLogMediumFitter = new TObjArray(36);
   fLogLongFitter = new TObjArray(36);

   for (UInt_t iSegment = 0; iSegment < 36; iSegment++) {   
      fShortFitter->Add(new TLinearFitter(6, "hyp5"));
      fMediumFitter->Add(new TLinearFitter(6, "hyp5"));
      fLongFitter->Add(new TLinearFitter(6, "hyp5"));
      
      fSqrtShortFitter->Add(new TLinearFitter(6, "hyp5"));
      fSqrtMediumFitter->Add(new TLinearFitter(6, "hyp5"));
      fSqrtLongFitter->Add(new TLinearFitter(6, "hyp5"));
      
      fLogShortFitter->Add(new TLinearFitter(6, "hyp5"));
      fLogMediumFitter->Add(new TLinearFitter(6, "hyp5"));
      fLogLongFitter->Add(new TLinearFitter(6, "hyp5"));
   
      ((TLinearFitter*)(fShortFitter->At(iSegment)))->StoreData(kFALSE);
      ((TLinearFitter*)(fMediumFitter->At(iSegment)))->StoreData(kFALSE);
      ((TLinearFitter*)(fLongFitter->At(iSegment)))->StoreData(kFALSE);
      ((TLinearFitter*)(fSqrtShortFitter->At(iSegment)))->StoreData(kFALSE);
      ((TLinearFitter*)(fSqrtMediumFitter->At(iSegment)))->StoreData(kFALSE);
      ((TLinearFitter*)(fSqrtLongFitter->At(iSegment)))->StoreData(kFALSE);
      ((TLinearFitter*)(fLogShortFitter->At(iSegment)))->StoreData(kFALSE);
      ((TLinearFitter*)(fLogMediumFitter->At(iSegment)))->StoreData(kFALSE);
      ((TLinearFitter*)(fLogLongFitter->At(iSegment)))->StoreData(kFALSE);
   }
   
   fTPCparam = new AliTPCParamSR();
 }   

AliTPCcalibTracksGain::~AliTPCcalibTracksGain() {
   //
   // destructor
   //
   
   //if (fDebugStream) delete fDebugStream;
   
   if (fShortFitter) { fShortFitter->Delete(); delete fShortFitter; }
   if (fMediumFitter) { fMediumFitter->Delete(); delete fMediumFitter; }
   if (fLongFitter) { fLongFitter->Delete(); delete fLongFitter; }

   if (fSqrtShortFitter) { fSqrtShortFitter->Delete(); delete fSqrtShortFitter; }
   if (fSqrtMediumFitter) { fSqrtMediumFitter->Delete(); delete fSqrtMediumFitter; }
   if (fSqrtLongFitter) { fSqrtLongFitter->Delete(); delete fSqrtLongFitter; }

   if (fLogShortFitter) { fLogShortFitter->Delete(); delete fLogShortFitter; }
   if (fLogMediumFitter) { fLogMediumFitter->Delete(); delete fLogMediumFitter; }
   if (fLogLongFitter) { fLogLongFitter->Delete(); delete fLogLongFitter; }
   //fDebugStream->GetFile()->Close();
   printf("Deleting debug stream\n");
   delete fDebugStream;
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
  Float_t mpt = track->GetSigned1Pt();
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
  if (!fDebugStream) fDebugStream = new TTreeSRedirector("TPCCalibTracksGain.root");

   for (Int_t iCluster = 0; iCluster < 159; iCluster++) {
      AliTPCclusterMI* cluster = seed->GetClusterPointer(iCluster);
      if (cluster) AddCluster(cluster);
   }
   DumpTrack(seed);
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

   if (!cluster) {
      Error("AddCluster", "Cluster not valid.");
      return;
   }
   
   Double_t xx[5];
   
   // this block is for using the center of the region with the same pad size as origin for the fit function,
   // comment it out if you want the origin at lx=ly=0.
   Float_t centerPad[3] = {0};
   AliTPCROC* tpcROC = AliTPCROC::Instance();
   Int_t padType = GetPadType(cluster->GetX());
   Int_t IOROC = (padType == 0) ? 0 : tpcROC->GetNInnerSector();
   //tpcROC->GetPositionLocal(IOROC, tpcROC->GetNRows(IOROC)/2, tpcROC->GetNPads(IOROC, tpcROC->GetNRows(IOROC)/2)/2, centerPad);  // use this instead of the switch statement if you want to calculate the center of the ROC and not the center of the regions with the same pad size
   switch (padType) {
      case kShortPads:
         tpcROC->GetPositionLocal(IOROC, tpcROC->GetNRows(IOROC)/2, tpcROC->GetNPads(IOROC, tpcROC->GetNRows(IOROC)/2)/2, centerPad);
         break;
      case kMediumPads:
         tpcROC->GetPositionLocal(IOROC, 64/2, tpcROC->GetNPads(IOROC, 64/2)/2, centerPad);
         break;
      case kLongPads:
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

   Int_t segment = cluster->GetDetector() % 36;
   Double_t q = ((Double_t)(cluster->GetMax()));  // note: no normalization to pad size!
   Double_t sqrtQ = TMath::Sqrt(q);
   Double_t logQ = fgkM * TMath::Log(1 + q / fgkM);
   if (padType == kShortPads) {
      ((TLinearFitter*)(fShortFitter->At(segment)))->AddPoint(xx, q);
      ((TLinearFitter*)(fSqrtShortFitter->At(segment)))->AddPoint(xx, sqrtQ);
      ((TLinearFitter*)(fLogShortFitter->At(segment)))->AddPoint(xx, logQ);
      fNShortClusters[segment]++;
   } else if (padType == kMediumPads) {
      ((TLinearFitter*)(fMediumFitter->At(segment)))->AddPoint(xx, q);
      ((TLinearFitter*)(fSqrtMediumFitter->At(segment)))->AddPoint(xx, sqrtQ);
      ((TLinearFitter*)(fLogMediumFitter->At(segment)))->AddPoint(xx, logQ);
      fNMediumClusters[segment]++;
   } else if (padType == kLongPads) {
      ((TLinearFitter*)(fLongFitter->At(segment)))->AddPoint(xx, q);
      ((TLinearFitter*)(fSqrtLongFitter->At(segment)))->AddPoint(xx, sqrtQ);
      ((TLinearFitter*)(fLogLongFitter->At(segment)))->AddPoint(xx, logQ);
      fNLongClusters[segment]++;
   }
}

Int_t AliTPCcalibTracksGain::Evaluate(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD* fitParam, TVectorD* fitError, Double_t* redChi2, Bool_t robust) {
   //
   // Evaluate the tracks for obtaining the calibration information.
   // segment specifies the segment (IROC & OROC) for which the fitter shall be evaluated (it can take values from
   // 0 to 35; segment == i means IROC# i and OROC# (i+36) are meant).
   // padType is one of kShortPads, kMediumPads, kLongPads.
   // fitType is one of kSimpleFitter, kSqrtFitter, kLogFitter.
   // fitParam should be a TVectorD object with 6 elements where the fit parameters will be written to,
   // the same is valid for fitError, which will contain the errors of course.
   // redChi2 is a pointer to an Int_t value which will contain the reduced chi^2 of the fit.
   // If fitParam, fitError or redChi2 is a null pointer, the corresponding value is not calculated.
   // robust specifies wether the fitter's robust fitting mode shall be used (use with caution, it takes looooong!)
   //

   TLinearFitter* fitter = GetFitter(segment, padType, fitType);
   Int_t NClusters = 0;
   switch (padType) {
      case kShortPads:
         NClusters = fNShortClusters[segment];
         break;
      case kMediumPads:
         NClusters = fNMediumClusters[segment];
         break;
      case kLongPads:
         NClusters = fNLongClusters[segment];
         break;
   }

   if (robust) fitter->EvalRobust();
   else fitter->Eval();
   
   if (redChi2) *redChi2 = fitter->GetChisquare()/(NClusters - 6);
   if (fitParam) fitter->GetParameters(*fitParam);
   if (fitError) {
      fitter->GetErrors(*fitError);
      *fitError *= (redChi2) ? (TMath::Sqrt(*redChi2)) : (TMath::Sqrt(fitter->GetChisquare()/(NClusters - 6)));
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
      case kShortPads:
         startRow = 0;
         endRow = ROCfitted->GetNrows();
         tpcROC->GetPositionLocal(sector, endRow/2, ROCfitted->GetNPads(endRow/2)/2, centerPad);
         break;
      case kMediumPads:
         startRow = 0;
         endRow = 64;
         tpcROC->GetPositionLocal(sector, endRow/2, ROCfitted->GetNPads(endRow/2)/2, centerPad);
         break;
      case kLongPads:
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
   if ((Int_t)(roc1->GetSector()) < fTPCparam->GetNInnerSector()) return 0;
   if ((Int_t)(roc2->GetSector()) < fTPCparam->GetNInnerSector()) return 0;
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

void AliTPCcalibTracksGain::GetParameters(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD &fitParam) {
   //
   // Puts the fit parameters for the specified segment (IROC & OROC), padType and fitType
   // into the fitParam TVectorD (which should contain 6 elements).
   // padType is one of kShortPads, kMediumPads, kLongPads. fitType is one of kSimpleFitter, kSqrtFitter, kLogFitter.
   // Note: The fitter has to be evaluated first!
   //

   GetFitter(segment, padType, fitType)->GetParameters(fitParam);
}

void AliTPCcalibTracksGain::GetErrors(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD &fitError) {
   //
   // Puts the fit parameter errors for the specified segment (IROC & OROC), padType and fitType
   // into the fitParam TVectorD (which should contain 6 elements).
   // padType is one of kShortPads, kMediumPads, kLongPads. fitType is one of kSimpleFitter, kSqrtFitter, kLogFitter.
   // Note: The fitter has to be evaluated first!
   //

   GetFitter(segment, padType, fitType)->GetErrors(fitError);
   fitError *= TMath::Sqrt(GetRedChi2(segment, padType, fitType));
}

Double_t AliTPCcalibTracksGain::GetRedChi2(UInt_t segment, UInt_t padType, UInt_t fitType) {
   //
   // Returns the reduced chi^2 value for the specified segment, padType and fitType.
   // padType is one of kShortPads, kMediumPads, kLongPads. fitType is one of kSimpleFitter, kSqrtFitter, kLogFitter.
   // Note: The fitter has to be evaluated first!
   //

   Int_t NClusters = 0;
   switch (padType) {
      case kShortPads:
         NClusters = fNShortClusters[segment];
         break;
      case kMediumPads:
         NClusters = fNMediumClusters[segment];
         break;
      case kLongPads:
         NClusters = fNLongClusters[segment];
         break;
   }
   return GetFitter(segment, padType, fitType)->GetChisquare()/(NClusters - 6);
}

void AliTPCcalibTracksGain::GetCovarianceMatrix(UInt_t segment, UInt_t padType, UInt_t fitType, TMatrixD& covMatrix) {
   //
   // Returns the covariance matrix for the specified segment, padType, fitType.
   // padType is one of kShortPads, kMediumPads, kLongPads. fitType is one of kSimpleFitter, kSqrtFitter, kLogFitter.
   //

   GetFitter(segment, padType, fitType)->GetCovarianceMatrix(covMatrix);
}

TLinearFitter* AliTPCcalibTracksGain::GetFitter(UInt_t segment, UInt_t padType, UInt_t fitType) {
   //
   // Returns the TLinearFitter object for the specified segment, padType, fitType.
   // padType is one of kShortPads, kMediumPads, kLongPads. fitType is one of kSimpleFitter, kSqrtFitter, kLogFitter.
   //
   
   if (segment >= 36 || padType > 2 || fitType > 2) {
      Error("GetFitter", "Index out of bounds.");
      return 0;
   }
   switch (padType) {
      case kShortPads:
         return (TLinearFitter*)((fitType == kSimpleFitter) ? fShortFitter->At(segment) : ((fitType == kSqrtFitter) ? fSqrtShortFitter->At(segment) : fLogShortFitter->At(segment)));
      case kMediumPads:
         return (TLinearFitter*)((fitType == kSimpleFitter) ? fMediumFitter->At(segment) : ((fitType == kSqrtFitter) ? fSqrtMediumFitter->At(segment) : fLogMediumFitter->At(segment)));
      case kLongPads:
         return (TLinearFitter*)((fitType == kSimpleFitter) ? fLongFitter->At(segment) : ((fitType == kSqrtFitter) ? fSqrtLongFitter->At(segment) : fLogLongFitter->At(segment)));
   }
   return 0;
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



void  AliTPCcalibTracksGain::DumpTrack(AliTPCseed * track){
  //
  //  Dump track information to the  stream
  //   
   
  
  (*fDebugStream)<<"Track"<<
    "Track.="<<track<<       // track information
    "\n";
  Int_t rows[200];
  for (Int_t ipad=0; ipad<3; ipad++){
    GetDedx(track,ipad, rows);
  }
  
}

Bool_t   AliTPCcalibTracksGain::GetDedx(AliTPCseed * track, Int_t padType, Int_t *rows){
  //
  // GetDedx for given sector for given track
  // padType - type of pads
  //
  Int_t firstRow=0, lastRow=0;
  Int_t minRow=100;
  Float_t xcenter=0;
  const Float_t ktany = TMath::Tan(TMath::DegToRad()*10);
  const Float_t kedgey =4.;
  if (padType==0){
    firstRow=0;
    lastRow= fTPCparam->GetNRowLow();
    xcenter= 108.47;
  }
  if (padType==1){
    firstRow= fTPCparam->GetNRowLow();
    lastRow= fTPCparam->GetNRowLow()+fTPCparam->GetNRowUp1();
    xcenter= 166.60;
  }
  if (padType==2){
    firstRow= fTPCparam->GetNRowLow()+fTPCparam->GetNRowUp1();
    lastRow= fTPCparam->GetNRowLow()+fTPCparam->GetNRowUp();
    xcenter =222.6;
  }
  minRow= (lastRow-firstRow)/2;
  //
  //
  Int_t nclusters=0;
  Int_t nclustersNE=0; // number of not edge clusters
  Int_t lastSector=-1;
  Float_t amplitudeQ[100];
  Float_t amplitudeM[100];
  Int_t   rowIn[100];
  Int_t   index[100];
  //
  static TLinearFitter fitY(2,"pol1");
  static TLinearFitter fitZ(2,"pol1");
  static TVectorD      parY(2);
  static TVectorD      parZ(2);
  fitY.ClearPoints();
  fitZ.ClearPoints();
  TVectorD meanPos(6);
  

  for (Int_t iCluster = firstRow; iCluster < lastRow; iCluster++) {
    AliTPCclusterMI* cluster = track->GetClusterPointer(iCluster);
    if (cluster) {
      Int_t detector = cluster->GetDetector() ;
      if (lastSector==-1) lastSector= detector;
      if (lastSector!=detector) continue;
      amplitudeQ[nclusters]=cluster->GetQ();
      amplitudeM[nclusters]=cluster->GetMax();
      rowIn[nclusters]=iCluster;
      nclusters++;
      Double_t dx=cluster->GetX()-xcenter;
      Double_t y=cluster->GetY();
      Double_t z=cluster->GetZ();      
      fitY.AddPoint(&dx, y);
      fitZ.AddPoint(&dx, z);
      meanPos[0]+=dx;
      meanPos[1]+=dx;
      meanPos[2]+=y;
      meanPos[3]+=y*y;
      meanPos[4]+=z;
      meanPos[5]+=z*z;
      if (TMath::Abs(cluster->GetY())<cluster->GetX()*ktany-kedgey) nclustersNE++;
    }
  }

  if (nclusters<minRow/2) return kFALSE;
  if (nclustersNE<minRow/2) return kFALSE;
  for (Int_t i=0;i<6;i++) meanPos[i]/=Double_t(nclusters);
  fitY.Eval();
  fitZ.Eval();
  fitY.GetParameters(parY);
  fitZ.GetParameters(parZ);
  //
  // calculate truncated mean
  //
  TMath::Sort(nclusters,amplitudeQ,index, kFALSE);

  TVectorD dedxQ(5);
  TVectorD dedxM(5);
  Float_t ndedx[5];
  for (Int_t i=0; i<5; i++){
    dedxQ[i]=0;
    dedxM[i]=0;					
    ndedx[i]=0;
  }
  //
  // dedx calculation
  //
  Int_t inonEdge=0;
  for (Int_t i=0; i<nclusters; i++){
    Int_t rowSorted = rowIn[index[i]]; 
    AliTPCclusterMI* cluster = track->GetClusterPointer(rowSorted);
    
    if (TMath::Abs(cluster->GetY())> cluster->GetX()*ktany-kedgey) continue;  //don't take edge clusters
    inonEdge++;
    if (inonEdge<nclustersNE*0.5) { 
      ndedx[0]++; 
      dedxQ[0]+=amplitudeQ[index[i]];
      dedxM[0]+=amplitudeM[index[i]];
    }
    if (inonEdge<nclustersNE*0.6) { 
      ndedx[1]++; 
      dedxQ[1]+=amplitudeQ[index[i]];
      dedxM[1]+=amplitudeM[index[i]];
    }
    if (inonEdge<nclustersNE*0.7) { 
      ndedx[2]++; 
      dedxQ[2]+=amplitudeQ[index[i]];
      dedxM[2]+=amplitudeM[index[i]];
    }
    if (inonEdge<nclustersNE*0.8) { 
      ndedx[3]++; 
      dedxQ[3]+=amplitudeQ[index[i]];
      dedxM[3]+=amplitudeM[index[i]];
    }
    if (inonEdge<nclustersNE*0.9) { 
      ndedx[4]++; 
      dedxQ[4]+=amplitudeQ[index[i]];
      dedxM[4]+=amplitudeM[index[i]];
    }
  }
  for (Int_t i=0; i<5; i++){
    dedxQ[i]/=ndedx[i];
    dedxM[i]/=ndedx[i];
  }

  inonEdge=0;
  for (Int_t i=0; i<nclusters; i++){
    Int_t rowSorted = rowIn[index[i]]; 
    AliTPCclusterMI* cluster = track->GetClusterPointer(rowSorted);
    if (!cluster) {
      printf("Problem\n");
      continue;
    }
    if (TMath::Abs(cluster->GetY())< cluster->GetX()*ktany-kedgey) inonEdge++;
    Float_t dedge    = cluster->GetX()*ktany-TMath::Abs(cluster->GetY());
    Float_t fraction = Float_t(i)/Float_t(nclusters);
    Float_t fraction2= Float_t(inonEdge)/Float_t(nclustersNE);
    Float_t momenta = track->GetP();
    Float_t mdedx   = track->GetdEdx();
    (*fDebugStream)<<"dEdx"<<
      "Cl.="<<cluster<<    //cluster of interest
      "P="<<momenta<<     // track momenta
      "dedx="<<mdedx<<    // mean dedx - corrected for angle
      "IPad="<<padType<<   // pad type 0..2
      "xc="<<xcenter<<     // x center of chamber
      "dedxQ.="<<&dedxQ<<  // dedxQ  - total charge
      "dedxM.="<<&dedxM<<  // dedxM  - maximal charge
      "fraction="<<fraction<<  // fraction - order in statistic (0,1)
      "fraction2="<<fraction2<<  // fraction - order in statistic (0,1)
      "dedge="<<dedge<<      // distance to thhe edge
      "parY.="<<&parY<<      // line fit
      "parZ.="<<&parZ<<      // line fit
      "meanPos.="<<&meanPos<< // mean position (dx, dx^2, y,y^2, z, z^2)
      "\n";

  }


}

