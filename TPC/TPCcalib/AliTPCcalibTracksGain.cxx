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

////////////////////////////////////////////////////////////////////////////
//                                                                        
//       === Calibration class for gain calibration using tracks ===
//
//    1) Genereal idea
//    ================
//       A 6-parametric parabolic function
//
//          G(x, y) = p0 + p1*x + p2*y + p3*x^2 + p4*y^2 + p5 * x*y
//
//       is fitted to the maximum charge values or total charge values of
//       all the clusters contained in the tracks that are added to this
//       object. This fit is performed for each read out chamber, in fact even
//       for each type of pad sizes (thus for one segment, which consists of
//       an IROC and an OROC, there are three fitters used, corresponding to
//       the three pad sizes). The coordinate origin is at the center of the
//       particular pad size region on each ROC.
//
//       Because of the Landau nature of the charge deposition we use
//       different "types" of fitters instead of one to minimize the effect
//       of the long Landau tail. The difference between the fitters is only
//       the charge value, that is put into them, i.e. the charge is subject
//       to a transformation. At this point we use three different fit types:
//
//          a) simple: the charge is put in as it is
//          b) sqrt:   the square root of the charge is put into the fitter
//          c) log:    fgkM * Log(1+q/fgkM) is put into the fitter, with
//                     q being the untransformed charge and fgkM=25
//
//       The results of the fits may be visualized and further used by
//       creating an AliTPCCalROC or AliTPCCalPad. You may specify to undo
//       the transformation and/or to normalize to the pad size.
//
//       Not every track you add to this object is actually used for
//       calibration. There are some cuts and conditions to exclude bad
//       tracks, e.g. a pt cut to cut out tracks with too much charge
//       deposition or a cut on edge clusters which are not fully
//       registered and don't give a usable signal.
//
//    2) Interface / usage
//    ====================
//       For each track to be added you need to call Process().
//       This method expects an AliTPCseed, which contains the necessary
//       cluster information. At the moment of writing this information
//       is stored in an AliESDfriend corresponding to an AliESD.
//       You may also call AddTrack() if you don't want the cuts and
//       other quality conditions to kick in (thus forcing the object to
//       accept the track) or AddCluster() for adding single clusters.
//       Call one of the Evaluate functions to evaluate the fitter(s) and
//       to retrieve the fit parameters, erros and so on. You can also
//       do this later on by using the different Getters.
//
//       The visualization methods CreateFitCalPad() and CreateFitCalROC()
//       are straight forward to use.
//
//       Note: If you plan to write this object to a ROOT file, make sure
//             you evaluate all the fitters *before* writing, because due
//             to a bug in the fitter component writing fitters doesn't
//             work properly (yet). Be aware that you cannot re-evaluate
//             the fitters after loading this object from file.
//             (This will be gone for a new ROOT version > v5-17-05)
//                                                                        
////////////////////////////////////////////////////////////////////////////

#include <TPDGCode.h>
#include <TStyle.h>
#include "TSystem.h"
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
#include <TCollection.h>
#include <TIterator.h>

//
// AliRoot includes
//
#include "AliMagF.h"
#include "AliMathBase.h"
//
#include "AliTPCROC.h"
#include "AliTPCParamSR.h"
#include "AliTPCCalROC.h"
#include "AliTPCCalPad.h"
//
#include "AliTracker.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h" 
#include "AliTPCseed.h"
#include "AliTPCclusterMI.h"
#include "AliTPCcalibTracksCuts.h"
#include "AliTPCFitPad.h"

// REMOVE ALL OF THIS
#include <TTree.h>
#include "AliESDEvent.h"

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

ClassImp(AliTPCcalibTracksGain)

const Bool_t   AliTPCcalibTracksGain::fgkUseTotalCharge = kTRUE;
const Double_t AliTPCcalibTracksGain::fgkM = 25.;
const char*    AliTPCcalibTracksGain::fgkDebugStreamFileName = "TPCCalibTracksGain.root";
AliTPCParamSR* AliTPCcalibTracksGain::fgTPCparam = new AliTPCParamSR();

AliTPCcalibTracksGain::AliTPCcalibTracksGain() :
   TNamed(),
   fDebugCalPadRaw(0),
   fDebugCalPadCorr(0),
   fDebugStream(0),
   fSimpleFitter(0),
   fSqrtFitter(0),
   fLogFitter(0),
   fSingleSectorFitter(0),
   fPrevIter(0),
   fDebugStreamPrefix(0),
   fCuts(0)
{
   //
   // Default constructor.
   //
}

AliTPCcalibTracksGain::AliTPCcalibTracksGain(const AliTPCcalibTracksGain& obj) :
   TNamed(obj)
{
   //
   // Copy constructor.
   //

   fDebugCalPadRaw = new AliTPCCalPad(*(obj.fDebugCalPadRaw));
   fDebugCalPadCorr = new AliTPCCalPad(*(obj.fDebugCalPadCorr));
   fSimpleFitter = new AliTPCFitPad(*(obj.fSimpleFitter));
   fSqrtFitter = new AliTPCFitPad(*(obj.fSqrtFitter));
   fLogFitter = new AliTPCFitPad(*(obj.fLogFitter));
   fSingleSectorFitter = new AliTPCFitPad(*(obj.fSingleSectorFitter));
   fPrevIter = new AliTPCcalibTracksGain(*(obj.fPrevIter));
   fDebugStreamPrefix = new TObjString(*(obj.fDebugStreamPrefix));
   fCuts = new AliTPCcalibTracksCuts(*(obj.fCuts));
}

AliTPCcalibTracksGain& AliTPCcalibTracksGain::operator=(const AliTPCcalibTracksGain& rhs) {
   //
   // Assignment operator.
   //

   if (this != &rhs) {
      TNamed::operator=(rhs);
      fDebugCalPadRaw = new AliTPCCalPad(*(rhs.fDebugCalPadRaw));
      fDebugCalPadCorr = new AliTPCCalPad(*(rhs.fDebugCalPadCorr));
      fSimpleFitter = new AliTPCFitPad(*(rhs.fSimpleFitter));
      fSqrtFitter = new AliTPCFitPad(*(rhs.fSqrtFitter));
      fLogFitter = new AliTPCFitPad(*(rhs.fLogFitter));
      fSingleSectorFitter = new AliTPCFitPad(*(rhs.fSingleSectorFitter));
      fPrevIter = new AliTPCcalibTracksGain(*(rhs.fPrevIter));
      fDebugStreamPrefix = new TObjString(*(rhs.fDebugStreamPrefix));
      fCuts = new AliTPCcalibTracksCuts(*(rhs.fCuts));
   }
   return *this;
}

AliTPCcalibTracksGain::AliTPCcalibTracksGain(const char* name, const char* title, AliTPCcalibTracksCuts* cuts, TNamed* debugStreamPrefix, AliTPCcalibTracksGain* prevIter) :
   TNamed(name, title),
   fDebugCalPadRaw(0),
   fDebugCalPadCorr(0),
   fDebugStream(0),
   fSimpleFitter(0),
   fSqrtFitter(0),
   fLogFitter(0),
   fSingleSectorFitter(0),
   fPrevIter(0),
   fDebugStreamPrefix(0),
   fCuts(0)
{
   //
   // Constructor.
   //
   
   //TH1::AddDirectory(kFALSE);
   G__SetCatchException(0);

   fCuts = cuts;
   if (debugStreamPrefix) fDebugStreamPrefix = new TObjString(debugStreamPrefix->GetTitle());
   fPrevIter = prevIter;

   fSimpleFitter = new AliTPCFitPad(7, "hyp6", "");
   fSqrtFitter   = new AliTPCFitPad(7, "hyp6", "");
   fLogFitter    = new AliTPCFitPad(7, "hyp6", "");
   fSingleSectorFitter = new AliTPCFitPad(7, "hyp6", ""); // just for debugging


   // just for debugging
   fTotalTracks     = 0;
   fAcceptedTracks  = 0;
   fDebugCalPadRaw  = new AliTPCCalPad("DebugCalPadRaw", "All clusters simply added up before correction");
   fDebugCalPadCorr = new AliTPCCalPad("DebugCalPadCorr", "All clusters simply added up after correction");
   
   // this will be gone for the a new ROOT version > v5-17-05
   for (UInt_t i = 0; i < 36; i++) {
      fNShortClusters[i]  = 0;
      fNMediumClusters[i] = 0;
      fNLongClusters[i]   = 0;
   }
 }

AliTPCcalibTracksGain::~AliTPCcalibTracksGain() {
   //
   // Destructor.
   //
  Info("Destructor","");
   if (fSimpleFitter) delete fSimpleFitter;
   if (fSqrtFitter) delete fSqrtFitter;
   if (fLogFitter) delete fLogFitter;
   if (fSingleSectorFitter) delete fSingleSectorFitter;


   if (fDebugStream) {
      //fDebugStream->GetFile()->Close();
      printf("Deleting debug stream object\n");
      delete fDebugStream;
   }

   if (fDebugStreamPrefix) delete fDebugStreamPrefix;

   if (fDebugCalPadRaw) delete fDebugCalPadRaw;
   if (fDebugCalPadCorr) delete fDebugCalPadCorr;
}

void AliTPCcalibTracksGain::Terminate(){
   //
   // Close Debug streamer
   //
   Evaluate();
   if (fDebugStream) {
     delete fDebugStream;
     fDebugStream = 0;
   }

   if (fDebugStreamPrefix) {
      TString debugStreamPrefix = fDebugStreamPrefix->GetString();
      TString destFile("");
      destFile += debugStreamPrefix;
      destFile += "/";
      destFile += gSystem->HostName();
      destFile += "_TPCCalibTracksGain.root";
      if (debugStreamPrefix.BeginsWith("root://")) {
         TFile::Cp("TPCCalibTracksGain.root", destFile.Data());
      } else {
         TString command("mv TPCCalibTracksGain.root ");
         command += destFile;
         gSystem->Exec(command.Data());
      }
      //char *prefix = "/d/alice11/miranov/simulHEAD0907/pp/calib/";
      //char command[4000];
      //sprintf(command,"mv TPCCalibTracksGain.root %s/%s_TPCCalibTracksGain.root", prefix, gSystem->HostName());
      //gSystem->Exec(command);
   }
}

void AliTPCcalibTracksGain::AddInfo(TChain* chain, char* debugStreamPrefix, char* prevIterFileName) {
   // 
   // Add some parameters from a previous run (AliTPCcalibTracksGain object contained
   // in root file fileName) to the chain.
   // Note: The parameters are *not* added to this class, you need to do it later by retrieving
   // the parameters from the chain and passing them to the constructor!
   //

   if (debugStreamPrefix) {
      TNamed* objDebugStreamPrefix = new TNamed("debugStreamPrefix", debugStreamPrefix);
      chain->GetUserInfo()->AddLast((TObject*)objDebugStreamPrefix);
   }
   
   if (prevIterFileName) {
      TFile paramFile(prevIterFileName);
      if (paramFile.IsZombie()) {
         printf("File %s not found. Continuing without z dependence parametrisation.\n", prevIterFileName);
         return;
      }
      
      AliTPCcalibTracksGain *prevIter = (AliTPCcalibTracksGain*)paramFile.Get("calibTracksGain");
      if (prevIter) {
         //TVectorD* param = new TVectorD(2);
         //chain->GetUserInfo()->AddLast((TObject*)param);
         chain->GetUserInfo()->AddLast((TObject*)prevIter);
      } else
         printf("No calibTracksGain object found. Continuing without z dependence parametrisation.\n");
   }
}

Int_t AliTPCcalibTracksGain::AcceptTrack(AliTPCseed* track) {
   //
   // Decides whether to accept a track or not.
   // Tracks are discarded, if due to edge effects, the number of clusters
   // is too low, the ratio of the number of clusters and the findable clusters is too low
   // or the transverse momentum is too low.
   // The corresponding cut values are specified in the fCuts member.
   //
   
  //  if (track->GetNumberOfClusters() < fCuts->GetMinClusters()) return 1;
//    if ((TMath::Abs(track->GetY() / track->GetX()) > fCuts->GetEdgeYXCutNoise())
//       && (TMath::Abs(track->GetTgl()) < fCuts->GetEdgeThetaCutNoise())) return 2;
//    if (track->GetNumberOfClusters() / (track->GetNFoundable()+1.) < fCuts->GetMinRatio()) return 3;
//    if (TMath::Abs(track->GetSigned1Pt()) > fCuts->GetMax1pt()) return 4;
   
   //if (track->GetPt() < 50.) return kFALSE;
   return 0;
}

/*Bool_t AliTPCcalibTracksGain::AcceptCluster(AliTPCclusterMI* cluster) {
   //
   // Decides whether to accept a cluster or not.
   //

   // cluster type < 0: edge cluster
   // cluster type = 0: "golden" (i.e. good) cluster
   // cluster type > 0: overlapping cluster
   if (cluster->GetType() != 0) { Info("AcceptCluster", "Cluster not accepted (type %d)", cluster->GetType()); return kFALSE;}

   // remove edge clusters
   
   
   return kTRUE;
}*/

void AliTPCcalibTracksGain::Process(AliTPCseed* seed) {
   //
   // Main method to be called when a new seed is supposed to be processed
   // and be used for gain calibration. Its quality is checked before it
   // is added.
   //
   
   fTotalTracks++;
   //fTrackPt->Fill(seed->GetSignedPt());
   Int_t status = AcceptTrack(seed);
   if (status != 0) { /*cout << "Track not accepted (reason " << status << ")" << endl;*/ return; }
   fAcceptedTracks++;
   AddTrack(seed);
}

Long64_t AliTPCcalibTracksGain::Merge(TCollection *list) {
   //
   // Merge() merges the results of all AliTPCcalibTracksGain objects contained in
   // list, thus allowing a distributed computation of several files, e.g. on PROOF.
   // The merged results are stored in the data members of the AliTPCcalibTracksGain
   // object used for calling the Merge method.
   // The return value is 0 /*the total number of tracks used for calibration*/ if the merge
   // is successful, otherwise it is -1.
   //

   if (!list || list->IsEmpty()) return -1;
   
   // reset the data members first
   //if (fSimpleFitter) delete fSimpleFitter;
   //if (fSqrtFitter)   delete fSqrtFitter;
   //if (fLogFitter)    delete fLogFitter;
   //if (fSingleSectorFitter) delete fSingleSectorFitter;  // just for debugging
   if (!fSimpleFitter) {
     fSimpleFitter = new AliTPCFitPad(7, "hyp6", "");
     fSqrtFitter   = new AliTPCFitPad(7, "hyp6", "");
     fLogFitter    = new AliTPCFitPad(7, "hyp6", "");
     fSingleSectorFitter = new AliTPCFitPad(7, "hyp6", "");  // just for debugging
   }

   // this will be gone for the a new ROOT version > v5-17-05
   for (UInt_t i = 0; i < 36; i++) {
      fNShortClusters[i]  = 0;
      fNMediumClusters[i] = 0;
      fNLongClusters[i]   = 0;
   }

   // just for debugging
   if (fDebugCalPadRaw)  delete fDebugCalPadRaw;
   if (fDebugCalPadCorr) delete fDebugCalPadCorr;
   fDebugCalPadRaw  = new AliTPCCalPad("DebugCalPadRaw", "All clusters simply added up before correction");
   fDebugCalPadCorr = new AliTPCCalPad("DebugCalPadCorr", "All clusters simply added up after correction");
   fTotalTracks     = 0;
   fAcceptedTracks  = 0;
   
   TIterator* iter = list->MakeIterator();
   AliTPCcalibTracksGain* cal = 0;
   
   while ((cal = (AliTPCcalibTracksGain*)iter->Next())) {
      if (!cal->InheritsFrom(AliTPCcalibTracksGain::Class())) {
         Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
         return -1;
      }
      Add(cal);
   }
   return 0;
}

void AliTPCcalibTracksGain::Add(AliTPCcalibTracksGain* cal) {
   //
   // Adds another AliTPCcalibTracksGain object to this object.
   //
   
   fSimpleFitter->Add(cal->fSimpleFitter);
   fSqrtFitter->Add(cal->fSqrtFitter);
   fLogFitter->Add(cal->fLogFitter);
   fSingleSectorFitter->Add(cal->fSingleSectorFitter);  // just for debugging

   // this will be gone for the a new ROOT version > v5-17-05
   for (UInt_t iSegment = 0; iSegment < 36; iSegment++) {
      fNShortClusters[iSegment] += cal->fNShortClusters[iSegment];
      fNMediumClusters[iSegment] += cal->fNMediumClusters[iSegment];
      fNLongClusters[iSegment] += cal->fNLongClusters[iSegment];
   }
   
   // just for debugging, remove me
   fTotalTracks += cal->fTotalTracks;
   fAcceptedTracks += cal->fAcceptedTracks;
   fDebugCalPadRaw->Add(cal->fDebugCalPadRaw);
   fDebugCalPadCorr->Add(cal->fDebugCalPadCorr);

   // Let's see later what to do with fCuts and fDebugStream and fDebugStreamPrefix
}

void AliTPCcalibTracksGain::AddTrack(AliTPCseed* seed) {
   //
   // The clusters making up the track (seed) are added to various fit functions.
   // See AddCluster(...) for more detail.
   //
   
   if (!fDebugStream) fDebugStream = new TTreeSRedirector(fgkDebugStreamFileName);
   DumpTrack(seed);

   /*AliTPCcalibTracksGain::PreProcess preProc(seed);
   for (Int_t iCluster = 0; iCluster < 159; iCluster++) {
      AliTPCclusterMI* cluster = seed->GetClusterPointer(iCluster);
      if (cluster && preProc.IsClusterAccepted(iCluster)) AddCluster(cluster, preProc);
   }*/
}
   
void AliTPCcalibTracksGain::AddCluster2(AliTPCclusterMI* cluster, Float_t momenta, Float_t mdedx, Int_t padType,
            Float_t xcenter, TVectorD dedxQ, TVectorD dedxM, Float_t fraction, Float_t fraction2, Float_t dedge,
            TVectorD parY, TVectorD parZ, TVectorD meanPos) {
   if (!cluster) {
      Error("AddCluster", "Cluster not valid.");
      return;
   }

   if (dedge < 3.) return;
   if (fraction2 > 0.7) return;

   //Int_t padType = GetPadType(cluster->GetX());
   Double_t xx[6];
   //Double_t centerPad[2] = {0};
   //AliTPCFitPad::GetPadRegionCenterLocal(padType, centerPad);
   //xx[0] = cluster->GetX() - centerPad[0];
   //xx[1] = cluster->GetY() - centerPad[1];
   xx[0] = cluster->GetX() - xcenter;
   xx[1] = cluster->GetY();
   xx[2] = xx[0] * xx[0];
   xx[3] = xx[1] * xx[1];
   xx[4] = xx[0] * xx[1];
   xx[5] = TMath::Abs(cluster->GetZ()) - TMath::Abs(meanPos[4]);

   Int_t segment = cluster->GetDetector() % 36;
   Double_t q = fgkUseTotalCharge ? ((Double_t)(cluster->GetQ())) : ((Double_t)(cluster->GetMax()));  // note: no normalization to pad size!

   // just for debugging
   Int_t row = 0;
   Int_t pad = 0;
   GetRowPad(cluster->GetX(), cluster->GetY(), row, pad);
   fDebugCalPadRaw->GetCalROC(cluster->GetDetector())->SetValue(row, pad, q + fDebugCalPadRaw->GetCalROC(cluster->GetDetector())->GetValue(row, pad));
   
   // correct charge by normalising to mean charge per track
   q /= dedxQ[2];

   // correct charge for dependency on z distance using previous iteration
   
   // just for debugging
   fDebugCalPadCorr->GetCalROC(cluster->GetDetector())->SetValue(row, pad, q + fDebugCalPadCorr->GetCalROC(cluster->GetDetector())->GetValue(row, pad));

   Double_t sqrtQ = TMath::Sqrt(q);
   Double_t logQ = fgkM * TMath::Log(1 + q / fgkM);
   fSimpleFitter->GetFitter(segment, padType)->AddPoint(xx, q);
   fSqrtFitter->GetFitter(segment, padType)->AddPoint(xx, sqrtQ);
   fLogFitter->GetFitter(segment, padType)->AddPoint(xx, logQ);
   fSingleSectorFitter->GetFitter(0, padType)->AddPoint(xx, q);  // just for debugging

   Double_t zz = TMath::Abs(cluster->GetZ()) - TMath::Abs(meanPos[4]);
   
   // this will be gone for the a new ROOT version > v5-17-05
   if (padType == kShortPads)
      fNShortClusters[segment]++;
   else if (padType == kMediumPads)
      fNMediumClusters[segment]++;
   else if (padType == kLongPads)
      fNLongClusters[segment]++;
}

void AliTPCcalibTracksGain::AddCluster(AliTPCclusterMI* cluster, AliTPCcalibTracksGain::PreProcess& preProc) {
   //
   // Adds cluster to the appropriate fitter for later analysis.
   // The charge used for the fit is the maximum charge for this specific cluster or the
   // accumulated charge per cluster, depending on the value of fgkUseTotalCharge.
   // Depending on the pad size where the cluster is registered, the value will be put in
   // the appropriate fitter. Furthermore, for each pad size three different types of fitters
   // are used. The fit functions are the same for all fitters (parabolic functions), but the value
   // added to each fitter is different. The simple fitter gets the charge plugged in as is, the sqrt fitter
   // gets the square root of the charge, and the log fitter gets fgkM*(1+q/fgkM), where q is the original charge
   // and fgkM==25.
   // The preProc object is used for passing information which was gained by preprocessing the seed to which the
   // cluster belongs. It is used for correcting the charge due to various effects, e.g. its dependence
   // on the track length over each pad.
   //

   if (!cluster) {
      Error("AddCluster", "Cluster not valid.");
      return;
   }
   
   Double_t xx[5];
   
   Int_t padType = GetPadType(cluster->GetX());
   Double_t centerPad[2] = {0};
   // this line is for using the center of the region with the same pad size as origin for the fit function,
   // comment out the appropriate lines in AliTPCCalPadRegion::GetPadRegionCenterLocal() and here if you
   // want the origin at lx=ly=0.
   AliTPCFitPad::GetPadRegionCenterLocal(padType, centerPad);
   
   xx[0] = cluster->GetX() - centerPad[0];
   xx[1] = cluster->GetY() - centerPad[1];
   
   //xx[0] = cluster->GetX();   // if you want the origin of the fit func at lx=ly=0
   //xx[1] = cluster->GetY();   // if you want the origin of the fit func at lx=ly=0
   xx[2] = xx[0] * xx[0];
   xx[3] = xx[1] * xx[1];
   xx[4] = xx[0] * xx[1];

   Int_t segment = cluster->GetDetector() % 36;
   Double_t q = fgkUseTotalCharge ? ((Double_t)(cluster->GetQ())) : ((Double_t)(cluster->GetMax()));  // note: no normalization to pad size!
   //cerr << "AngleTrackPadrow(" << segment << ", " << padType << ") == " << preProc.GetAngleTrackPadrow(segment, padType) * TMath::RadToDeg() << endl;
   //cerr << "AngleTrackBeam  (" << segment << ", " << padType << ") == " << preProc.GetAngleTrackBeam(segment, padType) * TMath:: RadToDeg() << endl;
   //cerr << "Correction factor == " << TMath::Abs(TMath::Sin(preProc.GetAngleTrackPadrow(segment, padType))*TMath::Sin(preProc.GetAngleTrackBeam(segment, padType)))/GetPadLength(cluster->GetX()) << endl;

   // just for debugging
   Int_t row = 0;
   Int_t pad = 0;
   GetRowPad(cluster->GetX(), cluster->GetY(), row, pad);
   fDebugCalPadRaw->GetCalROC(cluster->GetDetector())->SetValue(row, pad, q + fDebugCalPadRaw->GetCalROC(cluster->GetDetector())->GetValue(row, pad));
   
   // correct charge by normalising to mean charge per track
   //q = TMath::Abs(TMath::Sin(preProc.GetAngleTrackPadrow(segment, padType))*TMath::Sin(preProc.GetAngleTrackBeam(segment, padType)))*q/(GetPadLength(cluster->GetX()) /* * preProc.GetMeanCharge(segment, padType) */);
   //q *= TMath::Abs(TMath::Sin(preProc.GetAngleTrackPadrow(segment, padType))) / (/*GetPadLength(cluster->GetX()) * */preProc.GetMeanCharge(segment, padType));
   q /= preProc.GetMeanCharge(segment, padType);
   //q /= (1 + 0.36 * TMath::Abs(1/TMath::Tan(preProc.GetAngleTrackPadrow(segment, padType))));

   // correct charge for dependency on z distance using previous iteration
   
   // just for debugging
   fDebugCalPadCorr->GetCalROC(cluster->GetDetector())->SetValue(row, pad, q + fDebugCalPadCorr->GetCalROC(cluster->GetDetector())->GetValue(row, pad));

   Double_t sqrtQ = TMath::Sqrt(q);
   Double_t logQ = fgkM * TMath::Log(1 + q / fgkM);
   fSimpleFitter->GetFitter(segment, padType)->AddPoint(xx, q);
   fSqrtFitter->GetFitter(segment, padType)->AddPoint(xx, sqrtQ);
   fLogFitter->GetFitter(segment, padType)->AddPoint(xx, logQ);

   Double_t zz = TMath::Abs(cluster->GetZ()) - TMath::Abs(preProc.GetMeanZ(segment, padType));
   
   // this will be gone for the a new ROOT version > v5-17-05
   if (padType == kShortPads)
      fNShortClusters[segment]++;
   else if (padType == kMediumPads)
      fNMediumClusters[segment]++;
   else if (padType == kLongPads)
      fNLongClusters[segment]++;
}

void AliTPCcalibTracksGain::Evaluate(Bool_t robust, Double_t frac) {
   //
   // Evaluates all fitters contained in this object.
   // If the robust option is set to kTRUE a robust fit is performed with frac as
   // the minimal fraction of good points (see TLinearFitter::EvalRobust for details).
   // Beware: Robust fitting is much slower!
   //
   
   fSimpleFitter->Evaluate(robust, frac);
   fSqrtFitter->Evaluate(robust, frac);
   fLogFitter->Evaluate(robust, frac);
   fSingleSectorFitter->Evaluate(robust, frac);  // just for debugging
}

Int_t AliTPCcalibTracksGain::Evaluate(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD* fitParam, TVectorD* fitError, Double_t* redChi2, Bool_t robust, Double_t frac) {
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
   // If the robust option is set to kTRUE a robust fit is performed with frac as
   // the minimal fraction of good points (see TLinearFitter::EvalRobust for details).
   // Beware: Robust fitting is much slower!
   // Evaluate() returns the number of clusters in the specified padType.
   //

   TLinearFitter* fitter = GetFitter(segment, padType, fitType);
   Int_t NClusters = 0;
   // this will be gone for a new ROOT version > v5-17-05
   // and replaced by
   // UInt_t NClusters = fSimpleFitter->GetFitter(segment, padType)->GetNpoints();
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

   if (robust) fitter->EvalRobust(frac);
   else fitter->Eval();
   
   if (redChi2) *redChi2 = fitter->GetChisquare()/(NClusters - 6);
   if (fitParam) fitter->GetParameters(*fitParam);
   if (fitError) {
      fitter->GetErrors(*fitError);
      *fitError *= (redChi2) ? (TMath::Sqrt(*redChi2)) : (TMath::Sqrt(fitter->GetChisquare()/(NClusters - 6)));
   }

   return NClusters;
}

AliTPCCalPad* AliTPCcalibTracksGain::CreateFitCalPad(UInt_t fitType, Bool_t undoTransformation, Bool_t normalizeToPadSize) {
   TObjArray tpc(72);
   for (UInt_t iSector = 0; iSector < 72; iSector++)
      tpc.Add(CreateFitCalROC(iSector, fitType, undoTransformation, normalizeToPadSize));
   return new AliTPCCalPad(&tpc);
}

AliTPCCalROC* AliTPCcalibTracksGain::CreateFitCalROC(UInt_t sector, UInt_t fitType, Bool_t undoTransformation, Bool_t normalizeToPadSize) {
   TVectorD par(7);
   if (sector < 36) {
      GetParameters(sector % 36, 0, fitType, par);
      return CreateFitCalROC(sector, 0, par, fitType, undoTransformation, normalizeToPadSize);
   }
   else {
      GetParameters(sector % 36, 1, fitType, par);
      AliTPCCalROC* roc1 = CreateFitCalROC(sector, 1, par, fitType, undoTransformation, normalizeToPadSize);
      GetParameters(sector % 36, 2, fitType, par);
      AliTPCCalROC* roc2 = CreateFitCalROC(sector, 2, par, fitType, undoTransformation, normalizeToPadSize);
      AliTPCCalROC* roc3 = CreateCombinedCalROC(roc1, roc2);
      delete roc1;
      delete roc2;
      return roc3;
   }
}

AliTPCCalROC* AliTPCcalibTracksGain::CreateFitCalROC(UInt_t sector, UInt_t padType, TVectorD &fitParam, UInt_t fitType, Bool_t undoTransformation, Bool_t normalizeToPadSize) {
   //
   // This function is essentially a copy of AliTPCCalROC::CreateGlobalFitCalROC(...), with the
   // modifications, that the center of the region of same pad size is used as the origin
   // of the fit function instead of the center of the ROC.
   // The possibility of a linear fit is removed as well because it is not needed.
   // Only values for pads with the given pad size are calculated, the rest is 0.
   // Set undoTransformation for undoing the transformation that was applied to the
   // charge values before they were put into the fitter (thus allowing comparison to the original
   // charge values). For fitType use 0 for the simple fitter, 1 for the sqrt fitter, 2 for the log fitter.
   // If normalizeToPadSize is true, the values are normalized to the pad size.
   // Please be aware, that you even need to specify the fitType if you want to normalize to the pad size without
   // undoing the transformation (because normalizing involves undoing the trafo first, then normalizing, then
   // applying the trafo again).
   // Please note: The normalization to the pad size is a simple linear scaling with the pad length, which
   //              actually doesn't describe reality!
   //
   
   Float_t dlx, dly;
   Double_t centerPad[2] = {0};
   Float_t localXY[3] = {0};
   AliTPCROC* tpcROC = AliTPCROC::Instance();
   if ((padType == 0 && sector >= tpcROC->GetNInnerSector()) || (padType > 0 && sector < tpcROC->GetNInnerSector()) || sector >= tpcROC->GetNSector())
      return 0;
   AliTPCCalROC* ROCfitted = new AliTPCCalROC(sector);
   //tpcROC->GetPositionLocal(sector, ROCfitted->GetNrows()/2, ROCfitted->GetNPads(ROCfitted->GetNrows()/2)/2, centerPad);  // use this instead of the switch statement if you want to calculate the center of the ROC and not the center of the regions with the same pad size
   UInt_t startRow = 0;
   UInt_t endRow = 0;
   switch (padType) {
      case kShortPads:
         startRow = 0;
         endRow = ROCfitted->GetNrows();
         break;
      case kMediumPads:
         startRow = 0;
         endRow = 64;
         break;
      case kLongPads:
         startRow = 64;
         endRow = ROCfitted->GetNrows();
         break;
   }

   AliTPCFitPad::GetPadRegionCenterLocal(padType, centerPad);   
   Double_t value = 0;
   for (UInt_t irow = startRow; irow < endRow; irow++) {
      for (UInt_t ipad = 0; ipad < ROCfitted->GetNPads(irow); ipad++) {
         tpcROC->GetPositionLocal(sector, irow, ipad, localXY);   // calculate position localXY by pad and row number
         dlx = localXY[0] - centerPad[0];
         dly = localXY[1] - centerPad[1];
         value = fitParam[0] + fitParam[1]*dlx + fitParam[2]*dly + fitParam[3]*dlx*dlx + fitParam[4]*dly*dly + fitParam[5]*dlx*dly;

         // Let q' = value be the transformed value without any pad size corrections,
         // let T be the transformation and let l be the pad size
         //    1) don't undo transformation, don't normalize: return q'
         //    2) undo transformation,       don't normalize: return T^{-1} q'
         //    3) undo transformation,       normalize:       return (T^{-1} q') / l
         //    4) don't undo transformation, normalize:       return T((T^{-1} q') / l)
         if (!undoTransformation && !normalizeToPadSize) {/* value remains unchanged */}  // (1)
         else {                                                                           // (2), (3), (4)
            //calculate T^{-1}
            switch (fitType) {
               case  0: /* value remains unchanged */ break;
               case  1: value = value * value; break;
               case  2: value = (TMath::Exp(value / fgkM) - 1) * fgkM; break;
               default: Error("CreateFitCalROC", "Wrong fit type."); break;
            }
            if (normalizeToPadSize) value /= GetPadLength(localXY[0]);                    // (3)
         }
         if (!undoTransformation && normalizeToPadSize) {                                 // (4)
            // calculate T
            switch (fitType) {
               case  0: /* value remains unchanged */ break;
               case  1: value = TMath::Sqrt(value); break;
               case  2: value = fgkM * TMath::Log(1 + value / fgkM); break;
               default: Error("CreateFitCalROC", "Wrong fit type."); break;
            }
         }
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
   if ((Int_t)(roc1->GetSector()) < fgTPCparam->GetNInnerSector()) return 0;
   if ((Int_t)(roc2->GetSector()) < fgTPCparam->GetNInnerSector()) return 0;
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
   
   switch (fitType) {
      case kSimpleFitter:
         return fSimpleFitter->GetFitter(segment, padType);
      case kSqrtFitter:
         return fSqrtFitter->GetFitter(segment, padType);
      case kLogFitter:
         return fLogFitter->GetFitter(segment, padType);
      case 3:  // just for debugging
         return fSingleSectorFitter->GetFitter(0, padType);
   }
   return 0;
}

Double_t AliTPCcalibTracksGain::GetPadLength(Double_t lx) {
   //
   // The function returns 0.75 for an IROC, 1. for an OROC at medium pad size position,
   // 1.5 for an OROC at long pad size position, -1 if out of bounds.
   //

   Double_t irocLow = fgTPCparam->GetPadRowRadiiLow(0) - fgTPCparam->GetInnerPadPitchLength()/2;
   Double_t irocUp = fgTPCparam->GetPadRowRadiiLow(fgTPCparam->GetNRowLow()-1) + fgTPCparam->GetInnerPadPitchLength()/2;
   Double_t orocLow1 = fgTPCparam->GetPadRowRadiiUp(0) - fgTPCparam->GetOuter1PadPitchLength()/2;
   Double_t orocUp1 = fgTPCparam->GetPadRowRadiiUp(fgTPCparam->GetNRowUp1()-1) + fgTPCparam->GetOuter1PadPitchLength()/2;
   Double_t orocLow2 = fgTPCparam->GetPadRowRadiiUp(fgTPCparam->GetNRowUp1()) - fgTPCparam->GetOuter2PadPitchLength()/2;
   Double_t orocUp2 = fgTPCparam->GetPadRowRadiiUp(fgTPCparam->GetNRowUp()-1) + fgTPCparam->GetOuter2PadPitchLength()/2;
   
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

// ONLY FOR DEBUGGING PURPOSES - REMOVE ME WHEN NOT NEEDED ANYMORE
Bool_t AliTPCcalibTracksGain::GetRowPad(Double_t lx, Double_t ly, Int_t& row, Int_t& pad) {
   //
   // Calculate the row and pad number when the local coordinates are given.
   // Returns kFALSE if the position is out of range, otherwise return kTRUE.
   // WARNING: This function is preliminary and probably isn't very accurate!!
   //
   
   Double_t irocLow = fgTPCparam->GetPadRowRadiiLow(0) - fgTPCparam->GetInnerPadPitchLength()/2;
   //Double_t irocUp = fgTPCparam->GetPadRowRadiiLow(fgTPCparam->GetNRowLow()-1) + fgTPCparam->GetInnerPadPitchLength()/2;
   Double_t orocLow1 = fgTPCparam->GetPadRowRadiiUp(0) - fgTPCparam->GetOuter1PadPitchLength()/2;
   //Double_t orocUp1 = fgTPCparam->GetPadRowRadiiUp(fgTPCparam->GetNRowUp1()-1) + fgTPCparam->GetOuter1PadPitchLength()/2;
   Double_t orocLow2 = fgTPCparam->GetPadRowRadiiUp(fgTPCparam->GetNRowUp1()) - fgTPCparam->GetOuter2PadPitchLength()/2;
   //Double_t orocUp2 = fgTPCparam->GetPadRowRadiiUp(fgTPCparam->GetNRowUp()-1) + fgTPCparam->GetOuter2PadPitchLength()/2;

   if (GetPadType(lx) == 0) {
      row = (Int_t)((lx - irocLow) / fgTPCparam->GetInnerPadPitchLength());
      pad = (Int_t)((ly + fgTPCparam->GetYInner(row)) / fgTPCparam->GetInnerPadPitchWidth());
   } else if (GetPadType(lx) == 1) {
      row = (Int_t)((lx - orocLow1) / fgTPCparam->GetOuter1PadPitchLength());
      pad = (Int_t)((ly + fgTPCparam->GetYOuter(row)) / fgTPCparam->GetOuterPadPitchWidth());
   } else if (GetPadType(lx) == 2) {
      row = fgTPCparam->GetNRowUp1() + (Int_t)((lx - orocLow2) / fgTPCparam->GetOuter2PadPitchLength());
      pad = (Int_t)((ly + fgTPCparam->GetYOuter(row)) / fgTPCparam->GetOuterPadPitchWidth());
   }
   else return kFALSE;
   return kTRUE;
}

void AliTPCcalibTracksGain::DumpTrack(AliTPCseed* track) {
   //
   //  Dump track information to the debug stream
   //
   
   (*fDebugStream) << "Track" <<
      "Track.=" << track <<        // track information
      "\n";
   Int_t rows[200];
   for (Int_t ipad = 0; ipad < 3; ipad++) {
      GetDedx(track, ipad, rows);
   }
}

Bool_t AliTPCcalibTracksGain::GetDedx(AliTPCseed* track, Int_t padType, Int_t* rows) {
   //
   // GetDedx for given sector for given track
   // padType - type of pads
   //
   
   Int_t firstRow = 0, lastRow = 0;
   Int_t minRow = 100;
   Float_t xcenter = 0;
   const Float_t ktany = TMath::Tan(TMath::DegToRad() * 10);
   const Float_t kedgey = 4.;
   if (padType == 0) {
      firstRow = 0;
      lastRow = fgTPCparam->GetNRowLow();
      xcenter = 108.47;
   }
   if (padType == 1) {
      firstRow = fgTPCparam->GetNRowLow();
      lastRow = fgTPCparam->GetNRowLow() + fgTPCparam->GetNRowUp1();
      xcenter = 166.60;
   }
   if (padType == 2) {
      firstRow = fgTPCparam->GetNRowLow() + fgTPCparam->GetNRowUp1();
      lastRow = fgTPCparam->GetNRowLow() + fgTPCparam->GetNRowUp();
      xcenter = 222.6;
   }
   minRow = (lastRow - firstRow) / 2;
   //
   //
   Int_t nclusters = 0;
   Int_t nclustersNE = 0; // number of not edge clusters
   Int_t lastSector = -1;
   Float_t amplitudeQ[100];
   Float_t amplitudeM[100];
   Int_t rowIn[100];
   Int_t index[100];
   //
   static TLinearFitter fitY(2, "pol1");
   static TLinearFitter fitZ(2, "pol1");
   static TVectorD parY(2);
   static TVectorD parZ(2);
   fitY.ClearPoints();
   fitZ.ClearPoints();
   TVectorD meanPos(6);
   
   for (Int_t iCluster = firstRow; iCluster < lastRow; iCluster++) {
      AliTPCclusterMI* cluster = track->GetClusterPointer(iCluster);
      if (cluster) {
         Int_t detector = cluster->GetDetector() ;
         if (lastSector == -1) lastSector = detector;
         if (lastSector != detector) continue;
         amplitudeQ[nclusters] = cluster->GetQ();
         amplitudeM[nclusters] = cluster->GetMax();
         rowIn[nclusters] = iCluster;
         nclusters++;
         Double_t dx = cluster->GetX() - xcenter;
         Double_t y = cluster->GetY();
         Double_t z = cluster->GetZ();
         fitY.AddPoint(&dx, y);
         fitZ.AddPoint(&dx, z);
         meanPos[0] += dx;
         meanPos[1] += dx;
         meanPos[2] += y;
         meanPos[3] += y*y;
         meanPos[4] += z;
         meanPos[5] += z*z;
         if (TMath::Abs(cluster->GetY()) < cluster->GetX()*ktany - kedgey) nclustersNE++;
      }
   }
   
   if (nclusters < minRow / 2) return kFALSE;
   if (nclustersNE < minRow / 2) return kFALSE;
   for (Int_t i = 0; i < 6; i++) meanPos[i] /= Double_t(nclusters);
   fitY.Eval();
   fitZ.Eval();
   fitY.GetParameters(parY);
   fitZ.GetParameters(parZ);
   //
   // calculate truncated mean
   //
   TMath::Sort(nclusters, amplitudeQ, index, kFALSE);
   
   TVectorD dedxQ(5);
   TVectorD dedxM(5);
   Float_t ndedx[5];
   for (Int_t i = 0; i < 5; i++) {
      dedxQ[i] = 0;
      dedxM[i] = 0;
      ndedx[i] = 0;
   }
   //
   // dedx calculation
   //
   Int_t inonEdge = 0;
   for (Int_t i = 0; i < nclusters; i++) {
      Int_t rowSorted = rowIn[index[i]];
      AliTPCclusterMI* cluster = track->GetClusterPointer(rowSorted);
      
      if (TMath::Abs(cluster->GetY()) > cluster->GetX()*ktany - kedgey) continue;  //don't take edge clusters
      inonEdge++;
      if (inonEdge < nclustersNE * 0.5) {
         ndedx[0]++;
         dedxQ[0] += amplitudeQ[index[i]];
         dedxM[0] += amplitudeM[index[i]];
      }
      if (inonEdge < nclustersNE * 0.6) {
         ndedx[1]++;
         dedxQ[1] += amplitudeQ[index[i]];
         dedxM[1] += amplitudeM[index[i]];
      }
      if (inonEdge < nclustersNE * 0.7) {
         ndedx[2]++;
         dedxQ[2] += amplitudeQ[index[i]];
         dedxM[2] += amplitudeM[index[i]];
      }
      if (inonEdge < nclustersNE * 0.8) {
         ndedx[3]++;
         dedxQ[3] += amplitudeQ[index[i]];
         dedxM[3] += amplitudeM[index[i]];
      }
      if (inonEdge < nclustersNE * 0.9) {
         ndedx[4]++;
         dedxQ[4] += amplitudeQ[index[i]];
         dedxM[4] += amplitudeM[index[i]];
      }
   }
   for (Int_t i = 0; i < 5; i++) {
      dedxQ[i] /= ndedx[i];
      dedxM[i] /= ndedx[i];
   }
   
   inonEdge = 0;
   for (Int_t i = 0; i < nclusters; i++) {
      Int_t rowSorted = rowIn[index[i]];
      AliTPCclusterMI* cluster = track->GetClusterPointer(rowSorted);
      if (!cluster) {
         printf("Problem\n");
         continue;
      }
      if (TMath::Abs(cluster->GetY()) < cluster->GetX()*ktany - kedgey) inonEdge++;
      Float_t dedge = cluster->GetX()*ktany - TMath::Abs(cluster->GetY());
      Float_t fraction = Float_t(i) / Float_t(nclusters);
      Float_t fraction2 = Float_t(inonEdge) / Float_t(nclustersNE);
      Float_t momenta = track->GetP();
      Float_t mdedx = track->GetdEdx();

      AddCluster2(cluster, momenta, mdedx, padType, xcenter, dedxQ, dedxM, fraction, fraction2, dedge, parY, parZ, meanPos);

      (*fDebugStream) << "dEdx" <<
         "Cl.=" << cluster <<           // cluster of interest
         "P=" << momenta <<             // track momenta
         "dedx=" << mdedx <<            // mean dedx - corrected for angle
         "IPad=" << padType <<          // pad type 0..2
         "xc=" << xcenter <<            // x center of chamber
         "dedxQ.=" << &dedxQ <<         // dedxQ  - total charge
         "dedxM.=" << &dedxM <<         // dedxM  - maximal charge
         "fraction=" << fraction <<     // fraction - order in statistic (0,1)
         "fraction2=" << fraction2 <<   // fraction - order in statistic (0,1)
         "dedge=" << dedge <<           // distance to the edge
         "parY.=" << &parY <<           // line fit
         "parZ.=" << &parZ <<           // line fit
         "meanPos.=" << &meanPos <<     // mean position (dx, dx^2, y,y^2, z, z^2)
         "\n";
   }
   return kTRUE;
}

AliTPCcalibTracksGain::PreProcess::PreProcess(AliTPCseed* seed) :
   fSeed(seed)
{
   //
   // Constructor. Preprocesses the track contained in seed. After that
   // all relevant values gained from preprocessing can be accessed using
   // the getter methods.
   // For each pad region a line fit using the clusters of this region is made
   // for obtaining the angle between track and padrow.
   // The mean charge for the clusters of each pad region is calculated.
   // For all these calculations inappropriate clusters are not used
   // and are accordingly marked (which can be found out with IsClusterAccepted).
   //

   // initialize data members and constants
   fAngleTrackPadrow = new Double_t[AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes()];
   fAngleTrackBeam   = new Double_t[AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes()];
   fMeanCharge       = new Double_t[AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes()];
   fMeanZ            = new Double_t[AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes()];
   for (UInt_t i = 0; i < AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes(); i++) {
      fMeanCharge[i] = 0;
      fMeanZ[i] = 0;
   }
   for (UInt_t i = 0; i < 159; i++) fAcceptedClusters[i] = kFALSE;
   const Double_t kSmallNumber = 1E-20;
   const Double_t kTan10 = TMath::Tan(10*TMath::DegToRad());
   const Double_t kCos10 = TMath::Cos(10*TMath::DegToRad());
      
   Int_t nClusters[AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes()];                              // this will be gone for a new ROOT version > v5-17-05
   for (UInt_t i = 0; i < AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes(); i++) nClusters[i] = 0; // this will be gone for a new ROOT version > v5-17-05
   Int_t nAllClusters[AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes()];
   for (UInt_t i = 0; i < AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes(); i++) nAllClusters[i] = 0;
   
   Double_t clusterCharges[AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes()][100];
   Int_t    clusterIndices[AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes()][100];
   Int_t    clusterRows   [AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes()][100];
   for (Int_t iCluster = 0; iCluster < 159; iCluster++) {
      AliTPCclusterMI* cluster = fSeed->GetClusterPointer(iCluster);
      // check if cluster is good cluster
      if (cluster) {
         Int_t segment = cluster->GetDetector() % 36;
         Int_t padType = (iCluster < AliTPCcalibTracksGain::fgTPCparam->GetNRowLow()) ? 0 : ((iCluster < fgTPCparam->GetNRowLow() + fgTPCparam->GetNRowUp1()) ? 1 : 2);
         Int_t index = segment + AliTPCFitPad::GetNSegments() * padType;
         nAllClusters[index]++;

         // cluster type < 0: edge cluster
         // cluster type = 0: "golden" (i.e. good) cluster
         // cluster type > 0: overlapping cluster
         // if (cluster->GetType() != 0) continue;
         
         // add up z positions for later calculation of mean z position
         fMeanZ[index] += cluster->GetZ();
         
         // remove edge clusters (only keep those farther away than 3 cm from the edge)
         Double_t edgeDistance = (cluster->GetX() * kTan10 - TMath::Abs(cluster->GetY())); //* kCos10; // the cos is for the shortest distance
         if (edgeDistance < 3.) continue;
         
         // put cluster charges and their row position in arrays for each pad region and count them
         clusterCharges[index][nClusters[index]] = AliTPCcalibTracksGain::fgkUseTotalCharge ? cluster->GetQ() : cluster->GetMax();
         clusterRows[index][nClusters[index]] = iCluster;
         nClusters[index]++;
      }
   }

   // order clusters according to their charge (in each pad region)
   //for (UInt_t index = 0; index < AliTPCFitPad::GetNSegments() * AliTPCFitPad::GetNPadTypes(); index++)
   //   TMath::Sort(nClusters[index], clusterCharges[index], clusterIndices[index], kFALSE);

   AliTPCFitPad localXYfitters(2, "pol1", ""); // used for determining the angle between XY projection of track and Y axis (padrow) - local c.s.
   AliTPCFitPad localXZfitters(2, "pol1", ""); // used for determining the angle between XZ projection of track and Z axis (beam) - local c.s.
   
   Int_t minRows[3] = {0};
   minRows[0] = AliTPCcalibTracksGain::fgTPCparam->GetNRowLow();
   minRows[1] = AliTPCcalibTracksGain::fgTPCparam->GetNRowUp1();
   minRows[2] = AliTPCcalibTracksGain::fgTPCparam->GetNRowUp() - AliTPCcalibTracksGain::fgTPCparam->GetNRowUp1();
   // the minimum number of rows occupied by clusters in a pad region should be a quarter of the available rows
   for (Int_t i = 0; i < 3; minRows[i++] /= 4);
   
   for (UInt_t iSegment = 0; iSegment < AliTPCFitPad::GetNSegments(); iSegment++) {
      for (UInt_t iPadType = 0; iPadType < AliTPCFitPad::GetNPadTypes(); iPadType++) {
         Int_t index = iSegment + AliTPCFitPad::GetNSegments() * iPadType;
         
         // remove pad regions with too few rows occupied by clusters
         if (nAllClusters[index] < minRows[iPadType]) continue;
         // calculate mean z position
         fMeanZ[index] /= nAllClusters[index];
         // order clusters according to their charge (in each pad region)
         TMath::Sort(nClusters[index], clusterCharges[index], clusterIndices[index], kFALSE);
         
         for (Int_t i = 0; i < nClusters[index]; i++) {
            Int_t iCluster = clusterRows[index][clusterIndices[index][i]];
            AliTPCclusterMI* cluster = fSeed->GetClusterPointer(iCluster);
            // keep only the lower 70% of the clusters ordered according to their charge
            Double_t fraction = (Double_t)i / (Double_t)(nClusters[index]);
            if (fraction > 0.7) break;
                        
            fAcceptedClusters[iCluster] = kTRUE;
            // put cluster data into fitters (I should move this before the fraction cut, because the charge values don't matter - the problem is the number of points -> additional counter)
            Double_t x = cluster->GetX();             // maybe it's better for fitter stability to use difference from pad region center instead
            localXYfitters.GetFitter(iSegment, iPadType)->AddPoint(&x, cluster->GetY());
            localXZfitters.GetFitter(iSegment, iPadType)->AddPoint(&x, cluster->GetZ());
            // add up charges for later calculation of mean charges
            fMeanCharge[index] += AliTPCcalibTracksGain::fgkUseTotalCharge ? cluster->GetQ() : cluster->GetMax();
         }
         // calculate mean charges
         if (nClusters[index] != 0) fMeanCharge[index] /= (Double_t)(nClusters[index]);
      }
   }
      
   // evaluate fitters and set corresponding angles or default values, respectively
   for (UInt_t iSegment = 0; iSegment < AliTPCFitPad::GetNSegments(); iSegment++) {
      for (UInt_t iPadType = 0; iPadType < AliTPCFitPad::GetNPadTypes(); iPadType++) {
         UInt_t index = iSegment + AliTPCFitPad::GetNSegments() * iPadType;
         // evaluate XY fitters
         TLinearFitter* fitter = localXYfitters.GetFitterSimple(iSegment, iPadType);
         if (fitter && nClusters[index] >= 2 && fitter->Eval() == 0) {                                    // this will be gone for a new ROOT version > v5-17-05 and replaced by (fitter && fitter->GetNPoints() >= 2 && fitter->Eval() == 0)
            Double_t slope = fitter->GetParameter(1);
            if (TMath::Abs(slope) < kSmallNumber)
               fAngleTrackPadrow[index] = TMath::Sign(TMath::Pi()/2, slope); // default value for small slopes
            else
               fAngleTrackPadrow[index] = TMath::ATan(1/slope);
         }
         else  fAngleTrackPadrow[index] = TMath::Pi()/2; // default value if not enough data points in fitter
         // evaluate XZ fitters
         fitter = localXZfitters.GetFitterSimple(iSegment, iPadType);
         if (fitter && nClusters[index] >= 2 && fitter->Eval() == 0) {                                   // this will be gone for a new ROOT version > v5-17-05 and replaced by (fitter && fitter->GetNPoints() >= 2 && fitter->Eval() == 0)
            Double_t slope = fitter->GetParameter(1);
            if (TMath::Abs(slope) < kSmallNumber)
               fAngleTrackBeam[index] = TMath::Sign(TMath::Pi()/2, slope); // default value for small slopes
            else
               fAngleTrackBeam[index] = TMath::ATan(1/slope);
         }
         else  fAngleTrackBeam[index] = TMath::Pi()/2; // default value if not enough data points in fitter, should be pi/2 for cosmics, something like pi/4 for pp
      }
   }
}

AliTPCcalibTracksGain::PreProcess::~PreProcess() {
   //
   // Destructor.
   //
   
   delete[] fAngleTrackPadrow;
   delete[] fAngleTrackBeam;
   delete[] fMeanCharge;
   delete[] fMeanZ;
}

// REMOVE ME, just for debugging
void AliTPCcalibTracksGain::testSeed(char* file, Int_t entry, Int_t tr) {
   // read a single AliTPCseed from AliESDs.root and AliESDfriends.root (with new AliRoot)
   TFile f(file);
   TTree* tree = (TTree*)(f.Get("esdTree"));
   tree->SetBranchStatus("*",1);
   AliESDEvent* fESDevent = new AliESDEvent();
   fESDevent->ReadFromTree(tree);
   tree->GetEntry(entry);
   AliESDfriend* fESDfriend = (AliESDfriend*)fESDevent->FindListObject("AliESDfriend");
   tree->SetBranchAddress("ESDfriend.",&fESDfriend);
   fESDevent->SetESDfriend(fESDfriend);
   AliESDtrack *esdTrack = (AliESDtrack*) fESDevent->GetTrack(tr);
   AliESDfriendTrack *friendtrack = (AliESDfriendTrack*) esdTrack->GetFriendTrack();
   AliTPCseed *seed = 0; TObject *cobject = 0;
   for (Int_t i = 0; ; i++){ cobject = friendtrack->GetCalibObject(i); if (!cobject) break; seed = dynamic_cast<AliTPCseed*>(cobject); if (seed) break;}

   AliTPCcalibTracksGain::PreProcess preProc(seed);
   for (Int_t iCluster = 0; iCluster < 159; iCluster++) {
      AliTPCclusterMI* cluster = seed->GetClusterPointer(iCluster);
      if (cluster) {
         Int_t padType = AliTPCcalibTracksGain::GetPadType(cluster->GetX());
         Int_t segment = cluster->GetDetector() % 36;
         cout << "AngleTrackPadrow(" << segment << ", " << padType << ") == " << preProc.GetAngleTrackPadrow(segment, padType) * TMath::RadToDeg() << endl;
         cout << "AngleTrackBeam  (" << segment << ", " << padType << ") == " << preProc.GetAngleTrackBeam(segment, padType) * TMath:: RadToDeg() << endl;
         cout << "Correction factor == " << TMath::Abs(TMath::Sin(preProc.GetAngleTrackPadrow(segment, padType))*TMath::Sin(preProc.GetAngleTrackBeam(segment, padType)))/AliTPCcalibTracksGain::GetPadLength(cluster->GetX()) << endl;
      }
   }
}
