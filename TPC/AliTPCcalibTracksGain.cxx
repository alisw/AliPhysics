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
//
//  Gain calibration using tracks
//
//  The main goal:
//    1.) Inner TPC gain alignement - (parabolic) parameterization (inside of one sector) 
//    2.) Angular and z-position correction (parabolic) parameterization 
//    3.) Sector gain alignment  
//   
//  Following histograms are accumulated
//    a.) Simple 1D histograms  per chamber  
//    b.) Profile histograms per chamber - local x dependence    
//    c.) 2D Profile histograms  - local x - fi dependence
//
//  To get the gain map - the simple solution - use the histograms - is not enough
//  The resulting mean amplitude map depends strongly on the track topology
//  These dependence can be reduced, taking into account angular effect, and diffusion effect 
//  Using proper fit modeles
//
//     
//
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
//
//    In order to debug some numerical algorithm all data data which are used for
//    fitters can be stored in the debug streamers. In case of fitting roblems the 
//    errors and tendencies can be checked.
//
//    Debug Streamers:
//      
// 
//
//
//
////////////////////////////////////////////////////////////////////////////

/*
  .x ~/UliStyle.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTAT");
  gSystem->Load("libTPCcalib");

  TFile fcalib("CalibObjects.root");
  TObjArray * array = (TObjArray*)fcalib.Get("TPCCalib");
  AliTPCcalibTracksGain * gain = ( AliTPCcalibTracksGain *)array->FindObject("calibTracksGain");

  //
  // Angular and drift correction
  //
  AliTPCClusterParam *param = new AliTPCClusterParam;param->SetInstance(param);
  gain->UpdateClusterParam(param);
  TF2 fdrty("fdrty","AliTPCClusterParam::SQnorm(0,0,x,y,0)",0,1,0,1)

  //
  // Make visual Tree - compare with Kr calibration
  // 
  AliTPCCalPad * gain010 = gain->CreateFitCalPad(0,kTRUE,0); gain010->SetName("CGain010");
  AliTPCCalPad * gain110 = gain->CreateFitCalPad(1,kTRUE,0); gain110->SetName("CGain110");
  AliTPCCalPad * gain210 = gain->CreateFitCalPad(2,kTRUE,0); gain210->SetName("CGain210");
  TFile fkr("/u/miranov/GainMap.root");
  AliTPCCalPad *gainKr = fkr.Get("GainMap"); fkr->SetName("KrGain");
  //
  AliTPCPreprocessorOnline * preprocesor = new AliTPCPreprocessorOnline;
  preprocesor->AddComponent(gain010);
  preprocesor->AddComponent(gain110);
  preprocesor->AddComponent(gain210);
  preprocesor->AddComponent(gainKr);
  preprocesor->DumpToFile("cosmicGain.root");
  //
  //
  //
  // Simple session using the debug streamers
  //

  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+")
  AliXRDPROOFtoolkit tool;   

  TChain * chain0 = tool.MakeChain("gain.txt","dEdx",0,1000000);
  TChain * chain1 = tool.MakeChain("gain.txt","Track",0,1000000);
  TChain * chain2 = tool.MakeChain("gain.txt","TrackG",0,1000000);
  chain0->Lookup();
  chain1->Lookup();
  chain2->Lookup();

  chain2->SetAlias("k1","1/0.855");
  chain2->SetAlias("k0","1/0.9928");
  chain2->SetAlias("k2","1/1.152");
  
*/



#include "AliTPCcalibTracksGain.h"


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
#include <TH1.h>
#include <TH3F.h>
#include <TLinearFitter.h>
#include <TTreeStream.h>
#include <TFile.h>
#include <TCollection.h>
#include <TIterator.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProof.h>
#include <TStatToolkit.h>

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
#include "AliTPCClusterParam.h"
#include "AliTPCcalibDB.h"
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
#include "TStatToolkit.h"
#include "TString.h"
#include "TCut.h"

//
#include <TTree.h>
#include "AliESDEvent.h"

/*
  
TFile f("TPCCalibTracksGain.root")

gSystem->Load("libPWGPP.so")
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
  AliTPCcalibBase(),
  fCuts(0),            // cuts that are used for sieving the tracks used for calibration
  //
  // Fitters
  //
  fSimpleFitter(0),         // simple fitter for short pads
  fSqrtFitter(0),           // sqrt fitter for medium pads
  fLogFitter(0),            // log fitter for long pads
  
  fFitter0M(0),              // fitting of the atenuation, angular correction, and mean chamber gain
  fFitter1M(0),              // fitting of the atenuation, angular correction, and mean chamber gain
  fFitter2M(0),              // fitting of the atenuation, angular correction, and mean chamber gain
  fFitter0T(0),              // fitting of the atenuation, angular correction, and mean chamber gain
  fFitter1T(0),              // fitting of the atenuation, angular correction, and mean chamber gain
  fFitter2T(0),              // fitting of the atenuation, angular correction, and mean chamber gain
  //
  fDFitter0M(0),              // fitting of the atenuation, angular correction
  fDFitter1M(0),              // fitting of the atenuation, angular correction
  fDFitter2M(0),              // fitting of the atenuation, angular correction
  fDFitter0T(0),              // fitting of the atenuation, angular correction
  fDFitter1T(0),              // fitting of the atenuation, angular correction
  fDFitter2T(0),              // fitting of the atenuation, angular correction
  //
  fSingleSectorFitter(0),   // just for debugging
  //
  // Counters
  //
  fTotalTracks(0),         // just for debugging
  fAcceptedTracks(0)      // just for debugging

{
   //
   // Default constructor.
   //
}

AliTPCcalibTracksGain::AliTPCcalibTracksGain(const AliTPCcalibTracksGain& obj) :
  AliTPCcalibBase(obj),
  fCuts(obj.fCuts),            // cuts that are used for sieving the tracks used for calibration
  //
  // Fitters
  //
  fSimpleFitter(obj.fSimpleFitter),         // simple fitter for short pads
  fSqrtFitter(obj.fSqrtFitter),           // sqrt fitter for medium pads
  fLogFitter(obj.fLogFitter),            // log fitter for long pads
  fFitter0M(obj.fFitter0M),
  fFitter1M(obj.fFitter1M),
  fFitter2M(obj.fFitter2M),
  fFitter0T(obj.fFitter0T),
  fFitter1T(obj.fFitter1T),
  fFitter2T(obj.fFitter2T),
  //
  fDFitter0M(obj.fDFitter0M),
  fDFitter1M(obj.fDFitter1M),
  fDFitter2M(obj.fDFitter2M),
  fDFitter0T(obj.fDFitter0T),
  fDFitter1T(obj.fDFitter1T),
  fDFitter2T(obj.fDFitter2T),
  fSingleSectorFitter(obj.fSingleSectorFitter),   // just for debugging
  //
  // Counters
  //
  fTotalTracks(obj.fTotalTracks),         // just for debugging
  fAcceptedTracks(obj.fAcceptedTracks)      // just for debugging

{
   //
   // Copy constructor.
   //
}

AliTPCcalibTracksGain& AliTPCcalibTracksGain::operator=(const AliTPCcalibTracksGain& rhs) {
   //
   // Assignment operator.
   //

   if (this != &rhs) {
      TNamed::operator=(rhs);
      fSimpleFitter = new AliTPCFitPad(*(rhs.fSimpleFitter));
      fSqrtFitter = new AliTPCFitPad(*(rhs.fSqrtFitter));
      fLogFitter = new AliTPCFitPad(*(rhs.fLogFitter));
      fSingleSectorFitter = new AliTPCFitPad(*(rhs.fSingleSectorFitter));
      fCuts = new AliTPCcalibTracksCuts(*(rhs.fCuts));
   }
   return *this;
}

AliTPCcalibTracksGain::AliTPCcalibTracksGain(const char* name, const char* title, AliTPCcalibTracksCuts* cuts) :
  AliTPCcalibBase(),
  fCuts(0),            // cuts that are used for sieving the tracks used for calibration
  //
  // Fitters
  //
  fSimpleFitter(0),         // simple fitter for short pads
  fSqrtFitter(0),           // sqrt fitter for medium pads
  fLogFitter(0),            // log fitter for long pads
  fFitter0M(0),              // fitting of the atenuation, angular correction, and mean chamber gain
  fFitter1M(0),              // fitting of the atenuation, angular correction, and mean chamber gain
  fFitter2M(0),              // fitting of the atenuation, angular correction, and mean chamber gain
  fFitter0T(0),              // fitting of the atenuation, angular correction, and mean chamber gain
  fFitter1T(0),              // fitting of the atenuation, angular correction, and mean chamber gain
  fFitter2T(0),              // fitting of the atenuation, angular correction, and mean chamber gain
  //
  fDFitter0M(0),              // fitting of the atenuation, angular correction
  fDFitter1M(0),              // fitting of the atenuation, angular correction
  fDFitter2M(0),              // fitting of the atenuation, angular correction
  fDFitter0T(0),              // fitting of the atenuation, angular correction
  fDFitter1T(0),              // fitting of the atenuation, angular correction
  fDFitter2T(0),              // fitting of the atenuation, angular correction
  fSingleSectorFitter(0),   // just for debugging
  //
  // Counters
  //
  fTotalTracks(0),         // just for debugging
  fAcceptedTracks(0)      // just for debugging
  
{
   //
   // Constructor.
   //   
   this->SetNameTitle(name, title);
   fCuts = cuts;
   //
   // Fitter initialization
   //
   fSimpleFitter = new AliTPCFitPad(8, "hyp7", "");
   fSqrtFitter   = new AliTPCFitPad(8, "hyp7", "");
   fLogFitter    = new AliTPCFitPad(8, "hyp7", "");
   fSingleSectorFitter = new AliTPCFitPad(8, "hyp7", "");
   // 
   fFitter0M      = new TLinearFitter(45,"hyp44");
   fFitter1M      = new TLinearFitter(45,"hyp44");
   fFitter2M      = new TLinearFitter(45,"hyp44");
   fFitter0T      = new TLinearFitter(45,"hyp44");
   fFitter1T      = new TLinearFitter(45,"hyp44");
   fFitter2T      = new TLinearFitter(45,"hyp44");
   //
   fDFitter0M      = new TLinearFitter(10,"hyp9");
   fDFitter1M      = new TLinearFitter(10,"hyp9");
   fDFitter2M      = new TLinearFitter(10,"hyp9");
   fDFitter0T      = new TLinearFitter(10,"hyp9");
   fDFitter1T      = new TLinearFitter(10,"hyp9");
   fDFitter2T      = new TLinearFitter(10,"hyp9");
   //
   // 
   fFitter0M->StoreData(kFALSE);
   fFitter1M->StoreData(kFALSE);
   fFitter2M->StoreData(kFALSE);
   fFitter0T->StoreData(kFALSE);
   fFitter1T->StoreData(kFALSE);
   fFitter2T->StoreData(kFALSE);   
   //
   fDFitter0M->StoreData(kFALSE);
   fDFitter1M->StoreData(kFALSE);
   fDFitter2M->StoreData(kFALSE);
   fDFitter0T->StoreData(kFALSE);
   fDFitter1T->StoreData(kFALSE);
   fDFitter2T->StoreData(kFALSE);   
   //
   //
   // just for debugging -counters
   //
   fTotalTracks     = 0;
   fAcceptedTracks  = 0;
}

AliTPCcalibTracksGain::~AliTPCcalibTracksGain() {
   //
   // Destructor.
   //

   Info("Destructor",":");
   if (fSimpleFitter) delete fSimpleFitter;
   if (fSqrtFitter) delete fSqrtFitter;
   if (fLogFitter) delete fLogFitter;
   if (fSingleSectorFitter) delete fSingleSectorFitter;

}




void AliTPCcalibTracksGain::Terminate(){
   //
   // Evaluate fitters and close the debug stream.
   // Also move or copy the debug stream, if a debugStreamPrefix is provided.
   //

   Evaluate();
   AliTPCcalibBase::Terminate();
}



void AliTPCcalibTracksGain::Process(AliTPCseed* seed) {
   //
   // Main method to be called when a new seed is supposed to be processed
   // and be used for gain calibration. Its quality is checked before it
   // is added.
   //
   

   fTotalTracks++;
   if (!fCuts->AcceptTrack(seed)) return;
   //
   // reinint on proof
   //   if (gProof){
     static Bool_t doinit= kTRUE;
     if (doinit){
       fSimpleFitter = new AliTPCFitPad(8, "hyp7", "");
       fSqrtFitter   = new AliTPCFitPad(8, "hyp7", "");
       fLogFitter    = new AliTPCFitPad(8, "hyp7", "");
       fSingleSectorFitter = new AliTPCFitPad(8, "hyp7", "");
       // 
       fFitter0M      = new TLinearFitter(45,"hyp44");
       fFitter1M      = new TLinearFitter(45,"hyp44");
       fFitter2M      = new TLinearFitter(45,"hyp44");
       fFitter0T      = new TLinearFitter(45,"hyp44");
       fFitter1T      = new TLinearFitter(45,"hyp44");
       fFitter2T      = new TLinearFitter(45,"hyp44");
       //
       fDFitter0M      = new TLinearFitter(10,"hyp9");
       fDFitter1M      = new TLinearFitter(10,"hyp9");
       fDFitter2M      = new TLinearFitter(10,"hyp9");
       fDFitter0T      = new TLinearFitter(10,"hyp9");
       fDFitter1T      = new TLinearFitter(10,"hyp9");
       fDFitter2T      = new TLinearFitter(10,"hyp9");
       doinit=kFALSE;
     }
     //}

   fAcceptedTracks++;
   AddTrack(seed);
}

Long64_t AliTPCcalibTracksGain::Merge(TCollection *list) {
   //
   // Merge() merges the results of all AliTPCcalibTracksGain objects contained in
   // list, thus allowing a distributed computation of several files, e.g. on PROOF.
   // The merged results are merged with the data members of the AliTPCcalibTracksGain
   // object used for calling the Merge method.
   // The return value is 0 /*the total number of tracks used for calibration*/ if the merge
   // is successful, otherwise it is -1.
   //

   if (!list || list->IsEmpty()) return -1;
   
   if (!fSimpleFitter) fSimpleFitter = new AliTPCFitPad(8, "hyp7", "");
   if (!fSqrtFitter)   fSqrtFitter   = new AliTPCFitPad(8, "hyp7", "");
   if (!fLogFitter)    fLogFitter    = new AliTPCFitPad(8, "hyp7", "");
   if (!fSingleSectorFitter) fSingleSectorFitter = new AliTPCFitPad(8, "hyp7", "");


   
   TIterator* iter = list->MakeIterator();
   AliTPCcalibTracksGain* cal = 0;
   
   while ((cal = (AliTPCcalibTracksGain*)iter->Next())) {
      if (!cal->InheritsFrom(AliTPCcalibTracksGain::Class())) {
	//Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
         return -1;
      }
      
      Add(cal);
   }
   return 0;
}

Float_t  AliTPCcalibTracksGain::GetGain(AliTPCclusterMI *cluster){
  //
  // get gain
  //
  Float_t gainPad= 1;
  AliTPCCalPad * gainMap =  AliTPCcalibDB::Instance()->GetDedxGainFactor();
  if (gainMap) {
    AliTPCCalROC * roc = gainMap->GetCalROC(cluster->GetDetector());
    gainPad  = roc->GetValue(cluster->GetRow(), TMath::Nint(cluster->GetPad()));
  }
  return gainPad;
}

Float_t   AliTPCcalibTracksGain::GetMaxNorm(AliTPCclusterMI * cluster, Float_t ky, Float_t kz){
  //
  // Get normalized amplituded if the gain map provided
  //
  AliTPCClusterParam * parcl = AliTPCcalibDB::Instance()->GetClusterParam();
  Float_t maxNorm =
    parcl->QmaxCorrection(cluster->GetDetector(), cluster->GetRow(),cluster->GetPad(), 
			  cluster->GetTimeBin(),ky,kz,0.5,0.5,1.6);

  return GetGain(cluster)*maxNorm;
}


Float_t   AliTPCcalibTracksGain::GetQNorm(AliTPCclusterMI * cluster,  Float_t ky, Float_t kz){
  //
  // Get normalized amplituded if the gain map provided
  //
  AliTPCClusterParam * parcl = AliTPCcalibDB::Instance()->GetClusterParam();
  Float_t totNorm = parcl->QtotCorrection(cluster->GetDetector(), cluster->GetRow(),cluster->GetPad(), 			      cluster->GetTimeBin(),ky,kz,0.5,0.5,cluster->GetQ(),2.5,1.6);
  return GetGain(cluster)*totNorm;
}



void AliTPCcalibTracksGain::Add(AliTPCcalibTracksGain* cal) {
   //
   // Adds another AliTPCcalibTracksGain object to this object.
   //
   
  fSimpleFitter->Add(cal->fSimpleFitter);
  fSqrtFitter->Add(cal->fSqrtFitter);
  fLogFitter->Add(cal->fLogFitter);
  fSingleSectorFitter->Add(cal->fSingleSectorFitter);
  //
  //
  //
  if (cal->fFitter0M->GetNpoints()>0) fFitter0M->Add(cal->fFitter0M);
  if (cal->fFitter1M->GetNpoints()>0) fFitter1M->Add(cal->fFitter1M);
  if (cal->fFitter2M->GetNpoints()>0) fFitter2M->Add(cal->fFitter2M);
  if (cal->fFitter0T->GetNpoints()>0) fFitter0T->Add(cal->fFitter0T);
  if (cal->fFitter1T->GetNpoints()>0) fFitter1T->Add(cal->fFitter1T);
  if (cal->fFitter2T->GetNpoints()>0) fFitter2T->Add(cal->fFitter2T);
  //
  if (cal->fDFitter0M->GetNpoints()>0) fDFitter0M->Add(cal->fDFitter0M);
  if (cal->fDFitter1M->GetNpoints()>0) fDFitter1M->Add(cal->fDFitter1M);
  if (cal->fDFitter2M->GetNpoints()>0) fDFitter2M->Add(cal->fDFitter2M);
  if (cal->fDFitter0T->GetNpoints()>0) fDFitter0T->Add(cal->fDFitter0T);
  if (cal->fDFitter1T->GetNpoints()>0) fDFitter1T->Add(cal->fDFitter1T);
  if (cal->fDFitter2T->GetNpoints()>0) fDFitter2T->Add(cal->fDFitter2T);
   //
   
   // just for debugging, remove me
   fTotalTracks += cal->fTotalTracks;
   fAcceptedTracks += cal->fAcceptedTracks;

}

void AliTPCcalibTracksGain::AddTrack(AliTPCseed* seed) {
   //
   // The clusters making up the track (seed) are added to various fit functions.
   // See AddCluster(...) for more detail.
   //
   
   DumpTrack(seed);   
}
   



void AliTPCcalibTracksGain::AddCluster(AliTPCclusterMI* cluster, Float_t /*momenta*/, Float_t/* mdedx*/, Int_t padType,
				       Float_t xcenter, TVectorD& dedxQ, TVectorD& /*dedxM*/, Float_t /*fraction*/, Float_t fraction2, Float_t dedge,
				       TVectorD& parY, TVectorD& parZ, TVectorD& meanPos) {
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
   //
  Float_t kedge     = 3;
  Float_t kfraction = 0.7;
  Int_t   kinorm    = 2;


  // Where to put selection on threshold? 
  // Defined by the Q/dEdxT variable - see debug streamer:
  //
  // Debug stream variables:  (Where tu cut ?)
  // chain0->Draw("Cl.fQ/dedxQ.fElements[1]>>his(100,0,3)","fraction2<0.6&&dedge>3","",1000000);
  //         mean 1  sigma 0.25
  // chain0->Draw("Cl.fMax/dedxM.fElements[1]>>his(100,0,3)","fraction2<0.6&&dedge>3","",1000000)
  //         mean 1  sigma 0.25
  // chain0->Draw("Cl.fQ/dedxQ.fElements[2]>>his(100,0,3)","fraction2<0.7&&dedge>3","",1000000)
  //         mean 1 sigma 0.29
  // chain0->Draw("Cl.fMax/dedxM.fElements[2]>>his(100,0,3)","fraction2<0.7&&dedge>3","",1000000)
  //         mean 1 sigma 0.27
  // chain0->Draw("Cl.fQ/dedxQ.fElements[3]>>his(100,0,3)","fraction2<0.8&&dedge>3","",1000000)
  //         mean 1 sigma 0.32
  // 
  // chain0->Draw("Cl.fQ/dedxQ.fElements[4]>>his(100,0,3)","fraction2<0.9&&dedge>3","",1000000)
  //         mean 1 sigma 0.4

  // Fraction choosen 0.7

   if (!cluster) {
      Error("AddCluster", "Cluster not valid.");
      return;
   }

   if (dedge < kedge) return;
   if (fraction2 > kfraction) return;

   //Int_t padType = GetPadType(cluster->GetX());
   Double_t xx[7];
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
   xx[6] = xx[5] * xx[5];

   //
   // Update fitters
   //
   Int_t segment = cluster->GetDetector() % 36;
   Double_t q = fgkUseTotalCharge ? 
     ((Double_t)(cluster->GetQ()/GetQNorm(cluster,parY[1], parZ[1]))) : ((Double_t)(cluster->GetMax()/GetMaxNorm(cluster,parY[1], parZ[1])));  
      
   // correct charge by normalising to mean charge per track
   q /= dedxQ[kinorm];

   // just for debugging

   Double_t sqrtQ = TMath::Sqrt(q);
   Double_t logQ = fgkM * TMath::Log(1 + q / fgkM);
   TLinearFitter * fitter =0;
   //
   fitter = fSimpleFitter->GetFitter(segment, padType);
   fitter->AddPoint(xx, q);
   //
   fitter = fSqrtFitter->GetFitter(segment, padType);
   fitter->AddPoint(xx, sqrtQ);
   //
   fitter = fLogFitter->GetFitter(segment, padType);
   fitter->AddPoint(xx, logQ);
   //
   fitter=fSingleSectorFitter->GetFitter(0, padType);
   fitter->AddPoint(xx, q);

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
   fSingleSectorFitter->Evaluate(robust, frac);
   fFitter0M->Eval();
   fFitter1M->Eval();
   fFitter2M->Eval();
   fFitter0T->Eval();
   fFitter1T->Eval();
   fFitter2T->Eval();
   //
   fDFitter0M->Eval();
   fDFitter1M->Eval();
   fDFitter2M->Eval();
   fDFitter0T->Eval();
   fDFitter1T->Eval();
   fDFitter2T->Eval();
}

AliTPCCalPad* AliTPCcalibTracksGain::CreateFitCalPad(UInt_t fitType, Bool_t undoTransformation, Bool_t normalizeToPadSize) {
  //
  // Creates the calibration object AliTPCcalPad using fitted parameterization
  //
   TObjArray tpc(72);
   for (UInt_t iSector = 0; iSector < 72; iSector++)
      tpc.Add(CreateFitCalROC(iSector, fitType, undoTransformation, normalizeToPadSize));
   return new AliTPCCalPad(&tpc);
}

AliTPCCalROC* AliTPCcalibTracksGain::CreateFitCalROC(UInt_t sector, UInt_t fitType, Bool_t undoTransformation, Bool_t normalizeToPadSize) {
  //
  // Create the AliTPCCalROC with the values per pad
  // sector  - sector of interest 
  // fitType - 
  //

   TVectorD par(8);
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
   AliTPCCalROC* lROCfitted = new AliTPCCalROC(sector);
   //tpcROC->GetPositionLocal(sector, lROCfitted->GetNrows()/2, lROCfitted->GetNPads(lROCfitted->GetNrows()/2)/2, centerPad);  // use this instead of the switch statement if you want to calculate the center of the ROC and not the center of the regions with the same pad size
   UInt_t startRow = 0;
   UInt_t endRow = 0;
   switch (padType) {
      case kShortPads:
         startRow = 0;
         endRow = lROCfitted->GetNrows();
         break;
      case kMediumPads:
         startRow = 0;
         endRow = 64;
         break;
      case kLongPads:
         startRow = 64;
         endRow = lROCfitted->GetNrows();
         break;
   }

   AliTPCFitPad::GetPadRegionCenterLocal(padType, centerPad);   
   Double_t value = 0;
   for (UInt_t irow = startRow; irow < endRow; irow++) {
      for (UInt_t ipad = 0; ipad < lROCfitted->GetNPads(irow); ipad++) {
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
         lROCfitted->SetValue(irow, ipad, value);
      }
   }
   return lROCfitted;
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

Bool_t AliTPCcalibTracksGain::GetParameters(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD &fitParam) {
   //
   // Puts the fit parameters for the specified segment (IROC & OROC), padType and fitType
   // into the fitParam TVectorD (which should contain 8 elements).
   // padType is one of kShortPads, kMediumPads, kLongPads. fitType is one of kSimpleFitter, kSqrtFitter, kLogFitter.
   // Note: The fitter has to be evaluated first!
   //
  TLinearFitter * fitter = GetFitter(segment, padType, fitType);
  if (fitter){
    fitter->Eval();
    fitter->GetParameters(fitParam);
    return kTRUE;
  }else{
    Error("AliTPCcalibTracksGain::GetParameters","Fitter%d_%d_%d not available", segment, padType, fitType);
    return kFALSE;
  }
  return kFALSE;
}

void AliTPCcalibTracksGain::GetErrors(UInt_t segment, UInt_t padType, UInt_t fitType, TVectorD &fitError) {
   //
   // Puts the fit parameter errors for the specified segment (IROC & OROC), padType and fitType
   // into the fitError TVectorD (which should contain 8 elements).
   // padType is one of kShortPads, kMediumPads, kLongPads. fitType is one of kSimpleFitter, kSqrtFitter, kLogFitter.
   // Note: The fitter has to be evaluated first!
   //

   GetFitter(segment, padType, fitType)->GetErrors(fitError);
   //fitError *= TMath::Sqrt(GetRedChi2(segment, padType, fitType));
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
      case 3:
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

void AliTPCcalibTracksGain::DumpTrack(AliTPCseed* track) {
   //
   //  Dump track information to the debug stream
   //
   
   Int_t rows[200];
   Int_t sector[3];
   Int_t npoints[3];
   TVectorD dedxM[3];
   TVectorD dedxQ[3];
   TVectorD parY[3];
   TVectorD parZ[3];
   TVectorD meanPos[3];
   //
   Int_t count=0;
   for (Int_t ipad = 0; ipad < 3; ipad++) {
     dedxM[ipad].ResizeTo(5);
     dedxQ[ipad].ResizeTo(5);
     parY[ipad].ResizeTo(3);
     parZ[ipad].ResizeTo(3);
     meanPos[ipad].ResizeTo(6);
     Bool_t isOK = GetDedx(track, ipad, rows, sector[ipad], npoints[ipad], dedxM[ipad], dedxQ[ipad], parY[ipad], parZ[ipad], meanPos[ipad]);
     if (isOK) 
       AddTracklet(sector[ipad],ipad, dedxQ[ipad], dedxM[ipad], parY[ipad], parZ[ipad], meanPos[ipad] );
     if (isOK) count++;
   }

   TTreeSRedirector * cstream =  GetDebugStreamer();
   if (cstream){
     (*cstream) << "Track" <<
       "run="<<fRun<<              //  run number
       "event="<<fEvent<<          //  event number
       "time="<<fTime<<            //  time stamp of event
       "trigger="<<fTrigger<<      //  trigger
       "mag="<<fMagF<<             //  magnetic field	      
       "Track.=" << track <<        // track information
       "\n";
     //
     //
     //
     if ( GetStreamLevel()>1 && count>1){
       (*cstream) << "TrackG" <<
	 "run="<<fRun<<              //  run number
	 "event="<<fEvent<<          //  event number
	 "time="<<fTime<<            //  time stamp of event
	 "trigger="<<fTrigger<<      //  trigger
	 "mag="<<fMagF<<             //  magnetic field	      
	 "Track.=" << track <<        // track information
	 "ncount="<<count<<
	 // info for pad type 0
	 "sector0="<<sector[0]<<     
	 "npoints0="<<npoints[0]<<
	 "dedxM0.="<<&dedxM[0]<<
	 "dedxQ0.="<<&dedxQ[0]<<
	 "parY0.="<<&parY[0]<<
	 "parZ0.="<<&parZ[0]<<
	 "meanPos0.="<<&meanPos[0]<<
	 //
	 // info for pad type 1
	 "sector1="<<sector[1]<<     
	 "npoints1="<<npoints[1]<<
	 "dedxM1.="<<&dedxM[1]<<
	 "dedxQ1.="<<&dedxQ[1]<<
	 "parY1.="<<&parY[1]<<
	 "parZ1.="<<&parZ[1]<<
	 "meanPos1.="<<&meanPos[1]<<
	 //
	 // info for pad type 2
	 "sector2="<<sector[2]<<     
	 "npoints2="<<npoints[2]<<
	 "dedxM2.="<<&dedxM[2]<<
	 "dedxQ2.="<<&dedxQ[2]<<
	 "parY2.="<<&parY[2]<<
	 "parZ2.="<<&parZ[2]<<
	 "meanPos2.="<<&meanPos[2]<<
	 //       
	 "\n";
     }
   }
}

Bool_t AliTPCcalibTracksGain::GetDedx(AliTPCseed* track, Int_t padType, Int_t* /*rows*/,
				      Int_t &sector, Int_t& npoints, 
				      TVectorD &dedxM, TVectorD &dedxQ, 
				      TVectorD &parY, TVectorD &parZ, TVectorD&meanPos)
{
  //
  // GetDedx for given sector for given track
  // padType - type of pads
  //

  static TLinearFitter fitY(2, "pol1");
  static TLinearFitter fitZ(2, "pol1");
  fitY.StoreData(kFALSE);
  fitZ.StoreData(kFALSE);
  //
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
   //
   fitY.ClearPoints();
   fitZ.ClearPoints();
   
   for (Int_t iCluster = firstRow; iCluster < lastRow; iCluster++) {
      AliTPCclusterMI* cluster = track->GetClusterPointer(iCluster);
      if (cluster) {
         Int_t detector = cluster->GetDetector() ;
         if (lastSector == -1) lastSector = detector;
         if (lastSector != detector) continue;
         amplitudeQ[nclusters] = cluster->GetQ()/GetQNorm(cluster,parY[1], parZ[1]);
         amplitudeM[nclusters] = cluster->GetMax()/GetMaxNorm(cluster,parY[1], parZ[1]);
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
   //
   //
   //
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
   TTreeSRedirector * cstream =  GetDebugStreamer();   
   inonEdge = 0;
   Float_t momenta = track->GetP();
   Float_t mdedx = track->GetdEdx();
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

      AddCluster(cluster, momenta, mdedx, padType, xcenter, dedxQ, dedxM, fraction, fraction2, dedge, parY, parZ, meanPos);
      
      Double_t qNorm = GetQNorm(cluster,parY[1], parZ[1]);
      Double_t mNorm = GetMaxNorm(cluster,parY[1], parZ[1]);


      Float_t gain = GetGain(cluster);
      if (cstream) (*cstream) << "dEdx" <<
	"run="<<fRun<<              //  run number
	"event="<<fEvent<<          //  event number
	"time="<<fTime<<            //  time stamp of event
	"trigger="<<fTrigger<<      //  trigger
	"mag="<<fMagF<<             //  magnetic field	      

	"Cl.=" << cluster <<           // cluster of interest
	"gain="<<gain<<                // gain at cluster position
	"mNorm="<<mNorm<<              // Q max normalization
	"qNorm="<<qNorm<<              // Q tot normalization
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
   
   if (cstream) (*cstream) << "dEdxT" <<
     "run="<<fRun<<              //  run number
     "event="<<fEvent<<          //  event number
     "time="<<fTime<<            //  time stamp of event
     "trigger="<<fTrigger<<      //  trigger
     "mag="<<fMagF<<             //  magnetic field	      
     "P=" << momenta <<             // track momenta
     "npoints="<<inonEdge<<         // number of points
     "sector="<<lastSector<<        // sector number
     "dedx=" << mdedx <<            // mean dedx - corrected for angle
     "IPad=" << padType <<          // pad type 0..2
     "xc=" << xcenter <<            // x center of chamber
     "dedxQ.=" << &dedxQ <<         // dedxQ  - total charge
     "dedxM.=" << &dedxM <<         // dedxM  - maximal charge
     "parY.=" << &parY <<           // line fit
     "parZ.=" << &parZ <<           // line fit
     "meanPos.=" << &meanPos <<     // mean position (dx, dx^2, y,y^2, z, z^2)
     "\n";
   
   sector = lastSector;
   npoints = inonEdge;
   return kTRUE;
}

void AliTPCcalibTracksGain::AddTracklet(UInt_t sector, UInt_t padType,TVectorD &dedxQ, TVectorD &dedxM,TVectorD& parY, TVectorD& parZ, TVectorD& meanPos){
  //
  // Add measured point - dedx to the fitter
  //
  //
  //chain->SetAlias("dr","(250-abs(meanPos.fElements[4]))/250");
  //chain->SetAlias("tz","(0+abs(parZ.fElements[1]))");
  //chain->SetAlias("ty","(0+abs(parY.fElements[1]))");
  //chain->SetAlias("corrg","sqrt((1+ty^2)*(1+tz^2))");
  //expession fast - TString *strq0 = toolkit.FitPlane(chain,"dedxQ.fElements[2]","dr++ty++tz++dr*ty++dr*tz++ty*tz++ty^2++tz^2","IPad==0",chi2,npoints,param,covar,0,100000);

  Double_t xxx[100];
  //
  // z and angular part
  //
 
  xxx[0] = (250.-TMath::Abs(meanPos[4]))/250.;
  xxx[1] = TMath::Abs(parY[1]);
  xxx[2] = TMath::Abs(parZ[1]);
  xxx[3] = xxx[0]*xxx[1];
  xxx[4] = xxx[0]*xxx[2];
  xxx[5] = xxx[1]*xxx[2];
  xxx[6] = xxx[0]*xxx[0];
  xxx[7] = xxx[1]*xxx[1];
  xxx[8] = xxx[2]*xxx[2];
  //
  // chamber part
  //
  Int_t tsector = sector%36;
  for (Int_t i=0;i<35;i++){
    xxx[9+i]=(i==tsector)?1:0;
  }
  TLinearFitter *fitterM = fFitter0M;
  if (padType==1) fitterM=fFitter1M;
  if (padType==2) fitterM=fFitter2M;
  fitterM->AddPoint(xxx,dedxM[1]);
  //
  TLinearFitter *fitterT = fFitter0T;
  if (padType==1) fitterT = fFitter1T;
  if (padType==2) fitterT = fFitter2T;
  fitterT->AddPoint(xxx,dedxQ[1]);
  //
  TLinearFitter *dfitterM = fDFitter0M;
  if (padType==1) dfitterM=fDFitter1M;
  if (padType==2) dfitterM=fDFitter2M;
  dfitterM->AddPoint(xxx,dedxM[1]);
  //
  TLinearFitter *dfitterT = fDFitter0T;
  if (padType==1) dfitterT = fDFitter1T;
  if (padType==2) dfitterT = fDFitter2T;
  dfitterT->AddPoint(xxx,dedxQ[1]);
}


TGraph *AliTPCcalibTracksGain::CreateAmpGraph(Int_t ipad, Bool_t qmax){
  //
  // create the amplitude graph
  // The normalized amplitudes are extrapolated to the 0 angle (y,z)  and 0 drift length
  //
  
  TVectorD vec;
  if (qmax){
    if (ipad==0) fFitter0M->GetParameters(vec);
    if (ipad==1) fFitter1M->GetParameters(vec);
    if (ipad==2) fFitter2M->GetParameters(vec);
  }else{
    if (ipad==0) fFitter0T->GetParameters(vec);
    if (ipad==1) fFitter1T->GetParameters(vec);
    if (ipad==2) fFitter2T->GetParameters(vec);
  }
  
  Float_t amp[36];
  Float_t sec[36];
  for (Int_t i=0;i<35;i++){
    sec[i]=i;
    amp[i]=vec[10+i]+vec[0];
  }
  amp[35]=vec[0];
  Float_t mean = TMath::Mean(36,amp);
  for (Int_t i=0;i<36;i++){
    sec[i]=i;
    amp[i]=(amp[i]-mean)/mean;
  }
  TGraph *gr = new TGraph(36,sec,amp);
  return gr;
}


void   AliTPCcalibTracksGain::UpdateClusterParam(AliTPCClusterParam* clparam){
  //
  //   SetQ normalization parameters
  //
  //  void SetQnorm(Int_t ipad, Int_t itype,  TVectorD * norm); 

  TVectorD vec;
  
  //
  fDFitter0T->Eval();
  fDFitter1T->Eval();
  fDFitter2T->Eval();
  fDFitter0M->Eval();
  fDFitter1M->Eval();
  fDFitter2M->Eval();
  fDFitter0T->GetParameters(vec);
  clparam->SetQnorm(0,0,&vec);
  fDFitter1T->GetParameters(vec);
  clparam->SetQnorm(1,0,&vec);
  fDFitter2T->GetParameters(vec);
  clparam->SetQnorm(2,0,&vec);
  //
  fDFitter0M->GetParameters(vec);
  clparam->SetQnorm(0,1,&vec);
  fDFitter1M->GetParameters(vec);
  clparam->SetQnorm(1,1,&vec);
  fDFitter2M->GetParameters(vec);
  clparam->SetQnorm(2,1,&vec);
  //

}


void   AliTPCcalibTracksGain::Analyze(){

 Evaluate();

}


