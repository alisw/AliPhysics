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

/* $Id$ */

//_________________________________________________________________________
// Manager class for TOF reconstruction.
// 
//
//-- Authors: Bologna-ITEP-Salerno Group
//
// Description: Manager class for TOF reconstruction (derived from TTask)
// Summary of the main methods:
// - extraction of the TPC (assumed to be) reconstructed tracks 
//   comment: it has to me moved as soon as possible into a separate
//   class AliTOFTrackReader (K. Safarik suggestion)
// - geometrical propagation of the above tracks till TOF detector
// - matching of the tracks with the TOF signals
// 
// Remark: the GEANT3.21 geometry is used during the geometrical propagation
// of the tracks in order to know the current volume reached by the track.
//
//////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <stdlib.h>

#include <TBenchmark.h>
#include <TClonesArray.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TFolder.h>
#include <TGeant3.h>
#include <TNtuple.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTask.h>
#include <TTree.h>
#include <TVirtualMC.h>

#include "AliConst.h"
#include "AliTOFGeometry.h"
#include "AliDetector.h"
#include "AliHeader.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliTOF.h"
#include "AliTOFHitMap.h"
#include "AliTOFPad.h"
#include "AliTOFRecHit.h"
#include "AliTOFReconstructioner.h"
#include "AliTOFSDigit.h"
#include "AliTOFTrack.h"
#include "AliTOFhit.h"
#include "AliTOFv1.h"
#include "AliTOFv2.h"
#include "AliTOFv3.h"
#include "AliTOFv4.h"
#include "AliTOFv4T0.h"
#include "AliMC.h"

// #include "../TPC/AliTPC.h"
// AliTPChit class or somewhere
// this line has to be commented till TPC will provide fPx fPy fPz and fL in

ClassImp(AliTOFReconstructioner)

//____________________________________________________________________________ 
  AliTOFReconstructioner::AliTOFReconstructioner():TTask("AliTOFReconstructioner","") 
{
  // default ctor
  fNevents = 0 ; 
  foutputfile  = 0; 
  foutputntuple= 0;
  fZnoise  = 0;
  ftail    = 0;
}
           
//____________________________________________________________________________ 
  AliTOFReconstructioner::AliTOFReconstructioner(char* headerFile, Option_t* opt, char *RecFile ):TTask("AliTOFReconstructioner","") 
{
  //
  // ctor
  //
  fNevents = 0 ;     // Number of events to reconstruct, 0 means all evens in current file
  foutputfile  = 0; 
  foutputntuple= 0;
  fZnoise  = 0;
  ftail    = 0;

  Init(opt);

  // create output file
  if (RecFile){
    foutputfile= new TFile(RecFile,"RECREATE","root file for matching");
  } else {
    char outFileName[100];
    strcpy(outFileName,"match");
    strcat(outFileName,headerFile);
    foutputfile= new TFile(outFileName,"RECREATE","root file for matching");
  }
  
  // initialize the ALIROOT geometry 
  gAlice->Init();
  gAlice->Print(); 

  CreateNTuple();  

  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
}
//____________________________________________________________________________ 
void AliTOFReconstructioner::Init(Option_t* opt)
{
  // Initialize the AliTOFReconstructioner setting parameters for
  // reconstruction.
  // Option values: Pb-Pb for Pb-Pb events
  //                pp    for pp    events

  // set common parameters
  fdbg=1;
  fNevents    = 1;
  fFirstEvent = 1;
  fLastEvent  = 1;
  fTimeResolution =0.120;
  fpadefficiency  =0.99 ;
  fEdgeEffect     = 2   ;
  fEdgeTails      = 0   ;
  fHparameter     = 0.4 ;
  fH2parameter    = 0.15;
  fKparameter     = 0.5 ;
  fK2parameter    = 0.35;
  fEffCenter      = fpadefficiency;
  fEffBoundary    = 0.65;
  fEff2Boundary   = 0.90;
  fEff3Boundary   = 0.08;
  fResCenter      = 50. ;
  fResBoundary    = 70. ;
  fResSlope       = 40. ;
  fTimeWalkCenter = 0.  ;
  fTimeWalkBoundary=0.  ;
  fTimeWalkSlope  = 0.  ;
  fTimeDelayFlag  = 1   ;
  fPulseHeightSlope=2.0 ;
  fTimeDelaySlope =0.060;
  // was fMinimumCharge = TMath::Exp(fPulseHeightSlope*fKparameter/2.);
  fMinimumCharge = TMath::Exp(-fPulseHeightSlope*fHparameter);
  fChargeSmearing=0.0   ;
  fLogChargeSmearing=0.13;
  fTimeSmearing   =0.022;
  fAverageTimeFlag=0    ;
  fChargeFactorForMatching=1;
  fTrackingEfficiency=1.0; // 100% TPC tracking efficiency assumed
  fSigmavsp = 1.        ;
  fSigmaZ   = 0.        ;
  fSigmarphi= 0.        ;
  fSigmap   = 0.        ;
  fSigmaPhi = 0.        ;
  fSigmaTheta=0.        ;
  fField    = 0.2       ;
  // fRadLenTPC : 0.2 includes TRD / 0.03 TPC only
  fRadLenTPC=0.06        ; // last value
  fCorrectionTRD=0.     ;
  fLastTPCRow=111       ;
  fRadiusvtxBound=50.   ; // expressed in [cm]
  fStep     = 0.1       ; // expressed in [cm] step during propagation of the
                          // track inside TOF volumes 
  fMatchingStyle=2      ;
  /* previous values default
  fMaxPixels=70000      ;
  fMaxAllTracks=70000   ;
  fMaxTracks=15000      ;
  */
  fMaxPixels=165000      ;
  fMaxAllTracks=500000   ;
  fMaxTracks=15000      ;

  fMaxTOFHits=35000     ;
  fPBound      =0.0     ; // bending effect: P_t=0.3*z*B*R , z particle charge
  fNoiseSlope=20.       ;
  // set parameters as specified in opt
  //pp case
  if(strstr(opt,"pp")){
  fMaxTestTracks=500    ; 
  fNoise    = 26.       ;
  fNoiseMeanTof= 26.4       ; // to check
  }
  //Pb-Pb case
  if(strstr(opt,"Pb-Pb")){
  fMaxTestTracks=20     ;
  fNoise    = 9400.     ;
  fNoiseMeanTof= 26.4       ;
  }
}

//____________________________________________________________________________ 
  AliTOFReconstructioner::~AliTOFReconstructioner()
{
  //
  // dtor
  //

  if (foutputfile)
    {
      delete foutputfile;
      foutputfile = 0;
    }
  if (foutputntuple)
    {
      delete foutputntuple;
      foutputntuple = 0;
    }

  if (fZnoise)
    {
      delete fZnoise;
      fZnoise = 0;
    }

  if (ftail)
    {
      delete ftail;
      ftail = 0;
    }
}

//____________________________________________________________________________
void AliTOFReconstructioner::CreateNTuple()
{
  //
  // Create a Ntuple where information about reconstructed charged particles 
  // (both primaries and secondaries) are stored 
  // Variables: event ipart imam xvtx yvtx zvtx pxvtx pyvtx pzvtx time leng matc text mext
  // Meaning:
  // event - event number (0, 1, ...)
  // ipart - PDG code of particles 
  // imam  - PDG code for the parent
  // =0 for primary particle
  // xvtx  - x-coordinate of the vertex (cm)
  // yvtx  - y-coordinate of the vertex (cm)
  // zvtx  - z-coordinate of the vertex (cm)
  // pxvtx - x-coordinate of the momentum in the vertex (GeV)
  // pyvtx - y-coordinate of the momentum in the vertex (GeV)
  // pzvtx - z-coordinate of the momentum in the vertex (GeV)
  // time  - time of flight from TOF for given track (ps) - TOF time for the
  //         first TOF hit of the track
  // leng  - track length to the TOF pixel (cm), evaluate as a sum of the
  // track length from the track vertex to TPC and the average
  // length of the extrapolated track from TPC to TOF.
  // for the track without TOF hits leng=-abs(leng)
  // matc  - index of the (TPC track) - (TOF pixel) matching
  // =0 for tracks which are not tracks for matching, i.e. 
  // there is not hit on the TPC or Rvxt>200 cm
  // >0 for tracks with positive matching procedure:
  //   =1 or 2 for non-identified tracks:
  //     =1, if the corresponding pixel is not fired,
  //     =2, if the corresponding pixel is also matched to the 
  //         other track,
  //   =3 or 4 for identified tracks:
  //     =3, if identified with true time,
  //     =4, if identified with wrong time.
  // <0 for tracks with negative mathing procedure:
  //   =-1, if track do not reach the pixel plate (curved in the 
  //        magnetic field),
  //   =-2, if track is out of z-size of the TOF,
  //   =-3, if track is or into the RICH hole, or into the PHOS hole, or in the space between the plates,
  //   =-4, if track is into the dead space of the TOF.
  // text  - time of fligth from the matching procedure = time of the 
  //         pixel corresponding to the track (ps)
  //         =0 for the tracks with matc<=1
  // mext  - mass of the track from the matching procedure
  //           =p*sqrt(900*(text/leng)**2-1), if 900*(text/leng)**2-1>=0
  //           =-p*sqrt(abs(900*(text/leng)**2-1)), if 900*(text/leng)**2-1<0

  foutputntuple= new TNtuple("Ntuple","matching","event:ipart:imam:xvtx:yvtx:zvtx:pxvtx:pyvtx:pzvtx:time:leng:matc:text:mext",2000000); // buffersize set for 25 Pb-Pb events
}

//__________________________________________________________________
Double_t TimeWithTailR(Double_t* x, Double_t* par)
{
  // sigma - par[0], alpha - par[1], part - par[2]
  //  at x<part*sigma - gauss
  //  at x>part*sigma - TMath::Exp(-x/alpha)
  Float_t xx =x[0];
  Double_t f;
  if(xx<par[0]*par[2]) {
    f = TMath::Exp(-xx*xx/(2*par[0]*par[0]));
  } else {
    f = TMath::Exp(-(xx-par[0]*par[2])/par[1]-0.5*par[2]*par[2]);
  }
  return f;
}

//____________________________________________________________________________
void AliTOFReconstructioner::Exec(const char* datafile, Option_t *option) 
{ 
  //
  // Performs reconstruction for TOF detector
  // 
  gBenchmark->Start("TOFReconstruction");

  
  AliRunLoader *rl = AliRunLoader::Open(datafile);
  if (rl == 0x0)
   {
     Error("Exec","Can not open session for file %s",datafile);
     return;
   }
  // Get AliRun object from file or create it if not on file
  rl->LoadgAlice();
  gAlice = rl->GetAliRun();

  AliTOF* TOF = (AliTOF *) gAlice->GetDetector ("TOF");
  AliDetector* TPC = gAlice->GetDetector("TPC");

  if (!TOF) {
    Error("AliTOFReconstructioner","TOF not found");
    delete rl;
    return;
  }
  if (!TPC) {
    Error("AliTOFReconstructioner","TPC Detector not found");
    delete rl;
    return;
  }
  AliLoader* tpcloader = rl->GetLoader("TPCLoader");
  if (tpcloader == 0x0)
   {
    Error("AliTOFReconstructioner","Can not get TPC Loader from Run Loader.");
    delete rl;
    return;
   }

  AliLoader* tofloader = rl->GetLoader("TOFLoader");
  if (tofloader == 0x0)
   {
    Error("AliTOFReconstructioner","Can not get TOF Loader from Run Loader.");
    delete rl;
    return;
   }
  
  if (fEdgeTails) ftail = new TF1("tail",TimeWithTailR,-2,2,3);
  
  if (fNevents == 0) fNevents = rl->GetNumberOfEvents();
  // You have to set the number of event with the ad hoc setter
  // see testrecon.C
  if (rl->GetHeader() == 0x0) rl->LoadHeader();
  
  tofloader->LoadHits();
  tpcloader->LoadHits(); 
  
  for (Int_t ievent = 0; ievent < fNevents; ievent++) { // start loop on events
    rl->GetEvent(ievent);
    Int_t nparticles= rl->GetHeader()->GetNtrack();
    if (nparticles <= 0) return;

    TClonesArray* tofhits=0;
    TClonesArray* tpchits=0;

    if (TOF) tofhits = TOF->Hits();
    if (TPC) tpchits = TPC->Hits();

    TTree *TH = tofloader->TreeH();
    if (!TH) return;
    Int_t ntracks    = (Int_t) (TH->GetEntries()); // primary tracks
    cout << "number of primary tracked tracks in current event " << ntracks << endl; // number of primary tracked tracks
    // array declaration and initialization
    // TOF arrays
    //    Int_t mapPixels[AliTOFGeometry::NSectors()*AliTOFGeometry::NPlates()][AliTOFGeometry::NStripC()][AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX()];

    Int_t *** mapPixels = new Int_t**[AliTOFGeometry::NSectors()*AliTOFGeometry::NPlates()];
    for (Int_t i=0; i<AliTOFGeometry::NSectors()*AliTOFGeometry::NPlates(); i++) mapPixels[i] = new Int_t*[AliTOFGeometry::NStripC()];
    for (Int_t i=0; i<AliTOFGeometry::NSectors()*AliTOFGeometry::NPlates(); i++) {
      for (Int_t j=0; j<AliTOFGeometry::NStripC(); j++) {
        mapPixels[i][j]= new Int_t[AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX()];
      }
    }


    // initializing the previous array
    for (Int_t i=0;i<AliTOFGeometry::NSectors()*AliTOFGeometry::NPlates();i++) {
      for (Int_t j=0;j<AliTOFGeometry::NStripC();j++) {
        for (Int_t l=0;l<AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX();l++) {
          mapPixels[i][j][l]=0;
        }
      } 
    }

    Float_t * toftime = new Float_t[fMaxAllTracks]; 
    InitArray(toftime, fMaxAllTracks);
    AliTOFPad* pixelArray = new AliTOFPad[fMaxPixels];
    Int_t* iTOFpixel        = new Int_t[fMaxAllTracks];
    InitArray(iTOFpixel   , fMaxAllTracks);
    Int_t* kTOFhitFirst     = new Int_t[fMaxAllTracks];
    InitArray(kTOFhitFirst, fMaxAllTracks);
    AliTOFRecHit* hitArray  = new AliTOFRecHit[fMaxTOFHits];
    Int_t isHitOnFiredPad=0; // index used to fill hitArray (array used to store informations
                             // about pads that contains an hit)
    Int_t ntotFiredPads=0;   // index used to fill array -> total number of fired pads (at least one time)

    // TPC arrays
    AliTOFTrack* trackArray = new AliTOFTrack[fMaxTracks];
    Int_t * iparticle = new Int_t[fMaxAllTracks];
    InitArray(iparticle,fMaxAllTracks); 
    Int_t * iTrackPt  = new Int_t[fMaxTracks];
    InitArray(iTrackPt, fMaxTracks);  // array 
    Float_t * ptTrack = new Float_t[fMaxTracks];
    InitArray( ptTrack, fMaxTracks);  // array for selected track pt  
    Int_t   ntotTPCtracks=0; // total number of selected TPC tracks

    
    // reading TOF hits
    if(TOF) ReadTOFHits(ntracks, TH, tofhits, mapPixels, kTOFhitFirst, pixelArray, iTOFpixel, toftime, hitArray,isHitOnFiredPad,ntotFiredPads);
    cout << "isHitOnFiredPad " << isHitOnFiredPad << " for event " << ievent << endl;

    // start debug for adding noise
    // adding noise
    Int_t nHitsNoNoise=isHitOnFiredPad;

    
    if(fNoise) AddNoiseFromOuter(option,mapPixels,pixelArray,hitArray,isHitOnFiredPad,ntotFiredPads);
    cout << "ntotFiredPads after adding noise  " << ntotFiredPads   << " for event " << ievent << endl;
    // set the hitArray distance to nearest hit
    SetMinDistance(hitArray,nHitsNoNoise);

    // these lines has to be commented till TPC will provide fPx fPy fPz 
    // and fL in AliTPChit class
    // reading TPC hits
    /*
    if(TPC) ReadTPCHits(ntracks, TH, tpchits, iTrackPt, iparticle, ptTrack, trackArray,ntotTPCtracks);
    */
    
    // geometrical matching
    if(TOF && TPC) Matching(trackArray,hitArray,mapPixels,pixelArray,kTOFhitFirst,ntotFiredPads,iTrackPt,iTOFpixel,ntotTPCtracks);
    
    // fill ntuple with reconstructed particles from current event
    FillNtuple(ntracks,trackArray,hitArray,pixelArray,iTOFpixel,iparticle,toftime,ntotFiredPads,ntotTPCtracks);
    

    // free used memory
    delete [] toftime;
    delete [] pixelArray;
    delete [] iTOFpixel;
    delete [] kTOFhitFirst;
    delete [] hitArray;
    delete [] trackArray;
    delete [] iparticle;
    delete [] iTrackPt;
    delete [] ptTrack;

   for (Int_t i=0; i<AliTOFGeometry::NSectors()*AliTOFGeometry::NPlates(); i++) {
      for (Int_t j=0; j<AliTOFGeometry::NStripC(); j++) {
        delete [] mapPixels[i][j];
      }
    }
    for (Int_t i=0; i<AliTOFGeometry::NSectors()*AliTOFGeometry::NPlates(); i++) delete [] mapPixels[i];

    delete [] mapPixels;

  }//event loop

  // free used memory for ftail
  if (ftail)
    {
      delete ftail;
      ftail = 0;
    }

  // writing ntuple on output file
  foutputfile->cd();
  //foutputntuple->Write(0,TObject::kOverwrite);
  foutputntuple->Write();
  foutputfile->Write();
  foutputfile->Close();

  gBenchmark->Stop("TOFReconstruction");
  cout << "AliTOFReconstructioner:" << endl ;
  cout << "   took " << gBenchmark->GetCpuTime("TOFReconstruction") << " seconds in order to make the reconstruction for " <<  fNevents << " events " << endl;
  cout <<  gBenchmark->GetCpuTime("TOFReconstruction")/fNevents << " seconds per event " << endl ;
  cout << endl ;
  
}
 
//__________________________________________________________________
void AliTOFReconstructioner::SetRecFile(char * file )
{
  //
  // Set the file name for reconstruction output 
  //
  if(!fRecFile.IsNull())
    cout << "Changing destination file for TOF reconstruction from " <<(char *)fRecFile.Data() << " to " << file << endl ;
  fRecFile=file ;
}
//__________________________________________________________________
void AliTOFReconstructioner::Print(Option_t* /*option*/)const
{
  //
  // Print reconstruction output file name
  //
  cout << "------------------- "<< GetName() << " -------------" << endl ;
  if(fRecFile.IsNull())
    cout << " Writing reconstructed particles to file galice.root "<< endl ;
  else
    cout << "    Writing reconstructed particle to file  " << (char*) fRecFile.Data() << endl ;

}

//__________________________________________________________________
void AliTOFReconstructioner::PrintParameters()const
{
  //
  // Print parameters used for reconstruction
  //
  cout << " ------------------- "<< GetName() << " -------------" << endl ;
  cout << " Parameters used for TOF reconstruction " << endl ;
  //  Printing the parameters
  
  cout << " Number of events:                        " << fNevents << endl; 
  cout << " Recostruction from event                 "<< fFirstEvent << "  to event "<< fLastEvent << endl;
  cout << " TOF geometry parameters                  " << endl;
  cout << " Min. radius of the TOF (cm)              "<< AliTOFGeometry::Rmin() << endl;
  cout << " Max. radius of the TOF (cm)              "<< AliTOFGeometry::Rmax() << endl;
  cout << " Number of TOF geom. levels               "<< AliTOFGeometry::MaxTOFTree()<< endl;
  cout << " Number of TOF sectors                    "<< AliTOFGeometry::NSectors() << endl;
  cout << " Number of TOF modules                    "<< AliTOFGeometry::NPlates() << endl;
  cout << " Max. Number of strips in a module        "<< AliTOFGeometry::NStripC() << endl;
  cout << " Number of pads per strip                 "<< AliTOFGeometry::NpadX()*AliTOFGeometry::NpadZ() << endl;
  cout << " Number of strips in central module       "<< AliTOFGeometry::NStripA() << endl;
  cout << " Number of strips in intermediate modules "<< AliTOFGeometry::NStripB() << endl;
  cout << " Number of strips in outer modules        "<< AliTOFGeometry::NStripC() << endl;
  cout << " Number of MRPC in x strip direction      "<< AliTOFGeometry::NpadX()<< endl;
  cout << " Size of MRPC (cm) along X                "<< AliTOFGeometry::XPad()<< endl;
  cout << " Number of MRPC in z strip direction      "<< AliTOFGeometry::NpadZ()<<endl;
  cout << " Size of MRPC (cm) along Z                "<< AliTOFGeometry::ZPad()<<endl;
  cout << " Module Lengths (cm)" << endl;
  cout << " A Module: "<< AliTOFGeometry::ZlenA()<< "  B Modules: "<< AliTOFGeometry::ZlenB()<< "  C Modules: "<< AliTOFGeometry::ZlenC()<< endl;
  cout << " Inner radius of the TOF detector (cm): "<<AliTOFGeometry::Rmin() << endl;
  cout << " Outer radius of the TOF detector (cm): "<<AliTOFGeometry::Rmax() << endl;
  cout << " Max. half z-size of TOF (cm)         : "<<AliTOFGeometry::MaxhZtof() << endl;
  cout << " TOF Pad parameters   " << endl;
  cout << " Time Resolution (ns) "<< fTimeResolution <<" Pad Efficiency: "<< fpadefficiency << endl;
  cout << " Edge Effect option:  "<<  fEdgeEffect<< endl;

  cout << " Boundary Effect Simulation Parameters " << endl;
  cout << " Hparameter: "<< fHparameter<<"  H2parameter:"<< fH2parameter <<"  Kparameter:"<< fKparameter<<"  K2parameter: "<< fK2parameter << endl;
  cout << " Efficiency in the central region of the pad: "<< fEffCenter << endl;
  cout << " Efficiency at the boundary region of the pad: "<< fEffBoundary << endl;
  cout << " Efficiency value at H2parameter "<< fEff2Boundary << endl;
  cout << " Efficiency value at K2parameter "<< fEff3Boundary << endl;
  cout << " Resolution (ps) in the central region of the pad: "<< fResCenter << endl;
  cout << " Resolution (ps) at the boundary of the pad      : "<< fResBoundary << endl;
  cout << " Slope (ps/K) for neighbouring pad               : "<< fResSlope <<endl;
  cout << " Time walk (ps) in the central region of the pad : "<< fTimeWalkCenter << endl;
  cout << " Time walk (ps) at the boundary of the pad       : "<< fTimeWalkBoundary<< endl;
  cout << " Slope (ps/K) for neighbouring pad               : "<< fTimeWalkSlope<<endl;
  cout << " Pulse Heigth Simulation Parameters " << endl;
  cout << " Flag for delay due to the PulseHeightEffect: "<< fTimeDelayFlag <<endl;
  cout << " Pulse Height Slope                           : "<< fPulseHeightSlope<<endl;
  cout << " Time Delay Slope                             : "<< fTimeDelaySlope<<endl;
  cout << " Minimum charge amount which could be induced : "<< fMinimumCharge<<endl;
  cout << " Smearing in charge in (q1/q2) vs x plot      : "<< fChargeSmearing<<endl;
  cout << " Smearing in log of charge ratio              : "<< fLogChargeSmearing<<endl;
  cout << " Smearing in time in time vs log(q1/q2) plot  : "<< fTimeSmearing<<endl;
  cout << " Flag for average time                        : "<< fAverageTimeFlag<<endl;
  cout << " Charge factor flag for matching              : "<< fChargeFactorForMatching<<endl;
  cout << " Edge tails option                            : "<< fEdgeTails << endl;
  cout << " TPC tracking  parameters " << endl;
  cout << " TPC tracking efficiency                      : "<< fTrackingEfficiency<< endl;
  cout << " Sigma vs momentum dependency flag            : "<< fSigmavsp << endl;
  cout << " Space uncertainties (cm). sigma(z) (cm): "<< fSigmaZ << " sigma(R(phi)) (cm): "<< fSigmarphi << endl;
  cout << " Momentum uncertainties.   sigma(delta(P)/P): "<< fSigmap <<" sigma(phi) (rad): "<< fSigmaPhi <<" sigma(theta) (rad): "<< fSigmaTheta << endl;   
  cout << " Parameters for additional noise hits " << endl;
  cout << " Number of noise hits : " << fNoise <<" Slope parameter (ns) in the time distribution: " << fNoiseSlope << endl;
  cout << " Mean TOF for noise from outer regions (ns)" <<  fNoiseMeanTof << endl;
  cout << " Physical parameters " << endl;
  cout << " Magnetic Field (tesla)                   : "<< fField <<endl;
  cout << " Radiation length of the outer wall of TPC: "<< fRadLenTPC << endl;
  cout << " (TPC tracks)-(TOF pads) matching parameters " << endl;
  cout << " TRD Correction flag       : "<< fCorrectionTRD <<endl;
  cout << " Number of the last TPC row: "<< fLastTPCRow <<" Vertex radius (cm) for selected tracks: "<<fRadiusvtxBound<<endl;
  cout << " Max. number of test tracks: "<<fMaxTestTracks << endl;
  cout << " Space step (cm)           : "<< fStep <<endl;
  cout << " Matching style option     : "<< fMatchingStyle <<endl;
  cout << " Array parameters " << endl;
  cout << " Max.number of pads involved in the matching procedure: "<< fMaxPixels << endl;
  cout << " Max.number of TOF hits per event                     : "<< fMaxTOFHits<< endl;
  cout << " Max.number of tracks selected for matching           : "<< fMaxTracks << endl;
  cout << " Max.number of all tracks including the neutral ones  : "<< fMaxAllTracks<< endl;
  cout << " Debug Flag                                           : "<< fdbg << endl;
  cout << " Cut on momentum for selecting tracks                 : "<< fPBound << endl;
  
}

//__________________________________________________________________
void AliTOFReconstructioner::IsInsideThePad(TVirtualMC *vmc, Float_t x, Float_t y, Float_t z, Int_t *nGeom, Float_t& zPad, Float_t& xPad) 
{
  //   input: x,y,z - coordinates of a hit
  //   output: array  nGeom[]
  //          nGeom[0] - the TOF sector number, 1,2,...,18 along azimuthal direction starting from -90 deg.!!!
  //          nGeom[1] - the TOF module number, 1,2,3,4,5=C,B,A,B,C along z-direction
  //          nGeom[2] - the TOF strip  number, 1,2,... along z-direction
  //          nGeom[3] - the TOF padz  number,  1,2=NPZ across a strip
  //          nGeom[4] - the TOF padx  number,  1,2,...,48=NPX along a strip
  //          zPad, xPad - coordinates of the hit in the pad frame
  //  numbering is adopted for the version 3.05 of AliRoot
  //  example:
  //   from Hits: sec,pla,str,padz,padx=4,2,14,2,35
  //  Vol. n.0: ALIC, copy number 1
  //  Vol. n.1: B077, copy number 1
  //  Vol. n.2: B074, copy number 5
  //  Vol. n.3: BTO2, copy number 1
  //  Vol. n.4: FTOB, copy number 2
  //  Vol. n.5: FLTB, copy number 0
  //  Vol. n.6: FSTR, copy number 14
  //  Vol. n.7: FSEN, copy number 0
  //  Vol. n.8: FSEZ, copy number 2
  //  Vol. n.9: FSEX, copy number 35
  //  Vol. n.10: FPAD, copy number 0


  Float_t xTOF[3];
  Int_t sector=0,module=0,strip=0,padz=0,padx=0;
  Int_t i,numed,nLevel,copyNumber;
  Gcvolu_t* gcvolu;
  char name[5];
  name[4]=0;
  
  for (i=0; i<AliTOFGeometry::MaxTOFTree(); i++) nGeom[i]=0;
  zPad=100.;
  xPad=100.;
  
  xTOF[0]=x;
  xTOF[1]=y;
  xTOF[2]=z;
  
  TGeant3 * g3 = (TGeant3*) vmc;

  g3->Gmedia(xTOF, numed);
  gcvolu=g3->Gcvolu();
  nLevel=gcvolu->nlevel;
  if(fdbg) {
    for (Int_t i=0; i<nLevel; i++) {
      strncpy(name,(char*) (&gcvolu->names[i]),4);
      cout<<"Vol. n."<<i<<": "<<name<<", copy number "<<gcvolu->number[i]<<endl;
    }
  }
  if(nLevel>=2) {
    // sector type name: B071(1,2,...,10),B074(1,2,3,4,5-PHOS),B075(1,2,3-RICH)
    strncpy(name,(char*) (&gcvolu->names[2]),4);
    // volume copy: 1,2,...,10 for B071, 1,2,3,4,5 for B074, 1,2,3 for B075
    copyNumber=gcvolu->number[2];
   if(!strcmp(name,"B071")) {
     if (copyNumber>=6 && copyNumber<=8) {
       sector=copyNumber+10;
     } else if (copyNumber>=1 && copyNumber<=5){
       sector=copyNumber+7;
     } else {
       sector=copyNumber-8;
     }
   } else if(!strcmp(name,"B075")) {
     sector=copyNumber+12;
   } else if(!strcmp(name,"B074")) {
     if (copyNumber>=1 && copyNumber<=3){
       sector=copyNumber+4;
     } else {
       sector=copyNumber-1;
     }
   }
  }
  if(sector) {
    nGeom[0]=sector;
    if(nLevel>=4) {
      // we'll use the module value in z-direction:
      //                                    1    2    3    4    5
      // the module order in z-direction: FTOC,FTOB,FTOA,FTOB,FTOC
      // the module copy:                   2    2    0    1    1
      // module type name: FTOA, FTOB, FTOC
      strncpy(name,(char*) (&gcvolu->names[4]),4);
      // module copy:  
      copyNumber=gcvolu->number[4];
      if(!strcmp(name,"FTOC")) {
	if (copyNumber==2) {
	  module=1;
	} else {
	  module=5;
	}
      } else if(!strcmp(name,"FTOB")) {
	if (copyNumber==2) {
	  module=2;
	} else {
	  module=4;
	}
      } else if(!strcmp(name,"FTOA")) {
	module=3;
      }
    }
  }
  
  if(module) {
    nGeom[1]=module;
    if(nLevel>=6) {
      // strip type name: FSTR
      strncpy(name,(char*) (&gcvolu->names[6]),4);
      // strip copy:  
      copyNumber=gcvolu->number[6];
      if(!strcmp(name,"FSTR")) strip=copyNumber; 
    }
  }
  
  if(strip) {
    nGeom[2]=strip;
    if(nLevel>=8) {
      // padz type name: FSEZ
      strncpy(name,(char*) (&gcvolu->names[8]),4);
      // padz copy:  
      copyNumber=gcvolu->number[8];
      if(!strcmp(name,"FSEZ")) padz=copyNumber; 
    }
  }
  if(padz) {
    nGeom[3]=padz;
    if(nLevel>=9) {
      // padx type name: FSEX
      strncpy(name,(char*) (&gcvolu->names[9]),4);
      // padx copy:  
      copyNumber=gcvolu->number[9];
      if(!strcmp(name,"FSEX")) padx=copyNumber; 
    }
  }
  
  if(padx) {
    nGeom[4]=padx;
    zPad=gcvolu->glx[2];  // check here
    xPad=gcvolu->glx[0];  // check here
  }
  
  //   printf(" nGeom[0,1,2,3,4]=%i,%i,%i,%i,%i\n",nGeom[0],nGeom[1],nGeom[2],nGeom[3],nGeom[4]); 
}

//__________________________________________________________________
void AliTOFReconstructioner::EpMulScatt(Float_t& px, Float_t& py, Float_t& pz, Float_t& p, Float_t& theta)
{
  //   Momentum p  - before mult.scat.
  //   Momentum p2 - after mult.scat.
  //   THE0 - r.m.s. of deviation angle in plane
  //           (see RPP'96: Phys.Rev.D54 (1996) 134)
  
  Float_t pt,thex,they,tantx,tanty,p2px,p2py,p2pz,costhe,sinthe,cospsi,sinpsi,p2x,p2y,p2z,p2,g;
  
  pt=TMath::Sqrt(px*px+py*py);
  //   angles for p in the ' frame with Z'along p
  if(fMatchingStyle==1) {
    thex=theta*gRandom->Gaus();
    they=theta*gRandom->Gaus();
  } else {
    thex=3*(-theta+2*theta*gRandom->Rndm());
    they=3*(-theta+2*theta*gRandom->Rndm());
  }
  tantx=TMath::Tan(thex);
  tanty=TMath::Tan(they);
  
  //   p2p - p2 in the ' frame
  p2pz=p/TMath::Sqrt(1.+tantx*tantx+tanty*tanty);
  p2py=p2pz*tanty;
  p2px=p2pz*tantx;
  //   choose X'so that PHI=0 (see Il'in, Pozdnyak Analiticheskaya geometriya, 1968, c.88
  //   for Euler angles PSI, THETA (PHI=0)
  costhe=pz/p;
  sinthe=pt/p;
  cospsi=-py/pt;
  sinpsi=px/pt;
  //
  g=p2py*costhe-p2pz*sinthe;
  p2x=p2px*cospsi-g*sinpsi;
  p2y=p2px*sinpsi+g*cospsi;
  p2z=p2py*sinthe+p2pz*costhe;
  p2=TMath::Sqrt(p2x*p2x+p2y*p2y+p2z*p2z);
  
  //   Test angle
  g=(px*p2x+py*p2y+pz*p2z)/(p*p2);
  if(g>1) g=1;
  theta=TMath::ACos(g);
  px=p2x;
  py=p2y;
  pz=p2z;
  p=p2;
  
}

// std border effect algorithm
//__________________________________________________________________
void AliTOFReconstructioner::BorderEffect(Float_t z0, Float_t x0, Float_t geantTime, Int_t& nActivatedPads, Int_t& nFiredPads, Bool_t* isFired, Int_t* nPlace, Float_t* qInduced, Float_t* tofTime, Float_t& averageTime)
{
  // Input:  z0, x0 - hit position in the strip system (0,0 - center of the strip), cm
  //         geantTime - time generated by Geant, ns
  // Output: nActivatedPads - the number of pads activated by the hit (1 || 2 || 4)
  //         nFiredPads - the number of pads fired (really activated) by the hit (nFiredPads <= nActivatedPads)
  //         qInduced[iPad]- charge induced on pad, arb. units
  //                         this array is initialized at zero by the caller
  //         tofAfterSimul[iPad] - time calculated with edge effect algorithm, ns
  //                                   this array is initialized at zero by the caller
  //         averageTime - time given by pad hited by the Geant track taking into account the times (weighted) given by the pads fired for edge effect also.
  //                       The weight is given by the qInduced[iPad]/qCenterPad
  //                                   this variable is initialized at zero by the caller
  //         nPlace[iPad] - the number of the pad place, iPad = 0, 1, 2, 3
  //                                   this variable is initialized at zero by the caller
  //
  // Description of used variables:
  //         eff[iPad] - efficiency of the pad
  //         res[iPad] - resolution of the pad, ns
  //         timeWalk[iPad] - time walk of the pad, ns
  //         timeDelay[iPad] - time delay for neighbouring pad to hited pad, ns
  //         PadId[iPad] - Pad Identifier
  //                    E | F    -->   PadId[iPad] = 5 | 6
  //                    A | B    -->   PadId[iPad] = 1 | 2
  //                    C | D    -->   PadId[iPad] = 3 | 4
  //         nTail[iPad] - the tail number, = 1 for tailA, = 2 for tailB
  //         qCenterPad - charge extimated for each pad, arb. units
  //         weightsSum - sum of weights extimated for each pad fired, arb. units
  
  const Float_t kSigmaForTail[2] = {AliTOFGeometry::SigmaForTail1(),AliTOFGeometry::SigmaForTail2()}; //for tail                                                   
  Int_t iz = 0, ix = 0;
  Float_t dX = 0., dZ = 0., x = 0., z = 0.;
  Float_t h = fHparameter, h2 = fH2parameter, k = fKparameter, k2 = fK2parameter;
  Float_t effX = 0., effZ = 0., resX = 0., resZ = 0., timeWalkX = 0., timeWalkZ = 0.;
  Float_t logOfqInd = 0.;
  Float_t weightsSum = 0.;
  Int_t nTail[4]  = {0,0,0,0};
  Int_t padId[4]  = {0,0,0,0};
  Float_t eff[4]  = {0.,0.,0.,0.};
  Float_t res[4]  = {0.,0.,0.,0.};
  //  Float_t qCenterPad = fMinimumCharge * fMinimumCharge;
  Float_t qCenterPad = 1.;
  Float_t timeWalk[4]  = {0.,0.,0.,0.};
  Float_t timeDelay[4] = {0.,0.,0.,0.};
  
  nActivatedPads = 0;
  nFiredPads = 0;
  
  (z0 <= 0) ? iz = 0 : iz = 1;
  dZ = z0 + (0.5 * AliTOFGeometry::NpadZ() - iz - 0.5) * AliTOFGeometry::ZPad(); // hit position in the pad frame, (0,0) - center of the pad
  z = 0.5 * AliTOFGeometry::ZPad() - TMath::Abs(dZ);                               // variable for eff., res. and timeWalk. functions
  iz++;                                                                              // z row: 1, ..., AliTOFGeometry::NpadZ() = 2
  ix = (Int_t)((x0 + 0.5 * AliTOFGeometry::NpadX() * AliTOFGeometry::XPad()) / AliTOFGeometry::XPad());
  dX = x0 + (0.5 * AliTOFGeometry::NpadX() - ix - 0.5) * AliTOFGeometry::XPad(); // hit position in the pad frame, (0,0) - center of the pad
  x = 0.5 * AliTOFGeometry::XPad() - TMath::Abs(dX);                               // variable for eff., res. and timeWalk. functions;
  ix++;                                                                              // x row: 1, ..., AliTOFGeometry::NpadX() = 48
  
  ////// Pad A:
  nActivatedPads++;
  nPlace[nActivatedPads-1] = (iz - 1) * AliTOFGeometry::NpadX() + ix;
  qInduced[nActivatedPads-1] = qCenterPad;
  padId[nActivatedPads-1] = 1;
  
  if (fEdgeEffect == 0) {
    eff[nActivatedPads-1] = fEffCenter;
    if (gRandom->Rndm() < eff[nActivatedPads-1]) {
      nFiredPads = 1;
      res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + fResCenter * fResCenter); // 10400=30^2+20^2+40^2+50^2+50^2+50^2  ns;
      isFired[nActivatedPads-1] = kTRUE;
      tofTime[nActivatedPads-1] = gRandom->Gaus(geantTime + fTimeWalkCenter, res[0]);
      averageTime = tofTime[nActivatedPads-1];
    }
  } else {
     
    if(z < h) {
      if(z < h2) {
	effZ = fEffBoundary + (fEff2Boundary - fEffBoundary) * z / h2;
      } else {
	effZ = fEff2Boundary + (fEffCenter - fEff2Boundary) * (z - h2) / (h - h2);
      }
      resZ = fResBoundary + (fResCenter - fResBoundary) * z / h;
      timeWalkZ = fTimeWalkBoundary + (fTimeWalkCenter - fTimeWalkBoundary) * z / h;
      nTail[nActivatedPads-1] = 1;
    } else {
      effZ = fEffCenter;
      resZ = fResCenter;
      timeWalkZ = fTimeWalkCenter;
    }
    
    if(x < h) {
      if(x < h2) {
	effX = fEffBoundary + (fEff2Boundary - fEffBoundary) * x / h2;
      } else {
	effX = fEff2Boundary + (fEffCenter - fEff2Boundary) * (x - h2) / (h - h2);
      }
      resX = fResBoundary + (fResCenter - fResBoundary) * x / h;
      timeWalkX = fTimeWalkBoundary + (fTimeWalkCenter - fTimeWalkBoundary) * x / h;
      nTail[nActivatedPads-1] = 1;
    } else {
      effX = fEffCenter;
      resX = fResCenter;
      timeWalkX = fTimeWalkCenter;
    }
    
    (effZ<effX) ? eff[nActivatedPads-1] = effZ : eff[nActivatedPads-1] = effX;
    (resZ<resX) ? res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resX * resX) : res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resZ * resZ); // 10400=30^2+20^2+40^2+50^2+50^2+50^2  ns
    (timeWalkZ<timeWalkX) ? timeWalk[nActivatedPads-1] = 0.001 *  timeWalkZ : timeWalk[nActivatedPads-1] = 0.001 * timeWalkX; // ns


    ////// Pad B:
    if(z < k2) {
      effZ = fEffBoundary - (fEffBoundary - fEff3Boundary) * (z / k2);
    } else {
      effZ = fEff3Boundary * (k - z) / (k - k2);
    }
    resZ = fResBoundary + fResSlope * z / k;
    timeWalkZ = fTimeWalkBoundary + fTimeWalkSlope * z / k;
    
    if(z < k && z > 0) {
      if( (iz == 1 && dZ > 0) || (iz == 2 && dZ < 0) ) {
	nActivatedPads++;
	nPlace[nActivatedPads-1] = nPlace[0] + (3 - 2 * iz) * AliTOFGeometry::NpadX();
	eff[nActivatedPads-1] = effZ;
	res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resZ * resZ); // 10400=30^2+20^2+40^2+50^2+50^2+50^2 ns 
	timeWalk[nActivatedPads-1] = 0.001 * timeWalkZ; // ns
	nTail[nActivatedPads-1] = 2;
	if (fTimeDelayFlag) {
	  //	  qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * z / 2.);
	  //	  qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * z / 2.);
	  qInduced[nActivatedPads-1] = TMath::Exp(-fPulseHeightSlope * z);
	  logOfqInd = gRandom->Gaus(-fPulseHeightSlope * z, fLogChargeSmearing);
	  timeDelay[nActivatedPads-1] = gRandom->Gaus(-fTimeDelaySlope * logOfqInd, fTimeSmearing);
	} else {
	  timeDelay[nActivatedPads-1] = 0.;
	}
	padId[nActivatedPads-1] = 2;
      }
    }

    
    ////// Pad C, D, E, F:
    if(x < k2) {
      effX = fEffBoundary - (fEffBoundary - fEff3Boundary) * (x / k2);
    } else {
      effX = fEff3Boundary * (k - x) / (k - k2);
    }
    resX = fResBoundary + fResSlope*x/k;
    timeWalkX = fTimeWalkBoundary + fTimeWalkSlope*x/k;
    
    if(x < k && x > 0) {
      //   C:
      if(ix > 1 && dX < 0) {
	nActivatedPads++;
	nPlace[nActivatedPads-1] = nPlace[0] - 1;
	eff[nActivatedPads-1] = effX;
	res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resX * resX); // 10400=30^2+20^2+40^2+50^2+50^2+50^2 ns 
	timeWalk[nActivatedPads-1] = 0.001 * timeWalkX; // ns
	nTail[nActivatedPads-1] = 2;
	if (fTimeDelayFlag) {
	  //	  qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * x / 2.);
	  //	  qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * x / 2.);
	  qInduced[nActivatedPads-1] = TMath::Exp(-fPulseHeightSlope * x);
	  logOfqInd = gRandom->Gaus(-fPulseHeightSlope * x, fLogChargeSmearing);
	  timeDelay[nActivatedPads-1] = gRandom->Gaus(-fTimeDelaySlope * logOfqInd, fTimeSmearing);
	} else {
	  timeDelay[nActivatedPads-1] = 0.;
	}
	padId[nActivatedPads-1] = 3;

	//     D:
	if(z < k && z > 0) {
	  if( (iz == 1 && dZ > 0) || (iz == 2 && dZ < 0) ) {
	    nActivatedPads++;
	    nPlace[nActivatedPads-1] = nPlace[0] + (3 - 2 * iz) * AliTOFGeometry::NpadX() - 1;
	    eff[nActivatedPads-1] = effX * effZ;
	    (resZ<resX) ? res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resX * resX) : res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resZ * resZ); // 10400=30^2+20^2+40^2+50^2+50^2+50^2 ns
	    (timeWalkZ<timeWalkX) ? timeWalk[nActivatedPads-1] = 0.001 * timeWalkZ : timeWalk[nActivatedPads-1] = 0.001 * timeWalkX; // ns
	    
	    nTail[nActivatedPads-1] = 2;
	    if (fTimeDelayFlag) {
	      if (TMath::Abs(x) < TMath::Abs(z)) {
		//		qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * z / 2.);
		//		qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * z / 2.);
		qInduced[nActivatedPads-1] = TMath::Exp(-fPulseHeightSlope * z);
		logOfqInd = gRandom->Gaus(-fPulseHeightSlope * z, fLogChargeSmearing);
	      } else {
		//		qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * x / 2.);
		//		qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * x / 2.);
		qInduced[nActivatedPads-1] = TMath::Exp(-fPulseHeightSlope * x);
		logOfqInd = gRandom->Gaus(-fPulseHeightSlope * x, fLogChargeSmearing);
	      }
	      timeDelay[nActivatedPads-1] = gRandom->Gaus(-fTimeDelaySlope * logOfqInd, fTimeSmearing);
	    } else {
	      timeDelay[nActivatedPads-1] = 0.;
	    }
	    padId[nActivatedPads-1] = 4;
	  }
	}  // end D
      }  // end C
      
      //   E:
      if(ix < AliTOFGeometry::NpadX() && dX > 0) {
	nActivatedPads++;
	nPlace[nActivatedPads-1] = nPlace[0] + 1;
	eff[nActivatedPads-1] = effX;
	res[nActivatedPads-1] = 0.001 * (TMath::Sqrt(10400 + resX * resX)); // ns
	timeWalk[nActivatedPads-1] = 0.001 * timeWalkX; // ns
	nTail[nActivatedPads-1] = 2;
	if (fTimeDelayFlag) {
	  //	  qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * x / 2.);
	  //	  qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * x / 2.);
	  qInduced[nActivatedPads-1] = TMath::Exp(-fPulseHeightSlope * x);
	  logOfqInd = gRandom->Gaus(-fPulseHeightSlope * x, fLogChargeSmearing);
	  timeDelay[nActivatedPads-1] = gRandom->Gaus(-fTimeDelaySlope * logOfqInd, fTimeSmearing);
	} else {
	  timeDelay[nActivatedPads-1] = 0.;
	}
	padId[nActivatedPads-1] = 5;


	//     F:
	if(z < k && z > 0) {
	  if( (iz == 1 && dZ > 0) || (iz == 2 && dZ < 0) ) {
	    nActivatedPads++;
	    nPlace[nActivatedPads - 1] = nPlace[0] + (3 - 2 * iz) * AliTOFGeometry::NpadX() + 1;
	    eff[nActivatedPads - 1] = effX * effZ;
	    (resZ<resX) ? res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resX * resX) : res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resZ * resZ); // 10400=30^2+20^2+40^2+50^2+50^2+50^2 ns
	    (timeWalkZ<timeWalkX) ? timeWalk[nActivatedPads-1] = 0.001 * timeWalkZ : timeWalk[nActivatedPads-1] = 0.001*timeWalkX; // ns
	    nTail[nActivatedPads-1] = 2;
	    if (fTimeDelayFlag) {
	      if (TMath::Abs(x) < TMath::Abs(z)) {
		//		qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * z / 2.);
		//		qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * z / 2.);
		qInduced[nActivatedPads-1] = TMath::Exp(-fPulseHeightSlope * z);
		logOfqInd = gRandom->Gaus(-fPulseHeightSlope * z, fLogChargeSmearing);
	      } else {
		//		qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * x / 2.);
		//		qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * x / 2.);
		qInduced[nActivatedPads-1] = TMath::Exp(-fPulseHeightSlope * x);
		logOfqInd = gRandom->Gaus(-fPulseHeightSlope * x, fLogChargeSmearing);
	      }
	      timeDelay[nActivatedPads-1] = gRandom->Gaus(-fTimeDelaySlope * logOfqInd, fTimeSmearing);
	    } else {
	      timeDelay[nActivatedPads-1] = 0.;
	    }
	    padId[nActivatedPads-1] = 6;
	  }
	}  // end F
      }  // end E
    } // end if(x < k)


    for (Int_t iPad = 0; iPad < nActivatedPads; iPad++) {
      if (res[iPad] < fTimeResolution) res[iPad] = fTimeResolution;
      if(gRandom->Rndm() < eff[iPad]) {
	isFired[iPad] = kTRUE;
	nFiredPads++;
	if(fEdgeTails) {
	  if(nTail[iPad] == 0) {
	    tofTime[iPad] = gRandom->Gaus(geantTime + timeWalk[iPad] + timeDelay[iPad], res[iPad]);
	  } else {
	    ftail->SetParameters(res[iPad], 2. * res[iPad], kSigmaForTail[nTail[iPad]-1]);
	    Double_t timeAB = ftail->GetRandom();
	    tofTime[iPad] = geantTime + timeWalk[iPad] + timeDelay[iPad] + timeAB;
	  }
	} else {
	  tofTime[iPad] = gRandom->Gaus(geantTime + timeWalk[iPad] + timeDelay[iPad], res[iPad]);
	}
	if (fAverageTimeFlag) {
	  averageTime += tofTime[iPad] * qInduced[iPad];
	  weightsSum += qInduced[iPad];
	} else {
	  averageTime += tofTime[iPad];
	  weightsSum += 1.;
	}
      }
    }
    if (weightsSum!=0) averageTime /= weightsSum;
  } // end else (fEdgeEffect != 0)
}


/* new algorithm (to be checked)
//__________________________________________________________________
void AliTOFReconstructioner::BorderEffect(Float_t z0, Float_t x0, Float_t geantTime, Int_t& nActivatedPads, Int_t& nFiredPads, Bool_t* isFired, Int_t* nPlace, Float_t* qInduced, Float_t* tofTime, Float_t& averageTime)
{
  // Input:  z0, x0 - hit position in the strip system (0,0 - center of the strip), cm
  //         geantTime - time generated by Geant, ns
  // Output: nActivatedPads - the number of pads activated by the hit (1 || 2 || 4)
  //         nFiredPads - the number of pads fired (really activated) by the hit (nFiredPads <= nActivatedPads)
  //         qInduced[iPad]- charge induced on pad, arb. units
  //                         this array is initialized at zero by the caller
  //         tofAfterSimul[iPad] - time calculated with edge effect algorithm, ns
  //                                   this array is initialized at zero by the caller
  //         averageTime - time given by pad hited by the Geant track taking into account the times (weighted) given by the pads fired for edge effect also.
  //                       The weight is given by the qInduced[iPad]/qCenterPad
  //                                   this variable is initialized at zero by the caller
  //         nPlace[iPad] - the number of the pad place, iPad = 0, 1, 2, 3
  //                                   this variable is initialized at zero by the caller
  //
  // Description of used variables:
  //         eff[iPad] - efficiency of the pad
  //         res[iPad] - resolution of the pad, ns
  //         timeWalk[iPad] - time walk of the pad, ns
  //         timeDelay[iPad] - time delay for neighbouring pad to hited pad, ns
  //         PadId[iPad] - Pad Identifier
  //                    E | F    -->   PadId[iPad] = 5 | 6
  //                    A | B    -->   PadId[iPad] = 1 | 2
  //                    C | D    -->   PadId[iPad] = 3 | 4
  //         nTail[iPad] - the tail number, = 1 for tailA, = 2 for tailB
  //         qCenterPad - charge extimated for each pad, arb. units
  //         weightsSum - sum of weights extimated for each pad fired, arb. units
  
  const Float_t kSigmaForTail[2] = {AliTOFGeometry::SigmaForTail1(),AliTOFGeometry::SigmaForTail2()}; //for tail                                                   
  Int_t iz = 0, ix = 0;
  Float_t dX = 0., dZ = 0., x = 0., z = 0.;
  Float_t h = fHparameter, h2 = fH2parameter, k = fKparameter, k2 = fK2parameter;
  Float_t effX = 0., effZ = 0., resX = 0., resZ = 0., timeWalkX = 0., timeWalkZ = 0.;
  Float_t logOfqInd = 0.;
  Float_t weightsSum = 0.;
  Int_t nTail[4]  = {0,0,0,0};
  Int_t padId[4]  = {0,0,0,0};
  Float_t eff[4]  = {0.,0.,0.,0.};
  Float_t res[4]  = {0.,0.,0.,0.};
  Float_t qCenterPad = fMinimumCharge * fMinimumCharge;
  Float_t timeWalk[4]  = {0.,0.,0.,0.};
  Float_t timeDelay[4] = {0.,0.,0.,0.};
  
  nActivatedPads = 0;
  nFiredPads = 0;
  
  (z0 <= 0) ? iz = 0 : iz = 1;
  dZ = z0 + (0.5 * AliTOFGeometry::NpadZ() - iz - 0.5) * AliTOFGeometry::ZPad(); // hit position in the pad frame, (0,0) - center of the pad
  z = 0.5 * AliTOFGeometry::ZPad() - TMath::Abs(dZ);                               // variable for eff., res. and timeWalk. functions
  iz++;                                                                              // z row: 1, ..., AliTOFGeometry::NpadZ() = 2
  ix = (Int_t)((x0 + 0.5 * AliTOFGeometry::NpadX() * AliTOFGeometry::XPad()) / AliTOFGeometry::XPad());
  dX = x0 + (0.5 * AliTOFGeometry::NpadX() - ix - 0.5) * AliTOFGeometry::XPad(); // hit position in the pad frame, (0,0) - center of the pad
  x = 0.5 * AliTOFGeometry::XPad() - TMath::Abs(dX);                               // variable for eff., res. and timeWalk. functions;
  ix++;                                                                              // x row: 1, ..., AliTOFGeometry::NpadX() = 48
  
  ////// Pad A:
  nActivatedPads++;
  nPlace[nActivatedPads-1] = (iz - 1) * AliTOFGeometry::NpadX() + ix;
  qInduced[nActivatedPads-1] = qCenterPad;
  padId[nActivatedPads-1] = 1;
  
  if (fEdgeEffect == 0) {
    eff[nActivatedPads-1] = fEffCenter;
    if (gRandom->Rndm() < eff[nActivatedPads-1]) {
      nFiredPads = 1;
      res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + fResCenter * fResCenter); // 10400=30^2+20^2+40^2+50^2+50^2+50^2  ns;
      isFired[nActivatedPads-1] = kTRUE;
      tofTime[nActivatedPads-1] = gRandom->Gaus(geantTime + fTimeWalkCenter, res[0]);
      averageTime = tofTime[nActivatedPads-1];
    }
  } else {
     
    if(z < h) {
      if(z < h2) {
	effZ = fEffBoundary + (fEff2Boundary - fEffBoundary) * z / h2;
      } else {
	effZ = fEff2Boundary + (fEffCenter - fEff2Boundary) * (z - h2) / (h - h2);
      }
      resZ = fResBoundary + (fResCenter - fResBoundary) * z / h;
      timeWalkZ = fTimeWalkBoundary + (fTimeWalkCenter - fTimeWalkBoundary) * z / h;
      nTail[nActivatedPads-1] = 1;
    } else {
      effZ = fEffCenter;
      resZ = fResCenter;
      timeWalkZ = fTimeWalkCenter;
    }
    
    if(x < h) {
      if(x < h2) {
	effX = fEffBoundary + (fEff2Boundary - fEffBoundary) * x / h2;
      } else {
	effX = fEff2Boundary + (fEffCenter - fEff2Boundary) * (x - h2) / (h - h2);
      }
      resX = fResBoundary + (fResCenter - fResBoundary) * x / h;
      timeWalkX = fTimeWalkBoundary + (fTimeWalkCenter - fTimeWalkBoundary) * x / h;
      nTail[nActivatedPads-1] = 1;
    } else {
      effX = fEffCenter;
      resX = fResCenter;
      timeWalkX = fTimeWalkCenter;
    }
    
    (effZ<effX) ? eff[nActivatedPads-1] = effZ : eff[nActivatedPads-1] = effX;
    (resZ<resX) ? res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resX * resX) : res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resZ * resZ); // 10400=30^2+20^2+40^2+50^2+50^2+50^2  ns
    (timeWalkZ<timeWalkX) ? timeWalk[nActivatedPads-1] = 0.001 *  timeWalkZ : timeWalk[nActivatedPads-1] = 0.001 * timeWalkX; // ns


    ////// Pad B:
    if(z < k2) {
      effZ = fEffBoundary - (fEffBoundary - fEff3Boundary) * (z / k2);
    } else {
      effZ = fEff3Boundary * (k - z) / (k - k2);
    }
    resZ = fResBoundary + fResSlope * z / k;
    timeWalkZ = fTimeWalkBoundary + fTimeWalkSlope * z / k;
    
    if(z < k && z > 0) {
      if( (iz == 1 && dZ > 0) || (iz == 2 && dZ < 0) ) {
	nActivatedPads++;
	nPlace[nActivatedPads-1] = nPlace[0] + (3 - 2 * iz) * AliTOFGeometry::NpadX();
	eff[nActivatedPads-1] = effZ;
	res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resZ * resZ); // 10400=30^2+20^2+40^2+50^2+50^2+50^2 ns 
	timeWalk[nActivatedPads-1] = 0.001 * timeWalkZ; // ns
	nTail[nActivatedPads-1] = 2;
	if (fTimeDelayFlag) {
	  qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * z / 2.);
	  qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * z / 2.);
	  logOfqInd = gRandom->Gaus(-fPulseHeightSlope * z, fLogChargeSmearing);
	  timeDelay[nActivatedPads-1] = gRandom->Gaus(-fTimeDelaySlope * logOfqInd, fTimeSmearing);
	} else {
	  timeDelay[nActivatedPads-1] = 0.;
	}
	padId[nActivatedPads-1] = 2;
      }
    }

    
    ////// Pad C, D, E, F:
    if(x < k2) {
      effX = fEffBoundary - (fEffBoundary - fEff3Boundary) * (x / k2);
    } else {
      effX = fEff3Boundary * (k - x) / (k - k2);
    }
    resX = fResBoundary + fResSlope*x/k;
    timeWalkX = fTimeWalkBoundary + fTimeWalkSlope*x/k;
    
    if(x < k && x > 0) {
      //   C:
      if(ix > 1 && dX < 0) {
	nActivatedPads++;
	nPlace[nActivatedPads-1] = nPlace[0] - 1;
	eff[nActivatedPads-1] = effX;
	res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resX * resX); // 10400=30^2+20^2+40^2+50^2+50^2+50^2 ns 
	timeWalk[nActivatedPads-1] = 0.001 * timeWalkX; // ns
	nTail[nActivatedPads-1] = 2;
	if (fTimeDelayFlag) {
	  qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * x / 2.);
	  qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * x / 2.);
	  logOfqInd = gRandom->Gaus(-fPulseHeightSlope * x, fLogChargeSmearing);
	  timeDelay[nActivatedPads-1] = gRandom->Gaus(-fTimeDelaySlope * logOfqInd, fTimeSmearing);
	} else {
	  timeDelay[nActivatedPads-1] = 0.;
	}
	padId[nActivatedPads-1] = 3;

	//     D:
	if(z < k && z > 0) {
	  if( (iz == 1 && dZ > 0) || (iz == 2 && dZ < 0) ) {
	    nActivatedPads++;
	    nPlace[nActivatedPads-1] = nPlace[0] + (3 - 2 * iz) * AliTOFGeometry::NpadX() - 1;
	    eff[nActivatedPads-1] = effX * effZ;
	    (resZ<resX) ? res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resX * resX) : res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resZ * resZ); // 10400=30^2+20^2+40^2+50^2+50^2+50^2 ns
	    (timeWalkZ<timeWalkX) ? timeWalk[nActivatedPads-1] = 0.001 * timeWalkZ : timeWalk[nActivatedPads-1] = 0.001 * timeWalkX; // ns
	    
	    nTail[nActivatedPads-1] = 2;
	    if (fTimeDelayFlag) {
	      if (TMath::Abs(x) < TMath::Abs(z)) {
		qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * z / 2.);
		qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * z / 2.);
		logOfqInd = gRandom->Gaus(-fPulseHeightSlope * z, fLogChargeSmearing);
	      } else {
		qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * x / 2.);
		qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * x / 2.);
		logOfqInd = gRandom->Gaus(-fPulseHeightSlope * x, fLogChargeSmearing);
	      }
	      timeDelay[nActivatedPads-1] = gRandom->Gaus(-fTimeDelaySlope * logOfqInd, fTimeSmearing);
	    } else {
	      timeDelay[nActivatedPads-1] = 0.;
	    }
	    padId[nActivatedPads-1] = 4;
	  }
	}  // end D
      }  // end C
      
      //   E:
      if(ix < AliTOFGeometry::NpadX() && dX > 0) {
	nActivatedPads++;
	nPlace[nActivatedPads-1] = nPlace[0] + 1;
	eff[nActivatedPads-1] = effX;
	res[nActivatedPads-1] = 0.001 * (TMath::Sqrt(10400 + resX * resX)); // ns
	timeWalk[nActivatedPads-1] = 0.001 * timeWalkX; // ns
	nTail[nActivatedPads-1] = 2;
	if (fTimeDelayFlag) {
	  qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * x / 2.);
	  qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * x / 2.);
	  logOfqInd = gRandom->Gaus(-fPulseHeightSlope * x, fLogChargeSmearing);
	  timeDelay[nActivatedPads-1] = gRandom->Gaus(-fTimeDelaySlope * logOfqInd, fTimeSmearing);
	} else {
	  timeDelay[nActivatedPads-1] = 0.;
	}
	padId[nActivatedPads-1] = 5;


	//     F:
	if(z < k && z > 0) {
	  if( (iz == 1 && dZ > 0) || (iz == 2 && dZ < 0) ) {
	    nActivatedPads++;
	    nPlace[nActivatedPads - 1] = nPlace[0] + (3 - 2 * iz) * AliTOFGeometry::NpadX() + 1;
	    eff[nActivatedPads - 1] = effX * effZ;
	    (resZ<resX) ? res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resX * resX) : res[nActivatedPads-1] = 0.001 * TMath::Sqrt(10400 + resZ * resZ); // 10400=30^2+20^2+40^2+50^2+50^2+50^2 ns
	    (timeWalkZ<timeWalkX) ? timeWalk[nActivatedPads-1] = 0.001 * timeWalkZ : timeWalk[nActivatedPads-1] = 0.001*timeWalkX; // ns
	    nTail[nActivatedPads-1] = 2;
	    if (fTimeDelayFlag) {
	      if (TMath::Abs(x) < TMath::Abs(z)) {
		qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * z / 2.);
		qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * z / 2.);
		logOfqInd = gRandom->Gaus(-fPulseHeightSlope * z, fLogChargeSmearing);
	      } else {
		qInduced[0] = fMinimumCharge * TMath::Exp(fPulseHeightSlope * x / 2.);
		qInduced[nActivatedPads-1] = fMinimumCharge * TMath::Exp(-fPulseHeightSlope * x / 2.);
		logOfqInd = gRandom->Gaus(-fPulseHeightSlope * x, fLogChargeSmearing);
	      }
	      timeDelay[nActivatedPads-1] = gRandom->Gaus(-fTimeDelaySlope * logOfqInd, fTimeSmearing);
	    } else {
	      timeDelay[nActivatedPads-1] = 0.;
	    }
	    padId[nActivatedPads-1] = 6;
	  }
	}  // end F
      }  // end E
    } // end if(x < k)


    for (Int_t iPad = 0; iPad < nActivatedPads; iPad++) {
      if (res[iPad] < fTimeResolution) res[iPad] = fTimeResolution;
      if(gRandom->Rndm() < eff[iPad]) {
	isFired[iPad] = kTRUE;
	nFiredPads++;
	if(fEdgeTails) {
	  if(nTail[iPad] == 0) {
	    tofTime[iPad] = gRandom->Gaus(geantTime + timeWalk[iPad] + timeDelay[iPad], res[iPad]);
	  } else {
	    ftail->SetParameters(res[iPad], 2. * res[iPad], kSigmaForTail[nTail[iPad]-1]);
	    Double_t timeAB = ftail->GetRandom();
	    tofTime[iPad] = geantTime + timeWalk[iPad] + timeDelay[iPad] + timeAB;
	  }
	} else {
	  tofTime[iPad] = gRandom->Gaus(geantTime + timeWalk[iPad] + timeDelay[iPad], res[iPad]);
	}
	if (fAverageTimeFlag) {
	  averageTime += tofTime[iPad] * qInduced[iPad];
	  weightsSum += qInduced[iPad];
	} else {
	  averageTime += tofTime[iPad];
	  weightsSum += 1.;
	}
      }
    }
    if (weightsSum!=0) averageTime /= weightsSum;

  } // end else (fEdgeEffect != 0)
  
  //cout << "timedelay " << timeDelay[0] << endl;
  //cout << "timedelay " << timeDelay[1] << endl;
  //cout << "timedelay " << timeDelay[2] << endl;
  //cout << "timedelay " << timeDelay[3] << endl;
  
}
*/


//__________________________________________________________________
Int_t AliTOFReconstructioner::PDGtoGeantCode(Int_t pdgcode) 
{
  //
  // Gives the GEANT code from KF code of LUND JETSET
  //
  Int_t geantCode=0; // default value
  switch (pdgcode) {
  case 22:
    geantCode=1;            // GAMMA
    break ;
  case -11:
    geantCode=2;            // E+
    break ;
  case 11:
    geantCode=3;            // E-
    break ;
  case 12:
    geantCode=4;            // NUE
    break ;
  case 14:
    geantCode=4;            // NUMU
    break ;
  case -13:
    geantCode=5;            // MU+
    break ;
  case 13:
    geantCode=6;            // MU-
    break ;
  case 111:
    geantCode=7;            // PI0
    break ;
  case 211:
    geantCode=8;            // PI+
    break ;
  case -211:
    geantCode=9;            // PI-
    break ;
  case 130:
    geantCode=10;           // K_L0
    break ;
  case 321:
    geantCode=11;           // K+
    break ;
  case -321:
    geantCode=12;           // K-
    break ;
  case 2112:
    geantCode=13;           // N0
    break ;
  case 2212:
    geantCode=14;           // P+
    break ;
  case -2212:
    geantCode=15;           // P~-
    break ;
  case 310:
    geantCode=16;           // K_S0
    break ;
  case 221:
    geantCode=17;           // ETA
    break ;
  case 3122:
    geantCode=18;           // LAMBDA0
    break ;
  case 3222:
    geantCode=19;           // SIGMA+
    break ;
  case 3212:
    geantCode=20;           // SIGMA0
    break ;
  case 3112:
    geantCode=21;           // SIGMA-
    break ;
  case 3322:
    geantCode=22;           // XI0
    break ;
  case 3312:
    geantCode=23;           // XI-
    break ;
  case 3334:
    geantCode=24;           // OMEGA-
    break ;
  case -2112:
    geantCode=25;           // N~0
    break ;
  case -3122:
    geantCode=26;           // LAMBDA~0
    break ;
  case -3112:
    geantCode=27;           // SIGMA~+
    break ;
  case -3212:
    geantCode=28;           // SIGMA~0
    break ;
  case -3222:
    geantCode=29;           // SIGMA~-
    break ;
  case -3322:
    geantCode=30;           // XI~0
    break ;
  case -3312:
    geantCode=31;           // XI~+
    break ;
  case -3334:
    geantCode=32;           // OMEGA~+
    break ;
  case 223:
    geantCode=33;           // OMEGA(782)
    break ;
  case 333:
    geantCode=34;           // PHI(1020)
    break ;
  case 411:
    geantCode=35;           // D+
    break ;
  case -411:
    geantCode=36;           // D-
    break ;
  case 421:
    geantCode=37;           // D0
    break ;
  case -421:
    geantCode=38;           // D~0
    break ;
  case 431:
    geantCode=39;           // D_S+
    break ;
  case -431:
    geantCode=40;           // D_S~-
    break ;
  case 4122:
    geantCode=41;           // LAMBDA_C+
    break ;
  case 213:
    geantCode=42;           // RHP(770)+
    break ;
  case -213:
    geantCode=43;           // RHO(770)-
    break ;
  case 113:
    geantCode=44;           // RHO(770)0
    break ;
  default:
    geantCode=45;
    break;
  }

  return geantCode;
}

//__________________________________________________________________
Bool_t AliTOFReconstructioner::operator==( AliTOFReconstructioner const & tofrec)const
{
  // Equal operator.
  // Reconstructioners are equal if their parameters are equal

  // split the member variables in analogous categories
  
  // time resolution and edge effect parameters
  Bool_t dummy0=(fTimeResolution==tofrec.fTimeResolution)&&
	  (fpadefficiency==tofrec.fpadefficiency)&&
	  (fEdgeEffect==tofrec.fEdgeEffect)&&
	  (fEdgeTails==tofrec.fEdgeTails)&&
	  (fHparameter==tofrec.fHparameter)&&
	  (fH2parameter==tofrec.fH2parameter)&&
	  (fKparameter==tofrec.fKparameter)&&
	  (fK2parameter==tofrec.fK2parameter);
  
  // pad efficiency parameters
  Bool_t dummy1=(fEffCenter==tofrec.fEffCenter)&&
	 (fEffBoundary==tofrec.fEffBoundary)&&
	 (fEff2Boundary==tofrec.fEff2Boundary)&&
	 (fEff3Boundary==tofrec.fEff3Boundary)&&
	 (fResCenter==tofrec.fResCenter)&&
	 (fResBoundary==tofrec.fResBoundary)&&
	 (fResSlope==tofrec.fResSlope);

  // time walk parameters
  Bool_t dummy2=(fTimeWalkCenter==tofrec.fTimeWalkCenter)&&
	  (fTimeWalkBoundary==tofrec.fTimeWalkBoundary)&&
	  (fTimeWalkSlope==tofrec.fTimeWalkSlope)&&
	  (fTimeDelayFlag==tofrec.fTimeDelayFlag)&&
	  (fPulseHeightSlope==tofrec.fPulseHeightSlope)&&
	  (fTimeDelaySlope==tofrec.fTimeDelaySlope);

  // ADC-TDC correlation parameters
  Bool_t dummy3=(fMinimumCharge==tofrec.fMinimumCharge)&&
	  (fChargeSmearing==tofrec.fChargeSmearing )&&
	  (fLogChargeSmearing==tofrec.fLogChargeSmearing )&&
	  (fTimeSmearing==tofrec.fTimeSmearing )&&
	  (fAverageTimeFlag==tofrec.fAverageTimeFlag)&&
	  (fChargeFactorForMatching==tofrec.fChargeFactorForMatching)&&
	  (fMatchingStyle==tofrec.fMatchingStyle);
  
  Bool_t dummy4=(fTrackingEfficiency==tofrec.fTrackingEfficiency)&&
	  (fSigmavsp==tofrec.fSigmavsp)&&
	  (fSigmaZ==tofrec.fSigmaZ)&&
	  (fSigmarphi==tofrec.fSigmarphi)&&
	  (fSigmap==tofrec.fSigmap)&&
	  (fSigmaPhi==tofrec.fSigmaPhi)&&
	  (fSigmaTheta==tofrec.fSigmaTheta)&&
	  (fNoise==tofrec.fNoise)&&
	  (fNoiseSlope==tofrec.fNoiseSlope)&&
	  (fField==tofrec.fField)&&
	  (fRadLenTPC==tofrec.fRadLenTPC)&&
	  (fCorrectionTRD==tofrec.fCorrectionTRD)&&
	  (fLastTPCRow==tofrec.fLastTPCRow)&&
	  (fRadiusvtxBound==tofrec.fRadiusvtxBound)&&
	  (fMaxTestTracks==tofrec.fMaxTestTracks)&&
	  (fStep==tofrec.fStep)&&
	  (fMaxPixels==tofrec.fMaxPixels)&&
	  (fMaxAllTracks==tofrec.fMaxAllTracks)&&
	  (fMaxTracks==tofrec.fMaxTracks)&&
	  (fMaxTOFHits==tofrec.fMaxTOFHits)&&
	  (fPBound==tofrec.fPBound);

  if( dummy0 && dummy1 && dummy2 && dummy3 && dummy4)
    return kTRUE ;
  else
    return kFALSE ;

}
//____________________________________________________________________________ 
void AliTOFReconstructioner::UseHitsFrom(const char * filename)
{
  SetTitle(filename) ; 
}

//____________________________________________________________________________ 
void AliTOFReconstructioner::InitArray(Float_t array[], Int_t nlocations)
{
  //
  // Initialize the array of Float_t
  // 
  for (Int_t i = 0; i < nlocations; i++) {
    array[i]=0.;
  }				// end loop

}

//____________________________________________________________________________ 
void AliTOFReconstructioner::InitArray(Int_t array[], Int_t nlocations)
{
  //
  // Initialize the array of Int_t
  // 
  for (Int_t i = 0; i < nlocations; i++) {
    array[i]=0;
  }				// end loop

}


//____________________________________________________________________________
void AliTOFReconstructioner::ReadTOFHits(Int_t ntracks, TTree* treehits, 
		TClonesArray* tofhits, Int_t ***MapPixels, Int_t* kTOFhitFirst, 
		AliTOFPad* pixelArray , Int_t* iTOFpixel, Float_t* toftime, 
		AliTOFRecHit* hitArray, Int_t& isHitOnFiredPad, Int_t& ipixel)
{
  //
  // Read TOF hits for the current event and fill arrays
  // 
  // Start loop on primary tracks in the hits containers
  //
  // Noise meaning in ReadTOFHits: we use the word 'noise' in the  
  // following cases
  // - signals produced by secondary particles
  // - signals produced by the next hits (out of the first) of a given track
  //   (both primary and secondary)
  // - signals produced by edge effect


  TParticle *particle;
  Int_t nHitOutofTofVolumes; // number of hits out of TOF GEANT volumes (it happens in very
                             // few cases)
  Int_t * npixel = new Int_t[AliTOFGeometry::MaxTOFTree()]; // array used by TOFRecon for check on TOF geometry
  Int_t npions=0;    // number of pions for the current event
  Int_t nkaons=0;    // number of kaons for the current event
  Int_t nprotons=0;  // number of protons for the current event
  Int_t nelectrons=0;// number of electrons for the current event
  Int_t nmuons=0;    // number of muons for the current event
  Float_t tofpos[3];     // TOF hit position and GEANT time
  Float_t zPad,xPad; 
  Int_t nbytes = 0;
  Int_t ipart, nhits=0, nHitsFromPrimaries=0;
  Int_t ntotalTOFhits=0; // total number of TOF hits for the current event
  Int_t ipartLast=-1;    // last track identifier
  Int_t iFirstHit;       // flag to check if the current hit is the first hit on TOF for the 
                         // current track
  Int_t iNoiseHit=0;     // flag used to tag noise hits (the noise meaning is reported in the
                         // header of the ReadTOFHits method)
  Int_t nhitWithoutNoise;// number of hits not due to noise
  Int_t inoise=0,inoise2=0;
  Int_t nMultipleSignOnSamePad=0; // number of cases where a pad is fired more than one time
  Int_t nPixEdge=0;      // additional pads fired due to edge effect in ReadTOFHits (local var)
  // array used for counting different types of primary particles
  Int_t particleTypeGEANT[50]={0,4,4,0,5,5,0,3,3,0,
		   2,2,0,1,1,0,0,0,0,0,
		   0,0,0,0,0,0,0,0,0,0,
		   0,0,0,0,0,0,0,0,0,0,
		   0,0,0,0,0,0,0,0,0,0};
  Int_t particleType,particleInTOFtype[6][3];
  for (Int_t i=0;i<6;i++) {
    for (Int_t j=0;j<3;j++) {
      particleInTOFtype[i][j]=0;
    }
  }

  // speed-up the code
  treehits->SetBranchStatus("*",0); // switch off all branches
  treehits->SetBranchStatus("TOF*",1); // switch on only TOF

  for (Int_t track=0; track<ntracks;track++) { // starting loop on primary tracks for the current event

    gAlice->ResetHits();
    nbytes += treehits->GetEvent(track);
    nhits = tofhits->GetEntriesFast();

    ntotalTOFhits+=nhits;      

    // Start loop on hits connected to the current primary tracked particle
    // (including hits produced by secondary particles generaterd by the
    // current ptimary tracked particle)
    for (Int_t hit=0;hit<nhits;hit++) {
      AliTOFhit* tofHit   = (AliTOFhit*)tofhits->UncheckedAt(hit);
      ipart    = tofHit->GetTrack();
      if(ipart>=fMaxAllTracks) break;
      Float_t geantTime= tofHit->GetTof(); // it is given in [s]
      particle = (TParticle*)gAlice->GetMCApp()->Particle(ipart);

      Int_t pdgCode=particle->GetPdgCode();
      // Only high momentum tracks (see fPBound value)
      // momentum components at vertex
      Float_t pxvtx = particle->Px();
      Float_t pyvtx = particle->Py();
      Float_t pzvtx = particle->Pz();
      Float_t pvtx = TMath::Sqrt(pxvtx*pxvtx+pyvtx*pyvtx+pzvtx*pzvtx);
      if(pvtx>fPBound) {
	
	if(particle->GetFirstMother() < 0) nHitsFromPrimaries++; // count primaries

	// x and y coordinates of the particle production vertex
	Float_t vx = particle->Vx();
	Float_t vy = particle->Vy();
	Float_t vr = TMath::Sqrt(vx*vx+vy*vy); // cylindrical radius of the particle production vertex
	
	Float_t x = tofHit->X(); tofpos[0]=x;
	Float_t y = tofHit->Y(); tofpos[1]=y;
	Float_t z = tofHit->Z(); tofpos[2]=z;
	/* var used for QA
	Float_t tofradius = TMath::Sqrt(x*x+y*y); // radius cilindrical coordinate of the TOF hit
	*/
	// momentum components (cosine) when striking the TOF
	Float_t pxtof = tofHit->GetPx();
	Float_t pytof = tofHit->GetPy();
	Float_t pztof = tofHit->GetPz();
	// scalar product indicating the direction of the particle when striking the TOF 
	/* var used for QA
	// (>0 for outgoing particles)
	Float_t isGoingOut = (x*pxtof+y*pytof+z*pztof)/TMath::Sqrt(x*x+y*y+z*z);
	*/
	Float_t momtof = tofHit->GetMom();
	// now momentum components when striking the TOF
	pxtof *= momtof;
	pytof *= momtof;
	pztof *= momtof;
	particleType=particleTypeGEANT[PDGtoGeantCode(pdgCode)-1];
	if(particleType) {
	  particleInTOFtype[5][2]++;
	  particleInTOFtype[particleType-1][2]++;
	}
	iFirstHit=0;
	//  without noise hits
	
	if(ipart!=ipartLast) {
	  iFirstHit=1;
	  toftime[ipart]=geantTime;    //time [s]
	  // tofMom[ipart]=momtof;
	  ipartLast=ipart;
	  if(particle->GetFirstMother() < 0) {
	    Int_t abspdgCode=TMath::Abs(pdgCode);
	    switch (abspdgCode) {
	    case 211:
	      npions++;
	      break ;
	    case 321:
	      nkaons++;
	      break ;
	    case 2212:
	      nprotons++;
	      break ;
	    case 11:
	      nelectrons++;
	      break ;
	    case 13:
	      nmuons++;
	      break ;
	    }
	  }
	  if(vr>fRadiusvtxBound) {
	    if(particleType) { 
	      particleInTOFtype[5][1]++;
	      particleInTOFtype[particleType-1][1]++;
	    }
	    inoise++;
	    inoise2++;
	  } else {
	    if(particleType) {
	      particleInTOFtype[5][0]++;
	      particleInTOFtype[particleType-1][0]++;
	    }
	  }
	} else {
	  inoise++;
	  if(particleType) {
	    particleInTOFtype[5][1]++;
	    particleInTOFtype[particleType-1][1]++;
	  }
	} //end if(ipart!=ipartLast)

	IsInsideThePad(gMC,x,y,z,npixel,zPad,xPad);

	Int_t sec  = tofHit->GetSector();
	Int_t pla  = tofHit->GetPlate();
	Int_t str  = tofHit->GetStrip();
	if(sec!=npixel[0] || pla!=npixel[1] || str!=npixel[2]){// check on volume 
	  cout << "sector" << sec << " npixel[0] " << npixel[0] << endl;
	  cout << "plate " << pla << " npixel[1] " << npixel[1] << endl;
	  cout << "strip " << str << " npixel[2] " << npixel[2] << endl;
	} // close check on volume
	
	Int_t padz = tofHit->GetPadz();
	Int_t padx = tofHit->GetPadx();
	Float_t Zpad = tofHit->GetDz();
	Float_t Xpad = tofHit->GetDx();
	
	
	if (npixel[4]==0){
	  IsInsideThePad(gMC,x,y,z,npixel,zPad,xPad);
	  if (npixel[4]==0){	      
	    nHitOutofTofVolumes++;
	  }
	} else {
	  Float_t zStrip=AliTOFGeometry::ZPad()*(padz-0.5-0.5*AliTOFGeometry::NpadZ())+Zpad; 
	  if(padz!=npixel[3]) printf("            : Zpad=%f, padz=%i, npixel[3]=%i, zStrip=%f\n",Zpad,padz,npixel[3],zStrip);
	  Float_t xStrip=AliTOFGeometry::XPad()*(padx-0.5-0.5*AliTOFGeometry::NpadX())+Xpad;
	  
	  Int_t nPlace[4]={0,0,0,0};
	  nPlace[0]=(padz-1)*AliTOFGeometry::NpadX()+padx;
	  
	  Int_t   nActivatedPads=0;
	  Int_t   nFiredPads=0;
	  Bool_t  isFired[4]={kFALSE,kFALSE,kFALSE,kFALSE};
	  Float_t tofAfterSimul[4]={0.,0.,0.,0.};
	  Float_t qInduced[4]={0.,0.,0.,0.};
	  Float_t averageTime=0.;    


	  BorderEffect(zStrip,xStrip,geantTime*1.0e+09,nActivatedPads,nFiredPads,isFired,nPlace,qInduced,tofAfterSimul,averageTime); // simulate edge effect


	  if(nFiredPads) {
	    for(Int_t indexOfPad=0; indexOfPad<nActivatedPads; indexOfPad++) {
	      if(isFired[indexOfPad]){// the pad has fired
		if(indexOfPad==0) {// the hit belongs to a fired pad
		  isHitOnFiredPad++;
		  hitArray[isHitOnFiredPad-1].SetHit(ipart,pdgCode,tofpos,momtof,vr,iFirstHit);
		  iNoiseHit=0;

		  if(vr>fRadiusvtxBound || iFirstHit==0) iNoiseHit=1;
		  
		  hitArray[isHitOnFiredPad-1].SetNoise(iNoiseHit);
		  if(iFirstHit) kTOFhitFirst[ipart]=isHitOnFiredPad;	  

		}// close - the hit belongs to a fired pad
		
		Int_t iMapFirstIndex=AliTOFGeometry::NSectors()*(npixel[1]-1)+npixel[0]-1;
		Int_t iMapValue=MapPixels[iMapFirstIndex][npixel[2]-1][nPlace[indexOfPad]-1];
	
		if(iMapValue==0) {
		  ipixel++;
		  if(indexOfPad) {
		    iNoiseHit=1;
		    nPixEdge++;
		  } else { 
		    iTOFpixel[ipart]=ipixel;
		  }
		  
		  if(ipixel>fMaxPixels){ // check on the total number of activated pads
		    cout << "ipixel=" << ipixel << " > fMaxPixels=" << fMaxPixels << endl;
		    return;
		  } // close check on the number of activated pads
		  
		  MapPixels[iMapFirstIndex][npixel[2]-1][nPlace[indexOfPad]-1]=ipixel;
		  pixelArray[ipixel-1].SetGeom(npixel[0],npixel[1],npixel[2],nPlace[indexOfPad]);
		  pixelArray[ipixel-1].SetTrack(ipart);
		  if(iNoiseHit) {
		    pixelArray[ipixel-1].AddState(1);
		  } else {
		    if(tofAfterSimul[indexOfPad]<0) cout << "Time of Flight after detector simulation is negative" << endl;
		    pixelArray[ipixel-1].AddState(10);
		  }
		  
		  pixelArray[ipixel-1].SetTofChargeHit(tofAfterSimul[indexOfPad],qInduced[indexOfPad],geantTime*1.0e+09,isHitOnFiredPad);
		} else { //else if(iMapValue==0)
		  if(indexOfPad==0) iTOFpixel[ipart]=iMapValue;
		  nMultipleSignOnSamePad++;
		  
		  if(tofAfterSimul[indexOfPad] < pixelArray[iMapValue-1].GetRealTime() ) {
		    pixelArray[iMapValue-1].SetTrack(ipart);
		    //                   if(indexOfPad==0) pixelArray[iMapValue-1].SetTrack(ipart);
		    if(indexOfPad) iNoiseHit=1;
		    if(iNoiseHit) {
		      pixelArray[iMapValue-1].AddState(1);
		    } else {
		      pixelArray[iMapValue-1].AddState(10);
		    }
		    pixelArray[iMapValue-1].SetRealTime(tofAfterSimul[indexOfPad]);
		    pixelArray[iMapValue-1].SetGeantTime(geantTime*1.0e+09);
		    pixelArray[iMapValue-1].SetHit(isHitOnFiredPad);
		  } // close if(tofAfterSimul[indexOfPad] < pixelArray[iMapValue-1].GetTime() )
		}  //end of Pixel filling
	      }  // close if(isFired[indexOfPad])	
	    }  //end loop on activated pads indexOfPad
	  } // close if(nFiredPads)		  
	}  //end of hit with npixel[3]!=0
      }  //high momentum tracks
    }  //end on TOF hits
  }  //end on primary tracks
  
  
  if(fdbg) {
    cout << ntotalTOFhits << " - total number of TOF hits   " << nHitsFromPrimaries << " -  primary     " <<  endl; 
    cout << inoise << " - noise hits, " << inoise2<< " - first crossing of a track with Rvtx>" << fRadiusvtxBound << endl;
    //   cout << inoise << " - noise hits (" << 100.*inoise/ihit << " %), " << inoise2 
    //<< " - first crossing of a track with Rvtx>" << RVTXBOUND << endl;
    nhitWithoutNoise=isHitOnFiredPad;
    
    cout << ipixel << " fired pixels (" << nMultipleSignOnSamePad << " multiple fired pads, " << endl;
    //j << " fired by noise, " << j1 << " noise+track)" <<  endl;
    printf(" %i additional pads are fired due to edge effect\n",nPixEdge);
    cout << npions <<   "   primary pions     reached TOF" << endl;
    cout << nkaons <<   "   primary kaons     reached TOF" << endl;
    cout << nprotons << "   primary protons   reached TOF" << endl;
    cout << nelectrons<<"   primary electrons reached TOF" << endl;
    cout << nmuons <<   "   primary muons     reached TOF" << endl;
    cout << "number of TOF hits for different species: 1-p, 2-K, 3-pi, 4-e, 5-mu, 6-all" << endl;
    cout << "   first number - track hits, second - noise ones, third - all" << endl;
    for (Int_t i=0;i<6;i++) cout << i+1 << "  " << particleInTOFtype[i][0] << "  " << particleInTOFtype[i][1] << "  " << particleInTOFtype[i][2] << endl; 

    Int_t primaryReachedTOF[6];
    primaryReachedTOF[0]=npions;
    primaryReachedTOF[1]=nkaons;
    primaryReachedTOF[2]=nprotons;
    primaryReachedTOF[3]=nelectrons;
    primaryReachedTOF[4]=nmuons;
    primaryReachedTOF[5]=npions+nkaons+nprotons+nelectrons+nmuons;
    
    cout << " Reading TOF hits done" << endl;
  }

  delete [] npixel;
}

//____________________________________________________________________________
void AliTOFReconstructioner::AddNoiseFromOuter(Option_t *option, Int_t ***MapPixels, AliTOFPad* pixelArray , AliTOFRecHit* hitArray, Int_t& isHitOnFiredPad, Int_t& ipixel)
{
  //
  // Add noise hits from outer regions (forward and backward) according
  // to parameterized fZNoise distribution (to be used with events 
  // generated in the barrel region)

  Float_t * zLen = new Float_t[AliTOFGeometry::NPlates()+1];
  Float_t * zStrips = new Float_t[AliTOFGeometry::NPlates()];
  zStrips[0]=(Float_t) (AliTOFGeometry::NStripC());
  zStrips[1]=(Float_t) (AliTOFGeometry::NStripB());
  zStrips[2]=(Float_t) (AliTOFGeometry::NStripA());
  zStrips[3]=(Float_t) (AliTOFGeometry::NStripB());
  zStrips[4]=(Float_t) (AliTOFGeometry::NStripC());

  zLen[5]=AliTOFGeometry::ZlenA()*0.5+AliTOFGeometry::ZlenB()+AliTOFGeometry::ZlenC();
  zLen[4]=zLen[5]-AliTOFGeometry::ZlenC();
  zLen[3]=zLen[4]-AliTOFGeometry::ZlenB();
  zLen[2]=zLen[3]-AliTOFGeometry::ZlenA();
  zLen[1]=zLen[2]-AliTOFGeometry::ZlenB();
  zLen[0]=zLen[1]-AliTOFGeometry::ZlenC();

  
  Int_t isector;    // random sector number 
  Int_t iplate;     // random plate number
  Int_t istrip;     // random strip number in the plate
  Int_t ipadAlongX; // random pad number along x direction
  Int_t ipadAlongZ; // random pad number along z direction
  Int_t ipad;
  Int_t nPixEdge=0; // additional pads fired due to edge effect when adding noise from outer
                    // regions

  // x -> time of flight given in ns
  TF1 *noiseTof = new TF1("noiseTof","exp(-x/20)",0,100);

  if(strstr(option,"pp")){
    fZnoise = new TF1("fZnoise","257.8-0.178*x-0.000457*x*x",-AliTOFGeometry::MaxhZtof(),AliTOFGeometry::MaxhZtof());
  }
  if(strstr(option,"Pb-Pb")){
    fZnoise = new TF1("fZnoise","182.2-0.09179*x-0.0001931*x*x",-AliTOFGeometry::MaxhZtof(),AliTOFGeometry::MaxhZtof());
  }

  if(fNoise) {
    if(fdbg) cout << " Start adding additional noise  hits from outer regions" << endl;

    for(Int_t i=0;i<fNoise;i++) {

      isector=(Int_t) (AliTOFGeometry::NSectors()*gRandom->Rndm())+1; //the sector number
      //  non-flat z-distribution of additional hits
      Float_t zNoise=fZnoise->GetRandom();

      // holes for PHOS and HMPID
      if(((AliTOF *) gAlice->GetDetector("TOF"))->IsVersion()==2) {
	// to be checked the holes case
	if(isector>12 && isector<16) { // sectors 13,14,15 - RICH
	  do {
	    iplate=(Int_t) (AliTOFGeometry::NPlates()*gRandom->Rndm())+1;
	  } while (iplate==2 || iplate==3 || iplate==4);
	  //         } else if(isector>11 && isector<17) { // sectors 12,13,14,15,16 - PHOS
	} else if(isector>2 && isector<8) { // sectors 3,4,5,6,7 - PHOS
	  do {
	    iplate=(Int_t) (AliTOFGeometry::NPlates()*gRandom->Rndm())+1;
	  } while (iplate==3);
	} else {
	  iplate=(Int_t) (AliTOFGeometry::NPlates()*gRandom->Rndm())+1;
	}
      } else {
	iplate=0;
	do {
	  iplate++;
	} while(zNoise>zLen[iplate]);
      }
      // end of holes

      if(iplate<1 || iplate>5) {
	printf("  iplate<1 or iplate>5, iplate=%i\n",iplate);
	return; 
      }

      Float_t nStripes=0;
      if(iplate>1) {
	for (Int_t i=0;i<iplate-1;i++) {
	  nStripes += zStrips[i];
	}
      }

      istrip=(Int_t)((zNoise-zLen[iplate-1])/((zLen[iplate]-zLen[iplate-1])/zStrips[iplate-1])); //the strip number in the plate
      istrip++;

      ipadAlongX = (Int_t)(AliTOFGeometry::NpadX()*gRandom->Rndm())+1;
      ipadAlongZ = (Int_t)(AliTOFGeometry::NpadZ()*gRandom->Rndm())+1;
      ipad=(Int_t)(ipadAlongZ-1)*AliTOFGeometry::NpadX()+ipadAlongX;    //the pad number
      
      Float_t xStrip=(ipadAlongX-1)*AliTOFGeometry::XPad()+AliTOFGeometry::XPad()*gRandom->Rndm()-0.5*AliTOFGeometry::NpadX()*AliTOFGeometry::XPad();//x-coor.in the strip frame
      Float_t zStrip=(ipadAlongZ-1)*AliTOFGeometry::ZPad()+AliTOFGeometry::ZPad()*gRandom->Rndm()-0.5*AliTOFGeometry::NpadZ()*AliTOFGeometry::ZPad();//z-coor.in the strip frame 

      Int_t nPlace[4]={0,0,0,0};
      nPlace[0]=ipad;

      Int_t   nActivatedPads=0;
      Int_t   nFiredPads=0;
      Bool_t  isFired[4]={kFALSE,kFALSE,kFALSE,kFALSE};
      Float_t tofAfterSimul[4]={0.,0.,0.,0.};
      Float_t qInduced[4]={0.,0.,0.,0.};
      Float_t averageTime=0.;    
      Float_t toffornoise=10.+noiseTof->GetRandom(); // 10 ns offset + parameterization [ns]

      BorderEffect(zStrip,xStrip,toffornoise,nActivatedPads,nFiredPads,isFired,nPlace,qInduced,tofAfterSimul,averageTime); // simulate edge effect

      if(nFiredPads) {
	for(Int_t indexOfPad=0; indexOfPad<nActivatedPads; indexOfPad++) {
	  if(isFired[indexOfPad]){// the pad has fired

	    if(indexOfPad==0) {// the hit belongs to a fired pad
	      isHitOnFiredPad++;
	      hitArray[isHitOnFiredPad-1].SetX(0.);
	      hitArray[isHitOnFiredPad-1].SetY(0.);
	      hitArray[isHitOnFiredPad-1].SetZ(zNoise);
	      hitArray[isHitOnFiredPad-1].SetNoise(1);
	    } // close if(indexOfPad==0)

	    ipad = nPlace[indexOfPad];
	    
	    Int_t iMapValue=MapPixels[AliTOFGeometry::NSectors()*(iplate-1)+isector-1][istrip-1][ipad-1];
	    
	    if(iMapValue==0) {
	      ipixel++;
	      if(indexOfPad) nPixEdge++;
	      MapPixels[AliTOFGeometry::NSectors()*(iplate-1)+isector-1][istrip-1][ipad-1]=ipixel;
	      pixelArray[ipixel-1].SetGeom(isector,iplate,istrip,ipad);
	      pixelArray[ipixel-1].AddState(1);
	      pixelArray[ipixel-1].SetRealTime(tofAfterSimul[indexOfPad]);
	      pixelArray[ipixel-1].SetHit(isHitOnFiredPad);
	    } else if( tofAfterSimul[indexOfPad] < pixelArray[iMapValue-1].GetRealTime() ) {
	      pixelArray[iMapValue-1].SetTrack(-1);
	      pixelArray[iMapValue-1].AddState(1);
	      pixelArray[iMapValue-1].SetRealTime(tofAfterSimul[indexOfPad]);
	      pixelArray[iMapValue-1].SetHit(isHitOnFiredPad);
	    }  //end of if(iMapValue==0)
	    
	  }// close if(isFired[indexOfPad])
	}  //end loop on activated pads indexOfPad
      } // close if(nFiredPads)
    }  //end of NOISE cycle
  }

  // free used memory
  if (fZnoise)
    {
      delete fZnoise;
      fZnoise = 0;
    }

  if (noiseTof)
    {
      delete noiseTof;
       noiseTof = 0;
    }

  Int_t nNoiseSignals=0;
  Int_t nAll=0;
  for(Int_t idummy=1; idummy<ipixel+1; idummy++) {
    if(hitArray[pixelArray[idummy-1].GetHit()-1].GetNoise()==1) {
      nNoiseSignals++;
      if(pixelArray[idummy-1].GetState()>10) nAll++;
    }
  }

  if(fdbg) {
    cout << " after adding " << fNoise << " noise hits: " << ipixel << " fired pixels (" << nNoiseSignals << " fired by noise, " << nAll << " noise+track)" << endl;
    printf(" %i additional pixels are fired by noise due to edge effect\n",nPixEdge);
    cout << " End of adding additional noise hits from outer regions" << endl;
  }

  Float_t occupancy;
  // numberOfPads for AliTOFV4 (Full coverage) 
  // - to be upgraded checking the used TOF version -
  Float_t numberOfPads=AliTOFGeometry::NPadXSector()*AliTOFGeometry::NSectors();
  occupancy=100.*ipixel/numberOfPads;   // percentage of fired pads
  printf(" Overall TOF occupancy (percentage of fired pads after adding noise) = %f\n",occupancy); 
  delete [] zLen;
  delete [] zStrips;
  
}


//____________________________________________________________________________
void AliTOFReconstructioner::SetMinDistance(AliTOFRecHit* hitArray, Int_t ilastEntry)
{
  //
  // Set the distance to the nearest hit for hitArray
  // ilastEntry is the index of the last entry of hitArray

  // starting the setting for the distance to the nearest TOFhit (cm)
  for(Int_t i=0; i<ilastEntry; i++) {
    
    if(hitArray[i].GetFirst()==1 && hitArray[i].GetNoise()==0) { // select the first hit of the track 
      // hits are not due to noise
      Float_t minDistance=10000.,squareDistance; // current values of the (square) distance
      Int_t jAtMin=0;                            // index of the hit nearest to the i-th hit
      Float_t xhit=hitArray[i].X(); // x coordinate for i-th hit
      Float_t yhit=hitArray[i].Y(); // y coordinate for i-th hit
      Float_t zhit=hitArray[i].Z(); // z coordinate for i-th hit
      //  was    for(Int_t j=0; j<isHitOnFiredPad; j++) {
      for(Int_t j=0; j<ilastEntry; j++) {
	if(i!=j) {
	  squareDistance=(hitArray[j].X()-xhit)*(hitArray[j].X()-xhit)+
	    (hitArray[j].Y()-yhit)*(hitArray[j].Y()-yhit)+
	    (hitArray[j].Z()-zhit)*(hitArray[j].Z()-zhit);
	  if(squareDistance<minDistance) {
	    minDistance=squareDistance;
	    jAtMin=j;
	  }
	}
      }
      minDistance=TMath::Sqrt(minDistance);
      hitArray[i].SetRmin(minDistance);
      if(minDistance==0.) printf(" Rmin=0, i=%i, j=%i, x=%f,y=%f,z=%f\n",i,jAtMin,xhit,yhit,zhit);// it cannot happen
    }
  }

}

// these lines has to be commented till TPC will provide fPx fPy fPz 
// and fL in AliTPChit class
//____________________________________________________________________________ 
/*
void AliTOFReconstructioner::ReadTPCHits(Int_t ntracks, TTree* treehits, TClonesArray* tpchits, Int_t* iTrackPt, Int_t* iparticle, Float_t* ptTrack, AliTOFTrack* trackArray, Int_t& itrack)
{
  //
  // Read TPC hits for the current event
  // 
  TParticle *particle=0;
  Int_t npions=0;    // number of pions for the current event
  Int_t nkaons=0;    // number of kaons for the current event
  Int_t nprotons=0;  // number of protons for the current event
  Int_t nelectrons=0;// number of electrons for the current event
  Int_t nmuons=0;    // number of muons for the current event
  Int_t ntotalTPChits=0; // total number of TPC hits for the current event
  Int_t idummy=-1;       // dummy var used to count double hit TPC cases
  Int_t nTpcDoubleHitsLastRow=0; // number of double TPC hits in the last pad row
  Int_t nTpcHitsLastRow=0;       // number of TPC hits in the last pad row
  Float_t trdpos[2]={0.,0.};
  Float_t pos[3];               // TPC hit position
  Float_t mom[3]; // momentum components in the last TPC row
  Float_t pt=0., tpclen; // pt: transverse momentum in the last TPC row
  Int_t nbytes = 0;
  Int_t ipart=0, nhits=0, iprim=0;

  itrack=0; // itrack: total number of selected TPC tracks

  // speed-up the code
  treehits->SetBranchStatus("*",0); // switch off all branches
  treehits->SetBranchStatus("TPC*",1); // switch on only TPC

  for (Int_t track=0; track<ntracks;track++) {
    gAlice->ResetHits();
    nbytes += treehits->GetEvent(track);
    
    
    nhits = tpchits->GetEntriesFast();
    
    for (Int_t hit=0;hit<nhits;hit++) {
      ntotalTPChits++;
      AliTPChit* tpcHit = (AliTPChit*)tpchits->UncheckedAt(hit);
      Int_t row = tpcHit->fPadRow;
      ipart    = tpcHit->GetTrack();
      if(ipart>=fMaxAllTracks) break;
      particle = (TParticle*)gAlice->Particle(ipart);
      Int_t pdgCode=particle->GetPdgCode();
      // only high momentum tracks
      // momentum components at production vertex
      Float_t pxvtx = particle->Px();
      Float_t pyvtx = particle->Py();
      Float_t pzvtx = particle->Pz();
      Float_t pvtx = TMath::Sqrt(pxvtx*pxvtx+pyvtx*pyvtx+pzvtx*pzvtx);
      if(pvtx>fPBound && row == fLastTPCRow) {
	Float_t vx = particle->Vx();
	Float_t vy = particle->Vy();
	Float_t vr = TMath::Sqrt(vx*vx+vy*vy);
	Float_t x = tpcHit->X();
	Float_t y = tpcHit->Y();
	Float_t z = tpcHit->Z();
	pos[0]=x; pos[1]=y; pos[2]=z;
	
	Float_t pxtpc = tpcHit->fPx;
	Float_t pytpc = tpcHit->fPy;
	Float_t pztpc = tpcHit->fPz;
	mom[0]=pxtpc; mom[1]=pytpc; mom[2]=pztpc; 
	Float_t momtpc = TMath::Sqrt(pxtpc*pxtpc+pytpc*pytpc+pztpc*pztpc);
	
	if(x*pxtpc+y*pytpc>0) { // only tracks going out of TPC
	  
	  Float_t isoutgoing = x*pxtpc+y*pytpc+z*pztpc;
	  isoutgoing /= (momtpc*TMath::Sqrt(x*x+y*y+z*z));
	  tpclen = tpcHit->fL;
	  
	  
	  if(ipart!=idummy) {
	    if(particle->GetFirstMother() < 0) {
	      Int_t abspdgCode=TMath::Abs(pdgCode);
	      switch (abspdgCode) {
	      case 211:
		npions++;
		break ;
	      case 321:
		nkaons++;
		break ;
	      case 2212:
		nprotons++;
		break ;
	      case 11:
		nelectrons++;
		break ;
	      case 13:
		nmuons++;
		break ;
	      }
	    } // close if(particle->GetFirstMother() < 0)
	  } // close if(ipart!=idummy)
	  
	  if(gRandom->Rndm()<fTrackingEfficiency && vr<fRadiusvtxBound && ipart!=idummy) {
	    
	    itrack++;
	    if(particle->GetFirstMother() < 0) iprim++;
	    
	    if(itrack>fMaxTracks) {
	      cout << "itrack=" << itrack << " > MAXTRACKS=" << fMaxTracks << endl;
	      return;
	    } // close if(itrack>fMaxTracks)
	    
	    
	    iparticle[ipart]=itrack;
	    
	    trackArray[itrack-1].SetTrack(ipart,pvtx,pdgCode,tpclen,pos,mom,trdpos);
	    
	    pt=TMath::Sqrt(pxtpc*pxtpc+pytpc*pytpc); // pt: transverse momentum at TPC
	    // Filling iTrackPt[MAXTRACKS] by itrack ordering on Pt
	    if(itrack==1) {
	      iTrackPt[itrack-1]=itrack;
	      ptTrack[itrack-1]=pt;
	    } else {
	      for (Int_t i=0; i<itrack-1; i++) {
		if(pt>ptTrack[i]) {
		  for(Int_t j=i; j<itrack-1; j++) {
		    Int_t k=itrack-1+i-j;
		    iTrackPt[k]= iTrackPt[k-1];
		    ptTrack[k] = ptTrack[k-1];
		  }
		  iTrackPt[i]=itrack;
		  ptTrack[i]=pt;
		  break;
		}
		if(i==itrack-2) {
		  iTrackPt[itrack-1]=itrack;
		  ptTrack[itrack-1]=pt;
		}
	      }
	    }
	    
	  }  //end of itrack
	  if(vr>fRadiusvtxBound) nTpcHitsLastRow++;
	  if(ipart==idummy) nTpcDoubleHitsLastRow++;
	  idummy=ipart;
	}  // close if(x*px+y*py>0)
      }  // close if(pvtx>fPBound && row == fLastTPCRow)
    }  //end of hits  
  }  // close loop on tracks   
  
  
  if(fdbg) {
    cout << ntotalTPChits << " TPC hits in the last TPC row " << fLastTPCRow << endl;
    cout << "   " << nTpcHitsLastRow << " - hits with Rvtx>fRadiusvtxBound=" << fRadiusvtxBound << endl;
    cout << "   " << nTpcDoubleHitsLastRow << " double TPC hits" << endl;
    cout << itrack    << " - extracted TPC tracks   "     << iprim << " - primary" << endl;
    cout << npions    << " primary pions reached TPC"     << endl;
    cout << nkaons    << " primary kaons reached TPC"     << endl;
    cout << nprotons  << " primary protons reached TPC"   << endl;
    cout << nelectrons<< " primary electrons reached TPC" << endl;
    cout << nmuons    << " primary muons reached TPC"     << endl;
  } // if(fdbg)
  
  Int_t primaryInTPC[6]={0,0,0,0,0,0};
  primaryInTPC[0]=npions;
  primaryInTPC[1]=nkaons;
  primaryInTPC[2]=nprotons;
  primaryInTPC[3]=nelectrons;
  primaryInTPC[4]=nmuons;
  primaryInTPC[5]=npions+nkaons+nprotons+nelectrons+nmuons;
  
  if(fdbg) {
    printf("  contents of iTrackPt[MAXTRACKS],PtTrack[MAXTRACKS]\n");
    for (Int_t i=0; i<itrack; i++) {
      printf(" %i : iTrackPt=%i, PtTrack=%f\n",i+1,iTrackPt[i],ptTrack[i]); 
    }
    printf(" Check ordered transverse momentum array\n");
    for (Int_t i=itrack-1; i>=0; i--) {
      printf(" %i : iTrackPt=%i, PtTrack=%f\n",i+1,iTrackPt[i],ptTrack[i]); 
    }
  }// if(fdbg)
  
}
*/
//____________________________________________________________________________
void cylcor(Float_t& x, Float_t& y) {
  Float_t rho,phi;
  
  rho=TMath::Sqrt(x*x+y*y);
  phi=0.;
  if(TMath::Abs(x)>0. || TMath::Abs(y)>0.) phi=TMath::ATan2(y,x);
  if(phi<0.) phi=phi+2.*TMath::Pi();
  x=rho;
  y=phi;
  
}

//____________________________________________________________________________
void AliTOFReconstructioner::Matching(AliTOFTrack* trackArray, AliTOFRecHit* hitArray, Int_t ***mapPixels, AliTOFPad* pixelArray, Int_t* kTOFhitFirst, Int_t& ipixel, Int_t* iTrackPt, Int_t* iTOFpixel, Int_t ntotTpcTracks)
{
  Int_t TestTracks,iTestTrack,itest,wPixel=0,itestc;
  Int_t * ntest = new Int_t[fMaxTestTracks];
  Int_t * testPixel = new Int_t[fMaxTestTracks];
  Float_t wLength=0.,wRho=0.,wZ=0.;
  Float_t * testLength = new Float_t[fMaxTestTracks];
  Float_t * testRho = new Float_t[fMaxTestTracks];
  Float_t * testZ = new Float_t[fMaxTestTracks];
  Float_t weight;
  Float_t * testWeight = new Float_t[fMaxTestTracks];
  Float_t rotationFactor,phi0,coslam,sinlam,helixRadius,xHelixCenter,yHelixCenter,zHelixCenter,helixFactor;
  Int_t npixel[5],iMapValue,iwork1,iwork2,iwork3,iwork4,ihit=0;
  Int_t charge[48]={ 0, 1,-1, 0, 1,-1, 0, 1,-1, 0,
		     1,-1, 0, 1,-1, 0, 0, 0, 1, 0,
                     -1, 0,-1,-1, 0, 0,-1, 0, 1, 0,
		     1, 1, 0, 0, 1,-1, 0, 0, 1,-1,
		     1, 1,-1, 0, 1, 1, 2, 0};
  Float_t theta0,gpx,gpy,gpz,gp,gpt,gtheta,gx,gy,gz,gr,gxLast,gyLast,gzLast,chargeField;
  Float_t sumOfTheta=0.,weightTestTracksOutTof[4];
  Float_t s,ds,xRespectToHelixCenter,yRespectToHelixCenter,deltaRadius,fp,xp,yp,grho;
  Float_t mass,energy,g;
  Int_t itrack=0,itr,particleCharge,istep,iplate=0,iPadAlongX=0;  
  Int_t itra,t34=0,t32=0,t44=0,t43=0,t42=0;
  Int_t wstate=0,m2state=0,wPix;
  Int_t idelR=0,idelR1=0,idelR2=0,iRmin=0,iRmin1=0,iRmin2=0;
  Float_t massArray[50] = {0.0,0.00051,0.00051,0.0,0.1057,0.1057,0.135,0.1396,0.1396,0.4977,
		       0.4936,0.4936,0.9396,0.9383,0.9383,0.4977,0.5488,1.1156,1.1894,1.1926,1.1926,
		       1.3149,1.3213,1.6724,0.9396,1.1156,1.1894,1.1926,1.1974,1.3149,
		       0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Float_t delR;
  Float_t radius,area,normR,normS,cosAngl;
  Int_t iPlateFirst,iTestGmax=0;
  Int_t fstate,iPrintM1=0,iPrintM2=0;
  Float_t gxExtrap=0.,gyExtrap=0.,gzExtrap=0.;
  Float_t avSigZ=0,avSigRPHI=0,avSigP=0,avSigPHI=0,avSigTHETA=0;

  Float_t gxW,gyW,gzW;
  Float_t length0;
  Float_t snr=0;
  Int_t indexOfTestTrack;
  Float_t zPad,xPad;
  Int_t istate=0,imax=0,match,iMaxTestTracksOutTof=0,matchw;
  Float_t w,wmax=0.,inverseOfParticleSpeed,w2,smat[9],largestWeightTracksOutTof,sw;
  Float_t sumWeightTracksOutTof,sGeomWeigth;
  Int_t imatched;
  Int_t m10=0,m20=0,m22=0,m23=0;
  Int_t PRINT=0;
  TParticle *particle;

  Float_t time=0.;
  itr=ntotTpcTracks;
  printf(" itr=%i\n",itr);
  for (itra=1; itra<itr+1; itra++) {

    Int_t itrack=iTrackPt[itra-1];
    if(itrack==0) printf("  iTrackPt[itra-1]=0 for itra=%i\n",itra);
    Int_t ipart=trackArray[itrack-1].GetTrack(); 
    Float_t pvtx=trackArray[itrack-1].GetP();
    Int_t pdgCode=trackArray[itrack-1].GetPdgCode();
    Float_t tpclength=trackArray[itrack-1].GetlTPC();
    Float_t x=trackArray[itrack-1].GetRxTPC();
    Float_t y=trackArray[itrack-1].GetRyTPC();
    Float_t z=trackArray[itrack-1].GetRzTPC();
    /* vars used for QA
    Float_t RxTPC=x;
    Float_t RyTPC=y;
    Float_t RzTPC=z;
    */
    Float_t Wx=x;
    Float_t Wy=y;
    Float_t Wz=z;
    Float_t px=trackArray[itrack-1].GetPxTPC();
    Float_t py=trackArray[itrack-1].GetPyTPC();
    Float_t pz=trackArray[itrack-1].GetPzTPC();
    /* vars used for QA
    Float_t pxTPC=px;
    Float_t pyTPC=py;
    Float_t pzTPC=pz;
    */
    Float_t p = TMath::Sqrt(px*px+py*py+pz*pz);
    /* var used for QA
    Float_t pTPC=p;
    */
    Float_t rho = TMath::Sqrt(x*x+y*y);
    Float_t phi=0.;
    if(TMath::Abs(x)>0. || TMath::Abs(y)>0.) phi=TMath::ATan2(y,x);
    if(phi<0.) phi=phi+2.*TMath::Pi();
    /* var used for QA
    Float_t phiTPC=phi*kRaddeg;
    */
    if(fSigmavsp) {
      if(p==0) printf(" p=%f in g=0.022/p\n",p);
      g=0.022/p;
      avSigRPHI += g;      // (cm)
      if(rho==0) printf(" rho=%f in phi += g*gRandom->Gaus()/rho\n",rho);
      phi += g*gRandom->Gaus()/rho; 
    } else {
      if(rho==0) printf(" rho=%f in phi += (SIGMARPHI*gRandom->Gaus()/rho\n",rho);
      phi += (fSigmarphi*gRandom->Gaus()/rho);
    }
    x=rho*TMath::Cos(phi);
    y=rho*TMath::Sin(phi);
    /* var used for QA
    Float_t zTPC=z;
    */
    if(fSigmavsp) {
      if(p==0) printf(" p=%f in g=0.0275/p\n",p);
      g=0.0275/p;
      avSigZ += g;      // (cm)
      z += g*gRandom->Gaus();
    } else {
      z += fSigmaZ*gRandom->Gaus();
    }

    // smearing on TPC momentum

    {                                                                             
      Float_t pmom,phi,theta,arg;
      
      pmom=TMath::Sqrt(px*px+py*py+pz*pz);
      phi=0.;
      if(TMath::Abs(px)>0. || TMath::Abs(py)>0.) phi=TMath::ATan2(py,px);
      if(phi<0.) phi=phi+2*TMath::Pi();
      arg=1.;
      if(pmom>0.) arg=pz/pmom;
      theta=0.;
      if(TMath::Abs(arg)<=1.) theta=TMath::ACos(arg);
      
      if(fSigmavsp) {
        if(pmom<=0) printf(" pmom=%f in g = TMath::Abs(TMath::Log(pmom)/TMath::Log(10)+0.5)/0.7\n",pmom);
        g = TMath::Abs(TMath::Log(pmom)/TMath::Log(10)+0.5)/0.7;
        g = 0.01*(g*g*g+1.5)*1.24;
        avSigP += g;
        pmom *= (1+g*gRandom->Gaus());
	
        if(p<10) {
	  if(pmom<=0) printf(" pmom=%f in g = 1-TMath::Log(pmom)/TMath::Log(10)\n",pmom);
          g = 1-TMath::Log(pmom)/TMath::Log(10);
          g = 0.001*(g*g*g+0.3)*0.65;  // (radian)
	} else {
          g = 0.001*0.3*0.65;
	}
        avSigPHI += g;
        phi += g*gRandom->Gaus();
        avSigTHETA += g;
        theta += g*gRandom->Gaus();
        
      } else {
        pmom *= (1+fSigmap*gRandom->Gaus());
        phi += fSigmaPhi*gRandom->Gaus();
        theta += fSigmaTheta*gRandom->Gaus();
      }
      gxW=px;
      gyW=py;
      gzW=pz;
      
      px=pmom*TMath::Sin(theta)*TMath::Cos(phi);
      py=pmom*TMath::Sin(theta)*TMath::Sin(phi);
      pz=pmom*TMath::Cos(theta);

      
      if(x*px+y*py<=0) {
        x=Wx;
        y=Wy;
        z=Wz;
        px=gxW;
        py=gyW;
        pz=gzW;
      }// if(x*px+y*py<=0)
    }
    
    p = TMath::Sqrt(px*px+py*py+pz*pz);
    
    particleCharge=charge[PDGtoGeantCode(pdgCode)-1];
    mass=massArray[PDGtoGeantCode(pdgCode)-1];
    mass=massArray[8-1];       //we take pion mass for all tracks
    //             mass=massArray[14-1];       //here we take proton mass for all tracks
    energy=TMath::Sqrt(p*p+mass*mass);
    chargeField=particleCharge*fField;
    
    g=fRadLenTPC/( (x*px+y*py)/(rho*p) );
    
    if(g<=0) printf(" error, g<=0: g=%f, itra=%i, x,y,px,py=%f, %f, %f, %f\n",g,itra,x,y,px,py);
    
    theta0=13.6*0.001*TMath::Sqrt(g)*(1.+0.038*TMath::Log(g))*energy/(p*p);
 
    
    // start Loop on test tracks
    sumOfTheta=0.;
    for(Int_t i=0;i<4;i++) {
      weightTestTracksOutTof[i]=0.;
    }
    
    itest=0;
    for(Int_t i=0;i<fMaxTestTracks;i++) {
      ntest[i]=0;
      testPixel[i]=0;
      testLength[i]=0.;
      testRho[i]=0.;
      testZ[i]=0.;
      testWeight[i]=0.;
    }
    
    iPlateFirst=0;
    TestTracks=0;
    iTestTrack=0;
    iTestGmax=0;
    
    length0=0;
    
    for (indexOfTestTrack=0; indexOfTestTrack<fMaxTestTracks; indexOfTestTrack++) {

      iTestTrack++;
      gpx=px;
      gpy=py;
      gpz=pz;
      gp=p;
      if(indexOfTestTrack) {
	gtheta=theta0;
	EpMulScatt(gpx,gpy,gpz,gp,gtheta);
	
      } else {
	gtheta=0;
      }
      
      weight=TMath::Exp(-gtheta*gtheta/(2*theta0*theta0));
      sumOfTheta += gtheta;
      
      //    ==========================================================
      // Calculate crossing of the track in magnetic field with cylidrical surface
      // of radius RTOFINNER
      //   chargeField = qB, where q is a charge of a particle in units of e,
      //                     B is magnetic field in tesla
      //   see 3.3.1.1. in the book "Data analysis techniques for
      //   high-energy physics experiments", edited by M.Regler
      //   in Russian: "Metody analiza dannykh v fizicheskom eksperimente"
      //   Moskva, "Mir", 1993. ctr.306
      
      // Initial constants
      rotationFactor=1.;
      if(chargeField<0.) rotationFactor=-1.;
      rotationFactor=-rotationFactor;
      gpt=gpx;
      phi0=gpy;
      cylcor(gpt,phi0);
      phi0 -= rotationFactor*TMath::Pi()*0.5;
      //               phi0 -= h*PID2;
      coslam=gpt/gp;
      sinlam=gpz/gp;
      //      helixRadius=100.*gpt/TMath::Abs(0.299792458*chargeField);
      helixRadius=100.*gpt/TMath::Abs(AliTOFGeometry::SpeedOfLight()*chargeField);
      xHelixCenter=x-helixRadius*TMath::Cos(phi0);
      yHelixCenter=y-helixRadius*TMath::Sin(phi0);
      zHelixCenter=z;
      helixFactor=rotationFactor*coslam/helixRadius;
      
      //   Solves the equation f(s)=r(s)-RTOFINNER=0 by the Newton's method:
      //   snew=s-f/f'
      istep=0;
      s=AliTOFGeometry::Rmin()-TMath::Sqrt(x*x+y*y);;
      do {
	istep++;
	xRespectToHelixCenter=helixRadius*TMath::Cos(phi0+s*helixFactor);
	yRespectToHelixCenter=helixRadius*TMath::Sin(phi0+s*helixFactor);
	gx=xHelixCenter+xRespectToHelixCenter;
	gy=yHelixCenter+yRespectToHelixCenter;
	gr=TMath::Sqrt(gx*gx+gy*gy);
	deltaRadius=gr-AliTOFGeometry::Rmin();
	xp=-helixFactor*yRespectToHelixCenter;
	yp= helixFactor*xRespectToHelixCenter;
	fp=(gx*xp+gy*yp)/gr;
	ds=deltaRadius/fp;
	s -= ds;
	if(istep==20) {
	  istep=0;
	  break;
	}
      } while (TMath::Abs(ds)>0.01);
      
      
      if(istep==0) goto end;
      
      //   Steps along the circle till a pad
      wPixel=0;
      wLength=0.;
      iplate=0;
      iPadAlongX=0;
      grho=0.;
      ds=fStep;
      gxLast=xHelixCenter+helixRadius*TMath::Cos(phi0+s*helixFactor);
      gyLast=yHelixCenter+helixRadius*TMath::Sin(phi0+s*helixFactor);
      gzLast=zHelixCenter+s*sinlam;

      
      do {
	istep++;
	s += ds;
	gx=xHelixCenter+helixRadius*TMath::Cos(phi0+s*helixFactor);
	gy=yHelixCenter+helixRadius*TMath::Sin(phi0+s*helixFactor);
	gz=zHelixCenter+s*sinlam;
	rho=TMath::Sqrt(gx*gx+gy*gy);
	
	IsInsideThePad(gMC,gx,gy,gz,npixel,zPad,xPad);
	
	iplate += npixel[1];
	iPadAlongX += npixel[4];
	
	if(indexOfTestTrack==0 && iplate && iPlateFirst==0) {
	  iPlateFirst=1;
	  length0=s;

	  radius=s*3*theta0;
	  area=TMath::Pi()*radius*radius;
	  normR=TMath::Sqrt(gx*gx+gy*gy);
	  normS=TMath::Sqrt((gx-gxLast)*(gx-gxLast)+
		     (gy-gyLast)*(gy-gyLast)+
		     (gz-gzLast)*(gz-gzLast));

	  cosAngl=(gx*(gx-gxLast)+gy*(gy-gyLast))/(normR*normS);
	  if(cosAngl<0) printf(" cosAngl<0: gx=%f,gy=%f,  gxLast=%f,gyLast=%f,gzLast=%f\n",gx,gy,gxLast,gyLast,gzLast);

	  area /= cosAngl;
	  TestTracks=(Int_t) (2*area/(AliTOFGeometry::XPad() * AliTOFGeometry::ZPad()));

	  if(TestTracks<12) TestTracks=12;

	  // Angles of entering into the TOF plate

	  Int_t iZ=0;
	  if(TMath::Abs(gz)>300) {
	    iZ=4;
	  } else if(TMath::Abs(gz)>200) {
	    iZ=3;
	  } else if(TMath::Abs(gz)>100) {
	    iZ=2;
	  } else if(TMath::Abs(gz)>0) {
	    iZ=1;
	  }
	  
	  
	} // end of if(indexOfTestTrack==0 && iplate && iPlateFirst==0)

	
	if(npixel[4]>0) {

	  iwork1=npixel[0];
	  iwork2=npixel[1];
	  iwork3=npixel[2];
	  //		   iwork4=npixel[3];
	  iwork4=(npixel[3]-1)*AliTOFGeometry::NpadX()+npixel[4];

	  Int_t ifirstindex=AliTOFGeometry::NSectors()*(npixel[1]-1)+npixel[0];
	  iMapValue=mapPixels[ifirstindex-1][iwork3-1][iwork4-1];
	  if(iMapValue==0) {
	    ipixel++;
	    if(ipixel>fMaxPixels) {
	      cout << "ipixel=" << ipixel << " > MAXPIXELS=" << fMaxPixels << endl;
	      break;
	    }
	    mapPixels[ifirstindex-1][iwork3-1][iwork4-1]=ipixel;
	    pixelArray[ipixel-1].SetGeom(iwork1,iwork2,iwork3,iwork4);
	    iMapValue=ipixel;
	  }
	  
	  wPixel=iMapValue;
	  wLength=tpclength+s;
	  wRho=rho;
	  wZ=gz;
	  
	  ihit=kTOFhitFirst[ipart];
	  
	  if(ihit) {
	    if(indexOfTestTrack==0) {
	      {
		idelR++;
		delR=TMath::Sqrt((gx-hitArray[ihit-1].X())*(gx-hitArray[ihit-1].X())+
			  (gy-hitArray[ihit-1].Y())*(gy-hitArray[ihit-1].Y())+
			  (gz-hitArray[ihit-1].Z())*(gz-hitArray[ihit-1].Z()));

	      }
	      
	      if(delR>hitArray[ihit-1].GetRmin()) iRmin++;
	      gxExtrap=gx;
	      gyExtrap=gy;
	      gzExtrap=gz;
	    } else {
	      delR=TMath::Sqrt((gx-gxExtrap)*(gx-gxExtrap)+
			(gy-gyExtrap)*(gy-gyExtrap)+
			(gz-gzExtrap)*(gz-gzExtrap));
	    }
	  }  //end of if(ihit)
	  
	  break;
	  
	}  //end of npixel[4]
	
	if(rho<grho) {
	  istep=0;
	  break;
	}
	grho=rho;
	
	gxLast=gx;
	gyLast=gy;
	gzLast=gz;
	
      } while(rho<AliTOFGeometry::Rmax()); //end of do 

      
      if(istep>0) {
	if(iplate) {
	  if(iPadAlongX==0) {
	    istep=-3;            // holes in TOF
	  }
	} else {
	  if(TMath::Abs(gz)<AliTOFGeometry::MaxhZtof()) {
	    //                   if(TMath::Abs(gz)<MAXZTOF2) {
	    istep=-2;            // PHOS and RICH holes or holes in between TOF plates
	  } else {
	    istep=-1;            // out of TOF on z-size
	  }
	}
      }
      
      if(iPadAlongX>0) {
	if(itest==0) {
	  itest=1;
	  ntest[itest-1]=1;
	  testPixel[itest-1]=wPixel;
	  testLength[itest-1]=wLength;
	  testRho[itest-1]=wRho;
	  testZ[itest-1]=wZ;
	  testWeight[itest-1]=weight;
	} else {
	  Int_t k=0;
	  for(Int_t i=0;i<itest;i++) {
	    k=0;
	    if(testPixel[i]==wPixel) {
	      k=1;
	      ntest[i]++;
	      testLength[i] += wLength;
	      testRho[i] += wRho;
	      testZ[i] += wZ;
	      testWeight[i] += weight;
	      break;
	    }
	  }  //end for i
	  if(k==0) {
	    itest++;
	    ntest[itest-1]=1;
	    testPixel[itest-1]=wPixel;
	    testLength[itest-1]=wLength;
	    testRho[itest-1]=wRho;
	    testZ[itest-1]=wZ;
	    testWeight[itest-1]=weight;
	  }
	}
      }
      
    end: ;
      //   Statistics
      if(fMatchingStyle==1) {
	if(istep>-4 && istep<1) weightTestTracksOutTof[-istep] ++;
      } else {
	if(istep>-4 && istep<1) weightTestTracksOutTof[-istep] += weight;
      }
      
      if(fMatchingStyle==2) {
	if(indexOfTestTrack==0 && istep==0) break;
	if(indexOfTestTrack+1==TestTracks) break;
      }
      
    }  //end of indexOfTestTrack

    snr += (Float_t) (indexOfTestTrack+1);
    
    //   Search for the "hole" with the largest weigth
    largestWeightTracksOutTof=0.;
    sumWeightTracksOutTof=0.;
    for(Int_t i=0;i<4;i++) {
      w=weightTestTracksOutTof[i];
      sumWeightTracksOutTof += w;
      if(w>largestWeightTracksOutTof) {
	largestWeightTracksOutTof=w;
	iMaxTestTracksOutTof=i;
      }
    }
    
    itestc=itest;
    if(itest>0) {
      for(Int_t i=0;i<itest;i++) {
	testLength[i] /= ntest[i];
	testRho[i] /= ntest[i];
	testZ[i] /= ntest[i];
      }
      //   Search for the pixel with the largest weigth
      wmax=0.;
      wstate=0;
      sw=0;
      sGeomWeigth=0;
      for(Int_t i=0;i<itest;i++) {
	istate=pixelArray[testPixel[i]-1].GetState();
	fstate=0;
	if(istate>0) {
	  fstate=1;
	  wstate++;
	}
	if(fMatchingStyle==1) {
	  sGeomWeigth += ntest[i];
	  w=(fpadefficiency*fstate+(1.-fpadefficiency)*(1-fstate))*ntest[i];
	  if(pixelArray[testPixel[i]-1].GetTrackMatched()>0) w *= 0.1;
	} else {
	  sGeomWeigth += testWeight[i];
	  w=(fpadefficiency*fstate+(1.-fpadefficiency)*(1-fstate))*testWeight[i];
	  if(pixelArray[testPixel[i]-1].GetTrackMatched()>0) w *= 0.1;
	}
	
	// weighting according to the Pulse Height (we use the square of weight)
	// if (fChargeFactorForMatching) w *= (pixelArray[testPixel[i]-1].GetCharge())*(pixelArray[testPixel[i]-1].GetCharge());
	if (fChargeFactorForMatching && fstate==1) w *= (pixelArray[testPixel[i]-1].GetCharge())*(pixelArray[testPixel[i]-1].GetCharge());

	if(w>wmax) {
	  wmax=w;
	  imax=i;
	}
	sw += w;
      }
      wPixel=testPixel[imax];
      wLength=testLength[imax];
      istate=pixelArray[wPixel-1].GetState();
      
      //Choose the TOF dead space
      //               if(istate==0 && largestWeightTracksOutTof>wmax) {
      //               if(istate==0 && largestWeightTracksOutTof>=sw) {
      if(istate==0 && sumWeightTracksOutTof>sGeomWeigth) {
	itestc=itest;
	itest=0;
      }
    }
    
    if(itest>0) {
      
      //   Set for MyTrack: Pixel
      trackArray[itrack-1].SetPixel(wPixel);
      
      istate=pixelArray[wPixel-1].GetState();
      
      if(istate) {
	
	//   Set for MyTrack: Pixel, Length, TOF, MassTOF
	//fp
	//time=pixelArray[wPixel-1].GetTime();
	time=pixelArray[wPixel-1].GetRealTime();
	trackArray[itrack-1].SetLength(wLength);
	trackArray[itrack-1].SetTof(time);
	
	inverseOfParticleSpeed=time/wLength;
	//w=900.*inverseOfParticleSpeed*inverseOfParticleSpeed-1.;
	w=(100.*AliTOFGeometry::SpeedOfLight())*(100.*AliTOFGeometry::SpeedOfLight())*inverseOfParticleSpeed*inverseOfParticleSpeed-1.;
	w2=pvtx*pvtx;
	Float_t squareMass=w2*w;
	mass=TMath::Sqrt(TMath::Abs(squareMass));
	if(w<0.) mass=-mass;
	
	trackArray[itrack-1].SetMassTOF(mass);
	
	//   Set for MyTrack: Matching
	match=4;
	//                 if(ipart==pixelArray[wPixel-1].GetTrack()) match=3;
	if( (ipart==pixelArray[wPixel-1].GetTrack()) && hitArray[pixelArray[wPixel-1].GetHit()-1].GetNoise()==0)match=3;
	imatched=pixelArray[wPixel-1].GetTrackMatched();
	//   Set for TOFPixel the number of matched track
	pixelArray[wPixel-1].SetTrackMatched(itrack);
	
	if(imatched>0) {
	  matchw=trackArray[imatched-1].GetMatching();
	  if(match==3 && matchw==4) t34++;
	  if(match==3 && matchw==2) t32++;
	  if(match==4 && matchw==4) t44++;
	  if(match==4 && matchw==3) t43++;
	  if(match==4 && matchw==2) t42++;
	  if(iTOFpixel[ipart]==0 || iTOFpixel[trackArray[imatched-1].GetTrack()]==0) {
	    m20++;
	  } else if(iTOFpixel[ipart]==iTOFpixel[trackArray[imatched-1].GetTrack()]) {
	    m22++;
	  } else {
	    m23++;
	    wPix=iTOFpixel[ipart];
	    if(PRINT && iPrintM1==10 && iPrintM2<10) {
	      if(iPrintM2==0) {
		printf("*** test print for tracks matched with the pixel for with we had matched track\n");
	      }
	      iPrintM2++;
	      printf(" m=2: ipart=%i, pdgCode=%i, p=%f, theta0=%f, %i Pixel(LP=%i,SP=%i,P=%i) \n", 
		     ipart,pdgCode,p,theta0,wPix,
		     pixelArray[wPix-1].GetSector(),pixelArray[wPix-1].GetPlate(),pixelArray[wPix-1].GetPixel());
	      printf("      mat=%i, %i Pixel(LP=%i,SP=%i,P=%i), Test(n=%i,i=%i,w=%f,z=%f), wst=%i \n",
		     match,wPixel,
		     pixelArray[wPixel-1].GetSector(),pixelArray[wPixel-1].GetPlate(),pixelArray[wPixel-1].GetPixel(),
		     itest,imax,wmax,testZ[imax],wstate);
	      Int_t fstat,istat;
	      for(Int_t i=0;i<itest;i++) {
		wPix=testPixel[i];
		istat=pixelArray[wPix-1].GetState();
		fstat=0;
		if(istat>0) fstat=1;
		w=(fpadefficiency*fstat+(1.-fpadefficiency)*(1-fstat))*ntest[i];
		if(istat>0)
		  printf("                     %i: %i Pixel(LP=%i,SP=%i,P=%i), istat=%i, ntest=%i, w=%f\n",i+1,
			 wPix,pixelArray[wPix-1].GetSector(),pixelArray[wPix-1].GetPlate(),pixelArray[wPix-1].GetPixel(),
			 istat,ntest[i],w);
	      }
	      printf("      mat=%i, %i Pixel \n",matchw,trackArray[imatched-1].GetPad());
	    }
	  }
	  if(wstate>1) m2state++;
	  smat[matchw+4]--;
	  match=2;
	  trackArray[imatched-1].SetMatching(match);
	  smat[match+4]++;
	  
	}  // if(imatched>0)
	
      } else {  //else if(istate)
	
	match=1;
	if(iTOFpixel[ipart]==0) m10++;
	if(PRINT && iPrintM1<10) {
	  Int_t wPix;
	  wPix=iTOFpixel[ipart];
	  if(wPix) {
	    if(iPrintM1==0) {
	      printf("*** test print for tracks fired a pixel but matched with non-fired pixel\n");
	    }
	    iPrintM1++;
	    printf(" m=1: itra=%i,ipart=%i, pdgCode=%i, p=%f, theta0=%f, %i Pixel(LP=%i,SP=%i,P=%i) \n", 
		   itra,ipart,pdgCode,p,theta0,wPix,
		   pixelArray[wPix-1].GetSector(),pixelArray[wPix-1].GetPlate(),pixelArray[wPix-1].GetPixel());
	    printf("      mat=%i, %i Pixel(LP=%i,SP=%i,P=%i), Test(n=%i,i=%i,w=%f,z=%f), wst=%i \n",
		   match,wPixel,
		   pixelArray[wPixel-1].GetSector(),pixelArray[wPixel-1].GetPlate(),pixelArray[wPixel-1].GetPixel(),
		   itest,imax,wmax,testZ[imax],wstate);
	    
	  }
	}  //end if(PRINT && iPrintM1<10)
	
      }  //end if(istate)
      
    } else {
      match=-1-iMaxTestTracksOutTof;
      
    }  //end itest
    
    trackArray[itrack-1].SetMatching(match);
    //             if(iTestGmax==1) hMTT->Fill(match);
    smat[match+4]++;

    sumOfTheta /= iTestTrack;
    
    itest=itestc;
    
    //Test
    if(PRINT) {
      if(iTOFpixel[ipart] && match!=3) {
	particle = (TParticle*)gAlice->GetMCApp()->Particle(ipart);  //for V3.05

	printf("          ipixel=%i (Sector=%i, Plate=%i, Strip=%i, Pixel=%i), fired by %i track\n",iTOFpixel[ipart],
			pixelArray[iTOFpixel[ipart]-1].GetSector(),pixelArray[iTOFpixel[ipart]-1].GetPlate(),
			pixelArray[iTOFpixel[ipart]-1].GetStrip(),pixelArray[iTOFpixel[ipart]-1].GetPixel(),
			pixelArray[iTOFpixel[ipart]-1].GetTrack());   
	printf("     indexOfTestTrack=%i itest=%i weightTestTracksOutTof[4]=%f weightTestTracksOutTof[2]=%f weightTestTracksOutTof[1]=%f weightTestTracksOutTof[0]=%f\n",
			indexOfTestTrack,itest,weightTestTracksOutTof[3],weightTestTracksOutTof[2],
			weightTestTracksOutTof[1],weightTestTracksOutTof[0]);
	if(itest) {

	  printf("     take ipixel=%i (Sector=%i, Plate=%i, Strip=%i, Pixel=%i), (fired by %i track), match=%i\n",
			  wPixel,pixelArray[wPixel-1].GetSector(),pixelArray[wPixel-1].GetPlate(),pixelArray[wPixel-1].GetStrip(),
			  pixelArray[wPixel-1].GetPixel(),pixelArray[wPixel-1].GetTrack(),match);   
	}
      }
    }
    if(PRINT && itra<10 ) {

      if(itest) {
	cout << "      number of pixels with test tracks=" << itest << endl;
	for(Int_t i=0;i<itest;i++) {
	  cout << "      " << i+1 << "  tr.=" << ntest[i] << "  w=" << testWeight[i] << "  pix.= " << testPixel[i] << " (" << 
	    pixelArray[testPixel[i]-1].GetSector() << " " << " " << pixelArray[testPixel[i]-1].GetPlate() << " " << 
	    pixelArray[testPixel[i]-1].GetPixel() << " )" << "  l= " << testLength[i] << " sig=" << 
	    theta0*(testLength[i]-tpclength) << "  rho= " << testRho[i] << "  z= " << testZ[i] << endl;
	}
	cout << "      pixel=" << wPixel << "  state=" << istate << "  l=" << wLength << "  TOF=" << time << "  m=" << mass << "  match=" << match <<  endl;
	if(istate>0) cout << "      fired by track " << pixelArray[wPixel-1].GetTrack() << endl;
      }
    }
  }  //end of track
  

  if(itr) {
    printf(" %f probe tracks per 1 real track\n",snr/itr);   
    itrack=itr;
  }
  
  
  cout << ipixel << " - total number of TOF pixels after matching" << endl;
  w=iRmin;
  if(idelR!=0) {
    w /= idelR;
    printf(" %i tracks with delR, %f of them have delR>Rmin \n",idelR,w);
  }
  w=iRmin1;
  if(idelR1!=0) {
    w /= idelR1;
    printf(" %i tracks with delR1 (|z|<175), %f of them have delR>Rmin \n",idelR1,w);
  }
  w=iRmin2;
  if(idelR2!=0) {
    w /= idelR2;
    printf(" %i tracks with delR2 (|z|>175), %f of them have delR>Rmin \n",idelR2,w);
  }
  
  cout << " ********************  End of matching **********" << endl;
  delete [] ntest;
  delete [] testPixel;
  delete [] testLength;
  delete [] testRho;
  delete [] testZ;
  delete [] testWeight;
}

//____________________________________________________________________________
void AliTOFReconstructioner::FillNtuple(Int_t ntracks, AliTOFTrack* trackArray, AliTOFRecHit* hitArray, AliTOFPad* pixelArray, Int_t* iTOFpixel, Int_t* iparticle, Float_t* toftime, Int_t& ipixelLastEntry, Int_t itrack){
  
  // itrack : total number of TPC selected tracks
  // for the caller is ntotTPCtracks
  
  cout << " ********************  Start of searching non-matched fired pixels **********" << endl;
  const Int_t charge[48]={ 0, 1,-1, 0, 1,-1, 0, 1,-1, 0,
			   1,-1, 0, 1,-1, 0, 0, 0, 1, 0,
			   -1, 0,-1,-1, 0, 0,-1, 0, 1, 0,
			   1, 1, 0, 0, 1,-1, 0, 0, 1,-1,
			   1, 1,-1, 0, 1, 1, 2, 0};

  Int_t macthm1=0;
  Int_t macthm2=0;
  Int_t macthm3=0;
  Int_t macthm4=0;
  Int_t macth0=0;
  Int_t macth1=0;
  Int_t macth2=0;
  Int_t macth3=0;
  Int_t macth4=0;
  
  
  Float_t smat[9],smat0[9],smat1[9];
  for(Int_t i=0;i<9;i++) {
    smat[i]=0.;
    smat0[i]=0.;
    smat1[i]=0.;
  }
  
  Int_t nFiredPixelsNotMatchedWithTracks=0;
  Int_t istate;
  for (Int_t i=0; i<ipixelLastEntry; i++) {
    istate=pixelArray[i].GetState();
    if(istate==0) break;
    if(pixelArray[i].GetTrackMatched()==-1) nFiredPixelsNotMatchedWithTracks++;
  }
  printf("  %i fired pixels have not matched tracks\n",nFiredPixelsNotMatchedWithTracks);
  cout << " ********************  End of searching non-matched fired pixels **********" << endl;
  
  Int_t nTPCHitMissing=0;
  for(Int_t i=0; i<ipixelLastEntry; i++) {
    if(pixelArray[i].GetHit()>0) {
      if(hitArray[pixelArray[i].GetHit()-1].GetNoise()==0) {
	if(iparticle[pixelArray[i].GetTrack()]==0) nTPCHitMissing++;
      }
    }
  }
  printf("  %i pixels fired by track hit without a hit on the last layer of TPC\n",nTPCHitMissing);
  
  
  Int_t icharge=0;   // total number of charged particles
  Int_t iprim=0;     // number of primaries
  Int_t ipions=0;    // number of primary pions
  Int_t ikaons=0;    // number of primary kaons
  Int_t iprotons=0;  // number of primary protons
  Int_t ielectrons=0;// number of primary electrons
  Int_t imuons=0;    // number of primary muons
  Float_t particleTypeArray[6][5][2];
  
  for (Int_t index1=0;index1<6;index1++) {
    for (Int_t index2=0;index2<5;index2++) {
      for (Int_t index3=0;index3<2;index3++) {
	particleTypeArray[index1][index2][index3]=0.;
      }
    }
  }
  
  Int_t nTOFhitsWithNoTPCTracks=0; // to be moved later when used
  
  /*
  TObjArray *Particles = gAlice->Particles();
  Int_t numberOfParticles=Particles->GetEntries();
  cout << "numberOfParticles " << numberOfParticles << endl;
  // fpdbg
  if(numberOfParticles>fMaxAllTracks) numberOfParticles=fMaxAllTracks;
  */

  for (Int_t i=0; i<ntracks; i++) { // starting loop on all primaries charged particles for current event)

    /*
    cout << "particle " << i << endl;
    cout << "total " << numberOfParticles << endl;
    */
    TParticle *part = (TParticle *) gAlice->GetMCApp()->Particle(i);
    if(charge[PDGtoGeantCode(part->GetPdgCode())-1]) {
      icharge++;
      /*
      cout << "charged particles " << icharge << endl;
      */
      Int_t particleType=0;
      Int_t absPdgCode = TMath::Abs(part->GetPdgCode());
      switch (absPdgCode) {
      case 211:
	particleType=3;
	break ;
      case 321:
	particleType=2;
	break ;
      case 2212:
	particleType=1;
	break ;
      case 11:
	particleType=4;
	break ;
      case 13:
	particleType=5;
	break ;
      }
      
      if(part->GetFirstMother() < 0) {
	iprim++;
	switch (particleType) {
	case 1:
	  iprotons++;
	  break ;
	case 2:
	  ikaons++;
	  break ;
	case 3:
	  ipions++;
	  break ;
	case 4:
	  ielectrons++;
	  break ;
	case 5:
	  imuons++;
	  break ;
	}
      }
      
      Int_t match=0;
      Float_t wLength=-1.;
      Float_t time=-1.;
      Float_t mass=-1.;
      
      Int_t itr=iparticle[i]; // get the track number for the current charged particle
      
      if(iTOFpixel[i]>0 && itr==0) nTOFhitsWithNoTPCTracks++;
      
      if(itr) {
	match=trackArray[itr-1].GetMatching();
	//cout << "match " << match << endl;
	wLength=trackArray[itr-1].GetLength();
	//cout << "wLength " << wLength << endl;
	time=trackArray[itr-1].GetTof();
	mass=trackArray[itr-1].GetMassTOF();
	//cout << "mext " << mass << endl;
	//        if(PRINT && (i>789 && i<800) ) cout << i << " track:  l=" << wLength << "  TOF=" << time << "  m=" << mass << "  match=" << match <<  endl; 
	if(iTOFpixel[i]==0) {
	  smat0[match+4]++;
	  wLength=-wLength;
	}
      }
      Int_t ikparen=part->GetFirstMother();
      Int_t imam;
      if(ikparen<0) {
	imam=0;
      } else {
	imam=part->GetPdgCode();
      }
      
      Int_t evnumber=gAlice->GetEvNumber();
      if(match==-1) macthm1++;
      if(match==-2) macthm2++;
      if(match==-3) macthm3++;
      if(match==-4) macthm4++;
      if(match==0) macth0++;
      if(match==1) macth1++;
      if(match==2) macth2++;
      if(match==3) macth3++;
      if(match==4) macth4++;
      foutputntuple->Fill(evnumber,part->GetPdgCode(),imam,part->Vx(),part->Vy(),part->Vz(),part->Px(),part->Py(),part->Pz(),toftime[i],wLength,match,time,mass);
      
      
      
      // -----------------------------------------------------------
      // Filling 2 dimensional Histograms true time vs matched time
      // Filling 1 dimensional Histogram true time - matched time
      //
      // time              = time associated to the matched pad [ns]
      //                     it could be the average time of the cluster fired
      //
      // toftime[i]        = real time (including pulse height delays) [s]
      //
      //
      // if (time>=0) {
      // if (imam==0) TimeTrueMatched->Fill(time, toftime[i]*1E+09);
      // if (imam==0) DeltaTrueTimeMatched->Fill(time-toftime[i]*1E+09);
      // }
      //
      //---------------------------------------------------------------
      
      if(match==-4 || match>0) {
	Int_t matchW;
	matchW=match;
	if(match==-4) matchW=1;
	if(particleType) {
	  particleTypeArray[particleType-1][matchW-1][1]++;
	  particleTypeArray[5][matchW-1][1]++;
	  particleTypeArray[particleType-1][4][1]++;
	  particleTypeArray[5][4][1]++;
	  if(part->GetFirstMother() < 0) {
	    particleTypeArray[particleType-1][matchW-1][0]++;
	    particleTypeArray[5][matchW-1][0]++;
	    particleTypeArray[particleType-1][4][0]++;
	    particleTypeArray[5][4][0]++;
	    
	    // fill histos for QA
	    //if(particleType==3 && matchW==3) hPiWithTrueTime->Fill(sqrt((part->Px())*(part->Px())+(part->Py())*(part->Py())+(part->Pz())*(part->Pz())));
	    //if(particleType==2 && matchW==3) hKWithTrueTime->Fill(sqrt((part->Px())*(part->Px())+(part->Py())*(part->Py())+(part->Pz())*(part->Pz())));
	    //if(particleType==1 && matchW==3) hPWithTrueTime->Fill(sqrt((part->Px())*(part->Px())+(part->Py())*(part->Py())+(part->Pz())*(part->Pz())));
	    //
	    
	  } // close if(part->GetFirstMother() < 0)
	} // close if(particleType)
      } // close if(match==-4 || match>0)
    } // close if(charge[PDGtoGeantCode(part->GetPdgCode())-1])
  } // close for (Int_t i=0; i<ntracks; i++) {

  cout <<  " macthm1 " << macthm1 << endl;
  cout <<  " macthm2 " << macthm2 << endl;
  cout <<  " macthm3 " << macthm3 << endl;
  cout <<  " macthm4 " << macthm4 << endl;
  cout <<  " macth0 " << macth0 << endl;
  cout <<  " macth1 " << macth1 << endl;
  cout <<  " macth2 " << macth2 << endl;
  cout <<  " macth3 " << macth3 << endl;
  cout <<  " macth4 " << macth4 << endl;
  

  printf(" %i TOF hits have not TPC track\n",nTOFhitsWithNoTPCTracks);
  Int_t imatch=0;
  for(Int_t i=0;i<9;i++) {
    if(itrack) cout << "   " << smat[i]*100./itrack << " % of them (="<<smat[i]<<") have match=" << i-4 << "  " << smat0[i] << " have not TOF hits" << endl;
    if(i==0 || i>4) imatch += (Int_t) (smat[i]);
    
    //     cout << "   " << smat[i]*100./itrack << " % of them (="<<smat[i]<<") have match=" << i-4 << "  " << smat0[i] << " have not TOF hits" << "  " << smat1[i] << " have (r.p)<0 for first hit" << endl;
  }
  
  if(fdbg){
    /*
    cout << " nparticles = " << numberOfParticles << "  charged = " << icharge << "  prim.=" << iprim << endl;
    */
    cout << " nparticles = " << ntracks << "  charged = " << icharge << "  prim.=" << iprim << endl;
    cout << ipions << " - primary pions" << endl;
    cout << ikaons << " - primary kaons" << endl;
    cout << iprotons << " - primary protons" << endl;
    cout << ielectrons << " - primary electrons" << endl;
    cout << imuons << " - primary muons reached TPC" << endl;
    cout << " ********** " << imatch << " TPC tracks are matched with TOF pixels (incl.match=-4) **********" << endl;
  }
  
  /*
    Float_t PrimaryInBarrel[6],Acceptance[6];
    PrimaryInBarrel[0]=ipions;
    PrimaryInBarrel[1]=ikaons;
    PrimaryInBarrel[2]=iprotons;
    PrimaryInBarrel[3]=ielectrons;
    PrimaryInBarrel[4]=imuons;
    PrimaryInBarrel[5]=ipions+ikaons+iprotons+ielectrons+imuons;
    
    //   cout << "   TPC acceptance for the primary species: 1-p, 2-K, 3-pi, 4-e, 5-mu, 6-all" << endl; 
    for(Int_t i=0; i<6; i++) {
     Acceptance[i]=0.;
     if(PrimaryInBarrel[i]) Acceptance[i]=100.*PrimaryReachedTPC[i]/PrimaryInBarrel[i];
     //hTPCacceptance[i]->Fill(Acceptance[i]);
     //     printf(" species: %i    %f\n",i+1,Acceptance[i]);     
     }
     
     //   cout << "   TOF acceptance for the primary species: 1-p, 2-K, 3-pi, 4-e, 5-mu, 6-all" << endl; 
     for(Int_t i=0; i<6; i++) {
     Acceptance[i]=0.;
     if(PrimaryInBarrel[i]) Acceptance[i]=100.*PrimaryReachedTOF[i]/PrimaryInBarrel[i];
     //hTOFacceptance[i]->Fill(Acceptance[i]);
     //     printf(" species: %i    %f\n",i+1,Acceptance[i]);     
     }
     
   for (Int_t index1=0;index1<6;index1++) {
   for (Int_t index2=0;index2<4;index2++) {
   for (Int_t index3=0;index3<2;index3++) {
   if(particleTypeArray[index1][4][index3]) particleTypeArray[index1][index2][index3]=
						    100.*particleTypeArray[index1][index2][index3]/particleTypeArray[index1][4][index3]; 
						    }
     }
     }
     
   cout << "species: 1-p, 2-K, 3-pi, 4-e, 5-mu, 6-all" << endl; 
   cout << " matched pixels(%): 1-unfired 2-double 3-true 4-wrong 5-total number of tracks" << endl;
   
   cout << "  primary tracks:" << endl; 
   for (Int_t i=0;i<6;i++) {
     cout << i+1 << "  " << particleTypeArray[i][0][0] << "  " << particleTypeArray[i][1][0] << "  " << particleTypeArray[i][2][0] << "  " << particleTypeArray[i][3][0] << "  " << particleTypeArray[i][4][0] << endl; 
     }
     
     //   cout<<"      contam.for all prim.(%)="<<100*particleTypeArray[5][3][0]/(particleTypeArray[5][3][0]+particleTypeArray[5][2][0])<<endl;

     cout << "  all tracks:" << endl; 
     for (Int_t i=0;i<6;i++) {
     cout << i+1 << "  " << particleTypeArray[i][0][1] << "  " << particleTypeArray[i][1][1] << "  " << particleTypeArray[i][2][1] << "  " << particleTypeArray[i][3][1] << "  " << particleTypeArray[i][4][1] << endl; 
   } 
   
   //   cout<<"      contam.for all (%)="<<100*particleTypeArray[5][3][1]/(particleTypeArray[5][3][1]+particleTypeArray[5][2][1])<<endl;
   //  printf(" t34=%i, t32=%i, t44=%i, t43=%i, t42=%i\n",t34,t32,t44,t43,t42);
   //  printf(" m10=%f, m20=%f, m22=%f, m23=%f, m2state=%i\n",m10,m20,m22,m23,m2state);
  */
}
