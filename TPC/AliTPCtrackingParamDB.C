/****************************************************************************
 * This macro is used to create a DataBase for the TPC tracking             *
 * parameterization.                                                        * 
 * 1) the function CreateAllGeantTracks gives all tracks at the 1st hit of  *
 *    the TPC                                                               *
 * 2) the function TrackCompare compares them with track found by the       *
 *    Kalman filter for the same event and computes efficiency and          *
 *    resolution on the track parameters for the Kalman filter.             *
 * 3) the function BuildDataBase calls many functions of AliTPCtrackerParam:*
 *    - merge results from TrackCompare for many events and compute         *
 *      average efficiency.                                                 *
 *    - analyze the pulls of the covariance matrix                          *
 *    - analyze the dE/dx                                                   *
 *    - regularize the covariance matrix as a function of the momentum      *
 *    - write all the informations and the trees with regularized cov.      *
 *      matrices in the DataBase file.                                      *
 *                                                                          *
 *  Origin: A.Dainese, Padova, andrea.dainese@pd,infn.it                    * 
 *                                                                          *
 ****************************************************************************/

#ifndef __CINT__
#include "Riostream.h"
#include <TFile.h>
#include <TStopwatch.h>
#include <TObject.h>
#include "alles.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliMagF.h"
#include "AliModule.h"
#include "AliArrayI.h"
#include "AliDigits.h"
#include "AliITS.h"
#include "AliTPC.h"
#include "AliITSgeom.h"
#include "AliITSRecPoint.h"
#include "AliITSclusterV2.h"
#include "AliITSsimulationFastPoints.h"
#include "AliITStrackerV2.h"
#include "AliKalmanTrack.h"
#include "AliTPCtrackerParam.h"
#endif

void TPCParamTracks(const Char_t *galice,const Char_t *outname,const Int_t coll,const Double_t Bfield,Int_t n);
void CreateAllGeantTracks(const Char_t *galice, const Char_t *outname,const Int_t coll,const Double_t Bfield,Int_t n);
void TrackCompare(const Int_t coll,const Double_t Bfield,Int_t n);
void BuildDataBase(const Int_t coll,const Double_t Bfield);

void AliTPCtrackingParamDB(Int_t n=1) {

  const Char_t *name=" AliTPCtrackingParamTest";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);


  const Char_t *TPCtrkName="AliTPCtracksParam.root";
  const Char_t *TPCgeatrkName="AliTPCtracksGeant.root";
  const Char_t *galiceName="galice.root";

  // set here the code for the type of collision (needed for TPC tracking
  // parameterization). available collisions:
  //
  // coll = 0 ->   PbPb6000 (HIJING with b<2fm) 
  const Int_t    collcode = 0;  
  // set here the value of the magnetic field
  const Double_t BfieldValue = 0.4;

  AliKalmanTrack::SetConvConst(100/0.299792458/BfieldValue);

  // ********** Build TPC tracks with parameterization *********** // 
  TPCParamTracks(galiceName,TPCtrkName,collcode,BfieldValue,n);
 

  // ********** Build All True TPC tracks *********** //
  CreateAllGeantTracks(galiceName,TPCgeatrkName,collcode,BfieldValue,n);
    
  // ********** Compare Kalman tracks with tracks at 1st hit *********** //
  TrackCompare(collcode,BfieldValue);
  
   
  // ********** Merge files, compute pulls, and dEdx *********** //
  BuildDataBase(collcode,BfieldValue);
  



  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return;
}

//_____________________________________________________________________________
void TPCParamTracks(const Char_t *galice, const Char_t *outname,
		    const Int_t coll,const Double_t Bfield,Int_t n) {
//
// Ordinary TPC tracking parameterization
//
  cerr<<"\n*******************************************************************\n";

  const Char_t *name="TPCParamTracks";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);

  TFile *outfile=TFile::Open(outname,"recreate");
  TFile *infile =TFile::Open(galice);

  AliTPCtrackerParam tracker(coll,Bfield,n);
  tracker.BuildTPCtracks(infile,outfile);

  delete gAlice; gAlice=0;

  infile->Close();
  outfile->Close();

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return;
}
//_____________________________________________________________________________
void CreateAllGeantTracks(const Char_t *galice, const Char_t *outname,
			  const Int_t coll,const Double_t Bfield,Int_t n) {
//
// Get all tracks at TPC 1st hit w/o selection and smearing
//
  cerr<<"\n*******************************************************************\n";

  const Char_t *name="CreateAllGeantTracks";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);

  TFile *outfile=TFile::Open(outname,"recreate");
  TFile *infile =TFile::Open(galice);

  AliTPCtrackerParam tracker(coll,Bfield,n);
  tracker.AllGeantTracks(); // this is to switch-off selection and smearing
  tracker.BuildTPCtracks(infile,outfile);

  delete gAlice; gAlice=0;

  infile->Close();
  outfile->Close();

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return;
}
//_____________________________________________________________________________
void TrackCompare(const Int_t coll,const Double_t Bfield,Int_t n) {
//
// Compare Kalman tracks with tracks at TPC 1st hit
//
  cerr<<"\n*******************************************************************\n";

  const Char_t *name="TrackCompare";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);

  AliTPCtrackerParam tracker(coll,Bfield,n);
  tracker.CompareTPCtracks();

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return;
}
//_____________________________________________________________________________
void BuildDataBase(const Int_t coll,const Double_t Bfield) {
//
//
//
  cerr<<"\n*******************************************************************\n";

  AliTPCtrackerParam tracker(coll,Bfield);

  // Merge files with cov. matrix and compute average efficiencies
  cerr<<"\n   --- Merging Events ---\n\n";
  tracker.MergeEvents(1,5);
 
  // Compute the pulls for pions, kaons, electrons
  cerr<<"\n   --- Analyzing Pulls ---\n\n";
  tracker.AnalyzePulls("pulls.root");

  // Draw pulls and efficiencies  
  tracker.DrawPulls("CovMatrixDB_PbPb6000_B0.4T.root",211,0);
  tracker.DrawEffs("CovMatrixDB_PbPb6000_B0.4T.root",13);

  // Regularize the covariance matrix
  tracker.RegularizeCovMatrix("regPi.root",211);

  // Analyze the dE/dx
  tracker.AnalyzedEdx("dEdxPi.root",211);


  // Put everything together and create the DB file
  tracker.MakeDataBase();

  return;
}







