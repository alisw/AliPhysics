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
// This is a TTask for reconstruction V2 in TOF
// Description of the algorithm
//-- Author: F. Pierella | pierella@bo.infn.it
//////////////////////////////////////////////////////////////////////////////


#include <Riostream.h>
#include <stdlib.h>

#include <TBenchmark.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TFolder.h>
#include <TGeant3.h>
#include <TParticle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTask.h>
#include <TTree.h>
#include <TVirtualMC.h>

#include "../TPC/AliTPCtrack.h"
#include "../TRD/AliTRDtrack.h"
#include "AliDetector.h"
#include "AliHeader.h"
#include "AliKalmanTrack.h"
#include "AliLoader.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliTOF.h"
#include "AliTOFGeometry.h"
#include "AliTOFDigitMap.h"
#include "AliTOFHitMap.h"
#include "AliTOFReconstructionerV2.h"
#include "AliTOFTrackV2.h"
#include "AliTOFdigit.h"
#include "AliTOFhitT0.h"
#include "AliMC.h"

ClassImp(AliTOFReconstructionerV2)

//____________________________________________________________________________ 
  AliTOFReconstructionerV2::AliTOFReconstructionerV2():TTask("AliTOFReconstructionerV2","") 
{
  //
  // std ctor
  // set all member vars to zero
  //
  fdbg      =0;
  fDigitsMap=0x0;
  fField    =0;
  fNDummyTracks=0;
  fScaleSigmaFactor=0.;
  fStep     =0; 
  fTOFDigits=0x0;
  fTOFTracks=0x0;
  fTOFDigitsFile    ="digits.root";
  fTPCBackTracksFile="AliTPCBackTracks.root";
  fKalmanTree       =0x0;
  fBranchWithTracks =0x0;
}
           
//____________________________________________________________________________ 
AliTOFReconstructionerV2::AliTOFReconstructionerV2(char* tpcBackTracks, char* tofDigits):TTask("AliTOFReconstructionerV2","") 
{
  //
  // par ctor
  // default arguments are specified only in
  // the header file
  // Parameters:
  // tpcBackTracks -> file name with backpropagated tracks in TPC
  // tofDigits     -> file name with TOF digits
  //

  fdbg      =0;
  fDigitsMap=0x0;
  fField    =0.2;   // default value 0.2 [Tesla]
  fNDummyTracks=20; // by default 20 test tracks used
  fScaleSigmaFactor=1.;
  fStep     =0.01;  // [cm] 
  fTOFDigits=0x0;
  fTOFTracks=0x0;
  fTOFDigitsFile    =tofDigits;
  fTPCBackTracksFile=tpcBackTracks;
  fKalmanTree       =0x0;
  fBranchWithTracks =0x0;

  // initialize the G3 geometry 
  gAlice->Init();
  gAlice->Print(); 

  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
}

//____________________________________________________________________________
AliTOFReconstructionerV2::AliTOFReconstructionerV2(const AliTOFReconstructionerV2 & /*rec*/)
:TTask("AliTOFReconstructionerV2","")
{
  //
  // Dummy copy constructor
  // required by coding conventions
  ;
}



//____________________________________________________________________________ 
  AliTOFReconstructionerV2::~AliTOFReconstructionerV2()
{
  //
  // dtor
  // some delete has to be moved
  // 

  if (fDigitsMap)
    {
      delete fDigitsMap;
      fDigitsMap = 0;
    }

  if (fTOFDigits)
    {
      delete fTOFDigits;
      fTOFDigits = 0;
    }

  if (fTOFTracks)
    {
      delete fTOFTracks;
      fTOFTracks = 0;
    }

  if (fTOFDigitsFile)
    {
      delete fTOFDigitsFile;
    }

  if (fTPCBackTracksFile)
    {
      delete fTPCBackTracksFile;
    }

  if (fKalmanTree)
    {
      delete fKalmanTree;
      fKalmanTree = 0;
    }

  if (fBranchWithTracks)
    {
      delete fBranchWithTracks;
      fBranchWithTracks = 0;
    }
}


//____________________________________________________________________________
void AliTOFReconstructionerV2::Exec(Option_t* option) 
{ 
  //
  // Description of the algorithm:
  //
  //
  //
  //
  //
  //

  // load TOF digits and fill the digit map
  Int_t tofDigitsLoading=LoadTOFDigits();

  // load back-propagated tracks in TPC
  Int_t tpcTracksLoading=LoadTPCTracks();

  if(tofDigitsLoading || tpcTracksLoading) {
    cout<<" Couldn't start reconstruction V2. Exit."<<endl;
    exit(1);
  }
  // or load TRD reconstructed tracks
  // Int_t trdTracksLoading=LoadTRDTracks();

  // create a TObjArray to store reconstructed tracks
  // and reject fake tracks
  const Int_t maxRecTracks = 10000; // max number of reconstructed tracks
                                    // per event 

  const Float_t stripRegionHeight = 15.3; // [cm] height in radial direction 
                                          // of the volume where strips are placed
  TObjArray trackArray(maxRecTracks);

  // buffer for output
  fTOFTracks= new TClonesArray("AliTOFTrackV2");
  // create a reference to fill the TClonesArray
  TClonesArray &aTOFTracks = *fTOFTracks;

  const Int_t maxIndex = 100000; // max number of primary tracks to be analysed 
  // the index of the rtIndex array is the track label
  // the content of the array is -1 for fake tracks
  // and the track index for good tracks
  Float_t dEdXarray[maxRecTracks];
  Int_t rtIndex[maxIndex];
  for(Int_t i = 0; i < maxIndex; i++) rtIndex[i] = -1;

  AliKalmanTrack::SetConvConst(100/0.299792458/fField);

  Int_t nRecTracks = (Int_t) fKalmanTree->GetEntries();
  cout<<"Found "<<nRecTracks<<" entries in the track tree "<<endl;

  // load the tracks into the array
  for (Int_t i=0; i<nRecTracks; i++) {
    AliTPCtrack *iotrack=new AliTPCtrack();
    fBranchWithTracks->SetAddress(&iotrack);
    fKalmanTree->GetEvent(i);
    trackArray.AddLast(iotrack);
    Int_t trackLabel = iotrack->GetLabel();
    dEdXarray[i]=iotrack->GetdEdx(); // usefull for PID
    // 
    // start filling the TClonesArray of AliTOFTrackV2
    Float_t trdXYZ[3]={0.,0.,0.};
    Float_t trdPxPyPz[3]={0.,0.,0.};

    // tpc outer wall positions
    Double_t xk=iotrack->GetX();
    // get the running coordinates in the lrf
    // and alpha
    Double_t x=xk;
    Double_t y=iotrack->GetY();
    Double_t z=iotrack->GetZ();
    Double_t alpha=iotrack->GetAlpha();
    GetGlobalXYZ(alpha, x, y, z);
    Float_t tpcXYZ[3]={x,y,z};

    // momentum at the end of TPC
    Float_t lambda=TMath::ATan(iotrack->GetTgl());
    Float_t invpt=TMath::Abs(iotrack->Get1Pt());
    Float_t pt=-99.;
    if (invpt) pt= 1./invpt; // pt
    Float_t tpcmom=1./(invpt*TMath::Cos(lambda));
    Float_t pz=tpcmom*TMath::Sin(lambda);
    Float_t tpcPtPz[2]={pt,pz};

    Int_t matchingStatus=-1;
    if(trackLabel < 0) matchingStatus=0;
    new(aTOFTracks[i]) AliTOFTrackV2(trackLabel,matchingStatus,tpcmom,dEdXarray[i],tpcXYZ,tpcPtPz,trdXYZ,trdPxPyPz);
    //    printf("rt with %d clusters and label %d \n",
    //     iotrack->GetNumberOfClusters(), trackLabel);

    if(trackLabel < 0) continue;
    if(trackLabel >= maxIndex) continue;
    rtIndex[trackLabel] = i;
    delete iotrack;
  }

  if(strstr(option,"MC")) Comparison(rtIndex);

  // start loop on tracks
  // and backpropagate them from TPC to TOF
  // backpropagation is performed only
  // for good tracks (fake tracks rejected)
  AliTPCtrack *rt;

  for (Int_t i=0; i<nRecTracks; i++) {

    //******* tracking:  extract track coordinates, momentum, etc. 
    rt = (AliTPCtrack*) trackArray.UncheckedAt(i);
    // track length to be implemented
    // Double_t tr_length = rt->GetLength();

    Int_t tpcTrackLabel=rt->GetLabel();
    // reject fake tracks
    //if(tpcTrackLabel< 0) continue;

    // starting backpropagation to TOF
    // here we assume to have backpropagated tracks in TPC
    // starting point xk=246.055
    // starting back propagation
    // outer wall of the TPC
    Int_t iOuterTPCWall=rt->PropagateTo(261.53,40.,0.06124);
    // frame with air just to the beginning of the TOF
    Int_t iFrameWithAir=rt->PropagateTo(370.,36.66,1.2931e-3);
    // trough the wall of the TOF plate
    // thickness has changed according to the
    // last value
    Int_t iTOFWall=rt->PropagateTo(370.11,24.01,2.7);
    
    /*
      // outer wall of the TPC
      Int_t iOuterTPCWall=rt->PropagateTo(261.53,40.,0.06124);
      // air in between TPC and TRD
      Int_t iFrameWithAir=rt->PropagateTo(294.5,36.66,1.2931e-3);
      // TRD itself
      // mean density for the TRD calculated from
      // TRD Technical Design Report
      // page 11  -> global thickness
      // page 23  -> different main layers thickness (Radiator Air/ Drift Chamber Gas /G10)
      // page 139 -> material budget and radiation lengths
      Int_t iTRD=rt->PropagateTo(369.1,171.7,0.33);
      // air in between TRD and TOF
      Int_t iFrameWithAirbis=rt->PropagateTo(370.,36.66,1.2931e-3);
      // trough the wall of the TOF plate
      Int_t iTOFWall=rt->PropagateTo(370.11,24.01,2.7);
     */
 
    // select only cases when
    // backpropagation succeded
    // and particle is in the geometrical TOF acceptance along Z
    // |z|<380 [cm]  
    AliTOFTrackV2* oTOFtracks=(AliTOFTrackV2*)fTOFTracks->UncheckedAt(i);
    Bool_t outOfZacceptance=(rt->GetZ()<=380.);
    if(outOfZacceptance) oTOFtracks->SetMatchingStatus(-2);

    if(iOuterTPCWall==1 && iFrameWithAir==1 && iTOFWall==1 && (!outOfZacceptance)){
      Double_t cc[15];
      // get sigmaY and sigmaZ
      rt->GetExternalCovariance(cc);
      //Double_t sigmaY =TMath::Sqrt(cc[0]); // [cm]
      //Double_t sigmaZ =TMath::Sqrt(cc[2]); // [cm]

      // arrays used by the DigitFinder
      Int_t nSlot=1;
      TArrayI *secArray= new TArrayI(nSlot);
      TArrayI *plaArray= new TArrayI(nSlot);
      TArrayI *strArray= new TArrayI(nSlot);
      TArrayI *pdzArray= new TArrayI(nSlot);
      TArrayI *pdxArray= new TArrayI(nSlot);

      // make fNDummyTracks clones of the current backpropagated track
      // make a copy of the current track
      // smear according to the backpropagated area
      for (Int_t j=0; j<fNDummyTracks; i++) {
	AliTPCtrack *dummyrt=new AliTPCtrack(*rt);
	// get ylrf and zlrf
	//Double_t ylrf= dummyrt->GetY();  // P0
	//Double_t zlrf= dummyrt->GetZ();  // P1

	// smear according to sigmaY and sigmaZ
	//Double_t ylrfNew=gRandom->Gaus(ylrf,fScaleSigmaFactor*sigmaY);
	//Double_t zlrfNew=gRandom->Gaus(zlrf,fScaleSigmaFactor*sigmaZ);

	// set Y and Z accordingly
	// setter to be added in the class AliTPCtrack
	// here I have to modify the AliTPCtrack class
	// adding the setters for Y and Z
	//dummyrt->SetY(ylrfNew);
	//dummyrt->SetZ(zlrfNew);

	// start fine-backpropagation inside the TOF
	Bool_t padNotFound =kTRUE;
	Bool_t isInStripsRegion=kTRUE;
	Double_t xk=dummyrt->GetX();

	while (padNotFound && isInStripsRegion){
	  xk+=fStep;
	  // here we assume a frame with air
	  dummyrt->PropagateTo(xk,36.66,1.2931e-3);
	  // get the running coordinates in the lrf
	  // and alpha
	  Double_t x=xk;
	  Double_t y=dummyrt->GetY();
	  Double_t z=dummyrt->GetZ();
	  Double_t alpha=dummyrt->GetAlpha();
	  GetGlobalXYZ(alpha, x, y, z);
	  
	  // check if the point falls into a pad
	  // using the G3 geometry
	  Int_t* volumeID = new Int_t[AliTOFGeometry::MaxTOFTree()];
	  // volumeID[0] -> TOF Sector range [1-18]
	  // volumeID[1] -> TOF Plate  range [1- 5]
	  // volumeID[2] -> TOF Strip  max range [1-20]
	  // volumeID[3] -> TOF Pad along Z range [1- 2]
	  // volumeID[4] -> TOF Pad along X range [1-48]
	  
	  Float_t zInPadFrame=0.;
	  Float_t xInPadFrame=0.;
	  IsInsideThePad((Float_t)x,(Float_t)y,(Float_t)z, volumeID, zInPadFrame, xInPadFrame);
	  // adding protection versus wrong volume numbering
	  // to be released in the next release after debugging
	  if(volumeID[4]){
	    padNotFound=kFALSE;
	    nSlot++;
	    secArray->Set(nSlot-1);
	    plaArray->Set(nSlot-1);
	    strArray->Set(nSlot-1);
	    pdzArray->Set(nSlot-1);
	    pdxArray->Set(nSlot-1);
	    
	    (*secArray)[nSlot-1]=volumeID[0];
	    (*plaArray)[nSlot-1]=volumeID[1];
	    (*strArray)[nSlot-1]=volumeID[2];
	    (*pdzArray)[nSlot-1]=volumeID[3];
	    (*pdxArray)[nSlot-1]=volumeID[4];
	    
	  } // track falls into a pad volume
	  
	  delete [] volumeID;
	  
	  // check on xk to stop the fine-propagation
	  if(xk>=(370.+stripRegionHeight)) isInStripsRegion=kFALSE;
	  
	} // close the while for fine-propagation
	
	delete dummyrt;
      } // end loop on test tracks
      
      // start TOF digit finder
      Int_t assignedVol[5]={0,0,0,0,0};
      Float_t tdc=-1.;
      Int_t* digitTrackArray=0x0;
      Bool_t assignedDigit=DigitFinder(secArray, plaArray, strArray, pdzArray, pdxArray, assignedVol, digitTrackArray, tdc);


      if(assignedDigit){
	// fill the tree for tracks with time of flight
	// tof is given in tdc bin
	// conversion to [ns]
	Float_t binWidth=50.; // [ps]
	Float_t timeOfFlight=tdc*binWidth/1000.;

	// only the first track number contributing
	// to the assigned digit
	Int_t tofDigitTrackLabel=digitTrackArray[0];

	// matching status for the current track
	Int_t matching=4;
	if(tpcTrackLabel==digitTrackArray[0] || tpcTrackLabel==digitTrackArray[1] || tpcTrackLabel==digitTrackArray[0]) matching=3;
	oTOFtracks->UpdateTrack(tofDigitTrackLabel, matching, timeOfFlight);
      } else {
	// fill the TClonesArray for tracks with no time of flight
	Int_t matching=1;
	oTOFtracks->SetMatchingStatus(matching);
      }
      
      // delete used memory for tmp arrays used by DigitFinder
      delete secArray;
      delete plaArray;
      delete strArray;
      delete pdzArray;
      delete pdxArray;
      
    } // close the if for succeded backpropagation in TOF acceptance along z
    
  } // end loop on reconstructed tracks
  
  // free used memory for digitmap
  delete fDigitsMap;

  // save array with TOF tracks
  Int_t output=SaveTracks();
  if(output) cout << "Error writing TOF tracks " << endl;
}


//__________________________________________________________________
void AliTOFReconstructionerV2::Init(Option_t* /*opt*/)
{
  //
  // Initialize the AliTOFReconstructionerV2
  //
  //
}


//__________________________________________________________________
Int_t AliTOFReconstructionerV2::LoadTPCTracks()
{
  //
  // Connect the tree and the branch
  // with reconstructed tracks
  // backpropagated in the TPC

  gBenchmark->Start("LoadTPCTracks");

  TFile *kalFile    = TFile::Open(fTPCBackTracksFile.Data());
  if (!kalFile->IsOpen()) {cerr<<"Can't open AliTPCBackTracks.root !\n"; return 3;}

  // tracks from Kalman
  Int_t event=0;
  char treename[100]; sprintf(treename,"TreeT_TPCb_%d",event);
  fKalmanTree=(TTree*)kalFile->Get(treename);
  if (!fKalmanTree) {cerr<<"Can't get a tree with TPC back tracks !\n"; return 4;}

  // otherwise you get always 0 for 1/pt
  AliKalmanTrack::SetConvConst(100/0.299792458/fField);

  fBranchWithTracks=fKalmanTree->GetBranch("tracks");
  Int_t kalEntries =(Int_t)fKalmanTree->GetEntries();
  cout<<"Number of loaded Tracks :"<< kalEntries <<endl;

  gBenchmark->Stop("LoadTPCTracks");
  gBenchmark->Show("LoadTPCTracks");   
  return 0;
}

//__________________________________________________________________
Int_t AliTOFReconstructionerV2::LoadTRDTracks()
{
  //
  // Connect the tree and the branch
  // with reconstructed tracks in TRD

  gBenchmark->Start("LoadTRDTracks");

  Int_t nEvent = 0;
  const Int_t nPrimaries = 84210/16;
  const Int_t maxIndex = nPrimaries;
  Int_t rtIndex[maxIndex];

  TFile *tf=TFile::Open("AliTRDtracks.root");

  if (!tf->IsOpen()) {cerr<<"Can't open AliTRDtracks.root !\n"; return 3;}
  TObjArray tarray(2000);
  char   tname[100];
  sprintf(tname,"TRDb_%d",nEvent);     
  TTree *tracktree=(TTree*)tf->Get(tname);

  TBranch *tbranch=tracktree->GetBranch("tracks");

  Int_t nRecTracks = (Int_t) tracktree->GetEntries();
  cerr<<"Found "<<nRecTracks<<" entries in the track tree"<<endl;

  for (Int_t i=0; i<nRecTracks; i++) {
    AliTRDtrack *iotrack=new AliTRDtrack();
    tbranch->SetAddress(&iotrack);
    tracktree->GetEvent(i);
    tarray.AddLast(iotrack);
    Int_t trackLabel = iotrack->GetLabel();

    //    printf("rt with %d clusters and label %d \n",
    //     iotrack->GetNumberOfClusters(), trackLabel);

    if(trackLabel < 0) continue;
    if(trackLabel >= maxIndex) continue;
    rtIndex[trackLabel] = i;
  }
  tf->Close();                 
  gBenchmark->Stop("LoadTRDTracks");
  gBenchmark->Show("LoadTRDTracks");   
  return 0;

}

//__________________________________________________________________
Int_t AliTOFReconstructionerV2::LoadTOFDigits()
{
  //
  // Connect the TClonesArray with TOF
  // digits and fill the digit map
  // used by the DigitFinder

  Int_t rc=0;

  gBenchmark->Start("LoadTOFDigits");
  
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(fTOFDigitsFile.Data());
  if(file){
    cout<<"headerFile already open \n";
  }
  else {
    if(!file)file=TFile::Open(fTOFDigitsFile.Data());
  }
  
  // Get AliRun object from file
  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
  }
  
  Int_t iEvNum = 0;
  if (iEvNum == 0) iEvNum = (Int_t) gAlice->TreeE()->GetEntries();


  AliTOFdigit *tofdigit;

  AliTOF * tof = (AliTOF *) gAlice->GetDetector("TOF") ;

  if (!tof) {
    cout << "<LoadTOFDigits> No TOF detector found" << endl;
    rc = 2;
    return rc;
  }

  for (Int_t ievent = 0; ievent < iEvNum; ievent++) {

    gAlice->GetEvent(ievent) ;
    if(gAlice->TreeD()==0) {
      cout << "<LoadTOFDigits> No  TreeD found" << endl;
      rc = 4;
      return rc;
    }
    

    Int_t ndig;
    gAlice->ResetDigits();
    gAlice->TreeD()->GetEvent(ievent);
    fTOFDigits   = tof->Digits();
    
    ndig=fTOFDigits->GetEntries();

    // create the digit map
    fDigitsMap = new AliTOFDigitMap(fTOFDigits);

    
    cout << "<LoadTOFDigits> found " << ndig
	 << " TOF digits for event " << ievent << endl;
    
    for (Int_t k=0; k<ndig; k++) {
      tofdigit= (AliTOFdigit*) fTOFDigits->UncheckedAt(k);
      Float_t tdc=tofdigit->GetTdc();
      // adc value can be used for weighting
      //Float_t adc=tofdigit->GetAdc();

      // TOF digit volumes
      Int_t    vol[5];       // location for a digit
      Int_t sector    = tofdigit->GetSector(); // range [1-18]
      Int_t plate     = tofdigit->GetPlate();  // range [1- 5]
      Int_t strip     = tofdigit->GetStrip();  // range [1-20]
      Int_t padx      = tofdigit->GetPadx();   // range [1-48]
      Int_t padz      = tofdigit->GetPadz();   // range [1- 2]

      vol[0] = sector;
      vol[1] = plate;
      vol[2] = strip;
      vol[3] = padx;
      vol[4] = padz;

      // QA
      Bool_t isDigitBad = (sector<1 || sector>18 || plate<1 || plate >5 || padz<1 || padz>2 || padx<1 || padx>48);

      if (isDigitBad) {
	cout << "<LoadTOFDigits>  strange digit found" << endl;
	rc = 3;
	return rc;
      } // if (isDigitBad)

      // Fill the digit map checking if the location is already used
      // in this case we take the earliest signal
      if (fDigitsMap->TestHit(vol) != kEmpty) {
	// start comparison in between the 2 digit
	AliTOFdigit *dig = static_cast<AliTOFdigit*>(fDigitsMap->GetHit(vol));
	if(tdc < (dig->GetTdc())) fDigitsMap->SetHit(vol,k);
	// we can add also the check on adc value
	// by selecting the largest adc value
      } else {
	fDigitsMap->SetHit(vol,k);
      }
      // to be added protection versus 2-digit on the same pad
      // we have to have memory also of the second digit
      
    } // for (k=0; k<ndig; k++)
  
  } // end loop on events

  gBenchmark->Stop("LoadTOFDigits");
  gBenchmark->Show("LoadTOFDigits");   
  return rc;
}

//__________________________________________________________________
void AliTOFReconstructionerV2::IsInsideThePad(Float_t x, Float_t y, Float_t z, Int_t *nGeom, Float_t& zPad, Float_t& xPad) 
{
  //   input: x,y,z - coordinates of a point in the mrf [cm]
  //   output: array  nGeom[]
  //          nGeom[0] - the TOF sector number, 1,2,...,18 along azimuthal direction starting from -90 deg.
  //          nGeom[1] - the TOF module number, 1,2,3,4,5=C,B,A,B,C along z-direction
  //          nGeom[2] - the TOF strip  number, 1,2,... along z-direction
  //          nGeom[3] - the TOF padz  number,  1,2=NPZ across a strip
  //          nGeom[4] - the TOF padx  number,  1,2,...,48=NPX along a strip
  //          zPad, xPad - coordinates of the hit in the pad frame
  //  numbering is adopted for the version 3.08 of AliRoot
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
  
  TGeant3 * fG3Geom = (TGeant3*) gMC;

  fG3Geom->Gmedia(xTOF, numed);
  gcvolu=fG3Geom->Gcvolu();
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
    zPad=gcvolu->glx[2];
    xPad=gcvolu->glx[0];
  }
  
}

//__________________________________________________________________
void AliTOFReconstructionerV2::GetGlobalXYZ(Double_t alpha, Double_t& x, Double_t& y, Double_t& /*z*/)
{
  //
  // return the current running coordinates of 
  // the track in the global reference frame
  // x, y and z have to initialized to the
  // local frame coordinates by the caller
  // alpha is the alpha coordinate in the TPC Kalman
  // reference frame 

  // it take into account differences in between
  // TPC and TRD local coordinates frames
  if      (alpha < -TMath::Pi()) alpha += 2*TMath::Pi();
  else if (alpha >= TMath::Pi()) alpha -= 2*TMath::Pi();

  // conversion
  Double_t tmp=x*TMath::Cos(alpha) - y*TMath::Sin(alpha);
  y=x*TMath::Sin(alpha) + y*TMath::Cos(alpha);
  x=tmp;            
}

//__________________________________________________________________
Bool_t AliTOFReconstructionerV2::DigitFinder(TArrayI *secArray, TArrayI *plaArray, TArrayI *strArray, TArrayI *pdzArray, TArrayI *pdxArray, Int_t* assignedVol, Int_t* digitTrackArray, Float_t& tdc)
{
  //
  // Description
  // input: arrays with sectors, plates, strips, padz, padx
  // found during fine-propagation of probe tracks
  //
  // output kFALSE if signal is not found
  //        kTRUE  if signal is found
  // in this case the assignedVol array contains the digit volume numbers
  // and digitTrackArray the track numbers (max 3) contributing to the
  // digit

  Bool_t dummy=kFALSE;
  Int_t nFilledSlot=secArray->GetSize();

  // start loop
  Float_t maxWeight=-1.;
  Int_t   indexOfMaxWeight=-1;
  for (Int_t i = 0; i < nFilledSlot; i++) {
    Int_t    vol[5];       // location for a digit
    vol[0] = (*secArray)[i];
    vol[1] = (*plaArray)[i];
    vol[2] = (*strArray)[i];
    vol[3] = (*pdxArray)[i];
    vol[4] = (*pdzArray)[i];

    // check for digit in the current location
    if (fDigitsMap->TestHit(vol) != kEmpty) {

      AliTOFdigit *dig = static_cast<AliTOFdigit*>(fDigitsMap->GetHit(vol));
      Float_t adcWeight=dig->GetAdc();
      if(adcWeight > maxWeight){
	maxWeight=adcWeight;
	indexOfMaxWeight=i;
	tdc=dig->GetTdc();
	digitTrackArray=dig->GetTracks();
 
      } // if(adcWeight > maxWeight)
    } // close if (fDigitsMap->TestHit(vol) != kEmpty)

  } // end loop

  if(indexOfMaxWeight!=-1){
    assignedVol[0]=(*secArray)[indexOfMaxWeight];
    assignedVol[1]=(*plaArray)[indexOfMaxWeight];
    assignedVol[2]=(*strArray)[indexOfMaxWeight];
    assignedVol[3]=(*pdxArray)[indexOfMaxWeight];
    assignedVol[4]=(*pdzArray)[indexOfMaxWeight];
    dummy=kTRUE;
  }

  return dummy;
}

//__________________________________________________________________
Int_t AliTOFReconstructionerV2::SaveTracks(const Char_t *outname, const Int_t split)
{
  //
  // save reconstructed tracks into 
  // outname file
  //
  TDirectory *savedir=gDirectory;
  const Char_t *name="Writing Output";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);
  
  TFile *out=TFile::Open(outname,"RECREATE");
  if (!out->IsOpen()) {
    cerr<<"AliTOFReconstructionerV2::SaveTracks(): ";
    cerr<<"file for TOF tracks is not open !\n";
     return 2;
  }
  
  out->cd();
  TTree T("T","tree with TOF tracks");
  T.Branch("tracks",&fTOFTracks,256000,split);


  T.Fill();
  T.Write();
  savedir->cd();
  out->Close();
  gBenchmark->Stop(name);
  gBenchmark->Show(name);
  return 0;
}

//__________________________________________________________________
void AliTOFReconstructionerV2::Comparison(Int_t* rtIndex)
{
  //
  // perform MC comparison
  // used also for track length
  // for the time being
  // Connect the AliRoot file containing Geometry, Kine, Hits, and Digits


  Int_t nEvent = 0;
  Char_t *datafile = "galice.root";
  
  AliRunLoader *rl = AliRunLoader::Open(datafile);
  if (rl == 0x0)
   {
     Error("Exec","Can not open session for file %s",datafile);
     return;
   }
  // Get AliRun object from file or create it if not on file
  rl->LoadgAlice();
  gAlice = rl->GetAliRun();
  if (gAlice)
    cout << "AliRun object found on file" << endl;
  else
    gAlice = new AliRun("gAlice","Alice test program");


  AliTOF* TOF = (AliTOF *) gAlice->GetDetector ("TOF");

  
  
  // Import the Trees for the event nEvent in the file
  rl->GetEvent(nEvent);
  const Int_t nparticles = rl->GetNumberOfEvents();
  if (nparticles <= 0) return;

  AliLoader* tofloader = rl->GetLoader("TOFLoader");
  if (tofloader == 0x0)
   {
    Error("AliTOFReconstructioner","Can not get TOF Loader from Run Loader.");
    delete rl;
    return;
   }

  // Get pointers to Alice detectors and Hits containers
  tofloader->LoadHits();
  Int_t ntracks    = (Int_t) tofloader->TreeH()->GetEntries();     
  TOF->SetTreeAddress();
  // Start loop on tracks in the hits containers
  for (Int_t track=0; track < ntracks; track++) {
    
    if(TOF) {
      for(AliTOFhitT0* tofHit = (AliTOFhitT0*)TOF->FirstHit(track); 
	  tofHit; 
	  tofHit=(AliTOFhitT0*)TOF->NextHit()) {
	
         Int_t ipart    = tofHit->GetTrack();
         if(ipart >= 80000) continue;
         if(rtIndex[ipart] < 0) continue; 

	 TParticle *part = gAlice->GetMCApp()->Particle(ipart);
	 
	 // get first the pdg code
	 Int_t pdgCode=part->GetPdgCode();
	 
	 // then track length
	 Float_t trackLength=tofHit->GetLen(); // [cm]

	 // update the tof TClonesArray with TOF tracks
	 AliTOFTrackV2* oTOFtracks=(AliTOFTrackV2*)fTOFTracks->UncheckedAt(rtIndex[ipart]);
	 oTOFtracks->UpdateTrack(pdgCode,trackLength);

      } // loop on hits connected to the current track
    } // if(TOF)
  } // end loop on primary tracks
}

//__________________________________________________________________
Bool_t AliTOFReconstructionerV2::operator==(const AliTOFReconstructionerV2 & tofrecv2)const
{
  // Equal operator.
  // Reconstructioner are equal if their fField, fNDummyTracks, fScaleSigmaFactor and fStep are equal
 
  if( (fField==tofrecv2.fField)&&(fNDummyTracks==tofrecv2.fNDummyTracks)&&(fScaleSigmaFactor==tofrecv2.fScaleSigmaFactor)&&(fStep==tofrecv2.fStep))
    return kTRUE ;
  else
    return kFALSE ;
}
