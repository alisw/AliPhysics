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

/*
$Log$
Revision 1.39  2001/03/12 17:45:32  hristov
Changes needed on Sun with CC 5.0

Revision 1.38  2001/03/07 14:04:51  barbera
Some vector dimensions increased to cope with full events

Revision 1.37  2001/03/07 12:36:35  barbera
A change added in the tracking part to manage delta rays

Revision 1.36  2001/03/02 19:44:11  barbera
 modified to taking into account new version tracking v1

Revision 1.35  2001/02/28 18:16:46  mariana
Make the code compatible with the new AliRun

Revision 1.34  2001/02/11 15:51:39  mariana
Set protection in MakeBranch

Revision 1.33  2001/02/10 22:26:39  mariana
Move the initialization of the containers for raw clusters in MakeTreeC()

Revision 1.32  2001/02/08 23:55:31  nilsen
Removed fMajor/MinorVersion variables in favor of variables in derived classes.
Set arrays char *det[3] = {"SPD","SDD","SSD"} as const.

Revision 1.31  2001/02/02 23:57:28  nilsen
Added include file that are no londer included in AliITSgeom.h

Revision 1.30  2001/01/30 09:23:13  hristov
Streamers removed (R.Brun)

Revision 1.29  2001/01/26 20:01:09  hristov
Major upgrade of AliRoot code

Revision 1.28  2000/12/18 14:02:00  barbera
new version of the ITS tracking to take into account the new TPC track parametrization

Revision 1.27  2000/12/08 13:49:27  barbera
Hidden declaration in a for loop removed to be compliant with HP-UX compiler

Revision 1.26  2000/11/27 13:12:13  barbera
New version containing the files for tracking

Revision 1.25  2000/11/12 22:38:05  barbera
Added header file for the SPD Bari model

Revision 1.24  2000/10/09 22:18:12  barbera
Bug fixes from MAriana to le AliITStest.C run correctly

Revision 1.23  2000/10/05 20:47:42  nilsen
fixed dependencies of include files. Tryed but failed to get a root automaticly
generates streamer function to work. Modified SetDefaults.

Revision 1.9.2.15  2000/10/04 16:56:40  nilsen
Needed to include stdlib.h

=======
Revision 1.22  2000/10/04 19:45:52  barbera
Corrected by F. Carminati for v3.04

Revision 1.21  2000/10/02 21:28:08  fca
Removal of useless dependecies via forward declarations

Revision 1.20  2000/10/02 16:31:39  barbera
General code clean-up

Revision 1.9.2.14  2000/10/02 15:43:51  barbera
General code clean-up (e.g., printf -> cout)

Revision 1.19  2000/09/22 12:13:25  nilsen
Patches and updates for fixes to this and other routines.

Revision 1.18  2000/07/12 05:32:20  fca
Correcting several syntax problem with static members

Revision 1.17  2000/07/10 16:07:18  fca
Release version of ITS code

Revision 1.9.2.3  2000/02/02 13:42:09  barbera
fixed AliITS.cxx for new AliRun structure. Added ITS hits list to list of hits which will have their track numbers updated

Revision 1.9.2.2  2000/01/23 03:03:13  nilsen
//fixed FillModule. Removed fi(fabs(xl)<dx....

Revision 1.9.2.1  2000/01/12 19:03:32  nilsen
This is the version of the files after the merging done in December 1999.
See the ReadMe110100.txt file for details

Revision 1.9  1999/11/14 14:33:25  fca
Correct problems with distructors and pointers, thanks to I.Hrivnacova

Revision 1.8  1999/09/29 09:24:19  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//
//      An overview of the basic philosophy of the ITS code development
// and analysis is show in the figure below.
//Begin_Html
/*
<img src="picts/ITS/ITS_Analysis_schema.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>Roberto Barbera is in charge of the ITS Offline code (1999).
<a href="mailto:roberto.barbera@ct.infn.it">Roberto Barbera</a>.
</font>
<pre>
*/
//End_Html
//
//  AliITS. Inner Traking System base class.
//  This class contains the base procedures for the Inner Tracking System
//
//Begin_Html
/*
<img src="picts/ITS/AliITS_Class_Diagram.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This show the class diagram of the different elements that are part of
the AliITS class.
</font>
<pre>
*/
//End_Html
//
// Version: 0
// Written by Rene Brun, Federico Carminati, and Roberto Barbera
//
// Version: 1
// Modified and documented by Bjorn S. Nilsen
// July 11 1999
//
// Version: 2
// Modified and documented by A. Bologna
// October 18 1999
//
// AliITS is the general base class for the ITS. Also see AliDetector for
// futher information.
//
///////////////////////////////////////////////////////////////////////////////
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <TMath.h>
#include <TRandom.h>
#include <TBranch.h>
#include <TVector.h>
#include <TClonesArray.h>
#include <TROOT.h>
#include <TObjectTable.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TParticle.h>


#include "AliRun.h"
#include "AliITS.h"
#include "AliITSMap.h"
#include "AliITSDetType.h"
#include "AliITSClusterFinder.h"
//#include "AliITSsimulation.h"
#include "AliITSsimulationSPD.h"
#include "AliITSsimulationSDD.h"
#include "AliITSsimulationSSD.h"
#include "AliITSresponse.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSresponseSPD.h"
#include "AliITSresponseSPDbari.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSresponseSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSresponseSSD.h"
#include "AliITShit.h"
#include "AliITSgeom.h"
#include "AliITSdigit.h"
#include "AliITSmodule.h"
#include "AliITSRecPoint.h"
#include "AliITSRawCluster.h"
#include "AliMC.h"
#include "stdlib.h"

#include "AliITStrack.h"
#include "AliITSiotrack.h"
#include "AliITStracking.h"
#include "AliITSRad.h"   
#include "../TPC/AliTPC.h"
#include "../TPC/AliTPCParam.h"


ClassImp(AliITS)
 
//_____________________________________________________________________________
AliITS::AliITS() : AliDetector() {
  //
  // Default initialiser for ITS
  //     The default constructor of the AliITS class. In addition to
  // creating the AliITS class it zeros the variables fIshunt (a member
  // of AliDetector class), fEuclidOut, and fIdN, and zeros the pointers
  // fITSpoints, fIdSens, and fIdName. The AliDetector default constructor
  // is also called.
  //


  fIshunt     = 0;
  fEuclidOut  = 0;

  fNDetTypes = kNTYPES;
  fIdN        = 0;
  fIdName     = 0;
  fIdSens     = 0;
  fITSmodules = 0;
  //
  fDetTypes   = 0;
  //
  fDtype  = 0;
  fNdtype = 0;
  fCtype  = 0;
  fNctype = 0;
  fRecPoints = 0;
  fNRecPoints = 0;
  fTreeC = 0;
  //
  fITSgeom=0;
}

//_____________________________________________________________________________
AliITS::AliITS(const char *name, const char *title):AliDetector(name,title){
  //
  // Default initialiser for ITS
  //     The constructor of the AliITS class. In addition to creating the
  // AliITS class, it allocates memory for the TClonesArrays fHits and
  // fDigits, and for the TObjArray fITSpoints. It also zeros the variables
  // fIshunt (a member of AliDetector class), fEuclidOut, and fIdN, and zeros
  // the pointers fIdSens and fIdName. To help in displaying hits via the ROOT
  // macro display.C AliITS also sets the marker color to red. The variables
  // passes with this constructor, const char *name and *title, are used by
  // the constructor of AliDetector class. See AliDetector class for a
  // description of these parameters and its constructor functions.
  //


  fHits       = new TClonesArray("AliITShit", 1560);
  gAlice->AddHitList(fHits);

  fNDetTypes = kNTYPES;

  fNdtype = new Int_t[kNTYPES];
  fDtype = new TObjArray(kNTYPES);

  fNctype = new Int_t[kNTYPES];
  fCtype = new TObjArray(kNTYPES);


  fRecPoints = 0;
  fNRecPoints = 0;

  fTreeC = 0;

  fITSmodules = 0; 

  fIshunt     = 0;
  fEuclidOut  = 0;
  fIdN        = 0;
  fIdName     = 0;
  fIdSens     = 0;
 
  fDetTypes = new TObjArray(kNTYPES);  

  Int_t i;
  for(i=0;i<kNTYPES;i++) {
    (*fDetTypes)[i]=new AliITSDetType(); 
    fNdtype[i]=0;
    fNctype[i]=0;
   }
  //

  SetMarkerColor(kRed);

  fITSgeom=0;
}
//___________________________________________________________________________
AliITS::AliITS(AliITS &source){
  // copy constructor
  if(this==&source) return;
  Error("AliITS::Copy constructor",
  	"You are not allowed to make a copy of the AliITS");
  exit(1);
}
//____________________________________________________________________________
AliITS& AliITS::operator=(AliITS &source){
  // assignment operator
  if(this==&source) return *this;
  Error("AliITS::operator=",
  	"You are not allowed to make a copy of the AliITS");
  exit(1);
  return *this; //fake return
}
//____________________________________________________________________________
void AliITS::ClearModules(){
  //clear the modules TObjArray

  if(fITSmodules) fITSmodules->Delete();

}
//_____________________________________________________________________________
AliITS::~AliITS(){
  //
  // Default distructor for ITS
  //     The default destructor of the AliITS class. In addition to deleting
  // the AliITS class it deletes the memory pointed to by the fHits, fDigits,
  // fIdSens, fIdName, and fITSpoints.
  //


  delete fHits;
  delete fDigits;
  delete fRecPoints;
//  delete fIdName;        // TObjArray of TObjStrings
  if(fIdName!=0) delete[] fIdName;  // Array of TStrings
  if(fIdSens!=0) delete[] fIdSens;
  if(fITSmodules!=0) {
      this->ClearModules();
      delete fITSmodules;
  }// end if fITSmodules!=0

  //
  if(fDtype) {
    fDtype->Delete();
    delete fDtype;
  }
  delete [] fNdtype;
  if (fCtype) {
    fCtype->Delete();
    delete fCtype;
  }
  delete [] fNctype;
  //

  if (fDetTypes) {
    fDetTypes->Delete();
    delete fDetTypes;
  }

  if (fTreeC) delete fTreeC;

  if (fITSgeom) delete fITSgeom;

}

//___________________________________________
AliITSDetType* AliITS::DetType(Int_t id)
{
  //return pointer to id detector type
    return ((AliITSDetType*) (*fDetTypes)[id]);

}
//___________________________________________
void AliITS::SetClasses(Int_t id, const char *digit, const char *cluster)
{
  //set the digit and cluster classes to be used for the id detector type
    ((AliITSDetType*) (*fDetTypes)[id])->ClassNames(digit,cluster);

}
//___________________________________________
void AliITS::SetResponseModel(Int_t id, AliITSresponse *response)
{
  //set the response model for the id detector type

    ((AliITSDetType*) (*fDetTypes)[id])->ResponseModel(response);

}

//___________________________________________
void AliITS::SetSegmentationModel(Int_t id, AliITSsegmentation *seg)
{
  //set the segmentation model for the id detector type

    ((AliITSDetType*) (*fDetTypes)[id])->SegmentationModel(seg);

}

//___________________________________________
void AliITS::SetSimulationModel(Int_t id, AliITSsimulation *sim)
{
  //set the simulation model for the id detector type

   ((AliITSDetType*) (*fDetTypes)[id])->SimulationModel(sim);

}
//___________________________________________
void AliITS::SetReconstructionModel(Int_t id, AliITSClusterFinder *reconst)
{
  //set the cluster finder model for the id detector type

   ((AliITSDetType*) (*fDetTypes)[id])->ReconstructionModel(reconst);

}

//_____________________________________________________________________________
void AliITS::AddHit(Int_t track, Int_t *vol, Float_t *hits){
  //
  // Add an ITS hit
  //     The function to add information to the AliITShit class. See the
  // AliITShit class for a full description. This function allocates the
  // necessary new space for the hit information and passes the variable
  // track, and the pointers *vol and *hits to the AliITShit constructor
  // function.
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliITShit(fIshunt,track,vol,hits);
}
//_____________________________________________________________________________
void AliITS::AddRealDigit(Int_t id, Int_t *digits) 
{
  // add a real digit - as coming from data

  TClonesArray &ldigits = *((TClonesArray*)(*fDtype)[id]);
  new(ldigits[fNdtype[id]++]) AliITSdigit(digits);

}
//_____________________________________________________________________________
void AliITS::AddSimDigit(Int_t id, AliITSdigit *d) 
{

  // add a simulated digit

  TClonesArray &ldigits = *((TClonesArray*)(*fDtype)[id]);

  switch(id)
  {
  case 0:
     new(ldigits[fNdtype[id]++]) AliITSdigitSPD(*((AliITSdigitSPD*)d));
     break;
  case 1:
     new(ldigits[fNdtype[id]++]) AliITSdigitSDD(*((AliITSdigitSDD*)d));
     break;
  case 2:
     new(ldigits[fNdtype[id]++]) AliITSdigitSSD(*((AliITSdigitSSD*)d));
     break;
  }

}

//_____________________________________________________________________________
void AliITS::AddSimDigit(Int_t id,Float_t phys,Int_t *digits,Int_t *tracks,Int_t *hits,Float_t *charges){

  // add a simulated digit to the list

  TClonesArray &ldigits = *((TClonesArray*)(*fDtype)[id]);
  switch(id)
  {
  case 0:
     new(ldigits[fNdtype[id]++]) AliITSdigitSPD(digits,tracks,hits);
     break;
  case 1:
     new(ldigits[fNdtype[id]++]) AliITSdigitSDD(phys,digits,tracks,hits,charges);
     break;
  case 2:
     new(ldigits[fNdtype[id]++]) AliITSdigitSSD(digits,tracks,hits);
     break;
  }
 
}

//_____________________________________________________________________________
void AliITS::AddCluster(Int_t id, AliITSRawCluster *c) 
{

  // add a cluster to the list

  TClonesArray &lcl = *((TClonesArray*)(*fCtype)[id]);

  switch(id)
  {
  case 0:
     new(lcl[fNctype[id]++]) AliITSRawClusterSPD(*((AliITSRawClusterSPD*)c));
     break;
  case 1:
     new(lcl[fNctype[id]++]) AliITSRawClusterSDD(*((AliITSRawClusterSDD*)c));
     break;
  case 2:
     new(lcl[fNctype[id]++]) AliITSRawClusterSSD(*((AliITSRawClusterSSD*)c));
     break;
  }

}


//_____________________________________________________________________________
void AliITS::AddRecPoint(const AliITSRecPoint &r)
{
  //
  // Add a reconstructed space point to the list
  //
  TClonesArray &lrecp = *fRecPoints;
  new(lrecp[fNRecPoints++]) AliITSRecPoint(r);
}
 

//____________________________________________
void AliITS::ResetDigits()
{
    //
    // Reset number of digits and the digits array for the ITS detector
    //

    if (!fDtype) return;

    Int_t i;
    for (i=0;i<kNTYPES;i++ ) {
	if ((*fDtype)[i])    ((TClonesArray*)(*fDtype)[i])->Clear();
	if (fNdtype)  fNdtype[i]=0;
    }
}

//____________________________________________
void AliITS::ResetDigits(Int_t i)
{
    //
    // Reset number of digits and the digits array for this branch
    //
  if ((*fDtype)[i])    ((TClonesArray*)(*fDtype)[i])->Clear();
  if (fNdtype)  fNdtype[i]=0;
}


//____________________________________________
void AliITS::ResetClusters()
{
    //
    // Reset number of clusters and the clusters array for ITS
    //

    Int_t i;
    for (i=0;i<kNTYPES;i++ ) {
	if ((*fCtype)[i])    ((TClonesArray*)(*fCtype)[i])->Clear();
	if (fNctype)  fNctype[i]=0;
    }

}

//____________________________________________
void AliITS::ResetClusters(Int_t i)
{
    //
    // Reset number of clusters and the clusters array for this branch
    //
	if ((*fCtype)[i])    ((TClonesArray*)(*fCtype)[i])->Clear();
	if (fNctype)  fNctype[i]=0;

}


//____________________________________________
void AliITS::ResetRecPoints()
{
    //
    // Reset number of rec points and the rec points array 
    //
    if (fRecPoints) fRecPoints->Clear();
    fNRecPoints = 0;

}

//_____________________________________________________________________________
Int_t AliITS::DistancetoPrimitive(Int_t , Int_t ){
  //
  // Distance from mouse to ITS on the screen. Dummy routine
  //     A dummy routine used by the ROOT macro display.C to allow for the
  // use of the mouse (pointing device) in the macro. In general this should
  // never be called. If it is it returns the number 9999 for any value of
  // x and y.
  //
  return 9999;
}

//_____________________________________________________________________________
void AliITS::Init(){
  //
  // Initialise ITS after it has been built
  //     This routine initializes the AliITS class. It is intended to be called
  // from the Init function in AliITSv?. Besides displaying a banner
  // indicating that it has been called it initializes the array fIdSens
  // and sets the default segmentation, response, digit and raw cluster classes
  // Therefore it should be called after a call to CreateGeometry.
  //
  Int_t i;

//
  SetDefaults();
// Array of TStrings
  for(i=0;i<fIdN;i++) fIdSens[i] = gMC->VolId(fIdName[i]);
//
}

//_____________________________________________________________________________
void AliITS::SetDefaults()
{
  // sets the default segmentation, response, digit and raw cluster classes

  printf("SetDefaults\n");

  AliITSDetType *iDetType;


  //SPD 

  iDetType=DetType(0); 
  if (!iDetType->GetSegmentationModel()) {
    AliITSsegmentationSPD *seg0=new AliITSsegmentationSPD(fITSgeom);
    SetSegmentationModel(0,seg0); 
  }
  if (!iDetType->GetResponseModel()) {
     SetResponseModel(0,new AliITSresponseSPD()); 
  }
  // set digit and raw cluster classes to be used
  
  const char *kData0=(iDetType->GetResponseModel())->DataType();
  if (strstr(kData0,"real")) {
      iDetType->ClassNames("AliITSdigit","AliITSRawClusterSPD");
  } else iDetType->ClassNames("AliITSdigitSPD","AliITSRawClusterSPD");

  // SDD					  //
  iDetType=DetType(1); 
  if (!iDetType->GetResponseModel()) {
    SetResponseModel(1,new AliITSresponseSDD()); 
  }
  AliITSresponse *resp1=iDetType->GetResponseModel();
  if (!iDetType->GetSegmentationModel()) {
    AliITSsegmentationSDD *seg1=new AliITSsegmentationSDD(fITSgeom,resp1);
    SetSegmentationModel(1,seg1); 
  }
  const char *kData1=(iDetType->GetResponseModel())->DataType();
  const char *kopt=iDetType->GetResponseModel()->ZeroSuppOption();
  if ((!strstr(kopt,"2D")) && (!strstr(kopt,"1D")) || strstr(kData1,"real") ) {
      iDetType->ClassNames("AliITSdigit","AliITSRawClusterSDD");
  } else iDetType->ClassNames("AliITSdigitSDD","AliITSRawClusterSDD");

  // SSD
  iDetType=DetType(2); 
  if (!iDetType->GetSegmentationModel()) {
    AliITSsegmentationSSD *seg2=new AliITSsegmentationSSD(fITSgeom);
    SetSegmentationModel(2,seg2); 
  }
  if (!iDetType->GetResponseModel()) {
    SetResponseModel(2,new AliITSresponseSSD()); 
  }
  const char *kData2=(iDetType->GetResponseModel())->DataType();
  if (strstr(kData2,"real")) {
      iDetType->ClassNames("AliITSdigit","AliITSRawClusterSSD");
  } else iDetType->ClassNames("AliITSdigitSSD","AliITSRawClusterSSD");

  if (kNTYPES>3) {
    Warning("SetDefaults","Only the three basic detector types are initialised!");
  } 

}



//_____________________________________________________________________________
void AliITS::SetDefaultSimulation()
{
  // to be written

}
//_____________________________________________________________________________
void AliITS::SetDefaultClusterFinders()
{
  // to be written

}
//_____________________________________________________________________________

void AliITS::MakeTreeC(Option_t *option)
{
  // create a separate tree to store the clusters

  cout << "AliITS::MakeTreeC" << endl;

     const char *optC = strstr(option,"C");
     if (optC && !fTreeC) fTreeC = new TTree("TC","Clusters in ITS");
     else return;

     Int_t buffersize = 4000;
     char branchname[30];

     const char *det[3] = {"SPD","SDD","SSD"};

     char digclass[40];
     char clclass[40];

     // one branch for Clusters per type of detector
     Int_t i;
     for (i=0; i<kNTYPES ;i++) {
       AliITSDetType *iDetType=DetType(i); 
       iDetType->GetClassNames(digclass,clclass);
       // clusters
       (*fCtype)[i] = new TClonesArray(clclass,10000); 
        if (kNTYPES==3) sprintf(branchname,"%sClusters%s",GetName(),det[i]);
	else  sprintf(branchname,"%sClusters%d",GetName(),i+1);
	if (fCtype   && fTreeC) {
	   TreeC()->Branch(branchname,&((*fCtype)[i]), buffersize);
	   cout << "Making Branch " << branchname;
	   cout << " for Clusters of detector type " << i+1 << endl;
	}	
     }

}

//_____________________________________________________________________________
void AliITS::GetTreeC(Int_t event)
{

  cout << "AliITS::GetTreeC" << endl;

  // get the clusters tree for this event and set the branch address
    char treeName[20];
    char branchname[30];

    const char *det[3] = {"SPD","SDD","SSD"};

    ResetClusters();
    if (fTreeC) {
	  delete fTreeC;
    }

    sprintf(treeName,"TreeC%d",event);
    fTreeC = (TTree*)gDirectory->Get(treeName);

    TBranch *branch;
    if (fTreeC) {
        Int_t i;
	for (i=0; i<kNTYPES; i++) {
	   if (kNTYPES==3) sprintf(branchname,"%sClusters%s",GetName(),det[i]);
	   else  sprintf(branchname,"%sClusters%d",GetName(),i+1);
	   if (fCtype) {
		branch = fTreeC->GetBranch(branchname);
                if (branch) branch->SetAddress(&((*fCtype)[i]));
	    }
	}
    } else {
	Error("AliITS::GetTreeC",
		"cannot find Clusters Tree for event:%d\n",event);
    }

}
//_____________________________________________________________________________
void AliITS::MakeBranch(Option_t* option, char *file)
{
  //
  // Creates Tree branches for the ITS.
  //
  //
  Int_t buffersize = 4000;
  char branchname[30];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option,file);

  const char *cD = strstr(option,"D");
  const char *cR = strstr(option,"R");

  if (cD) {
  //
  // one branch for digits per type of detector
  //
   const char *det[3] = {"SPD","SDD","SSD"};

   char digclass[40];
   char clclass[40];

   Int_t i;
   for (i=0; i<kNTYPES ;i++) {
       AliITSDetType *iDetType=DetType(i); 
       iDetType->GetClassNames(digclass,clclass);
       // digits
       if(!((*fDtype)[i])) (*fDtype)[i] = new TClonesArray(digclass,10000);
       else ResetDigits(i);
   }

   for (i=0; i<kNTYPES ;i++) {
      if (kNTYPES==3) sprintf(branchname,"%sDigits%s",GetName(),det[i]);
      else  sprintf(branchname,"%sDigits%d",GetName(),i+1);      
      if (fDtype && gAlice->TreeD()) {
        gAlice->MakeBranchInTree(gAlice->TreeD(), 
                         branchname, &((*fDtype)[i]), buffersize, file);
	    cout << "Making Branch " << branchname;
	    cout << " for digits of type "<< i+1 << endl;
      }	
    }
  }

  if (cR) {
  //
  // only one branch for rec points for all detector types
  //
    sprintf(branchname,"%sRecPoints",GetName());
 
    if(!fRecPoints) fRecPoints=new TClonesArray("AliITSRecPoint",10000);

    if (fRecPoints && gAlice->TreeR()) {
      gAlice->MakeBranchInTree(gAlice->TreeR(), 
                               branchname, &fRecPoints, buffersize, file) ;
      cout << "Making Branch " << branchname;
      cout << " for reconstructed space points" << endl;
    }
  }	
}

//___________________________________________
void AliITS::SetTreeAddress()
{

  // Set branch address for the Trees.

  char branchname[30];
  AliDetector::SetTreeAddress();

  const char *det[3] = {"SPD","SDD","SSD"};

  TBranch *branch;
  TTree *treeD = gAlice->TreeD();
  TTree *treeR = gAlice->TreeR();

  Int_t i;
  if (treeD) {
      for (i=0; i<kNTYPES; i++) {
	  if (kNTYPES==3) sprintf(branchname,"%sDigits%s",GetName(),det[i]);
	  else  sprintf(branchname,"%sDigits%d",GetName(),i+1);
	  if (fDtype) {
	      branch = treeD->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fDtype)[i]));
	  }
      }
  }

 
  if (treeR) {
    sprintf(branchname,"%sRecPoints",GetName());
      branch = treeR->GetBranch(branchname);
      if (branch) branch->SetAddress(&fRecPoints);
  }
  

}

//____________________________________________________________________________
void AliITS::InitModules(Int_t size,Int_t &nmodules){

  //initialize the modules array

  if(fITSmodules){ 
      fITSmodules->Delete();
      delete fITSmodules;
  }

    Int_t nl,indexMAX,index;

    if(size<=0){ // default to using data stored in AliITSgeom
	if(fITSgeom==0) {
	    Error("AliITS::InitModules",
	    	"in AliITS::InitModule fITSgeom not defined\n");
	    return;
	} // end if fITSgeom==0
	nl = fITSgeom->GetNlayers();
	indexMAX = fITSgeom->GetModuleIndex(nl,fITSgeom->GetNladders(nl),
					    fITSgeom->GetNdetectors(nl))+1;
	nmodules = indexMAX;
	fITSmodules = new TObjArray(indexMAX);
	for(index=0;index<indexMAX;index++){
		fITSmodules->AddAt( new AliITSmodule(index),index);
	} // end for index
    }else{
	fITSmodules = new TObjArray(size);
	for(index=0;index<size;index++) {
	    fITSmodules->AddAt( new AliITSmodule(index),index);
	}

        nmodules = size;
    } // end i size<=0
}

//____________________________________________________________________________
void AliITS::FillModules(Int_t evnt,Int_t bgrev,Int_t nmodules,Option_t *option,Text_t *filename){

  // fill the modules with the sorted by module hits; add hits from background
  // if option=Add


    static TTree *trH1;                 //Tree with background hits
    static TClonesArray *fHits2;        //List of hits for one track only

    static Bool_t first=kTRUE;
    static TFile *file;
    const char *addBgr = strstr(option,"Add");


    if (addBgr ) {
	if(first) {
	    cout<<"filename "<<filename<<endl;
	    file=new TFile(filename);
	    cout<<"I have opened "<<filename<<" file "<<endl;
	    fHits2     = new TClonesArray("AliITShit",1000  );
	}	    
	first=kFALSE;
	file->cd();
	file->ls();
	// Get Hits Tree header from file
	if(fHits2) fHits2->Clear();
	if(trH1) delete trH1;
	trH1=0;
	
	char treeName[20];
	sprintf(treeName,"TreeH%d",bgrev);
	trH1 = (TTree*)gDirectory->Get(treeName);
        //printf("TrH1 %p of treename %s for event %d \n",trH1,treeName,bgrev);
	
	if (!trH1) {
	    Error("AliITS::FillModules",
	    "cannot find Hits Tree for event:%d\n",bgrev);
	}
	// Set branch addresses
	TBranch *branch;
	char branchname[20];
	sprintf(branchname,"%s",GetName());
	if (trH1 && fHits2) {
	    branch = trH1->GetBranch(branchname);
	    if (branch) branch->SetAddress(&fHits2);
	}

        // test
	//Int_t ntracks1 =(Int_t)TrH1->GetEntries();
	//printf("background - ntracks1 - %d\n",ntracks1);
   }

    //Int_t npart = gAlice->GetEvent(evnt);
    //if(npart<=0) return;
    TClonesArray *itsHits = this->Hits();
    Int_t lay,lad,det,index;
    AliITShit *itsHit=0;
    AliITSmodule *mod=0;

    TTree *iTH = gAlice->TreeH();
    Int_t ntracks =(Int_t) iTH->GetEntries();

    Int_t t,h;
    for(t=0; t<ntracks; t++){
	gAlice->ResetHits();
	iTH->GetEvent(t);
	Int_t nhits = itsHits->GetEntriesFast();
	//printf("nhits %d\n",nhits);
        if (!nhits) continue;
	for(h=0; h<nhits; h++){
	    itsHit = (AliITShit *)itsHits->UncheckedAt(h);
	    itsHit->GetDetectorID(lay,lad,det);
	    // temporarily index=det-1 !!!
	    if(fITSgeom) index = fITSgeom->GetModuleIndex(lay,lad,det);
	    else index=det-1;
	    //
	    mod = this->GetModule(index);
	    mod->AddHit(itsHit,t,h);
	} // end loop over hits 
    } // end loop over tracks

    // open the file with background
    
    if (addBgr ) {
          Int_t track,i;
          ntracks =(Int_t)trH1->GetEntries();
	    //printf("background - ntracks1 %d\n",ntracks);
	    //printf("background - Start loop over tracks \n");     
            //   Loop over tracks

	    for (track=0; track<ntracks; track++) {

		if (fHits2)       fHits2->Clear();
		trH1->GetEvent(track);
                //   Loop over hits
		for(i=0;i<fHits2->GetEntriesFast();++i) {

		    itsHit=(AliITShit*) (*fHits2)[i];
		    itsHit->GetDetectorID(lay,lad,det);
		    // temporarily index=det-1 !!!
		    if(fITSgeom) index = fITSgeom->GetModuleIndex(lay,lad,det);
		    else index=det-1;
		    //
		    mod = this->GetModule(index);
		    mod->AddHit(itsHit,track,i);
	       }  // end loop over hits
	    } // end loop over tracks

	    TTree *fAli=gAlice->TreeK();
            TFile *fileAli=0;
	    
	    if (fAli) fileAli =fAli->GetCurrentFile();
	    fileAli->cd();

    } // end if add

    //gObjectTable->Print();

}

//____________________________________________________________________________

void AliITS::SDigits2Digits()
{
  
  AliITSgeom *geom = GetITSgeom();
  
  // SPD
  AliITSDetType *iDetType;
  iDetType=DetType(0);
  AliITSsegmentationSPD *seg0=(AliITSsegmentationSPD*)iDetType->GetSegmentationModel();
  AliITSresponseSPD *res0 = (AliITSresponseSPD*)iDetType->GetResponseModel();
  AliITSsimulationSPD *sim0=new AliITSsimulationSPD(seg0,res0);
  SetSimulationModel(0,sim0);
  // test
  // printf("SPD dimensions %f %f \n",seg0->Dx(),seg0->Dz());
  // printf("SPD npixels %d %d \n",seg0->Npz(),seg0->Npx());
  // printf("SPD pitches %d %d \n",seg0->Dpz(0),seg0->Dpx(0));
  // end test
  // 
  // SDD
  //Set response functions
  // SDD compression param: 2 fDecrease, 2fTmin, 2fTmax or disable, 2 fTolerance

  iDetType=DetType(1);
  AliITSresponseSDD *res1 = (AliITSresponseSDD*)iDetType->GetResponseModel();
  if (!res1) {
    res1=new AliITSresponseSDD();
    SetResponseModel(1,res1);
  }
  Float_t noise, baseline;
  res1->GetNoiseParam(noise,baseline);
  Float_t noise_after_el =  res1->GetNoiseAfterElectronics();
  Float_t fCutAmp = baseline + 2.*noise_after_el;
  Int_t cp[8]={0,0,(int)fCutAmp,(int)fCutAmp,0,0,0,0}; //1D
  res1->SetCompressParam(cp);
  AliITSsegmentationSDD *seg1=(AliITSsegmentationSDD*)iDetType->GetSegmentationModel();
  if (!seg1) {
    seg1 = new AliITSsegmentationSDD(geom,res1);
    SetSegmentationModel(1,seg1);
  }
  AliITSsimulationSDD *sim1=new AliITSsimulationSDD(seg1,res1);
  SetSimulationModel(1,sim1);

  // SSD
  iDetType=DetType(2);
  AliITSsegmentationSSD *seg2=(AliITSsegmentationSSD*)iDetType->GetSegmentationModel();
  AliITSresponseSSD *res2 = (AliITSresponseSSD*)iDetType->GetResponseModel();
  res2->SetSigmaSpread(3.,2.);
  AliITSsimulationSSD *sim2=new AliITSsimulationSSD(seg2,res2);
  SetSimulationModel(2,sim2);
 
  cerr<<"Digitizing ITS...\n";
 
  TStopwatch timer;
  timer.Start();
  HitsToDigits(0,0,-1," ","All"," ");
  timer.Stop(); timer.Print();

  delete sim0;
  delete sim1;
  delete sim2;
}


//____________________________________________________________________________
void AliITS::HitsToDigits(Int_t evNumber,Int_t bgrev,Int_t size, Option_t *option, Option_t *opt,Text_t *filename)
{
    // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !
  
   // the condition below will disappear when the geom class will be
   // initialised for all versions - for the moment it is only for v5 !
   // 7 is the SDD beam test version  
   Int_t ver = this->IsVersion(); 
   if(ver!=5 && ver!=7) return; 

   const char *all = strstr(opt,"All");
   const char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),strstr(opt,"SSD")};

   Int_t nmodules;
   InitModules(size,nmodules); 
   FillModules(evNumber,bgrev,nmodules,option,filename);

   //TBranch *branch;
   AliITSsimulation* sim;
   //TObjArray *branches=gAlice->TreeD()->GetListOfBranches();
   AliITSgeom *geom = GetITSgeom();

   Int_t id,module;
   Int_t first,last;
   for (id=0;id<kNTYPES;id++) {
        if (!all && !det[id]) continue;
	//branch = (TBranch*)branches->UncheckedAt(id);
	AliITSDetType *iDetType=DetType(id); 
	sim = (AliITSsimulation*)iDetType->GetSimulationModel();
	if (!sim) {
           Error("HitsToDigits","The simulation class was not instantiated!");
           exit(1);
	   // or SetDefaultSimulation();
	}
	if(geom) {
	  first = geom->GetStartDet(id);
	  last = geom->GetLastDet(id);
	} else first=last=0;
	cout << "det type " << id << " first, last "<< first << last << endl;
	for(module=first;module<=last;module++) {
	    AliITSmodule *mod = (AliITSmodule *)fITSmodules->At(module);
	    sim->DigitiseModule(mod,module,evNumber);
	    // fills all branches - wasted disk space
	    gAlice->TreeD()->Fill(); 
	    ResetDigits();
	    // try and fill only the branch 
	    //branch->Fill();
	    //ResetDigits(id);
	} // loop over modules
   } // loop over detector types

   ClearModules();

   Int_t nentries=(Int_t)gAlice->TreeD()->GetEntries();
   cout << "nentries in TreeD" << nentries << endl;

   char hname[30];
   sprintf(hname,"TreeD%d",evNumber);
   gAlice->TreeD()->Write(hname,TObject::kOverwrite);
   // reset tree
   gAlice->TreeD()->Reset();

}


//____________________________________________________________________________
void AliITS::DigitsToRecPoints(Int_t evNumber,Int_t lastentry,Option_t *opt)
{
  // cluster finding and reconstruction of space points
  
   // the condition below will disappear when the geom class will be
   // initialised for all versions - for the moment it is only for v5 !
   // 7 is the SDD beam test version  
   Int_t ver = this->IsVersion(); 
   if(ver!=5) return; 

   const char *all = strstr(opt,"All");
   const char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),strstr(opt,"SSD")};

   static Bool_t first=kTRUE;
   if (!TreeC() && first) {
       MakeTreeC("C");
       first=kFALSE;
   }

   TTree *treeC=TreeC();
 

   //TBranch *branch;
   AliITSClusterFinder* rec;

   //TObjArray *branches=gAlice->TreeR()->GetListOfBranches();
   AliITSgeom *geom = GetITSgeom();

   Int_t id,module;
   for (id=0;id<kNTYPES;id++) {
        if (!all && !det[id]) continue;
	//branch = (TBranch*)branches->UncheckedAt(id);
	AliITSDetType *iDetType=DetType(id); 
	rec = (AliITSClusterFinder*)iDetType->GetReconstructionModel();
	if (!rec) {
           Error("DigitsToRecPoints","The cluster finder class was not instantiated!");
           exit(1);
	   // or SetDefaultClusterFinders();
	}
        TClonesArray *itsDigits  = this->DigitsAddress(id);

        Int_t first,last;
	if(geom) {
	  first = geom->GetStartDet(id);
	  last = geom->GetLastDet(id);
	} else first=last=0;
	//printf("first last %d %d\n",first,last);
	for(module=first;module<=last;module++) {
              this->ResetDigits();
              if (all) gAlice->TreeD()->GetEvent(lastentry+module);
	      else gAlice->TreeD()->GetEvent(lastentry+(module-first));
	      Int_t ndigits = itsDigits->GetEntriesFast();
	      if (ndigits) rec->FindRawClusters();
	      gAlice->TreeR()->Fill(); 
	      ResetRecPoints();
	      treeC->Fill();
              ResetClusters();
	      // try and fill only the branch 
	      //branch->Fill();
	      //ResetRecPoints(id);
	} // loop over modules
   } // loop over detector types


   Int_t nentries=(Int_t)gAlice->TreeR()->GetEntries();
   Int_t ncentries=(Int_t)treeC->GetEntries();
   cout << " nentries ncentries " << nentries << ncentries <<  endl;

   char hname[30];
   sprintf(hname,"TreeR%d",evNumber);
   gAlice->TreeR()->Write(hname,TObject::kOverwrite);
   // reset tree
   gAlice->TreeR()->Reset();

   sprintf(hname,"TreeC%d",evNumber);
   treeC->Write(hname,TObject::kOverwrite);
   treeC->Reset();
}


//____________________________________________________________________________
void AliITS::HitsToFastRecPoints(Int_t evNumber,Int_t bgrev,Int_t size,
Option_t *option,Option_t *opt,Text_t *filename)
{
    // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !
  

   // the condition below will disappear when the geom class will be
   // initialised for all versions - for the moment it is only for v5 !  
   Int_t ver = this->IsVersion(); 
   if(ver!=5) return; 

   const char *all = strstr(opt,"All");
   const char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),strstr(opt,"SSD")};

   Int_t nmodules;
   InitModules(size,nmodules);
   FillModules(evNumber,bgrev,nmodules,option,filename);


   AliITSsimulation* sim;
   AliITSgeom *geom = GetITSgeom();

   TRandom *random=new TRandom[9];
   random[0].SetSeed(111);
   random[1].SetSeed(222);
   random[2].SetSeed(333);		
   random[3].SetSeed(444);
   random[4].SetSeed(555);
   random[5].SetSeed(666);		
   random[6].SetSeed(777);
   random[7].SetSeed(888);
   random[8].SetSeed(999);		


   Int_t id,module;
   for (id=0;id<kNTYPES;id++) {
        if (!all && !det[id]) continue;
	AliITSDetType *iDetType=DetType(id); 
	sim = (AliITSsimulation*)iDetType->GetSimulationModel();
	if (!sim) {
           Error("HitsToFastPoints","The simulation class was not instantiated!");
           exit(1);
	   // or SetDefaultSimulation();
	}
	Int_t first = geom->GetStartDet(id);
	Int_t last = geom->GetLastDet(id);
	for(module=first;module<=last;module++) {
	    AliITSmodule *mod = (AliITSmodule *)fITSmodules->At(module);
	    sim->CreateFastRecPoints(mod,module,random);
	    gAlice->TreeR()->Fill(); 
	    ResetRecPoints();
	} // loop over modules
   } // loop over detector types


   ClearModules();

   //Int_t nentries=(Int_t)gAlice->TreeR()->GetEntries();

   char hname[30];
   sprintf(hname,"TreeR%d",evNumber);
   gAlice->TreeR()->Write(hname,TObject::kOverwrite);
   // reset tree
   gAlice->TreeR()->Reset();

   delete [] random;

}

//________________________________________________________________

AliITStrack  AliITS::Tracking(AliITStrack &track, AliITStrack *reference,TObjArray *fastpoints, Int_t
**vettid, Bool_t flagvert,  AliITSRad *rl  ) { 
										
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it 

  										    
  TList *list= new TList();   

  AliITStrack tr(track);
  
  list->AddLast(&tr);
  
  Double_t Pt=(tr).GetPt();
  //cout << "\n Pt = " << Pt <<"\n";

  AliITStracking obj(list, reference, this, fastpoints,TMath::Abs(Pt),vettid, flagvert, rl);
  list->Delete();
  delete list;

  Int_t itot=-1;
  TVector VecTotLabref(18);
  Int_t lay, k;
  for(lay=5; lay>=0; lay--) {
    TVector VecLabref(3); 
    VecLabref=(*reference).GetLabTrack(lay);
    Float_t ClustZ=(*reference).GetZclusterTrack( lay);   //modified il 5-3-2001	 
    for(k=0; k<3; k++) { // {itot++; VecTotLabref(itot)=VecLabref(k);}  // modified 5-3-2002
 
	Int_t lpp=(Int_t)VecLabref(k);
	if(lpp>=0) {
	  TParticle *p=(TParticle*) gAlice->Particle(lpp);
	  Int_t pcode=p->GetPdgCode();
	  if(pcode==11) VecLabref(k)=p->GetFirstMother();
	}    
    
        itot++; VecTotLabref(itot)=VecLabref(k);
        if(VecLabref(k)==0. && ClustZ == 0.) VecTotLabref(itot) =-3.;    
    }    
  }
  Long_t labref;
  Int_t freq;  
  (*reference).Search(VecTotLabref, labref, freq);
    
  if(freq < 5) labref=-labref;   	
  (*reference).SetLabel(labref);

  return *reference; 

}



//________________________________________________________________



void AliITS::DoTracking(Int_t evNumber, Int_t min_t, Int_t max_t, TFile *file, Bool_t flagvert) {

//   ex macro for tracking ITS

  printf("begin DoTracking - file %p\n",file);

  const char *pname="75x40_100x60";
  
 Int_t imax=200,jmax=450;
  AliITSRad *rl = new AliITSRad(imax,jmax);
  //cout<<" dopo costruttore AliITSRad\n"; getchar();
    
  struct GoodTrack {
    Int_t lab,code;
    Float_t px,py,pz,x,y,z,pxg,pyg,pzg,ptg;
    Bool_t flag;
  };
  

  gAlice->GetEvent(0);

  AliTPC *TPC=(AliTPC*)gAlice->GetDetector("TPC");
  AliTPCParam *digp = (AliTPCParam*)file->Get(pname);
  if (digp!=0) TPC->SetParam(digp);
  
  GoodTrack gt[15000];
  Int_t ngood=0;
  ifstream in("itsgood_tracks");

  cerr<<"Reading itsgood tracks...\n";
  while (in>>gt[ngood].lab>>gt[ngood].code
	  >>gt[ngood].px >>gt[ngood].py>>gt[ngood].pz
	  >>gt[ngood].x  >>gt[ngood].y >>gt[ngood].z
	  >>gt[ngood].pxg  >>gt[ngood].pyg >>gt[ngood].pzg
	  >>gt[ngood].ptg >>gt[ngood].flag) {
    ngood++;
    cerr<<ngood<<'\r';
    if (ngood==15000) {
      cerr<<"Too many good tracks !\n";
      break;
    }
  }
  if (!in.eof()) cerr<<"Read error (itsgood_tracks) !\n";
  
  
// Load tracks
  TFile *tf=TFile::Open("tpctracks.root");
  if (!tf->IsOpen()) {cerr<<"Can't open tpctracks.root !\n"; return ;}
  TObjArray tracks(200000);
  TTree *tracktree=(TTree*)tf->Get("TreeT");
  TBranch *tbranch=tracktree->GetBranch("tracks");
  Int_t nentr=(Int_t)tracktree->GetEntries();
  Int_t kk;
  for (kk=0; kk<nentr; kk++) {
    AliTPCtrack *iotrack=new AliTPCtrack;
    tbranch->SetAddress(&iotrack);
    tracktree->GetEvent(kk);
    tracks.AddLast(iotrack);
  }   
  tf->Close();


  Int_t nt = tracks.GetEntriesFast();
  cerr<<"Number of found tracks "<<nt<<endl;
  
  TVector DataOut(9);
  Int_t kkk=0;
  
  Double_t ptg=0.,pxg=0.,pyg=0.,pzg=0.;

  //////////////////////////////  good tracks definition in TPC  ////////////////////////////////
      
  ofstream out1 ("AliITSTrag.out");
  Int_t i;
  for (i=0; i<ngood; i++) out1 << gt[i].ptg << "\n";
  out1.close();


  TVector vec(5);
  TTree *TR=gAlice->TreeR();
  Int_t nent=(Int_t)TR->GetEntries();
  TClonesArray  *recPoints = RecPoints();
  Int_t numbpoints;
  Int_t totalpoints=0;
  Int_t *np = new Int_t[nent];
  Int_t **vettid = new Int_t* [nent];
  Int_t mod;
  for (mod=0; mod<nent; mod++) {
    vettid[mod]=0;
    this->ResetRecPoints();
    //gAlice->TreeR()->GetEvent(mod+1); //first entry in TreeR is empty
    gAlice->TreeR()->GetEvent(mod); //first entry in TreeR is empty
    numbpoints = recPoints->GetEntries();
    totalpoints+=numbpoints;
    np[mod] = numbpoints;
  //cout<<" mod = "<<mod<<"   numbpoints = "<<numbpoints<<"\n"; getchar();
    vettid[mod] = new Int_t[numbpoints];
    Int_t ii;
    for (ii=0;ii<numbpoints; ii++) *(vettid[mod]+ii)=0;
  }

  AliTPCtrack *track;

     
  if(min_t < 0) {min_t = 0; max_t = nt-1;}   

/*
  ///////////////////////////////// Definition of vertex end its error ////////////////////////////
  ////////////////////////// In the future it will be given by a method ///////////////////////////
  Double_t Vx=0.;
  Double_t Vy=0.;
  Double_t Vz=0.;
  
  Float_t sigmavx=0.0050;      // 50  microns
  Float_t sigmavy=0.0050;      // 50  microns
  Float_t sigmavz=0.010;       // 100 microns

  //Vx+=gRandom->Gaus(0,sigmavx);  Vy+=gRandom->Gaus(0,sigmavy);  Vz+=gRandom->Gaus(0,sigmavz);
  TVector vertex(3), ervertex(3)
  vertex(0)=Vx; vertex(1)=Vy; vertex(2)=Vz;
  ervertex(0)=sigmavx;  ervertex(1)=sigmavy;  ervertex(2)=sigmavz;
  /////////////////////////////////////////////////////////////////////////////////////////////////
*/      

  //TDirectory *savedir=gDirectory; 

  TTree tracktree1("TreeT","Tree with ITS tracks");
  AliITSiotrack *iotrack=0;
  tracktree1.Branch("ITStracks","AliITSiotrack",&iotrack,32000,0);

  ofstream out ("AliITSTra.out");
  
  Int_t j;       
  for (j=min_t; j<=max_t; j++) {     
    track=(AliTPCtrack*)tracks.UncheckedAt(j);
    Int_t flaglab=0;
    if (!track) continue;
    ////// elimination of not good tracks ////////////	 
    Int_t ilab=TMath::Abs(track->GetLabel());
    Int_t iii;
    for (iii=0;iii<ngood;iii++) {
	 //cout<<" ilab, gt[iii].lab = "<<ilab<<" "<<gt[iii].lab<<"\n"; getchar();
      if (ilab==gt[iii].lab) { 
	flaglab=1;
	ptg=gt[iii].ptg; 
	pxg=gt[iii].pxg;
	pyg=gt[iii].pyg;
	pzg=gt[iii].pzg;	
	break;
      }
    }
	 //cout<<" j flaglab =  " <<j<<" "<<flaglab<<"\n";  getchar();
    if (!flaglab) continue;
	 //cout<<" j =  " <<j<<"\n";  getchar();
  /*		
    ////// old propagation to the end of TPC //////////////       
    Double_t xk=76.;
    track->PropagateTo(xk);
    xk-=0.11;
    track->PropagateTo(xk,42.7,2.27); //C
    xk-=2.6;
    track->PropagateTo(xk,36.2,1.98e-3); //C02
    xk-=0.051;
    track->PropagateTo(xk,42.7,2.27); //C 
    /////////////////////////////////////////////////// 
	 */
	 	 
	 ////// new propagation to the end of TPC //////////////
    Double_t xk=77.415;
    track->PropagateTo(xk, 28.94, 1.204e-3);	 //Ne
	 xk -=0.01;
    track->PropagateTo(xk, 44.77, 1.71);	 //Tedlar
	 xk -=0.04;
    track->PropagateTo(xk, 44.86, 1.45);	 //Kevlar
	 xk -=2.0;
    track->PropagateTo(xk, 41.28, 0.029);	 //Nomex	 
    xk-=16;
    track->PropagateTo(xk,36.2,1.98e-3); //C02
	 xk -=0.01;
    track->PropagateTo(xk, 24.01, 2.7);	 //Al	 
	 xk -=0.01;
    track->PropagateTo(xk, 44.77, 1.71);	 //Tedlar
	 xk -=0.04;
    track->PropagateTo(xk, 44.86, 1.45);	 //Kevlar
	 xk -=0.5;
    track->PropagateTo(xk, 41.28, 0.029);	 //Nomex	 	 	 	 	 	 
	    	     
       /////////////////////////////////////////////////////////////// 	 		 
  
   ///////////////////////////////////////////////////////////////
    AliITStrack trackITS(*track);
    AliITStrack result(*track);
    AliITStrack primarytrack(*track); 
    
///////////////////////////////////////////////////////////////////////////////////////////////
	 TVector Vgeant(3);
	 Vgeant=result.GetVertex(); 
			  
  // Definition of Dv and Zv for vertex constraint	
     Double_t sigmaDv=0.0050;  Double_t sigmaZv=0.010;	
    //Double_t sigmaDv=0.0015;  Double_t sigmaZv=0.0015;				  
	Double_t uniform= gRandom->Uniform();
	Double_t signdv;
	if(uniform<=0.5) signdv=-1.;
	   else
		 signdv=1.;
	 
	Double_t Vr=TMath::Sqrt(Vgeant(0)*Vgeant(0)+ Vgeant(1)*Vgeant(1));
	  Double_t Dv=gRandom->Gaus(signdv*Vr,(Float_t)sigmaDv); 
    Double_t Zv=gRandom->Gaus(Vgeant(2),(Float_t)sigmaZv);
				
  //cout<<" Dv e Zv = "<<Dv<<" "<<Zv<<"\n";				
    trackITS.SetDv(Dv);  trackITS.SetZv(Zv);
    trackITS.SetsigmaDv(sigmaDv); trackITS.SetsigmaZv(sigmaZv); 
    result.SetDv(Dv);  result.SetZv(Zv);
    result.SetsigmaDv(sigmaDv); result.SetsigmaZv(sigmaZv);
    primarytrack.SetDv(Dv);  primarytrack.SetZv(Zv);
    primarytrack.SetsigmaDv(sigmaDv); primarytrack.SetsigmaZv(sigmaZv); 				 				

/////////////////////////////////////////////////////////////////////////////////////////////////  	 
	 	
    primarytrack.PrimaryTrack(rl);
    TVector  d2=primarytrack.Getd2();
    TVector  tgl2=primarytrack.Gettgl2();
    TVector  dtgl=primarytrack.Getdtgl();
    trackITS.Setd2(d2); trackITS.Settgl2(tgl2);  trackITS.Setdtgl(dtgl); 
    result.Setd2(d2); result.Settgl2(tgl2);  result.Setdtgl(dtgl); 	   
	 /*	          	 
    trackITS.SetVertex(vertex); trackITS.SetErrorVertex(ervertex);
    result.SetVertex(vertex);   result.SetErrorVertex(ervertex);   
    */                         
    Tracking(trackITS,&result,recPoints,vettid, flagvert,rl);  
	     
    // cout<<" progressive track number = "<<j<<"\r";
   // cout<<j<<"\r";
   // cout<<" progressive track number = "<<j<<"\n";
    Long_t labITS=result.GetLabel();
   // cout << " ITS track label = " << labITS << "\n"; 		    
    int lab=track->GetLabel();		    
   // cout << " TPC track label = " << lab <<"\n";
    //result.GetClusters(); getchar();  //to print the cluster matrix
	 
//propagation to vertex
	
    Double_t rbeam=3.;
     
    result.Propagation(rbeam);
       
    TMatrix *cov;
	 cov=&result.GetCMatrix();	 
    Double_t pt=TMath::Abs(result.GetPt());
    Double_t Dr=result.GetD();
    Double_t Z=result.GetZ();
    Double_t tgl=result.GetTgl();
    Double_t C=result.GetC();
    Double_t Cy=C/2.;
    Double_t Dz=Z-(tgl/Cy)*TMath::ASin(result.arga(rbeam));
	  Dz-=Vgeant(2);
	  
	 // cout<<" Dr e dz alla fine = "<<Dr<<" "<<Dz<<"\n"; getchar();
    Double_t phi=result.Getphi();
    Double_t phivertex = phi - TMath::ASin(result.argA(rbeam));
    Double_t duepi=2.*TMath::Pi();	 
    if(phivertex>duepi) phivertex-=duepi;
    if(phivertex<0.) phivertex+=duepi;
    Double_t Dtot=TMath::Sqrt(Dr*Dr+Dz*Dz);
	 
//////////////////////////////////////////////////////////////////////////////////////////    	
    Int_t NumofCluster, idmodule,idpoint;
    NumofCluster=result.GetNumClust();
    if(NumofCluster >=5)  {  	 


      AliITSiotrack outtrack;

      iotrack=&outtrack;

      iotrack->SetStatePhi(phi);
      iotrack->SetStateZ(Z);
      iotrack->SetStateD(Dr);
      iotrack->SetStateTgl(tgl);
      iotrack->SetStateC(C);
		Double_t radius=result.Getrtrack();
		iotrack->SetRadius(radius);
		Int_t charge;
		if(C>0.) charge=-1;  else charge=1;
		iotrack->SetCharge(charge);
		
		

      iotrack->SetCovMatrix(cov);         

      Double_t px=pt*TMath::Cos(phi);
      Double_t py=pt*TMath::Sin(phi);
      Double_t pz=pt*tgl;
		
      Double_t xtrack=Dr*TMath::Sin(phi);
      Double_t ytrack=Dr*TMath::Cos(phi);
      Double_t ztrack=Dz+Vgeant(2);


      iotrack->SetPx(px);
      iotrack->SetPy(py);
      iotrack->SetPz(pz);
      iotrack->SetX(xtrack);
      iotrack->SetY(ytrack);
      iotrack->SetZ(ztrack);
      iotrack->SetLabel(labITS);
		
      Int_t il;		
		for(il=0;il<6; il++){
		  iotrack->SetIdPoint(il,result.GetIdPoint(il));
		  iotrack->SetIdModule(il,result.GetIdModule(il));		
		}
      tracktree1.Fill();

   //cout<<" labITS = "<<labITS<<"\n";
	//cout<<" phi z Dr tgl C = "<<phi<<" "<<Z<<" "<<Dr<<" "<<tgl<<" "<<C<<"\n";  getchar();	   

     DataOut(kkk) = ptg; kkk++; DataOut(kkk)=labITS; kkk++; DataOut(kkk)=lab; kkk++;		

      for (il=0;il<6;il++) {
        idpoint=result.GetIdPoint(il);
        idmodule=result.GetIdModule(il);
	*(vettid[idmodule]+idpoint)=1; 
	iotrack->SetIdPoint(il,idpoint);
        iotrack->SetIdModule(il,idmodule);
      }
      
    //  cout<<"  +++++++++++++  pt e ptg = "<<pt<<" "<<ptg<<"  ++++++++++\n";
      Double_t difpt= (pt-ptg)/ptg*100.;	            	                
      DataOut(kkk)=difpt; kkk++;                                             
      Double_t lambdag=TMath::ATan(pzg/ptg);
      Double_t   lam=TMath::ATan(tgl);      
      Double_t diflam = (lam - lambdag)*1000.;
      DataOut(kkk) = diflam; kkk++;	    	  			    
      Double_t phig=TMath::ATan(pyg/pxg);
      Double_t phi=phivertex;  
      Double_t difphi = (phi - phig)*1000.;
      DataOut(kkk)=difphi; kkk++;
      DataOut(kkk)=Dtot*1.e4; kkk++;
      DataOut(kkk)=Dr*1.e4; kkk++;
      DataOut(kkk)=Dz*1.e4; kkk++;
      Int_t r;
      for (r=0; r<9; r++) { out<<DataOut(r)<<" ";}
      out<<"\n";
      kkk=0;  
		
	    
    } // end if on NumofCluster
  //gObjectTable->Print();    // stampa memoria     
  }  //  end for (int j=min_t; j<=max_t; j++)
  
  out.close();  
  
 
  static Bool_t first=kTRUE;
  static TFile *tfile;

	if(first) {
	    tfile=new TFile("itstracks.root","RECREATE");
	    //cout<<"I have opened itstracks.root file "<<endl;
	}	    
	first=kFALSE;
	tfile->cd();
	tfile->ls();

   char hname[30];
   sprintf(hname,"TreeT%d",evNumber);

  tracktree1.Write(hname);


  
	    TTree *fAli=gAlice->TreeK();
            TFile *fileAli=0;
	    
	    if (fAli) fileAli =fAli->GetCurrentFile();
	    fileAli->cd();
     
  ////////////////////////////////////////////////////////////////////////////////////////////////

  printf("delete vectors\n");
  if(np) delete [] np;
  if(vettid) delete [] vettid;
  
}
