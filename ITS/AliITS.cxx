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
Revision 1.57  2001/07/24 14:26:11  mariana
Introduce the function Digits2Reco() and write the defaults for simulation and reconstruction

Revision 1.56  2001/07/05 12:49:49  mariana
Temporary patches required by root.v3.01.05

Revision 1.55  2001/06/14 14:59:00  barbera
Tracking V1 decoupled from AliITS

Revision 1.54  2001/05/31 20:37:56  barbera
Bari/Salerno model set as defaault SPD simulation

Revision 1.53  2001/05/31 18:52:24 barbera 
Bari model becomes the default

Revision 1.53  2001/05/30 07:52:24  hristov
TPC and CONTAINERS included in the search path

Revision 1.52  2001/05/30 06:04:58  hristov
Changes made to be consitant with changes in TPC tracking classes (B.Nilsen)

Revision 1.51  2001/05/16 14:57:15  alibrary
New files for folders and Stack

Revision 1.50  2001/05/11 09:15:21  barbera
Corrected to make fast point creation working with PPR geometry

Revision 1.49  2001/05/11 07:37:49  hristov
Legacy lines commented

Revision 1.48  2001/05/10 18:14:25  barbera
A typo corrected

Revision 1.47  2001/05/10 17:55:59  barbera
Modified to create rec points also for PPR geometries

Revision 1.46  2001/05/10 00:05:28  nilsen
Allowed for HitsToDigits function to work with versions 5, 7, 8, and 9. This
should probably be cleaned up to only check to make sure that fITSgeom has
been properly defined.

Revision 1.45  2001/05/01 22:35:48  nilsen
Remove/commented a number of cout<< statements. and made change needed by
SSD code.

Revision 1.44  2001/04/26 22:44:01  nilsen
Removed dependence on layer 5/6 in AliITS::HitsToDigits. This will be
done properly in AliITSv???.cxx via SetDefaults.

Revision 1.43  2001/04/26 13:22:52  barbera
TMatrix and TVector elimininated to speed up the code

Revision 1.42  2001/04/25 21:55:12  barbera
Updated version to be compatible with actual verion of STEER and TPC

Revision 1.41  2001/04/21 15:16:51  barbera
Updated with the new SSD reconstruction code

Revision 1.40  2001/03/17 15:07:06  mariana
Update SDD response parameters

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



#include "AliMC.h"
#include "AliRun.h"
#include "AliHeader.h"

#include "AliITS.h"
#include "AliITSDetType.h"
#include "AliITSresponseSPD.h"
#include "AliITSresponseSDD.h"
#include "AliITSresponseSSD.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSsimulationSPD.h"
#include "AliITSsimulationSDD.h"
#include "AliITSsimulationSSD.h"
#include "AliITSClusterFinderSPD.h"
#include "AliITSClusterFinderSDD.h"
#include "AliITSClusterFinderSSD.h"

#include "AliITShit.h"
#include "AliITSgeom.h"
#include "AliITSdigit.h"
#include "AliITSmodule.h"
#include "AliITSRecPoint.h"
#include "AliITSRawCluster.h"


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



  fRecPoints=new TClonesArray("AliITSRecPoint",1000);
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
    fDetTypes->AddAt(new AliITSDetType(),i); 
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

  if(fDebug) printf("%s: SetDefaults\n",ClassName());

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
  // sets the default simulation


  AliITSDetType *iDetType;
  iDetType=DetType(0);
  if (!iDetType->GetSimulationModel()) {
      AliITSsegmentation *seg0=
                    (AliITSsegmentation*)iDetType->GetSegmentationModel();
      AliITSresponse *res0 = (AliITSresponse*)iDetType->GetResponseModel();
      AliITSsimulationSPD *sim0=new AliITSsimulationSPD(seg0,res0);
      SetSimulationModel(0,sim0);
  }
  iDetType=DetType(1);
  if (!iDetType->GetSimulationModel()) {
      AliITSsegmentation *seg1=
                    (AliITSsegmentation*)iDetType->GetSegmentationModel();
      AliITSresponse *res1 = (AliITSresponse*)iDetType->GetResponseModel();
      AliITSsimulationSDD *sim1=new AliITSsimulationSDD(seg1,res1);
      SetSimulationModel(1,sim1);
  }
  iDetType=DetType(2);
  if (!iDetType->GetSimulationModel()) {
      AliITSsegmentation *seg2=
                    (AliITSsegmentation*)iDetType->GetSegmentationModel();
      AliITSresponse *res2 = (AliITSresponse*)iDetType->GetResponseModel();
      AliITSsimulationSSD *sim2=new AliITSsimulationSSD(seg2,res2);
      SetSimulationModel(2,sim2);
  }


}
//_____________________________________________________________________________
void AliITS::SetDefaultClusterFinders()
{
  // sets the default cluster finders

  MakeTreeC();
  AliITSDetType *iDetType;
  iDetType=DetType(0);
  if (!iDetType->GetReconstructionModel()) {
      AliITSsegmentation *seg0=
                   (AliITSsegmentation*)iDetType->GetSegmentationModel();
      TClonesArray *dig0=DigitsAddress(0);
      TClonesArray *recp0=ClustersAddress(0);
      AliITSClusterFinderSPD *rec0=new AliITSClusterFinderSPD(seg0,dig0,recp0);
      SetReconstructionModel(0,rec0);
  }
  iDetType=DetType(1);
  if (!iDetType->GetReconstructionModel()) {
      AliITSsegmentation *seg1=
                     (AliITSsegmentation*)iDetType->GetSegmentationModel();
      AliITSresponse *res1 = (AliITSresponse*)iDetType->GetResponseModel();
      TClonesArray *dig1=DigitsAddress(1);
      TClonesArray *recp1=ClustersAddress(1);
      AliITSClusterFinderSDD *rec1=
                           new AliITSClusterFinderSDD(seg1,res1,dig1,recp1);

      SetReconstructionModel(1,rec1);
  }
  iDetType=DetType(2);
  if (!iDetType->GetReconstructionModel()) {
      AliITSsegmentation *seg2=
                    (AliITSsegmentation*)iDetType->GetSegmentationModel();

      TClonesArray *dig2=DigitsAddress(2);
      AliITSClusterFinderSSD *rec2= new AliITSClusterFinderSSD(seg2,dig2);
      SetReconstructionModel(2,rec2);
  }

}
//_____________________________________________________________________________

void AliITS::MakeTreeC(Option_t *option)
{
  // create a separate tree to store the clusters

//  cout << "AliITS::MakeTreeC" << endl;

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
        fCtype->AddAt(new TClonesArray(clclass,1000),i);
        if (kNTYPES==3) sprintf(branchname,"%sClusters%s",GetName(),det[i]);
	else  sprintf(branchname,"%sClusters%d",GetName(),i+1);
	if (fCtype   && fTreeC) {
	   TreeC()->Branch(branchname,&((*fCtype)[i]), buffersize);
//	   cout << "Making Branch " << branchname;
//	   cout << " for Clusters of detector type " << i+1 << endl;
	}	
     }

}

//_____________________________________________________________________________
void AliITS::GetTreeC(Int_t event)
{

//  cout << "AliITS::GetTreeC" << endl;

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
	char digclass[40];
	char clclass[40];
	for (i=0; i<kNTYPES; i++) {

	   AliITSDetType *iDetType=DetType(i); 
	   iDetType->GetClassNames(digclass,clclass);
	   // clusters
	   if(!(*fCtype)[i]) fCtype->AddAt(new TClonesArray(clclass,1000),i); 
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
void AliITS::MakeBranch(Option_t* option, const char *file)
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
       DetType(i)->GetClassNames(digclass,clclass);
       // digits
       if(!((*fDtype)[i])) fDtype->AddAt(new TClonesArray(digclass,1000),i);
       else ResetDigits(i);
   }

   for (i=0; i<kNTYPES ;i++) {
      if (kNTYPES==3) sprintf(branchname,"%sDigits%s",GetName(),det[i]);
      else  sprintf(branchname,"%sDigits%d",GetName(),i+1);      
      if (fDtype && gAlice->TreeD()) {
        MakeBranchInTree(gAlice->TreeD(), 
                         branchname, &((*fDtype)[i]), buffersize, file);
//	    cout << "Making Branch " << branchname;
//	    cout << " for digits of type "<< i+1 << endl;
      }	
    }
  }

  if (cR) {
  //
  // only one branch for rec points for all detector types
  //
    sprintf(branchname,"%sRecPoints",GetName());
 
    //if(!fRecPoints) fRecPoints=new TClonesArray("AliITSRecPoint",1000);

    if (fRecPoints && gAlice->TreeR()) {
      MakeBranchInTree(gAlice->TreeR(), 
                               branchname, &fRecPoints, buffersize, file) ;
//      cout << "Making Branch " << branchname;
//      cout << " for reconstructed space points" << endl;
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

  char digclass[40];
  char clclass[40];

  Int_t i;
  if (treeD) {
      for (i=0; i<kNTYPES; i++) {
	DetType(i)->GetClassNames(digclass,clclass);
	  // digits
        if(!((*fDtype)[i])) fDtype->AddAt(new TClonesArray(digclass,1000),i);
	else ResetDigits(i);

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
//	    cout<<"filename "<<filename<<endl;
	    file=new TFile(filename);
//	    cout<<"I have opened "<<filename<<" file "<<endl;
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

  cerr<<"Digitizing ITS...\n";

  TStopwatch timer;
  timer.Start();
  AliHeader *header=gAlice->GetHeader();
  HitsToDigits(header->GetEvent(),0,-1," ","All"," ");
  timer.Stop(); timer.Print();

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
   if(ver!=5 && ver!=7 && ver!=8 && ver!=9) return; 

   const char *all = strstr(opt,"All");
   const char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),strstr(opt,"SSD")};

   static Bool_t setDef=kTRUE;
   if (setDef) SetDefaultSimulation();
   setDef=kFALSE;


//   cout<<" 1 AliITS "<<endl;
   Int_t nmodules;
   InitModules(size,nmodules); 
//   cout<<" 2 AliITS "<<endl;
   FillModules(evNumber,bgrev,nmodules,option,filename);
//   cout<<" 3 AliITS "<<endl;

   //TBranch *branch;
   AliITSsimulation* sim;
   //TObjArray *branches=gAlice->TreeD()->GetListOfBranches();
   AliITSgeom *geom = GetITSgeom();

   Int_t id,module;
//   Int_t lay, lad, detect;
   Int_t first,last;
   for (id=0;id<kNTYPES;id++) {
        if (!all && !det[id]) continue;
	//branch = (TBranch*)branches->UncheckedAt(id);
	AliITSDetType *iDetType=DetType(id); 
	sim = (AliITSsimulation*)iDetType->GetSimulationModel();
	if(geom) {
	  first = geom->GetStartDet(id);
	  last = geom->GetLastDet(id);
	} else first=last=0;
	printf("first module - last module %d %d\n",first,last);
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

//   Int_t nentries=(Int_t)
   gAlice->TreeD()->GetEntries();
//   cout << "nentries in TreeD" << nentries << endl;

   char hname[30];
   sprintf(hname,"TreeD%d",evNumber);
   gAlice->TreeD()->Write(hname,TObject::kOverwrite);
   // reset tree
   gAlice->TreeD()->Reset();

}


//_____________________________________________________________________________
void AliITS::Digits2Reco()
{
  // find clusters and reconstruct space points

  AliHeader *header=gAlice->GetHeader();
  printf("header->GetEvent() %d\n",header->GetEvent());
  DigitsToRecPoints(header->GetEvent(),0,"All");

}
//____________________________________________________________________________
void AliITS::DigitsToRecPoints(Int_t evNumber,Int_t lastentry,Option_t *opt)
{
  // cluster finding and reconstruction of space points
  
   // the condition below will disappear when the geom class will be
   // initialised for all versions - for the moment it is only for v5 !
   // 7 is the SDD beam test version  
   Int_t ver = this->IsVersion(); 
   if(ver!=5 && ver!=8 && ver!=9) return;

   const char *all = strstr(opt,"All");
   const char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),strstr(opt,"SSD")};

   static Bool_t setRec=kTRUE;
   if (setRec) SetDefaultClusterFinders();
   setRec=kFALSE;


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
        TClonesArray *itsDigits  = this->DigitsAddress(id);
        Int_t first,last;
	if(geom) {
	  first = geom->GetStartDet(id);
	  last = geom->GetLastDet(id);
	} else first=last=0;
	printf("first module - last module %d %d\n",first,last);
	for(module=first;module<=last;module++) {
              this->ResetDigits();
              if (all) gAlice->TreeD()->GetEvent(lastentry+module);
	      else gAlice->TreeD()->GetEvent(lastentry+(module-first));
	      Int_t ndigits = itsDigits->GetEntriesFast();
	      if (ndigits) rec->FindRawClusters(module);
	      gAlice->TreeR()->Fill(); 
	      ResetRecPoints();
	      treeC->Fill();
              ResetClusters();
	      // try and fill only the branch 
	      //branch->Fill();
	      //ResetRecPoints(id);
	} // loop over modules
   } // loop over detector types


//   Int_t nentries=(Int_t)
   gAlice->TreeR()->GetEntries();
//   Int_t ncentries=(Int_t)
   treeC->GetEntries();
//   cout << " nentries ncentries " << nentries << ncentries <<  endl;

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
   if(ver!=5 && ver!=8 && ver!=9) return;


   const char *all = strstr(opt,"All");
   const char *det[3] ={strstr(opt,"SPD"),strstr(opt,"SDD"),strstr(opt,"SSD")};

   Int_t nmodules;
   InitModules(size,nmodules);
   FillModules(evNumber,bgrev,nmodules,option,filename);


   AliITSsimulation* sim;
   AliITSgeom *geom = GetITSgeom();
/* MI change 
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
   */

   Int_t id,module;
   for (id=0;id<kNTYPES;id++) {
        if (!all && !det[id]) continue;
	AliITSDetType *iDetType=DetType(id); 
	sim = (AliITSsimulation*)iDetType->GetSimulationModel();
	if (!sim) {
           Error("HitsToFastPoints",
		 "The simulation class was not instantiated!");
           exit(1);
	   // or SetDefaultSimulation();
	}

        Int_t first,last;
        if(geom) {
	  first = geom->GetStartDet(id);
	  last = geom->GetLastDet(id);
        } else first=last=0;
        printf("first module - last module %d %d\n",first,last);
	for(module=first;module<=last;module++) {
	    AliITSmodule *mod = (AliITSmodule *)fITSmodules->At(module);
        // sim->CreateFastRecPoints(mod,module,random);
        sim->CreateFastRecPoints(mod,module,gRandom); //MI change
	   
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

   //   delete [] random; //MI change 


}

