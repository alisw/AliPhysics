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
Revision 1.9.2.8  2000/06/12 18:05:59  barbera
fixed posible compilation errors on HP unix

Revision 1.9.2.7  2000/06/11 20:20:18  barbera
New AliITS base clase  for the new structure.

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
 
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TROOT.h>
#include <TObjectTable.h>



#include "AliRun.h"
#include "AliITS.h"
#include "AliITSMap.h"
#include "AliITSClusterFinder.h"
#include "AliITSsimulation.h"
#include "AliITSsimulationFastPoints.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSresponseSPD.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSresponseSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSresponseSSD.h"
//#include "AliITStrack.h"


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
  fIdN        = 0;
  fIdName     = 0;
  fIdSens     = 0;
  fITSmodules = 0;
  //
  fDetTypes   = 0;
  SetNDetTypes();
  //
  fDtype  = 0;
  fNdtype = 0;
  fCtype  = 0;
  fNctype = 0;
  fRecPoints = 0;
  fNRecPoints = 0;
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

  SetNDetTypes();

  fNdtype = 0;
  fDtype = 0;
  fCtype  = 0;
  fNctype = 0;

  fRecPoints = 0;
  fNRecPoints = 0;


  fITSmodules = 0; 

  fIshunt     = 0;
  fEuclidOut  = 0;
  fIdN        = 0;
  fIdName     = 0;
  fIdSens     = 0;
 
  fDetTypes = new TObjArray(fNDetTypes);  

  Int_t i;
  for(i=0;i<fNDetTypes;i++) {
    (*fDetTypes)[i]=new AliITSDetType(); 
   }
  //

  SetMarkerColor(kRed);

  fITSgeom=0;
}
//___________________________________________________________________________
AliITS::AliITS(AliITS &source){
  if(this==&source) return;
  printf("Error: You are not allowed to make a copy of the AliITS\n");
  exit(1);
}
//____________________________________________________________________________
AliITS& AliITS::operator=(AliITS &source){
  if(this==&source) return *this;
  printf("Error: You are not allowed to make a copy of the AliITS\n");
  exit(1);
}
//____________________________________________________________________________
void AliITS::ClearModules(){
  //clear the modules TObjArray
  Int_t i;

  if(fITSmodules!=0) {
	Int_t indSPD = fITSgeom->GetModuleIndex(2,fITSgeom->GetNladders(2),
						fITSgeom->GetNdetectors(2));
	Int_t indSDD = fITSgeom->GetModuleIndex(4,fITSgeom->GetNladders(4),
						fITSgeom->GetNdetectors(4));
      for(i=0;i<fITSmodules->GetEntriesFast();i++){
	    if(i<indSPD)
		delete (AliITSmodule *) fITSmodules->At(i);
	    else if(i<indSDD)
		delete (AliITSmodule *) fITSmodules->At(i);
	    else
		delete (AliITSmodule *) fITSmodules->At(i);
      } // end for i
  }// end if fITSmodules!=0

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
  if(fIdName!=0) delete[] fIdName;
  if(fIdSens!=0) delete[] fIdSens;
  if(fITSmodules!=0) {
      this->ClearModules();
      delete fITSmodules;
  }// end if fITSmodules!=0

  //
  Int_t i;
  if(fDtype) {
    for(i=0;i<fNDetTypes;i++) {
      delete (*fDtype)[i];
      fNdtype[i]=0;
    }
  }

  for(i=0;i<fNDetTypes;i++) {
      delete (*fCtype)[i];
      fNctype[i]=0;
  }

  //

  if (fDetTypes) {
    fDetTypes->Delete();
    delete fDetTypes;
  }

  if (fTreeC) delete fTreeC;
  
}

//___________________________________________
AliITSDetType* AliITS::DetType(Int_t id)
{
  //return pointer to id detector type
    return ((AliITSDetType*) (*fDetTypes)[id]);

}
//___________________________________________
void AliITS::SetClasses(Int_t id, TString digit, TString cluster)
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
void AliITS::AddDigit(Int_t id, AliITSdigit *d) 
{

  // add a simulated digit

  // should have ctors of type AliITSdigitSDD(const AliITSdigitSDD &)

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
void AliITS::AddDigit(Int_t id,Float_t phys,Int_t *digits,Int_t *tracks,Float_t *charges){

  // add a simulated digit to the list

  TClonesArray &ldigits = *((TClonesArray*)(*fDtype)[id]);
  switch(id)
  {
  case 0:
     new(ldigits[fNdtype[id]++]) AliITSdigitSPD(digits,tracks);
     break;
  case 1:
     new(ldigits[fNdtype[id]++]) AliITSdigitSDD(phys,digits,tracks,charges);
     break;
  case 2:
     new(ldigits[fNdtype[id]++]) AliITSdigitSSD(digits,tracks);
     break;
  }
 
}

//_____________________________________________________________________________
void AliITS::AddCluster(Int_t id, AliITSRawCluster *c) 
{

  // add a cluster to the list

  // should have ctors of type AliITSRawClusterSDD(const AliITSRawClusterSDD &)

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
    // Reset number of digits and the digits array for thE ITS detector
    //

    if (!fDtype) return;

    Int_t i;
    for(i=0;i<fNDetTypes;i++ ) {
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
    for(i=0;i<fNDetTypes;i++ ) {
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
    fNRecPoints = 0;
    if (fRecPoints) fRecPoints->Clear();

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


  SetDefaults();

  Int_t i;

  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" ITS_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  //
  for(i=0;i<fIdN;i++) fIdSens[i] = gMC->VolId(fIdName[i]);
  //
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

//_____________________________________________________________________________
void AliITS::SetDefaults()
{
  // sets the default segmentation, response, digit and raw cluster classes


  AliITSDetType *iDetType;

  //SPD 

  AliITSsegmentationSPD *seg0=new AliITSsegmentationSPD(fITSgeom);
  AliITSresponseSPD *resp0=new AliITSresponseSPD();
  iDetType=DetType(0); 
  if (!iDetType->GetSegmentationModel()) SetSegmentationModel(0,seg0); 
  if (!iDetType->GetResponseModel()) SetResponseModel(0,resp0); 
  // set digit and raw cluster classes to be used
  const char *kData=resp0->DataType();
  if (strstr(kData,"real")) {
      iDetType->ClassNames("AliITSdigit","AliITSRawClusterSPD");
  } else iDetType->ClassNames("AliITSdigitSPD","AliITSRawClusterSPD");

  // SDD					  //
  AliITSresponseSDD *resp1=new AliITSresponseSDD();
  AliITSsegmentationSDD *seg1=new AliITSsegmentationSDD(fITSgeom,resp1);
  iDetType=DetType(1); 
  if (!iDetType->GetSegmentationModel()) SetSegmentationModel(1,seg1); 
  if (!iDetType->GetResponseModel()) SetResponseModel(1,resp1); 
  kData=resp1->DataType();
  Option_t *opt=resp1->ZeroSuppOption();
  if ((!strstr(opt,"2D")) && (!strstr(opt,"1D")) || strstr(kData,"real") ) {
      iDetType->ClassNames("AliITSdigit","AliITSRawClusterSDD");
  } else iDetType->ClassNames("AliITSdigitSDD","AliITSRawClusterSDD");

  // SSD
  AliITSsegmentationSSD *seg2=new AliITSsegmentationSSD(fITSgeom);
  AliITSresponseSSD *resp2=new AliITSresponseSSD();
  iDetType=DetType(2); 
  if (!iDetType->GetSegmentationModel()) SetSegmentationModel(2,seg2); 
  if (!iDetType->GetResponseModel()) SetResponseModel(2,resp2); 
  kData=resp2->DataType();
  if (strstr(kData,"real")) {
      iDetType->ClassNames("AliITSdigit","AliITSRawClusterSSD");
  } else iDetType->ClassNames("AliITSdigitSSD","AliITSRawClusterSSD");

  if (fNDetTypes>3) {
    Warning("SetDefaults","Only the three basic detector types are initialised!");
  } 

}


//_____________________________________________________________________________

void AliITS::MakeTreeC(Option_t *option)
{
  // create a separate tree to store the clusters

     char *optC = strstr(option,"C");
     if (optC && !fTreeC) fTreeC = new TTree("TC","Clusters in ITS");

     Int_t buffersize = 4000;
     char branchname[30];

     char *det[3] = {"SPD","SDD","SSD"};

     // one branch for Clusters per type of detector
     Int_t i;
     for(i=0; i<fNDetTypes ;i++) {
        if (fNDetTypes==3) sprintf(branchname,"%sClusters%s",GetName(),det[i]);
	else  sprintf(branchname,"%sClusters%d",GetName(),i+1);
	if (fCtype   && fTreeC) {
	   TreeC()->Branch(branchname,&((*fCtype)[i]), buffersize);
	   printf("Making Branch %s for Clusters of detector type %d\n",branchname,i+1);
	}	
     }

}

//_____________________________________________________________________________
void AliITS::GetTreeC(Int_t event)
{

  // get the clusters tree for this event and set the branch address
    char treeName[20];
    char branchname[30];

    char *det[3] = {"SPD","SDD","SSD"};

    ResetClusters();
    if (fTreeC) {
	  delete fTreeC;
    }

    sprintf(treeName,"TreeC%d",event);
    fTreeC = (TTree*)gDirectory->Get(treeName);


    TBranch *branch;
    if (fTreeC) {
        Int_t i;
	for(i=0; i<fNDetTypes; i++) {
	    if (fNDetTypes==3) sprintf(branchname,"%sClusters%s",GetName(),det[i]);
	    else  sprintf(branchname,"%sClusters%d",GetName(),i+1);
	    if (fCtype) {
		branch = fTreeC->GetBranch(branchname);
                if (branch) branch->SetAddress(&((*fCtype)[i]));
	    }
	}
    } else {
	printf("ERROR: cannot find Clusters Tree for event:%d\n",event);
    }

}
//_____________________________________________________________________________
void AliITS::MakeBranch(Option_t* option){
  //
  // Creates Tree branches for the ITS.
  //


  Int_t buffersize = 4000;
  char branchname[30];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option);


// one branch for digits per type of detector
  
   char *det[3] = {"SPD","SDD","SSD"};

   fNdtype = new Int_t[fNDetTypes];
   fDtype = new TObjArray(fNDetTypes);

   fNctype = new Int_t[fNDetTypes];
   fCtype = new TObjArray(fNDetTypes);

   const char *kDigclass;
   const char *kClclass;

   Int_t i;
   for(i=0; i<fNDetTypes ;i++) {
       AliITSDetType *iDetType=DetType(i); 
       iDetType->GetClassNames(kDigclass,kClclass);
       //printf("i, digclass, recclass %d %s %s\n",i,kDigclass,kClclass); 
       // digits
       (*fDtype)[i] = new TClonesArray(kDigclass,100); 
       fNdtype[i]=0;
       // clusters
       (*fCtype)[i] = new TClonesArray(kClclass,100); 
       fNctype[i]=0;
   }


  for(i=0; i<fNDetTypes ;i++) {
      if (fNDetTypes==3) sprintf(branchname,"%sDigits%s",GetName(),det[i]);
      else  sprintf(branchname,"%sDigits%d",GetName(),i+1);
      
      if (fDtype   && gAlice->TreeD()) {
	  gAlice->TreeD()->Branch(branchname,&((*fDtype)[i]), buffersize);
	  printf("Making Branch %s for digits of type %d\n",branchname,i+1);
      }	
  }

  // only one branch for rec points for all detector types
  sprintf(branchname,"%sRecPoints",GetName());

  fRecPoints=new TClonesArray("AliITSRecPoint",1000);

  if (fRecPoints && gAlice->TreeR()) {
    gAlice->TreeR()->Branch(branchname,&fRecPoints, buffersize);
    printf("Making Branch %s for reconstructed space points\n",branchname);
  }	


}

//___________________________________________
void AliITS::SetTreeAddress()
{

  // Set branch address for the Trees.

  char branchname[30];
  AliDetector::SetTreeAddress();

  char *det[3] = {"SPD","SDD","SSD"};

  TBranch *branch;
  TTree *treeD = gAlice->TreeD();
  TTree *treeR = gAlice->TreeR();

  Int_t i;
  if (treeD) {
      for(i=0; i<fNDetTypes; i++) {
	  if (fNDetTypes==3) sprintf(branchname,"%sDigits%s",GetName(),det[i]);
	  else  sprintf(branchname,"%sDigits%d",GetName(),i+1);
	  if (fDtype) {
	      branch = treeD->GetBranch(branchname);
	      if (branch) branch->SetAddress(&((*fDtype)[i]));
	  }
      }
  }

 
  if (treeR) {
    sprintf(branchname,"%sRecPoints",GetName());
    if (fRecPoints) {
      branch = treeR->GetBranch(branchname);
      if (branch) branch->SetAddress(&fRecPoints);
    }
  }
  

}

//____________________________________________________________________________
void AliITS::InitModules(Int_t size,Int_t &nmodules){

  //initialize the modules array

    Int_t nl,indexMAX,index;
    Int_t indSPD,indSDD;

    if(size<=0){ // default to using data stored in AliITSgeom
	if(fITSgeom==0) {
	    printf("Error in AliITS::InitModule fITSgeom not defined\n");
	    return;
	} // end if fITSgeom==0
	nl = fITSgeom->GetNlayers();
	indexMAX = fITSgeom->GetModuleIndex(nl,fITSgeom->GetNladders(nl),
					    fITSgeom->GetNdetectors(nl))+1;
	nmodules = indexMAX;
	fITSmodules = new TObjArray(indexMAX);
	indSPD = fITSgeom->GetModuleIndex(2,fITSgeom->GetNladders(2),
					    fITSgeom->GetNdetectors(2));
	indSDD = fITSgeom->GetModuleIndex(4,fITSgeom->GetNladders(4),
					    fITSgeom->GetNdetectors(4));
	for(index=0;index<indexMAX;index++){
	    if(index<=indSPD)
		fITSmodules->AddAt( new AliITSmodule(index),index);
	    else if(index<=indSDD)
		fITSmodules->AddAt( new AliITSmodule(index),index);
	    else
	        fITSmodules->AddAt( new AliITSmodule(index),index);
	} // end for index
    }else{
	fITSmodules = new TObjArray(size);
        nmodules = size;
    } // end i size<=0
}

//____________________________________________________________________________
void AliITS::FillModules(Int_t evnt,Int_t bgrev,Int_t lastev,Int_t nmodules,Option_t *option,Text_t *filename){

  // fill the modules with the sorted by module hits; add hits from background
  // if option=Add

    static TTree *trH1;                 //Tree with background hits
    static TClonesArray *fHits2;        //List of hits for one track only

    static Bool_t first=kTRUE;
    static TFile *file;
    char *add = strstr(option,"Add");


    if (add ) {
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
	    printf("ERROR: cannot find Hits Tree for event:%d\n",bgrev);
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

    // store temporarily coordinates of signal particles - see where delete
    static Int_t *signal;
    if(!signal) signal=new Int_t[nmodules];
    memset(signal,0,sizeof(int)*nmodules);
    Float_t xhit[nmodules][4];
    Float_t yhit[nmodules][4];

    Int_t npart = gAlice->GetEvent(evnt);
    if(npart<=0) return;
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
        if (!nhits) continue;
       	// cout << nhits << " hits in track " << t << endl;
	for(h=0; h<nhits; h++){
	    itsHit = (AliITShit *)itsHits->UncheckedAt(h);
	    itsHit->GetDetectorID(lay,lad,det);
	    index = fITSgeom->GetModuleIndex(lay,lad,det);
	    mod = this->GetModule(index);
	    if (add) {
	        xhit[index][signal[index]]=itsHit->fX;
	        yhit[index][signal[index]]=itsHit->fY;
	        signal[index]++;
                if (signal[index] >4) 
                      printf("index,nsignal %d %d\n",index,signal[index]);
	    }
	    if(lay == 1 || lay == 2)
		mod->AddHit((AliITShit *) itsHit,t,h);
	    else if(lay == 3 || lay == 4)
		    mod->AddHit((AliITShit *) itsHit,t,h);
	    else if(lay == 5 || lay ==6)
		mod->AddHit((AliITShit *)itsHit,t,h);
	    else
		mod->AddHit(itsHit,t,h);

	} // end loop over hits 
    } // end loop over tracks

    // open the file with background
    
    if (add ) {
          Int_t track,i,isig;
          ntracks =(Int_t)trH1->GetEntries();
	    //printf("background - ntracks1 %d\n",ntracks);
	    //printf("background - Start loop over tracks \n");     
            //   Loop over tracks

	    for(track=0; track<ntracks; track++) {

		if (fHits2)       fHits2->Clear();
		trH1->GetEvent(track);
                //   Loop over hits
		for(i=0;i<fHits2->GetEntriesFast();++i) {

		    itsHit=(AliITShit*) (*fHits2)[i];
		    itsHit->GetDetectorID(lay,lad,det);
		    index = fITSgeom->GetModuleIndex(lay,lad,det);
		    mod = this->GetModule(index);

                    Float_t xbgr=itsHit->fX;
		    Float_t ybgr=itsHit->fY;
		    Float_t ebgr=itsHit->GetIonization();
		    Bool_t cond=kFALSE;
		    
		    for(isig =0; isig < signal[index]; isig++) {
			Float_t dist= 
                            (xbgr-xhit[index][isig])*(xbgr-xhit[index][isig])
			    +(ybgr-yhit[index][isig])*(ybgr-yhit[index][isig]);
			if (dist<0.2&& ebgr!=0) cond=kTRUE; // check this number for ITS!
		    }
		    if (!cond) continue;
		    if(lay == 1 || lay == 2)
		       mod->AddHit((AliITShit *) itsHit,track,i);
		    else if(lay == 3 || lay == 4)
		            mod->AddHit((AliITShit *) itsHit,track,i);
		    else if(lay == 5 || lay ==6)
		            mod->AddHit((AliITShit *)itsHit,track,i);
		    else
		       mod->AddHit(itsHit,track,i);


	       }  // end loop over hits
	    } // end loop over tracks

	    TTree *fAli=gAlice->TreeK();
            TFile *fileAli=0;
	    
	    if (fAli) fileAli =fAli->GetCurrentFile();
            //printf("fAli, file %p %p\n",fAli,file);
	    file->cd();

    } // end if add

    //if (evnt==lastev) {delete [] signal; delete signal;}

    //gObjectTable->Print();

}


//____________________________________________________________________________
void AliITS::HitsToDigits(Int_t evNumber,Int_t bgrev,Int_t lastev,Int_t size,
Option_t *option,Option_t *opt,Text_t *filename)
{
    // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !
  
   char *all = strstr(opt,"All");
   char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),strstr(opt,"SSD")};
   //printf("Det 1 2 3 %s %s %s \n",det[0],det[1],det[2]);

   Int_t nmodules;
   InitModules(size,nmodules); 
   FillModules(evNumber,bgrev,lastev,nmodules,option,filename);
   //printf("nmodules %d\n",nmodules);

   TBranch *branch;
   AliITSsimulation* sim;

   TObjArray *branches=gAlice->TreeD()->GetListOfBranches();
   AliITSgeom *geom = GetITSgeom();

   Int_t id,module;
   for(id=0;id<3;id++) {
        //printf("id %d All %s  det[id] %s \n",id,all,det[id]);
        if (!all && !det[id]) continue;
	branch = (TBranch*)branches->UncheckedAt(id);
	AliITSDetType *iDetType=DetType(id); 
	sim = (AliITSsimulation*)iDetType->GetSimulationModel();
	if (!sim) {
           Error("HitsToDigits","The simulation class was not instantiated!");
           exit(1);
	   // or SetDefaultSimulation(id,iDetType*);
	}
	Int_t first = geom->GetStartDet(id);
	Int_t last = geom->GetLastDet(id);
	//printf("det type %d first, last %d %d \n",id,first,last);
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

   //Int_t nentries=(Int_t)gAlice->TreeD()->GetEntries();

   char hname[30];
   sprintf(hname,"TreeD%d",evNumber);
   gAlice->TreeD()->Write(hname);
   // reset tree
   gAlice->TreeD()->Reset();

}


//____________________________________________________________________________
void AliITS::DigitsToRecPoints(Int_t evNumber,Int_t lastentry,Option_t *opt)
{
  
   char *all = strstr(opt,"All");
   char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),strstr(opt,"SSD")};

   static Bool_t first=kTRUE;
   if (first) {
       MakeTreeC("C");
       first=kFALSE;
   }
   TTree *iTC=TreeC();
 
   TBranch *branch;
   AliITSClusterFinder* rec;

   TObjArray *branches=gAlice->TreeR()->GetListOfBranches();
   AliITSgeom *geom = GetITSgeom();

   Int_t id,module;
   for(id=0;id<fNDetTypes;id++) {
        if (!all && !det[id]) continue;
	branch = (TBranch*)branches->UncheckedAt(id);
	AliITSDetType *iDetType=DetType(id); 
	rec = (AliITSClusterFinder*)iDetType->GetReconstructionModel();
	if (!rec) {
           Error("DigitsToRecPoints","The cluster finder class was not instantiated!");
           exit(1);
	   // or SetDefaultClusterFinder(id,iDetType*);
	}
        TClonesArray *itsDigits  = this->DigitsAddress(id);

	Int_t first = geom->GetStartDet(id);
	Int_t last = geom->GetLastDet(id);
	//printf("det type %d first, last %d %d \n",id,first,last);
	for(module=first;module<=last;module++) {
	      //printf("AliITS: module=%d\n",module);
              this->ResetDigits();
              if (all) gAlice->TreeD()->GetEvent(lastentry+module);
	      else gAlice->TreeD()->GetEvent(lastentry+(module-first));
	      Int_t ndigits = itsDigits->GetEntriesFast();
	      if (ndigits) rec->FindRawClusters();
	      gAlice->TreeR()->Fill(); 
	      ResetRecPoints();
	      iTC->Fill();
              ResetClusters();
	      // try and fill only the branch 
	      //branch->Fill();
	      //ResetRecPoints(id);
	} // loop over modules
   } // loop over detector types


   //Int_t nentries=(Int_t)gAlice->TreeR()->GetEntries();

   //Int_t ncentries=(Int_t)TC->GetEntries();

   char hname[30];
   sprintf(hname,"TreeR%d",evNumber);
   gAlice->TreeR()->Write(hname);
   // reset tree
   gAlice->TreeR()->Reset();

   sprintf(hname,"TreeC%d",evNumber);
   iTC->Write(hname);
   iTC->Reset();
}


//____________________________________________________________________________
void AliITS::HitsToFastRecPoints(Int_t evNumber,Int_t bgrev,Int_t lastev,Int_t size,
Option_t *option,Option_t *opt,Text_t *filename)
{
    // keep galice.root for signal and name differently the file for 
    // background when add! otherwise the track info for signal will be lost !
  
   char *all = strstr(opt,"All");
   char *det[3] = {strstr(opt,"SPD"),strstr(opt,"SDD"),strstr(opt,"SSD")};

   Int_t nmodules;
   InitModules(size,nmodules); 
   FillModules(evNumber,bgrev,lastev,nmodules,option,filename);

   static AliITSsimulationFastPoints* sim=0;
   if (!sim) sim = new AliITSsimulationFastPoints(); 


   AliITSgeom *geom = GetITSgeom();

   Int_t id,module;
   for(id=0;id<3;id++) {
        if (!all && !det[id]) continue;
	Int_t first = geom->GetStartDet(id);
	Int_t last = geom->GetLastDet(id);
	for(module=first;module<=last;module++) {
	    AliITSmodule *mod = (AliITSmodule *)fITSmodules->At(module);
	    sim->CreateFastRecPoints(mod);
	    gAlice->TreeR()->Fill(); 
	    ResetRecPoints();
	} // loop over modules
   } // loop over detector types


   ClearModules();

   //Int_t nentries=(Int_t)gAlice->TreeR()->GetEntries();

   char hname[30];
   sprintf(hname,"TreeR%d",evNumber);
   gAlice->TreeR()->Write(hname);
   // reset tree
   gAlice->TreeR()->Reset();

}

//____________________________________________________________________________
void AliITS::Streamer(TBuffer &R__b){
   // Stream an object of class AliITS.
    Int_t i,j,l;

    AliITSDetType          *det;
    AliITSresponse         *response;
    AliITSsegmentation     *seg;
    TClonesArray           *digitsaddress;
    TClonesArray           *clustaddr;


   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); 
      if (R__v == 1) {
	  AliDetector::Streamer(R__b);
	  R__b >> fITSgeom;
	  R__b >> fEuclidOut;
	  R__b >> fIdN;
	  if(fIdSens!=0) delete[] fIdSens;
	  if(fIdName!=0) delete[] fIdName;
	  fIdSens = new Int_t[fIdN];
	  fIdName = new char*[fIdN];
	  for(i=0;i<fIdN;i++) R__b >> fIdSens[i];
	  for(i=0;i<fIdN;i++){
	      R__b >> l;
	      fIdName[i] = new char[l+1]; // add room for null character.
	      for(j=0;j<l;j++) R__b >> fIdName[i][j];
	      fIdName[i][l] = '\0'; // Null terminate this string.
	  } // end for i
	  R__b >> fMajorVersion;
	  R__b >> fMinorVersion;
	  R__b >> fDtype;
	  R__b.ReadArray(fNdtype);
	  R__b >> fCtype;
	  R__b.ReadArray(fNctype);
	  R__b >> fDetTypes;
	  R__b >> fNDetTypes;
	  R__b >> fRecPoints;
	  R__b >> fNRecPoints;
          //  stream out only response and segmentation
          for(i=0;i<fNDetTypes;i++) {        
	       det=(AliITSDetType*)(*fDetTypes)[i];
	       det->Streamer(R__b);
	       response=det->GetResponseModel();
	       if (response) response->Streamer(R__b);	  
	       seg=det->GetSegmentationModel();
	       if (seg) seg->Streamer(R__b);	  
	       digitsaddress=(TClonesArray*) (*fDtype)[i];
	       digitsaddress->Streamer(R__b);
	       clustaddr=(TClonesArray*) (*fCtype)[i];
	       clustaddr->Streamer(R__b);
	  }

      } // end if (R__v)
   } else {
      R__b.WriteVersion(AliITS::IsA());
      AliDetector::Streamer(R__b);
      R__b << fITSgeom;
      R__b << fEuclidOut;
      R__b << fIdN;
      for(i=0;i<fIdN;i++) R__b <<fIdSens[i];
      for(i=0;i<fIdN;i++){
	  l = strlen(fIdName[i]);
	  R__b << l;
	  for(j=0;j<l;j++) R__b << fIdName[i][j];
      } // end for i
      R__b << fMajorVersion;
      R__b << fMinorVersion;
      R__b << fDtype;
      R__b.WriteArray(fNdtype, fNDetTypes);
      R__b << fCtype;
      R__b.WriteArray(fNctype, fNDetTypes);
      R__b << fDetTypes;
      R__b << fNDetTypes;
      R__b << fRecPoints;
      R__b << fNRecPoints;
      for(i=0;i<fNDetTypes;i++) {
	   det=(AliITSDetType*)(*fDetTypes)[i];
	   det->Streamer(R__b);
           response=det->GetResponseModel();
	   if(response) response->Streamer(R__b);
	   seg=det->GetSegmentationModel();
	   if (seg) seg->Streamer(R__b);	  
	   digitsaddress=(TClonesArray*) (*fDtype)[i];
	   digitsaddress->Streamer(R__b);
	   clustaddr=(TClonesArray*) (*fCtype)[i];
	   clustaddr->Streamer(R__b);
      }	  
   }


}

ClassImp(AliITSRecPoint)
ClassImp(AliITSsegmentation)
ClassImp(AliITSresponse)	
//ClassImp(AliITStrack)

