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
Revision 1.8  2000/10/20 06:24:40  fca
Put the PMD at the right position in the event display

Revision 1.7  2000/10/02 21:28:12  fca
Removal of useless dependecies via forward declarations

Revision 1.6  2000/01/19 17:17:06  fca
Introducing a list of lists of hits -- more hits allowed for detector now

Revision 1.5  1999/09/29 09:24:27  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//
//                                                                           //
//  Photon Multiplicity Detector                                             //
//  This class contains the basic functions for the Photon Multiplicity      //
//  Detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliPMDClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:sub@vecdec.veccal.ernet.in">Subhasis Chattopadhyay</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TBRIK.h>
#include <TNode.h>
#include <TTree.h>
#include <TGeometry.h>
#include <TClonesArray.h>

#include "AliPMD.h"
#include "AliRun.h"
#include "AliMC.h" 
#include "AliConst.h" 
#include "AliPMDRecPoint.h"
  
ClassImp(AliPMD)
 
//_____________________________________________________________________________
AliPMD::AliPMD()
{
  //
  // Default constructor
  //
  fIshunt = 0;
}
 
//_____________________________________________________________________________
AliPMD::AliPMD(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Default constructor
  //

  // 
  // Allocate the array of hits
  fHits   = new TClonesArray("AliPMDhit",  405);
  gAlice->AddHitList(fHits);

  fRecPoints  = new TClonesArray("AliPMDRecPoint",10000); 
  fNRecPoints = 0;
  

  fIshunt =  1;
  
  fPar[0] = 1;
  fPar[1] = 1;
  fPar[2] = 0.8;
  fPar[3] = 0.02;
  fIn[0]  = 6;
  fIn[1]  = 20;
  fIn[2]  = 600;
  fIn[3]  = 27;
  fIn[4]  = 27;
  fGeo[0] = 0;
  fGeo[1] = 0.2;
  fGeo[2] = 4;
  fPadSize[0] = 0.8;
  fPadSize[1] = 1.0;
  fPadSize[2] = 1.2;
  fPadSize[3] = 1.5;
}

AliPMD::~AliPMD()
{
  //
  // Default constructor
  //
    delete fRecPoints;
    fNRecPoints=0;
}

//_____________________________________________________________________________
void AliPMD::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a PMD hit
  //
  TClonesArray &lhits = *fHits;
  AliPMDhit *newcell, *curcell;
//    printf("PMD++ Adding energy %f, prim %d, vol %d %d %d %d\n",
//  	 hits[3],gAlice->GetPrimary(track-1),vol[0],vol[1],vol[2],vol[3]);
  newcell = new AliPMDhit(fIshunt, track, vol, hits);
  Int_t i;
  for (i=0; i<fNhits; i++) {
    //
    // See if this cell has already been hit
    curcell=(AliPMDhit*) lhits[i];
    if (*curcell==*newcell) {
//        printf("Cell with same numbers found\n") ; curcell->Print();
      *curcell = *curcell+*newcell;
//        printf("Cell after addition\n") ; curcell->Print();
      delete newcell;
      return;
    }
  }
  new(lhits[fNhits++]) AliPMDhit(newcell);
  delete newcell;
}
 
//_____________________________________________________________________________
void AliPMD::BuildGeometry()
{
  //
  // Build simple ROOT TNode geometry for event display
  //

  TNode *Node, *Top;
  const int kColorPMD  = kRed;

  //
  Top=gAlice->GetGeometry()->GetNode("alice");

  // PMD
  new TBRIK("S_PMD","PMD box","void",300,300,5);
  Top->cd();
  Node = new TNode("PMD","PMD","S_PMD",0,0,-600,"");
  Node->SetLineColor(kColorPMD);
  fNodes->Add(Node);
}

//_____________________________________________________________________________
Int_t AliPMD::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Distance from mouse to detector on the screen
  // dummy routine
  //
   return 9999;
}
 
//_____________________________________________________________________________
void AliPMD::SetPAR(Float_t p1, Float_t p2, Float_t p3,Float_t p4)
{
  //
  // Set PMD parameters
  //
  fPar[0] = p1;
  fPar[1] = p2;
  fPar[2] = p3;
  fPar[3] = p4;
}
 
//_____________________________________________________________________________
void AliPMD::SetIN(Float_t p1, Float_t p2, Float_t p3,Float_t p4,Float_t p5)
{
  //
  // Set PMD parameters
  //
  fIn[0] = p1;
  fIn[1] = p2;
  fIn[2] = p3;
  fIn[3] = p4;
  fIn[4] = p5;
}
 
//_____________________________________________________________________________
void AliPMD::SetGEO(Float_t p1, Float_t p2, Float_t p3)
{
  //
  // Set geometry parameters
  //
  fGeo[0] = p1;
  fGeo[1] = p2;
  fGeo[2] = p3;
}
 
//_____________________________________________________________________________
void AliPMD::SetPadSize(Float_t p1, Float_t p2, Float_t p3,Float_t p4)
{
  //
  // Set pad size
  //
  fPadSize[0] = p1;
  fPadSize[1] = p2;
  fPadSize[2] = p3;
  fPadSize[3] = p4;
}
 
//_____________________________________________________________________________
void AliPMD::StepManager()
{
  //
  // Called at every step in PMD
  //
}

void AliPMD::AddRecPoint(const AliPMDRecPoint &p)
{
    //
    // Add a PMD reconstructed hit to the list
    //
    TClonesArray &lrecpoints = *fRecPoints;
    new(lrecpoints[fNRecPoints++]) AliPMDRecPoint(p);
}

void AliPMD::MakeBranch(Option_t* option)
{
    // Create Tree branches for the PMD
    printf("Make Branch - TreeR address %p\n",gAlice->TreeR());
    
    const Int_t kBufferSize = 4000;
    char branchname[30];

    AliDetector::MakeBranch(option);

    sprintf(branchname,"%sRecPoints",GetName());
    if (fRecPoints   && gAlice->TreeR()) {
	gAlice->TreeR()->Branch(branchname, &fRecPoints, kBufferSize);
	printf("Making Branch %s for reconstructed hits\n",branchname);
    }	
}


void AliPMD::SetTreeAddress()
{
  // Set branch address for the TreeR
    char branchname[30];
    AliDetector::SetTreeAddress();

    TBranch *branch;
    TTree *treeR = gAlice->TreeR();

    sprintf(branchname,"%s",GetName());
    if (treeR && fRecPoints) {
	branch = treeR->GetBranch(branchname);
	if (branch) branch->SetAddress(&fRecPoints);
    }
}

void AliPMD::ResetHits()
{
  //
  // Reset number of hits and the hits array
  //
  fNRecPoints   = 0;
  if (fRecPoints)   fRecPoints->Clear();
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Photon Multiplicity Detector Version 1                                   //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliPMDv1Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



ClassImp(AliPMDhit)
  
//_____________________________________________________________________________
AliPMDhit::AliPMDhit(Int_t shunt,Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
  //
  // Add a PMD hit
  //
  Int_t i;
  for (i=0;i<5;i++) fVolume[i] = vol[i];
  fX=hits[0];
  fY=hits[1];
  fZ=hits[2];
  fEnergy=hits[3];
}
  
