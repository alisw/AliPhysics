///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Inner Traking System version 5                                           //
//  This class contains the base procedures for the Inner Tracking System    //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliITSv5Class.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:roberto.barbera@ct.infn.it">Roberto Barbera</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h> 

#include <TMath.h>
#include "AliITSv5.h"
#include "AliRun.h"
#include "stdlib.h"
#include "TSystem.h"
#include "TGeant3.h"

ClassImp(AliITSv5)
 
//_____________________________________________________________________________
AliITSv5::AliITSv5() : AliITS() 
{
  //
  // Default constructor for the ITS
  //
}
 
//_____________________________________________________________________________
AliITSv5::AliITSv5(const char *name, const char *title)
  : AliITS(name, title)
{
  //
  // Standard constructor for the ITS
  //
  fEuclidMaterial="$(ALICE_ROOT)/ITSgeometry_60.tme";
  fEuclidGeometry="$(ALICE_ROOT)/ITSgeometry_60.euc";
}

 
//_____________________________________________________________________________
void AliITSv5::CreateMaterials()
{
  //
  // Read materials for the ITS
  //
  FILE *file = fopen(fEuclidMaterial.Data(),"r");
  if(file) {
    fclose(file);
    gAlice->ReadEuclidMedia(fEuclidMaterial.Data(),2);
  } else {
    Error("CreateMaterials"," THE MEDIA FILE %s DOES NOT EXIST !",fEuclidMaterial.Data());
    exit(1);
  }
}

//_____________________________________________________________________________
void AliITSv5::CreateGeometry()
{
  //
  // Read geometry for the ITS
  //

  AliMC* pMC = AliMC::GetMC();
  
  char topvol[5];
  char *filtmp;
//
  filtmp=gSystem->ExpandPathName(fEuclidGeometry.Data());
  FILE *file = fopen(filtmp,"r");
  delete [] filtmp;
  if(file) {
    fclose(file);
    gAlice->ReadEuclid(fEuclidGeometry.Data(),2,topvol);
  } else {
    Error("CreateGeometry"," THE GEOM FILE %s DOES NOT EXIST !",fEuclidGeometry.Data());
    exit(1);
  }
  //
  // --- Place the ITS ghost volume ITSV in its mother volume (ALIC) and make it
  //     invisible
  //
  pMC->Gspos("ITSV",1,"ALIC",0,0,0,0,"ONLY");
  //
  // --- Outputs the geometry tree in the EUCLID/CAD format
  
    if (fEuclidOut) {
      pMC->WriteEuclid("ITSgeometry", "ITSV", 1, 5);
    }
}

//_____________________________________________________________________________
void AliITSv5::Init()
{
  //
  // Initialise the ITS after it has been created
  //
  AliITS::Init();
} 
 
//_____________________________________________________________________________
void AliITSv5::StepManager()
{
  //
  // Called for every step in the ITS
  //
  Int_t         copy, id;
  Int_t         copy1,copy2;
  Float_t       hits[7];
  Int_t         vol[3];
  Float_t       position[3];
  Float_t       momentum[4];
  TClonesArray &lhits = *fHits;
  AliMC* pMC = AliMC::GetMC();
  //
  if(pMC->TrackCharge() && pMC->Edep()) {
    //
    // Only entering charged tracks
    if((id=pMC->CurrentVol(0,copy))==fIdSens1) {  
      vol[0]=1;
      id=pMC->CurrentVolOff(0,0,copy);        //detector copy in the ladder = 1<->4  (ITS1)
      vol[1]=copy;
      pMC->CurrentVolOff(1,0,copy1);      //ladder copy in the module   = 1<->2  (I186)
      pMC->CurrentVolOff(2,0,copy2);      //module copy in the layer    = 1<->10 (I132)
      vol[2]=copy1+(copy2-1)*2;               //# of ladders in one module  = 2          
    } else if(id==fIdSens2) {
      vol[0]=2;
      id=pMC->CurrentVolOff(0,0,copy);        //detector copy in the ladder = 1<->4  (ITS2)        
      vol[1]=copy;
      pMC->CurrentVolOff(1,0,copy1);      //ladder copy in the module   = 1<->4  (I131)
      pMC->CurrentVolOff(2,0,copy2);      //module copy in the layer    = 1<->10 (I132)
      vol[2]=copy1+(copy2-1)*4;   	          //# of ladders in one module  = 4   
    } else if(id==fIdSens3) {
      vol[0]=3;
      id=pMC->CurrentVolOff(1,0,copy);        //detector copy in the ladder = 1<->5  (ITS3 is inside I314)
      vol[1]=copy;
      id=pMC->CurrentVolOff(2,0,copy);        //ladder copy in the layer    = 1<->12 (I316)
      vol[2]=copy;             
    } else if(id==fIdSens4) {
      vol[0]=4;
      id=pMC->CurrentVolOff(1,0,copy);	      //detector copy in the ladder = 1<->8  (ITS4 is inside I414)
      vol[1]=copy;
      id=pMC->CurrentVolOff(2,0,copy);        //ladder copy in the layer    = 1<->22 (I417)
      vol[2]=copy;               
    } else if(id==fIdSens5) {
      vol[0]=5;
      id=pMC->CurrentVolOff(1,0,copy);	      //detector copy in the ladder = 1<->23  (ITS5 is inside I562)	  
      vol[1]=copy;
      id=pMC->CurrentVolOff(2,0,copy);        //ladder copy in the layer    = 1<->34 (I565)
      vol[2]=copy;               
    } else if(id==fIdSens6) {
      vol[0]=6;
      id=pMC->CurrentVolOff(1,0,copy);	      //detector copy in the ladder = 1<->26  (ITS6 is inside I566)	  
      vol[1]=copy;
      id=pMC->CurrentVolOff(2,0,copy);        //ladder copy in the layer    = 1<->38 (I569)
      vol[2]=copy;                      
    } else return;
    pMC->TrackPosition(position);
    pMC->TrackMomentum(momentum);
    hits[0]=position[0];
    hits[1]=position[1];
    hits[2]=position[2];          
    hits[3]=momentum[0]*momentum[3];
    hits[4]=momentum[1]*momentum[3];
    hits[5]=momentum[2]*momentum[3];        
    hits[6]=pMC->Edep();
    new(lhits[fNhits++]) AliITShit(fIshunt,gAlice->CurrentTrack(),vol,hits);
  }      
}
