///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Forward Multiplicity Detector version 0                                  //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliFMDv0Class.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:Valeri.Kondratiev@cern.ch">Valeri Kondratiev</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliRun.h"
#include "AliFMDv0.h"
#include "AliMC.h"
#include "stdlib.h"
 
ClassImp(AliFMDv0)
 
//_____________________________________________________________________________
AliFMDv0::AliFMDv0()
{
  //
  // Default constructor for FMD version 0
  //
}
 
//_____________________________________________________________________________
AliFMDv0::AliFMDv0(const char *name, const char *title)
  : AliFMD(name,title)
{
  //
  // Standard constructor for FMD version 0
  //
  AliModule *start = gAlice->GetModule("START");
  if(start) {
    Error("ctor","This version of FMD is incompatible with START\n");
    exit(1);
  }
}
 
//___________________________________________
void AliFMDv0::CreateGeometry()
{
  //
  // Creation of the geometry of the FMD version 0
  //
  //Begin_Html
  /*
    <img src="picts/AliFMDv0Tree.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliFMDv0.gif">
  */
  //End_Html

  
  Float_t rout, z;
  Float_t par[3], rin;
  
  Int_t *idtmed = fIdtmed->GetArray()-899;    
  
  // ******************************************************** 
  //       DEFINE RIGHT DISK#3  OF FMD 
  // ******************************************************** 
  
  //       Define parameters 
  
  rin  = 4.5;
  rout = 10.5;
  z    = 85.;
  //       Ring #1 
  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1R3", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("R1R3", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");
  //       Ring #2 
  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2R3", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("R2R3", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");
  //       Ring #3 
  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3R3", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("R3R3", 1, "ALIC", 0., 0., z + 1.525, 0, "ONLY");
  // *********************************************************** 
  //       DEFINE LEFT DISK#3 OF FMD 
  // *********************************************************** 
  
  //       Define parameters 
  
  rin  = 4.5;
  rout = 10.5;
  z    = -85.;
  //       Ring #1 
  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1L3", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("R1L3", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //       Ring #2 
  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2L3", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("R2L3", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //       Ring #3 
  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3L3", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("R3L3", 1, "ALIC", 0., 0., z - 1.525, 0, "ONLY");
  // ******************************************************** 
  //       DEFINE RIGHT DISK#2  OF FMD 
  // ******************************************************** 
  
  //       Define parameters 
  
  rin  = 8.;
  rout = 14.;
  z    = 69.7;
  //       Ring #1 
  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1R2", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("R1R2", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");
  //       Ring #2 
  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2R2", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("R2R2", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");
  //       Ring #3 
  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3R2", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("R3R2", 1, "ALIC", 0., 0., z + 1.525, 0, "ONLY");
  // *********************************************************** 
  //       DEFINE LEFT DISK#2 OF FMD 
  // *********************************************************** 
  
  //       Define parameters 
  
  rin  = 8.;
  rout = 14.;
  z    = -69.7;
  //       Ring #1 
  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1L2", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("R1L2", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //       Ring #2 
  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2L2", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("R2L2", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //       Ring #3 
  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3L2", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("R3L2", 1, "ALIC", 0., 0., z - 1.525, 0, "ONLY");
  // ******************************************************** 
  //       DEFINE RIGHT DISK#1  OF FMD 
  // ******************************************************** 
  
  //       Define parameters 
  
  rin  = 8.;
  rout = 17.5;
  z    = 42.5;
  //       Ring #1 
  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1R1", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("R1R1", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");
  //       Ring #2 
  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2R1", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("R2R1", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");
  //       Ring #3 
  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3R1", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("R3R1", 1, "ALIC", 0., 0., z + 1.525, 0, "ONLY");
  // *********************************************************** 
  //       DEFINE LEFT DISK#1 OF FMD 
  // *********************************************************** 
  
  //       Define parameters 
  
  rin  = 8.;
  rout = 17.5;
  z    = -42.5;
  //       Ring #1 
  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1L1", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("R1L1", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //       Ring #2 
  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2L1", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("R2L1", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //       Ring #3 
  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3L1", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("R3L1", 1, "ALIC", 0., 0., z - 1.525, 0, "ONLY");
  // *********************************************************** 
  //       DEFINE LEFT DISK#4 OF FMD 
  // *********************************************************** 

  //       Define parameters 

  rin  = 4.2;
  rout = 13.;
  z    = -229.5;
  //       Ring #1 
  par[0] = rin;
  par[1] = rout;
  par[2] = 1.5;
  gMC->Gsvolu("R1L4", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("R1L4", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //       Ring #2 
  par[0] = rout;
  par[1] = rout + .65;
  par[2] = 1.5;
  gMC->Gsvolu("R2L4", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("R2L4", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //       Ring #3 
  par[0] = rout + .65;
  par[1] = rout + 5.65;
  par[2] = .025;
  gMC->Gsvolu("R3L4", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("R3L4", 1, "ALIC", 0., 0., z - 1.525, 0, "ONLY");
}

//_____________________________________________________________________________
void AliFMDv0::DrawModule()
{
  //
  // Draw a shaded view of the FMD version 0
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("R1R3","SEEN",1);
  gMC->Gsatt("R2R3","SEEN",1);
  gMC->Gsatt("R3R3","SEEN",1);
  gMC->Gsatt("R1L3","SEEN",1);
  gMC->Gsatt("R2L3","SEEN",1);
  gMC->Gsatt("R3L3","SEEN",1);
  gMC->Gsatt("R1R2","SEEN",1);
  gMC->Gsatt("R2R2","SEEN",1);
  gMC->Gsatt("R3R2","SEEN",1);
  gMC->Gsatt("R1L2","SEEN",1);
  gMC->Gsatt("R2L2","SEEN",1);
  gMC->Gsatt("R3L2","SEEN",1);
  gMC->Gsatt("R1R1","SEEN",1);
  gMC->Gsatt("R2R1","SEEN",1);
  gMC->Gsatt("R3R1","SEEN",1);
  gMC->Gsatt("R1L1","SEEN",1);
  gMC->Gsatt("R2L1","SEEN",1);
  gMC->Gsatt("R3L1","SEEN",1);
  gMC->Gsatt("R1L4","SEEN",1);
  gMC->Gsatt("R2L4","SEEN",1);
  gMC->Gsatt("R3L4","SEEN",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 6, 9, .08, .08);
  gMC->Gdhead(1111, "Forward Multiplicity Detector version 0");
  gMC->Gdman(13, 9, "MAN");
}

//_____________________________________________________________________________
void AliFMDv0::CreateMaterials()
{
  //
  // Create Materials for version 0 of FMD
  //

  //     Material for ring #1 
  Float_t ar1[8] = { 55.8,58.7,52.,47.9,16.,28.,207.2,27. };
  Float_t zr1[8] = { 26.,28.,24.,22.,8.,14.,82.,13. };
  Float_t wr1[8] = { .27,.081,.054,.045,.18,.25,.06,.06 };
  //     Material for ring #2 
  Float_t ar2[3] = { 55.8,27.,16. };
  Float_t zr2[3] = { 26.,13.,8. };
  Float_t wr2[3] = { .35,.34,.31 };
  //     Material for ring #3 
  Float_t ar3[3] = { 28.,27.,16. };
  Float_t zr3[3] = { 14.,13.,8. };
  Float_t wr3[3] = { .37,.33,.3 };
  
  Float_t epsil, stmin, deemax, tmaxfd, stemax;
  
  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  //     Ring #1 
  
  AliMixture(1, "FMD_R1$", ar1, zr1, 2.69, 8, wr1);
  
  //     Ring #2 
  
  AliMixture(2, "FMD_R2$", ar2, zr2, 2.63, 3, wr2);
  
  //     Ring #3 
  
  AliMixture(3, "FMD_R3$", ar3, zr3, 3.15, 3, wr3);
  // ******************************************************* 
  //     Defines tracking media parameters. 
  // ******************************************************* 
  epsil  = .001; // Tracking precision, DLS 
  stemax = -1.;  // Maximum displacement for multiple scattering 
  tmaxfd = -20.; // Maximum angle due to field deflection 
  deemax = -.3;  // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  // ******************************************************** 
  AliMedium(1, "FMD_R1_L3        ", 1, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(2, "FMD_R2_L3        ", 2, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(3, "FMD_R3_L3        ", 3, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
}


 
