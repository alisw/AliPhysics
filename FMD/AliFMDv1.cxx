///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Forward Multiplicity Detector Version 1                                  //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliFMDv1Class.gif">
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
#include "AliFMDv1.h"
#include "AliMC.h"
#include "AliConst.h"
 
 
ClassImp(AliFMDv1)
 
//_____________________________________________________________________________
AliFMDv1::AliFMDv1()
{
  //
  // Defautl constructor for FMD version 1
  //
}
 
//_____________________________________________________________________________
AliFMDv1::AliFMDv1(const char *name, const char *title)
  : AliFMD(name,title)
{
  //
  // Standard constructor for FMD version 1
  //
  AliModule *start = gAlice->GetModule("START");
  if(start) {
    Error("ctor","This version of FMD is incompatible with START\n");
    exit(1);
  }
}
 
//_____________________________________________________________________________
void AliFMDv1::CreateGeometry()
{
  //
  // Creation of the geometry of the FMD version 1
  //
  //Begin_Html
  /*
    <img src="picts/AliFMDv1Tree.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliFMDv1.gif">
  */
  //End_Html

  
  Float_t rout;
  Float_t z;
  Float_t h1, h2, t0, t1, t2;
  Float_t tt;
  Float_t par[3], rin;
  const Float_t v1 = 42.0;
  const Float_t v2 = 4.5;
  const Float_t v3 = 5.62;
  const Float_t v4 = 16.0;
  
  Int_t *idtmed = fIdtmed->GetArray()-899;    
  
  // ******************************************************** 
  //     DEFINE RIGHT DISK#3  OF FMD 
  // ******************************************************** 
  
  //     Define parameters 
  
  rin  = 4.5;
  rout = 10.5;
  z    = 85.;
  tt   = 2.5;
  h1   = TMath::Sqrt(rout * rout + v1/ (tt * kPI)) - rout;
  h2   = .4;
  t0   = v2/ (h1 * kPI * (h1 + rout * 2.));
  t1   = v3/ (h2 * kPI * (h2 + (h1 + rout) * 2.));
  t2   = v4/ (h2 * kPI * (h2 + (h1 + rout) * 2.));
  
  // *******Inner slice*********** 
  //     Inner steel wall 
  par[0] = rin - .02;
  par[1] = rin;
  par[2] = 1.5;
  gMC->Gsvolu("IWR3", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("IWR3", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");
  //     Front steel wall 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("FWR3", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("FWR3", 1, "ALIC", 0., 0., z + .01, 0, "ONLY");
  //     Rear steel wall 
  gMC->Gsvolu("RWR3", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("RWR3", 1, "ALIC", 0., 0., z + 2.99, 0, "ONLY");
  //     MCP 
  par[0] = rin;
  par[1] = rout;
  par[2] = .07;
  gMC->Gsvolu("MPR3", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("MPR3", 1, "ALIC", 0., 0., z + 1.57, 0, "ONLY");
  //     Silicon plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .019;
  gMC->Gsvolu("SPR3", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SPR3", 1, "ALIC", 0., 0., z + 1.719, 0, "ONLY");
  //     Summator plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("SMR3", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SMR3", 1, "ALIC", 0., 0., z + 2.01, 0, "ONLY");
  // *******Outer slice ******* 
  //     Ceramic plate 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = 1.25;
  gMC->Gsvolu("CPR3", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CPR3", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");
  //     Covar spring 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = t0 / 2.;
  gMC->Gsvolu("C1R3", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C1R3", 1, "ALIC", 0., 0., z + .25 - t0 / 2., 0, "ONLY");
  gMC->Gsvolu("C2R3", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C2R3", 1, "ALIC", 0., 0., z + 2.75 + t0 / 2., 0, "ONLY");
  //     Getter camera 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t1 / 2.;
  gMC->Gsvolu("GKR3", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("GKR3", 1, "ALIC", 0., 0., z + .25 + t1 / 2., 0, "ONLY");
  //     Ceramic slice 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t2 / 2.;
  gMC->Gsvolu("SCR3", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SCR3", 1, "ALIC", 0., 0., z + .25 + t1 + t2 / 2., 0, "ONLY");
  // ******Electronic slice ******* 
  //     Silicon ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("SER3", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SER3", 1, "ALIC", 0., 0., z + 1.57 - .025, 0, "ONLY");
  //     Ceramic ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("CER3", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CER3", 1, "ALIC", 0., 0., z + 1.58 + .025, 0, "ONLY");
  // *********************************************************** 
  //     DEFINE LEFT DISK#3 OF FMD 
  // *********************************************************** 
  
  //     Define parameters 
  
  rin  = 4.5;
  rout = 10.5;
  z    = -85.;
  tt   = 2.5;
  h1   = TMath::Sqrt(rout * rout + v1/ (tt * kPI)) - rout;
  h2   = .4;
  t0   = v2/ (h1 * kPI * (h1 + rout * 2.));
  t1   = v3/ (h2 * kPI * (h2 + (h1 + rout) * 2.));
  t2   = v4/  (h2 * kPI * (h2 + (h1 + rout) * 2.));
  
  // *******Inner slice*********** 
  //     Inner steel wall 
  par[0] = rin - .02;
  par[1] = rin;
  par[2] = 1.5;
  gMC->Gsvolu("IWL3", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("IWL3", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //     Front steel wall 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("FWL3", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("FWL3", 1, "ALIC", 0., 0., z - .01, 0, "ONLY");
  //     Rear steel wall 
  gMC->Gsvolu("RWL3", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("RWL3", 1, "ALIC", 0., 0., z - 2.99, 0, "ONLY");
  //     MCP 
  par[0] = rin;
  par[1] = rout;
  par[2] = .07;
  gMC->Gsvolu("MPL3", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("MPL3", 1, "ALIC", 0., 0., z - 1.57, 0, "ONLY");
  //     Silicon plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .019;
  gMC->Gsvolu("SPL3", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SPL3", 1, "ALIC", 0., 0., z - 1.719, 0, "ONLY");
  //     Summator plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("SML3", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SML3", 1, "ALIC", 0., 0., z - 2.01, 0, "ONLY");
  // *******Outer slice ******* 
  //     Ceramic plate 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = 1.25;
  gMC->Gsvolu("CPL3", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CPL3", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //     Covar spring 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = t0 / 2.;
  gMC->Gsvolu("C1L3", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C1L3", 1, "ALIC", 0., 0., z - .25 + t0 / 2., 0, "ONLY");
  gMC->Gsvolu("C2L3", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C2L3", 1, "ALIC", 0., 0., z - 2.75 - t0 / 2., 0, "ONLY");
  //     Getter camera 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t1 / 2.;
  gMC->Gsvolu("GKL3", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("GKL3", 1, "ALIC", 0., 0., z - .25 - t1 / 2., 0, "ONLY");
  //     Ceramic slice 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t2 / 2.;
  gMC->Gsvolu("SCL3", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SCL3", 1, "ALIC", 0., 0., z - .25 - t1 - t2 / 2., 0, "ONLY");
  // ******Electronic slice ******* 
  //     Silicon ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("SEL3", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SEL3", 1, "ALIC", 0., 0., z - 1.57 + .025, 0, "ONLY");
  //     Ceramic ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("CEL3", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CEL3", 1, "ALIC", 0., 0., z - 1.58 - .025, 0, "ONLY");
  // ******************************************************** 
  //     DEFINE RIGHT DISK#2  OF FMD 
  // ******************************************************** 
  
  //     Define parameters 
  
  rin  = 8.;
  rout = 14.;
  z    = 69.7;
  tt   = 2.5;
  h1   = TMath::Sqrt(rout * rout + v1/ (tt * kPI)) - rout;
  h2   = .4;
  t0   = v2/ (h1 * kPI * (h1 + rout * 2.));
  t1   = v3/ (h2 * kPI * (h2 + (h1 + rout) * 2.));
  t2   = v4/  (h2 * kPI * (h2 + (h1 + rout) * 2.));
  
  // *******Inner slice*********** 
  //     Inner steel wall 
  par[0] = rin - .02;
  par[1] = rin;
  par[2] = 1.5;
  gMC->Gsvolu("IWR2", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("IWR2", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");
  //     Front steel wall 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("FWR2", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("FWR2", 1, "ALIC", 0., 0., z + .01, 0, "ONLY");
  //     Rear steel wall 
  gMC->Gsvolu("RWR2", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("RWR2", 1, "ALIC", 0., 0., z + 2.99, 0, "ONLY");
  //     MCP 
  par[0] = rin;
  par[1] = rout;
  par[2] = .07;
  gMC->Gsvolu("MPR2", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("MPR2", 1, "ALIC", 0., 0., z + 1.57, 0, "ONLY");
  //     Silicon plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .019;
  gMC->Gsvolu("SPR2", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SPR2", 1, "ALIC", 0., 0., z + 1.719, 0, "ONLY");
  //     Summator plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("SMR2", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SMR2", 1, "ALIC", 0., 0., z + 2.01, 0, "ONLY");
  // *******Outer slice ******* 
  //     Ceramic plate 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = 1.25;
  gMC->Gsvolu("CPR2", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CPR2", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");
  //     Covar spring 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = t0 / 2.;
  gMC->Gsvolu("C1R2", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C1R2", 1, "ALIC", 0., 0., z + .25 - t0 / 2., 0, "ONLY");
  gMC->Gsvolu("C2R2", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C2R2", 1, "ALIC", 0., 0., z + 2.75 + t0 / 2., 0, "ONLY");
  //     Getter camera 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t1 / 2.;
  gMC->Gsvolu("GKR2", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("GKR2", 1, "ALIC", 0., 0., z + .25 + t1 / 2., 0, "ONLY");
  //     Ceramic slice 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t2 / 2.;
  gMC->Gsvolu("SCR2", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SCR2", 1, "ALIC", 0., 0., z + .25 + t1 + t2 / 2., 0, "ONLY");
  // ******Electronic slice ******* 
  //     Silicon ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("SER2", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SER2", 1, "ALIC", 0., 0., z + 1.57 - .025, 0, "ONLY");
  //     Ceramic ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("CER2", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CER2", 1, "ALIC", 0., 0., z + 1.58 + .025, 0, "ONLY");
  // *********************************************************** 
  //       DEFINE LEFT DISK#2 OF FMD 
  // *********************************************************** 
  
  //       Define parameters 
  
  rin  = 8.;
  rout = 14.;
  z    = -69.7;
  tt   = 2.5;
  h1   = TMath::Sqrt(rout * rout + v1/ (tt * kPI)) - rout;
  h2   = .4;
  t0   = v2/ (h1 * kPI * (h1 + rout * 2.));
  t1   = v3/ (h2 * kPI * (h2 + (h1 + rout) * 2.));
  t2   = v4/  (h2 * kPI * (h2 + (h1 + rout) * 2.));
  
  // *******Inner slice*********** 
  //     Inner steel wall 
  par[0] = rin - .02;
  par[1] = rin;
  par[2] = 1.5;
  gMC->Gsvolu("IWL2", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("IWL2", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //     Front steel wall 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("FWL2", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("FWL2", 1, "ALIC", 0., 0., z - .01, 0, "ONLY");
  //     Rear steel wall 
  gMC->Gsvolu("RWL2", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("RWL2", 1, "ALIC", 0., 0., z - 2.99, 0, "ONLY");
  //     MCP 
  par[0] = rin;
  par[1] = rout;
  par[2] = .07;
  gMC->Gsvolu("MPL2", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("MPL2", 1, "ALIC", 0., 0., z - 1.57, 0, "ONLY");
  //     Silicon plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .019;
  gMC->Gsvolu("SPL2", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SPL2", 1, "ALIC", 0., 0., z - 1.719, 0, "ONLY");
  //     Summator plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("SML2", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SML2", 1, "ALIC", 0., 0., z - 2.01, 0, "ONLY");
  // *******Outer slice ******* 
  //     Ceramic plate 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = 1.25;
  gMC->Gsvolu("CPL2", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CPL2", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //     Covar spring 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = t0 / 2.;
  gMC->Gsvolu("C1L2", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C1L2", 1, "ALIC", 0., 0., z - .25 + t0 / 2., 0, "ONLY");
  gMC->Gsvolu("C2L2", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C2L2", 1, "ALIC", 0., 0., z - 2.75 - t0 / 2., 0, "ONLY");
  //     Getter camera 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t1 / 2.;
  gMC->Gsvolu("GKL2", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("GKL2", 1, "ALIC", 0., 0., z - .25 - t1 / 2., 0, "ONLY");
  //     Ceramic slice 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t2 / 2.;
  gMC->Gsvolu("SCL2", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SCL2", 1, "ALIC", 0., 0., z - .25 - t1 - t2 / 2., 0, "ONLY");
  // ******Electronic slice ******* 
  //     Silicon ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("SEL2", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SEL2", 1, "ALIC", 0., 0., z - 1.57 + .025, 0, "ONLY");
  //     Ceramic ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("CEL2", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CEL2", 1, "ALIC", 0., 0., z - 1.58 - .025, 0, "ONLY");
  // ******************************************************** 
  //     DEFINE RIGHT DISK#1  OF FMD 
  // ******************************************************** 
  
  //     Define parameters 

  rin  = 8.;
  rout = 17.5;
  z    = 42.5;
  tt   = 2.5;
  h1   = TMath::Sqrt(rout * rout + v1/ (tt * kPI)) - rout;
  h2   = .4;
  t0   = v2/ (h1 * kPI * (h1 + rout * 2.));
  t1   = v3/ (h2 * kPI * (h2 + (h1 + rout) * 2.));
  t2   = v4/ (h2 * kPI * (h2 + (h1 + rout) * 2.));
  
  // *******Inner slice*********** 
  //     Inner steel wall 
  par[0] = rin - .02;
  par[1] = rin;
  par[2] = 1.5;
  gMC->Gsvolu("IWR1", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("IWR1", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");
  //     Front steel wall 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("FWR1", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("FWR1", 1, "ALIC", 0., 0., z + .01, 0, "ONLY");
  //     Rear steel wall 
  gMC->Gsvolu("RWR1", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("RWR1", 1, "ALIC", 0., 0., z + 2.99, 0, "ONLY");
  //     MCP 
  par[0] = rin;
  par[1] = rout;
  par[2] = .07;
  gMC->Gsvolu("MPR1", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("MPR1", 1, "ALIC", 0., 0., z + 1.57, 0, "ONLY");
  //     Silicon plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .019;
  gMC->Gsvolu("SPR1", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SPR1", 1, "ALIC", 0., 0., z + 1.719, 0, "ONLY");
  //     Summator plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("SMR1", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SMR1", 1, "ALIC", 0., 0., z + 2.01, 0, "ONLY");
  // *******Outer slice ******* 
  //     Ceramic plate 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = 1.25;
  gMC->Gsvolu("CPR1", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CPR1", 1, "ALIC", 0., 0., z + 1.5, 0, "ONLY");
  //     Covar spring 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = t0 / 2.;
  gMC->Gsvolu("C1R1", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C1R1", 1, "ALIC", 0., 0., z + .25 - t0 / 2., 0, "ONLY");
  gMC->Gsvolu("C2R1", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C2R1", 1, "ALIC", 0., 0., z + 2.75 + t0 / 2., 0, "ONLY");
  //     Getter camera 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t1 / 2.;
  gMC->Gsvolu("GKR1", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("GKR1", 1, "ALIC", 0., 0., z + .25 + t1 / 2., 0, "ONLY");
  //     Ceramic slice 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t2 / 2.;
  gMC->Gsvolu("SCR1", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SCR1", 1, "ALIC", 0., 0., z + .25 + t1 + t2 / 2., 0, "ONLY");
  // ******Electronic slice ******* 
  //     Silicon ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("SER1", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SER1", 1, "ALIC", 0., 0., z + 1.57 - .025, 0, "ONLY");
  //     Ceramic ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("CER1", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CER1", 1, "ALIC", 0., 0., z + 1.58 + .025, 0, "ONLY");
  // *********************************************************** 
  //     DEFINE LEFT DISK#1 OF FMD 
  // *********************************************************** 
  
  //     Define parameters 
  
  rin  = 8.;
  rout = 17.5;
  z    = -42.5;
  tt   = 2.5;
  h1   = TMath::Sqrt(rout * rout + v1/ (tt * kPI)) - rout;
  h2   = .4;
  t0   = v2/ (h1 * kPI * (h1 + rout * 2.));
  t1   = v3/ (h2 * kPI * (h2 + (h1 + rout) * 2.));
  t2   = v4/ (h2 * kPI * (h2 + (h1 + rout) * 2.));
  
  // *******Inner slice*********** 
  //     Inner steel wall 
  par[0] = rin - .02;
  par[1] = rin;
  par[2] = 1.5;
  gMC->Gsvolu("IWL1", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("IWL1", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //     Front steel wall 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("FWL1", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("FWL1", 1, "ALIC", 0., 0., z - .01, 0, "ONLY");
  //     Rear steel wall 
  gMC->Gsvolu("RWL1", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("RWL1", 1, "ALIC", 0., 0., z - 2.99, 0, "ONLY");
  //     MCP 
  par[0] = rin;
  par[1] = rout;
  par[2] = .07;
  gMC->Gsvolu("MPL1", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("MPL1", 1, "ALIC", 0., 0., z - 1.57, 0, "ONLY");
  //     Silicon plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .019;
  gMC->Gsvolu("SPL1", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SPL1", 1, "ALIC", 0., 0., z - 1.719, 0, "ONLY");
  //     Summator plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("SML1", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SML1", 1, "ALIC", 0., 0., z - 2.01, 0, "ONLY");
  // *******Outer slice ******* 
  //     Ceramic plate 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = 1.25;
  gMC->Gsvolu("CPL1", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CPL1", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //     Covar spring 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = t0 / 2.;
  gMC->Gsvolu("C1L1", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C1L1", 1, "ALIC", 0., 0., z - .25 + t0 / 2., 0, "ONLY");
  gMC->Gsvolu("C2L1", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C2L1", 1, "ALIC", 0., 0., z - 2.75 - t0 / 2., 0, "ONLY");
  //     Getter camera 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t1 / 2.;
  gMC->Gsvolu("GKL1", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("GKL1", 1, "ALIC", 0., 0., z - .25 - t1 / 2., 0, "ONLY");
  //     Ceramic slice 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t2 / 2.;
  gMC->Gsvolu("SCL1", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SCL1", 1, "ALIC", 0., 0., z - .25 - t1 - t2 / 2., 0, "ONLY");
  // ******Electronic slice ******* 
  //     Silicon ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("SEL1", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SEL1", 1, "ALIC", 0., 0., z - 1.57 + .025, 0, "ONLY");
  //     Ceramic ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("CEL1", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CEL1", 1, "ALIC", 0., 0., z - 1.58 - .025, 0, "ONLY");
  // *********************************************************** 
  //     DEFINE LEFT DISK#4 OF FMD 
  // *********************************************************** 
  
  //     Define parameters 
  
  rin  = 4.2;
  rout = 13.;
  z    = -229.5;
  tt   = 2.5;
  h1   = TMath::Sqrt(rout * rout + v1/ (tt * kPI)) - rout;
  h2   = .4;
  t0   = v2/ (h1 * kPI * (h1 + rout * 2.));
  t1   = v3/ (h2 * kPI * (h2 + (h1 + rout) * 2.));
  t2   = v4/ (h2 * kPI * (h2 + (h1 + rout) * 2.));
  
  // *******Inner slice*********** 
  //     Inner steel wall 
  par[0] = rin - .02;
  par[1] = rin;
  par[2] = 1.5;
  gMC->Gsvolu("IWL4", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("IWL4", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //     Front steel wall 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("FWL4", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("FWL4", 1, "ALIC", 0., 0., z - .01, 0, "ONLY");
  //     Rear steel wall 
  gMC->Gsvolu("RWL4", "TUBE", idtmed[899], par, 3);
  gMC->Gspos("RWL4", 1, "ALIC", 0., 0., z - 2.99, 0, "ONLY");
  //     MCP 
  par[0] = rin;
  par[1] = rout;
  par[2] = .07;
  gMC->Gsvolu("MPL4", "TUBE", idtmed[900], par, 3);
  gMC->Gspos("MPL4", 1, "ALIC", 0., 0., z - 1.57, 0, "ONLY");
  //     Silicon plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .019;
  gMC->Gsvolu("SPL4", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SPL4", 1, "ALIC", 0., 0., z - 1.719, 0, "ONLY");
  //     Summator plate 
  par[0] = rin;
  par[1] = rout;
  par[2] = .01;
  gMC->Gsvolu("SML4", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SML4", 1, "ALIC", 0., 0., z - 2.01, 0, "ONLY");
  // *******Outer slice ******* 
  //     Ceramic plate 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = 1.25;
  gMC->Gsvolu("CPL4", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CPL4", 1, "ALIC", 0., 0., z - 1.5, 0, "ONLY");
  //     Covar spring 
  par[0] = rout;
  par[1] = rout + h1;
  par[2] = t0 / 2.;
  gMC->Gsvolu("C1L4", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C1L4", 1, "ALIC", 0., 0., z - .25 + t0 / 2., 0, "ONLY");
  gMC->Gsvolu("C2L4", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("C2L4", 1, "ALIC", 0., 0., z - 2.75 - t0 / 2., 0, "ONLY");
  //       Getter camera 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t1 / 2.;
  gMC->Gsvolu("GKL4", "TUBE", idtmed[903], par, 3);
  gMC->Gspos("GKL4", 1, "ALIC", 0., 0., z - .25 - t1 / 2., 0, "ONLY");
  //     Ceramic slice 
  par[0] = rout + h1;
  par[1] = rout + h1 + h2;
  par[2] = t2 / 2.;
  gMC->Gsvolu("SCL4", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("SCL4", 1, "ALIC", 0., 0., z - .25 - t1 - t2 / 2., 0, "ONLY");
  // ******Electronic slice ******* 
  //     Silicon ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("SEL4", "TUBE", idtmed[901], par, 3);
  gMC->Gspos("SEL4", 1, "ALIC", 0., 0., z - 1.57 + .025, 0, "ONLY");
  //     Ceramic ring 
  par[0] = rout + h1 + h2;
  par[1] = rout + h1 + h2 + 5.;
  par[2] = .025;
  gMC->Gsvolu("CEL4", "TUBE", idtmed[902], par, 3);
  gMC->Gspos("CEL4", 1, "ALIC", 0., 0., z - 1.58 - .025, 0, "ONLY");
}
 
//_____________________________________________________________________________
void AliFMDv1::DrawModule()
{
  //
  // Draw a shaded view of the FMD version 1
  //

  
  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("IWR3","seen",1);
  gMC->Gsatt("FWR3","seen",1);
  gMC->Gsatt("RWR3","seen",1);
  gMC->Gsatt("MPR3","seen",1);
  gMC->Gsatt("SPR3","seen",1);
  gMC->Gsatt("SMR3","seen",1);
  gMC->Gsatt("CPR3","seen",1);
  gMC->Gsatt("C1R3","seen",1);
  gMC->Gsatt("C2R3","seen",1);
  gMC->Gsatt("GKR3","seen",1);
  gMC->Gsatt("SCR3","seen",1);
  gMC->Gsatt("SER3","seen",1);
  gMC->Gsatt("CER3","seen",1);
  gMC->Gsatt("IWL3","seen",1);
  gMC->Gsatt("FWL3","seen",1);
  gMC->Gsatt("RWL3","seen",1);
  gMC->Gsatt("MPL3","seen",1);
  gMC->Gsatt("SPL3","seen",1);
  gMC->Gsatt("SML3","seen",1);
  gMC->Gsatt("CPL3","seen",1);
  gMC->Gsatt("C1L3","seen",1);
  gMC->Gsatt("C2L3","seen",1);
  gMC->Gsatt("GKL3","seen",1);
  gMC->Gsatt("SCL3","seen",1);
  gMC->Gsatt("SEL3","seen",1);
  gMC->Gsatt("CEL3","seen",1);
  gMC->Gsatt("IWR2","seen",1);
  gMC->Gsatt("FWR2","seen",1);
  gMC->Gsatt("RWR2","seen",1);
  gMC->Gsatt("MPR2","seen",1);
  gMC->Gsatt("SPR2","seen",1);
  gMC->Gsatt("SMR2","seen",1);
  gMC->Gsatt("CPR2","seen",1);
  gMC->Gsatt("C1R2","seen",1);
  gMC->Gsatt("C2R2","seen",1);
  gMC->Gsatt("GKR2","seen",1);
  gMC->Gsatt("SCR2","seen",1);
  gMC->Gsatt("SER2","seen",1);
  gMC->Gsatt("CER2","seen",1);
  gMC->Gsatt("IWL2","seen",1);
  gMC->Gsatt("FWL2","seen",1);
  gMC->Gsatt("RWL2","seen",1);
  gMC->Gsatt("MPL2","seen",1);
  gMC->Gsatt("SPL2","seen",1);
  gMC->Gsatt("SML2","seen",1);
  gMC->Gsatt("CPL2","seen",1);
  gMC->Gsatt("C1L2","seen",1);
  gMC->Gsatt("C2L2","seen",1);
  gMC->Gsatt("GKL2","seen",1);
  gMC->Gsatt("SCL2","seen",1);
  gMC->Gsatt("SEL2","seen",1);
  gMC->Gsatt("CEL2","seen",1);
  gMC->Gsatt("IWR1","seen",1);
  gMC->Gsatt("FWR1","seen",1);
  gMC->Gsatt("RWR1","seen",1);
  gMC->Gsatt("MPR1","seen",1);
  gMC->Gsatt("SPR1","seen",1);
  gMC->Gsatt("SMR1","seen",1);
  gMC->Gsatt("CPR1","seen",1);
  gMC->Gsatt("C1R1","seen",1);
  gMC->Gsatt("C2R1","seen",1);
  gMC->Gsatt("GKR1","seen",1);
  gMC->Gsatt("SCR1","seen",1);
  gMC->Gsatt("SER1","seen",1);
  gMC->Gsatt("CER1","seen",1);
  gMC->Gsatt("IWL1","seen",1);
  gMC->Gsatt("FWL1","seen",1);
  gMC->Gsatt("RWL1","seen",1);
  gMC->Gsatt("MPL1","seen",1);
  gMC->Gsatt("SPL1","seen",1);
  gMC->Gsatt("SML1","seen",1);
  gMC->Gsatt("CPL1","seen",1);
  gMC->Gsatt("C1L1","seen",1);
  gMC->Gsatt("C2L1","seen",1);
  gMC->Gsatt("GKL1","seen",1);
  gMC->Gsatt("SCL1","seen",1);
  gMC->Gsatt("SEL1","seen",1);
  gMC->Gsatt("CEL1","seen",1);
  gMC->Gsatt("IWL4","seen",1);
  gMC->Gsatt("FWL4","seen",1);
  gMC->Gsatt("RWL4","seen",1);
  gMC->Gsatt("MPL4","seen",1);
  gMC->Gsatt("SPL4","seen",1);
  gMC->Gsatt("SML4","seen",1);
  gMC->Gsatt("CPL4","seen",1);
  gMC->Gsatt("C1L4","seen",1);
  gMC->Gsatt("C2L4","seen",1);
  gMC->Gsatt("GKL4","seen",1);
  gMC->Gsatt("SCL4","seen",1);
  gMC->Gsatt("SEL4","seen",1);
  gMC->Gsatt("CEL4","seen",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 6, 9, .08, .08);
  gMC->Gdhead(1111, "Forward Multiplicity Detector version 1");
  gMC->Gdman(13, 9, "MAN");
}

//_____________________________________________________________________________
void AliFMDv1::CreateMaterials()
{
  //
  // Create materials for version 1 of FMD
  //

  //     Steel for chamber 
  Float_t ast[4]  = { 55.847,58.71,51.996,47.9 };
  Float_t zst[4]  = { 26.,28.,24.,22. };
  Float_t wst[4]  = { .6,.18,.12,.1 };
  //     Lead glass for MCP 
  Float_t amcp[3] = { 15.9994,28.086,207.19 };
  Float_t zmcp[3] = { 8.,14.,82. };
  Float_t wmcp[3] = { .43,.375,.195 };
  //     Ceramic for plates 
  Float_t acer[2] = { 26.98,15.9994 };
  Float_t zcer[2] = { 13.,8. };
  Float_t wcer[2] = { .4,.6 };
  
  Float_t epsil, stmin, deemax, tmaxfd, stemax;
  
  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  //     Silicon 
  
  AliMaterial(2, "SI$", 28.086, 14., 2.33, 9.36, 45.5);
  
  //     Covar 
  
  AliMaterial(4, "FE$", 55.85, 26., 8.3, 1.67, 15.9);
  
  //     Steel 
  
  AliMixture(0, "FMD_STEEL$", ast, zst, 7.9, 4, wst);
  
  //     Lead glass 
  
  AliMixture(1, "FMD_MCP$", amcp, zmcp, 1.55, 3, wmcp);
  
  //     Ceramic 
  
  AliMixture(3, "FMD_CERAMIC$", acer, zcer, 3.96, -2, wcer);
  // ******************************************************* 
  //     Defines tracking media parameters. 
  // ******************************************************* 
  epsil  = .001; // Tracking precision, DLS 
  stemax = -1.;  // Maximum displacement for multiple scattering 
  tmaxfd = -20.; // Maximum angle due to field deflection 
  deemax = -.3;  // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  // ******************************************************** 
  AliMedium(0, "STEEL_L3          ", 0, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(1, "LEAD_GLASS_L3     ", 1, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(2, "SILICON_L3        ", 2, 1, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(3, "CERAMIC_L3        ", 3, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(4, "COVAR_L3          ", 4, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
}

