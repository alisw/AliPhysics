///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Zero Degree Calorimeter                                                  //
//  This class contains the basic functions for the Time Of Flight           //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliZDCClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:Eugenio.Scomparin@cern.ch">Eugenio Scomparin</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TBRIK.h>
#include <TNode.h>

#include "AliZDC.h"
#include "AliRun.h"
#include "AliCallf77.h"
#include "AliConst.h"
#include "AliMC.h"

#ifndef WIN32
# define zdc_init       zdc_init_
# define zdc_step       zdc_step_
# define zdc_setbeam    zdc_setbeam_
# define zdc_sethijing  zdc_sethijing_
# define zdc_setvenus   zdc_setvenus_
# define zdc_setkine    zdc_setkine_
#else
# define zdc_step       ZDC_STEP
# define zdc_setbeam    ZDC_SETBEAM
# define zdc_sethijing  ZDC_SETHIJING
# define zdc_setvenus   ZDC_SETVENUS
# define zdc_setkine    ZDC_SETKINE
#endif

extern "C" void type_of_call zdc_init();
extern "C" void type_of_call zdc_step();
extern "C" void type_of_call zdc_setbeam(Int_t beam, Float_t fx, Float_t fy,
					 Float_t sx, Float_t sy, Float_t div,
					 Float_t angle, Int_t cross);
extern "C" void type_of_call zdc_sethijing(Int_t hij, Int_t hijf, Int_t hijsp,
					   DEFCHARD DEFCHARL);
extern "C" void type_of_call zdc_setvenus(Int_t hiv, Int_t hivf, Int_t hivsp,
					  DEFCHARD DEFCHARL);
extern "C" void type_of_call zdc_setkine(Int_t code, Float_t pmom, Float_t cx,
					 Float_t cy, Float_t cz, Int_t type,
					 Int_t fermi);
 
ClassImp(AliZDC)
 
//_____________________________________________________________________________
AliZDC::AliZDC()
{
  //
  // Default constructor for the Zero Degree Calorimeter base class
  //
  fIshunt   = 0;
}
 
//_____________________________________________________________________________
AliZDC::AliZDC(const char *name, const char *title)
  : AliDetector(name,title)
{
  //
  // Standard constructor for the Zero Degree Calorimeter base class
  //

  //
  // Allocate the array of hits
  fHits   = new TClonesArray("AliZDChit",  405);
  
  fIshunt     =  1;
}
 
//_____________________________________________________________________________
void AliZDC::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a Zero Degree Calorimeter hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliZDChit(fIshunt,track,vol,hits);
}
 
//_____________________________________________________________________________
void AliZDC::BuildGeometry()
{
  //
  // Build the ROOT TNode geometry for event display 
  // in the Zero Degree Calorimeter
  // This routine is dummy for the moment
  //

  //  TNode *Node, *Top;
  //  TBRIK *brik;
  //  const int kColorZDC  = kRed;
  
  //
  // Top=gAlice->GetGeometry()->GetNode("alice");
  
  // ZDC
  /*
    brik = new TBRIK("S_ZDC","ZDC box","void",300,300,5);
    Top->cd();
    Node = new TNode("ZDC","ZDC","S_ZDC",0,0,600,"");
    Node->SetLineColor(kColorZDC);
    fNodes->Add(Node);
  */
}

//_____________________________________________________________________________
Int_t AliZDC::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Distance from the mouse to the Zero Degree Calorimeter
  // Dummy routine
  //
  return 9999;
}
 
//_____________________________________________________________________________
void AliZDC::SetBeam(Int_t beam, Float_t fx, Float_t fy, Float_t sx,
		     Float_t sy, Float_t div, Float_t angle, Int_t cross)
{
  //
  // Set beam characteristic
  // This routine has to be revised as it is disconnected from the
  // actual generation in this version of AliRoot
  //

  // beam  : 1 = gaussian beam
  //       : 2 = uniform beam
  // fx    : x-coordinate of beam offset
  // fy    : y-coordinate of beam offset
  // sx    : sigma-x of the beam (gaussian or uniform)
  // sy    : sigma-y of the beam (gaussian or uniform)
  // div   : divergency of the beam (32*10**-6 rad for LHC)
  // angle : beam crossing angle (100*10**-6 rad for LHC)
  // cross : 1 = horizontal beam crossing
  //       : 2 = vertical beam crossing
  zdc_setbeam(beam,fx,fy,sx,sy,div,angle,cross);
}

//_____________________________________________________________________________
void AliZDC::SetHijing(Int_t hij, Int_t hijf, Int_t hijsp, const char *file)
{
  //
  // Set the parameter for the HIJING generation
  // This routine has to be revised as it is disconnected from the
  // actual generation in this version of AliRoot
  //

  // HIJ  : 1 = read HIJING event file
  //      : 2 =  "     "      "    "    + debug
  // HIJF : event number of the first event to be read from file
  // HIJSP: 0 = read all particles
  //      : 1 = remove spectator nucleons
  zdc_sethijing(hij,hijf,hijsp, PASSCHARD(file) PASSCHARL(file));
}

//_____________________________________________________________________________
void AliZDC::SetVenus(Int_t hiv, Int_t hivf, Int_t hivsp, const char *file)
{
  //
  // Set the parameter for the VENUS generation
  // This routine has to be revised as it is disconnected from the
  // actual generation in this version of AliRoot
  //

  // HIV  : 1 = read VENUS event file
  //      : 2 =  "     "      "    "    + debug
  // HIVF : event number of the first event to be read from file
  // HIVSP: 0 = read all particles
  //      : 1 = remove spectator nucleons
  zdc_setvenus(hiv,hivf,hivsp, PASSCHARD(file) PASSCHARL(file));
}

//_____________________________________________________________________________
void AliZDC::SetKine(Int_t code, Float_t pmom, Float_t cx, Float_t cy,
		     Float_t cz, Int_t type, Int_t fermi)
{
  //
  // Set the parameter for the event generation
  // This routine has to be revised as it is disconnected from the
  // actual generation in this version of AliRoot
  //

  // code     : GEANT code of the test particle
  // pmom     : absolute value of particle momentum
  // cx,cy,cz : director cosines of the track (if type)
  // type     :  0 = take director cosines from cx,cy,cz
  //          : <>0 = pseudorapidity of the test particle
  // fermi    : 0 = no Fermi motion for the spectator nucleons
  //          : 1 = Fermi motion for the spectator nucleons
  zdc_setkine(code,pmom,cx,cy,cz,type,fermi);
}
 
//_____________________________________________________________________________
void AliZDC::StepManager()
{
  //
  // Routine called at every step in the Zero Degree Calorimeter
  // This is a simple interface to the FORTRAN routine
  // A step manager should be written
  //
  zdc_step();
}

 
ClassImp(AliZDCv1)
 

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Zero Degree Calorimeter version 1                                        //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliZDCv1Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
AliZDCv1::AliZDCv1() : AliZDC()
{
  //
  // Default constructor for Zero Degree Calorimeter
  //
}
 
//_____________________________________________________________________________
AliZDCv1::AliZDCv1(const char *name, const char *title)
  : AliZDC(name,title)
{
  //
  // Standard constructor for Zero Degree Calorimeter 
  //
}
 
//_____________________________________________________________________________
void AliZDCv1::CreateGeometry()
{
  //
  // Create the geometry for the Zero Degree Calorimeter version 1
  // -- Author :    E Scomparin
  //
  //Begin_Html
  /*
    <img src="picts/AliZDCv1.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliZDCv1Tree.gif">
  */
  //End_Html

  // The following variables were illegaly initialized in zdc_init.
  // These variables should become data members of this class
  // once zdc_init has been converted
  //* Initialize COMMON block ZDC_CGEOM
  //*
  
  const Int_t NZPTX=4;
  const Int_t NZPTY=1;
  const Int_t NZNTX=2;
  const Int_t NZNTY=2;
  
  Float_t HDZN[3]  = {4.0,4.0,50.0};
  Float_t HDZP[3]  = {10.0,6.0,75.0};
  // Coordinates of the center of the ZDC front face in the MRS
  Float_t ZNPOS[3] = {-0.5,0.,11613.};
  Float_t ZPPOS[3] = {-21.0,0.,11563.};
  Float_t FIZN[3]  = {0.,0.01825,50.0};
  Float_t FIZP[3]  = {0.,0.01825,75.0};
  Float_t GRZN[3]  = {0.025,0.025,50.0};
  Float_t GRZP[3]  = {0.040,0.040,75.0};
  Int_t NCEN[3]    = {11,11,0};
  Int_t NCEP[3]    = {10,10,0};
  
  Float_t angle;
  Float_t zq, conpar[9], tubpar[3];
  Int_t im1, im2;
  Float_t zd1, zd2;
  
  
  Int_t *idtmed = fIdtmed->GetArray()-799;
  
  // -- Mother of the ZDC 
  conpar[0] = 0.;
  conpar[1] = 360.;
  conpar[2] = 2.;
  conpar[3] = 1920.;
  conpar[4] = 0.;
  conpar[5] = 55.;
  conpar[6] = 13060.;
  conpar[7] = 0.;
  conpar[8] = 55.;
  gMC->Gsvolu("ZDC ", "PCON", idtmed[891], conpar, 9);
  gMC->Gspos("ZDC ", 1, "ALIC", 0., 0., 0., 0, "ONLY");
  // -- FIRST SECTION OF THE BEAM PIPE (from compensator dipole to 
  //    beginning of D1) 
  
  zd1 = 1920.;
  
  tubpar[0] = 2.3;
  tubpar[1] = 2.5;
  tubpar[2] = 1961.75;
  gMC->Gsvolu("P001", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P001", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  //-- SECOND SECTION OF THE BEAM PIPE (FROM THE END OF D1 TO THE BEGINNING OF
  //    D2) 
  
  zd1 = 6316.+472.5;
  
  tubpar[0] = 7.3/2.;
  tubpar[1] = 7.7/2.;
  tubpar[2] = 90.*0.5;
  gMC->Gsvolu("P002", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P002", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 7.3/2.;
  tubpar[1] = 7.7/2.;
  tubpar[2] = 10.*0.5;
  gMC->Gsvolu("P003", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P003", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  conpar[0] = 3.16/2.;
  conpar[1] = 7.3/2.;
  conpar[2] = 7.7/2.;
  conpar[3] = 9.8/2.;
  conpar[4] = 10.0/2.;
  gMC->Gsvolu("P004", "CONE", idtmed[851], conpar, 5);
  gMC->Gspos("P004", 1, "ZDC ", 0., 0., conpar[0] + zd1, 0, "ONLY");
  
  zd1 += conpar[0] * 2.;
  
  tubpar[0] = 9.8/2.;
  tubpar[1] = 10.0/2;
  tubpar[2] = 490./2.;
  gMC->Gsvolu("P005", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P005", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  conpar[0] = 30./2.;
  conpar[1] = 9.8/2.;
  conpar[2] = 10.0/2.;
  conpar[3] = 20.4/2.;
  conpar[4] = 20.6/2.;
  gMC->Gsvolu("P006", "CONE", idtmed[851], conpar, 5);
  gMC->Gspos("P006", 1, "ZDC ", 0., 0., conpar[0] + zd1, 0, "ONLY");
  
  zd1 += conpar[0] * 2.;
  
  tubpar[0] = 20.4/2.;
  tubpar[1] = 20.6/2.;
  tubpar[2] = 150./2.;
  gMC->Gsvolu("P007", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P007", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  conpar[0] = 13.6/2.;
  conpar[1] = 20.4/2.;
  conpar[2] = 20.6/2.;
  conpar[3] = 25.2/2.;
  conpar[4] = 25.4/2.;
  gMC->Gsvolu("P008", "CONE", idtmed[851], conpar, 5);
  gMC->Gspos("P008", 1, "ZDC ", 0., 0., conpar[0] + zd1, 0, "ONLY");
  
  zd1 += conpar[0] * 2.;
  
  tubpar[0] = 25.2/2.;
  tubpar[1] = 25.4/2.;
  tubpar[2] = 205.8/2.;
  gMC->Gsvolu("P009", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P009", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 43.8/2.;
  tubpar[1] = 44.0/2.;
  tubpar[2] = 500./2.;
  gMC->Gsvolu("P010", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P010", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 31.8/2.;
  tubpar[1] = 32.0/2.;
  tubpar[2] = 757.5/2.;
  gMC->Gsvolu("P011", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P011", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  conpar[0] = 22.7/2.;
  conpar[1] = 31.8/2.;
  conpar[2] = 32.0/2.;
  conpar[3] = 39.8/2.;
  conpar[4] = 40.0/2.;
  gMC->Gsvolu("P012", "CONE", idtmed[851], conpar, 5);
  gMC->Gspos("P012", 1, "ZDC ", 0., 0., conpar[0] + zd1, 0, "ONLY");
  
  zd1 += conpar[0] * 2.;
  
  tubpar[0] = 39.8/2.;
  tubpar[1] = 40.0/2.;
  tubpar[2] = 100./2.;
  gMC->Gsvolu("P013", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P013", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 39.8/2.;
  tubpar[1] = 40.0/2.;
  tubpar[2] = 600./2.;
  gMC->Gsvolu("P014", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P014", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  conpar[0] = 28.4/2.;
  conpar[1] = 39.8/2.;
  conpar[2] = 40.0/2.;
  conpar[3] = 49.8/2.;
  conpar[4] = 50.0/2.;
  gMC->Gsvolu("P015", "CONE", idtmed[851], conpar, 5);
  gMC->Gspos("P015", 1, "ZDC ", 0., 0., conpar[0] + zd1, 0, "ONLY");
  
  zd1 += conpar[0] * 2.;
  
  tubpar[0] = 49.8/2.;
  tubpar[1] = 50.0/2.;
  tubpar[2] = 100./2.;
  gMC->Gsvolu("P016", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P016", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 49.8/2.;
  tubpar[1] = 50.0/2.;
  tubpar[2] = 600./2.;
  gMC->Gsvolu("P017", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P017", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  conpar[0] = 28.4/2.;
  conpar[1] = 49.8/2.;
  conpar[2] = 50.0/2.;
  conpar[3] = 59.8/2.;
  conpar[4] = 60.0/2.;
  gMC->Gsvolu("P018", "CONE", idtmed[851], conpar, 5);
  gMC->Gspos("P018", 1, "ZDC ", 0., 0., conpar[0] + zd1, 0, "ONLY");
  
  zd1 += conpar[0] * 2.;
  
  tubpar[0] = 59.8/2.;
  tubpar[1] = 60.0/2.;
  tubpar[2] = 50./2.;
  gMC->Gsvolu("P019", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P019", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 59.8/2.;
  tubpar[1] = 60.0/2.;
  tubpar[2] = 800./2.;
  gMC->Gsvolu("P020", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P020", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 0.;
  tubpar[1] = 60.0/2.;
  tubpar[2] = 0.2/2.;
  gMC->Gsvolu("P021", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("P021", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  
  zd1 += tubpar[2] * 2.;
  
  tubpar[0] = 0.;
  tubpar[1] = 4.4/2.;
  tubpar[2] = 0.2/2.;
  gMC->Gsvolu("Q021", "TUBE", idtmed[889], tubpar, 3);
  tubpar[0] = 0.;
  tubpar[1] = 7.0/2.;
  tubpar[2] = 0.2/2.;
  gMC->Gsvolu("R021", "TUBE", idtmed[889], tubpar, 3);
  // -- POSITION Q021 INSIDE P021 
  gMC->Gspos("Q021", 1, "P021", -7.7, 0., 0., 0, "ONLY");
  // -- POSITION R020 INSIDE P020 
  gMC->Gspos("R021", 1, "P021", 7.7, 0., 0., 0, "ONLY");
  
  // -- BEAM PIPES BETWEEN END OF CONICAL PIPE AND BEGINNING OF D2 
  tubpar[0] = 4.0/2.;
  tubpar[1] = 4.4/2.;
  tubpar[2] = 645.*0.5;
  gMC->Gsvolu("P022", "TUBE", idtmed[851], tubpar, 3);
  tubpar[0] = 7.0/2.;
  tubpar[1] = 7.4/2.;
  tubpar[2] = 645.*0.5;
  gMC->Gsvolu("P023", "TUBE", idtmed[851], tubpar, 3);
  
  // -- ROTATE PIPES 
  AliMatrix(im1, 90.-0.071, 0., 90., 90., .071, 180.);
  angle = .071*kDegrad;
  gMC->Gspos("P022", 1, "ZDC ", TMath::Sin(angle) * 322.5 - 9.7 + 
	       TMath::Sin(angle) * 472.5, 0., tubpar[2] + zd1, im1, "ONLY");
  AliMatrix(im2, 90.+0.071, 0., 90., 90., .071, 0.);
  gMC->Gspos("P023", 1, "ZDC ", 9.7 - TMath::Sin(angle) * 322.5, 0., 
	       tubpar[2] + zd1, im2, "ONLY");
  
  // --  END OF BEAM PIPE VOLUME DEFINITION. MAGNET DEFINITION FOLLOWS 
  //     (LHC OPTICS 6) 
  
  // -- COMPENSATOR DIPOLE (MCBWA) 
  //     GAP (VACUUM WITH MAGNETIC FIELD) 
  
  tubpar[0] = 0.;
  tubpar[1] = 4.5;
  tubpar[2] = 190./2.;
  gMC->Gsvolu("MCBW", "TUBE", idtmed[890], tubpar, 3);
  gMC->Gspos("MCBW", 1, "ZDC ", 0., 0., tubpar[2] + 1920., 0, "ONLY");
  
  // --  YOKE (IRON WITHOUT MAGNETIC FIELD) 
  
  tubpar[0] = 4.5;
  tubpar[1] = 55.;
  tubpar[2] = 190./2.;
  gMC->Gsvolu("YMCB", "TUBE", idtmed[851], tubpar, 3);
  gMC->Gspos("YMCB", 1, "ZDC ", 0., 0., tubpar[2] + 1920., 0, "ONLY");
  
  // -- INNER TRIPLET 
  
  zq = 2300.;
  
  // -- DEFINE MQXL AND MQX QUADRUPOLE ELEMENT 
  
  //     MQXL 
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  
  tubpar[0] = 0.;
  tubpar[1] = 3.5;
  tubpar[2] = 630./2.;
  gMC->Gsvolu("MQXL", "TUBE", idtmed[890], tubpar, 3);
  
  // --  YOKE 
  
  tubpar[0] = 3.5;
  tubpar[1] = 22.;
  tubpar[2] = 630./2.;
  gMC->Gsvolu("YMQL", "TUBE", idtmed[851], tubpar, 3);
  
  gMC->Gspos("MQXL", 1, "ZDC ", 0., 0., tubpar[2] + zq, 0, "ONLY");
  gMC->Gspos("YMQL", 1, "ZDC ", 0., 0., tubpar[2] + zq, 0, "ONLY");
  
  gMC->Gspos("MQXL", 2, "ZDC ", 0., 0., tubpar[2] + zq + 2430., 0, "ONLY");
  gMC->Gspos("YMQL", 2, "ZDC ", 0., 0., tubpar[2] + zq + 2430., 0, "ONLY");
  
  // --  MQX 
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  
  tubpar[0] = 0.;
  tubpar[1] = 3.5;
  tubpar[2] = 550./2.;
  gMC->Gsvolu("MQX ", "TUBE", idtmed[890], tubpar, 3);
  
  // --  YOKE 
  
  tubpar[0] = 3.5;
  tubpar[1] = 22.;
  tubpar[2] = 550./2.;
  gMC->Gsvolu("YMQ ", "TUBE", idtmed[851], tubpar, 3);
  
  gMC->Gspos("MQX ", 1, "ZDC ", 0., 0., tubpar[2] + zq + 880.,  0, "ONLY");
  gMC->Gspos("YMQ ", 1, "ZDC ", 0., 0., tubpar[2] + zq + 880.,  0, "ONLY");
  
  gMC->Gspos("MQX ", 2, "ZDC ", 0., 0., tubpar[2] + zq + 1530., 0, "ONLY");
  gMC->Gspos("YMQ ", 2, "ZDC ", 0., 0., tubpar[2] + zq + 1530., 0, "ONLY");
  
  // -- SEPARATOR DIPOLE D1 
  
  zd1 = 5843.5;
  
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  
  tubpar[0] = 0.;
  tubpar[1] = 4.5;
  tubpar[2] = 945/2.;
  gMC->Gsvolu("D1  ", "TUBE", idtmed[890], tubpar, 3);
  
  // --  YOKE 
  
  tubpar[0] = 0.;
  tubpar[1] = 55.;
  tubpar[2] = 945/2.;
  gMC->Gsvolu("YD1 ", "TUBE", idtmed[851], tubpar, 3);
  
  gMC->Gspos("YD1 ", 1, "ZDC ", 0., 0., tubpar[2] + zd1, 0, "ONLY");
  gMC->Gspos("D1  ", 1, "YD1 ", 0., 0., 0., 0, "ONLY");
  
  // -- DIPOLE D2 
  
  zd2 = 12113.2;
  
  // --  GAP (VACUUM WITH MAGNETIC FIELD) 
  
  tubpar[0] = 0.;
  tubpar[1] = 4.5;
  tubpar[2] = 945./2.;
  gMC->Gsvolu("D2  ", "TUBE", idtmed[890], tubpar, 3);
  
  // --  YOKE 
  
  tubpar[0] = 0.;
  tubpar[1] = 55.;
  tubpar[2] = 945./2.;
  gMC->Gsvolu("YD2 ", "TUBE", idtmed[851], tubpar, 3);
  
  gMC->Gspos("YD2 ", 1, "ZDC ", 0., 0., tubpar[2] + zd2, 0, "ONLY");
  
  gMC->Gspos("D2  ", 1, "YD2 ", -9.7, 0., 0., 0, "ONLY");
  gMC->Gspos("D2  ", 2, "YD2 ",  9.7, 0., 0., 0, "ONLY");
  
  // -- END OF MAGNET DEFINITION 
  
  // ----------------- Hadronic calorimeters -------------------- * 
  
  // Neutron calorimeter 
  
  gMC->Gsvolu("ZNEU", "BOX ", idtmed[800], HDZN, 3); // Passive material  
  gMC->Gsvolu("ZNFI", "TUBE", idtmed[802], FIZN, 3); // Active material
  gMC->Gsvolu("ZNGR", "BOX ", idtmed[889], GRZN, 3); // Empty grooves 
  
  // Divide ZNEU in towers 
  // (for hits purposes) 
  
  gMC->Gsdvn("ZNTX", "ZNEU", NZNTX, 1); // x-tower 
  gMC->Gsdvn("ZN1 ", "ZNTX", NZNTY, 2); // y-tower
  
  // Divide ZNEU in minitowers 
  // (NCEN(1)= NUMBER OF FIBERS PER TOWER ALONG X-AXIS, 
  //  NCEN(2)= NUMBER OF FIBERS PER TOWER ALONG Y-AXIS) 
  // (one fiber per minitower) 
  
  gMC->Gsdvn("ZNSL", "ZN1 ", NCEN[1], 2); // Slices 
  gMC->Gsdvn("ZNST", "ZNSL", NCEN[0], 1); // Sticks
  
  // --- Position the empty grooves in the sticks 
  gMC->Gspos("ZNGR", 1, "ZNST", 0., 0., 0., 0, "ONLY");
  // --- Position the fibers in the grooves 
  gMC->Gspos("ZNFI", 1, "ZNGR", 0., 0., 0., 0, "ONLY");
  // --- Position the neutron calorimeter in ZDC 
  gMC->Gspos("ZNEU", 1, "ZDC ", ZNPOS[0], ZNPOS[1], ZNPOS[2] + HDZN[2], 0, "ONLY");
  
  // Proton calorimeter 
  
  gMC->Gsvolu("ZPRO", "BOX ", idtmed[801], HDZP, 3); // Passive material
  gMC->Gsvolu("ZPFI", "TUBE", idtmed[802], FIZP, 3); // Active material 
  gMC->Gsvolu("ZPGR", "BOX ", idtmed[889], GRZP, 3); // Empty grooves
  
  // Divide ZPRO in towers 
  // (for hits purposes) 
  
  gMC->Gsdvn("ZPTX", "ZPRO", NZPTX, 1); // x-tower 
  gMC->Gsdvn("ZP1 ", "ZPTX", NZPTY, 2); // y-tower
  
  
  // Divide ZPRO in minitowers 
  // (NCEP[0]= NUMBER OF FIBERS ALONG X-AXIS PER MINITOWER, 
  //  NCEP[1]= NUMBER OF FIBERS ALONG Y-AXIS PER MINITOWER) 
  // (one fiber per minitower) 
  
  gMC->Gsdvn("ZPSL", "ZP1 ", NCEP[1], 2); // Slices 
  gMC->Gsdvn("ZPST", "ZPSL", NCEP[0], 1); // Sticks
  
  // --- Position the empty grooves in the sticks 
  gMC->Gspos("ZPGR", 1, "ZPST", 0., 0., 0., 0, "ONLY");
  // --- Position the fibers in the grooves 
  gMC->Gspos("ZPFI", 1, "ZPGR", 0., 0., 0., 0, "ONLY");
  // --- Position the proton calorimeter in ZDC 
  gMC->Gspos("ZPRO", 1, "ZDC ", ZPPOS[0], ZPPOS[1], ZPPOS[2] + HDZP[2], 0, "ONLY");
  
}
 
//_____________________________________________________________________________
void AliZDCv1::DrawModule()
{
  //
  // Draw a shaded view of the Zero Degree Calorimeter version 1
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("ZDC","SEEN",0);
  gMC->Gsatt("P001","SEEN",1);
  gMC->Gsatt("P002","SEEN",1);
  gMC->Gsatt("P003","SEEN",1);
  gMC->Gsatt("P004","SEEN",1);
  gMC->Gsatt("P005","SEEN",1);
  gMC->Gsatt("P006","SEEN",1);
  gMC->Gsatt("P007","SEEN",1);
  gMC->Gsatt("P008","SEEN",1);
  gMC->Gsatt("P009","SEEN",1);
  gMC->Gsatt("P010","SEEN",1);
  gMC->Gsatt("P011","SEEN",1);
  gMC->Gsatt("P012","SEEN",1);
  gMC->Gsatt("P013","SEEN",1);
  gMC->Gsatt("P014","SEEN",1);
  gMC->Gsatt("P015","SEEN",1);
  gMC->Gsatt("P016","SEEN",1);
  gMC->Gsatt("P017","SEEN",1);
  gMC->Gsatt("P018","SEEN",1);
  gMC->Gsatt("P019","SEEN",1);
  gMC->Gsatt("P020","SEEN",1);
  gMC->Gsatt("P021","SEEN",1);
  gMC->Gsatt("Q021","SEEN",1);
  gMC->Gsatt("R021","SEEN",1);
  gMC->Gsatt("P022","SEEN",1);
  gMC->Gsatt("P023","SEEN",1);
  gMC->Gsatt("D1  ","SEEN",1);
  gMC->Gsatt("YD1 ","SEEN",1);
  gMC->Gsatt("D2  ","SEEN",1);
  gMC->Gsatt("YD2 ","SEEN",1);
  gMC->Gsatt("MCBW","SEEN",1);
  gMC->Gsatt("YMCB","SEEN",1);
  gMC->Gsatt("MQXL","SEEN",1);
  gMC->Gsatt("YMQL","SEEN",1);
  gMC->Gsatt("MQX","SEEN",1);
  gMC->Gsatt("YMQ","SEEN",1);
  gMC->Gsatt("D1","SEEN",1);
  gMC->Gsatt("YD1","SEEN",1);
  gMC->Gsatt("D2","SEEN",1);
  gMC->Gsatt("YD2","SEEN",1);
  gMC->Gsatt("ZNEU","SEEN",0);
  gMC->Gsatt("ZNFI","SEEN",0);
  gMC->Gsatt("ZNGR","SEEN",0);
  gMC->Gsatt("ZNTX","SEEN",0);
  gMC->Gsatt("ZN1 ","COLO",2); 
  gMC->Gsatt("ZN1 ","SEEN",1);
  gMC->Gsatt("ZNSL","SEEN",0);
  gMC->Gsatt("ZNST","SEEN",0);
  gMC->Gsatt("ZPRO","SEEN",0);
  gMC->Gsatt("ZPFI","SEEN",0);
  gMC->Gsatt("ZPGR","SEEN",0);
  gMC->Gsatt("ZPTX","SEEN",0);
  gMC->Gsatt("ZP1 ","SEEN",1);
  gMC->Gsatt("ZPSL","SEEN",0);
  gMC->Gsatt("ZPST","SEEN",0);
  
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 100, -100, 100, 12000, 16000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 488, 220, .07, .07);
  gMC->Gdhead(1111, "Zero Degree Calorimeter Version 1");
  gMC->Gdman(18, 4, "MAN");
}

//_____________________________________________________________________________
void AliZDCv1::CreateMaterials()
{
  //
  // Create Materials for the Zero Degree Calorimeter
  //
  // Origin    : E. Scomparin 
  
  Int_t *idtmed = fIdtmed->GetArray()-799;
  
  Float_t dens, ubuf[1], wmat[2];
  Int_t isvol_active;
  Float_t a[2];
  Int_t i;
  Float_t z[2], epsil=0.001, stmin=0.01;
  Int_t isvol;
  Float_t fieldm = gAlice->Field()->Max();
  Int_t inofld;
  Float_t deemax=-1;
  Float_t tmaxfd=gAlice->Field()->Max();
  Int_t isxfld = gAlice->Field()->Integ();
  Float_t stemax;
  
  // --- Store in UBUF r0 for nuclear radius calculation R=r0*A**1/3 
  
  // --- Tungsten 
  ubuf[0] = 1.11;
  AliMaterial(1, "TUNG", 183.85, 74., 19.3, .35, 10.3, ubuf, 1);
  
  // --- Brass (CuZn) 
  dens = 8.48;
  a[0] = 63.546;
  a[1] = 65.39;
  z[0] = 29.;
  z[1] = 30.;
  wmat[0] = .63;
  wmat[1] = .37;
  AliMixture(2, "BRASS               ", a, z, dens, 2, wmat);
  
  // --- SiO2 
  dens = 2.64;
  a[0] = 28.086;
  a[1] = 15.9994;
  z[0] = 14.;
  z[1] = 8.;
  wmat[0] = 1.;
  wmat[1] = 2.;
  AliMixture(3, "SIO2                ", a, z, dens, -2, wmat);
  
  // --- Lead 
  ubuf[0] = 1.12;
  AliMaterial(4, "LEAD", 207.19, 82., 11.35, .56, 18.5, ubuf, 1);
  
  // --- Copper 
  ubuf[0] = 1.1;
  AliMaterial(5, "COPP", 63.54, 29., 8.96, 1.4, 0., ubuf, 1);
  
  // --- Tantalum 
  ubuf[0] = 1.1;
  AliMaterial(6, "TANT", 180.95, 73., 16.65, .4, 11.9, ubuf, 1);
  
  // Steel still to be added 
  
  // --- Iron 
  ubuf[0] = 1.1;
  AliMaterial(52, "IRON", 55.85, 26., 7.87, 1.76, 0., ubuf, 1);
  
  // --- Vacuum (no magnetic field) 
  AliMaterial(90, "VOID", 1e-16, 1e-16, 1e-16, 1e16, 1e16, ubuf,0);
  
  // --- Vacuum (magnetic field) 
  AliMaterial(91, "VOIM", 1e-16, 1e-16, 1e-16, 1e16, 1e16, ubuf,0);
  
  // --- Air non magnetic 
  AliMaterial(92, "Air    $", 14.61, 7.3, .001205, 30420., 67500., ubuf, 0);
  
  // ---  Definition of tracking media: 
  
  // --- Tungsten = 801 ; 
  // --- Brass = 802 ; 
  // --- Fibers (SiO2) = 803 ; 
  // --- Lead = 804 ; 
  // --- Copper = 805 ; 
  // --- Tantalum = 806 ; 
  // --- Steel = 851 ; 
  // --- Iron = 852 ; 
  // --- Vacuum (no field) = 890 
  // --- Vacuum (with field) = 891 
  // --- Air   (no field) = 892 
  
  
  // --- Tracking media parameters 
  epsil  = .01;
  stemax = 1.;
  isvol  = 0;
  isvol_active = 1;
  inofld = 0;
  fieldm = 0.;
  
  AliMedium(1, "ZW", 1, isvol_active, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(2, "ZBRASS", 2, isvol_active, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(3, "ZSIO2", 3, isvol_active, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(4, "ZLEAD", 4, isvol_active, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(5, "ZCOPP", 5, isvol_active, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(6, "ZTANT", 6, isvol_active, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(52, "ZIRON", 52, isvol, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(90, "ZVOID", 90, isvol, inofld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(92, "Air", 92, 0, inofld, fieldm, tmaxfd, stemax,deemax, epsil, stmin);
  
  fieldm = 45.;
  //    AliMedium(91, "ZVOIM", 91, isvol, isxfld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(91, "ZVOIM", 91, isvol, isxfld, fieldm, tmaxfd, stemax, deemax, epsil, stmin);
  
  // Thresholds for showering in the ZDCs 
  
  i = 801;
  gMC->Gstpar(idtmed[i-1], "CUTGAM", .01);
  gMC->Gstpar(idtmed[i-1], "CUTELE", .01);
  gMC->Gstpar(idtmed[i-1], "CUTNEU", .1);
  gMC->Gstpar(idtmed[i-1], "CUTHAD", .1);
  i = 802;
  gMC->Gstpar(idtmed[i-1], "CUTGAM", .01);
  gMC->Gstpar(idtmed[i-1], "CUTELE", .01);
  gMC->Gstpar(idtmed[i-1], "CUTNEU", .1);
  gMC->Gstpar(idtmed[i-1], "CUTHAD", .1);
  
  // Avoid too detailed showering along the beam line 
  
  i = 852;
  gMC->Gstpar(idtmed[i-1], "CUTGAM", .1);
  gMC->Gstpar(idtmed[i-1], "CUTELE", .1);
  gMC->Gstpar(idtmed[i-1], "CUTNEU", 1.);
  gMC->Gstpar(idtmed[i-1], "CUTHAD", 1.);
  
  // Avoid interaction in fibers (only energy loss allowed) 
  i = 803;
  gMC->Gstpar(idtmed[i-1], "DCAY", 0.);
  gMC->Gstpar(idtmed[i-1], "MULS", 0.);
  gMC->Gstpar(idtmed[i-1], "PFIS", 0.);
  gMC->Gstpar(idtmed[i-1], "MUNU", 0.);
  gMC->Gstpar(idtmed[i-1], "LOSS", 1.);
  gMC->Gstpar(idtmed[i-1], "PHOT", 0.);
  gMC->Gstpar(idtmed[i-1], "COMP", 0.);
  gMC->Gstpar(idtmed[i-1], "PAIR", 0.);
  gMC->Gstpar(idtmed[i-1], "BREM", 0.);
  gMC->Gstpar(idtmed[i-1], "DRAY", 0.);
  gMC->Gstpar(idtmed[i-1], "ANNI", 0.);
  gMC->Gstpar(idtmed[i-1], "HADR", 0.);
}

ClassImp(AliZDChit)
  
//_____________________________________________________________________________
AliZDChit::AliZDChit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
  //
  // Add a Zero Degree Calorimeter hit
  //
  Int_t i;
  for (i=0;i<4;i++) fVolume[i] = vol[i];
  fX=hits[0];
  fY=hits[1];
  fZ=hits[2];
  fEnergy=hits[3];
}
