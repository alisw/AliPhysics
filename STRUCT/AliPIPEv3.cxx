///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Beam pipe class                                                          //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliPIPEClass.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliPIPEv3.h"
#include "AliRun.h"
#include "AliMC.h"
 
ClassImp(AliPIPEv3)
 
//_____________________________________________________________________________
AliPIPEv3::AliPIPEv3()
{
  //
  // Default constructor for beam pipe
  //
}
 
//_____________________________________________________________________________
AliPIPEv3::AliPIPEv3(const char *name, const char *title)
       : AliPIPE(name,title)
{
  //
  // Standard constructor for beam pipe
  //
}
 
//_____________________________________________________________________________
void AliPIPEv3::CreateGeometry()
{
  //
  // Create Beam Pipe geometry
  //
  //Begin_Html
  /*
    <img src="picts/AliPIPE.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliPIPETree.gif">
  */
  //End_Html

  AliMC* pMC = AliMC::GetMC();
  
  Float_t tpar[3], dzmo, zpos;
  Float_t bepar[3], alpar[3],sspar[3],flange[3],vacpar[3];
  Float_t bellow[3];
//  Float_t undul[3];
//  const Double_t z_flange = 150;
//for undulated structure
  char cn18[][5]={"CN01","CN02","CN03","CN04","CN05","CN06","CN07","CN08"};
  char cn48[][5]={"CN21","CN22","CN23","CN24","CN25","CN26","CN27","CN28"};
//  char undul[][5]={'BELO','UNDL'};
  Float_t zundul;
  Float_t rundul;
  Float_t pitch;
  Float_t thick;

  
  Int_t *idtmed = gAlice->Idtmed();
//     the mother of all beam pipes

  tpar[0] = 0.;
  tpar[1] = 10.;
  tpar[2] = 1400. / 2;
  dzmo = tpar[2];
  pMC->Gsvolu("QQMO", "TUBE", idtmed[2015], tpar, 3);

//	All beam pipe details as per the provisonal drawings given by Lars
//	Leistam on 31.5.99 
    
//     Beryllium  beam pipe, length 56.6 cm, centered at vertex 
  
  bepar[0]=0.0;
  bepar[1]=3.0;
  bepar[2]=28.3;
  zpos=0.0;
  vacpar[0]=0.0;
  vacpar[1]=2.9;
  vacpar[2]=bepar[2];
  //
  pMC->Gsvolu("QQBE", "TUBE", idtmed[2004], bepar, 3);
  pMC->Gsvolu("VAC1", "TUBE", idtmed[2015], vacpar, 3);
  pMC->Gspos("VAC1", 1, "QQBE", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("QQBE", 1, "QQMO", 0., 0., zpos, 0, "ONLY");
  
  // now beam pipes only in negative z-part for use in PMD.
 
  // SS Flange 4 cm thick, 5.8 cm ID, 6.3 cm OD
  flange[0]=0.0;
  flange[1]=3.15;
  flange[2]=2.0;
  zpos = zpos -bepar[2] - flange[2];
  vacpar[0]=0.0;
  vacpar[1]=2.9;
  vacpar[2]=flange[2];
  //
  pMC->Gsvolu("QFL1", "TUBE", idtmed[2018], flange, 3);
  pMC->Gsvolu("VAC2", "TUBE", idtmed[2015], vacpar, 3);
  pMC->Gspos("VAC2", 1, "QFL1", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("QFL1", 1, "QQMO", 0., 0., zpos, 0, "ONLY");
  
  // Aluminium alloy beam pipe, 1mm thick, 230 cm long
  alpar[0]=0.0;
  alpar[1]=3.0;
  alpar[2]=115.;
  zpos = zpos - flange[2] - alpar[2];

  vacpar[0]=0.0;
  vacpar[1]=2.9;
  vacpar[2]=alpar[2];
  pMC->Gsvolu("QQAL", "TUBE", idtmed[2003], alpar, 3);
  pMC->Gsvolu("VAC3", "TUBE", idtmed[2015], vacpar, 3);
  pMC->Gspos("VAC3", 1, "QQAL", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("QQAL", 1, "QQMO", 0., 0., zpos, 0, "ONLY");

 
  // SS tube 2.0 cm long, 0.8 mm thick, 5.96 cm OD

  sspar[0]=0.0;
  sspar[1]=2.98;
  sspar[2]=1.0;
  zpos = zpos - alpar[2] - sspar[2];

  vacpar[0]=0.0;
  vacpar[1]=2.9;
  vacpar[2]=sspar[2];
  pMC->Gsvolu("QSS1", "TUBE", idtmed[2018], sspar, 3);
  pMC->Gsvolu("VAC4", "TUBE", idtmed[2015], vacpar, 3);
  pMC->Gspos("VAC4", 1, "QSS1", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("QSS1", 1, "QQMO", 0., 0., zpos, 0, "ONLY");


 // SS Flange 3 cm thick 7.4 cm OD, 5.8 cm ID
   
  flange[0]=0.0;
  flange[1]=3.7;
  flange[2]=1.5;
  zpos = zpos - sspar[2] - flange[2];

  vacpar[0]=0.0;
  vacpar[1]=2.9;
  vacpar[2]=flange[2];
  pMC->Gsvolu("QFL2", "TUBE", idtmed[2018], flange, 3);
  pMC->Gsvolu("VAC5", "TUBE", idtmed[2015], vacpar, 3);
  pMC->Gspos("VAC5", 1, "QFL2", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("QFL2", 1, "QQMO", 0., 0., zpos, 0, "ONLY");


  // SS tube 4.0 cm long, 0.8 mm thick, 5.96 cm OD

  sspar[0]=0.0;
  sspar[1]=2.98;
  sspar[2]=2.0;
  zpos = zpos - flange[2] - sspar[2];

  vacpar[0]=0.0;
  vacpar[1]=2.9;
  vacpar[2]=sspar[2];
  pMC->Gsvolu("QSS2", "TUBE", idtmed[2018], sspar, 3);
  pMC->Gsvolu("VAC6", "TUBE", idtmed[2015], vacpar, 3);
  pMC->Gspos("VAC6", 1, "QSS2", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("QSS2", 1, "QQMO", 0., 0., zpos, 0, "ONLY");


  // *************
  // SS Bellow 8.4 cm long, 6.5 cm ID, 7.5 cm OD
  // 0.8 mm thick material, 0.3 cm pitch.
  // zundul=4.2, rundul=6.5, thick=0.08
  // **************
  pitch=0.3;
  thick=0.08;
  zundul=4.2;
  rundul=6.5;
  Undulation("BELO",pitch,thick,zundul,rundul,cn18);
//
  bellow[2]=zundul;
  zpos = zpos - sspar[2] - bellow[2];
  pMC->Gspos("BELO", 1, "QQMO", 0., 0., zpos, 0, "ONLY");

  // SS tube 20.0 cm long, 0.8 mm thick, 5.96 cm OD

  sspar[0]=0.0;
  sspar[1]=2.98;
  sspar[2]=10.0;
  zpos = zpos - bellow[2] - sspar[2];

  vacpar[0]=0.0;
  vacpar[1]=2.9;
  vacpar[2]=sspar[2];
  pMC->Gsvolu("QSS3", "TUBE", idtmed[2018], sspar, 3);
  pMC->Gsvolu("VAC7", "TUBE", idtmed[2015], vacpar, 3);
  pMC->Gspos("VAC7", 1, "QSS3", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("QSS3", 1, "QQMO", 0., 0., zpos, 0, "ONLY");

  // *************
  // SS Bellow 8.4 cm long, 6.5 cm ID, 7.5 cm OD
  // 0.8 mm thick material, 0.3 cm pitch.
  // **************
//  
  zpos = zpos - sspar[2] - bellow[2];
  pMC->Gspos("BELO", 2, "QQMO", 0., 0., zpos, 0, "ONLY");

  // SS tube 4.7 cm long, 0.8 mm thick, 

  sspar[0]=0.0;
  sspar[1]=2.98;
  sspar[2]=4.7/2.;
  zpos = zpos - bellow[2] - sspar[2];

  vacpar[0]=0.0;
  vacpar[1]=2.9;
  vacpar[2]=sspar[2];
  pMC->Gsvolu("QSS4", "TUBE", idtmed[2018], sspar, 3);
  pMC->Gsvolu("VAC8", "TUBE", idtmed[2015], vacpar, 3);
  pMC->Gspos("VAC8", 1, "QSS4", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("QSS4", 1, "QQMO", 0., 0., zpos, 0, "ONLY");

  // SS Flange 2.2 cm thick, ID=5.8 cm, OD=9.8 cm

  flange[0]=0.0;
  flange[1]=4.9;
  flange[2]=1.1;
  zpos = zpos - sspar[2] - flange[2];

  vacpar[0]=0.0;
  vacpar[1]=2.9;
  vacpar[2]=flange[2];
  pMC->Gsvolu("QFL3", "TUBE", idtmed[2018], flange, 3);
  pMC->Gsvolu("VAC9", "TUBE", idtmed[2015], vacpar, 3);
  pMC->Gspos("VAC9", 1, "QFL3", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("QFL3", 1, "QQMO", 0., 0., zpos, 0, "ONLY");

//Total of 3150 mm from vertex on the negative side upto this point.

// SS tube 20.0 cm long, 0.15 cm thick, 5.8 cm ID, to support vac. pump

  sspar[0]=0.0;
  sspar[1]=3.05;
  sspar[2]=10.0;
  zpos = zpos - flange[2] - sspar[2];

  vacpar[0]=0.0;
  vacpar[1]=2.9;
  vacpar[2]=sspar[2];
  pMC->Gsvolu("QSS5", "TUBE", idtmed[2018], sspar, 3);
  pMC->Gsvolu("VA10", "TUBE", idtmed[2015], vacpar, 3);
  pMC->Gspos("VA10", 1, "QSS5", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("QSS5", 1, "QQMO", 0., 0., zpos, 0, "ONLY");

// 
  // last item, undulated SS beam pipe, pitch=0.25, length= 342.0 cm
  // material thickness 0.015 cm, ID=6.0 cm,
  // zundul=171.0, thick=0.015, rundul=3.0
  pitch=0.25;
  thick=0.015;
  zundul=171;
  rundul=3.0;
  Undulation("UNDL",pitch,thick,zundul,rundul,cn48);
  //
  zpos = zpos - sspar[2] - zundul;
  pMC->Gspos("UNDL", 1, "QQMO", 0., 0., zpos, 0, "ONLY");
//
  pMC->Gspos("QQMO", 1, "ALIC", 0., 0., 0.1, 0, "ONLY");

// 	total of 6770 mm length upto this point, end of undulated beam
//	pipe section.

//	SS flange 22*2 mm thick


  flange[0]=0.0;
  flange[1]=6.3;
  flange[2]=2.2;
  zpos = zpos  - zundul - flange[2];

  vacpar[0]=0.0;
  vacpar[1]=2.9;
  vacpar[2]=flange[2];
  pMC->Gsvolu("QFL4", "TUBE", idtmed[2018], flange, 3);
  pMC->Gsvolu("VC11", "TUBE", idtmed[2015], vacpar, 3);
  pMC->Gspos("VC11", 1, "QFL4", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("QFL4", 1, "QQMO", 0., 0., zpos, 0, "ONLY");

}

//_____________________________________________________________________________
void AliPIPEv3::DrawModule()
{  
  //
  // Draw a shaded view of the Beam Pipe
  //

  AliMC* pMC = AliMC::GetMC();

  // Set everything unseen
  pMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  pMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  pMC->Gsatt("QQMO","seen",1);
  pMC->Gsatt("QQBE","seen",1);
  pMC->Gsatt("QFL1","seen",1);
  pMC->Gsatt("QQAL","seen",1);
  pMC->Gsatt("QSS1","seen",1);
  pMC->Gsatt("QFL2","seen",1);
  pMC->Gsatt("QSS2","seen",1);
  pMC->Gsatt("QSS3","seen",1);
  pMC->Gsatt("QSS4","seen",1);
  pMC->Gsatt("QFL3","seen",1);
  pMC->Gsatt("QSS5","seen",1);
  pMC->Gsatt("BELO","seen",1);
  pMC->Gsatt("UNDL","seen",1);
  //
  pMC->Gdopt("hide", "on");
  pMC->Gdopt("shad", "on");
  pMC->Gsatt("*", "fill", 7);
  pMC->SetClipBox(".");
  pMC->SetClipBox("*", 0, 3000, -3000, 3000, -6000, 6000);
  pMC->DefaultRange();
  pMC->Gdraw("alic", 40, 30, 0, 3, 5, .04, .04);
  pMC->Gdhead(1111, "Beam Pipe");
  pMC->Gdman(16, 6, "MAN");
  pMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliPIPEv3::CreateMaterials()
{
  //
  // Create materials for beam pipe
  //

  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  
  Float_t epsil, stmin, tmaxfd, deemax, stemax;
  
  //     STEEL 
  
  
  // --- Define the various materials for GEANT --- 
  AliMaterial(5, "BERILLIUM$", 9.01, 4., 1.848, 35.3, 36.7);
  AliMaterial(4, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 18.5);
  AliMaterial(16, "VACUUM$ ", 1e-16, 1e-16, 1e-16, 1e16, 1e16);
  AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
  AliMixture(19, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);
  
  // **************** 
  //     Defines tracking media parameters. 
  //     Les valeurs sont commentees pour laisser le defaut 
  //     a GEANT (version 3-21, page CONS200), f.m. 
  epsil  = .001;  // Tracking precision, 
  stemax = -1.;   // Maximum displacement for multiple scat 
  tmaxfd = -20.;  // Maximum angle due to field deflection 
  deemax = -.3;   // Maximum fractional energy loss, DLS 
  stmin  = -.8;
  
  //    Air 
  
  AliMedium(2015, "AIR_L3_US", 15, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Beryllium 
  
  AliMedium(2005, "BE_L3_US", 5, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);

  
    //    Aluminium 
  
  AliMedium(2004, "AL_L3_US", 4, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);

  //   Vacuum

  AliMedium(2016, "VA_L3_US", 16, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  //    Steel 
  
  AliMedium(2019, "ST_L3_US", 19, 0, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
}
//
void AliPIPEv3::Undulation(char *undul, Float_t pitch, Float_t thick,
                        Float_t zundul, Float_t rundul, char (*cone)[5])
{
  //
  // RUNDUL   : Internal radius of the undulated chamber
  // THICK    : material thickness
  // PITCH    : one-QUARTER wave of undulation (cm)
  // ZUNDUL   : half length (cm)
  //
  // The undulated structure is desgned as a superposition of eight CONES
  // of suitable sizes, where the inner/outer radius of the cone increases,
  // then decreases, each half of the wave is assumed to be a semicircle,
  // which allows to calculate the thickness and the radii of the cone, by
  // dividing the semicircle into 4 parts of equal arc length.
  // Thus apear the constants 0.293 and 0.707.
  //

  const Float_t const1 = .293;
  const Float_t const2 = .707;

  AliMC* pMC = AliMC::GetMC();

  // Local variables
  Int_t j, nwave;
  Float_t dcone1[5], dcone2[5], dcone3[5], dcone4[5], dcone5[5],
    dcone6[5], dcone7[5], dcone8[5];
  Float_t xc, yc, zc, dundul[3];
  Int_t *idtmed = gAlice->Idtmed();

  // Function Body

  dcone1[0] = const1 * pitch / 2;
  dcone1[1] = rundul;
  dcone1[2] = dcone1[1] + thick;
  dcone1[3] = dcone1[1] + const2 * pitch;
  dcone1[4] = dcone1[3] + thick;

  dcone2[0] = const2 * pitch / 2;
  dcone2[1] = dcone1[3];
  dcone2[2] = dcone1[4];
  dcone2[3] = dcone2[1] + const1 * pitch;
  dcone2[4] = dcone2[3] + thick;

  dcone3[0] = dcone2[0];
  dcone3[1] = dcone2[3];
  dcone3[2] = dcone2[4];
  dcone3[3] = dcone2[1];
  dcone3[4] = dcone2[2];

  dcone4[0] = dcone1[0];
  dcone4[1] = dcone1[3];
  dcone4[2] = dcone1[4];
  dcone4[3] = dcone1[1];
  dcone4[4] = dcone1[2];

  dcone5[0] = dcone1[0];
  dcone5[1] = dcone1[1] - thick;
  dcone5[2] = dcone1[1];
  dcone5[3] = dcone5[1] - const2 * pitch;
  dcone5[4] = dcone5[3] + thick;

  dcone6[0] = dcone2[0];
  dcone6[1] = dcone5[3];
  dcone6[2] = dcone5[4];
  dcone6[3] = dcone6[1] - const1 * pitch;
  dcone6[4] = dcone6[3] + thick;
  dcone7[0] = dcone6[0];
  dcone7[1] = dcone6[3];
  dcone7[2] = dcone6[4];
  dcone7[3] = dcone5[3];
  dcone7[4] = dcone5[4];

  dcone8[0] = dcone5[0];
  dcone8[1] = dcone7[3];
  dcone8[2] = dcone7[4];
  dcone8[3] = dcone5[1];
  dcone8[4] = dcone5[2];

  pMC->Gsvolu(cone[0], "CONE", idtmed[2015], dcone1, 5);
  pMC->Gsvolu(cone[1], "CONE", idtmed[2015], dcone2, 5);
  pMC->Gsvolu(cone[2], "CONE", idtmed[2015], dcone3, 5);
  pMC->Gsvolu(cone[3], "CONE", idtmed[2015], dcone4, 5);
  pMC->Gsvolu(cone[4], "CONE", idtmed[2015], dcone5, 5);
  pMC->Gsvolu(cone[5], "CONE", idtmed[2015], dcone6, 5);
  pMC->Gsvolu(cone[6], "CONE", idtmed[2015], dcone7, 5);
  pMC->Gsvolu(cone[7], "CONE", idtmed[2015], dcone8, 5);
  pMC->Gsatt(cone[0], "SEEN", 0);
  pMC->Gsatt(cone[1], "SEEN", 0);
  pMC->Gsatt(cone[2], "SEEN", 0);
  pMC->Gsatt(cone[3], "SEEN", 0);
  pMC->Gsatt(cone[4], "SEEN", 0);
  pMC->Gsatt(cone[5], "SEEN", 0);
  pMC->Gsatt(cone[6], "SEEN", 0);
  pMC->Gsatt(cone[7], "SEEN", 0);

  // DEFINE AN IMAGINARY TUBE VOLUME FOR UNDULATED CHAMBER, FILL WITH VACUUM

  nwave = Int_t (zundul / (pitch * 2) + .1);
  dundul[2] = pitch * 2 * nwave;
  dundul[1] = rundul + pitch + thick * 2;
  //
  dundul[0] = 1e-4;
  pMC->Gsvolu(undul, "TUBE", idtmed[2015], dundul, 3);

  xc = 0;
  yc = 0;
  zc = -dundul[2] + dcone1[0];
  for (j = 1; j <= nwave; ++j) {
    pMC->Gspos(cone[0], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone1[0] + dcone2[0];
    pMC->Gspos(cone[1], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone2[0] + dcone3[0];
    pMC->Gspos(cone[2], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone3[0] + dcone4[0];
    pMC->Gspos(cone[3], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone4[0] + dcone5[0];
    pMC->Gspos(cone[4], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone5[0] + dcone6[0];
    pMC->Gspos(cone[5], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone6[0] + dcone7[0];
    pMC->Gspos(cone[6], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone7[0] + dcone8[0];
    pMC->Gspos(cone[7], j, undul, xc, yc, zc, 0, "ONLY");
    zc = zc + dcone8[0] + dcone1[0];
  }
}
