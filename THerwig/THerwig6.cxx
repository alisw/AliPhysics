// definition of c++ Class THerwig6 to be used in ROOT
// this is a c++ interface to the F77 Herwig6 program
// author: j. g. contreras jgcn@moni.mda.cinvestav.mx
// date: december 22, 2000
//  Class THerwig6 is an interface to the Herwig program
//
// C-----------------------------------------------------------------------
// C                           H E R W I G
// C
// C            a Monte Carlo event generator for simulating
// C        +---------------------------------------------------+
// C        | Hadron Emission Reactions With Interfering Gluons |
// C        +---------------------------------------------------+
// C I.G. Knowles(*), G. Marchesini(+), M.H. Seymour($) and B.R. Webber(#)
// C-----------------------------------------------------------------------
// C with Minimal Supersymmetric Standard Model Matrix Elements by
// C                  S. Moretti($) and K. Odagiri($)
// C-----------------------------------------------------------------------
// C R parity violating Supersymmetric Decays and Matrix Elements by
// C                          P. Richardson(&)
// C-----------------------------------------------------------------------
// C matrix element corrections to top decay and Drell-Yan type processes
// C                         by G. Corcella(+)
// C-----------------------------------------------------------------------
// C Deep Inelastic Scattering and Heavy Flavour Electroproduction by
// C                  G. Abbiendi(@) and L. Stanco(%)
// C-----------------------------------------------------------------------
// C and Jet Photoproduction in Lepton-Hadron Collisions by J. Chyla(~)
// C-----------------------------------------------------------------------
// C(*)  Department of Physics & Astronomy, University of Edinburgh
// C(+)  Dipartimento di Fisica, Universita di Milano
// C($)  Rutherford Appleton Laboratory
// C(#)  Cavendish Laboratory, Cambridge
// C(&)  Department of Physics, University of Oxford
// C(@)  Dipartimento di Fisica, Universita di Bologna
// C(%)  Dipartimento di Fisica, Universita di Padova
// C(~)  Institute of Physics, Prague
// C-----------------------------------------------------------------------
// C                  Version 6.100 - 16th December 1999
// C-----------------------------------------------------------------------
// C Main reference:
// C    G.Marchesini,  B.R.Webber,  G.Abbiendi,  I.G.Knowles,  M.H.Seymour,
// C    and L.Stanco, Computer Physics Communications 67 (1992) 465.
// C-----------------------------------------------------------------------
// C Please send e-mail about  this program  to one of the  authors at the
// C following Internet addresses:
// C    I.Knowles@ed.ac.uk        Giuseppe.Marchesini@mi.infn.it
// C    M.Seymour@rl.ac.uk        webber@hep.phy.cam.ac.uk
// C-----------------------------------------------------------------------


#include "THerwig6.h"
#include "HCommon.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TObjArray.h"


ClassImp(THerwig6)

extern "C" {
  void   herwig6_open_fortran_file_ (int* lun, char* name, int);
  void   herwig6_close_fortran_file_(int* lun);
}


THerwig6::THerwig6() : TGenerator("Herwig6","Herwig6")
{

// THerwig6 constructor: creates a TClonesArray in which it will store all
// particles. Note that there may be only one functional THerwig6 object
// at a time, so it's not use to create more than one instance of it.

  delete fParticles; // was allocated as TObjArray in TGenerator
  fParticles = new TClonesArray("TParticle",50);

  // initialize common-blocks
 }

THerwig6::THerwig6(const THerwig6 & source): TGenerator(source)
{
    Fatal("THerwig6","Copy constructor not implemented yet");
}
//------------------------------------------------------------------------------
 THerwig6::~THerwig6()
 {
   // Destructor. The data members of TGenerator are delete by itself
 }

//______________________________________________________________________________
void THerwig6::GenerateEvent()
{

  //initialize event
  hwuine_();
  // generate hard subprocess
  hwepro_();
  // generate parton cascades
  hwbgen_();
  // do heavy objects decay
  hwdhob_();
  // do cluster formation
  hwcfor_();
  // do cluster decays
  hwcdec_();
  // do unstable particle decays
  hwdhad_();
  // do heavy flavor hadrons decay
  hwdhvy_();
  // add soft underlying event
  hwmevt_();
  // finish event
  hwufne_();
}
//______________________________________________________________________________
void THerwig6::OpenFortranFile(int lun, char* name) {
  herwig6_open_fortran_file_(&lun, name, strlen(name));
}

//______________________________________________________________________________
void THerwig6::CloseFortranFile(int lun) {
  herwig6_close_fortran_file_(&lun);
}

void THerwig6::Initialize(const char *beam, const char *target, double pbeam1, double pbeam2, int iproc)

{
  // perform the initialization for Herwig6
  // sets correct title.
  // after calling this method all parameters are set to their default
  // values. If you want to modify any parameter you have to set the new
  // value after calling Initialize and before PrepareRun.

   char  cbeam[8];
    strncpy(cbeam,beam,8);
   char  ctarget[8];
   strncpy(ctarget,target,8);
   printf("\n Initializing Herwig !! \n");
   if ( (!strncmp(beam, "E+"    ,2)) &&
        (!strncmp(beam, "E-"    ,2)) &&
        (!strncmp(beam, "MU+"   ,3)) &&
        (!strncmp(beam, "MU-"   ,3)) &&
        (!strncmp(beam, "NUE"   ,3)) &&
        (!strncmp(beam, "NUEB"  ,4)) &&
        (!strncmp(beam, "NUMU"  ,4)) &&
        (!strncmp(beam, "NMUB"  ,4)) &&
        (!strncmp(beam, "NTAU"  ,4)) &&
        (!strncmp(beam, "NTAB"  ,4)) &&
        (!strncmp(beam, "GAMA"  ,4)) &&
        (!strncmp(beam, "P       ",8)) &&
        (!strncmp(beam, "PBAR    ",8)) &&
        (!strncmp(beam, "N"     ,1)) &&
        (!strncmp(beam, "NBAR"  ,4)) &&
        (!strncmp(beam, "PI+"   ,3)) &&
        (!strncmp(beam, "PI-"   ,3)) ) {
      printf("WARNING! In THerwig6:Initialize():\n");
      printf(" specified beam=%s is unrecognized .\n",beam);
      printf(" resetting to \"P\" .");
      sprintf(cbeam,"P");
   }

   if ( (!strncmp(target, "E+"    ,2)) &&
        (!strncmp(target, "E-"    ,2)) &&
        (!strncmp(target, "MU+"   ,3)) &&
        (!strncmp(target, "MU-"   ,3)) &&
        (!strncmp(target, "NUE"   ,3)) &&
        (!strncmp(target, "NUEB"  ,4)) &&
        (!strncmp(target, "NUMU"  ,4)) &&
        (!strncmp(target, "NMUB"  ,4)) &&
        (!strncmp(target, "NTAU"  ,4)) &&
        (!strncmp(target, "NTAB"  ,4)) &&
        (!strncmp(target, "GAMA"  ,4)) &&
        (!strncmp(target, "P       ",8)) &&
        (!strncmp(target, "PBAR    ",8)) &&
        (!strncmp(target, "N"     ,1)) &&
        (!strncmp(target, "NBAR"  ,4)) &&
        (!strncmp(target, "PI+"   ,3)) &&
        (!strncmp(target, "PI-"   ,3)) ) {
      printf("WARNING! In THerwig6:Initialize():\n");
      printf(" specified target=%s is unrecognized .\n",target);
      printf(" resetting to \"P\" .");
      sprintf(ctarget,"P");
   }

   // initialization:
   // type of beams
   strncpy(HWBMCH.PART1,beam,8);
   strncpy(HWBMCH.PART2,target,8);
   // momentum of beams
   HWPROC.PBEAM1=pbeam1;
   HWPROC.PBEAM2=pbeam2;
   // process to generate
   HWPROC.IPROC=iproc;
   // not used in the class definition
   HWPROC.MAXEV=1;
   
   // reset all parameters
   hwigin_();

   // set correct title
   //char atitle[132];
   double win=pbeam1+pbeam2;
   printf("\n %s - %s at %g GeV \n",beam,target,win);
   //sprintf(atitle,"%s-%s at %g GeV",cbeam,ctarget,win);
   //SetTitle(atitle);
}

void THerwig6::InitializeJimmy(const char *beam, const char *target, double pbeam1, double pbeam2, int iproc)

{
  // perform the initialization for Herwig6
  // sets correct title.
  // after calling this method all parameters are set to their default
  // values. If you want to modify any parameter you have to set the new
  // value after calling Initialize and before PrepareRun.

   char  cbeam[8];
   strncpy(cbeam,beam,8);
   char  ctarget[8];
   strncpy(ctarget,target,8);
   printf("\n Initializing Herwig !! \n");
   if ( (!strncmp(beam, "E+"    ,2)) &&
        (!strncmp(beam, "E-"    ,2)) &&
        (!strncmp(beam, "MU+"   ,3)) &&
        (!strncmp(beam, "MU-"   ,3)) &&
        (!strncmp(beam, "NUE"   ,3)) &&
        (!strncmp(beam, "NUEB"  ,4)) &&
        (!strncmp(beam, "NUMU"  ,4)) &&
        (!strncmp(beam, "NMUB"  ,4)) &&
        (!strncmp(beam, "NTAU"  ,4)) &&
        (!strncmp(beam, "NTAB"  ,4)) &&
        (!strncmp(beam, "GAMA"  ,4)) &&
        (!strncmp(beam, "P       ",8)) &&
        (!strncmp(beam, "PBAR    ",8)) &&
        (!strncmp(beam, "N"     ,1)) &&
        (!strncmp(beam, "NBAR"  ,4)) &&
        (!strncmp(beam, "PI+"   ,3)) &&
        (!strncmp(beam, "PI-"   ,3)) ) {
      printf("WARNING! In THerwig6:Initialize():\n");
      printf(" specified beam=%s is unrecognized .\n",beam);
      printf(" resetting to \"P\" .");
      sprintf(cbeam,"P");
   }

   if ( (!strncmp(target, "E+"    ,2)) &&
        (!strncmp(target, "E-"    ,2)) &&
        (!strncmp(target, "MU+"   ,3)) &&
        (!strncmp(target, "MU-"   ,3)) &&
        (!strncmp(target, "NUE"   ,3)) &&
        (!strncmp(target, "NUEB"  ,4)) &&
        (!strncmp(target, "NUMU"  ,4)) &&
        (!strncmp(target, "NMUB"  ,4)) &&
        (!strncmp(target, "NTAU"  ,4)) &&
        (!strncmp(target, "NTAB"  ,4)) &&
        (!strncmp(target, "GAMA"  ,4)) &&
        (!strncmp(target, "P       ",8)) &&
        (!strncmp(target, "PBAR    ",8)) &&
        (!strncmp(target, "N"     ,1)) &&
        (!strncmp(target, "NBAR"  ,4)) &&
        (!strncmp(target, "PI+"   ,3)) &&
        (!strncmp(target, "PI-"   ,3)) ) {
      printf("WARNING! In THerwig6:Initialize():\n");
      printf(" specified target=%s is unrecognized .\n",target);
      printf(" resetting to \"P\" .");
      sprintf(ctarget,"P");
   }

   // initialization:
   // type of beams
   strncpy(HWBMCH.PART1,beam,8);
   strncpy(HWBMCH.PART2,target,8);
   // momentum of beams
   HWPROC.PBEAM1=pbeam1;
   HWPROC.PBEAM2=pbeam2;
   // process to generate
   HWPROC.IPROC=iproc;
   // not used in the class definition
   HWPROC.MAXEV=1;

   // reset all parameters
   hwigin_();
   // JIMMY initialization
   jimmin_();

   // set correct title
//   char atitle[132];
   double win=pbeam1+pbeam2;
   printf("\n %s - %s at %g GeV",beam,target,win);
//   sprintf(atitle,"%s-%s at %g GeV",cbeam,ctarget,win);
//   SetTitle(atitle);
}

void THerwig6::PrepareRun()
{
  // compute parameter dependent constants
  hwuinc_();
  // initialize elementary processes
  hweini_();
}

void THerwig6::PrepareRunJimmy()
{
  // compute parameter dependent constants
  hwuinc_();
  // initialize elementary processes
  hweini_();
  // more initializations for JIMMY
  jminit_();
}
//______________________________________________________________________________
TObjArray* THerwig6::ImportParticles(Option_t *option)
{
//
//  Default primary creation method. It reads the /HEPEVT/ common block which
//  has been filled by the GenerateEvent method. If the event generator does
//  not use the HEPEVT common block, This routine has to be overloaded by
//  the subclasses.
//  The default action is to store only the stable particles (ISTHEP = 1)
//  This can be demanded explicitly by setting the option = "Final"
//  If the option = "All", all the particles are stored.
//
  fParticles->Clear();
  Int_t numpart = HEPEVT.NHEP;
  TClonesArray &a = *((TClonesArray*)fParticles);
  if (!strcmp(option,"") || !strcmp(option,"Final")) {
    for (Int_t i = 0; i < numpart; i++) {
      if (HEPEVT.ISTHEP[i] == 1) {
//
//  Use the common block values for the TParticle constructor
//
        new(a[i]) TParticle(
                                   HEPEVT.IDHEP[i],
                                   HEPEVT.ISTHEP[i],
                                   HEPEVT.JMOHEP[i][0]-1,
                                   HEPEVT.JMOHEP[i][1]-1,
                                   HEPEVT.JDAHEP[i][0]-1,
                                   HEPEVT.JDAHEP[i][1]-1,

                                   HEPEVT.PHEP[i][0],
                                   HEPEVT.PHEP[i][1],
                                   HEPEVT.PHEP[i][2],
                                   HEPEVT.PHEP[i][3],
                                   HEPEVT.VHEP[i][0],
                                   HEPEVT.VHEP[i][1],
                                   HEPEVT.VHEP[i][2],
                                   HEPEVT.VHEP[i][3]);
        }
     }
  }
  else if (!strcmp(option,"All")) {
    for (Int_t i = 0; i < numpart; i++) {
      new(a[i]) TParticle(
                                   HEPEVT.IDHEP[i],
                                   HEPEVT.ISTHEP[i],
                                   HEPEVT.JMOHEP[i][0]-1,
                                   HEPEVT.JMOHEP[i][1]-1,
                                   HEPEVT.JDAHEP[i][0]-1,
                                   HEPEVT.JDAHEP[i][1]-1,

                                   HEPEVT.PHEP[i][0],
                                   HEPEVT.PHEP[i][1],
                                   HEPEVT.PHEP[i][2],
                                   HEPEVT.PHEP[i][3],
                                   HEPEVT.VHEP[i][0],
                                   HEPEVT.VHEP[i][1],
                                   HEPEVT.VHEP[i][2],
                                   HEPEVT.VHEP[i][3]);
    }
  }
  return fParticles;
}

//______________________________________________________________________________
Int_t THerwig6::ImportParticles(TClonesArray *particles, Option_t *option)
{
//
//  Default primary creation method. It reads the /HEPEVT/ common block which
//  has been filled by the GenerateEvent method. If the event generator does
//  not use the HEPEVT common block, This routine has to be overloaded by
//  the subclasses.
//  The function loops on the generated particles and store them in
//  the TClonesArray pointed by the argument particles.
//  The default action is to store only the stable particles (ISTHEP = 1)
//  This can be demanded explicitly by setting the option = "Final"
//  If the option = "All", all the particles are stored.
//
  if (particles == 0) return 0;
  TClonesArray &refParticles = *particles;
  refParticles.Clear();
  Int_t numpart = HEPEVT.NHEP;
  if (!strcmp(option,"") || !strcmp(option,"Final")) {
    for (Int_t i = 0; i < numpart; i++) {
      if (HEPEVT.ISTHEP[i] == 1) {
//
//  Use the common block values for the TParticle constructor
//
        new(refParticles[i]) TParticle(
                                   HEPEVT.IDHEP[i],
                                   HEPEVT.ISTHEP[i],
                                   HEPEVT.JMOHEP[i][0]-1,
                                   HEPEVT.JMOHEP[i][1]-1,
                                   HEPEVT.JDAHEP[i][0]-1,
                                   HEPEVT.JDAHEP[i][1]-1,

                                   HEPEVT.PHEP[i][0],
                                   HEPEVT.PHEP[i][1],
                                   HEPEVT.PHEP[i][2],
                                   HEPEVT.PHEP[i][3],
                                   HEPEVT.VHEP[i][0],
                                   HEPEVT.VHEP[i][1],
                                   HEPEVT.VHEP[i][2],
                                   HEPEVT.VHEP[i][3]);
        }
     }
  }
  else if (!strcmp(option,"All")) {
    for (Int_t i = 0; i< numpart; i++) {
      new(refParticles[i]) TParticle(
                                   HEPEVT.IDHEP[i],
                                   HEPEVT.ISTHEP[i],
                                   HEPEVT.JMOHEP[i][0]-1,
                                   HEPEVT.JMOHEP[i][1]-1,
                                   HEPEVT.JDAHEP[i][0]-1,
                                   HEPEVT.JDAHEP[i][1]-1,

                                   HEPEVT.PHEP[i][0],
                                   HEPEVT.PHEP[i][1],
                                   HEPEVT.PHEP[i][2],
                                   HEPEVT.PHEP[i][3],
                                   HEPEVT.VHEP[i][0],
                                   HEPEVT.VHEP[i][1],
                                   HEPEVT.VHEP[i][2],
                                   HEPEVT.VHEP[i][3]); // 
    }
  }
  return numpart;
}

void THerwig6::Hwigin()
{
  hwigin_();
}

void THerwig6::Hwuinc()
{
  hwuinc_();
}

void THerwig6::Hwusta(const char* name)

{
  hwusta_(name,8);
}

void THerwig6::Hweini()

{
  hweini_();
}

void THerwig6::Hwuine()

{
  hwuine_();
}

void THerwig6::Hwepro()

{
  hwepro_();
}

void THerwig6::Hwbgen()

{
  hwbgen_();
}

void THerwig6::Hwdhob()

{
  hwdhob_();
}

void THerwig6::Hwcfor()

{
  hwcfor_();
}

void THerwig6::Hwcdec()

{
  hwcdec_();
}

void THerwig6::Hwdhad()

{
  hwdhad_();
}

void THerwig6::Hwdhvy()

{
  hwdhvy_();
}

void THerwig6::Hwmevt()

{
  hwmevt_();
}

void THerwig6::Hwufne()

{
  hwufne_();
}

void THerwig6::Hwefin()

{
  hwefin_();
}

void THerwig6::Hwiodk(int iopt)

{
  hwiodk_(iopt);
}

void THerwig6::SetupTest()
{
  // exampe of running herwig and generating one event
  // after changing some options
  Initialize("P","PBAR",900.,900.,1500);
  // here you can set some parameters
  SetPTMIN(15.); // Min pt in hadronic jet production
  SetYJMIN(-4.); // Min jet rapidity
  SetYJMAX(4.);  // Max jet rapidity
  // after you set your wished parameters
  // herwig can do its work
  PrepareRun();
  int nEvToGenerate=1;
  for (int i=0;i<nEvToGenerate;i++)
    {
      GenerateEvent();
      // do your stuff. For ex:
      int nOfPar=GetNumberOfParticles(); // from TGenerator
      for (int j=0; j<nOfPar; j++)
	{
	  TParticle* p=GetParticle(j);
	  // here you do whatever you want with the particle
	  p->Print();
	};
    };
}

// Jimmy subroutines

void THerwig6::Jminit()
{
  jminit_();
}

void THerwig6::Jimmin()
{
  jimmin_();
}

void THerwig6::Jmefin()
{
  jmefin_();
}

void THerwig6::PrintEvt()
{
    hwuepr_();
    
}

  // acces to hep common block
int         THerwig6::GetNEVHEP        () const     { return HEPEVT.NEVHEP; }
int         THerwig6::GetNhep          () const     { return HEPEVT.NHEP; }
int         THerwig6::GetISTHEP    (int i)const     { return HEPEVT.ISTHEP[i-1]; }
int         THerwig6::GetIDHEP     (int i)const     { return HEPEVT.IDHEP[i-1]; }
int         THerwig6::GetJMOHEP (int i, int j) const
{ return HEPEVT.JMOHEP[i-1][j-1]; }
int         THerwig6::GetJDAHEP (int i, int j) const
{ return HEPEVT.JDAHEP[i-1][j-1]; }
double      THerwig6::GetPHEP   (int i, int j) const
{ return HEPEVT.PHEP[i-1][j-1]; }
double      THerwig6::GetVHEP   (int i, int j) const
{ return HEPEVT.VHEP[i-1][j-1]; }

// access to Herwig6 common-blocks
// WARNING: Some arrays start in 1, others in 0. Look up the manual!

// /HWBEAM/

int         THerwig6::GetIPART1        () const     { return HWBEAM.IPART1; }
int         THerwig6::GetIPART2        () const     { return HWBEAM.IPART2; }

// /HWBMCH/
char*       THerwig6::GetPART1         () const     { return HWBMCH.PART1; }
char*       THerwig6::GetPART2         () const     { return HWBMCH.PART2; }


// /HWPROC/
double      THerwig6::GetEBEAM1        () const     { return HWPROC.EBEAM1; }
double      THerwig6::GetEBEAM2        () const     { return HWPROC.EBEAM2; }
double      THerwig6::GetPBEAM1        () const     { return HWPROC.PBEAM1; }
double      THerwig6::GetPBEAM2        () const     { return HWPROC.PBEAM2; }
int         THerwig6::GetIPROC         () const     { return HWPROC.IPROC; }
int         THerwig6::GetMAXEV         () const     { return HWPROC.MAXEV; }

// /HWPRAM/
double      THerwig6::GetQCDLAM        () const     { return HWPRAM.QCDLAM; }
void        THerwig6::SetQCDLAM   (double q)        { HWPRAM.QCDLAM = q; }
double      THerwig6::GetVQCUT         () const     { return HWPRAM.VQCUT; }
void        THerwig6::SetVQCUT    (double v)        { HWPRAM.VQCUT = v; }
double      THerwig6::GetVGCUT         () const     { return HWPRAM.VGCUT; }
void        THerwig6::SetVGCUT    (double v)        { HWPRAM.VGCUT = v; }
double      THerwig6::GetVPCUT         () const     { return HWPRAM.VPCUT; }
void        THerwig6::SetVPCUT    (double v)        { HWPRAM.VPCUT = v; }
double      THerwig6::GetCLMAX         () const     { return HWPRAM.CLMAX; }
void        THerwig6::SetCLMAX    (double c)        { HWPRAM.CLMAX = c; }
double      THerwig6::GetCLPOW         () const     { return HWPRAM.CLPOW; }
void        THerwig6::SetCLPOW    (double c)        { HWPRAM.CLPOW = c; }
double      THerwig6::GetPSPLT    (int i) const     { return HWPRAM.PSPLT[i-1];}
void        THerwig6::SetPSPLT    (int i, double p) { HWPRAM.PSPLT[i-1] = p;}
double      THerwig6::GetQDIQK         () const     { return HWPRAM.QDIQK; }
void        THerwig6::SetQDIQK    (double q)        { HWPRAM.QDIQK = q; }
double      THerwig6::GetPDIQK         () const     { return HWPRAM.PDIQK; }
void        THerwig6::SetPDIQK    (double p)        { HWPRAM.PDIQK = p; }
double      THerwig6::GetQSPAC         () const     { return HWPRAM.QSPAC; }
void        THerwig6::SetQSPAC    (double q)        { HWPRAM.QSPAC = q; }
double      THerwig6::GetPTRMS         () const     { return HWPRAM.PTRMS; }
void        THerwig6::SetPTRMS    (double p)        { HWPRAM.PTRMS = p; }
double      THerwig6::GetENSOF         () const     { return HWPRAM.ENSOF; }
void        THerwig6::SetENSOF    (double e)        { HWPRAM.ENSOF = e; }
int         THerwig6::GetIPRINT        () const     { return HWPRAM.IPRINT; }
void        THerwig6::SetIPRINT   (int i)           { HWPRAM.IPRINT = i; }
int         THerwig6::GetMODPDF   (int i) const     { return HWPRAM.MODPDF[i-1];}
void        THerwig6::SetMODPDF   (int i, int j)  { HWPRAM.MODPDF[i-1] = j; }
int         THerwig6::GetNSTRU         () const     { return HWPRAM.NSTRU; }
void        THerwig6::SetNSTRU    (int i)          { HWPRAM.NSTRU = i; }

// /HWPRCH/
char*       THerwig6::GetAUTPDF     (int i)         { return HWPRCH.AUTPDF[i-1]; }
void        THerwig6::SetAUTPDF(int i,const char* s){ strncpy(HWPRCH.AUTPDF[i-1],s,20);}
char*       THerwig6::GetBDECAY        ()           { return HWPRCH.BDECAY; }

// /HWEVNT/
double      THerwig6::GetAVWGT         () const     { return HWEVNT.AVWGT; }
int         THerwig6::GetMAXPR         () const     { return HWEVNT.MAXPR; }
void        THerwig6::SetMAXPR    (int i)           { HWEVNT.MAXPR = i; }
int         THerwig6::GetMAXER         () const     { return HWEVNT.MAXER; }
void        THerwig6::SetMAXER    (int i)           { HWEVNT.MAXER = i; }
int         THerwig6::GetNRN      (int i) const     { return HWEVNT.NRN[i-1]; }
void        THerwig6::SetNRN    (int i, int j)      { HWEVNT.NRN[i-1] = j; }
double      THerwig6::GetEVWGT         () const     { return HWEVNT.EVWGT; }

int         THerwig6::GetIDHW     (int i) const     { return HWEVNT.IDHW[i]; }

int         THerwig6::GetIERROR        () const     { return HWEVNT.IERROR; }

// /HWHARD/
double      THerwig6::GetPTMIN         () const     { return HWHARD.PTMIN; }
void        THerwig6::SetPTMIN    (double d)        { HWHARD.PTMIN = d; }
double      THerwig6::GetPTMAX         () const     { return HWHARD.PTMAX; }
void        THerwig6::SetPTMAX    (double d)        { HWHARD.PTMAX = d; }
double      THerwig6::GetPTPOW         () const     { return HWHARD.PTPOW; }
void        THerwig6::SetPTPOW    (double d)        { HWHARD.PTPOW = d; }
double      THerwig6::GetYJMIN         () const     { return HWHARD.YJMIN; }
void        THerwig6::SetYJMIN    (double d)        { HWHARD.YJMIN = d; }
double      THerwig6::GetYJMAX         () const     { return HWHARD.YJMAX; }
void        THerwig6::SetYJMAX    (double d)        { HWHARD.YJMAX = d; }
double      THerwig6::GetQ2MIN         () const     { return HWHARD.Q2MIN; }
void        THerwig6::SetQ2MIN    (double d)        { HWHARD.Q2MIN = d; }
double      THerwig6::GetQ2MAX         () const     { return HWHARD.Q2MAX; }
void        THerwig6::SetQ2MAX    (double d)        { HWHARD.Q2MAX = d; }
double      THerwig6::GetYBMIN         () const     { return HWHARD.YBMIN; }
void        THerwig6::SetYBMIN    (double d)        { HWHARD.YBMIN = d; }
double      THerwig6::GetYBMAX         () const     { return HWHARD.YBMAX; }
void        THerwig6::SetYBMAX    (double d)        { HWHARD.YBMAX = d; }
double      THerwig6::GetZJMAX        ()  const     { return HWHARD.ZJMAX; }
void        THerwig6::SetZJMAX    (double d)        { HWHARD.ZJMAX = d; }

// /HWPROP/
double      THerwig6::GetRMASS      (int i) const   { return HWPROP.RMASS[i]; }
void        THerwig6::SetRMASS    (int i, double r) { HWPROP.RMASS[i] = r; }


void        THerwig6::GetRNAME (int i, char a[9])   { for (int j=0;j<8;j++) a[j] = HWUNAM.RNAME[i][j]; a[8] = '\0';}
