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
#include "TClonesArray.h"
#include "TParticle.h"
#include "TObjArray.h"


ClassImp(THerwig6)

extern "C" {
  void*  herwig6_common_block_address_(char*, int len);
  void   herwig6_open_fortran_file_ (int* lun, char* name, int);
  void   herwig6_close_fortran_file_(int* lun);
}


THerwig6::THerwig6() : TGenerator("Herwig6","Herwig6"),
  fHepevt((Hepevt_t*) herwig6_common_block_address_((char*)"HEPEVT",6)),
  fHwbeam((Hwbeam_t*) herwig6_common_block_address_((char*)"HWBEAM",6)),
  fHwbmch((Hwbmch_t*) herwig6_common_block_address_((char*)"HWBMCH",6)),
  fHwproc((Hwproc_t*) herwig6_common_block_address_((char*)"HWPROC",6)),
  fHwpram((Hwpram_t*) herwig6_common_block_address_((char*)"HWPRAM",6)),
  fHwprch((Hwprch_t*) herwig6_common_block_address_((char*)"HWPRCH",6)),
  fHwpart((Hwpart_t*) herwig6_common_block_address_((char*)"HWPART",6)),
  fHwparp((Hwparp_t*) herwig6_common_block_address_((char*)"HWPARP",6)),
  fHwbosc((Hwbosc_t*) herwig6_common_block_address_((char*)"HWBOSC",6)),
  fHwparc((Hwparc_t*) herwig6_common_block_address_((char*)"HWPARC",6)),
  fHwbrch((Hwbrch_t*) herwig6_common_block_address_((char*)"HWBRCH",6)),
  fHwevnt((Hwevnt_t*) herwig6_common_block_address_((char*)"HWEVNT",6)),
  fHwhard((Hwhard_t*) herwig6_common_block_address_((char*)"HWHARD",6)),
  fHwprop((Hwprop_t*) herwig6_common_block_address_((char*)"HWPROP",6)),
  fHwunam((Hwunam_t*) herwig6_common_block_address_((char*)"HWUNAM",6)),
  fHwupdt((Hwupdt_t*) herwig6_common_block_address_((char*)"HWUPDT",6)),
  fHwuwts((Hwuwts_t*) herwig6_common_block_address_((char*)"HWUWTS",6)),
  fHwuclu((Hwuclu_t*) herwig6_common_block_address_((char*)"HWUCLU",6)),
  fHwdist((Hwdist_t*) herwig6_common_block_address_((char*)"HWDIST",6)),
  fHwqdks((Hwqdks_t*) herwig6_common_block_address_((char*)"HWQDKS",6)),
  fHwusud((Hwusud_t*) herwig6_common_block_address_((char*)"HWUSUD",6)),
  fHwsusy((Hwsusy_t*) herwig6_common_block_address_((char*)"HWSUSY",6)),
  fHwrpar((Hwrpar_t*) herwig6_common_block_address_((char*)"HWRPAR",6)),
  fHwminb((Hwminb_t*) herwig6_common_block_address_((char*)"HWMINB",6)),
  fHwclus((Hwclus_t*) herwig6_common_block_address_((char*)"HWCLUS",6)),
  fHwgrav((Hwgrav_t*) herwig6_common_block_address_((char*)"HWGRAV",6)),
  fHw6202((Hw6202_t*) herwig6_common_block_address_((char*)"HW6202",6)),
  fHw6203((Hw6203_t*) herwig6_common_block_address_((char*)"HW6203",6)),
  fHw6300((Hw6300_t*) herwig6_common_block_address_((char*)"HW6300",6)),
  fHwpmrs((Hwpmrs_t*) herwig6_common_block_address_((char*)"HWPMRS",6)),
  fHwcirc((Hwcirc_t*) herwig6_common_block_address_((char*)"HWCIRC",6)),
  fHwdspb((Hwdspb_t*) herwig6_common_block_address_((char*)"HWDSPB",6)),
  fHwdsp2((Hwdsp2_t*) herwig6_common_block_address_((char*)"HWDSP2",6)),
  fHwdsp3((Hwdsp3_t*) herwig6_common_block_address_((char*)"HWDSP3",6)),
  fHwdsp4((Hwdsp4_t*) herwig6_common_block_address_((char*)"HWDSP4",6)),
  fHwdspn((Hwdspn_t*) herwig6_common_block_address_((char*)"HWDSPN",6)),
  fHwspin((Hwspin_t*) herwig6_common_block_address_((char*)"HWSPIN",6)),
  fHwstau((Hwstau_t*) herwig6_common_block_address_((char*)"HWSTAU",6)),
  fHwgupr((Hwgupr_t*) herwig6_common_block_address_((char*)"HWGUPR",6)),
  fHw6500((Hw6500_t*) herwig6_common_block_address_((char*)"HW6500",6)),
  fHw6504((Hw6504_t*) herwig6_common_block_address_((char*)"HW6504",6)),
  fHw6506((Hw6506_t*) herwig6_common_block_address_((char*)"HW6506",6))
{

// THerwig6 constructor: creates a TClonesArray in which it will store all
// particles. Note that there may be only one functional THerwig6 object
// at a time, so it's not use to create more than one instance of it.

  delete fParticles; // was allocated as TObjArray in TGenerator
  fParticles = new TClonesArray("TParticle",50);

  // initialize common-blocks
 }

THerwig6::THerwig6(const THerwig6 & source): TGenerator(source),
  fHepevt((Hepevt_t*) herwig6_common_block_address_((char*)"HEPEVT",6)),
  fHwbeam((Hwbeam_t*) herwig6_common_block_address_((char*)"HWBEAM",6)),
  fHwbmch((Hwbmch_t*) herwig6_common_block_address_((char*)"HWBMCH",6)),
  fHwproc((Hwproc_t*) herwig6_common_block_address_((char*)"HWPROC",6)),
  fHwpram((Hwpram_t*) herwig6_common_block_address_((char*)"HWPRAM",6)),
  fHwprch((Hwprch_t*) herwig6_common_block_address_((char*)"HWPRCH",6)),
  fHwpart((Hwpart_t*) herwig6_common_block_address_((char*)"HWPART",6)),
  fHwparp((Hwparp_t*) herwig6_common_block_address_((char*)"HWPARP",6)),
  fHwbosc((Hwbosc_t*) herwig6_common_block_address_((char*)"HWBOSC",6)),
  fHwparc((Hwparc_t*) herwig6_common_block_address_((char*)"HWPARC",6)),
  fHwbrch((Hwbrch_t*) herwig6_common_block_address_((char*)"HWBRCH",6)),
  fHwevnt((Hwevnt_t*) herwig6_common_block_address_((char*)"HWEVNT",6)),
  fHwhard((Hwhard_t*) herwig6_common_block_address_((char*)"HWHARD",6)),
  fHwprop((Hwprop_t*) herwig6_common_block_address_((char*)"HWPROP",6)),
  fHwunam((Hwunam_t*) herwig6_common_block_address_((char*)"HWUNAM",6)),
  fHwupdt((Hwupdt_t*) herwig6_common_block_address_((char*)"HWUPDT",6)),
  fHwuwts((Hwuwts_t*) herwig6_common_block_address_((char*)"HWUWTS",6)),
  fHwuclu((Hwuclu_t*) herwig6_common_block_address_((char*)"HWUCLU",6)),
  fHwdist((Hwdist_t*) herwig6_common_block_address_((char*)"HWDIST",6)),
  fHwqdks((Hwqdks_t*) herwig6_common_block_address_((char*)"HWQDKS",6)),
  fHwusud((Hwusud_t*) herwig6_common_block_address_((char*)"HWUSUD",6)),
  fHwsusy((Hwsusy_t*) herwig6_common_block_address_((char*)"HWSUSY",6)),
  fHwrpar((Hwrpar_t*) herwig6_common_block_address_((char*)"HWRPAR",6)),
  fHwminb((Hwminb_t*) herwig6_common_block_address_((char*)"HWMINB",6)),
  fHwclus((Hwclus_t*) herwig6_common_block_address_((char*)"HWCLUS",6)),
  fHwgrav((Hwgrav_t*) herwig6_common_block_address_((char*)"HWGRAV",6)),
  fHw6202((Hw6202_t*) herwig6_common_block_address_((char*)"HW6202",6)),
  fHw6203((Hw6203_t*) herwig6_common_block_address_((char*)"HW6203",6)),
  fHw6300((Hw6300_t*) herwig6_common_block_address_((char*)"HW6300",6)),
  fHwpmrs((Hwpmrs_t*) herwig6_common_block_address_((char*)"HWPMRS",6)),
  fHwcirc((Hwcirc_t*) herwig6_common_block_address_((char*)"HWCIRC",6)),
  fHwdspb((Hwdspb_t*) herwig6_common_block_address_((char*)"HWDSPB",6)),
  fHwdsp2((Hwdsp2_t*) herwig6_common_block_address_((char*)"HWDSP2",6)),
  fHwdsp3((Hwdsp3_t*) herwig6_common_block_address_((char*)"HWDSP3",6)),
  fHwdsp4((Hwdsp4_t*) herwig6_common_block_address_((char*)"HWDSP4",6)),
  fHwdspn((Hwdspn_t*) herwig6_common_block_address_((char*)"HWDSPN",6)),
  fHwspin((Hwspin_t*) herwig6_common_block_address_((char*)"HWSPIN",6)),
  fHwstau((Hwstau_t*) herwig6_common_block_address_((char*)"HWSTAU",6)),
  fHwgupr((Hwgupr_t*) herwig6_common_block_address_((char*)"HWGUPR",6)),
  fHw6500((Hw6500_t*) herwig6_common_block_address_((char*)"HW6500",6)),
  fHw6504((Hw6504_t*) herwig6_common_block_address_((char*)"HW6504",6)),
  fHw6506((Hw6506_t*) herwig6_common_block_address_((char*)"HW6506",6))
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
   strncpy(fHwbmch->PART1,beam,8);
   strncpy(fHwbmch->PART2,target,8);
   // momentum of beams
   fHwproc->PBEAM1=pbeam1;
   fHwproc->PBEAM2=pbeam2;
   // process to generate
   fHwproc->IPROC=iproc;
   // not used in the class definition
   fHwproc->MAXEV=1;

   // reset all parameters
   hwigin_();

   // set correct title
   char atitle[132];
   double win=pbeam1+pbeam2;
   printf("\n %s - %s at %g GeV",beam,target,win);
   sprintf(atitle,"%s-%s at %g GeV",cbeam,ctarget,win);
   SetTitle(atitle);
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
   strncpy(fHwbmch->PART1,beam,8);
   strncpy(fHwbmch->PART2,target,8);
   // momentum of beams
   fHwproc->PBEAM1=pbeam1;
   fHwproc->PBEAM2=pbeam2;
   // process to generate
   fHwproc->IPROC=iproc;
   // not used in the class definition
   fHwproc->MAXEV=1;

   // reset all parameters
   hwigin_();
   // JIMMY initialization
   jimmin_();

   // set correct title
   char atitle[132];
   double win=pbeam1+pbeam2;
   printf("\n %s - %s at %g GeV",beam,target,win);
   sprintf(atitle,"%s-%s at %g GeV",cbeam,ctarget,win);
   SetTitle(atitle);
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
  Int_t numpart = fHepevt->NHEP;
  TClonesArray &a = *((TClonesArray*)fParticles);
  if (!strcmp(option,"") || !strcmp(option,"Final")) {
    for (Int_t i = 0; i<=numpart; i++) {
      if (fHepevt->ISTHEP[i] == 1) {
//
//  Use the common block values for the TParticle constructor
//
        new(a[i]) TParticle(
                                   fHepevt->IDHEP[i],
                                   fHepevt->ISTHEP[i],
                                   fHepevt->JMOHEP[i][0]-1,
                                   fHepevt->JMOHEP[i][1]-1,
                                   fHepevt->JDAHEP[i][0]-1,
                                   fHepevt->JDAHEP[i][1]-1,

                                   fHepevt->PHEP[i][0],
                                   fHepevt->PHEP[i][1],
                                   fHepevt->PHEP[i][2],
                                   fHepevt->PHEP[i][3],
                                   fHepevt->VHEP[i][0],
                                   fHepevt->VHEP[i][1],
                                   fHepevt->VHEP[i][2],
                                   fHepevt->VHEP[i][3]);
        }
     }
  }
  else if (!strcmp(option,"All")) {
    for (Int_t i = 0; i<=numpart; i++) {
      new(a[i]) TParticle(
                                   fHepevt->IDHEP[i],
                                   fHepevt->ISTHEP[i],
                                   fHepevt->JMOHEP[i][0]-1,
                                   fHepevt->JMOHEP[i][1]-1,
                                   fHepevt->JDAHEP[i][0]-1,
                                   fHepevt->JDAHEP[i][1]-1,

                                   fHepevt->PHEP[i][0],
                                   fHepevt->PHEP[i][1],
                                   fHepevt->PHEP[i][2],
                                   fHepevt->PHEP[i][3],
                                   fHepevt->VHEP[i][0],
                                   fHepevt->VHEP[i][1],
                                   fHepevt->VHEP[i][2],
                                   fHepevt->VHEP[i][3]);
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
  Int_t numpart = fHepevt->NHEP;
  if (!strcmp(option,"") || !strcmp(option,"Final")) {
    for (Int_t i = 0; i< numpart; i++) {
      if (fHepevt->ISTHEP[i] == 1) {
//
//  Use the common block values for the TParticle constructor
//
        new(refParticles[i]) TParticle(
                                   fHepevt->IDHEP[i],
                                   fHepevt->ISTHEP[i],
                                   fHepevt->JMOHEP[i][0]-1,
                                   fHepevt->JMOHEP[i][1]-1,
                                   fHepevt->JDAHEP[i][0]-1,
                                   fHepevt->JDAHEP[i][1]-1,

                                   fHepevt->PHEP[i][0],
                                   fHepevt->PHEP[i][1],
                                   fHepevt->PHEP[i][2],
                                   fHepevt->PHEP[i][3],
                                   fHepevt->VHEP[i][0],
                                   fHepevt->VHEP[i][1],
                                   fHepevt->VHEP[i][2],
                                   fHepevt->VHEP[i][3]);
        }
     }
  }
  else if (!strcmp(option,"All")) {
    for (Int_t i = 0; i< numpart; i++) {
      new(refParticles[i]) TParticle(
                                   fHepevt->IDHEP[i],
                                   fHepevt->ISTHEP[i],
                                   fHepevt->JMOHEP[i][0]-1,
                                   fHepevt->JMOHEP[i][1]-1,
                                   fHepevt->JDAHEP[i][0]-1,
                                   fHepevt->JDAHEP[i][1]-1,

                                   fHepevt->PHEP[i][0],
                                   fHepevt->PHEP[i][1],
                                   fHepevt->PHEP[i][2],
                                   fHepevt->PHEP[i][3],
                                   fHepevt->VHEP[i][0],
                                   fHepevt->VHEP[i][1],
                                   fHepevt->VHEP[i][2],
                                   fHepevt->VHEP[i][3]);
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

void THerwig6::Hwusta(char* name)

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






