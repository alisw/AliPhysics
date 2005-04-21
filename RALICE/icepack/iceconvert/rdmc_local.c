
/* local general purpose routines for rdmc */
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "rdmc.h"
#include "rdmc_local.h"

char *rdmc_which_format(mcfile *fp){
  static char returnstring[RDMC_MAXLINE];

  switch (fp->format) {
#ifdef DUMAND_ASCII_F
  case DUMAND_ASCII_F:
    sprintf(returnstring, "dumand-siegmund");
    break;
#endif
#ifdef AMANDA_ASCII_F
  case AMANDA_ASCII_F:
    sprintf(returnstring, "amanda-f2000");
    break;
#endif
#ifdef UWI_ASCII_F
  case UWI_ASCII_F:
    sprintf(returnstring, "uwi-raven");
    break;
#endif
#ifdef BAIKAL_BIN_F
  case BAIKAL_BIN_F:
    sprintf(returnstring, "baikal-binary");
    break;
#endif
  default:
    sprintf(returnstring, "unknown");
    break;
  }

  return returnstring;

} /* rdmc_error() */



/***************************************************************************/
/* id tables */
/***************************************************************************/

const rdmc_idtable_t rdmc_pmt_idtable[] =
  {
    {  0, "q370" },
    { 10, "xp2600" },
    { 20, "r2018" },
    { 21, "r5212" },
    { 30, "emi" },
    {  0, "std" }, /* default */
    { RDMC_NA, NULL } /* end id */
  };

const rdmc_idtable_t rdmc_sphere_idtable[] =
  {
    { 0, "russ" },
    { 1, "naut" },
    { 2, "bent" },
    { 3, "bill" },
    {  0, "std" }, /* default */
    { RDMC_NA, NULL } /* end id */
  };

const rdmc_idtable_t rdmc_datatrans_idtable[] =
  {
    {  0, "coax" },
    {  1, "tp" },
    {  2, "opt" },
    {  3, "qt" },
    {  4, "tot" },
    { 10, "dig" },
    {  3, "std" }, /* default */
    { RDMC_NA, NULL } /* end id */
  };

const rdmc_idtable_t rdmc_detector_idtable[] =
{
  {  0, "unknown" },
  {  BAIKAL,  "baikal" }, /* and default baikal */
  {  NT_96,  "baikal-nt96" }, /* and default baikal */
  {  NT_36,  "baikal-nt36" },
  {  NT_36s, "baikal-nt36b" },
  {  NT_36s, "baikal-nt36'" },
  {  NT_72,  "baikal-nt72" },
  {  NT_144, "baikal-nt144" },
  {  NT_192, "baikal-nt192" },
  {  NT_200, "baikal-nt200" },
  {  AMANDA,    "amanda" },
  {  AMANDA_A,    "amanda-a" },
  {  AMANDA_B_4,  "amanda-b-4" },
  {  AMANDA_B_10, "amanda-b-10" },/* and default Amanda */
  {  AMANDA_B_11, "amanda-b-11" }, /* and default Amanda */
  {  AMANDA_B_13, "amanda-b-13" }, /* and default Amanda */
  {  AMANDA_II,   "amanda-ii-20" },
  {  AMANDA_II,   "amanda-ii" },
  {  AMANDA_KM3,  "icecube" },
  {  AMANDA_KM3,  "amanda-km3" },
  {  AMANDA, "neutrino_telescope" }, /* old f2000 from old rdmc */
  {  JULIA, "julia" },
  { RDMC_NA, NULL } /* end id */
};

const rdmc_idtable_t  rdmc_particle_idtable[] = {
  { 0        , "?" },
  { 0        , "unkwn" },
  { 0        , "unknown" },
  { RDMC_NA  , "?" },
  { MUON_PLUS,  "mu+"      },
  { MUON_MINUS, "mu-"      },
  { MUON,       "mu"       },
  { BREMS,      "brems"    },
  { BREMS,      "brehm"    },
  { DELTAE,     "delta"    },
  { PAIRPROD,   "epair"    },
  { MU_PAIR,    "mupair"   },
  { NUCL_INT,   "munu"    },
  { NUCL_INT,   "nucle"    },
  { HADRONS,    "hadr"     },
  { GAMMA,      "gamma"    },
  { E_PLUS,     "e+"       },
  { E_MINUS,    "e-"       },
  { NU,         "nu"    },
  { NU_MU,      "nu_mu"    },
  { NU_MU_BAR,  "~nu_mu"    },
  { NU_EL,      "nu_e"    },
  { NU_EL_BAR,  "~nu_e"    },
  { NU_TAU,     "nu_tau"    },
  { NU_TAU_BAR, "~nu_tau"    },
  { PI_0,       "pi0"      },
  { PI_PLUS,    "pi+"      },
  { PI_MINUS,   "pi-"      },
  { P_PLUS,     "p+"       },
  { P_MINUS,    "p-"       },
  { FIBERLASER, "flaser"   },
  { N2LASER,    "n2laser"  },
  { YAGLASER,   "yaglaser" },
  { 0        , "unkwn" },
  { 0        , "unknown" },
  { 0        , "?" },
  { RDMC_NA, NULL } /* end id */
};

