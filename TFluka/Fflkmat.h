#include "cfortran.h"
#include "Rtypes.h"

#include "Fdimpar.h"

extern "C" {
/*$ CREATE FLKMAT.ADD
*COPY FLKMAT
*
*=== Flkmat ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     Partial (some variables come from FLUKA87)                       *
*     Copyright (C) 1996-2005      by        Alfredo Ferrari           *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     FLuKa MATerial properties and atomic data                        *
*                                                                      *
*     Version for Fluka91/.../2005/...:                                *
*                                                                      *
*     Last change on  28-Apr-05    by  Alfredo Ferrari, INFN-Milan     *
*                                                                      *
*                                                                      *
*     This common contains the basic properties of the materials used  *
*     in the FLUKA run. Other properties are recorded in specialized   *
*     commons (ie for dE/dx etc)                                       *
*                                                                      *
*     Aocmbm(i) = Atomic density of the i_th material in barn^-1 cm^-1 *
*                 (Atoms Over Cm times Barn for Materials)             *
*     Eocmbm(i) = Electron density of the i_th material in barn^-1cm^-1*
*                 (Atoms Over Cm times Barn for Materials)             *
*       Amss(i) = Atomic weight (g/mole) of the i_th material          *
*     Amssem(i) = "Effective" i_th material atomic weight for the para-*
*                 metrized EM cascade                                  *
*        Rho(i) = Density of the i_th material                         *
*       Ztar(i) = Atomic number of the i_th material                   *
*     Zsqtar(i) = Squared atomic number of the i_th material           *
*     Ztarem(i) = "Effective" atomic number for the i_th material for  *
*                 the parametrized EM cascade                          *
*     Ainlng(i) = Inelastic scattering length of the i_th material     *
*                 for beam particles at the average beam energy in cm  *
*     Aellng(i) = Elastic scattering length of the i_th material for   *
*                 beam particles at average beam energy in cm          *
*      X0rad(i) = Radiation lengths of the materials in cm             *
*     Ainnth(i) = Inelastic scattering length of the i_th material     *
*                 for neutrons at threshold energy in cm               *
*     Medium(k) = Material number of the k_th region                   *
*     Mulflg(i) = Flags for multiple scattering options for the i_th   *
*                 material                                             *
*      Icomp(i) = Starting address in the Matnum array if the i_th     *
*                 material is a compound/mixture, 0 otherwise          *
*     Mssnum(i) = Mass number of the target nucleus for the i_th mater-*
*                 ial, if =< 0 it means that it is in the natural isot-*
*                 opic composition                                     *
*     Msindx(i) = Index for tabulations for the given isotope of the   *
*                 target nucleus (meaningful only for mssnum > 0)      *
*                 that it is in the natural isotopic composition       *
*     Lcmpnd(i) = logical flag for real compounds versus mixtures      *
*     Matnam(i) = Alphabetical name of the i_th material number        *
*        Nregs  = total number of regions                              *
*        Nregcg = total number of combinatorial geometry regions       *
*        Nmat   = total number of materials used in the problem        *
*        Mtbsnm = medium for which inelastic interaction biasing must  *
*                 be done                                              *
*                                                                      *
*                        Mxxmdf = maximum number of materials          *
*                        Mxxrgn = maximum number of regions            *
*                                                                      *
*----------------------------------------------------------------------*
*
      CHARACTER*8 MATNAM
      LOGICAL     LCMPND
      COMMON / FLKMAT / AOCMBM(MXXMDF), EOCMBM(MXXMDF), AMSS  (MXXMDF),
     &                  AMSSEM(MXXMDF), RHO   (MXXMDF), ZTAR  (MXXMDF),
     &                  ZTAREM(MXXMDF), ZSQTAR(MXXMDF), AINLNG(MXXMDF),
     &                  AELLNG(MXXMDF), X0RAD (MXXMDF), AINNTH(MXXMDF),
     &                  MEDIUM(MXXRGN), MULFLG(MXXMDF), ICOMP (MXXMDF),
     &                  MSSNUM(MXXMDF), MSINDX(MXXMDF), LCMPND(MXXMDF),
     &                  NREGS , NMAT  , MTBSNM, NREGCG
      COMMON / CHFLKM / MATNAM(MXXMDF)
      SAVE / FLKMAT /
      SAVE / CHFLKM /
*/

    typedef struct {
	Double_t aocmbm[mxxmdf];
	Double_t eocmbm[mxxmdf];
	Double_t amss  [mxxmdf];
	Double_t amssem[mxxmdf];
	Double_t rho   [mxxmdf];
	Double_t ztar  [mxxmdf];
	Double_t ztarem[mxxmdf];
	Double_t zsqtar[mxxmdf];
	Double_t ainlng[mxxmdf];
	Double_t aellng[mxxmdf];
	Double_t x0rad [mxxmdf];
	Double_t ainnth[mxxmdf];
	Int_t    medium[mxxrgn];
	Int_t    mulflg[mxxmdf];
	Int_t    icomp [mxxmdf];
	Int_t    mssnum[mxxmdf];
	Int_t    msindx[mxxmdf];
	Int_t    lcmpnd[mxxmdf];
	Int_t    nregs;
	Int_t    nmat;
	Int_t    mtbsnm;
	Int_t    nregcg;
    } flkmatCommon;
#define FLKMAT COMMON_BLOCK(FLKMAT,flkmat)
    COMMON_BLOCK_DEF(flkmatCommon, FLKMAT);
}
