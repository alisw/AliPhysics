#ifndef FBEAMCM_H
#define FBEAMCM_H 1

#include "cfortran.h"
#include "Rtypes.h"
extern "C" {

//*=== beam =============================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     CoMmon for BEAM properties:                                      *
//*                                                                      *
//*        Pbeam  = average beam particle momentum (GeV/c)               *
//*        Pbmmax = maximum momentum for which tabulations must be       *
//*                 generated (GeV/c)                                    *
//*        Dpbeam = beam momentum spread (GeV/c)                         *
//*        Divbm  = beam angular divergense (mrad)                       *
//*        Xspot  = beam width in (beam frame) x-direction (cm)          *
//*        Yspot  = beam width in (beam frame) y-direction (cm)          *
//*        Xbeam  = beam spot centre (geom frame) x-coordinate (cm)      *
//*        Ybeam  = beam spot centre (geom frame) y-coordinate (cm)      *
//*        Zbeam  = beam spot centre (geom frame) z-coordinate (cm)      *
//*        Ubeam  = beam direction cosine wrt the (beam frame) x-axis    *
//*        Vbeam  = beam direction cosine wrt the (beam frame) y-axis    *
//*        Wbeam  = beam direction cosine wrt the (beam frame) z-axis    *
//*        Ubmpol = beam polarization cosine wrt the (beam frame) x-axis *
//*        Vbmpol = beam polarization cosine wrt the (beam frame) y-axis *
//*        Wbmpol = beam polarization cosine wrt the (beam frame) z-axis *
//*        Polfra = polarization fraction                                *
//*        Rflood = emission radius for a uniform and isotropic source   *
//*                 or maximum radius for a cylindrical/spherical volume *
//*                 source                                               *
//*        Rvlmax = emission radius for a uniform and isotropic source   *
//*                 or maximum radius for a cylindrical/spherical volume *
//*                 source                                               *
//*        Rvlmin = minimum radius for a cylindrical/spherical volume    *
//*                 source                                               *
//*        Dxvlmx = maximum Dx for a cartesian volume source             *
//*                (particle emitted inside [Xina+Dxvlmn/2,Xina+Dxvlmx/2]*
//*                 and inside [Xina-Dxvlmx/2, Xina-Dxvlmn/2])           *
//*        Dxvlmn = minimum Dx for a cartesian volume source             *
//*        Dyvlmx = maximum Dy for a cartesian volume source             *
//*                (particle emitted inside [Yina+Dyvlmn/2,Yina+Dyvlmx/2]*
//*                 and inside [Yina-Dyvlmx/2, Yina-Dyvlmn/2])           *
//*        Dyvlmn = minimum Dy for a cartesian volume source             *
//*        Dzvlmx = maximum Dz for a cartesian/cylindrical volume source *
//*                (particle emitted inside [Zina+Dzvlmn/2,Zina+Dzvlmx/2]*
//*                 and inside [Zina-Dzvlmx/2, Zina-Dzvlmn/2])           *
//*        Dzvlmn = minimum Dz for a cartesian/cylindrical volume source *
//*        Ijbeam = beam particle type (see btype in /paprop/)           *
//*        Ijhion = heavy ion type if ijbeam = -2                        *
//*        Ldpgss = true for a gaussian momentum distribution of the     *
//*                 beam particles, false for a rectangular one          *
//*        Ldvgss = true for a gaussian angular divergence distribution  *
//*                 of the beam particles, false for a rectangular one   *
//*        Ldxgss = true for a gaussian spatial distribution of the beam *
//*                 spot in the x-direction, false for a rectangular one *
//*        Ldygss = true for a gaussian spatial distribution of the beam *
//*                 spot in the y-direction, false for a rectangular one *
//*        Beawei = weight of the beam particles                         *
//*        Lbeamc = flag for an annular beam                             *
//*        Lpperp = flag for polar. perp. to the beam direction          *
//*        Lpfrac = flag for interpreting the polar. fraction            *
//*   Bmaxis(j,i) = j_th component of the i_th axis used to define the   *
//*                 conventional x,y,z beam reference frame              *
//*!!!!! ATTENTION: in C++ it is the component bmaxis(i,j) !!!!!         *
//*        Lbaxis = logical flag for using a beam axis frame different   *
//*                 from the standard one                                *
//*        Lflood = logical flag for using a uniform and isotropic beam  *
//*                 source out of a sphere of radius Rflood              *
//*        Lvlcar = logical flag for using a cartesian   volume source   *
//*        Lvlcyl = logical flag for using a cylindrical volume source   *
//*        Lvlsph = logical flag for using a spherical   volume source   *
//*        Lsourc = logical flag for a user written source routine       *    
//*                                                                      *
//*----------------------------------------------------------------------*

typedef struct {
   Double_t pbeam;
   Double_t dpbeam;
   Double_t pbmmax;
   Double_t divbm;
   Double_t xspot;
   Double_t yspot;
   Double_t xbeam;
   Double_t ybeam;
   Double_t zbeam;
   Double_t ubeam;
   Double_t vbeam;
   Double_t wbeam;
   Double_t ubmpol;
   Double_t vbmpol;
   Double_t wbmpol;
   Double_t polfra;
   Double_t beawei;
   Double_t bmaxis[3][3];
   Double_t rvlmin;
   Double_t rvlmax;
   Double_t dxvlmn;
   Double_t dxvlmx;
   Double_t dyvlmn;
   Double_t dyvlmx;
   Double_t dzvlmn;
   Double_t dzvlmx;
   Int_t    ijbeam;
   Int_t    ijhion;
   Int_t    ldpgss;
   Int_t    ldvgss;
   Int_t    ldxgss;
   Int_t    ldygss;
   Int_t    lbeamc;
   Int_t    lpperp;
   Int_t    lpfrac;
   Int_t    lbaxis;
   Int_t    lflood;
   Int_t    lvlcar;
   Int_t    lvlcyl;
   Int_t    lvlsph;
   Int_t    lsourc;
} beamcmCommon;
#define BEAMCM COMMON_BLOCK(BEAMCM,beamcm)
COMMON_BLOCK_DEF(beamcmCommon,BEAMCM);
}

#endif
