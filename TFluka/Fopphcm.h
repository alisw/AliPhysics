#ifndef FOPPHCM_H
#define FOPPHCM_H 1

#include "Rtypes.h"
#include "cfortran.h"
#include "Fdimpar.h"

extern "C" {
/*
*=== Opphcm ===========================================================*
*
*----------------------------------------------------------------------*
*                                                                      *
*     OPtical PHoton CoMmon:                                           *
*                                                                      *
*     Created on 19 september 1997 by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on 11-jan-99     by    Alfredo Ferrari               *
*                                                                      *
*       Opphpr (ip,im) = ip_th optical property parameter of the im_th *
*                        material (non metal)                          *
*                       ip =                                           *
*                            1: refraction index                       *
*                            2: absorption  coeff. (cm^-1)             *
*                            3: diffusion   coeff. (cm^-1)             *
*                            4: refraction index 1st derivative        *
*                            5: absorption 1st derivative              *
*                            6: diffusion  1st derivative              *
*                            7: refraction index 2nd derivative        *
*                            8: absorption 2nd derivative              *
*                            9: diffusion  2nd derivative              *
*                           10: refraction index 3rd derivative        *
*                           11: absorption 3rd derivative              *
*                           12: diffusion  3rd derivative              *
*                        metal:                                        *
*                       ip =                                           *
*                            1: refraction index (not used)            *
*                            2: absorption  coeff. (cm^-1) (not used)  *
*                            3: 1 - reflectivity index                 *
*                            7: 1 - reflectivity index 1st derivative  *
*                            9: 1 - reflectivity index 2nd derivative  *
*                           12: 1 - reflectivity index 3rd derivative  *
*          Emncer (im) = minimum energy for Cerenkov photon production *
*                        for im_th medium                              *
*          Emxcer (im) = maximum energy for Cerenkov photon production *
*                        for im_th medium                              *
*          Rmxcer (im) = maximum refractive index in the energy range  *
*                        of interest for Cerenkov photon production    *
*                        for im_th medium                              *
*       Escint (je,im) = energy for je_th scintillation photon produc- *
*                        tion for im_th medium                         *
*       Fscint (je,im) = fraction of energy emitted as the je_th scin- *
*                        tillation photon energy for im_th medium      *
*       Sscint (je,im) = sensitivity for the je_th scintillation photon*
*                        energy for im_th medium                       *
*       Tscint (je,im) = time constant of  je_th scintillation photon  *
*                        production  for im_th medium                  *
*          Emntrd (im) = minimum energy for transition radiation photon*
*                        production for im_th medium                   *
*          Emxtrd (im) = maximum energy for transition radiation photon*
*                        production for im_th medium                   *
*          Wvmnop (im) = minimum wavelength for opt. photon transport  *
*                        for im_th medium (default: 250 nm)            *
*          Wvcnop (im) = central wavelength for opt. photon transport  *
*                        for im_th medium (default: 589 nm, Na D)      *
*          Wvmxop (im) = maximum wavelength for opt. photon transport  *
*                        for im_th medium (default: 600 nm)            *
*          Ommnop (im) = minimum 2pi x freq. for opt. photon transport *
*                        for im_th medium                              *
*          Omcnop (im) = central 2pi x freq. for opt. photon transport *
*                        for im_th medium                              *
*          Ommxop (im) = maximum 2pi x freq. for opt. photon transport *
*                        for im_th medium                              *
*               Wvmnsn = minimum wavelength for opt. photon sensitivity*
*                        (default:  25 nm)                             *
*               Wvcnsn = central wavelength for opt. photon sensitivity*
*                        (default: 589 nm, Na D)                       *
*               Wvmxsn = maximum wavelength for opt. photon sensitivity*
*                        for im_th medium (default: 6000 nm)           *
*               Ommnsn = minimum 2pi x freq. for opt. photon sensiti-  *
*                        vity                                          *
*               Omcnsn = central 2pi x freq. for opt. photon sensiti-  *
*                        vity                                          *
*               Ommxsn = maximum 2pi x freq. for opt. photon sensiti-  *
*                        vity                                          *
*               Opsnmx = maximum of optical photon sensitivity         *
*          Rghnss (ib) = Roughness parameter for ib_th material-to-ma- *
*                        terial boundary ib_th                         *
*          M1rghn (ib) = 1st material of ib_th material-to-material    *
*                        boundary                                      *
*          M2rghn (ib) = 2nd material of ib_th material-to-material    *
*                        boundary                                      *
*          M1rgbx (ix) = 1st region of ix_th region-to-region special  *
*                        boundary                                      *
*          M2rgbx (ix) = 2nd region of ix_th region-to-region special  *
*                        boundary                                      *
*          Lopprp (im) = logical flag for optical properties of im_th  *
*                        material                                      *
*          Lopmtl (im) = logical flag whether the im_th optical mate-  *
*                        rial is a metal or not                        *
*          Lwvopp (im) = logical flag whether optical properties of    *
*                        im_th material are expressed as a function of *
*                        wavelength (true) or 2pi x frequency (false). *
*                        By default it is true.                        *
*               Lwvops = logical flag whether optical photon sensiti-  *
*                        vities are expressed as a function of wave-   *
*                        length (true) or 2pi x frequency (false).     *
*                        By default it is true.                        *
*          Lcrnkv (im) = logical flag for Cerenkov photon production   *
*                        for im_th material                            *
*          Lscntl (im) = logical flag for scintillation photon produ-  *
*                        ction for im_th material                      *
*       Ltscnt (je,im) = logical flag  for time constant for je_th     *
*                        scintill photon production in im_th medium    *
*          Ltrrad (im) = logical flag for transition radiation photon  *
*                        production for im_th material                 *
*          Lopphp (im) = logical flag for transition radiation photon  *
*                        production for im_th material                 *
*               Nxoppb = number of material boundaries for which the   *
*                        roughness has been defined                    *
*               Nxopbx = number of region boundaries for which the     *
*                        special user routine ophbdx should be called  *
*                                                                      *
*----------------------------------------------------------------------*

      PARAMETER ( MXOPSN =  4 )
      PARAMETER ( MXOPPR = 12 )
      PARAMETER ( MXOPPB = 20 )
      PARAMETER ( MXOPBX = 40 )
      PARAMETER ( MXSCPH =  3 )
      PARAMETER ( WVMNTR = 250.D-07 )
      PARAMETER ( WVCNTR = 589.D-07 )
      PARAMETER ( WVMXTR = 600.D-07 )
      LOGICAL LOPPRP, LOPMTL, LWVOPP, LCRNKV, LTRRAD, LSCNTL, LOPPHP,
     &        LWVOPS, LTSCNT
      COMMON / OPPHCM /  WVMNSN, WVCNSN, WVMXSN, OMMNSN, OMCNSN, OMMXSN,
     &                   OPSNMX, OPSNPR (MXOPSN),OPPHPR (MXOPPR,MXXMDF),
     &                EMNCER (MXXMDF), EMXCER (MXXMDF), RMXCER (MXXMDF),
     &                EMNTRD (MXXMDF), EMXTRD (MXXMDF), WVMNOP (MXXMDF),
     &                WVMXOP (MXXMDF), WVCNOP (MXXMDF), OMMNOP (MXXMDF),
     &                OMMXOP (MXXMDF), OMCNOP (MXXMDF), RGHNSS (MXOPPB),
     &                ESCINT (MXSCPH,MXXMDF),    FSCINT (MXSCPH,MXXMDF),
     &                SSCINT (MXSCPH,MXXMDF), TSCINT(MXSCPH,MXXMDF),
     &                M1RGHN (MXOPPB), M2RGHN (MXOPPB), M1RGBX (MXOPBX),
     &                M2RGBX (MXOPBX), LOPPRP (MXXMDF), LOPMTL (MXXMDF),
     &                LWVOPP (MXXMDF), LCRNKV (MXXMDF), LSCNTL (MXXMDF),
     &                LTRRAD (MXXMDF), LOPPHP (MXXMDF), LWVOPS, NXOPPB,
     &                NXOPBX, LTSCNT (MXSCPH,MXXMDF)
*/
    const Int_t mxopsn =  4;
    const Int_t mxoppr = 12;
    const Int_t mxoppb = 20;
    const Int_t mxopbx = 40;
    const Int_t mxscph =  3;
    const Double_t wvmntr = 250.e-07;
    const Double_t wvcntr = 589.e-07;
    const Double_t wvmxtr = 600.e-07;

    typedef struct {
	Double_t wvmnsn;
	Double_t wvcnsn;
	Double_t wvmxsn;
	Double_t ommnsn;
	Double_t omcnsn;
	Double_t ommxsn;
	Double_t opsnmx;
	Double_t opsnpr [mxopsn];
	Double_t opphpr [mxxmdf][mxoppr];
	Double_t emncer [mxxmdf];
	Double_t emxcer [mxxmdf];
	Double_t rmxcer [mxxmdf];
	Double_t emntrd [mxxmdf];
	Double_t emxtrd [mxxmdf];
	Double_t wvmnop [mxxmdf];
	Double_t wvmxop [mxxmdf];
	Double_t wvcnop [mxxmdf];
	Double_t ommnop [mxxmdf];
	Double_t ommxop [mxxmdf];
	Double_t omcnop [mxxmdf];
	Double_t rghnss [mxoppb];
	Double_t escint [mxxmdf][mxscph];
	Double_t fscint [mxxmdf][mxscph];
	Double_t sscint [mxxmdf][mxscph];
	Double_t tscint [mxxmdf][mxscph];
	Int_t    m1rghn [mxoppb];
	Int_t    m2rghn [mxoppb];
	Int_t    m1rgbx [mxopbx];
	Int_t    m2rgbx [mxopbx];
	Int_t    lopprp [mxxmdf];
	Int_t    lopmtl [mxxmdf];
	Int_t    lwvopp [mxxmdf];
	Int_t    lcrnkv [mxxmdf];
	Int_t    lscntl [mxxmdf];
	Int_t    ltrrad [mxxmdf];
	Int_t    lopphp [mxxmdf];
	Int_t    lwvops;
	Int_t    nxoppb;
        Int_t    nxopbx;
	Int_t    ltscnt [mxxmdf][mxscph];
    } opphcmCommon;
#define OPPHCM COMMON_BLOCK(OPPHCM,opphcm)
COMMON_BLOCK_DEF(opphcmCommon,OPPHCM);
}
#endif
