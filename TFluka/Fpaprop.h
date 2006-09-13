#ifndef FPAPROP_H
#define FPAPROP_H 1

#include "Rtypes.h"
#include "cfortran.h"

#include "Fdimpar.h"

extern "C" {
  //*$ CREATE PAPROP.ADD
  //*COPY PAPROP
  //*
  //*=== paprop ===========================================================*
  //*
  //*----------------------------------------------------------------------*
  //*                                                                      *
  //*     PArticle PROPerties:                                             *
  //*                                                                      *
  //*     !!!!  N E W   V E R S I O N   F O R   FLUKA9x/FLUKA200x  !!!!    *
  //*                                                                      *
  //*     Created on    07 may 1991    by    Alfredo Ferrari & Paola Sala  *
  //*                                                   Infn - Milan       *
  //*                                                                      *
  //*     Last change on 03-jul-97     by    Alfredo Ferrari               *
  //*                                                                      *
  //*     Variable description:                                            *
  //*                                                                      *
  //*        btype (i) = literal name of the i_th particle                 *
  //*        am    (i) = i_th particle mass (GeV)                          *
  //*        amdisc(i) = i_th particle mass (GeV) for energy conservation  *
  //*                    purposes when discarded                           *
  //*        ichrge(i) = electric charge of the i_th particle              *
  //*        ibarch(i) = baryonic charge of the i_th particle              *
  //*        iscore(j) = id for the j_th scored distribution               *
  //*        gnname(k) = name of the k_th generalized particle type        * 2006.3
  //*        ijdisc(i) = flag for discarding the i_th particle type        *
  //*        tmnlf (i) = mean (not half!) life of the i_th particle (s)    *
  //*        biasdc(i) = decay biasing factor for the i_th particle        *
  //*        biasin(i) = inelastic interaction biasing factor for the i_th *
  //*                    particle                                          *
  //*        lhadro(i) = True if the i_th particle type is a hadron        *
  //*        jspinp(i) = i_th particle spin (in units of 1/2)              *
  //*        iparty(i) = i_th particle parity (when meaningful)            *
  //*        iparid(i) = particle type id flag for the i_th particle       *
  //*        lbsdcy(i) = logical flag for biased decay for the i_th parti- *
  //*                    cle: if .true. the biasing factor is used as an   *
  //*                    upper limit to the decay length                   *
  //*        lprbsd    = logical flag for biased decay: if .true. the bia- *
  //*                    sing factor is applied only to primaries          *
  //*        lprbsi    = logical flag for inelastic interaction biasing:   *
  //*                    if .true. the biasing factor is applied only to   *
  //*                    primaries                                         *
  //*        lsclwf    = logical flag for low energy neutron fission sco-  *
  //*                    ring                                              *
  //*        lscnbl    = logical flag for neutron balance scoring          *
  //*                                                                      *
  //*----------------------------------------------------------------------*


//    const Int_t mxgnpr =  33; // 2006.3
    const Int_t mxgnpr =  35;
    typedef struct {
	Double_t am[nallwp+7];         //(-6:NALLWP)
	Double_t amdisc[nallwp+7];     //(-6:NALLWP)
	Double_t tmnlf[nallwp+7];      //(-6:NALLWP)
	Double_t biasdc[nallwp+7];     //(-6:NALLWP)
	Double_t biasin[nallwp+7];     //(-6:NALLWP)
	Int_t    ichrge[nallwp+7];     //(-6:NALLWP)
	Int_t    ibarch[nallwp+7];     //(-6:NALLWP)
	Int_t    ijdisc[nallwp+7];     //(-6:NALLWP)
	Int_t    jspinp[nallwp+7];     //(-6:NALLWP)
	Int_t    iparty[nallwp+7];     //(-6:NALLWP)
	Int_t    iparid[nallwp+7];     //(-6:NALLWP)
	Int_t    lhadro[nallwp+7];     //(-6:NALLWP)
	Int_t    lbsdcy[nallwp+7];     //(-6:NALLWP)
	Int_t    iscore[12];
	Int_t    lprbsd;
	Int_t    lprbsi;
	Int_t    lsclwf;
	Int_t    lscnbl;
    } papropCommon;
#define PAPROP COMMON_BLOCK(PAPROP,paprop)
COMMON_BLOCK_DEF(papropCommon,PAPROP);

    typedef struct {
	Char_t   btype[nallwp+7][8];     //(-6:NALLWP)
//	Char_t   genpar[mxgnpr][8];          // 2006.3 
	Char_t   gnname[mxgnpr][8];          // 2006.3
    } chpprpCommon;
#define CHPPRP COMMON_BLOCK(CHPPRP,chpprp)
COMMON_BLOCK_DEF(chpprpCommon,CHPPRP);
}
#endif
