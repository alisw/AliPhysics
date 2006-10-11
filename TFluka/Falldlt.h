#ifndef FALLDLT_H
#define FALLDLT_H 

#include "cfortran.h"
#include "Rtypes.h"
extern "C" {
//*$ CREATE ALLDLT.ADD
//*COPY ALLDLT
//*
//*=== Alldlt ===========================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     Copyright (C) 2005-2006         by        Alfredo Ferrari        *
//*     All Rights Reserved.                                             *
//*                                                                      *
//*                                                                      *
//*     Include file: alldlt  (ALL DeLTas)  vv                           *
//*                                                                      *
//*     Created  on  10 october 2005    by        Alfredo Ferrari        *
//*                                                INFN - Milan          *
//*                                                                      *
//*     Last change on  19-feb-06       by        Alfredo Ferrari        *
//*                                                                      *
//*     Included in the following routines:                              *
//*                                                                      *
//*              blockmvax/bdtrns.f                                      *
//*              dedxmvax/dedxfl.f                                       *
//*              dedxmvax/enion.f                                        *
//*              dedxmvax/enionf.f                                       *
//*              emfmvax/ededxf.f                                        *
//*              emfmvax/emenio.f                                        *
//*              emfmvax/emfsco.f                                        *
//*              emfmvax/pdedxf.f                                        *
//*              kaskadmvax/kaskad.f                                     *
//*                                                                      *
//*           Talldl (m) = kinetic energy of the m_th primary electron   *
//*                        emitted during energy loss fluctuation pro-   *
//*                        cesses                                        *
//*               Tallmn = minimum energy of the recorded primary        *
//*                        electrons
//*       X/Y/Zalldl (m) = position coord. of the m_th primary electron  *
//*                        emitted during energy loss fluctuation pro-   *
//*                        cesses                                        *
//*
//*               Nalldl = number of recorded primary electrons          *
//*               Lalldl = logical flag for primary electrons recording  *
//*                                                                      *
//*----------------------------------------------------------------------*
//*
//      PARAMETER ( MXALLD = 5000 )
//*
//      LOGICAL LALLDL
//         COMMON / ALLDLT / TALLDL (MXALLD), XALLDL (MXALLD),
//     &                  YALLDL (MXALLD), ZALLDL (MXALLD),
//     &                  TALLMN, NALLDL, LALLDL
//      SAVE / ALLDLT /

    const Int_t mxalld = 5000;
    typedef struct {
	Double_t talldl[mxalld];
	Double_t xalldl[mxalld];
	Double_t yalldl[mxalld];
	Double_t zalldl[mxalld];
	Double_t tallmn;
	Int_t    nalldl;
	Int_t    lalldl;
    } alldltCommon;
#define ALLDLT COMMON_BLOCK(ALLDLT,alldlt)
    COMMON_BLOCK_DEF(alldltCommon,ALLDLT);
}

#endif
