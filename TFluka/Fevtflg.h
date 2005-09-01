extern "C" {
//*$ create evtflg.add
//*copy evtflg
//*
//*=== evtflg ===========================================================*
//*
//*----------------------------------------------------------------------*
//*                                                                      *
//*     event flags:                                                     *
//*                                                                      *
//*     created on    19 may 1998    by    alfredo ferrari & paola sala  *
//*                                                   infn - milan       *
//*                                                                      *
//*     last change on   13-aug-99   by    alfredo ferrari               *
//*                                                                      *
//*                                                                      *
//*----------------------------------------------------------------------*
//*
//*

typedef struct {
   Int_t    lelevt;
   Int_t    linevt;
   Int_t    ldecay;
   Int_t    ldltry;
   Int_t    lpairp;
   Int_t    lbrmsp;
   Int_t    lannrs;
   Int_t    lannfl;
   Int_t    lphoel;
   Int_t    lcmptn;
   Int_t    lcohsc;
   Int_t    llensc;
   Int_t    loppsc;
   Int_t    leldis;
    Int_t    lrdcay;
   Int_t    ntrcks;
} evtflgCommon;
#define EVTFLG COMMON_BLOCK(EVTFLG,evtflg)
COMMON_BLOCK_DEF(evtflgCommon,EVTFLG);
}
