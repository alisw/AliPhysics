/*
 *###################################################################
 *#        EPOS 1.67     K. WERNER, T. PIEROG, S. PORTEBOEUF.       #
 *#                      Contact: werner@subatech.in2p3.fr          #
 *###################################################################
 *
 * eposproc.h
 * 
 * Definitions of EPOS subroutines used in the interface.
 *
 *      Author: Piotr Ostrowski, postrow@if.pw.edu.pl
 */


#ifndef ROOT_Eposproc
#define ROOT_Eposproc

#ifndef WIN32
#define aaset aaset_
#define atitle atitle_
#define xiniall xiniall_
#define aread aread_
#define utpri utpri_
#define IniModel inimodel_
#define ainit ainit_
#define aseed aseed_
#define xana xana_
#define bstora bstora_
#define bstore bstore_
#define ustore ustore_
#define evgen evgen_
#define astati astati_
#define bfinal bfinal_
#define setinp setinp_

#define type_of_call

#else

#define aaset AASET
#define atitle ATITLE
#define xiniall XINIALL
#define aread AREAD
#define utpri UTPRI
#define IniModel INIMODEL
#define ainit AINIT
#define aseed ASEED
#define xana XANA
#define bstora BSTORA
#define bstore BSTORE
#define ustore USTORE
#define evgen EVGEN
#define astati ASTATI
#define bfinal BFINAL
#define setinp SETINP
#define type_of_call _stdcall

#endif

#include <TROOT.h>
#ifndef WIN32

extern "C" void type_of_call aaset(const Int_t &);
extern "C" void type_of_call atitle();
extern "C" void type_of_call xiniall();
extern "C" void type_of_call aread();
extern "C" void type_of_call utpri(const char *, Int_t &, Int_t &, const Int_t &, int);
extern "C" void type_of_call IniModel(Int_t &);
extern "C" void type_of_call ainit();
extern "C" void type_of_call aseed(const Int_t &);
extern "C" void type_of_call xana();
extern "C" void type_of_call bstora();
extern "C" void type_of_call bstore();
extern "C" void type_of_call ustore();
extern "C" void type_of_call evgen(Int_t &);
extern "C" void type_of_call astati();
extern "C" void type_of_call bfinal();
extern "C" void type_of_call setinp(const char *, Int_t &, int);
#else
#endif

#endif
