/********************************************************************
* TIsajetCint.h
********************************************************************/
#ifdef __CINT__
#error TIsajetCint.h/C is only for compilation. Abort cint.
#endif
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G__ANSIHEADER
#define G__DICTIONARY
#include "G__ci.h"
extern "C" {
extern void G__cpp_setup_tagtableTIsajetCint();
extern void G__cpp_setup_inheritanceTIsajetCint();
extern void G__cpp_setup_typetableTIsajetCint();
extern void G__cpp_setup_memvarTIsajetCint();
extern void G__cpp_setup_globalTIsajetCint();
extern void G__cpp_setup_memfuncTIsajetCint();
extern void G__cpp_setup_funcTIsajetCint();
extern void G__set_cpp_environmentTIsajetCint();
}


#include "TROOT.h"
#include "TMemberInspector.h"
#include "TIsajet.h"

#ifndef G__MEMFUNCBODY
#endif

extern G__linked_taginfo G__TIsajetCintLN_TClass;
extern G__linked_taginfo G__TIsajetCintLN_TObject;
extern G__linked_taginfo G__TIsajetCintLN_TNamed;
extern G__linked_taginfo G__TIsajetCintLN_TGenerator;
extern G__linked_taginfo G__TIsajetCintLN_AliRndm;
extern G__linked_taginfo G__TIsajetCintLN_TIsajet;

/* STUB derived class for protected member access */
