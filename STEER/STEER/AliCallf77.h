#ifndef AliCallf77_H
#define AliCallf77_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#ifndef WIN32
# define type_of_call
# define DEFCHARD     const char* 
# define DEFCHARL   , const int 
# define PASSCHARD(string) string 
# define PASSCHARL(string) , strlen(string) 
#else
# define type_of_call  _stdcall
# define DEFCHARD   const char* , const int        
# define DEFCHARL          
# define PASSCHARD(string) string, strlen(string) 
# define PASSCHARL(string) 
#endif
#endif
