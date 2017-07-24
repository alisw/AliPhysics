// -*- C++ -*-
#ifndef ALIADLOADER_H
#define ALIADLOADER_H

/////////////////////////////////////////////////////////////////////
//                                                                 //
// Base class for ADloaders.                                       //
// Loader provides base I/O facilities for standard data.          //
// Each detector has a loader data member.                         //
// Loader is always accessible via folder structure as well.       //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include "AliLoader.h"

class AliADLoader: public AliLoader {
public:
  AliADLoader();
  AliADLoader(const Char_t *name,const Char_t *topfoldername);
  AliADLoader(const Char_t *name,TFolder *topfolder);
  virtual ~AliADLoader() {};

  AliADLoader & operator=(const AliADLoader & ) {return *this;}

private:
  static const TString fgkDefaultHitsFileName;      //! Default Name for hit file
  static const TString fgkDefaultDigitsFileName;    //! Default Name for digit file

  ClassDef(AliADLoader,1);
};

#endif
