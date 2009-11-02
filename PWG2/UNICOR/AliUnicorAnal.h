#ifndef ALIUNICORANAL_H
#define ALIUNICORANAL_H

/* Copyright(c) 1998-2048, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// Author: Dariusz Miskowiec <mailto:d.miskowiec@gsi.de> 2007

//=============================================================================
// parent class of all analyzers
//=============================================================================

#include <TNamed.h>
#include <TObjArray.h>
#include <TDatabasePDG.h>

class TCollection;

//=============================================================================
class AliUnicorAnal : public TNamed {
   
 public:
  AliUnicorAnal(char *nam="anal");                                       // constructor
  virtual ~AliUnicorAnal() {}                                            // destructor
  Long64_t Merge(const TCollection * const list);                // sumup histograms
  void     Save(const char *outfil, const char *mode="update");  // save histograms 
  TObject *GetHist(int i)                                        {return fHistos.At(i);}

 protected:
  static TDatabasePDG fgPDG;  // particle database
  TObjArray fHistos;          // histograms

  ClassDef(AliUnicorAnal,1)
};
//=============================================================================
#endif
