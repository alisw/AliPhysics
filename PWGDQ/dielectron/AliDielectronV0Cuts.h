#ifndef ALIDIELECTRONV0CUTS_H
#define ALIDIELECTRONV0CUTS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronV0Cuts                        #
//#         Provide cuts for all variables handled in         #
//#           AliDielectronV0Manager                         #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################

#include <Rtypes.h>

#include <AliDielectronVarCuts.h>

class AliDielectronV0Cuts : public AliDielectronVarCuts {
public:

  AliDielectronV0Cuts();
  AliDielectronV0Cuts(const char* name, const char* title);
  virtual ~AliDielectronV0Cuts();
  //TODO: make copy constructor and assignment operator public

  //
  //Analysis cuts interface
  //
  virtual void Init();
  virtual Bool_t IsSelected(TObject* track);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}

private:

  TArrayC *fV0TrackArr;                        // array where TrackID corresponds to index

  AliDielectronV0Cuts(const AliDielectronV0Cuts &c);
  AliDielectronV0Cuts &operator=(const AliDielectronV0Cuts &c);

  ClassDef(AliDielectronV0Cuts,0)
};

#endif

