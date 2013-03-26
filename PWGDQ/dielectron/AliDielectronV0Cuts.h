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
  void InitEvent(AliVTrack *trk);
  Bool_t IsNewEvent(const AliVEvent *ev);
  virtual Bool_t IsSelected(TObject* track);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}
  void SetPdgCodes(Int_t mother, Int_t negDaughter, Int_t posDaughter) {fMotherPdg=mother; fNegPdg=negDaughter; fPosPdg=posDaughter;}
  void SetExcludeTracks(Bool_t exclude) {fExcludeTracks=exclude;}
  void SetDefaultPID(Int_t def) {fPID=def;}
  void ResetUniqueEventNumbers() { fOrbit=0; fPeriod=0; fBunchCross=0; }

private:

  TBits fV0TrackArr;                        // array with booleans where TrackID corresponds to bitnumber
  Bool_t fExcludeTracks;                       // cut logic: exclude or include tracks corresponding to a V0 candidate

  Int_t fMotherPdg;                         // target pdg code of the mother
  Int_t fNegPdg;                            // target pdg code of the negative daughter
  Int_t fPosPdg;                            // target pdg code of the positive daughter
  Int_t fPID;                               // default PID usage (see AliDielectronPID)

  // memebers needed to identify an event
  UInt_t fOrbit;                            // orbit number
  UInt_t fPeriod;                           // period number
  UShort_t fBunchCross;                     // bunch cross number

  AliDielectronV0Cuts(const AliDielectronV0Cuts &c);
  AliDielectronV0Cuts &operator=(const AliDielectronV0Cuts &c);

  ClassDef(AliDielectronV0Cuts,1)
};

#endif

