#ifndef ALIINFOFRAME_H
#define ALIINFOFRAME_H
/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////////
// ALICE INFO FRAME CLASS                                              //
// Author: Mayeul   ROUSSELET                                          //
// e-mail: Mayeul.Rousselet@cern.ch                                    //
// Last update:26/08/2003                                              //
/////////////////////////////////////////////////////////////////////////

#include <Rtypes.h>
#include <RQ_OBJECT.h>

class TGCompositeFrame;
class TGLabel;

class AliInfoFrame{
  //This class implements the info frame where the number of particles... are displayed

public:
	
 AliInfoFrame(TGCompositeFrame *p, UInt_t w, UInt_t h);
 virtual ~AliInfoFrame(void);
 
 void			AddLabel(const char *text, UInt_t options);
 TGCompositeFrame	*GetInfoFrame() const {return fMainFrame;};
 void 			Update();
 
private:

 TGCompositeFrame	*fMainFrame; // Main frame
 TGCompositeFrame	*fTitleFrame; // Title frame
 TGCompositeFrame	*fFiguresFrame; // Fugures frame
 TGLabel       		*fNbParticuleLabel; // Label for particle number
 TGLabel       	      	*fNbEventLabel; // Label for event number
 TGLabel       	       	*fNbHitsLabel; // Label for hits number
 TGLabel                *fNbClustersLabel; // Label for clusters number

 RQ_OBJECT("AliInfoFrame")

 ClassDef(AliInfoFrame,0);
};

#endif
