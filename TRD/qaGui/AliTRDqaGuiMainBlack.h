#ifndef ALITRDQAGUIMAINBLACK_H
#define ALITRDQAGUIMAINBLACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDqaGuiMainBlack.h 23387 2008-01-17 17:25:16Z cblume $ */

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
// of clusters. It lets display and browse throu histograms created by 
// the AliTRDQADataMakerRec run during the reconstruction 
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "TGFrame.h"

class TGWindow;
class TGTab;

class AliTRDqaGuiBlackSM;
class AliTRDqaGuiBlackChamber;  
class AliTRDqaGuiBlackError;
class AliTRDqaGuiBlackGlobal;
class AliTRDqaGuiBlackGTU;

class AliTRDqaGuiMainBlack : public TGMainFrame {

 public:

  AliTRDqaGuiMainBlack(TGWindow *parent);
  ~AliTRDqaGuiMainBlack() {}

  void SetQAFile(const char *file);

 protected:
  
  TGTab *fGTabPanel;                   // Something
   
  AliTRDqaGuiBlackSM       *fGSM;      // Something else
  AliTRDqaGuiBlackChamber  *fGChamber; // Something 
  AliTRDqaGuiBlackError    *fGError;  // somethig
  AliTRDqaGuiBlackGlobal   *fGGlobal;  // global view
  AliTRDqaGuiBlackGTU      *fGGTU;      // view status of optical links

 private:  
  AliTRDqaGuiMainBlack& operator = (const AliTRDqaGuiMainBlack& /*g*/) { return *this; };
  AliTRDqaGuiMainBlack(const AliTRDqaGuiMainBlack &);

  ClassDef(AliTRDqaGuiMainBlack,1)     // Does something

};

#endif
