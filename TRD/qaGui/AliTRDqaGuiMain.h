#ifndef ALITRDQAGUIMAIN_H
#define ALITRDQAGUIMAIN_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDqaGuiMain.h 23387 2008-01-17 17:25:16Z cblume $ */

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

class AliTRDqaGuiESDs;
class AliTRDqaGuiClustersStack;
class AliTRDqaGuiClustersSM;
class AliTRDqaGuiClusters;

class AliTRDqaGuiMain : public TGMainFrame {
 
 public:

  AliTRDqaGuiMain();
  AliTRDqaGuiMain(TGWindow *parent);
  ~AliTRDqaGuiMain() {}

  void SetQAFile(const char *file);
  
 protected:
  
  TGTab *fGTabPanel;                  // main tab panel
   
  AliTRDqaGuiClusters       *fGDet;   // panel with clusers
  AliTRDqaGuiClustersSM     *fGSM;    // panel with clusers
  AliTRDqaGuiClustersStack  *fGStack; // panel with clusers
  AliTRDqaGuiESDs *fGESDs[4];         // panel with ESDs

 private:
  AliTRDqaGuiMain& operator = (const AliTRDqaGuiMain& /*g*/) { return *this; };
  AliTRDqaGuiMain (const AliTRDqaGuiMain&);

  ClassDef(AliTRDqaGuiMain,1) //
};

#endif
