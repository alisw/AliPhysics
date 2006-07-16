#ifndef ALIANALYSISRLCONTAINER_H
#define ALIANALYSISRLCONTAINER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 31/05/2006

//==============================================================================
//   AliAnalysysRLContainer - Special container working with AliRunLoader
//        with trees
//==============================================================================

#ifndef ALIANALYSISDATQCONTAINER_H
#include "AliAnalysisDataContainer.h"
#endif

class TFile;
class TTree;
class AliRunLoader;
class AliHeader;
class AliStack;

class AliAnalysisRLContainer : public AliAnalysisDataContainer {

public:
   AliAnalysisRLContainer();
   AliAnalysisRLContainer(const char *name);
   virtual ~AliAnalysisRLContainer();
   
//protected:
   AliRunLoader             *GetRunLoader();
   AliHeader                *GetHeader();
   AliStack                 *GetStack();
   TTree                    *GetKinematics();
   AliESD                   *GetESD() const {return fESD;}
   virtual void              GetEntry(Long64_t ientry);
   virtual Bool_t            SetData(TObject *data, Option_t *option="");
   // Send a notify signal to the container
   virtual void              NotifyChange(ENotifyMessage type);

private:
   void                      DeleteKinematicsFile();
   void                      DeleteRunLoader();   

   AliRunLoader*             fRunLoader;    //! pointer to the RunLoader if galice.root was opened
   AliESD*                   fESD;        //! "ESD" branch in fChain
   TFile*                    fKineFile;   //! pointer to Kinematics.root if the file was opened
   Bool_t                    fKinematicsLoaded;    // determines if the stack is properly loaded (AliRunLoader::LoadKinematics() succeeded), this is needed because the GetStack returnes a invalid stack object when the function failed
   Bool_t                    fHeaderLoaded;        // determines if the header is properly loaded

   ClassDef(AliAnalysisRLContainer,1)  // Class describing a data container using AliRunLoader
};
#endif
