/**
 * \file AliEmcalSetupTask.h
 * \brief Declaration of the EMCAL setup task.
 * \author Constantin Loizides <>, Lawrence Berkeley National Laboratory
 */
#ifndef ALIEMCALSETUPTASK_H
#define ALIEMCALSETUPTASK_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TClonesArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

/**
 * \class AliEmcalSetupTask
 * \ingroup EMCALFWTASKS
 * \brief Simple task setting up connections to databases for the EMCAL train
 *
 * This class, as an analysis task, handles the setup of the connection to the
 * databases (OCDB and OADB) and the EMCAL geometry. For the geometry one can
 * either use the normal raw OCDB on alien, then this has to be specified, or
 * one uses a snapshot from AliPhysics providing a minimal set of information
 * necessary for the EMCAL train. The snapshot is the default method. As tender
 * task, this task is expected to be the first task in the order of task attached
 * to the train. Note that the setup is performed only for the first event, so
 * consequently the setup task cannot handle several runs.
 */
class AliEmcalSetupTask : public AliAnalysisTaskSE {
 public:
  AliEmcalSetupTask();
  AliEmcalSetupTask(const char *name);
  virtual ~AliEmcalSetupTask();

  /**
   * Set the path to the geometry file
   * @param n Name of the geometry file
   */
  void               SetGeoPath(const char *n)  { fGeoPath  = n; }
  /**
   * Switch odd handling of the OCDB
   * @param b If true the task doesn't handle the OCDB
   */
  void               SetNoOCDB(Bool_t b)        { fNoOCDB   = b; }
  /**
   * Set the path where to find the OADB
   * @param n Path of the OADB
   */
  void               SetOadbPath(const char *n) { fOadbPath = n; }
  /**
   * Set the path where to find the OCDB
   * @param n Path of the OCDB
   */
  void               SetOcdbPath(const char *n) { fOcdbPath = n; }
  /**
   * Define which detectors (entries) to handle from the OCDB
   * @param n List, spearated by whitespace, with detectors (entries)
   */
  void               SetObjs(const char *n)     { fObjs     = n; }

 protected:
  void               ConnectInputData(Option_t *option = "");
  void               UserExec(Option_t *option);
  void               Setup(Int_t runno);
  void               Terminate(Option_t *option);

  TString            fOcdbPath;        ///< path to ocdb (def=uselocal)
  TString            fOadbPath;        ///< path to oadb
  TString            fGeoPath;         ///< path to geometry
  TString            fObjs;            ///< string of objects for alignment to apply
  Bool_t             fNoOCDB;          ///< if true then do not mess with OCDB
  Bool_t             fIsInit;          //!<!=true then already initialized
  TString            fLocalOcdb;       //!<!directory path to local ocdb
  TString            fLocalOcdbStor;   //!<!storage path to local ocdb

 private:
  AliEmcalSetupTask(const AliEmcalSetupTask&);            // not implemented
  AliEmcalSetupTask &operator=(const AliEmcalSetupTask&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliEmcalSetupTask, 6); // Class to setup geometry for EMCal
  /// \endcond
};

#endif
