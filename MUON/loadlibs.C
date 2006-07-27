/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

// Macro which loads the libraries needed for simulation and reconstruction
// with MUON configuration macros
// Christian Finck

void loadlibs () 
{
  gSystem->Load("libPhysics");
  gSystem->Load("libmicrocern");
  gSystem->Load("libpdf");
  gSystem->Load("libpythia6");
  gSystem->Load("libEG");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libEGPythia6");

  gSystem->Load("libESD");
  gSystem->Load("libSTEER");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libRAWDatasim");
  gSystem->Load("libEVGEN");
  gSystem->Load("libFASTSIM");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libSTRUCT");
  gSystem->Load("libMUONmapping");
  gSystem->Load("libMUONgeometry");
  gSystem->Load("libMUONbase");
  gSystem->Load("libMUONraw");

  gSystem->Load("libMUONsim");
  new AliRun("gAlice","The ALICE Off-line Simulation Framework");

  gSystem->Load("libMUONrec");
}
