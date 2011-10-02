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

/// \ingroup macros
/// \file loadlibs.C
/// \brief Macro which loads the libraries needed for simulation and reconstruction
/// with MUON configuration macros
///
/// \author Christian Finck
///
/// New libraries list by Laurent Aphecetche

void loadlibs () 
{
  gSystem->Load("libVMC");
  gSystem->Load("libMinuit");
  gSystem->Load("libTree");
  gSystem->Load("libProofPlayer");
  gSystem->Load("libXMLParser");
  gSystem->Load("libPhysics");

  gSystem->Load("libSTEERBase"); 
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libCDB");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTEER");
  gSystem->Load("libANALYSISalice");
  
  gSystem->Load("libMUONcore");
  gSystem->Load("libMUONmapping");
  gSystem->Load("libMUONcalib");
  gSystem->Load("libMUONgeometry");
  gSystem->Load("libMUONtrigger");
  
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libMUONraw");
    
  gSystem->Load("libMUONbase");

  gSystem->Load("libMUONshuttle");

  gSystem->Load("libMUONrec");
  
  gSystem->Load("libRAWDatasim");
  gSystem->Load("libMUONsim");
  
  gSystem->Load("libMUONevaluation");
  gSystem->Load("libMUONcalign");
  
  gSystem->Load("libMUONgraphics");

  gSystem->Load("libEVGEN");
  gSystem->Load("libgeant321");
  gSystem->Load("libhijing");
  gSystem->Load("libFASTSIM");
  gSystem->Load("libTHijing");
  gSystem->Load("libpythia6");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libAliPythia6");
  
}
