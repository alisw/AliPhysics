/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/// \ingroup macros
/// \file rootlogon.C
/// \brief Macro which is run when starting Root in MUON
///
/// It loads the MUON libraries needed for simulation and reconstruction
/// and sets the include path. 
///
/// \author Laurent Aphecetche

{
  gSystem->SetIncludePath("-I${ALICE_ROOT}/include");
}
