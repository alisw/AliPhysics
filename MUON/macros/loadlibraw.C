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
/// \file loadlibraw.C
/// \brief Macro for loading libraries needed for reading & decoding rawdata
///
/// \author Christian Finck

void loadlibraw () 
{
  gSystem->Load("libMatrix");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  
  gSystem->Load("libESD");
  
  gSystem->Load("libHist");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  
  gSystem->Load("libMUONraw");
  
}
