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
//_________________________________________________________________________
// Macros performing the full reconstruction chain starting from Digits
// Use Case : 
//          root> .L Reconstruction.C++
//          root> rec("RE", "PHOS EMCAL") --> does the reconstruction for 
//                                            PHOS and EMCAL and stores the 
//                                            reconstructed particles in 
//                                            AliESDs.root
// author  : Yves Schutz (CERN/SUBATECH)
// February 2004
//_________________________________________________________________________
 
#include "AliReconstruction.h"
#include "TString.h"
#include "Riostream.h"
#include "AliPHOSGetter.h"
#include "AliEMCALGetter.h"

void reco(TString opt="TVRE", TString name="all", Bool_t debug="kFALSE") 
{
  AliReconstruction rec ; 
  if ( !opt.Contains("T") ) 
    rec.SetRunTracking("") ;
  else
    rec.SetRunTracking(name.Data()) ;

  if ( !opt.Contains("V") ) 
    rec.SetRunVertexFinder(kFALSE) ; 

  if ( opt.Contains("R") ) 
    rec.SetRunLocalReconstruction(name.Data()) ; 
  else 
    rec.SetRunLocalReconstruction("") ;

  if ( !opt.Contains("E") )
    rec.SetFillESD("") ; 
  else 
    rec.SetFillESD(name.Data()) ; 

  rec.Run() ;
}   
