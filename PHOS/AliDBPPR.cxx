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
// DB Class for table PPR:
//
//*-- Author: Yves Schutz (SUBATECH)
//////////////////////////////////////////////////////////////////////////////
 
#include <stdlib.h> 
#include <iostream.h> 
 
// --- ROOT system ---
#include "TSQLServer.h" 
#include "TSQLRow.h" 
#include "TSQLResult.h" 
#include "TDatime.h" 
 
// --- Standard library ---
 
// --- AliRoot header files ---  
#include "AliDBPPR.h" 
 
ClassImp(AliDBPPR)
//____________________________________________________________________________
AliDBPPR::AliDBPPR()
{   
  fCurrentEntry = 1 ; // first valid entry (0 = labels)
  fNfields = 6 ;
  fFields = new TString[6]; 
  fFields[0] = "RUN"; 
  fFields[1] = "EVENT"; 
  fFields[2] = "DATE"; 
  fFields[3] = "SIMULATION"; 
  fFields[4] = "DIGITIZATION"; 
  fFields[5] = "RECONSTRUCTION"; 
}   
//____________________________________________________________________________
AliDBPPR::~AliDBPPR()
{   
  delete[] fFields ;
}   
//____________________________________________________________________________
void AliDBPPR::GetEntry(Option_t * opt)
{   
   // Retrieves one single row from the table   //  opt = first : retrieves first entry 
   //  opt = last  : retrieves last entry 
   //  opt = next  : retrieves next to current entry 
  TSQLServer * mysql = TSQLServer::Connect("mysql://ccmysql.in2p3.fr:3306/alice", "schutz", "po2hgwy") ;
  TSQLResult * result = mysql->Query("SELECT * FROM PPR") ; 
  Int_t count = result->GetRowCount() ; 
  if ( !strcmp(opt, "first") ) fCurrentEntry = 1 ; 
  if ( !strcmp(opt, "last") ) fCurrentEntry = count ; 
  if ( fCurrentEntry > count ) fCurrentEntry = 1 ; 
  Int_t i; 
  Int_t end = fCurrentEntry ; 
  TSQLRow * row = 0 ; 
  for ( i = 0 ; i < end ; i++ ) {
     fCurrentEntry++ ; 
    row = result->Next() ; 
  }   
    fRUN = atoi(row->GetField(0)) ; 
    fEVENT = atoi(row->GetField(1)) ; 
    fDATE = TDatime::TDatime(row->GetField(2)) ; 
    fSIMULATION = row->GetField(3) ; 
    fDIGITIZATION = row->GetField(4) ; 
    fRECONSTRUCTION = row->GetField(5) ; 
}   
