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
// DB Class for table PPRS:
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
#include "AliDBPPRS.h" 
 
ClassImp(AliDBPPRS)
//____________________________________________________________________________
AliDBPPRS::AliDBPPRS()
{   
  fCurrentEntry = 1 ; // first valid entry (0 = labels)
  fNfields = 12 ;
  fFields = new TString[12]; 
  fFields[0] = "RUN"; 
  fFields[1] = "EVENT"; 
  fFields[2] = "WORKER"; 
  fFields[3] = "STATUS"; 
  fFields[4] = "DATE"; 
  fFields[5] = "STORAGE"; 
  fFields[6] = "ID"; 
  fFields[7] = "POS"; 
  fFields[8] = "SIZE"; 
  fFields[9] = "FTP"; 
  fFields[10] = "LOG"; 
  fFields[11] = "COMMENT"; 
}   
//____________________________________________________________________________
AliDBPPRS::~AliDBPPRS()
{   
  delete[] fFields ;
}   
//____________________________________________________________________________
void AliDBPPRS::GetEntry(Option_t * opt)
{   
   // Retrieves one single row from the table   //  opt = first : retrieves first entry 
   //  opt = last  : retrieves last entry 
   //  opt = next  : retrieves next to current entry 
  TSQLServer * mysql = TSQLServer::Connect("mysql://ccmysql.in2p3.fr:3306/alice", "schutz", "po2hgwy") ;
  TSQLResult * result = mysql->Query("SELECT * FROM PPRS") ; 
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
    fWORKER = row->GetField(2) ; 
    fSTATUS = row->GetField(3) ; 
    fDATE = TDatime::TDatime(row->GetField(4)) ; 
    fSTORAGE = row->GetField(5) ; 
    fID = row->GetField(6) ; 
    fPOS = atoi(row->GetField(7)) ; 
    fSIZE = atoi(row->GetField(8)) ; 
    fFTP = row->GetField(9) ; 
    fLOG = row->GetField(10) ; 
    fCOMMENT = row->GetField(11) ; 
}   
