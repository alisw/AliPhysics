/***************************************************************************
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

//-----------------------------------------------------//
//                                                     //
//           Date   : August 05 2003                   //
//  This does photon hadron discrimination on the      // 
//  of matching a PMD cluster with a CPV cluster       //
//                                                     //
//-----------------------------------------------------//

#include <Riostream.h>
#include <TMath.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TParticle.h>

#include "AliLog.h"
#include "AliPMDrecdata.h"
#include "AliPMDclupid.h"
#include "AliPMDDiscriminator.h"
#include "AliPMDEmpDiscriminator.h"

ClassImp(AliPMDEmpDiscriminator)

AliPMDEmpDiscriminator::AliPMDEmpDiscriminator()
{
//
// Default Constructor
//
}
// -----------------------------------------------------------------------
AliPMDEmpDiscriminator::~AliPMDEmpDiscriminator()
{
  // Destructor
}

// -----------------------------------------------------------------------

void AliPMDEmpDiscriminator::Discrimination(TObjArray *pmdcontin, TObjArray *pmdcontout)
{
  // Does Photon/Hadron discrimination
  // matching the clusters of CPV and PREshower plane
  //

  Int_t   det,smn, trno, trpid, mstatus;
  Float_t clusdata[7];

  AliPMDrecdata  *pmdcl    = 0;
  AliPMDclupid   *pmdclout = 0;

  Int_t nentries1 = pmdcontin->GetEntries();

  AliDebug(1,Form("Number of total clusters from CPV PRE = %d",nentries1));


  for (Int_t ient1 = 0; ient1 < nentries1; ient1++)
    {
      pmdcl = (AliPMDrecdata*)pmdcontin->UncheckedAt(ient1);
      
      det         = pmdcl->GetDetector();
      smn         = pmdcl->GetSMNumber();
      trno        = pmdcl->GetClusTrackNo();
      trpid       = pmdcl->GetClusTrackPid();
      clusdata[0] = pmdcl->GetClusX();
      clusdata[1] = pmdcl->GetClusY();
      clusdata[2] = pmdcl->GetClusADC();
      clusdata[3] = pmdcl->GetClusCells();
      clusdata[4] = pmdcl->GetClusSigmaX();
      clusdata[5] = pmdcl->GetClusSigmaY();

      if (det == 0 && clusdata[2] > 300.)
	{
	  clusdata[6] = 1;     // photon
	}
      else
	{
	  clusdata[6] = 8;     // hadron
	}

      mstatus = 0;             // at this moment matching is not done

      pmdclout = new AliPMDclupid(det,smn,trno,trpid,mstatus,clusdata);
      pmdcontout->Add(pmdclout);
      
    } // Entries of TObjArray loop

}
// -----------------------------------------------------------------------
