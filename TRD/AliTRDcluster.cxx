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

/*
$Log$
Revision 1.4.10.2  2002/07/24 10:09:30  alibrary
Updating VirtualMC

Revision 1.6  2002/06/12 09:54:35  cblume
Update of tracking code provided by Sergei

Revision 1.4  2001/05/07 08:08:05  cblume
Update of TRD code

Revision 1.3  2000/12/08 16:07:02  cblume
Update of the tracking by Sergei

Revision 1.2  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.2.1  2000/09/22 14:47:52  cblume
Add the tracking code

*/

#include "AliTRDcluster.h"
#include "AliTRDgeometry.h"
#include "AliTRDrecPoint.h"

ClassImp(AliTRDcluster)
 

//_____________________________________________________________________________
  AliTRDcluster::AliTRDcluster(const AliTRDrecPoint &p):AliCluster()
{
  //
  // Constructor from AliTRDrecPoint
  //

  fDetector   = p.GetDetector();
  fTimeBin    = p.GetLocalTimeBin();

  fTracks[0]  = p.GetTrackIndex(0);
  fTracks[1]  = p.GetTrackIndex(1);
  fTracks[2]  = p.GetTrackIndex(2);

  fQ          = p.GetEnergy();

  fY          = p.GetY();
  fZ          = p.GetZ();
  fSigmaY2    = p.GetSigmaY2();
  fSigmaZ2    = p.GetSigmaZ2();  

  fSigmaY2    = 0.2;
  fSigmaZ2    = 5.;  

}

//_____________________________________________________________________________
AliTRDcluster::AliTRDcluster(const AliTRDcluster &c):AliCluster()
{
  //
  // Copy constructor 
  //

  fTracks[0]  = c.GetLabel(0);
  fTracks[1]  = c.GetLabel(1);
  fTracks[2]  = c.GetLabel(2);

  fY          = c.GetY();
  fZ          = c.GetZ();
  fSigmaY2    = c.GetSigmaY2();
  fSigmaZ2    = c.GetSigmaZ2();  

  fDetector   = c.GetDetector();
  fTimeBin    = c.GetLocalTimeBin();
  fQ          = c.GetQ();

}

//_____________________________________________________________________________
void AliTRDcluster::AddTrackIndex(Int_t *track)
{
  //
  // Adds track index. Currently assumed that track is an array of
  // size 9, and up to 3 track indexes are stored in fTracks[3].
  // Indexes are sorted according to:
  //  1) index of max number of appearances is stored first
  //  2) if two or more indexes appear equal number of times, the lowest
  //     ones are stored first;
  //

  const Int_t size = 9;

  Int_t entries[size][2], i, j, index;

  Bool_t index_added;

  for (i=0; i<size; i++) {
    entries[i][0]=-1;
    entries[i][1]=0;
  }                                 

  for (Int_t k=0; k<size; k++) {
    index=track[k];
    index_added=kFALSE; j=0;
    if (index >= 0) {
      while ( (!index_added) && ( j < size ) ) {
        if ((entries[j][0]==index) || (entries[j][1]==0)) {
          entries[j][0]=index;
          entries[j][1]=entries[j][1]+1;
          index_added=kTRUE;
        }
        j++;
      }
    }
  }             

  // sort by number of appearances and index value
  Int_t swap=1, tmp0, tmp1;
  while ( swap > 0) {
    swap=0;
    for(i=0; i<(size-1); i++) {
      if ((entries[i][0] >= 0) && (entries[i+1][0] >= 0)) {
        if ((entries[i][1] < entries[i+1][1]) ||
            ((entries[i][1] == entries[i+1][1]) &&
             (entries[i][0] > entries[i+1][0]))) {
               tmp0=entries[i][0];
               tmp1=entries[i][1];
               entries[i][0]=entries[i+1][0];
               entries[i][1]=entries[i+1][1];
               entries[i+1][0]=tmp0;
               entries[i+1][1]=tmp1;
               swap++;
        }
      }
    }
  }               

  // set track indexes
  for(i=0; i<3; i++) SetLabel(entries[i][0],i);

  return;

}          

