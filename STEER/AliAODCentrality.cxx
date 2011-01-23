/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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


//-------------------------------------------------------------------------
//     AOD centrality class
//     Author: Alberica Toia, CERN, Alberica.Toia@cern.ch
//-------------------------------------------------------------------------

#include "AliAODCentrality.h"
#include "AliAODTrack.h"
#include "AliLog.h"

ClassImp(AliAODCentrality)


//______________________________________________________________________________

//______________________________________________________________________________
AliAODCentrality::AliAODCentrality(): TNamed(),
fxVertex          (0),
fyVertex          (0),
fzVertex          (0),
fVertexer3d       (0),
fbMC 		  (0),
fNpartTargMC	  (0),
fNpartProjMC	  (0),
fNNColl     	  (0),
fNNwColl    	  (0),
fNwNColl    	  (0),
fNwNwColl   	  (0),
fNTracklets 	  (0),
fNSingleClusters  (0),
fbZDC             (0),
fNpartZDC         (0),
fbZDCA            (0),
fNpartZDCA        (0), 
fbZDCC            (0),   
fNpartZDCC        (0),     
fESDFlag 	  (0),
fZNCEnergy	  (0),
fZPCEnergy	  (0),
fZNAEnergy	  (0),
fZPAEnergy	  (0),
fZEM1Energy	  (0),
fZEM2Energy	  (0),
fNTracks    	  (0),
fNPmdTracks 	  (0),
fMultV0A    	  (0),
fMultV0C    	  (0),
fMultFMDA    	  (0),   
fMultFMDC         (0)

{
  // constructor
  for (int i=0;i<6;i++) fNClusters[i]=0;
  for (int i=0;i<2;i++) fNChips[i]=0;
  for (int i=0;i<5;i++) fZNCtower[i]=0;
  for (int i=0;i<5;i++) fZPCtower[i]=0;
  for (int i=0;i<5;i++) fZNAtower[i]=0;
  for (int i=0;i<5;i++) fZPAtower[i]=0;
  for (int i=0;i<2;i++) fCentrZNC[i]=0;
  for (int i=0;i<2;i++) fCentrZNA[i]=0;
}

//______________________________________________________________________________
AliAODCentrality::~AliAODCentrality() 
{
  // Destructor
}

//______________________________________________________________________________
AliAODCentrality::AliAODCentrality(const AliAODCentrality& cnt) : TNamed(cnt),
fxVertex          (cnt.fxVertex   ),
fyVertex          (cnt.fyVertex   ),
fzVertex          (cnt.fzVertex   ),
fVertexer3d       (cnt.fVertexer3d),
fbMC 		  (cnt.fbMC ),
fNpartTargMC	  (cnt.fNpartTargMC),
fNpartProjMC	  (cnt.fNpartProjMC),
fNNColl     	  (cnt.fNNColl     ),
fNNwColl    	  (cnt.fNNwColl    ),
fNwNColl    	  (cnt.fNwNColl    ),
fNwNwColl   	  (cnt.fNwNwColl   ),
fNTracklets 	  (cnt.fNTracklets ),
fNSingleClusters  (cnt.fNSingleClusters),
fbZDC             (cnt.fbZDC           ),
fNpartZDC         (cnt.fNpartZDC       ),
fbZDCA            (cnt.fbZDCA          ),
fNpartZDCA        (cnt.fNpartZDCA     ), 
fbZDCC            (cnt.fbZDCC       ),   
fNpartZDCC        (cnt.fNpartZDCC ),     
fESDFlag 	  (cnt.fESDFlag ),
fZNCEnergy	  (cnt.fZNCEnergy),
fZPCEnergy	  (cnt.fZPCEnergy),
fZNAEnergy	  (cnt.fZNAEnergy),
fZPAEnergy	  (cnt.fZPAEnergy),
fZEM1Energy	  (cnt.fZEM1Energy),
fZEM2Energy	  (cnt.fZEM2Energy),
fNTracks    	  (cnt.fNTracks    ),
fNPmdTracks 	  (cnt.fNPmdTracks ),
fMultV0A    	  (cnt.fMultV0A    ),
fMultV0C    	  (cnt.fMultV0C    ),
fMultFMDA    	  (cnt.fMultFMDA ),   
fMultFMDC         (cnt.fMultFMDC )
{
  // Copy constructor.
  for (int i=0;i<6;i++) fNClusters[i] = cnt.fNClusters[i];
  for (int i=0;i<2;i++) fNChips[i]    = cnt.fNChips[i];
  for (int i=0;i<5;i++) fZNCtower[i]  = cnt.fZNCtower[i];
  for (int i=0;i<5;i++) fZPCtower[i]  = cnt.fZPCtower[i];
  for (int i=0;i<5;i++) fZNAtower[i]  = cnt.fZNAtower[i];
  for (int i=0;i<5;i++) fZPAtower[i]  = cnt.fZPAtower[i];
  for (int i=0;i<2;i++) fCentrZNC[i]  = cnt.fCentrZNC[i];
  for (int i=0;i<2;i++) fCentrZNA[i]  = cnt.fCentrZNA[i];

}

//______________________________________________________________________________
AliAODCentrality& AliAODCentrality::operator=(const AliAODCentrality& cnt) 
{
  // Assignment operator
  if (this != &cnt) {

    // name and type
    AliAODCentrality::operator=(cnt);
  }
  
  return *this;
}

//______________________________________________________________________________
void AliAODCentrality::Print(Option_t* /*option*/) const 
{
  // Print information of some data members

  printf("Centrality information:\n");
  printf("fNTracks    = %d\n",	  fNTracks    );
  printf("fNTracklets = %d\n",	  fNTracklets );
  printf("fMultV0A    = %e\n",	  fMultV0A    );
  printf("fMultV0C    = %e\n",	  fMultV0C    );
  printf("fMultFMDA   = %e\n", 	  fMultFMDA   );   
  printf("fMultFMDC   = %e\n",    fMultFMDC   );
}

