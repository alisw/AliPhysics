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
Revision 1.1  2000/06/15 15:47:48  morsch
Proposal for an event header class for generated events.

*/

// Event header base class for generator. 
// Stores as a minimum the date, run number, event number,
// number of particles produced  
// and the impact parameter.
// 
// Author: andreas.morsch@cern.ch

#include "AliGenEventHeader.h"
ClassImp(AliGenEventHeader)


//_____________________________________________________________________________
AliGenEventHeader::AliGenEventHeader(const char * name)
    :TNamed(name, "Event Header")
{
// Constructor
//    fDate=new TDatime;           
    fRunNumber=-1;      
    fEventNumber=-1;    
    fNProduced=-1;      
    fImpactParameter=-1.;
}




