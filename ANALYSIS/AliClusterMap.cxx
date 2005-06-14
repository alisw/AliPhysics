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

//_________________________________________________
///////////////////////////////////////////////////
//
// class AliClusterMap
//
// class that describes cluster occupation at TPC
// Each padraw has a corresponding bit in fPadRawMap
// 
//
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////////////////////

#include <TString.h>

#include "AliClusterMap.h"
#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliTPCtrack.h"
#include "AliVAODParticle.h"

const Int_t AliClusterMap::fNPadRows = 159;

ClassImp(AliClusterMap)

AliClusterMap::AliClusterMap():
 fPadRawMap(fNPadRows)
{
//ctor

}
/***********************************************************************/
AliClusterMap::AliClusterMap(AliESDtrack* track):
 fPadRawMap( (track)?track->GetTPCClusterMap():fNPadRows )
{
 //ctor
 
 StdoutToAliDebug(3,Print());

} 
/***********************************************************************/

AliClusterMap::AliClusterMap(AliTPCtrack* track):
 fPadRawMap(fNPadRows)
{
 //ctor
 
 AliDebug(10,"#####################################################################"); 
 if (track == 0x0)
  {
    Error("AliClusterMap","Pointer to TPC track is NULL");
    return;
  }
 Int_t prevrow = -1;
 Int_t i = 0;
 for ( ; i < track->GetNumberOfClusters(); i++)
  {
    Int_t idx = track->GetClusterIndex(i);
    Int_t sect = (idx&0xff000000)>>24;
    Int_t row = (idx&0x00ff0000)>>16;
    if (sect > 18) row +=63; //if it is outer sector, add number of inner sectors
    AliDebug(9,Form("Cl.idx is %d, sect %d, row %d",idx,sect,row));
      
    fPadRawMap.SetBitNumber(row,kTRUE);
    
    //Fill the gap between previous row and this row with 0 bits
    if (prevrow < 0) 
     {
       prevrow = row;//if previous bit was not assigned yet == this is the first one
     }
    else
     { //we don't know the order (inner to outer or reverse)
       //just to be save in case it is going to change
       Int_t n = 0, m = 0;
       if (prevrow < row)
        {
          n = prevrow;
          m = row;
        }
       else
        {
          n = row;
          m = prevrow;
        }
       for (Int_t j = n+1; j < m; j++)
        {
          fPadRawMap.SetBitNumber(j,kFALSE);
        }
       prevrow = row; 
     }
  }
  
 StdoutToAliDebug(3,Print());

}
/***********************************************************************/

void AliClusterMap::Print() const
{
//Prints the bit map 
  TString msg;
  for ( Int_t i = 0; i < fNPadRows; i++)
   {
     if ( fPadRawMap.TestBitNumber(i) )
      {
        msg+="1";
      }
     else
      {
        msg+="0";
      }
   }
  Info("AliClusterMap","BitMap is\n  %s",msg.Data());
  
}

/***********************************************************************/

Float_t AliClusterMap::GetOverlapFactor(const AliClusterMap& clmap) const
{
  //Returns quality factor FQ = Sum(An)/Sum(clusters)
  //      | -1; if both tracks have a cluster on padrow n
  //An = <  0; if neither track has a cluster on padrow n
  //     |  1; if only one trackhas a cluster on padrow n
  // Returned value ranges between 
  //  -0.5 (low probability that these tracks are a split track)
  //  and
  //   1.0 (high probability that these tracks are a split track)
  
  Int_t nh = 0;
  Int_t an = 0;
  for ( Int_t i = 0; i < fNPadRows; i++)
   {
     Bool_t x = HasClAtPadRow(i);
     Bool_t y = clmap.HasClAtPadRow(i);
     
     if (x && y)//both have clasters
      {
       an--;
       nh+=2;
      }
     else 
      {
       
       if (x || y)//only one have cluters
        {
          an++;
          nh++;
        }
      }
   }
  
  
  Float_t retval = 0.0;
  if (nh > 0) retval = ((Float_t)an)/((Float_t)nh);
  else Warning("GetOverlapFactor","Number of counted cluters is 0.");
  
  return retval;
}
