#include "AliHBTClusterMap.h"

//_________________________________________________
///////////////////////////////////////////////////
//
// class AliHBTClusterMap
//
// class that describes cluster occupation at TPC
// Each padraw has a corresponding bit in fPadRawMap
// 
//
// more info: http://aliweb.cern.ch/people/skowron/analyzer/index.html
// Piotr.Skowronski@cern.ch
//
///////////////////////////////////////////////////////////////////////////


#include "AliESDtrack.h"
#include "AliTPCtrack.h"
#include "AliHBTParticle.h"
#include <TString.h>
const Int_t AliHBTClusterMap::fNPadRows = 159;

AliHBTClusterMap::AliHBTClusterMap():
 fPadRawMap(fNPadRows)
{
//ctor

}
/***********************************************************************/
AliHBTClusterMap::AliHBTClusterMap(AliESDtrack* track):
 fPadRawMap(fNPadRows)
{
 //cotor
 track->Print();//to shut up compiler warning
 
} 
/***********************************************************************/

AliHBTClusterMap::AliHBTClusterMap(AliTPCtrack* track):
 fPadRawMap(fNPadRows)
{
 //cotor
 
 //Does not work since indeces in the claster index array 
 //in the TPC track does not correspond to the padraw segmatation
 if (AliHBTParticle::GetDebug() > 9) 
   Info("AliHBTClusterMap",
      "#####################################################################"); 
 if (track == 0x0)
  {
    Error("AliHBTClusterMap","Pointer to TPC track is NULL");
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
    if (AliHBTParticle::GetDebug() > 9)  
      Info("AliHBTClusterMap","Cl.idx is %d, sect %d, row %d",idx,sect,row);
      
    fPadRawMap.SetBitNumber(row,kTRUE);
    
    //Fill the gap between previous row and this row with 0 bits
    if (prevrow < 0) 
     {
       prevrow = row;//if previous bit was not assigned yet == this is the first one
     }
    else
     { //we don't know the order (inner to outer or reverse)
       //just to be save in case it is going change
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
  if (AliHBTParticle::GetDebug() > 2)
   {
     Print();
   } 
}
/***********************************************************************/

void AliHBTClusterMap::Print() const
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
  Info("AliHBTClusterMap","BitMap is\n  %s",msg.Data());
  
}

/***********************************************************************/

Float_t AliHBTClusterMap::GetOverlapFactor(const AliHBTClusterMap& clmap) const
{
  //Returns quality factor FQ = Sum(An)/Sum(clusters)
  //      | -1; if both tracks have a cluster on padrow n
  //An = <  0; if neither track has a cluster on padrow n
  //     |  1; if only one trackhas a cluster on padrow n
  // Returned value ranges between 
  //  -0.5 (low probability that these tracks are a split track)
  //  and
  //   1.0 (high probability that these tracks are a split track)
  TString msg1;
  TString msg2;
  
  Int_t nh = 0;
  Int_t an = 0;
  for ( Int_t i = 0; i < fNPadRows; i++)
   {
     Bool_t x = HasClAtPadRow(i);
     Bool_t y = clmap.HasClAtPadRow(i);
     
     if (x) msg1+="1";else msg1+="0";
     if (y) msg2+="1";else msg2+="0";
     
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
  
  if (AliHBTParticle::GetDebug() > 2)
   {
     Info("GetOverlapFactor","Splitting Quality Factor is %f. SumAn = %d, SumClusters %d",retval,an,nh); 
     if (retval == 1.0) 
      { 
        Print();
        Info("AliHBTClusterMap","BitMap is\n  %s\n",msg1.Data());
        clmap.Print(); 
        Info("AliHBTClusterMap","BitMap is\n  %s\n\n\n\n",msg2.Data());
      }
     if (retval == -.5) 
      { 
        Print();
        Info("AliHBTClusterMap","BitMap is\n  %s\n",msg1.Data());
        clmap.Print(); 
        Info("AliHBTClusterMap","BitMap is\n  %s\n\n\n\n",msg2.Data());
      }
   } 
 
  return retval;
}
