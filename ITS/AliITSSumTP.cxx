#include "AliITSSumTP.h"
#include "AliTrackPointArray.h"

///////////////////////////////////////////////////////////////////
//                                                               //
// Class for ITS trackpoints summary + some aux. info  )         //
// Author: Ruben Shahoian                                        //
//                                                               //
///////////////////////////////////////////////////////////////////

/* $Id$ */

ClassImp(AliITSSumTP)

//__________________________________________
AliITSSumTP::AliITSSumTP(const AliITSSumTP& src) : 
		   TObject(src), fTracks(src.fTracks.GetEntriesFast()), fVertex(src.GetVertex()), 
		   fNVars(src.fNVars), fCrvVars(0)
{
  // copy c-tor
  fCrvVars = new Double32_t[fNVars];
  TObjArray& arrSrc = src.GetTracks();
  for (int i=fNVars;i--;) fCrvVars[i] = src.fCrvVars[i];
  for (int i=arrSrc.GetEntriesFast();i--;) fTracks.AddAtAndExpand(arrSrc.UncheckedAt(i),i);
}

//__________________________________________
AliITSSumTP& AliITSSumTP::operator=(const AliITSSumTP& src)
{
  // assignment op-r
  if (this == &src) return *this;
  Reset();
  TObject::operator=(src);
  fVertex = src.GetVertex();
  fNVars = src.fNVars;
  fCrvVars = new Double32_t[fNVars];
  TObjArray& arrSrc = src.GetTracks();
  for (int i=fNVars;i--;) fCrvVars[i] = src.fCrvVars[i];
  for (int i=arrSrc.GetEntriesFast();i--;) fTracks.AddAtAndExpand(arrSrc.UncheckedAt(i),i);
  return *this;
}

//__________________________________________
void AliITSSumTP::BookNTracks(Int_t n)
{
  // book space for tracks info
  delete[] fCrvVars; 
  fNVars = n*kNVarPerTrack; 
  fCrvVars = fNVars>0 ? new Double32_t[fNVars] : 0;
  for (int i=fNVars;i--;) fCrvVars[i]=0; 
}

//__________________________________________
void AliITSSumTP::Reset()
{
  // reset object
  fTracks.Delete();
  delete[] fCrvVars; 
  fCrvVars = 0;
  fNVars = 0;
  SetUniqueID(0);
}
