
#include "flipPCB.h"

#include "AliMpMotif.h"
#include "AliMpMotifType.h"
#include "AliMpMotifPosition.h"
#include "AliMpExMap.h"
#include "AliMpPCB.h"
#include "AliMpConnection.h"
#include "AliMpSt345Reader.h"
#include "Riostream.h"
#include "TCanvas.h"
#include "AliMpVPainter.h"
#include "AliMpConnection.h"
#include "AliMpSlatMotifMap.h"

const char* NameIt(const TString& baseName)
{
  return Form("%sR",baseName.Data());
}

//______________________________________________________________________________
AliMpPCB* Duplicate(const AliMpPCB& src, AliMpSlatMotifMap& motifMap)
{
  AliMpPCB* dest = new AliMpPCB(&motifMap,NameIt(src.GetID()),
                                src.PadSizeX(),src.PadSizeY(),
                                src.DX()*2,src.DY()*2);
  
  for ( Int_t i = 0; i < src.GetSize(); ++i )
  {
    AliMpMotifPosition* srcMotifPosition = src.GetMotifPosition(i);
    AliMpMotifType* motifType = srcMotifPosition->GetMotif()->GetMotifType();
    AliMpMotifType* mt = motifMap.FindMotifType(motifType->GetID());
    if (!mt)
    {
      mt = static_cast<AliMpMotifType*>(motifType->Clone());
      motifMap.AddMotifType(mt);
    }
    dest->Add(mt,
              srcMotifPosition->GetLowIndicesLimit().GetFirst(),
              srcMotifPosition->GetLowIndicesLimit().GetSecond());
  }
  return dest;
}

//______________________________________________________________________________
AliMpConnection* getConnection(const AliMpPCB& pcb, Int_t ix, Int_t iy)
{
  AliMpMotifPosition* pos = pcb.FindMotifPosition(ix,iy);
  if (pos)
  {
    AliMpMotifType* type = pos->GetMotif()->GetMotifType();
    return type->FindConnectionByLocalIndices(AliMpIntPair(ix,iy)-pos->GetLowIndicesLimit());
  }
  return 0x0;
}
                                                      
//______________________________________________________________________________
AliMpPCB* flipX(const AliMpPCB& src,
                AliMpSlatMotifMap& destMotifs)
{
  AliMpPCB* dest = Duplicate(src,destMotifs);

  for ( Int_t iy = dest->Iymin(); iy <= dest->Iymax(); ++iy )
  {
    for ( Int_t ix = dest->Ixmin(); ix <= dest->Ixmax(); ++ix )
    {
      AliMpConnection* srcConnection = getConnection(src,ix,iy);
      Int_t dix = dest->Ixmax() - ix;
      AliMpConnection* destConnection = getConnection(*dest,dix,iy);
      if ( srcConnection && destConnection )
      {
        destConnection->SetGassiNum(srcConnection->GetGassiNum());
      }
      else
      {
        if ( srcConnection && !destConnection )
        {
          cout << Form("Got no connection in dest for ix,iy=%d,%d src "
                       "ix,iy=%d,%d in PCB %s",
                       dix,iy,ix,iy,dest->GetID()) << endl;
          dest->Print("M");
          cout << "---------- src=" << endl;
          src.Print("M");
          delete dest;
          return 0x0;
        }
      }
    }
  }
  return dest;
}

//______________________________________________________________________________
void flipPCB(const char* srcName)
{
  // flip PCB in X-direction

  AliMpSlatMotifMap* srcMotifs = new AliMpSlatMotifMap;
  AliMpSlatMotifMap* destMotifs = new AliMpSlatMotifMap;
  
  AliMpSt345Reader reader(*srcMotifs);

  AliMpPCB* src = reader.ReadPCB(srcName);
  
  if (!src) return;

  AliMpPCB* dest = flipX(*src,*destMotifs);
    
  if (!dest)
  {
    cout << "Flipping failed" << endl;
  }
  
  TCanvas* csrc = new TCanvas("src","src");
  csrc->Draw();
  AliMpVPainter* psrc = AliMpVPainter::CreatePainter(src);
  psrc->Draw("MZT");

  if ( dest ) 
  {
    TCanvas* cdest = new TCanvas("dest","dest");
    cdest->Draw();
    AliMpVPainter* pdest = AliMpVPainter::CreatePainter(dest);
    pdest->Draw("MZT");
    dest->Save();

  }
  
}