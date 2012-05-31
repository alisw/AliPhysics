// $Id$
// $MpId: testUpdateGlobalIndices.C,v 1.7 2005/08/24 08:53:27 ivana Exp $
//
// Tests updating global indices of motif positions from file.

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSector.h"
#include "AliMpSectorReader.h"
#include "AliMpMotifMap.h"
#include "AliMpVPainter.h"

#include <Riostream.h>
#include <TCanvas.h>

#endif

void testUpdateGlobalIndices()
{
  AliMpDataProcessor mp;
  AliMpDataMap* dataMap = mp.CreateDataMap("data");
  AliMpDataStreams dataStreams(dataMap);

  AliMpSectorReader r(dataStreams, AliMq::kStation1, AliMp::kNonBendingPlane);
  AliMpSector* sector = r.BuildSector();

  sector->GetMotifMap()->UpdateGlobalIndices("motif_map.dat");  
  
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(sector);
  TCanvas* canvas = new TCanvas();
  canvas->cd();
  painter->Draw("ZSSMI");
}          
 
