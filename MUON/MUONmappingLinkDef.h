/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $MpId: basicLinkDefIn.h,v 1.8 2006/01/11 10:09:18 ivana Exp $ 

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedef;

#pragma link C++ class  AliMpArea+;
#pragma link C++ class  AliMpConstants+;
#pragma link C++ class  AliMpIntPair+;
#pragma link C++ class  AliMpExMap-;
#pragma link C++ class  AliMpArrayI+;
#pragma link C++ class  AliMpStringObjMap+;
#pragma link C++ class  AliMpPad+;
#pragma link C++ class  AliMpPadPair+;
#pragma link C++ class  AliMpVIndexed+;
#pragma link C++ class  AliMpVSegmentation+;
#pragma link C++ class  AliMpVPadIterator+;
#pragma link C++ class  AliMpPadIteratorPtr+;
#pragma link C++ class  AliMpFiles!;
#pragma link C++ class  AliMpManuList;

#pragma link C++ namespace AliMp;
#pragma link C++ enum   AliMp::Direction;
#pragma link C++ enum   AliMp::XDirection;
#pragma link C++ enum   AliMp::PlaneType;
#pragma link C++ enum   AliMp::StationType;
#pragma link C++ enum   AliMp::CathodType;

#pragma link C++ function operator-(const AliMpIntPair& ,const AliMpIntPair& );
#pragma link C++ function operator+(const AliMpIntPair& ,const AliMpIntPair& );
#pragma link C++ function operator<<(ostream& ,const AliMpIntPair& );
#pragma link C++ function operator<<(ostream& ,const AliMpPad& );
#pragma link C++ function operator<<(ostream& ,const AliMpArea& );
#pragma link C++ function operator<(const AliMpPad& ,const AliMpPad& );

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $MpId: motifLinkDefIn.h,v 1.3 2005/08/24 09:53:33 ivana Exp $ 

#pragma link C++ class  AliMpMotifReader+;
#pragma link C++ class  AliMpMotifMap+;
#pragma link C++ class  AliMpVMotif+;
#pragma link C++ class  AliMpMotif+;
#pragma link C++ class  AliMpMotifSpecial+;
#pragma link C++ class  AliMpMotifType+;
#pragma link C++ class  AliMpMotifTypePadIterator+;
#pragma link C++ class  AliMpMotifPosition+;
#pragma link C++ class  AliMpMotifPositionPadIterator+;
#pragma link C++ class  AliMpConnection+;

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $MpId: sectorLinkDefIn.h,v 1.4 2005/08/24 10:07:31 ivana Exp $ 

#pragma link C++ class  AliMpSectorReader+;
#pragma link C++ class  AliMpSector+;
#pragma link C++ class  AliMpSectorPadIterator+;
#pragma link C++ class  AliMpSectorAreaHPadIterator+;
#pragma link C++ class  AliMpSectorAreaVPadIterator+;
#pragma link C++ class  AliMpSectorSegmentation+;
#pragma link C++ class  AliMpNeighboursPadIterator+;
#pragma link C++ class  AliMpZone+;
#pragma link C++ class  AliMpSubZone+;
#pragma link C++ class  AliMpRow+;
#pragma link C++ class  AliMpVRowSegment+;
#pragma link C++ class  AliMpVRowSegmentSpecial+;
#pragma link C++ class  AliMpRowSegment+;
#pragma link C++ class  AliMpRowSegmentLSpecial+;
#pragma link C++ class  AliMpRowSegmentRSpecial+;
#pragma link C++ class  AliMpPadRow+;
#pragma link C++ class  AliMpVPadRowSegment+;
#pragma link C++ class  AliMpPadRowLSegment+;
#pragma link C++ class  AliMpPadRowRSegment+;

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $MpId: slatLinkDefIn.h,v 1.4 2005/09/19 19:01:09 ivana Exp $ 

#pragma link C++ class  AliMpHelper+;
#pragma link C++ class  AliMpSt345Reader+;
#pragma link C++ class  AliMpSlat+;
#pragma link C++ class  AliMpSlatSegmentation+;
#pragma link C++ class  AliMpPCB+;
#pragma link C++ class  AliMpSlatPadIterator+;
#pragma link C++ class  AliMpPCBPadIterator+;
#pragma link C++ class  AliMpSlatMotifMap+;

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $MpId: triggerLinkDefIn.h,v 1.1 2005/10/28 15:41:21 ivana Exp $ 

#pragma link C++ class  AliMpTrigger+;
#pragma link C++ class  AliMpTriggerReader+;
#pragma link C++ class  AliMpTriggerSegmentation+;

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $MpId: managementLinkDefIn.h,v 1.1 2006/01/11 10:24:20 ivana Exp $ 

#pragma link C++ class  AliMpSegmentation+;
#pragma link C++ class  AliMpDetElement+;
#pragma link C++ class  AliMpDEStore+;
#pragma link C++ class  AliMpDEIterator+;
#pragma link C++ class  AliMpDEManager+;
#pragma link C++ class  AliMpBusPatch+;
#pragma link C++ class  AliMpDDL+;
#pragma link C++ class  AliMpDDLStore+;

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $MpId: graphicsLinkDefIn.h,v 1.4 2005/08/24 09:42:12 ivana Exp $ 

#pragma link C++ class  AliMpGraphContext+;
#pragma link C++ class  AliMpVPainter+;
#pragma link C++ class  AliMpMotifPainter+;
#pragma link C++ class  AliMpRowPainter+;
#pragma link C++ class  AliMpRowSegmentPainter+;
#pragma link C++ class  AliMpSectorPainter+;
#pragma link C++ class  AliMpSubZonePainter+;
#pragma link C++ class  AliMpZonePainter+;
#pragma link C++ class  AliMpSlatPainter+;
#pragma link C++ class  AliMpPCBPainter+;

