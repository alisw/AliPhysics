#ifdef __CINT__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
//#pragma link C++ enum   Cluster_t;

//#pragma link C++ global gITSdisplay;  // global used by AliITSdisplay

// Standard ITS classes 

#pragma link C++ class  AliITSRecPoint+;
#pragma link C++ class  AliITSclusterV2+; 
#pragma link C++ class  AliITSdigit+;
#pragma link C++ class  AliITSdigitSPD+;
#pragma link C++ class  AliITSdigitSDD+;
#pragma link C++ class  AliITSdigitSSD+;
#pragma link C++ class  AliITSTransientDigit+;
#pragma link C++ class  AliITSInitGeometry+;
#pragma link C++ class  AliITSgeom+;
#pragma link C++ class  AliITSgeomMatrix-;
#pragma link C++ class  AliITSgeomSPD+;
#pragma link C++ class  AliITSgeomSDD+;
#pragma link C++ class  AliITSgeomSSD+;
// Standard ITS detector class initilizers
#pragma link C++ class  AliITSgeomSPD300+;
#pragma link C++ class  AliITSgeomSPD425Short+;
#pragma link C++ class  AliITSgeomSPD425Long+;
#pragma link C++ class  AliITSgeomSDD256+;
#pragma link C++ class  AliITSgeomSDD300+;
#pragma link C++ class  AliITSgeomSSD175+;
#pragma link C++ class  AliITSgeomSSD275and75+;
#pragma link C++ class  AliITSgeomSSD75and275+;

#pragma link C++ class  AliITSMap+;
#pragma link C++ class  AliITSMapA1+;
#pragma link C++ class  AliITSMapA2+;
#pragma link C++ class  AliITSsegmentation+;
#pragma link C++ class  AliITSsegmentationSPD+;
#pragma link C++ class  AliITSsegmentationSDD+;
#pragma link C++ class  AliITSsegmentationSSD+;
#pragma link C++ class  AliITSresponse+;
#pragma link C++ class  AliITSCalibration+;
#pragma link C++ class  AliITSresponseSPD+;
#pragma link C++ class  AliITSresponseSDD+;
#pragma link C++ class  AliITSresponseSSD+;
#pragma link C++ class  AliITSCalibrationSPD+;
#pragma link C++ class  AliITSCalibrationSDD+;
#pragma link C++ class  AliITSCalibrationSSD+;
#pragma link C++ class  AliITSpList+;
#pragma link C++ class  AliITSpListItem+;

#pragma link C++ class  AliITSRawData+;
// These streamers must be formatted according to the raw data fromat
#pragma link C++ class  AliITSInStream+;
#pragma link C++ class  AliITSOutStream+;
//
// Raw data
#pragma link C++ class AliITSDDLRawData+;

#pragma link C++ class AliITSLoader+;

#pragma link C++ class AliITSRawStream+;
#pragma link C++ class AliITSRawStreamSDD+;
#pragma link C++ class AliITSRawStreamSDDv2+;
#pragma link C++ class AliITSRawStreamSDDv3+;
#pragma link C++ class AliITSRawStreamSPD+;
#pragma link C++ class AliITSRawStreamSSD+;
#pragma link C++ class AliITSRawStreamSSDv1+;
#pragma link C++ class AliITSEventHeader+;



#endif
