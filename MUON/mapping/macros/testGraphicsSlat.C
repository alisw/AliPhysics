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

// $Id$
// $MpId: $

//
// Test macro for drawing slat motifs with real contours
// Christian Finck, Subatech
//
#if !defined(__CINT__) || defined(__MAKECINT__)

// MUON includes
#include "AliMpStationType.h"
#include "AliMpStation12Type.h"
#include "AliMpPlaneType.h"
#include "AliMpDataProcessor.h"
#include "AliMpDataMap.h"
#include "AliMpDataStreams.h"
#include "AliMpSt345Reader.h"
#include "AliMpPCB.h"
#include "AliMpSlat.h"
#include "AliMpVPainter.h"
#include "AliMpMotifReader.h"
#include "AliMpMotifType.h"
#include "AliMpMotifPosition.h"
#include "AliMpMotif.h"
#include "AliMpSlatMotifMap.h"
#include "AliMpDataStreams.h"

#include "TVector2.h"
#include "TCanvas.h"

#endif

AliMpDataStreams& DataStreams() {
  static AliMpDataStreams* dataStreams = nullptr;
  if (!dataStreams) {
    AliMpDataProcessor mp;
    AliMpDataMap* dataMap = mp.CreateDataMap("data");
    dataStreams = new AliMpDataStreams(dataMap);
  }
  return *dataStreams;
}

AliMpSlatMotifMap* SlatMotifMap() {
  static AliMpSlatMotifMap* motifMap = new AliMpSlatMotifMap();
  return motifMap;
}

AliMpSt345Reader* SlatReader() {
  static AliMpSt345Reader* reader = new AliMpSt345Reader(SlatMotifMap());
  return reader;
}

void testGraphicsMotif(Option_t* motifType = "R43", const char* option = "PT", Double_t padsizex=5.0, Double_t padsizey=0.5)
{
  // Warning : this function leaks memory. But should be fine as only used 
  // interactively to check a few motifs at once...
  //
  AliMp::StationType station = AliMp::kStation345;
  AliMq::Station12Type station12 = AliMq::kNotSt12;
  AliMp::PlaneType plane = AliMp::kBendingPlane;
  
  AliMpMotifReader reader(station, station12, plane);
  AliMpMotifType* type = reader.BuildMotifType(DataStreams(),motifType);
  if (!type)
  {
    cerr << "Motif not found" << endl;
    return;
  }
  type->Print("G");
  AliMpMotif* motif = new AliMpMotif(motifType,type,padsizex,padsizey);
  AliMpMotifPosition* pos = new AliMpMotifPosition(0,motif,0,0);
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(pos);
  if (!painter)
  {
    cerr << "Could not get a painter !" << endl;
    return;
  }
  new TCanvas(motifType,motifType);
  painter->Draw(option);
}

void testGraphicsMotifs(const char* option = "PT", Double_t padsizex=5.0, Double_t padsizey=0.5)
{
  std::vector <std::string> motifs = {
    "A1", "A10", "A11", "A12", "A13", "A14", "A15", "A16", "A17", "A18", "A19", "A2", "A20", "A3", "A4", "A5", "A6",
    "A7", "A8", "A9", "C1", "C10", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "E1", "E10", "E11", "E12", "E13",
    "E14", "E15", "E16", "E17", "E18", "E19", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "I1", "L1", "L10", "L11",
    "L12", "L13", "L14", "L15", "L16", "L17", "L18", "L19", "L2", "L20", "L21", "L22", "L23", "L24", "L25", "L3", "L4",
    "L5", "L6", "L7", "L8", "L9", "O1", "O10", "O11", "O12", "O13", "O14", "O15", "O16", "O18", "O19", "O2", "O20",
    "O21", "O22", "O23", "O24", "O25", "O26", "O27", "O3", "O4", "O5", "O6", "O7", "O8", "O9", "P1", "P2", "P3", "P4",
    "Q1", "Q2", "Q3", "Q4", "R1", "R10", "R11", "R12", "R13", "R14", "R15", "R16", "R17", "R18", "R19", "R2", "R20",
    "R21", "R22", "R23", "R24", "R25", "R26", "R27", "R28", "R29", "R3", "R30", "R31", "R32", "R33", "R34", "R35",
    "R36", "R37", "R38", "R39", "R4", "R40", "R41", "R42", "R43", "R44", "R45", "R5", "R6", "R7", "R8", "R9", "S0",
    "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8"
  };


  for (auto&& m:motifs) {
    if (m[0]=='E')
    testGraphicsMotif(m.c_str(), option, padsizex, padsizey);
  }
  std::cout << motifs.size() << " motifs" << std::endl;
}

void testGraphicsSlat(const char* slatType,
                      AliMp::PlaneType planeType = AliMp::kBendingPlane,
                      Option_t* option = "PMCI",
                      Bool_t saveJPG = false)
{
  // P plane
  // M motif
  // P pad
  // I indices

  TCanvas* canvas = new TCanvas(slatType,slatType,10,10,1200,800);

  AliMpSlat* slat = SlatReader()->ReadSlat(DataStreams(),slatType, planeType);
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(slat);
  painter->Draw(option);

  if (saveJPG) {
    std::string jpgname = slatType;

    if (planeType == AliMp::kNonBendingPlane)
      jpgname += "_NonBending";
    else
      jpgname += "_Bending";
    jpgname += ".jpg";
    canvas->SaveAs(jpgname.c_str());
  }
}

void testGraphicsSlats(AliMp::PlaneType planeType = AliMp::kBendingPlane,
	              Option_t* option = "PMCI",
		      Bool_t saveJPG = false)
{

  // PMPT to get manu channels numbering
  
  const char* slatName[19] = {"122000SR1", "112200SR2", "122200S", "222000N", "220000N",
			  "122000NR1", "112200NR2", "122200N",
			  "122330N", "112233NR3", "112230N", "222330N", "223300N", "333000N", "330000N",
			  "112233N", "222333N", "223330N", "333300N"};
			  

  for (Int_t i = 0; i < 19; i++) {
    testGraphicsSlat(slatName[i], planeType, option, saveJPG);
  }

}

void testGraphicsPCB(const char* pcbName="S2B",
                      Option_t* option = "MZT",
                      Bool_t savePNG = false)
{
  
  TCanvas* c = new TCanvas(pcbName,pcbName,10,10,1200,800);     
    
  AliMpPCB* pcb = SlatReader()->ReadPCB(DataStreams(),pcbName);
  if (!pcb) 
  {
    cerr << "PCB " << pcbName << " does not exist" << endl;
    return;
  }
  
  AliMpVPainter* painter = AliMpVPainter::CreatePainter(pcb);
  painter->Draw(option);
    
  if (savePNG) c->Print(Form("%s-%s.png",pcbName,option));
  
}
