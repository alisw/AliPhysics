///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Forward Multiplicity Detector                                            //
//  This class contains the base procedures for the Forward Multiplicity     //
//  detector                                                                 //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliFMDClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:Valeri.Kondratiev@cern.ch">Valeri Kondratiev</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TTUBE.h>
#include <TNode.h>
#include "AliRun.h"
#include "AliFMD.h"
 
ClassImp(AliFMD)
 
//_____________________________________________________________________________
AliFMD::AliFMD(): AliDetector()
{
  //
  // Default constructor for class AliFMD
  //
  fIshunt   = 0;
}
 
//_____________________________________________________________________________
AliFMD::AliFMD(const char *name, const char *title)
       : AliDetector(name,title)
{
  //
  // Standard constructor for Forward Multiplicity Detector
  //
 
  //
  // Initialise Hit array
  fHits   = new TClonesArray("AliFMDhit",  405);
  
  fIshunt     =  0;
}
 
//_____________________________________________________________________________
void AliFMD::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a FMD hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliFMDhit(fIshunt,track,vol,hits);
}
 
//_____________________________________________________________________________
void AliFMD::BuildGeometry()
{
  //
  // Build simple ROOT TNode geometry for event display
  //
  TNode *Node, *Top;
  const int kColorFMD  = 7;
  //
  Top=gAlice->GetGeometry()->GetNode("alice");

  // FMD define the different volumes
  
  new TTUBE("S_FMD1","FMD sensitive volume 1","void",4.5,10.5,3);
  Top->cd();
  Node = new TNode("FMD1","FMD1","S_FMD1",0,0,86.5,"");
  Node->SetLineColor(kColorFMD);
  fNodes->Add(Node);
  
  new TTUBE("S_FMD2","FMD sensitive volume 2","void",4.5,10.5,3);
  Top->cd();
  Node = new TNode("FMD2","FMD2","S_FMD2",0,0,-86.5,"");
  Node->SetLineColor(kColorFMD);
  fNodes->Add(Node);
  
  new TTUBE("S_FMD3","FMD sensitive volume 3","void",8,14,3);
  Top->cd();
  Node = new TNode("FMD3","FMD3","S_FMD3",0,0,71.2,"");
  Node->SetLineColor(kColorFMD);
  fNodes->Add(Node);
  
  new TTUBE("S_FMD4","FMD sensitive volume 4","void",8,14,3);
  Top->cd();
  Node = new TNode("FMD4","FMD4","S_FMD4",0,0,-71.2,"");
  fNodes->Add(Node);
  Node->SetLineColor(kColorFMD);
  
  new TTUBE("S_FMD5","FMD sensitive volume 5","void",8,17.5,3);
  Top->cd();
  Node = new TNode("FMD5","FMD5","S_FMD5",0,0,44,"");
  Node->SetLineColor(kColorFMD);
  fNodes->Add(Node);
  
  new TTUBE("S_FMD6","FMD sensitive volume 6","void",8,17.5,3);
  Top->cd();
  Node = new TNode("FMD6","FMD6","S_FMD6",0,0,-44,"");
  Node->SetLineColor(kColorFMD);
  fNodes->Add(Node);
  
  new TTUBE("S_FMD7","FMD sensitive volume 7","void",4.2,13,3);
  Top->cd();
  Node = new TNode("FMD7","FMD7","S_FMD7",0,0,-231,"");
  Node->SetLineColor(kColorFMD);
  fNodes->Add(Node);
}
 
//_____________________________________________________________________________
Int_t AliFMD::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Calculate the distance from the mouse to the FMD on the screen
  // Dummy routine
  //
  return 9999;
}
 
//_____________________________________________________________________________
void AliFMD::StepManager()
{
  //
  // Called for each step in the FMD
  //
  
  Float_t       hits[3];
  Int_t         i,copy,vol[1];
  TClonesArray &lhits = *fHits;
  TLorentzVector p;
  
  gMC->CurrentVolID(copy);
  vol[0] = copy;
  gMC->TrackPosition(p);
  for(i=0;i<3;++i) hits[i]=p[i];
  new(lhits[fNhits++]) AliFMDhit(fIshunt,gAlice->CurrentTrack(),vol,hits);
}

//___________________________________________
void AliFMD::Init()
{
  //
  // Initialis the FMD after it has been built
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" FMD_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  // Here the FMD initialisation code (if any!)
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

 
ClassImp(AliFMDhit)
 
//_____________________________________________________________________________
AliFMDhit::AliFMDhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
  AliHit(shunt, track)
{
  //
  // Add a FMD hit
  //
  fVolume = vol[0];
  fX=hits[0];
  fY=hits[1];
  fZ=hits[2];
}
 
