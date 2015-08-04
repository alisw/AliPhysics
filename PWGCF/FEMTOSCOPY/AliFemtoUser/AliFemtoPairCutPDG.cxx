/////////////////////////////////////////////////////////////////////////////////////
//                                                                                 //
// AliFemtoPairCutPDG - a pair cut which checks if the sum of the transverse       //
// momenta of two particles fit within given range and provides cut for the        //
// pair of particles originated from the same mother particle                     //
// Authors: Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch     //
//          Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch //
//          Piotr Modzelewski pmodzele@cern.ch                                     //
//  	       				                                                           //
/////////////////////////////////////////////////////////////////////////////////////

#include "AliFemtoPairCutPDG.h"
#include <string>
#include <cstdio>
#include <TMath.h>
#include "TChain.h"
#include "TTree.h"

#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"

#include "AliStack.h"

#ifdef __ROOT__
ClassImp(AliFemtoPairCutPDG)
#endif

#define PIH 1.57079632679489656
#define PIT 6.28318530717958623

//__________________
AliFemtoPairCutPDG::AliFemtoPairCutPDG():
  AliFemtoPairCut(),
  fSumPtMin(0),
  fSumPtMax(10000),
  fPDG1(0),
  fPDG2(0),
  fNPairsFailed(0),
  fNPairsPassed(0)
{

}
//__________________
AliFemtoPairCutPDG::AliFemtoPairCutPDG(double lo, double hi, int pdg_1, int pdg_2):
  AliFemtoPairCut(),
  fSumPtMin(lo),
  fSumPtMax(hi),
  fPDG1(pdg_1),
  fPDG2(pdg_2),
  fNPairsFailed(0),
  fNPairsPassed(0)
{

}
//__________________
AliFemtoPairCutPDG::AliFemtoPairCutPDG(const AliFemtoPairCutPDG& c) : 
  AliFemtoPairCut(c),
  fSumPtMin(0),
  fSumPtMax(0),
  fPDG1(0),
  fPDG2(0),
  fNPairsFailed(0),
  fNPairsPassed(0)
{ 
  fSumPtMin = c.fSumPtMin;
  fSumPtMax = c.fSumPtMax;
  fPDG1 = c.fPDG1;
  fPDG2 = c.fPDG2;
}
AliFemtoPairCutPDG& AliFemtoPairCutPDG::operator=(const AliFemtoPairCutPDG& c)
{
  if (this != &c) {
    fSumPtMin = c.fSumPtMin;
    fSumPtMax = c.fSumPtMax;
    fPDG1 = c.fPDG1;
    fPDG2 = c.fPDG2;
  }

  return *this;

}

//__________________
AliFemtoPairCutPDG::~AliFemtoPairCutPDG(){
  /* no-op */
}
//__________________
bool AliFemtoPairCutPDG::Pass(const AliFemtoPair* pair){

  bool temp = true;

  double pt1 = pair->Track1()->Track()->Pt();
  double pt2 = pair->Track2()->Track()->Pt();
  double pt_sum = pt1 + pt2;

  if(pt_sum >= fSumPtMin && pt_sum <= fSumPtMax) temp = true;
  else temp = false;

  //find particles which originate from the same mother particle and eliminate them
  if (temp)
  {
    AliMCEventHandler* mctruth = (AliMCEventHandler*)((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
    AliStack* fStack = mctruth->MCEvent()->Stack();

    //starting to reading tracks
    int nofTracks=0;  //number of all tracks in MC event
    nofTracks=fStack->GetNtrack();
    for(int i=0;i<nofTracks;i++)
    {
      if(!fStack->IsPhysicalPrimary(i))
      {
        continue;
      }
      TParticle *kinetrack1 = fStack->Particle(i);
      for(int j=i+1;j<nofTracks;j++)
      {
        //take only primaries
        if(!fStack->IsPhysicalPrimary(j))
        {
          continue;
        }

        //getting next track
        TParticle *kinetrack2 = fStack->Particle(j);

        int pdgCode1 = fPDG1;
        int pdgCode2 = fPDG2;

        if((kinetrack1->GetPdgCode()==pdgCode1 && kinetrack2->GetPdgCode()==pdgCode2) || (kinetrack1->GetPdgCode()==pdgCode2 && kinetrack2->GetPdgCode()==pdgCode1)) 
        {
	       int mother1 = kinetrack1->GetFirstMother();
	       int mother2 = kinetrack2->GetFirstMother();

          if(mother1==mother2)
	       {
            temp = false;
            cout<<"PDG cut -> Pair rejected!"<<endl;
	       }
        }
      }
    } 
  }
   if(temp) fNPairsPassed++;
   else fNPairsFailed++;
   return temp;
}

//__________________
AliFemtoString AliFemtoPairCutPDG::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoPairCutPDG Pair Cut\n";  
  char ctemp[100];
  stemp += ctemp;
  snprintf(ctemp,100,"Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",(long int) fNPairsPassed,(long int) fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoPairCutPDG::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings = new TList();
  char buf[200];
  snprintf(buf, 200, "AliFemtoPairCutPDG.sumptmin=%f", fSumPtMin);
  snprintf(buf, 200, "AliFemtoPairCutPr.sumptmax=%f", fSumPtMax);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}

void AliFemtoPairCutPDG::SetMinSumPt(Double_t sumptmin)
{
  fSumPtMin = sumptmin;
}

void AliFemtoPairCutPDG::SetMaxSumPt(Double_t sumptmax)
{
  fSumPtMax = sumptmax;
}

void AliFemtoPairCutPDG::SetPDG1(Double_t pdg1)
{
  fPDG1 = pdg1;
}

void AliFemtoPairCutPDG::SetPDG2(Double_t pdg2)
{
  fPDG2 = pdg2;
}

