#include "AliHBTAnalysis.h"
//_________________________________________________________
////////////////////////////////////////////////////////////////////////////
//
// class AliHBTAnalysis
//
// Central Object Of HBTAnalyser: 
// This class performs main looping within HBT Analysis
// User must plug a reader of Type AliReader
// User plugs in coorelation and monitor functions
// as well as monitor functions
//
// HBT Analysis Tool, which is integral part of AliRoot,
// ALICE Off-Line framework:
//
// Piotr.Skowronski@cern.ch
// more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
//
////////////////////////////////////////////////////////////////////////////
//_________________________________________________________


#include <TSystem.h>
#include <TFile.h>

#include "AliAOD.h"
#include "AliAODParticle.h"
#include "AliAODPairCut.h"
#include "AliEventCut.h"

#include "AliEventBuffer.h"

#include "AliReader.h"
#include "AliHBTPair.h"
#include "AliHBTFunction.h"
#include "AliHBTMonitorFunction.h"
 

ClassImp(AliHBTAnalysis)

const UInt_t AliHBTAnalysis::fgkFctnArraySize = 100;
const UInt_t AliHBTAnalysis::fgkDefaultMixingInfo = 1000;
const Int_t  AliHBTAnalysis::fgkDefaultBufferSize = 5;

AliHBTAnalysis::AliHBTAnalysis():
  fProcEvent(0x0),
  fReader(0x0),
  fNTrackFunctions(0),
  fNParticleFunctions(0),
  fNParticleAndTrackFunctions(0),
  fNTrackMonitorFunctions(0),
  fNParticleMonitorFunctions(0), 
  fNParticleAndTrackMonitorFunctions(0),
  fTrackFunctions ( new AliHBTOnePairFctn* [fgkFctnArraySize]),
  fParticleFunctions ( new AliHBTOnePairFctn* [fgkFctnArraySize]),
  fParticleAndTrackFunctions ( new AliHBTTwoPairFctn* [fgkFctnArraySize]),
  fParticleMonitorFunctions ( new AliHBTMonOneParticleFctn* [fgkFctnArraySize]),    
  fTrackMonitorFunctions ( new AliHBTMonOneParticleFctn* [fgkFctnArraySize]),    
  fParticleAndTrackMonitorFunctions ( new AliHBTMonTwoParticleFctn* [fgkFctnArraySize]),    
  fBkgEventCut(0x0),
  fPartBuffer(0x0),
  fTrackBuffer(0x0),
  fBufferSize(2),
  fDisplayMixingInfo(fgkDefaultMixingInfo),
  fIsOwner(kFALSE),
  fProcessOption(kSimulatedAndReconstructed),
  fNoCorrfctns(kFALSE),
  fOutputFileName(0x0),
  fVertexX(0.0),
  fVertexY(0.0),
  fVertexZ(0.0)
 {
   //default constructor
   
 }
/*************************************************************************************/ 

AliHBTAnalysis::AliHBTAnalysis(const AliHBTAnalysis& in):
  AliAnalysis(in),
  fProcEvent(0x0),
  fReader(0x0),
  fNTrackFunctions(0),
  fNParticleFunctions(0),
  fNParticleAndTrackFunctions(0),
  fNTrackMonitorFunctions(0),
  fNParticleMonitorFunctions(0), 
  fNParticleAndTrackMonitorFunctions(0),
  fTrackFunctions(0x0),
  fParticleFunctions(0x0),
  fParticleAndTrackFunctions(0x0),
  fParticleMonitorFunctions(0x0),
  fTrackMonitorFunctions(0x0),
  fParticleAndTrackMonitorFunctions(0x0),
  fBkgEventCut(0x0),
  fPartBuffer(0x0),
  fTrackBuffer(0x0),
  fBufferSize(fgkDefaultBufferSize),
  fDisplayMixingInfo(fgkDefaultMixingInfo),
  fIsOwner(kFALSE),
  fProcessOption(kSimulatedAndReconstructed),
  fNoCorrfctns(kFALSE),
  fOutputFileName(0x0),
  fVertexX(0.0),
  fVertexY(0.0),
  fVertexZ(0.0)
 {
//copy constructor
   Fatal("AliHBTAnalysis(const AliHBTAnalysis&)","Sensless");
 }
/*************************************************************************************/ 
AliHBTAnalysis& AliHBTAnalysis::operator=(const AliHBTAnalysis& /*right*/)
 {
//operator =
   Fatal("AliHBTAnalysis(const AliHBTAnalysis&)","Sensless");
   return *this;
 }
/*************************************************************************************/ 
AliHBTAnalysis::~AliHBTAnalysis()
 {
 //destructor
 //note that we do not delete functions itself
 // they should be deleted by whom where created
 //we only store pointers, and we use "new" only for pointers array
    
   if (fIsOwner)
    {
      if (AliVAODParticle::GetDebug()>5)Info("~AliHBTAnalysis","Is Owner: Attempting to delete functions");
      DeleteFunctions();
      if (AliVAODParticle::GetDebug()>5)Info("~AliHBTAnalysis","Delete functions done");
    }
   delete [] fTrackFunctions;
   delete [] fParticleFunctions;
   delete [] fParticleAndTrackFunctions;
   
   delete [] fParticleMonitorFunctions; 
   delete [] fTrackMonitorFunctions; 
   delete [] fParticleAndTrackMonitorFunctions;

   delete fBkgEventCut;
   delete fOutputFileName;
 }

/*************************************************************************************/ 

Int_t AliHBTAnalysis::ProcessEvent(AliAOD* aodrec, AliAOD* aodsim)
{
  //Processes one event
  if (fProcEvent == 0x0)
   {
     Error("ProcessEvent","Analysis <<%s>> option not specified.",GetName());
     return 1;
   }
  if ( Pass(aodrec,aodsim) ) return 0;
  
  //Move event to the apparent vertex -> must be after the event cut
  //It is important for any cut that use any spacial coordiantes, 
  //f.g. anti-merging cut in order to preserve the same bahavior of variables (f.g. distance between tracks)
  Double_t dvx = 0, dvy = 0, dvz = 0;
  if (aodrec)
   {
     Double_t pvx,pvy,pvz;
     aodrec->GetPrimaryVertex(pvx,pvy,pvz);
     
     dvx = fVertexX - pvx;
     dvy = fVertexY - pvy;
     dvz = fVertexZ - pvz;
     aodrec->Move(dvx,dvy,dvz);
   }  
  
  Int_t result = (this->*fProcEvent)(aodrec,aodsim);

  if (aodrec) aodrec->Move(-dvx,-dvy,-dvz);//move event back to the same spacial coordinates
  
  return  result;
}
/*************************************************************************************/ 

Int_t AliHBTAnalysis::Finish()
{
//Finishes analysis
  WriteFunctions();
  return 0;
}
/*************************************************************************************/ 

void AliHBTAnalysis::DeleteFunctions()
{
 //Deletes all functions added to analysis
 
 UInt_t ii;
 for(ii = 0;ii<fNParticleFunctions;ii++)
  { 
    if (AliVAODParticle::GetDebug()>5)
     {
       Info("DeleteFunctions","Deleting ParticleFunction %#x",fParticleFunctions[ii]);
       Info("DeleteFunctions","Deleting ParticleFunction %s",fParticleFunctions[ii]->Name());
     }
    delete fParticleFunctions[ii];
  } 
 fNParticleFunctions = 0;
                
 for(ii = 0;ii<fNTrackFunctions;ii++)
  { 
    if (AliVAODParticle::GetDebug()>5)
     {
       Info("DeleteFunctions","Deleting TrackFunction %#x",fTrackFunctions[ii]);
       Info("DeleteFunctions","Deleting TrackFunction %s",fTrackFunctions[ii]->Name());
     }
    delete fTrackFunctions[ii];
  }  
 fNTrackFunctions = 0;
 
 for(ii = 0;ii<fNParticleAndTrackFunctions;ii++)
  { 
    if (AliVAODParticle::GetDebug()>5)
     {
       Info("DeleteFunctions","Deleting ParticleAndTrackFunction %#x",fParticleAndTrackFunctions[ii]);
       Info("DeleteFunctions","Deleting ParticleAndTrackFunction %s",fParticleAndTrackFunctions[ii]->Name());
     } 
    delete fParticleAndTrackFunctions[ii];
  }  
 fNParticleAndTrackFunctions = 0;
 
 for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
  { 
    if (AliVAODParticle::GetDebug()>5)
     {
       Info("DeleteFunctions","Deleting ParticleMonitorFunction %#x",fParticleMonitorFunctions[ii]);
       Info("DeleteFunctions","Deleting ParticleMonitorFunction %s",fParticleMonitorFunctions[ii]->Name());
     }
    delete fParticleMonitorFunctions[ii];
  } 
 fNParticleMonitorFunctions = 0;
   
 for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
  { 
    if (AliVAODParticle::GetDebug()>5)
     {
       Info("DeleteFunctions","Deleting TrackMonitorFunction %#x",fTrackMonitorFunctions[ii]);
       Info("DeleteFunctions","Deleting TrackMonitorFunction %s",fTrackMonitorFunctions[ii]->Name());
     }
    delete fTrackMonitorFunctions[ii];
  } 
 fNTrackMonitorFunctions = 0;
   
 for(ii = 0; ii<fNParticleAndTrackMonitorFunctions; ii++)
  { 
    if (AliVAODParticle::GetDebug()>5)
     {
      Info("DeleteFunctions","Deleting ParticleAndTrackMonitorFunction %#x",fParticleAndTrackMonitorFunctions[ii]);
      Info("DeleteFunctions","Deleting ParticleAndTrackMonitorFunction %s",fParticleAndTrackMonitorFunctions[ii]->Name());
     }
    delete fParticleAndTrackMonitorFunctions[ii];
  } 
 fNParticleAndTrackMonitorFunctions = 0;
 
}
/*************************************************************************************/ 

Int_t AliHBTAnalysis::Init()
{
//Initializeation method
//calls Init for all functions

 //Depending on option and pair cut assigns proper analysis method
 Bool_t nonid = IsNonIdentAnalysis();
 switch(fProcessOption)
  {
    case kReconstructed:
      if (nonid) fProcEvent = &AliHBTAnalysis::ProcessRecNonId;
      else fProcEvent = &AliHBTAnalysis::ProcessRec;
      SetCutsOnRec();
      break;

    case kSimulated:
      if (nonid) fProcEvent = &AliHBTAnalysis::ProcessSimNonId;
      else fProcEvent = &AliHBTAnalysis::ProcessSim;
      SetCutsOnSim();
      break;

    case kSimulatedAndReconstructed:
      if (nonid) fProcEvent = &AliHBTAnalysis::ProcessRecAndSimNonId;
      else fProcEvent = &AliHBTAnalysis::ProcessRecAndSim;
      break;
  }
  
 if (fPartBuffer == 0x0)  fPartBuffer = new AliEventBuffer (fBufferSize);
 else fPartBuffer->Reset();

 if (fTrackBuffer == 0x0)  fTrackBuffer = new AliEventBuffer (fBufferSize);
 else fTrackBuffer->Reset();


 fNoCorrfctns = (fNParticleFunctions == 0) && (fNTrackFunctions == 0) && (fNParticleAndTrackFunctions == 0);
 
 UInt_t ii;
 for(ii = 0;ii<fNParticleFunctions;ii++)
   fParticleFunctions[ii]->Init();
                
 for(ii = 0;ii<fNTrackFunctions;ii++)
   fTrackFunctions[ii]->Init();

 for(ii = 0;ii<fNParticleAndTrackFunctions;ii++)
   fParticleAndTrackFunctions[ii]->Init();
 
 for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
   fParticleMonitorFunctions[ii]->Init();
   
 for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
   fTrackMonitorFunctions[ii]->Init();
   
 for(ii = 0; ii<fNParticleAndTrackMonitorFunctions; ii++)
   fParticleAndTrackMonitorFunctions[ii]->Init();
   
 return 0;
}
/*************************************************************************************/ 

void AliHBTAnalysis::ResetFunctions()
{
//In case fOwner is true, deletes all functions
//in other case, just set number of analysis to 0
 if (fIsOwner) DeleteFunctions();
 else
  {
    fNParticleFunctions = 0;
    fNTrackFunctions = 0;
    fNParticleAndTrackFunctions = 0;
    fNParticleMonitorFunctions = 0;
    fNTrackMonitorFunctions = 0;
    fNParticleAndTrackMonitorFunctions = 0;
  }
}
/*************************************************************************************/ 
Int_t AliHBTAnalysis::ProcessRecAndSim(AliAOD* aodrec, AliAOD* aodsim)
{
//Does analysis for both tracks and particles
//mainly for resolution study and analysies with weighting algirithms
  
// cut on particles only -- why?
// - PID: when we make resolution analysis we want to take only tracks with correct PID
// We need cut on tracks because there are data characteristic
  
  AliVAODParticle * part1, * part2;
  AliVAODParticle * track1, * track2;
  
  AliAOD * trackEvent = aodrec, *partEvent = aodsim;
  AliAOD* trackEvent1 = new AliAOD();
  AliAOD* partEvent1 = new AliAOD();
  partEvent1->SetOwner(kTRUE);
  trackEvent1->SetOwner(kTRUE);

  AliAOD * trackEvent2,*partEvent2;
  
//  Int_t N1, N2, N=0; //number of particles in current event(we prcess two events in one time)
  
//  Int_t nev = fReader->GetNumberOfTrackEvents();
  static AliHBTPair tpair;
  static AliHBTPair ppair;

  AliHBTPair* trackpair = &tpair;
  AliHBTPair* partpair = &ppair;
 
  AliHBTPair * tmptrackpair;//temprary pointers to pairs
  AliHBTPair * tmppartpair;
  
  register UInt_t ii;
  
  

  if ( !partEvent || !trackEvent ) 
   {
     Error("ProcessRecAndSim","<<%s>> Can not get event",GetName());
     return 1;
   }

  if ( partEvent->GetNumberOfParticles() != trackEvent->GetNumberOfParticles() )
   {
     Error("ProcessRecAndSim",
           "Number of simulated particles (%d) not equal to number of reconstructed tracks (%d). Skipping Event.",
            partEvent->GetNumberOfParticles() , trackEvent->GetNumberOfParticles());
     return 2;
   }


  for (Int_t j = 0; j<partEvent->GetNumberOfParticles() ; j++)
   {
     /***************************************/
     /******   Looping same events   ********/
     /******   filling numerators    ********/
     /***************************************/
     if ( (j%fDisplayMixingInfo) == 0)
        Info("ProcessTracksAndParticles",
             "Mixing particle %d with particles from the same event",j);

     part1= partEvent->GetParticle(j);
     track1= trackEvent->GetParticle(j);

     Bool_t firstcut = (this->*fkPass1)(part1,track1);
     if (fBufferSize != 0) 
       if ( (firstcut == kFALSE) || ( (this->*fkPass2)(part1,track1) == kFALSE ) )
        {
          //accepted by any cut
          // we have to copy because reader keeps only one event

          partEvent1->AddParticle(new AliAODParticle(*part1));
          trackEvent1->AddParticle(new AliAODParticle(*track1));
        }

     if (firstcut) continue;

     for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
       fParticleMonitorFunctions[ii]->Process(part1);
     for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
       fTrackMonitorFunctions[ii]->Process(track1);
     for(ii = 0; ii<fNParticleAndTrackMonitorFunctions; ii++)
       fParticleAndTrackMonitorFunctions[ii]->Process(track1,part1);

     if (fNoCorrfctns) continue;

     for (Int_t k =j+1; k < partEvent->GetNumberOfParticles() ; k++)
      {
        part2= partEvent->GetParticle(k);
        if (part1->GetUID() == part2->GetUID()) continue;
        partpair->SetParticles(part1,part2);

        track2= trackEvent->GetParticle(k);
        trackpair->SetParticles(track1,track2);

        if( (this->*fkPass)(partpair,trackpair) ) //check pair cut 
          { //do not meets crietria of the pair cut, try with swapped pairs
            if( (this->*fkPass)((AliHBTPair*)partpair->GetSwappedPair(),(AliHBTPair*)trackpair->GetSwappedPair()) )
              continue; //swaped pairs do not meet criteria of pair cut as well, take next particle
            else 
             { //swaped pair meets all the criteria
               tmppartpair = (AliHBTPair*)partpair->GetSwappedPair();
               tmptrackpair = (AliHBTPair*)trackpair->GetSwappedPair();
             }
          }
        else
         {//meets criteria of the pair cut
           tmptrackpair = trackpair;
           tmppartpair = partpair;
         }

        for(ii = 0;ii<fNParticleFunctions;ii++)
               fParticleFunctions[ii]->ProcessSameEventParticles(tmppartpair);

        for(ii = 0;ii<fNTrackFunctions;ii++)
               fTrackFunctions[ii]->ProcessSameEventParticles(tmptrackpair);

        for(ii = 0;ii<fNParticleAndTrackFunctions;ii++)
               fParticleAndTrackFunctions[ii]->ProcessSameEventParticles(tmptrackpair,tmppartpair);
       //end of 2nd loop over particles from the same event  
      }//for (Int_t k =j+1; k < partEvent->GetNumberOfParticles() ; k++)

     /***************************************/
     /***** Filling denominators    *********/
     /***************************************/
     if (fBufferSize == 0) continue;

     fPartBuffer->ResetIter();
     fTrackBuffer->ResetIter();
     Int_t m = 0;
     while (( partEvent2 = fPartBuffer->Next() ))
      {
        trackEvent2 = fTrackBuffer->Next();

        m++;
        if ( (j%fDisplayMixingInfo) == 0)
           Info("ProcessTracksAndParticles",
                "Mixing particle %d from current event with particles from event %d",j,-m);

        for(Int_t l = 0; l<partEvent2->GetNumberOfParticles();l++)   //  ... on all particles
          {
            part2= partEvent2->GetParticle(l);
            partpair->SetParticles(part1,part2);

            track2= trackEvent2->GetParticle(l);
            trackpair->SetParticles(track1,track2);

            if( (this->*fkPass)(partpair,trackpair) ) //check pair cut
              { //do not meets crietria of the 
                if( (this->*fkPass)((AliHBTPair*)partpair->GetSwappedPair(),(AliHBTPair*)trackpair->GetSwappedPair()) )
                  continue;
                else 
                 {
                   tmppartpair = (AliHBTPair*)partpair->GetSwappedPair();
                   tmptrackpair = (AliHBTPair*)trackpair->GetSwappedPair();
                 }
              }
            else
             {//meets criteria of the pair cut
              tmptrackpair = trackpair;
              tmppartpair = partpair;
             }

            for(ii = 0;ii<fNParticleFunctions;ii++)
              fParticleFunctions[ii]->ProcessDiffEventParticles(tmppartpair);

            for(ii = 0;ii<fNTrackFunctions;ii++)
              fTrackFunctions[ii]->ProcessDiffEventParticles(tmptrackpair);

            for(ii = 0;ii<fNParticleAndTrackFunctions;ii++)
              fParticleAndTrackFunctions[ii]->ProcessDiffEventParticles(tmptrackpair,tmppartpair);
          }//for(Int_t l = 0; l<N2;l++)   //  ... on all particles

      }
    //end of loop over particles from first event
   }//for (Int_t j = 0; j<partEvent->GetNumberOfParticles() ; j++)
  delete fPartBuffer->Push(partEvent1);
  delete fTrackBuffer->Push(trackEvent1);
 //end of loop over events  
  return 0;
}
/*************************************************************************************/ 

Int_t AliHBTAnalysis::ProcessSim(AliAOD* /*aodrec*/, AliAOD* aodsim)
{
  //Does analysis of simulated data
  AliVAODParticle * part1, * part2;
  
  AliAOD* partEvent = aodsim;
  AliAOD* partEvent1 = new AliAOD();
  partEvent1->SetOwner(kTRUE);

  AliAOD* partEvent2;
  
  AliHBTPair ppair;

  AliHBTPair* partpair = &ppair;
 
  AliHBTPair * tmppartpair;
  
  register UInt_t ii;
  
  

  if ( !partEvent )
   {
     Error("ProcessRecAndSim","Can not get event");
     return 1;
   }


  for (Int_t j = 0; j<partEvent->GetNumberOfParticles() ; j++)
   {
     /***************************************/
     /******   Looping same events   ********/
     /******   filling numerators    ********/
     /***************************************/
     if ( (j%fDisplayMixingInfo) == 0)
        Info("ProcessTracksAndParticles",
             "Mixing particle %d with particles from the same event",j);

     part1= partEvent->GetParticle(j);

     Bool_t firstcut = fPairCut->GetFirstPartCut()->Pass(part1);

     if (fBufferSize != 0) 
       if ( (firstcut == kFALSE) || ( fPairCut->GetSecondPartCut()->Pass(part1) == kFALSE ) )
        {
          //accepted by any cut
          // we have to copy because reader keeps only one event

          partEvent1->AddParticle(new AliAODParticle(*part1));
        }

     if (firstcut) continue;

     for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
       fParticleMonitorFunctions[ii]->Process(part1);

     if ( fNParticleFunctions == 0 ) continue;

     for (Int_t k =j+1; k < partEvent->GetNumberOfParticles() ; k++)
      {
        part2= partEvent->GetParticle(k);
        if (part1->GetUID() == part2->GetUID()) continue;
        partpair->SetParticles(part1,part2);

           if(fPairCut->Pass(partpair)) //check pair cut
            { //do not meets crietria of the 
              if( fPairCut->Pass((AliHBTPair*)partpair->GetSwappedPair()) ) continue;
              else tmppartpair = (AliHBTPair*)partpair->GetSwappedPair();
            }
           else
            {
              tmppartpair = partpair;
            }

        for(ii = 0;ii<fNParticleFunctions;ii++)
               fParticleFunctions[ii]->ProcessSameEventParticles(tmppartpair);

       //end of 2nd loop over particles from the same event  
      }//for (Int_t k =j+1; k < partEvent->GetNumberOfParticles() ; k++)

     /***************************************/
     /***** Filling denominators    *********/
     /***************************************/
     if (fBufferSize == 0) continue;

     fPartBuffer->ResetIter();
         Int_t m = 0;
         while (( partEvent2 = fPartBuffer->Next() ))
          {
            m++;
            if ( (j%fDisplayMixingInfo) == 0)
               Info("ProcessParticles",
                    "Mixing particle %d from current event with particles from event %d",j,-m);
            for(Int_t l = 0; l<partEvent2->GetNumberOfParticles();l++)   //  ... on all particles
              {

                part2= partEvent2->GetParticle(l);
                partpair->SetParticles(part1,part2);

                if( fPairCut->Pass(partpair) ) //check pair cut
                  { //do not meets crietria of the 
                    if( fPairCut->Pass((AliHBTPair*)partpair->GetSwappedPair()) )
                      continue;
                    else 
                     {
                       tmppartpair = (AliHBTPair*)partpair->GetSwappedPair();
                     }
                  }
                else
                 {//meets criteria of the pair cut
                  tmppartpair = partpair;
                 }
                 
                for(ii = 0;ii<fNParticleFunctions;ii++)
                  fParticleFunctions[ii]->ProcessDiffEventParticles(tmppartpair);
                 
             }//for(Int_t l = 0; l<N2;l++)   //  ... on all particles
          }
       }
  delete fPartBuffer->Push(partEvent1);
 //end of loop over events  
  return 0;
}
/*************************************************************************************/ 
Int_t AliHBTAnalysis::ProcessRec(AliAOD* aodrec, AliAOD* /*aodsim*/)
{
  //Does analysis of reconstructed data
  AliVAODParticle * track1, * track2;
  
  AliAOD* trackEvent = aodrec;
  AliAOD* trackEvent1 = new AliAOD();
  trackEvent1->SetOwner(kTRUE);

  AliAOD* trackEvent2;
  
  AliHBTPair tpair;

  AliHBTPair* trackpair = &tpair;
 
  AliHBTPair * tmptrackpair;
  
  register UInt_t ii;
  

  if ( !trackEvent )
   {
     Error("ProcessRecAndSim","Can not get event");
     return 1;
   }


  for (Int_t j = 0; j<trackEvent->GetNumberOfParticles() ; j++)
   {
     /***************************************/
     /******   Looping same events   ********/
     /******   filling numerators    ********/
     /***************************************/
     if ( (j%fDisplayMixingInfo) == 0)
        Info("ProcessTracksAndParticles",
             "Mixing Particle %d with Particles from the same event",j);

     track1= trackEvent->GetParticle(j);

     Bool_t firstcut = fPairCut->GetFirstPartCut()->Pass(track1);

     if (fBufferSize != 0) 
       if ( (firstcut == kFALSE) || ( fPairCut->GetSecondPartCut()->Pass(track1) == kFALSE ) )
        {
          //accepted by any cut
          // we have to copy because reader keeps only one event

          trackEvent1->AddParticle(new AliAODParticle(*track1));
        }

     if (firstcut) continue;

     for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
       fParticleMonitorFunctions[ii]->Process(track1);

     if ( fNParticleFunctions == 0 ) continue;

     for (Int_t k =j+1; k < trackEvent->GetNumberOfParticles() ; k++)
      {
        track2= trackEvent->GetParticle(k);
        if (track1->GetUID() == track2->GetUID()) continue;
        trackpair->SetParticles(track1,track2);

           if(fPairCut->Pass(trackpair)) //check pair cut
            { //do not meets crietria of the 
              if( fPairCut->Pass((AliHBTPair*)trackpair->GetSwappedPair()) ) continue;
              else tmptrackpair = (AliHBTPair*)trackpair->GetSwappedPair();
            }
           else
            {
              tmptrackpair = trackpair;
            }

        for(ii = 0;ii<fNTrackFunctions;ii++)
               fParticleFunctions[ii]->ProcessSameEventParticles(tmptrackpair);

       //end of 2nd loop over Particles from the same event  
      }//for (Int_t k =j+1; k < trackEvent->GetNumberOfParticles() ; k++)

     /***************************************/
     /***** Filling denominators    *********/
     /***************************************/
     if (fBufferSize == 0) continue;

     fTrackBuffer->ResetIter();
         Int_t m = 0;
         while (( trackEvent2 = fTrackBuffer->Next() ))
          {
            m++;
            if ( (j%fDisplayMixingInfo) == 0)
               Info("ProcessParticles",
                    "Mixing Particle %d from current event with Particles from event %d",j,-m);
            for(Int_t l = 0; l<trackEvent2->GetNumberOfParticles();l++)   //  ... on all Particles
              {

                track2= trackEvent2->GetParticle(l);
                trackpair->SetParticles(track1,track2);

                if( fPairCut->Pass(trackpair) ) //check pair cut
                  { //do not meets crietria of the 
                    if( fPairCut->Pass((AliHBTPair*)trackpair->GetSwappedPair()) )
                      continue;
                    else 
                     {
                       tmptrackpair = (AliHBTPair*)trackpair->GetSwappedPair();
                     }
                  }
                else
                 {//meets criteria of the pair cut
                  tmptrackpair = trackpair;
                 }
                 
                for(ii = 0;ii<fNTrackFunctions;ii++)
                  fTrackFunctions[ii]->ProcessDiffEventParticles(tmptrackpair);
                 
             }//for(Int_t l = 0; l<N2;l++)   //  ... on all Particles
          }
       }
  delete fTrackBuffer->Push(trackEvent1);
 //end of loop over events  
  return 0;
}
/*************************************************************************************/ 
     
Int_t AliHBTAnalysis::ProcessRecAndSimNonId(AliAOD* aodrec, AliAOD* aodsim)
{
//Analyzes both reconstructed and simulated data
  if (aodrec == 0x0) 
   {
     Error("ProcessTracksAndParticlesNonIdentAnal","Reconstructed event is NULL");
     return 1;
   }

  if (aodsim == 0x0) 
   {
     Error("ProcessTracksAndParticlesNonIdentAnal","Simulated event is NULL");
     return 1;
   }

  if ( aodrec->GetNumberOfParticles() != aodsim->GetNumberOfParticles() )
   {
     Error("ProcessTracksAndParticlesNonIdentAnal",
           "Number of simulated particles (%d) not equal to number of reconstructed tracks (%d)",
           aodsim->GetNumberOfParticles() , aodrec->GetNumberOfParticles());
     return 2;
   }

 
  AliVAODParticle * part1, * part2;
  AliVAODParticle * track1, * track2;

  static AliAOD aodrec1;
  static AliAOD aodsim1;

  AliAOD * trackEvent1=&aodrec1,*partEvent1=&aodsim1;//Particle that passes first particle cut, this event
  AliAOD * trackEvent2=0x0,*partEvent2=0x0;//Particle that passes second particle cut, this event
  AliAOD * trackEvent3=0x0,*partEvent3=0x0;//Particle that passes second particle cut, events from buffer

  AliAOD* rawtrackEvent = aodrec;//this we get from Reader
  AliAOD* rawpartEvent = aodsim;//this we get from Reader

  static AliHBTPair tpair;
  static AliHBTPair ppair;
  
  AliHBTPair* trackpair = &tpair;
  AliHBTPair* partpair = &ppair;

  register UInt_t ii;
  
  
  trackEvent1 = new AliAOD();
  partEvent1 = new AliAOD();
  trackEvent1->SetOwner(kFALSE);
  partEvent1->SetOwner(kFALSE);;
  
  /********************************/
  /*      Filtering out           */
  /********************************/
  if ( ( (partEvent2==0x0) || (trackEvent2==0x0)) )//in case fBufferSize == 0 and pointers are created do not eneter
   {
     partEvent2  = new AliAOD();
     trackEvent2 = new AliAOD();
     partEvent2->SetOwner(kTRUE);
     trackEvent2->SetOwner(kTRUE);
   }

  FilterOut(partEvent1, partEvent2, rawpartEvent, trackEvent1, trackEvent2, rawtrackEvent);

  for (Int_t j = 0; j<partEvent1->GetNumberOfParticles() ; j++)
   {
     if ( (j%fDisplayMixingInfo) == 0) 
        Info("ProcessTracksAndParticlesNonIdentAnal",
             "Mixing particle %d from current event with particles from current event",j);

     part1= partEvent1->GetParticle(j);
     track1= trackEvent1->GetParticle(j);


     for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
       fParticleMonitorFunctions[ii]->Process(part1);
     for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
       fTrackMonitorFunctions[ii]->Process(track1);
     for(ii = 0; ii<fNParticleAndTrackMonitorFunctions; ii++)
       fParticleAndTrackMonitorFunctions[ii]->Process(track1,part1);

     if (fNoCorrfctns) continue;

     /***************************************/
     /******   filling numerators    ********/
     /****** (particles from event2) ********/
     /***************************************/

     for (Int_t k = 0; k < partEvent2->GetNumberOfParticles() ; k++) //partEvent1 and partEvent2 are particles from the same event but separated to two groups 
      {
        part2= partEvent2->GetParticle(k);
        if (part1->GetUID() == part2->GetUID()) continue;//this is the same particle but with different PID
        partpair->SetParticles(part1,part2);

        track2= trackEvent2->GetParticle(k);
        trackpair->SetParticles(track1,track2);

        if( (this->*fkPassPairProp)(partpair,trackpair) ) //check pair cut
         { //do not meets crietria of the pair cut
          continue; 
         }
        else
         {//meets criteria of the pair cut
          for(ii = 0;ii<fNParticleFunctions;ii++)
                 fParticleFunctions[ii]->ProcessSameEventParticles(partpair);

          for(ii = 0;ii<fNTrackFunctions;ii++)
                 fTrackFunctions[ii]->ProcessSameEventParticles(trackpair);

          for(ii = 0;ii<fNParticleAndTrackFunctions;ii++)
                 fParticleAndTrackFunctions[ii]->ProcessSameEventParticles(trackpair,partpair);
         }
       }

 if ( fBufferSize == 0) continue;//do not mix diff histograms
 /***************************************/
 /***** Filling denominators    *********/
 /***************************************/
 fPartBuffer->ResetIter();
 fTrackBuffer->ResetIter();

 Int_t nmonitor = 0;

 while ( (partEvent3 = fPartBuffer->Next() ) != 0x0)
  {
    trackEvent3 = fTrackBuffer->Next();

    if ( (j%fDisplayMixingInfo) == 0) 
      Info("ProcessTracksAndParticlesNonIdentAnal",
           "Mixing particle %d from current event with particles from event%d",j,-(++nmonitor));

    for (Int_t k = 0; k < partEvent3->GetNumberOfParticles() ; k++)
      {
        part2= partEvent3->GetParticle(k);
        partpair->SetParticles(part1,part2);

        track2= trackEvent3->GetParticle(k);
        trackpair->SetParticles(track1,track2);

        if( (this->*fkPassPairProp)(partpair,trackpair) ) //check pair cut
         { //do not meets crietria of the pair cut
          continue; 
         }
        else
         {//meets criteria of the pair cut
          UInt_t ii;
          for(ii = 0;ii<fNParticleFunctions;ii++)
                 fParticleFunctions[ii]->ProcessDiffEventParticles(partpair);

          for(ii = 0;ii<fNTrackFunctions;ii++)
                 fTrackFunctions[ii]->ProcessDiffEventParticles(trackpair);

          for(ii = 0;ii<fNParticleAndTrackFunctions;ii++)
                 fParticleAndTrackFunctions[ii]->ProcessDiffEventParticles(trackpair,partpair);
         }
       }// for particles event2
     }//while event2
   }//for over particles in event1

 delete fPartBuffer->Push(partEvent2);
 delete fTrackBuffer->Push(trackEvent2);
 
 return 0;
}
/*************************************************************************************/ 
Int_t AliHBTAnalysis::ProcessSimNonId(AliAOD* /*aodrec*/, AliAOD* aodsim)
{
//does analysis of simulated (MC) data in non-identical mode
//i.e. when particles selected by first part. cut are a disjunctive set than particles
//passed by the second part. cut
 if (aodsim == 0x0) 
  {
    return 1;
  }
 
 
  AliVAODParticle * part1, * part2;

  static AliAOD aodsim1;

  AliAOD* partEvent1=&aodsim1;//Particle that passes first particle cut, this event
  AliAOD* partEvent2=0x0;//Particle that passes second particle cut, this event
  AliAOD* partEvent3=0x0;//Particle that passes second particle cut, events from buffer

  AliAOD* rawpartEvent = aodsim;//this we get from Reader

  static AliHBTPair ppair;
  
  AliHBTPair* partpair = &ppair;

  register UInt_t ii;
  
  
  partEvent1 = new AliAOD();
  partEvent1->SetOwner(kFALSE);;
  

  /********************************/
  /*      Filtering out           */
  /********************************/
  if (partEvent2==0x0)//in case fBufferSize == 0 and pointers are created do not eneter
   {
     partEvent2  = new AliAOD();
     partEvent2->SetOwner(kTRUE);
   }

  FilterOut(partEvent1, partEvent2, rawpartEvent);

  for (Int_t j = 0; j<partEvent1->GetNumberOfParticles() ; j++)
   {
     if ( (j%fDisplayMixingInfo) == 0) 
        Info("ProcessParticlesNonIdentAnal",
             "Mixing particle %d from current event with particles from current event",j);

     part1= partEvent1->GetParticle(j);


     for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
       fParticleMonitorFunctions[ii]->Process(part1);

     if (fNParticleFunctions == 0) continue;

     /***************************************/
     /******   filling numerators    ********/
     /****** (particles from event2) ********/
     /***************************************/

     for (Int_t k = 0; k < partEvent2->GetNumberOfParticles() ; k++) //partEvent1 and partEvent2 are particles from the same event but separated to two groups 
      {
        part2= partEvent2->GetParticle(k);
        if (part1->GetUID() == part2->GetUID()) continue;//this is the same particle but with different PID
        partpair->SetParticles(part1,part2);


        if(fPairCut->PassPairProp(partpair) ) //check pair cut
          { //do not meets crietria of the pair cut
              continue; 
          }
        else
         {//meets criteria of the pair cut
          for(ii = 0;ii<fNParticleFunctions;ii++)
              fParticleFunctions[ii]->ProcessSameEventParticles(partpair);
         }
       }

 if ( fBufferSize == 0) continue;//do not mix diff histograms
 /***************************************/
 /***** Filling denominators    *********/
 /***************************************/
 fPartBuffer->ResetIter();

 Int_t nmonitor = 0;

 while ( (partEvent3 = fPartBuffer->Next() ) != 0x0)
  {

    if ( (j%fDisplayMixingInfo) == 0) 
      Info("ProcessParticlesNonIdentAnal",
           "Mixing particle %d from current event with particles from event%d",j,-(++nmonitor));

    for (Int_t k = 0; k < partEvent3->GetNumberOfParticles() ; k++)
      {
        part2= partEvent3->GetParticle(k);
        partpair->SetParticles(part1,part2);


        if(fPairCut->PassPairProp(partpair) ) //check pair cut
          { //do not meets crietria of the pair cut
              continue; 
          }
        else
         {//meets criteria of the pair cut
          for(ii = 0;ii<fNParticleFunctions;ii++)
           {
             fParticleFunctions[ii]->ProcessDiffEventParticles(partpair);
           }
         }
       }// for particles event2
     }//while event2
   }//for over particles in event1

 delete fPartBuffer->Push(partEvent2);
 
 return 0;
}
/*************************************************************************************/ 
Int_t AliHBTAnalysis::ProcessRecNonId(AliAOD* aodrec, AliAOD* /*aodsim*/)
{
//Analyzes both reconstructed and simulated data
 if (aodrec == 0x0) 
  {
    return 1;
  }
 
  AliVAODParticle * track1, * track2;

  static AliAOD aodrec1;

  AliAOD * trackEvent1=&aodrec1;//Particle that passes first particle cut, this event
  AliAOD * trackEvent2=0x0;//Particle that passes second particle cut, this event
  AliAOD * trackEvent3=0x0;//Particle that passes second particle cut, events from buffer

  AliAOD* rawtrackEvent = aodrec;//this we get from Reader

  static AliHBTPair tpair;
  
  AliHBTPair* trackpair = &tpair;

  register UInt_t ii;
  
  
  trackEvent1 = new AliAOD();
  trackEvent1->SetOwner(kFALSE);
  
  /********************************/
  /*      Filtering out           */
  /********************************/
  if ( trackEvent2==0x0 )//in case fBufferSize == 0 and pointers are created do not eneter
   {
     trackEvent2 = new AliAOD();
     trackEvent2->SetOwner(kTRUE);
   }

  FilterOut(trackEvent1, trackEvent2, rawtrackEvent);

  for (Int_t j = 0; j<trackEvent1->GetNumberOfParticles() ; j++)
   {
     if ( (j%fDisplayMixingInfo) == 0) 
        Info("ProcessTracksNonIdentAnal",
             "Mixing particle %d from current event with particles from current event",j);

     track1= trackEvent1->GetParticle(j);


     for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
       fTrackMonitorFunctions[ii]->Process(track1);

     if (fNTrackFunctions == 0x0) continue;

     /***************************************/
     /******   filling numerators    ********/
     /****** (particles from event2) ********/
     /***************************************/

     for (Int_t k = 0; k < trackEvent2->GetNumberOfParticles() ; k++) //partEvent1 and partEvent2 are particles from the same event but separated to two groups 
      {
        track2= trackEvent2->GetParticle(k);
        if (track1->GetUID() == track2->GetUID()) continue;//this is the same particle but with different PID
        trackpair->SetParticles(track1,track2);


        if( fPairCut->PassPairProp(trackpair)) //check pair cut
         { //do not meets crietria of the pair cut
           continue; 
         }
        else
         {//meets criteria of the pair cut
           UInt_t ii;
           for(ii = 0;ii<fNTrackFunctions;ii++)
                  fTrackFunctions[ii]->ProcessSameEventParticles(trackpair);
         }
       }

 if ( fBufferSize == 0) continue;//do not mix diff histograms
 /***************************************/
 /***** Filling denominators    *********/
 /***************************************/
 fTrackBuffer->ResetIter();

 Int_t nmonitor = 0;

 while ( (trackEvent3 = fTrackBuffer->Next() ) != 0x0)
  {
    if ( (j%fDisplayMixingInfo) == 0) 
      Info("ProcessTracksNonIdentAnal",
           "Mixing particle %d from current event with particles from event%d",j,-(++nmonitor));

    for (Int_t k = 0; k < trackEvent3->GetNumberOfParticles() ; k++)
      {
        track2= trackEvent3->GetParticle(k);
        trackpair->SetParticles(track1,track2);

        if( fPairCut->PassPairProp(trackpair)) //check pair cut
         { //do not meets crietria of the pair cut
           continue; 
         }
        else
         {//meets criteria of the pair cut
           for(ii = 0;ii<fNTrackFunctions;ii++)
               fTrackFunctions[ii]->ProcessDiffEventParticles(trackpair);
         }
       }// for particles event2
     }//while event2
   }//for over particles in event1

 delete fTrackBuffer->Push(trackEvent2);
 
 return 0;
}
/*************************************************************************************/ 

void AliHBTAnalysis::Process(Option_t* option)
{
 //default option  = "TracksAndParticles"
 //Main method of the HBT Analysis Package
 //It triggers reading with the global cut (default is an empty cut)
 //Than it checks options and data which are read
 //if everything is OK, then it calls one of the looping methods
 //depending on tfReaderhe option
 //These methods differs on what thay are looping on
 //
 //        METHOD                           OPTION
 //--------------------------------------------------------------------
 //ProcessTracksAndParticles    -     "TracksAndParticles"
 //     DEFAULT
 //     it loops over both, tracks(reconstructed) and particles(simulated)
 //     all function gethered in all 3 lists are called for each (double)pair
 //                             
 //ProcessTracks                -          "Tracks" 
 //     it loops only on tracks(reconstructed),
 //     functions ONLY from fTrackFunctions list are called
 //
 //ProcessParticles             -         "Particles"
 //     it loops only on particles(simulated),
 //     functions ONLY from fParticleAndTrackFunctions list are called
 //
 //
 if (!fReader) 
  {
   Error("Process","The reader is not set");
   return;
  }
 
 const char *oT = strstr(option,"Tracks");
 const char *oP = strstr(option,"Particles");
 
 Bool_t nonid = IsNonIdentAnalysis();
 
 Init();
 
 if(oT && oP)
  { 
    if (nonid) ProcessTracksAndParticlesNonIdentAnal();
    else ProcessTracksAndParticles();
    return;
  }

 if(oT)
  {
    if (nonid) ProcessTracksNonIdentAnal();
    else ProcessTracks();
    return;
  }
 
 if(oP)
  {
    if (nonid) ProcessParticlesNonIdentAnal();
    else ProcessParticles();
    return;
  }
 
}
/*************************************************************************************/ 

void AliHBTAnalysis::ProcessTracksAndParticles()
{
//Makes analysis for both tracks and particles
//mainly for resolution study and analysies with weighting algirithms
//In order to minimize calling AliRun::GetEvent (we need at one time particles from different events),
//the loops are splited
  
// cut on particles only -- why?
// - PID: when we make resolution analysis we want to take only tracks with correct PID
// We need cut on tracks because there are data characteristic to 
  
  AliAOD * trackEvent, *partEvent;

  fReader->Rewind();
  while (fReader->Next() == kFALSE)
    {
      partEvent = fReader->GetEventSim();
      trackEvent = fReader->GetEventRec();
      ProcessRecAndSim(trackEvent,partEvent);

    }//while (fReader->Next() == kFALSE)
    
} 
/*************************************************************************************/

void AliHBTAnalysis::ProcessTracks()
{
//In order to minimize calling AliRun::GetEvent (we need at one time particles from different events),
//the loops are splited
  AliAOD * trackEvent;
  fReader->Rewind();
  while (fReader->Next() == kFALSE)
    {
      trackEvent = fReader->GetEventRec();
      ProcessRec(trackEvent,0x0);
    }//while (fReader->Next() == kFALSE)
}

/*************************************************************************************/

void AliHBTAnalysis::ProcessParticles()
{
//In order to minimize calling AliRun::GetEvent (we need at one time particles from different events),
//the loops are splited
  AliAOD * partEvent;
  fReader->Rewind();
  while (fReader->Next() == kFALSE)
    {
      partEvent = fReader->GetEventSim();
      ProcessSim(0x0,partEvent);
    }//while (fReader->Next() == kFALSE)
}
/*************************************************************************************/

void AliHBTAnalysis::WriteFunctions()
{
//Calls Write for all defined functions in analysis
//== writes all results
  TFile* oututfile = 0x0;
  if (fOutputFileName)
   {
     oututfile = TFile::Open(*fOutputFileName,"update");
   }
  UInt_t ii;
  for(ii = 0;ii<fNParticleFunctions;ii++)
   {
     if (AliVAODParticle::GetDebug()>5)
      {
        Info("WriteFunctions","Writing ParticleFunction %#x",fParticleFunctions[ii]);
        Info("WriteFunctions","Writing ParticleFunction %s",fParticleFunctions[ii]->Name());
      }
     fParticleFunctions[ii]->Write();
   }
                 
  for(ii = 0;ii<fNTrackFunctions;ii++)
   {
     if (AliVAODParticle::GetDebug()>5)
      {
        Info("WriteFunctions","Writing TrackFunction %#x",fTrackFunctions[ii]);
        Info("WriteFunctions","Writing TrackFunction %s",fTrackFunctions[ii]->Name());
      }
     fTrackFunctions[ii]->Write();
   }
                 
  for(ii = 0;ii<fNParticleAndTrackFunctions;ii++)
   {
     if (AliVAODParticle::GetDebug()>5)
      {
        Info("WriteFunctions","Writing ParticleAndTrackFunction %#x",fParticleAndTrackFunctions[ii]);
        Info("WriteFunctions","Writing ParticleAndTrackFunction %s",fParticleAndTrackFunctions[ii]->Name());
      } 
     fParticleAndTrackFunctions[ii]->Write();
   }

  for(ii = 0;ii<fNParticleMonitorFunctions;ii++)
   {
     if (AliVAODParticle::GetDebug()>5)
      {
        Info("WriteFunctions","Writing ParticleMonitorFunction %#x",fParticleMonitorFunctions[ii]);
        Info("WriteFunctions","Writing ParticleMonitorFunction %s",fParticleMonitorFunctions[ii]->Name());
      }
     fParticleMonitorFunctions[ii]->Write();
   }

  for(ii = 0;ii<fNTrackMonitorFunctions;ii++)
   {
     if (AliVAODParticle::GetDebug()>5)
      {
        Info("WriteFunctions","Writing TrackMonitorFunction %#x",fTrackMonitorFunctions[ii]);
        Info("WriteFunctions","Writing TrackMonitorFunction %s",fTrackMonitorFunctions[ii]->Name());
      }
     fTrackMonitorFunctions[ii]->Write();
   }

  for(ii = 0;ii<fNParticleAndTrackMonitorFunctions;ii++)
   {
     if (AliVAODParticle::GetDebug()>5)
      {
       Info("WriteFunctions","Writing ParticleAndTrackMonitorFunction %#x",fParticleAndTrackMonitorFunctions[ii]);
       Info("WriteFunctions","Writing ParticleAndTrackMonitorFunction %s",fParticleAndTrackMonitorFunctions[ii]->Name());
      }
     fParticleAndTrackMonitorFunctions[ii]->Write();
   }
  delete oututfile; 
}
/*************************************************************************************/

void AliHBTAnalysis::SetOutputFileName(const char* fname)
{
  //Sets fiele name where to dump results, 
  //if not specified reults are written to gDirectory
  if (fname == 0x0)
   {
     delete fOutputFileName;
     fOutputFileName = 0x0;
     return;
   }
  if ( strcmp(fname,"") == 0 ) 
   {
     delete fOutputFileName;
     fOutputFileName = 0x0;
     return;
   }
  if (fOutputFileName == 0x0) fOutputFileName = new TString(fname);
  else *fOutputFileName = fname;
}
/*************************************************************************************/

void AliHBTAnalysis::SetGlobalPairCut(AliAODPairCut* cut)
{
//Sets the global cut
  if (cut == 0x0)
   {
     Error("AliHBTAnalysis::SetGlobalPairCut","Pointer is NULL. Ignoring");
   }
  delete fPairCut;
  fPairCut = (AliAODPairCut*)cut->Clone();
}

/*************************************************************************************/

void AliHBTAnalysis::AddTrackFunction(AliHBTOnePairFctn* f)
{
//Adds track function
  if (f == 0x0) return;
  if (fNTrackFunctions == fgkFctnArraySize)
   {
    Error("AliHBTAnalysis::AddTrackFunction","Can not add this function, not enough place in the array.");
   }
  fTrackFunctions[fNTrackFunctions] = f;
  fNTrackFunctions++;
}
/*************************************************************************************/ 

void AliHBTAnalysis::AddParticleFunction(AliHBTOnePairFctn* f)
{
//adds particle function
  if (f == 0x0) return;
  
  if (fNParticleFunctions == fgkFctnArraySize)
   {
    Error("AliHBTAnalysis::AddParticleFunction","Can not add this function, not enough place in the array.");
   }
  fParticleFunctions[fNParticleFunctions] = f;
  fNParticleFunctions++;
}
/*************************************************************************************/ 

void AliHBTAnalysis::AddParticleAndTrackFunction(AliHBTTwoPairFctn* f)
{
//add resolution function
  if (f == 0x0) return;
  if (fNParticleAndTrackFunctions == fgkFctnArraySize)
   {
    Error("AliHBTAnalysis::AddParticleAndTrackFunction","Can not add this function, not enough place in the array.");
   }  
  fParticleAndTrackFunctions[fNParticleAndTrackFunctions] = f;
  fNParticleAndTrackFunctions++;
}
/*************************************************************************************/ 

void AliHBTAnalysis::AddParticleMonitorFunction(AliHBTMonOneParticleFctn* f)
{
//add particle monitoring function
  if (f == 0x0) return;

  if (fNParticleMonitorFunctions == fgkFctnArraySize)
    {
      Error("AliHBTAnalysis::AddParticleMonitorFunction","Can not add this function, not enough place in the array.");
   }
  fParticleMonitorFunctions[fNParticleMonitorFunctions] = f;
  fNParticleMonitorFunctions++;
}
/*************************************************************************************/ 

void AliHBTAnalysis::AddTrackMonitorFunction(AliHBTMonOneParticleFctn* f)
{
//add track monitoring function
  if (f == 0x0) return;

  if (fNTrackMonitorFunctions == fgkFctnArraySize)
   {
    Error("AliHBTAnalysis::AddTrackMonitorFunction","Can not add this function, not enough place in the array.");
   }
  fTrackMonitorFunctions[fNTrackMonitorFunctions] = f;
  fNTrackMonitorFunctions++;
}
/*************************************************************************************/ 

void AliHBTAnalysis::AddParticleAndTrackMonitorFunction(AliHBTMonTwoParticleFctn* f)
{
//add resolution monitoring function
  if (f == 0x0) return;
  if (fNParticleAndTrackMonitorFunctions == fgkFctnArraySize)
    {
      Error("AliHBTAnalysis::AddParticleAndTrackMonitorFunction","Can not add this function, not enough place in the array.");
    }
  fParticleAndTrackMonitorFunctions[fNParticleAndTrackMonitorFunctions] = f;
  fNParticleAndTrackMonitorFunctions++;
}


/*************************************************************************************/ 
/*************************************************************************************/  

Bool_t AliHBTAnalysis::RunCoherencyCheck()
{
 //Checks if both HBTRuns are similar
 //return true if error found
 //if they seem to be OK return false
/* 

 Int_t i;  
 Info("RunCoherencyCheck","Checking HBT Runs Coherency");

//When we use non-buffering reader this is a big waste of time -> We need to read all data to check it
//and reader is implemented safe in this case anyway
// Info("RunCoherencyCheck","Number of events ...");
// if (fReader->GetNumberOfPartEvents() == fReader->GetNumberOfTrackEvents() ) //check whether there is the same  number of events
//  {
//    Info("RunCoherencyCheck","OK. %d found\n",fReader->GetNumberOfTrackEvents());
//  }
// else
//  { //if not the same - ERROR
//   Error("RunCoherencyCheck",
//         "Number of simulated events (%d) is not equal to number of reconstructed events(%d)",
//         fReader->GetNumberOfPartEvents(),fReader->GetNumberOfTrackEvents());
//   return kTRUE;
//  }
 
 Info("RunCoherencyCheck","Checking number of Particles AND Particles Types in each event ...");
 
 AliAOD *partEvent;
 AliAOD *trackEvent;
 for( i = 0; i<fReader->GetNumberOfTrackEvents();i++)
  {
    partEvent= fReader->GetEventSim(i); //gets the "ith" event 
    trackEvent = fReader->GetEventRec(i);
    
    if ( (partEvent == 0x0) && (partEvent == 0x0) ) continue;
    if ( (partEvent == 0x0) || (partEvent == 0x0) )
     {
       Error("RunCoherencyCheck",
             "One event is NULL and the other one not. Event Number %d",i);
       return kTRUE;    
     }
    
    if ( partEvent->GetNumberOfParticles() != trackEvent->GetNumberOfParticles() )
     {
       Error("RunCoherencyCheck",
             "Event %d: Number of simulated particles (%d) not equal to number of reconstructed tracks (%d)",
             i,partEvent->GetNumberOfParticles() , trackEvent->GetNumberOfParticles());
       return kTRUE;
     }
    else
     for (Int_t j = 0; j<partEvent->GetNumberOfParticles(); j++)
      {
        if( partEvent->GetParticle(j)->GetPdgCode() != trackEvent->GetParticle(j)->GetPdgCode() )
         {
           Error("RunCoherencyCheck",
                 "Event %d: Particle %d: PID of simulated particle (%d) not the same of reconstructed track (%d)",
                 i,j, partEvent->GetParticle(j)->GetPdgCode(),trackEvent->GetParticle(j)->GetPdgCode() );
           return kTRUE;
           
         }
      }
  }
 Info("RunCoherencyCheck","  Done");
 Info("RunCoherencyCheck","  Everything looks OK");
*/ 
 return kFALSE;
}

/*************************************************************************************/  
 
void AliHBTAnalysis::ProcessTracksAndParticlesNonIdentAnal()
{
//Performs analysis for both, tracks and particles
  AliAOD* rawtrackEvent, * rawpartEvent;
  fReader->Rewind();
 
  Info("ProcessTracksAndParticlesNonIdentAnal","**************************************");
  Info("ProcessTracksAndParticlesNonIdentAnal","*****  NON IDENT MODE ****************");
  Info("ProcessTracksAndParticlesNonIdentAnal","**************************************");

  for (Int_t i = 0;;i++)//infinite loop
    {
      if (fReader->Next()) break; //end when no more events available
      
      rawpartEvent  = fReader->GetEventSim();
      rawtrackEvent = fReader->GetEventRec();
      
      ProcessRecAndSimNonId(rawtrackEvent,rawpartEvent);
  }//end of loop over events (1)
}
/*************************************************************************************/  
 
void AliHBTAnalysis::ProcessTracksNonIdentAnal()
{
//Process Tracks only with non identical mode
  AliAOD * rawtrackEvent;
  fReader->Rewind();
  
  Info("ProcessTracksNonIdentAnal","**************************************");
  Info("ProcessTracksNonIdentAnal","*****  NON IDENT MODE ****************");
  Info("ProcessTracksNonIdentAnal","**************************************");
  
  for (Int_t i = 0;;i++)//infinite loop
    {
      if (fReader->Next()) break; //end when no more events available
      rawtrackEvent = fReader->GetEventRec();
      ProcessRecNonId(rawtrackEvent,0x0);
  }//end of loop over events (1)
}
/*************************************************************************************/  

void AliHBTAnalysis::ProcessParticlesNonIdentAnal()
{
//process paricles only with non identical mode
  AliAOD * rawpartEvent = 0x0;
  fReader->Rewind();
  
  Info("ProcessParticlesNonIdentAnal","**************************************");
  Info("ProcessParticlesNonIdentAnal","*****  NON IDENT MODE ****************");
  Info("ProcessParticlesNonIdentAnal","**************************************");

  for (Int_t i = 0;;i++)//infinite loop
    {
      if (fReader->Next()) break; //end when no more events available
      
      rawpartEvent  = fReader->GetEventSim();
      ProcessSimNonId(0x0,rawpartEvent);
  }//end of loop over events (1)
}

/*************************************************************************************/  
void AliHBTAnalysis::FilterOut(AliAOD* outpart1, AliAOD* outpart2, AliAOD* inpart,
                     AliAOD* outtrack1, AliAOD* outtrack2, AliAOD* intrack) const
{
 //Puts particles accepted as a first particle by global cut in out1
 //and as a second particle in out2

  AliVAODParticle* part, *track;

  outpart1->Reset();
  outpart2->Reset();
  outtrack1->Reset();
  outtrack2->Reset();
  
  Bool_t in1, in2;
  
  for (Int_t i = 0; i < inpart->GetNumberOfParticles(); i++)
   {
     in1 = in2 = kTRUE;
     part = inpart->GetParticle(i);
     track = intrack->GetParticle(i);
     
     if ( ((this->*fkPass1)(part,track))  ) in1 = kFALSE; //if part  is rejected by cut1, in1 is false
     if ( ((this->*fkPass2)(part,track))  ) in2 = kFALSE; //if part  is rejected by cut2, in2 is false
     
     if (gDebug)//to be removed in real analysis     
     if ( in1 && in2 ) //both cuts accepted, should never happen, just in case
      {
      //Particle accpted by both cuts
       Error("FilterOut","Particle accepted by both cuts");
       continue;
      }

     if (in1)
      {
        outpart1->AddParticle(part);
        outtrack1->AddParticle(track);
        continue;
      }
     
     if (in2)
      {
        outpart2->AddParticle(new AliAODParticle(*part));
        outtrack2->AddParticle(new AliAODParticle(*track));
        continue;
      }
   }
}
/*************************************************************************************/  

void AliHBTAnalysis::FilterOut(AliAOD* out1, AliAOD* out2, AliAOD* in) const
{
 //Puts particles accepted as a first particle by global cut in out1
 //and as a second particle in out2
  AliVAODParticle* part;
  
  out1->Reset();
  out2->Reset();
  
  AliAODParticleCut *cut1 = fPairCut->GetFirstPartCut();
  AliAODParticleCut *cut2 = fPairCut->GetSecondPartCut();
  
  Bool_t in1, in2;
  
  for (Int_t i = 0; i < in->GetNumberOfParticles(); i++)
   {
     in1 = in2 = kTRUE;
     part = in->GetParticle(i);
     
     if ( cut1->Pass(part) ) in1 = kFALSE; //if part is rejected by cut1, in1 is false
     if ( cut2->Pass(part) ) in2 = kFALSE; //if part is rejected by cut2, in2 is false
     
     if (gDebug)//to be removed in real analysis     
     if ( in1 && in2 ) //both cuts accepted, should never happen, just in case
      {
      //Particle accpted by both cuts
       Error("FilterOut","Particle accepted by both cuts");
       continue;
      }

     if (in1)
      { 
        out1->AddParticle(part);
        continue;
      }
     
     if (in2)
      {
        out2->AddParticle(part);
        continue;
      }
   }
}
/*************************************************************************************/ 

Bool_t AliHBTAnalysis::IsNonIdentAnalysis()
{
 //checks if it is possible to use special analysis for non identical particles
 //it means - in global pair cut first particle id is different than second one 
 //and both are different from 0
 //in the future is possible to perform more sophisticated check 
 //if cuts have excluding requirements
 
 if (fPairCut->IsEmpty()) 
   return kFALSE;
 
 if (fPairCut->GetFirstPartCut()->IsEmpty()) 
   return kFALSE;

 if (fPairCut->GetSecondPartCut()->IsEmpty()) 
   return kFALSE;
 
 Int_t id1 = fPairCut->GetFirstPartCut()->GetPID();
 Int_t id2 = fPairCut->GetSecondPartCut()->GetPID();

 if ( (id1==0) || (id2==0) || (id1==id2) ) 
   return kFALSE;
 
 return kTRUE;
}
/*************************************************************************************/ 

void AliHBTAnalysis::SetApparentVertex(Double_t x, Double_t y, Double_t z)
{ 
  //Sets apparent vertex
  // All events have to be moved to the same vertex position in order to
  // to be able to comare any space positions (f.g. anti-merging)
  // This method defines this position
  
  fVertexX = x;
  fVertexY = y;
  fVertexZ = z;
}
/*************************************************************************************/ 

void AliHBTAnalysis::PressAnyKey()
{ 
 //small utility function that helps to make comfortable macros
  char c;
  int nread = -1;
  fcntl(0,  F_SETFL, O_NONBLOCK);
  ::Info("","Press Any Key to continue ...");
  while (nread<1) 
   {
     nread = read(0, &c, 1);
     gSystem->ProcessEvents();
   }
}

/*************************************************************************************/ 

