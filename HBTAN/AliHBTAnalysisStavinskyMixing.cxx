#include "AliHBTAnalysisStavinskyMixing.h"
//_________________________________________________________
////////////////////////////////////////////////////////////////////////////
//
// class AliHBTAnalysisStavinskyMixing
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
 

ClassImp(AliHBTAnalysisStavinskyMixing)


Int_t AliHBTAnalysisStavinskyMixing::ProcessRecAndSim(AliAOD* aodrec, AliAOD* aodsim)
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

          partEvent1->AddParticle(part1);
          trackEvent1->AddParticle(track1);
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


        /***************************************/
        /***** Filling denominators    *********/
        /***************************************/

        tmppartpair->MirrorSecond();
        tmptrackpair->MirrorSecond();
        
        if((this->*fkPass)(tmppartpair,tmptrackpair) == kFALSE)
         {
            for(ii = 0;ii<fNParticleFunctions;ii++)
              fParticleFunctions[ii]->ProcessDiffEventParticles(tmppartpair);

            for(ii = 0;ii<fNTrackFunctions;ii++)
              fTrackFunctions[ii]->ProcessDiffEventParticles(tmptrackpair);

            for(ii = 0;ii<fNParticleAndTrackFunctions;ii++)
              fParticleAndTrackFunctions[ii]->ProcessDiffEventParticles(tmptrackpair,tmppartpair);

         } 
        
        tmppartpair->DeleteSecond();
        tmptrackpair->DeleteSecond();
        
       //end of 2nd loop over particles from the same event  
      }//for (Int_t k =j+1; k < partEvent->GetNumberOfParticles() ; k++)

    //end of loop over particles from first event
   }//for (Int_t j = 0; j<partEvent->GetNumberOfParticles() ; j++)
  delete fPartBuffer->Push(partEvent1);
  delete fTrackBuffer->Push(trackEvent1);
 //end of loop over events  
  return 0;
}
/*************************************************************************************/ 

Int_t AliHBTAnalysisStavinskyMixing::ProcessSim(AliAOD* /*aodrec*/, AliAOD* aodsim)
{
  //Does analysis of simulated data
  AliVAODParticle * part1, * part2;
  
  AliAOD* partEvent = aodsim;
  AliAOD* partEvent1 = new AliAOD();

  
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
     Warning("ProcessTracksAndParticles","Non-std MIXING");
     if ( (j%fDisplayMixingInfo) == 0)
        Info("ProcessTracksAndParticles",
             "Mixing particle %d with particles from the same event",j);

     part1= partEvent->GetParticle(j);

     Bool_t firstcut = fPairCut->GetFirstPartCut()->Rejected(part1);

     if (fBufferSize != 0) 
       if ( (firstcut == kFALSE) || ( fPairCut->GetSecondPartCut()->Rejected(part1) == kFALSE ) )
        {
          //accepted by any cut
          // we have to copy because reader keeps only one event

          partEvent1->AddParticle(part1);
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

           if(fPairCut->Rejected(partpair)) //check pair cut
            { //do not meets crietria of the 
              if( fPairCut->Rejected((AliHBTPair*)partpair->GetSwappedPair()) ) continue;
              else tmppartpair = (AliHBTPair*)partpair->GetSwappedPair();
            }
           else
            {
              tmppartpair = partpair;
            }

        for(ii = 0;ii<fNParticleFunctions;ii++)
               fParticleFunctions[ii]->ProcessSameEventParticles(tmppartpair);


        /***************************************/
        /***** Filling denominators    *********/
        /***************************************/
        
        tmppartpair->MirrorSecond();
        
        if(fPairCut->Rejected(partpair) == kFALSE)
         {
           for(ii = 0;ii<fNParticleFunctions;ii++)
             fParticleFunctions[ii]->ProcessDiffEventParticles(tmppartpair);
         }  
        tmppartpair->DeleteSecond();
        
       //end of 2nd loop over particles from the same event  
      }//for (Int_t k =j+1; k < partEvent->GetNumberOfParticles() ; k++)

       }
  delete fPartBuffer->Push(partEvent1);
 //end of loop over events  
  return 0;
}
/*************************************************************************************/ 
Int_t AliHBTAnalysisStavinskyMixing::ProcessRec(AliAOD* aodrec, AliAOD* /*aodsim*/)
{
  //Does analysis of reconstructed data
  AliVAODParticle * track1, * track2;
  
  AliAOD* trackEvent = aodrec;
  AliAOD* trackEvent1 = new AliAOD();

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

     Bool_t firstcut = fPairCut->GetFirstPartCut()->Rejected(track1);

     if (fBufferSize != 0) 
       if ( (firstcut == kFALSE) || ( fPairCut->GetSecondPartCut()->Rejected(track1) == kFALSE ) )
        {
          //accepted by any cut
          // we have to copy because reader keeps only one event

          trackEvent1->AddParticle(track1);
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

           if(fPairCut->Rejected(trackpair)) //check pair cut
            { //do not meets crietria of the 
              if( fPairCut->Rejected((AliHBTPair*)trackpair->GetSwappedPair()) ) continue;
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

                if( fPairCut->Rejected(trackpair) ) //check pair cut
                  { //do not meets crietria of the 
                    if( fPairCut->Rejected((AliHBTPair*)trackpair->GetSwappedPair()) )
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
     
Int_t AliHBTAnalysisStavinskyMixing::ProcessRecAndSimNonId(AliAOD* aodrec, AliAOD* aodsim)
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
  trackEvent1->Reset();
  partEvent1->Reset();
  AliAOD * trackEvent2=0x0,*partEvent2=0x0;//Particle that passes second particle cut, this event
  AliAOD * trackEvent3=0x0,*partEvent3=0x0;//Particle that passes second particle cut, events from buffer

  AliAOD* rawtrackEvent = aodrec;//this we get from Reader
  AliAOD* rawpartEvent = aodsim;//this we get from Reader

  static AliHBTPair tpair;
  static AliHBTPair ppair;
  
  AliHBTPair* trackpair = &tpair;
  AliHBTPair* partpair = &ppair;

  register UInt_t ii;
  
  /********************************/
  /*      Filtering out           */
  /********************************/
  if ( ( (partEvent2==0x0) || (trackEvent2==0x0)) )//in case fBufferSize == 0 and pointers are created do not eneter
   {
     partEvent2  = new AliAOD();
     trackEvent2 = new AliAOD();
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
Int_t AliHBTAnalysisStavinskyMixing::ProcessSimNonId(AliAOD* /*aodrec*/, AliAOD* aodsim)
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
  partEvent1->Reset();
  AliAOD* partEvent2=0x0;//Particle that passes second particle cut, this event
  AliAOD* partEvent3=0x0;//Particle that passes second particle cut, events from buffer

  AliAOD* rawpartEvent = aodsim;//this we get from Reader

  static AliHBTPair ppair;
  
  AliHBTPair* partpair = &ppair;

  register UInt_t ii;
  
  /********************************/
  /*      Filtering out           */
  /********************************/
  if (partEvent2==0x0)//in case fBufferSize == 0 and pointers are created do not eneter
   {
     partEvent2  = new AliAOD();
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
Int_t AliHBTAnalysisStavinskyMixing::ProcessRecNonId(AliAOD* aodrec, AliAOD* /*aodsim*/)
{
//Analyzes both reconstructed and simulated data
 if (aodrec == 0x0) 
  {
    return 1;
  }
 
  AliVAODParticle * track1, * track2;

  static AliAOD aodrec1;
  AliAOD * trackEvent1=&aodrec1;//Particle that passes first particle cut, this event
  trackEvent1->Reset();
  AliAOD * trackEvent2=0x0;//Particle that passes second particle cut, this event
  AliAOD * trackEvent3=0x0;//Particle that passes second particle cut, events from buffer
  AliAOD* rawtrackEvent = aodrec;//this we get from Reader

  static AliHBTPair tpair;
  
  AliHBTPair* trackpair = &tpair;

  register UInt_t ii;
  
  
  /********************************/
  /*      Filtering out           */
  /********************************/
  if ( trackEvent2==0x0 )//in case fBufferSize == 0 and pointers are created do not eneter
   {
     trackEvent2 = new AliAOD();
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
