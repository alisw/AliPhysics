#include "AliHBTAnalysis.h"
#include "AliHBTRun.h"
#include "AliHBTEvent.h"
#include "AliHBTReader.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"
#include "AliHBTPair.h"
#include "AliHBTPairCut.h"
#include "AliHBTFunction.h"
#include "AliHBTMonitorFunction.h"
#include "AliHBTEventBuffer.h"
 
#include <TBenchmark.h>
#include <TList.h>

//_________________________________________________________
///////////////////////////////////////////////////////////
//
//Central Object Of HBTAnalyser: 
//This class performs main looping within HBT Analysis
//User must plug a reader of Type AliHBTReader
//User plugs in coorelation and monitor functions
//as well as monitor functions
//
//HBT Analysis Tool, which is integral part of AliRoot,
//ALICE Off-Line framework:
//
//Piotr.Skowronski@cern.ch
//more info: http://alisoft.cern.ch/people/skowron/analyzer/index.html
//
//_________________________________________________________

ClassImp(AliHBTAnalysis)

const UInt_t AliHBTAnalysis::fgkFctnArraySize = 100;
const UInt_t AliHBTAnalysis::fgkDefaultMixingInfo = 1000;
const Int_t  AliHBTAnalysis::fgkDefaultBufferSize = 5;

AliHBTAnalysis::AliHBTAnalysis():
  fReader(0x0),
  fNTrackFunctions(0),
  fNParticleFunctions(0),
  fNParticleAndTrackFunctions(0),
  fNTrackMonitorFunctions(0),
  fNParticleMonitorFunctions(0), 
  fNParticleAndTrackMonitorFunctions(0),
  fBufferSize(2),
  fDisplayMixingInfo(fgkDefaultMixingInfo),
  fIsOwner(kFALSE)
 {
//default constructor
   fTrackFunctions = new AliHBTOnePairFctn* [fgkFctnArraySize];
   fParticleFunctions = new AliHBTOnePairFctn* [fgkFctnArraySize];
   fParticleAndTrackFunctions = new AliHBTTwoPairFctn* [fgkFctnArraySize];
   
   fParticleMonitorFunctions = new AliHBTMonOneParticleFctn* [fgkFctnArraySize];    
   fTrackMonitorFunctions = new AliHBTMonOneParticleFctn* [fgkFctnArraySize];    
   fParticleAndTrackMonitorFunctions = new AliHBTMonTwoParticleFctn* [fgkFctnArraySize];    

   fPairCut = new AliHBTEmptyPairCut();//empty cut - accepts all particles
 }
/*************************************************************************************/ 

AliHBTAnalysis::AliHBTAnalysis(const AliHBTAnalysis& in):
  TObject(in),
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
  fPairCut(0x0),
  fBufferSize(fgkDefaultBufferSize),
  fDisplayMixingInfo(fgkDefaultMixingInfo),
  fIsOwner(kFALSE)
 {
//copy constructor
   Fatal("AliHBTAnalysis(const AliHBTAnalysis&)","Sensless");
 }
/*************************************************************************************/ 
const AliHBTAnalysis& AliHBTAnalysis::operator=(const AliHBTAnalysis& /*right*/)
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
    DeleteFunctions();
   delete [] fTrackFunctions;
   delete [] fParticleFunctions;
   delete [] fParticleAndTrackFunctions;
   
   delete [] fParticleMonitorFunctions; 
   delete [] fTrackMonitorFunctions; 
   delete [] fParticleAndTrackMonitorFunctions; 

   delete fPairCut; // always have an copy of an object - we create we dstroy
 }

/*************************************************************************************/ 

void AliHBTAnalysis::DeleteFunctions()
{
 //Deletes all functions added to analysis
 UInt_t ii;
 for(ii = 0;ii<fNParticleFunctions;ii++)
   delete fParticleFunctions[ii];
 fNParticleFunctions = 0;
                
 for(ii = 0;ii<fNTrackFunctions;ii++)
   delete fTrackFunctions[ii];
 fNTrackFunctions = 0;
 
 for(ii = 0;ii<fNParticleAndTrackFunctions;ii++)
   delete fParticleAndTrackFunctions[ii];
 fNParticleAndTrackFunctions = 0;
 
 for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
   delete fParticleMonitorFunctions[ii];
 fNParticleMonitorFunctions = 0;
   
 for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
   delete fTrackMonitorFunctions[ii];
 fNTrackMonitorFunctions = 0;
   
 for(ii = 0; ii<fNParticleAndTrackMonitorFunctions; ii++)
   delete fParticleAndTrackMonitorFunctions[ii];
 fNParticleAndTrackMonitorFunctions = 0;
}
/*************************************************************************************/ 

void AliHBTAnalysis::Init()
{
//Initializeation method
//calls Init for all functions
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
//In order to minimize calling AliRun::GetEvent (we need at one time particles from different events),
//the loops are splited
  
  
  AliHBTParticle * part1, * part2;
  AliHBTParticle * track1, * track2;
  
  AliHBTEvent * trackEvent, *partEvent;
  AliHBTEvent * trackEvent1 = 0x0,*partEvent1 = 0x0;
  AliHBTEvent * trackEvent2,*partEvent2;
  
//  Int_t N1, N2, N=0; //number of particles in current event(we prcess two events in one time)
  
//  Int_t nev = fReader->GetNumberOfTrackEvents();
  AliHBTPair * trackpair = new AliHBTPair();
  AliHBTPair * partpair = new AliHBTPair();

  AliHBTPair * tmptrackpair;//temprary pointers to pairs
  AliHBTPair * tmppartpair;

  AliHBTEventBuffer partbuffer(fBufferSize);
  AliHBTEventBuffer trackbuffer(fBufferSize);
  
  
  register UInt_t ii;
  
  Bool_t nocorrfctns = (fNParticleFunctions == 0) && (fNTrackFunctions == 0) && (fNParticleAndTrackFunctions == 0);
  
  fReader->Rewind();
  Int_t i = -1;
  while (fReader->Next() == kFALSE)
    {
      i++;
      partEvent= fReader->GetParticleEvent();
      trackEvent = fReader->GetTrackEvent();
      
      if ( !partEvent || !trackEvent ) 
       {
         Error("ProcessTracksAndParticles","Can not get event");
         return;
       }
       
      if ( partEvent->GetNumberOfParticles() != trackEvent->GetNumberOfParticles() )
       {
         Fatal("ProcessTracksAndParticles",
               "Event %d: Number of simulated particles (%d) not equal to number of reconstructed tracks (%d)",
               i,partEvent->GetNumberOfParticles() , trackEvent->GetNumberOfParticles());
         return;
       }
      
      if(partEvent1 == 0x0) 
       {
         partEvent1 = new AliHBTEvent();
         partEvent1->SetOwner(kTRUE);

         trackEvent1 = new AliHBTEvent();
         trackEvent1->SetOwner(kTRUE);
       }
      else
       {
         partEvent1->Reset();
         trackEvent1->Reset();
       }
       
      for (Int_t j = 0; j<partEvent->GetNumberOfParticles() ; j++)
       {
         /***************************************/
         /******   Looping same events   ********/
         /******   filling numerators    ********/
         /***************************************/
         if ( (j%fDisplayMixingInfo) == 0)
            Info("ProcessTracksAndParticles",
                 "Mixing particle %d from event %d with particles from event %d",j,i,i);

         part1= partEvent->GetParticle(j);
         track1= trackEvent->GetParticle(j);
  
//PID imperfections ???
//         if( part1->GetPdgCode() != track1->GetPdgCode() )
//          {
//            Fatal("ProcessTracksAndParticles",
//                  "Event %d: Particle %d: PID of simulated particle (%d) not the same of reconstructed track (%d)",
//                  i,j, part1->GetPdgCode(),track1->GetPdgCode() );
//            return;
//          }
         
         Bool_t firstcut = fPairCut->GetFirstPartCut()->Pass(part1);
         
         if ( (firstcut == kFALSE) || (fPairCut->GetSecondPartCut()->Pass(part1) == kFALSE) )
          {
            //accepted by any cut
            // we have to copy because reader keeps only one event
            
            partEvent1->AddParticle(new AliHBTParticle(*part1));
            trackEvent1->AddParticle(new AliHBTParticle(*track1));
          }

         if (firstcut) continue;
         
         for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
           fParticleMonitorFunctions[ii]->Process(part1);
         for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
           fTrackMonitorFunctions[ii]->Process(track1);
         for(ii = 0; ii<fNParticleAndTrackMonitorFunctions; ii++)
           fParticleAndTrackMonitorFunctions[ii]->Process(track1,part1);

         if (nocorrfctns) continue;
        
         for (Int_t k =j+1; k < partEvent->GetNumberOfParticles() ; k++)
          {
            part2= partEvent->GetParticle(k);
            if (part1->GetUID() == part2->GetUID()) continue;
            partpair->SetParticles(part1,part2);
           
            track2= trackEvent->GetParticle(k);
            trackpair->SetParticles(track1,track2);

            if(fPairCut->Pass(partpair) ) //check pair cut 
              { //do not meets crietria of the pair cut, try with swapped pairs
                if( fPairCut->Pass(partpair->GetSwapedPair()) )
                  continue; //swaped pairs do not meet criteria of pair cut as well, take next particle
                else 
                 { //swaped pair meets all the criteria
                   tmppartpair = partpair->GetSwapedPair();
                   tmptrackpair = trackpair->GetSwapedPair();
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
         
         partbuffer.ResetIter();
         trackbuffer.ResetIter();
         Int_t m = 0;
         while (( partEvent2 = partbuffer.Next() ))
          {
            trackEvent2 = trackbuffer.Next();
            
            m++;
            if ( (j%fDisplayMixingInfo) == 0)
               Info("ProcessTracksAndParticles",
                    "Mixing particle %d from event %d with particles from event %d",j,i,i-m);
         
            for(Int_t l = 0; l<partEvent2->GetNumberOfParticles();l++)   //  ... on all particles
              {
                part2= partEvent2->GetParticle(l);
                partpair->SetParticles(part1,part2);

                track2= trackEvent2->GetParticle(l);
                trackpair->SetParticles(track1,track2);

                if( fPairCut->Pass(partpair) ) //check pair cut
                  { //do not meets crietria of the 
                    if( fPairCut->Pass(partpair->GetSwapedPair()) )
                      continue;
                    else 
                     {
                       tmppartpair = partpair->GetSwapedPair();
                       tmptrackpair = trackpair->GetSwapedPair();
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
      partEvent1  = partbuffer.Push(partEvent1);
      trackEvent1 = trackbuffer.Push(trackEvent1);
     //end of loop over events  
    }//while (fReader->Next() == kFALSE)
} 
/*************************************************************************************/

void AliHBTAnalysis::ProcessTracks()
{
//In order to minimize calling AliRun::GetEvent (we need at one time particles from different events),
//the loops are splited
  AliHBTParticle * track1, * track2;
  AliHBTEvent * trackEvent;
  AliHBTEvent * trackEvent1 = 0x0;
  AliHBTEvent * trackEvent2;

  register UInt_t ii;

  AliHBTPair * trackpair = new AliHBTPair();
  AliHBTPair * tmptrackpair; //temporary pointer 
  
  AliHBTEventBuffer trackbuffer(fBufferSize);
  
  fReader->Rewind();
  Int_t i = -1;
  while (fReader->Next() == kFALSE)
    {
      i++;
      trackEvent = fReader->GetTrackEvent();
      if (!trackEvent) continue;
      
      if(trackEvent1 == 0x0) 
       {
         trackEvent1 = new AliHBTEvent();
         trackEvent1->SetOwner(kTRUE);
       }
      else
       {
         trackEvent1->Reset();
       }
      
      for (Int_t j = 0; j<trackEvent->GetNumberOfParticles() ; j++)
       {
         /***************************************/
         /******   Looping same events   ********/
         /******   filling numerators    ********/
         /***************************************/
         if ( (j%fDisplayMixingInfo) == 0) 
               Info("ProcessTracks",
                    "Mixing particle %d from event %d with particles from event %d",j,i,i);
         
         track1= trackEvent->GetParticle(j);
         Bool_t firstcut = fPairCut->GetFirstPartCut()->Pass(track1);
         
         if ( (firstcut == kFALSE) || (fPairCut->GetSecondPartCut()->Pass(track1) == kFALSE) )
          {
            //accepted by any cut
            // we have to copy because reader keeps only one event
            trackEvent1->AddParticle(new AliHBTParticle(*track1));
          }

         if (firstcut) continue;

         for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
           fTrackMonitorFunctions[ii]->Process(track1);

         if ( fNTrackFunctions ==0 ) continue; 
        
         for (Int_t k =j+1; k < trackEvent->GetNumberOfParticles() ; k++)
          {
           track2= trackEvent->GetParticle(k);
           if (track1->GetUID() == track2->GetUID()) continue;

           trackpair->SetParticles(track1,track2);
           if(fPairCut->Pass(trackpair)) //check pair cut
            { //do not meets crietria of the 
              if( fPairCut->Pass(trackpair->GetSwapedPair()) ) continue;
              else tmptrackpair = trackpair->GetSwapedPair();
            }
           else
            {
              tmptrackpair = trackpair;
            }
           for(ii = 0;ii<fNTrackFunctions;ii++)
               fTrackFunctions[ii]->ProcessSameEventParticles(tmptrackpair);
          }
         /***************************************/
         /***** Filling denominators    *********/
         /***************************************/
          
         if (fBufferSize == 0) continue;
         
         trackbuffer.ResetIter();
         
         Int_t m = 0;
         while (( trackEvent2 = trackbuffer.Next() ))
          {
            m++;
            if ( (j%fDisplayMixingInfo) == 0)
               Info("ProcessTracks",
                    "Mixing track %d from event %d with tracks from event %d",j,i,i-m);
            for(Int_t l = 0; l<trackEvent2->GetNumberOfParticles();l++)   //  ... on all particles
              {

                track2= trackEvent2->GetParticle(l);
                trackpair->SetParticles(track1,track2);

                if( fPairCut->Pass(trackpair) ) //check pair cut
                  { //do not meets crietria of the 
                    if( fPairCut->Pass(trackpair->GetSwapedPair()) )
                      continue;
                    else 
                     {
                       tmptrackpair = trackpair->GetSwapedPair();
                     }
                  }
                else
                 {//meets criteria of the pair cut
                  tmptrackpair = trackpair;
                 }
                 
                for(ii = 0;ii<fNTrackFunctions;ii++)
                  fTrackFunctions[ii]->ProcessDiffEventParticles(tmptrackpair);
                 
             }//for(Int_t l = 0; l<N2;l++)   //  ... on all particles
            
          }
       }
      trackEvent1 = trackbuffer.Push(trackEvent1); 
    }//while (fReader->Next() == kFALSE)
}

/*************************************************************************************/

void AliHBTAnalysis::ProcessParticles()
{
//In order to minimize calling AliRun::GetEvent (we need at one time particles from different events),
//the loops are splited
  AliHBTParticle * part1, * part2;
  AliHBTEvent * partEvent;
  AliHBTEvent * partEvent1 = 0x0;
  AliHBTEvent * partEvent2;

  register UInt_t ii;

  AliHBTPair * partpair = new AliHBTPair();
  AliHBTPair * tmppartpair; //temporary pointer 
  
  AliHBTEventBuffer partbuffer(fBufferSize);
  
  fReader->Rewind();
  Int_t i = -1;
  while (fReader->Next() == kFALSE)
    {
      i++;
      partEvent = fReader->GetParticleEvent();
      if (!partEvent) continue;
      
      if(partEvent1 == 0x0) 
       {
         partEvent1 = new AliHBTEvent();
         partEvent1->SetOwner(kTRUE);
       }
      else
       {
         partEvent1->Reset();
       }
      
      for (Int_t j = 0; j<partEvent->GetNumberOfParticles() ; j++)
       {
         /***************************************/
         /******   Looping same events   ********/
         /******   filling numerators    ********/
         /***************************************/
         if ( (j%fDisplayMixingInfo) == 0) 
               Info("ProcessParticles",
                    "Mixing particle %d from event %d with particles from event %d",j,i,i);
         
         part1 = partEvent->GetParticle(j);
         Bool_t firstcut = fPairCut->GetFirstPartCut()->Pass(part1);
         
         if ( (firstcut == kFALSE) || (fPairCut->GetSecondPartCut()->Pass(part1) == kFALSE) )
          {
            //accepted by any cut
            // we have to copy because reader keeps only one event
            partEvent1->AddParticle(new AliHBTParticle(*part1));
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
              if( fPairCut->Pass(partpair->GetSwapedPair()) ) continue;
              else tmppartpair = partpair->GetSwapedPair();
            }
           else
            {
              tmppartpair = partpair;
            }
           for(ii = 0;ii<fNParticleFunctions;ii++)
               fParticleFunctions[ii]->ProcessSameEventParticles(tmppartpair);
          }
         /***************************************/
         /***** Filling denominators    *********/
         /***************************************/
          
         if (fBufferSize == 0) continue;
         
         partbuffer.ResetIter();
         
         Int_t m = 0;
         while (( partEvent2 = partbuffer.Next() ))
          {
            m++;
            if ( (j%fDisplayMixingInfo) == 0)
               Info("ProcessParticles",
                    "Mixing particle %d from event %d with particles from event %d",j,i,i-m);
            for(Int_t l = 0; l<partEvent2->GetNumberOfParticles();l++)   //  ... on all particles
              {

                part2= partEvent2->GetParticle(l);
                partpair->SetParticles(part1,part2);

                if( fPairCut->Pass(partpair) ) //check pair cut
                  { //do not meets crietria of the 
                    if( fPairCut->Pass(partpair->GetSwapedPair()) )
                      continue;
                    else 
                     {
                       tmppartpair = partpair->GetSwapedPair();
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
      partEvent1 = partbuffer.Push(partEvent1); 
    }//while (fReader->Next() == kFALSE)
}
/*************************************************************************************/

void AliHBTAnalysis::WriteFunctions()
{
//Calls Write for all defined functions in analysis
//== writes all results
  UInt_t ii;
  for(ii = 0;ii<fNParticleFunctions;ii++)
    fParticleFunctions[ii]->Write();
                 
  for(ii = 0;ii<fNTrackFunctions;ii++)
    fTrackFunctions[ii]->Write();
                 
  for(ii = 0;ii<fNParticleAndTrackFunctions;ii++)
   fParticleAndTrackFunctions[ii]->Write();

  for(ii = 0;ii<fNParticleMonitorFunctions;ii++)
    fParticleMonitorFunctions[ii]->Write();

  for(ii = 0;ii<fNTrackMonitorFunctions;ii++)
    fTrackMonitorFunctions[ii]->Write();

  for(ii = 0;ii<fNParticleAndTrackMonitorFunctions;ii++)
   fParticleAndTrackMonitorFunctions[ii]->Write();
}
/*************************************************************************************/

void AliHBTAnalysis::SetGlobalPairCut(AliHBTPairCut* cut)
{
//Sets the global cut
  if (cut == 0x0)
   {
     Error("AliHBTAnalysis::SetGlobalPairCut","Pointer is NULL. Ignoring");
   }
  delete fPairCut;
  fPairCut = (AliHBTPairCut*)cut->Clone();
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
 Int_t i;  
 Info("RunCoherencyCheck","Checking HBT Runs Coherency");

 Info("RunCoherencyCheck","Number of events ...");
 if (fReader->GetNumberOfPartEvents() == fReader->GetNumberOfTrackEvents() ) //check whether there is the same  number of events
  {
    Info("RunCoherencyCheck","OK. %d found\n",fReader->GetNumberOfTrackEvents());
  }
 else
  { //if not the same - ERROR
   Error("AliHBTAnalysis::RunCoherencyCheck()",
         "Number of simulated events (%d) is not equal to number of reconstructed events(%d)",
         fReader->GetNumberOfPartEvents(),fReader->GetNumberOfTrackEvents());
   return kTRUE;
  }
 
 Info("RunCoherencyCheck","Checking number of Particles AND Particles Types in each event ...");
 
 AliHBTEvent *partEvent;
 AliHBTEvent *trackEvent;
 for( i = 0; i<fReader->GetNumberOfTrackEvents();i++)
  {
    partEvent= fReader->GetParticleEvent(i); //gets the "ith" event 
    trackEvent = fReader->GetTrackEvent(i);
    
    if ( (partEvent == 0x0) && (partEvent == 0x0) ) continue;
    if ( (partEvent == 0x0) || (partEvent == 0x0) )
     {
       Error("AliHBTAnalysis::RunCoherencyCheck()",
             "One event is NULL and the other one not. Event Number %d",i);
       return kTRUE;    
     }
    
    if ( partEvent->GetNumberOfParticles() != trackEvent->GetNumberOfParticles() )
     {
       Error("AliHBTAnalysis::RunCoherencyCheck()",
             "Event %d: Number of simulated particles (%d) not equal to number of reconstructed tracks (%d)",
             i,partEvent->GetNumberOfParticles() , trackEvent->GetNumberOfParticles());
       return kTRUE;
     }
    else
     for (Int_t j = 0; j<partEvent->GetNumberOfParticles(); j++)
      {
        if( partEvent->GetParticle(j)->GetPdgCode() != trackEvent->GetParticle(j)->GetPdgCode() )
         {
           Error("AliHBTAnalysis::RunCoherencyCheck()",
                 "Event %d: Particle %d: PID of simulated particle (%d) not the same of reconstructed track (%d)",
                 i,j, partEvent->GetParticle(j)->GetPdgCode(),trackEvent->GetParticle(j)->GetPdgCode() );
           return kTRUE;
           
         }
      }
  }
 Info("RunCoherencyCheck","  Done");
 Info("RunCoherencyCheck","  Everything looks OK");
 return kFALSE;
}

/*************************************************************************************/  
 
void AliHBTAnalysis::ProcessTracksAndParticlesNonIdentAnal()
{
//Performs analysis for both, tracks and particles

  AliHBTParticle * part1, * part2;
  AliHBTParticle * track1, * track2;

  AliHBTEvent * trackEvent1=0x0,*partEvent1=0x0;
  AliHBTEvent * trackEvent2=0x0,*partEvent2=0x0;
  AliHBTEvent * trackEvent3=0x0,*partEvent3=0x0;

  AliHBTEvent * rawtrackEvent, *rawpartEvent;//this we get from Reader
  
  AliHBTPair * trackpair = new AliHBTPair();
  AliHBTPair * partpair = new AliHBTPair();

  AliHBTEventBuffer partbuffer(fBufferSize);
  AliHBTEventBuffer trackbuffer(fBufferSize);
  
  register UInt_t ii;

  trackEvent1 = new AliHBTEvent();
  partEvent1 = new AliHBTEvent();
  trackEvent1->SetOwner(kFALSE);
  partEvent1->SetOwner(kFALSE);;
  
  Bool_t nocorrfctns = (fNParticleFunctions == 0) && (fNTrackFunctions == 0) && (fNParticleAndTrackFunctions == 0);
  
  fReader->Rewind();
  
  Info("ProcessTracksAndParticlesNonIdentAnal","**************************************");
  Info("ProcessTracksAndParticlesNonIdentAnal","*****  NON IDENT MODE ****************");
  Info("ProcessTracksAndParticlesNonIdentAnal","**************************************");
  
  for (Int_t i = 0;;i++)//infinite loop
    {
      if (fReader->Next()) break; //end when no more events available
      
      rawpartEvent  = fReader->GetParticleEvent();
      rawtrackEvent = fReader->GetTrackEvent();
      if ( (rawpartEvent == 0x0) || (rawtrackEvent == 0x0) ) continue;//in case of any error
      
      if ( rawpartEvent->GetNumberOfParticles() != rawtrackEvent->GetNumberOfParticles() )
       {
         Fatal("ProcessTracksAndParticlesNonIdentAnal",
               "Event %d: Number of simulated particles (%d) not equal to number of reconstructed tracks (%d)",
               i,rawpartEvent->GetNumberOfParticles() , rawtrackEvent->GetNumberOfParticles());
         return;
       }
      
      /********************************/
      /*      Filtering out           */
      /********************************/
      if ( ( (partEvent2==0x0) || (trackEvent2==0x0)) )//in case fBufferSize == 0 and pointers are created do not eneter
       {
         partEvent2  = new AliHBTEvent();
         trackEvent2 = new AliHBTEvent();
         partEvent2->SetOwner(kTRUE);
         trackEvent2->SetOwner(kTRUE);
       }
       
      FilterOut(partEvent1, partEvent2, rawpartEvent, trackEvent1, trackEvent2, rawtrackEvent);
      
      for (Int_t j = 0; j<partEvent1->GetNumberOfParticles() ; j++)
       {
         if ( (j%fDisplayMixingInfo) == 0) 
            Info("ProcessTracksAndParticlesNonIdentAnal",
                 "Mixing particle %d from event %d with particles from event %d",j,i,i);

         part1= partEvent1->GetParticle(j);
         track1= trackEvent1->GetParticle(j);
  
//PID reconstruction imperfections
//         if( part1->GetPdgCode() != track1->GetPdgCode() )
//          {
//            Fatal("ProcessTracksAndParticlesNonIdentAnal",
//                  "Event %d: Particle %d: PID of simulated particle (%d) not the same of reconstructed track (%d)",
//                  i,j, part1->GetPdgCode(),track1->GetPdgCode() );
//            return;
//          }
         
         for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
           fParticleMonitorFunctions[ii]->Process(part1);
         for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
           fTrackMonitorFunctions[ii]->Process(track1);
         for(ii = 0; ii<fNParticleAndTrackMonitorFunctions; ii++)
           fParticleAndTrackMonitorFunctions[ii]->Process(track1,part1);

         if (nocorrfctns) continue;

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

            if( (fPairCut->PassPairProp(partpair)) ) //check pair cut
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
     partbuffer.ResetIter();
     trackbuffer.ResetIter();
     
     Int_t nmonitor = 0;
     
     while ( (partEvent3 = partbuffer.Next() ) != 0x0)
      {
        trackEvent3 = trackbuffer.Next();
        
        if ( (j%fDisplayMixingInfo) == 0) 
          Info("ProcessTracksAndParticlesNonIdentAnal",
               "Mixing particle %d from event %d with particles from event %d",j,i,i-(++nmonitor));
        
        for (Int_t k = 0; k < partEvent3->GetNumberOfParticles() ; k++)
          {
            part2= partEvent3->GetParticle(k);
            partpair->SetParticles(part1,part2);

            track2= trackEvent3->GetParticle(k);
            trackpair->SetParticles(track1,track2);

            if( (fPairCut->PassPairProp(partpair)) ) //check pair cut
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
     
     partEvent2 = partbuffer.Push(partEvent2);
     trackEvent2 = trackbuffer.Push(trackEvent2);

  }//end of loop over events (1)
  
 partbuffer.SetOwner(kTRUE);
 trackbuffer.SetOwner(kTRUE);
 
 delete partEvent1;
 delete trackEvent1;
 delete partEvent2;
 delete trackEvent2;
 delete partpair;
 delete trackpair;

}
/*************************************************************************************/  
 
void AliHBTAnalysis::ProcessTracksNonIdentAnal()
{
//Process Tracks only with non identical mode
  AliHBTParticle * track1, * track2;

  AliHBTEvent * trackEvent1=0x0;
  AliHBTEvent * trackEvent2=0x0;
  AliHBTEvent * trackEvent3=0x0;

  AliHBTEvent * rawtrackEvent;
  
  AliHBTPair * trackpair = new AliHBTPair();
  AliHBTEventBuffer trackbuffer(fBufferSize);

  register UInt_t ii;

  trackEvent1 = new AliHBTEvent();
  trackEvent1->SetOwner(kFALSE);
  
  fReader->Rewind();
  
  Info("ProcessTracksNonIdentAnal","**************************************");
  Info("ProcessTracksNonIdentAnal","*****  NON IDENT MODE ****************");
  Info("ProcessTracksNonIdentAnal","**************************************");
  
  
  for (Int_t i = 0;;i++)//infinite loop
    {
      if (fReader->Next()) break; //end when no more events available
      rawtrackEvent = fReader->GetTrackEvent();
      
      if (rawtrackEvent == 0x0)  continue;//in case of any error

      /********************************/
      /*      Filtering out           */
      /********************************/
      if ( trackEvent2==0x0 )//in case fBufferSize == 0 and pointers are created do not eneter
       {
         trackEvent2 = new AliHBTEvent();
         trackEvent2->SetOwner(kTRUE);
       }
       
      FilterOut(trackEvent1, trackEvent2, rawtrackEvent);
      
      for (Int_t j = 0; j<trackEvent1->GetNumberOfParticles() ; j++)
       {
         if ( (j%fDisplayMixingInfo) == 0) 
           Info("ProcessTracksNonIdentAnal",
                "Mixing particle %d from event %d with particles from event %d",j,i,i);

         track1= trackEvent1->GetParticle(j);

         for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
           fTrackMonitorFunctions[ii]->Process(track1);
         
         if (fNTrackFunctions == 0x0) continue;
         
         /***************************************/
         /******   filling numerators    ********/
         /****** (particles from event2) ********/
         /***************************************/
         for (Int_t k = 0; k < trackEvent2->GetNumberOfParticles() ; k++)
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
     trackbuffer.ResetIter();
     Int_t nmonitor = 0;
     
     while ( (trackEvent3 = trackbuffer.Next() ) != 0x0)
      {
        
        if ( (j%fDisplayMixingInfo) == 0) 
            Info("ProcessTracksNonIdentAnal",
                 "Mixing particle %d from event %d with particles from event %d",j,i,i-(++nmonitor));
        
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
     
     trackEvent2 = trackbuffer.Push(trackEvent2);

  }//end of loop over events (1)

 trackbuffer.SetOwner(kTRUE);
 delete trackpair;
 delete trackEvent1;
 delete trackEvent2;
}
/*************************************************************************************/  

void AliHBTAnalysis::ProcessParticlesNonIdentAnal()
{
//process paricles only with non identical mode
  AliHBTParticle * part1 = 0x0, * part2 = 0x0;

  AliHBTEvent * partEvent1 = 0x0;
  AliHBTEvent * partEvent2 = 0x0;
  AliHBTEvent * partEvent3 = 0x0;

  AliHBTEvent * rawpartEvent = 0x0;

  AliHBTPair * partpair = new AliHBTPair();
  AliHBTEventBuffer partbuffer(fBufferSize);

  register UInt_t ii;

  partEvent1 = new AliHBTEvent();
  partEvent1->SetOwner(kFALSE);
  
  fReader->Rewind();
  
  Info("ProcessParticlesNonIdentAnal","**************************************");
  Info("ProcessParticlesNonIdentAnal","*****  NON IDENT MODE ****************");
  Info("ProcessParticlesNonIdentAnal","**************************************");

  for (Int_t i = 0;;i++)//infinite loop
    {
      if (fReader->Next()) break; //end when no more events available
      
      rawpartEvent  = fReader->GetParticleEvent();
      if ( rawpartEvent == 0x0  ) continue;//in case of any error

      /********************************/
      /*      Filtering out           */
      /********************************/
      if (partEvent2==0x0)//in case fBufferSize == 0 and pointers are created do not eneter
       {
         partEvent2  = new AliHBTEvent();
         partEvent2->SetOwner(kTRUE);
       }
       
      FilterOut(partEvent1, partEvent2, rawpartEvent);
      
      for (Int_t j = 0; j<partEvent1->GetNumberOfParticles() ; j++)
       {
         if ( (j%fDisplayMixingInfo) == 0) 
            Info("ProcessParticlesNonIdentAnal",
                 "Mixing particle %d from event %d with particles from event %d",j,i,i);

         part1= partEvent1->GetParticle(j);

         for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
           fParticleMonitorFunctions[ii]->Process(part1);
         
         if (fNParticleFunctions == 0) continue;
         
         /***************************************/
         /******   filling numerators    ********/
         /****** (particles from event2) ********/
         /***************************************/
         for (Int_t k = 0; k < partEvent2->GetNumberOfParticles() ; k++)
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
     partbuffer.ResetIter();
     Int_t nmonitor = 0;

     while ( (partEvent3 = partbuffer.Next() ) != 0x0)
      {
        if ( (j%fDisplayMixingInfo) == 0) 
            Info("ProcessParticlesNonIdentAnal",
                 "Mixing particle %d from event %d with particles from event %d",j,i,i-(++nmonitor));

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
    partEvent2 = partbuffer.Push(partEvent2);
  }//end of loop over events (1)
  
 partbuffer.SetOwner(kTRUE);
 delete partpair;
 delete partEvent1;
 delete partEvent2;
}

/*************************************************************************************/  
void AliHBTAnalysis::FilterOut(AliHBTEvent* outpart1, AliHBTEvent* outpart2, AliHBTEvent* inpart,
                     AliHBTEvent* outtrack1, AliHBTEvent* outtrack2, AliHBTEvent* intrack)
{
 //Puts particles accepted as a first particle by global cut in out1
 //and as a second particle in out2

  AliHBTParticle* part, *track;

  outpart1->Reset();
  outpart2->Reset();
  outtrack1->Reset();
  outtrack2->Reset();
  
  AliHBTParticleCut *cut1 = fPairCut->GetFirstPartCut();
  AliHBTParticleCut *cut2 = fPairCut->GetSecondPartCut();
  
  Bool_t in1, in2;
  
  for (Int_t i = 0; i < inpart->GetNumberOfParticles(); i++)
   {
     in1 = in2 = kTRUE;
     part = inpart->GetParticle(i);
     track = intrack->GetParticle(i);
     
     if ( (cut1->Pass(part))  ) in1 = kFALSE; //if part  is rejected by cut1, in1 is false
     if ( (cut2->Pass(part))  ) in2 = kFALSE; //if part  is rejected by cut2, in2 is false
     
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
        outpart2->AddParticle(new AliHBTParticle(*part));
        outtrack2->AddParticle(new AliHBTParticle(*track));
        continue;
      }
   }
}
/*************************************************************************************/  

void AliHBTAnalysis::FilterOut(AliHBTEvent* out1, AliHBTEvent* out2, AliHBTEvent* in)
{
 //Puts particles accepted as a first particle by global cut in out1
 //and as a second particle in out2
  AliHBTParticle* part;
  
  out1->Reset();
  out2->Reset();
  
  AliHBTParticleCut *cut1 = fPairCut->GetFirstPartCut();
  AliHBTParticleCut *cut2 = fPairCut->GetSecondPartCut();
  
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
