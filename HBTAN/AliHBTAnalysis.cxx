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
 
#include <TBenchmark.h>
#include <TList.h>

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
   Fatal("AliHBTAnalysis(const AliHBTAnalysis&)","Sensless");
 }
/*************************************************************************************/ 
const AliHBTAnalysis& AliHBTAnalysis::operator=(const AliHBTAnalysis& right)
 {
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
 UInt_t ii;
 for(ii = 0;ii<fNParticleFunctions;ii++)
   delete fParticleFunctions[ii];
                
 for(ii = 0;ii<fNTrackFunctions;ii++)
   delete fTrackFunctions[ii];

 for(ii = 0;ii<fNParticleAndTrackFunctions;ii++)
   delete fParticleAndTrackFunctions[ii];
 
 for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
   delete fParticleMonitorFunctions[ii];
   
 for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
   delete fTrackMonitorFunctions[ii];
   
 for(ii = 0; ii<fNParticleAndTrackMonitorFunctions; ii++)
   delete fParticleAndTrackMonitorFunctions[ii];
}
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
 
 if(oT && oP)
  { 
    if (fReader->GetNumberOfPartEvents() <1)
     {
       Error("Process","There is no Particles. Maybe change the option?");
       return;
     }
    if (fReader->GetNumberOfTrackEvents() <1)
     {
       Error("Process","There is no Tracks. Maybe change the option?");
       return;
     }
    
    if ( RunCoherencyCheck() )
      {
        Error("Process",
              "Coherency check not passed. Maybe change the option?\n");
        return;
      }
    if (nonid) ProcessTracksAndParticlesNonIdentAnal();
    else ProcessTracksAndParticles();
    return;
  }

 if(oT)
  {
    if (fReader->GetNumberOfTrackEvents() <1)
     {
       Error("Process","There is no data to analyze.");
       return;
     }
    if (nonid) ProcessTracksNonIdentAnal();
    else ProcessTracks();
    return;
  }
 
 if(oP)
  {
    if (fReader->GetNumberOfPartEvents() <1)
     {
       Error("Process","There is no data to analyze.");
       return;
     }
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
  AliHBTEvent * trackEvent2,*partEvent2;
  
//  Int_t N1, N2, N=0; //number of particles in current event(we prcess two events in one time)
  
  Int_t Nev = fReader->GetNumberOfTrackEvents();
  
  /***************************************/
  /******   Looping same events   ********/
  /******   filling numerators    ********/
  /***************************************/
  AliHBTPair * trackpair = new AliHBTPair();
  AliHBTPair * partpair = new AliHBTPair();

  AliHBTPair * tmptrackpair;//temprary pointers to pairs
  AliHBTPair * tmppartpair;

  register UInt_t ii;
  
  for (Int_t i = 0;i<Nev;i++)
    {
      partEvent= fReader->GetParticleEvent(i);
      trackEvent = fReader->GetTrackEvent(i);
      
      if (!partEvent) continue;
      
      //N = 0;
      
      for (Int_t j = 0; j<partEvent->GetNumberOfParticles() ; j++)
       {

         if ( (j%fDisplayMixingInfo) == 0) 
            Info("ProcessTracksAndParticles",
                 "Mixing particle %d from event %d with particles from event %d",j,i,i);

         part1= partEvent->GetParticle(j);
         track1= trackEvent->GetParticle(j);       

         if (fPairCut->GetFirstPartCut()->Pass(part1)) continue;

         for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
           fParticleMonitorFunctions[ii]->ProcessSameEventParticles(part1);
         for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
           fTrackMonitorFunctions[ii]->ProcessSameEventParticles(track1);
         for(ii = 0; ii<fNParticleAndTrackMonitorFunctions; ii++)
           fParticleAndTrackMonitorFunctions[ii]->ProcessSameEventParticles(track1,part1);

         if ( (fNParticleFunctions == 0) && (fNTrackFunctions ==0) && (fNParticleAndTrackFunctions == 0))
           continue; 
        
         for (Int_t k =j+1; k < partEvent->GetNumberOfParticles() ; k++)
          {
            part2= partEvent->GetParticle(k);
            partpair->SetParticles(part1,part2);
           
            track2= trackEvent->GetParticle(k);       
            trackpair->SetParticles(track1,track2);

            if(fPairCut->Pass(partpair) ) //check pair cut
              { //do not meets crietria of the pair cut, try with swaped pairs
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
           }
       }
    }

  /***************************************/
  /***** Filling denominators    *********/
  /***************************************/
  if (fBufferSize != 0)
   for (Int_t i = 0;i<Nev-1;i++)   //In each event (but last) ....
    {
  
      if ((fNParticleFunctions == 0) && (fNTrackFunctions ==0) && (fNParticleAndTrackFunctions == 0))
        continue;  

      partEvent= fReader->GetParticleEvent(i);
      if (!partEvent) continue;
      
      trackEvent = fReader->GetTrackEvent(i); 
      
      for (Int_t j = 0; j< partEvent->GetNumberOfParticles(); j++) // ... Loop over all particles ...
       {
           part1= partEvent->GetParticle(j);
           track1= trackEvent->GetParticle(j);

           if (fPairCut->GetFirstPartCut()->Pass(part1)) continue;
 
           Int_t NNN;
           
           if ( ((i+fBufferSize) >= Nev) ||( fBufferSize < 0) ) //if buffer size is negative 
                                                                //or current event+buffersize is greater
                                                                //than max nuber of events
             {
               NNN = Nev; //set the max event number 
             }
           else 
             {
               NNN = i+fBufferSize; //set the current event number + fBufferSize
             }
 
           for (Int_t k = i+1; k<NNN;k++)  // ... Loop over next event
            {
             
             partEvent2= fReader->GetParticleEvent(k);
             if (!partEvent2) continue;
              
             trackEvent2 = fReader->GetTrackEvent(k);
             
             if ( (j%fDisplayMixingInfo) == 0) 
                 Info("ProcessTracksAndParticles",
                      "Mixing particle %d from event %d with particles from event %d",j,i,k);
            
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
            }//for (Int_t k = i+1; k<NNN;k++)  // ... Loop over next event
       }
    } 
  /***************************************/
} 
/*************************************************************************************/

void AliHBTAnalysis::ProcessTracks()
{
  //In order to minimize calling AliRun::GetEvent (we need at one time particles from different events),
//the loops are splited
  AliHBTParticle * track1, * track2;
  AliHBTEvent * trackEvent;
  AliHBTEvent * trackEvent2;

  UInt_t ii;
  Int_t Nev = fReader->GetNumberOfTrackEvents();

  AliHBTPair * trackpair = new AliHBTPair();
  AliHBTPair * tmptrackpair; //temporary pointer 
  
  /***************************************/
  /******   Looping same events   ********/
  /******   filling numerators    ********/
  /***************************************/
  for (Int_t i = 0;i<Nev;i++)
    {
      trackEvent = fReader->GetTrackEvent(i);
      if (!trackEvent) continue;
      
      for (Int_t j = 0; j<trackEvent->GetNumberOfParticles() ; j++)
       {
         if ( (j%fDisplayMixingInfo) == 0) 
               Info("ProcessTracks",
                    "Mixing particle %d from event %d with particles from event %d",j,i,i);
         
         track1= trackEvent->GetParticle(j);       
         if (fPairCut->GetFirstPartCut()->Pass(track1)) continue;

         for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
           fTrackMonitorFunctions[ii]->ProcessSameEventParticles(track1);

         if ( fNTrackFunctions ==0 )
           continue; 
        
         for (Int_t k =j+1; k < trackEvent->GetNumberOfParticles() ; k++)
          {
           track2= trackEvent->GetParticle(k);       
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
       }
    }

  /***************************************/
  /***** Filling diff histogram *********/
  /***************************************/
  if (fBufferSize != 0)
   for (Int_t i = 0;i<Nev-1;i++)   //In each event (but last) ....
    {
      if ( fNTrackFunctions ==0 )
        continue; 

      trackEvent = fReader->GetTrackEvent(i);
      if (!trackEvent) continue;

      for (Int_t j = 0; j< trackEvent->GetNumberOfParticles(); j++) // ... Loop over all particles ...
       {
         track1= trackEvent->GetParticle(j);
         if (fPairCut->GetFirstPartCut()->Pass(track1)) continue;

         Int_t NNN;
           
         if ( ((i+fBufferSize) >= Nev) ||( fBufferSize < 0) ) //if buffer size is negative 
                                                              //or current event+buffersize is greater
                                                              //than max nuber of events
          {
            NNN = Nev; //set the max event number 
          }
         else 
          {
            NNN = i+fBufferSize; //set the current event number + fBufferSize
          }
 
         for (Int_t k = i+1; k<NNN;k++)  // ... Loop over next event
          {
            trackEvent2 = fReader->GetTrackEvent(k);
            if (!trackEvent2) continue;
             
            if ( (j%fDisplayMixingInfo) == 0) 
                 Info("ProcessTracksAndParticles",
                      "Mixing particle %d from event %d with particles from event %d",j,i,k);
            
            for(Int_t l = 0; l<trackEvent2->GetNumberOfParticles();l++)   //  ... on all particles
             {
               track2= trackEvent2->GetParticle(l);
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
                   fTrackFunctions[ii]->ProcessDiffEventParticles(tmptrackpair);
               
             }//for(Int_t l = 0; l<N2;l++)   //  ... on all particles
          }//for (Int_t k = i+1; k<NNN;k++)  // ... Loop over next event
       }
    } 
  /***************************************/
}

/*************************************************************************************/

void AliHBTAnalysis::ProcessParticles()
{
//In order to minimize calling AliRun::GetEvent (we need at one time particles from different events),
//the loops are splited
  AliHBTParticle * part1, * part2;
  
  AliHBTEvent * partEvent;
  AliHBTEvent * partEvent2;

  AliHBTPair * partpair = new AliHBTPair();
  AliHBTPair * tmppartpair; //temporary pointer to the pair
  
  Int_t Nev = fReader->GetNumberOfPartEvents();
  
  /***************************************/
  /******   Looping same events   ********/
  /******   filling numerators    ********/
  /***************************************/
  for (Int_t i = 0;i<Nev;i++)
    {
      partEvent= fReader->GetParticleEvent(i);
      if (!partEvent) continue;
      //N = 0;
      
      for (Int_t j = 0; j<partEvent->GetNumberOfParticles() ; j++)
       {
         if ( (j%fDisplayMixingInfo) == 0) 
                 Info("ProcessParticles",
                      "Mixing particle %d from event %d with particles from event %d",j,i,i);

         part1= partEvent->GetParticle(j);
         if (fPairCut->GetFirstPartCut()->Pass(part1)) continue;
        
         UInt_t zz;
         for(zz = 0; zz<fNParticleMonitorFunctions; zz++)
           fParticleMonitorFunctions[zz]->ProcessSameEventParticles(part1);

         if ( fNParticleFunctions ==0 )
           continue; 

         for (Int_t k =j+1; k < partEvent->GetNumberOfParticles() ; k++)
          {
            part2= partEvent->GetParticle(k);
            partpair->SetParticles(part1,part2);
            
            if( fPairCut->Pass(partpair) ) //check pair cut
              { //do not meets crietria of the pair cut, try with swaped pairs
                if(  fPairCut->Pass(partpair->GetSwapedPair() )  ) 
                  continue; //swaped pairs do not meet criteria of pair cut as well, take next particle
                else 
                 { //swaped pair meets all the criteria
                   tmppartpair = partpair->GetSwapedPair();
                 }
              }
            else
              {
                tmppartpair = partpair;
              }
            
            UInt_t ii;
            for(ii = 0;ii<fNParticleFunctions;ii++)
                   fParticleFunctions[ii]->ProcessSameEventParticles(tmppartpair);
           }
       }
    }

  /***************************************/
  /***** Filling diff histogram *********/
  /***************************************/
  if (fBufferSize != 0) //less then 0 mix everything, == 0 do not mix denominator
   for (Int_t i = 0;i<Nev-1;i++)   //In each event (but last)....
    {
      if ( fNParticleFunctions ==0 )
        continue; 

      partEvent= fReader->GetParticleEvent(i);
      if (!partEvent) continue;

      for (Int_t j = 0; j< partEvent->GetNumberOfParticles(); j++) // ... Loop over all particles ...
       {
           part1= partEvent->GetParticle(j);
           if (fPairCut->GetFirstPartCut()->Pass(part1)) continue;
           Int_t NNN;

           if ( ((i+fBufferSize) >= Nev) ||( fBufferSize < 0) ) //if buffer size is negative 
                                                                //or current event+buffersize is greater
                                                                //than max nuber of events
            {
             NNN = Nev; //set the max event number 
            }
           else 
            {
             NNN = i+fBufferSize; //set the current event number + fBufferSize
            }
           
           for (Int_t k = i+1; k<NNN;k++)  // ... Loop over next event
            {
             
             partEvent2= fReader->GetParticleEvent(k);
             if (!partEvent2) continue;
             
             if ( (j%fDisplayMixingInfo) == 0) 
                Info("ProcessParticles",
                     "Mixing particle %d from event %d with particles from event %d",j,i,k);
            
             for(Int_t l = 0; l<partEvent2->GetNumberOfParticles();l++)   //  ... on all particles
              {
                part2= partEvent2->GetParticle(l);
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
                UInt_t ii;
                for(ii = 0;ii<fNParticleFunctions;ii++)
                  fParticleFunctions[ii]->ProcessDiffEventParticles(tmppartpair);
               
              }//for(Int_t l = 0; l<N2;l++)   //  ... on all particles
            }//for (Int_t k = i+1; k<NNN;k++)  // ... Loop over next event
       }
    } 
  /***************************************/

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

  AliHBTEvent * rawtrackEvent, *rawpartEvent;
  
  Int_t Nev = fReader->GetNumberOfTrackEvents();

  AliHBTPair * trackpair = new AliHBTPair();
  AliHBTPair * partpair = new AliHBTPair();

  TList tbuffer;
  TList pbuffer;
  Int_t ninbuffer = 0;
  UInt_t ii;

  trackEvent1 = new AliHBTEvent();
  partEvent1 = new AliHBTEvent();
  trackEvent1->SetOwner(kFALSE);
  partEvent1->SetOwner(kFALSE);;
  
  Info("ProcessTracksAndParticlesNonIdentAnal","**************************************");
  Info("ProcessTracksAndParticlesNonIdentAnal","*****  NON IDENT MODE ****************");
  Info("ProcessTracksAndParticlesNonIdentAnal","**************************************");
  
  for (Int_t i = 0;i<Nev;i++)
    {
      rawpartEvent  = fReader->GetParticleEvent(i);
      rawtrackEvent = fReader->GetTrackEvent(i);
      if ( (rawpartEvent == 0x0) || (rawtrackEvent == 0x0) ) continue;//in case of any error

      /********************************/
      /*      Filtering out           */
      /********************************/
      if ( (fBufferSize != 0) || ( (partEvent2==0x0)||(trackEvent2==0x0)) )//in case fBufferSize == 0 and pointers are created do not eneter
       if ((ninbuffer > fBufferSize) && (fBufferSize > 0))
        {//if we have in buffer desired number of events, use the last. If fBufferSize<0 mix all events for background
         partEvent2  = (AliHBTEvent*)pbuffer.Remove(pbuffer.Last()); //remove from the end to be reset, filled and put on the beginning
         trackEvent2 = (AliHBTEvent*)tbuffer.Remove(tbuffer.Last());
         ninbuffer--;
        }
       else
        {
         partEvent2  = new AliHBTEvent();
         trackEvent2 = new AliHBTEvent();
         partEvent2->SetOwner(kFALSE);
         trackEvent2->SetOwner(kFALSE);
        }
      FilterOut(partEvent1, partEvent2, rawpartEvent, trackEvent1, trackEvent2, rawtrackEvent);
      
      for (Int_t j = 0; j<partEvent1->GetNumberOfParticles() ; j++)
       {
         if ( (j%fDisplayMixingInfo) == 0) 
            Info("ProcessTracksAndParticlesNonIdentAnal",
                 "Mixing particle %d from event %d with particles from event %d",j,i,i);

         part1= partEvent1->GetParticle(j);
         track1= trackEvent1->GetParticle(j);

         for(ii = 0; ii<fNParticleMonitorFunctions; ii++)
           fParticleMonitorFunctions[ii]->ProcessSameEventParticles(part1);
         for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
           fTrackMonitorFunctions[ii]->ProcessSameEventParticles(track1);
         for(ii = 0; ii<fNParticleAndTrackMonitorFunctions; ii++)
           fParticleAndTrackMonitorFunctions[ii]->ProcessSameEventParticles(track1,part1);

         /***************************************/
         /******   filling numerators    ********/
         /****** (particles from event2) ********/
         /***************************************/
         for (Int_t k = 0; k < partEvent2->GetNumberOfParticles() ; k++)
          {
            part2= partEvent2->GetParticle(k);
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
     TIter piter(&pbuffer);
     TIter titer(&tbuffer);
     Int_t nmonitor = 0;
     
     while ( (partEvent3 = (AliHBTEvent*)piter.Next()) != 0x0)
      {
        trackEvent3 = (AliHBTEvent*)titer.Next();
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
     
     if ( fBufferSize == 0) continue;//do not mix diff histograms-> do not add to buffer list
     pbuffer.AddFirst(partEvent2);
     tbuffer.AddFirst(trackEvent2);
     partEvent2 = 0x0;
     trackEvent2 = 0x0;
     ninbuffer++;

  }//end of loop over events (1)
 pbuffer.SetOwner();  //to clean stored events by the way of buffer destruction
 tbuffer.SetOwner();  
 delete partEvent1;
 delete trackEvent1;
 delete partpair;
 delete trackpair;
 if ( fBufferSize == 0)
  {//in that case we did not add these events to list
    delete partEvent2;
    delete trackEvent2;
  }

}
/*************************************************************************************/  
 
void AliHBTAnalysis::ProcessTracksNonIdentAnal()
{
  AliHBTParticle * track1, * track2;

  AliHBTEvent * trackEvent1=0x0;
  AliHBTEvent * trackEvent2=0x0;
  AliHBTEvent * trackEvent3=0x0;

  AliHBTEvent * rawtrackEvent;
  
  Int_t Nev = fReader->GetNumberOfTrackEvents();

  AliHBTPair * trackpair = new AliHBTPair();

  TList tbuffer;
  Int_t ninbuffer = 0;
  register UInt_t ii;

  trackEvent1 = new AliHBTEvent();
  trackEvent1->SetOwner(kFALSE);
  
  Info("ProcessTracksNonIdentAnal","**************************************");
  Info("ProcessTracksNonIdentAnal","*****  NON IDENT MODE ****************");
  Info("ProcessTracksNonIdentAnal","**************************************");
  
  for (Int_t i = 0;i<Nev;i++)
    {
      rawtrackEvent = fReader->GetTrackEvent(i);
      if (rawtrackEvent == 0x0)  continue;//in case of any error

      /********************************/
      /*      Filtering out           */
      /********************************/
      if ( (fBufferSize != 0) || ( trackEvent2==0x0) )//in case fBufferSize == 0 and pointers are created do not eneter
       if ((ninbuffer > fBufferSize) && (fBufferSize > 0))
        {//if we have in buffer desired number of events, use the last. If fBufferSize<0 mix all events for background
         trackEvent2 = (AliHBTEvent*)tbuffer.Remove(tbuffer.Last());
         ninbuffer--;
        }
       else
        {
         trackEvent2 = new AliHBTEvent();
         trackEvent2->SetOwner(kFALSE);
        }
      FilterOut(trackEvent1, trackEvent2, rawtrackEvent);
      
      for (Int_t j = 0; j<trackEvent1->GetNumberOfParticles() ; j++)
       {
         if ( (j%fDisplayMixingInfo) == 0) 
           Info("ProcessTracksNonIdentAnal",
                "Mixing particle %d from event %d with particles from event %d",j,i,i);

         track1= trackEvent1->GetParticle(j);

         for(ii = 0; ii<fNTrackMonitorFunctions; ii++)
           fTrackMonitorFunctions[ii]->ProcessSameEventParticles(track1);

         /***************************************/
         /******   filling numerators    ********/
         /****** (particles from event2) ********/
         /***************************************/
         for (Int_t k = 0; k < trackEvent2->GetNumberOfParticles() ; k++)
          {
            track2= trackEvent2->GetParticle(k);
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
     TIter titer(&tbuffer);
     Int_t nmonitor = 0;
     
     while ( (trackEvent3 = (AliHBTEvent*)titer.Next()) != 0x0)
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
     
     if ( fBufferSize == 0) continue;//do not mix diff histograms     
     tbuffer.AddFirst(trackEvent2);
     trackEvent2 = 0x0;
     ninbuffer++;

  }//end of loop over events (1)

 delete trackpair;
 delete trackEvent1;
 if (fBufferSize == 0) delete trackEvent2;
 tbuffer.SetOwner();  
}
/*************************************************************************************/  

void AliHBTAnalysis::ProcessParticlesNonIdentAnal()
{
  AliHBTParticle * part1 = 0x0, * part2 = 0x0;

  AliHBTEvent * partEvent1 = 0x0;
  AliHBTEvent * partEvent2 = 0x0;
  AliHBTEvent * partEvent3 = 0x0;

  AliHBTEvent * rawpartEvent = 0x0;

  Int_t Nev = fReader->GetNumberOfPartEvents();

  AliHBTPair * partpair = new AliHBTPair();

  TList pbuffer;
  Int_t ninbuffer = 0;

  partEvent1 = new AliHBTEvent();
  partEvent1->SetOwner(kFALSE);;
  
  Info("ProcessParticlesNonIdentAnal","**************************************");
  Info("ProcessParticlesNonIdentAnal","*****  NON IDENT MODE ****************");
  Info("ProcessParticlesNonIdentAnal","**************************************");

  for (Int_t i = 0;i<Nev;i++)
    {
      rawpartEvent  = fReader->GetParticleEvent(i);
      if ( rawpartEvent == 0x0  ) continue;//in case of any error

      /********************************/
      /*      Filtering out           */
      /********************************/
      if ( (fBufferSize != 0) || ( partEvent2==0x0) )//in case fBufferSize == 0 and pointers are created do not eneter
       if ((ninbuffer > fBufferSize) && (fBufferSize > 0))
        {//if we have in buffer desired number of events, use the last. If fBufferSize<0 mix all events for background
         partEvent2  = (AliHBTEvent*)pbuffer.Remove(pbuffer.Last()); //remove from the end to be reset, filled and put on the beginning
         ninbuffer--; 
        }
       else
        {
         partEvent2  = new AliHBTEvent();
         partEvent2->SetOwner(kFALSE);
        }
      FilterOut(partEvent1, partEvent2, rawpartEvent);
      
      for (Int_t j = 0; j<partEvent1->GetNumberOfParticles() ; j++)
       {
         if ( (j%fDisplayMixingInfo) == 0) 
            Info("ProcessParticlesNonIdentAnal",
                 "Mixing particle %d from event %d with particles from event %d",j,i,i);

         part1= partEvent1->GetParticle(j);

         UInt_t zz;
         for(zz = 0; zz<fNParticleMonitorFunctions; zz++)
           fParticleMonitorFunctions[zz]->ProcessSameEventParticles(part1);

         /***************************************/
         /******   filling numerators    ********/
         /****** (particles from event2) ********/
         /***************************************/
         for (Int_t k = 0; k < partEvent2->GetNumberOfParticles() ; k++)
          {
            part2= partEvent2->GetParticle(k);
            partpair->SetParticles(part1,part2);

            if(fPairCut->PassPairProp(partpair) ) //check pair cut
              { //do not meets crietria of the pair cut
                  continue; 
              }
            else
             {//meets criteria of the pair cut
              UInt_t ii;
              for(ii = 0;ii<fNParticleFunctions;ii++)
                {
                  fParticleFunctions[ii]->ProcessSameEventParticles(partpair);
                }
             }
           }
     if ( fBufferSize == 0) continue;//do not mix diff histograms
     /***************************************/
     /***** Filling denominators    *********/
     /***************************************/
     TIter piter(&pbuffer);
     Int_t nmonitor = 0;

     while ( (partEvent3 = (AliHBTEvent*)piter.Next()) != 0x0)
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
              UInt_t ii;
              for(ii = 0;ii<fNParticleFunctions;ii++)
               {
                 fParticleFunctions[ii]->ProcessDiffEventParticles(partpair);
               }
             }
           }// for particles event2
         }//while event2
       }//for over particles in event1
     if ( fBufferSize == 0) continue;//do not mix diff histograms
     pbuffer.AddFirst(partEvent2);
     partEvent2 = 0x0; 
     ninbuffer++;

  }//end of loop over events (1)
 delete partpair;
 delete partEvent1;
 if ( fBufferSize == 0) delete partEvent2;
 pbuffer.SetOwner();//to delete stered events.
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
        outpart2->AddParticle(part);
        outtrack2->AddParticle(track);
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
