
#include "AliHBTAnalysis.h"

#include <iostream.h>

#include "AliHBTRun.h"
#include "AliHBTReader.h"
#include "AliHBTParticle.h"
#include "AliHBTParticleCut.h"
#include "AliHBTPair.h"
#include "AliHBTPairCut.h"
#include "AliHBTFunction.h"

#include <TList.h>



ClassImp(AliHBTAnalysis)

const UInt_t AliHBTAnalysis::fgkFctnArraySize = 100;
const Int_t AliHBTAnalysis::fgkHbtAnalyzeAll = 0;

AliHBTAnalysis::AliHBTAnalysis()
 {
   fReader = 0x0;
   
   fTrackFunctions = new AliHBTOnePairFctn* [fgkFctnArraySize];
   fParticleFunctions = new AliHBTOnePairFctn* [fgkFctnArraySize];
   fParticleAndTrackFunctions = new AliHBTTwoPairFctn* [fgkFctnArraySize];
   
   fNTrackFunctions = 0;
   fNParticleFunctions = 0;
   fNParticleAndTrackFunctions = 0;
   
   fPairCut = new AliHBTEmptyPairCut();//empty cut - accepts all particles
   
 }
/*************************************************************************************/ 

AliHBTAnalysis::~AliHBTAnalysis()
 {
 //destructor
 //note that we do not delete functions itself
 // they should be deleted by whom where created
 //we only store pointers, and we use "new" only for pointers array
   delete [] fTrackFunctions;
   delete [] fParticleFunctions;
   delete [] fParticleAndTrackFunctions;
   
   delete fPairCut; // always have an copy of an object - we create we dstroy
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
    ProcessTracksAndParticles();
    return;
  }

 if(oT)
  {
    ProcessTracks();
    return;
  }
 
 if(oP)
  {
    ProcessParticles();
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
  
  for (Int_t i = 0;i<Nev;i++)
    {
      partEvent= fReader->GetParticleEvent(i);
      trackEvent = fReader->GetTrackEvent(i);
      
      if (!partEvent) continue;
      
      //N = 0;
      
      for (Int_t j = 0; j<partEvent->GetNumberOfParticles() ; j++)
       {
         if ( (j%100) == 0) cout<<"Mixing particle "<<j<<" from event "<<i<<" with particles from event "<<i<<endl;

         part1= partEvent->GetParticle(j);
         track1= trackEvent->GetParticle(j);	
	 
         for (Int_t k =j+1; k < partEvent->GetNumberOfParticles() ; k++)
          {
            part2= partEvent->GetParticle(k);
            partpair->SetParticles(part1,part2);
           
            track2= trackEvent->GetParticle(k);	
            trackpair->SetParticles(track1,track2);

            if(fPairCut->Pass(partpair) || (fPairCut->Pass(trackpair))) //check pair cut
              { //do not meets crietria of the pair cut, try with swaped pairs
                if( ( fPairCut->Pass(partpair->GetSwapedPair()) ) || ( fPairCut->Pass(trackpair->GetSwapedPair()) ) ) 
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
            UInt_t ii;
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
  for (Int_t i = 0;i<Nev;i++)   //In each event ....
    {
      
      partEvent= fReader->GetParticleEvent(i);
      if (!partEvent) continue;
      
      trackEvent = fReader->GetTrackEvent(i); 
      
//      N=0;
      
      for (Int_t j = 0; j< partEvent->GetNumberOfParticles(); j++) // ... Loop over all particles ...
       {
//         if (N>MAXCOMB) break;
           
           part1= partEvent->GetParticle(j);

           track1= trackEvent->GetParticle(j);
 
//         for (Int_t k = i+1; k<Nev;k++)  //  ... Loop over all proceeding events ...
           Int_t NNN;
  
           if ( (i+2) < Nev) NNN = i+2;
           else NNN = Nev;
 
           for (Int_t k = i+1; k<NNN;k++)  // ... Loop over next event
            {
             
              partEvent2= fReader->GetParticleEvent(k);
              if (!partEvent2) continue;
              
              trackEvent2 = fReader->GetTrackEvent(k);
             
             if ( (j%100) == 0) cout<<"Mixing particle "<<j<<" from event "<<i<<" with particles from event "<<k<<endl;
	     
             for(Int_t l = 0; l<partEvent2->GetNumberOfParticles();l++)   //  ... on all particles
              {
               
                // if (N>MAXCOMB) break;
                
                part2= partEvent2->GetParticle(l);
                partpair->SetParticles(part1,part2);

                track2= trackEvent2->GetParticle(l);
                trackpair->SetParticles(track1,track2);

                if( (fPairCut->Pass(partpair)) || (fPairCut->Pass(trackpair)) ) //check pair cut
                  { //do not meets crietria of the 
                    if( ( fPairCut->Pass(partpair->GetSwapedPair()) ) || ( fPairCut->Pass(trackpair->GetSwapedPair()) ) )
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
                UInt_t ii;
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
  
//  Int_t N1, N2, N=0; //number of particles in current event(we prcess two events in one time)
  
  Int_t Nev = fReader->GetNumberOfTrackEvents();
  
  /***************************************/
  /******   Looping same events   ********/
  /******   filling numerators    ********/
  /***************************************/
  AliHBTPair * trackpair = new AliHBTPair();
  AliHBTPair * tmptrackpair; //temporary pointer 
  
  for (Int_t i = 0;i<Nev;i++)
    {
      trackEvent = fReader->GetTrackEvent(i);
      if (!trackEvent) continue;
      //N = 0;
      
      for (Int_t j = 0; j<trackEvent->GetNumberOfParticles() ; j++)
       {
         if ( (j%100) == 0) cout<<"Mixing particle "<<j<<" from event "<<i<<" with particles from event "<<i<<endl;

         track1= trackEvent->GetParticle(j);	
	 
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
            UInt_t ii;
                
            for(ii = 0;ii<fNTrackFunctions;ii++)
                   fTrackFunctions[ii]->ProcessSameEventParticles(tmptrackpair);
                 
           
           }
       }
    }

  /***************************************/
  /***** Filling diff histogram *********/
  /***************************************/
  for (Int_t i = 0;i<Nev;i++)   //In each event ....
    {
      trackEvent = fReader->GetTrackEvent(i);
      if (!trackEvent) continue;
//      N=0;
      
      for (Int_t j = 0; j< trackEvent->GetNumberOfParticles(); j++) // ... Loop over all particles ...
       {
//         if (N>MAXCOMB) break;
           
           track1= trackEvent->GetParticle(j);
 
//         for (Int_t k = i+1; k<Nev;k++)  //  ... Loop over all proceeding events ...
           Int_t NNN;
  
           if ( (i+2) < Nev) NNN = i+2;
           else NNN = Nev;
 
           for (Int_t k = i+1; k<NNN;k++)  // ... Loop over next event
            {
             
             trackEvent2 = fReader->GetTrackEvent(k);
             if (!trackEvent2) continue;
             
             if ( (j%100) == 0) cout<<"Mixing particle "<<j<<" from event "<<i<<" with particles from event "<<k<<endl;
	     
             for(Int_t l = 0; l<trackEvent2->GetNumberOfParticles();l++)   //  ... on all particles
              {
               
                // if (N>MAXCOMB) break;
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
                UInt_t ii;
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
  
//  Int_t N1, N2, N=0; //number of particles in current event(we prcess two events in one time)
  
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
         if ( (j%100) == 0) cout<<"Mixing particle "<<j<<" from event "<<i<<" with particles from event "<<i<<endl;

         part1= partEvent->GetParticle(j);
	 
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
  for (Int_t i = 0;i<Nev;i++)   //In each event ....
    {
      partEvent= fReader->GetParticleEvent(i);
      if (!partEvent) continue;

//      N=0;
      
      for (Int_t j = 0; j< partEvent->GetNumberOfParticles(); j++) // ... Loop over all particles ...
       {
//         if (N>MAXCOMB) break;
           
           part1= partEvent->GetParticle(j);
 
//         for (Int_t k = i+1; k<Nev;k++)  //  ... Loop over all proceeding events ...
           Int_t NNN;
  
           if ( (i+2) < Nev) NNN = i+2; //loop over next event
           else NNN = Nev;
 
           for (Int_t k = i+1; k<NNN;k++)  // ... Loop over next event
            {
             
             partEvent2= fReader->GetParticleEvent(k);
             if (!partEvent2) continue;
             
             if ( (j%100) == 0) cout<<"Mixing particle "<<j<<" from event "<<i<<" with particles from event "<<k<<endl;
	     
             for(Int_t l = 0; l<partEvent2->GetNumberOfParticles();l++)   //  ... on all particles
              {
               
                // if (N>MAXCOMB) break;
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
  UInt_t ii;
  for(ii = 0;ii<fNParticleFunctions;ii++)
    fParticleFunctions[ii]->Write();
                 
  for(ii = 0;ii<fNTrackFunctions;ii++)
    fTrackFunctions[ii]->Write();
                 
  for(ii = 0;ii<fNParticleAndTrackFunctions;ii++)
   fParticleAndTrackFunctions[ii]->Write();
}
/*************************************************************************************/

void AliHBTAnalysis::SetGlobalPairCut(AliHBTPairCut* cut)
{
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
  if (f == 0x0) return;
  
  if (fNParticleFunctions == fgkFctnArraySize)
   {
    Error("AliHBTAnalysis::AddParticleFunction","Can not add this function, not enough place in the array.");
   }
  fParticleFunctions[fNParticleFunctions] = f;
  fNParticleFunctions++;
  
  
}
void AliHBTAnalysis::AddParticleAndTrackFunction(AliHBTTwoPairFctn* f)
{
  if (f == 0x0) return;
  if (fNParticleAndTrackFunctions == fgkFctnArraySize)
   {
    Error("AliHBTAnalysis::AddParticleAndTrackFunction","Can not add this function, not enough place in the array.");
   }  
  fParticleAndTrackFunctions[fNParticleAndTrackFunctions] = f;
  fNParticleAndTrackFunctions++;
}


/*************************************************************************************/ 


/*************************************************************************************/  

Bool_t AliHBTAnalysis::RunCoherencyCheck()
{
 //Checks if both HBTRuns are similar
 //return true if error found
 //if they seem to be OK return false
 Int_t i;  
 cout<<"Checking HBT Runs Coherency"<<endl;
 
 cout<<"Number of events ..."; fflush(stdout);
 
 if (fReader->GetNumberOfPartEvents() == fReader->GetNumberOfTrackEvents() ) //check whether there is the same  number of events
  {
    cout<<"OK. "<<fReader->GetNumberOfTrackEvents()<<"  found"<<endl;
  }
 else
  { //if not the same - ERROR
   Error("AliHBTAnalysis::RunCoherencyCheck()",
         "Number of simulated events (%d) is not equal to number of reconstructed events(%d)",
         fReader->GetNumberOfPartEvents(),fReader->GetNumberOfTrackEvents());
   return kTRUE;
  }
 
 cout<<"Checking number of Particles AND Particles Types in each event ...";fflush(stdout);
 
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
 cout<<"  Done"<<endl;
 cout<<"  Everything looks OK"<<endl;
 return kFALSE;
}


/*************************************************************************************/  
 

/*************************************************************************************/  


