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

/* $Id$ */

//____________________________________________________________________
//////////////////////////////////////////////////////////////////////
//                                                                  //
// class AliReaderESD                                               //
//                                                                  //
// reader for ALICE Event Summary Data (ESD).                       //
//                                                                  //
// Piotr.Skowronski@cern.ch                                         //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TGliteXmlEventlist.h>
#include <TH1.h>
#include <TKey.h>
#include <TObjString.h>
#include <TPDGCode.h>
#include <TParticle.h>
#include <TString.h>
#include <TTree.h>

#include "AliAOD.h"
#include "AliAODParticle.h"
#include "AliAODParticleCut.h"
#include "AliAODRun.h"
#include "AliAnalysis.h"
#include "AliClusterMap.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliKalmanTrack.h"
#include "AliLog.h"
#include "AliReaderESD.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliStack.h"

ClassImp(AliReaderESD)

AliReaderESD::AliReaderESD(const Char_t* esdfilename, const Char_t* galfilename):
 fESDFileName(esdfilename),
 fGAlFileName(galfilename),
 fFile(0x0),
 fRunLoader(0x0),
 fKeyIterator(0x0),
 fReadSim(kFALSE),
 fCheckParticlePID(kFALSE),
 fReadMostProbableOnly(kFALSE),
 fNTrackPoints(0),
 fdR(0.0),
 fClusterMap(kFALSE),
 fITSTrackPoints(kFALSE),
 fITSTrackPointsType(AliTrackPoints::kITS),
 fMustTPC(kFALSE),
 fReadCentralBarrel(kTRUE),
 fReadMuon(kFALSE),
 fReadPHOS(kFALSE),
 fNTPCClustMin(0),
 fNTPCClustMax(1500),
 fTPCChi2PerClustMin(0.0),
 fTPCChi2PerClustMax(10e5),
 fChi2Min(0.0),
 fChi2Max(10e5),
 fC00Min(0.0),
 fC00Max(10e5),
 fC11Min(0.0),
 fC11Max(10e5),
 fC22Min(0.0),
 fC22Max(10e5),
 fC33Min(0.0),
 fC33Max(10e5),
 fC44Min(0.0),
 fC44Max(10e5),
 fTPCC00Min(0.0),
 fTPCC00Max(10e5),
 fTPCC11Min(0.0),
 fTPCC11Max(10e5),
 fTPCC22Min(0.0),
 fTPCC22Max(10e5),
 fTPCC33Min(0.0),
 fTPCC33Max(10e5),
 fTPCC44Min(0.0),
 fTPCC44Max(10e5)

{
  //cosntructor
  if ( ((Int_t)kNSpecies) != ((Int_t)AliPID::kSPECIES))
    Fatal("AliReaderESD","ESD defintions probobly changed. Ask Youra.");
}
/********************************************************************/
  
AliReaderESD::AliReaderESD(TObjArray* dirs,const Char_t* esdfilename, const Char_t* galfilename):
 AliReader(dirs), 
 fESDFileName(esdfilename),
 fGAlFileName(galfilename),
 fFile(0x0),
 fRunLoader(0x0),
 fKeyIterator(0x0),
 fReadSim(kFALSE),
 fCheckParticlePID(kFALSE),
 fReadMostProbableOnly(kFALSE),
 fNTrackPoints(0),
 fdR(0.0),
 fClusterMap(kFALSE),
 fITSTrackPoints(kFALSE),
 fITSTrackPointsType(AliTrackPoints::kITS),
 fMustTPC(kFALSE),
 fReadCentralBarrel(kTRUE),
 fReadMuon(kFALSE),
 fReadPHOS(kFALSE),
 fNTPCClustMin(0),
 fNTPCClustMax(150),
 fTPCChi2PerClustMin(0.0),
 fTPCChi2PerClustMax(10e5),
 fChi2Min(0.0),
 fChi2Max(10e5),
 fC00Min(0.0),
 fC00Max(10e5),
 fC11Min(0.0),
 fC11Max(10e5),
 fC22Min(0.0),
 fC22Max(10e5),
 fC33Min(0.0),
 fC33Max(10e5),
 fC44Min(0.0),
 fC44Max(10e5),
 fTPCC00Min(0.0),
 fTPCC00Max(10e5),
 fTPCC11Min(0.0),
 fTPCC11Max(10e5),
 fTPCC22Min(0.0),
 fTPCC22Max(10e5),
 fTPCC33Min(0.0),
 fTPCC33Max(10e5),
 fTPCC44Min(0.0),
 fTPCC44Max(10e5)
{
  //cosntructor
  if ( ((Int_t)kNSpecies) != ((Int_t)AliPID::kSPECIES))
    Fatal("AliReaderESD","ESD defintions probobly changed. Ask Youra.");
}
/********************************************************************/

AliReaderESD::~AliReaderESD()
{
 //desctructor
    delete fRunLoader;
    delete fKeyIterator;
    delete fFile;
}

/**********************************************************/
Int_t AliReaderESD::ReadNext()
{
//reads next event from fFile
  //fRunLoader is for reading Kine
  
  AliDebug(1,"Entered");
    
  if (fEventSim == 0x0)  fEventSim = new AliAOD();
  if (fEventRec == 0x0)  fEventRec = new AliAOD();
  
  fEventSim->Reset();
  fEventRec->Reset();
        
  do  //do{}while; is OK even if 0 dirs specified. In that case we try to read from "./"
    {
      if (fFile == 0x0)
       {
         fFile = OpenFile(fCurrentDir);//rl is opened here
         if (fFile == 0x0)
           {
             Error("ReadNext","Cannot get fFile for dir no. %d",fCurrentDir);
             fCurrentDir++;
             continue;
           }
         fCurrentEvent = 0;
       }
     TString esdname = "ESD";
     esdname+=fCurrentEvent;
     AliESD* esd = dynamic_cast<AliESD*>(fFile->Get(esdname));
     if (esd == 0x0)
      {
        AliDebug(3,Form("Can not find AliESD object named %s",esdname.Data()));
        fCurrentDir++;
        delete fFile;//we have to assume there is no more ESD objects in the fFile
        fFile = 0x0;
        delete fRunLoader;
        fRunLoader = 0x0;
        continue;
      }
      ReadESD(esd);
      
      fCurrentEvent++;
      fNEventsRead++;
      delete esd;
      return 0;//success -> read one event
    }while(fCurrentDir < GetNumberOfDirs());//end of loop over directories specified in fDirs Obj Array  
   
  return 1; //no more directories to read
}

/**********************************************************/
Int_t AliReaderESD::ReadESD(AliESD* esd)
{
//Reads esd data
  if (esd == 0x0)
   {
     Error("ReadESD","ESD is NULL");
     return 1;
   }
  
  // seperate each method
  if (fReadCentralBarrel) ReadESDCentral(esd);

  if (fReadMuon) ReadESDMuon(esd);

  if (fReadPHOS) ReadESDPHOS(esd);

  return 1;
}

/**********************************************************/
Int_t AliReaderESD::ReadESDCentral(AliESD* esd)
{
  //****** Tentative particle type "concentrations"
  static const Double_t concentr[5]={0.05, 0., 0.85, 0.10, 0.05};
  
  Double_t pidtable[kNSpecies];//array used for reading pid probabilities from ESD track
  Double_t w[kNSpecies];
  Double_t mom[3];//momentum
  Double_t pos[3];//position
  Double_t vertexpos[3];//vertex position
  //Reads one ESD

  TDatabasePDG* pdgdb = TDatabasePDG::Instance();
  if (pdgdb == 0x0)
   {
     Error("ReadESD","Can not get PDG Database Instance.");
     return 1;
   }
   
  Float_t mf = esd->GetMagneticField()/10.; //AliESD::GetMagnField returns mf in kG

  if ( (mf == 0.0) && ((fNTrackPoints > 0) || fITSTrackPoints) )
   {
      Error("ReadESD","Magnetic Field is 0 and Track Points Demended. Skipping to next event.");

   }

  if (fITSTrackPoints)
   {
     Info("ReadESD","Magnetic Field is %f",mf);
     AliKalmanTrack::SetMagneticField(mf);
   }
 
  AliStack* stack = 0x0;
  if (fReadSim && fRunLoader)
   {
     fRunLoader->GetEvent(fCurrentEvent);
     stack = fRunLoader->Stack();
   }

  const AliESDVertex* vertex = esd->GetVertex();
  if (vertex == 0x0)
   {
     Info("ReadESD","ESD returned NULL pointer to vertex - assuming (0.0,0.0,0.0)");
     vertexpos[0] = 0.0;
     vertexpos[1] = 0.0;
     vertexpos[2] = 0.0;
   }
  else
   {
     vertex->GetXYZ(vertexpos);
   }
   
  AliDebug(1,Form("Primary Vertex is (%f,%f,%f)",vertexpos[0],vertexpos[1],vertexpos[2]));
   
  Info("ReadESD","Reading Event %d",fCurrentEvent);

  Int_t ntr = esd->GetNumberOfTracks();
  Info("ReadESD","Found %d tracks.",ntr);
  for (Int_t i = 0;i<ntr; i++)
   {
     AliESDtrack *esdtrack = esd->GetTrack(i);
     if (esdtrack == 0x0)
      {
        Error("Next","Can not get track %d", i);
        continue;
      }

     //if (esdtrack->HasVertexParameters() == kFALSE) 
     if ((esdtrack->GetStatus() & AliESDtrack::kITSrefit) == kFALSE)
      {
        AliDebug(3,Form("Particle skipped: Data at vertex not available."));
        continue;
      }

     if (fMustTPC)
      {
       if ((esdtrack->GetStatus() & AliESDtrack::kTPCin) == kFALSE)
        {
          AliDebug(3,"Particle skipped: Was not reconstructed in TPC.");
          continue;
        }
      }     

     if ((esdtrack->GetStatus() & AliESDtrack::kESDpid) == kFALSE) 
      {
        AliDebug(3,"Particle skipped: PID BIT is not set.");
        continue;
      }


     Double_t extx;
     Double_t extp[5];
     esdtrack->GetConstrainedExternalParameters(extx,extp);
     if (extp[4] == 0.0)
      {
        AliDebug(3,"Track has 0 contrianed curvature -> Probobly parameters never updated. Skipping.");
        continue;
      } 
     esdtrack->GetESDpid(pidtable);
     esdtrack->GetConstrainedPxPyPz(mom);
     esdtrack->GetConstrainedXYZ(pos);
     pos[0] -= vertexpos[0];//we are interested only in relative position to Primary vertex at this point
     pos[1] -= vertexpos[1];
     pos[2] -= vertexpos[2];

     Int_t charge = (extp[4] > 0)?1:-1;//if curvature=charg/Pt is positive charge is positive

     //Particle from kinematics
     AliAODParticle* particle = 0;
     Bool_t keeppart = kFALSE;
     if ( fReadSim && stack  )
      {
        if (esdtrack->GetLabel() < 0) continue;//this is fake -  we are not able to match any track 
        TParticle *p = stack->Particle(esdtrack->GetLabel());
        if (p==0x0) 
         {
           Error("ReadNext","Can not find track with such label.");
           continue;
         }
         
        if (fCheckParticlePID)
         {
           if(Rejected(p->GetPdgCode())) 
            {
              AliDebug(6,Form("Simulated Particle PID (%d) did not pass the cut.",p->GetPdgCode()));
              continue; //check if we are intersted with particles of this type 
            }
         }  
//           if(p->GetPdgCode()<0) charge = -1;
        particle = new AliAODParticle(*p,i);

      }
      
     if(CheckTrack(esdtrack)) continue;
      
     //Here we apply Bayes' formula
     Double_t rc=0.;
     for (Int_t s=0; s<AliPID::kSPECIES; s++) rc+=concentr[s]*pidtable[s];
     if (rc==0.0) 
      {
        AliDebug(3,"Particle rejected since total bayessian PID probab. is zero.");
        continue;
      }

     for (Int_t s=0; s<AliPID::kSPECIES; s++) w[s]=concentr[s]*pidtable[s]/rc;

     if (AliDebugLevel() > 4)
      { 
        Info("ReadNext","###########################################################################");
        Info("ReadNext","Momentum: %f %f %f",mom[0],mom[1],mom[2]);
        Info("ReadNext","Position: %f %f %f",pos[0],pos[1],pos[2]);
        TString msg("Pid list got from track:");
        for (Int_t s = 0;s<kNSpecies;s++)
         {
           msg+="\n    ";
           msg+=s;
           msg+="(";
           msg+=charge*GetSpeciesPdgCode((ESpecies)s);
           msg+="): ";
           msg+=w[s];
           msg+=" (";
           msg+=pidtable[s];
           msg+=")";
         }
        Info("ReadNext","%s",msg.Data());
      }//if (AliDebugLevel()>4)

      AliTrackPoints* tpts = 0x0;
      if (fNTrackPoints > 0) 
       {
         tpts = new AliTrackPoints(fNTrackPoints,esdtrack,mf*10.0,fdR);
//         tpts->Move(-vertexpos[0],-vertexpos[1],-vertexpos[2]);
       }

      AliTrackPoints* itstpts = 0x0;
      if (fITSTrackPoints) 
       {
         itstpts = new AliTrackPoints(fITSTrackPointsType,esdtrack,mf*10.0);
//         itstpts->Move(-vertexpos[0],-vertexpos[1],-vertexpos[2]);
       }

      AliClusterMap* cmap = 0x0; 
      if (  fClusterMap ) 
       {
         cmap = new AliClusterMap(esdtrack);
       }
     
     //If this flag fReadMostProbableOnly is false the 
     //loop over species (see "LOOPPIDS") is over all possible PIDs
     //in other case the most probablle PID is searched
     //and the loop is limited to that PID only

     Int_t firstspecie = 0;
     Int_t lastspecie = kNSpecies;
     
     if (fReadMostProbableOnly)
      { 
        //find the most probable PID
        Int_t spec = 0;
        Float_t maxprob = w[0];
        for (Int_t s=1; s<AliPID::kSPECIES; s++) 
         {
           if (w[s]>maxprob)
            {
              maxprob = w[s];
              spec = s;
            }
         }
        firstspecie = spec;
        lastspecie = spec + 1;
      }
     
     for (Int_t s = firstspecie; s<lastspecie; s++)//LOOPPIDS
      {
        Int_t pdgcode = charge*GetSpeciesPdgCode((ESpecies)s);
        Float_t pp = w[s];
        if (pp == 0.0) 
         {
           AliDebug(6,Form("Probability of being PID %d is zero. Continuing.",pdgcode));
           continue;
         }

        if(Rejected(pdgcode)) 
         {
           AliDebug(6,Form("PID (%d) did not pass the cut.",pdgcode));
           continue; //check if we are intersted with particles of this type 
         }

        Double_t mass = pdgdb->GetParticle(pdgcode)->Mass();
        Double_t tEtot = TMath::Sqrt( mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2] + mass*mass);//total energy of the track

        AliAODParticle* track = new AliAODParticle(pdgcode, w[s],i, 
                                                   mom[0], mom[1], mom[2], tEtot,
                                                   pos[0], pos[1], pos[2], 0.);
        //copy probabilitis of other species (if not zero)
        for (Int_t k = 0; k<kNSpecies; k++)
         {
           if (k == s) continue;
           if (w[k] == 0.0) continue;
           track->SetPIDprobability(charge*GetSpeciesPdgCode( (ESpecies)k ),w[k]);
         }

        if(Rejected(track))//check if meets all criteria of any of our cuts
                       //if it does not delete it and take next good track
         { 
           AliDebug(5,"Track did not pass the cut");
           delete track;
           continue;
         }

        //Single Particle cuts on cluster map and track points rather do not have sense
        if (tpts)
         {
           track->SetTPCTrackPoints(tpts);
         }
         
        if (itstpts)
         {
           track->SetITSTrackPoints(itstpts); 
         }

        if (cmap) 
         { 
           track->SetClusterMap(cmap);
         }

        fEventRec->AddParticle(track);
        if (particle) fEventSim->AddParticle(particle);
        keeppart = kTRUE;

        if (AliDebugLevel() > 4 )
         {
           Info("ReadNext","\n\nAdding Particle with incarnation %d",pdgcode);
           track->Print();
           if (particle) particle->Print();
           Info("ReadNext","\n----------------------------------------------\n");
         }
      }//for (Int_t s = 0; s<kNSpecies; s++)

     if (keeppart == kFALSE) 
      {
        delete particle;//particle was not stored in event
        delete tpts;
        delete itstpts;
        delete cmap;
      }
     else
      {
        if ( fReadSim && stack  )
         {
           if (particle->P() < 0.00001)
            {
              Info("ReadNext","###################################");
              Info("ReadNext","###################################");
              Info("ReadNext","Track Label %d",esdtrack->GetLabel());
              TParticle *p = stack->Particle(esdtrack->GetLabel());
              Info("ReadNext","");
              p->Print();
              Info("ReadNext","");
              particle->Print();
            }
         }   
      } 

   }//for (Int_t i = 0;i<ntr; i++)  -- loop over tracks

  Info("ReadNext","Read %d tracks and %d particles from event %d (event %d in dir %d).",
         fEventRec->GetNumberOfParticles(), fEventSim->GetNumberOfParticles(),
         fNEventsRead,fCurrentEvent,fCurrentDir);
  fTrackCounter->Fill(fEventRec->GetNumberOfParticles());
  
  /******************************************************/
  /******     Setting glevet properties     *************/
  /******************************************************/
    
  if (fEventRec->GetNumberOfParticles() > 0)
   {
     fEventRec->SetPrimaryVertex(vertexpos[0],vertexpos[1],vertexpos[2]);
   }
  return 0;
}

/**********************************************************/
Int_t AliReaderESD::ReadESDMuon(AliESD* esd)
{

  Double_t vertexpos[3];//vertex position, assuming no secondary decay

  const AliESDVertex* vertex = esd->GetVertex();

  if (vertex == 0x0) {
    Info("ReadESD","ESD returned NULL pointer to vertex - assuming (0.0,0.0,0.0)");
    vertexpos[0] = 0.0;
    vertexpos[1] = 0.0;
    vertexpos[2] = 0.0;
  } else {
    vertex->GetXYZ(vertexpos);
  }

 Int_t nTracks = (Int_t)esd->GetNumberOfMuonTracks() ;

 AliDebug(1,Form("Reading Event %d \nFound %d tracks.",fCurrentEvent,nTracks));

 // settings
  Float_t Chi2Cut = 100.;
  Float_t PtCutMin = 1.;
  Float_t PtCutMax = 10000.;
  Float_t muonMass = 0.105658389;
  Int_t pdgcode = -13;
  Double_t thetaX, thetaY, pYZ;
  Double_t pxRec1, pyRec1, pzRec1, E1;
  Int_t charge;

  Int_t ntrackhits;
  Double_t fitfmin;

  TLorentzVector fV1;
  fEventRec->Reset();
  for (Int_t iTrack = 0; iTrack <  nTracks;  iTrack++) {

      AliESDMuonTrack* muonTrack = esd->GetMuonTrack(iTrack);

      thetaX = muonTrack->GetThetaX();
      thetaY = muonTrack->GetThetaY();

      pYZ     =  1./TMath::Abs(muonTrack->GetInverseBendingMomentum());
      pzRec1  = - pYZ / TMath::Sqrt(1.0 + TMath::Tan(thetaY)*TMath::Tan(thetaX));
      pxRec1  = pzRec1 * TMath::Tan(thetaX);
      pyRec1  = pzRec1 * TMath::Tan(thetaY);
      charge = Int_t(TMath::Sign(1.,muonTrack->GetInverseBendingMomentum()));
      E1 = TMath::Sqrt(muonMass * muonMass + pxRec1 * pxRec1 + pyRec1 * pyRec1 + pzRec1 * pzRec1);
      fV1.SetPxPyPzE(pxRec1, pyRec1, pzRec1, E1);

      ntrackhits = muonTrack->GetNHit();
      fitfmin    = muonTrack->GetChi2();

      // transverse momentum
      Float_t pt1 = fV1.Pt();

      // chi2 per d.o.f.
      Float_t ch1 =  fitfmin / (2.0 * ntrackhits - 5);

      if ((ch1 < Chi2Cut) && (pt1 > PtCutMin) && (pt1 < PtCutMax)) {
	AliAODParticle* track = new AliAODParticle(pdgcode*charge,1,iTrack, 
                                                   pxRec1, pyRec1,pzRec1, E1,
                                                   vertexpos[0], vertexpos[1], vertexpos[2], 0.);
        fEventRec->AddParticle(track);
      }
 
  }
  fTrackCounter->Fill(fEventRec->GetNumberOfParticles());
  return 0;
}

/**********************************************************/

void AliReaderESD::Rewind()
{
  //rewinds reading 
  //  delete fKeyIterator;
  delete fFile;
  fFile = 0x0;
  // fKeyIterator = 0x0;
  delete fRunLoader;
  fRunLoader = 0x0;
  fCurrentDir = 0;
  fNEventsRead = 0;
  if (fEventList) fEventList->Reset(); 
  if (fTrackCounter) fTrackCounter->Reset();
}
/**********************************************************/

TFile* AliReaderESD::OpenFile(Int_t n)
{
//opens fFile with  tree
 if (fEventList)
  { 
    if (fCurrentDir > n)
     {
       fEventList->Reset();
       fCurrentDir = 0;
     }
    
    while (fCurrentDir < n)
     {
       fEventList->Next();
       fCurrentDir++;
     }
    fEventList->Next();
  }
 
 const TString& dirname = GetDirName(n);
 if (dirname == "")
  {
   Error("OpenFiles","Can not get directory name");
   return 0x0;
  }
 TString filename;
 
 if (fEventList)
  {
    filename = fEventList->GetURL(fESDFileName);
  }
 else
  {
    filename = dirname +"/"+ fESDFileName;
  }
 Info("OpenFile","%s ==> %s",fESDFileName.Data(),filename.Data()); 
 TFile *ret = TFile::Open(filename.Data()); 

 if ( ret == 0x0)
  {
    Error("OpenFiles","Can't open fFile %s",filename.Data());
    return 0x0;
  }
 if (!ret->IsOpen())
  {
    Error("OpenFiles","Can't open fFile  %s",filename.Data());
    return 0x0;
  }
 
 if (fReadSim )
  {
    TString gafilename;
    if (fEventList)
     {
       gafilename = fEventList->GetURL(fGAlFileName);
     }
    else
     {
       gafilename = dirname +"/"+ fGAlFileName;
     }
   Info("OpenFile","%s ==> %s",fGAlFileName.Data(),gafilename.Data()); 

   fRunLoader = AliRunLoader::Open(gafilename);
   
   if (fRunLoader == 0x0)
    {
      Error("OpenFiles","Can't get RunLoader for directory %s",dirname.Data());
      delete ret;
      return 0x0;
    }
    
   fRunLoader->LoadHeader();
   
   if (fEventList)
    {
      TString kinefilename = fEventList->GetURL("Kinematics.root");
      fRunLoader->SetKineFileName(kinefilename);
      Info("OpenFile","%s ==> %s","Kinematics.root",kinefilename.Data()); 
    } 
   
   if (fRunLoader->LoadKinematics())
    {
      Error("Next","Error occured while loading kinematics.");
      delete fRunLoader;
      delete ret;
      return 0x0;
    }
  }
   
 return ret;
}
/**********************************************************/

Int_t AliReaderESD::GetSpeciesPdgCode(ESpecies spec)//skowron
{
  //returns pdg code from the PID index
  //ask jura about charge
  switch (spec)
   {
     case kESDElectron:
       return kPositron;
       break;
     case kESDMuon:
       return kMuonPlus;
       break;
     case kESDPion:
       return kPiPlus;
       break;
     case kESDKaon:
       return kKPlus;
       break;
     case kESDProton:
       return kProton;
       break;
     default:
       ::Warning("GetSpeciesPdgCode","Specie with number %d is not defined.",(Int_t)spec);
       break;
   }
  return 0;
}
/********************************************************************/
Bool_t AliReaderESD::CheckTrack(AliESDtrack* t) const
{
  //Performs check of the track
  
  if ( (t->GetConstrainedChi2() < fChi2Min) || (t->GetConstrainedChi2() > fChi2Max) ) return kTRUE;
  
  if ( (t->GetTPCclusters(0x0) < fNTPCClustMin) || (t->GetTPCclusters(0x0) > fNTPCClustMax) ) return kTRUE;

  if (t->GetTPCclusters(0x0) > 0)
   {
     Float_t chisqpercl = t->GetTPCchi2()/((Double_t)t->GetTPCclusters(0x0));
     if ( (chisqpercl < fTPCChi2PerClustMin) || (chisqpercl > fTPCChi2PerClustMax) ) return kTRUE;
   }

  Double_t cc[15];
  t->GetConstrainedExternalCovariance(cc);

  if ( (cc[0]  < fC00Min) || (cc[0]  > fC00Max) ) return kTRUE;
  if ( (cc[2]  < fC11Min) || (cc[2]  > fC11Max) ) return kTRUE;
  if ( (cc[5]  < fC22Min) || (cc[5]  > fC22Max) ) return kTRUE;
  if ( (cc[9]  < fC33Min) || (cc[9]  > fC33Max) ) return kTRUE;
  if ( (cc[14] < fC44Min) || (cc[14] > fC44Max) ) return kTRUE;


  t->GetInnerExternalCovariance(cc);

  if ( (cc[0]  < fTPCC00Min) || (cc[0]  > fTPCC00Max) ) return kTRUE;
  if ( (cc[2]  < fTPCC11Min) || (cc[2]  > fTPCC11Max) ) return kTRUE;
  if ( (cc[5]  < fTPCC22Min) || (cc[5]  > fTPCC22Max) ) return kTRUE;
  if ( (cc[9]  < fTPCC33Min) || (cc[9]  > fTPCC33Max) ) return kTRUE;
  if ( (cc[14] < fTPCC44Min) || (cc[14] > fTPCC44Max) ) return kTRUE;

  return kFALSE;

}
/********************************************************************/

void AliReaderESD::SetChi2Range(Float_t min, Float_t max)
{
  //sets range of Chi2 per Cluster
  fChi2Min = min;
  fChi2Max = max;
}
/********************************************************************/

void AliReaderESD::SetTPCNClustersRange(Int_t min,Int_t max)
{
 //sets range of Number Of Clusters that tracks have to have
 fNTPCClustMin = min;
 fNTPCClustMax = max;
}
/********************************************************************/

void AliReaderESD::SetTPCChi2PerCluserRange(Float_t min, Float_t max)
{
  //sets range of Chi2 per Cluster
  fTPCChi2PerClustMin = min;
  fTPCChi2PerClustMax = max;
}
/********************************************************************/

void AliReaderESD::SetC00Range(Float_t min, Float_t max)
{
 //Sets range of C00 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC00Min = min;
  fC00Max = max;
}
/********************************************************************/

void AliReaderESD::SetC11Range(Float_t min, Float_t max)
{
 //Sets range of C11 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC11Min = min;
  fC11Max = max;
}
/********************************************************************/

void AliReaderESD::SetC22Range(Float_t min, Float_t max)
{
 //Sets range of C22 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC22Min = min;
  fC22Max = max;
}
/********************************************************************/

void AliReaderESD::SetC33Range(Float_t min, Float_t max)
{
 //Sets range of C33 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC33Min = min;
  fC33Max = max;
}
/********************************************************************/

void AliReaderESD::SetC44Range(Float_t min, Float_t max)
{
 //Sets range of C44 parameter of covariance matrix of the track
 //it defines uncertainty of the momentum
  fC44Min = min;
  fC44Max = max;
}
