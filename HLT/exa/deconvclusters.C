//$Id$

/* Example of how to use the AliL3ClusterFitter to fit clusters
   to the track candidates given as a AliL3TrackArray.

   The path "path" should contain the link to the digitsfile, 
   and the directories called fitter (for the results) and hough
   (for the track files). For NEWIO, make sure that the file 
   TPC.Digits.root is in the path!

   Also provide the neccessary parameters in SetFitParameters.C.

   RUN with ALIROOT (not ROOT), no other method is 
   supported right now.
*/

void deconvclusters(char *path,int minslice=0,int maxslice=35,int nevent=1)
{
  
  AliL3Transform::Init(path,kTRUE);
  
  char filename[1024];
  AliL3FileHandler *file = new AliL3FileHandler();
  UInt_t ndigits=0;
  
  sprintf(filename,"%s/digitfile.root",path);
  
  file->SetAliInput(filename);
  int index=0;
  
  AliL3ClusterFitter *fitter = new AliL3ClusterFitter(path);
  
  char macroname[1024];
  gROOT->LoadMacro("SetFitParameters.C");
  SetFitParameters(fitter);
  
  int patch=-1;
  int rowrange[2] = {0,AliL3Transform::GetNRows()-1};
  
  for(int ev=0; ev<nevent; ev++)
    {
      fitter->LoadSeeds(rowrange,kFALSE,ev);//Takes input from global hough tracks
      
      for(int slice=minslice; slice<=maxslice; slice++)
	{
	  file->Init(slice,-1);
	  cout<<"Processing event "<<ev<<" slice "<<slice<<" patch "<<patch<<endl;
	  AliL3DigitRowData *digits = (AliL3DigitRowData*)file->AliAltroDigits2Memory(ndigits,ev);
	  
	  fitter->Init(slice,patch);
	  fitter->SetInputData(digits);
	  
	  fitter->FindClusters();
	  
	  fitter->WriteClusters();
	  
	  file->Free();
	  
	}
      
      //If you want a refit of the clusters;-------------------------
      AliL3Vertex vertex;
      AliL3TrackArray *tracks = fitter->GetSeeds(); //The seeds are the input tracks from circle HT
      AliL3Fitter *ft = new AliL3Fitter(&vertex,1);
      sprintf(filename,"%s/fitter/",path);
      ft->LoadClusters(filename,0,kTRUE);
      for(Int_t i=0; i<tracks->GetNTracks(); i++)
	{
	  track = tracks->GetCheckedTrack(i);
	  if(!track) continue;
	  if(track->GetNHits() < 40) continue;
	  ft->SortTrackClusters(track);
	  ft->FitHelix(track);
	  ft->UpdateTrack(track);
	}
      delete ft;
      //-------------------------------------------------------------
      
      fitter->WriteTracks(5); //Write the final tracks
      file->FreeDigitsTree();
    }
  delete fitter;
}


void deconvlocally(char *path,int minslice=0,int maxslice=17)
{
  
  AliL3Transform::Init(path,kTRUE);
  
  char filename[1024];
  AliL3FileHandler *file = new AliL3FileHandler();
  UInt_t ndigits=0;
  
  sprintf(filename,"%s/digitfile.root",path);
  
  file->SetAliInput(filename);
  int index=0;
  
  AliL3ClusterFitter *fitter = new AliL3ClusterFitter(path);
  
  char macroname[1024];
  sprintf(macroname,"%s/SetFitParameters.C",path);
  gROOT->LoadMacro(macroname);
  SetFitParameters(fitter);
  
  int patch=-1;
  
  for(int slice=minslice; slice<=maxslice; slice++)
    {
      file->Init(slice,patch);
      cout<<"Processing slice "<<slice<<" patch "<<patch<<endl;
      AliL3DigitRowData *digits = (AliL3DigitRowData*)file->AliAltroDigits2Memory(ndigits);
      
      fitter->Init(slice,patch);
      fitter->LoadLocalSegments();
      
      fitter->SetInputData(digits);
      
      fitter->FindClusters();
      fitter->WriteClusters();
      
      file->Free();
    
    }
  
  fitter->WriteTracks(5);
  delete fitter;
}

