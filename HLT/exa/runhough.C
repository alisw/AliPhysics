//$Id$

/**
   Run this macro for Hough track candidate finder
   (see steering class AliL3Hough).
   In argument path, you have to provide the path to 
   the directory where the data files should be located. 
   In case of reading from a rootfile, you have to
   make a symbolic link "digitfile.root", which points 
   to the rootfile containing AliROOT digits tree 
   and a symbolic link "alirunfile.root" pointing to a file 
   containing the ALIROOT geometry (TPC param).
   For NEWIO, make sure that the file TPC.Digits.root is in the path!
 
   Also provide the neccessary parameters in SetHoughParameters.C.

   RUN with ALIROOT (not ROOT) if using root files.

*/

void runhough(Char_t *path,Char_t *outpath,int s1=0,int s2=35,int nevent=1)
{

  Bool_t isinit=AliL3Transform::Init(path,kTRUE);
  if(!isinit){
    cerr << "Could not create transform settings, please check log for error messages!" << endl;
    return;
  }
  
  hough = new AliL3Hough();

  char macroname[1024];
  sprintf(macroname,"SetHoughParameters.C",path);
  gROOT->LoadMacro(macroname);
  SetHoughParameters(hough,path);

  for(int ev=0; ev<nevent; ev++)
    {
      for(int slice=s1; slice<=s2; slice++)
	{
	  cout<<"Processing slice "<<slice<<endl;
	  hough->ReadData(slice,ev);
	  hough->Transform();
	  hough->AddAllHistograms();
	  hough->FindTrackCandidates();
	  hough->AddTracks();
	}
      hough->WriteTracks(outpath);
    }
  
  delete hough;
  
}

