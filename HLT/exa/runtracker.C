// $Id$

/**
   Run this macro for cluster finder and track follower 
   (see steering class AliLevel3).
   In argument path, you have to provide the path to the directory 
   where the data files should be located. In case of reading from a rootfile, you have to
   make a symbolic link "digitfile.root", which points to the rootfile containing AliROOT 
   digits tree and a symbolic link "alirunfile.root" pointing to a file containing
   the ALIROOT geometry (TPC param). For NEWIO, make sure that the 
   file TPC.Digits.root is in the path!

   RUN with ALIROOT (not ROOT) if using root files.
*/

void runtracker(Int_t minslice=0,Int_t maxslice=35,Char_t* path="./",Int_t nevent=1,Char_t *opath="./")
{
  //Set your configuration here:
  AliLevel3::EFileType filetype=AliLevel3::kRoot; //Input is RLE binary files or rootfile.
  Bool_t pileup=kFALSE; //Assume input is pileup event = non RLE binary files.
  Int_t npatches = 1;   //Options; 1, 2 and 6.
  Char_t trackparams[] = "SetTrackingParameters_4000bf04.C"; //Set this to correspond 
                                                             //with mult. and BField
  
  //for aliroot the path should point to a file 
  //containing the tpc geometry called alirunfile.root
  //Bool_t isinit=AliL3Transform::Init(path,!binary);
  Bool_t isinit=AliL3Transform::Init(path,(filetype!=AliLevel3::kBinary));
  if(!isinit){
    cerr << "Could not create transform settings, please check log for error messages!" << endl;
    return;
  }

  for(Int_t ev=0; ev<nevent; ev++)
    {
      if(filetype==AliLevel3::kBinary)
	a = new AliLevel3();
      else 
	{
	  Char_t fname[1024];
	  if(filetype==AliLevel3::kRaw)
           sprintf(fname,"%s/raw.root",path);
	  else
           sprintf(fname,"%s/digitfile.root",path);
          a = new AliLevel3(fname);
	}
      
      a->Init(path,filetype,npatches);
      
      if(pileup)
	a->DoPileup();
      
      gROOT->LoadMacro(trackparams);
      SetTrackingParameters(a);

#if 0 /* little comments */
      phi_segments = 50;   //devide the space into phi_segments and eta_segments
      eta_segments = 100;  //to access the search of points to that area!
      trackletlength = 3;  //number of hits a tracklet has to have
      tracklength = 5;     //number of hits a track has to have 
      rowscopetracklet = 2;//search range of rows for a tracklet
      rowscopetrack = 2;   //search range of rows for a track
      min_pt_fit = 0;      
      maxangle = 1.31;     //maximum angle for the three point look ahead
      goodDist = 5;        //threshold distance between two hits when building tracklets
      maxphi=100;          //maximum phi difference for neighboring hits
      maxeta=100;          //maximum eta difference for neighboring hits
      hitChi2Cut = 100;    //maximum chi2 of added hit to track
      goodHitChi2 = 20;    //stop looking for next hit to add if chi2 is less than goodHitChi2
      trackChi2Cut = 50;   //maximum chi2 for track after final fit
#endif


      //a->DoRoi();    /*do region of interest*/
      //a->DoMc();     /*do monte carlo identification*/
      a->WriteFiles(opath); /*enable output*/
      
      a->ProcessEvent(minslice,maxslice);
      //a->DoBench("benchmark_0");

      delete a;
    }
}


