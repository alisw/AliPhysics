// $Id$

/* Example macro showing the usage of the 
   cluster finder storing the space points
   in an array. 
*/

void runcf(Char_t *path)
{
  AliL3Transform::Init(path,kTRUE);
  
  Char_t fname[1024];
  Char_t digitfile[1024];
  sprintf(digitfile,"%s/digitfile.root",path);
  
  
  AliL3MemHandler *memory = new AliL3MemHandler();
  AliL3MemHandler *out = new AliL3MemHandler();

  for(Int_t event=0; event<1; event++)
    {
      AliL3FileHandler *file = new AliL3FileHandler();
      file->SetAliInput(digitfile);
      for(Int_t slice=0; slice<=35; slice++)
	{
	  for(Int_t patch=0; patch<6; patch++)
	    {
	      cout<<"Processing event "<<event<<" slice "<<slice<<" patch "<<patch<<endl;
	      file->Init(slice,patch);
	      UInt_t ndigits=0;
	      UInt_t maxclusters=100000;
	      UInt_t pointsize = maxclusters*sizeof(AliL3SpacePointData);
	      
	      AliL3SpacePointData *points = (AliL3SpacePointData*)memory->Allocate(pointsize);
	      AliL3DigitRowData *digits = (AliL3DigitRowData*)file->AliAltroDigits2Memory(ndigits,event);
	      AliL3ClustFinderNew *cf = new AliL3ClustFinderNew();
	      //cf->SetMatchWidth(2);
	      cf->InitSlice(slice,patch,AliL3Transform::GetFirstRow(patch),AliL3Transform::GetLastRow(patch),maxclusters);
	      cf->SetSTDOutput(kTRUE);
	      cf->SetThreshold(5);
	      cf->SetDeconv(kTRUE);
	      cf->SetCalcErr(kTRUE);
	      cf->SetOutputArray(points);
	      cf->Read(ndigits,digits);
	      cf->ProcessDigits();
	      Int_t npoints = cf->GetNumberOfClusters();
	      
	      sprintf(fname,"%s/points_%d_%d_%d.raw",path,event,slice,patch);
	      
	      //Transform points into global system
	      //out->Transform(npoints,points,slice);
	      
	      out->SetBinaryOutput(fname);
	      out->Memory2Binary(npoints,points);
	      out->CloseBinaryOutput();
	      out->Free();

	      memory->Free();
	      file->Free();
	      delete cf;
	    }
	  
	}
      delete file;
    }

  delete memory;
  delete out;
}
