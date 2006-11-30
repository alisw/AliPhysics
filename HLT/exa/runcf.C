// $Id$

/* Example macro showing the usage of the 
   cluster finder storing the space points
   in an array. 
*/

void runcf(Char_t *path)
{
  AliHLTTransform::Init(path,kTRUE);
  
  Char_t fname[1024];
  Char_t digitfile[1024];
  sprintf(digitfile,"%s/digitfile.root",path);
  
  
  AliHLTMemHandler *memory = new AliHLTMemHandler();
  AliHLTMemHandler *out = new AliHLTMemHandler();

  for(Int_t event=0; event<1; event++)
    {
      AliHLTFileHandler *file = new AliHLTFileHandler();
      file->SetAliInput(digitfile);
      for(Int_t slice=0; slice<=35; slice++)
	{
	  for(Int_t patch=0; patch<6; patch++)
	    {
	      cout<<"Processing event "<<event<<" slice "<<slice<<" patch "<<patch<<endl;
	      file->Init(slice,patch);
	      UInt_t ndigits=0;
	      UInt_t maxclusters=100000;
	      UInt_t pointsize = maxclusters*sizeof(AliHLTSpacePointData);
	      
	      AliHLTSpacePointData *points = (AliHLTSpacePointData*)memory->Allocate(pointsize);
	      AliHLTDigitRowData *digits = (AliHLTDigitRowData*)file->AliAltroDigits2Memory(ndigits,event);
	      AliHLTClustFinderNew *cf = new AliHLTClustFinderNew();
	      //cf->SetMatchWidth(2);
	      cf->InitSlice(slice,patch,AliHLTTransform::GetFirstRow(patch),AliHLTTransform::GetLastRow(patch),maxclusters);
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
