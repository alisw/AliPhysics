// $Id$

void read(int min=0,int max=35)
{

  for(int slice=0; slice<35; slice++)
    {
      char fname[256];
      //sprintf(fname,"/prog/alice/data/Rawdata/PileUp/digits_%d_0.raw",slice);
      sprintf(fname,"/prog/alice/data/Rawdata/1_patch/pp/test_pileup/digits_%d_0.raw",slice);
      //sprintf(fname,"digits_%d_0.raw",slice);
      file = new AliL3FileHandler();
      if(!file->SetBinaryInput(fname))
	{
	  cerr<<"Error opening file "<<fname<<endl;
	  return;
	}
      int row[2]={0,175};
      file->Init(slice,0,row);
      UInt_t size;
      char name[256];
      AliL3DigitRowData *data = file->CompBinary2Memory(size);
      
      for(Int_t r=0; r<175; r++)
	{
	  UInt_t padrow=data->fRow;
	  AliL3DigitData *dPt = (AliL3DigitData*)data->fDigitData;
	  cout<<"padrow "<<padrow<<" ndigits "<<data->fNDigit<<endl;
	  
	  for(Int_t d=0; d<data->fNDigit; d++)
	    {
	      if(d>0 && dPt[d]->fPad == dPt[d-1]->fPad && dPt[d]->fTime == dPt[d-1]->fTime)
		cout<<"Slice "<<slice<<" padrow "<<padrow<<" pad "<<(int)dPt[d]->fPad<<" time "
		    <<(int)dPt[d]->fTime<<" charge "<<(int)dPt[d]->fCharge<<endl;
	    }
	  
	  file->UpdateRowPointer(data);
	}
      
      file->CloseBinaryInput();
      delete file;
    }
}
