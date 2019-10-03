
Int_t localMergeFiles(Char_t * outfile=NULL, Char_t * list=NULL)
{
  Int_t filesCounter=0;
  TFileMerger merger ; 
  if (!list) {
    printf("Invalid list of files given as input: nothing done\n");
    return 0;
  }	
  merger.OutputFile(outfile); 
  TString infile ; 
  FILE * files = fopen(list, "r") ; 
  while ( infile.Gets(files) ){
    if (merger.AddFile(infile)) filesCounter++;     
  } 
  printf("Number of files to be merged = %i\n",filesCounter);
  merger.Merge();  
  return 1; 
}
