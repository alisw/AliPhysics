/*
 *  MergeQAMerged.C
 *  
 *
 *  Created by schutz on 29/09/09.
 *  Copyright 2009 CERN. All rights reserved.
 *
 */
MergeQAMerged(Char_t * outfile, Char_t * list)
{
  TFileMerger merger ; 
  merger.OutputFile(outfile); 
  TString QAfile ; 
  FILE * QAfiles = fopen(list, "r") ; 
  while ( QAfile.Gets(QAfiles) ){
    merger.AddFile(QAfile) ;     
  } 
  merger.Merge() ;   
}

