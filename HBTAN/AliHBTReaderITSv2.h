#ifndef ALIHBTREADERITSV2_H
#define ALIHBTREADERITSV2_H

#include "AliHBTReader.h"

#include <TString.h>


class AliHBTReaderITSv2: public AliHBTReaderTPC
{
  public:    
    AliHBTReaderITSv2(const Char_t* trackfilename = "AliITStracksV2.root",
                    const Char_t* clusterfilename = "AliITSclustersV2.root",
	const Char_t* goodtracksfilename = "good_tracks_its",
	const Char_t* galicefilename = "");
	
    virtual ~AliHBTReaderITSv2();
    
    Int_t Read(AliHBTRun* particles, AliHBTRun *tracks);//reads tracks and particles and puts them in runs
    
  protected:    
  private:
  public:
    ClassDef(AliHBTReaderITSv2,1)
};
