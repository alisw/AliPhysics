Re-implementation of tracklet dNch/deta code 
============================================

The code in this sub-directory is a re-implementation of Ruben's code
for analysing SPD tracklet data for dNch/deta.  

The code differs from Ruben's mostly in organization e.g., the MC
part is split into a separate (derived) class, the output is organized
in to a hierarchy of containers, and the Terminate method tries to do
as many calculations as possible up front. 

Eventually, the code here can go into a binary compiled library. 

Content
-------

* `AliTrackletdNdetaPost.C` 
  Post processing class 
  
* `AliTrackletdNdetaTask.C`
  Analysis tasks 
  
* `AliTrackletdNdetaUtils.C`
  Common utilities used by tasks and post processing 

* `TrackletdNdetaTrain.C` 
  Train definition 

* `FixPaths.C`
  Fix header paths on AliEn
  
* `runTwo.sh` 
  Run pass on both data and simulation 
  
* `test.sh` 
  Run a test pass 
  
