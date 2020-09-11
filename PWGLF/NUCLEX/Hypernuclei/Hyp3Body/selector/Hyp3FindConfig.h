#ifndef HYP3FINDCONFIG_H
#define HYP3FINDCONFIG_H

const int kNBinCt = 10;
const int kNBinPt = 10;
//number of variables for the cuts
const int kNVar = 3;
//number of different cuts
const int kNCut = 5;
//kCuts[0][0][] must be the choosen cut
const float kCuts[kNVar][kNCut][3] = {{{3.5,3.5,3.5},{2.5,2.5,2.5},{4.5,4.5,4.5},{5.5,5.5,5.5},{6.5,6.5,6.5}},
{{100,100,100},{50,50,50},{75,75,75},{125,125,125},{150,150,150}},
{{0,0,0},{1,1,1},{2,2,2},{3,3,3},{4,4,4}}};


#endif
