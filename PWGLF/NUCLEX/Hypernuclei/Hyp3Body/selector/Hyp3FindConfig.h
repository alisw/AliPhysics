#ifndef HYP3FINDCONFIG_H
#define HYP3FINDCONFIG_H

const int kNVar = 3;
const int kNCuts[4]={5,5,5,5};
const float kNsigmaTPC[5][3] = {{3.5,3.5,3.5},{2.5,2.5,2.5},{4.5,4.5,4.5},{5.5,5.5,5.5},{6.5,6.5,6.5}};
const float kNsigmaTOF[5][3] = {{3.5,3.5,3.5},{2.5,2.5,2.5},{4.5,4.5,4.5},{5.5,5.5,5.5},{6.5,6.5,6.5}};

const int kNclusTPC[5][3] = {{100,100,100},{50,50,50},{75,75,75},{125,125,125},{150,150,150}};
const int kNclusITS[5][3] = {{0,0,0},{1,1,1},{2,2,2},{3,3,3},{4,4,4}};


#endif
