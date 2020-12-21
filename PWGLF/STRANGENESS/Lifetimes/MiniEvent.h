#ifndef Lifetimes_MiniEvent_h
#define Lifetimes_MiniEvent_h

namespace Lifetimes {

struct MiniEvent {
  Double32_t fXvtx;         //[-20,20,16]
  Double32_t fYvtx;         //[-20,20,16]
  Double32_t fZvtx;         //[-20,20,16]
  Double32_t fMultiplicity; //[-1,100,16]
};

}
#endif