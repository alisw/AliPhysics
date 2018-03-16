#include "AliMagFast.h"
#include "AliMagFastDip2k.h"
typedef unsigned short ushort;
const float Infinity = INFINITY;
AliMagFast::SegmentEnd_t dip2k_z_slices[] = {
};
AliMagFast::SegmentSearch_t dip2k_z_segments[] = {
};
AliMagFast::SegmentSearch_t dip2k_z = {-1, 0.0, -1760.0, dip2k_z_slices, dip2k_z_segments};
AliMagFast::ChebFormula_t dip2k_params[] = {
};
