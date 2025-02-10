#pragma once

#define RpSize 3
#define RpSizeGX 2
#define ParmSize 6
#define ParmSizeGX 10
#define fluxSize 3

#define D2(s1, i1, i2) ((i1)+(i2)*(s1))
#define D3(s1, s2, i1, i2, i3) ((i1)+((i2)+(i3)*(s2))*(s1))