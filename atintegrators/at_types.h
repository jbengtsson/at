#ifndef ATTYPES_H
#define ATTYPES_H

#ifndef OMP_PARTICLE_THRESHOLD
#define OMP_PARTICLE_THRESHOLD (10)
#endif

#define PS_DIM 6

#define C0     2.99792458e8


struct parameters {
  int    nturn;
  double RingLength;
  double T0;
};

#endif /*ATTYPES_H*/
