#ifndef __TRACK_H__
#define __TRACK_H__

/* element interface */


struct lattice {
  elem_type *next;
  int       N;
};

template<typename T>
void track_element(T *x, const elem_type *Elem);

#endif

