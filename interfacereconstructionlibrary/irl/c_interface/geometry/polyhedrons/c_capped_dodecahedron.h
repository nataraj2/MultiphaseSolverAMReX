// This file is part of the Interface Reconstruction Library (IRL),
// a library for interface reconstruction and computational geometry operations.
//
// Copyright (C) 2019 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef IRL_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_CAPPED_DODECAHEDRON_H_
#define IRL_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_CAPPED_DODECAHEDRON_H_

#include "irl/geometry/polyhedrons/capped_dodecahedron.h"

extern "C" {

struct c_CapDod {
  IRL::CappedDodecahedron *obj_ptr = nullptr;
};

void c_CapDod_new(c_CapDod *a_self);
void c_CapDod_delete(c_CapDod *a_self);
void c_CapDod_construct(c_CapDod *a_self, const double *a_dodecahedron);
double c_CapDod_calculateVolume(const c_CapDod *a_self);
void c_CapDod_calculateCentroid(c_CapDod *a_ptr, double *a_centroid);
void c_CapDod_adjustCapToMatchVolume(c_CapDod *a_self,
                                     const double *a_correct_signed_volume);
void c_CapDod_getBoundingPts(const c_CapDod *a_self, double *a_lower_pt,
                             double *a_upper_pt);
void c_CapDod_getPt(const c_CapDod *a_self, const int *a_index, double *a_pt);
}

#endif // SRC_C_INTERFACE_GEOMETRY_POLYHEDRONS_C_CAPPED_DODECAHEDRON_H_
