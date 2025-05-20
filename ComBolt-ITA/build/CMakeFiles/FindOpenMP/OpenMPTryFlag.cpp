// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.


#include <omp.h>
int main(void) {
#ifdef _OPENMP
  omp_get_max_threads();
  return 0;
#elif defined(__HIP_DEVICE_COMPILE__)
  return 0;
#else
  breaks_on_purpose
#endif
}
