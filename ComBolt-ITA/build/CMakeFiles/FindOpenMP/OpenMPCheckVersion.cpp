// SPDX-License-Identifier: MIT
// Copyright (c) 2025 Seyed Farid Taghavi, Technische Universität München – DFG Project (grant no. 517518417)
//
// This file is part of a DFG-funded research project (grant no. 517518417).
// See the top-level LICENSE file for license details.


#include <stdio.h>
#include <omp.h>
const char ompver_str[] = { 'I', 'N', 'F', 'O', ':', 'O', 'p', 'e', 'n', 'M',
                            'P', '-', 'd', 'a', 't', 'e', '[',
                            ('0' + ((_OPENMP/100000)%10)),
                            ('0' + ((_OPENMP/10000)%10)),
                            ('0' + ((_OPENMP/1000)%10)),
                            ('0' + ((_OPENMP/100)%10)),
                            ('0' + ((_OPENMP/10)%10)),
                            ('0' + ((_OPENMP/1)%10)),
                            ']', '\0' };
int main(void)
{
  puts(ompver_str);
  return 0;
}
