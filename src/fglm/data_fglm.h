/* This file is part of msolve.
 *
 * msolve is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * msolve is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with msolve.  If not, see <https://www.gnu.org/licenses/>
 *
 * Authors:
 * Jérémy Berthomieu
 * Christian Eder
 * Mohab Safey El Din */

#include<stdint.h>
#include <flint/flint.h>
#include <flint/longlong.h>
#include <flint/mpn_extras.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/ulong_extras.h>


typedef uint32_t szmat_t;
typedef uint32_t CF_t;
typedef uint64_t CF_l_t;
typedef uint32_t mod_t;
typedef int32_t nvars_t;
/* typedef __uint128_t CF_L_t; */

/**

M:=Matrix([[0, 1, 0, 0, 0],[ 0, 0, 0, 1, 0],[1, 2, 4, 5, 6],[0, 0, 0, 0, 1],[7, 8, 1, 9, 4]]);

   Une matrice de la forme

   0 1 0 0 0
   0 0 0 1 0
   1 2 4 5 6
   0 0 0 0 1
   7 8 1 9 4

   est code de la forme
   ncols  = 5
   nrows = 2
   dense_mat -> 1 2 4 5 6 7 8 1 9
   triv_idx -> 0 1 3 (lgr donnee par ncols-nrows)
   triv_pos -> 1 3 4
   dense_idx -> 2 4 (lgr donnee par nrows)
 **/

typedef struct{
  CF_t charac;
  szmat_t ncols; //dimension du quotient
  szmat_t nrows; //nbre de lignes non triviales
  CF_t *dense_mat; // matrice nrows lignes et ncols colonnes (elements donnes par lignes)
  szmat_t *triv_idx; //tableau d'indices des lignes ne contenant que des 0 et un 1
  szmat_t *triv_pos; //position des 1
  szmat_t *dense_idx; //position des lignes non triviales (qui constituent donc
                      //dense_mat)
  szmat_t *dst; //pour la gestion des lignes "denses" mais avec un bloc de zero a la fin
} sp_matfglm_t;

typedef struct{
  CF_t charac;
  szmat_t ncols; //dimension du sev du quotient
  szmat_t nrows; //nbre de lignes non triviales
  szmat_t nzero; //nbre de lignes nulles
  CF_t *dense_mat; // matrice nrows lignes et ncols colonnes (elements donnes par lignes)
  szmat_t *triv_idx; //tableau d'indices des lignes ne contenant que des 0 et un 1
  szmat_t *triv_pos; //position des 1
  szmat_t *dense_idx; //position des lignes non triviales (qui constituent donc
                      //dense_mat)
  szmat_t *zero_idx; //tableau d'indices des lignes ne contenant que des 0
  szmat_t *dst; //pour la gestion des lignes "denses" mais avec un bloc de zero a la fin
} sp_matfglmcol_t;


#ifndef ALIGNED32
#define ALIGNED32 __attribute__((aligned(32)))
#endif
typedef struct{
  CF_t *vecinit ALIGNED32; //vecteur initial
  CF_t *res ALIGNED32; //va contenir les termes de suites qui nous interessent
  CF_t *vecmult ALIGNED32; //utilise pour la multiplication
  CF_t *vvec ALIGNED32;
  CF_l_t *vec_cache ALIGNED32; //useless
  mp_limb_t *pts;
} fglm_data_t;

/* typedef struct { */
/*   slong npoints; */
/*   nmod_poly_t R0, R1; */
/*   nmod_poly_t V0, V1; */
/*   nmod_poly_t qt, rt; //temporaries */
/*   nmod_poly_t points; */
/* } nmod_berlekamp_massey_struct; */
/* typedef nmod_berlekamp_massey_struct nmod_berlekamp_massey_t[1]; */


typedef struct{
  nmod_berlekamp_massey_t BMS;
  nmod_poly_t Z1;
  nmod_poly_t Z2;
  nmod_poly_t rZ1;
  nmod_poly_t rZ2;
  nmod_poly_t A;
  nmod_poly_t B;
  nmod_poly_t V;
  nmod_poly_t param;
  nmod_poly_factor_t sqf;
} fglm_bms_data_t;


typedef struct{
  mp_limb_t charac;
  nvars_t nvars;
  nmod_poly_t elim;
  nmod_poly_t denom;
  nmod_poly_t *coords;
} param_t;


