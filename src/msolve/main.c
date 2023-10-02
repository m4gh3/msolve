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

#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpq_mpoly.h>

#include "libmsolve.c"

#define DEBUGGB 0
#define DEBUGBUILDMATRIX 0
#define IO_DEBUG 0


data_gens_ff_t *gens2msolve(fmpq_mpoly_t *polys, size_t n, char **vnames, fmpq_mpoly_ctx_t ctx )
{

	//slong nr_terms = 0;

	data_gens_ff_t *gens = (data_gens_ff_t *)(malloc(sizeof(data_gens_ff_t)));
	gens->cfs   = NULL;
  	gens->elim = 0;

	gens->nvars = fmpq_mpoly_ctx_nvars(ctx);
	gens->ngens = n;
	gens->nterms = 0;
	gens->field_char = 0;
	gens->vnames = NULL;
	gens->change_var_order = -1;
	gens->linear_form_base_coef = 0;
	gens->rand_linear = 0;
	gens->random_linear_form = malloc(sizeof(int32_t)*(gens->nvars));
	gens->elim = 0;

	gens->vnames = (char **)malloc( (gens->nvars) * (sizeof(char *)) );

	for(size_t i=0; i < gens->nvars; i++ )
	{
		/*size_t s = strlen(vnames[i])+1;
		gens->vnames[i] = (char *)malloc(s);
		strcpy(gens->vnames[i], vnames[i] );*/
		gens->vnames[i] = vnames[i];
	}

	gens->lens = (int32_t *)malloc((unsigned long)(gens->ngens) * sizeof(int32_t));

	ulong exps[gens->nvars];

	for(size_t i=0; i < n; i++ )
	{
		gens->lens[i] = fmpq_mpoly_length(polys[i], ctx );
		gens->nterms +=  gens->lens[i];
	}

	gens->exps = (int32_t *) malloc(gens->nterms * gens->nvars * sizeof(int32_t));
	gens->mpz_cfs = (mpz_t **)(malloc(sizeof(mpz_t *) * 2 * gens->nterms));

	for(int i=0; i < 2*gens->nterms; i++ )
	{
		gens->mpz_cfs[i] = (mpz_t *)malloc(sizeof(mpz_t));
		mpz_init(*(gens->mpz_cfs[i]));
	}

	for(size_t i=0,p=0; i < n; i++ )
	{
		fmpq_t c;
		fmpq_init(c);
		for(size_t j=0; j < gens->lens[i]; j++, p++ )
		{

			fmpq_mpoly_get_term_coeff_fmpq(c, polys[i], j, ctx );
			fmpq_get_mpz_frac(gens->mpz_cfs[2*p][0], gens->mpz_cfs[2*p+1][0], c );
			fmpq_mpoly_get_term_exp_ui(exps/*gens->exps + gens->nvars*p*/, polys[i], j, ctx );

			for(size_t k=0; k < gens->nvars; k++ )
				(gens->exps + gens->nvars*p)[k] = (int32_t) exps[k];

		}
	}

	return gens;

}


int msolve_from_fmpq_mpolys(fmpq_mpoly_t *polys, size_t n, char **vnames, fmpq_mpoly_ctx_t ctx )
{

    /* timinigs */
    double st0 = cputime();
    double rt0 = realtime();

    /**
      We get values from the command line.
     **/
    int32_t la_option             = 2; // by default
    int32_t use_signatures        = 0;
    int32_t nr_threads            = 1;
    int32_t info_level            = 0;
    int32_t initial_hts           = 17;
    int32_t max_pairs             = 0;
    int32_t elim_block_len        = 0;
    int32_t update_ht             = 0;
    int32_t generate_pbm          = 0;
    int32_t reduce_gb             = 1;
    int32_t print_gb              = 0;
    int32_t genericity_handling   = 2;
    int32_t saturate              = 0;
    int32_t colon                 = 0;
    int32_t normal_form           = 0;
    int32_t normal_form_matrix    = 0;
    int32_t is_gb                 = 0;
    int32_t get_param             = 0;
    int32_t precision             = 128;
    int32_t refine                = 0; /* not used at the moment */
    int32_t isolate               = 0; /* not used at the moment */

    files_gb *files = malloc(sizeof(files_gb));
    files->in_file = NULL;
    files->bin_file = NULL;
    files->out_file = NULL;
    files->bin_out_file = NULL;
    /**
       We get from files the requested data. 
    **/
    data_gens_ff_t *gens = gens2msolve(polys, n, vnames, ctx );//get_data_from_flint();
#ifdef IODEBUG
    display_gens(stdout, gens);
#endif
    
    /* data structures for parametrization */
    param_t *param  = NULL;
    mpz_param_t mpz_param;
    mpz_param_init(mpz_param);
    
    long nb_real_roots      = 0;
    interval *real_roots    = NULL;
    real_point_t *real_pts  = NULL;

    /* main msolve functionality */
    int ret = core_msolve(la_option, use_signatures, nr_threads, info_level,
                          initial_hts, max_pairs, elim_block_len, update_ht,
                          generate_pbm, reduce_gb, print_gb, get_param,
                          genericity_handling, saturate, colon, normal_form,
                          normal_form_matrix, is_gb, precision, 
                          files, gens,
            &param, &mpz_param, &nb_real_roots, &real_roots, &real_pts);

    /* free parametrization */
    free(param);
    mpz_param_clear(mpz_param);


    if (nb_real_roots > 0) {
        for(long i = 0; i < nb_real_roots; i++){
          real_point_clear(real_pts[i]);
          mpz_clear(real_roots[i].numer);
        }
        free(real_pts);
    }
    free(real_roots);

    /* timings */
    if (info_level > 0) {
        double st1 = cputime();
        double rt1 = realtime();
        fprintf(stderr, "-------------------------------------------------\
-----------------------------------\n");
        fprintf(stderr, "msolve overall time  %13.2f sec (elapsed) / %5.2f sec (cpu)\n",
                rt1-rt0, st1-st0);
        fprintf(stderr, "-------------------------------------------------\
-----------------------------------\n");
    }
    //free_data_gens(gens);
    /* for(long i = 0; i < gens->nvars; i++){
        free(gens->vnames[i]);
    }
    free(gens->vnames);*/
    free(gens->lens);
    free(gens->cfs);
    free(gens->exps);
    free(gens->random_linear_form);
    free(files);
    return ret;
}

int main()
{
	fmpq_mpoly_ctx_t ctx;
	fmpq_mpoly_ctx_init(ctx, 6, ORD_DEGREVLEX );

  data_gens_ff_t *gens = (data_gens_ff_t *)(malloc(sizeof(data_gens_ff_t)));

  gens->nvars = fmpq_mpoly_ctx_nvars(ctx);
  
  char **vnames = (char **)malloc((gens->nvars) * sizeof(char *));

  const char *ring_gen_names[] = {"w0", "w1", "w2", "w3", "l0", "l1" };
  for(int i=0; i < 6; i++ )
    vnames[i] = ring_gen_names[i];
  	
  fmpq_mpoly_t P[6];
  
  fmpq_mpoly_init(P[0], ctx );
  fmpq_mpoly_init(P[1], ctx );
  fmpq_mpoly_init(P[2], ctx );
  fmpq_mpoly_init(P[3], ctx );
  fmpq_mpoly_init(P[4], ctx );
  fmpq_mpoly_init(P[5], ctx );
  
  
  fmpq_mpoly_set_str_pretty(P[0], "8*w0*w2^2 + 16*w1*w2*w3 + 8*w0*w3^2 + 2*w0*l0 - 8*w3", ring_gen_names, ctx );
  fmpq_mpoly_set_str_pretty(P[1], "8*w1*w2^2 + 16*w0*w2*w3 + 8*w1*w3^2 + 2*w1*l0 - 8*w2", ring_gen_names, ctx );
  fmpq_mpoly_set_str_pretty(P[2], "8*w0^2*w2 + 8*w1^2*w2 + 16*w0*w1*w3 + 2*w2*l1 - 8*w1", ring_gen_names, ctx );
  fmpq_mpoly_set_str_pretty(P[3], "16*w0*w1*w2 + 8*w0^2*w3 + 8*w1^2*w3 + 2*w3*l1 - 8*w0", ring_gen_names, ctx );
  fmpq_mpoly_set_str_pretty(P[4], "w0^2 + w1^2 - 1", ring_gen_names, ctx );
  fmpq_mpoly_set_str_pretty(P[5], "w2^2 + w3^2 - 1", ring_gen_names, ctx );

  return msolve_from_fmpq_mpolys(P, 6, ring_gen_names, ctx );

}
