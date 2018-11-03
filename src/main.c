/*
 * main.c: entry point for meteodetect
 * Creation date: Fri Nov  2 22:45:20 2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>

#include <meteodetect.h>

#include <sigutils/iir.h>
#include <sigutils/ncqo.h>
#include <sigutils/log.h>
#include <sigutils/sampling.h>

#include <assert.h>

#define LPF1_CUTOFF 300. /* In Hz */
#define LPF2_CUTOFF 50.  /* In Hz */
#define POWER_RATIO (LPF2_CUTOFF / LPF1_CUTOFF)
#define POWER_THRESHOLD (2 * POWER_RATIO)

#define MIN_CHIRP_DURATION 0.07

struct meteo_detect {
  SUFLOAT fs;          /* Samp rate */
  SUSCOUNT n;
  su_iir_filt_t lpf1; /* Low pass filter 1. Used to detect noise power */
  su_iir_filt_t lpf2; /* Low pass filter 2. Used to isolate chirps */
  su_ncqo_t lo;
  SUFLOAT alpha;
  SUFLOAT n0; /* Noise power */
  SUFLOAT s0; /* Signal power */

  SUSCOUNT powerbuf_len;
  SUSCOUNT p;
  SUFLOAT *powerbuf;
  SUCOMPLEX *delay_line;

  SUFLOAT  energy_thres;
  SUBOOL in_chirp;
  SUSCOUNT chirp_len;
  SUSCOUNT chirp_remaining;
  SUCOMPLEX prev;

  FILE *ofp;
};

typedef struct meteo_detect meteo_detect_t;

void
meteo_detect_destroy(meteo_detect_t *detect)
{
  if (detect->ofp != NULL)
    fclose(detect->ofp);

  if (detect->powerbuf != NULL)
    free(detect->powerbuf);

  if (detect->delay_line != NULL)
    free(detect->delay_line);

  su_iir_filt_finalize(&detect->lpf1);
  su_iir_filt_finalize(&detect->lpf2);

  free(detect);
}

void
meteo_detect_feed(meteo_detect_t *md, SUCOMPLEX x)
{
  SUCOMPLEX y;
  SUCOMPLEX out;
  SUFLOAT ratio;
  SUFLOAT integral;
  SUBOOL valid_data;
  unsigned int i;
  SUSCOUNT seconds;

  y = su_iir_filt_feed(&md->lpf1, x * SU_C_CONJ(su_ncqo_read(&md->lo)));
  md->n0 += md->alpha * (SU_C_REAL(y * SU_C_CONJ(y)) - md->n0);

  y = su_iir_filt_feed(&md->lpf2, y);
  md->s0 += md->alpha * (SU_C_REAL(y * SU_C_CONJ(y)) - md->s0);

  ratio = md->s0 / md->n0;

  /* Put SNR value */
  md->powerbuf[md->p] = ratio;

  /* Keep demodulated sample in delay line */
  md->delay_line[md->p] = y * SU_C_CONJ(md->prev);

  if (++md->p == md->powerbuf_len)
    md->p = 0;

  /* md->p now points to the OLDEST sample */

  /* Compute integral */
  integral = 0;
  for (i = 0; i < md->powerbuf_len; ++i)
    integral += md->powerbuf[i];

  if (md->in_chirp) {
    if (integral < md->energy_thres) {
      seconds = SU_FLOOR((md->n - md->chirp_len) / md->fs);
      printf(
          "Chirp of length %5d detected (at %02d:%02d:%02d)\n",
          md->chirp_len,
          seconds / 3600,
          (seconds / 60) % 60,
          seconds % 60);
      md->in_chirp = SU_FALSE;
    } else {
      ++md->chirp_len;
    }
  } else {
    if (integral >= md->energy_thres) {
      md->in_chirp = SU_TRUE;
      md->chirp_len = md->chirp_remaining = md->powerbuf_len;
    }
  }

  /* If not inside chirp: don't save anything */
  if (!md->in_chirp && md->chirp_remaining > 0)
    --md->chirp_remaining;

  valid_data = md->chirp_remaining != 0;

  out = valid_data + (valid_data ? I * SU_C_ARG(md->delay_line[md->p]) : 0);

  fwrite(&out, sizeof(SUCOMPLEX), 1, md->ofp);

  md->prev = y;
  ++md->n;
}

meteo_detect_t *
meteo_detect_new(SUFLOAT fs, SUFLOAT fc, const char *ofile)
{
  meteo_detect_t *new = NULL;
  unsigned int len;

  SU_TRYCATCH(new = calloc(1, sizeof (meteo_detect_t)), goto fail);

  new->fs = fs;
  new->alpha = 1 - SU_EXP(-1. / (fs * MIN_CHIRP_DURATION));

  su_ncqo_init(&new->lo, SU_ABS2NORM_FREQ(fs, fc));

  SU_TRYCATCH(
      su_iir_bwlpf_init(&new->lpf1, 5, SU_ABS2NORM_FREQ(fs, LPF1_CUTOFF)),
      goto fail);

  SU_TRYCATCH(
      su_iir_bwlpf_init(&new->lpf2, 4, SU_ABS2NORM_FREQ(fs, LPF2_CUTOFF)),
      goto fail);

  new->powerbuf_len = SU_CEIL(fs * MIN_CHIRP_DURATION);
  new->energy_thres = POWER_THRESHOLD * new->powerbuf_len;

  SU_TRYCATCH(
      new->powerbuf   = calloc(sizeof(SUFLOAT), new->powerbuf_len),
      goto fail);

  SU_TRYCATCH(
      new->delay_line = calloc(sizeof(SUCOMPLEX), new->powerbuf_len),
      goto fail);

  SU_TRYCATCH(new->ofp = fopen(ofile, "wb"), goto fail);

  return new;

fail:
  if (new != NULL)
    meteo_detect_destroy(new);

  return new;
}

int
main (int argc, char *argv[], char *envp[])
{
  FILE *fp;
  meteo_detect_t *md;

  SUCOMPLEX sample;

  assert(sizeof(SUCOMPLEX) == sizeof(complex float));

  if (argc != 2) {
    fprintf(stderr, "Usage: %s <file>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  if ((fp = fopen(argv[1], "rb")) == NULL) {
    fprintf(stderr, "%s: cannot open `%s': %s\n", argv[0], argv[1], strerror(errno));
    exit(EXIT_FAILURE);
  }

  SU_TRYCATCH(
      md = meteo_detect_new(8000, 1000, "detect.raw"),
      exit(EXIT_FAILURE));

  while (fread(&sample, sizeof(sample), 1, fp) == 1)
    meteo_detect_feed(md, sample);

  meteo_detect_destroy(md);

  fclose(fp);

  
  return 0;
}

