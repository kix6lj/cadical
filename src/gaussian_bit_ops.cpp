#include "internal.hpp"
#include <immintrin.h>

#define WIDTH 256
//#define AVX512

#define DIV_CEIL(x, y) ((x) + (y)-1) / (y)

#define DIV_FLOOR(x, y) (x) / (y)

namespace CaDiCaL {
void GaussianDense::XOR(char *i, char *j) { XOR(i, i, j); }

void GaussianDense::XOR(char *r, char *i, char *j) {
  int num_bound = DIV_FLOOR(num_bytes_per_eqn, WIDTH >> 3);
  int bytes = WIDTH >> 3;
  char *offset_i = i;
  char *offset_j = j;
  char *offset_r = r;
  for (int i = 0; i < num_bound; ++i, offset_i += bytes, offset_j += bytes, offset_r += bytes) {
    __m256i coeff_i = _mm256_load_si256((__m256i *)offset_i);
    __m256i coeff_j = _mm256_load_si256((__m256i *)offset_j);
    __m256i result = _mm256_castps_si256(_mm256_xor_ps(
        _mm256_castsi256_ps(coeff_i), _mm256_castsi256_ps(coeff_j)));
    _mm256_storeu_si256((__m256i *)offset_r, result);
  }
}

void GaussianDense::OR(char *i, char *j) { OR(i, i, j); }

void GaussianDense::OR(char *r, char *i, char *j) {
  int num_bound = DIV_FLOOR(num_bytes_per_eqn, WIDTH >> 3);
  int bytes = WIDTH >> 3;
  char *offset_i = i;
  char *offset_j = j;
  char *offset_r = r;
  for (int i = 0; i < num_bound; ++i, offset_i += bytes, offset_j += bytes, offset_r += bytes) {
    __m256i coeff_i = _mm256_load_si256((__m256i *)offset_i);
    __m256i coeff_j = _mm256_load_si256((__m256i *)offset_j);
    __m256i result = _mm256_castps_si256(_mm256_or_ps(
        _mm256_castsi256_ps(coeff_i), _mm256_castsi256_ps(coeff_j)));
    _mm256_storeu_si256((__m256i *)offset_r, result);
  }
}

void GaussianDense::AND(char *i, char *j) { AND(i, i, j); }

void GaussianDense::AND(char *r, char *i, char *j) {
  int num_bound = DIV_FLOOR(num_bytes_per_eqn, WIDTH >> 3);
  int bytes = WIDTH >> 3;
  char *offset_i = i;
  char *offset_j = j;
  char *offset_r = r;
  for (int i = 0; i < num_bound; ++i, offset_i += bytes, offset_j += bytes, offset_r += bytes) {
    __m256i coeff_i = _mm256_load_si256((__m256i *)offset_i);
    __m256i coeff_j = _mm256_load_si256((__m256i *)offset_j);
    __m256i result = _mm256_castps_si256(_mm256_and_ps(
        _mm256_castsi256_ps(coeff_i), _mm256_castsi256_ps(coeff_j)));
    _mm256_storeu_si256((__m256i *)offset_r, result);
  }
}

void GaussianDense::IAND(char *i, char *j) { IAND(i, i, j); }
  
void GaussianDense::IAND(char *r, char *i, char *j) {
  int num_bound = DIV_FLOOR(num_bytes_per_eqn, WIDTH >> 3);
  int bytes = WIDTH >> 3;
  char *offset_i = i;
  char *offset_j = j;
  char *offset_r = r;
  __m256i mask = _mm256_set1_epi32(-1);
  for (int i = 0; i < num_bound; ++i, offset_i += bytes, offset_j += bytes, offset_r += bytes) {
    __m256i coeff_i = _mm256_load_si256((__m256i *)offset_i);
    __m256i coeff_j = _mm256_load_si256((__m256i *)offset_j);
    __m256i result = _mm256_castps_si256(
        _mm256_and_ps(_mm256_castsi256_ps(coeff_i),
                      _mm256_xor_ps(_mm256_castsi256_ps(coeff_j),
				    _mm256_castsi256_ps(mask))));
    _mm256_storeu_si256((__m256i *)offset_r, result);
  }
}

int GaussianDense::count(char *ptr) {
#ifdef AVX512
  __m256i sum = _mm256_setzero_si256();;
  char *offset = ptr;
  int num_bound = DIV_FLOOR(num_bytes_per_eqn, WIDTH >> 3);
  int bytes = WIDTH >> 3;

  for (int i = 0; i < num_bound; ++i, offset += bytes) {
    __m256i coeff = _mm256_load_si256((__m256i *)offset);
    __m256i num_ones = _mm256_popcnt_epi64(coeff);
    sum = _mm256_add_epi64(sum, num_ones);
  }

  int64_t *sum_ptr = (int64_t *)&sum;

  return (int)(sum_ptr[0] + sum_ptr[1] + sum_ptr[2] + sum_ptr[3]) - !!get_bit(ptr, 0);

#else

  int sum = 0;
  char *offset = ptr;
  int num_bound = DIV_FLOOR(num_bytes_per_eqn, 64 >> 3);
  int bytes = 64 >> 3;

  for (int i = 0; i < num_bound; ++i, offset += bytes) {
    auto coeff = *((int64_t *)offset);
    sum += __builtin_popcountll(coeff);
  }

  return sum - !!get_bit(ptr, 0);
#endif
}

int GaussianDense::msb_position(char *ptr) {
#ifdef AVX512
  int num_bound = DIV_FLOOR(num_bytes_per_eqn, WIDTH >> 3);
  int bytes = 64 >> 3;
  char *offset = ptr + (num_bound - 1) * bytes;

  int pos = 0;
  for (int i = num_bound - 1; i >= 0; --i, offset -= bytes) {
    __m256i coeff = _mm256_load_si256((__m256i *)offset);
    __m256i lzcnt = _mm256_lzcnt_epi64(coeff);
    int all_zero = _mm256_testz_si256(coeff, coeff);

    if (all_zero == 0) {
      int64_t *lz_ptr = (int64_t *)&lzcnt;
      pos += lz_ptr[3] + (lz_ptr[3] == 64 ? lz_ptr[2] : 0) +
             (lz_ptr[3] == 64 && lz_ptr[2] == 64 ? lz_ptr[1] : 0) +
             (lz_ptr[3] == 64 && lz_ptr[2] == 64 && lz_ptr[1] == 64 ? lz_ptr[0]
                                                                    : 0);
      break;
    }
    pos += WIDTH;
  }

  return (num_bytes_per_eqn << 3) - pos - 1;
#else
  for (int i = num_vars; i > 0; --i)
    if (ptr[i >> 3] & (1 << (i & 7)))
      return i;
  return 0;
#endif
}

int GaussianDense::lsb_position(char *ptr) {
  for (int i = 1; i <= num_vars; ++i)
    if (ptr[i >> 3] & (1 << (i & 7)))
      return i;
  return 0;
}

char GaussianDense::get_bit(char *ptr, int i) {
  return ptr[i >> 3] & (1 << (i & 7));
}

void GaussianDense::set_bit(char *ptr, int i) {
  ptr[i >> 3] |= (1 << (i & 7));
}

void GaussianDense::unset_bit(char *ptr, int i) {
  ptr[i >> 3] &= ~(1 << (i & 7));
}

void GaussianDense::print(char *ptr) {
  for (int i = 1; i <= num_vars; ++i)
    if (get_bit(ptr, i))
      printf("%d ", var_id[i - 1]);
  printf("\n");
}

}; // namespace CaDiCaL
