#include "matrix.h"

/* this type represents a square sub-matrix */
typedef struct {
  size_t f_row;// first row of the sub-matrix
  size_t f_col;// first column of the sub-matrix
  float **matrix;// pointer to the real matrix
  size_t n; // number of rows and columns of the
             // sub-matrix
} sub_matrix_type;

/*
 * this function performs the element-wise
 * subtraction of B from A and put the resulting
 * sub-matrix in C
 */
void sub_matrix_subtract(sub_matrix_type C,
                         const sub_matrix_type A,
                         const sub_matrix_type B)
{
    for (size_t y = 0;y < C.n;y++) {
      for (size_t x = 0;x < C.n;x++) {
          C.matrix[y+C.f_row][x+C.f_col] =
               A.matrix[y+A.f_row][x+A.f_col] -
               B.matrix[y+B.f_row][x+B.f_col];
      }
    }
}

/*
 * this function performs the element-wise
 * sum of A and B and put the resulting
 * sub-matrix in C
 */
void sub_matrix_sum(sub_matrix_type C,
                    const sub_matrix_type A,
                    const sub_matrix_type B)
{
    for (size_t y = 0;y < C.n;y++) {
      for (size_t x = 0;x < C.n;x++) {
          C.matrix[y+C.f_row][x+C.f_col] =
               A.matrix[y+A.f_row][x+A.f_col] +
               B.matrix[y+B.f_row][x+B.f_col];
      }
    }
}

/*
 * this function returns one of the four blocks
 * of the sub-matrix A.
 */
sub_matrix_type get_block(sub_matrix_type A,
                          const size_t i, const size_t j)
{
    sub_matrix_type block;

    block.n = A.n/2;
    block.f_row = A.f_row + i*block.n;
    block.f_col = A.f_col + j*block.n;
    block.matrix = A.matrix;

    return block;
}

/*
 * this function embeds a matrix in a sub-matrix of
 * the same size
 */
sub_matrix_type embed_matrix(float **A, const size_t n)
{
  sub_matrix_type sm = {0, 0, A, n};

  return sm;
}

/*
 * this function builds a matrix and embeds it in
 * a sub-matrix of the same size
 */
sub_matrix_type build_sub_matrix(const size_t n)
{
  return embed_matrix(allocate_matrix(n, n), n);
}

/*
 * this function implements the naive algorithm
 * for matrix multiplication between sub-matrixes.
 * The result is placed in the sub-matrix C.
 */
void naive_aux(sub_matrix_type C,
               const sub_matrix_type A,
               const sub_matrix_type B)
{
  for (size_t y = 0;y < A.n;y++) {
    for (size_t x = 0;x < A.n;x++) {
      float value = 0.0;
      for (size_t z = 0;z < A.n;z++) {
        value += (A.matrix[y + A.f_row][z + A.f_col]*
                  B.matrix[z + B.f_row][x + B.f_col]);
      }

      C.matrix[y + C.f_row][x + C.f_col] = value;
    }
  }
}

/*
 * this function implements the Strassen algorithm
 * for matrix multiplication between sub-matrixes.
 * The result is placed in the sub-matrix C.
 */
void strassen_aux(sub_matrix_type C,
                  const sub_matrix_type A,
                  const sub_matrix_type B)
{
    // This is the base case of the
    // recursive algorithm
    if (A.n <= (1<<5)) {
        naive_aux(C, A, B);

        return;
    }

    // initialize S's and P's matrices
    sub_matrix_type S[10];
    for (size_t i = 0;i < 10;i++) {
      S[i] = build_sub_matrix(A.n/2);
    }

    sub_matrix_type P[7];
    for (size_t i = 0;i < 7;i++) {
      P[i] = build_sub_matrix(A.n/2);
    }

    // build A's and B's blocks
    sub_matrix_type A11, A12, A21, A22;

    A11 = get_block(A, 0, 0);
    A12 = get_block(A, 0, 1);
    A21 = get_block(A, 1, 0);
    A22 = get_block(A, 1, 1);

    sub_matrix_type B11, B12, B21, B22;

    B11 = get_block(B, 0, 0);
    B12 = get_block(B, 0, 1);
    B21 = get_block(B, 1, 0);
    B22 = get_block(B, 1, 1);

    // compute the matrices S by
    // performing sums and subtractions

    // S1 = B12 - B22
    sub_matrix_subtract(S[0], B12, B22);

    // S2 = A11 + A12
    sub_matrix_sum(S[1], A11, A12);

    // S3 = A21 + A22
    sub_matrix_sum(S[2], A21, A22);

    // S4 = B21 - B11
    sub_matrix_subtract(S[3], B21, B11);

    // S5 = A11 + A22
    sub_matrix_sum(S[4], A11, A22);

    // S6 = B11 + B22
    sub_matrix_sum(S[5], B11, B22);

    // S7 = A12 - A22
    sub_matrix_subtract(S[6], A12, A22);

    // S8 = B21 + B22
    sub_matrix_sum(S[7], B21, B22);

    // S9 = A11 - A21
    sub_matrix_subtract(S[8], A11, A21);

    // S10 = B11 + B12
    sub_matrix_sum(S[9], B11, B12);

    // compute the matrices P by
    // calling recursively strassen_aux

    // P1 = A11 x S1
    strassen_aux(P[0], A11, S[0]);

    // P2 = S2 x B22
    strassen_aux(P[1], S[1], B22);

    // P3 = S3 x B11
    strassen_aux(P[2], S[2], B11);

    // P4 = A22 x S4
    strassen_aux(P[3], A22, S[3]);

    // P5 = S5 x S6
    strassen_aux(P[4], S[4], S[5]);

    // P6 = S7 x S8
    strassen_aux(P[5], S[6], S[7]);

    // P7 = S9 x S10
    strassen_aux(P[6], S[8], S[9]);

    // build C's blocks
    sub_matrix_type C11, C12, C21, C22;

    C11 = get_block(C, 0, 0);
    C12 = get_block(C, 0, 1);
    C21 = get_block(C, 1, 0);
    C22 = get_block(C, 1, 1);

    // compute values in C's blocks by
    // performing sums and subtractions

    // C11 = P5 + P4 - P2 + P6
    sub_matrix_sum(C11, P[4], P[3]);
    sub_matrix_subtract(C11, C11, P[1]);
    sub_matrix_sum(C11, C11, P[5]);

    // C12 = P1 + P2
    sub_matrix_sum(C12, P[0], P[1]);

    // C21 = P3 + P4
    sub_matrix_sum(C21, P[2], P[3]);

    // C22 = P5 + P1 - P3 - P7
    sub_matrix_sum(C22, P[4], P[0]);
    sub_matrix_subtract(C22, C22, P[2]);
    sub_matrix_subtract(C22, C22, P[6]);

    // free S's and P's matrices
    for (size_t i = 0;i < 10;i++) {
      deallocate_matrix(S[i].matrix, S[i].n);
    }

    for (size_t i = 0;i < 7;i++) {
      deallocate_matrix(P[i].matrix, P[i].n);
    }
}


/*
 * this functions is exclusively meant to provide an
 * easy to use API
 */
void strassen_matrix_multiplication(float **C, float const *const *const A,
                                    float const *const *const B, size_t n)
{

  strassen_aux(embed_matrix(C, n),
               embed_matrix((float **)A, n),
               embed_matrix((float **)B, n));

}

