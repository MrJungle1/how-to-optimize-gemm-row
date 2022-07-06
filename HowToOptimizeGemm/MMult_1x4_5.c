/* Create macros so that the matrices are stored in column-major order */

#define A(i,j) a[ (i)*lda + (j) ]
#define B(i,j) b[ (i)*ldb + (j) ]
#define C(i,j) c[ (i)*ldc + (j) ]


/* Routine for computing C = A * B + C */

// void AddDot( int, double *, double *, int, double * );
void AddDot1x4( int, double *, int,  double *, int, double *, int );


void MY_MMult( int m, int n, int k, double *a, int lda, 
                                    double *b, int ldb,
                                    double *c, int ldc )
{
  int i, j;

  for ( i=0; i<m; i+=4 ){        /* Loop over the columns of C */
    for ( j=0; j<n; j+=1 ){        /* Loop over the rows of C */
      /* Update the C( i,j ) with the inner product of the ith row of A
	 and the jth column of B */
      AddDot1x4( k, &A( i,0 ), lda, &B( 0,j ), ldb, &C( i,j ), ldc );

    }
  }
}

void AddDot1x4( int k, double *a, int lda,  double *b, int ldb, double *c, int ldc )
{

  /* So, this routine computes four elements of C: 

           C( 0, 0 ), C( 0, 1 ), C( 0, 2 ), C( 0, 3 ).  

     Notice that this routine is called with c = C( i, j ) in the
     previous routine, so these are actually the elements 

           C( i, j ), C( i, j+1 ), C( i, j+2 ), C( i, j+3 ) 
	  
     in the original matrix C.

     In this version, we "inline" AddDot */ 

  int p;

  // AddDot( k, &A( 0,0 ), &B( 0,0 ), ldb, &C( 0,0 ) );
  // AddDot( k, &A( 1,0 ), &B( 0,0 ), ldb, &C( 1,0 ) );
  // AddDot( k, &A( 2,0 ), &B( 0,0 ), ldb, &C( 2,0 ) );
  // AddDot( k, &A( 3,0 ), &B( 0,0 ), ldb, &C( 3,0 ) );
  for ( p=0; p<k; p++ ){
    C( 0, 0 ) += A( 0, p ) * B( p, 0 );   
    C( 1, 0 ) += A( 1, p ) * B( p, 0 );  
    C( 2, 0 ) += A( 2, p ) * B( p, 0 );   
    C( 3, 0 ) += A( 3, p ) * B( p, 0 );      
  }




}

// /* Create macro to let X( i ) equal the ith element of x */

// #define Y(i) y[ (i)*incy ]

// void AddDot( int k, double *x,   double *y, int incy, double *gamma )
// {
//   /* compute gamma := x' * y + gamma with vectors x and y of length n.

//      Here x starts at location x with increment (stride) incx and y starts at location y and has (implicit) stride of 1.
//   */
 
//   int p;

//   for ( p=0; p<k; p++ ){
//     *gamma += x[ p ] * Y( p );     
//   }
// }
