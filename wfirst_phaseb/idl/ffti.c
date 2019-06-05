#include <stdio.h>
#include "idl_export.h"
#include "mkl_dfti.h"
#include "mkl_service.h"

typedef struct {
    double re;
    double im;
} mkl_double_complex;

typedef struct {
    float re;
    float im;
} mkl_float_complex;

int ffti( int argc, void *argv[] )
{
    	DFTI_DESCRIPTOR_HANDLE Desc_Handle = 0;
	void *in; 
	IDL_LONG 	data_type, direction;
    	MKL_LONG 	nx, ny, status, nthreads, lengths[2];

	data_type = *(IDL_LONG *)argv[0];
        in = (void *)argv[1];
        direction = *(IDL_LONG *)argv[2];
        nx = (MKL_LONG)(*(IDL_LONG *)argv[3]);
        ny = (MKL_LONG)(*(IDL_LONG *)argv[4]);
        nthreads = (MKL_LONG)(*(IDL_LONG *)argv[5]);

	lengths[0] = nx;
	lengths[1] = ny;

  	// if ( nthreads != 0 )
	//	mkl_set_num_threads_local( nthreads );

	if ( data_type == 9 )
		status = DftiCreateDescriptor( &Desc_Handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, lengths );
	else
		status = DftiCreateDescriptor( &Desc_Handle, DFTI_SINGLE, DFTI_COMPLEX, 2, lengths );
	if ( !DftiErrorClass( status, DFTI_NO_ERROR ) )
	{
        	printf("FAIL : DftiCreateDescriptor( &Desc_Handle, DFTI_DOUBLE,...\n");
        	return 0;
    	}

	if ( direction == -1 ) 
  		status = DftiSetValue( Desc_Handle, DFTI_FORWARD_SCALE, 1.0/(nx*ny) );
	else
  		status = DftiSetValue( Desc_Handle, DFTI_BACKWARD_SCALE, 1.0 );
    	if ( !DftiErrorClass( status, DFTI_NO_ERROR ) )
	{
        	printf("FAIL : DftiSetValue(Desc_Handle, DFTI_BACKWARD_SCALE, Scale)\n"); 
		goto FREE_DESCRIPTOR;
    	}

	status = DftiSetValue( Desc_Handle, DFTI_THREAD_LIMIT, nthreads );
    	if ( !DftiErrorClass( status, DFTI_NO_ERROR ) )
	{
        	printf("FAIL : DftiSetValue(Desc_Handle, DFTI_THREAD_LIMIT, nthreads)\n"); 
		goto FREE_DESCRIPTOR;
    	}

	status = DftiCommitDescriptor( Desc_Handle );
    	if ( !DftiErrorClass( status, DFTI_NO_ERROR ) )
	{
        	printf("FAIL : DftiCommitDescriptor( Desc_Handle )\n"); 
		goto FREE_DESCRIPTOR;
    	}

	if ( direction == -1 )
	{
		if ( data_type == 8 )
			status = DftiComputeForward( Desc_Handle, (mkl_double_complex *)in );
		else
			status = DftiComputeForward( Desc_Handle, (mkl_float_complex *)in );
	}
	else
	{
		if ( data_type == 8 ) 
    			status = DftiComputeBackward( Desc_Handle, (mkl_double_complex *)in );
		else
    			status = DftiComputeBackward( Desc_Handle, (mkl_float_complex *)in );
	}
    	if ( !DftiErrorClass( status, DFTI_NO_ERROR ) )
	{
        	printf("FAIL : DftiComputeForward( Desc_Handle, x_in)\n"); 
		goto FREE_DESCRIPTOR;
    	}

FREE_DESCRIPTOR:
    	status = DftiFreeDescriptor( &Desc_Handle );
    	if( !DftiErrorClass( status, DFTI_NO_ERROR ) )
        	printf("FAIL : DftiFreeDescriptor(&Desc_Handle)\n");

	return( 1 );
}

