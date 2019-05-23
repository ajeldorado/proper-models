library_file = 'ffti'
input_file = 'ffti'
exported_routines = ['ffti']

os = strupcase(!version.os)

if ( os eq 'LINUX' ) then begin
	intel_dir = '/opt/intel/'
	mkl_dir = intel_dir + 'mkl/lib/intel64'
	lib_dir = intel_dir + 'lib/intel64'
	if ( file_search(lib_dir) eq '' ) then begin
		lib_dir = intel_dir + 'lib/intel64_lin'
		if ( file_search(lib_dir) eq '' ) then begin
			print, 'Cannot find library directory equivalent to /opt/intel/lib/intel64'
			print, 'Cannot compile Intel MKL interface'
			return 
		endif
	endif
	extra_lflags = '-L' + mkl_dir + ' -L' + lib_dir + ' -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lmkl_avx -liomp5 -lpthread -lm -ldl'
	extra_cflags = '-DMKL_ILP64 -m64 -I' + intel_dir + 'mkl/include'
   	make_dll, input_file, library_file, exported_routines, EXTRA_CFLAGS=extra_cflags, $
             INPUT_DIRECTORY='.', COMPILE_DIRECTORY='.', OUTPUT_DIRECTORY='.', EXTRA_LFLAGS=extra_lflags
endif else begin
	print, 'Unsupported OS.'
endelse

end
