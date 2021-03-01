      subroutine ncerror(status, fatal)
      integer*4 status, fatal
#ifdef NETCDF
      character*80 nf_strerror
 601  format(/,2x,a,/)
      write(6,601) nf_strerror(status)
#endif
      if(fatal.eq.1) stop
      end
