      subroutine xtcerror(message, fatal)
      character*40 message
      integer*4 fatal
 601  format(/'XTC Error: ',2x,a,/)
      write(6,601) message
      if(fatal.eq.1) stop
      end
