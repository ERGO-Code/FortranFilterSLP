      integer global_ifail
      integer global_int_arg
      double precision global_double_arg
      character*500 global_err_msg


      common /retc/ global_ifail,
     &              global_err_msg,
     &              global_int_arg(2),
     &              global_double_arg(4)


c     
cHM     Values for global ifail
cHM      0 ok
cHM      1 out of memory 1        No arguments
cHM      2 out of memory 2        2 x int arguments
cHM      3 presolve infeasible 1  2 x int + 2 x double
cHM      4 presolve infeasible 2  2 x int
cHM      5 presolve infeasible 3  2 x int + 4 x double
cHM      99 other                 No arguments
c
c     These are now suggested
c     Values for global ifail
c      0 ok
c      1 infeasible
c      20 out of passed memory
c      21 out of resources (dynamically allocated)
c      30 hardwired limit exceeded
c      50 bad usage
c      99 inconsistent state
c     
c     
