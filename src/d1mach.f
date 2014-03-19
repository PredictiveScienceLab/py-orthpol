c
c
      double precision function d1mach(i)
c
c  Double-precision machine constants
c
c  d1mach( 1) = b**(emin-1), the smallest positive magnitude.
c
c  d1mach( 2) = b**emax*(1 - b**(-t)), the largest magnitude.
c
c  d1mach( 3) = b**(-t), the smallest relative spacing.
c
c  d1mach( 4) = b**(1-t), the largest relative spacing.
c
c  d1mach( 5) = log10(b)
c
c  To alter this function for a particular environment,
c  the desired set of data statements should be activated by
c  removing the c from column 1.
c  On rare machines a static statement may need to be added.
c  (But probably more systems prohibit it than require it.)
c
c  For IEEE-arithmetic machines (binary standard), one of the second
c  two sets of constants below should be appropriate.
c
c  Where possible, decimal, octal or hexadecimal constants are used
c  to specify the constants exactly.  Sometimes this requires using
c  equivalent integer arrays.  If your compiler uses half-word
c  integers by default (sometimes called integer*2), you may need to
c  change integer to integer*4 or otherwise instruct your compiler
c  to use full-word integers in the next 5 declarations.
c
      integer small(2)
      integer large(2)
      integer right(2)
      integer diver(2)
      integer log10(2)
      integer sc
c
      double precision dmach(5)
c
      equivalence (dmach(1),small(1))
      equivalence (dmach(2),large(1))
      equivalence (dmach(3),right(1))
      equivalence (dmach(4),diver(1))
      equivalence (dmach(5),log10(1))
c
c     machine constants for cdc cyber 205 and eta-10.
c
c     data small(1) / x'9000400000000000' /
c     data small(2) / x'8000000000000000' /
c
c     data large(1) / x'6FFF7FFFFFFFFFFF' /
c     data large(2) / x'6FD07FFFFFFFFFFF' /
c
c     data right(1) / x'FF74400000000000' /
c     data right(2) / x'8000000000000000' /
c
c     data diver(1) / x'FF75400000000000' /
c     data diver(2) / x'8000000000000000' /
c
c     data log10(1) / x'FFD04D104D427DE7' /
c     data log10(2) / x'FFA17DE623E2566B' /, sc/987/
c     
c     machine constants for ieee arithmetic machines, such as the at&t
c     3b series and motorola 68000 based machines (e.g. sun 3 and at&t
c     pc 7300), in which the most significant byte is stored first.
c
       data small(1),small(2) /    1048576,          0 /
       data large(1),large(2) / 2146435071,         -1 /
       data right(1),right(2) / 1017118720,          0 /
       data diver(1),diver(2) / 1018167296,          0 /
       data log10(1),log10(2) / 1070810131, 1352628735 /, sc/987/
c
c     machine constants for ieee arithmetic machines and 8087-based
c     micros, such as the ibm pc and at&t 6300, in which the least
c     significant byte is stored first.
c
c      data small(1),small(2) /          0,    1048576 /
c      data large(1),large(2) /         -1, 2146435071 /
c      data right(1),right(2) /          0, 1017118720 /
c      data diver(1),diver(2) /          0, 1018167296 /
c      data log10(1),log10(2) / 1352628735, 1070810131 /, sc/987/
c
c     machine constants for amdahl machines.
c
c      data small(1),small(2) /    1048576,          0 /
c      data large(1),large(2) / 2147483647,         -1 /
c      data right(1),right(2) /  856686592,          0 /
c      data diver(1),diver(2) /  873463808,          0 /
c      data log10(1),log10(2) / 1091781651, 1352628735 /, sc/987/
c
c     machine constants for the burroughs 1700 system.
c
c      data small(1) / zc00800000 /
c      data small(2) / z000000000 /
c
c      data large(1) / zdffffffff /
c      data large(2) / zfffffffff /
c
c      data right(1) / zcc5800000 /
c      data right(2) / z000000000 /
c
c      data diver(1) / zcc6800000 /
c      data diver(2) / z000000000 /
c
c      data log10(1) / zd00e730e7 /
c      data log10(2) / zc77800dc0 /, sc/987/
c
c     machine constants for the burroughs 5700 system.
c
c      data small(1) / o1771000000000000 /
c      data small(2) / o0000000000000000 /
c
c      data large(1) / o0777777777777777 /
c      data large(2) / o0007777777777777 /
c
c      data right(1) / o1461000000000000 /
c      data right(2) / o0000000000000000 /
c
c      data diver(1) / o1451000000000000 /
c      data diver(2) / o0000000000000000 /
c
c      data log10(1) / o1157163034761674 /
c      data log10(2) / o0006677466732724 /, sc/987/
c
c     machine constants for the burroughs 6700/7700 systems.
c
c      data small(1) / o1771000000000000 /
c      data small(2) / o7770000000000000 /
c
c      data large(1) / o0777777777777777 /
c      data large(2) / o7777777777777777 /
c
c      data right(1) / o1461000000000000 /
c      data right(2) / o0000000000000000 /
c
c      data diver(1) / o1451000000000000 /
c      data diver(2) / o0000000000000000 /
c
c      data log10(1) / o1157163034761674 /
c      data log10(2) / o0006677466732724 /, sc/987/
c
c     machine constants for ftn4 on the cdc 6000/7000 series.
c
c      data small(1) / 00564000000000000000b /
c      data small(2) / 00000000000000000000b /
c
c      data large(1) / 37757777777777777777b /
c      data large(2) / 37157777777777777774b /
c
c      data right(1) / 15624000000000000000b /
c      data right(2) / 00000000000000000000b /
c
c      data diver(1) / 15634000000000000000b /
c      data diver(2) / 00000000000000000000b /
c
c      data log10(1) / 17164642023241175717b /
c      data log10(2) / 16367571421742254654b /, sc/987/
c
c     machine constants for ftn5 on the cdc 6000/7000 series.
c
c      data small(1) / o"00564000000000000000" /
c      data small(2) / o"00000000000000000000" /
c
c      data large(1) / o"37757777777777777777" /
c      data large(2) / o"37157777777777777774" /
c
c      data right(1) / o"15624000000000000000" /
c      data right(2) / o"00000000000000000000" /
c
c      data diver(1) / o"15634000000000000000" /
c      data diver(2) / o"00000000000000000000" /
c
c      data log10(1) / o"17164642023241175717" /
c      data log10(2) / o"16367571421742254654" /, sc/987/
c
c     machine constants for convex c-1
c
c      data small(1),small(2) / '00100000'x, '00000000'x /
c      data large(1),large(2) / '7fffffff'x, 'ffffffff'x /
c      data right(1),right(2) / '3cc00000'x, '00000000'x /
c      data diver(1),diver(2) / '3cd00000'x, '00000000'x /
c      data log10(1),log10(2) / '3ff34413'x, '509f79ff'x /, sc/987/
c
c     machine constants for the cray 1, xmp, 2, and 3.
c
c      data small(1) / 201354000000000000000b /
c      data small(2) / 000000000000000000000b /
c
c      data large(1) / 577767777777777777777b /
c      data large(2) / 000007777777777777776b /
c
c      data right(1) / 376434000000000000000b /
c      data right(2) / 000000000000000000000b /
c
c      data diver(1) / 376444000000000000000b /
c      data diver(2) / 000000000000000000000b /
c
c      data log10(1) / 377774642023241175717b /
c      data log10(2) / 000007571421742254654b /, sc/987/
c
c     machine constants for the data general eclipse s/200
c
c     small, large, right, diver, log10 should be declared
c     integer small(4), large(4), right(4), diver(4), log10(4)
c
c     note - it may be appropriate to include the following line -
c     static dmach(5)
c
c      data small/20k,3*0/,large/77777k,3*177777k/
c      data right/31420k,3*0/,diver/32020k,3*0/
c      data log10/40423k,42023k,50237k,74776k/, sc/987/
c
c     machine constants for the harris slash 6 and slash 7
c
c      data small(1),small(2) / '20000000, '00000201 /
c      data large(1),large(2) / '37777777, '37777577 /
c      data right(1),right(2) / '20000000, '00000333 /
c      data diver(1),diver(2) / '20000000, '00000334 /
c      data log10(1),log10(2) / '23210115, '10237777 /, sc/987/
c
c     machine constants for the honeywell dps 8/70 series.
c
c      data small(1),small(2) / o402400000000, o000000000000 /
c      data large(1),large(2) / o376777777777, o777777777777 /
c      data right(1),right(2) / o604400000000, o000000000000 /
c      data diver(1),diver(2) / o606400000000, o000000000000 /
c      data log10(1),log10(2) / o776464202324, o117571775714 /, sc/987/
c
c     machine constants for the ibm 360/370 series,
c     the xerox sigma 5/7/9 and the sel systems 85/86.
c
c      data small(1),small(2) / z00100000, z00000000 /
c      data large(1),large(2) / z7fffffff, zffffffff /
c      data right(1),right(2) / z33100000, z00000000 /
c      data diver(1),diver(2) / z34100000, z00000000 /
c      data log10(1),log10(2) / z41134413, z509f79ff /, sc/987/
c
c     machine constants for the interdata 8/32
c     with the unix system fortran 77 compiler.
c
c     for the interdata fortran vii compiler replace
c     the z's specifying hex constants with y's.
c
c      data small(1),small(2) / z'00100000', z'00000000' /
c      data large(1),large(2) / z'7effffff', z'ffffffff' /
c      data right(1),right(2) / z'33100000', z'00000000' /
c      data diver(1),diver(2) / z'34100000', z'00000000' /
c      data log10(1),log10(2) / z'41134413', z'509f79ff' /, sc/987/
c
c     machine constants for the pdp-10 (ka processor).
c
c      data small(1),small(2) / "033400000000, "000000000000 /
c      data large(1),large(2) / "377777777777, "344777777777 /
c      data right(1),right(2) / "113400000000, "000000000000 /
c      data diver(1),diver(2) / "114400000000, "000000000000 /
c      data log10(1),log10(2) / "177464202324, "144117571776 /, sc/987/
c
c     machine constants for the pdp-10 (ki processor).
c
c      data small(1),small(2) / "000400000000, "000000000000 /
c      data large(1),large(2) / "377777777777, "377777777777 /
c      data right(1),right(2) / "103400000000, "000000000000 /
c      data diver(1),diver(2) / "104400000000, "000000000000 /
c      data log10(1),log10(2) / "177464202324, "047674776746 /, sc/987/
c
c     machine constants for pdp-11 fortrans supporting
c     32-bit integers (expressed in integer and octal).
c
c      data small(1),small(2) /    8388608,           0 /
c      data large(1),large(2) / 2147483647,          -1 /
c      data right(1),right(2) /  612368384,           0 /
c      data diver(1),diver(2) /  620756992,           0 /
c      data log10(1),log10(2) / 1067065498, -2063872008 /, sc/987/
c
c      data small(1),small(2) / o00040000000, o00000000000 /
c      data large(1),large(2) / o17777777777, o37777777777 /
c      data right(1),right(2) / o04440000000, o00000000000 /
c      data diver(1),diver(2) / o04500000000, o00000000000 /
c      data log10(1),log10(2) / o07746420232, o20476747770 /, sc/987/
c
c     machine constants for pdp-11 fortrans supporting
c     16-bit integers (expressed in integer and octal).
c
c     small, large, right, diver, log10 should be declared
c     integer small(4), large(4), right(4), diver(4), log10(4)
c
c      data small(1),small(2) /    128,      0 /
c      data small(3),small(4) /      0,      0 /
c
c      data large(1),large(2) /  32767,     -1 /
c      data large(3),large(4) /     -1,     -1 /
c
c      data right(1),right(2) /   9344,      0 /
c      data right(3),right(4) /      0,      0 /
c
c      data diver(1),diver(2) /   9472,      0 /
c      data diver(3),diver(4) /      0,      0 /
c
c      data log10(1),log10(2) /  16282,   8346 /
c      data log10(3),log10(4) / -31493, -12296 /, sc/987/
c
c      data small(1),small(2) / o000200, o000000 /
c      data small(3),small(4) / o000000, o000000 /
c
c      data large(1),large(2) / o077777, o177777 /
c      data large(3),large(4) / o177777, o177777 /
c
c      data right(1),right(2) / o022200, o000000 /
c      data right(3),right(4) / o000000, o000000 /
c
c      data diver(1),diver(2) / o022400, o000000 /
c      data diver(3),diver(4) / o000000, o000000 /
c
c      data log10(1),log10(2) / o037632, o020232 /
c      data log10(3),log10(4) / o102373, o147770 /, sc/987/
c
c     machine constants for the prime 50 series systems
c     with 32-bit integers and 64v mode instructions,
c     supplied by igor bray.
c
c      data small(1),small(2) / :10000000000, :00000100001 /
c      data large(1),large(2) / :17777777777, :37777677775 /
c      data right(1),right(2) / :10000000000, :00000000122 /
c      data diver(1),diver(2) / :10000000000, :00000000123 /
c      data log10(1),log10(2) / :11504046501, :07674600177 /, sc/987/
c
c     machine constants for the sequent balance 8000
c
c      data small(1),small(2) / sh,   /
c      data large(1),large(2) / ,  fefffff /
c      data right(1),right(2) / sh,  ca00000 /
c      data diver(1),diver(2) / sh,  cb00000 /
c      data log10(1),log10(2) / f79ff,  fd34413 /, sc/987/
c
c     machine constants for the univac 1100 series.
c
c      data small(1),small(2) / o000040000000, o000000000000 /
c      data large(1),large(2) / o377777777777, o777777777777 /
c      data right(1),right(2) / o170540000000, o000000000000 /
c      data diver(1),diver(2) / o170640000000, o000000000000 /
c      data log10(1),log10(2) / o177746420232, o411757177572 /, sc/987/
c
c     machine constants for the vax unix f77 compiler
c
c      data small(1),small(2) /        128,           0 /
c      data large(1),large(2) /     -32769,          -1 /
c      data right(1),right(2) /       9344,           0 /
c      data diver(1),diver(2) /       9472,           0 /
c      data log10(1),log10(2) /  546979738,  -805796613 /, sc/987/
c
c     machine constants for the vax-11 with
c     fortran iv-plus compiler
c
c      data small(1),small(2) / z00000080, z00000000 /
c      data large(1),large(2) / zffff7fff, zffffffff /
c      data right(1),right(2) / z00002480, z00000000 /
c      data diver(1),diver(2) / z00002500, z00000000 /
c      data log10(1),log10(2) / z209a3f9a, zcff884fb /, sc/987/
c
c     machine constants for vax/vms version 2.2
c
c      data small(1),small(2) /       '80'x,        '0'x /
c      data large(1),large(2) / 'ffff7fff'x, 'ffffffff'x /
c      data right(1),right(2) /     '2480'x,        '0'x /
c      data diver(1),diver(2) /     '2500'x,        '0'x /
c      data log10(1),log10(2) / '209a3f9a'x, 'cff884fb'x /, sc/987/
c
c  ***  issue stop 779 if all data statements are commented...
      if (sc .ne. 987) stop 779
c  ***  issue stop 778 if all data statements are obviously wrong...
      if (dmach(4) .ge. 1.0d0) stop 778
      if (i .lt. 1  .or.  i .gt. 5) goto 999
      d1mach = dmach(i)
      return
  999 write(*,1999) i
 1999 format(' d1mach - i out of bounds',i10)
      stop
      end

