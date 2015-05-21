      program stubs
c gfortran -arch x86_64 -g -O0 -save-temps -o stubs stubs.F

      i = itestme()
      print *,"itestme first returned", i
      j = itestme()
      print *,"itestme next returned", j
c      print *,"itestme last returned",itestme()
      i = itest2()
      print *,"itest2 returned", i
      end

      integer function itestme()
      print *, "Hello, world"
      itestme = 1
      return
      end

      integer function itest2()
      itest2 = itestme()
      print *, "Hello, Earth"
      return
      end
