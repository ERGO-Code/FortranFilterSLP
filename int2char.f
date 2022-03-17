      character*4 function int2char ( int )

c     ... converts integer int into character string

      implicit none

      integer     int, int1, d(4), i
      character*4 char

c     ==========================   function body   =====================

c     ... find thousands, hundreds, tens and ones
      int1 = int
      d(1) = int1 / 1000
      int1 = int1 - d(1)*1000
      d(2) = int1 / 100
      int1 = int1 - d(2)*100
      d(3) = int1 / 10
      d(4) = int1 - d(3)*10

c     ... convert each into a character
      char = ' '
      do i=1,4
         if (d(i).eq.0) then
            char(i:i) = '0'
         elseif (d(i).eq.1) then
            char(i:i) = '1'
         elseif (d(i).eq.2) then
            char(i:i) = '2'
         elseif (d(i).eq.3) then
            char(i:i) = '3'
         elseif (d(i).eq.4) then
            char(i:i) = '4'
         elseif (d(i).eq.5) then
            char(i:i) = '5'
         elseif (d(i).eq.6) then
            char(i:i) = '6'
         elseif (d(i).eq.7) then
            char(i:i) = '7'
         elseif (d(i).eq.8) then
            char(i:i) = '8'
         elseif (d(i).eq.9) then
            char(i:i) = '9'
         endif
      enddo

c     ... save result in int2char
      int2char = char

      return
      end


