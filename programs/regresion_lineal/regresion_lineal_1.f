c
c     #########################################################################
c     ##                                                                     ##
c     ##   Programa egresion_lineal                                          ##
c     ##   input: channel 30 use to import data                              ##
c     ##   output: channel 40                                                ##
c     ##                                                                     ##
c     ##                                                                     ##
c     #########################################################################
c

      program regresion_lineal

      implicit double precision (a-h,o-z)
      implicit integer (i-n)

      character*100 filehb

      read(5,*)filehb

      
      open(30,file=filehb)

      sumx=0.d0
      sumy=0.d0
      sumx2=0.d0
      sumy2=0.d0
      sumxy=0.d0

      icon=0
      
 10   continue

      read(30,*,end=20)x,y

      icon=icon+1

      sumx=sumx+x
      sumy=sumy+y
      sumx2=sumx2+x**2
      sumy2=sumy2+y**2
      sumxy=sumxy+x*y
      
      goto 10

 20   continue

      if(icon.eq.1)goto 40

      sumx=sumx/dfloat(icon)
      sumy=sumy/dfloat(icon)
      sumx2=sumx2/dfloat(icon)
      sumy2=sumy2/dfloat(icon)
      sumxy=sumxy/dfloat(icon)

      a=(sumxy-sumx*sumy)/(sumx2-sumx**2)
      b=sumy-a*sumx
      r=(sumxy-sumx*sumy)/dsqrt((sumx2-sumx**2)*(sumy2-sumy**2))
      
      write(40,1000)a,b,r,icon
      
 1000 format(1x,'p=',1x,g13.5,1x,'o=',1x,g13.5,1x,'r=',1x,g13.5,
     $     x,'N=',1x,i8)
c      write(44,1030)textiop
c 1030 format(1x,a100)

 40         continue
            
      end
