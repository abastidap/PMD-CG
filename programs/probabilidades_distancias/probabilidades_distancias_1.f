c
c     #########################################################################
c     ##                                                                     ##
c     ##   Program probabilidades_distancias                                 ##
c     ##                                                                     ##
c     #########################################################################
c
c     Programa generar distancias aleatorias de los atomos del backbone a partir de distribuciones

     
      program probabilidades_distancias

      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)

      character*100 fichaux

      real*8,allocatable :: x(:),y(:)
      real*8,allocatable :: ytest(:)
      
      read(5,*)nfilas
      read(5,*)fichaux
      read(5,*)ndatos

      
      allocate(x(nfilas),y(nfilas))
      allocate(ytest(nfilas))


c     Lectura de histograma

      open(31,file=fichaux)
      
      do i=1,nfilas
         read(31,*)x(i),y(i)
c         print*,i,x(i),y(i)
      enddo


c     Normalizacion

      sum=0.d0
      do i=1,nfilas
         sum=sum+y(i)
      enddo
      
      do i=1,nfilas
         y(i)=y(i)/sum
      enddo

c     inicializacion del test

      do i=1,nfilas
         ytest(i)=0.d0
      enddo

      
c     generación de numeros aleatorios


      do m=1,ndatos
           
         call random_number(alea)

         sum=0.d0
         do i=1,nfilas
            sum=sum+y(i)
            if(sum.ge.alea)goto 100
         enddo
         
 100     continue
         
         xalea=x(i)
            
         ytest(i)=ytest(i)+1.d0
         
         write(40,200)xalea
 200     format(2(1x,f10.3))
         
      enddo

      close(40)
      
c     test

      sum=0.d0
      do i=1,nfilas
         sum=sum+ytest(i)
      enddo

      do i=1,nfilas
         ytest(i)=ytest(i)/sum
      enddo
  
      sx=0.d0
      sx2=0.d0
      sxy=0.d0
      sy=0.d0
      sy2=0.d0

      m=0
      do i=1,nfilas
         if(y(i).gt.0.d0.or.ytest(i).gt.0.d0)then
            m=m+1
            sx=sx+y(i)
            sx2=sx2+y(i)**2
            sy=sy+ytest(i)
            sy2=sy2+ytest(i)**2
            sxy=sxy+y(i)*ytest(i)
         endif
      enddo

      a=dfloat(m)
      r=(a*sxy-sx*sy)/dsqrt(a*sx2-sx**2)
     $     /dsqrt(a*sy2-sy**2)
      
      write(41,*)' Numero de valores = ',ndatos
      write(41,*)' r = ',r
         
      
      end      
