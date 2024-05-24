c
c     #########################################################################
c     ##                                                                     ##
c     ##   Program probabilidades_ramachandran                               ##
c     ##                                                                     ##
c     #########################################################################
c
c     Programa generar parejas de diedros probabilisticamente a partir de mapas de ramanchandran

c     v2 -> incluimos las probabilidades calculadas con p(2)
      
      program probabilidades_ramachandran

      implicit double precision (a-h,o-z)
      implicit integer*8 (i-n)

      character*100 fichaux,seq,taux
      character*1,allocatable :: aa(:)

      real*8,allocatable :: rama(:,:,:)
      real*8,allocatable :: test(:,:,:)
      
c     Res, reg. Centros y semianchos de las regiones de cada residuo
      real*8, allocatable :: phic(:,:),psic(:,:)
      real*8, allocatable :: dephi(:,:),depsi(:,:)

c     Res. Numero de regiones de un residuo
      integer, allocatable :: nreg(:)

c     Res, reg. Nombre de cada region de un residuo
      character*10, allocatable :: areg(:,:)

      real*8, allocatable :: ramaizq(:,:,:,:),ramader(:,:,:,:)

      real*8, allocatable :: phiaux(:),psiaux(:)
      
      read(5,*)seq
      read(5,*)n2d
      read(5,*)ndatos

      ndi=len(trim(seq))
      
      allocate(rama(ndi,n2d,n2d))
      allocate(test(ndi,n2d,n2d))

      allocate(aa(ndi))

      do i=1,ndi
         aa(i)=seq(i:i)
      enddo

c     Consideramos un maximo de 5 regiones

      nregmax=5
      
      allocate (nreg(ndi))
      allocate (areg(ndi,nregmax))
      allocate (phic(ndi,nregmax),psic(ndi,nregmax))
      allocate (dephi(ndi,nregmax),depsi(ndi,nregmax)) 
      allocate (ramaizq(ndi,nregmax,n2d,n2d))
      allocate (ramader(ndi,nregmax,n2d,n2d))
      allocate (phiaux(ndi),psiaux(ndi))

c     Lectura de las regiones de cada residuo 
      
      do i=1,ndi
         read(5,*)nreg(i)
         nreg(i)=nreg(i)+1
         areg(i,nreg(i))='Resto'
         
         do j=1,nreg(i)-1
            read(5,*)areg(i,j)
            read(5,*)phic(i,j),psic(i,j),dephi(i,j),depsi(i,j)                    
         enddo               
      enddo
      

c     Calculo con p(1) --------------------------------------------
      
c     lectura y renormalizacion a 1 de los ficheros

      taux='raman-'//trim(aa(1))//trim(aa(2))//trim(aa(3))//'-di1.dat'
      open(31,file=taux)

      do i=1,n2d
         read(31,*)(rama(1,i,j),j=1,n2d)
      enddo

      close(31)

      do n=2,ndi-1
         taux='raman-'//trim(aa(n-1))//trim(aa(n))//trim(aa(n+1))
     $        //'-di2.dat'
         open(31,file=taux)
         do i=1,n2d
            read(31,*)(rama(n,i,j),j=1,n2d)
         enddo
         close(31)
      enddo

      taux='raman-'//trim(aa(ndi-2))//trim(aa(ndi-1))
     $     //trim(aa(ndi))//'-di3.dat'
      open(31,file=taux)

      do i=1,n2d
         read(31,*)(rama(ndi,i,j),j=1,n2d)
      enddo

      close(31)

      
      do n=1,ndi
         
         sum=0.d0
         do i=1,n2d
            do j=1,n2d
               sum=sum+rama(n,i,j)
            enddo
         enddo

         do i=1,n2d
            do j=1,n2d
               rama(n,i,j)=rama(n,i,j)/sum
            enddo
         enddo   
         
      enddo

c     inicializacion del test

      do n=1,ndi
         do i=1,n2d
            do j=1,n2d
               test(n,i,j)=0.d0
            enddo
         enddo
      enddo

c     sorteo

      da=360.d0/dfloat(n2d)

      icon=0
      
      do m=1,ndatos
         do n=1,ndi

            icon=icon+1
            
            call random_number(alea)

            sum=0.d0
            do i=1,n2d
               do j=1,n2d
c$$$                  sum=sum+rama(n,i,j)
                  sum=sum+rama(n,j,i)
                  if(sum.ge.alea)goto 100
               enddo
            enddo

 100        continue
            
            phi=-180.d0+dfloat(i-1)*da
            psi=-180.d0+dfloat(j-1)*da
            
c$$$            ix=int(sngl((phi+180.d0)/da))+1
c$$$            iy=int(sngl((psi+180.d0)/da))+1         
c$$$            test(n,ix,iy)=test(n,ix,iy)+1.d0

            test(n,j,i)=test(n,j,i)+1.d0

c     sumamos un valor entre 0 y la resolucion de los mapas

            call random_number(alea)
            phi=phi+da*alea
            call random_number(alea)
            psi=psi+da*alea      
            
            write(50,200)phi,psi
 200        format(2(1x,f9.2))
            
         enddo        
      enddo

      close(50)
      
c     test

      do n=1,ndi

         do i=1,n2d
            do j=1,n2d
               test(n,i,j)=test(n,i,j)/dfloat(ndatos)
            enddo
         enddo

         sx=0.d0
         sx2=0.d0
         sxy=0.d0
         sy=0.d0
         sy2=0.d0

         m=0
         do i=1,n2d
            do j=1,n2d
               if(rama(n,i,j).gt.0.d0.or.test(n,i,j).gt.0.d0)then
                  m=m+1
                  sx=sx+rama(n,i,j)
                  sx2=sx2+rama(n,i,j)**2
                  sy=sy+test(n,i,j)
                  sy2=sy2+test(n,i,j)**2
                  sxy=sxy+rama(n,i,j)*test(n,i,j)
               endif
            enddo
         enddo

         a=dfloat(m)
         r=(a*sxy-sx*sy)/dsqrt(a*sx2-sx**2)
     $        /dsqrt(a*sy2-sy**2)

         write(40,*)' diedro = ',n
         write(40,*)' Numero de valores = ',ndatos
         write(40,*)' r = ',r
         
      enddo

c     Calculo con p(2) --------------------------------------------
      
c     lectura y renormalizacion a 1 de los ficheros

      do k=1,nreg(2)
         
         taux='raman-'//trim(aa(1))//trim(aa(2))//trim(aa(3))//
     $        '-di1-d-'//trim(areg(2,k))//'.dat'
         open(31,file=taux)

         do i=1,n2d
            read(31,*)(ramader(1,k,i,j),j=1,n2d)
         enddo

         close(31)
      enddo

      do n=2,ndi-1
         do k=1,nreg(n+1)
            taux='raman-'//trim(aa(n-1))//trim(aa(n))//trim(aa(n+1))
     $           //'-di2-d-'//trim(areg(n+1,k))//'.dat'
            open(31,file=taux)
            do i=1,n2d
               read(31,*)(ramader(n,k,i,j),j=1,n2d)
            enddo
            close(31)
         enddo

         do k=1,nreg(n-1)
            taux='raman-'//trim(aa(n-1))//trim(aa(n))//trim(aa(n+1))
     $           //'-di2-i-'//trim(areg(n-1,k))//'.dat'
            open(31,file=taux)
            do i=1,n2d
               read(31,*)(ramaizq(n,k,i,j),j=1,n2d)
            enddo
            close(31)
         enddo
         
      enddo

      do k=1,nreg(ndi-1)
         taux='raman-'//trim(aa(ndi-2))//trim(aa(ndi-1))
     $        //trim(aa(ndi))//'-di3-i-'//trim(areg(ndi-1,k))//'.dat'
         open(31,file=taux)

         do i=1,n2d
            read(31,*)(ramaizq(ndi,k,i,j),j=1,n2d)
         enddo

         close(31)
      enddo
      
      do n=1,ndi-1
         do k=1,nreg(n+1)
            
            sum=0.d0
            do i=1,n2d
               do j=1,n2d
                  sum=sum+ramader(n,k,i,j)
               enddo
            enddo

            do i=1,n2d
               do j=1,n2d
                  ramader(n,k,i,j)=ramader(n,k,i,j)/sum
               enddo
            enddo          
         enddo
      enddo

      do n=2,ndi
         do k=1,nreg(n-1)
            
            sum=0.d0
            do i=1,n2d
               do j=1,n2d
                  sum=sum+ramaizq(n,k,i,j)
               enddo
            enddo

            do i=1,n2d
               do j=1,n2d
                  ramaizq(n,k,i,j)=ramaizq(n,k,i,j)/sum
               enddo
            enddo          
         enddo
      enddo

c     sorteo condicionado a la izq
      
      do nn=1,ndatos

         call random_number(alea)

         sum=0.d0
         do i=1,n2d
            do j=1,n2d
c$$$               sum=sum+rama(1,i,j)
               sum=sum+rama(1,j,i)
               if(sum.ge.alea)goto 300
            enddo
         enddo

 300     continue
         
         phi=-180.d0+dfloat(i-1)*da
         psi=-180.d0+dfloat(j-1)*da
         
         call random_number(alea)
         phi=phi+da*alea
         call random_number(alea)
         psi=psi+da*alea      
         
         write(51,200)phi,psi
         
         do m=1,nreg(1)-1
            
            phi1d=phi-phic(1,m)
            psi1d=psi-psic(1,m)
            
            if(phi1d.gt.180.d0)phi1d=phi1d-360.d0
            if(psi1d.gt.180.d0)psi1d=psi1d-360.d0
            
            if(phi1d.lt.-180.d0)phi1d=phi1d+360.d0
            if(psi1d.lt.-180.d0)psi1d=psi1d+360.d0
            
            difphi=dabs(phi1d)
            difpsi=dabs(psi1d)
            
            if(difphi.le.dephi(1,m).and.difpsi.le.depsi(1,m))then
               ireg=m
               goto 400
            endif
            
         enddo

         ireg=nreg(1)
         
 400     continue

         do n=2,ndi

            call random_number(alea)

            sum=0.d0
            do i=1,n2d
               do j=1,n2d
c$$$                  sum=sum+ramaizq(n,ireg,i,j)
                  sum=sum+ramaizq(n,ireg,j,i)
                  if(sum.ge.alea)goto 500
               enddo
            enddo

 500        continue
            
            phi=-180.d0+dfloat(i-1)*da
            psi=-180.d0+dfloat(j-1)*da
            
            call random_number(alea)
            phi=phi+da*alea
            call random_number(alea)
            psi=psi+da*alea      
            
            write(51,200)phi,psi
            
            do m=1,nreg(n)-1
               
               phi1d=phi-phic(n,m)
               psi1d=psi-psic(n,m)
               
               if(phi1d.gt.180.d0)phi1d=phi1d-360.d0
               if(psi1d.gt.180.d0)psi1d=psi1d-360.d0
               
               if(phi1d.lt.-180.d0)phi1d=phi1d+360.d0
               if(psi1d.lt.-180.d0)psi1d=psi1d+360.d0
               
               difphi=dabs(phi1d)
               difpsi=dabs(psi1d)
               
               if(difphi.le.dephi(n,m).and.difpsi.le.depsi(n,m))then
                  ireg=m
                  goto 600
               endif
               
            enddo

            ireg=nreg(n)
            
 600        continue

         enddo


      enddo

      close(51)

      
c     sorteo condicionado a la der
      
      do nn=1,ndatos

         call random_number(alea)

         sum=0.d0
         do i=1,n2d
            do j=1,n2d
c$$$               sum=sum+rama(ndi,i,j)
               sum=sum+rama(ndi,j,i)
               if(sum.ge.alea)goto 310
            enddo
         enddo

 310     continue
         
         phi=-180.d0+dfloat(i-1)*da
         psi=-180.d0+dfloat(j-1)*da
         
         call random_number(alea)
         phi=phi+da*alea
         call random_number(alea)
         psi=psi+da*alea      

         phiaux(ndi)=phi
         psiaux(ndi)=psi
         
         do m=1,nreg(ndi)-1
            
            phi1d=phi-phic(ndi,m)
            psi1d=psi-psic(ndi,m)
            
            if(phi1d.gt.180.d0)phi1d=phi1d-360.d0
            if(psi1d.gt.180.d0)psi1d=psi1d-360.d0
            
            if(phi1d.lt.-180.d0)phi1d=phi1d+360.d0
            if(psi1d.lt.-180.d0)psi1d=psi1d+360.d0
            
            difphi=dabs(phi1d)
            difpsi=dabs(psi1d)
            
            if(difphi.le.dephi(ndi,m).and.difpsi.le.depsi(ndi,m))then
               ireg=m
               goto 410
            endif
            
         enddo

         ireg=nreg(ndi)
         
 410     continue

         do n=ndi-1,1,-1

            call random_number(alea)

            sum=0.d0
            do i=1,n2d
               do j=1,n2d
c$$$                  sum=sum+ramader(n,ireg,i,j)
                  sum=sum+ramader(n,ireg,j,i)
                  if(sum.ge.alea)goto 510
               enddo
            enddo

 510        continue
            
            phi=-180.d0+dfloat(i-1)*da
            psi=-180.d0+dfloat(j-1)*da
            
            call random_number(alea)
            phi=phi+da*alea
            call random_number(alea)
            psi=psi+da*alea      
            
            phiaux(n)=phi
            psiaux(n)=psi
            
            do m=1,nreg(n)-1
               
               phi1d=phi-phic(n,m)
               psi1d=psi-psic(n,m)
               
               if(phi1d.gt.180.d0)phi1d=phi1d-360.d0
               if(psi1d.gt.180.d0)psi1d=psi1d-360.d0
               
               if(phi1d.lt.-180.d0)phi1d=phi1d+360.d0
               if(psi1d.lt.-180.d0)psi1d=psi1d+360.d0
               
               difphi=dabs(phi1d)
               difpsi=dabs(psi1d)
               
               if(difphi.le.dephi(n,m).and.difpsi.le.depsi(n,m))then
                  ireg=m
                  goto 610
               endif
               
            enddo

            ireg=nreg(n)
            
 610        continue

         enddo

         do n=1,ndi
            write(52,200)phiaux(n),psiaux(n)
         enddo

         

      enddo
      
      close(52)
      
      end      
