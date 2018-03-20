      program kicky
C    
C     I am sorry that this is in fortran. 
C    This code reads in data about the globular clusters, then randomly 
C selects from the distribution of BH spins, mass ratios, ecc and 
C  runs through a 
C  chain of BH 25 mergers It calculates the recoil at each 
C   merger and compares it
C  to the escape velocity of a globular. You will likely only need lines 
C that have to do with the recoil equations. 
C They are scattered throughout -- lines 158-179,207-218,285-294,384-398,
C 415-441. At some point, I'll port the math to a python module. Promise.


      parameter(n1=501,n2=1001,n3=1000000)
      integer i,j,k,count,iseed,mbin,kkk,histo(2001),countesc
cc      real*8 ecc,q,probq,probecc,
      real*8 ecc1(n2+1),q,probm,probecc1(n2+1),allprob,eccmc,
     &     bigmass(n3),
     &     totprob,theta,testprob,minit,cprobeccbin,
     &     cprobecc(n2+1),
     &     minmass1, minmass2,maxmass,maxmass2,
     &     mm,ne,massran,mtest(n1),cprobm(n1)

      real*8 ecc2(n2+1),probecc2(n2+1),ecc3(n2+1),probecc3(n2+1),
     &       ecc4(n2+1),probecc4(n2+1),ecc5(n2+1),probecc5(n2+1),
     &       ecc6(n2+1),probecc6(n2+1),ecc7(n2+1),probecc7(n2+1),
     &       ecc8(n2+1),probecc8(n2+1),ecc9(n2+1),probecc9(n2+1),
     &       ecc(n2+1),probecc(n2+1),vescgc(139),cprobecc1(n2+1),
     &       cprobecc2(n2+1),cprobecc3(n2+1),cprobecc4(n2+1),
     &       cprobecc5(n2+1),cprobecc6(n2+1),cprobecc7(n2+1),
     &       cprobecc8(n2+1),cprobecc9(n2+1)

      real*8 A, H, Kcostheta, vtot(n3),B,vm,vperp,vpara,spin1para,
     &      spin2para,spin2perp,spin1perp,psi,pi,m1,m2,s1para,
     &      cos1angle,
     &      s2para(n3),cos2angle,bound,normal,vtotmax,probcollapse,
     &      smallmass,smallmasstemp,cumulative,eccprox(n3),mhisto

      real*8 Z1,Z2,Z3,uphi,jay,r,grade

      character*100 string

      rn3=n3

            do i=1,n2
               cprobecc1(i)=0
               cprobecc2(i)=0
               cprobecc3(i)=0
               cprobecc4(i)=0
               cprobecc5(i)=0
               cprobecc6(i)=0
               cprobecc7(i)=0
               cprobecc8(i)=0
               cprobecc9(i)=0
            enddo


      open(11,file='p0010.dat')
      cprobecc1(1)=0
      ecc1(1)=0
      do i=2,n2+1
         read(11,*) ecc1(i),probecc1(i)
         ecc1(i)=ecc1(i)+0.001
         cprobecc1(i)=cprobecc1(i-1)+probecc1(i)
      enddo
      close(11)

      open(12,file='p0020.dat')
      cprobecc2(1)=0
      ecc2(1)=0
      do i=2,n2+1
         read(12,*) ecc2(i),probecc2(i)
         ecc2(i)=ecc2(i)+0.001
         cprobecc2(i)=cprobecc2(i-1)+probecc2(i)
      enddo
      close(12)

      open(13,file='p0030.dat')
      cprobecc3(1)=0
      ecc3(1)=0

      do i=2,n2+1
         read(13,*) ecc3(i),probecc3(i)
         ecc3(i)=ecc3(i)+0.001
         cprobecc3(i)=cprobecc3(i-1)+probecc3(i)
      enddo
      close(13)

      open(14,file='p0050.dat')
      cprobecc4(1)=0
      ecc4(1)=0
      do i=2,n2+1
         read(14,*) ecc4(i),probecc4(i)
         ecc4(i)=ecc4(i)+0.001
         cprobecc4(i)=cprobecc4(i-1)+probecc4(i)
      enddo
      close(14)

      open(15,file='p0100.dat')
      cprobecc5(1)=0
      ecc5(1)=0
      do i=2,n2+1
         read(15,*) ecc5(i),probecc5(i)
         ecc5(i)=ecc5(i)+0.001
         cprobecc5(i)=cprobecc5(i-1)+probecc5(i)
      enddo
      close(15)

      open(16,file='p0200.dat')
      cprobecc6(1)=0
      ecc6(1)=0
      do i=2,n2+1
         read(16,*) ecc6(i),probecc6(i)
         ecc6(i)=ecc6(i)+0.001
         cprobecc6(i)=cprobecc6(i-1)+probecc6(i)
      enddo
      close(16)

      open(17,file='p0300.dat')
      cprobecc7(1)=0
      ecc7(1)=0
      do i=2,n2+1
         read(17,*) ecc7(i),probecc7(i)
         ecc7(i)=ecc7(i)+0.001
         cprobecc7(i)=cprobecc7(i-1)+probecc7(i)
      enddo
      close(17)

      open(22,file='p0500.dat')
      cprobecc8(1)=0
      ecc8(1)=0
      do i=2,n2+1
         read(22,*) ecc8(i),probecc8(i)
         ecc8(i)=ecc8(i)+0.001
         cprobecc8(i)=cprobecc8(i-1)+probecc8(i)
      enddo
      close(22)

      open(23,file='p1000.dat')
      cprobecc9(1)=0
      ecc9(1)=0
      do i=2,n2+1
         read(23,*) ecc9(i),probecc9(i)
         ecc9(i)=ecc9(i)+0.001
         cprobecc9(i)=cprobecc9(i-1)+probecc9(i)
      enddo
      close(23)





      open(11,file='escape.dat')
      do i=1,139
         read(11,*) rblah, vescgc(i)
         vescgc(i)=vescgc(i)*sqrt(2.0)
      enddo
      close(11)
      
      open(99,file='monte.mwgc.dat')
      open(9,file='monte.escape.mwgc.dat')
      open(10,file='vhisto.mwgc.kroupa.dat')

      A =  12000.0d0
      H = 7.3d3
      Kcostheta=6.0d4
      B=-0.93
      psi=45.0d0
      pi=3.1415926535
      iseed=-53562
      minmass1=0.0
      minmass2=20.0d0
      maxmass=80.0d0
      maxmass2=120.0



         do kkk=1,100
            rkkk=kkk
            if(kkk.eq.1)then
  
             minit=50.0
            else
               minit=rkkk*100.
            endif


            do i=1,n3
               eccprox(i)=0
            enddo



            do ll=1,139

            ne=1.0
            smallmass=10.
            allprob=1.0


            do j=1,25
            vtotmax=0
c            countesc=0.0
c               write(*,*)'j, minit, vesc', j,minit,vescgc(ll)
               bound=0
               if(j.eq.1)then
                  do k=1,n3
                     bigmass(k)=minit
                  enddo
               endif
               do k=1,n3
                  psi=360.0*ran2(iseed)

c         theta=2.*pi*ran2(iseed)
                  theta=360.0*ran2(iseed)
                  cos2angle=2.*ran2(iseed)-1.0
                  if(j.eq.1)then
                     s2para(k)=ran2(iseed)
c                    s2para(k)=2.*ran2(iseed)-1.0
c                    s2para(k)=0.001
c                    s2para(k)=0.99

                     mbin=0
                     do massran=minmass1,maxmass2,1.0
                        mbin=mbin+1.0
                        cprobm(mbin)=0
                     enddo
                     mbin=0
                     cumulative=0
                     do massran=minmass1,maxmass2,1.0
                        mbin=mbin+1.0
                        mtest(mbin)=massran

c               if(massran.le.0.08)then
c                  cprobm(mbin)=ne*0.254d0*(massran**(-0.3))
c               elseif((massran.gt.0.08).and.(massran.lt.0.5))then
c                  cprobm(mbin)=ne*0.254d0*(massran**(-1.3))
c               elseif((massran.ge.0.5).and.(massran.le.maxmass))then
c                  cprobm(mbin)=ne*0.254d0*(massran**(-2.3))
c               else
c                  cprobm(mbin)=0.0d0
c               endif


                        if((massran.gt.10).and.
     &                       (massran.lt.20))then
                           cprobm(mbin)=(1e-2)*(20.-massran)/20.0 
                        elseif((massran.ge.20).and.
     &                       (massran.le.80))then
                           cprobm(mbin)=5e-5
                        elseif((massran.ge.80).and.
     &                       (massran.le.120))then
                           cprobm(mbin)=(8e-5)*
     &                         (120.-massran)/360.0 + 1e-6
                        elseif((massran.ge.0).and.
     &                          (massran.le.10))then
                           cprobm(mbin)=(5e-3)*(massran)/10.0  
                        else

                        endif

                        cumulative=cumulative+cprobm(mbin)
c               write(*,*) mbin,cprobm(mbin),cumulative

                        cprobm(mbin)=cumulative
                     enddo
                     do i=1,mbin
                        cprobm(i)=cprobm(i)/cumulative
                     enddo

                  endif


                  testm=0
                  countm=1
                  do while(testm.eq.0)
                     cprobmbin=ran2(iseed)
                     do i=2,mbin
                        if((cprobmbin.gt.cprobm(i-1)).and.
     &                       (cprobmbin.le.cprobm(i)))then
                           smallmass=mtest(i-1)
                           if(smallmass.ge.10.0)then
                              testm=1.0
                           endif
                        endif
                     enddo
                     countm=countm+1
                  enddo

                  cos1angle=2.*ran2(iseed)-1.0
                  spin2para=(s2para(k))*cos2angle
                  spin2perp=(s2para(k))*sqrt(1.-cos2angle*cos2angle)
                  s1para=ran2(iseed)
c         s1para=2.*ran2(iseed)-1.0
c         s1para=0.001
                  spin1para=(s1para)*cos1angle
                  spin1perp=(s1para)*sqrt(1.-cos1angle*cos1angle)
                  q=smallmass/bigmass(k)

                  if(q.gt.1)then
                     q=bigmass(k)/smallmass
                     smallmasstemp=smallmass
                     smallmass=bigmass(k)
                     bigmass(k)=smallmasstemp
                  endif

c         write(*,*) k,smallmass
                  if((q.le.1.).and.(q.gt.0.5))then
                     do i=1,n2+1
                        cprobecc(i)=cprobecc1(i)
                        probecc(i)=probecc1(i)
                        ecc(i)=ecc1(i)
                     enddo
                  elseif((q.le.0.5).and.(q.gt.0.33))then
                     do i=1,n2+1
                        cprobecc(i)=cprobecc2(i)
                        probecc(i)=probecc2(i)
                        ecc(i)=ecc2(i)
                     enddo
                  elseif((q.le.0.33).and.(q.gt.0.2))then
                     do i=1,n2+1
                        cprobecc(i)=cprobecc3(i)
                        probecc(i)=probecc3(i)
                        ecc(i)=ecc3(i)
                     enddo
                  elseif((q.le.0.2).and.(q.gt.0.1))then
                     do i=1,n2+1
                        cprobecc(i)=cprobecc4(i)
                        probecc(i)=probecc4(i)
                        ecc(i)=ecc4(i)
                     enddo
                  elseif((q.le.0.1).and.(q.gt.0.05))then
                     do i=1,n2+1
                        cprobecc(i)=cprobecc5(i)
                        probecc(i)=probecc5(i)
                        ecc(i)=ecc5(i)
                     enddo
                  elseif((q.le.0.05).and.(q.gt.0.033))then
                     do i=1,n2+1
                        cprobecc(i)=cprobecc6(i)
                        probecc(i)=probecc6(i)
                        ecc(i)=ecc6(i)
                     enddo
                  elseif((q.le.0.033).and.(q.gt.0.02))then
                     do i=1,n2+1
                        cprobecc(i)=cprobecc7(i)
                        probecc(i)=probecc7(i)
                        ecc(i)=ecc7(i)
                     enddo
                  elseif((q.le.0.02).and.(q.ge.0.01))then
                     do i=1,n2+1                     
                        cprobecc(i)=cprobecc8(i)
                        probecc(i)=probecc8(i)
                        ecc(i)=ecc8(i)
                     enddo
                  elseif((q.le.0.01).and.(q.ge.0.0))then
                     do i=1,n2+1                     
                        cprobecc(i)=cprobecc9(i)
                        probecc(i)=probecc9(i)
                        ecc(i)=ecc9(i)
                     enddo
                  else
                   
                  endif
c                  do i=1,n2+1

c                  enddo

                  testprob=0
                  count=0
                  do while(testprob.eq.0)
                     cprobeccbin=ran2(iseed)
                     do i=1,n2
                        if((cprobeccbin.gt.cprobecc(i)).and.
     &                       (cprobeccbin.le.cprobecc(i+1)))then
c                  write(*,*) i,ecc(i),cprobecc(i)                        
                           eccmc=ecc(i+1)
                           testprob=probecc(i+1)
c     write(*,*) eccmc

                        endif
                     enddo
                     count=count+1
                  
                  enddo
c                     write(*,*) 'count',k, count,eccmc,testprob,
c     &                        cprobeccbin


                  vm=A*(((q*q*(1.-q))/(1.+q)**5.)*
     &                 (1.+B*(q/(1.+q)**2.0d0)))

                  vpara=Kcostheta*cos(pi*theta/180.)*
     &                 ((q*q/(1.+q)**5.))*
     &                 (spin2perp-q*spin1perp)
                  vperp=H*((q*q/(1.+q)**5.))*(spin2para-q*spin1para)
                  vtot(k)=sqrt((vm+
     &                 vperp*cos(psi*pi/180.0))**2.0+
     &                 (vperp*sin(psi*pi/180.0))**2.0 +
     &                 vpara*vpara)


                  vtot(k)=vtot(k)*(1.0d0+eccmc)
                  eccprox(k)=eccmc
                  if(abs(vtot(k)).gt.vtotmax)then
                     vtotmax=abs(vtot(k))
                  endif
                  if(abs(vtot(k)).le.vescgc(ll))then
                     bound=bound+1.0
                  endif

c         write(9,999) k,bigmass(k),smallmass,vtot(k),eccprox(k),
c     &              s2para(k),s1para,theta, cos2angle,cos1angle,psi
         
c                  if(vtot(k).ge.100)then
c         write(99,999) k,bigmass(k),smallmass,vtot(k),eccprox(k),
c     &              s2para(k),s1para,theta, cos2angle,cos1angle,psi
c                  endif

                  s2para(k)=abs(s2para(k))
                  jay=s2para(k)*bigmass(k)*bigmass(k)
                  mu=smallmass*bigmass(k)/(smallmass+bigmass(k))
c     if(k.eq.1)then
                  bigmass(k)=bigmass(k)+0.95*mu
c      endif

      
                  grade=s2para(k)/sqrt(s2para(k)*s2para(k))
                  Z1=(1.0+s2para(k))**(1.0/3.0)+
     &                 (1.0-s2para(k))**(1.0/3.0)
                  Z1=1.0+((1.0-s2para(k)*s2para(k))**(1.0/3.0))
     &                 *Z1
                  Z2=sqrt(3.0*s2para(k)*s2para(k)+Z1*Z1)
                  r=3.0+Z2-grade*sqrt((3.0-Z1)*(3.0+Z1+2.0*Z2))
                  uphi=sqrt(r)*(r*r-2.0*s2para(k)*
     &                 sqrt(r)+s2para(k)*s2para(k))
                  uphi=uphi/(r*sqrt(r*r-3.0*r+2.0*
     &                 s2para(k)*sqrt(r)))
                  uphi=uphi*bigmass(k)
                  jay=jay+0.75*(mu*uphi*cos2angle)
                  jay=sqrt(jay*jay+0.75*(mu*mu*uphi*uphi*
     &                 (1.0-cos2angle*cos2angle)))
                  s2para(k)=(jay/(bigmass(k)*bigmass(k)))
                  if(abs(s2para(k)).ge.1.0)then
                     s2para(k)=0.998
                  endif
c     if(s2para(k).lt.0)then
c     s2para(k)=0.001
c     endif


c     write(*,*) k,mu,bigmass(k),grade,Z1,Z2,r,uphi,jay,s2para(k),vtot
c                  if(mod(k,500).eq.0)then
c 

    
c                  if (vtot(k).gt.vescgc(ll))then
c                     countesc=countesc+1.0
c                     if(mod(countesc,1000).eq.0)then
c                     write(*,*) 'escape', k,vtot(k),vescgc(ll),
c     &                  bigmass(k),
c     &                    smallmass,eccprox(k)
c     endif
c               endif


c     endif
                  jay=0
                  uphi=0
                  Z1=0
                  Z2=0
                  Z3=0
                  r=0
                  mu=0

               enddo


               write(*,*) bound/k,vtotmax,bigmass(n3),s2para(n3),q
               allprob=(allprob*bound/k)
               write(*,*) j, allprob
               if(allprob.le.1e-10)then
                  goto 5000
               endif


            enddo
 5000       continue

            write(9,*) kkk,minit,bound/rn3,allprob,vescgc(ll),
     &                 bigmass(n3)
            write(*,*) 'kkk',kkk,minit,bound/rn3,allprob,vescgc(ll),
     &                 bigmass(n3)
         enddo
      enddo



      do mhisto=0,2000
         histo(mhisto)=0.1
      enddo


      mhisto=0
      do massran=0,2000
         mhisto=mhisto+1
         do k=1,n3
            if((vtot(k).lt.massran+1).and.
     &           (vtot(k).ge.massran))then
               histo(mhisto)=histo(mhisto)+1.0
            endif
         enddo
      enddo
      mhisto=0
      do massran=0,2000
         mhisto=mhisto+1
         write(10,*) massran,histo(mhisto)

      enddo
                  
 999   format(i9,1x,10(1x,e10.4))

       close(1)
       end


        FUNCTION ran2(idum)
        INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        REAL ran2,AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END



      
