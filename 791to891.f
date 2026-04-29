c     791to891.f: compile with ifort/gfortran 791to891.f -o 791to891.exe -> same fort.891
c---------------------------------------------------------------------------------------------

c     This file converts the.791 file into a .891 file which considers the clusters as formed when
c     they have a negative binding energy and the cluster nucleons are freezed out
c files 791(MST), fort.793(SACA) [identical for anti-clusters fort.781(MST), fort.783(SACA)]:
c     column 1:    current entry number of baryon in the recording vector [1,2,3,...Nparticles]
c     column 2:    'id' of baryons in QMD according to the QMD iso(i) vector - defined in iso.data [0=n, 1=p, 17=Lambda, 18=Sigma0]
c     column 3-5:  3-momentum of baryon in calculational A+A frame : (Px,Py,Pz) [GeV/c]
c     column 6-8:  position of baryon (X;Y;Z) in calculational A+A frame [fm]
c     column 9:    baryon mass [GeV] (i.e. mass of p or n or Lambda or Sigma0)
c     column 10:   number of cluster (Ncluster) to which this baryon belongs to - in MST(791) or SACA(793) mode
c     column 11:   size of cluster A_cluster to which this baryon belongs to - in MST(791) or SACA(793) mode
c     column 12:   original position "J" of baryon in the PHSD vector (integer number)
C                  needed for study of cluster stability and for sinchronization of cluster output with phsd.dat
c    *column 13:   MST mode: identical to column 10 (actually not needed);
C                  SACA mode: number of cluster (Ncluster) to which this baryon belongs to in the MST 'precalculation' for SACA
c    *column 14:   MST mode: identical to column 11 (actually not needed);
C                  SACA mode: size of cluster A to which this baryon belongs to in the MST 'precalculation' for SACA
C     column 15:   infor from PHQMD: ID(J,6) - if : =-1 - nucleon from projectile, =+1 - from target (wo collision);  = other - made collisions
C     column 16:   infor from PHQMD: ID(J,5) - production channel of baryon according to the PHSD list (cf. PHSD Manual)
C     column 17:   infor from PHQMD: P(4,J) - production time of nucleon in fm/c (needed for reconstruction of cluster stability)
C     column 18:   binding energy of cluster per baryon [GeV]
C     used format: format(1x,2I4,7(1XE11.4),2(1XI3),3(1XI10),2(1XI10),(2XE11.4))

c     changed in order to remove some bugs 4.12.25  

      character*130 filename91,filename93,filename95,filename97

      parameter (imax1=1000)    !number of particles
      parameter (imn=200)       !number of num
      parameter (imt=60)        !number of timesteps
      parameter (ims=1000)      !size of cluster
      parameter (nentries=1000000)

      real emax,aman,tstart,dtstep, ebt
      real emm(imax1,imt,imn)
      real time(imt),xxt(imax1,3,imt,imn),ppt(imax1,3,imt,imn)
      integer ipos(nentries,imt),idb(ims)
      real frou(nentries)
      integer IMST
      integer kclus1SMt(imax1,imt,imn),lsize1SM(imax1,imt,imn)
      integer kclus1SMt2(imax1,imt,imn),lsize1SM2(imax1,imt,imn)
      Integer IPHSD4FRIGA(imax1,imt,imn),It(imt,imn)
      integer isot(imax1,imt,imn)
      integer IRUN(imt,imn),ISUB(imt,imn)
      integer nparticlest(imt,imn),ipo(imax1)
      integer IDfriga1(imax1,imt,imn),IDfriga2(imax1,imt,imn)
      real pfriga(imax1,imt,imn),ebind(imax1,imt,imn)
      real*8 rcluster
      real pfriga_iphsd(nentries)
      integer iphsd_last,itstep,iStatbOrBig,ilbo,iub,num

      logical :: file_exists
      integer :: tageyuk1,tageasy1,tagepair1
      real :: vasy01,eta_pairing1
c...      real :: eratio_check

cccccccccccccccccccccccccccccccccccccccc
c  read the files 791 and/or 793
cccccccccccccccccccccccccccccccccccccccccc
c     variables which can (have to) be changed

      OPEN(4,FILE='parameters.txt',STATUS='OLD')
      READ(4,*) num             ! number of parallel ensembles
      READ(4,*) itstep          ! number of timesteps
      READ(4,*) tstart          ! initial time to start SACA output
      READ(4,*) dtstep          ! size of timestep
      READ(4,*) filename91      ! input filename
      READ(4,*) filename95      ! output filename
      CLOSE(4)

      itstep=itstep+1

      write(*,*)'NUM:                      ',num
      write(*,*)'number of timesteps + 1:  ',itstep
      write(*,*)'start time:               ',tstart
      write(*,*)'size of timestep:         ',dtstep
      write(*,*)'input filename:           ',filename91
      write(*,*)'output filename:          ',filename95

      ilbo=0
      iub=0
      emax=0.0                  ! upper limit of bindungsenergy per nucleon in GeV
      iFutOrItin=0              ! =0: stabilize all future timesteps in itin, =1: only stabilize current timestep itin
      aman=.938
cccccccccccccccccccccccccccccccccccccccc
        do 108 imu=ilbo,iub

        INQUIRE(FILE=filename91, EXIST=file_exists)
        write(*,*)filename91, file_exists
        if(file_exists.eqv..false.)goto 909
        if(file_exists.eqv..true.)open(791,file=filename91,status='old')
        write(*,*)filename91
        open(891,file=filename95,status='unknown')

         do 109 iru=1,1000
            aev=aev+1

            IRUN(:,:)=0
            ISUB(:,:)=0
            time(:)=0.
            it(:,:)=0
            nparticlest(:,:)=0
            isot(:,:,:)=0.
            ppt(:,:,:,:)=0.
            xxt(:,:,:,:)=0.
            emm(:,:,:)=0.
            kclus1SMt(:,:,:)=0
            lsize1SM(:,:,:)=0
            IPHSD4FRIGA(:,:,:)=0
            kclus1SMt2(:,:,:)=0
            lsize1SM2(:,:,:)=0
            IDfriga1(:,:,:)=0
            IDfriga2(:,:,:)=0
            pfriga(:,:,:)=0.
            ebind(:,:,:)=0.
            frou(:)=0.
            pfriga_iphsd(:)=200.0
            iphsd_last=0
            ebt=0.

         do 300 Itst=1,itstep
         do 301 inum=1,num
            read(791,9333,end=999) IRUN(itst,inum),ISUB(itst,inum),b,
     &   time(itst),
     &   it(itst,inum),npr,nzpr,nta,nzta,frepp,iqmdeos

                 if(it(itst,inum).ne.itst)write(*,*)'error timestep',
     &   it(itst,inum),itst

9333  format(2I5,2(2XF10.4),5I5,E12.5,I5)
      read(791,653)tageyuk1,tageasy1,tagepair1,vasy01,
     &eta_pairing1

653   format(1x,3I4,2f10.4)
654   format(3I12,f24.12)

      read(791,*) nparticlest(itst,inum),iMST,LSIZE_POT2,rcluster

      do 4415 i = 1, NPARTICLESt(itst,inum) ! loop over FRIGA clusters (max. number = number of particles which have participated to the clusterisation)
      read(791,569)k,isot(i,itst,inum),(ppt(i,j,itst,inum),j=1,3),
     &(xxt(i,j,itst,inum),j=1,3),emm(i,itst,inum),
     & kclus1SMt(i,itst,inum),
     & lsize1SM(i,itst,inum),
     & IPHSD4FRIGA(i,itst,inum),
     & kclus1SMt2(i,itst,inum),
     & lsize1SM2(i,itst,inum), IDfriga1(i,itst,inum),
     & IDfriga2(i,itst,inum),pfriga(i,itst,inum),ebind(i,itst,inum)

c********************************
4415  enddo

      write(*,*)iru, inum,itst,time(itst)

301   continue ! end num
300   continue ! end time steps

 569  format(1x,2I4,7(1XE16.9),2(1XI3),3(1XI10),2(1XI10),2(1XE11.4))  ! for PHSD-PHQMD: format extended for coordinates

c************************************

c     end of reading 791/793

c************************************
c
c     freeze out time for each baryon id
c
c************************************

         do 110 inum=1,num
           do ktstlast=1,nparticlest(itstep,inum)
              iphsd_last=IPHSD4FRIGA(ktstlast,itstep,inum)
c... check for too many nentries
c...              if(iphsd_last.gt.nentries) then
c...                write(*,*) 'increase nentries',iphsd_last
c...                stop
c...              endif
c... end check
              pfriga_iphsd(iphsd_last)=pfriga(ktstlast,itstep,inum)
           enddo
           do itin=1,itstep     ! time steps into the past
              do ktst2=1,nparticlest(itin,inum)
                 frou(IPHSD4FRIGA(ktst2,itin,inum))
     &                =pfriga_iphsd(IPHSD4FRIGA(ktst2,itin,inum))
              enddo
           enddo

c************************************
c
c     search for the future of clusters
c
c************************************
            ipos(:,:)=0

            do itin=1,itstep    ! time steps into the future
               do ktst9=1,nparticlest(itin,inum) !loop over baryons at that time
                  ipos(IPHSD4FRIGA(ktst9,itin,inum),itin)=ktst9 !to allow trace back when the cluster disappears
               enddo
            enddo

      do 114 itin=2,itstep ! time steps into the future
c     tref=tstart+dtstep*(itin-1) ! change first timestep
c     tref=tstart+dtstep*(itin-2) ! change first timestep corrected
         tref=time(itin-1)        ! change first timestep exact time (referst to same timestep as "corrected")
                  do 113 ktst1=1,nparticlest(itin,inum)!loop over baryons at that time
                  irr=0
      idn=IPHSD4FRIGA(ktst1,itin,inum)!to shorten lines
      ktst1o=ipos(idn,itin-1)
      if(ipos(idn,itin-1).eq.0)goto 113 !if new id do not check clusters
      ido=IPHSD4FRIGA(ipos(IPHSD4FRIGA(ktst1,itin,inum),
     &itin-1),itin-1,inum)!id time step before (check for algorithm)
        if(idn.ne.ido)!have not the same PHSD ID
     &Write(*,*)'err',ido,IPHSD4FRIGA(ktst1,itin,inum)
c************************************
      if(ipos(idn,itstep).eq.0)goto 123! leave if this PHSD ID does not exist in the last time step
c     if(pfriga(ipos(idn,itstep),itstep,inum).gt.10+itin*5)goto 123! not yet freezed out
      if(frou(idn).gt.tref)goto 123 ! not yet freezed out
c************************************

      kln=kclus1SMt(ktst1,itin,inum)! cluster to which baryon belongs at time step itin
      lsn=lsize1SM(ktst1,itin,inum) !size of the cluster to which baryon belongs
      Enew=ebind(ktst1,itin,inum)!binding energy of cluster at time step itin
      klo=kclus1SMt(ktst1o,itin-1,inum)
      lso=lsize1SM(ktst1o,itin-1,inum)
      Eold=ebind(ktst1o,itin-1,inum)
      irr=0
      if(lso.eq.1)goto 839
      if(lso.ne.lsn)goto 839
      lsn_content=0
      do ktst7=1,nparticlest(itin-1,inum) ! check if the cluster content did not change from itin-1 to itin
         if(klo.ne.kclus1SMt(ktst7,itin-1,inum)) goto 840 ! check if baryon is part of cluster in itin-1
         ktst7_new = ipos(IPHSD4FRIGA(ktst7,itin-1,inum),itin) ! position of ktst7 in itin
         if(kln.ne.kclus1SMt(ktst7_new,itin,inum))goto 840 ! check if baryon is part of cluster in itin
         lsn_content=lsn_content+1
 840  enddo
      if(eold/enew.lt.0.and.lso.eq.lsn_content)irr=1 !=1 if sign of binding energy has changed
839   continue
c********************************
c      find the nucleons which were part of the cluster at itin-1
c      and identify what is their fate at itin
      mm1=0
      ebt=0.
               
       do ktst4=1,nparticlest(itin-1,inum)
          if(lsize1SM(ktst4,itin-1,inum).eq.1)goto 777
          if(klo.ne.kclus1SMt(ktst4,itin-1,inum)) goto 777 ! selcts baryons in the same cluster as ktst1 at itin-1
          ktst4_new = ipos(IPHSD4FRIGA(ktst4,itin-1,inum),itin) ! position of ktst4 in itin
          ebt=ebt+ebind(ktst4_new,itin,inum) ! total binding energy of the baryons from the fomer cluster in itin-1 in itin
          if(lsize1SM(ktst4_new,itin,inum).eq.1)goto 777 ! checks if baryon is a single baroyn in itin
          mm1=mm1+1      
           
 777   enddo
 778   continue

      icl=0
      if(lso.gt.lsn)icl=1
      if(icl.eq.1.and.mm1.eq.0)icl=2 ! baryons from former cluster (itin-1) are all single baryons in itin
c      if(icl.eq.1.and.ebt.ge.0)icl=2 ! baryons from former cluster (itin-1) don't have a negative binding energy in itin (in total)

      if(icl.eq.2.or.irr.eq.1)then

c********************************
c     reparation of the file 791
      iclu=0
      isurv=0
      jk=1
      idb(:)=0
      ipo(:)=0
      do 214 ktst5=1,nparticlest(itin-1,inum) ! look for the particle with the same PHSD ID
         if(kclus1SMt(ktst1o,
     &itin-1,inum).ne.kclus1SMt(ktst5,itin-1,inum)) goto 214!baryons do not belong to the same cluster

         isurv=isurv+1

         if(frou(IPHSD4FRIGA(ktst5,itin-1,inum)).gt.tref)goto 214 ! all baryons for the disintegrating cluster have to have an earlier freeze out time
c         if(frou(IPHSD4FRIGA(ktst5,itin-1,inum)).eq.0.)goto 214
         if(ebind(ktst5,itin-1,inum).gt.emax)goto 214
         if(ipos(IPHSD4FRIGA(ktst5,itin-1,inum),itstep).eq.0)goto 214 ! leave if this PHSD ID does not exist in the last time step
         if(ipos(IPHSD4FRIGA(ktst5,itin-1,inum),itin).eq.0)goto 214 ! leave if this PHSD ID does not exist in the current time step
         
         iclu=iclu+1
         idb(jk)=IPHSD4FRIGA(ktst5,itin-1,inum)
         ipo(jk)=ktst5
         jk=jk+1                !jk-1=nucleons in the disintegration cluster with a freeze out < cluster disintegration time

 214  enddo
      if(isurv.ne.iclu)goto 113

      if(iFutOrItin.eq.0) itstepStab=itstep ! stabilize cluster in all future times
      if(iFutOrItin.eq.1) itstepStab=itin ! stabilize cluster only in itin
      do 115 itime=itin,itstepStab ! go through all future times
         numcl=0
         
         do  ktst6=1,nparticlest(itime,inum) !do loop to check the number of clusters
            numcl=max0(kclus1SMt(ktst6,itime,inum),numcl)
         enddo
         iff=0

      do 116 ktst3=1,nparticlest(itime,inum) ! look for the particles with the same PHSD ID
         do jkk=1,jk-1
            if(idb(jkk).eq.IPHSD4FRIGA(ktst3,itime,inum))then !identifies the baryon with the same ID as cluster baryons before decay at that time step
               iff=iff+1
               
         if(iFutOrItin.eq.0) then ! check if baryon is part of a different cluster in itime: if yes, correct the size of this cluster for the missing baryon

         ktst3_old=0
         ktst3_old=ipos(IPHSD4FRIGA(ktst3,itime,inum),itime-1)
         if(ktst3_old.eq.0)goto 118
            
         do 117 ktst7=1,nparticlest(itime,inum)
         if(kclus1SMt(ktst7,itime,inum).ne.kclus1SMt(ktst3,itime,inum)) ! check if baryon is part of the cluster in itime
     &           goto 117
         ktst7_old=ipos(IPHSD4FRIGA(ktst7,itime,inum),itime-1)
         if(kclus1SMt(ktst7_old,itime-1,inum).eq.
     &        kclus1SMt(ktst3_old,itime-1,inum))goto 117 ! check if baryon was part of the cluster already in itime-1
         lsize1SM(ktst7,itime,inum)=lsize1SM(ktst7,itime,inum)-1
 117  enddo
 118  continue

      endif                     ! end check if baryon is part of different cluster in itime

      lsize1SM(ktst3,itime,inum)=lso ! replace lsize there
      kclus1SMt(ktst3,itime,inum)=numcl+1 ! gives the cluster its kclus
      ebind(ktst3,itime,inum)=eold ! replace lsize there

      endif
      enddo                     !selected cluster nucleons
 116  enddo                     !baryons
 115  enddo                     !future times

      endif                     !end reparation-condition
c     end reparation of the file 791
c********************************

123   continue

c********************************

113   continue !end particle loop
114       continue !end time loop
1111   continue !end part loop
110   continue !end num loop
      do Itst=1,itstep
      do inum=1,num
      write(891,9333) IRUN(itst,inum),ISUB(itst,inum),b,time(itst),
     &   it(itst,inum),npr,nzpr,nta,nzta,frepp,iqmdeos
      write(891,653)tageyuk1,tageasy1,tagepair1,vasy01,
     &eta_pairing1
      write(891,654) nparticlest(itst,inum),iMST,LSIZE_POT2,rcluster
      do i = 1, NPARTICLESt(itst,inum)
      write(891,569)i,isot(i,itst,inum),(ppt(i,j,itst,inum),j=1,3),
     &(xxt(i,j,itst,inum),j=1,3),emm(i,itst,inum),
     & kclus1SMt(i,itst,inum),
     & lsize1SM(i,itst,inum),
     & IPHSD4FRIGA(i,itst,inum),
     & kclus1SMt(i,itst,inum),
     & lsize1SM(i,itst,inum), IDfriga1(i,itst,inum),
     & IDfriga2(i,itst,inum),pfriga(i,itst,inum),ebind(i,itst,inum)
           enddo
           enddo
           enddo

 109   enddo                    ! end of sequential runs in the same file
 999   continue
 108  enddo                     ! end file loop

 909  continue

      end
