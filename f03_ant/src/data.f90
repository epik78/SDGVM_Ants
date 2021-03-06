!> @brief Collection of subroutines to read in data.
!! @details
!! @author Mark Lomas
!! @date July 2016

module data

use real_precision
use dims
use func
use input_methods

implicit none

contains

!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!


!**********************************************************************!
!                                                                      !
!                          ex_clim :: data                             !
!                          ---------------                             !
!                                                                      !
! subroutine ex_clim(lat,lon,xlatf,xlatres,xlatresn,xlon0,     !
! xlonres,xlonresn,yr0,yrf,tmpv,humv,prcv,isite,year0,yearf,siteno,du) !
!                                                                      !                                                            !
!----------------------------------------------------------------------!
!> @brief Extract climate data from the climate database for the
!! nearest site.
!! @details ! Extract climate data from the climate database for the
!! nearest site to lat,lon, replace lat,lon with the nearest cell from
!! the database.
!!
!!            UNIX                DOS
!!
!!            ii(4)               ii(5)   for beginning of binary file
!!                                        records
!!            recl = 728          recl = 730     for binary climate
!!            recl = 577          recl = 578     for text map
!!
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine ex_clim(lat,lon,xlatf,xlatres,xlatresn,xlon0, &
 xlonres,xlonresn,yr0,yrf,tmpv,humv,prcv,isite,year0,yearf,siteno,du, &
 read_par)
!**********************************************************************!
real(dp) :: lat,lon,xlon0,xlatf,xlatres,xlonres,ans(12)
real(dp), dimension(500,12,31) :: swrv
real(dp) :: tmpv(300,12,31),humv(300,12,31),prcv(300,12,31)
integer :: year,year0,yearf,nrec,ncol,ans2(1000),siteno,i,du
integer :: nyears,yr0,mnth,day,isite,yrf,recl1,recl2
integer :: xlatresn,xlonresn,fno,read_par
character :: ii(4),jj(5),num*3
character(len=str_len) :: fname1,fname2,fname3,fname4
!----------------------------------------------------------------------!

if (du==1) then
  recl1 = 730
  recl2 = 6*xlonresn
else
  recl1 = 728
  recl2 = 6*xlonresn + 1
endif

nyears = yearf - year0 + 1

nrec = int((xlatf - lat)/xlatres + 1.0)
ncol = int((lon - xlon0)/xlonres + 1.0)

if ((nrec<=xlatresn).and.(ncol<=xlonresn)) then

lat = xlatf - xlatres/2.0 - real(nrec - 1)*xlatres
lon = xlon0 + xlonres/2.0 + real(ncol - 1)*xlonres

fno = 90

open(fno+1,file=trim(inp%dirs%climate)//'/maskmap.dat', &
 access='DIRECT',recl=recl2,form='formatted',status='OLD')

read(fno+1,'(96i6)',rec=nrec) (ans2(i),i=1,xlonresn)
close(fno+1)
siteno = ans2(ncol)

isite = 0
print*,siteno
if (siteno>0) then

  isite = 1

  write(num,'(i3.3)') (siteno-1)/100
  write(fname1,'(100a)') trim(inp%dirs%climate)//'/tmp_',num
  write(fname2,'(100a)') trim(inp%dirs%climate)//'/hum_',num
  write(fname3,'(100a)') trim(inp%dirs%climate)//'/prc_',num
  write(fname4,'(100a)') trim(inp%dirs%climate)//'/swr_',num

  siteno = mod(siteno-1,100) + 1

  open(fno+1,file=fname1,access='direct',recl=recl1,form='unformatted',status='old')
  open(fno+2,file=fname2,access='direct',recl=recl1,form='unformatted',status='old')
  open(fno+3,file=fname3,access='direct',recl=recl1,form='unformatted',status='old')
  open(fno+4,file=fname4,access='direct',recl=recl1,form='unformatted',status='old')

  if (du==1) then
    do year=yr0,yrf
      read(fno+1,rec=(siteno-1)*nyears+year-year0+1) jj, &
 ((tmpv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
      read(fno+2,rec=(siteno-1)*nyears+year-year0+1) jj, &
 ((humv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
      read(fno+3,rec=(siteno-1)*nyears+year-year0+1) jj, &
 ((prcv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
      do mnth=1,12
        do day=1,30
          prcv(year-yr0+1,mnth,day) = int(real(prcv(year-yr0+1,mnth,day))/10.0 + 0.5)
        enddo
      enddo
    enddo
  else
    do year=yr0,yrf
      read(fno+1,rec=(siteno-1)*nyears+year-year0+1) ii, &
 ((tmpv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
      read(fno+2,rec=(siteno-1)*nyears+year-year0+1) ii, &
 ((humv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
      read(fno+3,rec=(siteno-1)*nyears+year-year0+1) ii, &
 ((prcv(year-yr0+1,mnth,day),day=1,30),mnth=1,12)
      do mnth=1,12
        do day=1,30
          prcv(year-yr0+1,mnth,day) = int(real(prcv(year-yr0+1,mnth,day))/10.0 + 0.5)
        enddo
      enddo
    enddo
  endif

  close(fno+1)
  close(fno+2)
  close(fno+3)
  close(fno+4)

endif

else
  siteno = 0
endif

do mnth=1,12
  ans(mnth) = 0.0
  do day=1,30
    ans(mnth) = ans(mnth) + real(tmpv(1,mnth,day))/100.0
  enddo
enddo

end subroutine ex_clim





!**********************************************************************!
!                                                                      !
!                          ex_clim_site :: data                        !
!                          --------------------                        !
!                                                                      !
!                                                                      !
! subroutine ex_clim_site(yr0,yrf,tmpv,humv,prcv,year0,yearf)  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Extract climate data from the climate database for the
!! nearest site.
!! @details Extract climate data from the climate database for the
!! nearest site to lat,lon, replace lat,lon with the nearest cell from
!! the database.
!!
!!            UNIX                DOS
!!
!!            ii(4)               ii(5)   for beginning of binary file
!!                                        records
!!            recl = 728          recl = 730     for binary climate
!!            recl = 577          recl = 578     for text map
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine ex_clim_site(yr0,yrf,tmpv,humv,prcv,year0,yearf)
!----------------------------------------------------------------------!
real(dp) :: tmp,prc,hum
real(dp) :: tmpv(300,12,31),humv(300,12,31),prcv(300,12,31)
integer :: year,year0,yearf,yr0,mnth,day,yrf,iyear,imnth,iday
!**********************************************************************!

open(91,file=trim(inp%dirs%climate)//'/site.dat')

do year=year0,yearf
  do mnth=1,12
    do day=1,no_days(year,mnth,0)
      read(91,*) iyear,imnth,iday,tmp,prc,hum
      if ((iday/=day).or.(imnth/=mnth).or.(iyear/=year)) then
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,*) 'Error in climate data file',year,mnth,day
        stop
      endif

      if ((year<=yrf).and.(year>=yr0)) then
        tmpv(year-yr0+1,mnth,day) = tmp*100.0
        prcv(year-yr0+1,mnth,day) = prc*10.0
        humv(year-yr0+1,mnth,day) = hum*100.0
      endif
    enddo
  enddo
enddo

close(91)

end subroutine ex_clim_site





!**********************************************************************!
!                                                                      !
!                          ex_soil :: data                             !
!                          ---------------                             !
!                                                                      !
! subroutine ex_soil(lat,lon,sol_chr2,du,l_soil)                !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Extract % sand % silt bulk density and depth from soils
!! database.
!! @details Extract % sand % silt bulk density and depth from soils
!! database.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine ex_soil(lat,lon,sol_chr2,du,l_soil)
!----------------------------------------------------------------------!
real(dp) :: lat,lon,lon0,latf,latr,lonr,xlat,xlon,sol_chr2(10)
integer :: row,col,recn,recl1,du,latn,lonn,i,ii,jj,kode
integer, dimension(4,4) :: indx1,indx2,indx3,indx4,indx5
integer, dimension(4,4) :: indx6,indx7,indx8,indx9,indx10
real(dp),dimension(4,4) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8,xx9,xx10
real(dp) :: ynorm,xnorm,rrow,rcol
real(dp) :: ans
logical ::  l_soil(20)

if (du==1) then
! This works for ftn95
!  recl1 = 16+8*10
  recl1 = 16+10*10
else
  recl1 = 16+10*10+1
endif

open(99,file=trim(inp%dirs%soil)//'/readme.dat',status='old',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Soil data file does not exist.'
  write(*,'('' "'',A,''/readme.dat"'')') trim(inp%dirs%soil)
  stop
endif

read(99,*)
read(99,*) latf,lon0
read(99,*)
read(99,*) latr,lonr
read(99,*)
read(99,*) latn,lonn
close(99)

!----------------------------------------------------------------------!
! Find the real(dp) :: row col corresponding to lat and lon.                  !
!----------------------------------------------------------------------!
rrow = (latf - lat)/latr
rcol = (lon - lon0)/lonr

ynorm = rrow - real(int(rrow))
xnorm = rcol - real(int(rcol))
!----------------------------------------------------------------------!

open(99,file=trim(inp%dirs%soil)//'/data.dat',status='old', &
 form='formatted',access='direct',recl=recl1,iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Soil data-base.'
  write(*,*) 'File does not exist:',trim(inp%dirs%soil),'/data.dat'
  write(*,*) 'Or record length missmatch.'
  stop
endif

do ii=1,4
  do jj=1,4
    row = int(rrow)+jj-1
    col = int(rcol)+ii-1
    if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) then
      recn = (row-1)*lonn + col
      read(99,'(f7.3,f9.3,10f10.4)',rec=recn) xlat,xlon,(sol_chr2(i),i=1,10)
      xx1(ii,jj) = sol_chr2(1)
      xx2(ii,jj) = sol_chr2(2)
      xx3(ii,jj) = sol_chr2(3)
      xx4(ii,jj) = sol_chr2(4)
      xx5(ii,jj) = sol_chr2(5)
      xx6(ii,jj) = sol_chr2(6)
      xx7(ii,jj) = sol_chr2(7)
      xx8(ii,jj) = sol_chr2(8)
      xx9(ii,jj) = sol_chr2(9)
      xx10(ii,jj) = sol_chr2(10)
      if (sol_chr2(1)<0.0) then
        indx1(ii,jj) = 0
      else
        indx1(ii,jj) = 1
      endif
      if (sol_chr2(2)<0.0) then
        indx2(ii,jj) = 0
      else
        indx2(ii,jj) = 1
      endif
      if (sol_chr2(3)<0.0) then
        indx3(ii,jj) = 0
      else
        indx3(ii,jj) = 1
      endif
      if (sol_chr2(4)<0.0) then
        indx4(ii,jj) = 0
      else
        indx4(ii,jj) = 1
      endif
      if (sol_chr2(5)<0.0) then
        indx5(ii,jj) = 0
      else
        indx5(ii,jj) = 1
      endif
      if (sol_chr2(6)<0.0) then
        indx6(ii,jj) = 0
      else
        indx6(ii,jj) = 1
      endif
      if (sol_chr2(7)<0.0) then
        indx7(ii,jj) = 0
      else
        indx7(ii,jj) = 1
      endif
      if (sol_chr2(8)<0.0) then
        indx8(ii,jj) = 0
      else
        indx8(ii,jj) = 1
      endif
      IF (sol_chr2(9).LT.0.0d0) THEN
        indx9(ii,jj) = 0
      ELSE
        indx9(ii,jj) = 1
      ENDIF
      IF (sol_chr2(10).LT.0.0d0) THEN
        indx10(ii,jj) = 0
      ELSE
        indx10(ii,jj) = 1
      ENDIF
    else
      indx1(ii,jj) = -1
      indx2(ii,jj) = -1
      indx3(ii,jj) = -1
      indx4(ii,jj) = -1
      indx5(ii,jj) = -1
      indx6(ii,jj) = -1
      indx7(ii,jj) = -1
      indx8(ii,jj) = -1
      indx9(ii,jj) = -1
      indx10(ii,jj) = -1
    endif
  enddo
enddo

close(99)

if ((indx1(2,2)==1).or.(indx1(2,3)==1).or.(indx1(3,2)==1).or. &
 (indx1(3,3)==1)) then
  call bi_lin(xx1,indx1,xnorm,ynorm,ans)
  sol_chr2(1) = ans
  l_soil(1) = .true.
else
  l_soil(1) = .false.
endif
if ((indx2(2,2)==1).or.(indx2(2,3)==1).or.(indx2(3,2)==1).or. &
 (indx2(3,3)==1)) then
  call bi_lin(xx2,indx2,xnorm,ynorm,ans)
  sol_chr2(2) = ans
  l_soil(2) = .true.
else
  l_soil(2) = .false.
endif
if ((indx3(2,2)==1).or.(indx3(2,3)==1).or.(indx3(3,2)==1).or. &
 (indx3(3,3)==1)) then
  call bi_lin(xx3,indx3,xnorm,ynorm,ans)
  sol_chr2(3) = ans
  l_soil(3) = .true.
else
  l_soil(3) = .false.
endif
if ((indx4(2,2)==1).or.(indx4(2,3)==1).or.(indx4(3,2)==1).or. &
 (indx4(3,3)==1)) then
  call bi_lin(xx4,indx4,xnorm,ynorm,ans)
  sol_chr2(4) = ans
  l_soil(4) = .true.
else
  l_soil(4) = .false.
endif
if ((indx5(2,2)==1).or.(indx5(2,3)==1).or.(indx5(3,2)==1).or. &
 (indx5(3,3)==1)) then
  call bi_lin(xx5,indx5,xnorm,ynorm,ans)
  sol_chr2(5) = ans
  l_soil(5) = .true.
else
  l_soil(5) = .false.
endif
if ((indx6(2,2)==1).or.(indx6(2,3)==1).or.(indx6(3,2)==1).or. &
 (indx6(3,3)==1)) then
  call bi_lin(xx6,indx6,xnorm,ynorm,ans)
  sol_chr2(6) = ans
  l_soil(6) = .true.
else
  l_soil(6) = .false.
endif
if ((indx7(2,2)==1).or.(indx7(2,3)==1).or.(indx7(3,2)==1).or. &
 (indx7(3,3)==1)) then
  call bi_lin(xx7,indx7,xnorm,ynorm,ans)
  sol_chr2(7) = ans
  l_soil(7) = .true.
else
  l_soil(7) = .false.
endif
if ((indx8(2,2)==1).or.(indx8(2,3)==1).or.(indx8(3,2)==1).or. &
 (indx8(3,3)==1)) then
  call bi_lin(xx8,indx8,xnorm,ynorm,ans)
  sol_chr2(8) = ans
  l_soil(8) = .true.
else
  l_soil(8) = .false.
endif
IF ((indx9(2,2).EQ.1).OR.(indx9(2,3).EQ.1).OR.(indx9(3,2).EQ.1).or. &
 (indx9(3,3).EQ.1)) THEN
  CALL BI_LIN(xx9,indx9,xnorm,ynorm,ans)
  sol_chr2(9) = ans
  l_soil(9) = .true.
ELSE
  l_soil(9) = .false.
ENDIF
IF ((indx10(2,2).EQ.1).OR.(indx10(2,3).EQ.1).OR.(indx10(3,2).EQ.1).or. &
 (indx10(3,3).EQ.1)) THEN
  CALL BI_LIN(xx10,indx10,xnorm,ynorm,ans)
  sol_chr2(10) = ans
  l_soil(10) = .true.
ELSE
  l_soil(10) = .false.
ENDIF

end subroutine ex_soil





!**********************************************************************!
!                                                                      !
!                          ex_lu :: data                               !
!                          -------------                               !
!                                                                      !
! subroutine ex_lu(lat,lon,luse,yr0,yrf,du)                            !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Extract land use.
!! @details Extract land use from land use soils database for the
!! nearest site to lat,lon, replace lat,lon with the nearest cell from
!! the database.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine ex_lu(lat,lon,luse,yr0,yrf,du)
!**********************************************************************!
real(dp) :: lat,lon,lon0,latf,latr,lonr,xlat,xlon
integer :: i,n,j,du,latn,lonn,row,col,recn,kode
integer :: luse(max_years),yr0,yrf,rep,years(1000),lu(1000),nrecl
character(len=str_len) :: st1
!**********************************************************************!

open(99,file=trim(inp%dirs%land_use)//'/readme.dat')
read(99,*)
read(99,*) latf,lon0
read(99,*)
read(99,*) latr,lonr
read(99,*)
read(99,*) latn,lonn
read(99,*)
read(99,'(A)') st1
close(99)

n = n_fields(st1)
call ST2ARR(st1,years,1000,n)

if (du==1) then
  nrecl = 16+n*3
else
  nrecl = 16+n*3+1
endif

lu(1) = 0
open(99,file=trim(inp%dirs%land_use)//'/landuse.dat',status='old', &
 form='formatted',access='direct',recl=nrecl,iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'File does not exist:',trim(inp%dirs%land_use)//'/landuse.dat'
  stop
endif

row = int((latf - lat)/latr + 1.0)
col = int((lon - lon0)/lonr + 1.0)

recn = (row-1)*lonn + col
if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=recn) xlat,xlon,(lu(i),i=1,n)

if (lu(1)==0) then
  row = int((latf - (lat + latr))/latr + 1.0)
  col = int(((lon) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat - latr))/latr + 1.0)
  col = int(((lon) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat))/latr + 1.0)
  col = int(((lon + lonr) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat))/latr + 1.0)
  col = int(((lon - lonr) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat + latr))/latr + 1.0)
  col = int(((lon + lonr) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat - latr))/latr + 1.0)
  col = int(((lon + lonr) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat + latr))/latr + 1.0)
  col = int(((lon - lonr) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

if (lu(1)==0) then
  row = int((latf - (lat - latr))/latr + 1.0)
  col = int(((lon - lonr) - lon0)/lonr + 1.0)
  recn = (row-1)*lonn + col
  if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) &
 read(99,'(f7.3,f9.3,1000i3)',rec=int(min(recn,latn*lonn))) xlat,xlon,(lu(i),i=1,n)
endif

close(99)

rep = lu(1)
j = 1
do i=1,n
  if (yr0>=years(j+1)) then
    rep = lu(i)
    j = i
  endif
enddo

do i=1,yrf-yr0+1
  if ((i+yr0-1>=years(j+1)).and.(j<n)) then
    j = j+1
    rep = lu(j)
  endif
  luse(i) = rep
enddo

end subroutine ex_lu





!**********************************************************************!
!                                                                      !
!                          ex_clu :: data                              !
!                          --------------                              !
!                                                                      !
!  subroutine EX_CLU(lat,lon,nft,lutab,cluse,yr0,yrf,du,l_lu)   !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details  Extract land use from land use database for the
!! nearest site to lat,lon, replace lat,lon with the nearest cell from
!! the database. The data contains the percent (in byte format) for
!! each pft from 2 to nft.
!! First reads the readme.dat from the land use dataset with info such 
!! as resolution,years available and number of classes.Description of
!! classes and proportion assigned to model ft are ignored.
!! Land use must be written per class per year in vector format (i3)
!! 0-100 with 255 for water.Direction is West to East,North to South.
!! lutab(landcover classes,nft) is read from setup file.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine ex_clu(lat,lon,nft,lutab,cluse,yr0,yrf,du,l_lu)
!**********************************************************************!
real(dp) :: lat,lon,lon0,latf,latr,lonr,classprop(255)
real(dp) :: cluse(max_cohorts,max_years),lutab(255,100),ans
real(dp) :: ftprop(max_cohorts),rrow,rcol,xx(4,4),xnorm,ynorm
integer :: i,n,j,du,latn,lonn,row,col,recn,k,x,nft,ift
integer :: ii,jj,indx(4,4),yr0,yrf,years(1000),nrecl
integer :: classes(1000),nclasses,kode
character(len=str_len) :: st1,st2,st3
logical :: l_lu
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! read in the readme file 'readme.dat'.                                !
!----------------------------------------------------------------------!
open(99,file=trim(inp%dirs%land_use)//'/readme.dat',&
 status='old',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Land use file does not exist.'
  write(*,'('' "'',A,''/readme.dat"'')') trim(inp%dirs%land_use)
  stop
endif

read(99,*) st1
st2='CONTINUOUS'
if (stcmp(st1,st2)==0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'landuse is not a continuous field ?'
  write(*,*) 'readme.dat should begin with CONTINUOUS'
  stop
endif
read(99,*)
read(99,*) latf,lon0
read(99,*)
read(99,*) latr,lonr
read(99,*)
read(99,*) latn,lonn
read(99,*)
read(99,'(A)') st1
n = n_fields(st1)
call ST2ARR(st1,years,1000,n)
read(99,*)
read(99,'(A)') st1
close(99)
nclasses = n_fields(st1)
call ST2ARR(st1,classes,1000,nclasses)

!----------------------------------------------------------------------!
if (du==1) then
!  This works for ftn95
!  nrecl = 3
  nrecl = 5
else
  nrecl = 4
endif

if ((n>1).and.(yr0<years(1))) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Can''t start running in ',yr0,&
 ' since landuse map begin in ',years(1)
  stop
endif

!     look for the first year
j=1

 10   continue
if ((j<n).and.(years(j)<yr0)) then
   j = j + 1
   goto 10
endif

!----------------------------------------------------------------------!
! Find the real(dp) :: row col corresponding to lat and lon.                  !
!----------------------------------------------------------------------!
rrow = (latf - lat)/latr
rcol = (lon - lon0)/lonr

ynorm = rrow - real(int(rrow))
xnorm = rcol - real(int(rcol))
!----------------------------------------------------------------------!

j=1
do i=1,yrf-yr0+1
  if ((i==1).or.((i+yr0-1)==years(j))) then
    st2=in2st(years(j))
    call STRIPB(st2)
    j=j+1

    do k=1,nclasses
      classprop(classes(k)) = 0

      st3=in2st(classes(k))
      call STRIPB(st3)
      open(99,file= &
 trim(inp%dirs%land_use)//'/cont_lu-'//st3(1:blank(st3))//'-'//st2(1:4)//'.dat', & 
 status='old',form='formatted',access='direct',recl=nrecl,iostat=kode)
      if (kode/=0) then
        write(*,'('' PROGRAM TERMINATED'')')
        write(*,*) 'Land Use data-base.'
        write(*,*) 'File does not exist:'
        write(*,*) trim(inp%dirs%land_use)//'/cont_lu-',st3(1:blank(st3)),'-',st2(1:4),'.dat'
        stop
      endif

      do ii=1,4
        do jj=1,4
          row = int(rrow)+jj-1
          col = int(rcol)+ii-1
          if ((row>=1).and.(row<=latn).and.(col>=1).and.(col<=lonn)) then
            recn = (row-1)*lonn + col
            read(99,'(i3)',rec=recn) x 
            xx(ii,jj) = real(x)
            if (x<200) then
              indx(ii,jj) = 1
            else
              indx(ii,jj) = 0
            endif
          else
            indx(ii,jj) = -1
          endif
        enddo
      enddo

      call bi_lin(xx,indx,xnorm,ynorm,ans)
      
      x = int(ans+0.5)
      
      !This allows only specific crops
      !IF (k.NE.29) ans=0.
      
      classprop(classes(k)) = ans
      close(99)

    enddo ! end of loop over the classes
    
!----------------------------------------------------------------------!
! Now calculate the ftprop.
!----------------------------------------------------------------------!
    do ift=2,nft
      ftprop(ift)=0.0
      do k=1,nclasses
        ftprop(ift)=ftprop(ift)+lutab(classes(k),ift)*classprop(classes(k))/100.0
      enddo
    enddo

!----------------------------------------------------------------------!
! Calculate the bare soil.
!----------------------------------------------------------------------!
    if ((ftprop(2)<=100).and.(ftprop(2)>=0)) then
      ftprop(1)=100
      do ift=2,nft
        ftprop(1)=ftprop(1)-ftprop(ift)
      enddo
    endif

  endif ! finished reading

  do ift=1,nft
    cluse(ift,i) = ftprop(ift)
  enddo

enddo

!! TRENDY SDGVM method - assumes smoothed change in land-use between years specified in input dataset
!DO i=1,years(n)-yr0a
!  IF ( (i.EQ.1).OR.((i+yr0a-1).EQ.years(j1)) ) THEN
!    ij  = i
!    ij1 = ij + years(j1+1) - years(j1) 
!    j1  = j1 + 1
!  ELSE
!    DO ift=1,nft
!      cluse(ift,i) = cluse(ift,ij) + ( (real(i)-real(ij))/(real(ij1)-real(ij))* & 
! (cluse(ift,ij1) - cluse(ift,ij)) )
!    ENDDO
!  ENDIF
!ENDDO

if ((indx(2,2)==1).or.(indx(2,3)==1).or.(indx(3,2)==1).or.(indx(3,3)==1)) then
  l_lu = .true.
else
  l_lu = .false.
endif

end subroutine ex_clu





!**********************************************************************!
!                                                                      !
!                     lorc :: data                                     !
!                     ------------                                     !
!                                                                      !
! subroutine lorc(du,lat,lon,latdel,londel,xx)                  !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Determine if the site is a land site
!! @details Using a half an arc-second land sea mask. Determine whether
!! the site is land or sea based on the majority of the mask within
!! the limits of the gridcell.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine lorc(du,lat,lon,latdel,londel,xx)
!**********************************************************************!
logical :: xx
real(dp) :: lat,lon,latdel,londel,del,latf,lon0
character :: outc(6200)
integer :: i,j,k,x(7),ians(43300),sum1,n,col,row,check,nrecl,du,kode
!----------------------------------------------------------------------!

if (du==1) then
!  This works for ftn95
!  nrecl = 6172
  nrecl = 6172 + 2
else
  nrecl = 6172 + 1
endif

open(99,file=trim(inp%dirs%land_mask)//'/land_mask.dat', &
 form='formatted',recl=nrecl,access='direct',status='old',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Land sea mask.'
  write(*,*) 'Either the file doesn''t exist or there is a record length miss match.'
  write(*,*) 'Check that the correct DOS|UNIX switch is being used in the input file.'
  write(*,*) 'Land sea file :',trim(inp%dirs%land_mask)
  stop
endif

del = 1.0d0/60.0d0/2.0d0
latf = 90.0d0 - del/2.0d0
lon0 =-180.0d0 + del/2.0d0

col = int((lon - lon0)/del + 0.5d0)
row = int((latf - lat)/del + 0.5d0)

n = min((latdel/del-1.0d0)/2.0d0,(londel/del-1.0d0)/2.0d0)

!----------------------------------------------------------------------!
! Check nearest pixel for land.                                        !
!----------------------------------------------------------------------!
sum1 = 0
read(99,'(6172a)',rec=row) (outc(j),j=1,6172)

do j=1,6172
  call base72i(outc(j),x)
  do k=1,7
    ians(7*(j-1)+k) = x(k)
  enddo
enddo
if (ians(col)>0) sum1 = sum1 + 1
check = 1

!----------------------------------------------------------------------!
if (n>0) then
!----------------------------------------------------------------------!
! Check outward diagonals for land.                                    !
!----------------------------------------------------------------------!
  do i=1,n
    read(99,'(6172a)',rec=row+i) (outc(j),j=1,6172)
    do j=1,6172
      call base72i(outc(j),x)
      do k=1,7
        ians(7*(j-1)+k) = x(k)
      enddo
    enddo
    if (ians(col+i)>0) sum1 = sum1 + 1
    if (ians(col-i)>0) sum1 = sum1 + 1

    read(99,'(6172a)',rec=row-i) (outc(j),j=1,6172)
    do j=1,6172
      call base72i(outc(j),x)
      do k=1,7
        ians(7*(j-1)+k) = x(k)
      enddo
    enddo
    if (ians(col+i)>0) sum1 = sum1 + 1
    if (ians(col-i)>0) sum1 = sum1 + 1
    check = check + 4
  enddo
!----------------------------------------------------------------------!
endif

close(99)

if (real(sum1)>=(check+1)/2) then
  xx = .true.
else
  xx = .false.
endif
xx = .true.

end subroutine lorc





!**********************************************************************!
!                                                                      !
!                          base72I :: data                             !
!                          ---------------                             !
!                                                                      !
! subroutine base72i(c,x)                                              !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine base72i(c,x)
!**********************************************************************!
character :: c
integer :: x(7),i
!----------------------------------------------------------------------!

i=ichar(c)-100
x(1)=i/64
i=i-x(1)*64
x(2)=i/32
i=i-x(2)*32
x(3)=i/16
i=i-x(3)*16
x(4)=i/8
i=i-x(4)*8
x(5)=i/4
i=i-x(5)*4
x(6)=i/2
i=i-x(6)*2
x(7)=i

end subroutine base72i





!**********************************************************************!
!                                                                      !
!                          n7 :: data                                  !
!                          ----------                                  !
!                                                                      !
! Returns the value of a seven digit binary number given as a seven    !
! dimensional array of 0's and 1's                                     !
!                                                                      !
!    integer function n7(x)                                            !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
integer function n7(x)
!----------------------------------------------------------------------!
integer :: x(7)
!----------------------------------------------------------------------!

n7 = 64*x(1)+32*x(2)+16*x(3)+8*x(4)+4*x(5)+2*x(6)+x(7)
 
!**********************************************************************!
end function n7





!**********************************************************************!
!                                                                      !
!                          bi_lin :: data                              !
!                          --------------                              !
!                                                                      !
!                                                                      !
!  subroutine bi_lin(xx,indx,xnorm,ynorm,ans)                          !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Performs bilinear interpolation between four points.
!! @details Performs bilinear interpolation between four points, the
!! normalised
!! distances from the point (1,1) are given by 'xnorm' and 'ynorm'.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine bi_lin(xx,indx,xnorm,ynorm,ans)
!----------------------------------------------------------------------!
real(dp) :: xx(4,4),xnorm,ynorm,ans,av
integer :: indx(4,4),iav,ii,jj
!**********************************************************************!

!----------------------------------------------------------------------!
! Fill in averages if necessary.                                       !
!----------------------------------------------------------------------!
do ii=2,3
  do jj=2,3
    if (indx(ii,jj)/=1) then
      av = 0.0
      iav = 0
      if (indx(ii+1,jj)==1) then
        av = av + xx(ii+1,jj)
        iav = iav + 1
      endif
      if (indx(ii-1,jj)==1) then
        av = av + xx(ii-1,jj)
        iav = iav + 1
      endif
      if (indx(ii,jj+1)==1) then
        av = av + xx(ii,jj+1)
        iav = iav + 1
      endif
      if (indx(ii,jj-1)==1) then
        av = av + xx(ii,jj-1)
        iav = iav + 1
      endif
      if (indx(ii+1,jj+1)==1) then
        av = av + xx(ii+1,jj+1)
        iav = iav + 1
      endif
      if (indx(ii-1,jj-1)==1) then
        av = av + xx(ii-1,jj-1)
        iav = iav + 1
      endif
      if (indx(ii+1,jj-1)==1) then
        av = av + xx(ii+1,jj-1)
        iav = iav + 1
      endif
      if (indx(ii-1,jj+1)==1) then
        av = av + xx(ii-1,jj+1)
        iav = iav + 1
      endif
      if (iav>0) then
        xx(ii,jj) = av/real(iav)
      endif
    endif
  enddo
enddo

!----------------------------------------------------------------------!
! Bilinear interpolation.                                              !
!----------------------------------------------------------------------!
ans = xx(2,2)*(1.0-xnorm)*(1.0-ynorm) + xx(3,2)*xnorm*(1.0-ynorm) + &
 xx(2,3)*(1.0-xnorm)*ynorm + xx(3,3)*xnorm*ynorm

!----------------------------------------------------------------------!
! nearest pixel.                                                       !
!----------------------------------------------------------------------!
! ans = xx(int(xnorm+2.5),int(ynorm+2.5))
!----------------------------------------------------------------------!

end subroutine bi_lin





!**********************************************************************!
!                                                                      !
!                          readCO2 :: data                             !
!                          ---------------                             !
!                                                                      !
! subroutine readCO2(yr0,yrf,co2)                                      !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief Reads co2 values for each year
!! @details It will read from the co2 file the concentration vector 
!! with starting and end year as produced by the read_input_file sub.
!! If the co2 file has only one value then it will be used for all the
!! years meaning constant co2 concentrations.
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine readCO2(yr0,yrf,co2,daily_co2,nyears)
!----------------------------------------------------------------------!
real(dp) :: co2(max_years,12,31),ca
integer :: yr0,yrf,norecs,year,const,kode
integer :: yra,yrfa,prev_year,nyears
integer :: mnth,yr0a,daily_co2,day
logical :: co2spin
!**********************************************************************!

open(98,file=trim(inp%dirs%co2),status='OLD',iostat=kode)
if (kode/=0) then
  write(*,'('' PROGRAM TERMINATED'')')
  write(*,*) 'Co2 file does not exist.'
  write(*,'('' "'',A,''"'')') trim(inp%dirs%co2)
  stop
endif

      norecs = 0
      const  = 0
      co2spin = .FALSE.
      yrfa    = yrf
      if((inp%run%co2_constant.lt.0.0).and.(inp%run%spinup_length.gt.0)) then 
        if (inp%run%spinup_length.lt.nyears) then
          yr0a = yr0 - inp%run%spinup_length
        else
          co2spin = .TRUE.
        endif
      else 
        yr0a = yr0
      endif 

10    CONTINUE

!      print*, yr0, yrf

!      print*, daily_co2
      IF(daily_co2==1) THEN 
	READ(98,*,end=99) year,mnth,day,ca 
        if(co2spin.and.(norecs==0)) then
          yr0a = year
          yrfa = yr0a + inp%run%spinup_length - 1
        endif
	IF ((year==yr0a).and.(mnth==1).and.(day==1)) &
 prev_year = year-1
!     print*, prev_year, norecs      
        IF ((year.GE.yr0a).AND.(year.LE.yrfa)) THEN 
          IF(prev_year==(year-1)) norecs = norecs + 1 
          co2(norecs,mnth,day) = ca 
        ENDIF 
        prev_year = year 
      ELSE 
        READ(98,*,end=99) year,ca
!        print*, year,ca
        if(co2spin.and.(norecs==0)) then
          yr0a = year
          yrfa = yr0a + inp%run%spinup_length - 1
          print*, 'co2spin',yr0a,yrfa
        endif
	IF ((year.GE.yr0a).AND.(year.LE.yrfa)) THEN
          norecs = norecs + 1
          co2(norecs,:,:) = ca
        ENDIF
      ENDIF 
      const = const + 1
      GOTO 10
99    CONTINUE

      CLOSE(98)

      IF (const==1) THEN
	DO year=yr0,yrf
	  co2(year-yr0+1,:,:) = ca
	ENDDO
      ENDIF

end subroutine readCO2





!**********************************************************************!
!                                                                      !
!                         landuse1 :: data                             !
!                         ----------------                             !
!                                                                      !
!----------------------------------------------------------------------!
!> @brief
!! @details
!! @author Mark Lomas
!! @date Feb 2006
!----------------------------------------------------------------------!
subroutine landuse1(luse,yr0,yrf,fire_ant,harvest_ant)
!**********************************************************************!
integer :: yr0,yrf,year,early,rep,luse(max_years),i,use
logical :: fire_ant(max_years),harvest_ant(max_years)
!----------------------------------------------------------------------!

early = yrf+max_years

do i=1,max_years
  luse(i) = 0
enddo

do i=1,inp%land_use%n
  year = inp%land_use%year(i)
  if (year-yr0+1>0)  luse(year-yr0+1) = int(inp%land_use%map(i))
  if (year<early) rep = int(inp%land_use%map(i))
enddo

do i=1,max_years
  if (luse(i)>0) rep = luse(i)
  luse(i) = rep
enddo

return



!do i=1,max_years
!  luse(i) = 0
!enddo

!i = 1
!10    continue
!  read(fid_input,*,end=20) year,use
!  print*,use
!  if (year-yr0+1>0)  luse(year-yr0+1) = use
!  if (year<early) rep = use
!  i = i+1
!goto 10
!20    continue

!do i=1,max_years
!  if (luse(i)>0) rep = luse(i)
!  luse(i) = rep
!enddo


!return

! Anthony's updated version this seems to cause a problem.
!      early = yrf+max_years
!      luse(:) = 1000
      
!      i = 1
!11    CONTINUE
!        READ(98,*,end=21) year,use
!        IF (year-yr0+1.GT.0)  luse(year-yr0+1) = use
!        IF ((year.LT.early).AND.(use.GT.0)) rep = use
!        !early = year
!        i = i+1
!      GOTO 11
!21    CONTINUE

!      if(rep.LE.0) then
!        print*, 'ERROR:: Land use mapping equal to or below 0'
!        print*, 'you must specifiy at least one year with a cover type'
!      endif

!      DO i=1,max_years
!        If ((luse(i).GT.0).AND.(luse(i).LT.1000)) rep = luse(i)
!        If (luse(i).LT.0) fire_ant(i)    = .TRUE.
!        If (luse(i)==0) harvest_ant(i) = .TRUE.
!        luse(i) = rep
!      ENDDO

end subroutine landuse1

end module data
