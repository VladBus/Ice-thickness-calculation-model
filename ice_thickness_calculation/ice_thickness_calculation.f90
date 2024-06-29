program Ice_Thickness_Calculation

  implicit none

  real :: tt(5000), tr(5000), ta(15,30,8), gist(500), temp, tpl, tps, qw, dt, tti, roi, teta, hs
  real :: tsur, tt0, rkm, rdel, rku, ttor, rkmn, grad0, grad1, rj0, rj1
  integer :: np(5000), inp(5000), km, i, id, ih, kmn, k, j, ktor, idel, np_prev, n, im, press
  character*2 cim, cid
  character*1 cdec
  character*4 cnp(5000)
  character*80 temper_a_1, temper_a_2

  call DeleteOldFile

  open (1,file='.\\data\\thickness_ini\\thick-ini-38-62.txt',status='old') ! Модуль чтения толщин
    read (1,'(i8)') km              ! количество маркеров
    do i=1,km
      read (1,*) np(i),tt(i),tr(i)  ! номер маркера, начальная толщина маркера (м) и толщина торосов на этом маркере (м)
      cnp(i)='1001'
    enddo
  close (1)

  call DataFiles(temper_a_1, temper_a_2)

  open (1,file='.\\data\\temper_a_ini\\temper_a_2.txt',status='old') ! Модуль чтения температуры
    do im=1,8                                    ! цикл по месяцам
      do id=1,30                                 ! цикл по дням
        do ih=1,8                                ! цикл по часам
          read (1,*) ta(im,id,ih)                ! темп. воздуха (град)
        enddo
      enddo
    enddo
  close (1)

  open (3,file='.\\result\\press_type\\press-type.txt',status='unknown')  !открываем файл по торосам
  open (4, file= '.\\result\\mon_gst\\mon-all.gst',status='unknown') !открываем общий файл по гистограмме

  ! ---константы---
  tpl=2.22                              ! теплопроводность льда Вт/(м*К)
  tps=0.26                              ! теплопроводность снега Вт/(м*К)
  qw=0.              ! поток тепла от воды к нижней поверхности льда (Вт)
  dt=10800.                             ! шаг по времени (с)
  tti=3.35e5                    ! удельная теплота плавления льда (Дж/кг)
  roi=920.                              ! плотность льда (кг/м**3)
  teta=-1.5                             ! температура замерзания воды
  ! ---константы---

  do im=1,8                ! цикл по месяцам
    do id=1,30             ! цикл по дням
      do ih=1,8            ! цикл по часам внутри каждого дня

        if(im.eq.1) cim='10'
        if(im.eq.2) cim='11'
        if(im.eq.3) cim='12'
        if(im.eq.4) cim='01'
        if(im.eq.5) cim='02'
        if(im.eq.6) cim='03'
        if(im.eq.7) cim='04'
        if(im.eq.8) cim='05'

        if (id.lt.10) then
          write (cid,'(a1,i1)') '0',id
          else
          write (cid,'(i2)') id
        endif

        if(ta(im,id,ih).ge.teta) go to 1

        do i=1,km                              ! цикл по маркерам

          if (tt(i).gt.0.1) hs=0.1*tt(i)       ! hs - толщина снега (м)
          if (tt(i).le.0.1) hs=0

          tsur=(tt(i)*ta(im,id,ih))/(tt(i)+0.03) ! темп. снежной-ледяной поверхности (град)

          tt0=tt(i)
          tt(i)=-hs*tpl/tps-qw*dt/tti/roi+ &    
					& sqrt((hs*tpl/tps+qw*dt/tti/roi)**2+ &
					& tt0**2+2*tpl*abs(tsur-teta)*dt/roi/tti- &
					& 2*hs*tpl/tps*(qw*dt/roi/tti-tt0)) ! формула расчета толщины льда (м)
          
        enddo                                  ! конец цикла по маркерам
      1 end do             ! конец цикла по часам

      ! Начало блока ввода новых маркеров (разрывы) и исключения старых маркеров (сжатия и торошения)

      ! ---функция рандомайзера---
      call random_number (temp)

      if(temp.le.0.25) press=0
      if(temp.gt.0.25) press=1
      if(temp.gt.0.82) press=2
      if(temp.gt.0.95) press=3
      ! ---функция рандомайзера---

      write (3,'(2(2x,a2),f6.3,i4)') cim,cid,temp,press ! запись в файл press-type

      if(press.eq.0.or.press.eq.1) then ! условия для сжатий 0,1
        rkm=km
        if(press.eq.0) rkm=rkm+0.02*rkm
        if(press.eq.1) rkm=rkm+0.01*rkm
        kmn=nint(rkm)
        rdel=kmn-km

        do k=km+1,kmn
          rku=k-km
          if(rku.le.0.25*rdel) tt(k)=0.03
          if(rku.gt.0.25*rdel.and.rku.le.0.75*rdel) tt(k)=0.04
          if(rku.gt.0.75*rdel) tt(k)=0.05
          tr(k)=0
          np(k)=k
          cnp(k)=cim//cid
          inp(k)=(im-1)*30+id
        enddo

        km=kmn
      endif                       ! конец условий для сжатий 0,1

      if (press.eq.2.or.press.eq.3) then ! условия для сжатий 2,3
        ktor=1
        ttor=0
        rkm=km
        if(press.eq.2) rkm=rkm-0.04*rkm
        if(press.eq.3) rkm=rkm-0.08*rkm
        kmn=nint(rkm)
        idel=km-kmn

        do 2 k=1,km
          if(ktor.eq.idel) go to 2

          if (press.eq.2.and.tt(k).le.0.20.or.press.eq.3.and.tt(k).le.0.40) then
            np_prev=np(k)
            np(k)=0
            ktor=ktor+1
            ttor=ttor+tt(k)

            if (im.eq.7.and.id.eq.9) then
              write (6,101) im,id,temp,press,km,kmn,idel,k,ktor,tt(k),np_prev,np(k),cnp(k)
              101 format (' im,id,temp,press,km,kmn,idel,k,ktor,tt,np_pr,np,cnp',2i6,2f6.3,5i6,f12.4,2i6,3x,a4)
            endif
          endif
        2 end do

        kmn=km-ktor
        rkmn=float(kmn)
        ttor=ttor/rkmn

        open (1,file='.\\result\\markers\\marker'//cim//cid//'.tmp',status='unknown') !запись маркеров
          write (1,'(i8)') kmn

          do n=1,km

            if (np(n).ne.0) then
              tr(n)=tr(n)+ttor
              write (1,'(i8,2f12.4,3x,a4)') np(n),tt(n),tr(n),cnp(n),inp(n)
            endif
          enddo
        close (1) !конец записи маркеров

        do n=1,1000
          np(n)=0
          tt(n)=0
        enddo

        open (1,file='.\\result\\markers\\marker'//cim//cid//'.tmp',status='old') !читаем маркеры
          read (1,'(i8)') kmn

          do n=1,kmn
            read (1,'(i8,2f12.4,3x,a4)') np(n),tt(n),tr(n),cnp(n),inp(n)
          enddo

        close (1) !конец чтения маркеров

        km=kmn
      endif                        ! конец условий для сжатий 2,3

      ! Конец блока ввода новых маркеров (разрывы) и исключения старых маркеров (сжатия и торошения)

      continue

      if (id.eq.10.or.id.eq.20.or.id.eq.30) then  ! если конец декады, открываем файл и записываем промежуточный результат

        if(id.eq.10) cdec='1'
        if(id.eq.20) cdec='2'
        if(id.eq.30) cdec='3'

        open (1,file='.\\result\\mon_res\\mon'//cim//'-dec'//cdec//'.res',status='unknown')  ! файл с толщинами, торосами и датами образования маркеров на конец декады

          write (1,'(i8)') km

          do i=1,km
            write (1,'(i8,2f12.4,3x,a4,i6)') np(i),tt(i),tr(i),cnp(i), inp(i)
          enddo

        close (1) ! конец записи в файл с тол., тор. и датами образования маркеров на конец декады

        ! ---модуль расчета гистограммы---
        do j=1,21 ! вычисляем гистограмму распределения толщин, 0.1 см
          gist(j)=0
        enddo

        do i=1,km
          do j=1,201
            rj0=(j-1)*10
            rj1=j*10
            grad0=rj0/100
            grad1=rj1/100
            if(tt(i).ge.grad0.and.tt(i).lt.grad1) gist(j)=gist(j)+1
          enddo
        enddo

        rkm=km
        ! ---модуль расчета гистограммы---

        open (1,file='.\\result\\mon_gst\\mon'//cim//'-dec'//cdec//'.gst',status='unknown') ! файл с гистограммами толщин на конец декады

          do j=1,21
            rj0=(j-1)*10
            rj1=j*10
            grad0=rj0/100
            grad1=rj1/100
            gist(j)=100*gist(j)/rkm
            write (1,'(3f8.2)') grad0,grad1,gist(j)
          enddo

          write (4,'(2x,a2,2x,a1,21f8.2)')  cim,cdec,(gist(j),j=1,21)

        close (1) ! конец записи файла с гистограммами толщин на конец декады

      endif                ! конец записи промежуточного результата
    enddo                  ! конец цикла по дням
  enddo                    ! конец цикла по месяцам

  close (3)
  close (4)

end program Ice_Thickness_Calculation

subroutine DataFiles(temper_a_1, temper_a_2)

  implicit none

  real :: ta(150,40)
  integer :: m, kd, n, j, i
  character cmon*2, cgod*4, emp*1, fnam*60
  character*80 temper_a_1, temper_a_2

  open (2,file='.\\data\\temper_a_ini\\temper_a_2.txt',status='unknown')
  open (3,file='.\\data\\temper_a_ini\\temper_a_1.txt',status='unknown')

  do m=1,11  ! цикл по месяцам m - месяц
		
    if(m.eq.1) cmon='10'
    if(m.eq.2) cmon='11'
    if(m.eq.3) cmon='12'
    if(m.eq.4) cmon='01'
    if(m.eq.5) cmon='02'
    if(m.eq.6) cmon='03'
    if(m.eq.7) cmon='04'
    if(m.eq.8) cmon='05'
    if(m.eq.9) cmon='06'
    if(m.eq.10) cmon='07'
    if(m.eq.11) cmon='08'

    if(m.le.3) cgod='2022'
    if(m.gt.3) cgod='2023'

    if(m.eq.1) kd=31
    if(m.eq.2) kd=30
    if(m.eq.3) kd=31
    if(m.eq.4) kd=31
    if(m.eq.5) kd=28
    if(m.eq.6) kd=31
    if(m.eq.7) kd=30
    if(m.eq.8) kd=31
    if(m.eq.9) kd=30
    if(m.eq.10) kd=31
    if(m.eq.11) kd=31

    fnam='.\\data\\data_temper_a_ini\\'//cgod//cmon//'01a.tb'  !20221001a.tb

    open (1,file=fnam,status='old')

      do n=1,kd*8  ! в каждом дне по 8 значений или по значению в каждые 3 часа
        read (1,'(a1)') emp

        do j=1,37
          read (1,'(144f8.2)') (ta(i,j),i=1,144)
        enddo

        write (2,'(f8.2)') ta(29,6)
        write (3,'(f8.2)') ta(65,7)

        !do j=1,10
          !write (4, '()') (ta(i,j),i=25,35) !модуль для выбора поля данных
        !enddo

      enddo
    close (1)

  enddo

  close (2)
  close (3)

  return

end subroutine DataFiles

subroutine DeleteOldFile

  implicit none

  !Для Windows удаление старых файлов
  call system('del /f result\markers\*.tmp')
  call system('del /f result\mon_gst\*.gst')
  call system('del /f result\mon_res\*.res')
  call system('del /f result\press_type\*.txt')
  call system('del /f data\temper_a_ini\*.txt')

  !Для Unix удаление старых файлов
  !call system('rm -f result\markers\*.tmp')
  !call system('rm -f result\mon_gst\*.gst')
  !call system('rm -f result\moo_res\*.res')
  !call system('rm -f result\press_type\*.txt')
  !call system('rm -f data\temper_a_ini\*.txt')

  return

end subroutine DeleteOldFile
