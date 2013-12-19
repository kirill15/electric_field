        PROGRAM TXTTODAT
        
        CHARACTER(12) path
        
        dimension n(4)
        real*8 r, z
        real*8 f

        path = './for_telma/'
        
        m0=0
        m1=1
        
        open(2, file = 'rz.txt')
        open(3, file = path//'rz.dat', access = 'direct', recl = 16)     
        read(2, *) kuzlov
        do i = 1, kuzlov
          read(2, *) r, z
          write(3, rec = i) r, z
        enddo
        close(2)
        close(3)
       
        
        open(2, file = 'nvtr.txt')
        open(3, file = path//'nvtr.dat', access = 'direct', recl = 24)
        read(2, *) ktr
        do i = 1, ktr
          read(2, *) (n(j), j = 1, 4)
*          write(3, rec = i) (n(j),j=1, 4), m0, m1
          write(3, rec = i) n(3), n(4), n(1), n(2), m0, m1
        enddo
        close(2)
        close(3)
       
       
        open(2, file = 'nvkat.txt')
        open(3, file = path//'nvkat2d.dat', access = 'direct', recl = 4)
        do i = 1, ktr
          read(2, *) k
          write(3, rec = i) k
        enddo
        close(2)
        close(3)
        
        
        open(2, file = path//'inf2tr.dat')
        write(2, *) 'islau= 0 indku1= 0 indfro= 1'
        write(2, *) 'kuzlov= ', kuzlov
        write(2, *) 'ktr= ', ktr
        write(2, *) 'kt1= ', 0
        write(2, *) 'kreb2= 0  kreb3= 0'
        write(2, *) 'kisrr1= 2 kisrr2= 2 kisrr3= 2  kbrsr= 8'
        write(2, *) 'kreb4= 0'


        open(2, file = 'v.txt');
        open(3, file = path//'v2.dat', access = 'direct', recl=8)
        do i = 1, kuzlov
          read(2, *) f
          write(3, rec = i) f
        enddo
        close(2)
        close(3)
                
        END PROGRAM
