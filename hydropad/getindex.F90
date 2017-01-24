SUBROUTINE getindex(i,ch1,ch2,ch3)
!
        INTEGER :: n1,n2,n3,n2aux
        INTEGER, INTENT(IN) :: i
        CHARACTER*1, INTENT(OUT) :: ch1,ch2,ch3

        n1=i/100
        n2aux=i-n1*100
        n2=n2aux/10
        n3=n2aux-n2*10

        n1=48+n1
        n2=48+n2
        n3=48+n3

        ch1=char(n1)
        ch2=char(n2)
        ch3=char(n3)

END SUBROUTINE getindex

SUBROUTINE getindex10000(i,ch1,ch2,ch3,ch4,ch5)
!
	INTEGER :: n1,n2,n3,n2aux,n4,n5
	INTEGER, INTENT(IN) :: i
	CHARACTER*1, INTENT(OUT) :: ch1,ch2,ch3,ch4,ch5

	n1=i/10000
	n2aux=i-n1*10000
	n2=n2aux/1000
	n2aux=n2aux-n2*1000
	n3=n2aux/100
	n2aux=n2aux-n3*100
	n4=n2aux/10
	n5=n2aux-n4*10

        n1=48+n1
        n2=48+n2
        n3=48+n3
        n4=48+n4
        n5=48+n5

        ch1=char(n1)
        ch2=char(n2)
        ch3=char(n3)
        ch4=char(n4)
        ch5=char(n5)

END SUBROUTINE getindex10000
