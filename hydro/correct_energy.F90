SUBROUTINE correct_energy(eleft, eright, p3dnew

   do i=nmin,nmax
!CLA this introduces differences in parallel results
      eeleft = etotnew(i-1)
      if(i .EQ. nmin) eeleft = 0.0
      ee = etotnew(i)
      eeright = etotnew(i+1)
      if(i .EQ. nmax) eeright = 0.0
      ee2 = max(eeleft, ee, eeright)
      vvv = vxnew(i)*vxnew(i)+vynew(i)*vynew(i)+vznew(i)*vznew(i)
      eeaux = etotnew(i)-0.50*rhonew(i)*vvv
      ee1 = eeaux/ee2
      if (ee1 >= eta2) eintnew(i)=eeaux

      ee0 = eeaux/etotnew(i)
      if (ee0 >= eta1) then
         pnew(i) = (gamma-1)*eeaux
      else
         pnew(i) = (gamma-1)*eintnew(i)
         !write(*,*)"Correcting with internal energy!!!!!!!!!!!!!!!"
      endif
   enddo



END SUBROUTINE correct_energy
