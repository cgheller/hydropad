MODULE ccool
!
      IMPLICIT NONE
      SAVE

! curve di emissione che dipendono da rho
      real*8, dimension (:), allocatable :: trtable
      real*8, dimension (:), allocatable :: mrtable
      real*8, dimension (:), allocatable :: qrtable

! curve di emissione che dipendono da rho**2
      real*8, dimension (:), allocatable :: tr2table
      real*8, dimension (:), allocatable :: mr2table
      real*8, dimension (:), allocatable :: qr2table

! curve di emissione che dipendono da rho**3
      real*8, dimension (:), allocatable :: tr3table
      real*8, dimension (:), allocatable :: mr3table
      real*8, dimension (:), allocatable :: qr3table


!
      real*8 :: TT 
      real*8, parameter :: Hf=0.75, Hef=1.0-Hf, ff= 0.5*(1+Hf)
      real*8, parameter :: nb2na=1./(Hf+Hef/4.0)

!     Hf: frazione di barioni che forma atomi di idrogeno (n(H)/ntot)
!     Hef:frazione di barioni che forma atomi di elio
!     ff: numero di elettroni per barione presenti nell'universo 
!     nb2na: numero di barioni su numero di atomi
!
      real*8 ciH,crHp,ciHe,ciHep,crHep,crHepp
      real*8 fH, fHp, fHe, fHep, fHepp
      real*8 oldf, hem, heem

!     coefficienti di ionizzazione e ricombinazione
!      cix = coeff. di ionizzazione collisionale per x
!      crx = coeff. di ricombinazione per x

!     abbondanze relative
!      fH, fHp: idrogeno neutro e ionizzato
!      fHe, fHep ,fHepp: elio neutro, ion. una volta, ion. due volte

      real*8 nh2n,nhe2n,ne2n
      parameter(nh2n=Hf,nhe2n=Hef/4.0)

!     nh2n: numero di atomi di H/numero di barioni
!     nhe2n: numero di atomi di elio/numero di barioni
!     ne2n: numero di elettroni liberi per barione
!     queste quantita' rientrano nel calcolo dei coefficienti di
!     emissione che sono proporzionali alla densita' in numero degli
!     elettroni e alla densita' in numero della specie considerata.
!     Quindi, avendo a disposizione Yxi
!         Yxi= n(xi)/sum(tutti gli stati i della specie x)
!     si ottiene n(xi) come
!         n(xi)+ Yxi* sum(tutti gli stati..)
!     ora, sum(tutti gli stati) deve essere uguale a nh o nhe ovvero
!         nh= nh2n*n
!     dove n e' il numero di barioni:
!         nh= nh2n*(rho/mp)

      real*8 K2eV   
      parameter(K2eV=10000.0)

      real*8, parameter :: minimum=1.0e-50 
      real*8 mu,kbl,sectoy
      parameter(mu=0.56,kbl=1.38e-30,sectoy=3.171e-8)

      real*8 ecom,ebreem,eH,eHe
      real*8 ecih,ecihe,ecihep
      real*8 eceh,eceheI,eceheII
      real*8 erhp,erhep,erhepp,erdhe
      real*8 lambdar, lambdar2, lambdar3


!      integer npoints
!      parameter(npoints=800)
!      real*8 lambda(npoints)
!
END MODULE ccool
