!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

pure subroutine NNN(r,s,N,m,order)

use iso_fortran_env

implicit none

integer(int32),intent(in):: m,order
real(real64),intent(in):: r,s
real(real64),intent(out):: N(m)

real(real64):: Nmr,Nlr,Nrr,Nls,Nms,Nrs,Nlt,Nmt,Nrt

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{NNN}
!@@ Basis functions ${\bN}_i$. Spaces supported: $Q_1$, $Q_2$.
!==================================================================================================!

select case(order)

   case(1)

      N(1)=0.25*(1-r)*(1-s)
      N(2)=0.25*(1+r)*(1-s)
      N(3)=0.25*(1+r)*(1+s)
      N(4)=0.25*(1-r)*(1+s)

   case(2)

      Nlr=0.5d0*r*(r-1d0) ; Nls=0.5d0*s*(s-1d0)
      Nmr=(1d0-r**2)      ; Nms=(1d0-s**2)     
      Nrr=0.5d0*r*(r+1d0) ; Nrs=0.5d0*s*(s+1d0)
      N(1)= Nlr * Nls
      N(5)= Nmr * Nls
      N(2)= Nrr * Nls
      N(8)= Nlr * Nms
      N(9)= Nmr * Nms
      N(6)= Nrr * Nms
      N(4)= Nlr * Nrs
      N(7)= Nmr * Nrs
      N(3)= Nrr * Nrs

end select

end subroutine

!==================================================================================================!
!==================================================================================================!

subroutine dNNNdr(r,s,dNdr,m,order)

use iso_fortran_env

implicit none

integer(int32),intent(in):: m,order
real(real64),intent(in):: r,s
real(real64),intent(out):: dNdr(m)

real(real64) dNmr,dNlr,dNrr,Nls,Nms,Nrs

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{dNNNdr}
!@@ Basis functions derivatives $\partial{\bN}_i/\partial r$. Spaces supported: $Q_1$, $Q_2$.
!==================================================================================================!

select case(order)

   case(1)

      dNdr(1)=-0.25*(1-s)
      dNdr(2)=+0.25*(1-s)
      dNdr(3)=+0.25*(1+s)
      dNdr(4)=-0.25*(1+s)

   case(2)

      dNlr=r-0.5d0 ; Nls=0.5d0*s*(s-1d0)
      dNmr=-2d0*r  ; Nms=(1d0-s**2)     
      dNrr=r+0.5d0 ; Nrs=0.5d0*s*(s+1d0)
      dNdr(1)= dNlr * Nls
      dNdr(5)= dNmr * Nls
      dNdr(2)= dNrr * Nls
      dNdr(8)= dNlr * Nms
      dNdr(9)= dNmr * Nms
      dNdr(6)= dNrr * Nms
      dNdr(4)= dNlr * Nrs
      dNdr(7)= dNmr * Nrs
      dNdr(3)= dNrr * Nrs

end select 

end subroutine

!==================================================================================================!
!==================================================================================================!

pure subroutine dNNNds(r,s,dNds,m,order)

use iso_fortran_env

implicit none

integer(int32),intent(in):: m,order
real(real64),intent(in):: r,s
real(real64),intent(out):: dNds(m)

real(real64):: Nmr,Nlr,Nrr,dNls,dNms,dNrs

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{dNNNds}
!@@ Basis functions derivatives $\partial{\bN}_i/\partial s$. Spaces supported: $Q_1$, $Q_2$.
!==================================================================================================!

select case(order)

   case(1)

      dNds(1)=-0.25*(1-r)
      dNds(2)=-0.25*(1+r)
      dNds(3)=+0.25*(1+r)
      dNds(4)=+0.25*(1-r)

   case(2)

      Nlr=0.5d0*r*(r-1d0) ; dNls=s-0.5d0
      Nmr=(1d0-r**2)      ; dNms=-2d0*s 
      Nrr=0.5d0*r*(r+1d0) ; dNrs=s+0.5d0
      dNds(1)= Nlr * dNls
      dNds(5)= Nmr * dNls
      dNds(2)= Nrr * dNls
      dNds(8)= Nlr * dNms
      dNds(9)= Nmr * dNms
      dNds(6)= Nrr * dNms
      dNds(4)= Nlr * dNrs
      dNds(7)= Nmr * dNrs
      dNds(3)= Nrr * dNrs

end select

end subroutine

!==================================================================================================!
