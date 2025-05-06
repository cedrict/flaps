!==================================================================================================!
!==================================================================================================!
!                                                                                                  !
! FLAPS                                                                                            !
!                                                                                                  !
!==================================================================================================!
!==================================================================================================!

subroutine compute_dNdx_dNdz_at_qpt(iel,iq,dNdx,dNdz,jcob)

use module_core
use module_constants,only: eps
use module_mesh
use module_arrays 

implicit none

integer(int32),intent(in):: iq,iel
real(real64),intent(out):: dNdx(mV),dNdz(mV)
real(real64),intent(out):: jcob 

real(real64):: jcb(2,2),jcbi(2,2)

!==================================================================================================!
!==================================================================================================!
!@@ \subsection{compute\_dNdx\_dNdz\_at\_qpt.f90}
!@@ This subroutine computes $\partial{\bN}/\partial \{x,z\}$, 
!@@ at a quadrature point iq (between 1 and nqel)) passed as argument.
!@@ It also returns the determinant jcob of the Jacobian matrix.
!==================================================================================================!

jcb(1,1)=sum(dNVqdr(1:mV,iq)*mesh(iel)%xV(1:mV))
jcb(1,2)=sum(dNVqdr(1:mV,iq)*mesh(iel)%zV(1:mV))
jcb(2,1)=sum(dNVqds(1:mV,iq)*mesh(iel)%xV(1:mV))
jcb(2,2)=sum(dNVqds(1:mV,iq)*mesh(iel)%zV(1:mV))

jcob=jcb(1,1)*jcb(2,2)-jcb(1,2)*jcb(2,1)

if (jcob<eps) stop 'jcob=0'

jcbi(1,1)=    jcb(2,2) /jcob
jcbi(1,2)=  - jcb(1,2) /jcob
jcbi(2,1)=  - jcb(2,1) /jcob
jcbi(2,2)=    jcb(1,1) /jcob

dNdx(:)=jcbi(1,1)*dNVqdr(:,iq)+jcbi(1,2)*dNVqds(:,iq)
dNdz(:)=jcbi(2,1)*dNVqdr(:,iq)+jcbi(2,2)*dNVqds(:,iq)

end subroutine

!==================================================================================================!
!==================================================================================================!
