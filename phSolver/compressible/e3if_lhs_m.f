      module e3if_lhs_m
c
        use e3if_defs_m
c
        implicit none
c
      contains
c
      subroutine set_lhs_matrices
c
        integer :: n,p,q,r
c
        do q = 1,nflow
          do p = 1,nflow
c
            do n = 1,nshl0
c
              AiNa0(:,1,p,q) = Ai0(:,1,p,q) * shp0(:,n)
              AiNa0(:,2,p,q) = Ai0(:,2,p,q) * shp0(:,n)
              AiNa0(:,3,p,q) = Ai0(:,3,p,q) * shp0(:,n)
c
              KijNaj0(:,1,p,q) = Kij0(:,1,1,p,q)*shg0(:,n,1)
     &                         + Kij0(:,1,2,p,q)*shg0(:,n,2)
     &                         + Kij0(:,1,3,p,q)*shg0(:,n,3)
              KijNaj0(:,2,p,q) = Kij0(:,2,1,p,q)*shg0(:,n,1)
     &                         + Kij0(:,2,2,p,q)*shg0(:,n,2)
     &                         + Kij0(:,2,3,p,q)*shg0(:,n,3)
              KijNaj0(:,3,p,q) = Kij0(:,3,1,p,q)*shg0(:,n,1)
     &                         + Kij0(:,3,2,p,q)*shg0(:,n,2)
     &                         + Kij0(:,3,3,p,q)*shg0(:,n,3)
c
            enddo
c
            do n = 1,nshl1
c
              AiNa1(:,1,p,q) = Ai1(:,1,p,q) * shp1(:,n)
              AiNa1(:,2,p,q) = Ai1(:,2,p,q) * shp1(:,n)
              AiNa1(:,3,p,q) = Ai1(:,3,p,q) * shp1(:,n)
c
              KijNaj1(:,1,p,q) = Kij1(:,1,1,p,q)*shg1(:,n,1)
     &                         + Kij1(:,1,2,p,q)*shg1(:,n,2)
     &                         + Kij1(:,1,3,p,q)*shg1(:,n,3)
              KijNaj1(:,2,p,q) = Kij1(:,2,1,p,q)*shg1(:,n,1)
     &                         + Kij1(:,2,2,p,q)*shg1(:,n,2)
     &                         + Kij1(:,2,3,p,q)*shg1(:,n,3)
              KijNaj1(:,3,p,q) = Kij1(:,3,1,p,q)*shg1(:,n,1)
     &                         + Kij1(:,3,2,p,q)*shg1(:,n,2)
     &                         + Kij1(:,3,3,p,q)*shg1(:,n,3)
c
            enddo
c
          enddo
        enddo
c
        KijNajC0 = zero
        KijNajC1 = zero
c
        do q = 1,nflow
          do p = 1,nflow
            do r = 1,nflow
              KijNajC0(:,1,p,q) = KijNajC0(:,1,p,q) + KijNaj0(:,1,p,r)*cmtrx(:,r,q)
              KijNajC0(:,2,p,q) = KijNajC0(:,2,p,q) + KijNaj0(:,2,p,r)*cmtrx(:,r,q)
              KijNajC0(:,3,p,q) = KijNajC0(:,3,p,q) + KijNaj0(:,3,p,r)*cmtrx(:,r,q)
              KijNajC1(:,1,p,q) = KijNajC1(:,1,p,q) + KijNaj1(:,1,p,r)*cmtrx(:,r,q)
              KijNajC1(:,2,p,q) = KijNajC1(:,2,p,q) + KijNaj1(:,2,p,r)*cmtrx(:,r,q)
              KijNajC1(:,3,p,q) = KijNajC1(:,3,p,q) + KijNaj1(:,3,p,r)*cmtrx(:,r,q)
            enddo
          enddo
        enddo
c
      end subroutine set_lhs_matrices
c
      subroutine calc_egmass(egmass,AiNa1,KijNaj0,KijNaj1,KijNajC0,shp0,shp1,n0,n1,WdetJ,nshl0,nshl1)
c
        real*8, dimension(:,:,:), intent(out) :: egmass
        real*8, dimension(:,:,:,:), intent(in) :: AiNa1,KijNaj0,KijNaj1,KijNajC0
        real*8, dimension(:,:), intent(in) :: shp0,shp1,n0,n1
        real*8, dimension(:), intent(in) :: WdetJ
        integer, intent(in) :: nshl0,nshl1
c
        integer :: i,j,p,q,inode0,inode1,isd
c
c... loop through the columns
c
        do inode1 = 1, nshl1
          q = nflow * (inode1 - 1)
c
          do inode0 = 1, nshl0
            p = nflow * (inode0 - 1)
c
            do isd = 1,nsd
              egmass(:,p,q) = egmass(:,p,q) + (
     &              pt50 * shp0(:,inode0) * (AiNa1(:,isd,p,q)+KijNaj1(:,isd,p,q))*n0(:,isd)
c     &            + pt50 * s * KijNajC0(:,isd,p,q)*n1(:,isd)
c     &            + e*mu/h * ctc(:,p,q)*shp0(:,inode0)*n0(:,isd)*shp1(:,inode1)*n1(:,isd)
     &            )*WdetJ
            enddo
c
          enddo
c
        enddo
c
      end subroutine calc_egmass
c
      end module e3if_lhs_m
