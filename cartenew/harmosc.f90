module harmosc_solver
    real(kind=8),allocatable :: sol(:,:)

    contains
        subroutine solver(T,alphak,freq,re,im,Q)
            !$ use omp_lib
            implicit none
            integer :: K,n,nfreq,i,j
            real(kind=8) :: Q
            real(kind=8),dimension(:),intent(in) :: T,alphak,freq
            real(kind=8),dimension(:,:),intent(in) :: re,im
            real(kind=8), allocatable, dimension(:) ::  cphase,sphase,A,B
            !complex(kind=8), allocatable, dimension(:,:) :: tfsp
            !complex(kind=8) :: ij=(0.,1.)

            n=size(T)
            K=size(alphak)
            nfreq=size(freq)
            !print*,omp_get_num_procs(),omp_get_num_threads()
            allocate(sol(K,n))!,tfsp(K,nfreq))
            ! do i=1,K
            ! tfsp(i,:)=cmplx(re(i,:),im(i,:))/(alphak(i)**2-freq**2-ij*freq*alphak(i)/Q)
            ! enddo

            sol=0
             !$omp parallel default(shared) private(i,j,cphase,sphase)
            allocate(cphase(nfreq),sphase(nfreq))
             !$omp do schedule(dynamic,3)
             do i=1,n
                 cphase=cos(freq*T(i))
                 sphase=sin(freq*T(i))
                 do j=1,nfreq
                     sol(:,i)=sol(:,i)+re(:,j)*cphase(j)-im(:,j)*sphase(j)
                 enddo
             enddo
             !$omp end do
             deallocate(cphase,sphase)
             !$omp end parallel
            allocate(A(k),B(k))
            A=0
            B=0
            A=-sum(re,dim=2)
            B=sum(im*spread(freq,dim=1,ncopies=K),dim=2)/sqrt(alphak**2-alphak**2/Q/2._8)&
            +alphak/2/Q*A/sqrt(alphak**2-alphak**2/Q/2._8)
            !print*,freq
             !A(i)=-sum(re(i,:))
            !B=sum(im*spread(freq,dim=1,ncopies=K),dim=2)/alphak
            !do i=1,k
            !A(i)=-sum(re(i,:))
            !B=sum(im*spread(freq,dim=1,ncopies=K),dim=2)/alphak
            !A(i)=-sol(i,1)
            !B(i)=alphak(i)/2/Q*A(i)-(-0.5*sol(i,3)+2*sol(i,2)-3./2.*sol(i,1))/(T(2)-T(1))  
            ! A(i)=-sum(((alphak(i)**2-freq**2)*re(i,:)-freq*alphak(i)/Q*im(i,:))/&
            ! ((alphak(i)**2-freq**2)**2+freq**2*alphak(i)**2/Q**2))
            ! B(i)=sum(-(alphak(i)/2._8/Q/sqrt(alphak(i)**2-alphak(i)**2/Q**2/4._8))*&
            ! ((alphak(i)**2-freq**2)*re(i,:)-freq*alphak(i)/Q*im(i,:))/((alphak(i)**2-freq**2)**2+freq**2*alphak(i)**2/Q**2)-&
            ! freq/sqrt(alphak(i)**2-alphak(i)**2/Q**2/4._8)*&
            ! ((alphak(i)**2-freq**2)*im(i,:)-freq*alphak(i)/Q*re(i,:))/((alphak(i)**2-freq**2)**2+freq**2*alphak(i)**2/Q**2))
            !enddo
            do i=1,n
                    sol(:,i)=sol(:,i)+exp(-alphak/2/Q*T(i))*(A*cos(sqrt(alphak**2-alphak**2/Q/2._8)*T(i))&
                    +B*sin(sqrt(alphak**2-alphak**2/Q/2._8)*T(i)))
            enddo
        end subroutine
end module
