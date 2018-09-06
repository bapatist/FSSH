Module mod_afssh
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: mass_h=1.007825d0,kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=4184.d0
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! Potential
real*8 A_pot,B_pot,C_pot,D_pot,V_exothermicity
real*8 mass
real*8 gamma,temperature

!! Input/Output
integer n_E0
real*8 E0,log_E0_min,log_E0_max
real*8 x_start,x_left,x_right
real*8 l1,l2,r1,r2

!! Classical
integer nclass
real*8,allocatable :: x(:),v(:),acc(:)
real*8,allocatable :: x_old(:),v_old(:)
complex*16,allocatable :: delr(:,:,:),delp(:,:,:)

!! Quantum
integer nquant,state
real*8,allocatable :: si_adiab(:,:),V_k(:),d_ij(:,:,:),vdotd(:,:)
real*8,allocatable :: pot(:,:),force(:,:,:)
complex*16,allocatable :: ci(:)
real*8,allocatable :: si_adiab_old(:,:),ci_old(:)
complex*16,allocatable :: mat(:,:)

!! Evolution
integer n_traj,nsteps,nstep_write,traj_num
real*8 dtc,total_time,curr_time,prob_tun
real*8 energy,pot_en
real*8 ensq_avg,en_avg
integer flag_hop,flag_collapse,flag_terminate,flag_friction,flag_write,flag_frust

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
integer,allocatable :: seed(:)
integer nold
real*8 tim_tot,tim_diag,tim_cl,tim_wr_out,tim_evolve_qm,tim_evolve,tim_coll,tim_T_jk,tim_mat
integer clock_rate,clock_max

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  integer i,size_seed,seed2(2)
  real*8 rnd,c_0,c_e

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar
  nold=0

  open(10,file="AFSSH.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) n_E0
  read(10,*) log_E0_min
  read(10,*) log_E0_max
  read(10,*) N_traj
  read(10,*) dtc
  read(10,*) total_time
  read(10,*) nstep_write
  read(10,*) flag_write
  read(10,*) flag_frust
  read(10,*) flag_hop
  read(10,*) flag_collapse
  read(10,*) flag_friction
  read(10,*) flag_terminate
  read(10,*) nclass
  read(10,*) nquant
  read(10,*) mass
  read(10,*) A_pot
  read(10,*) B_pot
  read(10,*) C_pot
  read(10,*) D_pot
  read(10,*) V_exothermicity
  read(10,*) x_start
  read(10,*) x_left
  read(10,*) x_right
  read(10,*) seed2
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif
  !---------------------------------------------------------- 

  mass=mass*au2kg
  A_pot=A_pot*au2J
  B_pot=B_pot/au2m**2
  C_pot=C_pot*au2J
  D_pot=D_pot/au2m**2
  V_exothermicity=V_exothermicity*au2J

  nsteps=nint(total_time/dtc)+1


  !-----------------------------------------------------------------  
  allocate(x(nclass),v(nclass),acc(nclass))
  allocate(x_old(nclass),v_old(nclass))
  allocate(delr(nquant,nquant,nclass),delp(nquant,nquant,nclass))
  allocate(si_adiab(nquant,nquant),ci(nquant),V_k(nquant))
  allocate(pot(nquant,nquant),force(nquant,nquant,nclass))
  allocate(mat(nquant,nquant))
  allocate(d_ij(nquant,nquant,nclass),vdotd(nquant,nquant))
  allocate(si_adiab_old(nquant,nquant),ci_old(nquant))
  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)
  enddo
  call random_seed(put=seed)
  !-----------------------------------------------------------------  

  if(iflow==2) then
    open(10,file="ifolder.inp")
    read(10,*) ifolder
    close(10)
    N_traj=N_traj/iparallel
    if(ifolder>1) then
      do i=1,(ifolder-1)*N_traj
        seed=seed+1
        call init_cond
      enddo
      call random_seed(put=seed)
    endif
  else
    ifolder=1
  endif

  tim_tot=0.d0
  tim_diag=0.d0
  tim_cl=0.d0
  tim_wr_out=0.d0
  tim_evolve_qm=0.d0
  tim_evolve=0.d0
  tim_coll=0.d0
  tim_T_jk=0.d0
  tim_mat=0.d0

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer i,j,ieta
  integer t1,t2
  real*8 E0_min,E0_max

  call files(0)

  !call cpu_time(t1)
  call system_clock ( t1, clock_rate, clock_max )

  if(iflow==1.or.iflow==2) then
    do j=1,n_E0
      E0_min=dexp(log_E0_min)
      E0_max=dexp(log_E0_max)
      if(n_E0>1) E0=E0_min+(E0_max-E0_min)*(j-1)/dfloat(n_E0-1)
      if(n_E0==1) E0=E0_min
      E0=E0*au2J

      call initialize_averages

      do i=1,N_traj
        call init_cond
        call setup_evolve
        call evolve(nsteps)
        call average_end
      enddo
      call write_average
    enddo

    !call cpu_time(t2);tim_tot=tim_tot+t2-t1
    call system_clock ( t2, clock_rate, clock_max )
    tim_tot=tim_tot+real(t2-t1)/real(clock_rate)

    call files(1)
  endif

end subroutine main
!---------------------------------------------------------- 

subroutine files(flag)
  implicit none
  integer,intent(in)::flag

  if(flag==0) then
    open(10,file="output")
    open(11,file="output_cl")
    open(12,file="output_qm")
    open(13,file="output_dec")

    open(100,file="transmission.out")
  else
    write(10,*)
    write(10,*)"Total time=",tim_tot
    write(10,*)"Diag time=",tim_diag
    write(10,*)"Classical time=",tim_cl
    write(10,*)"Write Output time=",tim_wr_out
    write(10,*)"Evolve Quantum time=",tim_evolve_qm
    write(10,*)"Evolve time=",tim_evolve
    write(10,*)"Collapse time=",tim_coll
    write(10,*)"T_jk compute time=",tim_T_jk
    write(10,*)"mat compute time=",tim_mat
    close(10);close(11);close(12);close(13)
    close(100);close(101)
  endif

end subroutine files
!-----------------------------------------------------------------  

subroutine initialize_averages
  implicit none

  l1=0.d0;l2=0.d0
  r1=0.d0;r2=0.d0

end subroutine initialize_averages
!-----------------------------------------------------------------  

subroutine init_cond
  implicit none
  real*8 energy_0

  energy_0=E0

  x(1)=x_start
  state=1
  call evaluate_variables(0)
  v(1)=dsqrt(2/mass*(energy_0-V_k(1)))

  call evaluate_variables(0)
  call evaluate_variables(1)

  ci=0.d0;ci(1)=1.d0

end subroutine init_cond
!-----------------------------------------------------------------  

subroutine setup_evolve
  implicit none
  integer i,j
  real*8 rnd,sig_x,sig_p,su
  real*8 fac,beta,P_LZ,sgn

  delr=0.d0
  delp=0.d0

  call evaluate_variables(0)
  call evaluate_variables(1)
  call compute_mat_diab

  curr_time=0.d0
  en_avg=0.d0
  ensq_avg=0.d0

end subroutine setup_evolve
!-----------------------------------------------------------------  

subroutine evolve(nsteps)
  implicit none
  integer,intent(in) :: nsteps
  integer i,iterm,flag_coll
  real*8 t1,t2

  call cpu_time(t1)

  do i=1,nsteps
    call write_output(i,0)
    call average(i)
    call save_old_state
    call evolve_quantum
!    call evolve_moments(dtc)
    if(flag_hop==1)call hop
    if(flag_collapse==1)call collapse(dtc,flag_coll)
    if(flag_terminate==1) call traj_terminate(iterm)
      if(iterm==1) exit

    curr_time=curr_time+dtc
  enddo
  call write_output(1,1)

  call cpu_time(t2)
  tim_evolve=tim_evolve+t2-t1

end subroutine evolve
!-----------------------------------------------------------------  

subroutine average(i)
  implicit none
  integer,intent(in) :: i
  complex*16 ci_diab(nquant)
  integer j

  if(flag_write==1) then
    en_avg=en_avg+energy
    ensq_avg=ensq_avg+energy*energy
  endif

end subroutine average
!-----------------------------------------------------------------  

subroutine average_end
  implicit none

  if(x(1)>0.d0.and.state==1)r1=r1+1
  if(x(1)>0.d0.and.state==2)r2=r2+1
  if(x(1)<0.d0.and.state==1)l1=l1+1
  if(x(1)<0.d0.and.state==2)l2=l2+1

end subroutine average_end
!-----------------------------------------------------------------  

subroutine save_old_state
  implicit none

  x_old=x
  v_old=v
  ci_old=ci
  si_adiab_old=si_adiab

end subroutine save_old_state
!-----------------------------------------------------------------  

subroutine revert_state
  implicit none

  x=x_old
  v=v_old
  ci=ci_old
  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine revert_state
!-----------------------------------------------------------------  

subroutine evolve_quantum
  !! 4th order runge kutta
  implicit none
  complex*16,dimension(1:nquant):: k1,k2,k3,k4,ci_ex,ci_diab
  complex*16,dimension(2,nquant,nquant,nclass):: kd1,kd2,kd3,kd4,vec
  real*8 t1,t2

  call cpu_time(t1)

  ci_diab=matmul(si_adiab,ci)
  vec(1,:,:,:)=delr
  vec(2,:,:,:)=delp

  k1=matmul(mat,ci_diab)
  call compute_T_jk(kd1,vec)

  call evolve_classical(dtc/2.d0)
  call compute_mat_diab

  k2=matmul(mat,ci_diab+0.5*dtc*k1)
  k3=matmul(mat,ci_diab+0.5*dtc*k2)

  call compute_T_jk(kd2,vec+0.5*dtc*kd1)
  call compute_T_jk(kd3,vec+0.5*dtc*kd2)

  call evolve_classical(dtc/2.d0)
  call compute_mat_diab

  k4=matmul(mat,ci_diab+dtc*k3)
  call compute_T_jk(kd4,vec+dtc*kd3)

  ci_diab=ci_diab+dtc/6.d0*(k1+2*k2+2*k3+k4)

  vec=vec+dtc/6.d0*(kd1+2*kd2+2*kd3+kd4)

  ci=matmul(transpose(si_adiab),ci_diab)
  delr=vec(1,:,:,:)
  delp=vec(2,:,:,:)

  call cpu_time(t2)
  tim_evolve_qm=tim_evolve_qm+t2-t1

end subroutine evolve_quantum
!-----------------------------------------------------------------  

subroutine evolve_classical(dt)
  !! Velocity Verlet
  implicit none
  real*8,intent(in) :: dt
  real*8 acc_old(nclass)
  real*8 gama_dt,c0,c1,c2
  real*8 delta_r(nclass),delta_v(nclass)
  real*8 t1,t2

  if(flag_friction==0) then
    x=x+v*dt+0.5*acc*dt*dt
    acc_old=acc
    call evaluate_variables(0)
    v=v+0.5*dt*(acc_old+acc)
    call evaluate_variables(1)
  endif

  if(flag_friction==1) then
    gama_dt=gamma*dt
    c0=dexp(-gama_dt)
    c1=1.d0/gama_dt*(1.d0-c0)
    c2=1.d0/gama_dt*(1.d0-c1)

    call stochastic_force(delta_r,delta_v,dt)
    x=x+c1*dt*v+c2*dt*dt*acc+delta_r
    acc_old=acc
    call evaluate_variables(0)
    v=c0*v+(c1-c2)*dt*acc_old+c2*dt*acc+delta_v
    call evaluate_variables(1)

  endif


end subroutine evolve_classical
!-----------------------------------------------------------------  

subroutine evolve_moments(dt)
  !! Eq. 14-17 of JCP 137, 22A513
  implicit none
  real*8,intent(in)::dt
  integer j
  complex*16 mat(nquant,nquant,nclass),T_jk(nquant,nquant,nclass)

  call compute_T_jk_R(T_jk,delr,delp)
  mat=T_jk
  do j=1,nquant
    mat(j,j,:)=mat(j,j,:)-T_jk(state,state,:)
  enddo
  delr=delr+mat*dt

  call compute_T_jk_P(T_jk,delr,delp)
  mat=T_jk
  do j=1,nquant
    mat(j,j,:)=mat(j,j,:)-T_jk(state,state,:)
  enddo
  delp=delp+mat*dt

end subroutine evolve_moments
!-----------------------------------------------------------------  

subroutine compute_T_jk(T_jk,vec)
  implicit none
  complex*16,intent(in):: vec(2,nquant,nquant,nclass)
  complex*16,intent(out):: T_jk(2,nquant,nquant,nclass)
  complex*16 delr(nquant,nquant,nclass),delp(nquant,nquant,nclass)
  complex*16 Tr(nquant,nquant,nclass)
  complex*16 tmp1(nclass),tmp2(nclass)
  integer i
  real*8 t1,t2

  call cpu_time(t1)

  delr=vec(1,:,:,:)
  delp=vec(2,:,:,:)

  call compute_T_jk_R(Tr,delr,delp)
  T_jk(1,:,:,:)=Tr
  call compute_T_jk_P(Tr,delr,delp)
  T_jk(2,:,:,:)=Tr

  tmp1=T_jk(1,state,state,:)
  tmp2=T_jk(2,state,state,:)

  do i=1,nquant
    T_jk(1,i,i,:)=T_jk(1,i,i,:)-tmp1
    T_jk(2,i,i,:)=T_jk(2,i,i,:)-tmp2
  enddo

  call cpu_time(t2)
  tim_T_jk=tim_T_jk+t2-t1

end subroutine compute_T_jk
!-----------------------------------------------------------------  

subroutine compute_T_jk_R(T_jk,delr,delp)
  !! Eq. 14 of JCP 137, 22A513
  implicit none
  complex*16,intent(in) :: delr(nquant,nquant,nclass),delp(nquant,nquant,nclass)
  complex*16,intent(out) :: T_jk(nquant,nquant,nclass)
  integer i,j,k

  do i=1,nclass
    T_jk(:,:,i)=-iota/hbar*commute(pot,delr(:,:,i))+delp(:,:,i)/mass
    do j=1,nclass
      T_jk(:,:,i)=T_jk(:,:,i)-v(j)*commute(d_ij(:,:,j),delr(:,:,i))
    enddo
  enddo

end subroutine compute_T_jk_R
!-----------------------------------------------------------------  

subroutine compute_T_jk_P(T_jk,delr,delp)
  !! Eq. 16 of JCP 137, 22A513
  implicit none
  complex*16,intent(in) :: delr(nquant,nquant,nclass),delp(nquant,nquant,nclass)
  complex*16,intent(out) :: T_jk(nquant,nquant,nclass)
  complex*16 sigma(nquant,nquant)
  real*8 delF(nquant,nquant,nclass)
  integer i,j

  do i=1,nquant
    do j=1,nquant
      sigma(i,j)=ci(i)*dconjg(ci(j))
    enddo
  enddo

  delF=force
  do i=1,nquant
    delF(i,i,:)=delF(i,i,:)-force(state,state,:)
  enddo

  do i=1,nclass
    T_jk(:,:,i)=-iota/hbar*commute(pot,delp(:,:,i))
    T_jk(:,:,i)=T_jk(:,:,i)+0.5*(matmul(delF(:,:,i),sigma)+matmul(sigma,delF(:,:,i)))
    do j=1,nclass
!      T_jk(:,:,i)=T_jk(:,:,i)-v(j)*commute(d_ij(:,:,j),delp(:,:,i))
      T_jk(:,:,i)=T_jk(:,:,i)-v(j)*commute(d_ij(:,:,j),delp(:,:,i))
    enddo
  enddo

end subroutine compute_T_jk_P
!-----------------------------------------------------------------  

subroutine traj_terminate(iterm)
  implicit none
  integer,intent(out) :: iterm

  iterm=0
  if(x(1)>x_right.or.x(1)<x_left)iterm=1

end subroutine traj_terminate
!-----------------------------------------------------------------  

subroutine compute_mat_diab
  implicit none
  integer i,j
  real*8 t1,t2

  call cpu_time(t1)

  mat=0.d0
  do i=1,nquant
    do j=1,nquant
      mat(i,j)=-iota/hbar*sum(si_adiab(i,:)*si_adiab(j,:)*V_k(1:nquant))
    enddo
  enddo

  call cpu_time(t2)
  tim_mat=tim_mat+t2-t1

end subroutine compute_mat_diab
!-----------------------------------------------------------------  

subroutine hop
  implicit none
  integer j,state_tentative,ifrust
  real*8 g_kj,b_jk,a_kk
  real*8 rnd
  complex*16 a_jk

  a_kk=cdabs(ci(state))**2
  call random_number(rnd)
  do j=1,nquant
    if(j.ne.state) then
      a_jk=ci(j)*dconjg(ci(state))
      b_jk=-2*dble(dconjg(a_jk)*vdotd(j,state))
      g_kj=b_jk*dtc/a_kk
      if(g_kj<0.d0)g_kj=0.d0
      prob_tun=g_kj
      if(rnd<g_kj) then
        state_tentative=j
!        write(20,'(3f15.7$)') curr_time*1.d15,0.5*mass*v(1)**2/wave_to_J,energy/wave_to_J
        call velocity_adjust(state_tentative,ifrust)
!        write(20,'(f15.7,i5)') energy/wave_to_J,ifrust
        if(ifrust==0) exit
!        if(ifrust==1) write(21,*)curr_time*1.d15,0.5*mass*v(1)**2/wave_to_J,(V_k(2)-V_k(1))/wave_to_J
        if(ifrust==1)then
          if(flag_frust==0)then
            if(force(1,1,1)*force(2,2,1)<0.d0.and.v(1)*force(2,2,1)<0.d0) v(1)=-v(1)
          endif
        endif
      endif
    endif
  enddo

end subroutine hop
!-----------------------------------------------------------------  

subroutine velocity_adjust(state_tentative,ifrust)
  implicit none
  integer,intent(in)::state_tentative
  integer,intent(out)::ifrust
  real*8 gama,aa,bb,cc,discr,dir(nclass)
  integer i

  aa=0.5d0/mass*(sum(d_ij(state,state_tentative,:)**2))
  bb=vdotd(state,state_tentative)
  cc=V_k(state)-V_k(state_tentative)

  dir=d_ij(state,state_tentative,:)

  discr=bb**2+4*aa*cc
  if(discr<0.d0) then
    ifrust=1
  else
    ifrust=0
    if(bb>=0.d0) gama=(bb-dsqrt(discr))/(2*aa)
    if(bb<0.d0)  gama=(bb+dsqrt(discr))/(2*aa)
    v=v-gama*d_ij(state,state_tentative,:)/mass
    state=state_tentative

    call evaluate_variables(0)
    call evaluate_variables(1)

    delr=0.d0
    delp=0.d0
    !dir=1.d0
    !call adjust_delp(dir)

  endif

end subroutine velocity_adjust
!-----------------------------------------------------------------  

subroutine adjust_delp(dir)
  implicit none
  real*8,intent(in) :: dir(nclass)
  integer i,j,k
  real*8 aa,bb,cc
  real*8 sig,discr,root1,root2

  delp=0.d0

  do k=1,nquant
    sig=cdabs(ci(k))**2
    aa=sum(dir*dir)
    bb=2*mass*sum(v*dir)*sig
    cc=sig**2*(mass**2*sum(v*v)-(2*mass*(energy-V_k(k))))

    discr=bb**2-4*aa*cc

    if(discr<0.d0) then
      delp(k,k,:)=0.d0
    else
      root1=(-bb+dsqrt(discr))/(2*aa)
      root2=(-bb-dsqrt(discr))/(2*aa)
      if(dabs(root1)<dabs(root2))delp(k,k,:)=root1*dir
      if(dabs(root1)>=dabs(root2))delp(k,k,:)=root2*dir
    endif
  enddo

end subroutine adjust_delp
!-----------------------------------------------------------------  

subroutine collapse(dt,iflag_coll)
  implicit none
  real*8,intent(in) :: dt
  integer,intent(out) :: iflag_coll
  real*8 rnd,gama_collapse,gama_reset
  integer n,i,j

  i=state

  do n=1,nquant
    if(n.ne.state) then
      gama_reset=sum((force(n,n,:)-force(i,i,:))*dble(delr(n,n,:)))/(2*hbar)
      gama_collapse=gama_reset-2/hbar*cdabs(sum(force(i,n,:)*delr(n,n,:)))
      gama_collapse=gama_collapse*dt
      gama_reset=-gama_reset*dt
      call random_number(rnd)

      iflag_coll=0
      if(rnd<gama_collapse) then
        iflag_coll=1
        if(flag_collapse==1) then
          do j=1,nquant
            if(j.ne.n) ci(j)=ci(j)/dsqrt(1-cdabs(ci(n)**2))
          enddo
          ci(n)=0.d0
        endif
      endif
      if(rnd<gama_collapse.or.rnd<gama_reset) then
        if(flag_collapse==1) then
          do j=1,nquant
            delr(j,n,:)=0.d0;delr(n,j,:)=0.d0
            delp(j,n,:)=0.d0;delp(n,j,:)=0.d0
          enddo
        endif
      endif
    endif
  enddo

end subroutine collapse
!-----------------------------------------------------------------  

subroutine write_output(n,nflag)
  implicit none
  integer,intent(in)::nflag,n
  integer i
  real*8 t1,t2

  call cpu_time(t1)

  if(nflag==0) then
    if(flag_write==1) then
      if(mod(n,nstep_write)==1.or.nstep_write==1) then
        write(10,'(3f17.5,i5)')curr_time*1.d15,energy/wave_to_J,sum(cdabs(ci)**2),state
        write(11,'(f15.5$)')curr_time*1.d15
        write(12,'(5f15.5)')curr_time*1.d15,cdabs(ci)**2,cdabs(matmul(si_adiab,ci))**2
        write(13,'(9es15.5)')curr_time*1.d15,1.d10*delr(1,1,1),delr(2,2,1)*1.d10,delp(1,1,1)/mass,delp(2,2,1)/mass
        do i=1,1!nclass
          write(11,'(2f15.5$)')x(i)*1.d10,v(i)
        enddo
        write(11,*)!;write(12,*)
      endif
    endif
  endif

  if(nflag==1) then
    if(flag_write==1) then
      write(10,*)"traj num=",traj_num
      write(10,*)"standard deviation=",dsqrt((ensq_avg-en_avg**2/dfloat(nsteps))/dfloat(nsteps))/wave_to_J
      write(10,*)"ci**2=",sum(cdabs(ci)**2)
      write(10,*);write(10,*)
      write(11,*);write(11,*)
      write(12,*);write(12,*)
      write(13,*);write(13,*)
    endif
  endif

  call cpu_time(t2)
  tim_wr_out=tim_wr_out+t2-t1

end subroutine write_output
!-----------------------------------------------------------------  

subroutine write_average
  implicit none
  integer i
  real*8 nf

  nf=dfloat(n_traj)

  l1=l1/nf;l2=l2/nf
  r1=r1/nf;r2=r2/nf

  write(100,'(8f15.7)') dlog(E0/au2J),l1,l2,r1,r2

end subroutine write_average
!-----------------------------------------------------------------  

subroutine evaluate_variables(flag)
  implicit none
  integer,intent(in):: flag
  integer i,j

  if(flag==0) then
    call tise
  endif

  if(flag==1) then
    energy=pot_en+0.5*mass*sum(v*v)

    !! CHECK XXX
    vdotd=0.d0
    vdotd(1,2)=v(1)*d_ij(1,2,1)
    vdotd(2,1)=-vdotd(1,2)
    !! CHECK XXX
  endif

end subroutine evaluate_variables
!-----------------------------------------------------------------  

subroutine tise
  !! Output - pot_en,acc
  !! Output - V_k,d_ij
  implicit none
  integer i,j,k
  real*8 Hamil(nquant,nquant),ens(nquant),vect(nquant,nquant)
  real*8 delH_dels(nquant,nquant),delH_dels_ad(nquant,nquant)
  real*8 pot_cl,acc_cl(nclass),acc_qm(nclass)
  real*8 t1,t2

  call cpu_time(t1)

  call compute_potential(Hamil,delH_dels)
  call diag_2(Hamil,nquant,ens,vect,nquant,nold)

  si_adiab_old=si_adiab
  do i=1,nquant
    si_adiab(:,i)=vect(:,i)
    if(sum(si_adiab(:,i)*si_adiab_old(:,i))<0.d0)si_adiab(:,i)=-si_adiab(:,i)
  enddo

  V_k=ens
  delH_dels_ad=matmul(transpose(si_adiab),matmul(delH_dels,si_adiab))

  !! CHECK XXX !!
  d_ij=0.d0
  d_ij(1,2,1)=delH_dels_ad(1,2)/(V_k(2)-V_k(1))
  d_ij(2,1,1)=-d_ij(1,2,1)
  !! CHECK XXX !!

  call cpu_time(t2);tim_diag=tim_diag+(t2-t1)

  call cpu_time(t1)
  call potential_classical(pot_cl,acc_cl)
  acc_qm=0.d0
  acc_qm(1)=-1.d0/mass*delH_dels_ad(state,state)
  pot_en=pot_cl+ens(state)
  acc=acc_cl+acc_qm

  pot=0.d0
  force=0.d0
  do i=1,nquant
    pot(i,i)=pot_cl+ens(i)
    force(i,i,1)=mass*acc_cl(1)-delH_dels_ad(i,i)
    do j=1,nquant
      if(j.ne.i)force(i,j,:)=d_ij(i,j,:)*(V_k(i)-V_k(j))
    enddo
  enddo

  call cpu_time(t2);tim_cl=tim_cl+(t2-t1)

end subroutine tise
!-----------------------------------------------------------------  

subroutine potential_classical(pot_cl,acc_cl)
  implicit none
  real*8,intent(out) :: pot_cl,acc_cl(nclass)
  integer i
  real*8 tmp,mw2

  pot_cl=0.d0
  acc_cl=0.d0

end subroutine potential_classical
!-----------------------------------------------------------------  

subroutine compute_potential(H_diab,delV_dels)
  implicit none
  real*8,intent(out) :: H_diab(2,2),delV_dels(2,2)
  real*8 H1,H2,H12,dH1_ds,dH2_ds,dH12_ds
  real*8 tmp

  H1=0.d0
  dH1_ds=0.d0

  H2=-A_pot*dexp(-B_pot*x(1)**2) + V_exothermicity
  dH2_ds=2*A_pot*B_pot*x(1)*dexp(-B_pot*x(1)**2)

  H12=C_pot*dexp(-D_pot*x(1)**2)
  dH12_ds=-2*D_pot*x(1)*H12

  H_diab(1,1)=H1;  H_diab(1,2)=H12
  H_diab(2,1)=H12; H_diab(2,2)=H2

  delV_dels(1,1)=dH1_ds ; delV_dels(1,2)=dH12_ds
  delV_dels(2,1)=dH12_ds; delV_dels(2,2)=dH2_ds

end subroutine compute_potential
!-----------------------------------------------------------------  

subroutine check_acceleration
  implicit none
  integer i,nflag
  real*8 delx,en_old,acc_sav(nclass)
  real*8 q0,rnd

  delx=1.d-17
  state=1

  do i=1,nclass
    call random_number(rnd)
    x(i)=(rnd*2-1.d0)*1.d-10
  enddo

  call evaluate_variables(0)
  en_old=pot_en;acc_sav=acc

  write(6,*) "delx=",delx
  write(6,*)

  do i=1,nclass
    x(i)=x(i)+delx
    call evaluate_variables(0)
    acc(i)=-(pot_en-en_old)/delx/mass
    write(6,*)"Analytical acceleration =",acc_sav(i)
    write(6,*)"Numerical acceleration  =",acc(i)
    write(6,*)"Error =",(acc(i)-acc_sav(i))/acc(i)*100.d0
    write(6,*)
    x(i)=x(i)-delx
  enddo

  stop

end subroutine check_acceleration
!---------------------------------------------------------- 

function commute(A,B)
  real*8,intent(in) :: A(:,:)
  complex*16,intent(in) :: B(:,:)
  complex*16 commute(nquant,nquant)

  commute=matmul(A,B)-matmul(B,A)

end function commute
!-----------------------------------------------------------------  

subroutine stochastic_force(delr,delv,dt)
  implicit none
  real*8,intent(in)::dt
  real*8,intent(out) :: delr(nclass),delv(nclass)!f(nclass)
  integer i
  real*8 rnd1,rnd2,sig_r,sig_v,sig_rv,gdt

  gdt=gamma*dt
  sig_r=dt*dsqrt(kb*temperature/mass *1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
  sig_v=dsqrt(kb*temperature/mass*(1-dexp(-2*gdt)))
  sig_rv=(dt*kb*temperature/mass* 1.d0/gdt *(1-dexp(-gdt))**2)/(sig_r*sig_v)  !! correlation coeffecient

  do i=1,nclass
    call gaussian_random_number(rnd1)
    call gaussian_random_number(rnd2)
    delr(i)=sig_r*rnd1
    delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
  enddo

!  delr=delr-sum(delr)/dfloat(nclass)
!  delv=delv-sum(delv)/dfloat(nclass)

end subroutine stochastic_force
!-----------------------------------------------------------------  

subroutine draw_pes
  implicit none
  integer i
  real*8 H(2,2),dH(2,2)

  state=1
  open(10,file="pes.out")
  do i=1,5000
    x(1)=-2.d-10+4.d-10*i/5000.d0
    call compute_potential(H,dH)
    call evaluate_variables(0)
!    write(10,'(4f15.7)')x(1)/au2m,H(1,1)/au2J,H(2,2)/au2J,H(1,2)/au2J
    write(10,'(4f15.7)')x(1)*1.d10,V_k/wave_to_J,H(1,2)/wave_to_J
  enddo
  close(10)
  stop

end subroutine draw_pes
!-----------------------------------------------------------------  

subroutine diag_2(Hamil,nquant,ens,vect,m_values,nold)
  implicit none
  integer,intent(in) :: nquant,nold,m_values
  real*8,intent(in) :: Hamil(nquant,nquant)
  real*8,intent(out) :: ens(nquant)
  real*8,intent(out) :: vect(nquant,nquant)
  real*8 e1,e2,a,b,c,eps,c1,c2

  a=Hamil(1,1)
  c=Hamil(1,1)
  b=Hamil(1,2)

  eps=(c-a)/2.d0

  e1=(a+c)/2.d0 - dsqrt(eps**2+b**2)
  e2=(a+c)/2.d0 + dsqrt(eps**2+b**2)

  c1=b/dsqrt((e1+eps)**2+b**2)
  c2=(e1+eps)/dsqrt((e1+eps)**2+b**2)
  vect(1,1)=c1
  vect(2,1)=c2
  vect(1,2)=-c2
  vect(2,2)=c1

end subroutine diag_2
!-----------------------------------------------------------------  

subroutine gaussian_random_number(rnd)
  !! generates gaussian distribution with center 0, sigma 1
  !! q0+sig*rnd gives center=q0, sigma=sig
  implicit none
  real*8,intent(out)::rnd
  real*8 rnd1,rnd2,pi

  pi=dacos(-1.d0)

  call random_number(rnd1)
  call random_number(rnd2)
  rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!----------------------------------------------------------
End Module mod_AFSSH
