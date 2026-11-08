program upwind
    implicit none

    !variables
    integer i,xmax,t,tmax,h,it,j
    real pi,dt,c,uinit,uaft,u,dx
    dimension uinit(301),uaft(301)

    !fixed variable
    dx=0.01
    xmax=301
    pi=3.14159265359
    tmax=0.45

    !free variable
    dt=0.01
    u=0.09

    !create file for t=1s until t=10s
    open(11,file='upwind_1.plt')
    open(12,file='upwind_2.plt')
    open(13,file='upwind_3.plt')
    open(14,file='upwind_4.plt')
    open(15,file='upwind_5.plt')
    open(16,file='upwind_6.plt')
    open(17,file='upwind_7.plt')
    open(18,file='upwind_8.plt')
    open(19,file='upwind_9.plt')
    open(20,file='upwind_10.plt')

    !initial condition
    do i=1,xmax
        uinit(i)=0
        uaft(i)=0
    end do
    
    do i=1,11
        uinit(i)=sin(10*pi*(i-1)*dx)
    end do

    !time accumulation
    do 20 j=1,10
        tmax=100*j
    do 30 it=1,tmax

        !upwind
        do i=2,xmax-1
            c=u*dt/dx
            uaft(i)=uinit(i)-c*(uinit(i)-uinit(i-1))
        end do

        !finalize
        do i=1,xmax
            uinit(i)=uaft(i)
        end do

        30 continue
        
        !write to each created files
        h=j+10
        do i=1,xmax
            write(h,*)i,uaft(i)
        end do

        20 continue
end program upwind
