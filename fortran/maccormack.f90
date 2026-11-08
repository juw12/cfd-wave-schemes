program maccormack
    implicit none

    !variables
    integer i,xmax,t,tmax,h,it,j
    real pi,dt,c,uinit,uaft,up,u,dx
    dimension uinit(301),uaft(301),up(301)

    !fixed variable
    dx=0.01
    xmax=301
    pi=3.14159265359
    tmax=0.45

    !free variable
    dt=0.01
    u=0.09

    !create file for t=1s until t=10s
    open(11,file='maccormack_1.plt')
    open(12,file='maccormack_2.plt')
    open(13,file='maccormack_3.plt')
    open(14,file='maccormack_4.plt')
    open(15,file='maccormack_5.plt')
    open(16,file='maccormack_6.plt')
    open(17,file='maccormack_7.plt')
    open(18,file='maccormack_8.plt')
    open(19,file='maccormack_9.plt')
    open(20,file='maccormack_10.plt')

    !initial condition
    do i=1,xmax
        uinit(i)=0
        uaft(i)=0
        up(i)=0
    end do
    
    do i=1,11
        uinit(i)=sin(10*pi*(i-1)*dx)
    end do

    !time accumulation
    do 20 j=1,10
        tmax=100*j
    do 30 it=1,tmax

        !explicit
        !predictor
        do i=1,xmax-1
            c=u*dt/dx
            up(i)=uinit(i)-c*(uinit(i+1)-uinit(i))
        end do

        !corrector
        do i=2,xmax
            c=u*dt/dx
            uaft(i)=0.5*((uinit(i)+up(i))-c*(up(i)-up(i-1)))
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
end program maccormack