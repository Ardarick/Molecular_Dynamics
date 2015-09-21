
    program prog

    implicit none

    real h,m,sx,sy,sz,r0,a1,a0,p,xi,q,nav,kb,rs,ekin,etot,epot,ebound,er,v,kbt,betta,qu
    real dxt,dyt,dzt,drt,sigma,t,t0,Ron,Roff,obrez,proizv_obrez,pexp,qexp,vscal,sigma1,R,Rf, KineticEnergy
    real x(4000),xt(4000),y(4000),yt(4000),z(4000),zt(4000)
    real vx(4000),vxt(4000),vy(4000),vyt(4000),vz(4000),vzt(4000)
    real ax(4000),axt(4000),ay(4000),ayt(4000),az(4000),azt(4000)
    real Eba(4000),fx(4000),fy(4000),fz(4000)
    real dy(200000),dz(200000),dx(200000),dr(200000),pexp1(200000),qexp1(200000)
    real fr(200000),fb(200000),f(200000),w(200000),dw(200000)
    real Vmax,dVmax,vkv,Vkvmax,dVkvmax,mxx(20000),mxy(20000),mxz(20000),mxv(20000)
    integer i,j,k,natom,np,ip(200000),jp(200000),isrt(4000),nstat,flag,count,iscal,count1
    integer jx,jy,jz,jv

    open(8,file="Result.csv")
    !open(9,file="InputRelax.txt")
    open(9,file="input.txt")
    open(10,file="Output.csv")
    open(11,file="Temperature.csv")
    open(12,file="Velocity.csv")

    natom=2000
    nstat=0
    iscal=100
    
    do i=1,natom
        read(9,*) x(i),y(i),z(i),vx(i),vy(i),vz(i),isrt(i)
    end do

    do i=1,natom
        ax(i)=0
        ay(i)=0
        az(i)=0
    end do
    
    !====================Константы и параметры==============================================================================================
    h=1.0e-16
    nav=6.02e23
    kb=1.38e-23
    kbt=8.617e-5
    m=63.55
    m=m*1.e-3/nav
    t0=0.0
    sigma=1e10*sqrt(kb*t0/m)
    m=m/16
    Ron=6.0
    Roff=7.0
    sx=36.15
    sy=36.15
    sz=32.535
    a1=0.0
    a0=0.0854
    p=10.939
    r0=2.5563
    xi=1.2243
    q=2.2799
    count=0
    count1=0
    Vmax = -2e3
    dVmax = 0.2
    Vkvmax = 0.0
    dVkvmax = 1e11
    betta=1e12
    t=t0
    sigma1=0d0
    qu=1.0
    Rf=0.01

    call Maxwell(natom,vx,sigma)
    call Maxwell(natom,vy,sigma)
    call Maxwell(natom,vz,sigma)
    write (8,*) ";epot;ekin;etot"
    write (11,*) ";t"   
    !================================Начало эволюции системы=============================================================================================
    flag=1
    do k=1,100
        print *, k

        do i=1,np
            ip(i)=0
            jp(i)=0
        end do
        np=0

        do i=1,natom
            Eba(i)=0
            fx(i)=0
            fy(i)=0
            fz(i)=0
        end do
        
        call Verlett(natom,x,y,z,vx,vy,vz,ax,ay,az,h,sx,sy,sz)
        call Pairs(natom,x,y,z,sx,sy,sz,dr,dx,dy,dz,ip,jp,np,Roff)
        call Energy(natom,np,ip,jp,Ron,Roff,xi,q,r0,a1,a0,p,dr,Eba,er,ebound)
        call Force(np, ip, jp, Ron, Roff, q, xi, dr, r0, a1, p, a0, Eba, qexp1, pexp1, dx, dy, dz, fx, fy, fz)
        !==========================Расчёт ускорений и скоростей на шаге t+dt=============================================================================
        do i=1,natom
            ax(i)=fx(i)/m
            ay(i)=fy(i)/m
            az(i)=fz(i)/m

            vx(i)=vx(i)+0.5*(axt(i)+ax(i))*h
            vy(i)=vy(i)+0.5*(ayt(i)+ay(i))*h
            vz(i)=vz(i)+0.5*(azt(i)+az(i))*h

            axt(i)=ax(i)
            ayt(i)=ay(i)
            azt(i)=az(i)

        end do
        
        ekin = KineticEnergy(natom, vx, vy, vz)

        t=ekin*1e-20
        t=1.05e-25*t/2000/2
        t=t*2/3/1.38e-23
        ekin=0.5*ekin*m
        epot=er-ebound
        etot=epot+ekin

        write (11,*) ";",t        
        write (8,*) ";",epot,";",ekin,";",etot
        
        ekin=0
        ebound=0
        er=0
    end do
    
    write(10,*) ";x;y;z"
    write(12,*) ";vx;vy;vz"
    do i=1,natom
        write (10,*) ";",x(i),";",y(i),";",z(i)
        write (12,*) ";",vx(i),";",vy(i),";",vz(i)
    end do

    end program prog
    

    function obrez(r,Ron,Roff) result(res)
    real r,Ron,Roff,res

    if (Ron==Roff) then
        if (r<Roff) then
            res=1
        else
            res=0
        end if
    else
        if (r>Roff) then
            res=0
        end if

        if (r<=Ron) then
            res=1
        end if

        if ((r>Ron).and.(r<=Roff)) then
            res=(Roff**2-r**2)**2*(Roff**2-3*Ron**2+2*r**2)/(Roff**2-Ron**2)**3
        end if
    end if
    end function
    
    function KineticEnergy(natom, vx, vy, vz) result(Energy)
    real vx(4000), vy(4000), vz(4000), Energy
    integer natom
    
    Energy = 0
    
    do i = 1, natom
        Energy = Energy + vx(i) ** 2 + vy(i) ** 2 + vz(i) ** 2
    end do
    
    end function
    
    
    function proizv_obrez(r,Ron,Roff) result(res)
    real r,Ron,Roff,res

    if (Ron==Roff) then
        res=0
    else
        if ((r<=Ron).or.(r>=Roff)) then
            res=0
        else
            res=-12*r*(Roff**2-r**2)*(r**2-Ron**2)/(Roff**2-Ron**2)**3
        end if
    end if

    end function

    subroutine Verlett(natom,x,y,z,vx,vy,vz,ax,ay,az,h,sx,sy,sz)
    integer i,natom
    real x(4000),y(4000),z(4000),vx(4000),vy(4000),vz(4000)
    real ax(4000),ay(4000),az(4000),sx,sy,sz,h
    !==========================Интегрирование уравнений движения(расчёт координат)=======================================================================
    do i=1,natom
        x(i)=x(i)+vx(i)*h+0.5*ax(i)*h**2
        x(i)=x(i)-sx*nint(x(i)/sx-0.5)

        y(i)=y(i)+vy(i)*h+0.5*ay(i)*h**2
        y(i)=y(i)-sy*nint(y(i)/sy-0.5)

        z(i)=z(i)+vz(i)*h+0.5*az(i)*h**2
        z(i)=z(i)-sz*nint(z(i)/sz-0.5)
    end do
    return
    end

    subroutine Pairs(natom,x,y,z,sx,sy,sz,dr,dx,dy,dz,ip,jp,np,Roff)
    integer natom,np,i ,j
    real x(4000),y(4000),z(4000)
    real sx,sy,sz,Roff,drt
    real dx(200000),dy(200000),dz(200000),dr(200000)
    integer ip(200000),jp(200000)
    !================================Спаривание атомов=======================================================================================
    do i=1,natom-1
        do j=i+1,natom

            dxt=x(i)-x(j)
            dxt=dxt-sx*nint(dxt/sx)

            dyt=y(i)-y(j)
            dyt=dyt-sy*nint(dyt/sy)

            dzt=z(i)-z(j)
            dzt=dzt-sz*nint(dzt/sz)

            drt=dxt**2+dyt**2+dzt**2
            drt=sqrt(drt)

            if ((drt.lt.Roff).and.(drt.gt.1e-7)) then
                np=np+1
                ip(np)=i
                jp(np)=j
                dx(np)=dxt
                dy(np)=dyt
                dz(np)=dzt
                dr(np)=drt
            end if
        end do
    end do
    return
    end

    subroutine Energy(natom, np,ip,jp,Ron,Roff,xi,q,r0,a1,a0,p,dr,Eba,er,ebound)
    integer natom, np,i,ip(200000),jp(200000)
    real Ron,Roff,xi,q,r0,a1,a0,p,ebound,er
    real Eba(4000)
    real qexp1(200000),pexp1(200000),w(200000),dr(200000)
    !==========================Расчёт энергий=============================================================================================

    do i=1,np
        w(i) = obrez(dr(i),Ron,Roff)
        qexp1(i)=-(xi**2)*exp(-2*q*((dr(i)/r0)-1))
        pexp1(i)=(a1*(dr(i)/r0-1)+a0)*exp(-p*((dr(i)/r0)-1))
        Eba(ip(i))=Eba(ip(i))-qexp1(i)*w(i)
        Eba(jp(i))=Eba(jp(i))-qexp1(i)*w(i)
        er=er+2*pexp1(i)*w(i)
    end do

    do i=1,natom
        Eba(i)=sqrt(Eba(i))
        ebound=ebound+Eba(i)
    end do
    
    return
    end
    
    subroutine Force(np, ip, jp, Ron, Roff, q, xi, dr, r0, a1, p, a0, Eba, qexp1, pexp1, dx, dy, dz, fx, fy, fz)
    real dw(200000), qexp, pexp, f(200000)
    !=========================input=======================================================================================================
    integer np, ip(200000), jp(200000)
    real Ron, Roff, q, xi, r0, a1, p, a0
    real dr(200000), Eba(4000), qexp1(200000), pexp1(200000), dx(200000), dy(200000), dz(200000)
    !=========================output======================================================================================================
    real fx(4000), fy(4000), fz(4000)

    !==========================Расчёт сил=================================================================================================
    do i=1,np
        dw(i)=proizv_obrez(dr(i),Ron,Roff)
        qexp=-(q*xi**2)*exp(-2*q*((dr(i)/r0)-1))
        pexp=(a1*p*(dr(i)/r0-1)+a0*p-a1)*exp(-p*((dr(i)/r0)-1))

        f(i)=((1/Eba(ip(i))+1/Eba(jp(i)))*(qexp+qexp1(i)*dw(i)/2)+2*(pexp+pexp1(i)*dw(i)))/r0/dr(i)

        fx(ip(i))=fx(ip(i))+f(i)*dx(i)
        fx(jp(i))=fx(jp(i))-f(i)*dx(i)

        fy(ip(i))=fy(ip(i))+f(i)*dy(i)
        fy(jp(i))=fy(jp(i))-f(i)*dy(i)

        fz(ip(i))=fz(ip(i))+f(i)*dz(i)
        fz(jp(i))=fz(jp(i))-f(i)*dz(i)
    end do
    
    return
    end
    
    
    subroutine Maxwell(n,x,sigma)
    implicit real (a-h,o-z)
    dimension x(n)
    real sigma
    do i=1,n,2
        call random_number (u1)
        call random_number (u2)
        x(i)=sigma*cos(2*3.1415*u1)*sqrt(-2*log(u2))
        x(i+1)=sigma*sin(2*3.1415*u1)*sqrt(-2*log(u2))
    end do
    return
    end
