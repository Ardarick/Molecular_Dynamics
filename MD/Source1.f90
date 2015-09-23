
    program prog

    implicit none

    real h,m,sx,sy,sz,r0,a1,a0,p,xi,q,nav,kb,rs,ekin,etot,epot,ebound,er,v,kbt,tau,qu
    real dxt,dyt,dzt,drt,sigma,Ron,Roff,pexp,qexp,vscal,sigma1,R,Rf
    real kinetic_energy, cut_off, cut_off_div
    real x(4000),xt(4000),y(4000),yt(4000),z(4000),zt(4000)
    real vx(4000),vxt(4000),vy(4000),vyt(4000),vz(4000),vzt(4000)
    real ax(4000),axt(4000),ay(4000),ayt(4000),az(4000),azt(4000)
    real Eba(4000),fx(4000),fy(4000),fz(4000)
    real dy(200000),dz(200000),dx(200000),dr(200000),pexp1(200000),qexp1(200000)
    real fr(200000),fb(200000),f(200000),w(200000),dw(200000)
    real(8) t,t_ref
    integer i,j,k,natom,np,ip(200000),jp(200000),isrt(4000),nstat,count,iscal,count1

    open(8,file="Result.csv")
    !open(9,file="InputRelax.txt")
    open(9,file="input.txt")
    open(10,file="Output.csv")
    open(11,file="Temperature.csv")
    open(12,file="Velocity.csv")
    open(13, file = "Statistics.csv")

    natom = 2000
    nstat = 0
    iscal = 100
    
    do i = 1, natom
        read(9,*) x(i), y(i), z(i), vx(i), vy(i), vz(i), isrt(i)
    end do

    do i = 1, natom
        ax(i) = 0
        ay(i) = 0
        az(i) = 0
    end do
    
    !====================Константы и параметры==============================================================================================
    h = 1.0e-14
    nav = 6.02e23
    kb = 1.38e-23
    kbt = 8.617e-5
    m = 63.55 * 1.e-3 / nav
    t_ref = 300.0
    sigma = sqrt(kb * t_ref / m)
 
    Ron = 6.0
    Roff = 7.0
    sx = 36.15
    sy = 36.15
    sz = 32.535
    a1 = 0.0
    a0 = 0.0854
    p = 10.939
    r0 = 2.5563
    xi = 1.2243
    q = 2.2799
    tau = 1e-12
    sigma1 = 0d0
    qu = 1.0
    Rf = 0.01
    count1 = 0
    count =  0
    
    call maxwell(natom,vx,sigma)
    call maxwell(natom,vy,sigma)
    call maxwell(natom,vz,sigma)
    
    do i = 1, natom
        vx(i) = vx(i) * 1e10
        vy(i) = vy(i) * 1e10
        vz(i) = vz(i) * 1e10
    end do

    write (8,*) ";epot;ekin;etot"
    write (11,*) ";t"   
    write(13, *) ";vx;vy;vz"
    !================================Начало эволюции системы=============================================================================================

    do k=1,10000
        print *, '==============================================='
        print *, 'Step number', k

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
        
        call verlet_coords(natom, x, y, z, vx, vy, vz, ax, ay, az, h, sx, sy, sz)
        call pairs(natom, x, y, z, sx, sy, sz, dr, dx, dy, dz, ip, jp, np, Roff)
        call energy(natom, np, ip, jp, Ron, Roff, xi, q, r0, a1, a0, p, dr, Eba, er, ebound)
        call force(np, ip, jp, Ron, Roff, q, xi, dr, r0, a1, p, a0, Eba, qexp1, pexp1, dx, dy, dz, fx, fy, fz)
        call verlet_velocity(natom, m, fx, fy, fz, vx, vy, vz, h, axt, ayt, azt, ax, ay, az)

        
        count = count + 1
        if (count.eq.100) then
            count = 0
            do i=1,natom
                write (13,*) ";",vx(i),";",vy(i),";",vz(i)
            end do
        end if
        
        ekin = kinetic_energy(natom, m, vx, vy, vz)
       
        t = 2 * ekin / (3 * kbt) / natom
        print *, 'Temperature = ', t, ' K'
        
        epot=er-ebound
        etot=epot+ekin

        write (11,*) ";",t        
        !count1 = count1 + 1
        !if (count1.eq.iscal) then
            !count1 = 0
            !call thermo_scale(natom, vx, vy, vz, t, t_ref)
        !end if
        
        !call thermo_andersen(natom, vx, vy, vz, t_ref, Rf, m)
        call thermo_berendsen(natom, vx, vy, vz, t_ref, t, tau, h)
        
        
        print *, 'Potential energy = ', epot, ' eV'
        print *, 'Kinetic energy = ', ekin, ' eV'
        print *, 'Full energy = ', etot, ' eV'
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
    
    

    function cut_off(r,Ron,Roff) result(res)
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
    
    function kinetic_energy(natom, m, vx, vy, vz) result(Energy)
    real vx(4000), vy(4000), vz(4000), Energy, m
    integer natom
    
    Energy = 0
    
    do i = 1, natom
        Energy = Energy + vx(i) ** 2 + vy(i) ** 2 + vz(i) ** 2
    end do
    
    Energy = 0.5 * Energy * m / 16
    
    end function
    
    
    function cut_off_div(r,Ron,Roff) result(res)
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

    subroutine verlet_coords(natom,x,y,z,vx,vy,vz,ax,ay,az,h,sx,sy,sz)
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

    subroutine pairs(natom,x,y,z,sx,sy,sz,dr,dx,dy,dz,ip,jp,np,Roff)
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

    subroutine energy(natom, np,ip,jp,Ron,Roff,xi,q,r0,a1,a0,p,dr,Eba,er,ebound)
    integer natom, np,i,ip(200000),jp(200000)
    real Ron,Roff,xi,q,r0,a1,a0,p,ebound,er
    real Eba(4000)
    real qexp1(200000),pexp1(200000),w(200000),dr(200000)
    !==========================Расчёт энергий=============================================================================================

    do i=1,np
        w(i) = cut_off(dr(i),Ron,Roff)
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
    
    subroutine force(np, ip, jp, Ron, Roff, q, xi, dr, r0, a1, p, a0, Eba, qexp1, pexp1, dx, dy, dz, fx, fy, fz)
    real dw(200000), qexp, pexp, f(200000)
    !=========================input=======================================================================================================
    integer np, ip(200000), jp(200000)
    real Ron, Roff, q, xi, r0, a1, p, a0
    real dr(200000), Eba(4000), qexp1(200000), pexp1(200000), dx(200000), dy(200000), dz(200000)
    !=========================output======================================================================================================
    real fx(4000), fy(4000), fz(4000)

    !==========================Расчёт сил=================================================================================================
    do i=1,np
        dw(i)=cut_off_div(dr(i),Ron,Roff)
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
    
    subroutine verlet_velocity(natom, m, fx, fy, fz, vx, vy, vz, h, axt, ayt, azt, ax, ay, az)
    integer natom
    real h, m
    real fx(4000), fy(4000), fz(4000), vx(4000), vy(4000), vz(4000)
    real ax(4000), ay(4000), az(4000), axt(4000), ayt(4000), azt(4000)
    !==========================Расчёт ускорений и скоростей на шаге t+dt=========================================================================
    do i=1,natom
        ax(i)=fx(i)/m * 16
        ay(i)=fy(i)/m * 16
        az(i)=fz(i)/m * 16

        vx(i)=vx(i)+0.5*(axt(i)+ax(i))*h
        vy(i)=vy(i)+0.5*(ayt(i)+ay(i))*h
        vz(i)=vz(i)+0.5*(azt(i)+az(i))*h

        axt(i)=ax(i)
        ayt(i)=ay(i)
        azt(i)=az(i)

    end do
    
    return
    end
    
    
    subroutine maxwell(n,x,sigma)
    implicit real (a-h,o-z)
    dimension x(n)
    real sigma, pi
    pi = 3.1415
    do i=1,n,2
        call random_number (u1)
        call random_number (u2)
        x(i)=sigma*cos(2*pi*u1)*sqrt(-2*log(u2))
        x(i+1)=sigma*sin(2*pi*u1)*sqrt(-2*log(u2))
    end do
    return
    end
!======================================Термостаты===============================================================================================
    
    subroutine thermo_scale(natom, vx, vy, vz, t, t_ref)
    integer natom
    real vx(4000), vy(4000), vz(4000), scale
    real(8) t_ref, t
  

    scale = sqrt(t_ref / t)

    do i = 1, natom
        vx(i) = vx(i) * scale
        vy(i) = vy(i) * scale
        vz(i) = vz(i) * scale
    end do

    return
    end
    
    
    subroutine thermo_andersen(natom, vx, vy, vz, t_ref, Rf, m)
    integer natom
    real vx(4000), vy(4000), vz(4000), Rf, u1, u2, sigma, kb, R, m, pi
    real(8) t_ref
    kb = 1.38e-23
    pi = 3.1415
    sigma = sqrt(kb * t_ref / m)
    
    do i = 1, natom
        call random_number (R)
        if (R.lt.Rf) then
            call random_number (u1)
            call random_number (u2)
            vx(i)=sigma*cos(2*pi*u1)*sqrt(-2*log(u2)) * 1e10
            
            call random_number (u1)
            call random_number (u2)
            vy(i)=sigma*cos(2*pi*u1)*sqrt(-2*log(u2)) * 1e10
            
            call random_number (u1)
            call random_number (u2)
            vz(i)=sigma*cos(2*pi*u1)*sqrt(-2*log(u2)) * 1e10
        end if
    end do
    
    return
    end
    
    
    subroutine thermo_berendsen(natom, vx, vy, vz, t_ref, t, tau, h)
    integer natom
    real vx(4000), vy(4000), vz(4000), tau, h, lambda
    real(8) t, t_ref
    
    lambda = 1 + h / (2 * tau) * (t_ref / t - 1)
    
    do i = 1, natom
        vx(i) = vx(i) * lambda
        vy(i) = vy(i) * lambda
        vy(i) = vy(i) * lambda
    end do
    
    return
    end